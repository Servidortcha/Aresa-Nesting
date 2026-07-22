"""
dxf_io.py
Lector y escritor mínimo de archivos DXF, sin dependencias externas
(no requiere ezdxf, ya que no hay acceso a internet para instalarlo).

Soporta, en lectura, las entidades más comunes en piezas de corte láser:
  - LWPOLYLINE (con arcos vía "bulge")
  - POLYLINE / VERTEX (formato clásico)
  - LINE
  - CIRCLE
  - ARC
  - SPLINE — evaluada de verdad con el algoritmo de De Boor (no solo sus
    puntos de control), soporta curvas racionales (con pesos) y no
    racionales, usando el vector de nudos ("knots") del archivo.

Muchos programas de CAD (ej. exportaciones vistas de SolidWorks/Rhino/
AutoCAD) NO agrupan el contorno de una pieza en una sola polilínea
cerrada: lo parten en varias entidades sueltas (LINE + SPLINE + ARC) que
hay que unir por sus extremos para reconstruir el contorno cerrado. Este
módulo hace ese "encadenado" automáticamente.

En escritura, genera un DXF (formato R2000 / AC1015) válido con entidades
LWPOLYLINE, suficiente para abrir en LightBurn, AutoCAD, Illustrator, etc.

Limitaciones conocidas (documentadas para el usuario):
  - No soporta BLOCK/INSERT (piezas definidas como bloques anidados).
  - Se asume que cada DXF de entrada contiene UNA sola pieza (el contorno
    cerrado de mayor área = pieza; el resto de contornos cerrados = agujeros).
"""
import math


def _read_group_codes(path):
    """Lee un archivo DXF como lista de pares (codigo:int, valor:str)."""
    with open(path, "r", errors="ignore") as f:
        lines = [l.rstrip("\n").rstrip("\r") for l in f]
    pairs = []
    i = 0
    while i < len(lines) - 1:
        code_line = lines[i].strip()
        value_line = lines[i + 1]
        try:
            code = int(code_line)
        except ValueError:
            i += 2
            continue
        pairs.append((code, value_line))
        i += 2
    return pairs


def _iter_entities(pairs):
    """Devuelve lista de (tipo_entidad, lista_de_codigos) por cada entidad
    dentro de la seccion ENTITIES."""
    in_entities = False
    current_type = None
    current = []
    out = []
    for code, value in pairs:
        if code == 2 and value == "ENTITIES":
            in_entities = True
            continue
        if code == 0 and value == "ENDSEC":
            if in_entities and current_type is not None:
                out.append((current_type, current))
                current_type, current = None, []
            in_entities = False
            continue
        if not in_entities:
            continue
        if code == 0:
            if current_type is not None:
                out.append((current_type, current))
            current_type = value
            current = []
        else:
            current.append((code, value))
    if current_type is not None:
        out.append((current_type, current))
    return out


def _bulge_to_arc_points(p1, p2, bulge, segments=12):
    """Convierte un segmento con 'bulge' (arco en LWPOLYLINE) a una polilinea."""
    if abs(bulge) < 1e-9:
        return [p2]
    x1, y1 = p1
    x2, y2 = p2
    theta = 4 * math.atan(bulge)
    chord = math.hypot(x2 - x1, y2 - y1)
    if chord < 1e-9:
        return [p2]
    radius = chord / (2 * math.sin(theta / 2))
    mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
    dx, dy = x2 - x1, y2 - y1
    h = radius * math.cos(theta / 2)
    nx, ny = -dy / chord, dx / chord
    sign = 1 if bulge > 0 else -1
    cx = mx + sign * h * nx
    cy = my + sign * h * ny
    start_ang = math.atan2(y1 - cy, x1 - cx)
    pts = []
    n = max(2, int(abs(theta) / (math.pi / segments)) + 1)
    for k in range(1, n + 1):
        a = start_ang + theta * (k / n)
        pts.append((cx + radius * math.cos(a), cy + radius * math.sin(a)))
    return pts


# ---------------------------------------------------------------------------
# Evaluacion de curvas SPLINE (B-spline / NURBS) via algoritmo de De Boor
# ---------------------------------------------------------------------------

def _find_span(u, degree, knots, n_ctrl):
    if u >= knots[n_ctrl]:
        return n_ctrl - 1
    if u <= knots[degree]:
        return degree
    lo, hi = degree, n_ctrl
    mid = (lo + hi) // 2
    while u < knots[mid] or u >= knots[mid + 1]:
        if u < knots[mid]:
            hi = mid
        else:
            lo = mid
        mid = (lo + hi) // 2
    return mid


def _bspline_point(u, degree, knots, ctrl_pts, weights):
    n_ctrl = len(ctrl_pts)
    span = _find_span(u, degree, knots, n_ctrl)
    d = []
    for j in range(degree + 1):
        idx = span - degree + j
        w = weights[idx]
        d.append([ctrl_pts[idx][0] * w, ctrl_pts[idx][1] * w, w])
    for r in range(1, degree + 1):
        for j in range(degree, r - 1, -1):
            i = span - degree + j
            denom = knots[i + degree - r + 1] - knots[i]
            alpha = 0.0 if abs(denom) < 1e-12 else (u - knots[i]) / denom
            d[j][0] = (1 - alpha) * d[j - 1][0] + alpha * d[j][0]
            d[j][1] = (1 - alpha) * d[j - 1][1] + alpha * d[j][1]
            d[j][2] = (1 - alpha) * d[j - 1][2] + alpha * d[j][2]
    w = d[degree][2]
    if abs(w) < 1e-12:
        w = 1.0
    return (d[degree][0] / w, d[degree][1] / w)


def _sample_spline(degree, knots, ctrl_pts, weights, closed_flag, n_samples=40):
    n_ctrl = len(ctrl_pts)
    if closed_flag:
        # curva periodica cerrada: el vector de nudos ya viene preparado por
        # la mayoria de exportadores para recorrer todo el rango util.
        u0, u1 = knots[degree], knots[n_ctrl]
    else:
        u0, u1 = knots[degree], knots[n_ctrl]
    pts = []
    for k in range(n_samples + 1):
        u = u0 + (u1 - u0) * (k / n_samples)
        u = min(max(u, knots[0]), knots[-1] - 1e-9)
        pts.append(_bspline_point(u, degree, knots, ctrl_pts, weights))
    return pts


def _parse_spline(codes, n_samples=40):
    degree = 3
    knots = []
    ctrl_pts = []
    weights = []
    fit_pts = []
    flags = 0
    cur_ctrl = None
    cur_fit = None
    for code, value in codes:
        if code == 70:
            flags = int(float(value))
        elif code == 71:
            degree = int(float(value))
        elif code == 40:
            knots.append(float(value))
        elif code == 41:
            weights.append(float(value))
        elif code == 10:
            if cur_ctrl:
                ctrl_pts.append((cur_ctrl["x"], cur_ctrl["y"]))
            cur_ctrl = {"x": float(value)}
        elif code == 20 and cur_ctrl is not None:
            cur_ctrl["y"] = float(value)
        elif code == 11:
            if cur_fit:
                fit_pts.append((cur_fit["x"], cur_fit["y"]))
            cur_fit = {"x": float(value)}
        elif code == 21 and cur_fit is not None:
            cur_fit["y"] = float(value)
    if cur_ctrl:
        ctrl_pts.append((cur_ctrl["x"], cur_ctrl["y"]))
    if cur_fit:
        fit_pts.append((cur_fit["x"], cur_fit["y"]))

    closed = bool(flags & 1)

    if len(ctrl_pts) >= 2 and len(knots) == len(ctrl_pts) + degree + 1:
        if not weights:
            weights = [1.0] * len(ctrl_pts)
        try:
            return _sample_spline(degree, knots, ctrl_pts, weights, closed, n_samples)
        except Exception:
            pass
    # Respaldo: si no se pudo evaluar como NURBS (datos incompletos), usar
    # los puntos de ajuste (fit points) si vienen, o si no, los de control.
    if fit_pts:
        return fit_pts
    return ctrl_pts


# ---------------------------------------------------------------------------
# Extraccion de segmentos "crudos" (abiertos o cerrados) por entidad
# ---------------------------------------------------------------------------

def _entity_to_segment(etype, codes, arc_segments=16, spline_samples=40):
    """Devuelve (points, is_closed) o None si la entidad no aporta geometria
    de corte (texto, cotas, etc.)."""
    if etype == "LINE":
        x1 = y1 = x2 = y2 = 0.0
        for code, value in codes:
            if code == 10:
                x1 = float(value)
            elif code == 20:
                y1 = float(value)
            elif code == 11:
                x2 = float(value)
            elif code == 21:
                y2 = float(value)
        return [(x1, y1), (x2, y2)], False

    if etype == "CIRCLE":
        cx = cy = 0.0
        r = 0.0
        for code, value in codes:
            if code == 10:
                cx = float(value)
            elif code == 20:
                cy = float(value)
            elif code == 40:
                r = float(value)
        n = arc_segments * 4
        pts = [(cx + r * math.cos(2 * math.pi * k / n),
                cy + r * math.sin(2 * math.pi * k / n)) for k in range(n)]
        return pts, True

    if etype == "ARC":
        cx = cy = 0.0
        r = 0.0
        a1, a2 = 0.0, 360.0
        for code, value in codes:
            if code == 10:
                cx = float(value)
            elif code == 20:
                cy = float(value)
            elif code == 40:
                r = float(value)
            elif code == 50:
                a1 = float(value)
            elif code == 51:
                a2 = float(value)
        if a2 < a1:
            a2 += 360.0
        n = max(2, int(arc_segments * (a2 - a1) / 90.0) + 1)
        pts = []
        for k in range(n + 1):
            a = math.radians(a1 + (a2 - a1) * k / n)
            pts.append((cx + r * math.cos(a), cy + r * math.sin(a)))
        return pts, False

    if etype == "SPLINE":
        pts = _parse_spline(codes, n_samples=spline_samples)
        flags = 0
        for code, value in codes:
            if code == 70:
                flags = int(float(value))
        return pts, bool(flags & 1)

    if etype == "LWPOLYLINE":
        verts = []
        cur = {}
        closed = False
        for code, value in codes:
            if code == 70:
                closed = bool(int(float(value)) & 1)
            elif code == 10:
                if cur:
                    verts.append(cur)
                cur = {"x": float(value)}
            elif code == 20:
                cur["y"] = float(value)
            elif code == 42:
                cur["bulge"] = float(value)
        if cur:
            verts.append(cur)
        pts = []
        for idx, v in enumerate(verts):
            p1 = (v["x"], v["y"])
            pts.append(p1)
            bulge = v.get("bulge", 0.0)
            if bulge:
                nxt = None
                if idx + 1 < len(verts):
                    nxt = verts[idx + 1]
                elif closed:
                    nxt = verts[0]
                if nxt is not None:
                    p2 = (nxt["x"], nxt["y"])
                    pts.extend(_bulge_to_arc_points(p1, p2, bulge, arc_segments)[:-1])
        return pts, closed

    return None


def _entities_with_vertex_groups(entities):
    """Reconstruye POLYLINE clasicos (POLYLINE + VERTEX... + SEQEND) como
    una sola entidad 'sintetica' con sus puntos ya extraidos."""
    out = []
    i = 0
    while i < len(entities):
        etype, codes = entities[i]
        if etype == "POLYLINE":
            closed = False
            for code, value in codes:
                if code == 70:
                    closed = bool(int(float(value)) & 1)
            pts = []
            j = i + 1
            while j < len(entities) and entities[j][0] == "VERTEX":
                vx = vy = None
                for code, value in entities[j][1]:
                    if code == 10:
                        vx = float(value)
                    elif code == 20:
                        vy = float(value)
                if vx is not None and vy is not None:
                    pts.append((vx, vy))
                j += 1
            if j < len(entities) and entities[j][0] == "SEQEND":
                j += 1
            out.append(("_POLYLINE_RESOLVED", (pts, closed)))
            i = j
            continue
        out.append((etype, codes))
        i += 1
    return out


def _points_equal(a, b, tol):
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol


def _chain_open_segments(segments, tol):
    """Une segmentos abiertos (listas de puntos) por sus extremos hasta
    formar contornos cerrados. Devuelve (loops_cerrados, restos_abiertos)."""
    remaining = [list(s) for s in segments]
    closed_loops = []
    leftover_open = []

    while remaining:
        chain = remaining.pop(0)
        if len(chain) < 2:
            continue
        progressed = True
        while progressed:
            progressed = False
            if _points_equal(chain[0], chain[-1], tol) and len(chain) > 2:
                break
            start, end = chain[0], chain[-1]
            for idx in range(len(remaining)):
                other = remaining[idx]
                if _points_equal(end, other[0], tol):
                    chain = chain + other[1:]
                elif _points_equal(end, other[-1], tol):
                    chain = chain + list(reversed(other))[1:]
                elif _points_equal(start, other[-1], tol):
                    chain = other[:-1] + chain
                elif _points_equal(start, other[0], tol):
                    chain = list(reversed(other))[:-1] + chain
                else:
                    continue
                remaining.pop(idx)
                progressed = True
                break
        if _points_equal(chain[0], chain[-1], tol) and len(chain) > 2:
            closed_loops.append(chain[:-1])
        else:
            leftover_open.append(chain)
    return closed_loops, leftover_open


def _polygon_area(points):
    a = 0.0
    n = len(points)
    for i in range(n):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def load_dxf_piece(path, arc_segments=16, spline_samples=40, join_tolerance=0.01):
    """Carga un DXF y devuelve (outer, holes): el contorno principal como
    lista de puntos (x,y) cerrada, y una lista de contornos interiores
    (agujeros). Reconstruye automaticamente piezas cuyo contorno viene
    partido en varias entidades sueltas (LINE/SPLINE/ARC) uniendolas por
    sus extremos."""
    pairs = _read_group_codes(path)
    entities = _iter_entities(pairs)
    entities = _entities_with_vertex_groups(entities)

    closed_loops = []
    open_segments = []

    for etype, payload in entities:
        if etype == "_POLYLINE_RESOLVED":
            pts, closed = payload
            if len(pts) < 2:
                continue
            if closed or _points_equal(pts[0], pts[-1], join_tolerance):
                closed_loops.append(pts[:-1] if _points_equal(pts[0], pts[-1], join_tolerance) else pts)
            else:
                open_segments.append(pts)
            continue

        result = _entity_to_segment(etype, payload, arc_segments, spline_samples)
        if result is None:
            continue
        pts, closed = result
        if len(pts) < 2:
            continue
        if closed or _points_equal(pts[0], pts[-1], join_tolerance):
            closed_loops.append(pts[:-1] if _points_equal(pts[0], pts[-1], join_tolerance) else pts)
        else:
            open_segments.append(pts)

    joined_loops, leftover = _chain_open_segments(open_segments, join_tolerance)
    closed_loops.extend(joined_loops)

    closed_loops = [l for l in closed_loops if len(l) >= 3]

    if not closed_loops:
        raise ValueError(
            "El DXF no contiene ningun contorno cerrado reconocible (%s). "
            "Si el dibujo tiene tramos sueltos que no llegan a cerrar el "
            "contorno, revisa que no haya micro-huecos entre segmentos." % path
        )

    closed_loops.sort(key=_polygon_area, reverse=True)
    outer = closed_loops[0]
    holes = closed_loops[1:]
    warning = None
    if leftover:
        warning = (
            "%d tramo(s) sueltos no pudieron unirse a un contorno cerrado "
            "(posible micro-hueco entre segmentos); se ignoraron." % len(leftover)
        )
    return outer, holes, warning

# ---------------------------------------------------------------------------
# Escritura de DXF de salida (layout final con todas las piezas ubicadas)
# ---------------------------------------------------------------------------

_HEADER = """0
SECTION
2
HEADER
9
$INSUNITS
70
4
0
ENDSEC
0
SECTION
2
TABLES
0
TABLE
2
LAYER
0
LAYER
2
0
70
0
62
7
6
CONTINUOUS
0
LAYER
2
PIEZAS
70
0
62
1
6
CONTINUOUS
0
LAYER
2
CHAPA
70
0
62
5
6
CONTINUOUS
0
ENDTAB
0
ENDSEC
0
SECTION
2
ENTITIES
"""

_FOOTER = """0
ENDSEC
0
EOF
"""


def _polyline_entity(points, layer="PIEZAS", closed=True):
    lines = ["0", "LWPOLYLINE", "8", layer, "90", str(len(points)),
             "70", "1" if closed else "0"]
    for (x, y) in points:
        lines += ["10", f"{x:.6f}", "20", f"{y:.6f}"]
    return "\n".join(lines) + "\n"


def save_dxf_layout(path, placed_pieces, sheet_w, sheet_h, draw_sheet_outline=True):
    """
    placed_pieces: lista de dicts con:
        {"points": [(x,y), ...], "holes": [[(x,y),...], ...]}
      ya en coordenadas absolutas dentro de la chapa.
    """
    body = []
    if draw_sheet_outline:
        sheet_pts = [(0, 0), (sheet_w, 0), (sheet_w, sheet_h), (0, sheet_h)]
        body.append(_polyline_entity(sheet_pts, layer="CHAPA", closed=True))
    for piece in placed_pieces:
        body.append(_polyline_entity(piece["points"], layer="PIEZAS", closed=True))
        for hole in piece.get("holes", []):
            body.append(_polyline_entity(hole, layer="PIEZAS", closed=True))
    with open(path, "w") as f:
        f.write(_HEADER)
        f.write("".join(body))
        f.write(_FOOTER)
