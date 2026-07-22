"""
nesting.py
Algoritmo de nesting: ubica todas las piezas (con sus cantidades y rotacion
libre) dentro de una o mas chapas, minimizando desperdicio.

Estrategia:
  1. Cada pieza se rasteriza a una mascara binaria por cada angulo candidato
     (rotacion discretizada en pasos de `rotation_step_deg`; con un paso
     chico esto se aproxima a "rotacion libre").
  2. El espaciado/kerf configurado se aplica dilatando la mascara la mitad
     del espaciado (asi, cuando dos piezas quedan "tocandose" en la grilla,
     la distancia real entre sus contornos es el espaciado pedido).
  3. Para encontrar la posicion valida (sin colision) se usa correlacion
     2D via FFT (scipy.signal.fftconvolve) entre la ocupacion de la chapa
     y la mascara de la pieza: es mucho mas rapido que probar posicion por
     posicion con loops de Python.
  4. Heuristica de orden: piezas de mayor area primero (heuristica estandar
     de nesting: las piezas grandes son las que mas restringen el layout).
  5. Heuristica de posicion: entre todas las posiciones/angulos validos, se
     elige la mas "abajo-izquierda" (bottom-left) para compactar el layout.
"""
import math
import numpy as np
from scipy.signal import fftconvolve

import geometry as geo


def _build_piece_variants(points, holes, angles_deg, cell_mm, spacing_mm):
    """Para cada angulo candidato, devuelve la mascara rasterizada (con
    dilatacion por espaciado) junto con los puntos normalizados y el
    origen de la mascara en espacio 'normalizado' (bbox del poligono
    rotado empieza en (0,0))."""
    dilation_cells = max(0, int(round((spacing_mm / 2.0) / cell_mm)))
    variants = {}
    for ang in angles_deg:
        rotated = geo.rotate_points(points, ang)
        rotated_holes = [geo.rotate_points(h, ang) for h in (holes or [])]
        norm_pts, offset = geo.normalize_to_origin(rotated)
        norm_holes = [geo.translate_points(h, -offset[0], -offset[1]) for h in rotated_holes]
        mask, mask_origin = geo.rasterize_polygon(
            norm_pts, norm_holes, cell_mm, margin_cells=dilation_cells + 1)
        mask_d = geo.dilate_mask(mask, dilation_cells)
        variants[ang] = {
            "mask": mask_d,
            "shape": mask_d.shape,
            "mask_origin": mask_origin,
            "norm_pts": norm_pts,
            "norm_holes": norm_holes,
        }
    return variants


def _find_placement(sheet_occ, mask):
    mh, mw = mask.shape
    sh, sw = sheet_occ.shape
    if mh > sh or mw > sw:
        return None
    corr = fftconvolve(sheet_occ.astype(np.float32), mask[::-1, ::-1].astype(np.float32), mode="valid")
    valid = corr < 0.5
    if not valid.any():
        return None
    ys, xs = np.nonzero(valid)
    order = np.lexsort((xs, ys))  # ordena por y asc, luego x asc -> bottom-left
    best = order[0]
    return int(ys[best]), int(xs[best])


def auto_params(parts):
    """Elige automaticamente la resolucion de grilla y el paso de rotacion
    segun el tamano de las piezas y la cantidad total de copias, para no
    tener que pedirle estos parametros tecnicos al usuario."""
    total_qty = sum(int(p["qty"]) for p in parts) or 1
    min_dim = None
    for p in parts:
        minx, miny, maxx, maxy = geo.polygon_bounds(p["points"])
        d = min(maxx - minx, maxy - miny)
        if d > 0 and (min_dim is None or d < min_dim):
            min_dim = d
    if not min_dim:
        min_dim = 50.0
    cell_mm = max(1.0, min(5.0, min_dim / 12.0))

    if total_qty > 80:
        rotation_step_deg = 30.0
    elif total_qty > 30:
        rotation_step_deg = 20.0
    elif total_qty > 10:
        rotation_step_deg = 12.0
    else:
        rotation_step_deg = 6.0
    return round(cell_mm, 2), rotation_step_deg


def _fits_within_width(points, holes, sheet_w_mm, spacing_mm, cell_mm, angles):
    variants = _build_piece_variants(points, holes, angles, cell_mm, spacing_mm)
    w_cells_limit = max(1, int(math.floor(sheet_w_mm / cell_mm)))
    return any(v["shape"][1] <= w_cells_limit for v in variants.values())


def nest_strip(parts, sheet_w_mm, spacing_mm=3.0, cell_mm=None,
               allow_rotation=True, rotation_step_deg=None,
               max_length_mm=None, progress_cb=None):
    """Acomoda todas las piezas en UNA sola tira de ANCHO FIJO
    (`sheet_w_mm`, el ancho de chapa/rollo que se compra) calculando
    automaticamente el LARGO MINIMO necesario -- asi el area total de
    material usado queda minimizada, en vez de asumir una chapa de
    tamano fijo con desperdicio.

    Devuelve (sheet, no_ubicables):
      sheet: {"width", "height" (largo minimo calculado), "pieces",
              "used_area_mm2", "utilization_pct"} o None si no hay nada
              para ubicar.
      no_ubicables: piezas que nunca entran (mas anchas que sheet_w_mm
              en cualquier rotacion), independientemente del largo.
    """
    if cell_mm is None or rotation_step_deg is None:
        auto_cell, auto_rot = auto_params(parts)
        if cell_mm is None:
            cell_mm = auto_cell
        if rotation_step_deg is None:
            rotation_step_deg = auto_rot

    angles = [0.0] if not allow_rotation else list(np.arange(0, 360, rotation_step_deg))

    placeable_parts = []
    impossible = []
    for p in parts:
        if _fits_within_width(p["points"], p.get("holes", []), sheet_w_mm, spacing_mm, cell_mm, angles):
            placeable_parts.append(p)
        else:
            impossible.append({"name": p["name"]})

    if not placeable_parts:
        return None, impossible

    total_area = sum(geo.polygon_area(p["points"]) * int(p["qty"]) for p in placeable_parts)
    guess_h = max(sheet_w_mm * 0.25, (total_area * 1.6) / sheet_w_mm if sheet_w_mm else 500.0, 80.0)
    if max_length_mm:
        guess_h = min(guess_h, max_length_mm)

    # Limite de memoria: una grilla demasiado fina sobre una chapa grande
    # puede consumir mucha RAM (cada busqueda de posicion hace una
    # convolucion FFT del tamano de la chapa). Si la estimacion de largo
    # es grande, usamos una resolucion mas gruesa para que el numero total
    # de celdas se mantenga acotado, aunque eso implique perder algo de
    # precision en el acomodo.
    MAX_GRID_CELLS = 1_500_000
    estimated_h_for_sizing = max_length_mm or (guess_h * 3.0)
    min_cell_for_memory = math.sqrt((sheet_w_mm * estimated_h_for_sizing) / MAX_GRID_CELLS)
    if min_cell_for_memory > cell_mm:
        cell_mm = round(min_cell_for_memory, 2)

    sheets, unplaced = [], []
    attempts = 0
    variants_cache = {}
    while True:
        attempts += 1
        sheets, unplaced = nest(
            placeable_parts, sheet_w_mm, guess_h, spacing_mm=spacing_mm, cell_mm=cell_mm,
            allow_rotation=allow_rotation, rotation_step_deg=rotation_step_deg,
            progress_cb=progress_cb, variants_cache=variants_cache,
        )
        ok = (len(sheets) <= 1 and not unplaced)
        reached_cap = max_length_mm and guess_h >= max_length_mm
        if ok or attempts >= 6 or reached_cap:
            break
        guess_h = min(guess_h * 1.7, max_length_mm) if max_length_mm else guess_h * 1.7

    for u in unplaced:
        impossible.append({"name": u["name"]})

    if not sheets or not sheets[0]["pieces"]:
        return None, impossible

    sheet = sheets[0]
    all_pts = [pt for p in sheet["pieces"] for pt in (p["points"] + [c for h in p["holes"] for c in h])]
    max_y = max(pt[1] for pt in all_pts)
    used_length = min(guess_h, max_y + spacing_mm / 2.0)

    used_area = sum(p["area"] for p in sheet["pieces"])
    sheet_area = sheet_w_mm * used_length if used_length else 0.0
    result = {
        "width": sheet_w_mm,
        "height": round(used_length, 1),
        "pieces": sheet["pieces"],
        "used_area_mm2": used_area,
        "utilization_pct": round(100.0 * used_area / sheet_area, 2) if sheet_area else 0.0,
    }
    return result, impossible


def nest(parts, sheet_w_mm, sheet_h_mm, spacing_mm=3.0, cell_mm=4.0,
         allow_rotation=True, rotation_step_deg=15.0, progress_cb=None,
         variants_cache=None):
    """
    parts: lista de dicts:
        {"part_id": str, "name": str, "points": [(x,y),...],
         "holes": [[(x,y),...], ...], "qty": int}
    Devuelve: lista de "sheets", cada una:
        {"width": sheet_w_mm, "height": sheet_h_mm,
         "pieces": [{"part_id","name","angle","points","holes"}, ...],
         "used_area_mm2": float, "utilization_pct": float}
    y ademas una lista de piezas que no entraron en ninguna chapa
    (por ejemplo si una pieza es mas grande que la chapa).
    """
    angles = [0.0] if not allow_rotation else list(np.arange(0, 360, rotation_step_deg))

    sheet_w_cells = max(1, int(math.floor(sheet_w_mm / cell_mm)))
    sheet_h_cells = max(1, int(math.floor(sheet_h_mm / cell_mm)))

    instances = []
    for p in parts:
        area = geo.polygon_area(p["points"])
        for _ in range(int(p["qty"])):
            instances.append({
                "part_id": p["part_id"], "name": p["name"],
                "points": p["points"], "holes": p.get("holes", []),
                "area": area,
            })
    instances.sort(key=lambda x: x["area"], reverse=True)

    # Las copias de la misma pieza comparten la misma geometria: calculamos
    # las variantes (mascara rasterizada por angulo) UNA sola vez por
    # part_id, no por cada copia -- esto es la optimizacion mas importante
    # de rendimiento (evita recalcular lo mismo N veces por N copias).
    variants_cache = {} if variants_cache is None else variants_cache

    def get_variants(inst):
        key = inst["part_id"]
        if key not in variants_cache:
            variants_cache[key] = _build_piece_variants(
                inst["points"], inst["holes"], angles, cell_mm, spacing_mm)
        return variants_cache[key]

    sheets = []  # cada uno: {"occ": ndarray, "pieces": [...]}
    unplaced = []

    total = len(instances)
    for idx, inst in enumerate(instances):
        variants = get_variants(inst)

        # ¿la pieza entra, en algun angulo, dentro de las dimensiones de la chapa?
        fits_any_sheet_size = any(
            v["shape"][0] <= sheet_h_cells and v["shape"][1] <= sheet_w_cells
            for v in variants.values()
        )
        if not fits_any_sheet_size:
            unplaced.append(inst)
            if progress_cb:
                progress_cb(idx + 1, total, placed=False)
            continue

        placed = False
        for sheet in sheets:
            best_choice = None  # (score, ang, row, col)
            cur_max_y = sheet.get("max_y", 0)
            cur_max_x = sheet.get("max_x", 0)
            for ang, v in variants.items():
                if v["shape"][0] > sheet["occ"].shape[0] or v["shape"][1] > sheet["occ"].shape[1]:
                    continue
                pos = _find_placement(sheet["occ"], v["mask"])
                if pos is None:
                    continue
                row, col = pos
                mh, mw = v["shape"]
                new_max_y = max(cur_max_y, row + mh)
                new_max_x = max(cur_max_x, col + mw)
                # score principal: area total ocupada (minimizarla es lo que
                # de verdad ahorra material); empate: menor alto, luego bottom-left
                score = (new_max_y * new_max_x, new_max_y, new_max_x, row, col)
                if best_choice is None or score < best_choice[0]:
                    best_choice = (score, ang, row, col)
            if best_choice is not None:
                _, ang, row, col = best_choice
                v = variants[ang]
                mh, mw = v["shape"]
                sheet["occ"][row:row + mh, col:col + mw] |= v["mask"]
                sheet["max_y"] = max(cur_max_y, row + mh)
                sheet["max_x"] = max(cur_max_x, col + mw)
                dx = col * cell_mm - v["mask_origin"][0]
                dy = row * cell_mm - v["mask_origin"][1]
                abs_pts = geo.translate_points(v["norm_pts"], dx, dy)
                abs_holes = [geo.translate_points(h, dx, dy) for h in v["norm_holes"]]
                sheet["pieces"].append({
                    "part_id": inst["part_id"], "name": inst["name"],
                    "angle": ang, "points": abs_pts, "holes": abs_holes,
                    "area": inst["area"],
                })
                placed = True
                break

        if not placed:
            # nueva chapa
            occ = np.zeros((sheet_h_cells, sheet_w_cells), dtype=bool)
            best_choice = None
            for ang, v in variants.items():
                if v["shape"][0] > occ.shape[0] or v["shape"][1] > occ.shape[1]:
                    continue
                pos = _find_placement(occ, v["mask"])
                if pos is None:
                    continue
                row, col = pos
                mh, mw = v["shape"]
                new_max_y = row + mh
                new_max_x = col + mw
                score = (new_max_y * new_max_x, new_max_y, new_max_x, row, col)
                if best_choice is None or score < best_choice[0]:
                    best_choice = (score, ang, row, col)
            if best_choice is not None:
                _, ang, row, col = best_choice
                v = variants[ang]
                mh, mw = v["shape"]
                occ[row:row + mh, col:col + mw] |= v["mask"]
                dx = col * cell_mm - v["mask_origin"][0]
                dy = row * cell_mm - v["mask_origin"][1]
                abs_pts = geo.translate_points(v["norm_pts"], dx, dy)
                abs_holes = [geo.translate_points(h, dx, dy) for h in v["norm_holes"]]
                sheets.append({
                    "occ": occ,
                    "max_y": row + mh,
                    "max_x": col + mw,
                    "pieces": [{
                        "part_id": inst["part_id"], "name": inst["name"],
                        "angle": ang, "points": abs_pts, "holes": abs_holes,
                        "area": inst["area"],
                    }],
                })
            else:
                unplaced.append(inst)

        if progress_cb:
            progress_cb(idx + 1, total, placed=placed)

    result_sheets = []
    for sheet in sheets:
        used_area = sum(p["area"] for p in sheet["pieces"])
        sheet_area = sheet_w_mm * sheet_h_mm
        result_sheets.append({
            "width": sheet_w_mm,
            "height": sheet_h_mm,
            "pieces": sheet["pieces"],
            "used_area_mm2": used_area,
            "utilization_pct": round(100.0 * used_area / sheet_area, 2) if sheet_area else 0.0,
        })

    return result_sheets, unplaced
