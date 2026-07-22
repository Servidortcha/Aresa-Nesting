"""
geometry.py
Utilidades geometricas para el nesting: transformaciones de poligonos y
rasterizacion a grilla (bitmap) para deteccion de colisiones rapida con numpy.

Se eligio un enfoque de grilla (en vez de No-Fit Polygon exacto) porque:
  - Permite rotacion libre (cualquier angulo) sin geometria analitica compleja.
  - Es mucho mas simple de implementar correctamente sin librerias externas
    pesadas de geometria (no hay shapely/pyclipper disponibles, no hay acceso
    a internet para instalarlas).
  - Es suficientemente preciso para corte laser ajustando la resolucion de
    la grilla (configurable).

El test punto-en-poligono se hace con `matplotlib.path.Path.contains_points`,
vectorizado en C sobre toda la grilla de una vez (mucho mas rapido que un
scanline linea por linea en Python puro, que era el cuello de botella
original con piezas de muchos puntos como splines finas).

Contrapartida: piezas muy grandes o resoluciones muy finas consumen mas
memoria/tiempo. Se documenta como limitacion conocida.
"""
import math
import numpy as np
from matplotlib.path import Path


def polygon_area(points):
    a = 0.0
    n = len(points)
    for i in range(n):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def polygon_bounds(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return min(xs), min(ys), max(xs), max(ys)


def rotate_points(points, angle_deg, origin=(0, 0)):
    a = math.radians(angle_deg)
    ca, sa = math.cos(a), math.sin(a)
    ox, oy = origin
    out = []
    for x, y in points:
        x -= ox
        y -= oy
        out.append((x * ca - y * sa + ox, x * sa + y * ca + oy))
    return out


def translate_points(points, dx, dy):
    return [(x + dx, y + dy) for x, y in points]


def normalize_to_origin(points):
    """Traslada el poligono para que su bounding box empiece en (0,0)."""
    minx, miny, _, _ = polygon_bounds(points)
    return translate_points(points, -minx, -miny), (minx, miny)


def rasterize_polygon(points, holes, cell_mm, margin_cells=1):
    """
    Devuelve (mask_bool_2d, (w_cells, h_cells)) donde True = ocupado por la
    pieza (incluye el margen/espaciado ya aplicado por fuera si se agranda
    el poligono antes de llamar a esta funcion).
    Usa el algoritmo de punto-en-poligono por scanline, vectorizado con numpy.
    """
    minx, miny, maxx, maxy = polygon_bounds(points)
    w = int(math.ceil((maxx - minx) / cell_mm)) + 2 * margin_cells
    h = int(math.ceil((maxy - miny) / cell_mm)) + 2 * margin_cells
    w = max(w, 1)
    h = max(h, 1)

    mask = np.zeros((h, w), dtype=bool)
    # coordenadas del centro de cada celda en espacio real
    xs = (np.arange(w) - margin_cells + 0.5) * cell_mm + minx
    ys = (np.arange(h) - margin_cells + 0.5) * cell_mm + miny
    grid_x, grid_y = np.meshgrid(xs, ys)
    grid_pts = np.column_stack([grid_x.ravel(), grid_y.ravel()])

    def fill(poly, value):
        if len(poly) < 3:
            return
        path = Path(poly)
        inside = path.contains_points(grid_pts, radius=1e-9).reshape(h, w)
        mask[inside] = value

    fill(points, True)
    for hole in holes or []:
        fill(hole, False)

    return mask, (minx - margin_cells * cell_mm, miny - margin_cells * cell_mm)


def dilate_mask(mask, cells):
    """Dilata la mascara `cells` celdas (aprox. circular) para representar el
    espaciado/kerf entre piezas. Se implementa como una convolucion con un
    kernel circular via FFT (scipy), mucho mas rapido que iterar offset por
    offset en Python puro cuando `cells` crece."""
    if cells <= 0:
        return mask
    yy, xx = np.ogrid[-cells:cells + 1, -cells:cells + 1]
    kernel = (xx * xx + yy * yy) <= cells * cells
    from scipy.signal import fftconvolve
    conv = fftconvolve(mask.astype(np.float32), kernel.astype(np.float32), mode="same")
    return conv > 0.5
