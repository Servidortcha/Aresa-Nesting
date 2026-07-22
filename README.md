# Optimización de cortes (nesting) para láser

Web app que recibe archivos `.dxf` (una pieza por archivo), la cantidad de
copias de cada una, y arma el mejor acomodo posible sobre una o más chapas,
exportando el resultado como un nuevo `.dxf` listo para cortar.

## Instalación

```bash
cd nesting_app
python3 -m venv venv
source venv/bin/activate        # en Windows: venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

Abrí `http://localhost:5000` en el navegador.

## Cómo se usa

1. Arrastrá o seleccioná uno o varios `.dxf` (cada uno con **una sola pieza**).
2. Indicá cuántas copias de cada pieza necesitás.
3. Definí el **ancho de chapa** (el ancho fijo del material que comprás,
   en mm). El **largo lo calcula la app automáticamente** — es el mínimo
   necesario para que entren todas las piezas, así no desperdiciás
   material. Si tenés un largo máximo real (por ejemplo, el tamaño de la
   chapa que ya tenés), completalo en "Largo máximo"; si hace falta más
   espacio del que entra en esa chapa, la app arma varias chapas.
4. Definí el espaciado/kerf entre piezas (en mm) — configurable por corte.
5. Dejá tildado "Permitir rotación": el **ángulo de cada pieza se calcula
   solo**, probando varias orientaciones y quedándose con la que menos
   espacio ocupa — no hace falta elegir grados a mano.
6. Click en "Optimizar layout". Vas a ver una vista previa por chapa con el
   % de aprovechamiento de material, y un botón para descargar el `.dxf`
   resultante (o todas las chapas juntas en un `.zip` si hizo falta más de
   una).

Hay una sección "Opciones avanzadas" opcional (paso de rotación en grados y
resolución de grilla) por si querés forzar esos valores manualmente; dejarla
vacía usa el cálculo automático, que es lo recomendado.

## Cómo funciona el algoritmo (nesting.py)

No se usa fuerza bruta ni una librería de nesting comercial: es un
algoritmo propio pensado para andar sin dependencias pesadas (no hay acceso
a internet para instalar librerías como `ezdxf`, `shapely` o `pyclipper`).

1. Cada pieza se **rasteriza** a una grilla/bitmap (la resolución en mm por
   celda es configurable — más fino = más preciso pero más lento).
2. El espaciado pedido se aplica **dilatando** la máscara de cada pieza la
   mitad del espaciado, así cuando dos piezas "se tocan" en la grilla, la
   distancia real entre sus contornos es la pedida.
3. Para la rotación libre, se prueban ángulos discretos (paso configurable,
   por defecto cada 15°) — no es rotación 100% continua, pero con un paso
   fino se acerca bastante, ajustable según cuánto tiempo de cálculo estés
   dispuesto a esperar.
4. La búsqueda de la posición sin colisión se hace con **correlación 2D vía
   FFT** (`scipy.signal.fftconvolve`), mucho más rápida que probar
   posición por posición con loops.
5. Orden de acomodo: piezas de mayor área primero (heurística estándar en
   nesting — las piezas grandes son las que más condicionan el resultado).
6. Entre todas las posiciones y ángulos válidos para cada pieza, se elige
   la que **menos hace crecer el área total ocupada** hasta el momento
   (no solo "la primera que entra") — así el layout final tiende a quedar
   compacto en vez de esparcido con ángulos raros que desperdician espacio.
7. El ancho de chapa es fijo (el material que comprás) y el **largo se
   calcula automáticamente**: se prueba con una estimación inicial y se
   va ajustando hasta encontrar el mínimo largo que hace entrar todas las
   piezas. Si definiste un largo máximo, al llegarlo se abre una chapa
   nueva con las piezas que quedaron afuera.
8. Si una pieza no entra en ninguna chapa en ningún ángulo (más ancha que
   el ancho de chapa), queda reportada como "sin ubicar".

## Formato DXF: soporte y limitaciones

`dxf_io.py` es un lector/escritor propio (no usa librerías externas).

**Lectura soportada:**
- `LWPOLYLINE` (incluye arcos vía "bulge")
- `POLYLINE` / `VERTEX` (formato clásico)
- `LINE`, `CIRCLE`, `ARC`
- `SPLINE` — se evalúa la curva real (algoritmo de De Boor, con vector de
  nudos y pesos si es racional), no solo sus puntos de control.
- **Contornos partidos en varias entidades sueltas**: muchos CAD (ej.
  exports que vimos con LINE + SPLINE por separado, sin agrupar en una
  polilínea) no arman el contorno como una sola entidad. La app detecta
  automáticamente los tramos abiertos (LINE/SPLINE/ARC) y los **encadena
  por sus extremos** hasta formar los contornos cerrados de la pieza.
  La tolerancia de unión es de 0.01 mm por defecto (`join_tolerance` en
  `dxf_io.load_dxf_piece`); si tu CAD deja micro-huecos más grandes entre
  tramos, avisame para ajustarla.

**No soportado (limitación conocida):**
- `BLOCK` / `INSERT` (piezas armadas como bloques anidados) — exportá la
  pieza "explotada" (sin bloques) desde tu programa de CAD.
- Múltiples piezas dentro de un mismo archivo DXF: se asume **un DXF = una
  pieza**. Si tu archivo tiene varios contornos cerrados, se toma el de
  mayor área como pieza principal y el resto como agujeros internos (así
  es como se detectan automáticamente ranuras/recortes interiores).

**Escritura:** genera un DXF válido (formato R2000) con capas `PIEZAS`
(las piezas ubicadas) y `CHAPA` (el contorno de la chapa, de referencia).
Se probó abriendo en lectores DXF estándar; si tu software específico
(LightBurn, AutoCAD, etc.) tiene algún requisito particular de layer/units,
avisame y lo ajustamos.

## Estructura del proyecto

```
nesting_app/
├── app.py            # servidor Flask (rutas web + orquestación)
├── dxf_io.py         # lectura y escritura de DXF (sin dependencias)
├── geometry.py        # transformaciones y rasterización de polígonos
├── nesting.py         # algoritmo de nesting (FFT + heurística bottom-left)
├── templates/
│   └── index.html     # interfaz web (subida, parámetros, preview)
├── requirements.txt
└── uploads/ outputs/   # carpetas de trabajo (se crean archivos temporales)
```

## Posibles mejoras futuras (no incluidas en esta versión)

- Reconocer y **empaquetar mejor formas cóncavas complejas** anidando piezas
  chicas dentro de huecos de piezas grandes (nesting con agujeros
  aprovechables).
- Paso de rotación adaptativo (más fino solo cuando hace falta) para mejorar
  velocidad sin perder aprovechamiento.
- Guardar/recordar configuraciones de chapa y espaciado más usadas.
- Autenticación / multiusuario si esto se va a compartir en red.
