import os
import io
import json
import uuid
import zipfile
import traceback

from flask import Flask, request, render_template, jsonify, send_file, abort

import dxf_io
import nesting
import geometry as geo

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
UPLOAD_DIR = os.path.join(BASE_DIR, "uploads")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
os.makedirs(UPLOAD_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024  # 50 MB


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/optimize", methods=["POST"])
def optimize():
    try:
        sheet_w = float(request.form.get("sheet_w"))
    except (TypeError, ValueError):
        return jsonify({"error": "Falta o es inválido el ancho de chapa."}), 400

    sheet_h_raw = (request.form.get("sheet_h") or "").strip()
    max_length_mm = None
    if sheet_h_raw:
        try:
            v = float(sheet_h_raw)
            if v > 0:
                max_length_mm = v
        except ValueError:
            return jsonify({"error": "El largo máximo de chapa es inválido."}), 400

    try:
        spacing = float(request.form.get("spacing", 3))
    except (TypeError, ValueError):
        spacing = 3.0
    allow_rotation = request.form.get("allow_rotation", "1") == "1"

    # Parametros avanzados opcionales: si no vienen, se calculan solos
    # (resolucion de grilla y paso de rotacion) segun el tamano/cantidad
    # de piezas -- el usuario no necesita elegir "a mano" los grados.
    rotation_step = request.form.get("rotation_step", "").strip()
    cell_mm = request.form.get("cell_mm", "").strip()
    rotation_step = float(rotation_step) if rotation_step else None
    cell_mm = float(cell_mm) if cell_mm else None

    files = request.files.getlist("dxf_files")
    if not files:
        return jsonify({"error": "No se subió ningún archivo DXF."}), 400

    qtys_raw = request.form.get("quantities", "{}")
    try:
        qtys = json.loads(qtys_raw)
    except json.JSONDecodeError:
        qtys = {}

    job_id = uuid.uuid4().hex[:12]
    job_upload_dir = os.path.join(UPLOAD_DIR, job_id)
    os.makedirs(job_upload_dir, exist_ok=True)

    parts = []
    errors = []
    for f in files:
        fname = f.filename
        if not fname.lower().endswith(".dxf"):
            errors.append(f"'{fname}' no es un .dxf, se omitió.")
            continue
        save_path = os.path.join(job_upload_dir, fname)
        f.save(save_path)
        qty = int(qtys.get(fname, 1) or 1)
        try:
            outer, holes, warning = dxf_io.load_dxf_piece(save_path)
            if warning:
                errors.append(f"'{fname}': {warning}")
        except Exception as e:
            errors.append(f"No se pudo leer '{fname}': {e}")
            continue
        part_id = uuid.uuid4().hex[:8]
        parts.append({
            "part_id": part_id,
            "name": fname,
            "points": outer,
            "holes": holes,
            "qty": qty,
        })

    if not parts:
        return jsonify({"error": "Ningún archivo pudo procesarse.", "details": errors}), 400

    try:
        sheets = []
        remaining = parts
        impossible_names = set()
        for _ in range(20):  # tope de seguridad de chapas por corrida
            if not remaining:
                break
            sheet, impossible = nesting.nest_strip(
                remaining, sheet_w, spacing_mm=spacing,
                cell_mm=cell_mm, allow_rotation=allow_rotation,
                rotation_step_deg=rotation_step, max_length_mm=max_length_mm,
            )
            for imp in impossible:
                impossible_names.add(imp["name"])
            if sheet is None or not sheet["pieces"]:
                break
            sheets.append(sheet)

            placed_count = {}
            for p in sheet["pieces"]:
                placed_count[p["part_id"]] = placed_count.get(p["part_id"], 0) + 1
            next_remaining = []
            for p in remaining:
                used = placed_count.get(p["part_id"], 0)
                left_qty = p["qty"] - used
                if left_qty > 0 and p["name"] not in impossible_names:
                    next_remaining.append({**p, "qty": left_qty})
            remaining = next_remaining
            if max_length_mm is None:
                break  # sin tope de largo, una sola chapa alcanza siempre
        unplaced = [{"name": n} for n in impossible_names]
        for p in remaining:
            if p["name"] not in impossible_names:
                unplaced.append({"name": p["name"]})
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Error durante el nesting: {e}"}), 500

    job_output_dir = os.path.join(OUTPUT_DIR, job_id)
    os.makedirs(job_output_dir, exist_ok=True)
    sheet_files = []
    for i, sheet in enumerate(sheets):
        out_path = os.path.join(job_output_dir, f"chapa_{i + 1}.dxf")
        dxf_io.save_dxf_layout(out_path, sheet["pieces"], sheet["width"], sheet["height"])
        sheet_files.append(f"chapa_{i + 1}.dxf")

    response_sheets = []
    for i, sheet in enumerate(sheets):
        response_sheets.append({
            "index": i + 1,
            "width": sheet["width"],
            "height": sheet["height"],
            "utilization_pct": sheet["utilization_pct"],
            "file": sheet_files[i],
            "pieces": [
                {
                    "name": p["name"],
                    "angle": round(p["angle"], 1),
                    "points": p["points"],
                    "holes": p["holes"],
                } for p in sheet["pieces"]
            ],
        })

    return jsonify({
        "job_id": job_id,
        "sheets": response_sheets,
        "unplaced": unplaced,
        "warnings": errors,
    })


@app.route("/download/<job_id>/<filename>")
def download(job_id, filename):
    path = os.path.join(OUTPUT_DIR, job_id, filename)
    if not os.path.isfile(path):
        abort(404)
    return send_file(path, as_attachment=True, download_name=filename)


@app.route("/download/<job_id>/all")
def download_all(job_id):
    job_output_dir = os.path.join(OUTPUT_DIR, job_id)
    if not os.path.isdir(job_output_dir):
        abort(404)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for fname in sorted(os.listdir(job_output_dir)):
            zf.write(os.path.join(job_output_dir, fname), arcname=fname)
    buf.seek(0)
    return send_file(buf, as_attachment=True, download_name=f"layout_{job_id}.zip", mimetype="application/zip")


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
