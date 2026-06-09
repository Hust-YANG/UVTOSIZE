"""
MODS Flask framework - UVTOSIZE edition.
Adds /api/analyze (JSON) and /api/report (DOCX) alongside standard /api/run (PNG).
"""
import os
import io
import tempfile
import traceback

from flask import Flask, request, send_file, jsonify, render_template

import algo

app = Flask(__name__)
MAX_UPLOAD_BYTES = 200 * 1024 * 1024


@app.route("/")
def index():
    return render_template("tool.html")


@app.route("/api/run", methods=["POST"])
def run():
    """Standard MODS: upload file -> PNG plot."""
    f = request.files.get("file")
    if not f or not f.filename:
        return jsonify(error="No file uploaded"), 400

    fd, tmp_path = tempfile.mkstemp(suffix=os.path.splitext(f.filename)[1])
    os.close(fd)
    f.save(tmp_path)
    try:
        out_bytes, out_name, mime = algo.process(tmp_path, dict(request.form))
    except Exception as e:
        traceback.print_exc()
        return jsonify(error=str(e)), 500
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass

    return send_file(io.BytesIO(out_bytes), mimetype=mime,
                     as_attachment=True, download_name=out_name)


@app.route("/api/analyze", methods=["POST"])
def analyze():
    """Upload file -> full analysis JSON with plot data arrays."""
    f = request.files.get("file")
    if not f or not f.filename:
        return jsonify(error="No file uploaded"), 400

    fd, tmp_path = tempfile.mkstemp(suffix=os.path.splitext(f.filename)[1])
    os.close(fd)
    f.save(tmp_path)
    try:
        result = algo.process_json(tmp_path, dict(request.form))
    except Exception as e:
        traceback.print_exc()
        return jsonify(error=str(e)), 500
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass

    return jsonify(result)


@app.route("/api/report", methods=["POST"])
def report():
    """Upload file -> DOCX Word report."""
    f = request.files.get("file")
    if not f or not f.filename:
        return jsonify(error="No file uploaded"), 400

    fd, tmp_path = tempfile.mkstemp(suffix=os.path.splitext(f.filename)[1])
    os.close(fd)
    f.save(tmp_path)
    try:
        out_bytes, out_name, mime = algo.process_report(tmp_path, dict(request.form))
    except Exception as e:
        traceback.print_exc()
        return jsonify(error=str(e)), 500
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass

    return send_file(io.BytesIO(out_bytes), mimetype=mime,
                     as_attachment=True, download_name=out_name)
