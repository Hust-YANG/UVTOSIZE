#!/usr/bin/env python3
"""
UVTOSIZE Web Server — Flask REST API for quantum dot UV-Vis analysis.
"""

import sys
import os
import uuid
import shutil
import tempfile
from pathlib import Path
from datetime import datetime

# Add the skill scripts directory to Python path
SKILL_DIR = Path(__file__).resolve().parent.parent / ".claude" / "skills" / "UVTOSIZE" / "scripts"
sys.path.insert(0, str(SKILL_DIR))

from flask import Flask, request, jsonify, send_file, make_response
from flask_cors import CORS
from werkzeug.middleware.proxy_fix import ProxyFix

from uv_analysis import (
    analyze_to_dict, get_info,
    detect_qd_type_from_filename, infer_qd_type_from_spectrum,
    parse_uv_txt,
)

app = Flask(__name__, static_folder="static", static_url_path="")
app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1)

# CORS: restrict to production domain in prod, allow all in dev
ALLOWED_ORIGINS = os.environ.get("CORS_ORIGINS", "*")
CORS(app, resources={r"/api/*": {"origins": ALLOWED_ORIGINS.split(",") if ALLOWED_ORIGINS != "*" else "*"}})

# In-memory store for temp directories (token -> path)
TEMP_STORE = {}

# Maximum file size: 50 MB
MAX_FILE_SIZE = 50 * 1024 * 1024

# Auth
AUTH_REQUIRED = os.environ.get("UVTOSIZE_AUTH_REQUIRED", "0") == "1"


def check_auth():
    """Simple auth check."""
    if not AUTH_REQUIRED:
        return True
    token = request.headers.get("Authorization", "").replace("Bearer ", "")
    return bool(token)


def require_auth(f):
    """Decorator for endpoints requiring authentication."""
    from functools import wraps
    @wraps(f)
    def decorated(*args, **kwargs):
        if not check_auth():
            return jsonify({"error": "Authentication required", "code": "AUTH_REQUIRED"}), 401
        return f(*args, **kwargs)
    return decorated


def create_temp_dir():
    """Create a temp directory with a unique token."""
    token = uuid.uuid4().hex[:12]
    tmpdir = tempfile.mkdtemp(prefix=f"uvtosize_{token}_")
    TEMP_STORE[token] = tmpdir
    return token, tmpdir


def cleanup_old_temps():
    """Remove temp directories older than 2 hours."""
    import time
    now = time.time()
    expired = []
    for token, path in list(TEMP_STORE.items()):
        try:
            mtime = os.path.getmtime(path)
            if now - mtime > 7200:
                shutil.rmtree(path, ignore_errors=True)
                expired.append(token)
        except OSError:
            pass
    for token in expired:
        TEMP_STORE.pop(token, None)


def parse_axis_range(request, prefix):
    """Parse xmin/xmax or ymin/ymax from form data."""
    vmin = request.form.get(f"{prefix}min")
    vmax = request.form.get(f"{prefix}max")
    if vmin or vmax:
        return (float(vmin) if vmin else None, float(vmax) if vmax else None)
    return None


@app.route("/api/uvtosize/info", methods=["GET"])
def api_info():
    """Return formula metadata for all supported QD types."""
    try:
        info = get_info()
        return jsonify(info)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/uvtosize/analyze", methods=["POST"])
# @require_auth
def api_analyze():
    """Upload a UV-Vis .txt file and run full analysis."""
    # Validate file
    if "file" not in request.files:
        return jsonify({"error": "No file provided", "code": "NO_FILE"}), 400

    file = request.files["file"]
    if file.filename == "":
        return jsonify({"error": "No file selected", "code": "NO_FILE"}), 400

    if not file.filename.lower().endswith(".txt"):
        return jsonify({"error": "Only .txt files are supported", "code": "INVALID_FORMAT"}), 400

    file_content = file.read()
    if len(file_content) > MAX_FILE_SIZE:
        return jsonify({"error": f"File too large (max {MAX_FILE_SIZE // 1024 // 1024} MB)", "code": "FILE_TOO_LARGE"}), 400
    file.seek(0)

    # ---- Parse parameters ----
    qd_type = request.form.get("qd_type", "auto").strip().lower()
    if qd_type not in ("pbs", "pbse", "cds", "cdse", "auto"):
        return jsonify({"error": f"Invalid QD type: {qd_type}", "code": "INVALID_TYPE"}), 400

    try:
        smooth_window = int(request.form.get("smooth_window", 15))
        if smooth_window < 3 or smooth_window > 201 or smooth_window % 2 == 0:
            smooth_window = 15
    except ValueError:
        smooth_window = 15

    # Axis ranges
    x_range = parse_axis_range(request, "x")
    y_range = parse_axis_range(request, "y")

    # Output formats
    fmt_param = request.form.get("output_formats", "png").strip().lower()
    if fmt_param == "all":
        output_formats = ["png", "pdf", "svg"]
    else:
        output_formats = [f.strip() for f in fmt_param.split(",") if f.strip() in ("png", "pdf", "svg")]
        if not output_formats:
            output_formats = ["png"]

    # DPI
    try:
        dpi = int(request.form.get("dpi", 300))
        if dpi not in (150, 300, 600):
            dpi = 300
    except ValueError:
        dpi = 300

    # Figure width
    figure_width = request.form.get("figure_width", "double").strip().lower()
    if figure_width not in ("single", "double"):
        figure_width = "double"

    # Show annotation box
    show_annotation = request.form.get("show_annotation", "0") == "1"

    # Fit range (optional, auto-selected per QD type if not provided)
    fit_range = None
    fit_range_str = request.form.get("fit_range")
    if fit_range_str:
        try:
            fit_range = float(fit_range_str)
            if fit_range < 5 or fit_range > 300:
                fit_range = None
        except ValueError:
            fit_range = None

    # Baseline mode
    baseline_mode = request.form.get("baseline_mode", "auto").strip().lower()
    if baseline_mode not in ("auto", "constant", "linear", "exponential"):
        baseline_mode = "auto"

    # ---- Process ----
    token, tmpdir = create_temp_dir()

    try:
        safe_filename = file.filename.replace("\\", "_").replace("/", "_")
        tmp_filepath = os.path.join(tmpdir, safe_filename)
        file.save(tmp_filepath)

        # Handle QD type auto-detection
        if qd_type == "auto":
            detected = detect_qd_type_from_filename(file.filename)
            if detected:
                qd_type = detected
            else:
                rankings, _ = infer_qd_type_from_spectrum(tmp_filepath, smooth_window)
                best = rankings[0] if rankings else None
                if best and best["confidence"] in ("high", "medium"):
                    qd_type = best["qd_type"]
                else:
                    return jsonify({
                        "error": "Could not auto-detect QD type",
                        "code": "AUTO_DETECT_FAILED",
                        "rankings": rankings,
                        "message": "Please select the QD type manually."
                    }), 422

        # Run analysis with all parameters
        result = analyze_to_dict(
            tmp_filepath, qd_type, output_dir=tmpdir,
            smooth_window=smooth_window,
            x_range=x_range, y_range=y_range,
            output_formats=output_formats, dpi=dpi,
            figure_width=figure_width,
            show_annotation=show_annotation,
            fit_range_nm=fit_range,
            baseline_mode=baseline_mode,
        )

        # Add token-prefixed download URLs for all formats
        result["plot_url"] = ""
        result["plot_urls"] = {}
        plot_paths = result.get("plot_paths", {})
        for fmt, fpath in plot_paths.items():
            fname = os.path.basename(fpath)
            url = f"/api/uvtosize/download/{token}/{fname}"
            result["plot_urls"][fmt] = url
            if fmt == "png":
                result["plot_url"] = url

        doc_filename = os.path.basename(result.get("doc_path", ""))
        if doc_filename:
            result["docx_url"] = f"/api/uvtosize/download/{token}/{doc_filename}"

        result["token"] = token
        result["qd_type"] = qd_type

        if len(TEMP_STORE) > 50:
            cleanup_old_temps()

        return jsonify(result)

    except Exception as e:
        shutil.rmtree(tmpdir, ignore_errors=True)
        TEMP_STORE.pop(token, None)
        return jsonify({
            "error": f"Analysis failed: {str(e)}",
            "code": "ANALYSIS_ERROR"
        }), 500


@app.route("/api/uvtosize/autodetect", methods=["POST"])
# @require_auth
def api_autodetect():
    """Upload a file and get QD type rankings without full analysis."""
    if "file" not in request.files:
        return jsonify({"error": "No file provided", "code": "NO_FILE"}), 400

    file = request.files["file"]
    if file.filename == "":
        return jsonify({"error": "No file selected", "code": "NO_FILE"}), 400

    file_content = file.read()
    if len(file_content) > MAX_FILE_SIZE:
        return jsonify({"error": f"File too large (max {MAX_FILE_SIZE // 1024 // 1024} MB)", "code": "FILE_TOO_LARGE"}), 400
    file.seek(0)

    token, tmpdir = create_temp_dir()

    try:
        safe_filename = file.filename.replace("\\", "_").replace("/", "_")
        tmp_filepath = os.path.join(tmpdir, safe_filename)
        file.save(tmp_filepath)

        filename_hint = detect_qd_type_from_filename(file.filename)
        rankings, data = infer_qd_type_from_spectrum(tmp_filepath)

        rankings_out = []
        for r in rankings:
            rankings_out.append({
                "qd_type": r["qd_type"],
                "name": r["name"],
                "confidence": r["confidence"],
                "reason": r["reason"],
                "peak_wl_nm": float(r["peak_wl"]),
                "peak_eV": float(r["peak_eV"]),
                "margin_norm": float(r["margin_norm"]),
            })

        return jsonify({
            "filename_hint": filename_hint,
            "rankings": rankings_out,
            "peak_wl_nm": float(rankings[0]["peak_wl"]) if rankings else None,
            "data_points": len(data["wavelength"]) if data else 0,
            "wl_range": [float(data["wavelength"].min()), float(data["wavelength"].max())] if data else None,
        })

    except Exception as e:
        return jsonify({"error": f"Auto-detect failed: {str(e)}", "code": "DETECT_ERROR"}), 500


@app.route("/api/uvtosize/download/<token>/<filename>", methods=["GET"])
# @require_auth
def api_download(token, filename):
    """Download a generated artifact (PNG, PDF, SVG, or DOCX)."""
    tmpdir = TEMP_STORE.get(token)
    if not tmpdir:
        return jsonify({"error": "File not found or expired", "code": "NOT_FOUND"}), 404

    safe_name = os.path.basename(filename)
    filepath = os.path.join(tmpdir, safe_name)

    if not os.path.isfile(filepath):
        return jsonify({"error": "File not found", "code": "NOT_FOUND"}), 404

    # Determine MIME type
    ext = safe_name.lower().rsplit(".", 1)[-1] if "." in safe_name else ""
    mime_map = {"png": "image/png", "pdf": "application/pdf", "svg": "image/svg+xml", "docx": "application/vnd.openxmlformats-officedocument.wordprocessingml.document"}
    return send_file(filepath, as_attachment=True, download_name=safe_name,
                     mimetype=mime_map.get(ext, "application/octet-stream"))


@app.route("/api/uvtosize/health", methods=["GET"])
def api_health():
    """Health check endpoint."""
    return jsonify({
        "status": "ok",
        "timestamp": datetime.now().isoformat(),
        "version": "2.0.0",
    })


# Serve the tool page
@app.route("/")
@app.route("/uvtosize")
def serve_index():
    return send_file(os.path.join(app.static_folder, "uvtosize.html"))


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8765))
    debug = os.environ.get("FLASK_DEBUG", "1") == "1"
    print(f"Starting UVTOSIZE server on http://localhost:{port}")
    print(f"Auth required: {AUTH_REQUIRED}")
    app.run(host="0.0.0.0", port=port, debug=debug)
