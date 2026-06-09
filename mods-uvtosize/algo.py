"""
UVTOSIZE - Quantum Dot Size & Size Distribution from UV-Vis Spectra (MODS Tool)

Standard MODS interface:
    process(input_path, params) -> (output_bytes, output_filename, mimetype)

Extended interfaces (used by app.py for rich frontend):
    process_json(input_path, params) -> dict  (all results + plot data)
    process_report(input_path, params) -> (docx_bytes, filename, mimetype)
"""
import io
import sys
import os
import json
import tempfile
import zipfile
import numpy as np

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from uv_analysis import (
    parse_uv_txt,
    find_first_exciton_peak,
    calculate_size,
    fit_gaussian_to_peak,
    calculate_size_distribution,
    create_uv_plot,
    detect_qd_type_from_filename,
    infer_qd_type_from_spectrum,
    SIZE_FORMULAS,
    HWHM_CONSTRAINTS,
    DEFAULT_FIT_RANGE,
)


def _parse_params(params):
    """Parse and validate all form parameters. Returns a clean dict."""
    p = params or {}
    clean = {}

    # qd_type
    qd_type = p.get("qd_type", "auto").strip().lower()
    clean["qd_type"] = qd_type if qd_type in ("pbs", "pbse", "cds", "cdse", "auto") else "auto"

    # smooth_window
    try:
        sw = int(p.get("smooth_window", 15))
        clean["smooth_window"] = sw if (3 <= sw <= 201 and sw % 2 == 1) else 15
    except (ValueError, TypeError):
        clean["smooth_window"] = 15

    # fit_range
    fr = p.get("fit_range", "").strip()
    clean["fit_range"] = None
    if fr:
        try:
            v = float(fr)
            if 5 <= v <= 300:
                clean["fit_range"] = v
        except (ValueError, TypeError):
            pass

    # baseline_mode
    bm = p.get("baseline_mode", "auto").strip().lower()
    clean["baseline_mode"] = bm if bm in ("auto", "constant", "linear", "exponential") else "auto"

    # x_range / y_range
    for ax in ("x", "y"):
        clean[f"{ax}min"] = None
        clean[f"{ax}max"] = None
        vmin = p.get(f"{ax}min", "").strip()
        vmax = p.get(f"{ax}max", "").strip()
        if vmin:
            try:
                clean[f"{ax}min"] = float(vmin)
            except (ValueError, TypeError):
                pass
        if vmax:
            try:
                clean[f"{ax}max"] = float(vmax)
            except (ValueError, TypeError):
                pass

    # output_formats
    fmt = p.get("output_formats", "png").strip().lower()
    if fmt == "all":
        clean["output_formats"] = ["png", "pdf", "svg"]
    else:
        fmts = [f.strip() for f in fmt.split(",") if f.strip() in ("png", "pdf", "svg")]
        clean["output_formats"] = fmts if fmts else ["png"]

    # dpi
    try:
        dpi = int(p.get("dpi", 300))
        clean["dpi"] = dpi if dpi in (150, 300, 600) else 300
    except (ValueError, TypeError):
        clean["dpi"] = 300

    # figure_width
    fw = p.get("figure_width", "double").strip().lower()
    clean["figure_width"] = fw if fw in ("single", "double") else "double"

    # show_annotation
    clean["show_annotation"] = p.get("show_annotation", "0") == "1"

    return clean


def _run_analysis(input_path, params):
    """Core analysis: parse, detect, fit, calculate. Returns all results dict."""
    c = _parse_params(params)

    data = parse_uv_txt(input_path)
    if data is None or len(data.get("wavelength", [])) < 10:
        raise ValueError("Could not parse the UV-Vis file. Check format (need 2-column wavelength/absorbance data).")

    # QD type auto-detection
    qd_type = c["qd_type"]
    rankings = None
    if qd_type == "auto":
        detected = detect_qd_type_from_filename(os.path.basename(input_path))
        if detected:
            qd_type = detected
        else:
            rankings, _ = infer_qd_type_from_spectrum(input_path, c["smooth_window"])
            best = rankings[0] if rankings else None
            if best and best["confidence"] in ("high", "medium"):
                qd_type = best["qd_type"]
            else:
                raise ValueError(
                    "Could not auto-detect QD type. "
                    "Please select manually: PbS, PbSe, CdS, or CdSe."
                )

    # Peak detection
    peak_info = find_first_exciton_peak(
        data["wavelength"], data["absorbance"], qd_type, c["smooth_window"]
    )
    if peak_info is None:
        raise ValueError("Could not detect the first exciton peak.")

    # Size calculation
    size_nm, E0, warnings_list = calculate_size(qd_type, peak_info["peak_wavelength"])

    # Gaussian fit
    peak_idx = int(abs(data["wavelength"] - peak_info["peak_wavelength"]).argmin())
    gfit = fit_gaussian_to_peak(
        data["wavelength"], data["absorbance"], peak_idx,
        qd_type=qd_type, fit_range_nm=c["fit_range"], baseline_mode=c["baseline_mode"],
    )

    # Size distribution
    sigma_info = None
    if gfit["fit_success"]:
        sigma_info = calculate_size_distribution(
            qd_type, size_nm, gfit["hwhm_eV"], peak_info["peak_wavelength"]
        )

    # Build x_range / y_range tuples
    x_range = None
    if c["xmin"] is not None or c["xmax"] is not None:
        x_range = (c["xmin"], c["xmax"])
    y_range = None
    if c["ymin"] is not None or c["ymax"] is not None:
        y_range = (c["ymin"], c["ymax"])

    return {
        "data": data,
        "peak_info": peak_info,
        "qd_type": qd_type,
        "size_nm": size_nm,
        "gfit": gfit,
        "sigma_info": sigma_info,
        "warnings": warnings_list,
        "rankings": rankings,
        "params": c,
        "x_range": x_range,
        "y_range": y_range,
    }


def process(input_path, params=None):
    """Standard MODS interface: returns PNG plot bytes."""
    r = _run_analysis(input_path, params)
    c = r["params"]

    output_dir = tempfile.mkdtemp(prefix="uvtosize_")
    plot_paths = create_uv_plot(
        r["data"], r["peak_info"], r["qd_type"], r["size_nm"], output_dir,
        gaussian_fit=r["gfit"], sigma_info=r["sigma_info"],
        output_formats=["png"], dpi=c["dpi"],
        figure_width=c["figure_width"],
        show_annotation=c["show_annotation"],
        x_range=r["x_range"], y_range=r["y_range"],
    )

    png_path = plot_paths.get("png")
    if not png_path or not os.path.exists(png_path):
        raise RuntimeError("Failed to generate plot.")

    with open(png_path, "rb") as f:
        out = f.read()

    base = os.path.splitext(os.path.basename(input_path))[0]
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in base)
    name = "UVTOSIZE_{0}_{1}_d{2:.2f}nm.png".format(safe, r["qd_type"].upper(), r["size_nm"])
    return out, name, "image/png"


def process_json(input_path, params=None):
    """Extended interface: returns dict with all results + plot data arrays for charts."""
    r = _run_analysis(input_path, params)
    c = r["params"]
    gfit = r["gfit"]
    sigma = r["sigma_info"]
    peak = r["peak_info"]
    data = r["data"]
    formula_info = SIZE_FORMULAS.get(r["qd_type"].lower(), {})
    constraints = HWHM_CONSTRAINTS.get(r["qd_type"].lower(), {})

    # Build fit curve arrays for Plotly
    fit_wl = []
    fit_abs = []
    gauss_only_wl = []
    gauss_only_abs = []
    if gfit["fit_success"] and gfit.get("fit_wavelength") is not None:
        fit_wl = gfit["fit_wavelength"].tolist()
        fit_abs = gfit["fit_absorbance"].tolist()
        if gfit.get("gauss_only") is not None:
            gauss_only_wl = gfit["fit_wavelength"].tolist()
            gauss_only_abs = gfit["gauss_only"].tolist()

    # HWHM markers
    peak_wl = float(peak["peak_wavelength"])
    hwhm_left = None
    hwhm_right = None
    hwhm_half_max = None
    if gfit["fit_success"] and not np.isnan(gfit.get("hwhm_nm", float("nan"))):
        hwhm_left = peak_wl - float(gfit["hwhm_nm"])
        hwhm_right = peak_wl + float(gfit["hwhm_nm"])
        if gfit["fit_params"].get("half_max") is not None:
            hwhm_half_max = float(gfit["fit_params"]["half_max"])
        elif gfit["fit_params"]["A"] is not None:
            hwhm_half_max = float(gfit["fit_params"]["A"]) / 2.0

    result = {
        "success": True,
        "sample_name": data.get("sample_name", os.path.basename(input_path)),
        "qd_type": r["qd_type"],
        "qd_name": formula_info.get("name", r["qd_type"].upper()),
        "peak_wavelength_nm": peak_wl,
        "peak_energy_eV": float(peak["peak_energy_eV"]),
        "peak_absorbance": float(peak["peak_absorbance"]),
        "size_nm": float(r["size_nm"]),
        "formula_name": formula_info.get("reference", ""),
        "formula_str": formula_info.get("formula_str", ""),
        "valid_range_wl": constraints.get("wl_range", ""),
        "valid_range_size": "{0}-{1} nm".format(
            formula_info.get("size_min", "?"), formula_info.get("size_max", "?")
        ),
        # Gaussian fit
        "fit_success": gfit["fit_success"],
        "fit_r2": float(gfit["fit_r2"]) if gfit["fit_r2"] is not None else None,
        "hwhm_nm": float(gfit["hwhm_nm"]) if not np.isnan(gfit.get("hwhm_nm", float("nan"))) else None,
        "hwhm_eV": float(gfit["hwhm_eV"]) if not np.isnan(gfit.get("hwhm_eV", float("nan"))) else None,
        "fwhm_nm": float(gfit["fwhm_nm"]) if not np.isnan(gfit.get("fwhm_nm", float("nan"))) else None,
        "gaussian_sigma_nm": float(gfit["gaussian_sigma_nm"]) if not np.isnan(gfit.get("gaussian_sigma_nm", float("nan"))) else None,
        "fit_amplitude": float(gfit["fit_params"]["A"]) if gfit["fit_success"] else None,
        "fit_center": float(gfit["fit_params"]["x0"]) if gfit["fit_success"] else None,
        "baseline_at_peak": float(gfit["fit_params"]["baseline"]) if gfit["fit_success"] else None,
        "baseline_model": gfit.get("baseline_model"),
        "aic": float(gfit.get("aic", float("nan"))) if gfit.get("aic") is not None else None,
        # Size distribution
        "sigma_nm": float(sigma["sigma_nm"]) if sigma else None,
        "relative_sigma_percent": float(sigma["relative_sigma_percent"]) if sigma else None,
        "deriv_dE_dd": float(sigma["deriv_dE_dd"]) if sigma else None,
        # Warnings
        "warnings": r["warnings"],
        "hwhm_warning": gfit.get("hwhm_warning"),
        # Plot data arrays
        "wavelength": data["wavelength"].tolist(),
        "absorbance": data["absorbance"].tolist(),
        "fit_wl": fit_wl,
        "fit_abs": fit_abs,
        "gauss_only_wl": gauss_only_wl,
        "gauss_only_abs": gauss_only_abs,
        "hwhm_left_wl": hwhm_left,
        "hwhm_right_wl": hwhm_right,
        "hwhm_half_max_abs": hwhm_half_max,
        # Auto-detection rankings
        "rankings": None,
    }

    if r["rankings"]:
        result["rankings"] = []
        for rk in r["rankings"]:
            result["rankings"].append({
                "qd_type": rk["qd_type"],
                "name": rk["name"],
                "confidence": rk["confidence"],
                "reason": rk["reason"],
                "peak_wl_nm": float(rk["peak_wl"]),
                "peak_eV": float(rk["peak_eV"]),
                "margin_norm": float(rk["margin_norm"]),
            })

    return result


def process_report(input_path, params=None):
    """Extended interface: returns DOCX Word report bytes."""
    r = _run_analysis(input_path, params)
    c = r["params"]

    output_dir = tempfile.mkdtemp(prefix="uvtosize_")
    plot_paths = create_uv_plot(
        r["data"], r["peak_info"], r["qd_type"], r["size_nm"], output_dir,
        gaussian_fit=r["gfit"], sigma_info=r["sigma_info"],
        output_formats=["png"], dpi=c["dpi"],
        figure_width=c["figure_width"],
        show_annotation=c["show_annotation"],
        x_range=r["x_range"], y_range=r["y_range"],
    )

    png_path = plot_paths.get("png", "")
    from uv_analysis import create_word_report
    docx_path = create_word_report(
        r["data"], r["peak_info"], r["qd_type"], r["size_nm"],
        png_path, output_dir,
        warnings=r["warnings"], gaussian_fit=r["gfit"], sigma_info=r["sigma_info"],
    )

    if not docx_path or not os.path.exists(docx_path):
        raise RuntimeError("Failed to generate Word report.")

    with open(docx_path, "rb") as f:
        out = f.read()

    base = os.path.splitext(os.path.basename(input_path))[0]
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in base)
    name = "UVTOSIZE_Report_{0}_{1}.docx".format(safe, r["qd_type"].upper())
    return out, name, "application/vnd.openxmlformats-officedocument.wordprocessingml.document"
