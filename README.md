# UVTOSIZE — Quantum Dot Size & Size Distribution from UV-Vis Spectra

UV-Vis absorption spectroscopy analysis tool for semiconductor quantum dots (PbS, PbSe, CdS, CdSe). Automatically calculates nanoparticle diameter and size distribution from the first exciton peak.

## Features

- **Peak Detection** — Auto-identifies first exciton peak (1Se–1Sh) via Savitzky-Golay smoothing + scipy.signal.find_peaks
- **QD Type Auto-Detection** — Three-tier strategy: filename keywords → spectral inference → manual selection
- **Baseline Subtraction** — Three models (constant, linear, exponential) with AIC-based auto-selection
- **Two-Step Baseline** — Key innovation: baseline fitted from peak-excluded data, then pure Gaussian on subtracted spectrum (prevents CdS/CdSe absorption-edge overfitting)
- **Literature HWHM Constraints** — Prevents unphysical Gaussian fits using published HWHM ranges per QD type
- **Size Distribution** — Error propagation: σ_d = HWHM(E) / |dE/dd| (Nikolaev & Averkiev, 2009)
- **Publication-Quality Plots** — Wong (2011) colorblind-friendly palette, Times New Roman, inward ticks, multi-format export (PNG/PDF/SVG)
- **Word Report** — Auto-generated .docx with formulas, step-by-step calculations, and references
- **Web Interface** — Flask REST API + Plotly.js interactive charts


## Supported QD Materials

| Material | Sizing Formula | Reference | Valid Range |
|----------|---------------|-----------|-------------|
| PbS | E₀ = 0.41 + 1/(0.0252·d² + 0.283·d) | Moreels et al., *ACS Nano* 2009 | 2.3–10 nm |
| PbSe | E₀ = 0.278 + 1/(0.016·d² + 0.209·d + 0.45) | Moreels et al., *Chem. Mater.* 2007 | 2–20 nm |
| CdS | D = poly(λ) | Yu et al., *Chem. Mater.* 2003 | 1–6 nm |
| CdSe | D = poly(λ) | Yu et al., *Chem. Mater.* 2003 | 1.5–8 nm |

## Quick Start

### CLI
```bash
pip install numpy scipy matplotlib python-docx
python .claude/skills/UVTOSIZE/scripts/uv_analysis.py PbS.txt --type pbs --output results_PbS
```



## Project Structure

```
UVTOSIZE/
├── .claude/skills/UVTOSIZE/   # Claude Code skill definition
│   ├── SKILL.md                # Skill metadata & documentation
│   └── scripts/
│       └── uv_analysis.py      # Core analysis engine (~2400 lines)
├── web/                        # Flask web application
│   ├── server.py               # REST API (5 endpoints)
│   ├── static/uvtosize.html    # Plotly.js frontend
│   ├── wsgi.py                 # Gunicorn entry point
│   └── deploy/                 # Nginx, systemd, deploy script
├── mods-uvtosize/              # MODS tool platform package
│   ├── algo.py                 # process() interface
│   ├── app.py                  # Flask framework
│   ├── templates/tool.html     # Frontend (hust-mods.com theme)
│   └── sample/PbS.txt          # Test sample
├── PbS.txt                     # Sample data (PbS, d=6.08 nm)
├── CdS.txt                     # Sample data (CdS, d=4.53 nm)

```

## Key Design Decisions

### Two-Step Baseline Subtraction
Simultaneous fitting of Gaussian + exponential baseline fails for CdS/CdSe because the baseline is "squeezed" to near-zero at the peak to allow the Gaussian to match peak height. This causes severe left-side underfit (baseline cannot simultaneously capture the high absorption edge and stay low at the peak).

**Solution**: Fit baseline to peak-excluded data first → subtract → fit pure Gaussian. Improves CdS R² from 0.94 to 0.975.

### AIC-Based Model Selection with Physical Penalties
Three baseline models are evaluated, and AIC selects the best one. Physical penalties prevent unrealistic solutions:
- **Peak fraction penalty**: If Gaussian amplitude / total at peak < 0.3 → penalize
- **Negative baseline penalty**: Absorbance must be ≥ 0

### HWHM Constraints
Literature-derived HWHM ranges (10–350 meV for CdS, 10–200 meV for others) constrain the Gaussian fit bounds, preventing grossly unphysical results.

## References

- Moreels, I. et al. Size-Tunable PbS QDs. *ACS Nano* **2009**, 3, 3023–3030.
- Moreels, I. et al. Composition and Size-Dependent Extinction Coefficient of Colloidal PbSe QDs. *Chem. Mater.* **2007**, 19, 6101–6106.
- Yu, W. W. et al. Experimental Determination of the Extinction Coefficient of CdTe, CdSe, and CdS Nanocrystals. *Chem. Mater.* **2003**, 15, 2854–2860.
- Nikolaev, V. V. & Averkiev, N. S. Size distribution of QDs from absorption spectra. *Appl. Phys. Lett.* **2009**, 95, 263107.
- Wong, B. Points of view: Color blindness. *Nature Methods* **2011**, 8, 441.




