---
name: UVTOSIZE
description: >
  Analyze UV-Vis absorption spectra of quantum dots (PbS, PbSe, CdS, CdSe) to
  determine nanoparticle diameter and size distribution. Use this skill whenever
  the user mentions UV-Vis data, quantum dot size calculation, QD absorption
  spectra, UV peak analysis, size distribution (sigma), polydispersity, or wants
  to process .txt files from a UV-Vis spectrophotometer for quantum dots. Also
  trigger when the user mentions "UVTOSIZE", "量子点尺寸", "紫外吸收光谱", "UV图",
  "粒径分布", or asks to calculate nanoparticle diameter/polydispersity from
  optical spectra.
---

# UVTOSIZE — Quantum Dot Size & Size Distribution from UV-Vis Spectra

## Overview

This skill analyzes UV-Vis absorption spectroscopy data (`.txt` files) for
semiconductor quantum dots and calculates:

1. **Quantum dot diameter (d)** — using empirical sizing formulas from the
   Moreels group (PbS, PbSe) and Yu/Peng group (CdS, CdSe)
2. **Size distribution width (sigma)** — via Gaussian fitting of the first
   exciton peak and error propagation on the sizing curve

**Supported quantum dot materials:**

| Material | Sizing Formula | Reference | Valid Range |
|----------|---------------|-----------|-------------|
| PbS | E₀ = 0.41 + 1/(0.0252·d² + 0.283·d) | Moreels et al., *ACS Nano* 2009 | 2.3–10 nm |
| PbSe | E₀ = 0.278 + 1/(0.016·d² + 0.209·d + 0.45) | Moreels et al., *Chem. Mater.* 2007 | 2–20 nm |
| CdS | D = poly(λ) | Yu et al., *Chem. Mater.* 2003 | 1–6 nm |
| CdSe | D = poly(λ) | Yu et al., *Chem. Mater.* 2003 | 1.5–8 nm |

**Outputs produced:**
1. **UV-Vis spectrum plot** (PNG) — spectrum with Gaussian fit overlay, HWHM
   markers, peak annotation, size and sigma in legend
2. **Word document** (`.docx`) — complete analysis report with formulas,
   step-by-step calculations for both size and sigma, references

---

## Workflow

### Step 1: Gather information from the user

Determine the quantum dot type. The script supports three methods, tried in order:

**Method 1 — User prompt**: Look for clues in the user's message:
- e.g., "PbS量子点", "CdSe nanocrystals", "分析这个PbS数据"

**Method 2 — Filename detection** (`--type auto`): The script checks the filename
for keywords (pbs, pbse, cds, cdse). Examples:
- `500pbs 20mg ml.txt` → PbS
- `CdSe_sample1.txt` → CdSe

**Method 3 — Spectral inference** (`--type auto`, fallback): If the filename
gives no clue, the script analyzes the spectrum itself:
1. Finds the most prominent long-wavelength peak
2. Compares its position against each material's calibrated range
3. Ranks all four types by confidence (high / medium / low)
4. Auto-selects the best match if confidence is "high"
5. Asks the user to specify manually if all types have "low" confidence

| Peak Wavelength | PbS | PbSe | CdS | CdSe | Verdict |
|----------------|-----|------|-----|------|---------|
| 800–2500 nm | high | high | low | low | → PbS or PbSe (PbS preferred) |
| 2500–4000 nm | low | high | low | low | → PbSe |
| 360–500 nm | low | low | high | low | → CdS |
| 510–640 nm | low | low | low | high | → CdSe |
| <360 or >4000 nm | low | low | low | low | → Ask user |

Also identify the target files. Use `Glob` with `**/*.txt` to discover all
relevant files in the workspace.

### Step 2: Ensure dependencies are installed

```bash
pip install numpy scipy matplotlib python-docx
```

### Step 3: Run the analysis script

```bash
python "<skill_dir>/scripts/uv_analysis.py" "<data_file.txt>" --type <qd_type> --output "<output_dir>"
```

**Parameters:**
- `files` — one or more `.txt` file paths (supports glob patterns like `*.txt`)
- `--type` / `-t` — `pbs`, `pbse`, `cds`, `cdse`, or `auto`
- `--output` / `-o` — output directory (default: same directory as input file)
- `--smooth` / `-s` — Savitzky-Golay window size (default: 15, must be odd)

**Axis range control:**
- `--xmin` / `--xmax` — X-axis wavelength range (nm). Auto-detected if omitted.
- `--ymin` / `--ymax` — Y-axis absorbance range. Auto-detected if omitted.
- Either bound may be omitted independently (e.g. `--ymin 0` sets only the lower bound).

**Output format & quality:**
- `--format` / `-f` — Output format(s): `png`, `pdf`, `svg`, or comma-separated list (e.g. `png,pdf,svg`). Use `all` for all three. Default: `png`.
- `--dpi` — PNG resolution: `150`, `300`, or `600`. Default: `300`.
- `--figure-width` — Figure width preset: `single` (8.5 cm) or `double` (17 cm). Default: `double`.

### Step 4: Present results to the user

After the script completes, summarize:

```
## UVTOSIZE Analysis Results

**Sample:** 1550 (PbS)
**First Exciton Peak:** 1576.0 nm (0.787 eV)
**Diameter (d):** 6.08 nm
**HWHM:** 53.3 nm (26.8 meV) | Gaussian R² = 0.979
**Size Distribution (σ):** 0.32 nm (5.3%)

**Output files:**
- Plot: ./results/UV_spectrum_1550.png
- Report: ./results/UV_Analysis_Report_1550.docx

**Sizing Formula:** E₀ = 0.41 + 1/(0.0252·d² + 0.283·d)
  — Moreels et al., ACS Nano 2009
```

---

## How the Analysis Works

### 1. Peak Detection & Sizing

Standard Savitzky-Golay smoothing → `scipy.signal.find_peaks` → first exciton
peak = most prominent peak at longest wavelength (lowest energy, 1S(e)–1S(h)).

### 2. Gaussian Fitting & HWHM

A Gaussian function is fitted to the first exciton peak region:

```
A(λ) = A₀ · exp(-(λ - λ₀)² / (2σ_g²)) + baseline
```

From the fit, we extract:
- **HWHM** = σ_g · √(2·ln(2)) ≈ 1.1774 · σ_g (wavelength domain)
- Converted to energy domain via E = 1240/λ

### 3. Size Distribution (σ) from Inhomogeneous Broadening

The fundamental relationship (error propagation on sizing curve):

```
σ_d = HWHM(E) / |dE/dd|_{d=d_mean}
```

where `|dE/dd|` is the absolute derivative of the sizing curve evaluated at the
mean diameter.

**References:**
- Nikolaev & Averkiev, *Appl. Phys. Lett.* 95, 263107 (2009)
- Wu et al., *Appl. Phys. Lett.* 51, 710 (1987)

### 4. Sizing Curve Derivatives

| Material | dE/dd Formula |
|----------|--------------|
| PbS | dE/dd = -(0.0504·d + 0.283) / (0.0252·d² + 0.283·d)² |
| PbSe | dE/dd = -(0.032·d + 0.209) / (0.016·d² + 0.209·d + 0.45)² |
| CdS | dE/dD = (-1240/λ²) / (dD/dλ), where dD/dλ = -1.9956×10⁻⁷·λ² + 3.9114×10⁻⁴·λ - 9.2352×10⁻² |
| CdSe | dE/dD = (-1240/λ²) / (dD/dλ), where dD/dλ = 6.4488×10⁻⁹·λ³ - 7.9725×10⁻⁶·λ² + 3.2484×10⁻³·λ - 0.4277 |

### 5. Full Sizing Formulas

**PbS** (Moreels 2009):
```
E₀ = 0.41 + 1/(0.0252·d² + 0.283·d)
→ 0.0252·d² + 0.283·d - 1/(E₀-0.41) = 0 (quadratic, solve for d)
```

**PbSe** (Moreels 2007):
```
E₀ = 0.278 + 1/(0.016·d² + 0.209·d + 0.45)
→ 0.016·d² + 0.209·d + (0.45 - 1/(E₀-0.278)) = 0 (quadratic)
```

**CdS** (Yu 2003, λ in nm → D in nm):
```
D = -6.6521×10⁻⁸·λ³ + 1.9557×10⁻⁴·λ² - 9.2352×10⁻²·λ + 13.29
```

**CdSe** (Yu 2003, λ in nm → D in nm):
```
D = 1.6122×10⁻⁹·λ⁴ - 2.6575×10⁻⁶·λ³ + 1.6242×10⁻³·λ² - 0.4277·λ + 41.57
```

---

## Data File Format

The script handles UV-Vis `.txt` files with this typical format:

```
"Sample Name - RawData"
"Wavelength(nm)","Absorbance"
800.00,0.95910
802.00,0.95387
...
```

Auto-detects encoding (UTF-8, UTF-16 LE/BE, GBK), strips headers/quotes,
handles comma/tab/whitespace delimiters.

---

## Troubleshooting

**Gaussian fit fails (R² < 0.8)**: Peak may be asymmetric or overlapping with
other transitions. Try narrowing the fitting window; the script automatically
uses a 150 nm window around the peak center.

**HWHM negative**: Indicates the energy conversion is incorrect — this is a bug
(should be fixed in current version).

**Sigma unreasonably large (>20%)**: The sample may have genuine polydispersity,
or the first exciton peak may be overlapped with other transitions, artificially
broadening the Gaussian fit. Verify with TEM if possible.

**Peak outside calibrated range**: The script will print warnings. For PbS/PbSe,
the first exciton peak MUST be in the NIR region. If your measurement only goes
to 800 nm and you're analyzing PbS, the true first exciton peak is likely beyond
your measurement range.
