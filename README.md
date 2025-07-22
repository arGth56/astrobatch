# astrobatch

Advanced batch pipeline for professional astronomical image processing.

---

## 1. Overview
`astrobatch` automates the **entire nightly workflow** for small- to
medium-sized observatories:

1. Split raw FITS images into a calibrated folder tree
2. Calibrate & stack each exposure set via *Siril* (through **pySiril**)
3. Plate-solve every stacked image (Astrometry.net)
4. Upload to **STDWeb** for photometry & template subtraction
5. Retrieve and analyse the results, including transient matching

The pipeline is controlled by a **single CLI**:

```bash
astrobatch --split --calibrate --upload --analyse --night 2025-07-16
```

Any combination of steps can be run; inputs and outputs are reproducible
and fully logged.

---

## 2. Installation

### 2.1. System prerequisites

| Component | macOS (Homebrew) | Ubuntu / Debian |
|-----------|------------------|-----------------|
| Python ≥ 3.10 | `brew install python@3.11` | `sudo apt install python3.11 python3.11-venv` |
| Siril 1.2+ | `brew install --cask siril` | Official AppImage [link](https://siril.org/download/) |
| Astrometry.net | `brew install astrometry-net netpbm` | `sudo apt install astrometry.net netpbm` |
| Index files | download 4200-4210 series (see below) | same |

### 2.2. Python environment

```bash
# 1) create and activate a virtual-environment (any Python ≥3.10)
python3 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip

# 2) install runtime dependencies directly from the source tree
pip install -r requirements.txt

# 3) (optional, for development) install the package in *editable* mode so the
#    `astrobatch` command is available and live-reloads your changes
pip install -e .
```

> A future `pip install astrobatch` will be available once the project is
> published on PyPI. For now you must work from the cloned repository as shown
> above.

### 2.3. Astrometry index catalogues
Download the 4200–4210 series (2–8° FOV) into the directory listed in
`astrometry.cfg` (`/opt/homebrew/Cellar/astrometry-net/*/data` on macOS).

Example one-liner (4200 only):
```bash
cd /opt/homebrew/Cellar/astrometry-net/0.97/data
for n in {00..47}; do curl -O http://data.astrometry.net/4200/index-4200-$n.fits; done
```

---

## 3. Quick-start
```bash
# Split raw images (night automatically detected from path)
astrobatch --split --root /data/2025-07-16/LIGHT

# Calibrate via Siril (server mode)
astrobatch --calibrate --root /data/2025-07-16/LIGHT

# Upload to STDWeb using token in $STDWEB_API_TOKEN
astrobatch --upload --root /data/2025-07-16/LIGHT

# Analyse results (transient search, sky maps)
astrobatch --analyse --root /data/2025-07-16/LIGHT
```

---

## 3.1. Authentication / secrets

STDWeb requires an **API token** for uploads and automation.  Create one in your
STDWeb user profile (‣ *Settings → API Tokens → New*).  The pipeline expects it
in the environment variable `STDWEB_API_TOKEN`:

```bash
export STDWEB_API_TOKEN="<paste-your-token-here>"
```

For convenience you can store it in a `.env` file or in your shell’s start-up
script (`~/.zprofile`, `~/.bashrc`, …).  **.env files are already excluded from
Git via `.gitignore`, so the secret will never be committed.**

If the variable is missing `astrobatch` aborts with a clear error message.

---

## 4. Command-line reference

```
astrobatch [OPTIONS]

Options:
  --split           Organise raw FITS into folder tree
  --calibrate       Run Siril calibration/stacking
  --upload          Upload res.fit to STDWeb (photometry S/N=10, ZTF DR7)
  --analyse         Download & analyse STDWeb results (reports → Output/)
  --plate-solve     Plate-solve all res.fit that lack WCS headers
  --mosaic          Build a sky mosaic of the night (TAN projection)
  --start-server    Launch Flask review server for candidate labelling
  --port INT        Custom port for --start-server (default 5100)
  --night YYYY-MM-DD  Override night root directory
  --root PATH       Explicit data root (default derived from --night)
  --dry-run         Show what would be executed
  -v, --verbose     Increase logging
  -h, --help        Show this message and exit
```

---

## 5. Architecture

```
astrobatch/
├── spliter.py          # (legacy) processing engine – will be modularised
├── analyse.py          # thin wrapper around legacy analysis
├── db.py               # SQLite CandidateDB helper
├── server.py           # Flask API for manual labelling
├── models/             # ML models used by analyse
├── catalog/NGC.csv     # NGC galaxy catalogue
└── cli.py              # single entry-point (dispatches the above)

The refactor keeps legacy code operational inside the package while new
modules are extracted incrementally.

Each module is importable for custom pipelines; all functions carry full
type hints and numpy-style docstrings.

---

## 6. Development & contribution
1. Fork the repository and create a feature branch.
2. Follow the *pre-commit* hooks (`black`, `isort`, `flake8`).
3. Write unit tests under `tests/` (pytest).
4. Submit a pull request; CI must pass (GitHub Actions on macOS & Ubuntu).

---

## 7. License
MIT License © 2025 – contribute back improvements! 