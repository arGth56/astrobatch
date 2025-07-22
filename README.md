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
python3.11 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install astrobatch     # from source for now
```

**pySiril** is not on PyPI – install the wheel provided by the Siril
team:

```bash
pip install https://gitlab.com/-/project/20510105/uploads/8224707c29669f255ad43da3b93bc5ec/pysiril-0.0.15-py3-none-any.whl
```

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

# Upload to STDWeb using token in $STDWEB_TOKEN
astrobatch --upload --root /data/2025-07-16/LIGHT

# Analyse results (transient search, sky maps)
astrobatch --analyse --root /data/2025-07-16/LIGHT
```

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