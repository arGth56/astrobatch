# astrobatch

Advanced batch pipeline for professional astronomical image processing.

---

## 1. Overview
`astrobatch` automates the **entire nightly workflow** for small- to
medium-sized observatories:

1. Split raw FITS images into a calibrated folder tree
2. Calibrate & stack each exposure set via *Siril* (**CLI**)
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

Choose your platform:

| OS | Quick path |
|----|------------|
| **macOS** | `brew install --cask siril` then follow Python steps below (section “Siril & pySiril” has details). |
| **Windows** | Install Siril MSI, then continue with Python/pySiril (see platform notes). |
| **Debian / Ubuntu** | `chmod u+x scripts/install_debian.sh && ./scripts/install_debian.sh` – sets up everything automatically. |

If you prefer manual setup continue with the generic instructions below, then consult section 2.3.1 for per-OS specifics.

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

### 2.3. Siril (CLI)

Install Siril ≥1.2 (GUI or CLI); the pipeline uses `siril-cli` for head-less operation. If the binary is on `$PATH` you’re good to go. No Python wrapper is needed.

```bash
# Debian/Ubuntu:
sudo apt install siril siril-cli   # or download the AppImage

# macOS (Homebrew):
brew install --cask siril
```

If `siril-cli` isn’t found the calibration step is skipped with a clear error message.

#### 2.3.1 Platform notes

| Platform | Recommended install | Notes |
|----------|---------------------|-------|
| **macOS** | `brew install --cask siril`<br/>then install Siril | Home-brew ships the up-to-date build. |
| **Windows** | Download the MSI from the [official site](https://siril.org/download/) and run the installer.<br/>Then:<br/>`pip install pysiril-0.0.15-py3-none-any.whl` | pySiril looks for `C:\Program Files\Siril\bin\siril.exe` automatically. No extra display configuration needed. |
| **Linux (Debian/Ubuntu)** | **Option A (recommended)** – AppImage:<br/>`wget https://free-astro.org/download/siril-1.2.1-linux64.appimage -O /usr/local/bin/siril`<br/>`chmod +x /usr/local/bin/siril`<br/>Install pySiril wheel.<br/><br/>**Option B** – PPA / distro package:<br/>`sudo add-apt-repository ppa:lock042/siril && sudo apt update && sudo apt install siril siril-cli` | On head-less servers you either need:<br/>• `siril-cli` (no X needed), **or**<br/>• Xvfb: `sudo apt install xvfb` and run `xvfb-run python -m astrobatch.cli --calibrate …` |

> Tip: To make Xvfb permanent add `export DISPLAY=:99` in your shell profile and
> start `Xvfb :99 -screen 0 1024x768x24 -nolisten tcp -ac &` on boot (systemd
> service or rc.local).

### 2.4. Astrometry index catalogues

---

### 2.0. Debian/Ubuntu quick installer

If you’re on Debian or Ubuntu, the repository ships a convenience script that
sets up everything (system packages, virtual-env, Siril, pySiril):

```bash
chmod u+x scripts/install_debian.sh
./scripts/install_debian.sh
```

After it finishes:
```bash
source .venv/bin/activate
python -m astrobatch.cli --help
```

You can still follow the manual steps below if you prefer full control.

---

## 3. Quick-start
```bash
# Split raw images (night automatically detected from path)
astrobatch --split --root /data/2025-07-16/LIGHT

# Calibrate via Siril CLI (head-less servers)
astrobatch --calibrate --root /data/2025-07-16/LIGHT

# Upload to STDWeb (requires API token)
#
#   export STDWEB_API_TOKEN="<your-token>"
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