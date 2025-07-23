#!/usr/bin/env bash
# install_debian.sh – one-stop installer for astrobatch on Debian / Ubuntu
#
# This script:
#   1. Installs required system packages (Python, netpbm, etc.)
#   2. Creates a local virtual-environment (.venv) in the repo root
#   3. Installs Python dependencies from requirements.txt
#   4. Installs pySiril
#   5. Downloads the Siril AppImage if no Siril binary is present
#
# Run with:  bash scripts/install_debian.sh
# ------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"

appimage_url="https://free-astro.org/download/siril-1.2.1-linux64.appimage"
appimage_path="tools/siril.appimage"

# ------------------------------------------------------------------
# 1. System packages
# ------------------------------------------------------------------

echo "\n🔧 Installing system packages (sudo)…"
sudo apt update
sudo apt install -y python3 python3-venv python3-pip git wget curl netpbm

# Optional: astrometry.net (comment out if you already have it)
sudo apt install -y astrometry.net || true

# ------------------------------------------------------------------
# 2. Create / activate virtual-environment
# ------------------------------------------------------------------
if [[ ! -d .venv ]]; then
    python3 -m venv .venv
fi
source .venv/bin/activate
python -m pip install --upgrade pip

# ------------------------------------------------------------------
# 3. Python deps
# ------------------------------------------------------------------

pip install -r requirements.txt

# ------------------------------------------------------------------
# 4. pySiril (server-mode wrapper)
# ------------------------------------------------------------------

echo "\n📦 Installing pySiril wrapper …"
pip install https://gitlab.com/-/project/20510105/uploads/8224707c29669f255ad43da3b93bc5ec/pysiril-0.0.15-py3-none-any.whl

# ------------------------------------------------------------------
# 5. Siril binary (AppImage) – only if not already in PATH
# ------------------------------------------------------------------

if ! command -v siril &>/dev/null && ! command -v siril-cli &>/dev/null; then
    echo "\n⬇️  Downloading Siril AppImage …"
    mkdir -p tools
    curl -L "$appimage_url" -o "$appimage_path"
    chmod +x "$appimage_path"
    echo "   AppImage saved to $appimage_path"
    echo "   Creating symlink in /usr/local/bin (needs sudo) …"
    sudo ln -sf "$REPO_ROOT/$appimage_path" /usr/local/bin/siril
fi

# ------------------------------------------------------------------
# Done
# ------------------------------------------------------------------

echo -e "\n✅ All done!  Activate your environment with:\n   source .venv/bin/activate\nThen run:\n   python -m astrobatch.cli --help\n" 