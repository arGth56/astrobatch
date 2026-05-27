#!/usr/bin/env bash
# setup-astro-sync.sh — Install persistent CIFS mount + automount + sync service
#
# Run once with sudo from the repo root:
#   sudo bash server/setup-astro-sync.sh
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CREDS_DIR="/etc/astropc"
CREDS_FILE="$CREDS_DIR/credentials"

echo "=== AstroBatch sync setup ==="
echo "    repo: $REPO_DIR"

# ── 1. Credentials file ────────────────────────────────────────────────────
if [[ ! -f "$CREDS_FILE" ]]; then
  echo ""
  echo "Creating $CREDS_FILE …"
  mkdir -p "$CREDS_DIR"
  cat > "$CREDS_FILE" << 'EOF'
username=pyl
password=
EOF
  chmod 600 "$CREDS_FILE"
  chown root:root "$CREDS_FILE"
  echo "  Created. Edit $CREDS_FILE if the share ever needs a password."
else
  echo "  $CREDS_FILE already exists — leaving untouched."
fi

# ── 2. Mount point ─────────────────────────────────────────────────────────
if [[ ! -d /mnt/astropc ]]; then
  echo "Creating /mnt/astropc …"
  mkdir -p /mnt/astropc
fi

# ── 3. Install systemd units ───────────────────────────────────────────────
echo ""
echo "Installing systemd units …"

cp "$REPO_DIR/mnt-astropc.mount"     /etc/systemd/system/mnt-astropc.mount
cp "$REPO_DIR/mnt-astropc.automount" /etc/systemd/system/mnt-astropc.automount
cp "$REPO_DIR/server/astro-sync.service" /etc/systemd/system/astro-sync.service

chmod 644 /etc/systemd/system/mnt-astropc.mount
chmod 644 /etc/systemd/system/mnt-astropc.automount
chmod 644 /etc/systemd/system/astro-sync.service

# ── 4. Reload & enable ────────────────────────────────────────────────────
echo "Reloading systemd …"
systemctl daemon-reload

echo "Enabling units …"
systemctl enable mnt-astropc.automount
systemctl enable astro-sync.service

# ── 5. Start / restart ────────────────────────────────────────────────────
echo "Starting automount …"
systemctl start mnt-astropc.automount

# Restart sync service to pick up the new unit file
if systemctl is-active --quiet astro-sync.service; then
  echo "Restarting astro-sync.service …"
  systemctl restart astro-sync.service
else
  echo "Starting astro-sync.service …"
  systemctl start astro-sync.service
fi

# ── 6. Status ─────────────────────────────────────────────────────────────
echo ""
echo "=== Status ==="
systemctl status mnt-astropc.automount --no-pager -l | head -12
echo ""
systemctl status astro-sync.service --no-pager -l | head -12

echo ""
echo "✓ Done. The mount and sync service will start automatically on every boot."
echo ""
echo "  Useful commands:"
echo "    journalctl -fu astro-sync.service          # live sync log"
echo "    systemctl status mnt-astropc.automount      # mount health"
echo "    ls /mnt/astropc/                            # browse NINA PC"
