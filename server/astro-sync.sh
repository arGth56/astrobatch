#!/bin/bash
# astro-sync.sh — Mirror new FITS files from astropc USB disk to NAS
# Runs as a daemon: polls every 30s, rsyncs any new/updated files.

SRC="/mnt/astropc/"
DST="/mnt/nas/input/pyl/astro/input/"
LOG="/home/pyl/Documents/astrobatch/server/logs/astro-sync.log"
LOCKFILE="/tmp/astro-sync.lock"

log() {
  echo "$(date '+%Y-%m-%dT%H:%M:%S')  $*" | tee -a "$LOG"
}

# Single-instance guard
if [ -e "$LOCKFILE" ] && kill -0 "$(cat "$LOCKFILE")" 2>/dev/null; then
  log "[SKIP] Already running (pid $(cat "$LOCKFILE"))"
  exit 0
fi
echo $$ > "$LOCKFILE"
trap 'rm -f "$LOCKFILE"' EXIT

log "=== astro-sync started (pid $$) ==="
log "  src: $SRC"
log "  dst: $DST"

while true; do
  # Check both mounts are alive before doing anything
  if ! mountpoint -q /mnt/astropc; then
    log "[WARN] /mnt/astropc not mounted — skipping"
    sleep 30
    continue
  fi
  if ! mountpoint -q /mnt/nas/input; then
    log "[WARN] /mnt/nas/input not mounted — skipping"
    sleep 30
    continue
  fi

  # Rsync: copy new/updated files, preserve structure, skip Windows system folder
  RESULT=$(rsync \
    --archive \
    --update \
    --human-readable \
    --exclude 'System Volume Information' \
    --exclude '*.tmp' \
    --exclude '*.part' \
    --out-format='  [%t] %n (%b)' \
    "$SRC" "$DST" 2>&1)

  # Only log if files were actually transferred
  if echo "$RESULT" | grep -q '^\s*\['; then
    log "[SYNC]"
    echo "$RESULT" >> "$LOG"
  fi

  sleep 30
done
