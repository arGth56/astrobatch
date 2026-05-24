#!/bin/bash
# astro-sync.sh — Mirror new FITS files from astropc USB disk to NAS
# Syncs only the newest date folders (fast path); polls every 5s.

SRC="/mnt/astropc/"
DST="/mnt/nas/input/pyl/astro/input/"
LOG="/home/pyl/Documents/astrobatch/server/logs/astro-sync.log"
LOCKFILE="/tmp/astro-sync.lock"
POLL_SEC=5
RECENT_FOLDERS=2
FULL_SYNC_EVERY=720   # ~1 h at 5 s — catch older nights on D:

log() {
  echo "$(date '+%Y-%m-%dT%H:%M:%S')  $*" | tee -a "$LOG"
}

# YYYY-MM-DD folders on D:, newest first
recent_date_dirs() {
  ls -1 "$SRC" 2>/dev/null \
    | grep -E '^[0-9]{4}-[0-9]{2}-[0-9]{2}$' \
    | sort -r \
    | head -n "$RECENT_FOLDERS"
}

rsync_opts=(
  --archive
  --update
  --human-readable
  --exclude 'System Volume Information'
  --exclude '.Trash*'
  --exclude '.Spotlight*'
  --exclude '.fseventsd'
  --exclude '.*'
  --exclude '*.tmp'
  --exclude '*.part'
  --exclude '*.FW*'
  --out-format='  [%t] %n (%b)'
)

run_rsync() {
  local label="$1"
  local from="$2"
  local to="$3"
  rsync "${rsync_opts[@]}" "$from" "$to" 2>&1
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
log "  poll: ${POLL_SEC}s  recent folders: $RECENT_FOLDERS  full sync every: ${FULL_SYNC_EVERY} cycles"

cycle=0

while true; do
  cycle=$((cycle + 1))

  if ! mountpoint -q /mnt/astropc; then
    log "[WARN] /mnt/astropc not mounted — skipping"
    sleep "$POLL_SEC"
    continue
  fi
  if ! mountpoint -q /mnt/nas/input; then
    log "[WARN] /mnt/nas/input not mounted — skipping"
    sleep "$POLL_SEC"
    continue
  fi

  RESULT=""
  if [ $((cycle % FULL_SYNC_EVERY)) -eq 0 ]; then
    log "[FULL] periodic sync of entire D: tree"
    RESULT=$(run_rsync "full" "$SRC" "$DST")
  else
    mapfile -t DIRS < <(recent_date_dirs)
    if [ "${#DIRS[@]}" -eq 0 ]; then
      sleep "$POLL_SEC"
      continue
    fi
    for folder in "${DIRS[@]}"; do
      src_dir="${SRC%/}/$folder/"
      dst_dir="${DST%/}/$folder/"
      if [ ! -d "$src_dir" ]; then
        continue
      fi
      PART=$(run_rsync "$folder" "$src_dir" "$dst_dir")
      if [ -n "$PART" ]; then
        RESULT+="$PART"$'\n'
      fi
    done
  fi

  if echo "$RESULT" | grep -q '^\s*\['; then
    log "[SYNC]"
    echo "$RESULT" >> "$LOG"
  fi

  sleep "$POLL_SEC"
done
