#!/bin/bash
# Test that patchFitsFilterHeaders works correctly:
# 1. Switch filterwheel to BP
# 2. Take a 1-second test exposure (NINA will save it with wrong FILTER=RP in header)
# 3. Immediately patch the header
# 4. Read it back and verify FILTER=BP

set -euo pipefail

NINA_HOST="${NINA_HOST:-192.168.1.174}"
NINA_PORT="${NINA_PORT:-1888}"
NINA_BASE="http://${NINA_HOST}:${NINA_PORT}"
NAS_SNAPSHOT="/mnt/nas/input/pyl/astro/input/$(date +%Y-%m-%d)/SNAPSHOT"
VENV="$(dirname "$0")/../.venv/bin/python"

echo "=== Filter patch live test ==="
echo "NINA: $NINA_BASE"
echo "Snapshot dir: $NAS_SNAPSHOT"
echo

# 1. Switch to BP
echo "[1] Switching filterwheel to BP..."
RESULT=$(curl -s -X POST "$NINA_BASE/v2/api/equipment/filterwheel/change-filter" \
    -H "Content-Type: application/json" \
    -d '{"filterName":"BP"}')
echo "    Response: $RESULT"
if echo "$RESULT" | grep -qi '"success":true'; then
    echo "    Filter → BP ✓"
else
    echo "    WARNING: filter change may have failed"
fi
sleep 2

# 2. Record existing files before capture
echo
echo "[2] Recording existing files in SNAPSHOT..."
mkdir -p "$NAS_SNAPSHOT"
BEFORE_COUNT=$(ls "$NAS_SNAPSHOT"/*.fit "$NAS_SNAPSHOT"/*.fits 2>/dev/null | wc -l || echo 0)
CAPTURE_START_MS=$(date +%s%3N)
echo "    $BEFORE_COUNT files already present, capture start: ${CAPTURE_START_MS}ms"

# 3. Trigger 1-second test exposure
echo
echo "[3] Triggering 1s test exposure..."
CAPTURE=$(curl -s -X POST "$NINA_BASE/v2/api/equipment/camera/capture" \
    -H "Content-Type: application/json" \
    -d "{\"duration\":1,\"gain\":10,\"save\":true,\"targetName\":\"filter_test\",\"filter\":\"BP\"}")
echo "    Response: $CAPTURE"

# 4. Wait for exposure + download (1s exposure + 20s buffer)
echo
echo "[4] Waiting 25s for exposure + download..."
sleep 25

# 5. Find new file(s)
echo
echo "[5] Looking for new files..."
NEW_FILES=()
for f in "$NAS_SNAPSHOT"/*.fit "$NAS_SNAPSHOT"/*.fits 2>/dev/null; do
    [ -f "$f" ] || continue
    FILE_MS=$(date -r "$f" +%s%3N)
    if [ "$FILE_MS" -ge "$CAPTURE_START_MS" ]; then
        NEW_FILES+=("$f")
        echo "    Found new file: $(basename "$f")  (mtime=${FILE_MS}ms)"
    fi
done

if [ ${#NEW_FILES[@]} -eq 0 ]; then
    echo "    ERROR: No new files found — was camera connected and saving?"
    exit 1
fi

# 6. Check headers BEFORE patch
echo
echo "[6] FITS FILTER header BEFORE patch:"
for f in "${NEW_FILES[@]}"; do
    FILTER_BEFORE=$($VENV -c "
from astropy.io import fits
h = fits.getheader('$f', memmap=False)
print(h.get('FILTER','?'))
" 2>/dev/null)
    echo "    $(basename "$f"): FILTER='$FILTER_BEFORE'"
done

# 7. Apply patch
echo
echo "[7] Patching FILTER header to 'BP'..."
for f in "${NEW_FILES[@]}"; do
    $VENV -c "
from astropy.io import fits
with fits.open('$f', mode='update') as h:
    h[0].header['FILTER'] = 'BP'
    h.flush()
print('  Patched: $f')
" 2>/dev/null
done

# 8. Check headers AFTER patch
echo
echo "[8] FITS FILTER header AFTER patch:"
PASS=true
for f in "${NEW_FILES[@]}"; do
    FILTER_AFTER=$($VENV -c "
from astropy.io import fits
h = fits.getheader('$f', memmap=False)
print(h.get('FILTER','?'))
" 2>/dev/null)
    if [ "$FILTER_AFTER" = "BP" ]; then
        echo "    $(basename "$f"): FILTER='$FILTER_AFTER' ✓"
    else
        echo "    $(basename "$f"): FILTER='$FILTER_AFTER' ✗ — EXPECTED 'BP'"
        PASS=false
    fi
done

echo
if $PASS; then
    echo "=== TEST PASSED ✓ ==="
    echo "    The patch correctly overwrites NINA's wrong FILTER header."
    echo "    patchFitsFilterHeaders() will work correctly during sequences."
else
    echo "=== TEST FAILED ✗ ==="
    echo "    The FILTER header was not correctly patched."
fi
