#!/usr/bin/env python3
"""
Live test for the FITS FILTER header patch logic.

Workflow:
  1. Switch filterwheel to G filter
  2. Capture a 1-second test frame (NINA will tag it with wrong filter in header)
  3. Find the new file on the NAS
  4. Record the (wrong) FILTER header BEFORE patching
  5. Apply the patch (same logic as patchFitsFilterHeaders in server.js)
  6. Read the header back and verify it says G

Usage:
    cd /home/pyl/Documents/astrobatch
    .venv/bin/python tools/test_filter_patch.py
"""

import os
import sys
import time
from pathlib import Path
from datetime import datetime, timezone

import requests
from astropy.io import fits

NINA_HOST = os.environ.get("NINA_HOST", "192.168.1.174")
NINA_PORT = os.environ.get("NINA_PORT", "1888")
NINA_BASE = f"http://{NINA_HOST}:{NINA_PORT}"
NAS_INPUT  = os.environ.get("NAS_INPUT", "/mnt/nas/input/pyl/astro/input")

TEST_FILTER = "G"       # filter to test with
TEST_DURATION = 1       # seconds
TEST_GAIN = 10


def nina_get(path, **params):
    url = f"{NINA_BASE}{path}"
    r = requests.get(url, params=params, timeout=30)
    try:
        return r.json()
    except Exception:
        return {"raw": r.text, "status_code": r.status_code}


def all_snapshot_dirs():
    """Return all .../DATE/SNAPSHOT/ directories under NAS_INPUT."""
    base = Path(NAS_INPUT)
    dirs = []
    if not base.exists():
        return dirs
    for d in base.iterdir():
        snap = d / "SNAPSHOT"
        if d.is_dir() and snap.exists():
            dirs.append(snap)
    return dirs


def new_fits_files(since_ts):
    """Return .fit/.fits files in any SNAPSHOT dir modified at or after since_ts (epoch seconds)."""
    result = []
    for snap_dir in all_snapshot_dirs():
        for p in snap_dir.iterdir():
            if p.suffix.lower() not in (".fit", ".fits"):
                continue
            if p.stat().st_mtime >= since_ts:
                result.append(p)
    return sorted(result, key=lambda p: p.stat().st_mtime)


def patch_filter(fpath, filter_name):
    with fits.open(str(fpath), mode="update") as hdul:
        hdul[0].header["FILTER"] = filter_name
        hdul.flush()


def read_filter(fpath):
    return fits.getheader(str(fpath), memmap=False).get("FILTER", "?")


def main():
    print(f"=== Live filter-patch test ===")
    print(f"NINA:     {NINA_BASE}")
    print()

    # ── 1. Check NINA is reachable ────────────────────────────────────────────
    print("[1] Checking NINA connection...")
    info = nina_get("/v2/api/equipment/info")
    if "raw" in info:
        print(f"    NINA unreachable: {info.get('raw','')[:120]}")
        sys.exit(1)
    fw_info = info.get("Response", {}).get("FilterWheel", {})
    available = [f.get("Name") for f in fw_info.get("AvailableFilters", [])]
    print(f"    Connected ✓  Available filters: {available}")

    if TEST_FILTER not in available:
        print(f"    ERROR: filter '{TEST_FILTER}' not available (got {available})")
        sys.exit(1)

    # ── 2. Switch to test filter ──────────────────────────────────────────────
    print(f"\n[2] Switching filterwheel to {TEST_FILTER}...")
    r = nina_get("/v2/api/equipment/filterwheel/change-filter", filterName=TEST_FILTER)
    success = r.get("Success", r.get("success", False))
    if success:
        print(f"    Filter → {TEST_FILTER} ✓")
    else:
        print(f"    WARNING: {r}")
    time.sleep(2)

    # ── 3. Record capture start time ─────────────────────────────────────────
    capture_start = time.time()
    print(f"\n[3] Triggering {TEST_DURATION}s test exposure (gain {TEST_GAIN})...")
    cr = nina_get("/v2/api/equipment/camera/capture",
                  duration=TEST_DURATION, gain=TEST_GAIN,
                  save=True, targetName="FILTER_TEST", filter=TEST_FILTER)
    capture_ok = cr.get("Success", cr.get("success", False))
    print(f"    Capture response: Success={capture_ok}  Error={cr.get('Error','none')}")
    if not capture_ok:
        print("    WARNING: capture may have failed — waiting anyway")

    # ── 4. Wait for file to appear ────────────────────────────────────────────
    print(f"\n[4] Waiting up to 30s for file to appear in any SNAPSHOT dir ...")
    deadline = time.time() + 30
    new_files = []
    while time.time() < deadline:
        new_files = new_fits_files(capture_start)
        if new_files:
            break
        time.sleep(2)
        print("    ...", end="", flush=True)
    print()

    if not new_files:
        print("    ERROR: No new FITS file appeared — is camera saving to NAS?")
        sys.exit(1)

    f = new_files[-1]  # most recent
    print(f"    Found: {f.name}")

    # Wait until file size is stable (NINA has finished writing it)
    print("    Waiting for file to be fully written...", end="", flush=True)
    prev_size = -1
    for _ in range(15):
        size = f.stat().st_size
        if size == prev_size and size > 0:
            break
        prev_size = size
        time.sleep(2)
        print(".", end="", flush=True)
    print(f" {f.stat().st_size:,} bytes ✓")

    # ── 5. Read FILTER BEFORE patch ───────────────────────────────────────────
    before = read_filter(f)
    print(f"\n[5] FILTER header BEFORE patch: {before!r}")
    if before == TEST_FILTER:
        print(f"    Interesting! NINA already wrote the correct filter ({before}) — "
              "the bug may be fixed in your NINA version, but patching still won't hurt.")
    else:
        print(f"    Bug confirmed: NINA wrote '{before}' instead of '{TEST_FILTER}'")

    # ── 6. Apply patch ────────────────────────────────────────────────────────
    print(f"\n[6] Patching FILTER → {TEST_FILTER!r} ...")
    patch_filter(f, TEST_FILTER)
    print("    Patch applied ✓")

    # ── 7. Read FILTER AFTER patch ────────────────────────────────────────────
    after = read_filter(f)
    print(f"\n[7] FILTER header AFTER patch: {after!r}")

    # ── 8. Result ─────────────────────────────────────────────────────────────
    print()
    if after == TEST_FILTER:
        print(f"=== TEST PASSED ✓ ===")
        if before != TEST_FILTER:
            print(f"    NINA wrote '{before}' → patch corrected it to '{after}'.")
            print(f"    patchFitsFilterHeaders() will work correctly during sequences.")
        else:
            print(f"    NINA wrote the correct filter AND patch preserved it.")
    else:
        print(f"=== TEST FAILED ✗ ===")
        print(f"    Expected FILTER='{TEST_FILTER}', got '{after}' after patching.")
        sys.exit(1)


if __name__ == "__main__":
    main()
