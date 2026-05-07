#!/usr/bin/env python3
"""
verify_filter_labels.py

Cross-reference sequence logs (filter-change events) with raw FITS filenames
to determine what filter was actually active when each frame was taken.

Only covers nights where a seq-*.log exists. Nights without logs are skipped.

Usage:
  python3 verify_filter_labels.py --dry-run    # report mismatches only
  python3 verify_filter_labels.py --execute    # rename files + patch headers
"""

import argparse, re, sys, shutil, datetime
from pathlib import Path
from astropy.io import fits

LOG_DIR  = Path("/home/pyl/Documents/astrobatch/server/logs")
NAS_DIR  = Path("/mnt/nas/input/pyl/astro/input")
TMP      = "VFTMP"

# ID → canonical name (from NINA config confirmed in logs)
SLOT_NAME = {0: "RP", 1: "BP", 2: "G", 3: "dark"}
NAME_SLOT = {v: k for k, v in SLOT_NAME.items()}

# ── Parse logs ──────────────────────────────────────────────────────────────

RE_CHANGE  = re.compile(
    r'(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+Z).*'
    r'(?:Changing to filter|Switching to G filter)\s*[:\(]?\s*'
    r'(?P<name>[A-Za-z]+)\s*\(ID\s*(?P<id>\d+)\)')

def parse_logs():
    """Return sorted list of (datetime_utc, filter_name, slot_id)."""
    events = []
    for log_file in sorted(LOG_DIR.glob("seq-*.log")):
        for line in log_file.read_text(errors="replace").splitlines():
            m = RE_CHANGE.search(line)
            if m:
                ts = datetime.datetime.fromisoformat(
                    m.group(1).replace("Z", "+00:00"))
                name = m.group("name").strip()
                slot = int(m.group("id"))
                events.append((ts, name, slot))
    events.sort()
    return events

def active_filter_at(events, ts_naive):
    """Given a naive UTC datetime, return (name, slot) of the active filter."""
    ts = ts_naive.replace(tzinfo=datetime.timezone.utc)
    active = None
    for ev_ts, name, slot in events:
        if ev_ts <= ts:
            active = (name, slot)
        else:
            break
    return active  # None if before all log events

# ── FITS filename helpers ────────────────────────────────────────────────────

RE_FNAME = re.compile(
    r'^(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})_([A-Za-z]+)_(.+)$')

def parse_fits_name(name):
    """Return (date_str, time_str, filter_label, rest) or None."""
    m = RE_FNAME.match(Path(name).stem)
    if not m: return None
    return m.group(1), m.group(2), m.group(3), m.group(4)

def fits_utc(date_str, time_str):
    """Parse YYYY-MM-DD and HH-MM-SS → naive UTC datetime."""
    return datetime.datetime.strptime(
        f"{date_str}T{time_str.replace('-', ':')}", "%Y-%m-%dT%H:%M:%S")

def new_name(stem, old_filter, new_filter):
    return stem.replace(f"_{old_filter}_", f"_{new_filter}_", 1)

# ── Main ────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run",  action="store_true")
    p.add_argument("--execute",  action="store_true")
    args = p.parse_args()
    if not args.dry_run and not args.execute:
        p.error("Specify --dry-run or --execute")

    dry = args.dry_run
    print(f"Mode: {'DRY-RUN' if dry else 'EXECUTE'}\n")

    events = parse_logs()
    if not events:
        print("ERROR: no log events found"); sys.exit(1)
    first_ev = events[0][0]
    last_ev  = events[-1][0]
    print(f"Log coverage: {first_ev.date()} → {last_ev.date()}")
    print(f"Total filter-change events: {len(events)}\n")

    # Map slot_id → canonical name from NINA logs (cross-check)
    seen_slots = {}
    for _, name, slot in events:
        seen_slots.setdefault(slot, set()).add(name)
    print("Slot mapping seen in logs:")
    for slot in sorted(seen_slots):
        print(f"  ID {slot} → {seen_slots[slot]}")
    print()

    files = sorted(NAS_DIR.glob("*/SNAPSHOT/*.fits"))
    in_range = out_range = mismatch = correct = errors = fixed = 0

    changes = []  # (path, old_filter, new_filter)

    for fpath in files:
        parsed = parse_fits_name(fpath.name)
        if not parsed: continue
        date_s, time_s, label, rest = parsed

        # Skip non-BP/RP/G/dark files
        if label.upper() not in ("BP", "RP", "G", "DARK", "G", "dark"):
            continue

        ts = fits_utc(date_s, time_s)
        # Add 60s buffer: frame starts ~60s after filter confirmed
        ts_check = ts + datetime.timedelta(seconds=60)

        if ts_check < first_ev.replace(tzinfo=None):
            out_range += 1
            continue

        in_range += 1
        active = active_filter_at(events, ts_check)
        if active is None:
            print(f"  SKIP (before log start): {fpath.name}")
            continue

        log_filter, log_slot = active

        if label.upper() == log_filter.upper():
            correct += 1
        else:
            mismatch += 1
            print(f"  MISMATCH  {fpath.name}")
            print(f"    filename says '{label}'  but  log says '{log_filter}' (slot {log_slot})")
            changes.append((fpath, label, log_filter))

    print(f"\nSummary:")
    print(f"  Files in log date range : {in_range}")
    print(f"  Files outside log range  : {out_range} (skipped — no logs)")
    print(f"  Correct label            : {correct}")
    print(f"  Mismatched label         : {mismatch}")

    if not changes:
        print("\nNo mismatches — all labels match logs.")
        return

    if dry:
        print(f"\n[DRY-RUN] Would rename/patch {len(changes)} files.")
        return

    # Execute: 3-pass safe rename + header patch
    print(f"\nApplying {len(changes)} corrections...")
    for fpath, old_f, new_f in changes:
        tmp_name = fpath.parent / fpath.name.replace(f"_{old_f}_", f"_{TMP}_", 1)
        new_stem = new_name(fpath.stem, old_f, new_f)
        new_path = fpath.parent / (new_stem + fpath.suffix)
        try:
            fpath.rename(tmp_name)
            with fits.open(tmp_name, mode="update") as h:
                h[0].header["FILTER"] = new_f
                h.flush()
            tmp_name.rename(new_path)
            fixed += 1
            print(f"  FIXED: {fpath.name} → {new_path.name}")
        except Exception as e:
            print(f"  ERROR: {fpath.name}: {e}", file=sys.stderr)
            errors += 1

    print(f"\nDone. Fixed {fixed}, errors {errors}.")

if __name__ == "__main__":
    main()
