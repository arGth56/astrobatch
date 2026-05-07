#!/usr/bin/env python3
"""
Retroactively fix FILTER headers in NINA-saved FITS files.

NINA sometimes writes the wrong filter (from its internal sequence context) to all saved
frames, even when the physical filterwheel is changed between captures via the API.
This script reconstructs the correct filter assignment from capture timestamps.

Usage:
    python tools/fix_filter_headers.py --date 2026-03-21 --dry-run
    python tools/fix_filter_headers.py --date 2026-03-21
"""

import argparse
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from astropy.io import fits

NAS_INPUT = os.environ.get("NAS_INPUT", "/mnt/nas/input/pyl/astro/input")

# Sequence start cutoff (UTC). Files captured AFTER this time are from the 23:34 run
# and may need filter correction. Files before are treated as genuinely RP.
#
# Each target's post-cutoff frames are sorted by DATE-OBS; first N=G, next N=BP, last N=RP.
# Override via --filters and --count if your sequence used different settings.

def parse_args():
    p = argparse.ArgumentParser(description="Fix FITS FILTER headers for a night's SNAPSHOT data")
    p.add_argument("--date", required=True, help="NAS date folder, e.g. 2026-03-21")
    p.add_argument("--cutoff-utc", default=None,
                   help="UTC datetime after which to fix headers (ISO format, e.g. 2026-03-21T22:34:00). "
                        "Files before this are assumed to be correctly tagged.")
    p.add_argument("--filters", default="G,BP,RP",
                   help="Ordered comma-separated filter names used in the sequence (default: G,BP,RP)")
    p.add_argument("--count", type=int, default=10,
                   help="Number of frames per filter per target (default: 10)")
    p.add_argument("--dry-run", action="store_true",
                   help="Print what would change without writing")
    return p.parse_args()


def dateobs_to_utc(dateobs_str):
    """Parse FITS DATE-OBS string to UTC datetime."""
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S"):
        try:
            return datetime.strptime(dateobs_str[:26], fmt).replace(tzinfo=timezone.utc)
        except ValueError:
            continue
    return None


def main():
    args = parse_args()
    snapshot_dir = Path(NAS_INPUT) / args.date / "SNAPSHOT"
    if not snapshot_dir.exists():
        print(f"ERROR: {snapshot_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    filter_seq = [f.strip() for f in args.filters.split(",")]
    n_per_filter = args.count

    cutoff_dt = None
    if args.cutoff_utc:
        cutoff_dt = datetime.fromisoformat(args.cutoff_utc).replace(tzinfo=timezone.utc)

    # Load all FITS files and read headers
    print(f"Scanning {snapshot_dir} ...")
    records = []
    fits_files = sorted(snapshot_dir.glob("*.fit")) + sorted(snapshot_dir.glob("*.fits"))
    for fpath in fits_files:
        try:
            hdr = fits.getheader(str(fpath), memmap=False)
        except Exception as e:
            print(f"  SKIP {fpath.name}: cannot read header ({e})")
            continue
        dateobs = hdr.get("DATE-OBS", "")
        dt = dateobs_to_utc(dateobs)
        records.append({
            "path": fpath,
            "object": hdr.get("OBJECT", "?"),
            "filter_current": hdr.get("FILTER", "?"),
            "dateobs": dateobs,
            "dt": dt,
        })

    print(f"  Found {len(records)} FITS files\n")

    # Group by OBJECT
    from collections import defaultdict
    by_object = defaultdict(list)
    for r in records:
        by_object[r["object"]].append(r)

    changes = []
    for obj, recs in sorted(by_object.items()):
        recs.sort(key=lambda r: r["dt"] or datetime.min.replace(tzinfo=timezone.utc))

        # Split into: early run (before cutoff) and main sequence (after cutoff)
        if cutoff_dt:
            early = [r for r in recs if r["dt"] and r["dt"] < cutoff_dt]
            main_run = [r for r in recs if not r["dt"] or r["dt"] >= cutoff_dt]
        else:
            early = []
            main_run = recs

        print(f"Target: {obj!r}  total={len(recs)}  early={len(early)}  main_run={len(main_run)}")

        # Assign filters to main-run frames: first N→filter[0], next N→filter[1], etc.
        for idx, rec in enumerate(main_run):
            filter_slot = idx // n_per_filter
            if filter_slot < len(filter_seq):
                correct_filter = filter_seq[filter_slot]
            else:
                # Extra frames beyond expected count — tag as last filter
                correct_filter = filter_seq[-1]

            if rec["filter_current"] != correct_filter:
                changes.append((rec["path"], rec["filter_current"], correct_filter, rec["dateobs"]))
                print(f"  [{idx+1:3d}] {rec['path'].name[:40]:40}  "
                      f"{rec['filter_current']!r:6} → {correct_filter!r}  ({rec['dateobs'][:19]})")
            else:
                print(f"  [{idx+1:3d}] {rec['path'].name[:40]:40}  {rec['filter_current']!r:6}  OK")

        print()

    print(f"\n{'DRY RUN — ' if args.dry_run else ''}Changes to apply: {len(changes)}")

    if not changes:
        print("Nothing to fix.")
        return

    if args.dry_run:
        print("Re-run without --dry-run to apply.")
        return

    print("Applying...")
    ok = 0
    for fpath, old_f, new_f, _ in changes:
        try:
            # 1. Patch FITS header
            fits.setval(str(fpath), "FILTER", value=new_f)

            # 2. Rename file: replace filter token in filename
            # Format: YYYY-MM-DD_HH-MM-SS_FILTER_TEMP_DURATION_SEQ.ext
            stem = fpath.stem          # filename without extension
            ext  = fpath.suffix        # .fits / .fit
            parts = stem.split("_")
            if len(parts) >= 3 and parts[2] == old_f:
                parts[2] = new_f
                new_name = "_".join(parts) + ext
                new_path = fpath.parent / new_name
                fpath.rename(new_path)
                print(f"  FIXED {fpath.name} → {new_name}  (FILTER {old_f!r} → {new_f!r})")
            else:
                print(f"  FIXED header only: {fpath.name}  (FILTER {old_f!r} → {new_f!r})")
            ok += 1
        except Exception as e:
            print(f"  ERROR {fpath.name}: {e}", file=sys.stderr)

    print(f"\nDone. Fixed {ok}/{len(changes)} files.")


if __name__ == "__main__":
    main()
