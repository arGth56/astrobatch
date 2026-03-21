#!/usr/bin/env python3
"""
migrate_results_to_nas.py

One-off script: scan data/ for completed pipeline results (res.fit + res_preview.png),
move them to NAS_OUTPUT/<date>/<target>/<filter>/ and delete Siril intermediaries.

Usage (from repo root, venv active):
    python -m astrobatch.migrate_results_to_nas [--dry-run]
"""

import argparse
import os
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR     = PROJECT_ROOT / "data"
NAS_OUTPUT   = Path(os.environ.get("NAS_OUTPUT", "/mnt/nas/input/pyl/astro/output"))

INTERMEDIARY_PATTERNS = (
    "i_*.fit", "i_*.fits", "i_conversion.txt",
    "pp_i_*.fit", "pp_i_*.fits",
    "r_pp_i_*.fit", "r_pp_i_*.fits",
    "*.seq", "*.ssf", "*.lst",
)


def migrate_work_dir(work_dir: Path, nas_dest: Path, dry_run: bool) -> dict:
    stats = {"moved": 0, "cleaned": 0, "skipped": 0}

    for fname in ("res.fit", "res_preview.png"):
        src = work_dir / fname
        if not src.exists():
            continue
        dst = nas_dest / fname
        if dst.exists():
            print(f"    SKIP (already on NAS): {dst}")
            stats["skipped"] += 1
            continue
        print(f"    {'[dry] ' if dry_run else ''}MOVE {src} → {dst}")
        if not dry_run:
            nas_dest.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), dst)
        stats["moved"] += 1

    for pattern in INTERMEDIARY_PATTERNS:
        for f in work_dir.glob(pattern):
            print(f"    {'[dry] ' if dry_run else ''}DEL  {f.name}")
            if not dry_run:
                f.unlink(missing_ok=True)
            stats["cleaned"] += 1

    return stats


def main():
    parser = argparse.ArgumentParser(description="Migrate pipeline results to NAS")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without doing it")
    args = parser.parse_args()

    if not DATA_DIR.exists():
        print(f"data/ not found at {DATA_DIR} — nothing to do")
        return

    total_moved = total_cleaned = total_skipped = 0

    # Walk: data/<date>/<target>/<filter>/<exp>s/
    for date_dir in sorted(DATA_DIR.iterdir()):
        if not date_dir.is_dir():
            continue
        for target_dir in sorted(date_dir.iterdir()):
            if not target_dir.is_dir():
                continue
            for filt_dir in sorted(target_dir.iterdir()):
                if not filt_dir.is_dir():
                    continue
                for exp_dir in sorted(filt_dir.iterdir()):
                    if not exp_dir.is_dir():
                        continue

                    # Only process dirs that have a res.fit or intermediaries
                    has_result = (exp_dir / "res.fit").exists()
                    has_intermediaries = any(
                        list(exp_dir.glob(p)) for p in INTERMEDIARY_PATTERNS
                    )
                    if not (has_result or has_intermediaries):
                        continue

                    nas_dest = NAS_OUTPUT / date_dir.name / target_dir.name / filt_dir.name
                    print(f"\n{date_dir.name}/{target_dir.name}/{filt_dir.name}/{exp_dir.name}")

                    s = migrate_work_dir(exp_dir, nas_dest, args.dry_run)
                    total_moved   += s["moved"]
                    total_cleaned += s["cleaned"]
                    total_skipped += s["skipped"]

    print(f"\n{'[DRY RUN] ' if args.dry_run else ''}Done — "
          f"moved {total_moved} result files, "
          f"deleted {total_cleaned} intermediaries, "
          f"skipped {total_skipped} already on NAS")


if __name__ == "__main__":
    main()
