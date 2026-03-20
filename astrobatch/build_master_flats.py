#!/usr/bin/env python3
"""
build_master_flats.py

Scans NAS flat frames, groups by filter (from FITS header), stacks each group
into a master flat using Siril, and writes results to calib/flats/.

Filter name → master flat filename mapping (edit FILTER_MAP to customise):
  BP       → flat_Gbp.fit   (Gaia BP convention used in astrobatch)
  RP       → flat_Grp.fit
  G        → flat_G.fit
  Filter 5 → flat_Filter5.fit  (sanitised)

Usage:
    # From repo root, with venv active:
    python -m astrobatch.build_master_flats

    # Override NAS path:
    python -m astrobatch.build_master_flats --flat-input /custom/path

    # Dry run (show what would be done):
    python -m astrobatch.build_master_flats --dry-run
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from astropy.io import fits

# ── Paths ─────────────────────────────────────────────────────────────────────

PROJECT_ROOT    = Path(__file__).resolve().parent.parent
DEFAULT_INPUT   = Path("/mnt/nas/input/pyl/astro/input/FLAT")
CALIB_FLATS_DIR = PROJECT_ROOT / "calib" / "flats"
SIRIL_BIN       = "siril"

# ── Filter name → master flat stem ────────────────────────────────────────────
# Keys are UPPERCASE. Edit to match your naming convention.

FILTER_MAP: dict[str, str] = {
    "BP":       "Gbp",
    "RP":       "Grp",
    "G":        "G",
    "FILTER 5": "Filter5",
    "FILTER5":  "Filter5",
    "L":        "L",
    "R":        "R",
    "V":        "V",
    "B":        "B",
    "HA":       "Ha",
    "OIII":     "OIII",
    "SII":      "SII",
}


def sanitize(name: str) -> str:
    """Make a filter name safe for use as a filename."""
    return re.sub(r"[^A-Za-z0-9_-]", "", name.replace(" ", "_"))


def master_flat_name(filter_name: str) -> str:
    """Return the master flat filename stem for a given filter name."""
    key = filter_name.strip().upper()
    if key in FILTER_MAP:
        return f"flat_{FILTER_MAP[key]}.fit"
    return f"flat_{sanitize(filter_name)}.fit"


def get_filter(fits_path: Path) -> str | None:
    """Read the FILTER keyword from a FITS header."""
    try:
        return fits.getval(str(fits_path), "FILTER", ext=0)
    except Exception:
        return None


def scan_flats(input_dir: Path) -> dict[str, list[Path]]:
    """Group flat FITS files by filter name."""
    groups: dict[str, list[Path]] = defaultdict(list)
    for ext in ("*.fit", "*.fits", "*.FIT", "*.FITS"):
        for f in input_dir.glob(ext):
            filt = get_filter(f)
            if filt:
                groups[filt].append(f)
            else:
                print(f"  ⚠️  No FILTER keyword in {f.name} — skipping")
    return dict(groups)


def build_master_flat(
    filter_name: str,
    flat_files: list[Path],
    output_path: Path,
    dry_run: bool = False,
) -> bool:
    """Stack flat files with Siril and save master flat to output_path."""

    print(f"\n  Filter: {filter_name!r}  ({len(flat_files)} frames)")
    print(f"  Output: {output_path}")

    if dry_run:
        print("  [dry-run] skipping Siril execution")
        return True

    with tempfile.TemporaryDirectory(prefix=f"flat_{sanitize(filter_name)}_") as _tmp:
        tmpdir = Path(_tmp)

        # Copy flat frames into temp dir, renamed sequentially
        print(f"  Copying {len(flat_files)} files to temp dir...")
        for i, src in enumerate(sorted(flat_files)):
            shutil.copy2(src, tmpdir / f"flat_{i:04d}.fits")

        # Write Siril stacking script
        script = tmpdir / "stack_flats.ssf"
        script.write_text(
            f"requires 1.2.0\n"
            f"cd {tmpdir}\n"
            "convert flat\n"
            "stack flat rej w 3 3 -norm=mul\n"
        )

        # Run Siril
        print("  Running Siril...")
        result = subprocess.run(
            [SIRIL_BIN, "-s", str(script)],
            capture_output=True,
            text=True,
            timeout=300,
        )

        if result.returncode != 0:
            print(f"  ✗ Siril failed (exit {result.returncode})")
            print(result.stderr[-800:] if result.stderr else "(no stderr)")
            return False

        # Siril names the stacked output "flat_stacked.fit"
        candidates = sorted(
            list(tmpdir.glob("flat_stacked.fit"))
            + list(tmpdir.glob("*stacked*.fit"))
            + list(tmpdir.glob("stacked*.fit"))
        )

        if not candidates:
            print("  ✗ Could not find Siril output file")
            print("  Files in tmpdir:", [p.name for p in tmpdir.iterdir()])
            return False

        stacked = candidates[0]
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Back up existing master flat
        if output_path.exists():
            backup = output_path.with_suffix(".fit.bak")
            shutil.copy2(output_path, backup)
            print(f"  Backed up old master → {backup.name}")

        shutil.copy2(stacked, output_path)
        size_mb = output_path.stat().st_size / 1_048_576
        print(f"  ✓ Saved {output_path.name} ({size_mb:.1f} MB)")
        return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build master flat frames from NAS input")
    parser.add_argument(
        "--flat-input",
        default=str(DEFAULT_INPUT),
        help=f"Directory containing raw flat FITS files (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--calib-dir",
        default=str(CALIB_FLATS_DIR),
        help=f"Output directory for master flats (default: {CALIB_FLATS_DIR})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without running Siril",
    )
    parser.add_argument(
        "--filter",
        metavar="NAME",
        help="Process only this filter (case-insensitive)",
    )
    args = parser.parse_args(argv)

    flat_input  = Path(args.flat_input)
    calib_dir   = Path(args.calib_dir)

    if not flat_input.exists():
        print(f"✗ Flat input directory not found: {flat_input}")
        return 1

    print(f"Scanning flats in: {flat_input}")
    groups = scan_flats(flat_input)

    if not groups:
        print("No flat frames found.")
        return 0

    print(f"\nFound {sum(len(v) for v in groups.values())} flat files across {len(groups)} filters:")
    for filt, files in sorted(groups.items()):
        out = calib_dir / master_flat_name(filt)
        status = "✓ exists" if out.exists() else "new"
        print(f"  {filt:<15} {len(files):>3} frames  →  {out.name}  [{status}]")

    # Filter to single filter if requested
    if args.filter:
        target = args.filter.strip().upper()
        groups = {k: v for k, v in groups.items() if k.upper() == target}
        if not groups:
            print(f"\n✗ No frames found for filter: {args.filter!r}")
            return 1

    print()
    ok = 0
    fail = 0
    for filt, files in sorted(groups.items()):
        out = calib_dir / master_flat_name(filt)
        if build_master_flat(filt, files, out, dry_run=args.dry_run):
            ok += 1
        else:
            fail += 1

    print(f"\n{'='*50}")
    print(f"Done: {ok} master flat(s) built successfully, {fail} failed.")

    if ok and not args.dry_run:
        print(f"\nMaster flats in {calib_dir}:")
        for f in sorted(calib_dir.glob("flat_*.fit")):
            print(f"  {f.name}  ({f.stat().st_size/1_048_576:.1f} MB)")

    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
