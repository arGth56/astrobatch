#!/usr/bin/env python3
"""Aperture-photometry light curve for SN 2026fvx (wrapper around astrobatch.lightcurves)."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astrobatch.lightcurves import DEFAULT_DB, generate_for_target  # noqa: E402

TARGET = "sn_2026fvx"
OUTPUT = PROJECT_ROOT / "doc" / "sn2026fvx_aperture_lc.png"


def main() -> None:
    out = generate_for_target(TARGET, OUTPUT.parent, DEFAULT_DB)
    path = out.get("aperture")
    if path and Path(path) != OUTPUT:
        Path(path).rename(OUTPUT)
        print(f"Saved {OUTPUT}")
    elif path:
        print(f"Saved {path}")
    else:
        raise SystemExit("No aperture photometry found for sn_2026fvx")


if __name__ == "__main__":
    main()
