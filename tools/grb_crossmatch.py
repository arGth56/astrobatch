#!/usr/bin/env python3
"""
GRB 260402A optical counterpart search
=======================================
- Queries Gaia DR3 within 40" of the EP-FXT position
- Scans NAS FITS directory for new frames matching the target
- Source-extracts each frame with sep, matches against Gaia
- Flags any source inside the 20" error circle with no Gaia counterpart
  → candidate optical afterglow

Usage:
    python3 grb_crossmatch.py [--wait]   # --wait polls until frames appear
"""

import sys
import time
import glob
import argparse
import warnings
warnings.filterwarnings("ignore")

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

# ── GRB parameters ─────────────────────────────────────────────────────────
GRB_NAME   = "GRB 260402A"
GRB_RA     = 170.9607   # deg
GRB_DEC    = 63.6128    # deg
ERR_ARCSEC = 20.0       # EP FXT 1-sigma position error
SEARCH_RAD = 40.0       # arcsec — slightly wider than error circle for context

# ── Paths ───────────────────────────────────────────────────────────────────
NAS_ROOT   = "/mnt/nas/input/pyl/astro/input"
TARGET_KEY = "GRB"      # matched against FITS OBJECT header (case-insensitive)

# ── Detection thresholds ────────────────────────────────────────────────────
SEP_THRESH     = 3.0    # sep detection threshold (sigma above background)
GAIA_MATCH_RAD = 3.0    # arcsec — sources within this of a Gaia star are "known"

def query_gaia_dr3(ra, dec, radius_arcsec):
    """Return SkyCoord array of all Gaia DR3 sources within radius (via VizieR)."""
    try:
        from astroquery.vizier import Vizier
        v = Vizier(columns=["RA_ICRS", "DE_ICRS", "Gmag"],
                   row_limit=500, timeout=30)
        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame="icrs")
        result = v.query_region(coord,
                                radius=radius_arcsec * u.arcsec,
                                catalog="I/355/gaiadr3")
        if not result or len(result) == 0 or len(result[0]) == 0:
            return SkyCoord([], [], unit="deg"), []
        tbl = result[0]
        coords = SkyCoord(ra=tbl["RA_ICRS"].data*u.deg,
                          dec=tbl["DE_ICRS"].data*u.deg, frame="icrs")
        mags = list(tbl["Gmag"].data.filled(99.0))
        return coords, mags
    except Exception as e:
        print(f"  [Gaia] query failed: {e}")
        return SkyCoord([], [], unit="deg"), []

def find_grb_fits_files():
    """
    Find FITS files that are GRB 260402A frames.
    NINA puts all frames in DATE/SNAPSHOT/ with no target subfolder.
    We identify GRB frames by:
      1. OBJECT header contains GRB/260402, OR
      2. Filename timestamp >= 22:49 UTC today (sequence start time)
    """
    SEQ_START_UTC = "2026-04-02_22-49"  # sequence started at 22:49 UTC
    candidates = []
    snap_dir = f"{NAS_ROOT}/2026-04-02/SNAPSHOT"
    all_fits = (glob.glob(f"{snap_dir}/*.fits") +
                glob.glob(f"{snap_dir}/*.fit") +
                glob.glob(f"{NAS_ROOT}/**/*.fits", recursive=True) +
                glob.glob(f"{NAS_ROOT}/**/*.fit",  recursive=True))
    for fits_path in all_fits:
        fname = fits_path.split("/")[-1]
        # Match by filename timestamp
        if fname >= SEQ_START_UTC:
            candidates.append(fits_path)
            continue
        # Also check OBJECT header
        try:
            hdr = fits.getheader(fits_path)
            obj = str(hdr.get("OBJECT", "")).upper()
            if TARGET_KEY in obj or "260402" in obj:
                candidates.append(fits_path)
        except Exception:
            pass
    return sorted(set(candidates))

def extract_sources(fits_path):
    """Extract sources from a FITS image using sep. Returns (coords, fluxes)."""
    try:
        import sep
    except ImportError:
        print("  [sep] not installed — run: pip install sep")
        return SkyCoord([], [], unit="deg"), []

    try:
        with fits.open(fits_path) as hdul:
            # Find image HDU
            hdu = next((h for h in hdul if h.data is not None and h.data.ndim == 2), hdul[0])
            data = hdu.data.astype(np.float64)
            wcs  = WCS(hdu.header)

        # Background subtraction
        bkg  = sep.Background(data)
        data_sub = data - bkg

        # Source extraction
        objects = sep.extract(data_sub, SEP_THRESH, err=bkg.globalrms)

        if len(objects) == 0:
            return SkyCoord([], [], unit="deg"), []

        # Pixel → sky
        sky = wcs.pixel_to_world(objects["x"], objects["y"])
        fluxes = objects["flux"].tolist()
        return sky, fluxes

    except Exception as e:
        print(f"  [extract] {fits_path}: {e}")
        return SkyCoord([], [], unit="deg"), []

def analyse_frame(fits_path, grb_coord, gaia_coords, gaia_mags, frame_num):
    """Check one frame for a new source inside the error circle."""
    print(f"\n  Frame {frame_num}: {fits_path.split('/')[-1]}")
    sources, fluxes = extract_sources(fits_path)

    if len(sources) == 0:
        print("    No sources extracted (clouds? bad WCS?)")
        return None

    print(f"    {len(sources)} sources extracted")

    # Separate into: inside error circle / outside
    sep_from_grb = grb_coord.separation(sources).to(u.arcsec).value
    inside_mask  = sep_from_grb <= ERR_ARCSEC

    if not np.any(inside_mask):
        print(f"    No sources within {ERR_ARCSEC}\" of GRB position")
        return None

    inside_sources = sources[inside_mask]
    inside_seps    = sep_from_grb[inside_mask]
    inside_fluxes  = np.array(fluxes)[inside_mask]

    candidates = []
    for src, sep_val, flux in zip(inside_sources, inside_seps, inside_fluxes):
        # Check against Gaia
        if len(gaia_coords) > 0:
            gaia_sep = src.separation(gaia_coords).to(u.arcsec).value
            nearest_gaia = gaia_sep.min()
            nearest_mag  = gaia_mags[gaia_sep.argmin()] if gaia_mags else None
        else:
            nearest_gaia = 9999.0
            nearest_mag  = None

        is_new = nearest_gaia > GAIA_MATCH_RAD

        status = "*** NEW SOURCE ***" if is_new else f"Gaia match ({nearest_gaia:.1f}\" / G={nearest_mag:.1f})"
        print(f"    Source at offset {sep_val:.1f}\" from GRB: {status}")
        print(f"      RA={src.ra.deg:.5f}  Dec={src.dec.deg:.5f}  flux={flux:.0f}")

        if is_new:
            candidates.append({
                "ra": src.ra.deg, "dec": src.dec.deg,
                "sep_arcsec": sep_val, "flux": flux,
                "nearest_gaia_arcsec": nearest_gaia,
            })

    return candidates if candidates else None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wait", action="store_true", help="Poll until frames appear")
    args = parser.parse_args()

    grb_coord = SkyCoord(ra=GRB_RA*u.deg, dec=GRB_DEC*u.deg, frame="icrs")

    print(f"\n{'='*60}")
    print(f"  {GRB_NAME} optical counterpart search")
    print(f"  RA {GRB_RA}°  Dec {GRB_DEC}°  err={ERR_ARCSEC}\"")
    print(f"{'='*60}")

    # Query Gaia once
    print(f"\nQuerying Gaia DR3 within {SEARCH_RAD}\" ...")
    gaia_coords, gaia_mags = query_gaia_dr3(GRB_RA, GRB_DEC, SEARCH_RAD)
    print(f"  {len(gaia_coords)} Gaia sources found")
    if len(gaia_coords) > 0:
        for i, (gc, gm) in enumerate(zip(gaia_coords, gaia_mags)):
            sep_v = grb_coord.separation(gc).to(u.arcsec).value
            flag  = " ← inside error circle" if sep_v <= ERR_ARCSEC else ""
            print(f"  [{i+1}] G={gm:.2f}  sep={sep_v:.1f}\"{flag}")

    # Wait for / find FITS files
    seen = set()
    all_candidates = []
    print(f"\nWatching for FITS frames in {NAS_ROOT} ...")

    poll_count = 0
    while True:
        files = find_grb_fits_files()
        new_files = [f for f in files if f not in seen]

        for i, fp in enumerate(new_files):
            seen.add(fp)
            result = analyse_frame(fp, grb_coord, gaia_coords, gaia_mags, len(seen))
            if result:
                all_candidates.extend(result)

        if all_candidates:
            print(f"\n{'!'*60}")
            print(f"  CANDIDATE OPTICAL COUNTERPART DETECTED")
            print(f"{'!'*60}")
            for c in all_candidates:
                print(f"  RA={c['ra']:.5f}  Dec={c['dec']:.5f}")
                print(f"  Offset from GRB centre: {c['sep_arcsec']:.1f}\"")
                print(f"  Nearest Gaia source: {c['nearest_gaia_arcsec']:.1f}\" away (not in catalog)")
                print(f"  → Consistent with optical afterglow of {GRB_NAME}")

        if not args.wait:
            if not new_files:
                if not files:
                    print("  No GRB frames found yet. Run with --wait to poll continuously.")
            break

        poll_count += 1
        if poll_count % 10 == 0:
            print(f"  [{time.strftime('%H:%M:%S')}] waiting for frames ({len(seen)} processed so far)...")
        time.sleep(30)

if __name__ == "__main__":
    main()
