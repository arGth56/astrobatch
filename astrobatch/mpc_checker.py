"""
mpc_checker.py
==============
Query the IMCCE SkyBoT service for solar system objects in an image field
and annotate the stacked preview PNG with labeled markers.

Public API:
  run_mpc_check(work_dir, result_id, conn, jlog) -> list[dict]
"""

from __future__ import annotations

import json
import logging
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import requests
from astropy.io import fits as afits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
import warnings
warnings.filterwarnings("ignore", category=FITSFixedWarning)

if TYPE_CHECKING:
    from .processing_service import JobLogger

log = logging.getLogger(__name__)

# ── IMCCE SkyBoT ─────────────────────────────────────────────────────────────

SKYBOT_URL = "https://ssp.imcce.fr/webservices/skybot/api/conesearch.php"
SKYBOT_TIMEOUT = 20  # seconds

# SkyBoT "Class" values that indicate Near-Earth Objects / NEAs
# Format from v2 API: "Amor", "Apollo", "Aten", "IEO", or legacy codes
NEO_CODES = {"AMOR", "APOLLO", "ATEN", "IEO", "AMO", "APO", "ATE", "NEA", "NEO"}

# ── Visual style ──────────────────────────────────────────────────────────────

STYLE = {
    "asteroid": {"color": "#00e5ff",  "lw": 1.5, "radius_frac": 0.012},
    "neo":      {"color": "#ff9500",  "lw": 2.0, "radius_frac": 0.016},
    "uncertain":{"color": "#ff9500",  "lw": 1.0, "ls": "--"},
}


def _jd_from_fits(hdr) -> float | None:
    """Extract Julian Date from common FITS header keywords."""
    # Direct JD
    for key in ("JD", "JD-OBS", "JD_OBS"):
        if key in hdr:
            return float(hdr[key])
    # DATE-OBS  →  astropy Time
    for key in ("DATE-OBS", "DATE_OBS"):
        if key in hdr:
            try:
                t = Time(hdr[key], format="isot", scale="utc")
                return float(t.jd)
            except Exception:
                pass
    return None


def _wcs_center_and_fov(fits_path: Path):
    """
    Read WCS from the FITS file.
    Returns (ra_deg, dec_deg, fov_radius_deg, wcs_or_None, image_shape, hdr).

    First tries CRVAL1/2 (proper WCS from plate solve).
    Falls back to RA/DEC header + FOCALLEN/XPIXSZ for the FOV.
    Returns None only if no usable coordinates are found at all.
    """
    try:
        with afits.open(str(fits_path)) as hdul:
            hdu = next((h for h in hdul
                        if h.data is not None and h.data.ndim >= 2), hdul[0])
            hdr = hdu.header
            shape = hdu.data.shape[-2:]  # (ny, nx)
    except Exception as exc:
        log.debug("_wcs_center_and_fov open failed: %s", exc)
        return None

    ny, nx = shape

    # ── Try full WCS (plate solved) ───────────────────────────────────────────
    if "CRVAL1" in hdr and "CRVAL2" in hdr:
        try:
            wcs = WCS(hdr, naxis=2)
            ra, dec = float(hdr["CRVAL1"]), float(hdr["CRVAL2"])
            pscale_deg = float(np.sqrt(
                wcs.pixel_scale_matrix[0, 0] ** 2
                + wcs.pixel_scale_matrix[0, 1] ** 2
            ))
            half_diag = pscale_deg * np.sqrt(nx**2 + ny**2) / 2.0
            return ra, dec, half_diag, wcs, shape, hdr
        except Exception:
            pass  # fall through to header-based fallback

    # ── Fallback: RA/DEC header + instrument keywords ─────────────────────────
    ra  = hdr.get("RA")  or hdr.get("OBJCTRA")
    dec = hdr.get("DEC") or hdr.get("OBJCTDEC")
    if ra is None or dec is None:
        return None

    # Convert sexagesimal strings (OBJCTRA/OBJCTDEC) to decimal if needed
    try:
        ra_deg  = _header_to_deg(ra, is_ra=True)
        dec_deg = _header_to_deg(dec, is_ra=False)
    except Exception:
        return None

    # Estimate pixel scale from FOCALLEN (mm) + XPIXSZ (µm)
    fl_mm  = hdr.get("FOCALLEN")
    px_um  = hdr.get("XPIXSZ") or hdr.get("PIXSCALE")
    xbin   = max(1, int(hdr.get("XBINNING", 1) or 1))
    if fl_mm and px_um:
        pscale_arcsec = 206265.0 * (float(px_um) * 1e-3 * xbin) / float(fl_mm)
    else:
        pscale_arcsec = 1.034  # typical 750mm fl + 3.76µm pixel

    half_diag = pscale_arcsec / 3600.0 * np.sqrt(nx**2 + ny**2) / 2.0

    # Build a simple TAN WCS from the header info (CRPIX = image center)
    try:
        from astropy.wcs import WCS as _WCS
        w = _WCS(naxis=2)
        w.wcs.crpix = [nx / 2.0, ny / 2.0]
        w.wcs.crval = [ra_deg, dec_deg]
        w.wcs.cdelt = [-pscale_arcsec / 3600.0, pscale_arcsec / 3600.0]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.set()
        approx_wcs = w
    except Exception:
        approx_wcs = None

    return ra_deg, dec_deg, half_diag, approx_wcs, shape, hdr


def _header_to_deg(v, is_ra: bool) -> float:
    """Convert a FITS header RA/Dec value to decimal degrees.
    If the value is already numeric (decimal degrees), return it directly.
    Otherwise parse as sexagesimal (HMS for RA, DMS for Dec).
    """
    try:
        # If it's already a float/int, it's decimal degrees
        return float(v)
    except (TypeError, ValueError):
        pass
    # Must be a string in sexagesimal format
    return _sexagesimal_to_deg(str(v), is_ra=is_ra)


def _sexagesimal_to_deg(s: str, is_ra: bool) -> float:
    """Convert '12 14 55' or '12:14:55' (RA) or '+63 46 15' (Dec) to decimal degrees."""
    s = s.strip().replace(":", " ")
    sign = -1.0 if s.startswith("-") else 1.0
    parts = s.lstrip("+-").split()
    val = float(parts[0]) + float(parts[1]) / 60.0 + float(parts[2]) / 3600.0 if len(parts) >= 3 else float(parts[0])
    val *= sign
    return val * 15.0 if is_ra else val


def query_skybot(ra_deg: float, dec_deg: float,
                 sr_deg: float, epoch_jd: float) -> list[dict]:
    """
    Query IMCCE SkyBoT (ssp.imcce.fr v2) for solar system objects in a cone.

    Returns a list of dicts with keys:
      name, number, ra, dec, class, mag, uncertainty_arcsec
    """
    params = {
        "EPOCH":          f"{epoch_jd:.6f}",
        "RA":             f"{ra_deg:.6f}",
        "DEC":            f"{dec_deg:.6f}",
        "SR":             f"{sr_deg:.4f}",
        "RESPONSEFORMAT": "json",
        "OBSERVER":       "500",   # geocentric
    }
    try:
        resp = requests.get(SKYBOT_URL, params=params, timeout=SKYBOT_TIMEOUT)
        resp.raise_for_status()
    except Exception as exc:
        log.warning("SkyBoT query failed: %s", exc)
        return []

    try:
        data = resp.json()
    except Exception:
        return []

    # v2 API returns a JSON array directly
    if isinstance(data, list):
        raw_objects = data
    elif isinstance(data, dict):
        raw_objects = data.get("data", {}).get("sky-bodies", [])
        if not raw_objects:
            raw_objects = data.get("results", [])
    else:
        return []

    objects: list[dict] = []
    for obj in raw_objects:
        try:
            # v2 API keys: "Name", "RA (hms)", "DEC (dms)", "Class", "VMag (mag)", "Err (arcsec)"
            ra_val  = _parse_ra(obj.get("RA (hms)",  obj.get("RA",  obj.get("ra",  None))))
            dec_val = _parse_dec(obj.get("DEC (dms)", obj.get("DEC", obj.get("dec", None))))
            if ra_val is None or dec_val is None:
                continue
            objects.append({
                "name":               str(obj.get("Name", obj.get("name", "?"))),
                "number":             str(obj.get("Num", obj.get("Number", obj.get("number", "")))),
                "ra":                 ra_val,
                "dec":                dec_val,
                "class":              str(obj.get("Class", obj.get("SYSOBJTYPE", obj.get("class", "MBA")))).strip(),
                "mag":                _safe_float(obj.get("VMag (mag)", obj.get("V", obj.get("Mv", obj.get("mag"))))),
                "uncertainty_arcsec": _safe_float(obj.get("Err (arcsec)", obj.get("errRA", obj.get("d", obj.get("uncertainty"))))),
            })
        except Exception:
            continue

    return objects


def _parse_ra(v) -> float | None:
    """Parse RA from decimal degrees or HMS string '11 59 44.748'."""
    if v is None:
        return None
    try:
        return float(v)  # already decimal
    except (TypeError, ValueError):
        pass
    try:
        from astropy.coordinates import Angle
        import astropy.units as u
        return float(Angle(str(v), unit=u.hourangle).deg)
    except Exception:
        return None


def _parse_dec(v) -> float | None:
    """Parse Dec from decimal degrees or DMS string '+05 00 12.039'."""
    if v is None:
        return None
    try:
        return float(v)  # already decimal
    except (TypeError, ValueError):
        pass
    try:
        from astropy.coordinates import Angle
        import astropy.units as u
        return float(Angle(str(v), unit=u.deg).deg)
    except Exception:
        return None


def _safe_float(v) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def classify_neo(obj_class: str) -> bool:
    """Return True if the SkyBoT class corresponds to a NEO/NEA.

    v2 API uses human-readable names: 'Amor', 'Apollo', 'Aten', 'IEO'.
    Legacy codes: AMO, APO, ATE, IEO.
    """
    cls_upper = obj_class.upper().replace(" ", "").replace(">", "")
    return any(code in cls_upper for code in NEO_CODES)


def annotate_preview(
    res_fit: Path,
    preview_png: Path,
    objects: list[dict],
    output_png: Path,
    image_shape: tuple[int, int],
    wcs: WCS,
) -> bool:
    """
    Draw labeled circles on the preview PNG for each solar system object.
    Saves the annotated image to output_png.
    Returns True on success.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from PIL import Image

    try:
        img = Image.open(str(preview_png)).convert("RGB")
    except Exception as exc:
        log.warning("annotate_preview: cannot open preview %s: %s", preview_png, exc)
        return False

    # Preview may be downscaled vs. the FITS; compute the scale factor
    px_w, px_h = img.size           # preview pixel dims
    fits_h, fits_w = image_shape    # FITS pixel dims
    scale_x = px_w / fits_w
    scale_y = px_h / fits_h

    dpi = 96
    fig_w = px_w / dpi
    fig_h = px_h / dpi
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h), dpi=dpi)
    fig.subplots_adjust(0, 0, 1, 1)
    ax.axis("off")

    ax.imshow(np.array(img), origin="upper")

    radius_px = max(12, int(min(px_w, px_h) * STYLE["asteroid"]["radius_frac"]))

    for obj in objects:
        is_neo = classify_neo(obj["class"])
        style  = STYLE["neo"] if is_neo else STYLE["asteroid"]
        color  = style["color"]

        # World → FITS pixel → preview pixel
        try:
            fits_x, fits_y = wcs.world_to_pixel_values(obj["ra"], obj["dec"])
            # FITS pixel origin is (0,0) bottom-left; preview is flipped vertically
            px_x = fits_x * scale_x
            px_y = (fits_h - fits_y) * scale_y
        except Exception:
            continue

        # Skip if far off the image
        if not (
            -radius_px * 3 <= px_x <= px_w + radius_px * 3
            and -radius_px * 3 <= px_y <= px_h + radius_px * 3
        ):
            continue

        r = radius_px * (1.3 if is_neo else 1.0)
        circle = mpatches.Circle(
            (px_x, px_y), r,
            linewidth=style["lw"],
            edgecolor=color,
            facecolor="none",
            zorder=5,
        )
        ax.add_patch(circle)

        # Label: name + "(NEO)" suffix for near-earth objects
        label_parts = [obj["name"]]
        if is_neo:
            label_parts.append("(NEO)")
        if obj.get("mag") is not None:
            label_parts.append(f"V={obj['mag']:.1f}")
        label = " ".join(label_parts)

        ax.annotate(
            label,
            xy=(px_x, px_y),
            xytext=(r + 4, -r - 4),
            textcoords="offset points",
            fontsize=max(7, int(min(px_w, px_h) / 100)),
            color=color,
            fontweight="bold" if is_neo else "normal",
            zorder=6,
            bbox=dict(
                boxstyle="round,pad=0.2",
                facecolor="#000000",
                edgecolor=color,
                alpha=0.55,
                linewidth=0.8,
            ),
        )

        # Uncertainty circle for NEOs with large positional errors
        unc = obj.get("uncertainty_arcsec")
        if is_neo and unc and unc > 30:
            try:
                pscale_arcsec = abs(wcs.pixel_scale_matrix[0, 0]) * 3600.0
                unc_px = (unc / pscale_arcsec) * scale_x
                unc_circle = mpatches.Circle(
                    (px_x, px_y), unc_px,
                    linewidth=STYLE["uncertain"]["lw"],
                    edgecolor=STYLE["uncertain"]["color"],
                    facecolor="none",
                    linestyle=STYLE["uncertain"]["ls"],
                    alpha=0.5,
                    zorder=4,
                )
                ax.add_patch(unc_circle)
            except Exception:
                pass

    # Legend
    legend_patches = [
        mpatches.Patch(edgecolor=STYLE["asteroid"]["color"], facecolor="none",
                       linewidth=1.5, label="Asteroid"),
        mpatches.Patch(edgecolor=STYLE["neo"]["color"], facecolor="none",
                       linewidth=2.0, label="NEO / NEA"),
    ]
    ax.legend(
        handles=legend_patches,
        loc="lower left",
        fontsize=max(6, int(min(px_w, px_h) / 120)),
        facecolor="#000000cc",
        edgecolor="#444",
        labelcolor="white",
        framealpha=0.7,
    )

    ax.set_xlim(0, px_w)
    ax.set_ylim(px_h, 0)

    try:
        fig.savefig(str(output_png), dpi=dpi, bbox_inches="tight", pad_inches=0)
    except Exception as exc:
        log.warning("annotate_preview: save failed: %s", exc)
        plt.close(fig)
        return False

    plt.close(fig)
    return True


def run_mpc_check(
    work_dir: Path,
    res_fit: Path,
    result_id: int,
    db_path: Path,
    jlog,
) -> list[dict]:
    """
    Top-level function called from the pipeline after plate_solve().

    1. Reads WCS from res.fit
    2. Queries SkyBoT for solar system objects in the field
    3. Annotates res_preview.png → res_preview_mpc.png
    4. Stores the object list as JSON in pipeline_results.mpc_objects
       and the preview path in pipeline_results.mpc_preview

    Returns the list of found objects (may be empty).
    Always non-fatal: any error is logged and the function returns [].
    """
    jlog.info("Step MPC: Checking for solar system objects (SkyBoT)...")

    wcs_info = _wcs_center_and_fov(res_fit)
    if wcs_info is None:
        jlog.info("  MPC check skipped — no WCS in res.fit")
        return []

    ra, dec, half_diag_deg, wcs, shape, hdr = wcs_info
    epoch_jd = _jd_from_fits(hdr)
    if epoch_jd is None:
        jlog.info("  MPC check skipped — no DATE-OBS in FITS header")
        return []

    # Search radius: field diagonal + 10% margin so NEO uncertainty circles
    # that overlap the FOV edge are also returned
    sr = half_diag_deg * 1.1
    epoch_dt = Time(epoch_jd, format="jd", scale="utc").to_datetime(timezone=timezone.utc)
    jlog.info(
        f"  Field center RA={ra:.4f}° Dec={dec:.4f}°  "
        f"radius={sr*60:.2f}′  epoch={epoch_dt.strftime('%Y-%m-%d %H:%M UTC')}"
    )

    objects = query_skybot(ra, dec, sr, epoch_jd)

    if not objects:
        jlog.info("  SkyBoT: no solar system objects in field")
        _store_mpc(result_id, [], None, db_path)
        return []

    n_neo = sum(1 for o in objects if classify_neo(o["class"]))
    jlog.info(f"  SkyBoT: {len(objects)} object(s) found  ({n_neo} NEO/NEA)")
    for o in objects:
        neo_tag = " [NEO]" if classify_neo(o["class"]) else ""
        mag_str = f"  V={o['mag']:.1f}" if o.get("mag") is not None else ""
        unc_str = f"  ±{o['uncertainty_arcsec']:.0f}\"" if o.get("uncertainty_arcsec") else ""
        jlog.info(f"    {o['name']}{neo_tag}  cls={o['class']}{mag_str}{unc_str}")

    # Annotate preview
    preview_png = work_dir / "res_preview.png"
    mpc_png     = work_dir / "res_preview_mpc.png"
    mpc_png_path: str | None = None

    if preview_png.exists():
        ok = annotate_preview(res_fit, preview_png, objects, mpc_png, shape, wcs)
        if ok:
            mpc_png_path = str(mpc_png)
            jlog.info(f"  MPC preview saved: {mpc_png.name}")
        else:
            jlog.info("  MPC preview annotation failed (non-fatal)")
    else:
        jlog.info("  res_preview.png not found — skipping annotation")

    _store_mpc(result_id, objects, mpc_png_path, db_path)
    return objects


def _store_mpc(
    result_id: int,
    objects: list[dict],
    mpc_preview_path: str | None,
    db_path: Path,
) -> None:
    """Persist MPC results into pipeline_results."""
    try:
        with sqlite3.connect(str(db_path), timeout=10) as conn:
            conn.execute(
                "UPDATE pipeline_results SET mpc_objects=?, mpc_preview=? WHERE id=?",
                (json.dumps(objects), mpc_preview_path, result_id),
            )
    except Exception as exc:
        log.warning("_store_mpc: DB update failed: %s", exc)
