#!/usr/bin/env python3
"""
Build a per-target imagette in finding-chart style.

For each real target:
  - Find the latest date that has a res_preview.png for each filter
  - Crop each panel to the center 80%
  - Draw a standard research-paper reticle at the target position
    (4 tick marks around center, gap so ticks don't cover the target)
  - Add N/E compass and an arcminute scale bar
  - Stack filters vertically: BP / G / RP
  - Save to OUTPUT_DIR/<target>.png

Plate scale is derived from FITS header (FOCALLEN + XPIXSZ).
Target is always at image center (NINA always points at the target).
WCS is not required.

Usage:
    python3 tools/make_imagettes.py
"""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from collections import defaultdict
from astropy.io import fits as afits   # venv has astropy
import re, math

OUTPUT_ROOT = Path("/mnt/nas/input/pyl/astro/output")
OUTPUT_DIR  = OUTPUT_ROOT / "imagettes"

FILTER_ORDER    = ["G", "BP", "RP"]   # preference order for single-filter imagette
SKIP_RE         = re.compile(r"^(SN2026test|too-fermi-|2026fvx$)", re.I)

CROP_FRACTION   = 0.80   # keep central 80% of each axis
PANEL_W         = 600    # final width of each panel in px
HEADER_H        = 56     # title bar height
BAND_H          = 28     # per-filter label bar height
GAP             = 6      # gap between stacked panels

BG              = (12, 12, 20)
HEADER_BG       = (25, 25, 45)
BAND_BG         = (20, 20, 35)
TEXT_C          = (230, 230, 230)
DIM_C           = (140, 140, 160)
RETICLE_C       = (255, 255, 255)     # tick marks
SHADOW_C        = (0, 0, 0)
COMPASS_C       = (255, 220, 80)      # compass arrows
SCALEBAR_C      = (255, 255, 255)
ACCENT = {"BP": (100, 150, 255), "G": (110, 230, 110), "RP": (255, 110, 100)}

# Plate scale from NINA hardware (fallback if header missing)
DEFAULT_FOCALLEN_MM = 750.0
DEFAULT_PIXSZ_UM    = 3.76


def load_font(size: int, bold: bool = False):
    candidates = []
    if bold:
        candidates += [
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
            "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",
        ]
    candidates += [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
    ]
    for name in candidates:
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            pass
    return ImageFont.load_default()


def crop_center(img: Image.Image, frac: float) -> Image.Image:
    w, h = img.size
    nw, nh = int(w * frac), int(h * frac)
    left = (w - nw) // 2
    top  = (h - nh) // 2
    return img.crop((left, top, left + nw, top + nh))


def fit_to_width(img: Image.Image, target_w: int) -> Image.Image:
    w, h = img.size
    if w == target_w:
        return img
    nh = int(h * target_w / w)
    return img.resize((target_w, nh), Image.LANCZOS)


def plate_scale_arcsec_per_px(fits_path: Path) -> float:
    """
    Derive arcsec/pixel from FITS header.
    Returns the plate scale for the original (uncropped) pixel size.
    """
    try:
        hdr = afits.getheader(str(fits_path))
        fl  = float(hdr.get("FOCALLEN", DEFAULT_FOCALLEN_MM))   # mm
        ps  = float(hdr.get("XPIXSZ",   DEFAULT_PIXSZ_UM))     # μm
        if fl > 0 and ps > 0:
            return 206265.0 * (ps * 1e-3) / fl   # arcsec/px
    except Exception:
        pass
    return 206265.0 * (DEFAULT_PIXSZ_UM * 1e-3) / DEFAULT_FOCALLEN_MM


def get_target_coords(fits_path: Path) -> tuple[str, str] | tuple[None, None]:
    """
    Return (ra_str, dec_str) from FITS header, e.g. ("12h14m53s", "+63°47′14″").
    """
    try:
        hdr = afits.getheader(str(fits_path))
        ra  = hdr.get("OBJCTRA")   # e.g. "12 14 53"
        dec = hdr.get("OBJCTDEC")  # e.g. "+63 47 14"
        if ra and dec:
            rp = ra.strip().split()
            dp = dec.strip().split()
            ra_s  = f"{rp[0]}h{rp[1]}m{rp[2]}s" if len(rp) == 3 else ra.strip()
            dec_s = f"{dp[0]}°{dp[1]}′{dp[2]}″"  if len(dp) == 3 else dec.strip()
            return ra_s, dec_s
    except Exception:
        pass
    return None, None


def draw_reticle(draw: ImageDraw.ImageDraw, cx: int, cy: int, scale: float):
    """
    Draw 4 tick marks around (cx, cy).
    scale: display pixels per arcsec (for sizing the gap and tick length).
    gap = 15 arcsec, tick_len = 20 arcsec (reasonable for a typical SNe field).
    """
    gap = max(18, int(15 * scale))      # gap from center in display px
    tlen = max(14, int(20 * scale))     # tick length in display px
    lw = 2

    def line_with_shadow(xy0, xy1, width=lw):
        x0, y0, x1, y1 = *xy0, *xy1
        for dx, dy in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
            draw.line([(x0+dx, y0+dy), (x1+dx, y1+dy)], fill=SHADOW_C, width=width+2)
        draw.line([(x0, y0), (x1, y1)], fill=RETICLE_C, width=width)

    # Top tick
    line_with_shadow((cx, cy - gap - tlen), (cx, cy - gap))
    # Bottom tick
    line_with_shadow((cx, cy + gap), (cx, cy + gap + tlen))
    # Left tick
    line_with_shadow((cx - gap - tlen, cy), (cx - gap, cy))
    # Right tick
    line_with_shadow((cx + gap, cy), (cx + gap + tlen, cy))


def draw_compass(draw: ImageDraw.ImageDraw, x: int, y: int, font):
    """
    Draw N↑ E← compass at position (x, y) = top-left anchor.
    North is up, East is left (standard astronomical orientation).
    """
    arm = 28
    lw  = 2
    tip_n = (x + arm // 2, y)
    tip_e = (x, y + arm // 2)
    base  = (x + arm // 2, y + arm // 2)

    def arrow(p0, p1, label, lx, ly):
        for dx, dy in [(-1,-1),(1,-1),(-1,1),(1,1)]:
            draw.line([(p0[0]+dx, p0[1]+dy), (p1[0]+dx, p1[1]+dy)], fill=SHADOW_C, width=lw+2)
        draw.line([p0, p1], fill=COMPASS_C, width=lw)
        # Arrowhead (3 short lines)
        dx_full = p1[0] - p0[0]
        dy_full = p1[1] - p0[1]
        ln = math.hypot(dx_full, dy_full)
        if ln == 0:
            return
        nx, ny = dx_full / ln, dy_full / ln
        px, py = -ny, nx
        tip = p1
        ah = 6
        aw = 4
        left  = (int(tip[0] - nx * ah + px * aw), int(tip[1] - ny * ah + py * aw))
        right = (int(tip[0] - nx * ah - px * aw), int(tip[1] - ny * ah - py * aw))
        draw.polygon([tip, left, right], fill=COMPASS_C)
        draw.text((lx, ly), label, fill=COMPASS_C, font=font, anchor="mm")

    arrow(base, tip_n, "N", tip_n[0], tip_n[1] - 9)
    arrow(base, tip_e, "E", tip_e[0] - 9, tip_e[1])


def draw_scalebar(draw: ImageDraw.ImageDraw, x: int, y: int,
                  bar_px: int, label: str, font):
    """
    Draw a horizontal scale bar at (x, y) = right edge, centered vertically.
    """
    lw  = 3
    tick = 5
    x0, x1 = x - bar_px, x
    for dx, dy in [(-1,-1),(1,-1),(-1,1),(1,1)]:
        draw.line([(x0+dx, y+dy), (x1+dx, y+dy)], fill=SHADOW_C, width=lw+2)
        draw.line([(x0+dx, y-tick+dy), (x0+dx, y+tick+dy)], fill=SHADOW_C, width=lw)
        draw.line([(x1+dx, y-tick+dy), (x1+dx, y+tick+dy)], fill=SHADOW_C, width=lw)
    draw.line([(x0, y), (x1, y)], fill=SCALEBAR_C, width=lw)
    draw.line([(x0, y - tick), (x0, y + tick)], fill=SCALEBAR_C, width=lw)
    draw.line([(x1, y - tick), (x1, y + tick)], fill=SCALEBAR_C, width=lw)
    draw.text(((x0 + x1) // 2, y - tick - 6), label,
              fill=SCALEBAR_C, font=font, anchor="mb")


def discover() -> dict[str, dict[str, tuple[str, Path]]]:
    """
    Returns { target: { filter: (date_str, png_path) } }  — latest date per filter.
    Also collects the corresponding res.fit path.
    """
    result: dict[str, dict[str, tuple[str, Path, Path]]] = defaultdict(dict)
    for png in OUTPUT_ROOT.glob("*/*/*/res_preview.png"):
        parts   = png.parts
        date_str = parts[-4]
        target   = parts[-3]
        filt     = parts[-2]
        if SKIP_RE.match(target):
            continue
        fit = png.parent / "res.fit"
        existing = result[target].get(filt)
        if existing is None or date_str > existing[0]:
            result[target][filt] = (date_str, png, fit)
    return result


def make_imagette(target: str,
                  filters: dict[str, tuple[str, Path, Path]],
                  out_path: Path) -> None:

    # Pick the single best filter: G > BP > RP > first available
    filt = next((f for f in FILTER_ORDER if f in filters), next(iter(filters)))
    date_str, png_path, fit_path = filters[filt]

    try:
        img = Image.open(str(png_path)).convert("RGB")
    except Exception as e:
        print(f"  ⚠ Cannot open {png_path}: {e}")
        return

    arcsec_per_orig_px = plate_scale_arcsec_per_px(fit_path)
    orig_w = img.width
    arcsec_per_display_px = arcsec_per_orig_px * orig_w / PANEL_W

    img = crop_center(img, CROP_FRACTION)
    img = fit_to_width(img, PANEL_W)
    pw, ph = img.size

    ra_str, dec_str = get_target_coords(fit_path)

    panels = [{
        "filt":  filt,
        "date":  date_str,
        "img":   img,
        "scale": arcsec_per_display_px,
        "ra":    ra_str,
        "dec":   dec_str,
    }]

    if not panels:
        print(f"  ✗ No usable panels for {target}")
        return

    font_title  = load_font(22, bold=True)
    font_coord  = load_font(14)
    font_band   = load_font(14, bold=True)
    font_annot  = load_font(12)

    pw = panels[0]["img"].width
    ph = panels[0]["img"].height

    total_w = pw
    total_h = HEADER_H + len(panels) * (ph + BAND_H) + (len(panels) - 1) * GAP

    canvas = Image.new("RGB", (total_w, total_h), BG)
    draw   = ImageDraw.Draw(canvas)

    # ── Header bar ───────────────────────────────────────────────────────────────
    draw.rectangle([(0, 0), (total_w, HEADER_H)], fill=HEADER_BG)
    draw.text((total_w // 2, HEADER_H // 2 - 8),
              target.upper().replace("_", " "),
              fill=TEXT_C, font=font_title, anchor="mm")
    coord_line = ""
    if panels[0]["ra"]:
        coord_line = f"α {panels[0]['ra']}   δ {panels[0]['dec']}"
    if coord_line:
        draw.text((total_w // 2, HEADER_H - 12),
                  coord_line, fill=DIM_C, font=font_coord, anchor="mm")

    # ── Panels (stacked vertically) ───────────────────────────────────────────
    y_cursor = HEADER_H
    for idx, p in enumerate(panels):
        img_y = y_cursor

        # Paste image
        canvas.paste(p["img"], (0, img_y))

        # ── Reticle at center of panel ────────────────────────────────────────
        cx, cy = pw // 2, img_y + ph // 2
        scale_px_per_arcsec = 1.0 / p["scale"] if p["scale"] > 0 else 0.5
        draw_reticle(draw, cx, cy, scale_px_per_arcsec)

        # ── Compass (first panel only, top-left corner) ───────────────────────
        if idx == 0:
            draw_compass(draw, 14, img_y + 14, font_annot)

        # ── Scale bar (last panel only, bottom-right corner) ──────────────────
        if idx == len(panels) - 1:
            scale = p["scale"]     # arcsec per display px
            # pick a round arcminute label that fills ~15% of panel width
            for arcmin in [1, 2, 3, 5]:
                bar_px = int((arcmin * 60) / scale)
                if 0.10 * pw <= bar_px <= 0.25 * pw:
                    break
            if bar_px < 10:
                bar_px = 30
                arcmin = round(bar_px * scale / 60, 1)
            bar_x = pw - 14
            bar_y = img_y + ph - 22
            draw_scalebar(draw, bar_x, bar_y, bar_px, f"{arcmin}′", font_annot)

        # ── Band label bar ────────────────────────────────────────────────────
        band_y = img_y + ph
        draw.rectangle([(0, band_y), (total_w, band_y + BAND_H)], fill=BAND_BG)
        accent = ACCENT.get(p["filt"], DIM_C)
        # Colored left stripe
        draw.rectangle([(0, band_y), (4, band_y + BAND_H)], fill=accent)
        draw.text((16, band_y + BAND_H // 2),
                  p["filt"], fill=accent, font=font_band, anchor="lm")
        draw.text((pw - 8, band_y + BAND_H // 2),
                  p["date"], fill=DIM_C, font=font_annot, anchor="rm")

        y_cursor = band_y + BAND_H + GAP

    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(str(out_path), format="PNG", optimize=True)
    filters_str = " / ".join(p["filt"] for p in panels)
    print(f"  ✓  {target:30s}  [{filters_str}]  → {out_path.name}")


def main():
    targets = discover()
    print(f"Found {len(targets)} targets → {OUTPUT_DIR}/\n")
    for target, filters in sorted(targets.items()):
        out_path = OUTPUT_DIR / f"{target}.png"
        make_imagette(target, filters, out_path)
    print(f"\nDone — {len(targets)} imagettes in {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
