#!/usr/bin/env python3
"""
make_time_animation.py — Perfectly aligned time-lapse GIF.

Two-pass alignment:
  Pass 1: astroalign (handles rotation, translation, scale via star pattern matching)
  Pass 2: sub-pixel refinement via phase cross-correlation on aligned result
  
High-quality Lanczos interpolation for rotation/warping.
Histogram matching for consistent brightness.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from skimage.exposure import match_histograms
from skimage.transform import SimilarityTransform, warp
from skimage.registration import phase_cross_correlation
from scipy.ndimage import shift as ndi_shift
import astroalign as aa
import imageio.v3 as iio
import cv2

sys.stdout.reconfigure(line_buffering=True)


def load_fits(path):
    """Load full FITS image."""
    with fits.open(str(path), memmap=True) as hdul:
        data = hdul[0].data.astype(np.float64)
    data = np.nan_to_num(data, nan=0.0)
    return np.clip(data, 0, None)


def find_and_apply_transform(source, target):
    """
    Find transform using astroalign, apply with high-quality interpolation.
    Returns (aligned_image, success).
    """
    configs = [
        dict(detection_sigma=3, max_control_points=80, min_area=5),
        dict(detection_sigma=5, max_control_points=50, min_area=5),
        dict(detection_sigma=8, max_control_points=30, min_area=10),
    ]

    for cfg in configs:
        try:
            transf, (src_pts, dst_pts) = aa.find_transform(source, target, **cfg)
            # Apply with order=3 (bicubic) for better interpolation
            aligned = warp(source, inverse_map=transf.inverse,
                          output_shape=target.shape,
                          order=3, mode='constant', cval=0,
                          preserve_range=True)
            n_stars = len(src_pts)
            rot_deg = np.degrees(transf.rotation)
            return aligned, True, n_stars, rot_deg
        except Exception:
            continue
    return None, False, 0, 0


def subpixel_refine(source, target, upsample=100):
    """
    Sub-pixel alignment refinement using phase cross-correlation.
    Uses central 60% of the image (avoids edge artifacts from rotation).
    """
    h, w = source.shape
    margin = int(0.2 * min(h, w))
    src_crop = source[margin:h-margin, margin:w-margin]
    tgt_crop = target[margin:h-margin, margin:w-margin]

    # Mask out zeros (from rotation fill)
    src_crop = np.where(src_crop > 0, src_crop, np.median(src_crop[src_crop > 0]))
    tgt_crop = np.where(tgt_crop > 0, tgt_crop, np.median(tgt_crop[tgt_crop > 0]))

    shift_detected, error, diffphase = phase_cross_correlation(
        tgt_crop, src_crop, upsample_factor=upsample
    )

    if np.any(np.abs(shift_detected) > 5):
        # Something wrong, skip refinement
        return source, shift_detected

    refined = ndi_shift(source, shift_detected, order=3, mode='constant', cval=0)
    return refined, shift_detected


def crop_center(data, size):
    """Crop center of image."""
    h, w = data.shape
    y0 = (h - size) // 2
    x0 = (w - size) // 2
    return data[y0:y0+size, x0:x0+size]


def normalize_frame(data):
    """Normalize to uint8 using ZScale."""
    valid = data[np.isfinite(data) & (data > 0)]
    if len(valid) < 100:
        return np.zeros(data.shape, dtype=np.uint8)
    interval = ZScaleInterval(contrast=0.15)
    vmin, vmax = interval.get_limits(valid)
    stretched = np.clip((data - vmin) / (vmax - vmin), 0, 1)
    return (stretched * 255).astype(np.uint8)


def add_label(frame, text, filter_name=""):
    """Burn date + filter label."""
    h, w = frame.shape[:2]
    labeled = cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR)
    font = cv2.FONT_HERSHEY_SIMPLEX
    scale = 0.7
    thickness = 2
    cv2.putText(labeled, text, (12, 32), font, scale, (0, 0, 0), thickness+2, cv2.LINE_AA)
    cv2.putText(labeled, text, (10, 30), font, scale, (255, 255, 255), thickness, cv2.LINE_AA)
    if filter_name:
        (tw, _), _ = cv2.getTextSize(filter_name, font, scale, thickness)
        cv2.putText(labeled, filter_name, (w-tw-12, 32), font, scale, (0, 0, 0), thickness+2, cv2.LINE_AA)
        cv2.putText(labeled, filter_name, (w-tw-10, 30), font, scale, (200, 200, 255), thickness, cv2.LINE_AA)
    return cv2.cvtColor(labeled, cv2.COLOR_BGR2GRAY)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="sn_2026fvx")
    parser.add_argument("--base", default="/mnt/nas/input/pyl/astro/output")
    parser.add_argument("--output", default="/var/www/astrobatch/data/animations")
    parser.add_argument("--final-size", type=int, default=800)
    parser.add_argument("--fps", type=float, default=3.0)
    parser.add_argument("--filters", nargs="*", default=["G", "BP", "RP"])
    args = parser.parse_args()

    base = Path(args.base)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Discover frames
    filter_frames = defaultdict(list)
    for date_dir in sorted(base.iterdir()):
        if not date_dir.is_dir():
            continue
        target_dir = date_dir / args.target
        if not target_dir.exists():
            continue
        for filt in args.filters:
            res_path = target_dir / filt / "res.fit"
            if res_path.exists():
                filter_frames[filt].append((date_dir.name, res_path))

    if not filter_frames:
        print(f"No frames found for {args.target}")
        sys.exit(1)

    print(f"Target: {args.target}")
    for filt, frames in sorted(filter_frames.items()):
        print(f"  {filt}: {len(frames)} epochs ({frames[0][0]} -> {frames[-1][0]})")

    # Load FULL reference frame (first G) — use full 3008x3008 for best star matching
    ref_filter = "G" if "G" in filter_frames else list(filter_frames.keys())[0]
    ref_date, ref_path = filter_frames[ref_filter][0]
    ref_full = load_fits(ref_path)
    print(f"\n  Reference: {ref_date} ({ref_filter}), full {ref_full.shape[0]}x{ref_full.shape[1]}")
    print(f"  Final output: {args.final_size}x{args.final_size}px")
    print(f"  Method: astroalign (full frame) + sub-pixel phase cross-correlation\n")

    # Pre-crop reference to output region for sub-pixel refinement
    ref_cropped = crop_center(ref_full, args.final_size)

    for filt in args.filters:
        if filt not in filter_frames:
            continue

        frames = filter_frames[filt]
        print(f"=== {filt} filter: {len(frames)} frames ===")

        aligned_data = []
        dates_used = []

        for i, (date, fpath) in enumerate(frames):
            src_full = load_fits(fpath)

            if i == 0 and filt == ref_filter:
                # Reference frame — just crop
                cropped = crop_center(src_full, args.final_size)
                aligned_data.append(cropped)
                dates_used.append(date)
                print(f"  [{i+1}/{len(frames)}] {date} (reference)")
                continue

            # Pass 1: astroalign on FULL frames for robust transform detection
            aligned_full, ok, n_stars, rot = find_and_apply_transform(src_full, ref_full)
            if not ok:
                print(f"  [{i+1}/{len(frames)}] {date} SKIP (astroalign failed)")
                continue

            # Crop to output region
            aligned_cropped = crop_center(aligned_full, args.final_size)

            # Pass 2: sub-pixel refinement on cropped region
            refined, subpx_shift = subpixel_refine(aligned_cropped, ref_cropped)

            aligned_data.append(refined)
            dates_used.append(date)
            dy, dx = subpx_shift
            print(f"  [{i+1}/{len(frames)}] {date} rot={rot:+.2f}° stars={n_stars} refine=({dx:+.2f},{dy:+.2f})px")

        if len(aligned_data) < 2:
            print(f"  Not enough frames for {filt}")
            continue

        # Histogram match all to first frame
        print(f"  Histogram matching {len(aligned_data)} frames...")
        ref_hist = aligned_data[0]

        gif_frames = []
        for data, date in zip(aligned_data, dates_used):
            mask = (data > 0) & np.isfinite(data)
            if mask.sum() > 100:
                matched = match_histograms(data, ref_hist).astype(np.float64)
                matched[~mask] = 0
            else:
                matched = data
            normalized = normalize_frame(matched)
            labeled = add_label(normalized, date, filt)
            gif_frames.append(labeled)

        # Write GIF
        gif_path = output_dir / f"{args.target}_{filt}_timelapse.gif"
        duration_ms = int(1000 / args.fps)
        iio.imwrite(str(gif_path), gif_frames, duration=duration_ms, loop=0)
        size_mb = gif_path.stat().st_size / 1024 / 1024
        print(f"  -> {gif_path.name} ({len(gif_frames)} frames, {size_mb:.1f} MB)\n")

    print(f"Done! Animations in {output_dir}")


if __name__ == "__main__":
    main()
