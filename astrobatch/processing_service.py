#!/usr/bin/env python3
"""
processing_service.py

Flask microservice (port 5200) that runs the full astrobatch pipeline:
  scan → split → calibrate (Siril) → plate-solve → upload to STDWeb

Called by the Night Manager Node.js server when a pipeline job is triggered.

Endpoints:
  POST /process   { job_id, fits_dir, target }
  GET  /health
  GET  /jobs/<job_id>/log
"""

import logging
import os
import re
import shutil
import sqlite3
import subprocess
import threading
import queue
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import requests
from astropy.io import fits
from flask import Flask, jsonify, request
from flask_cors import CORS

# ── Paths ─────────────────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CALIB_DIR    = PROJECT_ROOT / "calib"
DATA_DIR     = Path(os.environ.get("DATA_DIR", PROJECT_ROOT / "data"))
LOG_DIR      = PROJECT_ROOT / "logs"
NAS_OUTPUT   = Path(os.environ.get("NAS_OUTPUT", "/mnt/nas/input/pyl/astro/output"))
DB_PATH      = Path(os.environ.get(
    "NIGHTMANAGER_DB",
    PROJECT_ROOT / "nightmanager.db"
))

LOG_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)

# ── STDWeb ────────────────────────────────────────────────────────────────────

STDWEB_URL   = os.environ.get("STDWEB_URL",   "http://86.253.141.183:7000")
STDWEB_TOKEN = os.environ.get("STDWEB_TOKEN", "1e296ddd6738af45467b7bc6558c00a9524447ab")
SIRIL_BIN    = os.environ.get("SIRIL_BIN",    "siril")

# IMX533 @ NINA gain=10 → ~2.8 e/ADU (14-bit ADC, scaled ×4 to 16-bit by NINA,
# then Siril normalises to [0,1] by dividing by 65535).
# Effective gain for STDWeb = real_gain × 2^14 = 2.8 × 16384 ≈ 45875
# Override with STDWEB_GAIN env var if camera/gain changes.
STDWEB_GAIN  = int(os.environ.get("STDWEB_GAIN", str(int(2.8 * 16384))))

# ── Filter name → flat lookup key mapping ─────────────────────────────────────
# Maps NINA filter name → key in astrobatch calib/flats/ discovery dict.
# calib/flats/flat_Gbp.fit  → key "GBP"
# calib/flats/flat_G.fit    → key "G"
# calib/flats/flat_Grp.fit  → key "GRP"

FILTER_FLAT_MAP: dict[str, str] = {
    "G":        "G",
    "BP":       "GBP",
    "RP":       "GRP",
    "FILTER 5": "FILTER5",
    "FILTER5":  "FILTER5",
}

# ── Logging ───────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("processing_service")

# Siril (AppImage GTK) needs a display even in pipe mode.
# On Xubuntu the local session is always at :0.
if not os.environ.get("DISPLAY"):
    os.environ["DISPLAY"] = ":0"
    log.info("DISPLAY not set — defaulting to :0 for Siril/GTK")

# ── Flask app ─────────────────────────────────────────────────────────────────

app = Flask(__name__)
CORS(app)

# ── Database helpers ──────────────────────────────────────────────────────────

def _db():
    conn = sqlite3.connect(str(DB_PATH), timeout=10)
    conn.row_factory = sqlite3.Row
    return conn


def _now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def update_job(job_id: int, status: str, error: str | None = None,
               stdweb_task_id: str | None = None, stdweb_url: str | None = None):
    try:
        with _db() as conn:
            fields = ["status = ?", "updated_at = ?"]
            values: list = [status, _now()]
            if error is not None:
                fields.append("error = ?");        values.append(error)
            if stdweb_task_id is not None:
                fields.append("stdweb_task_id = ?"); values.append(stdweb_task_id)
            if stdweb_url is not None:
                fields.append("stdweb_url = ?");    values.append(stdweb_url)
            values.append(job_id)
            conn.execute(
                f"UPDATE pipeline_jobs SET {', '.join(fields)} WHERE id = ?",
                values
            )
    except Exception as exc:
        log.error("DB update_job failed (job=%s status=%s): %s", job_id, status, exc)


def _ensure_phot_columns():
    """Add photometry columns to pipeline_results if they don't exist yet."""
    cols = [
        ("obs_date",   "TEXT"),
        ("mjd",        "REAL"),
        ("mag_ap",     "REAL"),
        ("magerr_ap",  "REAL"),
        ("mag_sub",    "REAL"),
        ("magerr_sub", "REAL"),
        ("mag_sub_ul", "REAL"),
    ]
    with _db() as conn:
        for col, typ in cols:
            try:
                conn.execute(f"ALTER TABLE pipeline_results ADD COLUMN {col} {typ}")
            except Exception:
                pass  # already exists

_ensure_phot_columns()


def create_result(job_id: int, filt: str, exposure: str, n_frames: int,
                  target: str = "", obs_date: str | None = None) -> int:
    """Insert a pipeline_results row for one filter/target. Returns the new row id."""
    try:
        with _db() as conn:
            cur = conn.execute(
                "INSERT INTO pipeline_results (job_id, filter, exposure, n_frames, status, target, obs_date) "
                "VALUES (?, ?, ?, ?, 'processing', ?, ?)",
                (job_id, filt, exposure, n_frames, target, obs_date)
            )
            return cur.lastrowid
    except Exception as exc:
        log.error("DB insert result failed: %s", exc)
        return -1


def update_result(result_id: int, status: str, error: str | None = None,
                  stdweb_task_id: str | None = None, stdweb_url: str | None = None,
                  frame_previews: list | None = None, stdweb_state: str | None = None,
                  obs_date: str | None = None, mjd: float | None = None,
                  mag_ap: float | None = None, magerr_ap: float | None = None,
                  mag_sub: float | None = None, magerr_sub: float | None = None,
                  mag_sub_ul: float | None = None):
    """Update a pipeline_results row."""
    import json as _json
    try:
        with _db() as conn:
            fields = ["status = ?", "updated_at = ?"]
            values: list = [status, _now()]
            if error is not None:
                fields.append("error = ?");          values.append(error)
            if stdweb_task_id is not None:
                fields.append("stdweb_task_id = ?"); values.append(stdweb_task_id)
            if stdweb_url is not None:
                fields.append("stdweb_url = ?");     values.append(stdweb_url)
            if frame_previews is not None:
                fields.append("frame_previews = ?"); values.append(_json.dumps(frame_previews))
            if stdweb_state is not None:
                fields.append("stdweb_state = ?");   values.append(stdweb_state)
            if obs_date is not None:
                fields.append("obs_date = ?");       values.append(obs_date)
            if mjd is not None:
                fields.append("mjd = ?");            values.append(mjd)
            if mag_ap is not None:
                fields.append("mag_ap = ?");         values.append(mag_ap)
            if magerr_ap is not None:
                fields.append("magerr_ap = ?");      values.append(magerr_ap)
            if mag_sub is not None:
                fields.append("mag_sub = ?");        values.append(mag_sub)
            if magerr_sub is not None:
                fields.append("magerr_sub = ?");     values.append(magerr_sub)
            if mag_sub_ul is not None:
                fields.append("mag_sub_ul = ?");     values.append(mag_sub_ul)
            values.append(result_id)
            conn.execute(
                f"UPDATE pipeline_results SET {', '.join(fields)} WHERE id = ?",
                values
            )
    except Exception as exc:
        log.error("DB update_result failed (result=%s status=%s): %s", result_id, status, exc)

# ── Job logger ────────────────────────────────────────────────────────────────

class JobLogger:
    """Writes timestamped log lines to a per-job file."""
    def __init__(self, job_id: int):
        self.path = LOG_DIR / f"job_{job_id}.log"
        self._fh  = open(self.path, "a", buffering=1)

    def info(self, msg: str):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        self._fh.write(line + "\n")
        log.info("[job] %s", msg)

    def error(self, msg: str):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] ✗ {msg}"
        self._fh.write(line + "\n")
        log.error("[job] %s", msg)

    def close(self):
        self._fh.close()

# ── Calibration frame discovery ───────────────────────────────────────────────

def discover_cal_frames() -> tuple[dict, dict]:
    """Return (darks, flats) dicts: exposure_str → path, filter_key → path."""
    darks: dict[str, str] = {}
    for p in (CALIB_DIR / "darks").glob("dark_*s.fit"):
        m = re.match(r"dark_(\d+)s\.fit", p.name, re.IGNORECASE)
        if m:
            darks[m.group(1)] = str(p)

    flats: dict[str, str] = {}
    for p in (CALIB_DIR / "flats").glob("flat_*.fit"):
        m = re.match(r"flat_([A-Za-z0-9]+)\.fit", p.name, re.IGNORECASE)
        if m:
            flats[m.group(1).upper()] = str(p)

    return darks, flats

# ── FITS scanning ─────────────────────────────────────────────────────────────

def _obj_slug(name: str) -> str:
    """Normalise an object name to a slug for comparison: 'SN 2026fvx' → 'sn_2026fvx'."""
    return name.strip().lower().replace(" ", "_")


def scan_fits(fits_dir: Path, jlog: JobLogger,
              selected_files: list[str] | None = None,
              object_filter: str | None = None) -> tuple[dict, dict]:
    """
    Scan fits_dir for FITS files grouped by (object_slug, filter, exp_str).
    If selected_files is provided, those paths are used directly — no rglob.
    If object_filter is provided (e.g. 'SN 2026fvx'), only groups whose
    object slug matches are returned — critical for SNAPSHOT dirs with
    multiple mixed targets.
    Returns: (groups, raw_names)
      groups    : { (slug, filter, exp_str): [Path, ...] }
      raw_names : { slug: original_OBJECT_header }  e.g. 'sn_2026fvx' → 'SN 2026fvx'
    """
    folder_fallback = fits_dir.name
    if folder_fallback.upper() == "SNAPSHOT":
        folder_fallback = "unknown"

    filter_slug = _obj_slug(object_filter) if object_filter else None

    groups: dict[tuple, list] = defaultdict(list)
    raw_names: dict[str, str] = {}       # slug → first raw OBJECT seen

    def _classify(f: Path):
        """Read headers and add f to the right group."""
        try:
            hdr     = fits.getheader(str(f), ext=0)
            filt    = str(hdr.get("FILTER", "Unknown")).strip()
            exp     = float(hdr.get("EXPTIME", hdr.get("EXPOSURE", 0)))
            exp_str = str(int(round(exp)))
            raw_obj = str(hdr.get("OBJECT", "")).strip()
            if not raw_obj or raw_obj.lower() in ("snapshot", "unknown", "none", ""):
                raw_obj = folder_fallback
            obj = _obj_slug(raw_obj)
            # Skip if we're filtering to a specific target and this file doesn't match
            if filter_slug and obj != filter_slug:
                return
            groups[(obj, filt, exp_str)].append(f)
            raw_names.setdefault(obj, raw_obj)   # keep first occurrence
        except Exception as exc:
            jlog.error(f"Cannot read {f.name}: {exc}")

    if selected_files:
        # Use the provided list directly — no rglob, no path matching
        candidates = [Path(p) for p in selected_files if p]
        jlog.info(f"  Using {len(candidates)} selected file(s) (skipping full directory scan)")
        for f in candidates:
            _classify(f)
    else:
        for ext in ("*.fit", "*.fits", "*.FIT", "*.FITS"):
            for f in fits_dir.rglob(ext):
                _classify(f)

    if filter_slug:
        jlog.info(f"  Object filter: '{object_filter}' (slug={filter_slug})")
    for (obj, filt, exp_str), files in groups.items():
        jlog.info(f"  Found {len(files)} frames  object={obj}  filter={filt}  exp={exp_str}s")

    return dict(groups), raw_names

# ── Work directory setup ──────────────────────────────────────────────────────

def prepare_work_dir(target: str, filt: str, exp_str: str, date_str: str,
                     files: list[Path], jlog: JobLogger,
                     clean_first: bool = False) -> Path:
    """Copy raw files into the work tree and return the folder path.

    When clean_first=True (re-process a subset), remove existing science FITS
    from the work dir so only the selected files are present for Siril.
    Uses up to 4 parallel copy threads to saturate NAS bandwidth.
    """
    work_folder = DATA_DIR / date_str / target / filt / f"{exp_str}s"
    work_folder.mkdir(parents=True, exist_ok=True)

    if clean_first:
        for old in work_folder.glob("*.fit*"):
            old.unlink(missing_ok=True)
        jlog.info(f"  Work dir cleaned for re-processing subset")

    existing = {f.name for f in work_folder.glob("*.fit*")}
    to_copy = [src for src in files if src.name not in existing]

    if not to_copy:
        jlog.info(f"  Work dir: {work_folder}  (all {len(files)} files already present)")
        return work_folder

    jlog.info(f"  Work dir: {work_folder}")
    jlog.info(f"  Copying {len(to_copy)} new file(s) of {len(files)} total…")

    done_count = 0
    lock = threading.Lock()

    def _copy_one(src: Path):
        nonlocal done_count
        dst = work_folder / src.name
        shutil.copy2(src, dst)
        with lock:
            done_count += 1
            if done_count % 10 == 0 or done_count == len(to_copy):
                jlog.info(f"    Copied {done_count}/{len(to_copy)} files…")

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(_copy_one, to_copy))

    jlog.info(f"  Copy done — {len(to_copy)} file(s) copied to {work_folder}")
    return work_folder

def finalise_work_dir(work_dir: Path, date_str: str, target: str,
                      filt: str, jlog: JobLogger) -> Path | None:
    """
    After a successful pipeline run:
      1. Move res.fit + res_preview.png to NAS_OUTPUT/<date>/<target>/<filter>/
      2. Delete all Siril intermediary files (i_*, pp_i_*, r_pp_i_*, *.seq, *.ssf)
         but keep the raw input .fits, _calib/, _previews/ subdirs intact.
    Returns the NAS destination folder, or None if the move failed.
    """
    nas_dest = NAS_OUTPUT / date_str / target / filt
    try:
        nas_dest.mkdir(parents=True, exist_ok=True)
    except Exception as exc:
        jlog.error(f"  Cannot create NAS output dir {nas_dest}: {exc}")
        return None

    # Move the two result files
    moved = 0
    for fname in ("res.fit", "res_preview.png"):
        src = work_dir / fname
        if src.exists():
            dst = nas_dest / fname
            shutil.move(str(src), dst)
            moved += 1
            jlog.info(f"  Moved {fname} → {dst}")

    if moved == 0:
        jlog.error("  No result files found to move")
        return None

    # Clean up Siril intermediaries (keep raw .fits and subdirs)
    cleanup_patterns = (
        "i_*.fit", "i_*.fits", "i_conversion.txt",
        "pp_i_*.fit", "pp_i_*.fits",
        "r_pp_i_*.fit", "r_pp_i_*.fits",
        "*.seq", "*.ssf", "*.lst",
    )
    removed = 0
    for pattern in cleanup_patterns:
        for f in work_dir.glob(pattern):
            f.unlink(missing_ok=True)
            removed += 1
    jlog.info(f"  Cleaned {removed} intermediary file(s) from work dir")
    return nas_dest


# ── Siril calibration + stacking ──────────────────────────────────────────────

def _run_siril_script(siril_exe: str, script_path: Path, folder: Path,
                      label: str, jlog: JobLogger) -> tuple[bool, bool]:
    """
    Execute a Siril script, stream output to jlog.
    Returns (success, no_stars) — no_stars=True when registration failed
    specifically because no stars were found in the reference frame.
    """
    jlog.info(f"  Siril script ({label}): {script_path}")
    jlog.info("  --- Siril output ---")
    no_stars = False
    try:
        env = os.environ.copy()
        env.setdefault("DISPLAY", ":0")
        proc = subprocess.Popen(
            [siril_exe, "-s", str(script_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
        )
        for line in proc.stdout:
            line = line.rstrip()
            if line:
                jlog.info(f"  [siril] {line}")
                if "not enough stars" in line.lower() or "found 0 stars" in line.lower():
                    no_stars = True
        proc.wait()
        jlog.info("  --- Siril end ---")
        if proc.returncode != 0:
            jlog.error(f"  Siril exited with code {proc.returncode}")
            return False, no_stars
    except Exception as exc:
        jlog.error(f"  Siril launch error: {exc}")
        return False, no_stars
    return True, False


def generate_frame_previews(folder: Path, jlog: JobLogger,
                            source_files: list[Path] | None = None,
                            max_frames: int = 6) -> list[dict]:
    """
    Generate PNG thumbnails for the first `max_frames` calibrated frames (pp_i_*.fit).
    Saved in _previews/ subfolder so Siril's `convert i` never picks them up.

    Returns a list of dicts: [{"preview": "<rel_data_url>", "source": "<abs_path>"}]
    `source` is the original NAS FITS path so the user can re-trigger with a subset.
    """
    import numpy as _np
    from PIL import Image as _Image

    previews_dir = folder / "_previews"
    previews_dir.mkdir(exist_ok=True)

    # pp_i_N maps positionally to source_files[N-1] (sorted order preserved by Siril)
    pp_files = sorted(folder.glob("pp_i_*.fit"))[:max_frames]
    if not pp_files:
        jlog.info("  No calibrated frames found for preview generation")
        return []

    # Siril numbers frames from the sorted list of original files it converted
    # source_files is that same sorted list; index them in order.
    src_list: list[Path | None] = (
        sorted(source_files)[:max_frames] if source_files else [None] * len(pp_files)
    )

    jlog.info(f"  Generating {len(pp_files)} frame preview(s)...")
    results: list[dict] = []
    for idx, fit_path in enumerate(pp_files, start=1):
        src = src_list[idx - 1] if idx - 1 < len(src_list) else None
        out = previews_dir / f"frame_{idx:03d}_preview.png"
        try:
            with fits.open(str(fit_path)) as hdul:
                data = hdul[0].data.astype(_np.float32)
            lo, hi = _np.percentile(data, (0.5, 99.5))
            if hi > lo:
                data = _np.clip((data - lo) / (hi - lo), 0, 1)
            else:
                data = _np.zeros_like(data)
            img = _Image.fromarray((data * 255).astype(_np.uint8))
            w, h = img.size
            if max(w, h) > 800:
                scale = 800 / max(w, h)
                img = img.resize((int(w * scale), int(h * scale)), _Image.LANCZOS)
            img.save(str(out), format="PNG", optimize=True)
            try:
                rel = str(out.relative_to(DATA_DIR))
            except ValueError:
                rel = str(out)
            entry: dict = {"preview": rel, "source": str(src) if src else ""}
            results.append(entry)
            jlog.info(f"  Frame {idx} preview → {out.name}  source={src.name if src else '?'}")
        except Exception as exc:
            jlog.info(f"  Frame {idx} preview failed (non-fatal): {exc}")

    return results


def run_siril(folder: Path, flat_path: str | None, dark_path: str | None,
              n_images: int, jlog: JobLogger,
              source_files: list[Path] | None = None) -> tuple[bool, list[dict]]:
    """
    Calibrate + stack with Siril.
    Returns (success, frame_previews).
    frame_previews is non-empty only when registration fails due to no stars —
    contains paths to individual calibrated frame previews for the user to diagnose.
    """
    jlog.info(f"  Running Siril ({n_images} frame{'s' if n_images != 1 else ''})...")

    siril_exe = shutil.which("siril-cli") or shutil.which(SIRIL_BIN) or shutil.which("siril")
    if not siril_exe:
        jlog.error("Siril executable not found in PATH")
        return False, []

    # Stage calibration masters in _calib/ so `convert i` never picks them up
    # as science frames (Siril converts ALL fits in the work dir).
    calib_dir = folder / "_calib"
    calib_dir.mkdir(exist_ok=True)

    def _stage_calib(src: str | None, stem: str) -> str | None:
        if not src:
            return None
        dst = calib_dir / (stem + Path(src).suffix)
        if not dst.exists():
            shutil.copy2(src, dst)
        return f"_calib/{stem}"

    flat_arg = _stage_calib(flat_path, "master_flat")
    dark_arg = _stage_calib(dark_path, "master_dark")

    def _clean_stale():
        # Clean Siril working files + any stale previews left in the root.
        # Subdirectories _calib/ and _previews/ are intentionally left intact.
        for pattern in ("i_*.fit", "i_*.fits", "pp_i_*.fit", "pp_i_*.fits",
                        "r_pp_i_*.fit", "r_pp_i_*.fits", "*.seq", "*.lst",
                        "res.fit", "res.fits", "res_preview.png", "_siril_script.ssf", "_siril_retry.ssf",
                        "frame_*_preview.png"):
            for stale in folder.glob(pattern):
                stale.unlink(missing_ok=True)

    _clean_stale()
    script_path = folder / "_siril_script.ssf"

    # ── Single-frame shortcut ─────────────────────────────────────────────────
    if n_images == 1:
        cal_cmd = "calibrate_single i_00001.fit"
        if flat_arg:
            cal_cmd += f" -flat={flat_arg}"
        if dark_arg:
            cal_cmd += f" -dark={dark_arg}"
        lines = ["requires 1.2.0", f'cd "{folder}"', "convert i -out=.",
                 cal_cmd, "load pp_i_00001", "save res"]
        script_path.write_text("\n".join(lines) + "\n")
        ok, _ = _run_siril_script(siril_exe, script_path, folder, "single-frame", jlog)
        if not ok:
            return False, []

    else:
        # ── Full pipeline: calibrate → register → stack ───────────────────────
        cal_cmd = "calibrate i"
        if flat_arg:
            cal_cmd += f" -flat={flat_arg}"
        if dark_arg:
            cal_cmd += f" -dark={dark_arg}"
        lines = [
            "requires 1.2.0", f'cd "{folder}"', "convert i -out=.",
            cal_cmd,
            "register pp_i -interp=cu",
            "stack r_pp_i rej 3 3 -norm=addscale -filter-fwhm=80% -filter-round=80%",
            "load r_pp_i_stacked",
            "save res",
        ]
        script_path.write_text("\n".join(lines) + "\n")
        ok, no_stars = _run_siril_script(siril_exe, script_path, folder, "calibrate+register+stack", jlog)

        if not ok:
            if no_stars:
                # Auto-drop the bad reference frame and retry until success or
                # too few frames remain.
                register_lines = [
                    "requires 1.2.0", f'cd "{folder}"',
                    "register pp_i -interp=cu",
                    "stack r_pp_i rej 3 3 -norm=addscale -filter-fwhm=80% -filter-round=80%",
                    "load r_pp_i_stacked",
                    "save res",
                ]
                retry_script = folder / "_siril_retry.ssf"
                while no_stars:
                    pp_frames = sorted(folder.glob("pp_i_*.fit"))
                    if len(pp_frames) < 3:
                        jlog.error(f"  Only {len(pp_frames)} frames left — giving up")
                        break
                    bad = pp_frames[0]
                    jlog.info(f"  Dropping frame with no stars: {bad.name}")
                    bad.unlink(missing_ok=True)
                    # Remove stale sequence files so Siril rescans
                    for seq in folder.glob("*.seq"):
                        seq.unlink(missing_ok=True)
                    retry_script.write_text("\n".join(register_lines) + "\n")
                    ok, no_stars = _run_siril_script(
                        siril_exe, retry_script, folder, "register+stack (retry)", jlog
                    )
                retry_script.unlink(missing_ok=True)
                if not ok:
                    jlog.error("  Registration still failed — generating frame previews for diagnosis")
                    previews = generate_frame_previews(folder, jlog, source_files=source_files)
                    return False, previews
            else:
                return False, []

    res = folder / "res.fit"
    if not res.exists():
        jlog.error("  Siril finished but res.fit not found")
        return False, []

    jlog.info(f"  res.fit created ({res.stat().st_size / 1_048_576:.1f} MB)")
    return True, []


def make_preview(fits_path: Path, jlog: JobLogger) -> Path | None:
    """
    Convert a FITS file to a PNG preview with zscale stretch.
    Returns the path to res_preview.png, or None on failure.
    """
    out = fits_path.parent / "res_preview.png"
    try:
        import numpy as np
        from astropy.io import fits as afits
        from astropy.visualization import ZScaleInterval
        from PIL import Image

        with afits.open(str(fits_path)) as hdul:
            data = hdul[0].data
            if data is None and len(hdul) > 1:
                data = hdul[1].data
            data = np.array(data, dtype=np.float32)

        # Collapse colour axis if present (e.g. 3×H×W → H×W via luminance)
        if data.ndim == 3:
            data = data.mean(axis=0)

        # ZScale stretch → [0, 1]
        interval = ZScaleInterval(contrast=0.25)
        vmin, vmax = interval.get_limits(data)
        data = np.clip((data - vmin) / max(vmax - vmin, 1e-10), 0, 1)

        # Flip vertically (FITS origin is bottom-left)
        data = np.flipud(data)

        # Scale to uint8 and save
        img = Image.fromarray((data * 255).astype(np.uint8), mode="L")

        # Resize for web: cap longest side at 1200 px
        max_side = 1200
        w, h = img.size
        if max(w, h) > max_side:
            scale = max_side / max(w, h)
            img = img.resize((int(w * scale), int(h * scale)), Image.LANCZOS)

        img.save(str(out), format="PNG", optimize=True)
        jlog.info(f"  Preview saved: {out.name} ({out.stat().st_size // 1024} KB)")
        return out
    except Exception as exc:
        jlog.info(f"  Preview generation failed (non-fatal): {exc}")
        return None


# ── Plate solving ─────────────────────────────────────────────────────────────

def plate_solve(res_fit: Path, jlog: JobLogger) -> bool:
    """Run astrometry.net solve-field on res.fit, write WCS into the file.
    Returns True if a valid WCS was obtained and written."""
    jlog.info("  Plate solving...")
    try:
        import sys
        if str(PROJECT_ROOT) not in sys.path:
            sys.path.insert(0, str(PROJECT_ROOT))
        from astrobatch.spliter import PlateSolver
        from astropy.io import fits as _fits

        solver = PlateSolver()
        result = solver.solve_field(str(res_fit))

        # solve_field returns (WCS, source_data) or (None, None) on failure
        solved_wcs = result[0] if isinstance(result, tuple) else result

        if solved_wcs is None:
            jlog.error("  Plate solve failed — no WCS solution")
            return False

        # Write the WCS header keywords back into res.fit in-place
        wcs_header = solved_wcs.to_header(relax=True)
        with _fits.open(str(res_fit), mode="update") as hdul:
            hdr = hdul[0].header
            # Remove old broken WCS keys before writing the new solution
            for key in ("CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
                        "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                        "CDELT1", "CDELT2", "CTYPE1", "CTYPE2",
                        "CUNIT1", "CUNIT2", "PC1_1", "PC1_2", "PC2_1", "PC2_2"):
                hdr.remove(key, ignore_missing=True, remove_all=True)
            hdr.update(wcs_header)
            hdul.flush()

        jlog.info(f"  WCS written to res.fit (CRVAL1={wcs_header.get('CRVAL1', '?'):.4f}, "
                  f"CRVAL2={wcs_header.get('CRVAL2', '?'):.4f})")
        return True
    except Exception as exc:
        jlog.error(f"  Plate solve error: {exc}")
        return False

# ── STDWeb upload + pipeline ───────────────────────────────────────────────────

def _read_fits_radec(fits_path: Path) -> tuple[float, float] | tuple[None, None]:
    """Read telescope pointing RA/DEC from FITS header. Returns (ra_deg, dec_deg) or (None, None)."""
    try:
        from astropy.io import fits as _fits
        h = _fits.getheader(str(fits_path))
        ra  = h.get("RA")  or h.get("CRVAL1") or h.get("OBJCTRA")
        dec = h.get("DEC") or h.get("CRVAL2") or h.get("OBJCTDEC")
        if ra is None or dec is None:
            return None, None
        # Convert sexagesimal strings to degrees if needed
        if isinstance(ra, str):
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            c = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
            return float(c.ra.deg), float(c.dec.deg)
        return float(ra), float(dec)
    except Exception:
        return None, None


def stdweb_upload(res_fit: Path, title: str, jlog: JobLogger,
                  target: str | None = None) -> str | None:
    """Upload res.fit to STDWeb. Returns task_id or None."""
    jlog.info(f"  Uploading to STDWeb: {title}" + (f" (target={target})" if target else ""))
    try:
        cmd = [
            "curl", "-s",
            "-H", f"Authorization: Token {STDWEB_TOKEN}",
            "-F", f"file=@{res_fit}",
            "-F", f"title={title}",
            "-F", f"gain={STDWEB_GAIN}",
            "-F", "do_inspect=true",
            "-F", "do_photometry=true",
        ]
        if target:
            cmd += ["-F", f"target={target}"]
        # Pass telescope pointing as blind-match hint so STDWeb doesn't need
        # to resolve the target name via external APIs (Fink/TNS SSL issues).
        ra, dec = _read_fits_radec(res_fit)
        if ra is not None:
            center_str = f"{ra:.6f} {dec:+.6f}"
            cmd += [
                "-F", f"blind_match_center={center_str}",
                "-F", "blind_match_sr0=1.0",
            ]
            jlog.info(f"  blind_match_center={center_str}")
        cmd.append(f"{STDWEB_URL}/api/tasks/upload/")
        # -w outputs HTTP status code after the body
        cmd_with_status = cmd[:-1] + ["-w", "\nHTTP_STATUS:%{http_code}", cmd[-1]]
        result = subprocess.run(cmd_with_status, capture_output=True, text=True, timeout=120)

        raw = result.stdout
        jlog.info(f"  STDWeb raw response: {raw[:400]}")

        import json
        # Strip the HTTP_STATUS trailer before parsing JSON
        body = raw.split("\nHTTP_STATUS:")[0].strip()
        data = json.loads(body)
        # Response: {"message":"...", "task": {"id": 3446, ...}}
        task_obj = data.get("task") or data
        task_id = str(task_obj.get("id") or task_obj.get("task_id") or "")
        if not task_id:
            jlog.error(f"Upload response has no task_id: {raw[:400]}")
            return None
        jlog.info(f"  Uploaded → task_id={task_id}")
        return task_id
    except Exception as exc:
        jlog.error(f"Upload error: {exc} | raw={result.stdout[:400] if 'result' in dir() else 'n/a'}")
        return None


def stdweb_action(task_id: str, action: str, jlog: JobLogger,
                  extra: dict | None = None) -> bool:
    """POST an action to STDWeb task."""
    payload = {"action": action}
    if action == "inspect":
        payload.update({"gain": STDWEB_GAIN, "bias": 0})
    elif action == "photometry":
        payload.update({"sn": 10, "catalog": "gaiadr2"})
    elif action == "subtract":
        payload.update({"template": os.environ.get("STDWEB_TEMPLATE", "ztf")})
    if extra:
        payload.update(extra)
    try:
        r = requests.post(
            f"{STDWEB_URL}/api/tasks/{task_id}/action/",
            headers={"Authorization": f"Token {STDWEB_TOKEN}"},
            json=payload,
            timeout=30,
        )
        jlog.info(f"  {action} → HTTP {r.status_code}")
        return r.status_code < 400
    except Exception as exc:
        jlog.error(f"STDWeb action {action} error: {exc}")
        return False


def wait_state(task_id: str, target_states: list[str],
               jlog: JobLogger, timeout: int = 300) -> str:
    """Poll STDWeb task until it reaches one of target_states. Returns final state."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            r = requests.get(
                f"{STDWEB_URL}/api/tasks/{task_id}/",
                headers={"Authorization": f"Token {STDWEB_TOKEN}"},
                timeout=15,
            )
            state = r.json().get("state", "unknown")
            if state in target_states:
                return state
        except Exception:
            pass
        time.sleep(10)
    return "timeout"

# ── Main pipeline ──────────────────────────────────────────────────────────────

def _detect_date(fits_path: Path) -> str:
    """Walk up the path looking for a YYYY-MM-DD folder name."""
    for part in reversed(fits_path.parts):
        if re.match(r"\d{4}-\d{2}-\d{2}", part):
            return part
    return datetime.now().strftime("%Y-%m-%d")



def run_pipeline(job_id: int, fits_dir: str, target: str,
                 selected_files: list[str] | None = None,
                 object_filter: str | None = None):
    """Full pipeline — runs in a background thread."""
    jlog = JobLogger(job_id)
    try:
        _run_pipeline(job_id, fits_dir, target, jlog,
                      selected_files=selected_files, object_filter=object_filter)
    except Exception as exc:
        log.exception("Unhandled exception in job %s", job_id)
        jlog.error(f"UNHANDLED EXCEPTION: {exc}")
        update_job(job_id, "error", error=str(exc))
    finally:
        jlog.close()


# Steps in order; used by resume to skip already-done work
PIPELINE_STEP_ORDER = ["calibrate", "solve", "upload", "inspect", "photometry", "subtraction"]


def resume_pipeline(result_id: int, from_step: str):
    """Resume a single pipeline result from a given step, in a background thread."""
    with sqlite3.connect(str(DB_PATH), timeout=10) as con:
        con.row_factory = sqlite3.Row
        cur = con.execute("SELECT * FROM pipeline_results WHERE id=?", (result_id,))
        row = cur.fetchone()
        result = dict(row) if row else None

    if not result:
        log.error("resume_pipeline: result %s not found", result_id)
        return

    job_id = result["job_id"]
    with sqlite3.connect(str(DB_PATH), timeout=10) as con:
        con.row_factory = sqlite3.Row
        cur = con.execute("SELECT * FROM pipeline_jobs WHERE id=?", (job_id,))
        row = cur.fetchone()
        job = dict(row) if row else None

    if not job:
        log.error("resume_pipeline: job %s not found", job_id)
        return

    jlog = JobLogger(job_id)
    try:
        _resume_pipeline(result, job, from_step, jlog)
    except Exception as exc:
        log.exception("Unhandled exception resuming result %s from %s", result_id, from_step)
        jlog.error(f"UNHANDLED EXCEPTION: {exc}")
        update_result(result_id, "error", error=str(exc))
    finally:
        jlog.close()


def _resume_pipeline(result: dict, job: dict, from_step: str, jlog: JobLogger):
    result_id = result["id"]
    job_id    = result["job_id"]
    filt      = result["filter"] or "?"
    exp_str   = result["exposure"] or "?"
    obj       = result["target"] or job["target"] or "unknown"
    date_str  = result.get("obs_date") or _detect_date(Path(job["fits_dir"]))

    jlog.info(f"=== Resuming result {result_id} from step '{from_step}' ===")
    jlog.info(f"  object={obj}  filter={filt}  exp={exp_str}s")

    # Clear previous error explicitly and reset job to running
    try:
        with _db() as conn:
            conn.execute("UPDATE pipeline_results SET status='processing', error=NULL, updated_at=? WHERE id=?",
                         (_now(), result_id))
            conn.execute("UPDATE pipeline_jobs SET status='uploading', error=NULL, updated_at=? WHERE id=?",
                         (_now(), job_id))
    except Exception as exc:
        log.error("resume clear error failed: %s", exc)

    step_idx = PIPELINE_STEP_ORDER.index(from_step) if from_step in PIPELINE_STEP_ORDER else 0

    # ── calibrate ─────────────────────────────────────────────────────────────
    if step_idx <= PIPELINE_STEP_ORDER.index("calibrate"):
        fits_path = Path(job["fits_dir"])
        groups, raw_names = scan_fits(fits_path, jlog, object_filter=obj)
        key = next((k for k in groups if k[0] == obj and k[1] == filt), None)
        if not key:
            update_result(result_id, "error", error="No matching FITS files found for re-calibration")
            update_job(job_id, "error", error="No matching FITS files found for re-calibration")
            return
        files = groups[key]

        darks, flats = discover_cal_frames()
        work_dir = prepare_work_dir(obj, filt, exp_str, date_str, files, jlog)

        flat_key  = FILTER_FLAT_MAP.get(filt.upper(), filt.upper())
        flat_path = flats.get(flat_key)
        dark_path = darks.get(exp_str)

        update_result(result_id, "calibrating")
        update_job(job_id, "calibrating")
        ok, frame_previews = run_siril(work_dir, flat_path, dark_path, len(files), jlog,
                                       source_files=files)
        if not ok:
            err_msg = "No stars found — see frame previews" if frame_previews else "Siril calibration failed"
            update_result(result_id, "error", error=err_msg,
                          frame_previews=frame_previews if frame_previews else None)
            update_job(job_id, "error", error=err_msg)
            return

        res_fit = work_dir / "res.fit"
        make_preview(res_fit, jlog)
        from_step = "solve"
        step_idx  = PIPELINE_STEP_ORDER.index("solve")
    else:
        # Work dir must already exist
        work_dir = DATA_DIR / date_str / obj / filt / f"{exp_str}s"
        res_fit  = work_dir / "res.fit"
        if not res_fit.exists():
            update_result(result_id, "error", error=f"Stacked file missing: {res_fit}")
            update_job(job_id, "error", error=f"Stacked file missing: {res_fit}")
            return

    # ── solve ─────────────────────────────────────────────────────────────────
    if step_idx <= PIPELINE_STEP_ORDER.index("solve"):
        update_result(result_id, "solving")
        update_job(job_id, "solving")
        jlog.info("Plate solving…")
        solve_ok = plate_solve(res_fit, jlog)
        if not solve_ok:
            jlog.info("  Plate solve failed — continuing without WCS")
        from_step = "upload"
        step_idx  = PIPELINE_STEP_ORDER.index("upload")

    # ── upload ────────────────────────────────────────────────────────────────
    if step_idx <= PIPELINE_STEP_ORDER.index("upload"):
        update_result(result_id, "uploading")
        update_job(job_id, "uploading")
        jlog.info("Uploading to STDWeb…")
        fits_path  = Path(job["fits_dir"])
        _, raw_names = scan_fits(fits_path, jlog, object_filter=obj)
        raw_obj     = raw_names.get(obj, obj.replace('_', ' ').title())
        target_clean = raw_obj.replace(' ', '')
        title        = f"{raw_obj} {filt} {exp_str}s {date_str}"
        jlog.info(f"  STDWeb target: '{target_clean}'")
        task_id = stdweb_upload(res_fit, title, jlog, target=target_clean)
        if not task_id:
            update_result(result_id, "error", error="STDWeb upload failed")
            update_job(job_id, "error", error="STDWeb upload failed")
            return
        task_url = f"{STDWEB_URL}/tasks/{task_id}"
        update_result(result_id, "uploaded", stdweb_task_id=task_id, stdweb_url=task_url)
        from_step = "inspect"
        step_idx  = PIPELINE_STEP_ORDER.index("inspect")
    else:
        task_id  = result.get("stdweb_task_id")
        task_url = result.get("stdweb_url") or (f"{STDWEB_URL}/tasks/{task_id}" if task_id else None)
        if not task_id:
            update_result(result_id, "error", error="No STDWeb task ID — re-run from upload")
            update_job(job_id, "error", error="No STDWeb task ID — re-run from upload")
            return

    # ── inspect ───────────────────────────────────────────────────────────────
    if step_idx <= PIPELINE_STEP_ORDER.index("inspect"):
        update_result(result_id, "inspecting")
        jlog.info("Triggering inspection…")
        stdweb_action(task_id, "inspect", jlog)
        state = wait_state(task_id, ["inspect_done", "inspect_failed", "failed", "error"], jlog, timeout=300)
        jlog.info(f"  Inspection state: {state}")
        if state in ("inspect_failed", "failed", "error", "timeout"):
            update_result(result_id, "error", error=f"Inspection {state}")
            update_job(job_id, "error", error=f"Inspection {state}")
            return
        from_step = "photometry"
        step_idx  = PIPELINE_STEP_ORDER.index("photometry")

    # ── photometry ────────────────────────────────────────────────────────────
    if step_idx <= PIPELINE_STEP_ORDER.index("photometry"):
        update_result(result_id, "photometry")
        jlog.info("Triggering photometry…")
        stdweb_action(task_id, "photometry", jlog)
        wait_state(task_id, ["photometry_done", "done", "completed", "failed", "error"], jlog, timeout=600)
        from_step = "subtraction"
        step_idx  = PIPELINE_STEP_ORDER.index("subtraction")

    # ── subtraction ───────────────────────────────────────────────────────────
    if step_idx <= PIPELINE_STEP_ORDER.index("subtraction"):
        update_result(result_id, "subtraction")
        jlog.info("Triggering template subtraction…")
        stdweb_action(task_id, "subtraction", jlog)
        wait_state(task_id, ["subtraction_done", "done", "completed", "failed", "error"], jlog, timeout=600)

    update_result(result_id, "done", stdweb_task_id=task_id, stdweb_url=task_url)
    update_job(job_id, "done")
    jlog.info(f"=== Resume of result {result_id} complete ===")


def _run_pipeline(job_id: int, fits_dir: str, target: str, jlog: JobLogger,
                  selected_files: list[str] | None = None,
                  object_filter: str | None = None):
    update_job(job_id, "running")
    jlog.info(f"=== Job {job_id} started ===")
    jlog.info(f"Input dir : {fits_dir}")
    if object_filter:
        jlog.info(f"Object    : {object_filter} (filtering to this target only)")
    if selected_files:
        jlog.info(f"Selected  : {len(selected_files)} specific file(s)")

    # Purge any stale result rows from a previous run of this job so we start
    # with a clean slate — avoids ghost errors in the UI from old failures.
    try:
        with _db() as conn:
            deleted = conn.execute(
                "DELETE FROM pipeline_results WHERE job_id=? AND status != 'done'",
                (job_id,)
            ).rowcount
        if deleted:
            jlog.info(f"Purged {deleted} stale result row(s) from previous run")
    except Exception as exc:
        jlog.warning(f"Could not purge old results: {exc}")

    fits_path = Path(fits_dir)
    date_str  = _detect_date(fits_path)

    darks, flats = discover_cal_frames()
    jlog.info(f"Darks: {list(darks.keys())}  Flats: {list(flats.keys())}")

    # ── Step 1: Scan ──────────────────────────────────────────────────────────
    update_job(job_id, "scanning")
    jlog.info("Step 1: Scanning FITS files...")
    groups, raw_names = scan_fits(fits_path, jlog, selected_files=selected_files,
                                  object_filter=object_filter)
    if not groups:
        update_job(job_id, "error", error="No valid FITS files found")
        jlog.error("No FITS files found")
        jlog.close()
        return

    success_count = 0
    total_groups = len(groups)
    objects = sorted({obj for (obj, *_) in groups})
    jlog.info(f"Objects found: {objects}  —  {total_groups} group(s)")

    for (obj, filt, exp_str), files in groups.items():
        jlog.info(f"\n--- Processing object={obj}  filter={filt}  exp={exp_str}s ({len(files)} frames) ---")

        # Create a per-result row immediately
        result_id = create_result(job_id, filt, exp_str, len(files), target=obj, obs_date=date_str)

        # ── Step 2: Copy to work dir ──────────────────────────────────────────
        update_job(job_id, "splitting")
        work_dir = prepare_work_dir(obj, filt, exp_str, date_str, files, jlog,
                                    clean_first=selected_files is not None)

        # ── Step 3: Calibrate + Stack ─────────────────────────────────────────
        update_job(job_id, "calibrating")
        jlog.info("Step 3: Calibrating + stacking...")

        flat_key = FILTER_FLAT_MAP.get(filt.upper(), filt.upper())
        flat_path = flats.get(flat_key)
        dark_path = darks.get(exp_str)

        if flat_path:
            jlog.info(f"  Using flat : {Path(flat_path).name}")
        else:
            jlog.info(f"  ⚠️  No flat for filter '{filt}' (key={flat_key}) — calibrating without")

        # Validate dark compatibility: compare GAIN with a science frame
        if dark_path and files:
            try:
                sci_hdr  = fits.getheader(str(files[0]))
                dark_hdr = fits.getheader(dark_path)
                sci_gain  = sci_hdr.get("GAIN")
                dark_gain = dark_hdr.get("GAIN")
                if sci_gain is not None and dark_gain is not None and abs(sci_gain - dark_gain) > max(5, sci_gain * 0.2):
                    jlog.info(f"  ⚠️  Dark GAIN={dark_gain} ≠ science GAIN={sci_gain} — skipping dark (incompatible)")
                    dark_path = None
                else:
                    jlog.info(f"  Using dark : {Path(dark_path).name}  (GAIN={dark_gain})")
            except Exception as _e:
                jlog.info(f"  Using dark : {Path(dark_path).name}")
        elif dark_path:
            jlog.info(f"  Using dark : {Path(dark_path).name}")
        else:
            jlog.info(f"  ⚠️  No dark for {exp_str}s — calibrating without")

        ok, frame_previews = run_siril(work_dir, flat_path, dark_path, len(files), jlog,
                                       source_files=files)
        if not ok:
            err_msg = "No stars found — see frame previews" if frame_previews else "Siril calibration failed"
            jlog.error(f"Calibration failed for {obj}/{filt}/{exp_str}s — skipping")
            update_result(result_id, "error", error=err_msg,
                          frame_previews=frame_previews if frame_previews else None)
            continue

        res_fit = work_dir / "res.fit"
        make_preview(res_fit, jlog)

        # ── Step 4: Plate solve ───────────────────────────────────────────────
        update_job(job_id, "solving")
        update_result(result_id, "solving")
        jlog.info("Step 4: Plate solving...")
        solve_ok = plate_solve(res_fit, jlog)
        if not solve_ok:
            jlog.info("  Plate solve failed — continuing without WCS")

        # ── Step 5: Upload to STDWeb ──────────────────────────────────────────
        update_job(job_id, "uploading")
        update_result(result_id, "uploading")
        jlog.info("Step 5: Uploading to STDWeb...")
        # Use original OBJECT header so STDWeb gets the proper TNS name
        # e.g. raw "SN 2026fvx" → target "SN2026fvx",  "AT 2026fuh" → "AT2026fuh"
        raw_obj = raw_names.get(obj, obj.replace('_', ' ').title())
        title = f"{raw_obj} {filt} {exp_str}s {date_str}"
        target_clean = raw_obj.replace(' ', '')   # "SN 2026fvx" → "SN2026fvx"
        jlog.info(f"  STDWeb target: '{target_clean}'")
        task_id = stdweb_upload(res_fit, title, jlog, target=target_clean)

        if not task_id:
            jlog.error("Upload failed")
            update_result(result_id, "error", error="STDWeb upload failed")
            continue

        task_url = f"{STDWEB_URL}/tasks/{task_id}"
        update_result(result_id, "uploaded", stdweb_task_id=task_id, stdweb_url=task_url)
        # Keep the job's stdweb fields pointing to the most recent upload
        update_job(job_id, "uploading", stdweb_task_id=task_id, stdweb_url=task_url)

        # ── Step 6: Inspect ───────────────────────────────────────────────────
        jlog.info("Step 6: Triggering inspection...")
        update_result(result_id, "inspecting")
        stdweb_action(task_id, "inspect", jlog)
        state = wait_state(task_id,
                           ["inspect_done", "inspect_failed", "failed", "error"],
                           jlog, timeout=300)
        jlog.info(f"  Inspection state: {state}")
        if state in ("inspect_failed", "failed", "error", "timeout"):
            jlog.error(f"  Inspection did not complete (state={state}) — aborting STDWeb steps")
            update_result(result_id, "error", error=f"Inspection {state}")
            continue

        # ── Step 7: Photometry ────────────────────────────────────────────────
        jlog.info("Step 7: Triggering photometry...")
        update_result(result_id, "photometry")
        stdweb_action(task_id, "photometry", jlog)
        state = wait_state(task_id,
                           ["photometry_done", "done", "completed", "failed", "error"],
                           jlog, timeout=600)
        jlog.info(f"  Photometry state: {state}")

        # ── Step 8: Template subtraction ──────────────────────────────────────
        jlog.info("Step 8: Triggering template subtraction...")
        update_result(result_id, "subtraction")
        stdweb_action(task_id, "subtraction", jlog)
        state = wait_state(task_id,
                           ["subtraction_done", "done", "completed", "failed", "error"],
                           jlog, timeout=600)
        jlog.info(f"  Subtraction state: {state}")

        update_result(result_id, "done", stdweb_task_id=task_id, stdweb_url=task_url)

        # ── Step 9: Move results to NAS, clean intermediaries ─────────────────
        jlog.info("Step 9: Archiving results to NAS...")
        finalise_work_dir(work_dir, date_str, obj, filt, jlog)

        success_count += 1

    # ── Done ──────────────────────────────────────────────────────────────────
    if success_count == total_groups and total_groups > 0:
        update_job(job_id, "done")
        jlog.info(f"\n=== Job {job_id} complete: {success_count}/{total_groups} group(s) processed ===")
    elif success_count > 0:
        update_job(job_id, "done", error=f"{total_groups - success_count} group(s) failed")
        jlog.info(f"\n=== Job {job_id} partial: {success_count}/{total_groups} succeeded ===")
    else:
        update_job(job_id, "error", error="All groups failed")
        jlog.error(f"\n=== Job {job_id} failed ===")

# ── Sequential job queue ──────────────────────────────────────────────────────
# A single worker thread drains this queue so jobs run one at a time.
# Each item is a (job_id, callable) tuple so we can deduplicate and track.
_job_queue: queue.Queue = queue.Queue()
_worker_lock  = threading.Lock()
_active_jobs: set[int] = set()      # job_ids currently queued OR running
_active_lock  = threading.Lock()
_current_job_id: int | None = None  # job_id currently executing


def _enqueue_job(job_id: int, task_fn) -> bool:
    """Add job_id to queue. Returns False (409) if already queued/running."""
    with _active_lock:
        if job_id in _active_jobs:
            return False
        _active_jobs.add(job_id)
    _job_queue.put((job_id, task_fn))
    return True


def _queue_worker():
    global _current_job_id
    log.info("Pipeline worker thread started")
    while True:
        try:
            item = _job_queue.get(timeout=5)
        except queue.Empty:
            continue
        job_id, task = item
        _current_job_id = job_id
        try:
            task()
        except BaseException as exc:
            log.exception("Unhandled exception in queue worker task: %s", exc)
        finally:
            _current_job_id = None
            with _active_lock:
                _active_jobs.discard(job_id)
            try:
                _job_queue.task_done()
            except Exception:
                pass
    log.warning("Pipeline worker thread exiting — watchdog will restart")


def _ensure_worker():
    """Start or restart the worker thread if it is not alive."""
    global _worker_thread
    with _worker_lock:
        if not _worker_thread.is_alive():
            log.warning("Pipeline worker thread is dead — restarting")
            _worker_thread = threading.Thread(
                target=_queue_worker, daemon=True, name="pipeline-worker"
            )
            _worker_thread.start()


def _watchdog():
    """Periodically check the worker thread is alive and restart if needed."""
    while True:
        time.sleep(10)
        _ensure_worker()


_worker_thread = threading.Thread(target=_queue_worker, daemon=True, name="pipeline-worker")
_worker_thread.start()
threading.Thread(target=_watchdog, daemon=True, name="pipeline-watchdog").start()


def _recover_queued_jobs():
    """On startup, re-enqueue any jobs left in an interrupted state.

    Includes 'running' (crashed mid-execution) and all intermediate step statuses.
    Deletes non-done pipeline_results before re-queuing to avoid duplicates.
    """
    RECOVERABLE = {
        'queued', 'running',
        'scanning', 'splitting', 'calibrating', 'solving',
        'uploading', 'inspecting', 'photometry', 'subtraction',
    }
    try:
        with _db() as conn:
            rows = conn.execute(
                "SELECT id, fits_dir, target, target_filter, status FROM pipeline_jobs "
                "WHERE status IN ({}) ORDER BY id ASC".format(
                    ",".join("?" * len(RECOVERABLE))),
                tuple(RECOVERABLE)
            ).fetchall()
        if not rows:
            return
        log.info("Recovering %d interrupted job(s) from DB", len(rows))
        for row in rows:
            job_id, fits_dir, target, target_filter, old_status = row
            if not fits_dir or not Path(fits_dir).exists():
                log.warning("Skipping job %s — fits_dir missing: %s", job_id, fits_dir)
                with _db() as conn:
                    conn.execute(
                        "UPDATE pipeline_jobs SET status='error', error='Directory not found on recovery' WHERE id=?",
                        (job_id,)
                    )
                continue
            # If every result is already done, just mark the job done — no re-run
            with _db() as conn:
                total = conn.execute(
                    "SELECT COUNT(*) FROM pipeline_results WHERE job_id=?", (job_id,)
                ).fetchone()[0]
                done  = conn.execute(
                    "SELECT COUNT(*) FROM pipeline_results WHERE job_id=? AND status='done'",
                    (job_id,)
                ).fetchone()[0]
            if total > 0 and done == total:
                log.info("  Skipping job %s — all %d result(s) already done", job_id, done)
                with _db() as conn:
                    conn.execute(
                        "UPDATE pipeline_jobs SET status='done', error=NULL WHERE id=?",
                        (job_id,)
                    )
                continue

            log.info("  Re-queuing job %s (%s) [was: %s]", job_id, target, old_status)
            # Delete non-done results to avoid duplicate rows on re-run
            with _db() as conn:
                conn.execute(
                    "DELETE FROM pipeline_results WHERE job_id=? AND status != 'done'",
                    (job_id,)
                )
                conn.execute(
                    "UPDATE pipeline_jobs SET status='queued', error=NULL WHERE id=?",
                    (job_id,)
                )
            is_snapshot = fits_dir and Path(fits_dir).name.upper() == "SNAPSHOT"
            obj_filter  = (target_filter or target) if is_snapshot else None
            _enqueue_job(int(job_id),
                         lambda jid=int(job_id), fd=fits_dir, t=target or "Unknown",
                                of=obj_filter:
                         run_pipeline(jid, fd, t, object_filter=of))
    except Exception as exc:
        log.error("Job recovery failed: %s", exc)


# Run recovery after a short delay so Flask has fully started
threading.Timer(2.0, _recover_queued_jobs).start()

# ── Flask routes ───────────────────────────────────────────────────────────────

@app.route("/health", methods=["GET"])
def health():
    db_ok = True
    try:
        with _db() as conn:
            conn.execute("SELECT 1").fetchone()
    except Exception:
        db_ok = False
    with _active_lock:
        active = list(_active_jobs)
    return jsonify({
        "status": "ok" if _worker_thread.is_alive() and db_ok else "degraded",
        "service": "processing_service",
        "queue_depth": _job_queue.qsize(),
        "current_job": _current_job_id,
        "active_jobs": active,
        "worker_alive": _worker_thread.is_alive(),
        "db_ok": db_ok,
    })


@app.route("/process", methods=["POST"])
def trigger():
    data           = request.json or {}
    job_id         = data.get("job_id")
    fits_dir       = data.get("fits_dir") or ""
    target         = data.get("target") or (Path(fits_dir).name if fits_dir else "") or "Unknown"
    selected_files = data.get("selected_files") or None  # list of abs paths, or None
    # object_filter: when set, only FITS files whose OBJECT header matches this
    # name are processed — prevents SNAPSHOT dirs from processing all targets.
    object_filter  = data.get("object_filter") or None

    if not job_id or not fits_dir:
        return jsonify({"success": False, "error": "job_id and fits_dir are required"}), 400

    if not Path(fits_dir).exists():
        return jsonify({"success": False, "error": f"Directory not found: {fits_dir}"}), 400

    enqueued = _enqueue_job(
        int(job_id),
        lambda jid=int(job_id), fd=fits_dir, t=target,
               sf=selected_files, of=object_filter:
        run_pipeline(jid, fd, t, selected_files=sf, object_filter=of)
    )
    if not enqueued:
        return jsonify({"success": True, "job_id": job_id,
                        "message": "Already active — skipped duplicate"}), 409

    queue_depth = _job_queue.qsize()
    msg = "Running" if _current_job_id is None else f"Queued (position {queue_depth})"
    return jsonify({"success": True, "job_id": job_id,
                    "message": msg, "queue_depth": queue_depth})


@app.route("/resume", methods=["POST"])
def resume():
    data       = request.json or {}
    result_id  = data.get("result_id")
    from_step  = data.get("from_step")

    if not result_id or not from_step:
        return jsonify({"success": False, "error": "result_id and from_step are required"}), 400
    if from_step not in PIPELINE_STEP_ORDER:
        return jsonify({"success": False, "error": f"Invalid step. Valid: {PIPELINE_STEP_ORDER}"}), 400

    _job_queue.put(lambda rid=int(result_id), fs=from_step:
                   resume_pipeline(rid, fs))

    queue_pos = _job_queue.qsize()
    return jsonify({"success": True, "result_id": result_id, "from_step": from_step,
                    "queue_position": queue_pos})


@app.route("/jobs/<int:job_id>/log", methods=["GET"])
def job_log(job_id: int):
    log_path = LOG_DIR / f"job_{job_id}.log"
    if not log_path.exists():
        return jsonify({"success": False, "error": "Log not found"}), 404
    return jsonify({"success": True, "log": log_path.read_text()})


# ── Background STDWeb poller ───────────────────────────────────────────────────
# Runs independently of pipeline threads — survives service restarts.
# Every 30 s it fetches the state of every pipeline_result that has a
# stdweb_task_id and hasn't reached a terminal status yet.

_STDWEB_TERMINAL = {
    "done", "completed",
    "subtraction_done",
    "inspect_failed", "failed", "error",
}

# Map raw STDWeb state → our pipeline_results.status
_STDWEB_STATE_MAP = {
    "inspect_done":     "photometry",
    "photometry_done":  "subtraction",
    "subtraction_done": "done",
    "done":             "done",
    "completed":        "done",
    "inspect_failed":   "error",
    "failed":           "error",
    "error":            "error",
}

def _fetch_stdweb_photometry(task_id: str) -> dict:
    """Scrape MJD and magnitude measurements from a STDWeb task page.

    Returns a dict with keys: filter, mjd, direct (mag/magerr), sub (mag/magerr),
    sub_ul (ul), target — any of which may be None/absent if not available.
    """
    try:
        task_json = requests.get(
            f"{STDWEB_URL}/api/tasks/{task_id}/",
            headers={"Authorization": f"Token {STDWEB_TOKEN}"},
            timeout=10,
        ).json()
    except Exception:
        task_json = {}

    try:
        html = requests.get(
            f"{STDWEB_URL}/tasks/{task_id}",
            headers={"Authorization": f"Token {STDWEB_TOKEN}"},
            timeout=20,
        ).text
    except Exception:
        html = ""

    filt   = task_json.get("config", {}).get("filter")
    target = task_json.get("config", {}).get("target")

    mjd_m     = re.search(r"MJD is ([\d.]+)", html)
    direct_m  = re.search(r"Primary target magnitude is \w+ = ([\d.]+) \+/- ([\d.]+)", html)
    sub_m     = re.search(r"Target magnitude is \w+ = ([\d.]+) \+/- ([\d.]+)", html)
    sub_ul_m  = re.search(r"Target magnitude upper limit is \w+ > ([\d.]+)", html) if not sub_m else None

    return {
        "filter":  filt,
        "target":  target,
        "mjd":     float(mjd_m.group(1))    if mjd_m    else None,
        "direct":  {"mag": float(direct_m.group(1)), "magerr": float(direct_m.group(2))} if direct_m else None,
        "sub":     {"mag": float(sub_m.group(1)),    "magerr": float(sub_m.group(2))}    if sub_m    else None,
        "sub_ul":  float(sub_ul_m.group(1)) if sub_ul_m else None,
    }


def _stdweb_poller_loop():
    """Background thread: polls STDWeb for all live tasks every 30 s.
    When a task reaches a terminal state it also fetches photometry measurements
    and persists them into pipeline_results so history is available offline.
    """
    while True:
        time.sleep(30)
        try:
            with _db() as conn:
                rows = conn.execute(
                    """SELECT id, stdweb_task_id, status, mjd
                       FROM pipeline_results
                       WHERE stdweb_task_id IS NOT NULL
                         AND status NOT IN ('done', 'error')"""
                ).fetchall()
            for (result_id, task_id, cur_status, stored_mjd) in rows:
                try:
                    r = requests.get(
                        f"{STDWEB_URL}/api/tasks/{task_id}/",
                        headers={"Authorization": f"Token {STDWEB_TOKEN}"},
                        timeout=10,
                    )
                    raw_state = r.json().get("state", "unknown")
                    new_status = _STDWEB_STATE_MAP.get(raw_state)
                    reaching_terminal = (
                        new_status in ("done", "error") and new_status != cur_status
                    )

                    if new_status and new_status != cur_status:
                        update_result(result_id, new_status, stdweb_state=raw_state,
                                      error=(raw_state if new_status == "error" else None))
                        log.info("[poller] result %s task %s: %s → %s (stdweb=%s)",
                                 result_id, task_id, cur_status, new_status, raw_state)
                    elif raw_state not in ("unknown",):
                        update_result(result_id, cur_status, stdweb_state=raw_state)

                    # Fetch and persist photometry when newly done (and not yet stored)
                    if raw_state == "subtraction_done" and stored_mjd is None:
                        try:
                            phot = _fetch_stdweb_photometry(task_id)
                            kwargs: dict = {}
                            if phot.get("mjd"):
                                kwargs["mjd"] = phot["mjd"]
                                # Derive obs_date from MJD (approx)
                                from astropy.time import Time as _ATime
                                kwargs["obs_date"] = _ATime(phot["mjd"], format="mjd").to_datetime().strftime("%Y-%m-%d")
                            if phot.get("direct"):
                                kwargs["mag_ap"]    = phot["direct"]["mag"]
                                kwargs["magerr_ap"] = phot["direct"]["magerr"]
                            if phot.get("sub"):
                                kwargs["mag_sub"]    = phot["sub"]["mag"]
                                kwargs["magerr_sub"] = phot["sub"]["magerr"]
                            elif phot.get("sub_ul") is not None:
                                kwargs["mag_sub_ul"] = phot["sub_ul"]
                            if kwargs:
                                update_result(result_id, new_status or cur_status,
                                              stdweb_state=raw_state, **kwargs)
                                log.info("[poller] stored photometry for result %s task %s: %s",
                                         result_id, task_id, kwargs)
                        except Exception as exc:
                            log.warning("[poller] photometry fetch failed result %s: %s", result_id, exc)
                except Exception as exc:
                    log.warning("[poller] result %s task %s: %s", result_id, task_id, exc)
        except Exception as exc:
            log.warning("[poller] DB error: %s", exc)


if __name__ == "__main__":
    port = int(os.environ.get("PROCESSING_PORT", 5200))
    log.info(f"Processing service starting on port {port}")
    log.info(f"DB path : {DB_PATH}")
    log.info(f"STDWeb  : {STDWEB_URL}")

    poller = threading.Thread(target=_stdweb_poller_loop, daemon=True, name="stdweb-poller")
    poller.start()
    log.info("STDWeb background poller started")

    app.run(host="0.0.0.0", port=port, debug=False)
