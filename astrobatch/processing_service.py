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
DATA_DIR     = PROJECT_ROOT / "data"
LOG_DIR      = PROJECT_ROOT / "logs"
NAS_OUTPUT   = Path(os.environ.get("NAS_OUTPUT", "/mnt/nas/input/pyl/astro/output"))
DB_PATH      = Path(os.environ.get(
    "NIGHTMANAGER_DB",
    PROJECT_ROOT / "server" / "nightmanager.db"
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


def create_result(job_id: int, filt: str, exposure: str, n_frames: int, target: str = "") -> int:
    """Insert a pipeline_results row for one filter/target. Returns the new row id."""
    try:
        with _db() as conn:
            cur = conn.execute(
                "INSERT INTO pipeline_results (job_id, filter, exposure, n_frames, status, target) "
                "VALUES (?, ?, ?, ?, 'processing', ?)",
                (job_id, filt, exposure, n_frames, target)
            )
            return cur.lastrowid
    except Exception as exc:
        log.error("DB insert result failed: %s", exc)
        return -1


def update_result(result_id: int, status: str, error: str | None = None,
                  stdweb_task_id: str | None = None, stdweb_url: str | None = None,
                  frame_previews: list | None = None):
    """Update a pipeline_results row."""
    import json as _json
    try:
        with _db() as conn:
            fields = ["status = ?", "updated_at = ?"]
            values: list = [status, _now()]
            if error is not None:
                fields.append("error = ?");         values.append(error)
            if stdweb_task_id is not None:
                fields.append("stdweb_task_id = ?"); values.append(stdweb_task_id)
            if stdweb_url is not None:
                fields.append("stdweb_url = ?");    values.append(stdweb_url)
            if frame_previews is not None:
                fields.append("frame_previews = ?"); values.append(_json.dumps(frame_previews))
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

def scan_fits(fits_dir: Path, jlog: JobLogger,
              selected_files: list[str] | None = None) -> dict:
    """
    Scan fits_dir for FITS files grouped by (object, filter, exp_str).
    If selected_files is provided, those paths are used directly — no rglob,
    no path-matching: avoids all NFS/SMB path normalisation issues.
    Returns: { (object_slug, filter, exp_str): [Path, ...] }
    """
    folder_fallback = fits_dir.name
    if folder_fallback.upper() == "SNAPSHOT":
        folder_fallback = "unknown"

    groups: dict[tuple, list] = defaultdict(list)

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
            obj = raw_obj.lower().replace(" ", "_")
            groups[(obj, filt, exp_str)].append(f)
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

    for (obj, filt, exp_str), files in groups.items():
        jlog.info(f"  Found {len(files)} frames  object={obj}  filter={filt}  exp={exp_str}s")

    return dict(groups)

# ── Work directory setup ──────────────────────────────────────────────────────

def prepare_work_dir(target: str, filt: str, exp_str: str, date_str: str,
                     files: list[Path], jlog: JobLogger,
                     clean_first: bool = False) -> Path:
    """Copy raw files into the work tree and return the folder path.

    When clean_first=True (re-process a subset), remove existing science FITS
    from the work dir so only the selected files are present for Siril.
    """
    work_folder = DATA_DIR / date_str / target / filt / f"{exp_str}s"
    work_folder.mkdir(parents=True, exist_ok=True)

    if clean_first:
        # Remove old science files; keep _calib/, _previews/ and other subdirs
        for old in work_folder.glob("*.fit*"):
            old.unlink(missing_ok=True)
        jlog.info(f"  Work dir cleaned for re-processing subset")

    existing = {f.name for f in work_folder.glob("*.fit*")}
    copied = 0
    for src in files:
        dst = work_folder / src.name
        if src.name not in existing:
            shutil.copy2(src, dst)
            copied += 1

    jlog.info(f"  Work dir: {work_folder}  ({copied} new files, {len(files)} total)")
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
            "register pp_i -interp=cu -framing=min",
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
                    "register pp_i -interp=cu -framing=min",
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
    """Run astrometry.net solve-field on res.fit. Returns True if WCS written."""
    jlog.info("  Plate solving...")
    try:
        import sys
        if str(PROJECT_ROOT) not in sys.path:
            sys.path.insert(0, str(PROJECT_ROOT))
        from astrobatch.spliter import PlateSolver
        solver = PlateSolver()
        wcs = solver.solve_field(str(res_fit))
        if wcs is None:
            jlog.error("Plate solve failed — no WCS solution")
            return False
        jlog.info("  WCS solution written to res.fit")
        return True
    except Exception as exc:
        jlog.error(f"Plate solve error: {exc}")
        return False

# ── STDWeb upload + pipeline ───────────────────────────────────────────────────

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
        payload.update({"template": "ps1"})
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
                 selected_files: list[str] | None = None):
    """Full pipeline — runs in a background thread."""
    jlog = JobLogger(job_id)
    try:
        _run_pipeline(job_id, fits_dir, target, jlog, selected_files=selected_files)
    except Exception as exc:
        log.exception("Unhandled exception in job %s", job_id)
        jlog.error(f"UNHANDLED EXCEPTION: {exc}")
        update_job(job_id, "error", error=str(exc))
    finally:
        jlog.close()


def _run_pipeline(job_id: int, fits_dir: str, target: str, jlog: JobLogger,
                  selected_files: list[str] | None = None):
    jlog.info(f"=== Job {job_id} started ===")
    jlog.info(f"Input dir : {fits_dir}")
    if selected_files:
        jlog.info(f"Selected  : {len(selected_files)} specific file(s)")

    fits_path = Path(fits_dir)
    date_str  = _detect_date(fits_path)

    darks, flats = discover_cal_frames()
    jlog.info(f"Darks: {list(darks.keys())}  Flats: {list(flats.keys())}")

    # ── Step 1: Scan ──────────────────────────────────────────────────────────
    update_job(job_id, "scanning")
    jlog.info("Step 1: Scanning FITS files...")
    groups = scan_fits(fits_path, jlog, selected_files=selected_files)
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
        result_id = create_result(job_id, filt, exp_str, len(files), target=obj)

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
        title = f"{obj.replace('_', ' ')} {filt} {exp_str}s {date_str}"
        target_clean = obj.replace('_', '').replace(' ', '')  # e.g. "AT2026FTZ"
        task_id = stdweb_upload(res_fit, title, jlog, target=target_clean)

        if not task_id:
            jlog.error("Upload failed")
            update_result(result_id, "error", error="STDWeb upload failed")
            continue

        task_url = f"{STDWEB_URL}/tasks/{task_id}/"
        update_result(result_id, "uploaded", stdweb_task_id=task_id, stdweb_url=task_url)
        # Keep the job's stdweb fields pointing to the most recent upload
        update_job(job_id, "uploading", stdweb_task_id=task_id, stdweb_url=task_url)

        # ── Step 6: Inspect ───────────────────────────────────────────────────
        jlog.info("Step 6: Triggering inspection...")
        update_result(result_id, "inspecting")
        stdweb_action(task_id, "inspect", jlog)
        state = wait_state(task_id, ["inspect_done", "failed", "error"], jlog, timeout=300)
        jlog.info(f"  Inspection state: {state}")
        if state in ("failed", "error", "timeout"):
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

# ── Flask routes ───────────────────────────────────────────────────────────────

@app.route("/health", methods=["GET"])
def health():
    return jsonify({"status": "ok", "service": "processing_service"})


@app.route("/process", methods=["POST"])
def trigger():
    data           = request.json or {}
    job_id         = data.get("job_id")
    fits_dir       = data.get("fits_dir") or ""
    target         = data.get("target") or (Path(fits_dir).name if fits_dir else "") or "Unknown"
    selected_files = data.get("selected_files") or None  # list of abs paths, or None

    if not job_id or not fits_dir:
        return jsonify({"success": False, "error": "job_id and fits_dir are required"}), 400

    if not Path(fits_dir).exists():
        return jsonify({"success": False, "error": f"Directory not found: {fits_dir}"}), 400

    thread = threading.Thread(
        target=run_pipeline,
        args=(int(job_id), fits_dir, target),
        kwargs={"selected_files": selected_files},
        daemon=True,
    )
    thread.start()

    return jsonify({"success": True, "job_id": job_id, "message": "Pipeline started"})


@app.route("/jobs/<int:job_id>/log", methods=["GET"])
def job_log(job_id: int):
    log_path = LOG_DIR / f"job_{job_id}.log"
    if not log_path.exists():
        return jsonify({"success": False, "error": "Log not found"}), 404
    return jsonify({"success": True, "log": log_path.read_text()})


if __name__ == "__main__":
    port = int(os.environ.get("PROCESSING_PORT", 5200))
    log.info(f"Processing service starting on port {port}")
    log.info(f"DB path : {DB_PATH}")
    log.info(f"STDWeb  : {STDWEB_URL}")
    app.run(host="0.0.0.0", port=port, debug=False)
