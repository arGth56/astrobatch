#!/usr/bin/env python3
"""
Analyse STDWeb transient-candidate pages for a range of task IDs.

Example usage:
    python analyse.py --begin 400 --end 450 --csv summary.csv

The script fetches each URL of the form
    http://51.178.73.109:7000/tasks/<TASK_ID>/candidates
parses the HTML table of candidates and prints a summary. If --csv is
provided, all candidate rows are written to that file for further
exploration (one row per candidate, with the originating task_id
included).
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Any

import requests
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
from astrobatch.db import CandidateDB
import json, math
import functools
import re
import os
import joblib

# Base URL of the STDWeb instance (updated – July 2025)
HOST = "http://86.253.141.183:7000"
# Candidate list page for a given task
BASE_URL = f"{HOST}/tasks/{{task_id}}/candidates"
HEADERS = {"User-Agent": "STDWeb-Analyzer/1.0 (+https://example.com)"}
MODEL_PATH = Path("candidate_model.json")
MODEL_PKL = Path("candidate_model.pkl")

HPM_CACHE: Dict[Tuple[int,int], bool] = {}

# cache for SIMBAD proximity queries
NEAR_CACHE = {}

# cache for inside galaxy queries
INSIDE_CACHE = {}

FEATURE_KEYS = [
    "mag","magerr","flux","fluxerr","fwhm","flux_radius",
    "ab_ratio","a","b","theta","bg","bg_fluxerr","bg_local",
    "xerr","yerr","number","mag_auto","isoarea_image",
    "mag_calib","mag_calib_err","flags","is_hpm",
    "near_star","near_galaxy","inside_galaxy",
    "fwhm_diff_psf","fwhm_sq"
]

def is_high_pm(ra: float, dec: float) -> bool:
    """Query SIMBAD for objects near (ra,dec); return True if any has type containing 'PM'."""
    if ra is None or dec is None:
        return False
    # cache key at 1 arcsec precision to reduce queries
    key = (int(ra*3600), int(dec*3600))
    if key in HPM_CACHE:
        return HPM_CACHE[key]
    try:
        url = (
            "https://simbad.cds.unistra.fr/simbad/sim-coo?"
            f"Coord={ra:.6f}+{dec:.6f}&Radius=5&Radius.unit=arcsec&output.format=ASCII"
        )
        resp = requests.get(url, timeout=10)
        if resp.status_code != 200:
            HPM_CACHE[key] = False
            return False
        for line in resp.text.splitlines():
            if line.startswith('#') or line.startswith('-'):
                continue
            parts = line.split('|')
            if len(parts) >= 4:
                typ = parts[3].strip()
                try:
                    dist = float(parts[1])  # dist(asec)
                except ValueError:
                    continue
                if dist <= 5 and 'PM' in typ:
                    HPM_CACHE[key] = True
                    return True
        HPM_CACHE[key] = False
        return False
    except Exception:
        HPM_CACHE[key] = False
        return False

def has_object_type(ra: float, dec: float, radius_arcsec: float, pattern: str) -> bool:
    """Query SIMBAD around (ra,dec) with given radius; return True if any object matches the desired pattern.

    The SIMBAD ASCII cone-search answer contains free-form text lines such as::

        Object 2MASS J10343999+2817489  ---  LM*  ---  OID=@9115777

    which do not contain pipe ('|') separators as originally assumed.  They are instead delimited by triple
    dashes.  We therefore match the full line with a regex instead of splitting on '|'.

    The *pattern* argument is user-friendly ("star" / "galaxy").  Internally we translate this into a regex
    that covers the various SIMBAD object type codes, e.g.  ``LM*`` (low-mass star), ``V*`` (variable star),
    ``PM*`` (high proper-motion star) or ``G``/``GiG`` for galaxies.
    """
    if ra is None or dec is None:
        return False

    wanted = pattern.lower()
    if wanted not in {"star", "galaxy"}:
        category_re = re.compile(pattern, re.I)
    else:
        category_re = None  # we will check explicitly by type code

    key = (int(ra * 3600), int(dec * 3600), int(radius_arcsec), pattern)
    if key in NEAR_CACHE:
        return NEAR_CACHE[key]

    try:
        url = (
            "https://simbad.cds.unistra.fr/simbad/sim-coo?"
            f"Coord={ra:.6f}+{dec:.6f}&Radius={radius_arcsec}&Radius.unit=arcsec&output.format=ASCII"
        )
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            NEAR_CACHE[key] = False
            return False

        for line in r.text.splitlines():
            ls = line.strip()
            if not ls or ls.startswith(("#", "-", "=")):
                continue
            if not ls.lower().startswith("object "):
                # Try summary line format: row|dist|identifier|typ|coords ...
                parts = ls.split("|")
                if len(parts) >= 4:
                    typ = parts[3].strip()
                else:
                    continue
            else:
                # Object line: split after first --- field to get type
                try:
                    typ = ls.split("---")[1].strip()
                except Exception:
                    typ = ""

            is_match = False
            if category_re is not None:
                is_match = bool(category_re.search(ls))
            else:
                if wanted == "star":
                    # Star types in SIMBAD end with '*'
                    is_match = typ.endswith("*") or typ.lower() == "star"
                else:  # galaxy
                    gal_codes = {"G", "GiG", "IG", "BCD", "LIN", "HII"}
                    typ_lower = typ.lower()
                    is_match = (
                        typ in gal_codes or
                        typ_lower.startswith(("g", "s", "sb", "sa", "sc", "sd", "sab", "sbc", "pec", "sy")) or
                        "gal" in typ_lower
                    )

            if is_match:
                NEAR_CACHE[key] = True
                # Optional debug trace
                if os.environ.get("ANALYSE_DEBUG"):
                    print(f"[has_object_type] match: '{ls}' -> {pattern}")
                return True

        NEAR_CACHE[key] = False
        if os.environ.get("ANALYSE_DEBUG"):
            print(f"[has_object_type] no {pattern} within {radius_arcsec} arcsec for RA={ra}, Dec={dec}")
        return False
    except Exception as exc:
        NEAR_CACHE[key] = False
        if os.environ.get("ANALYSE_DEBUG"):
            print(f"[has_object_type] error: {exc}")
        return False

def has_gaia_star(ra: float, dec: float, radius_arcsec: float = 5) -> bool:
    """Return True if Gaia DR3 has at least one source within radius_arcsec."""
    if ra is None or dec is None:
        return False
    key = (int(ra*3600), int(dec*3600), int(radius_arcsec), 'gaia')
    if key in NEAR_CACHE:
        return NEAR_CACHE[key]
    # Use ARI Heidelberg Gaia TAP sync service
    radius_deg = radius_arcsec / 3600.0
    adql = (
        "SELECT TOP 1 source_id FROM gaiadr3.gaia_source "
        "WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {ra}, {dec}, {radius_deg})) "
        "AND parallax_over_error >= 2"
    )
    url = "https://gaia.ari.uni-heidelberg.de/tap/sync"
    try:
        r = requests.post(
            url,
            data={
                'REQUEST': 'doQuery',
                'LANG': 'ADQL',
                'FORMAT': 'csv',
                'QUERY': adql,
            },
            timeout=10,
        )
        r.raise_for_status()
        # CSV header + possibly one data row
        lines = [l for l in r.text.splitlines() if l and not l.startswith('#')]
        is_star = len(lines) > 1  # header + at least one row
        NEAR_CACHE[key] = is_star
        return is_star
    except Exception:
        NEAR_CACHE[key] = False
        return False

def has_gaia_galaxy(ra: float, dec: float, radius_arcsec: float = 60) -> bool:
    """Return True if Gaia DR3 classifies at least one source as *galaxy* inside radius.

    We rely on the probability column `classprob_dsc_combmod_galaxy`.  A threshold
    of 0.5 (same as DSC's default winner-takes-all) is sufficient to catch most
    galaxies while keeping contamination low.
    """
    if ra is None or dec is None:
        return False

    key = (int(ra * 3600), int(dec * 3600), int(radius_arcsec), "gaia_gal")
    if key in NEAR_CACHE:
        return NEAR_CACHE[key]

    radius_deg = radius_arcsec / 3600.0
    adql = (
        "SELECT TOP 1 source_id "
        "FROM gaiadr3.gaia_source "
        "WHERE classprob_dsc_combmod_galaxy >= 0.5 "
        "AND 1 = CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {ra}, {dec}, {radius_deg})) "
    )

    url = "https://gaia.ari.uni-heidelberg.de/tap/sync"
    try:
        r = requests.post(
            url,
            data={
                "REQUEST": "doQuery",
                "LANG": "ADQL",
                "FORMAT": "csv",
                "QUERY": adql,
            },
            timeout=10,
        )
        r.raise_for_status()

        lines = [l for l in r.text.splitlines() if l and not l.startswith("#")]
        has_gal = len(lines) > 1  # header + at least one result row
        NEAR_CACHE[key] = has_gal
        return has_gal
    except Exception:
        NEAR_CACHE[key] = False
        return False

def fetch_candidates(task_id: int) -> Tuple[List[str], List[Dict[str, str]]]:
    """Return (header, rows) for the candidate table of a given task."""
    url = BASE_URL.format(task_id=task_id)
    try:
        r = requests.get(url, headers=HEADERS, timeout=15)
        if r.status_code != 200:
            print(f"⚠️  Task {task_id}: HTTP {r.status_code}")
            return [], []

        # Use pandas to parse; fallback to BeautifulSoup if pandas fails
        try:
            tables = pd.read_html(r.text)
            if not tables:
                raise ValueError("no tables")
            df = tables[0]
            header = [str(c) for c in df.columns]
            rows = df.astype(str).to_dict(orient="records")
            return header, rows
        except Exception:
            pass  # fallback

        soup = BeautifulSoup(r.text, "html.parser")
        tables = soup.find_all("table")
        if not tables:
            print(f"⚠️  Task {task_id}: no <table> found")
            return [], []

        rows = []
        header = []
        for table in tables:
            header_cells = table.find("tr").find_all(["th", "td"])
            header = [c.get_text(strip=True) or f"col{i}" for i, c in enumerate(header_cells)]
            for tr in table.find_all("tr")[1:]:
                cells = [td.get_text(strip=True) for td in tr.find_all("td")]
                if len(cells) != len(header):
                    continue
                rows.append({h: v for h, v in zip(header, cells)})
        return header, rows
    except Exception as exc:
        print(f"⚠️  Task {task_id}: failed to fetch/parse ({exc})")
        return [], []

def summarise(rows: List[Dict[str, str]]) -> str:
    """Return a one-line summary of candidate rows."""
    if not rows:
        return "0 candidates"
    mags = []
    for row in rows:
        m = row.get("mag") or row.get("MAG_AUTO")
        try:
            mags.append(float(m))
        except (TypeError, ValueError):
            pass
    if mags:
        return f"{len(rows)} candidates (brightest mag={min(mags):.2f})"
    return f"{len(rows)} candidates"

def write_csv(csv_path: Path, header: List[str], all_rows: List[Dict[str, str]]):
    """Write rows to CSV (append mode if file exists)."""
    write_header = not csv_path.exists()
    with csv_path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["task_id"] + header)
        if write_header:
            writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

def cast_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """Attempt to convert all columns to numeric when possible."""
    for col in df.columns:
        if col in ("task_id", "label"):
            continue
        df[col] = pd.to_numeric(df[col], errors="ignore")
    return df

def compare_features(df: pd.DataFrame):
    """Print simple statistics contrasting TP vs FP for numeric features."""
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    tp = df[df["label"] == "tp"]
    fp = df[df["label"] == "fp"]
    if tp.empty or fp.empty:
        print("Need non-empty TP and FP samples for comparison.")
        return

    print("\n=== Feature comparison (mean ± std) ===")
    rows = []
    for col in numeric_cols:
        rows.append((col, tp[col].mean(), fp[col].mean(), abs(tp[col].mean()-fp[col].mean())))
    rows.sort(key=lambda x: x[3], reverse=True)
    for col, tp_mean, fp_mean, diff in rows[:10]:
        print(f"{col:15s}  TP={tp_mean:.3g}  FP={fp_mean:.3g}  |Δ|={diff:.3g}")

# ---------- Simple rule-based filter & ranking ----------

def passes_rule(row: Dict[str, str]) -> Tuple[bool, float]:
    """Return (is_candidate, score) based on heuristic rules.

    score is larger for more likely real transients.
    """
    try:
        flags = row.get("flags", "0x0")
        if flags != "0x0":
            return False, 0.0

        mag = float(row.get("mag") or row.get("MAG_AUTO") or 0)
        fwhm = float(row.get("fwhm") or 0)
        fr   = float(row.get("FLUX_RADIUS") or 0)
        a = float(row.get("a") or 0)
        b = float(row.get("b") or 0)

        if mag <= 1:
            return False, 0.0
        if fwhm <= 2.5 or fr <= 1.4:
            return False, 0.0
        if b == 0:
            return False, 0.0
        ab_ratio = a / b if b else 99
        if not (1.0 < ab_ratio < 1.6):
            return False, 0.0

        # Score: higher mag (but capped), higher fwhm, roundness closeness to 1→1.3
        score = (mag - 1) + (fwhm - 2.5) + (fr - 1.4) - abs(ab_ratio - 1.2)
        return True, score
    except Exception:
        return False, 0.0

def find_cutout_url(task_id: int, row: Dict[str, str]) -> str:
    name = row.get("cutout_name")
    if not name:
        return ""
    # Construct full PNG URL (no timestamp param)
    return f"{HOST}/tasks/{task_id}/cutout/{name}?format=png"

def extract_features(row: Dict[str, str]) -> Dict[str, float]:
    def to_f(x):
        try:
            return float(x)
        except Exception:
            return 0.0

    def f(k):
        return to_f(row.get(k))

    mag = f("mag") or f("MAG_AUTO")
    a = f("a")
    b = f("b")
    ab_ratio = a / b if b else 0.0
    flags = 0 if row.get("flags") == "0x0" else 1
    # high proper motion flag
    if "is_hpm" in row:
        is_hpm = row["is_hpm"]
    else:
        try:
            is_hpm = 1.0 if is_high_pm(float(row.get("ra",0)), float(row.get("dec",0))) else 0.0
        except Exception:
            is_hpm = 0.0

    # Non-linear FWHM features
    psf_ref = 3.5  # reference PSF FWHM in pixels (adjust if needed)
    fwhm_diff_psf = abs(f("fwhm") - psf_ref)
    fwhm_sq = f("fwhm") * f("fwhm")

    ra_val = f("ra") or f("RA")
    dec_val = f("dec") or f("DEC")

    # near_star and near_galaxy flags – use stored value if available
    if row.get("near_star") not in (None, ""):
        near_star = float(row["near_star"])
    else:
        near_star = 1.0 if has_gaia_star(ra_val, dec_val, 5) else 0.0

    if row.get("near_galaxy") not in (None, ""):
        near_galaxy = float(row["near_galaxy"])
    else:
        near_galaxy = 1.0 if has_object_type(ra_val, dec_val, 60, "galaxy") else 0.0

    # inside_galaxy flag – prefer cached value
    if row.get("inside_galaxy") not in (None, ""):
        inside_galaxy = float(row["inside_galaxy"])
    else:
        inside_galaxy = 1.0 if is_inside_galaxy(ra_val, dec_val) else 0.0

    feat = {
        "mag": mag,
        "magerr": f("magerr"),
        "flux": f("flux"),
        "fluxerr": f("fluxerr"),
        "fwhm": f("fwhm"),
        "flux_radius": f("FLUX_RADIUS"),
        "ab_ratio": ab_ratio,
        "a": a,
        "b": b,
        "theta": f("theta"),
        "bg": f("bg"),
        "bg_fluxerr": f("bg_fluxerr"),
        "bg_local": f("bg_local"),
        "xerr": f("xerr"),
        "yerr": f("yerr"),
        "number": f("NUMBER"),
        "mag_auto": f("MAG_AUTO"),
        "isoarea_image": f("ISOAREA_IMAGE"),
        "mag_calib": f("mag_calib"),
        "mag_calib_err": f("mag_calib_err"),
        "flags": float(flags),
        "is_hpm": is_hpm,
        "near_star": near_star,
        "near_galaxy": near_galaxy,
        "inside_galaxy": inside_galaxy,
        "fwhm_diff_psf": fwhm_diff_psf,
        "fwhm_sq": fwhm_sq,
    }
    return feat

def load_model():
    if not MODEL_PKL.exists():
        return None
    return joblib.load(MODEL_PKL)

def predict_probability(feat: Dict[str, float], model):
    if model is None:
        return 0.0

    vec = np.array([feat.get(k, 0.0) for k in FEATURE_KEYS], dtype=float).reshape(1, -1)
    try:
        Z = model["pca"].transform(model["scaler"].transform(vec))
        prob = model["clf"].predict_proba(Z)[0, 1]
        return float(prob)
    except Exception:
        return 0.0

def _nearest_galaxy(ra: float, dec: float, search_radius_arcsec: float = 1200):
    """Return (dist_arcsec, identifier) of nearest galaxy within search radius, or None."""
    key = (int(ra*3600), int(dec*3600), int(search_radius_arcsec))
    if key in INSIDE_CACHE:
        return INSIDE_CACHE[key]
    url = (
        "https://simbad.cds.unistra.fr/simbad/sim-coo?"
        f"Coord={ra:.6f}+{dec:.6f}&Radius={search_radius_arcsec}&Radius.unit=arcsec&output.format=ASCII"
    )
    try:
        txt = requests.get(url, timeout=10).text
        min_dist = None
        gal_id = None
        for line in txt.splitlines():
            if line.startswith("#") or line.startswith("-"):
                continue
            parts = line.split("|")
            if len(parts) < 4:
                continue
            try:
                dist = float(parts[1].strip())
            except Exception:
                continue
            typ = parts[3].strip()
            if typ not in {"G", "GiG", "IG", "BCD", "LIN", "HII"} and not typ.startswith("S") and not typ.startswith("SB"):
                continue
            if min_dist is None or dist < min_dist:
                min_dist = dist
                gal_id = parts[2].strip()
        INSIDE_CACHE[key] = (min_dist, gal_id) if min_dist is not None else None
        return INSIDE_CACHE[key]
    except Exception:
        INSIDE_CACHE[key] = None
        return None

def is_inside_galaxy(ra: float, dec: float) -> bool:
    """Return True if (ra,dec) lies inside the angular size of nearest galaxy within 20'."""
    info = _nearest_galaxy(ra, dec)
    if not info:
        return False
    dist_arcsec, gal_id = info
    if gal_id is None:
        return False
    # Fetch galaxy details to get angular size
    try:
        url = f"https://simbad.cds.unistra.fr/simbad/sim-id?Ident={gal_id}&output.format=ASCII"
        txt = requests.get(url, timeout=10).text
        maj_arcsec = None
        for line in txt.splitlines():
            if line.startswith("Angular size"):
                # format: Angular size: 0.330 0.304  70 (NIR )
                tokens = line.split(":",1)[1].strip().split()
                if tokens:
                    try:
                        maj = float(tokens[0])
                        # if value <10 assume arcmin else arcsec? Here we assume arcmin.
                        maj_arcsec = maj * 60.0
                    except Exception:
                        pass
                break
        if maj_arcsec is None:
            # fallback threshold 20 arcsec
            maj_arcsec = 20.0
        return dist_arcsec <= maj_arcsec
    except Exception:
        return False

# ------------------------------------------------------------------
# Simple proximity-rule classifier
# ------------------------------------------------------------------

def rule_label(near_star: int, near_galaxy: int, inside_galaxy: int) -> str:
    """Return 'fp', 'tp' or 'unknown' from the three SIMBAD/Gaia flags.

    Logic:
      1. clear star residual  → FP  (near_star=1 & near_galaxy=0)
      2. clearly inside galaxy → TP  (inside_galaxy=1)
      3. otherwise               unknown – defer to ML model / human
    """
    try:
        ns = int(near_star)
        ng = int(near_galaxy)
        inside = int(inside_galaxy)
    except Exception:
        return "unknown"

    if ns == 1 and ng == 0:
        return "fp"
    elif inside == 1:
        return "tp"
    else:
        return "unknown"

def generate_html_report(html_path: Path, probable: List, args, current_task: int = None):
    """Generate HTML report with current progress."""
    if not probable:
        return
    
    import html
    
    # color scale based on probability 0→red, 1→green
    def prob_to_color(p: float) -> str:
        p_clamped = max(0.0, min(1.0, p))
        hue = p_clamped * 120  # 0:red 120:green
        return f"hsl({hue:.0f},70%,85%)"

    # Sort by probability if available otherwise by heuristic score
    probable_sorted = sorted(probable, key=lambda t: t[1].get('prob', 0.0), reverse=True)
    
    rows_html = []
    for rank, (sc, row, passed) in enumerate(probable_sorted, 1):
        tid = row["task_id"]
        link_page = f"{HOST}/tasks/{tid}/candidates"
        cutout = find_cutout_url(int(tid), row)
        mag = row.get("mag") or row.get("MAG_AUTO")
        prob_val = row.get('prob', 0.0)
        color = prob_to_color(prob_val)
        prob_cell = f"<td style='background:{color}'>{prob_val:.2f}</td>"
        label = row.get("label", "?")
        btns = (
            f"<button onclick=\"flag('{tid}','{row.get('cutout_name')}','tp')\">TP</button> "
            f"<button onclick=\"flag('{tid}','{row.get('cutout_name')}','fp')\">FP</button>"
        )
        ns = int(row.get('near_star', 0))
        ng = int(row.get('near_galaxy', 0))
        simbad_cell = ""
        if row.get("ra") and row.get("dec"):
            try:
                coord = f"{float(row.get('ra')):.6f}%20{float(row.get('dec')):.6f}"
                simbad_url = (
                    f"https://simbad.cds.unistra.fr/simbad/sim-fcoo?Coord={coord}&Radius=5&submit=submit+query"
                )
                style = " style='background:#faa'" if row.get("is_hpm") else ""
                simbad_cell = f"<td{style}><a href='{simbad_url}' target='_blank'>SIMBAD</a></td>"
            except Exception:
                simbad_cell = "<td></td>"

        if cutout:
            cutout_cell = (
                f"<td><a href='{cutout}' target='_blank' onmouseover=\"showImg('{cutout}')\" onmouseout=\"hideImg()\">cutout</a></td>"
            )
        else:
            cutout_cell = "<td></td>"

        rows_html.append(
            "<tr>"
            f"<td>{rank}</td>"
            f"<td><a href='{link_page}' target='_blank'>{tid}</a></td>"
            f"<td>{html.escape(mag or '')}</td>"
            f"{prob_cell}"
            f"<td>{row.get('fwhm','')}</td>"
            f"<td>{row.get('FLUX_RADIUS','')}</td>"
            f"<td>{row.get('a','')}</td>"
            f"<td>{row.get('b','')}</td>"
            f"<td>{row.get('flags','')}</td>"
            f"<td>{ns}</td><td>{ng}</td>"
            f"{simbad_cell}"
            f"{cutout_cell}"
            f"<td>{label}</td><td>{btns}</td>"
            "</tr>"
        )

    # Progress indicator
    if current_task:
        progress_info = f"Processing... Current task: {current_task} (last updated: {__import__('datetime').datetime.now().strftime('%H:%M:%S')})"
    else:
        progress_info = f"Analysis complete (finished at: {__import__('datetime').datetime.now().strftime('%H:%M:%S')})"

    html_content = f"""
<!DOCTYPE html>
<html><head>
  <meta charset='utf-8'/>
  <title>STDWeb candidate report</title>
  <style>
    table {{ border-collapse: collapse; }}
    td, th {{ border: 1px solid #ccc; padding: 4px 6px; }}
    th {{ background:#eee; cursor:pointer; user-select: none; }}
    th:hover {{ background:#ddd; }}
    .sort-asc::after {{ content: ' ▲'; }}
    .sort-desc::after {{ content: ' ▼'; }}
    .progress {{ background: #f0f0f0; padding: 10px; margin: 10px 0; border-radius: 5px; font-weight: bold; }}
  </style>
  <script>
    function sortTable(columnIndex) {{
      const table = document.querySelector('table');
      const tbody = table.querySelector('tbody');
      const rows = Array.from(tbody.querySelectorAll('tr'));
      const headers = table.querySelectorAll('th');
      const currentHeader = headers[columnIndex];
      
      // Determine sort direction based on current state
      let isAscending = true;
      if (currentHeader.classList.contains('sort-asc')) {{
        isAscending = false; // Switch to descending
      }} else if (currentHeader.classList.contains('sort-desc')) {{
        isAscending = true;  // Switch to ascending
      }}
      
      // Remove existing sort indicators from all headers
      headers.forEach(h => h.classList.remove('sort-asc', 'sort-desc'));
      
      // Sort rows
      rows.sort((a, b) => {{
        const aVal = a.children[columnIndex].textContent.trim();
        const bVal = b.children[columnIndex].textContent.trim();
        
        // Try to parse as numbers first
        const aNum = parseFloat(aVal);
        const bNum = parseFloat(bVal);
        
        let comparison;
        if (!isNaN(aNum) && !isNaN(bNum)) {{
          comparison = aNum - bNum;
        }} else {{
          comparison = aVal.localeCompare(bVal);
        }}
        
        return isAscending ? comparison : -comparison;
      }});
      
      // Update sort indicator
      currentHeader.classList.add(isAscending ? 'sort-asc' : 'sort-desc');
      
      // Re-append sorted rows
      rows.forEach(row => tbody.appendChild(row));
    }}

    function flag(task_id, cutout_name, label) {{
      const url = `http://localhost:5100/flag?task_id=${'{'}task_id{'}'}&cutout=${'{'}encodeURIComponent(cutout_name){'}'}&label=${'{'}label{'}'}`;
      fetch(url)
        .then(r => {{
          if (r.ok) {{
            alert(`Saved as ${{label}}`);
          }} else {{
            alert('Server error');
          }}
        }})
        .catch(err => alert('Fetch error'));
    }}

    function showImg(src) {{
      // Force HTTP scheme to avoid Django dev-server 400 errors when the browser
      // automatically probes HTTPS (HTTPS-first mode / speculative preconnect).
      if (src.startsWith('https://')) {{
        src = 'http://' + src.substring(8);
      }}
      const img = document.getElementById('previewImg');
      img.src = src;
      img.style.display = 'block';
    }}
    function hideImg() {{
      const img = document.getElementById('previewImg');
      img.style.display = 'none';
    }}
    
    // Add click listeners to headers when page loads
    document.addEventListener('DOMContentLoaded', function() {{
      const headers = document.querySelectorAll('th');
      headers.forEach((header, index) => {{
        header.addEventListener('click', () => sortTable(index));
      }});
    }});
  </script>
</head><body>
  <h2>Candidates from tasks {args.begin} – {args.end}</h2>
  <div class="progress">{progress_info}</div>
  <p>{len(probable_sorted)} candidates listed; color scale indicates probability (green = more likely real).</p>
  <p><em>Click any column header to sort (ascending/descending)</em></p>
  <table>
    <thead><tr><th>#</th><th>Task</th><th>Mag</th><th>Prob</th><th>FWHM</th><th>Flux_Radius</th><th>a</th><th>b</th><th>Flags</th><th>NearStar</th><th>NearGal</th><th>SIMBAD</th><th>Cutout</th><th>Label</th><th>Actions</th></tr></thead>
    <tbody>
      {''.join(rows_html)}
    </tbody>
  </table>

  <img id="previewImg" style="display:none; position:fixed; top:10px; right:10px; max-width:350px; border:1px solid #888; background:#fff; z-index:999;" />
</body></html>
"""
    html_path.write_text(html_content, encoding="utf-8")

def main():
    parser = argparse.ArgumentParser(description="Analyse STDWeb transient candidates over a task range.")
    parser.add_argument("--begin", type=int, help="First task ID (inclusive)")
    parser.add_argument("--end", type=int, help="Last task ID (inclusive)")
    parser.add_argument("--csv", type=str, help="Optional CSV output file to store all candidates")
    parser.add_argument("--tp", type=str, help="Comma-separated list of task IDs known to be true positives")
    parser.add_argument("--fp", type=str, help="Comma-separated list of task IDs known to be false positives")
    parser.add_argument("--top", type=int, default=10, help="Show N best candidates according to heuristic")
    parser.add_argument("--html", type=str, help="Optional HTML report file to generate with sortable table")
    parser.add_argument("--db", type=str, default="candidates.sqlite", help="SQLite database path")
    parser.add_argument("--set-label", type=str, choices=["tp","fp","unknown"], help="Set label for a specific candidate")
    parser.add_argument("--task", type=int, help="Task id for --set-label mode")
    parser.add_argument("--cutout", type=str, help="Cutout name for --set-label mode")
    parser.add_argument("--retrain-model", action="store_true", help="Retrain logistic model from labelled data")
    parser.add_argument("--host", type=str, help="Base URL of the STDWeb instance (e.g. http://127.0.0.1:7000)")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Determine HOST dynamically (CLI overrides env var, which overrides default)
    # ------------------------------------------------------------------
    global HOST, BASE_URL
    if args.host:
        HOST = args.host.rstrip("/")  # strip trailing slash if any
    elif os.getenv("STDWEB_HOST"):
        HOST = os.getenv("STDWEB_HOST").rstrip("/")
    # Recompute BASE_URL with possibly updated HOST
    BASE_URL = f"{HOST}/tasks/{{task_id}}/candidates"

    if args.begin is not None and args.end is not None and args.begin > args.end:
        parser.error("--begin must be <= --end")

    all_rows: List[Dict[str, str]] = []
    stored_header: List[str] = []
    probable: List[Tuple[float, Dict[str, str], bool]] = []  # (score, row, passed)

    # ------------------------------------------------------------------
    # Database setup
    # ------------------------------------------------------------------
    db = CandidateDB(Path(args.db))

    # Handle explicit label setting
    if args.set_label:
        if args.task is None or args.cutout is None:
            parser.error("--set-label requires --task and --cutout")
        # Ensure the candidate row exists; if not, insert a placeholder
        if db.get_label(args.task, args.cutout) is None:
            db.upsert({"task_id": args.task, "cutout": args.cutout})
        db.set_label(args.task, args.cutout, args.set_label)
        print(f"Label for task {args.task} cutout {args.cutout} set to {args.set_label}")
        return

    # Require begin/end for normal analysis; not needed for --set-label or --retrain-model
    if not args.retrain_model and (args.begin is None or args.end is None):
        parser.error("--begin and --end required for analysis mode (omit only when --set-label or --retrain-model)")

    tp_ids = set(int(t) for t in args.tp.split(",")) if args.tp else set()
    fp_ids = set(int(t) for t in args.fp.split(",")) if args.fp else set()

    if args.retrain_model:
        # collect labeled data
        rows = list(db.fetch_all("label IN ('tp','fp')"))
        if len(rows) < 6:
            print("Not enough labelled data to train (need at least a few TPs and FPs).")
            return
        import numpy as np
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA
        from sklearn.linear_model import LogisticRegressionCV

        X = []
        y = []
        for r in rows:
            rd = dict(r)
            feat = extract_features(rd)
            vec = [feat.get(k,0.0) for k in FEATURE_KEYS]
            if any((v is None) or math.isnan(v) for v in vec):
                continue  # skip rows with missing data
            X.append(vec)
            y.append(1 if rd.get("label") == "tp" else 0)
        if not X:
            print("No valid rows without NaNs; cannot train model.")
            return

        X = np.asarray(X, dtype=float)
        scaler = StandardScaler().fit(X)
        Xs = scaler.transform(X)

        pca = PCA(n_components=0.95, svd_solver="full").fit(Xs)
        Z = pca.transform(Xs)

        clf = LogisticRegressionCV(cv=5, class_weight="balanced", max_iter=2000)
        clf.fit(Z, y)

        joblib.dump({"scaler": scaler, "pca": pca, "clf": clf}, MODEL_PKL)
        print("Model retrained and saved to", MODEL_PKL)
        return

    for task_id in range(args.begin, args.end + 1):
        header, rows = fetch_candidates(task_id)
        print(f"Task {task_id}: {summarise(rows)}")
        if rows:
            # Ensure every row has consistent header and include task_id
            for r in rows:
                r["task_id"] = task_id
                if task_id in tp_ids:
                    r["label"] = "tp"
                elif task_id in fp_ids:
                    r["label"] = "fp"

                # Determine high proper motion flag
                try:
                    ra_val = float(r.get("ra") or r.get("RA") or 0)
                    dec_val = float(r.get("dec") or r.get("DEC") or 0)
                except Exception:
                    ra_val = dec_val = 0.0
                r["is_hpm"] = 1.0 if is_high_pm(ra_val, dec_val) else 0.0

                # ------------------------------------------------------------------
                # Feature extraction (inc. SIMBAD proximity) & model probability
                # ------------------------------------------------------------------
                feat = extract_features(r)

                # Predict probability with the ML model (no hand-coded rules)
                model = load_model()
                if model:
                    r["prob"] = predict_probability(feat, model)

                passed, sc = passes_rule(r)
                r["score"] = sc

                # Persist SIMBAD proximity flags back onto the row so that they are
                # later written to the SQLite database via cand_dict.
                r["near_star"] = feat.get("near_star", 0)
                r["near_galaxy"] = feat.get("near_galaxy", 0)
                r["inside_galaxy"] = feat.get("inside_galaxy", 0)

                # Upsert into DB
                cand_dict = {
                    "task_id": int(task_id),
                    "cutout": r.get("cutout_name", ""),
                    "ra": float(r.get("ra", 0) or 0),
                    "dec": float(r.get("dec", 0) or 0),
                    "mag": float(r.get("mag") or r.get("MAG_AUTO") or 0),
                    "score": sc,
                    "flags": r.get("flags"),
                    "fwhm": float(r.get("fwhm", 0) or 0),
                    "flux_radius": float(r.get("FLUX_RADIUS", 0) or 0),
                    "a": float(r.get("a", 0) or 0),
                    "b": float(r.get("b", 0) or 0),
                    "is_hpm": int(r.get("is_hpm", 0)),
                    "magerr": float(r.get("magerr", 0) or 0),
                    "flux": float(r.get("flux", 0) or 0),
                    "fluxerr": float(r.get("fluxerr", 0) or 0),
                    "x": float(r.get("x", 0) or 0),
                    "y": float(r.get("y", 0) or 0),
                    "xerr": float(r.get("xerr", 0) or 0),
                    "yerr": float(r.get("yerr", 0) or 0),
                    "theta": float(r.get("theta", 0) or 0),
                    "bg": float(r.get("bg", 0) or 0),
                    "number": int(r.get("NUMBER", 0) or 0),
                    "mag_auto": float(r.get("MAG_AUTO", 0) or 0),
                    "isoarea_image": float(r.get("ISOAREA_IMAGE", 0) or 0),
                    "bg_fluxerr": float(r.get("bg_fluxerr", 0) or 0),
                    "bg_local": float(r.get("bg_local", 0) or 0),
                    "mag_calib": float(r.get("mag_calib", 0) or 0),
                    "mag_calib_err": float(r.get("mag_calib_err", 0) or 0),
                    "mag_filter_name": r.get("mag_filter_name"),
                    "mag_color_name": r.get("mag_color_name"),
                    "mag_color_term": r.get("mag_color_term"),
                    "near_star": float(r.get("near_star", 0) or 0),
                    "near_galaxy": float(r.get("near_galaxy", 0) or 0),
                    "inside_galaxy": float(r.get("inside_galaxy", 0) or 0),
                    "fwhm_diff_psf": float(r.get("fwhm_diff_psf", 0) or 0),
                    "fwhm_sq": float(r.get("fwhm_sq", 0) or 0),
                }
                if cand_dict["cutout"]:
                    db.upsert(cand_dict)

                # Fetch label if exists
                label = db.get_label(int(task_id), r.get("cutout_name", ""))
                if label:
                    r["label"] = label

                probable.append((sc, r, passed))
            all_rows.extend(rows)
            if not stored_header:
                stored_header = header
        
        # Update HTML report after each task if requested
        if args.html and probable:
            generate_html_report(Path(args.html).expanduser(), probable, args, current_task=task_id)
            print(f"📝 Updated HTML report ({len(probable)} candidates so far)")
    
    # Final HTML report update (remove "processing" indicator)
    if args.html and probable:
        generate_html_report(Path(args.html).expanduser(), probable, args)
        print(f"📝 Final HTML report written to {args.html}")

    if args.csv and all_rows:
        csv_path = Path(args.csv).expanduser()
        write_csv(csv_path, stored_header, all_rows)
        print(f"📄 Wrote {len(all_rows)} rows to {csv_path}")

    # Analyse if labels provided
    if any(r.get("label") for r in all_rows):
        df = pd.DataFrame(all_rows)
        df = cast_numeric(df)
        compare_features(df)

    # Show top candidates
    if probable:
        # Sort by probability if available otherwise by heuristic score
        probable.sort(key=lambda t: t[1].get('prob', 0.0), reverse=True)
        print(f"\n=== Top {min(args.top, len(probable))} candidates (rule pass marked ✔) ===")
        for rank, (sc, row, passed) in enumerate(probable[: args.top], 1):
            tid = row["task_id"]
            mag = row.get("mag") or row.get("MAG_AUTO")
            link = find_cutout_url(int(tid), row)
            flag = "✔" if passed else "✖"
            if 'prob' in row:
                prob_str = f" p={row['prob']:.2f}"
            else:
                prob_str = ""
            print(f"#{rank:2d} {flag} task {tid} mag={mag} score={sc:.2f}{prob_str} a={row.get('a')} b={row.get('b')} fwhm={row.get('fwhm')} → {link}")

    # HTML report generation is now handled incrementally in the loop above

    if not all_rows:
        print("No candidates found in the specified range.")

if __name__ == "__main__":
    main() 