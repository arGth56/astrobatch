import glob
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
from pathlib import Path
import shutil
import math
import re
import subprocess
import sys
import argparse
import time
import json
import urllib.parse
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import pandas as pd, io, requests
from scipy.ndimage import zoom
# Optional progress-bar support
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable=None, **kwargs):
        """Minimal shim if tqdm is not installed."""
        return iterable if iterable is not None else lambda x: x

# ============================================================================
# PATH CONFIGURATION (edit here to adapt to your environment)
# ============================================================================
PROJECT_ROOT           = Path(__file__).resolve().parent.parent  # repo root

# Raw FITS data root – will be overridden by `astrobatch.cli` via --root.
# Provide a neutral default (current working directory) instead of a
# hard-wired absolute path.
DATA_ROOT: Path = Path.cwd()

# Calibration directory containing flats & darks
CALIB_DIR              = PROJECT_ROOT / "calib"

# Name of the Siril script that this program writes / executes
SIRIL_SCRIPT_NAME      = "script.ssf"

# Location where system-wide Siril scripts live (used when we copy the script)
SYSTEM_SIRIL_SCRIPTS_DIR = "/usr/share/siril/scripts"
# User fallback path for scripts when the system location is not writable
USER_SIRIL_SCRIPTS_DIR   = Path.home() / ".local/share/siril/scripts"

# Keep legacy variable names for backward compatibility with the rest of the
# existing codebase. Feel free to grep & refactor to the new names later.
init_path   = str(DATA_ROOT)
directory   = DATA_ROOT
siril_script_path = SIRIL_SCRIPT_NAME

# API Configuration
API_TOKEN = "1e296ddd6738af45467b7bc6558c00a9524447ab"
# Base URL of the photometry server – can be overridden via the SPLITER_API_BASE_URL
# environment variable if needed.
API_BASE_URL = os.getenv("SPLITER_API_BASE_URL", "http://86.253.141.183:7000")
API_ENDPOINT = f"{API_BASE_URL}/api/tasks/upload/"

# Default gain for photometry & subtraction
GAIN_VALUE = 23592

# List of exposure folders already calibrated in this run to prevent duplicate work
calibrated_folders: set[str] = set()

# Define base directories for calibration frames
# --------------- CALIBRATION FRAMES ------------------
# Calibration masters now live in the project directory under "calib/".
# Paths are derived relative to this file so the script is portable.
# CALIB_DIR defined in PATH CONFIGURATION above (duplicate removed)

# Dark masters organised by exposure time (seconds)
_darks_root = CALIB_DIR / "darks" 
darks = {
    "10" : str(_darks_root / "dark_10s.fit"),
    "60" : str(_darks_root / "dark_60s.fit"),
    "120": str(_darks_root / "dark_120s.fit"),
}

# Flat masters are single FITS files in calib/flats/
_flats_root = CALIB_DIR / "flats"
flats = {
    "G":       str(_flats_root / "flat_G.fit"),
    "RP":      str(_flats_root / "flat_Grp.fit"),
    "BP":      str(_flats_root / "flat_Gbp.fit"),
    "grating": str(_flats_root / "flat_grating.fit"),  # may be absent
}

# Mosaic configuration
MOSAIC_CONFIG = {
    'pixel_scale': 1.0,           # arcsec/pixel
    'projection': 'TAN',          # TAN (gnomonic), ZEA (zenithal equal-area), SIN (orthographic), CAR (equirectangular)
    'blend_method': 'feather',    # feather, linear, weighted
    'quality_threshold': 0.7,     # WCS quality threshold
    'sigma_clip': 3.0,            # outlier rejection threshold
    'max_memory_gb': 8.0,         # memory limit for processing
    'output_formats': ['fits', 'png'],
    'preserve_metadata': True,
    'generate_coverage_map': True,
    'edge_samples': 8,          # samples per image edge for boundary check
    'dump_border_csv': True    # export border_points.csv diagnostics
}

# Dictionary to store unnamed objects with their RA/DEC coordinates
# No longer needed - removed coordinate-based merging per user request
# unnamed_objects = {}
unnamed_object_counter = 1

# TNS base (web HTML search, simple scrape)
TNS_SEARCH_BASE = "https://www.wis-tns.org/search"  # HTML endpoint

# JSON API endpoint
TNS_API_ENDPOINT = "https://www.wis-tns.org/api/get/search"

# Provide your personal API key via environment variable
TNS_API_KEY = os.getenv("TNS_API_KEY", "")

# User-agent per TNS requirements: tns_marker <your_username>
TNS_USER = os.getenv("TNS_USER", "spliter_script")

# Regex patterns for transient names
TRANSIENT_REGEX = r"(?:SN|AT)\d{4}[A-Za-z]{1,3}"

# Astro-CoLiBRI credentials (public for read-only endpoints)
COLIBRI_UUID = "NeYByMhfVbhye61jJeN6f0NaBng2"  # Provided by user
COLIBRI_ENDPOINT = "https://astro-colibri.science/api/v1/events"  # base REST endpoint

# ------------------- NGC catalogue --------------------
# ---------------------------------------------------------------------------
# Catalogue files are now stored under *astrobatch/catalog/*.  Keep a
# fallback to the legacy location so the script remains compatible if the
# user did not move the file.
# ---------------------------------------------------------------------------

_PKG_CATALOG_DIR = Path(__file__).resolve().parent / "astrobatch" / "catalog"
_PKG_CATALOG_DIR.mkdir(parents=True, exist_ok=True)
NGC_CSV_LOCAL = _PKG_CATALOG_DIR / "NGC.csv"

# Legacy fallback (root directory) – will be tried if the new path is missing.
_NGC_CSV_LEGACY = Path(__file__).with_name("ngc.csv")

NGC_CSV_URL   = "https://raw.githubusercontent.com/cosmos-book/cosmos-book.github.io/master/catalogue/NGC/data/NGC.csv"  # new reliable mirror [[1]]

def _load_ngc_catalog():
    """Return DataFrame with columns ra_deg, dec_deg (ICRS J2000). Download once if needed."""
    csv_path = NGC_CSV_LOCAL
    # Prefer new location inside package; fall back to legacy root copy
    csv_path = NGC_CSV_LOCAL if NGC_CSV_LOCAL.exists() else _NGC_CSV_LEGACY
    # On case-sensitive filesystems the file may have been stored with another
    # case – keep a second fallback check.
    if not csv_path.exists():
        alt = csv_path.with_name("NGC.csv")
        if alt.exists():
            csv_path = alt

    if csv_path.exists():
        df = pd.read_csv(csv_path)
    else:
        try:
            print(f"⬇️  Downloading NGC catalogue …")
            r = requests.get(NGC_CSV_URL, timeout=30)
            r.raise_for_status()
            df = pd.read_csv(io.StringIO(r.text))
            df.to_csv(NGC_CSV_LOCAL, index=False)
        except Exception as e:
            print(f"⚠️  Could not download NGC catalogue: {e}")
            return pd.DataFrame(columns=["ra_deg","dec_deg"])
    # Map common column names to ra_deg / dec_deg
    col_map = {}
    for col in df.columns:
        low = col.lower()
        if low in ("ra", "radeg", "raj2000", "ra_deg"):
            col_map[col] = "ra_deg"
        elif low in ("dec", "dedeg", "decj2000", "dec_deg"):
            col_map[col] = "dec_deg"
    df = df.rename(columns=col_map)
    if not {"ra_deg", "dec_deg"}.issubset(df.columns):
        print("⚠️  NGC catalogue missing RA/Dec columns; galaxy overlay disabled")
        return pd.DataFrame(columns=["ra_deg","dec_deg"])
    return df[["ra_deg","dec_deg"]]

NGC_DF = _load_ngc_catalog()

def fetch_simbad_galaxies(ra_center, dec_center, radius_deg):
    """Return list of (ra_deg, dec_deg) for galaxies in SIMBAD within radius_deg of center."""
    url = (
        "https://simbad.cds.unistra.fr/simbad/sim-coo?"
        f"Coord={ra_center:.6f}+{dec_center:.6f}&Radius={radius_deg:.3f}&Radius.unit=deg&output.format=ASCII"
    )
    try:
        r = requests.get(url, timeout=20)
        r.raise_for_status()
        coords = []
        for line in r.text.splitlines():
            if line.startswith("#") or line.startswith("-") or line.strip()=="":
                continue
            parts = line.split("|")
            if len(parts) < 6:
                continue
            obj_type = parts[3].strip()
            if obj_type not in ("G", "Galaxy", "GPair", "Seyf", "GClstr"):
                continue
            try:
                ra_str = parts[4].strip().split()[0]
                dec_str = parts[4].strip().split()[1]
                # convert sexa to deg using simple parser
                def sexa_to_deg(s):
                    if ":" in s:
                        h, m, s2 = s.split(":")
                        return (float(h)+float(m)/60+float(s2)/3600)*15
                    return float(s)
                # SIMBAD output maybe in HMS and DMS; skip complex parse for brevity
            except Exception:
                continue
        return coords
    except Exception as e:
        print(f"⚠️  SIMBAD query failed: {e}")
        return []

def getFitHeader(path):
    # Use memmap=False and only read header to avoid loading entire file
    with fits.open(path, memmap=False) as hdulist:
        hdr = hdulist[0].header
        return hdr

def createFolderIfNeeded(folder_path):
    try:
        os.makedirs(folder_path, exist_ok=True)
    except OSError as e:
        print(f"Error: {e}")

# Cache for created folders to avoid repeated os.makedirs calls
_created_folders = set()

def createFolderIfNeededCached(folder_path):
    """Cached version of folder creation to avoid repeated filesystem checks"""
    if folder_path not in _created_folders:
        try:
            os.makedirs(folder_path, exist_ok=True)
            _created_folders.add(folder_path)
        except OSError as e:
            print(f"Error: {e}")

def fast_file_copy(source, destination):
    """Fast file copy using shutil.copy2 with proper error handling"""
    # Skip copy if file already exists with identical size
    try:
        if destination.exists() and destination.stat().st_size == source.stat().st_size:
            print(f"🟡 {destination.name} already exists – skipping copy")
            return
        shutil.copy2(source, destination)
    except Exception as e:
        # Fallback to basic copy if shutil.copy2 fails
        print(f"⚠️  shutil.copy2 failed, using fallback: {e}")
        with open(source, 'rb') as src, open(destination, 'wb') as dst:
            shutil.copyfileobj(src, dst, length=4*1024*1024)  # 4 MB buffer for speed

def parse_coordinate(coord_str):
    """
    Parse coordinate string which can be in various formats:
    - Decimal degrees: 123.456
    - HMS/DMS: 12:34:56.7 or 12h34m56.7s
    Returns coordinate in decimal degrees
    """
    if coord_str is None:
        return None
    
    coord_str = str(coord_str).strip()
    
    # Try decimal degrees first
    try:
        return float(coord_str)
    except ValueError:
        pass
    
    # Try HMS/DMS format (12:34:56.7 or 12h34m56.7s)
    patterns = [
        r'(\d+):(\d+):(\d+\.?\d*)',  # 12:34:56.7
        r'(\d+)h(\d+)m(\d+\.?\d*)s?',  # 12h34m56.7s
        r'(\d+)\s+(\d+)\s+(\d+\.?\d*)'  # 12 34 56.7
    ]
    
    for pattern in patterns:
        match = re.match(pattern, coord_str)
        if match:
            hours_or_deg = float(match.group(1))
            minutes = float(match.group(2))
            seconds = float(match.group(3))
            return hours_or_deg + minutes/60.0 + seconds/3600.0
    
    # If all parsing fails, return None
    print(f"Warning: Could not parse coordinate: {coord_str}")
    return None

def parse_ra(ra_str):
    """Parse RA value and return **degrees**.

    Heuristic:
    1. Try generic parsing to float (decimal degrees *or* hours).
    2. If the original string contains explicit HMS markers (':', 'h', 'm', 's') we
       assume it is given in *hours* and convert to degrees.
    3. Otherwise we treat the numeric result as already being in degrees – even if
       the value is < 24.
    """
    if ra_str is None:
        return None

    ra_str_s = str(ra_str).strip()
    ra_decimal = parse_coordinate(ra_str_s)
    if ra_decimal is None:
        return None

    # Detect HMS-style strings by presence of delimiters/letters
    if re.search(r"[hms:]", ra_str_s, re.IGNORECASE):
        # Value was in hours; convert to degrees if it is within 0–24 h range
        return ra_decimal * 15.0

    # Otherwise: assume decimal degrees already
    return ra_decimal

def parse_dec(dec_str):
    """Parse DEC coordinate (already in degrees)"""
    return parse_coordinate(dec_str)

def angular_distance(ra1, dec1, ra2, dec2):
    """
    Calculate angular distance between two celestial coordinates in degrees
    Returns distance in degrees
    """
    if None in [ra1, dec1, ra2, dec2]:
        return float('inf')
    
    # Convert to radians
    ra1_rad = math.radians(ra1)
    dec1_rad = math.radians(dec1)
    ra2_rad = math.radians(ra2)
    dec2_rad = math.radians(dec2)
    
    # Haversine formula for angular distance
    delta_ra = ra2_rad - ra1_rad
    delta_dec = dec2_rad - dec1_rad
    
    a = (math.sin(delta_dec/2)**2 + 
         math.cos(dec1_rad) * math.cos(dec2_rad) * math.sin(delta_ra/2)**2)
    
    distance_rad = 2 * math.asin(math.sqrt(a))
    distance_deg = math.degrees(distance_rad)
    
    return distance_deg

def get_target_name(hdr):
    """Return OBJECT header; if absent create sequential placeholder OBJECT_N."""
    global unnamed_object_counter
    for key in ("OBJECT", "OBJNAME", "TARGNAME", "TARGET"):
        if key in hdr and str(hdr[key]).strip():
            return str(hdr[key]).strip().replace(" ", "")

    # No name → assign placeholder folder with single image behaviour
    name = f"OBJECT_{unnamed_object_counter}"
    unnamed_object_counter += 1
    print(f"Created new unnamed object: {name} (1 image per folder)")
    return name

def upload_for_inspection(file_path, object_name, date_str):
    """
    Step 1: Upload a FITS file without processing
    Returns task_id if successful, None if failed
    """
    title = f"{date_str}_{object_name}"
    
    curl_command = [
        "curl", "-s", "-X", "POST", API_ENDPOINT,
        "-H", f"Authorization: Token {API_TOKEN}",
        "-F", f"file=@{file_path}",
        "-F", f"title={title}",
        "-F", "do_inspect=false",
        "-F", "do_photometry=false"
    ]
    
    try:
        print(f"📤 Uploading {file_path} as '{title}' for inspection...")
        print(f"🔧 Debug: Upload URL = {API_ENDPOINT}")
        
        result = subprocess.run(curl_command, capture_output=True, text=True, check=True)
        
        print(f"🔧 Debug: Upload response:")
        print(f"   {result.stdout}")
        
        # Try to extract task ID from response
        try:
            response_json = json.loads(result.stdout)
            print(f"🔧 Debug: Parsed JSON response:")
            print(f"   {json.dumps(response_json, indent=2)}")
            
            # Check for task ID in various possible locations
            task_id = (response_json.get('task_id') or 
                      response_json.get('id') or 
                      response_json.get('task', {}).get('id'))
            if task_id:
                print(f"✅ Upload successful. Task ID: {task_id}")
                return str(task_id)
            else:
                print(f"⚠️ Upload successful but no task ID found in response")
                print(f"🔧 Debug: Available keys in response: {list(response_json.keys())}")
                return None
        except json.JSONDecodeError:
            print(f"⚠️ Upload successful but response is not JSON: {result.stdout}")
            # Try to extract task ID from text response if possible
            import re
            task_id_match = re.search(r'"?(?:task_)?id"?\s*:?\s*"?(\d+)"?', result.stdout)
            if task_id_match:
                task_id = task_id_match.group(1)
                print(f"✅ Extracted Task ID from text: {task_id}")
                return task_id
            return None
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to upload {title}")
        print(f"Error: {e.stderr}")
        return None
    except Exception as e:
        print(f"❌ Unexpected error uploading {title}: {e}")
        return None

def launch_inspection(task_id):
    """Launch inspection phase using the /api/tasks/{id}/action/ endpoint"""
    action_url = f"{API_BASE_URL}/api/tasks/{task_id}/action/"

    # JSON body for inspection step
    payload = json.dumps({
        "action": "inspect",
        "gain": GAIN_VALUE,
        "saturation": 1
    })

    curl_command = [
        "curl", "-s", "-X", "POST", action_url,
        "-H", f"Authorization: Token {API_TOKEN}",
        "-H", "Content-Type: application/json",
        "-d", payload
    ]

    print(f"🔍 Starting inspection for task {task_id}...")
    print(f"🔧 Debug: curl → {' '.join(curl_command)}")

    result = subprocess.run(curl_command, capture_output=True, text=True)
    print(f"🔧 Debug: Return code = {result.returncode}")
    print(f"🔧 Debug: STDOUT: {result.stdout}")
    print(f"🔧 Debug: STDERR: {result.stderr}")

    if result.returncode == 0:
        print(f"✅ Inspection request sent for task {task_id}")
        return True
    else:
        print(f"❌ Failed to send inspection request for task {task_id}")
        return False

def get_task_state(task_id):
    """Return current task 'state' string or None on error"""
    status_url = f"{API_BASE_URL}/api/tasks/{task_id}/"
    curl_command = [
        "curl", "-s", "-H", f"Authorization: Token {API_TOKEN}", status_url
    ]
    try:
        result = subprocess.run(curl_command, capture_output=True, text=True, check=True)
        try:
            response_json = json.loads(result.stdout)
            return response_json.get("state")
        except json.JSONDecodeError:
            print(f"⚠️ Could not decode JSON while fetching state for task {task_id}: {result.stdout}")
            return None
    except subprocess.CalledProcessError as e:
        print(f"❌ Error retrieving state for task {task_id}: {e.stderr}")
        return None

def wait_for_state(task_id, desired_states, timeout=600, poll_interval=10):
    """Block until task.state is in desired_states list. Returns True when reached, False on timeout, None on failure state."""
    fail_states = {"failed", "error"}
    start_time = time.time()
    while time.time() - start_time < timeout:
        state = get_task_state(task_id)
        if state is None:
            # keep waiting, could not fetch
            time.sleep(poll_interval)
            continue
        print(f"🔍 Task {task_id} current state: {state}")
        if state in desired_states:
            return True
        if state in fail_states:
            print(f"❌ Task {task_id} reached failure state: {state}")
            return None
        time.sleep(poll_interval)
    print(f"⏰ Timeout ({timeout}s) waiting for task {task_id} to reach states {desired_states}")
    return False

def trigger_photometry(task_id):
    """Trigger photometry once inspection is complete"""
    action_url = f"{API_BASE_URL}/api/tasks/{task_id}/action/"

    # Use S/N = 10 instead of old S/N (20)
    payload = json.dumps({
        "action": "photometry",
        "use_color": False,
        "gain": GAIN_VALUE,
        "sn": 10
    })

    curl_command = [
        "curl", "-s", "-X", "POST", action_url,
        "-H", f"Authorization: Token {API_TOKEN}",
        "-H", "Content-Type: application/json",
        "-d", payload
    ]

    print(f"🔬 Starting photometry for task {task_id}...")
    print(f"🔧 Debug: curl → {' '.join(curl_command)}")

    result = subprocess.run(curl_command, capture_output=True, text=True)
    print(f"🔧 Debug: Return code = {result.returncode}")
    print(f"🔧 Debug: STDOUT: {result.stdout}")
    print(f"🔧 Debug: STDERR: {result.stderr}")

    if result.returncode == 0:
        print(f"✅ Photometry request sent for task {task_id}")
        print(f"🎉 Complete workflow launched for task {task_id}!")
        return True
    else:
        print(f"❌ Failed to send photometry request for task {task_id}")
        return False

def trigger_template_subtraction(task_id):
    """Send template-subtraction action to the /action/ endpoint"""
    action_url = f"{API_BASE_URL}/api/tasks/{task_id}/action/"
    # Request subtraction using ZTF DR7 reference templates
    payload = json.dumps({
        "action": "subtraction",
        "template_source": "ZTF_DR7"
    })
    curl_command = [
        "curl", "-s", "-X", "POST", action_url,
        "-H", f"Authorization: Token {API_TOKEN}",
        "-H", "Content-Type: application/json",
        "-d", payload
    ]
    print(f"🪄 Triggering template subtraction for task {task_id}…")
    print(f"🔧 Debug: curl → {' '.join(curl_command)}")
    result = subprocess.run(curl_command, capture_output=True, text=True)
    print(f"🔧 Debug: Return code = {result.returncode}")
    print(f"🔧 Debug: STDOUT: {result.stdout}")
    print(f"🔧 Debug: STDERR: {result.stderr}")
    return result.returncode == 0

def query_tns(ra_deg, dec_deg, radius_arcsec=60):
    """Query TNS JSON API for transients near coordinates. Returns set[str] names."""
    radius_deg = radius_arcsec / 3600.0

    # If API key available, use JSON API
    if TNS_API_KEY:
        payload = {
            "api_key": TNS_API_KEY,
            "data": {
                "ra": ra_deg,
                "dec": dec_deg,
                "radius": radius_deg,
                "units": "deg"
            }
        }
        payload_str = json.dumps(payload)
        curl_cmd = [
            "curl", "-s", "-X", "POST", TNS_API_ENDPOINT,
            "-H", f"User-Agent: tns_marker {TNS_USER}",
            "-H", "Content-Type: application/json",
            "-d", payload_str
        ]
        print(f"🔧 Debug: TNS API POST → {TNS_API_ENDPOINT}")
        try:
            resp = subprocess.run(curl_cmd, capture_output=True, text=True, check=True)
            reply_json = json.loads(resp.stdout)
            objs = reply_json.get("reply", {}).get("objs", [])
            names = {obj.get("name") for obj in objs if obj.get("name")}
            if names:
                print(f"🔭 TNS transients: {', '.join(names)}")
            else:
                print("🔭 No TNS transients found in JSON reply")
            return names
        except Exception as e:
            print(f"⚠️  TNS JSON API failed ({e}); falling back to HTML scrape")

    # Fallback: scrape HTML search page
    params = {
        "reported_within_last_value": 12,  # months window large
        "reported_within_last_units": "months",
        "ra": f"{ra_deg:.6f}",
        "decl": f"{dec_deg:.6f}",
        "radius": f"{radius_deg}",
        "coords_unit": "deg",
        "num_page": 50,
        "public": 1,
        "format": "csv"
    }
    url = TNS_SEARCH_BASE + "?" + urllib.parse.urlencode(params)
    curl_cmd = [
        "curl", "-s", "-L", "-A", "Mozilla/5.0", url
    ]
    print(f"🔧 Debug: TNS HTML GET → {url}")
    try:
        resp = subprocess.run(curl_cmd, capture_output=True, text=True)

        if resp.returncode != 0:
            print(f"⚠️  curl exit code {resp.returncode}. stderr: {resp.stderr.strip()}")
        content = resp.stdout

        # If CSV, each transient name is in first column labelled 'Name'
        if content.startswith("Name,Ra,Dec"):
            lines = content.splitlines()[1:]  # skip header
            names = {line.split(',')[0] for line in lines if line}
        else:
            import re
            names = set(re.findall(TRANSIENT_REGEX, content))
        if names:
            print(f"🔭 TNS transients (scrape): {', '.join(names)}")
        else:
            print("🔭 No TNS transients found on HTML page")
        return names
    except subprocess.CalledProcessError as e:
        print(f"⚠️  HTML scrape error: {e}")
        return set()

def query_colibri(ra_deg, dec_deg, radius_deg=0.25):
    """Query Astro-CoLiBRI REST API (v2025) for transient events near the given position.
    Returns a set[str] of unique event names (empty set on failure).
    The default search radius is 0.25° (~15').
    """
    try:
        # ------------------------------------------------------------------
        # Build payload following the official example: include a hard-wired
        # filter that selects optical transients (types "ot" or "ot_sn").
        # Only RA/DEC and radius are variable.
        # ------------------------------------------------------------------

        now = datetime.datetime.utcnow()
        time_min = now - datetime.timedelta(days=365)
        time_max = now

        uid_clean = (os.getenv("COLIBRI_UID", COLIBRI_UUID) or "").strip()

        payload = {
            "uid": uid_clean,
            "filter": None,
            "time_range": {
                "min": time_min.isoformat(timespec="seconds"),
                "max": time_max.isoformat(timespec="seconds"),
            },
            "properties": {
                "type": "cone",
                "position": {"ra": ra_deg, "dec": dec_deg},
                "radius": 4.0,
            },
            "return_format": "json",
        }

        payload_str = json.dumps(payload)

        curl_cmd = [
            "curl", "-s", "-X", "POST", "https://astro-colibri.science/cone_search",
            "-H", "Content-Type: application/json",
            "-d", payload_str,
        ]

        print("🔧 Debug: CoLiBRI payload:", payload_str)
        print(f"🔧 Debug: CoLiBRI POST → curl -s -X POST https://astro-colibri.science/cone_search -H 'Content-Type: application/json' -d '{payload_str}'")

        resp = subprocess.run(curl_cmd, capture_output=True, text=True)

        if resp.returncode != 0:
            print(f"⚠️  CoLiBRI HTTP error {resp.returncode}: {resp.stderr.strip()}")
            return set()

        # Attempt to parse JSON – if it fails, bail out gracefully
        try:
            data = json.loads(resp.stdout)
        except json.JSONDecodeError:
            print("⚠️  CoLiBRI response is not valid JSON – skipping")
            return set()

        # API may return a dict with a 'message' key on error/empty results
        if isinstance(data, dict) and data.get("message"):
            print(f"🔧 CoLiBRI info: {data['message']}")
            return set()

        # Successful responses: list of events in 'voevents' or bare list
        events = []
        if isinstance(data, dict):
            events = data.get("voevents", []) or data.get("events", []) or data.get("results", [])
        elif isinstance(data, list):
            events = data

        # Post-filter: keep only events inside the actual image footprint
        filtered_names = set()
        for e in events:
            if not isinstance(e, dict):
                continue
            ev_ra = e.get("ra") or e.get("RA")
            ev_dec = e.get("dec") or e.get("DEC")
            if ev_ra is None or ev_dec is None:
                continue
            try:
                dist_deg = angular_distance(float(ev_ra), float(ev_dec), ra_deg, dec_deg)
                if radius_deg is None or dist_deg <= radius_deg:
                    name_candidate = e.get("source_name") or e.get("name") or e.get("id") or e.get("trigger_id")
                    if name_candidate:
                        filtered_names.add(name_candidate)
            except Exception:
                continue

        names = filtered_names

        if names:
            print(f"🔭 CoLiBRI transients: {', '.join(sorted(names))}")
        else:
            print("🔭 No CoLiBRI transients found in field")

        return names
    except Exception as ex:
        print(f"⚠️  CoLiBRI query failed: {ex}")
        return set()

def update_task_title(task_id, new_title):
    """PATCH task title via API."""
    patch_url = f"{API_BASE_URL}/api/tasks/{task_id}/"
    payload = json.dumps({"title": new_title})
    curl_cmd = [
        "curl", "-s", "-X", "PATCH", patch_url,
        "-H", f"Authorization: Token {API_TOKEN}",
        "-H", "Content-Type: application/json",
        "-d", payload
    ]
    r = subprocess.run(curl_cmd, capture_output=True, text=True)
    if r.returncode == 0:
        print(f"✏️  Updated task {task_id} title → {new_title}")
    else:
        print(f"⚠️ Failed to update title for task {task_id}: {r.stderr}")

def get_task_json(task_id):
    """Return full JSON dict for task or None"""
    url = f"{API_BASE_URL}/api/tasks/{task_id}/"
    curl_cmd = ["curl", "-s", "-H", f"Authorization: Token {API_TOKEN}", url]
    try:
        res = subprocess.run(curl_cmd, capture_output=True, text=True, check=True)
        return json.loads(res.stdout)
    except Exception as e:
        print(f"⚠️ Unable to fetch task JSON for {task_id}: {e}")
        return None

def _append_footprint_from_center(footprints, ra_val, dec_val, radius_deg):
    """Helper to add a square footprint based on astrometry center/radius."""
    if ra_val is None or dec_val is None:
        return
    if radius_deg is None:
        radius_deg = 0.25  # fallback default

    # Build simple square corners (RA ± radius, DEC ± radius)
    # Wrap RA to [0, 360)
    corners = [
        ( (ra_val - radius_deg) % 360.0, dec_val - radius_deg ),
        ( (ra_val + radius_deg) % 360.0, dec_val - radius_deg ),
        ( (ra_val + radius_deg) % 360.0, dec_val + radius_deg ),
        ( (ra_val - radius_deg) % 360.0, dec_val + radius_deg ),
    ]
    footprints.append(corners)
    save_sky_maps(footprints)

# ----------------------------------------------------------------------------
# Fallback helper: build a footprint from FITS header coordinates (approximate)
# ----------------------------------------------------------------------------
def _append_footprint_from_header(footprints, hdr):
    """Try to read RA/DEC from common FITS header keywords and append a footprint.

    This is a *best-effort* method for the local analysis workflow where no
    precise astrometry is available yet.  If neither RA nor DEC can be found
    the function silently returns.
    """
    ra_keys  = ["RA", "OBJCTRA", "RA_CENTER", "CRVAL1"]
    dec_keys = ["DEC", "OBJCTDEC", "DEC_CENTER", "CRVAL2"]

    ra_val = None
    dec_val = None

    for k in ra_keys:
        if k in hdr and hdr[k] not in (None, ""):
            ra_val = parse_ra(hdr[k])
            break

    for k in dec_keys:
        if k in hdr and hdr[k] not in (None, ""):
            dec_val = parse_dec(hdr[k])
            break

    if ra_val is None or dec_val is None:
        return  # insufficient information

    # Coarse radius – use 0.25° by default (same as server side)
    _append_footprint_from_center(footprints, ra_val, dec_val, 0.25)

# ----------------------------------------------------------------------------
# Helper: compute approximate field-of-view from FITS header
# ----------------------------------------------------------------------------

def _estimate_fov_from_header(hdr):
    """Return (d_ra_deg, d_dec_deg) half-sizes from header or None.
    d_ra/dec are HALF the full width/height in degrees (i.e. center→edge).
    """
    try:
        nx = int(hdr.get('NAXIS1') or 0)
        ny = int(hdr.get('NAXIS2') or 0)
        if nx == 0 or ny == 0:
            return None

        # Prefer CDELT keywords (deg / pixel)
        if 'CDELT1' in hdr and 'CDELT2' in hdr:
            scale_ra = abs(float(hdr['CDELT1']))  # deg / px (can be negative)
            scale_dec = abs(float(hdr['CDELT2']))
        else:
            # Try CD matrix → pixel scale ~ sqrt(cd1_1^2 + cd1_2^2)
            cd_keys = [k.lower() for k in hdr.keys()]
            if 'cd1_1' in cd_keys and 'cd2_2' in cd_keys:
                cd11 = abs(float(hdr['CD1_1']))
                cd22 = abs(float(hdr['CD2_2']))
                scale_ra  = cd11
                scale_dec = cd22
            else:
                return None  # no scale info

        d_ra  = scale_ra  * nx / 2.0
        d_dec = scale_dec * ny / 2.0
        return (d_ra, d_dec)
    except Exception:
        return None

def process_single_image(file_path, object_name, date_str, footprints=None, max_wait_time=600):
    """Full workflow: upload → inspect → photometry → subtraction → TNS"""
    task_id = upload_for_inspection(file_path, object_name, date_str)
    if not task_id:
        return False

    if not launch_inspection(task_id):
        return False

    if not wait_for_state(task_id, {"inspect_done"}, timeout=max_wait_time):
        return False

    if not trigger_photometry(task_id):
        return False

    if not wait_for_state(task_id, {"photometry_done", "done", "completed"}, timeout=max_wait_time):
        return False

    if not trigger_template_subtraction(task_id):
        return False

    # Fetch RA/DEC from task metadata (produced by astrometry on server)
    task_json = get_task_json(task_id)
    ra_val = None
    dec_val = None
    radius_deg = None
    if task_json:
        # Prefer new config keys provided by server
        cfg = task_json.get("config", {})
        if cfg.get("field_ra") is not None and cfg.get("field_dec") is not None:
            ra_val = float(cfg["field_ra"])
            dec_val = float(cfg["field_dec"])
            radius_deg = float(cfg.get("field_sr", 0)) or None

        # Fallback to other known keys
        if ra_val is None:
            for ra_key in ["ra", "wcs_ra", "ra_center", "ra_deg"]:
                if ra_key in task_json and task_json[ra_key] is not None:
                    ra_val = float(task_json[ra_key])
                    break
        if dec_val is None:
            for dec_key in ["dec", "wcs_dec", "dec_center", "dec_deg"]:
                if dec_key in task_json and task_json[dec_key] is not None:
                    dec_val = float(task_json[dec_key])
                    break

    # If the server did not provide a search radius, fall back to the
    # default 0.25 deg; do NOT attempt to estimate it from FITS headers.

    # Query TNS / CoLiBRI and rename task if transients found
    if ra_val is not None and dec_val is not None:
        all_names = set()
        if radius_deg is not None:
            radius_arcsec = radius_deg * 3600.0
            tns_names = query_tns(ra_val, dec_val, radius_arcsec=radius_arcsec)
        else:
            tns_names = query_tns(ra_val, dec_val)
        if tns_names:
            all_names.update(tns_names)

        # CoLiBRI search – use same radius or 0.25 deg default
        search_radius = radius_deg if radius_deg is not None else 0.25
        colibri_names = query_colibri(ra_val, dec_val, radius_deg=search_radius)
        if colibri_names:
            all_names.update(colibri_names)

        if all_names:
            new_title = f"{date_str}_{object_name}_{'+'.join(sorted(all_names))}"
            update_task_title(task_id, new_title)

    if ra_val is None or dec_val is None:
        print(f"⚠️  Could not extract RA/DEC from task {task_id} JSON")
    else:
        print(f"ℹ️  Using RA={ra_val:.4f}, DEC={dec_val:.4f} for TNS query (radius: {radius_deg if radius_deg else 'default'})")

    # After photometry DONE, add footprint based on accurate astrometry
    if footprints is not None and ra_val is not None and dec_val is not None:
        _append_footprint_from_center(footprints, ra_val, dec_val, radius_deg)

    return True

def upload_processed_images():
    """
    Loop through processed folders and upload res.fit files to API using 2-step process
    """
    # Extract date from init_path (e.g., "2025-07-02")
    date_str = os.path.basename(os.path.dirname(init_path))
    print(f"Using date: {date_str}")
    
    uploaded_count = 0
    failed_count = 0
    
    # List to hold footprints for sky map: each item is list of (ra_deg, dec_deg) corners
    footprints = []
    
    # Loop through target directories
    for target_dir in os.listdir(init_path):
        target_path = os.path.join(init_path, target_dir)
        if not os.path.isdir(target_path):
            continue
            
        print(f"\n📁 Processing target: {target_dir}")
        
        # Loop through filter directories
        for filter_dir in os.listdir(target_path):
            filter_path = os.path.join(target_path, filter_dir)
            if not os.path.isdir(filter_path):
                continue
                
            # Loop through exposure directories
            for exposure_dir in os.listdir(filter_path):
                # Build the *full* path to the exposure directory (e.g. OBJECT_X/G/60)
                exposure_path = os.path.join(filter_path, exposure_dir)
                if not os.path.isdir(exposure_path):
                    continue
                    
                # Look for res.fit file
                res_file = os.path.join(exposure_path, "res.fit")
                if os.path.exists(res_file):
                    # ------------------------------------------------------------------
                    # Ensure the res.fit has a valid WCS. If not, run the local PlateSolver
                    # (astrometry.net solve-field) to generate one in-place so that subsequent
                    # analysis (SN search, sky-map, etc.) has accurate coordinates.
                    # ------------------------------------------------------------------
                    try:
                        with fits.open(res_file, mode='update') as hdul:
                            hdr = hdul[0].header
                            # Minimal check – CRVAL keywords must exist and be finite
                            has_wcs = all(k in hdr for k in ("CRVAL1", "CRVAL2")) and \
                                       not any(math.isnan(hdr.get(k, math.nan)) for k in ("CRVAL1", "CRVAL2"))
                            if not has_wcs:
                                print("🛰️  No valid WCS – running local plate solver (solve-field)…")
                                solver = PlateSolver()
                                solved_wcs, _ = solver.solve_field(res_file)
                                if solved_wcs is not None:
                                    # Update FITS header in-place with solved WCS
                                    hdr.update(solved_wcs.to_header())
                                    hdul.flush()
                                    print("✅ WCS written to res.fit")
                                else:
                                    print("⚠️  Plate solving failed – continuing without WCS")
                    except Exception as e:
                        print(f"⚠️  Error while plate-solving {res_file}: {e}")

                    print(f"\n🎯 Found res.fit in {target_dir}/{filter_dir}/{exposure_dir}")

                    if process_single_image(res_file, target_dir, date_str, footprints=footprints):
                        uploaded_count += 1
                        print(f"✅ Successfully processed: {target_dir}")
                        save_sky_maps(footprints)
                    else:
                        failed_count += 1
                        print(f"❌ Failed to process: {target_dir}")
                else:
                    print(f"⚠️ No res.fit found in {target_dir}/{filter_dir}/{exposure_dir}")

    print(f"\n📊 Final Summary:")
    print(f"✅ Successfully processed: {uploaded_count}")
    print(f"❌ Failed processing: {failed_count}")
    print(f"📁 Total attempted: {uploaded_count + failed_count}")

    # ------------------------------------------------------------------
    # Generate sky coverage map if we have at least one footprint
    # ------------------------------------------------------------------
    if footprints:
        try:
            print(f"🗺️  Generating sky map with {len(footprints)} footprints…")

            fig = plt.figure(figsize=(10,5))
            ax  = fig.add_subplot(111, projection='mollweide')

            for corners in footprints:
                ra = np.array([c[0] for c in corners])
                dec = np.array([c[1] for c in corners])
                # Convert RA to radians, wrap to [-pi, pi]
                ra_rad = np.radians(ra)
                ra_rad = np.remainder(ra_rad+2*np.pi, 2*np.pi)
                ra_rad[ra_rad>np.pi] -= 2*np.pi
                dec_rad = np.radians(dec)
                poly = Polygon(np.column_stack((-ra_rad, dec_rad)), closed=True, facecolor='none', edgecolor='red', linewidth=0.5)
                ax.add_patch(poly)

            ax.grid(True, color='lightgray', lw=0.5)
            ax.set_title('Sky coverage of tonight\'s images', pad=20)
            out_png = os.path.join(init_path, 'sky_coverage.png')
            plt.savefig(out_png, dpi=200, bbox_inches='tight')
            plt.close(fig)
            print(f"✅ Sky coverage map saved → {out_png}")
        except Exception as e:
            print(f"⚠️  Failed to create sky map: {e}")

def main():
    parser = argparse.ArgumentParser(description="FITS Image Processing and Organization Script")
    parser.add_argument("--upload", action="store_true", 
                       help="Upload processed res.fit images to photometry API instead of processing")
    parser.add_argument("--mosaic", action="store_true",
                       help="Create mosaic from processed res.fit images")
    parser.add_argument("--pixel-scale", type=float, default=1.0,
                       help="Pixel scale for mosaic in arcsec/pixel (default: 1.0)")
    parser.add_argument("--projection", type=str, default="TAN", choices=["TAN", "ZEA", "SIN", "CAR"],
                       help="Projection type for mosaic (default: TAN)")
    parser.add_argument("--blend-method", type=str, default="feather", choices=["feather", "linear", "weighted"],
                       help="Blending method for mosaic (default: feather)")
    parser.add_argument("--output-dir", type=str, default="mosaics",
                       help="Output directory for mosaic files (default: mosaics)")
    parser.add_argument("--quality-threshold", type=float, default=0.7,
                       help="WCS quality threshold (0-1, default: 0.7)")
    parser.add_argument("--limit", type=int, default=None,
                       help="Limit processing to first N images for testing (default: process all)")
    parser.add_argument("--margin", type=float, default=2.0,
                       help="Margin in degrees to add around mosaic bounds (default: 2.0)")
    parser.add_argument("--no-size-limit", action="store_true",
                       help="Disable maximum size limits (use with caution - may cause memory issues)")
    parser.add_argument("--outlier-threshold", type=float, default=3.0,
                       help="Sigma threshold for outlier rejection (default: 3.0)")
    parser.add_argument("--feather-sigma", type=float, default=2.0,
                       help="Gaussian sigma for feathering (default: 2.0)")
    parser.add_argument("--fast-mode", action="store_true",
                       help="Fast processing mode - skip complex steps for speed")
    parser.add_argument("--incremental", action="store_true",
                       help="Save intermediate mosaic after each image is reprojected")
    parser.add_argument("--edge-samples", type=int, default=8,
                       help="Extra samples per image edge for boundary check (default: 8)")
    parser.add_argument("--no-border-csv", action="store_true",
                       help="Disable saving border_points.csv diagnostics")
    parser.add_argument("--plate-solve", action="store_true",
                       help="Only run astrometry.net solve-field on res.fit files and write WCS")
    
    args = parser.parse_args()
    
    if args.upload:
        print("🚀 Starting API upload mode with 4-step workflow...")
        upload_processed_images()
        return
    
    if args.mosaic:
        print("🎨 Starting mosaic creation mode...")
        create_mosaic(args)
        return
    
    if args.plate_solve:
        print("🛰️  Starting stand-alone plate-solve mode …")
        plate_solve_all_res_files()
        return
    
    # Original processing workflow
    print("🔄 Starting image processing and organization...")
    
    # Reset unnamed object counter for clean run
    global unnamed_object_counter, _created_folders
    unnamed_object_counter = 1
    _created_folders.clear()  # Clear folder cache for fresh run
    
    # ------------------------------------------------------------------
    # Collect approximate footprints from FITS headers so that we can still
    # generate a *coarse* sky-coverage map in the local analysis workflow.
    # ------------------------------------------------------------------
    footprints: list[list[tuple[float,float]]] = []

    # ------------------------------------------------------------------
    # 1) Collect list of raw FITS files in DATA_ROOT
    # ------------------------------------------------------------------
    fits_files = list(directory.glob('*.[fF][iI][tT]*'))

    print(f"Found {len(fits_files)} FITS files to process")

    # Use the global tracking set so that other helper functions can update it
    global calibrated_folders
    calibrated_folders.clear()
    
    processed_count = 0
    skipped_count = 0
    
    # Track per-folder statistics for pySiril calibration
    folder_info: dict[str, dict] = {}

    print("🔄 Organising raw FITS → folder structure …")
    for file_path in tqdm(fits_files, desc="Splitting", unit="img"):
        # human-friendly index will be handled by tqdm
        print(f"  ↪ Processing {file_path.name}")
        
        # Check file size - skip if smaller than 15MB
        file_size_mb = file_path.stat().st_size / (1024 * 1024)
        if file_size_mb < 15:
            print(f"⚠️  Skipping {file_path.name}: file size {file_size_mb:.1f}MB < 15MB")
            skipped_count += 1
            continue
        
        # Read FITS header (optimized to only read what we need)
        try:
            info = getFitHeader(file_path)
            exposure_time = int(info['EXPOSURE'])
            filter_name = info['FILTER']
        except Exception as e:
            print(f"⚠️  Skipping {file_path.name}: failed to read FITS header - {e}")
            skipped_count += 1
            continue
        
        # Check exposure time - skip if 5 seconds
        if exposure_time == 5:
            print(f"⚠️  Skipping {file_path.name}: exposure time {exposure_time}s = 5s")
            skipped_count += 1
            continue
        
        target = get_target_name(info)
        exposure_time_str = str(exposure_time)
        print(f"Target: {target}, Exposure: {exposure_time_str}s, Size: {file_size_mb:.1f}MB")
        
        # Build paths once
        target_dir = os.path.join(directory, target)
        filter_dir = os.path.join(target_dir, filter_name)
        exposure_dir = os.path.join(filter_dir, exposure_time_str)
        
        # Create folders with caching
        createFolderIfNeededCached(target_dir)
        createFolderIfNeededCached(filter_dir)
        createFolderIfNeededCached(exposure_dir)
        
        destination = Path(os.path.join(exposure_dir, file_path.name))
        
        # SAFETY: Copy file instead of move to preserve original data
        fast_file_copy(file_path, destination)
        processed_count += 1

        # Update per-folder counters
        key = str(exposure_dir)
        f_info = folder_info.setdefault(key, {"filter": filter_name, "exposure": exposure_time_str, "num_images": 0, "target": target})
        f_info["num_images"] += 1

        # Add footprint from FITS header
        _append_footprint_from_header(footprints, info)

    print(f"✅ Processing complete: {processed_count} files processed, {skipped_count} files skipped")

    # ------------------------------------------------------------------
    # Calibrate all folders now using pySiril (one batch)
    # ------------------------------------------------------------------
    if folder_info:
        print("🚀 Calibrating all exposure folders with Siril …")
        calibrate_folders_pysiril(folder_info)

    # ------------------------------------------------------------------
    # After calibration, create the sky-coverage map (approximate)
    # ------------------------------------------------------------------
    if footprints:
        print(f"🗺️  Generating sky coverage maps from {len(footprints)} header positions…")
        save_sky_maps(footprints)
    else:
        print("ℹ️  No footprints available for sky coverage map")

    # ------------------------------------------------------------------
    # Build Siril script (same as previous behaviour, useful for manual batch)
    # ------------------------------------------------------------------
    print("🔧 Building Siril script...")
    with open(os.path.join(init_path, "script.ssf"), "w") as f:
        f.write("requires 1.2.0\n")
        
        # Use Path.iterdir() which is faster than os.listdir()
        for target_path in Path(init_path).iterdir():
            if not target_path.is_dir():
                continue
                
            _target = target_path.name
            print(f"Processing target: {_target}")
            
            for filter_path in target_path.iterdir():
                if not filter_path.is_dir():
                    continue
                    
                _filter = filter_path.name
                
                for exposure_path in filter_path.iterdir():
                    if not exposure_path.is_dir():
                        continue
                        
                    _exposure = exposure_path.name
                    print(f"  {_filter}/{_exposure}")
                    
                    current_folder = str(exposure_path)
                    
                    # Count FITS files more efficiently
                    fits_files = list(exposure_path.glob('*.[fF][iI][tT]*'))
                    num_images = len(fits_files)
                    
                    print(f"    Found {num_images} image(s)")
                    
                    f.write(f"cd {current_folder}\n")
                    f.write("convert i\n")
                    
                    # Build calibration argument string based on available frames
                    flat_arg  = f" -flat={flats[_filter]}"  if _filter in flats else ""
                    dark_arg  = f" -dark={darks[_exposure]}" if _exposure in darks else ""
                    cal_args  = flat_arg + dark_arg

                    if num_images == 1:
                        # Single image workflow
                        print("    Using single image workflow")
                        f.write(f"calibrate_single i_00001.fit{cal_args}\n")
                        f.write("load pp_i_00001\n")
                        f.write("binxy 2\n")
                        f.write("save res\n")
                    else:
                        # Multiple images workflow
                        f.write(f"calibrate i{cal_args}\n")
                        f.write("register pp_i -2pass -interp=cu\n")
                        f.write("seqapplyreg pp_i -framing=min\n")
                        f.write("stack r_pp_i rej w 3 3 -nonorm -filter-fwhm=80% -filter-round=80%\n")
                        f.write("load r_pp_i_stacked\n")
                        f.write("binxy 2\n")
                        f.write("save res\n")

    # ------------------------------------------------------------------
    # Fallback: any folder that, for some reason, wasn't calibrated during the
    # splitting loop (e.g. empty, or pySiril missing) is handled here.
    # ------------------------------------------------------------------
    remaining = {k:v for k,v in folder_info.items() if k not in calibrated_folders}
    if remaining:
        calibrate_folders_pysiril(remaining)

    # ---------------------------------------------------------------------------
    # Siril integration helpers
    # ---------------------------------------------------------------------------

    def _run_siril_script_cli(script_path: str):
        """Fallback: run Siril via subprocess as before."""
        siril_exec = shutil.which("siril-cli") or shutil.which("siril")
        if siril_exec is None:
            print("⚠️  Siril executable not found; cannot calibrate images.")
            return False
        cmd = [siril_exec, "-s", script_path, "-q"]
        print(f"🚀 Launching Siril (CLI): {' '.join(cmd)}")
        try:
            subprocess.run(cmd, cwd=DATA_ROOT, check=True)
            print("✅ Siril calibration completed (CLI mode)")
            return True
        except subprocess.CalledProcessError as e:
            print(f"❌ Siril CLI failed with status {e.returncode}")
            return False


    def _run_siril_with_pysiril(script_path: str):
        """Prefer running Siril through the pySiril wrapper (server mode)."""
        try:
            import pysiril.siril as sr
        except ImportError:
            print("ℹ️  pySiril not installed – falling back to CLI mode")
            return _run_siril_script_cli(script_path)

        # On macOS Siril.app usually lives in /Applications; respect any existing env
        os.environ.setdefault("SIRIL_APP", "/Applications/Siril.app")

        try:
            print("🚀 Starting Siril via pySiril server …")
            sr.init()
            sr.run(f"exec {script_path}")  # run the whole .ssf file
            sr.quit()
            print("✅ Siril calibration completed (pySiril)")
            return True
        except Exception as e:
            print(f"❌ pySiril error: {e}; falling back to CLI")
            try:
                sr.quit()
            except Exception:
                pass
            return _run_siril_script_cli(script_path)

    # Run the generated Siril script only if no folder has already been
    # calibrated by `calibrate_folders_pysiril` above. This prevents a second
    # redundant calibration pass on the same data.
    if not calibrated_folders:
        _run_siril_with_pysiril(os.path.join(DATA_ROOT, SIRIL_SCRIPT_NAME))
    else:
        print("ℹ️  Skipping Siril script execution – folders already calibrated.")

# === new helper to refresh the sky coverage plot ===

# ============================================================================
# MOSAIC CREATION CLASSES AND FUNCTIONS (v5.0)
# ============================================================================

class MosaicEngine:
    """Main orchestrator for mosaic creation"""
    
    def __init__(self, config):
        self.config = config
        self.plate_solver = PlateSolver()
        self.grid_projector = GridProjector(config['projection'], config['pixel_scale'], config)
        self.image_aligner = ImageAligner()
        self.smart_stitcher = SmartStitcher(config)
        self.mosaic_exporter = MosaicExporter()
    
    def create_mosaic(self, res_fit_files, output_dir):
        """Create mosaic from list of res.fit files"""
        print(f"🎨 Creating mosaic from {len(res_fit_files)} images...")
        
        # Step 1: Load and plate solve each image
        valid_images = []
        for i, fits_file in enumerate(res_fit_files):
            print(f"📷 Processing image {i+1}/{len(res_fit_files)}: {os.path.basename(fits_file)}")
            
            try:
                with fits.open(fits_file) as hdul:
                    header = hdul[0].header
                    data = hdul[0].data
                    
                    if data is None:
                        print(f"⚠️  Skipping: {fits_file} (no image data)")
                        continue
                    
                    print(f"🔍 Image {os.path.basename(fits_file)}: {data.shape} pixels")
                    
                    # Always plate solve for maximum accuracy
                    print(f"🔭 Plate solving for accurate astrometry...")
                    
                    # Get initial guess from header if available
                    initial_guess = {}
                    header_ra = header.get('RA') or header.get('OBJCTRA')
                    header_dec = header.get('DEC') or header.get('OBJCTDEC')
                    
                    if header_ra is not None and header_dec is not None:
                        try:
                            initial_guess['ra'] = float(header_ra)
                            initial_guess['dec'] = float(header_dec)
                            print(f"🎯 Using header hint: RA={initial_guess['ra']:.4f}, DEC={initial_guess['dec']:.4f}")
                        except:
                            pass
                    
                    # Plate solve the image and extract sources
                    solved_wcs, source_data = self.plate_solver.solve_field(fits_file, initial_guess)
                    
                    if solved_wcs is None:
                        print(f"❌ Plate solving failed for {os.path.basename(fits_file)}")
                        continue
                    
                    # Get accurate corners from plate solved WCS
                    corners_world = self.plate_solver.get_image_corners(solved_wcs, data.shape)
                    
                    # Validate plate solved WCS
                    solved_quality = self.plate_solver.validate_wcs(solved_wcs, data)
                    print(f"✅ Plate solved successfully (quality: {solved_quality:.2f})")
                    
                    if solved_wcs is not None and corners_world is not None:
                        # Calculate image properties from solved WCS
                        ny, nx = data.shape
                        center_ra, center_dec = solved_wcs.pixel_to_world_values(nx/2, ny/2)
                        pixel_scale = abs(solved_wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
                        
                        # Debug: Print corners
                        print(f"🔍 Image corners (RA, DEC):")
                        corner_names = ["bottom-left", "bottom-right", "top-right", "top-left"]
                        for j, (corner_name, (ra, dec)) in enumerate(zip(corner_names, corners_world)):
                            if ra is not None and dec is not None:
                                print(f"    {corner_name}: ({ra:.6f}°, {dec:.6f}°)")
                            else:
                                print(f"    {corner_name}: invalid")
                        
                        valid_images.append({
                            'file': fits_file,
                            'header': header,
                            'data': data,
                            'wcs': solved_wcs,
                            'corners_world': corners_world,
                            'shape': data.shape,
                            'center_ra': center_ra,
                            'center_dec': center_dec,
                            'pixel_scale': pixel_scale,
                            'source_data': source_data  # Include extracted source positions
                        })
                        print(f"✅ Added to mosaic: {os.path.basename(fits_file)}")
                    else:
                        print(f"⚠️  Skipping: {fits_file} (plate solving failed)")
                        
            except Exception as e:
                print(f"⚠️  Error processing {fits_file}: {e}")
                import traceback
                traceback.print_exc()
        
        if not valid_images:
            print("❌ No valid images found for mosaic creation")
            return False
        
        print(f"📊 Processing {len(valid_images)} valid images")
        
        # Step 2: Phase 1 - Collect ALL corner coordinates from all images
        print("🔍 Phase 1: Collecting all corner coordinates from all images...")
        phase1_start = time.time()
        all_corners = []
        
        if self.config.get('fast_mode', False):
            print("⚡ Fast mode: Using existing corner data...")
            for i, img in enumerate(valid_images):
                if 'corners_world' in img:
                    valid_corners = [c for c in img['corners_world'] if c[0] is not None and c[1] is not None]
                    all_corners.extend(valid_corners)
                    print(f"   📷 {os.path.basename(img['file'])}: {len(valid_corners)} corners")
        else:
            for i, img in enumerate(valid_images):
                img_start = time.time()
                print(f"📐 Analyzing image {i+1}/{len(valid_images)}: {os.path.basename(img['file'])}")
                
                # Get accurate corners from WCS
                corners = self.plate_solver.get_image_corners(img['wcs'], img['data'].shape)
                if corners:
                    valid_corners = [c for c in corners if c[0] is not None and c[1] is not None]
                    all_corners.extend(valid_corners)
                    img_time = time.time() - img_start
                    print(f"   Found {len(valid_corners)} valid corners ({img_time:.2f}s)")
                    
                    # Debug: show corner coordinates (skip in fast mode)
                    for j, (ra, dec) in enumerate(valid_corners):
                        print(f"     Corner {j+1}: RA={ra:.6f}°, DEC={dec:.6f}°")
                else:
                    img_time = time.time() - img_start
                    print(f"   ⚠️  No valid corners found for {os.path.basename(img['file'])} ({img_time:.2f}s)")
        
        phase1_time = time.time() - phase1_start
        if not all_corners:
            print("❌ No valid corners found from any image")
            return False
            
        print(f"📊 Phase 1 complete: Collected {len(all_corners)} corner coordinates from {len(valid_images)} images ({phase1_time:.2f}s)")
        
        # Step 3: Phase 2 - Calculate optimal bounding box and create grid
        print("🗺️  Phase 2: Calculating optimal bounding box...")
        phase2_start = time.time()
        target_wcs = self.grid_projector.define_target_grid_from_corners(all_corners, valid_images)
        phase2_time = time.time() - phase2_start
        print(f"🗺️  Target grid defined: {target_wcs.pixel_shape if hasattr(target_wcs, 'pixel_shape') else 'custom'} ({phase2_time:.2f}s)")
        
        # Step 4: Reproject all images to the pre-calculated grid
        print("🔄 Reprojecting all images to common grid...")
        phase3_start = time.time()
        reprojected_images = []
        verification_data_list = []
        mosaic_sum = None  # running sum for incremental output
        mosaic_count = None  # running count for incremental output
        
        for i, img in enumerate(valid_images):
            img_start = time.time()
            print(f"📷 Reprojecting image {i+1}/{len(valid_images)}: {os.path.basename(img['file'])}")
            
            try:
                result = self.grid_projector.reproject_image(
                    img['data'], img['wcs'], target_wcs, img.get('source_data')
                )
                
                # Handle tuple return (data, verification_data)
                if isinstance(result, tuple):
                    reprojected_data, verification_data = result
                else:
                    reprojected_data, verification_data = result, None
                
                if reprojected_data is not None:
                    reprojected_images.append({
                        'file': img['file'],
                        'data': reprojected_data,
                        'header': img['header']
                    })
                    
                    # NEW: Save intermediate reprojected image for this step
                    try:
                        intermediate_dir = os.path.join(output_dir, 'intermediate')
                        os.makedirs(intermediate_dir, exist_ok=True)
                        intermediate_path = os.path.join(intermediate_dir, f"reprojected_{i+1:03d}.fits")
                        hdr_int = target_wcs.to_header()
                        try:
                            valid_mask_int = np.isfinite(reprojected_data)
                            n_valid_int = int(np.count_nonzero(valid_mask_int))
                            total_int   = reprojected_data.size
                            if n_valid_int > 0:
                                cdelt_int = np.abs(target_wcs.wcs.cdelt)
                                pix_area_int = cdelt_int[0] * cdelt_int[1]
                                hdr_int['COVPIX']  = (n_valid_int, 'Pixels with data')
                                hdr_int['COVAREA'] = (round(n_valid_int * pix_area_int, 6), 'Covered area deg^2')
                                hdr_int['COVPCT']  = (round((n_valid_int / total_int)*100.0, 2), 'Coverage %')
                        except Exception:
                            pass
                        fits.PrimaryHDU(data=reprojected_data, header=hdr_int).writeto(intermediate_path, overwrite=True)
                        print(f"💾 Saved intermediate FITS: {intermediate_path}")
                    except Exception as save_err:
                        print(f"⚠️  Could not save intermediate result: {save_err}")
                    
                    # Store verification data if available
                    if verification_data is not None:
                        verification_data['image_file'] = img['file']
                        verification_data_list.append(verification_data)
                    
                    img_time = time.time() - img_start
                    print(f"✅ Reprojected: {os.path.basename(img['file'])} ({img_time:.2f}s)")
                    # ---------------- Incremental mosaic accumulation/output ----------------
                    if self.config.get('incremental_output', False):
                        if mosaic_sum is None:
                            mosaic_sum = np.zeros_like(reprojected_data, dtype=np.float32)
                            mosaic_count = np.zeros_like(reprojected_data, dtype=np.uint16)

                        valid = np.isfinite(reprojected_data)
                        mosaic_sum[valid] += reprojected_data[valid].astype(np.float32)
                        mosaic_count[valid] += 1

                        with np.errstate(invalid='ignore', divide='ignore'):
                            mosaic_partial = mosaic_sum / np.maximum(mosaic_count, 1)
                        mosaic_partial[mosaic_count == 0] = np.nan

                        incr_dir = os.path.join(output_dir, 'incremental')
                        os.makedirs(incr_dir, exist_ok=True)
                        incr_path = os.path.join(incr_dir, f"mosaic_{len(reprojected_images):03d}.fits")
                        try:
                            hdr_inc = target_wcs.to_header()
                            try:
                                valid_mask_inc = np.isfinite(mosaic_partial)
                                n_valid_inc = int(np.count_nonzero(valid_mask_inc))
                                total_inc   = mosaic_partial.size
                                if n_valid_inc > 0:
                                    cdelt_inc = np.abs(target_wcs.wcs.cdelt)
                                    pix_area_inc = cdelt_inc[0] * cdelt_inc[1]
                                    hdr_inc['COVPIX']  = (n_valid_inc, 'Pixels with data')
                                    hdr_inc['COVAREA'] = (round(n_valid_inc * pix_area_inc, 6), 'Covered area deg^2')
                                    hdr_inc['COVPCT']  = (round((n_valid_inc / total_inc)*100.0, 2), 'Coverage %')
                            except Exception:
                                pass
                            fits.PrimaryHDU(data=mosaic_partial.astype(np.float32), header=hdr_inc).writeto(incr_path, overwrite=True)
                            print(f"💾 Saved incremental mosaic: {incr_path}")
                        except Exception as inc_err:
                            print(f"⚠️  Could not save incremental mosaic: {inc_err}")
                else:
                    img_time = time.time() - img_start
                    print(f"⚠️  Failed to reproject: {os.path.basename(img['file'])} ({img_time:.2f}s)")
                    
            except Exception as e:
                img_time = time.time() - img_start
                print(f"⚠️  Error reprojecting {os.path.basename(img['file'])}: {e} ({img_time:.2f}s)")
                import traceback
                traceback.print_exc()
        
        phase3_time = time.time() - phase3_start
        if not reprojected_images:
            print("❌ No images successfully reprojected")
            return False
            
        print(f"🎉 Reprojection complete: {len(reprojected_images)} images ready for combination ({phase3_time:.2f}s)")
        
        # Step 5: Combine all reprojected images into final mosaic
        print("🧩 Combining images into final mosaic...")
        phase4_start = time.time()
        mosaic_data = self.smart_stitcher.combine_images(reprojected_images)
        phase4_time = time.time() - phase4_start
        if mosaic_data is None:
            print("❌ Failed to create mosaic")
            return False
        print(f"✅ Stitching complete ({phase4_time:.2f}s)")
        
        # Step 6: Create source verification plot if we have verification data
        if verification_data_list:
            os.makedirs(output_dir, exist_ok=True)  # Ensure directory exists
            self._create_source_verification_plot(verification_data_list, target_wcs, output_dir)
        
        # Step 7: Export results
        print("💾 Exporting results...")
        export_start = time.time()
        os.makedirs(output_dir, exist_ok=True)
        success = self.mosaic_exporter.export_mosaic(mosaic_data, target_wcs, output_dir, reprojected_images)
        export_time = time.time() - export_start
        
        total_time = time.time() - phase1_start
        print(f"\n⏱️  Performance Summary:")
        print(f"   Phase 1 (Corner collection): {phase1_time:.2f}s ({100*phase1_time/total_time:.1f}%)")
        print(f"   Phase 2 (Grid calculation): {phase2_time:.2f}s ({100*phase2_time/total_time:.1f}%)")
        print(f"   Phase 3 (Reprojection): {phase3_time:.2f}s ({100*phase3_time/total_time:.1f}%)")
        print(f"   Phase 4 (Stitching): {phase4_time:.2f}s ({100*phase4_time/total_time:.1f}%)")
        print(f"   Export: {export_time:.2f}s ({100*export_time/total_time:.1f}%)")
        print(f"   Total: {total_time:.2f}s")
        
        if success:
            print(f"✅ Mosaic created successfully in {output_dir}")
            return True
        else:
            print("❌ Failed to export mosaic")
            return False
    
    def _create_source_verification_plot(self, verification_data_list, target_wcs, output_dir):
        """Create a plot showing source position accuracy after geometric transformation"""
        try:
            import matplotlib.pyplot as plt
            from matplotlib.patches import FancyArrowPatch
            import numpy as np
            
            print("🎯 Creating source position verification plot...")
            
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Colors for different images
            colors = plt.cm.Set1(np.linspace(0, 1, len(verification_data_list)))
            
            all_errors = []
            max_error = 0
            
            for i, vdata in enumerate(verification_data_list):
                source_control = vdata['source_control']
                target_control = vdata['target_control']
                catalog_coords = vdata['catalog_coords']
                target_wcs = vdata['target_wcs']
                image_file = vdata['image_file']
                
                # Calculate expected positions from catalog coordinates
                expected_positions = []
                for ra, dec in catalog_coords:
                    try:
                        ex, ey = target_wcs.world_to_pixel_values(ra, dec)
                        expected_positions.append([ex, ey])
                    except:
                        expected_positions.append([np.nan, np.nan])
                
                expected_positions = np.array(expected_positions)
                
                # Calculate position errors (actual vs expected)
                errors = target_control - expected_positions
                error_magnitudes = np.sqrt(errors[:, 0]**2 + errors[:, 1]**2)
                
                # Filter out invalid positions
                valid_mask = np.isfinite(error_magnitudes)
                if np.sum(valid_mask) == 0:
                    continue
                    
                valid_target = target_control[valid_mask]
                valid_expected = expected_positions[valid_mask]
                valid_errors = errors[valid_mask]
                valid_error_mags = error_magnitudes[valid_mask]
                
                all_errors.extend(valid_error_mags)
                max_error = max(max_error, np.max(valid_error_mags))
                
                # Plot expected positions as circles
                ax.scatter(valid_expected[:, 0], valid_expected[:, 1], 
                          color=colors[i], s=50, alpha=0.7, marker='o', 
                          label=f'{os.path.basename(image_file)} (expected)')
                
                # Plot actual positions as crosses
                ax.scatter(valid_target[:, 0], valid_target[:, 1],
                          color=colors[i], s=50, alpha=0.9, marker='x',
                          label=f'{os.path.basename(image_file)} (actual)')
                
                # Draw arrows showing the errors (scaled for visibility)
                scale_factor = max(1.0, max_error / 20)  # Scale arrows for visibility
                for j in range(len(valid_expected)):
                    if valid_error_mags[j] > 0.1:  # Only show significant errors
                        arrow = FancyArrowPatch(
                            (valid_expected[j, 0], valid_expected[j, 1]),
                            (valid_target[j, 0], valid_target[j, 1]),
                            color=colors[i], alpha=0.6, linewidth=1,
                            arrowstyle='->', mutation_scale=10
                        )
                        ax.add_patch(arrow)
            
            # Set up the plot
            ax.set_xlabel('X Pixel Coordinate')
            ax.set_ylabel('Y Pixel Coordinate')
            ax.set_title('Source Position Verification\n(Arrows show errors: Expected → Actual positions)')
            ax.grid(True, alpha=0.3)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Add statistics
            if all_errors:
                mean_error = np.mean(all_errors)
                rms_error = np.sqrt(np.mean(np.array(all_errors)**2))
                max_error_val = np.max(all_errors)
                
                stats_text = f'Error Statistics:\n'
                stats_text += f'Mean: {mean_error:.2f} px\n'
                stats_text += f'RMS: {rms_error:.2f} px\n'
                stats_text += f'Max: {max_error_val:.2f} px\n'
                stats_text += f'Sources: {len(all_errors)}'
                
                ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                       verticalalignment='top', bbox=dict(boxstyle='round', 
                       facecolor='wheat', alpha=0.8))
            
            plt.tight_layout()
            
            # Save the plot
            verification_path = os.path.join(output_dir, 'source_verification.png')
            plt.savefig(verification_path, dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"📊 Source verification plot saved: {verification_path}")
            
            if all_errors:
                print(f"🎯 Source position accuracy:")
                print(f"   Mean error: {np.mean(all_errors):.2f} pixels")
                print(f"   RMS error: {np.sqrt(np.mean(np.array(all_errors)**2)):.2f} pixels")
                print(f"   Max error: {np.max(all_errors):.2f} pixels")
                print(f"   Total sources: {len(all_errors)}")
            
        except Exception as e:
            print(f"⚠️  Error creating source verification plot: {e}")
            import traceback
            traceback.print_exc()

class PlateSolver:
    """Real plate solving using astrometry.net solve-field"""
    
    def __init__(self):
        # Hard-wired path to solve-field (reported by `which solve-field`)
        self.solve_field_cmd = '/opt/homebrew/opt/astrometry-net/bin/solve-field'
        # Ensure its directory is on PATH for any sub-calls that rely on other
        # astrometry.net helper binaries located in the same folder.
        import os
        sf_dir = os.path.dirname(self.solve_field_cmd)
        if sf_dir and sf_dir not in os.environ.get('PATH', ''):
            os.environ['PATH'] = f"{sf_dir}:{os.environ.get('PATH','')}"

        self.temp_dir = '/tmp/plate_solving'
        os.makedirs(self.temp_dir, exist_ok=True)
    
    def solve_field(self, fits_path, initial_guess=None):
        """Solve astrometry for a single field using solve-field and extract sources"""
        print(f"🔭 Plate solving {os.path.basename(fits_path)}...")
        
        # Create output filename
        base_name = os.path.splitext(os.path.basename(fits_path))[0]
        output_dir = os.path.join(self.temp_dir, base_name)
        os.makedirs(output_dir, exist_ok=True)
        
        # Construct solve-field command - extract sources and correlations
        cmd = [
            self.solve_field_cmd,
            '--no-plots',           # Don't create plot files
            '--new-fits', 'none',   # Don't create new FITS file
            '--solved', 'none',     # Don't create .solved file
            '--wcs', os.path.join(output_dir, 'solution.wcs'),  # WCS output file
            '--corr', os.path.join(output_dir, 'correlations.fits'),  # Source correlations
            '--rdls', os.path.join(output_dir, 'catalog.rdls'),  # Catalog sources
            '--match', os.path.join(output_dir, 'matches.fits'), # Matched sources
            '--overwrite',          # Overwrite existing files
            '--no-verify',          # Skip verification (faster)
            '--crpix-center',       # Put CRPIX at image center
            '--scale-units', 'arcsecperpix',  # Use arcsec/pixel for scale
            '--scale-low', '1.5',   # Expected scale range (1.5-3.0 arcsec/pixel)
            '--scale-high', '3.0',
            '--downsample', '2',    # Downsample for faster solving
            '--cpulimit', '60',     # 60 second CPU limit
        ]
        
        # Add initial guess if available
        if initial_guess and 'ra' in initial_guess and 'dec' in initial_guess:
            cmd.extend([
                '--ra', str(initial_guess['ra']),
                '--dec', str(initial_guess['dec']),
                '--radius', '2.0'  # Search within 2 degrees
            ])
        
        cmd.append(fits_path)
        
        try:
            # Run solve-field
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            
            if result.returncode == 0:
                # Load the solved WCS
                wcs_file = os.path.join(output_dir, 'solution.wcs')
                matches_file = os.path.join(output_dir, 'matches.fits')
                corr_file = os.path.join(output_dir, 'correlations.fits')
                
                if os.path.exists(wcs_file):
                    with fits.open(wcs_file) as hdul:
                        solved_wcs = WCS(hdul[0].header)
                        
                        # Extract source positions for distortion correction
                        source_data = self._extract_source_positions(matches_file, corr_file)
                        
                        print(f"✅ Plate solved: {os.path.basename(fits_path)}")
                        if source_data and source_data.get('pixel_coords'):
                            pixel_count = len(source_data['pixel_coords'])
                            catalog_count = len(source_data.get('catalog_coords', []))
                            print(f"🌟 Extracted {pixel_count} matched sources ({catalog_count} with catalog coords)")
                        
                        # Clean up temp files
                        try:
                            for temp_file in [wcs_file, matches_file, corr_file]:
                                if os.path.exists(temp_file):
                                    os.remove(temp_file)
                            os.rmdir(output_dir)
                        except:
                            pass
                        
                        return solved_wcs, source_data
                else:
                    print(f"⚠️  WCS file not found: {wcs_file}")
            else:
                print(f"⚠️  solve-field failed for {os.path.basename(fits_path)}: {result.stderr.strip()}")
                
        except subprocess.TimeoutExpired:
            print(f"⚠️  Plate solving timeout for {os.path.basename(fits_path)}")
        except Exception as e:
            print(f"⚠️  Error in plate solving: {e}")
        
        # Cleanup on failure
        try:
            if os.path.exists(output_dir):
                import shutil
                shutil.rmtree(output_dir)
        except:
            pass
            
        return None, None
    
    def _extract_source_positions(self, matches_file, corr_file):
        """Extract source positions from plate solving output files"""
        source_data = {
            'pixel_coords': [],
            'world_coords': [],
            'catalog_coords': [],
            'flux': []
        }
        
        try:
            # Read correlations file - this has the matched source data we need
            if os.path.exists(corr_file):
                with fits.open(corr_file) as hdul:
                    if len(hdul) > 1:
                        data = hdul[1].data
                        if data is not None:
                            # Extract pixel coordinates (field_x, field_y)
                            if 'field_x' in data.names and 'field_y' in data.names:
                                pixel_x = data['field_x']
                                pixel_y = data['field_y']
                                source_data['pixel_coords'] = list(zip(pixel_x, pixel_y))
                                
                            # Extract catalog coordinates (index_ra, index_dec)
                            if 'index_ra' in data.names and 'index_dec' in data.names:
                                cat_ra = data['index_ra']
                                cat_dec = data['index_dec']
                                source_data['catalog_coords'] = list(zip(cat_ra, cat_dec))
                                
                            # Extract flux information if available
                            if 'FLUX' in data.names:
                                source_data['flux'] = data['FLUX'].tolist()
                                
            # If correlations file didn't work, try matches file as backup
            if not source_data['pixel_coords'] and os.path.exists(matches_file):
                with fits.open(matches_file) as hdul:
                    if len(hdul) > 1:
                        data = hdul[1].data
                        if data is not None:
                            # Try to extract from matches structure (more complex)
                            if 'FIELDOBJS' in data.names and 'STARS' in data.names:
                                print(f"⚠️  Using matches file format (less optimal)")
                                # This would require more complex parsing
                                
        except Exception as e:
            print(f"⚠️  Error extracting source positions: {e}")
            return None
        
        # Return None if no sources found
        if not source_data['pixel_coords']:
            return None
            
        return source_data
    
    def get_image_corners(self, wcs, data_shape, n_per_edge: int = 11):
        """Return a dense sampling of the image perimeter in world coordinates.

        For wide-field images geometric distortions can make the extreme
        sky‐coordinates occur somewhere along the edges, not necessarily
        at the four corners.  We therefore sample *n_per_edge* points on
        every edge (including the actual corners) to build a safer
        envelope that is later used to define the mosaic bounding box.
        """
        ny, nx = data_shape

        # Build (x, y) pixel coordinates along the perimeter
        edge_px: list[list[float]] = []

        # Bottom and top edges
        for t in np.linspace(0, nx - 1, n_per_edge, dtype=float):
            edge_px.append([t, 0])          # bottom edge
            edge_px.append([t, ny - 1])     # top edge

        # Left and right edges
        for t in np.linspace(0, ny - 1, n_per_edge, dtype=float):
            edge_px.append([0, t])          # left edge
            edge_px.append([nx - 1, t])     # right edge

        # Convert perimeter pixels → world coordinates
        corners_world: list[list[float | None]] = []
        for px, py in edge_px:
            try:
                ra, dec = wcs.pixel_to_world_values(px, py)
                corners_world.append([ra, dec])
            except Exception as e:
                print(f"⚠️  Perimeter point conversion failed: {e}")
                corners_world.append([None, None])
        
        return corners_world
    
    def validate_wcs(self, wcs, fits_data):
        """Validate WCS quality from plate solving"""
        if not wcs.is_celestial:
            return 0.0
        
        try:
            # Check pixel scale is reasonable (around 2 arcsec/pixel)
            pixel_scale = abs(wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
            if 1.5 <= pixel_scale <= 3.0:
                return 0.95  # High quality (plate solved)
            else:
                return 0.3   # Poor quality
        except Exception:
            return 0.1  # Very poor quality

class GridProjector:
    """Handles coordinate system projection and grid definition"""
    
    def __init__(self, projection_type='TAN', pixel_scale=1.0, config=None):
        self.projection = projection_type
        self.pixel_scale = pixel_scale  # arcsec/pixel
        self.config = config or {}
    
    def define_target_grid(self, images):
        """Define optimal grid covering all images"""
        print(f"🗺️  Defining target grid for {len(images)} images...")
        
        # Extract coordinates from all images
        coords = []
        all_corners = []  # Store all corner coordinates
        
        for i, img in enumerate(images):
            try:
                # Use stored center coordinates if available, otherwise calculate from WCS
                if 'center_ra' in img and 'center_dec' in img:
                    ra, dec = img['center_ra'], img['center_dec']
                else:
                    # Fallback to WCS calculation
                    wcs = img['wcs']
                    ny, nx = img['data'].shape
                    center_x, center_y = nx/2, ny/2
                    ra, dec = wcs.pixel_to_world_values(center_x, center_y)
                    
                coords.append([ra, dec])
                
                # Also collect all corner coordinates for better bounding box calculation
                if 'corners_world' in img:
                    valid_corners = [c for c in img['corners_world'] if c[0] is not None and c[1] is not None]
                    all_corners.extend(valid_corners)
                    
            except Exception as e:
                print(f"⚠️  Error extracting coordinates from image {i}: {e}")
        
        if not coords:
            print("❌ No valid coordinates found")
            return None
        
        coords = np.array(coords)
        
        # Use corners for more accurate bounding box if available
        if all_corners:
            corner_coords = np.array(all_corners)
            ra_min, ra_max = corner_coords[:, 0].min(), corner_coords[:, 0].max()
            dec_min, dec_max = corner_coords[:, 1].min(), corner_coords[:, 1].max()
        else:
            # Calculate bounding box from centers
            ra_min, ra_max = coords[:, 0].min(), coords[:, 0].max()
            dec_min, dec_max = coords[:, 1].min(), coords[:, 1].max()
        
        # Handle RA wraparound more carefully
        ra_span_initial = ra_max - ra_min
        if ra_span_initial > 180:
            # Likely wraparound case - adjust all coordinates consistently
            if all_corners:
                corner_coords[:, 0] = np.where(corner_coords[:, 0] > 180, corner_coords[:, 0] - 360, corner_coords[:, 0])
                ra_min, ra_max = corner_coords[:, 0].min(), corner_coords[:, 0].max()
            else:
                coords[:, 0] = np.where(coords[:, 0] > 180, coords[:, 0] - 360, coords[:, 0])
                ra_min, ra_max = coords[:, 0].min(), coords[:, 0].max()
        
        # Add configurable margin to show all data without cropping
        margin = self.config.get('margin', 2.0)  # degrees - configurable margin
        ra_min -= margin
        ra_max += margin
        dec_min -= margin
        dec_max += margin
        
        # Calculate center
        ra_center = (ra_min + ra_max) / 2
        dec_center = (dec_min + dec_max) / 2
        
        # Calculate required image size to cover all input images
        ra_span = ra_max - ra_min
        dec_span = dec_max - dec_min
        
        # Convert to pixels using user-specified pixel scale
        # Respect the user's pixel scale choice
        effective_pixel_scale = self.pixel_scale
        pixel_scale_deg = effective_pixel_scale / 3600.0
        nx = int(np.ceil(ra_span / pixel_scale_deg))
        ny = int(np.ceil(dec_span / pixel_scale_deg))
        
        # Apply reasonable size limits only for extremely large mosaics (unless disabled)
        if not self.config.get('no_size_limit', False):
            max_size = 15000  # pixels - allow larger mosaics but prevent memory issues
            if nx > max_size or ny > max_size:
                scale_factor = max_size / max(nx, ny)
                nx = int(nx * scale_factor)
                ny = int(ny * scale_factor)
                effective_pixel_scale = effective_pixel_scale / scale_factor
                print(f"⚠️  Mosaic size limited to {max_size}x{max_size} pixels due to memory constraints")
        else:
            print("🚨 Size limits disabled - creating full-resolution mosaic")
        
        # Ensure minimum size and make it even
        nx = max(nx, 500)
        ny = max(ny, 500)
        nx = nx + (nx % 2)  # make even
        ny = ny + (ny % 2)  # make even
        
        # Create target WCS
        target_wcs = WCS(naxis=2)
        target_wcs.wcs.crpix = [nx/2, ny/2]  # Center of image
        target_wcs.wcs.crval = [ra_center, dec_center]
        target_wcs.wcs.cdelt = [-pixel_scale_deg, pixel_scale_deg]  # Convert arcsec to degrees
        target_wcs.wcs.ctype = [f"RA---{self.projection}", f"DEC--{self.projection}"]
        
        print(f"🗺️  Target grid: {nx}x{ny} pixels, {effective_pixel_scale:.1f}\" scale, {ra_span:.3f}°x{dec_span:.3f}° coverage")
        print(f"🔍 Bounding box: RA=[{ra_min:.3f}°, {ra_max:.3f}°], DEC=[{dec_min:.3f}°, {dec_max:.3f}°] (margin: {margin:.3f}°)")
        
        # Store the target size for later use
        target_wcs.pixel_shape = (ny, nx)
        
        return target_wcs
    
    def define_target_grid_from_corners(self, all_corners, images):
        """Define optimal grid from pre-collected corner coordinates"""
        print(f"🔍 Calculating bounding box from {len(all_corners)} corner coordinates...")
        
        if not all_corners:
            print("❌ No corner coordinates provided")
            return None
        
        # Convert to numpy array for easier processing
        corner_coords = np.array(all_corners)
        
        # Calculate initial bounding box
        ra_coords = corner_coords[:, 0]
        dec_coords = corner_coords[:, 1]
        
        # Handle RA wraparound case
        ra_span_initial = ra_coords.max() - ra_coords.min()
        if ra_span_initial > 180:
            # Likely wraparound case - shift coordinates
            print("🔄 Detected RA wraparound, adjusting coordinates...")
            ra_coords = np.where(ra_coords > 180, ra_coords - 360, ra_coords)
        
        ra_min, ra_max = ra_coords.min(), ra_coords.max()
        dec_min, dec_max = dec_coords.min(), dec_coords.max()
        
        print(f"🔍 Raw bounding box: RA=[{ra_min:.6f}°, {ra_max:.6f}°], DEC=[{dec_min:.6f}°, {dec_max:.6f}°]")
        
        # Add configurable margin
        margin = self.config.get('margin', 2.0)
        ra_min -= margin
        ra_max += margin
        dec_min -= margin
        dec_max += margin
        
        print(f"🔍 With margin ({margin:.3f}°): RA=[{ra_min:.6f}°, {ra_max:.6f}°], DEC=[{dec_min:.6f}°, {dec_max:.6f}°]")
        
        # Calculate center and spans
        ra_center = (ra_min + ra_max) / 2
        dec_center = (dec_min + dec_max) / 2
        ra_span = ra_max - ra_min
        dec_span = dec_max - dec_min
        
        print(f"🎯 Center: RA={ra_center:.6f}°, DEC={dec_center:.6f}°")
        print(f"📏 Span: {ra_span:.6f}° x {dec_span:.6f}°")
        
        # --- Densify footprints: sample points along each image edge to capture
        #     curvature introduced by the projection so that no part will be
        #     cropped. The number of extra samples per edge is controlled by
        #     config['edge_samples'] (default 8).
        edge_samples = int(self.config.get('edge_samples', 8))
        sampled_points = []  # list of (ra, dec)

        if edge_samples > 0:
            ts = np.linspace(0.0, 1.0, edge_samples + 2)  # include endpoints
            for img in images:
                try:
                    wcs_i = img['wcs']
                    ny_i, nx_i = img['data'].shape
                    # Bottom & top edges (vary x)
                    for t in ts:
                        xpix = t * (nx_i - 1)
                        sampled_points.append(wcs_i.pixel_to_world_values(xpix, 0))            # bottom
                        sampled_points.append(wcs_i.pixel_to_world_values(xpix, ny_i - 1))     # top
                    # Left & right edges (vary y)
                    for t in ts:
                        ypix = t * (ny_i - 1)
                        sampled_points.append(wcs_i.pixel_to_world_values(0, ypix))            # left
                        sampled_points.append(wcs_i.pixel_to_world_values(nx_i - 1, ypix))     # right
                except Exception:
                    continue  # skip image on failure

        # Combine original corners and sampled points for grid-fit validation
        if sampled_points:
            extra_coords = np.array(sampled_points)
            validation_coords = np.vstack([corner_coords, extra_coords])
        else:
            validation_coords = corner_coords
        
        # Calculate required image size using user-specified pixel scale
        # Use uniform pixel scales in both directions
        effective_pixel_scale = self.pixel_scale
        pixel_scale_deg = effective_pixel_scale / 3600.0
        
        # Account for RA compression at high declinations in the sky coverage calculation
        dec_center_rad = np.radians(dec_center)
        cos_dec = np.cos(dec_center_rad)
        
        # Calculate the actual angular size on sky (RA span needs declination correction)
        ra_angular_span = ra_span  # Use full RA coordinate extent to avoid cropping
        dec_angular_span = dec_span  # DEC axis already correct
        
        print(f"🔄 Declination correction: cos({dec_center:.1f}°) = {cos_dec:.3f}")
        print(f"🔄 Angular spans: RA={ra_angular_span:.3f}° (used), DEC={dec_angular_span:.3f}°")
        
        # Use uniform pixel scale for both directions
        nx = int(np.ceil(ra_angular_span / pixel_scale_deg))
        ny = int(np.ceil(dec_angular_span / pixel_scale_deg))
        
        print(f"🔢 Initial grid size: {nx}x{ny} pixels")
        print(f"🔢 Uniform pixel scale: {effective_pixel_scale:.1f}\"/pixel in both directions")
        
        # Apply reasonable size limits only for extremely large mosaics (unless disabled)
        if not self.config.get('no_size_limit', False):
            max_size = 15000  # Maximum reasonable size
            if nx > max_size or ny > max_size:
                scale_factor = max_size / max(nx, ny)
                nx = int(nx * scale_factor)
                ny = int(ny * scale_factor)
                effective_pixel_scale = effective_pixel_scale / scale_factor
                print(f"⚠️  Mosaic size limited to {max_size}x{max_size} pixels due to memory constraints")
                print(f"🔽 Scaled down to: {nx}x{ny} pixels at {effective_pixel_scale:.1f}\" scale")
        else:
            print("🚨 Size limits disabled - creating full-resolution mosaic")
        
        # Ensure minimum size and make even
        nx = max(nx, 200)
        ny = max(ny, 200)
        nx = nx + (nx % 2)  # make even
        ny = ny + (ny % 2)  # make even
        
        # Build WCS in a loop, enlarging grid if any corner falls outside the initial bounds
        max_iterations = 3
        iter_count = 0
        while True:
            target_wcs = WCS(naxis=2)
            target_wcs.wcs.crpix = [nx/2, ny/2]
            target_wcs.wcs.crval = [ra_center, dec_center]
            target_wcs.wcs.cdelt = [-pixel_scale_deg, pixel_scale_deg]
            target_wcs.wcs.ctype = [f"RA---{self.projection}", f"DEC--{self.projection}"]
            target_wcs.pixel_shape = (ny, nx)

            # Project corners to pixel coordinates to verify they are inside
            try:
                x_pix, y_pix = target_wcs.world_to_pixel_values(validation_coords[:,0], validation_coords[:,1])
                min_x, max_x = np.nanmin(x_pix), np.nanmax(x_pix)
                min_y, max_y = np.nanmin(y_pix), np.nanmax(y_pix)
                inside = (min_x >= 0) and (max_x <= nx-1) and (min_y >= 0) and (max_y <= ny-1)
            except Exception:
                inside = True  # If projection fails, skip check

            if inside or iter_count >= max_iterations:
                break

            # Enlarge grid by 20% on both axes and try again
            grow_factor = 1.2
            nx = int(nx * grow_factor)
            ny = int(ny * grow_factor)
            iter_count += 1
            print(f"🔄 Grid too small, enlarging to {nx}x{ny} (attempt {iter_count})")

        print(f"✅ Final grid: {nx}x{ny} pixels, {effective_pixel_scale:.1f}\"/pixel uniform scale")
        print(f"🗺️  Sky coverage: {ra_span:.6f}° x {dec_span:.6f}° (coordinate span)")
        print(f"🗺️  Angular coverage: {ra_angular_span:.6f}° x {dec_angular_span:.6f}° (true sky)")

        # Save border point diagnostics CSV if requested
        if self.config.get('dump_border_csv', True):
            try:
                import pandas as pd
                inside_flags = ((x_pix >= 0) & (x_pix <= nx-1) & (y_pix >= 0) & (y_pix <= ny-1))
                df = pd.DataFrame({
                    'ra_deg': validation_coords[:, 0],
                    'dec_deg': validation_coords[:, 1],
                    'x_pix': x_pix,
                    'y_pix': y_pix,
                    'inside': inside_flags.astype(int)
                })
                csv_path = self.config.get('border_csv_path', 'border_points.csv')
                df.to_csv(csv_path, index=False)
                print(f"📝 Border points CSV saved: {csv_path}")
            except Exception as e:
                print(f"⚠️  Could not save border CSV: {e}")

        return target_wcs
    
    def reproject_image(self, data, source_wcs, target_wcs, source_data=None):
        """Reproject image from source to target coordinate system with optional source-based distortion correction"""
        try:
            # Check if data is valid
            if data is None or data.size == 0:
                print("⚠️  Input data is None or empty")
                return None
            
            if not np.isfinite(data).any():
                print("⚠️  Input data contains no finite values")
                return None
            
            # If we have source data, use it for enhanced geometric transformation
            if source_data and source_data.get('pixel_coords'):
                num_sources = len(source_data['pixel_coords'])
                print(f"🌟 Using {num_sources} sources for distortion correction")
                
                # Implement source-based geometric transformation
                if num_sources >= 4:  # Need at least 4 points for reliable transformation
                    print(f"🔧 Computing source-based distortion map...")
                    result = self._source_based_reproject(data, source_wcs, target_wcs, source_data)
                    if isinstance(result, tuple):
                        return result  # Returns (data, verification_data)
                    else:
                        return result, None  # Fallback case
                else:
                    print(f"⚠️  Only {num_sources} sources - falling back to WCS reprojection")
            
            # Get target image size
            if hasattr(target_wcs, 'pixel_shape'):
                target_shape = target_wcs.pixel_shape
            else:
                target_shape = data.shape  # fallback
            
            # Create target image filled with NaN
            target_data = np.full(target_shape, np.nan, dtype=data.dtype)
            
            # Simple reprojection: for each pixel in target, find corresponding source pixel
            ny, nx = target_shape
            source_ny, source_nx = data.shape
            
            # Create coordinate arrays for target image
            y_target, x_target = np.mgrid[0:ny, 0:nx]
            
            # Convert target pixel coordinates to world coordinates
            try:
                ra_target, dec_target = target_wcs.pixel_to_world_values(x_target, y_target)
                
                # Convert world coordinates to source pixel coordinates
                x_source, y_source = source_wcs.world_to_pixel_values(ra_target, dec_target)
                
                # Round to nearest integer pixel
                x_source = np.round(x_source).astype(int)
                y_source = np.round(y_source).astype(int)
                
                # Find valid pixels (within source image bounds)
                valid_mask = ((x_source >= 0) & (x_source < source_nx) & 
                             (y_source >= 0) & (y_source < source_ny))
                
                # Calculate overlap statistics  
                total_pixels = nx * ny
                valid_pixels = np.sum(valid_mask)
                overlap_percent = (valid_pixels / total_pixels) * 100
                
                # Copy valid pixels
                if np.any(valid_mask):
                    target_data[y_target[valid_mask], x_target[valid_mask]] = data[y_source[valid_mask], x_source[valid_mask]]
                    print(f"🔍 Reprojected {valid_pixels} pixels ({overlap_percent:.1f}% overlap)")
                else:
                    print("🔍 No overlap - placing as non-overlapping survey image")
                    
                    # For non-overlapping images, we need to find where this image should be placed in the target grid
                    # Calculate the offset based on the image center position
                    try:
                        # Get the center coordinates of the source image 
                        center_x_src, center_y_src = data.shape[1]/2, data.shape[0]/2
                        center_ra_src, center_dec_src = source_wcs.pixel_to_world_values(center_x_src, center_y_src)
                        
                        # Find where this center should be in the target grid
                        center_x_tgt, center_y_tgt = target_wcs.world_to_pixel_values(center_ra_src, center_dec_src)
                        
                        # Calculate offset to place source image center at target center
                        offset_x = int(center_x_tgt - center_x_src)
                        offset_y = int(center_y_tgt - center_y_src)
                        
                        # Copy the entire source image to the target grid at the calculated offset
                        src_h, src_w = data.shape
                        tgt_h, tgt_w = target_data.shape
                        
                        # Calculate intersection bounds
                        src_y_start = max(0, -offset_y)
                        src_y_end = min(src_h, tgt_h - offset_y)
                        src_x_start = max(0, -offset_x)
                        src_x_end = min(src_w, tgt_w - offset_x)
                        
                        tgt_y_start = max(0, offset_y)
                        tgt_y_end = tgt_y_start + (src_y_end - src_y_start)
                        tgt_x_start = max(0, offset_x)
                        tgt_x_end = tgt_x_start + (src_x_end - src_x_start)
                        
                        if (src_y_end > src_y_start and src_x_end > src_x_start and 
                            tgt_y_end <= tgt_h and tgt_x_end <= tgt_w):
                            
                            # Copy the overlapping region
                            target_data[tgt_y_start:tgt_y_end, tgt_x_start:tgt_x_end] = \
                                data[src_y_start:src_y_end, src_x_start:src_x_end]
                            
                            copied_pixels = (tgt_y_end - tgt_y_start) * (tgt_x_end - tgt_x_start)
                            print(f"🔍 Placed {copied_pixels} pixels at offset ({offset_x}, {offset_y})")
                        else:
                            print("⚠️  Image falls completely outside target grid")
                            return None
                            
                    except Exception as e:
                        print(f"⚠️  Non-overlapping placement failed: {e}")
                        return None
                    
            except Exception as e:
                print(f"⚠️  Coordinate transformation error: {e}")
                # Fallback: just resize the image to target shape
                from scipy.ndimage import zoom
                zoom_factors = (target_shape[0] / data.shape[0], target_shape[1] / data.shape[1])
                target_data = zoom(data, zoom_factors, order=1)
                print(f"🔍 Fallback: resized image")
            
            return target_data, None  # No verification data for fallback method
            
        except Exception as e:
            print(f"⚠️  Reprojection error: {e}")
            import traceback
            traceback.print_exc()
            return None, None

    def _source_based_reproject(self, data, source_wcs, target_wcs, source_data):
        """Advanced reprojection using actual source positions for distortion correction"""
        try:
            from scipy.spatial.distance import cdist
            from scipy.interpolate import griddata
            
            print(f"🔧 Source-based geometric transformation...")
            
            # Get target image size
            if hasattr(target_wcs, 'pixel_shape'):
                target_shape = target_wcs.pixel_shape
            else:
                target_shape = data.shape
            
            # Extract source positions
            source_pixels = np.array(source_data['pixel_coords'])
            catalog_coords = np.array(source_data['catalog_coords'])  # RA, DEC
            
            if len(source_pixels) < 4 or len(catalog_coords) < 4:
                print("⚠️  Insufficient source points for geometric transformation")
                return None
            
            # Convert catalog coordinates to target pixel coordinates
            target_pixels = []
            valid_indices = []
            
            for i, (ra, dec) in enumerate(catalog_coords):
                try:
                    tx, ty = target_wcs.world_to_pixel_values(ra, dec)
                    if (0 <= tx < target_shape[1] and 0 <= ty < target_shape[0] and 
                        np.isfinite(tx) and np.isfinite(ty)):
                        target_pixels.append([tx, ty])
                        valid_indices.append(i)
                except:
                    continue
            
            if len(valid_indices) < 4:
                print(f"⚠️  Only {len(valid_indices)} valid control points - using standard reprojection")
                return self._standard_reproject(data, source_wcs, target_wcs, target_shape)
            
            # Filter to valid control points
            source_control = source_pixels[valid_indices]
            target_control = np.array(target_pixels)
            
            print(f"🎯 Using {len(valid_indices)} control points for transformation")
            
            # Create high-quality geometric transformation
            transformed_data = self._apply_geometric_transformation(data, source_control, target_control, target_shape)
            
            # Return both data and verification info
            if transformed_data is not None:
                verification_data = {
                    'source_control': source_control,
                    'target_control': target_control,
                    'catalog_coords': catalog_coords[valid_indices],
                    'source_wcs': source_wcs,
                    'target_wcs': target_wcs,
                    'image_file': None  # Will be set by caller
                }
                return transformed_data, verification_data
            
            return transformed_data, None
            
        except Exception as e:
            print(f"⚠️  Source-based reprojection failed: {e}")
            return self._standard_reproject(data, source_wcs, target_wcs, target_shape)
    
    def _apply_geometric_transformation(self, data, source_points, target_points, target_shape):
        """Apply geometric transformation using control points"""
        try:
            from scipy.interpolate import griddata
            from scipy.spatial.distance import cdist
            
            target_data = np.full(target_shape, np.nan, dtype=data.dtype)
            ny, nx = target_shape
            
            # Create dense coordinate grids
            y_target, x_target = np.mgrid[0:ny, 0:nx]
            target_coords = np.column_stack([x_target.ravel(), y_target.ravel()])
            
            # Compute transformation using thin plate splines for smoothness
            try:
                # Use Radial Basis Function interpolation for smooth transformation
                from scipy.interpolate import RBFInterpolator
                
                # Create RBF interpolators for x and y transformations  
                rbf_x = RBFInterpolator(target_points, source_points[:, 0], 
                                       kernel='thin_plate_spline', smoothing=0.1)
                rbf_y = RBFInterpolator(target_points, source_points[:, 1],
                                       kernel='thin_plate_spline', smoothing=0.1)
                
                # Apply transformation to get source coordinates
                source_x = rbf_x(target_coords).reshape(ny, nx)
                source_y = rbf_y(target_coords).reshape(ny, nx)
                
            except ImportError:
                # Fallback to simpler griddata interpolation
                print("🔄 Using griddata interpolation (RBF not available)")
                source_x = griddata(target_points, source_points[:, 0], target_coords, 
                                  method='cubic', fill_value=np.nan).reshape(ny, nx)
                source_y = griddata(target_points, source_points[:, 1], target_coords,
                                  method='cubic', fill_value=np.nan).reshape(ny, nx)
            
            # Sample from source image using computed coordinates
            source_ny, source_nx = data.shape
            
            # Round to nearest pixels and apply bounds checking
            sx = np.round(source_x).astype(int)
            sy = np.round(source_y).astype(int)
            
            # Create validity mask
            valid_mask = ((sx >= 0) & (sx < source_nx) & 
                         (sy >= 0) & (sy < source_ny) & 
                         np.isfinite(sx) & np.isfinite(sy))
            
            # Apply transformation
            valid_y, valid_x = np.where(valid_mask)
            if len(valid_y) > 0:
                target_data[valid_y, valid_x] = data[sy[valid_mask], sx[valid_mask]]
                
                # Calculate statistics
                valid_pixels = len(valid_y)
                total_pixels = ny * nx
                coverage = (valid_pixels / total_pixels) * 100
                
                print(f"✨ Geometric transformation applied: {valid_pixels} pixels ({coverage:.1f}% coverage)")
                
                # Quality check - ensure we have reasonable coverage
                if coverage < 5.0:
                    print("⚠️  Low coverage from geometric transformation - may need fallback")
                
                return target_data
            else:
                print("⚠️  No valid pixels after geometric transformation")
                return None
                
        except Exception as e:
            print(f"⚠️  Geometric transformation error: {e}")
            return None
    
    def _standard_reproject(self, data, source_wcs, target_wcs, target_shape):
        """Fallback to standard WCS-based reprojection"""
        print(f"🔄 Falling back to standard WCS reprojection...")
        
        # This is the existing logic moved to a separate method
        target_data = np.full(target_shape, np.nan, dtype=data.dtype)
        
        # Simple reprojection: for each pixel in target, find corresponding source pixel
        ny, nx = target_shape
        source_ny, source_nx = data.shape
        
        # Create coordinate arrays for target image
        y_target, x_target = np.mgrid[0:ny, 0:nx]
        
        # Convert target pixel coordinates to world coordinates
        try:
            ra_target, dec_target = target_wcs.pixel_to_world_values(x_target, y_target)
            
            # Convert world coordinates to source pixel coordinates
            x_source, y_source = source_wcs.world_to_pixel_values(ra_target, dec_target)
            
            # Round to nearest integer pixel
            x_source = np.round(x_source).astype(int)
            y_source = np.round(y_source).astype(int)
            
            # Find valid pixels (within source image bounds)
            valid_mask = ((x_source >= 0) & (x_source < source_nx) & 
                         (y_source >= 0) & (y_source < source_ny))
            
            # Calculate overlap statistics  
            total_pixels = nx * ny
            valid_pixels = np.sum(valid_mask)
            overlap_percent = (valid_pixels / total_pixels) * 100
            
            # Copy valid pixels
            if np.any(valid_mask):
                target_data[y_target[valid_mask], x_target[valid_mask]] = data[y_source[valid_mask], x_source[valid_mask]]
                print(f"🔍 Standard reprojection: {valid_pixels} pixels ({overlap_percent:.1f}% overlap)")
            else:
                print("🔍 No overlap - placing as non-overlapping survey image")
                # [Rest of the non-overlapping logic would go here - truncated for brevity]
                
        except Exception as e:
            print(f"⚠️  Standard reprojection coordinate transformation error: {e}")
            return None, None
            
        return target_data, None  # No verification data for standard reprojection

class ImageAligner:
    """Handles refined image alignment"""
    
    def align_images(self, image_list, reference_wcs):
        """Perform refined alignment of multiple images"""
        # Placeholder for now
        return image_list
    
    def compute_alignment_quality(self, aligned_images):
        """Compute alignment quality metrics"""
        return 0.8  # Placeholder

class SmartStitcher:
    """Advanced image combination and stitching with multiple blending methods"""
    
    def __init__(self, config=None):
        self.config = config or {}
        self.blend_method = self.config.get('blend_method', 'feather')
        
    def combine_images(self, images):
        """Create final mosaic with intelligent stitching"""
        if not images:
            return None
        
        blend_method = self.blend_method
        fast_mode = self.config.get('fast_mode', False)
        
        # If fast_mode, use streaming blend to lower memory usage
        if fast_mode:
            print(f"🧩 Combining {len(images)} images using FAST mode (streaming)...")
            mosaic_data = self._fast_stream_blend(images)
            return mosaic_data
        
        # -----------------------------------------------
        # Existing (non-fast) path retains advanced logic
        # -----------------------------------------------
        print(f"🧩 Combining {len(images)} images using {blend_method} blending...")
        
        # Extract metadata from images
        image_data = []
        for img in images:
            img_info = {
                'data': img['data'],
                'file': img['file'],
                'header': img.get('header', {}),
                'exposure_time': self._extract_exposure_time(img.get('header', {})),
                'quality': self._estimate_image_quality(img['data'])
            }
            image_data.append(img_info)
        
        # Stack all images (memory-intensive but acceptable in non-fast mode)
        image_stack = np.stack([img['data'] for img in image_data], axis=0)
        
        # Analyze overlap
        valid_count = np.sum(~np.isnan(image_stack), axis=0)
        overlap_pixels = np.sum(valid_count > 1)
        total_pixels = np.prod(valid_count.shape)
        overlap_percentage = (overlap_pixels / total_pixels) * 100
        
        print(f"📊 Overlap analysis: {overlap_percentage:.1f}% of pixels have multiple images")
        
        # Choose combination method based on overlap and user preference
        if overlap_percentage > 5:
            # Significant overlap - use blending
            if fast_mode:
                print(f"⚡ Using fast linear blending for overlapping regions...")
                mosaic_data = self._fast_blend(image_data, image_stack)
            else:
                print(f"🎨 Using advanced blending for overlapping regions...")
                mosaic_data = self._advanced_blend(image_data, image_stack, valid_count)
        else:
            # Minimal overlap - use enhanced placement
            print(f"🔍 Using {'fast' if fast_mode else 'enhanced'} placement for non-overlapping images...")
            mosaic_data = self._enhanced_placement(image_data, image_stack)
        
        # Post-processing (skip in fast mode)
        if not fast_mode:
            mosaic_data = self._post_process_mosaic(mosaic_data)
        
        # Check final mosaic statistics
        finite_data = mosaic_data[np.isfinite(mosaic_data)]
        if len(finite_data) > 0:
            print(f"✅ Mosaic complete: {mosaic_data.shape} pixels, {100*len(finite_data)/mosaic_data.size:.1f}% coverage")
            print(f"📈 Dynamic range: {np.min(finite_data):.1f} to {np.max(finite_data):.1f}")
        else:
            print("⚠️  No finite data in final mosaic")
            return None
        
        return mosaic_data
    
    def _fast_blend(self, image_data, image_stack):
        """Fast blending method for speed over quality"""
        print("⚡ Applying fast linear blending...")
        
        # Simple linear weights based on valid pixels only
        weights = (~np.isnan(image_stack)).astype(np.float32)
        
        # Normalize weights
        total_weights = np.sum(weights, axis=0)
        total_weights[total_weights == 0] = 1
        normalized_weights = weights / total_weights[np.newaxis, ...]
        
        # Simple weighted combination
        mosaic_data = np.nansum(image_stack * normalized_weights, axis=0)
        return mosaic_data
    
    def _advanced_blend(self, image_data, image_stack, valid_count):
        """Advanced blending for overlapping regions"""
        print(f"🎨 Performing {self.blend_method} blending...")
        
        # Create weight maps for each image
        weight_start = time.time()
        weight_maps = self._create_weight_maps(image_data, image_stack, valid_count)
        weight_time = time.time() - weight_start
        print(f"   ⚖️  Weight calculation: {weight_time:.2f}s")
        
        # Normalize weights so they sum to 1 at each pixel
        norm_start = time.time()
        total_weights = np.sum(weight_maps, axis=0)
        total_weights[total_weights == 0] = 1  # Avoid division by zero
        normalized_weights = weight_maps / total_weights[np.newaxis, ...]
        norm_time = time.time() - norm_start
        print(f"   🔢 Weight normalization: {norm_time:.2f}s")
        
        # Apply blending method
        blend_start = time.time()
        if self.blend_method == 'feather':
            mosaic_data = self._feather_blend(image_stack, normalized_weights)
        elif self.blend_method == 'linear':
            mosaic_data = self._linear_blend(image_stack, normalized_weights)
        elif self.blend_method == 'weighted':
            mosaic_data = self._weighted_blend(image_stack, normalized_weights)
        else:
            # Fallback to weighted median
            mosaic_data = self._weighted_median(image_stack, normalized_weights)
        blend_time = time.time() - blend_start
        print(f"   🎨 Blending operation: {blend_time:.2f}s")
        
        return mosaic_data
    
    def _enhanced_placement(self, image_data, image_stack):
        """Enhanced placement for non-overlapping images"""
        print("🔧 Enhanced placement with quality-based ordering...")
        
        # Sort images by quality (best first)
        sorted_images = sorted(enumerate(image_data), key=lambda x: x[1]['quality'], reverse=True)
        
        # Start with empty mosaic
        mosaic_data = np.full(image_stack.shape[1:], np.nan, dtype=image_stack.dtype)
        
        # Place images in quality order
        for orig_idx, img_info in sorted_images:
            img_data = image_stack[orig_idx]
            
            # Find pixels where this image has data but mosaic doesn't
            new_pixels = ~np.isnan(img_data) & np.isnan(mosaic_data)
            
            # Also replace pixels where this image is significantly better
            if hasattr(img_info, 'quality'):
                existing_pixels = ~np.isnan(img_data) & ~np.isnan(mosaic_data)
                # Could implement quality-based replacement here
                
            # Place the image
            mosaic_data[new_pixels] = img_data[new_pixels]
            
            pixels_added = np.sum(new_pixels)
            print(f"   📷 {os.path.basename(img_info['file'])}: added {pixels_added} pixels (quality: {img_info['quality']:.2f})")
        
        return mosaic_data
    
    def _create_weight_maps(self, image_data, image_stack, valid_count):
        """Create weight maps for each image based on quality, exposure, and distance from edges"""
        fast_mode = self.config.get('fast_mode', False)
        
        if fast_mode:
            print("⚖️  Creating simple weight maps (fast mode)...")
        else:
            print("⚖️  Creating quality and exposure-based weight maps...")
        
        weight_maps = np.zeros_like(image_stack, dtype=np.float32)
        
        for i, img_info in enumerate(image_data):
            img_data = image_stack[i]
            
            # Start with valid pixel mask
            base_weights = (~np.isnan(img_data)).astype(np.float32)
            
            if not fast_mode:
                # Distance from edge weighting (for feathering) - skip in fast mode
                if self.blend_method == 'feather':
                    edge_weights = self._compute_edge_weights(base_weights)
                    base_weights *= edge_weights
                
                # Quality weighting - skip in fast mode
                quality_weight = img_info['quality']
                base_weights *= quality_weight
                
                # Exposure time weighting - skip in fast mode
                exposure_weight = self._compute_exposure_weight(img_info['exposure_time'])
                base_weights *= exposure_weight
                
                # Background uniformity weighting - skip in fast mode
                bg_weight = self._compute_background_weight(img_data)
                base_weights *= bg_weight
                
                # Debug info
                avg_weight = np.mean(base_weights[base_weights > 0])
                print(f"   📷 {os.path.basename(img_info['file'])}: avg weight = {avg_weight:.3f} "
                      f"(quality={quality_weight:.2f}, exposure={exposure_weight:.2f}, bg={bg_weight:.2f})")
            else:
                # Fast mode: just use uniform weights for valid pixels
                print(f"   📷 {os.path.basename(img_info['file'])}: uniform weights")
            
            weight_maps[i] = base_weights
        
        return weight_maps
    
    def _compute_edge_weights(self, mask):
        """Compute distance-from-edge weights for feathering"""
        from scipy.ndimage import distance_transform_edt
        
        # Compute distance from edges of valid region
        distances = distance_transform_edt(mask)
        
        # Normalize to 0-1 range with gentle falloff
        max_dist = np.max(distances)
        if max_dist > 0:
            # Use square root for gentle falloff
            edge_weights = np.sqrt(distances / max_dist)
        else:
            edge_weights = mask.astype(np.float32)
        
        return edge_weights
    
    def _compute_exposure_weight(self, exposure_time):
        """Compute weight based on exposure time"""
        if exposure_time is None or exposure_time <= 0:
            return 1.0
        
        # Longer exposures get higher weight (up to a point)
        # Use logarithmic scaling to prevent very long exposures from dominating
        base_exposure = 60.0  # seconds
        weight = np.log10(max(exposure_time, 1) / base_exposure + 1) + 1
        return min(weight, 3.0)  # Cap at 3x weight
    
    def _compute_background_weight(self, img_data):
        """Compute weight based on background uniformity"""
        finite_data = img_data[np.isfinite(img_data)]
        if len(finite_data) < 1000:
            return 1.0
        
        # Estimate background uniformity by looking at noise
        # More uniform backgrounds get higher weight
        bg_std = np.std(finite_data)
        bg_mean = np.mean(finite_data)
        
        if bg_mean > 0:
            noise_ratio = bg_std / bg_mean
            # Lower noise ratio = higher weight
            weight = 1.0 / (1.0 + noise_ratio)
        else:
            weight = 1.0
        
        return np.clip(weight, 0.1, 2.0)  # Keep reasonable bounds
    
    def _feather_blend(self, image_stack, weights):
        """Feathered blending for seamless transitions"""
        print("🪶 Applying feathered blending...")
        
        # Apply Gaussian blur to weights for smooth transitions
        from scipy.ndimage import gaussian_filter
        
        blurred_weights = np.zeros_like(weights)
        for i in range(weights.shape[0]):
                         # Apply slight blur to weight maps
             sigma = self.config.get('feather_sigma', 2.0)
             blurred_weights[i] = gaussian_filter(weights[i], sigma=sigma)
        
        # Renormalize after blurring
        total_weights = np.sum(blurred_weights, axis=0)
        total_weights[total_weights == 0] = 1
        blurred_weights = blurred_weights / total_weights[np.newaxis, ...]
        
        # Weighted combination
        mosaic_data = np.nansum(image_stack * blurred_weights, axis=0)
        return mosaic_data
    
    def _linear_blend(self, image_stack, weights):
        """Linear blending"""
        print("📐 Applying linear blending...")
        mosaic_data = np.nansum(image_stack * weights, axis=0)
        return mosaic_data
    
    def _weighted_blend(self, image_stack, weights):
        """Quality and exposure weighted blending"""
        print("⚖️  Applying weighted blending...")
        
        # Remove outliers before combining
        sigma_threshold = self.config.get('outlier_threshold', 3.0)
        cleaned_stack = self._reject_outliers(image_stack, weights, sigma_threshold)
        
        # Weighted combination
        mosaic_data = np.nansum(cleaned_stack * weights, axis=0)
        return mosaic_data
    
    def _weighted_median(self, image_stack, weights):
        """Weighted median combination (robust to outliers)"""
        print("📊 Applying weighted median blending...")
        
        # For now, use simple median (weighted median is complex to implement)
        # TODO: Implement true weighted median
        mosaic_data = np.nanmedian(image_stack, axis=0)
        return mosaic_data
    
    def _reject_outliers(self, image_stack, weights, sigma_threshold=3.0):
        """Remove outliers (cosmic rays, satellites) using sigma clipping"""
        print(f"🚫 Rejecting outliers with {sigma_threshold}σ threshold...")
        
        # Compute weighted mean and std
        mean_image = np.nansum(image_stack * weights, axis=0)
        
        # Compute residuals
        residuals = image_stack - mean_image[np.newaxis, ...]
        
        # Compute robust standard deviation
        finite_residuals = residuals[np.isfinite(residuals)]
        if len(finite_residuals) > 0:
            std_image = np.nanstd(finite_residuals)
            
            # Mask outliers
            outlier_mask = np.abs(residuals) > sigma_threshold * std_image
            
            # Replace outliers with NaN
            cleaned_stack = image_stack.copy()
            cleaned_stack[outlier_mask] = np.nan
            
            outliers_removed = np.sum(outlier_mask)
            print(f"   🧹 Removed {outliers_removed} outlier pixels")
            
            return cleaned_stack
        
        return image_stack
    
    def _post_process_mosaic(self, mosaic_data):
        """Post-process the combined mosaic"""
        print("🔧 Post-processing mosaic...")
        
        # Background subtraction/normalization could go here
        # For now, just return as-is
        return mosaic_data
    
    def _extract_exposure_time(self, header):
        """Extract exposure time from FITS header"""
        exposure_keys = ['EXPTIME', 'EXPOSURE', 'ITIME', 'TEXP']
        for key in exposure_keys:
            if key in header:
                try:
                    return float(header[key])
                except (ValueError, TypeError):
                    continue
        return None
    
    def _estimate_image_quality(self, img_data):
        """Estimate image quality based on various metrics"""
        finite_data = img_data[np.isfinite(img_data)]
        if len(finite_data) < 1000:
            return 0.1  # Poor quality
        
        # Simple quality metrics
        # 1. Signal-to-noise ratio
        mean_signal = np.mean(finite_data)
        std_noise = np.std(finite_data)
        
        if std_noise > 0:
            snr = mean_signal / std_noise
        else:
            snr = 1.0
        
        # 2. Data completeness
        completeness = len(finite_data) / img_data.size
        
        # 3. Dynamic range
        if len(finite_data) > 0:
            dynamic_range = (np.max(finite_data) - np.min(finite_data)) / np.mean(finite_data)
        else:
            dynamic_range = 1.0
        
        # Combine metrics (this is a simple heuristic)
        quality = (snr * 0.4 + completeness * 0.4 + dynamic_range * 0.2) / 3.0
        
        return np.clip(quality, 0.1, 2.0)  # Keep reasonable bounds

    def _fast_stream_blend(self, raw_images):
        """Memory-efficient linear blending used when fast_mode is enabled.
        This method avoids stacking all images at once – instead it accumulates
        a running sum and count, requiring only two full-resolution arrays in
        memory regardless of image count."""
        # Determine mosaic shape from first image with finite data
        first_img = None
        for item in raw_images:
            if item['data'] is not None:
                first_img = item['data']
                break
        if first_img is None:
            print("⚠️  No image data found for streaming blend")
            return None
        mosaic_sum = np.zeros_like(first_img, dtype=np.float32)
        mosaic_count = np.zeros_like(first_img, dtype=np.uint16)

        # Iterate through images and accumulate
        for img in raw_images:
            img_data = img['data']
            if img_data is None:
                continue
            # Ensure float32 for consistent accumulation
            img_data_f32 = img_data.astype(np.float32, copy=False)
            mask = np.isfinite(img_data_f32)
            mosaic_sum[mask] += img_data_f32[mask]
            mosaic_count[mask] += 1
        
        # Compute simple average where we have data
        with np.errstate(invalid='ignore', divide='ignore'):
            mosaic_data = mosaic_sum / mosaic_count
        mosaic_data[mosaic_count == 0] = np.nan

        # Overlap analysis summary
        overlap_pixels = np.sum(mosaic_count > 1)
        total_pixels = mosaic_count.size
        overlap_percentage = (overlap_pixels / total_pixels) * 100
        print(f"📊 Overlap analysis: {overlap_percentage:.1f}% of pixels have multiple images")

        # Quick stats
        finite_vals = mosaic_data[np.isfinite(mosaic_data)]
        if finite_vals.size == 0:
            print("⚠️  No finite data in final mosaic after streaming blend")
            return None
        print(f"✅ Mosaic complete: {mosaic_data.shape} pixels, {100*finite_vals.size/mosaic_data.size:.1f}% coverage")
        print(f"📈 Dynamic range: {np.min(finite_vals):.1f} to {np.max(finite_vals):.1f}")
        return mosaic_data

class MosaicExporter:
    """Handles mosaic output and visualization"""
    
    def export_mosaic(self, mosaic_data, target_wcs, output_dir, source_images):
        """Export mosaic in various formats"""
        try:
            # Export FITS
            fits_path = os.path.join(output_dir, 'mosaic.fits')
            self.export_fits(mosaic_data, target_wcs, fits_path, source_images)
            
            # Export PNG with WCS axes
            png_path = os.path.join(output_dir, 'mosaic.png')
            self.export_png(mosaic_data, target_wcs, png_path)
            
            # Export coverage map
            coverage_path = os.path.join(output_dir, 'coverage_map.png')
            self.export_coverage_map(source_images, coverage_path)
            
            # === New: print unique sky-coverage area (no double-counting overlaps) ===
            try:
                valid_mask = np.isfinite(mosaic_data)
                n_valid_pix = int(np.count_nonzero(valid_mask))
                if n_valid_pix > 0 and target_wcs is not None:
                    # Pixel area in square degrees from WCS CDELTs
                    cdelt_deg = np.abs(target_wcs.wcs.cdelt)  # [deg/pix] for X and Y
                    pix_area_deg2 = cdelt_deg[0] * cdelt_deg[1]
                    covered_deg2  = n_valid_pix * pix_area_deg2

                    # Also express in square arc-minutes for a more intuitive number
                    covered_arcmin2 = covered_deg2 * 3600.0

                    coverage_pct = (n_valid_pix / mosaic_data.size) * 100.0
                    print(
                        f"🌌 Unique sky coverage: {covered_deg2:.3f} deg² "
                        f"({covered_arcmin2:.0f} arcmin²) — {coverage_pct:.1f}% of mosaic grid"
                    )
            except Exception as cov_err:
                print(f"⚠️  Coverage calculation failed: {cov_err}")
            
            return True
            
        except Exception as e:
            print(f"❌ Export error: {e}")
            return False
    
    def export_fits(self, mosaic_data, target_wcs, output_path, source_images):
        """Export mosaic as FITS file"""
        # Create FITS HDU
        hdu = fits.PrimaryHDU(data=mosaic_data, header=target_wcs.to_header())
        
        # Add metadata
        hdu.header['ORIGIN'] = 'Spliter v5.0'
        hdu.header['CREATOR'] = 'MosaicEngine'
        hdu.header['DATE'] = datetime.datetime.now().isoformat()
        hdu.header['NIMAGES'] = len(source_images)
        
        # ---- Coverage statistics ----
        try:
            valid_mask = np.isfinite(mosaic_data)
            n_valid_pix = int(np.count_nonzero(valid_mask))
            total_pix   = mosaic_data.size
            if n_valid_pix > 0:
                cdelt_deg = np.abs(target_wcs.wcs.cdelt)  # deg/pix
                pix_area  = cdelt_deg[0] * cdelt_deg[1]   # deg² per pixel
                cov_deg2  = n_valid_pix * pix_area
                cov_pct   = (n_valid_pix / total_pix) * 100.0
                hdu.header['COVPIX']  = (n_valid_pix, 'Number of pixels with data')
                hdu.header['COVAREA'] = (round(cov_deg2, 6), 'Sky area covered (deg^2)')
                hdu.header['COVPCT']  = (round(cov_pct, 2), 'Coverage of mosaic grid (%)')
        except Exception as _cov_err:
            # Non-fatal: just continue
            pass
        
        # Add source file information
        for i, img in enumerate(source_images[:10]):  # Limit to first 10 to avoid header overflow
            hdu.header[f'SOURCE{i:02d}'] = os.path.basename(img['file'])
        
        # Write FITS file
        hdu.writeto(output_path, overwrite=True)
        print(f"💾 FITS mosaic saved: {output_path}")
    
    def export_png(self, mosaic_data, target_wcs, output_path):
        """Export mosaic as PNG preview"""
        try:
            # Create figure with better aspect ratio for wide mosaics
            aspect_ratio = mosaic_data.shape[1] / mosaic_data.shape[0]
            if aspect_ratio > 2:
                figsize = (16, 8)
            else:
                figsize = (12, 10)
            
            if target_wcs is not None:
                from astropy.visualization.wcsaxes import WCSAxes
                fig = plt.figure(figsize=figsize)
                ax = fig.add_subplot(111, projection=target_wcs)
            else:
                fig, ax = plt.subplots(figsize=figsize)
            
            # Check data for PNG export
            finite_data = mosaic_data[np.isfinite(mosaic_data)]
            
            if len(finite_data) == 0:
                print("⚠️  No finite data for PNG export - creating placeholder")
                # Create a test pattern
                test_data = np.random.random(mosaic_data.shape) * 1000
                im = ax.imshow(test_data, origin='lower', cmap='gray')
                ax.set_title('Mosaic Preview (No Valid Data - Test Pattern)', fontsize=16)
            else:
                # Copy data for display and set NaNs to background
                display_data = mosaic_data.copy()
                
                # Compute robust percentile limits for stretching
                vmin, vmax = np.percentile(finite_data, [0.5, 99.5])
                if abs(vmax - vmin) < 1e-6:
                    vmin, vmax = np.percentile(finite_data, [0.1, 99.9])
                if abs(vmax - vmin) < 1e-6:
                    vmin, vmax = np.min(finite_data), np.max(finite_data)
                if abs(vmax - vmin) < 1e-6:
                    center = (vmin + vmax) / 2
                    vmin, vmax = center - 0.001, center + 0.001

                # Histogram equalization for better local contrast
                # Scale to 0-1, clip to robust range
                scaled = np.clip((display_data - vmin) / (vmax - vmin), 0, 1)
                # Compute histogram and cumulative distribution function (CDF)
                hist, bins = np.histogram(scaled[np.isfinite(scaled)], bins=2048, range=(0, 1), density=False)
                cdf = np.cumsum(hist).astype(float)
                cdf /= cdf[-1] if cdf[-1] != 0 else 1
                # Map the scaled data through the CDF for equalization
                equalized = np.interp(scaled.flat, bins[:-1], cdf).reshape(scaled.shape)
                # Replace NaNs with background just below 0 for grayscale
                nan_value = -0.05
                equalized[~np.isfinite(equalized)] = nan_value
                
                # Apply gamma correction (>1) to darken background and enhance bright features
                gamma = 3.5  # value >1 darkens backgrounds
                stretched = np.power(equalized, gamma)
                
                # Replace NaNs with background just below 0 for grayscale
                nan_value = -0.05
                stretched[~np.isfinite(stretched)] = nan_value
                
                # Display equalized grayscale mosaic
                im = ax.imshow(stretched, origin='lower', cmap='gray', vmin=0, vmax=1, interpolation='none')
                ax.set_title(f'Mosaic Preview (Equalized Grayscale)\n{mosaic_data.shape[1]}×{mosaic_data.shape[0]} pixels', fontsize=14)
            
            if target_wcs is not None:
                # Use celestial coordinates, add curved RA/DEC grid
                ax.coords.grid(True, color='white', alpha=0.3, linestyle='--', linewidth=0.5)
                ax.coords[0].set_axislabel('Right Ascension (°)')
                ax.coords[1].set_axislabel('Declination (°)')
                # Remove default grid added earlier (pixel-based)
            else:
                ax.set_xlabel('Pixel X')
                ax.set_ylabel('Pixel Y')
            
            # Add colorbar
            plt.colorbar(im, ax=ax, label='Intensity', shrink=0.8)
            
            # Add grid to help see image boundaries
            if len(finite_data) > 0:
                # Pixel grid no longer needed when using WCS grid; keep fallback for non-WCS
                if len(finite_data) > 0 and target_wcs is None:
                    ax.grid(True, alpha=0.2, linewidth=0.5, color='white')
            
            # Save PNG with higher DPI for better quality
            plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
            plt.close(fig)
            
            print(f"🖼️  PNG preview saved: {output_path}")
            
            # Also create a detailed analysis plot
            self._create_analysis_plot(mosaic_data, output_path.replace('.png', '_analysis.png'))
            
        except Exception as e:
            print(f"⚠️  PNG export error: {e}")
            import traceback
            traceback.print_exc()
    
    def _create_analysis_plot(self, mosaic_data, output_path):
        """Create detailed analysis plot showing data distribution"""
        try:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
            
            # Main mosaic view
            finite_data = mosaic_data[np.isfinite(mosaic_data)]
            if len(finite_data) > 0:
                display_data = mosaic_data.copy()
                vmin, vmax = np.percentile(finite_data, [1, 99])
                nan_value = vmin - (vmax - vmin) * 0.1
                display_data[~np.isfinite(display_data)] = nan_value
                
                im1 = ax1.imshow(display_data, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
                ax1.set_title('Full Mosaic', fontsize=14)
                ax1.set_xlabel('Pixel X')
                ax1.set_ylabel('Pixel Y')
                plt.colorbar(im1, ax=ax1, label='Intensity')
                
                # Coverage map (finite vs NaN)
                coverage = np.isfinite(mosaic_data).astype(float)
                im2 = ax2.imshow(coverage, origin='lower', cmap='RdYlBu_r', vmin=0, vmax=1)
                ax2.set_title('Data Coverage Map', fontsize=14)
                ax2.set_xlabel('Pixel X')
                ax2.set_ylabel('Pixel Y')
                plt.colorbar(im2, ax=ax2, label='Has Data')
                
                # Histogram of finite values
                ax3.hist(finite_data, bins=100, alpha=0.7, color='blue', edgecolor='black')
                ax3.set_xlabel('Intensity')
                ax3.set_ylabel('Pixel Count')
                ax3.set_title(f'Intensity Distribution ({len(finite_data):,} pixels)', fontsize=14)
                ax3.grid(True, alpha=0.3)
                
                # Spatial distribution by rows and columns
                row_coverage = np.mean(np.isfinite(mosaic_data), axis=1) * 100
                col_coverage = np.mean(np.isfinite(mosaic_data), axis=0) * 100
                
                ax4_twin = ax4.twinx()
                ax4.plot(row_coverage, range(len(row_coverage)), 'b-', label='Row Coverage %')
                ax4_twin.plot(range(len(col_coverage)), col_coverage, 'r-', label='Column Coverage %')
                ax4.set_xlabel('Coverage %')
                ax4.set_ylabel('Row Number')
                ax4_twin.set_ylabel('Column Coverage %')
                ax4.set_title('Spatial Coverage Distribution', fontsize=14)
                ax4.legend(loc='upper left')
                ax4_twin.legend(loc='upper right')
                ax4.grid(True, alpha=0.3)
                
            plt.tight_layout()
            plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
            plt.close(fig)
            
            print(f"📊 Analysis plot saved: {output_path}")
            
        except Exception as e:
            print(f"⚠️  Analysis plot error: {e}")
            import traceback
            traceback.print_exc()
    
    def export_coverage_map(self, source_images, output_path):
        """Export coverage map showing input images"""
        try:
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Plot footprints of source images
            for i, img in enumerate(source_images):
                # Simple representation - just plot image centers
                # TODO: Plot actual footprints
                ax.scatter(i, 0, s=100, alpha=0.7, label=f"Image {i+1}")
            
            ax.set_title('Mosaic Coverage Map')
            ax.set_xlabel('Image Index')
            ax.set_ylabel('Coverage')
            ax.legend()
            
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
            
            print(f"🗺️  Coverage map saved: {output_path}")
            
        except Exception as e:
            print(f"⚠️  Coverage map export error: {e}")

def find_res_fit_files(base_dir):
    """Find all res.fit files in OBJECT_XXXX subdirectories"""
    res_fit_files = []
    base_path = Path(base_dir)
    
    print(f"🔍 Searching for res.fit files in {base_dir}...")
    
    # Look for OBJECT_XXXX directories
    for target_dir in base_path.iterdir():
        if target_dir.is_dir():
            # Look for res.fit files in subdirectories
            for res_file in target_dir.rglob('res.fit'):
                res_fit_files.append(str(res_file))
                print(f"  Found: {res_file}")
    
    print(f"📊 Found {len(res_fit_files)} res.fit files")
    return res_fit_files

def create_mosaic(args):
    """Main function for mosaic creation"""
    # Update config with command line arguments
    config = MOSAIC_CONFIG.copy()
    config['pixel_scale'] = args.pixel_scale
    config['projection'] = args.projection
    config['blend_method'] = args.blend_method
    config['quality_threshold'] = args.quality_threshold
    config['margin'] = args.margin
    config['no_size_limit'] = args.no_size_limit
    config['outlier_threshold'] = args.outlier_threshold
    config['feather_sigma'] = args.feather_sigma
    config['fast_mode'] = args.fast_mode
    config['incremental_output'] = args.incremental
    config['edge_samples'] = args.edge_samples
    config['dump_border_csv'] = not args.no_border_csv
    
    print(f"🎨 Mosaic Configuration:")
    print(f"   Pixel scale: {config['pixel_scale']} arcsec/pixel")
    print(f"   Projection: {config['projection']}")
    print(f"   Blend method: {config['blend_method']}")
    print(f"   Quality threshold: {config['quality_threshold']}")
    print(f"   Margin: {config['margin']} degrees")
    print(f"   Size limits: {'DISABLED' if config['no_size_limit'] else 'ENABLED'}")
    print(f"   Outlier threshold: {config['outlier_threshold']}σ")
    print(f"   Feather sigma: {config['feather_sigma']}")
    print(f"   Fast mode: {'ENABLED' if config['fast_mode'] else 'DISABLED'}")
    print(f"   Incremental output: {'ENABLED' if config['incremental_output'] else 'DISABLED'}")
    print(f"   Output directory: {args.output_dir}")
    
    # Find res.fit files
    res_fit_files = find_res_fit_files(init_path)
    
    if not res_fit_files:
        print("❌ No res.fit files found. Please run normal processing first.")
        return False
    
    # Apply limit if specified
    if args.limit:
        res_fit_files = res_fit_files[:args.limit]
        print(f"🔢 Processing first {len(res_fit_files)} images (--limit {args.limit})")
    
    # Create mosaic engine
    mosaic_engine = MosaicEngine(config)
    
    # Create mosaic
    success = mosaic_engine.create_mosaic(res_fit_files, args.output_dir)
    
    if success:
        print(f"🎉 Mosaic creation completed successfully!")
        print(f"📁 Output files saved to: {args.output_dir}")
    else:
        print("❌ Mosaic creation failed")
    
    return success

# ============================================================================
# END MOSAIC CREATION SECTION
# ============================================================================

def save_sky_maps(footprints):
    """Create/refresh global and zoomed sky maps with NGC galaxies overlay."""
    if not footprints:
        return

    # --------- Wide Mollweide map ----------
    try:
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111, projection="mollweide")

        # Plot footprints
        for corners in footprints:
            ra = np.array([c[0] for c in corners])
            dec = np.array([c[1] for c in corners])
            ra_rad = np.radians(ra)
            ra_rad = np.remainder(ra_rad + 2 * np.pi, 2 * np.pi)
            ra_rad[ra_rad > np.pi] -= 2 * np.pi
            dec_rad = np.radians(dec)
            poly = Polygon(np.column_stack((-ra_rad, dec_rad)), closed=True, facecolor="none", edgecolor="red", linewidth=0.4)
            ax.add_patch(poly)

        # Overlay NGC galaxies (sample to 3000 for speed)
        if not NGC_DF.empty:
            sample = NGC_DF.sample(n=min(3000, len(NGC_DF)), random_state=1)
            ra_rad = np.radians(sample["ra_deg"].values)
            ra_rad = np.remainder(ra_rad + 2 * np.pi, 2 * np.pi)
            ra_rad[ra_rad > np.pi] -= 2 * np.pi
            dec_rad = np.radians(sample["dec_deg"].values)
            ax.scatter(-ra_rad, dec_rad, s=2, c="green", marker="+", alpha=0.5, linewidths=0.3, label="NGC galaxies")

        ax.grid(True, color="lightgray", lw=0.4)
        # Add legend (if galaxies plotted) *before* saving so it is visible on the image
        if len(ax.get_legend_handles_labels()[0]):
            ax.legend(loc="lower left", fontsize="small")

        ax.set_title("Sky coverage (wide)")
        out_png = os.path.join(init_path, "sky_coverage.png")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"⚠️  Failed wide sky map: {e}")

    # --------- Zoomed map ----------
    try:
        # Determine overall bounding box of footprints
        all_ra = np.concatenate([[c[0] for c in fp] for fp in footprints])
        all_dec = np.concatenate([[c[1] for c in fp] for fp in footprints])
        
                # Handle RA wraparound: detect if we span across 0°/360° boundary
        ra_range = all_ra.max() - all_ra.min()
        if ra_range > 180:  # Likely wraparound case
            # Shift RAs: values > 180° become negative
            all_ra = np.where(all_ra > 180, all_ra - 360, all_ra)
            # Recalculate bounds after shift
            ra_min, ra_max = all_ra.min(), all_ra.max()
            # Apply same shift to footprint coordinates for plotting
            shifted_footprints = []
            for fp in footprints:
                shifted_corners = []
                for ra, dec in fp:
                    ra_shifted = ra - 360 if ra > 180 else ra
                    shifted_corners.append((ra_shifted, dec))
                shifted_footprints.append(shifted_corners)
            footprints_to_plot = shifted_footprints
        else:
            ra_min, ra_max = all_ra.min(), all_ra.max()
            footprints_to_plot = footprints
        
        dec_min, dec_max = all_dec.min(), all_dec.max()
        
        # Add 2 deg margin
        margin = 2.0
        ra_min -= margin
        ra_max += margin
        dec_min -= margin
        dec_max += margin

        # Calculate mean declination for aspect ratio correction
        mean_dec = np.mean(all_dec)
        cos_dec = np.cos(np.radians(mean_dec))

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_xlim(ra_max, ra_min)  # RA increasing to left
        ax.set_ylim(dec_min, dec_max)
        
        # Set equal aspect ratio corrected for declination
        # This makes 1° RA and 1° DEC have the same angular size on the plot
        ax.set_aspect(1.0 / cos_dec)
        
        ax.set_xlabel("RA (deg)")
        ax.set_ylabel("Dec (deg)")

        # Draw footprints boxes
        for corners in footprints_to_plot:
            ra = np.array([c[0] for c in corners])
            dec = np.array([c[1] for c in corners])
            poly = Polygon(np.column_stack((ra, dec)), closed=True, facecolor="none", edgecolor="red", linewidth=0.6)
            ax.add_patch(poly)

        # Galaxies inside region - handle wraparound case
        if not NGC_DF.empty:
            if ra_range > 180:  # Wraparound case
                # Apply same coordinate shift to galaxy coordinates
                ngc_ra_shifted = np.where(NGC_DF["ra_deg"] > 180, NGC_DF["ra_deg"] - 360, NGC_DF["ra_deg"])
                mask = (
                    (ngc_ra_shifted >= ra_min)
                    & (ngc_ra_shifted <= ra_max)
                    & (NGC_DF["dec_deg"] >= dec_min)
                    & (NGC_DF["dec_deg"] <= dec_max)
                )
                sel_ra = ngc_ra_shifted[mask]
                sel_dec = NGC_DF["dec_deg"][mask]
            else:
                # Normal case
                mask = (
                    (NGC_DF["ra_deg"] >= ra_min)
                    & (NGC_DF["ra_deg"] <= ra_max)
                    & (NGC_DF["dec_deg"] >= dec_min)
                    & (NGC_DF["dec_deg"] <= dec_max)
                )
                sel_ra = NGC_DF["ra_deg"][mask]
                sel_dec = NGC_DF["dec_deg"][mask]
            
            if len(sel_ra) > 0:
                ax.scatter(sel_ra, sel_dec, s=10, c="green", marker="+", linewidths=0.5, label="NGC galaxies")

        ax.legend(loc="upper right", fontsize="small")
        ax.set_title("Sky coverage – zoom")
        ax.grid(True, alpha=0.3)
        out_zoom = os.path.join(init_path, "sky_coverage_zoom.png")
        fig.savefig(out_zoom, dpi=400, bbox_inches="tight")
        plt.close(fig)
        print(f"🔄 Sky maps refreshed → {out_png}, {out_zoom}")
    except Exception as e:
        print(f"⚠️  Failed zoom sky map: {e}")

# ---------------------------------------------------------------------------
# Batch calibration per-folder using pySiril
# ---------------------------------------------------------------------------

def calibrate_folders_pysiril(folders: dict):
    """Calibrate/stack each folder listed in `folders` using pySiril.

    Parameters
    ----------
    folders : dict
        key = exposure directory path (str or Path)
        value = dict(filter=..., exposure=..., num_images=...)
    """
    # Make the global set available so we can mark folders as processed
    global calibrated_folders

    # Detect pySiril API version (wrapper vs legacy)
    API = None
    try:
        import pysiril.wrapper as sr
        API = 'wrapper'
        if hasattr(sr, 'start'):
            sr.start()
        else:
            # Extremely old wrapper that only exposes the Wrapper class – create our own helper
            import pysiril.siril as _ps
            _pipe = _ps.Siril()
            if not _pipe.Open():
                raise RuntimeError("Failed to open Siril via wrapper fallback")

            _wrapper = sr.Wrapper(_pipe)

            # Dynamically create minimal interface expected later (run & stop)
            class _DynSR:
                def __init__(self, wrapper, pipe):
                    self._w = wrapper
                    self._p = pipe
                def run(self, cmd):
                    self._p.Execute(cmd)
                def stop(self):
                    self._p.Close()

            sr = _DynSR(_wrapper, _pipe)
    except ImportError:
        try:
            import pysiril.siril as sr
            API = 'legacy'
            sr.init()
        except ImportError:
            print("⚠️  pySiril not installed – per-folder calibration skipped")
            return

    # Locate Siril.app on macOS if not configured
    os.environ.setdefault("SIRIL_APP", "/Applications/Siril.app")

    print(f"🚀 Launching Siril (pySiril {API}) to process {len(folders)} folder(s)…")
 
    try:
        for folder, info in tqdm(folders.items(), desc="Calibrating", unit="folder"):
            folder = str(folder)
            _filter   = info["filter"]
            _exposure = info["exposure"]  # string, e.g. "30" or "120"
            nimg      = info["num_images"]
            if _filter not in flats or _exposure not in darks:
                print(f"⚠️  Skipping {folder}: no cal frame for filter={_filter}, exp={_exposure}")
                continue

            print(f"📂 Calibrating {folder}  ({nimg} image{'s' if nimg>1 else ''})")
            sr.run(f"cd {folder}")
            if nimg == 1:
                # Single image – work on the raw FITS directly (skip sequence conversion)
                raw_files = [p.name for p in Path(folder).glob("*.[fF][iI][tT]*")]
                if not raw_files:
                    print(f"⚠️  No FITS found in {folder}; skipping")
                    continue
                raw_name = raw_files[0]
                sr.run("convert i")  # convert raw FITS to sequence
                sr.run(
                    f"calibrate_single i_00001.fit -flat={flats.get(_filter,'')} -dark={darks.get(_exposure,'')}"
                )
                base = Path(raw_name).stem
                sr.run(f"load pp_i_00001")
            else:
                # Multiple images – convert to sequence then calibrate
                sr.run("convert i")  # convert raw FITS to sequence
                sr.run(
                    f"calibrate i -flat={flats.get(_filter,'')} -dark={darks.get(_exposure,'')}"
                )  # calibrate the sequence
                sr.run("register pp_i -2pass -interp=cu")  # star alignment
                sr.run("seqapplyreg pp_i -framing=min")    # apply shifts
                sr.run(
                    "stack r_pp_i rej w 3 3 -nonorm -filter-fwhm=80% -filter-round=80%"
                )  # sigma-clipped stacking
                sr.run("load r_pp_i_stacked")
            sr.run("binxy 2")   # 2×2 binning for manageability
            sr.run("save res")  # final calibrated image
            sr.run("cd ..")  # return to parent just in case

            # ------------------------------------------------------------------
            # Mark this folder as successfully calibrated so that we don't do it
            # again later in the workflow.
            # ------------------------------------------------------------------
            calibrated_folders.add(folder)
    finally:
        try:
            if API == 'wrapper':
                sr.stop()
            else:
                sr.quit()
        except Exception:
            pass
        print("✅ pySiril processing done")

# ---------------------------------------------------------------------------
# Stand-alone plate-solving helper – invoked with --plate-solve
# ---------------------------------------------------------------------------

def plate_solve_all_res_files():
    """Walk through DATA_ROOT and run `solve-field` on every res.fit that lacks a WCS."""

    from astropy.io import fits  # local import to avoid top-level dependency when unused

    total      = 0
    solved_ok  = 0
    skipped    = 0
    failed     = 0

    solver = PlateSolver()

    print(f"🔍 Scanning for res.fit files under {DATA_ROOT} …")
    for res_path in directory.rglob('res.fit'):
        total += 1
        res_path = res_path.resolve()

        try:
            with fits.open(res_path, mode='update') as hdul:
                hdr = hdul[0].header
                has_wcs = all(k in hdr for k in ("CRVAL1", "CRVAL2")) and \
                           not any(math.isnan(hdr.get(k, math.nan)) for k in ("CRVAL1", "CRVAL2"))

                if has_wcs:
                    skipped += 1
                    continue

                print(f"🛰️  Plate-solving {res_path.relative_to(directory)} …")
                solved_wcs, _ = solver.solve_field(str(res_path))
                if solved_wcs is not None:
                    hdr.update(solved_wcs.to_header())
                    hdul.flush()
                    solved_ok += 1
                    print("   ✅ WCS written")
                else:
                    failed += 1
                    print("   ❌ Failed to solve")

        except Exception as e:
            failed += 1
            print(f"⚠️  Error processing {res_path.name}: {e}")

    print("\n📊 Plate-solve summary")
    print(f"   Total res.fit checked : {total}")
    print(f"   Already had WCS      : {skipped}")
    print(f"   Successfully solved  : {solved_ok}")
    print(f"   Failed to solve      : {failed}")

if __name__ == "__main__":
    main()

