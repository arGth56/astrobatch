# Processing Pipeline

The pipeline takes raw FITS files from NAS and produces photometry measurements via STDWeb. It is implemented in `astrobatch/processing_service.py` and runs as a Flask service on port 5200.

## Overview

```
NAS input FITS
    │
    ▼  1. Scan & group
    │     Group by (OBJECT, FILTER, EXPTIME)
    ▼  2. Copy to work dir
    │     /mnt/usb/astrobatch_data/<date>/<target>/<filter>/<exp>s/
    ▼  3. Siril calibration + stacking
    │     convert → calibrate (darks/flats) → register → stack → res.fit
    ▼  4. Plate solve
    │     solve-field (astrometry.net) → WCS into res.fit
    ▼  5. STDWeb upload
    │     POST multipart res.fit with gain, target name, blind match hint
    ▼  6. STDWeb processing
    │     inspect → photometry → template subtraction
    ▼  7. Archive
    │     res.fit + res_preview.png → NAS output
    ▼
Magnitudes stored in nightmanager.db
```

## Step Details

### 1. Scan

The service reads the FITS directory and groups files by:
- `OBJECT` header (slugified)
- `FILTER` header
- `EXPTIME` header (rounded to nearest second)

For **SNAPSHOT** directories, each subdirectory is treated as one target group. A `selected_files` override can restrict which files are processed.

Special case: if only one frame is found for a group, it is processed as a single (no stacking).

### 2. Copy to Work Directory

Files are copied in parallel to:
```
DATA_DIR/<date>/<target>/<filter>/<exptime>s/
```
`DATA_DIR` defaults to `/mnt/usb/astrobatch_data/` (set via env var). Files already present (by name) are skipped.

### 3. Siril Calibration + Stacking

Siril is invoked via CLI (`siril-cli` or `siril -s`).

**Calibration frames** are auto-discovered from:
```
calib/
  dark_<exptime>s.fit (or .fits)
  flat_<filter>.fit
```

**Script sequence:**
```
convert i_      # convert to Siril internal format
calibrate i_    # apply dark/flat (if available)
register pp_i_  # star alignment
stack r_pp_i_   # Kappa-Sigma stack → res.fit
```

**Single frame:** `calibrate_single` is used instead of register+stack.

**Retry on failure:** if stacking fails due to "no stars", drops one frame and retries (up to 3 times).

**Output:** `res.fit` in the work directory.

### 4. Plate Solve

`PlateSolver` from `astrobatch/spliter.py` calls `solve-field` (astrometry.net local install):

```
solve-field --scale-units arcsecperpix \
            --scale-low 0.8 --scale-high 1.3 \
            --ra <CRVAL1> --dec <CRVAL2> \
            --radius 5 \
            res.fit
```

Scale range `0.8–1.3 arcsec/px` matches the primary camera setup (~1.044 arcsec/px).

On success, WCS is written into `res.fit` headers. Failure is non-fatal — upload proceeds without WCS.

### 5. STDWeb Upload

`res.fit` is uploaded to STDWeb via multipart POST:

```
POST /api/tasks/upload/
  file: res.fit
  gain: <STDWEB_GAIN>
  do_inspect: true
  do_photometry: true
  target_name: <OBJECT>
  blind_match_center_ra: <CRVAL1>   (if WCS present)
  blind_match_center_dec: <CRVAL2>
```

Returns a `task_id` stored in `pipeline_results.stdweb_task_id`.

### 6. STDWeb Processing

STDWeb processes the image through sequential steps. The service polls task state every ~30s:

| STDWeb state | Internal status |
|---|---|
| `uploaded` | `uploaded` |
| `inspect_running` | `inspecting` |
| `inspect_done` | `photometry` |
| `photometry_done` | `subtraction` |
| `subtraction_done` | `done` |
| `failed` / `error` | `error` |

**Inspect** — source detection, background estimation, PSF fitting, cosmic ray masking.

**Photometry** — aperture photometry on all detected sources; primary target matched by coordinates.

**Template subtraction** — subtracts a reference image (default: ZTF DR7 template) to isolate new sources.

### 7. Archive

On success:
```
res.fit        → NAS_OUTPUT/<date>/<target>/<filter>/res.fit
res_preview.png → NAS_OUTPUT/<date>/<target>/<filter>/res_preview.png
```
Siril intermediate files (work directory) are deleted to free USB disk space.

### Photometry Storage

After `subtraction_done`, the STDWeb poller scrapes the task for:
- `mjd` — Modified Julian Date
- `mag_ap`, `magerr_ap` — aperture photometry
- `mag_sub`, `magerr_sub` — subtraction photometry
- `mag_sub_ul` — upper limit flag (boolean)

Stored in `pipeline_results` columns; visible in the Pipeline and History tabs.

---

## Job Status Flow

```
queued → scanning → calibrating → stacking → solving →
uploading → inspecting → photometry → subtraction → done
                                                  ↘ error
```

Each status maps to a pipeline step. A job can be resumed from any step using `POST /resume`.

---

## Resume Steps

If a job fails mid-pipeline, it can be resumed from any step:

| `from_step` value | Resumes from |
|---|---|
| `calibrate` | Siril calibration |
| `solve` | Plate solving only |
| `upload` | STDWeb upload |
| `inspect` | Trigger STDWeb inspect |
| `photometry` | Trigger STDWeb photometry |
| `subtraction` | Trigger STDWeb subtraction |

---

## Job Queue

The processing service runs one job at a time.

- Jobs are enqueued by `POST /process` from Night Manager
- **Deduplication**: if a `job_id` is already active or queued, the request returns `409 Conflict` (idempotent)
- **Startup recovery**: on service start, all `queued` or `running` jobs in the DB are automatically re-enqueued (non-done results are deleted first for a clean retry)
- **Watchdog**: a background thread restarts the worker if it dies

---

## Per-group Result Cleanup

When a job is re-queued (restart), only results with `status != 'done'` are deleted. Groups that already completed successfully are **skipped** on re-run — their `done` results survive.

---

## Disk Layout

```
DATA_DIR/                          # USB disk work area
  2026-03-28/
    at_2026gze/
      G/
        30s/
          i_00001.fit              # Siril converted
          pp_i_00001.fit           # calibrated
          r_pp_i_stack.fit         # registered+stacked
          res.fit                  # final stacked

NAS_INPUT/                         # NAS input (read-only)
  2026-03-28/
    SNAPSHOT/
      AT2026gze/
        AT2026gze_G_30s_001.fits

NAS_OUTPUT/                        # NAS output (archived results)
  2026-03-28/
    at_2026gze/
      G/
        res.fit
        res_preview.png
```
