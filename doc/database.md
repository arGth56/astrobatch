# Database Schema

AstroBatch uses a single SQLite database (`nightmanager.db`) shared between the Night Manager and the Processing Service.

---

## `pipeline_jobs`

One row per processing job triggered from the Pipeline tab.

| Column | Type | Description |
|--------|------|-------------|
| `id` | INTEGER PK | Auto-increment job ID |
| `target` | TEXT | Target name (e.g. `AT 2026gze`) |
| `filter` | TEXT | Comma-separated filter list (legacy single-job field) |
| `exposure` | REAL | Exposure time in seconds |
| `fits_dir` | TEXT | NAS input directory that was processed |
| `status` | TEXT | Current status (see flow below) |
| `stdweb_task_id` | TEXT | STDWeb task ID (legacy top-level field) |
| `stdweb_url` | TEXT | STDWeb task URL |
| `error` | TEXT | Error message if failed |
| `created_at` | TEXT | ISO timestamp |
| `updated_at` | TEXT | ISO timestamp |
| `target_filter` | TEXT | Combined key (migration column) |

**Status flow:**
```
queued → scanning → calibrating → stacking → solving →
uploading → inspecting → photometry → subtraction → done | error
```

---

## `pipeline_results`

One row per (job, filter, exposure) group. A job with 3 filters produces 3 result rows.

| Column | Type | Description |
|--------|------|-------------|
| `id` | INTEGER PK | Auto-increment result ID |
| `job_id` | INTEGER FK | → `pipeline_jobs.id` |
| `target` | TEXT | Target name |
| `filter` | TEXT | Filter name (e.g. `G`) |
| `exposure` | REAL | Exposure in seconds |
| `n_frames` | INTEGER | Number of frames stacked |
| `status` | TEXT | Same status flow as jobs |
| `stdweb_task_id` | TEXT | STDWeb task ID for this result |
| `stdweb_url` | TEXT | STDWeb task URL |
| `stdweb_state` | TEXT | Raw STDWeb state string |
| `error` | TEXT | Error message if failed |
| `obs_date` | TEXT | Observation date (from FITS DATE-OBS) |
| `mjd` | REAL | Modified Julian Date |
| `mag_ap` | REAL | Aperture photometry magnitude |
| `magerr_ap` | REAL | Aperture magnitude uncertainty |
| `mag_sub` | REAL | Subtraction photometry magnitude |
| `magerr_sub` | REAL | Subtraction magnitude uncertainty |
| `mag_sub_ul` | INTEGER | 1 if `mag_sub` is an upper limit |
| `frame_previews` | TEXT | JSON array of preview image paths |
| `created_at` | TEXT | ISO timestamp |
| `updated_at` | TEXT | ISO timestamp |

---

## `seq_queue`

Imaging sequence queue — targets to observe tonight.

| Column | Type | Description |
|--------|------|-------------|
| `position` | INTEGER | Sort order (0-based) |
| `name` | TEXT | Target name |
| `ra` | TEXT | RA in HMS format |
| `dec` | TEXT | Dec in DMS format |
| `ra_deg` | REAL | RA in decimal degrees |
| `dec_deg` | REAL | Dec in decimal degrees |
| `done` | INTEGER | 1 when the target has been observed |
| `added_at` | TEXT | ISO timestamp |

---

## `targets`

Saved target library (persists across nights).

| Column | Type | Description |
|--------|------|-------------|
| `id` | INTEGER PK | Auto-increment |
| `name` | TEXT | Target name |
| `ra` | TEXT | RA in HMS |
| `dec` | TEXT | Dec in DMS |
| `ra_deg` | REAL | RA decimal degrees |
| `dec_deg` | REAL | Dec decimal degrees |
| `type` | TEXT | Object type (SN, AT, Nova, etc.) |
| `source` | TEXT | Where it came from (`tns`, `manual`) |
| `notes` | TEXT | Free text |
| `created_at` | TEXT | ISO timestamp |

---

## `ocs_history`

Weather and observatory condition history, stored every 10 minutes.

| Column | Type | Description |
|--------|------|-------------|
| `id` | INTEGER PK | Auto-increment |
| `ts` | TEXT | UTC timestamp `YYYY-MM-DD HH:MM:SS` — unique index after first backfill |
| `roof` | TEXT | Roof state (`Open`, `Closed`, `Opening`, `Closing`) |
| `safe` | TEXT | Safety state (`SAFE`, `UNSAFE`) |
| `rain` | TEXT | Rain state (`Rain`, `Dry`) |
| `temp` | REAL | Outdoor temperature (°C) |
| `humidity` | REAL | Relative humidity (%) |
| `pressure` | REAL | Barometric pressure (hPa) |
| `sky` | REAL | Sky quality (mag/arcsec²) |
| `ir_sky` | REAL | IR thermal camera average sky temperature (°C, from MLX90621) |

**Note:** Backfilled rows (from `POST /api/ocs/backfill`) will not have `roof`, `safe`, or `rain` values (OCS does not include those in its chart data).

---

## `candidates` (separate DB)

Used by the CLI candidate labelling server (`astrobatch/server.py`). Stored in a separate SQLite file managed by `astrobatch/db.py`.

| Column | Description |
|--------|-------------|
| `id` | Auto-increment |
| `target` | Target name |
| `filter` | Filter |
| `mjd` | MJD |
| `mag`, `magerr` | Photometry |
| `ra`, `dec` | Coordinates |
| `cutout_*` | Preview image paths |
| `label` | Human review label (`real`, `bogus`, `unclear`) |
| `created_at` | Timestamp |

---

## Querying Useful Data

```sql
-- Latest photometry per target
SELECT target, filter, mjd, mag_ap, magerr_ap, mag_sub, magerr_sub
FROM pipeline_results
WHERE status = 'done'
ORDER BY mjd DESC;

-- Jobs still stuck
SELECT id, target, status, updated_at
FROM pipeline_jobs
WHERE status NOT IN ('done', 'error')
ORDER BY created_at;

-- Weather last 24h
SELECT ts, roof, rain, temp, humidity, sky
FROM ocs_history
WHERE ts >= datetime('now', '-24 hours')
ORDER BY ts;
```
