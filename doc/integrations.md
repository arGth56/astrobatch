# External Integrations

## NINA (Nighttime Imaging 'N' Astronomy)

NINA runs on a Windows PC and controls the imaging hardware. AstroBatch communicates with NINA via its built-in HTTP API.

**Base URL:** `http://<NINA_HOST>:<NINA_PORT>/v2/api/`

**Key endpoints used:**

| Endpoint | Usage |
|----------|-------|
| `GET /equipment/*` | Device connection status |
| `POST /equipment/camera/connect` | Connect camera |
| `POST /equipment/mount/slew` | Slew to RA/Dec |
| `POST /equipment/mount/tracking` | Enable/disable tracking |
| `POST /equipment/mount/park` | Park mount |
| `POST /equipment/mount/unpark` | Unpark |
| `POST /equipment/camera/capture` | Take an exposure |
| `POST /equipment/camera/cool` | Set cooling target temperature |
| `POST /equipment/camera/abort-exposure` | Stop current exposure |
| `POST /equipment/filterwheel/change-filter` | Select filter by ID |
| `POST /equipment/focuser/auto-focus` | Run autofocus routine |
| `POST /equipment/focuser/move` | Move focuser by steps |
| `POST /equipment/guider/start` | Start PHD2 guiding |
| `POST /equipment/guider/stop` | Stop guiding |

**Long exposures** use a custom polling helper (`callNinaLong`) that waits with a configurable timeout and retries on transient network errors.

---

## OCS (Observatory Control System)

The OCS is a small embedded controller (OCS 3.12k) at `192.168.1.220` that controls the roof/dome and reads weather sensors.

### Status polling
```
GET http://192.168.1.220/index.txt
```
Returns `key|value` lines:
```
roof_sta|Closed
stat_safe|UNSAFE
wea_rain|Rain
wea_temp| 10.6 °C
wea_humd|  65.8 %
wea_pres|  1033 mb
wea_sq| 5.4 mpsas
```

### Roof control
```
GET http://192.168.1.220/index-ajax-get.txt?cmd=roof_open
GET http://192.168.1.220/index-ajax-get.txt?cmd=roof_close
GET http://192.168.1.220/index-ajax-get.txt?cmd=roof_stop
```

### History scraping (backfill)
The OCS stores 60-min, 24h, and 48h weather history in Chart.js data embedded in HTML pages:

```
GET /weatherpage.htm              → recent 60 min (temp, pressure, humidity)
GET /weatherpage.htm?chart=last24 → last 24 hours
GET /weatherpage.htm?chart=last48 → last 48 hours
GET /skypage.htm                  → recent 60 min (sky quality)
GET /skypage.htm?chart=last48     → last 48 hours
```

Data points are `{x: hoursAgo, y: value}` format, parsed by the backfill endpoint.

---

## STDWeb

STDWeb is a photometric pipeline web service that performs source extraction, aperture photometry, and image subtraction.

**Base URL:** configured via `STDWEB_URL` env var.
**Auth:** API token via `STDWEB_TOKEN` (sent as `Authorization: Token ...` header).

### Upload
```
POST /api/tasks/upload/
Content-Type: multipart/form-data

file: <res.fit binary>
gain: 1.0
do_inspect: true
do_photometry: true
target_name: AT 2026gze
blind_match_center_ra: 131.524   (optional, from WCS)
blind_match_center_dec: 10.795
```
Returns `{ id: "3548" }` — the task ID.

### State polling
```
GET /api/tasks/<task_id>/
```
Returns task state string. Key states:

| STDWeb state | Meaning |
|---|---|
| `uploaded` | File received |
| `inspect_running` | Source detection in progress |
| `inspect_done` | Inspection complete |
| `photometry_running` | Photometry in progress |
| `photometry_done` | Photometry complete |
| `subtraction_running` | Template subtraction in progress |
| `subtraction_done` | All done |
| `failed` / `error` | Pipeline error |

### Actions
```
POST /api/tasks/<task_id>/action/
{ "action": "inspect" }
{ "action": "photometry" }
{ "action": "subtract", "template": "ztf" }
```

### Photometry scraping
Magnitudes are scraped from the STDWeb task HTML page since the API doesn't expose them directly. The scraper extracts:
- Aperture magnitude and error
- Subtraction magnitude and error (or upper limit)
- MJD of observation

---

## Astrometry.net (solve-field)

Local installation of astrometry.net provides plate solving.

**Binary:** `solve-field` (path configurable via `SIRIL_BIN` or auto-detected).

**Call parameters for this setup:**
```bash
solve-field \
  --scale-units arcsecperpix \
  --scale-low 0.8 \
  --scale-high 1.3 \
  --ra <CRVAL1> --dec <CRVAL2> \
  --radius 5 \
  --no-plots \
  --overwrite \
  res.fit
```

Scale range `0.8–1.3 arcsec/px` is calibrated for the primary camera (~1.044 arcsec/px actual).

On success: WCS headers (`CRVAL1`, `CRVAL2`, `CD1_1`, etc.) are written into `res.fit`. If solving fails, the result is uploaded to STDWeb with blind matching enabled.

---

## Siril

Siril is the calibration and stacking engine.

**Binary:** `siril-cli` (or `siril -s` for GUI-less operation).
**Required:** `DISPLAY` env var set (e.g. `:0`) on Linux.

**Typical script sequence:**
```
requires 1.2.0
cd /work/dir/

convert i_       # load raw FITS as Siril sequence
calibrate i_ \
  -dark=dark_30s \
  -flat=flat_G \
  -cc=dark       # calibration

register pp_i_   # star alignment (multi-frame only)
stack r_pp_i_ \
  rej 3 3 \
  type=ks        # Kappa-Sigma stack

save res         # output as res.fit
```

**Single frame:** uses `calibrate_single` instead of register+stack.

**Calibration frames** are auto-discovered from the `calib/` directory (naming: `dark_<exptime>s.fit`, `flat_<filter>.fit`).

---

## TNS (Transient Name Server)

Used to look up transient target coordinates and classifications.

**API:** `https://www.wis-tns.org/api/get/object`

**Auth:** Bot credentials (API key, bot ID, bot name) — set as env vars or entered in the UI per-session.

**Usage:** Type a target name (e.g. `2026gze`), AstroBatch calls TNS and returns RA/Dec, type, discovery date, redshift, host galaxy, and classification reports.

If the TNS API returns no coordinates, AstroBatch falls back to scraping the TNS object HTML page.

---

## AstroColibri

Optional enrichment for TNS targets.

**API:** `https://astro-colibri.science/event?trigger_id=TNS<name>`

Returns additional context: multi-messenger alerts, GCN notices, related events.

---

## 7Timer (Weather Forecast)

Astronomical weather forecast in the Night Plan tab.

**API:** `http://www.7timer.info/bin/astro.php?lon=<lon>&lat=<lat>&ac=0&unit=metric&output=json`

Returns hourly forecast for: cloud cover, seeing, transparency, humidity, wind speed. Displayed as a Gantt-style chart overlay in the Night Plan tab.
