# Configuration

## Night Manager (Node.js)

Set via environment variables before starting `node server.js`.

| Variable | Default | Description |
|----------|---------|-------------|
| `PORT` | `3000` | HTTP port for the web UI and API |
| `NINA_HOST` | `localhost` | Default NINA host (can be overridden per-request) |
| `NINA_PORT` | `1888` | Default NINA port |
| `NINA_PROTOCOL` | `http` | `http` or `https` |
| `OCS_HOST` | `192.168.1.220` | OCS observatory controller IP |
| `NAS_WATCH_PATH` | `/mnt/nas/input/pyl/astro/input` | Root of NAS input FITS tree |
| `NAS_OUTPUT` | `/mnt/nas/input/pyl/astro/output` | Root of NAS output (archived results) |
| `DB_FILE` | `../nightmanager.db` (relative to `server/`) | SQLite database path |
| `PROCESSING_SERVICE_URL` | `http://127.0.0.1:5200` | Base URL of Python processing service |
| `PROC_URL` | `http://127.0.0.1:5200` | Same, used for resume endpoint |
| `TNS_API_KEY` | — | TNS bot API key (optional; can be set per-request from UI) |
| `TNS_BOT_ID` | — | TNS bot ID |
| `TNS_BOT_NAME` | — | TNS bot name |
| `DEBUG` | — | Set to `1` to log OCS polling errors |

---

## Processing Service (Python)

Set via environment variables before starting `python -m astrobatch.processing_service`.

| Variable | Default | Description |
|----------|---------|-------------|
| `NIGHTMANAGER_DB` | `../nightmanager.db` | SQLite database path — **must match Night Manager's DB** |
| `DATA_DIR` | `/mnt/usb/astrobatch_data` | Work directory for Siril (use fast local disk, not NAS) |
| `NAS_OUTPUT` | `/mnt/nas/input/pyl/astro/output` | Where to archive results |
| `STDWEB_URL` | — | STDWeb base URL (e.g. `https://stdweb.example.com`) |
| `STDWEB_TOKEN` | — | STDWeb API token |
| `STDWEB_GAIN` | — | Detector gain in e-/ADU passed to STDWeb |
| `STDWEB_TEMPLATE` | `ztf` | Template source for subtraction (`ztf`, `ps1`, etc.) |
| `SIRIL_BIN` | `siril-cli` | Path to Siril binary |
| `DISPLAY` | `:0` | Required for Siril on Linux |
| `PROCESSING_PORT` | `5200` | Flask port |

---

## Critical: DB Path Must Match

Both services share **one SQLite database**. If they point to different files, jobs will appear stuck forever (the classic bug).

Recommended setup:
```
DB path (Night Manager):    /home/pyl/Documents/astrobatch/nightmanager.db
NIGHTMANAGER_DB (Python):   /home/pyl/Documents/astrobatch/nightmanager.db
```

---

## Systemd Services

Two service unit files are included in the repo root.

### `nina-control.service` — Night Manager

```ini
[Unit]
Description=AstroBatch Night Manager
After=network.target

[Service]
WorkingDirectory=/home/pyl/Documents/astrobatch/server
ExecStart=/usr/bin/node server.js
Restart=on-failure
Environment=PORT=3000
Environment=NAS_WATCH_PATH=/mnt/nas/input/pyl/astro/input
Environment=NAS_OUTPUT=/mnt/nas/input/pyl/astro/output
Environment=DB_FILE=/home/pyl/Documents/astrobatch/nightmanager.db
Environment=OCS_HOST=192.168.1.220

[Install]
WantedBy=multi-user.target
```

### `astrobatch-processing.service` — Processing Service

```ini
[Unit]
Description=AstroBatch Processing Service
After=network.target

[Service]
WorkingDirectory=/home/pyl/Documents/astrobatch
ExecStart=/home/pyl/Documents/astrobatch/.venv/bin/python -m astrobatch.processing_service
Restart=on-failure
Environment=NIGHTMANAGER_DB=/home/pyl/Documents/astrobatch/nightmanager.db
Environment=DATA_DIR=/mnt/usb/astrobatch_data
Environment=DISPLAY=:0
Environment=STDWEB_URL=https://...
Environment=STDWEB_TOKEN=...
Environment=STDWEB_GAIN=...

[Install]
WantedBy=multi-user.target
```

Install with:
```bash
sudo cp nina-control.service /etc/systemd/system/
sudo cp astrobatch-processing.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable nina-control astrobatch-processing
sudo systemctl start nina-control astrobatch-processing
```

---

## Calibration Frames

Flat and dark frames are auto-discovered from the `calib/` directory at the project root.

**Naming convention:**
```
calib/
  dark_30s.fit      # dark for 30-second exposures
  dark_120s.fit     # dark for 120-second exposures
  flat_G.fit        # flat for G filter
  flat_RP.fit       # flat for RP filter
  flat_BP.fit       # flat for BP filter
```

File extensions `.fit` and `.fits` are both accepted. If no calibration frame is found for a group, that step is skipped (Siril processes without calibration).

---

## NAS Directory Layout

```
NAS_WATCH_PATH/
  2026-03-28/
    SNAPSHOT/
      AT2026gze/          ← target folder
        AT2026gze_G_30s_001.fits
        AT2026gze_G_30s_002.fits
        AT2026gze_RP_120s_001.fits
      SN2026fvx/
        ...
    (other subdirs)

NAS_OUTPUT/
  2026-03-28/
    at_2026gze/
      G/
        res.fit
        res_preview.png
      RP/
        res.fit
        res_preview.png
```

NINA is configured to save FITS directly to `NAS_WATCH_PATH/<date>/SNAPSHOT/<target>/`. The Night Manager patches `FILTER` and `OBJECT` headers on NAS after each capture and fixes filename tokens.

---

## Time Accuracy

PC system time must be within **100ms of UTC** for accurate OCS polling timestamps. Use NTP:

```bash
sudo timedatectl set-ntp true
timedatectl status   # verify NTP synchronized: yes
```
