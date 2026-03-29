# Architecture

## Overview

```
Browser (UI)
    │
    ▼
┌──────────────────────────────────────────┐
│  Night Manager  (Node.js, port 3000)     │
│  server/server.js                        │
│                                          │
│  • REST API for UI                       │
│  • NINA proxy                            │
│  • OCS proxy + 10-min weather poller     │
│  • Sequence engine                       │
│  • SQLite R/W (nightmanager.db)          │
│  • Job dispatcher → processing service   │
└─────────────────┬────────────────────────┘
                  │  HTTP (port 5200)
                  ▼
┌──────────────────────────────────────────┐
│  Processing Service  (Python/Flask)      │
│  astrobatch/processing_service.py        │
│                                          │
│  • Single-worker job queue               │
│  • Siril calibration + stacking          │
│  • Plate solving (solve-field)           │
│  • STDWeb upload + polling               │
│  • SQLite R/W (nightmanager.db)          │
│  • Results archived to NAS output        │
└──────────────────────────────────────────┘

External systems:
  NINA (Windows PC)       → port configurable (HTTP)
  OCS (192.168.1.220)     → HTTP /index.txt, /weatherpage.htm
  NAS                     → /mnt/nas/input/pyl/astro/
  STDWeb                  → https://stdweb.… (API + scrape)
  Astrometry.net          → solve-field (local install)
  Siril                   → siril-cli (local install)
```

## Components

### Night Manager (`server/server.js`)

Node.js / Express application. Single process, runs on the observatory machine.

**Responsibilities:**
- Serves the web UI (`server/public/`)
- Proxies all NINA HTTP API calls
- Controls OCS roof/dome via HTTP
- Polls OCS weather every 10 minutes → stores in `ocs_history`
- Runs automated multi-target imaging sequences
- Manages the pipeline job queue (dispatches to Python service)
- Reads/writes shared SQLite database

**Database:** `nightmanager.db` (default: project root)

### Processing Service (`astrobatch/processing_service.py`)

Flask application, background worker thread.

**Responsibilities:**
- Receives jobs via `POST /process`
- Runs one job at a time (single-worker queue)
- Full pipeline: scan → copy → Siril → plate solve → STDWeb → archive
- Recovers stuck jobs on startup
- Polls STDWeb every 30s for pending tasks
- Reads/writes shared SQLite database

**Database:** same `nightmanager.db` (path via `NIGHTMANAGER_DB` env)

### Web UI (`server/public/`)

Static files served by Night Manager.

| File | Role |
|------|------|
| `index.html` | Shell, tab structure, all HTML |
| `app.js` | All interactivity, API calls, polling, chart rendering |
| `astronomy.js` | Alt/Az computation, sky dome canvas, altitude curves |
| `nightplan.js` | Night plan scheduling, 7Timer forecast, Gantt chart |
| `styles.css` | Styling |

## Data Flow

### Imaging night
```
1. Operator opens UI → connects devices (NINA + OCS)
2. Searches TNS for targets → adds to Night Plan queue
3. Starts sequence → Night Manager drives NINA:
   cool camera → slew → plate solve → guide → autofocus → capture
4. FITS files land on NAS input path (NINA saves there)
5. Night Manager patches FILTER/OBJECT headers on NAS
```

### Pipeline (next morning)
```
1. Operator opens Pipeline tab → picks date → triggers job
2. Night Manager calls POST /process on Python service
3. Python service:
   - Scans NAS input FITS, groups by (target, filter, exposure)
   - Copies to work dir on USB disk (DATA_DIR)
   - Runs Siril: convert → calibrate → register → stack → res.fit
   - Runs astrometry.net: adds WCS to res.fit
   - Uploads res.fit to STDWeb
   - STDWeb: inspect → photometry → template subtraction
   - Archives res.fit + preview PNG to NAS output
4. Operator views results in Pipeline tab → photometry magnitudes
```
