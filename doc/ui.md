# Web Dashboard

Access at `http://<observatory-machine>:3000`

## Connection Bar (top)

Always visible. Before doing anything, connect NINA and OCS.

| Field | Description |
|-------|-------------|
| NINA Host | IP of the Windows PC running NINA (default: `localhost`) |
| Port | NINA HTTP API port (default: `1888`) |
| Protocol | `http` or `https` |
| OCS Host | OCS observatory controller IP (default: `192.168.1.220`) |
| Connect | Tests NINA connection and fetches device map |

Status indicators show NINA connection state and OCS reachability.

---

## Tab: Devices

Shows the status of all connected equipment.

**NINA devices** (polled from NINA):
- Mount — connection state, tracking, parked/unparked
- Camera — connection, temperature, cooling state
- Filter Wheel — current filter
- Focuser — position, temperature

**OCS row** — roof status, safety flag, rain, temperature (live from OCS).

**Connect controls:**
- Connect All — connects every device at once
- Individual connect buttons per device

---

## Tab: Actions

Manual one-shot controls. Use for testing or manual operation.

### Mount
| Button | Effect |
|--------|--------|
| Park | Sends park command to NINA |
| Unpark | Unparks and enables tracking |
| Home | Finds home position |

### Camera
- **Duration** — exposure length (seconds)
- **Gain** — camera gain
- **Save path** — optional NAS path to save FITS
- **Capture** — triggers one exposure, shows result

### Filter Wheel
- Dropdown of available filters → **Change** applies it

### Focuser
- **Steps** — number of steps to move
- **In / Out** — relative move buttons

---

## Tab: Dome

Observatory weather and roof control.

### Roof controls
| Button | Effect |
|--------|--------|
| Open Roof | Sends open command to OCS |
| Close Roof | Sends close command |
| Stop | Emergency stop |

### Live weather (from OCS)
Current: temperature, humidity, pressure, rain state, safety flag, sky quality.

### 48h Weather Timeline (Chart.js)

Displays stored `ocs_history` as a timeline graph:

| Series | Color | Y axis |
|--------|-------|--------|
| Temperature (°C) | Amber | Left |
| Humidity (%) | Sky blue | Right (0–100) |
| Pressure (hPa) | Purple | Right (auto) |
| Sky Quality (mag/″²) | Cyan | Left (shared) |
| Roof open | Green tint | Background band |
| Rain | Indigo bar | Bottom 20% strip |

**Controls:**
- **Time range dropdown** — Last 12h / 24h / 48h / 72h
- **Refresh** — reload from DB, reset zoom
- **Poll Now** — force an immediate OCS reading and store it
- **Backfill 60 min** — scrapes OCS embedded charts (recent 60 min + last 24h + last 48h) and inserts missing rows into DB. Run once to seed history.

**Zoom:** Click and drag on the chart to zoom into a time range. Refresh to reset.

---

## Tab: Target

Search for transient targets and build the observation list.

### TNS Search
1. Enter TNS bot credentials (API key, bot ID, bot name) — stored in browser
2. Type target name (e.g. `2026gze` → looks up `AT 2026gze`)
3. Result shows: RA/Dec, type, discovery date, redshift, host galaxy
4. Optional: **AstroColibri enrichment** — fetches additional context

### Sky Visibility
- Observer coordinates (lat/lon/elevation) used to compute:
  - Current altitude/azimuth
  - **Sky dome** — 2D canvas showing target position
  - **Altitude curve** — rises/sets, twilight windows

### Target table
Shows all matched targets with coordinates, altitude, rise/set times.

**Actions per target:**
- **Add to Night Plan** — queues for automated sequence
- **Save** — saves to local `targets` table for reuse

---

## Tab: Night Plan

Automated multi-target sequence planner and runner.

### Sequence parameters
| Parameter | Description |
|-----------|-------------|
| Exposure | Seconds per frame |
| Gain | Camera gain |
| Frames per filter | Number of exposures per filter per target |
| Filters | Comma-separated list (e.g. `G,RP,BP`) |
| Plate-solve & center | Enable slew+solve loop before imaging |
| Per-frame check | Verify each frame is on-target after capture |

### Target queue
- Ordered list of targets to observe
- **Drag to reorder** — drag handles to change priority
- **Delete** — remove a target
- **Add Wait** — insert a timed pause between targets
- **Optimize order** — sorts by altitude/visibility window

### Altitude Gantt
Visual timeline showing when each target is above the horizon tonight. Helps spot scheduling conflicts.

### 7Timer forecast
Astronomical weather forecast for tonight at the observer location.

### Sequence controls
| Button | Action |
|--------|--------|
| Run | Starts the sequence from the first pending target |
| Abort | Gracefully stops after current frame |
| Restart | Resets all targets to pending and restarts |
| Clear | Empties the queue |
| Reset AF | Forces autofocus on next target regardless of timing |
| Next (manual mode) | Confirms next step when manual mode is on |

### Sequence log
Live log of all sequence steps, NINA responses, errors. Retained until next run.

---

## Tab: Pipeline

Trigger and monitor processing jobs.

### Trigger a job
1. **Pick date** — dropdown of dates found on NAS input path
2. **Pick target folder** — lists SNAPSHOT subdirectories for that date
3. **Trigger** — sends job to Python processing service

### Job table
| Column | Meaning |
|--------|---------|
| ID | Job number |
| Target | Object name |
| Status | `queued` → `scanning` → `calibrating` → `stacking` → `solving` → `uploading` → `inspecting` → `photometry` → `subtraction` → `done` / `error` |
| Groups | Number of (filter, exposure) groups found |
| Results | Per-group result rows with status and STDWeb links |

**Per-job actions:**
- **Log** — live log from the processing service
- **Rerun** — reset and reprocess
- **Delete** — remove job and results

**Per-result actions:**
- **View on STDWeb** — opens the task in STDWeb
- **Photometry** — shows stored magnitude measurements
- **Retry step** — resume from a specific pipeline step
- **Download** — download the archived `res.fit`

### Photometry modal
For each completed result, shows:
- Aperture magnitude (`mag_ap`) ± error
- Subtraction magnitude (`mag_sub`) ± error or upper limit
- MJD of observation
- Copy as CSV for reporting

---

## Tab: History

Search all completed pipeline results.

**Search filters:**
- Target name (partial match)
- Date range
- Filter

**Result table:**
- Target, filter, exposure, date, MJD, magnitudes
- Link to STDWeb task
- Copy all rows as CSV (for ATel / TNS reporting)
