# API Reference

## Night Manager (Node.js, port 3000)

All endpoints return JSON. POST/PATCH endpoints expect `Content-Type: application/json`.

---

### Config

#### `GET /api/config/defaults`
Returns default configuration: NINA defaults, OCS host, whether TNS credentials are configured server-side.

---

### NINA

All NINA endpoints accept `{ ninaHost, ninaPort, ninaProtocol }` in the body (fallback to server defaults).

#### `POST /api/nina/test`
Test connection to NINA. Returns device map and equipment info.

#### `POST /api/nina/devices/status`
Refresh all device statuses from NINA.

#### `POST /api/nina/devices/connect`
Connect equipment.
```json
{ "device": "all" }          // connect everything
{ "device": "camera" }       // connect one device
// device: all | mount | camera | filterwheel | focuser
```

#### `POST /api/nina/actions/mount`
```json
{ "command": "park" }        // park | unpark | home
```

#### `POST /api/nina/actions/camera/capture`
```json
{
  "duration": 30,
  "gain": 100,
  "savePath": "/mnt/nas/...",   // optional
  "targetName": "AT2026gze"     // optional, for FITS header
}
```

#### `POST /api/nina/actions/filterwheel/change`
```json
{ "filterId": 2 }
```

#### `POST /api/nina/actions/focuser/move-relative`
```json
{ "direction": "in", "steps": 500 }
// direction: in | out
```

---

### Sequence

#### `GET /api/sequence/state`
Returns full sequence state:
```json
{
  "running": false,
  "currentTarget": null,
  "queue": [...],
  "log": [...],
  "manualMode": false,
  "lastAFTime": null
}
```

#### `POST /api/sequence/queue`
Add a target to the queue.
```json
{
  "name": "AT 2026gze",
  "raDeg": 131.524,
  "decDeg": 10.795,
  "ra": "08h46m05s",     // optional human-readable
  "dec": "+10°47'42\""
}
```

#### `DELETE /api/sequence/queue/:index`
Remove target at position `index`.

#### `POST /api/sequence/run`
Start the sequence.
```json
{
  "ninaHost": "192.168.1.100",
  "duration": 120,
  "gain": 10,
  "count": 10,
  "filters": ["G", "RP", "BP"],
  "solveEnabled": true,
  "manualMode": false
}
```

#### `POST /api/sequence/abort`
Gracefully abort after the current frame.

#### `POST /api/sequence/next`
In manual mode: confirm and proceed to next step.

#### `POST /api/sequence/manual-mode`
```json
{ "enabled": true }
```

#### `POST /api/sequence/clear`
Empty the queue.

#### `POST /api/sequence/reset-af`
Force autofocus on the next target.

#### `POST /api/sequence/restart`
Reset all targets to pending and restart from the beginning.

#### `POST /api/sequence/reorder`
Reorder the queue.
```json
{ "order": ["AT 2026gze", "SN 2026fvx", "AT 2026hbm"] }
```

---

### Targets

#### `GET /api/targets`
List all saved targets.

#### `POST /api/targets`
Save a target.
```json
{
  "name": "AT 2026gze",
  "ra": "08h46m05s",
  "dec": "+10°47'42\"",
  "ra_deg": 131.524,
  "dec_deg": 10.795,
  "type": "SN",
  "source": "tns",
  "notes": ""
}
```

#### `DELETE /api/targets/:id`
Delete a saved target.

#### `POST /api/target/tns`
Look up a target on TNS.
```json
{
  "name": "2026gze",
  "tnsApiKey": "...",
  "tnsBotId": "...",
  "tnsBotName": "..."
}
```
Returns TNS object data including RA/Dec, type, classification, redshift.

#### `GET /api/target/astrocolibri?name=AT2026gze`
Fetch AstroColibri event data for a target.

---

### OCS

#### `POST /api/ocs/status`
Get current OCS status (roof, weather, safety).
```json
{ "ocsHost": "192.168.1.220" }
```

#### `POST /api/ocs/roof`
Control the roof.
```json
{ "ocsHost": "192.168.1.220", "command": "open" }
// command: open | close | stop
```

#### `GET /api/ocs/history?hours=48`
Get stored weather history. `hours` can be 1–168 (max 7 days).

Returns:
```json
{
  "success": true,
  "history": [
    {
      "id": 42,
      "ts": "2026-03-28 08:30:00",
      "roof": "Closed",
      "safe": "UNSAFE",
      "rain": "Rain",
      "temp": 10.6,
      "humidity": 65.8,
      "pressure": 1033.0,
      "sky": 5.4
    }
  ],
  "hours": 48
}
```

#### `POST /api/ocs/poll-now`
Force an immediate OCS reading and store it.
```json
{ "ocsHost": "192.168.1.220" }
```

#### `POST /api/ocs/backfill`
Scrape OCS embedded chart history (recent 60 min + last 24h + last 48h) and insert missing rows.
```json
{ "ocsHost": "192.168.1.220" }
```
Returns:
```json
{ "success": true, "inserted": 48, "total": 73 }
```

---

### Pipeline

#### `GET /api/pipeline/jobs`
List all jobs with nested results.
```json
{
  "jobs": [
    {
      "id": 99,
      "target": "AT 2026gze",
      "status": "done",
      "fits_dir": "/mnt/nas/.../SNAPSHOT/AT2026gze",
      "created_at": "2026-03-28T...",
      "results": [
        {
          "id": 76,
          "filter": "G",
          "exposure": 30,
          "n_frames": 40,
          "status": "done",
          "stdweb_task_id": "3548",
          "mag_ap": 18.23,
          "magerr_ap": 0.05,
          "mag_sub": 18.31,
          "magerr_sub": 0.07
        }
      ]
    }
  ]
}
```

#### `DELETE /api/pipeline/jobs/:id`
Delete a job and all its results.

#### `DELETE /api/pipeline/results/:id`
Delete a single result row.

#### `POST /api/pipeline/trigger`
Start a new processing job.
```json
{
  "fits_dir": "/mnt/nas/input/pyl/astro/input/2026-03-28/SNAPSHOT/AT2026gze",
  "target": "AT 2026gze",
  "selected_files": ["file1.fits", "file2.fits"]  // optional
}
```

#### `POST /api/pipeline/rerun/:id`
Reset job `id` and reprocess from scratch.

#### `POST /api/pipeline/results/:id/retry-step`
Resume a result from a specific step.
```json
{ "from_step": "upload" }
// from_step: calibrate | solve | upload | inspect | photometry | subtraction
```

#### `GET /api/pipeline/nas-dates`
List available dates on NAS input path.

#### `GET /api/pipeline/nas-dates/:date`
List target folders for a given date (e.g. `2026-03-28`).

#### `GET /api/pipeline/job/:id/log`
Get the processing log for job `id`. Proxied from Python service.

#### `GET /api/pipeline/history`
Search completed results.

Query params: `target`, `filter`, `date_from`, `date_to`.

#### `GET /api/stdweb/task/:task_id/photometry`
Get photometry for a STDWeb task. Returns cached DB values or live scrape.

---

## Processing Service (Python/Flask, port 5200)

#### `GET /health`
```json
{
  "status": "ok",
  "worker_alive": true,
  "current_job": 99,
  "active_jobs": [99],
  "queue_depth": 2,
  "db_ok": true
}
```

#### `POST /process`
Enqueue a processing job.
```json
{
  "job_id": 99,
  "fits_dir": "/mnt/nas/...",
  "target": "AT 2026gze",
  "selected_files": [...],     // optional
  "object_filter": "AT2026gze" // optional, for SNAPSHOT dirs
}
```
Returns `409` if `job_id` is already active.

#### `POST /resume`
Resume a result from a specific step.
```json
{
  "result_id": 76,
  "from_step": "upload"
}
```

#### `GET /jobs/<job_id>/log`
Stream the processing log file for job `job_id`.
