# Automated Imaging Sequence

The Night Plan sequence engine is implemented in `server/server.js` (`runSequence` function). It drives NINA to image a queue of targets fully automatically.

## Prerequisites

Before starting the sequence:
1. All devices connected (Devices tab)
2. Targets added to Night Plan queue
3. Sequence parameters configured (exposure, gain, frames, filters)
4. Weather safe (OCS showing SAFE)

## Sequence Flow

### Startup
1. Camera cooling — cools to target temperature (~−5°C), polls until stable
2. Manual confirm (if manual mode enabled)

### Per Target Loop

For each pending target in queue order:

#### 1. Pre-checks
- Camera must be idle (not capturing)
- Target horizon check: altitude must be > 5° using observer lat/lon

#### 2. Unpark & tracking
- Unpark mount
- Enable tracking

#### 3. Slew to target
- Slew to target RA/Dec
- Wait for slew to complete

#### 4. Plate-solve & center *(if `solveEnabled`)*
- Take a short exposure
- Solve the field (via NINA)
- If offset > threshold: correct and re-slew
- Repeat until centered

#### 5. Start guiding
- Start PHD2 guiding via NINA
- Wait for guiding to settle

#### 6. Autofocus *(if > 1h since last AF or reset requested)*
- Trigger NINA autofocus routine
- Record autofocus time on success

#### 7. Capture loop — per filter

For each filter in the filter list:
1. Switch filter via NINA
2. For N frames:
   - Trigger exposure (gain, duration)
   - NINA saves FITS to NAS
   - Night Manager patches `FILTER` and `OBJECT` FITS headers on NAS
   - Fix NINA filename tokens if needed
3. Manual confirm between filters (if manual mode)

#### 8. Target complete
- Mark target as `done` in queue

### End of Queue
- Park mount (home position)

---

## Error Handling

### Per-frame network errors
If a NINA API call fails with a network error (e.g. brief disconnect):
- Retry the frame once after a 15-second wait
- If retry fails, log and skip to next frame

### Per-target errors
If a target fails (non-network error):
- Log the error
- Home mount
- Continue with next target
- Up to N retries per target before giving up

### Camera busy (409)
If NINA returns "camera busy":
- Wait 30 seconds and retry the capture
- Repeats up to 3 times

---

## Manual Mode

When manual mode is enabled, the sequence pauses at key steps and waits for the operator to click **Next** in the UI:

- After camera cooling
- After slewing to target (before guiding)
- Between filters
- After autofocus

Useful for supervised runs or when testing a new setup.

---

## Per-Frame Verification

*(Partially implemented)*

When enabled, each captured frame is plate-solved after capture. If the solved position differs from the target by more than a threshold (arcminutes), the sequence logs a warning. This detects mount drift, guiding failures, or incorrect pointing.

---

## Autofocus Logic

- Autofocus runs automatically on the first target of the night
- Subsequent targets skip autofocus if the last AF was < 1 hour ago
- Temperature changes can also trigger AF (if NINA's built-in trigger is enabled)
- **Reset AF** button forces AF on the next target regardless of timing

---

## Filter Order

Filters are imaged in the order specified in the sequence parameters (e.g. `G, RP, BP`).

All N frames of one filter are taken before moving to the next. This minimises filter wheel wear and avoids frequent refocusing.

---

## Sequence Parameters Reference

| Parameter | Description | Example |
|-----------|-------------|---------|
| `duration` | Seconds per exposure | `120` |
| `gain` | Camera gain | `10` |
| `count` | Frames per filter per target | `10` |
| `filters` | Filter list (comma-separated) | `G,RP,BP` |
| `solveEnabled` | Plate-solve before imaging | `true` |
| `manualMode` | Pause at each step | `false` |

---

## State Polling

The Night Plan UI polls `GET /api/sequence/state` every ~2 seconds while the tab is active. The state response includes:

```json
{
  "running": true,
  "currentTarget": "AT 2026gze",
  "currentFilter": "G",
  "currentFrame": 7,
  "totalFrames": 10,
  "queue": [...],
  "log": ["[22:01:15] Slewing to AT 2026gze", "..."],
  "manualMode": false,
  "lastAFTime": "2026-03-28T22:00:00.000Z"
}
```

The log is displayed live in the Night Plan tab and retained until the next sequence run.
