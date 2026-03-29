# CLI Tools

The CLI provides standalone pipeline operations without the web UI, useful for batch reprocessing, testing, or running on a different machine.

## Entry Point

```bash
python -m astrobatch.cli [options]
```

Or via the `astrobatch` command if installed as a package.

## Commands

### Split & group FITS
```bash
python -m astrobatch.cli --split --input /path/to/fits/dir
```
Groups FITS files by `(OBJECT, FILTER, EXPTIME)` and creates a structured directory tree.

### Calibrate + Stack
```bash
python -m astrobatch.cli --calibrate --input /path/to/fits/dir
```
Runs Siril calibration and stacking on each group. Produces `res.fit` per group.

### Plate Solve
```bash
python -m astrobatch.cli --plate-solve --input /path/to/res.fit
```
Runs `solve-field` and writes WCS into the FITS file.

### Upload to STDWeb
```bash
python -m astrobatch.cli --upload --input /path/to/res.fit \
  --target "AT 2026gze"
```
Uploads `res.fit` to STDWeb. Requires `STDWEB_URL` and `STDWEB_TOKEN` env vars.

### Analyse (photometry)
```bash
python -m astrobatch.cli --analyse --task-id 3548
```
Polls STDWeb task and fetches photometry results.

### Mosaic
```bash
python -m astrobatch.cli --mosaic --input /path/to/results/dir
```
Combines multiple `res.fit` files into a mosaic. Output goes to `ROOT/mosaics/`.

### Full pipeline
```bash
python -m astrobatch.cli --split --calibrate --upload --analyse \
  --input /mnt/nas/input/.../SNAPSHOT/AT2026gze \
  --target "AT 2026gze"
```
Chains all steps in order.

### Start candidate labelling server
```bash
python -m astrobatch.cli --start-server
```
Starts a Flask app on ~port 5100 for labelling candidate detections as `real` / `bogus` / `unclear`.

## Environment Variables (CLI)

| Variable | Description |
|----------|-------------|
| `STDWEB_API_TOKEN` | STDWeb API token (note: different name from service's `STDWEB_TOKEN`) |
| `STDWEB_URL` | STDWeb base URL |
| `DATA_DIR` | Work directory for Siril |
| `SIRIL_BIN` | Path to `siril-cli` |
| `DISPLAY` | X display for Siril (Linux) |

## Key Modules

| File | Purpose |
|------|---------|
| `astrobatch/cli.py` | CLI entry point and argument parsing |
| `astrobatch/spliter.py` | Core engine: split, Siril scripts, upload, plate solve, TNS helpers, mosaic |
| `astrobatch/analyse.py` | STDWeb analysis wrapper |
| `astrobatch/db.py` | Candidate database (separate from nightmanager.db) |
| `astrobatch/server.py` | Candidate labelling Flask server |
| `astrobatch/build_master_flats.py` | Utility to build master flat frames |

## Candidate Labelling

After photometric analysis, candidates (potential transients) can be reviewed via a web interface:

```bash
python -m astrobatch.cli --start-server
# Open http://localhost:5100 in browser
```

Shows cutout images for each candidate. Operator labels each as `real`, `bogus`, or `unclear`. Labels are stored in the `candidates` SQLite table for later analysis.
