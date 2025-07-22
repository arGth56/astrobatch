"""astrobatch.analyse – temporary wrapper around legacy analyse.py

This keeps the external interface unchanged while we progressively
move code out of the monolithic top-level script.  All symbols from
`analyse.py` are re-exported so existing imports continue to work.
"""
from importlib import import_module as _imp
from pathlib import Path as _Path

# ---------------------------------------------------------------------------
# Redirect model & database files into the package directory so users do not
# clutter the project root.  We patch the legacy module *after* import so that
# any subsequent access picks up the new locations.
# ---------------------------------------------------------------------------

_PKG_DIR = _Path(__file__).resolve().parent
_MODEL_DIR = _PKG_DIR / "models"
_MODEL_DIR.mkdir(exist_ok=True)

_legacy = _imp("astrobatch.legacy_analyse")

globals().update(_legacy.__dict__)

# Redirect model file constants inside the legacy module
try:
    _legacy.MODEL_PATH = _MODEL_DIR / "candidate_model.json"
    _legacy.MODEL_PKL = _MODEL_DIR / "candidate_model.pkl"
except Exception:
    pass

# Fallback: if analyse.py defines a CLI default for --db, rewrite it to point
# inside the package (this is done by mutating argparse defaults after parser
# creation – the legacy script sets it via parser.add_argument(... default=...))
try:
    # The parser is created inside main(); we cannot access it easily. Instead
    # we override os.environ so that astrobatch.server and callers pick up the
    # relocated file.
    import os as _os
    _os.environ.setdefault("ASTROBATCH_CAND_DB", str(_PKG_DIR / "candidates.sqlite"))
except Exception:
    pass

# Provide a thin main() wrapper that injects --db <path> when user did not
# specify one.  This keeps CLI identical while transparently using the
# relocated database.

def main(argv=None):  # type: ignore
    import sys as _sys
    argv_in = list(_sys.argv[1:] if argv is None else argv)
    if "--db" not in argv_in and not any(arg.startswith("--db=") for arg in argv_in):
        # Prefer .sqlite, fall back to .db if it exists
        db_path = _PKG_DIR / "candidates.sqlite"
        if not db_path.exists():
            alt = _PKG_DIR / "candidates.db"
            if alt.exists():
                db_path = alt
        argv_in.extend(["--db", str(db_path)])

    # ------------------------------------------------------------------
    # Default locations for CSV/HTML reports (Output/)
    # ------------------------------------------------------------------
    from datetime import datetime as _dt
    _out_dir = (_Path.cwd() / "Output")
    _out_dir.mkdir(exist_ok=True)

    def _has_flag(flag: str) -> bool:
        return any(arg == flag or arg.startswith(f"{flag}=") for arg in argv_in)

    stamp = _dt.utcnow().strftime("%Y%m%d_%H%M%S")
    if not _has_flag("--csv"):
        argv_in.extend(["--csv", str(_out_dir / f"candidates_{stamp}.csv")])
    if not _has_flag("--html"):
        argv_in.extend(["--html", str(_out_dir / f"report_{stamp}.html")])

    _sys.argv = ["analyse"] + argv_in
    return _legacy.main()

__all__ = [k for k in globals().keys() if not k.startswith('_')] 