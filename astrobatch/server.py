"""Light-weight REST service exposing the *candidates* SQLite database.

Run via::

    python -m astrobatch.cli --start-server [--port 5100]

It mirrors the original `server.py` functionality but is namespaced
inside the *astrobatch* package and imports `CandidateDB` from
``astrobatch.db``.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Iterable

from flask import Flask, request, jsonify
from flask_cors import CORS

from .db import CandidateDB

DEFAULT_PORT = int(os.environ.get("PORT", 5100))
# Default DB path inside package
_PKG_DIR = Path(__file__).resolve().parent
_DB_PATH = Path(os.environ.get("ASTROBATCH_CAND_DB") or (_PKG_DIR / "candidates.sqlite"))

# Backwards-compat: if the .sqlite file does not exist but a legacy
# `candidates.db` does, use that instead.
if not _DB_PATH.exists():
    legacy_db = _PKG_DIR / "candidates.db"
    if legacy_db.exists():
        _DB_PATH = legacy_db

# Global objects ------------------------------------------------------------
db = CandidateDB(_DB_PATH)
app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})

# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.route("/flag", methods=["POST", "GET"])
def flag_candidate():
    """Mark a candidate row with a human label (tp/fp/unknown)."""
    try:
        task_id = int(request.values.get("task_id"))
    except (TypeError, ValueError):
        return jsonify(error="task_id must be integer"), 400
    cutout = request.values.get("cutout")
    label = request.values.get("label")
    if not cutout or not label:
        return jsonify(error="cutout and label required"), 400
    db.upsert({"task_id": task_id, "cutout": cutout})
    db.set_label(task_id, cutout, label)
    return jsonify(ok=True)


@app.route("/candidate", methods=["POST"])
def upsert_candidate():
    """Insert or update a single candidate row via JSON payload."""
    data: Dict[str, Any] | None = request.get_json(silent=True)
    if not data:
        return jsonify(error="JSON body required"), 400
    if not {"task_id", "cutout"}.issubset(data):
        return jsonify(error="task_id and cutout required"), 400
    db.upsert(data)
    return jsonify(ok=True)


@app.route("/candidates", methods=["GET"])
def list_candidates():
    rows = [dict(r) for r in db.fetch_all()]
    return jsonify(rows)


# ---------------------------------------------------------------------------
# Entrypoint callable from CLI
# ---------------------------------------------------------------------------

def start_server(port: int = DEFAULT_PORT) -> None:
    """Blocking call to launch the Flask development server."""
    app.run(host="0.0.0.0", port=port, debug=True, use_reloader=False) 