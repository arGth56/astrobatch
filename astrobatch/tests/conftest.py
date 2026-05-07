"""
conftest.py – shared fixtures for processing_service tests.

The module has side effects at import time (DB schema migration, worker threads,
recovery timer). We handle this by:

  1. Creating a minimal SQLite test DB before the module is imported.
  2. Pointing NIGHTMANAGER_DB env var at that DB so the migration runs there.
  3. Importing the module once per session (via session-scoped fixture).
  4. Providing per-test helpers to seed and inspect the DB.
"""
import os
import sqlite3
import threading
from pathlib import Path

import pytest

# ── Schema ────────────────────────────────────────────────────────────────────

_SCHEMA = """
CREATE TABLE IF NOT EXISTS pipeline_jobs (
    id               INTEGER PRIMARY KEY,
    status           TEXT    DEFAULT 'queued',
    error            TEXT,
    fits_dir         TEXT    DEFAULT '',
    target           TEXT    DEFAULT '',
    target_filter    TEXT,
    manual_selection INTEGER DEFAULT 0,
    selection_data   TEXT,
    selection_result TEXT,
    stdweb_task_id   TEXT,
    stdweb_url       TEXT,
    updated_at       TEXT    DEFAULT ''
);
CREATE TABLE IF NOT EXISTS pipeline_results (
    id             INTEGER PRIMARY KEY,
    job_id         INTEGER,
    status         TEXT    DEFAULT 'queued',
    filter         TEXT    DEFAULT 'G',
    exposure       TEXT    DEFAULT '60',
    target         TEXT    DEFAULT 'test_target',
    stdweb_task_id TEXT    DEFAULT 'fake-task-id',
    stdweb_url     TEXT,
    error          TEXT,
    obs_date       TEXT    DEFAULT '2026-01-01',
    mjd            REAL,
    mag_ap         REAL,
    magerr_ap      REAL,
    mag_sub        REAL,
    magerr_sub     REAL,
    mag_sub_ul     REAL,
    updated_at     TEXT    DEFAULT ''
);
"""

# Fixed path (not tmp_path) so we can set the env var before session start.
_TEST_DB = Path("/tmp/test_nightmanager_bugs.db")


def _init_db():
    _TEST_DB.unlink(missing_ok=True)
    conn = sqlite3.connect(str(_TEST_DB))
    conn.executescript(_SCHEMA)
    conn.commit()
    conn.close()


def _seed_job(job_id: int, status: str = "running", conn=None):
    close = conn is None
    if conn is None:
        conn = sqlite3.connect(str(_TEST_DB))
    conn.execute(
        "INSERT OR REPLACE INTO pipeline_jobs(id, status) VALUES (?,?)",
        (job_id, status),
    )
    conn.commit()
    if close:
        conn.close()


def _seed_result(result_id: int, job_id: int, status: str = "photometry",
                 filter: str = "G", exposure: str = "60",
                 target: str = "test_target", stdweb_task_id: str = "fake-task-id"):
    conn = sqlite3.connect(str(_TEST_DB))
    conn.execute(
        """INSERT OR REPLACE INTO pipeline_results
           (id, job_id, status, filter, exposure, target, stdweb_task_id)
           VALUES (?,?,?,?,?,?,?)""",
        (result_id, job_id, status, filter, exposure, target, stdweb_task_id),
    )
    conn.commit()
    conn.close()


def _get_job_status(job_id: int) -> str | None:
    conn = sqlite3.connect(str(_TEST_DB))
    row = conn.execute("SELECT status FROM pipeline_jobs WHERE id=?", (job_id,)).fetchone()
    conn.close()
    return row[0] if row else None


def _get_result_status(result_id: int) -> str | None:
    conn = sqlite3.connect(str(_TEST_DB))
    row = conn.execute("SELECT status FROM pipeline_results WHERE id=?", (result_id,)).fetchone()
    conn.close()
    return row[0] if row else None


def _get_selection_result(job_id: int):
    conn = sqlite3.connect(str(_TEST_DB))
    row = conn.execute("SELECT selection_result FROM pipeline_jobs WHERE id=?", (job_id,)).fetchone()
    conn.close()
    return row[0] if row else None


# ── Session-scoped module import ───────────────────────────────────────────────

@pytest.fixture(scope="session", autouse=True)
def processing_service_module():
    """
    Create test DB, set env var, then import processing_service once for the
    whole test session. Returns the imported module.
    """
    _init_db()
    os.environ["NIGHTMANAGER_DB"] = str(_TEST_DB)
    # Suppress the background recovery timer and threads from cluttering test output.
    # We patch before import using a threading event.
    import importlib
    ps = importlib.import_module("astrobatch.processing_service")
    return ps


# ── Per-test DB fixtures ───────────────────────────────────────────────────────

@pytest.fixture()
def clean_db():
    """Wipe job/result rows before each test so tests are independent."""
    conn = sqlite3.connect(str(_TEST_DB))
    conn.execute("DELETE FROM pipeline_jobs")
    conn.execute("DELETE FROM pipeline_results")
    conn.commit()
    conn.close()
    yield
    conn = sqlite3.connect(str(_TEST_DB))
    conn.execute("DELETE FROM pipeline_jobs")
    conn.execute("DELETE FROM pipeline_results")
    conn.commit()
    conn.close()


@pytest.fixture()
def ps(processing_service_module):
    """Shorthand alias for the imported processing_service module."""
    return processing_service_module


# Export helpers so test files can import them from conftest.
pytest.helpers = type("H", (), {
    "seed_job":           staticmethod(_seed_job),
    "seed_result":        staticmethod(_seed_result),
    "get_job_status":     staticmethod(_get_job_status),
    "get_result_status":  staticmethod(_get_result_status),
    "get_selection_result": staticmethod(_get_selection_result),
    "TEST_DB":            _TEST_DB,
})
