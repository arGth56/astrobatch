"""Lightweight SQLite wrapper for candidate objects.

This is a direct migration of the legacy CandidateDB previously located
in `candidate_db.py`.  No behavioural changes – only the import path is
new (``from astrobatch.db import CandidateDB``).
"""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Dict, Any, Optional, Iterable

SCHEMA = """
CREATE TABLE IF NOT EXISTS candidates (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    task_id INTEGER NOT NULL,
    cutout TEXT,
    ra REAL,
    dec REAL,
    mag REAL,
    magerr REAL,
    flux REAL,
    fluxerr REAL,
    x REAL,
    y REAL,
    xerr REAL,
    yerr REAL,
    score REAL,
    flags TEXT,
    fwhm REAL,
    flux_radius REAL,
    theta REAL,
    bg REAL,
    number INTEGER,
    mag_auto REAL,
    isoarea_image REAL,
    bg_fluxerr REAL,
    bg_local REAL,
    mag_calib REAL,
    mag_calib_err REAL,
    mag_limit REAL,
    mag_filter_name TEXT,
    mag_color_name TEXT,
    mag_color_term TEXT,
    near_star INTEGER,
    near_galaxy INTEGER,
    a REAL,
    b REAL,
    is_hpm INTEGER,
    label TEXT DEFAULT 'unknown',
    inside_galaxy INTEGER,
    peak_snr REAL,
    flux_pos REAL,
    flux_neg REAL,
    asymmetry REAL,
    footprint_area REAL,
    UNIQUE(task_id, cutout)
);
"""


class CandidateDB:
    """Simple SQLite backend for STDWeb candidate records."""

    def __init__(self, db_path: Path | str):
        self.db_path = Path(db_path)
        self.conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
        self._ensure_schema()

    # ------------------------------------------------------------------
    # Schema handling
    # ------------------------------------------------------------------

    def _ensure_schema(self) -> None:
        """Create table & add missing columns when upgrading."""
        self.conn.execute(SCHEMA)
        self.conn.commit()

        required_cols: dict[str, str] = {
            "a": "REAL",
            "b": "REAL",
            "is_hpm": "INTEGER",
            "magerr": "REAL",
            "flux": "REAL",
            "fluxerr": "REAL",
            "x": "REAL",
            "y": "REAL",
            "xerr": "REAL",
            "yerr": "REAL",
            "theta": "REAL",
            "bg": "REAL",
            "number": "INTEGER",
            "mag_auto": "REAL",
            "isoarea_image": "REAL",
            "bg_fluxerr": "REAL",
            "bg_local": "REAL",
            "mag_calib": "REAL",
            "mag_calib_err": "REAL",
            "mag_limit": "REAL",
            "mag_filter_name": "TEXT",
            "mag_color_name": "TEXT",
            "mag_color_term": "TEXT",
            "near_star": "INTEGER",
            "near_galaxy": "INTEGER",
            "inside_galaxy": "INTEGER",
            "peak_snr": "REAL",
            "flux_pos": "REAL",
            "flux_neg": "REAL",
            "asymmetry": "REAL",
            "footprint_area": "REAL",
        }

        existing = {row[1] for row in self.conn.execute("PRAGMA table_info(candidates)").fetchall()}
        for col, typ in required_cols.items():
            if col not in existing:
                self.conn.execute(f"ALTER TABLE candidates ADD COLUMN {col} {typ}")
        self.conn.commit()

    # ------------------------------------------------------------------
    # CRUD helpers
    # ------------------------------------------------------------------

    def upsert(self, cand: Dict[str, Any]) -> None:
        """Insert new candidate or update fields if already present."""
        cols = [
            "task_id",
            "cutout",
            "ra",
            "dec",
            "mag",
            "magerr",
            "flux",
            "fluxerr",
            "x",
            "y",
            "xerr",
            "yerr",
            "score",
            "flags",
            "fwhm",
            "flux_radius",
            "theta",
            "bg",
            "number",
            "mag_auto",
            "isoarea_image",
            "bg_fluxerr",
            "bg_local",
            "mag_calib",
            "mag_calib_err",
            "mag_limit",
            "mag_filter_name",
            "mag_color_name",
            "mag_color_term",
            "near_star",
            "near_galaxy",
            "inside_galaxy",
            "a",
            "b",
            "is_hpm",
            "peak_snr",
            "flux_pos",
            "flux_neg",
            "asymmetry",
            "footprint_area",
        ]
        placeholders = ",".join(["?"] * len(cols))
        values = [cand.get(c) for c in cols]
        self.conn.execute(
            f"INSERT OR IGNORE INTO candidates ({','.join(cols)}) VALUES ({placeholders})",
            values,
        )

        # Overwrite mutable fields (score, photometry, etc.)
        update_sql = (
            "UPDATE candidates SET "
            "ra=?, dec=?, mag=?, magerr=?, flux=?, fluxerr=?, x=?, y=?, xerr=?, yerr=?, "
            "score=?, flags=?, fwhm=?, flux_radius=?, theta=?, bg=?, number=?, mag_auto=?, isoarea_image=?, bg_fluxerr=?, "
            "bg_local=?, mag_calib=?, mag_calib_err=?, mag_limit=?, mag_filter_name=?, mag_color_name=?, mag_color_term=?, "
            "near_star=?, near_galaxy=?, inside_galaxy=?, a=?, b=?, is_hpm=?, peak_snr=?, flux_pos=?, flux_neg=?, asymmetry=?, footprint_area=? "
            "WHERE task_id=? AND cutout=?"
        )
        self.conn.execute(
            update_sql,
            (
                cand.get("ra"),
                cand.get("dec"),
                cand.get("mag"),
                cand.get("magerr"),
                cand.get("flux"),
                cand.get("fluxerr"),
                cand.get("x"),
                cand.get("y"),
                cand.get("xerr"),
                cand.get("yerr"),
                cand.get("score"),
                cand.get("flags"),
                cand.get("fwhm"),
                cand.get("flux_radius"),
                cand.get("theta"),
                cand.get("bg"),
                cand.get("number"),
                cand.get("mag_auto"),
                cand.get("isoarea_image"),
                cand.get("bg_fluxerr"),
                cand.get("bg_local"),
                cand.get("mag_calib"),
                cand.get("mag_calib_err"),
                cand.get("mag_limit"),
                cand.get("mag_filter_name"),
                cand.get("mag_color_name"),
                cand.get("mag_color_term"),
                cand.get("near_star"),
                cand.get("near_galaxy"),
                cand.get("inside_galaxy"),
                cand.get("a"),
                cand.get("b"),
                cand.get("is_hpm"),
                cand.get("peak_snr"),
                cand.get("flux_pos"),
                cand.get("flux_neg"),
                cand.get("asymmetry"),
                cand.get("footprint_area"),
                cand.get("task_id"),
                cand.get("cutout"),
            ),
        )
        self.conn.commit()

    def set_label(self, task_id: int, cutout: str, label: str) -> None:
        self.conn.execute(
            "UPDATE candidates SET label=? WHERE task_id=? AND cutout=?",
            (label, task_id, cutout),
        )
        self.conn.commit()

    def fetch_all(self, where: Optional[str] = None) -> Iterable[sqlite3.Row]:
        sql = "SELECT * FROM candidates"
        if where:
            sql += " WHERE " + where
        return self.conn.execute(sql)

    def get_label(self, task_id: int, cutout: str) -> Optional[str]:
        row = self.conn.execute(
            "SELECT label FROM candidates WHERE task_id=? AND cutout=?",
            (task_id, cutout),
        ).fetchone()
        return row["label"] if row else None 