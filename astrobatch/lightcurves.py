"""Generate aperture and subtraction light curves from nightmanager.db."""

from __future__ import annotations

import json
import logging
import os
import sqlite3
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Literal

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = Path(os.environ.get("NIGHTMANAGER_DB", PROJECT_ROOT / "nightmanager.db"))
DEFAULT_OUTPUT = Path(
    os.environ.get("LIGHTCURVE_OUTPUT", PROJECT_ROOT / "data" / "lightcurves")
)

FILT_COLORS = {"G": "#2ecc71", "BP": "#3498db", "RP": "#e74c3c"}
FILT_MARKERS = {"G": "o", "BP": "s", "RP": "^"}
FILTER_ORDER = ["G", "BP", "RP"]

PhotMode = Literal["aperture", "subtraction"]


@dataclass(frozen=True)
class PhotPoint:
    mjd: float
    mag: float
    magerr: float
    obs_date: str
    is_upper_limit: bool = False


def target_display_name(slug: str) -> str:
    """sn_2026fvx → SN 2026fvx."""
    s = slug.replace("_", " ")
    if s.lower().startswith("sn "):
        return "SN" + s[2:]
    if s.lower().startswith("at "):
        return "AT" + s[2:]
    return s.upper() if s.islower() else s


def _db_path(db_path: Path | None) -> Path:
    path = db_path or DEFAULT_DB
    if not path.exists():
        raise FileNotFoundError(f"Database not found: {path}")
    return path


def recent_targets(
    db_path: Path | None = None,
    *,
    hours: float = 24,
) -> list[str]:
    """Targets with pipeline_results updated in the last *hours* (done, with photometry)."""
    path = _db_path(db_path)
    cutoff = (datetime.now() - timedelta(hours=hours)).strftime("%Y-%m-%d %H:%M:%S")
    conn = sqlite3.connect(path)
    try:
        rows = conn.execute(
            """
            SELECT DISTINCT target
            FROM pipeline_results
            WHERE status = 'done'
              AND updated_at >= ?
              AND (
                    mag_ap IS NOT NULL
                    OR mag_sub IS NOT NULL
                    OR mag_sub_ul IS NOT NULL
                  )
            ORDER BY target
            """,
            (cutoff,),
        ).fetchall()
    finally:
        conn.close()
    return [r[0] for r in rows if r[0]]


def load_points(
    target: str,
    mode: PhotMode,
    db_path: Path | None = None,
) -> dict[str, list[PhotPoint]]:
    """Load deduplicated photometry per filter for one target."""
    path = _db_path(db_path)
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    try:
        if mode == "aperture":
            rows = conn.execute(
                """
                SELECT filter, mjd, mag_ap AS mag, magerr_ap AS magerr,
                       obs_date, NULL AS mag_ul
                FROM pipeline_results
                WHERE target = ?
                  AND status = 'done'
                  AND mag_ap IS NOT NULL
                  AND mjd IS NOT NULL
                ORDER BY mjd, filter
                """,
                (target,),
            ).fetchall()
        else:
            rows = conn.execute(
                """
                SELECT filter, mjd, mag_sub AS mag, magerr_sub AS magerr,
                       obs_date, mag_sub_ul AS mag_ul
                FROM pipeline_results
                WHERE target = ?
                  AND status = 'done'
                  AND mjd IS NOT NULL
                  AND (mag_sub IS NOT NULL OR mag_sub_ul IS NOT NULL)
                ORDER BY mjd, filter
                """,
                (target,),
            ).fetchall()
    finally:
        conn.close()

    groups: dict[tuple[str, float], list] = defaultdict(list)
    for r in rows:
        groups[(r["filter"], round(r["mjd"], 4))].append(r)

    by_filter: dict[str, list[PhotPoint]] = {f: [] for f in FILTER_ORDER}
    for (filt, mjd_key), grp in sorted(groups.items()):
        if filt not in by_filter:
            continue
        ul_rows = [x for x in grp if x["mag"] is None and x["mag_ul"] is not None]
        det_rows = [x for x in grp if x["mag"] is not None]

        if det_rows:
            errs = [max(x["magerr"] or 0.01, 1e-3) for x in det_rows]
            w = 1.0 / np.array(errs) ** 2
            mag = float(np.average([x["mag"] for x in det_rows], weights=w))
            err = float(np.sqrt(1.0 / w.sum()))
            by_filter[filt].append(
                PhotPoint(
                    mjd=mjd_key,
                    mag=mag,
                    magerr=err,
                    obs_date=det_rows[0]["obs_date"] or "",
                    is_upper_limit=False,
                )
            )
        elif ul_rows and mode == "subtraction":
            ul = ul_rows[0]
            by_filter[filt].append(
                PhotPoint(
                    mjd=mjd_key,
                    mag=float(ul["mag_ul"]),
                    magerr=0.0,
                    obs_date=ul["obs_date"] or "",
                    is_upper_limit=True,
                )
            )
    return by_filter


def _mjd_ref(by_filter: dict[str, list[PhotPoint]]) -> float:
    all_mjd = [p.mjd for pts in by_filter.values() for p in pts]
    if not all_mjd:
        return 61100.0
    return float(int(min(all_mjd) / 50) * 50)


def plot_lightcurve(
    by_filter: dict[str, list[PhotPoint]],
    *,
    title: str,
    ylabel: str,
    output_path: Path,
    mode: PhotMode,
) -> bool:
    """Save a light-curve PNG. Returns False if there is no data to plot."""
    has_data = any(by_filter.get(f) for f in FILTER_ORDER)
    if not has_data:
        return False

    mjd_ref = _mjd_ref(by_filter)
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=120)
    fig.patch.set_facecolor("#0d1117")
    ax.set_facecolor("#161b22")

    for filt in FILTER_ORDER:
        pts = by_filter.get(filt, [])
        if not pts:
            continue
        det = [p for p in pts if not p.is_upper_limit]
        ul = [p for p in pts if p.is_upper_limit]

        if det:
            mjd = np.array([p.mjd for p in det]) - mjd_ref
            mag = np.array([p.mag for p in det])
            err = np.array([p.magerr for p in det])
            ax.errorbar(
                mjd,
                mag,
                yerr=err,
                fmt=FILT_MARKERS[filt],
                color=FILT_COLORS[filt],
                markeredgecolor="white",
                markeredgewidth=0.4,
                markersize=6,
                capsize=2.5,
                elinewidth=1,
                label=filt,
                zorder=5,
            )
            ax.plot(
                mjd,
                mag,
                "-",
                color=FILT_COLORS[filt],
                alpha=0.25,
                linewidth=0.8,
                zorder=4,
            )

        if ul and mode == "subtraction":
            mjd_ul = np.array([p.mjd for p in ul]) - mjd_ref
            mag_ul = np.array([p.mag for p in ul])
            ax.plot(
                mjd_ul,
                mag_ul,
                "v",
                color=FILT_COLORS[filt],
                markersize=7,
                markeredgecolor="white",
                markeredgewidth=0.4,
                label=f"{filt} (UL)" if not det else None,
                zorder=5,
            )

    ax.invert_yaxis()
    ax.set_xlabel(f"MJD − {mjd_ref:.0f}", color="#c9d1d9", fontsize=11)
    ax.set_ylabel(ylabel, color="#c9d1d9", fontsize=11)
    ax.set_title(title, color="#f0f6fc", fontsize=13, fontweight="bold")
    ax.tick_params(colors="#8b949e")
    for spine in ax.spines.values():
        spine.set_color("#30363d")
    ax.grid(True, alpha=0.15, color="#8b949e")
    ax.legend(title="Filter", framealpha=0.9, loc="lower left")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, facecolor=fig.get_facecolor())
    plt.close(fig)
    return True


def generate_for_target(
    target: str,
    output_dir: Path,
    db_path: Path | None = None,
) -> dict[str, str | None]:
    """Write aperture and subtraction PNGs for *target*. Returns output paths."""
    name = target_display_name(target)
    out: dict[str, str | None] = {"aperture": None, "subtraction": None}

    ap_data = load_points(target, "aperture", db_path)
    ap_path = output_dir / f"{target}_aperture.png"
    if plot_lightcurve(
        ap_data,
        title=f"{name} — aperture photometry",
        ylabel="Gaia magnitude (aperture)",
        output_path=ap_path,
        mode="aperture",
    ):
        out["aperture"] = str(ap_path)

    sub_data = load_points(target, "subtraction", db_path)
    sub_path = output_dir / f"{target}_subtraction.png"
    if plot_lightcurve(
        sub_data,
        title=f"{name} — template subtraction",
        ylabel="Gaia magnitude (subtraction)",
        output_path=sub_path,
        mode="subtraction",
    ):
        out["subtraction"] = str(sub_path)

    return out


def run_daily(
    db_path: Path | None = None,
    output_root: Path | None = None,
    *,
    hours: float = 24,
    targets: list[str] | None = None,
) -> dict:
    """Generate light curves for recently processed targets."""
    root = output_root or DEFAULT_OUTPUT
    day_dir = root / datetime.now().strftime("%Y-%m-%d")
    day_dir.mkdir(parents=True, exist_ok=True)

    if targets is None:
        targets = recent_targets(db_path, hours=hours)

    summary: dict = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "hours": hours,
        "output_dir": str(day_dir),
        "targets": [],
    }

    if not targets:
        log.info("No targets processed in the last %.0f h", hours)
        summary["message"] = "no targets"
        _write_index(day_dir, summary)
        return summary

    for target in targets:
        paths = generate_for_target(target, day_dir, db_path)
        entry = {"target": target, "display": target_display_name(target), **paths}
        summary["targets"].append(entry)
        log.info(
            "%s: aperture=%s subtraction=%s",
            target,
            "yes" if paths["aperture"] else "skip",
            "yes" if paths["subtraction"] else "skip",
        )

    _write_index(day_dir, summary)

    if os.environ.get("NOTION_TOKEN") or os.environ.get("NOTION_API_KEY"):
        try:
            from astrobatch.notion_upload import upload_daily_summary

            notion = upload_daily_summary(summary)
            if notion:
                summary["notion"] = notion
                _write_index(day_dir, summary)
        except Exception as exc:
            log.error("Notion upload failed (light curves saved locally): %s", exc)
            summary["notion_error"] = str(exc)
            _write_index(day_dir, summary)

    return summary


def _write_index(day_dir: Path, summary: dict) -> None:
    index_path = day_dir / "index.json"
    index_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
