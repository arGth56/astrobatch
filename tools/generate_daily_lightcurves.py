#!/usr/bin/env python3
"""Generate light curves for targets processed in the last 24 hours.

Intended to run daily at noon via systemd timer (astrobatch-lightcurves.timer).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def _load_notion_env() -> None:
    """Load config/notion.env if present (gitignored)."""
    import os

    env_path = PROJECT_ROOT / "config" / "notion.env"
    if not env_path.is_file():
        return
    for line in env_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        os.environ.setdefault(key.strip(), val.strip())

from astrobatch.lightcurves import (  # noqa: E402
    DEFAULT_DB,
    DEFAULT_OUTPUT,
    generate_for_target,
    recent_targets,
    run_daily,
    target_display_name,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


def main() -> int:
    _load_notion_env()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hours", type=float, default=24, help="Look-back window (default: 24)")
    p.add_argument("--db", type=Path, default=DEFAULT_DB, help="nightmanager.db path")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output root directory")
    p.add_argument(
        "--target",
        action="append",
        dest="targets",
        metavar="SLUG",
        help="Only these target(s), e.g. sn_2026fvx (skips recent-target filter)",
    )
    p.add_argument(
        "--list-recent",
        action="store_true",
        help="Print targets updated in the last --hours and exit",
    )
    p.add_argument(
        "--no-notion",
        action="store_true",
        help="Skip Notion upload even if NOTION_TOKEN is set",
    )
    p.add_argument(
        "--notion-only",
        action="store_true",
        help="Upload existing index.json for today (no regeneration)",
    )
    args = p.parse_args()

    if args.notion_only:
        import json
        from datetime import datetime

        from astrobatch.notion_upload import upload_daily_summary

        day_dir = args.output / datetime.now().strftime("%Y-%m-%d")
        index = day_dir / "index.json"
        if not index.is_file():
            log.error("No index at %s", index)
            return 1
        summary = json.loads(index.read_text(encoding="utf-8"))
        if args.no_notion:
            return 0
        upload_daily_summary(summary)
        return 0

    if args.list_recent:
        for t in recent_targets(args.db, hours=args.hours):
            print(t)
        return 0

    if args.no_notion:
        import os

        os.environ.pop("NOTION_TOKEN", None)
        os.environ.pop("NOTION_API_KEY", None)

    if args.targets:
        from datetime import datetime

        day_dir = args.output / datetime.now().strftime("%Y-%m-%d")
        day_dir.mkdir(parents=True, exist_ok=True)
        summary = {
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "output_dir": str(day_dir),
            "targets": [],
        }
        for target in args.targets:
            paths = generate_for_target(target, day_dir, args.db)
            entry = {"target": target, "display": target_display_name(target), **paths}
            summary["targets"].append(entry)
            log.info("%s → %s", target, paths)
        if not args.no_notion:
            import os

            if os.environ.get("NOTION_TOKEN") or os.environ.get("NOTION_API_KEY"):
                from astrobatch.notion_upload import upload_daily_summary

                upload_daily_summary(summary)
        return 0

    summary = run_daily(args.db, args.output, hours=args.hours)
    n = len(summary.get("targets", []))
    log.info("Done: %d target(s) → %s", n, summary.get("output_dir"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
