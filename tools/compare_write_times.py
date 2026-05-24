#!/usr/bin/env python3
"""Compare FITS write timing: NAS vs astropc (D:) using mtime and filename timestamps."""
import re
import statistics
import sys
from datetime import datetime
from pathlib import Path

FNAME_RE = re.compile(
    r"^(?P<date>\d{4}-\d{2}-\d{2})_(?P<time>\d{2}-\d{2}-\d{2})_(?P<filt>[^_]+)_.*?(?P<exp>\d+(?:\.\d+)?)s_\d+\.fits$"
)


def parse_name(name: str):
    m = FNAME_RE.match(name)
    if not m:
        return None
    d, t, filt, exp = m.group("date", "time", "filt", "exp")
    ts = datetime.strptime(f"{d} {t.replace('-', ':')}", "%Y-%m-%d %H:%M:%S").timestamp()
    return ts, filt, float(exp)


def load_snap(root: str):
    snap = Path(root)
    if not snap.is_dir():
        return []
    rows = []
    for f in snap.glob("*.fits"):
        p = parse_name(f.name)
        if not p:
            continue
        fts, filt, exp = p
        st = f.stat()
        rows.append(
            {
                "name": f.name,
                "fts": fts,
                "filt": filt,
                "exp": exp,
                "mtime": st.st_mtime,
                "size": st.st_size,
            }
        )
    rows.sort(key=lambda r: r["fts"])
    return rows


def find_blocks(rows, min_n=10, max_gap=90):
    blocks = []
    cur = []
    for r in rows:
        if not cur:
            cur = [r]
            continue
        gap_name = r["fts"] - cur[-1]["fts"]
        same_filt = r["filt"] == cur[-1]["filt"]
        if same_filt and 25 <= gap_name <= max_gap:
            cur.append(r)
        else:
            if len(cur) >= min_n:
                blocks.append(cur)
            cur = [r]
    if len(cur) >= min_n:
        blocks.append(cur)
    return blocks


def block_metrics(block):
    name_gaps = [block[i]["fts"] - block[i - 1]["fts"] for i in range(1, len(block))]
    mtime_gaps = [block[i]["mtime"] - block[i - 1]["mtime"] for i in range(1, len(block))]
    tight = [(ng, mg) for ng, mg in zip(name_gaps, mtime_gaps) if 25 <= ng <= 90]
    if len(tight) < 8:
        return None
    ng = [x[0] for x in tight]
    mg = [x[1] for x in tight]
    return {
        "n": len(block),
        "filt": block[0]["filt"],
        "exp": block[0]["exp"],
        "start": datetime.fromtimestamp(block[0]["fts"]),
        "end": datetime.fromtimestamp(block[-1]["fts"]),
        "name_gap_med": statistics.median(ng),
        "mtime_gap_med": statistics.median(mg),
        "mtime_gap_mean": statistics.mean(mg),
        "mtime_gap_p90": sorted(mg)[max(0, int(0.9 * len(mg)) - 1)],
    }


def compare_mounts(date):
    nas = load_snap(f"/mnt/nas/input/pyl/astro/input/{date}/SNAPSHOT")
    pc = load_snap(f"/mnt/astropc/{date}/SNAPSHOT")
    nas_by = {r["name"]: r for r in nas}
    pc_by = {r["name"]: r for r in pc}
    common = sorted(set(nas_by) & set(pc_by))
    lags = [nas_by[n]["mtime"] - pc_by[n]["mtime"] for n in common]
    return nas, pc, lags


def main():
    dates = sys.argv[1:] or ["2026-05-21", "2026-05-23"]
    for date in dates:
        print("=" * 72)
        print(date)
        nas, pc, lags = compare_mounts(date)
        print(f"  NAS files: {len(nas)}   astropc: {len(pc)}   matched: {len(lags)}")
        if lags:
            print(
                f"  NAS-astropc mtime lag: med={statistics.median(lags):.1f}s "
                f"mean={statistics.mean(lags):.1f}s p90={sorted(lags)[int(0.9*len(lags))-1]:.1f}s "
                f"min={min(lags):.1f}s max={max(lags):.1f}s"
            )
        for label, rows in [("NAS", nas), ("astropc (D:)", pc)]:
            metrics = [m for b in find_blocks(rows) if (m := block_metrics(b))]
            print(f"  {label}: {len(metrics)} consecutive blocks")
            for m in sorted(metrics, key=lambda x: -x["n"])[:4]:
                print(
                    f"    {m['filt']:4s} n={m['n']:3d} exp={m['exp']:.0f}s  "
                    f"fname_delta={m['name_gap_med']:.1f}s  mtime_delta med={m['mtime_gap_med']:.2f}s "
                    f"mean={m['mtime_gap_mean']:.2f}s p90={m['mtime_gap_p90']:.2f}s  "
                    f"{m['start'].strftime('%H:%M:%S')}-{m['end'].strftime('%H:%M:%S')}"
                )


if __name__ == "__main__":
    main()
