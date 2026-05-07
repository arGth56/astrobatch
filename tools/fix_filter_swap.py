#!/usr/bin/env python3
"""
fix_filter_swap.py — correct BP/RP filter mislabelling in raw FITS files.

Filter wheel was wired wrong:
  Position 1 labelled BP  ->  physically RP
  Position 2 labelled RP  ->  physically BP
  Position 3 G             ->  correct
  Position 4 Dark          ->  correct

Usage:
  python3 fix_filter_swap.py --dry-run          # print changes, touch nothing
  python3 fix_filter_swap.py --execute          # apply all changes
  python3 fix_filter_swap.py --execute --date 2026-03-31  # single date only
  python3 fix_filter_swap.py --execute --step raw         # one step only
"""

import argparse, sys, shutil, datetime
from pathlib import Path

NAS_INPUT  = Path("/mnt/nas/input/pyl/astro/input")
NAS_OUTPUT = Path("/mnt/nas/input/pyl/astro/output")
FLAT_DIR   = Path("/home/pyl/Documents/astrobatch/calib/flats")
DB_PATH    = Path("/home/pyl/Documents/astrobatch/nightmanager.db")
TMP        = "FSTMP"

def swapped(s):
    if "_BP_" in s: return s.replace("_BP_", "_RP_", 1)
    if "_RP_" in s: return s.replace("_RP_", "_BP_", 1)
    return s

def to_tmp(s):
    if "_BP_" in s: return s.replace("_BP_", f"_{TMP}_", 1)
    if "_RP_" in s: return s.replace("_RP_", f"_{TMP}_", 1)
    return s

def find_raw(date_filter):
    for p in sorted(NAS_INPUT.glob("*/SNAPSHOT/*.fits")):
        if date_filter and date_filter not in str(p): continue
        n = p.name
        if TMP in n: continue
        if "_BP_" in n or "_RP_" in n:
            yield p

def patch_header(path):
    from astropy.io import fits
    with fits.open(path, mode="update") as hdul:
        old = str(hdul[0].header.get("FILTER", "")).strip()
        if old == "BP":   hdul[0].header["FILTER"] = "RP"
        elif old == "RP": hdul[0].header["FILTER"] = "BP"
        hdul.flush()

def fix_raw(dry, date_filter):
    files = list(find_raw(date_filter))
    bp = [p for p in files if "_BP_" in p.name]
    rp = [p for p in files if "_RP_" in p.name]
    print(f"\n{'[DRY] ' if dry else ''}=== Step 1: Raw FITS  BP:{len(bp)}  RP:{len(rp)} ===")
    errors = 0

    def mv(src, dst, patch=False):
        nonlocal errors
        if dry:
            print(f"  RENAME  {src.name}  ->  {dst.name}")
            return
        try:
            src.rename(dst)
            if patch: patch_header(dst)
        except Exception as e:
            print(f"  ERROR {src}: {e}", file=sys.stderr); errors += 1

    print("  Pass A: BP -> TMP")
    for p in bp:  mv(p, p.parent / to_tmp(p.name))
    print("  Pass B: RP -> BP")
    for p in rp:  mv(p, p.parent / swapped(p.name), patch=True)
    print("  Pass C: TMP -> RP")
    for p in bp:  mv(p.parent / to_tmp(p.name), p.parent / swapped(p.name), patch=True)

    if dry: print(f"  Would rename {len(files)} files + patch FILTER headers.")
    else:   print(f"  Done. {len(files)} files renamed, {errors} errors.")

def fix_flats(dry):
    print(f"\n{'[DRY] ' if dry else ''}=== Step 2: Master flats ===")
    for ext in ("", ".bak"):
        gbp = FLAT_DIR / f"flat_Gbp.fit{ext}"
        grp = FLAT_DIR / f"flat_Grp.fit{ext}"
        tmp = FLAT_DIR / f"flat_Gbp_FSTMP.fit{ext}"
        if not gbp.exists() and not grp.exists():
            print(f"  Skip: {gbp.name} / {grp.name}")
            continue
        if dry:
            print(f"  RENAME  {gbp.name}  <->  {grp.name}")
            continue
        if gbp.exists(): gbp.rename(tmp)
        if grp.exists(): grp.rename(gbp)
        if tmp.exists(): tmp.rename(grp)
        print(f"  Swapped: {gbp.name} <-> {grp.name}")

def delete_stacks(dry, date_filter):
    dirs = sorted(NAS_OUTPUT.glob("*/*/BP")) + sorted(NAS_OUTPUT.glob("*/*/RP"))
    if date_filter: dirs = [d for d in dirs if date_filter in str(d)]
    print(f"\n{'[DRY] ' if dry else ''}=== Step 3: Delete stacked BP/RP outputs ({len(dirs)} dirs) ===")
    for d in dirs:
        if dry: print(f"  DELETE  {d}")
        else:
            shutil.rmtree(d, ignore_errors=True)
            print(f"  Deleted {d}")
    if not dry: print(f"  Removed {len(dirs)} directories.")

def fix_db(dry, date_filter):
    import sqlite3
    print(f"\n{'[DRY] ' if dry else ''}=== Step 4: Database ===")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    where = f"AND obs_date='{date_filter}'" if date_filter else ""
    c.execute(f"SELECT filter,count(*) FROM pipeline_results WHERE filter IN ('BP','RP') {where} GROUP BY filter")
    counts = dict(c.fetchall())
    print(f"  BP rows: {counts.get('BP',0)}  RP rows: {counts.get('RP',0)}")
    if dry:
        print("  Would swap BP<->RP and reset all to status='pending'.")
        conn.close(); return
    ts  = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = DB_PATH.parent / f"nightmanager_{ts}.db.bak"
    shutil.copy2(DB_PATH, bak)
    print(f"  DB backup -> {bak.name}")
    c.execute(f"UPDATE pipeline_results SET filter='BP_TMP' WHERE filter='BP' {where}")
    c.execute(f"UPDATE pipeline_results SET filter='BP'     WHERE filter='RP'  {where}")
    c.execute(f"UPDATE pipeline_results SET filter='RP'     WHERE filter='BP_TMP' {where}")
    c.execute(f"""UPDATE pipeline_results
                  SET status='pending', error=NULL,
                      stdweb_task_id=NULL, stdweb_url=NULL, stdweb_state=NULL,
                      mag_ap=NULL, magerr_ap=NULL,
                      mag_sub=NULL, magerr_sub=NULL, mag_sub_ul=NULL, mjd=NULL
                  WHERE filter IN ('BP','RP') {where}""")
    conn.commit(); conn.close()
    print("  Done.")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run",  action="store_true")
    p.add_argument("--execute",  action="store_true")
    p.add_argument("--date",     default=None)
    p.add_argument("--step",     default="all", choices=["all","raw","flats","stacks","db"])
    args = p.parse_args()
    if not args.dry_run and not args.execute:
        p.error("Specify --dry-run or --execute")
    dry = args.dry_run; date = args.date; step = args.step
    print(f"Mode={'DRY-RUN' if dry else 'EXECUTE'}"
          + (f"  date={date}" if date else "")
          + (f"  step={step}" if step != "all" else ""))
    if step in ("all","raw"):               fix_raw(dry, date)
    if step in ("all","flats") and not date: fix_flats(dry)
    if step in ("all","stacks"):            delete_stacks(dry, date)
    if step in ("all","db"):                fix_db(dry, date)
    print("\n" + ("Dry-run complete -- nothing changed." if dry else "All done."))

if __name__ == "__main__":
    main()
