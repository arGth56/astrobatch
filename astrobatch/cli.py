"""Command-line interface for the *astrobatch* pipeline.

This small wrapper preserves the original functions that still live in
`spliter.py`, providing a clean entry-point:

    python -m astrobatch.cli --split --calibrate --upload --analyse \
        --root /data/2025-07-16/LIGHT

Once the old monolithic script is refactored into package modules these
imports will be updated accordingly.
"""
from __future__ import annotations

import argparse
import importlib
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------
# Lazy import helper – avoids importing heavy dependencies if the user
# only wants a subset of the pipeline.
# ---------------------------------------------------------------------

def _lazy_import(module_name: str):
    return importlib.import_module(module_name)


def _main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="astrobatch",
        description="Professional nightly pipeline – split, calibrate, upload, analyse",
    )
    parser.add_argument("--root", type=Path, required=True, help="Night LIGHT directory root")
    parser.add_argument("--split", action="store_true", help="Split raw FITS into folder tree")
    parser.add_argument("--calibrate", action="store_true", help="Calibrate/stack via Siril (pySiril)")
    parser.add_argument("--process", action="store_true", help="Shortcut for --split --calibrate")
    parser.add_argument("--upload", action="store_true", help="Upload stacked images to STDWeb")
    parser.add_argument("--mosaic", action="store_true", help="Create sky mosaic from res.fit files in night root")
    parser.add_argument("--analyse", action="store_true", help="Analyse STDWeb results & transients")
    parser.add_argument("--plate-solve", action="store_true", help="Run stand-alone plate solving on res.fit")
    parser.add_argument("--dry-run", action="store_true", help="Show steps without executing")
    parser.add_argument("--start-server", action="store_true", help="Launch Flask candidate review server")
    parser.add_argument("--port", type=int, default=None, help="Port for --start-server (default 5100 or $PORT)")
    args, extra = parser.parse_known_args(argv)

    # --process implies both split and calibrate
    if args.process:
        args.split = True
        args.calibrate = True

    # ------------------------------------------------------------------
    # Always use CLI-only calibration (pySiril disabled).  The environment
    # flag is picked up inside spliter.py.
    # ------------------------------------------------------------------
    os.environ["ASTROBATCH_CLI_ONLY"] = "1"

    # Ensure at least one action selected
    if not any((args.split, args.calibrate, args.upload, args.analyse, args.plate_solve, args.start_server, args.mosaic)):
        parser.error("Select at least one action: --split / --calibrate / --upload / --analyse / --plate-solve / --mosaic / --start-server")

    # Import processing module (now inside package)
    spliter = _lazy_import("astrobatch.spliter")

    # Override DATA_ROOT dynamically
    if hasattr(spliter, "DATA_ROOT"):
        spliter.DATA_ROOT = Path(args.root)
        spliter.init_path = str(args.root)
        spliter.directory = Path(args.root)

    if args.dry_run:
        print("ℹ️  DRY-RUN – following actions would be executed:")

    # Start server early if requested (it blocks)
    if args.start_server:
        if args.dry_run:
            print("• would start candidate review server")
        else:
            from importlib import import_module
            srv = import_module("astrobatch.server")
            port = args.port or srv.DEFAULT_PORT
            print(f"🚀 Starting candidate review server on port {port} …")
            srv.start_server(port=port)
            return  # server is blocking, nothing else executes

    if args.split:
        if args.dry_run:
            print("• split raw FITS")
        else:
            import sys as _sys
            _old = _sys.argv[:]
            _sys.argv = ["spliter"]
            try:
                import os as _os
                if not args.calibrate:
                    _os.environ["ASTROBATCH_SKIP_CALIBRATE"] = "1"
                spliter.main()
                _os.environ.pop("ASTROBATCH_SKIP_CALIBRATE", None)
            finally:
                _sys.argv = _old

    if args.calibrate:
        # ------------------------------------------------------------------
        # Early check: ensure either pySiril is installed or a usable siril-cli
        # binary is available.  This prevents cryptic errors on systems that
        # ship an old Siril version.
        # ------------------------------------------------------------------
        import shutil, importlib.util

        has_pysiril = importlib.util.find_spec("pysiril") is not None
        has_siril_cli = shutil.which("siril-cli") is not None or shutil.which("siril") is not None

        if not has_pysiril and not has_siril_cli:
            print("⚠️  Neither pySiril nor siril-cli found – skipping calibration.\n"
                  "    Install Siril ≥1.2 and pySiril, or run without --calibrate.")
        else:
            if args.dry_run:
                print("• calibrate via Siril")
            else:
                spliter.calibrate_folders_pysiril({})  # fallback to internal logic

    if args.plate_solve:
        if args.dry_run:
            print("• plate solve all res.fit")
        else:
            spliter.plate_solve_all_res_files()

    if args.upload:
        if args.dry_run:
            print("• upload to STDWeb")
        else:
            spliter.upload_processed_images()

    if args.mosaic:
        if args.dry_run:
            print("• create sky mosaic")
        else:
            # Build a minimal arg namespace with defaults identical to those in spliter's parser
            from types import SimpleNamespace
            mosaic_args = SimpleNamespace(
                pixel_scale=1.0,
                projection="TAN",
                blend_method="feather",
                output_dir=str(Path(args.root) / "mosaics"),
                quality_threshold=0.7,
                limit=None,
                margin=2.0,
                no_size_limit=False,
                outlier_threshold=3.0,
                feather_sigma=2.0,
                fast_mode=False,
                incremental=False,
                edge_samples=8,
                no_border_csv=False,
            )
            spliter.create_mosaic(mosaic_args)

    if args.analyse:
        if args.dry_run:
            print("• analyse STDWeb results")
        else:
            try:
                analyse = _lazy_import("astrobatch.analyse")
                analyse.main(extra)
            except ModuleNotFoundError:
                print("⚠️  analyse.py not found – skipping analysis step")


if __name__ == "__main__":
    _main() 