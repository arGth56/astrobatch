import types
import importlib
from pathlib import Path

import astrobatch.cli as cli


def test_cli_upload_triggers_upload(monkeypatch, tmp_path):
    """Verify that running the CLI with --upload calls spliter.upload_processed_images
    without importing the heavy astropy stack or triggering network requests.
    """
    # ------------------------------------------------------------------
    # 1. Prepare a fake NIGHT/LIGHT directory (structure is irrelevant for the
    #    dry-run, but create it for completeness).
    # ------------------------------------------------------------------
    night_root: Path = tmp_path / "2025-07-21" / "LIGHT"
    night_root.mkdir(parents=True)

    # ------------------------------------------------------------------
    # 2. Build a minimal fake replacement for the `astrobatch.spliter` module.
    #    We only provide the attributes that the CLI touches.
    # ------------------------------------------------------------------
    fake_spliter = types.ModuleType("astrobatch.spliter")
    called_flag = {"upload": False}

    def dummy_upload_processed_images():
        called_flag["upload"] = True

    # Stubbed no-op functions/attributes
    fake_spliter.upload_processed_images = dummy_upload_processed_images
    fake_spliter.main = lambda: None
    fake_spliter.calibrate_folders_pysiril = lambda _: None
    fake_spliter.plate_solve_all_res_files = lambda: None
    fake_spliter.DATA_ROOT = night_root
    fake_spliter.init_path = str(night_root)
    fake_spliter.directory = night_root

    # ------------------------------------------------------------------
    # 3. Monkey-patch importlib so that the CLI receives our fake module when
    #    it lazily imports "astrobatch.spliter".
    # ------------------------------------------------------------------
    original_import_module = importlib.import_module

    def fake_import(name, *args, **kwargs):
        if name == "astrobatch.spliter":
            return fake_spliter
        return original_import_module(name, *args, **kwargs)

    monkeypatch.setattr(importlib, "import_module", fake_import)

    # ------------------------------------------------------------------
    # 4. Execute the CLI. We omit --dry-run because we *expect* the real call
    #    to happen (but against the stub, so it's safe).
    # ------------------------------------------------------------------
    cli._main(["--upload", "--root", str(night_root)])

    # ------------------------------------------------------------------
    # 5. Assert that our stub was invoked.
    # ------------------------------------------------------------------
    assert called_flag["upload"], "upload_processed_images() was not triggered" 