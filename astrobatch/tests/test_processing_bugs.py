"""
test_processing_bugs.py
=======================

Tests for three confirmed bugs in processing_service.py.

Each test class:
  - States the bug clearly.
  - Has tests that FAIL before the fix and PASS after.
  - Includes a regression test to ensure correct paths still work.

────────────────────────────────────────────────────────────────────────────────
B1  wait_state return value ignored in photometry / subtraction steps
    → STDWeb job failures are silently recorded as "done"
    File: _resume_pipeline(), lines ~1151 and ~1160
    Fix:  capture return value; if "failed"/"error"/"timeout" → mark error, return

B2  confirm_selection has no job-status guard
    → any job's selection_result can be overwritten regardless of its state
    File: POST /jobs/<id>/selection, line ~1743
    Fix:  check status == "selection_pending"; return 409 otherwise

B3  /resume bypasses _active_jobs deduplication
    → calling /resume twice for the same result spawns duplicate workers
    File: POST /resume, line ~1693
    Fix:  maintain a _active_resumes set and return 409 on duplicate
────────────────────────────────────────────────────────────────────────────────
"""
import json
import sqlite3
import time
from unittest.mock import MagicMock, patch

import pytest

# Helpers are registered on pytest.helpers in conftest.py
from astrobatch.tests.conftest import (
    _seed_job,
    _seed_result,
    _get_job_status,
    _get_result_status,
    _get_selection_result,
    _TEST_DB,
)


# ─────────────────────────────────────────────────────────────────────────────
# B1 – wait_state return value is ignored for photometry / subtraction
# ─────────────────────────────────────────────────────────────────────────────

class TestB1WaitStateIgnored:
    """
    BUG
    ---
    _resume_pipeline() calls wait_state() for both the photometry and
    subtraction STDWeb steps but discards the return value. The pipeline then
    unconditionally calls:
        update_result(result_id, "done")
        update_job(job_id, "done")
    This means a failed/timed-out STDWeb step is recorded as a success.

    EXPECTED BEHAVIOUR (after fix)
    --------------------------------
    When wait_state() returns "failed", "error", or "timeout":
        • result.status must be "error"  (NOT "done")
        • job.status   must be "error"   (NOT "done")
        • the pipeline function must return early (not continue to next step)

    TEST SETUP NOTE
    ---------------
    _resume_pipeline() checks for res.fit (the stacked output) before reaching
    the STDWeb steps. We create a minimal fake res.fit in tmp_path so that
    guard passes, and then control the STDWeb behaviour via wait_state mocks.
    """

    # ── helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    def _make_wait_state_returning(state: str):
        """Drop-in replacement for wait_state() that immediately returns `state`."""
        def _mock(task_id, target_states, jlog, timeout=300):
            return state
        return _mock

    @staticmethod
    def _make_requests_mock(task_state: str = "failed"):
        """requests.get mock that reports `task_state` for /api/tasks/{id}/."""
        def _get(url, *args, **kwargs):
            resp = MagicMock()
            resp.json.return_value = {
                "state": task_state,
                "config": {"filter": "G", "target": "test_target"},
                "result": {},
            }
            resp.text = ""
            return resp
        return _get

    @staticmethod
    def _setup_fake_res_fit(tmp_path, ps):
        """
        Create the directory structure that _resume_pipeline expects to exist
        when resuming from 'photometry':

            DATA_DIR / obs_date / target / filter / exposures / res.fit

        Returns (data_dir, obs_date, target, filt, exposure) for DB seeding.
        """
        obs_date = "2026-01-01"
        target   = "test_target"
        filt     = "G"
        exposure = "60"

        work_dir = tmp_path / "data" / obs_date / target / filt / f"{exposure}s"
        work_dir.mkdir(parents=True)
        (work_dir / "res.fit").write_bytes(b"FAKE")   # minimal placeholder

        # Redirect the module's DATA_DIR to our temp directory
        ps.DATA_DIR = tmp_path / "data"

        return obs_date, target, filt, exposure

    # ── tests ─────────────────────────────────────────────────────────────────

    def test_photometry_failure_marks_job_error(self, ps, clean_db, tmp_path):
        """
        STDWeb photometry returns 'failed'.
        Expected: result='error', job='error'.
        Current (bug): result='done', job='done'.
        """
        obs_date, target, filt, exp = self._setup_fake_res_fit(tmp_path, ps)

        _seed_job(job_id=101, status="running")
        _seed_result(result_id=1001, job_id=101, status="photometry",
                     filter=filt, exposure=exp, target=target)

        with patch.object(ps, "wait_state", self._make_wait_state_returning("failed")), \
             patch("requests.get",  self._make_requests_mock("failed")), \
             patch("requests.post", MagicMock(
                 return_value=MagicMock(json=MagicMock(return_value={}))
             )), \
             patch("time.sleep", lambda _: None):
            ps.resume_pipeline(result_id=1001, from_step="photometry")

        result_status = _get_result_status(1001)
        job_status    = _get_job_status(101)

        assert result_status == "error", (
            f"[B1] Photometry failure: result.status should be 'error', got '{result_status}'. "
            "Fix: capture wait_state return value and call update_result('error') when failed."
        )
        assert job_status == "error", (
            f"[B1] Photometry failure: job.status should be 'error', got '{job_status}'."
        )

    def test_subtraction_failure_marks_job_error(self, ps, clean_db, tmp_path):
        """
        STDWeb subtraction returns 'failed'.
        Expected: result='error', job='error'.
        Current (bug): result='done', job='done'.
        """
        obs_date, target, filt, exp = self._setup_fake_res_fit(tmp_path, ps)

        _seed_job(job_id=102, status="running")
        _seed_result(result_id=1002, job_id=102, status="subtraction",
                     filter=filt, exposure=exp, target=target)

        with patch.object(ps, "wait_state", self._make_wait_state_returning("failed")), \
             patch("requests.get",  self._make_requests_mock("failed")), \
             patch("requests.post", MagicMock(
                 return_value=MagicMock(json=MagicMock(return_value={}))
             )), \
             patch("time.sleep", lambda _: None):
            ps.resume_pipeline(result_id=1002, from_step="subtraction")

        assert _get_result_status(1002) == "error", (
            "[B1] Subtraction failure: result.status should be 'error'. "
            "Fix: check wait_state return value for subtraction step too."
        )
        assert _get_job_status(102) == "error", (
            "[B1] Subtraction failure: job.status should be 'error'."
        )

    def test_photometry_timeout_marks_job_error(self, ps, clean_db, tmp_path):
        """
        STDWeb photometry times out (wait_state returns 'timeout').
        Expected: result='error', job='error'.
        """
        obs_date, target, filt, exp = self._setup_fake_res_fit(tmp_path, ps)

        _seed_job(job_id=103, status="running")
        _seed_result(result_id=1003, job_id=103, status="photometry",
                     filter=filt, exposure=exp, target=target)

        with patch.object(ps, "wait_state", self._make_wait_state_returning("timeout")), \
             patch("requests.get",  self._make_requests_mock("unknown")), \
             patch("requests.post", MagicMock(
                 return_value=MagicMock(json=MagicMock(return_value={}))
             )), \
             patch("time.sleep", lambda _: None):
            ps.resume_pipeline(result_id=1003, from_step="photometry")

        assert _get_result_status(1003) == "error", "[B1] Timeout should mark result as error."
        assert _get_job_status(103)    == "error", "[B1] Timeout should mark job as error."

    def test_success_still_marks_done(self, ps, clean_db, tmp_path):
        """
        Regression: when all STDWeb steps succeed, job must still reach 'done'.
        """
        obs_date, target, filt, exp = self._setup_fake_res_fit(tmp_path, ps)

        _seed_job(job_id=104, status="running")
        _seed_result(result_id=1004, job_id=104, status="photometry",
                     filter=filt, exposure=exp, target=target)

        call_count = [0]
        def _success_wait(task_id, target_states, jlog, timeout=300):
            call_count[0] += 1
            return "photometry_done" if call_count[0] == 1 else "subtraction_done"

        with patch.object(ps, "wait_state", _success_wait), \
             patch("requests.get",  self._make_requests_mock("subtraction_done")), \
             patch("requests.post", MagicMock(
                 return_value=MagicMock(json=MagicMock(return_value={}))
             )), \
             patch("time.sleep", lambda _: None):
            ps.resume_pipeline(result_id=1004, from_step="photometry")

        assert _get_result_status(1004) == "done", "[B1 regression] Success path: result must be 'done'."
        assert _get_job_status(104)     == "done", "[B1 regression] Success path: job must be 'done'."


# ─────────────────────────────────────────────────────────────────────────────
# B2 – confirm_selection has no job-status guard
# ─────────────────────────────────────────────────────────────────────────────

class TestB2ConfirmSelectionGuard:
    """
    BUG
    ---
    POST /jobs/<id>/selection writes selection_result for ANY job, regardless
    of its current status. A job in 'running', 'calibrating', 'error', etc.
    will have its selection_result overwritten, potentially unblocking a stale
    poll loop if the service is restarted at just the wrong moment.

    EXPECTED BEHAVIOUR (after fix)
    --------------------------------
    • Status 200 only when job.status == 'selection_pending'
    • Status 409 (Conflict) with success=False when status is anything else
    • Status 404 when the job does not exist
    """

    def test_rejects_running_job(self, ps, clean_db):
        """
        A 'running' job must not accept a selection confirmation.
        Expected: HTTP 409, success=False.
        Current (bug): HTTP 200, success=True.
        """
        _seed_job(job_id=201, status="running")

        client = ps.app.test_client()
        resp = client.post(
            "/jobs/201/selection",
            json={"files": ["/tmp/frame1.fits"]},
        )

        assert resp.status_code == 409, (
            f"[B2] Running job: expected HTTP 409, got {resp.status_code}. "
            "Fix: check status == 'selection_pending' before writing selection_result."
        )
        data = resp.get_json()
        assert data.get("success") is False, "[B2] Response must have success=False."

    def test_rejects_error_job(self, ps, clean_db):
        """
        An 'error' job must not accept a selection confirmation.
        """
        _seed_job(job_id=202, status="error")

        client = ps.app.test_client()
        resp = client.post("/jobs/202/selection", json={"files": []})

        assert resp.status_code == 409, (
            f"[B2] Error job: expected HTTP 409, got {resp.status_code}."
        )

    def test_rejects_done_job(self, ps, clean_db):
        """
        A 'done' job must not accept a selection confirmation.
        """
        _seed_job(job_id=203, status="done")

        client = ps.app.test_client()
        resp = client.post("/jobs/203/selection", json={"files": []})

        assert resp.status_code == 409, (
            f"[B2] Done job: expected HTTP 409, got {resp.status_code}."
        )

    def test_accepts_selection_pending_job(self, ps, clean_db):
        """
        A 'selection_pending' job must accept the confirmation.
        """
        _seed_job(job_id=204, status="selection_pending")

        client = ps.app.test_client()
        resp = client.post(
            "/jobs/204/selection",
            json={"files": ["/tmp/frame1.fits", "/tmp/frame2.fits"]},
        )

        assert resp.status_code == 200, (
            f"[B2 regression] Pending job: expected HTTP 200, got {resp.status_code}."
        )
        data = resp.get_json()
        assert data.get("success") is True
        assert data.get("files_kept") == 2

    def test_returns_404_for_unknown_job(self, ps, clean_db):
        """
        A non-existent job must return 404, not 200 or 500.
        """
        client = ps.app.test_client()
        resp = client.post("/jobs/99999/selection", json={"files": []})

        assert resp.status_code == 404, (
            f"[B2] Non-existent job: expected HTTP 404, got {resp.status_code}."
        )

    def test_does_not_overwrite_selection_result_on_wrong_status(self, ps, clean_db):
        """
        After a rejected call, the DB row must be unchanged (no selection_result
        written for a job that was not in selection_pending state).
        """
        _seed_job(job_id=205, status="calibrating")

        client = ps.app.test_client()
        client.post("/jobs/205/selection", json={"files": ["/tmp/bad.fits"]})

        stored = _get_selection_result(205)
        assert stored is None, (
            f"[B2] selection_result must not be written for non-pending job, "
            f"but got: {stored}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# B3 – /resume bypasses _active_jobs deduplication
# ─────────────────────────────────────────────────────────────────────────────

class TestB3ResumeDuplicate:
    """
    BUG
    ---
    POST /resume calls _job_queue.put() directly instead of _enqueue_job().
    This bypasses the _active_jobs set, so two rapid POST /resume calls for the
    same result_id both succeed (200) and both tasks run concurrently.

    EXPECTED BEHAVIOUR (after fix)
    --------------------------------
    • First POST /resume → 200
    • Second POST /resume with same result_id → 409 (already active)
    • Calling with a different result_id while first is active → 200
    """

    def test_second_resume_same_result_returns_409(self, ps, clean_db):
        """
        Two consecutive /resume calls for the same result_id must not both
        succeed. The second must return 409.
        Current (bug): both return 200, two workers are spawned.
        """
        # Prevent the worker thread from actually processing the task
        # by using a mock queue.put (the task lambda is never executed).
        original_put = ps._job_queue.put
        ps._job_queue.put = MagicMock()

        try:
            client = ps.app.test_client()

            resp1 = client.post(
                "/resume",
                json={"result_id": "301", "from_step": "photometry"},
            )
            resp2 = client.post(
                "/resume",
                json={"result_id": "301", "from_step": "photometry"},
            )
        finally:
            ps._job_queue.put = original_put
            # Clean up any tracking state added by the fix
            if hasattr(ps, "_active_resumes"):
                ps._active_resumes.discard(301)

        assert resp1.status_code == 200, (
            f"[B3] First /resume should succeed with 200, got {resp1.status_code}."
        )
        assert resp2.status_code == 409, (
            f"[B3] Second /resume for same result_id should return 409, "
            f"got {resp2.status_code}. "
            "Fix: maintain _active_resumes set; use it in POST /resume the same "
            "way _active_jobs is used in _enqueue_job()."
        )

    def test_different_result_ids_both_succeed(self, ps, clean_db):
        """
        Regression: two different result_ids must both be accepted.
        """
        original_put = ps._job_queue.put
        ps._job_queue.put = MagicMock()

        try:
            client = ps.app.test_client()
            resp_a = client.post("/resume", json={"result_id": "401", "from_step": "photometry"})
            resp_b = client.post("/resume", json={"result_id": "402", "from_step": "photometry"})
        finally:
            ps._job_queue.put = original_put
            if hasattr(ps, "_active_resumes"):
                ps._active_resumes.discard(401)
                ps._active_resumes.discard(402)

        assert resp_a.status_code == 200, "[B3 regression] First result_id must get 200."
        assert resp_b.status_code == 200, "[B3 regression] Different result_id must also get 200."

    def test_resume_rejects_invalid_step(self, ps, clean_db):
        """
        /resume must still validate from_step (pre-existing behaviour must not
        be broken by the deduplication fix).
        """
        client = ps.app.test_client()
        resp = client.post(
            "/resume",
            json={"result_id": "501", "from_step": "nonexistent_step"},
        )
        assert resp.status_code == 400, (
            f"[B3 regression] Invalid step should give 400, got {resp.status_code}."
        )
        assert resp.get_json().get("success") is False
