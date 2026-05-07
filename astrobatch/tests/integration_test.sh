#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# integration_test.sh
#
# Layer-2 integration tests: runs against a LIVE processing service.
# Each test issues a real HTTP request and validates the JSON response.
#
# Requires:
#   - Processing service running (default: http://localhost:5200)
#   - sqlite3 CLI available
#   - jq available  (apt install jq)
#
# Usage:
#   # Against the default local service:
#   bash astrobatch/tests/integration_test.sh
#
#   # Against a custom URL:
#   SVC_URL=http://localhost:5200 bash astrobatch/tests/integration_test.sh
#
#   # To keep the test DB rows for inspection afterwards:
#   KEEP_ROWS=1 bash astrobatch/tests/integration_test.sh
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SVC_URL="${SVC_URL:-http://localhost:5200}"
DB_PATH="${NIGHTMANAGER_DB:-$(dirname "$0")/../../nightmanager.db}"
# Use absolute path
DB_PATH="$(realpath "$DB_PATH")"
KEEP_ROWS="${KEEP_ROWS:-0}"

# ── Colours ───────────────────────────────────────────────────────────────────
GREEN="\033[0;32m"; RED="\033[0;31m"; YELLOW="\033[1;33m"; RESET="\033[0m"
PASS=0; FAIL=0; SKIP=0

pass() { echo -e "${GREEN}  ✓ PASS${RESET}  $1"; ((PASS++)); }
fail() { echo -e "${RED}  ✗ FAIL${RESET}  $1"; ((FAIL++)); }
skip() { echo -e "${YELLOW}  ↷ SKIP${RESET}  $1"; ((SKIP++)); }
section() { echo -e "\n${YELLOW}── $1 ──${RESET}"; }

# ── Prerequisites ─────────────────────────────────────────────────────────────
echo "Integration tests → $SVC_URL"
echo "DB: $DB_PATH"

if ! command -v jq &>/dev/null; then
  echo "ERROR: jq not found. Install with: sudo apt install jq"; exit 1
fi

# Check service is up
if ! curl -sf "$SVC_URL/health" >/dev/null 2>&1; then
  echo "ERROR: Processing service not reachable at $SVC_URL/health"
  echo "Start it with:  NIGHTMANAGER_DB=\$DB_PATH python -m astrobatch.processing_service"
  exit 1
fi

# ── Helpers ───────────────────────────────────────────────────────────────────

# Insert a job row directly via sqlite3
db_insert_job() {
  local id="$1" status="$2"
  sqlite3 "$DB_PATH" "INSERT OR REPLACE INTO pipeline_jobs(id, status) VALUES ($id, '$status');"
}

db_get_job_status() {
  sqlite3 "$DB_PATH" "SELECT status FROM pipeline_jobs WHERE id=$1;" 2>/dev/null
}

db_get_selection_result() {
  sqlite3 "$DB_PATH" "SELECT selection_result FROM pipeline_jobs WHERE id=$1;" 2>/dev/null
}

db_cleanup() {
  if [[ "$KEEP_ROWS" == "0" ]]; then
    sqlite3 "$DB_PATH" "DELETE FROM pipeline_jobs WHERE id IN ($1);" 2>/dev/null || true
  fi
}

http_status() {
  # Returns HTTP status code only
  curl -s -o /dev/null -w "%{http_code}" \
    -X POST -H "Content-Type: application/json" \
    -d "$2" "$SVC_URL$1"
}

http_json() {
  # Returns response body
  curl -s -X POST -H "Content-Type: application/json" \
    -d "$2" "$SVC_URL$1"
}

# ─────────────────────────────────────────────────────────────────────────────
section "B2 – confirm_selection status guard"
# ─────────────────────────────────────────────────────────────────────────────
# These tests do NOT require the pipeline to be running — they just check
# the HTTP layer.

# B2-1: 'running' job must be rejected (expected: 409 after fix; 200 = bug still present)
db_insert_job 8001 "running"
code=$(http_status "/jobs/8001/selection" '{"files":[]}')
if [[ "$code" == "409" ]]; then
  pass "B2-1: running job → 409 (guard in place)"
elif [[ "$code" == "200" ]]; then
  fail "B2-1: running job → 200 (BUG: no status guard, fix pending)"
else
  fail "B2-1: running job → unexpected HTTP $code"
fi
db_cleanup "8001"

# B2-2: 'error' job must be rejected
db_insert_job 8002 "error"
code=$(http_status "/jobs/8002/selection" '{"files":[]}')
if [[ "$code" == "409" ]]; then
  pass "B2-2: error job → 409"
elif [[ "$code" == "200" ]]; then
  fail "B2-2: error job → 200 (BUG)"
else
  fail "B2-2: error job → unexpected HTTP $code"
fi
db_cleanup "8002"

# B2-3: 'selection_pending' job must be accepted
db_insert_job 8003 "selection_pending"
code=$(http_status "/jobs/8003/selection" '{"files":["/tmp/a.fits"]}')
if [[ "$code" == "200" ]]; then
  pass "B2-3: selection_pending job → 200 (accepted)"
else
  fail "B2-3: selection_pending job → unexpected HTTP $code (expected 200)"
fi
db_cleanup "8003"

# B2-4: non-existent job must be 404
code=$(http_status "/jobs/99999/selection" '{"files":[]}')
if [[ "$code" == "404" ]]; then
  pass "B2-4: non-existent job → 404"
else
  fail "B2-4: non-existent job → unexpected HTTP $code (expected 404)"
fi

# B2-5: verify DB is NOT written for rejected jobs
db_insert_job 8005 "calibrating"
http_json "/jobs/8005/selection" '{"files":["/tmp/sneaky.fits"]}' >/dev/null
stored=$(db_get_selection_result 8005)
if [[ -z "$stored" ]]; then
  pass "B2-5: no selection_result written to DB for rejected job"
else
  fail "B2-5: selection_result was written for non-pending job: $stored"
fi
db_cleanup "8005"

# ─────────────────────────────────────────────────────────────────────────────
section "B3 – /resume deduplication"
# ─────────────────────────────────────────────────────────────────────────────
# We use a result_id range (9000+) that should never exist in the real DB.
# The requests won't actually process anything meaningful; we only care about
# the HTTP response code indicating whether deduplication is in place.

# B3-1: two identical /resume calls — second must be 409
resp1=$(http_json "/resume" '{"result_id":"9001","from_step":"photometry"}')
code1=$(echo "$resp1" | jq -r '.success // "null"')

resp2=$(http_json "/resume" '{"result_id":"9001","from_step":"photometry"}')
code2_http=$(http_status "/resume" '{"result_id":"9001","from_step":"photometry"}')

# Give the worker a moment to pick up the first task
sleep 1

# Re-issue the second call now (after the worker may have started)
resp2b_http=$(http_status "/resume" '{"result_id":"9001","from_step":"photometry"}')

if [[ "$resp2b_http" == "409" ]]; then
  pass "B3-1: second /resume for same result_id → 409 (dedup working)"
elif [[ "$resp2b_http" == "200" ]]; then
  fail "B3-1: second /resume → 200 (BUG: no deduplication)"
else
  fail "B3-1: second /resume → unexpected HTTP $resp2b_http"
fi

# B3-2: different result_ids must both succeed
resp_a=$(http_status "/resume" '{"result_id":"9101","from_step":"photometry"}')
resp_b=$(http_status "/resume" '{"result_id":"9102","from_step":"photometry"}')
if [[ "$resp_a" == "200" && "$resp_b" == "200" ]]; then
  pass "B3-2: different result_ids both → 200 (no false positives)"
else
  fail "B3-2: different result_ids → $resp_a / $resp_b (expected 200/200)"
fi

# B3-3: invalid step must return 400 (regression — existing validation)
code_bad=$(http_status "/resume" '{"result_id":"9200","from_step":"bogus_step"}')
if [[ "$code_bad" == "400" ]]; then
  pass "B3-3: invalid from_step → 400 (validation intact)"
else
  fail "B3-3: invalid from_step → $code_bad (expected 400)"
fi

# ─────────────────────────────────────────────────────────────────────────────
section "B1 – wait_state ignored  (manual verification only)"
# ─────────────────────────────────────────────────────────────────────────────
# B1 cannot be fully tested via curl because it requires the processing service
# to receive a "failed" response from STDWeb, which we cannot fake without
# modifying the service. Instead we provide a manual verification procedure.
echo ""
echo "  B1 cannot be automated with curl alone."
echo "  Use the pytest unit tests instead:"
echo ""
echo "    cd /home/pyl/Documents/astrobatch"
echo "    .venv/bin/python -m pytest astrobatch/tests/test_processing_bugs.py::TestB1WaitStateIgnored -v"
echo ""
echo "  To verify the bug manually (BEFORE fix):"
echo "    1. Trigger a /resume for a result whose STDWeb task is already in 'failed' state."
echo "    2. Check: sqlite3 nightmanager.db 'SELECT status FROM pipeline_results WHERE id=<id>;'"
echo "    3. Expected (after fix): 'error'. Current (bug): 'done'."
skip "B1: manual/pytest only — see instructions above"

# ─────────────────────────────────────────────────────────────────────────────
section "Summary"
# ─────────────────────────────────────────────────────────────────────────────
TOTAL=$((PASS + FAIL + SKIP))
echo ""
echo "  Passed: $PASS / $TOTAL"
[[ $FAIL -gt 0 ]] && echo -e "  ${RED}Failed: $FAIL / $TOTAL${RESET}"
[[ $SKIP -gt 0 ]] && echo -e "  ${YELLOW}Skipped: $SKIP / $TOTAL${RESET}"
echo ""

if [[ $FAIL -gt 0 ]]; then
  echo "  Tests marked FAIL show bugs that are NOT yet fixed."
  echo "  Tests marked BUG in parentheses are expected to fail before the fix."
  exit 1
fi
exit 0
