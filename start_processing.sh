#!/bin/bash
# Start the Python processing service (port 5200)
# Usage: ./start_processing.sh
# Or: source /home/pyl/Documents/astrobatch/.venv/bin/activate && python3 -m astrobatch.processing_service

set -e
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV="$PROJECT_DIR/.venv"

export NIGHTMANAGER_DB="$PROJECT_DIR/server/nightmanager.db"
export STDWEB_URL="${STDWEB_URL:-http://86.253.141.183:7000}"
export STDWEB_TOKEN="${STDWEB_TOKEN:-1e296ddd6738af45467b7bc6558c00a9524447ab}"
export SIRIL_BIN="${SIRIL_BIN:-siril}"
export PROCESSING_PORT="${PROCESSING_PORT:-5200}"

# Load .env from server/ if present
ENV_FILE="$PROJECT_DIR/server/.env"
if [ -f "$ENV_FILE" ]; then
  set -o allexport
  # shellcheck disable=SC1090
  source "$ENV_FILE"
  set +o allexport
fi

exec "$VENV/bin/python3" -m astrobatch.processing_service
