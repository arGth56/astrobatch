#!/usr/bin/env bash
# Beta server — port 3001, isolated DB, separate processing service on port 5201
set -e
cd "$(dirname "$0")"

# Kill any existing beta instances
pkill -f "PROCESSING_PORT=5201" 2>/dev/null || true
pkill -f "PORT=3001" 2>/dev/null || true
sleep 1

# Start beta processing service on port 5201
PROCESSING_PORT=5201 nohup .venv/bin/python astrobatch/processing_service.py \
  >> /tmp/proc_svc_beta.log 2>&1 &
echo "Beta processing service PID=$!"

sleep 2

# Start beta Node server on port 3001
PORT=3001 \
DB_FILE=nightmanager_beta.db \
PROCESSING_SERVICE_URL=http://127.0.0.1:5201 \
nohup node server/server.js >> /tmp/node_server_beta.log 2>&1 &
echo "Beta server PID=$!"

echo "Beta env starting at http://localhost:3001"
