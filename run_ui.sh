#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ADDRESS="${ONT_UI_ADDRESS:-0.0.0.0}"
PORT="${ONT_UI_PORT:-8501}"
MAX_UPLOAD_MB="${ONT_UI_MAX_UPLOAD_MB:-102400}"

if ! command -v streamlit >/dev/null 2>&1; then
    echo "streamlit was not found. Activate/update the NGS_env conda environment." >&2
    exit 1
fi

exec streamlit run "$SCRIPT_DIR/app.py" \
    --server.address "$ADDRESS" \
    --server.port "$PORT" \
    --server.maxUploadSize "$MAX_UPLOAD_MB" \
    --server.headless true \
    --server.fileWatcherType none \
    --server.runOnSave false \
    --browser.gatherUsageStats false
