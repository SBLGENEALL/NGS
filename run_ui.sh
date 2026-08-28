#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ADDRESS="${ONT_UI_ADDRESS:-0.0.0.0}"
PORT="${ONT_UI_PORT:-8501}"

if ! command -v streamlit >/dev/null 2>&1; then
    echo "streamlit was not found. Activate/update the NGS_env conda environment." >&2
    exit 1
fi

exec streamlit run "$SCRIPT_DIR/app.py" \
    --server.address "$ADDRESS" \
    --server.port "$PORT" \
    --server.headless true \
    --browser.gatherUsageStats false
