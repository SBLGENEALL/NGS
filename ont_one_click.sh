#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${ONT_PROJECT_DIR:-$SCRIPT_DIR}"
CONDA_BASE="${CONDA_BASE:-/home/mcet/anaconda3}"
CONDA_ENV="${ONT_CONDA_ENV:-$HOME/conda_envs/NGS_ONT_env}"
UI_ADDRESS="${ONT_UI_ADDRESS:-127.0.0.1}"
UI_PORT="${ONT_UI_PORT:-8502}"
SERVER_ROOT="${ONT_SERVER_ROOT:-/data}"
SERVER_START="${ONT_SERVER_START:-$SERVER_ROOT}"
LOG_ROOT="${ONT_LOG_ROOT:-$PROJECT_DIR/usage_logs}"
LAUNCH_USER="$(printf '%s' "${*:-unknown}" | tr -cd '[:alnum:]_.@ -')"
LAUNCHED_AT="$(TZ=Asia/Seoul date '+%Y-%m-%dT%H:%M:%S%z')"

if [[ ! -d "$PROJECT_DIR" ]]; then
    echo "ONT project folder was not found: $PROJECT_DIR" >&2
    exit 1
fi
if [[ ! -x "$CONDA_ENV/bin/streamlit" ]]; then
    echo "NGS_ONT_env does not contain Streamlit: $CONDA_ENV" >&2
    exit 1
fi

if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1091
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
else
    export PATH="$CONDA_ENV/bin:$PATH"
    export CONDA_PREFIX="$CONDA_ENV"
fi

cd "$PROJECT_DIR"
chmod +x launch_ui.sh run_ui.sh run_pipeline.sh
mkdir -p "$LOG_ROOT"
printf '%s\tUI_LAUNCH\t%s\n' "$LAUNCHED_AT" "${LAUNCH_USER:-unknown}" \
    >> "$LOG_ROOT/ui_access.log"

ONT_CONDA_ENV="$CONDA_ENV" \
ONT_UI_ADDRESS="$UI_ADDRESS" \
ONT_UI_PORT="$UI_PORT" \
ONT_SERVER_ROOT="$SERVER_ROOT" \
ONT_SERVER_START="$SERVER_START" \
ONT_LOG_ROOT="$LOG_ROOT" \
./launch_ui.sh --restart

echo "ONT Plasmid Analyzer is ready on the SSH tunnel port."
