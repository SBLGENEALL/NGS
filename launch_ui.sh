#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${ONT_UI_PORT:-8502}"
ADDRESS="${ONT_UI_ADDRESS:-0.0.0.0}"
MODE="${1:-foreground}"

find_environment() {
    local candidate
    if [[ -n "${ONT_CONDA_ENV:-}" ]]; then
        candidate="$ONT_CONDA_ENV"
        if [[ -x "$candidate/bin/streamlit" ]]; then
            echo "$candidate"
            return 0
        fi
        echo "ONT_CONDA_ENV does not contain streamlit: $candidate" >&2
        return 1
    fi
    for candidate in \
        /home/MCET03/conda_envs/NGS_ONT_env \
        /home/MCET03/conda_envs/NGS_env \
        "$HOME/conda_envs/NGS_ONT_env" \
        "$HOME/conda_envs/NGS_env" \
        "$HOME/miniconda3/envs/NGS_ONT_env" \
        "$HOME/miniconda3/envs/NGS_env" \
        "$HOME/anaconda3/envs/NGS_ONT_env" \
        "$HOME/anaconda3/envs/NGS_env" \
        "$HOME/miniconda3/envs/ONT_UI_TEST" \
        "$HOME/anaconda3/envs/ONT_UI_TEST"; do
        if [[ -x "$candidate/bin/streamlit" ]]; then
            echo "$candidate"
            return 0
        fi
    done
    if command -v streamlit >/dev/null 2>&1; then
        dirname "$(dirname "$(command -v streamlit)")"
        return 0
    fi
    echo "Could not find NGS_ONT_env, NGS_env, or ONT_UI_TEST with streamlit installed." >&2
    echo "Set the environment once, for example:" >&2
    echo "  ONT_CONDA_ENV=/home/MCET03/conda_envs/NGS_ONT_env ./launch_ui.sh" >&2
    return 1
}

port_is_open() {
    timeout 1 bash -c "</dev/tcp/127.0.0.1/$PORT" >/dev/null 2>&1
}

open_browser() {
    if command -v xdg-open >/dev/null 2>&1; then
        xdg-open "http://127.0.0.1:$PORT" >/dev/null 2>&1 || true
    fi
}

ENV_PATH="$(find_environment)"
export PATH="$ENV_PATH/bin:$PATH"
export CONDA_PREFIX="$ENV_PATH"

if port_is_open; then
    echo "ONT Plasmid Analyzer is already running: http://127.0.0.1:$PORT"
    open_browser
    exit 0
fi

if [[ "$MODE" == "--background" ]]; then
    nohup env ONT_UI_ADDRESS="$ADDRESS" ONT_UI_PORT="$PORT" \
        "$SCRIPT_DIR/run_ui.sh" > "$SCRIPT_DIR/ui_launcher.log" 2>&1 &
    echo "$!" > "$SCRIPT_DIR/.ont_ui.pid"
    echo "Starting ONT Plasmid Analyzer: http://127.0.0.1:$PORT"
    exit 0
fi

( sleep 3; open_browser ) &
exec env ONT_UI_ADDRESS="$ADDRESS" ONT_UI_PORT="$PORT" "$SCRIPT_DIR/run_ui.sh"
