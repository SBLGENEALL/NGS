#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${ONT_UI_PORT:-8502}"
ADDRESS="${ONT_UI_ADDRESS:-0.0.0.0}"
MODE="${1:-foreground}"
SERVER_ROOT="${ONT_SERVER_ROOT:-/data}"
SERVER_START="${ONT_SERVER_START:-$SERVER_ROOT}"
LOG_ROOT="${ONT_LOG_ROOT:-$SCRIPT_DIR/usage_logs}"
UI_RUN_ROOT="${ONT_UI_RUN_ROOT:-$SCRIPT_DIR/ui_runs}"

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

stop_running_ui() {
    local pid=""
    if [[ -f "$SCRIPT_DIR/.ont_ui.pid" ]]; then
        pid="$(tr -cd '0-9' < "$SCRIPT_DIR/.ont_ui.pid")"
        if [[ -n "$pid" ]] && ! kill -0 "$pid" 2>/dev/null; then
            pid=""
        fi
    fi
    if [[ -z "$pid" ]] && command -v ss >/dev/null 2>&1; then
        pid="$(
            ss -ltnp "sport = :$PORT" 2>/dev/null \
                | sed -n 's/.*pid=\([0-9]*\).*/\1/p' \
                | head -1
        )"
    fi
    if [[ -n "$pid" ]]; then
        echo "Stopping ONT Plasmid Analyzer (PID $pid)..."
        kill "$pid"
        for _ in {1..20}; do
            port_is_open || break
            sleep 0.25
        done
    elif port_is_open; then
        echo "Port $PORT is in use, but its process could not be identified." >&2
        return 1
    fi
    rm -f "$SCRIPT_DIR/.ont_ui.pid"
}

ENV_PATH="$(find_environment)"
export PATH="$ENV_PATH/bin:$PATH"
export CONDA_PREFIX="$ENV_PATH"

if [[ "$MODE" == "--restart" ]]; then
    stop_running_ui
    MODE="--background"
fi

if port_is_open; then
    echo "ONT Plasmid Analyzer is already running: http://127.0.0.1:$PORT"
    open_browser
    exit 0
fi

if [[ "$MODE" == "--background" ]]; then
    nohup env ONT_UI_ADDRESS="$ADDRESS" ONT_UI_PORT="$PORT" \
        ONT_SERVER_ROOT="$SERVER_ROOT" ONT_SERVER_START="$SERVER_START" \
        ONT_LOG_ROOT="$LOG_ROOT" ONT_UI_RUN_ROOT="$UI_RUN_ROOT" \
        "$SCRIPT_DIR/run_ui.sh" > "$SCRIPT_DIR/ui_launcher.log" 2>&1 &
    echo "$!" > "$SCRIPT_DIR/.ont_ui.pid"
    echo "Starting ONT Plasmid Analyzer: http://127.0.0.1:$PORT"
    exit 0
fi

( sleep 3; open_browser ) &
exec env ONT_UI_ADDRESS="$ADDRESS" ONT_UI_PORT="$PORT" \
    ONT_SERVER_ROOT="$SERVER_ROOT" ONT_SERVER_START="$SERVER_START" \
    ONT_LOG_ROOT="$LOG_ROOT" ONT_UI_RUN_ROOT="$UI_RUN_ROOT" \
    "$SCRIPT_DIR/run_ui.sh"
