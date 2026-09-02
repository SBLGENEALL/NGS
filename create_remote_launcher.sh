#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SERVER_HOST="${1:-}"
OUTPUT="$SCRIPT_DIR/Launch_ONT_UI_Remote.cmd"

if [[ -z "$SERVER_HOST" ]]; then
    echo "Usage: $0 <Linux-server-IP>" >&2
    exit 1
fi
if [[ ! "$SERVER_HOST" =~ ^[0-9A-Za-z.-]+$ ]]; then
    echo "Invalid server address: $SERVER_HOST" >&2
    exit 1
fi

cat > "$OUTPUT" <<EOF
@echo off
setlocal
set "SERVER_HOST=$SERVER_HOST"
set "SERVER_USER=MCET03"
set "REMOTE_PROJECT=/data/user/MCET03/04_ONT/NGS_ONT"
set "REMOTE_ENV=/home/MCET03/conda_envs/NGS_ONT_env"
set "UI_PORT=8502"

start "ONT Server Connection - keep this window open" cmd /k ^
ssh -L %UI_PORT%:127.0.0.1:%UI_PORT% %SERVER_USER%@%SERVER_HOST% ^
"cd %REMOTE_PROJECT% && chmod +x launch_ui.sh run_ui.sh && ONT_CONDA_ENV=%REMOTE_ENV% ONT_UI_ADDRESS=127.0.0.1 ONT_UI_PORT=%UI_PORT% ./launch_ui.sh --background; tail -f /dev/null"

timeout /t 5 /nobreak >nul
start "" "http://localhost:%UI_PORT%"
endlocal
EOF

echo "Created: $OUTPUT"
echo "Download this file to the Windows desktop and double-click it."
