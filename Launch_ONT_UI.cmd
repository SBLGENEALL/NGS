@echo off
setlocal
set "WSL_PROJECT=~/NGS_ONT_batch"
set "UI_PORT=8502"

wsl.exe bash -lc "cd %WSL_PROJECT% && chmod +x launch_ui.sh run_ui.sh && ONT_UI_PORT=%UI_PORT% ./launch_ui.sh --background"
if errorlevel 1 (
    echo.
    echo Failed to start ONT Plasmid Analyzer.
    echo Check the WSL_PROJECT path and launch environment.
    pause
    exit /b 1
)

timeout /t 3 /nobreak >nul
start "" "http://localhost:%UI_PORT%"
endlocal
