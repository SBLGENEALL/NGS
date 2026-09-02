@echo off
setlocal EnableExtensions
title ONT Plasmid Analyzer

set "SERVER_HOST=182.198.164.21"
set "SERVER_USER=MCET03"
set "LOCAL_PORT=18502"
set "REMOTE_PORT=8502"
set "REMOTE_SCRIPT=/data/user/MCET03/04_ONT/NGS_ONT/ont_one_click.sh"
set "UI_URL=http://127.0.0.1:%LOCAL_PORT%"
set "MOBA_PATH="

rem If the tunnel and UI are already healthy, do not touch MobaXterm.
call :ui_ready
if not errorlevel 1 goto :open_browser

call :find_mobaxterm
if not defined MOBA_PATH (
    echo.
    echo MobaXterm.exe could not be found automatically.
    echo Put this BAT file beside MobaXterm.exe, or install MobaXterm normally.
    echo.
    pause
    exit /b 1
)

rem -newtab reuses a running MobaXterm instance, or starts one if needed.
start "" "%MOBA_PATH%" -newtab "ssh -o ExitOnForwardFailure=yes -L %LOCAL_PORT%:127.0.0.1:%REMOTE_PORT% %SERVER_USER%@%SERVER_HOST% 'bash %REMOTE_SCRIPT% %USERNAME%; exec bash -l'"

rem Give MobaXterm a moment, then wait only until Streamlit actually responds.
ping 127.0.0.1 -n 1 -w 500 >nul
powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "$deadline=(Get-Date).AddSeconds(120); do { try { $r=Invoke-WebRequest -UseBasicParsing -Uri '%UI_URL%' -TimeoutSec 1; if ($r.StatusCode -eq 200) { exit 0 } } catch {}; Start-Sleep -Milliseconds 300 } while ((Get-Date) -lt $deadline); exit 1"

if errorlevel 1 (
    echo.
    echo The SSH window opened, but the UI did not respond within 120 seconds.
    echo Complete the SSH password or host-key prompt in MobaXterm, then run this file again.
    echo.
    pause
    exit /b 1
)

:open_browser
start "" "%UI_URL%"
exit /b 0

:ui_ready
powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "try { $r=Invoke-WebRequest -UseBasicParsing -Uri '%UI_URL%' -TimeoutSec 1; if ($r.StatusCode -eq 200) { exit 0 } } catch {}; exit 1" >nul 2>&1
exit /b %errorlevel%

:find_mobaxterm
for /f "delims=" %%I in ('where MobaXterm.exe 2^>nul') do if not defined MOBA_PATH set "MOBA_PATH=%%~fI"
for %%I in (
    "%~dp0MobaXterm.exe"
    "%~dp0MobaXterm_Personal*.exe"
    "%ProgramFiles%\Mobatek\MobaXterm\MobaXterm*.exe"
    "%ProgramFiles%\MobaXterm\MobaXterm*.exe"
    "%ProgramFiles(x86)%\Mobatek\MobaXterm\MobaXterm*.exe"
    "%LOCALAPPDATA%\Programs\MobaXterm\MobaXterm*.exe"
    "%USERPROFILE%\Desktop\MobaXterm*.exe"
    "%USERPROFILE%\Downloads\MobaXterm*.exe"
) do if not defined MOBA_PATH if exist "%%~fI" set "MOBA_PATH=%%~fI"
if defined MOBA_PATH exit /b 0

for %%R in ("%ProgramFiles%\Mobatek" "%LOCALAPPDATA%\Programs" "%USERPROFILE%\Desktop" "%USERPROFILE%\Downloads") do (
    if exist "%%~R" for /f "delims=" %%I in ('where /r "%%~R" MobaXterm*.exe 2^>nul') do if not defined MOBA_PATH set "MOBA_PATH=%%~fI"
)
exit /b 0
