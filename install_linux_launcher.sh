#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$HOME/.local/share/applications"
DESKTOP_FILE="$APP_DIR/ont-plasmid-analyzer.desktop"
ICON="applications-science"

for candidate in "$SCRIPT_DIR/branding_sidebar_logo.png" "$SCRIPT_DIR/branding_logo.png"; do
    if [[ -f "$candidate" ]]; then
        ICON="$candidate"
        break
    fi
done

mkdir -p "$APP_DIR"
cat > "$DESKTOP_FILE" <<EOF
[Desktop Entry]
Type=Application
Name=ONT Plasmid Analyzer
Comment=Reference-based ONT plasmid variant analysis
Exec=$SCRIPT_DIR/launch_ui.sh
Icon=$ICON
Terminal=true
Categories=Science;Biology;
StartupNotify=true
EOF
chmod +x "$DESKTOP_FILE" "$SCRIPT_DIR/launch_ui.sh" "$SCRIPT_DIR/run_ui.sh"

if [[ -d "$HOME/Desktop" ]]; then
    cp "$DESKTOP_FILE" "$HOME/Desktop/ONT_Plasmid_Analyzer.desktop"
    chmod +x "$HOME/Desktop/ONT_Plasmid_Analyzer.desktop"
    echo "Desktop icon created: $HOME/Desktop/ONT_Plasmid_Analyzer.desktop"
else
    echo "Application launcher created: $DESKTOP_FILE"
fi
