#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ORGANIZATION="${1:-}"
DISTRIBUTOR="${2:-Jongin Baek}"

if [[ -z "$ORGANIZATION" ]]; then
    echo "Usage: $0 <organization> [distributed-by]" >&2
    exit 1
fi

python3 - "$SCRIPT_DIR/branding.local.json" "$ORGANIZATION" "$DISTRIBUTOR" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
values = {
    "organization": sys.argv[2],
    "distributed_by": sys.argv[3],
}
path.write_text(
    json.dumps(values, ensure_ascii=False, indent=2) + "\n",
    encoding="utf-8",
)
print(f"Created local branding: {path}")
PY
