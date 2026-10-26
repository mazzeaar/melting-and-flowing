#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$ROOT/scripts/tidy.sh" "$@"
"$ROOT/scripts/format.sh"

echo "All checks complete."
