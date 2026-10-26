#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

command -v clang-format >/dev/null || {
    echo "Error: clang-format not found" >&2
    exit 1
}

find "$ROOT/src" -type f \( \
    -name '*.cpp' -o \
    -name '*.cc' -o \
    -name '*.c' -o \
    -name '*.hpp' -o \
    -name '*.hh' -o \
    -name '*.h' \
    \) -print0 |
    xargs -0 clang-format -i --style=file

echo "Formatting complete."
