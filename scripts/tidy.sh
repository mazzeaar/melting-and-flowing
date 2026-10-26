#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <build-directory>" >&2
    exit 1
fi

BUILD_DIR="$(cd "$ROOT" && realpath "$1")"

command -v clang-tidy >/dev/null || {
    echo "Error: clang-tidy not found" >&2
    exit 1
}

if [[ ! -f "$BUILD_DIR/compile_commands.json" ]]; then
    echo "Error: missing $BUILD_DIR/compile_commands.json" >&2
    exit 1
fi

find "$ROOT/src" -type f \( \
    -name '*.cpp' -o \
    -name '*.cc' -o \
    -name '*.c' \
    \) -print0 |
    xargs -0 clang-tidy \
        -p "$BUILD_DIR" \
        --config-file="$ROOT/.clang-tidy" \
        --header-filter="^$ROOT/src/"

echo "Clang-tidy complete."
