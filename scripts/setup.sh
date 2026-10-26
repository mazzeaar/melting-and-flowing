#!/bin/sh

ensure_brew_installed()
{
if command -v brew >/dev/null 2>&1; then
    echo 0
else
    echo "brew not installed, please install brew first"
    return 1
fi
}

ensure_ninja_installed() {
  if command -v ninja >/dev/null 2>&1; then
    return 0
  else
    echo "Installing ninja..."
    brew install ninja
    return 1
  fi
}

init_submodules() {
    echo "Initializing git submodules..."
    git submodule update --init --recursive
}

if [ $(ensure_brew_installed) -eq 0 ]; then
    init_submodules
    ensure_ninja_installed
    echo "[Setup Done]\n"
    echo "Start running the application using"
    echo "cmake --preset [release|debug]"
    echo "ninja -C build/release && ./build/release/fluid_sim"
fi
