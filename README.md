# Melting and Flowing

Our implementation of the 2002 paper [Melting and Flowing](https://faculty.cc.gatech.edu/~turk/my_papers/melt.pdf).

## Setup
To set up the project you can simply run
```bash
./scripts/setup.sh
```

## todo
- user needs ninja, cmake, some cpp23 compiler
- how to clone/init submodules
- presets, compile, run

## Getting started
'''
cmake --preset release # or debug
ninja -C release && ./build/release/sim
'''

## Benchmark
Compile and run can be done with:
```bash
cmake --build --preset release --target fluid_benchmark
./build/release/benchmark/fluid_benchmark benchmark/configs/{BENCH_TYPE}.toml
```

Plots can be generated like:
```bash
python benchmark/plot.py {BENCH_TYPE}
```

Currently there is only one benchmark: 'basic'

## Integrate with your tools

## Collaborate with your team

