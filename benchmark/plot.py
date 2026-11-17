#!/usr/bin/env python3
"""Plot local solver benchmark history.

Usage: python benchmark/plot.py basic
Requires: matplotlib (python -m pip install matplotlib)
Reads: benchmark/results/basic/*.csv
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Works without a display on Windows, macOS and Linux.
import matplotlib.pyplot as plt


RESULTS = Path(__file__).resolve().parent / "results"


def load(path: Path):
    with path.open(newline="", encoding="utf-8") as file:
        rows = list(csv.DictReader(file))
    if not rows or "step" not in rows[0] or "total_ms" not in rows[0]:
        raise ValueError(f"Invalid benchmark CSV: {path}")
    return {
        key: [float(row[key]) for row in rows]
        for key in rows[0]
        if key and key != "step"
    } | {"step": [int(row["step"]) for row in rows]}


def main():
    parser = argparse.ArgumentParser(description="Plot solver benchmark results")
    parser.add_argument("config", help="Benchmark name, e.g. basic")
    args = parser.parse_args()

    directory = RESULTS / args.config
    files = sorted(directory.glob("*.csv"))
    if not files:
        parser.error(f"No CSV results found in {directory}")

    runs = [(path, load(path)) for path in files]
    latest_file, latest = runs[-1]
    previous = runs[-2][1] if len(runs) > 1 else None

    # Plot each recorded operation for latest run, with previous run overlaid.
    metrics = [key for key in latest if key.endswith("_ms")]
    fig, axes = plt.subplots(len(metrics), 1, figsize=(11, max(3, len(metrics) * 2.3)),
                             sharex=True, squeeze=False)
    for ax, key in zip(axes[:, 0], metrics):
        ax.plot(latest["step"], latest[key], label="Latest")
        if previous is not None and key in previous:
            ax.plot(runs[-2][1]["step"], previous[key], "--", label="Previous")
        ax.set_ylabel(f"{key}\nms")
        ax.grid(alpha=0.3)
        ax.legend()
    axes[-1, 0].set_xlabel("Solver step")
    fig.suptitle(f"{args.config}: latest vs previous")
    fig.tight_layout()
    comparison = directory / "comparison.png"
    fig.savefig(comparison, dpi=150)
    plt.close(fig)

    # Track mean full-step duration across every saved benchmark run.
    labels = [p.stem for p, _ in runs]
    means = [sum(data["total_ms"]) / len(data["total_ms"]) for _, data in runs]
    fig, ax = plt.subplots(figsize=(max(8, len(runs) * 0.7), 4))
    ax.plot(range(len(runs)), means, marker="o")
    ax.set_xticks(range(len(runs)), labels, rotation=45, ha="right")
    ax.set_ylabel("Mean solver step (ms)")
    ax.set_title(f"{args.config}: benchmark history")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    history = directory / "history.png"
    fig.savefig(history, dpi=150)
    plt.close(fig)

    print(f"Latest: {latest_file.name} ({len(runs)} saved runs)")
    print(f"Saved {comparison} and {history}")


if __name__ == "__main__":
    main()
