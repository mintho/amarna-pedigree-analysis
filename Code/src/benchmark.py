#!/usr/bin/env python3
"""Benchmark exact likelihood engines on a run configuration."""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

try:
    from .compare import evaluate_with_engine
except ImportError:
    from compare import evaluate_with_engine  # type: ignore[no-redef]


def timed_engine(run_config: Path, engine: str, repeats: int) -> dict[str, Any]:
    durations: list[float] = []
    last_results: list[dict[str, Any]] = []

    for _ in range(repeats):
        start = time.perf_counter()
        last_results = evaluate_with_engine(run_config, engine)
        durations.append(time.perf_counter() - start)

    return {
        "engine": engine,
        "repeats": repeats,
        "best_seconds": min(durations),
        "mean_seconds": sum(durations) / len(durations),
        "worst_seconds": max(durations),
        "top_log_likelihood": last_results[0]["log_likelihood"] if last_results else None,
        "models": len(last_results),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark likelihood engines")
    parser.add_argument(
        "run_config",
        nargs="?",
        default="Data/models/runs/focused_yuya_thuya_tiye_trio_l1_l2_run.json",
        help="Path to a run configuration JSON file",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Number of repeats per engine",
    )
    parser.add_argument(
        "--engine",
        choices=["both", "factors", "enumeration"],
        default="both",
        help="Which engine to benchmark",
    )
    args = parser.parse_args()

    if args.repeats < 1:
        print("BENCHMARK FAILED: --repeats must be >= 1")
        return 1

    run_config = Path(args.run_config).resolve()
    engines = ["factors", "enumeration"] if args.engine == "both" else [args.engine]

    print("=== Likelihood Engine Benchmark ===")
    print(f"Run: {run_config}")
    print(f"Repeats: {args.repeats}")
    print("")
    print(f"{'Engine':<12} {'Best s':>10} {'Mean s':>10} {'Worst s':>10} {'Top logL':>14} {'Models':>6}")

    rows = []
    for engine in engines:
        try:
            row = timed_engine(run_config, engine, args.repeats)
        except Exception as exc:
            print(f"{engine:<12} FAILED: {exc}")
            return 1
        rows.append(row)
        print(
            f"{row['engine']:<12} "
            f"{row['best_seconds']:>10.6f} "
            f"{row['mean_seconds']:>10.6f} "
            f"{row['worst_seconds']:>10.6f} "
            f"{row['top_log_likelihood']:>14.6f} "
            f"{row['models']:>6}"
        )

    if len(rows) == 2:
        by_engine = {row["engine"]: row for row in rows}
        enum = by_engine["enumeration"]["mean_seconds"]
        factors = by_engine["factors"]["mean_seconds"]
        speedup = enum / factors if factors > 0 else float("inf")
        delta = abs(
            by_engine["enumeration"]["top_log_likelihood"]
            - by_engine["factors"]["top_log_likelihood"]
        )
        print("")
        print(f"Mean speedup factors vs enumeration: {speedup:.2f}x")
        print(f"Top logL absolute difference: {delta:.12g}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
