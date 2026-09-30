#!/usr/bin/env python3
"""Aggregate sensitivity matrices into a stability report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def structural_model_key(row: dict[str, Any]) -> str:
    return f"{row['pedigree']} / {row['identity_hypothesis']}"


def matrix_point_best(point: dict[str, Any]) -> dict[str, Any]:
    if not point.get("results"):
        raise ValueError(f"{point.get('matrix_key', 'matrix point')}: no results")
    best = point["results"][0]
    return {
        "model": structural_model_key(best),
        "weight_percent": best["weight_percent"],
        "log_likelihood": best["log_likelihood"],
    }


def changed_axes(point: dict[str, Any], matrix: dict[str, Any]) -> list[str]:
    axes = []
    if point["founder_scenario_id"] != matrix["baseline_founder_scenario_id"]:
        axes.append("founder_scenario")
    if point["background_error"] != matrix["baseline_background_error"]:
        axes.append("background_error")
    if point["dropout"] != matrix["baseline_dropout"]:
        axes.append("dropout")
    return axes


def summarize_matrix(path: Path) -> dict[str, Any]:
    matrix = load_json(path)
    baseline_best = structural_model_key(matrix["baseline"][0])
    baseline_best_weight = matrix["baseline"][0]["weight_percent"]
    points = matrix.get("matrix_points", [])
    if not points:
        raise ValueError(f"{path}: matrix_points must not be empty")

    best_counts: dict[str, int] = {}
    rank_flip_points = []
    max_abs_weight_change = 0.0
    for point in points:
        best = matrix_point_best(point)
        best_counts[best["model"]] = best_counts.get(best["model"], 0) + 1
        max_abs_weight_change = max(
            max_abs_weight_change,
            max(
                (
                    abs(change["weight_percent_change"])
                    for change in point.get("changes_vs_baseline", [])
                ),
                default=0.0,
            ),
        )
        if best["model"] != baseline_best:
            rank_flip_points.append(
                {
                    "matrix_key": point["matrix_key"],
                    "founder_scenario_id": point["founder_scenario_id"],
                    "background_error": point["background_error"],
                    "dropout": point["dropout"],
                    "best_model": best["model"],
                    "best_weight_percent": best["weight_percent"],
                    "changed_axes": changed_axes(point, matrix),
                }
            )

    total_points = len(points)
    best_fractions = {
        model: count / total_points
        for model, count in sorted(best_counts.items(), key=lambda item: item[0])
    }
    baseline_best_count = best_counts.get(baseline_best, 0)

    return {
        "matrix_path": str(path),
        "run_config": matrix["run_config"],
        "engine": matrix["engine"],
        "loci": matrix["loci"],
        "matrix_dimensions": {
            "founder_scenarios": len(matrix["founder_scenarios"]),
            "background_error_values": len(matrix["background_error_values"]),
            "dropout_values": len(matrix["dropout_values"]),
            "total_points": total_points,
        },
        "baseline_best_model": baseline_best,
        "baseline_best_weight_percent": baseline_best_weight,
        "baseline_best_count": baseline_best_count,
        "baseline_best_fraction": baseline_best_count / total_points,
        "rank_flip_count": len(rank_flip_points),
        "rank_flip_fraction": len(rank_flip_points) / total_points,
        "best_model_fractions": best_fractions,
        "max_abs_weight_percent_change": max_abs_weight_change,
        "rank_flip_points": rank_flip_points,
        "summary": matrix["summary"],
    }


def build_report(paths: list[Path]) -> dict[str, Any]:
    runs = [summarize_matrix(path) for path in paths]
    return {
        "schema_version": 1,
        "matrix_files": [str(path) for path in paths],
        "runs": runs,
    }


def markdown_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Stability Report",
        "",
        "## Run Summary",
        "",
        "| Run | Baseline best | Baseline best in matrix | Rank flips | Max weight change |",
        "|---|---|---:|---:|---:|",
    ]
    for run in payload["runs"]:
        lines.append(
            f"| `{Path(run['run_config']).name}` | "
            f"`{run['baseline_best_model']}` | "
            f"{run['baseline_best_fraction'] * 100.0:.1f}% | "
            f"{run['rank_flip_count']}/{run['matrix_dimensions']['total_points']} | "
            f"{run['max_abs_weight_percent_change']:.3f} pp |"
        )

    lines.extend(["", "## Best-Model Fractions", ""])
    for run in payload["runs"]:
        lines.extend(
            [
                f"### `{Path(run['run_config']).name}`",
                "",
                "| Model | Matrix share as rank 1 |",
                "|---|---:|",
            ]
        )
        for model, fraction in run["best_model_fractions"].items():
            lines.append(f"| `{model}` | {fraction * 100.0:.1f}% |")
        lines.append("")

    lines.extend(["", "## Rank-Flip Points", ""])
    for run in payload["runs"]:
        lines.extend(
            [
                f"### `{Path(run['run_config']).name}`",
                "",
                "| Founder scenario | Background error | Dropout | New rank 1 model | Rank 1 weight | Changed axes |",
                "|---|---:|---:|---|---:|---|",
            ]
        )
        if not run["rank_flip_points"]:
            lines.append("| n/a | n/a | n/a | n/a | n/a | n/a |")
        for point in run["rank_flip_points"]:
            lines.append(
                f"| `{point['founder_scenario_id']}` | "
                f"{point['background_error']:.6g} | "
                f"{point['dropout']:.6g} | "
                f"`{point['best_model']}` | "
                f"{point['best_weight_percent']:.3f}% | "
                f"{', '.join(point['changed_axes']) or 'baseline'} |"
            )
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_report(output_dir: Path, payload: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "stability_report.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (output_dir / "stability_report.md").write_text(
        markdown_report(payload),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize sensitivity matrices into a stability report"
    )
    parser.add_argument("matrix_json", nargs="+")
    parser.add_argument(
        "--out_dir",
        default="reproduced/stability_report",
    )
    args = parser.parse_args()

    try:
        payload = build_report([Path(path) for path in args.matrix_json])
        write_report(Path(args.out_dir), payload)
    except Exception as exc:
        print(f"STABILITY REPORT FAILED: {exc}")
        return 1

    print("=== Stability Report ===")
    print(f"Results: {Path(args.out_dir) / 'stability_report.json'}")
    print(f"Report: {Path(args.out_dir) / 'stability_report.md'}")
    for run in payload["runs"]:
        print(
            f"{Path(run['run_config']).name}: "
            f"baseline best in "
            f"{run['baseline_best_fraction'] * 100.0:.1f}% of matrix, "
            f"{run['rank_flip_count']} rank flips"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
