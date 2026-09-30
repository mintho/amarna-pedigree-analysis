#!/usr/bin/env python3
"""Classify model-family robustness from a stability report."""

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


def classify_run(run: dict[str, Any]) -> str:
    """Classify robustness as stable, assumption-dependent, or unstable."""
    if run["rank_flip_count"] == 0 and run["baseline_best_fraction"] == 1.0:
        return "stabil"
    if run["baseline_best_fraction"] >= 0.5:
        return "annahmeabhaengig"
    return "instabil"


def main_driver_axis(run: dict[str, Any]) -> str:
    axis_counts: dict[str, int] = {}
    for point in run.get("rank_flip_points", []):
        for axis in point.get("changed_axes", []):
            axis_counts[axis] = axis_counts.get(axis, 0) + 1
    if not axis_counts:
        return "keine"
    return max(axis_counts.items(), key=lambda item: (item[1], item[0]))[0]


def conclusion_text(run: dict[str, Any], classification: str) -> str:
    baseline = run["baseline_best_model"]
    fraction = run["baseline_best_fraction"] * 100.0
    flips = run["rank_flip_count"]
    total = run["matrix_dimensions"]["total_points"]
    if classification == "stabil":
        return (
            f"`{baseline}` bleibt in allen {total} Matrixpunkten Rang 1. "
            "Das Ergebnis ist innerhalb des gerechneten Annahmeraums stabil."
        )
    if classification == "annahmeabhaengig":
        axis = main_driver_axis(run)
        return (
            f"`{baseline}` bleibt in {fraction:.1f}% der Matrixpunkte Rang 1, "
            f"aber {flips} von {total} Punkten zeigen Rangwechsel. "
            f"Die Rangwechsel haengen vor allem an der Achse `{axis}`."
        )
    return (
        f"`{baseline}` bleibt nur in {fraction:.1f}% der Matrixpunkte Rang 1. "
        "Das Ergebnis ist im gerechneten Annahmeraum instabil."
    )


def classify_stability_report(path: Path) -> dict[str, Any]:
    stability = load_json(path)
    conclusions = []
    for run in stability.get("runs", []):
        classification = classify_run(run)
        conclusions.append(
            {
                "run_config": run["run_config"],
                "run_name": Path(run["run_config"]).name,
                "classification": classification,
                "baseline_best_model": run["baseline_best_model"],
                "baseline_best_fraction": run["baseline_best_fraction"],
                "rank_flip_count": run["rank_flip_count"],
                "total_matrix_points": run["matrix_dimensions"]["total_points"],
                "max_abs_weight_percent_change": run["max_abs_weight_percent_change"],
                "main_driver_axis": main_driver_axis(run),
                "conclusion": conclusion_text(run, classification),
            }
        )
    return {
        "schema_version": 1,
        "stability_report": str(path),
        "classification_rules": {
            "stabil": "rank_flip_count == 0 and baseline_best_fraction == 1.0",
            "annahmeabhaengig": (
                "rank flips occur, but baseline_best_fraction >= 0.5"
            ),
            "instabil": "baseline_best_fraction < 0.5",
        },
        "conclusions": conclusions,
    }


def markdown_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Robustness Conclusions",
        "",
        f"Stability report: `{payload['stability_report']}`",
        "",
        "## Classification Rules",
        "",
        "| Label | Rule |",
        "|---|---|",
    ]
    for label, rule in payload["classification_rules"].items():
        lines.append(f"| `{label}` | {rule} |")

    lines.extend(
        [
            "",
            "## Conclusions",
            "",
            "| Run | Classification | Baseline best in matrix | Rank flips | Main driver |",
            "|---|---|---:|---:|---|",
        ]
    )
    for row in payload["conclusions"]:
        lines.append(
            f"| `{row['run_name']}` | `{row['classification']}` | "
            f"{row['baseline_best_fraction'] * 100.0:.1f}% | "
            f"{row['rank_flip_count']}/{row['total_matrix_points']} | "
            f"`{row['main_driver_axis']}` |"
        )

    lines.extend(["", "## Interpretation", ""])
    for row in payload["conclusions"]:
        lines.extend(
            [
                f"### `{row['run_name']}`",
                "",
                row["conclusion"],
                "",
            ]
        )

    return "\n".join(lines).rstrip() + "\n"


def write_report(output_dir: Path, payload: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "robustness_conclusions.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (output_dir / "robustness_conclusions.md").write_text(
        markdown_report(payload),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Classify robustness from a stability report"
    )
    parser.add_argument("stability_report_json")
    parser.add_argument(
        "--out_dir",
        default="reproduced/robustness_conclusions",
    )
    args = parser.parse_args()

    try:
        payload = classify_stability_report(Path(args.stability_report_json))
        write_report(Path(args.out_dir), payload)
    except Exception as exc:
        print(f"ROBUSTNESS CONCLUSIONS FAILED: {exc}")
        return 1

    print("=== Robustness Conclusions ===")
    print(f"Results: {Path(args.out_dir) / 'robustness_conclusions.json'}")
    print(f"Report: {Path(args.out_dir) / 'robustness_conclusions.md'}")
    for row in payload["conclusions"]:
        print(f"{row['run_name']}: {row['classification']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
