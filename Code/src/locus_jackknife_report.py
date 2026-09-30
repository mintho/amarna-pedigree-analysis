#!/usr/bin/env python3
"""Summarize locus-jackknife results into JSON and Markdown reports."""

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


def model_key(row: dict[str, Any]) -> str:
    return (
        f"{row['pedigree']} / "
        f"{row['founder_scenario']} / "
        f"{row['identity_hypothesis']}"
    )


def best_model_summary(row: dict[str, Any]) -> dict[str, Any]:
    best = row["results"][0]
    return {
        "model": model_key(best),
        "weight_percent": best["weight_percent"],
        "log_likelihood": best["log_likelihood"],
    }


def summarize_leave_one_out(row: dict[str, Any], baseline_best: str) -> dict[str, Any]:
    changes = row["changes"]
    rank_changes = [change for change in changes if change["rank_change"] != 0]
    max_abs_weight_change = max(
        (abs(change["weight_percent_change"]) for change in changes),
        default=0.0,
    )
    largest_weight_changes = [
        change
        for change in changes
        if abs(change["weight_percent_change"]) == max_abs_weight_change
    ]
    best = best_model_summary(row)
    return {
        "omitted_locus": row["omitted_locus"],
        "retained_loci": row["retained_loci"],
        "best_model": best["model"],
        "best_weight_percent": best["weight_percent"],
        "rank_flip": best["model"] != baseline_best,
        "rank_changes": rank_changes,
        "max_abs_weight_percent_change": max_abs_weight_change,
        "largest_weight_changes": largest_weight_changes,
    }


def summarize_jackknife(path: Path) -> dict[str, Any]:
    data = load_json(path)
    baseline_best_row = data["baseline"][0]
    baseline_best = model_key(baseline_best_row)
    leave_one_out = [
        summarize_leave_one_out(row, baseline_best)
        for row in data.get("leave_one_out", [])
    ]
    driving_loci = sorted(
        leave_one_out,
        key=lambda row: (
            row["rank_flip"],
            row["max_abs_weight_percent_change"],
        ),
        reverse=True,
    )
    return {
        "jackknife_path": str(path),
        "run_config": data["run_config"],
        "engine": data["engine"],
        "loci": data["loci"],
        "baseline_best_model": baseline_best,
        "baseline_best_weight_percent": baseline_best_row["weight_percent"],
        "rank_flip_count": sum(1 for row in leave_one_out if row["rank_flip"]),
        "max_abs_weight_percent_change": max(
            (row["max_abs_weight_percent_change"] for row in leave_one_out),
            default=0.0,
        ),
        "driving_loci": driving_loci,
        "leave_one_out": leave_one_out,
    }


def build_report(paths: list[Path]) -> dict[str, Any]:
    summaries = [summarize_jackknife(path) for path in paths]
    return {
        "schema_version": 1,
        "jackknife_files": [str(path) for path in paths],
        "runs": summaries,
    }


def markdown_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Locus Jackknife Report",
        "",
        "## Run Summary",
        "",
        "| Run | Baseline best | Baseline weight | Rank flips | Max weight change | Driving locus |",
        "|---|---|---:|---:|---:|---|",
    ]
    for run in payload["runs"]:
        driving = run["driving_loci"][0] if run["driving_loci"] else None
        driving_locus = driving["omitted_locus"] if driving else "n/a"
        lines.append(
            f"| `{Path(run['run_config']).name}` | "
            f"`{run['baseline_best_model']}` | "
            f"{run['baseline_best_weight_percent']:.3f}% | "
            f"{run['rank_flip_count']} | "
            f"{run['max_abs_weight_percent_change']:.3f} pp | "
            f"{driving_locus} |"
        )

    lines.extend(["", "## Leave-One-Out Details", ""])
    for run in payload["runs"]:
        lines.extend(
            [
                f"### `{Path(run['run_config']).name}`",
                "",
                "| Omitted locus | Best model after omission | Best weight | Rank flip | Max weight change | Largest affected model |",
                "|---|---|---:|---|---:|---|",
            ]
        )
        for row in run["leave_one_out"]:
            affected = row["largest_weight_changes"][0] if row["largest_weight_changes"] else None
            affected_model = affected["model"] if affected else "n/a"
            lines.append(
                f"| {row['omitted_locus']} | "
                f"`{row['best_model']}` | "
                f"{row['best_weight_percent']:.3f}% | "
                f"{'yes' if row['rank_flip'] else 'no'} | "
                f"{row['max_abs_weight_percent_change']:.3f} pp | "
                f"`{affected_model}` |"
            )
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_report(output_dir: Path, payload: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "jackknife_report.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (output_dir / "jackknife_report.md").write_text(
        markdown_report(payload),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize one or more locus-jackknife JSON files"
    )
    parser.add_argument("jackknife_json", nargs="+")
    parser.add_argument(
        "--out_dir",
        default="reproduced/locus_jackknife_report",
    )
    args = parser.parse_args()

    try:
        payload = build_report([Path(path) for path in args.jackknife_json])
        write_report(Path(args.out_dir), payload)
    except Exception as exc:
        print(f"JACKKNIFE REPORT FAILED: {exc}")
        return 1

    print("=== Locus Jackknife Report ===")
    print(f"Results: {Path(args.out_dir) / 'jackknife_report.json'}")
    print(f"Report: {Path(args.out_dir) / 'jackknife_report.md'}")
    for run in payload["runs"]:
        print(
            f"{Path(run['run_config']).name}: "
            f"{run['rank_flip_count']} rank flips, "
            f"max weight change {run['max_abs_weight_percent_change']:.3f} pp"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
