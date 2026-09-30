#!/usr/bin/env python3
"""Summarize primary all-locus rankings after reruns."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


SCENARIOS = ("a", "b", "c", "d")


def load_results(result_dir: Path) -> list[dict[str, Any]]:
    path = result_dir / "results.json"
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload["results"]


def model_label(row: dict[str, Any]) -> str:
    return row["pedigree"].replace("_all_loci", "")


def bayes_factor(delta_log_likelihood: float) -> float:
    return math.exp(-delta_log_likelihood)


def locus_map(row: dict[str, Any]) -> dict[str, float]:
    return {
        locus["locus"]: float(locus["log_likelihood"])
        for locus in row.get("loci", [])
    }


def build_summary(input_dir: Path = Path("Results/Primary_All_Locus")) -> dict[str, Any]:
    scenarios = []
    for scenario in SCENARIOS:
        generative_results = load_results(
            input_dir / f"scenario_{scenario}"
        )
        best = generative_results[0]
        second = generative_results[1]
        nearest = generative_results[2]
        best_loci = locus_map(best)
        nearest_loci = locus_map(nearest)
        per_locus_vs_nearest = [
            {
                "locus": locus,
                "best_log_likelihood": best_loci[locus],
                "nearest_log_likelihood": nearest_loci[locus],
                "delta_log_likelihood": best_loci[locus] - nearest_loci[locus],
            }
            for locus in sorted(best_loci)
        ]
        scenarios.append(
            {
                "scenario": scenario,
                "best_model": model_label(best),
                "second_model": model_label(second),
                "tie_delta_log_likelihood": (
                    float(best["log_likelihood"]) - float(second["log_likelihood"])
                ),
                "nearest_non_tied_model": model_label(nearest),
                "nearest_non_tied_delta_log_likelihood": float(
                    nearest["delta_log_likelihood"]
                ),
                "nearest_non_tied_bayes_factor": bayes_factor(
                    float(nearest["delta_log_likelihood"])
                ),
                "generative_results": generative_results,
                "per_locus_vs_nearest": per_locus_vs_nearest,
            }
        )
    return {"schema_version": 1, "scenarios": scenarios}


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def ranking_rows(summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for scenario in summary["scenarios"]:
        for row in scenario["generative_results"]:
            rows.append(
                {
                    "scenario": scenario["scenario"],
                    "rank": row["rank"],
                    "model": model_label(row),
                    "log_likelihood": row["log_likelihood"],
                    "delta_log_likelihood": row["delta_log_likelihood"],
                    "weight_percent": row["weight_percent"],
                    "observation_model_type": row["observation_model"]["type"],
                }
            )
    return rows


def structural_label(label: str) -> str:
    if label.startswith("hawass_v2_") or label.startswith("belmonte_v1_"):
        return "shared_kv55_kv35yl_tutankhamun_structure"
    return label


def structural_ranking_rows(summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for scenario in summary["scenarios"]:
        collapsed: dict[str, dict[str, Any]] = {}
        for row in scenario["generative_results"]:
            label = model_label(row)
            structure = structural_label(label)
            if structure not in collapsed:
                collapsed[structure] = {
                    "scenario": scenario["scenario"],
                    "structure": structure,
                    "represented_by": label,
                    "log_likelihood": float(row["log_likelihood"]),
                    "observation_model_type": row["observation_model"]["type"],
                }
            elif label not in collapsed[structure]["represented_by"].split(" / "):
                collapsed[structure]["represented_by"] += f" / {label}"

        ranked = sorted(collapsed.values(), key=lambda item: item["log_likelihood"], reverse=True)
        best = ranked[0]["log_likelihood"]
        denom = sum(math.exp(item["log_likelihood"] - best) for item in ranked)
        for rank, item in enumerate(ranked, start=1):
            delta = item["log_likelihood"] - best
            rows.append(
                {
                    **item,
                    "rank": rank,
                    "delta_log_likelihood": delta,
                    "weight_percent": 100.0 * math.exp(delta) / denom,
                }
            )
    return rows


def per_locus_rows(summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for scenario in summary["scenarios"]:
        for row in scenario["per_locus_vs_nearest"]:
            rows.append(
                {
                    "scenario": scenario["scenario"],
                    "nearest_non_tied_model": scenario["nearest_non_tied_model"],
                    **row,
                }
            )
    return rows


def markdown_report(summary: dict[str, Any]) -> str:
    lines = [
        "# Primary All-Locus Results: Generative STR Observation Model",
        "",
        "These tables summarize the primary all-locus calculations under the "
        "normalized `generative_str` observation model.",
        "",
        "Hawass-derived and Belmonte are also reported as one unique scored "
        "structure for likelihood-weight normalization, because they assign the "
        "same sampled profiles to the same inheritance positions.",
        "",
        "## Scenario Summary",
        "",
        "| Scenario | Tied best readings | Tie ΔlogL | Nearest non-tied | ΔlogL | Support ratio |",
        "|---|---|---:|---|---:|---:|",
    ]
    for scenario in summary["scenarios"]:
        tied = f"{scenario['best_model']} / {scenario['second_model']}"
        ratio_label = "BF" if scenario["scenario"] in {"c", "d"} else "LR"
        lines.append(
            f"| {scenario['scenario']} | {tied} | "
            f"{scenario['tie_delta_log_likelihood']:.6g} | "
            f"{scenario['nearest_non_tied_model']} | "
            f"{scenario['nearest_non_tied_delta_log_likelihood']:.6f} | "
            f"{ratio_label} {scenario['nearest_non_tied_bayes_factor']:.6g} |"
        )

    lines.extend(
        [
            "",
            "## Unique Scored-Structure Rankings",
            "",
            "| Scenario | Rank | Structure | Represented by | logL | ΔlogL | Weight |",
            "|---|---:|---|---|---:|---:|---:|",
        ]
    )
    for row in structural_ranking_rows(summary):
        lines.append(
            f"| {row['scenario']} | {row['rank']} | {row['structure']} | "
            f"{row['represented_by']} | {row['log_likelihood']:.6f} | "
            f"{row['delta_log_likelihood']:.6f} | {row['weight_percent']:.3f}% |"
        )

    lines.extend(
        [
            "",
            "## Labelled Family Outputs Before Structural Collapse",
            "",
            "These rows preserve the machine-readable labelled family outputs. Hawass-derived "
            "and Belmonte are not independent scored structures; the structural weights above "
            "are the weights used for manuscript interpretation.",
            "",
            "| Scenario | Rank | Model | logL | ΔlogL | Labelled weight |",
            "|---|---:|---|---:|---:|---:|",
        ]
    )
    for row in ranking_rows(summary):
        lines.append(
            f"| {row['scenario']} | {row['rank']} | {row['model']} | "
            f"{row['log_likelihood']:.6f} | {row['delta_log_likelihood']:.6f} | "
            f"{row['weight_percent']:.3f}% |"
        )

    lines.extend(
        [
            "",
            "## Per-Locus Contribution Against Nearest Non-Tied Model",
            "",
            "Positive values favour the tied best structure over the nearest non-tied alternative.",
            "",
            "| Scenario | Nearest non-tied model | Locus | ΔlogL |",
            "|---|---|---|---:|",
        ]
    )
    for row in per_locus_rows(summary):
        lines.append(
            f"| {row['scenario']} | {row['nearest_non_tied_model']} | "
            f"{row['locus']} | {row['delta_log_likelihood']:.6f} |"
        )

    return "\n".join(lines).rstrip() + "\n"


def write_outputs(summary: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (output_dir / "report.md").write_text(markdown_report(summary), encoding="utf-8")
    write_csv(
        output_dir / "ranking_summary.csv",
        ranking_rows(summary),
        [
            "scenario",
            "rank",
            "model",
            "log_likelihood",
            "delta_log_likelihood",
            "weight_percent",
            "observation_model_type",
        ],
    )
    write_csv(
        output_dir / "structural_ranking_summary.csv",
        structural_ranking_rows(summary),
        [
            "scenario",
            "rank",
            "structure",
            "represented_by",
            "log_likelihood",
            "delta_log_likelihood",
            "weight_percent",
            "observation_model_type",
        ],
    )
    write_csv(
        output_dir / "per_locus_vs_nearest.csv",
        per_locus_rows(summary),
        [
            "scenario",
            "nearest_non_tied_model",
            "locus",
            "best_log_likelihood",
            "nearest_log_likelihood",
            "delta_log_likelihood",
        ],
    )

def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize primary all-locus outputs")
    parser.add_argument("--input-dir", type=Path, default=Path("Results/Primary_All_Locus"))
    parser.add_argument(
        "--out-dir",
        default="reproduced/Primary_All_Locus/summary",
        help="Output directory for summary tables",
    )
    args = parser.parse_args()

    try:
        summary = build_summary(args.input_dir)
        write_outputs(summary, Path(args.out_dir))
    except Exception as exc:
        print(f"SUMMARY FAILED: {exc}")
        return 1

    print(f"Wrote primary all-locus summary to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
