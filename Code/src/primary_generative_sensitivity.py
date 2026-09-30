#!/usr/bin/env python3
"""Full primary-comparison sensitivity analyses for the generative STR observation model."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from .locus_jackknife import collect_run_loci, evaluate_models_for_loci, load_json
    from .validate import validate_run
except ImportError:
    from locus_jackknife import collect_run_loci, evaluate_models_for_loci, load_json  # type: ignore[no-redef]
    from validate import validate_run  # type: ignore[no-redef]


SCENARIO_RUNS = {
    "a": Path("Data/models/runs/scenario_a_primary_models_all_loci_generative_str.json"),
    "b": Path("Data/models/runs/scenario_b_primary_models_all_loci_generative_str.json"),
    "c": Path("Data/models/runs/scenario_c_primary_models_all_loci_generative_str.json"),
    "d": Path("Data/models/runs/scenario_d_primary_models_all_loci_generative_str.json"),
}
BACKGROUND_ERRORS = (1e-4, 1e-3, 1e-2)
DROPOUTS = (0.1, 0.2, 0.3)
BASELINE_BACKGROUND_ERROR = 1e-3
BASELINE_DROPOUT = 0.2


def model_label(row: dict[str, Any]) -> str:
    return row["pedigree"].replace("_all_loci", "")


def expected_tied_labels(scenario: str) -> set[str]:
    return {
        f"hawass_v2_scenario_{scenario}",
        f"belmonte_v1_scenario_{scenario}",
    }


def run_with_observation_parameters(
    *,
    run_path: Path,
    run: dict[str, Any],
    loci: list[str],
    background_error: float,
    dropout: float,
    engine: str,
) -> list[dict[str, Any]]:
    varied = copy.deepcopy(run)
    varied["observation_model"]["background_error"] = background_error
    varied["observation_model"]["dropout"] = dropout
    return evaluate_models_for_loci(
        run_path=run_path,
        run=varied,
        loci=loci,
        engine=engine,
    )


def tied_structure_status(scenario: str, results: list[dict[str, Any]]) -> dict[str, Any]:
    top_two = results[:2]
    top_labels = {model_label(row) for row in top_two}
    expected = expected_tied_labels(scenario)
    tie_delta = (
        float(top_two[0]["log_likelihood"]) - float(top_two[1]["log_likelihood"])
        if len(top_two) == 2
        else float("nan")
    )
    nearest = results[2] if len(results) > 2 else None
    return {
        "top_labels": sorted(top_labels),
        "expected_top_labels": sorted(expected),
        "shared_structure_top": top_labels == expected and abs(tie_delta) < 1e-9,
        "tie_delta_log_likelihood": tie_delta,
        "nearest_non_tied_model": model_label(nearest) if nearest else "",
        "nearest_non_tied_delta_log_likelihood": (
            float(nearest["delta_log_likelihood"]) if nearest else float("nan")
        ),
        "nearest_non_tied_bayes_factor": (
            math.exp(-float(nearest["delta_log_likelihood"])) if nearest else float("nan")
        ),
    }


def compact_ranking(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "rank": row["rank"],
            "model": model_label(row),
            "log_likelihood": row["log_likelihood"],
            "delta_log_likelihood": row["delta_log_likelihood"],
            "weight_percent": row["weight_percent"],
        }
        for row in results
    ]


def run_sensitivity(engine: str = "factors") -> dict[str, Any]:
    scenario_payloads = []
    for scenario, relative_run_path in SCENARIO_RUNS.items():
        run_path = relative_run_path.resolve()
        validate_run(run_path)
        run = load_json(run_path)
        loci = collect_run_loci(run_path, run)

        parameter_grid = []
        for background_error in BACKGROUND_ERRORS:
            for dropout in DROPOUTS:
                print(
                    f"scenario {scenario}: background_error={background_error:g}, "
                    f"dropout={dropout:g}",
                    flush=True,
                )
                results = run_with_observation_parameters(
                    run_path=run_path,
                    run=run,
                    loci=loci,
                    background_error=background_error,
                    dropout=dropout,
                    engine=engine,
                )
                parameter_grid.append(
                    {
                        "background_error": background_error,
                        "dropout": dropout,
                        **tied_structure_status(scenario, results),
                        "ranking": compact_ranking(results),
                    }
                )

        leave_one_locus = []
        for omitted_locus in loci:
            retained_loci = [locus for locus in loci if locus != omitted_locus]
            print(f"scenario {scenario}: omit {omitted_locus}", flush=True)
            results = run_with_observation_parameters(
                run_path=run_path,
                run=run,
                loci=retained_loci,
                background_error=BASELINE_BACKGROUND_ERROR,
                dropout=BASELINE_DROPOUT,
                engine=engine,
            )
            leave_one_locus.append(
                {
                    "omitted_locus": omitted_locus,
                    "retained_loci": retained_loci,
                    **tied_structure_status(scenario, results),
                    "ranking": compact_ranking(results),
                }
            )

        scenario_payloads.append(
            {
                "scenario": scenario,
                "run_config": str(run_path),
                "loci": loci,
                "parameter_grid": parameter_grid,
                "leave_one_locus": leave_one_locus,
            }
        )

    return {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "engine": engine,
        "baseline": {
            "background_error": BASELINE_BACKGROUND_ERROR,
            "dropout": BASELINE_DROPOUT,
        },
        "background_error_values": list(BACKGROUND_ERRORS),
        "dropout_values": list(DROPOUTS),
        "scenarios": scenario_payloads,
    }


def parameter_rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for scenario in payload["scenarios"]:
        for point in scenario["parameter_grid"]:
            rows.append(
                {
                    "scenario": scenario["scenario"],
                    "background_error": point["background_error"],
                    "dropout": point["dropout"],
                    "shared_structure_top": point["shared_structure_top"],
                    "tie_delta_log_likelihood": point["tie_delta_log_likelihood"],
                    "nearest_non_tied_model": point["nearest_non_tied_model"],
                    "nearest_non_tied_delta_log_likelihood": point[
                        "nearest_non_tied_delta_log_likelihood"
                    ],
                    "nearest_non_tied_bayes_factor": point[
                        "nearest_non_tied_bayes_factor"
                    ],
                }
            )
    return rows


def leave_one_rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for scenario in payload["scenarios"]:
        for point in scenario["leave_one_locus"]:
            rows.append(
                {
                    "scenario": scenario["scenario"],
                    "omitted_locus": point["omitted_locus"],
                    "shared_structure_top": point["shared_structure_top"],
                    "tie_delta_log_likelihood": point["tie_delta_log_likelihood"],
                    "nearest_non_tied_model": point["nearest_non_tied_model"],
                    "nearest_non_tied_delta_log_likelihood": point[
                        "nearest_non_tied_delta_log_likelihood"
                    ],
                    "nearest_non_tied_bayes_factor": point[
                        "nearest_non_tied_bayes_factor"
                    ],
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def all_points_stable(payload: dict[str, Any]) -> bool:
    return all(row["shared_structure_top"] for row in parameter_rows(payload)) and all(
        row["shared_structure_top"] for row in leave_one_rows(payload)
    )


def markdown_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Full primary-comparison Sensitivity: Generative STR Observation Model",
        "",
        f"Engine: `{payload['engine']}`",
        f"Baseline background error: `{payload['baseline']['background_error']}`",
        f"Baseline dropout: `{payload['baseline']['dropout']}`",
        f"Overall stable tied structure: **{str(all_points_stable(payload)).upper()}**",
        "",
        "## Observation-Parameter Grid",
        "",
        "| Scenario | background_error | dropout | Tied structure top? | Nearest non-tied | ΔlogL | BF tied:nearest |",
        "|---|---:|---:|---|---|---:|---:|",
    ]
    for row in parameter_rows(payload):
        lines.append(
            f"| {row['scenario']} | {row['background_error']:.6g} | "
            f"{row['dropout']:.6g} | {row['shared_structure_top']} | "
            f"{row['nearest_non_tied_model']} | "
            f"{row['nearest_non_tied_delta_log_likelihood']:.6f} | "
            f"{row['nearest_non_tied_bayes_factor']:.6g} |"
        )

    lines.extend(
        [
            "",
            "## Leave-One-Locus Full-Ranking Reruns",
            "",
            "| Scenario | Omitted locus | Tied structure top? | Nearest non-tied | ΔlogL | BF tied:nearest |",
            "|---|---|---|---|---:|---:|",
        ]
    )
    for row in leave_one_rows(payload):
        lines.append(
            f"| {row['scenario']} | {row['omitted_locus']} | "
            f"{row['shared_structure_top']} | {row['nearest_non_tied_model']} | "
            f"{row['nearest_non_tied_delta_log_likelihood']:.6f} | "
            f"{row['nearest_non_tied_bayes_factor']:.6g} |"
        )
    return "\n".join(lines).rstrip() + "\n"


def write_outputs(payload: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "sensitivity.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (output_dir / "report.md").write_text(markdown_report(payload), encoding="utf-8")
    write_csv(
        output_dir / "parameter_grid.csv",
        parameter_rows(payload),
        [
            "scenario",
            "background_error",
            "dropout",
            "shared_structure_top",
            "tie_delta_log_likelihood",
            "nearest_non_tied_model",
            "nearest_non_tied_delta_log_likelihood",
            "nearest_non_tied_bayes_factor",
        ],
    )
    write_csv(
        output_dir / "leave_one_locus.csv",
        leave_one_rows(payload),
        [
            "scenario",
            "omitted_locus",
            "shared_structure_top",
            "tie_delta_log_likelihood",
            "nearest_non_tied_model",
            "nearest_non_tied_delta_log_likelihood",
            "nearest_non_tied_bayes_factor",
        ],
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run full primary generative observation-model sensitivity analyses"
    )
    parser.add_argument("--engine", choices=["factors", "enumeration"], default="factors")
    parser.add_argument(
        "--out-dir",
        default="reproduced/primary_generative_sensitivity",
        help="Output directory for sensitivity results",
    )
    args = parser.parse_args()

    try:
        payload = run_sensitivity(engine=args.engine)
        write_outputs(payload, Path(args.out_dir))
    except Exception as exc:
        print(f"SENSITIVITY FAILED: {exc}")
        return 1

    print(f"Wrote sensitivity outputs to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
