#!/usr/bin/env python3
"""Sensitivity matrix over founder priors, background_error, and dropout."""

from __future__ import annotations

import argparse
import copy
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from .compare import add_model_comparison
    from .background_error_sensitivity import collect_common_loci, parse_background_error_values
    from .founder_prior_sensitivity import (
        base_model_specs_for_founder_variation,
        candidate_founder_refs,
        founder_scenario_label,
        model_changes_vs_baseline,
        normalize_candidate_ref,
        structural_model_key,
    )
    from .likelihood_enumeration import load_json
    from .locus_jackknife import evaluate_model_for_loci
    from .validate import validate_run
except ImportError:
    from compare import add_model_comparison  # type: ignore[no-redef]
    from background_error_sensitivity import (  # type: ignore[no-redef]
        collect_common_loci,
        parse_background_error_values,
    )
    from founder_prior_sensitivity import (  # type: ignore[no-redef]
        base_model_specs_for_founder_variation,
        candidate_founder_refs,
        founder_scenario_label,
        model_changes_vs_baseline,
        normalize_candidate_ref,
        structural_model_key,
    )
    from likelihood_enumeration import load_json  # type: ignore[no-redef]
    from locus_jackknife import evaluate_model_for_loci  # type: ignore[no-redef]
    from validate import validate_run  # type: ignore[no-redef]


DEFAULT_DROPOUT_VALUES = [0.0, 0.05, 0.1, 0.2, 0.35]


def validate_rate(value: float, name: str) -> float:
    rate = float(value)
    if rate < 0.0 or rate >= 1.0:
        raise ValueError(f"{name} values must satisfy 0 <= {name} < 1")
    return rate


def unique_sorted_rates(values: list[float], name: str) -> list[float]:
    clean = [validate_rate(value, name) for value in values]
    return sorted(set(clean))


def parse_dropout_values(raw: str | None) -> list[float]:
    if raw is None or not raw.strip():
        return list(DEFAULT_DROPOUT_VALUES)
    values = []
    for item in raw.split(","):
        stripped = item.strip()
        if not stripped:
            continue
        values.append(float(stripped))
    if not values:
        raise ValueError("at least one dropout value is required")
    return values


def run_matrix_point(
    *,
    run_path: Path,
    run: dict[str, Any],
    loci: list[str],
    founder_ref: str,
    background_error: float,
    dropout: float,
    engine: str,
) -> list[dict[str, Any]]:
    varied_run = copy.deepcopy(run)
    varied_run.setdefault("observation_model", {})["background_error"] = validate_rate(
        background_error,
        "background_error",
    )
    varied_run.setdefault("observation_model", {})["dropout"] = validate_rate(
        dropout,
        "dropout",
    )

    specs = []
    for spec in base_model_specs_for_founder_variation(varied_run):
        varied_spec = dict(spec)
        varied_spec["founder_scenario"] = founder_ref
        specs.append(varied_spec)

    results = [
        evaluate_model_for_loci(
            run_path=run_path,
            run=varied_run,
            spec=spec,
            loci=loci,
            engine=engine,
        )
        for spec in specs
    ]
    return add_model_comparison(results)


def summarize_matrix(
    *,
    baseline: list[dict[str, Any]],
    matrix_points: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    baseline_by_model = {structural_model_key(row): row for row in baseline}
    summary = []
    for key, base in baseline_by_model.items():
        ranks = []
        weights = []
        log_likelihoods = []
        best_count = 0
        best_points = []
        for point in matrix_points:
            rows = {structural_model_key(row): row for row in point["results"]}
            row = rows[key]
            ranks.append(row["rank"])
            weights.append(row["weight_percent"])
            log_likelihoods.append(row["log_likelihood"])
            if point["results"] and structural_model_key(point["results"][0]) == key:
                best_count += 1
                best_points.append(point["matrix_key"])
        summary.append(
            {
                "model": key,
                "baseline_rank": base["rank"],
                "baseline_weight_percent": base["weight_percent"],
                "baseline_log_likelihood": base["log_likelihood"],
                "best_count": best_count,
                "best_fraction": best_count / (len(matrix_points) or 1),
                "best_points": best_points,
                "rank_min": min(ranks),
                "rank_max": max(ranks),
                "weight_percent_min": min(weights),
                "weight_percent_max": max(weights),
                "max_abs_weight_percent_change": max(
                    abs(weight - base["weight_percent"]) for weight in weights
                ),
                "log_likelihood_min": min(log_likelihoods),
                "log_likelihood_max": max(log_likelihoods),
            }
        )
    return sorted(summary, key=lambda row: row["baseline_rank"])


def run_sensitivity_matrix(
    *,
    run_config: Path,
    founder_scenarios: list[str] | None = None,
    baseline_founder_scenario: str | None = None,
    background_error_values: list[float] | None = None,
    dropout_values: list[float] | None = None,
    engine: str = "factors",
) -> dict[str, Any]:
    run_path = run_config.resolve()
    validate_run(run_path)
    run = load_json(run_path)
    loci = collect_common_loci(run_path, run)

    refs = candidate_founder_refs(
        run_path=run_path,
        run=run,
        founder_scenarios=founder_scenarios,
    )
    baseline_ref = (
        normalize_candidate_ref(run_path, baseline_founder_scenario)
        if baseline_founder_scenario
        else refs[0]
    )
    if baseline_ref not in refs:
        refs = [baseline_ref] + refs

    baseline_background_error = validate_rate(float(run["observation_model"]["background_error"]), "background_error")
    baseline_dropout = validate_rate(float(run["observation_model"]["dropout"]), "dropout")
    background_errors = unique_sorted_rates(
        list(background_error_values or parse_background_error_values(None)) + [baseline_background_error],
        "background_error",
    )
    dropouts = unique_sorted_rates(
        list(dropout_values or DEFAULT_DROPOUT_VALUES) + [baseline_dropout],
        "dropout",
    )

    baseline = run_matrix_point(
        run_path=run_path,
        run=run,
        loci=loci,
        founder_ref=baseline_ref,
        background_error=baseline_background_error,
        dropout=baseline_dropout,
        engine=engine,
    )

    matrix_points = []
    for founder_ref in refs:
        label = founder_scenario_label(run_path, founder_ref)
        for background_error in background_errors:
            for dropout in dropouts:
                results = run_matrix_point(
                    run_path=run_path,
                    run=run,
                    loci=loci,
                    founder_ref=founder_ref,
                    background_error=background_error,
                    dropout=dropout,
                    engine=engine,
                )
                matrix_key = (
                    f"{label['id']} | background_error={background_error:.6g} | "
                    f"dropout={dropout:.6g}"
                )
                matrix_points.append(
                    {
                        "matrix_key": matrix_key,
                        "founder_scenario_ref": founder_ref,
                        "founder_scenario_path": label["path"],
                        "founder_scenario_id": label["id"],
                        "background_error": background_error,
                        "dropout": dropout,
                        "results": results,
                        "changes_vs_baseline": model_changes_vs_baseline(
                            baseline,
                            results,
                        ),
                    }
                )

    return {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_config": str(run_path),
        "engine": engine,
        "baseline_founder_scenario_ref": baseline_ref,
        "baseline_founder_scenario_id": founder_scenario_label(run_path, baseline_ref)["id"],
        "baseline_background_error": baseline_background_error,
        "baseline_dropout": baseline_dropout,
        "founder_scenarios": [founder_scenario_label(run_path, ref) for ref in refs],
        "background_error_values": background_errors,
        "dropout_values": dropouts,
        "loci": loci,
        "baseline": baseline,
        "matrix_points": matrix_points,
        "summary": summarize_matrix(
            baseline=baseline,
            matrix_points=matrix_points,
        ),
    }


def markdown_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Sensitivity Matrix Report",
        "",
        f"Run configuration: `{payload['run_config']}`",
        f"Engine: `{payload['engine']}`",
        f"Baseline founder scenario: `{payload['baseline_founder_scenario_id']}`",
        f"Baseline background error: `{payload['baseline_background_error']}`",
        f"Baseline dropout: `{payload['baseline_dropout']}`",
        "",
        "## Summary",
        "",
        "| Model | Baseline rank | Best fraction | Rank range | Weight range | Max weight change |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in payload["summary"]:
        lines.append(
            f"| {row['model']} | {row['baseline_rank']} | "
            f"{row['best_fraction'] * 100.0:.1f}% | "
            f"{row['rank_min']}..{row['rank_max']} | "
            f"{row['weight_percent_min']:.3f}%..{row['weight_percent_max']:.3f}% | "
            f"{row['max_abs_weight_percent_change']:.3f} pp |"
        )

    lines.extend(
        [
            "",
            "## Matrix Overview",
            "",
            "| Founder scenario | Background error | Dropout | Rank 1 model | Rank 1 weight |",
            "|---|---:|---:|---|---:|",
        ]
    )
    for point in payload["matrix_points"]:
        best = point["results"][0]
        lines.append(
            f"| {point['founder_scenario_id']} | {point['background_error']:.6g} | "
            f"{point['dropout']:.6g} | {structural_model_key(best)} | "
            f"{best['weight_percent']:.3f}% |"
        )

    lines.extend(["", "## Model Weights", ""])
    for point in payload["matrix_points"]:
        lines.append(f"### {point['matrix_key']}")
        lines.append("")
        lines.append("| Rank | Model | logL | Weight | Delta weight vs baseline |")
        lines.append("|---:|---|---:|---:|---:|")
        changes = {
            row["model"]: row["weight_percent_change"]
            for row in point["changes_vs_baseline"]
        }
        for row in point["results"]:
            key = structural_model_key(row)
            lines.append(
                f"| {row['rank']} | {key} | {row['log_likelihood']:.6f} | "
                f"{row['weight_percent']:.3f}% | {changes[key]:.3f} pp |"
            )
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_matrix_json(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def write_matrix_report(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(markdown_report(payload), encoding="utf-8")


def default_output_dir(run_config: Path) -> Path:
    return Path("reproduced") / f"{run_config.stem}_sensitivity_matrix"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate pedigree likelihood sensitivity over founder priors, "
            "background_error values, and dropout assumptions"
        )
    )
    parser.add_argument("run_config", help="Path to a run configuration JSON file")
    parser.add_argument(
        "--founder_scenarios",
        nargs="+",
        default=None,
        help=(
            "Alternative founder scenario JSON files. If omitted, the founder "
            "scenarios already referenced by the run are used."
        ),
    )
    parser.add_argument(
        "--baseline_founder_scenario",
        default=None,
        help="Founder scenario used as baseline. Defaults to the first candidate.",
    )
    parser.add_argument(
        "--background_errors",
        default=None,
        help="Comma-separated background_error grid.",
    )
    parser.add_argument(
        "--dropouts",
        default=None,
        help=(
            "Comma-separated dropout grid. Defaults to "
            f"{','.join(str(value) for value in DEFAULT_DROPOUT_VALUES)}"
        ),
    )
    parser.add_argument(
        "--engine",
        choices=["factors", "enumeration"],
        default="factors",
    )
    parser.add_argument(
        "--out_dir",
        default=None,
        help="Output directory. Defaults to reproduced/<run_stem>_sensitivity_matrix",
    )
    args = parser.parse_args()

    run_config = Path(args.run_config)
    out_dir = Path(args.out_dir) if args.out_dir else default_output_dir(run_config)
    try:
        payload = run_sensitivity_matrix(
            run_config=run_config,
            founder_scenarios=args.founder_scenarios,
            baseline_founder_scenario=args.baseline_founder_scenario,
            background_error_values=parse_background_error_values(args.background_errors),
            dropout_values=parse_dropout_values(args.dropouts),
            engine=args.engine,
        )
        write_matrix_json(out_dir / "matrix.json", payload)
        write_matrix_report(out_dir / "matrix.md", payload)
    except Exception as exc:
        print(f"SENSITIVITY MATRIX FAILED: {exc}")
        return 1

    print("=== Sensitivity Matrix ===")
    print(f"Run: {payload['run_config']}")
    print(f"Engine: {payload['engine']}")
    print(f"Founder scenarios: {len(payload['founder_scenarios'])}")
    print(f"Background error values: {payload['background_error_values']}")
    print(f"Dropout values: {payload['dropout_values']}")
    print(f"Results: {out_dir / 'matrix.json'}")
    print(f"Report: {out_dir / 'matrix.md'}")
    for row in payload["summary"]:
        print(
            f"{row['model']}: best {row['best_fraction'] * 100.0:.1f}%, "
            f"rank {row['rank_min']}..{row['rank_max']}, "
            f"weight {row['weight_percent_min']:.3f}%.."
            f"{row['weight_percent_max']:.3f}%"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
