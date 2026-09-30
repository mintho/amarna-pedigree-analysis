#!/usr/bin/env python3
"""Sensitivity analysis over observation background-error rates."""

from __future__ import annotations

import argparse
import copy
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from .locus_jackknife import (
        evaluate_models_for_loci,
        model_changes,
        model_key,
    )
    from .likelihood_enumeration import load_json, model_specs_from_run, resolve_ref
    from .validate import validate_run
except ImportError:
    from locus_jackknife import (  # type: ignore[no-redef]
        evaluate_models_for_loci,
        model_changes,
        model_key,
    )
    from likelihood_enumeration import (  # type: ignore[no-redef]
        load_json,
        model_specs_from_run,
        resolve_ref,
    )
    from validate import validate_run  # type: ignore[no-redef]


DEFAULT_BACKGROUND_ERROR_VALUES = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2]


def validate_background_error(value: float) -> float:
    background_error = float(value)
    if background_error < 0.0 or background_error >= 1.0:
        raise ValueError("background_error values must satisfy 0 <= background_error < 1")
    return background_error


def unique_sorted_background_errors(values: list[float]) -> list[float]:
    clean = [validate_background_error(value) for value in values]
    return sorted(set(clean))


def parse_background_error_values(raw: str | None) -> list[float]:
    if raw is None or not raw.strip():
        return list(DEFAULT_BACKGROUND_ERROR_VALUES)
    values = []
    for item in raw.split(","):
        stripped = item.strip()
        if not stripped:
            continue
        values.append(float(stripped))
    if not values:
        raise ValueError("at least one background_error value is required")
    return values


def collect_common_loci(run_path: Path, run: dict[str, Any]) -> list[str]:
    loci_by_model: list[tuple[str, ...]] = []
    for spec in model_specs_from_run(run):
        observations_ref = spec.get("observations", run.get("observations"))
        observations = load_json(resolve_ref(run_path, observations_ref))
        loci = tuple(observations.get("loci", []))
        if not loci:
            raise ValueError(f"{observations_ref}: no loci available")
        loci_by_model.append(loci)

    first = loci_by_model[0]
    for loci in loci_by_model[1:]:
        if loci != first:
            raise ValueError(
                "error-rate sensitivity requires all model specs to use the "
                "same locus list"
            )
    return list(first)


def run_with_background_error(
    *,
    run_path: Path,
    run: dict[str, Any],
    loci: list[str],
    background_error: float,
    engine: str,
) -> list[dict[str, Any]]:
    varied_run = copy.deepcopy(run)
    varied_run.setdefault("observation_model", {})["background_error"] = validate_background_error(background_error)
    return evaluate_models_for_loci(
        run_path=run_path,
        run=varied_run,
        loci=loci,
        engine=engine,
    )


def summarize_sensitivity(
    *,
    baseline: list[dict[str, Any]],
    sensitivity_points: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    baseline_by_model = {model_key(row): row for row in baseline}
    summary = []
    for key, base in baseline_by_model.items():
        ranks = []
        weights = []
        log_likelihoods = []
        best_count = 0
        for point in sensitivity_points:
            rows = {model_key(row): row for row in point["results"]}
            row = rows[key]
            ranks.append(row["rank"])
            weights.append(row["weight_percent"])
            log_likelihoods.append(row["log_likelihood"])
            if point["results"] and model_key(point["results"][0]) == key:
                best_count += 1
        summary.append(
            {
                "model": key,
                "baseline_rank": base["rank"],
                "baseline_weight_percent": base["weight_percent"],
                "baseline_log_likelihood": base["log_likelihood"],
                "best_count": best_count,
                "best_fraction": best_count / (len(sensitivity_points) or 1),
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


def run_background_error_sensitivity(
    *,
    run_config: Path,
    background_error_values: list[float] | None = None,
    engine: str = "factors",
) -> dict[str, Any]:
    run_path = run_config.resolve()
    validate_run(run_path)
    run = load_json(run_path)
    baseline_background_error = validate_background_error(float(run["observation_model"]["background_error"]))
    background_errors = unique_sorted_background_errors(
        list(background_error_values or DEFAULT_BACKGROUND_ERROR_VALUES) + [baseline_background_error]
    )
    loci = collect_common_loci(run_path, run)
    baseline = run_with_background_error(
        run_path=run_path,
        run=run,
        loci=loci,
        background_error=baseline_background_error,
        engine=engine,
    )

    sensitivity_points = []
    for background_error in background_errors:
        results = run_with_background_error(
            run_path=run_path,
            run=run,
            loci=loci,
            background_error=background_error,
            engine=engine,
        )
        sensitivity_points.append(
            {
                "background_error": background_error,
                "results": results,
                "changes_vs_baseline": model_changes(baseline, results),
            }
        )

    return {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_config": str(run_path),
        "engine": engine,
        "baseline_background_error": baseline_background_error,
        "background_error_values": background_errors,
        "loci": loci,
        "baseline": baseline,
        "sensitivity_points": sensitivity_points,
        "summary": summarize_sensitivity(
            baseline=baseline,
            sensitivity_points=sensitivity_points,
        ),
    }


def markdown_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Error-Rate Sensitivity Report",
        "",
        f"Run configuration: `{payload['run_config']}`",
        f"Engine: `{payload['engine']}`",
        f"Baseline background error: `{payload['baseline_background_error']}`",
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
            "## Background error Grid",
            "",
            "| Background error | Rank 1 model | Rank 1 weight |",
            "|---:|---|---:|",
        ]
    )
    for point in payload["sensitivity_points"]:
        best = point["results"][0]
        lines.append(
            f"| {point['background_error']:.6g} | {model_key(best)} | "
            f"{best['weight_percent']:.3f}% |"
        )

    lines.extend(["", "## Model Weights", ""])
    for point in payload["sensitivity_points"]:
        lines.append(f"### Background error {point['background_error']:.6g}")
        lines.append("")
        lines.append("| Rank | Model | logL | Weight | Δweight vs baseline |")
        lines.append("|---:|---|---:|---:|---:|")
        changes = {
            row["model"]: row["weight_percent_change"]
            for row in point["changes_vs_baseline"]
        }
        for row in point["results"]:
            key = model_key(row)
            lines.append(
                f"| {row['rank']} | {key} | {row['log_likelihood']:.6f} | "
                f"{row['weight_percent']:.3f}% | {changes[key]:.3f} pp |"
            )
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_sensitivity_json(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def write_sensitivity_report(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(markdown_report(payload), encoding="utf-8")


def default_output_dir(run_config: Path) -> Path:
    return Path("reproduced") / f"{run_config.stem}_background_error_sensitivity"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate pedigree likelihood sensitivity over background_error values"
    )
    parser.add_argument("run_config", help="Path to a run configuration JSON file")
    parser.add_argument(
        "--background_errors",
        default=None,
        help=(
            "Comma-separated background_error grid. Defaults to "
            f"{','.join(str(value) for value in DEFAULT_BACKGROUND_ERROR_VALUES)}"
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
        help="Output directory. Defaults to reproduced/<run_stem>_background_error_sensitivity",
    )
    args = parser.parse_args()

    run_config = Path(args.run_config)
    out_dir = Path(args.out_dir) if args.out_dir else default_output_dir(run_config)
    try:
        payload = run_background_error_sensitivity(
            run_config=run_config,
            background_error_values=parse_background_error_values(args.background_errors),
            engine=args.engine,
        )
        write_sensitivity_json(out_dir / "sensitivity.json", payload)
        write_sensitivity_report(out_dir / "sensitivity.md", payload)
    except Exception as exc:
        print(f"ERROR-RATE SENSITIVITY FAILED: {exc}")
        return 1

    print("=== Error-Rate Sensitivity ===")
    print(f"Run: {payload['run_config']}")
    print(f"Engine: {payload['engine']}")
    print(f"Baseline background error: {payload['baseline_background_error']}")
    print(f"Results: {out_dir / 'sensitivity.json'}")
    print(f"Report: {out_dir / 'sensitivity.md'}")
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
