#!/usr/bin/env python3
"""Sensitivity analysis over alternative founder prior scenarios."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from .compare import add_model_comparison
    from .background_error_sensitivity import collect_common_loci
    from .likelihood_enumeration import load_json, model_specs_from_run, resolve_ref
    from .locus_jackknife import evaluate_model_for_loci
    from .validate import validate_run
except ImportError:
    from compare import add_model_comparison  # type: ignore[no-redef]
    from background_error_sensitivity import collect_common_loci  # type: ignore[no-redef]
    from likelihood_enumeration import (  # type: ignore[no-redef]
        load_json,
        model_specs_from_run,
        resolve_ref,
    )
    from locus_jackknife import evaluate_model_for_loci  # type: ignore[no-redef]
    from validate import validate_run  # type: ignore[no-redef]


def structural_model_key(row: dict[str, Any]) -> str:
    return f"{row['pedigree']} / {row['identity_hypothesis']}"


def normalize_candidate_ref(run_path: Path, raw_ref: str) -> str:
    candidate = Path(raw_ref)
    if candidate.is_absolute():
        return str(candidate)
    if candidate.exists():
        return str(candidate.resolve())
    return raw_ref


def default_founder_refs(run: dict[str, Any]) -> list[str]:
    if "models" in run:
        refs = []
        for spec in model_specs_from_run(run):
            ref = spec["founder_scenario"]
            if ref not in refs:
                refs.append(ref)
        return refs
    return list(run["founder_scenarios"])


def base_model_specs_for_founder_variation(run: dict[str, Any]) -> list[dict[str, str]]:
    """Return model specs with founder scenarios removed from the comparison axis."""
    if "models" in run:
        seen = set()
        specs = []
        for spec in model_specs_from_run(run):
            varied_spec = dict(spec)
            varied_spec.pop("founder_scenario", None)
            key = tuple(sorted(varied_spec.items()))
            if key not in seen:
                seen.add(key)
                specs.append(varied_spec)
        return specs

    return [
        {
            "pedigree": pedigree,
            "identity_hypothesis": identity_hypothesis,
        }
        for pedigree in run["pedigrees"]
        for identity_hypothesis in run["identity_hypotheses"]
    ]


def candidate_founder_refs(
    *,
    run_path: Path,
    run: dict[str, Any],
    founder_scenarios: list[str] | None,
) -> list[str]:
    refs = founder_scenarios or default_founder_refs(run)
    if not refs:
        raise ValueError("at least one founder scenario is required")
    normalized = [normalize_candidate_ref(run_path, ref) for ref in refs]
    unique = []
    for ref in normalized:
        if ref not in unique:
            unique.append(ref)
    return unique


def run_with_founder_scenario(
    *,
    run_path: Path,
    run: dict[str, Any],
    loci: list[str],
    founder_ref: str,
    engine: str,
) -> list[dict[str, Any]]:
    specs = []
    for spec in base_model_specs_for_founder_variation(run):
        varied_spec = dict(spec)
        varied_spec["founder_scenario"] = founder_ref
        specs.append(varied_spec)

    results = [
        evaluate_model_for_loci(
            run_path=run_path,
            run=run,
            spec=spec,
            loci=loci,
            engine=engine,
        )
        for spec in specs
    ]
    return add_model_comparison(results)


def founder_scenario_label(run_path: Path, founder_ref: str) -> dict[str, str]:
    path = resolve_ref(run_path, founder_ref)
    scenario = load_json(path)
    return {
        "ref": founder_ref,
        "path": str(path),
        "id": scenario["id"],
        "description": scenario.get("description", ""),
    }


def model_changes_vs_baseline(
    baseline: list[dict[str, Any]],
    results: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    baseline_by_key = {structural_model_key(row): row for row in baseline}
    changes = []
    for row in results:
        key = structural_model_key(row)
        base = baseline_by_key[key]
        changes.append(
            {
                "model": key,
                "baseline_founder_scenario": base["founder_scenario"],
                "candidate_founder_scenario": row["founder_scenario"],
                "baseline_rank": base["rank"],
                "candidate_rank": row["rank"],
                "rank_change": row["rank"] - base["rank"],
                "baseline_weight_percent": base["weight_percent"],
                "candidate_weight_percent": row["weight_percent"],
                "weight_percent_change": (
                    row["weight_percent"] - base["weight_percent"]
                ),
                "baseline_log_likelihood": base["log_likelihood"],
                "candidate_log_likelihood": row["log_likelihood"],
                "log_likelihood_change": (
                    row["log_likelihood"] - base["log_likelihood"]
                ),
            }
        )
    return sorted(
        changes,
        key=lambda row: (
            abs(row["rank_change"]),
            abs(row["weight_percent_change"]),
        ),
        reverse=True,
    )


def summarize_sensitivity(
    *,
    baseline: list[dict[str, Any]],
    sensitivity_points: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    baseline_by_model = {structural_model_key(row): row for row in baseline}
    summary = []
    for key, base in baseline_by_model.items():
        ranks = []
        weights = []
        log_likelihoods = []
        best_count = 0
        for point in sensitivity_points:
            rows = {structural_model_key(row): row for row in point["results"]}
            row = rows[key]
            ranks.append(row["rank"])
            weights.append(row["weight_percent"])
            log_likelihoods.append(row["log_likelihood"])
            if point["results"] and structural_model_key(point["results"][0]) == key:
                best_count += 1
        summary.append(
            {
                "model": key,
                "baseline_founder_scenario": base["founder_scenario"],
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


def run_founder_prior_sensitivity(
    *,
    run_config: Path,
    founder_scenarios: list[str] | None = None,
    baseline_founder_scenario: str | None = None,
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

    baseline = run_with_founder_scenario(
        run_path=run_path,
        run=run,
        loci=loci,
        founder_ref=baseline_ref,
        engine=engine,
    )

    sensitivity_points = []
    for ref in refs:
        results = run_with_founder_scenario(
            run_path=run_path,
            run=run,
            loci=loci,
            founder_ref=ref,
            engine=engine,
        )
        label = founder_scenario_label(run_path, ref)
        sensitivity_points.append(
            {
                "founder_scenario_ref": ref,
                "founder_scenario_path": label["path"],
                "founder_scenario_id": label["id"],
                "description": label["description"],
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
        "founder_scenarios": [
            founder_scenario_label(run_path, ref) for ref in refs
        ],
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
        "# Founder-Prior Sensitivity Report",
        "",
        f"Run configuration: `{payload['run_config']}`",
        f"Engine: `{payload['engine']}`",
        f"Baseline founder scenario: `{payload['baseline_founder_scenario_id']}`",
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
            "## Founder Scenarios",
            "",
            "| Founder scenario | Rank 1 model | Rank 1 weight |",
            "|---|---|---:|",
        ]
    )
    for point in payload["sensitivity_points"]:
        best = point["results"][0]
        lines.append(
            f"| {point['founder_scenario_id']} | "
            f"{structural_model_key(best)} | {best['weight_percent']:.3f}% |"
        )

    lines.extend(["", "## Model Weights", ""])
    for point in payload["sensitivity_points"]:
        lines.append(f"### {point['founder_scenario_id']}")
        lines.append("")
        lines.append("| Rank | Model | logL | Weight | Δweight vs baseline |")
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
    return Path("reproduced") / f"{run_config.stem}_founder_prior_sensitivity"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate pedigree likelihood sensitivity over founder priors"
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
        "--engine",
        choices=["factors", "enumeration"],
        default="factors",
    )
    parser.add_argument(
        "--out_dir",
        default=None,
        help="Output directory. Defaults to reproduced/<run_stem>_founder_prior_sensitivity",
    )
    args = parser.parse_args()

    run_config = Path(args.run_config)
    out_dir = Path(args.out_dir) if args.out_dir else default_output_dir(run_config)
    try:
        payload = run_founder_prior_sensitivity(
            run_config=run_config,
            founder_scenarios=args.founder_scenarios,
            baseline_founder_scenario=args.baseline_founder_scenario,
            engine=args.engine,
        )
        write_sensitivity_json(out_dir / "sensitivity.json", payload)
        write_sensitivity_report(out_dir / "sensitivity.md", payload)
    except Exception as exc:
        print(f"FOUNDER-PRIOR SENSITIVITY FAILED: {exc}")
        return 1

    print("=== Founder-Prior Sensitivity ===")
    print(f"Run: {payload['run_config']}")
    print(f"Engine: {payload['engine']}")
    print(f"Baseline founder scenario: {payload['baseline_founder_scenario_id']}")
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
