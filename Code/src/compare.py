#!/usr/bin/env python3
"""Compare likelihood results and write reproducible reports."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from .likelihood_enumeration import evaluate_run
    from .likelihood_factors import evaluate_run_factor_graph
except ImportError:
    from likelihood_enumeration import evaluate_run  # type: ignore[no-redef]
    from likelihood_factors import evaluate_run_factor_graph  # type: ignore[no-redef]


def add_model_comparison(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Add rank, delta-logL, and normalized likelihood weight to results."""
    if not results:
        return []

    ranked = sorted(results, key=lambda row: row["log_likelihood"], reverse=True)
    best_log_likelihood = ranked[0]["log_likelihood"]
    structures = {}
    for row in ranked:
        key = row.get("structure_id", row["pedigree"])
        if key in structures and not math.isclose(structures[key], row["log_likelihood"], abs_tol=1e-9, rel_tol=0):
            raise ValueError(f"Different likelihoods assigned to shared structure {key}")
        structures[key] = row["log_likelihood"]
    structure_ranks = {key: rank for rank, key in enumerate(structures, 1)}
    weight_total = sum(math.exp(value - best_log_likelihood) for value in structures.values()) or 1.0

    compared: list[dict[str, Any]] = []
    for row in ranked:
        key = row.get("structure_id", row["pedigree"])
        raw_weight = math.exp(row["log_likelihood"] - best_log_likelihood)
        enriched = dict(row)
        enriched["structure_id"] = key
        enriched["weight_scope"] = "unique_scored_structure"
        enriched["rank"] = structure_ranks[key]
        enriched["delta_log_likelihood"] = row["log_likelihood"] - best_log_likelihood
        enriched["weight"] = raw_weight / weight_total
        enriched["weight_percent"] = 100.0 * enriched["weight"]
        compared.append(enriched)

    return compared


def bayes_factor(log_likelihood_a: float, log_likelihood_b: float) -> float:
    """Return BF(a:b) from two log-likelihoods."""
    return math.exp(log_likelihood_a - log_likelihood_b)


def write_results_json(
    output_path: Path,
    *,
    run_config: Path,
    engine: str,
    compared_results: list[dict[str, Any]],
) -> None:
    payload = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_config": str(run_config),
        "engine": engine,
        "normalisation": "equal prior mass per unique structure_id; historical aliases share a structural weight",
        "results": compared_results,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def markdown_report(
    *,
    run_config: Path,
    engine: str,
    compared_results: list[dict[str, Any]],
) -> str:
    lines = [
        "# Pedigree Likelihood Report",
        "",
        f"Run configuration: `{run_config}`",
        f"Engine: `{engine}`",
    ]
    if compared_results and compared_results[0].get("observation_model"):
        lines.append(
            "Observation model: "
            f"`{json.dumps(compared_results[0]['observation_model'], sort_keys=True)}`"
        )
    lines.extend(
        [
            "",
            "## Ranking",
            "",
            "| Rank | Model | logL | ΔlogL | Weight |",
            "|---:|---|---:|---:|---:|",
        ]
    )

    for row in compared_results:
        model = (
            f"{row['pedigree']} / "
            f"{row['founder_scenario']} / "
            f"{row['identity_hypothesis']}"
        )
        lines.append(
            f"| {row['rank']} | {model} | "
            f"{row['log_likelihood']:.6f} | "
            f"{row['delta_log_likelihood']:.6f} | "
            f"{row['weight_percent']:.3f}% |"
        )

    if len(compared_results) >= 2:
        best = compared_results[0]
        lines.extend(
            [
                "",
                "## Support Ratios (LR for sensitivity settings; BF for population-prior comparisons)",
                "",
                "| Comparison | Ratio |",
                "|---|---:|",
            ]
        )
        for row in compared_results[1:]:
            model = (
                f"{best['pedigree']} / {best['founder_scenario']} / "
                f"{best['identity_hypothesis']} vs. "
                f"{row['pedigree']} / {row['founder_scenario']} / "
                f"{row['identity_hypothesis']}"
            )
            bf = bayes_factor(best["log_likelihood"], row["log_likelihood"])
            lines.append(f"| {model} | {bf:.6g} |")

    lines.extend(
        [
            "",
            "## Locus Details",
            "",
        ]
    )
    for row in compared_results:
        model = (
            f"{row['pedigree']} / "
            f"{row['founder_scenario']} / "
            f"{row['identity_hypothesis']}"
        )
        lines.append(f"### Rank {row['rank']}: {model}")
        lines.append("")
        lines.append("| Locus | logL | Alleles | Genotype states | Work summary |")
        lines.append("|---|---:|---|---:|---:|")
        for locus in row["loci"]:
            alleles = ", ".join(locus.get("alleles", [])) or "n/a"
            if "assignments_nonzero" in locus:
                work_summary = (
                    f"{locus['assignments_nonzero']}/"
                    f"{locus['assignments_checked']} assignments"
                )
            elif "max_factor_rows" in locus:
                work_summary = f"max factor rows {locus['max_factor_rows']}"
            else:
                work_summary = "n/a"
            lines.append(
                f"| {locus['locus']} | {locus['log_likelihood']:.6f} | "
                f"{alleles} | {locus['genotype_states']} | "
                f"{work_summary} |"
            )
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_report_md(
    output_path: Path,
    *,
    run_config: Path,
    engine: str,
    compared_results: list[dict[str, Any]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        markdown_report(
            run_config=run_config,
            engine=engine,
            compared_results=compared_results,
        ),
        encoding="utf-8",
    )


def default_output_dir(run_config: Path) -> Path:
    primaryd = run_config.stem
    return Path("reproduced") / primaryd


def evaluate_with_engine(run_config: Path, engine: str) -> list[dict[str, Any]]:
    if engine == "factors":
        return evaluate_run_factor_graph(run_config)
    if engine == "enumeration":
        return evaluate_run(run_config)
    raise ValueError(f"unknown engine: {engine}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate and compare exact pedigree likelihoods"
    )
    parser.add_argument(
        "run_config",
        nargs="?",
        default="Data/models/runs/scenario_c_primary_models_all_loci_generative_str.json",
        help="Path to a run configuration JSON file",
    )
    parser.add_argument(
        "--out_dir",
        default=None,
        help="Output directory. Defaults to reproduced/<run_config_stem>",
    )
    parser.add_argument(
        "--engine",
        choices=["factors", "enumeration"],
        default="factors",
        help="Likelihood engine. 'factors' is exact variable elimination and is the default.",
    )
    args = parser.parse_args()

    run_config = Path(args.run_config).resolve()
    out_dir = Path(args.out_dir) if args.out_dir else default_output_dir(Path(args.run_config))

    try:
        compared_results = add_model_comparison(evaluate_with_engine(run_config, args.engine))
        write_results_json(
            out_dir / "results.json",
            run_config=run_config,
            engine=args.engine,
            compared_results=compared_results,
        )
        write_report_md(
            out_dir / "report.md",
            run_config=run_config,
            engine=args.engine,
            compared_results=compared_results,
        )
    except Exception as exc:
        print(f"COMPARE FAILED: {exc}")
        return 1

    print("=== Model Comparison ===")
    print(f"Run: {run_config}")
    print(f"Engine: {args.engine}")
    print(f"Results: {out_dir / 'results.json'}")
    print(f"Report: {out_dir / 'report.md'}")
    print(f"{'Rank':>4}  {'logL':>14}  {'ΔlogL':>10}  {'Weight':>9}  Model")
    for row in compared_results:
        model = (
            f"{row['pedigree']} / "
            f"{row['founder_scenario']} / "
            f"{row['identity_hypothesis']}"
        )
        print(
            f"{row['rank']:>4}  {row['log_likelihood']:>14.6f}  "
            f"{row['delta_log_likelihood']:>10.6f}  "
            f"{row['weight_percent']:>8.3f}%  {model}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
