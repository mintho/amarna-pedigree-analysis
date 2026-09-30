#!/usr/bin/env python3
"""Leave-one-locus-out robustness analysis for pedigree likelihood runs."""

from __future__ import annotations

import argparse
import copy
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from .compare import add_model_comparison
    from .likelihood_enumeration import (
        load_json,
        load_observation_model,
        model_specs_from_run,
        pedigree_log_likelihood,
        resolve_ref,
        selected_observation_model_summary,
    )
    from .likelihood_factors import pedigree_log_likelihood_factor_graph
    from .validate import validate_run
except ImportError:
    from compare import add_model_comparison  # type: ignore[no-redef]
    from likelihood_enumeration import (  # type: ignore[no-redef]
        load_json,
        load_observation_model,
        model_specs_from_run,
        pedigree_log_likelihood,
        resolve_ref,
        selected_observation_model_summary,
    )
    from likelihood_factors import pedigree_log_likelihood_factor_graph  # type: ignore[no-redef]
    from validate import validate_run  # type: ignore[no-redef]


def model_key(row: dict[str, Any]) -> str:
    return (
        f"{row['pedigree']} / "
        f"{row['founder_scenario']} / "
        f"{row['identity_hypothesis']}"
    )


def subset_observations(observations: dict[str, Any], loci: list[str]) -> dict[str, Any]:
    subset = copy.deepcopy(observations)
    subset["loci"] = list(loci)
    locus_set = set(loci)
    subset["observations"] = {
        pid: {
            locus: value
            for locus, value in by_locus.items()
            if locus in locus_set
        }
        for pid, by_locus in observations.get("observations", {}).items()
    }
    return subset


def subset_founder_scenario(
    founder_scenario: dict[str, Any],
    loci: list[str],
) -> dict[str, Any]:
    subset = copy.deepcopy(founder_scenario)
    locus_set = set(loci)
    subset["allele_frequencies"] = {
        locus: frequencies
        for locus, frequencies in founder_scenario.get("allele_frequencies", {}).items()
        if locus in locus_set
    }
    subset["explicit_founder_priors"] = {
        pid: {
            locus: dist
            for locus, dist in by_locus.items()
            if locus in locus_set
        }
        for pid, by_locus in founder_scenario.get("explicit_founder_priors", {}).items()
        if any(locus in locus_set for locus in by_locus)
    }
    return subset


def collect_run_loci(run_path: Path, run: dict[str, Any]) -> list[str]:
    """Return the common locus list used by all model specs in a run."""
    loci_by_model: list[tuple[str, ...]] = []
    for spec in model_specs_from_run(run):
        observations_ref = spec.get("observations", run.get("observations"))
        observations = load_json(resolve_ref(run_path, observations_ref))
        loci = tuple(observations.get("loci", []))
        if not loci:
            raise ValueError(f"{observations_ref}: no loci available for jackknife")
        loci_by_model.append(loci)

    first = loci_by_model[0]
    for loci in loci_by_model[1:]:
        if loci != first:
            raise ValueError(
                "locus jackknife requires all model specs to use the same locus list"
            )
    if len(first) < 2:
        raise ValueError("locus jackknife requires at least two loci")
    return list(first)


def evaluate_model_for_loci(
    *,
    run_path: Path,
    run: dict[str, Any],
    spec: dict[str, str],
    loci: list[str],
    engine: str,
) -> dict[str, Any]:
    persons = load_json(resolve_ref(run_path, spec.get("persons", run.get("persons"))))
    observations = subset_observations(
        load_json(resolve_ref(run_path, spec.get("observations", run.get("observations")))),
        loci,
    )
    pedigree = load_json(resolve_ref(run_path, spec["pedigree"]))
    founder_scenario = subset_founder_scenario(
        load_json(resolve_ref(run_path, spec["founder_scenario"])),
        loci,
    )
    identity_hypothesis = load_json(resolve_ref(run_path, spec["identity_hypothesis"]))
    background_error = float(run.get("error_model", {}).get("background_error", 1e-3))
    error_model = run.get("error_model")
    observation_model = load_observation_model(run_path, run)

    if engine == "factors":
        log_likelihood, locus_details = pedigree_log_likelihood_factor_graph(
            persons=persons,
            observation_data=observations,
            pedigree=pedigree,
            founder_scenario=founder_scenario,
            identity_hypothesis=identity_hypothesis,
            background_error=background_error,
            error_model=error_model,
            observation_model=observation_model,
        )
    elif engine == "enumeration":
        log_likelihood, locus_details = pedigree_log_likelihood(
            persons=persons,
            observation_data=observations,
            pedigree=pedigree,
            founder_scenario=founder_scenario,
            identity_hypothesis=identity_hypothesis,
            background_error=background_error,
            error_model=error_model,
            observation_model=observation_model,
        )
    else:
        raise ValueError(f"unknown engine: {engine}")

    return {
        "pedigree": pedigree["id"],
        "structure_id": spec.get("structure_id", pedigree["id"]),
        "founder_scenario": founder_scenario["id"],
        "identity_hypothesis": identity_hypothesis["id"],
        "log_likelihood": log_likelihood,
        "observation_model": selected_observation_model_summary(
            error_model=error_model,
            observation_model=observation_model,
        ),
        "loci": locus_details,
    }


def evaluate_models_for_loci(
    *,
    run_path: Path,
    run: dict[str, Any],
    loci: list[str],
    engine: str,
) -> list[dict[str, Any]]:
    return add_model_comparison(
        [
            evaluate_model_for_loci(
                run_path=run_path,
                run=run,
                spec=spec,
                loci=loci,
                engine=engine,
            )
            for spec in model_specs_from_run(run)
        ]
    )


def model_changes(
    baseline: list[dict[str, Any]],
    jackknife_results: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    baseline_by_key = {model_key(row): row for row in baseline}
    changes = []
    for row in jackknife_results:
        key = model_key(row)
        base = baseline_by_key[key]
        changes.append(
            {
                "model": key,
                "baseline_rank": base["rank"],
                "jackknife_rank": row["rank"],
                "rank_change": row["rank"] - base["rank"],
                "baseline_weight_percent": base["weight_percent"],
                "jackknife_weight_percent": row["weight_percent"],
                "weight_percent_change": (
                    row["weight_percent"] - base["weight_percent"]
                ),
                "baseline_log_likelihood": base["log_likelihood"],
                "jackknife_log_likelihood": row["log_likelihood"],
                "log_likelihood_change": (
                    row["log_likelihood"] - base["log_likelihood"]
                ),
            }
        )
    return sorted(changes, key=lambda row: (abs(row["rank_change"]), abs(row["weight_percent_change"])), reverse=True)


def run_locus_jackknife(
    *,
    run_config: Path,
    engine: str = "factors",
) -> dict[str, Any]:
    run_path = run_config.resolve()
    validate_run(run_path)
    run = load_json(run_path)
    loci = collect_run_loci(run_path, run)
    baseline = evaluate_models_for_loci(
        run_path=run_path,
        run=run,
        loci=loci,
        engine=engine,
    )

    leave_one_out = []
    for omitted_locus in loci:
        retained_loci = [locus for locus in loci if locus != omitted_locus]
        results = evaluate_models_for_loci(
            run_path=run_path,
            run=run,
            loci=retained_loci,
            engine=engine,
        )
        leave_one_out.append(
            {
                "omitted_locus": omitted_locus,
                "retained_loci": retained_loci,
                "results": results,
                "changes": model_changes(baseline, results),
            }
        )

    return {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_config": str(run_path),
        "engine": engine,
        "loci": loci,
        "baseline": baseline,
        "leave_one_out": leave_one_out,
    }


def default_output_dir(run_config: Path) -> Path:
    return Path("reproduced") / f"{run_config.stem}_locus_jackknife"


def write_jackknife_json(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Leave-one-locus-out jackknife for pedigree likelihood runs"
    )
    parser.add_argument("run_config", help="Path to a run configuration JSON file")
    parser.add_argument(
        "--engine",
        choices=["factors", "enumeration"],
        default="factors",
    )
    parser.add_argument(
        "--out_dir",
        default=None,
        help="Output directory. Defaults to reproduced/<run_stem>_locus_jackknife",
    )
    args = parser.parse_args()

    run_config = Path(args.run_config)
    out_dir = Path(args.out_dir) if args.out_dir else default_output_dir(run_config)
    try:
        payload = run_locus_jackknife(run_config=run_config, engine=args.engine)
        write_jackknife_json(out_dir / "jackknife.json", payload)
    except Exception as exc:
        print(f"JACKKNIFE FAILED: {exc}")
        return 1

    print("=== Locus Jackknife ===")
    print(f"Run: {payload['run_config']}")
    print(f"Engine: {payload['engine']}")
    print(f"Results: {out_dir / 'jackknife.json'}")
    print(f"Baseline best: {model_key(payload['baseline'][0])}")
    for row in payload["leave_one_out"]:
        print(
            f"omit {row['omitted_locus']}: best {model_key(row['results'][0])} "
            f"({row['results'][0]['weight_percent']:.3f}%)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
