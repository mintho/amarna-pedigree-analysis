#!/usr/bin/env python3
"""Exact full-pedigree likelihood by genotype enumeration.

This is the deliberately simple reference implementation. It sums over all
latent genotypes for each locus and therefore scales poorly for large
pedigrees. Later factor-graph code should be checked against this module.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import math
from pathlib import Path
from typing import Any

try:  # Package import when used by tests.
    from .genetics import (
        Genotype,
        founder_prior_from_allele_frequencies,
        genotype_prior_from_mapping,
        mendelian_child_distribution,
        observation_likelihood_from_observation_model,
        parse_genotype_key,
        unordered_genotypes,
    )
    from .validate import validate_run
except ImportError:  # Script import when run as `python3 src/...`.
    from genetics import (  # type: ignore[no-redef]
        Genotype,
        founder_prior_from_allele_frequencies,
        genotype_prior_from_mapping,
        mendelian_child_distribution,
        observation_likelihood_from_observation_model,
        parse_genotype_key,
        unordered_genotypes,
    )
    from validate import validate_run  # type: ignore[no-redef]


Relationship = dict[str, str]


class LikelihoodError(Exception):
    """Raised when likelihood inputs cannot be evaluated."""


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise LikelihoodError(f"{path} must contain a JSON object")
    return data


def resolve_ref(base_file: Path, ref: str) -> Path:
    return (base_file.parent / ref).resolve()


def load_observation_model(run_path: Path, run: dict[str, Any]) -> dict[str, Any]:
    model = copy.deepcopy(run.get("observation_model"))
    if model is None:
        raise LikelihoodError("observation_model is required")
    if "background_frequency_file" in model:
        source = load_json(resolve_ref(run_path, model["background_frequency_file"]))
        model["background_allele_frequencies"] = source["allele_frequencies"]
    return model


def background_frequencies_for_locus(
    observation_model: dict[str, Any], locus: str, founder_frequencies: dict[str, float],
) -> dict[str, float]:
    # Freeze the background independently of varied founder support.
    if "background_allele_frequencies" not in observation_model:
        return founder_frequencies
    frequencies = observation_model["background_allele_frequencies"]
    if locus not in frequencies:
        raise LikelihoodError(f"Missing fixed background frequencies for {locus}")
    return frequencies[locus]


def logsumexp(values: list[float]) -> float:
    if not values:
        return float("-inf")
    max_value = max(values)
    if max_value == float("-inf"):
        return float("-inf")
    return max_value + math.log(sum(math.exp(value - max_value) for value in values))


def relationship_parent_map(relationships: list[Relationship]) -> dict[str, tuple[str, ...]]:
    parents: dict[str, tuple[str, ...]] = {}
    for rel in relationships:
        if "father" in rel and "mother" in rel:
            parents[rel["child"]] = (rel["father"], rel["mother"])
        elif "parent" in rel:
            parents[rel["child"]] = (rel["parent"],)
    return parents


def collect_locus_alleles(
    locus: str,
    person_ids: list[str],
    observations: dict[str, Any],
    allele_frequencies: dict[str, float],
    explicit_founder_priors: dict[str, Any] | None = None,
) -> list[str]:
    alleles = set(allele_frequencies)
    for pid in person_ids:
        obs = observations.get(pid, {}).get(locus, {})
        for allele in obs.get("alleles", []) or []:
            alleles.add(allele)
        explicit = (explicit_founder_priors or {}).get(pid, {}).get(locus, {})
        for key in explicit:
            alleles.update(parse_genotype_key(key))
    if not alleles:
        raise LikelihoodError(f"{locus}: no alleles available")
    return sorted(alleles)


def founder_prior_for_person(
    founder: str,
    locus: str,
    allele_frequencies: dict[str, float],
    explicit_founder_priors: dict[str, Any] | None = None,
) -> dict[Genotype, float]:
    explicit = (explicit_founder_priors or {}).get(founder, {}).get(locus)
    if explicit:
        return genotype_prior_from_mapping(explicit)
    return founder_prior_from_allele_frequencies(allele_frequencies)


def relationship_person_ids(relationships: list[Relationship]) -> set[str]:
    result: set[str] = set()
    for rel in relationships:
        result.add(rel["child"])
        if "father" in rel and "mother" in rel:
            result.update([rel["father"], rel["mother"]])
        elif "parent" in rel:
            result.add(rel["parent"])
    return result


def single_parent_child_probability(
    *,
    parent_genotype: Genotype,
    child_genotype: Genotype,
    allele_frequencies: dict[str, float],
) -> float:
    """Return P(child genotype | one known parent, unknown founder co-parent)."""
    unknown_parent_prior = founder_prior_from_allele_frequencies(allele_frequencies)
    probability = 0.0
    for unknown_genotype, unknown_probability in unknown_parent_prior.items():
        child_dist = mendelian_child_distribution(parent_genotype, unknown_genotype)
        probability += unknown_probability * child_dist.get(child_genotype, 0.0)
    return probability


def merge_observation(
    *,
    existing: dict[str, Any],
    incoming: dict[str, Any],
    context: str,
) -> dict[str, Any]:
    existing_alleles = list(existing.get("alleles", []) or [])
    incoming_alleles = list(incoming.get("alleles", []) or [])
    if not incoming_alleles:
        return existing
    if not existing_alleles:
        return dict(incoming)
    if sorted(existing_alleles) != sorted(incoming_alleles):
        raise LikelihoodError(
            f"{context}: conflicting observations {existing_alleles!r} vs {incoming_alleles!r}"
        )
    return existing


def effective_observation_error_model(
    global_error_model: dict[str, Any] | None,
    observation: dict[str, Any],
) -> dict[str, Any] | None:
    """Return run-level error settings with observation-level overrides applied."""
    local_error_model = observation.get("error_model") if observation else None
    if local_error_model is None:
        return global_error_model
    if not isinstance(local_error_model, dict):
        raise LikelihoodError("observation error_model must be an object")

    effective = dict(global_error_model or {})
    effective.update(local_error_model)
    return effective


def selected_observation_model_summary(
    *,
    error_model: dict[str, Any] | None = None,
    observation_model: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return the run-level observation model used for audit output."""
    if observation_model is None:
        raise LikelihoodError("observation_model is required")
    summary = {"type": observation_model.get("type")}
    for key in ("dropout", "background_error", "background_source", "background_frequency_file", "missing"):
        if key in observation_model:
            summary[key] = observation_model[key]
    return summary


def apply_identity_hypothesis(
    *,
    person_ids: list[str],
    relationships: list[Relationship],
    observations: dict[str, Any],
    identity_hypothesis: dict[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    """Apply sample-to-person identity assignments to observations.

    An assignment {"sample": "SampleA", "person": "PersonB"} means that the
    STR observations recorded under the sample node are evidence on the latent
    genotype of the named person. If the sample is not used in any relationship,
    it is removed from the variable set after its observations are moved.
    """
    assignments = identity_hypothesis.get("assignments", [])
    if not assignments:
        return list(person_ids), copy.deepcopy(observations)

    effective_person_ids = list(person_ids)
    effective_observations = copy.deepcopy(observations)
    relationship_nodes = relationship_person_ids(relationships)

    for assignment in assignments:
        sample = assignment["sample"]
        person = assignment["person"]
        if sample == person:
            continue
        if sample in relationship_nodes:
            raise LikelihoodError(
                f"identity sample {sample!r} is used in pedigree relationships; "
                "non-trivial identity merging for relationship nodes is not supported yet"
            )

        sample_observations = effective_observations.get(sample, {})
        person_observations = effective_observations.setdefault(person, {})
        for locus, sample_observation in sample_observations.items():
            person_observations[locus] = merge_observation(
                existing=person_observations.get(locus, {}),
                incoming=sample_observation,
                context=f"identity {sample}->{person} at {locus}",
            )

        effective_observations.pop(sample, None)
        if sample in effective_person_ids:
            effective_person_ids.remove(sample)

    return effective_person_ids, effective_observations


def locus_log_likelihood(
    *,
    person_ids: list[str],
    relationships: list[Relationship],
    observations: dict[str, Any],
    locus: str,
    allele_frequencies: dict[str, float],
    explicit_founder_priors: dict[str, Any] | None = None,
    observation_model: dict[str, Any] | None = None,
) -> tuple[float, dict[str, Any]]:
    """Compute exact log P(observations at one locus | pedigree)."""
    parent_map = relationship_parent_map(relationships)
    founders = [pid for pid in person_ids if pid not in parent_map]
    alleles = collect_locus_alleles(
        locus,
        person_ids,
        observations,
        allele_frequencies,
        explicit_founder_priors,
    )
    genotype_states = unordered_genotypes(alleles)
    founder_priors = {
        founder: founder_prior_for_person(
            founder,
            locus,
            allele_frequencies,
            explicit_founder_priors,
        )
        for founder in founders
    }

    log_terms: list[float] = []
    assignments_checked = 0
    assignments_nonzero = 0

    for genotype_tuple in itertools.product(genotype_states, repeat=len(person_ids)):
        assignments_checked += 1
        genotypes = dict(zip(person_ids, genotype_tuple))
        logp = 0.0
        impossible = False

        for founder in founders:
            prior = founder_priors[founder].get(genotypes[founder], 0.0)
            if prior <= 0.0:
                impossible = True
                break
            logp += math.log(prior)
        if impossible:
            continue

        for rel in relationships:
            child_gt = genotypes[rel["child"]]
            if "father" in rel and "mother" in rel:
                father_gt = genotypes[rel["father"]]
                mother_gt = genotypes[rel["mother"]]
                transmission = mendelian_child_distribution(father_gt, mother_gt)
                probability = transmission.get(child_gt, 0.0)
            else:
                probability = single_parent_child_probability(
                    parent_genotype=genotypes[rel["parent"]],
                    child_genotype=child_gt,
                    allele_frequencies=allele_frequencies,
                )
            if probability <= 0.0:
                impossible = True
                break
            logp += math.log(probability)
        if impossible:
            continue

        for pid in person_ids:
            observation = observations.get(pid, {}).get(locus, {})
            observed_alleles = observation.get("alleles", [])
            observation_error_model = effective_observation_error_model(
                observation_model,
                observation,
            )
            if observation_model is None:
                raise LikelihoodError("observation_model is required")
            probability = observation_likelihood_from_observation_model(
                observed_alleles,
                genotypes[pid],
                observation_error_model or observation_model,
                background_frequencies_for_locus(observation_model, locus, allele_frequencies),
            )
            if probability <= 0.0:
                impossible = True
                break
            logp += math.log(probability)
        if impossible:
            continue

        assignments_nonzero += 1
        log_terms.append(logp)

    log_likelihood = logsumexp(log_terms)
    details = {
        "locus": locus,
        "alleles": alleles,
        "observation_model": selected_observation_model_summary(
            observation_model=observation_model,
        ),
        "genotype_states": len(genotype_states),
        "persons": len(person_ids),
        "founders": founders,
        "assignments_checked": assignments_checked,
        "assignments_nonzero": assignments_nonzero,
        "log_likelihood": log_likelihood,
    }
    return log_likelihood, details


def pedigree_log_likelihood(
    *,
    persons: dict[str, Any],
    observation_data: dict[str, Any],
    pedigree: dict[str, Any],
    founder_scenario: dict[str, Any],
    identity_hypothesis: dict[str, Any],
    observation_model: dict[str, Any] | None = None,
) -> tuple[float, list[dict[str, Any]]]:
    """Compute exact full-pedigree log-likelihood over all loci."""
    person_ids = [person["id"] for person in persons["persons"]]
    relationships = pedigree.get("relationships", [])
    observations = observation_data.get("observations", {})
    person_ids, observations = apply_identity_hypothesis(
        person_ids=person_ids,
        relationships=relationships,
        observations=observations,
        identity_hypothesis=identity_hypothesis,
    )
    loci = observation_data.get("loci", [])
    allele_frequencies_by_locus = founder_scenario.get("allele_frequencies", {})
    explicit_founder_priors = founder_scenario.get("explicit_founder_priors", {})

    locus_details: list[dict[str, Any]] = []
    total = 0.0
    for locus in loci:
        allele_frequencies = allele_frequencies_by_locus[locus]
        log_likelihood, details = locus_log_likelihood(
            person_ids=person_ids,
            relationships=relationships,
            observations=observations,
            locus=locus,
            allele_frequencies=allele_frequencies,
            explicit_founder_priors=explicit_founder_priors,
            observation_model=observation_model,
        )
        total += log_likelihood
        locus_details.append(details)
    return total, locus_details


def model_specs_from_run(run: dict[str, Any]) -> list[dict[str, str]]:
    """Return explicit model specs or expand cartesian-product run fields."""
    if "models" in run:
        return [
            {
                "pedigree": model["pedigree"],
                "founder_scenario": model["founder_scenario"],
                "identity_hypothesis": model["identity_hypothesis"],
                **({"persons": model["persons"]} if "persons" in model else {}),
                **({"observations": model["observations"]} if "observations" in model else {}),
                **({"structure_id": model["structure_id"]} if "structure_id" in model else {}),
            }
            for model in run["models"]
        ]

    return [
        {
            "pedigree": pedigree,
            "founder_scenario": founder_scenario,
            "identity_hypothesis": identity_hypothesis,
        }
        for pedigree in run["pedigrees"]
        for founder_scenario in run["founder_scenarios"]
        for identity_hypothesis in run["identity_hypotheses"]
    ]


def evaluate_run(run_path: Path) -> list[dict[str, Any]]:
    """Validate and evaluate every model combination referenced by a run file."""
    validate_run(run_path)
    run = load_json(run_path)

    specs = model_specs_from_run(run)
    observation_model = load_observation_model(run_path, run)
    observation_model_summary = selected_observation_model_summary(
        observation_model=observation_model,
    )

    results: list[dict[str, Any]] = []
    for spec in specs:
        persons = load_json(resolve_ref(run_path, spec.get("persons", run.get("persons"))))
        observations = load_json(resolve_ref(run_path, spec.get("observations", run.get("observations"))))
        pedigree = load_json(resolve_ref(run_path, spec["pedigree"]))
        founder_scenario = load_json(resolve_ref(run_path, spec["founder_scenario"]))
        identity_hypothesis = load_json(resolve_ref(run_path, spec["identity_hypothesis"]))
        log_likelihood, locus_details = pedigree_log_likelihood(
            persons=persons,
            observation_data=observations,
            pedigree=pedigree,
            founder_scenario=founder_scenario,
            identity_hypothesis=identity_hypothesis,
            observation_model=observation_model,
        )
        results.append(
            {
                "pedigree": pedigree["id"],
                "structure_id": spec.get("structure_id", pedigree["id"]),
                "founder_scenario": founder_scenario["id"],
                "identity_hypothesis": identity_hypothesis["id"],
                "log_likelihood": log_likelihood,
                "observation_model": observation_model_summary,
                "loci": locus_details,
            }
        )
    results.sort(key=lambda row: row["log_likelihood"], reverse=True)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Exact full-pedigree likelihood by genotype enumeration"
    )
    parser.add_argument(
        "run_config",
        nargs="?",
        default="Data/models/runs/scenario_c_primary_models_all_loci_generative_str.json",
        help="Path to a run configuration JSON file",
    )
    args = parser.parse_args()

    run_path = Path(args.run_config).resolve()
    try:
        results = evaluate_run(run_path)
    except Exception as exc:
        print(f"LIKELIHOOD FAILED: {exc}")
        return 1

    print("=== Exact Full-Pedigree Likelihoods ===")
    print(f"Run: {run_path}")
    print(f"{'Rank':>4}  {'logL':>14}  Pedigree / Founder / Identity")
    for rank, result in enumerate(results, start=1):
        print(
            f"{rank:>4}  {result['log_likelihood']:>14.6f}  "
            f"{result['pedigree']} / "
            f"{result['founder_scenario']} / "
            f"{result['identity_hypothesis']}"
        )
        for locus in result["loci"]:
            print(
                f"      {locus['locus']}: logL={locus['log_likelihood']:.6f}, "
                f"states={locus['genotype_states']}, "
                f"assignments={locus['assignments_nonzero']}/"
                f"{locus['assignments_checked']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
