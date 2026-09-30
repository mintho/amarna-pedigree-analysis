#!/usr/bin/env python3
"""Validate run configuration and referenced input files.

This module intentionally does not compute likelihoods. It checks that a run is
well-defined before any model code consumes it.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


class ValidationError(Exception):
    """Raised when input data fails validation."""


def load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except FileNotFoundError as exc:
        raise ValidationError(f"File not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValidationError(f"Invalid JSON in {path}: {exc}") from exc

    if not isinstance(data, dict):
        raise ValidationError(f"{path} must contain a JSON object")
    return data


def require_schema(data: dict[str, Any], path: Path, expected: int = 1) -> None:
    version = data.get("schema_version")
    if version != expected:
        raise ValidationError(
            f"{path}: schema_version must be {expected}, got {version!r}"
        )


def resolve_ref(base_file: Path, ref: str) -> Path:
    return (base_file.parent / ref).resolve()


def assert_probability(value: Any, context: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ValidationError(f"{context}: probability must be a number")
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValidationError(f"{context}: probability must be finite and >= 0")
    return value


def assert_distribution(dist: Any, context: str, tol: float = 1e-6) -> None:
    if not isinstance(dist, dict) or not dist:
        raise ValidationError(f"{context}: distribution must be a non-empty object")

    total = 0.0
    for key, value in dist.items():
        if not isinstance(key, str) or not key:
            raise ValidationError(f"{context}: distribution keys must be strings")
        total += assert_probability(value, f"{context}.{key}")

    if abs(total - 1.0) > tol:
        raise ValidationError(f"{context}: probabilities must sum to 1, got {total}")


def assert_genotype_distribution(dist: Any, context: str, tol: float = 1e-6) -> None:
    assert_distribution(dist, context, tol=tol)
    for key in dist:
        parts = key.split("/") if isinstance(key, str) else []
        if len(parts) != 2 or not parts[0] or not parts[1]:
            raise ValidationError(
                f"{context}: genotype keys must have form 'allele_a/allele_b', got {key!r}"
            )


def validate_error_model_rates(
    error_model: Any,
    context: str,
    required_keys: tuple[str, ...] = (),
) -> None:
    if not isinstance(error_model, dict):
        raise ValidationError(f"{context}: error_model must be an object")
    allowed_keys = {"background_error", "dropout", "miscall", "contamination", "background_error"}
    for key in required_keys:
        if key not in error_model:
            raise ValidationError(f"{context}: error_model.{key} is required")
    for key, value in error_model.items():
        if key not in allowed_keys:
            continue
        probability = assert_probability(value, f"{context}: error_model.{key}")
        if probability >= 1.0:
            raise ValidationError(f"{context}: error_model.{key} must be < 1")


def validate_persons(data: dict[str, Any], path: Path) -> set[str]:
    require_schema(data, path)
    persons = data.get("persons")
    if not isinstance(persons, list) or not persons:
        raise ValidationError(f"{path}: persons must be a non-empty list")

    ids: set[str] = set()
    allowed_kinds = {"historical", "sample", "latent"}
    for index, person in enumerate(persons):
        context = f"{path}: persons[{index}]"
        if not isinstance(person, dict):
            raise ValidationError(f"{context} must be an object")

        pid = person.get("id")
        if not isinstance(pid, str) or not pid:
            raise ValidationError(f"{context}.id must be a non-empty string")
        if pid in ids:
            raise ValidationError(f"{path}: duplicate person id {pid!r}")
        ids.add(pid)

        label = person.get("label")
        if label is not None and not isinstance(label, str):
            raise ValidationError(f"{context}.label must be a string")

        kind = person.get("kind")
        if kind not in allowed_kinds:
            raise ValidationError(
                f"{context}.kind must be one of {sorted(allowed_kinds)}, got {kind!r}"
            )

    return ids


def validate_observations(
    data: dict[str, Any], path: Path, person_ids: set[str]
) -> set[str]:
    require_schema(data, path)

    loci = data.get("loci")
    if not isinstance(loci, list) or not loci:
        raise ValidationError(f"{path}: loci must be a non-empty list")
    if any(not isinstance(locus, str) or not locus for locus in loci):
        raise ValidationError(f"{path}: every locus must be a non-empty string")
    if len(set(loci)) != len(loci):
        raise ValidationError(f"{path}: loci must be unique")
    locus_set = set(loci)

    observations = data.get("observations")
    if not isinstance(observations, dict):
        raise ValidationError(f"{path}: observations must be an object")

    for pid, by_locus in observations.items():
        if pid not in person_ids:
            raise ValidationError(f"{path}: observation references unknown person {pid!r}")
        if not isinstance(by_locus, dict):
            raise ValidationError(f"{path}: observations.{pid} must be an object")

        for locus, obs in by_locus.items():
            context = f"{path}: observations.{pid}.{locus}"
            if locus not in locus_set:
                raise ValidationError(f"{context}: unknown locus")
            if not isinstance(obs, dict):
                raise ValidationError(f"{context} must be an object")

            alleles = obs.get("alleles")
            if not isinstance(alleles, list) or len(alleles) > 2:
                raise ValidationError(f"{context}.alleles must be a list of length 0..2")
            for allele in alleles:
                if not isinstance(allele, str) or not allele:
                    raise ValidationError(
                        f"{context}.alleles must contain non-empty strings"
                    )

            source = obs.get("source")
            if source is not None and not isinstance(source, str):
                raise ValidationError(f"{context}.source must be a string")
            if "error_model" in obs:
                validate_error_model_rates(
                    obs["error_model"],
                    context,
                )

    return locus_set


def validate_pedigree(
    data: dict[str, Any], path: Path, person_ids: set[str]
) -> set[str]:
    require_schema(data, path)
    model_id = data.get("id")
    if not isinstance(model_id, str) or not model_id:
        raise ValidationError(f"{path}: id must be a non-empty string")

    relationships = data.get("relationships")
    if not isinstance(relationships, list):
        raise ValidationError(f"{path}: relationships must be a list")

    children: set[str] = set()
    for index, rel in enumerate(relationships):
        context = f"{path}: relationships[{index}]"
        if not isinstance(rel, dict):
            raise ValidationError(f"{context} must be an object")

        child = rel.get("child")
        father = rel.get("father")
        mother = rel.get("mother")
        parent = rel.get("parent")

        has_full_parent_pair = father is not None or mother is not None
        has_single_parent = parent is not None
        if has_full_parent_pair == has_single_parent:
            raise ValidationError(
                f"{context}: relationship must define either father/mother or parent"
            )
        if has_full_parent_pair and (father is None or mother is None):
            raise ValidationError(f"{context}: father and mother must be provided together")

        roles = [("child", child)]
        if has_full_parent_pair:
            roles.extend([("father", father), ("mother", mother)])
        else:
            roles.append(("parent", parent))

        for role, pid in roles:
            if not isinstance(pid, str) or not pid:
                raise ValidationError(f"{context}.{role} must be a non-empty string")
            if pid not in person_ids:
                raise ValidationError(f"{context}.{role} references unknown person {pid!r}")

        if child in children:
            raise ValidationError(f"{path}: child {child!r} has multiple parent pairs")
        if has_full_parent_pair:
            if child == father or child == mother:
                raise ValidationError(f"{context}: child cannot be its own parent")
            if father == mother:
                raise ValidationError(f"{context}: father and mother must differ")
        elif child == parent:
            raise ValidationError(f"{context}: child cannot be its own parent")
        children.add(child)

    return children


def validate_founder_scenario(
    data: dict[str, Any], path: Path, locus_set: set[str], person_ids: set[str]
) -> None:
    require_schema(data, path)
    scenario_id = data.get("id")
    if not isinstance(scenario_id, str) or not scenario_id:
        raise ValidationError(f"{path}: id must be a non-empty string")

    allele_freqs = data.get("allele_frequencies")
    if not isinstance(allele_freqs, dict):
        raise ValidationError(f"{path}: allele_frequencies must be an object")
    missing = sorted(locus_set - set(allele_freqs))
    if missing:
        raise ValidationError(f"{path}: missing allele frequencies for loci {missing}")
    unknown = sorted(set(allele_freqs) - locus_set)
    if unknown:
        raise ValidationError(f"{path}: allele frequencies contain unknown loci {unknown}")

    for locus, dist in allele_freqs.items():
        assert_distribution(dist, f"{path}: allele_frequencies.{locus}")

    explicit = data.get("explicit_founder_priors", {})
    if not isinstance(explicit, dict):
        raise ValidationError(f"{path}: explicit_founder_priors must be an object")

    for pid, by_locus in explicit.items():
        if pid not in person_ids:
            raise ValidationError(
                f"{path}: explicit prior references unknown person {pid!r}"
            )
        if not isinstance(by_locus, dict):
            raise ValidationError(f"{path}: explicit_founder_priors.{pid} must be object")
        for locus, dist in by_locus.items():
            if locus not in locus_set:
                raise ValidationError(
                    f"{path}: explicit prior for {pid} references unknown locus {locus!r}"
                )
            assert_genotype_distribution(
                dist, f"{path}: explicit_founder_priors.{pid}.{locus}"
            )


def validate_identity_hypothesis(
    data: dict[str, Any], path: Path, person_ids: set[str]
) -> None:
    require_schema(data, path)
    hypothesis_id = data.get("id")
    if not isinstance(hypothesis_id, str) or not hypothesis_id:
        raise ValidationError(f"{path}: id must be a non-empty string")

    assignments = data.get("assignments")
    if not isinstance(assignments, list):
        raise ValidationError(f"{path}: assignments must be a list")

    seen_samples: set[str] = set()
    for index, assignment in enumerate(assignments):
        context = f"{path}: assignments[{index}]"
        if not isinstance(assignment, dict):
            raise ValidationError(f"{context} must be an object")
        sample = assignment.get("sample")
        person = assignment.get("person")
        for role, pid in (("sample", sample), ("person", person)):
            if not isinstance(pid, str) or not pid:
                raise ValidationError(f"{context}.{role} must be a non-empty string")
            if pid not in person_ids:
                raise ValidationError(f"{context}.{role} references unknown person {pid!r}")
        if sample in seen_samples:
            raise ValidationError(f"{path}: sample {sample!r} assigned more than once")
        seen_samples.add(sample)


def validate_run(path: Path) -> list[str]:
    run = load_json(path)
    require_schema(run, path)

    primaryd = run.get("id")
    if not isinstance(primaryd, str) or not primaryd:
        raise ValidationError(f"{path}: id must be a non-empty string")

    if "models" in run and ("persons" not in run or "observations" not in run):
        messages = []
        for index, spec in enumerate(validate_model_specs(run, path)):
            persons_path = resolve_ref(path, require_string(spec, "persons", path))
            observations_path = resolve_ref(path, require_string(spec, "observations", path))
            model_messages = validate_model_references(
                path=path,
                persons_path=persons_path,
                observations_path=observations_path,
                pedigree_refs=[spec["pedigree"]],
                founder_refs=[spec["founder_scenario"]],
                identity_refs=[spec["identity_hypothesis"]],
            )
            messages.append(f"model[{index}] {spec.get('label') or spec['pedigree']}: ok")
            messages.extend(f"  {message}" for message in model_messages)
        validate_run_observation_settings(run, path)
        messages.append(f"observation model type: {run_observation_model_type(run)}")
        messages.append("observation model: ok")
        return messages

    persons_path = resolve_ref(path, require_string(run, "persons", path))
    observations_path = resolve_ref(path, require_string(run, "observations", path))

    if "models" in run:
        model_specs = validate_model_specs(run, path)
        pedigree_refs = [spec["pedigree"] for spec in model_specs]
        founder_refs = [spec["founder_scenario"] for spec in model_specs]
        identity_refs = [spec["identity_hypothesis"] for spec in model_specs]
    else:
        pedigree_refs = require_string_list(run, "pedigrees", path)
        founder_refs = require_string_list(run, "founder_scenarios", path)
        identity_refs = require_string_list(run, "identity_hypotheses", path)

    messages = validate_model_references(
        path=path,
        persons_path=persons_path,
        observations_path=observations_path,
        pedigree_refs=pedigree_refs,
        founder_refs=founder_refs,
        identity_refs=identity_refs,
    )

    validate_run_observation_settings(run, path)
    messages.append(f"observation model type: {run_observation_model_type(run)}")
    messages.append("observation model: ok")

    return messages


def validate_model_references(
    *,
    path: Path,
    persons_path: Path,
    observations_path: Path,
    pedigree_refs: list[str],
    founder_refs: list[str],
    identity_refs: list[str],
) -> list[str]:
    persons = load_json(persons_path)
    person_ids = validate_persons(persons, persons_path)

    observations = load_json(observations_path)
    locus_set = validate_observations(observations, observations_path, person_ids)

    messages = [
        f"persons: {len(person_ids)}",
        f"loci: {len(locus_set)}",
    ]

    for ref in pedigree_refs:
        pedigree_path = resolve_ref(path, ref)
        pedigree = load_json(pedigree_path)
        children = validate_pedigree(pedigree, pedigree_path, person_ids)
        messages.append(f"pedigree {pedigree.get('id')}: {len(children)} relationships")

    for ref in founder_refs:
        scenario_path = resolve_ref(path, ref)
        scenario = load_json(scenario_path)
        validate_founder_scenario(scenario, scenario_path, locus_set, person_ids)
        messages.append(f"founder scenario {scenario.get('id')}: ok")

    for ref in identity_refs:
        identity_path = resolve_ref(path, ref)
        identity = load_json(identity_path)
        validate_identity_hypothesis(identity, identity_path, person_ids)
        messages.append(f"identity hypothesis {identity.get('id')}: ok")

    return messages


def require_string(data: dict[str, Any], key: str, path: Path) -> str:
    value = data.get(key)
    if not isinstance(value, str) or not value:
        raise ValidationError(f"{path}: {key} must be a non-empty string")
    return value


def require_string_list(data: dict[str, Any], key: str, path: Path) -> list[str]:
    value = data.get(key)
    if not isinstance(value, list) or not value:
        raise ValidationError(f"{path}: {key} must be a non-empty list")
    for index, item in enumerate(value):
        if not isinstance(item, str) or not item:
            raise ValidationError(f"{path}: {key}[{index}] must be a non-empty string")
    return value


def validate_model_specs(run: dict[str, Any], path: Path) -> list[dict[str, str]]:
    models = run.get("models")
    if not isinstance(models, list) or not models:
        raise ValidationError(f"{path}: models must be a non-empty list")

    result: list[dict[str, str]] = []
    for index, model in enumerate(models):
        context = f"{path}: models[{index}]"
        if not isinstance(model, dict):
            raise ValidationError(f"{context} must be an object")
        spec = {
            "pedigree": require_string(model, "pedigree", Path(context)),
            "founder_scenario": require_string(model, "founder_scenario", Path(context)),
            "identity_hypothesis": require_string(model, "identity_hypothesis", Path(context)),
        }
        if "persons" in model:
            spec["persons"] = require_string(model, "persons", Path(context))
        if "observations" in model:
            spec["observations"] = require_string(model, "observations", Path(context))
        label = model.get("label")
        if label is not None and not isinstance(label, str):
            raise ValidationError(f"{context}.label must be a string")
        if label is not None:
            spec["label"] = label
        result.append(spec)
    return result


def validate_error_model(error_model: Any, path: Path) -> None:
    validate_error_model_rates(
        error_model,
        str(path),
        required_keys=("background_error", "dropout"),
    )


def validate_observation_model(observation_model: Any, path: Path) -> None:
    if not isinstance(observation_model, dict):
        raise ValidationError(f"{path}: observation_model must be an object")
    model_type = observation_model.get("type")
    if model_type != "generative_str":
        raise ValidationError(
            f"{path}: observation_model.type must be 'generative_str', "
            f"got {model_type!r}"
        )

    for key in ("dropout", "background_error"):
        if key not in observation_model:
            raise ValidationError(f"{path}: observation_model.{key} is required")
        probability = assert_probability(
            observation_model[key],
            f"{path}: observation_model.{key}",
        )
        if probability >= 1.0:
            raise ValidationError(f"{path}: observation_model.{key} must be < 1")

    background_source = observation_model.get("background_source")
    if background_source not in (None, "run_locus_allele_frequencies", "fixed_reference_frequencies"):
        raise ValidationError(
            f"{path}: observation_model.background_source must be "
            "'run_locus_allele_frequencies' or 'fixed_reference_frequencies'"
        )
    if background_source == "fixed_reference_frequencies":
        ref = observation_model.get("background_frequency_file")
        if not isinstance(ref, str):
            raise ValidationError(f"{path}: fixed background requires background_frequency_file")
        source = load_json(resolve_ref(path, ref))
        for locus, frequencies in source.get("allele_frequencies", {}).items():
            if any(not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0
                   for value in frequencies.values()) or not math.isclose(sum(frequencies.values()), 1.0, abs_tol=1e-9):
                raise ValidationError(f"{path}: invalid fixed background for {locus}")

    missing = observation_model.get("missing")
    if missing is not None and missing != "condition_on_call_available":
        raise ValidationError(
            f"{path}: observation_model.missing must be "
            "'condition_on_call_available'"
        )


def validate_run_observation_settings(run: dict[str, Any], path: Path) -> None:
    has_observation_model = "observation_model" in run
    if has_observation_model:
        validate_observation_model(run["observation_model"], path)
        return
    raise ValidationError(f"{path}: observation_model is required")


def run_observation_model_type(run: dict[str, Any]) -> str:
    if "observation_model" in run:
        observation_model = run["observation_model"]
        if isinstance(observation_model, dict):
            return str(observation_model.get("type", "missing_observation_model_type"))
        return "invalid_observation_model"
    return "missing_observation_model"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate a pedigree likelihood run")
    parser.add_argument(
        "run_config",
        nargs="?",
        default="Data/models/runs/scenario_c_primary_models_all_loci_generative_str.json",
        help="Path to a run configuration JSON file",
    )
    args = parser.parse_args()

    path = Path(args.run_config).resolve()
    try:
        messages = validate_run(path)
    except ValidationError as exc:
        print(f"VALIDATION FAILED: {exc}")
        return 1

    print(f"VALIDATION OK: {path}")
    for message in messages:
        print(f"- {message}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
