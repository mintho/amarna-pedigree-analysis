#!/usr/bin/env python3
"""Exact full-pedigree likelihood by variable elimination."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

try:
    from .genetics import (
        Genotype,
        mendelian_child_distribution,
        observation_likelihood_from_observation_model,
        unordered_genotypes,
    )
    from .likelihood_enumeration import (
        Relationship,
        background_frequencies_for_locus,
        load_observation_model,
        collect_locus_alleles,
        effective_observation_error_model,
        apply_identity_hypothesis,
        founder_prior_for_person,
        model_specs_from_run,
        relationship_parent_map,
        resolve_ref,
        selected_observation_model_summary,
        single_parent_child_probability,
    )
    from .validate import validate_run
except ImportError:
    from genetics import (  # type: ignore[no-redef]
        Genotype,
        mendelian_child_distribution,
        observation_likelihood_from_observation_model,
        unordered_genotypes,
    )
    from likelihood_enumeration import (  # type: ignore[no-redef]
        Relationship,
        background_frequencies_for_locus,
        load_observation_model,
        collect_locus_alleles,
        effective_observation_error_model,
        apply_identity_hypothesis,
        founder_prior_for_person,
        model_specs_from_run,
        relationship_parent_map,
        resolve_ref,
        selected_observation_model_summary,
        single_parent_child_probability,
    )
    from validate import validate_run  # type: ignore[no-redef]


@dataclass(frozen=True)
class Factor:
    variables: tuple[str, ...]
    table: dict[tuple[Genotype, ...], float]


class FactorLikelihoodError(Exception):
    """Raised when factor likelihood inputs cannot be evaluated."""


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise FactorLikelihoodError(f"{path} must contain a JSON object")
    return data


def make_unary_factor(variable: str, distribution: dict[Genotype, float]) -> Factor:
    return Factor(
        variables=(variable,),
        table={(genotype,): probability for genotype, probability in distribution.items()},
    )


def make_transmission_factor(
    father: str,
    mother: str,
    child: str,
    genotype_states: list[Genotype],
) -> Factor:
    """Build sparse P(child genotype | father genotype, mother genotype).

    The underlying Mendelian distributions are cached in genetics.py, so this
    remains exact while avoiding repeated recomputation across loci/models.
    """
    table: dict[tuple[Genotype, ...], float] = {}
    for father_gt in genotype_states:
        for mother_gt in genotype_states:
            child_dist = mendelian_child_distribution(father_gt, mother_gt)
            for child_gt, probability in child_dist.items():
                if probability > 0.0:
                    table[(father_gt, mother_gt, child_gt)] = probability
    return Factor(variables=(father, mother, child), table=table)


def make_single_parent_transmission_factor(
    parent: str,
    child: str,
    genotype_states: list[Genotype],
    allele_frequencies: dict[str, float],
) -> Factor:
    """Build sparse P(child genotype | known parent genotype)."""
    table: dict[tuple[Genotype, ...], float] = {}
    for parent_gt in genotype_states:
        for child_gt in genotype_states:
            probability = single_parent_child_probability(
                parent_genotype=parent_gt,
                child_genotype=child_gt,
                allele_frequencies=allele_frequencies,
            )
            if probability > 0.0:
                table[(parent_gt, child_gt)] = probability
    return Factor(variables=(parent, child), table=table)


def multiply_factors(left: Factor, right: Factor) -> Factor:
    """Multiply two sparse factors."""
    result_variables = tuple(dict.fromkeys(left.variables + right.variables))
    if not left.table or not right.table:
        return Factor(variables=result_variables, table={})

    left_positions = [result_variables.index(variable) for variable in left.variables]
    right_positions = [result_variables.index(variable) for variable in right.variables]
    common = [variable for variable in left.variables if variable in right.variables]
    left_common_positions = [left.variables.index(variable) for variable in common]
    right_common_positions = [right.variables.index(variable) for variable in common]

    right_index: dict[tuple[Genotype, ...], list[tuple[tuple[Genotype, ...], float]]] = {}
    for assignment, probability in right.table.items():
        key = tuple(assignment[position] for position in right_common_positions)
        right_index.setdefault(key, []).append((assignment, probability))

    result_table: dict[tuple[Genotype, ...], float] = {}
    for left_assignment, left_probability in left.table.items():
        key = tuple(left_assignment[position] for position in left_common_positions)
        for right_assignment, right_probability in right_index.get(key, []):
            result_assignment: list[Genotype | None] = [None] * len(result_variables)
            for position, value in zip(left_positions, left_assignment):
                result_assignment[position] = value
            for position, value in zip(right_positions, right_assignment):
                result_assignment[position] = value
            final_assignment = tuple(value for value in result_assignment if value is not None)
            result_table[final_assignment] = (
                result_table.get(final_assignment, 0.0)
                + left_probability * right_probability
            )

    return Factor(variables=result_variables, table=result_table)


def multiply_all(factors: Iterable[Factor]) -> Factor:
    iterator = iter(factors)
    try:
        product = next(iterator)
    except StopIteration:
        return Factor(variables=(), table={(): 1.0})
    for factor in iterator:
        product = multiply_factors(product, factor)
    return product


def sum_out(factor: Factor, variable: str) -> Factor:
    if variable not in factor.variables:
        return factor
    variable_index = factor.variables.index(variable)
    result_variables = tuple(v for v in factor.variables if v != variable)
    result_table: dict[tuple[Genotype, ...], float] = {}

    for assignment, probability in factor.table.items():
        reduced = tuple(
            value for index, value in enumerate(assignment) if index != variable_index
        )
        result_table[reduced] = result_table.get(reduced, 0.0) + probability

    return Factor(variables=result_variables, table=result_table)


def choose_elimination_variable(factors: list[Factor], variables: list[str]) -> str:
    """Choose a variable with the smallest current product scope."""
    best_variable = variables[0]
    best_score: tuple[int, int] | None = None
    for variable in variables:
        involved = [factor for factor in factors if variable in factor.variables]
        scope = set()
        row_sum = 0
        for factor in involved:
            scope.update(factor.variables)
            row_sum += len(factor.table)
        score = (len(scope), row_sum)
        if best_score is None or score < best_score:
            best_variable = variable
            best_score = score
    return best_variable


def variable_elimination_sparse(factors: list[Factor], variables: list[str]) -> tuple[float, dict[str, Any]]:
    """Eliminate all variables and return the scalar probability."""
    working = list(factors)
    remaining_variables = list(variables)
    max_factor_rows = max((len(factor.table) for factor in working), default=0)
    elimination_steps: list[dict[str, Any]] = []

    while remaining_variables:
        variable = choose_elimination_variable(working, remaining_variables)
        involved = [factor for factor in working if variable in factor.variables]
        uninvolved = [factor for factor in working if variable not in factor.variables]
        product = multiply_all(involved)
        summed = sum_out(product, variable)
        max_factor_rows = max(max_factor_rows, len(product.table), len(summed.table))
        elimination_steps.append(
            {
                "variable": variable,
                "involved_factors": len(involved),
                "product_variables": list(product.variables),
                "product_rows": len(product.table),
                "result_variables": list(summed.variables),
                "result_rows": len(summed.table),
            }
        )
        working = uninvolved
        if summed.table:
            working.append(summed)
        remaining_variables.remove(variable)

    final = multiply_all(working)
    if final.variables:
        raise FactorLikelihoodError("variable elimination ended with non-scalar factor")
    return final.table.get((), 0.0), {
        "initial_factors": len(factors),
        "max_factor_rows": max_factor_rows,
        "elimination_steps": elimination_steps,
    }


def variable_elimination(factors: list[Factor], variables: list[str]) -> tuple[float, dict[str, Any]]:
    """Exact contraction; NumPy is an optional acceleration of the same factors."""
    try:
        import numpy as np
    except ImportError:
        return variable_elimination_sparse(factors, variables)
    states = sorted({state for factor in factors for assignment in factor.table for state in assignment})
    state_index = {state: i for i, state in enumerate(states)}
    variable_index = {variable: i for i, variable in enumerate(variables)}
    working = []
    for factor in factors:
        table = np.zeros((len(states),) * len(factor.variables), dtype=float)
        for assignment, value in factor.table.items():
            table[tuple(state_index[state] for state in assignment)] = value
        working.append((factor.variables, table))
    remaining = list(variables)
    steps = []
    max_rows = max((len(factor.table) for factor in factors), default=0)
    while remaining:
        def score(variable):
            involved = [(scope, table) for scope, table in working if variable in scope]
            return (len({node for scope, _ in involved for node in scope}),
                    sum(table.size for _, table in involved))
        variable = min(remaining, key=score)
        involved = [(scope, table) for scope, table in working if variable in scope]
        working = [(scope, table) for scope, table in working if variable not in scope]
        union = tuple(dict.fromkeys(node for scope, _ in involved for node in scope))
        output_scope = tuple(node for node in union if node != variable)
        operands = []
        for scope, table in involved:
            operands.extend((table, [variable_index[node] for node in scope]))
        operands.append([variable_index[node] for node in output_scope])
        contracted = np.einsum(*operands, optimize=True) if involved else np.array(float(len(states)))
        rows = int(np.count_nonzero(contracted))
        max_rows = max(max_rows, rows)
        steps.append({'variable': variable, 'involved_factors': len(involved),
                      'result_variables': list(output_scope), 'result_rows': rows})
        working.append((output_scope, contracted))
        remaining.remove(variable)
    probability = math.prod(float(table) for _, table in working)
    return probability, {'initial_factors': len(factors), 'max_factor_rows': max_rows,
                         'elimination_steps': steps, 'backend': 'numpy_exact_contraction',
                         'factor_size_definition': 'Maximum stored initial or contracted nonzero rows; internal contraction buffers are not counted.'}


def build_locus_factors(
    *,
    person_ids: list[str],
    relationships: list[Relationship],
    observations: dict[str, Any],
    locus: str,
    allele_frequencies: dict[str, float],
    background_error: float,
    explicit_founder_priors: dict[str, Any] | None = None,
    error_model: dict[str, Any] | None = None,
    observation_model: dict[str, Any] | None = None,
) -> tuple[list[Factor], list[Genotype], list[str]]:
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

    factors: list[Factor] = []
    for founder in founders:
        founder_prior = founder_prior_for_person(
            founder,
            locus,
            allele_frequencies,
            explicit_founder_priors,
        )
        factors.append(make_unary_factor(founder, founder_prior))

    for rel in relationships:
        if "father" in rel and "mother" in rel:
            factors.append(
                make_transmission_factor(
                    rel["father"],
                    rel["mother"],
                    rel["child"],
                    genotype_states,
                )
            )
        else:
            factors.append(
                make_single_parent_transmission_factor(
                    rel["parent"],
                    rel["child"],
                    genotype_states,
                    allele_frequencies,
                )
            )

    for pid in person_ids:
        observation = observations.get(pid, {}).get(locus, {})
        observed_alleles = observation.get("alleles", [])
        observation_error_model = effective_observation_error_model(
            observation_model or error_model,
            observation,
        )
        if observation_model is None:
            raise FactorLikelihoodError("observation_model is required")
        obs_dist = {
            genotype: observation_likelihood_from_observation_model(
                observed_alleles,
                genotype,
                observation_error_model or observation_model,
                background_frequencies_for_locus(observation_model, locus, allele_frequencies),
            )
            for genotype in genotype_states
        }
        factors.append(make_unary_factor(pid, obs_dist))

    return factors, genotype_states, founders


def locus_likelihood_factor_graph(
    *,
    person_ids: list[str],
    relationships: list[Relationship],
    observations: dict[str, Any],
    locus: str,
    allele_frequencies: dict[str, float],
    background_error: float,
    explicit_founder_priors: dict[str, Any] | None = None,
    error_model: dict[str, Any] | None = None,
    observation_model: dict[str, Any] | None = None,
) -> tuple[float, dict[str, Any]]:
    factors, genotype_states, founders = build_locus_factors(
        person_ids=person_ids,
        relationships=relationships,
        observations=observations,
        locus=locus,
        allele_frequencies=allele_frequencies,
        explicit_founder_priors=explicit_founder_priors,
        background_error=background_error,
        error_model=error_model,
        observation_model=observation_model,
    )
    probability, ve_details = variable_elimination(factors, person_ids)
    alleles = collect_locus_alleles(
        locus,
        person_ids,
        observations,
        allele_frequencies,
        explicit_founder_priors,
    )
    return probability, {
        "locus": locus,
        "alleles": alleles,
        "observation_model": selected_observation_model_summary(
            error_model=error_model,
            observation_model=observation_model,
        ),
        "genotype_states": len(genotype_states),
        "persons": len(person_ids),
        "founders": founders,
        "probability": probability,
        "log_likelihood": math.log(probability) if probability > 0.0 else float("-inf"),
        **ve_details,
    }


def pedigree_log_likelihood_factor_graph(
    *,
    persons: dict[str, Any],
    observation_data: dict[str, Any],
    pedigree: dict[str, Any],
    founder_scenario: dict[str, Any],
    identity_hypothesis: dict[str, Any],
    background_error: float,
    error_model: dict[str, Any] | None = None,
    observation_model: dict[str, Any] | None = None,
) -> tuple[float, list[dict[str, Any]]]:
    person_ids = [person["id"] for person in persons["persons"]]
    relationships = pedigree.get("relationships", [])
    observations = observation_data.get("observations", {})
    try:
        person_ids, observations = apply_identity_hypothesis(
            person_ids=person_ids,
            relationships=relationships,
            observations=observations,
            identity_hypothesis=identity_hypothesis,
        )
    except Exception as exc:
        raise FactorLikelihoodError(str(exc)) from exc
    loci = observation_data.get("loci", [])
    allele_frequencies_by_locus = founder_scenario.get("allele_frequencies", {})
    explicit_founder_priors = founder_scenario.get("explicit_founder_priors", {})

    total_log_likelihood = 0.0
    locus_details: list[dict[str, Any]] = []
    for locus in loci:
        probability, details = locus_likelihood_factor_graph(
            person_ids=person_ids,
            relationships=relationships,
            observations=observations,
            locus=locus,
            allele_frequencies=allele_frequencies_by_locus[locus],
            explicit_founder_priors=explicit_founder_priors,
            background_error=background_error,
            error_model=error_model,
            observation_model=observation_model,
        )
        total_log_likelihood += math.log(probability) if probability > 0.0 else float("-inf")
        locus_details.append(details)
    return total_log_likelihood, locus_details


def evaluate_run_factor_graph(run_path: Path) -> list[dict[str, Any]]:
    validate_run(run_path)
    run = load_json(run_path)

    specs = model_specs_from_run(run)
    background_error = float(run.get("error_model", {}).get("background_error", 1e-3))
    error_model = run.get("error_model")
    observation_model = load_observation_model(run_path, run)
    observation_model_summary = selected_observation_model_summary(
        error_model=error_model,
        observation_model=observation_model,
    )

    results: list[dict[str, Any]] = []
    for spec in specs:
        persons = load_json(resolve_ref(run_path, spec.get("persons", run.get("persons"))))
        observations = load_json(resolve_ref(run_path, spec.get("observations", run.get("observations"))))
        pedigree = load_json(resolve_ref(run_path, spec["pedigree"]))
        founder_scenario = load_json(resolve_ref(run_path, spec["founder_scenario"]))
        identity_hypothesis = load_json(resolve_ref(run_path, spec["identity_hypothesis"]))
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
    return sorted(results, key=lambda row: row["log_likelihood"], reverse=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Exact full-pedigree likelihood by variable elimination"
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
        results = evaluate_run_factor_graph(run_path)
    except Exception as exc:
        print(f"FACTOR LIKELIHOOD FAILED: {exc}")
        return 1

    print("=== Exact Full-Pedigree Likelihoods: Variable Elimination ===")
    print(f"Run: {run_path}")
    print(f"{'Rank':>4}  {'logL':>14}  Pedigree / Founder / Identity")
    for rank, result in enumerate(results, start=1):
        print(
            f"{rank:>4}  {result['log_likelihood']:>14.6f}  "
            f"{result['pedigree']} / {result['founder_scenario']} / "
            f"{result['identity_hypothesis']}"
        )
        for locus in result["loci"]:
            print(
                f"      {locus['locus']}: p={locus['probability']:.8g}, "
                f"states={locus['genotype_states']}, "
                f"max_factor_rows={locus['max_factor_rows']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
