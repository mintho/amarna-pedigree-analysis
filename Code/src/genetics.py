#!/usr/bin/env python3
"""Core genetic primitives for STR pedigree likelihoods."""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from math import isfinite
from typing import Mapping

Allele = str
Genotype = tuple[Allele, Allele]
Distribution = dict[Genotype, float]


class GeneticsError(ValueError):
    """Raised when genetic inputs are malformed."""


def canonical_genotype(allele_a: Allele, allele_b: Allele) -> Genotype:
    """Return an unordered genotype represented as a sorted 2-tuple."""
    if not isinstance(allele_a, str) or not allele_a:
        raise GeneticsError("allele_a must be a non-empty string")
    if not isinstance(allele_b, str) or not allele_b:
        raise GeneticsError("allele_b must be a non-empty string")
    return tuple(sorted((allele_a, allele_b)))  # type: ignore[return-value]


def genotype_key(genotype: Genotype) -> str:
    """Serialize a genotype as a stable JSON key."""
    gt = canonical_genotype(genotype[0], genotype[1])
    return f"{gt[0]}/{gt[1]}"


def parse_genotype_key(key: str) -> Genotype:
    """Parse a genotype key of the form 'allele_a/allele_b'."""
    if not isinstance(key, str) or "/" not in key:
        raise GeneticsError(f"invalid genotype key: {key!r}")
    parts = key.split("/")
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise GeneticsError(f"invalid genotype key: {key!r}")
    return canonical_genotype(parts[0], parts[1])


def normalize_distribution(dist: Mapping[Genotype, float]) -> Distribution:
    """Validate and normalize a genotype distribution."""
    if not dist:
        raise GeneticsError("distribution must not be empty")

    clean: Distribution = {}
    total = 0.0
    for genotype, probability in dist.items():
        if len(genotype) != 2:
            raise GeneticsError(f"invalid genotype key: {genotype!r}")
        gt = canonical_genotype(genotype[0], genotype[1])
        if not isinstance(probability, (int, float)) or not isfinite(probability):
            raise GeneticsError(f"invalid probability for {genotype!r}: {probability!r}")
        if probability < 0.0:
            raise GeneticsError(f"negative probability for {genotype!r}")
        clean[gt] = clean.get(gt, 0.0) + float(probability)
        total += float(probability)

    if total <= 0.0:
        raise GeneticsError("distribution total must be positive")

    return {genotype: probability / total for genotype, probability in clean.items()}


def genotype_prior_from_mapping(dist: Mapping[str, float]) -> Distribution:
    """Build a genotype distribution from JSON genotype-key probabilities."""
    parsed: dict[Genotype, float] = {}
    for key, probability in dist.items():
        genotype = parse_genotype_key(key)
        parsed[genotype] = parsed.get(genotype, 0.0) + probability
    return normalize_distribution(parsed)


def unordered_genotypes(alleles: list[Allele] | tuple[Allele, ...]) -> list[Genotype]:
    """Return all unordered diploid genotypes over the given alleles."""
    unique = sorted(set(alleles))
    if not unique:
        raise GeneticsError("alleles must not be empty")
    for allele in unique:
        if not isinstance(allele, str) or not allele:
            raise GeneticsError("alleles must contain non-empty strings")

    genotypes: list[Genotype] = []
    for i, allele_a in enumerate(unique):
        for allele_b in unique[i:]:
            genotypes.append(canonical_genotype(allele_a, allele_b))
    return genotypes


def founder_prior_from_allele_frequencies(
    allele_frequencies: Mapping[Allele, float],
) -> Distribution:
    """Build an unordered genotype prior under Hardy-Weinberg proportions."""
    if not allele_frequencies:
        raise GeneticsError("allele_frequencies must not be empty")

    total = 0.0
    freqs: dict[Allele, float] = {}
    for allele, probability in allele_frequencies.items():
        if not isinstance(allele, str) or not allele:
            raise GeneticsError("allele labels must be non-empty strings")
        if not isinstance(probability, (int, float)) or not isfinite(probability):
            raise GeneticsError(f"invalid allele frequency for {allele!r}")
        if probability < 0.0:
            raise GeneticsError(f"negative allele frequency for {allele!r}")
        freqs[allele] = float(probability)
        total += float(probability)

    if total <= 0.0:
        raise GeneticsError("allele frequencies must have positive total")

    freqs = {allele: probability / total for allele, probability in freqs.items()}
    alleles = sorted(freqs)
    prior: Distribution = {}

    for i, allele_a in enumerate(alleles):
        pa = freqs[allele_a]
        prior[canonical_genotype(allele_a, allele_a)] = pa * pa
        for allele_b in alleles[i + 1 :]:
            pb = freqs[allele_b]
            prior[canonical_genotype(allele_a, allele_b)] = 2.0 * pa * pb

    return normalize_distribution(prior)


def gamete_distribution(parent_genotype: Genotype) -> dict[Allele, float]:
    """Return allele transmission probabilities for one parent genotype."""
    genotype = canonical_genotype(parent_genotype[0], parent_genotype[1])
    return _gamete_distribution_cached(genotype)


@lru_cache(maxsize=None)
def _gamete_distribution_cached(genotype: Genotype) -> dict[Allele, float]:
    allele_a, allele_b = genotype
    if allele_a == allele_b:
        return {allele_a: 1.0}
    return {allele_a: 0.5, allele_b: 0.5}


def mendelian_child_distribution(
    father_genotype: Genotype,
    mother_genotype: Genotype,
) -> Distribution:
    """Return P(child genotype | father genotype, mother genotype)."""
    father = canonical_genotype(father_genotype[0], father_genotype[1])
    mother = canonical_genotype(mother_genotype[0], mother_genotype[1])
    return _mendelian_child_distribution_cached(father, mother)


@lru_cache(maxsize=None)
def _mendelian_child_distribution_cached(
    father_genotype: Genotype,
    mother_genotype: Genotype,
) -> Distribution:
    father_gametes = gamete_distribution(father_genotype)
    mother_gametes = gamete_distribution(mother_genotype)
    child_dist: defaultdict[Genotype, float] = defaultdict(float)

    for paternal, paternal_prob in father_gametes.items():
        for maternal, maternal_prob in mother_gametes.items():
            child_dist[canonical_genotype(paternal, maternal)] += (
                paternal_prob * maternal_prob
            )

    return normalize_distribution(child_dist)


def clear_genetics_caches() -> None:
    """Clear internal caches used by genetic primitive functions."""
    _gamete_distribution_cached.cache_clear()
    _mendelian_child_distribution_cached.cache_clear()


def genetics_cache_info() -> dict[str, object]:
    """Return cache statistics for diagnostics/tests."""
    return {
        "gamete_distribution": _gamete_distribution_cached.cache_info(),
        "mendelian_child_distribution": _mendelian_child_distribution_cached.cache_info(),
    }


def validate_rate(value: float, name: str) -> float:
    if not isinstance(value, (int, float)) or not isfinite(value):
        raise GeneticsError(f"{name} must be finite")
    rate = float(value)
    if rate < 0.0 or rate >= 1.0:
        raise GeneticsError(f"{name} must satisfy 0 <= {name} < 1")
    return rate


def dropout_observation_likelihood(
    observed_alleles: list[Allele] | tuple[Allele, ...],
    true_genotype: Genotype,
    dropout: float,
) -> float:
    """Return P(observation | true genotype) under independent allelic dropout.

    Empty observations remain uninformative because current data files use an
    empty allele list to mean "no usable observation", not a confirmed blank
    electropherogram.
    """
    dropout = validate_rate(dropout, "dropout")
    alleles = list(observed_alleles)
    if len(alleles) > 2:
        raise GeneticsError("observed_alleles must have length 0..2")
    for allele in alleles:
        if not isinstance(allele, str) or not allele:
            raise GeneticsError("observed alleles must be non-empty strings")
    if not alleles:
        return 1.0

    genotype = canonical_genotype(true_genotype[0], true_genotype[1])
    observed = (
        canonical_genotype(alleles[0], alleles[0])
        if len(alleles) == 1
        else canonical_genotype(alleles[0], alleles[1])
    )

    probability = 0.0
    copy_alleles = list(genotype)
    for keep_first in (False, True):
        for keep_second in (False, True):
            kept = []
            pattern_probability = 1.0
            for keep, allele in zip((keep_first, keep_second), copy_alleles):
                if keep:
                    kept.append(allele)
                    pattern_probability *= 1.0 - dropout
                else:
                    pattern_probability *= dropout
            if len(kept) != len(alleles):
                continue
            if not kept:
                continue
            kept_observation = (
                canonical_genotype(kept[0], kept[0])
                if len(kept) == 1
                else canonical_genotype(kept[0], kept[1])
            )
            if kept_observation == observed:
                probability += pattern_probability
    return probability


def authentic_scored_observation_likelihood(
    observed_alleles: list[Allele] | tuple[Allele, ...],
    true_genotype: Genotype,
    dropout: float,
) -> float:
    """Return normalized P(scored observation | true genotype, dropout).

    Empty observations are not scored calls in the current data model and
    should be handled by the caller as uninformative missing data.
    """
    dropout = validate_rate(dropout, "dropout")
    alleles = list(observed_alleles)
    if not alleles:
        raise GeneticsError("scored observation must contain one or two alleles")
    raw = dropout_observation_likelihood(alleles, true_genotype, dropout)
    denominator = 1.0 - dropout * dropout
    if denominator <= 0.0:
        raise GeneticsError("dropout is too close to 1 for scored-call normalization")
    return raw / denominator


def background_call_distribution_probability(
    observed_alleles: list[Allele] | tuple[Allele, ...],
    allele_frequencies: Mapping[Allele, float],
    dropout: float,
) -> float:
    """Return model-independent background P(scored observation).

    The background component draws a genotype from locus-level allele
    frequencies and then applies the same normalized scored-call dropout model.
    """
    if not observed_alleles:
        raise GeneticsError("background call probability requires a scored call")
    genotype_prior = founder_prior_from_allele_frequencies(allele_frequencies)
    probability = 0.0
    for genotype, genotype_probability in genotype_prior.items():
        probability += genotype_probability * authentic_scored_observation_likelihood(
            observed_alleles,
            genotype,
            dropout,
        )
    return probability


def generative_str_observation_likelihood(
    observed_alleles: list[Allele] | tuple[Allele, ...],
    true_genotype: Genotype,
    allele_frequencies: Mapping[Allele, float],
    dropout: float,
    background_error: float,
) -> float:
    """Return normalized STR observation probability.

    Empty allele lists mean "no usable call" in the current data files and
    contribute a constant factor of 1.0.
    """
    alleles = list(observed_alleles)
    if len(alleles) > 2:
        raise GeneticsError("observed_alleles must have length 0..2")
    for allele in alleles:
        if not isinstance(allele, str) or not allele:
            raise GeneticsError("observed alleles must be non-empty strings")
    if not alleles:
        return 1.0

    dropout = validate_rate(dropout, "dropout")
    background_error = validate_rate(background_error, "background_error")
    authentic = authentic_scored_observation_likelihood(
        alleles,
        true_genotype,
        dropout,
    )
    background = background_call_distribution_probability(
        alleles,
        allele_frequencies,
        dropout,
    )
    return (1.0 - background_error) * authentic + background_error * background


def combined_background_error(error_model: Mapping[str, float]) -> float:
    """Return a single background-error mixture rate from model fields."""
    if "background_error" in error_model:
        return validate_rate(float(error_model["background_error"]), "background_error")
    miscall = validate_rate(float(error_model.get("miscall", 0.0)), "miscall")
    contamination = validate_rate(
        float(error_model.get("contamination", 0.0)),
        "contamination",
    )
    return 1.0 - (1.0 - miscall) * (1.0 - contamination)


def observation_likelihood_from_observation_model(
    observed_alleles: list[Allele] | tuple[Allele, ...],
    true_genotype: Genotype,
    observation_model: Mapping[str, float | str],
    allele_frequencies: Mapping[Allele, float],
) -> float:
    """Return observation probability under an explicit observation model."""
    model_type = observation_model.get("type")
    if model_type != "generative_str":
        raise GeneticsError(f"unsupported observation model type: {model_type!r}")

    dropout = validate_rate(float(observation_model.get("dropout", 0.0)), "dropout")
    background_error = combined_background_error(observation_model)  # type: ignore[arg-type]
    return generative_str_observation_likelihood(
        observed_alleles,
        true_genotype,
        allele_frequencies,
        dropout,
        background_error,
    )
