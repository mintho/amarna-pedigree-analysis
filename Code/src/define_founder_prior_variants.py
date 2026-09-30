#!/usr/bin/env python3
"""Generate named founder-prior variants from a baseline scenario."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_BROAD_UNIFORM_BLEND = 0.5
DEFAULT_ENDOGAMOUS_OTHER_CAP = 0.01
DEFAULT_MIN_ALLELE_FREQUENCY = 1e-6


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def normalize(dist: dict[str, float]) -> dict[str, float]:
    total = sum(dist.values())
    if total <= 0.0:
        raise ValueError("cannot normalize empty or zero distribution")
    return {key: value / total for key, value in sorted(dist.items()) if value > 0.0}


def observed_alleles_by_locus(observations: dict[str, Any] | None) -> dict[str, set[str]]:
    if not observations:
        return {}
    result: dict[str, set[str]] = {}
    for by_locus in observations.get("observations", {}).values():
        if not isinstance(by_locus, dict):
            continue
        for locus, observation in by_locus.items():
            if not isinstance(observation, dict):
                continue
            alleles = observation.get("alleles", [])
            if not isinstance(alleles, list):
                continue
            result.setdefault(locus, set()).update(
                str(allele) for allele in alleles if isinstance(allele, str) and allele
            )
    return result


def smooth_towards_uniform(
    dist: dict[str, float],
    *,
    extra_alleles: set[str] | None = None,
    blend: float = DEFAULT_BROAD_UNIFORM_BLEND,
) -> dict[str, float]:
    if blend < 0.0 or blend > 1.0:
        raise ValueError("blend must be in [0, 1]")
    normalized = normalize({key: float(value) for key, value in dist.items()})
    keys = set(normalized)
    keys.update(extra_alleles or set())
    if not keys:
        raise ValueError("cannot smooth empty distribution")
    uniform = 1.0 / len(keys)
    return normalize(
        {
            key: (1.0 - blend) * normalized.get(key, 0.0) + blend * uniform
            for key in keys
        }
    )


def restrict_to_observed_pool(
    dist: dict[str, float],
    *,
    observed_alleles: set[str] | None = None,
    other_cap: float = DEFAULT_ENDOGAMOUS_OTHER_CAP,
    min_allele_frequency: float = DEFAULT_MIN_ALLELE_FREQUENCY,
) -> dict[str, float]:
    if other_cap < 0.0 or other_cap >= 1.0:
        raise ValueError("other_cap must satisfy 0 <= other_cap < 1")
    if min_allele_frequency < 0.0:
        raise ValueError("min_allele_frequency must be >= 0")

    normalized = normalize({key: float(value) for key, value in dist.items()})
    observed = set(observed_alleles or set())
    allele_keys = {key for key in normalized if key != "other"}
    retained = observed & allele_keys
    if not retained:
        retained = allele_keys
    if not retained:
        retained = observed
    if not retained:
        raise ValueError("restricted pool has no alleles")

    other_mass = min(normalized.get("other", 0.0), other_cap)
    retained_mass = 1.0 - other_mass
    retained_dist = normalize(
        {
            allele: max(normalized.get(allele, 0.0), min_allele_frequency)
            for allele in retained
        }
    )
    result = {
        allele: probability * retained_mass
        for allele, probability in retained_dist.items()
    }
    if other_mass > 0.0:
        result["other"] = other_mass
    return normalize(result)


def variant_payload(
    base: dict[str, Any],
    *,
    variant_id: str,
    variant_type: str,
    description: str,
    allele_frequencies: dict[str, dict[str, float]],
    explicit_founder_priors: dict[str, Any],
    parameters: dict[str, Any],
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "id": variant_id,
        "description": description,
        "variant_type": variant_type,
        "base_founder_scenario": base["id"],
        "allele_frequency_source": base.get("allele_frequency_source", base["id"]),
        "allele_frequency_parameters": {
            **dict(base.get("allele_frequency_parameters", {})),
            "variant_parameters": parameters,
        },
        "allele_frequencies": allele_frequencies,
        "explicit_founder_priors": explicit_founder_priors,
    }


def build_founder_prior_variants(
    base: dict[str, Any],
    *,
    observed_alleles: dict[str, set[str]] | None = None,
    id_prefix: str | None = None,
    broad_uniform_blend: float = DEFAULT_BROAD_UNIFORM_BLEND,
    endogamous_other_cap: float = DEFAULT_ENDOGAMOUS_OTHER_CAP,
    min_allele_frequency: float = DEFAULT_MIN_ALLELE_FREQUENCY,
) -> dict[str, dict[str, Any]]:
    base_id = str(id_prefix or base["id"])
    allele_frequencies = base.get("allele_frequencies", {})
    if not isinstance(allele_frequencies, dict) or not allele_frequencies:
        raise ValueError("base founder scenario must contain allele_frequencies")
    observed = observed_alleles or {}

    broad = {
        locus: smooth_towards_uniform(
            {str(allele): float(probability) for allele, probability in dist.items()},
            extra_alleles=observed.get(locus, set()),
            blend=broad_uniform_blend,
        )
        for locus, dist in allele_frequencies.items()
    }
    endogamous = {
        locus: restrict_to_observed_pool(
            {str(allele): float(probability) for allele, probability in dist.items()},
            observed_alleles=observed.get(locus, set()),
            other_cap=endogamous_other_cap,
            min_allele_frequency=min_allele_frequency,
        )
        for locus, dist in allele_frequencies.items()
    }

    return {
        "broad_empirical_priors": variant_payload(
            base,
            variant_id=f"{base_id}_broad_empirical_priors",
            variant_type="broad_empirical_priors",
            description=(
                f"Broad empirical founder priors derived from {base['id']} by "
                "shrinking each locus toward a uniform allele distribution. "
                "Explicit person-specific founder genotype priors are removed."
            ),
            allele_frequencies=broad,
            explicit_founder_priors={},
            parameters={
                "uniform_blend": broad_uniform_blend,
                "explicit_founder_priors_policy": "removed",
            },
        ),
        "restricted_endogamous_priors": variant_payload(
            base,
            variant_id=f"{base_id}_restricted_endogamous_priors",
            variant_type="restricted_endogamous_priors",
            description=(
                f"Restricted endogamous founder priors derived from {base['id']} "
                "by concentrating each locus on the observed/listed allele pool "
                "and capping 'other' mass. Explicit person-specific founder "
                "genotype priors are removed."
            ),
            allele_frequencies=endogamous,
            explicit_founder_priors={},
            parameters={
                "other_cap": endogamous_other_cap,
                "min_allele_frequency": min_allele_frequency,
                "explicit_founder_priors_policy": "removed",
            },
        ),
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def write_variant_files(
    *,
    variants: dict[str, dict[str, Any]],
    output_dir: Path,
) -> dict[str, str]:
    written = {}
    for variant_name, payload in variants.items():
        path = output_dir / f"{payload['id']}.json"
        write_json(path, payload)
        written[variant_name] = str(path)
    return written


def manifest_payload(
    *,
    base_path: Path,
    observations_path: Path | None,
    variants: dict[str, dict[str, Any]],
    written_files: dict[str, str],
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "base_founder_scenario": str(base_path),
        "observations": str(observations_path) if observations_path else None,
        "variant_order": list(variants),
        "variants": {
            name: {
                "id": payload["id"],
                "variant_type": payload["variant_type"],
                "path": written_files[name],
                "description": payload["description"],
                "parameters": payload.get("allele_frequency_parameters", {}).get(
                    "variant_parameters", {}
                ),
            }
            for name, payload in variants.items()
        },
    }


def default_manifest_path(base_path: Path, output_dir: Path) -> Path:
    return output_dir / f"{base_path.stem}_variant_set.json"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate founder-prior variant JSON files"
    )
    parser.add_argument("base_founder_scenario")
    parser.add_argument(
        "--observations",
        default=None,
        help="Optional observation JSON used to define the restricted allele pool.",
    )
    parser.add_argument(
        "--output_dir",
        default="reproduced/Founder_Priors",
    )
    parser.add_argument(
        "--id_prefix",
        default=None,
        help="Optional id prefix. Defaults to the base founder scenario id.",
    )
    parser.add_argument("--broad_uniform_blend", type=float, default=DEFAULT_BROAD_UNIFORM_BLEND)
    parser.add_argument("--endogamous_other_cap", type=float, default=DEFAULT_ENDOGAMOUS_OTHER_CAP)
    parser.add_argument("--min_allele_frequency", type=float, default=DEFAULT_MIN_ALLELE_FREQUENCY)
    parser.add_argument("--manifest", default=None)
    args = parser.parse_args()

    base_path = Path(args.base_founder_scenario)
    observations_path = Path(args.observations) if args.observations else None
    output_dir = Path(args.output_dir)
    manifest_path = (
        Path(args.manifest)
        if args.manifest
        else default_manifest_path(base_path, output_dir)
    )

    try:
        observations = load_json(observations_path) if observations_path else None
        variants = build_founder_prior_variants(
            load_json(base_path),
            observed_alleles=observed_alleles_by_locus(observations),
            id_prefix=args.id_prefix,
            broad_uniform_blend=args.broad_uniform_blend,
            endogamous_other_cap=args.endogamous_other_cap,
            min_allele_frequency=args.min_allele_frequency,
        )
        written_files = write_variant_files(variants=variants, output_dir=output_dir)
        write_json(
            manifest_path,
            manifest_payload(
                base_path=base_path,
                observations_path=observations_path,
                variants=variants,
                written_files=written_files,
            ),
        )
    except Exception as exc:
        print(f"FOUNDER PRIOR VARIANT GENERATION FAILED: {exc}")
        return 1

    print("=== Founder-Prior Variants ===")
    print(f"Base: {base_path}")
    print(f"Manifest: {manifest_path}")
    for name, path in written_files.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
