#!/usr/bin/env python3
"""Build usable color-coded founder scenarios from Duesseldorf STR mappings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


DEFAULT_ALLELE_FLOOR = 1e-5


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def distribution_from_mapping(
    *,
    base_distribution: dict[str, float],
    locus_mapping: dict[str, Any],
    allele_floor: float,
) -> tuple[dict[str, float], dict[str, Any]]:
    if allele_floor < 0.0:
        raise ValueError("allele_floor must be >= 0")

    result: dict[str, float] = {}
    provenance: dict[str, Any] = {}
    for color in base_distribution:
        if color == "other":
            continue
        entry = locus_mapping.get(color)
        if not isinstance(entry, dict):
            raise ValueError(f"missing mapping for color {color!r}")
        if entry.get("requires_manual_review"):
            raise ValueError(f"mapping for color {color!r} still requires review")
        selected_frequency = entry.get("selected_frequency")
        if selected_frequency is None:
            raise ValueError(f"mapping for color {color!r} has no selected_frequency")
        raw_frequency = float(selected_frequency)
        operational_frequency = max(raw_frequency, allele_floor)
        result[color] = operational_frequency
        provenance[color] = {
            "selected_str_allele": entry.get("selected_str_allele"),
            "source_frequency": raw_frequency,
            "operational_frequency": operational_frequency,
            "selection_status": entry.get("selection_status"),
            "manual_override": entry.get("manual_override"),
        }

    used = sum(result.values())
    if used > 1.0 + 1e-9:
        raise ValueError(f"mapped allele frequencies exceed 1: {used}")
    result["other"] = max(0.0, 1.0 - used)
    return result, provenance


def build_duesseldorf_founder_scenario(
    *,
    base_founder_scenario: dict[str, Any],
    color_mapping: dict[str, Any],
    scenario_id: str | None = None,
    allele_floor: float | None = None,
    preserve_explicit_founder_priors: bool = False,
) -> dict[str, Any]:
    base_id = base_founder_scenario["id"]
    parameters = dict(base_founder_scenario.get("allele_frequency_parameters", {}))
    effective_floor = (
        DEFAULT_ALLELE_FLOOR
        if allele_floor is None
        else float(allele_floor)
    )
    if allele_floor is None and "allele_floor_background_error0" in parameters:
        effective_floor = float(parameters["allele_floor_background_error0"])

    allele_frequencies: dict[str, dict[str, float]] = {}
    provenance: dict[str, Any] = {}
    for locus, base_distribution in base_founder_scenario.get("allele_frequencies", {}).items():
        locus_mapping = (
            color_mapping.get("loci", {})
            .get(locus, {})
            .get("color_to_str_allele", {})
        )
        if not locus_mapping:
            raise ValueError(f"missing color mapping for {locus}")
        allele_frequencies[locus], provenance[locus] = distribution_from_mapping(
            base_distribution={
                str(color): float(probability)
                for color, probability in base_distribution.items()
            },
            locus_mapping=locus_mapping,
            allele_floor=effective_floor,
        )

    return {
        "schema_version": 1,
        "id": scenario_id or f"{base_id}_duesseldorf_color_priors",
        "description": (
            f"Color-coded Duesseldorf founder prior scenario derived from {base_id}. "
            "Colors are preserved as model alleles; frequencies are taken from the "
            "verified locus-specific Duesseldorf STR mapping."
        ),
        "variant_type": "duesseldorf_color_priors",
        "base_founder_scenario": base_id,
        "allele_frequency_source": "data/allele_frequencies.duesseldorf_egypt.json",
        "color_mapping_source": "data/color_to_str_allele.duesseldorf_egypt.json",
        "allele_frequency_parameters": {
            **parameters,
            "variant_parameters": {
                "allele_floor_background_error0": effective_floor,
                "other_policy": (
                    "For each locus, mapped color frequencies are copied from the "
                    "selected Duesseldorf STR allele frequencies; raw zero "
                    "frequencies are raised to allele_floor_background_error0; remaining "
                    "mass is stored as color allele 'other'."
                ),
                "explicit_founder_priors_policy": (
                    "preserved"
                    if preserve_explicit_founder_priors
                    else "removed_population_prior_variant"
                ),
            },
        },
        "allele_frequencies": allele_frequencies,
        "color_to_str_allele_provenance": provenance,
        "explicit_founder_priors": (
            base_founder_scenario.get("explicit_founder_priors", {})
            if preserve_explicit_founder_priors
            else {}
        ),
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def default_output_path(base_path: Path) -> Path:
    return Path("reproduced/Founder_Priors") / f"{base_path.stem}_reference_priors.json"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a usable color-coded Duesseldorf founder scenario"
    )
    parser.add_argument("base_founder_scenario")
    parser.add_argument(
        "--color_mapping",
        default="Data/data/color_to_str_allele.duesseldorf_egypt.json",
    )
    parser.add_argument("--output", default=None)
    parser.add_argument("--id", default=None)
    parser.add_argument("--allele_floor", type=float, default=None)
    parser.add_argument(
        "--preserve_explicit_founder_priors",
        action="store_true",
        help="Preserve explicit person-specific founder priors from the base scenario.",
    )
    args = parser.parse_args()

    base_path = Path(args.base_founder_scenario)
    output = Path(args.output) if args.output else default_output_path(base_path)
    try:
        payload = build_duesseldorf_founder_scenario(
            base_founder_scenario=load_json(base_path),
            color_mapping=load_json(Path(args.color_mapping)),
            scenario_id=args.id,
            allele_floor=args.allele_floor,
            preserve_explicit_founder_priors=args.preserve_explicit_founder_priors,
        )
        write_json(output, payload)
    except Exception as exc:
        print(f"DUESSELDORF FOUNDER SCENARIO BUILD FAILED: {exc}")
        return 1

    print("=== Duesseldorf Founder Scenario ===")
    print(f"Base: {base_path}")
    print(f"Output: {output}")
    print(f"ID: {payload['id']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
