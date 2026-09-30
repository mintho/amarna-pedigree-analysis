#!/usr/bin/env python3
"""Build locus-specific color-to-STR-allele mapping from the Duesseldorf audit."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def label_contains_length(label: str, length: str) -> bool:
    parts = label.replace("<", "").replace(">", "").split("/")
    return any(part == length for part in parts)


def frequency_for_allele(
    frequencies: dict[str, float],
    allele_label: str | None,
) -> float | None:
    if allele_label is None:
        return None
    return frequencies.get(allele_label)


def override_for_color(
    manual_overrides: dict[str, Any] | None,
    locus: str,
    color: str,
) -> dict[str, Any] | None:
    if not manual_overrides:
        return None
    overrides = manual_overrides.get("overrides", manual_overrides)
    by_locus = overrides.get(locus, {})
    if not isinstance(by_locus, dict):
        return None
    override = by_locus.get(color)
    return override if isinstance(override, dict) else None


def select_mapping(row: dict[str, Any]) -> dict[str, Any]:
    paper_length = str(row["paper_length"])
    paper_prior = float(row["paper_prior_raw"])
    candidates = list(row.get("duesseldorf_alleles_with_same_frequency", []) or [])
    source_at_length = row.get("duesseldorf_frequency_at_paper_length")

    if row["status"] == "no_egypt_population_found":
        return {
            "selected_str_allele": None,
            "selection_status": "unavailable_no_egypt_population",
            "confidence": "unavailable",
            "requires_manual_review": True,
        }

    exact_or_grouped = [
        candidate for candidate in candidates if label_contains_length(candidate, paper_length)
    ]
    if source_at_length is not None and abs(float(source_at_length) - paper_prior) <= 5e-5:
        return {
            "selected_str_allele": paper_length,
            "selection_status": "paper_length_frequency_match",
            "confidence": "high",
            "requires_manual_review": False,
        }

    if len(exact_or_grouped) == 1:
        return {
            "selected_str_allele": exact_or_grouped[0],
            "selection_status": "paper_length_matches_grouped_allele_label",
            "confidence": "medium" if paper_prior == 0.0 else "high",
            "requires_manual_review": paper_prior == 0.0 and len(candidates) > 1,
        }

    if len(candidates) == 1:
        return {
            "selected_str_allele": candidates[0],
            "selection_status": "unique_frequency_match",
            "confidence": "high" if paper_prior > 0.0 else "medium",
            "requires_manual_review": False,
        }

    if candidates:
        return {
            "selected_str_allele": None,
            "selection_status": (
                "ambiguous_zero_frequency_match"
                if paper_prior == 0.0
                else "ambiguous_equal_frequency_match"
            ),
            "confidence": "low",
            "requires_manual_review": True,
        }

    return {
        "selected_str_allele": None,
        "selection_status": "no_frequency_match",
        "confidence": "low",
        "requires_manual_review": True,
    }


def build_color_str_mapping(
    duesseldorf_payload: dict[str, Any],
    *,
    manual_overrides: dict[str, Any] | None = None,
) -> dict[str, Any]:
    loci = {}
    review_items = []
    for locus, info in duesseldorf_payload["loci"].items():
        frequencies = {
            allele: float(probability)
            for allele, probability in info.get("duesseldorf_allele_frequencies", {}).items()
        }
        color_map = {}
        for row in info.get("paper_color_mapping_audit", []):
            selected = select_mapping(row)
            manual_override = override_for_color(
                manual_overrides,
                locus,
                row["color"],
            )
            if manual_override:
                selected = {
                    "selected_str_allele": manual_override["selected_str_allele"],
                    "selection_status": "manual_color_legend_override",
                    "confidence": "manual",
                    "requires_manual_review": False,
                }
            selected_frequency = frequency_for_allele(
                frequencies,
                selected["selected_str_allele"],
            )
            entry = {
                "color": row["color"],
                "paper_length_original": row["paper_length"],
                "paper_prior_raw": row["paper_prior_raw"],
                "duesseldorf_frequency_at_paper_length": row.get(
                    "duesseldorf_frequency_at_paper_length"
                ),
                "candidate_str_alleles": row.get(
                    "duesseldorf_alleles_with_same_frequency", []
                ),
                "selected_str_allele": selected["selected_str_allele"],
                "selected_frequency": selected_frequency,
                "selection_status": selected["selection_status"],
                "confidence": selected["confidence"],
                "requires_manual_review": selected["requires_manual_review"],
                "manual_override": manual_override or None,
            }
            color_map[row["color"]] = entry
            if entry["requires_manual_review"]:
                review_items.append(
                    {
                        "locus": locus,
                        "marker": info["marker"],
                        "color": row["color"],
                        "paper_length_original": row["paper_length"],
                        "paper_prior_raw": row["paper_prior_raw"],
                        "candidate_str_alleles": entry["candidate_str_alleles"],
                        "selection_status": entry["selection_status"],
                    }
                )

        loci[locus] = {
            "marker": info["marker"],
            "selected_population": info.get("selected_population"),
            "selected_sample_size": info.get("selected_sample_size"),
            "archive_url": info.get("archive_url"),
            "color_to_str_allele": color_map,
        }

    return {
        "schema_version": 1,
        "id": "duesseldorf_egypt_color_to_str_allele_mapping",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "description": (
            "Locus-specific mapping from Amarna diagram colour labels to "
            "Duesseldorf STR allele labels. Colors are retained as operational "
            "and visualization labels."
        ),
        "source": {
            "duesseldorf_egypt_priors": "data/allele_frequencies.duesseldorf_egypt.json",
            "audit_report": "reproduced/duesseldorf_egypt_prior_audit/audit.md",
        },
        "policy": {
            "color_scope": "locus-specific",
            "colors_are_preserved": True,
            "manual_review_required_before_drop_in_prior_replacement": bool(review_items),
            "manual_overrides_applied": manual_overrides is not None,
        },
        "loci": loci,
        "manual_review_items": review_items,
    }


def markdown_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Color To STR Allele Mapping",
        "",
        "Colors are retained as locus-specific operational labels. The mapped STR "
        "allele is the Duesseldorf source-table allele whose frequency explains "
        "the Paper prior value, where this can be inferred automatically.",
        "",
        "## Summary",
        "",
        "| Locus | Marker | Population | Colors | Manual review |",
        "|---|---|---|---:|---:|",
    ]
    for locus, info in payload["loci"].items():
        entries = list(info["color_to_str_allele"].values())
        review_count = sum(1 for entry in entries if entry["requires_manual_review"])
        lines.append(
            f"| {locus} | {info['marker']} | {info['selected_population'] or 'n/a'} | "
            f"{len(entries)} | {review_count} |"
        )

    lines.extend(
        [
            "",
            "## Mapping",
            "",
            "| Locus | Color | Paper length | Paper prior | Selected STR allele | Source frequency | Status | Review |",
            "|---|---|---:|---:|---|---:|---|---:|",
        ]
    )
    for locus, info in payload["loci"].items():
        for entry in info["color_to_str_allele"].values():
            frequency = entry["selected_frequency"]
            frequency_text = "" if frequency is None else f"{frequency:.6g}"
            selected = entry["selected_str_allele"] or ""
            review = "yes" if entry["requires_manual_review"] else "no"
            lines.append(
                f"| {locus} | {entry['color']} | {entry['paper_length_original']} | "
                f"{entry['paper_prior_raw']:.6g} | {selected} | {frequency_text} | "
                f"{entry['selection_status']} | {review} |"
            )

    if payload["manual_review_items"]:
        lines.extend(
            [
                "",
                "## Manual Review Items",
                "",
                "| Locus | Marker | Color | Paper length | Paper prior | Candidates | Reason |",
                "|---|---|---|---:|---:|---|---|",
            ]
        )
        for item in payload["manual_review_items"]:
            lines.append(
                f"| {item['locus']} | {item['marker']} | {item['color']} | "
                f"{item['paper_length_original']} | {item['paper_prior_raw']:.6g} | "
                f"{', '.join(item['candidate_str_alleles'])} | "
                f"{item['selection_status']} |"
            )

    return "\n".join(lines).rstrip() + "\n"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build locus-specific color-to-STR-allele mapping"
    )
    parser.add_argument(
        "--duesseldorf_priors",
        default="Data/data/allele_frequencies.duesseldorf_egypt.json",
    )
    parser.add_argument(
        "--output",
        default="reproduced/Frequency_prior_audit/color_to_str_allele.duesseldorf_egypt.json",
    )
    parser.add_argument(
        "--report",
        default="reproduced/duesseldorf_egypt_prior_audit/color_mapping.md",
    )
    parser.add_argument(
        "--manual_overrides",
        default="Data/data/color_to_str_allele.manual_overrides.json",
        help="Optional manual override JSON. Used if the file exists.",
    )
    args = parser.parse_args()

    try:
        override_path = Path(args.manual_overrides) if args.manual_overrides else None
        manual_overrides = (
            load_json(override_path)
            if override_path is not None and override_path.exists()
            else None
        )
        payload = build_color_str_mapping(
            load_json(Path(args.duesseldorf_priors)),
            manual_overrides=manual_overrides,
        )
        write_json(Path(args.output), payload)
        report = Path(args.report)
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(markdown_report(payload), encoding="utf-8")
    except Exception as exc:
        print(f"COLOR MAPPING BUILD FAILED: {exc}")
        return 1

    print("=== Color To STR Allele Mapping ===")
    print(f"Output: {args.output}")
    print(f"Report: {args.report}")
    print(f"Manual review items: {len(payload['manual_review_items'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
