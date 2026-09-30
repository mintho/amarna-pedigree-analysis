#!/usr/bin/env python3
"""Audit primary evidence-set consistency before likelihood reruns."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

try:
    from .validate import validate_run
except ImportError:
    from validate import validate_run  # type: ignore[no-redef]


class EvidenceAuditError(Exception):
    """Raised when a primary evidence audit fails."""


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise EvidenceAuditError(f"{path} must contain a JSON object")
    return data


def resolve_ref(base_file: Path, ref: str) -> Path:
    return (base_file.parent / ref).resolve()


def profile_tuple(observation: dict[str, Any], loci: list[str]) -> tuple[tuple[str, tuple[str, ...]], ...]:
    return tuple(
        (
            locus,
            tuple(sorted((observation.get(locus, {}) or {}).get("alleles", []) or [])),
        )
        for locus in loci
    )


def has_scored_call(profile: tuple[tuple[str, tuple[str, ...]], ...]) -> bool:
    return any(alleles for _, alleles in profile)


def profile_hash(profile: tuple[tuple[str, tuple[str, ...]], ...]) -> str:
    payload = json.dumps(profile, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def relationship_nodes(pedigree: dict[str, Any]) -> set[str]:
    nodes: set[str] = set()
    for relationship in pedigree.get("relationships", []):
        nodes.add(relationship["child"])
        if "father" in relationship and "mother" in relationship:
            nodes.add(relationship["father"])
            nodes.add(relationship["mother"])
        else:
            nodes.add(relationship["parent"])
    return nodes


def fetal_mothers(pedigree: dict[str, Any]) -> dict[str, str]:
    mothers: dict[str, str] = {}
    for relationship in pedigree.get("relationships", []):
        child = relationship.get("child")
        if child in {"Fetus1", "Fetus2"}:
            mother = relationship.get("mother")
            if not isinstance(mother, str) or not mother:
                raise EvidenceAuditError(f"{pedigree.get('id')}: {child} has no mother")
            mothers[child] = mother
    missing = {"Fetus1", "Fetus2"} - set(mothers)
    if missing:
        raise EvidenceAuditError(
            f"{pedigree.get('id')}: missing fetal-mother links for {sorted(missing)}"
        )
    return mothers


def model_evidence_summary(run_path: Path, model: dict[str, str]) -> dict[str, Any]:
    persons = load_json(resolve_ref(run_path, model["persons"]))
    observations = load_json(resolve_ref(run_path, model["observations"]))
    pedigree = load_json(resolve_ref(run_path, model["pedigree"]))

    person_ids = {person["id"] for person in persons.get("persons", [])}
    nodes = relationship_nodes(pedigree)
    loci = observations.get("loci", [])
    observed_profiles = observations.get("observations", {})

    scored_profiles: dict[str, Any] = {}
    for person_id, by_locus in observed_profiles.items():
        profile = profile_tuple(by_locus, loci)
        if has_scored_call(profile):
            scored_profiles[person_id] = {
                "profile_hash": profile_hash(profile),
                "profile": profile,
            }

    mothers = fetal_mothers(pedigree)
    mother_ids = sorted(set(mothers.values()))
    if len(mother_ids) != 1:
        raise EvidenceAuditError(
            f"{model.get('label')}: Fetus1 and Fetus2 have different mothers: {mothers}"
        )
    fetal_mother = mother_ids[0]
    fetal_mother_profile = scored_profiles.get(fetal_mother)
    if fetal_mother_profile is None:
        raise EvidenceAuditError(
            f"{model.get('label')}: fetal mother {fetal_mother!r} has no scored STR profile"
        )

    kv21b_mentions = sorted(
        person_id
        for person_id in set(person_ids) | set(observed_profiles) | nodes
        if person_id == "KV21B"
    )

    return {
        "label": model.get("label") or model["pedigree"],
        "pedigree": pedigree.get("id"),
        "observations": model["observations"],
        "loci": loci,
        "scored_profile_count": len(scored_profiles),
        "scored_labels": sorted(scored_profiles),
        "scored_profile_hashes": sorted(
            item["profile_hash"] for item in scored_profiles.values()
        ),
        "fetal_mother": fetal_mother,
        "fetal_mother_profile_hash": fetal_mother_profile["profile_hash"],
        "kv21b_mentions": kv21b_mentions,
    }


def audit_run(run_config: Path) -> dict[str, Any]:
    run_path = run_config.resolve()
    validate_run(run_path)
    run = load_json(run_path)
    observation_model = run.get("observation_model")
    if not isinstance(observation_model, dict):
        raise EvidenceAuditError(f"{run_path}: missing explicit observation_model")
    if observation_model.get("type") != "generative_str":
        raise EvidenceAuditError(
            f"{run_path}: observation_model.type must be generative_str"
        )

    model_summaries = [
        model_evidence_summary(run_path, model)
        for model in run.get("models", [])
    ]
    if not model_summaries:
        raise EvidenceAuditError(f"{run_path}: no explicit models to audit")

    reference_hashes = model_summaries[0]["scored_profile_hashes"]
    reference_mother_hash = model_summaries[0]["fetal_mother_profile_hash"]
    failures: list[str] = []
    for summary in model_summaries:
        if summary["scored_profile_hashes"] != reference_hashes:
            failures.append(f"{summary['label']}: scored profile set differs")
        if summary["kv21b_mentions"]:
            failures.append(
                f"{summary['label']}: KV21B appears in {summary['kv21b_mentions']}"
            )
        if summary["fetal_mother_profile_hash"] != reference_mother_hash:
            failures.append(
                f"{summary['label']}: fetal mother profile differs from KV21A reference"
            )

    return {
        "run_config": str(run_path),
        "primaryd": run.get("id"),
        "observation_model": observation_model,
        "status": "fail" if failures else "pass",
        "failures": failures,
        "scored_profile_hashes": reference_hashes,
        "models": model_summaries,
    }


def audit_runs(run_configs: list[Path]) -> dict[str, Any]:
    audits = [audit_run(path) for path in run_configs]
    return {
        "schema_version": 1,
        "status": "fail" if any(audit["status"] != "pass" for audit in audits) else "pass",
        "runs": audits,
    }


def markdown_report(payload: dict[str, Any]) -> str:
    lines = [
        "# primary Evidence-Set Audit",
        "",
        f"Overall status: **{payload['status'].upper()}**",
        "",
        "This audit checks the primary generative STR run files. It verifies "
        "that each run uses the explicit `generative_str` observation model, that "
        "all model families within a scenario use the same scored STR profile set, "
        "that KV21B is not present in the scored primary evidence set, and that the "
        "maternal profile attached to the KV62 fetuses is the same KV21A profile "
        "across model families even when historical labels differ.",
        "",
    ]
    for run in payload["runs"]:
        lines.extend(
            [
                f"## {run['primaryd']}",
                "",
                f"Status: **{run['status'].upper()}**",
                "",
                "Observation model:",
                "",
                f"```json\n{json.dumps(run['observation_model'], indent=2)}\n```",
                "",
            ]
        )
        if run["failures"]:
            lines.extend(["Failures:", ""])
            lines.extend(f"- {failure}" for failure in run["failures"])
            lines.append("")
        lines.extend(
            [
                "| Model | Scored labels | Fetal mother | KV21B present? |",
                "|---|---|---|---|",
            ]
        )
        for model in run["models"]:
            labels = ", ".join(model["scored_labels"])
            kv21b = "yes" if model["kv21b_mentions"] else "no"
            lines.append(
                f"| {model['label']} | {labels} | "
                f"{model['fetal_mother']} "
                f"(`{model['fetal_mother_profile_hash']}`) | {kv21b} |"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit primary evidence sets")
    parser.add_argument("run_config", nargs="+", help="primary config JSON files")
    parser.add_argument("--json-out", default=None, help="Optional JSON output path")
    parser.add_argument("--md-out", default=None, help="Optional Markdown output path")
    args = parser.parse_args()

    try:
        payload = audit_runs([Path(path) for path in args.run_config])
    except Exception as exc:
        print(f"AUDIT FAILED: {exc}")
        return 1

    if args.json_out:
        Path(args.json_out).write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    if args.md_out:
        Path(args.md_out).write_text(markdown_report(payload), encoding="utf-8")

    print(f"primary evidence audit: {payload['status'].upper()}")
    for run in payload["runs"]:
        print(f"- {run['primaryd']}: {run['status'].upper()}")
    return 0 if payload["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
