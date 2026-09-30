#!/usr/bin/env python3
"""Extract Egypt STR allele-frequency tables from the archived Duesseldorf database."""

from __future__ import annotations

import argparse
import json
import re
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Any


ARCHIVE_BASE = (
    "https://web.archive.org/web/20120920010159/"
    "http://www.uni-duesseldorf.de/WWW/MedFak/Serology/DNA-Systeme/"
)

LOCUS_PAGES = {
    "Locus1": {"marker": "D13S317", "page": "d13s317.htm"},
    "Locus2": {"marker": "D7S820", "page": "D7S820.html"},
    "Locus3": {
        "marker": "D2S1338",
        "page": "D2S1338.html",
        "population_selection_order": ["Israel (Jews)"],
        "population_selection_note": (
            "Paper page 6 states that Egyptian-specific D2S1338 frequencies "
            "were unavailable and modern Israeli frequencies from the same "
            "database were used."
        ),
    },
    "Locus4": {"marker": "D21S11", "page": "d21s11.html"},
    "Locus5": {"marker": "D16S539", "page": "D16S539.html"},
    "Locus6": {"marker": "D18S51", "page": "D18S51.html"},
    "Locus7": {"marker": "CSF1PO", "page": "csf.html"},
    "Locus8": {"marker": "FGA", "page": "fga.htm"},
}

PREFERRED_POPULATIONS = [
    "Egypt (pooled)",
    "Egypt",
    "Egypt (Cairo area)",
    "Egypt (Caucasians from the Cairo area)",
    "Egypt (Central, El-Minia)",
    "Egypt (El-Minia)",
    "Egypt (El-Minia City)",
]

SUMMARY_ROW_LABELS = {
    "s",
    "Σ",
    "σ",
    "sum",
    "h(exp)",
    "pe",
    "pi",
    "pd",
    "pic",
    "hwe",
}


@dataclass(frozen=True)
class PopulationRecord:
    population: str
    reference: str
    sample_size: int | None
    frequencies: dict[str, float]


class HtmlTableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.tables: list[list[list[str]]] = []
        self._table_stack: list[list[list[str]] | None] = []
        self._current_table: list[list[str]] | None = None
        self._current_row: list[str] | None = None
        self._current_cell: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        tag = tag.lower()
        if tag == "table":
            self._table_stack.append(self._current_table)
            self._current_table = []
        elif tag == "tr" and self._current_table is not None:
            self._current_row = []
        elif tag in {"td", "th"} and self._current_row is not None:
            self._current_cell = []

    def handle_endtag(self, tag: str) -> None:
        tag = tag.lower()
        if tag in {"td", "th"} and self._current_cell is not None:
            text = " ".join("".join(self._current_cell).split())
            self._current_row.append(text)
            self._current_cell = None
        elif tag == "tr" and self._current_row is not None:
            if any(cell for cell in self._current_row):
                assert self._current_table is not None
                self._current_table.append(self._current_row)
            self._current_row = None
        elif tag == "table" and self._current_table is not None:
            self.tables.append(self._current_table)
            self._current_table = self._table_stack.pop() if self._table_stack else None

    def handle_data(self, data: str) -> None:
        if self._current_cell is not None:
            self._current_cell.append(data)


def parse_float(value: str) -> float | None:
    try:
        return float(value.strip().replace(",", "."))
    except ValueError:
        return None


def parse_int(value: str) -> int | None:
    value = value.strip()
    if not value:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def normalize_allele_label(label: str) -> str:
    return label.strip().replace(" ", "")


def records_from_row_table(table: list[list[str]]) -> list[PopulationRecord]:
    header = table[0]
    if len(header) < 4 or header[:3] != ["Population", "Ref.", "n"]:
        return []

    allele_labels = [normalize_allele_label(label) for label in header[3:]]
    records = []
    for row in table[1:]:
        if len(row) < len(header) or not row[0] or row[0] in {"Europe", "Africa", "Asia"}:
            continue
        frequencies = {}
        for allele, value in zip(allele_labels, row[3:]):
            probability = parse_float(value)
            if probability is not None:
                frequencies[allele] = probability
        if frequencies:
            records.append(
                PopulationRecord(
                    population=row[0],
                    reference=row[1] if len(row) > 1 else "",
                    sample_size=parse_int(row[2]) if len(row) > 2 else None,
                    frequencies=frequencies,
                )
            )
    return records


def records_from_column_table(table: list[list[str]]) -> list[PopulationRecord]:
    population_row_index = next(
        (index for index, row in enumerate(table) if row and row[0] == "Population"),
        None,
    )
    if population_row_index is None:
        return []

    population_row = table[population_row_index]
    ref_row = table[population_row_index + 1] if len(table) > population_row_index + 1 else []
    n_row = table[population_row_index + 2] if len(table) > population_row_index + 2 else []
    allele_start = next(
        (
            index + 1
            for index, row in enumerate(table[population_row_index:], start=population_row_index)
            if row and row[0] == "Alleles"
        ),
        population_row_index + 3,
    )

    records = []
    for column_index, population in enumerate(population_row[1:], start=1):
        if not population:
            continue
        frequencies = {}
        for row in table[allele_start:]:
            if not row:
                continue
            label = normalize_allele_label(row[0])
            if not label or label.lower() in SUMMARY_ROW_LABELS:
                break
            if column_index >= len(row):
                continue
            probability = parse_float(row[column_index])
            if probability is not None:
                frequencies[label] = probability
        if frequencies:
            records.append(
                PopulationRecord(
                    population=population,
                    reference=ref_row[column_index] if column_index < len(ref_row) else "",
                    sample_size=parse_int(n_row[column_index]) if column_index < len(n_row) else None,
                    frequencies=frequencies,
                )
            )
    return records


def population_records_from_html(html_text: str) -> list[PopulationRecord]:
    parser = HtmlTableParser()
    parser.feed(html_text)
    records = []
    for table in parser.tables:
        records.extend(records_from_row_table(table))
        records.extend(records_from_column_table(table))
    return records


def select_egypt_record(records: list[PopulationRecord]) -> PopulationRecord | None:
    return select_population_record(records, PREFERRED_POPULATIONS, required_term="Egypt")


def select_population_record(
    records: list[PopulationRecord],
    preferred_populations: list[str],
    *,
    required_term: str | None = None,
) -> PopulationRecord | None:
    candidates = [
        record
        for record in records
        if required_term is None or required_term in record.population
    ]
    for preferred in preferred_populations:
        for record in candidates:
            if record.population == preferred:
                return record
    return candidates[0] if candidates else None


def fetch_text(url: str, cache_dir: Path | None = None) -> str:
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / re.sub(r"[^A-Za-z0-9_.-]+", "_", url)
        if cache_path.exists():
            return cache_path.read_text(encoding="latin-1", errors="replace")
    with urllib.request.urlopen(url, timeout=45) as response:
        raw = response.read()
    text = raw.decode("latin-1", errors="replace")
    if cache_dir is not None:
        cache_path.write_text(text, encoding="latin-1")
    return text


def allele_value_matches(
    frequencies: dict[str, float],
    value: float,
    *,
    tolerance: float = 5e-5,
) -> list[str]:
    return [
        allele
        for allele, probability in frequencies.items()
        if abs(probability - value) <= tolerance
    ]


def paper_prior_audit(
    *,
    paper_locus: dict[str, Any],
    record: PopulationRecord | None,
) -> list[dict[str, Any]]:
    rows = []
    for entry in paper_locus.get("alleles", []):
        length = str(entry["length"])
        raw = float(entry["prior_raw"])
        if record is None:
            rows.append(
                {
                    "color": entry["color"],
                    "paper_length": length,
                    "paper_prior_raw": raw,
                    "status": "no_egypt_population_found",
                }
            )
            continue
        source_at_length = record.frequencies.get(length)
        value_matches = allele_value_matches(record.frequencies, raw)
        if source_at_length is not None and abs(source_at_length - raw) <= 5e-5:
            status = "matches_paper_length"
        elif value_matches:
            status = "matches_different_duesseldorf_allele"
        else:
            status = "no_exact_value_match"
        rows.append(
            {
                "color": entry["color"],
                "paper_length": length,
                "paper_prior_raw": raw,
                "duesseldorf_frequency_at_paper_length": source_at_length,
                "duesseldorf_alleles_with_same_frequency": value_matches,
                "status": status,
            }
        )
    return rows


def build_duesseldorf_payload(
    *,
    paper_payload: dict[str, Any],
    cache_dir: Path | None = None,
) -> dict[str, Any]:
    loci = {}
    allele_frequencies = {}
    audit = {}

    for locus, source in LOCUS_PAGES.items():
        url = ARCHIVE_BASE + source["page"]
        html_text = fetch_text(url, cache_dir=cache_dir)
        records = population_records_from_html(html_text)
        selection_order = source.get("population_selection_order", PREFERRED_POPULATIONS)
        required_term = None if "population_selection_order" in source else "Egypt"
        selected = select_population_record(
            records,
            selection_order,
            required_term=required_term,
        )
        paper_locus = paper_payload["loci"][locus]
        audit_rows = paper_prior_audit(paper_locus=paper_locus, record=selected)
        audit[locus] = audit_rows

        loci[locus] = {
            "label": paper_locus.get("label"),
            "marker": source["marker"],
            "archive_url": url,
            "selected_population": selected.population if selected else None,
            "selected_reference": selected.reference if selected else None,
            "selected_sample_size": selected.sample_size if selected else None,
            "population_selection_order": selection_order,
            "population_selection_note": source.get(
                "population_selection_note",
                "Egyptian frequencies selected where available.",
            ),
            "available_egypt_populations": [
                {
                    "population": record.population,
                    "reference": record.reference,
                    "sample_size": record.sample_size,
                }
                for record in records
                if "Egypt" in record.population
            ],
            "duesseldorf_allele_frequencies": selected.frequencies if selected else {},
            "paper_color_mapping_audit": audit_rows,
        }
        allele_frequencies[locus] = selected.frequencies if selected else {}

    return {
        "schema_version": 1,
        "id": "duesseldorf_egypt_archive_2012",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "description": (
            "STR allele frequencies extracted from the archived "
            "Duesseldorf/Huckenbeck-Scheil database. Egyptian populations are "
            "used where available; Locus3/D2S1338 uses the Israeli fallback "
            "specified by the paper. Frequencies are stored by original STR "
            "allele labels, not by the diagram colour labels used in the Amarna "
            "scripts."
        ),
        "source": {
            "database_url": (
                "https://web.archive.org/web/20120920010159/"
                "http://www.uni-duesseldorf.de/WWW/MedFak/Serology/database.html"
            ),
            "authors": "Wolfgang Huckenbeck and Hans-Georg Scheil",
            "institution": "Universitaetsklinikum Duesseldorf",
            "archive_timestamp": "2012-09-20 01:01:59 UTC",
        },
        "parameters": {
            "population_selection_order": PREFERRED_POPULATIONS,
            "locus_specific_population_selection": {
                locus: source.get("population_selection_order", PREFERRED_POPULATIONS)
                for locus, source in LOCUS_PAGES.items()
            },
            "paper_match_tolerance": 5e-5,
        },
        "loci": loci,
        "allele_frequencies": allele_frequencies,
        "audit": audit,
    }


def markdown_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Duesseldorf Prior Audit",
        "",
        f"Source: `{payload['source']['database_url']}`",
        "",
        "## Selected Populations",
        "",
        "| Locus | Marker | Selected population | n | Audit status |",
        "|---|---|---|---:|---|",
    ]
    for locus, info in payload["loci"].items():
        statuses = [row["status"] for row in info["paper_color_mapping_audit"]]
        status_counts = ", ".join(
            f"{status}: {statuses.count(status)}" for status in sorted(set(statuses))
        )
        n = info["selected_sample_size"] if info["selected_sample_size"] is not None else ""
        lines.append(
            f"| {locus} | {info['marker']} | "
            f"{info['selected_population'] or 'n/a'} | {n} | {status_counts} |"
        )

    lines.extend(
        [
            "",
            "## Paper Color Mapping Audit",
            "",
            "| Locus | Color | Paper length | Paper prior | Duesseldorf at paper length | Same-frequency Duesseldorf alleles | Status |",
            "|---|---|---:|---:|---:|---|---|",
        ]
    )
    for locus, info in payload["loci"].items():
        for row in info["paper_color_mapping_audit"]:
            at_length = row.get("duesseldorf_frequency_at_paper_length")
            at_length_text = "" if at_length is None else f"{at_length:.6g}"
            matches = ", ".join(row.get("duesseldorf_alleles_with_same_frequency", []))
            lines.append(
                f"| {locus} | {row['color']} | {row['paper_length']} | "
                f"{row['paper_prior_raw']:.6g} | {at_length_text} | "
                f"{matches} | {row['status']} |"
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The extracted Duesseldorf frequencies use original STR allele labels. "
            "The Amarna scripts currently use diagram colour labels. Several Paper "
            "priors match Duesseldorf frequencies at a different allele label than "
            "the stored `paper_length`; therefore this file is an audited source "
            "dataset, not yet a drop-in replacement for the operational color-coded "
            "founder scenarios. For Locus3/D2S1338 the selected population is the "
            "paper-specified Israeli fallback, not Egypt.",
        ]
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
        description="Extract STR priors from the archived Duesseldorf database"
    )
    parser.add_argument(
        "--paper_priors",
        default="Data/data/allele_frequencies.reference.json",
    )
    parser.add_argument(
        "--output",
        default="reproduced/Frequency_prior_audit/allele_frequencies.duesseldorf_egypt.json",
    )
    parser.add_argument(
        "--report",
        default="reproduced/duesseldorf_egypt_prior_audit/audit.md",
    )
    parser.add_argument(
        "--cache_dir",
        default="Data/data/source_cache/duesseldorf_archive",
    )
    args = parser.parse_args()

    try:
        paper_payload = json.loads(Path(args.paper_priors).read_text(encoding="utf-8"))
        payload = build_duesseldorf_payload(
            paper_payload=paper_payload,
            cache_dir=Path(args.cache_dir) if args.cache_dir else None,
        )
        output = Path(args.output)
        report = Path(args.report)
        write_json(output, payload)
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(markdown_report(payload), encoding="utf-8")
    except Exception as exc:
        print(f"DUESSELDORF PRIOR EXTRACTION FAILED: {exc}")
        return 1

    print("=== Duesseldorf Egypt Prior Extraction ===")
    print(f"Output: {output}")
    print(f"Report: {report}")
    for locus, info in payload["loci"].items():
        print(
            f"{locus} {info['marker']}: "
            f"{info['selected_population'] or 'no Egypt population'}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
