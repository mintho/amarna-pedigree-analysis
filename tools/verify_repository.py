#!/usr/bin/env python3
"""Check configuration closure, archived rankings and distinct founder settings."""

import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'Code/src'))
from summarize_primary_results import build_summary, structural_ranking_rows
from result_provenance import verify as verify_provenance
from repository_files import retained_files


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify():
    configs = sorted((ROOT / 'Data/models/runs').glob('*.json'))
    require(len(configs) == 7, 'Expected four primary and three focused configurations')
    referenced = set()
    for path in configs:
        run = json.loads(path.read_text())
        reference = run['observation_model'].get('background_frequency_file')
        require(isinstance(reference, str) and (path.parent / reference).resolve().is_file(),
                f'Missing fixed background: {path.name}')
        for model in run.get('models', [run]):
            for key in ('persons', 'observations', 'pedigree', 'founder_scenario', 'identity_hypothesis',
                        'pedigrees', 'founder_scenarios', 'identity_hypotheses'):
                refs = model.get(key, [])
                for ref in ([refs] if isinstance(refs, str) else refs):
                    target = (path.parent / ref).resolve()
                    require(target.is_relative_to(ROOT) and target.is_file(), f'Missing input: {target}')
                    referenced.add(target)
    for stem in ('focused_yuya_thuya_tiye_trio_l1_l2', 'focused_amenhotep_tiye_children_l1_l2',
                 'focused_yuya_kv55_grandparent_l1_l2'):
        signatures = set()
        for suffix in ('', '_broad_empirical_priors', '_restricted_endogamous_priors'):
            path = ROOT / 'Data/models/founder_scenarios' / f'{stem}_founders{suffix}.json'
            data = json.loads(path.read_text())
            referenced.add(path.resolve())
            signature = json.dumps({key: data.get(key, {}) for key in
                                    ('allele_frequencies', 'explicit_founder_priors')}, sort_keys=True)
            signatures.add(signature)
        require(len(signatures) == 3, f'Duplicate analytical founder setting: {stem}')
        matrix = json.loads((ROOT / 'Results/Focused_Robustness' / stem / 'matrix.json').read_text())
        require(len(matrix['matrix_points']) == 27, f'Incorrect focused grid: {stem}')
        require(all(row['results'][0]['rank'] == 1 for row in matrix['matrix_points']), 'Invalid focused ranking')
        jackknife = json.loads((ROOT / 'Results/Focused_Robustness' / stem / 'jackknife.json').read_text())
        require({row['omitted_locus'] for row in jackknife['leave_one_out']} == {'Locus1', 'Locus2'},
                f'Incorrect focused omissions: {stem}')
    for directory in ('pedigrees', 'founder_scenarios', 'identity_hypotheses'):
        supplied = {p.resolve() for p in (ROOT / 'Data/models' / directory).glob('*.json')}
        needed = {p for p in referenced if p.parent == ROOT / 'Data/models' / directory}
        require(supplied == needed, f'Unused or missing {directory}: {supplied ^ needed}')
    summary = build_summary(ROOT / 'Results/Primary_All_Locus')
    rows = structural_ranking_rows(summary)
    with (ROOT / 'Results/Primary_All_Locus/summary/structural_ranking_summary.csv').open() as handle:
        archived = list(csv.DictReader(handle))
    require(len(rows) == len(archived) == 20, 'Expected five structures per scenario')
    for actual, expected in zip(rows, archived):
        require(actual['structure'] == expected['structure'], 'Structure order differs')
        for key in ('log_likelihood', 'delta_log_likelihood', 'weight_percent'):
            require(math.isclose(actual[key], float(expected[key]), abs_tol=1e-8, rel_tol=1e-9),
                    f'Archived {key} differs')
    for scenario in summary['scenarios']:
        require(abs(scenario['tie_delta_log_likelihood']) < 1e-9, 'Historical readings no longer tied')
        require(len(scenario['generative_results']) == 6, 'Expected six labelled model outputs')
        require(math.isclose(sum(r['weight_percent'] for r in rows
                                 if r['scenario'] == scenario['scenario']), 100.0, abs_tol=1e-8),
                'Structural weights do not sum to 100%')
    for name, count in (('parameter_grid.csv', 36), ('leave_one_locus.csv', 32)):
        with (ROOT / 'Results/Primary_All_Locus_Sensitivity' / name).open() as handle:
            require(len(list(csv.DictReader(handle))) == count, f'Incorrect count: {name}')
    for doc in (ROOT / 'README.md', ROOT / 'LICENSING.md', ROOT / 'Code/README.md',
                ROOT / 'Data/README.md', ROOT / 'Data/docs/FREQUENCY_SOURCES.md', ROOT / 'Results/README.md'):
        for link in re.findall(r'\]\(([^)]+)\)', doc.read_text()):
            if link.startswith(('http:', 'https:', '#')):
                continue
            target = (doc.parent / unquote(link.split('#')[0])).resolve()
            require(target.is_relative_to(ROOT) and target.exists(), f'Broken documentation link: {link}')
    manifest = ROOT / 'MANIFEST.sha256'
    if manifest.exists():
        listed = set()
        for line in manifest.read_text().splitlines():
            digest, relative = line.split('  ', 1)
            listed.add(relative)
            require(hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest,
                    f'Integrity mismatch: {relative}')
        actual = {str(p.relative_to(ROOT)) for p in retained_files(ROOT) if p != manifest}
        require(listed == actual, 'Manifest inventory differs from the retained files')
        require(not any(Path(relative).suffix.lower() in {'.html', '.htm', '.pdf', '.png'}
                        for relative in listed), 'Publication assets or cached webpages in bundle')
    provenance = ROOT / 'Results/Validation/calculation_provenance.json'
    if provenance.exists():
        verify_provenance()
    print('PASS: seven configurations; input closure; three distinct focused priors each;')
    print('five-structure weights; 36/32 sensitivity rows; documentation links; portable file hashes.')


if __name__ == '__main__':
    verify()
