#!/usr/bin/env python3
"""Reproduce the documented configurations without overwriting retained results."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'Code/src'
OUT = ROOT / 'reproduced'
FOCUSED = (
    'focused_yuya_thuya_tiye_trio_l1_l2',
    'focused_amenhotep_tiye_children_l1_l2',
    'focused_yuya_kv55_grandparent_l1_l2',
)


def run(script, *args):
    command = [sys.executable, str(SRC / script), *map(str, args)]
    print('+ ' + ' '.join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def verify_csv(actual, expected):
    with actual.open(newline='') as handle:
        a = list(csv.DictReader(handle))
    with expected.open(newline='') as handle:
        b = list(csv.DictReader(handle))
    if len(a) != len(b):
        raise ValueError(f'Row count differs: {actual}')
    for i, (left, right) in enumerate(zip(a, b), 1):
        if left.keys() != right.keys():
            raise ValueError(f'Columns differ: {actual}')
        for key in left:
            try:
                matches = math.isclose(float(left[key]), float(right[key]), rel_tol=1e-9, abs_tol=1e-8)
            except ValueError:
                matches = left[key] == right[key]
            if not matches:
                raise ValueError(f'{actual}, row {i}, {key}: {left[key]} != {right[key]}')
    print(f'PASS: {actual.name}, {len(a)} rows match archive', flush=True)


def check():
    env = dict(os.environ, PYTHONPATH=str(ROOT / 'Code'))
    tested = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', 'Code/tests', '-v'],
                            cwd=ROOT, env=env, text=True, capture_output=True)
    validation = OUT / 'validation'
    validation.mkdir(parents=True, exist_ok=True)
    (validation / 'unit_tests.txt').write_text(tested.stdout + tested.stderr)
    print(tested.stdout + tested.stderr, flush=True)
    tested.check_returncode()
    subprocess.run([sys.executable, str(ROOT / 'tools/verify_repository.py')], cwd=ROOT, check=True)
    for config in sorted((ROOT / 'Data/models/runs').glob('*.json')):
        run('validate.py', config)
    run('audit_primary_evidence.py', *(ROOT / 'Data/models/runs' /
        f'scenario_{s}_primary_models_all_loci_generative_str.json' for s in 'abcd'),
        '--json-out', OUT / 'validation/primary_evidence_audit.json',
        '--md-out', OUT / 'validation/PRIMARY_EVIDENCE_AUDIT.md')


def primary():
    output = OUT / 'Primary_All_Locus'
    for s in 'abcd':
        run('compare.py', ROOT / 'Data/models/runs' /
            f'scenario_{s}_primary_models_all_loci_generative_str.json',
            '--engine', 'factors', '--out_dir', output / f'scenario_{s}')
    summary(output)


def summary(input_dir):
    output = OUT / 'Primary_All_Locus/summary'
    run('summarize_primary_results.py', '--input-dir', input_dir, '--out-dir', output)
    if not REFRESH:
        for name in ('structural_ranking_summary.csv', 'ranking_summary.csv', 'per_locus_vs_nearest.csv'):
            verify_csv(output / name, ROOT / 'Results/Primary_All_Locus/summary' / name)


def focused():
    output = OUT / 'Focused_Robustness'
    matrices, jackknives = [], []
    for stem in FOCUSED:
        config = ROOT / 'Data/models/runs' / f'{stem}_run.json'
        founders = [ROOT / 'Data/models/founder_scenarios' / f'{stem}_founders{s}.json'
                    for s in ('', '_broad_empirical_priors', '_restricted_endogamous_priors')]
        signatures = {json.dumps({key: data.get(key, {}) for key in
                                  ('allele_frequencies', 'explicit_founder_priors')}, sort_keys=True)
                      for data in (json.loads(path.read_text()) for path in founders)}
        if len(signatures) != 3:
            raise ValueError(f'{stem}: expected three distinct founder priors, found {len(signatures)}; '
                             'resolve input provenance before reproducing the published grid')
        destination = output / stem
        run('sensitivity_matrix.py', config, '--founder_scenarios', *founders,
            '--background_errors', '0.0001,0.001,0.01', '--dropouts', '0,0.1,0.2',
            '--engine', 'factors', '--out_dir', destination)
        run('locus_jackknife.py', config, '--engine', 'factors', '--out_dir', destination)
        matrices.append(destination / 'matrix.json')
        jackknives.append(destination / 'jackknife.json')
    run('stability_report.py', *matrices, '--out_dir', output)
    run('locus_jackknife_report.py', *jackknives, '--out_dir', output)
    run('robustness_conclusions.py', output / 'stability_report.json', '--out_dir', output)
    if REFRESH:
        return
    # Compare numerical Markdown tables, excluding path-dependent headings/provenance.
    for new, old in (('stability_report.md', 'stability_report.md'),
                     ('jackknife_report.md', 'jackknife_report.md'),
                     ('robustness_conclusions.md', 'robustness_conclusions.md')):
        table = lambda p: [line for line in p.read_text().splitlines() if line.startswith('|')]
        if table(output / new) != table(ROOT / 'Results/Focused_Robustness' / old):
            raise ValueError(f'Focused tables differ: {new}')
        print(f'PASS: {new} tables match archive', flush=True)


def sensitivity():
    output = OUT / 'Primary_All_Locus_Sensitivity'
    run('primary_generative_sensitivity.py', '--engine', 'factors', '--out-dir', output)
    if not REFRESH:
        for name in ('parameter_grid.csv', 'leave_one_locus.csv'):
            verify_csv(output / name, ROOT / 'Results/Primary_All_Locus_Sensitivity' / name)


def priors():
    output = OUT / 'Frequency_prior_audit'
    run('extract_duesseldorf_egypt_priors.py',
        '--paper_priors', ROOT / 'Data/data/allele_frequencies.reference.json',
        '--output', output / 'allele_frequencies.duesseldorf_egypt.json',
        '--report', output / 'audit.md',
        '--cache_dir', OUT / 'source_cache/duesseldorf_archive')
    run('build_color_str_mapping.py',
        '--duesseldorf_priors', output / 'allele_frequencies.duesseldorf_egypt.json',
        '--manual_overrides', ROOT / 'Data/data/color_to_str_allele.manual_overrides.json',
        '--output', output / 'color_to_str_allele.duesseldorf_egypt.json',
        '--report', output / 'color_mapping.md')


def main():
    global REFRESH
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task', choices=('check', 'summary', 'primary', 'focused', 'sensitivity', 'priors', 'all'))
    parser.add_argument('--refresh', action='store_true', help='Generate new outputs without comparing saved results; never overwrites Results.')
    args = parser.parse_args()
    REFRESH = args.refresh
    os.chdir(ROOT)
    if args.task == 'summary':
        summary(ROOT / 'Results/Primary_All_Locus')
    elif args.task == 'all':
        for task in (check, priors, primary, focused, sensitivity):
            task()
        subprocess.run([sys.executable, str(ROOT / 'tools/validate_engines.py')], cwd=ROOT, check=True)
    else:
        globals()[args.task]()
    print(f'Completed {args.task}. Generated outputs: {OUT}', flush=True)


if __name__ == '__main__':
    main()
