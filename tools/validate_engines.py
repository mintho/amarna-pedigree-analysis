#!/usr/bin/env python3
"""Compare exact implementations on current inputs and store validation results."""

import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'Code'))
from src import likelihood_factors
from src.likelihood_enumeration import model_specs_from_run
from src.locus_jackknife import evaluate_model_for_loci
from src.compare import evaluate_with_engine


def main():
    checks = []
    run_path = ROOT / 'Data/models/runs/scenario_c_primary_models_all_loci_generative_str.json'
    run = json.loads(run_path.read_text())
    accelerated = likelihood_factors.variable_elimination
    for spec in model_specs_from_run(run):
        arguments = dict(run_path=run_path, run=run, spec=spec, loci=['Locus2'], engine='factors')
        fast = evaluate_model_for_loci(**arguments)['log_likelihood']
        try:
            likelihood_factors.variable_elimination = likelihood_factors.variable_elimination_sparse
            reference = evaluate_model_for_loci(**arguments)['log_likelihood']
        finally:
            likelihood_factors.variable_elimination = accelerated
        difference = abs(fast - reference)
        assert difference < 1e-9, (spec, difference)
        checks.append({'comparison': 'full_pedigree_L2_sparse_vs_accelerated', 'pedigree': spec['pedigree'],
                       'accelerated_log_likelihood': fast, 'sparse_log_likelihood': reference,
                       'absolute_difference': difference})
    path = ROOT / 'Data/models/runs/focused_yuya_thuya_tiye_trio_l1_l2_run.json'
    factors = evaluate_with_engine(path, 'factors')
    enumeration = {row['pedigree']: row for row in evaluate_with_engine(path, 'enumeration')}
    for row in factors:
        reference = enumeration[row['pedigree']]['log_likelihood']
        difference = abs(row['log_likelihood'] - reference)
        assert difference < 1e-9
        checks.append({'comparison': 'focused_trio_enumeration_vs_factors', 'pedigree': row['pedigree'],
                       'factor_log_likelihood': row['log_likelihood'], 'enumeration_log_likelihood': reference,
                       'absolute_difference': difference})
    target = ROOT / 'reproduced/validation/engine_agreement.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps({'checks': checks, 'tolerance': 1e-9}, indent=2) + '\n')
    print(f'PASS: {len(checks)} exact-engine comparisons; maximum difference {max(c["absolute_difference"] for c in checks):.3g}')


if __name__ == '__main__':
    main()
