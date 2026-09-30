#!/usr/bin/env python3
"""Record or verify the exact code, input and numerical-result file hashes."""

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / 'Results/Validation/calculation_provenance.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory():
    inputs = sorted((ROOT / 'Data').rglob('*.json'))
    inputs = [p for p in inputs if 'docs' not in p.relative_to(ROOT).parts]
    code = sorted((ROOT / 'Code/src').glob('*.py'))
    results = []
    for directory in ('Primary_All_Locus', 'Primary_All_Locus_Sensitivity', 'Focused_Robustness'):
        results.extend(p for p in (ROOT / 'Results' / directory).rglob('*')
                       if p.suffix in ('.json', '.csv'))
    return {group: {str(p.relative_to(ROOT)): digest(p) for p in paths}
            for group, paths in (('code', code), ('inputs', inputs), ('results', sorted(results)))}


def verify():
    saved = json.loads(RECORD.read_text())
    actual = inventory()
    for group, files in actual.items():
        if saved[group] != files:
            changed = sorted(key for key in files.keys() | saved[group].keys()
                             if files.get(key) != saved[group].get(key))
            raise ValueError(f'Calculation provenance differs ({group}): {changed}. '
                             'Reproduce and review the calculations before updating the record.')
    print('PASS: calculation code, complete JSON inputs and numerical results match their provenance record.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-verified-calculation', action='store_true',
                        help='Explicitly record a newly reproduced, reviewed calculation; does not run it.')
    args = parser.parse_args()
    if args.record_verified_calculation:
        payload = {'schema_version': 1,
                   'description': 'Integrity record for the reviewed calculation; hashes are not a substitute for reproduction.',
                   **inventory()}
        RECORD.parent.mkdir(parents=True, exist_ok=True)
        RECORD.write_text(json.dumps(payload, indent=2) + '\n')
    verify()


if __name__ == '__main__':
    main()
