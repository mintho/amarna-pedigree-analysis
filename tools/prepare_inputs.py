#!/usr/bin/env python3
"""Synchronise operational founder frequencies with the verified reference key."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'Data'


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def prepare():
    mapping = json.loads((DATA / 'data/color_to_str_allele.duesseldorf_egypt.json').read_text())
    reference_path = DATA / 'data/allele_frequencies.reference.json'
    reference = json.loads(reference_path.read_text())
    reported_key = json.loads((DATA / 'data/color_to_reported_str_allele.json').read_text())['loci']
    floor = float(reference['parameters']['allele_floor_epsilon0'])
    frequencies = {}
    for locus, row in reference['loci'].items():
        values = {}
        for allele in row['alleles']:
            entry = mapping['loci'][locus]['color_to_str_allele'][allele['color']]
            if entry.get('requires_manual_review'):
                raise ValueError(f'Unverified mapping: {locus} {allele["color"]}')
            allele['length'] = reported_key.get(locus, {}).get(allele['color'], entry['selected_str_allele'])
            allele['source_table_allele'] = entry['selected_str_allele']
            allele['prior_raw'] = entry['selected_frequency']
            allele['prior_operational'] = max(float(entry['selected_frequency']), floor)
            values[allele['color']] = allele['prior_operational']
        if sum(values.values()) > 1:
            raise ValueError(f'Frequencies exceed one: {locus}')
        values['other'] = 1 - sum(values.values())
        row['other_prior_operational'] = values['other']
        frequencies[locus] = values
    reference['allele_frequencies'] = frequencies
    reference['id'] = 'reference_allele_frequencies'
    reference['description'] = 'Operational reference frequencies with verified locus-specific repeat-count coding and the documented zero-frequency floor.'
    reference['source'] = {'frequency_data': 'allele_frequencies.duesseldorf_egypt.json',
                           'verified_key': 'color_to_str_allele.duesseldorf_egypt.json',
                           'archive_cache': 'source_cache/duesseldorf_archive'}
    write(reference_path, reference)
    write(DATA / 'data/background_allele_frequencies.json', {
        'schema_version': 1, 'id': 'fixed_reference_background',
        'description': 'Fixed reference distribution for background calls in all comparisons. Founder sensitivity variants never modify this distribution.',
        'allele_frequencies': frequencies,
    })
    for config_path in sorted((DATA / 'models/runs').glob('*.json')):
        config = json.loads(config_path.read_text())
        config['observation_model']['background_source'] = 'fixed_reference_frequencies'
        config['observation_model']['background_frequency_file'] = '../../data/background_allele_frequencies.json'
        for spec in config.get('models', []):
            if spec['label'].startswith(('Hawass_', 'Belmonte_')):
                spec['structure_id'] = 'shared_kv55_kv35yl_tutankhamun_structure'
            path = (config_path.parent / spec['founder_scenario']).resolve()
            founder = json.loads(path.read_text())
            founder['allele_frequencies'] = frequencies
            founder['allele_frequency_source'] = 'data/allele_frequencies.reference.json'
            write(path, founder)
        write(config_path, config)
    print('Prepared eight reference distributions, fixed background and 24 primary founder files.')


if __name__ == '__main__':
    prepare()
