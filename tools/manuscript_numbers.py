#!/usr/bin/env python3
"""Generate bilingual manuscript table rows and numerical macros from Results."""

import csv
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'Results'
OUT = RESULTS / 'Manuscript_Tables'
FOCUSED = ('focused_yuya_thuya_tiye_trio_l1_l2',
           'focused_amenhotep_tiye_children_l1_l2',
           'focused_yuya_kv55_grandparent_l1_l2')


def load(path):
    return json.loads(path.read_text())


def csv_rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def fmt(value, precision, lang):
    text = f'{float(value):.{precision}f}'
    return text.replace('.', ',') if lang == 'de' else text


def integer(value, lang):
    return f'{round(value):,}'.replace(',', r'\,') if lang == 'de' else f'{round(value):,}'


def write_rows(lang, name, rows):
    target = OUT / lang / name
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text('\n'.join(' & '.join(row) + r' \\' for row in rows) + '\n')


def main():
    summary = load(RESULTS / 'Primary_All_Locus/summary/summary.json')
    ranking = csv_rows(RESULTS / 'Primary_All_Locus/summary/structural_ranking_summary.csv')
    loci = csv_rows(RESULTS / 'Primary_All_Locus/summary/per_locus_vs_nearest.csv')
    jack = csv_rows(RESULTS / 'Primary_All_Locus_Sensitivity/leave_one_locus.csv')
    weakest = min(jack, key=lambda row: float(row['nearest_non_tied_bayes_factor']))
    facts = {'primary': [], 'focused': [], 'weakest_primary_omission': weakest}
    for scenario in summary['scenarios']:
        best = next(row for row in ranking if row['scenario'] == scenario['scenario'] and row['rank'] == '1')
        facts['primary'].append({'scenario': scenario['scenario'], 'log_likelihood': float(best['log_likelihood']),
                                 'weight_percent': float(best['weight_percent']),
                                 'gap': -scenario['nearest_non_tied_delta_log_likelihood'],
                                 'ratio': scenario['nearest_non_tied_bayes_factor']})
    for stem in FOCUSED:
        matrix = load(RESULTS / 'Focused_Robustness' / stem / 'matrix.json')
        omissions = load(RESULTS / 'Focused_Robustness' / stem / 'jackknife.json')['leave_one_out']
        baseline = matrix['baseline'][0]
        weights = [next(row['weight_percent'] for row in point['results']
                        if row['pedigree'] == baseline['pedigree']) for point in matrix['matrix_points']]
        omit = {row['omitted_locus']: next(result['weight_percent'] for result in row['results']
                if result['pedigree'] == baseline['pedigree']) for row in omissions}
        assert len(weights) == 27 and all(point['results'][0]['pedigree'] == baseline['pedigree']
                                         for point in matrix['matrix_points'])
        facts['focused'].append({'excerpt': stem, 'baseline_weight': baseline['weight_percent'],
                                 'grid_min': min(weights), 'grid_max': max(weights), 'omissions': omit})
    for lang in ('en', 'de'):
        primary_rows = []
        for row in facts['primary']:
            shared = r'Shared Hawass-derived/\newline Belmonte structure' if lang == 'en' else r'Gemeinsame Hawass-derived/\newline Belmonte-Struktur'
            logl = fmt(row['log_likelihood'], 6, lang)
            gap = fmt(-row['gap'], 6, lang)
            if lang == 'de':
                logl = '$' + logl.replace(',', '{,}') + '$'
                gap = '$' + gap.replace(',', '{,}') + '$'
            ratio = f'{row["ratio"]:,.2f}' if lang == 'en' else f'{row["ratio"]:,.2f}'.replace(',', '@').replace('.', ',').replace('@', r'\,')
            tag = 'LR' if row['scenario'] in 'ab' else 'BF'
            primary_rows.append([row['scenario'], shared, logl, fmt(row['weight_percent'], 3, lang) + r'\%',
                                 rf'Phizackerley-derived; \(\Delta\)logL = {gap};\newline {tag} = {ratio}'])
        write_rows(lang, 'primary_summary.tex', primary_rows)
        alternatives = []
        interpretations = {
            'phizackerley': ('Strongest non-tied alternative; the additional KV21B placement is not scored.', 'Stärkste genetisch verschiedene Alternative; die zusätzliche KV21B-Position wird nicht bewertet.'),
            'variant': ('Exploratory later-generation alternative; especially weak in c--d.', 'Explorative Alternative einer späteren Generation; besonders schwach in c--d.'),
            'dodson': ('Derived Dodson-family structure; consistently weak across settings.', 'Abgeleitete Dodson-Struktur; unter allen Einstellungen durchgehend schwach.'),
            'tawfik': ('Published later-generation alternative; consistently weak across settings.', 'Publizierte Alternative einer späteren Generation; unter allen Einstellungen durchgehend schwach.'),
        }
        for key, label in (('phizackerley', 'Phizackerley-derived'), ('variant', 'Variant E1' if lang == 'en' else 'Variante E1'),
                           ('dodson', 'Dodson-derived'), ('tawfik', 'Tawfik')):
            rows = [row for row in ranking if row['structure'].startswith(key)]
            gaps = [-float(row['delta_log_likelihood']) for row in rows]
            alternatives.append([label, fmt(min(gaps), 2, lang) + '--' + fmt(max(gaps), 2, lang),
                                 fmt(max(float(row['weight_percent']) for row in rows), 3, lang) + r'\%',
                                 interpretations[key][lang == 'de']])
        write_rows(lang, 'primary_alternatives.tex', alternatives)
        locus_rows = []
        markers = ('D13S317', 'D7S820', 'D2S1338', 'D21S11', 'D16S539', 'D18S51', 'CSF1PO', 'FGA')
        for i, marker in enumerate(markers, 1):
            values = [float(row['delta_log_likelihood']) for row in loci if row['locus'] == f'Locus{i}']
            assert min(values) > 0
            interpretation = 'Favours the shared structure in all settings.' if lang == 'en' else 'Bevorzugt die gemeinsame Struktur unter allen Einstellungen.'
            locus_rows.append([f'L{i} / {marker}', fmt(statistics.mean(values), 2, lang),
                               fmt(min(values), 2, lang) + '--' + fmt(max(values), 2, lang), interpretation])
        write_rows(lang, 'per_locus.tex', locus_rows)
        labels = ('Yuya--Thuya--Tiye trio', r'Amenhotep III--Tiye \(\rightarrow\) KV55/KV35YL', 'Yuya-to-KV55 grandparent link') if lang == 'en' else (
            'Trio Juja--Tuja--Teje', r'Amenhotep III.--Teje \(\rightarrow\) KV55/KV35YL', 'Juja--KV55-Großelternverbindung')
        focused_rows = []
        detail_rows = []
        for label, row in zip(labels, facts['focused']):
            weight = fmt(row['omissions']['Locus2'], 3, lang)
            driver = rf'No rank flips; without L2 the relationship weight is {weight}\%.' if lang == 'en' else rf'Kein Rangwechsel; ohne L2 beträgt das Beziehungsgewicht {weight}\%.'
            focused_rows.append([label, 'Stable' if lang == 'en' else 'Stabil', fmt(100, 1, lang) + r'\%', '0/27', driver])
            detail_rows.append([label, fmt(row['baseline_weight'], 3, lang),
                                fmt(row['grid_min'], 3, lang) + '--' + fmt(row['grid_max'], 3, lang),
                                fmt(row['omissions']['Locus1'], 3, lang), weight])
        write_rows(lang, 'focused_classes.tex', focused_rows)
        if lang == 'en':
            write_rows(lang, 'focused_weights.tex', detail_rows)
        for scenario in ('abcd' if lang == 'en' else ''):
            rows = [row for row in ranking if row['scenario'] == scenario]
            labels_by_key = {'shared': 'Shared Hawass-derived/Belmonte structure', 'phizackerley': 'Phizackerley-derived',
                             'dodson': 'Dodson-derived', 'tawfik': 'Tawfik', 'variant': 'Variant E1'}
            write_rows(lang, f'supplement_ranking_{scenario}.tex', [[row['rank'], next(value for key, value in labels_by_key.items()
                        if row['structure'].startswith(key)), fmt(row['log_likelihood'], 6, lang),
                        fmt(row['weight_percent'], 3, lang) + r'\%'] for row in rows])
        macros = {}
        for row in facts['primary']:
            macros['PrimaryRatio' + row['scenario'].upper()] = integer(row['ratio'], lang)
        gaps = [row['gap'] for row in facts['primary']]
        macros['PrimaryGapMin'] = fmt(min(gaps), 2, lang).replace(',', '{,}')
        macros['PrimaryGapMax'] = fmt(max(gaps), 2, lang).replace(',', '{,}')
        macros['WeakestScenario'] = weakest['scenario']
        macros['WeakestLocus'] = weakest['omitted_locus'].replace('Locus', 'L')
        macros['WeakestGap'] = fmt(-float(weakest['nearest_non_tied_delta_log_likelihood']), 2, lang).replace(',', '{,}')
        macros['WeakestRatio'] = integer(float(weakest['nearest_non_tied_bayes_factor']), lang)
        for name, row in zip(('Trio', 'Core', 'Grandparent'), facts['focused']):
            macros[name + 'WithoutLTwo'] = fmt(row['omissions']['Locus2'], 3, lang)
        (OUT / lang / 'numerical_macros.tex').write_text('\n'.join('\\newcommand{\\' + key + '}{' + value + '}' for key, value in macros.items()) + '\n')
    (OUT / 'numerical_facts.json').write_text(json.dumps(facts, indent=2) + '\n')
    print(json.dumps(facts, indent=2))


if __name__ == '__main__':
    main()
