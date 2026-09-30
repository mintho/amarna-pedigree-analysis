import math
from collections import Counter
import json
from pathlib import Path
import unittest

from src.summarize_primary_results import build_summary, structural_ranking_rows
from src.compare import add_model_comparison

ROOT = Path(__file__).resolve().parents[2]


class RepositoryTests(unittest.TestCase):
    def test_all_primary_profiles_match_reported_repeat_counts(self):
        data = ROOT / 'Data/data'
        original = json.loads((data / 'observed_str.reported_repeat_counts.json').read_text())['observations']
        coding = json.loads((data / 'color_to_reported_str_allele.json').read_text())['loci']
        loci = [f'Locus{i}' for i in range(1, 9)]
        expected = Counter(tuple(tuple(sorted(profile[locus])) for locus in loci)
                           for sample, profile in original.items() if sample != 'KV21B')
        for path in data.glob('observed_str.*scenario_*_all_loci.json'):
            encoded = json.loads(path.read_text())['observations']
            actual = Counter(tuple(tuple(sorted(coding[locus][token] for token in profile[locus]['alleles']))
                                   for locus in loci) for profile in encoded.values()
                             if any(profile[locus]['alleles'] for locus in loci))
            self.assertEqual(actual, expected, path.name)

    def test_historical_aliases_have_one_prior_mass(self):
        rows = add_model_comparison([
            {'pedigree': 'reading_a', 'structure_id': 'shared', 'log_likelihood': 0},
            {'pedigree': 'reading_b', 'structure_id': 'shared', 'log_likelihood': 0},
            {'pedigree': 'alternative', 'log_likelihood': -math.log(9)},
        ])
        self.assertAlmostEqual(rows[0]['weight_percent'], 90)
        self.assertAlmostEqual(rows[1]['weight_percent'], 90)
        self.assertAlmostEqual(rows[2]['weight_percent'], 10)
        self.assertEqual([row['rank'] for row in rows], [1, 1, 2])

    def test_shared_structure_with_different_likelihood_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Different likelihoods'):
            add_model_comparison([
                {'pedigree': 'a', 'structure_id': 'shared', 'log_likelihood': 0},
                {'pedigree': 'b', 'structure_id': 'shared', 'log_likelihood': -1},
            ])

    def test_unique_structural_normalisation(self):
        summary = build_summary(ROOT / 'Results/Primary_All_Locus')
        rows = structural_ranking_rows(summary)
        self.assertEqual(len(rows), 20)
        for scenario in 'abcd':
            subset = [row for row in rows if row['scenario'] == scenario]
            self.assertEqual(len(subset), 5)
            self.assertAlmostEqual(sum(row['weight_percent'] for row in subset), 100)
            self.assertEqual(subset[0]['structure'], 'shared_kv55_kv35yl_tutankhamun_structure')

    def test_reported_population_prior_support(self):
        summary = build_summary(ROOT / 'Results/Primary_All_Locus')
        for scenario in summary['scenarios']:
            self.assertAlmostEqual(scenario['tie_delta_log_likelihood'], 0)
            gap = scenario['nearest_non_tied_delta_log_likelihood']
            self.assertLess(gap, 0)
            self.assertAlmostEqual(math.exp(-gap), scenario['nearest_non_tied_bayes_factor'])
