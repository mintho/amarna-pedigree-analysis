import json
import math
from pathlib import Path
import unittest

from src.genetics import observation_likelihood_from_observation_model
from src.likelihood_enumeration import (
    background_frequencies_for_locus, locus_log_likelihood, load_observation_model,
)
from src.likelihood_factors import (
    build_locus_factors, locus_likelihood_factor_graph,
    variable_elimination, variable_elimination_sparse,
)

ROOT = Path(__file__).resolve().parents[2]


class FixedBackgroundTests(unittest.TestCase):
    def test_founder_changes_do_not_change_observation_probability(self):
        model = {'type': 'generative_str', 'dropout': 0.2, 'background_error': 0.01,
                 'background_allele_frequencies': {'Locus1': {'A': 0.3, 'B': 0.7}}}
        probabilities = []
        for founder in ({'A': 0.9, 'B': 0.1}, {'A': 0.1, 'B': 0.9}):
            q = background_frequencies_for_locus(model, 'Locus1', founder)
            probabilities.append(observation_likelihood_from_observation_model(['B'], ('A', 'A'), model, q))
        self.assertEqual(probabilities[0], probabilities[1])

    def test_missing_fixed_background_locus_is_rejected(self):
        with self.assertRaisesRegex(Exception, 'Missing fixed background'):
            background_frequencies_for_locus({'background_allele_frequencies': {}}, 'Locus1', {'A': 1})

    def test_exact_engines_agree_with_distinct_background_and_founder(self):
        inputs = dict(person_ids=['P', 'Q', 'C'],
                      relationships=[{'child': 'C', 'father': 'P', 'mother': 'Q'}],
                      observations={'C': {'Locus1': {'alleles': ['A', 'B']}}},
                      locus='Locus1', allele_frequencies={'A': 0.8, 'B': 0.2},
                      observation_model={'type': 'generative_str', 'dropout': 0.2,
                                         'background_error': 0.01,
                                         'background_allele_frequencies': {'Locus1': {'A': 0.3, 'B': 0.7}}})
        logp, _ = locus_log_likelihood(**inputs)
        probability, _ = locus_likelihood_factor_graph(**inputs, background_error=0.01)
        self.assertAlmostEqual(logp, math.log(probability), places=12)

    def test_current_runs_use_same_fixed_reference(self):
        backgrounds = []
        for path in sorted((ROOT / 'Data/models/runs').glob('*.json')):
            run = json.loads(path.read_text())
            backgrounds.append(load_observation_model(path, run)['background_allele_frequencies'])
        self.assertEqual(len(backgrounds), 7)
        self.assertTrue(all(q == backgrounds[0] for q in backgrounds))
        for path in (ROOT / 'Data/models/founder_scenarios').glob('*scenario_*_all_loci.json'):
            self.assertEqual(json.loads(path.read_text())['allele_frequencies'], backgrounds[0])

    def test_accelerated_contraction_matches_sparse_reference(self):
        for dropout in (0.0, 0.2, 0.3):
            inputs = dict(person_ids=['G', 'P', 'Q', 'C'],
                          relationships=[{'child': 'P', 'parent': 'G'},
                                         {'child': 'C', 'father': 'P', 'mother': 'Q'}],
                          observations={'G': {'Locus1': {'alleles': ['A']}},
                                        'C': {'Locus1': {'alleles': ['A', 'B']}}},
                          locus='Locus1', allele_frequencies={'A': 0.8, 'B': 0.2},
                          background_error=0.01,
                          observation_model={'type': 'generative_str', 'dropout': dropout,
                                             'background_error': 0.01,
                                             'background_allele_frequencies': {'Locus1': {'A': 0.3, 'B': 0.7}}})
            factors, _, _ = build_locus_factors(**inputs)
            accelerated, _ = variable_elimination(factors, inputs['person_ids'])
            reference, _ = variable_elimination_sparse(factors, inputs['person_ids'])
            enumeration, _ = locus_log_likelihood(**{k: v for k, v in inputs.items() if k != 'background_error'})
            self.assertAlmostEqual(accelerated, reference, places=13)
            if accelerated == 0:
                self.assertEqual(enumeration, float('-inf'))
            else:
                self.assertAlmostEqual(math.log(accelerated), enumeration, places=12)
