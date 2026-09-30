import unittest

from src.genetics import (
    authentic_scored_observation_likelihood,
    background_call_distribution_probability,
    canonical_genotype,
    dropout_observation_likelihood,
    founder_prior_from_allele_frequencies,
    gamete_distribution,
    generative_str_observation_likelihood,
    genotype_key,
    genotype_prior_from_mapping,
    mendelian_child_distribution,
    observation_likelihood_from_observation_model,
    parse_genotype_key,
    unordered_genotypes,
)


class GeneticsTests(unittest.TestCase):
    def assertDistributionAlmostEqual(self, actual, expected):
        self.assertEqual(set(actual), set(expected))
        for key, value in expected.items():
            self.assertAlmostEqual(actual[key], value)

    def test_canonical_genotype_sorts_alleles(self):
        self.assertEqual(canonical_genotype("B", "A"), ("A", "B"))

    def test_genotype_key_roundtrip(self):
        self.assertEqual(genotype_key(("B", "A")), "A/B")
        self.assertEqual(parse_genotype_key("B/A"), ("A", "B"))

    def test_genotype_prior_from_mapping(self):
        prior = genotype_prior_from_mapping({"B/A": 0.25, "A/B": 0.25, "B/B": 0.5})
        self.assertDistributionAlmostEqual(
            prior,
            {
                ("A", "B"): 0.5,
                ("B", "B"): 0.5,
            },
        )

    def test_unordered_genotypes(self):
        self.assertEqual(
            unordered_genotypes(["B", "A"]),
            [("A", "A"), ("A", "B"), ("B", "B")],
        )

    def test_founder_prior_from_equal_allele_frequencies(self):
        prior = founder_prior_from_allele_frequencies({"A": 0.5, "B": 0.5})
        self.assertDistributionAlmostEqual(
            prior,
            {
                ("A", "A"): 0.25,
                ("A", "B"): 0.5,
                ("B", "B"): 0.25,
            },
        )

    def test_gamete_distribution(self):
        self.assertDistributionAlmostEqual(gamete_distribution(("A", "A")), {"A": 1.0})
        self.assertDistributionAlmostEqual(
            gamete_distribution(("A", "B")),
            {"A": 0.5, "B": 0.5},
        )

    def test_mendelian_aa_x_bb(self):
        child = mendelian_child_distribution(("A", "A"), ("B", "B"))
        self.assertDistributionAlmostEqual(child, {("A", "B"): 1.0})

    def test_mendelian_ab_x_ab(self):
        child = mendelian_child_distribution(("A", "B"), ("A", "B"))
        self.assertDistributionAlmostEqual(
            child,
            {
                ("A", "A"): 0.25,
                ("A", "B"): 0.5,
                ("B", "B"): 0.25,
            },
        )

    def test_dropout_observation_likelihood_partial_heterozygote(self):
        self.assertAlmostEqual(
            dropout_observation_likelihood(["A"], ("A", "B"), 0.2),
            0.16,
        )

    def test_dropout_observation_likelihood_partial_homozygote(self):
        self.assertAlmostEqual(
            dropout_observation_likelihood(["A"], ("A", "A"), 0.2),
            0.32,
        )

    def test_dropout_observation_likelihood_full_heterozygote(self):
        self.assertAlmostEqual(
            dropout_observation_likelihood(["B", "A"], ("A", "B"), 0.2),
            0.64,
        )

    def test_dropout_observation_likelihood_unknown_stays_uninformative(self):
        self.assertEqual(dropout_observation_likelihood([], ("A", "B"), 0.2), 1.0)

    def scored_states(self, alleles):
        states = [[allele] for allele in alleles]
        for i, allele_a in enumerate(alleles):
            for allele_b in alleles[i:]:
                states.append([allele_a, allele_b])
        return states

    def test_authentic_scored_observation_distribution_normalizes(self):
        alleles = ["A", "B", "other"]
        dropout = 0.2
        total = sum(
            authentic_scored_observation_likelihood(state, ("A", "B"), dropout)
            for state in self.scored_states(alleles)
        )
        self.assertAlmostEqual(total, 1.0)

    def test_authentic_scored_observation_distribution_normalizes_homozygote(self):
        alleles = ["A", "B", "other"]
        dropout = 0.2
        total = sum(
            authentic_scored_observation_likelihood(state, ("A", "A"), dropout)
            for state in self.scored_states(alleles)
        )
        self.assertAlmostEqual(total, 1.0)

    def test_authentic_scored_partial_distinguishes_homozygote_and_heterozygote(self):
        dropout = 0.2
        denominator = 1.0 - dropout * dropout
        self.assertAlmostEqual(
            authentic_scored_observation_likelihood(["A"], ("A", "A"), dropout),
            2.0 * dropout * (1.0 - dropout) / denominator,
        )
        self.assertAlmostEqual(
            authentic_scored_observation_likelihood(["A"], ("A", "B"), dropout),
            dropout * (1.0 - dropout) / denominator,
        )
        self.assertAlmostEqual(
            authentic_scored_observation_likelihood(["B"], ("A", "B"), dropout),
            dropout * (1.0 - dropout) / denominator,
        )
        self.assertAlmostEqual(
            authentic_scored_observation_likelihood(["A", "B"], ("A", "B"), dropout),
            (1.0 - dropout) ** 2 / denominator,
        )

    def test_background_call_distribution_normalizes(self):
        allele_frequencies = {"A": 0.7, "B": 0.2, "other": 0.1}
        alleles = sorted(allele_frequencies)
        dropout = 0.2
        total = sum(
            background_call_distribution_probability(state, allele_frequencies, dropout)
            for state in self.scored_states(alleles)
        )
        self.assertAlmostEqual(total, 1.0)

    def test_background_call_distribution_uses_allele_frequencies(self):
        allele_frequencies = {"A": 0.7, "B": 0.2, "other": 0.1}
        dropout = 0.2
        probability = background_call_distribution_probability(
            ["A"],
            allele_frequencies,
            dropout,
        )
        # AA contributes 0.49 * 1/3; AB contributes 0.28 * 1/6;
        # A/other contributes 0.14 * 1/6.
        self.assertAlmostEqual(probability, 0.23333333333333334)

    def test_generative_str_missing_is_uninformative(self):
        self.assertEqual(
            generative_str_observation_likelihood(
                [],
                ("A", "B"),
                {"A": 0.5, "B": 0.5},
                dropout=0.2,
                background_error=0.001,
            ),
            1.0,
        )

    def test_generative_str_impossible_authentic_call_uses_background_error(self):
        allele_frequencies = {"A": 0.45, "B": 0.45, "C": 0.1}
        dropout = 0.2
        background_error = 0.01
        expected_background = background_call_distribution_probability(
            ["C"],
            allele_frequencies,
            dropout,
        )
        self.assertGreater(expected_background, 0.0)
        self.assertEqual(
            authentic_scored_observation_likelihood(["C"], ("A", "B"), dropout),
            0.0,
        )
        self.assertAlmostEqual(
            generative_str_observation_likelihood(
                ["C"],
                ("A", "B"),
                allele_frequencies,
                dropout=dropout,
                background_error=background_error,
            ),
            background_error * expected_background,
        )

    def test_generative_str_other_participates_in_background_normalization(self):
        allele_frequencies = {"A": 0.7, "B": 0.2, "other": 0.1}
        dropout = 0.2
        self.assertGreater(
            background_call_distribution_probability(
                ["other"],
                allele_frequencies,
                dropout,
            ),
            0.0,
        )
        self.assertGreater(
            background_call_distribution_probability(
                ["A", "other"],
                allele_frequencies,
                dropout,
            ),
            0.0,
        )

    def test_observation_model_dispatcher_uses_generative_str(self):
        probability = observation_likelihood_from_observation_model(
            ["A"],
            ("A", "B"),
            {
                "type": "generative_str",
                "dropout": 0.2,
                "background_error": 0.0,
                "background_source": "run_locus_allele_frequencies",
            },
            {"A": 0.5, "B": 0.5},
        )
        self.assertAlmostEqual(probability, 1.0 / 6.0)

if __name__ == "__main__":
    unittest.main()
