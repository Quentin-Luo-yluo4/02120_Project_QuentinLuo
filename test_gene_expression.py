"""Tests using the special cases from the project proposal. Run: python3 -m unittest -v"""

import unittest

from gene_expression import (
    ACTIVATOR,
    REPRESSOR,
    hill,
    mrna_derivative,
    mrna_steady_state,
    protein_derivative,
    protein_steady_state,
    transcription_rate,
)

K, N = 10, 2


class TestTranscriptionRate(unittest.TestCase):
    def test_repressor_values(self):
        self.assertAlmostEqual(transcription_rate([(0, K, N, REPRESSOR)]), 1.0)
        self.assertAlmostEqual(transcription_rate([(K, K, N, REPRESSOR)]), 0.5005)
        self.assertAlmostEqual(transcription_rate([(1000, K, N, REPRESSOR)]), 0.0011, places=4)

    def test_activator_values(self):
        self.assertAlmostEqual(transcription_rate([(0, K, N, ACTIVATOR)]), 0.001)
        self.assertAlmostEqual(transcription_rate([(K, K, N, ACTIVATOR)]), 0.5005)
        self.assertAlmostEqual(transcription_rate([(1000, K, N, ACTIVATOR)]), 1.0, places=3)

    def test_multiple_factors_multiply(self):
        # activator and repressor both at x = K: product of h = 0.25
        factors = [(K, K, N, ACTIVATOR), (K, K, N, REPRESSOR)]
        self.assertAlmostEqual(transcription_rate(factors), 0.001 + 0.999 * 0.25)

    def test_rejects_invalid_inputs(self):
        for kind in (REPRESSOR, ACTIVATOR):
            with self.assertRaises(ValueError):
                hill(-1, K, N, kind)
            with self.assertRaises(ValueError):
                hill(1, 0, N, kind)
            with self.assertRaises(ValueError):
                hill(1, K, -2, kind)
        with self.assertRaises(ValueError):
            hill(1, K, N, "enhancer")
        with self.assertRaises(ValueError):
            transcription_rate([(1, K, N, REPRESSOR)], r_max=1, r_basal=1)
        with self.assertRaises(ValueError):
            transcription_rate([(1, K, N, REPRESSOR)], r_max=-1, r_basal=-2)
        with self.assertRaises(ValueError):
            transcription_rate([])


class TestMrnaAndProtein(unittest.TestCase):
    r, gamma_m, beta, gamma_p = 1, 0.5, 1, 0.5

    def test_steady_states(self):
        m = mrna_steady_state(self.r, self.gamma_m)
        self.assertAlmostEqual(m, 2)
        self.assertAlmostEqual(protein_steady_state(m, self.beta, self.gamma_p), 4)

    def test_derivatives_zero_at_steady_state(self):
        self.assertAlmostEqual(mrna_derivative(2, self.r, self.gamma_m), 0)
        self.assertAlmostEqual(protein_derivative(4, 2, self.beta, self.gamma_p), 0)

    def test_beta_zero(self):
        m = mrna_steady_state(self.r, self.gamma_m)
        self.assertAlmostEqual(m, 2)
        self.assertAlmostEqual(protein_steady_state(m, 0, self.gamma_p), 0)

    def test_rejects_negative_rates(self):
        with self.assertRaises(ValueError):
            mrna_derivative(1, -1, self.gamma_m)
        with self.assertRaises(ValueError):
            protein_derivative(1, 1, -1, self.gamma_p)
        with self.assertRaises(ValueError):
            mrna_steady_state(self.r, 0)


if __name__ == "__main__":
    unittest.main()
