"""Tests for gillespie.py against the analytic steady state. Run: python3 -m unittest -v"""

import unittest

import numpy as np

from gillespie import gillespie, steady_state_stats

R, GAMMA_M, BETA, GAMMA_P = 1, 0.5, 1, 0.5
T_MAX, T_BURN, SEED = 20000, 100, 1


class TestGillespie(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.stats = steady_state_stats(R, GAMMA_M, BETA, GAMMA_P, T_MAX, T_BURN, seed=SEED)

    def test_mrna_matches_analytic(self):
        self.assertAlmostEqual(self.stats["mean_m"], R / GAMMA_M, delta=0.05 * R / GAMMA_M)  # 2
        self.assertAlmostEqual(self.stats["fano_m"], 1, delta=0.1)

    def test_protein_matches_analytic(self):
        expected_mean = R * BETA / (GAMMA_M * GAMMA_P)   # 4
        expected_fano = 1 + BETA / (GAMMA_M + GAMMA_P)   # 2
        self.assertAlmostEqual(self.stats["mean_p"], expected_mean, delta=0.05 * expected_mean)
        self.assertAlmostEqual(self.stats["fano_p"], expected_fano, delta=0.1 * expected_fano)

    def test_trajectory_shape_and_values(self):
        times, ms, ps = gillespie(R, GAMMA_M, BETA, GAMMA_P, 100, seed=SEED, n_points=501)
        self.assertEqual(len(times), 501)
        self.assertEqual(times[0], 0)
        self.assertEqual(times[-1], 100)
        self.assertEqual((ms[0], ps[0]), (0, 0))
        self.assertTrue((ms >= 0).all() and (ps >= 0).all())

    def test_trajectory_starts_from_initial_counts(self):
        _, ms, ps = gillespie(R, GAMMA_M, BETA, GAMMA_P, 10, m0=5, p0=20, seed=SEED)
        self.assertEqual((ms[0], ps[0]), (5, 20))

    def test_same_seed_same_trajectory(self):
        a = gillespie(R, GAMMA_M, BETA, GAMMA_P, 100, seed=7)
        b = gillespie(R, GAMMA_M, BETA, GAMMA_P, 100, seed=7)
        for x, y in zip(a, b):
            np.testing.assert_array_equal(x, y)

    def test_beta_zero_gives_no_protein(self):
        _, _, ps = gillespie(R, GAMMA_M, 0, GAMMA_P, 1000, seed=SEED)
        self.assertTrue((ps == 0).all())

    def test_no_reactions_possible(self):
        times, ms, ps = gillespie(0, GAMMA_M, BETA, GAMMA_P, 10, n_points=3)
        np.testing.assert_array_equal(times, [0, 5, 10])
        np.testing.assert_array_equal(ms, [0, 0, 0])
        np.testing.assert_array_equal(ps, [0, 0, 0])

    def test_negative_feedback_lowers_protein(self):
        # protein represses its own transcription with K = 2, n = 2, r_basal = 0.001
        s = steady_state_stats(R, GAMMA_M, BETA, GAMMA_P, 5000, T_BURN, seed=SEED,
                               autorepression=(2, 2, 0.001))
        self.assertLess(s["mean_p"], R * BETA / (GAMMA_M * GAMMA_P))

    def test_rejects_invalid_inputs(self):
        with self.assertRaises(ValueError):
            gillespie(R, -1, BETA, GAMMA_P, 10)
        with self.assertRaises(ValueError):
            gillespie(R, GAMMA_M, BETA, GAMMA_P, 0)
        with self.assertRaises(ValueError):
            gillespie(-1, GAMMA_M, BETA, GAMMA_P, 10)
        with self.assertRaises(ValueError):
            gillespie(R, GAMMA_M, BETA, GAMMA_P, 10, n_points=1)
        with self.assertRaises(ValueError):
            steady_state_stats(R, GAMMA_M, BETA, GAMMA_P, 10, t_burn=10)
        with self.assertRaises(ValueError):
            gillespie(R, GAMMA_M, BETA, GAMMA_P, 10, autorepression=(0, 2, 0.001))
        with self.assertRaises(ValueError):
            gillespie(R, GAMMA_M, BETA, GAMMA_P, 10, autorepression=(2, 2, 1))


if __name__ == "__main__":
    unittest.main()
