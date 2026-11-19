"""Tests for noise_theory.py. Run: python3 -m unittest -v test_noise_theory

Eq. 3 is checked against the exact linear-noise result, computed independently here
by solving the Lyapunov equation J C + C J^T + D = 0 with numpy (no scipy needed).
"""

import math
import unittest

import numpy as np

from gillespie import gillespie
from noise_theory import eq1, eq3, feedback_rate, feedback_steady_state, moments, predict

LN2 = math.log(2)
GM, GP = LN2 / 120, LN2 / 3600   # the paper's mRNA and protein half-lives


def lyapunov_fano(r, gamma_m, beta, gamma_p, K, n, r_basal):
    """Protein Fano factor of the feedback model linearized around its steady state."""
    p = feedback_steady_state(r, gamma_m, beta, gamma_p, K, n, r_basal)
    m = p * gamma_p / beta
    h = 1e-6 * p
    slope = (feedback_rate(p + h, r, K, n, r_basal) - feedback_rate(p - h, r, K, n, r_basal)) / (2 * h)
    J = np.array([[-gamma_m, slope], [beta, -gamma_p]])
    D = np.diag([2 * gamma_m * m, 2 * gamma_p * p])  # creation = decay at steady state
    # vec(J C + C J^T) = (I (x) J + J (x) I) vec(C)
    A = np.kron(np.eye(2), J) + np.kron(J, np.eye(2))
    C = np.linalg.solve(A, -D.reshape(-1)).reshape(2, 2)
    return C[1, 1] / p


class TestEq1(unittest.TestCase):

    def test_paper_values(self):
        mean, fano = eq1(0.01, GM, 20 * GM, GP)        # Fig. 1 caption: k_R, b = 20
        self.assertAlmostEqual(mean, 0.01 * 20 / GP)
        self.assertAlmostEqual(fano, 20 / (1 + GP / GM) + 1)
        self.assertAlmostEqual(fano, 20.3548, places=4)

    def test_fano_depends_on_burst_size_not_transcription(self):
        _, f1 = eq1(0.01, GM, 20 * GM, GP)
        _, f2 = eq1(0.08, GM, 20 * GM, GP)
        _, f3 = eq1(0.04, GM, 5 * GM, GP)
        self.assertAlmostEqual(f1, f2)
        self.assertLess(f3, f1)


class TestEq3(unittest.TestCase):

    def test_matches_linear_noise_theory(self):
        r_max = 1200 * GP / 10
        for K, n, r_basal in [(100, 2, 0), (800, 2, 0), (800, 10, 0), (800, 2, 0.005), (300, 4, 0.01)]:
            with self.subTest(K=K, n=n, r_basal=r_basal):
                _, fano = eq3(r_max, GM, 10 * GM, GP, K, n, r_basal)
                self.assertAlmostEqual(fano, lyapunov_fano(r_max, GM, 10 * GM, GP, K, n, r_basal), places=5)

    def test_steady_state_is_self_consistent(self):
        r_max = 1200 * GP / 10
        p, _ = eq3(r_max, GM, 10 * GM, GP, 800, 2)
        self.assertAlmostEqual(p, feedback_rate(p, r_max, 800, 2, 0) * 10 / GP, places=6)
        self.assertAlmostEqual(p, 689.0, places=0)

    def test_weak_feedback_approaches_eq1(self):
        r_max = 1200 * GP / 10
        _, fano = eq3(r_max, GM, 10 * GM, GP, K=1e9, n=2)
        self.assertAlmostEqual(fano, eq1(r_max, GM, 10 * GM, GP)[1], places=6)

    def test_feedback_lowers_noise_at_the_same_mean(self):
        r_max = 1200 * GP / 10
        p, fano_fb = eq3(r_max, GM, 10 * GM, GP, 800, 2)
        mean_nf, fano_nf = eq1(p * GP / 10, GM, 10 * GM, GP)
        self.assertAlmostEqual(mean_nf, p)
        self.assertLess(fano_fb, fano_nf)


class TestMoments(unittest.TestCase):
    """The paper's time-dependent theory (mean and variance over time)."""

    def test_reproduces_fig2_exact_values(self):
        # exact Fano factors printed by replicate_fig2.py (scipy LSODA, rtol 1e-10)
        expected = {120: {100: 7.93, 500: 23.94, 1000: 30.31, 2000: 31.70, 3000: 30.50, 5000: 27.65},
                    24: {100: 24.68, 500: 36.42, 1000: 36.28, 2000: 34.08, 3000: 31.90, 5000: 28.44}}
        times = np.linspace(0, 5000, 101)
        for half_life, points in expected.items():
            gm = LN2 / half_life
            _, mean_p, _, var_p = moments(0.01, gm, 20 * gm, GP, times)
            for t, fano in points.items():
                with self.subTest(half_life=half_life, t=t):
                    k = int(t / 50)
                    self.assertAlmostEqual(var_p[k] / mean_p[k], fano, delta=0.005)

    def test_mean_matches_closed_form(self):
        gm, gp, r, beta = 0.5, 0.05, 1.0, 10.0
        t = np.linspace(0, 100, 201)
        mean_m, mean_p, _, _ = moments(r, gm, beta, gp, t)
        np.testing.assert_allclose(mean_m, r / gm * (1 - np.exp(-gm * t)), rtol=1e-6, atol=1e-9)
        exact_p = r * beta / (gm * gp) * (1 - (gp * np.exp(-gm * t) - gm * np.exp(-gp * t)) / (gp - gm))
        np.testing.assert_allclose(mean_p, exact_p, rtol=1e-6, atol=1e-6)

    def test_long_time_limit_is_eq1(self):
        t = np.linspace(0, 400, 101)   # 20 protein lifetimes
        mean_m, mean_p, var_m, var_p = moments(1.0, 0.5, 10.0, 0.05, t)
        mean, fano = eq1(1.0, 0.5, 10.0, 0.05)
        self.assertAlmostEqual(mean_p[-1], mean, delta=1e-3 * mean)
        self.assertAlmostEqual(var_p[-1] / mean_p[-1], fano, delta=1e-3 * fano)
        self.assertAlmostEqual(var_m[-1] / mean_m[-1], 1.0, places=6)  # mRNA is Poisson

    def test_starts_with_no_spread(self):
        mean_m, mean_p, var_m, var_p = moments(1.0, 0.5, 10.0, 0.05, np.linspace(0, 10, 11), m0=3, p0=50)
        self.assertEqual((mean_m[0], mean_p[0], var_m[0], var_p[0]), (3.0, 50.0, 0.0, 0.0))

    def test_agrees_with_many_simulated_cells(self):
        # 2000 independent cells at t = 20, while the noise is still building up
        r, gm, beta, gp, t_end, cells = 1.0, 0.5, 10.0, 0.05, 20.0, 2000
        ps = np.array([gillespie(r, gm, beta, gp, t_end, seed=k, n_points=2)[2][-1] for k in range(cells)])
        _, mean_p, _, var_p = moments(r, gm, beta, gp, np.array([0.0, t_end]))
        self.assertAlmostEqual(ps.mean(), mean_p[-1], delta=4 * ps.std() / math.sqrt(cells))
        self.assertAlmostEqual(ps.var(), var_p[-1], delta=0.12 * var_p[-1])  # SE of a variance ~ 3%


class TestPredict(unittest.TestCase):

    def test_dispatch_and_inapplicable_cases(self):
        self.assertEqual(predict(0.01, GM, 20 * GM, GP), eq1(0.01, GM, 20 * GM, GP))
        self.assertEqual(predict(0.02, GM, 10 * GM, GP, (800, 2, 0)), eq3(0.02, GM, 10 * GM, GP, 800, 2, 0))
        for args in [(0, GM, 1, GP), (1, 0, 1, GP), (1, GM, 0, GP), (1, GM, 1, 0)]:
            self.assertIsNone(predict(*args))


if __name__ == "__main__":
    unittest.main()
