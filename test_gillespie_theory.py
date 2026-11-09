"""Check that gillespie.py produces the results theory predicts.

Each test compares the simulation with an exact mathematical answer:
  1. Pure production: the count after time T is Poisson with mean r*T.
  2. Pure decay: the average count falls as N * exp(-gamma * t).
  3. Approach to steady state: the average over many cells follows the ODE
     solution from gene_expression.py's model.
  4. Steady state for several parameter sets: means and Fano factors match
     the analytic formulas, including a high-noise (bursty) case.
  5. The whole mRNA distribution at steady state is Poisson(r / gamma_m).

The simulation is random, so each check allows a margin of a few standard
errors (or a small fixed tolerance). Fixed seeds keep the results repeatable.

Run: python3 -m unittest test_gillespie_theory -v
"""

import math
import unittest

import numpy as np

from gillespie import gillespie, steady_state_stats

SEED = 2026
N_CELLS = 2000   # independent trajectories for ensemble tests
Z = 4            # allowed difference, in standard errors


def run_cells(n_cells, *args, **kwargs):
    """Simulate n_cells independent trajectories; return times, ms, ps (cells x times)."""
    ms, ps = [], []
    for i in range(n_cells):
        times, m, p = gillespie(*args, seed=SEED + i, **kwargs)
        ms.append(m)
        ps.append(p)
    return times, np.array(ms), np.array(ps)


def ode_mean(t, r, gamma_m, beta, gamma_p):
    """Exact ODE solution from m(0) = p(0) = 0 (requires gamma_m != gamma_p)."""
    m = r / gamma_m * (1 - np.exp(-gamma_m * t))
    p = r * beta / gamma_m * ((1 - np.exp(-gamma_p * t)) / gamma_p
                              - (np.exp(-gamma_m * t) - np.exp(-gamma_p * t)) / (gamma_p - gamma_m))
    return m, p


class TestAgainstTheory(unittest.TestCase):
    def assert_close_to_mean(self, samples, expected, label):
        """Sample average within Z standard errors of the expected value."""
        se = samples.std(ddof=1) / math.sqrt(len(samples))
        margin = Z * se + 1e-9
        self.assertLess(abs(samples.mean() - expected), margin,
                        f"{label}: simulated {samples.mean():.4f}, expected {expected:.4f}")

    def test_pure_production_is_poisson(self):
        # only transcription can happen, so m(T) ~ Poisson(r * T): mean = variance = r*T
        r, T = 3.0, 5.0
        _, ms, _ = run_cells(N_CELLS, r, 0, 0, 0, T, n_points=2)
        final = ms[:, -1].astype(float)
        self.assert_close_to_mean(final, r * T, "mean of m(T)")
        self.assertAlmostEqual(final.var(ddof=1) / final.mean(), 1.0, delta=0.1)

    def test_pure_decay_is_exponential(self):
        # only mRNA decay can happen: average m(t) = N * exp(-gamma_m * t)
        N, gamma_m = 50, 0.3
        times, ms, _ = run_cells(N_CELLS, 0, gamma_m, 0, 0, 10, m0=N, n_points=11)
        for k, t in enumerate(times):
            with self.subTest(t=t):
                self.assert_close_to_mean(ms[:, k].astype(float), N * math.exp(-gamma_m * t),
                                          f"mean m at t={t}")

    def test_average_follows_ode(self):
        # starting from zero, the average over many cells follows the ODE solution
        r, gamma_m, beta, gamma_p = 2.0, 1.0, 3.0, 0.2
        times, ms, ps = run_cells(N_CELLS, r, gamma_m, beta, gamma_p, 20, n_points=11)
        m_ode, p_ode = ode_mean(times, r, gamma_m, beta, gamma_p)
        for k, t in enumerate(times[1:], start=1):
            with self.subTest(t=t):
                self.assert_close_to_mean(ms[:, k].astype(float), m_ode[k], f"mean m at t={t}")
                self.assert_close_to_mean(ps[:, k].astype(float), p_ode[k], f"mean p at t={t}")

    def test_steady_state_for_several_parameter_sets(self):
        cases = [
            # r,   gamma_m, beta, gamma_p
            (1.0,  0.5,     1.0,  0.5),    # the default example
            (5.0,  1.0,     2.0,  0.1),    # high expression: mean p = 100
            (0.5,  1.0,     10.0, 0.1),    # rare, bursty mRNA: Fano p about 10
        ]
        for r, gamma_m, beta, gamma_p in cases:
            with self.subTest(r=r, gamma_m=gamma_m, beta=beta, gamma_p=gamma_p):
                s = steady_state_stats(r, gamma_m, beta, gamma_p, t_max=50000,
                                       t_burn=200, seed=SEED)
                mean_m = r / gamma_m
                mean_p = r * beta / (gamma_m * gamma_p)
                fano_p = 1 + beta / (gamma_m + gamma_p)
                self.assertAlmostEqual(s["mean_m"], mean_m, delta=0.05 * mean_m)
                self.assertAlmostEqual(s["fano_m"], 1.0, delta=0.1)
                self.assertAlmostEqual(s["mean_p"], mean_p, delta=0.08 * mean_p)
                self.assertAlmostEqual(s["fano_p"], fano_p, delta=0.15 * fano_p)

    def test_mrna_distribution_is_poisson(self):
        # at steady state, P(m = k) = e^-mu * mu^k / k!  with mu = r / gamma_m
        r, gamma_m, beta, gamma_p = 3.0, 1.0, 1.0, 0.5
        mu = r / gamma_m
        _, ms, _ = gillespie(r, gamma_m, beta, gamma_p, t_max=100000,
                             seed=SEED, n_points=100001)
        ms = ms[100:]  # skip the first 100 time units (burn-in)
        for k in range(10):
            with self.subTest(k=k):
                simulated = np.mean(ms == k)
                expected = math.exp(-mu) * mu ** k / math.factorial(k)
                self.assertAlmostEqual(simulated, expected, delta=0.01,
                                       msg=f"P(m={k}): simulated {simulated:.4f}, "
                                           f"expected {expected:.4f}")


if __name__ == "__main__":
    unittest.main()
