"""Gillespie stochastic simulation of the mRNA/protein model in gene_expression.py.

Reactions (m = mRNA count, p = protein count):
    1. transcription     0 -> m      rate r            (or r(p) with autorepression)
    2. mRNA decay        m -> 0      rate gamma_m * m
    3. translation       m -> m + p  rate beta * m
    4. protein decay     p -> 0      rate gamma_p * p

Expected steady state (constant r):
    mean m = r / gamma_m,                    Fano m = 1
    mean p = r * beta / (gamma_m * gamma_p), Fano p = 1 + beta / (gamma_m + gamma_p)

Speed: the simulation loop is compiled with numba when it is installed (~60x
faster after a one-time compile of about 1 s, cached on disk); otherwise the
same code runs as plain Python.
"""

import math
import random

import numpy as np

try:
    from numba import njit
except ImportError:  # plain Python fallback
    def njit(**options):
        return lambda f: f


def _check_nonnegative(name, value):
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value}")


def _check_positive(name, value):
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value}")


@njit(cache=True)
def _simulate(r, feedback, K, n, r_basal, gamma_m, beta, gamma_p,
              m0, p0, t_burn, grid, seed):
    """Run one trajectory up to grid[-1].

    Records the state at each time in grid, and accumulates time-weighted
    sums after t_burn so mean and variance need no second pass.
    """
    if seed >= 0:
        random.seed(seed)
    t_max = grid[-1]
    ms = np.empty(len(grid), np.int64)
    ps = np.empty(len(grid), np.int64)
    k = 0  # next grid index to fill
    t = 0.0
    m = m0
    p = p0
    total = sum_m = sum_m2 = sum_p = sum_p2 = 0.0
    Kn = K ** n

    while True:
        if feedback:
            a1 = r_basal + (r - r_basal) * Kn / (Kn + p ** n)
        else:
            a1 = r
        a2 = gamma_m * m
        a3 = beta * m
        a4 = gamma_p * p
        a_total = a1 + a2 + a3 + a4

        # time of next reaction: exponential with rate a_total
        if a_total > 0:
            t_next = t - math.log(1.0 - random.random()) / a_total
        else:
            t_next = t_max  # nothing can happen any more
        if t_next > t_max:
            t_next = t_max

        # the current state holds on [t, t_next)
        while k < len(grid) and grid[k] < t_next:
            ms[k] = m
            ps[k] = p
            k += 1
        if t_next > t_burn:
            dt = t_next - max(t, t_burn)
            total += dt
            sum_m += m * dt
            sum_m2 += m * m * dt
            sum_p += p * dt
            sum_p2 += p * p * dt

        t = t_next
        if t >= t_max:
            break

        # pick which reaction happens, weighted by its propensity
        u = random.random() * a_total
        if u < a1:
            m += 1
        elif u < a1 + a2:
            m -= 1
        elif u < a1 + a2 + a3:
            p += 1
        else:
            p -= 1

    while k < len(grid):  # grid points at exactly t_max
        ms[k] = m
        ps[k] = p
        k += 1
    sums = np.array([total, sum_m, sum_m2, sum_p, sum_p2])
    return ms, ps, sums


def _run(r, gamma_m, beta, gamma_p, grid, m0, p0, seed, autorepression, t_burn):
    for name, value in (("r", r), ("gamma_m", gamma_m), ("beta", beta),
                        ("gamma_p", gamma_p), ("m0", m0), ("p0", p0)):
        _check_nonnegative(name, value)
    _check_positive("t_max", grid[-1])
    if not 0 <= t_burn < grid[-1]:
        raise ValueError(f"t_burn must be in [0, t_max), got {t_burn}")

    if autorepression is None:
        feedback, K, n, r_basal = False, 1.0, 1.0, 0.0
    else:
        K, n, r_basal = autorepression
        _check_positive("K", K)
        _check_positive("n", n)
        _check_nonnegative("r_basal", r_basal)
        if r_basal >= r:
            raise ValueError(f"r_basal ({r_basal}) must be less than r ({r})")
        feedback = True

    return _simulate(float(r), feedback, float(K), float(n), float(r_basal),
                     float(gamma_m), float(beta), float(gamma_p),
                     int(m0), int(p0), float(t_burn), grid,
                     -1 if seed is None else int(seed))


def gillespie(r, gamma_m, beta, gamma_p, t_max, m0=0, p0=0, seed=None,
              autorepression=None, n_points=1001):
    """Simulate one trajectory and sample it at n_points evenly spaced times.

    r: transcription rate. With autorepression = (K, n, r_basal), the protein
       represses its own gene: r(p) = r_basal + (r - r_basal) * K^n / (K^n + p^n),
       the Level 1 repressor formula with r as r_max.
    seed: set for reproducible runs.

    Returns (times, ms, ps) as numpy arrays of length n_points.
    """
    if n_points < 2:
        raise ValueError(f"n_points must be at least 2, got {n_points}")
    _check_positive("t_max", t_max)
    grid = np.linspace(0.0, t_max, n_points)
    ms, ps, _ = _run(r, gamma_m, beta, gamma_p, grid, m0, p0, seed, autorepression, 0.0)
    return grid, ms, ps


def steady_state_stats(r, gamma_m, beta, gamma_p, t_max, t_burn=0.0, m0=0, p0=0,
                       seed=None, autorepression=None):
    """Time-weighted mean, variance, and Fano factor (variance / mean) of m and p
    over [t_burn, t_max] of one long trajectory.

    Each state is weighted by how long the system stayed in it, which is the
    right way to average a Gillespie trajectory.
    """
    _check_positive("t_max", t_max)
    grid = np.array([float(t_max)])
    _, _, sums = _run(r, gamma_m, beta, gamma_p, grid, m0, p0, seed, autorepression, t_burn)
    total, sum_m, sum_m2, sum_p, sum_p2 = sums

    stats = {}
    for name, s, s2 in (("m", sum_m, sum_m2), ("p", sum_p, sum_p2)):
        mean = s / total
        variance = s2 / total - mean ** 2
        stats[f"mean_{name}"] = mean
        stats[f"var_{name}"] = variance
        stats[f"fano_{name}"] = variance / mean if mean > 0 else float("nan")
    return stats


if __name__ == "__main__":
    import time

    r, gamma_m, beta, gamma_p = 1, 0.5, 1, 0.5
    start = time.perf_counter()
    steady_state_stats(r, gamma_m, beta, gamma_p, t_max=10, seed=1)
    print(f"first call (numba compile or cache load): {time.perf_counter() - start:.3f} s")

    start = time.perf_counter()
    s = steady_state_stats(r, gamma_m, beta, gamma_p, t_max=20000, t_burn=100, seed=1)
    print(f"t_max = 20000 run: {time.perf_counter() - start:.4f} s")
    print(f"mRNA:    mean = {s['mean_m']:.3f} (expected {r / gamma_m:.3f}), "
          f"Fano = {s['fano_m']:.3f} (expected 1.000)")
    print(f"protein: mean = {s['mean_p']:.3f} (expected {r * beta / (gamma_m * gamma_p):.3f}), "
          f"Fano = {s['fano_p']:.3f} (expected {1 + beta / (gamma_m + gamma_p):.3f})")
