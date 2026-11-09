"""Gillespie simulation of a gene cascade (Thattai & van Oudenaarden 2001, Fig. 5).

Gene 1 (the regulator) is transcribed at a constant rate k_R1. Each later gene i is
activated by the protein of gene i-1 with a hyperbolic (n = 1) activator function:

    k_Ri = k_R_max * (p_{i-1} / K_d) / (1 + p_{i-1} / K_d)

Every gene has the same four reactions as gillespie.py: make mRNA, destroy mRNA,
make protein (rate k_P * m), destroy protein. Protein i only affects gene i+1.

Speed: compiled with numba when installed, like gillespie.py.
"""

import math
import random

import numpy as np

try:
    from numba import njit
except ImportError:  # plain Python fallback
    def njit(**options):
        return lambda f: f


@njit(cache=True)
def _simulate_cascade(k_r1, k_r_max, k_d, k_p, gamma_m, gamma_p, m0, p0, t_max, t_burn, seed):
    """One long trajectory; returns time-weighted sums of p and p^2 for every gene."""
    random.seed(seed)
    n = len(k_p)
    m = m0.copy()
    p = p0.copy()
    a = np.zeros(4 * n)
    total = 0.0
    sum_p = np.zeros(n)
    sum_p2 = np.zeros(n)
    t = 0.0
    while True:
        a_total = 0.0
        for i in range(n):
            if i == 0:
                tx = k_r1
            else:
                x = p[i - 1] / k_d
                tx = k_r_max * x / (1.0 + x)
            a[4 * i] = tx
            a[4 * i + 1] = gamma_m * m[i]
            a[4 * i + 2] = k_p[i] * m[i]
            a[4 * i + 3] = gamma_p * p[i]
            a_total += tx + a[4 * i + 1] + a[4 * i + 2] + a[4 * i + 3]

        t_next = t - math.log(1.0 - random.random()) / a_total
        if t_next > t_max:
            t_next = t_max
        if t_next > t_burn:
            dt = t_next - max(t, t_burn)
            total += dt
            for i in range(n):
                sum_p[i] += p[i] * dt
                sum_p2[i] += p[i] * p[i] * dt
        t = t_next
        if t >= t_max:
            break

        u = random.random() * a_total
        j = 0
        c = a[0]
        while c < u and j < 4 * n - 1:
            j += 1
            c += a[j]
        gene, kind = j // 4, j % 4
        if kind == 0:
            m[gene] += 1
        elif kind == 1:
            m[gene] -= 1
        elif kind == 2:
            p[gene] += 1
        else:
            p[gene] -= 1
    return total, sum_p, sum_p2


def cascade_stats(k_r1, k_r_max, k_d, b, mrna_half_life, protein_half_life,
                  t_max, t_burn, seed=0):
    """Time-weighted mean and Fano factor of each protein in the cascade.

    b: list of burst sizes, one per gene (k_P = b * gamma_m).
    Half-lives are in seconds. Starts near the deterministic steady state.
    Returns (means, fanos) as numpy arrays, one entry per gene.
    """
    gamma_m = math.log(2) / mrna_half_life
    gamma_p = math.log(2) / protein_half_life
    b = np.asarray(b, dtype=float)
    k_p = b * gamma_m

    # deterministic steady state, used as the starting point
    n = len(b)
    p_star = np.zeros(n)
    m_star = np.zeros(n)
    for i in range(n):
        tx = k_r1 if i == 0 else k_r_max * (p_star[i - 1] / k_d) / (1 + p_star[i - 1] / k_d)
        m_star[i] = tx / gamma_m
        p_star[i] = k_p[i] * m_star[i] / gamma_p
    m0 = np.round(m_star).astype(np.int64)
    p0 = np.round(p_star).astype(np.int64)

    total, s, s2 = _simulate_cascade(float(k_r1), float(k_r_max), float(k_d), k_p,
                                     gamma_m, gamma_p, m0, p0, float(t_max), float(t_burn), int(seed))
    mean = s / total
    var = s2 / total - mean ** 2
    return mean, var / mean
