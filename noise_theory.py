"""The paper's analytic predictions for protein noise (Thattai & van Oudenaarden 2001),
written in the same parameters as gillespie.py:

    r        transcription rate (with feedback: the unrepressed rate r_max)
    gamma_m  mRNA decay rate          beta     translation rate
    gamma_p  protein decay rate       b = beta / gamma_m  (burst size)

Eq. 1 (no feedback):   mean p = r * b / gamma_p
                       Fano   = b / (1 + eta) + 1,          eta = gamma_p / gamma_m

Eq. 3 (protein represses its own gene, r(p) = r_basal + (r - r_basal) K^n / (K^n + p^n)):
    linearize r(p) ~ k0 - k1 * p around the steady state p*, phi = k1 / gamma_p, then
                       mean p = p*   (the self-consistent level p* = r(p*) * b / gamma_p)
                       Fano   = (1 - phi) / (1 + b * phi) * b / (1 + eta) + 1
Eq. 3 is exact only for small fluctuations; for steep feedback (large n) or strong
repression the simulation is somewhat noisier than Eq. 3 (see replicate_figures.py).

Uses numpy but not scipy (scipy is not installed on Vercel).
"""

import math

import numpy as np


def eq1(r, gamma_m, beta, gamma_p):
    """Mean and Fano factor of protein without feedback (paper's Eq. 1)."""
    b = beta / gamma_m
    eta = gamma_p / gamma_m
    return r * b / gamma_p, b / (1 + eta) + 1


def feedback_rate(p, r, K, n, r_basal):
    """Transcription rate when the protein represses its own gene (as in gillespie.py)."""
    return r_basal + (r - r_basal) / (1 + (p / K) ** n)


def feedback_steady_state(r, gamma_m, beta, gamma_p, K, n, r_basal):
    """Solve p = r(p) * b / gamma_p by bisection. r(p) falls as p rises, so the
    root is unique and lies between r_basal * b / gamma_p and r * b / gamma_p."""
    scale = beta / (gamma_m * gamma_p)
    lo, hi = r_basal * scale, r * scale
    for _ in range(200):
        mid = (lo + hi) / 2
        if feedback_rate(mid, r, K, n, r_basal) * scale > mid:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def eq3(r, gamma_m, beta, gamma_p, K, n, r_basal=0.0):
    """Mean and Fano factor of protein with negative autoregulation (paper's Eq. 3)."""
    b = beta / gamma_m
    eta = gamma_p / gamma_m
    p = feedback_steady_state(r, gamma_m, beta, gamma_p, K, n, r_basal)
    if p > 0:
        x = (p / K) ** n
        k1 = (r - r_basal) * n * x / (p * (1 + x) ** 2)   # -dr/dp at p
    else:
        k1 = 0.0
    phi = k1 / gamma_p
    return p, (1 - phi) / (1 + b * phi) * b / (1 + eta) + 1


def predict(r, gamma_m, beta, gamma_p, autorepression=None):
    """Paper's (mean, Fano) for these parameters; autorepression = (K, n, r_basal) or None.
    Returns None when the formulas don't apply (a decay rate of 0, or no protein)."""
    if gamma_m <= 0 or gamma_p <= 0 or r <= 0 or beta <= 0:
        return None
    if autorepression is None:
        return eq1(r, gamma_m, beta, gamma_p)
    K, n, r_basal = autorepression
    return eq3(r, gamma_m, beta, gamma_p, K, n, r_basal)


# ---------- time course: the paper's moment equations ----------

def moment_derivatives(y, r, gamma_m, beta, gamma_p):
    """Exact equations for the means and (co)variances of the linear model
    (no feedback), from the master equation; same as exact_fano in replicate_fig2.py.
    y = (mean m, mean p, var m, cov(m, p), var p)."""
    m, p, vm, c, vp = y
    return (r - gamma_m * m,
            beta * m - gamma_p * p,
            r + gamma_m * m - 2 * gamma_m * vm,
            beta * vm - (gamma_m + gamma_p) * c,
            beta * m + gamma_p * p + 2 * beta * c - 2 * gamma_p * vp)


def moments(r, gamma_m, beta, gamma_p, times, m0=0, p0=0):
    """Thattai et al.'s predicted mean and variance of m and p at each time in
    `times` (evenly spaced, starting at 0), for a cell that starts with exactly
    m0 mRNAs and p0 proteins (so both variances start at 0).

    The equations are linear, y' = M y with y = (m, p, vm, c, vp, 1), so one RK4
    step is a fixed matrix S. Taking ~20 steps per fastest decay time, the update
    over one gap between output times is S^k, computed once and reused.
    Returns four lists: mean_m, mean_p, var_m, var_p.
    """
    if len(times) < 2:
        raise ValueError("need at least 2 time points")
    # With r = 0 the equations are linear in y, so column i of M is the derivative
    # at the i-th unit vector; r only adds constants (the derivative at y = 0).
    M = np.zeros((6, 6))
    for i, unit in enumerate(np.eye(5)):
        M[:5, i] = moment_derivatives(unit, 0.0, gamma_m, beta, gamma_p)
    M[:5, 5] = moment_derivatives((0.0,) * 5, r, gamma_m, beta, gamma_p)

    interval = times[1] - times[0]
    fastest = 2 * max(gamma_m, gamma_p, 1e-12)     # variances decay at up to 2 * gamma
    steps = max(1, math.ceil(interval * fastest * 20))
    hM = M * (interval / steps)
    S = np.eye(6) + hM @ (np.eye(6) + hM @ (np.eye(6) / 2 + hM @ (np.eye(6) / 6 + hM / 24)))
    step = np.linalg.matrix_power(S, steps)

    y = np.array([float(m0), float(p0), 0.0, 0.0, 0.0, 1.0])
    out = np.empty((len(times), 6))
    out[0] = y
    for k in range(1, len(times)):
        y = step @ y
        out[k] = y
    # tiny negative values from rounding at t = 0 would break sqrt
    return (out[:, 0].tolist(), out[:, 1].tolist(),
            np.maximum(out[:, 2], 0.0).tolist(), np.maximum(out[:, 4], 0.0).tolist())
