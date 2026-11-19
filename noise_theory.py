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

Only numpy-free standard library math, so it runs on Vercel without scipy.
"""

import math


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
