"""Basic calculation functions for a three-level prokaryotic gene expression model.

Level 1 (TF regulation):
    r(x) = r_basal + (r_max - r_basal) * prod_i h_i(x_i)
    h_i(x_i) = K_i^n_i / (K_i^n_i + x_i^n_i)   repressor
    h_i(x_i) = x_i^n_i / (K_i^n_i + x_i^n_i)   activator

Level 2 (mRNA):     dm/dt = r(x) - gamma_m * m
Level 3 (protein):  dp/dt = beta * m - gamma_p * p
"""

R_MAX = 1.0
R_BASAL = 0.001

REPRESSOR = "repressor"
ACTIVATOR = "activator"


def _check_nonnegative(name, value):
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value}")


def _check_positive(name, value):
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value}")


def _check_hill_inputs(x, K, n):
    _check_nonnegative("x", x)
    _check_positive("K", K)
    _check_positive("n", n)


# ---------- Level 1: transcription factor regulation ----------

def hill_repressor(x, K, n):
    """Fraction of maximal activity left when a repressor is at concentration x."""
    _check_hill_inputs(x, K, n)
    Kn = K ** n
    return Kn / (Kn + x ** n)


def hill_activator(x, K, n):
    """Fraction of maximal activity reached when an activator is at concentration x."""
    _check_hill_inputs(x, K, n)
    xn = x ** n
    return xn / (K ** n + xn)


def hill(x, K, n, kind):
    """Dispatch to the repressor or activator Hill function."""
    if kind == REPRESSOR:
        return hill_repressor(x, K, n)
    if kind == ACTIVATOR:
        return hill_activator(x, K, n)
    raise ValueError(f"kind must be '{REPRESSOR}' or '{ACTIVATOR}', got {kind!r}")


def transcription_rate(factors, r_max=R_MAX, r_basal=R_BASAL):
    """Transcription rate r(x) for one or more transcription factors.

    factors: list of (x, K, n, kind) tuples, one per transcription factor,
             where kind is "repressor" or "activator".
    """
    _check_nonnegative("r_max", r_max)
    _check_nonnegative("r_basal", r_basal)
    if r_basal >= r_max:
        raise ValueError(f"r_basal ({r_basal}) must be less than r_max ({r_max})")
    if len(factors) == 0:
        raise ValueError("at least one transcription factor is required")

    product = 1.0
    for x, K, n, kind in factors:
        product *= hill(x, K, n, kind)
    return r_basal + (r_max - r_basal) * product


# ---------- Level 2: mRNA ----------

def mrna_derivative(m, r, gamma_m):
    """dm/dt = r - gamma_m * m"""
    _check_nonnegative("m", m)
    _check_nonnegative("r", r)
    _check_nonnegative("gamma_m", gamma_m)
    return r - gamma_m * m


def mrna_steady_state(r, gamma_m):
    """m* = r / gamma_m (where dm/dt = 0)."""
    _check_nonnegative("r", r)
    _check_positive("gamma_m", gamma_m)
    return r / gamma_m


# ---------- Level 3: protein ----------

def protein_derivative(p, m, beta, gamma_p):
    """dp/dt = beta * m - gamma_p * p"""
    _check_nonnegative("p", p)
    _check_nonnegative("m", m)
    _check_nonnegative("beta", beta)
    _check_nonnegative("gamma_p", gamma_p)
    return beta * m - gamma_p * p


def protein_steady_state(m, beta, gamma_p):
    """p* = beta * m / gamma_p (where dp/dt = 0)."""
    _check_nonnegative("m", m)
    _check_nonnegative("beta", beta)
    _check_positive("gamma_p", gamma_p)
    return beta * m / gamma_p


if __name__ == "__main__":
    K, n = 10, 2
    for x in (0, K, 1000):
        rep = transcription_rate([(x, K, n, REPRESSOR)])
        act = transcription_rate([(x, K, n, ACTIVATOR)])
        print(f"x = {x:>4}: repressor r = {rep:.4f}, activator r = {act:.4f}")

    r, gamma_m, beta, gamma_p = 1, 0.5, 1, 0.5
    m = mrna_steady_state(r, gamma_m)
    p = protein_steady_state(m, beta, gamma_p)
    print(f"steady state: m = {m}, p = {p}")
    print(f"dm/dt at m*: {mrna_derivative(m, r, gamma_m)}, dp/dt at p*: {protein_derivative(p, m, beta, gamma_p)}")
