"""Replicate Fig. 2 of Thattai & van Oudenaarden 2001 (PNAS 98:8614): transient noise
after a single gene is switched on.

Paper's parameters (Fig. 2 caption):
    protein half-life 1 h, burst size b = 20, transcription rate k_R = 0.01 per s,
    mRNA half-life 120 s, 24 s and 0 s  (eta = gamma_P / gamma_R = 1/30, 1/150, 0).
The gene is switched on at t = 0 with no mRNA and no protein, and the Fano factor
(variance / mean of protein number across cells) is followed for 5000 s.

Paper's notation -> gillespie.py arguments:
    k_R -> r,  gamma_R -> gamma_m,  k_P = b * gamma_R -> beta,  gamma_P -> gamma_p

What this script compares:
    simulation : Fano factor across N_CELLS independent runs of gillespie()
    exact      : the model's mean and variance solved exactly over time
                 (the paper's "exact analytic results" for finite eta)
    Eq. 2      : the paper's formula for eta -> 0
    Eq. 1      : the steady-state value b / (1 + eta) + 1

Run: python3 replicate_fig2.py      (about a minute; saves fig2_replication.png)
"""

import math

import numpy as np
from scipy.integrate import solve_ivp

from gillespie import gillespie

K_R = 0.01                      # transcripts per second
B = 20                          # proteins per transcript
GAMMA_P = math.log(2) / 3600    # protein half-life 1 h
MRNA_HALF_LIVES = [120, 24]     # seconds; eta = 1/30 and 1/150 (0 s can't be simulated)
T_MAX = 5000
N_POINTS = 101                  # one snapshot every 50 s
N_CELLS = 10000                 # the paper used 10,000 trials per data point
N_GROUPS = 10                   # split cells into groups to get error bars
SEED = 2001


def rates(mrna_half_life):
    gamma_m = math.log(2) / mrna_half_life
    return K_R, gamma_m, B * gamma_m, GAMMA_P


def exact_fano(times, r, gamma_m, beta, gamma_p):
    """Solve the moment equations of the linear model exactly (from zero molecules)."""
    def deriv(t, y):
        m, p, vm, c, vp = y
        return [r - gamma_m * m,
                beta * m - gamma_p * p,
                r + gamma_m * m - 2 * gamma_m * vm,
                beta * vm - (gamma_m + gamma_p) * c,
                beta * m + gamma_p * p + 2 * beta * c - 2 * gamma_p * vp]
    sol = solve_ivp(deriv, (0, times[-1]), [0, 0, 0, 0, 0], t_eval=times,
                    rtol=1e-10, atol=1e-12, method="LSODA")
    p, vp = sol.y[1], sol.y[4]
    return np.divide(vp, p, out=np.full_like(p, np.nan), where=p > 1e-9)


def eq2_fano(times):
    """Paper's Eq. 2 (eta -> 0)."""
    x = np.exp(-GAMMA_P * times)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (1 - x ** 2) / (1 - x) * B + 1


def simulated_fano(r, gamma_m, beta, gamma_p):
    """Fano factor across cells at each snapshot time, with a standard error."""
    proteins = np.empty((N_CELLS, N_POINTS))
    for i in range(N_CELLS):
        times, _, ps = gillespie(r, gamma_m, beta, gamma_p, T_MAX, seed=SEED + i, n_points=N_POINTS)
        proteins[i] = ps
    mean = proteins.mean(axis=0)
    fano = np.divide(proteins.var(axis=0, ddof=1), mean, out=np.full(N_POINTS, np.nan), where=mean > 0)
    groups = proteins.reshape(N_GROUPS, -1, N_POINTS)
    g_mean = groups.mean(axis=1)
    g_fano = np.divide(groups.var(axis=1, ddof=1), g_mean,
                       out=np.full_like(g_mean, np.nan), where=g_mean > 0)
    se = g_fano.std(axis=0, ddof=1) / math.sqrt(N_GROUPS)
    return times, fano, se


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    dense = np.linspace(0, T_MAX, 1001)
    ax.plot(dense[1:], eq2_fano(dense[1:]), color="black", lw=1.2, label="η = 0 (paper Eq. 2)")

    colors = {120: "tab:blue", 24: "tab:red"}
    for half_life in MRNA_HALF_LIVES:
        r, gamma_m, beta, gamma_p = rates(half_life)
        eta = gamma_p / gamma_m
        steady = B / (1 + eta) + 1
        times, fano, se = simulated_fano(r, gamma_m, beta, gamma_p)
        exact = exact_fano(times, r, gamma_m, beta, gamma_p)
        exact_dense = exact_fano(dense, r, gamma_m, beta, gamma_p)

        k_sim, k_ex = np.nanargmax(fano), np.nanargmax(exact_dense)
        ok = np.abs(fano[1:] - exact[1:]) <= 3 * se[1:]
        print(f"\nmRNA half-life {half_life} s  (eta = 1/{round(1 / eta)})")
        print(f"  steady state (Eq. 1):      {steady:6.2f}")
        print(f"  peak, exact theory:        {exact_dense[k_ex]:6.2f} at t = {dense[k_ex]:.0f} s"
              f"  -> {exact_dense[k_ex] / steady:.2f} x steady state")
        print(f"  peak, simulation:          {fano[k_sim]:6.2f} ± {se[k_sim]:.2f} at t = {times[k_sim]:.0f} s"
              f"  -> {fano[k_sim] / steady:.2f} x steady state")
        print(f"  simulation within 3 SE of exact theory at {ok.sum()}/{ok.size} time points")
        for t in (100, 500, 1000, 2000, 3000, 5000):
            j = int(round(t / (T_MAX / (N_POINTS - 1))))
            print(f"    t = {t:>4} s: simulated {fano[j]:6.2f} ± {se[j]:.2f}   exact {exact[j]:6.2f}")

        ax.plot(dense[1:], exact_dense[1:], color=colors[half_life], lw=1.2,
                label=f"η = 1/{round(1 / eta)} exact")
        ax.errorbar(times[1::2], fano[1::2], yerr=se[1::2], fmt="o", ms=3.5,
                    color=colors[half_life], label=f"η = 1/{round(1 / eta)} Gillespie ({N_CELLS:,} cells)")

    ax.axhline(B + 1, color="gray", ls=":", lw=1, label="steady state (Eq. 1, η = 0)")
    ax.set_xlabel("t (seconds)")
    ax.set_ylabel("Fano factor  δp² / ⟨p⟩")
    ax.set_xlim(0, T_MAX)
    ax.set_ylim(0, 45)
    ax.set_title("Replication of Thattai & van Oudenaarden 2001, Fig. 2")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig("fig2_replication.png", dpi=150)
    print("\nsaved fig2_replication.png")


if __name__ == "__main__":
    main()
