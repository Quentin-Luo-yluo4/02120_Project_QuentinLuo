"""Replicate Figs. 1, 3, 4a-b and 5 of Thattai & van Oudenaarden 2001 (PNAS 98:8614).
(Fig. 2 is in replicate_fig2.py.)

For each figure, the paper's parameters come from its caption. Each simulated point
is compared with the paper's formula and counts as a match if it lies within
3 standard errors (SE) of it. SEs come from several independent runs per point.

Paper's notation -> gillespie.py arguments:
    k_R -> r,  gamma_R -> gamma_m,  k_P = b * gamma_R -> beta,  gamma_P -> gamma_p
    half-life h -> decay rate ln(2) / h

Run: python3 replicate_figures.py        (a few minutes; plots go to figures/)
"""

import math
from pathlib import Path

import numpy as np
from scipy.linalg import solve_continuous_lyapunov

from cascade import cascade_stats
from gillespie import gillespie, steady_state_stats

LN2 = math.log(2)
OUT = Path("figures")
RUNS = 8              # independent runs per point (for the SE)
RUN_LENGTH = 6e6      # seconds of simulated time per run, after burn-in
Z = 3                 # a point matches if |simulation - theory| <= Z * SE


# ---------- the paper's formulas ----------

def eq1(k_r, b, gamma_m, gamma_p):
    """Eq. 1, single gene: mean and Fano factor of protein."""
    eta = gamma_p / gamma_m
    return k_r * b / gamma_p, b / (1 + eta) + 1


def hill_rep(p, k_r_max, k_d, n):
    return k_r_max / (1 + (p / k_d) ** n)


def feedback_fixed_point(k_r_max, k_d, n, b, gamma_p):
    """Solve p = k_R(p) * b / gamma_P by bisection (k_R(p) is a Hill repressor)."""
    lo, hi = 0.0, k_r_max * b / gamma_p
    for _ in range(200):
        mid = (lo + hi) / 2
        if hill_rep(mid, k_r_max, k_d, n) * b / gamma_p > mid:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def eq3(k_r_max, k_d, n, b, gamma_m, gamma_p):
    """Eq. 3, autoregulation: linearize k_R(p) ~ k0 - k1 p around the mean."""
    p = feedback_fixed_point(k_r_max, k_d, n, b, gamma_p)
    x = (p / k_d) ** n
    k1 = k_r_max * n * x / (p * (1 + x) ** 2) if p > 0 else 0.0
    k0 = hill_rep(p, k_r_max, k_d, n) + k1 * p
    phi = k1 / gamma_p
    eta = gamma_p / gamma_m
    mean = k0 * b / gamma_p / (1 + b * phi)
    fano = (1 - phi) / (1 + b * phi) * b / (1 + eta) + 1
    return mean, fano, phi


def linear_noise_fano(jacobian, means):
    """Exact steady-state Fano factors of a linear(ized) network (paper's Appendix).

    jacobian: d(rate of change)/d(species); means: steady-state molecule numbers.
    Each species' noise source is 2 * (decay rate) * mean, because creation = decay
    at steady state. Returns variance / mean for every species.
    """
    decay = -np.diag(jacobian)
    noise = np.diag(2 * decay * means)
    cov = solve_continuous_lyapunov(jacobian, -noise)
    return np.diag(cov) / means


# ---------- simulation helper ----------

def simulate(r, gamma_m, beta, gamma_p, autorepression=None, seed0=0):
    """Mean and Fano factor of protein from RUNS independent long runs, with SE."""
    burn = 10 * LN2 / gamma_p  # paper: steady state after 10 protein half-lives
    means, fanos = [], []
    for k in range(RUNS):
        s = steady_state_stats(r, gamma_m, beta, gamma_p, t_max=burn + RUN_LENGTH,
                               t_burn=burn, seed=seed0 + k, autorepression=autorepression)
        means.append(s["mean_p"])
        fanos.append(s["fano_p"])
    se = np.std(fanos, ddof=1) / math.sqrt(RUNS)
    return np.mean(means), np.mean(fanos), se


class Scorecard:
    def __init__(self):
        self.rows = []

    def add(self, figure, label, theory, sim, se):
        ok = abs(sim - theory) <= Z * se
        self.rows.append((figure, label, theory, sim, se, ok))
        print(f"  {label:<28} theory {theory:7.2f}   simulation {sim:7.2f} ± {se:4.2f}   {'match' if ok else 'MISS'}")
        return ok

    def summary(self):
        print("\n=== Summary ===")
        for fig in dict.fromkeys(r[0] for r in self.rows):
            rows = [r for r in self.rows if r[0] == fig]
            print(f"  {fig}: {sum(r[5] for r in rows)}/{len(rows)} points within {Z} SE of the paper's formula")


# ---------- figures ----------

def fig1b(plt):
    """Fig. 1b: protein arrives in bursts (one burst per mRNA)."""
    gm, gp = LN2 / 120, LN2 / 3600
    times, ms, ps = gillespie(0.01, gm, 20 * gm, gp, t_max=3000, seed=4, n_points=3001)
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.plot(times, ps, color="tab:purple", lw=1)
    ax.set_xlabel("t (seconds)")
    ax.set_ylabel("protein p")
    ax.set_title("Fig. 1b: protein is made in bursts (b = 20)")
    fig.tight_layout()
    fig.savefig(OUT / "fig1b_bursts.png", dpi=150)
    plt.close(fig)


def fig1c(plt, card):
    """Fig. 1c: Fano factor vs mean, varying b, k_R or protein half-life."""
    print("\nFig. 1c  (mRNA half-life 2 min; base b = 20, k_R = 0.01/s, protein half-life 1 h)")
    gm = LN2 / 120
    series = {
        "b varied 5-40": [(0.01, b, 3600) for b in range(5, 41, 5)],
        "k_R varied 0.0025-0.02": [(k, 20, 3600) for k in np.arange(0.0025, 0.02001, 0.0025)],
        "protein half-life 15 min-2 h": [(0.01, 20, h * 60) for h in range(15, 121, 15)],
    }
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for (name, points), marker in zip(series.items(), "o^s"):
        print(f" {name}")
        xs, ys, es, txs, tys = [], [], [], [], []
        for i, (k_r, b, hl) in enumerate(points):
            gp = LN2 / hl
            mean_t, fano_t = eq1(k_r, b, gm, gp)
            mean_s, fano_s, se = simulate(k_r, gm, b * gm, gp, seed0=100 * i)
            card.add("Fig. 1c", f"<p>={mean_t:6.0f}", fano_t, fano_s, se)
            xs.append(mean_s), ys.append(fano_s), es.append(se), txs.append(mean_t), tys.append(fano_t)
        ax.plot(txs, tys, color="black", lw=1)
        ax.errorbar(xs, ys, yerr=es, fmt=marker, ms=5, mfc="none", label=name)
    ax.axhline(1, color="gray", ls=":", label="Poisson")
    ax.set_xlabel("mean protein ⟨p⟩")
    ax.set_ylabel("Fano factor δp²/⟨p⟩")
    ax.set_xlim(0, 2200)
    ax.set_ylim(0, 45)
    ax.set_title("Fig. 1c: lines = Eq. 1, symbols = your Gillespie simulation")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig1c_mean_vs_noise.png", dpi=150)
    plt.close(fig)


def fig3b(plt, card):
    """Fig. 3b: negative feedback lowers noise. Vary K_d (n = 2) or n (K_d = 800).

    Expect some MISSes for steep feedback (n >= 6) and K_d = 100: Eq. 3 linearizes the
    Hill curve, which the paper says is exact only "in the limit of small fluctuations".
    The simulation is slightly noisier there (5-15%). Making the system 10x larger
    (k_R_max and K_d both x10, so Eq. 3 is unchanged) shrinks the gap about 10-fold
    (n = 20: +15.4% -> +1.7%), which confirms it's the approximation, not a bug.
    """
    print("\nFig. 3b  (protein half-life 1 h, mRNA half-life 2 min, b = 10, unrepressed mean 1200)")
    gm, gp, b = LN2 / 120, LN2 / 3600, 10
    k_r_max = 1200 * gp / b
    print(f"  -> k_R_max = 1200 * gamma_P / b = {k_r_max:.5f} per s")

    # check that our reading of Eq. 3 agrees with the exact linear-noise solution
    worst = 0.0
    for k_d, n in [(100, 2), (800, 2), (800, 10), (800, 20)]:
        mean, fano, phi = eq3(k_r_max, k_d, n, b, gm, gp)
        k1 = phi * gp
        jac = np.array([[-gm, -k1], [b * gm, -gp]])
        lna = linear_noise_fano(jac, np.array([mean * gp / (b * gm), mean]))[1]
        worst = max(worst, abs(lna - fano))
    print(f"  (Eq. 3 vs exact linear-noise theory: largest difference {worst:.2e})")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    kd_points = list(range(100, 2001, 100)) + [3000, 4000, 5000]
    n_points = list(range(0, 21))
    for name, params, marker in [("K_d varied (n = 2)", [(kd, 2) for kd in kd_points], "^"),
                                 ("n varied (K_d = 800)", [(800, n) for n in n_points], "o")]:
        print(f" {name}")
        xs, ys, es, txs, tys = [], [], [], [], []
        for i, (k_d, n) in enumerate(params):
            if n == 0:   # n = 0: constant rate k_R_max / 2 (no feedback)
                mean_t, fano_t = eq1(k_r_max / 2, b, gm, gp)
                mean_s, fano_s, se = simulate(k_r_max / 2, gm, b * gm, gp, seed0=500 + 100 * i)
            else:
                mean_t, fano_t, _ = eq3(k_r_max, k_d, n, b, gm, gp)
                mean_s, fano_s, se = simulate(k_r_max, gm, b * gm, gp, autorepression=(k_d, n, 0.0),
                                              seed0=500 + 100 * i)
            card.add("Fig. 3b", f"K_d={k_d:<5} n={n:<2} <p>={mean_t:5.0f}", fano_t, fano_s, se)
            xs.append(mean_s), ys.append(fano_s), es.append(se), txs.append(mean_t), tys.append(fano_t)
        ax.plot(txs, tys, color="black", lw=1)
        ax.errorbar(xs, ys, yerr=es, fmt=marker, ms=5, mfc="none", label=name)
    mean_u, fano_u = eq1(k_r_max, b, gm, gp)
    ax.plot([mean_u], [fano_u], "k*", ms=10, label="unrepressed (K_d → ∞)")
    ax.set_xlabel("mean protein ⟨p⟩")
    ax.set_ylabel("Fano factor δp²/⟨p⟩")
    ax.set_xlim(0, 1400)
    ax.set_ylim(0, 12)
    ax.set_title("Fig. 3b: lines = Eq. 3, symbols = your Gillespie simulation")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig3b_feedback.png", dpi=150)
    plt.close(fig)


def fig3a_4ab(plt):
    """Figs. 3a / 4a-b: protein histograms without and with feedback, with the
    analytic Gaussian (mean and variance from Eq. 1 / Eq. 3) on top."""
    print("\nFigs. 3a / 4a-b  (histograms; feedback K_d = 800, n = 2 — the caption doesn't say, see notes)")
    gm, gp, b = LN2 / 120, LN2 / 3600, 10
    k_r_max = 1200 * gp / b
    fig, ax = plt.subplots(figsize=(7, 4))
    p_axis = np.linspace(0, 1800, 400)
    for label, autorep, color, theory in [
        ("unrepressed", None, "tab:gray", eq1(k_r_max, b, gm, gp)),
        ("repressed (K_d = 800, n = 2)", (800, 2, 0.0), "tab:green", eq3(k_r_max, 800, 2, b, gm, gp)[:2]),
    ]:
        _, _, ps = gillespie(k_r_max, gm, b * gm, gp, t_max=5e7, seed=11, autorepression=autorep, n_points=10001)
        ps = ps[100:]  # drop the first 100 snapshots (burn-in)
        ax.hist(ps, bins=np.arange(0, 1801, 25), density=True, alpha=0.45, color=color, label=f"{label}: simulation")
        mean, fano = theory
        sd = math.sqrt(fano * mean)
        ax.plot(p_axis, np.exp(-(p_axis - mean) ** 2 / (2 * sd ** 2)) / (sd * math.sqrt(2 * math.pi)),
                color=color, lw=2, label=f"{label}: Gaussian from theory")
        print(f"  {label:<30} theory mean {mean:6.0f} sd {sd:5.1f}   simulation mean {ps.mean():6.0f} sd {ps.std():5.1f}")
    ax.set_xlabel("protein p")
    ax.set_ylabel("fraction of time")
    ax.set_title("Figs. 3a / 4a-b: feedback lowers and narrows the distribution")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig3a_4ab_histograms.png", dpi=150)
    plt.close(fig)


def cascade_theory(b_list, k_d=500):
    """Fig. 5 analytic result: linear-noise solution of the cascade."""
    gm, gp = LN2 / 120, LN2 / 3600
    n = len(b_list)
    k_r_max = 2 * k_d * gp / b_list[1] if n > 1 else 0
    k_r1 = k_d * gp / b_list[0]
    slope = k_r_max / (4 * k_d)   # d k_R / d p at p = K_d for the hyperbolic activator
    jac = np.zeros((2 * n, 2 * n))
    means = np.zeros(2 * n)
    for i in range(n):
        mi, pi = 2 * i, 2 * i + 1
        jac[mi, mi], jac[pi, pi], jac[pi, mi] = -gm, -gp, b_list[i] * gm
        if i > 0:
            jac[mi, 2 * (i - 1) + 1] = slope
        k_r = k_r1 if i == 0 else k_r_max / 2
        means[mi], means[pi] = k_r / gm, k_d
    return linear_noise_fano(jac, means)[1::2], k_r1, k_r_max


def fig5(plt, card):
    """Fig. 5: noise propagates down a 6-gene cascade (all means fixed at K_d = 500)."""
    print("\nFig. 5  (6 genes; protein half-life 1 h, mRNA half-life 2 min, K_d = 500, downstream b = 20)")
    k_d, runs, length = 500, 4, 2e7
    burn = 10 * 3600

    def sim(b_list, seed0):
        _, k_r1, k_r_max = cascade_theory(b_list, k_d)
        fanos = np.array([cascade_stats(k_r1, k_r_max, k_d, b_list, 120, 3600, burn + length, burn,
                                        seed=seed0 + k)[1] for k in range(runs)])
        return fanos.mean(0), fanos.std(0, ddof=1) / math.sqrt(runs)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10, 4))
    print(" Fig. 5a: regulator burst size b1 varied")
    b1_values = [1, 2, 5, 10, 15, 20]
    theory_a, sim_a, se_a = [], [], []
    for j, b1 in enumerate(b1_values):
        b_list = [b1] + [20] * 5
        th, _, _ = cascade_theory(b_list, k_d)
        sf, se = sim(b_list, 9000 + 10 * j)
        for g in (1, 2, 3):
            card.add("Fig. 5", f"b1={b1:<2} gene P{g + 1}", th[g], sf[g], se[g])
        theory_a.append(th), sim_a.append(sf), se_a.append(se)
    theory_a, sim_a, se_a = map(np.array, (theory_a, sim_a, se_a))
    dense_b1 = np.linspace(1, 20, 40)
    dense = np.array([cascade_theory([b1] + [20] * 5, k_d)[0] for b1 in dense_b1])
    for g, color in zip((1, 2, 3), ("tab:blue", "tab:orange", "tab:green")):
        ax_a.plot(dense[:, 0], dense[:, g], color=color, lw=1)
        ax_a.errorbar(sim_a[:, 0], sim_a[:, g], yerr=se_a[:, g], fmt="o", ms=4, color=color, label=f"P{g + 1}")
    ax_a.axhline(20 / (1 + (LN2 / 3600) / (LN2 / 120)) + 1, color="gray", ls=":", label="no input noise (Eq. 1)")
    ax_a.set_xlabel("Fano factor of regulator P1")
    ax_a.set_ylabel("Fano factor downstream")
    ax_a.set_title("Fig. 5a")
    ax_a.legend(fontsize=8, frameon=False)

    print(" Fig. 5b: b = 20 for all six genes")
    b_list = [20] * 6
    th, _, _ = cascade_theory(b_list, k_d)
    sf, se = sim(b_list, 9900)
    for g in range(6):
        card.add("Fig. 5", f"all b=20, gene P{g + 1}", th[g], sf[g], se[g])
    idx = np.arange(1, 7)
    ax_b.plot(idx, th, "k-", lw=1, label="analytic")
    ax_b.errorbar(idx, sf, yerr=se, fmt="o", ms=5, label="Gillespie")
    ax_b.set_xlabel("protein index i")
    ax_b.set_ylabel("Fano factor of P_i")
    ax_b.set_title("Fig. 5b")
    ax_b.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_cascade.png", dpi=150)
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUT.mkdir(exist_ok=True)
    card = Scorecard()
    fig1b(plt)
    fig1c(plt, card)
    fig3b(plt, card)
    fig3a_4ab(plt)
    fig5(plt, card)
    card.summary()
    print(f"\nplots saved in {OUT}/")


if __name__ == "__main__":
    main()
