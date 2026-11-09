"""Try your own parameter values: does the Gillespie simulation match the paper's formula?

Give the values in the paper's notation (Thattai & van Oudenaarden 2001):

    --kR     transcription rate, transcripts per second           (paper: k_R)
    --b      proteins made per mRNA, the burst size               (paper: b = k_P / gamma_R)
    --mrna-half-life     mRNA half-life in seconds                (gamma_R = ln2 / half-life)
    --protein-half-life  protein half-life in seconds             (gamma_P = ln2 / half-life)

Optional:
    --Kd K --n N   add negative feedback: k_R becomes kR / (1 + (p/K)^N)   (paper Eq. 3)
    --induction    also check the noise after switching the gene on       (paper Eq. 2 / Fig. 2)
    --runs, --length, --seed   how much simulation to do

Examples:
    python3 try_parameters.py --kR 0.01 --b 20 --mrna-half-life 120 --protein-half-life 3600
    python3 try_parameters.py --kR 0.0231 --b 10 --mrna-half-life 120 --protein-half-life 3600 --Kd 800 --n 2
    python3 try_parameters.py --kR 0.01 --b 20 --mrna-half-life 24 --protein-half-life 3600 --induction
"""

import argparse
import math

import numpy as np

from gillespie import gillespie, steady_state_stats
from replicate_fig2 import exact_fano

LN2 = math.log(2)
Z = 3  # a result matches if it is within 3 standard errors of the formula


def theory(k_r, b, gamma_m, gamma_p, k_d=None, n=None):
    """The paper's prediction: Eq. 1 without feedback, Eq. 3 with it."""
    eta = gamma_p / gamma_m
    if k_d is None:
        return k_r * b / gamma_p, b / (1 + eta) + 1, "Eq. 1"
    # mean: solve p = k_R(p) * b / gamma_P (bisection)
    lo, hi = 0.0, k_r * b / gamma_p
    for _ in range(200):
        mid = (lo + hi) / 2
        if k_r / (1 + (mid / k_d) ** n) * b / gamma_p > mid:
            lo = mid
        else:
            hi = mid
    p = (lo + hi) / 2
    x = (p / k_d) ** n
    phi = k_r * n * x / (p * (1 + x) ** 2) / gamma_p   # feedback strength k1 / gamma_P
    return p, (1 - phi) / (1 + b * phi) * b / (1 + eta) + 1, "Eq. 3"


def verdict(name, predicted, simulated, se):
    ok = abs(simulated - predicted) <= Z * se
    diff = 100 * (simulated - predicted) / predicted
    print(f"  {name:<14} formula {predicted:9.2f}   simulation {simulated:9.2f} ± {se:.2f}"
          f"   ({diff:+.1f}%)   {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kR", type=float, required=True)
    ap.add_argument("--b", type=float, required=True)
    ap.add_argument("--mrna-half-life", type=float, required=True)
    ap.add_argument("--protein-half-life", type=float, required=True)
    ap.add_argument("--Kd", type=float)
    ap.add_argument("--n", type=float)
    ap.add_argument("--induction", action="store_true")
    ap.add_argument("--runs", type=int, default=8)
    ap.add_argument("--length", type=float, default=None,
                    help="simulated seconds per run (default: 2000 protein half-lives)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if (a.Kd is None) != (a.n is None):
        ap.error("give both --Kd and --n for feedback, or neither")
    for name, v in [("kR", a.kR), ("b", a.b), ("mrna-half-life", a.mrna_half_life),
                    ("protein-half-life", a.protein_half_life)]:
        if v <= 0:
            ap.error(f"--{name} must be positive")
    if a.runs < 2:
        ap.error("--runs must be at least 2 to estimate an error")

    # 1. convert the paper's notation to gillespie.py's arguments
    gamma_m = LN2 / a.mrna_half_life
    gamma_p = LN2 / a.protein_half_life
    beta = a.b * gamma_m
    autorep = None if a.Kd is None else (a.Kd, a.n, 0.0)
    print("Your values, converted to gillespie.py arguments:")
    print(f"  r = kR = {a.kR:g}   gamma_m = ln2/{a.mrna_half_life:g} = {gamma_m:.5g}"
          f"   beta = b*gamma_m = {beta:.5g}   gamma_p = ln2/{a.protein_half_life:g} = {gamma_p:.5g}")
    if autorep:
        print(f"  autorepression = (K={a.Kd:g}, n={a.n:g}, r_basal=0)")

    # 2. what the paper's formula predicts
    mean_t, fano_t, eq = theory(a.kR, a.b, gamma_m, gamma_p, a.Kd, a.n)
    print(f"\nPaper's prediction ({eq}): mean protein {mean_t:.1f}, Fano factor {fano_t:.2f}")

    # warn if the run would be very slow (about one event per protein made or destroyed)
    length = a.length or 2000 * a.protein_half_life
    events = a.runs * length * (2 * mean_t * gamma_p + 2 * a.kR)
    if events > 3e9:
        print(f"  (about {events:.1e} reactions to simulate: this may take minutes; lower --length or --runs)")

    # 3. simulate: several independent long runs, then average
    burn = 10 * a.protein_half_life  # the paper: steady state after 10 protein half-lives
    means, fanos = [], []
    for k in range(a.runs):
        s = steady_state_stats(a.kR, gamma_m, beta, gamma_p, t_max=burn + length, t_burn=burn,
                               seed=a.seed + k, autorepression=autorep)
        means.append(s["mean_p"])
        fanos.append(s["fano_p"])
    se = lambda xs: np.std(xs, ddof=1) / math.sqrt(len(xs))
    print(f"\nSimulation: {a.runs} runs x {length:.3g} s each. Steady-state check (pass = within {Z} SE):")
    results = [verdict("mean protein", mean_t, np.mean(means), se(means)),
               verdict("Fano factor", fano_t, np.mean(fanos), se(fanos))]

    # 4. optional: noise right after switching the gene on (no feedback only)
    if a.induction:
        if autorep:
            print("\n--induction compares with the no-feedback theory, so it can't be used with --Kd/--n.")
        else:
            cells, t_end = 4000, 5 * a.protein_half_life
            times = np.linspace(0, t_end, 101)
            ps = np.array([gillespie(a.kR, gamma_m, beta, gamma_p, t_end, seed=a.seed + i, n_points=101)[2]
                           for i in range(cells)])
            mean = ps.mean(0)
            fano = np.divide(ps.var(0, ddof=1), mean, out=np.full(101, np.nan), where=mean > 0)
            g = ps.reshape(8, -1, 101)
            gm_ = g.mean(1)
            g_f = np.divide(g.var(1, ddof=1), gm_, out=np.full_like(gm_, np.nan), where=gm_ > 0)
            f_se = g_f.std(0, ddof=1) / math.sqrt(8)
            exact = exact_fano(times, a.kR, gamma_m, beta, gamma_p)
            k = int(np.nanargmax(exact))
            print(f"\nInduction check ({cells} cells, from zero molecules):")
            print(f"  exact theory peaks at t = {times[k]:.0f} s, {exact[k] / fano_t:.2f} x the steady state")
            results.append(verdict("peak Fano", exact[k], fano[k], f_se[k]))
            inside = np.abs(fano[1:] - exact[1:]) <= Z * f_se[1:]
            print(f"  whole curve: {inside.sum()}/{inside.size} time points within {Z} SE")

    if all(results):
        print("\nOverall: everything matches the paper's formula.")
    else:
        print("\nOverall: something doesn't match. Things to check:\n"
              "  - rerun with a different --seed: about 1 check in 100 fails by chance at 3 SE\n"
              "  - run longer (--length): too-short runs underestimate the noise\n"
              "  - with feedback, Eq. 3 is an approximation; it gets worse for very steep (large n)\n"
              "    or very strong feedback, where the formula itself may be off")


if __name__ == "__main__":
    main()
