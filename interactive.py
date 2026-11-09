"""Interactive command-line front end for gene_expression.py. Run: python3 interactive.py

The user enters one or more transcription factors (x, K, n, repressor/activator)
and the mRNA/protein rates; the program prints r(x) and the steady-state m and p.
Pressing Enter at a prompt that shows [default] accepts the default.
"""

from gene_expression import (
    ACTIVATOR,
    R_BASAL,
    R_MAX,
    REPRESSOR,
    hill,
    mrna_steady_state,
    protein_steady_state,
    transcription_rate,
)


def ask_float(prompt, default=None, check=None, input_fn=input):
    """Ask until the user types a number that passes check (a function raising ValueError)."""
    suffix = f" [{default}]" if default is not None else ""
    while True:
        text = input_fn(f"{prompt}{suffix}: ").strip()
        if text == "" and default is not None:
            return default
        try:
            value = float(text)
            if check is not None:
                check(value)
            return value
        except ValueError as err:
            print(f"  Invalid input: {err}")


def ask_kind(prompt, input_fn=input):
    """Ask until the user types repressor/activator (or r/a)."""
    while True:
        text = input_fn(f"{prompt} (r = repressor, a = activator): ").strip().lower()
        if text in ("r", REPRESSOR):
            return REPRESSOR
        if text in ("a", ACTIVATOR):
            return ACTIVATOR
        print("  Invalid input: type 'r' or 'a'")


def ask_int(prompt, minimum, input_fn=input):
    while True:
        text = input_fn(f"{prompt}: ").strip()
        try:
            value = int(text)
        except ValueError:
            print("  Invalid input: enter a whole number")
            continue
        if value < minimum:
            print(f"  Invalid input: must be at least {minimum}")
            continue
        return value


def _nonnegative(value):
    if value < 0:
        raise ValueError("must be non-negative")


def _positive(value):
    if value <= 0:
        raise ValueError("must be positive")


def ask_factor(i, input_fn=input):
    """Ask for one transcription factor; returns an (x, K, n, kind) tuple."""
    print(f"\nTranscription factor {i}")
    kind = ask_kind("  Type", input_fn=input_fn)
    K = ask_float("  K (half-maximal concentration)", check=_positive, input_fn=input_fn)
    n = ask_float("  n (Hill coefficient)", check=_positive, input_fn=input_fn)
    x = ask_float("  x (concentration)", check=_nonnegative, input_fn=input_fn)
    hill(x, K, n, kind)  # final check with the model's own validation
    return (x, K, n, kind)


def ask_rate_limits(input_fn=input):
    """Ask for r_max and r_basal, re-asking both until r_basal < r_max."""
    while True:
        r_max = ask_float("r_max", default=R_MAX, check=_nonnegative, input_fn=input_fn)
        r_basal = ask_float("r_basal", default=R_BASAL, check=_nonnegative, input_fn=input_fn)
        if r_basal < r_max:
            return r_max, r_basal
        print(f"  Invalid input: r_basal ({r_basal}) must be less than r_max ({r_max})")


def run(input_fn=input):
    print("=== Gene expression simulator ===")
    print("\nLevel 1: transcription regulation")
    r_max, r_basal = ask_rate_limits(input_fn=input_fn)
    count = ask_int("How many transcription factors", minimum=1, input_fn=input_fn)
    factors = [ask_factor(i + 1, input_fn=input_fn) for i in range(count)]
    r = transcription_rate(factors, r_max=r_max, r_basal=r_basal)

    print("\nLevels 2-3: mRNA and protein")
    gamma_m = ask_float("gamma_m (mRNA degradation rate)", check=_positive, input_fn=input_fn)
    beta = ask_float("beta (translation rate)", check=_nonnegative, input_fn=input_fn)
    gamma_p = ask_float("gamma_p (protein degradation rate)", check=_positive, input_fn=input_fn)
    m = mrna_steady_state(r, gamma_m)
    p = protein_steady_state(m, beta, gamma_p)

    print("\n=== Results ===")
    for i, (x, K, n, kind) in enumerate(factors, start=1):
        print(f"TF {i}: {kind:<9} x = {x:g}, K = {K:g}, n = {n:g}, h = {hill(x, K, n, kind):.4f}")
    print(f"Transcription rate r(x) = {r:.4f}")
    print(f"Steady-state mRNA    m* = {m:.4f}")
    print(f"Steady-state protein p* = {p:.4f}")
    return r, m, p


if __name__ == "__main__":
    try:
        run()
    except (KeyboardInterrupt, EOFError):
        print("\nExiting.")
