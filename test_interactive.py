"""Tests for interactive.py, feeding scripted answers in place of keyboard input."""

import io
import unittest
from contextlib import redirect_stdout

from interactive import run


def scripted(answers):
    """Return an input function that replays answers in order."""
    it = iter(answers)
    return lambda prompt="": next(it)


def run_quietly(answers):
    with redirect_stdout(io.StringIO()):
        return run(input_fn=scripted(answers))


class TestInteractive(unittest.TestCase):
    def test_proposal_repressor_case(self):
        # defaults for r_max/r_basal, one repressor with K=10, n=2, x=10
        r, m, p = run_quietly(["", "", "1", "r", "10", "2", "10", "0.5", "1", "0.5"])
        self.assertAlmostEqual(r, 0.5005)
        self.assertAlmostEqual(m, 0.5005 / 0.5)
        self.assertAlmostEqual(p, 0.5005 / 0.5 / 0.5)

    def test_proposal_steady_state_case(self):
        # r_basal = 0 and activator saturated so r ~ 1; gamma_m=0.5, beta=1, gamma_p=0.5 -> m=2, p=4
        r, m, p = run_quietly(["1", "0", "1", "a", "1", "1", "1e9", "0.5", "1", "0.5"])
        self.assertAlmostEqual(m, 2, places=6)
        self.assertAlmostEqual(p, 4, places=6)

    def test_two_factors(self):
        r, _, _ = run_quietly(["", "", "2", "a", "10", "2", "10", "r", "10", "2", "10", "1", "1", "1"])
        self.assertAlmostEqual(r, 0.001 + 0.999 * 0.25)

    def test_reprompts_on_bad_input(self):
        answers = [
            "", "",
            "zero", "0", "1",          # count: not a number, too small, then 1
            "enhancer", "r",           # kind: invalid, then repressor
            "-5", "0", "10",           # K: negative, zero, then 10
            "abc", "2",                # n
            "-1", "0",                 # x
            "0", "0.5", "1", "0.5",    # gamma_m: zero rejected, then 0.5
        ]
        r, _, _ = run_quietly(answers)
        self.assertAlmostEqual(r, 1.0)

    def test_reprompts_when_basal_not_below_max(self):
        r, _, _ = run_quietly(["1", "2", "1", "0.001", "1", "r", "10", "2", "0", "1", "1", "1"])
        self.assertAlmostEqual(r, 1.0)


if __name__ == "__main__":
    unittest.main()
