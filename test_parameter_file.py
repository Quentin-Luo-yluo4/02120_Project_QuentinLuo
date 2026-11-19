"""Tests for parameter_file.py. Run: python3 -m unittest -v test_parameter_file"""

import os
import tempfile
import unittest

from noise_theory import predict
from parameter_file import parse_parameters, read_parameter_file, uses_feedback

FILES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "parameter_files")


class TestParseParameters(unittest.TestCase):

    def test_all_names(self):
        text = "r = 1\ngamma_m = 0.5\nbeta = 10\ngamma_p = 0.05\nt_max = 300\nK = 200\nn = 2\nr_basal = 0\n"
        self.assertEqual(parse_parameters(text), {
            "r": 1.0, "gamma_m": 0.5, "beta": 10.0, "gamma_p": 0.05, "t_max": 300.0,
            "K": 200.0, "n": 2.0, "r_basal": 0.0})

    def test_spacing_comments_blank_lines_and_windows_line_endings(self):
        text = "# my run\r\n\r\n  r=2.5   # per minute\r\nbeta =1e1\r\n\t\r\n"
        self.assertEqual(parse_parameters(text), {"r": 2.5, "beta": 10.0})

    def test_partial_file(self):
        self.assertEqual(parse_parameters("beta = 20"), {"beta": 20.0})

    def test_feedback_detection(self):
        self.assertIsNone(uses_feedback(parse_parameters("r = 1\nbeta = 2")))  # file doesn't say
        self.assertTrue(uses_feedback(parse_parameters("r = 1\nK = 100")))
        self.assertTrue(uses_feedback(parse_parameters("feedback = 1")))
        self.assertFalse(uses_feedback(parse_parameters("feedback = 0\nK = 100")))  # explicit wins

    def test_errors_name_the_line(self):
        cases = {
            "no equals sign": ("r = 1\nbeta 10", "line 2: expected 'name = value'"),
            "unknown name": ("r = 1\n\ngamma = 3", "line 3: unknown parameter 'gamma'"),
            "wrong case": ("R = 1", "line 1: unknown parameter 'R'"),
            "duplicate": ("r = 1\nr = 2", "line 2: r is given more than once"),
            "not a number": ("beta = ten", "line 1: value of beta is not a number"),
            "missing value": ("beta =", "line 1: value of beta is not a number"),
            "nan": ("r = nan", "line 1: value of r must be a finite number"),
            "infinity": ("t_max = inf", "line 1: value of t_max must be a finite number"),
            "empty file": ("\n# only a comment\n", "no parameters found"),
            "feedback not 0/1": ("feedback = 2", "line 1: feedback must be 0 (off) or 1 (on)"),
        }
        for label, (text, message) in cases.items():
            with self.subTest(label):
                with self.assertRaises(ValueError) as ctx:
                    parse_parameters(text)
                self.assertIn(message, str(ctx.exception))

    def test_read_file_and_example_file(self):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write("r = 3\ngamma_p = 0.1\n")
        try:
            self.assertEqual(read_parameter_file(f.name), {"r": 3.0, "gamma_p": 0.1})
        finally:
            os.remove(f.name)
        example = read_parameter_file(os.path.join(FILES, "example_parameters.txt"))
        self.assertEqual(len(example), 9)
        self.assertTrue(uses_feedback(example))


class TestPaperFiles(unittest.TestCase):
    """Each paper file gives the mean and Fano factor its comments promise."""

    # file -> (feedback on?, mean, Fano) as written in the file's comments
    EXPECTED = {
        "paper_parameters.txt": (False, 1039, 20.35),
        "paper_small_bursts.txt": (False, 1039, 5.84),
        "paper_more_transcription.txt": (False, 2077, 20.35),
        "paper_feedback_parameters.txt": (True, 689, 5.78),
        "paper_no_feedback_same_mean.txt": (False, 689, 10.68),
    }

    def test_every_paper_file_is_listed(self):
        on_disk = {f for f in os.listdir(FILES) if f.startswith("paper_")}
        self.assertEqual(on_disk, set(self.EXPECTED))

    def test_predictions_match_comments(self):
        for name, (feedback, mean, fano) in self.EXPECTED.items():
            with self.subTest(name):
                v = read_parameter_file(os.path.join(FILES, name))
                self.assertIs(uses_feedback(v), feedback)
                autorep = (v["K"], v["n"], v["r_basal"]) if feedback else None
                got_mean, got_fano = predict(v["r"], v["gamma_m"], v["beta"], v["gamma_p"], autorep)
                self.assertAlmostEqual(got_mean, mean, delta=0.5)
                self.assertAlmostEqual(got_fano, fano, delta=0.005)
                self.assertEqual(v["t_max"], 36000)


if __name__ == "__main__":
    unittest.main()
