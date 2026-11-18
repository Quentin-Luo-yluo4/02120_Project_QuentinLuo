"""Tests for app.py (the web API). Run: python3 -m unittest -v test_app

Three groups:
  TestConnection  - each endpoint returns exactly what calling the model module
                    directly returns, so the API is wired to the real code.
  TestBehaviour   - answers are right, bad input gives a clear 400/422, limits hold.
  TestLiveServer  - real HTTP requests to a running server. By default it starts
                    uvicorn locally; set APP_URL to test a deployed site instead:
                        APP_URL=https://your-app.vercel.app python3 -m unittest -v test_app

Needs: fastapi, httpx, uvicorn (pip install -r requirements.txt httpx uvicorn).
"""

import json
import os
import socket
import subprocess
import sys
import time
import unittest
import urllib.error
import urllib.request

import numpy as np
from fastapi.testclient import TestClient

from app import MAX_TRAJECTORY_POINTS, app
from cascade import cascade_stats
from gene_expression import (
    ACTIVATOR,
    REPRESSOR,
    hill,
    mrna_steady_state,
    protein_steady_state,
    transcription_rate,
)
from gillespie import gillespie, steady_state_stats

client = TestClient(app)

R, GAMMA_M, BETA, GAMMA_P = 1, 0.5, 1, 0.5
FACTORS = [
    {"x": 2.0, "K": 1.0, "n": 2.0, "kind": REPRESSOR},
    {"x": 5.0, "K": 3.0, "n": 1.5, "kind": ACTIVATOR},
]


def as_tuples(factors):
    return [(f["x"], f["K"], f["n"], f["kind"]) for f in factors]


def post_ok(test, path, body):
    response = client.post(path, json=body)
    test.assertEqual(response.status_code, 200, response.text)
    return response.json()


class TestConnection(unittest.TestCase):
    """API output == direct call into the existing modules."""

    def test_health(self):
        response = client.get("/api/health")
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["model_check"], 1.0)  # repressor at x = 0 -> r = r_max
        self.assertEqual(data["kinds"], [REPRESSOR, ACTIVATOR])

    def test_steady_state_matches_gene_expression(self):
        # distinct rates, so a swapped argument (e.g. gamma_m <-> gamma_p) is caught
        gamma_m, beta, gamma_p = 0.5, 2.0, 0.2
        body = {"factors": FACTORS, "r_max": 1.0, "r_basal": 0.01,
                "gamma_m": gamma_m, "beta": beta, "gamma_p": gamma_p}
        data = post_ok(self, "/api/steady-state", body)

        factors = as_tuples(FACTORS)
        r = transcription_rate(factors, r_max=1.0, r_basal=0.01)
        m = mrna_steady_state(r, gamma_m)
        p = protein_steady_state(m, beta, gamma_p)
        self.assertEqual(data["hill"], [hill(*f) for f in factors])
        self.assertEqual(data["r"], r)
        self.assertEqual(data["m"], m)
        self.assertEqual(data["p"], p)

    def test_rate_curve_matches_gene_expression(self):
        body = {"factors": FACTORS, "sweep": 1, "x_min": 0, "x_max": 10, "n_points": 50}
        data = post_ok(self, "/api/rate-curve", body)

        x = np.linspace(0, 10, 50)
        factors = as_tuples(FACTORS)
        factors[1] = (x, 3.0, 1.5, ACTIVATOR)
        np.testing.assert_allclose(data["x"], x)
        np.testing.assert_allclose(data["r"], transcription_rate(factors), rtol=1e-12)

    def test_gillespie_matches_module_with_same_seed(self):
        body = {"r": R, "gamma_m": GAMMA_M, "beta": BETA, "gamma_p": GAMMA_P,
                "t_max": 100, "seed": 7, "n_points": 201}
        data = post_ok(self, "/api/gillespie", body)

        times, ms, ps = gillespie(R, GAMMA_M, BETA, GAMMA_P, 100, seed=7, n_points=201)
        np.testing.assert_allclose(data["t"], times)
        self.assertEqual(data["m"], ms.tolist())
        self.assertEqual(data["p"], ps.tolist())

    def test_gillespie_autorepression_matches_module(self):
        body = {"r": 5, "gamma_m": GAMMA_M, "beta": 2, "gamma_p": GAMMA_P, "t_max": 50,
                "seed": 3, "n_points": 101,
                "autorepression": {"K": 10, "n": 2, "r_basal": 0.1}}
        data = post_ok(self, "/api/gillespie", body)

        _, ms, ps = gillespie(5, GAMMA_M, 2, GAMMA_P, 50, seed=3, n_points=101,
                              autorepression=(10, 2, 0.1))
        self.assertEqual(data["m"], ms.tolist())
        self.assertEqual(data["p"], ps.tolist())

    def test_gillespie_stats_matches_module(self):
        body = {"r": R, "gamma_m": GAMMA_M, "beta": BETA, "gamma_p": GAMMA_P,
                "t_max": 2000, "t_burn": 50, "seed": 1}
        data = post_ok(self, "/api/gillespie-stats", body)

        stats = steady_state_stats(R, GAMMA_M, BETA, GAMMA_P, 2000, t_burn=50, seed=1)
        for key, value in stats.items():
            self.assertAlmostEqual(data["simulated"][key], value, places=12)

    def test_cascade_matches_module(self):
        body = {"k_r1": 0.01, "k_r_max": 0.02, "k_d": 50, "b": [5, 5],
                "mrna_half_life": 120, "protein_half_life": 3600,
                "t_max": 20000, "t_burn": 1000, "seed": 2}
        data = post_ok(self, "/api/cascade", body)

        means, fanos = cascade_stats(0.01, 0.02, 50, [5, 5], 120, 3600, 20000, 1000, seed=2)
        np.testing.assert_allclose(data["mean"], means, rtol=1e-12)
        np.testing.assert_allclose(data["fano"], fanos, rtol=1e-12)


class TestBehaviour(unittest.TestCase):
    """Correct answers, error handling, and limits."""

    # ----- correctness -----

    def test_repressor_curve_limits(self):
        # r(0) = r_max, r(K) = midpoint, r(large x) -> r_basal
        body = {"factors": [{"x": 0, "K": 2, "n": 4, "kind": REPRESSOR}],
                "x_min": 0, "x_max": 2, "n_points": 3, "r_max": 1.0, "r_basal": 0.0}
        data = post_ok(self, "/api/rate-curve", body)
        self.assertAlmostEqual(data["r"][0], 1.0)
        self.assertAlmostEqual(data["r"][2], 0.5)

    def test_log_x_curve(self):
        body = {"factors": FACTORS[:1], "x_min": 0.01, "x_max": 100,
                "n_points": 5, "log_x": True}
        data = post_ok(self, "/api/rate-curve", body)
        np.testing.assert_allclose(data["x"], [0.01, 0.1, 1, 10, 100])

    def test_gillespie_stats_close_to_theory(self):
        body = {"r": R, "gamma_m": GAMMA_M, "beta": BETA, "gamma_p": GAMMA_P,
                "t_max": 20000, "t_burn": 100, "seed": 1}
        data = post_ok(self, "/api/gillespie-stats", body)
        sim, theory = data["simulated"], data["theory"]
        self.assertEqual(theory, {"mean_m": 2.0, "fano_m": 1.0, "mean_p": 4.0, "fano_p": 2.0})
        for key in theory:
            self.assertAlmostEqual(sim[key], theory[key], delta=0.1 * theory[key], msg=key)

    def test_no_theory_with_autorepression(self):
        body = {"r": 5, "gamma_m": GAMMA_M, "beta": 2, "gamma_p": GAMMA_P, "t_max": 100,
                "seed": 1, "autorepression": {"K": 10, "n": 2, "r_basal": 0.1}}
        self.assertNotIn("theory", post_ok(self, "/api/gillespie-stats", body))

    def test_trajectory_counts_are_nonnegative_integers(self):
        body = {"r": R, "gamma_m": GAMMA_M, "beta": BETA, "gamma_p": GAMMA_P,
                "t_max": 200, "seed": 5}
        data = post_ok(self, "/api/gillespie", body)
        self.assertEqual(len(data["t"]), 1001)  # default n_points
        for values in (data["m"], data["p"]):
            self.assertTrue(all(v >= 0 and v == int(v) for v in values))

    def test_nan_becomes_null(self):
        # r = 0 from an empty start: nothing ever happens, mean 0 -> Fano is NaN,
        # which is not valid JSON; the API must send null instead.
        body = {"r": 0, "gamma_m": GAMMA_M, "beta": BETA, "gamma_p": GAMMA_P, "t_max": 10}
        response = client.post("/api/gillespie-stats", json=body)
        self.assertEqual(response.status_code, 200)
        self.assertIsNone(response.json()["simulated"]["fano_m"])
        json.loads(response.text)  # strict JSON parse succeeds

    # ----- invalid values: the model's own checks -> 400 with a message -----

    def assert_400(self, path, body, words):
        response = client.post(path, json=body)
        self.assertEqual(response.status_code, 400, response.text)
        self.assertIn(words, response.json()["error"])

    def test_negative_hill_coefficient(self):
        bad = [{"x": 1, "K": 1, "n": -1, "kind": REPRESSOR}]
        self.assert_400("/api/steady-state",
                        {"factors": bad, "gamma_m": 1, "beta": 1, "gamma_p": 1},
                        "n must be positive")

    def test_basal_above_max(self):
        self.assert_400("/api/steady-state",
                        {"factors": FACTORS, "r_max": 0.1, "r_basal": 0.5,
                         "gamma_m": 1, "beta": 1, "gamma_p": 1},
                        "must be less than r_max")

    def test_zero_degradation_rate(self):
        self.assert_400("/api/steady-state",
                        {"factors": FACTORS, "gamma_m": 0, "beta": 1, "gamma_p": 1},
                        "gamma_m must be positive")

    def test_bad_axis_ranges(self):
        self.assert_400("/api/rate-curve", {"factors": FACTORS, "x_min": 5, "x_max": 1},
                        "x_min (5.0) must be less than x_max")
        self.assert_400("/api/rate-curve",
                        {"factors": FACTORS, "x_min": 0, "x_max": 1, "log_x": True},
                        "positive for a log x-axis")
        self.assert_400("/api/rate-curve", {"factors": FACTORS, "sweep": 2, "x_max": 1},
                        "sweep must be less than")

    def test_gillespie_bad_values(self):
        base = {"r": R, "gamma_m": GAMMA_M, "beta": BETA, "gamma_p": GAMMA_P}
        self.assert_400("/api/gillespie", {**base, "t_max": -5}, "t_max must be positive")
        self.assert_400("/api/gillespie-stats", {**base, "t_max": 10, "t_burn": 20},
                        "t_burn must be in")
        self.assert_400("/api/gillespie",
                        {**base, "r": 1, "t_max": 10,
                         "autorepression": {"K": 1, "n": 1, "r_basal": 2}},
                        "r_basal (2.0) must be less than r")

    def test_cascade_bad_values(self):
        base = {"k_r1": 0.01, "k_r_max": 0.02, "k_d": 50, "b": [5],
                "mrna_half_life": 120, "protein_half_life": 3600, "t_max": 100, "t_burn": 0}
        self.assert_400("/api/cascade", {**base, "k_r1": 0}, "k_r1 must be positive")
        self.assert_400("/api/cascade", {**base, "b": [5, -1]}, "burst sizes")
        self.assert_400("/api/cascade", {**base, "t_burn": 100}, "t_burn must be in")

    # ----- malformed requests: rejected before reaching the model -> 422 -----

    def test_malformed_requests(self):
        cases = {
            "missing field": {"factors": FACTORS, "beta": 1, "gamma_p": 1},
            "unknown kind": {"factors": [{"x": 1, "K": 1, "n": 1, "kind": "enhancer"}],
                             "gamma_m": 1, "beta": 1, "gamma_p": 1},
            "text for a number": {"factors": FACTORS, "gamma_m": "fast",
                                  "beta": 1, "gamma_p": 1},
            "no factors": {"factors": [], "gamma_m": 1, "beta": 1, "gamma_p": 1},
        }
        for label, body in cases.items():
            with self.subTest(label):
                self.assertEqual(client.post("/api/steady-state", json=body).status_code, 422)

    def test_point_cap(self):
        body = {"r": R, "gamma_m": GAMMA_M, "beta": BETA, "gamma_p": GAMMA_P,
                "t_max": 10, "n_points": MAX_TRAJECTORY_POINTS + 1}
        self.assertEqual(client.post("/api/gillespie", json=body).status_code, 422)

    def test_home_redirects_to_docs(self):
        response = client.get("/", follow_redirects=False)
        self.assertEqual(response.status_code, 307)
        self.assertEqual(response.headers["location"], "/docs")
        self.assertEqual(client.get("/docs").status_code, 200)

    def test_unknown_route_and_wrong_method(self):
        self.assertEqual(client.get("/api/nope").status_code, 404)
        self.assertEqual(client.get("/api/steady-state").status_code, 405)

    def test_all_endpoints_listed(self):
        paths = set(client.get("/openapi.json").json()["paths"])
        self.assertEqual(paths, {"/api/health", "/api/steady-state", "/api/rate-curve",
                                 "/api/gillespie", "/api/gillespie-stats", "/api/cascade"})


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class TestLiveServer(unittest.TestCase):
    """Real HTTP over the network: a local uvicorn server, or APP_URL if set."""

    @classmethod
    def setUpClass(cls):
        cls.server = None
        cls.base = os.environ.get("APP_URL", "").rstrip("/")
        if cls.base:
            return
        port = _free_port()
        cls.base = f"http://127.0.0.1:{port}"
        cls.server = subprocess.Popen(
            [sys.executable, "-m", "uvicorn", "app:app", "--port", str(port)],
            cwd=os.path.dirname(os.path.abspath(__file__)),
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        deadline = time.time() + 30
        while time.time() < deadline:
            try:
                urllib.request.urlopen(cls.base + "/api/health", timeout=1)
                return
            except (urllib.error.URLError, ConnectionError):
                time.sleep(0.2)
        cls.server.kill()
        raise RuntimeError("uvicorn did not start within 30 s")

    @classmethod
    def tearDownClass(cls):
        if cls.server is not None:
            cls.server.terminate()
            cls.server.wait(timeout=10)

    def request(self, path, body=None):
        data = None if body is None else json.dumps(body).encode()
        req = urllib.request.Request(self.base + path, data=data,
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=60) as response:
                return response.status, json.loads(response.read())
        except urllib.error.HTTPError as err:
            return err.code, json.loads(err.read())

    def test_health_over_http(self):
        status, data = self.request("/api/health")
        self.assertEqual(status, 200)
        self.assertEqual(data["status"], "ok")

    def test_steady_state_over_http(self):
        status, data = self.request("/api/steady-state", {
            "factors": [{"x": 0, "K": 1, "n": 1, "kind": REPRESSOR}],
            "r_max": 1, "r_basal": 0, "gamma_m": 0.5, "beta": 2, "gamma_p": 0.25})
        self.assertEqual(status, 200)
        self.assertEqual((data["r"], data["m"], data["p"]), (1.0, 2.0, 16.0))

    def test_gillespie_over_http(self):
        # also exercises numba on the server (compile or cache load)
        status, data = self.request("/api/gillespie", {
            "r": 1, "gamma_m": 0.5, "beta": 1, "gamma_p": 0.5, "t_max": 50,
            "seed": 1, "n_points": 51})
        self.assertEqual(status, 200)
        self.assertEqual(len(data["p"]), 51)

    def test_error_over_http(self):
        status, data = self.request("/api/steady-state", {
            "factors": [{"x": 1, "K": -1, "n": 1, "kind": REPRESSOR}],
            "gamma_m": 1, "beta": 1, "gamma_p": 1})
        self.assertEqual(status, 400)
        self.assertIn("K must be positive", data["error"])


if __name__ == "__main__":
    unittest.main()
