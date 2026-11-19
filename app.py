"""Web API for the gene expression model. Every calculation is done by the existing
modules (gene_expression.py, gillespie.py, cascade.py); this file only turns web
requests into function calls and the results into JSON.

Run locally:   uvicorn app:app --reload      then open http://127.0.0.1:8000/docs
Deploy:        Vercel detects the FastAPI `app` object automatically.

Endpoints (all return JSON):
    GET  /api/health             is the server up, and can it reach the model code?
    POST /api/steady-state       Levels 1-3: h_i, r(x), m*, p*   (same as interactive.py)
    POST /api/rate-curve         r(x) while one TF's x is swept   (for plotting)
    POST /api/gillespie          one stochastic trajectory m(t), p(t)
    POST /api/time-course        trajectory + predicted m(t), p(t) and spread (for the web page)
    POST /api/parse-parameters   read "name = value" lines from an uploaded text file
    POST /api/noise-check        measured protein mean/Fano (several long runs) vs the paper
    POST /api/gillespie-stats    time-weighted mean / variance / Fano of m and p
    POST /api/cascade            mean and Fano of each protein in a gene cascade

Bad values (e.g. n <= 0) are rejected by the model's own checks and come back as
HTTP 400 with {"error": "<message>"}.
"""

import math
import os

# Vercel's file system is read-only except /tmp, and numba's cache=True writes a
# cache file next to the code. Point it at /tmp before numba is imported.
if os.environ.get("VERCEL"):
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")

from typing import List, Literal, Optional  # noqa: E402  (must follow the env setup)

import numpy as np  # noqa: E402
from fastapi import FastAPI, Request  # noqa: E402
from fastapi.responses import HTMLResponse, JSONResponse  # noqa: E402
from pydantic import BaseModel, Field  # noqa: E402

from cascade import cascade_stats  # noqa: E402
from gene_expression import (  # noqa: E402
    ACTIVATOR,
    R_BASAL,
    R_MAX,
    REPRESSOR,
    hill,
    mrna_derivative,
    mrna_steady_state,
    protein_derivative,
    protein_steady_state,
    transcription_rate,
)
from gillespie import gillespie, steady_state_stats  # noqa: E402
from noise_theory import moments, predict  # noqa: E402
from parameter_file import parse_parameters, uses_feedback  # noqa: E402

# Caps that keep one request inside Vercel's time and 4.5 MB response limits.
MAX_CURVE_POINTS = 2000
MAX_TRAJECTORY_POINTS = 5000
MAX_CASCADE_GENES = 10
MAX_ODE_STEPS = 50_000
MAX_PARAMETER_FILE_CHARS = 10_000

# /api/noise-check: like replicate_figures.py, several independent long runs give a
# mean and a standard error (SE); a result "matches" if within MATCH_SE of the paper.
NOISE_RUNS = 8
NOISE_BURN_LIFETIMES = 10      # paper: steady state after ~10 protein half-lives
NOISE_MAX_LIFETIMES = 2000     # length of each run, in protein lifetimes 1/gamma_p ...
NOISE_EVENT_BUDGET = 3e7       # ... shortened so all runs together stay under ~0.5 s
NOISE_MIN_LIFETIMES = 20       # below this the estimate is meaningless; refuse
MATCH_SE = 3

PAGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static", "index.html")

app = FastAPI(title="Gene expression simulator")


@app.exception_handler(ValueError)
async def model_error(request: Request, err: ValueError):
    """The model raises ValueError for invalid parameters; report it as a 400."""
    return JSONResponse(status_code=400, content={"error": str(err)})


def _num(value):
    """Plain float for JSON; NaN/inf (not valid JSON) become null."""
    value = float(value)
    return value if math.isfinite(value) else None


def _nums(array):
    return [_num(v) for v in np.asarray(array, dtype=float)]


# ---------- request bodies ----------

class Factor(BaseModel):
    x: float
    K: float
    n: float
    kind: Literal["repressor", "activator"]

    def as_tuple(self):
        return (self.x, self.K, self.n, self.kind)


class SteadyStateRequest(BaseModel):
    factors: List[Factor] = Field(min_length=1)
    r_max: float = R_MAX
    r_basal: float = R_BASAL
    gamma_m: float
    beta: float
    gamma_p: float


class RateCurveRequest(BaseModel):
    factors: List[Factor] = Field(min_length=1)
    sweep: int = Field(0, ge=0, description="index of the factor whose x is swept")
    x_min: float = 0.0
    x_max: float
    n_points: int = Field(500, ge=2, le=MAX_CURVE_POINTS)
    log_x: bool = False
    r_max: float = R_MAX
    r_basal: float = R_BASAL


class Autorepression(BaseModel):
    K: float
    n: float
    r_basal: float


class GillespieRequest(BaseModel):
    r: float
    gamma_m: float
    beta: float
    gamma_p: float
    t_max: float
    m0: int = 0
    p0: int = 0
    seed: Optional[int] = None
    autorepression: Optional[Autorepression] = None
    n_points: int = Field(1001, ge=2, le=MAX_TRAJECTORY_POINTS)


class GillespieStatsRequest(BaseModel):
    r: float
    gamma_m: float
    beta: float
    gamma_p: float
    t_max: float
    t_burn: float = 0.0
    m0: int = 0
    p0: int = 0
    seed: Optional[int] = None
    autorepression: Optional[Autorepression] = None


class CascadeRequest(BaseModel):
    k_r1: float
    k_r_max: float
    k_d: float
    b: List[float] = Field(min_length=1, max_length=MAX_CASCADE_GENES)
    mrna_half_life: float
    protein_half_life: float
    t_max: float
    t_burn: float
    seed: int = 0


class NoiseCheckRequest(BaseModel):
    r: float
    gamma_m: float
    beta: float
    gamma_p: float
    autorepression: Optional[Autorepression] = None
    seed: int = 0


class ParameterFileRequest(BaseModel):
    text: str = Field(max_length=MAX_PARAMETER_FILE_CHARS)


def _autorepression_tuple(a):
    return None if a is None else (a.K, a.n, a.r_basal)


def _theory(r, gamma_m, beta, gamma_p):
    """Analytic steady state for constant r (see gillespie.py docstring)."""
    return {
        "mean_m": _num(r / gamma_m),
        "fano_m": 1.0,
        "mean_p": _num(r * beta / (gamma_m * gamma_p)),
        "fano_p": _num(1 + beta / (gamma_m + gamma_p)),
    }


# ---------- endpoints ----------

@app.get("/", include_in_schema=False)
def home():
    """The interactive graphing page (static/index.html)."""
    with open(PAGE, encoding="utf-8") as f:
        return HTMLResponse(f.read())


@app.get("/api/health")
def health():
    """Checks that the model modules are importable and give a known answer."""
    r = transcription_rate([(0.0, 1.0, 1.0, REPRESSOR)], r_max=1.0, r_basal=0.0)
    return {"status": "ok", "model_check": _num(r), "kinds": [REPRESSOR, ACTIVATOR]}


@app.post("/api/steady-state")
def steady_state(req: SteadyStateRequest):
    factors = [f.as_tuple() for f in req.factors]
    r = transcription_rate(factors, r_max=req.r_max, r_basal=req.r_basal)
    m = mrna_steady_state(r, req.gamma_m)
    p = protein_steady_state(m, req.beta, req.gamma_p)
    return {
        "hill": [_num(hill(*f)) for f in factors],
        "r": _num(r),
        "m": _num(m),
        "p": _num(p),
    }


@app.post("/api/rate-curve")
def rate_curve(req: RateCurveRequest):
    if req.sweep >= len(req.factors):
        raise ValueError(f"sweep must be less than the number of factors ({len(req.factors)})")
    if req.x_min >= req.x_max:
        raise ValueError(f"x_min ({req.x_min}) must be less than x_max ({req.x_max})")
    if req.log_x:
        if req.x_min <= 0:
            raise ValueError("x_min must be positive for a log x-axis")
        x = np.logspace(np.log10(req.x_min), np.log10(req.x_max), req.n_points)
    else:
        x = np.linspace(req.x_min, req.x_max, req.n_points)

    factors = [f.as_tuple() for f in req.factors]
    _, K, n, kind = factors[req.sweep]
    factors[req.sweep] = (x, K, n, kind)  # the model accepts an array for x
    r = transcription_rate(factors, r_max=req.r_max, r_basal=req.r_basal)
    return {"x": _nums(x), "r": _nums(r)}


@app.post("/api/gillespie")
def gillespie_trajectory(req: GillespieRequest):
    times, ms, ps = gillespie(req.r, req.gamma_m, req.beta, req.gamma_p, req.t_max,
                              m0=req.m0, p0=req.p0, seed=req.seed,
                              autorepression=_autorepression_tuple(req.autorepression),
                              n_points=req.n_points)
    return {"t": _nums(times), "m": _nums(ms), "p": _nums(ps)}


def _deterministic(r, gamma_m, beta, gamma_p, times, m0, p0, autorepression):
    """Integrate dm/dt and dp/dt (gene_expression.py) with RK4, recording at times.

    With autorepression the rate is the Level 1 repressor formula with x = p:
    r(p) = transcription_rate([(p, K, n, repressor)], r_max=r, r_basal=r_basal).
    """
    def rate(p):
        if autorepression is None:
            return r
        return transcription_rate([(p, autorepression.K, autorepression.n, REPRESSOR)],
                                  r_max=r, r_basal=autorepression.r_basal)

    def deriv(m, p):
        m, p = max(m, 0.0), max(p, 0.0)
        return (mrna_derivative(m, rate(p), gamma_m),
                protein_derivative(p, m, beta, gamma_p))

    interval = times[1] - times[0]
    # about 20 steps per decay time keeps RK4 accurate; cap the total work
    per_interval = math.ceil(interval * max(gamma_m, gamma_p, 1e-12) * 20)
    per_interval = max(1, min(per_interval, MAX_ODE_STEPS // (len(times) - 1)))
    h = interval / per_interval
    m, p = float(m0), float(p0)
    ms, ps = [m], [p]
    for _ in range(len(times) - 1):
        for _ in range(per_interval):
            k1 = deriv(m, p)
            k2 = deriv(m + h / 2 * k1[0], p + h / 2 * k1[1])
            k3 = deriv(m + h / 2 * k2[0], p + h / 2 * k2[1])
            k4 = deriv(m + h * k3[0], p + h * k3[1])
            m += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            p += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        ms.append(m)
        ps.append(p)
    return ms, ps


@app.post("/api/time-course")
def time_course(req: GillespieRequest):
    """Everything the graphing page needs: one stochastic run plus the predicted
    curves from the same starting point.

    Without feedback the prediction is Thattai et al.'s theory: the exact mean and
    standard deviation (sd_m, sd_p) of m and p over time, from the moment equations
    (curve = "paper"). With feedback the paper only gives the steady state (Eq. 3),
    so det_m / det_p are the model's deterministic equations instead and there is no
    spread (curve = "deterministic")."""
    times, ms, ps = gillespie(req.r, req.gamma_m, req.beta, req.gamma_p, req.t_max,
                              m0=req.m0, p0=req.p0, seed=req.seed,
                              autorepression=_autorepression_tuple(req.autorepression),
                              n_points=req.n_points)
    if req.autorepression is None:
        det_m, det_p, var_m, var_p = moments(req.r, req.gamma_m, req.beta, req.gamma_p,
                                             times, req.m0, req.p0)
        curve, sd_m, sd_p = "paper", _nums(np.sqrt(var_m)), _nums(np.sqrt(var_p))
    else:
        det_m, det_p = _deterministic(req.r, req.gamma_m, req.beta, req.gamma_p, times,
                                      req.m0, req.p0, req.autorepression)
        curve, sd_m, sd_p = "deterministic", None, None
    result = {"t": _nums(times), "m": _nums(ms), "p": _nums(ps),
              "det_m": _nums(det_m), "det_p": _nums(det_p),
              "sd_m": sd_m, "sd_p": sd_p, "curve": curve}
    if req.autorepression is None and req.gamma_m > 0 and req.gamma_p > 0:
        result["theory"] = _theory(req.r, req.gamma_m, req.beta, req.gamma_p)
    return result


@app.post("/api/parse-parameters")
def parse_parameter_file(req: ParameterFileRequest):
    """The page sends the uploaded file's text; the rules are in parameter_file.py."""
    values = parse_parameters(req.text)
    return {"values": values, "feedback": uses_feedback(values)}


@app.post("/api/noise-check")
def noise_check(req: NoiseCheckRequest):
    """Measure the steady-state protein mean and Fano factor from NOISE_RUNS long
    runs and compare them with the paper's Eq. 1 (or Eq. 3 with feedback)."""
    autorep = _autorepression_tuple(req.autorepression)
    theory = predict(req.r, req.gamma_m, req.beta, req.gamma_p, autorep)
    if theory is None:
        raise ValueError("the paper's formulas need r, beta, gamma_m and gamma_p all > 0")

    # events per protein lifetime ~ 2 * r * (1 + b) / gamma_p (every molecule made
    # is also destroyed); pick the longest run that fits the budget
    b = req.beta / req.gamma_m
    events_per_lifetime = 2 * req.r * (1 + b) / req.gamma_p
    lifetimes = min(NOISE_MAX_LIFETIMES,
                    NOISE_EVENT_BUDGET / (NOISE_RUNS * max(events_per_lifetime, 1e-12)))
    if lifetimes < NOISE_MIN_LIFETIMES:
        raise ValueError("these rates are too fast to measure the noise quickly; "
                         "try smaller r or beta")
    t_burn = NOISE_BURN_LIFETIMES / req.gamma_p
    t_max = t_burn + lifetimes / req.gamma_p

    means, fanos = [], []
    for k in range(NOISE_RUNS):
        s = steady_state_stats(req.r, req.gamma_m, req.beta, req.gamma_p, t_max,
                               t_burn=t_burn, seed=req.seed + k, autorepression=autorep)
        means.append(s["mean_p"])
        fanos.append(s["fano_p"])

    def summary(values, predicted):
        mean = float(np.mean(values))
        se = float(np.std(values, ddof=1) / math.sqrt(len(values)))
        return {"simulated": _num(mean), "se": _num(se), "paper": _num(predicted),
                "match": bool(abs(mean - predicted) <= MATCH_SE * se)}

    return {"mean": summary(means, theory[0]),
            "fano": summary(fanos, theory[1]),
            "formula": "Eq. 1" if autorep is None else "Eq. 3",
            "runs": NOISE_RUNS, "lifetimes_per_run": _num(lifetimes)}


@app.post("/api/gillespie-stats")
def gillespie_stats(req: GillespieStatsRequest):
    stats = steady_state_stats(req.r, req.gamma_m, req.beta, req.gamma_p, req.t_max,
                               t_burn=req.t_burn, m0=req.m0, p0=req.p0, seed=req.seed,
                               autorepression=_autorepression_tuple(req.autorepression))
    result = {"simulated": {k: _num(v) for k, v in stats.items()}}
    if req.autorepression is None:  # the analytic formulas only hold for constant r
        result["theory"] = _theory(req.r, req.gamma_m, req.beta, req.gamma_p)
    return result


@app.post("/api/cascade")
def cascade(req: CascadeRequest):
    # cascade.py does not validate its inputs, so check them here. k_r1 must be
    # positive: with k_r1 = 0 every reaction rate is 0 and the simulation divides by 0.
    for name in ("k_r1", "k_d", "mrna_half_life", "protein_half_life", "t_max"):
        if getattr(req, name) <= 0:
            raise ValueError(f"{name} must be positive, got {getattr(req, name)}")
    if req.k_r_max < 0:
        raise ValueError(f"k_r_max must be non-negative, got {req.k_r_max}")
    if min(req.b) < 0:
        raise ValueError(f"burst sizes b must be non-negative, got {req.b}")
    if not 0 <= req.t_burn < req.t_max:
        raise ValueError(f"t_burn must be in [0, t_max), got {req.t_burn}")
    means, fanos = cascade_stats(req.k_r1, req.k_r_max, req.k_d, req.b,
                                 req.mrna_half_life, req.protein_half_life,
                                 req.t_max, req.t_burn, seed=req.seed)
    return {"mean": _nums(means), "fano": _nums(fanos)}
