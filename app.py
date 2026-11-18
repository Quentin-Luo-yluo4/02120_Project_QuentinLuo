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
from fastapi.responses import JSONResponse, RedirectResponse  # noqa: E402
from pydantic import BaseModel, Field  # noqa: E402

from cascade import cascade_stats  # noqa: E402
from gene_expression import (  # noqa: E402
    ACTIVATOR,
    R_BASAL,
    R_MAX,
    REPRESSOR,
    hill,
    mrna_steady_state,
    protein_steady_state,
    transcription_rate,
)
from gillespie import gillespie, steady_state_stats  # noqa: E402

# Caps that keep one request inside Vercel's time and 4.5 MB response limits.
MAX_CURVE_POINTS = 2000
MAX_TRAJECTORY_POINTS = 5000
MAX_CASCADE_GENES = 10

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
    """Until the web page exists, send visitors to the interactive API docs."""
    return RedirectResponse("/docs")


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
