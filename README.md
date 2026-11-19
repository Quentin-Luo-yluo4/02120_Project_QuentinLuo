# Noise in Gene Expression: replicating Thattai & van Oudenaarden (2001)

**02-120 Programming for Scientists, final project, Quentin Luo**

**Live website:** https://02120-project-quentinluo.vercel.app

This project rebuilds, in Python, the model of gene expression noise from

> M. Thattai and A. van Oudenaarden, *Intrinsic noise in gene regulatory networks*,
> PNAS 98:8614–8619 (2001).

A gene makes mRNA, and each mRNA makes protein. Every one of these events is random, so
two identical cells end up with different numbers of protein molecules. The paper
predicts how big that cell-to-cell difference (the "noise") is and what controls it.
This project simulates the process one molecule at a time with the Gillespie algorithm,
checks the results against the paper's formulas, and puts everything on an interactive
website.

---

## The website

Open https://02120-project-quentinluo.vercel.app. Nothing to install.

It works like a graphing calculator: type a parameter or drag its slider, and the graph
redraws immediately. There is no Submit button.

| Feature | What it does |
|---|---|
| Parameter boxes and sliders | transcription rate `r`, mRNA decay `γₘ`, translation rate `β`, protein decay `γₚ`, simulated time `t`; optional negative feedback (`K`, `n`, `r₀`) |
| Graph | one stochastic (Gillespie) run of protein and mRNA over time |
| **Show Thattai et al. predicted curve** (checkbox) | dashed line = the paper's predicted average; shaded band = predicted average ± 1 standard deviation, which changes over time (including the noise overshoot after the gene switches on, Fig. 2) |
| Paper's numbers (under the graph) | steady-state mean and Fano factor from the paper's Eq. 1, or Eq. 3 with feedback |
| Measured noise (under the graph) | 8 long simulations give the protein mean and Fano factor with a standard error, marked **✓ match** if within 3 standard errors of the paper |
| Keep this run | freezes the current curve so the next run is drawn on the same axes (up to 4) |
| Load from text file | fills in all parameters from a `.txt` file (format below) |

### Reproducing the paper's conclusions on the website

The files in [`parameter_files/`](parameter_files/) use the paper's parameters (from its
figure captions; time in seconds). For each conclusion: load the first file, press
**Keep this run**, then load the second.

| Paper's conclusion | First file | Second file | What you should see |
|---|---|---|---|
| Noise is set by the burst size b (proteins per mRNA), not by the amount of protein | `paper_parameters.txt` (b = 20) | `paper_small_bursts.txt` (b = 5) | same mean (1039), Fano factor 20.35 → 5.84, much narrower curve |
| Transcription rate changes the amount of protein, not the Fano factor | `paper_parameters.txt` | `paper_more_transcription.txt` | mean doubles (1039 → 2077), Fano factor stays 20.35 |
| Negative feedback reduces noise | `paper_no_feedback_same_mean.txt` | `paper_feedback_parameters.txt` | same mean (689), Fano factor 10.68 → 5.78 |

Measured on the live site, all five files match the paper:

| File | Fano factor, simulated | Paper |
|---|---|---|
| `paper_parameters.txt` | 19.85 ± 0.17 | 20.35 (Eq. 1) |
| `paper_small_bursts.txt` | 5.82 ± 0.07 | 5.84 (Eq. 1) |
| `paper_more_transcription.txt` | 20.01 ± 0.38 | 20.35 (Eq. 1) |
| `paper_no_feedback_same_mean.txt` | 10.68 ± 0.12 | 10.68 (Eq. 1) |
| `paper_feedback_parameters.txt` | 5.87 ± 0.05 | 5.78 (Eq. 3) |

### Parameter file format

One `name = value` per line. `#` starts a comment, blank lines are ignored, and any
parameter left out keeps its current value on the page.

```
r = 0.01
gamma_m = 0.0057762265
beta = 0.1155245301
gamma_p = 0.0001925409
t_max = 36000
feedback = 0          # 1 = protein represses its own gene (then also give K, n, r_basal)
```

Mistakes are reported with the line number, e.g. `line 3: unknown parameter 'gamma'`.
The same reader works from Python: `parameter_file.read_parameter_file("my_file.txt")`.

### How the website is built

```
browser (static/index.html, Plotly.js)
   │  sends the parameters as JSON
   ▼
FastAPI server (app.py), hosted on Vercel
   │  calls the project's own model code; app.py only converts requests and results
   ▼
gillespie.py        stochastic simulation (numba-compiled)
gene_expression.py  rate equations
noise_theory.py     the paper's formulas: Eq. 1, Eq. 3, and the moment equations over time
parameter_file.py   reads uploaded parameter files
```

Vercel is connected to this GitHub repository and redeploys automatically on every push
to `main`, so the website always runs the code in this repository.

API endpoints (interactive documentation at
[/docs](https://02120-project-quentinluo.vercel.app/docs)):

| Endpoint | Returns |
|---|---|
| `POST /api/time-course` | one simulated run plus the predicted mean and ± 1 SD curves (used by the graph) |
| `POST /api/noise-check` | measured mean and Fano factor from 8 long runs vs the paper's formula |
| `POST /api/parse-parameters` | the values read from an uploaded parameter file |
| `POST /api/steady-state` | transcription rate and steady-state mRNA and protein for given transcription factors |
| `POST /api/rate-curve` | transcription rate as one transcription factor's concentration is swept |
| `POST /api/gillespie`, `/api/gillespie-stats` | a raw trajectory, or its time-averaged statistics |
| `POST /api/cascade` | mean and Fano factor of each gene in a cascade (Fig. 5) |
| `GET /api/health` | whether the server is running |

---

## Replicating the paper's figures

Two scripts compare the simulation with the paper's formulas, using the parameters from
the figure captions. A point counts as a **match** if the simulation is within 3 standard
errors of the formula. The plots are in [`figures/`](figures/) and
[`fig2_replication.png`](fig2_replication.png): black lines are the paper's formulas,
symbols with error bars are this project's simulation.

| Figure | What it shows | Result |
|---|---|---|
| Fig. 1c | noise vs mean as burst size, transcription rate and protein lifetime vary | 24 / 24 points match |
| Fig. 2 | noise over time after a gene switches on, including the overshoot | 100 / 100 time points match, for both mRNA half-lives |
| Figs. 3a, 4a–b | protein distributions with and without feedback | mean 1199 vs 1200 and SD 111.8 vs 113.2 (no feedback); 688 vs 689 and 62.2 vs 63.1 (feedback) |
| Fig. 3b | feedback lowers noise | 28 / 44 points match (see below) |
| Fig. 5 | noise passed along a 6-gene cascade | 24 / 24 points match |

**Fig. 3b misses.** All 16 misses are for very steep feedback (Hill coefficient n ≥ 6) or
very strong repression (K_d = 100), where the simulation is 3–14 % noisier than Eq. 3.
The paper derives Eq. 3 by approximating the feedback curve with a straight line, which it
says is exact only "in the limit of small fluctuations". The steeper the feedback, the
worse that approximation, which matches where the gap grows. The script also confirms
that its implementation of Eq. 3 agrees with an independent exact linear calculation to
10⁻¹⁵, so the gap comes from the approximation, not from the code.

---

## Running it yourself

Requires Python 3.9 or newer.

```bash
pip install -r requirements.txt httpx uvicorn   # website + tests
pip install scipy matplotlib                    # for replicate_*.py and try_parameters.py
```

| Task | Command |
|---|---|
| Run the website locally | `uvicorn app:app --reload`, then open http://127.0.0.1:8000 |
| Run all tests (93 tests, a few seconds) | `python3 -m unittest` |
| Test the deployed website | `APP_URL=https://02120-project-quentinluo.vercel.app python3 -m unittest -v test_app.TestLiveServer` |
| Replicate Figs. 1, 3, 4, 5 (a few minutes) | `python3 replicate_figures.py` |
| Replicate Fig. 2 (about a minute) | `python3 replicate_fig2.py` |
| Check your own parameters against the paper | `python3 try_parameters.py --kR 0.01 --b 20 --mrna-half-life 120 --protein-half-life 3600` |
| Command-line calculator (Levels 1–3) | `python3 interactive.py` |

## Files

| File | Purpose |
|---|---|
| `gene_expression.py` | Level 1 transcription regulation (Hill functions), Level 2–3 mRNA and protein rate equations |
| `gillespie.py` | Gillespie simulation of one gene, optionally with negative autoregulation |
| `cascade.py` | Gillespie simulation of a gene cascade (Fig. 5) |
| `noise_theory.py` | the paper's predictions: Eq. 1, Eq. 3 and the time-dependent moment equations |
| `parameter_file.py` | reads `name = value` parameter files |
| `app.py` | the web server (FastAPI) |
| `static/index.html` | the web page |
| `parameter_files/` | the paper's parameters and the comparison files |
| `replicate_figures.py`, `replicate_fig2.py` | figure replication against the paper's formulas |
| `try_parameters.py`, `interactive.py` | command-line tools |
| `test_*.py` | unit tests (`unittest`) |
| `Proposal_Proj_QuentinLuo_yluo4.pdf` | project proposal |
| `journal/` | build journal (see below) |

## Use of AI tools

Parts of this project, including the website and its deployment, were built with the help
of an AI coding assistant (Claude Code). Every working session is logged in
[`journal/`](journal/), with readable versions in [`journal/readable/`](journal/readable/):
the requests made, the code written, and the decisions taken. The git history shows how
the project was built step by step.
