---
name: deployment-status-2026-11-18
description: Where the Vercel deployment stood at end of 2026-11-18 session — API committed locally but not pushed; terminal has no GitHub credentials
metadata:
  node_type: memory
  type: project
  originSessionId: 02458a38-db29-449d-a6f7-9fc6cf0ee127
  modified: 2026-11-18T23:05:28.034Z
---

State at end of session 2026-11-18:
- UPDATE (2026-11-18, later session): Vercel project CREATED via REST API using ~/.vercel-token: name `02120-project-quentinluo` (id prj_jm1kqivnEzZfQ9DCXO6wjs7gI3aR), framework fastapi, linked to GitHub repo (repoId 1438925220) so pushes to main auto-deploy; env NUMBA_CACHE_DIR=/tmp/numba_cache set. First prod deploy READY at https://02120-project-quentinluo.vercel.app ; all 4 TestLiveServer tests pass against it. Items below about the import being pending are superseded.
- Front end: user wants a Desmos-style page (type params, graph redraws live; no submit button, no axis-picker form). Built static/index.html (served at `/` by app.py) + POST /api/time-course (Gillespie trajectory + RK4 deterministic curve + theory). Tests pass locally (34). Pushed as d48ef24 and live at / (5 live tests pass).
- Pushed as 5a455e9 and LIVE (user tested locally first, then approved push): parameter_file.py + /api/parse-parameters + "Load from text file" button; noise_theory.py (Eq.1/Eq.3, no scipy) + /api/noise-check (8 runs, match if within 3 SE); compare mode ("Keep this run"); parameter_files/ with paper_*.txt comparison files (feedback = 0/1 line). 59 tests pass.
- Pushed as a91bf70 and LIVE (2026-11-19): "Show Thattai et al. predicted curve" checkbox; noise_theory.moments() solves the paper's 5 moment equations (matches replicate_fig2 exact values) -> dashed mean + shaded ±1 SD band; with feedback the dashed line is labelled deterministic (no paper formula), no band. 65 tests pass.
- Replication scripts need scipy: .venv lacks it; run with default python3 (/opt/miniconda3).
- UPDATE (later 2026-11-18): `git status -sb` shows `main...origin/main` with no ahead count, so the API commit b34eeb0 has since been pushed. Credentials issue resolved: `git push` from Claude's terminal works (pushed 9b8edf9). `__pycache__/` was untracked in that commit (open question settled).
- Vercel import NOT done yet. Offered 3 routes: (1, recommended) user creates a 7-day Vercel token and saves it to `~/.vercel-token` (chmod 600) via Terminal `read -s`, never pasting it in chat (journal is pushed to GitHub); then Claude creates the project + deploy via Vercel REST API with curl; (2) Claude in Chrome; (3) user clicks Import at vercel.com/new. User hadn't chosen yet. Check `ls ~/.vercel-token` next session.
- Memory is mirrored to repo `claude_memory/` and auto-committed/pushed by `.claude/hooks/sync_memory.py` on SessionEnd.
- Commit "Add FastAPI web API and tests for Vercel deployment" (app.py, test_app.py, requirements.txt, .gitignore) exists locally on `main` but was **not pushed** — `git push` fails here: osxkeychain helper has no github.com entry and there are no SSH keys. Verify with `git status -sb` before assuming it's still unpushed.
- User was told to create a fine-grained PAT (Contents read/write, this repo only) and push once from the macOS Terminal app (username Quentin-Luo-yluo4) so the keychain stores it; or use VS Code / GitHub Desktop. User signs into GitHub with Google, so has no GitHub password.
- After push: user imports the repo at vercel.com/new (FastAPI auto-detected), optionally sets env NUMBA_CACHE_DIR=/tmp/numba_cache, then pastes the URL so live tests run via `APP_URL=<url> .venv/bin/python -m unittest -v test_app.TestLiveServer`.
- Pending: build the front-end page (inputs, axis controls, Plotly) — `/` currently redirects to `/docs`. Open question: user hasn't decided whether to untrack `__pycache__/` (tracked in git, modified by numba on every test run).
- Local env: `.venv` uses system Python 3.9.4 (only Python installed; no Node, Vercel CLI, brew, or gh).

Relates to [[vercel-deployment-requirement]].

- 2026-11-19: README.md written (site features, paper comparisons, replication results, run instructions, AI-use note pointing to journal/). User may add a site screenshot (e.g. figures/website.png) for the README — not done yet.
- Note: pushing to main auto-deploys to Vercel; the .git pre-commit hook runs save_journal.py and adds journal/ to every commit; run save_journal.py by hand with </dev/null or it hangs on stdin.
