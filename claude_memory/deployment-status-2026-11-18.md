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
- UPDATE (later 2026-11-18): `git status -sb` shows `main...origin/main` with no ahead count, so the API commit b34eeb0 has since been pushed. Credentials issue resolved: `git push` from Claude's terminal works (pushed 9b8edf9). `__pycache__/` was untracked in that commit (open question settled).
- Vercel import NOT done yet. Offered 3 routes: (1, recommended) user creates a 7-day Vercel token and saves it to `~/.vercel-token` (chmod 600) via Terminal `read -s`, never pasting it in chat (journal is pushed to GitHub); then Claude creates the project + deploy via Vercel REST API with curl; (2) Claude in Chrome; (3) user clicks Import at vercel.com/new. User hadn't chosen yet. Check `ls ~/.vercel-token` next session.
- Memory is mirrored to repo `claude_memory/` and auto-committed/pushed by `.claude/hooks/sync_memory.py` on SessionEnd.
- Commit "Add FastAPI web API and tests for Vercel deployment" (app.py, test_app.py, requirements.txt, .gitignore) exists locally on `main` but was **not pushed** — `git push` fails here: osxkeychain helper has no github.com entry and there are no SSH keys. Verify with `git status -sb` before assuming it's still unpushed.
- User was told to create a fine-grained PAT (Contents read/write, this repo only) and push once from the macOS Terminal app (username Quentin-Luo-yluo4) so the keychain stores it; or use VS Code / GitHub Desktop. User signs into GitHub with Google, so has no GitHub password.
- After push: user imports the repo at vercel.com/new (FastAPI auto-detected), optionally sets env NUMBA_CACHE_DIR=/tmp/numba_cache, then pastes the URL so live tests run via `APP_URL=<url> .venv/bin/python -m unittest -v test_app.TestLiveServer`.
- Pending: build the front-end page (inputs, axis controls, Plotly) — `/` currently redirects to `/docs`. Open question: user hasn't decided whether to untrack `__pycache__/` (tracked in git, modified by numba on every test run).
- Local env: `.venv` uses system Python 3.9.4 (only Python installed; no Node, Vercel CLI, brew, or gh).

Relates to [[vercel-deployment-requirement]].
