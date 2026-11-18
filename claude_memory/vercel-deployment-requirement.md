---
name: vercel-deployment-requirement
description: Instructor requires the 02-120 project web app be deployed on Vercel; user wants local-first dev then deploy in under an hour
metadata:
  node_type: memory
  type: project
  originSessionId: 02458a38-db29-449d-a6f7-9fc6cf0ee127
  modified: 2026-11-18T23:05:21.251Z
---

The course instructor requires the project's interactive web app to be hosted on **Vercel** (user's account: vercel.com/quentin-luo-yluo4, signed in via GitHub). The user wants to do almost all coding locally and move it onto the website in under an hour.

**Why:** Streamlit / matplotlib-widget UIs were ruled out — Streamlit can't run on Vercel's serverless functions, and matplotlib widgets don't work in a browser.

**How to apply:** Architecture chosen (2026-11-18): FastAPI `app.py` at repo root (Vercel zero-config detection) calling the existing model modules, plus a static HTML page with Plotly.js for inputs and x/y axis controls (next step, not built yet). Keep simulations within Vercel limits (Hobby 300 s duration, 4.5 MB responses); numba cache must go to /tmp on Vercel. See [[deployment-status-2026-11-18]] and [[user-profile]].
