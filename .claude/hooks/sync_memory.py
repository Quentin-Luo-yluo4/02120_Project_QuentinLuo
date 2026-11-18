"""Mirror Claude Code's project memory into claude_memory/ and push it to GitHub.

Runs at the end of every Claude Code session (SessionEnd hook). It only ever
stages claude_memory/ and journal/, so other uncommitted work is left alone.
Run it by hand with: python3 .claude/hooks/sync_memory.py
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def project_dir() -> Path:
    env_dir = os.environ.get("CLAUDE_PROJECT_DIR")
    if env_dir:
        return Path(env_dir)
    return Path(__file__).resolve().parent.parent.parent


def memory_dir(proj: Path) -> Path:
    # Claude Code names the folder after the project path, with every
    # character other than letters and digits replaced by "-".
    slug = "".join(c if c.isascii() and c.isalnum() else "-" for c in str(proj))
    return Path.home() / ".claude" / "projects" / slug / "memory"


def git(proj: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=proj, capture_output=True, text=True, timeout=60)


def main() -> int:
    proj = project_dir()
    src = memory_dir(proj)
    if not src.is_dir():
        print(f"sync_memory: no memory folder at {src}")
        return 0

    dest = proj / "claude_memory"
    if dest.exists():
        shutil.rmtree(dest)  # so deleted memories disappear from the repo too
    shutil.copytree(src, dest)

    git(proj, "add", "-A", "--", "claude_memory", "journal")
    if git(proj, "diff", "--cached", "--quiet").returncode == 0:
        print("sync_memory: nothing new to commit")
        return 0

    msg = "Sync Claude memory and build journal\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
    commit = git(proj, "commit", "-m", msg)
    if commit.returncode != 0:
        print("sync_memory: commit failed:", commit.stderr.strip())
        return 0
    push = git(proj, "push", "origin", "HEAD")
    print("sync_memory: pushed" if push.returncode == 0 else "sync_memory: push failed: " + push.stderr.strip())
    return 0


if __name__ == "__main__":
    sys.exit(main())
