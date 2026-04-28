#!/usr/bin/env python3
"""release.py — cut a new SpliceCraft release.

Usage:  ./release.py 0.4.0

Pure-Python replacement for the old release.sh. Same five steps, same
ordering, same abort-on-failure semantics:

  1. Bumps the version in BOTH pyproject.toml and splicecraft.py.
  2. Runs the test suite (aborts on any failure).
  3. Builds the package locally (aborts if twine check fails).
  4. Commits the bump, tags v<version>, and pushes.
  5. GitHub Actions (publish.yml) builds + uploads to PyPI from the tag.

Prereqs (one-time, see README):

  - Trusted Publishing configured at pypi.org for this project.
  - git clean working tree (commit everything else first).

Aborts at the first error — every shelled-out step uses ``check=True`` so
a failed test / build / push raises ``CalledProcessError`` and the script
exits non-zero before tagging or pushing. Mirrors the old script's
``set -euo pipefail`` guarantee.

The whole thing takes ~2 minutes from run to "live on PyPI".
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
PYPROJECT = REPO_ROOT / "pyproject.toml"
SPLICECRAFT = REPO_ROOT / "splicecraft.py"

# Accept canonical X.Y.Z plus PEP-440 suffixes (rc1, post1, dev0, …).
_VERSION_RE = re.compile(r"^\d+\.\d+\.\d+([a-z0-9.-]*)?$")

# In-file regexes for the version bump. Both files live at the repo
# root; the regexes target the canonical lines that release.sh used to
# rewrite with `sed -i -E`.
_PYPROJECT_VERSION_RE = re.compile(
    r'^(version\s*=\s*)"[^"]+"', re.MULTILINE,
)
_SPLICECRAFT_VERSION_RE = re.compile(
    r'^(__version__\s*=\s*)"[^"]+"', re.MULTILINE,
)


def _heading(msg: str) -> None:
    print(f"── {msg} ──", flush=True)


def _die(msg: str, code: int = 1) -> None:
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(code)


def _run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    """Run *cmd* with ``check=True`` so any non-zero exit aborts the
    release. We deliberately do not capture stdout/stderr — the user
    needs to see pytest / twine / git output as it streams."""
    return subprocess.run(cmd, check=True, **kwargs)


def _ensure_clean_tree() -> None:
    """Refuse to proceed if there are uncommitted changes — same
    invariant release.sh enforced via ``git status --porcelain``."""
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        check=True, capture_output=True, text=True,
    )
    if result.stdout.strip():
        print("Error: working tree is dirty. Commit or stash first.",
              file=sys.stderr)
        subprocess.run(["git", "status", "--short"], check=False)
        sys.exit(1)


def _ensure_tag_unused(version: str) -> None:
    """Refuse to proceed if the target tag already exists locally."""
    rev = subprocess.run(
        ["git", "rev-parse", f"v{version}"],
        capture_output=True, text=True,
    )
    if rev.returncode == 0:
        _die(f"tag v{version} already exists.")


def _bump_version_in_file(path: Path, pattern: re.Pattern[str],
                          new_version: str, label: str) -> None:
    """Substitute the canonical version line in *path* via regex. Raises
    if zero substitutions happened — protects against the file's
    formatting drifting (e.g., quotes changing) without us noticing."""
    text = path.read_text(encoding="utf-8")
    new_text, n = pattern.subn(rf'\1"{new_version}"', text)
    if n == 0:
        _die(f"failed to find version line in {label} ({path.name}).")
    if n > 1:
        # More than one match would mean we're rewriting something
        # unintended; bail before clobbering it.
        _die(f"found {n} version-line matches in {label}; refusing to "
             f"bump (please simplify {path.name}).")
    path.write_text(new_text, encoding="utf-8")


def _verify_bump(path: Path, new_version: str, var_name: str) -> None:
    """Read the file back and confirm the new version line is exactly
    where we expect it. Catches odd encoding / line-ending issues that
    would otherwise let a silently-failed bump through."""
    text = path.read_text(encoding="utf-8")
    expected = re.compile(
        rf'^{re.escape(var_name)}\s*=\s*"{re.escape(new_version)}"',
        re.MULTILINE,
    )
    if not expected.search(text):
        _die(f"failed to update {path.name} (expected "
             f'{var_name} = "{new_version}").')


def _clean_build_artifacts() -> None:
    """``rm -rf dist/ build/ *.egg-info`` in pure Python."""
    for d in ("dist", "build"):
        shutil.rmtree(REPO_ROOT / d, ignore_errors=True)
    for egg in REPO_ROOT.glob("*.egg-info"):
        shutil.rmtree(egg, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Cut a new SpliceCraft release.",
    )
    parser.add_argument(
        "version",
        help="New version (X.Y.Z, optionally with PEP-440 suffix).",
    )
    args = parser.parse_args(argv)
    new_version = args.version

    if not _VERSION_RE.match(new_version):
        _die(f"version must look like X.Y.Z (got {new_version!r}).")

    _ensure_clean_tree()
    _ensure_tag_unused(new_version)

    _heading(f"Bumping version to {new_version}")
    _bump_version_in_file(
        PYPROJECT, _PYPROJECT_VERSION_RE, new_version, "pyproject.toml",
    )
    _bump_version_in_file(
        SPLICECRAFT, _SPLICECRAFT_VERSION_RE, new_version, "splicecraft.py",
    )
    _verify_bump(PYPROJECT,   new_version, "version")
    _verify_bump(SPLICECRAFT, new_version, "__version__")

    _heading("Running test suite")
    _run([sys.executable, "-m", "pytest", "-q", "--tb=short"])

    _heading("Building sdist + wheel")
    _clean_build_artifacts()
    _run([sys.executable, "-m", "build"])

    _heading("Verifying package metadata")
    # Pass the dist files explicitly so we don't depend on shell glob
    # expansion. ``twine check dist/*`` would fail under a literal-glob
    # interpreter (Windows cmd, restricted shells); enumerate via Path.
    dist_files = sorted(str(p) for p in (REPO_ROOT / "dist").iterdir())
    if not dist_files:
        _die("no build artifacts found in dist/.")
    _run([sys.executable, "-m", "twine", "check", *dist_files])

    _heading("Committing + tagging + pushing")
    _run(["git", "add", "pyproject.toml", "splicecraft.py"])
    _run(["git", "commit", "-m", f"Release v{new_version}"])
    _run(["git", "tag", f"v{new_version}"])
    _run(["git", "push", "origin", "master"])
    _run(["git", "push", "origin", f"v{new_version}"])

    print()
    print("═" * 61)
    print(f" Release v{new_version} pushed.")
    print(" GitHub Actions will publish to PyPI in ~2 minutes.")
    print(" Watch:  https://github.com/Binomica-Labs/SpliceCraft/actions")
    print(" Verify: https://pypi.org/project/splicecraft/")
    print("═" * 61)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except subprocess.CalledProcessError as exc:
        # Subprocess output already streamed; surface the failing
        # command + exit code so the user knows where the release
        # halted.
        cmd = " ".join(exc.cmd) if isinstance(exc.cmd, list) else str(exc.cmd)
        print(f"\nRelease aborted — `{cmd}` exited with code {exc.returncode}.",
              file=sys.stderr)
        sys.exit(exc.returncode or 1)
    except KeyboardInterrupt:
        print("\nRelease cancelled.", file=sys.stderr)
        sys.exit(130)
