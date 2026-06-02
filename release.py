#!/usr/bin/env python3
"""release.py — cut a new SpliceCraft release.

Usage:  ./release.py 0.8.7
        ./release.py --bioconda-only           # re-submit current recipe

Pure-Python replacement for the old release.sh. Same seven steps, same
ordering, same abort-on-failure semantics:

  1. Bumps the version in BOTH pyproject.toml and splicecraft.py.
  2. Runs the test suite (aborts on any failure).
  3. Builds the package locally (aborts if twine check fails).
  4. Syncs `conda-recipe/meta.yaml` (version + sha256 + run-deps) from
     the just-built sdist + pyproject — keeps the in-repo recipe in
     lockstep with PyPI so a bioconda PR is one click away.
  5. Bundles ALL accumulated working-tree changes (version bump +
     anything else the user has iterated on since the last release)
     into a single ``Release v<version>`` commit, tags v<version>,
     and pushes — then creates a GitHub Release for the tag (changelog
     section as the body + the built wheel/sdist as assets; optional +
     non-fatal, skipped cleanly when `gh` is absent/unauthed).
  6. GitHub Actions (publish.yml) builds + uploads to PyPI from the tag.
  7. Polls PyPI until the sdist is hosted, then forks
     bioconda-recipes (if needed) and opens a PR with the recipe
     update so end users can `conda install -c bioconda splicecraft`.

Sweep #18 (2026-05-21): the working tree no longer needs to be
clean before release. release.py prints a summary of pending
changes, then ``git add -A``'s them into the release commit. This
keeps the iteration cycle tight — the user accumulates work in the
tree and ships it all at release time instead of paying a commit
toll on every change. `.gitignore` is still in effect, so build
artefacts / caches / logs stay out of the release commit.

Prereqs (one-time, see README):

  - Trusted Publishing configured at pypi.org for this project.
  - `gh` CLI authenticated (`gh auth status`) — needed for the
    bioconda fork + PR step. Skipped with a note if missing.

Aborts at the first error — every shelled-out step uses ``check=True`` so
a failed test / build / push raises ``CalledProcessError`` and the script
exits non-zero before tagging or pushing. Mirrors the old script's
``set -euo pipefail`` guarantee.

The whole thing takes ~2 minutes from run to "live on PyPI".
"""
from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import NoReturn

REPO_ROOT = Path(__file__).resolve().parent
PYPROJECT = REPO_ROOT / "pyproject.toml"
SPLICECRAFT = REPO_ROOT / "splicecraft.py"
CHANGELOG = REPO_ROOT / "CHANGELOG.md"
CONDA_RECIPE = REPO_ROOT / "conda-recipe" / "meta.yaml"

# Upstream bioconda repository — the canonical channel every
# `conda install -c bioconda <pkg>` user pulls from. The release flow
# forks this into the maintainer's account, pushes a branch with the
# updated recipe, and opens a PR against this upstream repo. Once
# merged by a bioconda maintainer, bioconda's CI builds the package
# and uploads to anaconda.org/bioconda — that's how the recipe
# actually reaches end users.
BIOCONDA_UPSTREAM = "bioconda/bioconda-recipes"
BIOCONDA_RECIPE_SUBPATH = "recipes/splicecraft/meta.yaml"
PYPI_SDIST_URL_TPL = (
    "https://pypi.io/packages/source/s/splicecraft/"
    "splicecraft-{version}.tar.gz"
)

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

# conda-recipe/meta.yaml uses jinja2-style `{% set version = "..." %}`
# at the top and a `sha256: <hex>` line under `source:`. The run-deps
# block under `requirements:` is a YAML list we rewrite wholesale from
# pyproject's `[project] dependencies` so the two never drift.
_CONDA_VERSION_RE = re.compile(
    r'(\{%\s*set\s+version\s*=\s*)"[^"]+"(\s*%\})'
)
_CONDA_SHA256_RE = re.compile(
    r'^(\s*sha256:\s*)[0-9a-fA-F]{64}', re.MULTILINE,
)
# `  run:` followed by one or more `    - <name> ...` list items, up to
# the next less-indented line. Captures the whole block for replacement.
_CONDA_RUN_BLOCK_RE = re.compile(
    r"^  run:\n(?:    -.*\n)+", re.MULTILINE,
)
# Match a single "package>=X.Y.Z" or "package" line inside pyproject's
# `dependencies = [ ... ]` array (PEP 621).
_PYPROJECT_DEPS_BLOCK_RE = re.compile(
    r"^dependencies\s*=\s*\[\s*\n(.*?)^\]",
    re.MULTILINE | re.DOTALL,
)
_PEP440_SPEC_RE = re.compile(r"^([A-Za-z0-9._-]+)\s*(.*)$")


def _heading(msg: str) -> None:
    print(f"── {msg} ──", flush=True)


def _die(msg: str, code: int = 1) -> NoReturn:
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(code)


def _run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    """Run *cmd* with ``check=True`` so any non-zero exit aborts the
    release. We deliberately do not capture stdout/stderr — the user
    needs to see pytest / twine / git output as it streams."""
    return subprocess.run(cmd, check=True, **kwargs)


def _summarize_pending_changes() -> None:
    """Sweep #18 (2026-05-21): release.py now BUNDLES accumulated
    work into the release commit instead of refusing to proceed.
    Pre-sweep `_ensure_clean_tree` rejected any uncommitted changes;
    the new flow lets the user accumulate iterative work in the
    working tree between releases and ship it in a single
    ``Release vX.Y.Z`` commit. This function is informational only:
    print the pending changes so the user can confirm visually
    before the build/test/commit/push pipeline proceeds.
    """
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        check=True, capture_output=True, text=True,
    )
    pending = result.stdout.strip()
    if not pending:
        return
    print("─" * 61)
    print(" Pending changes that will be bundled into this release:")
    print("─" * 61)
    subprocess.run(["git", "status", "--short"], check=False)
    print()


# Pre-sweep #18 alias — kept around for any external invocation that
# imports release.py and calls the old name. New code should call
# ``_summarize_pending_changes`` directly.
_ensure_clean_tree = _summarize_pending_changes


def _ensure_tag_unused(version: str) -> None:
    """Refuse to proceed if the target tag already exists locally."""
    rev = subprocess.run(
        ["git", "rev-parse", f"v{version}"],
        capture_output=True, text=True,
    )
    if rev.returncode == 0:
        _die(f"tag v{version} already exists.")


def _previous_release_ref() -> str | None:
    """Return the most recent SpliceCraft release tag (``vX.Y.Z``) for
    bounding ``git log`` when drafting a changelog section. Falls back
    to ``None`` if no such tag exists (first-ever release on a fresh
    clone), in which case the caller emits a placeholder rather than
    dumping the whole history.
    """
    result = subprocess.run(
        ["git", "tag", "--list", "v*", "--sort=-v:refname"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        return None
    for line in result.stdout.splitlines():
        tag = line.strip()
        if tag:
            return tag
    return None


def _commits_since(ref: "str | None") -> list[str]:
    """Return the subject line of every non-merge commit since *ref*
    (or every commit on the branch if *ref* is None / unknown).
    Strips out the bookkeeping commits ("Release v…", standalone
    "Changelog: …") that would clutter a user-facing summary.
    """
    if ref is None:
        rng = ["HEAD"]
    else:
        rng = [f"{ref}..HEAD"]
    result = subprocess.run(
        ["git", "log", "--no-merges", "--pretty=format:%s", *rng],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        return []
    out: list[str] = []
    for ln in result.stdout.splitlines():
        s = ln.strip()
        if not s:
            continue
        low = s.lower()
        if low.startswith("release v"):
            continue
        if low.startswith("changelog:"):
            continue
        out.append(s)
    return out


def _draft_changelog_section(version: str) -> str:
    """Build a fresh ``## [<version>]`` section from the commits landed
    since the previous release tag. Used by `_ensure_changelog_entry`
    when no entry has been hand-written yet — guarantees the What's
    New modal renders an up-to-date brief on every release, never a
    stale one.

    Commit subjects in this repo are descriptive enough to make a
    reasonable user-facing summary; the maintainer is free to write a
    richer entry by hand before running release.py and the auto-draft
    will not overwrite it.
    """
    import datetime as _dt
    today = _dt.date.today().isoformat()
    prev = _previous_release_ref()
    commits = _commits_since(prev)
    if commits:
        bullets = "\n".join(f"* {c}" for c in commits)
        provenance = (
            f"_(auto-generated from commits since {prev})_"
            if prev is not None else
            "_(auto-generated changelog)_"
        )
        body = f"{provenance}\n\n{bullets}\n"
    else:
        # No commits since the last tag (e.g. a re-release of the same
        # tree, or release.py is invoked twice). Still emit a heading
        # so the modal has something to display rather than rendering
        # the previous version as if it were current.
        body = (
            "_(auto-generated changelog — no notable commits found "
            "since the previous release)_\n"
        )
    return (
        f"## [{version}] — {today}\n\n"
        f"{body}\n---\n\n"
    )


def _insert_changelog_section(section: str) -> None:
    """Insert *section* at the top of ``CHANGELOG.md``, just below the
    file's header + first ``---`` separator. Preserves all existing
    entries verbatim. Caller has already confirmed the new version's
    heading is absent.
    """
    text = CHANGELOG.read_text(encoding="utf-8")
    # Canonical file shape:
    #   # SpliceCraft Changelog
    #   <blank>
    #   ---
    #   <blank>
    #   ## [X.Y.Z] — …
    # New sections insert between the leading `---` and the most
    # recent `## [` heading. If the anchor is missing (mangled file),
    # fall back to prepending immediately after the H1.
    marker = "\n---\n\n"
    idx = text.find(marker)
    if idx == -1:
        prefix, _, rest = text.partition("\n")
        new_text = f"{prefix}\n\n{section}{rest}"
    else:
        split_at = idx + len(marker)
        new_text = text[:split_at] + section + text[split_at:]
    CHANGELOG.write_text(new_text, encoding="utf-8")


_AUTO_STUB_MARKER = (
    "_(auto-generated changelog — no notable commits found "
    "since the previous release)_"
)


def _promote_unreleased_to_version(version: str) -> bool:
    """Sweep #36 (2026-05-27): if the maintainer has been keeping a
    `## [unreleased]` block updated with the real changelog content,
    relabel that heading to `## [<version>] — <today>` so the rich
    content becomes the new version's body. Returns True if a
    promotion happened, False otherwise.

    Why this exists: the prior flow auto-generated a fresh section
    from `git log` commit subjects, but the user often accumulates
    bullet-pointed `[unreleased]` content as they go and forgets to
    hand-copy it into the new version heading before releasing.
    Result: 0.9.27 / 0.9.28 / 0.9.29 all shipped with the
    "no notable commits found" auto-stub even though substantial
    work had landed. Auto-promotion makes the `[unreleased]` pocket
    the source of truth.

    Only the FIRST `## [unreleased]` heading is promoted — any
    later ones (legacy duplicates from manual edits) stay put so
    the maintainer can fold them in by hand.
    """
    import datetime as _dt
    text = CHANGELOG.read_text(encoding="utf-8")
    # Anchor: heading starts at col 0, matches `## [unreleased]`
    # case-insensitively to allow `[Unreleased]` too.
    m = re.search(
        r"^## \[unreleased\][^\n]*$",
        text, flags=re.IGNORECASE | re.MULTILINE,
    )
    if m is None:
        return False
    today = _dt.date.today().isoformat()
    new_heading = f"## [{version}] — {today}"
    new_text = text[:m.start()] + new_heading + text[m.end():]
    CHANGELOG.write_text(new_text, encoding="utf-8")
    return True


def _changelog_section_body(version: str) -> str:
    """Return the body of `## [<version>]` (everything between that
    heading and the next `## [` heading) so callers can inspect it.
    Empty string if the heading isn't present.
    """
    if not CHANGELOG.is_file():
        return ""
    text = CHANGELOG.read_text(encoding="utf-8")
    head_re = re.compile(
        rf"^## \[{re.escape(version)}\][^\n]*$",
        re.MULTILINE,
    )
    m = head_re.search(text)
    if not m:
        return ""
    tail = text[m.end():]
    next_head = re.search(r"^## \[", tail, flags=re.MULTILINE)
    return tail[:next_head.start()] if next_head else tail


def _section_has_real_content(body: str) -> bool:
    """Return True if *body* contains real changelog content (real
    bullet points or sub-headings), False if it's empty / blank /
    just the auto-stub marker. Drives the empty-stub refusal in
    `_ensure_changelog_entry`.

    Heuristic: strip the heading-trailing `---`, blank lines, and
    the auto-stub marker; if anything substantive remains, accept.
    `_draft_changelog_section`'s commit-bullet body has real `*`
    bullets and passes; the no-notable-commits body strips to
    empty and fails.
    """
    text = body
    text = text.replace(_AUTO_STUB_MARKER, "")
    text = re.sub(r"^---\s*$", "", text, flags=re.MULTILINE)
    text = text.strip()
    return bool(text)


def _ensure_changelog_entry(version: str) -> None:
    """Make sure `CHANGELOG.md` carries a heading for the new version.

    Resolution order (sweep #36, 2026-05-27):
      1. If `## [<version>]` already exists, leave it alone — the
         maintainer hand-wrote it.
      2. Else, if `## [unreleased]` exists, relabel it. The pocket
         is the source of truth for "what shipped this release".
      3. Else, auto-generate a section from `git log
         <previous tag>..HEAD` so the in-app What's New modal at
         least lists the commits since the previous release.

    After resolution, asserts the section body has real content
    (anything beyond the no-notable-commits auto-stub). If it
    doesn't, aborts with a clear message asking the maintainer to
    populate the section first — the "0.9.29 shipped with an
    empty stub" failure mode triggered this guard.

    Missing `CHANGELOG.md` falls through with a friendlier
    message than ``FileNotFoundError``. The resulting file is
    added to the release commit downstream (see `add_targets` in
    `main()`), so the file on disk and the modal stay in lockstep.
    """
    if not CHANGELOG.is_file():
        _die(f"{CHANGELOG.name} not found at {CHANGELOG}. "
             "Add it before releasing.")
    text = CHANGELOG.read_text(encoding="utf-8")
    needle = f"## [{version}]"
    if needle in text:
        if not _section_has_real_content(
            _changelog_section_body(version)
        ):
            _die(
                f"{CHANGELOG.name} has a hand-written `{needle}` "
                f"heading but its body contains no real content "
                f"(only the no-notable-commits stub or blanks). "
                f"Populate it with the actual bug-fix / feature "
                f"bullets before re-running release.py."
            )
        return
    if _promote_unreleased_to_version(version):
        print(
            f"  ↳ promoted `## [unreleased]` → `{needle}` in "
            f"{CHANGELOG.name} (the pocket had queued content)"
        )
    else:
        section = _draft_changelog_section(version)
        _insert_changelog_section(section)
        print(
            f"  ↳ auto-generated `{needle}` section in "
            f"{CHANGELOG.name} (no hand-written entry, no "
            f"`[unreleased]` pocket)"
        )
    if not _section_has_real_content(
        _changelog_section_body(version)
    ):
        _die(
            f"{CHANGELOG.name} section for `{needle}` ended up "
            f"empty (no `[unreleased]` content, no commits since "
            f"the previous tag). Either: (a) hand-write the "
            f"section with bullets describing what's in this "
            f"release, or (b) add bullets under a "
            f"`## [unreleased]` heading and re-run release.py. "
            f"Refusing to ship a What's New modal with no real "
            f"content."
        )


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


def _read_pyproject_runtime_deps() -> list[str]:
    """Extract the [project] dependencies array from pyproject.toml in
    declaration order. Returns the raw PEP-440 specifiers (e.g.
    ``"textual>=8.2.6"``).

    Regex-based to avoid a tomllib (3.11+) / tomli dependency — the
    pyproject's dependency block uses one quoted string per line with
    comments and blank lines interspersed, all of which we skip.
    """
    text = PYPROJECT.read_text(encoding="utf-8")
    m = _PYPROJECT_DEPS_BLOCK_RE.search(text)
    if m is None:
        _die("could not find `dependencies = [...]` in pyproject.toml.")
    deps: list[str] = []
    for line in m.group(1).splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        # Strip trailing comma + trailing inline comment.
        no_comment = stripped.split("#", 1)[0].strip().rstrip(",").strip()
        m2 = re.match(r'^"([^"]+)"$', no_comment)
        if m2:
            deps.append(m2.group(1))
    if not deps:
        _die("pyproject's dependencies array parsed empty — refusing "
             "to wipe the conda recipe's run-deps.")
    return deps


def _pep440_to_conda(spec: str) -> str:
    """Translate a single PEP-440 spec to conda's `pkg >=X.Y.Z` form.

    Conda's recipe parser wants whitespace between the package name and
    the version constraint; pyproject (PEP 508) accepts both styles but
    we normalise here so the recipe stays consistent regardless of how
    pyproject is written.
    """
    m = _PEP440_SPEC_RE.match(spec.strip())
    if not m:
        return spec
    name, rest = m.group(1), m.group(2).strip()
    return f"{name} {rest}".strip() if rest else name


def _sha256_of(path: Path) -> str:
    """SHA-256 hex digest of a file, streamed in 64 KB chunks so we
    don't load multi-MB sdists fully into RAM."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(64 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _sync_conda_recipe(new_version: str) -> None:
    """Rewrite ``conda-recipe/meta.yaml`` to match the just-built
    release: bump version, regenerate sha256 from the local sdist
    (hatchling produces deterministic sdists, so this matches what
    PyPI will host), and rewrite the ``run:`` block from
    pyproject.toml's runtime deps.

    The recipe is the reference copy that lives in-repo; the canonical
    bioconda copy at ``bioconda-recipes/recipes/splicecraft/`` still
    needs a PR with the same bump, but having the local recipe always
    pre-bumped means the PR is a one-shot copy-paste.

    No-op if the recipe file is missing (older checkouts predate it).
    """
    if not CONDA_RECIPE.is_file():
        print(f"Note: {CONDA_RECIPE.relative_to(REPO_ROOT)} not present; "
              "skipping conda recipe sync.")
        return

    sdists = sorted((REPO_ROOT / "dist").glob("*.tar.gz"))
    if not sdists:
        _die("no sdist found in dist/ — conda sha256 cannot be computed.")
    if len(sdists) > 1:
        _die(f"expected exactly one sdist in dist/, found {len(sdists)}.")
    sdist = sdists[0]
    sha256 = _sha256_of(sdist)

    text = CONDA_RECIPE.read_text(encoding="utf-8")

    text, n = _CONDA_VERSION_RE.subn(rf'\1"{new_version}"\2', text)
    if n != 1:
        _die(f"conda recipe: expected 1 `{{% set version = ... %}}` "
             f"line, got {n}.")

    text, n = _CONDA_SHA256_RE.subn(rf"\g<1>{sha256}", text)
    if n != 1:
        _die(f"conda recipe: expected 1 sha256 line, got {n}.")

    pep_deps = _read_pyproject_runtime_deps()
    # `python >=3.10` is always the first run-dep — keep it pinned to
    # the value in pyproject.toml's `requires-python` rather than
    # hard-coding so a future Python-floor bump stays in sync. Parsing
    # the floor again here avoids a separate constant.
    py_floor_m = re.search(
        r'^requires-python\s*=\s*"\s*([^"]+?)\s*"',
        PYPROJECT.read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    py_floor = (
        py_floor_m.group(1).strip() if py_floor_m is not None else ">=3.10"
    )
    run_lines = [f"    - python {py_floor}"]
    run_lines.extend(f"    - {_pep440_to_conda(d)}" for d in pep_deps)
    new_run_block = "  run:\n" + "\n".join(run_lines) + "\n"

    text, n = _CONDA_RUN_BLOCK_RE.subn(new_run_block, text, count=1)
    if n != 1:
        _die("conda recipe: failed to rewrite the `run:` deps block "
             "(expected `  run:` followed by `    - …` lines).")

    CONDA_RECIPE.write_text(text, encoding="utf-8")
    print(f"  conda recipe → version {new_version}, sha256 {sha256[:16]}…, "
          f"{len(pep_deps)} runtime deps from pyproject")


def _wait_for_pypi(new_version: str, timeout_s: int = 360) -> bool:
    """Poll PyPI until the new sdist is downloadable. Bioconda's bot
    tries to fetch the sdist immediately when the PR opens; if PyPI's
    Trusted Publishing pipeline hasn't finished yet, the first CI run
    fails and the user has to re-trigger it manually. Polling here is
    a 1-3 minute wait that avoids the dance entirely.

    Returns True on success, False on timeout (in which case the
    caller opens the PR anyway and the bot's automatic retry will
    eventually pick up the sdist).
    """
    import time
    url = PYPI_SDIST_URL_TPL.format(version=new_version)
    deadline = time.time() + timeout_s
    print(f"  Polling {url}")
    while time.time() < deadline:
        result = subprocess.run(
            ["curl", "-sIfL", "-o", "/dev/null", url],
            capture_output=True,
        )
        if result.returncode == 0:
            print("  sdist live on PyPI.")
            return True
        time.sleep(10)
    print(f"  Timed out after {timeout_s}s. Bioconda PR will open "
          "anyway; their bot retries until the sdist is reachable.")
    return False


def _gh_logged_in_user() -> "str | None":
    """Return the GitHub username the local `gh` CLI is authenticated
    as, or None if `gh` isn't installed / authed (caller skips PR)."""
    if not shutil.which("gh"):
        return None
    result = subprocess.run(
        ["gh", "api", "user", "--jq", ".login"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return None
    name = result.stdout.strip()
    return name or None


def _ensure_bioconda_fork(owner: str) -> None:
    """Ensure `<owner>/bioconda-recipes` exists. Forks the upstream if
    it doesn't. Idempotent — `gh repo fork` is a no-op when the fork
    already exists."""
    result = subprocess.run(
        ["gh", "repo", "view", f"{owner}/bioconda-recipes"],
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        return
    print(f"  Forking {BIOCONDA_UPSTREAM} → {owner}/bioconda-recipes")
    _run(["gh", "repo", "fork", BIOCONDA_UPSTREAM,
          "--clone=false", "--default-branch-only"])


def _version_tuple(v: str) -> tuple[int, ...]:
    """Parse 'X.Y.Z[.suffix]' into a comparable tuple. Non-numeric
    suffix segments contribute their leading digits only (so 0.9.0rc1
    sorts as (0, 9, 0, 1)) — good enough for the supersede check in
    `_close_superseded_bioconda_prs`; the rich PEP-440 ordering isn't
    needed because we only compare same-tool branches."""
    parts: list[int] = []
    for piece in v.split("."):
        try:
            parts.append(int(piece))
        except ValueError:
            digits = "".join(ch for ch in piece if ch.isdigit())
            parts.append(int(digits) if digits else 0)
    return tuple(parts)


def _close_superseded_bioconda_prs(owner: str, new_version: str) -> None:
    """Close any open splicecraft PRs in bioconda-recipes whose branch
    targets a version OLDER than `new_version`. Same-version PRs are
    left alone — `_submit_bioconda_pr`'s existing check handles those.

    Why: bioconda reviewers are volunteers and a multi-PR stack for
    the same software with sequential versions creates ambiguity ("is
    one of these the canonical ask, or did the maintainer change
    their mind?"). Closing supersedes with an explicit pointer makes
    the queue obvious. Without this guard, every `release.py` invocation
    piles on a new PR while leaving prior ones open — exactly what
    happened during the 0.8.9 → 0.9.2 stretch where four "first
    submission" PRs sat green for days while the maintainer waited
    for the queue to clarify itself.

    Branch convention: `splicecraft-<version>` (see
    `_submit_bioconda_pr` below — both producer + this consumer).
    Failures are non-fatal (`check=False`); the new PR opens regardless.
    """
    import json as _json
    result = subprocess.run(
        ["gh", "pr", "list",
         "--repo", BIOCONDA_UPSTREAM,
         "--author", owner,
         "--state", "open",
         "--search", "splicecraft in:title",
         "--json", "number,headRefName,url",
         "--limit", "20"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return
    try:
        prs = _json.loads(result.stdout)
    except _json.JSONDecodeError:
        return
    if not isinstance(prs, list):
        return
    new_tuple = _version_tuple(new_version)
    for pr in prs:
        branch = pr.get("headRefName", "") or ""
        if not branch.startswith("splicecraft-"):
            continue
        pr_version = branch[len("splicecraft-"):]
        if _version_tuple(pr_version) >= new_tuple:
            # Equal or newer — leave alone. Same-version is handled
            # by the existing branch-match check inside
            # `_submit_bioconda_pr`; newer would mean something
            # unexpected (manual PR ahead of release.py) and we'd
            # rather not auto-close that.
            continue
        print(f"  Closing superseded bioconda PR #{pr['number']} "
              f"(v{pr_version} < v{new_version})")
        comment_body = (
            f"Superseded by the upcoming PR for v{new_version}. "
            f"Closing to clear the queue so a maintainer only needs "
            f"to review the latest — the recipe content for each "
            f"interim version is functionally identical apart from "
            f"version + sha256."
        )
        # Comment THEN close, both `check=False` — leaving them
        # non-fatal so a transient gh / API hiccup doesn't block the
        # release. Worst case: an older PR stays open and we manually
        # close it later (the same situation as before this guard).
        subprocess.run(
            ["gh", "pr", "comment", str(pr["number"]),
             "--repo", BIOCONDA_UPSTREAM,
             "--body", comment_body],
            check=False, capture_output=True,
        )
        subprocess.run(
            ["gh", "pr", "close", str(pr["number"]),
             "--repo", BIOCONDA_UPSTREAM],
            check=False, capture_output=True,
        )


def _submit_bioconda_pr(new_version: str) -> None:
    """Open (or update) a bioconda PR with the current recipe.

    Steps:
      1. Verify `gh` is authed; skip with a note if not.
      2. Fork bioconda/bioconda-recipes if the maintainer's account
         doesn't already host a fork.
      3. Close any open splicecraft PRs at an older version so the
         queue surfaces only the canonical ask
         (`_close_superseded_bioconda_prs`).
      4. Shallow-clone the fork to a temp dir, fast-forward its
         master to upstream (so our branch is fresh).
      5. Drop in `recipes/splicecraft/meta.yaml` from our in-repo
         recipe (which `_sync_conda_recipe` just refreshed).
      6. Branch, commit, force-push to the fork (force-with-lease so
         a stale branch from a previous failed attempt doesn't block
         us; per-version branches make collisions impossible across
         releases).
      7. Open the PR against the upstream master. If one for this
         version already exists (re-run), print the URL and skip.

    First-time submissions (no `recipes/splicecraft/` directory in
    bioconda yet) and updates use the same flow; bioconda's bot
    handles both transparently.
    """
    if not CONDA_RECIPE.is_file():
        print(f"Note: {CONDA_RECIPE.relative_to(REPO_ROOT)} not present; "
              "skipping bioconda PR.")
        return

    owner = _gh_logged_in_user()
    if owner is None:
        print("Note: `gh` CLI not installed or not authenticated; "
              "skipping bioconda PR. Run `gh auth login` and rerun "
              "with --bioconda-only to publish to bioconda.")
        return

    print(f"  gh authenticated as {owner}")
    _ensure_bioconda_fork(owner)
    _close_superseded_bioconda_prs(owner, new_version)

    import tempfile
    with tempfile.TemporaryDirectory(prefix="bioconda-pr-") as tmp:
        clone = Path(tmp) / "bioconda-recipes"
        # Shallow clone keeps the bioconda-recipes mega-repo manageable
        # (it has 10K+ recipes and a multi-GB history).
        _run([
            "git", "clone", "--depth=1", "--single-branch",
            "--branch=master",
            f"https://github.com/{owner}/bioconda-recipes.git",
            str(clone),
        ])
        # Pull the latest upstream master into the shallow clone so
        # our branch starts from current head, not the fork's
        # potentially-stale snapshot.
        _run(["git", "-C", str(clone), "remote", "add",
              "upstream", f"https://github.com/{BIOCONDA_UPSTREAM}.git"])
        _run(["git", "-C", str(clone), "fetch", "--depth=1",
              "upstream", "master"])
        _run(["git", "-C", str(clone), "reset", "--hard",
              "upstream/master"])

        recipe_dir = clone / "recipes" / "splicecraft"
        is_first_submission = not recipe_dir.exists()
        recipe_dir.mkdir(parents=True, exist_ok=True)
        target_meta = recipe_dir / "meta.yaml"
        shutil.copy2(CONDA_RECIPE, target_meta)

        branch = f"splicecraft-{new_version}"
        _run(["git", "-C", str(clone), "checkout", "-b", branch])
        _run(["git", "-C", str(clone), "add",
              "recipes/splicecraft/meta.yaml"])

        commit_msg = (
            f"splicecraft: add v{new_version} (first submission)"
            if is_first_submission
            else f"splicecraft: bump to v{new_version}"
        )
        # Configure git identity for the commit in the temp clone so
        # the commit doesn't fail on a clean machine without a global
        # gitconfig. Use the same identity gh is authed as.
        _run(["git", "-C", str(clone), "config", "user.name",  owner])
        _run(["git", "-C", str(clone), "config", "user.email",
              f"{owner}@users.noreply.github.com"])
        _run(["git", "-C", str(clone), "commit", "-m", commit_msg])
        # `--force-with-lease` rather than `--force`: safer if someone
        # else is collaborating on the same branch (won't happen for
        # auto-named per-version branches, but the habit is cheap).
        _run(["git", "-C", str(clone), "push",
              "--force-with-lease", "origin", branch])

        # Check whether a PR for this version already exists. If yes,
        # skip the create call and just print the URL — a re-run
        # shouldn't open duplicate PRs.
        pr_check = subprocess.run(
            ["gh", "pr", "list",
             "--repo", BIOCONDA_UPSTREAM,
             "--head", f"{owner}:{branch}",
             "--state", "open",
             "--json", "url,number",
             "--limit", "1"],
            capture_output=True, text=True, check=False,
        )
        if pr_check.returncode == 0 and pr_check.stdout.strip() not in (
            "", "[]"
        ):
            # A PR already exists for this branch — surface it.
            import json as _json
            try:
                existing = _json.loads(pr_check.stdout)
                if existing:
                    print(f"  Existing PR: {existing[0]['url']}")
                    return
            except _json.JSONDecodeError:
                pass

        pr_title = commit_msg
        pr_body = (
            f"Updates the SpliceCraft recipe to v{new_version}.\n\n"
            f"- PyPI: https://pypi.org/project/splicecraft/{new_version}/\n"
            f"- Upstream: https://github.com/Binomica-Labs/SpliceCraft\n"
            f"- Tag: https://github.com/Binomica-Labs/SpliceCraft/"
            f"releases/tag/v{new_version}\n\n"
            "Recipe synced from the in-repo `conda-recipe/meta.yaml` "
            "via `release.py`: version + sha256 (regenerated from the "
            "live sdist) + runtime deps (rewritten from `pyproject.toml`'s "
            "`[project] dependencies`).\n"
        )
        if is_first_submission:
            pr_body += (
                "\n_This is the first bioconda submission for "
                "splicecraft._ Recipe maintainer: @"
                + owner
                + ". Happy to address any bot feedback."
            )

        create_result = subprocess.run(
            ["gh", "pr", "create",
             "--repo", BIOCONDA_UPSTREAM,
             "--head", f"{owner}:{branch}",
             "--base", "master",
             "--title", pr_title,
             "--body", pr_body],
            check=False,
        )
        if create_result.returncode != 0:
            print("  Note: `gh pr create` exited non-zero. The fork "
                  "branch is pushed; open the PR manually at "
                  f"https://github.com/{BIOCONDA_UPSTREAM}/compare/"
                  f"master...{owner}:bioconda-recipes:{branch}")


def _current_version() -> str:
    """Read the canonical version from pyproject.toml. Used by the
    `--bioconda-only` re-submit path which skips the bump step."""
    text = PYPROJECT.read_text(encoding="utf-8")
    m = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if m is None:
        _die("could not read current version from pyproject.toml.")
    return m.group(1)


def _create_github_release(version: str) -> None:
    """Create a GitHub Release for ``v<version>`` — the matching
    ``CHANGELOG.md`` section as the body, plus the built wheel + sdist
    from ``dist/`` as downloadable assets.

    Optional + NON-FATAL: the tag is already pushed and PyPI is already
    publishing by the time we get here, so a Release is a nice-to-have
    on top. Skips with a printed note (never aborts) when ``gh`` is
    missing / unauthenticated, or when the Release already exists — so a
    box without GitHub CLI still ships to PyPI exactly as before.

    One Release per version tag is the intended GitHub workflow (it
    surfaces the changelog on the Releases page + notifies watchers);
    this just stops them being created by hand.
    """
    import shutil
    tag = f"v{version}"
    if shutil.which("gh") is None:
        print(f"  gh CLI not found — skipping GitHub Release for {tag} "
              f"(create later with: gh release create {tag}).")
        return
    if subprocess.run(["gh", "auth", "status"],
                      capture_output=True, text=True).returncode != 0:
        print(f"  gh not authenticated — skipping GitHub Release for {tag}.")
        return
    # Idempotent: a re-run (or a hand-made Release) is left untouched
    # rather than erroring out the already-completed release.
    if subprocess.run(["gh", "release", "view", tag],
                      capture_output=True, text=True).returncode == 0:
        print(f"  GitHub Release {tag} already exists — leaving it.")
        return
    # Body = the `## [<version>]` changelog section, minus a trailing
    # `---` separator so the Release doesn't end on a stray rule.
    notes = re.sub(r"\n*-{3,}\s*$", "",
                   _changelog_section_body(version).strip()).strip()
    if not notes:
        notes = f"SpliceCraft {tag}."
    dist = REPO_ROOT / "dist"
    assets = (sorted(str(p) for p in dist.iterdir() if p.is_file())
              if dist.is_dir() else [])
    result = subprocess.run(
        ["gh", "release", "create", tag,
         "--title", f"SpliceCraft {tag}",
         "--notes", notes, "--latest", *assets],
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        url = (result.stdout.strip().splitlines() or [""])[-1].strip()
        print(f"  GitHub Release {tag} created"
              + (f" (+{len(assets)} asset(s))" if assets else "")
              + (f": {url}" if url else "."))
    else:
        print("  GitHub Release create failed (non-fatal — the tag is "
              f"pushed + PyPI is publishing): "
              f"{result.stderr.strip() or result.returncode}")
        print(f"  Create it later with: gh release create {tag} …")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Cut a new SpliceCraft release.",
    )
    parser.add_argument(
        "version",
        nargs="?",
        help="New version (X.Y.Z, optionally with PEP-440 suffix). "
             "Omit when using --bioconda-only.",
    )
    parser.add_argument(
        "--bioconda-only",
        action="store_true",
        help="Skip every PyPI step and just open / refresh the "
             "bioconda PR against the upstream `bioconda-recipes` "
             "repo using the current `conda-recipe/meta.yaml`. Use "
             "after a PyPI release that pre-dated this automation, "
             "or to retry the PR after addressing bot feedback.",
    )
    args = parser.parse_args(argv)

    if args.bioconda_only:
        if args.version is not None:
            _die("--bioconda-only does not accept a version argument; "
                 "it operates on the recipe currently in the tree.")
        new_version = _current_version()
        _heading(f"Bioconda-only re-submit for v{new_version}")
        _wait_for_pypi(new_version)
        _heading("Opening / refreshing bioconda PR")
        _submit_bioconda_pr(new_version)
        return 0

    if args.version is None:
        _die("version argument required (or pass --bioconda-only).")
    new_version = args.version

    if not _VERSION_RE.match(new_version):
        _die(f"version must look like X.Y.Z (got {new_version!r}).")

    _summarize_pending_changes()
    _ensure_tag_unused(new_version)
    _ensure_changelog_entry(new_version)

    _heading(f"Bumping version to {new_version}")
    _bump_version_in_file(
        PYPROJECT, _PYPROJECT_VERSION_RE, new_version, "pyproject.toml",
    )
    _bump_version_in_file(
        SPLICECRAFT, _SPLICECRAFT_VERSION_RE, new_version, "splicecraft.py",
    )
    _verify_bump(PYPROJECT,   new_version, "version")
    _verify_bump(SPLICECRAFT, new_version, "__version__")

    # Run ruff BEFORE pytest — sweep #16 (2026-05-21) added this so a
    # lint-failing commit can't slip through release and turn the CI
    # badge red. The same `ruff check` runs in `.github/workflows/
    # test.yml`'s `lint` job; matching the local invocation means a
    # release that passes here also passes CI's lint gate. Pyright is
    # NOT bundled into release.py because it shells out to a separate
    # JS-pinned binary that's slow + flaky offline — it's a CI-only
    # gate (`PYRIGHT_PYTHON_FORCE_VERSION=latest` in test.yml).
    _heading("Running ruff lint")
    _run(["ruff", "check", "."])

    _heading("Running test suite")
    # Parallel via pytest-xdist; previously serial took ~13 min, -n auto
    # cuts that to ~5 min on an 8-core box. Tests are isolated by the
    # autouse `_protect_user_data` fixture so cross-worker collisions
    # are impossible.
    _run([sys.executable, "-m", "pytest", "-n", "auto", "-q", "--tb=short"])

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

    # 2026-05-27: conda recipe sync + bioconda PR removed from the
    # default release flow at user request — too-frequent bioconda
    # submissions were drawing reviewer complaints. The recipe in
    # `conda-recipe/meta.yaml` stays in the tree for the rare manual
    # re-submission, but is no longer touched per release. The
    # `--bioconda-only` flag still works for explicit re-submission
    # via `./release.py --bioconda-only`.

    _heading("Committing + tagging + pushing")
    # Sweep #18 (2026-05-21): bundle ALL accumulated working-tree
    # changes (tracked + new untracked) into the release commit, not
    # just the version-bump files. The user's workflow is "iterate
    # in the working tree, only commit at release time", so the
    # release commit needs to pick up every modification + every
    # new file the user authored since the previous tag. `git add
    # -A` is bounded by `.gitignore` so build artifacts, caches,
    # logs, etc. stay out.
    _run(["git", "add", "-A"])
    _run(["git", "commit", "-m", f"Release v{new_version}"])
    _run(["git", "tag", f"v{new_version}"])
    _run(["git", "push", "origin", "master"])
    _run(["git", "push", "origin", f"v{new_version}"])

    # One GitHub Release per version tag (changelog notes + dist assets).
    # Optional + non-fatal — see `_create_github_release`.
    _heading("Creating GitHub Release")
    _create_github_release(new_version)

    print()
    print("═" * 61)
    print(f" Release v{new_version} pushed.")
    print(" GitHub Actions will publish to PyPI in ~2 minutes.")
    print(" Watch:   https://github.com/Binomica-Labs/SpliceCraft/actions")
    print(" Verify:  https://pypi.org/project/splicecraft/")
    print(f" Release: https://github.com/Binomica-Labs/SpliceCraft/"
          f"releases/tag/v{new_version}")
    print(" (Bioconda PR step skipped — run `./release.py "
          "--bioconda-only` to re-submit the recipe manually.)")
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
