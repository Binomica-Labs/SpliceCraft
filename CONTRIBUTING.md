# Contributing to SpliceCraft

Thanks for considering a contribution. SpliceCraft is built and used
by a practicing bioengineer; the bar is "would a senior reviewer
merge this on the first look." This document walks through what that
means in practice.

> If you're an AI agent landing here, read `CLAUDE.md` first.
> The 10 numbered sacred invariants are not negotiable, and most of
> the friction-causing edge cases are documented there (and in the
> `[PIT-]`/`[INV-]` catalog in `docs/invariants.md`).

## Table of contents

- [Before you open a PR](#before-you-open-a-pr)
- [Local setup](#local-setup)
- [Hub and layered siblings](#hub-and-layered-siblings)
- [Sacred invariants](#sacred-invariants)
- [Testing](#testing)
- [Style](#style)
- [Logging](#logging)
- [Security-sensitive code](#security-sensitive-code)
- [Documentation](#documentation)
- [Commits, releases, and version bumps](#commits-releases-and-version-bumps)
- [Code review](#code-review)

## Before you open a PR

1. **Read `CLAUDE.md`.** It is the project handoff document and the
   source of truth for invariants, conventions, and known pitfalls.
2. **Read the related test file.** Tests are 1:1 named after the
   subsystem they cover (`test_dna_sanity.py`, `test_commercialsaas_io.py`,
   `test_agent_api.py`, etc.).
3. **For non-trivial changes, open an issue first.** Five-line bug
   fixes can go straight to PR. Anything that touches biology
   primitives, persistence, or the agent-API benefits from a design
   conversation up front.

## Local setup

```bash
git clone https://github.com/Binomica-Labs/SpliceCraft.git
cd SpliceCraft
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

That installs the runtime deps plus `pytest`, `pytest-xdist`,
`pytest-asyncio`, `hypothesis`, `build`, and `twine`. The CI image
adds `ruff`, `pyright`, and `coverage` (already in `[dev]` extras);
to mirror CI's `pyright` exactly, ALSO export the latest-pyright
pin in your shell:

```bash
export PYRIGHT_PYTHON_FORCE_VERSION=latest
ruff check .
pyright
```

The `pyright` Python package is a launcher that downloads the
actual binary at first run; without `PYRIGHT_PYTHON_FORCE_VERSION`
the launcher caches whatever release was current when it was
installed, which drifts away from CI's version over time. CI sets
the env var explicitly (`.github/workflows/test.yml`).

`release.py` goes one step further: right before the lint gate it
runs `pip install --upgrade pyright`, so the release is cleared by
the current checker and the launcher stops printing "a new pyright
version is available" on every run — a warning you learn to skim
past is a warning that hides the next real one. The upgrade is
non-fatal (offline releases still work) and skippable with
`SPLICECRAFT_SKIP_PYRIGHT_UPGRADE=1`.

Inner-loop sanity check (under 2 seconds):

```bash
python3 -m pytest tests/test_dna_sanity.py
```

Full suite (~5 minutes on 8 cores):

```bash
python3 -m pytest -n auto -q
```

## Hub and layered siblings

The application is a **hub + layered-siblings** layout: `splicecraft.py`
is the ~99k-line hub (the application core — `PlasmidApp`, the big-3
panels, the data-safety policy, the `.dna` blob store, the agent-API
server, and the app-coupled modals), and the cleanly-separable layers
live in 25 flat `splicecraft_*.py` siblings the hub imports and
re-exports — including the stdlib-only `splicecraft_cli.py`
agent-client sidecar. The hub stays greppable as one totally-ordered file; each
sibling is a bounded, independently-loadable unit so a context-limited
model can hold one whole. The Phase-D modularization that carved the
~131k-line monolith down to ~99k is essentially complete — see
`docs/architecture.md` and `CLAUDE.md` for the full sibling map.

**A new sibling — or a move into an existing one — must pass three guards:**

1. **No upward imports / no cycle, and it's packaged.**
   `tests/test_import_layers.py` classifies every sibling into a layer
   (L0→L7), forbids upward imports, and verifies the file is in
   `pyproject.toml`'s wheel **and** sdist lists (miss either and the
   released wheel breaks at the re-import step).
2. **The public surface stays byte-for-byte identical.** Re-export every
   moved name so `import splicecraft as sc; sc.<name>` resolves
   unchanged; `tests/public_surface_baseline.json` enforces it (only a
   deliberate migration to `_state` drops a name).
3. **Shared mutable state lives in `_state`**, accessed by attribute
   (`_state.X`, never `from splicecraft_state import X` — that binds a
   stale copy). `tests/test_state_module.py` guards the
   single-source-of-truth contract. A sibling reaches hub-pinned
   dependencies through `_state` hooks, never by importing the hub
   (which would be a cycle).

Everything else must still pass without modification, including the
byte-for-byte assertions in `test_commercialsaas_io.py`. What stays
hub-side: `PlasmidApp`, the big-3 panels (`PlasmidMap` / `SequencePanel`
/ `LibraryPanel`), their anchored modals, and the blob writers — the
God-class application core, which would need in-App decomposition rather
than a lift.

## Sacred invariants

`CLAUDE.md` lists 10 numbered sacred invariants (biology correctness +
the save chokepoint); the wider catalog — known pitfalls, persistence /
concurrency / UI contracts — lives in `docs/invariants.md` as
`[PIT-]`/`[INV-]` tags. Each has at least one regression test. The
classes most often touched:

- **Biology correctness** (sacred #1–#6, #8, #9): palindromic-enzyme
  scanning, IUPAC reverse-complement, wrap-feature math, restriction-
  cloning, Type IIS digest, exact-match annotation transfer.
- **Data safety** (sacred #7, #10 + the `[INV-]` data-safety tags):
  atomic saves, the deepcopy-on-read+save cache contract, four-layer
  backup, symlink refusal, pre-update snapshots.
- **Concurrency & UI/UX** (the `[PIT-]`/`[INV-]` catalog): worker
  cancellation contract, RLock on saves, stale-record guard, undo
  deepcopy, modal Ctrl+Z opt-out, natural-sort row-mapping symmetry.

**Touching invariant code without updating its test, or weakening
an invariant without explicit discussion, is grounds for the PR
being closed.**

## Testing

| Suite                                  | Use                                              | Runtime |
|----------------------------------------|--------------------------------------------------|---------|
| `tests/test_dna_sanity.py`             | Inner loop while iterating on biology code       | < 2 s   |
| `tests/test_commercialsaas_io.py`      | When touching `.dna` reader / writer             | ~30 s   |
| `tests/test_agent_api.py`              | When touching `_h_*` endpoints                   | ~45 s   |
| `tests/test_smoke.py`                  | End-to-end + update / restore flows              | ~90 s   |
| `python3 -m pytest -n auto -q`         | Before commit, before release                    | ~5 min  |

**Conventions:**

- Cross-validate biology against Biopython where applicable.
- No network calls from tests. Use the autouse `_protect_user_data`
  fixture in `tests/conftest.py` — it monkeypatches every `_*_FILE`
  path, so production data is never touched.
- Async tests: `async with app.run_test(size=...)` with a double
  `await pilot.pause()` for `call_after_refresh` to settle.
- Regression guards cite the date in their docstring
  (`# Regression guard for 2026-MM-DD fix`).
- Add new persisted libraries to `_protect_user_data` in
  `tests/conftest.py` AND to `_check_data_files`. Cover corruption
  recovery in `test_data_safety.py`.
- New modals add a row to `test_modal_boundaries.py::_MODAL_CASES`
  (every modal must fit 160×48).

## Style

- Python 3.10+ syntax. `from __future__ import annotations` is fine.
- Type hints encouraged but not strictly enforced (86%+ today).
  Strict pyright is enabled on `splicecraft.py` and
  `splicecraft_cli.py` — tests are excluded.
- **No bare `except`.** Use narrow exception types
  (`NoMatches`, `ET.ParseError`, `(OSError, json.JSONDecodeError)`,
  etc.). The one exception: `except Exception as exc` is reserved
  for `@work` thread bodies, and must `_log.exception(...)` the
  failure.
- **User-facing errors:** `self.notify(...)` or
  `Static.update("[red]...[/]")`. Never raw tracebacks. Diagnostic
  detail goes to `_log.exception`.
- Comments are for *why*, not *what*. Don't restate code in prose.
- Don't add scaffolding for hypothetical future requirements.

## Logging

Structured event logging (invariant #42) is how field bugs get
triaged. Use `_log_event(event, **fields)` for any state transition
worth a log; use `@_action_log("name")` on `action_*` methods; use
`@_timed("path")` for heavy ops.

Event-name convention:

- `app.<area>.<verb>` — user actions
- `op.<area>.<verb>` — operations
- `<noun>.<verb>` — state transitions (`save.ok`, `record.loaded`,
  `migration.step`)

**Never log sequence content.** `_repr_for_log` exists to truncate
and tag any value that might carry DNA bases.

## Security-sensitive code

If your change touches any of the surfaces listed in `SECURITY.md`:

- XML parsing → must route through `_safe_xml_parse`.
- New JSON file → must save through `_safe_save_json` (atomic,
  backed up, snapshot-eligible).
- New HTTP fetch → must follow the `resp.read(MAX + 1)` + bail
  pattern; add a `_*_MAX_RESPONSE_BYTES` constant.
- New agent-API endpoint → must use `_check_agent_read_dir` /
  `_check_agent_write_path` for any filesystem touch; must use
  `_agent_save_or_500(...)` for any save call.
- New path-derived filename → must sanitize via `_sanitize_path` or
  `_dna_sidecar_path`.

When in doubt, copy the pattern from an adjacent endpoint.

## Documentation

- **User-visible changes update the README in the same PR.** Stale
  README is a release-blocker.
- **New invariants get a numbered entry in `CLAUDE.md`.**
- **Changelog entry in `CHANGELOG.md`** under the unreleased section
  (`release.py` rotates it on tag).
- **Public endpoints get docstrings.** Internal helpers — only if the
  *why* is non-obvious.

## Commits, releases, and version bumps

- One concept per commit. `release.py` rotates the changelog and
  builds atomically, so commits inside a release window are easier
  to bisect when they're focused.
- Version bumps are usually patch (`Z` in `X.Y.Z`). Minor bumps
  signal a UX-visible change; major bumps are reserved for the
  v1.0.0 milestone (see `V1_GATE.md`).
- **Do not push or release without explicit maintainer approval.**
  Contributors land via PR; the maintainer cuts the release.

## Code review

What a maintainer looks for, in order:

1. Does it break any of the 41 invariants?
2. Does it pass `python3 -m pytest -n auto -q`?
3. Does it pass `pyright` strict?
4. Does it pass `ruff check`?
5. Is the diff focused? (Cleanup belongs in its own commit.)
6. Does the README still reflect reality?
7. Does the commit message explain *why*?
8. Is there a regression test for any bug being fixed?

Thanks again. Quality contributions make the next user's bench day
better.
