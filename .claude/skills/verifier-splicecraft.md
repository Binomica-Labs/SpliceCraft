---
name: verifier-splicecraft
description: Sandboxed evidence-capture harness for verifying SpliceCraft changes. Drives the Textual TUI via Pilot, captures SVG screenshots, and pins the data dir away from the user's real `~/.local/share/splicecraft/` so a probe can never touch live state. Use this when verifying any change that reaches the TUI, an agent endpoint, or a persisted-data file.
---

# Verifier — SpliceCraft

The surface is the **Textual TUI** and the **agent-API endpoints**.
Tests cover helper functions; this verifier covers what the user (or
an agent) actually drives through the real event loop.

## Sacred: sandbox the data dir BEFORE the import

SpliceCraft's `_DATA_DIR` is set at module-import time from
`XDG_DATA_HOME` (or `~/.local/share/`). If you `import splicecraft`
without first redirecting `XDG_DATA_HOME`, the probe touches the
user's real plasmid library, primer collection, and `.bak` history.

A probe that touches real data is **not a probe**, it's a malfunction.

```python
import os, sys, tempfile
_sandbox = tempfile.mkdtemp(prefix="sc-verify-")
os.environ["XDG_DATA_HOME"] = _sandbox
os.environ.setdefault("SPLICECRAFT_SKIP_LOCK", "1")
sys.path.insert(0, "/home/seb/SpliceCraft")
import splicecraft as sc

# `_DATA_DIR` lives in `_state` since the modularization — `sc._DATA_DIR`
# raises AttributeError, so an assert written against it never guards
# anything. Read it through `_state`, which is the single source of truth.
assert "sc-verify-" in str(sc._state._DATA_DIR), \
    f"sandbox not active — _DATA_DIR={sc._state._DATA_DIR}"

# L2 chokepoint (2026-05-22): authorise writes against the sandbox.
# Without this, every `_save_*` raises RuntimeError. The helper
# itself refuses if the data dir isn't under the OS tempdir, so
# you can't bypass it without first sandboxing.
sc._authorize_writes_for_sandbox(sc._state._DATA_DIR)
```

After the probe finishes, `shutil.rmtree(_sandbox, ignore_errors=True)`
or just leave it for the OS to clean — they're in `/tmp` either way.

### Seeding the sandbox: JSON only, NEVER the blob store

An empty sandbox is the default and the right one for most probes. When a
probe genuinely needs the real library to reproduce something, copy the
**JSON files only** — `collections.json`, `library.json`, `primers.json`,
`parts_bin.json`, `features.json`, … — and leave the blobs where they are.

`plasmid_blobs/` and `dna_originals/` are hundreds of megabytes across
thousands of files. `/tmp` is a **7.5 GB tmpfs**, i.e. RAM: a handful of
probe runs that each copy the blob store fill it, and the failure is not a
clean "disk full" but `EDQUOT` out of unrelated commands — the whole session
stops working (field report 2026-11-08). Blobs are also content-addressed and
read-only in every read path, so copying them buys nothing.

```python
import shutil, pathlib
_real = pathlib.Path.home() / ".local/share/splicecraft"
_dst  = pathlib.Path(_sandbox) / "splicecraft"
_dst.mkdir(parents=True, exist_ok=True)
for _f in _real.glob("*.json"):          # JSON ONLY — no dirs, no .bak
    shutil.copy2(_f, _dst / _f.name)
```

A probe that needs real blob CONTENT should point the app at the real blob
dir read-only (`sc._state._DNA_ORIGINALS_DIR`) rather than copy it, and must
not write there — the L2 chokepoint refuses writes outside the sandbox
anyway. If a probe run is interrupted, delete **only your own**
`/tmp/sc-verify-*` directories; other agents and other projects keep scratch
under `/tmp` too.

## Drive the app

Two surfaces — pick the one that reaches your change.

### Surface A — Pilot (Textual's headless event loop)

Pilot is the same event-loop Textual uses for real keystrokes; it's
not a unit-test mock. Anything queryable in production is queryable
here. Use it for:
- Menu / modal / action verification
- Multi-screen navigation flows
- Anything that hits the reactive system

```python
import asyncio

async def probe():
    sc.PlasmidApp._skip_seed = True
    sc.PlasmidApp._skip_snapshot = True
    sc.PlasmidApp._skip_update_check = True
    sc.PlasmidApp._preload_demo_record = None  # NOT False — that's a bool

    app = sc.PlasmidApp()
    async with app.run_test(size=(200, 50)) as pilot:
        # Settle initial mount. Nested compose paths (e.g. the
        # SequencingScreen → TabbedContent → TabPane → sub-TabbedContent
        # → TabPane structure) need 4-6 pauses, not 2.
        for _ in range(6):
            await pilot.pause()

        # Dismiss any startup modal — What's New surfaces on a fresh
        # data dir; it eats the first keypress otherwise.
        while len(app.screen_stack) > 1:
            app.pop_screen()
            for _ in range(2):
                await pilot.pause()

        # ── Your probe steps here ──────────────────────────────────
        # Action invocation (skips menu — use when testing the action):
        app.action_open_sequencing()
        for _ in range(6):
            await pilot.pause()
        assert type(app.screen).__name__ == "SequencingScreen"

        # Keyboard binding (use when testing the binding itself):
        await pilot.press("alt+x")
        for _ in range(6):
            await pilot.pause()
        assert type(app.screen).__name__ == "ExperimentProjectsPickerModal"

asyncio.run(probe())
```

**Settle rules:**
- 4 pauses: simple screen push
- 6 pauses: TabbedContent with sub-tabs (Sequencing, Synthesis)
- 6+ pauses: anything that fires a `@work(thread=True)` worker

**Query rules:**
- `app.screen` returns the topmost screen — always wait pauses before
  querying it
- `app.screen.query_one("#widget-id", WidgetType)` — same syntax as
  production code
- Inner `TabPane`s mount AFTER the outer `TabbedContent`, so flip
  `tabs.active = "tab-id"` then pause again before querying inside

### Surface B — tmux (real terminal output)

Use when you need to see actual terminal rendering, key sequences,
or a specific terminal capability quirk.

```bash
tmux -L sc-verify new-session -d -s sc -x 200 -y 50 'XDG_DATA_HOME=$(mktemp -d) SPLICECRAFT_SKIP_LOCK=1 python3 /home/seb/SpliceCraft/splicecraft.py'
sleep 4
tmux -L sc-verify send-keys -t sc Enter   # dismiss splash
sleep 2
tmux -L sc-verify capture-pane -t sc -p > /tmp/sc-pane.txt

# Send Alt+letter via xterm escape (Textual reads these as M-letter):
tmux -L sc-verify send-keys -t sc M-f     # Alt+F → File menu

# Cleanup:
tmux -L sc-verify kill-server
```

The `-L sc-verify` socket name isolates from any other tmux session.

## Capture evidence

```python
# Per-touchpoint SVG. Textual's export_screenshot returns SVG text.
with open("/tmp/sc-snap-<touchpoint>.svg", "w") as f:
    f.write(app.export_screenshot(title="<touchpoint label>"))
```

Reviewers open the SVG in a browser; SpliceCraft's braille-dot maps
render correctly because Textual's screenshot includes the full font
information.

For tmux verifies, `capture-pane -p > /tmp/<file>.txt`.

## Agent-endpoint quirk

`_h_load_file` and similar agent handlers call `app.call_from_thread`
which **must run from a non-app thread**. They work via the HTTP
server in production. From Pilot (single-thread test loop) they raise
`RuntimeError: call_from_thread must run in a different thread`.

**Workaround:** drive the underlying helper directly. The agent
handler is mostly payload validation + dispatch; the actual work is
done by helpers like `_gff3_path_to_record`, `_fasta_path_to_record`,
`_ab1_path_to_record`. Tests in `tests/test_agent_api.py` cover the
HTTP path; your verifier covers the helper.

```python
# DON'T (raises RuntimeError):
res = sc._h_load_file(app, {"path": "/tmp/x.gff3"})

# DO:
rec = sc._gff3_path_to_record("/tmp/x.gff3")
```

## What FAIL looks like in this codebase

- `AttributeError: 'bool' object has no attribute 'id'` from `_load_demo`
  → you set `_preload_demo_record = False` instead of `None`. The
  attr expects a SeqRecord or None.
- `No nodes match '#<id>'` → you queried before the screen finished
  mounting. Add more `pilot.pause()` calls.
- Modal didn't dismiss what you expected → a startup modal (What's New,
  Update Available, Splash) is on top. Loop `pop_screen` until
  `screen_stack == 1`.
- The screen looks right but the user's real `~/.local/share/splicecraft`
  shows up in the trace → **sandbox failed**. Stop, fix the
  `XDG_DATA_HOME` order, never proceed.
- Test suite passes but UI broken → the change touches mount-time
  code paths the tests skip (e.g., `_skip_snapshot = True` in tests
  flips False in `main()`).

## Useful entry points

| Touchpoint | Action / binding | Verify by |
|---|---|---|
| File → Open | `action_open_file()` / Ctrl+O | `app.screen` == `OpenFileModal` |
| File → Save | `action_save()` / Ctrl+S | dirty flag clears, source path set |
| File → Export collection (bulk) | `action_export_collection()` | `app.screen` == `BulkExportCollectionModal` |
| Sequencing → Sanger tab | `action_open_sequencing()` + flip `#seq-tabs.active = "seq-tab-sanger"` | `#sanger-tree` mounts |
| Synthesis | `action_open_synthesis()` / Alt+Y | `app.screen` == `SynthesisScreen` |
| Experiments | `action_open_experiments()` / Alt+X | `app.screen` == `ExperimentProjectsPickerModal` |
| Settings → Restore from backup | `app.push_screen(sc.RestoreFromBackupModal())` | `#restore-target` Select mounts |
| Agent endpoint (read) | `sc._h_<name>(app, payload)` | returns dict on success, `(dict, code)` tuple on error |
| Agent endpoint (write) | drive the underlying helper directly | see "Agent-endpoint quirk" above |

## Pre-flight before reporting

- [ ] Sandboxed (`assert "sc-verify-" in str(sc._state._DATA_DIR)`)
- [ ] No `~/.local/share/splicecraft/` paths in captured output
- [ ] Each touchpoint has at least one PASS/FAIL line with the actual
      screen / widget name as evidence
- [ ] At least one 🔍 probe (off the happy path) per non-trivial change
- [ ] SVG / pane captures saved to `/tmp/sc-snap-*` or referenced
      inline

Then write the report per the parent `verify` skill's format.
