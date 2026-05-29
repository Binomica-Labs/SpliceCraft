# OS + Terminal Compatibility Matrix

SpliceCraft is a Textual TUI; it inherits Textual's cross-platform
support and adds a few platform-specific extras (data-dir locking,
clipboard, signal handling, file paths). This page lists every
tested OS × terminal pairing, known issues, and the env-var fixes
the launcher recommends when a capability is missing.

Last update: 2026-05-19 (SpliceCraft 0.9.7+).

---

## Quick reference

| OS / WSL | Terminal | Status | Notes |
|---|---|---|---|
| **Linux** | gnome-terminal | ✅ Fully supported | Reference platform; CI runs on Ubuntu |
| Linux | kitty | ✅ Fully supported | True-color, mouse, OSC 52 all work |
| Linux | alacritty | ✅ Fully supported | Same as kitty |
| Linux | xterm | ⚠ Limited | OSC 52 clipboard needs `allowWindowOps`; mouse drag-select may need `xtermMouseProtocol` set |
| Linux | tmux / screen | ✅ Works under | Pass `-e` (env), use `tmux`'s `set -g allow-passthrough on` for OSC 52 |
| **macOS** | Terminal.app | ✅ Fully supported | macOS 11+; older macOS may lack true-color |
| macOS | iTerm2 | ✅ Fully supported | Preferred for OSC 52 reliability |
| macOS | tmux on macOS | ✅ Works under | Same OSC 52 caveat as Linux |
| **WSL** | Windows Terminal | ✅ Fully supported | Pillow clipboard image grab disabled (Linux side) — use file picker |
| WSL | VS Code integrated | ✅ Works | VS Code's terminal handles OSC 52 |
| WSL | tmux inside WSL | ✅ Works under | Same OSC 52 caveat |
| **Windows native** | Windows Terminal | ✅ Supported | Pillow clipboard image grab available. Local HMMscan unavailable (`pyhmmer` is POSIX-only — use WSL2); BLASTN/BLASTP + everything else work natively |
| Windows native | ConPTY (cmd, PowerShell) | ⚠ Limited | Some Textual mouse modes unreliable; recommend Windows Terminal |
| Windows native | conhost.exe (legacy) | ❌ Unsupported | No true-color, no mouse — Textual will refuse to launch |

---

## Required terminal features

SpliceCraft hard-requires:

* **UTF-8 stdout encoding** — the plasmid map uses U+2800-U+28FF braille
  dots. Non-UTF-8 terminals will render the map as gibberish. The
  launcher refuses to start with a clear error message in this case.
  Fix: `export PYTHONIOENCODING=utf-8` or use a UTF-8 locale
  (`LANG=C.UTF-8` / `LANG=en_US.UTF-8`).
* **256-color ANSI** — feature colors degrade visibly without it.
  Every modern terminal supports this; only ancient `vt100` doesn't.
* **Mouse events** — required for click-to-select, drag-select on
  parts bin, hover tooltips on map.
* **Alternate screen buffer** — Textual uses it for the full-screen
  TUI; without it the terminal scrollback gets polluted.

SpliceCraft soft-requires (degrades on absence):

* **True-color (24-bit) ANSI** — used for ApEinfo feature colors; on
  256-color-only terminals colors approximate to nearest palette
  entry.
* **OSC 52 clipboard** — copy-to-clipboard falls back through a
  3-tier chain: Textual's API → OSC 52 escape → temp file in
  `<DATA_DIR>/clipboard/`. The last tier always works.

---

## Platform-specific implementations

### Data directory lock

* **POSIX** (Linux / macOS / WSL): `fcntl.flock(LOCK_EX | LOCK_NB)`
* **Windows native**: `msvcrt.locking(LK_NBLCK, 1)`

Either failure (module missing on a stripped-down interpreter) skips
locking with a warning — the cache coherence guarantee weakens but
SpliceCraft still launches.

Stale lock detection:
* POSIX: reads `/proc/<pid>/cmdline` to confirm PID hasn't been
  recycled to an unrelated process (long-uptime systems)
* macOS / Windows: relies on `os.kill(pid, 0)` for liveness check;
  PID-recycle to unrelated process can't be detected without psutil

### Clipboard

| Operation | Linux / WSL | macOS | Windows |
|---|---|---|---|
| Copy text (Ctrl+C selection) | Textual → OSC 52 → file fallback | Textual → OSC 52 → file fallback | Textual → OSC 52 (via `CONOUT$`) → file fallback |
| Paste image (Experiments tab) | ❌ Button disabled | ✅ via Pillow `ImageGrab` | ✅ via Pillow `ImageGrab` |

Linux/WSL has no pure-Python clipboard image API; users must drop
the file via the file picker. Pillow on Win/Mac handles the
clipboard bitmap directly.

### Signals

* `SIGUSR1` faulthandler stack-dump handler installed only when
  `signal.SIGUSR1` exists (POSIX). Windows lacks it — gracefully
  skipped.
* `SIGTERM` / `SIGHUP` translated to `KeyboardInterrupt` when
  available. Windows lacks `SIGHUP`; the `getattr(_signal, name,
  None)` check filters it out.

### File permissions

* `os.chmod(path, 0o600)` for logs / bundles is skipped on Windows
  (the call is a no-op there but we avoid the syscall for clarity).

### Subprocess invocation

`subprocess.run` is always called with `cmd` as a **list** (never
`shell=True`) so Windows quoting rules and POSIX globbing don't
diverge. The update install path (`splicecraft update`) follows
this convention end-to-end.

### Optional Python dependencies

| Dep | Purpose | Degradation if missing |
|---|---|---|
| `Pillow>=10.0` | Experiments image attach, clipboard grab | Image-attach button warns; clipboard grab disabled |
| `pyspellchecker>=0.8.0` | Experiments F7 spellcheck | F7 silently no-ops |
| `primer3-py` | Primer Tm calculation | Falls back to 2+4 GC rule (less accurate) |
| `rich-pixels>=3.0.0` | Image rendering in Experiments | Image renders as placeholder |
| `pyhmmer` | In-process HMMscan | HMMscan button warns |

`Pillow`, `pyspellchecker`, `primer3-py`, and `rich-pixels` are
unconditional hard deps in `pyproject.toml`, so a normal `pip install
splicecraft` brings them on every platform. **`pyhmmer` is a hard dep
only on POSIX** (`sys_platform != 'win32'`): it ships no Windows wheels
and HMMER's C core doesn't build on native Windows, so it is
intentionally omitted there. On native Windows, HMMscan is unavailable
(the button toasts + explains the WSL2 path) while BLASTN/BLASTP fall
back to the pure-Python engine; WSL2 reports as Linux and gets pyhmmer
normally. The runtime checks also cover other edge cases (e.g. `pip
install --no-deps` for an offline minimal install).

---

## Known issues

### Windows ConPTY / legacy `conhost.exe`

* Mouse drag-select on parts bin may not register all events under
  ConPTY. Use Windows Terminal instead.
* OSC 52 clipboard via `CONOUT$` works on Windows Terminal but not
  on `conhost.exe` (legacy console). Tier-3 file fallback kicks in.

### WSL → host clipboard

* Text copy via OSC 52 forwards through the WSL→Windows Terminal
  pipeline correctly.
* Image paste is disabled (Linux side has no pure-Python API);
  use the file picker after saving the screenshot from the Windows
  side.

### macOS Terminal.app vs iTerm2

* macOS Terminal.app's OSC 52 support is gated behind
  "Allow Mouse Reporting" + "Send escape sequences to clipboard"
  in Preferences. iTerm2 supports OSC 52 out of the box. If text
  copy doesn't reach the system clipboard, tier-3 file fallback
  still works — check `~/.local/share/splicecraft/clipboard/`.

### tmux / screen

OSC 52 passthrough requires:

```tmux
set -g allow-passthrough on
set -g set-clipboard on
```

Without these, OSC 52 escapes get filtered and tier-3 file fallback
takes over.

---

## Launch-time capability probe

On every launch, `_log_terminal_capabilities()` runs and emits a
`startup.terminal_capabilities` structured event with:

```json
{
  "encoding":       "utf-8",
  "is_tty":         true,
  "platform":       "linux",
  "blocking_count": 0,
  "warning_count":  0
}
```

Blocking failures (non-UTF-8 encoding) print to stderr and abort
launch with exit 1. Warnings (missing optional deps) log only and
allow launch to proceed.

Grep your log:

```bash
grep "startup.terminal_capabilities" ~/.local/share/splicecraft/logs/splicecraft.log
```

---

## Reporting compatibility issues

If SpliceCraft refuses to launch or renders badly on your terminal:

1. Run `splicecraft logs --bundle --out splicecraft-bundle.zip`
2. Open an issue at github.com/Binomica-Labs/SpliceCraft including
   the bundle + your terminal name / version + OS version.
3. The bundle includes the last 5 UI snapshots and sanitised
   settings; sequence content NEVER ships (privacy invariant #38).
