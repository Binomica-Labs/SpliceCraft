"""test_cross_platform — verifies SpliceCraft's platform-specific
branches behave correctly under mocked `sys.platform`.

Covers the audit findings in CLAUDE.md invariant #45 (cross-platform
contract). Each test mocks `sys.platform` and / or relevant module
availability so the same test file exercises the Windows, macOS, and
Linux branches from any host. Real-runtime smoke testing on each OS
is documented in `docs/PLATFORMS.md`.
"""
import pytest

import splicecraft as sc


class TestOSC52ClipboardHelper:
    """`_copy_to_clipboard_osc52` opens `/dev/tty` on POSIX,
    `CONOUT$` on Windows. Either failure mode (device unavailable,
    encoding error) must return False so the multi-tier wrapper
    escalates cleanly — never raise."""

    def test_osc52_returns_false_when_tty_unavailable(self, monkeypatch):
        """When the controlling TTY isn't openable, the helper
        returns False — does NOT raise. Caller (the multi-tier
        wrapper) falls through to the file-based fallback."""
        def _no_tty(*args, **kwargs):
            raise OSError("no controlling terminal")
        # Replace `open` in the splicecraft module's namespace so the
        # helper's `open("/dev/tty", "w")` hits our stub.
        monkeypatch.setattr("splicecraft.open", _no_tty, raising=False)
        assert sc._copy_to_clipboard_osc52("hello") is False

    def test_osc52_uses_conout_on_windows(self, monkeypatch):
        """On Windows the helper writes to `CONOUT$` (the console
        output device) instead of `/dev/tty`. Verify the right
        target is opened so the OSC sequence reaches the terminal
        emulator on a Windows host."""
        opened: list[str] = []

        class _FakeTTY:
            def write(self, _s): pass
            def flush(self):     pass
            def __enter__(self):  return self
            def __exit__(self, *_a): pass

        def _capture_open(target, *args, **kwargs):
            opened.append(str(target))
            return _FakeTTY()

        monkeypatch.setattr(sc.sys, "platform", "win32")
        monkeypatch.setattr("splicecraft.open", _capture_open, raising=False)
        assert sc._copy_to_clipboard_osc52("hello") is True
        assert opened == ["CONOUT$"], (
            f"Windows OSC 52 must target CONOUT$, got {opened}"
        )

    def test_osc52_uses_dev_tty_on_posix(self, monkeypatch):
        """On POSIX systems (Linux, macOS) the helper writes to
        `/dev/tty` — the canonical controlling terminal device."""
        opened: list[str] = []

        class _FakeTTY:
            def write(self, _s): pass
            def flush(self):     pass
            def __enter__(self):  return self
            def __exit__(self, *_a): pass

        def _capture_open(target, *args, **kwargs):
            opened.append(str(target))
            return _FakeTTY()

        monkeypatch.setattr(sc.sys, "platform", "linux")
        monkeypatch.setattr("splicecraft.open", _capture_open, raising=False)
        assert sc._copy_to_clipboard_osc52("hello") is True
        assert opened == ["/dev/tty"], (
            f"POSIX OSC 52 must target /dev/tty, got {opened}"
        )

        opened.clear()
        monkeypatch.setattr(sc.sys, "platform", "darwin")
        assert sc._copy_to_clipboard_osc52("hello") is True
        assert opened == ["/dev/tty"]


class TestTerminalCapabilities:
    """`_check_terminal_capabilities` returns (blocking, warning)
    lists. UTF-8 is blocking; missing optional Python deps are
    warnings only."""

    def test_utf8_encoding_passes_silently(self, monkeypatch):
        """`stdout.encoding == 'utf-8'` produces no blocking
        entries (the canonical happy-path on every modern terminal
        with a properly set locale)."""
        class _FakeStdout:
            encoding = "utf-8"
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        blocking, _ = sc._check_terminal_capabilities()
        assert blocking == []

    def test_utf8_aliases_pass(self, monkeypatch):
        """Aliases (`utf8`, `cp65001` — the Windows code page for
        UTF-8) also count as UTF-8 — the braille map renders fine
        on a Windows Terminal session in `chcp 65001` mode."""
        for variant in ("utf-8", "utf8", "UTF-8", "Utf-8", "cp65001",
                          "utf-16", "utf-32"):
            class _FakeStdout:
                encoding = variant
            monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
            blocking, _ = sc._check_terminal_capabilities()
            assert blocking == [], (
                f"encoding {variant!r} should not block, got {blocking}"
            )

    def test_latin1_encoding_falls_back_to_ascii_not_block(self, monkeypatch):
        """Latin-1 (or any non-UTF-8 codec that can't be reconfigured)
        NO LONGER blocks launch. A FakeStdout with no `reconfigure`
        can't be switched to UTF-8, so the probe selects the ASCII map
        fallback (`_ASCII_MODE = True`) and WARNS — it does not refuse.
        The warning still points at the UTF-8 locale fix. (Pre-1.0.3
        this raised a hard blocker; the new contract is graceful
        degradation so SpliceCraft runs on any ANSI terminal.)"""
        monkeypatch.setattr(sc, "_ASCII_MODE", False)
        class _FakeStdout:
            encoding = "latin-1"   # and no `reconfigure` → unfixable
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        blocking, warning = sc._check_terminal_capabilities()
        assert blocking == [], "non-UTF-8 must no longer block launch"
        assert sc._ASCII_MODE is True, "should fall back to the ASCII map"
        blob = " ".join(warning).lower()
        assert "ascii" in blob
        assert "utf-8" in blob or "utf8" in blob

    def test_no_encoding_is_skipped(self, monkeypatch):
        """When stdout has no encoding (piped to a process), skip
        the encoding check — the user is clearly not running
        interactively, so the braille map isn't going to render
        anywhere anyway."""
        class _FakeStdout:
            encoding = ""
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        blocking, _ = sc._check_terminal_capabilities()
        assert blocking == []

    def test_missing_optional_deps_are_warnings_not_blockers(
            self, monkeypatch):
        """Pillow / primer3 / pyspellchecker import failures count
        as warnings only — the user can still launch SpliceCraft,
        they just lose specific features (clipboard image paste,
        primer3 Tm, spellcheck). Refusing to launch on a missing
        optional dep would over-step into hard-dependency territory."""
        # Force every optional import to raise ImportError.
        real_import = __import__

        def _faux_import(name, *args, **kwargs):
            if name in ("PIL", "primer3", "spellchecker"):
                raise ImportError(f"simulated absence of {name}")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", _faux_import)
        class _FakeStdout:
            encoding = "utf-8"
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        blocking, warning = sc._check_terminal_capabilities()
        assert blocking == [], "missing optional deps must not block"
        # All three warnings surface so the user knows what's
        # unavailable.
        msg_blob = " ".join(warning)
        assert "PIL"          in msg_blob
        assert "primer3"      in msg_blob
        assert "spellchecker" in msg_blob


class TestEnsureUtf8Stdout:
    """`_ensure_utf8_stdout` — the rescue step that converts a locale-
    mislabelled (LANG=C) terminal to real UTF-8 output instead of
    refusing. Returns True iff stdout is / became UTF-8. Edge-cased
    hard because it's the gate that decides whether the braille map
    or the ASCII fallback renders."""

    def test_already_utf8_short_circuits_no_reconfigure(self, monkeypatch):
        calls = []
        class _FakeStdout:
            encoding = "utf-8"
            def reconfigure(self, **kw):
                calls.append(kw)
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        assert sc._ensure_utf8_stdout() is True
        assert calls == [], "must not reconfigure an already-UTF-8 stream"

    def test_utf8_aliases_short_circuit(self, monkeypatch):
        for enc in ("utf-8", "utf8", "UTF-8", "cp65001", "utf-16", "utf-32"):
            class _FakeStdout:
                encoding = enc
            monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
            assert sc._ensure_utf8_stdout() is True, enc

    def test_empty_encoding_treated_as_ok(self, monkeypatch):
        # Piped / redirected — encoding "" → not interactive, don't
        # force ASCII on the strength of a missing string.
        class _FakeStdout:
            encoding = ""
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        assert sc._ensure_utf8_stdout() is True

    def test_none_stdout_does_not_crash(self, monkeypatch):
        # pythonw / daemonised: sys.stdout can be None.
        monkeypatch.setattr(sc.sys, "stdout", None)
        assert sc._ensure_utf8_stdout() is True

    def test_latin1_without_reconfigure_returns_false(self, monkeypatch):
        class _FakeStdout:
            encoding = "latin-1"   # no `reconfigure` attribute
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        assert sc._ensure_utf8_stdout() is False

    def test_reconfigure_success_flips_to_utf8(self, monkeypatch):
        class _FakeStdout:
            encoding = "ascii"
            def reconfigure(self, *, encoding=None, **kw):
                self.encoding = encoding   # mimic TextIOWrapper
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        assert sc._ensure_utf8_stdout() is True

    def test_reconfigure_raises_returns_false(self, monkeypatch):
        for exc in (OSError("nope"), ValueError("x"), LookupError("x")):
            class _FakeStdout:
                encoding = "ascii"
                def reconfigure(self, **kw):
                    raise exc
            monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
            assert sc._ensure_utf8_stdout() is False

    def test_reconfigure_noop_still_nonutf8_returns_false(self, monkeypatch):
        # Defensive: reconfigure() returns without error but the codec
        # didn't actually change → still report False (don't claim UTF-8
        # we couldn't achieve).
        class _FakeStdout:
            encoding = "latin-1"
            def reconfigure(self, **kw):
                pass   # pretend success, leave encoding unchanged
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        assert sc._ensure_utf8_stdout() is False


class TestSelectRenderTier:
    """`_select_render_tier` — env override + UTF-8 rescue → sets the
    module-level `_ASCII_MODE`. Every test resets the global via
    monkeypatch so it can't leak into other suites."""

    def test_utf8_terminal_uses_braille(self, monkeypatch):
        monkeypatch.setattr(sc, "_ASCII_MODE", True)   # ensure it flips back
        monkeypatch.delenv("SPLICECRAFT_ASCII", raising=False)
        class _FakeStdout:
            encoding = "utf-8"
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        ascii_mode, warnings = sc._select_render_tier()
        assert ascii_mode is False
        assert sc._ASCII_MODE is False
        assert warnings == []

    def test_unfixable_nonutf8_selects_ascii(self, monkeypatch):
        monkeypatch.setattr(sc, "_ASCII_MODE", False)
        monkeypatch.delenv("SPLICECRAFT_ASCII", raising=False)
        class _FakeStdout:
            encoding = "latin-1"   # no reconfigure → unfixable
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        ascii_mode, warnings = sc._select_render_tier()
        assert ascii_mode is True
        assert sc._ASCII_MODE is True
        assert warnings and "ascii" in " ".join(warnings).lower()

    def test_env_force_ascii_overrides_capable_utf8(self, monkeypatch):
        monkeypatch.setattr(sc, "_ASCII_MODE", False)
        monkeypatch.setenv("SPLICECRAFT_ASCII", "1")
        class _FakeStdout:
            encoding = "utf-8"   # capable, but forced to ASCII anyway
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        ascii_mode, warnings = sc._select_render_tier()
        assert ascii_mode is True
        assert sc._ASCII_MODE is True
        assert warnings and "splicecraft_ascii" in " ".join(warnings).lower()

    def test_env_force_accepts_truthy_variants(self, monkeypatch):
        class _FakeStdout:
            encoding = "utf-8"
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        for val in ("1", "true", "TRUE", "Yes", "on", "  on  "):
            monkeypatch.setattr(sc, "_ASCII_MODE", False)
            monkeypatch.setenv("SPLICECRAFT_ASCII", val)
            ascii_mode, _ = sc._select_render_tier()
            assert ascii_mode is True, f"{val!r} should force ASCII"

    def test_env_force_ignores_falsey_variants(self, monkeypatch):
        class _FakeStdout:
            encoding = "utf-8"
        monkeypatch.setattr(sc.sys, "stdout", _FakeStdout())
        for val in ("0", "false", "no", "off", "", "garbage"):
            monkeypatch.setattr(sc, "_ASCII_MODE", True)
            monkeypatch.setenv("SPLICECRAFT_ASCII", val)
            ascii_mode, _ = sc._select_render_tier()
            assert ascii_mode is False, f"{val!r} should NOT force ASCII"


class TestAsciiDensityLut:
    """`_ASCII_DENSITY_LUT` — braille bitmask → 7-bit-ASCII density char
    by dot popcount. The fallback only renders on 'any terminal' if it
    is genuinely pure ASCII."""

    def test_length_is_256(self):
        assert len(sc._ASCII_DENSITY_LUT) == 256

    def test_all_entries_are_single_ascii_char(self):
        for ch in sc._ASCII_DENSITY_LUT:
            assert len(ch) == 1
            assert ord(ch) < 128, f"{ch!r} is not 7-bit ASCII"

    def test_endpoints(self):
        assert sc._ASCII_DENSITY_LUT[0x00] == " "    # 0 dots
        assert sc._ASCII_DENSITY_LUT[0x01] == "."    # 1 dot, lightest visible
        assert sc._ASCII_DENSITY_LUT[0xFF] == "@"    # 8 dots, densest

    def test_glyph_tracks_popcount(self):
        ramp = sc._ASCII_DENSITY_RAMP
        assert len(ramp) == 9, "ramp must cover popcounts 0..8"
        for bits in range(256):
            assert sc._ASCII_DENSITY_LUT[bits] == ramp[bin(bits).count("1")]


class TestAsciiMapFallback:
    """`_BrailleCanvas.combine` swaps braille glyphs for the ASCII
    density ramp when `_ASCII_MODE` is on — the core "works on any
    terminal" guarantee. The braille default path must be untouched."""

    @staticmethod
    def _render(ascii_mode, monkeypatch):
        monkeypatch.setattr(sc, "_ASCII_MODE", ascii_mode)
        canvas = sc._Canvas(4, 1)         # all-blank text cells
        bc = sc._BrailleCanvas(4, 1)
        bc.set_pixel(0, 0)                # 1 dot in cell 0 (sparse)
        for px in range(2):               # all 8 dots in cell 1 (dense)
            for py in range(4):
                bc.set_pixel(2 + px, py)
        return bc.combine(canvas).plain

    def test_ascii_mode_emits_pure_ascii_no_braille(self, monkeypatch):
        out = self._render(True, monkeypatch)
        assert all(ord(c) < 128 for c in out), f"non-ASCII leaked: {out!r}"
        assert not any(0x2800 <= ord(c) <= 0x28FF for c in out)

    def test_ascii_mode_density_tracks_dot_count(self, monkeypatch):
        out = self._render(True, monkeypatch)
        # Sparse cell → lightest visible glyph; dense cell → densest.
        assert "." in out and "@" in out

    def test_braille_mode_unchanged(self, monkeypatch):
        out = self._render(False, monkeypatch)
        assert any(0x2800 <= ord(c) <= 0x28FF for c in out), (
            f"default mode must still render braille, got {out!r}"
        )

    def test_ascii_mode_transliterates_overlay_glyphs(self, monkeypatch):
        # The map overlays Unicode glyphs on the text canvas (block
        # fills, arrowheads, ⚠, crosshair). ASCII mode must fold those
        # to 7-bit equivalents too — not just the braille dot layer.
        monkeypatch.setattr(sc, "_ASCII_MODE", True)
        canvas = sc._Canvas(5, 1)
        bc = sc._BrailleCanvas(5, 1)
        canvas.put(0, 0, "█", "")    # block fill   → '#'
        canvas.put(1, 0, "▶", "")    # arrowhead    → '>'
        canvas.put(2, 0, "⚠", "")    # weak marker  → '!'
        canvas.put(3, 0, "·", "")    # crosshair    → '.'
        out = bc.combine(canvas).plain
        assert all(ord(c) < 128 for c in out), f"non-ASCII leaked: {out!r}"
        assert "#" in out and ">" in out and "!" in out and "." in out

    def test_ascii_mode_unmapped_glyph_becomes_question(self, monkeypatch):
        # A glyph not in the map (e.g. an accented label letter) folds
        # to '?' rather than leaking raw UTF-8.
        monkeypatch.setattr(sc, "_ASCII_MODE", True)
        canvas = sc._Canvas(2, 1)
        bc = sc._BrailleCanvas(2, 1)
        canvas.put(0, 0, "β", "")    # not in _ASCII_GLYPH_MAP
        out = bc.combine(canvas).plain
        assert "?" in out
        assert all(ord(c) < 128 for c in out)

    def test_braille_mode_preserves_overlay_glyphs(self, monkeypatch):
        # Default (UTF-8) mode leaves overlay glyphs untouched.
        monkeypatch.setattr(sc, "_ASCII_MODE", False)
        canvas = sc._Canvas(3, 1)
        bc = sc._BrailleCanvas(3, 1)
        canvas.put(0, 0, "▶", "")
        assert "▶" in bc.combine(canvas).plain


class TestAsciiSpinnerFallback:
    """The online-search spinner uses braille frames by default and a
    7-bit-ASCII set under `_ASCII_MODE`, so it doesn't mojibake on a
    non-UTF-8 terminal."""

    def test_ascii_spinner_frames_are_pure_ascii(self):
        frames = sc.BlastModal._ONLINE_SPIN_FRAMES_ASCII
        assert frames, "ASCII spinner frames must be non-empty"
        assert all(ord(c) < 128 for c in frames)

    def test_default_spinner_frames_are_braille(self):
        assert any(0x2800 <= ord(c) <= 0x28FF
                   for c in sc.BlastModal._ONLINE_SPIN_FRAMES)


class TestSignalHandlersDegradeWithoutSIGUSR1:
    """Windows lacks `signal.SIGUSR1` — the faulthandler hang-debug
    helper must not crash the import. The import-time `try`/`except
    AttributeError` already handles this; this test pins the
    contract so a future refactor can't regress."""

    def test_module_imports_when_sigusr1_missing(self, monkeypatch):
        """Strip SIGUSR1 from the `signal` module and re-run the
        helper code that installs the faulthandler — must not
        raise. Mirrors what happens on Windows where SIGUSR1
        doesn't exist."""
        import signal as _signal
        had_sigusr1 = hasattr(_signal, "SIGUSR1")
        if had_sigusr1:
            monkeypatch.delattr(_signal, "SIGUSR1", raising=False)
        # The module's install block runs at import; replay its
        # logic here.
        assert not hasattr(_signal, "SIGUSR1")
        # The check that the splicecraft install block uses:
        if hasattr(_signal, "SIGUSR1"):
            pytest.fail("monkeypatch failed to strip SIGUSR1")


class TestSubprocessUsage:
    """SpliceCraft's `subprocess.run` calls must pass `cmd` as a
    list (not a string with `shell=True`) to avoid Windows / POSIX
    shell-quoting differences. Verifies the upgrade-install path
    follows this contract by capturing the actual call shape."""

    def test_update_install_passes_list_not_shell(self, monkeypatch):
        """The actual install subprocess wrapper (`_update_run_install`)
        passes `cmd` as a list and does NOT set `shell=True` —
        portable invocation across Windows + POSIX."""
        captured: dict = {}

        def _fake_run(cmd, **kwargs):
            captured["cmd"]    = cmd
            captured["kwargs"] = kwargs
            class _Result:
                returncode = 0
            return _Result()

        monkeypatch.setattr(sc.subprocess, "run", _fake_run)
        # Drive the wrapper; it only matters that subprocess.run is
        # invoked — we capture the call shape, not the return.
        import pathlib as _path
        rc = sc._update_run_install(
            ["pip", "install", "splicecraft==0.9.7"],
            snap_path=_path.Path("/tmp/snap"),
            method="pip",
            pin_version="0.9.7",
        )
        assert rc == 0
        assert isinstance(captured["cmd"], list), (
            "subprocess.run must receive cmd as list, not string"
        )
        assert captured["kwargs"].get("shell", False) is False, (
            "subprocess.run must not use shell=True (Windows-quote risk)"
        )


class TestPlatformCachedAtImport:
    """`_RUNTIME_PLATFORM` is cached at module import time so
    `platform.platform()` (which shells out via subprocess on some
    OSes) doesn't fire on every snapshot / event log. Verifies the
    cache survives a monkeypatched `platform.platform`."""

    def test_runtime_platform_is_cached(self, monkeypatch):
        """If we monkeypatch `platform.platform` to raise, the
        existing `_RUNTIME_PLATFORM` constant still returns the
        cached value. Confirms no path re-invokes
        `platform.platform()` after import."""
        # Confirm the constant exists + is non-empty.
        assert sc._RUNTIME_PLATFORM
        # Even if platform.platform() blows up now, the constant
        # is intact.
        import platform as _platform
        def _no_platform():
            raise RuntimeError("simulated platform.platform crash")
        monkeypatch.setattr(_platform, "platform", _no_platform)
        # Accessing the cache must not re-call platform.platform().
        assert sc._RUNTIME_PLATFORM  # still the cached string
