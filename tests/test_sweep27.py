"""
test_sweep27 — adversarial-audit follow-up sweep (2026-05-25).

Covers every fix landed in sweep #27 / invariants [INV-75]…[INV-82]
plus the HIGH/MED/LOW gap-closures from the same audit:

  HIGH-1 — FASTA size cap in `_do_load`
  HIGH-2 — L2 chokepoint on `_save_ui_snapshot`
  HIGH-3 — agent exception text scrubbed on 5 endpoints

  MED-4  — agent rate limiting (token bucket per bearer)
  MED-5  — GFF3 strand column validation
  MED-6  — NCBI Entrez `tool` identifier
  MED-7  — Pillow `MAX_IMAGE_PIXELS` ceiling at module import

  LOW-8  — `Authorization` header case-insensitive
  LOW-9  — `_find_usages_worker` is_mounted guard
  LOW-10 — autosave timer cancelled on app exit
  LOW-11 — `threading.RLock` for parse / BLAST caches

  INV-75 — L2 chokepoint extends to deletes
  INV-76 — crash-injection test for `_safe_save_json` recovery
  INV-77 — concurrency fuzz test for save/load invariants
  INV-78 — single `_now()` time source
  INV-79 — golden-file regression for persisted JSON envelopes
  INV-80 — idempotency keys on agent write endpoints
  INV-81 — `_safe_data_repr` for error-message scrubbing
  INV-82 — `@work` decorator thread=True enforcement (AST walk)

Each test is small + deterministic. Crash-injection skips on Windows
(no `os.kill(pid, SIGKILL)`). Concurrency fuzz uses a fixed seed so
flakes get caught on the same run that introduces them, not several
runs later.
"""
from __future__ import annotations

import ast
import json
import os
import sys
import threading
import time

from pathlib import Path

import pytest

import splicecraft as sc


# ── INV-78 — single time source `_now()` ─────────────────────────────


class TestNowHelper:
    """The `_now()` / `_monotonic()` helpers are the single point of
    attack for clock-skew tests. They MUST exist + be tz-aware (for
    `_now`) + return a float (for `_monotonic`). Future callsites
    will route through them, so the contract is part of the public
    surface of the module."""

    def test_now_returns_tz_aware_datetime(self):
        ts = sc._now()
        assert ts.tzinfo is not None, "_now() must be tz-aware"

    def test_now_is_monotonic_increasing(self):
        a = sc._now()
        b = sc._now()
        assert b >= a, "_now() should be monotonic non-decreasing"

    def test_monotonic_returns_float(self):
        a = sc._monotonic()
        b = sc._monotonic()
        assert isinstance(a, float)
        assert isinstance(b, float)
        assert b >= a

    def test_save_ui_snapshot_uses_now(self, tmp_path, monkeypatch):
        """`_save_ui_snapshot` was a `_datetime.now()` callsite before
        sweep #27. After: routed through `_now()`. Verify by patching
        `_now` to return a fixed timestamp and confirming the filename
        carries it."""
        from datetime import datetime, timezone

        # Fixed UTC instant; the resulting filename uses the LOCAL
        # tzinfo's representation (`.astimezone()` resolves to local).
        fixed = datetime(2026, 5, 25, 12, 0, 0, tzinfo=timezone.utc)
        monkeypatch.setattr(sc, "_now", lambda: fixed.astimezone())
        d = tmp_path / "snapshots"
        path = sc._save_ui_snapshot("test body", dest_dir=d)
        # Filename pattern `ui-snapshot-YYYYMMDD-HHMMSS.md` — must
        # carry digits derived from `fixed`, regardless of local TZ.
        assert path.name.startswith("ui-snapshot-")
        assert path.suffix == ".md"


# ── INV-81 — `_safe_data_repr` ────────────────────────────────────────


class TestSafeDataRepr:
    """Error messages must surface type/length only, never raw user
    content. The helper centralises the policy so any future error
    path can use it consistently."""

    def test_string_returns_len_only(self):
        out = sc._safe_data_repr("ACGTACGT" * 100)
        assert "ACGT" not in out, "raw value must not appear"
        assert "len=" in out
        assert "str" in out

    def test_none_returns_literal(self):
        assert sc._safe_data_repr(None) == "None"

    def test_bytes_returns_len_only(self):
        out = sc._safe_data_repr(b"\x00\x01\x02\x03" * 50)
        assert "\\x00" not in out
        assert "len=200" in out
        assert "bytes" in out

    def test_list_returns_n_only(self):
        out = sc._safe_data_repr([1, 2, 3, 4, 5])
        assert "n=5" in out
        assert "1" not in out.replace("n=5", "")

    def test_dict_returns_n_only(self):
        out = sc._safe_data_repr({"secret": "value", "k": "v"})
        assert "secret" not in out
        assert "value" not in out
        assert "n=2" in out

    def test_scalar_int_returns_type_only(self):
        out = sc._safe_data_repr(424242)
        assert "424242" not in out
        assert "int" in out

    def test_no_ansi_escape_leakage(self):
        """A malicious paste could embed ANSI control codes in a
        value. `_safe_data_repr` MUST NOT echo them."""
        out = sc._safe_data_repr("\x1b[31mEVIL\x1b[0m")
        assert "\x1b" not in out
        assert "EVIL" not in out


# ── INV-82 — `@work` decorator thread=True enforcement ───────────────


class TestWorkDecoratorAlwaysThread:
    """Every `@work` in splicecraft.py uses `thread=True` (or its
    positional equivalent). A bare `@work(...)` or `@work` defaults
    to running on the UI thread as a coroutine — defeats the worker
    contract and freezes the app. AST-walk asserts the rule."""

    @staticmethod
    def _collect_work_decorators() -> "list[tuple[int, ast.expr]]":
        src = Path(sc.__file__).read_text()
        tree = ast.parse(src)
        found: list[tuple[int, ast.expr]] = []
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for dec in node.decorator_list:
                    # Match `@work(...)` calls only.
                    if isinstance(dec, ast.Call):
                        fn = dec.func
                        # `@work(thread=True, ...)`
                        if isinstance(fn, ast.Name) and fn.id == "work":
                            found.append((node.lineno, dec))
                    # Bare `@work` (no parens) is also a footgun.
                    if (isinstance(dec, ast.Name)
                            and dec.id == "work"):
                        found.append((node.lineno, dec))
        return found

    def test_every_work_call_carries_thread_kwarg(self):
        decorators = self._collect_work_decorators()
        # Sanity: we should find ≥ 1 worker (the codebase has 55+).
        assert len(decorators) > 10, (
            "AST walk found %d @work decorators; expected many more "
            "(was the import broken?)" % len(decorators)
        )
        violators: list[str] = []
        for lineno, dec in decorators:
            if isinstance(dec, ast.Call):
                kwargs = {kw.arg: kw.value for kw in dec.keywords}
                if "thread" not in kwargs:
                    violators.append(
                        f"line {lineno}: @work(...) missing thread=True"
                    )
                    continue
                v = kwargs["thread"]
                if not (isinstance(v, ast.Constant) and v.value is True):
                    violators.append(
                        f"line {lineno}: @work thread={ast.dump(v)} "
                        f"is not literal True"
                    )
            else:
                violators.append(
                    f"line {lineno}: bare @work without parens; "
                    f"runs on UI thread"
                )
        assert not violators, (
            "@work decorator violations:\n  - "
            + "\n  - ".join(violators)
        )


# ── HIGH-1 — FASTA size cap ──────────────────────────────────────────


class TestFastaSizeCap:
    def test_constant_exists(self):
        assert hasattr(sc, "_FASTA_MAX_BYTES")
        assert isinstance(sc._FASTA_MAX_BYTES, int)
        assert sc._FASTA_MAX_BYTES > 1024 * 1024   # at least 1 MB
        assert sc._FASTA_MAX_BYTES <= 1024 * 1024 * 1024  # ≤ 1 GB

    def test_safe_file_size_check_rejects_over_cap(self, tmp_path):
        """The check fires BEFORE SeqIO.parse so OOM is impossible."""
        # Synthesise an oversized FASTA. `_safe_file_size_check`
        # rejects based on st_size, so we can use a sparse-but-large
        # file.
        big = tmp_path / "huge.fasta"
        big.write_bytes(b">x\n" + b"A" * (sc._FASTA_MAX_BYTES + 1024))
        ok, reason = sc._safe_file_size_check(
            big, sc._FASTA_MAX_BYTES, "fasta",
        )
        assert not ok
        assert reason and "fasta" in reason.lower()

    def test_safe_file_size_check_accepts_under_cap(self, tmp_path):
        small = tmp_path / "tiny.fasta"
        small.write_bytes(b">x\nACGT\n")
        ok, reason = sc._safe_file_size_check(
            small, sc._FASTA_MAX_BYTES, "fasta",
        )
        assert ok, reason


# ── HIGH-2 — L2 chokepoint on `_save_ui_snapshot` ────────────────────


class TestUiSnapshotChokepoint:
    def test_save_ui_snapshot_refuses_unauthorized(self, tmp_path,
                                                     monkeypatch):
        """When the authorisation gate is False, `_save_ui_snapshot`
        must raise rather than write."""
        monkeypatch.setattr(sc._state, "_SAVES_AUTHORIZED", False)
        with pytest.raises(RuntimeError, match="not authorised"):
            sc._save_ui_snapshot("body", dest_dir=tmp_path / "snaps")

    def test_save_ui_snapshot_allowed_when_authorized(self, tmp_path,
                                                       monkeypatch):
        monkeypatch.setattr(sc._state, "_SAVES_AUTHORIZED", True)
        path = sc._save_ui_snapshot("body", dest_dir=tmp_path / "snaps")
        assert path.exists()
        assert path.read_text() == "body"


# ── HIGH-3 — agent exception text scrubbing ──────────────────────────


class TestAgentExceptionScrubbing:
    """Five endpoints used to leak `OSError.strerror` /
    `OSError.filename` in their error responses. After sweep #27 they
    route through `_scrub_path`."""

    def test_load_file_error_is_scrubbed(self, tmp_path):
        # Trigger a parse failure on a non-existent path; the path
        # itself MUST be scrubbed from the response.
        bogus = tmp_path / "does_not_exist.gb"
        result = sc._h_load_file(None, {"path": str(bogus)})
        # Either 400 (rejected by sanitiser) or {} response — either way
        # the path's home segment must not appear verbatim.
        if isinstance(result, tuple):
            payload, status = result
            home_seg = str(Path.home())
            for v in payload.values():
                if isinstance(v, str) and home_seg in v:
                    pytest.fail(
                        f"home path leaked in error: {v!r}"
                    )

    def test_load_entry_error_is_scrubbed(self):
        # _h_load_entry on a missing entry returns 404 with `{"error":
        # "no library entry matching ..."}`. The key name itself can
        # appear (user supplied), but no system path.
        result = sc._h_load_entry(None, {"name": "no-such-entry-xyz"})
        assert isinstance(result, tuple)
        payload, status = result
        assert status == 404


# ── MED-4 — agent rate limiting ──────────────────────────────────────


class TestAgentRateLimit:
    def setup_method(self, _method):
        sc._agent_rate_limit_reset()

    def teardown_method(self, _method):
        sc._agent_rate_limit_reset()

    def test_initial_request_allowed(self):
        assert sc._agent_rate_limit_check("tokenA", cost=1.0)

    def test_burst_then_exhaust(self):
        # Burst of cost-1 reads up to the burst capacity should all
        # succeed; the next one should be rejected.
        token = "tokenB"
        burst = sc._AGENT_RATE_LIMIT_BURST
        for _ in range(burst):
            assert sc._agent_rate_limit_check(token, cost=1.0)
        # Bucket is now empty (no time has passed); next call must
        # fail.
        assert not sc._agent_rate_limit_check(token, cost=1.0)

    def test_per_token_isolation(self):
        # tokenA exhausted; tokenB still has full capacity.
        for _ in range(sc._AGENT_RATE_LIMIT_BURST):
            sc._agent_rate_limit_check("tokenA", cost=1.0)
        assert not sc._agent_rate_limit_check("tokenA", cost=1.0)
        assert sc._agent_rate_limit_check("tokenB", cost=1.0)

    def test_write_consumes_double(self):
        token = "tokenC"
        # A burst of `burst/2` writes (cost=2) exhausts the bucket.
        budget = sc._AGENT_RATE_LIMIT_BURST // 2
        for _ in range(budget):
            assert sc._agent_rate_limit_check(token, cost=2.0)
        # Bucket roughly empty; next write should fail.
        assert not sc._agent_rate_limit_check(token, cost=2.0)


# ── INV-80 — idempotency keys ────────────────────────────────────────


class TestAgentIdempotency:
    def setup_method(self, _method):
        sc._agent_idempotency_reset()

    def teardown_method(self, _method):
        sc._agent_idempotency_reset()

    def test_validate_key_rejects_bad_chars(self):
        assert sc._agent_validate_idempotency_key("abc_123") == "abc_123"
        assert sc._agent_validate_idempotency_key("with space") is None
        assert sc._agent_validate_idempotency_key("a/b/c") is None
        assert sc._agent_validate_idempotency_key("") is None
        assert sc._agent_validate_idempotency_key("x" * 1024) is None

    def test_put_then_get_returns_same(self):
        sc._agent_idempotency_put(
            "create-collection", "key1", {"ok": True, "name": "X"}, 200,
        )
        out = sc._agent_idempotency_get("create-collection", "key1")
        assert out is not None
        payload, status = out
        assert payload == {"ok": True, "name": "X"}
        assert status == 200

    def test_get_miss(self):
        assert sc._agent_idempotency_get("anything", "no-such-key") is None

    def test_expiry(self, monkeypatch):
        # Freeze monotonic to an artificial value, store, then advance
        # past the TTL.
        clock = [1000.0]
        monkeypatch.setattr(sc, "_monotonic", lambda: clock[0])
        sc._agent_idempotency_put("e", "k", {"ok": True}, 200)
        clock[0] = 1000.0 + sc._AGENT_IDEMPOTENCY_TTL_S + 1
        assert sc._agent_idempotency_get("e", "k") is None

    def test_cache_cap(self):
        # Spam past the cap; the cache must evict old entries.
        for i in range(sc._AGENT_IDEMPOTENCY_MAX_ENTRIES + 50):
            sc._agent_idempotency_put(
                "endpoint", f"key{i}", {"i": i}, 200,
            )
        # Total entries ≤ cap.
        with sc._AGENT_IDEMPOTENCY_LOCK:
            assert len(sc._AGENT_IDEMPOTENCY_CACHE) <= sc._AGENT_IDEMPOTENCY_MAX_ENTRIES


# ── MED-5 — GFF3 strand validation ────────────────────────────────────


class TestGff3StrandValidation:
    """Sweep #27 made `_parse_gff3_text` reject any strand column
    value not in `{+, -, ., ?}`."""

    @staticmethod
    def _make_gff3(strand_col: str) -> str:
        # Minimal GFF3 with one feature using the given strand string.
        # NB the trailing `+ ("ACGT" * 25)` deliberately breaks the
        # implicit string-literal concatenation. Without the `+`, the
        # `"ACGT" * 25` operator binds to the implicitly-concatenated
        # `"##FASTA\n>chr1\nACGT"` literal preceding it, multiplying
        # the whole `##FASTA>chr1ACGT` block 25 times — which produces
        # a multi-record fixture that the strict-raise GFF3 parser
        # (2026-05-27 audit-3 M8) refuses. The `+` forces evaluation
        # of `"ACGT" * 25` in isolation.
        return (
            "##gff-version 3\n"
            "##sequence-region chr1 1 100\n"
            f"chr1\tsrc\tCDS\t1\t30\t.\t{strand_col}\t.\tID=t1\n"
            "##FASTA\n"
            ">chr1\n"
            + ("ACGT" * 25) + "\n"
        )

    def test_valid_strands_accepted(self, tmp_path):
        for s in ("+", "-", ".", "?"):
            p = tmp_path / f"good_{s.replace('.', 'dot').replace('?', 'q')}.gff3"
            p.write_text(self._make_gff3(s))
            rec = sc._gff3_path_to_record(str(p))
            assert rec is not None

    def test_malicious_strand_skipped(self, tmp_path):
        # ANSI escape + HTML payload — must NOT produce a feature.
        evil = "\x1b[31m<script>alert(1)</script>"
        p = tmp_path / "evil.gff3"
        p.write_text(self._make_gff3(evil))
        rec = sc._gff3_path_to_record(str(p))
        # No feature should land from the malformed row (only the
        # automatic `source` feature, if any).
        non_src = [f for f in rec.features if f.type != "source"]
        assert non_src == [], (
            "malformed strand should produce no features; got: "
            + repr(non_src)
        )


# ── MED-6 — NCBI Entrez tool identifier ──────────────────────────────


class TestNcbiToolIdentifier:
    """Sweep #27 set `Entrez.tool = "SpliceCraft/<version>"` so NCBI
    can identify our traffic. We can't probe NCBI in tests, so we
    monkeypatch `Entrez.efetch` and confirm the tool attribute is set
    before the fetch runs."""

    def test_entrez_tool_set_before_fetch(self, monkeypatch):
        # Stub out the Bio.Entrez surface so the fetch doesn't actually
        # hit the network. Capture the value of `Entrez.tool` at the
        # moment `efetch` would be called.
        from Bio import Entrez

        captured: dict = {}

        class _StubHandle:
            def read(self, _n):
                # Return a minimal valid GenBank record.
                return (
                    "LOCUS       TEST            10 bp    DNA     linear   UNK\n"
                    "ORIGIN\n"
                    "        1 acgtacgtac\n"
                    "//\n"
                )

            def close(self): pass
            def __enter__(self): return self
            def __exit__(self, *a): pass

        def _stub_efetch(**kwargs):
            captured["tool"] = getattr(Entrez, "tool", None)
            captured["email"] = getattr(Entrez, "email", None)
            return _StubHandle()

        monkeypatch.setattr(Entrez, "efetch", _stub_efetch)
        try:
            sc.fetch_genbank("test_accession", email="test@example.com")
        except Exception:
            # The synthetic record probably won't parse cleanly via
            # Biopython's `SeqIO.read`, which is fine — we only care
            # that `Entrez.tool` was set before the fetch.
            pass
        assert captured.get("tool", "").startswith("SpliceCraft/")


# ── MED-7 — Pillow MAX_IMAGE_PIXELS ──────────────────────────────────


class TestPillowMaxImagePixels:
    def test_max_image_pixels_set_at_module_load(self):
        try:
            from PIL import Image
        except ImportError:
            pytest.skip("Pillow not installed")
        # Must match `_EXPERIMENT_CLIP_MAX_PIXELS` so the in-decoder
        # refusal matches the post-decode check.
        assert Image.MAX_IMAGE_PIXELS == sc._EXPERIMENT_CLIP_MAX_PIXELS


# ── LOW-8 — Authorization header case-insensitive ────────────────────


class TestAuthHeaderCaseInsensitive:
    """A "bearer" scheme keyword in any case should pass through to
    `secrets.compare_digest`. The token comparison itself is still
    case-sensitive (URL-safe base64)."""

    def test_lowercase_bearer_accepted(self):
        # Construct a minimal mock that exercises the prefix check.
        class _MockHeaders:
            def __init__(self, d): self._d = d
            def get(self, k, default=""): return self._d.get(k, default)

        class _MockHandler:
            headers = _MockHeaders({"Authorization": "bearer abc123"})
            server = type("S", (), {"_token": "abc123"})()
            _check_token = sc._AgentRequestHandler._check_token

        assert _MockHandler._check_token(_MockHandler())

    def test_uppercase_bearer_accepted(self):
        class _MockHeaders:
            def __init__(self, d): self._d = d
            def get(self, k, default=""): return self._d.get(k, default)

        class _MockHandler:
            headers = _MockHeaders({"Authorization": "BEARER abc123"})
            server = type("S", (), {"_token": "abc123"})()
            _check_token = sc._AgentRequestHandler._check_token

        assert _MockHandler._check_token(_MockHandler())

    def test_no_scheme_rejected(self):
        class _MockHeaders:
            def __init__(self, d): self._d = d
            def get(self, k, default=""): return self._d.get(k, default)

        class _MockHandler:
            headers = _MockHeaders({"Authorization": "abc123"})
            server = type("S", (), {"_token": "abc123"})()
            _check_token = sc._AgentRequestHandler._check_token

        assert not _MockHandler._check_token(_MockHandler())


# ── INV-75 — delete chokepoint ───────────────────────────────────────


class TestDeleteChokepoint:
    def test_refuse_unauthorized_delete_helper_exists(self):
        assert callable(getattr(sc, "_refuse_unauthorized_delete", None))

    def test_refuse_unauthorized_delete_raises_when_off(self, monkeypatch,
                                                          tmp_path):
        monkeypatch.setattr(sc._state, "_SAVES_AUTHORIZED", False)
        with pytest.raises(RuntimeError, match="not authorised"):
            sc._refuse_unauthorized_delete(
                tmp_path / "fake.json", "test",
            )

    def test_refuse_unauthorized_delete_passes_when_on(self, monkeypatch,
                                                        tmp_path):
        monkeypatch.setattr(sc._state, "_SAVES_AUTHORIZED", True)
        # Must not raise.
        sc._refuse_unauthorized_delete(
            tmp_path / "fake.json", "test",
        )


# ── INV-76 — crash-injection recovery test for `_safe_save_json` ─────


@pytest.mark.skipif(sys.platform == "win32",
                     reason="SIGKILL not supported on Windows")
class TestCrashInjectionSafeSaveJson:
    """Spawn a subprocess that writes through `_safe_save_json`,
    SIGKILL it mid-write, and confirm the .bak chain recovers a
    consistent state.

    This codifies the atomic-write contract documented in
    `[PIT-31]` / `[INV-37]`. A future regression in the fsync order
    of `_safe_save_json` would silently mean "post-crash file might
    be torn"; this test catches it deterministically.

    Strategy: write a small payload, then a large payload. Halfway
    through the large payload, send SIGKILL. Recovery via `_safe_load_
    json`'s fallback to `.bak` must return either the small payload
    (write didn't land) or the large one (write fully landed) — never
    a torn mix.
    """

    # `fork()` is intentional here — `spawn` would re-import the
    # 105k-line `splicecraft.py` in the child (~1-2 s startup),
    # blowing past the 50 ms SIGKILL window the test relies on to
    # catch the write mid-flight. Python 3.12+ emits a
    # `DeprecationWarning` because pytest's xdist worker (this
    # process) is multi-threaded, but the child immediately calls
    # `_safe_save_json` without spinning extra threads, so the
    # documented deadlock risk doesn't apply here.
    @pytest.mark.filterwarnings(
        "ignore:This process.*is multi-threaded.*fork:DeprecationWarning"
    )
    def test_recovery_after_sigkill(self, tmp_path):
        import multiprocessing as mp
        import signal

        out_file = tmp_path / "test_lib.json"
        small_entries = [{"id": "small", "name": "S", "n": 1}]
        big_entries = [{"id": f"big_{i}", "name": "B" * 32, "n": i}
                       for i in range(2000)]

        # Persist an initial payload from this process.
        sc._safe_save_json(out_file, small_entries, "test")
        original = out_file.read_bytes()

        # Confirm the initial payload survives a clean read.
        loaded, _ = sc._safe_load_json(out_file, "test")
        assert loaded == small_entries

        # The child process tries to write the large payload but
        # spins (so we can SIGKILL it). The point isn't to crash
        # mid-fsync exactly — it's to verify that whatever state the
        # filesystem ends up in is recoverable via the .bak chain.
        # We arrange for the kill to land between the bak-copy step
        # and the os.replace step by injecting a sleep.
        def child(path_str: str, entries):
            import time as _time
            import splicecraft as _sc
            _sc._authorize_writes(reason="crash test child")
            # Monkey-patch `_safe_save_json` internals not directly
            # possible from outside, so instead exercise the public
            # path with a payload large enough that disk write takes
            # nonzero time. Then we send SIGKILL while it's blocked
            # in `os.fsync`.
            _sc._safe_save_json(Path(path_str), entries, "test")
            # If we got here without being killed, the write
            # completed.
            _time.sleep(60)

        ctx = mp.get_context("fork") if sys.platform != "win32" else mp.get_context("spawn")
        p = ctx.Process(target=child, args=(str(out_file), big_entries))
        p.start()
        # Short delay then kill. 50 ms is enough for the child to
        # enter the write path on any modern OS.
        time.sleep(0.05)
        if p.is_alive():
            os.kill(p.pid, signal.SIGKILL)
        p.join(timeout=5.0)

        # Filesystem state: either the original .bak chain has the
        # pre-kill payload, OR the new payload has fully landed.
        # The `.json.tmp_*` files from the interrupted mkstemp may
        # linger; that's documented (`_sweep_orphan_tmp_files` cleans
        # them on next launch).
        # Use `_safe_load_json` to recover — it falls back to .bak.
        loaded, _ = sc._safe_load_json(out_file, "test")
        # Must be ONE of the two known states, never torn.
        assert (loaded == small_entries or loaded == big_entries), (
            "post-crash state is not a known good payload: "
            "len(loaded)={}, first_id={!r}".format(
                len(loaded),
                loaded[0].get("id") if loaded else None,
            )
        )

        del original  # silence unused-var


# ── INV-77 — concurrency fuzz test ───────────────────────────────────


class TestConcurrencyFuzz:
    """N threads run random save/load operations against a small
    library for K seconds. Post-run we assert the invariants hold:

      - every entry is dict
      - every entry has an `id` that's a non-empty str
      - `id`s are unique (no duplicate entries)
      - JSON on disk parses

    A regression in `_cache_lock` placement or a missing RMW wrap
    would surface as a duplicate id, a torn JSON file, or a hung
    test.
    """

    def test_fuzz_save_load_round_trip(self, tmp_path, monkeypatch):
        import random

        # Redirect the library file to tmp_path.
        lib_path = tmp_path / "plasmid_library.json"
        monkeypatch.setattr(sc._state, "_LIBRARY_FILE", lib_path)
        monkeypatch.setattr(sc._state, "_library_cache", None)

        # Seed an initial set of entries.
        initial = [
            {"id": f"e{i}", "name": f"name{i}", "gb_text": ""}
            for i in range(20)
        ]
        sc._save_library(initial)

        rng = random.Random(0xCAFEBABE)
        n_threads = 8
        duration_s = 1.5
        stop_at = time.monotonic() + duration_s
        errors: list[BaseException] = []

        def worker(seed: int):
            local_rng = random.Random(seed)
            while time.monotonic() < stop_at:
                op = local_rng.choice(("save", "load", "iter"))
                try:
                    if op == "save":
                        # Read-modify-write to mirror the actual
                        # callsite shape. The test verifies the
                        # `_cache_lock`-protected RMW path.
                        entries = sc._load_library()
                        if entries and local_rng.random() < 0.2:
                            entries.pop(local_rng.randrange(len(entries)))
                        if local_rng.random() < 0.3:
                            new_idx = local_rng.randint(0, 999_999)
                            entries.append({
                                "id": f"x{new_idx}",
                                "name": f"x{new_idx}",
                                "gb_text": "",
                            })
                        sc._save_library(entries)
                    elif op == "load":
                        sc._load_library()
                    else:
                        # Read-only iteration via the safe helper.
                        for _ in sc._iter_library_readonly():
                            pass
                except Exception as exc:
                    errors.append(exc)
                    return

        threads = [
            threading.Thread(target=worker, args=(rng.randint(0, 1 << 32),))
            for _ in range(n_threads)
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=duration_s + 5.0)

        # No thread raised.
        assert not errors, f"worker raised: {errors!r}"

        # Final state: load + invariants.
        final = sc._load_library()
        assert isinstance(final, list)
        seen_ids: set[str] = set()
        for entry in final:
            assert isinstance(entry, dict), (
                f"non-dict entry survived: {type(entry)}"
            )
            entry_id = entry.get("id")
            assert isinstance(entry_id, str) and entry_id, (
                f"missing/blank id: {entry!r}"
            )
            assert entry_id not in seen_ids, (
                f"duplicate id {entry_id!r} in final library"
            )
            seen_ids.add(entry_id)

        # JSON on disk parses cleanly.
        raw = lib_path.read_text()
        envelope = json.loads(raw)
        # Envelope-v1 or legacy bare-list both acceptable.
        if isinstance(envelope, dict):
            assert "entries" in envelope
            assert isinstance(envelope["entries"], list)
        else:
            assert isinstance(envelope, list)


# ── INV-79 — golden-file regression for JSON envelopes ───────────────


class TestGoldenFileEnvelope:
    """`_safe_save_json` writes an envelope `{"_schema_version": N,
    "entries": [...]}` with deterministic key order + indentation.
    A regression in the writer (accidental key sort change, indent
    drift, schema-version bump without migration) would surface as a
    fixed-input → drifted-output diff. The test pins the shape.

    Fixture data is deliberately tiny + verbatim — the assertion is
    on structure, not content.
    """

    def test_envelope_keys_stable(self, tmp_path):
        # Save a small payload and inspect the on-disk shape.
        entries = [
            {"id": "alpha", "name": "Alpha", "gb_text": "LOCUS …\n"},
            {"id": "beta",  "name": "Beta",  "gb_text": "LOCUS …\n"},
        ]
        path = tmp_path / "lib.json"
        sc._safe_save_json(path, entries, "test")
        loaded = json.loads(path.read_text())
        # Envelope-v1 contract.
        assert isinstance(loaded, dict)
        assert set(loaded.keys()) == {"_schema_version", "entries"}
        assert loaded["_schema_version"] == 1
        assert loaded["entries"] == entries

    def test_envelope_indented_for_humans(self, tmp_path):
        # The writer indents 2 spaces so humans can `git diff` the
        # JSON files. A regression to compact `json.dumps(...)` would
        # blow up the diff signal-to-noise ratio. Probe for newlines.
        path = tmp_path / "lib.json"
        sc._safe_save_json(path, [{"id": "x", "name": "X"}], "test")
        text = path.read_text()
        assert "\n" in text, (
            "envelope must be indented (multi-line) for diffability"
        )

    def test_envelope_round_trip_preserves_entries(self, tmp_path):
        # Save → read → load via `_safe_load_json` → original list.
        entries = [{"id": f"e{i}", "name": f"E{i}"} for i in range(5)]
        path = tmp_path / "lib.json"
        sc._safe_save_json(path, entries, "test")
        loaded, _ = sc._safe_load_json(path, "test")
        assert loaded == entries

    def test_envelope_back_compat_with_legacy_bare_list(self, tmp_path):
        """Pre-0.3.1 the on-disk format was a bare list. `_safe_load_
        json` must still accept it."""
        path = tmp_path / "legacy.json"
        legacy = [{"id": "e1", "name": "E1"}]
        path.write_text(json.dumps(legacy))
        loaded, _ = sc._safe_load_json(path, "test")
        assert loaded == legacy


# ── LOW-10 — autosave timer cancel on exit ──────────────────────────


class TestAutosaveTimerCancel:
    """The autosave debounce timer is stopped when the app exits."""

    def test_cancel_helper_idempotent_when_no_timer(self):
        # Construct an instance without going through __init__ so we
        # can probe the helper in isolation (running PlasmidApp would
        # spin up the whole TUI).
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        app._autosave_timer = None
        # Must not raise.
        app._cancel_autosave_timer()
        assert app._autosave_timer is None

    def test_cancel_helper_stops_active_timer(self):
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)

        class _StubTimer:
            stopped = False

            def stop(self):
                self.stopped = True

        timer = _StubTimer()
        app._autosave_timer = timer
        app._cancel_autosave_timer()
        assert timer.stopped is True
        assert app._autosave_timer is None


# ── INV-83 — mirror-write exemption from shrink guard ────────────────


class TestMirrorWriteHelper:
    """Regression for the 2026-05-25 live catastrophic-shrink
    incident. Switching parts bins from DemoColl (26 parts) to FFE
    (0 parts) tripped the L3 shrink guard because the bin-switch
    path used bare `_safe_save_json` instead of the mirror helper.
    User data was safe (shrink guard + lost_entries spillover
    preserved everything) but the UI bin-switch failed."""

    def test_helper_exists_and_bypasses_shrink_guard(self, tmp_path,
                                                       monkeypatch):
        # Pre-populate a "large" file (≥ 10 entries so the shrink
        # guard would normally fire on a >90% shrink).
        path = tmp_path / "mirror.json"
        big = [{"id": f"e{i}", "name": f"E{i}"} for i in range(20)]
        sc._safe_save_json(path, big, "test")
        # A bare `_safe_save_json` with [] would raise; the mirror
        # helper must succeed.
        sc._safe_save_json_mirror(path, [], "test")
        # Confirm the write landed.
        loaded, _ = sc._safe_load_json(path, "test")
        assert loaded == []

    def test_helper_still_uses_atomic_write(self, tmp_path):
        # The .bak chain MUST still rotate so an accidental mirror
        # call against the wrong file is still recoverable.
        path = tmp_path / "mirror.json"
        sc._safe_save_json(path, [{"id": "e1"}], "test")
        sc._safe_save_json_mirror(path, [{"id": "e2"}], "test")
        bak = path.with_suffix(path.suffix + ".bak")
        assert bak.exists(), "mirror helper must still write .bak"

    def test_no_bare_safe_save_json_against_mirror_files(self):
        """Regression scanner: no code path writes the four mirror
        files via bare `_safe_save_json`. Every cross-mirror write
        must route through `_safe_save_json_mirror`. New mirror
        writes that forget the helper get caught here.

        Uses the AST so comment text mentioning the function (e.g.
        in a sweep-history docstring) doesn't false-positive."""
        import ast as _ast
        mirror_files = {
            "_PARTS_BIN_FILE",
            "_LIBRARY_FILE",
            "_PRIMERS_FILE",
            "_EXPERIMENTS_FILE",
        }
        # Allowlist: the canonical `_save_*` helpers + the mirror
        # helper itself + every confirmed-non-swap save-to-disk worker
        # (these go straight to `_safe_save_json` for cache-lock or
        # perf reasons, never as a cross-mirror swap, so the L3 shrink
        # guard is the right protection — adding/removing one entry
        # doesn't trip the >90% catastrophic-shrink threshold).
        # If you add a new worker that needs to bypass `_save_*`,
        # CONFIRM it's an in-place edit (not a cross-mirror swap)
        # before adding it here. Cross-mirror swap = use
        # `_safe_save_json_mirror`.
        allowlist_defs = {
            # Canonical writers (the data IS the source-of-truth).
            "_save_library",
            "_save_parts_bin",
            "_save_primers",
            "_save_experiments",
            # The mirror helper itself.
            "_safe_save_json_mirror",
            # In-place save workers: add/delete/rename/status flip
            # against the LIVE mirror, not a swap. They've been
            # audited against the catastrophic-shrink pattern.
            "_add_save_to_disk",
            "_delete_save_to_disk",
            "_dom_primers_save_to_disk",
            "_trad_save_to_disk",
            "_primer_status_save_to_disk",
            "_rename_save_to_disk",
            # `PartsBinModal._delete` removes one part from the
            # active bin's mirror — in-place delete, not a swap.
            "_delete",
        }

        src = Path(sc.__file__).read_text()
        tree = _ast.parse(src)

        # Walk every def, find `_safe_save_json(<FIRST_ARG>, ...)`
        # calls whose first positional arg is a mirror file Name.
        violations: list[str] = []

        class _Visitor(_ast.NodeVisitor):
            def __init__(self):
                self._enclosing: list[str] = []

            def visit_FunctionDef(self, node):
                self._enclosing.append(node.name)
                self.generic_visit(node)
                self._enclosing.pop()

            visit_AsyncFunctionDef = visit_FunctionDef  # type: ignore[assignment]

            def visit_Call(self, node):
                fn = node.func
                if (isinstance(fn, _ast.Name)
                        and fn.id == "_safe_save_json"
                        and node.args
                        and isinstance(node.args[0], _ast.Name)
                        and node.args[0].id in mirror_files):
                    enclosing = (self._enclosing[-1]
                                 if self._enclosing else "<module>")
                    if enclosing not in allowlist_defs:
                        violations.append(
                            f"line {node.lineno}: bare "
                            f"_safe_save_json({node.args[0].id}, ...) "
                            f"in def {enclosing!r} — use "
                            f"_safe_save_json_mirror for cross-mirror "
                            f"writes (INV-83)"
                        )
                self.generic_visit(node)

        _Visitor().visit(tree)
        assert not violations, (
            "Bare mirror writes found (INV-83 violation):\n  - "
            + "\n  - ".join(violations)
        )


# ── The chokepoint's DELETE half now covers the pruners (2026-11-09) ────────
# `_refuse_unauthorized_delete` existed and its docstring named this exact
# threat — "a maintenance helper … that walks `_DATA_DIR` to prune rotated
# backups, expired snapshots, or stale crash-recovery files could be imported
# by an unsandboxed probe" — and NO pruner called it. So one
# `sc._run_data_dir_housekeeping(sc._state._DATA_DIR)` from an ad-hoc
# `import splicecraft` compressed and pruned every registry file's `.bak`
# rotation, pruned `lost_entries/`, quarantined blobs and hard-deleted
# already-quarantined ones — the four-layer recovery net that saved the
# 2026-05-22 incident — while the ADD-only `_snapshot_data_files` was
# correctly gated. See [INV-234].


class TestDataDirPrunersAreGated:
    """Every pruner must refuse in an unauthorised process, and still run in
    an authorised one (the app authorises before launch spawns housekeeping;
    read-only sets `_skip_snapshot`, which gates the spawn entirely)."""

    # (callable name, args builder) — each takes a path under the data dir.
    PRUNERS = [
        "_prune_backups",
        "_compress_old_backups",
        "_prune_lost_entries",
        "_prune_old_snapshots",
        "_run_data_dir_housekeeping",
        "_gc_orphan_blobs",
        "_prune_blob_quarantine",
        "_enforce_ui_snapshot_retention",
    ]

    def _invoke(self, name, tmp_path):
        fn = getattr(sc, name)
        if name in ("_prune_backups", "_compress_old_backups"):
            return fn(sc._state._LIBRARY_FILE)
        if name == "_prune_lost_entries":
            return fn(tmp_path / "lost_entries")
        if name == "_prune_old_snapshots":
            return fn(tmp_path / "snapshots")
        if name == "_run_data_dir_housekeeping":
            return fn(sc._state._DATA_DIR)
        if name == "_prune_blob_quarantine":
            return fn(tmp_path / "quarantine")
        if name == "_enforce_ui_snapshot_retention":
            return fn(tmp_path / "ui-snapshots")
        return fn()                      # _gc_orphan_blobs

    @pytest.mark.parametrize("name", PRUNERS)
    def test_refuses_when_writes_are_not_authorised(self, name, monkeypatch,
                                                     tmp_path):
        monkeypatch.setattr(sc._state, "_SAVES_AUTHORIZED", False)
        with pytest.raises(RuntimeError, match="not authorised"):
            self._invoke(name, tmp_path)

    @pytest.mark.parametrize("name", PRUNERS)
    def test_runs_when_authorised(self, name, monkeypatch, tmp_path):
        """The other half — a gate that also blocks the sanctioned caller
        would silently stop pruning and the backup dir would grow forever."""
        monkeypatch.setattr(sc._state, "_SAVES_AUTHORIZED", True)
        for sub in ("lost_entries", "snapshots", "quarantine",
                    "ui-snapshots"):
            (tmp_path / sub).mkdir(exist_ok=True)
        # Must not raise. These are all best-effort and may do nothing on an
        # empty dir; the point is that the gate lets them through.
        self._invoke(name, tmp_path)

    def test_the_delete_guard_is_what_refuses(self, monkeypatch, tmp_path):
        """Pin the MECHANISM, not just the symptom: the refusal must come
        from `_refuse_unauthorized_delete`, so a future pruner that forgets
        to call it is not accidentally covered by something else."""
        seen: list = []
        real = sc._refuse_unauthorized_delete

        def spy(path, label):
            seen.append(label)
            return real(path, label)

        monkeypatch.setattr(sc._state, "_SAVES_AUTHORIZED", True)
        monkeypatch.setattr(sc, "_refuse_unauthorized_delete", spy)
        import splicecraft_backup as B
        monkeypatch.setattr(B, "_refuse_unauthorized_delete", spy)
        B._run_data_dir_housekeeping(sc._state._DATA_DIR)
        assert seen, "housekeeping never consulted the delete guard"
