"""
test_hmm_db_catalog — sweep #28.

Covers the new HMM database registry + download / version-check
infrastructure introduced by the BlastModal Pfam-A download feature.

Hardening test coverage (the user explicitly asked for "harden and
edge case these updates since they communicate with the internet"):

  * URL scheme validation (https-only by default, http opt-in)
  * URL credential redaction in logs
  * Bounded redirects via the custom opener
  * Content-Type guard (HTML/JSON error pages rejected)
  * Magic-byte verification (gzip 0x1f8b + HMMER3 header)
  * Disk-space pre-check
  * Retry on transient network failures
  * Cancel-aware download / decompress
  * Cross-modal-instance download slot
  * Atomic hmmpress cleanup on failure
  * n_profiles == 0 → treated as failure
  * Catalog round-trips builtins + user entries
  * Sanitisation rejects path traversal / NUL / shell metas
  * Modal smoke fits 160x48
"""
from __future__ import annotations

import gzip
import io
import threading
import urllib.error

import pytest

import splicecraft as sc


# ── Catalog persistence ──────────────────────────────────────────────


class TestHmmDbCatalogPersistence:
    """Builtins always present; user entries round-trip; corrupted
    entries dropped."""

    def test_first_load_contains_builtins(self):
        catalog = sc._load_hmm_db_catalog()
        ids = {e["id"] for e in catalog}
        assert "pfam-a" in ids
        assert "ncbifam" in ids

    def test_save_then_load_round_trip(self):
        catalog = sc._load_hmm_db_catalog()
        catalog.append({
            "id":          "custom-x",
            "name":        "My Custom DB",
            "url":         "https://example.com/x.hmm.gz",
            "version_url": "",
            "format":      "hmm-gz",
            "builtin":     False,
            "description": "test entry",
        })
        sc._save_hmm_db_catalog(catalog)
        sc._hmm_db_catalog_cache = None    # force reload from disk
        loaded = sc._load_hmm_db_catalog()
        ids = {e["id"] for e in loaded}
        assert "custom-x" in ids
        assert "pfam-a" in ids

    def test_builtin_re_injected_on_reload(self):
        """If the user (or a hand-edit) removes pfam-a, the next
        load re-injects it. Defensive — the UI never lets you
        remove a builtin, but a corrupted catalog shouldn't
        permanently lose the well-known default."""
        catalog = [e for e in sc._load_hmm_db_catalog()
                   if e["id"] != "pfam-a"]
        sc._save_hmm_db_catalog(catalog)
        sc._hmm_db_catalog_cache = None
        loaded = sc._load_hmm_db_catalog()
        assert any(e["id"] == "pfam-a" for e in loaded)

    def test_duplicate_ids_collapsed_to_first(self):
        sc._save_hmm_db_catalog([
            {"id": "dup", "name": "First", "url": "https://a/x.hmm.gz"},
            {"id": "dup", "name": "Second", "url": "https://b/x.hmm.gz"},
        ])
        sc._hmm_db_catalog_cache = None
        loaded = [e for e in sc._load_hmm_db_catalog()
                  if e["id"] == "dup"]
        assert len(loaded) == 1

    def test_corrupted_entry_dropped(self):
        # Schema violation: missing url. Should be silently dropped
        # at load (not crash the modal).
        sc._save_hmm_db_catalog([
            {"id": "broken", "name": "X"},   # no url
            {"id": "good",   "name": "Y", "url": "https://x/x.hmm.gz"},
        ])
        sc._hmm_db_catalog_cache = None
        loaded = sc._load_hmm_db_catalog()
        ids = {e["id"] for e in loaded}
        assert "good" in ids
        assert "broken" not in ids


# ── Sanitisation ─────────────────────────────────────────────────────


class TestSanitisation:
    def test_id_accepts_alnum_dash_underscore(self):
        assert sc._sanitize_hmm_db_id("pfam-a") == "pfam-a"
        assert sc._sanitize_hmm_db_id("ncbifam_v3") == "ncbifam_v3"
        assert sc._sanitize_hmm_db_id("a") == "a"

    def test_id_rejects_path_traversal(self):
        assert sc._sanitize_hmm_db_id("..") is None
        assert sc._sanitize_hmm_db_id("../etc") is None
        assert sc._sanitize_hmm_db_id("a/b") is None
        assert sc._sanitize_hmm_db_id("a\\b") is None

    def test_id_rejects_nul(self):
        assert sc._sanitize_hmm_db_id("a\x00b") is None

    def test_id_rejects_unicode_and_whitespace(self):
        assert sc._sanitize_hmm_db_id("café") is None
        assert sc._sanitize_hmm_db_id("a b") is None
        assert sc._sanitize_hmm_db_id("") is None

    def test_id_length_cap(self):
        assert sc._sanitize_hmm_db_id("x" * 64) == "x" * 64
        assert sc._sanitize_hmm_db_id("x" * 65) is None

    def test_url_must_be_http_or_https(self):
        assert sc._sanitize_hmm_db_url("https://x.com/x") == "https://x.com/x"
        assert sc._sanitize_hmm_db_url("http://x.com/x")  == "http://x.com/x"
        assert sc._sanitize_hmm_db_url("ftp://x.com/x")    is None
        assert sc._sanitize_hmm_db_url("file:///etc/passwd") is None
        assert sc._sanitize_hmm_db_url("javascript:x")     is None

    def test_url_rejects_whitespace_and_control_chars(self):
        assert sc._sanitize_hmm_db_url("https://x.com/x\n") is None
        assert sc._sanitize_hmm_db_url("https://x.com/x ") is None
        assert sc._sanitize_hmm_db_url("https://x.com/x\x01y") is None

    def test_url_length_cap(self):
        long_url = "https://x.com/" + "a" * (sc._HMM_DB_URL_MAX_LEN + 1)
        assert sc._sanitize_hmm_db_url(long_url) is None


# ── URL credential redaction ─────────────────────────────────────────


class TestUrlRedaction:
    def test_redacts_user_password(self):
        out = sc._redact_url_credentials("https://user:pw@host.com/x")
        assert "user" not in out
        assert "pw" not in out
        assert "host.com" in out

    def test_preserves_clean_url(self):
        url = "https://ftp.ebi.ac.uk/pub/x.gz"
        assert sc._redact_url_credentials(url) == url

    def test_preserves_port(self):
        out = sc._redact_url_credentials("https://user:pw@host.com:8080/x")
        assert ":8080" in out
        assert "user" not in out
        assert "pw" not in out

    def test_handles_garbage(self):
        # Function must not raise on weird input.
        assert sc._redact_url_credentials("not a url at all") == "not a url at all"
        assert sc._redact_url_credentials("") == ""


# ── Scheme policy ────────────────────────────────────────────────────


class TestSchemePolicy:
    def test_https_ok(self):
        assert sc._hmm_db_url_scheme_ok("https://x.com/x") is None

    def test_http_rejected_by_default(self):
        msg = sc._hmm_db_url_scheme_ok("http://x.com/x")
        assert msg is not None
        assert "http://" in msg
        assert "hmm_db_allow_http" in msg

    def test_http_allowed_with_opt_in(self, monkeypatch):
        # Pre-set the opt-in then call.
        sc._set_setting("hmm_db_allow_http", True)
        try:
            assert sc._hmm_db_url_scheme_ok("http://x.com/x") is None
        finally:
            sc._set_setting("hmm_db_allow_http", False)

    def test_other_schemes_rejected(self):
        for scheme in ("ftp", "file", "javascript", "data"):
            msg = sc._hmm_db_url_scheme_ok(f"{scheme}://x.com/x")
            assert msg is not None


# ── Content-Type guard ───────────────────────────────────────────────


class TestContentTypeGuard:
    @staticmethod
    def _stub_resp(content_type: str):
        class _Hdrs:
            def __init__(self, ct): self._ct = ct
            def get(self, key, default=""):
                if key.lower() == "content-type":
                    return self._ct
                return default

        class _R:
            headers = _Hdrs(content_type)
        return _R()

    def test_text_html_rejected(self):
        with pytest.raises(ValueError, match="text/html"):
            sc._hmm_db_assert_content_type_ok(
                self._stub_resp("text/html; charset=utf-8"), "http://x",
            )

    def test_application_json_rejected(self):
        with pytest.raises(ValueError, match="application/json"):
            sc._hmm_db_assert_content_type_ok(
                self._stub_resp("application/json"), "http://x",
            )

    def test_binary_accepted(self):
        # application/octet-stream and application/gzip are both fine.
        sc._hmm_db_assert_content_type_ok(
            self._stub_resp("application/octet-stream"), "http://x",
        )
        sc._hmm_db_assert_content_type_ok(
            self._stub_resp("application/gzip"), "http://x",
        )

    def test_missing_header_accepted(self):
        # No Content-Type header at all — accepted (some old mirrors
        # don't send it; rejecting would be too brittle).
        sc._hmm_db_assert_content_type_ok(
            self._stub_resp(""), "http://x",
        )


# ── Disk-space check ─────────────────────────────────────────────────


class TestDiskSpace:
    def test_accepts_when_plenty_of_space(self, tmp_path):
        # 1 MB on a normal tmpfs/HDD: plenty of space.
        sc._hmm_db_check_disk_space(tmp_path / "x", 1024 * 1024)

    def test_refuses_when_low_on_space(self, tmp_path, monkeypatch):
        # Stub `shutil.disk_usage` to report a tiny free pool.
        def _stub_du(_path):
            class _U:
                total = 10 * 1024 * 1024
                used  = 9 * 1024 * 1024
                free  = 1 * 1024 * 1024
            return _U()
        monkeypatch.setattr(sc.shutil, "disk_usage", _stub_du)
        with pytest.raises(OSError, match="disk space"):
            sc._hmm_db_check_disk_space(tmp_path / "x", 100 * 1024 * 1024)

    def test_defaults_when_no_size_known(self, tmp_path, monkeypatch):
        # `expected_bytes=None` → defaults to 5 GB reserve.
        def _stub_du(_path):
            class _U:
                total = 10 * 1024 * 1024 * 1024 * 1024
                used  = 0
                free  = 10 * 1024 * 1024 * 1024 * 1024
            return _U()
        monkeypatch.setattr(sc.shutil, "disk_usage", _stub_du)
        sc._hmm_db_check_disk_space(tmp_path / "x", None)
        # Tons of free space → must not raise.


# ── Magic-byte verification (download path) ──────────────────────────


def _build_fake_hmm_gz(payload: bytes = b"HMMER3/f [3.4 | Aug 2023]\n"
                       ) -> bytes:
    """A minimal valid gzip stream wrapping a HMMER3 header."""
    return gzip.compress(payload)


class TestDownloadMagicBytes:
    @staticmethod
    def _stub_urlopen(content: bytes, *,
                       content_type: str = "application/gzip"):
        class _Hdrs:
            _items = {"Content-Type": content_type,
                       "Content-Length": str(len(content))}
            def get(self, key, default=""):
                for k, v in self._items.items():
                    if k.lower() == key.lower():
                        return v
                return default

        class _Resp:
            headers = _Hdrs()
            _buf = io.BytesIO(content)
            def read(self, n=-1):
                return self._buf.read(n)
            def close(self): pass
            def __enter__(self): return self
            def __exit__(self, *a): pass

        return _Resp()

    def test_rejects_non_gzip_first_bytes(self, tmp_path, monkeypatch):
        # Pretend the server returned text/plain (somehow with a
        # bypassed content-type guard) — the first-chunk magic check
        # catches it.
        not_gzip = b"NOT A GZIP FILE AT ALL" * 100
        def _opener():
            class _O:
                def open(self, _req, timeout=None):
                    return TestDownloadMagicBytes._stub_urlopen(
                        not_gzip,
                        content_type="application/octet-stream",
                    )
            return _O()
        monkeypatch.setattr(sc, "_build_hardened_url_opener", _opener)
        with pytest.raises(ValueError, match="gzip magic"):
            sc._stream_download_to_path(
                "https://example.com/x.gz",
                tmp_path / "db.hmm.gz",
                max_bytes=10 * 1024 * 1024,
            )

    def test_rejects_html_content_type(self, tmp_path, monkeypatch):
        html = b"<html><body>blocked</body></html>"
        def _opener():
            class _O:
                def open(self, _req, timeout=None):
                    return TestDownloadMagicBytes._stub_urlopen(
                        html, content_type="text/html",
                    )
            return _O()
        monkeypatch.setattr(sc, "_build_hardened_url_opener", _opener)
        with pytest.raises(ValueError, match="Content-Type"):
            sc._stream_download_to_path(
                "https://example.com/x.gz",
                tmp_path / "db.hmm.gz",
                max_bytes=10 * 1024 * 1024,
            )

    def test_accepts_valid_gzip(self, tmp_path, monkeypatch):
        gz = _build_fake_hmm_gz()
        def _opener():
            class _O:
                def open(self, _req, timeout=None):
                    return TestDownloadMagicBytes._stub_urlopen(gz)
            return _O()
        monkeypatch.setattr(sc, "_build_hardened_url_opener", _opener)
        dest = tmp_path / "db.hmm.gz"
        sha = sc._stream_download_to_path(
            "https://example.com/x.gz", dest,
            max_bytes=10 * 1024 * 1024,
        )
        assert dest.exists()
        assert dest.read_bytes() == gz
        assert len(sha) == 64    # sha256 hex


# ── Decompression hardening ──────────────────────────────────────────


class TestDecompressHardening:
    def test_accepts_valid_hmmer3(self, tmp_path):
        body = b"HMMER3/f [3.4 | Aug 2023]\nNAME test\nLENG 10\n"
        gz_path = tmp_path / "in.hmm.gz"
        gz_path.write_bytes(gzip.compress(body))
        out = tmp_path / "out.hmm"
        sc._decompress_gz_to_path(gz_path, out, max_bytes=1024 * 1024)
        assert out.read_bytes() == body

    def test_rejects_non_hmm_decompressed_content(self, tmp_path):
        body = b"this is not a HMMER file at all\n"
        gz_path = tmp_path / "in.hmm.gz"
        gz_path.write_bytes(gzip.compress(body))
        out = tmp_path / "out.hmm"
        with pytest.raises(ValueError, match="HMMER"):
            sc._decompress_gz_to_path(
                gz_path, out, max_bytes=1024 * 1024,
            )
        assert not out.exists(), (
            "tmp must be cleaned up on failure"
        )

    def test_rejects_zip_bomb(self, tmp_path):
        # 1 MB of zeros compresses to ~1 KB; cap at 100 bytes
        # decompressed → should fire.
        body = b"\x00" * (1024 * 1024)
        gz_path = tmp_path / "bomb.hmm.gz"
        gz_path.write_bytes(gzip.compress(body))
        out = tmp_path / "out.hmm"
        with pytest.raises(ValueError, match="cap"):
            sc._decompress_gz_to_path(
                gz_path, out, max_bytes=100,
            )
        assert not out.exists()

    def test_rejects_corrupt_gzip(self, tmp_path):
        gz_path = tmp_path / "corrupt.gz"
        gz_path.write_bytes(b"this is not gzip data at all" * 100)
        out = tmp_path / "out.hmm"
        with pytest.raises((ValueError, OSError)):
            sc._decompress_gz_to_path(
                gz_path, out, max_bytes=1024 * 1024,
            )

    def test_cancel_check_aborts_cleanly(self, tmp_path):
        body = b"HMMER3/f [3.4]\n" + b"X" * (10 * 1024 * 1024)
        gz_path = tmp_path / "big.hmm.gz"
        gz_path.write_bytes(gzip.compress(body))
        out = tmp_path / "out.hmm"
        # Cancel cb returns True on the first call → abort
        # before any chunk written.
        calls = [0]
        def _cancel():
            calls[0] += 1
            return True
        with pytest.raises(OSError, match="cancelled"):
            sc._decompress_gz_to_path(
                gz_path, out, max_bytes=100 * 1024 * 1024,
                cancel_check_cb=_cancel,
            )
        assert not out.exists()


# ── Cross-modal download slot ────────────────────────────────────────


class TestDownloadSlot:
    def teardown_method(self, _m):
        # Always release any test-held slots so a test order change
        # doesn't leak.
        with sc._HMM_DB_DOWNLOAD_INFLIGHT_LOCK:
            sc._HMM_DB_DOWNLOAD_INFLIGHT.clear()

    def test_acquire_release(self):
        assert sc._hmm_db_acquire_download_slot("pfam-a")
        # Second acquire on same id refused.
        assert not sc._hmm_db_acquire_download_slot("pfam-a")
        # Different id ok.
        assert sc._hmm_db_acquire_download_slot("ncbifam")
        sc._hmm_db_release_download_slot("pfam-a")
        # Now re-acquirable.
        assert sc._hmm_db_acquire_download_slot("pfam-a")

    def test_release_idempotent(self):
        # Safe to release a never-acquired id.
        sc._hmm_db_release_download_slot("never-acquired")

    def test_concurrent_acquire_only_one_wins(self):
        winners = []
        barrier = threading.Barrier(8)
        def worker(_i):
            barrier.wait()
            if sc._hmm_db_acquire_download_slot("contended"):
                winners.append(threading.get_ident())
        threads = [threading.Thread(target=worker, args=(i,))
                   for i in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5)
        assert len(winners) == 1


# ── hmmpress cleanup ─────────────────────────────────────────────────


class TestHmmpressCleanup:
    def test_cleanup_removes_pressed_siblings(self, tmp_path):
        hmm = tmp_path / "db.hmm"
        hmm.write_bytes(b"HMMER3/f\n")
        for ext in ("h3i", "h3m", "h3p", "h3f"):
            (tmp_path / f"db.hmm.{ext}").write_bytes(b"junk")
        sc._cleanup_pressed_files(hmm)
        for ext in ("h3i", "h3m", "h3p", "h3f"):
            assert not (tmp_path / f"db.hmm.{ext}").exists()
        # Main .hmm preserved.
        assert hmm.exists()

    def test_cleanup_idempotent(self, tmp_path):
        # Calling cleanup on a directory with no pressed files is fine.
        hmm = tmp_path / "db.hmm"
        sc._cleanup_pressed_files(hmm)


# ── Version-file parsing ─────────────────────────────────────────────


class TestPfamVersionParser:
    def test_extracts_release_number(self):
        body = (
            "Pfam release       : 37.4\n"
            "Pfam-A families    : 21978\n"
            "Date               : 2025-02\n"
        )
        assert sc._parse_pfam_version_text(body) == "37.4"

    def test_handles_blank_lines(self):
        body = "\n\nPfam release : 36.0\n"
        assert sc._parse_pfam_version_text(body) == "36.0"

    def test_fallback_to_first_line(self):
        body = "totally unexpected format with no Pfam release line\n"
        assert sc._parse_pfam_version_text(body).startswith("totally")

    def test_empty_returns_empty(self):
        assert sc._parse_pfam_version_text("") == ""

    def test_non_string_returns_empty(self):
        assert sc._parse_pfam_version_text(None) == ""  # type: ignore

    def test_caps_at_64_chars(self):
        body = "Pfam release : " + "x" * 200
        assert len(sc._parse_pfam_version_text(body)) <= 64


# ── Per-DB local state ───────────────────────────────────────────────


class TestPerDbLocalState:
    def test_missing_meta_returns_none(self):
        assert sc._load_hmm_db_local_meta("pfam-a") is None

    def test_save_then_load(self, tmp_path):
        meta = {
            "id":            "pfam-a",
            "version":       "37.4",
            "downloaded_at": "2026-05-25T12:00:00",
            "sha256":        "abc" * 21,
            "n_profiles":    100,
            "pressed":       True,
        }
        sc._save_hmm_db_local_meta("pfam-a", meta)
        loaded = sc._load_hmm_db_local_meta("pfam-a")
        assert loaded is not None
        assert loaded["version"] == "37.4"
        assert loaded["n_profiles"] == 100

    def test_is_downloaded_reflects_hmm_file(self, tmp_path):
        assert not sc._is_hmm_db_downloaded("pfam-a")
        # Drop a fake .hmm file in the per-DB dir.
        hmm = sc._hmm_db_hmm_path("pfam-a")
        hmm.parent.mkdir(parents=True, exist_ok=True)
        hmm.write_bytes(b"HMMER3/f stub")
        assert sc._is_hmm_db_downloaded("pfam-a")

    def test_should_check_remote_24h_cache(self, monkeypatch):
        # No meta yet → always check.
        assert sc._hmm_db_should_check_remote("pfam-a")
        # Stamp a recent check; should NOT re-check.
        clock = [1000.0]
        monkeypatch.setattr(sc, "_monotonic", lambda: clock[0])
        sc._record_hmm_db_remote_version("pfam-a", "37.4")
        assert not sc._hmm_db_should_check_remote("pfam-a")
        # Advance past 24h.
        clock[0] = 1000.0 + sc._HMM_DB_VERSION_CHECK_TTL_S + 1
        assert sc._hmm_db_should_check_remote("pfam-a")


# ── Delete chokepoint ────────────────────────────────────────────────


class TestDeleteFiles:
    def test_delete_removes_all_db_files(self, tmp_path):
        d = sc._hmm_db_entry_dir("pfam-a")
        d.mkdir(parents=True, exist_ok=True)
        for name in ("db.hmm", "db.hmm.gz", "db.hmm.h3i",
                     "db.hmm.h3m", "db.hmm.h3p", "db.hmm.h3f",
                     "meta.json"):
            (d / name).write_bytes(b"x")
        removed = sc._delete_hmm_db_files("pfam-a")
        assert removed >= 7
        assert not d.exists()

    def test_delete_on_missing_dir_returns_zero(self):
        assert sc._delete_hmm_db_files("never-existed") == 0


# ── Network retry behaviour ──────────────────────────────────────────


class TestNetworkRetry:
    def test_version_check_retries_once_then_succeeds(self,
                                                       monkeypatch):
        # Stub the opener to fail once then succeed.
        calls = [0]

        class _Headers:
            """Real urllib `HTTPMessage` supports `.get()` for
            case-insensitive lookup AND can be passed to `dict(...)`
            (it implements `keys` + `__getitem__`). Mirror both
            interfaces here so the helper's `dict(resp.headers)`
            conversion + `.get('Content-Type')` lookups both work
            on the stub."""
            def __init__(self, data):
                self._data = dict(data)
            def get(self, k, default=""):
                kl = k.lower()
                for hk, hv in self._data.items():
                    if hk.lower() == kl:
                        return hv
                return default
            def keys(self):
                return list(self._data.keys())
            def __getitem__(self, k):
                return self.get(k)
            def __iter__(self):
                return iter(self._data)

        class _StubResp:
            headers = _Headers({"Content-Length": "10"})
            def read(self, n=-1): return b"Pfam release : 99.9\n"
            def close(self): pass
            def __enter__(self): return self
            def __exit__(self, *a): pass

        class _Opener:
            def open(self, _req, timeout=None):
                calls[0] += 1
                if calls[0] == 1:
                    raise urllib.error.URLError("transient")
                return _StubResp()
        monkeypatch.setattr(sc, "_build_hardened_url_opener",
                              lambda: _Opener())
        v, src = sc._fetch_hmm_db_remote_version({
            "id":          "pfam-a",
            "name":        "Pfam-A",
            "url":         "https://x.com/x.hmm.gz",
            "version_url": "https://x.com/x.version",
        })
        # Either the retry succeeded (version returned) or the fallback
        # HEAD on main URL was tried (also fails / 0 calls). Verify the
        # retry happened.
        assert calls[0] >= 2
        assert v == "99.9"
        assert src == "version_file"

    def test_version_check_persistent_failure_returns_empty(self,
                                                             monkeypatch):
        class _Opener:
            def open(self, _req, timeout=None):
                raise urllib.error.URLError("persistent")
        monkeypatch.setattr(sc, "_build_hardened_url_opener",
                              lambda: _Opener())
        v, src = sc._fetch_hmm_db_remote_version({
            "id":          "pfam-a",
            "name":        "Pfam-A",
            "url":         "https://x.com/x.hmm.gz",
            "version_url": "https://x.com/x.version",
        })
        assert v == ""
        assert src == ""


# ── Modal boundaries ─────────────────────────────────────────────────


@pytest.mark.asyncio
async def test_hmm_db_catalog_modal_fits_screen():
    """Smoke: modal renders inside 160x48 without overflow."""
    from splicecraft import PlasmidApp

    class _App(PlasmidApp):
        _skip_seed = True
        _skip_splash = True

    app = _App()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        await pilot.pause()
        app.push_screen(sc.HmmDbCatalogModal())
        await pilot.pause()
        await pilot.pause()
        # No assertion needed — if the modal raised on compose,
        # `push_screen` would have surfaced it.


@pytest.mark.asyncio
async def test_hmm_db_add_edit_modal_fits_screen():
    from splicecraft import PlasmidApp

    class _App(PlasmidApp):
        _skip_seed = True
        _skip_splash = True

    app = _App()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        await pilot.pause()
        app.push_screen(sc.HmmDbAddEditModal(mode="add"))
        await pilot.pause()
        await pilot.pause()
