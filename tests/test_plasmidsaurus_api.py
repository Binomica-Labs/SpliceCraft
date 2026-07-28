"""Plasmidsaurus REST API client + agent endpoints + GUI fetch modal
(sweep #29).

The network is fully mocked at ``splicecraft_search._build_hardened_url_opener``
— no real egress ever happens. Tests are auto-sandboxed by the conftest
``_protect_user_data`` fixture, so the library-import paths write to a
throwaway data dir (the suite is also authorised for the L2 chokepoint).
"""
from __future__ import annotations

import io
import json
import zipfile

import pytest

import splicecraft as sc
import splicecraft_search as _search
import splicecraft_fileio as _fileio


# ── response / opener stubs ──────────────────────────────────────────────────
def _resp(body: bytes, content_type: str = "application/json"):
    class _H:
        _d = {"Content-Type": content_type, "Content-Length": str(len(body))}

        def get(self, k, default=""):
            for kk, v in self._d.items():
                if kk.lower() == k.lower():
                    return v
            return default

    class _R:
        headers = _H()

        def __init__(self):
            self._b = io.BytesIO(body)

        def read(self, n=-1):
            return self._b.read(n)

        def close(self):
            pass

    return _R()


def _mkgb(rid: str = "samp") -> str:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    r = SeqRecord(
        Seq("ATGAAACGCATTAGCACCACCATTACCACCACCATCGGTACCTAA" * 3),
        id=rid, name=rid)
    r.annotations["molecule_type"] = "DNA"
    r.annotations["topology"] = "circular"
    return sc._record_to_gb_text(r)


def _zip_bytes(names=("my_sample_1.gbk", "my_sample_2.gbk")) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for i, n in enumerate(names):
            zf.writestr(n, _mkgb(f"s{i}"))
    return buf.getvalue()


class _Router:
    """Fake hardened opener; routes by URL to canned responses."""

    def __init__(self, *, token="tok", items=None,
                 link="https://files.example/x.zip", zipbytes=None):
        self.token = token
        self.items = items if items is not None else []
        self.link = link
        self.zipbytes = zipbytes if zipbytes is not None else _zip_bytes()
        self.calls: list[str] = []

    def open(self, req, timeout=None):
        url = req.get_full_url()
        self.calls.append(url)
        if url.endswith("/oauth/token"):
            assert (req.get_header("Authorization") or "").startswith("Basic ")
            return _resp(json.dumps({"access_token": self.token}).encode())
        if "/api/items" in url:
            return _resp(json.dumps(self.items).encode())
        if url.endswith("/results"):
            return _resp(json.dumps({"link": self.link}).encode())
        if url == self.link:
            return _resp(self.zipbytes, content_type="application/zip")
        raise AssertionError("unexpected URL " + url)


@pytest.fixture
def use_router(monkeypatch):
    def _install(router: _Router) -> _Router:
        monkeypatch.setattr(_search, "_build_hardened_url_opener",
                            lambda: router)
        return router
    return _install


# A realistic listing: one downloadable run, one shipping label (complete but
# undated, no results), one canceled order.
_ORDERS = [
    {"code": "AAAAAA", "status": "complete", "order_name": "run one",
     "product_name": "plasmid_low_copy", "quantity": 2,
     "done_date": "2026-06-01T10:00:00+00:00", "gross": 30.0},
    {"code": "BBBBBB", "status": "complete", "order_name": "label",
     "product_name": "ups_shipping_label", "quantity": 1,
     "done_date": None, "gross": 12.0},
    {"code": "CCCCCC", "status": "canceled", "order_name": "dropped",
     "product_name": "plasmid", "quantity": 1,
     "done_date": "2026-05-01T10:00:00+00:00", "gross": 0.0},
]


@pytest.fixture(autouse=True)
def _clear_ps_order_cache():
    """The GUI order cache is module-level; a leak across tests would let one
    test serve another's listing."""
    blank = {"at": 0.0, "items": [], "token": "", "key": ""}
    sc._PS_ORDERS_CACHE.update(blank)
    yield
    sc._PS_ORDERS_CACHE.update(blank)


# ── item-code sanitiser ──────────────────────────────────────────────────────
class TestItemCodeSanitizer:
    def test_valid_uppercased(self):
        assert sc._sanitize_plasmidsaurus_item_code("abc123") == "ABC123"
        assert sc._sanitize_plasmidsaurus_item_code(" ABC123 ") == "ABC123"

    def test_rejects_bad_shapes(self):
        for bad in ("ABC12", "ABCDEFG", "AB C12", "ABC123/../x",
                    "../secr", "http://x", "", "ABC-12", None, 123):
            assert sc._sanitize_plasmidsaurus_item_code(bad) is None, bad


# ── credentials (env-first, then settings) ───────────────────────────────────
class TestCredentials:
    def test_env_first(self, monkeypatch):
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "envid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "envsec")
        sc._set_setting("plasmidsaurus_client_id", "setid")
        sc._set_setting("plasmidsaurus_client_secret", "setsec")
        assert sc._plasmidsaurus_credentials() == ("envid", "envsec")

    def test_settings_fallback(self, monkeypatch):
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_ID", raising=False)
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_SECRET", raising=False)
        sc._set_setting("plasmidsaurus_client_id", "setid")
        sc._set_setting("plasmidsaurus_client_secret", "setsec")
        assert sc._plasmidsaurus_credentials() == ("setid", "setsec")

    def test_partial_is_none(self, monkeypatch):
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_ID", raising=False)
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_SECRET", raising=False)
        sc._set_setting("plasmidsaurus_client_id", "only-id")
        sc._set_setting("plasmidsaurus_client_secret", "")
        cid, sec = sc._plasmidsaurus_credentials()
        assert cid == "only-id" and sec is None



class TestCredentialShapeHint:
    """A rejection should say WHAT is wrong when the shape gives it away.

    2026-07-28 (user report): the Settings fields held a Plasmidsaurus
    *website login* — email in Client ID, account password in Client Secret —
    and the modal reported only a bare `HTTP 401`. An email address sitting in
    a field that wants 32 hex characters is diagnosable at a glance, so the
    message now says so.
    """

    ID, SEC = "a" * 32, "b" * 64

    def test_a_correct_pair_produces_no_hint(self):
        # Append-safe: "" must be the answer whenever nothing looks wrong,
        # or every successful-shape rejection grows a bogus tail.
        assert sc._plasmidsaurus_credential_hint(self.ID, self.SEC) == ""
        assert sc._plasmidsaurus_credential_hint(self.ID.upper(), self.SEC) == ""

    def test_an_email_in_the_client_id_is_named(self):
        h = sc._plasmidsaurus_credential_hint("seb@example.org", "hunter2")
        assert "email address" in h and "website login" in h
        # The wrong-KIND diagnosis wins outright — listing lengths underneath
        # would bury it.
        assert "characters, expected" not in h

    def test_wrong_length_or_non_hex_is_named_per_field(self):
        h = sc._plasmidsaurus_credential_hint("abc", self.SEC)
        assert "Client ID" in h and "Client Secret" not in h
        h = sc._plasmidsaurus_credential_hint(self.ID, "not-hex-" + "z" * 56)
        assert "Client Secret" in h and "isn't hex" in h
        h = sc._plasmidsaurus_credential_hint("a" * 30, "b" * 60)
        assert "30 characters, expected 32" in h
        assert "60 characters, expected 64" in h

    def test_a_missing_half_is_not_reported_as_malformed(self):
        # `_plasmidsaurus_credentials` already refuses a missing half with its
        # own message; complaining about "0 characters" here would double up.
        assert sc._plasmidsaurus_credential_hint("", self.SEC) == ""
        assert sc._plasmidsaurus_credential_hint(self.ID, "") == ""

    def test_the_hint_never_echoes_a_credential(self):
        """It reaches the UI and the log — lengths and character classes only."""
        secret = "sk-DO-NOT-ECHO-ME-abc123"
        h = sc._plasmidsaurus_credential_hint("seb@example.org", secret)
        assert secret not in h
        h = sc._plasmidsaurus_credential_hint("zz" * 20, secret)
        assert secret not in h and "zz" * 20 not in h

    @staticmethod
    def _raise_status(monkeypatch, code: int):
        import urllib.error

        class _O:
            def open(self, req, timeout=None):
                raise urllib.error.HTTPError(
                    req.get_full_url(), code, "err", {}, None)
        monkeypatch.setattr(_search, "_build_hardened_url_opener", lambda: _O())

    def test_a_401_carries_the_hint(self, monkeypatch):
        """The whole point: the message the modal shows must contain it."""
        self._raise_status(monkeypatch, 401)
        with pytest.raises(OSError) as ei:
            sc._plasmidsaurus_oauth_token("seb@example.org", "hunter2")
        msg = str(ei.value)
        assert "HTTP 401" in msg and "email address" in msg

    def test_hostile_input_can_never_raise(self):
        """It runs INSIDE the `except HTTPError` block that reports the
        rejection — anything it raises replaces a clear "HTTP 401" with a
        traceback, i.e. the diagnostic destroys the diagnosis. Probed with
        non-str types, which every one of these used to raise on."""
        for cid, sec in ((None, None), (12345678, "b" * 64),
                         ("a" * 32, b"b" * 64), (["a"], "b" * 64),
                         ("a" * 32, {"k": 1}), (True, False),
                         ("   ", "\t\n"), (object(), object())):
            assert sc._plasmidsaurus_credential_hint(cid, sec) == ""

    def test_a_value_cannot_forge_a_log_line_or_emit_escapes(self):
        """The hint reaches both the log and the terminal, and it is built
        from an attacker-influenced-ish value (whatever is in Settings)."""
        for bad in ("b" * 32 + "\nFAKE LOG LINE", "\x1b[31mred", "\r\n" * 5):
            h = sc._plasmidsaurus_credential_hint("a" * 32, bad)
            assert "\n" not in h and "\r" not in h and "\x1b" not in h
            assert bad not in h

    def test_an_oversized_value_reports_a_length_not_the_value(self):
        h = sc._plasmidsaurus_credential_hint("a" * 32, "b" * 1_000_000)
        assert "1000000 characters" in h and len(h) < 200

    def test_the_agent_endpoint_surfaces_the_hint_too(self, monkeypatch):
        """One placement, every surface: the GUI modal, the download path and
        the agent endpoints all obtain their token through
        `_plasmidsaurus_oauth_token`, so the hint rides along. Pinned rather
        than asserted, since a future caller could bypass it."""
        import splicecraft_agent as _ag
        # The moved handler resolves its deps in the SIBLING's namespace.
        monkeypatch.setattr(_ag, "_plasmidsaurus_credentials",
                            lambda: ("seb@example.org", "hunter2"))
        self._raise_status(monkeypatch, 401)
        handler = sc._state._AGENT_HANDLERS["plasmidsaurus-items"][0]
        body, status = handler(None, {})
        assert status == 502
        assert "email address" in body["error"]

    def test_a_429_is_still_not_blamed_on_the_credentials(self, monkeypatch):
        """The pre-existing carve-out must survive: a rate limit says nothing
        about the key, so it must not pick up a shape hint either."""
        self._raise_status(monkeypatch, 429)
        with pytest.raises(OSError) as ei:
            sc._plasmidsaurus_oauth_token("seb@example.org", "hunter2")
        assert "email address" not in str(ei.value)


# ── OAuth + JSON API ─────────────────────────────────────────────────────────
class TestApiClient:
    def test_token_happy(self, use_router):
        use_router(_Router(token="abc123"))
        assert sc._plasmidsaurus_oauth_token("cid", "sec") == "abc123"

    def test_token_bad_credentials_raises_oserror(self, monkeypatch):
        import urllib.error

        class _O:
            def open(self, req, timeout=None):
                raise urllib.error.HTTPError(
                    req.get_full_url(), 401, "Unauthorized", {}, None)
        monkeypatch.setattr(_search, "_build_hardened_url_opener", lambda: _O())
        with pytest.raises(OSError, match="credential"):
            sc._plasmidsaurus_oauth_token("cid", "bad")

    def test_token_requires_both(self):
        with pytest.raises(ValueError):
            sc._plasmidsaurus_oauth_token("", "sec")

    def test_api_get_404_raises(self, monkeypatch):
        import urllib.error

        class _O:
            def open(self, req, timeout=None):
                raise urllib.error.HTTPError(
                    req.get_full_url(), 404, "Not Found", {}, None)
        monkeypatch.setattr(_search, "_build_hardened_url_opener", lambda: _O())
        with pytest.raises(OSError, match="404"):
            sc._plasmidsaurus_api_get("/api/item/ABCDEF", "tok")

    def test_api_get_oversize_raises(self, monkeypatch):
        big = b"[" + b"0," * 200 + b"0]"

        class _O:
            def open(self, req, timeout=None):
                return _resp(big)
        monkeypatch.setattr(_search, "_build_hardened_url_opener", lambda: _O())
        with pytest.raises(ValueError, match="too large"):
            sc._plasmidsaurus_api_get("/api/items", "tok", max_bytes=8)

    def test_list_items(self, use_router):
        use_router(_Router(items=[{"code": "ABCDEF", "status": "complete"}]))
        items = sc._plasmidsaurus_list_items("tok")
        assert items and items[0]["code"] == "ABCDEF"

    def test_result_link(self, use_router):
        use_router(_Router(link="https://files/x.zip"))
        assert sc._plasmidsaurus_result_link("tok", "ABCDEF") == \
            "https://files/x.zip"

    def test_result_link_bad_kind(self):
        with pytest.raises(ValueError, match="kind"):
            sc._plasmidsaurus_result_link("tok", "ABCDEF", kind="genome")

    def test_result_link_missing_link(self, monkeypatch):
        class _O:
            def open(self, req, timeout=None):
                return _resp(b"{}")
        monkeypatch.setattr(_search, "_build_hardened_url_opener", lambda: _O())
        with pytest.raises(ValueError, match="no results download link"):
            sc._plasmidsaurus_result_link("tok", "ABCDEF")


# ── live-API behaviours (found against a real account, 2026-07-27) ───────────
class TestLiveApiQuirks:
    """Server behaviour the mocked happy path never exercised, all of it
    confirmed against a live Plasmidsaurus account: `/api/items` silently
    truncating at its 100-item default (which hid 118 of 218 real orders),
    the 10-request/minute rate limit, and the shipping-label / canceled
    orders that sit in the listing marked `complete` with nothing to fetch.
    """

    def test_list_items_sends_explicit_limit(self, use_router):
        r = use_router(_Router(items=[]))
        sc._plasmidsaurus_list_items("tok")
        assert f"limit={sc._PLASMIDSAURUS_ITEMS_LIMIT}" in r.calls[-1]

    def test_list_items_shared_keeps_both_params(self, use_router):
        r = use_router(_Router(items=[]))
        sc._plasmidsaurus_list_items("tok", shared=True)
        assert "shared=true" in r.calls[-1]
        assert f"limit={sc._PLASMIDSAURUS_ITEMS_LIMIT}" in r.calls[-1]

    def test_list_items_clamps_limit_to_server_max(self, use_router):
        # The server 400s a limit above 1000 rather than capping it, so an
        # over-large request must be clamped client-side.
        r = use_router(_Router(items=[]))
        sc._plasmidsaurus_list_items("tok", limit=99999)
        assert f"limit={sc._PLASMIDSAURUS_ITEMS_LIMIT}" in r.calls[-1]
        sc._plasmidsaurus_list_items("tok", limit=0)
        assert "limit=1" in r.calls[-1]

    @staticmethod
    def _raise_status(monkeypatch, code: int):
        import urllib.error

        class _O:
            def open(self, req, timeout=None):
                raise urllib.error.HTTPError(
                    req.get_full_url(), code, "err", {}, None)
        monkeypatch.setattr(_search, "_build_hardened_url_opener", lambda: _O())

    def test_token_429_says_rate_limit_not_bad_credentials(self, monkeypatch):
        # 429 shares the 4xx space with the credential rejections; reporting
        # it as "check your client ID" sends the user debugging a good key.
        self._raise_status(monkeypatch, 429)
        with pytest.raises(OSError) as ei:
            sc._plasmidsaurus_oauth_token("cid", "sec")
        assert "rate limit" in str(ei.value).lower()
        assert "credential" not in str(ei.value).lower()
        assert getattr(ei.value, "http_status", None) == 429

    def test_api_get_429_says_rate_limit(self, monkeypatch):
        self._raise_status(monkeypatch, 429)
        with pytest.raises(OSError, match="rate limit"):
            sc._plasmidsaurus_api_get("/api/items", "tok")

    def test_api_get_tags_http_status(self, monkeypatch):
        self._raise_status(monkeypatch, 404)
        with pytest.raises(OSError) as ei:
            sc._plasmidsaurus_api_get("/api/item/ABCDEF/results", "tok")
        assert getattr(ei.value, "http_status", None) == 404

    def test_result_link_404_names_the_real_cause(self, monkeypatch):
        self._raise_status(monkeypatch, 404)
        with pytest.raises(OSError) as ei:
            sc._plasmidsaurus_result_link("tok", "ABCDEF")
        msg = str(ei.value).lower()
        assert "no results to download" in msg
        assert "canceled" in msg and "shipping" in msg

    def test_result_link_non_404_passes_through(self, monkeypatch):
        self._raise_status(monkeypatch, 403)
        with pytest.raises(OSError, match="not authorised"):
            sc._plasmidsaurus_result_link("tok", "ABCDEF")

    @pytest.mark.parametrize("item,expected", [
        ({"status": "complete", "done_date": "2026-06-01",
          "product_name": "plasmid_low_copy"}, True),
        # Shipping label: marked complete, sorts to the top, 404s on /results.
        ({"status": "complete", "done_date": None,
          "product_name": "ups_shipping_label"}, False),
        ({"status": "canceled", "done_date": "2026-06-01",
          "product_name": "plasmid_low_copy"}, False),
        ({"status": "complete", "done_date": None,
          "product_name": "plasmid_low_copy"}, False),
        ({"status": "COMPLETE", "done_date": "2026-06-01",
          "product_name": "Plasmid_Low_Copy"}, True),
        ({}, False),
        (None, False),
        ("not-a-dict", False),
    ])
    def test_has_results_heuristic(self, item, expected):
        assert sc._plasmidsaurus_item_has_results(item) is expected


# ── zip download (PK magic + caps + content-type) ────────────────────────────
class TestDownloadZip:
    def _opener_returning(self, monkeypatch, body, ctype="application/zip"):
        class _O:
            def open(self, req, timeout=None):
                return _resp(body, content_type=ctype)
        monkeypatch.setattr(_search, "_build_hardened_url_opener", lambda: _O())

    def test_accepts_real_zip(self, tmp_path, monkeypatch):
        zb = _zip_bytes()
        self._opener_returning(monkeypatch, zb)
        dest = tmp_path / "ABCDEF_results.zip"
        sha = sc._plasmidsaurus_download_zip(
            "https://files/x.zip", dest, max_bytes=10 * 1024 * 1024)
        assert dest.exists() and dest.read_bytes() == zb and len(sha) == 64

    def test_rejects_non_zip(self, tmp_path, monkeypatch):
        self._opener_returning(monkeypatch, b"NOT A ZIP" * 50,
                               ctype="application/octet-stream")
        with pytest.raises(ValueError, match="zip"):
            sc._plasmidsaurus_download_zip(
                "https://files/x.zip", tmp_path / "x.zip",
                max_bytes=10 * 1024 * 1024)
        assert not (tmp_path / "x.zip").exists()   # nothing left on disk

    def test_rejects_html_error_page(self, tmp_path, monkeypatch):
        self._opener_returning(monkeypatch, b"<html>blocked</html>",
                               ctype="text/html")
        with pytest.raises(ValueError, match="Content-Type"):
            sc._plasmidsaurus_download_zip(
                "https://files/x.zip", tmp_path / "x.zip",
                max_bytes=10 * 1024 * 1024)

    def test_enforces_cap(self, tmp_path, monkeypatch):
        self._opener_returning(monkeypatch, _zip_bytes())
        with pytest.raises(ValueError, match="cap"):
            sc._plasmidsaurus_download_zip(
                "https://files/x.zip", tmp_path / "x.zip", max_bytes=8)

    @pytest.mark.parametrize("code", [403, 404])
    def test_expired_presigned_link_says_so(self, tmp_path, monkeypatch, code):
        # The results URL is a short-lived pre-signed S3 link; the object
        # store answers an expired signature with 403/404, which as a bare
        # status reads like an account-permission problem.
        import urllib.error

        class _O:
            def open(self, req, timeout=None):
                raise urllib.error.HTTPError(
                    req.get_full_url(), code, "err", {}, None)
        monkeypatch.setattr(_search, "_build_hardened_url_opener", lambda: _O())
        with pytest.raises(OSError) as ei:
            sc._plasmidsaurus_download_zip(
                "https://files/x.zip", tmp_path / "x.zip",
                max_bytes=10 * 1024 * 1024)
        msg = str(ei.value).lower()
        assert "pre-signed" in msg and "short-lived" in msg
        assert getattr(ei.value, "http_status", None) == code

    def test_rejects_non_https(self, tmp_path):
        with pytest.raises(ValueError, match="HTTPS"):
            sc._plasmidsaurus_download_zip(
                "http://files/x.zip", tmp_path / "x.zip",
                max_bytes=10 * 1024 * 1024)


# ── zip → library entries ────────────────────────────────────────────────────
class TestZipToEntries:
    def test_builds_entries(self, tmp_path):
        zp = tmp_path / "run.zip"
        zp.write_bytes(_zip_bytes(("alpha.gbk", "beta.gbk")))
        entries, warnings = _fileio._plasmidsaurus_zip_to_entries(
            zp, run_id="ABCDEF")
        assert warnings == []
        assert {e["name"] for e in entries} == {"alpha", "beta"}
        for e in entries:
            assert e["source"].startswith("plasmidsaurus:ABCDEF:")
            assert e["status"] == "" and e["gb_text"] and e["size"] > 0
            assert sc._gb_text_to_record(e["gb_text"]) is not None

    def test_bad_member_becomes_warning(self, tmp_path):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as zf:
            zf.writestr("good.gbk", _mkgb("g"))
            zf.writestr("broken.gbk", "this is not genbank at all")
        zp = tmp_path / "run.zip"
        zp.write_bytes(buf.getvalue())
        entries, warnings = _fileio._plasmidsaurus_zip_to_entries(zp)
        assert [e["name"] for e in entries] == ["good"]
        assert len(warnings) == 1 and "broken.gbk" in warnings[0]

    def test_non_zip_raises(self, tmp_path):
        bad = tmp_path / "x.zip"
        bad.write_bytes(b"not a zip")
        with pytest.raises(ValueError):
            _fileio._plasmidsaurus_zip_to_entries(bad)


# ── fetch orchestration ──────────────────────────────────────────────────────
class TestFetchItemZip:
    def test_end_to_end(self, tmp_path, use_router):
        use_router(_Router())
        out = sc._plasmidsaurus_fetch_item_zip(
            "abcdef", tmp_path, client_id="cid", client_secret="sec")
        assert out.name == "ABCDEF_results.zip" and out.exists()

    def test_bad_code_raises(self, tmp_path):
        with pytest.raises(ValueError, match="item code"):
            sc._plasmidsaurus_fetch_item_zip(
                "bad/x", tmp_path, client_id="c", client_secret="s")

    def test_no_credentials_raises(self, tmp_path, monkeypatch):
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_ID", raising=False)
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_SECRET", raising=False)
        sc._set_setting("plasmidsaurus_client_id", "")
        sc._set_setting("plasmidsaurus_client_secret", "")
        with pytest.raises(ValueError, match="credential"):
            sc._plasmidsaurus_fetch_item_zip("ABCDEF", tmp_path)


# ── agent endpoints ──────────────────────────────────────────────────────────
class TestAgentEndpoints:
    def _items(self):
        return sc._state._AGENT_HANDLERS["plasmidsaurus-items"][0]

    def _download(self):
        return sc._state._AGENT_HANDLERS["download-plasmidsaurus"][0]

    def test_registered(self):
        H = sc._state._AGENT_HANDLERS
        assert H["plasmidsaurus-items"][1] is False     # read
        assert H["download-plasmidsaurus"][1] is True   # write

    def test_items_no_creds_400(self, monkeypatch):
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_ID", raising=False)
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_SECRET", raising=False)
        sc._set_setting("plasmidsaurus_client_id", "")
        sc._set_setting("plasmidsaurus_client_secret", "")
        payload, status = self._items()(None, {})
        assert status == 400 and "credential" in payload["error"].lower()

    def test_items_happy(self, monkeypatch, use_router):
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        use_router(_Router(items=[{"code": "ABCDEF", "status": "complete",
                                   "product_name": "plasmid_high_copy",
                                   "quantity": 2, "done_date": "2026-06-01",
                                   "gross": 30.0}]))
        r = self._items()(None, {})
        assert r["ok"] and r["count"] == 1
        assert r["items"][0]["code"] == "ABCDEF"

    def test_items_projects_order_name_and_flags_has_results(
            self, monkeypatch, use_router):
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        use_router(_Router(items=[
            {"code": "AAAAAA", "status": "complete", "order_name": "my run",
             "product_name": "plasmid_low_copy", "quantity": 2,
             "done_date": "2026-06-01", "gross": 30.0},
            {"code": "BBBBBB", "status": "complete", "order_name": "label",
             "product_name": "ups_shipping_label", "quantity": 1,
             "done_date": None, "gross": 12.0},
        ]))
        r = self._items()(None, {})
        assert r["count"] == 2 and r["downloadable"] == 1
        assert r["items"][0]["order_name"] == "my run"
        assert r["items"][0]["has_results"] is True
        # The shipping label is KEPT in the listing (dropping rows would break
        # the agent-API compat promise) but flagged so a picker can skip it.
        assert r["items"][1]["has_results"] is False
        assert "truncated" not in r

    def test_items_flags_truncation_at_the_server_ceiling(
            self, monkeypatch, use_router):
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        n = sc._PLASMIDSAURUS_ITEMS_LIMIT
        use_router(_Router(items=[{"code": f"{i:06d}"} for i in range(n)]))
        r = self._items()(None, {})
        assert r["count"] == n and r["truncated"] is True
        assert "pagination" in r["warning"]

    def test_truncation_measured_on_the_raw_server_list(
            self, monkeypatch, use_router):
        # A junk row drops out of the projection, so counting the PROJECTED
        # rows would put a genuinely-full page under the ceiling and retract
        # the warning on exactly the listing that needed it.
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        n = sc._PLASMIDSAURUS_ITEMS_LIMIT
        items = [{"code": f"{i:06d}"} for i in range(n - 2)] + ["junk", None]
        use_router(_Router(items=items))
        r = self._items()(None, {})
        assert r["count"] == n - 2          # junk rows dropped from output
        assert r["truncated"] is True       # but the page WAS full

    def test_download_bad_code_400(self):
        payload, status = self._download()(None, {"item_code": "bad/x"})
        assert status == 400 and "item_code" in payload["error"]

    def test_download_non_results_kind_400(self):
        payload, status = self._download()(
            None, {"item_code": "ABCDEF", "kind": "reads"})
        assert status == 400 and "results" in payload["error"].lower()

    def test_download_no_creds_400(self, monkeypatch):
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_ID", raising=False)
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_SECRET", raising=False)
        sc._set_setting("plasmidsaurus_client_id", "")
        sc._set_setting("plasmidsaurus_client_secret", "")
        payload, status = self._download()(None, {"item_code": "ABCDEF"})
        assert status == 400 and "credential" in payload["error"].lower()

    def test_download_imports_into_library(self, monkeypatch, use_router):
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        use_router(_Router(zipbytes=_zip_bytes(("foo.gbk", "bar.gbk"))))
        before = len(sc._load_library())
        r = self._download()(None, {"item_code": "abcdef"})
        assert r["ok"] and r["n_added"] == 2
        assert {a["name"] for a in r["added"]} == {"foo", "bar"}
        lib = sc._load_library()
        assert len(lib) == before + 2
        tagged = [e for e in lib if e.get("source", "").startswith(
            "plasmidsaurus:ABCDEF:")]
        assert len(tagged) == 2
        # Re-import APPENDS (never overwrites / drops) → +2 more.
        r2 = self._download()(None, {"item_code": "abcdef"})
        assert r2["n_added"] == 2 and len(sc._load_library()) == before + 4
        # …and every library id stays UNIQUE. Library `id` is the canonical
        # key (delete-by-id filters on it), so a duplicate would make a later
        # single-entry delete nuke both copies — a data-loss class bug.
        ids = [e["id"] for e in sc._load_library() if "id" in e]
        assert len(ids) == len(set(ids)), "duplicate library ids after re-import"


# ── secret never reaches the log / event stream ──────────────────────────────
class TestSecretRedaction:
    def test_secret_redacted(self, monkeypatch):
        import splicecraft_dataaccess as _da
        evs = []
        monkeypatch.setattr(_da, "_log_event",
                            lambda ev, **kw: evs.append((ev, kw)))
        import logging
        records = []

        class _Cap(logging.Handler):
            def emit(self, record):
                records.append(record.getMessage())
        cap = _Cap()
        sc._log.addHandler(cap)
        sc._log.setLevel(logging.DEBUG)
        try:
            secret = "sk-DO-NOT-LOG-ME-123abc"
            sc._set_setting("plasmidsaurus_client_secret", secret)
            blob = "\n".join(records) + repr(evs)
            assert secret not in blob, "SECRET LEAKED"
            assert any(kw.get("key") == "plasmidsaurus_client_secret"
                       and kw.get("value") == "<redacted>" for _, kw in evs)
        finally:
            sc._log.removeHandler(cap)


# ── secret stays OUT of the agent settings surface ───────────────────────────
class TestSecretExcludedFromAgentSurface:
    """The CHANGELOG promises a remote agent "can neither read it back nor
    change it." Pin that: the credential keys are absent from the settings
    allowlist, so get-settings omits the secret and set-setting refuses it."""

    def test_not_in_allowlist(self):
        assert "plasmidsaurus_client_secret" not in sc._AGENT_SETTINGS_ALLOWLIST
        assert "plasmidsaurus_client_id" not in sc._AGENT_SETTINGS_ALLOWLIST

    def test_get_settings_omits_secret(self):
        get = sc._state._AGENT_HANDLERS["get-settings"][0]
        r = get(None, {})
        assert "plasmidsaurus_client_secret" not in r["settings"]

    def test_set_setting_refuses_secret(self):
        setter = sc._state._AGENT_HANDLERS["set-setting"][0]
        payload, status = setter(
            None, {"key": "plasmidsaurus_client_secret", "value": "leak"})
        assert status == 400 and "unknown setting" in payload["error"]
        assert sc._get_setting("plasmidsaurus_client_secret", "") != "leak"


# ── GUI: SettingsModal credential fields ─────────────────────────────────────
class TestSettingsModalCreds:
    async def test_save_and_clear(self, monkeypatch):
        from textual.widgets import Input
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_ID", raising=False)
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_SECRET", raising=False)
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.SettingsModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#set-ps-id", Input).value = "my-id"
            modal.query_one("#set-ps-secret", Input).value = "my-secret"
            modal._ps_save_creds(None)
            await pilot.pause()
            assert sc._get_setting("plasmidsaurus_client_id") == "my-id"
            assert sc._get_setting("plasmidsaurus_client_secret") == "my-secret"
            modal._ps_clear_creds(None)
            await pilot.pause()
            assert sc._get_setting("plasmidsaurus_client_id") == ""
            assert sc._get_setting("plasmidsaurus_client_secret") == ""
            app.exit()


# ── GUI: PlasmidsaurusFetchModal ─────────────────────────────────────────────
class TestOrderRows:
    """`_plasmidsaurus_order_rows` — pure, so the selection logic is testable
    without driving a DataTable (headless DataTable has its own quirks)."""

    def test_hides_orders_with_nothing_to_download(self):
        rows = sc._plasmidsaurus_order_rows(_ORDERS)
        assert [r["code"] for r in rows] == ["AAAAAA"]

    def test_include_empty_shows_all_and_marks_them(self):
        rows = sc._plasmidsaurus_order_rows(_ORDERS, include_empty=True)
        assert [r["code"] for r in rows] == ["AAAAAA", "BBBBBB", "CCCCCC"]
        assert [r["has_results"] for r in rows] == [True, False, False]

    def test_date_uses_the_universal_slash_free_form(self):
        rows = sc._plasmidsaurus_order_rows(_ORDERS)
        assert rows[0]["date"] == "JUN 1 2026"

    def test_undated_order_shows_a_dash_not_a_bogus_date(self):
        rows = sc._plasmidsaurus_order_rows(_ORDERS, include_empty=True)
        assert rows[1]["date"] == "—"

    def test_product_underscores_become_spaces(self):
        rows = sc._plasmidsaurus_order_rows(_ORDERS)
        assert rows[0]["product"] == "plasmid low copy"

    def test_unusable_codes_are_dropped(self):
        bad = [{"code": "toolongcode", "status": "complete",
                "done_date": "2026-06-01", "product_name": "plasmid"},
               {"code": None, "status": "complete",
                "done_date": "2026-06-01", "product_name": "plasmid"}]
        assert sc._plasmidsaurus_order_rows(bad, include_empty=True) == []

    def test_order_name_is_sanitised(self):
        # Order names are remote input rendered into a terminal.
        evil = [{"code": "AAAAAA", "status": "complete",
                 "order_name": "boom\x1b[31mred", "product_name": "plasmid",
                 "done_date": "2026-06-01", "quantity": 1}]
        rows = sc._plasmidsaurus_order_rows(evil)
        assert "\x1b" not in rows[0]["order"]

    def test_duplicate_codes_are_collapsed(self):
        """The code is the DataTable row key and Textual raises DuplicateKey
        on a repeat, so a duplicated server row must not reach the table."""
        dupes = _ORDERS + [dict(_ORDERS[0], order_name="same code again")]
        rows = sc._plasmidsaurus_order_rows(dupes, include_empty=True)
        codes = [r["code"] for r in rows]
        assert len(codes) == len(set(codes))
        # Newest-first: the first occurrence is the one kept.
        assert rows[0]["order"] == "run one"

    def test_non_dict_rows_are_skipped(self):
        assert sc._plasmidsaurus_order_rows(
            ["junk", None, 7], include_empty=True) == []

    def test_missing_quantity_shows_a_dash(self):
        one = [{"code": "AAAAAA", "status": "complete",
                "done_date": "2026-06-01", "product_name": "plasmid"}]
        assert sc._plasmidsaurus_order_rows(one)[0]["qty"] == "—"


class TestOrdersCache:
    """The account is capped at 10 requests/minute, so the listing is cached."""

    def test_second_call_costs_no_requests(self, use_router):
        r = use_router(_Router(items=_ORDERS))
        items_a, tok_a = sc._plasmidsaurus_orders("cid", "sec")
        spent = len(r.calls)
        items_b, tok_b = sc._plasmidsaurus_orders("cid", "sec")
        assert len(r.calls) == spent          # nothing new went out
        assert tok_a == tok_b
        assert [i["code"] for i in items_a] == [i["code"] for i in items_b]

    def test_force_refetches(self, use_router):
        r = use_router(_Router(items=_ORDERS))
        sc._plasmidsaurus_orders("cid", "sec")
        spent = len(r.calls)
        sc._plasmidsaurus_orders("cid", "sec", force=True)
        assert len(r.calls) > spent

    def test_different_credentials_miss_the_cache(self, use_router):
        # Otherwise changing the key in Settings would keep showing the
        # previous account's orders.
        r = use_router(_Router(items=_ORDERS))
        sc._plasmidsaurus_orders("cid", "sec")
        spent = len(r.calls)
        sc._plasmidsaurus_orders("other", "creds")
        assert len(r.calls) > spent

    def test_expired_entry_refetches(self, use_router):
        r = use_router(_Router(items=_ORDERS))
        sc._plasmidsaurus_orders("cid", "sec")
        spent = len(r.calls)
        sc._PS_ORDERS_CACHE["at"] = (
            sc._monotonic() - sc._PS_ORDERS_CACHE_TTL_S - 1.0)
        sc._plasmidsaurus_orders("cid", "sec")
        assert len(r.calls) > spent

    def test_raw_credentials_are_never_stored(self, use_router):
        use_router(_Router(items=_ORDERS))
        sc._plasmidsaurus_orders("cid-secret-value", "sec-secret-value")
        blob = repr(sc._PS_ORDERS_CACHE)
        assert "cid-secret-value" not in blob
        assert "sec-secret-value" not in blob


class TestFetchModal:
    async def test_invalid_code_shows_error(self, monkeypatch, use_router):
        from textual.widgets import Input, Static
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        # Credentials present means the modal auto-enumerates on mount, so the
        # opener MUST be stubbed or this test would reach the real API.
        use_router(_Router(items=_ORDERS))
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#ps-fetch-code", Input).value = "xx"
            modal._start()
            await pilot.pause()
            status = str(modal.query_one("#ps-fetch-status", Static).render())
            assert "valid" in status.lower() and modal._busy is False
            app.exit()

    async def test_done_callback_refreshes_and_reports(self):
        from textual.widgets import Static
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            entries = [{"name": "s1"}, {"name": "s2"}]
            modal._busy = True
            modal._fetch_done(("ok", entries, ["skipped.gbk: bad"]))
            await pilot.pause()
            status = str(modal.query_one("#ps-fetch-status", Static).render())
            assert "Imported 2" in status and "1 skipped" in status
            assert modal._busy is False
            app.exit()

    async def test_open_enumerates_orders(self, monkeypatch, use_router):
        from textual.widgets import DataTable, Static
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        use_router(_Router(items=_ORDERS))
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            # Drive the callback directly rather than racing the thread worker.
            modal._list_done(("ok", _ORDERS, "tok"))
            await pilot.pause()
            tbl = modal.query_one("#ps-fetch-table", DataTable)
            assert tbl.row_count == 1          # only the downloadable order
            assert modal._token == "tok"       # reused for the download
            status = str(modal.query_one("#ps-fetch-status", Static).render())
            assert "1 of 3" in status
            app.exit()

    async def test_show_empty_checkbox_reveals_the_rest(self, monkeypatch,
                                                        use_router):
        from textual.widgets import Checkbox, DataTable
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        use_router(_Router(items=_ORDERS))
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            modal._list_done(("ok", _ORDERS, "tok"))
            await pilot.pause()
            modal.query_one("#ps-show-empty", Checkbox).value = True
            await pilot.pause()
            assert modal.query_one("#ps-fetch-table", DataTable).row_count == 3
            app.exit()

    async def test_selecting_a_no_results_row_refuses_before_requesting(
            self, monkeypatch, use_router):
        from textual.widgets import Checkbox, Static
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        use_router(_Router(items=_ORDERS))
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            modal._list_done(("ok", _ORDERS, "tok"))
            modal.query_one("#ps-show-empty", Checkbox).value = True
            await pilot.pause()

            class _Ev:                      # only `.row_key.value` is read
                class row_key:              # noqa: N801
                    value = "BBBBBB"        # the shipping label
            modal._row_selected(_Ev())
            await pilot.pause()
            status = str(modal.query_one("#ps-fetch-status", Static).render())
            assert "nothing to download" in status
            assert modal._busy is False     # no download was started
            app.exit()

    async def test_listing_does_not_steal_focus_from_the_code_box(
            self, monkeypatch, use_router):
        from textual.widgets import DataTable, Input
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        use_router(_Router(items=_ORDERS))
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            # User starts typing a shared code before the listing lands.
            box = modal.query_one("#ps-fetch-code", Input)
            box.focus()
            await pilot.pause()
            modal._list_done(("ok", _ORDERS, "tok"))
            await pilot.pause()
            assert box.has_focus, "listing yanked focus out of the code box"
            app.exit()

    async def test_listing_focuses_the_table_when_code_box_is_idle(
            self, monkeypatch, use_router):
        from textual.widgets import DataTable
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        use_router(_Router(items=_ORDERS))
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            modal._list_done(("ok", _ORDERS, "tok"))
            await pilot.pause()
            assert modal.query_one("#ps-fetch-table", DataTable).has_focus
            app.exit()

    async def test_listing_failure_is_reported_not_swallowed(self, monkeypatch,
                                                             use_router):
        from textual.widgets import Static
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_ID", "cid")
        monkeypatch.setenv("PLASMIDSAURUS_CLIENT_SECRET", "sec")
        use_router(_Router(items=_ORDERS))
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            modal._list_done(("err", "rate limit reached", ""))
            await pilot.pause()
            status = str(modal.query_one("#ps-fetch-status", Static).render())
            assert "rate limit" in status and modal._listing is False
            app.exit()

    async def test_no_credentials_does_not_enumerate(self, monkeypatch):
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_ID", raising=False)
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_SECRET", raising=False)
        sc._set_setting("plasmidsaurus_client_id", "")
        sc._set_setting("plasmidsaurus_client_secret", "")
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            # No creds must mean no request attempt at all.
            assert modal._listing is False and modal._items == []
            app.exit()

    async def test_creds_hint_when_unset(self, monkeypatch):
        from textual.widgets import Static
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_ID", raising=False)
        monkeypatch.delenv("PLASMIDSAURUS_CLIENT_SECRET", raising=False)
        sc._set_setting("plasmidsaurus_client_id", "")
        sc._set_setting("plasmidsaurus_client_secret", "")
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PlasmidsaurusFetchModal())
            await pilot.pause()
            modal = app.screen
            hint = str(modal.query_one("#ps-fetch-creds", Static).render())
            assert "No API credentials" in hint
            app.exit()
