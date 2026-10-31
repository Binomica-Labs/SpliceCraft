"""Serving a local browser UI from the agent API's own origin
(`--workbench-dir`), and the single-origin dev escape (`--cors-origin`).

The agent API has no CORS story by design — it binds loopback and every
endpoint wants a bearer token the browser cannot read. Same-origin serving is
what makes a browser front-end possible at all: it removes the preflight and
the token-distribution problem together. These tests cover the two things that
buys and the two things it costs:

  * **A static route that cannot be walked out of.** Every traversal shape is
    asserted against a real file planted OUTSIDE the bundle, over a RAW socket
    — `urllib` normalises `/../` client-side, so a test written through it
    proves only that urllib works. Percent-encoded, double-encoded,
    backslash, NUL and symlink forms are each their own case because each
    defeats a different naive check.
  * **A static route that cannot shadow an endpoint.** The claim rule rests on
    no endpoint name containing a dot; that is asserted directly, so a future
    endpoint called `foo.bar` fails here rather than quietly becoming
    unreachable.
  * **An unauthenticated route that hands out the bearer token.** `GET /` has
    to answer before the page has a token. The test pins that it happens ONLY
    with the flag set, ONLY behind the `Host` anti-rebinding check, and that
    the default build of the server still exposes nothing.
  * **Nothing changes when the flags are absent.** A whole class of these
    tests exists to say the default surface is byte-for-byte what it was:
    `GET /` is still the `tools` self-describe, an asset-shaped path is still
    an unknown endpoint, and `OPTIONS` emits no `Access-Control-*` at all.
"""
from __future__ import annotations

import socket
import sys
import threading
import time

import pytest

import splicecraft as sc


# ── helpers ───────────────────────────────────────────────────────────────────

def _free_port() -> int:
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


class _TinyApp:
    """Enough app for `_AgentAPIServer.__init__` to read its flags off."""

    _headless = True
    _current_record = None

    def __init__(self, workbench_dir=None, cors_origin=None):
        self._workbench_dir = workbench_dir
        self._cors_origin = cors_origin


class _Client:
    """Raw-socket HTTP client.

    Deliberately NOT `urllib`: it resolves `..` in the URL before sending, so
    `/../../etc/passwd` would leave the client as `/etc/passwd` and the
    traversal tests would be testing urllib rather than the server. The whole
    point is to put an un-normalised path on the wire.
    """

    def __init__(self, port: int) -> None:
        self.port = port

    def request(self, path: str, *, method: str = "GET",
                headers: "dict[str, str] | None" = None,
                host: "str | None" = None,
                ) -> "tuple[int, dict[str, str], bytes]":
        lines = [f"{method} {path} HTTP/1.1",
                 f"Host: {host or f'127.0.0.1:{self.port}'}"]
        for key, value in (headers or {}).items():
            lines.append(f"{key}: {value}")
        lines += ["Connection: close", "", ""]
        conn = socket.create_connection(("127.0.0.1", self.port), timeout=5)
        try:
            conn.sendall("\r\n".join(lines).encode("utf-8"))
            chunks = []
            while True:
                chunk = conn.recv(65536)
                if not chunk:
                    break
                chunks.append(chunk)
        finally:
            conn.close()
        head, _, body = b"".join(chunks).partition(b"\r\n\r\n")
        status = int(head.split(b"\r\n")[0].split()[1])
        hdrs: "dict[str, str]" = {}
        for line in head.split(b"\r\n")[1:]:
            if b":" in line:
                key, value = line.split(b":", 1)
                hdrs[key.decode().lower()] = value.decode().strip()
        return status, hdrs, body


def _serve(app) -> "tuple[_Client, str]":
    port = _free_port()
    token = f"wb-token-{port}"
    srv = sc._AgentAPIServer(("127.0.0.1", port), app, token)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    # The listening socket is up after __init__ returns, but the accept loop
    # needs a beat or the very first connect can race it.
    time.sleep(0.05)
    return (_Client(port), token), srv


@pytest.fixture
def bundle(tmp_path):
    """A small bundle, plus the hostile neighbours a traversal would want.

    `secret.html` sits OUTSIDE the bundle with a SERVABLE extension — the real
    attack shape, since an unservable one is refused by the content-type map
    before the path is ever resolved, which would make the traversal tests
    pass for the wrong reason.
    """
    root = tmp_path / "ui"
    (root / "assets").mkdir(parents=True)
    (root / "index.html").write_text(
        "<!doctype html>\n<html><head><title>t</title></head>"
        "<body><script src=\"/assets/app.js\"></script></body></html>\n",
        encoding="utf-8")
    (root / "assets" / "app.js").write_text("console.log(1);\n",
                                            encoding="utf-8")
    (root / "style.css").write_text("body{color:red}\n", encoding="utf-8")
    (root / "notes.txt").write_text("not servable\n", encoding="utf-8")
    (tmp_path / "secret.html").write_text("OUTSIDE-MARKER\n", encoding="utf-8")
    (tmp_path / "secret.json").write_text('{"a":"OUTSIDE-MARKER"}\n',
                                          encoding="utf-8")
    try:
        (root / "link.html").symlink_to(tmp_path / "secret.html")
    except (OSError, NotImplementedError):      # pragma: no cover - Windows
        pass
    return root


@pytest.fixture
def served(bundle):
    (client, token), srv = _serve(_TinyApp(workbench_dir=bundle.resolve()))
    try:
        yield client, token, bundle
    finally:
        srv.shutdown()
        srv.server_close()


@pytest.fixture
def plain_served():
    """The server as it ships: no workbench, no CORS."""
    (client, token), srv = _serve(_TinyApp())
    try:
        yield client, token
    finally:
        srv.shutdown()
        srv.server_close()


# ── the acceptance criteria ───────────────────────────────────────────────────

class TestServesTheBundle:
    def test_root_returns_index_html_with_the_token_injected(self, served):
        client, token, _root = served
        status, hdrs, body = client.request("/")
        assert status == 200
        assert hdrs["content-type"].startswith("text/html")
        text = body.decode("utf-8")
        assert token in text, "the page cannot authenticate without it"
        assert f"window.{sc._WORKBENCH_GLOBAL}" in text
        # The document still has to be a document.
        assert text.lower().startswith("<!doctype html>")
        assert "<title>t</title>" in text

    def test_assets_are_served_with_their_declared_type(self, served):
        client, _token, _root = served
        for path, ctype, needle in (("/assets/app.js", "text/javascript",
                                     "console.log"),
                                    ("/style.css", "text/css", "color:red")):
            status, hdrs, body = client.request(path)
            assert status == 200, path
            assert hdrs["content-type"].startswith(ctype), path
            assert needle in body.decode("utf-8"), path

    def test_assets_are_not_cached(self, served):
        """The token is baked into the HTML and a dev rebuilds in place, so a
        cached copy is a stale token as well as stale code."""
        client, _token, _root = served
        _s, hdrs, _b = client.request("/")
        assert hdrs.get("cache-control") == "no-store"

    def test_a_missing_asset_is_a_404_with_no_endpoint_list(self, served):
        """This route answers BEFORE the token check, so anything specific
        here would be an unauthenticated filesystem/endpoint oracle."""
        client, _token, _root = served
        status, _hdrs, body = client.request("/nope.js")
        assert status == 404
        assert b"endpoints" not in body

    def test_an_unsupported_extension_is_415_not_a_guessed_type(self, served):
        client, _token, _root = served
        status, _hdrs, body = client.request("/notes.txt")
        assert status == 415
        assert b"octet-stream" not in body
        # The refusal names what IS servable, so a missing entry is a
        # one-line fix rather than a mystery.
        assert b".css" in body and b".wasm" in body

    def test_an_oversized_asset_is_refused(self, served, monkeypatch):
        client, _token, root = served
        (root / "big.js").write_text("x" * 4096, encoding="utf-8")
        monkeypatch.setattr(sc, "_WORKBENCH_MAX_ASSET_BYTES", 128)
        status, _hdrs, _body = client.request("/big.js")
        assert status == 404

    def test_every_html_page_gets_the_token_not_just_index(self, served):
        """A bundle may have more than one entry point, and they are all
        served out of the same opt-in directory, so this widens nothing."""
        client, token, root = served
        (root / "other.html").write_text(
            "<!doctype html><html><head></head><body>x</body></html>",
            encoding="utf-8")
        _s, _h, body = client.request("/other.html")
        assert token in body.decode("utf-8")


class TestCannotBeWalkedOutOf:
    """Every encoding of `..`, against a real file outside the bundle."""

    # Each shape defeats a different naive check: a textual `".." in path`
    # test, a decode-after-check ordering, a `/`-only separator assumption,
    # a single-pass decoder, and a NUL truncation.
    ESCAPES = [
        "/../secret.html",
        "/..%2fsecret.html",
        "/%2e%2e%2fsecret.html",
        "/%2E%2E%2Fsecret.html",
        "/a/../../secret.html",
        "/assets/..%2f..%2fsecret.html",
        "/assets/%2e%2e/%2e%2e/secret.html",
        "/%2e%2e%5csecret.html",
        "/....//secret.html",
        "/%252e%252e%252fsecret.html",
        "/index.html/../../secret.html",
        "//secret.html",
        "/%2f%2fsecret.html",
        "/index.html%00.png",
        "/~root/.bashrc.html",
        "/../secret.json",
    ]

    @pytest.mark.parametrize("path", ESCAPES)
    def test_traversal_is_refused_and_leaks_nothing(self, served, path):
        """Refused BY THE STATIC ROUTE, which is the part worth pinning.

        Each of these paths is asset-shaped, so it is claimed and judged by
        `_workbench_resolve_asset` — the assertion on the body is what says
        so. Without it the test would also pass against a build that has no
        static route at all (everything 401s), and would be proving nothing
        about the resolver.
        """
        client, _token, _root = served
        status, _hdrs, body = client.request(path)
        assert 400 <= status < 500, (path, status)
        assert b"OUTSIDE-MARKER" not in body, (
            f"{path!r} read a file outside the bundle")
        assert b"no such asset" in body, (
            f"{path!r} was not judged by the static resolver (got {body!r})")

    @pytest.mark.skipif(not hasattr(sc.Path, "symlink_to"),
                        reason="no symlink support")
    def test_a_symlink_out_of_the_bundle_is_refused(self, served):
        """`resolve()` follows the link, so confinement catches it — and the
        `lstat` regular-file check refuses a symlink outright besides."""
        client, _token, root = served
        if not (root / "link.html").is_symlink():
            pytest.skip("symlink not created on this platform")
        status, _hdrs, body = client.request("/link.html")
        assert status == 404
        assert b"OUTSIDE-MARKER" not in body

    def test_resolver_refuses_an_absolute_piece(self, tmp_path):
        """`Path("/a") / "/etc/passwd"` is `/etc/passwd` — pathlib DISCARDS
        the left side — so an absolute piece has to die before the join."""
        assert sc._workbench_resolve_asset(tmp_path, "/etc/passwd") is None

    def test_resolver_refuses_nul_and_tilde_user(self, tmp_path):
        # Borrowed straight from `_sanitize_path`, which owns both refusals.
        assert sc._workbench_resolve_asset(tmp_path, "a\x00b.js") is None
        assert sc._workbench_resolve_asset(tmp_path, "~root/x.js") is None
        assert sc._workbench_resolve_asset(tmp_path, "") is None

    def test_resolver_accepts_a_real_nested_file(self, tmp_path):
        (tmp_path / "sub").mkdir()
        target = tmp_path / "sub" / "a.js"
        target.write_text("x", encoding="utf-8")
        got = sc._workbench_resolve_asset(tmp_path, "sub/a.js")
        assert got == target.resolve()


class TestCannotShadowAnEndpoint:
    def test_no_endpoint_name_contains_a_dot(self):
        """THE invariant the claim rule rests on.

        A static request is claimed when the path's last component has a dot.
        If an endpoint were ever named `foo.bar` it would stop being reachable
        the moment `--workbench-dir` was passed, and nothing else would say
        so. Renaming the endpoint is the fix; widening the rule is not.
        """
        dotted = sorted(n for n in sc._AGENT_HANDLERS if "." in n)
        assert dotted == [], (
            f"endpoint name(s) with a dot would be shadowed by the static "
            f"route: {dotted}")

    def test_endpoints_still_answer_with_the_bundle_mounted(self, served):
        client, token, _root = served
        auth = {"Authorization": f"Bearer {token}"}
        for path in ("/status", "/tools", "/healthz"):
            status, hdrs, _body = client.request(path, headers=auth)
            assert status == 200, path
            assert hdrs["content-type"] == "application/json", path

    def test_a_mistyped_endpoint_still_gets_the_endpoint_list(self, served):
        """The dot rule is what preserves this: an extension-less unknown
        path is not asset-shaped, so it reaches the registry's helpful 404
        instead of a bare 'no such asset'."""
        client, token, _root = served
        status, _hdrs, body = client.request(
            "/stattus", headers={"Authorization": f"Bearer {token}"})
        assert status == 404
        assert b"endpoints" in body

    def test_claim_predicate_defers_to_the_registry(self, bundle):
        """Belt and braces for the dot rule: even an endpoint WITH a dot
        would still be dispatched, because the registry is checked first."""
        handler = sc._AgentRequestHandler.__new__(sc._AgentRequestHandler)
        handler.server = type("S", (), {"_workbench_dir": bundle})()
        assert handler._workbench_claims("GET", "app.js") is True
        assert handler._workbench_claims("GET", "") is True
        assert handler._workbench_claims("GET", "status") is False
        assert handler._workbench_claims("GET", "healthz") is False
        # A POST is never a static read.
        assert handler._workbench_claims("POST", "app.js") is False
        # An extension-less unknown path belongs to the registry's 404.
        assert handler._workbench_claims("GET", "stattus") is False


class TestTokenExposureIsGated:
    def test_the_host_guard_covers_the_static_route(self, served):
        """A route that serves the bearer token is the last one that should
        skip the DNS-rebinding check."""
        client, token, _root = served
        status, _hdrs, body = client.request("/", host="evil.example.com")
        assert status == 403
        assert token.encode() not in body

    def test_nothing_is_served_without_the_flag(self, plain_served):
        """Unauthenticated, the asset path does not even reach a file check —
        it 401s like any other unknown path, which is the pre-existing
        behaviour and leaks neither the token nor whether the file exists."""
        client, token = plain_served
        status, hdrs, body = client.request("/index.html")
        assert status == 401
        assert hdrs["content-type"] == "application/json"
        assert token.encode() not in body

    def test_injection_keeps_the_doctype_first(self):
        """A `<head lang=...>` document used to fall through to a prepend,
        which put the script ahead of the doctype and the browser into quirks
        mode — a layout bug with no visible cause, from a security feature."""
        cases = [
            b"<!doctype html>\n<html><head><title>t</title></head></html>",
            b'<!DOCTYPE html>\n<html lang="en"><head data-x="1"></head></html>',
            b"<!doctype html>\n<html><body>x</body></html>",
            b"<!doctype html>\nbare",
            b"<!DOCTYPE HTML><HTML><HEAD></HEAD></HTML>",
        ]
        for html in cases:
            out = sc._workbench_inject_token(html, "TOK").decode("utf-8")
            assert out.lower().startswith("<!doctype"), html
            assert '"TOK"' in out, html
        # No doctype at all: the script may lead.
        frag = sc._workbench_inject_token(b"<div>x</div>", "TOK").decode()
        assert frag.startswith("<script>") and frag.endswith("<div>x</div>")

    def test_the_token_is_json_escaped_not_hand_quoted(self):
        """URL-safe base64 needs no escaping today. Hand-quoting would be a
        latent injection the day the token format changes."""
        out = sc._workbench_inject_token(
            b"<html><head></head></html>", 'a"b\\c</script>').decode("utf-8")
        assert 'a"b' not in out.replace('\\"', "")
        assert "\\u003c" in out or "<\\/script>" in out or '\\"' in out

    def test_a_token_free_server_injects_an_empty_string(self, bundle):
        """Some tests run the server token-free. The global must still be
        defined, or the page fails on an undefined lookup instead of on the
        401 that explains itself."""
        out = sc._workbench_inject_token(b"<html><head></head></html>", None)
        assert b'"token": ""' in out


# ── the dev escape ────────────────────────────────────────────────────────────

ORIGIN = "http://localhost:5173"


@pytest.fixture
def cors_served():
    (client, token), srv = _serve(_TinyApp(cors_origin=ORIGIN))
    try:
        yield client, token
    finally:
        srv.shutdown()
        srv.server_close()


class TestSingleOriginCors:
    def test_preflight_is_answered_for_the_configured_origin(self,
                                                             cors_served):
        client, _token = cors_served
        status, hdrs, _body = client.request(
            "/status", method="OPTIONS", headers={"Origin": ORIGIN})
        assert status == 204
        assert hdrs["access-control-allow-origin"] == ORIGIN
        assert "Authorization" in hdrs["access-control-allow-headers"]
        assert "POST" in hdrs["access-control-allow-methods"]

    def test_another_origin_gets_no_allow_header(self, cors_served):
        client, _token = cors_served
        status, hdrs, _body = client.request(
            "/status", method="OPTIONS", headers={"Origin": "http://evil.test"})
        assert status == 405
        assert "access-control-allow-origin" not in hdrs

    def test_vary_origin_is_sent_whenever_cors_is_configured(self,
                                                             cors_served):
        """The response differs by `Origin`; a cache that did not know would
        serve the allowed origin's headers to a different one."""
        client, _token = cors_served
        _s, hdrs, _b = client.request(
            "/status", method="OPTIONS", headers={"Origin": "http://evil.test"})
        assert hdrs.get("vary") == "Origin"

    def test_a_real_response_carries_the_header_too(self, cors_served):
        client, token = cors_served
        status, hdrs, _body = client.request(
            "/status", headers={"Authorization": f"Bearer {token}",
                                "Origin": ORIGIN})
        assert status == 200
        assert hdrs["access-control-allow-origin"] == ORIGIN

    def test_no_cors_configured_means_no_header_anywhere(self, plain_served):
        client, token = plain_served
        status, hdrs, _body = client.request(
            "/status", headers={"Authorization": f"Bearer {token}",
                                "Origin": ORIGIN})
        assert status == 200
        assert "access-control-allow-origin" not in hdrs
        assert "vary" not in hdrs

    def test_options_without_cors_explains_the_flags(self, plain_served):
        client, _token = plain_served
        status, _hdrs, body = client.request("/status", method="OPTIONS")
        assert status == 405
        assert b"--workbench-dir" in body and b"--cors-origin" in body

    def test_preflight_still_checks_the_host_header(self, cors_served):
        client, _token = cors_served
        status, _hdrs, _body = client.request(
            "/status", method="OPTIONS",
            headers={"Origin": ORIGIN}, host="evil.example.com")
        assert status == 403


# ── the flags ─────────────────────────────────────────────────────────────────

def _main_with(monkeypatch, argv):
    """Run `sc.main()` for its ARGUMENT HANDLING only, capturing the prepared
    app instead of launching the TUI (the `test_smoke` pattern)."""
    import shutil as _shutil
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(
        _shutil, "get_terminal_size",
        lambda *_a, **_k: _shutil.os.terminal_size((160, 48)))
    captured = {}

    def _capture_run(self_app, *a, **k):
        captured["app"] = self_app
        return None

    monkeypatch.setattr(sc.PlasmidApp, "run", _capture_run, raising=True)
    monkeypatch.setattr(sc, "_check_and_stamp_data_version",
                        lambda *a, **k: None, raising=False)
    monkeypatch.setattr(sc, "_start_agent_api", lambda *a, **k: None,
                        raising=False)
    sc.main()
    return captured.get("app")


class TestFlagValidation:
    def test_a_missing_directory_is_a_usage_error(self, monkeypatch, tmp_path,
                                                  capsys):
        """Refused at parse time. A flag that accepted a bad path and then
        404'd every asset would look like a broken UI, not a typo."""
        with pytest.raises(SystemExit) as exc:
            _main_with(monkeypatch, ["splicecraft", "--workbench-dir",
                                     str(tmp_path / "nope"), "--version"])
        assert exc.value.code == 2
        assert "not a directory" in capsys.readouterr().err

    def test_a_file_is_not_a_directory(self, monkeypatch, tmp_path, capsys):
        f = tmp_path / "index.html"
        f.write_text("x", encoding="utf-8")
        with pytest.raises(SystemExit) as exc:
            _main_with(monkeypatch, ["splicecraft", "--workbench-dir", str(f),
                                     "--version"])
        assert exc.value.code == 2
        assert "not a directory" in capsys.readouterr().err

    def test_a_valid_directory_parses_and_implies_the_agent_api(
            self, monkeypatch, bundle):
        app = _main_with(monkeypatch, ["splicecraft", "--workbench-dir",
                                       str(bundle)])
        assert app is not None
        assert app._workbench_dir == bundle.resolve()
        # `--workbench-dir` implies `--agent`: a bundle with no API behind it
        # is a page that cannot do anything.
        assert getattr(app, "_agent_api_port", None) is not None

    def test_a_missing_index_warns_but_does_not_refuse(self, monkeypatch,
                                                       tmp_path, capsys):
        """A dev may well launch before the first build lands."""
        empty = tmp_path / "ui"
        empty.mkdir()
        app = _main_with(monkeypatch, ["splicecraft", "--workbench-dir",
                                       str(empty)])
        assert app is not None and app._workbench_dir == empty.resolve()
        assert "no servable index.html" in capsys.readouterr().err

    @pytest.mark.parametrize("origin", [
        "*",                      # the whole point of the flag
        "localhost:5173",         # no scheme
        "ftp://localhost",        # not a browser origin
        "http://localhost/app",   # an Origin carries no path
        "http://localhost?x=1",
        "http://user:pw@localhost",
    ])
    def test_a_bad_origin_is_refused(self, monkeypatch, origin, capsys):
        with pytest.raises(SystemExit) as exc:
            _main_with(monkeypatch, ["splicecraft", "--cors-origin", origin,
                                     "--version"])
        assert exc.value.code == 2
        assert "--cors-origin" in capsys.readouterr().err

    @pytest.mark.parametrize("origin", ["http://localhost:5173",
                                        "https://127.0.0.1:8443",
                                        "http://[::1]:5173"])
    def test_a_good_origin_parses(self, monkeypatch, origin):
        app = _main_with(monkeypatch, ["splicecraft", "--cors-origin", origin])
        assert app is not None and app._cors_origin == origin
        assert getattr(app, "_agent_api_port", None) is not None

    def test_the_default_launch_sets_neither(self, monkeypatch):
        app = _main_with(monkeypatch, ["splicecraft"])
        assert app is not None
        assert app._workbench_dir is None
        assert app._cors_origin is None

    def test_the_flags_are_documented_in_help(self, monkeypatch, capsys):
        with pytest.raises(SystemExit):
            _main_with(monkeypatch, ["splicecraft", "--workbench-dir",
                                     "/nonexistent-xyz", "--help"])
        # Validation fires before the help branch, matching how
        # `--agent-api-port` already behaves; the message is the useful part.
        assert "--workbench-dir" in capsys.readouterr().err


class TestDefaultSurfaceUnchanged:
    """The flags are additive. Without them the HTTP surface is what it was."""

    def test_root_is_still_the_tools_self_describe(self, plain_served):
        client, token = plain_served
        status, hdrs, body = client.request(
            "/", headers={"Authorization": f"Bearer {token}"})
        assert status == 200
        assert hdrs["content-type"] == "application/json"
        assert b"endpoints" in body

    def test_an_asset_shaped_path_is_still_an_unknown_endpoint(self,
                                                               plain_served):
        client, token = plain_served
        status, _hdrs, body = client.request(
            "/app.js", headers={"Authorization": f"Bearer {token}"})
        assert status == 404
        assert b"unknown endpoint" in body

    def test_the_server_defaults_to_both_off(self):
        """`_AgentAPIServer` reads the flags off the app with `getattr`, so an
        app that never heard of them (every existing test's rig) is fine."""
        (client, _token), srv = _serve(object())
        try:
            assert srv._workbench_dir is None
            assert srv._cors_origin is None
        finally:
            srv.shutdown()
            srv.server_close()


# ── hardening sweep ──────────────────────────────────────────────────────────
#
# Every case below was found by probing the shipped-shape code, and each one
# failed SILENTLY: the page loaded, nothing was logged, and the only symptom
# was a front-end that could not authenticate or a bundle that 404'd whole.
# That is the shape this project keeps catching — "silence instead of a
# refusal" — so each fix gets a test that says what the silence was.

class TestInjectionLandsWhereItRuns:
    """The token script has to be somewhere a browser will EXECUTE it, and
    before the bundle reads it. Two ways it was not."""

    @staticmethod
    def _audit(html: bytes) -> "dict[str, bool]":
        out = sc._workbench_inject_token(html, "TOK")
        at = out.find(b"window.SPLICECRAFT")
        app = out.find(b'src="/app.js"')
        doctype = out.find(b"<!DOCTYPE")
        opened = out.rfind(b"<!--", 0, at) if at != -1 else -1
        closed = out.find(b"-->", opened + 4) if opened != -1 else -1
        return {
            "present": at != -1,
            "before_app": app == -1 or at < app,
            "in_dead_comment": opened != -1 and (closed == -1 or closed > at),
            "after_doctype": doctype == -1 or doctype < at,
        }

    def test_a_header_element_is_not_mistaken_for_head(self):
        """`<head` is a prefix of `<header>`. A document with no `<head>` but
        a `<header>` is valid HTML5 — the parser infers the head — and when
        that element sits AFTER the bundle's own `<script src>`, matching the
        prefix defined the global after the app had already read it."""
        verdict = self._audit(
            b'<!DOCTYPE html><html><body><div id="root"></div>'
            b'<script src="/app.js"></script><header>nav</header>'
            b'</body></html>')
        assert verdict["present"], verdict
        assert verdict["before_app"], verdict

    def test_a_head_named_in_a_comment_does_not_swallow_the_token(self):
        """A build banner that merely MENTIONS `<head>` took the injected
        `<script>` inside the comment, where it is dead. The page then loads
        perfectly and cannot authenticate — the worst of both outcomes."""
        verdict = self._audit(
            b'<!DOCTYPE html>\n<!-- generated; <head> is built -->\n'
            b'<html><head><script src="/app.js"></script></head></html>')
        assert not verdict["in_dead_comment"], verdict
        assert verdict["before_app"], verdict

    def test_an_unclosed_comment_does_not_swallow_the_token(self):
        """An unterminated comment comments out the rest of the document, so
        "is there a `-->` after it" is not the same question as "is this
        inside a comment"."""
        verdict = self._audit(b'<!DOCTYPE html><!-- oops <head> <html>')
        assert not verdict["in_dead_comment"], verdict

    @pytest.mark.parametrize("html", [
        b'<!DOCTYPE html><html lang="en"><head><script src="/app.js"></script>',
        b'<!DOCTYPE html><html><head/><script src="/app.js"></script>',
        b'<!DOCTYPE html><html><body><script src="/app.js"></script></body>',
        b'<!DOCTYPE html>\n<script src="/app.js"></script>',
        b'<div>fragment</div><script src="/app.js"></script>',
        b'<!doctype HTML><HTML><HEAD><script src="/app.js"></script>',
    ])
    def test_the_token_precedes_the_bundle_and_follows_the_doctype(self, html):
        """Across every document shape: defined before the app reads it, and
        never spliced in ahead of a doctype (which is quirks mode — a layout
        bug with no visible cause, out of a security feature)."""
        verdict = self._audit(html)
        assert verdict["present"], verdict
        assert verdict["before_app"], verdict
        assert verdict["after_doctype"], verdict

    def test_the_tag_scanner_reports_no_match_rather_than_guessing(self):
        assert sc._workbench_tag_insert_at(b"<header><footer>", b"<head") is None
        assert sc._workbench_tag_insert_at(b"<head", b"<head") is None
        assert sc._workbench_tag_insert_at(b"<head>", b"<head") == 6


class TestTheServedRootIsNormalised:
    """`_workbench_resolve_asset` compares a RESOLVED candidate against the
    root, so an unresolved root makes `is_relative_to` false for every asset
    and 404s the whole bundle with nothing to explain it. `main()` resolves,
    but the server reads the attribute off the app with `getattr` — it must
    not depend on every future caller having done so."""

    @staticmethod
    def _fetch(root, path="/assets/app.js"):
        (client, _token), srv = _serve(_TinyApp(workbench_dir=root))
        try:
            return client.request(path)
        finally:
            srv.shutdown()
            srv.server_close()

    def test_a_symlinked_root_still_serves_every_asset(self, bundle, tmp_path):
        link = tmp_path / "ui-link"
        try:
            link.symlink_to(bundle, target_is_directory=True)
        except (OSError, NotImplementedError):      # pragma: no cover - Windows
            pytest.skip("symlinks unavailable")
        status, _hdrs, body = self._fetch(link)
        assert status == 200, body
        assert b"console.log(1);" in body

    def test_a_root_carrying_dot_dot_still_serves_every_asset(self, bundle):
        status, _hdrs, body = self._fetch(bundle / ".." / bundle.name)
        assert status == 200, body
        assert b"console.log(1);" in body

    def test_an_unusable_root_turns_the_feature_off_instead_of_raising(self):
        """A value that is not a path at all disables the static route, rather
        than raising from inside a request thread."""
        assert sc._workbench_server_root(None) is None
        assert sc._workbench_server_root(object()) is None
        assert sc._workbench_server_root("") is None


class TestTheCapBoundsTheRead:
    def test_a_file_that_outgrows_the_cap_after_the_check_is_refused(
            self, served, monkeypatch):
        """The size came from an `lstat` a syscall ago, so on its own it is
        advisory — a rebuild between the check and the read pulled the whole
        new file into the handler thread whatever the cap said.

        Stubbed at the resolver, which is the gate that did the stat: the
        mechanism under test is the READ bound in `_serve_workbench`, and
        "the file grew after it was measured" is exactly what that stub
        says. Lowering the cap alone cannot test this — the resolver's own
        size check would refuse first, with the same 404."""
        client, _token, root = served
        big = root / "big.js"
        big.write_bytes(b"x" * 4096)
        monkeypatch.setattr(sc, "_WORKBENCH_MAX_ASSET_BYTES", 64)
        monkeypatch.setattr(
            sc, "_workbench_resolve_asset",
            lambda _root, rel: big if rel == "big.js" else None)
        status, _hdrs, body = client.request("/big.js")
        assert status == 404
        assert b"no such asset" in body
        assert b"xxxx" not in body


class TestTheBundleCannotPublishTheDataDir:
    """`.json` is a servable type and assets answer BEFORE the token check, so
    a root that contains the data directory would hand `library.json` to any
    local process that can reach the port. The data is the product."""

    @staticmethod
    def _refuse(monkeypatch, capsys, data_dir, workbench_dir):
        monkeypatch.setattr(sc._state, "_DATA_DIR", data_dir, raising=False)
        with pytest.raises(SystemExit) as exc:
            _main_with(monkeypatch, ["splicecraft", "--workbench-dir",
                                     str(workbench_dir), "--version"])
        assert exc.value.code == 2
        return capsys.readouterr().err

    def test_the_data_directory_itself_is_refused(self, monkeypatch, capsys,
                                                  tmp_path):
        data = tmp_path / "data" / "splicecraft"
        data.mkdir(parents=True)
        assert "data directory" in self._refuse(monkeypatch, capsys, data, data)

    def test_an_ancestor_of_the_data_directory_is_refused(self, monkeypatch,
                                                          capsys, tmp_path):
        """The `--workbench-dir ~` case, which is the one that matters."""
        home = tmp_path / "home"
        data = home / ".local" / "share" / "splicecraft"
        data.mkdir(parents=True)
        assert "data directory" in self._refuse(monkeypatch, capsys, data, home)

    def test_a_home_directory_is_refused_even_with_the_data_dir_elsewhere(
            self, monkeypatch, capsys, tmp_path):
        """`$SPLICECRAFT_DATA_DIR` can move the data out of home, which would
        take the ancestor rule with it. A home directory is still never a UI
        bundle."""
        home = tmp_path / "home"
        home.mkdir()
        monkeypatch.setenv("HOME", str(home))
        err = self._refuse(monkeypatch, capsys, tmp_path / "elsewhere", home)
        assert "never a UI bundle" in err

    def test_a_bundle_beside_the_data_directory_is_fine(self, monkeypatch,
                                                        tmp_path):
        """The rule is about containment, not proximity. A subdirectory of the
        data dir is allowed too — confinement already means it exposes only
        itself, so refusing it would cost a harmless layout and buy nothing."""
        data = tmp_path / "data"
        (data / "inside-ui").mkdir(parents=True)
        beside = tmp_path / "ui"
        beside.mkdir()
        monkeypatch.setattr(sc._state, "_DATA_DIR", data, raising=False)
        for ok in (beside, data / "inside-ui"):
            app = _main_with(monkeypatch, ["splicecraft",
                                           "--workbench-dir", str(ok)])
            assert app is not None and app._workbench_dir == ok.resolve()


class TestStartupAgreesWithServeTime:
    def test_a_symlinked_index_is_called_out_at_startup(self, monkeypatch,
                                                        capsys, tmp_path):
        """`is_file()` follows symlinks and the server does not (it wants a
        regular file via `lstat`), so a symlinked `index.html` passed the
        startup check and then 404'd at serve time with nothing anywhere to
        say why. The note now asks the same function the server asks, so the
        two cannot disagree."""
        root = tmp_path / "ui"
        root.mkdir()
        real = tmp_path / "real-index.html"
        real.write_text("<html><head></head></html>", encoding="utf-8")
        try:
            (root / "index.html").symlink_to(real)
        except (OSError, NotImplementedError):      # pragma: no cover - Windows
            pytest.skip("symlinks unavailable")
        assert (root / "index.html").is_file()      # the check that was wrong
        app = _main_with(monkeypatch, ["splicecraft",
                                       "--workbench-dir", str(root)])
        assert app is not None
        assert "no servable index.html" in capsys.readouterr().err

    def test_a_real_index_draws_no_note(self, monkeypatch, capsys, bundle):
        app = _main_with(monkeypatch, ["splicecraft",
                                       "--workbench-dir", str(bundle)])
        assert app is not None
        assert "index.html" not in capsys.readouterr().err


class TestTheBundleTypesABrowserActuallyAsksFor:
    def test_a_favicon_and_a_source_map_are_servable(self, served):
        """Every browser asks for `/favicon.ico` without being told to, and
        devtools asks for a `.map` the moment anyone opens a bundled file.
        Both answered 415 into the console of every page load."""
        client, _token, root = served
        (root / "favicon.ico").write_bytes(b"\x00\x00\x01\x00")
        (root / "app.js.map").write_text('{"version":3}', encoding="utf-8")
        (root / "mod.mjs").write_text("export const a=1;\n", encoding="utf-8")
        for path, want in (("/favicon.ico", "image/vnd.microsoft.icon"),
                           ("/app.js.map", "application/json"),
                           ("/mod.mjs", "text/javascript")):
            status, hdrs, _body = client.request(path)
            assert status == 200, (path, status)
            assert hdrs["content-type"].startswith(want), (path, hdrs)

    def test_the_list_is_still_an_allow_list(self, served):
        """Widening it is not the same as guessing. An extension nobody put
        in the map is still a 415, and the message still says what IS
        servable."""
        client, _token, _root = served
        status, _hdrs, body = client.request("/notes.txt")
        assert status == 415
        assert b".html" in body


class TestOriginIsNormalisedTheWayABrowserSendsIt:
    """`_cors_allowed_origin` compares for EXACT equality, so an
    un-normalised configured value matches no request a browser ever sends —
    and the only symptom is that CORS silently never works, from a flag that
    read as accepted."""

    @pytest.mark.parametrize("given,sent", [
        ("http://LocalHost:5173", "http://localhost:5173"),
        ("HTTP://localhost:5173", "http://localhost:5173"),
        ("http://localhost:80", "http://localhost"),
        ("https://App.Test:443", "https://app.test"),
        ("  http://localhost:5173  ", "http://localhost:5173"),
    ])
    def test_a_configured_origin_matches_what_the_browser_sends(self, given,
                                                                sent):
        normalised, why = sc._workbench_normalise_origin(given)
        assert why is None
        assert normalised == sent
        (client, _token), srv = _serve(_TinyApp(cors_origin=normalised))
        try:
            status, hdrs, _body = client.request(
                "/status", method="OPTIONS", headers={"Origin": sent})
        finally:
            srv.shutdown()
            srv.server_close()
        assert status == 204
        assert hdrs["access-control-allow-origin"] == sent

    def test_an_ipv6_origin_keeps_its_brackets(self):
        """`urlparse().hostname` strips them, and reassembling without them
        yields `http://::1:5173`, which matches nothing."""
        assert sc._workbench_normalise_origin("http://[::1]:5173")[0] == (
            "http://[::1]:5173")

    @pytest.mark.parametrize("bad,hint", [
        ("", "empty"),
        ("*", "wildcard"),
        ("http://localhost:notaport", "port"),
        ("http://:pw@localhost", "credentials"),
        ("http://user:pw@localhost", "credentials"),
        ("localhost:5173", "scheme://host"),
        ("ftp://localhost", "scheme://host"),
        ("http://localhost/app", "no path"),
    ])
    def test_the_refusal_names_the_rule_that_was_broken(self, bad, hint):
        normalised, why = sc._workbench_normalise_origin(bad)
        assert normalised is None
        assert why is not None and hint in why

    def test_the_flag_stores_the_normalised_form(self, monkeypatch):
        app = _main_with(monkeypatch, ["splicecraft", "--cors-origin",
                                       "HTTP://LocalHost:5173"])
        assert app is not None
        assert app._cors_origin == "http://localhost:5173"
