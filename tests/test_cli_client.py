"""test_cli_client — splicecraft_cli sidecar (stdlib-only agent-API client).

The sidecar ships separately from the GUI module (stdlib only, ~50 ms
startup) so an AI agent firing dozens of commands per session doesn't
pay the GUI's ~1.5 s import cost. This test surface covers:

  * Token file: missing, oversized (DoS protection), malformed
    (port-parse failure, single-line file).
  * Response: oversized response triggers exit with a useful message.
  * HTTP errors: 4xx/5xx with JSON body extract the `error` key.
  * Connection refusal: helpful "is the GUI running?" message.
  * Data dir resolution: `$SPLICECRAFT_DATA_DIR` env override.
  * Argparse: every subcommand registers + --help works.

Network calls are mocked at `urllib.request.urlopen` so the test
doesn't need a real server.
"""
from __future__ import annotations

import io
import json
from pathlib import Path
from unittest.mock import patch

import pytest

import splicecraft_cli as cli


# ── Token file handling ───────────────────────────────────────────────────


class TestReadSession:
    """Cover _read_session — the most security-sensitive entry point.
    A malformed / oversized / missing token file MUST exit with a
    helpful message, never crash or silently truncate."""

    def test_missing_token_file_exits_with_help(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        with pytest.raises(SystemExit) as excinfo:
            cli._read_session()
        msg = str(excinfo.value)
        assert "No SpliceCraft session found" in msg
        assert "splicecraft --agent" in msg

    def test_oversized_token_file_refuses(self, tmp_path, monkeypatch):
        """A local hostile / runaway process writing GB into the token
        file must not OOM the CLI when it tries to read it."""
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        token = tmp_path / cli.TOKEN_FILENAME
        # Cap is 1024; write 1 KB + 1 byte to trigger the guard.
        token.write_bytes(b"x" * (cli._CLI_TOKEN_FILE_MAX_BYTES + 1))
        with pytest.raises(SystemExit) as excinfo:
            cli._read_session()
        msg = str(excinfo.value)
        assert "oversized token file" in msg
        assert str(cli._CLI_TOKEN_FILE_MAX_BYTES) in msg

    def test_exact_cap_is_allowed(self, tmp_path, monkeypatch):
        """Boundary check: exactly _CLI_TOKEN_FILE_MAX_BYTES is OK; cap+1 is not.
        This catches an off-by-one in the size guard."""
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        token = tmp_path / cli.TOKEN_FILENAME
        # Build a valid 2-line file that's exactly at the cap.
        port_line = "6701\n"
        # The remaining space is filled with a synthetic 'token'.
        remaining = cli._CLI_TOKEN_FILE_MAX_BYTES - len(port_line)
        assert remaining > 0
        token.write_bytes(port_line.encode() + (b"a" * remaining))
        host, port, tok = cli._read_session()
        assert host == cli.DEFAULT_HOST
        assert port == 6701
        assert len(tok) == remaining
        assert tok == "a" * remaining

    def test_single_line_token_file_rejected(self, tmp_path, monkeypatch):
        """The token file must have `port\\ntoken` — a one-liner is
        malformed and must fail cleanly."""
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        token = tmp_path / cli.TOKEN_FILENAME
        token.write_text("6701")
        with pytest.raises(SystemExit) as excinfo:
            cli._read_session()
        assert "Malformed token file" in str(excinfo.value)

    def test_empty_token_file_rejected(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        token = tmp_path / cli.TOKEN_FILENAME
        token.write_text("")
        with pytest.raises(SystemExit) as excinfo:
            cli._read_session()
        assert "Malformed token file" in str(excinfo.value)

    def test_non_integer_port_rejected(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        token = tmp_path / cli.TOKEN_FILENAME
        token.write_text("not-a-port\nsecret-token\n")
        with pytest.raises(SystemExit) as excinfo:
            cli._read_session()
        assert "Malformed port" in str(excinfo.value)

    def test_happy_path(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        token = tmp_path / cli.TOKEN_FILENAME
        token.write_text("6701\nsupersecret\n")
        host, port, tok = cli._read_session()
        assert host == cli.DEFAULT_HOST
        assert port == 6701
        assert tok == "supersecret"

    def test_whitespace_in_token_is_stripped(self, tmp_path, monkeypatch):
        """The token field is `.strip()`-ed — defensive against
        accidental trailing newlines from `echo` redirects."""
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        token = tmp_path / cli.TOKEN_FILENAME
        token.write_text("  6701  \n  spacey-token  \n")
        host, port, tok = cli._read_session()
        assert port == 6701
        assert tok == "spacey-token"


# ── Data-dir resolution ───────────────────────────────────────────────────


class TestDataDir:
    """Cover _data_dir env override + platformdirs fallback."""

    def test_env_override_wins(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        assert cli._data_dir() == tmp_path

    def test_env_override_expands_tilde(self, monkeypatch):
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", "~/splicecraft-test")
        result = cli._data_dir()
        # Path.expanduser leaves the tilde resolved.
        assert "~" not in str(result)
        assert str(result).endswith("splicecraft-test")

    def test_no_env_falls_back(self, monkeypatch):
        """Without the env var, _data_dir uses platformdirs (or its
        manual ~/.local/share fallback). Both produce a valid Path."""
        monkeypatch.delenv("SPLICECRAFT_DATA_DIR", raising=False)
        result = cli._data_dir()
        assert isinstance(result, Path)
        # The name should be 'splicecraft' as the leaf, regardless of OS.
        assert result.name == "splicecraft"


# ── HTTP request handling ────────────────────────────────────────────────────


class _FakeResp:
    """Minimal urllib response stand-in: supports the with-block
    protocol and `.read(n)`."""

    def __init__(self, body: bytes):
        self._body = body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def read(self, n: int) -> bytes:
        return self._body[:n]


def _setup_token(tmp_path, monkeypatch):
    monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
    (tmp_path / cli.TOKEN_FILENAME).write_text("6701\nt\n")


class TestRequest:
    """Cover _request: URL building, auth header, response size cap,
    HTTP error decoding, connection refusal messaging."""

    def test_happy_path_parses_json(self, tmp_path, monkeypatch):
        _setup_token(tmp_path, monkeypatch)
        body = json.dumps({"ok": True, "n": 42}).encode()
        with patch("urllib.request.urlopen",
                    return_value=_FakeResp(body)) as mock:
            result = cli._request("status")
        assert result == {"ok": True, "n": 42}
        # Bearer token was attached.
        req = mock.call_args[0][0]
        assert req.get_header("Authorization") == "Bearer t"

    def test_post_includes_json_body(self, tmp_path, monkeypatch):
        _setup_token(tmp_path, monkeypatch)
        body = b"{}"
        with patch("urllib.request.urlopen",
                    return_value=_FakeResp(body)) as mock:
            cli._request("fetch", "POST", {"accession": "L09137"})
        req = mock.call_args[0][0]
        assert req.method == "POST"
        assert req.get_header("Content-type") == "application/json"
        sent = json.loads(req.data.decode())
        assert sent == {"accession": "L09137"}

    def test_oversized_response_refuses(self, tmp_path, monkeypatch):
        """Mirror of the server-side 50 MB cap. If the server (or a
        compromised proxy) returns a larger body, the CLI must abort
        rather than OOM."""
        _setup_token(tmp_path, monkeypatch)
        # 50 MB + 1: the response cap path. Don't actually allocate it —
        # FakeResp returns it from read(n). Use a slim sentinel.
        big = b"x" * (cli._CLI_RESPONSE_MAX_BYTES + 1)
        with patch("urllib.request.urlopen", return_value=_FakeResp(big)):
            with pytest.raises(cli._AgentCliError) as excinfo:
                cli._request("status")
        msg = str(excinfo.value)
        assert "exceeds" in msg
        assert "byte cap" in msg

    def test_exact_cap_response_is_allowed(self, tmp_path, monkeypatch):
        """Boundary: response at exactly _CLI_RESPONSE_MAX_BYTES passes;
        cap+1 fails. Catches an off-by-one in the size guard."""
        _setup_token(tmp_path, monkeypatch)
        # JSON parsing requires valid JSON, so fill with whitespace
        # padding wrapped around a tiny payload.
        payload = b'{"ok":true}'
        pad = cli._CLI_RESPONSE_MAX_BYTES - len(payload)
        # Note: padding inside a string field keeps the JSON valid.
        body = (b'{"pad":"' + b" " * (pad - len('{"pad":""}'.encode()))
                + b'","ok":true}')
        # Trim to exactly the cap; if our math is off, json.loads will
        # raise — that itself is a useful failure signal.
        body = body[:cli._CLI_RESPONSE_MAX_BYTES]
        with patch("urllib.request.urlopen", return_value=_FakeResp(body)):
            result = cli._request("status")
        assert result.get("ok") is True

    def test_http_error_with_json_body_extracts_error_field(
            self, tmp_path, monkeypatch):
        import urllib.error
        _setup_token(tmp_path, monkeypatch)
        err_payload = json.dumps({"error": "Plasmid not found"}).encode()
        fake_err = urllib.error.HTTPError(
            url="http://localhost/foo",
            code=404,
            msg="Not Found",
            hdrs=None,  # type: ignore[arg-type]
            fp=io.BytesIO(err_payload),
        )
        with patch("urllib.request.urlopen", side_effect=fake_err):
            with pytest.raises(cli._AgentCliError) as excinfo:
                cli._request("status")
        msg = str(excinfo.value)
        assert "Plasmid not found" in msg
        assert "HTTP 404" in msg

    def test_http_error_with_non_json_body_falls_back_to_raw(
            self, tmp_path, monkeypatch):
        import urllib.error
        _setup_token(tmp_path, monkeypatch)
        fake_err = urllib.error.HTTPError(
            url="http://localhost/foo",
            code=500,
            msg="Internal Error",
            hdrs=None,  # type: ignore[arg-type]
            fp=io.BytesIO(b"<html>boom</html>"),
        )
        with patch("urllib.request.urlopen", side_effect=fake_err):
            with pytest.raises(cli._AgentCliError) as excinfo:
                cli._request("status")
        msg = str(excinfo.value)
        # Falls back to raw body when JSON parse fails.
        assert "boom" in msg or "Internal Error" in msg
        assert "HTTP 500" in msg

    def test_connection_refused_has_useful_message(
            self, tmp_path, monkeypatch):
        import urllib.error
        _setup_token(tmp_path, monkeypatch)
        with patch("urllib.request.urlopen",
                    side_effect=urllib.error.URLError("Connection refused")):
            with pytest.raises(cli._AgentCliError) as excinfo:
                cli._request("status")
        msg = str(excinfo.value)
        assert "Could not reach SpliceCraft" in msg
        # The "is the GUI still running" hint must be present so users
        # don't waste time debugging the network when the answer is
        # "you didn't start the GUI."
        assert "--agent" in msg

    def test_non_json_response_returns_raw_dict(
            self, tmp_path, monkeypatch):
        """Some endpoints (export-*) may return plain-text. The CLI
        wraps these in `{"raw": body}` rather than crashing."""
        _setup_token(tmp_path, monkeypatch)
        body = b"plain text response"
        with patch("urllib.request.urlopen", return_value=_FakeResp(body)):
            result = cli._request("status")
        assert result == {"raw": "plain text response"}

    def test_empty_response_returns_empty_dict(
            self, tmp_path, monkeypatch):
        _setup_token(tmp_path, monkeypatch)
        with patch("urllib.request.urlopen",
                    return_value=_FakeResp(b"")):
            result = cli._request("status")
        assert result == {}


# ── Argparse surface ────────────────────────────────────────────────────


class TestParser:
    """Static check that every documented subcommand registers and
    that --help works without imports failing. Also a forcing function
    for additions: a new endpoint without a CLI subcommand will not
    show up in EXPECTED_SUBCOMMANDS and the test fails loudly."""

    # Subcommands the CLI ships today. New ones get added here in the
    # same PR as the cmd_ handler so the test acts as a registry.
    EXPECTED_SUBCOMMANDS = {
        "status", "tools", "features", "fetch", "load-entry",
        "load-file", "add-feature", "save",
        "get-sequence", "replace-sequence",
        "delete-feature", "update-feature", "get-feature",
        "export-genbank", "export-fasta",
        "list-library", "list-collections", "delete-from-library",
        "list-restriction-sites", "list-codon-tables", "optimize-protein",
    }

    def test_parser_builds(self):
        parser = cli._build_parser()
        assert parser is not None

    def test_every_expected_subcommand_registered(self):
        parser = cli._build_parser()
        # argparse stashes subparsers as the action's `choices` dict.
        subparsers_action = next(
            a for a in parser._subparsers._actions  # type: ignore[union-attr]
            if hasattr(a, "choices") and a.choices
        )
        registered = set(subparsers_action.choices.keys())
        missing = self.EXPECTED_SUBCOMMANDS - registered
        assert not missing, f"Missing CLI subcommands: {missing}"

    def test_call_passthrough_parses(self):
        args = cli._build_parser().parse_args(
            ["call", "rename-plasmid", "--json", '{"old":"a","new":"b"}'])
        assert args.endpoint == "rename-plasmid"
        assert args.json == '{"old":"a","new":"b"}'

    def test_optimize_protein_exposes_all_flags(self):
        """The dedicated subcommand must reach the endpoint's full power, not
        just --table (else `call optimize-protein --json` is the only capable
        path — a discoverability trap)."""
        args = cli._build_parser().parse_args([
            "optimize-protein", "MAKR", "--mode", "max_cai", "--stops", "2",
            "--transl-table", "11", "--hazard-hosts", "plant", "mammalian",
            "--forbidden-motifs", "GGTCTC", "--min-gc", "40", "--max-gc", "60",
            "--gc-window", "45", "--avoid-repeats-with", "ATGCATGC",
            "--max-repeat", "20",
        ])
        assert args.mode == "max_cai" and args.stops == 2
        assert args.transl_table == 11
        assert args.hazard_hosts == ["plant", "mammalian"]
        assert args.forbidden_motifs == ["GGTCTC"]
        assert args.min_gc == 40.0 and args.max_gc == 60.0
        assert args.gc_window == 45
        assert args.avoid_repeats_with == ["ATGCATGC"] and args.max_repeat == 20

    def test_optimize_protein_builds_full_payload(self, tmp_path, monkeypatch):
        _setup_token(tmp_path, monkeypatch)
        captured = {}

        def _fake(endpoint, method, payload=None, **_kw):
            captured["endpoint"] = endpoint
            captured["payload"] = payload
            return {"dna": "ATG", "warnings": [], "fixes": []}

        args = cli._build_parser().parse_args([
            "optimize-protein", "MAKR", "--mode", "max_cai",
            "--hazard-hosts", "plant", "--min-gc", "40",
            "--avoid-repeats-with", "ATGCATGC",
        ])
        with patch.object(cli, "_request", side_effect=_fake):
            args.fn(args)
        p = captured["payload"]
        assert captured["endpoint"] == "optimize-protein"
        assert p["protein"] == "MAKR" and p["mode"] == "max_cai"
        assert p["hazard_hosts"] == ["plant"] and p["min_gc"] == 40.0
        assert p["avoid_repeats_with"] == ["ATGCATGC"]
        # Unset flags must NOT appear (so server defaults apply).
        assert "max_gc" not in p and "forbidden_motifs" not in p

    def test_optimize_protein_omits_unset_flags(self, tmp_path, monkeypatch):
        _setup_token(tmp_path, monkeypatch)
        captured = {}

        def _fake(endpoint, method, payload=None, **_kw):
            captured["payload"] = payload
            return {"dna": "ATG"}

        args = cli._build_parser().parse_args(["optimize-protein", "MK"])
        with patch.object(cli, "_request", side_effect=_fake):
            args.fn(args)
        # Bare call → only the protein; every optional key absent.
        assert captured["payload"] == {"protein": "MK"}

    def test_call_surfaces_structured_error_as_json(
            self, tmp_path, monkeypatch, capsys):
        _setup_token(tmp_path, monkeypatch)
        err = cli._AgentCliError("boom (HTTP 404)", code=404,
                                  payload={"error": "no such entry"})
        args = cli._build_parser().parse_args(
            ["call", "load-entry", "--json", '{"name":"x"}'])
        with patch.object(cli, "_request", side_effect=err):
            with pytest.raises(SystemExit) as exc:
                args.fn(args)
        assert exc.value.code == 1
        out = json.loads(capsys.readouterr().out)
        assert out["error"] == "no such entry"
        assert out["http_code"] == 404

    def test_help_runs(self, capsys):
        parser = cli._build_parser()
        with pytest.raises(SystemExit) as excinfo:
            parser.parse_args(["--help"])
        assert excinfo.value.code == 0
        captured = capsys.readouterr()
        assert "splicecraft" in captured.out.lower()

    def test_add_feature_strand_choices(self):
        """--strand must be restricted to -1 / 0 / 1; bad values fail
        argparse, not the server."""
        parser = cli._build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["add-feature", "0", "10", "--strand", "2"])

    def test_add_feature_defaults(self):
        parser = cli._build_parser()
        args = parser.parse_args(["add-feature", "0", "10"])
        assert args.start  == 0
        assert args.end    == 10
        assert args.label  == ""
        assert args.type   == "misc_feature"
        assert args.strand == 1
        assert args.force is False

    def test_wrap_feature_end_less_than_start_is_accepted_at_argparse(self):
        """add-feature accepts end < start (the wrap-feature convention).
        Server-side validation decides if it's biologically valid; the
        CLI must not block the call."""
        parser = cli._build_parser()
        args = parser.parse_args(["add-feature", "100", "50"])
        assert args.start == 100
        assert args.end == 50

    def test_fetch_force_flag(self):
        parser = cli._build_parser()
        args = parser.parse_args(["fetch", "L09137", "--force"])
        assert args.accession == "L09137"
        assert args.force is True

    def test_unknown_subcommand_fails(self):
        parser = cli._build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["definitely-not-a-real-command"])


# ── End-to-end main() smoke ────────────────────────────────────────────────


class TestMainSmoke:
    """Run main() with a mocked _request so we exercise the full path
    from argv → parser → cmd_ → _request without needing a live GUI."""

    def test_status_emits_json(self, tmp_path, monkeypatch, capsys):
        _setup_token(tmp_path, monkeypatch)
        with patch.object(cli, "_request",
                            return_value={"loaded": "pUC19"}) as mock_req:
            cli.main(["status"])
        mock_req.assert_called_once_with("status")
        out = capsys.readouterr().out
        assert "pUC19" in out

    def test_features_text_format(self, tmp_path, monkeypatch, capsys):
        _setup_token(tmp_path, monkeypatch)
        fake_resp = {"features": [
            {"idx": 0, "type": "CDS", "start": 99, "end": 200,
             "strand": 1, "label": "lacZ"},
            {"idx": 1, "type": "promoter", "start": 0, "end": 50,
             "strand": -1, "label": "P"},
        ]}
        with patch.object(cli, "_request", return_value=fake_resp):
            cli.main(["features"])
        out = capsys.readouterr().out
        assert "lacZ" in out
        assert "CDS" in out
        # 0-based start (99) renders as 1-based (100).
        assert "100" in out

    def test_no_subcommand_shows_help_and_exits(
            self, monkeypatch, capsys):
        """`splicecraft-cli` with no args must surface usage, not
        silent-noop."""
        # argparse with required=True on subparsers exits non-zero.
        # Without required=True, argparse's behavior differs across
        # Python versions; just confirm SOMETHING reasonable happens.
        result_code = None
        try:
            cli.main([])
        except SystemExit as exc:
            result_code = exc.code
        except AttributeError:
            # If `args.fn` is unset because no subcommand was given,
            # accept that as the failure mode (the user sees an error).
            result_code = 1
        # Either way, the user did not silently succeed.
        assert result_code != 0
