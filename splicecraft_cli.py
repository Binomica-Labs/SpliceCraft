#!/usr/bin/env python3
"""splicecraft-cli — drive a running SpliceCraft GUI session from any
external CLI agent (Claude Code, Cursor, aider, hand-rolled scripts).

Connects to the localhost JSON API exposed by `splicecraft --agent`
(alias of the original `--agent-api`). Reads connection details
(port + bearer token) from
``~/.local/share/splicecraft/agent_token`` (or
``$SPLICECRAFT_DATA_DIR/agent_token`` when overridden), so the running
GUI is always the destination — no flag-fiddling required.

Stdlib-only by design: imports complete in ~50 ms (vs ~1.5 s for the
GUI module), so an AI agent firing dozens of commands per session
doesn't pay startup cost on every call.

Examples::

    splicecraft-cli status
    splicecraft-cli features
    splicecraft-cli fetch L09137
    splicecraft-cli load-entry pUC19
    splicecraft-cli add-feature 100 200 --label lacZ --type CDS --strand 1
    splicecraft-cli save

Use ``splicecraft-cli tools`` to list every endpoint the running
session exposes (handy when wiring a new tool surface for the agent).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 6701
TOKEN_FILENAME = "agent_token"


def _data_dir() -> Path:
    """Resolve the SpliceCraft user-data directory the same way the
    GUI does — env override first, then the platform default. We keep
    a hand-rolled fallback so the CLI doesn't pull in `platformdirs`
    just for one path lookup (keeps imports stdlib-only)."""
    override = os.environ.get("SPLICECRAFT_DATA_DIR")
    # Sweep #9 (2026-05-19): strip + truthiness check so a stray
    # `SPLICECRAFT_DATA_DIR=" "` (shell scripting bug) doesn't land
    # us on a literal whitespace-named relative dir.
    if override and override.strip():
        return Path(override.strip()).expanduser()
    try:
        from platformdirs import user_data_dir   # noqa: WPS433 (lazy)
        return Path(user_data_dir("splicecraft", appauthor=False,
                                    roaming=False))
    except ImportError:
        return Path.home() / ".local" / "share" / "splicecraft"


def _token_file() -> Path:
    return _data_dir() / TOKEN_FILENAME


_CLI_TOKEN_FILE_MAX_BYTES = 1024


def _read_session() -> tuple[str, int, str]:
    """Return `(host, port, token)` from the running session's token
    file. Exits with a helpful message if no session is up.

    Capped at 1 KB defensively. A local-attacker-with-write or a
    misbehaving co-resident process can fill the token file; reading
    unbounded would let a hostile fill DoS the CLI."""
    f = _token_file()
    if not f.exists():
        sys.exit(
            f"No SpliceCraft session found.\n"
            f"  Expected token file: {f}\n"
            f"  Start the GUI with: splicecraft --agent"
        )
    # Sweep #25 (2026-05-23): `lstat` not `stat` — `f.stat()`
    # follows symlinks. A pre-placed symlink at the token path
    # pointing to (e.g.) `/etc/passwd` would be silently followed,
    # read up to 1 KB, and stderr-leak via the malformed-port
    # error message ("Malformed port in <token-file>: 'root:x:0...'").
    # Niche but trivially fixable.
    import stat as _stat
    try:
        st = f.lstat()
        if _stat.S_ISLNK(st.st_mode):
            sys.exit(
                f"Refusing to read symlinked token file {f}. "
                f"Restart the GUI to regenerate."
            )
        size = st.st_size
    except OSError as exc:
        sys.exit(f"Could not stat token file {f}: {exc}")
    if size > _CLI_TOKEN_FILE_MAX_BYTES:
        sys.exit(
            f"Refusing to read oversized token file {f} "
            f"({size:,} bytes > {_CLI_TOKEN_FILE_MAX_BYTES}-byte cap). "
            f"Restart the GUI to regenerate."
        )
    lines = f.read_text(encoding="utf-8").strip().splitlines()
    if len(lines) < 2:
        sys.exit(
            f"Malformed token file at {f} "
            f"(expected `port\\ntoken`)."
        )
    try:
        port = int(lines[0].strip())
    except ValueError:
        sys.exit(f"Malformed port in {f}: {lines[0]!r}")
    return DEFAULT_HOST, port, lines[1].strip()


_CLI_RESPONSE_MAX_BYTES = 50 * 1024 * 1024


class _AgentCliError(Exception):
    """Raised by `_request` on an unrecoverable call failure (HTTP 4xx /
    5xx, transport error, oversize response). CATCHABLE (snag #6) — a
    4xx no longer hard-kills a batch that imports `_request`. `main()`
    turns an uncaught one into a clean `sys.exit` (preserving the CLI's
    message-and-non-zero-exit UX); `call` catches it to emit the
    structured server error as JSON. Carries the HTTP `code` (None for
    transport failures) and the parsed error `payload` (None if the body
    wasn't a JSON object)."""

    def __init__(self, message, *, code=None, payload=None):
        super().__init__(message)
        self.code = code
        self.payload = payload


def _request(endpoint: str, method: str = "GET",
              payload: "dict | None" = None,
              timeout: float = 30.0) -> dict:
    host, port, token = _read_session()
    url = f"http://{host}:{port}/{endpoint}"
    data = (json.dumps(payload or {}).encode("utf-8")
            if method == "POST" else None)
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {token}")
    if data is not None:
        req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            # Cap the read size symmetric with the server's cap; a
            # compromised / buggy server returning an unbounded stream
            # would otherwise OOM the CLI. Read 1 byte over the cap so
            # we can detect oversize and abort with a useful message
            # rather than silently truncating.
            raw = resp.read(_CLI_RESPONSE_MAX_BYTES + 1)
            if len(raw) > _CLI_RESPONSE_MAX_BYTES:
                raise _AgentCliError(
                    f"Error: response from SpliceCraft exceeds "
                    f"{_CLI_RESPONSE_MAX_BYTES:,}-byte cap "
                    f"(endpoint={endpoint!r}). Refusing to read."
                )
            body = raw.decode("utf-8")
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8") if exc.fp else ""
        err_payload: "dict | None" = None
        try:
            parsed = json.loads(body) if body else {}
            if isinstance(parsed, dict):
                err_payload = parsed
                msg = parsed.get("error", body or exc.reason)
            else:
                msg = body or exc.reason
        except json.JSONDecodeError:
            msg = body or exc.reason
        # Snag #6: raise a CATCHABLE error instead of a hard `sys.exit`,
        # so a 4xx no longer kills an importing batch mid-run. `main()`
        # turns an uncaught one into a clean exit; `call` catches it to
        # surface the structured error as JSON.
        raise _AgentCliError(f"Error: {msg} (HTTP {exc.code})",
                              code=exc.code, payload=err_payload)
    except urllib.error.URLError as exc:
        raise _AgentCliError(
            f"Could not reach SpliceCraft at {host}:{port} ({exc.reason}).\n"
            f"  Is the GUI still running with --agent-api?"
        )
    try:
        return json.loads(body) if body else {}
    except json.JSONDecodeError:
        return {"raw": body}


# ── Subcommand handlers ────────────────────────────────────────────────────────
# Each subcommand maps 1:1 to a `/<endpoint>` on the server. Keep these
# thin — the server is the source of truth for validation / error messages.

def _emit_json(obj) -> None:
    print(json.dumps(obj, indent=2, default=str))


def _call_pick_method(endpoint: str, method: str, payload: "dict | None",
                       *, allow_upgrade: bool) -> dict:
    """Issue the request, retrying once with the method the server names.

    `call` has to guess a method: with no `--json` there is no body to
    infer POST from, so a body-less WRITE endpoint (`add-current-to-library`,
    `save`) guessed GET and came back 405 "requires POST" — leaving the
    caller to re-run by hand. The server's 405 body carries the authoritative
    `allow`, so honour it instead of guessing twice. Reads pay nothing: the
    retry only fires on the 405 that would otherwise have been a hard failure.

    An explicit `--method` is never overridden — the caller asked for that
    method and gets the real 405 back.
    """
    try:
        return _request(endpoint, method=method, payload=payload)
    except _AgentCliError as exc:
        allow = (exc.payload or {}).get("allow") if isinstance(
            exc.payload, dict) else None
        if (allow_upgrade and exc.code == 405 and method == "GET"
                and str(allow or "").upper() == "POST"):
            return _request(endpoint, method="POST", payload=payload or {})
        raise


def cmd_call(args) -> None:
    """Generic passthrough: call ANY agent endpoint by name (snag #1).

    Reuses the same auth / host / port plumbing as the blessed
    subcommands, so an agent integrating SpliceCraft no longer has to
    hand-roll HTTP against the private `_request`. `tools` lists every
    endpoint and its `doc_full` body schema."""
    payload = None
    if args.json:
        try:
            payload = json.loads(args.json)
        except json.JSONDecodeError as exc:
            sys.exit(f"Error: --json is not valid JSON: {exc}")
        if not isinstance(payload, dict):
            sys.exit("Error: --json must be a JSON object")
    # Default to POST when a body is supplied (writes require POST), else
    # GET — matching the server's read/write method split.
    explicit_method = args.method is not None
    method = (args.method or ("POST" if payload is not None else "GET")).upper()
    try:
        result = _call_pick_method(args.endpoint, method, payload,
                                    allow_upgrade=not explicit_method)
    except _AgentCliError as exc:
        # Surface the structured server error as JSON (then exit
        # non-zero) rather than a bare message, so a calling script can
        # parse the failure instead of scraping stderr.
        out: dict = (dict(exc.payload) if isinstance(exc.payload, dict)
                     else {"error": str(exc)})
        out["http_code"] = exc.code
        _emit_json(out)
        raise SystemExit(1)
    _emit_json(result)


def cmd_status(args) -> None:
    _emit_json(_request("status"))


def cmd_tools(args) -> None:
    result = _request("tools")
    if args.json:
        _emit_json(result)
        return
    for ep in result.get("endpoints", []):
        flag = "WRITE" if ep.get("write") else "READ "
        print(f"  {flag}  {ep['name']:14}  {ep.get('doc', '')}")


def cmd_features(args) -> None:
    result = _request("features")
    feats = result.get("features", [])
    if args.json:
        _emit_json(feats)
        return
    if not feats:
        print("(no features)")
        return
    for f in feats:
        strand = ("+" if f["strand"] == 1
                   else "-" if f["strand"] == -1 else ".")
        label = f.get("label") or ""
        print(
            f"  [{f['idx']:3}] {f.get('type','?'):14} "
            f"{f['start'] + 1:>7,}..{f['end']:<7,} {strand}  {label}"
        )


def cmd_fetch(args) -> None:
    payload = {"accession": args.accession}
    if args.force:
        payload["force"] = True
    _emit_json(_request("fetch", "POST", payload, timeout=60))


def cmd_load_entry(args) -> None:
    payload = {"name": args.name}
    if args.force:
        payload["force"] = True
    _emit_json(_request("load-entry", "POST", payload))


def cmd_load_file(args) -> None:
    payload = {"path": args.path}
    if args.force:
        payload["force"] = True
    # Server-side parse on a chromosome can take 10s+; bump timeout.
    _emit_json(_request("load-file", "POST", payload, timeout=120))


def cmd_add_feature(args) -> None:
    payload = {
        "start":  args.start,
        "end":    args.end,
        "label":  args.label,
        "type":   args.type,
        "strand": args.strand,
    }
    if args.force:
        payload["force"] = True
    _emit_json(_request("add-feature", "POST", payload))


def cmd_save(args) -> None:
    _emit_json(_request("save", "POST"))


# ── Tier 1: sequence + feature CRUD ────────────────────────────────────────────


def cmd_get_sequence(args) -> None:
    payload = {"start": args.start, "end": args.end, "bottom": args.bottom}
    result = _request("get-sequence", "POST", payload)
    if args.json:
        _emit_json(result)
        return
    seq = result.get("seq", "")
    print(seq)


def cmd_replace_sequence(args) -> None:
    payload = {"start": args.start, "end": args.end, "bases": args.bases}
    if args.force:
        payload["force"] = True
    _emit_json(_request("replace-sequence", "POST", payload))


def cmd_delete_feature(args) -> None:
    payload = {"idx": args.idx}
    if args.force:
        payload["force"] = True
    _emit_json(_request("delete-feature", "POST", payload))


def cmd_update_feature(args) -> None:
    payload = {"idx": args.idx}
    if args.label  is not None: payload["label"]  = args.label
    if args.type   is not None: payload["type"]   = args.type
    if args.strand is not None: payload["strand"] = args.strand
    if args.force:
        payload["force"] = True
    _emit_json(_request("update-feature", "POST", payload))


def cmd_get_feature(args) -> None:
    _emit_json(_request("get-feature", "POST", {"idx": args.idx}))


def cmd_export_genbank(args) -> None:
    payload = {"path": args.path}
    if args.force:
        payload["force"] = True
    _emit_json(_request("export-genbank", "POST", payload))


def cmd_export_fasta(args) -> None:
    payload = {"path": args.path}
    if args.force:
        payload["force"] = True
    _emit_json(_request("export-fasta", "POST", payload))


# ── Tier 2: library + collections ──────────────────────────────────────────────


def cmd_list_library(args) -> None:
    result = _request("list-library")
    if args.json:
        _emit_json(result)
        return
    entries = result.get("library", [])
    if not entries:
        print("(empty library)")
        return
    for e in entries:
        print(
            f"  {e.get('name','?'):24}  {e.get('length',0):>7,} bp  "
            f"{e.get('n_features',0):>3} feat  {e.get('topology','')}"
        )


def cmd_list_collections(args) -> None:
    result = _request("list-collections")
    if args.json:
        _emit_json(result)
        return
    active = result.get("active") or "(none)"
    print(f"  active: {active}")
    for c in result.get("collections", []):
        marker = "*" if c.get("name") == active else " "
        print(f"  {marker} {c.get('name','?'):24}  "
              f"{c.get('n_plasmids',0):>3} plasmids")


def cmd_delete_from_library(args) -> None:
    payload = {"name": args.name}
    if args.force:
        payload["force"] = True
    _emit_json(_request("delete-from-library", "POST", payload))


# ── Tier 3: cloning / design helpers ───────────────────────────────────────────


def cmd_list_restriction_sites(args) -> None:
    payload: dict = {}
    if args.enzymes:
        payload["enzymes"] = args.enzymes
    if args.min_length is not None:
        payload["min_length"] = args.min_length
    if args.unique_only:
        payload["unique_only"] = True
    result = _request("list-restriction-sites", "POST", payload)
    if args.json:
        _emit_json(result)
        return
    sites = result.get("sites", [])
    if not sites:
        print("(no sites)")
        return
    for s in sites:
        strand = "+" if s.get("strand", 1) == 1 else "-"
        print(f"  {s.get('enzyme','?'):14}  "
              f"{(s.get('start') or 0)+1:>7,}..{s.get('end','?'):<7}  "
              f"{strand}  cut@{s.get('cut_bp', '?')}")


def cmd_list_codon_tables(args) -> None:
    result = _request("list-codon-tables")
    if args.json:
        _emit_json(result)
        return
    for t in result.get("tables", []):
        print(f"  {t.get('taxid',''):>6}  {t.get('source','?'):8}  "
              f"{t.get('name','?')}")


def cmd_optimize_protein(args) -> None:
    payload = {"protein": args.protein}
    if args.table:
        payload["table"] = args.table
    result = _request("optimize-protein", "POST", payload)
    if args.json:
        _emit_json(result)
        return
    print(result.get("dna", ""))


# ── Argparse wiring ────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="splicecraft-cli",
        description=(__doc__ or "").splitlines()[0] if __doc__ else "",
    )
    sub = parser.add_subparsers(dest="cmd", required=True,
                                  metavar="COMMAND")

    p_status = sub.add_parser("status",
                                help="Show what's loaded + dirty flag.")
    p_status.set_defaults(fn=cmd_status)

    p_tools = sub.add_parser("tools",
                                help="List all endpoints the session exposes.")
    p_tools.add_argument("--json", action="store_true",
                          help="Emit JSON instead of a table")
    p_tools.set_defaults(fn=cmd_tools)

    p_feat = sub.add_parser("features",
                              help="List features on the loaded record.")
    p_feat.add_argument("--json", action="store_true",
                         help="Emit JSON instead of a table")
    p_feat.set_defaults(fn=cmd_features)

    p_fetch = sub.add_parser("fetch",
                               help="Fetch a GenBank record from NCBI.")
    p_fetch.add_argument("accession",
                          help="GenBank accession (e.g. L09137).")
    p_fetch.add_argument("--force", action="store_true",
                          help="Override unsaved-changes guard.")
    p_fetch.set_defaults(fn=cmd_fetch)

    p_load = sub.add_parser("load-entry",
                              help="Load a plasmid library entry by name.")
    p_load.add_argument("name", help="Library entry name or id.")
    p_load.add_argument("--force", action="store_true",
                         help="Override unsaved-changes guard.")
    p_load.set_defaults(fn=cmd_load_entry)

    p_loadf = sub.add_parser(
        "load-file",
        help="Load a .gb / .gbk / .dna file from a server-side path "
             "(bypasses the JSON-body size cap; works for chromosome-"
             "scale files).",
    )
    p_loadf.add_argument("path", help="Server-side file path.")
    p_loadf.add_argument("--force", action="store_true",
                          help="Override unsaved-changes guard.")
    p_loadf.set_defaults(fn=cmd_load_file)

    p_add = sub.add_parser(
        "add-feature",
        help="Add a feature to the loaded record.",
    )
    p_add.add_argument("start", type=int, help="0-based start bp.")
    p_add.add_argument("end",   type=int,
                        help="0-based end bp (exclusive). For wrap "
                             "features, pass end < start.")
    p_add.add_argument("--label", default="",
                        help="Feature label (qualifier).")
    p_add.add_argument("--type", default="misc_feature",
                        help="GenBank feature type (CDS, promoter, …).")
    p_add.add_argument("--strand", type=int, default=1,
                        choices=[-1, 0, 1],
                        help="1=forward (default), -1=reverse, 0=both.")
    p_add.add_argument("--force", action="store_true",
                        help="Override unsaved-changes guard.")
    p_add.set_defaults(fn=cmd_add_feature)

    p_save = sub.add_parser("save",
                              help="Save the loaded record (file + library).")
    p_save.set_defaults(fn=cmd_save)

    # ── Tier 1 ────────────────────────────────────────────────────────────

    p_getseq = sub.add_parser(
        "get-sequence",
        help="Extract DNA from a bp range (forward or --bottom strand).",
    )
    p_getseq.add_argument("start", type=int, help="0-based start bp.")
    p_getseq.add_argument("end",   type=int, help="0-based end bp (exclusive).")
    p_getseq.add_argument("--bottom", action="store_true",
                            help="Return reverse-complement (bottom strand 5'→3').")
    p_getseq.add_argument("--json", action="store_true",
                            help="Emit JSON instead of plain seq.")
    p_getseq.set_defaults(fn=cmd_get_sequence)

    p_repseq = sub.add_parser(
        "replace-sequence",
        help="Replace bp range with new bases (mutagenesis).",
    )
    p_repseq.add_argument("start", type=int)
    p_repseq.add_argument("end",   type=int)
    p_repseq.add_argument("bases",
                            help="New bases (IUPAC ACGTNRYWSMKBDHV; will be uppercased).")
    p_repseq.add_argument("--force", action="store_true",
                            help="Override unsaved-changes guard.")
    p_repseq.set_defaults(fn=cmd_replace_sequence)

    p_delfeat = sub.add_parser(
        "delete-feature", help="Delete the feature at the given index.",
    )
    p_delfeat.add_argument("idx", type=int)
    p_delfeat.add_argument("--force", action="store_true",
                             help="Override unsaved-changes guard.")
    p_delfeat.set_defaults(fn=cmd_delete_feature)

    p_updfeat = sub.add_parser(
        "update-feature",
        help="Update label / type / strand of the feature at idx.",
    )
    p_updfeat.add_argument("idx", type=int)
    p_updfeat.add_argument("--label", default=None)
    p_updfeat.add_argument("--type",  default=None,
                             help="GenBank feature type (CDS, promoter, …).")
    p_updfeat.add_argument("--strand", type=int, default=None,
                             choices=[-1, 0, 1])
    p_updfeat.add_argument("--force", action="store_true")
    p_updfeat.set_defaults(fn=cmd_update_feature)

    p_getfeat = sub.add_parser(
        "get-feature", help="Detail of one feature (idx, qualifiers, …).",
    )
    p_getfeat.add_argument("idx", type=int)
    p_getfeat.set_defaults(fn=cmd_get_feature)

    p_expgb = sub.add_parser(
        "export-genbank", help="Write the loaded record to PATH as GenBank.",
    )
    p_expgb.add_argument("path", help="Output path (.gb / .gbk).")
    p_expgb.add_argument("--force", action="store_true",
                           help="Override unsaved-changes guard.")
    p_expgb.set_defaults(fn=cmd_export_genbank)

    p_expfa = sub.add_parser(
        "export-fasta", help="Write the loaded record's seq to PATH as FASTA.",
    )
    p_expfa.add_argument("path", help="Output path (.fa / .fasta / .fna).")
    p_expfa.add_argument("--force", action="store_true")
    p_expfa.set_defaults(fn=cmd_export_fasta)

    # ── Tier 2 ────────────────────────────────────────────────────────────

    p_lib = sub.add_parser(
        "list-library", help="List saved plasmid library entries.",
    )
    p_lib.add_argument("--json", action="store_true")
    p_lib.set_defaults(fn=cmd_list_library)

    p_col = sub.add_parser(
        "list-collections", help="List collections + active one.",
    )
    p_col.add_argument("--json", action="store_true")
    p_col.set_defaults(fn=cmd_list_collections)

    p_dellib = sub.add_parser(
        "delete-from-library", help="Remove a plasmid library entry by name.",
    )
    p_dellib.add_argument("name")
    p_dellib.add_argument("--force", action="store_true")
    p_dellib.set_defaults(fn=cmd_delete_from_library)

    # ── Tier 3 ────────────────────────────────────────────────────────────

    p_re = sub.add_parser(
        "list-restriction-sites",
        help="Scan record for restriction sites (default: NEB catalog, "
             "min recognition 4 bp).",
    )
    p_re.add_argument("--enzymes", nargs="*",
                       help="Limit to these enzyme names.")
    p_re.add_argument("--min-length", type=int, default=None,
                       dest="min_length",
                       help="Minimum recognition site length (default 4).")
    p_re.add_argument("--unique-only", action="store_true",
                       dest="unique_only",
                       help="Only return enzymes that cut once.")
    p_re.add_argument("--json", action="store_true")
    p_re.set_defaults(fn=cmd_list_restriction_sites)

    p_codons = sub.add_parser(
        "list-codon-tables", help="List available codon usage tables.",
    )
    p_codons.add_argument("--json", action="store_true")
    p_codons.set_defaults(fn=cmd_list_codon_tables)

    p_harm = sub.add_parser(
        "optimize-protein",
        help="Codon-optimize an AA sequence to DNA "
             "(default table: E. coli K12, taxid 83333).",
    )
    p_harm.add_argument("protein", help="1-letter AA sequence.")
    p_harm.add_argument("--table", default=None,
                          help="Codon-table taxid (see list-codon-tables).")
    p_harm.add_argument("--json", action="store_true")
    p_harm.set_defaults(fn=cmd_optimize_protein)

    p_call = sub.add_parser(
        "call",
        help="Call ANY agent endpoint by name (generic passthrough).",
    )
    p_call.add_argument("endpoint",
                         help="Endpoint name (see `tools`), e.g. "
                              "rename-plasmid.")
    p_call.add_argument("--method", choices=["GET", "POST", "get", "post"],
                         default=None,
                         help="HTTP method. Default: POST when --json is "
                              "given, else GET auto-upgraded to POST if the "
                              "endpoint is a mutation. Setting this "
                              "explicitly disables the auto-upgrade.")
    p_call.add_argument("--json", default=None, metavar="JSON",
                         help="Request body as a JSON object, e.g. "
                              "'{\"old\": \"a\", \"new\": \"b\"}'.")
    p_call.set_defaults(fn=cmd_call)

    return parser


def main(argv=None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        args.fn(args)
    except _AgentCliError as exc:
        # A call failure a subcommand didn't handle itself becomes a
        # clean message + non-zero exit (back-compat with the old
        # `sys.exit`), not a traceback. `call` handles its own.
        sys.exit(str(exc))


if __name__ == "__main__":
    main()
