"""splicecraft_search — online NCBI search operations (Phase D, layer L1).

The online network-search operations that resolve user queries against external
services, lifted out of the hub: the NCBI taxonomy name->taxid search
(_ncbi_taxid_search + _ncbi_prep_term + _ncbi_taxid_search_terms) AND the **online
search engines** — NCBI BLAST (URL API: _ncbi_blast_online / _ncbi_blast_parse_xml
/ _ncbi_blast_delete_rid) + EBI HMMER web (_hmmer_web_hmmscan / _hmmer_web_parse_json)
with their shared submit->poll->fetch->parse infra (_online_http, _online_clean_query,
_program_query_kind, _online_max_query_len, _ncbi_blast_db_for, _online_safe_*, the
_OnlineSearchCancelled exception + per-engine threading.Event cancel flags, and the
_NCBI_BLAST_*/_HMMER_*/_ONLINE_* consts). Builds on the shared SSRF-hardened primitives
in splicecraft_net (_build_hardened_url_opener / _NCBI_TIMEOUT_S) + the XML-security
parse in splicecraft_util, gates egress through the fail-closed
_state._demo_block_network_hook, stamps the User-Agent from _state._sc_version (= hub
__version__), and never logs sequence content ([INV-38]). The BlastModal / HmmerModal
UI screens (which run these in @work(thread=True) workers) STAY hub-side.

ALSO owns the **HMM-database downloader** (sweep #28) — the download / version-check /
gz-decompress / hmmpress / delete pipeline for HMMER3 databases (Pfam-A, NCBIfam, custom
URLs): _hmm_db_perform_download + _stream_download_to_path / _decompress_gz_to_path /
_hmmpress_db / _delete_hmm_db_files / _fetch_hmm_db_remote_version + the entry-dir/meta
helpers + the download-slot guard + the _HMM_DB_* / _GZIP_MAGIC / _HMMER3_MAGIC_PREFIXES
consts. Writes land under _state._HMM_DATABASES_DIR, gated by the persistence L2
chokepoint (_refuse_unauthorized_write/_delete + atomic _safe_save/load_json); the
cross-modal download slot + lock live in _state (_HMM_DB_DOWNLOAD_INFLIGHT[_LOCK]). This
part imports the same-layer L1 siblings dataaccess (_get_setting / _sanitize_hmm_db_id /
catalog) + persistence (cycle-free — neither imports search). The HMMscan *search* engine
(_hmmscan_run) + its pyhmmer helpers (_BLASTP_QUERY_ALPHABET / _pyhmmer_alignment_identity)
STAY hub-side (_sanitize_path is a shared util L0 helper, not pyhmmer-specific). Re-exported
by the hub so sc.<name> + every call site resolves unchanged.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import threading
from pathlib import Path
from typing import Any as _Any, Callable as _Callable

import splicecraft_state as _state
from splicecraft_logging import _log, _log_event
from splicecraft_util import (
    _CONTROL_CHARS_RE, _safe_xml_parse, _strip_fasta_headers, _now, _monotonic,
)
from splicecraft_net import (
    _NCBI_MAX_RESPONSE_BYTES, _NCBI_TIMEOUT_S, _build_hardened_url_opener,
    _HMM_DB_RETRY_BACKOFF_S, _hmm_db_assert_content_type_ok, _redact_url_credentials,
)
# Data-safety: the HMM-DB downloader writes multi-GB databases under
# `_state._HMM_DATABASES_DIR`; every write/delete routes through the persistence
# L2 chokepoint (`_refuse_unauthorized_*`) + atomic `_safe_save/load_json`.
from splicecraft_persistence import (
    _fsync_parent_dir, _refuse_unauthorized_write, _refuse_unauthorized_delete,
    _safe_save_json, _safe_load_json,
)
from splicecraft_dataaccess import _get_setting, _sanitize_hmm_db_id


def _ncbi_prep_term(query: str) -> str:
    """Turn a user query into an NCBI Entrez term.

    * Single token (typed-as-genus): combine an exact-taxon subtree search
      restricted to species rank with a wildcard prefix search via OR, so
      typing 'Escherichia' returns every Escherichia species (subtree hit)
      AND typing a partial like 'Escher' still matches via the wildcard.
    * Multi-word query: append '*' to the trailing token so 'Homo sapien'
      matches 'Homo sapiens' etc.
    * User-supplied wildcards or field tags pass through untouched.
    """
    q = (query or "").strip()
    if not q or "*" in q or "[" in q:
        return q
    tokens = q.split()
    if len(tokens) == 1:
        t = tokens[0]
        return f"({t}[Subtree] AND species[Rank]) OR {t}*"
    tokens[-1] = tokens[-1] + "*"
    return " ".join(tokens)


def _ncbi_taxid_search_terms(query: str) -> list[str]:
    """Cascading NCBI-taxonomy search terms, strictest first — mirrors the
    sister project ScriptoScope's `genbank_search` strategy (see
    `[SISTER]` in docs/invariants.md). The caller tries each term in
    order and stops at the first that returns hits, so an imprecise or
    partial species name still surfaces *related* taxa instead of a
    dead-end "no results".

    * Strategy 1 is `_ncbi_prep_term` — the precise-ish query (genus
      subtree + prefix wildcard for a single token; trailing wildcard
      for a multi-word name). The overwhelming majority of real queries
      resolve here, so this is identical to the pre-cascade behaviour.
    * For a MULTI-word query that comes back empty, two broader rounds
      kick in: every token prefix-wildcarded and AND-joined (tolerates a
      partial genus AS WELL AS a partial species — "Sacchar cerev"),
      then OR-joined (ANY token matches, so a stray or mistyped word
      still finds relatives — e.g. a genus typo "Homon sapiens" is
      rescued by the species epithet).

    A single-token query gets only strategy 1: its broader forms are
    already subsumed by the `OR {t}*` clause `_ncbi_prep_term` emits, so
    re-querying would just repeat the same search. A query the user has
    already steered with a `*` wildcard or a `[Field]` tag passes through
    verbatim — never second-guess an explicit Entrez query.
    """
    q = (query or "").strip()
    if not q:
        return []
    first = _ncbi_prep_term(q)
    if "*" in q or "[" in q:
        return [first]
    tokens = q.split()
    if len(tokens) < 2:
        return [first]
    terms = [first]
    for broadened in (
        " AND ".join(f"{tok}*" for tok in tokens),
        " OR ".join(f"{tok}*" for tok in tokens),
    ):
        if broadened not in terms:
            terms.append(broadened)
    return terms


def _ncbi_taxid_search(query: str, retmax: int = 200,
                      timeout: float = 15.0) -> tuple:
    """Search NCBI taxonomy for candidates matching `query`. Returns
    (hits, total_count, status_message) where each hit is
    {'taxid': str, 'name': str}. Names come from a batched esummary call
    (one round-trip for up to `retmax` ids).

    Imprecise / partial queries are handled by a cascading search
    (`_ncbi_taxid_search_terms`, mirroring ScriptoScope's `genbank_search`):
    a strict query is tried first, and if it comes back EMPTY the search
    is automatically broadened — every token prefix-wildcarded, then
    OR-joined — so a related-but-not-exact species name still surfaces
    candidates. A network / parse / oversize failure aborts the cascade
    immediately (an empty result, by contrast, just falls through to the
    next, broader strategy). Pure network — run from a worker."""
    _state._demo_block_network_hook("NCBI taxonomy search")
    import urllib.parse
    import urllib.request
    import urllib.error as _urllib_error   # sweep #25 — narrow excepts
    import xml.etree.ElementTree as ET
    q = (query or "").strip()
    if not q:
        return [], 0, "Empty query"
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    ids: list[str] = []
    total = 0
    broadened = False
    # Route NCBI taxonomy traffic through the shared hardened opener (verifying
    # SSL context, bounded redirects, https-downgrade refusal, and the private-
    # IP/SSRF host filter) instead of a raw urlopen — parity with every other
    # network path. The egress guard above already fired; the opener's own
    # guard is a harmless no-op here.
    opener = _build_hardened_url_opener()
    for idx, term in enumerate(_ncbi_taxid_search_terms(q)):
        params = urllib.parse.urlencode({
            "db": "taxonomy", "term": term,
            "retmax": str(retmax), "retmode": "xml",
        })
        try:
            req = urllib.request.Request(
                f"{base}/esearch.fcgi?{params}",
                headers={"User-Agent": "SpliceCraft/1.0"})
            with opener.open(req, timeout=timeout) as r:
                raw = r.read(_NCBI_MAX_RESPONSE_BYTES + 1)
            if len(raw) > _NCBI_MAX_RESPONSE_BYTES:
                _log.warning("NCBI esearch response exceeded %d bytes; truncating",
                              _NCBI_MAX_RESPONSE_BYTES)
                return [], 0, "NCBI returned an oversized response"
            xml_data = raw.decode("utf-8", errors="replace")
        except (OSError, _urllib_error.URLError) as exc:
            # Sweep #25 (2026-05-23): narrowed from bare `Exception` —
            # network failures land as `URLError` (timeout, DNS, refused)
            # or `OSError` (socket-level). Other exceptions (KeyError,
            # AttributeError) here would be real bugs and should
            # propagate to surface in the log + dialog. A network error
            # aborts the cascade (matches ScriptoScope's `?` propagation)
            # rather than masking an outage behind broader retries.
            _log.exception("NCBI esearch failed for %r", q)
            return [], 0, f"Network error: {exc}"
        try:
            # allow_dtd=True: NCBI eutils responses open with an EXTERNAL
            # DTD reference — `<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD
            # esearch 20060628//EN" ".../esearch.dtd">` — which has no
            # internal subset, so it stays XXE / billion-laughs safe (expat
            # never fetches the external DTD since Py 3.7.1, and there are no
            # entities to expand). Without this flag the strict parse rejects
            # EVERY real NCBI response with "XML contains DTD/ENTITY —
            # refusing to parse", which surfaced as the taxon-search "XML
            # parser error". Mirrors the BLAST-XML caller. See [PIT-19].
            root = _safe_xml_parse(xml_data, allow_dtd=True)
        except ET.ParseError as exc:
            return [], 0, f"Could not parse NCBI response: {exc}"
        cur_ids = [e.text for e in root.findall(".//Id") if e.text]
        if cur_ids:
            ids = cur_ids
            count_elem = root.find(".//Count")
            try:
                total = (int(count_elem.text)
                         if count_elem is not None and count_elem.text
                         else len(ids))
            except ValueError:
                total = len(ids)
            broadened = idx > 0
            break
    if not ids:
        return [], 0, f"No NCBI taxonomy entry for '{q}'"
    # Batched esummary: one round-trip for all retrieved ids
    names_by_id: dict[str, str] = {}
    try:
        sparams = urllib.parse.urlencode({
            "db": "taxonomy", "id": ",".join(ids), "retmode": "xml",
        })
        req = urllib.request.Request(f"{base}/esummary.fcgi?{sparams}",
                                     headers={"User-Agent": "SpliceCraft/1.0"})
        with opener.open(req, timeout=timeout) as r:
            sraw = r.read(_NCBI_MAX_RESPONSE_BYTES + 1)
        if len(sraw) > _NCBI_MAX_RESPONSE_BYTES:
            _log.warning("NCBI esummary response exceeded %d bytes; ignoring",
                          _NCBI_MAX_RESPONSE_BYTES)
            sxml = ""
        else:
            sxml = sraw.decode("utf-8", errors="replace")
        # allow_dtd=True — esummary carries the same external-DTD DOCTYPE
        # as esearch (`<!DOCTYPE eSummaryResult PUBLIC ... esummary-v1.dtd>`);
        # see the esearch parse above for why this is safe.
        sroot = _safe_xml_parse(sxml, allow_dtd=True) if sxml else None
        for doc in (sroot.findall(".//DocSum") if sroot is not None else []):
            did_el = doc.find("Id")
            if did_el is None or not did_el.text:
                continue
            did = did_el.text
            for item in doc.findall("Item"):
                if item.get("Name") == "ScientificName" and item.text:
                    names_by_id[did] = item.text
                    break
    except (OSError, _urllib_error.URLError, ET.ParseError):
        # Sweep #25 (2026-05-23): narrowed from bare `Exception`.
        _log.exception("NCBI esummary failed for ids %s", ids[:3])
    hits = [{"taxid": tid,
             "name":  names_by_id.get(tid, f"(taxid {tid})")}
            for tid in ids]
    msg = f"{total} hit(s) for '{q}'"
    if total > len(hits):
        msg = f"Showing {len(hits)} of {total} hits for '{q}' (refine to narrow)"
    if broadened:
        # No exact match — surfaced related taxa via the broadened
        # cascade. Tell the user so a loose hit isn't mistaken for a
        # precise one.
        msg += " · broadened to related names"
    return hits, total, msg


# ═══ Online search engines (NCBI BLAST URL API + EBI HMMER web) — from hub ═══
# submit -> poll -> fetch -> parse; run in BlastModal/HmmerModal workers (UI stays
# hub-side). Cancel via the per-engine threading.Event. __version__ -> _state._sc_version.

# ── Online search engines (NCBI BLAST URL API + EBI HMMER web) ───────────────
#
# Borrowed + generalised from the sister project ScriptoScope
# (`/home/seb/proteoscope/scriptoscope.py::ncbi_blastp`) per the
# "[RECIPE] borrow before respinning" playbook. Both engines follow the
# same submit → poll → fetch → parse shape and run inside BlastModal's
# `@work(thread=True)` workers so the UI keeps redrawing during the
# (potentially multi-minute) round trip. A `threading.Event` per engine
# lets the user cancel mid-flight; the NCBI side also releases the
# server-side RID so an abandoned job doesn't keep running unattended.
#
# Hardening (per the new-feature sweep convention):
#   * resp.read(MAX + 1) + bail-if-exceeded — never raw .read() ([PIT-20]).
#   * NCBI XML routed through _safe_xml_parse ([PIT-19], defangs DOCTYPE /
#     billion-laughs / unbounded nesting).
#   * HMMER JSON parsed defensively (.get() chains, several shapes) since
#     the P7Hit object is `additionalProperties: true` in the EBI schema.
#   * Friendly RuntimeErrors for timeout / unreachable host — never a raw
#     traceback to the user.
#   * Never logs sequence content ([INV-38]) — lengths / program only.

_NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
_HMMER_WEB_SUBMIT_URL = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan"
_HMMER_WEB_RESULT_URL = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/"

# Response caps ([PIT-20]). BLAST XML for a 50-hit list runs into the
# megabytes; Pfam hmmscan JSON is smaller. Both refuse a pathological
# multi-hundred-MB body (compromised / misconfigured server / MITM).
_NCBI_BLAST_MAX_RESPONSE_BYTES = 48 * 1024 * 1024
_HMMER_WEB_MAX_RESPONSE_BYTES = 24 * 1024 * 1024

# Poll cadence + overall ceiling shared by both engines. 10 s honours
# NCBI's URL-API politeness floor (don't poll a single RID more than
# ~once / 10 s); ScriptoScope used 5 s but we err polite. 300 s overall
# matches the crib.
_ONLINE_POLL_INTERVAL_S = 10
_ONLINE_MAX_WAIT_S = 300

# Programs offered by the Online tab's dropdown. The trailing "hmmscan"
# routes to the EBI HMMER engine; every other value hits NCBI BLAST.
_ONLINE_BLAST_PROGRAMS: "tuple[tuple[str, str], ...]" = (
    ("blastn   (DNA/RNA → nucleotide nt)",        "blastn"),
    ("blastp   (protein → protein nr)",           "blastp"),
    ("blastx   (translated DNA → protein nr)",    "blastx"),
    ("tblastn  (protein → translated nt)",        "tblastn"),
    ("tblastx  (translated DNA → translated nt)", "tblastx"),
    ("hmmscan  (protein → Pfam, online)",         "hmmscan"),
)

# Programs that take a nucleotide query (so a pasted RNA can be U→T
# normalised before submit) vs a protein query.
_ONLINE_NUCLEOTIDE_QUERY = frozenset({"blastn", "blastx", "tblastx"})
_ONLINE_PROTEIN_QUERY = frozenset({"blastp", "tblastn", "hmmscan"})


def _program_query_kind(program: str) -> str:
    """"nt" if the program takes a nucleotide query, else "protein".
    Spans both tabs (local blastn/blastp/hmmscan + online's five BLAST
    programs). Used to decide whether switching program invalidates the
    current query (nt↔protein) and should clear it."""
    return "nt" if program in _ONLINE_NUCLEOTIDE_QUERY else "protein"

# NCBI's hard query-length limits (chars). Submitting past these makes
# the server reject the job ("query too long"), so we refuse the search
# client-side with a clear message rather than truncate (a silently
# trimmed BLAST query returns misleading hits). Per NCBI's published
# limits: 1,000,000 for nucleotide queries; 100,000 for protein. The
# EBI HMMER hmmscan service takes a single protein — 100,000 aa is far
# beyond any real Pfam query but still caps a pathological paste.
_NCBI_BLAST_MAX_QUERY = {
    "blastn": 1_000_000, "blastx": 1_000_000, "tblastx": 1_000_000,
    "blastp": 100_000, "tblastn": 100_000,
    "hmmscan": 100_000,
}


def _online_max_query_len(program: str) -> int:
    return _NCBI_BLAST_MAX_QUERY.get(program, 100_000)

# EBI HMMER job-status buckets (case-folded). Anything non-empty that
# isn't a terminal-success / terminal-error state counts as "still
# running" and keeps the poll loop going.
_HMMER_DONE_STATUS = frozenset(
    {"SUCCESS", "DONE", "COMPLETE", "COMPLETED", "FINISHED", "OK"})
_HMMER_ERROR_STATUS = frozenset(
    {"ERROR", "ERR", "FAILURE", "FAILED", "FAIL"})


class _OnlineSearchCancelled(Exception):
    """Raised inside an online-search worker when the user cancels."""


# One cancel flag per engine. Set by the UI thread (Cancel button),
# polled by the worker via `Event.wait(timeout=…)` so cancellation is
# near-instant even mid-sleep between polls.
_ncbi_blast_cancel = threading.Event()
_hmmer_web_cancel = threading.Event()


def _ncbi_blast_db_for(program: str) -> str:
    """Default NCBI database for a BLAST program: a protein DB (`nr`) when
    the *subject* is protein (blastp / blastx), a nucleotide DB (`nt`)
    otherwise (blastn / tblastn / tblastx)."""
    return "nr" if program in ("blastp", "blastx") else "nt"


def _online_safe_float(val: "_Any") -> "float | None":
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def _online_safe_int(val: "_Any") -> "int | None":
    try:
        return int(str(val).strip())
    except (TypeError, ValueError):
        return None


def _online_clean_query(raw: str, program: str) -> str:
    """Strip FASTA headers + whitespace, upper-case, and (for nucleotide-
    query programs) normalise RNA U→T. Does NOT truncate — the caller
    enforces `_online_max_query_len(program)` and refuses oversize input
    so a trimmed query can't return misleading hits."""
    seq = _strip_fasta_headers(raw or "")
    seq = re.sub(r"\s+", "", seq).upper()
    if program in _ONLINE_NUCLEOTIDE_QUERY:
        seq = seq.replace("U", "T")
    return seq


def _online_http(url: str, *, data: "bytes | None" = None,
                 headers: "dict | None" = None,
                 timeout: int = _NCBI_TIMEOUT_S,
                 max_bytes: int = _NCBI_BLAST_MAX_RESPONSE_BYTES) -> str:
    """POST (``data`` given) / GET (``data`` None) ``url`` and return the
    decoded body. Translates network failures into friendly RuntimeErrors
    and enforces the response-size cap ([PIT-20]). HTTP error statuses
    surface as ``RuntimeError`` carrying the code so the poll loop can tell
    a transient 404-not-ready from a fatal 400."""
    import socket
    import urllib.error
    import urllib.request
    hdrs = {"User-Agent": f"SpliceCraft/{_state._sc_version}"}
    if headers:
        hdrs.update(headers)
    req = urllib.request.Request(url, data=data, headers=hdrs)
    # Sweep #30 (2026-05-28): route through the hardened opener so the
    # online BLAST/HMMER path gets the same network discipline as the
    # HMM-DB downloader + PyPI check — explicit verifying SSL context,
    # bounded redirects, and refusal of an https->http downgrade on
    # redirect. Pre-fix this used the default global opener, which would
    # follow a MITM/compromised-CDN redirect to plaintext and read the
    # RID/XML over a tamperable channel (feeding the result table). [INV-85]
    opener = _build_hardened_url_opener()
    try:
        with opener.open(req, timeout=timeout) as resp:
            raw = resp.read(max_bytes + 1)
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code} from {url}") from exc
    except urllib.error.URLError as exc:
        reason = getattr(exc, "reason", exc)
        if isinstance(reason, socket.timeout):
            raise RuntimeError(
                "Connection timed out — check your internet connection."
            ) from exc
        raise RuntimeError(f"Cannot reach server: {reason}") from exc
    except socket.timeout as exc:
        raise RuntimeError(
            "Connection timed out — check your internet connection."
        ) from exc
    if len(raw) > max_bytes:
        raise RuntimeError(
            f"Response exceeded the {max_bytes // (1024 * 1024)} MB cap — "
            f"refusing to load it. Narrow the query or hit list."
        )
    return raw.decode("utf-8", "replace")


def _ncbi_blast_delete_rid(rid: str) -> None:
    """Best-effort release of a server-side BLAST job so an abandoned /
    cancelled search doesn't keep running. Failures are swallowed — the
    job times out on NCBI's side anyway."""
    from urllib.parse import urlencode
    try:
        _online_http(
            _NCBI_BLAST_URL,
            data=urlencode({"CMD": "Delete", "RID": rid}).encode(),
            timeout=10,
        )
    except Exception:
        pass


def _ncbi_blast_online(query: str, program: str, database: str,
                       max_hits: int,
                       progress_cb: "_Callable[[str], None] | None" = None,
                       cancel_event: "threading.Event | None" = None
                       ) -> "list[dict]":
    """Run a remote NCBI BLAST via the public URL API (CMD=Put → RID →
    poll SearchInfo → CMD=Get XML → parse). Returns hit dicts shaped for
    ``BlastModal._online_render_blast``. Raises ``_OnlineSearchCancelled``
    if ``cancel_event`` fires, ``RuntimeError`` on network / server error.
    Borrowed from ScriptoScope's ``ncbi_blastp`` ([RECIPE])."""
    from urllib.parse import urlencode
    cancel = cancel_event or threading.Event()
    if progress_cb:
        progress_cb(f"Submitting {program} to NCBI…")
    put = urlencode({
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": query,
        "HITLIST_SIZE": str(max_hits),
        "FORMAT_TYPE": "XML",
    }).encode()
    put_text = _online_http(
        _NCBI_BLAST_URL, data=put,
        max_bytes=_NCBI_BLAST_MAX_RESPONSE_BYTES)
    rid = ""
    for line in put_text.splitlines():
        s = line.strip()
        if s.startswith("RID = "):
            rid = s.split("=", 1)[1].strip()
            break
    if not rid:
        raise RuntimeError(
            "NCBI did not return a job id (RID) — the query may be invalid "
            "for this program, or NCBI is rejecting traffic right now.")
    try:
        elapsed = 0
        checks = 0
        consecutive_fail = 0
        while elapsed < _ONLINE_MAX_WAIT_S:
            if cancel.wait(timeout=_ONLINE_POLL_INTERVAL_S):
                raise _OnlineSearchCancelled()
            elapsed += _ONLINE_POLL_INTERVAL_S
            checks += 1
            if progress_cb:
                progress_cb(
                    f"Waiting for NCBI {program} results "
                    f"(checked {checks}×)…")
            try:
                status_text = _online_http(
                    _NCBI_BLAST_URL,
                    data=urlencode({
                        "CMD": "Get",
                        "FORMAT_OBJECT": "SearchInfo",
                        "RID": rid,
                    }).encode(),
                    max_bytes=_NCBI_BLAST_MAX_RESPONSE_BYTES)
                consecutive_fail = 0
            except (RuntimeError, OSError) as exc:
                # Transient NCBI 5xx / network blip while the job is still
                # queued server-side: tolerate a few in a row before giving
                # up (mirrors the EBI HMMER poll loop). Pre-fix a single 502
                # aborted the search AND leaked the RID — the delete only ran
                # on the cancel path.
                consecutive_fail += 1
                if consecutive_fail >= 5:
                    raise RuntimeError(
                        f"NCBI BLAST polling failed {consecutive_fail}× in a "
                        f"row: {exc}") from exc
                _log.debug("NCBI poll transient failure %d/5: %s",
                            consecutive_fail, exc)
                continue
            if "Status=WAITING" in status_text:
                continue
            if "Status=FAILED" in status_text:
                raise RuntimeError(
                    "NCBI BLAST job failed — the server errored on the "
                    "query (check the sequence matches the program).")
            if "Status=UNKNOWN" in status_text:
                raise RuntimeError(
                    "NCBI BLAST job expired or is unknown — please retry.")
            if "Status=READY" in status_text:
                break
        else:
            raise RuntimeError(
                f"NCBI BLAST timed out after {_ONLINE_MAX_WAIT_S}s — the "
                f"query may be too large; try a shorter region.")
    except BaseException:
        # Any abort after the RID was issued (cancel, timeout, FAILED /
        # UNKNOWN, repeated poll failures) must delete the server-side job —
        # pre-fix only cancel did, leaking RIDs on every other error path.
        # Success (READY) breaks WITHOUT raising, so the RID survives for the
        # result fetch below.
        _ncbi_blast_delete_rid(rid)
        raise
    if progress_cb:
        progress_cb("Downloading NCBI results…")
    xml_text = _online_http(
        _NCBI_BLAST_URL,
        data=urlencode({
            "CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid,
        }).encode(),
        timeout=60, max_bytes=_NCBI_BLAST_MAX_RESPONSE_BYTES)
    return _ncbi_blast_parse_xml(xml_text, max_hits)


def _ncbi_blast_parse_xml(xml_text: str, max_hits: int) -> "list[dict]":
    """Parse NCBI BLAST XML into hit dicts. Routes through _safe_xml_parse
    ([PIT-19]); one row per hit (first/best HSP, mirroring the crib)."""
    import xml.etree.ElementTree as ET
    try:
        # NCBI BLAST XML opens with an external <!DOCTYPE … BlastOutput.dtd>
        # — allowed here (expat never fetches it); internal subsets remain
        # refused so this stays XXE / billion-laughs safe.
        root = _safe_xml_parse(xml_text, allow_dtd=True)
    except ET.ParseError as exc:
        raise RuntimeError(
            f"NCBI returned XML SpliceCraft couldn't parse: {exc}") from exc
    hits: "list[dict]" = []
    for hit in root.iter("Hit"):
        # Strip terminal control bytes (\x1b / OSC etc.) from server-supplied
        # text: the detail pane escapes Rich markup but NOT raw ESC, so a
        # hostile Hit_def could otherwise inject terminal escape sequences.
        hit_def = _CONTROL_CHARS_RE.sub(
            "", hit.findtext("Hit_def") or "").strip()
        acc = (_CONTROL_CHARS_RE.sub(
                   "", hit.findtext("Hit_accession") or "").strip()
               or _CONTROL_CHARS_RE.sub(
                   "", hit.findtext("Hit_id") or "").strip() or "?")
        hsp = hit.find(".//Hsp")
        if hsp is None:
            continue
        identity = _online_safe_int(hsp.findtext("Hsp_identity"))
        align_len = _online_safe_int(hsp.findtext("Hsp_align-len"))
        pct = (round(identity / align_len * 100, 1)
               if identity is not None and align_len else None)
        hits.append({
            "accession": acc,
            "description": hit_def,
            "identity_pct": pct,
            "aln_len": align_len,
            "evalue": _online_safe_float(hsp.findtext("Hsp_evalue")),
            "bit_score": _online_safe_float(hsp.findtext("Hsp_bit-score")),
            "q_start": _online_safe_int(hsp.findtext("Hsp_query-from")),
            "q_end": _online_safe_int(hsp.findtext("Hsp_query-to")),
            "s_start": _online_safe_int(hsp.findtext("Hsp_hit-from")),
            "s_end": _online_safe_int(hsp.findtext("Hsp_hit-to")),
        })
        if len(hits) >= max_hits:
            break
    return hits


def _hmmer_web_parse_json(obj: "_Any", max_hits: int) -> "list[dict]":
    """Pull Pfam hits out of an EBI HMMER result body. Defensive across
    `result.hits`, `results.hits`, and top-level `hits`.

    The human-readable family name + description live in the hit's
    ``metadata`` sub-object (``identifier`` / ``description``), NOT in the
    top-level ``name`` (an internal numeric id) or ``desc`` (always null in
    the live v1 API). We read metadata first and fall back to the
    top-level / legacy keys so canned fixtures without metadata still
    parse. ``clan`` / ``type`` / ``external_link`` enrich the detail pane;
    ``included`` flags whether the hit cleared Pfam's inclusion (gathering)
    threshold vs. merely the reporting threshold."""
    if not isinstance(obj, dict):
        return []
    hits = None
    for container_key in ("result", "results"):
        container = obj.get(container_key)
        if isinstance(container, dict) and isinstance(
                container.get("hits"), list):
            hits = container["hits"]
            break
    if hits is None and isinstance(obj.get("hits"), list):
        hits = obj["hits"]
    if not isinstance(hits, list):
        return []
    out: "list[dict]" = []
    for h in hits:
        if not isinstance(h, dict):
            continue
        md = h.get("metadata")
        if not isinstance(md, dict):
            md = {}
        # Strip terminal control bytes from server-supplied strings (the
        # detail pane escapes Rich markup but NOT raw ESC / OSC sequences).
        acc = _CONTROL_CHARS_RE.sub("", str(
            h.get("acc") or md.get("accession")
            or h.get("accession") or "?"))
        # `metadata.identifier` is the Pfam family name (e.g. "Pkinase");
        # the top-level `name` is an internal numeric id — avoid it unless
        # nothing better exists.
        name = _CONTROL_CHARS_RE.sub("", str(
            md.get("identifier") or md.get("id")
            or h.get("name") or acc))
        desc = _CONTROL_CHARS_RE.sub("", str(
            md.get("description") or h.get("desc")
            or h.get("description") or ""))
        ndom = _online_safe_int(
            h.get("ndom") if h.get("ndom") is not None
            else h.get("nincluded") if h.get("nincluded") is not None
            else h.get("nreported"))
        if ndom is None:
            doms = h.get("domains")
            ndom = len(doms) if isinstance(doms, list) else 0
        out.append({
            "acc": acc,
            "name": name,
            "description": desc,
            "evalue": _online_safe_float(
                h.get("evalue") if h.get("evalue") is not None
                else h.get("eval")),
            "bit_score": _online_safe_float(
                h.get("score") if h.get("score") is not None
                else h.get("bitscore")),
            "n_dom": ndom,
            # Sweep #30 (2026-05-28): strip control bytes here too — these
            # three EBI-supplied fields also reach the markup=True detail
            # pane (clan/type/link) and were the gap the acc/name/desc
            # strips above left open. [INV-85]
            "clan": _CONTROL_CHARS_RE.sub("", str(md.get("clan") or "")),
            "type": _CONTROL_CHARS_RE.sub("", str(md.get("type") or "")),
            "link": _CONTROL_CHARS_RE.sub(
                "", str(md.get("external_link") or "")),
            "included": bool(h.get("is_included", True)),
        })
        if len(out) >= max_hits:
            break
    return out


def _hmmer_web_hmmscan(protein: str, max_hits: int,
                       progress_cb: "_Callable[[str], None] | None" = None,
                       cancel_event: "threading.Event | None" = None
                       ) -> "list[dict]":
    """Run a remote hmmscan vs Pfam via the EBI HMMER web API (POST search
    → poll result/{id} until SUCCESS → parse). Returns Pfam hit dicts.
    Raises ``_OnlineSearchCancelled`` on cancel, ``RuntimeError`` on
    network / server error."""
    cancel = cancel_event or threading.Event()
    if progress_cb:
        progress_cb("Submitting hmmscan to EBI HMMER…")
    body = json.dumps({"input": protein, "database": "pfam"}).encode()
    submit_text = _online_http(
        _HMMER_WEB_SUBMIT_URL, data=body,
        headers={"Content-Type": "application/json",
                 "Accept": "application/json"},
        max_bytes=_HMMER_WEB_MAX_RESPONSE_BYTES)
    try:
        submit = json.loads(submit_text)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            "EBI HMMER returned an unparseable submit response.") from exc
    # Some deployments return hits inline on submit — short-circuit.
    inline = _hmmer_web_parse_json(submit, max_hits)
    if inline:
        return inline
    job_id = ""
    if isinstance(submit, dict):
        job_id = str(submit.get("id") or submit.get("uuid") or "").strip()
    if not job_id:
        raise RuntimeError("EBI HMMER did not return a job id.")
    result_url = _HMMER_WEB_RESULT_URL + job_id
    elapsed = 0
    checks = 0
    consecutive_fail = 0
    while elapsed < _ONLINE_MAX_WAIT_S:
        if cancel.wait(timeout=_ONLINE_POLL_INTERVAL_S):
            raise _OnlineSearchCancelled()
        elapsed += _ONLINE_POLL_INTERVAL_S
        checks += 1
        if progress_cb:
            progress_cb(
                f"Waiting for Pfam hmmscan results (checked {checks}×)…")
        try:
            rtext = _online_http(
                result_url,
                headers={"Accept": "application/json"},
                max_bytes=_HMMER_WEB_MAX_RESPONSE_BYTES)
            obj = json.loads(rtext)
            consecutive_fail = 0
        except (RuntimeError, json.JSONDecodeError):
            # 404-not-ready / transient flake — keep polling, but don't
            # spin forever against a persistent failure.
            consecutive_fail += 1
            if consecutive_fail >= 5:
                raise RuntimeError(
                    "EBI HMMER kept returning errors — try again later.")
            continue
        status = str(
            (obj.get("status") if isinstance(obj, dict) else "") or ""
        ).upper()
        if status in _HMMER_ERROR_STATUS:
            raise RuntimeError("EBI HMMER reported the job failed.")
        if status and status not in _HMMER_DONE_STATUS:
            continue  # PENDING / RUNNING / QUEUED — still cooking
        return _hmmer_web_parse_json(obj, max_hits)
    raise RuntimeError(
        f"EBI HMMER timed out after {_ONLINE_MAX_WAIT_S}s — try again "
        f"or use the Local tab with a downloaded Pfam database.")


# ── HMM database downloader (sweep #28) — relocated from the hub (Phase D) ──
# The download / version-check / decompress / hmmpress / delete pipeline for
# HMMER3 databases (Pfam-A, NCBIfam, custom URLs). Writes land under
# `_state._HMM_DATABASES_DIR` and are L2-chokepoint-guarded
# (`_refuse_unauthorized_write` / `_refuse_unauthorized_delete`). The HMMscan
# *search* engine (`_hmmscan_run`) + the Blast/Hmmer UI modals stay hub-side.

# Per-download caps. The 4 GB ceiling is well above the largest
# realistic single-HMM database (NCBIfam compressed is ~600 MB,
# uncompressed ~3.5 GB; Pfam-A compressed ~300 MB, uncompressed ~2 GB).
# A 100 GB typo'd URL (or a hostile mirror serving /dev/zero) bails
# at the cap instead of OOM'ing the disk.
_HMM_DB_DOWNLOAD_MAX_BYTES = 4 * 1024 * 1024 * 1024
# Version-file cap is tiny — Pfam.version.gz is ~150 bytes.
_HMM_DB_VERSION_MAX_BYTES = 64 * 1024

# Re-check the remote version no more than once every 24h. Pfam
# releases happen every 6-12 months; NCBIfam every few weeks. The
# cache is keyed by entry id so a manual "Check now" on one DB
# doesn't reset the timer on the others.
_HMM_DB_VERSION_CHECK_TTL_S = 24 * 3600


def _hmm_db_entry_dir(entry_id: str) -> Path:
    """Per-DB directory: `<DATA_DIR>/hmm_databases/<id>/`. The id has
    already been sanitised via `_sanitize_hmm_db_id` before any call
    site reaches here; we still re-sanitise defensively so a future
    caller that forgets the sanitisation can't talk us into writing
    outside the databases dir."""
    safe = _sanitize_hmm_db_id(entry_id) or "_unknown"
    return _state._HMM_DATABASES_DIR / safe


def _hmm_db_meta_path(entry_id: str) -> Path:
    """Per-DB metadata file: download timestamp, sha256, n_profiles,
    pressed-state. Routed through `_safe_save_json` so it gets the
    same atomic-write + backup chain as every other persisted file."""
    return _hmm_db_entry_dir(entry_id) / "meta.json"


def _hmm_db_hmm_path(entry_id: str) -> Path:
    """The decompressed HMM file. pyhmmer's hmmscan reads this (or
    the `.h3*` sibling index files if pressed)."""
    return _hmm_db_entry_dir(entry_id) / "db.hmm"






# ── HMM DB local state ────────────────────────────────────────────────
# Each downloaded DB gets a `meta.json` sibling tracking:
#   * `version`         — remote release identifier (e.g. "Pfam 37.4")
#                         or HTTP Last-Modified fallback for DBs
#                         without an explicit version file
#   * `downloaded_at`   — ISO-8601, tz-aware
#   * `sha256`          — checksum of the source `.hmm.gz` (so a
#                         truncated download is detectable)
#   * `n_profiles`      — count of HMMs in the file (smoke check)
#   * `pressed`         — bool: are the `.h3i/.h3m/.h3p/.h3f`
#                         siblings present?
#   * `last_remote_version_check_ts` — wall-clock monotonic seconds
#                         for the 24h auto-check cache
#   * `last_remote_version` — what the remote said last time we polled


def _load_hmm_db_local_meta(entry_id: str) -> "dict | None":
    """Read the per-DB metadata. Returns None if no download has
    happened yet OR if the file is corrupted (we don't try to
    recover; the UI will show "not downloaded" and a fresh download
    overwrites)."""
    path = _hmm_db_meta_path(entry_id)
    if not path.exists():
        return None
    try:
        entries, _warning = _safe_load_json(path, "HMM DB meta")
    except (OSError, ValueError):
        return None
    if not entries:
        return None
    meta = entries[0]
    return meta if isinstance(meta, dict) else None


def _save_hmm_db_local_meta(entry_id: str, meta: dict) -> None:
    """Persist per-DB metadata. The entry's directory is created
    on demand (mkdir parents=True, exist_ok=True). The L2 chokepoint
    fires via `_safe_save_json` so an unsandboxed probe can't talk
    its way past the gate."""
    d = _hmm_db_entry_dir(entry_id)
    d.mkdir(parents=True, exist_ok=True)
    path = _hmm_db_meta_path(entry_id)
    _safe_save_json(path, [meta], "HMM DB meta")


def _is_hmm_db_downloaded(entry_id: str) -> bool:
    """True iff the decompressed HMM file is present on disk. Doesn't
    require the pressed `.h3*` files (those are an optimisation; a
    raw `.hmm` is still usable, just slower)."""
    return _hmm_db_hmm_path(entry_id).exists()


def _hmm_db_pressed(entry_id: str) -> bool:
    """True iff `pyhmmer.hmmer.hmmpress` has been run (all four
    index files present). pyhmmer can hmmscan an unpressed HMM,
    just much slower on Pfam-scale DBs."""
    base = _hmm_db_hmm_path(entry_id)
    return all((base.parent / f"{base.name}.{ext}").exists()
               for ext in ("h3i", "h3m", "h3p", "h3f"))


# ── Network hardening (sweep #28) ────────────────────────────────────
# These helpers are shared by every internet-facing HMM helper below
# (`_fetch_hmm_db_remote_version`, `_stream_download_to_path`) so the
# hardening lands in one place. Rationale per defense layer:
#
#   * **HTTPS only by default** — every builtin URL is https; a user
#     adding a custom URL gets a debug log + the download proceeds
#     (we don't refuse http outright because some legacy mirrors only
#     serve over plain HTTP). The setting `hmm_db_allow_http` (off by
#     default) gates this; absent it, http URLs are rejected.
#   * **Explicit SSL context** — `ssl.create_default_context()` so
#     cert validation uses the system trust store. Pre-sweep `urllib.
#     request.urlopen` defaulted to whatever Python's compiled-in
#     defaults were; on a system with a broken CA bundle this could
#     silently fall through to no validation. Explicit context makes
#     the validation intent obvious + greppable.
#   * **Bounded redirects (5 max)** — urllib's default is 10. We tighten
#     to 5 so a hostile mirror can't bounce us through an arbitrarily
#     long chain of redirects (each one is a network round-trip + a
#     fresh chance for partial-response weirdness).
#   * **Content-Type guard** — refuse responses whose Content-Type is
#     `text/html` or `application/json`. Both indicate the server
#     returned an error page (CDN block, captcha, "site under
#     maintenance" wrapper) instead of the binary payload we asked
#     for. Pre-sweep we'd happily save 12 KB of HTML as `db.hmm.gz`
#     and bewilder the user with a "not a gzip file" error from
#     hmmpress later.
#   * **Magic-byte verification** — gunzip output starts with the
#     ASCII tag `HMMER3/f` (or `HMMER2.0`, rare). We sniff the first
#     line and reject if missing — catches a download that was the
#     right size but the wrong file (CDN substitution, mirror drift,
#     redirect to a `README.txt`).
#   * **n_profiles ≥ 1** — after hmmpress, if the press returns 0
#     profiles, the file was syntactically valid HMMER3 but empty;
#     refuse + clean up. Real Pfam-A has ~21000 profiles; NCBIfam ~14000.
#   * **Retry on transient failures** — one retry with 250 ms backoff
#     for `URLError` / `socket.timeout` (mirrors `[INV-37]`'s
#     `_fetch_latest_pypi_version` + `fetch_genbank`). Persistent
#     errors surface immediately.
#   * **Disk-space pre-check** — `shutil.disk_usage` before any
#     download. Refuse if free < expected × 2.5 (the ×2.5 reserves
#     headroom for decompression + hmmpress on top of the .gz size).
#   * **Bounded URL in logs** — `_redact_url_credentials` strips any
#     embedded `user:pass@host` from a logged URL (a paranoid user
#     could paste a URL with creds; we don't want them in the
#     diagnostic bundle).
#   * **Mid-flight cancellation** — every long phase (download chunk
#     loop, decompress chunk loop) calls `cancel_check_cb()` between
#     iterations. The worker's cb peeks `is_mounted` so a closed
#     modal aborts cleanly instead of running to completion against
#     a no-longer-visible UI.
#   * **Global cross-modal download lock** — `_HMM_DB_DOWNLOAD_INFLIGHT`
#     dict-of-id-→-bool gates concurrent downloads of the same DB
#     even across re-opened modal instances. Each modal still has its
#     own `_active_downloads` set for per-modal UI state.
#   * **Atomic hmmpress cleanup** — `_hmmpress_db` removes any
#     partial `.h3*` siblings on failure so a half-pressed state
#     doesn't get treated as "ready" by `_hmm_db_pressed`.


_HMM_DB_NETWORK_TIMEOUT_S = 60      # connect + read; per-op
# Disk-space reserve multiplier: total free space must be ≥ expected
# download size × this. Covers decompression (~7×) + hmmpress (~0.3×)
# headroom on top of the .gz. 2.5 is comfortable for compressed HMM
# files; for an unknown-size download we use a fixed 5 GB reserve.
_HMM_DB_DISK_RESERVE_MULTIPLIER = 2.5
_HMM_DB_DISK_RESERVE_DEFAULT_BYTES = 5 * 1024 * 1024 * 1024

def _hmm_db_acquire_download_slot(entry_id: str) -> bool:
    """Reserve the download slot for `entry_id`. Returns True if
    acquired (caller must release), False if another flow holds it."""
    with _state._HMM_DB_DOWNLOAD_INFLIGHT_LOCK:
        if entry_id in _state._HMM_DB_DOWNLOAD_INFLIGHT:
            return False
        _state._HMM_DB_DOWNLOAD_INFLIGHT.add(entry_id)
        return True


def _hmm_db_release_download_slot(entry_id: str) -> None:
    """Release the download slot. Idempotent — safe in `finally`."""
    with _state._HMM_DB_DOWNLOAD_INFLIGHT_LOCK:
        _state._HMM_DB_DOWNLOAD_INFLIGHT.discard(entry_id)


def _hmm_db_check_disk_space(target: Path,
                              expected_bytes: "int | None") -> None:
    """Raise OSError if there's insufficient free space at `target`'s
    filesystem. `expected_bytes` is the Content-Length of the planned
    download (None → use `_HMM_DB_DISK_RESERVE_DEFAULT_BYTES`)."""
    needed = (int(expected_bytes * _HMM_DB_DISK_RESERVE_MULTIPLIER)
              if expected_bytes
              else _HMM_DB_DISK_RESERVE_DEFAULT_BYTES)
    # Walk up to an existing ancestor for the disk_usage probe; the
    # entry dir itself may not exist yet.
    probe = target if target.exists() else target.parent
    while probe and not probe.exists():
        probe = probe.parent
    if probe is None:
        return  # disk-space check is best-effort
    try:
        usage = shutil.disk_usage(str(probe))
    except OSError:
        return
    if usage.free < needed:
        raise OSError(
            f"insufficient disk space at {probe}: need ≥ "
            f"{needed // (1024*1024):,} MB free (have "
            f"{usage.free // (1024*1024):,} MB). The download is "
            f"refused before any bytes hit disk."
        )


def _hmm_db_user_agent() -> str:
    """Identify our traffic to EBI / NCBI per their crawler policies.
    Mirrors the NCBI Entrez `tool` identifier added in sweep #27."""
    return f"SpliceCraft/{_state._sc_version} (HMM database client)"


def _hmm_db_url_scheme_ok(url: str) -> "str | None":
    """Return None if the URL passes the http/https policy, else a
    user-facing error string. http URLs are rejected by default
    (any custom-URL builtin override should be https); enable via
    `settings.hmm_db_allow_http = true` for the rare legacy mirror."""
    if not isinstance(url, str) or not url:
        return "URL missing"
    if url.startswith("https://"):
        return None
    if url.startswith("http://"):
        if _get_setting("hmm_db_allow_http", False):
            return None
        return ("URL uses plain http:// — refused for safety. Enable "
                "the `hmm_db_allow_http` setting if you really need "
                "to download from an http-only mirror.")
    return "URL must start with https:// (or http:// with opt-in)"


def _parse_pfam_version_text(text: str) -> str:
    """Pfam's `Pfam.version` file contains lines like:

        Pfam release       : 37.4
        Pfam-A families    : 21978
        Date               : 2025-02
        Based on UniProtKB : 2024_06

    Return the bare release identifier (`"37.4"`) or the full first
    non-empty line if the format ever drifts. We never raise on a
    malformed version file — degraded "version present but I can't
    parse it" is better than the modal blocking on an EBI hiccup."""
    if not isinstance(text, str):
        return ""
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        # Match `Pfam release : 37.4` (case-insensitive, lenient
        # whitespace + colon).
        low = line.lower()
        if low.startswith("pfam release"):
            _, _, rest = line.partition(":")
            v = rest.strip()
            if v:
                return v[:64]
        # Fallback for hypothetical future formats: take the first
        # numeric token that looks like a version.
    # Last resort: first non-empty trimmed line, capped.
    for line in text.splitlines():
        s = line.strip()
        if s:
            return s[:64]
    return ""


def _fetch_hmm_db_remote_version(entry: dict) -> "tuple[str, str]":
    """Return `(version_string, source)` for the remote DB. `source`
    is one of `"version_file"` (parsed from a dedicated version URL),
    `"last_modified"` (HTTP HEAD on the main URL), or `""` on failure.

    Hardened (sweep #28):
      * HTTPS scheme check (bypass via `hmm_db_allow_http` setting).
      * Hardened opener (bounded redirects, explicit SSL context).
      * One retry with 250 ms backoff on transient network failures
        (mirrors `[INV-37]`).
      * Cap-bounded read (`_HMM_DB_VERSION_MAX_BYTES = 64 KB`).
      * URL credentials redacted in any log line.

    Best-effort: every error path returns `("", "")` rather than
    raising, because the UI should degrade to "couldn't check"
    rather than failing the modal mount."""
    import gzip
    import socket
    import time as _time
    import urllib.error

    opener = _build_hardened_url_opener()
    entry_id = entry.get("id") or "?"
    version_url = (entry.get("version_url") or "").strip()

    def _redacted(u: str) -> str:
        return _redact_url_credentials(u)[:200]

    def _try_get(url: str,
                 *, method: str,
                 max_bytes: int) -> "tuple[bytes, dict] | None":
        """One-shot fetch with retry. Returns (raw, headers) or None."""
        from urllib.request import Request as _Req
        scheme_err = _hmm_db_url_scheme_ok(url)
        if scheme_err is not None:
            _log.warning(
                "HMM DB version check %r refused: %s (%s)",
                entry_id, scheme_err, _redacted(url),
            )
            return None
        last_exc: "BaseException | None" = None
        for attempt in range(2):
            try:
                req = _Req(
                    url, method=method,
                    headers={"User-Agent": _hmm_db_user_agent()},
                )
                with opener.open(
                    req, timeout=_HMM_DB_NETWORK_TIMEOUT_S,
                ) as resp:
                    # Refuse an HTML / JSON / XML error page served with a 200
                    # at the version URL — pre-fix a CDN interstitial flowed
                    # into gzip.decompress / the version parser and surfaced
                    # its first line as the "remote version", mislabelling the
                    # DB as out-of-date. HEAD has no body to check.
                    if method == "GET":
                        _hmm_db_assert_content_type_ok(resp, url)
                    raw = resp.read(max_bytes + 1) if method == "GET" else b""
                    if method == "GET" and len(raw) > max_bytes:
                        _log.warning(
                            "HMM DB version file at %s exceeds %d byte "
                            "cap; aborting check.",
                            _redacted(url), max_bytes,
                        )
                        return None
                    return (raw, dict(resp.headers))
            except (urllib.error.URLError, socket.timeout, OSError,
                    ValueError) as exc:
                last_exc = exc
                if attempt == 0:
                    _log.debug(
                        "HMM DB version probe attempt 1 failed for %s "
                        "(%s); retrying after %.1fs",
                        _redacted(url), exc, _HMM_DB_RETRY_BACKOFF_S,
                    )
                    _log_event(
                        "net.retry", endpoint="hmm_db_version",
                        url=_redacted(url),
                        exc_type=type(exc).__name__,
                    )
                    _time.sleep(_HMM_DB_RETRY_BACKOFF_S)
                    continue
                _log.warning(
                    "HMM DB version check %r failed (%s); url=%s",
                    entry_id, exc, _redacted(url),
                )
                return None
        _ = last_exc   # silence "unused" — already logged
        return None

    if version_url:
        out = _try_get(version_url, method="GET",
                        max_bytes=_HMM_DB_VERSION_MAX_BYTES)
        if out is not None:
            raw, _hdrs = out
            try:
                if version_url.endswith(".gz"):
                    # Sweep #30 (2026-05-28): bound the DECOMPRESSED size.
                    # `raw` is 64 KB-capped, but gzip can expand ~1000x, so
                    # a crafted bomb would balloon to tens of MB in RAM on
                    # every (24h-gated) modal mount. The real version file
                    # is ~150 bytes; stream-decompress with the same cap and
                    # bail on overflow rather than one-shot gzip.decompress
                    # (which would inflate the whole bomb). [INV-84]
                    import io as _io
                    with gzip.GzipFile(fileobj=_io.BytesIO(raw)) as _gz:
                        _dec = _gz.read(_HMM_DB_VERSION_MAX_BYTES + 1)
                    if len(_dec) > _HMM_DB_VERSION_MAX_BYTES:
                        _log.warning(
                            "HMM DB version file at %s decompresses past the "
                            "%d-byte cap; ignoring (possible gzip bomb).",
                            _redacted(version_url),
                            _HMM_DB_VERSION_MAX_BYTES,
                        )
                        text = ""
                    else:
                        text = _dec.decode("utf-8", errors="replace")
                else:
                    text = raw.decode("utf-8", errors="replace")
            except (OSError, ValueError, EOFError):
                text = raw.decode("utf-8", errors="replace")
            parsed = _parse_pfam_version_text(text)
            if parsed:
                return (parsed, "version_file")

    main_url = (entry.get("url") or "").strip()
    if not main_url:
        return ("", "")
    out = _try_get(main_url, method="HEAD", max_bytes=0)
    if out is not None:
        _, hdrs = out
        last_mod = (hdrs.get("Last-Modified")
                    or hdrs.get("last-modified") or "")
        if last_mod:
            return (last_mod.strip()[:64], "last_modified")
    return ("", "")


def _hmm_db_local_version(entry_id: str) -> str:
    """Return the locally-stored version string, or "" if no download
    has happened."""
    meta = _load_hmm_db_local_meta(entry_id)
    if meta is None:
        return ""
    v = meta.get("version")
    return v if isinstance(v, str) else ""


def _hmm_db_should_check_remote(entry_id: str) -> bool:
    """24h cache: if the last remote check was less than `_HMM_DB_
    VERSION_CHECK_TTL_S` ago, skip the network round-trip."""
    meta = _load_hmm_db_local_meta(entry_id)
    if meta is None:
        return True
    last = meta.get("last_remote_version_check_ts")
    if not isinstance(last, (int, float)):
        return True
    return (_monotonic() - last) > _HMM_DB_VERSION_CHECK_TTL_S


def _record_hmm_db_remote_version(entry_id: str,
                                    remote_version: str) -> None:
    """Stamp the remote version + check timestamp into the per-DB
    meta. Creates the meta if no download has happened yet (so a
    pure version-check user gets the 24h cache benefit even before
    they download anything)."""
    meta = _load_hmm_db_local_meta(entry_id) or {
        "id":             entry_id,
        "version":        "",
        "downloaded_at":  "",
        "sha256":         "",
        "n_profiles":     0,
        "pressed":        False,
    }
    meta["last_remote_version"] = (
        remote_version or ""
    )[:64]
    meta["last_remote_version_check_ts"] = _monotonic()
    try:
        _save_hmm_db_local_meta(entry_id, meta)
    except (OSError, RuntimeError) as exc:
        _log.warning(
            "HMM DB version stamp save failed for %r: %s",
            entry_id, exc,
        )


# ── Download worker ──────────────────────────────────────────────────


_GZIP_MAGIC = b"\x1f\x8b"


def _stream_download_to_path(url: str, dest: Path,
                              *,
                              max_bytes: int,
                              progress_cb=None,
                              cancel_check_cb=None,
                              chunk_size: int = 64 * 1024) -> str:
    """Stream `url` to `dest` atomically. Returns the SHA-256 hex.

    Hardened (sweep #28):
      * HTTPS-only by default; http requires `hmm_db_allow_http` opt-in.
      * Hardened opener (bounded redirects, explicit SSL context).
      * One retry with 250 ms backoff on transient network failures.
      * Content-Type guard refuses HTML/JSON error pages.
      * Disk-space pre-check (`needed = Content-Length × 2.5` or 5 GB).
      * `cancel_check_cb()` polled between chunks for cooperative abort
        (returns True to abort).
      * Gzip magic bytes verified on the first read so a malformed
        download fails fast instead of after writing GB.
      * Atomic write: tmp → fsync → os.replace → fsync_parent.
      * L2 chokepoint via `_refuse_unauthorized_write`.

    Raises:
      ValueError on cap exceedance / bad content-type / bad magic.
      OSError on network / disk failure / cancellation.
    """
    import hashlib
    import socket
    import time as _time
    import urllib.error

    redacted = _redact_url_credentials(url)

    scheme_err = _hmm_db_url_scheme_ok(url)
    if scheme_err is not None:
        raise ValueError(scheme_err)

    dest.parent.mkdir(parents=True, exist_ok=True)
    _refuse_unauthorized_write(dest, "hmm db download")

    tmp = dest.with_name(dest.name + ".download_tmp")
    try:
        tmp.unlink()
    except FileNotFoundError:
        pass

    opener = _build_hardened_url_opener()
    from urllib.request import Request as _Req

    def _open_with_retry():
        """One retry with 250 ms backoff, only on transient errors
        (URLError, socket.timeout). Permanent errors (HTTPError 4xx,
        ValueError) raise immediately."""
        last_exc: "BaseException | None" = None
        for attempt in range(2):
            try:
                req = _Req(
                    url,
                    headers={"User-Agent": _hmm_db_user_agent()},
                )
                return opener.open(
                    req, timeout=_HMM_DB_NETWORK_TIMEOUT_S,
                )
            except urllib.error.HTTPError:
                # 4xx/5xx is a permanent server response; don't retry
                # (we'd just hit the same response).
                raise
            except (urllib.error.URLError, socket.timeout) as exc:
                last_exc = exc
                if attempt == 0:
                    _log.debug(
                        "HMM DB download attempt 1 failed for %s (%s); "
                        "retrying after %.1fs",
                        redacted, exc, _HMM_DB_RETRY_BACKOFF_S,
                    )
                    _log_event(
                        "net.retry", endpoint="hmm_db_download",
                        url=redacted, exc_type=type(exc).__name__,
                    )
                    _time.sleep(_HMM_DB_RETRY_BACKOFF_S)
                    continue
                raise
        assert last_exc is not None
        raise last_exc  # unreachable; for the type checker

    sha = hashlib.sha256()
    bytes_so_far = 0
    cancelled = False
    saw_magic = False
    try:
        resp = _open_with_retry()
        try:
            _hmm_db_assert_content_type_ok(resp, url)
            total_header = resp.headers.get("Content-Length")
            total: "int | None"
            try:
                total = int(total_header) if total_header else None
            except (TypeError, ValueError):
                total = None
            if total is not None and total > max_bytes:
                raise ValueError(
                    f"refusing download: server reports "
                    f"{total:,} bytes > cap {max_bytes:,}"
                )
            # Disk-space pre-check (after we know total, before we
            # commit any bytes to disk).
            _hmm_db_check_disk_space(dest, total)
            first_chunk = True
            with open(tmp, "wb") as fh:
                while True:
                    if (cancel_check_cb is not None
                            and bool(cancel_check_cb())):
                        cancelled = True
                        raise OSError("download cancelled by user")
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    bytes_so_far += len(chunk)
                    if bytes_so_far > max_bytes:
                        raise ValueError(
                            f"download exceeded {max_bytes:,} byte cap "
                            f"({bytes_so_far:,} read) — likely zip-bomb "
                            f"or wrong URL"
                        )
                    if first_chunk:
                        first_chunk = False
                        if not chunk.startswith(_GZIP_MAGIC):
                            raise ValueError(
                                "first bytes don't match gzip magic "
                                "(0x1f8b) — server didn't return a "
                                "gzipped HMM file"
                            )
                        saw_magic = True
                    sha.update(chunk)
                    fh.write(chunk)
                    if progress_cb is not None:
                        try:
                            progress_cb(bytes_so_far, total)
                        except Exception:
                            _log.exception(
                                "HMM DB download: progress callback raised",
                            )
                fh.flush()
                try:
                    os.fsync(fh.fileno())
                except OSError:
                    pass
        finally:
            try:
                resp.close()
            except Exception:
                pass
        if not saw_magic:
            raise ValueError(
                "download empty — server returned 0 bytes"
            )
        os.replace(str(tmp), str(dest))
        _fsync_parent_dir(dest)
    except Exception:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass
        if cancelled:
            _log_event(
                "hmm_db.download.cancelled",
                url=redacted, bytes_read=bytes_so_far,
            )
        raise
    return sha.hexdigest()


_HMMER3_MAGIC_PREFIXES: tuple[bytes, ...] = (b"HMMER3/", b"HMMER2.0")


def _decompress_gz_to_path(src: Path, dest: Path,
                            *,
                            max_bytes: int,
                            cancel_check_cb=None) -> None:
    """Streaming gunzip with size cap + HMMER3 magic verification.

    Hardened (sweep #28):
      * `max_bytes` bounds the DECOMPRESSED output (legitimate Pfam-A
        compression ratio is ~7:1; we cap at 4 GB for the .hmm).
      * HMMER3 magic check on first chunk — rejects a file that
        decompresses fine but isn't an HMM (CDN-served README, etc.).
      * `cancel_check_cb()` polled between chunks.
      * Disk-space pre-check via `_hmm_db_check_disk_space` (uses
        compressed file size × 10 as a conservative reserve).
      * L2 chokepoint via `_refuse_unauthorized_write`.

    Raises ValueError on cap exceedance / bad magic / corrupt gzip;
    OSError on disk failure or cancellation."""
    import gzip

    _refuse_unauthorized_write(dest, "hmm db gunzip")
    # Reserve roughly 10× the .gz size (Pfam-A is ~7:1 compressed;
    # rounded up for safety + hmmpress headroom).
    try:
        src_size = src.stat().st_size
        _hmm_db_check_disk_space(dest, src_size * 10)
    except OSError:
        # Disk check itself can fail (e.g. on a path that doesn't
        # exist yet); re-raise.
        raise

    tmp = dest.with_name(dest.name + ".gz_tmp")
    try:
        tmp.unlink()
    except FileNotFoundError:
        pass
    written = 0
    cancelled = False
    saw_magic = False
    try:
        with gzip.open(str(src), "rb") as fh_in, open(tmp, "wb") as fh_out:
            first = True
            while True:
                if (cancel_check_cb is not None
                        and bool(cancel_check_cb())):
                    cancelled = True
                    raise OSError("decompression cancelled by user")
                try:
                    chunk = fh_in.read(1024 * 1024)
                except (OSError, EOFError) as exc:
                    raise ValueError(
                        f"corrupt gzip data after {written:,} bytes: "
                        f"{exc}"
                    ) from exc
                if not chunk:
                    break
                written += len(chunk)
                if written > max_bytes:
                    raise ValueError(
                        f"decompressed output exceeded {max_bytes:,} "
                        f"byte cap (likely zip-bomb / bad upload)"
                    )
                if first:
                    first = False
                    if not any(chunk.startswith(m)
                               for m in _HMMER3_MAGIC_PREFIXES):
                        raise ValueError(
                            "decompressed file isn't a HMMER profile "
                            "database — expected first bytes to start "
                            "with 'HMMER3/' or 'HMMER2.0'"
                        )
                    saw_magic = True
                fh_out.write(chunk)
            fh_out.flush()
            try:
                os.fsync(fh_out.fileno())
            except OSError:
                pass
        if not saw_magic:
            raise ValueError(
                "decompressed file is empty — source .hmm.gz holds no "
                "HMMER content"
            )
        os.replace(str(tmp), str(dest))
        _fsync_parent_dir(dest)
    except Exception:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass
        if cancelled:
            _log_event(
                "hmm_db.decompress.cancelled",
                src=str(src.name), bytes_written=written,
            )
        raise


def _cleanup_pressed_files(hmm_path: Path) -> None:
    """Remove any partial `.h3i/.h3m/.h3p/.h3f` siblings of `hmm_path`.
    Called from `_hmmpress_db`'s failure path so a half-pressed state
    doesn't get treated as "ready" by `_hmm_db_pressed`."""
    for ext in ("h3i", "h3m", "h3p", "h3f"):
        sibling = hmm_path.parent / f"{hmm_path.name}.{ext}"
        _refuse_unauthorized_delete(sibling, "HMM DB pressed index")
        try:
            sibling.unlink()
        except FileNotFoundError:
            pass
        except OSError as exc:
            _log.warning(
                "HMM DB cleanup: couldn't remove %s (%s)",
                sibling, exc,
            )


def _hmmpress_db(hmm_path: Path) -> "tuple[int, bool]":
    """Run `pyhmmer.hmmer.hmmpress` on `hmm_path` to create the
    four `.h3*` index files pyhmmer needs for fast hmmscan.

    Returns `(n_profiles, success_bool)`. `success_bool=False` is
    surfaced separately from raising because hmmpress failing is
    not fatal — pyhmmer can still scan an unpressed `.hmm`, just
    much slower. The user gets a notify either way.

    Hardened (sweep #28):
      * Pyhmmer absence treated cleanly (pressed=False, n=0).
      * Half-pressed leftovers cleaned on any failure via
        `_cleanup_pressed_files` so `_hmm_db_pressed` doesn't
        report success for a corrupted index.
      * Validates pressed files exist after hmmpress returns; if any
        of the four are missing, treats it as failure.
      * n_profiles == 0 → returns (0, False) since an "empty" press
        is functionally identical to an absent index.
    """
    n: int = 0
    if not _state._PYHMMER_AVAILABLE:
        return (0, False)
    try:
        import pyhmmer
    except ImportError:
        return (0, False)
    try:
        # pyhmmer's HMMFile yields HMMs lazily; pass the file
        # handle directly to hmmpress so the press is streaming.
        with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
            n = pyhmmer.hmmer.hmmpress(hf, str(hmm_path))
    except Exception as exc:
        _log.exception(
            "HMM DB hmmpress failed for %s (%s); will fall back to "
            "unpressed scan",
            hmm_path, exc,
        )
        _cleanup_pressed_files(hmm_path)
        return (0, False)
    # Post-press verification: all 4 index files present + n_profiles ≥ 1.
    pressed_ok = all(
        (hmm_path.parent / f"{hmm_path.name}.{ext}").exists()
        for ext in ("h3i", "h3m", "h3p", "h3f")
    )
    n_int = int(n) if isinstance(n, (int, float)) else 0
    if not pressed_ok or n_int < 1:
        _log.warning(
            "HMM DB hmmpress %s returned n=%d but index files "
            "incomplete (pressed_ok=%s); cleaning up.",
            hmm_path, n_int, pressed_ok,
        )
        _cleanup_pressed_files(hmm_path)
        return (n_int, False)
    return (n_int, True)


def _delete_hmm_db_files(entry_id: str) -> int:
    """Remove every file under `<DATA_DIR>/hmm_databases/<id>/`.
    Returns the count of files removed. Used by the catalog modal's
    "Delete download" button and by the entry-delete flow.

    Routed through `_refuse_unauthorized_delete` per [INV-75] so an
    unsandboxed probe can't wipe a real DB by guessing the id."""
    d = _hmm_db_entry_dir(entry_id)
    if not d.exists():
        return 0
    removed = 0
    for p in sorted(d.glob("*")):
        if not p.is_file():
            continue
        try:
            _refuse_unauthorized_delete(p, "hmm db file")
            p.unlink()
            removed += 1
        except (OSError, RuntimeError) as exc:
            _log.warning(
                "HMM DB delete: failed on %s (%s)", p, exc,
            )
    try:
        d.rmdir()
    except OSError:
        pass
    return removed


def _hmm_db_perform_download(
    entry: dict,
    *,
    progress_cb: "_Callable[[int, int | None], None] | None" = None,
    status_cb: "_Callable[[str], None] | None" = None,
    cancel_check_cb: "_Callable[[], bool] | None" = None,
) -> dict:
    """Download → decompress → hmmpress → stamp meta for one catalog entry.

    The HEADLESS core shared by `HmmDbCatalogModal._download_worker` (GUI)
    and the `download-hmm-database` agent endpoint, so an agent-triggered
    download lands byte-identical on disk to a user-triggered one — same
    `db.hmm(.gz)` + `.h3*` siblings, same `meta.json`, same catalog state
    ("sits nicely as if the user downloaded it themselves").

    `entry` is a catalog dict with at least ``id`` + ``url``. The optional
    callbacks are pure UI side-channels — all no-ops when None (the agent
    path passes none): ``progress_cb(done, total)`` during the byte
    stream, ``status_cb(phase)`` where *phase* is ``"decompress"`` /
    ``"hmmpress"``, and ``cancel_check_cb()`` polled between chunks
    (return True to abort).

    Returns ``{id, n_profiles, pressed, sha256, version, bytes}``. Raises
    OSError / ValueError / RuntimeError on any failure (the caller logs +
    notifies). Does NOT manage the cross-modal `_HMM_DB_DOWNLOAD_INFLIGHT`
    slot — the caller (GUI worker / agent endpoint) owns acquire+release
    so each can report "already running" in its own idiom.

    Writes go through the L2-gated `_save_hmm_db_local_meta`, so this only
    runs from a sanctioned caller (GUI app, agent server, sandboxed
    verifier); an unsandboxed probe still trips the chokepoint."""
    eid = entry["id"]
    url = entry["url"]
    redacted = _redact_url_credentials(url)
    dest_gz = _hmm_db_entry_dir(eid) / "db.hmm.gz"
    dest_hmm = _hmm_db_hmm_path(eid)

    # Phase 1: stream download (content-type + gzip/HMMER3 magic + size +
    # disk-space + redirect/scheme guards all live inside the helper).
    sha = _stream_download_to_path(
        url, dest_gz,
        max_bytes=_HMM_DB_DOWNLOAD_MAX_BYTES,
        progress_cb=progress_cb,
        cancel_check_cb=cancel_check_cb,
    )
    # Phase 2: decompress (gzip-bomb-bounded, cancel-aware).
    if status_cb:
        status_cb("decompress")
    _decompress_gz_to_path(
        dest_gz, dest_hmm,
        max_bytes=_HMM_DB_DOWNLOAD_MAX_BYTES,
        cancel_check_cb=cancel_check_cb,
    )
    # Phase 3: hmmpress. n_profiles == 0 means the bytes landed but the
    # file isn't a usable HMM DB — treat as failure (don't record success
    # for a DB hmmscan can't open).
    if status_cb:
        status_cb("hmmpress")
    n_profiles, pressed = _hmmpress_db(dest_hmm)
    if n_profiles < 1:
        raise ValueError(
            "post-press validation: 0 profiles parsed from the "
            "downloaded file. Source URL likely served the "
            "wrong file or a corrupted upload."
        )
    # Phase 4: stamp meta (mirrors the GUI worker exactly).
    existing = _load_hmm_db_local_meta(eid) or {}
    existing.update({
        "id":            eid,
        "version":       _hmm_db_local_version(eid) or "downloaded",
        "downloaded_at": _now().isoformat(timespec="seconds"),
        "sha256":        sha,
        "source_url":    redacted,
        "n_profiles":    n_profiles,
        "pressed":       pressed,
    })
    remote = (existing.get("last_remote_version") or "").strip()
    if remote:
        existing["version"] = remote
    _save_hmm_db_local_meta(eid, existing)
    return {
        "id":         eid,
        "n_profiles": n_profiles,
        "pressed":    pressed,
        "sha256":     sha,
        "version":    existing["version"],
        "bytes":      dest_gz.stat().st_size if dest_gz.exists() else 0,
    }


# ===========================================================================
# Plasmidsaurus REST API client (sweep #29) — fetch sequencing results by item
# code over the official OAuth2 client-credentials API
# (github.com/plasmidsaurus/api_docs). The downloaded `<code>_results.zip`
# flows straight into the existing `[SUB-plasmidsaurus]` zip importer in
# splicecraft_fileio (so an API-fetched run lands identically to one the user
# downloaded from their dashboard and picked by hand).
#
# Hardened like the rest of the network layer: egress-gated (the hardened
# opener calls `_state._demo_block_network_hook`), HTTPS-only with bounded
# redirects + no http downgrade, credentials redacted in logs, JSON/zip size
# capped, and the disk write goes through the L2 `_refuse_unauthorized_write`
# chokepoint. The item code is validated against the published `^[A-Z0-9]{6}$`
# shape BEFORE it's interpolated into any URL. Credentials resolve env-first
# (PLASMIDSAURUS_CLIENT_ID / _SECRET) then settings (plaintext — the Settings
# UI flags that). Kept as a DEDICATED downloader rather than generalising the
# HMM-DB `_stream_download_to_path` so the shipped HMM path carries no risk.
# ===========================================================================

_PLASMIDSAURUS_API_URL = "https://app.plasmidsaurus.com"
_PLASMIDSAURUS_ITEM_CODE_RE = re.compile(r"[A-Z0-9]{6}")
_PLASMIDSAURUS_RESULT_KINDS: frozenset[str] = frozenset({
    "results", "reads", "pod5"})
_PLASMIDSAURUS_NETWORK_TIMEOUT_S = 120         # connect + read; per request
_PLASMIDSAURUS_TOKEN_MAX_BYTES = 64 * 1024     # OAuth token JSON
_PLASMIDSAURUS_API_MAX_BYTES = 16 * 1024 * 1024   # items/samples JSON listing
# Results-zip cap a touch above the importer's 500 MB so the clearer
# "parser refused" message wins over a raw byte-cap abort for a borderline run.
_PLASMIDSAURUS_DOWNLOAD_MAX_BYTES = 600 * 1024 * 1024
# `/api/items` DEFAULTS to 100 items and silently truncates — an account with
# 218 orders looked like it only had its most recent 100 (one year of a
# four-year history) until this was found against a live account (2026-07-27).
# `?limit=N` is honoured up to 1000 ("Limit parameter can't be more than a
# 1000."); there is NO pagination — `offset` / `page` / `skip` / `per_page` /
# `cursor` are all accepted-and-ignored, returning the identical first page. So
# `limit` is the only lever, and past 1000 items the tail is unreachable via
# this API (callers surface `_PLASMIDSAURUS_ITEMS_TRUNCATED_HINT` at the cap).
_PLASMIDSAURUS_ITEMS_LIMIT = 1000
_PLASMIDSAURUS_ITEMS_TRUNCATED_HINT = (
    f"listing hit the server's {_PLASMIDSAURUS_ITEMS_LIMIT}-item ceiling; "
    f"older orders exist but the API offers no pagination to reach them")
# Live-observed rate limit, shared across `/oauth/token` + `/api/*`: exceeding
# it returns `429 {"message": "10 per 1 minute"}`. Tight enough that a single
# GUI fetch (token + link + download) is already 3 of the 10.
_PLASMIDSAURUS_RATE_LIMIT_MSG = (
    "Plasmidsaurus rate limit reached (10 requests per minute) — wait about "
    "a minute and try again.")
# Order types that never carry a results zip. Shipping labels are ~40% of a
# typical listing, are marked `status: complete` with a null `done_date`, and
# sort to the TOP — so "the newest complete item" is usually a shipping label,
# whose `/results` 404s. Used by `_plasmidsaurus_item_has_results` to flag them
# rather than hide them (dropping rows would break the agent-API compat
# promise for anything already consuming the listing).
_PLASMIDSAURUS_NON_RESULT_PRODUCTS: frozenset[str] = frozenset({
    "ups_shipping_label"})


def _plasmidsaurus_user_agent() -> str:
    return f"SpliceCraft/{_state._sc_version} (Plasmidsaurus API client)"


def _plasmidsaurus_oserror(msg: str, status: "int | None") -> OSError:
    """Build the OSError callers see, tagging it with the HTTP status so a
    caller can branch on it (`getattr(exc, "http_status", None)`) instead of
    string-sniffing the message."""
    err = OSError(msg)
    err.http_status = status        # type: ignore[attr-defined]
    return err


def _plasmidsaurus_item_has_results(item: "dict | None") -> bool:
    """Best-effort 'this order has a downloadable results zip'.

    NOT authoritative — the API exposes no such field, so this is the
    heuristic that matches live data: a real completed run is `status:
    complete` AND carries a `done_date` AND isn't one of the known
    non-sequencing products. Shipping labels fail on the last two, canceled
    orders on the first. Callers use it to steer a picker away from the
    dead rows; the download path still has to handle a 404.

    Both error directions are tolerable because the row is listed either
    way: a run that is complete but not yet dated reads False for that
    window (the user can still feed its code to the download), and the
    product list is a snapshot that a new non-sequencing order type could
    outrun — which is why `done_date` is checked too rather than trusting
    `_PLASMIDSAURUS_NON_RESULT_PRODUCTS` alone to stay current."""
    if not isinstance(item, dict):
        return False
    if str(item.get("status") or "").strip().lower() != "complete":
        return False
    if not item.get("done_date"):
        return False
    product = str(item.get("product_name") or "").strip().lower()
    return product not in _PLASMIDSAURUS_NON_RESULT_PRODUCTS


def _plasmidsaurus_credentials() -> "tuple[str | None, str | None]":
    """Resolve the OAuth Client ID + Secret env-first, then settings.

    Env vars (``PLASMIDSAURUS_CLIENT_ID`` / ``PLASMIDSAURUS_CLIENT_SECRET``)
    win — they match Plasmidsaurus's own official scripts and keep the secret
    off disk. The settings fallback (plaintext in the data dir; the Settings
    UI flags this) is the convenience path for users who'd rather not export
    env vars. Returns ``(None, None)`` for any missing half so callers refuse
    with a clear 'set credentials' message instead of attempting a request."""
    cid = (os.environ.get("PLASMIDSAURUS_CLIENT_ID") or "").strip()
    sec = (os.environ.get("PLASMIDSAURUS_CLIENT_SECRET") or "").strip()
    if not cid:
        cid = str(_get_setting("plasmidsaurus_client_id", "") or "").strip()
    if not sec:
        sec = str(_get_setting("plasmidsaurus_client_secret", "") or "").strip()
    return (cid or None, sec or None)


def _sanitize_plasmidsaurus_item_code(code: "str | None") -> "str | None":
    """Validate/normalise an item code to the published ``^[A-Z0-9]{6}$``
    shape (upper-cased). Returns None on any deviation so callers refuse it
    BEFORE it's interpolated into a URL path (defends against
    ``ABC123/../../etc`` and friends). Type-strict: non-string → None."""
    if not isinstance(code, str):
        return None
    c = code.strip().upper()
    if _PLASMIDSAURUS_ITEM_CODE_RE.fullmatch(c):
        return c
    return None


def _plasmidsaurus_oauth_token(client_id: str, client_secret: str,
                                *, timeout: "float | None" = None) -> str:
    """Redeem a short-lived Bearer token via the OAuth2 client-credentials
    grant: POST ``grant_type=client_credentials&scope=item:read`` to
    ``/oauth/token`` with HTTP Basic auth. Raises OSError on network failure
    or rejected credentials (400/401/403), ValueError on a malformed token
    body."""
    import base64
    import json
    import urllib.error
    import urllib.parse
    import urllib.request

    if not client_id or not client_secret:
        raise ValueError("missing Plasmidsaurus client_id / client_secret")
    if timeout is None:
        timeout = _PLASMIDSAURUS_NETWORK_TIMEOUT_S

    basic = base64.b64encode(
        f"{client_id}:{client_secret}".encode("utf-8")).decode("ascii")
    body = urllib.parse.urlencode({
        "grant_type": "client_credentials",
        "scope": "item:read",
    }).encode("ascii")
    opener = _build_hardened_url_opener()
    req = urllib.request.Request(
        f"{_PLASMIDSAURUS_API_URL}/oauth/token",
        data=body,
        headers={
            "Authorization": f"Basic {basic}",
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "application/json",
            "User-Agent": _plasmidsaurus_user_agent(),
        },
    )
    try:
        resp = opener.open(req, timeout=timeout)
    except urllib.error.HTTPError as exc:
        # 429 first: it shares the 4xx space but says nothing about the
        # credentials, and telling a rate-limited user to "check your
        # client ID" sends them off debugging a key that's perfectly fine.
        if exc.code == 429:
            raise _plasmidsaurus_oserror(
                _PLASMIDSAURUS_RATE_LIMIT_MSG, 429) from exc
        if exc.code in (400, 401, 403):
            raise _plasmidsaurus_oserror(
                "Plasmidsaurus rejected the API credentials "
                f"(HTTP {exc.code}). Check PLASMIDSAURUS_CLIENT_ID / "
                "PLASMIDSAURUS_CLIENT_SECRET (or the Settings values).",
                exc.code) from exc
        raise _plasmidsaurus_oserror(
            f"Plasmidsaurus token request failed: HTTP {exc.code}",
            exc.code) from exc
    except urllib.error.URLError as exc:
        raise OSError(
            f"Plasmidsaurus token request failed: {exc.reason}") from exc
    try:
        raw = resp.read(_PLASMIDSAURUS_TOKEN_MAX_BYTES + 1)
    finally:
        try:
            resp.close()
        except Exception:
            pass
    if len(raw) > _PLASMIDSAURUS_TOKEN_MAX_BYTES:
        raise ValueError("Plasmidsaurus token response too large")
    try:
        token = json.loads(raw.decode("utf-8", "replace")).get("access_token")
    except (ValueError, AttributeError) as exc:
        raise ValueError(
            "Plasmidsaurus token response wasn't valid JSON") from exc
    if not isinstance(token, str) or not token:
        raise ValueError("Plasmidsaurus token response had no access_token")
    return token


def _plasmidsaurus_api_get(path: str, token: str,
                            *, timeout: "float | None" = None,
                            max_bytes: "int | None" = None) -> "_Any":
    """GET ``{API_URL}{path}`` with the Bearer token; return parsed JSON.
    ``path`` MUST be server-constructed (callers build it from a sanitised
    item code), never raw user input. Raises OSError on network / HTTP
    failure, ValueError on an oversized or malformed JSON body."""
    import json
    import urllib.error
    import urllib.request

    if timeout is None:
        timeout = _PLASMIDSAURUS_NETWORK_TIMEOUT_S
    if max_bytes is None:
        max_bytes = _PLASMIDSAURUS_API_MAX_BYTES
    opener = _build_hardened_url_opener()
    req = urllib.request.Request(
        f"{_PLASMIDSAURUS_API_URL}{path}",
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/json",
            "User-Agent": _plasmidsaurus_user_agent(),
        },
    )
    try:
        resp = opener.open(req, timeout=timeout)
    except urllib.error.HTTPError as exc:
        if exc.code == 429:
            raise _plasmidsaurus_oserror(
                _PLASMIDSAURUS_RATE_LIMIT_MSG, 429) from exc
        if exc.code == 404:
            raise _plasmidsaurus_oserror(
                "Plasmidsaurus: item or results not found (HTTP 404)",
                404) from exc
        if exc.code in (401, 403):
            raise _plasmidsaurus_oserror(
                f"Plasmidsaurus: not authorised for this item "
                f"(HTTP {exc.code})", exc.code) from exc
        raise _plasmidsaurus_oserror(
            f"Plasmidsaurus API GET failed: HTTP {exc.code}", exc.code) from exc
    except urllib.error.URLError as exc:
        raise OSError(f"Plasmidsaurus API GET failed: {exc.reason}") from exc
    try:
        raw = resp.read(max_bytes + 1)
    finally:
        try:
            resp.close()
        except Exception:
            pass
    if len(raw) > max_bytes:
        raise ValueError("Plasmidsaurus API response too large")
    try:
        return json.loads(raw.decode("utf-8", "replace"))
    except ValueError as exc:
        raise ValueError(
            "Plasmidsaurus API response wasn't valid JSON") from exc


def _plasmidsaurus_list_items(token: str,
                               *, shared: bool = False,
                               limit: "int | None" = None,
                               timeout: "float | None" = None) -> "list":
    """Return the caller's items (most-recent first). ``shared=True`` returns
    items shared with the caller by other users. Always a list (empty if the
    account has none / the API returns a non-list).

    ALWAYS sends an explicit ``limit`` — the server's own default is 100 and
    it truncates silently, which hid 118 of a real account's 218 orders (75 of
    them actual sequencing runs) behind a listing that looked complete.
    ``limit`` is clamped to the server's 1000 ceiling; a request for more is
    rejected with a 400 rather than capped server-side."""
    if limit is None:
        limit = _PLASMIDSAURUS_ITEMS_LIMIT
    limit = max(1, min(int(limit), _PLASMIDSAURUS_ITEMS_LIMIT))
    path = (f"/api/items?shared=true&limit={limit}" if shared
            else f"/api/items?limit={limit}")
    data = _plasmidsaurus_api_get(path, token, timeout=timeout)
    return data if isinstance(data, list) else []


def _plasmidsaurus_result_link(token: str, item_code: str,
                                *, kind: str = "results",
                                timeout: "float | None" = None) -> str:
    """Return the (pre-signed, short-lived) download URL for an item's
    results / reads / pod5 zip. ``item_code`` MUST already be sanitised."""
    if kind not in _PLASMIDSAURUS_RESULT_KINDS:
        raise ValueError(f"unknown Plasmidsaurus result kind: {kind!r}")
    try:
        data = _plasmidsaurus_api_get(
            f"/api/item/{item_code}/{kind}", token, timeout=timeout)
    except OSError as exc:
        # A bare "not found" reads as a bad item code, but the common cause is
        # an order that simply has no results: shipping-label orders and
        # canceled orders both live in the listing (the former marked
        # `complete`!) and both 404 here. Name them so the user stops
        # re-checking a code that's perfectly valid.
        if getattr(exc, "http_status", None) == 404:
            raise _plasmidsaurus_oserror(
                f"Plasmidsaurus item {item_code} has no {kind} to download. "
                f"Shipping-label orders, canceled orders and runs still in "
                f"progress all appear in the item list but have no results.",
                404) from exc
        raise
    link = data.get("link") if isinstance(data, dict) else None
    if not isinstance(link, str) or not link:
        raise ValueError(
            f"Plasmidsaurus returned no {kind} download link for {item_code}")
    return link


def _plasmidsaurus_download_zip(url: str, dest: Path,
                                 *, max_bytes: int,
                                 progress_cb=None,
                                 cancel_check_cb=None,
                                 chunk_size: int = 64 * 1024) -> str:
    """Stream a Plasmidsaurus pre-signed zip URL to ``dest`` atomically;
    return the SHA-256 hex. Independently hardened (mirrors the HMM-DB
    downloader's safety machinery, kept separate so the shipped HMM path
    carries no Plasmidsaurus risk):

      * Hardened opener (HTTPS-only, bounded redirects, no http downgrade).
      * Content-Type guard rejects HTML/JSON error interstitials.
      * ZIP magic (``PK\\x03\\x04``) verified on the first chunk — fail fast
        if the link served an error page or the wrong file.
      * ``max_bytes`` cap (Content-Length pre-check + running tally) —
        zip-bomb / wrong-URL guard.
      * Disk-space pre-check (shared ``_hmm_db_check_disk_space``).
      * Atomic write: tmp → fsync → os.replace → fsync_parent.
      * L2 chokepoint via ``_refuse_unauthorized_write``.
      * ``cancel_check_cb()`` polled between chunks (return True to abort).

    Raises ValueError (cap / bad magic / bad content-type) or OSError
    (network / disk / cancellation)."""
    import hashlib
    import socket
    import urllib.error
    import urllib.request

    redacted = _redact_url_credentials(url)
    if not url.lower().startswith("https://"):
        raise ValueError("refusing non-HTTPS Plasmidsaurus download URL")

    dest.parent.mkdir(parents=True, exist_ok=True)
    _refuse_unauthorized_write(dest, "plasmidsaurus download")
    tmp = dest.with_name(dest.name + ".download_tmp")
    try:
        tmp.unlink()
    except FileNotFoundError:
        pass

    opener = _build_hardened_url_opener()
    req = urllib.request.Request(
        url, headers={"User-Agent": _plasmidsaurus_user_agent()})
    try:
        resp = opener.open(req, timeout=_PLASMIDSAURUS_NETWORK_TIMEOUT_S)
    except urllib.error.HTTPError as exc:
        # The results URL is a short-lived pre-signed S3 link, not an API
        # endpoint — the object store answers an expired or already-consumed
        # signature with 403/404, which as a bare status reads like "your
        # account can't see this run". Re-fetching the link is the fix, so
        # say that rather than sending the user at their credentials.
        if exc.code in (403, 404):
            raise _plasmidsaurus_oserror(
                f"Plasmidsaurus download link rejected (HTTP {exc.code}) — "
                f"these links are pre-signed and short-lived, so request the "
                f"download again to get a fresh one.", exc.code) from exc
        raise _plasmidsaurus_oserror(
            f"Plasmidsaurus download failed: HTTP {exc.code}",
            exc.code) from exc
    except (urllib.error.URLError, socket.timeout) as exc:
        raise OSError(f"Plasmidsaurus download failed: {exc}") from exc

    sha = hashlib.sha256()
    bytes_so_far = 0
    cancelled = False
    saw_bytes = False
    try:
        try:
            _hmm_db_assert_content_type_ok(resp, url)
            total_header = resp.headers.get("Content-Length")
            total: "int | None"
            try:
                total = int(total_header) if total_header else None
            except (TypeError, ValueError):
                total = None
            if total is not None and total > max_bytes:
                raise ValueError(
                    f"refusing download: server reports {total:,} bytes "
                    f"> cap {max_bytes:,}")
            _hmm_db_check_disk_space(dest, total)
            first_chunk = True
            with open(tmp, "wb") as fh:
                while True:
                    if (cancel_check_cb is not None
                            and bool(cancel_check_cb())):
                        cancelled = True
                        raise OSError("download cancelled by user")
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    bytes_so_far += len(chunk)
                    if bytes_so_far > max_bytes:
                        raise ValueError(
                            f"download exceeded {max_bytes:,} byte cap "
                            f"({bytes_so_far:,} read) — likely wrong URL "
                            f"or zip bomb")
                    if first_chunk:
                        first_chunk = False
                        if not chunk.startswith(b"PK\x03\x04"):
                            raise ValueError(
                                "first bytes aren't a zip (PK 03 04) — the "
                                "link didn't serve a results archive")
                        saw_bytes = True
                    sha.update(chunk)
                    fh.write(chunk)
                    if progress_cb is not None:
                        try:
                            progress_cb(bytes_so_far, total)
                        except Exception:
                            _log.exception(
                                "Plasmidsaurus download: progress cb raised")
                fh.flush()
                try:
                    os.fsync(fh.fileno())
                except OSError:
                    pass
        finally:
            try:
                resp.close()
            except Exception:
                pass
        if not saw_bytes:
            raise ValueError("download empty — server returned 0 bytes")
        os.replace(str(tmp), str(dest))
        _fsync_parent_dir(dest)
    except Exception:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass
        if cancelled:
            _log_event("plasmidsaurus.download.cancelled",
                       url=redacted, bytes_read=bytes_so_far)
        raise
    _log_event("plasmidsaurus.download.ok", url=redacted, bytes=bytes_so_far)
    return sha.hexdigest()


def _plasmidsaurus_fetch_item_zip(item_code: str, dest_dir: "str | Path",
                                   *, kind: str = "results",
                                   client_id: "str | None" = None,
                                   client_secret: "str | None" = None,
                                   token: "str | None" = None,
                                   progress_cb=None,
                                   cancel_check_cb=None) -> Path:
    """End-to-end: resolve creds → token → results link → stream the zip to
    ``dest_dir/<code>_<kind>.zip``. Returns the written path. The HEADLESS
    core shared by the agent endpoint and the GUI fetch worker, so both land
    the same bytes on disk. Raises ValueError (bad item code / kind / cred
    gap) or OSError (network / disk)."""
    code = _sanitize_plasmidsaurus_item_code(item_code)
    if code is None:
        raise ValueError("invalid item code (expected 6 chars, A-Z / 0-9)")
    if kind not in _PLASMIDSAURUS_RESULT_KINDS:
        raise ValueError(f"unknown result kind: {kind!r}")
    if token is None:
        if not client_id or not client_secret:
            cid, sec = _plasmidsaurus_credentials()
            client_id = client_id or cid
            client_secret = client_secret or sec
        if not client_id or not client_secret:
            raise ValueError(
                "no Plasmidsaurus API credentials — set "
                "PLASMIDSAURUS_CLIENT_ID / PLASMIDSAURUS_CLIENT_SECRET "
                "or add them in Settings")
        token = _plasmidsaurus_oauth_token(client_id, client_secret)
    link = _plasmidsaurus_result_link(token, code, kind=kind)
    dest = Path(dest_dir) / f"{code}_{kind}.zip"
    _plasmidsaurus_download_zip(
        link, dest, max_bytes=_PLASMIDSAURUS_DOWNLOAD_MAX_BYTES,
        progress_cb=progress_cb, cancel_check_cb=cancel_check_cb)
    return dest


# ── Online reference-database lookups (Babs / agent-API `*-search`) ─────────
# DISTINCT from the BLAST/HMMER egress path above: these send only a short
# QUERY STRING (a name / search term) to public databases — never the user's
# sequence. Armed by the SEPARATE `allow_online_lookups` setting (see the
# agent-side `_online_lookups_armed`). Every fetch still routes through the
# shared hardened opener (`_online_http` → `_build_hardened_url_opener`): SSRF
# host filter, bounded redirects, https-only, and the fail-closed DEMO_MODE
# egress gate. Never log query/response content ([INV-38]).
_ONLINE_LOOKUP_MAX_HITS = 25              # polite hit cap per public service
_ONLINE_LOOKUP_TIMEOUT_S = _NCBI_TIMEOUT_S
_ONLINE_LOOKUP_MAX_RESPONSE_BYTES = _NCBI_MAX_RESPONSE_BYTES   # 4 MB
_ONLINE_LOOKUP_QUERY_MAX = 400            # reject absurd query strings
# Google Patents' XHR and DuckDuckGo's HTML endpoint refuse a non-browser UA
# (they answer a "Sorry…" block page instead). Used ONLY for those keyless
# best-effort scrapers; the official JSON APIs get the honest SpliceCraft UA.
_ONLINE_LOOKUP_BROWSER_UA = (
    "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0"
)
# FPbase stores oligomerization + switching as short codes; map the common
# ones to readable words, pass an unknown code through unchanged.
_FPBASE_AGG = {"m": "monomer", "d": "dimer", "td": "tandem dimer",
               "wd": "weak dimer", "t": "tetramer"}
_FPBASE_SWITCH = {"b": "basic", "pa": "photoactivatable",
                  "ps": "photoswitchable", "pc": "photoconvertible",
                  "tf": "timer", "o": "other", "c": "chromoprotein"}


def _lookup_clean_html(s: "_Any") -> str:
    """Strip HTML tags + unescape entities from a snippet/abstract fragment
    (search APIs and scrapers return `<b>`-marked-up text). Non-string → ''."""
    import html as _html
    if not isinstance(s, str):
        return ""
    return _html.unescape(re.sub(r"<[^>]+>", "", s)).strip()


def _online_lookup_json(url: str, *, headers: "dict | None" = None,
                        browser_ua: bool = False) -> "_Any":
    """GET `url` through the hardened opener and parse JSON. `browser_ua`
    swaps in a browser User-Agent for the endpoints that refuse ours. Raises
    RuntimeError on a network/HTTP failure or a non-JSON body — the agent
    handlers map that to a 502."""
    hdrs: "dict[str, str]" = {"Accept": "application/json"}
    if browser_ua:
        hdrs["User-Agent"] = _ONLINE_LOOKUP_BROWSER_UA
    if headers:
        hdrs.update(headers)
    body = _online_http(url, headers=hdrs, timeout=_ONLINE_LOOKUP_TIMEOUT_S,
                        max_bytes=_ONLINE_LOOKUP_MAX_RESPONSE_BYTES)
    try:
        return json.loads(body)
    except ValueError as exc:
        raise RuntimeError(f"non-JSON response from server: {exc}")


def _fpbase_search(query: str,
                   max_hits: int = _ONLINE_LOOKUP_MAX_HITS) -> "list[dict]":
    """FPbase (fpbase.org) fluorescent-protein lookup by name (substring).
    Returns compact records: spectra (ex/em maxima, QY, extinction,
    brightness), oligomerization, switch type, GenBank/UniProt xrefs, the
    aa sequence, and the FPbase URL."""
    import urllib.parse as _up
    qs = _up.urlencode({"format": "json", "name__icontains": query})
    data = _online_lookup_json(f"https://www.fpbase.org/api/proteins/?{qs}")
    rows = data.get("results", data) if isinstance(data, dict) else data
    if not isinstance(rows, list):
        return []
    out: "list[dict]" = []
    for p in rows[:max_hits]:
        if not isinstance(p, dict):
            continue
        name = (p.get("name") or "").strip()
        slug = (p.get("slug") or "").strip()
        if not name:
            continue
        states = p.get("states") or []
        default = next((s for s in states if isinstance(s, dict)
                        and s.get("slug") == p.get("default_state")), None)
        if default is None:
            default = states[0] if states and isinstance(states[0], dict) else {}
        agg = p.get("agg")
        switch = p.get("switch_type")
        out.append({
            "name": name,
            "ex_max": default.get("ex_max"),
            "em_max": default.get("em_max"),
            "quantum_yield": default.get("qy"),
            "ext_coeff": default.get("ext_coeff"),
            "brightness": default.get("brightness"),
            "oligomerization": (_FPBASE_AGG.get(agg, agg)
                                if isinstance(agg, str) else agg),
            "switch_type": (_FPBASE_SWITCH.get(switch, switch)
                            if isinstance(switch, str) else switch),
            "genbank": p.get("genbank") or None,
            "uniprot": p.get("uniprot") or None,
            "seq": p.get("seq") or None,
            "url": (f"https://www.fpbase.org/protein/{slug}/"
                    if slug else "https://www.fpbase.org/"),
        })
    return out


def _uniprot_name(entry: dict) -> str:
    desc = entry.get("proteinDescription") or {}
    for slot in ("recommendedName", "submissionNames"):
        v = desc.get(slot)
        if isinstance(v, list):
            v = v[0] if v else None
        if isinstance(v, dict):
            fn = (v.get("fullName") or {}).get("value")
            if fn:
                return fn
    return entry.get("uniProtkbId") or entry.get("primaryAccession") or "protein"


def _uniprot_function(entry: dict) -> "str | None":
    for c in (entry.get("comments") or []):
        if isinstance(c, dict) and c.get("commentType") == "FUNCTION":
            vals = [str(t.get("value")) for t in (c.get("texts") or [])
                    if isinstance(t, dict) and t.get("value")]
            if vals:
                s = " ".join(vals)
                return s[:600] + ("…" if len(s) > 600 else "")
    return None


def _uniprot_search(query: str,
                    max_hits: int = _ONLINE_LOOKUP_MAX_HITS) -> "list[dict]":
    """UniProtKB protein lookup. Returns accession, protein name, organism,
    the curated FUNCTION comment, keywords, and reviewed (Swiss-Prot) flag."""
    import urllib.parse as _up
    qs = _up.urlencode({
        "query": query, "format": "json", "size": max_hits,
        "fields": "accession,id,protein_name,organism_name,cc_function,keyword",
    })
    data = _online_lookup_json(f"https://rest.uniprot.org/uniprotkb/search?{qs}")
    results = data.get("results") if isinstance(data, dict) else None
    if not isinstance(results, list):
        return []
    out: "list[dict]" = []
    for e in results[:max_hits]:
        if not isinstance(e, dict):
            continue
        acc = (e.get("primaryAccession") or "").strip()
        et = e.get("entryType") or ""
        out.append({
            "accession": acc,
            "name": _uniprot_name(e),
            "organism": (e.get("organism") or {}).get("scientificName"),
            "function": _uniprot_function(e),
            "keywords": [k.get("name") for k in (e.get("keywords") or [])
                         if isinstance(k, dict) and k.get("name")][:12],
            "reviewed": "Swiss-Prot" in et,
            "url": (f"https://www.uniprot.org/uniprotkb/{acc}"
                    if acc else "https://www.uniprot.org/"),
        })
    return out


def _europepmc_search(query: str,
                      max_hits: int = _ONLINE_LOOKUP_MAX_HITS) -> "list[dict]":
    """Europe PMC literature search (PubMed + PMC + Agricola + preprints).
    Returns title, authors, journal, year, DOI/PMID/PMCID, a truncated
    abstract, open-access flag, and a resolvable URL."""
    import urllib.parse as _up
    qs = _up.urlencode({"query": query, "format": "json",
                        "pageSize": max_hits, "resultType": "core"})
    data = _online_lookup_json(
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{qs}")
    results = ((data.get("resultList") or {}).get("result")
               if isinstance(data, dict) else None)
    if not isinstance(results, list):
        return []
    out: "list[dict]" = []
    for r in results[:max_hits]:
        if not isinstance(r, dict):
            continue
        doi, pmid, pmcid = r.get("doi"), r.get("pmid"), r.get("pmcid")
        if doi:
            url = f"https://doi.org/{doi}"
        elif pmid:
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        elif pmcid:
            url = f"https://europepmc.org/article/PMC/{pmcid}"
        else:
            url = "https://europepmc.org/"
        abstract = r.get("abstractText") or ""
        out.append({
            "title": _lookup_clean_html(r.get("title")),
            "authors": r.get("authorString"),
            "journal": r.get("journalTitle") or r.get("source"),
            "year": r.get("pubYear"),
            "doi": doi, "pmid": pmid, "pmcid": pmcid,
            "abstract": ((_lookup_clean_html(abstract)[:800]
                          + ("…" if len(abstract) > 800 else ""))
                         if abstract else None),
            "is_open_access": r.get("isOpenAccess") == "Y",
            "url": url,
        })
    return out


def _ncbi_db_search(query: str, db: str = "nucleotide",
                    max_hits: int = 10) -> "tuple[list[dict], int | None]":
    """NCBI Entrez term search over `nucleotide` or `protein` (esearch →
    esummary). Returns (records, total_matches); each record carries
    accession, title, organism, length, and an ncbi.nlm.nih.gov URL — the
    accession feeds the existing `fetch` endpoint to pull the full record."""
    import urllib.parse as _up
    eutils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    es = eutils + "esearch.fcgi?" + _up.urlencode(
        {"db": db, "term": query, "retmode": "json", "retmax": max_hits})
    data = _online_lookup_json(es)
    er = data.get("esearchresult") or {} if isinstance(data, dict) else {}
    ids = [i for i in (er.get("idlist") or []) if isinstance(i, str)]
    cnt = er.get("count")
    total: "int | None" = (int(cnt) if isinstance(cnt, str)
                           and cnt.isdigit() else None)
    if not ids:
        return [], total
    su = eutils + "esummary.fcgi?" + _up.urlencode(
        {"db": db, "id": ",".join(ids), "retmode": "json"})
    sdata = _online_lookup_json(su)
    res = (sdata.get("result") or {}) if isinstance(sdata, dict) else {}
    is_prot = db == "protein"
    out: "list[dict]" = []
    for uid in (res.get("uids") or ids):
        e = res.get(uid)
        if not isinstance(e, dict):
            continue
        acc = (e.get("accessionversion") or e.get("caption") or "").strip()
        out.append({
            "accession": acc,
            "title": (e.get("title") or "").strip(),
            "organism": e.get("organism"),
            "length": e.get("slen"),
            "moltype": e.get("moltype") or e.get("biomol"),
            "update_date": e.get("updatedate"),
            "url": (f"https://www.ncbi.nlm.nih.gov/"
                    f"{'protein' if is_prot else 'nuccore'}/{acc}"
                    if acc else "https://www.ncbi.nlm.nih.gov/"),
        })
    return out, total


def _wikipedia_search(query: str,
                      max_hits: int = _ONLINE_LOOKUP_MAX_HITS) -> "list[dict]":
    """Wikipedia (MediaWiki) search. One search call for the hit list, then a
    bounded second call for the lead-section plaintext of the top results;
    hits beyond that carry the search snippet. Returns title, summary, URL."""
    import urllib.parse as _up
    api = "https://en.wikipedia.org/w/api.php?"
    s = api + _up.urlencode({
        "action": "query", "list": "search", "srsearch": query,
        "srlimit": max_hits, "srprop": "snippet", "format": "json"})
    data = _online_lookup_json(s)
    hits = ((data.get("query") or {}).get("search")
            if isinstance(data, dict) else None)
    if not isinstance(hits, list) or not hits:
        return []
    top = [str(h.get("pageid")) for h in hits[:5]
           if isinstance(h, dict) and h.get("pageid")]
    extracts: "dict[str, str]" = {}
    if top:
        ex = api + _up.urlencode({
            "action": "query", "prop": "extracts", "exintro": 1,
            "explaintext": 1, "redirects": 1, "format": "json",
            "pageids": "|".join(top)})
        try:
            pages = (((_online_lookup_json(ex).get("query") or {})
                      .get("pages")) or {})
            for pid, pg in pages.items():
                if isinstance(pg, dict) and pg.get("extract"):
                    extracts[str(pg.get("pageid") or pid)] = pg["extract"]
        except RuntimeError:
            pass    # extracts are best-effort enrichment; snippets still ship
    out: "list[dict]" = []
    for h in hits[:max_hits]:
        if not isinstance(h, dict):
            continue
        pid = h.get("pageid")
        ext = extracts.get(str(pid))
        summary = ((ext[:600] + ("…" if len(ext) > 600 else "")) if ext
                   else _lookup_clean_html(h.get("snippet")))
        out.append({
            "title": (h.get("title") or "").strip(),
            "summary": summary,
            "url": (f"https://en.wikipedia.org/?curid={pid}"
                    if pid else "https://en.wikipedia.org/"),
        })
    return out


def _ddg_parse_html(html_text: str) -> "list[dict]":
    """Parse the DuckDuckGo HTML-endpoint result page into {title, url,
    snippet} dicts. Raises RuntimeError if DDG served a block/anomaly page
    (keyless search is rate-limited)."""
    import urllib.parse as _up
    if "result__a" not in html_text:
        low = html_text.lower()
        if ("anomaly" in low or "blocked" in low
                or "<title>sorry" in low[:600]):
            raise RuntimeError(
                "DuckDuckGo blocked the request (keyless search is "
                "rate-limited). Try again shortly, or set a Brave Search "
                "API key in Settings for robust web search.")
        return []
    anchors = re.findall(
        r'result__a"[^>]*href="([^"]+)"[^>]*>(.*?)</a>', html_text, re.DOTALL)
    snippets = re.findall(
        r'result__snippet"[^>]*>(.*?)</a>', html_text, re.DOTALL)
    out: "list[dict]" = []
    for i, (href, title) in enumerate(anchors):
        url = href
        if "uddg=" in href:
            try:
                probe = ("https:" + href) if href.startswith("//") else href
                q = _up.parse_qs(_up.urlsplit(probe).query)
                if q.get("uddg"):
                    url = q["uddg"][0]
            except Exception:
                pass
        out.append({
            "title": _lookup_clean_html(title),
            "url": url,
            "snippet": _lookup_clean_html(snippets[i]) if i < len(snippets) else "",
        })
    return out


def _web_search(query: str, max_hits: int = 10,
                brave_key: str = "") -> "list[dict]":
    """General web search. With a Brave Search API key → the official Brave
    API (robust). Without → DuckDuckGo's HTML endpoint (keyless, best-effort,
    rate-limited). Returns {title, url, snippet}."""
    import urllib.parse as _up
    if brave_key:
        qs = _up.urlencode({"q": query, "count": max_hits})
        data = _online_lookup_json(
            f"https://api.search.brave.com/res/v1/web/search?{qs}",
            headers={"X-Subscription-Token": brave_key})
        results = ((data.get("web") or {}).get("results")
                   if isinstance(data, dict) else None) or []
        out: "list[dict]" = []
        for r in results[:max_hits]:
            if not isinstance(r, dict):
                continue
            out.append({"title": _lookup_clean_html(r.get("title")),
                        "url": r.get("url"),
                        "snippet": _lookup_clean_html(r.get("description"))})
        return out
    body = _online_http(
        "https://html.duckduckgo.com/html/?" + _up.urlencode({"q": query}),
        headers={"User-Agent": _ONLINE_LOOKUP_BROWSER_UA},
        timeout=_ONLINE_LOOKUP_TIMEOUT_S,
        max_bytes=_ONLINE_LOOKUP_MAX_RESPONSE_BYTES)
    return _ddg_parse_html(body)[:max_hits]


# ── read-url: fetch ONE page → readable text ─────────────────────────────────
# Babs' "open a search result and actually read it" primitive (pairs with
# web-search). SSRF-hardened through the same `_build_hardened_url_opener` as
# every other egress path — public http(s) hosts only, on the initial URL AND
# every redirect hop. It is a READER, not a browser: no JavaScript is run.
_READ_URL_MAX_URL_LEN = 2048
_READ_URL_DEFAULT_CHARS = 20000       # default returned-text budget
_READ_URL_MAX_CHARS = 100000          # hard ceiling (keep _h_read_url doc in sync)
# Content types we can turn into text. Anything else (PDF, images, octet-stream)
# is refused rather than decoded into mojibake; `text/*` is allowed broadly.
_READ_URL_TEXT_TYPES = frozenset({
    "text/html", "application/xhtml+xml", "application/xml", "text/xml",
    "text/plain",
})
# Tags whose TEXT content is dropped entirely (scripts/styles/non-content).
_HTML_DROP_TAGS = frozenset({
    "script", "style", "noscript", "template", "svg", "head", "iframe",
    "object", "canvas", "math",
})
# Tags that imply a line break around their content (readability).
_HTML_BLOCK_TAGS = frozenset({
    "p", "div", "br", "li", "tr", "section", "article", "header", "footer",
    "nav", "aside", "main", "h1", "h2", "h3", "h4", "h5", "h6", "ul", "ol",
    "table", "thead", "tbody", "blockquote", "pre", "hr", "figure",
    "figcaption", "dd", "dt", "dl", "form",
})


def _html_to_text(html_text: str) -> "tuple[str, str]":
    """Extract ``(title, readable_text)`` from an HTML document using ONLY the
    stdlib parser — script/style/head text dropped, block tags become line
    breaks, whitespace collapsed, entities decoded. NO JavaScript is executed
    (this is a reader, not a browser). Pure + deterministic → unit-tested."""
    from html.parser import HTMLParser

    class _Extractor(HTMLParser):
        def __init__(self):
            super().__init__(convert_charrefs=True)
            self.parts: "list[str]" = []
            self.title_parts: "list[str]" = []
            self._drop = 0
            self._in_title = False

        def handle_starttag(self, tag, attrs):
            if tag in _HTML_DROP_TAGS:
                self._drop += 1
            elif tag == "title":
                self._in_title = True
            elif tag in _HTML_BLOCK_TAGS:
                self.parts.append("\n")

        def handle_startendtag(self, tag, attrs):
            if tag in _HTML_BLOCK_TAGS:
                self.parts.append("\n")

        def handle_endtag(self, tag):
            if tag in _HTML_DROP_TAGS:
                self._drop = max(0, self._drop - 1)
            elif tag == "title":
                self._in_title = False
            elif tag in _HTML_BLOCK_TAGS:
                self.parts.append("\n")

        def handle_data(self, data):
            # Title first: it lives inside <head> (a drop tag) but we want it.
            if self._in_title:
                self.title_parts.append(data)
                return
            if self._drop:
                return
            # Collapse whitespace WITHIN the text node (incl. source newlines,
            # which are just whitespace in HTML) so only the "\n" markers we
            # inject for block tags become real line breaks.
            self.parts.append(re.sub(r"\s+", " ", data))

    p = _Extractor()
    try:
        p.feed(html_text)
        p.close()
    except Exception:
        # A malformed document must not sink the fetch — keep what parsed.
        pass
    title = re.sub(r"\s+", " ", "".join(p.title_parts)).strip()
    out_lines: "list[str]" = []
    blank = False
    for ln in "".join(p.parts).split("\n"):
        ln = re.sub(r"[^\S\n]+", " ", ln).strip()   # collapse intra-line spaces
        if ln:
            out_lines.append(ln)
            blank = False
        elif not blank:
            out_lines.append("")                    # one blank between paragraphs
            blank = True
    return title, "\n".join(out_lines).strip()


def _maybe_decompress(raw: bytes, encoding: str, cap: int) -> bytes:
    """Decompress a gzip / deflate response body. Unlike the fixed-domain API
    lookups, `read-url` fetches arbitrary pages and a CDN may compress the body
    even though we never advertise ``Accept-Encoding`` — so decode it here or
    the charset step would return mojibake. Output is bounded to 4× the byte cap
    so a decompression bomb can't blow memory; anything we can't decode
    (brotli / unknown / corrupt) is passed through unchanged for the
    replace-decode to handle."""
    if not encoding or encoding == "identity":
        return raw
    import zlib
    limit = cap * 4
    try:
        if encoding in ("gzip", "x-gzip"):
            return zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw, limit)
        if encoding == "deflate":                      # zlib-wrapped or raw
            try:
                return zlib.decompressobj().decompress(raw, limit)
            except zlib.error:
                return zlib.decompressobj(-zlib.MAX_WBITS).decompress(raw, limit)
    except (zlib.error, OSError, EOFError):
        return raw
    return raw                                          # br / unknown → passthrough


def _fetch_page_raw(url: str) -> "tuple[str, str, str]":
    """GET one http(s) ``url`` through the SSRF-hardened opener and return
    ``(final_url, content_type, decoded_body)``. Split out from `_read_url` so
    the extraction layer is unit-testable without network. Enforces the shared
    response-size cap, decompresses gzip/deflate, and decodes per the response
    charset (utf-8 fallback). Raises RuntimeError on any network / HTTP failure
    (incl. a refused non-public host — the opener blocks SSRF before the socket
    opens)."""
    import socket
    import urllib.error
    import urllib.request
    cap = _ONLINE_LOOKUP_MAX_RESPONSE_BYTES
    req = urllib.request.Request(url, headers={
        "User-Agent": _ONLINE_LOOKUP_BROWSER_UA,
        "Accept": "text/html,application/xhtml+xml,text/plain;q=0.9,*/*;q=0.5",
    })
    opener = _build_hardened_url_opener()
    try:
        with opener.open(req, timeout=_ONLINE_LOOKUP_TIMEOUT_S) as resp:
            ctype = (resp.headers.get("Content-Type") or "").split(
                ";", 1)[0].strip().lower()
            enc = (resp.headers.get("Content-Encoding") or "").strip().lower()
            raw = resp.read(cap + 1)
            charset = resp.headers.get_content_charset() or "utf-8"
            final_url = resp.geturl() or url
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code} fetching the page") from exc
    except urllib.error.URLError as exc:
        reason = getattr(exc, "reason", exc)
        if isinstance(reason, socket.timeout):
            raise RuntimeError("Connection timed out fetching the page.") from exc
        raise RuntimeError(f"cannot fetch the page: {reason}") from exc
    except socket.timeout as exc:
        raise RuntimeError("Connection timed out fetching the page.") from exc
    if len(raw) > cap:
        raise RuntimeError(
            f"page exceeded the {cap // (1024 * 1024)} MB fetch cap — "
            f"refusing it.")
    raw = _maybe_decompress(raw, enc, cap)
    try:
        body = raw.decode(charset, "replace")
    except (LookupError, TypeError):
        body = raw.decode("utf-8", "replace")
    return final_url, ctype, body


def _read_url(url: str, max_chars: int = _READ_URL_DEFAULT_CHARS) -> dict:
    """Fetch ONE http(s) page and return its readable text (HTML → text; NO
    JavaScript). Normalises a bare host to https, rejects non-http(s) schemes,
    SSRF-hardened via `_fetch_page_raw`, refuses binary content types, and caps
    the returned text at `max_chars`. Returns
    ``{url, title, content_type, text, chars, truncated}``. Raises RuntimeError
    on a network/HTTP failure, a bad scheme, or an unsupported content type."""
    if not isinstance(url, str) or not url.strip():
        raise RuntimeError("missing URL")
    url = url.strip()
    if url.startswith("//"):
        url = "https:" + url                         # protocol-relative → https
    m = re.match(r"(?i)^([a-z][a-z0-9+.\-]*):", url)
    if m:
        if m.group(1).lower() not in ("http", "https"):
            raise RuntimeError(
                f"unsupported URL scheme {m.group(1).lower()!r} — read-url "
                f"fetches http/https pages only.")
    else:
        url = "https://" + url                       # bare host/path → assume https
    if len(url) > _READ_URL_MAX_URL_LEN:
        raise RuntimeError(f"URL too long (max {_READ_URL_MAX_URL_LEN} chars).")
    try:
        max_chars = int(max_chars)
    except (TypeError, ValueError):
        max_chars = _READ_URL_DEFAULT_CHARS
    max_chars = max(500, min(max_chars, _READ_URL_MAX_CHARS))

    final_url, ctype, body = _fetch_page_raw(url)
    if ctype and ctype not in _READ_URL_TEXT_TYPES and not ctype.startswith(
            "text/"):
        raise RuntimeError(
            f"unsupported content type {ctype!r} — read-url returns page TEXT "
            f"(HTML / plain text), not binary files (PDF, images, downloads).")
    if ctype == "text/plain":
        title, text = "", body
    else:
        title, text = _html_to_text(body)
    text = text.strip()
    truncated = len(text) > max_chars
    if truncated:
        text = text[:max_chars].rstrip() + "\n…[truncated]"
    return {"url": final_url, "title": title, "content_type": ctype or None,
            "text": text, "chars": len(text), "truncated": truncated}


def _google_patents_parse(data: "_Any",
                          max_hits: int) -> "tuple[list[dict], int | None]":
    res = (data or {}).get("results") or {} if isinstance(data, dict) else {}
    total = res.get("total_num_results")
    out: "list[dict]" = []
    for cl in (res.get("cluster") or []):
        for item in (cl.get("result") or []):
            if not isinstance(item, dict):
                continue
            p = item.get("patent") or {}
            num = (p.get("publication_number") or "").strip()
            pid = (item.get("id") or "").strip()          # "patent/XX/en"
            if not num and pid.startswith("patent/"):
                parts = pid.split("/")
                num = parts[1] if len(parts) > 1 else ""
            if not num:
                continue
            out.append({
                "number": num,
                "title": _lookup_clean_html(p.get("title")),
                "snippet": _lookup_clean_html(p.get("snippet")),
                "assignee": _lookup_clean_html(p.get("assignee")) or None,
                "inventor": _lookup_clean_html(p.get("inventor")) or None,
                "priority_date": p.get("priority_date") or None,
                "grant_date": p.get("grant_date") or None,
                "filing_date": p.get("filing_date") or None,
                "url": f"https://patents.google.com/patent/{num}/en",
            })
            if len(out) >= max_hits:
                return out, total
    return out, total


def _patentsview_search(query: str, max_hits: int,
                        key: str) -> "tuple[list[dict], int | None]":
    import urllib.parse as _up
    qs = _up.urlencode({
        "q": json.dumps({"_text_any": {"patent_title": query}}),
        "f": json.dumps(["patent_id", "patent_title", "patent_date",
                         "patent_abstract", "assignees.assignee_organization"]),
        "o": json.dumps({"size": max_hits}),
    })
    data = _online_lookup_json(
        f"https://search.patentsview.org/api/v1/patent/?{qs}",
        headers={"X-Api-Key": key})
    patents = (data.get("patents")
               if isinstance(data, dict) else None) or []
    total = (data.get("total_hits") or data.get("count")
             if isinstance(data, dict) else None)
    out: "list[dict]" = []
    for p in patents[:max_hits]:
        if not isinstance(p, dict):
            continue
        num = (p.get("patent_id") or "").strip()
        assignees = p.get("assignees") or []
        org = None
        if assignees and isinstance(assignees[0], dict):
            org = assignees[0].get("assignee_organization")
        abstract = p.get("patent_abstract") or ""
        out.append({
            "number": num,
            "title": _lookup_clean_html(p.get("patent_title")),
            "snippet": ((_lookup_clean_html(abstract)[:600]
                         + ("…" if len(abstract) > 600 else ""))
                        if abstract else ""),
            "assignee": org,
            "grant_date": p.get("patent_date") or None,
            "url": (f"https://patents.google.com/patent/US{num}/en"
                    if num else "https://patentsview.org/"),
        })
    return out, total


def _patent_search(query: str, max_hits: int = 10,
                   patentsview_key: str = "") -> "tuple[list[dict], int | None]":
    """Patent search. With a PatentsView API key → the official PatentsView
    API (robust, US patents). Without → Google Patents' XHR endpoint (keyless,
    best-effort, rate-limited). Returns (records, total_matches)."""
    import urllib.parse as _up
    if patentsview_key:
        return _patentsview_search(query, max_hits, patentsview_key)
    url = ("https://patents.google.com/xhr/query?url=q%3D"
           + _up.quote(query) + "&exp=")
    body = _online_http(url, headers={"User-Agent": _ONLINE_LOOKUP_BROWSER_UA},
                        timeout=_ONLINE_LOOKUP_TIMEOUT_S,
                        max_bytes=_ONLINE_LOOKUP_MAX_RESPONSE_BYTES)
    try:
        data = json.loads(body)
    except ValueError:
        raise RuntimeError(
            "Google Patents blocked the request (keyless search is "
            "rate-limited). Try again shortly, or set a PatentsView API key "
            "in Settings for robust patent search.")
    return _google_patents_parse(data, max_hits)
