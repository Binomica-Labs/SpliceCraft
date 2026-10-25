"""Logging primitives for SpliceCraft (layer 0).

Extracted from the monolith so any sibling can log. Holds the `splicecraft`
logger (`_log`), the surrogate/session filters, the structured-event emitter
`_log_event`, and `_repr_for_log`. Imports nothing from the rest of the package.

DELIBERATELY LEFT IN THE HUB (do not move here):
  * the `_DATA_DIR`-dependent file-handler wiring (`_default_log_path` /
    `_LOG_PATH` + the `RotatingFileHandler` setup) -- it reads the data dir;
    the hub configures THIS logger after importing it.
  * the `_action_log` / `_timed` decorators -- they call `_log_event`
    internally, and tests monkeypatch `sc._log_event`; keeping the decorators
    in the hub makes that call resolve the patchable re-export rather than a
    private copy here. See tests/test_logging.py.
"""
from __future__ import annotations

import functools
import json
import logging
import os
import time as _time
import re
import uuid as _uuid


_SESSION_ID = _uuid.uuid4().hex[:8]

class _SessionFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        record.session = _SESSION_ID
        return True


# Control bytes (C0 incl. ESC, DEL, the C1 block incl. 8-bit CSI 0x9b)
# in a logged value would let a crafted name / filename / accession
# smuggle a terminal-escape sequence into the log file — which then
# fires when the user `cat`s the log or pastes it into a bug report
# (or hands it to an LLM scanning it). Escape them to a visible `\xNN`
# form on the logger so no callsite can leak one regardless of how it
# logs. Cheap: the search short-circuits when no control byte is present.
_LOG_CTRL_RE = re.compile(r"[\x00-\x1f\x7f-\x9f]")


class _SurrogateScrubFilter(logging.Filter):
    """A lone surrogate (e.g. what `os.fsdecode` yields for an undecodable
    filename under a non-UTF-8 locale) cannot be UTF-8-encoded. Beyond
    crashing a UTF-8 stream, it crashes **pytest-xdist's** worker→controller
    report serialization — `execnet` raises `UnicodeEncodeError` and the whole
    CI run `INTERNALERROR`s (pytest 9.1 serialises the captured record, which
    the rotating handler's `errors=` escaping happens too late to save). Scrub
    every record's message + args to a safe backslash-escaped form HERE, on the
    logger, so no surrogate can leave this logger regardless of handler config
    or test-capture path. Also escapes control bytes (incl. ESC) so a logged
    name can't smuggle a terminal-escape into the log file. Sweep 2026-06-14."""

    @staticmethod
    def _scrub(value: object) -> object:
        if isinstance(value, str):
            try:
                value.encode("utf-8")
            except UnicodeEncodeError:
                # `backslashreplace` escapes ONLY what UTF-8 cannot encode
                # (the lone surrogates); everything else stays UTF-8, so the
                # result is decoded as UTF-8. Decoding as ASCII raised on any
                # other non-ASCII character in the same string (`Gel é` +
                # a surrogate) — out of the filter and out of the `_log` call,
                # turning a SUCCESSFUL agent write into a 500 whose retry made
                # a duplicate (audit 2026-10-23).
                value = value.encode("utf-8", "backslashreplace").decode("utf-8")
            if _LOG_CTRL_RE.search(value):
                value = _LOG_CTRL_RE.sub(
                    lambda m: "\\x%02x" % ord(m.group()), value)
        return value

    def filter(self, record: logging.LogRecord) -> bool:
        # A logging filter must never raise: the exception would surface at
        # whatever call site happened to log — usually after the real work
        # was done.
        try:
            record.msg = self._scrub(record.msg)
            if isinstance(record.args, tuple):
                record.args = tuple(self._scrub(a) for a in record.args)
            elif isinstance(record.args, dict):
                record.args = {k: self._scrub(v)
                               for k, v in record.args.items()}
        except Exception:                              # noqa: BLE001
            record.msg, record.args = ascii(record.msg), ()
        return True

_log = logging.getLogger("splicecraft")
# DEBUG-level logs are off by default — bump via SPLICECRAFT_DEBUG=1
# so a user reproducing a bug can capture cache-invalidation / retry /
# pre-flight events that normally don't surface. Network retry events,
# cache-clear events, and other diagnostic-only signals key off this.
if os.environ.get("SPLICECRAFT_DEBUG", "").strip().lower() in (
    "1", "true", "yes"
):
    _log.setLevel(logging.DEBUG)
else:
    _log.setLevel(logging.INFO)
_log.propagate = False
# Scrub lone surrogates on EVERY record before any handler / test-capture sees
# them — keeps pytest-xdist's report serialization from crashing CI (see class).
_log.addFilter(_SurrogateScrubFilter())


def _repr_for_log(value, max_len: int = 100) -> str:
    """Compact, bounded `repr()` for log lines. Long sequences,
    nested dicts, and accidentally-huge values are truncated so a
    single chatty log line can't blow out the rotation window. The
    prefix preserves the type ("[123 items]") so a grepper can still
    find the relevant entry. Plasmid sequences are NOT logged (callers
    avoid passing them) but this guard keeps an accidental long-string
    log call from leaking content."""
    try:
        if isinstance(value, str):
            if len(value) > max_len:
                return f"{value[:max_len // 2]!r}…[{len(value)} chars]"
            return repr(value)
        if isinstance(value, (list, tuple, set)):
            n = len(value)
            if n > 10:
                return f"<{type(value).__name__} of {n} items>"
            return repr(value)[:max_len]
        if isinstance(value, dict):
            if len(value) > 10:
                return f"<dict with {len(value)} keys>"
            return repr(value)[:max_len]
        out = repr(value)
        if len(out) > max_len:
            return out[:max_len] + f"…[{len(out)} chars]"
        return out
    except Exception:
        return "<unrepr-able>"


def _log_event(event: str, *, _stacklevel: int = 2, **fields) -> None:
    """One-line AI-parseable structured event for diagnostic logs.

    Output line shape (JSON payload after the prefix):

        2026-05-15 12:34:56,789 [a3f2c1d8] INFO  splicecraft.action_save:53321 event save.ok {"rec":"pUC19","path":"/tmp/x.gb"}

    The event payload is JSON so any downstream parser (jq, Python
    json.loads, an LLM scanning the log on the user's behalf) can
    extract every field unambiguously without regex tricks against
    embedded whitespace or quotes. `stacklevel=2` makes the
    `funcName:lineno` prefix point at the originating caller, not at
    `_log_event` itself — so even an empty-field event names the
    method that emitted it.

    Sacred invariant: every string value over 200 chars is truncated
    before encoding, so a caller that accidentally passes raw
    sequence content can't leak bases into the log file. Non-scalar
    values pass through `_repr_for_log` which truncates and adds a
    size hint.

    Use this at click handlers, key actions, save/load/annotate
    boundaries — anywhere a user-visible state change happens. The
    output goes to the rotating log file at INFO level so a user
    pasting their log into a bug report tells the reader (or an
    AI assistant) exactly what was being done when the symptom hit.

    Event-name conventions (see also CLAUDE.md invariant #43):
      * `app.<area>.<verb>`  — user actions
            e.g. `app.save.trigger`, `app.library.add`,
            `app.feature.add`, `app.click_debug.toggle`.
      * `op.<area>.<verb>`   — heavy ops emitted by `@_timed`
            e.g. `op.fetch_genbank`, `op.gibson_simulate`,
            `op.pairwise_align`, `op.blast_search`, `op.hmmscan`,
            `op.annotation_transfer`.
      * `<noun>.<verb>`      — state changes
            e.g. `save.ok` / `save.failed`, `record.loaded`,
            `undo.trigger` / `undo.refused` / `undo.empty`,
            `redo.trigger` / `redo.refused` / `redo.empty`,
            `lock.acquired` / `lock.contended` / `lock.stale` /
            `lock.released`, `shutdown.drain.ok` /
            `shutdown.drain.timeout`, `agent.write.ok` /
            `agent.write.failed`, `agent.error`.

    Performance: short-circuits to a no-op when the logger isn't
    INFO-enabled, so every callsite pays one `isEnabledFor` check
    (~100 ns) in the happy path. The `logging` machinery is
    thread-safe, so this is safe to call from `@work` workers.
    """
    if not _log.isEnabledFor(logging.INFO):
        return
    if not fields:
        _log.info("event %s", event, stacklevel=_stacklevel)
        return
    safe: dict = {}
    for k, v in fields.items():
        if isinstance(v, str):
            safe[k] = v if len(v) <= 200 else (
                v[:100] + f"…[+{len(v) - 100}]"
            )
        elif isinstance(v, (int, float, bool)) or v is None:
            safe[k] = v
        elif isinstance(v, (list, tuple, dict)):
            # Lists / tuples / dicts of primitives JSON-encode as
            # arrays / objects (better for downstream parsers than
            # a stringified repr). Bail to bounded repr only when
            # the structure is too big to fit on one log line OR
            # carries un-JSON-encodable values.
            try:
                encoded = json.dumps(v, default=str,
                                      ensure_ascii=False)
                if len(encoded) <= 300:
                    safe[k] = list(v) if isinstance(v, tuple) else v
                else:
                    safe[k] = _repr_for_log(v, max_len=200)
            except (TypeError, ValueError):
                safe[k] = _repr_for_log(v, max_len=200)
        elif isinstance(v, (bytes, bytearray)):
            # Raw byte blobs (e.g. a .dna file body) are NEVER
            # rendered — only a size tag.
            safe[k] = f"<{type(v).__name__} len={len(v)}>"
        else:
            # Catch BioPython sequence-bearing classes by name
            # (avoids importing BioPython here, which would
            # already be in sys.modules anyway). `repr(Seq)` /
            # `repr(SeqRecord)` embed the first ~55 bases — the
            # sacred invariant says we must never log sequence
            # content, so render an opaque tag instead.
            cls_name = type(v).__name__
            if cls_name in ("SeqRecord", "Seq", "MutableSeq"):
                rid = getattr(v, "id", None)
                try:
                    sz = len(v)
                except (TypeError, ValueError):
                    sz = None
                parts = [f"<{cls_name}"]
                if rid:
                    parts.append(f"id={rid}")
                if sz is not None:
                    parts.append(f"len={sz}")
                safe[k] = " ".join(parts) + ">"
            else:
                # Path, datetime, arbitrary objects — bounded
                # repr. Final defence against an unexpected type
                # that happens to embed sequence in its repr.
                safe[k] = _repr_for_log(v, max_len=200)
    try:
        payload = json.dumps(safe, default=str,
                              separators=(",", ":"),
                              ensure_ascii=False)
    except (TypeError, ValueError):
        # Final-resort fallback for non-JSON-encodable values that
        # also slip past `_repr_for_log`. Should be unreachable.
        payload = repr(safe)[:500]
    _log.info("event %s %s", event, payload, stacklevel=_stacklevel)


# ── Decorators: structured action logging + perf timing (moved from hub, Phase D)
# `_log_event`/`_log` resolve in THIS module now, so event-capture tests patch
# `splicecraft_logging._log_event` (the patchable source moved WITH the decorators;
# see tests/test_logging.py::test_action_log_routes_through_patchable_log_event).
def _action_log(event_name: str):
    """Decorator that emits an INFO event at action-method entry
    plus the active record id (when present) so the log line tells
    the reader what the user was looking at, not just what they
    pressed.

    Applied to ``action_*`` methods on `PlasmidApp` so the per-action
    boilerplate stays a one-line decorator instead of an explicit
    `_log_event` call at the head of every body. Any exception in
    the logging path is swallowed — logging must never break the
    underlying action.

    Use as::

        @_action_log("app.fetch")
        def action_fetch(self):
            ...
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            # Level-gate FIRST so an INFO-suppressed run pays one
            # attribute lookup per action invocation (no getattr chain,
            # no dict build, no JSON encode). The same guard inside
            # `_log_event` covers the explicit-call path; this one
            # short-circuits the decorator wrapper too.
            if _log.isEnabledFor(logging.INFO):
                try:
                    ctx: dict = {}
                    rec = getattr(self, "_current_record", None)
                    if rec is not None:
                        rid = getattr(rec, "id", None) or getattr(
                            rec, "name", None
                        )
                        if rid:
                            ctx["rec"] = rid
                    # _stacklevel=3 so the logger's funcName:lineno
                    # points at the wrapped action method, not at this
                    # wrapper. Without it every action shows up in the
                    # log as "_action_log:531" — useless for tracing.
                    _log_event(event_name, _stacklevel=3, **ctx)
                except Exception:  # noqa: BLE001 — logging must never raise
                    # INV-73 (2026-05-25): the silent-swallow used to
                    # mean a malformed _log_event payload caused the
                    # entire log entry to vanish with no trace. Now
                    # we leave a debug-level breadcrumb on the std
                    # logger so a power user grepping the log can
                    # spot the regression — the action still runs
                    # so user-facing behaviour is untouched.
                    try:
                        _log.debug(
                            "_action_log decorator: event emit "
                            "failed for %s; action still ran",
                            event_name,
                        )
                    except Exception:  # noqa: BLE001 — last-resort
                        pass
            return func(self, *args, **kwargs)
        return wrapper
    return decorator


def _timed(path: str, threshold_ms: float = 0.0):
    """Decorator counterpart to `_log_timing` — wraps a whole function
    body in the same start/elapsed harness. Use on top-level heavy
    operations (NCBI fetch, primer3, BLAST, Gibson simulate) where
    every call is potentially diagnostic-worthy.

    Default `threshold_ms=0` means every call emits an event — useful
    for known-heavy paths where you always want a timestamp + duration
    in the log. Bump to a non-zero value to silence fast cases on
    paths that mostly return quickly.

    Naming: pass a `path` like ``"op.fetch_genbank"`` or
    ``"op.gibson_simulate"`` so the AI-parser can group all slow
    events by operation type.

    Performance: same level-gate short-circuit as `_log_event` — when
    INFO is suppressed the wrapper still measures (`perf_counter` is
    sub-microsecond) but `_log_event` early-returns before any JSON
    encode. Net cost per call when INFO is off: ~300 ns.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            t0 = _time.perf_counter()
            try:
                return func(*args, **kwargs)
            finally:
                dt_ms = (_time.perf_counter() - t0) * 1000
                if dt_ms >= threshold_ms:
                    # _stacklevel=3 so the funcName:lineno prefix
                    # points at the wrapped function, not at this
                    # wrapper closure.
                    _log_event("op.timed", _stacklevel=3,
                                path=path,
                                elapsed_ms=round(dt_ms, 1))
        return wrapper
    return decorator
