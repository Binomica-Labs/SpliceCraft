"""splicecraft_babs — the BABS chat engine (L1).

Pure, app-free core for the in-app **BABS** chat tab: a faithful reimplementation
of the Babs terminal-chat UX (persona, streaming answers, ``<think>``-stripping,
markdown-lite rendering, the ❤ context lifebar, slash commands) talking
**directly to a local Ollama server over stdlib HTTP** — no ``ollama`` /
``chromadb`` dependency — plus a model browser that lists installed Ollama
models, searches the HuggingFace Hub (GGUF) and pulls a chosen model on demand.

Why a sibling, not hub code: every network / protocol / markdown helper here is
app-free and unit-testable in isolation. The Textual UI (``BabsScreen`` + the
model-picker modal) stays hub-side because it is app-coupled (``push_screen``,
workers, ``call_from_thread``). This module imports only L0 siblings
(``net`` / ``util`` / ``logging`` / ``state``) and never reaches the hub.

SECURITY / hardening
--------------------
* Ollama is a **local** service (default ``127.0.0.1:11434``), so its calls
  deliberately bypass the public-host SSRF assertion in ``splicecraft_net``
  (which *refuses* loopback). They are still bounded on every axis: explicit
  connect/read timeouts, a hard response-size cap, a per-line cap and a total
  cap on streamed NDJSON, and a cooperative cancel ``Event``.
* The HuggingFace search is a **public** fetch and therefore goes through the
  shared hardened opener (``_build_hardened_url_opener``) — same SSRF / redirect
  / size discipline as every other outbound SpliceCraft fetch.
* The whole feature is gated OFF in web-demo mode at the screen level
  (``_demo_web_refuse``), so none of these sockets open on the public demo.
* No data-dir writes happen here — the engine is read/compute/network only.
"""
from __future__ import annotations

import functools
import json
import math
import os
import re
import urllib.error
import urllib.parse
import urllib.request

from rich.markup import escape as _rich_escape

import splicecraft_logging as _logging
import splicecraft_net as _net
import splicecraft_state as _state

_log = _logging._log


# ── Defaults / tunables ────────────────────────────────────────────────────────
# Mirror Babs' bb_config defaults so a user who knows Babs feels at home: same
# model floor, same deterministic temperature, same generous num_ctx.
DEFAULT_OLLAMA_HOST = "http://127.0.0.1:11434"
DEFAULT_CHAT_MODEL = "qwen2.5:7b"   # == bb_config.CHAT_MODEL default
DEFAULT_NUM_CTX = 4096              # == bb_config.ANSWER_NUM_CTX default
# Conservative ceiling when AUTO-sizing num_ctx to a model's real window: lift
# the working window above the 4096 floor (so the ❤ lifebar and chat memory
# reflect a large-context model) without requesting a giant KV cache that could
# balloon RAM on modest hardware.
MAX_AUTO_NUM_CTX = 16384
DEFAULT_TEMP = 0.0                  # == bb_config.ANSWER_TEMP default (deterministic)

# Timeouts (seconds). The *connect/short-call* budget is tight (a missing Ollama
# should fail fast, not hang the worker); chat + pull are long because CPU
# inference and multi-GB downloads legitimately take minutes.
_SHORT_TIMEOUT = 8       # /api/tags, /api/show, connectivity ping
_CHAT_TIMEOUT = 600      # == bb_config.CHAT_TIMEOUT; a wedged model can't hang forever
_PULL_TIMEOUT = 900      # a big model layer over a slow link
_HF_TIMEOUT = 20

# Response-size caps (bytes). Defend against a misbehaving / hijacked endpoint
# streaming gigabytes at us (mirrors the _*_MAX_RESPONSE_BYTES discipline in
# splicecraft_net). Generous vs. real payloads: /api/tags for 100 models is tens
# of KB; one HF search page is a few hundred KB.
_TAGS_MAX_BYTES = 8 * 1024 * 1024
_SHOW_MAX_BYTES = 8 * 1024 * 1024
_HF_MAX_BYTES = 8 * 1024 * 1024
_STREAM_MAX_TOTAL_BYTES = 256 * 1024 * 1024   # a whole chat/pull stream
_STREAM_MAX_LINE_BYTES = 4 * 1024 * 1024      # one NDJSON line

# Input guard: clamp a pasted prompt so a megabyte blob can't blow num_ctx (and
# the UI TextArea) — mirrors the modal paste caps. The model would silently
# evict context past num_ctx anyway.
MAX_PROMPT_CHARS = 32000

# A reasonable ceiling on how many models a search returns to the UI table.
HF_SEARCH_LIMIT = 30


# ── The Babs persona (system prompt) ───────────────────────────────────────────
# Mirrors bb_config / rag_bot.SYSTEM's persona verbatim where it still applies,
# with the RAG-only clauses ("answer ONLY from the retrieved context", "cite by
# bracket number") removed — there is no corpus here, so demanding citations
# would make every answer apologise. The injection-safety spirit is kept. The
# user can edit this live with `/system`; it is persisted by the screen.
BABS_SYSTEM = (
    "You are Babs, an expert assistant on plant tissue culture and plant "
    "bioengineering: in vitro regeneration and callus culture; transformation by "
    "Agrobacterium, biolistics and silicon-carbide whisker; molecular cloning and "
    "DNA assembly (Golden Gate, MoClo, GoldenBraid, Gibson, traditional cloning); "
    "genome integration (transposases such as Himar1, recombinases, homologous "
    "recombination) and CRISPR/Cas editing; and the mechanisms and optimization of "
    "plant growth regulators. You are embedded inside SpliceCraft, a terminal "
    "plasmid-design workbench, as its chat assistant. Give clear, practical, "
    "well-organized answers. This is a direct conversation (no document retrieval), "
    "so answer from your own knowledge and be honest about its limits: if you are "
    "unsure, or a number/protocol needs a primary source, say so plainly rather "
    "than inventing specifics."
)


# ── Errors ─────────────────────────────────────────────────────────────────────
class BabsError(Exception):
    """Any BABS-engine failure with a user-facing, actionable message."""


class OllamaUnavailable(BabsError):
    """Ollama isn't reachable (not installed / not running / wrong host)."""


_OLLAMA_START_HINT = (
    "Ollama isn't reachable. Install it from https://ollama.com, then start it "
    "with `ollama serve` (it usually runs on 127.0.0.1:11434). Set a different "
    "endpoint with $SPLICECRAFT_OLLAMA_HOST or $OLLAMA_HOST."
)


def _conn_error(exc: Exception) -> OllamaUnavailable:
    """Map a urllib failure talking to Ollama into an actionable message."""
    reason = getattr(exc, "reason", None) or exc
    return OllamaUnavailable(f"{_OLLAMA_START_HINT}\n  (reason: {reason})")


def _http_error_message(exc) -> str:
    """Pull a clean message out of an Ollama HTTP error RESPONSE (the server WAS
    reached but returned a non-2xx — e.g. a 404 'model not found' when chatting
    with an un-pulled model). Reads the body BOUNDED so a giant error page can't
    flood the UI, and prefers the JSON ``error`` field Ollama returns."""
    body = ""
    try:
        body = exc.read(8192).decode("utf-8", "replace")
    except Exception:
        pass
    detail = ""
    try:
        parsed = json.loads(body)
        if isinstance(parsed, dict):
            detail = str(parsed.get("error") or "").strip()
    except Exception:
        pass
    if not detail:
        detail = body.strip()[:200] or "no detail"
    return f"Ollama error (HTTP {getattr(exc, 'code', '?')}): {detail}"


# ── Host resolution ────────────────────────────────────────────────────────────
def ollama_base() -> str:
    """Resolve the Ollama base URL from the environment, defaulting to localhost.

    Accepts ``host:port``, ``http://host:port`` or ``https://host:port`` in
    ``$SPLICECRAFT_OLLAMA_HOST`` / ``$OLLAMA_HOST`` (the latter is Ollama's own
    convention). A malformed value degrades to the safe loopback default rather
    than raising at call time."""
    raw = (os.getenv("SPLICECRAFT_OLLAMA_HOST")
           or os.getenv("OLLAMA_HOST") or "").strip()
    if not raw:
        return DEFAULT_OLLAMA_HOST
    if "://" not in raw:
        raw = "http://" + raw
    try:
        p = urllib.parse.urlparse(raw)
        # `.port` raises ValueError on an out-of-range / non-numeric port
        # (e.g. "host:999999" or "host:abc"); read it INSIDE the guard so a
        # malformed $OLLAMA_HOST / $SPLICECRAFT_OLLAMA_HOST degrades to the
        # loopback default instead of raising from every Babs call
        # (ping / list / chat / pull / delete), per this function's contract.
        if p.scheme not in ("http", "https") or not p.hostname:
            return DEFAULT_OLLAMA_HOST
        port = p.port or (443 if p.scheme == "https" else 11434)
        host = p.hostname
    except ValueError:
        return DEFAULT_OLLAMA_HOST
    if ":" in host:                 # an IPv6 literal needs brackets in a URL
        host = f"[{host}]"
    return f"{p.scheme}://{host}:{port}"


def _user_agent() -> str:
    return f"SpliceCraft/{_state._sc_version or '?'} (BABS chat client)"


class _NoRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Refuse to follow 3xx redirects on the Ollama opener.

    The Ollama REST API never legitimately redirects. Following one blindly (what
    the default opener does) would let a hostile or MITM'd *remote* endpoint —
    only reachable when ``$OLLAMA_HOST`` points off-loopback — bounce the client
    into the internal network (a classic SSRF pivot to ``169.254.169.254`` /
    intranet hosts). The default-loopback config is unaffected. Surfaced as an
    ``HTTPError`` so the existing error mapping shows a clean message."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: D401
        raise urllib.error.HTTPError(
            req.full_url, code,
            f"Ollama endpoint attempted a redirect to {newurl!r} — refused "
            f"(possible SSRF). Point $OLLAMA_HOST straight at the server.",
            headers, fp)


@functools.lru_cache(maxsize=1)
def _ollama_opener() -> "urllib.request.OpenerDirector":
    """A urllib opener for the local Ollama server identical to the default
    EXCEPT it refuses redirects (see :class:`_NoRedirectHandler`). Cached — one
    opener for the process. Localhost still works; only off-loopback redirect
    chains are blocked."""
    return urllib.request.build_opener(_NoRedirectHandler())


# ── Low-level HTTP (localhost-allowed; bounded) ────────────────────────────────
def _request_json(path: str, *, method: str = "GET", payload: "dict | None" = None,
                  timeout: float = _SHORT_TIMEOUT, max_bytes: int = _TAGS_MAX_BYTES) -> dict:
    """One bounded JSON call to the local Ollama server. Uses the *default*
    urllib opener on purpose (NOT the public-host-asserting hardened opener,
    which refuses loopback). Size-capped via ``read(max_bytes + 1)`` per
    [PIT-20]. Raises ``OllamaUnavailable`` on a transport error."""
    url = ollama_base() + path
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    headers = {"User-Agent": _user_agent(), "Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with _ollama_opener().open(req, timeout=timeout) as resp:
            raw = resp.read(max_bytes + 1)
    except urllib.error.HTTPError as exc:
        raise BabsError(_http_error_message(exc))
    except (urllib.error.URLError, OSError) as exc:
        raise _conn_error(exc)
    if len(raw) > max_bytes:
        raise BabsError(f"Ollama response exceeded {max_bytes // (1024*1024)} MB — aborted.")
    try:
        return json.loads(raw.decode("utf-8", "replace"))
    except ValueError as exc:
        raise BabsError(f"Ollama returned a non-JSON response: {exc}")


def _iter_ndjson(resp, *, cancel, max_total: int, max_line: int):
    """Yield parsed objects from a newline-delimited-JSON HTTP stream, bounded on
    BOTH total bytes and single-line bytes so a runaway / malicious stream can't
    OOM the box. Checks the cancel ``Event`` between chunks (a blocked socket
    read is interrupted separately by closing the response). Unparseable lines
    are skipped, never fatal."""
    total = 0
    buf = b""
    while True:
        if cancel is not None and cancel.is_set():
            return
        try:
            chunk = resp.read(16384)
        except (urllib.error.URLError, OSError):
            # Socket closed (often our own cancel) — stop cleanly.
            return
        if not chunk:
            break
        total += len(chunk)
        if total > max_total:
            raise BabsError("Ollama stream exceeded its size cap — aborted.")
        buf += chunk
        while b"\n" in buf:
            line, buf = buf.split(b"\n", 1)
            if len(line) > max_line:
                raise BabsError("Ollama stream emitted an oversized line — aborted.")
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line.decode("utf-8", "replace"))
            except ValueError:
                continue
            if isinstance(obj, dict):     # Ollama NDJSON is always objects; ignore stray scalars/arrays
                yield obj
        if len(buf) > max_line:
            raise BabsError("Ollama stream emitted an oversized line — aborted.")
    tail = buf.strip()
    if tail:
        try:
            obj = json.loads(tail.decode("utf-8", "replace"))
        except ValueError:
            obj = None
        if isinstance(obj, dict):
            yield obj


def _open_stream(path: str, payload: dict, *, timeout: float, register=None):
    """POST a streaming request to Ollama and return the open response. ``register``
    (if given) receives the raw response object so the caller can ``.close()`` it
    from another thread to abort a blocked read (responsive cancel)."""
    url = ollama_base() + path
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=data, method="POST",
        headers={"User-Agent": _user_agent(), "Content-Type": "application/json",
                 "Accept": "application/x-ndjson"},
    )
    try:
        resp = _ollama_opener().open(req, timeout=timeout)
    except urllib.error.HTTPError as exc:
        raise BabsError(_http_error_message(exc))
    except (urllib.error.URLError, OSError) as exc:
        raise _conn_error(exc)
    if register is not None:
        try:
            register(resp)
        except Exception:   # pragma: no cover - registration must never break the stream
            pass
    return resp


# ── Ollama API ─────────────────────────────────────────────────────────────────
def list_installed(*, timeout: float = _SHORT_TIMEOUT) -> "list[dict]":
    """Models already pulled into Ollama (``GET /api/tags``), normalized to
    ``{name, size, family, params, quant, modified}``. Raises
    ``OllamaUnavailable`` if Ollama isn't up."""
    data = _request_json("/api/tags", timeout=timeout, max_bytes=_TAGS_MAX_BYTES)
    out = []
    for m in (data.get("models") or []):
        if not isinstance(m, dict):
            continue
        det = m.get("details") or {}
        name = m.get("name") or m.get("model") or ""
        if not name:
            continue
        out.append({
            "name": name,
            "size": m.get("size") or 0,
            "family": (det.get("family") or "").strip(),
            "params": (det.get("parameter_size") or "").strip(),
            "quant": (det.get("quantization_level") or "").strip(),
            "modified": (m.get("modified_at") or "").strip(),
        })
    out.sort(key=lambda d: d["name"].lower())
    return out


def show(name: str, *, timeout: float = _SHORT_TIMEOUT) -> dict:
    """Best-effort model metadata (``POST /api/show``). Returns ``{}`` on any
    failure — metadata is a nice-to-have, never load-bearing."""
    name = (name or "").strip()
    if not name:
        return {}
    try:
        return _request_json("/api/show", method="POST",
                             payload={"model": name, "name": name},
                             timeout=timeout, max_bytes=_SHOW_MAX_BYTES)
    except BabsError:
        return {}


def context_length_from_show(meta: "dict | None") -> "int | None":
    """Pull the model's real context window from a ``/api/show`` payload, or
    ``None``. Ollama nests it as ``model_info["<arch>.context_length"]`` (e.g.
    ``qwen2.context_length``); scan for the architecture-prefixed key so no arch
    is hard-coded. Lets the UI size the ❤ lifebar to the ACTUAL window instead of
    the 4096 default (a 128k model otherwise reads as nearly-full far too soon)."""
    if not isinstance(meta, dict):
        return None
    info = meta.get("model_info")
    if isinstance(info, dict):
        for k, v in info.items():
            if isinstance(k, str) and k.endswith(".context_length"):
                n = _to_int(v)
                if n and n > 0:
                    return n
    return None


def ping(*, timeout: float = _SHORT_TIMEOUT) -> bool:
    """True iff the Ollama server answers. Never raises."""
    try:
        list_installed(timeout=timeout)
        return True
    except BabsError:
        return False


def model_is_installed(name: str, installed: "list[dict] | None" = None) -> bool:
    """Match a model name against the installed set with Ollama's exact
    semantics: a tagless name (``qwen2.5``) means ``qwen2.5:latest``, so it
    counts as installed only if THAT specific tag is present — NOT merely because
    some other tag (``qwen2.5:7b``) of the same family exists (pulling it would
    fetch a different artifact). An explicit tag must match exactly."""
    name = (name or "").strip()
    if not name:
        return False
    try:
        names = {m["name"] for m in (installed if installed is not None else list_installed())}
    except BabsError:
        return False
    if name in names:
        return True
    if ":" not in name:
        return f"{name}:latest" in names
    return False


def chat_stream(model: str, messages: "list[dict]", *, tools: "list[dict] | None" = None,
                options: "dict | None" = None,
                think: "bool | None" = None, cancel=None, register=None,
                timeout: float = _CHAT_TIMEOUT):
    """Stream a chat completion (``POST /api/chat``). Yields
    ``{"content": str, "thinking": str, "tool_calls": list, "done": bool,
    "done_reason": str|None, "error": str|None}`` per chunk — ``content`` is
    visible answer text, ``thinking`` is the model's separate reasoning channel
    (newer Ollama) which the UI hides by default, ``tool_calls`` is any
    function-call(s) the model emitted this chunk (present only when ``tools``
    is supplied AND the model decides to call one).

    ``tools`` (optional) is the Ollama function-tool manifest — a list of
    ``{"type":"function","function":{name, description, parameters}}`` — that
    powers the autonomous Babs tool-loop. Streamed tool calls require a
    tool-capable model (e.g. qwen2.5) and Ollama ≳0.4; older servers simply
    never populate ``tool_calls`` and the loop degrades to plain chat.

    ``cancel`` (Event) + ``register`` (receives the response for ``.close()``)
    give a responsive stop."""
    payload: dict = {"model": model, "messages": messages, "stream": True,
                     "options": options or {"temperature": DEFAULT_TEMP, "num_ctx": DEFAULT_NUM_CTX}}
    if tools:
        payload["tools"] = tools
    if think is not None:
        payload["think"] = bool(think)
    resp = _open_stream("/api/chat", payload, timeout=timeout, register=register)
    try:
        for obj in _iter_ndjson(resp, cancel=cancel,
                                max_total=_STREAM_MAX_TOTAL_BYTES,
                                max_line=_STREAM_MAX_LINE_BYTES):
            msg = obj.get("message") or {}
            tcs = msg.get("tool_calls")
            yield {
                "content": msg.get("content") or "",
                "thinking": msg.get("thinking") or "",
                "tool_calls": tcs if isinstance(tcs, list) else [],
                "done": bool(obj.get("done")),
                "done_reason": obj.get("done_reason"),
                "error": obj.get("error"),
            }
    finally:
        try:
            resp.close()
        except Exception:
            pass


def pull_stream(name: str, *, cancel=None, register=None, timeout: float = _PULL_TIMEOUT):
    """Stream a model download (``POST /api/pull``). Yields the raw Ollama
    progress objects — ``{"status", "digest", "total", "completed", "error"}``.
    Works for both registry names (``qwen2.5:7b``) and HuggingFace GGUF repos
    (``hf.co/<owner>/<repo>[:quant]``)."""
    payload = {"model": name, "name": name, "stream": True}
    resp = _open_stream("/api/pull", payload, timeout=timeout, register=register)
    try:
        for obj in _iter_ndjson(resp, cancel=cancel,
                                max_total=_STREAM_MAX_TOTAL_BYTES,
                                max_line=_STREAM_MAX_LINE_BYTES):
            yield obj
    finally:
        try:
            resp.close()
        except Exception:
            pass


def delete_model(name: str, *, timeout: float = _SHORT_TIMEOUT) -> bool:
    """Uninstall a model from the local Ollama store (``DELETE /api/delete``).

    Ollama answers ``200`` with an EMPTY body on success, so this does NOT route
    through ``_request_json`` (which expects JSON and would raise on the empty
    body). Returns True on success. Raises ``OllamaUnavailable`` if the server
    is unreachable, ``BabsError`` on an HTTP error (404 = model not installed)."""
    name = (name or "").strip()
    if not name:
        raise BabsError("no model name given")
    url = ollama_base() + "/api/delete"
    # Send both keys: newer Ollama wants "model", older wants "name".
    data = json.dumps({"model": name, "name": name}).encode("utf-8")
    req = urllib.request.Request(
        url, data=data, method="DELETE",
        headers={"User-Agent": _user_agent(), "Content-Type": "application/json"})
    try:
        with _ollama_opener().open(req, timeout=timeout) as resp:
            resp.read(1024)        # drain (bounded); body is empty on success
    except urllib.error.HTTPError as exc:
        raise BabsError(_http_error_message(exc))
    except (urllib.error.URLError, OSError) as exc:
        raise _conn_error(exc)
    return True


# ── HuggingFace Hub search (public; hardened opener) ───────────────────────────
def hf_search_gguf(query: str, *, limit: int = HF_SEARCH_LIMIT, opener=None,
                   timeout: float = _HF_TIMEOUT) -> "list[dict]":
    """Search the HuggingFace Hub for GGUF model repos (the ones Ollama can run
    via ``ollama pull hf.co/<id>``). Official JSON API, fetched through the shared
    SSRF-hardened opener. ``opener`` is injectable for tests. Returns
    ``[{id, downloads, likes, gated, pull}]`` sorted by downloads; ``[]`` on a
    blank query."""
    query = (query or "").strip()
    if not query:
        return []
    limit = max(1, min(int(limit), 100))
    qs = urllib.parse.urlencode({
        "search": query, "filter": "gguf", "sort": "downloads",
        "direction": "-1", "limit": str(limit), "full": "false",
    })
    url = f"https://huggingface.co/api/models?{qs}"
    opener = opener or _net._build_hardened_url_opener()
    req = urllib.request.Request(
        url, headers={"User-Agent": _user_agent(), "Accept": "application/json"})
    try:
        with opener.open(req, timeout=timeout) as resp:
            raw = resp.read(_HF_MAX_BYTES + 1)
    except (urllib.error.URLError, OSError) as exc:
        raise BabsError(f"HuggingFace search failed: {getattr(exc, 'reason', exc)}")
    if len(raw) > _HF_MAX_BYTES:
        raise BabsError("HuggingFace response exceeded its size cap — aborted.")
    try:
        data = json.loads(raw.decode("utf-8", "replace"))
    except ValueError as exc:
        raise BabsError(f"HuggingFace returned a non-JSON response: {exc}")
    if not isinstance(data, list):
        return []
    out = []
    for m in data:
        if not isinstance(m, dict):
            continue
        mid = m.get("id") or m.get("modelId") or ""
        if not isinstance(mid, str):     # hostile/malformed JSON: id not a string
            continue
        mid = mid.strip()
        if not mid:
            continue
        out.append({
            "id": mid,
            "downloads": _to_int(m.get("downloads")) or 0,
            "likes": _to_int(m.get("likes")) or 0,
            "gated": bool(m.get("gated")),
            "pull": f"hf.co/{mid}",
        })
    out.sort(key=lambda d: d["downloads"], reverse=True)
    return out


# ── Reasoning-channel handling (mirror bb_config.strip_think + rag_bot stream) ──
_THINK_BLOCK = re.compile(r"<think>.*?</think>\s*", re.DOTALL | re.IGNORECASE)
_THINK_OPEN, _THINK_CLOSE = "<think>", "</think>"
# Case-insensitive locators for the streaming stripper (a model may emit
# <THINK>/<Think>). Regex search keeps correct indices on the ORIGINAL buffer
# regardless of any Unicode elsewhere in it (str.lower() can change length).
_THINK_OPEN_RE = re.compile(re.escape(_THINK_OPEN), re.IGNORECASE)
_THINK_CLOSE_RE = re.compile(re.escape(_THINK_CLOSE), re.IGNORECASE)


def strip_think(text: str) -> str:
    """Drop ``<think>...</think>`` reasoning a model may inline, tolerating an
    unclosed tag (generation cut off mid-thought). Byte-for-byte port of
    ``bb_config.strip_think`` minus the JSON code-fence unwrap (irrelevant to a
    prose chat). No-op for a normal model."""
    if not text:
        return text
    if "<think>" in text.lower():
        text = _THINK_BLOCK.sub("", text)
        low = text.lower()
        if "<think>" in low and "</think>" not in low:
            text = text[:low.index("<think>")]
        text = text.strip()
    return text


def _held_suffix(buf: str, tag: str) -> int:
    """Length of the longest suffix of ``buf`` that is a proper prefix of ``tag``
    (compared case-insensitively; ``tag`` is lowercase) — bytes held back because
    they might begin a tag split across two chunks."""
    for k in range(min(len(buf), len(tag) - 1), 0, -1):
        if buf[-k:].lower() == tag[:k]:
            return k
    return 0


class ThinkStripper:
    """Streaming twin of :func:`strip_think`: feed answer deltas, get back only
    the *visible* text with any inline ``<think>...</think>`` span swallowed
    token-by-token — even when a tag marker is split across chunks. A normal
    model (no such block) passes through byte-for-byte. Port of rag_bot's
    ``_visible_stream`` as a stateful feeder (Ollama's separate ``thinking``
    field is handled upstream; this is the belt-and-suspenders for models that
    still inline the tag in ``content``)."""

    def __init__(self) -> None:
        self.buf = ""
        self.in_think = False

    def feed(self, text: str) -> str:
        if not text:
            return ""
        self.buf += text
        out: "list[str]" = []
        while self.buf:
            if self.in_think:
                m = _THINK_CLOSE_RE.search(self.buf)
                if m is None:
                    hold = _held_suffix(self.buf, _THINK_CLOSE)
                    self.buf = self.buf[len(self.buf) - hold:] if hold else ""
                    break
                self.buf = self.buf[m.end():].lstrip()
                self.in_think = False
            else:
                m = _THINK_OPEN_RE.search(self.buf)
                if m is None:
                    hold = _held_suffix(self.buf, _THINK_OPEN)
                    visible = self.buf[:len(self.buf) - hold] if hold else self.buf
                    if visible:
                        out.append(visible)
                    self.buf = self.buf[len(self.buf) - hold:] if hold else ""
                    break
                if m.start():
                    out.append(self.buf[:m.start()])
                self.buf = self.buf[m.end():]
                self.in_think = True
        return "".join(out)

    def flush(self) -> str:
        """Emit any held tail that turned out not to be a tag start."""
        tail = self.buf if not self.in_think else ""
        self.buf = ""
        return tail


# ── Context lifebar (mirror rag_bot) ───────────────────────────────────────────
_TOKENISH_RE = re.compile(r"\w+|[^\w\s]")
_CTX_SAFETY = 0.92      # never let a turn's prompt exceed this fraction of num_ctx
_ANSWER_RESERVE = 512   # tokens kept free for the reply when sizing memory budget


@functools.lru_cache(maxsize=2048)
def est_tokens(text: str) -> int:
    """Conservative BPE-token estimate without a model call (port of rag_bot's
    ``_est_tokens``): bias HIGH (undercounting overflows num_ctx). lru_cache keyed
    on the exact string makes re-renders free and self-staleguarding."""
    if not text:
        return 1
    by_char = len(text) / 3.6
    by_piece = len(_TOKENISH_RE.findall(text))
    return max(1, int(round(max(by_char, by_piece))))


def history_tokens(history: "list[dict] | None") -> int:
    return sum(est_tokens(m.get("content", "")) for m in (history or []))


def lifebar(history: "list[dict] | None", reserve_tokens: int,
            num_ctx: int = DEFAULT_NUM_CTX, width: int = 20) -> dict:
    """Compute the ❤ context gauge. 'Memory budget' is the window minus a reserve
    for the system prompt and room to answer, so an empty chat reads ~100% and it
    depletes as turns accumulate. Returns render-ready fields (the UI paints the
    bar; we don't emit ANSI here)."""
    mem_budget = max(256, int(num_ctx * _CTX_SAFETY) - reserve_tokens - _ANSWER_RESERVE)
    used = history_tokens(history)
    remaining = max(0.0, min(1.0, 1 - used / mem_budget))
    filled = int(round(remaining * width))
    level = "ok" if remaining > 0.5 else ("warn" if remaining > 0.2 else "crit")
    return {
        "remaining": remaining,
        "pct": int(round(remaining * 100)),
        "filled": filled,
        "width": width,
        "turns": len(history or []) // 2,
        "level": level,
    }


def trim_history(history: "list[dict]", reserve_tokens: int,
                 num_ctx: int = DEFAULT_NUM_CTX, query_tokens: int = 0) -> bool:
    """Drop the oldest user+assistant pairs in place until this turn fits the
    window (port of rag_bot's hard-trim backstop). Returns True if it trimmed.

    ``reserve_tokens`` is the system-prompt budget and ``query_tokens`` THIS
    turn's user message — both are counted (alongside ``_ANSWER_RESERVE`` room
    for the reply) so the assembled prompt can't overflow ``num_ctx`` once the
    query is appended. Pre-fix the query was omitted, so a large paste pushed the
    total past the window and llama.cpp silently evicted the FRONT — i.e. the
    ``BABS_SYSTEM`` persona. The ceiling now mirrors :func:`lifebar` exactly."""
    ceiling = max(256, int(num_ctx * _CTX_SAFETY) - _ANSWER_RESERVE)
    trimmed = False
    while (len(history) >= 2
           and reserve_tokens + query_tokens + history_tokens(history) > ceiling):
        del history[0:2]
        trimmed = True
    return trimmed


def build_messages(system: str, history: "list[dict] | None", query: str) -> "list[dict]":
    """Assemble the Ollama ``messages`` array: system + prior turns + this query."""
    msgs = [{"role": "system", "content": system or BABS_SYSTEM}]
    if history:
        msgs.extend(history)
    msgs.append({"role": "user", "content": query})
    return msgs


# ── Markdown-lite → Rich markup (for the transcript Static) ─────────────────────
_BOLD_RE = re.compile(r"\*\*(.+?)\*\*")
_CODE_RE = re.compile(r"`([^`]+)`")
_HEADING_RE = re.compile(r"^(#{1,6})\s+(.*)$")
_BULLET_RE = re.compile(r"^(\s*)[-*+]\s+(.*)$")


def _inline_md(s: str) -> str:
    """Escape Rich metacharacters in a line, then re-introduce **bold** and
    `code`. Order matters: escape FIRST (so user ``[1]`` shows literally), then
    add real markup tags (``**`` / `` ` `` are untouched by rich.escape)."""
    s = _rich_escape(s)
    s = _BOLD_RE.sub(r"[b]\1[/b]", s)
    s = _CODE_RE.sub(r"[cyan]\1[/cyan]", s)
    return s


def md_to_rich(text: str) -> str:
    """Render a markdown-lite answer as a Rich-markup string for a ``Static`` —
    **bold**, `code`, ``` fenced blocks ```, # headings and - bullets. Modest by
    design (no per-token re-parse cost, fully under our control, unit-testable);
    everything else passes through escaped. Safe to call on partial/streaming
    text (an unterminated fence just styles the tail as code)."""
    out: "list[str]" = []
    in_fence = False
    for line in (text or "").split("\n"):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            out.append("[dim cyan]" + _rich_escape(line) + "[/dim cyan]")
            continue
        mh = _HEADING_RE.match(line)
        if mh:
            out.append("[b]" + _inline_md(mh.group(2)) + "[/b]")
            continue
        mb = _BULLET_RE.match(line)
        if mb:
            out.append(mb.group(1) + "  • " + _inline_md(mb.group(2)))
            continue
        out.append(_inline_md(line))
    return "\n".join(out)


# ── Slash commands ─────────────────────────────────────────────────────────────
# The applicable subset of Babs' REPL commands (the corpus-only ones — /explore,
# /sources, /jobs, /get, /update, /kill, /n, /rerank, /stats — don't exist
# without a corpus). Kept in Babs' help layout so the UX reads the same.
HELP_COMMANDS = [
    ("/help  /?", "show this help"),
    ("/models", "browse + pick a model — installed · HuggingFace GGUF (pull on demand)"),
    ("/model [name]", "show the active model, or switch to <name> (pulls if needed)"),
    ("/system [text]", "show the system prompt (persona), or set it to <text>"),
    ("/recall [query]", "toggle corpus recall, or search your corpus and show the passages"),
    ("/ingest <url>", "read a page and KEEP it — adds it to Babs' knowledge corpus"),
    ("/learn <topic>", "grow the corpus about a topic (focused crawl; watch it in the Learn tab)"),
    ("/remember <fact>", "save a durable fact Babs recalls in every future session"),
    ("/memory", "show what Babs remembers (persistent notes)"),
    ("/forget <slug>", "remove one memory (/memory lists the slugs)"),
    ("/temp [0-2]", "show or set the answer temperature (0 = deterministic)"),
    ("/think [on|off]", "show or hide the model's reasoning (<think>) — default off"),
    ("/context  /ctx", "show the context lifebar (chat memory left before /reset)"),
    ("/retry", "regenerate the last answer"),
    ("/reset", "forget the conversation so far (clear chat memory)"),
    ("/clear", "clear the transcript on screen (keeps memory)"),
    ("/agent", "toggle agent mode — let Babs drive SpliceCraft (call its endpoints)"),
    ("/autonomy [ask|auto|readonly|off]", "set the write policy in agent mode"),
    ("/agentmodel [auto|chat|name]", "which model runs agent/tool turns (default: a fast tool-capable one)"),
    ("/exit  /quit  /q", "close the BABS tab"),
]


def parse_command(raw: str) -> "tuple[str, str]":
    """Split a ``/command rest`` line into ``(cmd_lower, rest)``. A bare exit word
    (``exit``/``quit``/``q``) is normalized to its slash form so the dispatcher
    has one path."""
    body = (raw or "").strip()
    if body.lower() in ("exit", "quit", "bye", "q"):
        return "exit", ""
    body = body.lstrip("/").strip()
    cmd, _, rest = body.partition(" ")
    return cmd.lower(), rest.strip()


def toggle_value(arg: str, current: bool) -> bool:
    """`on|off|1|0|...` sets explicitly; anything else flips ``current``."""
    a = (arg or "").strip().lower()
    if a in ("on", "1", "true", "yes", "y"):
        return True
    if a in ("off", "0", "false", "no", "n"):
        return False
    return not current


# ── Display helpers ─────────────────────────────────────────────────────────────
def fmt_size(nbytes: "int | float") -> str:
    """Human file size for the model table ('4.4 GB')."""
    try:
        n = float(nbytes)
    except (TypeError, ValueError, OverflowError):
        # OverflowError: a hijacked endpoint can deliver `size` as a
        # 300+-digit JSON integer, which `float()` rejects with OverflowError.
        return "?"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit in ("B", "KB") else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def clamp_prompt(text: str, *, num_ctx: "int | None" = None,
                 reserve_tokens: int = 0) -> "tuple[str, bool]":
    """Clamp an over-long prompt; returns ``(text, truncated)``.

    Always enforces the ``MAX_PROMPT_CHARS`` paste-bomb guard. When ``num_ctx``
    is given it ALSO clamps the query so it can't, on its own, exceed the room
    left in the window (``num_ctx*_CTX_SAFETY − reserve_tokens − _ANSWER_RESERVE``)
    — otherwise a single huge message forces llama.cpp to evict the persona even
    after history is trimmed. The token→char trim chops the tail until it fits
    (``est_tokens`` is lru-cached, so the few iterations are cheap)."""
    text = text or ""
    truncated = False
    if len(text) > MAX_PROMPT_CHARS:
        text = text[:MAX_PROMPT_CHARS]
        truncated = True
    if num_ctx:
        budget = max(128, int(num_ctx * _CTX_SAFETY) - reserve_tokens - _ANSWER_RESERVE)
        while text and est_tokens(text) > budget:
            text = text[:max(1, int(len(text) * 0.9))]
            truncated = True
    return text, truncated


def _to_float(x: object) -> "float | None":
    """Best-effort float coercion → ``None`` for None / non-numeric / a bad
    numeric string / a non-finite (inf, nan), so a malformed Ollama or HF value
    never raises and never propagates a bogus inf/nan into size / speed /
    fraction math. Type-clean stand-in for ``try: float(x) except ...`` (pyright
    models ``float()``'s argument strictly)."""
    if isinstance(x, (int, float, str)):
        try:
            f = float(x)
        except (ValueError, OverflowError):
            # OverflowError: a huge (300+-digit) JSON integer from a
            # misbehaving Ollama/HF endpoint can't be cast to float.
            return None
        return f if math.isfinite(f) else None
    return None


def _to_int(x: object) -> "int | None":
    """Best-effort int coercion via :func:`_to_float` (so it inherits the
    inf/nan rejection). ``None`` for non-numeric — used to tolerate HuggingFace
    JSON whose ``downloads``/``likes`` may arrive as strings or junk."""
    f = _to_float(x)
    return None if f is None else int(f)


def pull_progress_fraction(obj: dict) -> "float | None":
    """Extract a 0..1 download fraction from an Ollama pull progress object, or
    None when the phase has no measurable total (e.g. 'pulling manifest')."""
    total = _to_float(obj.get("total"))
    completed = _to_float(obj.get("completed"))
    if total is None or completed is None or total <= 0:
        return None
    return max(0.0, min(1.0, completed / total))


def fmt_speed(bps: "int | float | None") -> str:
    """Human download rate ('18.4 MB/s'); '' when not yet measurable so the UI
    can simply omit it on the first ('pulling manifest') ticks."""
    n = _to_float(bps)
    if n is None or not (n > 0):
        return ""
    return f"{fmt_size(n)}/s"


class PullMeter:
    """Turns the raw Ollama ``/api/pull`` progress stream into a human readout
    with a smoothed download speed, so a user can SEE a multi-GB pull is real
    and ongoing rather than staring at a bar that only ticks a percent.

    PURE + deterministic: the caller passes a monotonic timestamp on every
    ``update`` (the single time source, INV-78), so tests drive it with a
    synthetic clock and need no network or real clock. One meter per pull.

    Speed is averaged over a trailing ``window`` of wall-clock seconds within
    the CURRENT layer; the baseline resets when Ollama moves to a new layer
    (``digest`` changes) or the byte count jumps backwards, so a multi-layer
    pull never reports a bogus negative or spiked rate."""

    def __init__(self, *, window: float = 3.0) -> None:
        self._window = float(window)
        self._digest: "str | None" = None
        self._samples: "list[tuple[float, float]]" = []   # (t, completed) this layer

    def update(self, obj: dict, t: float) -> dict:
        """Feed one Ollama progress object + a monotonic timestamp. Returns
        ``{status, fraction, completed, total, bps, text}`` — ``bps`` is None
        until two byte samples in one layer make a rate measurable, ``text`` is
        the ready-to-display status line (status · bytes · percent · speed)."""
        status = (obj.get("status") or "").strip()
        frac = pull_progress_fraction(obj)
        completed = _to_float(obj.get("completed"))
        total = _to_float(obj.get("total"))
        digest = obj.get("digest")

        bps = None
        if completed is not None and total and total > 0:
            # New layer, or bytes went backwards → start a fresh speed baseline.
            if digest != self._digest or (self._samples and completed < self._samples[-1][1]):
                self._digest = digest
                self._samples = []
            self._samples.append((t, completed))
            # Keep only the trailing window (but always ≥2 samples to measure a rate).
            cutoff = t - self._window
            while len(self._samples) > 2 and self._samples[0][0] < cutoff:
                self._samples.pop(0)
            if len(self._samples) >= 2:
                dt = self._samples[-1][0] - self._samples[0][0]
                dc = self._samples[-1][1] - self._samples[0][1]
                if dt > 0 and dc >= 0:
                    bps = dc / dt
        # else: non-byte phase (manifest / verifying / writing manifest / success)
        # — leave the bar where it is and show no speed.

        text = self._compose(status, frac, completed, total, bps)
        return {"status": status, "fraction": frac, "completed": completed,
                "total": total, "bps": bps, "text": text}

    @staticmethod
    def _compose(status, frac, completed, total, bps) -> str:
        parts = []
        if status:
            parts.append(status)
        if completed is not None and total:
            parts.append(f"{fmt_size(completed)}/{fmt_size(total)}")
        if frac is not None:
            parts.append(f"{int(frac * 100)}%")
        spd = fmt_speed(bps)
        if spd:
            parts.append(spd)
        return "  ".join(parts) if parts else "…"
