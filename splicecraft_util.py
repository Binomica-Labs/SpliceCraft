"""splicecraft_util — pure cross-cutting helpers (Phase D, layer L0).

Domain-agnostic utilities with no SpliceCraft dependencies (stdlib only): natural
sort, label sanitising, a DataTable cursor-row-key reader, file-extension
predicates, and an export-path splitter. Extracted from the hub so the modal /
screen / widget siblings can import them instead of calling them bare on the hub
(which would be an import cycle). Re-exported by the hub so `sc.<name>` + every
existing call site resolves unchanged.
"""
from __future__ import annotations

import re
import platform
import functools as _functools
import time as _time_mod
from datetime import datetime as _datetime
from typing import Any as _Any
from pathlib import Path

from splicecraft_logging import _log, _log_event


# Snapshot the runtime platform string ONCE at import (INV-36). On some OSes
# `platform.platform()` shells out via `subprocess.run`; tests that monkeypatch
# `subprocess.run` (e.g. the `splicecraft update` command capture) would
# otherwise re-trigger it and crash on the mock. Read this constant rather than
# re-invoking `platform.platform()` in any hot path. Lives in L0 util so the hub
# and the backup sibling share one cached value; re-exported as `sc._RUNTIME_PLATFORM`.
try:
    _RUNTIME_PLATFORM: str = platform.platform()
except Exception:  # pragma: no cover — defensive against weird platforms
    _RUNTIME_PLATFORM = "?"


def _sanitize_path(p: "str | None") -> "Path | None":
    """Resolve a user-supplied path with `~` expansion. Returns None
    on empty / None / non-string input. Doesn't reject any specific
    paths — the user owns their filesystem and the agent runs with
    their permissions, so OS-level file mode controls trust. We just
    normalize so a relative path lands at CWD rather than a surprise
    location, and reject non-string types so a JSON ``{"path": {}}``
    can't smuggle a dict's str() into a Path constructor.

    Refuses ``~user`` syntax (where ``user`` is anything other than
    the empty current user). `Path.expanduser` will happily resolve
    ``~daemon/.bashrc`` to /var/lib/daemon/... and similar — combined
    with `_h_load_file`'s 404 vs 400 distinction this would be a
    user-enumeration oracle for an agent caller. Only ``~`` (current
    user, no trailing identifier) is allowed.
    """
    if not isinstance(p, str) or not p:
        return None
    if p.startswith("~") and len(p) > 1 and p[1] not in ("/", "\\"):
        return None
    return Path(p).expanduser()


_NATURAL_SORT_RE = re.compile(r"(\d+)")

_NATURAL_SORT_KEY_CACHE: "dict[str, tuple]" = {}

_NATURAL_SORT_CACHE_CAP: int = 4096

_FASTA_EXTS: frozenset[str] = frozenset({
    ".fa", ".fasta", ".fna", ".ffn", ".frn", ".fas", ".mpfa", ".faa",
})

_SEQ_ZIP_EXTS: frozenset[str] = frozenset({".zip"})

_CONTROL_CHARS_RE = re.compile(
    # C0 + DEL/C1 + surrogates + line/para separators AND the Unicode
    # DIRECTIONAL format controls - bidi embeddings/overrides (U+202A-202E),
    # isolates (U+2066-2069), LRM/RLM (U+200E/200F), BOM/ZWNBSP (U+FEFF). A
    # name carrying U+202E (RIGHT-TO-LEFT OVERRIDE) visually reorders DataTable
    # cells + `.gb` qualifiers (Trojan-Source-style spoofing) and has no
    # legitimate use in a plasmid / feature name. ZWNJ/ZWJ (U+200C/200D) are
    # deliberately LEFT intact - they legitimately join scripts + emoji.
    r"[\x00-\x1f\x7f-\x9f\u200e\u200f\u2028-\u202e\u2066-\u2069\ufeff"
    r"\ud800-\udfff]+"
)


def _natural_sort_key(s: str) -> tuple:
    """Return a tuple suitable for natural / human-friendly sorting.

    Splits `s` into alternating text and integer runs so that
    `pBin2` sorts before `pBin10`, instead of the lexicographic
    `pBin1 < pBin10 < pBin2 < pBin20`. Tuples carry `(0, str)` for
    text runs and `(1, int)` for digit runs so Python's tuple
    comparison never tries to order an `int` against a `str` (which
    raises in Py3) — that mixed comparison is the gotcha that bites
    the naïve `(int_or_str, ...)` formulation when a name starts
    with a digit (like `"5kb_backbone"` vs `"pBin1"`).

    Cached: identical input strings (the common case across re-sort
    calls when filter keystrokes only change the predicate, not the
    underlying library) reuse the previously computed tuple.
    """
    cached = _NATURAL_SORT_KEY_CACHE.get(s)
    if cached is not None:
        return cached
    out: list[tuple[int, "int | str"]] = []
    for part in _NATURAL_SORT_RE.split(s.lower()):
        if not part:
            continue
        if part.isdigit():
            # `\d+` (str-mode) captures ASCII + Unicode Nd digits, all
            # int()-able. But a *text* run can ALSO be `.isdigit()`-True for
            # No-category chars (superscripts ²³⁴, category No) that `int()`
            # REJECTS — a plasmid / collection / file name like "10²" would
            # otherwise raise ValueError and take down every sort routed
            # through this key (library panel, pickers, file browser). Guard it.
            try:
                out.append((1, int(part)))
            except ValueError:
                out.append((0, part))
        else:
            out.append((0, part))
    result = tuple(out)
    if len(_NATURAL_SORT_KEY_CACHE) >= _NATURAL_SORT_CACHE_CAP:
        # Bounded cache: drop oldest insertion when at cap. dict
        # preserves insertion order, so popping the first key is the
        # cheapest FIFO eviction we can do without a full LRU dance.
        try:
            first = next(iter(_NATURAL_SORT_KEY_CACHE))
            del _NATURAL_SORT_KEY_CACHE[first]
        except StopIteration:
            pass
    _NATURAL_SORT_KEY_CACHE[s] = result
    return result


def _sanitize_label(s: "str | None", *, max_len: int = 200) -> str:
    """Clean a feature label / qualifier value: strip control chars,
    collapse to single line, trim, cap length. Empty / None / non-
    string input → empty string (callers decide the default).

    Type-strict: a dict / list / int payload value is treated as
    "missing" rather than coerced via ``str()`` — coercion would
    silently accept a JSON ``{"name": {"x": 1}}`` and store
    ``"{'x': 1}"`` as the label. Wrong by design.
    """
    if not isinstance(s, str) or not s:
        return ""
    s = _CONTROL_CHARS_RE.sub("", s).strip()
    return s[:max_len]


def _cursor_row_key(table) -> "str | None":
    """Return the value of a DataTable's cursor row key, or None when
    the table is empty or the cursor is out of bounds.

    Centralises the boilerplate `list(t.rows.keys())` + bounds-check
    that was open-coded at ~10 sites (library panel buttons, picker
    modals, primer table). Always pair with the empty-table branch in
    the caller — this helper is read-only.
    """
    if table.row_count == 0:
        return None
    row_keys = list(table.rows.keys())
    if not (0 <= table.cursor_row < len(row_keys)):
        return None
    rk = row_keys[table.cursor_row]
    return rk.value if rk else None


def _is_fasta_path(path) -> bool:
    """True if ``path`` looks like a FASTA file by extension. Accepts
    anything with a ``suffix`` attribute (``pathlib.Path`` or ``DirEntry``)
    or a plain string."""
    try:
        suffix = getattr(path, "suffix", None)
        if suffix is None:
            suffix = Path(str(path)).suffix
    except Exception:
        return False
    return suffix.lower() in _FASTA_EXTS


def _is_seq_zip_path(path) -> bool:
    """True if `path` looks like a sequencing-data archive by extension.
    Currently a `.zip` check; if vendors start shipping `.tar.gz` /
    `.7z` we add them here."""
    try:
        suffix = getattr(path, "suffix", None)
        if suffix is None:
            suffix = Path(str(path)).suffix
    except Exception:
        return False
    return suffix.lower() in _SEQ_ZIP_EXTS


def _split_default_export_path(default_path: str, fallback_filename: str
                                  ) -> "tuple[str, str]":
    """Split `default_path` into ``(parent_dir, filename)`` for the
    "save as" modals. Falls back to `Path.home()` for the dir and
    `fallback_filename` for the filename when either component is
    missing or unreadable. Centralised so each export modal handles
    a missing default the same way."""
    try:
        p = Path(default_path).expanduser() if default_path else None
    except (OSError, ValueError):
        p = None
    if p is not None and p.name:
        parent = p.parent
        filename = p.name
    else:
        parent = Path.home()
        filename = fallback_filename
    try:
        if not parent.is_dir():
            parent = Path.home()
    except OSError:
        parent = Path.home()
    return (str(parent), filename)


# ── Save-failure notifier + display formatters (moved from hub, Phase D) ────
def _notify_save_failure(app, label: str, exc: BaseException,
                          *, severity: str = "error") -> None:
    """Surface a save failure to the user via the app's notify channel.
    Wraps every `_save_*` call site that needs to recover gracefully
    when disk-full / RO-mount / EACCES bubbles out of `_safe_save_json`.

    Per sacred invariant #7, `_safe_save_json` re-raises on failure so
    callers can notify; this helper makes the notify pattern uniform
    rather than each call site composing its own message.

    `app` may be `None` (test contexts that don't run an App). Falls
    back to logging only.

    Emits a structured `save.failed` event so AI parsers of bug-report
    log dumps can correlate the failure target + exception class
    without regex-scraping the human-readable `Save failed for X` log
    line. Sweep #5 — pre-fix, 30+ save sites went through this helper
    silently from a structured-event perspective.
    """
    _log.exception("Save failed for %s", label)
    _log_event("save.failed", target=label,
                exc_type=type(exc).__name__,
                exc_msg=str(exc))
    msg = f"{label} save failed: {exc}"
    if app is None:
        return
    try:
        app.notify(msg, severity=severity, timeout=12)
    except Exception:
        # If notify itself fails (no app, no screen mounted yet) we've
        # already logged via _log.exception — that's enough.
        pass


def _format_identity_pct(pct: "float | int | None", *,
                         decimals: int = 1) -> str:
    """Format an alignment identity percentage for a compact table cell
    such that a value below 100% NEVER renders as ``"100%"``.

    The naïve ``f"{v:.1f}%"`` rounds 99.99% up to ``"100.0%"`` — which
    then reads as a perfect alignment even though `_identity_pct_color`
    (strict ``>= 100.0`` → light-blue) correctly keeps the cell *green*.
    The user flagged exactly that contradiction ("says 100% but it's
    green") for a one-bp mismatch in an 18 kb plasmid. Here, when
    rounding at ``decimals`` places would land on "100", precision is
    escalated one place at a time (capped at 4) until the rendered
    number is strictly < 100, so the same alignment shows e.g.
    ``"99.99%"`` instead. A genuine 100.0 (``n_matches == aligned_cols``)
    still renders the clean ``"100%"`` — no decimals — so a true perfect
    match is visually distinct from a near-perfect one at a glance.

    A value pathologically close to but below 100 (e.g. 99.999999, which
    rounds to "100.0000" even at 4 places) renders ``"<100%"`` rather
    than implying perfection. Non-numeric / None → ``"—"`` (no tier
    implied — matches `_identity_pct_color`'s neutral handling).
    """
    if pct is None:
        return "—"
    try:
        v = float(pct)
    except (TypeError, ValueError):
        return "—"
    # Use the SAME strict ``>= 100.0`` boundary as `_identity_pct_color`
    # so the number and the colour can never disagree about perfection.
    if v >= 100.0:
        return "100%"
    d = max(0, int(decimals))
    for places in range(d, 5):
        s = f"{v:.{places}f}"
        try:
            shown = float(s)
        except ValueError:
            break
        if shown < 100.0:
            return f"{s}%"
    return "<100%"


def _sanitize_plasmid_name(raw: str, *,
                            fallback: str = "assembly",
                            max_len: int = 60) -> str:
    """Clean a user-entered plasmid name for safe storage in the
    library / parts bin and for use as a SeqRecord ``id`` / ``name``.

    Strips control chars (including NUL — would break C-string-style
    handling in downstream tools), trims whitespace, and truncates
    to ``max_len`` chars to keep DataTable rows from blowing the
    row width. Empty / whitespace-only input falls back to
    ``fallback`` so the caller never gets a zero-length string.

    Forbidden character set is conservative: only printable ASCII +
    Unicode letters / digits / a small punctuation set (``_-+. ·:``).
    Colons are kept because the constructor uses them in source
    annotations (``constructor:gid:role``); slashes and backslashes
    are dropped because they look like paths and tools commonly
    interpret them as such.
    """
    if not isinstance(raw, str):
        raw = str(raw or "")
    # Whitespace control chars (\t \n \r \v \f) become spaces FIRST so
    # they don't silently fuse adjacent words after the control-strip
    # pass. e.g. ``"foo\tbar"`` becomes ``"foo bar"``, not ``"foobar"``.
    for ch in "\t\n\r\v\f":
        raw = raw.replace(ch, " ")
    # Drop NUL + remaining C0 control chars (\x00–\x1F + \x7F) — these
    # never belong in a user-facing identifier and break naive
    # C-string handling in some downstream tools.
    cleaned = "".join(
        ch for ch in raw
        if (ord(ch) >= 0x20 and ord(ch) != 0x7F)
    )
    # Drop path-like separators outright — a name like
    # ``../../etc/passwd`` becomes ``etc passwd`` after this filter,
    # so even a malicious agent prompt can't escape into a file path
    # via the library-save flow downstream.
    for ch in "/\\":
        cleaned = cleaned.replace(ch, " ")
    # Normalise whitespace runs to single spaces; lots of tools
    # render multiple spaces awkwardly in TUI tables.
    cleaned = " ".join(cleaned.split())
    cleaned = cleaned.strip()
    if not cleaned:
        return fallback
    if len(cleaned) > max_len:
        cleaned = cleaned[:max_len].rstrip()
    return cleaned or fallback


# ── More pure helpers (moved from hub, Phase D) ─────────────────────────────
def _feature_location(start: int, end: int, total: int, strand):
    """Build the BioPython location for a dict-model feature `(start, end,
    strand)` on a molecule of `total` bp — the INVERSE of `_feat_bounds`.

    `end < start` is the dict convention for an origin-spanning (wrap)
    feature, which becomes a two-part `CompoundLocation`. Everything else
    is a plain `FeatureLocation`. Returns `None` for a degenerate span
    (`end == start`, which would be zero-length or the whole molecule with
    no break point — callers already refuse those).

    **Why this exists — the part ORDER of a wrap on the minus strand.**
    BioPython concatenates a `CompoundLocation`'s parts in STORED order when
    it extracts, and writes `complement(join(...))` by reversing them, so the
    stored order *is* the 5'→3' reading order. For a plus-strand wrap that is
    `[(start, total), (0, end)]` — read the tail, then the head. For a MINUS
    strand it is the other way round: the feature reads from `end - 1` down to
    `start`, so the head piece `(0, end)` comes FIRST. Emitting the ascending
    order for both (which every wrap-building call site used to do) is
    invisible inside SpliceCraft — `_feat_bounds` recovers `(start, end)`
    from either order, and `_translate_cds` works off those numbers — but it
    exports `complement(join(1..end, start+1..total))`, which BioPython,
    SnapGene, Benchling and NCBI all read as `rc(tail) + rc(head)`: the right
    bases in the wrong order, i.e. a CYCLICALLY ROTATED protein for any
    reverse-strand CDS that crosses the origin. Route every wrap through here
    so there is one definition of "which half comes first".

    `strand` is passed through untouched (BioPython wants ±1 / 0 / None);
    SpliceCraft's `2` = double-stranded convention must be mapped by the
    caller before it gets here.
    """
    from Bio.SeqFeature import CompoundLocation, FeatureLocation
    start, end, total = int(start), int(end), int(total)
    if end > start:
        return FeatureLocation(start, end, strand=strand)
    if end == 0 and start < total:
        # Ends exactly ON the origin: a wrap with an empty head half.
        # `CompoundLocation` refuses a zero-length part, and the tail alone
        # already describes every base.
        return FeatureLocation(start, total, strand=strand)
    if end == start:
        return None
    tail = FeatureLocation(start, total, strand=strand)
    head = FeatureLocation(0, end, strand=strand)
    return CompoundLocation([head, tail] if strand == -1 else [tail, head])


def _feat_bounds(feat, total: int) -> "tuple[int, int, int] | None":
    """Wrap-aware extraction of `(start, end, strand)` from a Biopython
    `SeqFeature`. The returned `(start, end)` follows the dict-feature
    convention: `end < start` signals an origin-spanning wrap; otherwise
    `end > start`. Returns `None` if the location has non-integer coords
    (UnknownPosition / BetweenPosition).

    For a `CompoundLocation` of exactly two parts whose outer bounds are
    `[0, ..)` and `[.., total)`, re-encodes as `(tail_start, head_end)`
    so callers can slice with `_slice_circular` and length with `_feat_len`.
    Other compound shapes flatten to outer bounds.

    Callers that read `int(feat.location.start)` / `int(feat.location.end)`
    directly silently flatten wrap features (Biopython returns `min(part.start)`
    for a CompoundLocation), so any code that later does `seq[s:e]` returns
    the BACKBONE GAP rather than the feature. Always route through this
    helper instead. See sacred invariant #9.
    """
    loc = getattr(feat, "location", None)
    if loc is None:
        return None
    # Preserve `loc.strand == None` (BioPython's "no strand info") as
    # 0 (arrowless) rather than coercing to 1. See `PlasmidMap._parse`
    # for the same fix and rationale.
    try:
        _raw_strand = getattr(loc, "strand", None)
        strand = int(_raw_strand) if _raw_strand is not None else 0
    except (TypeError, ValueError):
        strand = 0
    try:
        from Bio.SeqFeature import CompoundLocation
    except ImportError:
        CompoundLocation = None
    if CompoundLocation is not None and isinstance(loc, CompoundLocation):
        try:
            parts = sorted(loc.parts, key=lambda p: int(p.start))
            if (
                total > 0 and len(parts) == 2
                and int(parts[0].start) == 0
                and int(parts[-1].end) == total
                and int(parts[0].end) < int(parts[-1].start)
            ):
                # Origin wrap → (tail_start, head_end) so end < start.
                return int(parts[-1].start), int(parts[0].end), strand
            # Other compound shapes: outer bounds, lossy but oriented.
            return int(parts[0].start), int(parts[-1].end), strand
        except (TypeError, ValueError):
            return None
    try:
        return int(loc.start), int(loc.end), strand
    except (TypeError, ValueError):
        return None


def _name_modal_result(result: "_Any",
                       default_collection: str) -> "tuple[str, str] | None":
    """Normalise a `NamePlasmidModal` dismiss payload into
    ``(name, collection)`` — or ``None`` for cancel / empty.

    Accepts BOTH the collection-mode dict ``{"name", "collection"}`` and a
    bare name ``str`` (legacy callers + direct-dismiss tests), so the
    universal save callbacks stay robust regardless of how the modal was
    dismissed. ``collection`` defaults to ``default_collection`` when the
    payload doesn't carry one."""
    if isinstance(result, dict):
        nm = (result.get("name") or "").strip()
        if not nm:
            return None
        coll = (result.get("collection") or "").strip() or default_collection
        return (nm, coll)
    if isinstance(result, str):
        nm = result.strip()
        return (nm, default_collection) if nm else None
    return None


# ── Export/collection/color/DNA-normalise pure helpers (moved, Phase D) ─────
_IUPAC_NUC_CHARS = frozenset("ACGTUMRWSYKVHDBN")


_IUPAC_NUC_PATTERN = re.compile(r"^[ACGTUMRWSYKVHDBN]+$")


_FASTA_HEADER_PATTERN = re.compile(r"^>[^\n]*\n?", re.MULTILINE)


_SCRUB_PATTERN     = re.compile(r"[\s\d]+")


_HEX_COLOR_RE = re.compile(r"^#[0-9A-Fa-f]{6}$|^#[0-9A-Fa-f]{3}$")


_MAX_COLLECTION_NAME_LEN = 200


def _normalize_dna_for_align(seq: str) -> str:
    """Scrub FASTA header lines + whitespace + digits from ``seq``,
    uppercase, map ``U`` → ``T`` (so RNA pasted into a DNA field
    aligns instead of every base mismatching), and validate that every
    remaining character is an IUPAC nucleotide code
    (``ACGTMRWSYKVHDBN`` after the U-mapping).

    Covers the common bad-paste failure modes upstream of the C-loop:

      * FASTA pasted as-is — ``>name`` line gets stripped along with
        the embedded newlines, leaving only the sequence body.
      * GenBank ``ORIGIN`` block — leading bp position numbers, line
        wraps, spaces between every 10 bp.
      * Protein pasted in a DNA-only field — would have ``EFILPQ`` etc.
        Pre-fix Biopython chewed through these and produced a
        "no alignment" result with no clue why.
      * RNA consensus pasted into a DNA target — every ``U`` used to
        mismatch ``T`` for ~0% identity with no hint. Now mapped to
        ``T`` so the alignment is biologically meaningful.

    Returns the cleaned string. Raises ``ValueError`` on a foreign
    character (with the offending char(s) named in the message so
    the user can find them in the source).
    """
    if not seq:
        return ""
    # Two-pass scrub: FASTA headers first (line-anchored), then
    # whitespace/digits across the whole string. Order matters —
    # stripping `\n` first would let the header text leak into the
    # body and trip the IUPAC check on a leading "MYPLASMID" etc.
    s = _FASTA_HEADER_PATTERN.sub("", seq)
    s = _SCRUB_PATTERN.sub("", s).upper()
    if not s:
        return ""
    # 2026-05-27: silent RNA→DNA. The IUPAC alphabet still admits
    # ``U`` syntactically (frozenset includes it) but the aligner +
    # state classification work on DNA bases — leaving ``U`` would
    # mismatch every paired ``T``. The remap is post-uppercase so
    # both ``u`` and ``U`` in the source collapse to ``T``.
    if "U" in s:
        s = s.replace("U", "T")
    if not _IUPAC_NUC_PATTERN.match(s):
        bad = sorted(set(s) - _IUPAC_NUC_CHARS)
        raise ValueError(
            f"sequence contains non-IUPAC nucleotide character(s): "
            f"{', '.join(repr(c) for c in bad[:6])}"
            f"{' (truncated)' if len(bad) > 6 else ''}"
        )
    return s


def _safe_color_for_picker(raw) -> "str | None":
    """Filter a raw color value to something ColorPickerModal can
    safely consume. The picker assigns the value to
    `styles.background` which raises `StyleValueError` on palette
    references like `color(39)` and on malformed hex strings like
    `#OLDCOL`. Members / library entries CAN carry these
    (canvas `_parse` stamps palette refs at load time, hand-
    edited `.gb` files can carry garbage in `ApEinfo_fwdcolor`),
    so we normalize on the way INTO the picker — None means
    "Auto" / no starting color, and the picker presents its full
    palette fresh.

    Validation: must match `#RGB` or `#RRGGBB`. Anything else
    (palette refs, named colours, mangled hex) → None.

    Sweep #30 (2026-05-26 hardening): defends against the
    `StyleValueError: Invalid color value '...'` crash reported
    when opening per-row color picker on a feature that had a
    palette-ref color or a malformed hex."""
    if not isinstance(raw, str):
        return None
    raw = raw.strip()
    if not raw:
        return None
    if _HEX_COLOR_RE.match(raw) is None:
        return None
    return raw


def _check_export_extension(path: Path, allowed: "tuple[str, ...]",
                              fmt: str) -> "str | None":
    """Enforce an extension whitelist on agent export targets. Without
    this an agent can write `/home/user/.bashrc` as GenBank text (which
    starts with `LOCUS` — not executable but visually hostile / footgun-
    y) or write a `.sh` extension that the user later double-clicks by
    accident. Matches the GUI ExportModal's "save as <FMT>" behaviour
    where the user can't pick an arbitrary extension."""
    suffix = path.suffix.lower()
    if suffix in allowed:
        return None
    return (
        f"refusing to write {fmt} to {path.name!r}: extension must be "
        f"one of {allowed}"
    )


def _normalize_collection_name(s: "str | None") -> "str | None":
    """Trim, strip control chars, cap length, reject blank. Returns
    None on empty input so the caller can 400 the request."""
    name = _sanitize_label(s, max_len=_MAX_COLLECTION_NAME_LEN)
    return name or None


# ── Feature-label + note-sanitise pure helpers (moved from hub, Phase D) ────
_FEAT_LABEL_DISPLAY_MAX = 28


_NOTE_CTRL_RE = re.compile(
    r"[\x00-\x08\x0b-\x1f\x7f-\x9f\u2028\u2029\ud800-\udfff]+"
)


def _feat_label_full(feat) -> str:
    """Canonical, UNtruncated feature label — the value the feature data
    model (`PlasmidMap._parse`) stores and that cloning / export / the
    agent API carry forward. Whitespace is collapsed and a generic
    "unnamed protein product; <gene>" note is reduced to <gene>.

    DISPLAY callers truncate to their own width (the map + seq-panel via
    `_feat_decorated_label`, the sidebar via its own `[:14]`); `_feat_label`
    is the `_FEAT_LABEL_DISPLAY_MAX`-char display wrapper. Baking that cap
    in here used to corrupt persisted data — a lifted operon's six genes
    ("…; luxC" … "…; luxG") all collapsed to the identical
    "unnamed protein product; lux" because the gene letter sits at char 29,
    past the cap, so the cloned plasmid stored six indistinguishable CDS
    bars (user-reported)."""
    for q in ("label", "gene", "product", "standard_name", "note", "bound_moiety"):
        if q in feat.qualifiers:
            v = feat.qualifiers[q]
            # Biopython normally wraps qualifier values in a 1+ element
            # list, but malformed GenBank files can produce empty lists
            # or bare strings. Guard both.
            if isinstance(v, list):
                if not v:
                    continue
                s = v[0]
            else:
                s = v
            if not isinstance(s, str):
                continue
            # Collapse whitespace characters (newline, tab, vertical tab)
            # into single spaces so a multi-line /note="…" qualifier
            # doesn't break the sidebar row or clobber the map label.
            # Then strip and fall through if the result is empty.
            s = " ".join(s.split())
            if s:
                return _surface_placeholder_gene(s)
    return feat.type


def _feat_label(feat) -> str:
    """Display label: `_feat_label_full` truncated to the map / sidebar
    width cap (`_FEAT_LABEL_DISPLAY_MAX`). Tested at 28 chars by
    `test_real_label_truncates_at_28`."""
    return _feat_label_full(feat)[:_FEAT_LABEL_DISPLAY_MAX]


def _sanitize_note(s: "str | None", *, max_len: int = 8000) -> str:
    """Clean a feature ``/note`` body: strip dangerous control bytes
    (preserves `\\t` and `\\n` so multi-paragraph Markdown survives),
    cap at `max_len` characters, trim trailing whitespace.

    Type-strict like `_sanitize_label` — a JSON dict / int payload
    becomes empty rather than `str()`-coerced. Empty / None / non-
    string input → empty string. The 8 KB cap matches typical
    GenBank ``/note`` conventions and prevents adversarial / accidental
    pasted blobs from bloating `.gb` exports or stalling the Markdown
    parser. Callers split on blank-line paragraphs after sanitizing,
    so the cap applies to the combined note body, not per-line.
    """
    if not isinstance(s, str) or not s:
        return ""
    s = _NOTE_CTRL_RE.sub("", s).rstrip()
    return s[:max_len]


# ── Placeholder-gene helper (completes the _feat_label closure, Phase D) ────
_GENERIC_PRODUCT_PLACEHOLDERS: tuple[str, ...] = (
    "unnamed protein product",
)


def _surface_placeholder_gene(s: str) -> str:
    """Reduce a generic ``"unnamed protein product; luxC"`` placeholder to
    its trailing identifier (``"luxC"``). A bare placeholder with nothing
    meaningful after it is returned unchanged. Case-insensitive match on
    the placeholder; the surfaced remainder keeps its original case."""
    low = s.lower()
    for ph in _GENERIC_PRODUCT_PLACEHOLDERS:
        if low.startswith(ph):
            rest = s[len(ph):].lstrip(" \t;:,-/|")
            return rest if rest else s
    return s


# ── Primer-Tm + single-record-pick pure helpers (moved, Phase D) ────────────
@_functools.lru_cache(maxsize=512)
def _primer_tm_safe(seq: str) -> "float | None":
    """Memoized, defensive primer3 Tm calculation. Returns None if
    primer3 is unavailable, the seq is outside the calc's useful
    range (5..200 bp), or the underlying call raises. Caller
    renders a `—` placeholder for None.

    Cached because `PrimerEditModal._seq_changed` repaints the
    stats line on every keystroke — a user typing/backspacing
    re-hits the same intermediate strings, and the nearest-
    neighbor thermodynamics is the slow part of the path.
    """
    s = (seq or "").upper()
    if not (5 <= len(s) <= 200):
        return None
    try:
        import primer3
        return float(primer3.calc_tm(s))
    except (ImportError, OSError, ValueError, RuntimeError, TypeError):
        return None


def _pick_single_record(records: list, source: str):
    """Given a list of SeqRecords, return the single one if there's exactly
    one, else raise ValueError with a user-friendly message. Used by both
    NCBI fetch and file load so the error text is consistent.
    """
    if not records:
        raise ValueError(
            f"{source} contained no GenBank records. Is it a valid .gb/.gbk file?"
        )
    if len(records) > 1:
        ids = ", ".join(r.id for r in records[:3])
        more = f" (and {len(records) - 3} more)" if len(records) > 3 else ""
        raise ValueError(
            f"{source} contains {len(records)} records — SpliceCraft loads "
            f"one plasmid at a time. Split the file or extract a single "
            f"record first (found: {ids}{more})."
        )
    return records[0]


# ── Identity-color + path-scrub pure helpers (moved from hub, Phase D) ──────
def _scrub_path(text: str) -> str:
    """Replace home-directory paths with `~` so a snapshot/bundle is
    safe to share without exposing the user's username. Best-effort:
    if `Path.home()` is unset (rare), the input is returned as-is.

    Always returns a string. Multi-platform — handles Linux/Mac
    (`/home/<user>`, `/Users/<user>`) and Windows
    (`C:\\Users\\<user>`, `%USERPROFILE%`).
    """
    if not isinstance(text, str):
        return str(text)
    try:
        home = str(Path.home())
    except (RuntimeError, OSError):
        return text
    if not home:
        return text
    out = text.replace(home, "~")
    # Cross-OS: if the user runs the bundle command on a different
    # filesystem layout from where the log was written (rare but
    # possible: WSL / network home dirs), the literal /home/<user>
    # pattern may still appear. Fall back to a regex scrub.
    try:
        out = re.sub(r"(?:/home|/Users)/[A-Za-z0-9_.\-]+", "~", out)
        # Windows: handle both forward and back slashes.
        out = re.sub(
            r"[A-Za-z]:[\\/]Users[\\/][A-Za-z0-9_.\-]+",
            "~",
            out,
        )
    except re.error:
        pass
    return out


def _identity_pct_color(pct: "float | int | None") -> str:
    """Map an alignment identity percentage to a Rich color name for
    table cells. Tiers picked 2026-05-27 to match the user's
    sequencing-QC color grading:

      * **100.0% (strict)** → ``bright_cyan`` (light blue). No
        rounding-up — a 99.999% identity does NOT promote to light
        blue; only a literal 100.0 from the pairwise aligner counts.
      * **>= 90%** → ``green``
      * **>= 80%** → ``yellow``
      * **>= 51%** → ``dark_orange``
      * **>= 11%** → ``red``
      * **<= 10%** → ``grey50`` (gray; "barely aligned" indicator)
      * Non-numeric / None → ``white`` (neutral; surface raw to avoid
        accidentally implying a quality tier the value doesn't carry)

    Used by `AlignmentManagerModal._repopulate` and
    `VerificationReportModal._add_row`. Pure function — testable
    without spinning up a screen.
    """
    if pct is None:
        return "white"
    try:
        v = float(pct)
    except (TypeError, ValueError):
        return "white"
    # Use `>= 100.0` (strict). pairwise-align computes
    # identity_pct = 100.0 * n_matches / aligned_len, which lands
    # on exactly 100.0 only when n_matches == aligned_len; any
    # mismatch / gap pushes below 100.0 even by a fractional bp,
    # so this guard isn't subject to float drift in practice.
    if v >= 100.0:
        return "bright_cyan"
    if v >= 90.0:
        return "green"
    if v >= 80.0:
        return "yellow"
    if v >= 51.0:
        return "dark_orange"
    if v >= 11.0:
        return "red"
    return "grey50"


# ── Status/fuzzy/gel/feat sanitisers (moved from hub, Phase D) ──────────────
_GEL_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._\-]{0,63}$")


_PLASMID_STATUS_VALUES: tuple[str, ...] = (
    "DESIGNING", "CLONING", "SEQUENCING", "VERIFIED", "ERROR",
)


def _fuzzy_match(query: str, name: str) -> bool:
    """Case-insensitive subsequence match: True if every char of `query`
    appears in `name` in order (not necessarily contiguous). Empty query
    matches everything. Used by LibraryPanel's search filter and by
    `_search_collections_library` (cross-collection search).

    Performance: O(len(query) * len(name)) worst case via repeated
    `str.find`. Early-rejects when ``len(query) > len(name)`` (no
    subsequence can be longer than its container) — saves the
    full lower() + scan on huge libraries where most names are
    shorter than a typical fuzzy query. The lower() call is the hot
    spot for very long names; the early reject runs on the original
    strings before any allocation.
    """
    if not query:
        return True
    if not name:
        return False
    # Subsequence membership requires `len(query) <= len(name)`. Skip
    # the lower() + scan entirely when impossible. This is the hot
    # exit on huge libraries where typical names are short and the
    # user's filter has 5+ chars.
    if len(query) > len(name):
        return False
    q, n = query.lower(), name.lower()
    i = 0
    for ch in q:
        i = n.find(ch, i)
        if i < 0:
            return False
        i += 1
    return True


def _sanitize_gel_id(raw: object) -> "str | None":
    """Mirror `_sanitize_experiment_id` — a gel id is filesystem-safe
    (it shows up in tag tokens + JSON path keys + log events).
    Rejects empty / NUL / `..` / `/` / `\\` / >64 chars."""
    if not isinstance(raw, str) or not raw:
        return None
    if "\x00" in raw or ".." in raw or "/" in raw or "\\" in raw:
        return None
    if not _GEL_ID_RE.match(raw):
        return None
    return raw


def _sanitize_feat_type(s: "str | None", *,
                         max_len: int = 50,
                         default: str = "misc_feature") -> str:
    """Clean a GenBank feature-type string. Empty / None / non-string
    input → ``default``. INSDC types are short alpha words with
    underscore — non-strict here (we accept anything printable) since
    users may legitimately add custom types not in the curated list."""
    if not isinstance(s, str) or not s:
        return default
    s = _CONTROL_CHARS_RE.sub("", s).strip()
    return (s or default)[:max_len]


def _sanitize_plasmid_status(s) -> str:
    """Type-strict acceptance of the four canonical workflow status
    strings. Anything else — including None, empty string, dicts,
    leading/trailing whitespace variants of recognised values, or
    case mismatches — collapses to empty (the "no status" sentinel).

    Strict matching is intentional: a hand-edited library JSON with
    `"status": "Designing"` (mixed case) silently degrades to no
    status rather than being normalised, because we want the on-
    disk value to round-trip through this function exactly. Callers
    that want lenient parsing should upper-case + strip themselves
    before passing through.
    """
    if not isinstance(s, str):
        return ""
    return s if s in _PLASMID_STATUS_VALUES else ""


# ── Group-member validation (moved from hub, Phase D) ──────────────────────
_MAX_GROUP_MEMBERS = 64


_MAX_GROUP_LABEL_LEN = 200


_MAX_GROUP_COLOR_LEN = 32


def _validate_group_members(
    members: list,
    sequence_len: int,
) -> "list[dict]":
    """Normalise + validate a group entry's `members` list. Returns
    a fresh list of well-formed member dicts in `rel_start` order
    (stable sort, so duplicate-start members preserve input order).

    Raises `ValueError` for unrecoverable shape errors so the save
    / load paths surface a clear error to the user rather than
    persisting a half-broken entry that later trips the annotate
    path with a NoneType / KeyError deep in the rendering call
    stack.

    Validation rules:
      * `members` must be a non-empty list of dicts.
      * Each member needs int `rel_start`, `rel_end` with
        `0 <= rel_start < rel_end <= sequence_len`. The closed-open
        half-interval matches the existing single-feature
        convention used by `_annotate_with_feature_impl`.
      * `feature_type` defaults to `"misc_feature"` if missing /
        empty (matches `_annotate_with_feature_impl`).
      * `strand` defaults to 1, clamped to `{-1, 0, 1, 2}`.
      * `color` defaults to `None` (Auto / palette fallback);
        non-hex strings are passed through unmodified for
        forward-compat (Rich may accept named colours).
      * `qualifiers` defaults to `{}`; non-dict values are
        replaced with `{}` rather than raising.
      * `label` defaults to `""`; non-str values are stringified.
      * `description` defaults to `""`; same coercion.

    Members CAN overlap (rare but legitimate — e.g. a parent CDS
    with a sub-feature like an active site marker) and CAN leave
    gaps (= unannotated bases inside the group's sequence). The
    validator deliberately does NOT enforce tiling so the user
    can model real biology."""
    if not isinstance(members, list) or not members:
        raise ValueError(
            "group entry must have a non-empty `members` list"
        )
    if len(members) > _MAX_GROUP_MEMBERS:
        # Hardening (2026-05-26): defends against pasted JSON /
        # textarea content with thousands of rows. A real Golden
        # Gate adapter has 4 members; biology rarely exceeds 8.
        # `_MAX_GROUP_MEMBERS = 64` leaves comfortable headroom.
        raise ValueError(
            f"too many members ({len(members)} > "
            f"{_MAX_GROUP_MEMBERS}) — split into multiple groups "
            "or raise `_MAX_GROUP_MEMBERS` if you really need "
            "this many sub-features"
        )
    if not isinstance(sequence_len, int) or sequence_len < 1:
        raise ValueError(
            f"group sequence_len must be a positive int "
            f"(got {sequence_len!r})"
        )
    out: list[dict] = []
    for i, m in enumerate(members):
        if not isinstance(m, dict):
            raise ValueError(
                f"member {i} must be a dict (got {type(m).__name__})"
            )
        rs_raw = m.get("rel_start")
        re_raw = m.get("rel_end")
        if rs_raw is None or re_raw is None:
            raise ValueError(
                f"member {i}: rel_start / rel_end must be ints "
                f"(missing)"
            )
        try:
            rs = int(rs_raw)
            re = int(re_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"member {i}: rel_start / rel_end must be ints "
                f"({exc})"
            )
        if not (0 <= rs < re <= sequence_len):
            raise ValueError(
                f"member {i}: rel_start={rs}, rel_end={re} not in "
                f"[0, {sequence_len}] half-open"
            )
        raw_strand = m.get("strand", 1)
        try:
            strand = int(raw_strand)
        except (TypeError, ValueError):
            strand = 1
        if strand not in (-1, 0, 1, 2):
            strand = 1
        # Sweep #29 hardening (2026-05-26): every user-supplied
        # string field is scrubbed of control chars + length-
        # capped. Paste of a binary blob or ANSI-escape-laden
        # text can't smuggle escape sequences into Rich Text /
        # `.gb` qualifier values via the group library route.
        raw_ftype = m.get("feature_type") or ""
        if not isinstance(raw_ftype, str):
            raw_ftype = str(raw_ftype) if raw_ftype is not None else ""
        ftype = _sanitize_feat_type(raw_ftype) or "misc_feature"
        color = m.get("color")
        if not isinstance(color, str):
            color = None
        else:
            color = _CONTROL_CHARS_RE.sub("", color).strip()
            if not color:
                color = None
            elif len(color) > _MAX_GROUP_COLOR_LEN:
                # Suspiciously long "colour" → drop, palette
                # fallback applies at render time.
                color = None
        raw_label = m.get("label")
        if not isinstance(raw_label, str):
            raw_label = str(raw_label) if raw_label is not None else ""
        label = _sanitize_label(raw_label, max_len=_MAX_GROUP_LABEL_LEN)
        raw_desc = m.get("description")
        if not isinstance(raw_desc, str):
            raw_desc = str(raw_desc) if raw_desc is not None else ""
        desc = _sanitize_note(raw_desc)
        quals = m.get("qualifiers")
        if not isinstance(quals, dict):
            quals = {}
        # Sanitize qualifier values too — they round-trip to
        # `.gb` exports, so a control char in a value would
        # corrupt the file. Keys are not user-input in practice
        # (they come from GenBank spec / the library JSON)
        # but we cap their length defensively.
        clean_quals: dict = {}
        for qk, qv in quals.items():
            if not isinstance(qk, str):
                continue
            qk_clean = _CONTROL_CHARS_RE.sub("", qk).strip()
            if not qk_clean or len(qk_clean) > 64:
                continue
            if isinstance(qv, (list, tuple)):
                clean_vals = []
                for v in qv:
                    if isinstance(v, str):
                        cv = _CONTROL_CHARS_RE.sub("", v)
                        if len(cv) > 8192:
                            cv = cv[:8192]
                        clean_vals.append(cv)
                clean_quals[qk_clean] = clean_vals
            elif isinstance(qv, str):
                cv = _CONTROL_CHARS_RE.sub("", qv)
                if len(cv) > 8192:
                    cv = cv[:8192]
                clean_quals[qk_clean] = [cv]
        out.append({
            "rel_start":    rs,
            "rel_end":      re,
            "feature_type": ftype,
            "label":        label,
            "color":        color,
            "strand":       strand,
            "qualifiers":   clean_quals,
            "description":  desc,
        })
    out.sort(key=lambda m: (m["rel_start"], m["rel_end"]))
    return out


# ── Markdown-escape + circular-topology pure helpers (moved, Phase D) ───────
def _gb_text_is_circular(gb_text: "str | None") -> bool:
    """Cheap topology read from a GenBank record's LOCUS line, for the
    origin-history helpers. Returns False only when the LOCUS line
    explicitly says ``linear`` (PCR amplicons, synthesis fragments);
    defaults to circular for the common plasmid case + unmarked text."""
    if not gb_text:
        return True
    return "linear" not in gb_text.split("\n", 1)[0].lower()


def _esc_md(s: str) -> str:
    """Local helper — escape Rich markup metachars so a hostile / odd
    file path or error string can't inject styling into Static updates
    in this modal. Equivalent to `rich.markup.escape`."""
    from rich.markup import escape
    return escape(s)


def _now_iso() -> str:
    """ISO-8601 timestamp with local-timezone offset, e.g.
    `2026-05-18T14:30:00-07:00`. Used for entry `created_at` /
    `updated_at`. Timezone-aware so a user roaming across timezones
    still gets unambiguous timestamps."""
    return _datetime.now().astimezone().isoformat(timespec="seconds")


# ── Feature-type default colours (relocated from widgets, fileio .dna prereq) ──
# Pure data: feature_type -> hex colour. Read by widgets `_resolve_feature_color`,
# a feature-colour modal, and the .dna writer's feature-packet builder. Lives at
# L0 so the L2 fileio .dna codec can reach it without an upward widgets import.
_DEFAULT_TYPE_COLORS: dict[str, str] = {
    "CDS":             "#FFA500",
    "gene":            "#FFD700",
    "mRNA":            "#FFA07A",
    "tRNA":            "#FF69B4",
    "rRNA":            "#FF1493",
    "ncRNA":           "#DA70D6",
    "misc_RNA":        "#BA55D3",
    "promoter":        "#00CED1",
    "terminator":      "#DC143C",
    "RBS":             "#00FF7F",
    "polyA_signal":    "#FF6347",
    "regulatory":      "#7FFFD4",
    "5'UTR":           "#87CEEB",
    "3'UTR":           "#4682B4",
    "intron":          "#A9A9A9",
    "exon":            "#90EE90",
    "operon":          "#DDA0DD",
    "primer_bind":     "#00BFFF",
    "protein_bind":    "#F08080",
    "misc_binding":    "#FF8C00",
    "repeat_region":   "#CD853F",
    "LTR":             "#8B4513",
    "mobile_element":  "#8B008B",
    "rep_origin":      "#9370DB",
    "oriT":            "#BA55D3",
    "sig_peptide":     "#ADFF2F",
    "mat_peptide":     "#9ACD32",
    "transit_peptide": "#7CFC00",
    "propeptide":      "#6B8E23",
    "misc_feature":    "#20B2AA",
    "misc_recomb":     "#48D1CC",
    "stem_loop":       "#FF4500",
    "variation":       "#800080",
}


# ── XML-security parse primitive (relocated from hub, [PIT-19]) ──────────────
# Pure stdlib (lazy io + xml.etree). Defangs billion-laughs / XXE / deep-nesting.
# Shared by NCBI fetch, BLAST XML, .dna history XML — an L0 security leaf.
def _safe_xml_parse(xml_data: str, *, allow_dtd: bool = False):
    """Parse XML with defense against billion-laughs / XXE tricks.

    Python's stdlib ET (expat) already refuses to fetch external entities
    since 3.7.1, so the remaining attack surface is DTD-declared entity
    expansion. Pre-fix we substring-matched `<!doctype` / `<!entity` only
    in the FIRST 4096 chars — exploitable with a long whitespace /
    comment prefix that pushes the DOCTYPE past the window.

    Replaced with a streaming prologue scan: skip leading whitespace +
    comments + processing instructions and inspect the first non-trivial
    declaration. If it's a DOCTYPE / ENTITY, refuse before handing to
    expat. NCBI / Kazusa / .dna history XML have no legitimate DTD, so
    this never false-positives.

    ``allow_dtd=True`` opts a caller into permitting an **external** DTD
    reference — e.g. NCBI BLAST's ``FORMAT_TYPE=XML`` output, which opens
    with ``<!DOCTYPE BlastOutput PUBLIC ... NCBI_BlastOutput.dtd>``. A
    standalone ``<!ENTITY>`` and any DOCTYPE carrying an internal subset
    (``[ … ]`` — where billion-laughs entity definitions live) are still
    refused. expat never fetches the external DTD (no network since
    3.7.1) and there are no entities to expand, so this stays XXE-safe.
    Default ``False`` keeps every existing caller's strict behaviour.
    """
    import io
    import xml.etree.ElementTree as ET
    # Streaming prologue scan: advance past whitespace, comments, and
    # processing instructions to find the first non-trivial declaration.
    i = 0
    n = len(xml_data)
    while i < n:
        c = xml_data[i]
        if c in " \t\r\n":
            i += 1
            continue
        # Comment: <!-- ... -->
        if xml_data.startswith("<!--", i):
            end = xml_data.find("-->", i + 4)
            if end == -1:
                # Unterminated comment — defer to expat for the parse
                # error so the user sees the proper diagnostic.
                break
            i = end + 3
            continue
        # Processing instruction: <?xml version="1.0"?> etc.
        if xml_data.startswith("<?", i):
            end = xml_data.find("?>", i + 2)
            if end == -1:
                break
            i = end + 2
            continue
        # DOCTYPE / ENTITY declaration — refuse (or, with allow_dtd, permit
        # only an external-DTD reference). Case-insensitive check on the
        # next ~16 chars (case-folding the WHOLE document is expensive on
        # multi-MB inputs).
        head = xml_data[i:i + 16].lower()
        if head.startswith("<!entity"):
            raise ET.ParseError(
                "XML contains a standalone ENTITY declaration — refusing"
            )
        if head.startswith("<!doctype"):
            if not allow_dtd:
                raise ET.ParseError(
                    "XML contains DTD/ENTITY — refusing to parse"
                )
            # Permit an external-DTD reference but refuse any internal
            # subset (`[ … ]`), which is where entity-expansion attacks live.
            # Find the DOCTYPE's terminating `>` QUOTE-AWARELY: a SYSTEM/PUBLIC
            # literal can legally contain a `>` (e.g. `SYSTEM "a>b"`), so a naive
            # first-`>` scan would stop INSIDE the literal and miss a trailing
            # `[ … ]` subset — a billion-laughs bypass. Walk char-by-char,
            # skipping quoted spans, and stop at the first UNquoted `[` (subset →
            # refuse) or `>` (clean external ref → skip past it).
            j = i + len("<!doctype")
            quote = ""
            close = -1
            saw_bracket = False
            n_xml = len(xml_data)
            while j < n_xml:
                ch = xml_data[j]
                if quote:
                    if ch == quote:
                        quote = ""
                elif ch in ("'", '"'):
                    quote = ch
                elif ch == "[":
                    saw_bracket = True
                    break
                elif ch == ">":
                    close = j
                    break
                j += 1
            if saw_bracket:
                raise ET.ParseError(
                    "XML DOCTYPE has an internal subset — refusing to parse"
                )
            if close == -1:
                break  # malformed / unterminated — let expat surface it
            i = close + 1
            continue
        # Reached a normal start tag — safe to hand off to expat.
        break
    # 2026-05-27 (audit-3 M6): cap nesting depth via iterparse so a
    # 10k-deep `<Node><Node>...</Node></Node>` chain can't blow the
    # Python stack in downstream consumers (HistoryTree walkers,
    # qualifier readers). 256 levels is more than any legitimate
    # plasmid annotation tree ever produces.
    #
    # Parse ONCE: iterparse builds the tree while we check depth, so we capture
    # and RETURN its root instead of re-running `ET.fromstring` over the same
    # document. The old form parsed twice — a full iterparse tree for the depth
    # scan, then a second full `fromstring` build — doubling CPU and peak memory
    # on multi-MB XML (the amplifier behind the .dna packet-parse OOM).
    _MAX_DEPTH = 256
    depth = 0
    root = None
    for event, elem in ET.iterparse(
            io.StringIO(xml_data), events=("start", "end"),
    ):
        if event == "start":
            depth += 1
            if depth > _MAX_DEPTH:
                raise ET.ParseError(
                    f"XML nesting exceeds {_MAX_DEPTH} levels "
                    f"(unbounded recursion guard)"
                )
        else:  # "end"
            depth -= 1
            if depth == 0:
                root = elem   # the document root's end event fires last
    if root is None:
        # Empty / element-less input — defer to fromstring for the canonical
        # "no element found" ParseError so callers see the same diagnostic.
        return ET.fromstring(xml_data)
    return root


# ── Windows-reserved-filename guard (relocated from hub, blob/export prereq) ──
# Shared by every sanitiser that emits an on-disk filename (_dna_sidecar_path
# in fileio, _safe_export_filename + _safe_snapshot_token in the hub) — an L0
# filename-safety leaf so the L2 blob store can reach it without importing up.
# Windows reserved device names (case-insensitive). NTFS refuses to
# open any file whose stem matches one of these, so every sanitiser
# that produces an on-disk filename must rewrite a matching stem.
# Hoisted to module scope so the three sanitisers (`_dna_sidecar_path`,
# `_safe_export_filename`, `_safe_snapshot_token`) share one set and a
# new sanitiser added later can't accidentally drift to a stale list.
_WIN_RESERVED_FILENAMES: frozenset = frozenset({
    "con", "prn", "aux", "nul",
    *(f"com{i}" for i in range(1, 10)),
    *(f"lpt{i}" for i in range(1, 10)),
})


def _is_windows_reserved_stem(name: str) -> bool:
    """True if `name`'s first dot-separated segment is a Windows
    reserved device name (CON, PRN, AUX, NUL, COM1-9, LPT1-9). Case-
    insensitive; matches the Win32 namespace rules."""
    if not name:
        return False
    return name.split(".")[0].lower() in _WIN_RESERVED_FILENAMES


def _strip_fasta_headers(text: str) -> str:
    """Drop FASTA header lines (lines whose first non-whitespace char is
    ``>``) so a paste of a ``>id\\nATG…`` blob becomes just the body.
    Preserves the rest of the text — `_detect_query_program` then
    handles whitespace + alphabet.

    Tolerant to leading whitespace on the header so a copy-pasted
    FASTA from a wrapped editor or e-mail still gets cleaned. A line
    that's just ``>`` with nothing after it is also dropped (defends
    against misformatted multi-FASTA pastes).
    """
    if ">" not in (text or ""):
        return text or ""
    out = []
    for line in (text or "").splitlines():
        if line.lstrip().startswith(">"):
            continue
        out.append(line)
    return "\n".join(out)


# ── Single time source (INV-78, sweep #27) ────────────────────────────
# Every wall-clock and monotonic-clock read in the codebase routes through
# these two helpers (relocated here from the hub during the modularization so
# the siblings can import them at L0 instead of calling them bare on the hub).
# Reasons: (1) a monkeypatch on `_now` / `_monotonic` covers every callsite at
# once; (2) `_now()` is always tz-aware, avoiding silent DST off-by-hours bugs;
# (3) a future "freeze the clock for this transaction" need lands in one place.
def _now() -> _datetime:
    """Current wall-clock as a tz-aware datetime in the user's local
    zone. Use this everywhere instead of `datetime.now()`."""
    return _datetime.now().astimezone()


def _monotonic() -> float:
    """Monotonic clock reading in seconds. Use for elapsed-time
    measurements; `_now()` for absolute timestamps."""
    return _time_mod.monotonic()


# ── Citation metadata (Zenodo DOI) ────────────────────────────────────
# SpliceCraft is archived on Zenodo: every GitHub Release is snapshotted and
# stamped with its own version DOI, and one CONCEPT DOI (below) always resolves
# to the newest archived version. The concept DOI is what a paper should cite
# when the exact version doesn't matter; `--citation` prints the running
# version alongside it so a methods section can be precise.
#
# This is the CONCEPT DOI, minted with the first archived release (v1.2.53,
# 2026-08-15) and shared by every version since; the per-version DOIs sit one
# id higher and upward (v1.2.53 is 10.5281/zenodo.21960401) and are listed on
# the record. The formatters below fall back to the GitHub URL when this is
# empty rather than printing a placeholder that looks like a real identifier.
# tests/test_citation.py enforces that README.md, CITATION.cff, and
# docs/citation.md all carry the SAME string (all-or-nothing — a half-wired
# DOI can't ship).
_ZENODO_CONCEPT_DOI = "10.5281/zenodo.21960400"

# GitHub's numeric repository id for Binomica-Labs/SpliceCraft. Zenodo's badge
# endpoints are keyed on it (`zenodo.org/badge/<id>.svg` renders the badge,
# `zenodo.org/badge/latestdoi/<id>` 302s to the newest version DOI), which is
# how the README badge can be wired up before any DOI exists.
_ZENODO_GITHUB_REPO_ID = "1190666059"

_CITATION_TITLE = (
    "SpliceCraft: a terminal-native plasmid workbench for molecular cloning"
)
# BibTeX styles that title-case (plain / unsrt / abbrv run the title through
# `"t" change.case$`) lowercase every letter but the first, rendering
# "SpliceCraft" as "Splicecraft" in the bibliography. Wrapping the program name
# in its own brace group pins the capitalisation; the braces themselves do not
# render. Derived from the title above so the two can never drift.
_CITATION_TITLE_BIBTEX = _CITATION_TITLE.replace("SpliceCraft", "{SpliceCraft}", 1)

_CITATION_AUTHOR_APA = "Cocioba, S."
_CITATION_AUTHOR_BIBTEX = "Cocioba, Sebastian"
_CITATION_REPO_URL = "https://github.com/Binomica-Labs/SpliceCraft"


def _citation_doi_url() -> str:
    """Resolver URL for the concept DOI, or "" when none is minted yet."""
    return f"https://doi.org/{_ZENODO_CONCEPT_DOI}" if _ZENODO_CONCEPT_DOI else ""


def _citation_year(released: str) -> str:
    """Publication year for a citation, taken from a ``YYYY-MM-DD`` release
    date. Falls back to the current year when the stamp is missing or
    malformed — a citation with a wrong-ish year still resolves; one with a
    crashed formatter doesn't."""
    # `str(...)` so a caller handing over a `date` (or anything else with a
    # sane `__str__`) degrades to the fallback instead of raising mid-citation.
    # ASCII-only: `"\u0662\u0660\u0662\u0666".isdigit()` is True, and a year in
    # Arabic-Indic digits would sail into a bibliography unnoticed.
    head = str(released or "").strip()[:4]
    if len(head) == 4 and head.isascii() and head.isdigit():
        return head
    return str(_now().year)


def _citation_apa(version: str, released: str) -> str:
    """One-line APA-7 software reference for *version*."""
    where = _citation_doi_url() or _CITATION_REPO_URL
    publisher = "Zenodo. " if _ZENODO_CONCEPT_DOI else ""
    return (
        f"{_CITATION_AUTHOR_APA} ({_citation_year(released)}). "
        f"{_CITATION_TITLE} (Version {version}) [Computer software]. "
        f"{publisher}{where}"
    )


def _citation_bibtex(version: str, released: str) -> str:
    """BibTeX ``@software`` entry for *version*."""
    lines = [
        "@software{splicecraft,",
        f"  author    = {{{_CITATION_AUTHOR_BIBTEX}}},",
        f"  title     = {{{_CITATION_TITLE_BIBTEX}}},",
        f"  version   = {{{version}}},",
        f"  year      = {{{_citation_year(released)}}},",
    ]
    if _ZENODO_CONCEPT_DOI:
        lines.append("  publisher = {Zenodo},")
        lines.append(f"  doi       = {{{_ZENODO_CONCEPT_DOI}}},")
        lines.append(f"  url       = {{{_citation_doi_url()}}}")
    else:
        lines.append("  publisher = {GitHub},")
        lines.append(f"  url       = {{{_CITATION_REPO_URL}}}")
    lines.append("}")
    return "\n".join(lines)


def _citation_text(version: str, released: str) -> str:
    """The full block printed by ``splicecraft --citation`` and shown in the
    in-app help: a plain-text APA reference plus a BibTeX entry, both pinned to
    the version that is actually running."""
    out = [
        f"SpliceCraft {version} — how to cite",
        "",
        _citation_apa(version, released),
        "",
        _citation_bibtex(version, released),
        "",
    ]
    if _ZENODO_CONCEPT_DOI:
        out.append(
            f"The DOI above ({_ZENODO_CONCEPT_DOI}) is the concept DOI and always "
            f"resolves to the\nnewest release. To cite this exact version instead, "
            f"open the record and pick\nv{version} from its version list."
        )
    else:
        out.append(
            "A Zenodo DOI is being minted for the next release; until then cite "
            "the\nrepository URL above."
        )
    return "\n".join(out)
