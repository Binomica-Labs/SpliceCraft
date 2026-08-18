"""splicecraft_record — GenBank <-> SeqRecord serialization (Phase D, layer L1).

The hub's record serialization core, extracted so the modal/screen siblings that
parse/serialise records can import it instead of calling it bare on the hub. This
is MISSION-CRITICAL (INV-98: the GenBank LOCUS must never carry the human/display
name; the .dna + library round-trips depend on byte-faithful parse/serialise), so
it moved as one self-contained unit with its own LRU parse cache. Biopython is
imported lazily inside the two entry points (`from Bio import SeqIO`) exactly as
in the hub. The SpliceCraft version stamped into the GenBank COMMENT is read from
`_state._sc_version` (the hub sets it from `__version__`, which must stay in
splicecraft.py for release.py's bump regex). Re-exported by the hub so `sc.<name>`
+ every existing call site resolves unchanged.
"""
from __future__ import annotations

import re
import threading
from collections import OrderedDict as _OD
from copy import copy as _shallow_copy
from io import StringIO

import splicecraft_state as _state


_SC_STRAND_QUAL = "SpliceCraft_strand"

_GB_TEXT_MAX_BYTES = 64 * 1024 * 1024

_GB_PARSE_CACHE: "_OD" = _OD()

_GB_PARSE_CACHE_MAX = 16

_GB_PARSE_CACHE_LOCK = threading.RLock()

_GB_LOCUS_NAME_MAX = 28  # NCBI relaxed LOCUS name length (spec is 16)


def _normalize_primer_seq(raw: object) -> str:
    """Canonical primer-sequence sanitiser: whitespace-strip + uppercase.
    THE single function every primer string routes through before it is
    rendered, stored, indexed, or compared.

    A primer is a pure DNA oligo — every space / tab / newline in it is an
    artifact, never data. Sources that inject one:
      * BioPython's GenBank wrap-rejoin of a long ``/primer_seq`` (a value
        wrapped at the ~58-col boundary reloads as ``…GGAA TGAT``) — the
        original "phantom gap in the primer" report;
      * a hand-edited ``.gb`` / ``.dna`` qualifier or a pasted oligo with a
        stray space / trailing newline;
      * a multi-line qualifier from any external importer.
    Drawn verbatim, any of these becomes a phantom one-column gap that shifts
    the primer's bases off the template it anneals to — catastrophic for a tool
    whose primers must DISPLAY exactly where they bind
    ([[project_primer_design_catastrophic]]). Routing every read-for-render
    (`PlasmidMap._parse`), every load (`_gb_text_to_record` →
    `_repair_wrapped_primer_seqs`), every qualifier write, the primer editor,
    and the usage index (`_index_primer_usage_in_collections` /
    `_find_primer_plasmid_usages`) through here means no workflow or insertion
    point can leak a dirty primer onto the canvas. Non-string input → ``""``.
    """
    if not isinstance(raw, str):
        return ""
    return re.sub(r"\s+", "", raw).upper()


def _arrowless_encode_features(features):
    """Return a features list in which every strand-0 / strand-None feature
    carries ``SpliceCraft_strand=["none"]`` (so the GenBank writer's plain
    location round-trips as arrowless, not forward). Pure: the input features
    are never mutated — only the tagged ones are replaced with shallow copies
    bearing a fresh qualifiers dict. Features that already carry a
    ``SpliceCraft_strand`` value (the double-strand ``["double"]`` marker) and
    the ``source`` feature are left untouched. Returns the SAME list object
    when nothing needs encoding (common case → zero allocation)."""
    feats = features or []
    out: list = []
    changed = False
    for f in feats:
        loc = getattr(f, "location", None)
        strand = getattr(loc, "strand", None) if loc is not None else None
        quals = getattr(f, "qualifiers", None) or {}
        if (strand in (0, None) and getattr(f, "type", "") != "source"
                and _SC_STRAND_QUAL not in quals):
            nf = _shallow_copy(f)
            nf.qualifiers = {**quals, _SC_STRAND_QUAL: ["none"]}
            out.append(nf)
            changed = True
        else:
            out.append(f)
    return out if changed else features


def _arrowless_decode_features(record):
    """Inverse of `_arrowless_encode_features`, applied to a freshly-parsed
    record we own: any feature tagged ``SpliceCraft_strand=["none"]`` has its
    location strand restored to 0 (arrowless) and the marker stripped, so the
    live record is both faithful (location.strand == 0) and clean (no leftover
    marker to confuse a later edit or re-export). The strand setter propagates
    to every part of a compound/wrap location. Mutates in place; returns the
    record."""
    for f in getattr(record, "features", None) or []:
        quals = getattr(f, "qualifiers", None)
        if not isinstance(quals, dict):
            continue
        q = quals.get(_SC_STRAND_QUAL)
        if isinstance(q, list) and q and str(q[0]).strip().lower() == "none":
            loc = getattr(f, "location", None)
            if loc is not None:
                try:
                    loc.strand = 0
                except (AttributeError, ValueError, TypeError):
                    pass
            quals.pop(_SC_STRAND_QUAL, None)
    return record


def _repair_wrapped_primer_seqs(record):
    """Strip embedded whitespace from every ``/primer_seq`` qualifier on a
    freshly-parsed record.

    A primer sequence is pure DNA — it never contains spaces. But BioPython's
    GenBank writer WRAPS a long qualifier value across lines, and its parser
    rejoins the continuation line with a SPACE (the right thing for free-text
    ``/note``, the wrong thing for a sequence). So any primer longer than the
    GenBank value-wrap width (~58 chars) round-trips through ``.gb`` with a
    literal space jammed mid-sequence:

        /primer_seq="GCGCCGTCTCAAATGAATAAATGTATTCCAATGATAATTAATGGAA
                     TGAT"           →  "…GGAA TGAT"

    The seq panel faithfully drew that space as a one-column gap, shifting the
    primer's 3' bases off the template (user-reported "gap in the primer that
    isn't a real mismatch"). Stripping on load keeps the in-memory primer clean
    for the bound-bar renderer, the primer editor, the binding re-derivation,
    and any re-export. Idempotent: a re-save re-wraps the value, the next load
    strips it again. Mutates in place; returns the record.
    [primer display is catastrophic-class — project_primer_design_catastrophic]
    """
    for f in getattr(record, "features", None) or []:
        if getattr(f, "type", "") != "primer_bind":
            continue
        quals = getattr(f, "qualifiers", None)
        if not isinstance(quals, dict):
            continue
        ps = quals.get("primer_seq")
        if isinstance(ps, list) and ps:
            cleaned = [_normalize_primer_seq(v) for v in ps]
            if cleaned != [str(v) for v in ps]:
                quals["primer_seq"] = cleaned
    return record


def _split_multiline_qualifiers(features):
    """Return a features list in which no qualifier VALUE contains an embedded
    newline.

    GenBank cannot encode a newline INSIDE a single qualifier value: BioPython's
    writer emits the raw newline with the continuation text flush against
    column 0, and BioPython's OWN parser then rejects that on reload ("Problem
    with '<type>' feature"). Because every plasmid persists as its GenBank text,
    a feature ``/note`` typed with a single Enter (the everyday multi-line note)
    would serialise to an UNPARSEABLE record — the saved entry could never be
    loaded again.

    Split every newline-bearing value into one entry per non-blank line, which
    the GenBank "repeat the qualifier" convention round-trips cleanly (the note
    editor already stores multi-paragraph notes as exactly such a list). This is
    the universal guard: every serialisation path (feature editor, `.dna` / GFF3
    import, agent API) funnels through `_record_to_gb_text`.

    Pure / copy-on-write: a feature is only replaced (shallow copy + fresh
    qualifiers dict) when it actually carries a multi-line value, so the caller's
    record is never mutated and the common single-line case allocates nothing
    (returns the SAME list object)."""
    feats = features or []

    def _has_multiline(quals):
        for v in quals.values():
            if isinstance(v, str) and "\n" in v:
                return True
            if isinstance(v, list) and any(
                    isinstance(x, str) and "\n" in x for x in v):
                return True
        return False

    def _split_value(v):
        # Normalise CR / CRLF to LF, split on LF, drop blank lines (collapsing
        # paragraph gaps exactly like the editor's blank-line split). Fall back
        # to a single empty entry so an all-whitespace value still emits a
        # parseable qualifier instead of vanishing.
        s = v.replace("\r\n", "\n").replace("\r", "\n")
        return [ln for ln in s.split("\n") if ln.strip()] or [""]

    out: list = []
    changed = False
    for f in feats:
        quals = getattr(f, "qualifiers", None)
        if not isinstance(quals, dict) or not _has_multiline(quals):
            out.append(f)
            continue
        new_quals: dict = {}
        for k, v in quals.items():
            if isinstance(v, str) and "\n" in v:
                new_quals[k] = _split_value(v)
            elif isinstance(v, list):
                expanded: list = []
                for x in v:
                    if isinstance(x, str) and "\n" in x:
                        expanded.extend(_split_value(x))
                    else:
                        expanded.append(x)
                new_quals[k] = expanded
            else:
                new_quals[k] = v
        nf = _shallow_copy(f)
        nf.qualifiers = new_quals
        out.append(nf)
        changed = True
    return out if changed else features


# SC-E: the COMMENT marker that carries the human display name through a
# gb_text / file round-trip (the LOCUS line can't, and the library entry's
# name field doesn't survive an export). Written by `_record_to_gb_text`,
# read back by `_restore_display_name_from_comment`.
_DISPLAY_NAME_MARKER = "SpliceCraft-name:"


def _restore_display_name_from_comment(rec) -> None:
    """In-place: if the parsed record's COMMENT carries a SpliceCraft display-
    name marker, stamp it onto ``rec._tui_display_name`` — UNLESS the caller
    already set one (e.g. load-file from a hyphen-rich filename, which is
    equally authoritative and should win when present). The marker is the
    last line of our COMMENT stamp, so `(.+)` to the line end captures it;
    a name long enough to have wrapped degrades to its first physical line
    rather than being lost."""
    if getattr(rec, "_tui_display_name", None):
        return
    comment = (getattr(rec, "annotations", None) or {}).get("comment", "")
    if isinstance(comment, (list, tuple)):
        comment = "\n".join(str(x) for x in comment)
    m = re.search(re.escape(_DISPLAY_NAME_MARKER) + r"\s*(.+)", str(comment or ""))
    if not m:
        return
    name = m.group(1).splitlines()[0].strip()
    if name:
        try:
            rec._tui_display_name = name  # type: ignore[attr-defined]
        except Exception:
            pass


def _topology_from_gb_text(gb_text, default: str = "circular") -> str:
    """Read the molecule topology out of a GenBank LOCUS line.

    LOCUS carries topology as its own FIELD, between the molecule type and
    the division::

        LOCUS       pUC19    2686 bp    DNA     circular  SYN  27-OCT-2024

    so it is matched as a whole TOKEN, on the LOCUS line only. Three call
    sites used to ask `"linear" in gb_text[:200].lower()` instead, which also
    reads DEFINITION / ACCESSION / the plasmid's own NAME and answers with
    whatever word happens to appear there: an entry named `pLinear2`, or a
    circular vector whose DEFINITION says "linearized with EcoRI", reported
    the wrong topology — and topology is not cosmetic, it decides whether the
    origin-wrap scan runs at all (sacred invariant #6) and what the Kind
    column and the agent API tell a caller.

    ``default`` is returned when the text is empty, does not start with a
    LOCUS line, or carries neither token — the historical behaviour at each
    call site (mostly ``"circular"``: nearly every stored entry is a plasmid).
    """
    if not isinstance(gb_text, str) or not gb_text:
        return default
    line = gb_text[:200].lstrip().split("\n", 1)[0]
    if line[:5].upper() != "LOCUS":
        return default
    tokens = {t.lower() for t in line.split()}
    if "linear" in tokens:
        return "linear"
    if "circular" in tokens:
        return "circular"
    return default


def _record_to_gb_text(record) -> str:
    """Serialize a SeqRecord to GenBank format text.

    Biopython's genbank writer requires `molecule_type` in annotations
    — if the record came from elsewhere and doesn't have it, default to
    "DNA" rather than crashing. The fill-in happens on a shallow
    SeqRecord copy so the caller's record is never mutated (avoids
    subtle races with concurrent readers and surprise side effects for
    callers that compare records by annotation contents).

    Arrowless (strand 0) features are tagged with `SpliceCraft_strand=["none"]`
    here so they survive the round-trip (`_arrowless_encode_features`); the
    parse side restores them. No-op when nothing is arrowless.
    """
    from Bio import SeqIO
    anns = dict(getattr(record, "annotations", None) or {})
    anns.setdefault("molecule_type", "DNA")
    # Provenance: stamp which SpliceCraft version + date first wrote this file
    # into the GenBank COMMENT, so "which version made this?" is answerable
    # from the file alone (a legacy fragment's origin was otherwise
    # unknowable). Preserve an EXISTING stamp — the date reflects CREATION,
    # not the most recent re-save — and never duplicate it on round-trips.
    _prov = anns.get("comment", "")
    if isinstance(_prov, (list, tuple)):
        _prov = "\n".join(str(x) for x in _prov)
    _prov = str(_prov or "")
    if "Created by SpliceCraft v" not in _prov:
        from datetime import date
        _stamp = (f"Created by SpliceCraft v{_state._sc_version} "
                  f"on {date.today().isoformat()}")
        anns["comment"] = f"{_prov}\n{_stamp}".strip() if _prov else _stamp
    # SC-E: persist the human display name into the COMMENT so it survives a
    # gb_text / file round-trip. The LOCUS can't hold spaces / hyphens / >16
    # chars, and the library entry's `name` field is GONE the moment you
    # export to a `.gb` — so without this an exported construct re-imports
    # under the mangled LOCUS and needs a manual `rename-plasmid` every time.
    # The parse side (`_restore_display_name_from_comment`) reads it back onto
    # `_tui_display_name`. Idempotent: never re-stamps on a round-trip, and
    # only stamps when the display name actually differs from the LOCUS-safe
    # name (no point otherwise). [INV-98]
    _disp = getattr(record, "_tui_display_name", None)
    if isinstance(_disp, str):
        _disp_clean = _disp.replace("\n", " ").replace("\r", " ").strip()
        _cur = anns.get("comment", "")
        if isinstance(_cur, (list, tuple)):
            _cur = "\n".join(str(x) for x in _cur)
        _cur = str(_cur or "")
        _locus_form = re.sub(r"\s+", "_",
                             str(getattr(record, "name", "") or "").strip())
        if (_disp_clean and _disp_clean != _locus_form
                and _DISPLAY_NAME_MARKER not in _cur):
            _nm = f"{_DISPLAY_NAME_MARKER} {_disp_clean}"
            anns["comment"] = f"{_cur}\n{_nm}".strip() if _cur else _nm
    rec = _shallow_copy(record)
    rec.annotations = anns
    # The GenBank LOCUS line forbids whitespace and caps length — Biopython
    # raises "Invalid whitespace in '…' for LOCUS line" otherwise, which then
    # breaks EVERY save + autosave (user-reported after a scrub Add-to-Map: a
    # spaced display name had leaked into `record.name`). The human/display
    # name lives in `_tui_display_name` + the library entry, NOT the LOCUS, so
    # sanitise the LOCUS here on the COPY (never mutating the caller's record):
    # collapse whitespace to `_`, fall back to the id, cap to the INSDC max.
    _loc = re.sub(r"\s+", "_", str(getattr(rec, "name", "") or "").strip())
    if not _loc:
        _loc = re.sub(r"\s+", "_", str(getattr(rec, "id", "") or "").strip())
    rec.name = _loc[:_GB_LOCUS_NAME_MAX] or "PLASMID"
    rec.features = _split_multiline_qualifiers(
        _arrowless_encode_features(getattr(record, "features", None)))
    buf = StringIO()
    SeqIO.write(rec, buf, "genbank")
    return buf.getvalue()


def _clone_cached_record(rec):
    """Isolation copy for the `_GB_PARSE_CACHE` hit/store paths — the cheap,
    provably-safe replacement for a full `deepcopy` (~7-8× faster; 57 ms → 8 ms
    on a 5 Mb / 5000-feature chromosome, the per-hit cost every load-entry /
    diff / part-classify used to pay).

    The contract is pitfall #17: a caller may mutate the returned record freely
    without poisoning the cache. Callers only ever (a) reassign attributes
    (rebinds on THIS copy), (b) add/remove features (we hand back a fresh list),
    or (c) edit qualifiers (fresh dict + fresh value lists). They NEVER mutate a
    location, the `Seq`, `sub_features`, `letter_annotations`, or an annotation
    VALUE in place — verified by a codebase-wide grep audit (2026-06-23) and by
    sacred invariant #9 (feature edits rebuild whole records via
    `_rebuild_record_with_edit`). So those immutable-in-practice objects are
    shared by reference; only the mutable CONTAINERS are duplicated. If a future
    caller ever mutates a shared object in place, `test_record_cache.py`'s
    isolation test is the tripwire."""
    out = _shallow_copy(rec)                       # shares .seq (immutable Seq)
    out.annotations = dict(getattr(rec, "annotations", {}) or {})
    out.dbxrefs = list(getattr(rec, "dbxrefs", []) or [])
    la = getattr(rec, "letter_annotations", None)
    if la:
        # Restricted dict — assign via the public setter semantics (plain copy
        # of the per-letter arrays; they're never mutated in place, only read).
        out.letter_annotations = dict(la)
    new_feats = []
    for f in getattr(rec, "features", None) or []:
        f2 = _shallow_copy(f)                      # shares .location (immutable)
        f2.qualifiers = {k: list(v) for k, v in (f.qualifiers or {}).items()}
        new_feats.append(f2)
    out.features = new_feats
    return out


def _gb_cache_key(text: str) -> str:
    """128-bit content key for `_GB_PARSE_CACHE`. Plain `hash(text)` is only
    64-bit, so a collision between two DIFFERENT GenBank texts would return the
    WRONG cached SeqRecord from this mission-critical serialization path.
    blake2b-128 makes that collision astronomically impossible without the
    memory cost of storing every source text or a per-hit full-string compare."""
    import hashlib
    return hashlib.blake2b(
        text.encode("utf-8", "surrogatepass"), digest_size=16,
    ).hexdigest()


def _gb_text_to_record(text: str, *, cache: bool = True):
    """Parse GenBank format text back to a SeqRecord.

    Defence-in-depth: cap input length at 64 MB before handing to
    BioPython's parser. Library entries are gated through
    `_safe_load_json`'s 1 GB cap and zip extracts through the 50 MB
    member cap, but the GenBank parser itself is unbounded internally;
    a single record line that happens to slip past upstream caps would
    otherwise allocate intermediate parser objects without ceiling.

    Results are LRU-cached (`_GB_PARSE_CACHE`, capped at
    `_GB_PARSE_CACHE_MAX`) keyed on a 128-bit content hash. Returned records are
    isolation copies (`_clone_cached_record` — shares the immutable
    `Seq`/locations, duplicates the mutable feature/qualifier containers)
    so callers can mutate freely (pitfall #17 contract). Empty-string /
    oversize inputs bypass the cache.

    Sweep #26 (2026-05-23) — pass ``cache=False`` for one-shot batch
    parses (Plasmidsaurus zip ingest, bulk-import folder walk) where
    the cache would absorb the entire batch's parsed records before
    eviction. A 50-sample run × 5 MB assemblies = ~250 MB cache
    pressure on a single click; the records are consumed immediately
    so caching them has no payoff.
    """
    if not text:
        # Explicit, clean error instead of handing "" to SeqIO.read (whose
        # "No records found in handle" ValueError is opaque). Callers already
        # pre-check or catch ValueError; this just makes the message legible.
        raise ValueError("empty GenBank text")
    if len(text) > _GB_TEXT_MAX_BYTES:
        raise ValueError(
            f"GenBank text too large to parse "
            f"({len(text):,} bytes > "
            f"{_GB_TEXT_MAX_BYTES:,} cap)"
        )
    if cache:
        key = _gb_cache_key(text)
        with _GB_PARSE_CACHE_LOCK:
            hit = _GB_PARSE_CACHE.get(key)
            if hit is not None:
                _GB_PARSE_CACHE.move_to_end(key)
                return _clone_cached_record(hit)
    from Bio import SeqIO
    rec = SeqIO.read(StringIO(text), "genbank")
    # SC-E: restore the human display name stamped into the COMMENT so it
    # survives the round-trip (the LOCUS is the underscored/truncated form).
    # Done before caching so every caller — load-entry, add-current, diff —
    # sees the real name without a manual rename.
    _restore_display_name_from_comment(rec)
    # Restore arrowless (strand 0) features tagged on the write side before
    # the result is cached, so every caller sees a faithful + clean record.
    _arrowless_decode_features(rec)
    # Undo BioPython's wrap-rejoin space corruption of long `/primer_seq`
    # qualifiers BEFORE caching, so the renderer never paints a phantom gap.
    _repair_wrapped_primer_seqs(rec)
    if cache:
        with _GB_PARSE_CACHE_LOCK:
            # Store an isolation copy so the caller (who owns `rec`) can keep
            # mutating it without poisoning the cache; return `rec` itself.
            # Recompute the key here (rather than reuse the lookup block's
            # `key`) so it's unconditionally bound in this separate `if cache:`
            # block — miss-path only, so the extra hash is dwarfed by the parse.
            _GB_PARSE_CACHE[_gb_cache_key(text)] = _clone_cached_record(rec)
            while len(_GB_PARSE_CACHE) > _GB_PARSE_CACHE_MAX:
                _GB_PARSE_CACHE.popitem(last=False)
    return rec
