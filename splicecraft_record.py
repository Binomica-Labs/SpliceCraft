"""splicecraft_record — GenBank <-> SeqRecord serialization (Phase D, layer L1).

The hub's record serialization core, extracted so the modal/screen siblings that
parse/serialise records can import it instead of calling it bare on the hub. This
is MISSION-CRITICAL (INV-98: the GenBank LOCUS must never carry the human/display
name; the .dna + library round-trips depend on byte-faithful parse/serialise), so
it moved as one self-contained unit with its own LRU parse cache. Biopython is
imported lazily inside the two entry points (the writer subclass and
`SeqIO.read`) exactly as in the hub. The SpliceCraft version stamped into the GenBank COMMENT is read from
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
from splicecraft_logging import _log


_SC_STRAND_QUAL = "SpliceCraft_strand"

_GB_TEXT_MAX_BYTES = 64 * 1024 * 1024

_GB_PARSE_CACHE: "_OD" = _OD()

_GB_PARSE_CACHE_MAX = 16

_GB_PARSE_CACHE_LOCK = threading.RLock()

_GB_LOCUS_NAME_MAX = 28  # NCBI relaxed LOCUS name length (spec is 16)

# The classic 80-column LOCUS line packs the name and the sequence length into
# one 28-character field (name left-justified from column 13, length
# right-justified ending at column 40), which is why NCBI's relaxation of the
# 16-char name limit only stretches so far: as soon as
# ``len(name) + len(str(bp)) > 27`` Biopython gives up on the fixed columns and
# emits ``<name> <bp> bp ...`` instead. That line is LONGER than 80 characters
# and every field after the length — units, molecule type, topology, division,
# date — lands in the wrong column, so any reader that slices the LOCUS line
# positionally (and plenty still do) reads garbage for the topology.
# `_locus_name_cap` is what keeps us on the standard layout.
_GB_LOCUS_NAME_FIELD = 28


def _locus_name_cap(seq_len: int) -> int:
    """Longest LOCUS name that still yields a column-conformant LOCUS line
    for a record of ``seq_len`` bases.

    Never returns less than 8 — a pathological length (>19 digits) would
    otherwise shrink the name to nothing, and at that point the file is
    already outside anything the format contemplates.
    """
    try:
        digits = len(str(max(0, int(seq_len))))
    except (TypeError, ValueError):
        digits = 11
    return max(8, min(_GB_LOCUS_NAME_MAX, _GB_LOCUS_NAME_FIELD - 1 - digits))


# INSDC restricts the LOCUS name to letters, digits and underscore (NCBI also
# tolerates `.` and `-` in practice). Anything else — the `+` in a construct
# called "DEMO 33 MOD CDS+REPORTER", a `%` or `#` from a pasted name — is
# folded to `_`. Only whitespace used to be folded, so a name with a `+` in it
# put a `+` straight into the LOCUS.
_LOCUS_BAD_CHARS_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def _locus_name_for(record, *, cap: "int | None" = None) -> str:
    """THE LOCUS name for `record`: sanitised to the INSDC character set and
    truncated to the length that keeps the 80-column layout intact.

    One function so the value the writer stamps and the value the display-name
    COMMENT marker is compared against can never drift apart — when they did,
    a name the LOCUS could not hold went unstamped and its truncation became
    permanent. [INV-98]
    """
    if cap is None:
        try:
            cap = _locus_name_cap(len(record.seq))
        except (AttributeError, TypeError):
            cap = _GB_LOCUS_NAME_MAX
    raw = str(getattr(record, "name", "") or "").strip()
    if not raw:
        raw = str(getattr(record, "id", "") or "").strip()
    return _LOCUS_BAD_CHARS_RE.sub("_", raw)[:cap] or "PLASMID"


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
        # A CompoundLocation whose PARTS disagree on strand also reports
        # `.strand is None` — Biopython's way of saying "mixed", not
        # "arrowless". Tagging it drove the decoder to force strand 0 onto
        # every part, so a legal `join(11..50,complement(61..90))` came back
        # with its second part no longer reverse-complemented and its
        # extracted DNA silently changed. Only a location that is arrowless
        # in EVERY part is arrowless.
        if strand is None and loc is not None:
            parts = getattr(loc, "parts", None)
            if parts and not all(
                    getattr(p, "strand", None) in (0, None) for p in parts):
                out.append(f)
                continue
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
    record we own: a feature tagged ``SpliceCraft_strand=["none"]`` has its
    location strand restored to 0 (arrowless) and the marker stripped, so the
    live record is both faithful (location.strand == 0) and clean (no leftover
    marker to confuse a later edit or re-export). The strand setter propagates
    to every part of a compound/wrap location. Mutates in place; returns the
    record.

    ``SpliceCraft_strand=["double"]`` (the SpliceCraft-only "both arrows"
    convention) gets its strand restored to 0 as well but **keeps** the
    marker — that is what `PlasmidMap._parse` reads to promote the render
    back to strand 2. Without this the strand half of the round-trip was
    lossy in a way the marker hid: the writer emitted a plain location for
    the record's BioPython strand 0/None, the reader handed back `1`, and a
    double-stranded feature came home pointing FORWARD while still carrying
    a "double" marker. Because the marker drives the render, the map looked
    right and only `_export_genbank_to_path`'s signature guard noticed —
    refusing to export any record holding a double-stranded feature (agent
    field report 2026-09-12, found while tracing the strand-0 export trap)."""
    for f in getattr(record, "features", None) or []:
        quals = getattr(f, "qualifiers", None)
        if not isinstance(quals, dict):
            continue
        q = quals.get(_SC_STRAND_QUAL)
        if not (isinstance(q, list) and q):
            continue
        marker = str(q[0]).strip().lower()
        if marker not in ("none", "double"):
            continue
        loc = getattr(f, "location", None)
        if loc is not None:
            try:
                loc.strand = 0
            except (AttributeError, ValueError, TypeError):
                pass
        if marker == "none":
            # "double" is a live render instruction, not an encode artefact.
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
        # ANY feature type, not just `primer_bind`. SpliceCraft only writes
        # `/primer_seq` onto a primer_bind, but a `.dna` / GenBank file from
        # elsewhere can carry it anywhere, and a primer sequence is a primer
        # sequence wherever it sits — whitespace in it is an artefact either
        # way. Keying on the type left those unrepaired, which (once the
        # writer started wrapping long `/primer_seq` values at 80 columns)
        # made the export's round-trip guard refuse the file.
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
    LAST line of our COMMENT stamp, so everything after it — across however
    many physical lines Biopython's 68-column COMMENT wrapping produced — is
    the name, rejoined on single spaces. Keeping only the first physical line
    silently truncated every display name past ~50 characters, and the next
    save then re-stamped the truncation, making it permanent."""
    if getattr(rec, "_tui_display_name", None):
        return
    comment = (getattr(rec, "annotations", None) or {}).get("comment", "")
    if isinstance(comment, (list, tuple)):
        comment = "\n".join(str(x) for x in comment)
    # `[ \t]*` not `\s*`: `\s` matches a NEWLINE, so a marker with an empty
    # value swallowed the following comment line and restored it as the name.
    m = re.search(re.escape(_DISPLAY_NAME_MARKER) + r"[ \t]*(.+)",
                  str(comment or ""), re.S)
    if not m:
        return
    # The marker is the LAST line of our COMMENT stamp, so everything after it
    # belongs to the name. Biopython wraps COMMENT at 68 columns, so a name
    # with spaces arrives split across physical lines; keeping only the first
    # of them silently truncated every display name past ~50 characters (and
    # then re-stamped the truncation, making it permanent).
    name = " ".join(m.group(1).split())
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


def _backfill_topology(rec, gb_text: str) -> None:
    """In-place: set ``rec.annotations["topology"]`` from the LOCUS line when
    Biopython left it unset.

    Biopython reads the LOCUS topology POSITIONALLY (columns 56-63). Plenty of
    real files — ApE's, older SnapGene exports, anything hand-edited — put the
    same `circular` token in a slightly different column, and Biopython then
    silently reports no topology at all. That is not cosmetic: the export path
    defaults a topology-less record to **linear**, so loading one of those
    files and exporting it turned a circular plasmid into a linear fragment for
    every downstream tool, while SpliceCraft's own map (which reads the LOCUS
    line by token, via `_topology_from_gb_text`) kept drawing it circular.

    Only ever ADDS a value — an explicit annotation always wins, and a LOCUS
    line carrying neither token is left alone rather than guessed at.
    """
    anns = getattr(rec, "annotations", None)
    if not isinstance(anns, dict) or anns.get("topology"):
        return
    topo = _topology_from_gb_text(gb_text, default="")
    if topo:
        anns["topology"] = topo


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
    # `_tui_display_name`. The marker simply MIRRORS `_tui_display_name`: it is
    # rewritten when that changes (so a rename reaches the file), dropped when
    # the LOCUS can carry the name unchanged, and left alone when the record
    # has no display name at all — which keeps a round-trip idempotent. [INV-98]
    try:
        _locus_cap = _locus_name_cap(len(record.seq))
    except (AttributeError, TypeError):
        _locus_cap = _GB_LOCUS_NAME_MAX
    _disp = getattr(record, "_tui_display_name", None)
    if isinstance(_disp, str):
        _disp_clean = _disp.replace("\n", " ").replace("\r", " ").strip()
        _cur = anns.get("comment", "")
        if isinstance(_cur, (list, tuple)):
            _cur = "\n".join(str(x) for x in _cur)
        _cur = str(_cur or "")
        # Compare against the FINAL locus form — sanitised AND truncated — not
        # the raw `record.name`: a display name the LOCUS cannot hold verbatim
        # but which happens to equal `record.name` used to skip the stamp, so
        # the mangling the writer then applied was permanent and
        # unrecoverable. [INV-98]
        _locus_form = _locus_name_for(record, cap=_locus_cap)
        # REWRITE the marker rather than skip when one is already there. The
        # old "never re-stamp" rule kept the round-trip idempotent but made a
        # RENAME invisible to the file: the entry had been saved once with
        # `SpliceCraft-name: My Plasmid`, the user renamed it, and the marker
        # still said the old name — so exporting (or re-importing) the entry
        # brought the old name back. The marker is now simply whatever
        # `_tui_display_name` says, and is dropped when the LOCUS can carry
        # the name unchanged. Still idempotent: re-serialising a record whose
        # marker already matches produces the same text.
        _before, _sep, _ = _cur.partition(_DISPLAY_NAME_MARKER)
        _head = _before.rstrip() if _sep else _cur
        if _disp_clean and _disp_clean != _locus_form:
            _nm = f"{_DISPLAY_NAME_MARKER} {_disp_clean}"
            anns["comment"] = f"{_head}\n{_nm}".strip() if _head else _nm
        elif _sep:
            anns["comment"] = _head
    rec = _shallow_copy(record)
    rec.annotations = anns
    # The GenBank LOCUS line forbids whitespace and caps length — Biopython
    # raises "Invalid whitespace in '…' for LOCUS line" otherwise, which then
    # breaks EVERY save + autosave (user-reported after a scrub Add-to-Map: a
    # spaced display name had leaked into `record.name`). The human/display
    # name lives in `_tui_display_name` + the library entry, NOT the LOCUS, so
    # sanitise the LOCUS here on the COPY (never mutating the caller's record)
    # via `_locus_name_for`, which folds whitespace AND every other character
    # outside the INSDC set, falls back to the id, and applies the DYNAMIC cap
    # (`_locus_name_cap`): the name and the sequence length share one
    # 28-column field, so how much name fits depends on how many digits the
    # length needs. Overflowing it makes Biopython abandon the fixed-column
    # layout entirely — see `_locus_name_cap`.
    rec.name = _locus_name_for(rec, cap=_locus_cap)
    rec.features = _split_multiline_qualifiers(
        _arrowless_encode_features(getattr(record, "features", None)))
    buf = StringIO()
    _unbreakable_genbank_writer(buf).write_file([rec])
    return buf.getvalue()


# Qualifiers a long value may be HARD-BROKEN in, because whitespace in them is
# never data:
#
#   * `/translation` — the case the standard itself names. NCBI wraps every
#     protein translation mid-token, and Biopython's parser lists it in
#     `FeatureValueCleaner.keys_to_process` precisely so the continuation
#     rejoins with no separator. Refusing to break it is what produced
#     7,132-character lines for SARS-CoV-2 ORF1ab.
#   * `/primer_seq` — a vendor extension (not INSDC), and its value is a bare
#     nucleotide string. `_normalize_primer_seq` strips whitespace from it on
#     every read, so the break is lossless here; anything else reading it as a
#     sequence will strip too. A 60-mer Gibson or mutagenesis primer overflows
#     the line otherwise, which is far likelier to break a reader than a space
#     in an oligo is to mislead one.
#
# Every OTHER qualifier rejoins with a space that was never in the data, which
# is why `_unbreakable_genbank_writer` refuses to break those at all — see its
# docstring for that trade.
_HARD_BREAKABLE_QUALS: frozenset = frozenset({"translation", "primer_seq"})

# Back-compat alias for the name this set shipped under.
_NO_SPACE_REJOIN_QUALS = _HARD_BREAKABLE_QUALS


def _unbreakable_genbank_writer(handle):
    """A `GenBankWriter` that never HARD-BREAKS a whitespace-free qualifier
    value, built lazily so Biopython stays a deferred import.

    Biopython wraps a long qualifier at column 80. When the value contains a
    space it breaks there, and its own parser rejoins the continuation on that
    space — lossless. When the value contains NO space in reach it breaks
    mid-token anyway, and the parser rejoins with a space that was never
    there. The result is a silent, permanent mutation of user annotation on
    every load: a 51-character label came back with a trailing space, and
    `Construct_optimized_for_expression_v2_final_long_name` came back with a
    space inside it. Because the library persists every plasmid as GenBank
    text, this fired on ordinary saves, and `_export_genbank_to_path`'s
    round-trip guard then REFUSED to export the record at all.

    This writer instead breaks after the offending token (or, failing that,
    emits the rest of the line whole). The lines it produces are longer than
    80 columns, which INSDC discourages and every parser in practice accepts —
    a strictly better trade than corrupting the value.

    EXCEPT for `/translation` and its kin (`_NO_SPACE_REJOIN_QUALS`), where the
    rejoin inserts nothing, so the stock hard break is both lossless AND what
    NCBI itself emits. Routing those through this writer's "never break" path
    was a real regression: every record carrying a CDS translation — i.e. every
    record NCBI serves — exported with the whole protein on one line, up to
    7,132 characters for SARS-CoV-2 ORF1ab. `test_genbank_io.py` pins the
    80-column result on a real translation."""
    from Bio.SeqIO.InsdcIO import GenBankWriter

    class _Writer(GenBankWriter):
        def _write_feature_qualifier(self, key, value=None, quote=None):
            if value is None or not isinstance(value, str):
                return super()._write_feature_qualifier(key, value, quote)
            if key in _HARD_BREAKABLE_QUALS:
                # Lossless to hard-break — let the stock 80-column wrapper run.
                return super()._write_feature_qualifier(key, value, quote)
            if (quote is None and key in self.FTQUAL_NO_QUOTE
                    and any(c.isspace() for c in value)):
                # `/codon_start`, `/transl_table` and friends are written
                # UNQUOTED because their values are bare tokens. When one of
                # them somehow holds whitespace (a hand-edited file, a `.dna`
                # import, an agent write), the unquoted form wraps into a
                # continuation line the reader rejoins with a newline — so the
                # value came back different and the export's round-trip guard
                # refused the whole file. Quoting it costs nothing for a
                # well-formed value (there is no whitespace to trigger this)
                # and makes a malformed one survive instead of blocking the
                # export.
                quote = True
            if " " in value:
                # Breakable the normal way, but a single over-long token
                # inside it can still trip the hard break — fall through to
                # the safe path only when one actually would.
                longest = max((len(t) for t in value.split(" ")), default=0)
                if longest + self.QUALIFIER_INDENT + len(key) + 4 <= self.MAX_WIDTH:
                    return super()._write_feature_qualifier(key, value, quote)
            esc = value.replace('"', '""')
            if quote is None:
                quote = not (isinstance(value, int) or key in self.FTQUAL_NO_QUOTE)
            line = (f'{self.QUALIFIER_INDENT_STR}/{key}="{esc}"' if quote
                    else f"{self.QUALIFIER_INDENT_STR}/{key}={esc}")
            while line.lstrip():
                if len(line) <= self.MAX_WIDTH:
                    self.handle.write(line + "\n")
                    return
                idx = line.rfind(" ", self.QUALIFIER_INDENT + 1, self.MAX_WIDTH + 1)
                if idx <= self.QUALIFIER_INDENT:
                    # No break point in reach — take the NEXT one instead of
                    # splitting the token, and emit the remainder whole when
                    # there is none.
                    idx = line.find(" ", self.MAX_WIDTH)
                    if idx == -1:
                        # Deliberate, documented spec deviation. Log it so an
                        # over-wide line in an exported file is traceable to
                        # the qualifier that forced it rather than looking
                        # like corruption.
                        _log.warning(
                            "GenBank: /%s value has no break point in reach; "
                            "emitting a %d-column line (breaking it would add "
                            "a space the reader cannot tell from real data)",
                            key, len(line),
                        )
                        self.handle.write(line + "\n")
                        return
                self.handle.write(line[:idx] + "\n")
                line = self.QUALIFIER_INDENT_STR + line[idx:].lstrip()

    return _Writer(handle)


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
    _backfill_topology(rec, text)
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
