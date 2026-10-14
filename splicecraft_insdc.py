"""splicecraft_insdc — SpliceCraft's own INSDC flat-file reader (GenBank).

Layer L0: pure text -> record parsing. No app state, no data dir, no network.
It builds the same Biopython MODEL objects SpliceCraft works in everywhere
(`SeqRecord`, `SeqFeature`, `SimpleLocation` / `CompoundLocation` and the
position classes, `Reference`), imported lazily, but it does its own READING.

Why own the reader [INV-219]: Biopython's GenBank scanner finds everything by
column and by exact match, and real files from other tools break both. The
format scoreboard (`scripts/format_scoreboard.py`, [INV-218]) measured what
that costs: a LOCUS line with a spaced name or a field one column over is
refused outright; a tab or a shifted line in the feature table becomes a
garbage feature (`/organism="Cloning` as a feature KEY) or an AssertionError;
a sequence line without position numbers silently loses its first ten bases;
an origin-spanning feature on a record not marked circular gets location
None, and the next save then crashes.

The contract, enforced by `tests/test_insdc_reader.py`:

* **On every file Biopython reads correctly, this reads the same record** —
  name, id, description, dbxrefs, every annotation (references included),
  every feature (type, location object for location object, qualifiers in
  order) and the sequence. Measured on the scoreboard corpus: 313 NCBI
  records, Biopython's own GenBank fixtures, and a 1,303-entry plasmid
  library. Biopython's quirks are kept where a file depends on them (a
  wrapped value rejoins with ONE space; `KEYWORDS .` is `[""]`).
* **Where Biopython refuses, crashes or reads wrong, this reads what the
  file means** and says so with an `InsdcParserWarning` — it does not refuse
  a file it can make sense of, and it does not invent features. Every such
  repair is listed in `_REPAIRS` and counted on the record
  (`rec._insdc_repairs`) so the app can tell the user.
* **Never a location of None.** A location that cannot be read at all drops
  its feature, with a warning naming it — a None location crashed every
  save of the record downstream.

`parse_genbank(text)` returns every record in a text; `read_genbank(text)`
exactly one. Input is TEXT (the caller decodes: BOM, CP1252 and line endings
are `splicecraft_fileio._read_text_tolerant`'s job).
"""
from __future__ import annotations

import re
import warnings
from typing import Any

__all__ = ["InsdcParserWarning", "parse_genbank", "read_genbank",
           "genbank_record_count", "parse_location", "repair_note"]


class InsdcParserWarning(UserWarning):
    """Something in a flat file the reader had to repair or could not use."""


# Repairs counted on `rec._insdc_repairs` — what the app may tell the user.
_REPAIRS = (
    "origin_spanning_not_circular",   # start..end with start > end, record not
                                      # marked circular: kept as two joined parts
    "unreadable_location",            # feature dropped: location unreadable
    "negative_position",              # a position before base 1: kept as `<1`
    "locus_line",                     # LOCUS fields read by what they are
    "feature_table_layout",           # tabs / shifted lines read by content
    "unnumbered_sequence",            # sequence lines without position numbers
    "unclosed_quote",                 # a quoted value that never closed
)


def _warn(msg: str) -> None:
    warnings.warn(msg, InsdcParserWarning, stacklevel=3)


def _isnum(s: str) -> bool:
    """ASCII digits only: `str.isdigit` also says yes to `²`, which `int`
    then refuses."""
    return s.isascii() and s.isdigit()


# A real location nests three or four deep (`complement(join(complement(..`).
# Past this the reader gives up on the location instead of recursing into a
# RecursionError, which is not a ValueError and got past every caller.
_MAX_LOCATION_DEPTH = 32


def _location_too_deep(text: str) -> bool:
    if text.count("(") <= _MAX_LOCATION_DEPTH:
        return False
    depth = 0
    for ch in text:
        if ch == "(":
            depth += 1
            if depth > _MAX_LOCATION_DEPTH:
                return True
        elif ch == ")":
            depth -= 1
    return False


# ── the Biopython model, imported once on first use ─────────────────────────

_M: dict = {}


def _model() -> dict:
    if not _M:
        from Bio.Seq import Seq
        from Bio.SeqFeature import (
            AfterPosition, BeforePosition, CompoundLocation, ExactPosition,
            OneOfPosition, Reference, SeqFeature, SimpleLocation,
            UncertainPosition, UnknownPosition, WithinPosition,
        )
        from Bio.SeqRecord import SeqRecord
        _M.update(
            Seq=Seq, SeqRecord=SeqRecord, SeqFeature=SeqFeature,
            SimpleLocation=SimpleLocation, CompoundLocation=CompoundLocation,
            ExactPosition=ExactPosition, BeforePosition=BeforePosition,
            AfterPosition=AfterPosition, WithinPosition=WithinPosition,
            OneOfPosition=OneOfPosition, UnknownPosition=UnknownPosition,
            UncertainPosition=UncertainPosition, Reference=Reference)
    return _M


# ── locations ───────────────────────────────────────────────────────────────

class _LocationError(ValueError):
    pass


_SOLO_RE = re.compile(r"[<>]?\d+")
_WITHIN_RE = re.compile(r"\((\d+)\.(\d+)\)")
_ONEOF_RE = re.compile(r"one-of\((\d+(?:,\d+)+)\)")
_BETWEEN_RE = re.compile(r"(\d+)\^(\d+)")
_REF_RE = re.compile(r"[a-zA-Z][a-zA-Z0-9_.|]*[a-zA-Z0-9]?")
_POS_TOKEN_RE = re.compile(
    r"[<>]?-?\d+|\(\d+\.\d+\)|one-of\(\d+(?:,\d+)+\)|\?\d*")


def _position(text: str, offset: int):
    """One position of a range, as `Bio.SeqFeature.Position.fromstring`
    builds it (``offset`` -1 for a start, 0 for an end)."""
    M = _model()
    if text == "?":
        return M["UnknownPosition"]()
    c = text[0]
    if c == "?":
        return M["UncertainPosition"](int(text[1:]) + offset)
    if c == "<":
        return M["BeforePosition"](int(text[1:]) + offset)
    if c == ">":
        return M["AfterPosition"](int(text[1:]) + offset)
    if c == "(":
        m = _WITHIN_RE.fullmatch(text)
        if m is None:
            raise _LocationError(text)
        s, e = int(m.group(1)) + offset, int(m.group(2)) + offset
        return M["WithinPosition"](s if offset == -1 else e, left=s, right=e)
    if c == "o":
        m = _ONEOF_RE.fullmatch(text)
        if m is None:
            raise _LocationError(text)
        parts = [M["ExactPosition"](int(p) + offset)
                 for p in m.group(1).split(",")]
        default = (min(int(p) for p in parts) if offset == -1
                   else max(int(p) for p in parts))
        return M["OneOfPosition"](default, choices=parts)
    return M["ExactPosition"](int(text) + offset)


def _split_top(text: str) -> "list[str]":
    """``text`` split on the commas outside parentheses."""
    if "(" not in text:
        return text.split(",")
    parts, depth, start = [], 0, 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "," and depth == 0:
            parts.append(text[start:i])
            start = i + 1
    parts.append(text[start:])
    return parts


class _LocCtx:
    """What one location parse needs to know about its record, and what it
    had to repair."""
    __slots__ = ("length", "circular", "stranded", "repairs", "where")

    def __init__(self, length, circular, stranded, where=""):
        self.length = length
        self.circular = circular
        self.stranded = stranded
        self.repairs: "list[str]" = []
        self.where = where


def _simple(text: str, ctx: _LocCtx, strand) -> list:
    """One range (or point) as a list of SimpleLocation parts — two parts for
    an origin-spanning range — as `SimpleLocation.fromstring` reads it.
    ``strand`` is the strand to give it; a `complement(...)` part is -1."""
    M = _model()
    SimpleLocation = M["SimpleLocation"]
    if text.startswith("complement(") and text.endswith(")"):
        # A complement INSIDE a join: `join(complement(1..5),9..12)`. Always
        # the minus strand, even on a protein record (Biopython's reading).
        return _simple(text[11:-1], ctx, -1)
    # Fast path, Biopython's own: a plain `s..e`.
    dots = text.find("..")
    if dots > 0 and text[0].isdigit():
        a, b = text[:dots], text[dots + 2:]
        if _isnum(a) and _isnum(b):
            s, e = int(a) - 1, int(b)
            if 0 <= s < e:
                return [SimpleLocation(s, e, strand)]
    ref = None
    colon = text.find(":")
    if colon > 0:
        ref, text = text[:colon], text[colon + 1:]
        if not _REF_RE.fullmatch(ref) or ":" in text:
            raise _LocationError(text)
        if text.startswith("complement(") and text.endswith(")"):
            parts = _simple(text[11:-1], ctx, -1)
            for p in parts:
                p.ref = ref
            return parts
    if text.startswith("bond(") and text.endswith(")"):
        _warn("Dropping bond qualifier in feature location")
        text = text[5:-1]
    m = _BETWEEN_RE.fullmatch(text)
    if m is not None:
        s, e = int(m.group(1)), int(m.group(2))
        if s + 1 == e or (s == ctx.length and e == 1):
            pos = M["ExactPosition"](s)
            return [SimpleLocation(pos, pos, strand, ref=ref)]
        raise _LocationError(text)
    if ".." in text:
        a, _sep, b = text.partition("..")
        if not (_POS_TOKEN_RE.fullmatch(a) and _POS_TOKEN_RE.fullmatch(b)):
            raise _LocationError(text)
        neg = False
        if a.lstrip("<>").startswith("-"):
            # A position before base 1 — the feature runs off the 5' end of
            # the sequence, which INSDC writes `<1`. Biopython drops the whole
            # location ("negative starting position").
            a = "<1"
            neg = True
        s_pos = _position(a, -1)
        e_pos = _position(b, 0)
        if neg:
            ctx.repairs.append("negative_position")
            _warn(f"{ctx.where}: a position before base 1 in {text!r} — "
                  f"kept as a feature running off the start (<1)")
        if (isinstance(s_pos, int) and isinstance(e_pos, int)
                and e_pos <= s_pos and ctx.length and s_pos < ctx.length):
            # Spans the origin. On a circle that is the INSDC meaning, and
            # Biopython repairs it the same way. On a record NOT marked
            # circular Biopython gives up (location None — and the next save
            # crashed); the two parts are still exactly the bases the file
            # names, so keep them, and say so.
            if not ctx.circular:
                ctx.repairs.append("origin_spanning_not_circular")
                _warn(f"{ctx.where}: {text!r} runs across the origin of a "
                      f"record not marked circular — kept as two joined parts")
            else:
                _warn(f"Attempting to fix invalid location {text!r} as it "
                      f"looks like incorrect origin wrapping. Please fix input "
                      f"file, this could have unintended behavior.")
            f1 = SimpleLocation(s_pos, ctx.length, strand, ref=ref)
            f2 = SimpleLocation(0, e_pos, strand, ref=ref)
            return [f2, f1] if strand == -1 else [f1, f2]
        if isinstance(s_pos, int) and int(s_pos) < 0:
            raise _LocationError(text)
        try:
            return [SimpleLocation(s_pos, e_pos, strand, ref=ref)]
        except ValueError:
            raise _LocationError(text) from None
    if _SOLO_RE.fullmatch(text) or (text[:1] == "?"
                                     and _POS_TOKEN_RE.fullmatch(text)):
        s_pos = _position(text, -1)
        e_pos = _position(text, 0)
        if isinstance(s_pos, int) and int(s_pos) < 0:
            raise _LocationError(text)
        return [SimpleLocation(s_pos, e_pos, strand, ref=ref)]
    raise _LocationError(text)


def _compound(text: str, ctx: _LocCtx, strand) -> "tuple[list, str]":
    """``(parts, operator)`` of a `join(...)` / `order(...)` / `bond(...)`,
    or ``(parts, "")`` for a plain range."""
    for op in ("join", "order", "bond"):
        if text.startswith(op + "(") and text.endswith(")"):
            pieces = _split_top(text[len(op) + 1:-1])
            if any(not p for p in pieces):
                _warn("Dropping trailing comma in malformed feature location")
                pieces = [p for p in pieces if p]
            parts: list = []
            for p in pieces:
                if p.startswith(("join(", "order(", "bond(",
                                 "complement(join(", "complement(order(")):
                    # Nested operators are illegal INSDC (Biopython: location
                    # None); read what the nesting means.
                    _warn(f"{ctx.where}: nested operator in {text!r} — "
                          f"read as one {op}")
                    parts.extend(_flatten(p, ctx, strand))
                    continue
                parts.extend(_simple(p, ctx, strand))
            return parts, op
    return _simple(text, ctx, strand), ""


def _flatten(text: str, ctx: _LocCtx, strand) -> list:
    outer = text.startswith("complement(") and text.endswith(")")
    parts, _op = _compound(text[11:-1] if outer else text, ctx, strand)
    if outer:
        for p in parts:
            p.strand = 1 if p.strand == -1 else -1
        parts.reverse()
    return parts


def _location(text: str, ctx: _LocCtx):
    """A whitespace-free INSDC location string -> Biopython location, the
    way `Location.fromstring` builds it. Raises `_LocationError`."""
    outer = text.startswith("complement(")
    if outer:
        if not text.endswith(")"):
            raise _LocationError(text)
        text = text[11:-1]
    parts, op = _compound(text, ctx, 1 if ctx.stranded else None)
    if not parts:
        raise _LocationError(text)
    if outer:
        # `complement(join(a,b))` lists the parts forward; the feature reads
        # them back to front, each on the minus strand. A part complemented
        # twice reads forward (Biopython refuses it: "double complement").
        if any(p.strand == -1 for p in parts):
            _warn(f"{ctx.where}: a part complemented twice — read forward")
        for p in parts:
            p.strand = 1 if p.strand == -1 else -1
        parts.reverse()
    if len(parts) == 1:
        return parts[0]
    return _model()["CompoundLocation"](parts, operator=op or "join")


def parse_location(text: str, length: "int | None" = None,
                   circular: bool = False, stranded: bool = True):
    """Public face of the location reader: an INSDC location string to a
    Biopython location object (whitespace ignored). Raises ValueError for a
    location it cannot read."""
    ctx = _LocCtx(length, circular, stranded, where="location")
    loc_text = "".join(text.split())
    if _location_too_deep(loc_text):
        raise ValueError(f"cannot read location {text[:80]!r}: nested too deeply")
    try:
        return _location(loc_text, ctx)
    except (_LocationError, ValueError, IndexError, RecursionError) as exc:
        raise ValueError(f"cannot read location {text[:200]!r}") from exc


# ── LOCUS line ──────────────────────────────────────────────────────────────

_UNITS = (" bp ", " aa ", " rc ")
_TOPOLOGIES = ("linear", "circular")
_DATE_RE = re.compile(r"(\d{1,2})-([A-Za-z]{3})-(\d{4})")
_DIVISION_RE = re.compile(r"[A-Za-z]{3}")
_MOL_RE = re.compile(r"(?:ss-|ds-|ms-)?[A-Za-z]*(?:DNA|RNA|dna|rna)")


class _Rec:
    """What the scanner gathers for one record before it becomes one."""
    __slots__ = ("name", "size", "seq_type", "anns", "dbxrefs", "description",
                 "id", "refs", "cur_ref", "features", "repairs")

    def __init__(self):
        self.name: "str | None" = None
        self.size: "int | None" = None
        self.seq_type = ""
        self.anns: dict = {}
        self.dbxrefs: list = []
        self.description = ""
        self.id: "str | None" = None
        # None until the first REFERENCE: a record with none has no
        # `references` annotation at all (Biopython's reading).
        self.refs: "list | None" = None
        self.cur_ref: "Any" = None        # the Reference being filled in
        self.features: list = []
        self.repairs: "list[str]" = []

    def close_reference(self) -> None:
        """File the reference being filled in, if there is one."""
        if self.cur_ref is not None:
            if self.refs is None:
                self.refs = []
            self.refs.append(self.cur_ref)
            self.cur_ref = None

    def molecule_type(self, mol: str) -> None:
        if not mol:
            return
        if "circular" in mol or "linear" in mol:
            raise ValueError(f"Molecule type {mol!r} should not include topology")
        if mol[-3:].upper() in ("DNA", "RNA") and not mol[-3:].isupper():
            _warn(f"Non-upper case molecule type in LOCUS line: {mol}")
        self.anns["molecule_type"] = mol

    def topology(self, topo: str) -> None:
        if topo:
            if topo not in _TOPOLOGIES:
                raise ValueError(f"Unexpected topology {topo!r} should be "
                                 f"linear or circular")
            self.anns["topology"] = topo


def _locus_columns_old(line: str, r: _Rec) -> bool:
    """Biopython's 'old' fixed-column LOCUS layout. False when the line is
    not that layout (or fails its checks)."""
    if not (line[29:33] in _UNITS and line[55:62] == "       "):
        return False
    if (line[41:42] != " " or line[42:51].strip() not in ("", *_TOPOLOGIES)
            or line[51:52] != " "):
        return False
    name_len = " ".join(line[12:29].split()).split(" ")
    if len(name_len) != 2 or not _isnum(name_len[1]):
        return False
    parse_date = False
    if line[62:73].strip():
        parse_date = line[64:65] == "-" and line[68:69] == "-"
        if not parse_date:
            _warn("LOCUS line does not contain - at the expected position "
                  "in date:\n" + line)
    name, length = name_len
    if len(name) > 16:
        _warn("GenBank LOCUS line identifier over 16 characters")
    r.name = name
    r.size = int(length)
    if line[33:51].strip() == "" and line[29:33] == " aa ":
        r.seq_type = "PROTEIN"
    else:
        r.seq_type = line[33:51].strip()
    r.molecule_type(line[33:41].strip())
    r.topology(line[42:51].strip())
    r.anns["data_file_division"] = line[52:55]
    if parse_date:
        r.anns["date"] = line[62:73]
    return True


def _locus_columns_new(line: str, r: _Rec) -> bool:
    """Biopython's 'new' fixed-column LOCUS layout (the one NCBI writes).
    False when the line is not that layout or fails its checks — a JBEI-style
    line with no division and date fails them in Biopython too, after it pads
    the line, and is read by tokens instead."""
    if not (line[40:44] in _UNITS and line[54:64].strip() in ("", *_TOPOLOGIES)):
        return False
    if len(line) < 79:
        line = line + " " * (79 - len(line))
    if (line[44:47] not in ("   ", "ss-", "ds-", "ms-")
            or not (line[47:54].strip() == ""
                    or "DNA" in line[47:54].strip().upper()
                    or "RNA" in line[47:54].strip().upper())
            or line[54:55] != " "
            or line[55:63].strip() not in ("", *_TOPOLOGIES)
            or line[63:64] != " " or line[67:68] != " "):
        return False
    name_len = " ".join(line[12:40].split()).split(" ")
    if len(name_len) != 2 or not _isnum(name_len[1]):
        return False
    parse_date = False
    if line[68:79].strip():
        parse_date = line[70:71] == "-" and line[74:75] == "-"
        if not parse_date:
            _warn("LOCUS line does not contain - at the expected position "
                  "in date:\n" + line)
    r.name = name_len[0]
    r.size = int(name_len[1])
    if line[44:54].strip() == "" and line[40:44] == " aa ":
        r.seq_type = ("PROTEIN " + line[54:63]).strip()
    else:
        r.seq_type = line[44:63].strip()
    r.molecule_type(line[44:54].strip())
    r.topology(line[55:63].strip())
    if line[64:76].strip():
        r.anns["data_file_division"] = line[64:67]
    if parse_date:
        r.anns["date"] = line[68:79]
    return True


def _bio_token_reading(toks: "list[str]") -> "tuple | None":
    """``(name, size, topology, division, date)`` as Biopython's own
    token-based LOCUS branches read ``toks`` (None where it refuses)."""
    n = len(toks)
    if n == 8 and toks[3] in ("aa", "bp") and toks[5] in _TOPOLOGIES:
        return toks[1], toks[2], toks[5], toks[6], toks[7]
    if n == 7 and toks[3] in ("aa", "bp"):
        return toks[1], toks[2], "", toks[5], toks[6]   # topology lost
    if n >= 4 and toks[3] in ("aa", "bp"):
        return toks[1], toks[2], "", "", ""
    if n >= 4 and toks[-1] in ("aa", "bp"):
        return " ".join(toks[1:-2]), toks[-2], "", "", ""
    return None


def _locus_tokens(line: str, r: _Rec) -> None:
    """A LOCUS line read by what its fields ARE rather than where they sit —
    for every layout neither fixed-column one fits: a spaced name, a field
    one column over, a missing division or date, fields in another order.
    Counted as a repair (and warned about) only where this reading differs
    from what Biopython's own token rules make of the line."""
    all_toks = line.split()
    toks = all_toks[1:]
    unit = next((i for i in range(len(toks) - 1, 0, -1)
                 if toks[i].lower() in ("bp", "aa", "rc")
                 and _isnum(toks[i - 1])), None)
    if unit is None:
        if toks:
            r.name = " ".join(toks)
        else:
            _warn(f"Minimal LOCUS line found - is this correct?\n:{line!r}")
        return
    r.name = " ".join(toks[:unit - 1]) or None
    r.size = int(toks[unit - 1])
    mol = topo = div = date = ""
    dropped = []
    for t in toks[unit + 1:]:
        if t.lower() in _TOPOLOGIES and not topo:
            topo = t.lower()
        elif not mol and len(t) <= 7 and _MOL_RE.fullmatch(t):
            mol = t
        elif (not div and not date and _DIVISION_RE.fullmatch(t)
              and not _MOL_RE.fullmatch(t)):
            div = t.upper()
        elif not date and (m := _DATE_RE.fullmatch(t)):
            date = f"{int(m.group(1)):02d}-{m.group(2).upper()}-{m.group(3)}"
        else:
            dropped.append(t)
    if toks[unit].lower() == "aa" and not mol:
        r.seq_type = "PROTEIN"
    else:
        r.seq_type = (mol + " " + topo).strip()
    # The molecule type is then DERIVED from it at the end (`ds-DNA` ->
    # `DNA`), as Biopython does for every LOCUS line it reads by tokens.
    r.topology(topo)
    if div:
        r.anns["data_file_division"] = div
    if date:
        r.anns["date"] = date
    if _bio_token_reading(all_toks) != (r.name, str(r.size), topo, div, date):
        r.repairs.append("locus_line")
        _warn(f"LOCUS line read by its fields, not its columns:\n{line!r}")
    if dropped:
        _warn(f"LOCUS line: field(s) {' '.join(dropped)!r} not recognised")


def _locus(line: str, r: _Rec) -> None:
    if _locus_columns_old(line, r) or _locus_columns_new(line, r):
        return
    _locus_tokens(line, r)


# ── header ──────────────────────────────────────────────────────────────────

_HEADER_KEYS = {
    "DEFINITION", "ACCESSION", "NID", "PID", "DBSOURCE", "KEYWORDS",
    "SEGMENT", "SOURCE", "AUTHORS", "CONSRTM", "PROJECT", "TITLE", "JOURNAL",
    "MEDLINE", "PUBMED", "REMARK", "VERSION", "DBLINK", "REFERENCE",
    "ORGANISM", "COMMENT", "PRIMARY",
}
_SEQUENCE_HEADERS = ("CONTIG", "ORIGIN", "BASE COUNT", "WGS", "TSA", "TLS")
_LINEAGE_ROOTS = ("Bacteria.", "Archaea.", "Eukaryota.", "Unclassified.",
                  "Viruses.", "cellular organisms.", "other sequences.",
                  "unclassified sequences.")
_SC_START, _SC_END, _SC_DELIM = "-START##", "-END##", "::"
_SC_KEY_RE = re.compile(r"([^#]+)" + re.escape(_SC_START) + "$")
_SC_ITEM_RE = re.compile(r"(.+?)\s*::\s*(.*)")


def _split_accessions(s: str) -> list:
    return [x for x in s.replace("\n", " ").replace(";", " ").split() if x]


def _split_keywords(s: str) -> list:
    if s in ("", "."):
        s = ""
    elif s[-1] == ".":
        s = s[:-1]
    return [x.strip() for x in s.split(";")]


def _split_taxonomy(s: str) -> list:
    if not s or s == ".":
        return []
    if s[-1] == ".":
        s = s[:-1]
    return [item.strip() for tax in s.split(";") for item in tax.split("\n")
            if item]


def _header_keyword(line: str) -> "tuple[str, str]":
    """``(keyword, data)`` of a header line. The standard layout pads the
    keyword to column 12, and is read exactly as Biopython reads it; a line
    that does not (`DEFINITION A thing`) still gives its keyword, where
    Biopython ignores the line and loses its text."""
    kw = line[:12].strip()
    if kw in _HEADER_KEYS or not kw:
        return kw, line[12:].strip()
    tok = line.split(None, 1)
    if tok and tok[0] in _HEADER_KEYS:
        return tok[0], (tok[1].strip() if len(tok) > 1 else "")
    return kw, line[12:].strip()


def _is_cont(line: str) -> bool:
    """A header continuation line: no keyword. The standard is a blank first
    twelve columns; a line indented 6 to 11 and carrying no keyword still
    continues the field above it (a lineage line two columns short).

    Not one indented less: that is where a keyword sits (`  ORGANISM`,
    `   PUBMED`) and where a FEATURE line sits (`     source`). Read as a
    continuation, a misspelt sub-keyword or a feature table whose FEATURES
    line is missing poured into the field above it — the whole table into
    COMMENT. Biopython skips such a line, and so does this."""
    if line[:12] == "            ":
        return True
    if not line or line[0] != " ":
        return False
    body = line.lstrip(" ")
    if len(line) - len(body) < 6:
        return False
    tok = body.split(None, 1)
    return bool(tok) and tok[0] not in _HEADER_KEYS


def _cont_data(line: str) -> str:
    """A continuation line's text: from column 13 on, as Biopython takes it
    (keeping any further indent), or stripped when it is indented less."""
    return line[12:] if line[:12] == "            " else line.strip()


def _header(lines: "list[str]", r: _Rec) -> None:
    """The lines between LOCUS and FEATURES, as Biopython's
    `GenBankScanner._feed_header_lines` reads them."""
    lines = [x for x in lines if x]
    n = len(lines)
    i = 0
    anns = r.anns
    while i < n:
        line = lines[i]
        kw, data = _header_keyword(line)
        i += 1
        if kw == "VERSION":
            while "  " in data:
                data = data.replace("  ", " ")
            if " GI:" in data:
                _version(data.split(" GI:")[0], r)
                anns["gi"] = data.split(" GI:")[1]
            else:
                _version(data, r)
        elif kw == "DBLINK":
            cur = data.strip()
            while i < n and _is_cont(lines[i]):
                nxt = lines[i]
                i += 1
                if nxt.count(":") == 0:
                    cur += " " + nxt.strip()
                else:
                    _dblink(cur.strip(), r)
                    cur = nxt
            _dblink(cur.strip(), r)
        elif kw == "REFERENCE":
            data = data.strip()
            while i < n and _is_cont(lines[i]):
                data += " " + _cont_data(lines[i])
                i += 1
            while "  " in data:
                data = data.replace("  ", " ")
            r.close_reference()
            if r.refs is None:
                r.refs = []
            r.cur_ref = _model()["Reference"]()
            if " " in data:
                _reference_bases(data[data.find(" ") + 1:], r)
        elif kw == "ORGANISM":
            organism, lineage = data, ""
            while i < n and _is_cont(lines[i]):
                cont = _cont_data(lines[i])
                i += 1
                if (lineage or ";" in cont
                        or cont.strip() in _LINEAGE_ROOTS):
                    lineage += " " + cont
                elif cont.strip() == ".":
                    pass
                else:
                    organism += " " + cont.strip()
            anns["organism"] = organism
            tax = _split_taxonomy(lineage.strip())
            if "taxonomy" in anns:
                anns["taxonomy"].extend(tax)
            else:
                anns["taxonomy"] = tax
        elif kw == "COMMENT":
            _comment(line, lines, i, n, r)
            while i < n and _is_cont(lines[i]):
                i += 1
        elif kw in _SIMPLE_HEADERS:
            while i < n and _is_cont(lines[i]):
                data += " " + _cont_data(lines[i])
                i += 1
            if kw == "DEFINITION" and data.endswith("."):
                data = data[:-1]
            _SIMPLE_HEADERS[kw](data, r)
        else:
            # PRIMARY and anything unknown: skip it and its continuation.
            while i < n and _is_cont(lines[i]):
                i += 1


def _comment(line: str, lines: "list[str]", i: int, n: int, r: _Rec) -> None:
    if line[:12].strip() == "COMMENT":
        data = line[12:]
    else:
        tok = line.split(None, 1)
        data = tok[1] if len(tok) > 1 else ""
    comments: list = []
    structured: dict = {}
    key = None
    m = _SC_KEY_RE.search(data)
    if m is not None:
        key = m.group(1)
    else:
        comments.append(data)
    last = None
    while i < n and _is_cont(lines[i]):
        raw = lines[i]
        data = _cont_data(raw)
        i += 1
        if _SC_START in data:
            m = _SC_KEY_RE.search(data)
            if m is not None:
                key = m.group(1)
            else:
                comments.append(data)
        elif key is not None and _SC_DELIM in data:
            m = _SC_ITEM_RE.search(data)
            if m is not None:
                # Each value as its lines, joined once at the end: a value
                # continued over many lines must not cost quadratic time.
                structured.setdefault(key, {})[m.group(1)] = [m.group(2)]
                last = m.group(1)
        elif key is not None and _SC_END not in data:
            if last is None or last not in structured.get(key, {}):
                # (Also a continuation in a block reopened after another one,
                # whose last item belongs to that other block: Biopython
                # raises KeyError there.)
                _warn("Structured comment not parsed on malformed header "
                      f"line: {raw}")
                continue
            structured[key][last].append(raw.strip())
        elif _SC_END in data:
            key = None
        else:
            comments.append(data)
    if comments:
        if "comment" in r.anns:
            r.anns["comment"] += "\n" + "\n".join(comments)
        else:
            r.anns["comment"] = "\n".join(comments)
    if structured:
        r.anns["structured_comment"] = {
            k: {item: " ".join(parts) for item, parts in d.items()}
            for k, d in structured.items()}


def _version(v: str, r: _Rec) -> None:
    if v.count(".") == 1 and _isnum(v.split(".")[1]):
        _accession(v.split(".")[0], r)
        r.anns["sequence_version"] = int(v.split(".")[1])
    elif v:
        r.id = v


def _accession(s: str, r: _Rec) -> None:
    accs = _split_accessions(s)
    if "accessions" in r.anns:
        for a in accs:
            if a not in r.anns["accessions"]:
                r.anns["accessions"].append(a)
    elif accs:
        # An EMPTY accession line records nothing (Biopython's reading).
        r.anns["accessions"] = accs
    if not r.id and accs:
        r.id = r.anns["accessions"][0]


def _dblink(s: str, r: _Rec) -> None:
    while ": " in s:
        s = s.replace(": ", ":")
    if s.strip() not in r.dbxrefs:
        r.dbxrefs.append(s.strip())


def _reference_bases(content: str, r: _Rec) -> None:
    if not content.endswith(")") or not content.startswith("("):
        _warn(f"REFERENCE line not understood: {content!r}")
        return
    info = content[1:-1]
    locs = []
    try:
        if "bases" in info and "to" in info:
            locs = _ref_locations(info[5:])
        elif "residues" in info and "to" in info:
            locs = _ref_locations(info[info.find("residues") + 9:])
        elif not (info == "sites" or info.strip() == "bases"):
            raise ValueError(info)
    except ValueError:
        _warn(f"REFERENCE base range not understood: {content!r}")
        locs = []
    r.cur_ref.location = locs


def _ref_locations(s: str) -> list:
    SL = _model()["SimpleLocation"]
    out = []
    for part in s.split(";"):
        a, b = part.split("to")
        out.append(SL(int(a.strip()) - 1, int(b.strip())))
    return out


def _ref_attr(name: str):
    def _set(content: str, r: _Rec) -> None:
        ref = r.cur_ref
        if ref is None:
            if name == "title":
                _warn("GenBank TITLE line without REFERENCE line.")
            return
        cur = getattr(ref, name)
        setattr(ref, name, cur + " " + content if cur else content)
    return _set


def _ref_id(name: str):
    def _set(content: str, r: _Rec) -> None:
        if r.cur_ref is not None:
            setattr(r.cur_ref, name, content)
    return _set


def _set_ann(key: str, fn=None):
    def _set(content: str, r: _Rec) -> None:
        r.anns[key] = fn(content) if fn else content
    return _set


def _definition(s: str, r: _Rec) -> None:
    r.description = r.description + " " + s if r.description else s


def _keywords(s: str, r: _Rec) -> None:
    if "keywords" in r.anns:
        r.anns["keywords"].extend(_split_keywords(s))
    else:
        r.anns["keywords"] = _split_keywords(s)


def _source(s: str, r: _Rec) -> None:
    r.anns["source"] = s[:-1] if s and s[-1] == "." else s


def _project(s: str, r: _Rec) -> None:
    r.dbxrefs.extend(p for p in s.replace("GenomeProject:", "Project:").split()
                     if p)


_SIMPLE_HEADERS = {
    "DEFINITION": _definition,
    "ACCESSION": _accession,
    "NID": _set_ann("nid"),
    "PID": _set_ann("pid"),
    "DBSOURCE": _set_ann("db_source", str.rstrip),
    "KEYWORDS": _keywords,
    "SEGMENT": _set_ann("segment"),
    "SOURCE": _source,
    "AUTHORS": _ref_attr("authors"),
    "CONSRTM": _ref_attr("consrtm"),
    "PROJECT": _project,
    "TITLE": _ref_attr("title"),
    "JOURNAL": _ref_attr("journal"),
    "MEDLINE": _ref_id("medline_id"),
    "PUBMED": _ref_id("pubmed_id"),
    "REMARK": _ref_attr("comment"),
}


# ── feature table ───────────────────────────────────────────────────────────

_KEY_CHARS_RE = re.compile(r"[A-Za-z0-9_'*-]+")
_LOC_START_RE = re.compile(
    r"(?:complement|join|order|bond|one-of)\(|[<>?(]?\d|[A-Za-z][\w.|]*:\S")
_KNOWN_KEYS = frozenset("""
-10_signal -35_signal 3'clip 3'UTR 5'clip 5'UTR assembly_gap attenuator
C_region CAAT_signal CDS centromere conflict D-loop D_segment enhancer exon gap
GC_signal gene iDNA intron J_segment LTR mat_peptide misc_binding
misc_difference misc_feature misc_recomb misc_RNA misc_signal misc_structure
mobile_element modified_base mRNA ncRNA N_region old_sequence operon oriT
polyA_signal polyA_site precursor_RNA prim_transcript primer_bind promoter
propeptide protein_bind RBS regulatory repeat_region repeat_unit rep_origin
rRNA S_region satellite scRNA sig_peptide snoRNA snRNA source stem_loop STS
TATA_signal telomere terminator tmRNA transit_peptide tRNA unsure V_region
V_segment variation Protein Region Site Bond SecStr Het Non_std_residue
""".split())


class _Feat:
    __slots__ = ("key", "loc", "quals")

    def __init__(self, key: str, loc: str):
        self.key = key
        # The location's lines, joined once the feature is complete: text
        # accumulated with `+=` costs quadratic time on a long location.
        self.loc: "list[str]" = [loc] if loc else []
        # One `[key, lines | None, state]` per qualifier: `lines` is the
        # value's text line by line (None for a flag such as `/pseudo`);
        # `state` is True while its quotes are still open, "unclosed" when
        # they never closed.
        self.quals: list = []


def _feature_table(lines: "list[str]", r: _Rec) -> "list[_Feat]":
    """Group the feature table into features — by what each line IS rather
    than by column. On a file laid out to the standard (keys from column 6,
    everything else from column 22) the grouping is exactly Biopython's;
    tabs, a shifted line or a mis-indented qualifier no longer turn into a
    garbage feature key (`/organism="Cloning`) or a crash.

    A line is read as a new feature when it is indented like one (nearer
    this table's own key column than its qualifier column) and starts with
    a word a feature key can be; as a qualifier when it starts with `/`;
    otherwise it continues the location or the last qualifier. Inside a
    quoted value every line continues the value, unless it is plainly the
    next feature (a known key followed by a location) — a quote that never
    closes then ends there, with a warning, instead of swallowing the rest
    of the table.

    The key column is learnt from the first line that IS a feature, never
    from a stray qualifier above it: a table whose first key line was lost
    used to take its qualifier column for the key column, and every wrapped
    line after that became a feature of its own, taking the text it held
    with it."""
    feats: "list[_Feat]" = []
    cur = None
    k_ind = q_ind = None
    open_q = None                      # the qualifier still inside its quotes
    orphan_open = False                # inside the quotes of a stray qualifier
    in_header = False                  # inside header text within the table
    differs = False                    # read differently from fixed columns
    for raw in lines:
        if not raw:
            continue
        tabs = "\t" in raw
        if tabs:
            raw = raw.expandtabs(8)
        s = raw.lstrip()
        if not s:
            continue
        indent = len(raw) - len(s)
        col_key = tabs or bool(raw[2:21].strip())   # the fixed-column reading
        if open_q is not None:
            if not (indent < 21 and _is_feature_start(s)):
                open_q[1].append(s)
                if s[-1] == '"':
                    open_q[2] = False
                    open_q = None
                differs = differs or col_key
                continue
            r.repairs.append("unclosed_quote")
            _warn(f"/{open_q[0]} never closes its quotes — ended at the next "
                  f"feature")
            open_q[2] = "unclosed"
            open_q = None
        if orphan_open:
            if not (indent < 21 and _is_feature_start(s)):
                orphan_open = s[-1] != '"'
                continue
            orphan_open = False
        if in_header:
            # Header continuation lines sit at column 13; a line nearer the
            # key column shaped like `key location` is the table again —
            # whatever its key (a custom one too, not only a known one).
            if not (indent < 12 and _looks_like_feature_line(s)):
                continue
            in_header = False
        if s[0] != "/":
            # Junk in front of a qualifier (`"   /codon_start=1`, a stray
            # character where the indent was): the qualifier is the line.
            # Biopython reads the qualifier column and never sees the junk.
            sp = s.find(" ")
            if sp > 0 and s[sp:].lstrip().startswith("/") and (
                    cur is not None or not _KEY_CHARS_RE.fullmatch(s[:sp])):
                _warn(f"Text in front of a qualifier dropped: {s[:sp]!r}")
                s = s[sp:].lstrip()
                # Its column is where the qualifier starts, not where the
                # junk did: the first qualifier line sets the column every
                # later line is judged by.
                indent = len(raw) - len(s)
                differs = True
        if s[0] == "/":
            if cur is None:
                _warn(f"Text before the first feature in the feature table: "
                      f"{s[:80]!r}")
                orphan_open = bool(_qualifier_start(s)[2])
                continue
            if q_ind is None:
                q_ind = indent
            q = _qualifier_start(s)
            cur.quals.append(q)
            if q[2]:
                open_q = q
            differs = differs or col_key
            continue
        sp = s.find(" ")
        tok = s if sp < 0 else s[:sp]
        if cur is None:
            # The first feature: any key-like word not in the qualifier
            # column (a table may be shifted or tab-indented as a whole).
            is_key = indent < 21 and _KEY_CHARS_RE.fullmatch(tok) is not None
        else:
            k = k_ind if k_ind is not None else 5   # set by the first feature
            threshold = (k + (q_ind if q_ind is not None else k + 16)) // 2 + 1
            is_key = indent < threshold and _KEY_CHARS_RE.fullmatch(tok) is not None
        if is_key and indent == 0 and tok in _HEADER_KEYS:
            # A header keyword inside the table (`COMMENT` written after the
            # features): it and its continuation lines are text, never a
            # feature, up to the next line that plainly starts one.
            _warn(f"Header line inside the feature table skipped: {s[:80]!r}")
            in_header = True
            continue
        if is_key:
            if k_ind is None:
                k_ind = indent
            differs = differs or not col_key
            cur = _Feat(tok, "" if sp < 0 else s[sp:].strip())
            feats.append(cur)
            continue
        differs = differs or col_key
        if cur is None:
            _warn(f"Text before the first feature in the feature table: {s[:80]!r}")
            continue
        if not cur.quals:
            cur.loc.append(s)
            continue
        last = cur.quals[-1]
        if last[1] is None:
            _warn(f"/{last[0]} takes no value, but text follows it: {s[:80]!r}")
            last[1] = [s]
        else:
            last[1].append(s)
    if open_q is not None:
        r.repairs.append("unclosed_quote")
        _warn(f"/{open_q[0]} never closes its quotes")
        open_q[2] = "unclosed"
    if differs:
        r.repairs.append("feature_table_layout")
        _warn("Feature table not in the standard columns — read by layout")
    return feats


def _is_feature_start(s: str) -> bool:
    """Plainly the next feature: a KNOWN key followed by a location. Strict on
    purpose — it is what ends an unclosed quote, where a wrongly early end
    would cut a value short."""
    tok = s.split(None, 1)
    return (len(tok) == 2 and tok[0] in _KNOWN_KEYS
            and _LOC_START_RE.match(tok[1]) is not None)


def _looks_like_feature_line(s: str) -> bool:
    """Any key-shaped word followed by a location (custom keys included)."""
    tok = s.split(None, 1)
    return (len(tok) == 2 and _KEY_CHARS_RE.fullmatch(tok[0]) is not None
            and _LOC_START_RE.match(tok[1]) is not None)


def _qualifier_start(s: str) -> list:
    """``[key, [value] | None, still_open]`` for a line starting ``/``,
    exactly as Biopython's `InsdcScanner.parse_feature` reads one — except
    where the `=` is missing or a space sits before it (`/product "text…`,
    `/note ="x"`): a qualifier name never holds a space, so the first word is
    the name and the rest its value. Biopython makes the whole line one
    garbage name (`product "text…`), then refuses the file at its
    continuation."""
    i = s.find("=")
    sp = s.find(" ")
    if i == -1 or 0 <= sp < i:
        if sp < 0:
            return [s[1:], None, False]
        key, value = s[1:sp], s[sp + 1:].lstrip()
        if value.startswith("="):
            value = value[1:].lstrip()
        _warn(f"No '=' straight after /{key} — read the rest of the line as "
              f"its value")
        if not value or value == '"' or value[0] != '"':
            return [key, [value], False]
        return [key, [value], value[-1] != '"']
    key, value = s[1:i], s[i + 1:]
    if value.startswith(" ") and value.lstrip().startswith('"'):
        _warn("White space after equals in qualifier")
        value = value.lstrip()
    if not value or value == '"' or value[0] != '"':
        return [key, [value], False]
    return [key, [value], value[-1] != '"']


def _qualifier_value(key: str, value: str) -> str:
    """A qualifier's text as Biopython's `_FeatureConsumer.feature_qualifier`
    leaves it: line breaks to spaces, enclosing quotes off, `""` to `"`,
    every whitespace out of a `/translation`."""
    value = value.replace("\n", " ")
    if len(value) > 1 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]
    if '"' in value:
        if re.search(r'[^"]"[^"]|^"[^"]|[^"]"$', value):
            _warn('The NCBI states double-quote characters like " should be '
                  'escaped as "" (two double - quotes), but here it was not: '
                  f"{value!r}")
        value = value.replace('""', '"')
    if key == "translation":
        value = "".join(value.split())
    return value


def _build_features(feats: "list[_Feat]", r: _Rec) -> list:
    M = _model()
    SeqFeature = M["SeqFeature"]
    circular = "circular" in r.anns.get("topology", "").lower()
    ctx = _LocCtx(r.size, circular, "PROTEIN" not in r.seq_type.upper())
    out = []
    for f in feats:
        loc_text = "".join("".join(f.loc).split())
        if "replace" in loc_text:
            comma = loc_text.find(",")
            loc_text = loc_text[8:comma]
        ctx.where = f"{f.key} feature {loc_text[:200]!r}"
        try:
            if _location_too_deep(loc_text):
                raise _LocationError("nested too deeply")
            loc = _location(loc_text, ctx)
        except (_LocationError, ValueError, IndexError, RecursionError):
            r.repairs.append("unreadable_location")
            _warn(f"{f.key} feature: location {loc_text[:200]!r} cannot be "
                  f"read — the feature is left out")
            continue
        feat = SeqFeature(loc, type=f.key)
        quals: dict = {}
        for key, lines, state in f.quals:
            if lines is None:
                if key not in quals:
                    quals[key] = [""]
                continue
            value = "\n".join(lines)
            if state == "unclosed" and value[:1] == '"':
                value = value[1:]
            v = _qualifier_value(key, value)
            if key in quals:
                quals[key].append(v)
            else:
                quals[key] = [v]
        feat.qualifiers = quals
        out.append(feat)
    r.repairs.extend(ctx.repairs)
    return out


# ── footer + sequence ───────────────────────────────────────────────────────

_SEQ_DELETE = b"0123456789 \t\n\r\x0b\x0c"
_SEQ_DELETE_STR = str.maketrans("", "", "0123456789 \t\n\r\x0b\x0c")


def _bases(block: str) -> str:
    """The bases of an ORIGIN block: everything but the position numbers and
    the white space. (Biopython cuts the first ten COLUMNS of every line, so a
    line written without numbers lost its first ten bases.)

    A second ORIGIN or BASE COUNT line inside the block is not bases. A
    no-break or thin space (pasted from a web page or a word processor) is
    white space like any other; any other non-ASCII character is not a
    sequence letter, and is refused by name — `Seq` cannot hold it, and its
    own error ("'ascii' codec can't encode…") said nothing useful."""
    if "\nORIGIN" in block or "\nBASE COUNT" in block:
        block = "\n".join(ln for ln in block.split("\n")
                          if not ln.startswith(("ORIGIN", "BASE COUNT")))
    try:
        return block.encode("ascii").translate(None, _SEQ_DELETE).decode(
            "ascii").upper()
    except UnicodeEncodeError:
        pass
    kept = "".join(ch for ch in block.translate(_SEQ_DELETE_STR)
                   if not (ch.isspace() or ch.isdigit()))
    bad = next((ch for ch in kept if not ch.isascii()), None)
    if bad is not None:
        raise ValueError(f"the sequence holds {bad!r}, which is not a "
                         f"sequence letter — the file is damaged")
    return kept.upper()


def _misc(lines: "list[str]", r: _Rec) -> None:
    """CONTIG / WGS / TSA / TLS / BASE COUNT lines between the feature table
    and the sequence."""
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        i += 1
        if line.startswith("TLS "):
            r.anns["tls"] = line[3:].strip().split("-")
        elif line.startswith("TSA "):
            r.anns["tsa"] = line[3:].strip().split("-")
        elif line.startswith("WGS_SCAFLD"):
            r.anns.setdefault("wgs_scafld", []).append(
                line[10:].strip().split("-"))
        elif line.startswith("WGS "):
            r.anns["wgs"] = line[3:].strip().split("-")
        elif line.startswith("CONTIG"):
            contig = line[6:].strip()
            while i < n and lines[i][:12] == "            ":
                contig += lines[i][12:].rstrip()
                i += 1
            r.anns["contig"] = contig


# ── records ─────────────────────────────────────────────────────────────────

def _record_starts(text: str) -> "list[int]":
    starts = []
    i = text.find("LOCUS")
    while i >= 0:
        if (i == 0 or text[i - 1] == "\n") and (
                i + 5 >= len(text) or text[i + 5] in " \t\n"):
            starts.append(i)
        i = text.find("LOCUS", i + 5)
    return starts


def _normalise(text: str) -> str:
    """Text as the reader takes it: no BOM in front, `\\n` line ends. The file
    loaders decode with `_read_text_tolerant`, which does both; this covers
    TEXT handed in directly (an agent `gb_text`, a paste), where a BOM hid the
    one record and a CR-only text was one enormous LOCUS line."""
    if text[:1] == "\ufeff":
        text = text[1:]
    if "\r" in text:
        if "\n" in text:
            # CRLF line ends. A lone CR INSIDE a line is data — a `.dna` note
            # typed with one ("first\rsecond") — and Biopython reads it as
            # data too; turned into a line break, it changed the note.
            text = text.replace("\r\n", "\n")
        else:
            text = text.replace("\r", "\n")      # CR-only (classic Mac)
    return text


def _record_spans(text: str) -> "list[tuple[int, int]]":
    """``(start, end)`` of each record in normalised ``text``. A LOCUS line
    with nothing under it before the next LOCUS line (a duplicated or
    stray first line) starts no record of its own: it used to turn a
    one-plasmid file into a two-record file whose first record had no
    sequence."""
    starts = _record_starts(text)
    spans = [(s, starts[k + 1] if k + 1 < len(starts) else len(text))
             for k, s in enumerate(starts)]
    if len(spans) < 2:
        return spans
    kept = []
    for s, e in spans:
        nl = text.find("\n", s, e)
        if nl < 0 or not text[nl:e].strip():
            _warn(f"A LOCUS line with no record under it was skipped: "
                  f"{text[s:nl if nl >= 0 else e][:80]!r}")
            continue
        kept.append((s, e))
    return kept or spans[-1:]


def genbank_record_count(text: str) -> int:
    """How many records `parse_genbank` would find in ``text``, without
    reading them."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", InsdcParserWarning)
        return len(_record_spans(_normalise(text)))


def _is_seq_header(line: str) -> bool:
    return (line[:12].rstrip() in _SEQUENCE_HEADERS
            or line.startswith(("ORIGIN", "CONTIG", "BASE COUNT")))


def _one(text: str, start: int, end: int):
    """Parse the record in ``text[start:end]``."""
    M = _model()
    r = _Rec()
    # The `//` that ends this record.
    term = text.find("\n//", start, end)
    body_end = term if term >= 0 else end
    if term < 0:
        _warn("Premature end of file in sequence data")
    # Where the sequence starts: the ORIGIN line. Everything before it is
    # split into lines; the bases never are.
    origin = text.find("\nORIGIN", start, body_end)
    head_end = origin if origin >= 0 else body_end
    lines = text[start:head_end].split("\n")
    # NOT right-stripped: the fixed-column layouts are told apart by the
    # spaces a writer pads the line with.
    _locus(lines[0].rstrip("\r"), r)
    hdr: "list[str]" = []
    ft: "list[str]" = []
    misc: "list[str]" = []
    where = hdr
    for raw in lines[1:]:
        line = raw.rstrip()
        if where is hdr:
            if line[:8].upper() == "FEATURES" and (
                    len(line) == 8 or line[8] in " \t"):
                # Any spacing before `Location/Qualifiers`, any case: the
                # line that starts the table is the table's start.
                where = ft
                continue
            if _is_seq_header(line):
                where = misc
                misc.append(line)
                continue
            hdr.append(line)
        elif where is ft:
            if line and line[0] != " " and line[0] != "\t" and _is_seq_header(
                    line):
                where = misc
                misc.append(line)
                continue
            if line == "//":
                break
            ft.append(line)
        else:
            misc.append(line)
    _header(hdr, r)
    r.close_reference()
    if r.refs is not None:
        r.anns["references"] = r.refs
    feats = _feature_table(ft, r)
    seq = ""
    if origin >= 0:
        o_end = text.find("\n", origin + 1, body_end)
        o_end = body_end if o_end < 0 else o_end
        block = text[o_end:body_end]
        contig_at = block.find("\nCONTIG")
        if contig_at >= 0:
            misc.extend(block[contig_at + 1:].split("\n"))
            block = block[:contig_at]
        seq = _bases(block)
        if seq and _has_unnumbered_lines(block):
            r.repairs.append("unnumbered_sequence")
            _warn("Sequence lines without position numbers — every base read")
    _misc(misc, r)
    if term < 0 and origin < 0 and not any(
            k in r.anns for k in ("contig", "wgs", "wgs_scafld", "tsa", "tls")):
        # No sequence lines of any kind (ORIGIN, CONTIG, WGS / TSA / TLS) and
        # no `//`: the text stops before the record gets to its sequence. Biopython refuses it ("Premature end of line
        # during features table"); read anyway, it came out as a CONTIG-style
        # record of the declared length with no bases, and the error the user
        # then saw blamed a CONTIG record that was never there.
        named = f" {r.name}" if r.name else ""
        raise ValueError(
            f"the record{named} ends before its sequence (no ORIGIN line and "
            f"no //) — the file looks cut short")
    r.features = _build_features(feats, r)
    return _finish(r, seq, M)


def _has_unnumbered_lines(block: str) -> bool:
    """Does any sequence line lack its position number? Checked on the first
    few lines only: a writer numbers every line or none."""
    n = 0
    for ln in block.split("\n", 6)[:6]:
        s = ln.strip()
        if not s:
            continue
        n += 1
        if not s[0].isdigit():
            return True
    return False


def _finish(r: _Rec, seq: str, M: dict):
    SeqRecord = M["SeqRecord"]
    Seq = M["Seq"]
    rec = SeqRecord(None, id=None)
    rec.name = r.name if r.name else "<unknown name>"
    rec.description = r.description
    rec.dbxrefs = r.dbxrefs
    rid = r.id
    if not rid:
        rid = r.name if r.name else "<unknown id>"
    elif rid.count(".") == 0 and "sequence_version" in r.anns:
        rid += ".%i" % r.anns["sequence_version"]
    rec.id = rid
    if r.size is not None and seq and r.size != len(seq):
        _warn(f"Expected sequence length {r.size}, found {len(seq)} ({rid}).")
    st = r.seq_type.upper()
    mol = None
    if r.seq_type:
        if "DNA" in st or "MRNA" in st:
            mol = "DNA"
        elif "RNA" in st:
            mol = "RNA"
        elif "PROTEIN" in st or r.seq_type == "PRT":
            mol = "protein"
        elif r.seq_type not in ("circular", "linear", "unspecified"):
            _warn(f"Could not determine molecule_type for seq_type "
                  f"{r.seq_type!r}")
    if mol is not None:
        r.anns["molecule_type"] = r.anns.get("molecule_type", mol)
    rec.annotations = r.anns
    rec.features = r.features
    if not seq and r.size:
        rec.seq = Seq(None, length=r.size)
    else:
        rec.seq = Seq(seq)
    counts: dict = {}
    for k in r.repairs:
        counts[k] = counts.get(k, 0) + 1
    rec._insdc_repairs = counts       # type: ignore[attr-defined]
    return rec


def repair_note(rec) -> "str | None":
    """One sentence for the user about a repair that changed what they get —
    None when nothing worth interrupting them for happened. Layout repairs
    (LOCUS fields, a tab-indented table, unnumbered sequence lines, an
    unclosed quote) lose nothing and are only logged; a feature left out, or
    kept across the origin of a molecule the file calls linear, is worth a
    look before the file is trusted."""
    rep = getattr(rec, "_insdc_repairs", None) or {}
    parts = []
    n = int(rep.get("unreadable_location", 0) or 0)
    if n:
        parts.append(f"{n} feature{'s' if n != 1 else ''} had a location that "
                     f"could not be read and {'were' if n != 1 else 'was'} "
                     f"left out")
    n = int(rep.get("origin_spanning_not_circular", 0) or 0)
    if n:
        parts.append(f"{n} feature{'s' if n != 1 else ''} "
                     f"run{'' if n != 1 else 's'} across the origin, but the "
                     f"file marks the molecule linear")
    if not parts:
        return None
    return "; ".join(parts) + " (details in the log)."


def parse_genbank(text: str) -> list:
    """Every GenBank record in ``text``, in order. A text with no LOCUS line
    has no records (an empty list, not an error)."""
    if not isinstance(text, str):
        raise TypeError("parse_genbank wants text")
    text = _normalise(text)
    return [_one(text, s, e) for s, e in _record_spans(text)]


def read_genbank(text: str):
    """The ONE record in ``text``; ValueError for none or several."""
    recs = parse_genbank(text)
    if not recs:
        raise ValueError("No records found in handle")
    if len(recs) > 1:
        raise ValueError("More than one record found in handle")
    return recs[0]
