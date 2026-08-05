"""CommercialSaaS .dna history model + provenance date helpers (layer 2).

`_CommercialSaaSHistoryNode` (one node in the .dna <HistoryTree> lineage) plus
the pure history helpers it and the codec share: `_HISTORY_MONTHS`,
`_history_now_str` (machine storage stamp), `_history_human_dt` (the universal
slash-free "JUN 9 2026 14:30" render), and `_coerce_int_or_zero`.

Self-contained: ElementTree and datetime are imported function-locally; depends
on no other splicecraft module. The .dna codec FUNCTIONS that build and
serialise these trees stay in the hub (they couple to file I/O).
"""
from __future__ import annotations

from rich.text import Text

import re


# Hardcoded (NOT `strftime('%b')`, which is locale-dependent) so the History
# tab's dates render the same "JUN" everywhere.
_HISTORY_MONTHS: tuple[str, ...] = (
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
    "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
)


def _history_now_str() -> str:
    """Current LOCAL date+time as a compact, sortable storage stamp
    (``YYYY-MM-DDTHH:MM``) for a history node's ``date`` attribute. Stored
    machine-clean and rendered to the universal, slash-free "JUN 9 2026 14:30"
    by `_history_human_dt` (no MM/DD-vs-DD/MM ambiguity — user request
    2026-06-10). Stamped by the fresh-build callers, mirroring how
    `_build_commercialsaas_notes_packet` calls `datetime.now()` itself rather
    than baking non-determinism into the `_CommercialSaaSHistoryNode.new`
    constructor."""
    from datetime import datetime
    return datetime.now().strftime("%Y-%m-%dT%H:%M")


def _history_human_dt(stamp: str) -> str:
    """Render a stored history timestamp as the universal "JUN 9 2026 14:30"
    (uppercase 3-letter month, no leading-zero day, 24-hour clock). Accepts
    our ISO storage form AND the CommercialSaaS Notes ``YYYY.MM.DD`` date, so
    an imported file-level date renders too. Returns "" on empty/unparseable
    input — an undated node (reconstructed lineage / starting material whose
    real date is unknown) shows no date rather than a bogus one."""
    # `str(...)` + clamp: the `date` attr is always a string off XML, but a
    # hostile/garbage history_xml could carry an over-long or odd value, and a
    # direct caller might pass a non-string. The regex is start-anchored with
    # bounded quantifiers (no ReDoS); the clamp just caps the strip/scan cost.
    s = str(stamp or "").strip()[:40]
    if not s:
        return ""
    m = re.match(r"(\d{4})[-./](\d{1,2})[-./](\d{1,2})"
                 r"(?:[T ](\d{1,2}):(\d{2}))?", s)
    if not m:
        return ""
    year, mon, day = int(m.group(1)), int(m.group(2)), int(m.group(3))
    if not (1 <= mon <= 12 and 1 <= day <= 31):
        return ""
    out = f"{_HISTORY_MONTHS[mon - 1]} {day} {year}"
    if m.group(4) is not None:
        hh, mm = int(m.group(4)), int(m.group(5))
        if 0 <= hh <= 23 and 0 <= mm <= 59:   # drop a garbage time, keep date
            out += f" {hh:02d}:{mm:02d}"
    return out


# Node attributes with a typed getter on `_CommercialSaaSHistoryNode`.
# Anything NOT here surfaces through `extra_attributes` → the detail pane's
# "Other recorded fields" block, so an attribute this codebase has never seen
# is still SHOWN rather than silently dropped. Adding a typed getter means
# adding its attribute name here (otherwise the field renders twice).
_HISTORY_MODELLED_ATTRS: "frozenset[str]" = frozenset({
    "name", "type", "seqLen", "strandedness", "ID", "circular",
    "operation", "date", "resurrectable",
    "upstreamModification", "downstreamModification",
    "upstreamStickiness", "downstreamStickiness", "countBasesFrom",
    "customMapLabel", "useCustomMapLabel", "hideHistory",
})

# Child element tags with a dedicated renderer (or deliberately reported as a
# count). Anything else lands in `other_children`. `Features`/`HistoryColors`/
# `StrandColors` are the cosmetic annotation+colour snapshots — real, but bulk
# (~95% of a history's bytes), so they're summarised, not enumerated.
_HISTORY_MODELLED_CHILDREN: "frozenset[str]" = frozenset({
    "Node", "InputSummary", "RegeneratedSite", "Oligo", "Parameter",
    "Primers", "Features", "HistoryColors", "StrandColors",
    "UpstreamModification", "DownstreamModification",
})


class _CommercialSaaSHistoryNode:
    """One node in the CommercialSaaS `<HistoryTree>`.

    Wraps an `xml.etree.ElementTree.Element` so unknown attributes
    and child elements survive a parse → modify → serialise cycle.
    Typed properties expose the well-known fields; raw access is
    available via ``self.element`` for anything else.

    Construction modes:
      * ``_CommercialSaaSHistoryNode(element)`` — wrap an existing element
        (used by the parser).
      * ``_CommercialSaaSHistoryNode.new(name=..., seq_len=..., …)`` —
        create a fresh node with the canonical attributes set,
        ready to attach to a parent or use as the tree root.
    """

    __slots__ = ("element", "_sibling_elements")

    def __init__(self, element) -> None:
        self.element = element
        # 2026-05-27 (audit-3 M5): optional sibling top-level
        # `<Node>` elements preserved by `_parse_commercialsaas_history`
        # when the source file carried more than one. Empty list by
        # default; the serializer round-trips them back into the
        # output `<HistoryTree>`.
        self._sibling_elements: list = []

    @classmethod
    def new(cls, *, name: str, seq_len: int, circular: bool,
              operation: str, node_id: int = 0,
              strandedness: str = "double",
              date: str = "") -> "_CommercialSaaSHistoryNode":
        import xml.etree.ElementTree as _ET
        el = _ET.Element("Node")
        el.set("name", str(name))
        el.set("type", "DNA")
        el.set("seqLen", str(int(seq_len)))
        el.set("strandedness", str(strandedness))
        el.set("ID", str(int(node_id)))
        el.set("circular", "1" if circular else "0")
        el.set("operation", str(operation))
        # Optional creation timestamp (our ISO storage form, see
        # `_history_now_str`): fresh-build callers stamp it so the History tab
        # shows WHEN each step happened; reconstructed lineage / starting
        # material passes "" (real date unknown). Round-trips like any other
        # attr — and SnapGene/CommercialSaaS tolerates the extra attribute.
        if date:
            el.set("date", str(date))
        return cls(el)

    # ── Typed getters ────────────────────────────────────────────

    @property
    def name(self) -> str:
        return self.element.get("name", "") or ""

    @property
    def operation(self) -> str:
        return self.element.get("operation", "") or ""

    @property
    def date(self) -> str:
        """Creation timestamp for this step — our ISO storage form, or a
        CommercialSaaS Notes-style ``YYYY.MM.DD`` date carried in on import.
        Empty for reconstructed lineage / starting material with no recorded
        date; `_history_human_dt` renders it to "JUN 9 2026 14:30"."""
        return self.element.get("date", "") or ""

    @property
    def seq_len(self) -> int:
        try:
            return int(self.element.get("seqLen", "0"))
        except (TypeError, ValueError):
            return 0

    @property
    def circular(self) -> bool:
        return self.element.get("circular") == "1"

    @property
    def node_id(self) -> int:
        try:
            return int(self.element.get("ID", "0"))
        except (TypeError, ValueError):
            return 0

    @property
    def resurrectable(self) -> bool:
        """CommercialSaaS marks parent nodes as ``resurrectable="1"`` when
        the original fragment can be re-extracted from the history
        (i.e., a downstream user could reconstruct the parent
        plasmid). Defaults to False for the top-level result node.

        Any NON-ZERO value counts. Real files also carry
        ``resurrectable="2"``; the old ``== "1"`` test read those as
        False, hiding the badge on genuinely resurrectable ancestors."""
        return _coerce_int_or_zero(self.element.get("resurrectable")) != 0

    @property
    def parents(self) -> "list[_CommercialSaaSHistoryNode]":
        """Direct parent fragments (one level down). Use
        :meth:`walk` for a full traversal."""
        return [_CommercialSaaSHistoryNode(c)
                for c in self.element.findall("Node")]

    @property
    def regenerated_sites(self) -> "list[dict]":
        """Restriction sites that the cloning operation preserved
        or recreated. Each entry is a plain dict ``{name, pos,
        siteCount}`` so callers don't have to wrap individual
        elements."""
        out: list[dict] = []
        for el in self.element.findall("RegeneratedSite"):
            out.append({
                "name":      el.get("name", ""),
                "pos":       _coerce_int_or_zero(el.get("pos")),
                "siteCount": _coerce_int_or_zero(el.get("siteCount")),
            })
        return out

    @property
    def input_summaries(self) -> "list[dict]":
        out: list[dict] = []
        for el in self.element.findall("InputSummary"):
            out.append({
                "manipulation": el.get("manipulation", ""),
                "name1":        el.get("name1", ""),
                "name2":        el.get("name2", ""),
                "val1":         _coerce_int_or_zero(el.get("val1")),
                "val2":         _coerce_int_or_zero(el.get("val2")),
                "siteCount1":   _coerce_int_or_zero(el.get("siteCount1")),
                "siteCount2":   _coerce_int_or_zero(el.get("siteCount2")),
            })
        return out

    @property
    def oligos(self) -> "list[dict]":
        """PCR primers recorded on an ``amplifyFragment`` history node.
        CommercialSaaS stores them as ``<Oligo name= sequence= …>``
        children — the amplify ``InputSummary`` itself carries no
        name1/name2 (only the val1/val2 amplified-region coordinates),
        which is why an un-parsed amplify step reads as a bare "amplify
        ? ↔ ?". Returns ``{name, sequence, description}`` per oligo so
        the detail pane can name the primers that made the product."""
        out: list[dict] = []
        for el in self.element.findall("Oligo"):
            out.append({
                "name":        el.get("name", "") or "",
                "sequence":    el.get("sequence", "") or "",
                "description": el.get("description", "") or "",
            })
        return out

    @property
    def parameters(self) -> "list[dict]":
        """``<Parameter name= val=>`` children — the free-form knobs a
        manipulation was run with (real files record e.g.
        ``polymerase = "polymerase that creates blunt ends"`` on every
        PCR node). Returns ``{name, value}`` dicts. Purely informational,
        but it IS a recorded property of the step, so the detail pane
        lists it rather than dropping it on the floor."""
        return [{"name":  el.get("name", "") or "",
                 "value": el.get("val", "") or ""}
                for el in self.element.findall("Parameter")]

    @property
    def primer_details(self) -> "list[dict]":
        """The `<Primers>/<Primer>` block — richer than `oligos`.

        `<Oligo>` (on the amplify node) names the primers; the SEPARATE
        `<Primers>` block (on the TEMPLATE node) records where each one
        actually annealed. Per primer:

          * ``name`` / ``sequence`` / ``description``
          * ``binding_sites`` — ``{location, strand, annealed_bases, tm,
            flap}``, where ``flap`` is the 5′ tail that did NOT hybridise
            (reconstructed from the `<Component>` children: a component
            carrying ``hybridizedRange`` is annealed, one without is
            flap). That flap is the added restriction site / overhang,
            i.e. exactly what the cloning step engineered in.

        This is the single largest block of recorded detail the viewer
        used to discard (908 binding sites across the reference
        library)."""
        out: list[dict] = []
        for block in self.element.findall("Primers"):
            for el in block.findall("Primer"):
                sites: list[dict] = []
                for bs in el.findall("BindingSite"):
                    flap_parts: list[str] = []
                    for comp in bs.findall("Component"):
                        if comp.get("hybridizedRange"):
                            continue          # annealed, not flap
                        flap_parts.append(comp.get("bases", "") or "")
                    sites.append({
                        "location":       bs.get("location", "") or "",
                        "strand":         _coerce_int_or_zero(
                            bs.get("boundStrand")),
                        "annealed_bases": bs.get("annealedBases", "") or "",
                        "tm":             bs.get("meltingTemperature", "") or "",
                        "flap":           "".join(flap_parts),
                    })
                out.append({
                    "name":          el.get("name", "") or "",
                    "sequence":      el.get("sequence", "") or "",
                    "description":   el.get("description", "") or "",
                    "binding_sites": sites,
                })
        return out

    @property
    def hybridization_params(self) -> "dict[str, str]":
        """``<Primers>/<HybridizationParams>`` — the annealing stringency
        the binding sites above were computed under (min continuous
        match, mismatch allowance, min Tm, …). Empty dict when absent."""
        for block in self.element.findall("Primers"):
            for el in block.findall("HybridizationParams"):
                return {k: str(v) for k, v in el.attrib.items()}
        return {}

    @property
    def end_modifications(self) -> "dict[str, str]":
        """5′/3′ end chemistry recorded on the node — ``upstreamModification``
        / ``downstreamModification`` attributes (``"Unmodified"``,
        ``"5' phosphorylated"``, …), plus any `<UpstreamModification>` /
        `<DownstreamModification>` child element. Empty dict when the
        node records none. Real, benchwork-relevant state: whether a
        fragment carried a 5′ phosphate decides whether it ligates."""
        out: "dict[str, str]" = {}
        for attr, key in (("upstreamModification", "upstream"),
                          ("downstreamModification", "downstream")):
            v = (self.element.get(attr) or "").strip()
            if v:
                out[key] = v
        for tag, key in (("UpstreamModification", "upstream"),
                         ("DownstreamModification", "downstream")):
            for el in self.element.findall(tag):
                v = (el.get("name") or el.get("val")
                     or (el.text or "")).strip()
                if v:
                    out[key] = v
        return out

    @property
    def sticky_ends(self) -> "dict[str, str]":
        """``upstreamStickiness`` / ``downstreamStickiness`` (+
        ``countBasesFrom``) — the overhang state of a LINEAR fragment
        node. Empty dict when absent (every circular node)."""
        out: "dict[str, str]" = {}
        for attr, key in (("upstreamStickiness", "upstream"),
                          ("downstreamStickiness", "downstream"),
                          ("countBasesFrom", "count_bases_from")):
            v = (self.element.get(attr) or "").strip()
            if v:
                out[key] = v
        return out

    @property
    def custom_map_label(self) -> str:
        """The user's own label for this node, when the source file set
        ``useCustomMapLabel="1"`` + ``customMapLabel="…"``. Returned only
        when actually in use, so the caller can show it ALONGSIDE the
        real name rather than silently substituting it."""
        if _coerce_int_or_zero(self.element.get("useCustomMapLabel")) == 0:
            return ""
        return (self.element.get("customMapLabel") or "").strip()

    @property
    def source_collapsed(self) -> bool:
        """``hideHistory="1"`` — the source application had this branch
        COLLAPSED in its own history view. Purely a display preference
        there, but worth stating: it marks lineage the originating tool
        was hiding, which SpliceCraft shows anyway."""
        return _coerce_int_or_zero(self.element.get("hideHistory")) != 0

    @property
    def feature_snapshot_count(self) -> int:
        """How many `<Feature>` elements this node's `<Features>` snapshot
        carries — the annotation state of that intermediate at the time
        the step ran. The features themselves are cosmetic bulk (~95% of
        a real history's bytes), so the viewer reports the COUNT rather
        than re-rendering a whole annotation track per node."""
        n = 0
        for block in self.element.findall("Features"):
            n += len(block.findall("Feature"))
        return n

    @property
    def extra_attributes(self) -> "dict[str, str]":
        """Every node attribute WITHOUT a typed getter above. The
        no-info-is-lost backstop: a field this codebase has never seen —
        a newer writer's addition, a vendor extension — still reaches the
        detail pane instead of being invisible. Empty for the files we
        fully model."""
        return {k: str(v) for k, v in self.element.attrib.items()
                if k not in _HISTORY_MODELLED_ATTRS}

    @property
    def other_children(self) -> "list[tuple[str, int]]":
        """``[(tag, count), …]`` for child elements with no dedicated
        renderer — the element-level twin of `extra_attributes`. Cosmetic
        colour/annotation blocks are excluded (they're reported via
        `feature_snapshot_count` instead of counted as unknowns)."""
        seen: "dict[str, int]" = {}
        for el in self.element:
            tag = str(el.tag)
            if tag in _HISTORY_MODELLED_CHILDREN:
                continue
            seen[tag] = seen.get(tag, 0) + 1
        return sorted(seen.items())

    # ── Mutation ─────────────────────────────────────────────────

    def add_parent(self, parent: "_CommercialSaaSHistoryNode") -> None:
        """Attach ``parent`` as a child of this node. CommercialSaaS's
        convention: parent fragments hang off the result node so
        the tree reads "result → parents → grandparents …" as you
        descend. Caller is responsible for setting `parent.node_id`
        to a value unique within the tree."""
        self.element.append(parent.element)

    def add_regenerated_site(self, name: str, pos: int,
                              site_count: int = 1) -> None:
        import xml.etree.ElementTree as _ET
        el = _ET.SubElement(self.element, "RegeneratedSite")
        el.set("name", str(name))
        el.set("pos", str(int(pos)))
        el.set("siteCount", str(int(site_count)))

    def add_input_summary(self, *, manipulation: str,
                            name1: str = "", name2: str = "",
                            val1: int = 0, val2: int = 0,
                            site_count1: int = 1,
                            site_count2: int = 1) -> None:
        import xml.etree.ElementTree as _ET
        el = _ET.SubElement(self.element, "InputSummary")
        el.set("manipulation", str(manipulation))
        el.set("name1", str(name1))
        el.set("name2", str(name2))
        el.set("val1", str(int(val1)))
        el.set("val2", str(int(val2)))
        el.set("siteCount1", str(int(site_count1)))
        el.set("siteCount2", str(int(site_count2)))

    def add_oligo(self, *, name: str, sequence: str,
                   description: str = "") -> None:
        """Record a PCR primer on this node as an `<Oligo>` child —
        SnapGene/CommercialSaaS's representation. Lets SpliceCraft-
        generated PCR history carry the SAME primer detail a `.dna`
        import does (`oligos` reads them back; the History detail's
        Primers block renders them) — harmonised history regardless of
        whether the plasmid was imported or built de-novo in SpliceCraft
        (user request 2026-06-01)."""
        import xml.etree.ElementTree as _ET
        el = _ET.SubElement(self.element, "Oligo")
        el.set("name", str(name))
        el.set("sequence", str(sequence))
        if description:
            el.set("description", str(description))

    def add_parameter(self, *, name: str, value: str) -> None:
        """Record a `<Parameter name= val=>` on this node — the element
        an imported `.dna` uses for a step's run conditions. Lets
        SpliceCraft-generated history carry the same detail an import
        does (`parameters` reads it back), closing the last "not yet
        generated" gap noted in `[INV-93]`."""
        import xml.etree.ElementTree as _ET
        el = _ET.SubElement(self.element, "Parameter")
        el.set("name", str(name))
        el.set("val", str(value))

    # ── Traversal ────────────────────────────────────────────────

    def walk(self):
        """Pre-order depth-first traversal: yields self, then each
        child node's traversal. Useful for "find every node with
        operation X" / "count parent fragments" / etc. Returns a
        generator of `_CommercialSaaSHistoryNode`.

        Iterative (stack-based) so a hostile `.dna` file with a 1000+
        deep nested `<Node><Node>...` history can't trip the CPython
        recursion limit. The total node count is still bounded by the
        LZMA-decompression cap (`_COMMERCIALSAAS_HISTORY_MAX_XML`).
        """
        stack: list = [self]
        while stack:
            node = stack.pop()
            yield node
            # Push children in reverse so the first child pops next,
            # preserving the original pre-order traversal sequence.
            stack.extend(reversed(node.parents))


def _coerce_int_or_zero(s) -> int:
    """Best-effort int coercion; falls back to 0. Used by history
    attribute getters where a malformed XML attribute shouldn't
    blow up the whole tree traversal."""
    try:
        return int(s) if s is not None else 0
    except (TypeError, ValueError):
        return 0


# ── History-viewer presentation: tree labels, protocol table, detail lines (Phase D)
_HISTORY_LABEL_NAME_MAX: int = 40


_HISTORY_DETAIL_LIST_MAX: int = 12


_HISTORY_PROTOCOL_INPUT_MAX: int = 8   # parts listed per step before "+N more"


# EVERY operation string real `.dna` files carry, mapped to a bench verb.
# Derived from a full scan of the reference library (930 `.dna` files / 1184
# stored plasmids), not guessed: an unmapped op used to render as raw
# camelCase, so `goldenGateAssembly` — the 4th most common manipulation
# there — read as machine noise in the middle of a protocol.
# An op NOT listed here still passes through verbatim (never dropped).
_HISTORY_OP_FRIENDLY: "dict[str, str]" = {
    # assembly / joining
    "insertFragment":   "assemble",
    "insertFragments":  "assemble",
    "goldenGateAssembly": "Golden Gate",
    "gibsonAssembly":   "Gibson",
    "hifiAssembly":     "HiFi assembly",
    "insert":           "insert",
    "replace":          "replace",
    # amplification
    "amplifyFragment":  "PCR",
    # editing
    "editSequence":     "edit",
    "primerDirectedMutagenesis": "mutagenesis",
    "insertCodons":     "insert codons",
    "insertRestrictionSite": "insert site",
    "remove":           "delete",
    "editDNAEnds":      "edit ends",
    "changeMethylation": "methylation",
    # topology / framing
    "setOrigin":        "set origin",
    "flip":             "reverse-complement",
    "circularize":      "circularize",
    "linearize":        "linearize",
    # provenance roots
    "newFileFromSelection": "extract region",
    "importFile":       "import",
    "createDocument":   "create",
    # NB: `truncated` is deliberately absent. The marker node is NAMED
    # "(earlier history truncated)", so a friendly verb would just repeat
    # the name in the same row; the raw "truncated" passthrough reads as a
    # short tag beside it.
}


_HISTORY_OP_GLYPH: "dict[str, str]" = {
    "assemble":           "⊕",
    "Golden Gate":        "⊕",
    "Gibson":             "⊕",
    "HiFi assembly":      "⊕",
    "insert":             "⊕",
    "insert codons":      "⊕",
    "insert site":        "⊕",
    "replace":            "⇄",
    "delete":             "⊖",
    "extract region":     "⊖",
    "PCR":                "⇉",
    "edit":               "✎",
    "mutagenesis":        "✎",
    "edit ends":          "✎",
    "methylation":        "✎",
    "set origin":         "↻",
    "circularize":        "↻",
    "linearize":          "↧",
    "reverse-complement": "⇋",
    "import":             "↦",
    "create":             "✦",
}


def _history_clean_name(name: str) -> str:
    """Strip the cosmetic ``.dna`` suffix history nodes carry (CommercialSaaS
    convention — see `_build_history_for_product`) for a cleaner row.
    Empty / whitespace name → ``(unnamed)`` so a row never collapses to
    blank. Does NOT escape — callers escape at render time."""
    n = (name or "").strip()
    if n[-4:].lower() == ".dna":
        n = n[:-4].strip()
    return n or "(unnamed)"


def _history_size_label(bp: int) -> str:
    """Compact length for a tree row: ``712 bp`` · ``13.9 kb`` ·
    ``1.20 Mb``. The exact base count still shows in the detail pane —
    the tree just wants something short."""
    try:
        b = int(bp)
    except (TypeError, ValueError):
        b = 0
    if b < 0:
        b = 0
    if b < 1_000:
        return f"{b} bp"
    if b < 1_000_000:
        return f"{b / 1_000:.1f} kb"
    return f"{b / 1_000_000:.2f} Mb"


_HISTORY_OP_SENTINELS: "frozenset[str]" = frozenset(
    {"invalid", "unknown", "none", "unspecified"})


def _history_op_label(op: str) -> str:
    """Friendly verb for an operation string; unknown ops pass through
    verbatim. Empty op — OR a CommercialSaaS sentinel like ``invalid``
    (its placeholder for a base/starting sequence with no recorded
    operation) — returns "" so the caller omits the tag instead of
    printing a literal "invalid"."""
    raw = (op or "").strip()
    if not raw or raw.lower() in _HISTORY_OP_SENTINELS:
        return ""
    return _HISTORY_OP_FRIENDLY.get(raw, raw)


# Manipulations whose coordinate is a single POSITION, not a range. Driven by
# the manipulation name rather than "which value happens to be zero", because a
# genuine range can legitimately start at base 0 — `amplify 0–542` is a real
# region, `insertAt 4707–0` is one insertion point with an unused second slot.
# Rendering the latter as a range produced the nonsense "renumberBase 610–0".
_HISTORY_POINT_MANIPULATIONS: "frozenset[str]" = frozenset(
    {"insertAt", "renumberBase", "cutAt"})

# Phrasing per manipulation. Counts are occurrences across the reference
# library, so the common cases read naturally rather than as raw field names.
_HISTORY_COORD_PHRASE: "dict[str, str]" = {
    "amplify":      "region {a}–{b}",      # 480
    # Named rather than left to the default "{manip} {a}–{b}", which
    # printed the verb twice: "… · replace   replace 0–1,986".
    "replace":      "bases {a}–{b}",       # 1,991
    "remove":       "bases {a}–{b}",       # 208
    "select":       "selection {a}–{b}",   # 172
    "editSequence": "{a} → {b} bp",
    "insertAt":     "at {a}",              # 494
    "renumberBase": "origin → {a}",        # 124
    "cutAt":        "cut at {a}",
}


# A hostile `.dna` can carry thousands of RegeneratedSite / BindingSite
# children on one node. The pane shows a dozen and the header shows a count,
# so past this the extra work buys nothing — and the checker runs inside
# `compose()`, where a slow render is a hung screen.
_HISTORY_WARN_MAX: int = 50

# Bound the WORK, not just the output. The warning cap above stops at 50
# reported problems, but a node whose 900 binding sites all match produces
# zero warnings and still costs 900 scans — so cap the claims inspected too.
_HISTORY_CHECK_MAX_ITEMS: int = 200

# Above this, run only the arithmetic checks. Each scan is O(len(seq)), so a
# 4.5 Mb chromosome with 20 binding sites measured ~460 ms — and this runs
# inside `compose()`, where that is a visible stall. Chromosome-scale records
# aren't what construction-history checking is for; the recorded-length checks
# still run on them.
_HISTORY_CHECK_MAX_BP: int = 2_000_000


def _history_node_warnings(node: "_CommercialSaaSHistoryNode",
                            seq: str = "", *, enzymes=None) -> "list[str]":
    """Claims this node makes that the evidence contradicts. Plain text;
    empty list when everything checks out.

    STRICTLY READ-ONLY. It reports a disagreement between the record and
    the molecule — it never "repairs" either, because both are the user's
    data: the history is a faithful account of what some tool actually
    did, and the sequence is the plasmid. Rewriting either to silence a
    warning would destroy the evidence that they diverged, which is the
    one thing the user needs to see.

    ``seq`` must be THIS node's own molecule (callers pass it only for the
    node they can identify with certainty — in the viewers, the loaded
    record's root). Without it, only the checks that need nothing but the
    recorded lengths run, and those are safe on every node.

    ``enzymes`` is a ``{name: (site, …)}`` catalog; passed in rather than
    imported so this module stays free of the data layer.

    Every check here was validated against a 1,184-plasmid reference
    library, and only checks that held for ~99% of real records are
    included — a warning that fires on correct data is worse than no
    warning at all.
    """
    out: "list[str]" = []
    try:
        return _history_node_warnings_impl(node, seq, enzymes, out)
    except Exception:
        # This runs inside `compose()` on both viewers. A malformed value
        # anywhere in a hostile `.dna` must degrade to "no warnings", never
        # to a History screen that won't render — the user still needs to
        # SEE their lineage. Partial results are kept: a warning already
        # collected is still true.
        import logging
        logging.getLogger("splicecraft").debug(
            "history consistency check failed for %r; showing %d partial "
            "warning(s)", getattr(node, "name", "?"), len(out), exc_info=True)
        return out


def _history_node_warnings_impl(node, seq, enzymes, out) -> "list[str]":
    """Body of `_history_node_warnings`; see there for the contract. Split
    out so the wrapper can keep partial results on an unexpected value."""
    from splicecraft_biology import _iupac_pattern, _rc
    op = (node.operation or "").strip()
    parents = node.parents
    nlen = int(node.seq_len or 0)

    # ── recorded-length arithmetic (no sequence required) ────────────
    if parents:
        plen = int(parents[0].seq_len or 0)
        if op in ("setOrigin", "flip", "circularize", "linearize",
                  "changeMethylation") and plen and nlen and plen != nlen:
            out.append(
                f"{_history_op_label(op) or op} should not change length, "
                f"but the parent is {plen:,} bp and this is {nlen:,} bp")
        if op == "remove" and plen and nlen:
            for sm in node.input_summaries:
                if str(sm.get("manipulation") or "") != "remove":
                    continue
                # 0-based INCLUSIVE span — calibrated against the reference
                # library, where 193/193 removals balance on this rule.
                span = (_coerce_int_or_zero(sm.get("val2"))
                        - _coerce_int_or_zero(sm.get("val1")) + 1)
                if span > 0 and plen - span != nlen:
                    out.append(
                        f"deletion of {span:,} bp from {plen:,} bp should "
                        f"leave {plen - span:,} bp, but this is {nlen:,} bp")
                break

    n = len(seq)
    if (not seq or len(out) >= _HISTORY_WARN_MAX
            or n > _HISTORY_CHECK_MAX_BP):
        return out

    # ── claims checkable against the actual DNA ──────────────────────
    circ = bool(node.circular)
    seen = 0                       # claims inspected, capped independently
    # Built ONCE, and LAZILY — doubling the sequence inside the
    # per-binding-site loop turned a detail-pane render into hundreds of MB
    # of copying, and building it eagerly wasted a full copy on every node
    # that records no binding sites at all (most of them).
    _hay: "list[str]" = []

    def hay_of() -> str:
        if not _hay:
            _hay.append((seq + seq) if circ else seq)
        return _hay[0]

    for s in (node.regenerated_sites or []):
        if len(out) >= _HISTORY_WARN_MAX or seen >= _HISTORY_CHECK_MAX_ITEMS:
            break
        seen += 1
        nm = str(s.get("name") or "").strip()[:_HISTORY_MANIP_LABEL_MAX]
        pos = _coerce_int_or_zero(s.get("pos"))
        # pos == 0 is SpliceCraft's OWN "this is the enzyme the assembly
        # used" marker (`_build_history_for_assembly`, [INV-93]) — not a
        # claim that a site survives in the product. Type IIS assembly
        # consumes its sites, so checking those would warn on 126 correct
        # records in the reference library.
        if not nm or nm.startswith("(") or pos == 0:
            continue
        spec = (enzymes or {}).get(nm)
        site = str(spec[0]).upper() if spec else ""
        # An empty recognition site would compile to a pattern that matches
        # at every offset — "always present", i.e. the check silently
        # stops checking. Skip instead.
        if not site or len(site) > n:
            continue
        # Circular: extend by site_len-1 so a site straddling the origin is
        # found once (sacred invariant #6's wrap-scan rule), not twice.
        window = (seq + seq[:len(site) - 1]) if circ else seq
        if (_iupac_pattern(site).search(window)
                or _iupac_pattern(_rc(site)).search(window)):
            continue
        out.append(
            f"records a regenerated {nm} site at {pos:,}, but {nm} "
            f"({site}) does not occur anywhere in this molecule")

    for pr in (node.primer_details or []):
        if len(out) >= _HISTORY_WARN_MAX or seen >= _HISTORY_CHECK_MAX_ITEMS:
            break
        pname = str(pr.get("name") or "primer")[:_HISTORY_LABEL_NAME_MAX]
        for bs in (pr.get("binding_sites") or []):
            if (len(out) >= _HISTORY_WARN_MAX
                    or seen >= _HISTORY_CHECK_MAX_ITEMS):
                break
            seen += 1
            ab = str(bs.get("annealed_bases") or "").upper()
            loc = str(bs.get("location") or "")
            # A recorded "annealed" run longer than the molecule can't be
            # a real binding site; scanning for it is pure cost.
            if not ab or "-" not in loc or len(ab) > n:
                continue
            try:
                start = int(loc.split("-", 1)[0])
            except ValueError:
                continue
            # A negative offset would index from the END of the string and
            # quietly compare the WRONG bases — a false pass or a false
            # warning, both worse than skipping. Out-of-range is just a
            # miss, which the "anneals elsewhere" branch reports correctly.
            if start < 0 or start > 2 * n:
                continue
            want = ab if not _coerce_int_or_zero(bs.get("strand")) else _rc(ab)
            hay = hay_of()
            if hay[start:start + len(want)] == want:
                continue
            if want in hay:
                out.append(
                    f"{pname} is recorded at {start:,} but anneals elsewhere "
                    f"— the molecule was re-origined after this was written, "
                    f"so the position is stale (the primer still binds)")
            else:
                out.append(
                    f"{pname} does not bind this molecule at all — the "
                    f"sequence has changed since this step was recorded")
    return out


def _history_consistency_summary(warnings: "list[str]") -> str:
    """One-line markup for a viewer header. "" when there's nothing to
    flag, so a clean history shows no chrome at all."""
    n = len(warnings or [])
    if not n:
        return ""
    return (f"[yellow]⚠ {n} recorded detail"
            f"{'' if n == 1 else 's'} the sequence doesn't support[/yellow]")


_HISTORY_MANIP_LABEL_MAX: int = 32


def _history_coord_phrase(manip, val1, val2) -> str:
    """Render one summary's coordinates — ``region 469–969``,
    ``bases 450–458``, ``origin → 610``. Returns "" when there's nothing
    to say.

    THE single formatter: both the tree/protocol tag
    (`_history_manipulation_detail`) and the detail pane's Inputs block go
    through it, so a point manipulation can't render as a range in one
    place and a position in the other.

    ``manip`` is deliberately untyped: every caller reads it straight off
    an XML attribute dict, so it arrives as ``str | None`` (or anything a
    hostile file put there). It is normalised on the first line.

    SACRED: ``.format()`` is applied ONLY to the templates in
    `_HISTORY_COORD_PHRASE`, never to one built from ``manip``. The
    manipulation name comes straight off attacker-controllable `.dna` XML,
    and a value like ``{evil}`` folded into a format string raises
    ``KeyError`` — which would take out the tree label AND the protocol
    pane for that plasmid. The unknown-manipulation case is concatenated,
    and the label is length-capped so a 10 KB attribute can't push the
    row off-screen."""
    m = str(manip or "").strip()[:_HISTORY_MANIP_LABEL_MAX]
    v1 = _coerce_int_or_zero(val1)
    v2 = _coerce_int_or_zero(val2)
    if not (v1 or v2):
        return ""
    tmpl = _HISTORY_COORD_PHRASE.get(m)          # trusted templates only
    if m in _HISTORY_POINT_MANIPULATIONS:
        pos = v1 or v2
        if tmpl:
            return tmpl.format(a=f"{pos:,}", b=f"{pos:,}")
        return f"{m} {pos:,}"
    if tmpl:
        return tmpl.format(a=f"{v1:,}", b=f"{v2:,}")
    return f"{m or 'at'} {v1:,}–{v2:,}"


def _history_manipulation_detail(node: "_CommercialSaaSHistoryNode") -> str:
    """WHERE on the molecule a step acted, distilled from its
    `<InputSummary>` coordinates — ``region 469–969`` for a PCR,
    ``bases 450–458`` for a deletion, ``origin → 610`` for a re-origin.

    Those val1/val2 numbers were parsed but never rendered anywhere: the
    protocol showed WHAT combined, never WHERE it landed. Returns "" when
    the summaries carry no usable coordinates (an assembly records
    enzyme names in name1/name2 instead, and the ``✂`` tag covers that).
    Plain text — the caller escapes.

    PERF: iterates the `<InputSummary>` ELEMENTS rather than the
    `input_summaries` property, which materialises a dict per summary.
    `_history_tree_label` calls this for EVERY node, and the tree populate
    walks up to `_HISTORY_NODE_MAX_NODES` (100,000) of them — so a dict
    build here is paid 100,000 times to produce one short string. Measured
    at 5.7 µs of a 19.9 µs label before this change."""
    for el in node.element.findall("InputSummary"):
        # An assembly summary's val1/val2 are cut positions on two
        # DIFFERENT molecules — a range across them is meaningless.
        if (el.get("name1") or "") or (el.get("name2") or ""):
            continue
        phrase = _history_coord_phrase(
            el.get("manipulation"), el.get("val1"), el.get("val2"))
        if phrase:
            return phrase
    return ""


def _history_node_signature(node: "_CommercialSaaSHistoryNode") -> str:
    """Identity key for collapsing repeated ancestor subtrees in the
    viewer. Two nodes with the same cleaned name + length + operation
    are treated as the same plasmid — repeats come from the same parent
    entry's inherited ``history_xml``, so the subtree IS identical.

    Purely cosmetic: a false match merely renders the 2nd occurrence as
    a ``↳ … (shown above)`` reference instead of redrawing it. It never
    touches stored history."""
    return (f"{_history_clean_name(node.name)}\x00{int(node.seq_len)}"
            f"\x00{(node.operation or '').strip()}")


def _history_tree_label(node: "_CommercialSaaSHistoryNode") -> str:
    """De-noised one-line tree label: ``name   size   ⊕ op`` (markup
    string; render via ``Text.from_markup``).

    The old format crammed name · N bp · circular · operation onto
    every row; the bp count + topology + raw op string repeated on all
    of them. Now the row leads with the (``.dna``-stripped) name + a
    compact size, drops "circular" (the common case — only "linear" is
    flagged), and shows a friendly op verb. Full bp / topology /
    strandedness move to the detail pane (`_history_detail_lines`).

    Names + ops are Rich-escaped (XML attrs can legally contain ``[``)
    and truncated so a hostile value can't push the column off-screen."""
    from rich.markup import escape as _esc
    raw = _history_clean_name(node.name)
    if len(raw) > _HISTORY_LABEL_NAME_MAX:
        raw = raw[: _HISTORY_LABEL_NAME_MAX - 1] + "…"
    name = _esc(raw)
    bits = [f"[b]{name}[/b]", f"[dim]{_history_size_label(node.seq_len)}[/dim]"]
    if not node.circular:
        bits.append("[dim]linear[/dim]")
    # Every node that records a REAL operation is tagged — including a
    # leaf. The old rule ("only nodes with parents") existed to stop raw
    # starting material reading "⊕ assemble", but `_history_op_label`
    # already blanks CommercialSaaS's `invalid` sentinel, which is what
    # those leaves actually carry. Requiring parents on top of that hid
    # genuine root actions — an `importFile` / `createDocument` /
    # `newFileFromSelection` leaf showed as a bare name with no hint of
    # how the molecule entered the project.
    op = _history_op_label(node.operation)
    if op:
        glyph = _HISTORY_OP_GLYPH.get(op, "")
        tag = f"{glyph} {_esc(op)}".strip()
        bits.append(f"[cyan]{tag}[/cyan]")
    # WHERE the step acted (PCR region / deleted span / edited length),
    # straight from the InputSummary coordinates.
    where = _history_manipulation_detail(node)
    if where:
        bits.append(f"[dim]{_esc(where)}[/dim]")
    # Date+time of the step, in the universal slash-free "JUN 9 2026 14:30"
    # form, right alongside the action — undated reconstructed-lineage rows
    # just omit it (no bogus date).
    when = _history_human_dt(node.date)
    if when:
        bits.append(f"[dim]· {when}[/dim]")
    return "   ".join(bits)


def _history_reference_label(node: "_CommercialSaaSHistoryNode") -> str:
    """Compact label for the 2nd+ occurrence of a repeated ancestor —
    ``↳ name  size  (shown above)``. No operation tag (the canonical
    occurrence carries the detail); the marker tells the user the full
    subtree lives elsewhere in the tree rather than being drawn again."""
    from rich.markup import escape as _esc
    raw = _history_clean_name(node.name)
    if len(raw) > _HISTORY_LABEL_NAME_MAX:
        raw = raw[: _HISTORY_LABEL_NAME_MAX - 1] + "…"
    name = _esc(raw)
    return (f"[dim]↳[/dim] {name}   [dim]{_history_size_label(node.seq_len)}"
            f"[/dim]   [dim i](shown above)[/dim i]")


def _history_name_key(name: str) -> str:
    """Comparison key for "is this the same plasmid name?" — case-folded
    with ``_`` treated as a space and runs of whitespace collapsed.

    The underscore rule matters: a caller may hand us the LOCUS-safe id
    (``FFE_102_TU_35sPPDK``) where the history root holds the display name
    (``FFE 102 TU 35sPPDK``). Those are the SAME plasmid — reporting them
    as a rename would put a "built as …" banner on a plasmid that was
    never renamed ([INV-98] is about keeping those two forms distinct in
    STORAGE; here we must see through the difference)."""
    return " ".join(str(name or "").replace("_", " ").split()).casefold()


def _history_former_name(root: "_CommercialSaaSHistoryNode",
                          entry_name: str) -> str:
    """The name this molecule was BUILT UNDER, when it differs from what
    the library calls it today — else "".

    History roots keep the filename the source recorded at build time, so
    a plasmid renamed afterwards opens its History on a root reading
    ``AB000000`` for an entry the user knows as ``pDemo42``. That's
    accurate provenance, but silently confusing: 177 of the reference
    library's 1,184 plasmids are in exactly this state. Surfacing it as
    "built as …" turns a puzzling mismatch into a useful fact."""
    was = _history_clean_name(root.name).strip()
    now = str(entry_name or "").strip()
    if not was or not now or was == "(unnamed)":
        return ""
    if _history_name_key(was) == _history_name_key(now):
        return ""
    # Capped like a tree label. This lands in a header Static above the
    # tree; uncapped, a 20,000-character name (hostile file, or just a very
    # long one) wraps across the pane and pushes the lineage off-screen.
    if len(was) > _HISTORY_LABEL_NAME_MAX:
        was = was[:_HISTORY_LABEL_NAME_MAX - 1] + "…"
    return was


def _history_renumber_node_ids(root: "_CommercialSaaSHistoryNode") -> int:
    """Give every node in a tree a unique ``ID``, in pre-order from 0.
    Returns the count assigned.

    Nesting an imported parent's history under a new root merges two
    independently-numbered trees, so IDs collide (found on 10 plasmids in
    the reference library, e.g. IDs ``[3, 1, 2, 10, 4, 11]`` all
    duplicated in one tree). Nothing keys off the ID today, which is why
    it went unnoticed — but a duplicate ID in a lineage record is wrong
    on its face, and any consumer that indexes by it silently merges
    unrelated steps. Iterative; mutates in place."""
    n = 0
    stack: list = [root.element]
    while stack:
        el = stack.pop()
        el.set("ID", str(n))
        n += 1
        # Reverse so children are numbered in document order.
        stack.extend(reversed(el.findall("Node")))
    return n


def _history_populate_tree(tree, root: "_CommercialSaaSHistoryNode",
                            node_by_id: dict) -> bool:
    """Fill a Textual ``Tree`` from a history root, de-noised:

      * **Progressive disclosure** — only the product (the root) opens
        by default; deeper generations start collapsed so a multi-step
        build doesn't dump its whole fractal on open.
      * **Dedup repeated ancestors** — a subtree seen before renders as
        a single ``↳ … (shown above)`` leaf, so a backbone reused
        across N branches is drawn once, not N times (the dominant
        noise on real GB/MoClo lineages).

    Iterative DFS with the shared depth + node caps (mirrors
    `_history_node_to_dict` / `_CommercialSaaSHistoryNode.walk` — a
    hostile deeply-nested ``.dna`` can't trip the recursion limit).
    Populates ``node_by_id`` ``{textual_node_id: hist_node}`` for the
    detail pane. Returns ``True`` if a cap truncated the render."""
    seen: "set[str]" = set()
    truncated = False
    n_seen = 0
    # Frame: (textual_parent_node, hist_node, depth). Pre-order; push
    # children reversed so the first parent pops next.
    stack: list = [(tree.root, root, 0)]
    while stack:
        parent_tnode, hist, depth = stack.pop()
        if n_seen >= _HISTORY_NODE_MAX_NODES or depth >= _HISTORY_NODE_MAX_DEPTH:
            truncated = True
            continue
        sig = _history_node_signature(hist)
        is_ref = sig in seen
        if is_ref:
            label = _history_reference_label(hist)
        else:
            seen.add(sig)
            label = _history_tree_label(hist)
        child = parent_tnode.add(Text.from_markup(label), expand=(depth == 0))
        node_by_id[child.id] = hist
        n_seen += 1
        if is_ref:
            # Reference occurrence — don't redraw the (identical) subtree.
            continue
        for p in reversed(hist.parents):
            stack.append((child, p, depth + 1))
    return truncated


def _history_step_from_node(node: "_CommercialSaaSHistoryNode") -> dict:
    """Distil one build STEP from a history node that has parents.
    ``backbone`` is the entry-vector / acceptor; ``inputs`` the parts
    that went into it; ``enzyme`` the Type IIS / RE that joined them.

    Backbone detection: our builders add the acceptor first AND record
    it as ``InputSummary.name1`` — prefer that match, else fall back to
    the first parent. Imported CommercialSaaS history that doesn't follow
    the convention degrades to "first parent is the backbone", which is
    cosmetic-only (the protocol line may swap which input is labelled
    the acceptor)."""
    parents = node.parents
    sites = node.regenerated_sites
    # Collect EVERY enzyme — a restriction digest uses ≥2, and the
    # builders record one RegeneratedSite per enzyme. Pre-fix this took
    # only sites[0], so a KpnI + XbaI double digest showed just "✂ KpnI"
    # in the protocol (user report 2026-06-01). Dedup, order-preserve,
    # and drop the "(reverse insert)" orientation sentinel.
    enzymes: "list[str]" = []
    _seen_enz: set = set()
    for _s in sites:
        _nm = str(_s.get("name") or "").strip()
        if not _nm or _nm.startswith("(") or _nm in _seen_enz:
            continue
        _seen_enz.add(_nm)
        enzymes.append(_nm)
    enzyme = enzymes[0] if enzymes else ""   # back-compat single-enzyme field
    summaries = node.input_summaries
    backbone_label = (_history_clean_name(summaries[0].get("name1") or "")
                      if summaries else "")
    backbone = ""
    inputs: "list[str]" = []
    if parents:
        names = [_history_clean_name(p.name) for p in parents]
        if len(names) == 1:
            inputs = names
        else:
            bi = 0
            if backbone_label:
                for i, nm in enumerate(names):
                    if nm == backbone_label:
                        bi = i
                        break
            backbone = names[bi]
            inputs = [nm for i, nm in enumerate(names) if i != bi]
    return {
        "product":  _history_clean_name(node.name),
        "op":       _history_op_label(node.operation),
        "enzyme":   enzyme,        # first enzyme (back-compat)
        "enzymes":  enzymes,       # ALL enzymes in the digest
        "backbone": backbone,
        "inputs":   inputs,
        "seq_len":  int(node.seq_len),
        "circular": bool(node.circular),
        # WHERE the manipulation acted (PCR region / deleted span / …).
        "where":    _history_manipulation_detail(node),
        # Run conditions the source file recorded for the step (e.g.
        # "polymerase = polymerase that creates blunt ends").
        "parameters": node.parameters,
    }


def _history_build_steps(root: "_CommercialSaaSHistoryNode") -> "list[dict]":
    """Flatten a history tree into an ordered, de-duplicated list of
    build STEPS for the protocol view. A step is any node with parents
    (an assembly / cloning op); raw leaf inputs are not steps. A reused
    sub-assembly collapses to ONE step. Ordered earliest-first (deepest
    in the tree) → final product last, so it reads like a bench recipe.

    The protocol naturally sidesteps the combinatorial leaf-duplication
    that bloats the tree: duplication lives in the leaf ancestors, and
    leaves are never steps. Iterative + capped (hostile-XML safe)."""
    by_sig: "dict[str, dict]" = {}
    depth_of: "dict[str, int]" = {}
    order: "list[str]" = []
    stack: list = [(root, 0)]
    n_seen = 0
    while stack:
        node, depth = stack.pop()
        if n_seen >= _HISTORY_NODE_MAX_NODES or depth >= _HISTORY_NODE_MAX_DEPTH:
            break
        n_seen += 1
        parents = node.parents
        if parents:
            sig = _history_node_signature(node)
            # A reused sub-assembly sorts to its EARLIEST (deepest) use.
            depth_of[sig] = max(depth_of.get(sig, 0), depth)
            if sig not in by_sig:
                by_sig[sig] = _history_step_from_node(node)
                order.append(sig)
        # Push reversed so siblings pop in natural tree order — keeps
        # the protocol reading vector-then-parts, first-listed-first.
        for p in reversed(parents):
            stack.append((p, depth + 1))
    # Deepest first (earliest build) → product (depth 0) last. Stable on
    # insertion order for ties (siblings keep tree order).
    order.sort(key=lambda s: depth_of.get(s, 0), reverse=True)
    return [by_sig[s] for s in order]


def _history_protocol_step_cells(
        root: "_CommercialSaaSHistoryNode") -> "list[tuple[str, str]]":
    """``[(number, content_markup), …]`` — one pair per de-duplicated
    build step, earliest first. ``number`` is like ``"2."``; ``content``
    is the recipe body read left → right:

        <op>  <ingredients> into <backbone>  →  <product>   ✂ <enzymes>

    The forward arrow points AT the product (ingredients → result), the
    acceptor/vector is tagged ``into``, ``✂`` marks the enzymes, and an
    in-place edit (single input == product) collapses to one labelled
    product. A bare record with no steps returns a single
    ``("", <placeholder>)`` pair. Splitting the number from the content
    lets `_history_protocol_renderable`'s table hang-indent WRAPPED lines
    under the content rather than under the next step's number (user
    nitpick 2026-06-01)."""
    from rich.markup import escape as _esc
    steps = _history_build_steps(root)
    if not steps:
        # No step nodes — but the root still records HOW the molecule
        # entered the project (imported / created / extracted). Name it
        # rather than the old flat "no construction steps".
        root_op = _history_op_label(root.operation)
        if root_op:
            return [("", f"[dim]Single record — "
                          f"[/dim][cyan]{_esc(root_op)}[/cyan]"
                          f"[dim], no further steps recorded.[/dim]")]
        return [("",
                 "[dim]Single record — no construction steps recorded.[/dim]")]
    cells: "list[tuple[str, str]]" = []
    for i, s in enumerate(steps, 1):
        product = _esc(s["product"])
        ins = s["inputs"]
        backbone = s["backbone"]
        verb = _esc(s["op"]) if s.get("op") else ""
        enzymes = s.get("enzymes") or ([s["enzyme"]] if s.get("enzyme") else [])
        enz_tag = (
            f"   [magenta]✂ {' + '.join(_esc(e) for e in enzymes)}[/magenta]"
            if enzymes else ""
        )
        # Coordinates of the manipulation, when it recorded any — the
        # difference between "PCR template → product" and "PCR template →
        # product  region 469–969".
        where = s.get("where") or ""
        where_tag = f"   [dim]{_esc(where)}[/dim]" if where else ""
        # In-place edit: a single input that IS the product reads
        # redundantly as "X → X" — show one labelled product instead,
        # tagged with what was actually DONE to it. A flat "· edited" here
        # collapsed set-origin, delete, insert and replace into one word,
        # so a plasmid's last four steps all read identically.
        if len(ins) == 1 and not backbone and ins[0] == s["product"]:
            cells.append(
                (f"{i}.",
                 f"[b]{product}[/b]  [dim]· {verb or 'edited'}[/dim]"
                 f"{where_tag}{enz_tag}"))
            continue
        # Ingredients = inputs (+ "into <backbone>" for the acceptor).
        chunks: "list[str]" = []
        if ins:
            shown = ins[:_HISTORY_PROTOCOL_INPUT_MAX]
            joined = " + ".join(_esc(x) for x in shown)
            extra = len(ins) - len(shown)
            if extra > 0:
                joined += f" [dim]+{extra} more[/dim]"
            chunks.append(joined)
        if backbone:
            lead = "[dim]into[/dim] " if chunks else ""
            chunks.append(f"{lead}[cyan]{_esc(backbone)}[/cyan]")
        prefix = f"[dim]{verb}[/dim]  " if verb else ""
        if chunks:
            content = (f"{prefix}{' '.join(chunks)}  [b]→[/b]  "
                       f"[b]{product}[/b]{where_tag}{enz_tag}")
        else:
            content = f"{prefix}[b]{product}[/b]{where_tag}{enz_tag}"
        cells.append((f"{i}.", content))
    return cells


_HISTORY_PROTOCOL_LEGEND = (
    "[dim]ingredients  [b]→[/b]  product"
    "       [cyan]into[/cyan] = acceptor / backbone"
    "       [magenta]✂[/magenta] = enzymes[/dim]"
)


def _history_protocol_renderable(root: "_CommercialSaaSHistoryNode"):
    """Rich renderable for a viewer's protocol pane: the symbol legend
    above a borderless 2-column table (right-justified step-number gutter
    | recipe content). The content column wraps WITHIN its own width, so
    a long step's wrapped tail hangs-indents under the content instead of
    dropping back under the next step's number (user nitpick 2026-06-01).
    Falls back to a bare placeholder Text when there are no steps."""
    from rich.table import Table
    from rich.text import Text
    from rich.console import Group
    cells = _history_protocol_step_cells(root)
    if len(cells) == 1 and cells[0][0] == "":
        return Text.from_markup(cells[0][1])          # placeholder, no legend
    gutter = max((len(num) for num, _ in cells), default=2)
    table = Table(show_header=False, box=None, pad_edge=False,
                   padding=(0, 1, 0, 0), expand=True)
    table.add_column(justify="right", no_wrap=True, width=gutter,
                      style="bold")
    table.add_column(ratio=1, overflow="fold")
    for num, content in cells:
        table.add_row(num, Text.from_markup(content))
    return Group(Text.from_markup(_HISTORY_PROTOCOL_LEGEND), Text(""), table)


def _history_detail_lines(hist: "_CommercialSaaSHistoryNode",
                           seq: str = "", *, enzymes=None) -> "list[str]":
    """Full detail block for the selected history node — shared by both
    viewers so the modal and the full screen never drift. Every dynamic
    field is Rich-escaped and long lists are capped.

    ``seq`` (this node's own molecule) and ``enzymes`` enable the
    consistency block: claims the record makes that the DNA contradicts.
    Omit them and the block is skipped — callers pass a sequence only for
    a node they can identify with certainty."""
    from rich.markup import escape as _esc
    warnings = _history_node_warnings(hist, seq, enzymes=enzymes)
    name_disp = (f"[b]{_esc(_history_clean_name(hist.name))}[/b]"
                 if hist.name else "[dim](unnamed)[/]")
    op_raw = hist.operation
    op_friendly = _history_op_label(op_raw)
    if op_friendly:
        # Don't echo the raw op when the friendly label IS the raw op
        # (e.g. "createDocument") — that just prints "X (X)".
        if op_friendly == op_raw:
            op_disp = f"[cyan]{_esc(op_friendly)}[/]"
        else:
            op_disp = (f"[cyan]{_esc(op_friendly)}[/]  "
                       f"[dim]({_esc(op_raw)})[/]")
    else:
        # Empty, or a CommercialSaaS sentinel ("invalid" / "unknown")
        # for a base / starting sequence with no recorded operation —
        # don't echo the literal sentinel.
        op_disp = "[dim](no operation recorded)[/]"
    strandedness = hist.element.get("strandedness") or "?"
    lines: "list[str]" = [name_disp]
    # The source file's own map label for this node, when it set one.
    # Shown ALONGSIDE the real name, never instead of it.
    custom = hist.custom_map_label
    if custom:
        lines.append(f"[dim]map label:[/] {_esc(custom)}")
    # Straight under the name, because it changes how everything below it
    # should be read. Stated as a disagreement, never as a repair: the
    # record and the molecule are both the user's data and neither is
    # touched to make this go away.
    if warnings:
        lines.append("")
        lines.append("[yellow]⚠ Not supported by the sequence[/yellow]")
        for w in warnings[:_HISTORY_DETAIL_LIST_MAX]:
            lines.append(f"  [yellow]{_esc(w)}[/yellow]")
        if len(warnings) > _HISTORY_DETAIL_LIST_MAX:
            lines.append(
                f"  [dim]… (+{len(warnings) - _HISTORY_DETAIL_LIST_MAX} more)[/]")
    lines.append("")
    lines.append("[b]Properties[/]")
    lines.append(f"  Length:        {hist.seq_len:,} bp")
    lines.append("  Topology:      "
                 f"{'circular' if hist.circular else 'linear'}")
    lines.append(f"  Strandedness:  {_esc(strandedness)}")
    lines.append(f"  Node ID:       {hist.node_id}")
    # End chemistry + overhangs — whether a fragment carried a 5′
    # phosphate, and what its sticky ends were. Recorded on ~2,450 nodes
    # in the reference library and never surfaced before.
    ends = hist.end_modifications
    for key, label in (("upstream", "5′ end:"), ("downstream", "3′ end:")):
        if ends.get(key):
            lines.append(f"  {label:<14s} {_esc(ends[key])}")
    sticky = hist.sticky_ends
    for key, label in (("upstream", "5′ overhang:"),
                       ("downstream", "3′ overhang:"),
                       ("count_bases_from", "Numbering:")):
        if sticky.get(key):
            lines.append(f"  {label:<14s} {_esc(sticky[key])}")
    lines.append("")
    lines.append("[b]Operation[/]")
    lines.append(f"  {op_disp}")
    if hist.resurrectable:
        lines.append("  [green]✓ resurrectable[/] "
                     "[dim](parent can be re-extracted)[/]")
    if hist.source_collapsed:
        lines.append("  [dim]· collapsed in the source file's own "
                     "history view[/]")
    # NB: the manipulation coordinates are deliberately NOT repeated here.
    # `_history_manipulation_detail` only fires for a summary with no
    # name1/name2, and the Inputs block below renders exactly those as
    # "<manip>  (region v1–v2)". The tree label and the protocol — where
    # the coordinates were genuinely missing — carry them instead.
    when = _history_human_dt(hist.date)
    if when:
        lines.append("")
        lines.append("[b]Date[/]")
        lines.append(f"  {when}")
    # Run conditions the source recorded for this step — polymerase
    # choice and friends. Present on every real PCR node; previously
    # parsed by nothing and shown nowhere.
    params = hist.parameters
    if params:
        lines.append("")
        lines.append("[b]Conditions[/]")
        for p in params[:_HISTORY_DETAIL_LIST_MAX]:
            nm = _esc(str(p.get("name") or "(unnamed)"))
            val = _esc(str(p.get("value") or ""))
            lines.append(f"  {nm}: {val}" if val else f"  {nm}")
        if len(params) > _HISTORY_DETAIL_LIST_MAX:
            lines.append(
                f"  [dim]… (+{len(params) - _HISTORY_DETAIL_LIST_MAX} more)[/]")
    sites = hist.regenerated_sites
    if sites:
        lines.append("")
        lines.append("[b]Regenerated sites[/]")
        shown = sites[:_HISTORY_DETAIL_LIST_MAX]
        joined = ", ".join(
            f"{_esc(str(s['name']) or '(unnamed)')}@{s['pos']}"
            for s in shown
        )
        if len(sites) > _HISTORY_DETAIL_LIST_MAX:
            joined += f", … (+{len(sites) - _HISTORY_DETAIL_LIST_MAX} more)"
        lines.append(f"  {joined}")
    sums = hist.input_summaries
    if sums:
        lines.append("")
        lines.append("[b]Inputs[/]")
        for sm in sums[:_HISTORY_DETAIL_LIST_MAX]:
            manip = _esc(str(sm.get('manipulation') or "(unknown)"))
            n1 = str(sm.get('name1') or "")
            n2 = str(sm.get('name2') or "")
            if n1 or n2:
                lines.append(
                    f"  {manip}  ({_esc(n1 or '?')} ↔ {_esc(n2 or '?')})")
            else:
                # No name pair — e.g. an `amplify` (PCR) step, which
                # records its detail as val1/val2 (the amplified region)
                # + <Oligo> primers (the Primers block below), not
                # name1/name2. Don't render a bare "? ↔ ?".
                #
                # Routed through the shared `_history_coord_phrase` — this
                # block used to hardcode "(region v1–v2)", which turned a
                # single-position manipulation into the nonsense
                # "renumberBase (region 6,491–0)".
                phrase = _history_coord_phrase(
                    sm.get('manipulation'), sm.get('val1'), sm.get('val2'))
                suffix = f"  [dim]({_esc(phrase)})[/dim]" if phrase else ""
                lines.append(f"  {manip}{suffix}")
        if len(sums) > _HISTORY_DETAIL_LIST_MAX:
            lines.append(
                f"  [dim]… (+{len(sums) - _HISTORY_DETAIL_LIST_MAX} more)[/]")
    # Primers — CommercialSaaS records PCR oligos as <Oligo> children on
    # an amplify node; surface name + sequence so a PCR step shows WHICH
    # primers made it (the detail was in the .dna all along) instead of a
    # bare "amplify".
    oligos = hist.oligos
    if oligos:
        lines.append("")
        lines.append("[b]Primers[/]")
        for o in oligos[:_HISTORY_DETAIL_LIST_MAX]:
            nm = _esc(o.get("name") or "(unnamed)")
            seq = o.get("sequence") or ""
            seq_disp = (_esc(seq[:40]) + ("…" if len(seq) > 40 else "")
                        if seq else "")
            lines.append(
                f"  {nm}" + (f"   [dim]{seq_disp}[/dim]" if seq_disp else ""))
        if len(oligos) > _HISTORY_DETAIL_LIST_MAX:
            lines.append(
                f"  [dim]… (+{len(oligos) - _HISTORY_DETAIL_LIST_MAX} more)[/]")
    # Where those primers ANNEALED on this molecule — the `<Primers>`
    # block, which the source records on the TEMPLATE node (the `<Oligo>`
    # children above only name them, on the amplify node). Per site:
    # position, strand, Tm, and the 5′ flap that did NOT anneal — the
    # engineered tail (restriction site + overhang) the cloning step
    # added. The single largest block of recorded detail the viewer used
    # to discard.
    primers = hist.primer_details
    if primers:
        lines.append("")
        lines.append("[b]Primer binding[/]")
        hp = hist.hybridization_params
        if hp:
            tm_min = hp.get("minMeltingTemperature", "")
            match_min = hp.get("minContinuousMatchLen", "")
            crit = []
            if match_min:
                crit.append(f"≥{_esc(match_min)} bp match")
            if tm_min:
                crit.append(f"Tm ≥ {_esc(tm_min)}°C")
            if _coerce_int_or_zero(hp.get("allowMismatch")):
                crit.append("mismatches allowed")
            if crit:
                lines.append(f"  [dim]criteria: {', '.join(crit)}[/dim]")
        for p in primers[:_HISTORY_DETAIL_LIST_MAX]:
            nm = _esc(str(p.get("name") or "(unnamed)"))
            lines.append(f"  {nm}")
            for bs in (p.get("binding_sites") or [])[:4]:
                loc = _esc(str(bs.get("location") or "?"))
                arrow = "←" if _coerce_int_or_zero(bs.get("strand")) else "→"
                tm = str(bs.get("tm") or "")
                tm_disp = f"  Tm {_esc(tm)}°C" if tm else ""
                lines.append(f"    {arrow} {loc}{tm_disp}")
                flap = str(bs.get("flap") or "")
                if flap:
                    lines.append(
                        f"      [dim]5′ tail[/dim] {_esc(flap[:48])}"
                        + ("…" if len(flap) > 48 else ""))
        if len(primers) > _HISTORY_DETAIL_LIST_MAX:
            lines.append(
                f"  [dim]… (+{len(primers) - _HISTORY_DETAIL_LIST_MAX} more)[/]")
    # The annotation state this intermediate carried when the step ran.
    # Reported as a count, not re-rendered — the features are ~95% of a
    # real history's bytes and the map already draws the current ones.
    n_feats = hist.feature_snapshot_count
    if n_feats:
        lines.append("")
        lines.append("[b]Annotations at this step[/]")
        lines.append(f"  {n_feats:,} feature{'' if n_feats == 1 else 's'} "
                     "[dim]recorded on this intermediate[/dim]")
    # Backstop: anything this codebase does NOT model. An attribute or
    # child element from a newer writer, or a vendor extension, still
    # SHOWS here instead of being silently invisible. Normally empty.
    extras = hist.extra_attributes
    others = hist.other_children
    if extras or others:
        lines.append("")
        lines.append("[b]Other recorded fields[/]")
        for k in sorted(extras)[:_HISTORY_DETAIL_LIST_MAX]:
            lines.append(f"  {_esc(k)}: {_esc(str(extras[k])[:80])}")
        for tag, count in others[:_HISTORY_DETAIL_LIST_MAX]:
            suffix = f" ×{count}" if count > 1 else ""
            lines.append(f"  [dim]<{_esc(tag)}>{suffix}[/dim]")
    parents = hist.parents
    lines.append("")
    if parents:
        lines.append(f"[b]Parents ({len(parents)})[/]")
        shown_p = parents[:_HISTORY_DETAIL_LIST_MAX]
        joined = ", ".join(
            _esc(_history_clean_name(p.name)) for p in shown_p
        )
        if len(parents) > _HISTORY_DETAIL_LIST_MAX:
            joined += f", … (+{len(parents) - _HISTORY_DETAIL_LIST_MAX} more)"
        lines.append(f"  {joined}")
    else:
        lines.append("[dim](leaf — no recorded parents)[/]")
    return lines


_HISTORY_NODE_MAX_DEPTH = 500


_HISTORY_NODE_MAX_NODES = 100_000
