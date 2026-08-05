"""splicecraft_cloning — construction simulation (Phase D, layer L3).

The "simulate the real steps" construction helpers ([INV-127]): build real
amplicons (`_simulate_primed_amplicon`), assemble real cloned plasmids via
digest+ligation (`_simulate_cloned_plasmid`), the pUPD2 backbone stub, overhang
fusion, the Commercial-SaaS `.dna` history serialisation, and the **Gibson assembly
simulator** ([INV-85/86]): `_simulate_gibson_assembly` + its `_gibson_*` helpers
(overlap detect, body-length validate, product build, feature shift + origin-wrap
merge) and `_gibson_record_from_result`; plus the **traditional (restriction) cloning
simulator** ([INV-127]): the ligation primitives (`_ends_compatible`/`_ligate_fragments`
/`_close_circular`) and the cut-paste sims (`_simulate_traditional_cloning`(`_multi`),
`_classify_junction`, `_annotate_scars_on_product`, `_rc_fragment`,
`_label_disrupted_split_features`), with the fragment-prep that feeds them
(`_make_synthetic_fragment`, `_excise_pcr_insert`, `_excise_fragment_pair`,
`_enzyme_is_type_iis`); and the **Golden-Braid (BsaI Type IIS) fragment scrub**
([INV-127]): `_scrub_gb_design` + its `_scrub_gb_*` helpers + `_assemble_scrub_amplicons_real`
+ the `_SCRUB_GB_*` constants — cures sites by splitting at each cluster into BsaI-tailed
PCR fragments that Golden-Gate reassemble into the cured plasmid (real digest+ligate,
verified); and the **Golden-Braid / MoClo domestication primer designers**:
`_design_gb_primers` (single-part domestication, with internal-Type-IIS synonymous
repair) + `_design_operon_soe_primers` (whole-operon SOE) + the grammar-position /
forbidden-site / pair-name helpers (`_grammar_position_by_type`/`_position_overhangs`/
`_tu_overhangs`, `_gb_find_forbidden_hits`, `_gb_binding_region_advisory`,
`_dom_primer_pair_names`). Extracted so the cloning modal/screen siblings can import them.
Layer L3: imports state(L0), biology(L0), codon(L2), primer(L2), dataaccess(L1),
record(L1), history(L2), logging(L0); used by the modals (L4). The enzyme catalog is
reached via `_state._all_enzymes_hook` (it reads dataaccess but stays hub-side); the GB
scrub + designers reuse primer's `_scrub_design`/binding helpers + codon's `_codon_fix_*`
(all L2, downward — neither imports cloning). Re-exported by the hub so every call site
resolves unchanged.
"""
from __future__ import annotations

import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # annotation-only; the real Bio import is lazy, inside the fn
    from Bio.SeqRecord import SeqRecord

import splicecraft_state as _state
from splicecraft_biology import (
    _digest_with_enzymes, _enzyme_cuts, _forbidden_hit_set, _fragments_from_cuts,
    _iupac_pattern, _rc, _slice_circular,
)
from splicecraft_codon import _codon_fix_mutation_positions, _codon_fix_sites  # L2
from splicecraft_primer import (   # L2 (downward; primer never imports cloning)
    _PRIMER_MAX_OLIGO_LEN, _SCRUB_DEFAULT_ENZYMES, _SCRUB_PRIMER_FOOTPRINT,
    _binding_max_len, _circ_extract, _mut_gc_pct, _mut_translate, _pick_binding_region,
    _primer_tm, _scrub_cluster_edits, _scrub_cluster_span, _scrub_design,
    _scrub_resolve_sites,
)
from splicecraft_dataaccess import (
    _BUILTIN_GRAMMARS, _GB_CODING_PART_TYPES, _GB_DOMESTICATION_FORBIDDEN,
    _GB_L0_ENZYME_SITE, _GB_PAD, _GB_SPACER, _get_entry_vector,
)
from splicecraft_history import (
    _CommercialSaaSHistoryNode, _coerce_int_or_zero, _history_now_str,
)
from splicecraft_logging import _log, _timed
from splicecraft_record import _gb_text_to_record, _normalize_primer_seq


def _serialize_commercialsaas_history(root: "_CommercialSaaSHistoryNode | None"
                                  ) -> str:
    """Serialise a `_CommercialSaaSHistoryNode` back to UTF-8 XML text
    suitable for passing to `_pack_commercialsaas_history_payload`. Empty
    / None root yields an empty `<HistoryTree/>` document so we
    can still emit a valid history packet (e.g., to mark a
    construction with no parental input)."""
    import xml.etree.ElementTree as _ET
    tree_root = _ET.Element("HistoryTree")
    if root is not None:
        tree_root.append(root.element)
        # 2026-05-27 (audit-3 M5): also re-emit any sibling top-level
        # Node elements the parser stashed. Preserves round-trip for
        # rare multi-tree history files.
        for sibling in getattr(root, "_sibling_elements", []) or []:
            tree_root.append(sibling)
    body = _ET.tostring(tree_root, encoding="unicode")
    # CommercialSaaS's emitted XML always opens with the standard
    # declaration. Re-add it on the way out for compatibility — Python's
    # `tostring(encoding="unicode")` skips the declaration by default.
    return f'<?xml version="1.0" encoding="UTF-8"?>{body}'


def _build_amplicon_history_xml(*, name: str, seq_len: int,
                                 fwd_seq: str, rev_seq: str,
                                 start_1based: int, end_1based: int,
                                 fwd_name: str = "forward primer",
                                 rev_name: str = "reverse primer",
                                 ) -> "str | None":
    """Build a `<HistoryTree>` for a PCR amplicon made IN SpliceCraft
    that matches the element set a `.dna` import carries: an
    ``amplifyFragment`` node + an ``amplify`` InputSummary (val1/val2 =
    the amplified region) + the two primers as `<Oligo>` children. So a
    de-novo amplicon's History shows the same Primers / region detail as
    an imported one — harmonised history regardless of origin (user
    request 2026-06-01). Best-effort; returns None on failure."""
    try:
        root = _CommercialSaaSHistoryNode.new(
            name=str(name or "amplicon") + ".dna",
            seq_len=_coerce_int_or_zero(seq_len),
            circular=False, operation="amplifyFragment", node_id=0,
            date=_history_now_str(),
        )
        root.add_input_summary(
            manipulation="amplify",
            val1=_coerce_int_or_zero(start_1based),
            val2=_coerce_int_or_zero(end_1based),
        )
        if fwd_seq:
            root.add_oligo(name=fwd_name, sequence=str(fwd_seq))
        if rev_seq:
            root.add_oligo(name=rev_name, sequence=str(rev_seq))
        return _serialize_commercialsaas_history(root)
    except Exception:
        _log.debug("amplicon history: build failed for %r", name,
                    exc_info=True)
        return None


def _atg_offset_for_part(part_oh5: str, part_type: str) -> int:
    """Return the number of bases the 5' boundary of a coding-part
    feature should extend upstream to include its embedded start codon.

    GB 2.0 puts the ATG start codon inside the AATG fusion overhang
    (Pos 12→13 boundary): AATG = A + ATG, where the A is the spacer
    base and ATG is the start codon. The domesticator's forward primer
    encodes AATG as the part's 5' fusion overhang and PCR-binds at
    codon 2 of the source CDS — so the L0 part's body sequence
    (`part["sequence"]`) starts at codon 2, NOT at the ATG.

    When the part is assembled into an L1 plasmid the upstream LINK
    contributes its AATG (Pos 12 oh3) which fuses with the part's
    AATG oh5; the resulting cloned sequence reads
    ``...LINK-body...A + ATG + [codon2]...`` and the ATG sits in the
    last 3 nt of the fusion overhang. A feature annotation that only
    spans [body-start, body-end) would therefore drop the user's
    start codon — visibly broken in the L1 plasmid map.

    This helper returns ``3`` when the part needs the upstream
    extension (any coding part type with ``oh5 == AATG``: Signal
    peptide, CDS, CDS-NS, CDS-NS (CT)), and ``0`` otherwise. Returning
    a numeric offset rather than a bool lets callers do
    ``feature_start -= _atg_offset_for_part(oh5, ptype)`` without
    branching on every feature loop iteration.

    Regression guard for 2026-05-10: user reported "CDS's cloned seem
    to lose the annotation of their ATG because it also occupies the
    AATG overhang."
    """
    if not isinstance(part_oh5, str) or not isinstance(part_type, str):
        return 0
    oh = part_oh5.upper()
    # Coding types: known GB types + the MoClo-equivalent
    # ("CDS" / "C-tag") declared in the MoClo Plant grammar.
    # Custom grammars that introduce new translational part type
    # names should add them to `_GB_CODING_PART_TYPES` (the name
    # is historical — the set covers BOTH GB and MoClo by union).
    if part_type not in _GB_CODING_PART_TYPES:
        return 0
    # Detect ATG-fusion: any 4+ nt 5' overhang whose LAST 3 nt
    # are "ATG" carries the start codon (GB: AATG; hypothetical
    # NATG variants in custom grammars also qualify). MoClo
    # Plant's CDS oh5 is "AGGT" — last 3 = "GGT" ≠ ATG, so this
    # naturally returns 0 there. The `>= 4` guard enforces the
    # standard Type IIS overhang width (4 nt) and rejects a
    # degenerate 3-nt "ATG" overhang which isn't biologically
    # standard and would otherwise erroneously trigger the skip.
    if oh.endswith("ATG") and len(oh) >= 4:
        return 3
    return 0


def _fuse_overhang_body(oh5: str, body: str, part_type: str) -> str:
    """Join a 5' fusion overhang to a coding-part body, COLLAPSING the
    start-codon overlap.

    When the overhang embeds the start codon (Golden Braid / MoClo ``AATG``
    = ``A`` spacer + ``ATG`` start; any 4+ nt overhang ending in ATG) AND the
    stored body still carries that same ATG (the Domesticator-saved
    convention), drop the body's redundant leading ATG so the fused sequence
    reads ``AATG[codon2…]`` instead of ``AATG·ATG[codon2…]`` — a duplicated
    start codon that frameshifts the rest of the ORF (the 2026-05-30 'double
    ATG' bug). The designed forward primer already binds at codon 2, so this
    makes the simulated assembly match what the primers produce on the bench.

    Grammar-agnostic via `_atg_offset_for_part` — a no-op for MoClo-Plant
    ``AGGT``, non-coding parts, or a body already trimmed to codon 2 — so it
    is safe and idempotent at every site that assembles a part body behind
    its 5' overhang (Golden Braid AND MoClo)."""
    off = _atg_offset_for_part(oh5, part_type)
    if off and str(body)[:3].upper() == "ATG":
        body = body[off:]
    return oh5 + body


def _build_pupd2_backbone_stub(seed: int = 0xBACDBAC0, length: int = 420) -> str:
    """Return a deterministic ACGT string, free of BsaI/Esp3I/BsmBI sites on
    both strands, for use as a pUPD2-shaped placeholder backbone.

    Deterministic because the same insert must produce the same cloned
    sequence across sessions (otherwise the "Copy Cloned Sequence" output
    would silently drift). Seeded with a fixed constant.
    """
    import random as _random_mod
    rng = _random_mod.Random(seed)
    bases = [rng.choice("ACGT") for _ in range(length)]
    # Scrub both strands — the linear backbone becomes part of a circular
    # product, so a top-strand CGTCTC (Esp3I/BsmBI) and a bottom-strand
    # GAGACG are biologically equivalent and both must be absent.
    forbidden = ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG")
    i = 0
    while i <= length - 6:
        window = "".join(bases[i:i + 6])
        if window in forbidden:
            # Flip the middle base to something that can't re-hit any
            # forbidden site; ACGT minus the current base leaves 3 choices.
            middle = i + 3
            current = bases[middle]
            for replacement in "ACGT":
                if replacement != current:
                    bases[middle] = replacement
                    break
            # Rewind a bit to catch any new site created at the boundary.
            i = max(0, i - 5)
            continue
        i += 1
    return "".join(bases)


_PUPD2_BACKBONE_STUB: str = _build_pupd2_backbone_stub()


def _simulate_primed_amplicon(
    insert: str, oh5: str, oh3: str,
    grammar: "dict | None" = None,
    part_type: str = "",
    entry_overhangs: "tuple[str, str] | None" = None,
) -> str:
    """PCR amplicon top strand (5'→3'), as it would run on a pre-digest gel.

    Structure:  [pad] [enzyme site] [spacer] [ent5] [oh5] [insert] [oh3] [ent3]
                [rc(spacer)] [rc(enzyme site)] [rc(pad)]

    Matches the primer geometry in :func:`_design_gb_primers`. Defaults to
    Golden Braid L0 (Esp3I); pass ``grammar`` to use a different cloning
    grammar's enzyme/pad/spacer (e.g., MoClo Plant uses BsaI). Used by
    both DomesticatorModal (active grammar at design time) and
    PartsBinModal "Copy Primed Sequence" (the part's stored grammar).

    ``entry_overhangs`` = ``(external_oh5, external_oh3)`` nests a two-tier
    Golden-Braid UPD layout: the external entry pair (e.g. CTCG/TGAG) wraps the
    category overhangs, exactly as `_design_gb_primers` builds the primers, so
    the simulated amplicon equals what those primers produce on the bench.
    ``None`` — or a pair equal to ``(oh5, oh3)`` — keeps the classic one-tier
    amplicon (byte-identical to the pre-nesting output).
    """
    g = grammar if isinstance(grammar, dict) else _BUILTIN_GRAMMARS["gb_l0"]
    pad    = g.get("pad",    _GB_PAD)
    site   = g.get("site",   _GB_L0_ENZYME_SITE)
    spacer = g.get("spacer", _GB_SPACER)
    left_tail  = pad + site + spacer
    right_tail = _rc(spacer) + _rc(site) + _rc(pad)
    ent5 = ent3 = ""
    if isinstance(entry_overhangs, (tuple, list)) and len(entry_overhangs) == 2:
        e5 = str(entry_overhangs[0] or "").upper()
        e3 = str(entry_overhangs[1] or "").upper()
        if len(e5) == 4 and len(e3) == 4 \
                and (e5, e3) != ((oh5 or "").upper(), (oh3 or "").upper()):
            ent5, ent3 = e5, e3
    # `_fuse_overhang_body` collapses the AATG-CDS start-codon overlap so the
    # amplicon matches the designed primers (which bind at codon 2) — no
    # duplicated ATG. Grammar-agnostic (GB + MoClo); no-op otherwise. The
    # external entry overhangs (`ent5`/`ent3`) wrap that fused body.
    return (left_tail + ent5 + _fuse_overhang_body(oh5, insert, part_type)
            + oh3 + ent3 + right_tail)


def _simulate_cloned_plasmid(insert: str, oh5: str, oh3: str,
                             part_type: str = "") -> str:
    """Simulated cloned circular plasmid, linearised at the 5' overhang.

    After the cloning grammar's enzyme cuts both the amplicon and the
    backbone, the insert fragment carries `oh5…oh3` on its 4-nt sticky
    ends and ligates into the backbone in a single orientation. The
    circular product, read starting at `oh5`, is:

        [oh5] [insert] [oh3] [backbone_body]

    The backbone here is `_PUPD2_BACKBONE_STUB` — a scrubbed placeholder
    that contains no BsaI/Esp3I sites on either strand, so the simulated
    plasmid is guaranteed not to re-cut in either L0 or L1 assembly.
    """
    return (_fuse_overhang_body(oh5, insert, part_type)
            + oh3 + _PUPD2_BACKBONE_STUB)


# ── Gibson assembly ───────────────────────────────────────────────────────────
#
# Gibson chemistry: 5' exonuclease chews back the 5' end of each fragment,
# exposing a 3' single-stranded tail. Tails that share a complementary
# homology region anneal; DNA polymerase fills gaps; DNA ligase seals
# nicks. Net result: adjacent fragments are joined seamlessly and the
# homology region appears ONCE in the product.
#
# The simulator here treats each fragment as a top-strand DNA sequence
# (5' → 3') and finds the longest exact-match suffix-of-A / prefix-of-B
# overlap at each junction. The user is responsible for designing
# primer tails that produce these overlaps — the simulator doesn't
# pretend to extend short overlaps via heuristic search. Below
# `_GIBSON_MIN_OVERLAP_BP` (15 bp by default, Gibson's commonly cited
# floor), the junction is rejected so the user can fix the design
# before committing the assembled product to the library.

_GIBSON_MIN_OVERLAP_BP = 15
_GIBSON_MAX_OVERLAP_BP = 200   # cap the suffix/prefix probe; longer is
                                # unrealistic for Gibson primer tails


def _gibson_overlap_len(a_seq: str, b_seq: str, *,
                          min_overlap: int = _GIBSON_MIN_OVERLAP_BP,
                          max_overlap: int = _GIBSON_MAX_OVERLAP_BP,
                          ) -> int:
    """Length of the longest exact-match overlap between `a_seq`'s 3'
    end and `b_seq`'s 5' end, in `[min_overlap, max_overlap]`.
    Returns 0 if no overlap of at least `min_overlap` bp matches.

    Comparison is case-insensitive. Bigger overlap preferred —
    biologically, the homology arm the user designed is the longest
    exact match, not a degenerate short one inside it.

    When `a_seq` is identical to `b_seq` (e.g. the n=1 circular
    self-circularisation probe, where the simulator passes the same
    fragment as both a and b), the probe caps at `len(a) - 1` so the
    trivial whole-string match is skipped — without the cap, the user's
    intended short homology arm at the fragment ends would be masked
    by the always-matching full string.
    """
    if min_overlap <= 0:
        min_overlap = 1
    a = a_seq.upper()
    b = b_seq.upper()
    # Whole-string match is degenerate when the two sides are the same
    # sequence — cap one shorter so the probe finds a real arm overlap.
    # For distinct sequences a full-length match is biologically legal
    # (one fragment is a prefix/suffix of the other); the downstream
    # body-length validation in `_simulate_gibson_assembly` decides.
    full_match_safe = (a != b)
    max_check = min(max_overlap, len(a), len(b))
    if not full_match_safe:
        max_check = min(max_check, len(a) - 1)
    if max_check < min_overlap:
        return 0
    for k in range(max_check, min_overlap - 1, -1):
        if a[-k:] == b[:k]:
            return k
    return 0


@_timed("op.gibson_simulate")
# ── _simulate_gibson_assembly helpers ─────────────────────────────────────
#
# The Gibson simulator was 370 lines pre-refactor. The math is
# biology-critical (junction-overlap detection, body-length validation,
# wrap-aware feature shifting, wrap-half re-merge), so the helpers are
# defined as pure functions with explicit `(input) → output` signatures
# — no shared closure state, no hidden dependencies. Behaviour is bit-
# identical to the pre-refactor version, tested by tests/test_gibson.py
# (47 cases) which passed before AND after the extraction.


def _gibson_failure(circular: bool, errors: list[str],
                       overlaps: "list[dict] | None" = None,
                       warnings: "list[str] | None" = None) -> dict:
    """Standard failure-shaped result dict for the Gibson pipeline.

    Centralising it ensures every short-circuit returns the same shape
    (success=False, empty product_seq, empty features) — a property the
    UI consumer (`GibsonAssemblyPane`) relies on."""
    return {
        "success":     False,
        "product_seq": "",
        "circular":    circular,
        "features":    [],
        "overlaps":    list(overlaps or []),
        "errors":      list(errors),
        "warnings":    list(warnings or []),
    }


def _gibson_normalize_fragments(
    fragments: list[dict],
) -> "tuple[list[dict] | None, str | None]":
    """Uppercase + whitespace-strip each fragment's sequence.

    Returns ``(norm_fragments, None)`` on success or ``(None, error_msg)``
    on the two short-circuit failure modes: empty input list, or any
    fragment that isn't a dict.

    The normalised list is a parallel structure to `fragments` —
    name / sequence / features — with the sequence cleaned. We
    intentionally don't mutate the caller's dicts."""
    if not fragments:
        return None, "No fragments supplied."
    norm_fragments: list[dict] = []
    for f in fragments:
        if not isinstance(f, dict):
            return None, "Each fragment must be a dict."
        raw = str(f.get("sequence") or "")
        # Sweep #30 (2026-05-28): map RNA U->T to match the rest of the
        # app's normalisation. A fragment pasted in RNA notation would
        # otherwise fail overlap detection against a T-notation neighbour
        # (and a downstream _rc would mangle the U, which has no entry in
        # the complement table). Main sequence-load paths already do this.
        cleaned = "".join(
            ch for ch in raw.upper().replace("U", "T") if not ch.isspace()
        )
        norm_fragments.append({
            "name":     str(f.get("name") or "?"),
            "sequence": cleaned,
            "features": list(f.get("features") or []),
        })
    return norm_fragments, None


def _gibson_detect_overlaps(
    norm_fragments: list[dict], *, min_overlap: int, circular: bool,
) -> "tuple[list[dict], list[int], list[str]]":
    """For each junction (consecutive pair plus wrap when circular),
    detect the longest exact-match suffix/prefix overlap.

    Returns ``(overlaps, overlap_lens, errors)``:
      * `overlaps`: one dict per junction (junction, from, to, length,
        seq, ok, is_wrap, rc_hint) — the UI shows the full chain even
        on partial failure, so this is populated regardless of success.
      * `overlap_lens`: parallel list of int lengths for the build step.
      * `errors`: human-readable junction-failure messages.
    """
    n = len(norm_fragments)
    n_junctions = n if circular else n - 1
    overlaps: list[dict] = []
    overlap_lens: list[int] = []
    errors: list[str] = []
    for i in range(n_junctions):
        a = norm_fragments[i]
        b = norm_fragments[(i + 1) % n]
        a_seq = a["sequence"]
        b_seq = b["sequence"]
        k = _gibson_overlap_len(a_seq, b_seq, min_overlap=min_overlap)
        ok = k >= min_overlap
        rc_hint = _gibson_rc_hint(a, b, min_overlap=min_overlap) if not ok else ""
        overlaps.append({
            "junction": i + 1,
            "from":     a["name"],
            "to":       b["name"],
            "length":   k,
            "seq":      a_seq[-k:] if k else "",
            "ok":       ok,
            "is_wrap":  (circular and i == n - 1),
            "rc_hint":  rc_hint,
        })
        if not ok:
            errors.append(
                f"Junction {i+1} ({a['name']!r} → {b['name']!r}): "
                f"no overlap ≥ {min_overlap} bp." + rc_hint
            )
        overlap_lens.append(k)
    return overlaps, overlap_lens, errors


def _gibson_rc_hint(a: dict, b: dict, *, min_overlap: int) -> str:
    """Reverse-orientation diagnostic. If the forward overlap failed,
    probe RC(b) and RC(a) at a lower threshold (min of `min_overlap`
    and 10 bp so the user sees the hint even with a high min). Returns
    `""` if neither RC orientation yields a plausible overlap.

    Returning a string so the caller can concat into the error message
    without an extra branch — keeps the failure path linear."""
    a_seq = a["sequence"]
    b_seq = b["sequence"]
    probe_min = min(min_overlap, 10)
    k_b_rc = _gibson_overlap_len(a_seq, _rc(b_seq), min_overlap=probe_min)
    k_a_rc = _gibson_overlap_len(_rc(a_seq), b_seq, min_overlap=probe_min)
    if k_b_rc >= probe_min and k_b_rc >= k_a_rc:
        return (
            f" — but reverse-complement of {b['name']!r} "
            f"yields a {k_b_rc} bp overlap; "
            f"did you mean to flip {b['name']!r}?"
        )
    if k_a_rc >= probe_min:
        return (
            f" — but reverse-complement of {a['name']!r} "
            f"yields a {k_a_rc} bp overlap; "
            f"did you mean to flip {a['name']!r}?"
        )
    return ""


def _gibson_validate_body_lengths(
    norm_fragments: list[dict], overlap_lens: list[int], *, circular: bool,
) -> list[str]:
    """Check that every fragment has body bases left after its homology
    arms. A fragment fully consumed by overlap(s) is biologically
    redundant — its bases are entirely supplied by adjacent fragments.

    Returns a list of error messages; empty list means OK.

    Three sub-cases:
      * fragments[i>0]: leading overlap from junction i-1.
      * circular n==1: self-circularisation uses overlap_lens[0] as a
        single wrap arm.
      * circular n>1: last fragment has BOTH leading + trailing arms.
    """
    n = len(norm_fragments)
    errors: list[str] = []
    for i in range(1, n):
        oh_lead = overlap_lens[i - 1]
        frag_len = len(norm_fragments[i]["sequence"])
        if oh_lead >= frag_len:
            errors.append(
                f"Fragment {norm_fragments[i]['name']!r} is consumed "
                f"by its leading {oh_lead} bp overlap "
                f"(fragment is {frag_len} bp). "
                f"Use a longer fragment or shorter overlap."
            )
    if circular and not errors:
        if n == 1:
            wrap_oh = overlap_lens[0]
            frag_len = len(norm_fragments[0]["sequence"])
            if wrap_oh >= frag_len:
                errors.append(
                    f"Fragment {norm_fragments[0]['name']!r} is fully "
                    f"consumed by its self-circularisation overlap "
                    f"({wrap_oh} ≥ {frag_len} bp)."
                )
        else:
            last_lead  = overlap_lens[n - 2]
            last_trail = overlap_lens[n - 1]
            last_len = len(norm_fragments[-1]["sequence"])
            if last_lead + last_trail >= last_len:
                errors.append(
                    f"Fragment {norm_fragments[-1]['name']!r} is fully "
                    f"consumed by its homology arms ({last_lead} + "
                    f"{last_trail} ≥ {last_len} bp). Pick a longer "
                    f"fragment or shorter overlaps."
                )
            # Fragment 0 also carries BOTH arms in a circular assembly: a
            # leading wrap arm (junction (n-1)->0) and a trailing arm
            # (junction 0->1). The i>0 loop above never validates it.
            first_lead  = overlap_lens[n - 1]
            first_trail = overlap_lens[0]
            first_len = len(norm_fragments[0]["sequence"])
            if first_lead + first_trail >= first_len:
                errors.append(
                    f"Fragment {norm_fragments[0]['name']!r} is fully "
                    f"consumed by its homology arms ({first_lead} + "
                    f"{first_trail} ≥ {first_len} bp). Pick a longer "
                    f"fragment or shorter overlaps."
                )
    return errors


def _gibson_short_fragment_warnings(
    norm_fragments: list[dict], *, min_overlap: int,
) -> list[str]:
    """Soft warning: fragments < 3× min_overlap are legal but unusual,
    flag for awareness without blocking the save."""
    warnings: list[str] = []
    for f in norm_fragments:
        if 0 < len(f["sequence"]) < 3 * min_overlap:
            warnings.append(
                f"Fragment {f['name']!r} is short "
                f"({len(f['sequence'])} bp) relative to the "
                f"{min_overlap} bp homology arms — assembly may "
                f"be hard to confirm by gel."
            )
    return warnings


def _gibson_build_product(
    norm_fragments: list[dict], overlap_lens: list[int], *, circular: bool,
) -> "tuple[str, list[int]]":
    """Concatenate fragment bodies into the product sequence, tracking
    each fragment's offset in product coordinates. The trailing copy
    of each overlap is dropped (`frag[i+1][oh_lead:]`).

    Returns ``(product_seq, offsets)`` where offsets[i] is the product
    coord at which fragments[i]'s local-pos 0 lands. For wrap fragments
    this offset can be negative-equivalent (modulo product_len)."""
    n = len(norm_fragments)
    seq_parts: list[str] = []
    offsets: list[int] = []
    first_seq = norm_fragments[0]["sequence"]
    seq_parts.append(first_seq)
    offsets.append(0)
    cursor = len(first_seq)
    for i in range(1, n):
        oh_lead = overlap_lens[i - 1]
        frag_seq = norm_fragments[i]["sequence"]
        body = frag_seq[oh_lead:]
        seq_parts.append(body)
        # Fragment i's local-pos 0 maps to product pos (cursor - oh_lead)
        # — the leading overlap bases already exist as the previous
        # fragment's tail.
        offsets.append(cursor - oh_lead)
        cursor += len(body)
    if circular:
        if n == 1:
            # Self-circularisation: drop the trailing wrap overlap from
            # the only fragment. seq_parts has one entry == fragments[0].
            wrap_oh = overlap_lens[0]
            if wrap_oh > 0:
                seq_parts[0] = seq_parts[0][:-wrap_oh]
        else:
            wrap_oh = overlap_lens[n - 1]
            if wrap_oh > 0:
                # Drop the last fragment's trailing-overlap bases —
                # they equal the first fragment's leading bases, which
                # are already at the start of the product.
                seq_parts[-1] = seq_parts[-1][:-wrap_oh]
    product_seq = "".join(seq_parts)
    return product_seq, offsets


def _gibson_shift_features(
    norm_fragments: list[dict], overlap_lens: list[int],
    offsets: list[int], product_len: int, *, circular: bool,
) -> list[dict]:
    """Shift each fragment's features into product coords. Features
    wholly inside a fragment's leading-overlap region (i > 0) are
    skipped — the preceding fragment already supplies the same bases
    at the same product coords, so emitting them again duplicates the
    annotation. For circular, features that straddle the wrap junction
    become wrap features (`end < start`) per the dict-feature
    convention; `_feat_len` / `_bp_in` handle the resulting topology."""
    shifted: list[dict] = []
    for i, f_dict in enumerate(norm_fragments):
        offset = offsets[i]
        oh_lead = overlap_lens[i - 1] if i > 0 else 0
        for feat in (f_dict.get("features") or []):
            if not isinstance(feat, dict):
                continue
            if str(feat.get("type") or "") == "source":
                continue
            try:
                s = int(feat.get("start", 0))
                e = int(feat.get("end",   0))
            except (TypeError, ValueError):
                continue
            if e <= s:
                continue
            if i > 0 and e <= oh_lead:
                # Feature lies entirely in the leading-overlap region
                # — preceding fragment already supplies its annotation
                # at the same product coords. Skip to avoid a duplicate.
                continue
            new_s = offset + s
            new_e = offset + e
            span = new_e - new_s
            if circular and product_len > 0:
                # Wrap math: modulo into product coords. The span
                # (linear length) is invariant under shift, so we
                # decide product topology from `span` not from the
                # `ms <=> me` ordering — which is ambiguous when
                # `ms == me` mod product_len.
                ms = new_s % product_len
                me_raw = new_e % product_len
                # `me == 0` with span > 0 means the feature ends
                # exactly at the product's wrap point — keep it as
                # `product_len` so the linear-form expression matches.
                me = product_len if (me_raw == 0 and span > 0) else me_raw
                if span >= product_len:
                    new_s, new_e = 0, product_len
                else:
                    new_s, new_e = ms, me
            else:
                # Linear product. Negative `new_s` means the feature's
                # 5' edge sits before the product start — the only way
                # this happens is if a middle fragment is shorter than
                # its lead+trail arms (a pathological design the
                # simulator doesn't reject up front; see the body
                # length check above which only catches single-arm
                # exhaustion). Skip rather than clamp: clamping would
                # silently shift the feature's biological coordinates,
                # whereas skip is honest about the lost annotation.
                if new_s < 0:
                    continue
                if new_e > product_len:
                    new_e = product_len
                if new_e <= new_s:
                    continue
            shifted.append({**feat, "start": new_s, "end": new_e})
    return shifted


def _gibson_merge_wrap_halves(
    shifted: list[dict], product_len: int, *, circular: bool,
) -> list[dict]:
    """Re-merge `_wrap_pair`-tagged halves: when a circular source
    plasmid's wrap feature was split by `_record_features`, both halves
    carry the same ``_wrap_pair`` id + ``_wrap_role`` (head/tail). After
    shifting, if the two halves are still adjacent across the product's
    wrap (head.start == 0 AND tail.end == product_len) we collapse them
    back into one wrap feature. For linear products or non-adjacent
    halves the split is preserved — the biological feature was severed
    by the assembly geometry."""
    out_feats: list[dict] = []
    if not (circular and product_len > 0):
        # Linear product: just strip wrap-pair sentinels and pass through.
        for f in shifted:
            out_feats.append({
                k: v for k, v in f.items()
                if not k.startswith("_wrap_")
            })
        return out_feats
    # Index pairs in one pass so we can match heads to tails O(1) below.
    pair_index: dict[str, list[dict]] = {}
    for f in shifted:
        pid = f.get("_wrap_pair")
        if pid:
            pair_index.setdefault(pid, []).append(f)
    handled_pairs: set[str] = set()
    for f in shifted:
        pid = f.get("_wrap_pair")
        if pid and pid not in handled_pairs:
            halves = pair_index.get(pid) or []
            if len(halves) == 2:
                head = next((h for h in halves
                             if h.get("_wrap_role") == "head"), None)
                tail = next((h for h in halves
                             if h.get("_wrap_role") == "tail"), None)
                if (head is not None and tail is not None
                        and head["start"] == 0
                        and tail["end"] == product_len):
                    merged = {k: v for k, v in tail.items()
                              if not k.startswith("_wrap_")}
                    merged["start"] = tail["start"]
                    merged["end"]   = head["end"]
                    out_feats.append(merged)
                    handled_pairs.add(pid)
                    continue
        # Either no wrap pair, pair already handled, or halves not
        # adjacent — strip wrap sentinels and pass through individually.
        if pid in handled_pairs:
            continue
        out_feats.append({
            k: v for k, v in f.items()
            if not k.startswith("_wrap_")
        })
    return out_feats


def _simulate_gibson_assembly(fragments: list[dict], *,
                                min_overlap: int = _GIBSON_MIN_OVERLAP_BP,
                                circular: bool = True,
                              ) -> dict:
    """Simulate a Gibson assembly of N linear top-strand fragments.

    Each ``fragments[i]`` dict shape:
        {
          "name":     str,
          "sequence": str,            # linear DNA 5' → 3'
          "features": list[dict],     # optional, in fragment-local coords
        }

    For each junction (consecutive pair plus the wrap junction when
    ``circular=True``) the longest exact-match suffix/prefix overlap
    is detected. Any junction below ``min_overlap`` bp fails the
    assembly. The product sequence has each overlap appearing once
    (the trailing copy is dropped: ``frag[i] + frag[i+1][overlap:]``).
    Features are shifted into product coordinates; features wholly
    inside a fragment's leading-overlap region are skipped (they're
    already represented by the preceding fragment's trailing copy).

    Returns ``{success, product_seq, circular, features, overlaps,
    errors, warnings}`` — see UI consumer (``GibsonAssemblyPane``)
    for the rendering convention. ``overlaps`` always has one entry
    per junction so the UI can show the full chain even on partial
    failure.

    Pre-refactor this was a 370-line monolith; the per-stage logic
    now lives in `_gibson_*` helpers (normalise → detect overlaps →
    validate body lengths → build product seq → shift features →
    re-merge wrap halves). Behaviour is bit-identical; tested by
    tests/test_gibson.py (47 cases).
    """
    norm_fragments, norm_err = _gibson_normalize_fragments(fragments)
    if norm_err is not None:
        return _gibson_failure(circular, [norm_err])
    assert norm_fragments is not None  # err None ⇒ list present

    # Reject zero-length fragments up front — they can't carry homology.
    for f in norm_fragments:
        if not f["sequence"]:
            return _gibson_failure(
                circular, [f"Fragment {f['name']!r} has no sequence."]
            )

    if len(norm_fragments) == 1 and not circular:
        # A single linear fragment doesn't need Gibson — pass through.
        f = norm_fragments[0]
        return {
            "success":     True,
            "product_seq": f["sequence"],
            "circular":    False,
            "features":    list(f["features"]),
            "overlaps":    [],
            "errors":      [],
            "warnings":    [
                "Single linear fragment — no Gibson junctions to "
                "validate. Product is the fragment as supplied.",
            ],
        }

    overlaps, overlap_lens, junction_errors = _gibson_detect_overlaps(
        norm_fragments, min_overlap=min_overlap, circular=circular,
    )
    if junction_errors:
        return _gibson_failure(circular, junction_errors, overlaps=overlaps)

    body_errors = _gibson_validate_body_lengths(
        norm_fragments, overlap_lens, circular=circular,
    )
    if body_errors:
        return _gibson_failure(circular, body_errors, overlaps=overlaps)

    warnings = _gibson_short_fragment_warnings(
        norm_fragments, min_overlap=min_overlap,
    )

    product_seq, offsets = _gibson_build_product(
        norm_fragments, overlap_lens, circular=circular,
    )
    product_len = len(product_seq)

    shifted = _gibson_shift_features(
        norm_fragments, overlap_lens, offsets, product_len,
        circular=circular,
    )
    out_feats = _gibson_merge_wrap_halves(
        shifted, product_len, circular=circular,
    )

    return {
        "success":     True,
        "product_seq": product_seq,
        "circular":    circular,
        "features":    out_feats,
        "overlaps":    overlaps,
        "errors":      [],
        "warnings":    warnings,
    }


def _gibson_record_from_result(result: dict, *, name: str) -> "SeqRecord | None":
    """Build a SeqRecord from a successful ``_simulate_gibson_assembly``
    result. Returns ``None`` when ``result["success"] is False``.

    Wrap features (``end < start``) become ``CompoundLocation`` per the
    GenBank wrap convention. Linear features land as ``FeatureLocation``.
    Strand defaults to +1 when missing.
    """
    if not result or not result.get("success"):
        return None
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import (
        SeqFeature, FeatureLocation, CompoundLocation,
    )
    seq_str = str(result.get("product_seq") or "")
    n = len(seq_str)
    safe_name = re.sub(r"[^A-Za-z0-9_]+", "_", name or "gibson") or "gibson"
    topology = "circular" if result.get("circular") else "linear"
    rec = SeqRecord(
        Seq(seq_str),
        id=safe_name[:16] or "gibson",
        name=safe_name[:16] or "gibson",
        description=name or "Gibson assembly",
        annotations={
            "molecule_type": "DNA",
            "topology":      topology,
        },
    )
    for f in (result.get("features") or []):
        if not isinstance(f, dict):
            continue
        try:
            s = int(f.get("start", 0))
            e = int(f.get("end",   0))
            strand = int(f.get("strand", 1) or 1)
        except (TypeError, ValueError):
            continue
        if s == e or n == 0:
            continue
        ftype = str(f.get("type") or "misc_feature")
        quals: dict = {}
        label = f.get("label")
        if label:
            quals["label"] = [str(label)]
        color = f.get("color")
        if color:
            quals["ApEinfo_fwdcolor"] = [str(color)]
            quals["ApEinfo_revcolor"] = [str(color)]
        note = f.get("note")
        if note:
            quals["note"] = [str(note)]
        # Re-emit the primer sequence so an inherited primer_bind renders
        # bound bases + 5' flap on the Gibson product (not a plain bar) —
        # `_record_features` carried it through the assembly.
        ps = f.get("primer_seq")
        if ps:
            quals["primer_seq"] = [_normalize_primer_seq(ps)]
        if e > s:
            loc = FeatureLocation(s, e, strand=strand)
        else:
            # Wrap feature: (s, n) + (0, e)
            loc = CompoundLocation([
                FeatureLocation(s, n, strand=strand),
                FeatureLocation(0, e, strand=strand),
            ])
        rec.features.append(SeqFeature(loc, type=ftype, qualifiers=quals))
    return rec


# ═══ Traditional (restriction) cloning simulation — moved from the hub ═══════
# Ligation primitives + the cut-paste cloning sims ([INV-127]: the design IS
# the product). Enzyme catalog via _state._all_enzymes_hook.

def _label_disrupted_split_features(features: "list[dict]",
                                     enzymes: "list[str] | None" = None) -> None:
    """In-place: a feature split across a cloning cut (tagged ``_split`` head/
    tail by `_split_features_at_cuts`) had a cut site land INSIDE it — e.g.
    cloning into lacZα's MCS knocks the gene out. Mark each surviving half as
    disrupted: append `` (disrupted)`` to its label (once) and record a note
    naming the cut enzyme(s), so the product shows the broken feature as two
    flanking pieces instead of an intact one. Render-only — coords untouched.
    Idempotent (skips a half already marked); leaves un-split features (the
    carried-over insert annotations) alone."""
    enz = sorted({e for e in (enzymes or []) if e})
    where = f" ({'/'.join(enz)} cut site inside it)" if enz else ""
    for f in features:
        # head/tail = split across fragments; whole = cut(s) inside it but the
        # remnant stayed in one fragment (an excised middle, e.g. a lacZ MCS);
        # mid = this whole fragment lies INSIDE the feature, so it was cut on
        # BOTH sides — disrupted by definition (the excised lacZ-MCS stuffer
        # itself, which endpoint-slotting used to drop entirely).
        if (f.get("_split") not in ("head", "tail", "whole", "mid")
                or f.get("_disrupted")):
            continue
        base = str(f.get("label") or f.get("type") or "feature")
        if "(disrupted)" not in base:
            f["label"] = f"{base} (disrupted)"
        tag = f"Disrupted by the cloning insertion{where}."
        note = str(f.get("note") or "").strip()
        f["note"] = f"{note}; {tag}" if note else tag
        f["_disrupted"] = True


def _ends_compatible(end_a: dict, end_b: dict) -> bool:
    """Return True if two fragment edges can ligate. Same kind +
    matching overhang sequence (both stored top-strand-canonical). A
    `linear` edge never ligates."""
    ka, kb = end_a.get("kind"), end_b.get("kind")
    if ka == "linear" or kb == "linear":
        return False
    if ka != kb:
        return False
    return end_a.get("overhang_seq", "") == end_b.get("overhang_seq", "")


def _ligate_fragments(a: dict, b: dict) -> "dict | None":
    """Ligate `a.right` to `b.left`. Returns merged fragment (linear),
    or None if the overhangs are incompatible. The merged top strand
    is `a.top_seq + b.top_seq` (overhang bases live in whichever piece's
    top strand carried them — see the canonicalisation comment above).
    Features from `b` are shifted by `len(a.top_seq)`."""
    if not _ends_compatible(a["right"], b["left"]):
        return None
    shift = len(a["top_seq"])
    merged_feats = list(a["features"])
    for f in b["features"]:
        merged_feats.append({
            **f,
            "start": int(f.get("start", 0)) + shift,
            "end":   int(f.get("end",   0)) + shift,
        })
    return {
        "top_seq":      a["top_seq"] + b["top_seq"],
        "left":         a["left"],
        "right":        b["right"],
        "features":     merged_feats,
        "source_label": (f"{a['source_label']}+{b['source_label']}"
                         if a["source_label"] or b["source_label"] else ""),
    }


def _close_circular(frag: dict) -> "dict | None":
    """Close a linear fragment into a circle by ligating its right + left
    edges. Returns ``{top_seq, features, source_label, circular: True}``
    or None if the ends don't match."""
    if not _ends_compatible(frag["right"], frag["left"]):
        return None
    return {
        "top_seq":      frag["top_seq"],
        "features":     [dict(f) for f in frag["features"]],
        "source_label": frag["source_label"],
        "circular":     True,
    }


@_timed("op.simulate_traditional_cloning")
def _simulate_traditional_cloning(insert_frag: dict,
                                    vector_frag: dict,
                                   ) -> dict:
    """Try to ligate `insert_frag` into `vector_frag` in both possible
    orientations, returning a result dict:

      ``{
          "forward":  {"top_seq", "features", "compatible": bool},
          "reverse":  {"top_seq", "features", "compatible": bool},
          "warnings": [str, ...],
          "errors":   [str, ...],
        }``

    `forward` = vector + insert as supplied; `reverse` = vector + RC of
    insert. `compatible` is True when the overhangs actually permit
    that orientation; False means the orientation is rendered for
    reference (the insert can be flipped post-hoc) but isn't reachable
    by canonical ligation chemistry. When BOTH orientations are
    compatible (common with palindromic single-enzyme cuts), a
    warning calls out the ambiguity."""
    warnings: list[str] = []
    errors:   list[str] = []
    insert_rc = _rc_fragment(insert_frag)

    fwd_linear = _ligate_fragments(vector_frag, insert_frag)
    rev_linear = _ligate_fragments(vector_frag, insert_rc)

    fwd_compat = fwd_linear is not None and \
        _ends_compatible(fwd_linear["right"], fwd_linear["left"])
    rev_compat = rev_linear is not None and \
        _ends_compatible(rev_linear["right"], rev_linear["left"])

    fwd_seq = (fwd_linear["top_seq"] if fwd_linear is not None
               else vector_frag["top_seq"] + insert_frag["top_seq"])
    rev_seq = (rev_linear["top_seq"] if rev_linear is not None
               else vector_frag["top_seq"] + insert_rc["top_seq"])

    def _shift_feats(base: dict, shift: int) -> list[dict]:
        return [{**f,
                  "start": int(f.get("start", 0)) + shift,
                  "end":   int(f.get("end",   0)) + shift}
                 for f in base["features"]]

    # Cloning enzymes at the junctions — name the cut in the "(disrupted)"
    # note for any vector feature a cut site landed inside.
    _junction_enz = [e for e in (
        insert_frag.get("left",  {}).get("enzyme"),
        insert_frag.get("right", {}).get("enzyme"),
        vector_frag.get("left",  {}).get("enzyme"),
        vector_frag.get("right", {}).get("enzyme"),
    ) if e]
    # Fresh per-orientation copies of the vector features (don't mutate the
    # caller's fragment) so each can be labelled independently.
    fwd_feats = [dict(f) for f in vector_frag["features"]] + _shift_feats(
        insert_frag, len(vector_frag["top_seq"]))
    rev_feats = [dict(f) for f in vector_frag["features"]] + _shift_feats(
        insert_rc,   len(vector_frag["top_seq"]))
    # A cut site that fell inside a vector feature (e.g. cloning into lacZα's
    # MCS) split it into two halves — surface that as "(disrupted)".
    _label_disrupted_split_features(fwd_feats, _junction_enz)
    _label_disrupted_split_features(rev_feats, _junction_enz)

    if fwd_compat and rev_compat:
        warnings.append(
            "Ambiguous orientation: both forward and reverse ligation "
            "are chemically compatible. The cloning reaction will yield "
            "a mixture; pick by sequencing.")
    elif fwd_compat and not rev_compat:
        warnings.append(
            "Directional cloning: only the forward orientation is "
            "biologically achievable. Reverse is rendered for reference "
            "but cannot ligate.")
    elif rev_compat and not fwd_compat:
        warnings.append(
            "Directional cloning: only the reverse orientation is "
            "biologically achievable. Forward is rendered for reference "
            "but cannot ligate.")
    else:
        errors.append(
            "Neither orientation has matching overhangs at both junctions. "
            "Check that the insert and vector were cut with the same "
            "enzyme(s).")

    return {
        "forward": {"top_seq": fwd_seq, "features": fwd_feats,
                     "compatible": fwd_compat},
        "reverse": {"top_seq": rev_seq, "features": rev_feats,
                     "compatible": rev_compat},
        "warnings": warnings,
        "errors":   errors,
    }


def _classify_junction(left_enz: str, right_enz: str,
                          context_top: str,
                          *, context_left_offset: int = 6) -> dict:
    """Classify a ligation junction by checking whether the parent
    enzymes can still re-cut the joint sequence.

    ``context_top`` is the ~12 bp window straddling the junction
    (default 6 bp on each side; `context_left_offset` says how many
    bases of `context_top` are upstream of the cut). Returns:

      {
        "scar":         bool,    # True when NEITHER parent enzyme
                                 # recognises the joint — irreversible
                                 # (BioBrick-style idempotent assembly:
                                 # SpeI A/CTAGT + XbaI T/CTAGA → ACTAGA,
                                 # neither cuttable).
        "re_cuttable":  list[str],  # parent enzyme names whose
                                    # recognition site IS still
                                    # present at the junction.
        "label":        str,     # human-readable badge for warnings /
                                 # feature annotations.
      }

    Generic to all enzymes (palindromic + asymmetric + Type IIS):
    matches each parent enzyme's recognition site (via IUPAC pattern
    + reverse-complement) against the junction context window, so
    BamHI/BamHI (re-cuttable), BamHI/BglII (scar), and
    BsaI-Type-IIS junctions all classify correctly.
    """
    enzymes = []
    if left_enz:
        enzymes.append(left_enz)
    if right_enz and right_enz != left_enz:
        enzymes.append(right_enz)
    catalog = _state._all_enzymes_hook()
    re_cuttable: list[str] = []
    ctx = context_top.upper()
    for ename in enzymes:
        spec = catalog.get(ename)
        if spec is None:
            continue
        site = spec[0].upper()
        if not site:
            continue
        pat = _iupac_pattern(site)
        if pat.search(ctx):
            re_cuttable.append(ename)
            continue
        # Check the reverse complement too — asymmetric / Type IIS
        # enzymes can bind either strand and re-cut from the other
        # side. Palindromic sites are their own RC so a double-match
        # is fine.
        rc_pat = _iupac_pattern(_rc(site))
        if rc_pat.search(ctx):
            re_cuttable.append(ename)
    if re_cuttable:
        return {
            "scar":        False,
            "re_cuttable": re_cuttable,
            "label":       (f"{'/'.join(re_cuttable)} "
                              f"re-cuttable junction"),
        }
    # No parent enzyme site survives → idempotent scar.
    pair = f"{left_enz}/{right_enz}" if left_enz != right_enz else left_enz
    return {
        "scar":        True,
        "re_cuttable": [],
        "label":       f"{pair} scar (uncuttable by parent enzymes)",
    }


@_timed("op.simulate_traditional_cloning_multi")
def _simulate_traditional_cloning_multi(insert_frags: list[dict],
                                          vector_frag: dict) -> dict:
    """N-way version of `_simulate_traditional_cloning`. Pre-chains
    ``insert_frags`` in lane order — each adjacent pair must have
    matching sticky ends — then delegates to the 2-fragment engine
    for the final vector + chained-insert ligation. For N=1 this is
    bit-identical to calling the 2-fragment engine directly.

    Forward = inserts in lane order; Reverse = the entire chain RC'd
    (the 2-fragment engine handles the flip). Per-junction error
    messages name the failing pair by its ``source_label`` so the
    user can diagnose which sticky-end pair didn't match.

    Returns the same shape as `_simulate_traditional_cloning`:
    ``{"forward": {...}, "reverse": {...}, "warnings": [...],
       "errors": [...]}``.
    """
    if not insert_frags:
        empty = {"top_seq": vector_frag.get("top_seq", ""),
                 "features": [], "compatible": False}
        return {
            "forward":  empty,
            "reverse":  empty,
            "warnings": [],
            "errors":   ["No insert fragments queued for ligation."],
        }
    if len(insert_frags) == 1:
        result = _simulate_traditional_cloning(insert_frags[0],
                                                 vector_frag)
        _annotate_scars_on_product(result, insert_frags, vector_frag)
        return result
    chained = insert_frags[0]
    # Record junction info as we chain so the scar annotator below
    # can locate each junction in the final product.
    junction_info: list[dict] = []
    for i in range(1, len(insert_frags)):
        left_enz_at_junc  = chained["right"].get("enzyme") or ""
        right_enz_at_junc = insert_frags[i]["left"].get("enzyme") or ""
        junction_pos = len(chained["top_seq"])  # bp position in chain
        nxt = _ligate_fragments(chained, insert_frags[i])
        if nxt is None:
            a_label = chained.get("source_label") or f"fragment {i}"
            b_label = (insert_frags[i].get("source_label")
                        or f"fragment {i + 1}")
            empty = {"top_seq": "", "features": [], "compatible": False}
            return {
                "forward":  empty,
                "reverse":  empty,
                "warnings": [],
                "errors": [
                    f"Junction {i} → {i + 1}: sticky ends incompatible "
                    f"between {a_label!r} (3' end) and {b_label!r} "
                    f"(5' end). Check that adjacent fragments share an "
                    f"enzyme at the matching cut."
                ],
            }
        junction_info.append({
            "label":       f"insert {i} ↔ insert {i + 1}",
            "left_enz":    left_enz_at_junc,
            "right_enz":   right_enz_at_junc,
            "pos_in_chain": junction_pos,
        })
        chained = nxt
    result = _simulate_traditional_cloning(chained, vector_frag)
    _annotate_scars_on_product(result, insert_frags, vector_frag,
                                  internal_junctions=junction_info)
    return result


def _annotate_scars_on_product(
    result: dict,
    insert_frags: list[dict],
    vector_frag: dict,
    *,
    internal_junctions: "list[dict] | None" = None,
) -> None:
    """Walk every ligation junction in the simulated product and
    emit (a) a warning per junction describing whether it's
    re-cuttable or an idempotent scar (BioBrick-style), and (b) a
    ``misc_feature`` at the junction position labelled with the
    same. The annotations land on BOTH the forward and reverse
    orientation products so the user's saved plasmid carries the
    scar info regardless of which orientation they keep.

    ``internal_junctions`` is a list of ``{label, left_enz,
    right_enz, pos_in_chain}`` for the N-1 insert↔insert junctions
    when N inserts were chained pre-vector-ligation. The vector↔
    chain junctions are computed here from the parent frags'
    end-enzyme metadata. Mutates `result["warnings"]` and
    `result["forward"/"reverse"]["features"]` in place."""
    warnings: list[str] = result.setdefault("warnings", [])
    # Chained-insert sequence (everything except the vector).
    chain_top = "".join(f.get("top_seq", "") for f in insert_frags)
    chain_len = len(chain_top)
    vec_len   = len(vector_frag.get("top_seq", ""))
    # Per-orientation junctions are built lazily inside
    # `_build_orient_junctions(reverse, total)` below — forward and
    # reverse need different position math (the close-junction wraps
    # at the product's actual length, not at `vec_len + chain_len`).
    _ = chain_len  # kept for future re-use; no longer needed here

    def _build_orient_junctions(reverse: bool,
                                   total: int) -> list[dict]:
        """Per-orientation junction list. In REVERSE the insert chain
        is RC'd, so each chain end uses the OPPOSITE end's enzyme
        (RC swaps left↔right metadata). Positions also need
        re-derivation: the close lives at the product's wrap point
        (total), not at `vec_len + forward_chain_len`."""
        out: list[dict] = []
        if not reverse:
            # Forward = original layout.
            out.append({
                "label":     "vector ↔ insert 1",
                "left_enz":  vector_frag["right"].get("enzyme") or "",
                "right_enz": insert_frags[0]["left"].get("enzyme") or "",
                "pos":       vec_len,
            })
            for j in internal_junctions or []:
                out.append({
                    "label":     j["label"],
                    "left_enz":  j["left_enz"],
                    "right_enz": j["right_enz"],
                    "pos":       vec_len + j["pos_in_chain"],
                })
            out.append({
                "label":     f"insert {len(insert_frags)} ↔ vector",
                "left_enz":  insert_frags[-1]["right"].get("enzyme") or "",
                "right_enz": vector_frag["left"].get("enzyme") or "",
                "pos":       total,
            })
            return out
        # Reverse — chain order is REVERSED (insert N comes first,
        # insert 1 last) and each insert's left/right enzymes swap.
        n_inserts = len(insert_frags)
        # vec.right ↔ RC(chain_first).left = chain_first was
        # insert_frags[-1] (reversed order), so its left after RC
        # corresponds to original right.
        first_in_rc_chain = insert_frags[-1]
        out.append({
            "label":     f"vector ↔ insert {n_inserts}",
            "left_enz":  vector_frag["right"].get("enzyme") or "",
            "right_enz": first_in_rc_chain["right"].get("enzyme") or "",
            "pos":       vec_len,
        })
        # Internal junctions (only when N > 1): the chain runs in
        # reverse order, and each insert-to-insert junction uses the
        # RC'd ends. For an internal junction `insert i ↔ insert
        # i+1` in forward, the reverse equivalent is `insert (n - i)
        # ↔ insert (n - i - 1)` with swapped enzymes.
        if internal_junctions:
            # Walk forward chain junctions backwards.
            for jf in reversed(internal_junctions):
                out.append({
                    "label":     jf["label"] + " (reverse)",
                    # Enzymes swap because chain is RC'd
                    "left_enz":  jf["right_enz"],
                    "right_enz": jf["left_enz"],
                    "pos":       vec_len + (sum(
                        len(f.get("top_seq", ""))
                        for f in insert_frags
                    ) - jf["pos_in_chain"]),
                })
        # Closing junction (chain last → vec). chain_last was
        # insert_frags[0] in original order.
        last_in_rc_chain = insert_frags[0]
        out.append({
            "label":     "insert 1 ↔ vector",
            "left_enz":  last_in_rc_chain["left"].get("enzyme") or "",
            "right_enz": vector_frag["left"].get("enzyme") or "",
            "pos":       total,
        })
        return out

    def _annotate_orient(prod: dict, *, reverse: bool) -> None:
        if not prod.get("compatible", False):
            return
        top = prod.get("top_seq", "")
        if not top:
            return
        feats = prod.setdefault("features", [])
        total = len(top)
        for j in _build_orient_junctions(reverse, total):
            pos = j["pos"]
            # Circular wrap window — the closing junction sits at
            # `pos = total`. The context must straddle the linear-
            # string boundary or the regenerated parent recognition
            # site disappears.
            if pos >= total or pos == 0:
                pre = top[max(0, total - 6):total]
                post = top[:min(6, total)]
                context = pre + post
                ctx_left_offset = len(pre)
            else:
                window_l = max(0, pos - 6)
                window_r = min(total, pos + 6)
                context = top[window_l:window_r]
                ctx_left_offset = pos - window_l
            cls = _classify_junction(
                j["left_enz"], j["right_enz"], context,
                context_left_offset=ctx_left_offset,
            )
            warnings.append(f"Junction {j['label']}: {cls['label']}")
            # Tag the 4 bp ligation OVERHANG — light-blue, arrowless (strand
            # 0) — instead of labelling the junction a "LIGATION SCAR" (the
            # user wanted scars left as-is in the sequence, not annotated; the
            # re-cuttable / scar classification still rides the warnings
            # above). The ORIGIN junction's overhang straddles the
            # linearisation point, so tag it as a wrap feature (end < start →
            # CompoundLocation on save) — a full 4 bp, not the 2 bp head a
            # flat [0,2) clamp gives (review F6).
            if pos >= total or pos == 0:
                if total >= 4:
                    feats.append({
                        "start":  total - 2,
                        "end":    2,
                        "strand": 0,
                        "type":   "misc_feature",
                        "label":  (top[total - 2:total] + top[:2]).upper()
                                  or "overhang",
                        "color":  "#ADD8E6",
                    })
                continue
            feat_s = max(0, pos - 2)
            feat_e = min(total, pos + 2)
            if feat_e > feat_s:
                feats.append({
                    "start":  feat_s,
                    "end":    feat_e,
                    "strand": 0,
                    "type":   "misc_feature",
                    "label":  top[feat_s:feat_e].upper() or "overhang",
                    "color":  "#ADD8E6",
                })

    _annotate_orient(result.get("forward", {}), reverse=False)
    _annotate_orient(result.get("reverse", {}), reverse=True)


def _rc_fragment(frag: dict) -> dict:
    """Reverse-complement a Fragment, swapping its left/right ends and
    flipping its feature coordinates. The overhang sequences of the
    swapped ends are themselves reverse-complemented (the strand that
    sticks out is now the other strand). 5' overhangs stay 5' (the
    "5'-protruding" geometry is preserved across the flip — only the
    bases change).

    Convention-aware ``top_seq`` reconstruction (fix for the
    EcoRI+KpnI reverse-orientation scar-detection bug, 2026-05-23):
    excised fragments include overhang bases in ``top_seq`` at ends
    where the overhang is on the TOP strand (5' at left, 3' at
    right); synthetic fragments do not. After RC, the strand the
    overhang sits on flips, so the included/excluded bases must move
    accordingly — naive ``_rc(top_seq)`` produces junk at the
    junction for excise-convention fragments. The heuristic below
    detects per-end whether the convention is excise (overhang bases
    present in top_seq) by comparing the prefix/suffix of top_seq to
    the overhang_seq, and rebuilds the new top with the right
    strand-side overhang inclusion."""
    n = len(frag["top_seq"])
    top = frag["top_seq"]
    left  = frag["left"]
    right = frag["right"]
    left_oh  = left.get("overhang_seq", "") or ""
    right_oh = right.get("overhang_seq", "") or ""
    left_kind  = left.get("kind", "")
    right_kind = right.get("kind", "")
    # Convention detection: excise fragments include overhang bases
    # in top_seq at ends where the overhang is on the top strand
    # (5' at left OR 3' at right). Synthetic fragments
    # (`_make_synthetic_fragment`) never include them. Heuristic:
    # if top_seq's prefix/suffix matches the overhang at any on-top
    # end, the fragment is excise; otherwise synthetic. When neither
    # end is on-top (3' at left + 5' at right), heuristic can't
    # tell — default to excise (the common case from
    # `_excise_fragment_pair`).
    # Per-end positive checks: top_seq prefix/suffix matches the
    # overhang where the overhang sits on the TOP strand. A match is
    # a strong excise indicator; an end that COULD be on-top but
    # doesn't match is a strong synthetic indicator.
    can_check_left  = (left_kind == "5'" and bool(left_oh))
    can_check_right = (right_kind == "3'" and bool(right_oh))
    excise_match_left  = (can_check_left
        and top[:len(left_oh)].upper() == left_oh.upper())
    excise_match_right = (can_check_right
        and top[n - len(right_oh):].upper() == right_oh.upper())
    synth_match_left  = can_check_left  and not excise_match_left
    synth_match_right = can_check_right and not excise_match_right
    if excise_match_left or excise_match_right:
        is_excise = True
    elif synth_match_left or synth_match_right:
        is_excise = False
    else:
        # No on-top ends exist (both 3'-at-left or 5'-at-right) —
        # can't detect from top_seq. Default to excise (the common
        # case for `_excise_fragment_pair` output).
        is_excise = True
    if not is_excise:
        # Synthetic convention — preserve the pre-fix behaviour
        # (naive RC of top_seq, swap ends, no overhang-side
        # adjustment). The synthetic ligation path has its own
        # quirks but that's a separate bug.
        new_top = _rc(top)
        left_strip = 0
        right_strip = 0
        new_left_extra_len = 0
        new_right_extra_len = 0
    else:
        left_strip = len(left_oh) if excise_match_left else 0
        right_strip = len(right_oh) if excise_match_right else 0
        core_top = top[left_strip:n - right_strip] if right_strip else \
                   top[left_strip:]
        rc_core = _rc(core_top)
        # After RC, old.right (overhang on bot if 5') contributes
        # bases at new.left's top — prepend RC'd overhang. Old.left
        # (overhang on bot if 3') contributes at new.right's top —
        # append RC'd overhang.
        new_left_extra = (_rc(right_oh)
            if right_kind == "5'" and right_oh and right_strip == 0
            else "")
        new_right_extra = (_rc(left_oh)
            if left_kind == "3'" and left_oh and left_strip == 0
            else "")
        new_top = new_left_extra + rc_core + new_right_extra
        new_left_extra_len = len(new_left_extra)
        new_right_extra_len = len(new_right_extra)
    new_left  = {
        "overhang_seq": _rc(right_oh) if right_oh else "",
        "kind":         right_kind,
        "enzyme":       right.get("enzyme", ""),
    }
    new_right = {
        "overhang_seq": _rc(left_oh) if left_oh else "",
        "kind":         left_kind,
        "enzyme":       left.get("enzyme", ""),
    }
    # Feature coords flip relative to the OLD top, then translate
    # into the new top's frame via the strip/extra adjustments.
    new_n = len(new_top)
    _ = new_right_extra_len  # avoid unused-variable lint
    flipped_feats: list[dict] = []
    for f in frag["features"]:
        fs = int(f.get("start", 0))
        fe = int(f.get("end",   0))
        # Map old end → new start, old start → new end. Clamp to
        # the new top's length so out-of-range features (those in
        # the stripped-off old overhang region) collapse to a valid
        # zero-length slice rather than negative coords.
        new_start_raw = (n - fe) - left_strip + new_left_extra_len
        new_end_raw   = (n - fs) - left_strip + new_left_extra_len
        new_f = dict(f)
        new_f["start"]  = max(0, min(new_n, new_start_raw))
        new_f["end"]    = max(0, min(new_n, new_end_raw))
        new_f["strand"] = -int(f.get("strand", 1) or 0) or 0
        flipped_feats.append(new_f)
    return {
        "top_seq":      new_top,
        "left":         new_left,
        "right":        new_right,
        "features":     flipped_feats,
        "source_label": frag["source_label"],
    }


# ═══ Fragment prep for cloning — moved from the hub ═════════════════════════
# Build/excise the fragments the cloning sims consume (synthetic PCR-product
# fragments + digest-and-purify the insert). Enzyme catalog via the hook.

def _make_synthetic_fragment(seq: str, *, enz_left: str, enz_right: str,
                               source_label: str = "",
                               features: "list[dict] | None" = None,
                              ) -> dict:
    """Synthesise a Fragment from a linear DNA sequence by stamping
    canonical sticky ends from the chosen enzymes. The simulator doesn't
    pretend to find recognition sites in the input — it just synthesises the
    overhangs the chosen enzymes would leave.

    NOTE — schematic only. This keeps the WHOLE input as `top_seq` (it does
    NOT cut), so it is used for the Constructor's puzzle-piece T-view (which
    only needs the enzyme overhangs, fed a dummy sequence) and as a low-level
    overhang primitive in the unit tests. The REAL insert excision —
    digesting a PCR product and purifying away the primer pad / off-cut bases
    so they never reach the clone — goes through `_excise_pcr_insert`. Don't
    route a real cloning insert through here: its untrimmed `top_seq` would
    carry the pad + outside-the-cut site bases into the ligation product.

    Both enzymes must resolve via `_state._all_enzymes_hook()` (built-in NEB ∪
    user-added custom enzymes). Raises ValueError on unknown names.
    Features (if any) are carried through verbatim in `top_seq`-local
    coordinates — NOT shifted here."""
    catalog = _state._all_enzymes_hook()
    if enz_left not in catalog:
        raise ValueError(f"unknown enzyme: {enz_left!r}")
    if enz_right not in catalog:
        raise ValueError(f"unknown enzyme: {enz_right!r}")
    site_l, fwd_l, rev_l = catalog[enz_left]
    site_r, fwd_r, rev_r = catalog[enz_right]
    site_l_u = site_l.upper()
    site_r_u = site_r.upper()
    site_l_len = len(site_l_u)
    site_r_len = len(site_r_u)
    # Canonical overhang for `enz_left`: top-strand bases of the
    # recognition site between min(fwd, rev) and max(fwd, rev). This
    # is what the user effectively "buys" by adding the recognition
    # site as a primer tail and then digesting.
    lo_l, hi_l = min(fwd_l, rev_l), max(fwd_l, rev_l)
    lo_r, hi_r = min(fwd_r, rev_r), max(fwd_r, rev_r)
    # Type IIS enzymes (BsaI, Esp3I, BsmBI, etc.) cut OUTSIDE their
    # recognition site, so `hi > site_len` and the overhang depends
    # on bases the user supplies AFTER the recognition — not on
    # the recognition itself. The synthetic-fragment model can't
    # know what those bases are, so it'd produce an empty overhang
    # and silently fail to ligate. Reject upfront with a message
    # pointing at the literal-digest mode (a) which DOES handle
    # Type IIS correctly because it operates on the actual sequence.
    if hi_l > site_l_len or lo_l < 0:
        raise ValueError(
            f"{enz_left!r} cuts outside its recognition site (Type IIS); "
            f"use 'From plasmid' mode for literal digest"
        )
    if hi_r > site_r_len or lo_r < 0:
        raise ValueError(
            f"{enz_right!r} cuts outside its recognition site (Type IIS); "
            f"use 'From plasmid' mode for literal digest"
        )
    overhang_l = site_l_u[lo_l:hi_l]
    overhang_r = site_r_u[lo_r:hi_r]
    kind_l = ("blunt" if fwd_l == rev_l
              else "5'" if fwd_l < rev_l else "3'")
    kind_r = ("blunt" if fwd_r == rev_r
              else "5'" if fwd_r < rev_r else "3'")
    return {
        "top_seq":      seq.upper(),
        "left":  {"overhang_seq": overhang_l, "kind": kind_l,
                  "enzyme": enz_left},
        "right": {"overhang_seq": overhang_r, "kind": kind_r,
                  "enzyme": enz_right},
        "features":     [dict(f) for f in (features or [])],
        "source_label": source_label,
    }


def _enzyme_is_type_iis(name: str) -> bool:
    """True when ``name`` cuts OUTSIDE its recognition site (Type IIS —
    BsaI / BsmBI / BbsI / SapI / Esp3I …). The add-cut-sites cloning model
    only places the recognition site in the primer tail, so a Type IIS
    enzyme (whose overhang depends on bases AFTER the site) can't be used
    that way — mirrors the reject in `_make_synthetic_fragment`."""
    ent = _state._all_enzymes_hook().get(name)
    if not ent or len(ent) < 3:
        return False
    try:
        site = ent[0]
        lo = min(int(ent[1]), int(ent[2]))
        hi = max(int(ent[1]), int(ent[2]))
    except (TypeError, ValueError):
        return False
    return hi > len(site) or lo < 0


@_timed("op.excise_pcr_insert", threshold_ms=25)
def _excise_pcr_insert(seq: str, enz_left: str, enz_right: str, *,
                        features: "list[dict] | None" = None,
                        source_label: str = "",
                        pad: str = "GCGC",
                       ) -> "tuple[dict | None, str | None]":
    """Digest a PCR product and return ONLY the purified insert — the bench
    "cut, then gel-purify away the primer off-cuts" step. Returns
    ``(fragment_or_None, error_or_None)``.

    This is the biologically-faithful replacement for handing a PCR / feature
    insert to `_make_synthetic_fragment`, which kept the WHOLE input as the
    fragment's `top_seq` — so the GCGC primer pad + the recognition-site bases
    OUTSIDE the cut leaked into the ligated clone ("extra bases in the
    product"). At the bench those off-cuts are tiny end pieces you purify away;
    only the middle fragment carries into the ligation.

    Two input shapes, one identical result:

      * **Full PCR product** already carrying the ``[pad][site]…[site][pad]``
        primer tails — what the Clone-region (Alt+Shift+P) workflow builds, or
        a user pasting their real amplicon. A literal linear digest drops the
        flanking pad/off-cut pieces and keeps the single middle fragment.
      * **Bare region** the user amplified (no tails — the Constructor's manual
        "Paste DNA" donor). We wrap it in the SAME ``[pad][site]…[site][pad]``
        the tailed primers would add, then digest — so the purified insert is
        byte-identical to the full-product case.

    Either way NONE of the primer-added pad / outside-the-cut site bases
    survive into the fragment; the retained POST-cut site bases (which reform
    the recognition site when it ligates) DO stay, exactly as at the bench. The
    overhang sequence + kind come straight from the real cut (`_fragments_from_
    cuts`), so 5'/3'/blunt, palindromic and asymmetric cutters all land right
    and the fragment ligates against a literally-digested vector without the
    convention mismatch the synthetic model had.

    Features (amplicon-local coords) are slotted onto the insert by
    `_fragments_from_cuts`; on the wrap path they're first offset past the
    synthesised ``[pad][site]`` lead so they still land on the payload.

    Type IIS is refused (the recognition-only tail can't encode an outside
    cut — mirrors `_make_synthetic_fragment`); use 'From plasmid' literal
    digest for those."""
    catalog = _state._all_enzymes_hook()
    for e in (enz_left, enz_right):
        if e not in catalog:
            return None, f"Unknown enzyme: {e!r}"
        if _enzyme_is_type_iis(e):
            return None, (f"{e} is Type IIS (cuts outside its recognition "
                          "site) — not supported by add-cut-sites cloning; use "
                          "'From plasmid' for a literal digest.")
    cleaned = "".join(c for c in (seq or "").upper()
                       if c in "ACGTRYWSMKBDHVN")
    if not cleaned:
        return None, "Empty PCR-product sequence."

    def _interior(amp: str, feats: "list[dict] | None") -> list[dict]:
        # Literal LINEAR digest, then keep fragments cut on BOTH ends (a real
        # enzyme cut, not the molecule's own linear terminus). The pad/off-cut
        # end pieces always carry one `linear` edge, so they fall away; the
        # purified insert is the piece flanked by two cuts.
        cuts = _enzyme_cuts(amp, [enz_left, enz_right], circular=False)
        frags = _fragments_from_cuts(amp, cuts, circular=False,
                                       features=feats,
                                       source_label=source_label)
        return [f for f in frags
                if f["left"].get("kind") != "linear"
                and f["right"].get("kind") != "linear"]

    def _internal_site_err() -> str:
        return (f"The {enz_left}"
                + (f"/{enz_right}" if enz_right != enz_left else "")
                + " recognition site occurs inside the region you're cloning — "
                "the amplicon would be cut internally and the clone would fail. "
                "Pick an enzyme whose site isn't in your insert.")

    # 1) Treat the input as a FULL PCR product first (tails already present).
    interior = _interior(cleaned, features)
    if len(interior) == 1:
        return interior[0], None
    if len(interior) > 1:
        return None, _internal_site_err()

    # 2) No fragment flanked by two cuts → the input is the BARE region (no
    #    tails). Wrap it the way tailed primers would and digest that.
    site_l = catalog[enz_left][0].upper()
    site_r = catalog[enz_right][0].upper()
    lead = len(pad) + len(site_l)
    wrapped = pad + site_l + cleaned + site_r + _rc(pad)
    w_feats = ([{**f,
                 "start": int(f.get("start", 0)) + lead,
                 "end":   int(f.get("end",   0)) + lead}
                for f in features]
               if features else None)
    interior = _interior(wrapped, w_feats)
    if len(interior) == 1:
        return interior[0], None
    if len(interior) > 1:
        return None, _internal_site_err()
    return None, (f"Couldn't cut the PCR product with {enz_left} / "
                  f"{enz_right} — check the recognition sites.")


@_timed("op.excise_fragment_pair", threshold_ms=25)
def _excise_fragment_pair(seq: str, enzyme_names: list[str], *,
                           circular: bool = True,
                           features: "list[dict] | None" = None,
                           source_label: str = "",
                          ) -> tuple[list[dict], "dict | None"]:
    """Mode (a) helper for both insert and vector: digest `seq` with
    the given enzymes; return ``(fragments, error_dict_or_none)``. If
    the digest produces exactly 2 fragments (typical 1-cut + 1-cut
    case on a circular plasmid), the caller picks one as the insert
    or vector. Otherwise an error dict is returned describing the
    cut count problem.

    The error dict has shape ``{"error": str, "cuts": [{enzyme, top}]}``
    so the UI can surface it clearly."""
    cuts = _enzyme_cuts(seq, enzyme_names, circular=circular)
    n_cuts = len(cuts)
    if n_cuts == 0:
        return [], {"error": (f"no cut sites found for "
                                f"{', '.join(enzyme_names) or '(none)'}"),
                     "cuts": []}
    # Categorise per enzyme so the error message is specific.
    per_enzyme: dict[str, int] = {}
    for c in cuts:
        per_enzyme[c["enzyme"]] = per_enzyme.get(c["enzyme"], 0) + 1
    fragments = _fragments_from_cuts(seq, cuts, circular=circular,
                                       features=features,
                                       source_label=source_label)
    if circular and n_cuts < 2:
        return fragments, {
            "error":
                (f"need ≥2 cuts to excise an insert; "
                 f"got {n_cuts} on a circular plasmid"),
            "cuts": [{"enzyme": c["enzyme"], "top": c["top"]} for c in cuts],
        }
    # Exactly-2 cuts is the only ligation-compatible case for the
    # 2-fragment "insert + vector" model. ≥3 cuts on a circular plasmid
    # produces N fragments and the caller can't unambiguously pick the
    # "insert" vs the "vector" — surface the ambiguity here so future
    # callers can't quietly take fragments[0:2] and ship a wrong product.
    if circular and n_cuts > 2:
        per_enz_str = ", ".join(f"{e}×{n}" for e, n in per_enzyme.items())
        return fragments, {
            "error":
                (f"got {n_cuts} cut sites ({per_enz_str}); "
                 f"need exactly 2 for unambiguous excise. "
                 f"Pick a different enzyme pair."),
            "cuts": [{"enzyme": c["enzyme"], "top": c["top"]} for c in cuts],
        }
    return fragments, None


# ═══ Golden-Braid (BsaI Type IIS) fragment-based scrub — moved from the hub ══
# Cures sites by splitting at each cluster into BsaI-tailed PCR fragments that
# Golden-Gate reassemble into the cured plasmid (real digest+ligate, verified).
# Uses primer's _scrub_design + the now-local _ligate_fragments/_close_circular.

# ── Scrub: Golden Braid (BsaI Type IIS) fragment-based curing ───────────────
#
# An alternative re-circularization to the QuikChange path above. Instead of
# one whole-plasmid amplicon that self-anneals, the plasmid is split at each
# cure cluster into PCR fragments; each fragment's primers carry a BsaI
# (GGTCTC) 5' tail + the NATIVE 4 nt junction overhang, and a Golden Gate /
# Golden Braid (BsaI) reaction reassembles the fragments seamlessly — the
# only net change vs the original is the cured sites.
#
# Because BsaI is the ASSEMBLY enzyme, EVERY BsaI site in the plasmid (not
# just user-listed ones) must be cured — a retained BsaI site is cut
# mid-reaction and scrambles the product. That is the "an unwanted site IS
# BsaI" edge case: BsaI is force-added to the cure set, and a BsaI site that
# can't be silently removed makes the whole Golden Braid design fail (you
# can't Golden-Gate around a site the enzyme still cuts). The cure is sliced
# into the primer ends from the CURED template, so each amplicon is
# internally BsaI-clean; a digest+ligate simulation then proves the
# reassembled product equals the cured plasmid (catastrophic-class: the
# design IS the product). Primer-tail layout mirrors the Domesticator
# (PAD + GGTCTC + spacer + native binding); the binding's first 4 nt ARE the
# fusion overhang, so reassembly restores native sequence except the cures.

_SCRUB_GB_ENZYME    = "BsaI"
_SCRUB_GB_SITE      = "GGTCTC"
_SCRUB_GB_PAD       = "GCGC"        # 5' pad → efficient terminal cutting
_SCRUB_GB_SPACER    = "A"           # 1 nt between recognition and the overhang
_SCRUB_GB_TARGET_TM = 60.0          # NN target for each fragment's annealing arm
_SCRUB_GB_BIND_MIN  = 18
# Grow the binding to reach `_SCRUB_GB_TARGET_TM` while keeping the TOTAL
# scrub oligo (tail + binding) within `_PRIMER_MAX_OLIGO_LEN` — same
# oligo-length budget as every other primer path. 50 − 11 nt tail = 39.
_SCRUB_GB_BIND_MAX  = _PRIMER_MAX_OLIGO_LEN - (
    len(_SCRUB_GB_PAD) + len(_SCRUB_GB_SITE) + len(_SCRUB_GB_SPACER))
_SCRUB_GB_MIN_FRAG  = 50            # shorter than this is impractical to PCR + gel
_SCRUB_GB_OH_SLIDE  = 10            # ± window to slide a junction for a clean overhang
_SCRUB_GB_CLUSTER_FOOTPRINT = 18    # tighter than QuikChange: both arms must reach


def _scrub_gb_overhang_ok(oh: str) -> bool:
    """A 4 nt Golden Gate overhang is usable iff it is exactly 4 ACGT bases,
    NOT palindromic (a palindromic overhang ligates to its own RC → a
    fragment flips or concatemerizes), and not a single base ×4 (AAAA/TTTT/…
    ligate promiscuously and misassemble). These are the standard Golden Gate
    fidelity rules — the difference between a clean 4 nt junction and a
    scrambled assembly."""
    if len(oh) != 4 or any(c not in "ACGT" for c in oh):
        return False
    if oh == _rc(oh):                  # palindrome
        return False
    if len(set(oh)) == 1:              # AAAA / CCCC / GGGG / TTTT
        return False
    return True


def _scrub_gb_cluster_center(positions: list, n: int) -> int:
    """Midpoint coordinate of a cure cluster (wrap-aware), where the junction
    cut is centred so both flanking primers reach every cure in the cluster."""
    start, end = _scrub_cluster_span(positions, n)   # end < start ⇒ wraps origin
    span = (end - start) % n
    return (start + span // 2) % n


def _scrub_gb_pick_cuts(cured: str, clusters: list, n: int) -> tuple:
    """One cut position per cluster, each yielding a clean, mutually-unique
    4 nt Golden Gate overhang `cured[cut:cut+4]`. Searches outward from each
    cluster centre within ±_SCRUB_GB_OH_SLIDE. An overhang AND its reverse
    complement are both reserved once used, so no two junctions can ligate
    cross-wise. Returns `(cuts_sorted, None)` or `(None, reason)` when a
    cluster has no usable overhang window (e.g. an AT-homopolymer junction)."""
    used: set = set()
    raw: list = []
    for c in clusters:
        positions = c["positions"] if isinstance(c, dict) else c
        center = _scrub_gb_cluster_center(positions, n)
        chosen: "int | None" = None
        for d in range(0, _SCRUB_GB_OH_SLIDE + 1):
            for p in dict.fromkeys(((center + d) % n, (center - d) % n)):
                oh = _circ_extract(cured, p, 4, n)
                if not _scrub_gb_overhang_ok(oh):
                    continue
                if oh in used or _rc(oh) in used:
                    continue
                chosen = p
                used.add(oh)
                used.add(_rc(oh))
                break
            if chosen is not None:
                break
        if chosen is None:
            return None, (
                f"no unique non-palindromic BsaI overhang within "
                f"{_SCRUB_GB_OH_SLIDE} bp of the cure near {center} bp — "
                "the QuikChange method has no overhang constraint")
        raw.append(chosen)
    return sorted(set(raw)), None


def _scrub_gb_fragment(cured: str, cut_l: int, cut_r: int, n: int, *,
                       single: bool, idx: int) -> dict:
    """Design ONE Golden Braid fragment + its BsaI-tailed primer pair. The
    fragment spans junction `cut_l`→`cut_r` (the whole plasmid when `single`,
    i.e. a lone junction whose fragment self-circularises). Both primers are
    sliced from the CURED template, so the cure rides in the primer end and
    the binding's first 4 nt are the native fusion overhang."""
    span = n if single else (cut_r - cut_l) % n
    body_len = span + 4                 # include the right overhang's 4 nt
    tail = _SCRUB_GB_PAD + _SCRUB_GB_SITE + _SCRUB_GB_SPACER
    # Forward primer: native binding STARTING at the left junction; its first
    # 4 nt become the 5' overhang after BsaI cuts 1 nt past the recognition.
    fwd_full_bind = _circ_extract(cured, cut_l, _SCRUB_GB_BIND_MAX, n)
    fwd_bind, fwd_tm = _pick_binding_region(
        fwd_full_bind, _SCRUB_GB_TARGET_TM,
        _SCRUB_GB_BIND_MIN, _SCRUB_GB_BIND_MAX)
    fwd_full = tail + fwd_bind
    # Reverse primer: binding ENDS at cut_r+4 (so its last template base is the
    # 4 nt right overhang). rc() of that window is the bottom-strand 5'→3'
    # primer; its first 4 nt = rc(right overhang). Prefix-of-rc trick: growing
    # the picked length extends the 3' anchor leftward while the 5' overhang
    # stays fixed (same as `_pick_binding_region`'s prefix selection).
    rev_win = _circ_extract(cured, (cut_r + 4 - _SCRUB_GB_BIND_MAX) % n,
                            _SCRUB_GB_BIND_MAX, n)
    rev_bind, rev_tm = _pick_binding_region(
        _rc(rev_win), _SCRUB_GB_TARGET_TM,
        _SCRUB_GB_BIND_MIN, _SCRUB_GB_BIND_MAX)
    rev_full = tail + rev_bind
    return {
        "index": idx, "cut_l": cut_l, "cut_r": cut_r, "span": body_len,
        "oh_left": _circ_extract(cured, cut_l, 4, n),
        "oh_right": _circ_extract(cured, cut_r, 4, n),
        "fwd_seq": fwd_full, "rev_seq": rev_full,
        "fwd_bind_len": len(fwd_bind), "rev_bind_len": len(rev_bind),
        "fwd_tm": round(fwd_tm, 1), "rev_tm": round(rev_tm, 1),
        "fwd_gc": round(_mut_gc_pct(fwd_full), 1),
        "rev_gc": round(_mut_gc_pct(rev_full), 1),
    }


def _scrub_gb_build_amplicons(orig: str, cured: str, frags: list, n: int, *,
                              single: bool) -> list:
    """Reconstruct each Golden Braid fragment's REAL PCR product:
    ``tail + body + _rc(tail)`` — ``tail`` is the BsaI cassette
    (``_SCRUB_GB_PAD + _SCRUB_GB_SITE + _SCRUB_GB_SPACER``) and ``body`` carries
    the CURED base inside each primer's binding footprint but the ORIGINAL-
    template base in the PCR-copied middle (a cure outside primer reach would NOT
    make it into the product). Shared by `_scrub_gb_verify` (the design-time
    proof) and the Apply-cure save, which digests + ligates these EXACT amplicons
    so the saved intermediates and the assembled product are the same molecules.
    Returns ``[{index, cut_l, body, amplicon}, …]`` in fragment order.

    NOTE: not `_simulate_pcr` — that is exact-match only and explicitly does not
    handle 5' tailed primers; this is the faithful tailed-mutagenic-PCR product."""
    tail = _SCRUB_GB_PAD + _SCRUB_GB_SITE + _SCRUB_GB_SPACER
    out: list = []
    for fr in frags:
        cut_l, cut_r = fr["cut_l"], fr["cut_r"]
        span = n if single else (cut_r - cut_l) % n
        body_len = span + 4              # include the right overhang's 4 nt
        fwd_len, rev_len = fr["fwd_bind_len"], fr["rev_bind_len"]
        chars: list = []
        for i in range(body_len):
            g = (cut_l + i) % n
            if i < fwd_len or i >= body_len - rev_len:
                chars.append(cured[g])   # inside a primer footprint → cured
            else:
                chars.append(orig[g])    # PCR-copied middle → original template
        body = "".join(chars)
        out.append({"index": fr["index"], "cut_l": cut_l, "body": body,
                    "amplicon": tail + body + _rc(tail)})
    return out


def _scrub_gb_verify(orig: str, cured: str, frags: list, n: int, *,
                     single: bool) -> tuple:
    """Catastrophic-class proof. Build each amplicon AS REAL PCR PRODUCES IT
    (`_scrub_gb_build_amplicons`) then: (1) confirm exactly two BsaI sites per
    amplicon (the two tails, no internal site), (2) release each fragment
    body, (3) ligate by chaining native overhangs, and (4) assert the
    reassembled circle equals the CURED plasmid with no forbidden site left.
    Returns `(ok, errors)`."""
    errors: list = []
    site_pat = _iupac_pattern(_SCRUB_GB_SITE)
    rc_pat = _iupac_pattern(_rc(_SCRUB_GB_SITE))
    amps = _scrub_gb_build_amplicons(orig, cured, frags, n, single=single)
    bodies: list = []
    for a in amps:
        amplicon = a["amplicon"]
        n_sites = len(site_pat.findall(amplicon)) + len(rc_pat.findall(amplicon))
        if n_sites != 2:
            errors.append(
                f"fragment {a['index'] + 1}: amplicon carries {n_sites} BsaI "
                "site(s) (expected exactly 2 — the two tails); an internal "
                "site slipped through and would be cut mid-assembly")
        bodies.append((a["cut_l"], a["body"]))
    if not bodies:
        return False, ["no fragments to verify"]
    bodies.sort(key=lambda t: t[0])
    # Each body ends with the next fragment's left overhang (4 nt shared at the
    # junction); drop it so each base appears once in the reassembled circle.
    product = "".join(b[:-4] for _, b in bodies)
    expect = _circ_extract(cured, bodies[0][0], n, n)
    if len(product) != n or product != expect:
        errors.append(
            "re-assembled product does not equal the cured plasmid — a cure "
            "fell outside primer reach of its junction, or a junction is "
            "mis-placed")
    # Final whole-circle scan: the product must be BsaI-CLEAN (wrap-aware).
    # Every BsaI site was force-cured and the assembly cuts off the added
    # tails, so a residual BsaI would mean the reaction re-cuts the product
    # (fratricide). NON-assembly unwanted sites that couldn't be cured stay
    # reported in `sites_skipped` (exactly like QuikChange) — they leave the
    # product unchanged at that site but must NOT fail the whole design.
    asm = (_SCRUB_GB_SITE, _rc(_SCRUB_GB_SITE))
    if _forbidden_hit_set(product + product[:len(_SCRUB_GB_SITE) - 1], asm):
        errors.append("re-assembled product still carries a BsaI site")
    return (not errors), errors


def _assemble_scrub_amplicons_real(amplicon_specs: list, *,
                                   enzyme: str = _SCRUB_GB_ENZYME
                                   ) -> "str | None":
    """The REAL Golden Braid assembly behind Apply cure (vs the design-time
    proof in `_scrub_gb_verify`, which string-chains bodies): genuinely
    BsaI-digest each amplicon (`_digest_with_enzymes`), keep the released insert
    body (the ONLY fragment with a sticky overhang on BOTH edges — each end stub
    has one ``linear`` edge), chain the inserts by matching native 4 nt overhangs
    (`_ligate_fragments`), and close the circle (`_close_circular`) — the same
    primitives the Constructor's Golden Gate assembly uses.

    `amplicon_specs` is ``[{cut_l, amplicon}, …]``; ``cut_l`` orders the inserts
    around the plasmid. Returns the circular product's top-strand sequence, or
    None if any digest / overhang-mismatch / closure step fails — the caller then
    ABORTS the save rather than persist a mis-assembled plasmid."""
    inserts: list = []
    for spec in sorted(amplicon_specs, key=lambda s: s["cut_l"]):
        try:
            frags = _digest_with_enzymes(spec["amplicon"], [enzyme],
                                         circular=False)
        except Exception:
            _log.exception("scrub assemble: BsaI digest failed")
            return None
        body = [f for f in frags
                if f.get("left", {}).get("kind") != "linear"
                and f.get("right", {}).get("kind") != "linear"]
        if len(body) != 1:
            _log.error("scrub assemble: amplicon released %d insert(s) "
                       "(expected exactly 1 between the two BsaI cuts)",
                       len(body))
            return None
        inserts.append(body[0])
    if not inserts:
        return None
    chained = inserts[0]
    for nxt in inserts[1:]:
        chained = _ligate_fragments(chained, nxt)
        if chained is None:
            _log.error("scrub assemble: overhang mismatch chaining fragments")
            return None
    closed = _close_circular(chained)
    if closed is None:
        _log.error("scrub assemble: terminal overhangs don't close the circle")
        return None
    return str(closed.get("top_seq") or "") or None


# ═══ Golden Gate / MoClo (Type IIS) assembly — SC-H ════════════════════════════
# Overhang-directed N-part assembly: digest every part + the destination vector
# with the Type IIS enzyme (BsaI / BsmBI / BbsI / SapI / Esp3I), keep each
# released body (cut on BOTH ends), then chain them by their native 4 nt
# overhangs into a circle. Unlike the scrub assembler above the parts are NOT
# pre-ordered — the overhangs alone dictate the order (the whole point of Golden
# Gate). [INV-127: the design IS the product — a real digest + ligation.]


def _gg_released_bodies(seq: str, enzyme: str, *, circular: bool) -> list[dict]:
    """Digest `seq` with `enzyme` and return only the fragments cut on BOTH
    ends — the released part body / vector backbone with two sticky overhangs
    (the off-cut end stubs each keep one ``linear`` molecule-terminus edge and
    fall away). Mirrors `_assemble_scrub_amplicons_real`'s body selection."""
    frags = _digest_with_enzymes(seq, [enzyme], circular=circular)
    return [f for f in frags
            if f.get("left", {}).get("kind") != "linear"
            and f.get("right", {}).get("kind") != "linear"]


def _gg_greedy_chain(start: dict, parts: list[dict]) -> "list[dict] | None":
    """Order `parts` after `start` by overhang matching: repeatedly append the
    part whose LEFT end ligates the chain's current RIGHT end (`_ligate_fragments`
    returns non-None). Returns the ordered chain using EVERY part, or None if it
    ever gets stuck. A valid Golden Gate design has unique overhangs, so exactly
    one part matches at each step (deterministic); a tie/no-match means a
    mis-designed overhang set."""
    chain = [start]
    pool = list(parts)
    while pool:
        cur = chain[-1]
        idx = next((i for i, f in enumerate(pool)
                    if _ligate_fragments(cur, f) is not None), None)
        if idx is None:
            return None
        chain.append(pool.pop(idx))
    return chain


def _simulate_golden_gate(part_seqs: list[str], vector_seq: str, *,
                          enzyme: str = _SCRUB_GB_ENZYME) -> dict:
    """Simulate a Golden Gate / MoClo (Type IIS) one-pot assembly.

    Digests each part in `part_seqs` and the destination `vector_seq` with
    `enzyme`, then chains the released fragments by matching 4 nt overhangs
    into a circular product. The parts may be supplied in ANY order — the
    overhangs determine the assembly order. Returns:

    ``{ok, product_seq, length, circular, n_parts, enzyme, order,
       junctions:[{overhang}], n_residual_sites, warnings, errors}``

    `ok:false` (with `errors`) when the enzyme isn't Type IIS, a part doesn't
    release exactly one fragment, or the overhangs don't chain into a closed
    circle. Fidelity `warnings` flag a non-unique junction overhang (ambiguous
    assembly) or a residual enzyme site in the product (it would be re-cut)."""
    warnings: list[str] = []
    if not _enzyme_is_type_iis(enzyme):
        return {"ok": False, "errors": [
            f"{enzyme!r} is not a Type IIS enzyme — Golden Gate / MoClo needs "
            f"one that cuts OUTSIDE its recognition site (BsaI / BsmBI / BbsI "
            f"/ SapI / Esp3I)."], "warnings": []}
    if not part_seqs:
        return {"ok": False, "errors": ["no parts supplied"], "warnings": []}

    part_frags: list[dict] = []
    for i, ps in enumerate(part_seqs):
        try:
            bodies = _gg_released_bodies((ps or "").upper(), enzyme,
                                         circular=False)
        except Exception as exc:               # pragma: no cover - defensive
            return {"ok": False, "warnings": [],
                    "errors": [f"part {i + 1}: digest failed: {exc}"]}
        if len(bodies) != 1:
            return {"ok": False, "warnings": [], "errors": [
                f"part {i + 1} released {len(bodies)} fragment(s) when cut "
                f"with {enzyme} (expected exactly 1 between two {enzyme} "
                f"sites — check the part's flanking sites / orientation)."]}
        part_frags.append(bodies[0])

    try:
        vec_bodies = _gg_released_bodies((vector_seq or "").upper(), enzyme,
                                         circular=True)
    except Exception as exc:                   # pragma: no cover - defensive
        return {"ok": False, "warnings": [],
                "errors": [f"vector: digest failed: {exc}"]}
    if not vec_bodies:
        return {"ok": False, "warnings": [], "errors": [
            f"vector released no fragment when cut with {enzyme} — it needs "
            f"two {enzyme} sites flanking the dropout."]}

    # Try each vector fragment as the backbone seed (the dropout/stuffer
    # candidate won't chain — its overhangs don't match the parts).
    chosen: "tuple[list[dict], dict] | None" = None
    for vstart in vec_bodies:
        chain = _gg_greedy_chain(vstart, part_frags)
        if chain is None:
            continue
        ligated = chain[0]
        good = True
        for f in chain[1:]:
            ligated = _ligate_fragments(ligated, f)
            if ligated is None:
                good = False
                break
        if not good:
            continue
        assert ligated is not None         # `good` ⇒ no ligation returned None
        closed = _close_circular(ligated)
        if closed is not None:
            chosen = (chain, closed)
            break
    if chosen is None:
        return {"ok": False, "warnings": [], "errors": [
            "the overhangs don't chain every part + the vector into a closed "
            "circle — check that adjacent parts share a 4 nt overhang and the "
            "vector's two overhangs match the assembly's ends."]}

    chain, closed = chosen
    product = str(closed.get("top_seq") or "")
    # Fidelity: every junction overhang must be DISTINCT, or the reaction can
    # mis-assemble (two junctions with the same overhang are interchangeable).
    junctions = [f.get("right", {}).get("overhang_seq", "") for f in chain]
    if len(set(junctions)) != len(junctions):
        warnings.append(
            "non-unique junction overhang(s) — the assembly is ambiguous and "
            "can mis-assemble; redesign to all-distinct 4 nt overhangs.")
    # Residual enzyme sites in the product would be re-cut in the one-pot rxn.
    try:
        residual = _enzyme_cuts(product, [enzyme], circular=True)
    except Exception:                          # pragma: no cover - defensive
        residual = []
    if residual:
        warnings.append(
            f"{len(residual)} residual {enzyme} site(s) in the product — it "
            f"would be re-cut during the one-pot reaction; domesticate them "
            f"out of the parts first.")
    return {
        "ok":              True,
        "product_seq":     product,
        "length":          len(product),
        "circular":        True,
        "n_parts":         len(part_frags),
        "enzyme":          enzyme,
        "order":           [f.get("source_label") or f"frag{i}"
                            for i, f in enumerate(chain)],
        "junctions":       [{"overhang": j} for j in junctions],
        "n_residual_sites": len(residual),
        "warnings":        warnings,
        "errors":          [],
    }


def _scrub_gb_design(seq: str, feats: "list | None" = None, enzymes=None, *,
                     circular: bool = True, codon_raw: "dict | None" = None
                     ) -> dict:
    """Plan a Golden Braid (BsaI Type IIS) fragment cure of `seq`. Cures the
    selected unwanted sites AND every BsaI site (the assembly enzyme), splits
    the plasmid at each cure cluster into PCR fragments with BsaI-tailed
    primers carrying the native 4 nt junction overhangs, and verifies a
    digest+ligate reassembles the cured plasmid seamlessly. Pure biology —
    designs the primers but mutates nothing on disk.

    `codon_raw` biases coding cures toward host-frequent synonyms (tie-break
    only — synonymy is guaranteed independently). Returns a result dict with
    `ok`, `cured_seq`, `fragments` (each with its primer pair + overhangs),
    `sites_removed` / `sites_skipped`, `verified`, `warnings`, and `errors`."""
    seq = (seq or "").upper()
    n = len(seq)
    # The assembly enzyme MUST be in the cure set — a retained BsaI site is cut
    # during the Golden Gate reaction and scrambles the product. Force-add it
    # so BsaI is always cured, never merely tolerated ("unwanted site IS BsaI").
    base = list(enzymes) if enzymes is not None else list(_SCRUB_DEFAULT_ENZYMES)
    if _SCRUB_GB_ENZYME not in base:
        base.append(_SCRUB_GB_ENZYME)
    plan = _scrub_design(seq, feats, base, circular=circular,
                         codon_raw=codon_raw)
    result: dict = {
        "ok": True, "method": "golden_braid", "enzyme": _SCRUB_GB_ENZYME,
        "orig_seq": seq, "cured_seq": plan.get("cured_seq", seq),
        "edits": list(plan.get("edits", [])),
        "sites_removed": list(plan.get("sites_removed", [])),
        "sites_skipped": list(plan.get("sites_skipped", [])),
        "fragments": [], "n_fragments": 0, "verified": False,
        "warnings": list(plan.get("warnings", [])), "errors": [],
    }
    if not seq:
        result["ok"] = False
        result["errors"].append("No sequence loaded.")
        return result
    if not plan.get("ok", False):
        result["ok"] = False
        result["errors"].append("Curing failed; nothing to assemble.")
        return result
    # Edge case — a BsaI site that can't be silently cured is FATAL for Golden
    # Braid: you cannot Golden-Gate around a site the assembly enzyme cuts.
    bsai_skipped = [s for s in result["sites_skipped"]
                    if s.get("enzyme") == _SCRUB_GB_ENZYME]
    if bsai_skipped:
        result["ok"] = False
        spots = ", ".join(str(s.get("pos")) for s in bsai_skipped)
        result["errors"].append(
            f"{len(bsai_skipped)} BsaI site(s) (at {spots}) can't be silently "
            "removed, but BsaI is the assembly enzyme — it would cut the "
            "fragments mid-reaction. Golden Braid curing is impossible here; "
            "use QuikChange, or supply a codon table so the coding site has a "
            "synonymous alternative.")
        return result
    cured = result["cured_seq"]
    if not result["edits"]:
        result["warnings"].append(
            "No sites needed curing — nothing to fragment.")
        return result
    # Junctions sit at cure clusters (tighter footprint than QuikChange so both
    # flanking primer arms reach every cure in the cluster).
    gb_clusters = _scrub_cluster_edits(
        [e["pos"] for e in result["edits"]], n,
        footprint=_SCRUB_GB_CLUSTER_FOOTPRINT)
    cuts, reason = _scrub_gb_pick_cuts(cured, gb_clusters, n)
    if not cuts:                       # None (no clean overhang) or empty
        result["ok"] = False
        result["errors"].append(reason or "could not place any junction")
        return result
    single = (len(cuts) == 1)
    frags: list = []
    for k in range(len(cuts)):
        cut_l = cuts[k]
        cut_r = cuts[0] if single else cuts[(k + 1) % len(cuts)]
        frags.append(_scrub_gb_fragment(cured, cut_l, cut_r, n,
                                        single=single, idx=k))
    short = [fr for fr in frags if fr["span"] < _SCRUB_GB_MIN_FRAG]
    if short:
        result["ok"] = False
        result["errors"].append(
            f"{len(short)} fragment(s) shorter than {_SCRUB_GB_MIN_FRAG} bp — "
            "the junctions are too close to PCR + gel-purify reliably; use "
            "QuikChange or cure fewer sites.")
        result["fragments"] = frags
        result["n_fragments"] = len(frags)
        return result
    ok, errors = _scrub_gb_verify(seq, cured, frags, n, single=single)
    result["fragments"] = frags
    result["n_fragments"] = len(frags)
    result["verified"] = ok
    if not ok:
        result["ok"] = False
        result["errors"].extend(errors)
    return result


# ═══ Golden-Braid / MoClo domestication primer designers — moved from the hub ═
# _design_gb_primers (single-part domestication) + _design_operon_soe_primers
# (whole-operon SOE) + the grammar-position / forbidden-site / pair-name helpers.
# Reuse primer L2 (designers + _scrub_design) + codon L2 (_codon_fix_*) + the
# GB grammar data (dataaccess) + _atg_offset_for_part (local).

# Backbone-marker detection for the dropout picker in the resolver below.
# cloning (L3) can't import splicecraft_seqanalysis._fragment_has_backbone_marker
# (L3→L3 cycle: seqanalysis imports cloning's _excise_fragment_pair), so a
# compact copy of its constants lives here — keep in sync with
# splicecraft_seqanalysis._BACKBONE_{FEATURE_TYPES,LABEL_KEYWORDS}.
_ACCEPTOR_BACKBONE_FEATURE_TYPES = frozenset({"rep_origin", "oriT"})
_ACCEPTOR_BACKBONE_LABEL_KEYWORDS = (
    "ori", "rep_origin", "ampr", "kanr", "specr", "specinomycin",
    "spectinomycin", "cmr", "chloramphenicol", "tetr", "tetracyclin",
    "carbr", "carbenicillin",
    # bare "selection"/"antibiotic" excluded on purpose — see the rationale on
    # splicecraft_seqanalysis._BACKBONE_LABEL_KEYWORDS (a counter-selection
    # dropout labeled "...selection..." must NOT read as backbone).
)
# Whole-token keywords — mirrors seqanalysis `_BACKBONE_TOKEN_KEYWORDS`. These
# are too short to substring-match ("erm" is inside every "Terminator"), and
# are consulted ONLY on gene-bearing feature types — a `regulatory` or
# `primer_bind` label mentioning a marker names something NEAR the gene, not
# the gene. See the seqanalysis constant for the worked examples.
_ACCEPTOR_BACKBONE_TOKEN_FEATURE_TYPES = frozenset({"cds", "gene"})
_ACCEPTOR_BACKBONE_TOKEN_KEYWORDS = frozenset({
    "erm", "ermb", "ermc", "ermam", "erma", "ermg",
    "bla", "cat", "aada", "aac", "aph",
    "neo", "neor", "nptii", "kan",
    "smr", "hyg", "hygr", "hph",
    "zeo", "zeor", "ble", "gmr", "gent",
    "tetm", "teta", "tetl",
    "puror", "pac", "sh", "nat1",
})


def _frag_carries_backbone_marker(frag: dict) -> bool:
    """True iff ``frag``'s features include an origin / resistance marker — the
    signal that identifies the vector BACKBONE (never the dropout an insert
    replaces). Mirrors splicecraft_seqanalysis._fragment_has_backbone_marker."""
    for f in (frag.get("features") or []):
        if not isinstance(f, dict):
            continue
        ftype = str(f.get("type") or "").lower()
        if ftype in _ACCEPTOR_BACKBONE_FEATURE_TYPES:
            return True
        label = str(f.get("label") or "").lower()
        if not label:
            continue
        if any(kw in label for kw in _ACCEPTOR_BACKBONE_LABEL_KEYWORDS):
            return True
        if ftype in _ACCEPTOR_BACKBONE_TOKEN_FEATURE_TYPES:
            tokens = {t for t in re.split(r"[^a-z0-9]+", label) if t}
            if tokens & _ACCEPTOR_BACKBONE_TOKEN_KEYWORDS:
                return True
    return False


def _acceptor_vector_features(rec) -> "list[dict]":
    """Minimal feature dicts (type/label/start/end) from a parsed entry-vector
    record, so the digest can tag which fragment carries a backbone marker.
    Skips features whose coords don't resolve; wrap/compound features flatten to
    their span (fine — only used to identify the backbone fragment)."""
    out: list[dict] = []
    for f in getattr(rec, "features", None) or []:
        try:
            loc = getattr(f, "location", None)
            s, e = int(loc.start), int(loc.end)   # type: ignore[union-attr]
        except Exception:
            continue
        if e <= s:
            continue
        try:
            q = getattr(f, "qualifiers", {}) or {}
            lbl = (q.get("label") or q.get("gene") or q.get("product")
                   or [""])[0]
        except Exception:
            lbl = ""
        out.append({"type": str(getattr(f, "type", "") or ""),
                    "label": str(lbl or ""), "start": s, "end": e})
    return out


def _entry_vector_acceptor_overhangs(
    grammar_id: str, role: str = "", enzyme: str = "",
) -> "tuple[str, str] | None":
    """The 4-nt overhang pair a domestication fragment must PRESENT to drop
    into the CONFIGURED entry vector for ``(grammar_id, role)``.

    Digests that vector with ``enzyme`` (the grammar's Type IIS enzyme) and
    returns the dropout/stuffer fragment's ``(oh5, oh3)`` — the acceptor
    overhangs the incoming part inherits. For a two-tier Golden-Braid UPD
    vector these are the EXTERNAL entry overhangs (e.g. ``CTCG`` / ``TGAG``),
    which differ from the grammar's INTERNAL category overhangs (``AATG`` /
    ``GCTT`` for a CDS): the fragment enters the vector on the external pair
    (Esp3I/BsmBI) and later releases the L0 part on the internal pair (BsaI),
    which the vector's own nested sites regenerate. `_design_gb_primers` nests
    the category overhangs inside this pair when they differ.

    Returns ``None`` — and the caller keeps the classic one-tier layout
    (category overhang IS the entry overhang, standard MoClo pUPD2) — when no
    vector is configured, it carries no sequence, or it doesn't cut to a clean
    2-fragment acceptor. Never raises: a misconfigured vector degrades to the
    safe one-tier default rather than aborting the design."""
    if not enzyme:
        return None
    try:
        ev = _get_entry_vector(grammar_id, role or "")
    except Exception:
        _log.exception("_entry_vector_acceptor_overhangs: EV lookup failed")
        return None
    if not isinstance(ev, dict):
        return None
    gb = ev.get("gb_text") or ""
    if not gb:
        return None
    try:
        rec = _gb_text_to_record(gb)
        seq = str(getattr(rec, "seq", "") or "").upper()
        if not seq:
            return None
        frags, err = _excise_fragment_pair(
            seq, [enzyme], circular=True,
            features=_acceptor_vector_features(rec))
    except Exception:
        _log.exception("_entry_vector_acceptor_overhangs: EV digest failed")
        return None
    if err is not None or len(frags) != 2:
        return None
    # Identify the DROPOUT (the fragment an insert replaces) by the ABSENCE of a
    # backbone marker — NOT by size, which inverts the overhang order the moment
    # a visible-marker / counterselection dropout outgrows the backbone (never
    # assume the smaller frag is the payload). Fall back to smallest only when
    # the marker signal is ambiguous (both or neither fragment annotated).
    non_backbone = [f for f in frags if not _frag_carries_backbone_marker(f)]
    stuffer = (non_backbone[0] if len(non_backbone) == 1
               else min(frags, key=lambda f: len(f.get("top_seq", ""))))
    oh5 = ((stuffer.get("left") or {}).get("overhang_seq") or "").upper()
    oh3 = ((stuffer.get("right") or {}).get("overhang_seq") or "").upper()
    if len(oh5) == 4 and len(oh3) == 4:
        return (oh5, oh3)
    return None


def _dom_primer_pair_names(part_name: str, idx: int = 1) -> "tuple[str, str]":
    """The canonical ``(fwd, rev)`` primer NAMES for a domestication pair:
    ``{part_name}-DOM-{idx}-F`` / ``-R``. Single source of the naming
    convention shared by `_save_primers_to_library` (the primer LIBRARY
    entry), the `primer_bind` feature LABELS on the fragment + clone (via
    `_part_primer_labels`), and the construction-history ``<Oligo>`` names
    — so every surface shows the SAME name for a given primer (user
    request 2026-06-05: "all unified"). ``idx`` is the 1-indexed pair
    number within a run, matching `_save_primers_to_library`'s loop."""
    base = str(part_name or "").strip() or "part"
    return f"{base}-DOM-{idx}-F", f"{base}-DOM-{idx}-R"


def _grammar_position_by_type(grammar: dict, ptype: str) -> "dict | None":
    """Helper: find the position spec for a given part type within a
    grammar. ``None`` if the grammar doesn't define that type — which
    e.g. means CDS-NS isn't a valid pick under MoClo Plant."""
    for pos in grammar.get("positions", []):
        if pos.get("type") == ptype:
            return pos
    return None


def _gb_find_forbidden_hits(
    seq: str,
    sites: "dict[str, str] | None" = None,
) -> list[tuple[str, str, int]]:
    """Return ``(enzyme_name, site_found, position)`` for **every** internal
    Type IIS recognition in *seq*, on both forward and reverse strands.

    The ``sites`` map (``{enzyme_name: recognition}``) defaults to the
    Golden Braid L0 forbidden set (Esp3I + BsaI). Pass a different
    grammar's ``forbidden_sites`` to scan against MoClo (BsaI + BpiI),
    a custom grammar, or any other Type IIS combination. Returns every
    occurrence — not just the first per enzyme. Critical for accurate
    reporting when an insert contains multiple sites: the user must
    know about all of them before paying for a gBlock synthesis.
    Results are sorted by position to aid downstream reporting.
    """
    if sites is None:
        sites = _GB_DOMESTICATION_FORBIDDEN
    out: list[tuple[str, str, int]] = []
    for name, site in sites.items():
        if not isinstance(site, str) or not site:
            continue
        rc = _rc(site)
        needles = [site] if rc == site else [site, rc]
        for needle in needles:
            start = 0
            while True:
                pos = seq.find(needle, start)
                if pos == -1:
                    break
                out.append((name, needle, pos))
                start = pos + 1
    out.sort(key=lambda t: (t[2], t[0], t[1]))
    return out


def _gb_binding_region_advisory(
    mutations: list[str],
    insert_len: int,
    fwd_bind_len: int,
    rev_bind_len: int,
    fwd_skip: int = 0,
) -> list[dict]:
    """Return one entry per mutation that lands inside a primer binding
    window. Each entry is ``{"text", "region", "codon_start"}`` where
    ``region`` is ``"fwd"`` or ``"rev"``. Empty list when every mutation
    is safely inside the amplicon's interior.

    Why this matters: if a mutation falls inside the 5′ or 3′ binding
    window, the PCR primer won't bind perfectly to the user's original
    plasmid template — they must order the *mutated* insert as a gBlock
    and use that as the PCR template (or redesign around the site).
    """
    if not mutations or insert_len <= 0:
        return []
    positions = _codon_fix_mutation_positions(mutations)
    # The forward primer binds `insert[fwd_skip : fwd_skip+fwd_bind_len]` — an
    # ATG-fusion part (GB L0 `AATG`) sets fwd_skip=3, so the real binding does
    # NOT start at 0. Anchoring the window at 0 (the old behaviour) both
    # over-flagged mutations in the skipped overhang region AND missed the one
    # at the 3′ terminus of the binding (codon in [fwd_bind_len, fwd_bind_len+3)),
    # silently omitting the "order this as a gBlock" warning for the common
    # Golden-Braid L0 CDS path.
    fwd_lo = max(0, fwd_skip)                    # [fwd_lo, fwd_hi) covers fwd binding
    fwd_hi = fwd_lo + max(0, fwd_bind_len)
    rev_lo = insert_len - max(0, rev_bind_len)  # [rev_lo, insert_len) covers rev binding
    out: list[dict] = []
    for text, codon_start in zip(mutations, positions):
        if codon_start < 0:
            continue
        codon_end = codon_start + 3
        # A 3-nt codon overlaps the fwd window if it intersects [fwd_lo, fwd_hi).
        in_fwd = codon_start < fwd_hi and codon_end > fwd_lo
        # Overlaps the rev window if any nt is in [rev_lo, insert_len).
        in_rev = codon_end > rev_lo
        if in_fwd:
            out.append({"text": text, "region": "fwd",
                        "codon_start": codon_start})
        if in_rev:
            out.append({"text": text, "region": "rev",
                        "codon_start": codon_start})
    return out


@_timed("op.primer3.gb_design")
def _design_gb_primers(
    template_seq: str,
    start: int,
    end: int,
    part_type: str,
    target_tm: float = 60.0,
    codon_raw: "dict | None" = None,
    grammar: "dict | None" = None,
    entry_overhangs: "tuple[str, str] | None" = None,
) -> dict:
    """Design modular cloning domestication primers for a template region.

    Defaults to Golden Braid L0 (Esp3I, GGAG/TGAC/AATG/GCTT/CGCT
    overhangs); pass ``grammar`` to use a different cloning grammar
    (MoClo Plant, custom user-defined). The amplified product, after
    digestion with the grammar's enzyme, carries the 4-nt overhangs
    associated with ``part_type`` in that grammar's position table.

    Primer structure (5'→3'):

        Forward: [pad] [enzyme site] [spacer] [5' overhang]    [binding →]
        Reverse: [pad] [enzyme site] [spacer] [RC 3' overhang] [← binding RC]

    Pass ``entry_overhangs`` = ``(external_oh5, external_oh3)`` to nest a
    two-tier Golden-Braid layout: the grammar's category overhangs (``oh5`` /
    ``oh3``) move INSIDE the given external pair, so the Type IIS cut presents
    the external overhangs (which enter a UPD/pUPD entry vector) while the
    category overhangs survive to be released one level up. When ``None`` — or
    equal to the category pair — the classic one-tier layout is emitted
    unchanged (external overhang IS the category overhang, standard MoClo).
    Resolve the pair from the configured entry vector with
    :func:`_entry_vector_acceptor_overhangs`. The result reports the fragment's
    actual cut overhangs under ``entry_oh5`` / ``entry_oh3``.

    When *codon_raw* (a ``{codon: (aa, count)}`` dict from the codon-table
    registry) is supplied and *part_type* is a coding type (CDS / CDS-NS /
    C-tag), internal BsaI or Esp3I sites are silently repaired by
    substituting synonymous codons before the primers are designed — the
    returned ``insert_seq`` is then the *mutated* sequence (which is what
    the user should order as a gBlock / synthetic fragment for PCR).
    The list of substitutions made is returned under ``mutations``.

    Returns a dict with keys: part_type, position, oh5, oh3, insert_seq,
    fwd_binding, rev_binding, fwd_full, rev_full, fwd_tm, rev_tm,
    amplicon_len, mutations, and a ``pairs`` list. ``pairs`` holds one
    dict per amplicon — each with ``fwd_full``, ``rev_full``, ``fwd_tm``,
    ``rev_tm``, ``fwd_pos``, ``rev_pos``, ``fwd_binding``, ``rev_binding``,
    ``amplicon_len``. Callers that only need the first pair can continue
    to read the legacy top-level keys; multi-pair savers should iterate
    ``pairs`` (future SOE-PCR splitting for non-repairable internal sites
    will extend this list beyond one pair).
    """
    g = grammar if isinstance(grammar, dict) else _BUILTIN_GRAMMARS["gb_l0"]
    pos_spec = _grammar_position_by_type(g, part_type)
    if pos_spec is None:
        return {
            "error": f"Part type {part_type!r} is not defined in grammar "
                     f"{g.get('name', '?')}. Available types: "
                     f"{', '.join(p.get('type', '?') for p in g.get('positions', []))}.",
            "mutations": [],
        }
    pos_label, oh5, oh3 = pos_spec.get("name", "?"), pos_spec.get("oh5", ""), pos_spec.get("oh3", "")
    forbidden_sites = g.get("forbidden_sites", _GB_DOMESTICATION_FORBIDDEN) or _GB_DOMESTICATION_FORBIDDEN
    coding_types = set(g.get("coding_types", []) or _GB_CODING_PART_TYPES)
    enzyme_pad = g.get("pad", _GB_PAD)
    enzyme_site = g.get("site", _GB_L0_ENZYME_SITE)
    enzyme_spacer = g.get("spacer", _GB_SPACER)

    # Two-tier Golden-Braid nesting. `oh5`/`oh3` are the grammar's INTERNAL
    # category overhangs (AATG/GCTT for a CDS — what BsaI releases at L1). When
    # the destination entry vector exposes DIFFERENT external acceptor overhangs
    # (`entry_overhangs`, e.g. CTCG/TGAG for a UPD vector — resolved by the
    # caller via `_entry_vector_acceptor_overhangs`), the Type IIS cut must
    # present THOSE, with the category pair nested one layer in. `ent5`/`ent3`
    # are that external layer (empty ⇒ classic one-tier: category == entry, the
    # standard MoClo pUPD2 layout, byte-identical to the pre-nesting design).
    # `entry_oh5`/`entry_oh3` report the fragment's ACTUAL enzyme-cut overhangs
    # (what physically ligates into the acceptor) — the entry-vector check reads
    # these, not the category pair.
    ent5 = ent3 = ""
    entry_oh5, entry_oh3 = oh5, oh3
    if isinstance(entry_overhangs, (tuple, list)) and len(entry_overhangs) == 2:
        e5 = str(entry_overhangs[0] or "").upper()
        e3 = str(entry_overhangs[1] or "").upper()
        if len(e5) == 4 and len(e3) == 4 and (e5, e3) != (oh5.upper(), oh3.upper()):
            ent5, ent3 = e5, e3
            entry_oh5, entry_oh3 = e5, e3

    total  = len(template_seq)
    insert = _slice_circular(template_seq.upper(), start, end)
    wraps  = end < start

    # Need at least 18 bp to pick a proper binding region — otherwise
    # _pick_binding_region returns the whole (too-short) insert with Tm=0.
    if len(insert) < 18:
        return {
            "error": f"Cloning region is too short ({len(insert)} bp). "
                     f"Select at least 18 bp (recommended 25+ bp for a "
                     f"robust binding region).",
            "mutations": [],
        }

    # Internal Type IIS check. The grammar's `forbidden_sites` lists
    # every recognition that must be absent from the final part — for
    # GB L0 that's Esp3I (current cut) + BsaI (next-level reuse); for
    # MoClo Plant, BsaI + BpiI. Coding parts can be repaired via
    # synonymous codon substitution; non-coding parts have no reading
    # frame so internal sites must be fixed manually.
    mutations: list[str] = []
    initial_hits = _gb_find_forbidden_hits(insert, sites=forbidden_sites)
    if initial_hits:
        hit_str = ", ".join(
            f"{name} {site} at +{pos + 1}"
            for name, site, pos in initial_hits
        )
        can_attempt_fix = (
            part_type in coding_types
            and bool(codon_raw)
            and len(insert) % 3 == 0
        )
        if can_attempt_fix and codon_raw:
            # Inline `codon_raw` narrowing for pyright: the `bool(codon_raw)`
            # inside `can_attempt_fix` is not visible across the assignment,
            # so re-asserting truthiness here lets the type checker see
            # `codon_raw` as non-None at the call site below.
            protein = _mut_translate(insert)
            # 2026-05-27 (audit-5 H3): validate the reading frame
            # before codon repair. `_mut_translate` stops on the
            # first stop codon, so if the user's CDS selection is
            # off-frame (codon_start != 1, partial CDS, or an
            # off-by-one in the selected region) the protein comes
            # out short and `_codon_fix_sites` will silently
            # substitute synonyms in the WRONG codon table for every
            # position — silently corrupting the synthesized part.
            # Require the translated protein to cover ≥ 90 % of the
            # insert's codons (allows for a single mid-CDS stop in
            # error-tolerant cases but catches off-by-1 / off-by-2).
            expected_codons = len(insert) // 3
            if expected_codons and len(protein) < int(expected_codons * 0.9):
                return {
                    "error": (
                        f"CDS reading-frame validation failed: translated "
                        f"protein is {len(protein)} aa but the {len(insert)} "
                        f"bp insert should encode ~{expected_codons} aa. "
                        f"The selection is likely off-frame "
                        f"(check `codon_start` qualifier, partial CDS, or "
                        f"adjust selection boundaries to align with codon 1)."
                    ),
                }
            if protein:
                # 2026-05-27 (audit-5 H2): the insert is the user's
                # raw CDS region — NOT a `_codon_optimize` output
                # with an appended stop. Pass `has_appended_stop=False`
                # so the boundary check doesn't silently skip the
                # last 1-2 codons (leaving forbidden sites overlapping
                # the C-terminus unfixed).
                fixed_insert, mutations = _codon_fix_sites(
                    insert, protein, codon_raw,
                    sites=forbidden_sites,
                    has_appended_stop=False,
                )
                remaining = _gb_find_forbidden_hits(
                    fixed_insert, sites=forbidden_sites,
                )
                if remaining:
                    remain_str = ", ".join(
                        f"{name} {site} at +{pos + 1}"
                        for name, site, pos in remaining
                    )
                    return {
                        "error": f"Internal Type IIS site(s) remain after "
                                 f"silent-mutation attempt ({remain_str}). "
                                 f"The sites overlap codons with no "
                                 f"synonymous alternative in this codon "
                                 f"table — pick a different region or "
                                 f"redesign.",
                        "mutations": mutations,
                    }
                insert = fixed_insert
            else:
                return {
                    "error": f"Internal Type IIS site(s) found ({hit_str}) "
                             f"but the insert could not be translated for "
                             f"silent mutation — pick a different region.",
                    "mutations": [],
                }
        else:
            reasons: list[str] = []
            if part_type not in coding_types:
                reasons.append(f"{part_type} is non-coding")
            else:
                if not codon_raw:
                    reasons.append("no codon table selected")
                if len(insert) % 3 != 0:
                    reasons.append(f"insert length {len(insert)} bp is "
                                   f"not a multiple of 3")
            extra = f" ({'; '.join(reasons)})" if reasons else ""
            return {
                "error": f"Internal Type IIS site(s) found: {hit_str}. "
                         f"Silent-mutation repair unavailable{extra}. "
                         f"Pick a different region or redesign.",
                "mutations": [],
            }

    # Junction-inclusive self-cut check: the bare-insert scan above cannot see a
    # Type IIS site that FORMS across a fusion boundary — an insert ending
    # ...GGTCT plus a CGCT 3' overhang spells GGTCTC (BsaI); an AATG overhang
    # plus a GTCTC insert start spells GGTCTC — which would silently self-cut
    # during the domestication / L1 digest. The insert body is internally clean
    # here (started clean, or was repaired above), so rebuild the exact cloned
    # region the real amplicon produces (same `_simulate_primed_amplicon` +
    # `_tail` trim as the synthesis path) and scan the whole thing; any hit is a
    # boundary-straddling site. Refuse rather than hand back a self-cutting part.
    _ent = (ent5, ent3) if (ent5 and ent3) else None
    try:
        _fused = _simulate_primed_amplicon(
            insert, oh5, oh3, grammar=g, part_type=part_type,
            entry_overhangs=_ent,
        )
        _tail_len = len(enzyme_pad) + len(enzyme_site) + len(enzyme_spacer)
        _region = (_fused[_tail_len:len(_fused) - _tail_len]
                   if len(_fused) > 2 * _tail_len else _fused)
        _junction_hits = _gb_find_forbidden_hits(_region, sites=forbidden_sites)
    except Exception:
        _log.exception("design_gb_primers: junction self-cut scan failed")
        _junction_hits = []
    if _junction_hits:
        _jstr = ", ".join(f"{n} {s} at +{p + 1}"
                          for n, s, p in _junction_hits)
        return {
            "error": f"A Type IIS site forms across a fusion junction in the "
                     f"cloned part ({_jstr}) — it would self-cut during "
                     f"domestication / L1 even though the insert body is clean. "
                     f"Shift the selected region a few bp (so the overhang↔insert "
                     f"junction no longer spells the site) or pick a different "
                     f"boundary.",
            "mutations": mutations,
        }

    # CDS ATG-fusion rule (regression guard 2026-05-21):
    # When the 5' overhang carries the CDS start codon (e.g.
    # AATG = 'A' spacer + 'ATG' start codon in GB L0), the
    # forward primer must NOT re-include the CDS's own first
    # 3 bp or the assembled L1 reads `...AATG ATG codon2...` —
    # duplicated start codon that frameshifts the rest of the
    # ORF. `_atg_offset_for_part` is the canonical helper that
    # encodes which (overhang, part_type) pairs trigger this
    # skip; it returns 3 for GB AATG+coding and 0 for MoClo
    # Plant AGGT (no embedded ATG) or any other custom grammar
    # whose CDS overhang doesn't carry the start codon. Extend
    # `_atg_offset_for_part` to cover new ATG-carrying
    # overhangs introduced by future grammars rather than
    # hard-coding the rule here.
    # NOTE: NO symmetric `oh3` skip. The GCTT (GB) / GCTT
    # (MoClo Plant C-tag) 3' overhang does NOT embed a stop
    # codon — the stop codon lives in the user's CDS body OR
    # in the downstream LINK's body. Stripping the last 3 bp
    # of the insert would silently drop the user's real stop.
    fwd_skip = _atg_offset_for_part(oh5, part_type)
    fwd_insert = insert[fwd_skip:] if fwd_skip else insert
    # Assemble the 5' tails FIRST so each binding region can be capped to keep
    # the TOTAL oligo within `_PRIMER_MAX_OLIGO_LEN` (binding grows to reach
    # `target_tm` but never past 50 − len(tail)) — a low-GC arm reaches
    # ~60 °C without blowing the oligo-synthesis budget.
    fwd_tail = enzyme_pad + enzyme_site + enzyme_spacer + ent5 + oh5
    rev_tail = enzyme_pad + enzyme_site + enzyme_spacer + _rc(ent3) + _rc(oh3)
    fwd_bind, fwd_tm = _pick_binding_region(
        fwd_insert, target_tm, max_len=_binding_max_len(len(fwd_tail)))
    rev_bind, rev_tm = _pick_binding_region(
        _rc(insert), target_tm, max_len=_binding_max_len(len(rev_tail)))
    rev_skip = 0   # kept as a binding so the pos calc below works

    fwd_full = fwd_tail + fwd_bind
    rev_full = rev_tail + rev_bind

    # Amplicon = pad + Esp3I + spacer + oh + insert + oh_rc + spacer_rc
    #          + Esp3I_rc + pad_rc
    amplicon_len = len(fwd_tail) + len(insert) + len(rev_tail)
    # `_simulate_primed_amplicon` fuses the duplicated start codon (the AATG
    # overhang's ATG + the insert's leading ATG collapse via
    # `_fuse_overhang_body`), so the real amplicon is `fwd_skip` bp shorter
    # whenever that fusion fires. Mirror its exact condition so the displayed
    # length matches the fragment actually built + saved.
    if fwd_skip and insert[:3].upper() == "ATG":
        amplicon_len -= fwd_skip

    # Positions of the primer binding regions on the TEMPLATE (not the
    # full amplicon). The forward primer binds the top strand at the
    # start of the insert; the reverse primer binds the bottom strand at
    # the end of the insert (positions are reported in forward-strand
    # coordinates). Save-to-library needs these to add primer_bind
    # features to the map. For wrap regions, compute positions with
    # modular arithmetic so they land on the real plasmid coordinates.
    # AATG/GCTT skips shift the binding start/end by the skipped
    # codon so the primer_bind features land on the actual primed
    # bases (codon 2 .. last-but-one codon), not the duplicated
    # start/stop.
    if wraps:
        fwd_pos = ((start + fwd_skip) % total,
                   (start + fwd_skip + len(fwd_bind)) % total)
        rev_pos = ((end - rev_skip - len(rev_bind)) % total,
                   (end - rev_skip) % total)
    else:
        fwd_pos = (start + fwd_skip,
                   start + fwd_skip + len(fwd_bind))
        rev_pos = (end - rev_skip - len(rev_bind),
                   end - rev_skip)

    pair = {
        "fwd_full":     fwd_full,
        "rev_full":     rev_full,
        "fwd_binding":  fwd_bind,
        "rev_binding":  rev_bind,
        "fwd_tm":       round(fwd_tm, 1),
        "rev_tm":       round(rev_tm, 1),
        "fwd_pos":      fwd_pos,
        "rev_pos":      rev_pos,
        "amplicon_len": amplicon_len,
    }
    # Binding-region advisory: flag any silent mutation that lands inside
    # the forward or reverse primer binding window. When non-empty, the
    # user must order the mutated insert as a gBlock and use that — not
    # the original template — as the PCR template.
    binding_region_mutations = _gb_binding_region_advisory(
        mutations, len(insert), len(fwd_bind), len(rev_bind),
        fwd_skip=fwd_skip,
    )
    return {
        "part_type":    part_type,
        "position":     pos_label,
        "oh5":          oh5,
        "oh3":          oh3,
        # The fragment's ACTUAL Type IIS-released overhangs (the domestication
        # junction that ligates into the entry vector). Equal to oh5/oh3 for a
        # classic one-tier design; the EXTERNAL entry pair (e.g. CTCG/TGAG) when
        # nested into a two-tier Golden-Braid UPD vector. Callers validate
        # entry-vector compatibility against THESE, not the category overhangs.
        "entry_oh5":    entry_oh5,
        "entry_oh3":    entry_oh3,
        "insert_seq":   insert,
        "mutations":    mutations,
        "binding_region_mutations": binding_region_mutations,
        "pairs":        [pair],
        # Segment lengths for the results painter: GB primers
        # assemble as `pad + enzyme_site + spacer + overhang +
        # binding`. Exposed so `_show_result` can highlight the
        # RE recognition seq in blue and the overhang/spacer as
        # padding rather than leaving the whole tail gray.
        "enzyme_site":   enzyme_site,
        "enzyme_pad":    enzyme_pad,
        "enzyme_spacer": enzyme_spacer,
        # The binding-portion Tm this design aimed for. The achieved
        # fwd_tm / rev_tm are the closest reachable within the allowed
        # binding window (`_pick_binding_region`); `_show_result` flags any arm
        # whose achieved Tm lands off this target so the next-best choice
        # is visible to the user.
        "target_tm":     float(target_tm),
        # Legacy top-level mirror of pairs[0] for callers (cloning simulator,
        # PrimerDesignScreen) that don't iterate the list yet.
        **pair,
    }


def _released_insert_from_fragment(frag_seq: str, enzyme: str) -> dict:
    """Digest a LINEAR synthesised fragment and return the piece the enzyme
    releases — the one cut on BOTH ends, i.e. what actually enters a vector.

    Refuses loudly rather than guessing: a fragment that isn't a properly
    wrapped L0 fragment (no sites, one site, or extra internal sites) would
    otherwise yield a plausible-looking part that cannot assemble. Raises
    ValueError with a message naming the actual cause."""
    seq = "".join(c for c in (frag_seq or "").upper()
                  if c in "ACGTRYSWKMBDHVN")
    if not seq:
        raise ValueError("the fragment has no sequence")
    if not enzyme:
        raise ValueError("this grammar has no Type IIS enzyme configured")
    if enzyme not in _state._all_enzymes_hook():
        # Without this the digest simply finds no sites and the fragment gets
        # blamed — sending the user to rebuild a perfectly good fragment when
        # the fault is a custom grammar naming an enzyme nothing can cut with.
        raise ValueError(
            f"this grammar's enzyme {enzyme!r} isn't in the enzyme catalog — "
            f"fix the grammar (or add the enzyme under Settings → Enzymes); "
            f"the fragment is not the problem.")
    pieces = _digest_with_enzymes(seq, [enzyme], circular=False)
    cut_both = [p for p in pieces
                if p["left"].get("kind") != "linear"
                and p["right"].get("kind") != "linear"]
    if not cut_both:
        n_cuts = sum(1 for p in pieces
                     if p["right"].get("kind") != "linear")
        if n_cuts == 0:
            raise ValueError(
                f"no {enzyme} site found — this is a plain fragment, not a "
                f"wrapped L0 fragment. Build one with Synthesis → L0 Fragment.")
        raise ValueError(
            f"only one {enzyme} cut — an L0 fragment needs the enzyme on BOTH "
            f"ends so the insert can be released.")
    if len(cut_both) > 1:
        raise ValueError(
            f"{len(cut_both) + 1} {enzyme} cuts — the insert carries an extra "
            f"internal site, so the digest would fragment it. Scrub the site "
            f"(Optimize) and rebuild the fragment.")
    return cut_both[0]


def _clone_syn_fragment_into_entry_vector(
        frag_seq: str, vector_seq: str, *, grammar: dict,
        frag_features: "list[dict] | None" = None,
        vector_features: "list[dict] | None" = None,
) -> dict:
    """Simulate cloning a SYNTHESISED L0 fragment into its entry vector.

    The real bench steps, in order ([INV-127] — simulate, never hand-build the
    product): digest the fragment and the vector with the grammar's Type IIS
    enzyme, ligate the released insert into the surviving backbone, close the
    circle. The fragment's own layout is never reverse-engineered — the
    overhangs come from the digest, so the two-tier nesting
    (`_build_synthesis_l0_fragment`'s entry overhangs wrapping the category
    pair) needs no special handling and an AATG overlapping a body's start
    codon can't be miscounted.

    Deliberately does NOT report the part's fusion overhangs: reading them
    back off the finished plasmid is the very inference [INV-172] exists to
    forbid — a two-tier fragment answers with the ENTRY pair. The caller takes
    them from the grammar's position table instead.

    Returns ``{sequence, features, insert_len, enzyme, vector_len}``.
    Raises ValueError, message-first, on any step that can't happen."""
    enzyme = str(grammar.get("enzyme") or "")
    insert = _released_insert_from_fragment(frag_seq, enzyme)
    if frag_features:
        # Re-digest carrying the features so they ride onto the insert in
        # fragment-local coords instead of being dropped.
        seq = "".join(c for c in (frag_seq or "").upper()
                      if c in "ACGTRYSWKMBDHVN")
        for p in _digest_with_enzymes(seq, [enzyme], circular=False,
                                      features=list(frag_features)):
            if (p["left"].get("kind") != "linear"
                    and p["right"].get("kind") != "linear"):
                insert = p
                break
    vseq = "".join(c for c in (vector_seq or "").upper()
                   if c in "ACGTRYSWKMBDHVN")
    if not vseq:
        raise ValueError("the entry vector has no sequence")
    vpieces = _digest_with_enzymes(vseq, [enzyme], circular=True,
                                   features=list(vector_features or []))
    if len(vpieces) < 2:
        raise ValueError(
            f"the entry vector has no {enzyme} dropout to replace — it must "
            f"be cut twice to open a slot for the insert.")
    # The dropout is whichever piece presents the SAME overhang pair as the
    # insert; everything else chains into the backbone. Mirrors
    # `_clone_part_into_entry_vector`'s handling of vectors that picked up an
    # extra site: any number of survivors ligate end to end.
    want = (insert["left"].get("overhang_seq", ""),
            insert["right"].get("overhang_seq", ""))
    keep = [p for p in vpieces
            if (p["left"].get("overhang_seq", ""),
                p["right"].get("overhang_seq", "")) != want]
    if len(keep) == len(vpieces):
        raise ValueError(
            f"the insert's {want[0]}/{want[1]} overhangs match nothing the "
            f"vector releases — wrong entry vector for this grammar, or the "
            f"fragment was built against a different one.")
    if not keep:
        # Every piece looked like the dropout, which only happens on a
        # degenerate vector whose two cuts leave the SAME overhang on both
        # ends. Ligating nothing would close the insert on itself and report a
        # backbone-free "plasmid" as a success.
        raise ValueError(
            f"the entry vector leaves {want[0]}/{want[1]} at both cuts, so "
            f"nothing survives as a backbone — it can't take a directional "
            f"insert.")
    chain = insert
    for p in keep:
        merged = _ligate_fragments(chain, p)
        if merged is None:
            raise ValueError(
                "the vector backbone pieces don't ligate onto the insert — "
                "check that the fragment and vector share a grammar.")
        chain = merged
    closed = _close_circular(chain)
    if closed is None:
        raise ValueError(
            "the assembly won't close into a circle — the free ends don't "
            "match.")
    return {
        "sequence":   closed["top_seq"],
        "features":   closed["features"],
        "insert_len": len(insert["top_seq"]),
        "enzyme":     enzyme,
        "vector_len": len(vseq),
    }


def _grammar_position_for_type(grammar: dict, part_type: str) -> "dict | None":
    """Case-tolerant lookup of a grammar's position row for ``part_type``.

    Defers to the existing exact-match `_grammar_position_by_type` first, so
    there is ONE definition of what a position is; the fallback only relaxes
    case, which the agent API needs (a caller sending `"cds"` means `CDS`).
    Deliberately does NOT match on the position's `name` — names like
    `Pos 3-4` could alias another row's type."""
    hit = _grammar_position_by_type(grammar, part_type)
    if hit is not None:
        return hit
    want = str(part_type or "").strip().lower()
    if not want:
        return None
    for p in (grammar.get("positions") or []):
        if isinstance(p, dict) and str(p.get("type") or "").strip().lower() == want:
            return p
    return None


def _l0_part_from_syn_fragment(
        frag_seq: str, vector_seq: str, *, grammar: dict, part_type: str,
        name: str = "", frag_features: "list[dict] | None" = None,
        vector_features: "list[dict] | None" = None) -> dict:
    """Turn a synthesised L0 fragment into a Level-0 PART of ``grammar``.

    The caller chooses the part type; its fusion overhangs come from the
    grammar's own position table, never from inferring them off the sequence.
    Two guards make "wrong grammar / wrong part type" fail loudly instead of
    producing a part that looks right and cannot assemble:

      * the part type must BE a position in this grammar (so a MoClo type
        can't be filed under a Golden Braid grammar, and vice versa);
      * the overhangs that type declares must actually be present at the ends
        of the released insert — picking `CDS` for a fragment built as a
        promoter would otherwise store `AATG/GCTT` against a body that
        carries `GGAG/AATG`.

    The stored ``sequence`` is the body BETWEEN those overhangs, never the
    flanked region: the L0 storage contract holds `oh5`/`oh3` separately and
    the assembly chain re-adds them (`_assembly_fragment_strip_oh5`,
    `_clone_part_into_entry_vector`), so storing them inline duplicates 8 bp
    per part — a frameshift in a CDS. For a two-tier fragment the entry
    overhangs come off too; they belong to the vector junction, not the part.

    Returns a parts-bin-ready dict plus the cloned plasmid under
    ``cloned_seq`` / ``cloned_features``. Raises ValueError, message first."""
    pos = _grammar_position_for_type(grammar, part_type)
    if pos is None:
        known = ", ".join(
            str(p.get("type") or p.get("name") or "?")
            for p in (grammar.get("positions") or []) if isinstance(p, dict)
        ) or "none"
        raise ValueError(
            f"{grammar.get('name') or grammar.get('id') or 'this grammar'} has "
            f"no {part_type!r} position — it defines: {known}.")
    oh5 = str(pos.get("oh5") or "").upper()
    oh3 = str(pos.get("oh3") or "").upper()
    if len(oh5) != 4 or len(oh3) != 4:
        raise ValueError(
            f"the {part_type!r} position has no usable 4-nt overhang pair.")
    enzyme = str(grammar.get("enzyme") or "")
    insert = _released_insert_from_fragment(frag_seq, enzyme)
    # The region that lands in the vector is the insert's top strand plus the
    # overhang carried by the NEXT piece (see `_fragments_from_cuts`: top_seq
    # holds the left overhang, not the right).
    region = insert["top_seq"] + insert["right"].get("overhang_seq", "")
    # Locate the CATEGORY overhangs at the two EXACT offsets the layout
    # allows — flush (one-tier) or one 4-nt entry overhang in (two-tier). A
    # substring search over the same window would pass a body that merely
    # happens to begin with the other type's overhang, filing a part whose
    # declared ends aren't its real ones. The offsets double as the trim: what
    # sits outside the category pair is vector junction, not part.
    if region[:4] == oh5:
        off5 = 0
    elif region[4:8] == oh5:
        off5 = 4
    else:
        raise ValueError(
            f"this fragment doesn't start with the {oh5} overhang that "
            f"{part_type!r} requires — it was built as a different part type.")
    if region[-4:] == oh3:
        off3 = 0
    elif region[-8:-4] == oh3:
        off3 = 4
    else:
        raise ValueError(
            f"this fragment doesn't end with the {oh3} overhang that "
            f"{part_type!r} requires — it was built as a different part type.")
    body = region[off5 + 4:len(region) - off3 - 4]
    if not body:
        raise ValueError(
            f"there is nothing between the {oh5} and {oh3} overhangs — this "
            f"fragment carries no insert.")
    cloned = _clone_syn_fragment_into_entry_vector(
        frag_seq, vector_seq, grammar=grammar,
        frag_features=frag_features, vector_features=vector_features)
    return {
        "name":       (name or "").strip() or "part",
        "type":       str(pos.get("type") or part_type),
        "position":   str(pos.get("name") or pos.get("type") or part_type),
        "oh5":        oh5,
        "oh3":        oh3,
        "sequence":   body,
        "nested":     bool(off5 or off3),
        "grammar":    str(grammar.get("id") or ""),
        "level":      0,
        "enzyme":     enzyme,
        "cloned_seq": cloned["sequence"],
        "cloned_features": cloned["features"],
        "insert_len": cloned["insert_len"],
    }


def _build_synthesis_l0_fragment(
    sequence: str, oh5: str, oh3: str, *,
    grammar: "dict | None" = None,
    part_type: str = "",
    entry_overhangs: "tuple[str, str] | None" = None,
) -> dict:
    """Wrap ``sequence`` in the grammar's Type IIS ends + fusion overhangs so the
    directly-SYNTHESISED fragment (order it as a gBlock — no PCR) becomes a
    proper Level-0 part once cloned into a UPD/pUPD entry vector.

    Shared by the Synthesis-tab "L0 Fragment" button and the ``design-
    synthesis-fragment`` agent endpoint, so a human and an agent/Babs produce
    the SAME fragment. Two-tier aware: pass the vector's EXTERNAL acceptor
    overhangs (resolve with :func:`_entry_vector_acceptor_overhangs`) as
    ``entry_overhangs`` and the category pair (``oh5``/``oh3``) nests inside
    them; ``None`` — or a pair equal to the category — yields the classic
    one-tier fragment.

    Returns ``{fragment, length, oh5, oh3, entry_oh5, entry_oh3, part_type,
    enzyme, nested, internal_sites}`` where ``entry_oh5``/``entry_oh3`` are the
    fragment's ACTUAL Type IIS-cut overhangs (what enters the vector) and
    ``internal_sites`` lists any grammar-forbidden recognition site FOUND in the
    insert (each ``{enzyme, site, pos}``, 0-based) — those must be scrubbed (e.g.
    via codon Optimize) before ordering or the fragment self-cuts."""
    g = grammar if isinstance(grammar, dict) else _BUILTIN_GRAMMARS["gb_l0"]
    seq = "".join(c for c in (sequence or "").upper() if c in "ACGTRYSWKMBDHVN")
    cat5 = (oh5 or "").upper()
    cat3 = (oh3 or "").upper()
    ext: "tuple[str, str] | None" = None
    if isinstance(entry_overhangs, (tuple, list)) and len(entry_overhangs) == 2:
        e5 = str(entry_overhangs[0] or "").upper()
        e3 = str(entry_overhangs[1] or "").upper()
        if len(e5) == 4 and len(e3) == 4 and (e5, e3) != (cat5, cat3):
            ext = (e5, e3)
    fragment = _simulate_primed_amplicon(
        seq, cat5, cat3, grammar=g, part_type=part_type, entry_overhangs=ext)
    forbidden = g.get("forbidden_sites") or {}
    # Scan the CLONED region (everything between the two intended Type IIS
    # tails), NOT the raw insert — a forbidden site can straddle a fusion
    # junction (e.g. the AATG overhang + a GTCTC insert start = a BsaI GGTCTC)
    # and would silently self-cut during domestication / L1. `_tail` is the
    # pad+site+spacer machinery on each end; positions are reported
    # FRAGMENT-relative so the warning locates the site in the built fragment.
    _tail = (len(g.get("pad", _GB_PAD)) + len(g.get("site", _GB_L0_ENZYME_SITE))
             + len(g.get("spacer", _GB_SPACER)))
    _region = (fragment[_tail:len(fragment) - _tail]
               if len(fragment) > 2 * _tail else fragment)
    hits = ([(n, s, p + _tail)
             for (n, s, p) in _gb_find_forbidden_hits(_region, sites=forbidden)]
            if forbidden else [])
    return {
        "fragment":       fragment,
        "length":         len(fragment),
        "oh5":            cat5,
        "oh3":            cat3,
        "entry_oh5":      ext[0] if ext else cat5,
        "entry_oh3":      ext[1] if ext else cat3,
        "part_type":      part_type,
        "enzyme":         g.get("enzyme") or "",
        "nested":         ext is not None,
        "internal_sites": [{"enzyme": n, "site": s, "pos": p}
                           for (n, s, p) in hits],
    }


# Standard column oligo synthesis tops out ~60-100 nt; a mutagenic SOE window
# spanning several clustered edits can exceed that. We WARN (never truncate —
# dropping cured bases would be catastrophic) so the user orders an ultramer or
# splits the edits.
_SOE_PRIMER_WARN_LEN = 100


def _design_operon_soe_primers(
    operon_seq: str,
    feats: "list[dict] | None",
    grammar: dict,
    *,
    manual_edits: "list[dict] | None" = None,
    extra_enzymes: "list[str] | None" = None,
    codon_raw: "dict | None" = None,
    target_tm: float = 60.0,
    overlap_arm: int = 18,
) -> dict:
    """Design an SOE-PCR primer set that LIFTS a native operon off its template,
    cures every grammar-forbidden Type IIS site, and adds the grammar's TU
    cassette to the ends — ready to clone into the grammar's entry vector.

    ``extra_enzymes`` (names, e.g. ``["EcoRI", "KpnI"]``) are cured ALONGSIDE the
    grammar's Type IIS set — resolved against the merged enzyme catalog and
    treated identically (synonymous inside CDS, flagged outside), so the lifted
    operon can be scrubbed of downstream cloning sites too, not just the
    assembly enzymes.

    The user amplifies NATIVE DNA (no synthesis), so every cure must be carried
    by a primer. Curing = synonymous substitutions inside CDS features
    (`_scrub_design`, multi-CDS, reverse-strand- and stop-safe) PLUS any
    ``manual_edits`` ({pos, to}) the user marked for non-coding sites. Edits are
    clustered (`_SCRUB_PRIMER_FOOTPRINT`); each cluster becomes an SOE junction
    with a complementary mutagenic primer pair carrying the cured bases, and the
    two outermost primers carry ``pad + site + spacer + TU-overhang + binding``.

    Returns one of:
      * ``{ok: True, cured_seq, edits, primers:[{name, seq, kind, tm, ...}],
            n_clusters, tu_overhangs, amplicon_len, warnings}``
      * ``{needs_manual: True, sites_skipped:[{enzyme, site, pos, in_cds}],
            cured_seq, edits, warnings}`` — forbidden sites remain (non-coding,
            or a CDS site with no synonymous cure); caller collects
            ``manual_edits`` for the ``in_cds=False`` ones and re-runs.
      * ``{error: str}``

    CATASTROPHIC-CLASS GATE: a cure is only real if a primer carries it. Before
    returning ``ok`` every edit position is checked to fall inside a
    primer-covered window (a flank binding region or a mutagenic window);
    otherwise PCR off the native template would silently leave that site intact,
    so the design is REFUSED ([[project_primer_design_catastrophic]])."""
    operon_seq = (operon_seq or "").upper()
    n = len(operon_seq)
    if n < 40:
        return {"error": f"Operon is too short ({n} bp) for SOE domestication."}
    forbidden = dict(grammar.get("forbidden_sites") or _GB_DOMESTICATION_FORBIDDEN)
    if extra_enzymes:
        # Merge user-requested cloning enzymes (EcoRI / KpnI / …) into the
        # forbidden set, resolved against the merged built-in + custom catalog.
        forbidden.update(_scrub_resolve_sites(extra_enzymes))
    enzymes = list(forbidden.keys())

    # 1. Auto-cure synonymous CDS sites via the Scrub engine (strand/stop-safe,
    #    multi-CDS; it skips non-coding sites — those need manual override).
    scrub = _scrub_design(operon_seq, feats, enzymes=enzymes,
                          circular=False, codon_raw=codon_raw)
    if not scrub.get("ok", True):
        return {"error": "; ".join(scrub.get("warnings") or ["cure failed"])}
    # Auto-apply ONLY synonymous cures whose changed base sits INSIDE a CDS
    # feature — the rule is "mutate inside CDS, flag outside" (the Scrub engine
    # will freely substitute a non-coding base, but a stray change could break
    # an RBS / regulator, so we leave those for manual override). A site that
    # needs a non-coding edit is therefore left intact → flagged below.
    cds_spans = [(int(f["start"]), int(f["end"]))
                 for f in (feats or [])
                 if str(f.get("type", "")).upper() == "CDS"
                 and f.get("start") is not None and f.get("end") is not None]

    def _in_cds(p: int) -> bool:
        return any(s <= p < e for s, e in cds_spans)

    working = list(operon_seq)
    edits: "list[dict]" = []
    for e in scrub.get("edits", []):
        p = int(e["pos"])
        if not _in_cds(p):
            continue
        working[p] = str(e["to"]).upper()
        edits.append({"pos": p, "to": str(e["to"]).upper(),
                      "frm": operon_seq[p], "region": "CDS"})

    # 2. Apply user manual edits (single-base substitutions the user marked for
    #    non-coding sites — explicitly authorized, so they bypass the CDS rule).
    for m in (manual_edits or []):
        try:
            p, b = int(m["pos"]), str(m["to"]).upper()
        except (KeyError, TypeError, ValueError):
            continue
        if 0 <= p < n and b in "ACGT":
            working[p] = b
            edits.append({"pos": p, "to": b, "frm": operon_seq[p],
                          "region": "manual"})
    cured = "".join(working)
    if len(cured) != n:
        return {"error": "internal: cure changed the sequence length"}

    # 3. Any forbidden site still present? → caller must supply manual edits.
    remaining = _gb_find_forbidden_hits(cured, sites=forbidden)
    if remaining:
        flagged = [{"enzyme": nm, "site": st, "pos": ps,
                    "in_cds": _in_cds(ps)}
                   for nm, st, ps in remaining]
        return {"needs_manual": True, "sites_skipped": flagged,
                "cured_seq": cured, "edits": edits,
                "warnings": scrub.get("warnings", [])}

    # 4. Cluster edits (linear — no origin merge) into SOE junctions.
    positions = sorted(e["pos"] for e in edits)
    clusters: "list[list[int]]" = []
    if positions:
        clusters = [[positions[0]]]
        for p in positions[1:]:
            if p - clusters[-1][-1] <= _SCRUB_PRIMER_FOOTPRINT:
                clusters[-1].append(p)
            else:
                clusters.append([p])

    # 5. Flanking primers carry the grammar cassette. An operon clones as a
    #    CDS-EQUIVALENT L0 part, so use the OPERON position's overhangs (AATG
    #    carrying the first gene's start codon → GCTT before a terminator),
    #    falling back to CDS, then the TU boundary.
    pad = grammar.get("pad", _GB_PAD)
    site = grammar.get("site", _GB_L0_ENZYME_SITE)
    spacer = grammar.get("spacer", _GB_SPACER)
    oh = _grammar_position_overhangs(grammar, ("OPERON", "CDS"))
    op_oh5, op_oh3 = oh if oh else _grammar_tu_overhangs(grammar)
    if not op_oh5 or not op_oh3:                       # custom grammar w/o any
        op_oh5, op_oh3 = "AATG", "GCTT"
    fwd_tail = pad + site + spacer + op_oh5
    rev_tail = pad + site + spacer + _rc(op_oh3)
    # ATG-fusion: when the 5' overhang carries the start codon (AATG) AND the
    # operon really begins with the first gene's ATG, the forward primer binds
    # at codon 2 so the assembled L1 reads ...AATG + codon2... (no duplicated
    # ATG) — exactly the CDS convention (`_atg_offset_for_part`).
    fwd_skip = (_atg_offset_for_part(op_oh5, "OPERON")
                if cured[:3] == "ATG" else 0)
    fwd_bind, fwd_tm = _pick_binding_region(
        cured[fwd_skip:], target_tm, max_len=_binding_max_len(len(fwd_tail)))
    rev_bind, rev_tm = _pick_binding_region(
        _rc(cured), target_tm, max_len=_binding_max_len(len(rev_tail)))
    primers: "list[dict]" = []
    fwd_name, rev_name = _dom_primer_pair_names("operon", 1)
    primers.append({"name": fwd_name, "seq": fwd_tail + fwd_bind,
                    "kind": "flank-fwd", "tm": round(fwd_tm, 1),
                    "covers": [fwd_skip, fwd_skip + len(fwd_bind)]})
    primers.append({"name": rev_name, "seq": rev_tail + rev_bind,
                    "kind": "flank-rev", "tm": round(rev_tm, 1),
                    "covers": [n - len(rev_bind), n]})
    # Primer-covered windows: the two flank binding regions sit at the ends
    # (the forward one starts past the skipped start codon).
    cover: "list[tuple[int, int]]" = [(fwd_skip, fwd_skip + len(fwd_bind)),
                                      (n - len(rev_bind), n)]

    # 6. Internal mutagenic pairs — one per cluster, carrying the cured bases.
    #    Window = cluster span ± overlap_arm (clamped); the fwd primer IS the
    #    cured window, the rev its reverse-complement → a fully-complementary
    #    SOE overlap that anneals during the assembly PCR.
    warnings = list(scrub.get("warnings") or [])
    for j, cl in enumerate(clusters, start=2):
        ws = max(0, cl[0] - overlap_arm)
        we = min(n, cl[-1] + 1 + overlap_arm)
        if we - ws < 2 * overlap_arm:                 # near an end → widen in
            if ws == 0:
                we = min(n, 2 * overlap_arm)
            else:
                ws = max(0, n - 2 * overlap_arm if we == n else we - 2 * overlap_arm)
        win = cured[ws:we]
        if len(win) > _SOE_PRIMER_WARN_LEN:
            warnings.append(
                f"cluster {j - 1}: SOE mutagenic primer is {len(win)} nt "
                f"(> {_SOE_PRIMER_WARN_LEN}) — near/over the standard oligo-"
                f"synthesis limit; order as an ultramer or split the edits.")
        fj, rj = _dom_primer_pair_names("operon", j)
        primers.append({"name": fj, "seq": win, "kind": "soe-fwd",
                        "tm": round(_primer_tm(win) or 0.0, 1),
                        "covers": [ws, we]})
        primers.append({"name": rj, "seq": _rc(win), "kind": "soe-rev",
                        "tm": round(_primer_tm(_rc(win)) or 0.0, 1),
                        "covers": [ws, we]})
        cover.append((ws, we))

    # 7. CATASTROPHIC-CLASS GATE — every cure must ride a primer.
    for e in edits:
        if not any(s <= e["pos"] < q for s, q in cover):
            return {"error": (
                f"cure at +{e['pos'] + 1} ({e.get('frm')}→{e['to']}) "
                "falls outside every primer — the SOE set would not introduce "
                "it. Aborting (catastrophic-class primer safety).")}
    if _gb_find_forbidden_hits(cured, sites=forbidden):
        return {"error": "internal: forbidden site survived the cure"}

    return {"ok": True, "cured_seq": cured, "edits": edits, "primers": primers,
            "n_clusters": len(clusters), "overhangs": [op_oh5, op_oh3],
            "fwd_skip": fwd_skip,
            "amplicon_len": len(fwd_tail) + (n - fwd_skip) + len(rev_tail),
            "warnings": warnings}


def _grammar_tu_overhangs(grammar: dict) -> tuple[str, str]:
    """Return ``(tu_start, tu_end)`` — the boundary overhangs of a
    full TU under this grammar. ``tu_start`` is the first position's
    `oh5` (Promoter side); ``tu_end`` is the last position's `oh3`
    (Terminator side)."""
    positions = grammar.get("positions") or []
    if not positions:
        return ("", "")
    return (
        str(positions[0].get("oh5") or ""),
        str(positions[-1].get("oh3") or ""),
    )


def _grammar_position_overhangs(grammar: dict,
                                type_names: "tuple[str, ...]"
                                ) -> "tuple[str, str] | None":
    """Return ``(oh5, oh3)`` of the FIRST of ``type_names`` present in this
    grammar's position list (matched by part-type name), or ``None`` if none
    are defined. Lets the Native Operon designer prefer the ``OPERON`` position
    (CDS overhangs), fall back to ``CDS``, then to the TU boundary."""
    by_type = {str(p.get("type") or ""): p
               for p in (grammar.get("positions") or [])}
    for nm in type_names:
        p = by_type.get(nm)
        if p and p.get("oh5") and p.get("oh3"):
            return str(p["oh5"]), str(p["oh3"])
    return None


# ── PCR backend caps (Phase D) — shared by the hub PCR simulator (_simulate_pcr,
# which stays hub-side via its _state hook) + the agent PCR endpoint's input
# validation; relocated here so both resolve them (the early cloning re-export
# covers the hub default-arg uses).
_PCR_MIN_PRIMER_LEN     = 10       # primers shorter than this can't anneal


_PCR_MAX_PRIMER_LEN     = 80       # absurdly long primer = user error


_PCR_DEFAULT_MAX_AMPLICON = 20_000   # bp — function / agent default ceiling


# UI default for the Simulator's "Max amplicon" box: 500 bp is a common
# amplicon size, a friendlier starting point than the 20 kb long-PCR
# ceiling. UI-only — the function / agent default above is unchanged.
_PCR_UI_DEFAULT_MAX_AMPLICON = 500


_PCR_AMPLICON_HARD_CAP  = 100_000  # bp — safety cap regardless of UI input


_PCR_MAX_AMPLICONS      = 50       # cap on result count (a mispriming primer


                                   # on a repetitive template can yield 1000s)
_PCR_MAX_TEMPLATE_BP    = 5_000_000  # 5 Mb — above this we skip the run rather


                                     # than freeze the UI on chromosome-scale
                                     # inputs (genome chunks via FASTA import
                                     # routinely break this threshold).
# Pathological case: a 10-bp ACGT-only primer on a 5 Mb template can yield
# ≈ 4,768 expected hits at random; an A-rich primer on an A-rich tract can
# yield orders of magnitude more. The fwd × rev double-loop is O(N²); cap
# either side at this many positions and refuse — surfacing a clearer error
# than a multi-second UI freeze on a pure-A primer.
_PCR_MAX_PRIMER_HITS    = 5_000
