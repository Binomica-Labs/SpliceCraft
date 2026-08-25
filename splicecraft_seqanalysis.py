"""splicecraft_seqanalysis — sequence analysis & part classification (Phase D, layer L3).

Read-only analysis of a loaded sequence: the six-frame ORF finder (`_find_orfs`)
and the Golden-Braid / MoClo part classifier (`_classify_part_from_plasmid` — which
cloning grammar + position + level a circular plasmid represents) plus its helpers
(`_check_vector_match` / `_vector_half_top_seq`, `_grammar_acceptor_tu_pairs` /
`_grammar_canonical_overhangs`, `_fragment_has_backbone_marker`, `_part_level_label`).

Layer L3 (same as cloning): the classifier digests the candidate plasmid via cloning's
`_excise_fragment_pair` + reads `_grammar_tu_overhangs`, so it sits at cloning's layer —
cloning never imports this module, so no cycle. Other deps are strictly lower: codon L2
(`_CODON_TABLE` / `_STOP_CODONS`), record L1 (`_gb_text_to_record`), dataaccess L1
(`_all_grammars` / `_load_entry_vectors`), biology L0 (`_rc`), logging L0.

The two classification LRU caches live in `_state` (`_VECTOR_MATCH_CACHE` /
`_ACCEPTOR_TU_PAIRS_CACHE`) so the hub-side `_after_entry_vectors_save` bust
(`globals().get(...).clear()`) and this sibling's reads/writes hit the SAME dict —
read them as `_state._VECTOR_MATCH_CACHE`, never by a stale by-value import. The caps
(`_VECTOR_MATCH_CACHE_MAX` / `_ACCEPTOR_TU_PAIRS_CACHE_MAX`) are plain consts here.
Re-exported by the hub so `sc.<name>` + every call site resolves unchanged.
"""
from __future__ import annotations

import re

import splicecraft_state as _state
from splicecraft_biology import (_rc, _enzyme_cuts, _feat_len,
                                 _slice_circular)
from splicecraft_splice import (_splice_scan, _splice_model_summary)
from splicecraft_codon import _CODON_TABLE, _STOP_CODONS
from splicecraft_record import _gb_text_to_record
from splicecraft_dataaccess import _all_grammars, _load_entry_vectors
from splicecraft_cloning import _excise_fragment_pair, _grammar_tu_overhangs
from splicecraft_logging import _log, _log_event


def _find_orfs(seq: str, *,
               min_aa: int = 30,
               include_alt_starts: bool = False,
               circular: bool = True) -> list[dict]:
    """Six-frame ORF scan. Returns
    ``[{start, end, strand, length_aa, nt_len, exceeds_one_lap, aa_seq},
    ...]`` sorted by length descending. Wrap-aware on circular plasmids: an
    ORF crossing the origin is reported with ``end < start``, matching the
    wrap-feature convention `_feat_len` / `_bp_in` already implement.

    ``nt_len`` is the exact coding length in bases INCLUDING the stop codon.
    Read it (or ``length_aa``) for length — never derive length from
    ``(start, end)``, which on a circle cannot express an ORF of a full lap
    or more. Such an ORF (a frame with no in-frame stop for the whole
    molecule — reachable on short synthetic constructs) is reported with
    ``exceeds_one_lap=True`` and its span pinned to the near-full circle.

    `min_aa` excludes the stop codon (so ``min_aa=30`` ⇒ ORFs ≥ 30
    coded residues, i.e. ≥ 93 bp including the trailing stop).
    `include_alt_starts=True` adds GTG and TTG to the start-codon
    set — useful for bacterial genomes; off by default since most
    plasmid CDSes use ATG.
    """
    n = len(seq)
    if n < 6:
        return []
    seq_u = seq.upper()
    starts = {"ATG"}
    if include_alt_starts:
        starts |= {"GTG", "TTG"}

    orfs: list[dict] = []

    for strand in (1, -1):
        if strand == 1:
            scan_seq = (seq_u + seq_u) if circular else seq_u
        else:
            rc_seq = _rc(seq_u)
            scan_seq = (rc_seq + rc_seq) if circular else rc_seq
        scan_n = len(scan_seq)

        for frame in range(3):
            current_start = -1
            i = frame
            while i + 3 <= scan_n:
                codon = scan_seq[i:i+3]
                if current_start < 0:
                    if codon in starts:
                        current_start = i
                else:
                    if codon in _STOP_CODONS:
                        aa_len = (i - current_start) // 3
                        # Drop ORFs whose start codon falls in the
                        # second copy of the doubled scan — they're
                        # duplicates of one we already found from the
                        # first copy.
                        if (circular and current_start >= n) \
                                or aa_len < min_aa:
                            current_start = -1
                            i += 3
                            continue
                        nt_seq = scan_seq[current_start:i + 3]  # incl. stop
                        aa_seq = "".join(
                            _CODON_TABLE.get(nt_seq[k:k+3], "?")
                            for k in range(0, len(nt_seq), 3)
                        )
                        if strand == 1:
                            o_s = current_start
                            o_e = i + 3
                            if circular and o_e > n:
                                o_e -= n   # wrap: end < start
                        else:
                            p_rc = current_start
                            e_rc = i + 3
                            if circular:
                                o_s = (n - e_rc) % n
                                o_e = (n - p_rc) % n
                            else:
                                o_s = n - e_rc
                                o_e = n - p_rc
                        # Coding length in bases, INCLUDING the stop codon —
                        # reported explicitly so no consumer ever has to infer
                        # it from `(start, end)`, which on a circle cannot.
                        nt_len = len(nt_seq)
                        # An ORF that runs a full lap or more (a frame with no
                        # in-frame stop for the whole molecule) has NO (start,
                        # end) pair on an n-bp circle that can express it:
                        # `(o_s + nt_len) % n` lands back inside the ORF, so the
                        # arc drawn from it UNDER-states the ORF by whole turns
                        # — an exactly-one-lap ORF even collapses to `o_e ==
                        # o_s`, which `_feat_len` reads as length 0. Pin the
                        # span to the near-full circle (the honest "it goes all
                        # the way round" arc) and flag it; `length_aa` /
                        # `nt_len` / `aa_seq` stay exact either way.
                        over_lap = circular and nt_len >= n
                        if over_lap:
                            o_e = (o_s - 1) % n
                        orfs.append({
                            "start":     o_s,
                            "end":       o_e,
                            "strand":    strand,
                            "length_aa": aa_len,
                            "nt_len":    nt_len,
                            "exceeds_one_lap": over_lap,
                            "aa_seq":    aa_seq,
                        })
                        current_start = -1
                i += 3

    # Dedupe identical (start, end, strand) tuples — can happen when the
    # doubled-scan cycles past the origin and re-finds the same ORF, or
    # when alt-starts inside an ATG-bounded region produce nested hits
    # that then collapse onto the same boundary.
    orfs.sort(key=lambda o: o["length_aa"], reverse=True)
    seen: set[tuple] = set()
    unique: list[dict] = []
    for o in orfs:
        k = (o["start"], o["end"], o["strand"])
        if k in seen:
            continue
        seen.add(k)
        unique.append(o)
    return unique


# Backbone-marker keywords used by `_pick_insert_fragment` to tell the
# bacterial / replication-machinery half of a digested plasmid apart
# from the cloned insert. Match is case-insensitive substring on
# feature labels + types — covers the common annotation styles
# (`rep_origin`, `Ori*`, `AmpR`, `KanR`, `cat`, `Spec`, etc.) without
# trying to be exhaustive: false negatives just trigger the size
# fallback, which is fine.
_BACKBONE_FEATURE_TYPES: frozenset[str] = frozenset({
    "rep_origin", "oriT",
})


_BACKBONE_LABEL_KEYWORDS: tuple[str, ...] = (
    "ori", "rep_origin",
    "ampr", "kanr", "specr", "specinomycin", "spectinomycin",
    "cmr", "chloramphenicol",
    "carbr", "carbenicillin",
    # NB: "tetr" / "tetracyclin" USED to be here and are now whole-token +
    # gene-bearing-type keywords (`_BACKBONE_TOKEN_KEYWORDS`) — see the tet
    # note there. Every OTHER keyword here is long enough to survive as a
    # word-start match.
    # NB: bare "selection"/"antibiotic" deliberately excluded — a ccdB/sacB
    # counter-selection DROPOUT is often labeled "...selection cassette", and
    # matching it would tag the dropout as backbone, collapsing the marker-
    # absence dropout test back to the size heuristic it was meant to replace.
    # Specific resistance genes + origins only. Keep in sync with cloning's
    # `_ACCEPTOR_BACKBONE_LABEL_KEYWORDS`.
)


# Suffixes that must NOT follow a keyword: the English words and unrelated
# biology that keyword happens to OPEN. A word-start match alone doesn't stop
# these — "orientation" and "original" both start with a real keyword.
#   ori  → orient(ation), origina(l) …  but "origin" / "origins" must match
_BACKBONE_LABEL_TRAPS: "dict[str, tuple[str, ...]]" = {
    "ori": ("ent", "gina"),
}


def _label_keyword_re(
    keywords: "tuple[str, ...]",
    traps: "dict[str, tuple[str, ...]]",
) -> "re.Pattern[str]":
    r"""Compile marker keywords into one alternation that matches only at a
    WORD START — never mid-word — with a per-keyword suffix blocklist.

    THE bug this exists to prevent (field report 2026-08-17, agent-driven
    Golden Gate): the keywords used to be matched as bare case-insensitive
    substrings, and `"ori"` is three letters. It is the tail of `EcoRI` (also
    `EcoRI site`, `EcoRI-F`), of `Aequorea victoria`, of `a priori` — labels
    that sit on PAYLOAD fragments every day. A reporter's alpha/omega digest
    had an insert half annotated with nothing but an `EcoRI` site, and it read
    as backbone. Once BOTH halves look like backbone, every consumer degrades
    differently and all of them are wrong: `_pick_insert_fragment` /
    `_pick_backbone_fragment` empty their candidate list and collapse to the
    raw SIZE heuristic the marker test exists to replace
    (`feedback_never_assume_smaller_frag_is_payload`);
    `_classify_part_from_plasmid` pass-4 requires `sum(marked) == 1` and bails
    to "no detectable grammar"; `_grammar_acceptor_tu_pairs` and cloning's
    `_entry_vector_acceptor_overhangs` fall back to "smallest", which INVERTS
    the `(oh5, oh3)` acceptor pair whenever a dropout outgrew its backbone.

    Two rules, both learned from real annotation text:

    * **Word start, not whole word.** `(?<![^\W\d_])` = "not preceded by a
      letter" (Unicode-aware — `[^\W\d_]` is the word class minus digits and
      underscore). A matching RIGHT boundary would be wrong: `oriT`, `oriV`,
      `oriP`, `ori2`, `oriR6K`, `origin` and `AmpR promoter` all continue past
      the keyword. Digits and separators are NOT letters, so `pMB1ori` and
      `rep_origin` still match.
    * **Suffix traps** for the words a keyword opens (see
      `_BACKBONE_LABEL_TRAPS`), because a word-start match can't tell `origin`
      from `original`.

    Keywords too short to survive even this (`cat`, `bla`, `tetr`, …) don't
    belong here at all — they go in `_BACKBONE_TOKEN_KEYWORDS`, which demands
    a whole token AND a gene-bearing feature type.
    """
    parts: list[str] = []
    for kw in keywords:
        if not kw:
            continue          # an empty branch would match EVERY label
        bad = traps.get(kw)
        parts.append(re.escape(kw)
                     + (f"(?!{'|'.join(b for b in bad if b)})" if bad else ""))
    if not parts:
        # `(?!)` never matches. An empty alternation would compile to a
        # pattern that matches the empty string at every position, i.e. every
        # fragment would read as backbone — fail CLOSED instead.
        return re.compile(r"(?!)")
    return re.compile(r"(?<![^\W\d_])(?:" + "|".join(parts) + r")",
                      re.IGNORECASE)


_BACKBONE_LABEL_RE = _label_keyword_re(_BACKBONE_LABEL_KEYWORDS,
                                       _BACKBONE_LABEL_TRAPS)


# Resistance genes whose names are too SHORT to substring-match safely.
# Checked as whole tokens (split on non-alphanumerics), never as substrings:
#   * "erm" inside "T-erm-inator" / "T-erm 908" — terminators are on nearly
#     every construct, so substring "erm" would mark practically any fragment
#     as backbone and silently disable the whole heuristic.
#   * "cat" inside "category", "bla" inside "blast", "neo" inside "neomycin
#     phosphotransferase-like" — the same false-positive class that pushed
#     `_detect_selection_marker` to token matching in sweep #34.
# ermB/ermC are the reason this exists (found 2026-08-05): a very common
# Gram-positive selection marker that BOTH this predicate and
# `_SELECTION_MARKER_KEYWORDS` (which only lists "ermr") failed to recognise,
# so a shuttle-vector backbone carrying it read as marker-free and the
# insert/backbone pick silently fell back to fragment size.
#
# Only consulted on feature types that actually ENCODE a gene. Verified
# against real annotated vectors: a genuine marker is `gene`/`CDS`, whereas
# the tokens leak into other types as references TO a gene rather than the
# gene itself —
#     regulatory   /label="as_annotated_pXX-ermB"   (a promoter, annotated
#                                                    after another plasmid)
#     primer_bind  /label="ermB-CLO-F"              (an oligo that anneals
#                                                    near ermB)
# Both tokenise to include "ermb" and would mark a payload fragment as
# backbone, which is worse than not knowing ermB at all: it makes BOTH halves
# look like backbone and collapses the pick back to raw fragment size. The
# long-form substring keywords above are unrestricted (unchanged behaviour) —
# they're specific enough not to need this guard.
_BACKBONE_TOKEN_FEATURE_TYPES: frozenset[str] = frozenset({"cds", "gene"})

_BACKBONE_TOKEN_KEYWORDS: frozenset[str] = frozenset({
    "erm", "ermb", "ermc", "ermam", "erma", "ermg",      # erythromycin / MLS
    "bla", "cat", "aada", "aac", "aph",                  # amp / cm / spec / gent
    "neo", "neor", "nptii", "kan",                       # kanamycin / neomycin
    "smr", "hyg", "hygr", "hph",                         # spec / hygromycin
    "zeo", "zeor", "ble", "gmr", "gent",                 # zeocin / gentamicin
    "tetm", "teta", "tetl",                              # tetracycline
    "tetr", "tetracyclin", "tetracycline",               # ← tet, see below
    "puror", "pac", "sh", "nat1",                        # puromycin / nourseo.
})
# The tet family is here rather than in `_BACKBONE_LABEL_KEYWORDS` for BOTH
# reasons this tier exists (moved 2026-08-18, adversarial corpus sweep):
#   * substring `"tetr"` is inside `tetr·atricopeptide repeat`, `tetr·amer`
#     -ization domain, `tetr·aspanin`, `Tetr·ahymena` ribozyme, `tetr·aloop` —
#     none of which have anything to do with tetracycline, and the first two
#     are among the most-annotated protein domains there are;
#   * tet is the ONE marker family whose repressor and operator are standard
#     PAYLOAD regulatory parts — `tetO`, `TRE`, `tetracycline-responsive
#     promoter`, `tetR binding site`, primers named `TetR-F`. A tet mention
#     outside a `CDS`/`gene` names the regulatory element, not the resistance
#     gene, which is exactly the distinction this tier already draws for ermB.
# A real `CDS`/`gene` labelled `TetR` or `tetracycline resistance protein`
# still matches, on the token, as before.


def _label_tokens(label: str) -> "set[str]":
    """Lowercased alphanumeric tokens of a feature label. Mirrors
    `_detect_selection_marker`'s tokeniser so both marker checks agree on
    what counts as a whole word."""
    return {t.lower() for t in re.split(r"[^A-Za-z0-9]+", label) if t}


def _feature_is_backbone_marker(ftype, label) -> bool:
    """Does ONE feature's ``(type, label)`` read as a backbone marker?

    THE matching rule, in one place. `_fragment_has_backbone_marker` and
    `_fragment_backbone_marker_labels` used to hand-keep two copies of it and
    the copies DID drift: the evidence list clamped labels to 80 characters
    BEFORE matching, so a marker past character 80 made the predicate say True
    while the list that is supposed to be its evidence came back empty (found
    + fixed 2026-08-18). One rule, no drift.

    Three tiers, narrowest first:
      1. the feature TYPE is `rep_origin` / `oriT` — unambiguous;
      2. a `_BACKBONE_LABEL_KEYWORDS` keyword at a WORD START, minus its
         suffix traps (`_label_keyword_re` carries the worked examples);
      3. a `_BACKBONE_TOKEN_KEYWORDS` name as a WHOLE TOKEN, and only on a
         gene-bearing feature type — the short names (`cat`, `bla`, `tetr`)
         that even a word-start match can't make safe.
    """
    ft = str(ftype or "").lower()
    if ft in _BACKBONE_FEATURE_TYPES:
        return True
    lbl = str(label or "")
    if not lbl:
        return False
    low = lbl.lower()
    # Cheap pre-filter, then the real rule. `_BACKBONE_LABEL_RE` only ever
    # matches where the bare keyword occurs (it adds boundary + trap
    # restrictions, never new matches), so a plain substring scan is a SOUND
    # gate — and C-level `in` beats running an alternation with a lookbehind
    # over every label. Measured on a 20-feature fragment 2026-08-18: 44 µs
    # regex-always → 17 µs gated, i.e. back to the pre-boundary cost.
    # `test_word_start_rule_implies_substring` pins the soundness.
    if (any(kw in low for kw in _BACKBONE_LABEL_KEYWORDS)
            and _BACKBONE_LABEL_RE.search(lbl)):
        return True
    return (ft in _BACKBONE_TOKEN_FEATURE_TYPES
            and bool(_label_tokens(lbl) & _BACKBONE_TOKEN_KEYWORDS))


def _fragment_features(frag) -> list:
    """``frag``'s feature dicts, or ``[]`` for anything malformed.

    Both marker helpers are called on digest output, on parsed-GenBank
    output, and (through the agent API) on caller-supplied JSON, so a
    non-dict fragment or a `features` value that isn't a list must return
    "no markers" rather than raising out of a Constructor worker."""
    if not isinstance(frag, dict):
        return []
    feats = frag.get("features")
    # The list case is returned as-is, NOT copied: callers only iterate, and
    # these helpers sit in the classifier's inner loop. A tuple is converted
    # (rare, and it keeps the return type honest).
    if isinstance(feats, list):
        return feats
    return list(feats) if isinstance(feats, tuple) else []


def _fragment_has_backbone_marker(frag: dict) -> bool:
    """Return True iff ``frag``'s features include a typical
    bacterial-backbone marker (origin of replication or antibiotic
    resistance) — see `_feature_is_backbone_marker` for the match rule.

    Used by `_pick_insert_fragment` to avoid the "smallest fragment
    is the dropout" heuristic — that rule breaks the moment a
    stacked-TU/MOD insert outgrows its carrier vector. Looking for
    the ORIGIN/SELECTION markers is reliable because real Golden
    Braid / MoClo entry vectors annotate them, and the L0 parts
    chained INTO an insert never do.

    Biased toward false NEGATIVES: a missed marker falls back to the size
    heuristic and logs it, while a false positive silently picks the wrong
    half of a digest (see `_label_keyword_re`)."""
    for f in _fragment_features(frag):
        if isinstance(f, dict) and _feature_is_backbone_marker(
                f.get("type"), f.get("label")):
            return True
    return False


def _fragment_backbone_marker_labels(frag: dict) -> "list[str]":
    """The feature names that make `_fragment_has_backbone_marker` say True,
    de-duplicated in first-seen order (empty list when it says False).

    Same rule — literally the same function — so the two can never disagree;
    this only reports WHICH features tripped it. Lets a caller say "carries
    AmpR, Ori*", which reads as a diagnosis, instead of "carries a backbone
    marker", which reads as noise the user has no way to act on."""
    out: list[str] = []
    for f in _fragment_features(frag):
        if not isinstance(f, dict):
            continue
        if not _feature_is_backbone_marker(f.get("type"), f.get("label")):
            continue
        # Clamp for DISPLAY only (never before matching): a label is
        # user-controlled and can be arbitrarily long — imported GenBank
        # happily carries multi-line /note text — and these names go straight
        # into a one-line warning.
        name = (str(f.get("label") or "").strip()[:80]
                or str(f.get("type") or "").lower())
        if name and name not in out:
            out.append(name)
    return out


def _ev_frag_input_features(record) -> "list[dict]":
    """Marshal a parsed entry-vector record's features into the
    ``[{start, end, type, label}]`` shape `_excise_fragment_pair` splits onto
    its fragments, so `_fragment_has_backbone_marker` can tell the backbone half
    from the dropout by ORIGIN/RESISTANCE markers instead of by size (which
    inverts when a vector's dropout cassette outgrows its backbone)."""
    feats: "list[dict]" = []
    for f in (getattr(record, "features", None) or []):
        try:
            start = int(f.location.start)
            end = int(f.location.end)
        except (AttributeError, TypeError, ValueError):
            continue
        q = getattr(f, "qualifiers", None) or {}
        label = ""
        for k in ("label", "gene", "product", "note"):
            v = q.get(k)
            if v:
                label = str(v[0] if isinstance(v, (list, tuple)) else v)
                break
        feats.append({
            "start": start, "end": end,
            "type": str(getattr(f, "type", "") or ""),
            "label": label,
        })
    return feats


_VECTOR_MATCH_CACHE_MAX = 64


def _vector_half_top_seq(ev_gb: str, enzyme: str) -> "str | None":
    """Cached helper: return the entry vector's vector-half top_seq
    after digesting with ``enzyme``. Returns None on parse / digest
    failure or when the EV is non-circular (digest needs a ring).
    Cache is keyed by ``(ev_gb, enzyme)`` so a saved EV's digest is
    re-used across `_check_vector_match` calls within a classification
    run.

    Non-circular EVs are explicitly skipped (returns None) — without
    the topology check, a linearised EV file would dispatch through
    `_excise_fragment_pair(circular=True)` and produce nonsense
    fragments. Better to fail closed than report a phantom match.

    Sweep #25 (2026-05-23): key is `(hash(ev_gb), enzyme)` not
    `(ev_gb, enzyme)`. Pre-fix the key held the full gb_text string
    (potentially multi-MB); at 64 entries × 5 MB worst-case
    a steady-state cache could pin ~320 MB just in tuple keys.
    `hash(str)` is deterministic per-process so cache hits work
    correctly within a session. Collisions are functionally
    indistinguishable from a miss (we'd just recompute) — Python
    `hash(str)` collisions are vanishingly rare across the small
    enzyme inputs we feed.
    """
    key = (hash(ev_gb), enzyme)
    with _state._COMPUTE_CACHE_LOCK:
        if key in _state._VECTOR_MATCH_CACHE:
            return _state._VECTOR_MATCH_CACHE[key]
    result: "str | None" = None
    try:
        ev_rec = _gb_text_to_record(ev_gb)
        topology = (
            getattr(ev_rec, "annotations", {}) or {}
        ).get("topology", "")
        if str(topology).lower() != "circular":
            # Refuse linearised entry vectors — `_excise_fragment_pair`
            # with `circular=True` on a linear sequence misreports
            # cuts at the joined ends.
            _log.debug(
                "_vector_half_top_seq: skipping non-circular EV "
                "(topology=%r)", topology,
            )
        else:
            ev_seq = str(ev_rec.seq).upper()
            ev_frags, ev_err = _excise_fragment_pair(
                ev_seq, [enzyme], circular=True,
                features=_ev_frag_input_features(ev_rec),
            )
            if ev_err is None and len(ev_frags) == 2:
                # The VECTOR half carries the backbone marker (ori/resistance);
                # prefer it by marker, not size — size inverts when a vector's
                # dropout outgrows its backbone. Fall back to the larger fragment
                # when the marker signal is absent/ambiguous.
                marked = [f for f in ev_frags
                          if _fragment_has_backbone_marker(f)]
                ev_vector_half = (
                    marked[0] if len(marked) == 1
                    else max(ev_frags,
                             key=lambda f: len(f.get("top_seq") or ""))
                )
                result = (
                    ev_vector_half.get("top_seq") or ""
                ).upper() or None
    except Exception:
        _log.debug(
            "_vector_half_top_seq: digest failed for enzyme %r", enzyme,
        )
    # Bounded LRU-ish: drop one arbitrary entry when over cap. Insertion
    # order eviction (Python 3.7+ dict preserves insertion order) gives
    # us FIFO behaviour without an OrderedDict import.
    with _state._COMPUTE_CACHE_LOCK:
        if len(_state._VECTOR_MATCH_CACHE) >= _VECTOR_MATCH_CACHE_MAX:
            try:
                _state._VECTOR_MATCH_CACHE.pop(next(iter(_state._VECTOR_MATCH_CACHE)))
            except (StopIteration, KeyError):
                pass
        _state._VECTOR_MATCH_CACHE[key] = result
    return result


def _check_vector_match(
    gid: str, enzyme: str, user_vector_frag: dict,
) -> "dict | None":
    """Compare the user's vector half (the larger digest fragment)
    against every entry vector configured for grammar ``gid``. Returns
    ``{role, name, matches: True}`` for the first vector half that's
    rotationally identical (in either orientation), or ``None`` when
    no entry vector is configured for ``gid`` / no match is found.

    The comparison uses the same ``enzyme`` that produced the user's
    digest so the two vector halves are directly comparable. Rotation-
    invariance is handled by checking whether the user's top_seq is a
    substring of the entry vector's doubled top_seq — fragments from
    rotationally equivalent rings have the same content but may start
    at different positions in the linearised representation. The same
    test is run against the EV's reverse-complement-doubled seq so an
    RC-saved user plasmid (the same biological ring written with the
    other strand on top) still matches.

    Used by `_classify_part_from_plasmid` to surface "this plasmid was
    cloned into your configured Alpha1 entry vector" so Load Part can
    confirm the user's expected destination matches before saving.
    """
    user_vec_seq = (user_vector_frag.get("top_seq") or "").upper()
    if not user_vec_seq:
        return None
    for ev in _load_entry_vectors():
        if ev.get("grammar_id") != gid:
            continue
        ev_gb = ev.get("gb_text") or ""
        if not ev_gb:
            continue
        ev_vec_seq = _vector_half_top_seq(ev_gb, enzyme)
        if not ev_vec_seq or len(user_vec_seq) != len(ev_vec_seq):
            continue
        if user_vec_seq in (ev_vec_seq + ev_vec_seq):
            return {
                "role":    ev.get("role") or "",
                "name":    ev.get("name", "?"),
                "matches": True,
            }
        # Reverse-strand orientation: the user may have saved the
        # plasmid with the other strand on top, which makes the
        # digest's top_seq the RC of the canonical vector half. Try
        # the RC-doubled form too — same rotation-invariant check.
        try:
            rc_doubled = _rc(ev_vec_seq)
            if rc_doubled and user_vec_seq in (rc_doubled + rc_doubled):
                return {
                    "role":    ev.get("role") or "",
                    "name":    ev.get("name", "?"),
                    "matches": True,
                }
        except Exception:
            _log.debug(
                "_check_vector_match: RC fallback failed for ev %r",
                ev.get("name"),
            )
            continue
    return None


_ACCEPTOR_TU_PAIRS_CACHE_MAX = 64


# ── Synthesis / assembly readiness lint ─────────────────────────────────────
# The "is it safe to order + assemble?" pre-flight: internal Type IIS sites
# (kill Golden Gate / MoClo), GC extremes + homopolymer / tandem runs (kill
# gene synthesis), degenerate bases, and optional ORF integrity. Read-only;
# composes `_enzyme_cuts` (biology) + `_find_orfs` with a few local scans.
_SYNTHESIS_FORBIDDEN_ENZYMES = ("BsaI", "BsmBI", "BbsI", "SapI", "Esp3I")


def _lint_homopolymer_runs(seq: str, min_len: int
                           ) -> "list[tuple[int, str, int]]":
    """``[(start, base, length), …]`` for each single-base run ≥ ``min_len``."""
    out: "list[tuple[int, str, int]]" = []
    n, i = len(seq), 0
    while i < n:
        j = i + 1
        while j < n and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            out.append((i, seq[i], j - i))
        i = j
    return out


def _lint_tandem_repeats(seq: str, *, unit_min: int, unit_max: int,
                         min_copies: int) -> "list[tuple[int, str, int]]":
    """``[(start, unit, copies), …]`` for tandem k-mer repeats (``unit_min ≤ k
    ≤ unit_max``) of ≥ ``min_copies`` copies. Greedy + non-overlapping;
    single-base units are skipped (homopolymers are reported separately)."""
    out: "list[tuple[int, str, int]]" = []
    n, i = len(seq), 0
    while i < n:
        best = None
        for u in range(unit_min, unit_max + 1):
            if i + u > n:
                break
            unit = seq[i:i + u]
            if len(set(unit)) == 1:
                continue
            copies, j = 1, i + u
            while j + u <= n and seq[j:j + u] == unit:
                copies, j = copies + 1, j + u
            if copies >= min_copies and (best is None
                                          or copies * u > best[2] * len(best[1])):
                best = (i, unit, copies)
        if best is not None:
            out.append(best)
            i += best[2] * len(best[1])
        else:
            i += 1
    return out


def _lint_gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def _lint_worst_gc_window(seq: str, window: int) -> "tuple[int, float] | None":
    """``(start, gc_fraction)`` of the window whose GC is furthest from 0.5, or
    None when the sequence is shorter than ``window``. O(n) sliding sum."""
    n = len(seq)
    if window <= 0 or n < window:
        return None
    cur = sum(1 for c in seq[:window] if c in "GC")
    worst_i, worst_gc = 0, cur / window
    for i in range(1, n - window + 1):
        if seq[i - 1] in "GC":
            cur -= 1
        if seq[i + window - 1] in "GC":
            cur += 1
        frac = cur / window
        if abs(frac - 0.5) > abs(worst_gc - 0.5):
            worst_i, worst_gc = i, frac
    return worst_i, worst_gc


def _synthesis_lint(seq, *, circular: bool = True, expect_cds: bool = False,
                    forbidden_enzymes=_SYNTHESIS_FORBIDDEN_ENZYMES,
                    min_gc: float = 0.30, max_gc: float = 0.70,
                    gc_window: int = 50, gc_window_lo: float = 0.20,
                    gc_window_hi: float = 0.80, homopolymer_max: int = 9,
                    repeat_unit_max: int = 4, repeat_min_copies: int = 5,
                    max_warnings_per_kind: int = 20) -> dict:
    """Synthesis / assembly readiness lint. Returns
    ``{"score": int, "warnings": [...], "stats": {...}}``. Read-only.

    ``warnings`` are ``{level, kind, message, start?, end?}`` with
    ``level ∈ {"error","warn","info"}``. Checks: internal Type IIS sites
    (assembly killers → error), extreme overall + windowed GC, long homopolymer
    runs, tandem repeats, degenerate/non-ACGT bases, and — when ``expect_cds``
    — a clean full-length ORF. ``score`` starts at 100; each finding deducts by
    severity (error −20, warn −8, info −2, floored at 0)."""
    raw = (seq or "").upper()
    s = "".join(c for c in raw if c in "ACGT")
    warnings: "list[dict]" = []

    def add(level, kind, message, start=None, end=None):
        w = {"level": level, "kind": kind, "message": message}
        if start is not None:
            w["start"] = int(start)
        if end is not None:
            w["end"] = int(end)
        warnings.append(w)

    n_degenerate = len(raw) - len(s)
    if n_degenerate:
        add("warn", "degenerate_bases",
            f"{n_degenerate} non-ACGT / ambiguous base(s) — most synthesis "
            "vendors reject degenerate positions")

    catalog = _state._all_enzymes_hook()
    for enz in forbidden_enzymes:
        if enz not in catalog:
            continue
        try:
            cuts = _enzyme_cuts(s, [enz], circular=circular)
        except Exception:
            continue
        for c in cuts[:max_warnings_per_kind]:
            add("error", "type_iis_site",
                f"internal {enz} site — cut during Golden Gate / MoClo assembly "
                "(domesticate it out)", start=c.get("top"))
        if len(cuts) > max_warnings_per_kind:
            add("info", "type_iis_site",
                f"…and {len(cuts) - max_warnings_per_kind} more {enz} site(s)")

    gc = _lint_gc(s)
    if s and (gc < min_gc or gc > max_gc):
        add("warn", "gc_overall",
            f"overall GC {gc*100:.0f}% is outside the {min_gc*100:.0f}–"
            f"{max_gc*100:.0f}% synthesis-friendly range")
    ww = _lint_worst_gc_window(s, gc_window)
    if ww is not None and (ww[1] < gc_window_lo or ww[1] > gc_window_hi):
        add("warn", "gc_window",
            f"a {gc_window} bp window at +{ww[0]} is {ww[1]*100:.0f}% GC — local "
            "extremes hurt synthesis + PCR", start=ww[0], end=ww[0] + gc_window)

    for start, base, length in _lint_homopolymer_runs(
            s, homopolymer_max)[:max_warnings_per_kind]:
        add("warn", "homopolymer",
            f"{length}×{base} homopolymer run — synthesis-hard (slippage)",
            start=start, end=start + length)

    for start, unit, copies in _lint_tandem_repeats(
            s, unit_min=2, unit_max=repeat_unit_max,
            min_copies=repeat_min_copies)[:max_warnings_per_kind]:
        add("warn", "tandem_repeat",
            f"{copies}×({unit}) tandem repeat — synthesis + assembly hazard",
            start=start, end=start + copies * len(unit))

    orf_count = 0
    if s:
        try:
            orfs = _find_orfs(s, circular=circular)
        except Exception:
            orfs = []
        orf_count = len(orfs)
        if expect_cds:
            def _span(o):
                # Read the reported coding length, NOT the coordinate pair:
                # on a circle `(start, end)` cannot express an ORF of a full
                # lap or more, which is exactly the "spans ~the whole
                # sequence" case this check is looking for. Coord fallback
                # keeps a hand-built ORF dict working.
                nt = o.get("nt_len")
                if nt is not None:
                    return int(nt)
                st, en = int(o["start"]), int(o["end"])
                return (len(s) - st) + en if en < st else en - st
            if not any(_span(o) >= 0.9 * len(s) for o in orfs):
                add("error", "orf_integrity",
                    "no full-length ORF spanning ~the whole sequence — a "
                    "premature stop or frame shift would truncate the protein")

    penalty = sum({"error": 20, "warn": 8, "info": 2}.get(w["level"], 0)
                  for w in warnings)
    return {
        "score": max(0, 100 - penalty),
        "warnings": warnings,
        "stats": {"length": len(s), "gc": round(gc, 4),
                  "degenerate_bases": n_degenerate, "orf_count": orf_count,
                  "circular": bool(circular)},
    }


def _grammar_acceptor_tu_pairs(
    grammar_id: str, enzyme: str,
) -> "list[tuple[str, str, str, str]]":
    """Return ``[(role, ev_name, oh5, oh3), ...]`` — the stuffer's
    overhang pair released by digesting each configured entry vector
    for ``grammar_id`` with ``enzyme``.

    The L0 → L1 entry vector for Golden Braid has four canonical
    roles (Alpha1 / Alpha2 / Omega1 / Omega2), each receiving a TU
    in a different orientation. When BsaI digests the assembled TU
    plasmid (a TU INSIDE one of those acceptors), the released
    insert carries that acceptor's specific overhang pair — which
    is the SAME pair as the stuffer's overhangs in the empty
    acceptor itself. Pre-2026-05-13 the classifier only knew the
    canonical (Promoter.oh5, Terminator.oh3) pair via
    `_grammar_tu_overhangs`, so a TU in Alpha2 / Omega1 / Omega2
    silently failed to classify (the bug a user hit on DEMO-25 in
    alpha-2).

    Singleton entry vectors (role == "") are skipped — those are L0
    acceptors (pUPD2 et al.), not TU acceptors, and the L0 position
    table check upstream covers them.

    Result is cached per ``(grammar_id, enzyme)`` and invalidated by
    `_save_entry_vectors`. A multi-select Load Part batch that
    classifies 50 plasmids against 2 grammars × 2 enzymes was
    re-digesting every configured EV 200× without this cache.

    Failures (no gb_text, parse error, digest yields ≠ 2 fragments)
    are logged at warning level so a misconfigured EV surfaces in
    the diagnostic bundle — the user otherwise sees their TU silently
    fall through to None classification with no hint why.
    """
    cache_key = (grammar_id, enzyme)
    with _state._COMPUTE_CACHE_LOCK:
        cached = _state._ACCEPTOR_TU_PAIRS_CACHE.get(cache_key)
        if cached is not None:
            return list(cached)
    out: "list[tuple[str, str, str, str]]" = []
    for ev in _load_entry_vectors():
        if ev.get("grammar_id") != grammar_id:
            continue
        role = ev.get("role") or ""
        if not role:
            continue   # singleton L0 vector — not a TU acceptor
        ev_gb = ev.get("gb_text") or ""
        ev_name = ev.get("name") or "?"
        if not ev_gb:
            _log.warning(
                "_grammar_acceptor_tu_pairs: entry vector %r "
                "(role=%r, grammar=%r) has no gb_text — skipping",
                ev_name, role, grammar_id,
            )
            continue
        try:
            record = _gb_text_to_record(ev_gb)
            ev_seq = str(getattr(record, "seq", "") or "").upper()
            if not ev_seq:
                _log.warning(
                    "_grammar_acceptor_tu_pairs: entry vector %r "
                    "parsed to empty sequence — skipping",
                    ev_name,
                )
                continue
            frags, err = _excise_fragment_pair(
                ev_seq, [enzyme], circular=True,
                features=_ev_frag_input_features(record),
            )
        except Exception:
            _log.exception(
                "_grammar_acceptor_tu_pairs: digest failed for "
                "ev=%r role=%r enzyme=%r", ev_name, role, enzyme,
            )
            continue
        if err is not None or len(frags) != 2:
            _log.info(
                "_grammar_acceptor_tu_pairs: ev=%r role=%r digest "
                "with %r produced %d fragment(s) (err=%r) — skipped",
                ev_name, role, enzyme, len(frags), err,
            )
            continue
        # Stuffer/dropout = the placeholder (lacZα / ccdB) that the assembled
        # TU replaces; its 5'/3' overhangs ARE this acceptor's TU-boundary
        # overhangs. Identify it by marker-ABSENCE (not "smaller") so an
        # acceptor whose dropout outgrew its backbone still yields the right
        # overhangs; fall back to the smaller fragment when marker-ambiguous.
        non_marker = [f for f in frags
                      if not _fragment_has_backbone_marker(f)]
        insert = (non_marker[0] if len(non_marker) == 1
                  else min(frags, key=lambda f: len(f.get("top_seq", ""))))
        oh5 = (insert.get("left")  or {}).get(
            "overhang_seq", "",
        ).upper()
        oh3 = (insert.get("right") or {}).get(
            "overhang_seq", "",
        ).upper()
        if oh5 and oh3:
            out.append((role, ev_name, oh5, oh3))
    # Sweep #26: FIFO-evict oldest entry if at cap.
    with _state._COMPUTE_CACHE_LOCK:
        if len(_state._ACCEPTOR_TU_PAIRS_CACHE) >= _ACCEPTOR_TU_PAIRS_CACHE_MAX:
            try:
                _state._ACCEPTOR_TU_PAIRS_CACHE.pop(next(iter(_state._ACCEPTOR_TU_PAIRS_CACHE)))
            except (StopIteration, KeyError):
                pass
        _state._ACCEPTOR_TU_PAIRS_CACHE[cache_key] = list(out)
    return out


def _classify_part_from_plasmid(
    seq: str,
    *,
    circular: bool,
    features: "list[dict] | None" = None,
) -> "dict | None":
    """Identify which cloning grammar + position + level a circular
    plasmid holds, by trying each grammar's primary and secondary
    Type IIS enzymes and matching the released fragment's overhangs
    against either the grammar's L0 position table or its TU boundary
    overhangs (`_grammar_tu_overhangs`).

    SACRED INVARIANT: the (oh5, oh3) overhang pair released by the
    digest is the **only** input used to determine position type.
    The classifier never overrides this from feature labels, plasmid
    name, source filename, or any other heuristic — the user's
    biological molecule has ONE legal position per overhang pair, so
    the lookup is unambiguous and the code path must reflect that.
    Callers can re-tag manually via the Parts Bin Edit modal if they
    really need to, but the classifier itself stays mechanical.

    Detection cases (per grammar, in registry order, for each
    enzyme; **both fragments** are tried in each pass so a library
    entry without backbone-marker annotations still classifies
    correctly when the insert outgrew the backbone):

      * Either fragment's overhangs match an L0 position → ``level=0``
        (L0 part, regardless of which side of the cycle produced it).
      * Either fragment's overhangs match the canonical TU boundary
        OR a configured entry vector's stuffer pair → ``level=1`` (TU).

    Pre-2026-05-13 the classifier picked a single "insert" fragment
    via `_pick_insert_fragment` and inferred level from enzyme
    parity ("primary release ⇒ MOD"). Both heuristics broke for
    real-world plasmids:

      * `_pick_insert_fragment` falls through to "smallest fragment"
        when no `rep_origin`/`selection_marker` features are present;
        for a TU whose body is larger than its alpha backbone
        (common for any cassette ~2 kbp+) the backbone half got
        picked as the "insert" and its mirrored overhangs matched
        nothing. DEMO 26 in the DemoColl collection: 3250 bp body with
        the correct (GGAG, GTCA) overhangs, 1850 bp backbone with
        the mirrored (GTCA, GGAG) — the body was discarded.

      * Enzyme parity assumed the splicecraft convention
        (Esp3I = primary = L0 release, BsaI = secondary = L1
        release). The actual pDGB1 / GB 2.0 convention used by
        Sarrion-Perdigones 2013 and the user's DemoColl collection has
        these REVERSED (BsaI = L0 release from pUPD2, Esp3I = L1
        release from α-vectors). Under that convention, the user's
        TUs were being mis-classified as MOD (level=2).

    Auto MOD (level=2) detection from overhangs alone is unreliable
    across both conventions — a fragment with TU-boundary overhangs
    could be a TU or a MOD depending on which enzyme cycle the lab
    uses, and we can't tell from the overhangs alone. TUs that are
    actually MODs in the user's lab can be re-tagged via the
    Parts Bin Edit modal.

    Returns ``None`` when no grammar gives a clean (exactly 2-fragment)
    digest with recognised overhangs. Otherwise returns
    ``{grammar_id, grammar_name, level, position, insert, vector,
       release_enzyme, entry_vector}``.

    Linear records skip the digest — a linear "part" can't be cleanly
    excised, and the overhangs in the .gb annotation are the source
    of truth in that case.
    """
    if not seq or not circular:
        return None
    grammars = _all_grammars()
    for gid, grammar in grammars.items():
        primary   = grammar.get("enzyme")
        secondary = (grammar.get("level_up_enzyme")
                     or grammar.get("enzyme"))
        for enzyme_role, enzyme in (
            ("primary",   primary),
            ("secondary", secondary),
        ):
            if not isinstance(enzyme, str) or not enzyme:
                continue
            # Skip the secondary pass entirely when it's the same as
            # the primary (custom grammar that omits `level_up_enzyme`)
            # — we'd just re-run the same digest with the same outcome.
            if enzyme_role == "secondary" and enzyme == primary:
                continue
            try:
                frags, err = _excise_fragment_pair(
                    seq, [enzyme], circular=True,
                    features=features,
                    source_label=grammar.get("name", gid),
                )
            except Exception:
                _log.exception(
                    "_classify_part_from_plasmid: %s digest failed for %s",
                    enzyme_role, gid,
                )
                continue
            if err is not None or len(frags) != 2:
                continue
            # Try BOTH fragments — sized assumptions about which half
            # is the "insert" break when the user's library entries
            # don't carry rep_origin/selection-marker annotations
            # AND the assembled cassette is larger than its carrier
            # backbone (DEMO 26 family, 2026-05-13). For each fragment,
            # match its overhangs against the L0 position table, the
            # canonical TU boundary, and the per-acceptor stuffer
            # pairs; the FIRST match wins, with the other fragment
            # becoming the vector half by elimination.
            for insert_idx, insert in enumerate(frags):
                vector = frags[1 - insert_idx]
                oh5 = (insert.get("left")  or {}).get(
                    "overhang_seq", "",
                ).upper()
                oh3 = (insert.get("right") or {}).get(
                    "overhang_seq", "",
                ).upper()
                if not oh5 or not oh3:
                    continue
                # First check: the L0 position table. Both pre- and
                # post-cloning L0 parts land here — only the enzyme
                # differs depending on the lab's GB convention.
                for pos in (grammar.get("positions") or []):
                    pos_oh5 = str(pos.get("oh5", "")).upper()
                    pos_oh3 = str(pos.get("oh3", "")).upper()
                    if pos_oh5 == oh5 and pos_oh3 == oh3:
                        return {
                            "grammar_id":     gid,
                            "grammar_name":   grammar.get("name", gid),
                            "level":          0,
                            "position":       dict(pos),
                            "insert":         insert,
                            "vector":         vector,
                            "release_enzyme": enzyme,
                            "entry_vector":   _check_vector_match(
                                gid, enzyme, vector,
                            ),
                        }
                # Second check: canonical TU boundary overhangs
                # (Pos 1's oh5 + last position's oh3). Matches a TU
                # assembled into the canonical L1 acceptor — the
                # orientation that lines up with the grammar's L0
                # positions.
                tu_start, tu_end = _grammar_tu_overhangs(grammar)
                if (tu_start and tu_end
                        and tu_start.upper() == oh5
                        and tu_end.upper()   == oh3):
                    position = {
                        "name":  _part_level_label(1),
                        "type":  _part_level_label(1),
                        "oh5":   tu_start,
                        "oh3":   tu_end,
                        "color": "white",
                    }
                    return {
                        "grammar_id":     gid,
                        "grammar_name":   grammar.get("name", gid),
                        "level":          1,
                        "position":       position,
                        "insert":         insert,
                        "vector":         vector,
                        "release_enzyme": enzyme,
                        "entry_vector":   _check_vector_match(
                            gid, enzyme, vector,
                        ),
                    }
                # Third check: per-acceptor TU boundary overhangs.
                # Each configured entry vector's stuffer carries a
                # specific (oh5, oh3) pair; a TU assembled into
                # that acceptor will release with the same pair.
                # DEMO-25 in alpha-2 hit this pass (2026-05-13).
                for role, ev_name, acc_oh5, acc_oh3 in (
                    _grammar_acceptor_tu_pairs(gid, enzyme)
                ):
                    if acc_oh5 == oh5 and acc_oh3 == oh3:
                        label = _part_level_label(1)
                        position = {
                            "name":  f"{label} ({role})",
                            "type":  label,
                            "oh5":   acc_oh5,
                            "oh3":   acc_oh3,
                            "color": "white",
                        }
                        return {
                            "grammar_id":     gid,
                            "grammar_name":   grammar.get("name", gid),
                            "level":          1,
                            "position":       position,
                            "insert":         insert,
                            "vector":         vector,
                            "release_enzyme": enzyme,
                            "entry_vector":   _check_vector_match(
                                gid, enzyme, vector,
                            ),
                        }
            # Pass 4 (2026-05-20): lenient TU detection for L1
            # acceptors that release outside the canonical (Pos 1
            # oh5, Pos N oh3) pair but use overhangs from the
            # grammar's own canonical alphabet. Caught by user
            # report on the DemoColl collection DEMO 25-31 (TUx1/TUx2
            # acceptors) — TUx1 releases the TU with (GGAG, GTCA);
            # TUx2 with (GTCA, CGCT). Both pairs are valid GB 2.0
            # overhangs (GTCA = RC(TGAC), a Pos 1b operator
            # overhang in the GB 2.0 expanded grammar) but neither
            # matches the strict TU boundary and neither lab has
            # entry vectors configured in `entry_vectors.json`.
            #
            # Pick the TU candidate via `_pick_insert_fragment`'s
            # backbone-marker exclusion — NEVER by fragment size.
            # Size is unreliable for L1+ TUs (an assembled cassette
            # commonly outgrows the alpha-vector backbone), and
            # the user has explicitly called out this assumption
            # as a class of bugs (2026-05-20).
            #
            # Skip pass-4 when no fragment carries a clear backbone
            # marker — pass-4 cannot safely commit to one fragment
            # being the TU without that biology signal, and a
            # wrong-tag is worse than the existing "no detectable
            # grammar" outcome.
            backbone_marked = [
                _fragment_has_backbone_marker(f) for f in frags
            ]
            if sum(backbone_marked) != 1:
                # Either both have markers (annotation noise) or
                # neither does (un-annotated entry). Either way,
                # pass-4 can't safely commit. Fall through to the
                # next grammar/enzyme.
                continue
            tu_idx = 0 if not backbone_marked[0] else 1
            tu_candidate = frags[tu_idx]
            vector_candidate = frags[1 - tu_idx]
            oh5 = (tu_candidate.get("left") or {}).get(
                "overhang_seq", "",
            ).upper()
            oh3 = (tu_candidate.get("right") or {}).get(
                "overhang_seq", "",
            ).upper()
            if not oh5 or not oh3:
                continue
            canonical = _grammar_canonical_overhangs(grammar)
            if oh5 not in canonical or oh3 not in canonical:
                continue
            position = {
                "name":  f"TU ({oh5}→{oh3})",
                "type":  _part_level_label(1),
                "oh5":   oh5,
                "oh3":   oh3,
                "color": "white",
            }
            _log_event(
                "classify.lenient_tu",
                grammar=gid, enzyme=enzyme,
                oh5=oh5, oh3=oh3,
            )
            return {
                "grammar_id":     gid,
                "grammar_name":   grammar.get("name", gid),
                "level":          1,
                "position":       position,
                "insert":         tu_candidate,
                "vector":         vector_candidate,
                "release_enzyme": enzyme,
                "entry_vector":   _check_vector_match(
                    gid, enzyme, vector_candidate,
                ),
                "lenient":        True,
            }
    return None


def _grammar_canonical_overhangs(grammar: dict) -> "set[str]":
    """Set of every 4-bp overhang the grammar's position table knows,
    plus each one's reverse complement.

    Pre-2026-05-20 the classifier hard-coded the canonical TU
    boundary as ``(positions[0].oh5, positions[-1].oh3)`` — for
    `gb_l0` that's ``(GGAG, CGCT)``. Real Golden Braid α/Ω
    acceptors release TUs with NON-canonical pairs drawn from the
    same overhang alphabet (e.g. pDGB3 α1 releases with
    ``(GGAG, GTCA)`` — the second overhang is RC of TGAC, a Pos 1b
    operator from the GB 2.0 expanded grammar). This helper
    exposes the full alphabet so the lenient pass-4 of
    `_classify_part_from_plasmid` can recognise these without
    requiring the user to pre-configure entry vectors.

    Includes reverse complements because Type IIS overhangs read
    differently depending on which strand the cut surface is being
    described from — a released fragment whose 5' overhang on the
    top strand is ``GTCA`` has ``TGAC`` on the bottom strand, and
    both are equally "this overhang".
    """
    out: "set[str]" = set()
    for pos in (grammar.get("positions") or []):
        if not isinstance(pos, dict):
            continue
        for key in ("oh5", "oh3"):
            v = str(pos.get(key, "") or "").upper()
            if v and len(v) == 4 and all(c in "ACGT" for c in v):
                out.add(v)
                out.add(_rc(v))
    return out


def _part_level_label(level: int) -> str:
    """Map an integer level to the user-facing label used in the Parts
    Bin tabs and notify strings: L0 → 'L0', 1 → 'TU', ≥2 → 'MOD'."""
    if level <= 0:
        return "L0"
    if level == 1:
        return "TU"
    return "MOD"


# ── Spliced-transcript prediction ──────────────────────────────────────────
#
# The read-direction counterpart to `splicecraft_cassette`'s assembler: given a
# plasmid that ALREADY carries a transcription unit, reconstruct the mature
# mRNA and score translation initiation on it.
#
# It exists because of a specific, repeatable wrong answer. When a construct
# carries a genomic 5'UTR intron, the biology depends on the SPLICED leader,
# but every check written against the plasmid sequence sees the UNSPLICED one.
# An upstream-ATG screen run that way flags intronic ATGs that splicing
# removes — false positives that look exactly like real uORF hazards, and that
# each cost a round of investigation. Reconstructing the mature message by hand
# per script is how that mistake gets made repeatedly, so it lives here once.
#
# Introns come from ANNOTATION — `intron` features, or the gaps between the
# parts of a spliced CDS/mRNA location — never from de novo prediction.
# Guessing where a spliceosome cuts and then reporting the guess as the message
# would be the same class of failure this module was written to remove. The
# trained PWM in `splicecraft_splice` is still used, but only to flag CRYPTIC
# sites as an advisory alongside the answer, never to define the transcript.

_TX_MAX_UORFS = 200
_TX_MAX_SPLICE_HITS = 50


def _kozak_context(mrna: str, atg_pos: int) -> dict:
    """Score the translation-initiation context around `atg_pos` in `mrna`.

    Returns ``{context, minus3, plus4, strength, score, truncated, model}``.
    `context` is the 10 nt window ``[-6 .. +4]`` with the A of the ATG at
    offset 6, so it reads the way the literature prints it.

    **This is the positional consensus rule, not a trained matrix** — the two
    positions that carry most of the signal are a purine at −3 and a G at +4,
    and `score` just counts how many of those two are satisfied (0–2 →
    weak / adequate / strong). It is reported as `model` on every result so a
    caller never mistakes it for the calibrated PWM that
    `splicecraft_splice` uses for splice sites. A context that runs off the
    5' end of the message is marked `truncated` — a 5'UTR shorter than 6 nt
    is itself a finding, and padding it would invent bases."""
    n = len(mrna)
    lo = atg_pos - 6
    truncated = lo < 0 or atg_pos + 4 > n
    context = mrna[max(0, lo):atg_pos + 4]
    minus3 = mrna[atg_pos - 3] if atg_pos >= 3 else ""
    plus4 = mrna[atg_pos + 3] if atg_pos + 4 <= n else ""
    score = int(minus3 in ("A", "G")) + int(plus4 == "G")
    return {
        "context":   context,
        "minus3":    minus3,
        "plus4":     plus4,
        "score":     score,
        "strength":  ("strong" if score == 2
                      else "adequate" if score == 1 else "weak"),
        "truncated": truncated,
        "model":     "positional -3/+4 consensus (not a trained PWM)",
    }


def _tx_map_span(gs: int, ge: int, a: int, total: int,
                 strand: int, pre_len: int) -> "tuple[int, int] | None":
    """Map a genomic half-open span onto pre-mRNA coordinates.

    `a` is the plus-strand genomic start of the transcription unit and
    `pre_len` its length; `strand` is the unit's orientation. Returns None
    when the span doesn't lie inside the unit — a partial overlap included,
    because half of a feature mapped into a transcript is a coordinate lie,
    not a partial answer.

    On a minus-strand unit the span is FLIPPED as well as shifted: the
    transcript reads 3'→5' along the plus strand, so its first base is the
    unit's LAST plus-strand base. Getting this backwards is the classic
    origin-and-strand fault, so it is derived once here and tested against a
    sequence whose content identifies its own coordinates."""
    if total <= 0 or pre_len <= 0:
        return None
    off_s = (gs - a) % total
    off_e = (ge - a) % total
    # A span ending exactly at the unit's end wraps its offset to 0.
    if off_e == 0 and off_s != 0:
        off_e = pre_len
    if off_s > pre_len or off_e > pre_len or off_s >= off_e:
        return None
    if strand >= 0:
        return (off_s, off_e)
    return (pre_len - off_e, pre_len - off_s)


def _tx_feature_parts(feat: dict) -> "list[tuple[int, int]]":
    """A feature's location parts as ``[(start, end), ...]``.

    A spliced CDS/mRNA arrives as a compound location, and the gaps BETWEEN
    its parts are the introns — the standard GenBank way to say so. A feature
    with no `parts` is a single span."""
    raw = feat.get("parts")
    out: "list[tuple[int, int]]" = []
    if isinstance(raw, (list, tuple)):
        for p in raw:
            try:
                if isinstance(p, dict):
                    out.append((int(p["start"]), int(p["end"])))
                else:
                    out.append((int(p[0]), int(p[1])))
            except (KeyError, TypeError, ValueError, IndexError):
                continue
    if not out:
        try:
            out = [(int(feat["start"]), int(feat["end"]))]
        except (KeyError, TypeError, ValueError):
            return []
    return out


def _tx_usable_features(features, warnings: "list[str]") -> "list[dict]":
    """Feature dicts whose coordinates are actually usable.

    A record can carry a feature with a non-integer or missing position — a
    hand-edited GenBank, a `<1..>500` fuzzy location, or a caller-built list.
    Reading one used to raise `ValueError` out of the middle of the transcript
    walk, which reached the endpoint as an opaque 500. It is not worth failing
    the whole unit over one bad feature, but it IS worth saying so: a silently
    shorter feature list could be the reason the promoter went missing."""
    out: "list[dict]" = []
    dropped = 0
    for f in (features or []):
        if not isinstance(f, dict):
            dropped += 1
            continue
        try:
            int(f["start"]); int(f["end"])
            int(f.get("strand", 1) or 1)
        except (KeyError, TypeError, ValueError):
            dropped += 1
            continue
        out.append(f)
    if dropped:
        warnings.append(
            f"{dropped} feature(s) skipped — their coordinates aren't whole "
            "numbers (a fuzzy or hand-edited location)")
    return out


def _feat_label_or(feat, fallback: str = "unlabelled") -> str:
    """A feature dict's display label, or `fallback`. Tolerates a non-dict so
    it can be used inside a message built before the union is narrowed."""
    if not isinstance(feat, dict):
        return fallback
    return str(feat.get("label") or "") or fallback


def _tx_pick_feature(features, want, types, *, strand=None):
    """Resolve one feature by explicit label / index, else by type.

    `want` may be an int index into `features`, a label string, or None.
    Returns ``(feature, error_message)``. An explicit request that doesn't
    resolve is an ERROR, never a silent fall-back to a type guess — picking a
    different promoter than the one asked for and reporting a transcript for
    it is exactly the kind of confident wrong answer this module avoids."""
    if isinstance(want, bool):
        return None, "expected a feature index or label"
    if isinstance(want, int):
        if not (0 <= want < len(features)):
            return None, f"feature index {want} out of range"
        return features[want], None
    if isinstance(want, str) and want.strip():
        needle = want.strip().lower()
        hits = [f for f in features
                if str(f.get("label") or "").lower() == needle]
        if not hits:
            hits = [f for f in features
                    if needle in str(f.get("label") or "").lower()]
        if not hits:
            return None, f"no feature labelled {want!r}"
        if len(hits) > 1:
            return None, (f"{len(hits)} features match {want!r} — "
                          "pass an index instead")
        return hits[0], None
    cands = [f for f in features
             if str(f.get("type") or "").lower() in types
             and (strand is None or int(f.get("strand", 1) or 1) == strand)]
    if not cands:
        return None, None                  # caller decides if that's fatal
    return cands, None


def _tx_translate(mrna: str, start: int) -> "tuple[str, int]":
    """Translate `mrna` from `start` to the first in-frame stop.

    Returns ``(protein, stop_pos)`` where `stop_pos` is the offset of the
    stop codon, or ``-1`` when the frame runs off the 3' end without one
    (which is itself a finding — an unterminated ORF in a mature message
    means the annotation and the sequence disagree)."""
    aa: "list[str]" = []
    i = start
    n = len(mrna)
    while i + 3 <= n:
        codon = mrna[i:i + 3]
        if codon in _STOP_CODONS:
            return "".join(aa), i
        aa.append(_CODON_TABLE.get(codon, "X"))
        i += 3
    return "".join(aa), -1


def _predict_transcript(seq: str, features: "list[dict]", *,
                        circular: bool = True,
                        promoter=None, cds=None, terminator=None,
                        tx_start=None, tx_end=None,
                        check_splice: bool = True,
                        min_uorf_aa: int = 1) -> dict:
    """Reconstruct the mature mRNA of one annotated transcription unit and
    score translation initiation on it.

    `features` are dicts in `list-features` shape — ``{start, end, strand,
    type, label}`` — optionally carrying ``parts`` (a compound location's
    exons). `promoter` / `cds` / `terminator` each take a feature index, a
    label, or None to pick by type. `tx_start` / `tx_end` override the unit's
    genomic bounds directly (plus-strand coordinates) when the record has no
    promoter or terminator annotated.

    Returns a dict with ``ok``, the ``unit`` that was read, the ``pre_mrna``
    and ``mature_mrna``, the ``exons`` / ``introns`` that separate them, the
    ``five_utr`` / ``cds`` / ``three_utr`` blocks, the ``kozak`` context, the
    ``uorfs`` found upstream of the start codon, the ``removed_by_splicing``
    list, an advisory ``splice`` scan and any ``warnings``.

    **`removed_by_splicing` is the point of the whole function.** It lists
    every upstream ATG that IS in the unspliced sequence and is NOT in the
    mature message. Those are precisely the entries a uATG screen run against
    the plasmid reports as hazards and which the spliceosome deletes before
    a ribosome ever sees them.

    The transcript is taken to start at the END of the promoter (a promoter is
    bound, not transcribed — the same reasoning `_assemble_expression_cassette`
    uses for its splice scan) and to run to the end of the terminator. Introns
    are read from annotation only; see the module note above.
    """
    seq_u = (seq or "").upper()
    total = len(seq_u)
    warnings: "list[str]" = []
    if total == 0:
        return {"ok": False, "error": "empty sequence", "warnings": warnings}
    feats = _tx_usable_features(features, warnings)

    # ── 1. The coding feature fixes the unit's strand. ────────────────────
    picked, err = _tx_pick_feature(feats, cds, {"cds"})
    if err:
        return {"ok": False, "error": f"cds: {err}", "warnings": warnings}
    if isinstance(picked, list):
        if not picked:
            return {"ok": False,
                    "error": "no CDS feature on this record — pass 'cds' with "
                             "a label or index",
                    "warnings": warnings}
        # Longest by the wrap-aware length, never by `end - start`.
        picked = max(picked, key=lambda f: _feat_len(
            int(f["start"]) % total, int(f["end"]) % total, total))
        if len([f for f in feats
                if str(f.get("type") or "").lower() == "cds"]) > 1:
            warnings.append(
                "more than one CDS on this record; took the longest "
                f"({_feat_label_or(picked)}) — pass "
                "'cds' to choose")
    if not isinstance(picked, dict):        # unreachable; narrows the union
        return {"ok": False, "error": "could not resolve a CDS feature",
                "warnings": warnings}
    cds_feat: dict = picked
    strand = 1 if int(cds_feat.get("strand", 1) or 1) >= 0 else -1

    # ── 2. Bound the unit. ────────────────────────────────────────────────
    def _anchor(want, types, which, default_edge):
        f, e = _tx_pick_feature(feats, want, types, strand=strand)
        if e:
            return None, f"{which}: {e}"
        if isinstance(f, list):
            if not f:
                return None, None
            # Nearest one on the correct side of the CDS, measured AROUND the
            # circle in transcript direction — not by raw coordinate order,
            # which is meaningless across the origin.
            def _dist(g):
                cs = int(cds_feat["start"]) % total
                ce = int(cds_feat["end"]) % total
                if which == "promoter":
                    return ((cs - int(g["end"]) % total) % total if strand >= 0
                            else (int(g["start"]) % total - ce) % total)
                return ((int(g["end"]) % total - ce) % total if strand >= 0
                        else (cs - int(g["start"]) % total) % total)
            f = min(f, key=_dist)
        return f, None

    prom, e = _anchor(promoter, {"promoter"}, "promoter", None)
    if e:
        return {"ok": False, "error": e, "warnings": warnings}
    term, e = _anchor(terminator, {"terminator"}, "terminator", None)
    if e:
        return {"ok": False, "error": e, "warnings": warnings}

    def _coord(v):
        try:
            return int(v) % total
        except (TypeError, ValueError):
            return None

    if tx_start is not None or tx_end is not None:
        a_opt, b_opt = _coord(tx_start), _coord(tx_end)
        if a_opt is None or b_opt is None:
            return {"ok": False,
                    "error": "tx_start and tx_end must BOTH be integers when "
                             "either is given",
                    "warnings": warnings}
    elif prom is None or term is None:
        return {"ok": False, "error": _tx_missing_anchor(prom, term),
                "warnings": warnings}
    elif strand >= 0:
        a_opt, b_opt = _coord(prom.get("end")), _coord(term.get("end"))
    else:
        a_opt, b_opt = _coord(term.get("start")), _coord(prom.get("start"))
    if a_opt is None or b_opt is None:
        return {"ok": False,
                "error": "the promoter/terminator coordinates are not "
                         "integers",
                "warnings": warnings}
    a, b = a_opt, b_opt
    pre_len = _feat_len(a, b, total)
    if pre_len <= 0:
        return {"ok": False,
                "error": "the transcription unit is empty — check the "
                         "promoter/terminator order and strand",
                "warnings": warnings}
    if not circular and b < a:
        return {"ok": False,
                "error": "the transcription unit crosses the origin, but this "
                         "record is LINEAR",
                "warnings": warnings}
    pre_plus = _slice_circular(seq_u, a, b)
    pre = pre_plus if strand >= 0 else _rc(pre_plus)

    def _to_tx(gs, ge):
        cs, ce = _coord(gs), _coord(ge)
        if cs is None or ce is None:
            return None
        return _tx_map_span(cs, ce, a, total, strand, pre_len)

    cds_parts_tx = []
    for ps, pe in _tx_feature_parts(cds_feat):
        m = _to_tx(ps, pe)
        if m is not None:
            cds_parts_tx.append(m)
    if not cds_parts_tx:
        return {"ok": False,
                "error": "the CDS does not lie inside the transcription unit "
                         "— check the promoter/terminator choice and strand",
                "warnings": warnings}
    cds_parts_tx.sort()

    # ── 3. Introns: annotated only. ───────────────────────────────────────
    introns: "list[list]" = []
    for f in feats:
        if str(f.get("type") or "").lower() != "intron":
            continue
        m = _to_tx(f.get("start"), f.get("end"))
        if m is None:
            continue
        introns.append([m[0], m[1], str(f.get("label") or "")])
    # Gaps BETWEEN a spliced location's parts, computed in transcript
    # coordinates so an origin-crossing CDS needs no special case.
    for (s1, e1), (s2, _e2) in zip(cds_parts_tx, cds_parts_tx[1:]):
        if s2 > e1:
            introns.append([e1, s2, "CDS location gap"])
    introns.sort()
    merged: "list[list]" = []
    for iv in introns:
        if merged and iv[0] < merged[-1][1]:
            if iv[1] > merged[-1][1]:
                merged[-1][1] = iv[1]
                merged[-1][2] = (merged[-1][2] + " + " + iv[2]).strip(" +")
            continue
        merged.append(list(iv))
    exons: "list[tuple[int, int]]" = []
    cur = 0
    for s, e, _lbl in merged:
        if s > cur:
            exons.append((cur, s))
        cur = max(cur, e)
    if cur < pre_len:
        exons.append((cur, pre_len))
    mature = "".join(pre[s:e] for s, e in exons)

    def _mature_pos(local):
        """pre-mRNA offset → mature offset, or None when it's intronic."""
        out = 0
        for s, e in exons:
            if local < s:
                return None
            if local < e:
                return out + (local - s)
            out += e - s
        return out if local >= pre_len else None

    def _mature_end(local):
        """Half-open END mapping — an end sitting on an exon boundary belongs
        to the exon it closes, not to the intron that follows it."""
        if local <= 0:
            return 0
        p = _mature_pos(local - 1)
        return None if p is None else p + 1

    cds_tx_start = cds_parts_tx[0][0]
    cds_tx_end = max(e for _s, e in cds_parts_tx)
    cds_start_m = _mature_pos(cds_tx_start)
    cds_end_m = _mature_end(cds_tx_end)
    if cds_start_m is None or cds_end_m is None:
        return {"ok": False,
                "error": "the CDS start or end falls INSIDE an annotated "
                         "intron — the annotation is inconsistent",
                "warnings": warnings}
    # ── 4. Upstream ATGs: which survive splicing, and which don't. ───────
    uorfs = []
    for p in range(0, max(0, cds_start_m)):
        if mature[p:p + 3] != "ATG":
            continue
        protein, stop = _tx_translate(mature, p)
        if stop == -1:
            kind = "no_in_frame_stop"
        elif stop + 3 <= cds_start_m:
            kind = "upstream"                    # a true uORF: starts and
        elif (cds_start_m - p) % 3 == 0:         # stops before the CDS
            kind = "in_frame_extension"          # N-terminal extension
        else:
            kind = "out_of_frame_overlap"        # the worst kind
        if len(protein) < min_uorf_aa:
            continue
        uorfs.append({
            "position":        p,
            "relative_to_cds": p - cds_start_m,
            "kind":            kind,
            "length_aa":       len(protein),
            "protein":         protein[:200],
            "stop_position":   stop,
            "kozak":           _kozak_context(mature, p),
        })
        if len(uorfs) >= _TX_MAX_UORFS:
            warnings.append(
                f"more than {_TX_MAX_UORFS} upstream ATGs — list truncated")
            break

    # THE point of this function: upstream ATGs the plasmid sequence shows and
    # the mature message does not. A uATG screen run on the unspliced sequence
    # reports every one of these as a hazard; none of them reaches a ribosome.
    removed = []
    for p in range(0, min(cds_tx_start, max(0, len(pre) - 2))):
        if pre[p:p + 3] != "ATG":
            continue
        if _mature_pos(p) is not None:
            continue
        holder = next((lbl for s, e, lbl in merged if s <= p < e), "")
        removed.append({
            "pre_mrna_position": p,
            "genomic":           _tx_to_genomic(p, a, total, strand, pre_len),
            "intron":            holder or "annotated intron",
        })
        if len(removed) >= _TX_MAX_UORFS:
            break

    # ── 5. Advisory cryptic-splice scan over the PRE-mRNA. ────────────────
    # The spliceosome acts on the unspliced message, so that is what gets
    # scanned; sites coinciding with an ANNOTATED boundary are marked as such
    # rather than reported as cryptic. Never used to alter the transcript.
    splice: dict = {"skipped": True, "reason": "not requested"}
    if check_splice:
        annotated_edges = {s for s, _e, _l in merged} | {e for _s, e, _l in merged}
        try:
            hits = _splice_scan(pre, kind="both")
            for h in hits:
                h["annotated"] = (h.get("position") in annotated_edges
                                  or h.get("position", -9) + 2 in annotated_edges)
            cryptic = [h for h in hits if not h["annotated"]]
            splice = {
                "skipped":  False,
                "model":    _splice_model_summary().get("source", "?"),
                "scanned":  "pre_mrna",
                "n_sites":  len(hits),
                "n_cryptic": len(cryptic),
                "sites":    cryptic[:_TX_MAX_SPLICE_HITS],
            }
            if cryptic:
                warnings.append(
                    f"{len(cryptic)} cryptic splice site(s) in the pre-mRNA "
                    "score as strongly as real sites — they compete with the "
                    "annotated ones")
        except (RuntimeError, ValueError) as exc:
            # An untrained model must read as NOT CHECKED, never as clean.
            splice = {"skipped": True, "reason": str(exc)}
            warnings.append("cryptic-splice scan SKIPPED — " + str(exc))

    protein, stop = _tx_translate(mature, cds_start_m)
    if stop == -1:
        warnings.append(
            "the CDS has no in-frame stop inside the mature message — the "
            "annotation and the spliced sequence disagree")
    elif stop + 3 != cds_end_m:
        warnings.append(
            f"the first in-frame stop is at mature position {stop}, but the "
            f"CDS annotation ends at {cds_end_m} — splicing may have shifted "
            "the frame, or the annotation is stale")

    return {
        "ok": True,
        "unit": {
            "strand":      strand,
            "genomic":     {"start": a, "end": b, "length": pre_len},
            "wraps_origin": bool(b < a or (b == a and pre_len == total)),
            "promoter":    (prom or {}).get("label") or None,
            "cds":         cds_feat.get("label") or None,
            "terminator":  (term or {}).get("label") or None,
        },
        "pre_mrna":    pre,
        "mature_mrna": mature,
        "spliced":     bool(merged),
        "exons":   [{"start": s, "end": e, "length": e - s} for s, e in exons],
        "introns": [{"start": s, "end": e, "length": e - s, "label": lbl,
                     "donor": pre[s:s + 2], "acceptor": pre[max(0, e - 2):e],
                     "canonical": pre[s:s + 2] == "GT" and pre[e - 2:e] == "AG"}
                    for s, e, lbl in merged],
        "five_utr":  {"seq": mature[:cds_start_m], "length": cds_start_m},
        "cds":       {"start": cds_start_m, "end": cds_end_m,
                      "length": cds_end_m - cds_start_m,
                      "protein": protein},
        "three_utr": {"seq": mature[cds_end_m:],
                      "length": max(0, len(mature) - cds_end_m)},
        "kozak":     _kozak_context(mature, cds_start_m),
        "uorfs":     uorfs,
        "removed_by_splicing": removed,
        "splice":    splice,
        "warnings":  warnings,
    }


def _tx_to_genomic(local: int, a: int, total: int,
                   strand: int, pre_len: int) -> int:
    """Transcript offset → plus-strand genomic coordinate. The inverse of
    `_tx_map_span` for a single position, so a hazard reported in transcript
    space can be pointed at on the map."""
    if total <= 0:
        return 0
    off = local if strand >= 0 else (pre_len - 1 - local)
    return (a + off) % total


def _tx_missing_anchor(prom, term) -> str:
    missing = [n for n, v in (("promoter", prom), ("terminator", term))
               if v is None]
    return (f"no {' or '.join(missing)} feature on the same strand as the CDS "
            "— annotate one, or pass tx_start + tx_end to bound the unit "
            "directly")
