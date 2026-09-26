"""splicecraft_regulatory — bacterial regulatory-element scanners (L1).

Two scanners that answer "where does transcription START and where does it
STOP", on the plasmid's own forward-strand coordinates:

  * ``_scan_promoters``    — sigma-70 (-35 / spacer / -10) promoter search.
  * ``_scan_terminators``  — intrinsic (Rho-independent) terminator search:
                             a stem-loop followed by a U-tract.

Both are wrap-aware on a circular molecule and both report FORWARD-strand
coordinates for hits on either strand (sacred invariant #2 — a reverse hit at
``p`` is stored at ``p``, not mirrored). They feed ``_map_transcription``
(seqanalysis L3), which walks the circle from every promoter to every CDS and
needs the two halves in one coordinate frame to do it.

WHAT THESE DELIBERATELY DO NOT DO
---------------------------------
Neither returns a predicted transcription RATE, and neither claims to.

``_score_promoter`` returns a composite of NAMED, separately-reported
components — how well each hexamer matches the sigma-70 consensus, weighted by
which positions are actually conserved, plus a spacer-geometry term — and the
result carries ``model: "sigma70 consensus composite (relative, not
calibrated)"``. It is built for RANKING candidate promoters against each other
on one molecule, which is the question that gets asked ("does this insert carry
a cryptic promoter, and is it stronger than the others?"). It is NOT a
promoter-strength prediction, there is no trained model behind it, and a number
that looked like one would be worse than no number at all. This is the same
discipline ``splicecraft_crispr._score_guide`` applies to on-target efficiency
and ``splicecraft_splice`` applies to GT/AG regexes — see ``[INV-198]``.

Terminators DO carry a real thermodynamic number: ``dg`` comes from the
Turner-2004 nearest-neighbour fold in biology L0, not from a heuristic. The
``efficiency`` field beside it is an explicit COARSE TIER
(``strong``/``moderate``/``weak``) over the geometry, never a percentage.

THE BUG THIS MODULE EXISTS TO NOT HAVE
--------------------------------------
A hand-rolled terminator scan reported 22 terminators where there was one,
because every overlapping stem/loop REGISTER of the same hairpin scored as its
own hit. ``_scan_terminators`` is organised around that: candidates are
generated per U-tract, the best register per hairpin wins, and overlapping
survivors are merged. One hairpin is one record — asserted by
``test_regulatory.py::TestNonRedundancy``.
"""
from __future__ import annotations

from splicecraft_biology import _rc, _rna_fold
from splicecraft_logging import _log

# ── sigma-70 promoter geometry ────────────────────────────────────────
# The consensus hexamers, and per-position weights reflecting which bases are
# actually conserved across characterised E. coli promoters. The -35 TTG and
# the -10 invariant T (-12) plus its 3'-most T (-7) carry most of the
# recognition; the middle positions of both boxes drift widely in real
# promoters. Weights are RELATIVE importance for ranking, not frequencies —
# see the module docstring on why no frequency matrix is asserted here.
_P35_CONSENSUS = "TTGACA"
_P35_WEIGHTS = (1.0, 1.0, 1.0, 0.6, 0.4, 0.5)
_P10_CONSENSUS = "TATAAT"
_P10_WEIGHTS = (1.0, 0.8, 0.5, 0.5, 0.4, 1.0)
# Spacer between the boxes. 17 bp is optimal; 15-19 is the usable range, and
# outside it the two boxes cannot present to one sigma factor at once.
_SPACER_MIN = 15
_SPACER_MAX = 19
_SPACER_OPT = 17
_SPACER_PENALTY_PER_BP = 0.06
# Extended -10: a TG immediately 5' of the -10 hexamer can substitute for a
# weak -35 entirely, so a scan that only scored the two hexamers would miss a
# whole class of real promoters.
_EXT10_MOTIF = "TG"
_EXT10_BONUS = 0.10
# Transcription start site, as an offset from the -10 hexamer's 5' END (its
# first base, 0-based). The hexamer sits at -12..-7 and +1 is the start, so from
# its first base the start is 12 bases along — equivalently the ~7 nt downstream
# of its 3' end that the literature quotes. Written from the 5' end because that
# is the index the scan has in hand (`j`).
#
# This was 7 — the 3'-end spacing applied to the 5'-end index — which put every
# reported TSS 5 bp early, while the scan WINDOW below was sized correctly for
# 13. The two halves disagreed for as long as both existed (audit 2026-09-22).
_TSS_OFFSET_FROM_MINUS10 = 12
# Default reporting floor, CALIBRATED rather than guessed — an arbitrary floor
# is how a scan reports 245 "promoters" in 2.7 kb and buries the one that
# matters. Two measurements set it:
#
#   * Background, over 200 kb of random equal-base DNA, one strand:
#       score >= 0.55 -> 1 hit per     21 bp   (unusable)
#       score >= 0.70 -> 1 hit per    327 bp
#       score >= 0.75 -> 1 hit per    966 bp   <- default
#       score >= 0.80 -> 1 hit per  3,125 bp
#       score >= 0.90 -> 1 hit per 50,000 bp
#   * Real promoters, scored from this repo's own verified preset catalogue
#     (`splicecraft_presets`, sequences extracted from public GenBank records):
#       trc 1.000, J23119 1.000, tac 0.940, pTetL 0.804, EM7 0.800,
#       AmpR 0.774, CmR 0.687; median across 23 presets 0.774.
#
# 0.75 sits just under the median real promoter and accepts ~1 chance hit per
# kb per strand. Callers wanting only unambiguous promoters should pass
# `min_score=0.90`; callers hunting a weak cryptic promoter should lower it and
# read the numbers above as the price.
_PROMOTER_MIN_SCORE = 0.75
# Named tiers over the same calibration. Reported ALONGSIDE the score so a
# caller has a coarse answer that does not invite over-reading of 3 decimals.
_PROMOTER_STRONG = 0.90
_PROMOTER_MODERATE = 0.80

# ── Intrinsic terminator geometry ─────────────────────────────────────
_TERM_MIN_STEM = 5          # bp of duplex per side
_TERM_MAX_STEM = 20
_TERM_MIN_LOOP = 3          # nt — below this the strand cannot turn
_TERM_MAX_LOOP = 10
_TERM_MIN_U_TRACT = 4       # T's on the coding strand, 3' of the hairpin
_TERM_U_TRACT_WINDOW = 9    # how far past the stem the U-tract is measured
_TERM_U_TRACT_MAX_GAP = 1   # non-T bases tolerated inside the tract
# How far upstream of a U-tract a hairpin may close. A real intrinsic
# terminator's stem sits immediately 5' of the tract; a wide window would let
# an unrelated hairpin elsewhere claim the tract.
_TERM_UPSTREAM_WINDOW = 60
# Scoring tiers, CALIBRATED against this repo's own verified terminator presets
# (sequences extracted from public GenBank records) and a random background.
# Coarse ON PURPOSE (see module docstring) — a percentage here would be an
# invented number.
#
#   Real terminators, best hit per preset:
#       rrnB T1      dg -31.3  stem 20  U 6
#       T7 Tphi      dg -19.7  stem 13  U 6
#       BBa_B1008    dg -17.9  stem  8  U 5
#       T7Te         dg -14.2  stem  8  U 5
#       lambda t0    dg -12.1  stem  8  U 9
#       U6 (pol III) dg -10.0  stem  6  U 5
#     (CYC1/NOS/ADH1 are eukaryotic and are NOT intrinsic terminators; they
#      score at or above -6.4 and are correctly not called.)
#
#   Random equal-base DNA, one strand, per 100 kb:
#       dg <=  -5  -> 1 hit per  1,250 bp
#       dg <=  -8  -> 1 hit per  4,166 bp   <- reporting floor
#       dg <= -12  -> 1 hit per 25,000 bp   <- "strong" tier
#
# The two thresholds are picked for ASYMMETRIC cost. In `map-transcription` a
# spurious terminator reports a gene as insulated when it is not — the exact
# false-negative that cost this project months — so `strong`, the tier that
# calls a read-through BLOCKED, is set where random sequence essentially never
# reaches: every real bacterial intrinsic terminator above clears -12 with
# margin, and chance produces one per 25 kb.
_TERM_STRONG_DG = -12.0     # kcal/mol, hairpin fold
_TERM_MODERATE_DG = -8.0
_TERM_STRONG_U = 5
# Reporting floor. A candidate whose stem does not fold to a usefully negative
# free energy is not a terminator — bases can be written as complementary and
# still not form under any condition the cell offers. `None` (unfoldable, e.g.
# over the fold length cap) is KEPT, because unknown is not the same as
# disproved.
_TERM_MAX_DG = -8.0

_COMPLEMENT_PAIRS = frozenset({
    ("A", "T"), ("T", "A"), ("G", "C"), ("C", "G"),
    # G:U wobble — real in the RNA hairpin, and excluding it under-calls
    # GC-rich stems that genuinely form.
    ("G", "T"), ("T", "G"),
})


def _clean(seq: str) -> str:
    """Upper-case A/C/G/T view of `seq`; every other symbol becomes ``N``.

    Ambiguity codes are NOT expanded. A promoter or terminator call resting on
    a base the record does not actually know is a guess presented as a finding,
    and an ``N`` simply fails to match or pair, which is the honest outcome.
    """
    out = []
    for ch in (seq or "").upper():
        out.append(ch if ch in "ACGT" else "N")
    return "".join(out)


def _hexamer_score(window: str, consensus: str, weights) -> float:
    """Weighted fraction of `consensus` matched by `window`, in ``[0, 1]``."""
    if len(window) != len(consensus):
        return 0.0
    hit = sum(w for ch, cons, w in zip(window, consensus, weights)
              if ch == cons)
    total = sum(weights)
    return (hit / total) if total else 0.0


def _scan_window(seq: str, circular: bool, extra: int) -> str:
    """The string to scan so a circular molecule's origin-spanning hits are
    found: the sequence plus its own first `extra` bases.

    Hits are reported ``% len(seq)``, so an element that straddles bp 0 lands
    at its real forward-strand start rather than being lost — the same
    wrap-scan shape sacred invariant #6 mandates for restriction sites.

    Suffix-only, which is correct for a scanner that reads FORWARD from its
    anchor (the promoter scanner anchors on the -35 box and everything else
    lies downstream). A scanner that reads BACKWARD needs `_scan_window_padded`
    instead — see the note there.
    """
    if not circular or not seq or extra <= 0:
        return seq
    extra = min(extra, len(seq))
    return seq + seq[:extra]


def _scan_window_padded(seq: str, circular: bool, pre: int,
                        post: int) -> "tuple[str, int]":
    """``(scan, offset)`` — `seq` padded on BOTH sides with its own bases.

    Returns the padded string and the index within it that corresponds to
    ``seq[0]``, so a hit at scan index ``x`` is at ``(x - offset) % len(seq)``.

    Needed by any scan that searches UPSTREAM of its anchor. The terminator
    scanner anchors on the U-tract and then looks back for the hairpin, so on a
    circular molecule a terminator whose stem sits just before bp 0 and whose
    tract sits just after it has its 5' half at the far END of the string. With
    a suffix-only window the backward search clamps at index 0 and either finds
    a truncated stem or misses the terminator entirely — which is exactly what
    happened, and it was invisible until the molecule was rotated.
    """
    if not circular or not seq:
        return seq, 0
    n = len(seq)
    pre = min(max(0, pre), n)
    post = min(max(0, post), n)
    return (seq[n - pre:] + seq + seq[:post]), pre


# ──────────────────────────────────────────────────────────────────────
#  Promoters
# ──────────────────────────────────────────────────────────────────────

def _score_promoter(minus35: str, spacer_len: int, minus10: str,
                    ext10: str = "") -> "tuple[float, dict]":
    """Composite score in ``[0, 1]`` plus its named components.

    RELATIVE, not calibrated — see the module docstring. The components are
    returned alongside so a caller can rank on one of them instead, or show its
    working, rather than trusting one opaque number.
    """
    s35 = _hexamer_score(minus35, _P35_CONSENSUS, _P35_WEIGHTS)
    s10 = _hexamer_score(minus10, _P10_CONSENSUS, _P10_WEIGHTS)
    spacer_pen = abs(int(spacer_len) - _SPACER_OPT) * _SPACER_PENALTY_PER_BP
    ext = _EXT10_BONUS if ext10 == _EXT10_MOTIF else 0.0
    # The -10 carries more weight than the -35: sigma-70 makes its
    # sequence-specific, strand-opening contact there, and an extended -10 can
    # drive a promoter with essentially no -35 at all.
    raw = (0.40 * s35) + (0.60 * s10) + ext - spacer_pen
    score = max(0.0, min(1.0, raw))
    return score, {
        "minus35_match": round(s35, 3),
        "minus10_match": round(s10, 3),
        "spacer_penalty": round(spacer_pen, 3),
        "extended_minus10": bool(ext),
    }


def _scan_promoters_one_strand(seq: str, circular: bool,
                               min_score: float) -> "list[dict]":
    """sigma-70 hits on the GIVEN strand, in that string's coordinates."""
    n = len(seq)
    if n < (6 + _SPACER_MIN + 6):
        return []
    # -35 hexamer + spacer + everything from the -10's first base through the
    # TSS inclusive. Derived from the offset so the two cannot drift apart.
    span = 6 + _SPACER_MAX + _TSS_OFFSET_FROM_MINUS10 + 1
    scan = _scan_window(seq, circular, span)
    limit = n if circular else len(scan)
    out: "list[dict]" = []
    for i in range(0, limit):
        w35 = scan[i:i + 6]
        if len(w35) < 6:
            break
        s35 = _hexamer_score(w35, _P35_CONSENSUS, _P35_WEIGHTS)
        # Cheap rejection: even a perfect -10 with the extended motif and no
        # spacer penalty cannot reach `min_score` from here, so the inner
        # spacer loop would be wasted work on most of the molecule.
        if (0.40 * s35) + 0.60 + _EXT10_BONUS < min_score:
            continue
        best = None
        for spacer in range(_SPACER_MIN, _SPACER_MAX + 1):
            j = i + 6 + spacer
            w10 = scan[j:j + 6]
            if len(w10) < 6:
                continue
            ext10 = scan[j - 2:j] if j >= 2 else ""
            score, comps = _score_promoter(w35, spacer, w10, ext10)
            if score < min_score:
                continue
            if best is None or score > best[0]:
                best = (score, spacer, w10, j, comps)
        if best is None:
            continue
        score, spacer, w10, j, comps = best
        # One record per -35 anchor: the spacer that scores best IS the
        # register the polymerase would use, and emitting all five would be
        # the same redundancy bug the terminator scanner exists to avoid.
        tss = j + _TSS_OFFSET_FROM_MINUS10
        # `% n` only on a CIRCLE. A LINEAR molecule has no way round: the -10
        # can fit while the start site lies past the last base, and wrapping
        # put that TSS at bp 0 — a promoter at the 3' end then "drove" the
        # genes at the 5' end (audit 2026-09-24). There it is None: the
        # promoter's own coordinates stand, no start site does.
        if circular:
            out.append({
                "start": i % n,
                "end": (j + 6) % n,
                "minus35": w35,
                "minus35_start": i % n,
                "spacer": spacer,
                "minus10": w10,
                "minus10_start": j % n,
                "tss": tss % n,
                "score": round(score, 3),
                "components": comps,
            })
        else:
            out.append({
                "start": i,
                "end": j + 6,
                "minus35": w35,
                "minus35_start": i,
                "spacer": spacer,
                "minus10": w10,
                "minus10_start": j,
                "tss": tss if tss < n else None,
                "score": round(score, 3),
                "components": comps,
            })
    return out


def _scan_promoters(seq: str, *, circular: bool = True,
                    both_strands: bool = True,
                    min_score: "float | None" = None) -> "list[dict]":
    """Every sigma-70 promoter candidate on `seq`, best-scoring first.

    Each hit: ``{start, end, strand, minus35, minus35_start, spacer, minus10,
    minus10_start, tss, score, components, model}`` — all coordinates
    FORWARD-strand, 0-based, ``% len(seq)`` so a hit spanning the origin lands
    at its real start.

    ``strand`` is ``1`` or ``-1``. For a ``-1`` hit every coordinate is still
    forward-strand, and ``tss`` is the base transcription would start from,
    which lies at a LOWER coordinate than ``start`` because the promoter fires
    leftwards. That asymmetry is deliberate and is what lets
    ``_map_transcription`` walk both directions in one frame.

    ``score`` is relative and uncalibrated — see the module docstring. Sorting
    by it to ask "which of these inserts carries the strongest cryptic
    promoter" is exactly what it is for; reading it as a transcription rate is
    not.
    """
    s = _clean(seq)
    if not s:
        return []
    floor = (_PROMOTER_MIN_SCORE if min_score is None
             else max(0.0, min(1.0, float(min_score))))
    n = len(s)
    hits: "list[dict]" = [dict(h, strand=1) for h in
                          _scan_promoters_one_strand(s, circular, floor)]
    if both_strands:
        rcs = _rc(s)
        for h in _scan_promoters_one_strand(rcs, circular, floor):
            # Map back to forward coordinates. A feature occupying
            # [a, b) of the reverse-complement string occupies
            # [n - b, n - a) of the forward strand.
            def fwd(pos_in_rc: int) -> int:
                return (n - 1 - (pos_in_rc % n)) % n
            end = fwd(h["minus35_start"]) + 1
            hits.append({
                **h,
                "strand": -1,
                # The element's forward-strand extent: its RC-string end maps
                # to the lower forward coordinate.
                "start": fwd(h["minus10_start"] + 5),
                "end": end % n if circular else end,
                "minus35_start": fwd(h["minus35_start"]),
                "minus10_start": fwd(h["minus10_start"]),
                "tss": fwd(h["tss"]) if h["tss"] is not None else None,
            })
    for h in hits:
        h["model"] = ("sigma70 consensus composite "
                      "(relative, not calibrated)")
        h["strength"] = ("strong" if h["score"] >= _PROMOTER_STRONG
                         else "moderate" if h["score"] >= _PROMOTER_MODERATE
                         else "weak")
    hits.sort(key=lambda h: (-h["score"], h["start"], h["strand"]))
    return hits


# ──────────────────────────────────────────────────────────────────────
#  Intrinsic terminators
# ──────────────────────────────────────────────────────────────────────

def _u_tract_length(seq: str, start: int) -> int:
    """Length of the T-run beginning at `start`, tolerating
    ``_TERM_U_TRACT_MAX_GAP`` interruptions inside it.

    Real U-tracts are not pure. Requiring a perfect run misses most
    characterised terminators; allowing unlimited gaps makes every AT-rich
    stretch a tract. One gap is the documented middle.
    """
    gaps = 0
    length = 0
    for k in range(start, min(len(seq), start + _TERM_U_TRACT_WINDOW)):
        if seq[k] == "T":
            length = k - start + 1
        else:
            gaps += 1
            if gaps > _TERM_U_TRACT_MAX_GAP:
                break
    return length


def _best_hairpin_before(seq: str, tract_start: int) -> "tuple | None":
    """The best stem-loop closing immediately 5' of `tract_start`.

    Returns ``(stem_len, loop_len, stem5_start, stem3_end)`` in `seq`
    coordinates, or None. "Best" is the longest stem; ties break to the
    tighter loop. Searching per U-tract — rather than enumerating every
    stem/loop register in the molecule — is what makes one hairpin produce one
    record instead of one per register.
    """
    best = None
    lo = max(0, tract_start - _TERM_UPSTREAM_WINDOW)
    # The stem's 3' side must end at (or within a couple of bases of) the
    # tract: an intrinsic terminator's hairpin sits immediately before the
    # U-tract, and a gap means the two are unrelated elements.
    for stem3_end in range(tract_start, max(lo, tract_start - 3) - 1, -1):
        for stem in range(_TERM_MAX_STEM, _TERM_MIN_STEM - 1, -1):
            for loop in range(_TERM_MIN_LOOP, _TERM_MAX_LOOP + 1):
                stem3_start = stem3_end - stem
                loop_start = stem3_start - loop
                stem5_start = loop_start - stem
                if stem5_start < lo:
                    continue
                left = seq[stem5_start:stem5_start + stem]
                right = seq[stem3_start:stem3_end]
                if len(left) != stem or len(right) != stem:
                    continue
                if "N" in left or "N" in right:
                    continue
                # `right` read 3'->5' pairs with `left` read 5'->3'.
                ok = all((a, b) in _COMPLEMENT_PAIRS
                         for a, b in zip(left, reversed(right)))
                if not ok:
                    continue
                cand = (stem, -loop, stem5_start, stem3_end)
                if best is None or cand > best:
                    best = cand
                # Longest stem wins for this (stem3_end, loop); no shorter
                # stem at the same geometry can beat it.
                break
    if best is None:
        return None
    stem, neg_loop, stem5_start, stem3_end = best
    return (stem, -neg_loop, stem5_start, stem3_end)


def _terminator_tier(dg: "float | None", u_tract: int) -> str:
    """Coarse strength tier. Named, never a percentage — the module docstring
    says why. Both halves must be present for ``strong``: a deep hairpin with
    a two-base tract does not terminate, and neither does a long tract with no
    hairpin."""
    if dg is not None and dg <= _TERM_STRONG_DG and u_tract >= _TERM_STRONG_U:
        return "strong"
    if dg is not None and dg <= _TERM_MODERATE_DG:
        return "moderate"
    return "weak"


def _scan_terminators_one_strand(seq: str, circular: bool, min_stem: int,
                                 min_u_tract: int) -> "list[dict]":
    """Intrinsic terminators on the GIVEN strand, in that string's
    coordinates. One record per hairpin."""
    n = len(seq)
    if n < (min_stem * 2 + _TERM_MIN_LOOP + min_u_tract):
        return []
    scan, off = _scan_window_padded(seq, circular, _TERM_UPSTREAM_WINDOW,
                                    _TERM_U_TRACT_WINDOW)
    # Anchor positions to sweep: every base of the molecule exactly once on a
    # circle (offset into the padded string), or the whole string when linear.
    t_range = range(off, off + n) if circular else range(0, len(scan))
    out: "list[dict]" = []
    for t in t_range:
        u_len = _u_tract_length(scan, t)
        if u_len < min_u_tract:
            continue
        # Only consider a tract START — otherwise every offset into the same
        # T-run generates its own candidate.
        if t > 0 and scan[t - 1] == "T":
            continue
        hp = _best_hairpin_before(scan, t)
        if hp is None:
            continue
        stem, loop, stem5_start, stem3_end = hp
        if stem < min_stem:
            continue
        hairpin = scan[stem5_start:stem3_end]
        dg = None
        try:
            _struct, dg = _rna_fold(hairpin)
        except ValueError:
            # Over the fold cap or ambiguous — geometry still stands, and a
            # missing dg is reported as None rather than as 0.0, which would
            # read as "folds to nothing".
            dg = None
        except Exception:
            _log.exception("terminator scan: fold failed at %d", stem5_start)
            dg = None
        if dg is not None and dg > _TERM_MAX_DG:
            # INCLUSIVE, as `_TERM_MAX_DG`'s own calibration comment documents
            # ("dg <= -8 -> reporting floor"). The `>=` this replaced dropped a
            # hairpin folding at exactly the floor (audit 2026-09-22).
            continue
        # An exclusive END wraps only on a circle: on a linear molecule a
        # hairpin ending at the last base ends at n, and `% n` reported it
        # as ending at 0 (round-2 hardening, 2026-09-25 — the promoter scan
        # already had this rule).
        out.append({
            "start": (stem5_start - off) % n,
            "end": ((t + u_len - off) % n if circular else t + u_len - off),
            "hairpin_start": (stem5_start - off) % n,
            "hairpin_end": ((stem3_end - off) % n if circular
                            else stem3_end - off),
            "stem_len": stem,
            "loop_len": loop,
            "u_tract": u_len,
            "u_tract_start": (t - off) % n,
            "dg": (round(dg, 2) if dg is not None else None),
            "efficiency": _terminator_tier(dg, u_len),
            # Overlap is judged as an ARC — a (start, length) pair on the
            # molecule — not as a padded-coordinate interval. Two registers of
            # one origin-spanning hairpin can be discovered in different padded
            # frames (one anchored in the left pad, one in the right), and a
            # plain `a0 < b1 and b0 < a1` test then reads them as disjoint: the
            # same hairpin came back as two records at exactly the rotations
            # that put it across the origin (audit 2026-09-22).
            "_span": ((stem5_start - off) % n, (t + u_len) - stem5_start),
        })
    return _merge_overlapping(out, total=n, circular=circular)


def _arcs_overlap(a_start: int, a_len: int, b_start: int, b_len: int,
                  total: int) -> bool:
    """Do two arcs of a circle of `total` bp share a base?

    Each arc is a START and a LENGTH, so an arc that runs past the origin needs
    no special case — which is the whole reason the spans are carried this way.
    Two arcs overlap when either one's start lies inside the other.
    """
    if a_len <= 0 or b_len <= 0 or total <= 0:
        return False
    if a_len >= total or b_len >= total:
        return True
    return (((b_start - a_start) % total) < a_len
            or ((a_start - b_start) % total) < b_len)


def _merge_overlapping(hits: "list[dict]", *, total: int = 0,
                       circular: bool = False) -> "list[dict]":
    """Collapse candidates whose spans overlap, keeping the best.

    THE point of this module's terminator half. Two U-tracts a base apart, or
    a hairpin that can close in two registers, are ONE terminator; reporting
    each separately is how a hand-rolled scan turned one hairpin into 22.
    Ranked by stem length, then tract length, then depth of fold.

    ``_span`` is ``(start, length)`` on the molecule. On a circle the overlap
    is judged as an ARC (`_arcs_overlap`): a hairpin that straddles the origin
    is one record whatever rotation the file happens to be stored in.
    """
    if not hits:
        return []
    ordered = sorted(hits, key=lambda h: (
        -h["stem_len"], -h["u_tract"], h["dg"] if h["dg"] is not None else 0.0))
    kept: "list[dict]" = []
    for h in ordered:
        a0, alen = h["_span"]
        clash = False
        for k in kept:
            b0, blen = k["_span"]
            if circular and total > 0:
                if _arcs_overlap(a0, alen, b0, blen, total):
                    clash = True
                    break
            elif a0 < b0 + blen and b0 < a0 + alen:
                clash = True
                break
        if not clash:
            kept.append(h)
    for k in kept:
        k.pop("_span", None)
    kept.sort(key=lambda h: h["start"])
    return kept


def _scan_terminators(seq: str, *, circular: bool = True,
                      both_strands: bool = True,
                      min_stem: "int | None" = None,
                      min_u_tract: "int | None" = None) -> "list[dict]":
    """Every intrinsic (Rho-independent) terminator candidate on `seq`.

    Each hit: ``{start, end, strand, hairpin_start, hairpin_end, stem_len,
    loop_len, u_tract, u_tract_start, dg, efficiency}`` — FORWARD-strand,
    0-based coordinates for hits on either strand.

    ONE RECORD PER HAIRPIN. Overlapping registers of the same stem-loop, and
    U-tracts that belong to the same terminator, are merged to the best
    candidate; see ``_merge_overlapping``.

    ``dg`` is a real Turner-2004 fold of the hairpin (kcal/mol, ``None`` if it
    could not be folded). ``efficiency`` is a coarse named tier over the
    geometry — ``strong`` / ``moderate`` / ``weak`` — never a percentage.
    Rho-DEPENDENT termination is not modelled at all: it has no sequence
    signature to scan for, and silently omitting it is honest in a way that
    scoring it would not be.
    """
    s = _clean(seq)
    if not s:
        return []
    stem_floor = _TERM_MIN_STEM if min_stem is None else max(2, int(min_stem))
    u_floor = (_TERM_MIN_U_TRACT if min_u_tract is None
               else max(1, int(min_u_tract)))
    n = len(s)
    hits: "list[dict]" = [
        dict(h, strand=1) for h in
        _scan_terminators_one_strand(s, circular, stem_floor, u_floor)]
    if both_strands:
        rcs = _rc(s)
        for h in _scan_terminators_one_strand(rcs, circular, stem_floor,
                                              u_floor):
            def fwd(pos_in_rc: int) -> int:
                return (n - 1 - (pos_in_rc % n)) % n
            _wrap = (lambda x: x % n) if circular else (lambda x: x)
            hits.append({
                **h,
                "strand": -1,
                "start": fwd((h["u_tract_start"] + h["u_tract"] - 1)),
                "end": _wrap(fwd(h["hairpin_start"]) + 1),
                "hairpin_start": fwd(h["hairpin_end"] - 1),
                "hairpin_end": _wrap(fwd(h["hairpin_start"]) + 1),
                "u_tract_start": fwd(h["u_tract_start"] + h["u_tract"] - 1),
            })
    hits.sort(key=lambda h: (h["start"], -h["stem_len"]))
    return hits
