"""splicecraft_splice — plant splice-site scoring (layer L0).

Answers one question: **would a plant's spliceosome plausibly use this GT as a
donor, or this AG as an acceptor?** That matters when designing an intronless
transgene, because a codon-optimized CDS can accidentally spell a strong
splice site, and the plant will then splice a chunk out of the message. The
protein disappears and the DNA looks perfect — one of the nastier silent
failure modes in plant expression.

**Why not a regex.** `GT` and `AG` each occur roughly every 16 bp by chance, so
a pattern match flags essentially the entire CDS and forces meaningless edits.
The dinucleotide is necessary but nowhere near sufficient: real sites carry
information across a ~9 nt (donor) or ~23 nt (acceptor) window. Scoring that
window is what separates "a GT" from "a usable donor".

**The model.** Position weight matrices, scored as a log-odds ratio of two
EMPIRICAL models:

    score = sum_i log2( P_real(base_i at position i) / P_decoy(base_i at i) )

The null model is not uniform base composition. It is real GT/AG dinucleotides
that sit *inside* verified introns and are therefore demonstrably NOT used as
boundaries — exactly the population the scrubber is trying to discriminate
against. A composition-based background would make every GT look informative
and inflate every score.

Positive: a score above 0 means the window looks more like a real site than
like an unused one.

**Clade-matched training matters.** The shipped matrices are trained on
verified introns from a reference plant genome in the same clade as the
intended host. Intron recognition is NOT uniform across plants — some clades
lean hard on intron UA-richness, others tolerate GC-richer introns — so a
matrix trained on a mismatched clade carries the wrong bias into the design
and will mis-rank sites. Retrain with `_splice_train` if a better-matched set
arrives.

**Thresholds are FALSE-POSITIVE-anchored,** and that is the single most
important thing to know before reading a result. The cutoff is set where the
DECOY rate hits `fpr_target` (0.01) — about one chance hit per 100 unused
GT/AG dinucleotides — and the sensitivity that buys is then MEASURED and
recorded, not chosen:

    donor     threshold 3.766   sensitivity 0.4593   decoy FPR 0.00998
    acceptor  threshold 4.439   sensitivity 0.3325   decoy FPR 0.0092

So **a clean scan is not a clean sequence.** At these cutoffs the model flags
fewer than half of genuine donors and a third of genuine acceptors; a real
plant intron's own 3' site scoring 3.2 against the 4.439 bar is the expected
behaviour of a 0.33-sensitivity test, not a defect. `0 sites found` means "none
that clears a 1%-false-positive bar", which is a much weaker statement than
"nothing a spliceosome would use". Read `_splice_scan_basis()` alongside any
count — it ships the sensitivity WITH the result so a caller cannot quote the
count without the detection power — and pass an explicit lower
`donor_threshold` / `acceptor_threshold` (or a `margin`, for the codon
scrubber) when a miss is costlier than a needless synonymous swap. It usually
is: a missed site silently deletes protein, a false alarm costs one codon.

An earlier version of this paragraph claimed the cutoff was the 5th percentile
of real-site scores, i.e. 95% sensitivity. It never was — that is ~2x the real
donor figure and ~3x the acceptor one, and since this is the sentence a reader
uses to decide what "0 sites" means, the overstatement made clean-looking
scans out of sequences carrying usable sites (field report 2026-11-08).

Layer L0: imports only logging. Pure functions over strings.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from math import log2

from splicecraft_logging import _log

# Window geometry, matching the MaxEntScan convention so scores are comparable
# to published work.
#   donor    9 nt = exon[-3:] + intron[:6],  GT at window offsets 3..4
#   acceptor 23 nt = intron[-20:] + exon[:3], AG at window offsets 18..19
_SPLICE_DONOR_LEN = 9
_SPLICE_DONOR_GT = 3
_SPLICE_ACCEPTOR_LEN = 23
_SPLICE_ACCEPTOR_AG = 18

_SPLICE_BASES = "ACGT"

# Longest scoring window of either kind. A change of a few bases can only move
# the score of a window that OVERLAPS those bases, so this is the half-width a
# caller needs to rescan locally instead of rescanning the whole sequence
# (`_codon_fix_splice_sites` leans on exactly that).
_SPLICE_MAX_WINDOW = max(_SPLICE_DONOR_LEN, _SPLICE_ACCEPTOR_LEN)

# Pseudocount as a FRACTION of the sample, not a flat count. This matters more
# than it looks. The real and decoy sets differ in size by ~18x, so a flat
# Laplace +1 smooths them by different amounts: at an invariant position (the
# GT / AG anchor, where every sequence in BOTH classes carries the same base) a
# base that never occurs would score
#     log2( (1/N_real) / (1/N_decoy) ) = log2(N_decoy/N_real) ~ +4.2 bits
# — i.e. an IMPOSSIBLE base would outscore the correct one. Smoothing each
# class by the same proportion makes the log-odds at an invariant position
# exactly 0, which is the truthful answer: a position that is constant in both
# populations cannot discriminate between them.
_SPLICE_PSEUDO_FRAC = 0.001

# ── Trained matrices ───────────────────────────────────────────────────────
# Generated by `_splice_train` (see outputs/tools/train_splice_model.py in the
# thaumatin project for the reproducible pipeline). Each is a list of
# per-position {base: log2-odds} dicts. Provenance + calibration live in
# _SPLICE_MODEL_META so a stale matrix is visible rather than silent.
_SPLICE_DONOR_PWM: "list[dict[str, float]]" = []
_SPLICE_ACCEPTOR_PWM: "list[dict[str, float]]" = []
_SPLICE_MODEL_META: dict = {}


def _splice_count_matrix(windows, width: int) -> "list[dict[str, float]]":
    """Pseudocounted base frequencies per position."""
    cols = [{b: 0.0 for b in _SPLICE_BASES} for _ in range(width)]
    used = 0
    for w in windows:
        if len(w) != width:
            continue
        w = w.upper()
        if set(w) - set(_SPLICE_BASES):
            continue
        used += 1
        for i, base in enumerate(w):
            cols[i][base] += 1.0
    if not used:
        raise ValueError(f"no usable training windows of width {width}")
    # Proportional smoothing — see _SPLICE_PSEUDO_FRAC.
    pseudo = max(_SPLICE_PSEUDO_FRAC * used, 1e-9) / len(_SPLICE_BASES)
    for col in cols:
        for base in col:
            col[base] += pseudo
    for col in cols:
        total = sum(col.values())
        for base in col:
            col[base] /= total
    return cols


def _splice_train(real, decoy, width: int) -> "list[dict[str, float]]":
    """Build a per-position log2-odds matrix from real vs decoy windows.

    Both models are estimated from data; the decoy set defines the null. See
    the module docstring for why that beats a composition background."""
    p_real = _splice_count_matrix(real, width)
    p_decoy = _splice_count_matrix(decoy, width)
    return [{b: log2(p_real[i][b] / p_decoy[i][b]) for b in _SPLICE_BASES}
            for i in range(width)]


def _splice_score_window(window: str, pwm: "list[dict[str, float]]") -> "float | None":
    """Log2-odds score for one window. None when it cannot be scored (wrong
    length, or an ambiguity code) — never a silent 0.0, which would read as
    'perfectly neutral site' rather than 'no answer'."""
    if not pwm or len(window) != len(pwm):
        return None
    window = window.upper()
    total = 0.0
    for i, base in enumerate(window):
        col = pwm[i]
        if base not in col:
            return None
        total += col[base]
    return total


def _splice_score_donor(seq: str, gt_pos: int) -> "float | None":
    """Score the donor whose GT begins at `gt_pos` (0-based, forward strand).

    Returns None when there is no `GT` there. The anchor is deliberately NOT
    part of the score: the decoy population is GT-anchored too, so those two
    positions carry zero discriminative information and contribute 0 bits. The
    model's job is to rank GTs against each other, not to decide whether a
    dinucleotide is a GT — so that check lives here, as a precondition. Without
    it a `CAG|ATAAGT` window would score identically to `CAG|GTAAGT` and a
    caller scoring arbitrary offsets would get a confident number for a site
    that cannot splice."""
    if seq[gt_pos:gt_pos + 2].upper() != "GT":
        return None
    start = gt_pos - _SPLICE_DONOR_GT
    if start < 0 or start + _SPLICE_DONOR_LEN > len(seq):
        return None
    return _splice_score_window(seq[start:start + _SPLICE_DONOR_LEN],
                                _SPLICE_DONOR_PWM)


def _splice_score_acceptor(seq: str, ag_pos: int) -> "float | None":
    """Score the acceptor whose AG begins at `ag_pos` (0-based, forward).

    Returns None when there is no `AG` there — same precondition rationale as
    `_splice_score_donor`."""
    if seq[ag_pos:ag_pos + 2].upper() != "AG":
        return None
    start = ag_pos - _SPLICE_ACCEPTOR_AG
    if start < 0 or start + _SPLICE_ACCEPTOR_LEN > len(seq):
        return None
    return _splice_score_window(seq[start:start + _SPLICE_ACCEPTOR_LEN],
                                _SPLICE_ACCEPTOR_PWM)


def _splice_threshold(kind: str) -> "float | None":
    meta = _SPLICE_MODEL_META.get(kind) or {}
    return meta.get("threshold")


def _splice_scan_basis(*, donor_threshold: "float | None" = None,
                       acceptor_threshold: "float | None" = None) -> dict:
    """The detection power behind a scan, for shipping WITH its result.

    ``{donor: {threshold, calibrated_threshold, sensitivity, decoy_fpr},
    acceptor: {...}, fpr_target, note}``.

    This exists because a bare count is unreadable. The cutoffs are anchored
    at a 1% decoy false-positive rate, which buys ~0.46 donor / ~0.33 acceptor
    sensitivity, so `n_sites: 0` means "nothing clearing a 1%-FPR bar" and NOT
    "nothing a spliceosome would use" — the distinction a reader cannot make
    from the number alone. Every consumer that reports a site count reports
    this beside it.

    `sensitivity` is reported as `None` when the caller overrode the cutoff:
    the recorded figure was measured at the calibrated threshold and does not
    describe a different one. Saying nothing is correct; re-using the number
    would be a measurement claim nobody made."""
    out: dict = {}
    for kind, override in (("donor", donor_threshold),
                           ("acceptor", acceptor_threshold)):
        meta = _SPLICE_MODEL_META.get(kind) or {}
        cal = meta.get("threshold")
        used = cal if override is None else float(override)
        same = override is None or (cal is not None and float(override) == cal)
        out[kind] = {
            "threshold": used,
            "calibrated_threshold": cal,
            "sensitivity": meta.get("sensitivity") if same else None,
            "decoy_fpr": meta.get("decoy_fpr") if same else None,
        }
    out["fpr_target"] = _SPLICE_MODEL_META.get("fpr_target")
    out["note"] = (
        "thresholds are anchored at a ~1% decoy false-positive rate, NOT at a "
        "sensitivity target: a count of 0 means 'no site clears a 1%-FPR bar', "
        "not 'no usable site'. Fewer than half of real donors and a third of "
        "real acceptors clear it.")
    return out


def _splice_scan(seq: str, *, kind: str = "both",
                 donor_threshold: "float | None" = None,
                 acceptor_threshold: "float | None" = None,
                 both_strands: bool = False) -> "list[dict]":
    """Every cryptic splice site in `seq` scoring at or above threshold.

    Returns ``[{kind, position, score, threshold, above_calibrated, window,
    strand}]`` sorted by descending score. `position` is the 0-based offset of
    the GT (donor) or AG (acceptor) dinucleotide on the FORWARD strand.

    `above_calibrated` is the hit's standing against the MODEL's own cutoff,
    independent of the `threshold` this scan ran at, and is computed from the
    unrounded score (`None` when the model records no cutoff). It is what a
    relaxed scan — the pairing stage — uses to say which ends are genuinely
    sub-threshold; re-deriving it from the rounded `score` gets borderline
    sites wrong.

    `both_strands` is off by default and should usually stay off: an mRNA is
    single-stranded and only the sense strand is spliced, so a reverse-strand
    "site" in a CDS is not a splicing hazard. It is available for scanning
    regions whose transcribed orientation is unknown or bidirectional.

    Raises RuntimeError when the model has not been trained — a scan that
    silently returned [] would read as 'clean' when it means 'not checked'."""
    if not _SPLICE_DONOR_PWM or not _SPLICE_ACCEPTOR_PWM:
        raise RuntimeError(
            "splice model is not trained — _SPLICE_DONOR_PWM/_SPLICE_ACCEPTOR_PWM "
            "are empty. An untrained scan cannot report 'no sites found'.")
    seq = (seq or "").upper()
    if kind not in ("both", "donor", "acceptor"):
        raise ValueError("kind must be 'both', 'donor' or 'acceptor'")
    if donor_threshold is None:
        donor_threshold = _splice_threshold("donor")
    if acceptor_threshold is None:
        acceptor_threshold = _splice_threshold("acceptor")
    # `_splice_threshold` returns None when the metadata carries no calibrated
    # cutoff. That is a MODEL fault, not a scan input fault, and it would
    # otherwise reach the `score >= threshold` comparisons below and die on
    # `float >= None` — an opaque TypeError from deep inside the loop. Same
    # class as the untrained-PWM guard above: refuse loudly and say why.
    if donor_threshold is None or acceptor_threshold is None:
        raise RuntimeError(
            "splice model has no calibrated threshold — the model metadata is "
            "missing 'threshold' for donor and/or acceptor. Retrain the model, "
            "or pass explicit donor_threshold / acceptor_threshold.")

    strands = [("+", seq)]
    if both_strands:
        comp = str.maketrans("ACGT", "TGCA")
        strands.append(("-", seq.translate(comp)[::-1]))

    # Each hit records whether it clears the CALIBRATED cutoff, computed here
    # from the UNROUNDED score. A caller scanning at a relaxed threshold (the
    # pairing stage) must not re-derive this from the reported `score`: that
    # one is rounded to 3 dp, so a site at 3.76599 reads as exactly 3.766 and
    # compares as >= the 3.766 bar it actually fails. Two halves answering the
    # same question differently is the bug class this whole module keeps
    # tripping over, so the scan answers it ONCE.
    _cal = {"donor": _splice_threshold("donor"),
            "acceptor": _splice_threshold("acceptor")}

    def _above(kind: str, score: float) -> "bool | None":
        cal = _cal.get(kind)
        return None if cal is None else score >= cal

    hits: "list[dict]" = []
    for strand, s in strands:
        n = len(s)
        for i in range(n - 1):
            pair = s[i:i + 2]
            if kind in ("both", "donor") and pair == "GT":
                score = _splice_score_donor(s, i)
                if score is not None and score >= donor_threshold:
                    pos = i if strand == "+" else n - i - 2
                    hits.append({"kind": "donor", "position": pos,
                                 "score": round(score, 3),
                                 "threshold": donor_threshold,
                                 "above_calibrated": _above("donor", score),
                                 "window": s[i - _SPLICE_DONOR_GT:
                                             i - _SPLICE_DONOR_GT + _SPLICE_DONOR_LEN],
                                 "strand": strand})
            if kind in ("both", "acceptor") and pair == "AG":
                score = _splice_score_acceptor(s, i)
                if score is not None and score >= acceptor_threshold:
                    pos = i if strand == "+" else n - i - 2
                    hits.append({"kind": "acceptor", "position": pos,
                                 "score": round(score, 3),
                                 "threshold": acceptor_threshold,
                                 "above_calibrated": _above("acceptor", score),
                                 "window": s[i - _SPLICE_ACCEPTOR_AG:
                                             i - _SPLICE_ACCEPTOR_AG + _SPLICE_ACCEPTOR_LEN],
                                 "strand": strand})
    hits.sort(key=lambda h: -h["score"])
    return hits


def _splice_window_span(kind: str, position: int) -> "tuple[int, int]":
    """``(start, end)`` of the scoring window behind a `_splice_scan` hit.

    `position` is the hit's GT / AG offset; the window reaches BACK over the
    exon for a donor and FORWARD into the exon for an acceptor, so neither
    span is centred on the dinucleotide. Exposed because a caller editing the
    sequence (synonymous-codon scrubbing) has to know which bases carry a
    site's score — deriving it from the window STRING works only until a
    sequence repeats itself."""
    if kind == "donor":
        s = position - _SPLICE_DONOR_GT
        return (s, s + _SPLICE_DONOR_LEN)
    if kind == "acceptor":
        s = position - _SPLICE_ACCEPTOR_AG
        return (s, s + _SPLICE_ACCEPTOR_LEN)
    raise ValueError("kind must be 'donor' or 'acceptor'")


# ── Pair-stage calibration ────────────────────────────────────────────────
# The pairing stage deliberately scans LOWER than the reporting stage. A pair
# carries two independent pieces of evidence, so requiring each end to clear a
# cutoff sized for a LONE site discards the shape that actually bites: a
# mediocre donor plus a strong acceptor. Measured on a real case (field report
# 2026-11-08) a donor at 3.596 against the 3.766 bar paired with an acceptor at
# 5.135 — combined 8.73, comfortably over the bar for a pair — and the default
# scan could not see it at all.
#
# Two numbers, neither invented:
#
#   * `_SPLICE_PAIR_FLOOR = 0.0` is the log-odds break-even the module
#     docstring already defines ("more like a real site than like an unused
#     one"). MEASURED: the gated pair count SATURATES at or above this floor —
#     identical output at floors +0.5, 0.0 and -1.0 across GC 0.35/0.45/0.55/
#     0.65 — so nothing below 0.0 can reach the combined gate and a lower
#     floor only costs scan time. It is a compute bound, not a judgement.
#   * the combined gate is `donor_threshold + acceptor_threshold`, i.e. a pair
#     must carry as much total evidence as two individually-passing sites.
#     Derived from the two calibrated cutoffs, so it moves with a retrained
#     model and adds no third magic number.
#
# Flagging burden, MEASURED against the old both-ends-pass rule: on 15.5 kb of
# real GenBank-extracted coding sequence, 0.06 -> 0.19 pairs per kb; on random
# DNA at GC 0.45, 0.6 -> 5.2 per kb (random sequence carries far more chance
# pyrimidine tracts than real CDS, which is what the acceptor matrix keys on).
# The rule is a strict SUPERSET of the old one — every pair the old default
# reported still qualifies, verified over 12 x 6 kb in test_splice.py — so this
# can only add, never silence. The asymmetry is deliberate and matches
# `splicecraft_regulatory`: a MISSED pair silently deletes protein from the
# product, a false alarm costs one synonymous codon.
_SPLICE_PAIR_FLOOR = 0.0


def _splice_pair_risk(seq: str, *, pair_floor: "float | None" = None,
                      combined_threshold: "float | None" = None,
                      annotated_edges: "set | None" = None,
                      **kw) -> "list[dict]":
    """Donor/acceptor PAIRS that could excise a real cryptic intron.

    A lone strong donor is much less dangerous than a donor followed, at a
    plausible intron distance, by a strong acceptor — that pair is what
    actually removes sequence from the message. Plant introns are short
    (training-set median ~150 nt), so the window checked is 60–3000 nt.

    Because a pair carries two independent pieces of evidence, the SITES are
    found at a relaxed `pair_floor` (default `_SPLICE_PAIR_FLOOR`) and the PAIR
    is then gated on `combined_threshold` (default: the sum of the two
    calibrated single-site cutoffs) rather than on two separate cutoffs. See
    the block comment above for why, and for the measured burden. An explicit
    `donor_threshold` / `acceptor_threshold` in `kw` still wins — it replaces
    the floor for that kind — so a caller can tighten or loosen either end.

    Returns ``[{donor, acceptor, intron_len, combined_score, in_frame_loss,
    both_above_threshold, sub_threshold_ends, annotated_ends}]`` sorted by
    combined score. `in_frame_loss` flags an excision that is a multiple of 3:
    those are the truly insidious ones, since they delete residues without
    frameshifting and the product may still be detected by an antibody while
    being biologically wrong.

    `both_above_threshold` / `sub_threshold_ends` say WHICH ends cleared their
    own calibrated cutoff, so the stricter old answer is still recoverable
    from the richer one and a relaxed pair is never mistaken for a confirmed
    one.

    `annotated_edges` is a set of positions that are REAL annotated intron
    boundaries; `annotated_ends` then lists which end of each pair sits on one
    (`[]`, `["donor"]`, `["acceptor"]` or both). A pair with one annotated end
    is a different animal from a fully cryptic one — it is usually the real
    donor reaching past its own partner because the model missed the real
    acceptor (a 0.33-sensitivity test misses two in three) — and without the
    key nothing distinguishes them. The caller owns the annotation because
    this layer knows nothing about features.

    Splicing is STRAND-SPECIFIC, so a donor is only ever paired with an
    acceptor on its OWN strand. `_splice_scan` reports both strands in
    FORWARD coordinates, so on the minus strand the intron runs right to
    left and the acceptor sits at the LOWER coordinate — its span is
    measured the other way round. Measuring every combination in forward
    coordinates (the pre-fix behaviour, reachable only via
    ``both_strands=True``) invented plus-donor/minus-acceptor "introns"
    that cannot splice, and discarded every genuine minus-strand pair,
    whose forward-coordinate span is negative."""
    floor = _SPLICE_PAIR_FLOOR if pair_floor is None else float(pair_floor)
    cal_d, cal_a = _splice_threshold("donor"), _splice_threshold("acceptor")
    # `_splice_scan` stamps `above_calibrated` per hit, so the pair stage never
    # compares the rounded `score` against a cutoff itself.
    # An explicit per-kind threshold wins over the floor: a caller asking for a
    # specific cutoff must get it, not a relaxed one.
    kw.setdefault("donor_threshold", floor)
    kw.setdefault("acceptor_threshold", floor)
    if combined_threshold is None:
        if cal_d is None or cal_a is None:
            # Same class as `_splice_scan`'s guard: an uncalibrated model must
            # refuse, never quietly gate on a number it does not have.
            raise RuntimeError(
                "splice model has no calibrated threshold — cannot derive the "
                "pair gate. Retrain the model, or pass combined_threshold.")
        gate = cal_d + cal_a
    else:
        gate = float(combined_threshold)
        # A NaN gate compares False against every score, so EVERY pair would
        # be dropped and the answer would read as "no excisable intron" —
        # silence standing in for a refusal, which is the failure mode this
        # module keeps being audited for.
        if gate != gate:
            raise ValueError("combined_threshold must be a number, not NaN")
    if floor != floor:
        raise ValueError("pair_floor must be a number, not NaN")

    edges = annotated_edges or set()

    def _ann(h) -> bool:
        """Does this hit sit on an annotated boundary? Both the dinucleotide
        start and its far base count, matching `_predict_transcript`."""
        if h.get("annotated"):
            return True
        pos = h.get("position")
        return pos in edges or (pos is not None and pos + 2 in edges)

    hits = _splice_scan(seq, **kw)
    pairs = []
    for strand in ("+", "-"):
        donors = sorted([h for h in hits
                         if h["kind"] == "donor" and h.get("strand", "+") == strand],
                        key=lambda h: h["position"])
        acceptors = sorted([h for h in hits
                            if h["kind"] == "acceptor"
                            and h.get("strand", "+") == strand],
                           key=lambda h: h["position"])
        # Only acceptors inside the 60..3000 nt intron window can pair with a
        # given donor, and `acceptors` is sorted by position, so the candidates
        # are a CONTIGUOUS slice — found by bisect instead of walked. Without
        # this the loop is O(donors x acceptors) over the whole molecule: a
        # 1 Mbp AT-rich sequence at the relaxed pair floor put 8,352 donors
        # against 6,346 acceptors, i.e. 53 million combinations and 6.6 s, all
        # but ~0.5% of them rejected by the span test. Relaxing the floor (the
        # 2026-11-08 field report) multiplied the hit counts ~10x each, which
        # is what turned a tolerable quadratic into a reachable one — the
        # endpoint takes a sequence straight from the caller. Same output,
        # asserted against the exhaustive walk in test_splice.py.
        apos = [a["position"] for a in acceptors]
        for d in donors:
            dp = d["position"]
            if strand == "+":
                lo, hi = dp + 60 - 2, dp + 3000 - 2
            else:
                lo, hi = dp + 2 - 3000, dp + 2 - 60
            for a in acceptors[bisect_left(apos, lo):bisect_right(apos, hi)]:
                span = (a["position"] + 2 - dp if strand == "+"
                        else dp + 2 - a["position"])
                if not (60 <= span <= 3000):
                    continue
                combined = d["score"] + a["score"]
                if combined < gate:
                    continue
                sub = [k for k, h in (("donor", d), ("acceptor", a))
                       if h.get("above_calibrated") is False]
                pairs.append({
                    "donor": d, "acceptor": a, "intron_len": span,
                    "combined_score": round(combined, 3),
                    "in_frame_loss": span % 3 == 0,
                    "both_above_threshold": not sub,
                    "sub_threshold_ends": sub,
                    "annotated_ends": [k for k, h in (("donor", d),
                                                      ("acceptor", a))
                                       if _ann(h)],
                })
    pairs.sort(key=lambda p: -p["combined_score"])
    return pairs


def _splice_model_summary() -> dict:
    """What model is loaded, and how it was calibrated."""
    return {
        "trained": bool(_SPLICE_DONOR_PWM and _SPLICE_ACCEPTOR_PWM),
        "donor_width": len(_SPLICE_DONOR_PWM),
        "acceptor_width": len(_SPLICE_ACCEPTOR_PWM),
        **_SPLICE_MODEL_META,
    }


def _splice_load_model(donor_pwm, acceptor_pwm, meta) -> None:
    """Install a trained model (used by the vendored block below and by tests
    that want a deliberately tiny model)."""
    global _SPLICE_DONOR_PWM, _SPLICE_ACCEPTOR_PWM, _SPLICE_MODEL_META
    _SPLICE_DONOR_PWM = donor_pwm
    _SPLICE_ACCEPTOR_PWM = acceptor_pwm
    _SPLICE_MODEL_META = meta
    _log.info("splice model loaded: donor %d nt, acceptor %d nt, source %s",
              len(donor_pwm), len(acceptor_pwm), meta.get("source", "?"))


def _splice_consensus(pwm: "list[dict[str, float]]", *,
                      min_bits: float = 0.15) -> str:
    """Consensus encoded by the matrix — uppercase where the position carries
    real discriminative signal, `.` where it does not.

    **The GT / AG anchor positions read as `.`, and that is correct.** The
    decoy set is itself GT- and AG-anchored (unused GT/AG dinucleotides inside
    real introns), so at those positions both classes are constant, the
    log-odds is ~0 for every base, and `max()` would return an arbitrary
    letter. Printing that letter as if it were a consensus is actively
    misleading — it looks like the model failed to learn `GT` when in fact
    `GT` carries no information *relative to the population being
    discriminated against*.

    What the model must get right is the INFORMATIVE positions. For a correctly
    trained donor that is `CAG` on the exon side and `AAGT` at intron
    +3..+6 (the `CAG|GTAAGT` consensus); for the acceptor, a long pyrimidine
    tract and `G` at exon +1. Those are what `test_splice.py` asserts."""
    out = []
    for col in pwm:
        best = max(col, key=lambda b: col[b])
        spread = col[best] - min(col.values())
        out.append(best if spread >= min_bits else ".")
    return "".join(out)


def _splice_position_bits(pwm: "list[dict[str, float]]") -> "list[float]":
    """Per-position discriminative range (max - min log-odds), in bits.

    Near-zero identifies the anchor positions; the profile as a whole is what
    to inspect when a retrained model behaves oddly."""
    return [round(max(c.values()) - min(c.values()), 3) for c in pwm]


# ── Load the vendored model at import. Optional by design: the scoring code is
# useful with a caller-supplied model (tests, a retrained species-specific
# matrix), and a missing vendored block must not make the module unimportable.
# `_splice_scan` raises rather than returning [] when nothing is loaded, so an
# absent model can never be mistaken for a clean result.
try:                                    # pragma: no cover - import-time wiring
    from splicecraft_splice_model import (
        ACCEPTOR_PWM as _VENDORED_ACCEPTOR,
        DONOR_PWM as _VENDORED_DONOR,
        SPLICE_MODEL_META as _VENDORED_META,
    )
except ImportError:                     # pragma: no cover
    _log.warning("no vendored splice model found — _splice_scan will raise "
                 "until a model is loaded via _splice_load_model")
else:                                   # pragma: no cover
    _splice_load_model(_VENDORED_DONOR, _VENDORED_ACCEPTOR, _VENDORED_META)
