"""splicecraft_kozak — plant translation-initiation context scoring (layer L0).

Answers one question: **of the ATGs a scanning ribosome meets on this message,
which ones does a plant actually initiate at?** That is the question behind
every upstream-ORF finding. A 5'UTR ATG with a strong context captures
ribosomes and starves the real CDS; the same ATG in a poor context is mostly
read through and harmless. Reporting both as "an upstream ATG" is the
difference between a hazard and a curiosity.

**Why not just the consensus rule.** The shipped rule — a purine at −3 and a G
at +4, which is Kozak's 1987 mammalian result — does transfer to plants, and
it is kept and reported here as the interpretable component. But it has two
limits that matter for the job:

  * it is **binary twice over**, so it produces three buckets. 38% of real
    plant starts and 13% of upstream ATGs both land in the top bucket, and
    within a bucket the rule cannot order anything. Asked "which of my four
    upstream ATGs matters most", it has no answer.
  * it calls 10% of real plant start codons `weak`, because a real start can
    carry neither feature and still be the start.

The trained matrix measured here separates real starts from the upstream ATGs
they compete with at **AUC 0.766** on held-out genes, against **0.712** for
the rule (ΔAUC +0.054, 95% CI [+0.035, +0.074], 2,000 gene-resamples), and it
still separates them INSIDE each of the rule's three buckets (AUC 0.62 / 0.68
/ 0.66) — which is the capability the rule structurally lacks. It beats the
rule in six of seven species held out one at a time, and ties in the seventh.

**The model.** One position weight matrix over a −12..+4 window, scored as a
log-odds ratio of two EMPIRICAL distributions, exactly as
`splicecraft_splice`::

    score = sum_i log2( P_start(base_i at i) / P_uATG(base_i at i) )

The null is not base composition. It is the real 5'UTR upstream ATGs of the
same transcripts — the population a uORF screen has to rank. A composition
null would make every A-rich leader look informative.

**Three things this deliberately does NOT do**, each with its killing number
in `KOZAK_MODEL_META["refused"]`:

  * **No window past +4.** Extending to +9 raises AUC to 0.791, and all of
    that gain is the CODING/UTR base-composition shift — a real start has
    coding sequence downstream and a 5'UTR ATG does not. It discriminates
    starts beautifully while saying nothing about initiation, so it is
    refused. (The in-CDS decoy population, coding on both sides, is scored
    separately for exactly this reason: AUC 0.718 there, so the advantage is
    not the shortcut.)
  * **No initiation RATE, and no yield number.** This ranks contexts against
    each other and reports a percentile against real plant starts. A number
    that looked like an efficiency and was not would be worse than none —
    the same discipline as `_score_guide` and `splicecraft_regulettuce`.
  * **No partial windows.** An ATG with fewer than 12 nt of 5' context is NOT
    scored. A shorter window sums fewer log-odds terms, so its score is not
    comparable to the calibrated distribution; a `None` that says why beats a
    number that silently means something else. In real plant 5'UTRs this
    costs 4.2% of upstream ATGs.

Layer L0: imports only logging. Pure functions over strings.
"""
from __future__ import annotations

from math import log2

from splicecraft_logging import _log

# Window geometry. The A of the ATG sits at offset `_KOZAK_UP`, so the three
# invariant ATG columns are `_KOZAK_UP .. _KOZAK_UP + 2` and +4 is
# `_KOZAK_UP + 3`. Those three columns MUST carry ~0 bits (both populations
# are ATG-anchored); the calibration script fails loudly if they do not,
# because a non-zero there means the geometry or the smoothing slipped.
_KOZAK_UP = 12
_KOZAK_DOWN = 4
_KOZAK_WIDTH = _KOZAK_UP + _KOZAK_DOWN

# The window `context` is REPORTED in, unchanged from before the matrix
# landed: −6..+4 is how the literature prints an initiation context, and it is
# what a human reads in a result. It is deliberately NOT the scored window.
_KOZAK_PRINT_UP = 6
_KOZAK_PRINT_DOWN = 4

_KOZAK_BASES = "ACGT"

# Proportional smoothing — see `splicecraft_splice._SPLICE_PSEUDO_FRAC` for why
# a flat Laplace +1 over two populations of different size makes an impossible
# base outscore the right one at an invariant position.
_KOZAK_PSEUDO_FRAC = 0.001

_KOZAK_PWM: "list[dict[str, float]]" = []
_KOZAK_QUANTILES: "list[float]" = []
_KOZAK_MODEL_META: dict = {}


def _kozak_count_matrix(windows, width: int) -> "list[dict[str, float]]":
    """Pseudocounted base frequencies per position."""
    cols = [{b: 0.0 for b in _KOZAK_BASES} for _ in range(width)]
    used = 0
    for w in windows:
        if len(w) != width:
            continue
        w = w.upper()
        if set(w) - set(_KOZAK_BASES):
            continue
        used += 1
        for i, base in enumerate(w):
            cols[i][base] += 1.0
    if not used:
        raise ValueError(f"no usable training windows of width {width}")
    pseudo = max(_KOZAK_PSEUDO_FRAC * used, 1e-9) / len(_KOZAK_BASES)
    for col in cols:
        for base in col:
            col[base] += pseudo
        total = sum(col.values())
        for base in col:
            col[base] /= total
    return cols


def _kozak_train(real, decoy, width: int = _KOZAK_WIDTH
                 ) -> "list[dict[str, float]]":
    """Per-position log2-odds matrix from real starts vs competing ATGs.

    Exposed so a user with a clade-matched set can retrain — the same escape
    hatch `_splice_train` provides. `scripts/calibrate_kozak.py` is the
    reproducible pipeline that produced the vendored matrix."""
    p_real = _kozak_count_matrix(real, width)
    p_decoy = _kozak_count_matrix(decoy, width)
    return [{b: log2(p_real[i][b] / p_decoy[i][b]) for b in _KOZAK_BASES}
            for i in range(width)]


def _kozak_score_window(window: str,
                        pwm: "list[dict[str, float]] | None" = None
                        ) -> "float | None":
    """Log2-odds score for one window, or None when it cannot be scored.

    None — never 0.0 — for a wrong length, an ambiguity code, or no loaded
    model: 0.0 is a real score meaning 'as likely a start as an upstream ATG',
    and a caller cannot tell it from 'no answer'."""
    pwm = _KOZAK_PWM if pwm is None else pwm
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


def _kozak_score_at(seq: str, atg_pos: int,
                    pwm: "list[dict[str, float]] | None" = None
                    ) -> "float | None":
    """Score the initiation context of the ATG at `atg_pos` in `seq`.

    Returns None when there is no ATG there, or when the −12..+4 window runs
    off either end. The ATG is a PRECONDITION, not part of the score: both
    trained populations are ATG-anchored, so those three columns carry no
    information and a caller scoring an arbitrary offset would otherwise get a
    confident number for a codon the ribosome cannot initiate at. Same
    rationale as `_splice_score_donor`."""
    if seq[atg_pos:atg_pos + 3].upper() != "ATG":
        return None
    start = atg_pos - _KOZAK_UP
    if start < 0 or atg_pos + _KOZAK_DOWN > len(seq):
        return None
    return _kozak_score_window(seq[start:start + _KOZAK_WIDTH], pwm)


def _kozak_percentile(score: "float | None") -> "float | None":
    """Where `score` falls in the distribution of REAL plant start codons.

    ``62.0`` means "a stronger context than 62% of real plant start codons".
    That is the only calibrated statement this model can make, and it is the
    one worth showing a user: it is a comparison to measured biology rather
    than to an arbitrary cutoff. Linearly interpolated between the vendored
    held-out quantiles; None when there is no score or no model.

    A model carrying FEWER THAN TWO quantiles has no distribution to place a
    score in, so it gets None rather than a number: there is nothing to
    interpolate between, and the one-quantile case used to divide by zero.
    `_kozak_load_model` accepts any quantile list (a test's deliberately tiny
    model is the point of that escape hatch), so the guard belongs here."""
    if score is None or len(_KOZAK_QUANTILES or ()) < 2:
        return None
    q = _KOZAK_QUANTILES
    step = 100.0 / (len(q) - 1)
    if score <= q[0]:
        return 0.0
    if score >= q[-1]:
        return 100.0
    for i in range(len(q) - 1):
        lo, hi = q[i], q[i + 1]
        if lo <= score < hi:
            # A flat run of quantiles means many real starts share this score;
            # report the bottom of the run rather than dividing by zero.
            frac = 0.0 if hi == lo else (score - lo) / (hi - lo)
            return round((i + frac) * step, 1)
    return 100.0


def _kozak_tier(score: "float | None"
                ) -> "tuple[str | None, float | None]":
    """(tier, measured upstream-ATG rate at that tier) for a score.

    The cuts are anchored on the real-start distribution, not chosen: `strong`
    is the median real plant start (and 15% of upstream ATGs reach it),
    `adequate` is the 5th percentile (76% of upstream ATGs reach it — which is
    itself the honest message that this is a weak-to-moderate discriminator,
    not a filter).

    The rate is None for `weak`, which is not an oversight: the measured
    numbers are "fraction at or above this cut", and `weak` is the band below
    the lowest cut, so it HAS no such number. A derived one would read as
    measured."""
    tiers = (_KOZAK_MODEL_META or {}).get("tiers") or {}
    if score is None or not tiers:
        return (None, None)
    for name in ("strong", "adequate"):
        t = tiers.get(name) or {}
        cut = t.get("cut")
        if cut is not None and score >= cut:
            return (name, t.get("utr_uatgs_at_or_above"))
    return ("weak", None)


def _kozak_context(mrna: str, atg_pos: int) -> dict:
    """Score the translation-initiation context around `atg_pos` in `mrna`.

    Returns the long-standing keys — ``{context, minus3, plus4, score,
    strength, truncated, model}`` — plus ``rule_model`` and a nested ``pwm``
    block carrying the trained matrix's answer.

    `context` is the −6..+4 window with the A of the ATG at offset 6, so it
    reads the way the literature prints it. `score` / `strength` remain **the
    positional consensus RULE** (a purine at −3, a G at +4, 0–2 → weak /
    adequate / strong), kept because it is interpretable and because it is
    what every other tool reports.

    ``pwm`` is the calibrated answer and the one to prefer::

        {score, percentile, tier, uatg_rate, window, note}

    `percentile` is the readable form: "stronger than N% of real plant start
    codons". `pwm.score` is None, with `note` saying why, when there is no ATG
    at `atg_pos`, when fewer than 12 nt of 5' context exist, or when no model
    is loaded — the two halves of this dict can disagree, and each carries its
    own basis so a reader is never left guessing which one answered.

    A context that runs off the 5' end of the message is marked `truncated` —
    a 5'UTR shorter than 6 nt is itself a finding, and padding it would invent
    bases.

    `mrna` is normalised to upper-case DNA (``U`` → ``T``) before anything is
    read off it, so the two halves cannot disagree about CASE or ALPHABET —
    only about biology. The matrix half always up-cased its window, so a
    lower-case message used to come back `weak` from the rule and `strong`
    from the matrix for no reason a reader could act on; and the parameter is
    named `mrna`, so a message spelled in RNA has to score the same as the
    same message spelled in DNA. `context` is therefore reported as DNA, as
    every other sequence this app prints is.

    `atg_pos` outside ``0 <= atg_pos < len(mrna)`` raises ``IndexError``. It
    used to read `plus4` from the FAR END of the string for a negative index
    (Python's own wrap) and report a confident `adequate` for a position the
    message does not have, while a too-large index raised from a bare slice
    deep inside. One check, one exception, no invented base."""
    seq = (mrna or "").upper().replace("U", "T")
    n = len(seq)
    if not 0 <= atg_pos < n:
        raise IndexError(
            f"atg_pos {atg_pos} is outside the message (0..{n - 1})")
    lo = atg_pos - _KOZAK_PRINT_UP
    truncated = lo < 0 or atg_pos + _KOZAK_PRINT_DOWN > n
    context = seq[max(0, lo):atg_pos + _KOZAK_PRINT_DOWN]
    minus3 = seq[atg_pos - 3] if atg_pos >= 3 else ""
    plus4 = seq[atg_pos + 3] if atg_pos + 4 <= n else ""
    score = int(minus3 in ("A", "G")) + int(plus4 == "G")

    win_start = atg_pos - _KOZAK_UP
    window = (seq[win_start:win_start + _KOZAK_WIDTH]
              if win_start >= 0 and atg_pos + _KOZAK_DOWN <= n else "")
    pwm_score = _kozak_score_at(seq, atg_pos) if window else None
    if pwm_score is None:
        if not _KOZAK_PWM:
            note = "no initiation model loaded"
        elif not window:
            note = (f"needs {_KOZAK_UP} nt of 5' context and "
                    f"{_KOZAK_DOWN} nt downstream; a partial window is not "
                    "comparable to the calibrated distribution")
        elif seq[atg_pos:atg_pos + 3] != "ATG":
            note = f"no ATG at this position ({seq[atg_pos:atg_pos + 3]!r})"
        elif len(window) != len(_KOZAK_PWM):
            # A caller-installed model of a different width. Named, because
            # the fall-through below used to blame the SEQUENCE for a
            # mismatch that is the MODEL's — the one note a reader cannot
            # debug from, since the bases are plainly ACGT.
            note = (f"loaded model is {len(_KOZAK_PWM)} nt wide, not "
                    f"{len(window)} — retrain or reload it")
        else:
            note = "window carries a base outside ACGT"
    else:
        note = ""
    tier, uatg_rate = _kozak_tier(pwm_score)
    return {
        "context":   context,
        "minus3":    minus3,
        "plus4":     plus4,
        "score":     score,
        "strength":  ("strong" if score == 2
                      else "adequate" if score == 1 else "weak"),
        "truncated": truncated,
        "model":     (f"trained plant initiation PWM "
                      f"{(_KOZAK_MODEL_META or {}).get('window', '?')} "
                      f"(held-out AUC "
                      f"{(_KOZAK_MODEL_META or {}).get('auc_vs_utr_uatg', '?')}"
                      f" vs 5'UTR upstream ATGs), reported beside the "
                      f"positional consensus rule"
                      if _KOZAK_PWM else
                      "no trained model loaded — rule only"),
        "rule_model": ("positional -3 purine / +4 G consensus "
                       "(Kozak 1987, mammalian; measured to transfer to "
                       "plants at AUC "
                       f"{(_KOZAK_MODEL_META or {}).get('rule_auc_vs_utr_uatg', '?')})"),
        "pwm": {
            "score":      None if pwm_score is None else round(pwm_score, 3),
            "percentile": _kozak_percentile(pwm_score),
            "tier":       tier,
            "uatg_rate":  uatg_rate,
            "window":     window,
            "note":       note,
        },
    }


def _kozak_model_summary() -> dict:
    """What initiation model is loaded, and how it was calibrated.

    Surfaced to the user and to the agent API for the same reason
    `_splice_model_summary` and `_regulettuce_calibration` are: a threshold
    folded into a source comment cannot be checked."""
    return {
        "trained": bool(_KOZAK_PWM),
        "width": len(_KOZAK_PWM),
        "n_quantiles": len(_KOZAK_QUANTILES),
        **_KOZAK_MODEL_META,
    }


def _kozak_load_model(pwm, quantiles, meta) -> None:
    """Install a trained model (used by the vendored block below, and by tests
    that want a deliberately tiny one)."""
    global _KOZAK_PWM, _KOZAK_QUANTILES, _KOZAK_MODEL_META
    _KOZAK_PWM = pwm
    _KOZAK_QUANTILES = quantiles
    _KOZAK_MODEL_META = meta
    _log.info("initiation model loaded: %d nt window, source %s",
              len(pwm), (meta or {}).get("source", "?"))


# ── Load the vendored model at import. Optional by design, as for the splice
# matrices: the scoring code is useful with a caller-supplied model, and a
# missing vendored block must not make the module unimportable. Nothing here
# returns a silent 0.0 in that state — `pwm.score` is None and `pwm.note` says
# the model is missing, while the consensus rule keeps answering.
try:                                    # pragma: no cover - import-time wiring
    from splicecraft_kozak_model import (
        INITIATION_PWM as _VENDORED_PWM,
        KOZAK_MODEL_META as _VENDORED_META,
        REAL_START_QUANTILES as _VENDORED_QUANTILES,
    )
except ImportError:                     # pragma: no cover
    _log.warning("no vendored initiation model found — _kozak_context will "
                 "report the consensus rule only")
else:                                   # pragma: no cover
    _kozak_load_model(_VENDORED_PWM, _VENDORED_QUANTILES, _VENDORED_META)
