"""splicecraft_regulettuce — plant cis-regulatory annotation (L1).

The **Regulettuce** tab's engine. Where `splicecraft_regulatory` answers "where
does a *bacterium* start and stop transcribing", this answers the questions a
*plant* construct raises:

  * ``_scan_cis_elements``   — which curated plant cis-regulatory elements does
                               this sequence carry, and how many would chance
                               alone put there?
  * ``_scan_tata_boxes``     — is there a TATA box where a plant core promoter
                               would put one, RELATIVE TO AN ANCHOR?
  * ``_scan_3prime_risk``    — does a stretch of this construct read like a
                               plant 3'-end / terminator region, i.e. like
                               somewhere the message would be cut and
                               polyadenylated?

All three are wrap-aware on a circular molecule and all three report
FORWARD-strand coordinates for hits on either strand (sacred invariant #2 — a
reverse hit at ``p`` is stored at ``p``, not mirrored).

HOW THE NUMBERS IN THIS MODULE WERE OBTAINED
--------------------------------------------
Every threshold, enrichment and rate here was MEASURED, with the measurement
recorded beside the constant it set. Nothing was taken from a review article
and nothing was typed from memory — the same discipline `splicecraft_presets`
applies to sequences and `splicecraft_regulatory` applies to its thresholds.
Reproduce the whole set with ``scripts/calibrate_regulettuce.py``.

The reference sets, all fetched from public NCBI records:

  * **3,089 verified plant transcription start sites.** TAIR/RefSeq-annotated
    ``mRNA`` features from 8 Mb of *Arabidopsis thaliana* chromosomes 1, 3 and
    5, keeping ONLY transcripts whose ``mRNA`` start lies strictly upstream of
    their CDS start — i.e. that carry an annotated 5'UTR, so the 5' end is a
    transcription start and not a start codon. **That filter is not optional.**
    The first attempt at this set used ``gene``/``mRNA`` features from
    single-gene GenBank submissions, and 92% of those "TSSs" sat exactly on an
    ATG: such submissions annotate the gene from the start codon. Measuring
    promoter structure against them would have profiled the wrong 400 bp
    entirely. QC that caught it is now a release test: ``ATG`` at offset 0
    occurs in 1.5% of the kept windows, against 1.6% expected by chance.
  * **600 real plant poly(A) sites.** RefSeq mRNAs of *A. thaliana* and *Oryza
    sativa*; the 3' end of such a record IS the annotated cleavage site. A
    trailing poly(A) tail, where the record carries one, is trimmed first.
  * **12,000 composition-matched null windows** drawn from the same 8 Mb of
    genome, plus 600 internal-transcript windows as a coding-sequence control.
    The null is REAL PLANT DNA, not shuffled bases: *Arabidopsis* genomic DNA
    is itself ~64% AT, and a uniform-composition background would make every
    AT-rich hexamer look informative.

WHAT THIS MODULE DELIBERATELY DOES NOT DO
-----------------------------------------
Three things a "plant regulatory annotator" is expected to do are absent, each
because it was built, measured, and failed. The numbers are kept here so the
next person does not spend the week rediscovering them.

**1. It does not predict poly(A) cleavage sites.** A position-specific log-odds
model over the 60 nt ending at 600 real cleavage sites — the exact method
`splicecraft_splice` uses successfully for splice donors — earns no positional
resolution whatsoever. Sliding it over held-out 3' regions put its argmax
within ±20 nt of the true site 27.2% of the time, against 29.1% by chance:
**below chance.** Its window AUC (0.864) was itself beaten by a model that
knows nothing but base composition (0.873). The plant poly(A) signal is
genuinely not a positional motif: AATAAA appears in only 28% of real 3' ends
and its median position is −64, not the −25 the mammalian model would predict.
``_scan_3prime_risk`` therefore reports a REGION and no site, and carries no
``cleavage_site`` field at all — a null one would still invite a feature arrow.

**2. It does not call plant promoters or score their strength.** A composite of
TATA presence, pyrimidine content, Initiator geometry and local GC reached AUC
0.646 against the composition-matched null. At the only background rate worth
reporting at (1 chance hit per 95 windows) its sensitivity was 5.0%. For
contrast, `splicecraft_regulatory`'s sigma-70 composite separates real
bacterial promoters cleanly, and `splicecraft_splice` scores AUC 0.962. Plant
Pol II promoters mostly do not carry a scannable consensus, which is precisely
why that field moved to large DNA language models. A promoter caller at AUC
0.646 dressed as a finding would be worse than none.

**3. It does not scan for the Y-patch.** The pyrimidine patch is described as a
plant-specific core-promoter element, and it does not survive a
composition-matched control: a ≥8 nt C/T run appears in 40.6% of real promoter
windows and 38.6% of null windows — **1.05×**. At ≥10 nt it reaches 1.24×,
still inside what `_scan_cis_elements` labels ``weak``. There is no Y-patch
scanner here because measuring one produced no signal to report.

WHY A BARE TATA SCAN IS NOT OFFERED EITHER
------------------------------------------
``TATAWA`` occurs once per 429 bp of real plant genomic DNA, one strand
(measured over 8.0 Mb). On a 10 kb plasmid, both strands, chance alone puts
~47 of them there. Position is the entire signal: against verified TSSs the
hits peak sharply in **−40..−31** (281 hits in that 10 bp bin against ~75 in
bins 100 bp away), and 10.7% of real promoters carry one in −40..−20 versus
3.4% of null windows — a real 3.17× enrichment that exists ONLY relative to a
known start. So ``_scan_tata_boxes`` REQUIRES anchors and refuses to run
without them, rather than returning the 47.

Layer L1: imports biology L0 (`_rc`, `_iupac_pattern`, `_iter_match_starts`),
util L0 and logging L0, plus four wrap-scan helpers from the sibling
`splicecraft_regulatory` (same layer, no cycle — regulatory imports neither
this module nor anything above L0). Pure functions over strings; no app, no
network, no data-dir writes.
"""
from __future__ import annotations

from splicecraft_biology import _iter_match_starts, _iupac_pattern, _rc
from splicecraft_logging import _log
# Same-layer reuse rather than four more copies of the wrap-scan contract.
# `_clean`'s "every non-ACGT base becomes N" rule in particular must stay ONE
# definition: an element called on a base the record does not know is a guess
# presented as a finding, and both modules refuse it the same way.
from splicecraft_regulatory import _clean, _scan_window

# ──────────────────────────────────────────────────────────────────────
#  The curated plant cis-element catalogue
# ──────────────────────────────────────────────────────────────────────
# Each entry carries the measurement that decides whether finding it means
# anything. `promoter_rate` is the fraction of the 3,089 verified plant
# promoter windows (400 bp upstream of a verified TSS) that contain the motif
# on EITHER strand; `null_rate` is the same fraction over 12,000
# composition-matched genomic windows of the same length; `enrichment` is their
# ratio. A transcription factor binds double-stranded DNA, so both strands
# count as one presence — which is also why `as-1`, `TGACG-motif` and
# `CGTCA-motif` are ONE entry here: TGACG reverse-complements to CGTCA, and
# catalogues that list them separately report one site three times.
#
# `informativeness` is derived from those two numbers, never asserted:
#   enriched      enrichment >= 1.25 — finding one is worth reporting
#   weak          1.10 <= enrichment < 1.25
#   uninformative 0.95 <= enrichment < 1.10 — present about as often by chance
#   ubiquitous    null_rate >= 0.80 — in most random DNA; presence says nothing
#   depleted      enrichment < 0.95 — real promoters AVOID it
#
# The catalogue deliberately keeps the uninformative, ubiquitous and depleted
# entries instead of dropping them. People come looking for `GT1-motif`
# because other tools report it, and the useful answer is not silence: it is
# "that motif is in 94.8% of random 400 bp windows, so your promoter carrying
# one is not evidence of anything".
_CIS_ELEMENTS: "tuple[dict, ...]" = (
    {"name": "site II / TCP",      "consensus": "TGGGCY",
     "responds_to": "TCP factors; ribosomal protein + cell-cycle genes",
     "promoter_rate": 0.276, "null_rate": 0.141},
    {"name": "SURE",               "consensus": "AATAGAAAA",
     "responds_to": "sugar",
     "promoter_rate": 0.043, "null_rate": 0.022},
    {"name": "G-box",              "consensus": "CACGTG",
     "responds_to": "light, ABA; bZIP factors",
     "promoter_rate": 0.088, "null_rate": 0.046},
    {"name": "Box 4",              "consensus": "ATTAAT",
     "responds_to": "light",
     "promoter_rate": 0.374, "null_rate": 0.228},
    {"name": "evening element",    "consensus": "AAAATATCT",
     "responds_to": "circadian clock (evening phase)",
     "aliases": ("EE",),
     "promoter_rate": 0.029, "null_rate": 0.018},
    {"name": "ABRE",               "consensus": "ACGTGG",
     "responds_to": "abscisic acid; overlaps the G-box family",
     "aliases": ("ABA-responsive element",),
     "promoter_rate": 0.114, "null_rate": 0.071},
    {"name": "O2-site",            "consensus": "GATGAYRTGG",
     "responds_to": "Opaque-2 bZIP; zein / seed storage",
     "promoter_rate": 0.009, "null_rate": 0.006},
    {"name": "telo-box",           "consensus": "AAACCCTAA",
     "responds_to": "cell-cycle / proliferating tissue",
     "promoter_rate": 0.042, "null_rate": 0.030},
    {"name": "ABRE core",          "consensus": "ACGTG",
     "responds_to": "abscisic acid (the 5 nt core of ABRE / G-box)",
     "promoter_rate": 0.323, "null_rate": 0.237},
    {"name": "GARE",               "consensus": "TAACAAA",
     "responds_to": "gibberellin",
     "aliases": ("GARE-motif",),
     "promoter_rate": 0.235, "null_rate": 0.176},
    {"name": "RY-element",         "consensus": "CATGCATG",
     "responds_to": "seed-specific expression; B3 factors",
     "promoter_rate": 0.014, "null_rate": 0.011},
    {"name": "CE3",                "consensus": "GCGTGTC",
     "responds_to": "ABA coupling element (works with an ABRE)",
     "aliases": ("ABRE-like",),
     "promoter_rate": 0.011, "null_rate": 0.009},
    {"name": "DRE / CRT",          "consensus": "RCCGAC",
     "responds_to": "drought, cold; DREB/CBF factors",
     "aliases": ("DRE", "C-repeat"),
     "promoter_rate": 0.121, "null_rate": 0.099},
    {"name": "W-box",              "consensus": "TTGACY",
     "responds_to": "WRKY factors; defence",
     "promoter_rate": 0.384, "null_rate": 0.349},
    {"name": "LTRE",               "consensus": "CCGAAA",
     "responds_to": "low temperature",
     "promoter_rate": 0.131, "null_rate": 0.120},
    {"name": "CCAAT-box",          "consensus": "CCAAT",
     "responds_to": "NF-Y / HAP complex",
     "aliases": ("HAP2/3/5",),
     "promoter_rate": 0.629, "null_rate": 0.585},
    {"name": "GCC-box",            "consensus": "GCCGCC",
     "responds_to": "ethylene; ERF factors",
     "promoter_rate": 0.045, "null_rate": 0.043},
    {"name": "W-box core",         "consensus": "TTGAC",
     "responds_to": "WRKY factors (the 5 nt core of the W-box)",
     "promoter_rate": 0.594, "null_rate": 0.567},
    {"name": "MYB core",           "consensus": "CNGTTA",
     "responds_to": "MYB factors",
     "promoter_rate": 0.393, "null_rate": 0.379},
    {"name": "as-1",               "consensus": "TGACG",
     "responds_to": "salicylic acid, auxin, methyl jasmonate; TGA factors",
     # ONE site, three catalogue names: TGACG is the reverse complement of
     # CGTCA, so a both-strands scan finds the same base pairs under every
     # name. Measured identically (0.293 / 0.284) because they ARE identical.
     "aliases": ("TGA-element", "TGACG-motif", "CGTCA-motif"),
     "promoter_rate": 0.293, "null_rate": 0.284},
    {"name": "TCT-motif",          "consensus": "TCTTAC",
     "responds_to": "light",
     "promoter_rate": 0.217, "null_rate": 0.211},
    {"name": "GT1-motif",          "consensus": "GRWAAW",
     "responds_to": "light; GT-1 factors",
     "promoter_rate": 0.976, "null_rate": 0.948},
    {"name": "ARR1-binding site",  "consensus": "NGATT",
     "responds_to": "cytokinin; type-B ARR factors",
     "promoter_rate": 0.990, "null_rate": 0.989},
    {"name": "HSE",                "consensus": "AGAANNTTC",
     "responds_to": "heat shock; HSF factors",
     "promoter_rate": 0.098, "null_rate": 0.104},
    {"name": "MYC / E-box",        "consensus": "CANNTG",
     "responds_to": "bHLH factors; cold, jasmonate",
     "promoter_rate": 0.733, "null_rate": 0.785},
    {"name": "I-box",              "consensus": "GATAAG",
     "responds_to": "light; GATA factors",
     "promoter_rate": 0.177, "null_rate": 0.191},
    {"name": "MBS",                "consensus": "CNGTTR",
     "responds_to": "drought; MYB binding site",
     "aliases": ("MYB binding site",),
     "promoter_rate": 0.571, "null_rate": 0.623},
    {"name": "CAT-box",            "consensus": "GCCACT",
     "responds_to": "meristem expression",
     "promoter_rate": 0.078, "null_rate": 0.091},
    {"name": "circadian element",  "consensus": "CAANNNNATC",
     "responds_to": "circadian clock",
     "promoter_rate": 0.264, "null_rate": 0.310},
    {"name": "SORLIP1",            "consensus": "GCCAC",
     "responds_to": "phytochrome A light signalling",
     "promoter_rate": 0.238, "null_rate": 0.281},
    {"name": "GC-box",             "consensus": "GGGCGG",
     "responds_to": "Sp1-type GC-rich element",
     "aliases": ("Sp1",),
     "promoter_rate": 0.015, "null_rate": 0.018},
    {"name": "AuxRE",              "consensus": "TGTCTC",
     "responds_to": "auxin; ARF factors",
     "aliases": ("TGA-box",),
     "promoter_rate": 0.164, "null_rate": 0.208},
    {"name": "TC-rich repeat",     "consensus": "ATTTCTTCA",
     "responds_to": "defence and stress response",
     "promoter_rate": 0.013, "null_rate": 0.019},
    {"name": "GA-motif",           "consensus": "GGAGGA",
     "responds_to": "GA-rich element",
     "promoter_rate": 0.118, "null_rate": 0.185},
)

# Boundaries for the `informativeness` classes above.
_CIS_ENRICHED_MIN = 1.25
_CIS_WEAK_MIN = 1.10
_CIS_UNINFORMATIVE_MIN = 0.95
_CIS_UBIQUITOUS_NULL = 0.80
# The reference-set sizes behind `promoter_rate` / `null_rate`, reported in
# every scan's `model` string so a caller can weigh them.
_CIS_N_PROMOTERS = 3089
_CIS_N_NULL = 12000
_CIS_WINDOW_BP = 400
# A scan of a chromosome against 35 motifs on both strands can produce
# millions of records. The cap is on HITS, applied per motif so one ubiquitous
# element cannot starve the rest, and the summary still counts everything.
_CIS_MAX_HITS_PER_ELEMENT = 2000

# ── TATA box ──────────────────────────────────────────────────────────
# Consensus and geometry, both measured against the 3,089 verified TSSs.
_TATA_CONSENSUS = "TATAWA"
# The peak 10 bp bin of TATAWA hits relative to a verified TSS was -40..-31
# (281 hits, against ~75 per 10 bp in bins 100+ bp upstream). The reported
# window is widened to the shoulder of that peak: -40..-20 holds 10.7% of real
# promoters' hits against a 3.4% null rate (3.17x), and -60..-10 holds 19.0%.
# `_TATA_PEAK` is the tight claim, `_TATA_WINDOW` what is still worth showing.
_TATA_PEAK = (-40, -31)
_TATA_WINDOW = (-60, -10)
_TATA_PEAK_ENRICHMENT = 3.17
# Measured over 8.0 Mb of A. thaliana genomic DNA, one strand. Quoted in the
# output so "I found 47 TATA boxes in my plasmid" is read correctly.
_TATA_GENOMIC_BP_PER_HIT = 429

# ── 3'-end-like region ────────────────────────────────────────────────
# Window length matching the span over which the real signal was measured; the
# plant 3' region is diffuse (see the docstring), so a shorter window scores
# noise and a longer one blurs past the element.
_3P_WINDOW = 160
_3P_STEP = 20
# score = 0.70 * AT fraction + 0.30 * best T fraction over any 20 nt
# These are the only two features that measured as carrying signal. Each was
# checked alone on the forward strand before being combined:
#     AT fraction alone      AUC 0.903 vs coding, 0.584 vs genomic
#     best T20 alone         AUC 0.834 vs coding, 0.639 vs genomic
# Features that were tried and added NOTHING, so are absent: A-richness in
# -40..-8, presence of an AATAAA-like hexamer anywhere (AUC 0.650), presence
# of one in the canonical -45..-6 window (0.580), longest T-run.
_3P_AT_WEIGHT = 0.70
_3P_T20_WEIGHT = 0.30
_3P_T20_LEN = 20
# The SHIPPED composite, scored as `_scan_3prime_risk` scores it (both strands,
# 160 nt window): n = 600 real plant 3' ends, 600 internal-transcript (coding)
# windows, 6,000 plant genomic windows. Means 0.648 / 0.526 / 0.630.
_3P_AUC_VS_CODING = 0.895
_3P_AUC_VS_GENOMIC = 0.575
# Reporting floor, from the measured operating curve of that same scorer:
#     score >= 0.60  sensitivity 0.808  coding FPR 0.1633
#     score >= 0.62  sensitivity 0.700  coding FPR 0.0817
#     score >= 0.64  sensitivity 0.573  coding FPR 0.0433   <- default
#     score >= 0.66  sensitivity 0.415  coding FPR 0.0250
#     score >= 0.68  sensitivity 0.270  coding FPR 0.0117   <- "strong" tier
#     score >= 0.70  sensitivity 0.185  coding FPR 0.0083
# 0.64 puts roughly two chance flags per kb of coding sequence scanned (43
# windows at step 20 x 0.043), which merge into one or two regions. 0.66 would
# halve that and cost a quarter of the sensitivity; the trade is set this way
# because a missed cryptic 3' end costs a construct and a false flag costs a
# glance. A caller hunting a weak one should lower it and read the FPR as the
# price.
#
# Genomic FPR stays high at every threshold (0.45 at the default) and that is
# not a defect: plant intergenic DNA is AT-rich, so this separates a 3' end
# from CODING sequence and cannot separate one from plant DNA in general.
# Hence the claim in `_scan_3prime_risk`'s docstring, and hence why a hit
# inside a backbone means much less than one inside a CDS being designed.
_3P_MIN_SCORE = 0.64
_3P_STRONG_SCORE = 0.68
_3P_SENSITIVITY_AT_DEFAULT = 0.573
_3P_CODING_FPR_AT_DEFAULT = 0.043


def _cis_informativeness(promoter_rate: float, null_rate: float) -> str:
    """Which of the five classes an element's two measured rates put it in.

    Derived, never stored — so a corrected measurement cannot disagree with a
    hand-written label."""
    if null_rate >= _CIS_UBIQUITOUS_NULL:
        return "ubiquitous"
    ratio = (promoter_rate / null_rate) if null_rate else float("inf")
    if ratio >= _CIS_ENRICHED_MIN:
        return "enriched"
    if ratio >= _CIS_WEAK_MIN:
        return "weak"
    if ratio >= _CIS_UNINFORMATIVE_MIN:
        return "uninformative"
    return "depleted"


def _cis_element_catalogue() -> "list[dict]":
    """Every curated plant cis-element, most informative first.

    Each entry: ``{name, consensus, length, responds_to, aliases,
    promoter_rate, null_rate, enrichment, informativeness, palindromic}``.

    ``enrichment`` and ``informativeness`` are the measured verdict on whether
    finding this element means anything — see the module docstring for the
    reference sets. Copies are returned, so a caller mutating a row cannot
    corrupt the catalogue."""
    out: "list[dict]" = []
    for e in _CIS_ELEMENTS:
        pr = float(e["promoter_rate"])
        nr = float(e["null_rate"])
        cons = str(e["consensus"])
        out.append({
            "name": str(e["name"]),
            "consensus": cons,
            "length": len(cons),
            "responds_to": str(e["responds_to"]),
            "aliases": list(e.get("aliases") or ()),
            "promoter_rate": pr,
            "null_rate": nr,
            "enrichment": (round(pr / nr, 2) if nr else None),
            "informativeness": _cis_informativeness(pr, nr),
            # A palindromic consensus is found once, not twice, by a
            # both-strands scan — the same rule sacred invariant #1 applies to
            # palindromic restriction enzymes.
            "palindromic": _rc(cons) == cons,
        })
    out.sort(key=lambda d: (-(d["enrichment"] or 0.0), d["name"]))
    return out


def _find_cis_element(name: str) -> "dict | None":
    """Catalogue entry matching `name` or any of its aliases, case-folded."""
    key = (name or "").strip().casefold()
    if not key:
        return None
    for row in _cis_element_catalogue():
        if row["name"].casefold() == key:
            return row
        if any(a.casefold() == key for a in row["aliases"]):
            return row
    return None


# ──────────────────────────────────────────────────────────────────────
#  Cis-element scan
# ──────────────────────────────────────────────────────────────────────

def _element_hits_one_pattern(s: str, consensus: str, circular: bool,
                              limit: int) -> "list[tuple[int, str]]":
    """``(start, matched_bases)`` for every occurrence of `consensus` in `s`.

    Overlapping occurrences are kept: `_iter_match_starts` advances one base
    past each match's start rather than past its end, so a tandem motif such
    as ``ACGTGTGTG`` against ``GTG`` yields all three. `re.finditer` would drop
    two of them — the same undercount sacred invariant #1's territory fixed
    for bordered restriction sites.
    """
    n = len(s)
    if not n:
        return []
    pat = _iupac_pattern(consensus)
    width = len(consensus)
    # Suffix padding, which is right for a scan that reads FORWARD from its
    # anchor; `% n` then lands an origin-spanning motif on its real start.
    scan = _scan_window(s, circular, width - 1)
    out: "list[tuple[int, str]]" = []
    for i in _iter_match_starts(pat, scan):
        if circular:
            if i >= n:
                break
        elif i + width > n:
            break
        out.append((i % n, scan[i:i + width]))
        if len(out) >= limit:
            break
    return out


def _resolve_names(names: "list[str] | None") -> "set[str] | None":
    """Canonical catalogue names for `names`, or None for "all of them".

    Raises ``ValueError`` on anything the catalogue does not know. Shared by
    `_scan_cis_elements` and `_cis_element_summary` so the two cannot disagree
    about a typo: the summary used to resolve an unknown name to ``None``,
    build a filter set containing it, and return an EMPTY summary while the
    scan raised — two halves of one answer behaving differently on the same
    input.
    """
    if names is None:
        return None
    out: "set[str]" = set()
    for nm in names:
        row = _find_cis_element(str(nm))
        if row is None:
            raise ValueError(f"unknown cis-element {nm!r}")
        out.add(str(row["name"]))
    return out


def _scan_cis_elements(seq: str, *, circular: bool = True,
                       both_strands: bool = True,
                       names: "list[str] | None" = None,
                       min_enrichment: "float | None" = None,
                       include_uninformative: bool = True,
                       max_hits_per_element: "int | None" = None,
                       ) -> "list[dict]":
    """Every curated plant cis-element occurrence in `seq`.

    Each hit: ``{name, consensus, start, end, strand, match, length,
    enrichment, informativeness, responds_to, model}`` — FORWARD-strand,
    0-based coordinates, ``% len(seq)`` on a circle so an origin-spanning motif
    lands on its real start.

    ``strand`` is ``1`` (motif reads on the given strand), ``-1`` (reads on the
    complement) or **``0``** — the same bases read as the motif either way,
    which is every palindromic consensus and any degenerate one whose reverse
    complement also matches there. One stretch of DNA is ONE record: emitting a
    forward and a reverse hit over the same bases is how a scan turns one site
    into two, and a transcription factor binds the duplex, not a strand.

    ``include_uninformative=False`` drops the ``uninformative`` / ``ubiquitous``
    / ``depleted`` classes. It defaults to KEEPING them, because the useful
    answer to "does my promoter have a GT1-motif" is not silence but the
    measured 94.8% null rate that comes attached to the hit. Use
    ``min_enrichment`` for a numeric floor instead.

    ``max_hits_per_element`` bounds what is DISPLAYED, PER ELEMENT PER STRAND
    — a non-palindromic motif can therefore return up to ``2 * cap`` records,
    which is the deliberate shape: the cap exists so one ubiquitous element
    cannot starve the rest of the catalogue, and halving a strand's budget to
    hit an exact total would hide hits on one strand while showing them on the
    other. `_cis_element_summary`'s arithmetic is never capped.

    Finding an element here is a statement about the SEQUENCE, not a prediction
    that the factor binds: occupancy depends on chromatin, tissue and the rest
    of the promoter, none of which a string carries. Pair every scan with
    `_cis_element_summary`, which says how many of these chance alone explains.
    """
    s = _clean(seq)
    if not s:
        return []
    cap = (_CIS_MAX_HITS_PER_ELEMENT if max_hits_per_element is None
           else max(1, int(max_hits_per_element)))
    wanted = _resolve_names(names or None)
    floor = None if min_enrichment is None else float(min_enrichment)
    n = len(s)
    hits: "list[dict]" = []
    for row in _cis_element_catalogue():
        if wanted is not None and row["name"] not in wanted:
            continue
        if floor is not None and (row["enrichment"] or 0.0) < floor:
            continue
        if not include_uninformative and row["informativeness"] not in (
                "enriched", "weak"):
            continue
        cons = row["consensus"]
        width = row["length"]
        # span -> strand, so the same bases cannot be emitted twice.
        seen: "dict[int, int]" = {}
        match_at: "dict[int, str]" = {}
        # A palindromic consensus reads as the motif on BOTH strands wherever
        # it occurs, so its hits are strand 0 by construction and the reverse
        # scan below is skipped — the same "scan forward only" rule sacred
        # invariant #1 applies to palindromic restriction enzymes, for the same
        # reason: scanning both strands would double-count every site.
        fwd_strand = 0 if row["palindromic"] else 1
        for start, matched in _element_hits_one_pattern(s, cons, circular, cap):
            seen[start] = fwd_strand
            match_at[start] = matched
        if both_strands and not row["palindromic"]:
            # Scan for the reverse complement of the CONSENSUS on the forward
            # string rather than reverse-complementing the whole sequence:
            # the hit is then already in forward coordinates and no mapping
            # step can get the off-by-one wrong (the same choice
            # `splicecraft_crispr` makes for reverse PAMs).
            for start, matched in _element_hits_one_pattern(
                    s, _rc(cons), circular, cap):
                seen[start] = 0 if start in seen else -1
                match_at.setdefault(start, matched)
        if len(seen) >= cap:
            # The cap trims what is DISPLAYED, never `_cis_element_summary`'s
            # arithmetic. Logged so a short hit list on a large molecule is
            # explainable rather than mysterious.
            _log.debug("cis-element scan: %s hit the %d-hit display cap on "
                       "%d bp", row["name"], cap, n)
        for start in sorted(seen):
            end = start + width
            hits.append({
                "name": row["name"],
                "consensus": cons,
                "length": width,
                "start": start,
                "end": (end % n) if circular else end,
                "strand": seen[start],
                "match": match_at[start],
                "enrichment": row["enrichment"],
                "informativeness": row["informativeness"],
                "responds_to": row["responds_to"],
            })
    model = (f"curated plant cis-element consensus match; enrichment measured "
             f"on {_CIS_N_PROMOTERS} verified plant promoters vs "
             f"{_CIS_N_NULL} composition-matched genomic windows "
             f"({_CIS_WINDOW_BP} bp). A match is sequence, not occupancy.")
    for h in hits:
        h["model"] = model
    hits.sort(key=lambda h: (h["start"], h["name"]))
    return hits


def _composition(s: str) -> "dict[str, float]":
    """Base frequencies of `s` over A/C/G/T, ignoring every other symbol.

    The null for "how many of these motifs would chance put here" must use the
    QUERIED sequence's own composition. A uniform 25% null makes every AT-rich
    motif look enriched in AT-rich DNA, which is most plant DNA and most
    plasmid backbones."""
    counts = {b: s.count(b) for b in "ACGT"}
    total = sum(counts.values())
    if not total:
        return {b: 0.25 for b in "ACGT"}
    return {b: counts[b] / total for b in "ACGT"}


_IUPAC_SETS = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _expected_occurrences(consensus: str, comp: "dict[str, float]",
                          positions: int) -> float:
    """Expected count of `consensus` in `positions` independent sites.

    Zeroth-order model: the per-site probability is the product over the
    consensus's positions of the summed frequency of the bases that position
    accepts. Deliberately NOT a Markov model — a higher-order null would be
    more accurate and would also stop the number being explainable in one
    line, and its job here is to put a believable order of magnitude beside
    the observed count, not to be a p-value."""
    p = 1.0
    for ch in consensus.upper():
        allowed = _IUPAC_SETS.get(ch)
        if not allowed:
            return 0.0
        p *= sum(comp.get(b, 0.0) for b in allowed)
    return p * max(0, positions)


def _cis_element_summary(seq: str, *, circular: bool = True,
                         both_strands: bool = True,
                         names: "list[str] | None" = None,
                         include_uninformative: bool = True,
                         hits: "list[dict] | None" = None,
                         ) -> "list[dict]":
    """Observed vs chance-expected count per cis-element for `seq`.

    Each row: ``{name, consensus, observed, expected, ratio, promoter_rate,
    null_rate, enrichment, informativeness, verdict}``.

    `hits` lets a caller that already ran `_scan_cis_elements` hand the list
    back instead of paying for a second scan — but it MUST be an uncapped scan
    with the same `names` / `both_strands` / `include_uninformative`, or
    `observed` understates and every ratio below it is wrong.

    THE point of this module's cis half. A 6 bp motif occurs by chance roughly
    every 2 kb per strand, so a 10 kb plasmid shows a handful of every element
    in the catalogue and a bare hit list reads as a finding when it is
    arithmetic. ``ratio`` is observed/expected on this sequence's own base
    composition; ``verdict`` is a coarse named tier over it — never a p-value,
    because the zeroth-order null behind it does not earn one.
    """
    s = _clean(seq)
    if hits is None:
        # UNCAPPED on purpose: a capped `observed` would understate the count
        # and make a ubiquitous element look enriched. The cap exists to bound
        # what gets DISPLAYED, never the arithmetic.
        hits = _scan_cis_elements(seq, circular=circular,
                                  both_strands=both_strands, names=names,
                                  include_uninformative=include_uninformative,
                                  max_hits_per_element=10 ** 9)
    observed: "dict[str, int]" = {}
    for h in hits:
        observed[h["name"]] = observed.get(h["name"], 0) + 1
    comp = _composition(s)
    n = len(s)
    wanted = _resolve_names(names or None)
    rows: "list[dict]" = []
    for row in _cis_element_catalogue():
        if wanted is not None and row["name"] not in wanted:
            continue
        if not include_uninformative and row["informativeness"] not in (
                "enriched", "weak"):
            continue
        width = row["length"]
        sites = n if circular else max(0, n - width + 1)
        exp = _expected_occurrences(row["consensus"], comp, sites)
        if both_strands and not row["palindromic"]:
            exp += _expected_occurrences(_rc(row["consensus"]), comp, sites)
        obs = observed.get(str(row["name"]), 0)
        ratio = (obs / exp) if exp > 0 else None
        if ratio is None:
            verdict = "unknown"
        elif obs == 0:
            verdict = "absent"
        elif ratio >= 2.0:
            verdict = "over-represented"
        elif ratio >= 1.25:
            verdict = "slightly over-represented"
        elif ratio >= 0.75:
            verdict = "as expected by chance"
        else:
            verdict = "under-represented"
        rows.append({
            "name": row["name"],
            "consensus": row["consensus"],
            "observed": obs,
            "expected": round(exp, 2),
            "ratio": (round(ratio, 2) if ratio is not None else None),
            "promoter_rate": row["promoter_rate"],
            "null_rate": row["null_rate"],
            "enrichment": row["enrichment"],
            "informativeness": row["informativeness"],
            "verdict": verdict,
        })
    rows.sort(key=lambda r: (-r["observed"], r["name"]))
    return rows


# ──────────────────────────────────────────────────────────────────────
#  TATA box, relative to an anchor
# ──────────────────────────────────────────────────────────────────────

def _signed_offset(pos: int, anchor: int, total: int, circular: bool) -> int:
    """`pos - anchor`, as the SHORTER way round when the molecule is circular.

    Without this a TATA box 35 bp upstream of a start site that sits near bp 0
    reports an offset of ``total - 35`` — a huge positive number — and falls
    outside every window, so the one promoter that straddles the origin is the
    one the scan misses. Same class of bug as sacred invariant #5's wrap
    midpoint."""
    d = pos - anchor
    if not circular or total <= 0:
        return d
    d %= total
    return d - total if d > total // 2 else d


def _normalise_anchors(anchors, total: int, *,
                       circular: bool = True) -> "list[dict]":
    """Accept ``int`` or ``{position, strand?, label?}`` and return dicts.

    An int means a start site on the forward strand, which is the common case
    (a CDS start read off the map). A reverse-strand gene must say so: its
    promoter lies at HIGHER forward coordinates, and defaulting it to forward
    would search the transcript body instead of the promoter.

    **`position` is only wrapped on a CIRCLE.** There, ``% total`` is the
    coordinate system: bp −40 and bp ``total − 40`` are the same base. On a
    LINEAR molecule they are not, and wrapping an out-of-range anchor invents
    a start site — an anchor typed as 99999 on a 1,206 bp linear sequence used
    to come back as bp 640 and report a TATA box 40 bp upstream of a base the
    molecule does not have. Out of range on a linear molecule raises
    ``ValueError`` instead."""
    out: "list[dict]" = []
    for raw in (anchors or ()):
        if isinstance(raw, dict):
            pos = raw.get("position")
            if pos is None:
                raise ValueError("anchor dict needs a 'position'")
            pos = int(pos)
            strand = int(raw.get("strand", 1) or 1)
            label = str(raw.get("label") or "")
        else:
            pos = int(raw)
            strand = 1
            label = ""
        if strand not in (1, -1):
            raise ValueError(f"anchor strand must be 1 or -1, got {strand!r}")
        if total > 0:
            if circular:
                pos %= total
            elif not 0 <= pos < total:
                raise ValueError(
                    f"anchor position {pos} is outside a linear {total} bp "
                    f"molecule (0..{total - 1}); a circular one would wrap it")
        out.append({"position": pos, "strand": strand, "label": label})
    return out


def _scan_tata_boxes(seq: str, *, anchors, circular: bool = True,
                     window: "tuple[int, int] | None" = None,
                     ) -> "list[dict]":
    """TATA-box candidates positioned relative to each anchor in `anchors`.

    `anchors` are the start sites to measure from — a transcription start if
    one is known, otherwise a CDS start, as ``int`` (forward strand) or
    ``{"position": int, "strand": 1|-1, "label": str}``. **Required.** It
    raises ``ValueError`` on an empty list rather than scanning the whole
    molecule, because an unanchored TATA scan is noise: ``TATAWA`` occurs once
    per 429 bp of real plant genomic DNA, so chance alone puts ~47 of them on a
    10 kb plasmid read on both strands. Position is the entire signal — see the
    module docstring for the measured −40..−31 peak.

    Each hit: ``{start, end, strand, match, anchor, anchor_label, offset,
    in_peak, tier, model}``. Coordinates are FORWARD-strand and 0-based;
    ``offset`` is in the ANCHOR's frame (negative = upstream of it), so a
    reverse-strand anchor's upstream hits lie at higher forward coordinates and
    still report a negative offset.

    ``tier`` is ``peak`` inside −40..−31, ``window`` inside −60..−10, and
    nothing else is returned. No score: the measurement supports a positional
    claim and a 3.17× enrichment, not a strength.

    An anchor outside a LINEAR molecule raises ``ValueError`` rather than
    wrapping to a base that exists — see `_normalise_anchors`. A malformed
    `window` raises ``ValueError`` too, so an agent caller gets a 400 with a
    readable message instead of a 500.
    """
    # Arguments are checked BEFORE the sequence is looked at: the
    # anchors-required rule is this scan's whole contract, and an empty or
    # all-ambiguous `seq` used to return `[]` from the line below and skip it,
    # so `_scan_tata_boxes("", anchors=[])` answered "no TATA boxes" to a call
    # that should have been refused outright.
    s = _clean(seq)
    anc = _normalise_anchors(anchors, len(s), circular=circular)
    if not anc:
        raise ValueError(
            "_scan_tata_boxes requires at least one anchor (a transcription "
            "or CDS start). An unanchored TATAWA scan returns chance hits: "
            f"one per {_TATA_GENOMIC_BP_PER_HIT} bp of real plant DNA.")
    if window is None:
        lo, hi = _TATA_WINDOW
    else:
        try:
            lo, hi = (int(window[0]), int(window[1]))
        except (TypeError, ValueError, IndexError, KeyError) as exc:
            # ValueError for every malformed window, never an IndexError from
            # a short tuple or a bare int() message about a stray character:
            # the callers above translate a ValueError into a 400 and let
            # anything else reach the 500 path.
            raise ValueError(
                f"window must be a (low, high) pair of offsets relative to "
                f"the anchor, got {window!r}") from exc
    if lo > hi:
        raise ValueError(f"window must be (low, high), got {(lo, hi)!r}")
    if not s:
        return []
    n = len(s)
    width = len(_TATA_CONSENSUS)
    fwd = _element_hits_one_pattern(s, _TATA_CONSENSUS, circular, 10 ** 9)
    rev = _element_hits_one_pattern(s, _rc(_TATA_CONSENSUS), circular, 10 ** 9)
    out: "list[dict]" = []
    for a in anc:
        pool = fwd if a["strand"] == 1 else rev
        for start, matched in pool:
            if a["strand"] == 1:
                off = _signed_offset(start, a["position"], n, circular)
            else:
                # The motif's 5' base in the GENE's frame is its highest
                # forward coordinate, so the offset is measured from there.
                off = _signed_offset(a["position"], start + width - 1, n,
                                     circular)
            if not (lo <= off <= hi):
                continue
            out.append({
                "start": start,
                "end": (start + width) % n if circular else start + width,
                "strand": a["strand"],
                "match": matched,
                "anchor": a["position"],
                "anchor_label": a["label"],
                "offset": off,
                "in_peak": _TATA_PEAK[0] <= off <= _TATA_PEAK[1],
                "tier": ("peak" if _TATA_PEAK[0] <= off <= _TATA_PEAK[1]
                         else "window"),
                "model": (
                    f"TATAWA positioned against a start site; hits peak at "
                    f"{_TATA_PEAK[0]}..{_TATA_PEAK[1]} on "
                    f"{_CIS_N_PROMOTERS} verified plant promoters "
                    f"({_TATA_PEAK_ENRICHMENT}x enriched over "
                    f"composition-matched genomic windows). Positional claim "
                    f"only -- not a promoter strength."),
            })
    out.sort(key=lambda h: (h["anchor"], -h["offset"]))
    return out


# ──────────────────────────────────────────────────────────────────────
#  3'-end-like regions
# ──────────────────────────────────────────────────────────────────────

def _best_t_fraction(w: str, k: int = _3P_T20_LEN) -> float:
    """Highest T fraction over any `k`-mer of `w`, by rolling count."""
    if len(w) < k or k <= 0:
        return 0.0
    run = w[:k].count("T")
    best = run
    for i in range(k, len(w)):
        run += (w[i] == "T") - (w[i - k] == "T")
        if run > best:
            best = run
    return best / k


def _3prime_score(w: str) -> "tuple[float, float, float]":
    """``(score, at_fraction, best_t20)`` for one window.

    The two components are the only features that measured as carrying signal
    (see `_3P_AT_WEIGHT`'s comment for the per-feature AUCs), and they are
    returned alongside the composite so the user can see that AT fraction is
    doing most of the work. That is not a weakness to hide — it is the finding:
    the plant 3' region is a compositional zone, not a motif."""
    n = len(w)
    if not n:
        return 0.0, 0.0, 0.0
    acgt = sum(w.count(b) for b in "ACGT")
    at = ((w.count("A") + w.count("T")) / acgt) if acgt else 0.0
    t20 = _best_t_fraction(w)
    return (_3P_AT_WEIGHT * at + _3P_T20_WEIGHT * t20), at, t20


def _merge_arcs(spans: "list[tuple[int, int]]", total: int,
                circular: bool) -> "list[tuple[int, int]]":
    """Union overlapping/abutting ``(start, length)`` arcs into maximal runs.

    UNION, not pick-the-best. A 400 bp AT-rich stretch flagged by 13
    overlapping windows is ONE region; keeping only the best window would
    report a 160 bp slice of it and hide the extent, which is the number the
    user needs to decide whether to recode. (`splicecraft_regulatory`'s
    `_merge_overlapping` picks a winner instead, because there one hairpin
    really is one record — different question, different merge.)

    Merged in UNWRAPPED scan coordinates, then the first and last runs are
    joined if they meet across the origin — which is why a region spanning bp 0
    comes back as one arc rather than two stubs at either end.
    """
    if not spans:
        return []
    ivs = sorted((s, s + ln) for s, ln in spans if ln > 0)
    merged: "list[list[int]]" = []
    for a, b in ivs:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    if circular and total > 0 and len(merged) > 1:
        first, last = merged[0], merged[-1]
        # The scan pads past the end, so a run that began before bp 0 shows up
        # as a run reaching past `total`.
        if last[1] % total >= first[0] and last[1] > total:
            first[0] = last[0] - total
            merged.pop()
    out: "list[tuple[int, int]]" = []
    for a, b in merged:
        length = min(b - a, total) if total else b - a
        start = (a % total) if (circular and total) else a
        out.append((start, length))
    return out


def _scan_3prime_risk(seq: str, *, circular: bool = True,
                      window: "int | None" = None, step: "int | None" = None,
                      min_score: "float | None" = None,
                      both_strands: bool = True) -> "list[dict]":
    """Regions of `seq` whose composition reads like a plant 3'-end region.

    The question this answers: *is there a stretch of my construct that a plant
    would treat as the end of the message?* A cryptic 3' end inside a transgene
    truncates the transcript, the protein disappears, and the DNA sequences
    perfectly — the same silent failure mode `splicecraft_splice` exists for.

    Each hit: ``{start, end, length, score, at_fraction, best_t_fraction,
    strand, tier, model}`` — FORWARD-strand, 0-based, wrap-aware.
    ``tier`` is ``strong`` / ``moderate``.

    **There is no cleavage-site field, and that is the finding, not an
    omission.** A position-specific model trained on 600 real plant poly(A)
    sites located them WORSE than chance (±20 nt: 27.2% against 29.1%), so this
    reports an extent and nothing finer. See the module docstring.

    **Region boundaries are quantised to `step`.** Windows start at multiples
    of `step` from bp 0, so rotating a circular molecule can move a reported
    boundary by up to `step - 1` bases. Detection and extent are
    rotation-invariant; the exact edge is not, and cannot be: there is no base
    at which "3'-end-like" begins — this is a compositional zone, not a motif
    with ends. Pass ``step=1`` for base-resolution, rotation-exact boundaries
    at `step`-fold the work.

    Calibration: AUC 0.895 against coding sequence, 0.575 against plant
    genomic DNA. The honest claim is "this does not look like coding sequence",
    not "this does not look like plant DNA" — plant intergenic DNA is AT-rich
    too, so a hit inside a backbone region means much less than one inside a
    CDS you are designing.
    """
    s = _clean(seq)
    if not s:
        return []
    win = _3P_WINDOW if window is None else max(20, int(window))
    stp = _3P_STEP if step is None else max(1, int(step))
    floor = _3P_MIN_SCORE if min_score is None else float(min_score)
    n = len(s)
    if n < win:
        return []
    scan = _scan_window(s, circular, win - 1)
    last = n if circular else n - win + 1
    flagged: "list[tuple[int, int]]" = []
    best_at: "dict[int, tuple[float, float, float, int]]" = {}
    for i in range(0, last, stp):
        w = scan[i:i + win]
        if len(w) < win:
            break
        score, at, t20 = _3prime_score(w)
        strand = 1
        if both_strands:
            # AT fraction is strand-symmetric but T-richness is not: a 3' end
            # on the reverse strand reads A-rich here. Scoring the complement
            # too is what stops a reverse-oriented cassette's terminator from
            # being invisible.
            r_score, r_at, r_t20 = _3prime_score(_rc(w))
            if r_score > score:
                score, at, t20, strand = r_score, r_at, r_t20, -1
        if score < floor:
            continue
        flagged.append((i, win))
        best_at[i] = (score, at, t20, strand)
    out: "list[dict]" = []
    for start, length in _merge_arcs(flagged, n, circular):
        # The merged region's score is the best window inside it: the user is
        # deciding whether to recode this stretch, and the peak is what makes
        # that call.
        members = [v for i, v in best_at.items()
                   if _in_arc(i, start, length, n, circular)]
        if not members:
            continue
        score, at, t20, strand = max(members, key=lambda m: m[0])
        out.append({
            "start": start,
            "end": ((start + length) % n) if circular else start + length,
            "length": length,
            "score": round(score, 3),
            "at_fraction": round(at, 3),
            "best_t_fraction": round(t20, 3),
            "strand": strand,
            "tier": ("strong" if score >= _3P_STRONG_SCORE else "moderate"),
            "model": (
                f"plant 3'-end-like composition "
                f"({_3P_AT_WEIGHT:g}*AT + {_3P_T20_WEIGHT:g}*best T"
                f"{_3P_T20_LEN}); AUC {_3P_AUC_VS_CODING} vs coding, "
                f"{_3P_AUC_VS_GENOMIC} vs plant genomic. Region only -- "
                f"no cleavage site is predicted."),
        })
    out.sort(key=lambda h: (-h["score"], h["start"]))
    return out


def _in_arc(pos: int, start: int, length: int, total: int,
            circular: bool) -> bool:
    """Is scan index `pos` inside the arc ``(start, length)``?"""
    if not circular or total <= 0:
        return start <= pos < start + length
    return ((pos - start) % total) < length


def _regulettuce_calibration() -> dict:
    """Every measurement behind this module, and every model it refuses.

    Surfaced to the user (and to the agent API) so the numbers are not folded
    into source comments only — the same job `_splice_model_summary` does for
    the splice matrices. ``not_modelled`` is deliberately part of the payload:
    an annotator's refusals are as informative as its findings."""
    cat = _cis_element_catalogue()
    by_class: "dict[str, int]" = {}
    for row in cat:
        by_class[row["informativeness"]] = by_class.get(
            row["informativeness"], 0) + 1
    return {
        "reference_sets": {
            "plant_promoters": {
                "n": _CIS_N_PROMOTERS,
                "window_bp": _CIS_WINDOW_BP,
                "source": ("TAIR/RefSeq mRNA features over 8 Mb of "
                           "Arabidopsis thaliana chromosomes 1, 3 and 5, "
                           "keeping only transcripts with an annotated 5'UTR "
                           "so the 5' end is a transcription start and not a "
                           "start codon"),
            },
            "null_windows": {
                "n": _CIS_N_NULL,
                "source": ("composition-matched windows from the same 8 Mb of "
                           "genome -- real plant DNA, not shuffled bases"),
            },
            "polya_sites": {
                "n": 600,
                "source": ("3' ends of Arabidopsis thaliana and Oryza sativa "
                           "RefSeq mRNAs, poly(A) tail trimmed"),
            },
        },
        "cis_elements": {
            "total": len(cat),
            "by_informativeness": by_class,
            "enriched_min": _CIS_ENRICHED_MIN,
        },
        "tata": {
            "consensus": _TATA_CONSENSUS,
            "peak_window": list(_TATA_PEAK),
            "reported_window": list(_TATA_WINDOW),
            "peak_enrichment": _TATA_PEAK_ENRICHMENT,
            "genomic_bp_per_hit": _TATA_GENOMIC_BP_PER_HIT,
            "requires_anchor": True,
        },
        "three_prime": {
            "window": _3P_WINDOW,
            "step": _3P_STEP,
            "weights": {"at_fraction": _3P_AT_WEIGHT,
                        "best_t_fraction": _3P_T20_WEIGHT},
            "auc_vs_coding": _3P_AUC_VS_CODING,
            "auc_vs_genomic": _3P_AUC_VS_GENOMIC,
            "default_threshold": _3P_MIN_SCORE,
            "sensitivity_at_default": _3P_SENSITIVITY_AT_DEFAULT,
            "coding_fpr_at_default": _3P_CODING_FPR_AT_DEFAULT,
            "predicts_cleavage_site": False,
        },
        "not_modelled": [
            {"what": "poly(A) cleavage-site position",
             "why": ("a position-specific log-odds model over 600 real sites "
                     "located them below chance: argmax within +-20 nt in "
                     "27.2% of held-out 3' regions against 29.1% expected, "
                     "and its window AUC (0.864) lost to a composition-only "
                     "model (0.873)")},
            {"what": "plant promoter calling / strength",
             "why": ("a TATA + pyrimidine + Initiator + GC composite reached "
                     "AUC 0.646 against the composition-matched null; at a "
                     "background of 1 hit per 95 windows its sensitivity was "
                     "5.0%")},
            {"what": "Y-patch (plant pyrimidine patch)",
             "why": ("a >=8 nt C/T run appears in 40.6% of real promoter "
                     "windows and 38.6% of null windows -- 1.05x, no signal "
                     "to report")},
        ],
        "reproduce_with": "scripts/calibrate_regulettuce.py",
    }
