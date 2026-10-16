"""
test_introns — intron awareness through the whole annotation chain.

Guards `[INV-221]`. Before it, the engine's flat feature dict carried only
``(start, end)``, so a SPLICED feature was reduced to its outer bounds the
moment it entered a digest → ligate: a two-exon CDS came out of
``carry_annotations`` annotated straight through its own intron. That is not
a cosmetic error — the carried annotation translated to
``MAEAKGFEVSLI*FFQRLDLKAW*`` where the parent reads ``MAEAKGFERLDLKAW*``,
i.e. a premature stop, presented as the product's own CDS.

The contract these tests pin:

  * ``start``/``end`` keep their old meaning — outer traversal bounds, with
    ``end < start`` for an origin wrap.
  * an optional ``parts`` list carries the exons in forward-traversal order,
    and is present ONLY when the feature genuinely has introns (exonic length
    strictly less than the traversal span). **Absent means contiguous.**
  * the oracle is always the BASES or the PROTEIN, never the coordinates —
    a feature's identity is what it extracts, which is the lesson that found
    the bug in the first place.

Biopython's own ``SeqFeature.extract()`` is used as the outside reference
wherever a reading is checked.
"""
from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import (CompoundLocation, FeatureLocation, SeqFeature)

import splicecraft as sc
import splicecraft_agent as sca
import splicecraft_biology as scb
import splicecraft_cloning as scc
import splicecraft_dataaccess as scd
import splicecraft_util as su


TOTAL = 400
_RNG = random.Random(99)
SEQ = "".join(_RNG.choice("ACGT") for _ in range(TOTAL))


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _mk(spans, strand, total=TOTAL, ftype="CDS", label="t"):
    """A feature stored the way a GenBank file stores it: parts in READING
    order, which for the minus strand is the reverse of traversal order."""
    locs = [FeatureLocation(a, b, strand=strand) for a, b in spans]
    loc = (locs[0] if len(locs) == 1
           else CompoundLocation(list(reversed(locs)) if strand == -1
                                 else locs))
    return SeqFeature(loc, type=ftype, qualifiers={"label": [label]})


def _expect(spans, strand, seq=SEQ):
    """Independent oracle: the spans concatenated in traversal order."""
    fwd = "".join(seq[a:b] for a, b in spans)
    return _rc(fwd) if strand == -1 else fwd


# Every geometry that matters, as span lists in forward-traversal order.
GEOMETRIES = [
    ("contiguous",            [(50, 110)],                          False),
    ("spliced-2",             [(50, 80), (120, 160)],               True),
    ("spliced-3",             [(50, 80), (120, 160), (200, 215)],   True),
    ("spliced-4", [(10, 25), (50, 80), (120, 160), (300, 321)],     True),
    ("tiny-exons",            [(50, 53), (60, 63)],                 True),
    ("adjacent-exons",        [(50, 80), (80, 110)],                False),
    ("wrap",                  [(360, TOTAL), (0, 40)],              False),
    ("wrap+intron-in-head",   [(360, TOTAL), (20, 40)],             True),
    ("wrap+intron-in-tail",   [(300, 330), (360, TOTAL), (0, 40)],  True),
    ("wrap+2-introns",
     [(300, 330), (360, TOTAL), (0, 20), (40, 60)],                 True),
    ("ends-on-origin",        [(360, TOTAL)],                       False),
    ("starts-at-origin",      [(0, 40)],                            False),
    ("whole-molecule",        [(0, TOTAL)],                         False),
]


# ═══════════════════════════════════════════════════════════════════════════
# The L0 contract
# ═══════════════════════════════════════════════════════════════════════════

def test_util_span_len_matches_feat_len():
    """`_span_len_wrapaware` MIRRORS biology's `_feat_len` because util
    cannot import biology (that would be a cycle). Drift guard — if the two
    ever disagree the wrap arithmetic has forked."""
    for total in (1, 2, 7, 400, 2701):
        for start in range(0, total, max(1, total // 11)):
            for end in range(0, total, max(1, total // 13)):
                assert (su._span_len_wrapaware(start, end, total)
                        == sc._feat_len(start, end, total)), (
                    start, end, total)


@pytest.mark.parametrize("name,spans,spliced",
                         [(n, s, sp) for n, s, sp in GEOMETRIES])
@pytest.mark.parametrize("strand", (1, -1))
def test_feat_spans_reports_introns_only_when_present(name, spans, spliced,
                                                      strand):
    """`parts` is emitted exactly when the feature has introns — decided by
    MEASUREMENT (exonic < span), not by `len(parts) > 1`, because an origin
    wrap is a two-part location too (sacred #9)."""
    got = su._feat_spans(_mk(spans, strand), TOTAL, circular=True)
    assert got is not None, name
    _s, _e, _st, parts = got
    assert bool(parts) is spliced, f"{name}/{strand}: parts={parts}"
    if spliced:
        assert parts is not None
        exonic = sum(b - a for a, b in parts)
        assert exonic < su._span_len_wrapaware(_s, _e, TOTAL)


@pytest.mark.parametrize("name,spans,_spliced", GEOMETRIES)
@pytest.mark.parametrize("strand", (1, -1))
def test_spans_round_trip_the_bases(name, spans, _spliced, strand):
    """`_feat_spans` → `_feature_location_spans` → the SAME bases Biopython
    extracts. The coordinates being tidy is not the point."""
    want = _expect(spans, strand)
    feat = _mk(spans, strand)
    assert str(feat.extract(Seq(SEQ))).upper() == want, "oracle disagrees"
    s, e, st, parts = su._feat_spans(feat, TOTAL, circular=True)
    loc = (su._feature_location_spans(parts, TOTAL, st) if parts
           else su._feature_location(s, e, TOTAL, st))
    assert loc is not None
    assert str(SeqFeature(loc, type="CDS").extract(Seq(SEQ))).upper() == want


def test_feature_location_spans_stores_minus_strand_exons_descending():
    """Stored order IS reading order ([INV-182]): a minus-strand multi-exon
    location must store its parts DESCENDING, or every other tool reads the
    exons in the wrong order — a cyclically rotated protein."""
    loc = su._feature_location_spans([[50, 80], [120, 160]], TOTAL, -1)
    starts = [int(p.start) for p in loc.parts]
    assert starts == [120, 50], starts
    fwd = su._feature_location_spans([[50, 80], [120, 160]], TOTAL, 1)
    assert [int(p.start) for p in fwd.parts] == [50, 120]


def test_dict_spans_order_distinguishes_a_wrap_from_a_whole_span_gene():
    """Span ORDER is information. `[[0,40],[360,400]]` is either one 80 bp
    feature crossing the origin or a two-exon gene with a 320 bp intron, and
    only the traversal order says which — so `_feat_dict_set_spans` must not
    sort."""
    wrap = su._feat_dict_set_spans({}, [[360, TOTAL], [0, 40]], TOTAL)
    assert (wrap["start"], wrap["end"]) == (360, 40)
    assert "parts" not in wrap, "a plain wrap must not claim an intron"

    gene = su._feat_dict_set_spans({}, [[0, 40], [360, TOTAL]], TOTAL)
    assert (gene["start"], gene["end"]) == (0, TOTAL)
    assert gene["parts"] == [[0, 40], [360, TOTAL]]


def test_dict_set_spans_drops_parts_when_exons_become_contiguous():
    """ABSENT is the contract for "no introns"; a one-span `parts` left
    behind would be a second, disagreeing truth."""
    fd = su._feat_dict_set_spans({"parts": [[10, 20], [20, 30]]},
                                 [[10, 20], [20, 30]], TOTAL)
    assert (fd["start"], fd["end"]) == (10, 30)
    assert "parts" not in fd


@pytest.mark.parametrize("bad", [
    "nonsense", 5, {}, [], [[1, 2]], [[1]], [["a", "b"], ["c", "d"]],
    [[5, 5], [9, 9]], [None, None],
])
def test_dict_spans_never_raises_on_malformed_parts(bad):
    """A hand-edited GenBank or a caller-built list must degrade to the
    contiguous reading, never raise out of a digest and never return []."""
    fd = {"start": 10, "end": 50, "parts": bad}
    got = su._feat_dict_spans(fd, TOTAL)
    assert got, "must never return an empty span list"
    assert got == [[10, 50]]


def test_location_for_dict_falls_back_when_parts_unusable():
    assert su._feature_location_for_dict(
        {"start": 10, "end": 50, "parts": "junk"}, TOTAL, 1) is not None
    loc = su._feature_location_for_dict(
        {"start": 10, "end": 50, "parts": [[10, 20], [30, 50]]}, TOTAL, 1)
    assert len(loc.parts) == 2


# ═══════════════════════════════════════════════════════════════════════════
# The engine: a spliced gene keeps its exons through a real reaction
# ═══════════════════════════════════════════════════════════════════════════

EX1 = "ATGGCTGAAGCTAAGGGTTTCGAA"
INTRON = "GTAAGTTTGATCTAATTTTTTCAG"      # carries stops, so a read-through
EX2 = "CGTCTGGATCTGAAAGCTTGGTAA"          # is unmistakable in the protein
PROTEIN = "MAEAKGFERLDLKAW*"


@pytest.fixture(scope="module")
def spliced_plasmid():
    """A circular plasmid carrying one two-exon CDS, plus its coordinates."""
    rng = random.Random(7)
    pad = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    pre, post = pad(200), pad(600)
    seq = pre + EX1 + INTRON + EX2 + post
    a = len(pre)
    b = a + len(EX1)
    c = b + len(INTRON)
    d = c + len(EX2)
    rec = SeqRecord(Seq(seq), id="P", name="P",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"})
    rec.features.append(SeqFeature(
        CompoundLocation([FeatureLocation(a, b, 1), FeatureLocation(c, d, 1)]),
        type="CDS", qualifiers={"label": ["myGene"]}))
    assert str(rec.features[0].extract(rec.seq).translate()) == PROTEIN
    return rec, seq, (a, b, c, d)


def _protein_of(feat_dict, prod_seq):
    """Translate what a product's feature dict actually annotates."""
    spans = su._feat_dict_spans(feat_dict, len(prod_seq))
    fwd = "".join(prod_seq[s:e] for s, e in spans)
    dna = sc._rc(fwd) if int(feat_dict.get("strand", 1) or 1) == -1 else fwd
    return str(Seq(dna).translate())


def test_carry_dicts_emit_the_exons(spliced_plasmid):
    rec, seq, (a, b, c, d) = spliced_plasmid
    dicts = sca._agent_carry_feature_dicts(rec)
    gene = next(f for f in dicts if f["label"] == "myGene")
    assert gene["parts"] == [[a, b], [c, d]]
    exonic = sum(e - s for s, e in su._feat_dict_spans(gene, len(seq)))
    assert exonic == len(EX1) + len(EX2)


def test_digest_then_religate_keeps_the_protein(spliced_plasmid):
    """The whole point: cut the plasmid outside the gene, put it back, and
    the gene must still translate to the same protein."""
    rec, seq, (a, b, c, d) = spliced_plasmid
    dicts = sca._agent_carry_feature_dicts(rec)
    enz = next((e for e in ("EcoRI", "BamHI", "HindIII", "PstI", "SalI",
                            "XbaI", "NcoI", "KpnI", "SacI", "XhoI")
                if len(sc._enzyme_cuts(seq, [e], circular=True)) == 1
                and not (a <= sc._enzyme_cuts(seq, [e],
                                              circular=True)[0]["top"] < d)),
               None)
    assert enz is not None, "fixture needs a single cutter outside the gene"
    frags = sc._digest_with_enzymes(seq, [enz], circular=True,
                                    features=dicts)
    closed = sc._close_circular(frags[0])
    assert closed is not None
    feats = scc._rejoin_origin_split_features(closed["features"],
                                              len(closed["top_seq"]))
    feats = scc._rejoin_cut_split_features(feats, len(closed["top_seq"]),
                                           top_seq=closed["top_seq"])
    gene = next(f for f in feats if "myGene" in str(f.get("label")))
    assert gene.get("parts"), "the exon list did not survive the digest"
    assert _protein_of(gene, closed["top_seq"]) == PROTEIN


def test_golden_gate_one_pot_keeps_the_protein(spliced_plasmid):
    rec, seq, (a, b, c, d) = spliced_plasmid
    rng = random.Random(21)
    pad = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    payload = seq[a - 20:d + 20]
    part = ("AAAA" + "GGTCTC" + "A" + "AATG" + payload + "GCTT" + "T"
            + "GAGACC" + "AAAA")
    g = part.find(EX1)                      # PART coordinates, per contract
    assert g > 0
    pf = [{"start": g, "end": g + (d - a), "strand": 1, "type": "CDS",
           "label": "myGene",
           "parts": [[g, g + (b - a)], [g + (c - a), g + (d - a)]]}]
    acc = ("AATG" + "T" + "GAGACC" + pad(300) + "GGTCTC" + "A" + "GCTT"
           + pad(1500))
    r = scc._simulate_golden_gate([part], acc, enzyme="BsaI",
                                  part_features=[pf])
    assert r.get("ok"), r.get("errors")
    gene = next(f for f in r["features"] if "myGene" in str(f.get("label")))
    assert gene.get("parts"), "exons lost in the one-pot assembly"
    assert _protein_of(gene, r["product_seq"]) == PROTEIN


def test_product_record_writes_a_join(spliced_plasmid):
    """The exon list has to reach the FILE, or the next reader flattens it
    again. A `join()` read back by Biopython must give the same protein."""
    rec, seq, (a, b, c, d) = spliced_plasmid
    fd = {"start": a, "end": d, "strand": 1, "type": "CDS", "label": "g",
          "parts": [[a, b], [c, d]]}
    loc = su._feature_location_for_dict(fd, len(seq), 1)
    assert len(loc.parts) == 2
    assert str(SeqFeature(loc, type="CDS").extract(Seq(seq)).translate()) \
        == PROTEIN


def test_rc_fragment_mirrors_the_exons():
    """A reverse-orientation clone flips the molecule: the exons ride the
    same affine map as the bounds AND their order inverts."""
    seq = SEQ[:120]
    frag = {"top_seq": seq,
            "left":  {"overhang_seq": "AATT", "kind": "5'", "enzyme": "E"},
            "right": {"overhang_seq": "AATT", "kind": "5'", "enzyme": "E"},
            "features": [{"start": 10, "end": 70, "strand": 1,
                          "type": "CDS", "label": "g",
                          "parts": [[10, 30], [50, 70]]}],
            "source_label": ""}
    flipped = scc._rc_fragment(frag)
    g = flipped["features"][0]
    assert g["strand"] == -1
    before = "".join(seq[s:e] for s, e in [[10, 30], [50, 70]])
    after_spans = su._feat_dict_spans(g, len(flipped["top_seq"]))
    after = "".join(flipped["top_seq"][s:e] for s, e in after_spans)
    assert after == _rc(before), (before, after)


# ═══════════════════════════════════════════════════════════════════════════
# The three adjacent defects this pass also fixed
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("name,spans,strand", [
    ("plus-spliced",  [(10, 20), (30, 40)],        1),
    ("plus-wrap",     [(290, 300), (0, 10)],       1),
    ("minus-spliced", [(10, 20), (30, 40)],       -1),
    ("minus-wrap",    [(290, 300), (0, 10)],      -1),
])
def test_feature_library_export_matches_biopython(name, spans, strand):
    """`_extract_feature_entries_from_record` feeds "Import from library
    plasmid", so a wrong reading is SAVED into the user's feature library.
    It read `loc.parts` in STORED order, which on the minus strand is
    DESCENDING — the right bases in the wrong block order, i.e. a cyclically
    rotated gene. Oracle: Biopython's own extract()."""
    seq = SEQ[:300]
    rec = SeqRecord(Seq(seq), id="X", name="X",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"})
    rec.features = [_mk(spans, strand, total=300, label=name)]
    want = str(rec.features[0].extract(rec.seq)).upper()
    got = scd._extract_feature_entries_from_record(rec)[0]["sequence"].upper()
    assert got == want, f"{name}: {got} != {want}"


def test_gibson_reports_a_dropped_malformed_feature():
    """A Gibson fragment is LINEAR, so `end <= start` cannot be a wrap — it
    is malformed input. Dropping it is right; dropping it in SILENCE meant a
    saved product lost the annotation with nothing in the response."""
    ov = "GGTTCCAAGGTTCCAAGGTTCCAAG"
    rng = random.Random(3)
    core = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    a = ov + core(200) + ov
    b = ov + core(200) + ov
    frags = [{"sequence": a, "name": "A", "features": [
                 {"start": len(a) - 10, "end": 5, "strand": 1,
                  "type": "misc_feature", "label": "wrapped"},
                 {"start": 40, "end": 80, "strand": 1,
                  "type": "misc_feature", "label": "plain"}]},
             {"sequence": b, "name": "B", "features": []}]
    r = scc._simulate_gibson_assembly(frags, circular=True)
    labels = [f.get("label") for f in (r.get("features") or [])]
    assert "plain" in labels
    assert "wrapped" not in labels
    assert any("wrapped" in w for w in (r.get("warnings") or [])), \
        "the drop must be reported, not silent"


def test_protein_tab_excludes_features_outside_a_wrapping_cds():
    """The synthesis protein tab read the CDS's FLATTENED bounds, so for a
    wrap CDS the "inside the CDS" test spanned the whole molecule and an
    unrelated backbone feature was imported as a motif at an invented
    residue."""
    seq = SEQ[:600]
    rec = SeqRecord(Seq(seq), id="P", name="P",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"})
    cds = SeqFeature(CompoundLocation([FeatureLocation(540, 600, 1),
                                       FeatureLocation(0, 60, 1)]),
                     type="CDS", qualifiers={"label": ["gene"],
                                             "translation": ["M" * 40]})
    def _mf(s, e, nm):
        return SeqFeature(FeatureLocation(s, e, 1), type="misc_feature",
                          qualifiers={"label": [nm], "note": [
                              "splicecraft-aa-feature-type=Tag"]})
    rec.features = [cds, _mf(540, 570, "realTag"), _mf(300, 330, "far")]
    got = {f["label"]: (f["start"], f["end"])
           for f in sc.SynthesisScreen._extract_aa_feats_from_record(rec)}
    assert got.get("realTag") == (0, 10), got
    assert "far" not in got, "an off-CDS feature was imported as a motif"


def test_protein_tab_numbers_residues_through_the_intron():
    """A motif on exon 2 is numbered from the MESSAGE, not the DNA."""
    seq = SEQ[:600]
    rec = SeqRecord(Seq(seq), id="P", name="P",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"})
    rec.features = [
        SeqFeature(CompoundLocation([FeatureLocation(0, 60, 1),
                                     FeatureLocation(300, 360, 1)]),
                   type="CDS", qualifiers={"label": ["g"],
                                           "translation": ["M" * 40]}),
        SeqFeature(FeatureLocation(300, 330, 1), type="misc_feature",
                   qualifiers={"label": ["exon2Tag"], "note": [
                       "splicecraft-aa-feature-type=Tag"]}),
        SeqFeature(FeatureLocation(120, 150, 1), type="misc_feature",
                   qualifiers={"label": ["intronic"], "note": [
                       "splicecraft-aa-feature-type=Tag"]}),
    ]
    got = {f["label"]: (f["start"], f["end"])
           for f in sc.SynthesisScreen._extract_aa_feats_from_record(rec)}
    assert got.get("exon2Tag") == (20, 30), got
    assert "intronic" not in got, "an intronic motif has no residue number"


# ═══════════════════════════════════════════════════════════════════════════
# Robustness — the digest must survive anything a caller hands it
# ═══════════════════════════════════════════════════════════════════════════

def test_digest_ignores_non_dict_feature_entries():
    """Handing ONE feature list where a list PER PART is wanted used to
    surface as `'str' object has no attribute 'get'` from inside a digest —
    an AttributeError the `(TypeError, ValueError)` guards never caught."""
    feats = ["nonsense", 42, None, {"start": 10, "end": 40, "strand": 1,
                                    "type": "misc_feature", "label": "ok"}]
    frags = sc._digest_with_enzymes(SEQ, ["EcoRI"], circular=True,
                                    features=feats)
    labels = [f.get("label") for fr in frags for f in fr["features"]]
    assert "ok" in labels


def test_split_features_at_cuts_tolerates_malformed_parts():
    out = scb._split_features_at_cuts(
        [{"start": 10, "end": 300, "strand": 1, "type": "CDS",
          "label": "g", "parts": [["x", "y"], [1]]}],
        TOTAL, [100, 200], circular=True)
    assert out, "the feature must still be slotted"
    for pieces in out.values():
        for p in pieces:
            assert "parts" not in p or len(p["parts"]) > 1


def test_overlapping_parts_read_as_contiguous():
    """Malformed overlapping exons measure exonic >= span, so they read as
    contiguous — the outer bounds, which is what `_feat_bounds` always gave.
    Deliberately fail-safe rather than inventing a structure."""
    feat = _mk([(50, 120), (60, 140)], 1)
    s, e, st, parts = su._feat_spans(feat, TOTAL, circular=True)
    assert parts is None


@pytest.mark.parametrize("name,spans,spliced", GEOMETRIES)
@pytest.mark.parametrize("strand", (1, -1))
def test_agent_feature_parts_agrees_with_feat_spans(name, spans, spliced,
                                                    strand):
    """Two expressions of one rule — "exons only when genuinely
    discontiguous" — must never drift. `_agent_feature_parts` (hub, serves
    `find-feature` / `_agent_feature_detail`) decides it by a CONTIGUITY
    walk; `_feat_spans` (util L0, serves the cloning engine) decides it by
    MEASURING exonic length against the traversal span. They have to agree
    on every geometry, or the same feature reads as spliced through one door
    and contiguous through the other."""
    rec = SeqRecord(Seq(SEQ), id="X", name="X",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"})
    feat = _mk(spans, strand)
    rec.features = [feat]
    got = su._feat_spans(feat, TOTAL, circular=True)
    assert got is not None
    _s, _e, _st, parts = got
    row = {"wraps": _e < _s}
    hub = sc._agent_feature_parts(feat, row, rec)
    assert bool(hub) is bool(parts) is spliced, (
        f"{name}/{strand}: hub={hub} util={parts}")
    if spliced:
        assert [list(p) for p in (hub or [])] == [list(p) for p in
                                                  (parts or [])]


@pytest.mark.parametrize("strand", (1, -1))
def test_new_plasmid_can_author_a_spliced_feature(strand):
    """The WRITE direction: `new-plasmid` takes `parts` so a gene with
    introns can be authored from outside the GUI, and it must read back with
    its exons — not as one block over them."""
    bases = SEQ
    spans = [[50, 80], [120, 160]]
    rec = sc.NewPlasmidModal._build_record(
        None, "spliced test", bases, True,
        [{"start": 50, "end": 160, "strand": strand, "name": "g",
          "parts": spans}])
    feat = next(f for f in rec.features if f.type != "source")
    assert len(feat.location.parts) == 2
    want = _expect([(50, 80), (120, 160)], strand)
    assert str(feat.extract(rec.seq)).upper() == want
    # and the documented `end > n` wrap convention still works beside it
    rec2 = sc.NewPlasmidModal._build_record(
        None, "wrap test", bases, True,
        [{"start": TOTAL - 20, "end": TOTAL + 20, "strand": 1, "name": "w"}])
    w = next(f for f in rec2.features if f.type != "source")
    assert len(w.location.parts) == 2


@pytest.mark.parametrize("strand", (1, -1))
def test_fuzz_spans_round_trip(strand):
    """Random exon sets, including orders that cross the origin mid
    traversal, must round-trip their bases."""
    rng = random.Random(4242 + strand)
    for _ in range(600):
        n_ex = rng.randint(1, 5)
        pts = sorted(rng.sample(range(TOTAL + 1), n_ex * 2))
        spans = [(pts[2 * i], pts[2 * i + 1]) for i in range(n_ex)]
        spans = [(a, b) for a, b in spans if b > a]
        if not spans:
            continue
        k = rng.randrange(len(spans))
        spans = spans[k:] + spans[:k]
        want = _expect(spans, strand)
        loc = su._feature_location_spans(list(spans), TOTAL, strand)
        assert loc is not None
        assert str(SeqFeature(loc, type="CDS").extract(Seq(SEQ))).upper() \
            == want, spans
        s, e, st, parts = su._feat_spans(
            SeqFeature(loc, type="CDS"), TOTAL, circular=True)
        back = (su._feature_location_spans(parts, TOTAL, st) if parts
                else su._feature_location(s, e, TOTAL, st))
        assert back is not None
        assert str(SeqFeature(back, type="CDS").extract(Seq(SEQ))).upper() \
            == want, (spans, parts)
