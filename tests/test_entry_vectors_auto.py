"""Entry-vector auto-detection + EntryVectorsModal.

Detection engine: `_detect_entry_vector_role` digests a candidate
plasmid with both grammar enzymes, identifies the stuffer fragment via
backbone-marker exclusion (NEVER size — feedback
`feedback_never_assume_smaller_frag_is_payload`), and resolves
(inner-enzyme, outer-overhang-pair) → role via the Golden Braid binary-
assembly architecture.

Auto-bind helper: `_auto_bind_entry_vectors_from_entries` runs
detection across every grammar for every entry and binds new roles,
without clobbering existing user bindings.

Modal: `EntryVectorsModal` exposes one row per role with Pick/Clear/
Auto-detect controls. Replaces the single-slot widget in
`GrammarEditorModal`.
"""
import pytest

import splicecraft as sc

pytestmark = [pytest.mark.usefixtures("_protect_user_data")]


# ── Synthetic-acceptor builder ─────────────────────────────────────────────

def _build_acceptor(
    *, oh5_inner: str, oh3_inner: str,
    oh5_outer: str, oh3_outer: str,
    inner_enzyme_site: str = "GGTCTC",    # BsaI forward
    outer_enzyme_site: str = "CGTCTC",    # Esp3I forward
    backbone_body: str = "ATCGATCGAT" * 100,
    stuffer_body: str  = "TTGGAACCAA" * 20,
):
    """Build a circular acceptor where:

      * The INNER enzyme (default BsaI) cuts to release a stuffer with
        `(oh5_inner, oh3_inner)`.
      * The OUTER enzyme (default Esp3I) cuts at outer sites to release
        a larger fragment containing the inner sites + stuffer, with
        outer overhangs `(oh5_outer, oh3_outer)`.

    The backbone carries a `rep_origin` feature so
    `_fragment_has_backbone_marker` can identify it.

    Layout (linearised for clarity, actually circular):

        [backbone+rep_origin] —[outer-rc]——[inner-fwd]·oh5_inner·[stuffer]·oh3_inner·[inner-rc]——[outer-fwd]——
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    # Forward Esp3I-type site: ENZ + N + overhang + payload
    # Reverse site:           payload + overhang + N + RC(ENZ)
    inner_fwd = inner_enzyme_site + "A"
    inner_rc  = "T" + sc._rc(inner_enzyme_site)
    outer_fwd = outer_enzyme_site + "A"
    outer_rc  = "T" + sc._rc(outer_enzyme_site)

    # Inner cassette (between outer cuts) contains the stuffer +
    # flanking BsaI sites.
    inner_cassette = (
        inner_fwd + oh5_inner + stuffer_body + oh3_inner + inner_rc
    )
    # Outer "release" region surrounds the inner cassette. We arrange
    # outer Esp3I sites OUTSIDE the inner cassette so that an Esp3I
    # digest of the full ring releases the [inner_cassette] piece
    # with outer overhangs (oh5_outer, oh3_outer).
    inner_with_outer_flanks = (
        outer_fwd + oh5_outer + inner_cassette + oh3_outer + outer_rc
    )
    seq = inner_with_outer_flanks + backbone_body
    backbone_start = len(inner_with_outer_flanks)
    rec = SeqRecord(
        Seq(seq), id="acceptor", name="acceptor",
        annotations={"topology": "circular", "molecule_type": "DNA"},
    )
    rec.features.append(SeqFeature(
        FeatureLocation(backbone_start + 100, backbone_start + 500),
        type="rep_origin",
        qualifiers={"label": ["bla"]},
    ))
    return rec


# ── Detection engine: canonical α/Ω matches ────────────────────────────────

def test_alpha1_detected_strict():
    """α1: inner BsaI release = (GGAG, CGCT); outer Esp3I = (GGAG, GTCA)."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
        inner_enzyme_site="GGTCTC",   # BsaI
        outer_enzyme_site="CGTCTC",   # Esp3I
    )
    result = sc._detect_entry_vector_role(rec, gb_l0)
    assert result == ("Alpha1", "strict"), (
        f"Expected Alpha1/strict, got {result}"
    )


def test_alpha2_detected_strict():
    """α2: inner BsaI = (GGAG, CGCT); outer Esp3I = (GTCA, CGCT)."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GTCA", oh3_outer="CGCT",
        inner_enzyme_site="GGTCTC",
        outer_enzyme_site="CGTCTC",
    )
    result = sc._detect_entry_vector_role(rec, gb_l0)
    assert result == ("Alpha2", "strict")


def test_omega1_detected_strict():
    """Ω1: inner Esp3I = (GGAG, CGCT); outer BsaI = (GGAG, GTCA)."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
        inner_enzyme_site="CGTCTC",   # Esp3I inner = Ω
        outer_enzyme_site="GGTCTC",   # BsaI outer = Ω
    )
    result = sc._detect_entry_vector_role(rec, gb_l0)
    assert result == ("Omega1", "strict")


def test_omega2_detected_strict():
    """Ω2: inner Esp3I = (GGAG, CGCT); outer BsaI = (GTCA, CGCT)."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GTCA", oh3_outer="CGCT",
        inner_enzyme_site="CGTCTC",
        outer_enzyme_site="GGTCTC",
    )
    result = sc._detect_entry_vector_role(rec, gb_l0)
    assert result == ("Omega2", "strict")


# ── Detection: UPD weak match ──────────────────────────────────────────────

def test_upd_style_singleton_detected_weak():
    """A plasmid with non-canonical stuffer overhangs (not L0 positions,
    not in the canonical alphabet) but still digesting cleanly should
    detect as a singleton L0 donor ("" role, weak)."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    # `TGAG`/`CTCA` matches FFE1 UPD's BsaI digest — non-canonical.
    rec = _build_acceptor(
        oh5_inner="TGAG", oh3_inner="CTCA",
        oh5_outer="CTCG", oh3_outer="TGAG",
        inner_enzyme_site="GGTCTC",
        outer_enzyme_site="CGTCTC",
    )
    result = sc._detect_entry_vector_role(rec, gb_l0)
    assert result is not None
    role, conf = result
    assert role == ""
    assert conf == "weak"


# ── Detection: rejection cases ─────────────────────────────────────────────

def test_l0_part_not_detected_as_acceptor():
    """An L0 part has overhangs matching an L0 position (e.g. Promoter
    = (GGAG, AATG)). The detection engine must REJECT it — it's a
    part, not an empty acceptor."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="AATG",     # Pos 1 Promoter
        oh5_outer="GGAG", oh3_outer="AATG",
        inner_enzyme_site="GGTCTC",
        outer_enzyme_site="CGTCTC",
    )
    result = sc._detect_entry_vector_role(rec, gb_l0)
    assert result is None, (
        "L0 part overhangs must not be classified as a UPD donor"
    )


def test_tu_plasmid_not_detected_as_acceptor():
    """A TU plasmid's stuffer overhangs come from the canonical
    alphabet (matches `_classify_part_from_plasmid` pass-4). The
    detector must REJECT — TUs aren't empty acceptors."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    # Match DEMO 26 TUx1 release: Esp3I gives (GGAG, GTCA).
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="GTCA",
        oh5_outer="GGAG", oh3_outer="GTCA",
        inner_enzyme_site="GGTCTC",
        outer_enzyme_site="CGTCTC",
    )
    result = sc._detect_entry_vector_role(rec, gb_l0)
    assert result is None, "TU must not be classified as an acceptor"


def test_linear_record_not_detected():
    """Linear records can't be acceptors (no second cut, can't
    release a stuffer cleanly). Must return None."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
    )
    rec.annotations["topology"] = "linear"
    result = sc._detect_entry_vector_role(rec, gb_l0)
    assert result is None


def test_wrap_feature_backbone_marker_still_detected():
    """Sweep #35 (2026-05-26): a backbone marker that spans the origin
    (CompoundLocation: `join(tail..end, 0..head)`) used to silently
    flatten to `(0, total)` because `_detect_entry_vector_role`
    called `int(loc.start)` / `int(loc.end)` directly. BioPython
    returns `min(parts.start)=0` and `max(parts.end)=total` for a
    two-part wrap, so the marker appeared to cover the whole
    plasmid → `_fragment_has_backbone_marker` fired on BOTH digest
    fragments → `sum(marked) != 1` ambiguous → detection skipped.
    Post-fix the marker routes through `_feat_bounds`, which encodes
    a wrap as `end < start` and slots correctly into one fragment.
    """
    from Bio.SeqFeature import (SeqFeature, FeatureLocation,
                                  CompoundLocation)
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    gb_l0 = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
        inner_enzyme_site="GGTCTC",
        outer_enzyme_site="CGTCTC",
    )
    # Sanity baseline: the linear (non-wrap) version classifies cleanly.
    baseline = sc._detect_entry_vector_role(rec, gb_l0)
    assert baseline == ("Alpha1", "strict"), (
        f"Sanity check failed before wrap: got {baseline}"
    )
    # Rotate the sequence so the backbone's rep_origin straddles the
    # new origin (position 0). The biology is identical (circular
    # plasmids are rotation-invariant) but the rep_origin now lives
    # in a wrap-encoded CompoundLocation. Pick the rotation point
    # inside the rep_origin so the feature genuinely crosses 0.
    orig = next(
        (f for f in rec.features if f.type == "rep_origin"),
        None,
    )
    assert orig is not None, "_build_acceptor must place a rep_origin"
    orig_start = int(orig.location.start)
    orig_end   = int(orig.location.end)
    orig_mid   = (orig_start + orig_end) // 2
    total      = len(rec.seq)
    rot        = orig_mid
    rotated_seq = str(rec.seq)[rot:] + str(rec.seq)[:rot]
    new_rec = SeqRecord(
        Seq(rotated_seq), id=rec.id, name=rec.name,
        annotations=dict(rec.annotations),
    )
    # Shift each feature by `-rot`; features whose original interval
    # crossed `rot` become wrap CompoundLocations.
    for f in rec.features:
        fs, fe = int(f.location.start), int(f.location.end)
        new_fs = (fs - rot) % total
        new_fe = (fe - rot) % total
        if new_fs < new_fe:
            new_loc = FeatureLocation(new_fs, new_fe,
                                       strand=f.location.strand)
        else:
            # Wrap: tail [new_fs, total) + head [0, new_fe).
            new_loc = CompoundLocation([
                FeatureLocation(new_fs, total,
                                 strand=f.location.strand),
                FeatureLocation(0, new_fe,
                                 strand=f.location.strand),
            ])
        new_rec.features.append(SeqFeature(
            new_loc, type=f.type,
            qualifiers=dict(f.qualifiers),
        ))
    # Confirm the rep_origin is now a wrap (CompoundLocation) — if
    # it isn't, the test isn't exercising the bug we care about.
    new_rep = next(
        (f for f in new_rec.features if f.type == "rep_origin"),
        None,
    )
    assert new_rep is not None
    assert isinstance(new_rep.location, CompoundLocation), (
        "Test rotation did not produce a wrap-spanning rep_origin"
    )
    result = sc._detect_entry_vector_role(new_rec, gb_l0)
    assert result == ("Alpha1", "strict"), (
        f"Wrap-spanning backbone marker must still resolve to "
        f"Alpha1/strict; got {result}"
    )


def test_record_without_backbone_marker_rejected():
    """Sacred — backbone-marker exclusion is how the stuffer is
    identified. Without a marker, the detector can't safely pick
    a stuffer (size is forbidden). Returns None."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
    )
    # Strip the rep_origin feature.
    rec.features = []
    result = sc._detect_entry_vector_role(rec, gb_l0)
    assert result is None


# ── Auto-bind helper ───────────────────────────────────────────────────────

def _entry_from_record(rec, eid: str = "test") -> dict:
    """Wrap a SeqRecord into the library-entry dict shape."""
    return {
        "id":      eid,
        "name":    eid,
        "size":    len(rec.seq),
        "gb_text": sc._record_to_gb_text(rec),
        "n_feats": len(rec.features),
    }


def test_auto_bind_fills_all_alpha_omega_roles():
    """Feed all 4 acceptors + UPD-style donor to the bulk auto-bind;
    verify all 5 slots get filled with the correct role."""
    # Build 5 acceptors matching the FFE1-5 fingerprint.
    upd = _build_acceptor(
        oh5_inner="TGAG", oh3_inner="CTCA",
        oh5_outer="CTCG", oh3_outer="TGAG",
        inner_enzyme_site="GGTCTC", outer_enzyme_site="CGTCTC",
    )
    a1 = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
        inner_enzyme_site="GGTCTC", outer_enzyme_site="CGTCTC",
    )
    a2 = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GTCA", oh3_outer="CGCT",
        inner_enzyme_site="GGTCTC", outer_enzyme_site="CGTCTC",
    )
    o1 = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
        inner_enzyme_site="CGTCTC", outer_enzyme_site="GGTCTC",
    )
    o2 = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GTCA", oh3_outer="CGCT",
        inner_enzyme_site="CGTCTC", outer_enzyme_site="GGTCTC",
    )
    entries = [
        _entry_from_record(upd, "upd"),
        _entry_from_record(a1, "a1"),
        _entry_from_record(a2, "a2"),
        _entry_from_record(o1, "o1"),
        _entry_from_record(o2, "o2"),
    ]
    msg = sc._auto_bind_entry_vectors_from_entries(entries)
    assert msg, "Auto-bind should return a summary message"
    assert "Alpha1" in msg
    assert "Alpha2" in msg
    assert "Omega1" in msg
    assert "Omega2" in msg
    assert "UPD" in msg
    # Verify each binding lands in entry_vectors.json under the right
    # (grammar_id, role) key.
    assert sc._get_entry_vector("gb_l0", "")["id"] == "upd"
    assert sc._get_entry_vector("gb_l0", "Alpha1")["id"] == "a1"
    assert sc._get_entry_vector("gb_l0", "Alpha2")["id"] == "a2"
    assert sc._get_entry_vector("gb_l0", "Omega1")["id"] == "o1"
    assert sc._get_entry_vector("gb_l0", "Omega2")["id"] == "o2"


def test_auto_bind_does_not_clobber_existing_bindings():
    """Existing user-set bindings are sacred — the auto-bind must
    only fill in gaps, never replace a pre-existing binding."""
    # User has Alpha1 manually bound to "manual".
    sc._set_entry_vector("gb_l0", {
        "name": "manual", "size": 0, "gb_text": "", "id": "manual",
    }, "Alpha1")
    # Then auto-bind sees an Alpha1 candidate.
    a1 = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
        inner_enzyme_site="GGTCTC", outer_enzyme_site="CGTCTC",
    )
    sc._auto_bind_entry_vectors_from_entries(
        [_entry_from_record(a1, "auto")]
    )
    # The original "manual" binding survives — "auto" is skipped.
    bound = sc._get_entry_vector("gb_l0", "Alpha1")
    assert bound is not None
    assert bound["id"] == "manual", (
        "Auto-bind clobbered an existing user binding"
    )


def test_auto_bind_strict_wins_over_weak():
    """If two entries detect for the same role with different
    confidences, the strict match wins."""
    gb_l0 = sc._all_grammars()["gb_l0"]
    weak_upd = _build_acceptor(
        oh5_inner="TGAG", oh3_inner="CTCA",
        oh5_outer="CTCG", oh3_outer="TGAG",
    )
    strict_a1 = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
    )
    # Verify the synthetic acceptors detect with the expected
    # confidence — otherwise this test is asserting wrong state.
    assert sc._detect_entry_vector_role(weak_upd, gb_l0) == ("", "weak")
    assert sc._detect_entry_vector_role(strict_a1, gb_l0) == (
        "Alpha1", "strict",
    )
    # Submit weak first, strict second.
    entries = [
        _entry_from_record(weak_upd, "weak"),
        _entry_from_record(strict_a1, "strict"),
    ]
    sc._auto_bind_entry_vectors_from_entries(entries)
    # UPD goes to weak (only candidate); Alpha1 goes to strict.
    assert sc._get_entry_vector("gb_l0", "")["id"] == "weak"
    assert sc._get_entry_vector("gb_l0", "Alpha1")["id"] == "strict"


def test_auto_bind_empty_input_returns_empty_string():
    """No entries → no bindings → empty summary string."""
    assert sc._auto_bind_entry_vectors_from_entries([]) == ""


# ── Modal: smoke test ──────────────────────────────────────────────────────

@pytest.mark.asyncio
async def test_entry_vectors_modal_lists_all_gb_l0_roles():
    """Modal renders five rows for gb_l0 (UPD + α1/α2/Ω1/Ω2)."""
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.EntryVectorsModal("gb_l0")
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import DataTable
        t = modal.query_one("#ev-table", DataTable)
        # 5 rows = UPD + 4 named roles in `_CONSTRUCTOR_BACKBONES["gb_l0"]`.
        assert t.row_count == 5


@pytest.mark.asyncio
async def test_entry_vectors_modal_hint_uses_bound_vector_marker():
    """When a vector is bound to a role, the Hint column must show
    the antibiotic detected from the vector's gb_text — NOT the
    role's canonical `_CONSTRUCTOR_BACKBONES["selection"]` default.

    Regression (2026-05-22): a user binding an AmpR-bearing custom
    α-vector to Alpha1 still saw "Spectinomycin" in the hint column
    because `_refresh_table` rendered the convention default
    unchanged. Detection-from-bound now overrides it.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from io import StringIO
    from Bio import SeqIO
    vec_rec = SeqRecord(Seq("ATGC" * 25), id="amp_alpha",
                         name="amp_alpha")
    vec_rec.annotations["molecule_type"] = "DNA"
    vec_rec.annotations["topology"]      = "circular"
    vec_rec.features.append(SeqFeature(
        FeatureLocation(0, 30), type="CDS",
        qualifiers={"label": ["AmpR"]},
    ))
    buf = StringIO()
    SeqIO.write([vec_rec], buf, "genbank")
    sc._set_entry_vector(
        "gb_l0",
        {"name": "amp_alpha", "size": 100, "gb_text": buf.getvalue()},
        role="Alpha1",
    )
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.EntryVectorsModal("gb_l0")
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import DataTable
        t = modal.query_one("#ev-table", DataTable)
        # Find the Alpha1 row and inspect its Hint cell.
        alpha1_row = None
        for row_key in t.rows:
            if str(row_key.value) == "Alpha1":
                alpha1_row = row_key
                break
        assert alpha1_row is not None, "Alpha1 row missing"
        cells = t.get_row(alpha1_row)
        # Hint is the third column. The rendered Text object stringifies
        # to its plain content — substring-check for the detected marker.
        hint_text = str(cells[2])
        assert "Ampicillin" in hint_text, (
            f"expected 'Ampicillin' in hint, got {hint_text!r}"
        )
        assert "Spectinomycin" not in hint_text, (
            f"convention default leaked through: {hint_text!r}"
        )


@pytest.mark.asyncio
async def test_entry_vectors_modal_hint_omits_marker_when_unbound():
    """An unbound role's hint must NOT carry a hardcoded antibiotic.
    The 2026-05-22 overhaul removed `_CONSTRUCTOR_BACKBONES`'
    "selection" defaults; an unbound row shows only the slot
    descriptor so the user isn't misled into thinking the role
    enforces a particular marker.
    """
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.EntryVectorsModal("gb_l0")
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import DataTable
        t = modal.query_one("#ev-table", DataTable)
        alpha1_row = None
        for row_key in t.rows:
            if str(row_key.value) == "Alpha1":
                alpha1_row = row_key
                break
        assert alpha1_row is not None
        cells = t.get_row(alpha1_row)
        hint_text = str(cells[2])
        # No hardcoded antibiotic in the hint when unbound.
        for forbidden in ("Spectinomycin", "Kanamycin", "Ampicillin"):
            assert forbidden not in hint_text, (
                f"hardcoded {forbidden!r} leaked into unbound hint: "
                f"{hint_text!r}"
            )
        # The slot descriptor still surfaces so the user knows
        # which role this row represents.
        assert "slot" in hint_text.lower()


class TestMarkerWarnings:
    """`EntryVectorsModal._marker_warnings` enforces two rules from
    the Golden Braid iteration protocol:

      1. Intra-pair: slots within a family (α1/α2, Ω1/Ω2, …) must
         share an antibiotic.
      2. Cross-family: distinct families must use distinct
         antibiotics, so each iteration cycle's bench selection
         is unambiguous.

    Partial bindings (fewer than 2 bound + detected per family,
    or only one family bound) skip the corresponding check.
    """

    def test_no_warnings_when_consistent_alpha_amp_omega_kan(self):
        family_markers = {
            "Alpha": [("Alpha1", "Ampicillin"), ("Alpha2", "Ampicillin")],
            "Omega": [("Omega1", "Kanamycin"), ("Omega2", "Kanamycin")],
        }
        assert sc.EntryVectorsModal._marker_warnings(family_markers) == []

    def test_warns_on_alpha_pair_mismatch(self):
        family_markers = {
            "Alpha": [("Alpha1", "Ampicillin"), ("Alpha2", "Spectinomycin")],
        }
        warnings = sc.EntryVectorsModal._marker_warnings(family_markers)
        assert len(warnings) == 1
        assert "Alpha" in warnings[0]
        assert "Ampicillin" in warnings[0]
        assert "Spectinomycin" in warnings[0]
        assert "pair mismatch" in warnings[0].lower()

    def test_warns_on_omega_pair_mismatch(self):
        family_markers = {
            "Omega": [("Omega1", "Kanamycin"), ("Omega2", "Hygromycin")],
        }
        warnings = sc.EntryVectorsModal._marker_warnings(family_markers)
        assert len(warnings) == 1
        assert "Omega" in warnings[0]

    def test_warns_on_alpha_omega_collision(self):
        family_markers = {
            "Alpha": [("Alpha1", "Spectinomycin"), ("Alpha2", "Spectinomycin")],
            "Omega": [("Omega1", "Spectinomycin"), ("Omega2", "Spectinomycin")],
        }
        warnings = sc.EntryVectorsModal._marker_warnings(family_markers)
        # No pair mismatch, but α and Ω share Spectinomycin.
        assert len(warnings) == 1
        assert "Spectinomycin" in warnings[0]
        assert "Alpha" in warnings[0] and "Omega" in warnings[0]

    def test_warns_on_pair_mismatch_and_collision_together(self):
        family_markers = {
            "Alpha": [("Alpha1", "Ampicillin"), ("Alpha2", "Kanamycin")],
            "Omega": [("Omega1", "Kanamycin"), ("Omega2", "Kanamycin")],
        }
        warnings = sc.EntryVectorsModal._marker_warnings(family_markers)
        # 1 pair mismatch (Alpha) + 1 collision (Kanamycin shared
        # with Omega via Alpha2).
        assert len(warnings) == 2

    def test_no_warning_when_only_one_alpha_bound(self):
        """Single bound role in a family can't be checked for
        intra-pair mismatch — user is mid-config, not in error."""
        family_markers = {
            "Alpha": [("Alpha1", "Ampicillin")],
        }
        assert sc.EntryVectorsModal._marker_warnings(family_markers) == []

    def test_no_warning_when_only_one_family_bound(self):
        family_markers = {
            "Alpha": [("Alpha1", "Ampicillin"), ("Alpha2", "Ampicillin")],
        }
        assert sc.EntryVectorsModal._marker_warnings(family_markers) == []

    def test_empty_input_yields_no_warnings(self):
        assert sc.EntryVectorsModal._marker_warnings({}) == []


@pytest.mark.asyncio
async def test_entry_vectors_modal_status_shows_pair_mismatch():
    """End-to-end: binding two Alpha vectors with different markers
    surfaces a pair-mismatch warning in the modal's status line."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from io import StringIO
    from Bio import SeqIO

    def _gb_with(marker_label: str) -> str:
        r = SeqRecord(Seq("ATGC" * 25), id="v", name="v")
        r.annotations["molecule_type"] = "DNA"
        r.annotations["topology"]      = "circular"
        r.features.append(SeqFeature(
            FeatureLocation(0, 30), type="CDS",
            qualifiers={"label": [marker_label]},
        ))
        buf = StringIO()
        SeqIO.write([r], buf, "genbank")
        return buf.getvalue()

    sc._set_entry_vector(
        "gb_l0",
        {"name": "amp_a1", "size": 100, "gb_text": _gb_with("AmpR")},
        role="Alpha1",
    )
    sc._set_entry_vector(
        "gb_l0",
        {"name": "spec_a2", "size": 100, "gb_text": _gb_with("SpecR")},
        role="Alpha2",
    )
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.EntryVectorsModal("gb_l0")
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import Static
        status = modal.query_one("#ev-status", Static)
        text = str(status.render())
        assert "pair mismatch" in text.lower()


@pytest.mark.asyncio
async def test_entry_vectors_modal_blocks_undo():
    """Modal sets `_blocks_undo = True` so app-level Ctrl+Z doesn't
    fire underneath while the user mutates persistent state."""
    modal = sc.EntryVectorsModal("gb_l0")
    assert modal._blocks_undo is True


@pytest.mark.asyncio
async def test_entry_vectors_modal_default_focus_on_table():
    """Default focus is the role table so arrow keys + Enter drive
    the picks. Esc → close."""
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.EntryVectorsModal("gb_l0")
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import DataTable
        t = modal.query_one("#ev-table", DataTable)
        assert t.has_focus


def test_settings_modal_has_entry_vectors():
    """The Entry Vectors launcher must be reachable from the Settings
    surface. As of 2026-05-22 the menu was consolidated into a single
    `SettingsModal`; the Entry Vectors entry lives there as a button
    (`#set-entry-vectors`) plus the unchanged `action_open_entry_vectors`
    method on `PlasmidApp`. Greppable so future refactors trip the test."""
    with open(sc.__file__, encoding="utf-8") as f:
        src = f.read()
    assert "action_open_entry_vectors" in src
    assert "set-entry-vectors" in src
    assert "Entry Vectors" in src


# ── Modal auto-detect performance regression guards ────────────────────────

def test_auto_detect_runs_in_worker():
    """Regression guard for 2026-05-22 perf fix.

    `EntryVectorsModal._auto_detect_worker` must be decorated with
    `@work(thread=True, ...)` so a 200–500 plasmid active collection
    doesn't freeze the UI for seconds. Pre-fix the detection loop
    ran inline in `_auto_btn` on the UI thread."""
    worker = sc.EntryVectorsModal._auto_detect_worker
    # Textual's @work decorator preserves the original callable on
    # `_textual_worker_function` or wraps it; either way the wrapper
    # exposes `_thread` / `group` metadata. Sanity check: the method
    # is not the same as a plain function — it's wrapped.
    assert callable(worker)
    # Check by source: the decorator line is present immediately above
    # the def. (Whitebox is acceptable here per the project pattern:
    # `test_settings_menu_has_entry_vectors` uses the same approach.)
    with open(sc.__file__, encoding="utf-8") as f:
        src = f.read()
    assert (
        'group="ev-auto-detect"' in src
    ), "auto-detect worker must use exclusive group 'ev-auto-detect'"
    # Confirm the worker method is decorated with @work + thread=True
    # immediately above its def line.
    import re
    m = re.search(
        r'@work\(thread=True[^)]*group="ev-auto-detect"[^)]*\)\s+'
        r'def _auto_detect_worker',
        src,
    )
    assert m, "@work decorator missing on _auto_detect_worker"


def test_detect_cache_hit_skips_recompute(monkeypatch):
    """Regression guard for 2026-05-22 cache layer.

    Second call with the same (gb_text, grammar_id) must NOT call
    `_detect_entry_vector_role` again — that's the whole point of
    the cache (avoids 4–30 s re-digest cost on re-clicks)."""
    # Clear cache to start clean.
    sc._clear_entry_vector_detect_cache()
    grammar = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
    )
    from io import StringIO
    from Bio import SeqIO
    buf = StringIO()
    SeqIO.write([rec], buf, "genbank")
    gb_text = buf.getvalue()

    # First call: cold miss → runs detection.
    result1 = sc._detect_entry_vector_role_cached(gb_text, grammar)
    assert result1 is not None
    assert result1[0] == "Alpha1"

    # Patch the underlying detector to a tripwire — any further call
    # means the cache missed.
    called = []
    real = sc._detect_entry_vector_role

    def _tripwire(record, g):
        called.append((record, g))
        return real(record, g)

    monkeypatch.setattr(sc, "_detect_entry_vector_role", _tripwire)

    # Second call: cache hit → tripwire must NOT fire.
    result2 = sc._detect_entry_vector_role_cached(gb_text, grammar)
    assert result2 == result1
    assert called == [], (
        f"cache hit should have skipped recompute, but detector "
        f"was called {len(called)} time(s)"
    )


def test_detect_cache_invalidates_on_grammar_save():
    """Regression guard for 2026-05-22 cache layer.

    `_save_custom_grammars` must clear the EV detection cache so an
    enzyme change in a user grammar doesn't return a stale role. Same
    invalidation pattern as `_blast_clear_cache` from `_save_collections`
    ([PIT-16])."""
    sc._clear_entry_vector_detect_cache()
    grammar = sc._all_grammars()["gb_l0"]
    rec = _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
    )
    from io import StringIO
    from Bio import SeqIO
    buf = StringIO()
    SeqIO.write([rec], buf, "genbank")
    gb_text = buf.getvalue()

    # Warm the cache.
    sc._detect_entry_vector_role_cached(gb_text, grammar)
    assert len(sc._ENTRY_VECTOR_DETECT_CACHE) >= 1

    # Save grammars → cache should be empty.
    sc._save_custom_grammars([])
    assert len(sc._ENTRY_VECTOR_DETECT_CACHE) == 0, (
        "_save_custom_grammars must clear _ENTRY_VECTOR_DETECT_CACHE — "
        "stale results would otherwise survive grammar enzyme edits."
    )


def test_detect_cache_caches_parse_failures():
    """Regression guard for 2026-05-22 cache layer.

    A gb_text that fails to parse should cache `None` so a re-click
    doesn't re-pay the (failing) parse attempt. Without this, malformed
    entries in a 500-plasmid collection would re-parse on every
    Auto-detect click."""
    sc._clear_entry_vector_detect_cache()
    grammar = sc._all_grammars()["gb_l0"]
    junk = "not a genbank file at all"
    result1 = sc._detect_entry_vector_role_cached(junk, grammar)
    assert result1 is None
    assert len(sc._ENTRY_VECTOR_DETECT_CACHE) == 1, (
        "parse failure should be cached as None — otherwise re-click "
        "re-pays the failing parse"
    )


def test_auto_detect_uses_batched_save():
    """Regression guard for 2026-05-22 perf fix.

    The worker must call `_set_entry_vectors_batch` (one save) not
    N × `_set_entry_vector` (N saves). Each `_set_entry_vector`
    triggers `_safe_save_json` + `.bak` + fsync, so a 5-role grammar
    pre-fix did 5 disk round-trips. Now one."""
    with open(sc.__file__, encoding="utf-8") as f:
        src = f.read()
    # Find the worker body and assert the batch helper is called
    # inside it.
    worker_start = src.find("def _auto_detect_worker")
    assert worker_start > 0, "worker method not found"
    # Find the next method def after _auto_detect_worker.
    next_def = src.find("\n    def ", worker_start + 10)
    worker_body = src[worker_start:next_def]
    assert "_set_entry_vectors_batch" in worker_body, (
        "worker must call _set_entry_vectors_batch (one save), "
        "not per-role _set_entry_vector"
    )
    # And the per-role helper should NOT be inside the worker body —
    # if it slips back in, we're back to N saves.
    assert "_set_entry_vector(" not in worker_body, (
        "_set_entry_vector(per-role) leaked back into the worker; "
        "use _set_entry_vectors_batch for the bulk write"
    )


# ══════════════════════════════════════════════════════════════════════════
#  Field report 2026-11-08 — "is this actually an acceptor?", asked out loud
# ══════════════════════════════════════════════════════════════════════════
#
# Two α/Ω entry vectors in daily use carry FOUR sites of the assembly enzyme:
# a second cassette in the backbone brought its own pair, and one of the
# surplus fragments comes out with AATG/GCTT ends — the exact overhang pair an
# incoming L0 CDS part presents, so it competes for BOTH junctions of every
# assembly that vector is used in. Nothing said so, because nothing was asking
# the question: `_detect_entry_vector_role` returned None, the overhang lookup
# turned its own refusal into "no vector configured", and `_synthesis_lint`
# with the grammar's enzymes as `forbidden_enzymes` scored a WORKING acceptor
# 20/100 because an acceptor is SUPPOSED to carry those sites.


def _build_four_site_acceptor(part_oh5="AATG", part_oh3="GCTT"):
    """An α1-shaped acceptor carrying a SECOND inner-enzyme pair in its
    backbone, releasing a fragment whose ends are a real grammar part pair.
    The synthetic equivalent of the reported vectors."""
    bsa_fwd, bsa_rc = "GGTCTC" + "A", "T" + sc._rc("GGTCTC")
    extra = (bsa_fwd + part_oh5 + "CCTTGGAACC" * 30 + part_oh3 + bsa_rc)
    return _build_acceptor(
        oh5_inner="GGAG", oh3_inner="CGCT",
        oh5_outer="GGAG", oh3_outer="GTCA",
        inner_enzyme_site="GGTCTC", outer_enzyme_site="CGTCTC",
        backbone_body="ATCGATCGAT" * 60 + extra + "ATCGATCGAT" * 60,
    )



def _gb_text(rec) -> str:
    """GenBank text for a built record — what an entry-vector binding stores."""
    from io import StringIO
    from Bio import SeqIO
    buf = StringIO()
    rec.annotations.setdefault("molecule_type", "DNA")
    SeqIO.write(rec, buf, "genbank")
    return buf.getvalue()


class TestAcceptorVectorCheck:
    def _check(self, rec):
        return sc._check_acceptor_vector(
            str(rec.seq), circular=True,
            grammar=sc._all_grammars()["gb_l0"],
            features=sc._acceptor_vector_features(rec))

    def test_a_clean_alpha1_passes_on_both_enzymes(self):
        rec = _build_acceptor(
            oh5_inner="GGAG", oh3_inner="CGCT",
            oh5_outer="GGAG", oh3_outer="GTCA")
        res = self._check(rec)
        assert res["ok"] is True
        assert sorted(res["usable_enzymes"]) == ["BsaI", "Esp3I"]
        assert res["warnings"] == []

    def test_the_four_site_vector_is_refused_for_that_enzyme(self):
        res = self._check(_build_four_site_acceptor())
        row = next(r for r in res["enzymes"] if r["enzyme"] == "BsaI")
        assert row["n_sites"] == 4
        assert row["usable"] is False
        assert "exactly 2" in row["note"]

    def test_and_the_competing_fragment_is_named_with_its_overhangs(self):
        """The finding the endpoint exists for: not "digest failed" but WHICH
        fragment impersonates a part, and at which junction."""
        res = self._check(_build_four_site_acceptor())
        row = next(r for r in res["enzymes"] if r["enzyme"] == "BsaI")
        comp = [x for x in row["extra"] if x["collides_with"]]
        assert comp, row["extra"]
        assert any(set(x["collides_with"]) >= {"AATG", "GCTT"} for x in comp)
        assert any("COMPETES with a real part at BOTH junctions" in w
                   for w in res["warnings"]), res["warnings"]

    def test_the_other_enzyme_still_works_and_says_so(self):
        """A vector can be a perfect Esp3I acceptor and a hazard under BsaI;
        `ok` must not flatten that."""
        res = self._check(_build_four_site_acceptor())
        assert res["usable_enzymes"] == ["Esp3I"]
        assert res["ok"] is True

    def test_the_dropout_is_found_by_marker_exclusion_not_size(self):
        """Sacred — `feedback_never_assume_smaller_frag_is_payload`. The
        stuffer here is BIGGER than the backbone."""
        rec = _build_acceptor(
            oh5_inner="GGAG", oh3_inner="CGCT",
            oh5_outer="GGAG", oh3_outer="GTCA",
            backbone_body="ATCGATCGAT" * 70,
            stuffer_body="TTGGAACCAA" * 200)
        res = self._check(rec)
        row = next(r for r in res["enzymes"] if r["enzyme"] == "BsaI")
        assert row["usable"] is True
        assert row["dropout"]["has_backbone_marker"] is False
        assert row["dropout"]["length"] > max(
            f["length"] for f in row["fragments"]
            if f is not row["dropout"])

    def test_a_single_site_is_named_as_a_linearisation(self):
        seq = ("GGTCTCA" + "GGAG" + "TTGGAACCAA" * 40 + "ATCGATCGAT" * 40)
        res = sc._check_acceptor_vector(
            seq, circular=True, grammar=sc._all_grammars()["gb_l0"])
        row = next(r for r in res["enzymes"] if r["enzyme"] == "BsaI")
        assert row["n_sites"] == 1 and row["usable"] is False
        assert "linearises" in row["note"]

    def test_a_linear_record_is_refused_with_a_reason(self):
        res = sc._check_acceptor_vector(
            "ACGT" * 100, circular=False,
            grammar=sc._all_grammars()["gb_l0"])
        assert res["ok"] is False
        assert any("LINEAR" in e for e in res["errors"])

    def test_an_enzyme_that_does_not_cut_is_not_silently_fine(self):
        rec = _build_acceptor(
            oh5_inner="GGAG", oh3_inner="CGCT",
            oh5_outer="GGAG", oh3_outer="GTCA",
            outer_enzyme_site="GGTCTC")       # no Esp3I anywhere
        res = self._check(rec)
        row = next(r for r in res["enzymes"] if r["enzyme"] == "Esp3I")
        assert row["usable"] is False
        assert "does not cut" in row["note"]


class TestSynthesisLintAcceptorMode:
    def test_forbidden_mode_punishes_a_working_acceptor(self):
        """The premise — and it is CORRECT behaviour for a part, which is why
        the fix is a separate mode and not a changed rule."""
        rec = _build_acceptor(oh5_inner="GGAG", oh3_inner="CGCT",
                              oh5_outer="GGAG", oh3_outer="GTCA")
        rep = sc._synthesis_lint(str(rec.seq), circular=True,
                                 forbidden_enzymes=("BsaI", "Esp3I"))
        assert any(w["kind"] == "type_iis_site" for w in rep["warnings"])

    def test_acceptor_mode_expects_two_sites_and_deducts_nothing(self):
        rec = _build_acceptor(oh5_inner="GGAG", oh3_inner="CGCT",
                              oh5_outer="GGAG", oh3_outer="GTCA")
        rep = sc._synthesis_lint(str(rec.seq), circular=True,
                                 acceptor_enzymes=("BsaI", "Esp3I"))
        assert rep["stats"]["acceptor_sites"] == {"BsaI": 2, "Esp3I": 2}
        assert not [w for w in rep["warnings"]
                    if w["kind"] in ("acceptor_site", "type_iis_site")]

    def test_acceptor_mode_errors_on_four_sites(self):
        rec = _build_four_site_acceptor()
        rep = sc._synthesis_lint(str(rec.seq), circular=True,
                                 acceptor_enzymes=("BsaI", "Esp3I"))
        bad = [w for w in rep["warnings"] if w["kind"] == "acceptor_site"]
        assert len(bad) == 1 and bad[0]["level"] == "error"
        assert "exactly" in bad[0]["message"]
        assert rep["stats"]["acceptor_sites"]["BsaI"] == 4

    def test_an_acceptor_enzyme_is_not_also_scanned_as_forbidden(self):
        rec = _build_acceptor(oh5_inner="GGAG", oh3_inner="CGCT",
                              oh5_outer="GGAG", oh3_outer="GTCA")
        rep = sc._synthesis_lint(str(rec.seq), circular=True,
                                 forbidden_enzymes=("BsaI", "Esp3I"),
                                 acceptor_enzymes=("BsaI", "Esp3I"))
        assert not [w for w in rep["warnings"]
                    if w["kind"] == "type_iis_site"]

    def test_an_isoschizomer_of_an_acceptor_enzyme_is_excused_too(self):
        """Esp3I and BsmBI are the same CGTCTC site cut the same way. Excusing
        one by name while flagging the other described ONE physical site two
        contradictory ways."""
        rec = _build_acceptor(oh5_inner="GGAG", oh3_inner="CGCT",
                              oh5_outer="GGAG", oh3_outer="GTCA")
        assert sc._enzyme_signature("BsmBI") == sc._enzyme_signature("Esp3I")
        rep = sc._synthesis_lint(str(rec.seq), circular=True,
                                 forbidden_enzymes=("BsmBI",),
                                 acceptor_enzymes=("Esp3I",))
        assert not [w for w in rep["warnings"]
                    if w["kind"] == "type_iis_site"], rep["warnings"]


class TestWhichEnzymeReleasesWhichLevel:
    """The Golden-Braid convention, pinned against REAL vectors.

    `_enzyme_for_level_up(gb_l0, 0)` returns Esp3I — the DOMESTICATION enzyme,
    because that is what `gb_l0.enzyme` means — while a real pDGB3-convention
    α acceptor is assembled with BsaI and released with Esp3I. The code below
    `_resolve_acceptor_role` has always had it right and a comment above it
    once had it backwards; these assertions are against the plasmids, not
    against either one, so the next reader cannot "fix" the code to match the
    parity (field report 2026-11-08).
    """

    def test_alpha1_fixture_detects_as_alpha1(self):
        rec = sc.load_genbank("tests/FFE 2 ENTRY A1.dna")
        assert sc._detect_entry_vector_role(
            rec, sc._all_grammars()["gb_l0"]) == ("Alpha1", "strict")

    def test_omega1_fixture_detects_as_omega1(self):
        rec = sc.load_genbank("tests/FFE 4 ENTRY O1.dna")
        assert sc._detect_entry_vector_role(
            rec, sc._all_grammars()["gb_l0"]) == ("Omega1", "strict")

    def test_an_alpha_acceptor_takes_its_part_in_on_the_level_up_enzyme(self):
        """BsaI — `level_up_enzyme` — is what releases the canonical TU
        boundary stuffer from α, i.e. what an L0 part is assembled with."""
        gb = sc._all_grammars()["gb_l0"]
        rec = sc.load_genbank("tests/FFE 2 ENTRY A1.dna")
        cls = sc._classify_part_from_plasmid(str(rec.seq), circular=True)
        assert cls is not None
        assert (cls["position"]["oh5"], cls["position"]["oh3"]) == \
            sc._grammar_tu_overhangs(gb)
        assert cls["release_enzyme"] == gb["level_up_enzyme"] == "BsaI"

    def test_an_omega_acceptor_takes_its_tu_in_on_the_primary_enzyme(self):
        gb = sc._all_grammars()["gb_l0"]
        rec = sc.load_genbank("tests/FFE 4 ENTRY O1.dna")
        cls = sc._classify_part_from_plasmid(str(rec.seq), circular=True)
        assert cls is not None
        assert cls["release_enzyme"] == gb["enzyme"] == "Esp3I"

    def test_the_level_up_parity_disagrees_and_the_fallback_covers_it(self):
        """Not a bug report — a pinned fact. If someone "simplifies"
        `_assembly_fragment_enzyme_candidates` down to the parity's single
        enzyme, every L0→α assembly on a real GoldenBraid library breaks."""
        gb = sc._all_grammars()["gb_l0"]
        assert sc._enzyme_for_level_up(gb, 0) == "Esp3I"      # the parity
        cands = sc._assembly_fragment_enzyme_candidates(gb, "Esp3I")
        assert "BsaI" in cands, cands


class TestAcceptorQcInTheModal:
    """The warning has to reach the screen where the user binds the vector —
    the report's whole complaint is that nothing on any surface said it."""

    def test_a_clean_vector_produces_no_line(self):
        rec = _build_acceptor(oh5_inner="GGAG", oh3_inner="CGCT",
                              oh5_outer="GGAG", oh3_outer="GTCA")
        gb_text = _gb_text(rec)
        assert sc.EntryVectorsModal._acceptor_qc_lines(
            sc._all_grammars()["gb_l0"], "Alpha1", gb_text) == []

    def test_the_four_site_vector_produces_one_compact_line(self):
        gb_text = _gb_text(_build_four_site_acceptor())
        lines = sc.EntryVectorsModal._acceptor_qc_lines(
            sc._all_grammars()["gb_l0"], "Alpha1", gb_text)
        assert len(lines) == 1, lines
        line = lines[0]
        assert line.startswith("Alpha1: 4 BsaI sites")
        assert "competes with a real part" in line
        assert "AATG/GCTT" in line

    def test_an_oversized_vector_is_skipped_not_digested_on_the_ui_thread(self):
        """`_refresh_table` runs this once per bound vector on every refresh.
        A chromosome bound by mistake must not freeze the modal."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        big = SeqRecord(Seq("ATCG" * 60_000), id="big", name="big",
                        annotations={"topology": "circular",
                                     "molecule_type": "DNA"})
        assert sc.EntryVectorsModal._acceptor_qc_lines(
            sc._all_grammars()["gb_l0"], "Alpha1", _gb_text(big)) == []

    def test_unparseable_gb_text_is_not_fatal(self):
        """A QC that cannot run must never stop the user configuring their
        vectors."""
        assert sc.EntryVectorsModal._acceptor_qc_lines(
            sc._all_grammars()["gb_l0"], "Alpha1", "not genbank at all") == []
        assert sc.EntryVectorsModal._acceptor_qc_lines(
            sc._all_grammars()["gb_l0"], "Alpha1", "") == []


@pytest.mark.asyncio
async def test_entry_vectors_modal_flags_a_bound_four_site_vector():
    gb_text = _gb_text(_build_four_site_acceptor())
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        sc._set_entry_vector("gb_l0", {
            "name": "four-site α1", "size": len(gb_text),
            "source": "test", "id": "foursite", "gb_text": gb_text},
            "Alpha1")
        modal = sc.EntryVectorsModal("gb_l0")
        app.push_screen(modal)
        for _ in range(4):
            await pilot.pause()
        from textual.widgets import Static, DataTable
        status = modal.query_one("#ev-status", Static)
        # Textual keeps a Static's content on different attributes across
        # versions; read whichever holds it rather than pinning one.
        text = str(getattr(status, "_content", "")
                   or getattr(status, "_renderable", "")
                   or status.render())
        assert "4 BsaI sites" in text, text
        assert "competes with a real part" in text
        # and the row itself is marked, so the table alone is not reassuring
        t = modal.query_one("#ev-table", DataTable)
        row = t.get_row_at(1)
        assert "⚠" in str(row[1])


class TestAcceptorCheckIsBounded:
    """A 60-site plasmid answered with a 16 KB response, almost all of it
    fragments nobody will read."""

    SEQ = "GGTCTCAGGAGCCTTAA" * 60

    def _row(self):
        res = sc._check_acceptor_vector(
            self.SEQ, circular=True, grammar=sc._all_grammars()["gb_l0"])
        return next(r for r in res["enzymes"] if r["enzyme"] == "BsaI")

    def test_the_lists_are_capped(self):
        row = self._row()
        cap = sc._ACCEPTOR_MAX_FRAGMENTS
        assert len(row["fragments"]) == cap
        assert len(row["extra"]) == cap
        assert len(row["cut_positions"]) == cap
        assert row["cut_positions_truncated"] is True

    def test_but_the_counts_are_exact(self):
        """A truncated answer must not read as a complete one."""
        row = self._row()
        assert row["n_sites"] == 60
        assert row["n_fragments"] == 60
        assert row["n_extra"] == 60

    def test_the_response_stays_small(self):
        import json
        res = sc._check_acceptor_vector(
            self.SEQ, circular=True, grammar=sc._all_grammars()["gb_l0"])
        assert len(json.dumps(res, default=str)) < 12_000


def test_check_acceptor_vector_refuses_a_chromosome_scale_record():
    """An acceptor is a plasmid; a chromosome would digest into thousands of
    fragments nobody can act on."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq("ATCG" * 60_000), id="big", name="big",
                    annotations={"topology": "circular",
                                 "molecule_type": "DNA"})
    app = type("A", (), {"_current_record": rec})()
    body, code = sc._h_check_acceptor_vector(app, {})
    assert code == 413 and "capped at" in body["error"]
