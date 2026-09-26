"""
test_gibson_e4 — E4: Gibson from-plasmid linearization + homology-arm designer.

Linearization: a circular plasmid source opens at a user-chosen bp
("Linearize at" → `_rotate_seq_record`), so its Gibson fragment ends sit
where the user wants (not always base 0). Features straddling the cut are
re-framed wrap-correctly by the shared rotation helper.

Homology arms: "Design overlaps" (`_design_homology_arms`) appends a 5'
arm = the upstream fragment's 3'-terminal `min_overlap` bases to each
downstream fragment, so every junction reaches `min_overlap` and the
product is the fragments joined seamlessly (the arm collapses back into
the upstream's own 3' end). Idempotent — re-run after a reorder.

The last test verifies the whole flow end-to-end on a real circular
plasmid (`_make_demo_record`): linearize → add insert → design overlaps →
the simulator assembles a correct circular product.
"""
from __future__ import annotations

import pytest
from textual.widgets import DataTable, Input, RadioButton

import splicecraft as sc
from tests.test_smoke import _build_app, TERMINAL_SIZE


async def _pane(app, pilot):
    modal = sc.ConstructorModal()
    app.push_screen(modal)
    await pilot.pause()
    await pilot.pause(0.05)
    return modal.query_one("#ctor-gib-pane", sc.GibsonAssemblyPane)


def _frag(name, seq, features=None):
    return {"name": name, "sequence": seq,
            "features": list(features or []), "source": f"test:{name}"}


def _save_one(entry_id, name, rec):
    sc._save_library([{
        "id": entry_id, "name": name,
        "size": len(rec.seq), "gb_text": sc._record_to_gb_text(rec),
    }])


class TestGibsonLinearization:

    @pytest.mark.asyncio
    async def test_plasmid_source_linearizes_at_position(
            self, tiny_record, isolated_library):
        _save_one("circ1", "CircPlasmid", tiny_record)
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await _pane(app, pilot)
            pane._populate_library_table()
            await pilot.pause()
            pane.query_one("#gib-source-table", DataTable).move_cursor(row=0)
            pane.query_one("#gib-linearize-at", Input).value = "10"
            pane._on_add(None)
            await pilot.pause()
            assert len(pane._lane) == 1
            full = str(tiny_record.seq).upper()
            assert pane._lane[0]["sequence"] == full[10:] + full[:10]
            assert "@10" in pane._lane[0]["source"]

    @pytest.mark.asyncio
    async def test_linearize_zero_is_unrotated(
            self, tiny_record, isolated_library):
        _save_one("c", "C", tiny_record)
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await _pane(app, pilot)
            pane._populate_library_table()
            await pilot.pause()
            pane.query_one("#gib-source-table", DataTable).move_cursor(row=0)
            pane.query_one("#gib-linearize-at", Input).value = "0"
            pane._on_add(None)
            await pilot.pause()
            assert pane._lane[0]["sequence"] == str(tiny_record.seq).upper()

    @pytest.mark.asyncio
    async def test_linearize_past_end_rejected(
            self, tiny_record, isolated_library):
        n = len(tiny_record.seq)
        _save_one("c", "C", tiny_record)
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await _pane(app, pilot)
            pane._populate_library_table()
            await pilot.pause()
            pane.query_one("#gib-source-table", DataTable).move_cursor(row=0)
            pane.query_one("#gib-linearize-at", Input).value = str(n + 5)
            pane._on_add(None)   # ValueError caught in _on_add → notify + no-op
            await pilot.pause()
            assert len(pane._lane) == 0

    @pytest.mark.asyncio
    async def test_linear_source_ignores_linearize(
            self, tiny_record, isolated_library):
        # A linear source has fixed ends — 'Linearize at' is a no-op (used
        # as-is, no rotation, no @N tag) and must not crash.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        lin_rec = SeqRecord(Seq("ACGTACGTACGTACGTACGTAAAA"),
                            id="LinFrag", name="LinFrag")
        lin_rec.annotations["molecule_type"] = "DNA"
        lin_rec.annotations["topology"] = "linear"
        _save_one("linf", "LinFrag", lin_rec)
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await _pane(app, pilot)
            pane._populate_library_table()
            await pilot.pause()
            pane.query_one("#gib-source-table", DataTable).move_cursor(row=0)
            pane.query_one("#gib-linearize-at", Input).value = "5"
            pane._on_add(None)
            await pilot.pause()
            assert len(pane._lane) == 1
            assert pane._lane[0]["sequence"] == "ACGTACGTACGTACGTACGTAAAA"
            assert "@" not in pane._lane[0]["source"]


class TestGibsonArmDesigner:

    async def _mk(self, app, pilot, min_oh="15", circular=False):
        pane = await _pane(app, pilot)
        # Select the target radio (a RadioSet won't let its only pressed
        # button be turned OFF, so we press the OTHER one and let the set
        # settle exclusivity on the next pump cycle).
        target = "gib-topo-circular" if circular else "gib-topo-linear"
        pane.query_one(f"#{target}", RadioButton).value = True
        pane.query_one("#gib-min-overlap", Input).value = min_oh
        await pilot.pause()
        await pilot.pause()
        assert pane._is_circular() is circular   # topology actually applied
        return pane

    @pytest.mark.asyncio
    async def test_arms_make_linear_junction_assemble(
            self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await self._mk(app, pilot)
            a = "GGATCCACAGTGACTTGGCATATCAGTCGAC"         # 31 bp
            b = "CCTTAAGGACATCATTGCCTAGGAATTCACG"          # 30 bp, no overlap
            pane._lane = [_frag("A", a), _frag("B", b)]
            armed, already, skipped, _long = pane._design_homology_arms()
            assert armed == 1 and not skipped
            r = sc._simulate_gibson_assembly(pane._lane, min_overlap=15,
                                              circular=False)
            assert r["success"] is True, r
            b_armed = pane._lane[1]["sequence"]
            ov = sc._gibson_overlap_len(a, b_armed, min_overlap=15)
            assert ov >= 15
            # Seamless: product = upstream + (downstream minus the overlap).
            assert r["product_seq"].upper() == (a + b_armed[ov:]).upper()

    @pytest.mark.asyncio
    async def test_three_fragment_chain_assembles(
            self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await self._mk(app, pilot)
            a = "GGATCCACAGTGACTTGGCATATCAGTCGAC"
            b = "TTAACCGGACATGCATTGCCTAGGAATTCAA"
            c = "CACACAGTGTGTAAAACCCCGGGGTTTTACG"
            pane._lane = [_frag("A", a), _frag("B", b), _frag("C", c)]
            armed, already, skipped, _long = pane._design_homology_arms()
            assert armed == 2 and not skipped   # 2 junctions in a 3-frag chain
            r = sc._simulate_gibson_assembly(pane._lane, min_overlap=15,
                                              circular=False)
            assert r["success"] is True, r

    @pytest.mark.asyncio
    async def test_design_is_idempotent(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await self._mk(app, pilot)
            pane._lane = [_frag("A", "ACGT" * 12), _frag("B", "TTGG" * 12)]
            pane._design_homology_arms()
            after1 = [f["sequence"] for f in pane._lane]
            pane._design_homology_arms()
            after2 = [f["sequence"] for f in pane._lane]
            assert after1 == after2   # no double-stacking of arms

    @pytest.mark.asyncio
    async def test_features_shift_with_arm(
            self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await self._mk(app, pilot)
            feat = {"start": 0, "end": 10, "strand": 1,
                    "type": "CDS", "label": "g"}
            pane._lane = [_frag("A", "ACGT" * 12), _frag("B", "TTGG" * 12,
                                                          [feat])]
            pane._design_homology_arms()
            assert pane._lane[1]["_gib_arm5_len"] == 15
            assert pane._lane[1]["features"][0]["start"] == 15
            assert pane._lane[1]["features"][0]["end"] == 25

    @pytest.mark.asyncio
    async def test_single_fragment_is_noop(
            self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await self._mk(app, pilot)
            pane._lane = [_frag("A", "ACGT" * 12)]
            assert pane._design_homology_arms() is None

    @pytest.mark.asyncio
    async def test_upstream_too_short_skipped_not_crash(
            self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await self._mk(app, pilot)
            pane._lane = [_frag("short", "ACGTACGT"), _frag("B", "TTGG" * 12)]
            armed, already, skipped, _long = pane._design_homology_arms()
            assert armed == 0 and skipped == ["B"]


class TestGibsonLayout:

    @pytest.mark.asyncio
    async def test_linearize_row_has_bounded_height(
            self, tiny_record, isolated_library):
        # Regression: the E4 "Linearize at" row shipped (v1.2.24) without a CSS
        # height rule, so it ballooned to fill the source box (~23 rows) and
        # collapsed the assembly lane + action buttons below the fold. It must
        # stay a compact fixed row so the pane keeps a usable layout.
        from textual.widgets import TabbedContent
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            modal = sc.ConstructorModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            modal.query_one(TabbedContent).active = "ctor-tab-gibson"
            await pilot.pause()
            await pilot.pause(0.05)
            pane = modal.query_one("#ctor-gib-pane", sc.GibsonAssemblyPane)
            row = pane.query_one("#gib-linearize-row")
            h = row.region.height
            assert 0 < h <= 4, f"linearize row height = {h} (expected ~3)"
            # The assembly lane must be laid out with real height (not
            # squeezed to nothing by a runaway sibling).
            lane = pane.query_one("#gib-lane-host")
            assert lane.region.height >= 6, \
                f"lane host height = {lane.region.height} (too small)"


class TestGibsonE4RealPlasmid:

    @pytest.mark.asyncio
    async def test_real_plasmid_linearize_insert_assemble(
            self, isolated_library):
        demo = sc._make_demo_record()
        plen = len(demo.seq)
        _save_one("demo", "DemoPlasmid", demo)
        app = sc.PlasmidApp()
        app._preload_record = demo
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pane = await _pane(app, pilot)
            pane._populate_library_table()
            await pilot.pause()
            # Linearize the backbone at bp 100.
            pane.query_one("#gib-source-table", DataTable).move_cursor(row=0)
            pane.query_one("#gib-linearize-at", Input).value = "100"
            pane._on_add(None)
            await pilot.pause()
            assert len(pane._lane) == 1
            assert len(pane._lane[0]["sequence"]) == plen
            # Add a benign insert + design overlaps (circular).
            insert = "ATG" + "GCTAGCACC" * 8 + "TAA"
            pane._lane.append(_frag("insert", insert))
            pane.query_one("#gib-topo-circular", RadioButton).value = True
            pane.query_one("#gib-min-overlap", Input).value = "20"
            armed, already, skipped, _long = pane._design_homology_arms()
            assert not skipped
            r = sc._simulate_gibson_assembly(pane._lane, min_overlap=20,
                                              circular=True)
            assert r["success"] is True, r
            # Seamless circular product: backbone + insert, each shared arm
            # appearing exactly once → plen + len(insert).
            assert len(r["product_seq"]) == plen + len(insert)


class TestGibsonArmHardening:
    """Audit fixes 2026-07-16: F4 (arm design survives non-int feature coords)
    and F1 (a fragment consumed by its arms is reported invalid, not success)."""

    async def _mk(self, app, pilot, min_oh="15", circular=False):
        pane = await _pane(app, pilot)
        target = "gib-topo-circular" if circular else "gib-topo-linear"
        pane.query_one(f"#{target}", RadioButton).value = True
        pane.query_one("#gib-min-overlap", Input).value = min_oh
        await pilot.pause()
        await pilot.pause()
        return pane

    @pytest.mark.asyncio
    async def test_arm_design_survives_non_int_feature_coord(
            self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            pane = await self._mk(app, pilot, min_oh="15", circular=False)
            bad = _frag("A", "ACGTACGT" * 10,
                        [{"start": "x", "end": "y", "label": "bad"}])
            pane._lane = [bad, _frag("B", "TTTTGGGG" * 10)]
            armed, already, skipped, _long = pane._design_homology_arms()   # must not raise
            assert isinstance(armed, int)

    @pytest.mark.asyncio
    async def test_consumed_fragment_reported_invalid(
            self, tiny_record, isolated_library, monkeypatch):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            pane = await self._mk(app, pilot, min_oh="20", circular=True)
            # 2-fragment circular: a 20-bp insert (= min_overlap) is armed on
            # BOTH ends (leading + wrap) and fully consumed → assembly invalid.
            pane._lane = [_frag("BB", "ACGTACGT" * 10),          # 80 bp
                          _frag("INS", "GGGGTTTTAAAACCCCGGGG")]  # 20 bp
            msgs: list = []
            monkeypatch.setattr(pane.app, "notify",
                                lambda *a, **k: msgs.append((a, k)))
            pane._on_design_arms(None)
            for _ in range(3):
                await pilot.pause()
            # the assembly is genuinely invalid …
            r = sc._simulate_gibson_assembly(pane._lane, min_overlap=20,
                                             circular=True)
            assert not (r.get("success") and not r.get("errors"))
            # … and the notify says so instead of a bare success
            joined = " ".join(str(x) for m in msgs for x in m[0])
            assert "still invalid" in joined, msgs
