"""test_mapimage — circular plasmid-map image export (SVG + PNG).

Covers the pure renderer `splicecraft_mapimage`:
  * colour conversion (`_to_rgb`) across every accepted token form,
  * geometry (angle convention, wrap-aware arc membership, lane stacking),
  * SVG well-formedness + transparent-background behaviour,
  * PNG decodes at the right size/mode, transparent → alpha-0 corner,
  * the record extractors + high-level `export_plasmid_map` file write,
  * graceful handling of empty / wrap-crossing / label-dense plasmids.

The renderer is app-free; we import `splicecraft as sc` only for the
`_protect_user_data` sandbox, the demo record, and the live restriction
catalog. `_FALLBACK_PALETTE` is asserted to mirror the widget palette so
off-screen exports colour features exactly like the on-screen map.
"""
from __future__ import annotations

import math
from io import BytesIO
from xml.etree import ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from PIL import Image

import splicecraft as sc
import splicecraft_mapimage as mi
import splicecraft_widgets as _w
from tests._pilot_helpers import wait_for_modal, wait_for_state, wait_for_widget


def _rec(seq: str, feats=None):
    r = SeqRecord(Seq(seq), name="pTest", id="pTest")
    for (s, e, strand, ftype, label, color) in (feats or []):
        q = {"label": [label]}
        if color:
            q["ApEinfo_fwdcolor"] = [color]
        r.features.append(
            SeqFeature(FeatureLocation(s, e, strand=strand), type=ftype,
                       qualifiers=q)
        )
    return r


def _entry(name, *, gb_text=None):
    """A minimal library-entry dict (id/name/gb_text) for bulk-export tests."""
    if gb_text is None:
        gb_text = sc._record_to_gb_text(sc._make_demo_record())
    return {"id": name, "name": name, "gb_text": gb_text}


# ── Palette drift guard ──────────────────────────────────────────────────────

def test_palette_matches_widgets():
    # Off-screen exports must colour features identically to the live map.
    assert mi._FALLBACK_PALETTE == _w._FEATURE_PALETTE


# ── Colour conversion ────────────────────────────────────────────────────────

class TestToRgb:
    def test_hex6(self):
        assert mi._to_rgb("#FF8000") == (255, 128, 0)

    def test_hex3_expands(self):
        assert mi._to_rgb("#0a0") == (0, 170, 0)

    def test_xterm_color_paren(self):
        rgb = mi._to_rgb("color(196)")   # bright red in the cube
        assert rgb == (255, 0, 0)

    def test_bare_xterm_index(self):
        assert mi._to_rgb("46") == mi._xterm_to_rgb(46)

    def test_named(self):
        assert mi._to_rgb("blue") == (0, 0, 255)

    def test_rgb_func(self):
        assert mi._to_rgb("rgb(10,20,30)") == (10, 20, 30)

    def test_garbage_is_grey(self):
        assert mi._to_rgb("not-a-color") == (150, 150, 150)
        assert mi._to_rgb(None) == (150, 150, 150)
        assert mi._to_rgb("#zzzzzz") == (150, 150, 150)

    def test_every_palette_entry_converts(self):
        # No palette colour may fall through to the grey error sentinel.
        for c in mi._FALLBACK_PALETTE:
            assert mi._to_rgb(c) != (150, 150, 150)

    def test_ensure_readable_darkens_light(self):
        light = (255, 255, 200)
        assert mi._luminance(mi._ensure_readable(light)) < mi._luminance(light)
        dark = (30, 30, 30)
        assert mi._ensure_readable(dark) == dark


# ── Geometry ─────────────────────────────────────────────────────────────────

class TestGeometry:
    def test_origin_at_top(self):
        assert mi._angle(0, 1000, 0) == pytest.approx(-math.pi / 2)
        # Origin offset moves the reference bp to the top.
        assert mi._angle(250, 1000, 250) == pytest.approx(-math.pi / 2)

    def test_clockwise_quarter_turn(self):
        # A quarter of the plasmid ahead of the origin → right of centre (0 rad).
        assert mi._angle(250, 1000, 0) == pytest.approx(0.0, abs=1e-9)
        # Half-way → bottom (+pi/2).
        assert mi._angle(500, 1000, 0) == pytest.approx(math.pi / 2)

    def test_in_arc_no_wrap(self):
        assert mi._in_arc(50, 10, 100, 1000)
        assert not mi._in_arc(200, 10, 100, 1000)

    def test_in_arc_wraps_origin(self):
        # Arc [950, 1050) wraps to [950,1000)+[0,50).
        assert mi._in_arc(20, 950, 100, 1000)
        assert mi._in_arc(970, 950, 100, 1000)
        assert not mi._in_arc(500, 950, 100, 1000)

    def test_lanes_disjoint_share_lane0(self):
        feats = [{"start": 0, "end": 100}, {"start": 500, "end": 600}]
        lanes = mi._assign_lanes(feats, 1000)
        assert set(lanes.values()) == {0}

    def test_lanes_overlap_stack(self):
        feats = [{"start": 0, "end": 400}, {"start": 100, "end": 300}]
        lanes = mi._assign_lanes(feats, 1000)
        assert lanes[0] != lanes[1]

    def test_lanes_wrap_overlap_stack(self):
        # One arc crosses the origin and overlaps a feature at bp 10.
        feats = [{"start": 950, "end": 50}, {"start": 0, "end": 100}]
        lanes = mi._assign_lanes(feats, 1000)
        assert lanes[0] != lanes[1]

    def test_lanes_are_capped(self):
        # A huge overlapping pile-up must not spawn unbounded lanes (O(F²)).
        feats = [{"start": i, "end": i + 900} for i in range(200)]
        lanes = mi._assign_lanes(feats, 1000)
        assert max(lanes.values()) < mi._MAP_MAX_LANES

    def test_band_polygon_arrowhead_counts(self):
        pts_f = mi._band_polygon(100, 100, 30, 6, 0.0, 1.0, 1)
        pts_r = mi._band_polygon(100, 100, 30, 6, 0.0, 1.0, -1)
        pts_0 = mi._band_polygon(100, 100, 30, 6, 0.0, 1.0, 0)
        # All produce a closed, non-degenerate polygon.
        for pts in (pts_f, pts_r, pts_0):
            assert len(pts) >= 4
            assert all(len(p) == 2 for p in pts)


# ── SVG back-end ─────────────────────────────────────────────────────────────

class TestSvg:
    def _demo_svg(self, **kw):
        feats, total = mi._map_feats_from_record(sc._make_demo_record())
        return mi.render_plasmid_map_svg(feats, total, title="pDemo",
                                         size=800, **kw)

    def test_wellformed(self):
        root = ET.fromstring(self._demo_svg())
        assert root.tag.endswith("svg")
        assert root.attrib["width"] == "800"
        kinds = {c.tag.split("}")[-1] for c in root}
        assert "circle" in kinds       # backbone
        assert "polygon" in kinds      # feature arcs
        assert "text" in kinds         # labels / title

    def test_opaque_has_bg_rect(self):
        svg = self._demo_svg(transparent=False)
        assert "<rect" in svg and mi._BG_COLOR in svg

    def test_transparent_has_no_bg_rect(self):
        svg = self._demo_svg(transparent=True)
        assert "<rect" not in svg

    def test_title_and_bp_present(self):
        svg = self._demo_svg()
        assert "pDemo" in svg
        assert "bp" in svg

    def test_label_markup_escaped(self):
        rec = _rec("ACGT" * 200, [(10, 200, 1, "CDS", "a<b>&c", "#3366CC")])
        feats, total = mi._map_feats_from_record(rec)
        svg = mi.render_plasmid_map_svg(feats, total, title="x", size=600)
        ET.fromstring(svg)                       # must still parse
        assert "a<b>&c" not in svg               # raw angle brackets escaped
        assert "&lt;b&gt;" in svg


# ── PNG back-end ─────────────────────────────────────────────────────────────

class TestPng:
    def _demo_png(self, **kw):
        feats, total = mi._map_feats_from_record(sc._make_demo_record())
        return mi.render_plasmid_map_png(feats, total, title="pDemo",
                                         size=600, **kw)

    def test_decodes_at_size(self):
        img = Image.open(BytesIO(self._demo_png()))
        assert img.size == (600, 600)

    def test_opaque_is_rgb(self):
        img = Image.open(BytesIO(self._demo_png(transparent=False)))
        assert img.mode == "RGB"

    def test_transparent_corner_alpha_zero(self):
        img = Image.open(BytesIO(self._demo_png(transparent=True)))
        assert img.mode == "RGBA"
        assert img.getpixel((0, 0))[3] == 0     # top-left is fully transparent

    def test_supersample_1_still_valid(self):
        img = Image.open(BytesIO(self._demo_png(supersample=1)))
        assert img.size == (600, 600)


# ── Extractors ───────────────────────────────────────────────────────────────

class TestExtractors:
    def test_feats_shape(self):
        feats, total = mi._map_feats_from_record(sc._make_demo_record())
        assert total > 0
        assert feats
        for f in feats:
            assert {"start", "end", "strand", "color", "label", "type"} <= set(f)

    def test_source_feature_skipped(self):
        rec = _rec("ACGT" * 100)
        rec.features.append(SeqFeature(FeatureLocation(0, 400), type="source"))
        rec.features.append(
            SeqFeature(FeatureLocation(10, 50, strand=1), type="CDS",
                       qualifiers={"label": ["x"]}))
        feats, _ = mi._map_feats_from_record(rec)
        assert all(f["type"] != "source" for f in feats)
        assert len(feats) == 1

    def test_sites_from_record_recut_only(self):
        # pUC19-like: use the demo record + the live catalog.
        sites = mi._map_sites_from_record(sc._make_demo_record())
        assert isinstance(sites, list)
        assert all(s.get("type") == "recut" for s in sites)

    def test_qualifier_color_honoured(self):
        rec = _rec("ACGT" * 100, [(10, 50, 1, "CDS", "gene", "#AB12CD")])
        feats, _ = mi._map_feats_from_record(rec)
        assert feats[0]["color"] == "#AB12CD"


# ── Robustness ───────────────────────────────────────────────────────────────

class TestRobustness:
    def test_empty_record_no_crash(self):
        rec = SeqRecord(Seq(""), name="empty")
        feats, total = mi._map_feats_from_record(rec)
        assert total == 0
        svg = mi.render_plasmid_map_svg(feats, total, title="empty", size=400)
        ET.fromstring(svg)
        png = mi.render_plasmid_map_png(feats, total, title="empty", size=400)
        assert Image.open(BytesIO(png)).size == (400, 400)

    def test_no_feature_backbone_only(self):
        rec = _rec("ACGT" * 250)
        feats, total = mi._map_feats_from_record(rec)
        assert feats == []
        svg = mi.render_plasmid_map_svg(feats, total, title="bare", size=400)
        assert "circle" in svg              # backbone still drawn

    def test_wrap_crossing_feature(self):
        # Biopython models an origin-crossing feature as a join() of two
        # parts: [950:1000] + [0:50]. The extractor must yield end < start.
        rec = SeqRecord(Seq("ACGT" * 250), name="pWrap", id="pWrap")  # 1000 bp
        loc = CompoundLocation([FeatureLocation(950, 1000, strand=1),
                                FeatureLocation(0, 50, strand=1)])
        rec.features.append(
            SeqFeature(loc, type="CDS",
                       qualifiers={"label": ["wrapper"],
                                   "ApEinfo_fwdcolor": ["#20C020"]}))
        feats, total = mi._map_feats_from_record(rec)
        assert feats[0]["end"] < feats[0]["start"]      # origin-crossing arc
        svg = mi.render_plasmid_map_svg(feats, total, size=400)
        ET.fromstring(svg)

    def test_size_clamped(self):
        assert mi._clamp_size(10) == mi._MAP_IMAGE_MIN_SIZE
        assert mi._clamp_size(99999) == mi._MAP_IMAGE_MAX_SIZE
        assert mi._clamp_size("bad") == mi._MAP_IMAGE_DEFAULT_SIZE

    def test_label_dense_no_crash(self):
        # Many overlapping labelled features exercise the cap + de-collision.
        feats = [{"start": i * 3, "end": i * 3 + 200, "strand": 1,
                  "color": "#3399FF", "label": f"feat{i}", "type": "CDS"}
                 for i in range(120)]
        svg = mi.render_plasmid_map_svg(feats, 1000, title="dense", size=800)
        ET.fromstring(svg)


# ── High-level file export ───────────────────────────────────────────────────

class TestExport:
    def test_png_file_written(self, tmp_path):
        out = tmp_path / "map.png"
        summary = mi.export_plasmid_map(out, record=sc._make_demo_record(),
                                        fmt="png", title="pDemo", size=500)
        assert out.exists()
        assert summary["fmt"] == "png"
        assert summary["bytes"] == out.stat().st_size
        assert Image.open(out).size == (500, 500)

    def test_svg_file_written(self, tmp_path):
        out = tmp_path / "map.svg"
        mi.export_plasmid_map(out, record=sc._make_demo_record(), fmt="svg",
                              title="pDemo", size=500, transparent=True)
        assert out.exists()
        ET.fromstring(out.read_text())

    def test_feats_total_path(self, tmp_path):
        # Live-map path: pass parsed feats + total directly.
        feats, total = mi._map_feats_from_record(sc._make_demo_record())
        out = tmp_path / "m.png"
        summary = mi.export_plasmid_map(out, feats=feats, total=total,
                                        fmt="png", size=400)
        assert summary["features"] == len(feats)
        assert out.exists()

    def test_unknown_format_rejected(self, tmp_path):
        with pytest.raises(ValueError):
            mi.export_plasmid_map(tmp_path / "x.bmp", record=sc._make_demo_record(),
                                  fmt="bmp")

    def test_missing_inputs_rejected(self, tmp_path):
        with pytest.raises(ValueError):
            mi.export_plasmid_map(tmp_path / "x.png", fmt="png")

    def test_atomic_write_no_tempfile_left(self, tmp_path):
        out = tmp_path / "map.png"
        mi.export_plasmid_map(out, record=sc._make_demo_record(), fmt="png",
                              size=400)
        leftovers = [p.name for p in tmp_path.iterdir()
                     if p.name.startswith(".sc-map-")]
        assert leftovers == []


# ── UI integration (modal + hotkey action) ──────────────────────────────────

class TestModalIntegration:
    """Drive the real MapImageExportModal + action through a headless app so
    the hotkey → modal → worker → renderer → file path is exercised end-to-end."""

    async def test_modal_exports_png(self, tmp_path):
        from textual.widgets import Input
        app = sc.PlasmidApp()
        async with app.run_test(size=(170, 50)) as pilot:
            await pilot.pause()
            feats, total = mi._map_feats_from_record(sc._make_demo_record())
            modal = sc.MapImageExportModal(
                feats=feats, total=total, title="pDemo",
                default_path=str(tmp_path / "m.png"),
                record=sc._make_demo_record())
            app.push_screen(modal)
            await wait_for_widget(pilot, "#mapimg-filename")
            modal._selected_dir = str(tmp_path)
            modal.query_one("#mapimg-filename", Input).value = "out.png"
            modal._do_export()
            out = tmp_path / "out.png"
            await wait_for_state(pilot, lambda: out.exists(), what="png written")
            assert Image.open(out).size == (1400, 1400)

    async def test_modal_exports_svg_transparent(self, tmp_path):
        from textual.widgets import Input, Select, Checkbox
        app = sc.PlasmidApp()
        async with app.run_test(size=(170, 50)) as pilot:
            await pilot.pause()
            feats, total = mi._map_feats_from_record(sc._make_demo_record())
            modal = sc.MapImageExportModal(
                feats=feats, total=total, title="pDemo",
                default_path=str(tmp_path / "m.png"))
            app.push_screen(modal)
            await wait_for_widget(pilot, "#mapimg-fmt")
            modal.query_one("#mapimg-fmt", Select).value = "svg"
            modal.query_one("#mapimg-size", Select).value = "800"
            modal.query_one("#mapimg-transparent", Checkbox).value = True
            modal._selected_dir = str(tmp_path)
            modal.query_one("#mapimg-filename", Input).value = "out.svg"
            modal._do_export()
            out = tmp_path / "out.svg"
            await wait_for_state(pilot, lambda: out.exists(), what="svg written")
            body = out.read_text()
            ET.fromstring(body)
            assert "<rect" not in body        # transparent → no backing rect

    async def test_extension_swaps_with_format(self, tmp_path):
        from textual.widgets import Input, Select
        app = sc.PlasmidApp()
        async with app.run_test(size=(170, 50)) as pilot:
            await pilot.pause()
            feats, total = mi._map_feats_from_record(sc._make_demo_record())
            modal = sc.MapImageExportModal(feats=feats, total=total,
                                           title="pDemo")
            app.push_screen(modal)
            await wait_for_widget(pilot, "#mapimg-filename")
            modal.query_one("#mapimg-filename", Input).value = "thing.png"
            modal.query_one("#mapimg-fmt", Select).value = "svg"
            await pilot.pause()
            assert modal.query_one("#mapimg-filename", Input).value == "thing.svg"

    async def test_action_pushes_modal_for_loaded_plasmid(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(170, 50)) as pilot:
            await pilot.pause()
            app._apply_record(sc._make_demo_record())
            await pilot.pause()
            app.action_export_map_image()
            await wait_for_modal(pilot, sc.MapImageExportModal)
            modal = app.screen
            assert modal._total > 0 and modal._feats

    async def test_action_no_record_warns_not_crash(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(170, 50)) as pilot:
            await pilot.pause()
            app._current_record = None
            app.action_export_map_image()      # must not raise / push a modal
            await pilot.pause()
            assert not isinstance(app.screen, sc.MapImageExportModal)


# ── Bulk / mass export (library marks + collection) ─────────────────────────

class TestBulkExportMapImages:
    """`_bulk_export_map_images` — the entry-list → folder image exporter that
    powers the library's 'export marked plasmids as images' flow."""

    def test_png_writes_all(self, tmp_path):
        res = sc._bulk_export_map_images(
            [_entry("pA"), _entry("pB")], tmp_path, "png", size=400)
        assert res["total"] == 2
        assert len(res["written"]) == 2 and not res["failures"]
        assert (tmp_path / "pA.png").exists()
        assert Image.open(tmp_path / "pB.png").size == (400, 400)

    def test_svg_writes(self, tmp_path):
        sc._bulk_export_map_images([_entry("pX")], tmp_path, "svg", size=400)
        ET.fromstring((tmp_path / "pX.svg").read_text())

    def test_bad_entry_recorded_not_raised(self, tmp_path):
        entries = [_entry("good"), {"id": "bad", "name": "bad", "gb_text": ""}]
        res = sc._bulk_export_map_images(entries, tmp_path, "png", size=400)
        assert len(res["written"]) == 1 and len(res["failures"]) == 1
        assert res["failures"][0]["name"] == "bad"
        assert (tmp_path / "good.png").exists()

    def test_filename_collision_deduped(self, tmp_path):
        sc._bulk_export_map_images([_entry("dup"), _entry("dup")],
                                   tmp_path, "png", size=400)
        names = sorted(p.name for p in tmp_path.glob("*.png"))
        assert names == ["dup.png", "dup_2.png"]

    def test_progress_callback_fires(self, tmp_path):
        calls = []
        sc._bulk_export_map_images(
            [_entry("p1"), _entry("p2")], tmp_path, "png",
            progress_cb=lambda i, t, n, ok: calls.append((i, t, ok)), size=400)
        assert calls[-1] == (2, 2, True)

    def test_render_opts_passed_through(self, tmp_path):
        # transparent → RGBA PNG with alpha-0 corner.
        sc._bulk_export_map_images([_entry("pT")], tmp_path, "png",
                                   size=400, transparent=True)
        img = Image.open(tmp_path / "pT.png")
        assert img.mode == "RGBA" and img.getpixel((0, 0))[3] == 0

    def test_unknown_fmt_rejected(self, tmp_path):
        with pytest.raises(ValueError):
            sc._bulk_export_map_images([_entry("p")], tmp_path, "bmp")

    def test_bulk_export_collection_supports_png_svg(self):
        # The collection bulk exporter learned png/svg via _BULK_EXPORT_FORMATS.
        assert "png" in sc._BULK_EXPORT_FORMATS
        assert "svg" in sc._BULK_EXPORT_FORMATS
        assert sc._BULK_EXPORT_FORMATS["png"][0] == "png"


class TestMarkedExportModal:
    async def test_marked_modal_exports_folder(self, tmp_path):
        from textual.widgets import Select
        app = sc.PlasmidApp()
        async with app.run_test(size=(170, 50)) as pilot:
            await pilot.pause()
            entries = [_entry("pM1"), _entry("pM2")]
            modal = sc.MarkedMapImageExportModal(entries,
                                                 default_dir=str(tmp_path))
            app.push_screen(modal)
            await wait_for_widget(pilot, "#btn-markexp")
            modal.query_one("#markexp-size", Select).value = "800"
            modal._selected_dir = str(tmp_path)
            modal._do_export(None)
            await wait_for_state(
                pilot,
                lambda: (tmp_path / "pM1.png").exists()
                and (tmp_path / "pM2.png").exists(),
                what="bulk pngs written")
            assert Image.open(tmp_path / "pM1.png").size == (800, 800)


# ── Label layout: on-band names, radial placement, title, padding ────────────

class TestLayout:
    """The 2026-07 relayout: names on the arc when they fit, external labels on
    an outer ellipse with non-crossing leaders, a title that drops below the
    circle when too long, and no text clipping the canvas."""

    def _texts(self, rec, *, title="", size=1000, sites=False):
        feats, total = mi._map_feats_from_record(rec)
        return mi._build_primitives(feats, total, title=title, size=size,
                                    show_sites=sites)["texts"]

    def test_every_text_carries_rotation(self):
        texts = self._texts(sc._make_demo_record(), title="pDemo")
        assert texts and all(len(t) == 8 for t in texts)

    def test_long_arc_short_name_is_on_band_rotated(self):
        # A feature spanning most of the plasmid with a short name → drawn ON
        # the band (rotated), never as an external horizontal label.
        rec = _rec("ACGT" * 250, [(50, 900, 1, "CDS", "bigGene", "#3366CC")])
        texts = self._texts(rec, size=1000)
        matches = [t for t in texts if t[2] == "bigGene"]
        assert matches, "feature label missing entirely"
        assert all(t[7] != 0 for t in matches)     # rotated to the arc tangent

    def test_short_tight_feature_is_external_horizontal(self):
        # A tiny feature can't hold its name → external label, horizontal.
        rec = _rec("ACGT" * 250, [(500, 510, 1, "CDS",
                                   "TinyFeatureLabelHere", "#3366CC")])
        texts = self._texts(rec, size=1000)
        m = [t for t in texts if t[2] == "TinyFeatureLabelHere"]
        assert m and all(t[7] == 0 for t in m)      # horizontal external label

    def test_title_centres_when_it_fits(self):
        texts = self._texts(sc._make_demo_record(), title="pX", size=1000)
        t = [x for x in texts if x[2] == "pX"][0]
        assert abs(t[1] - 500) < 0.15 * 1000        # near the centre (cy=500)

    def test_long_title_moves_below_the_circle(self):
        long_t = "LONG SPEC pDemoVector-Pdemo-GeneX-TermDemo-99"
        texts = self._texts(sc._make_demo_record(), title=long_t, size=1000)
        t = [x for x in texts if x[2] == long_t]
        assert t, "title text missing"
        assert t[0][1] > 500 + mi._R_BACKBONE * 1000   # below the ring bottom

    def test_no_label_clips_the_canvas(self):
        # Ample padding: every label (with its estimated width) stays inside the
        # canvas on a dense plasmid full of long names.
        specs = [(i * 30, i * 30 + 18, (1 if i % 2 else -1), "misc_feature",
                  f"very-long-feature-name-{i:02d}", "#8833AA")
                 for i in range(28)]
        size = 1400
        texts = self._texts(_rec("ACGT" * 500, specs), title="dense", size=size)
        for x, y, text, fs, _c, anchor, _w, _r in texts:
            assert 0 <= y <= size, f"{text!r} y={y} off-canvas"
            w = len(text) * fs * mi._CHAR_W
            if anchor == "start":
                assert x + w <= size + 3, f"{text!r} clips right edge"
            elif anchor == "end":
                assert x - w >= -3, f"{text!r} clips left edge"

    def test_rotated_on_band_label_renders_both_backends(self):
        # The on-band (rotated) label path has its own PNG code — a temp text
        # tile that is rotated + alpha-composited (`_draw_rotated_text`). Render
        # both back-ends for a big-arc/short-name plasmid so that whole path is
        # exercised, not just the primitive tuple.
        rec = _rec("ACGT" * 250, [(50, 900, 1, "CDS", "bigGene", "#3366CC")])
        feats, total = mi._map_feats_from_record(rec)
        svg = mi.render_plasmid_map_svg(feats, total, title="pRot", size=600)
        ET.fromstring(svg)
        assert "rotate(" in svg                       # tangent-rotated name
        png = mi.render_plasmid_map_png(feats, total, title="pRot", size=600)
        assert Image.open(BytesIO(png)).size == (600, 600)

    def test_dense_long_names_render_without_clip_or_crash(self):
        # Drive the real emit path (not just `_build_primitives`) on a dense
        # plasmid of long names + many cut sites: adaptive font, per-side
        # de-collision, and the downsample cap must all survive rendering.
        specs = [(i * 30, i * 30 + 18, (1 if i % 2 else -1), "misc_feature",
                  f"very-long-feature-name-{i:02d}", "#8833AA")
                 for i in range(40)]
        feats, total = mi._map_feats_from_record(_rec("ACGT" * 500, specs))
        sites = [{"type": "recut", "start": i * 61, "label": f"BsaI-{i}"}
                 for i in range(20)]
        ET.fromstring(mi.render_plasmid_map_svg(feats, total, sites=sites,
                                                title="dense", size=1200))
        png = mi.render_plasmid_map_png(feats, total, sites=sites,
                                        title="dense", size=1200)
        assert Image.open(BytesIO(png)).size == (1200, 1200)

    def test_long_title_below_renders(self):
        # The title-below-the-circle branch must render in both back-ends.
        long_t = "LONG SPEC pDemoVector-Pdemo-GeneX-TermDemo-99-extralong"
        feats, total = mi._map_feats_from_record(sc._make_demo_record())
        svg = mi.render_plasmid_map_svg(feats, total, title=long_t, size=700)
        ET.fromstring(svg)
        assert "pDemoVector" in svg
        assert Image.open(BytesIO(
            mi.render_plasmid_map_png(feats, total, title=long_t,
                                      size=700))).size == (700, 700)


# ── audit 2026-10-23: lane packing is O(F log F) and arcs are capped ────────

def _ref_assign_lanes(feats, total):
    """The pre-2026-10-23 pairwise algorithm, frozen as the reference."""
    import splicecraft_biology as B
    import splicecraft_mapimage as M

    def in_arc(x, start, length):
        return ((x - start) % total) < length

    order = sorted(range(len(feats)),
                   key=lambda i: B._feat_len(feats[i]["start"],
                                             feats[i]["end"], total),
                   reverse=True)
    lanes, result = [], {}
    for i in order:
        f = feats[i]
        s = int(f["start"])
        ln = B._feat_len(s, f["end"], total)
        placed = False
        for lane_idx, occ in enumerate(lanes):
            if all(not (in_arc(s, os_, ol) or in_arc(os_, s, ln))
                   for os_, ol in occ):
                occ.append((s, ln))
                result[i] = lane_idx
                placed = True
                break
        if not placed:
            if len(lanes) >= M._MAP_MAX_LANES:
                lanes[-1].append((s, ln))
                result[i] = len(lanes) - 1
            else:
                lanes.append([(s, ln)])
                result[i] = len(lanes) - 1
    return result


def test_lane_assignment_matches_the_pairwise_reference():
    import random
    import splicecraft_mapimage as M
    rng = random.Random(20261024)
    for trial in range(400):
        total = rng.choice([50, 400, 3000, 10_000])
        feats = []
        for _ in range(rng.randint(0, 120)):
            s = rng.randrange(-total, 2 * total)
            r = rng.random()
            if r < 0.08:
                e = s                                      # zero length
            elif r < 0.12:
                e = s + total * rng.randint(1, 3)          # full turn(s)
            elif r < 0.25:
                e = rng.randrange(0, max(1, s % total))    # origin wrap
            else:
                e = s + rng.randint(1, max(1, total // rng.choice([2, 5, 50])))
            feats.append({"start": s, "end": e})
        assert M._assign_lanes(feats, total) == _ref_assign_lanes(feats, total), trial


def test_lane_assignment_is_not_pairwise():
    """8,000 non-overlapping features: the pairwise scan took 6 s here (and
    21 s at 16k); the bisect packer takes ~15 ms. 100x headroom for a busy
    machine, still 4x under what the old algorithm needs."""
    import time
    import splicecraft_mapimage as M
    n = 8000
    feats = [{"start": i * 100, "end": i * 100 + 60} for i in range(n)]
    t0 = time.perf_counter()
    lanes = M._assign_lanes(feats, n * 100)
    dt = time.perf_counter() - t0
    assert set(lanes.values()) == {0}
    assert dt < 1.5, f"lane assignment took {dt:.2f} s for {n} features"


def test_a_feature_far_past_the_end_draws_one_turn():
    # 300 turns: the old renderer sampled every one of them (72k points here;
    # `end=1e9` was 48 million points and the 11.7 GB OOM). Kept small so the
    # test fails fast — not by exhausting memory — on code without the cap.
    import splicecraft_mapimage as M
    feats = [{"start": 10, "end": 10 + 5000 * 300, "strand": 1, "label": "x",
              "color": "color(40)"}]
    prims = M._build_primitives(feats, 5000, size=400)
    n_points = sum(len(poly) for poly, *_rest in prims["polys"])
    assert n_points < 2000, n_points
