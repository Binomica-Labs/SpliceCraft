"""Speed-ups must not change answers ([INV-214]).

Every optimisation in that pass replaced a hot path with a faster one that is
supposed to compute the SAME thing. These tests hold each one to its
definition: the fast form against the slow form it replaced, over inputs
chosen to reach the corners (IUPAC codes, case, characters outside the
alphabet, empty input), plus a golden digest of the aligner's output — whose
tie-breaking is load-bearing ([INV-212], [INV-213]) and must only ever change
on purpose.
"""
from __future__ import annotations

import hashlib
import itertools
import os
import random

import pytest

import splicecraft as sc
import splicecraft_biology as bio


# ══════════════════════════════════════════════════════════════════════════
# IUPAC compatibility fast path
# ══════════════════════════════════════════════════════════════════════════

def _slow_compatible(a, b):
    """The definition, exactly as it was before the lookup table."""
    if not a or not b:
        return False
    sa = bio._IUPAC_BASE_SET.get(a.upper())
    sb = bio._IUPAC_BASE_SET.get(b.upper())
    if sa is None or sb is None:
        return False
    return bool(sa & sb)


def test_iupac_fast_path_equals_the_definition():
    # Both cases, the gap glyph, junk, multi-character and empty input, and
    # two characters whose `.upper()` IS an ASCII IUPAC letter (long s -> S,
    # dotless i -> I): the table covers only the ASCII letters, so those must
    # still reach the slow path and answer the same.
    alphabet = (list("ACGTUMRWSYKVHDBNacgtumrwsykvhdbn-*. XZxz")
                + ["", "AC", "ac", "ſ", "ı", "ﬀ", "K"])
    bad = [(a, b) for a, b in itertools.product(alphabet, repeat=2)
           if _slow_compatible(a, b) != sc._iupac_compatible(a, b)]
    assert not bad, bad[:10]


# ══════════════════════════════════════════════════════════════════════════
# Myers pattern bit-vectors
# ══════════════════════════════════════════════════════════════════════════

def _slow_peq(pattern):
    peq = {}
    for j, ch in enumerate(pattern):
        bit = 1 << j
        for c in sc._MYERS_COMPAT_CHARS.get(ch, ()):
            peq[c] = peq.get(c, 0) | bit
    return peq


@pytest.mark.parametrize("alph", ["ACGT", "ACGTRYSWKMBDHVN", "ACGTN-*x"])
def test_peq_bitvectors_equal_the_per_position_definition(alph):
    rng = random.Random(7)
    for _ in range(300):
        pat = "".join(rng.choice(alph) for _ in range(rng.randint(0, 200)))
        assert sc._myers_build_peq(pat) == _slow_peq(pat), pat


# ══════════════════════════════════════════════════════════════════════════
# The aligner's output, pinned
# ══════════════════════════════════════════════════════════════════════════

# sha256 over 240 seeded alignments (see `_aligner_digest`). Computed on
# v1.3.6, before any of the speed-ups, and unchanged by them. If this fails,
# the engine now prefers a different one of several equal-cost alignments:
# reported variant positions move. Change it only deliberately.
_ALIGNER_DIGEST = (
    "1f0ea69ecdb2c04c48d154dce8f1c9d5a0c75ee27e9a654026cdde5a466eac5e")


def _aligner_digest(n: int = 240) -> str:
    rng = random.Random(2026)
    h = hashlib.sha256()
    for k in range(n):
        alph = "ACGTRYSWKMBDHVN" if k % 4 == 0 else "ACGT"
        if k % 3 == 0:
            a = "".join(rng.choice(alph) for _ in range(rng.randint(0, 80)))
            b = "".join(rng.choice(alph) for _ in range(rng.randint(0, 80)))
        else:
            base = "".join(rng.choice("ACGT")
                           for _ in range(rng.randint(100, 1200)))
            out = []
            for c in base:
                r = rng.random()
                if r < 0.01:
                    continue
                if r < 0.02:
                    out.append(rng.choice(alph))
                out.append(rng.choice(alph) if rng.random() < 0.01 else c)
            a, b = "".join(out), base
        ga, gb = sc._myers_align_global(a, b)
        h.update(ga.encode() + b"|" + gb.encode() + b"\n")
    return h.hexdigest()


def test_the_aligner_returns_the_same_alignments():
    assert _aligner_digest() == _ALIGNER_DIGEST


# ══════════════════════════════════════════════════════════════════════════
# Blob store: parallel reads, the ref cache, the mtime rule
# ══════════════════════════════════════════════════════════════════════════

GB = "LOCUS       test 20 bp DNA circular\nORIGIN\n        1 acgtacgtac gtacgtacgt\n//\n"


def test_blob_read_many_equals_reading_one_at_a_time():
    refs = [sc._blob_write(GB + "x" * i) for i in range(40)]
    missing = "0" * 64
    corrupt = sc._blob_write(GB + "corrupt-me")
    sc._blob_path(corrupt).write_text("not what the name says")
    sc._BLOB_VERIFIED_PATHS.discard(str(sc._blob_path(corrupt)))
    want = {r: sc._blob_read(r) for r in refs + [missing, corrupt]}
    got = sc._blob_read_many(refs + refs[:5] + [missing, corrupt])   # dupes
    assert got == want
    assert got[missing] is None and got[corrupt] is None


def test_rehydrating_a_library_reads_each_blob_once(monkeypatch):
    ref = sc._blob_write(GB)
    entries = [{"id": f"p{i}", "gb_ref": ref} for i in range(50)]
    calls = []
    real = sc._blob_read
    monkeypatch.setattr(sc, "_blob_read", lambda r: calls.append(r) or real(r))
    out = sc._rehydrate_collections([{"name": "A", "plasmids": entries[:25]},
                                     {"name": "B", "plasmids": entries[25:]}])
    assert calls == [ref]
    assert all(p["gb_text"] == GB for c in out for p in c["plasmids"])


def test_saving_a_loaded_text_does_not_rehash_it(monkeypatch):
    ref = sc._blob_write(GB)
    text = sc._blob_read(ref)                      # the text a load hands out
    hashed = []
    real = sc._blob_hash_bytes
    monkeypatch.setattr(sc, "_blob_hash_bytes",
                        lambda d: hashed.append(len(d)) or real(d))
    assert sc._blob_write(text) == ref
    assert hashed == []                            # reused, not re-hashed
    # An equal text that is a DIFFERENT object is hashed as before.
    assert sc._blob_write("".join(list(text))) == ref
    assert hashed


def test_a_cached_ref_still_catches_a_damaged_blob():
    text = sc._blob_read(sc._blob_write(GB))
    path = sc._blob_path(sc._blob_hash(GB))
    path.write_text(GB[:10])                       # truncated on disk
    assert sc._blob_write(text) == sc._blob_hash(GB)
    assert path.read_text() == GB                  # rewritten from the text


def test_reuse_refreshes_only_a_blob_near_the_gc_window():
    ref = sc._blob_write(GB)
    path = sc._blob_path(ref)
    st = path.stat()
    young = st.st_mtime - 60                       # well inside half the window
    os.utime(path, (young, young))
    sc._blob_write(GB)
    assert path.stat().st_mtime == pytest.approx(young, abs=1.0)
    old = st.st_mtime - sc._BLOB_GC_GRACE_SECONDS  # past half the window
    os.utime(path, (old, old))
    sc._blob_write(GB)
    assert path.stat().st_mtime > old + sc._BLOB_GC_GRACE_SECONDS / 2


# ══════════════════════════════════════════════════════════════════════════
# _safe_save_json: an unchanged file is left alone
# ══════════════════════════════════════════════════════════════════════════

def test_an_unchanged_save_writes_nothing(tmp_path):
    p = tmp_path / "t.json"
    sc._safe_save_json(p, [{"id": "a"}], "t")
    sc._safe_save_json(p, [{"id": "b"}], "t")      # changed: backs up
    before = sorted(x.name for x in tmp_path.iterdir())
    st = p.stat()
    sc._safe_save_json(p, [{"id": "b"}], "t")      # identical
    assert sorted(x.name for x in tmp_path.iterdir()) == before
    assert p.stat().st_mtime_ns == st.st_mtime_ns
    assert p.stat().st_ino == st.st_ino


def test_a_changed_save_still_writes_and_backs_up(tmp_path):
    p = tmp_path / "t.json"
    sc._safe_save_json(p, [{"id": "a"}], "t")
    sc._safe_save_json(p, [{"id": "b"}], "t")
    n = len(list(tmp_path.glob("t.json.bak.*")))
    sc._safe_save_json(p, [{"id": "c"}], "t")
    assert len(list(tmp_path.glob("t.json.bak.*"))) > n
    import json
    assert json.loads(p.read_text())["entries"] == [{"id": "c"}]


# ══════════════════════════════════════════════════════════════════════════
# Rotation picker: a full tie goes to the alignment with fewer gaps
# ══════════════════════════════════════════════════════════════════════════

def test_a_partial_read_just_across_the_origin_is_placed_whole():
    """1 kb read crossing bp 0 by 3 bp. The query-rotated candidate ties the
    target-rotated one exactly on matches and identity, but threads the 3 bp
    through a chance match with a 7 bp gap; list order picked it, and a
    perfect read came back with a deletion (seed 1 of the sweep that found
    it; 14 of 96 such reads were wrong on v1.3.6)."""
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(6000))
    read = ref[len(ref) - 997:] + ref[:3]
    r = sc._state._AGENT_HANDLERS["verify-against-reads"][0](
        None, {"reference": ref, "reads": [read], "circular": True,
               "read_quality": [[40] * len(read)]})
    assert r["read_quality"][0]["verdict"] == "clean", r["read_quality"]
    assert r["consensus"]["variants"] == []


# ══════════════════════════════════════════════════════════════════════════
# Sequence panel: the window + padding paints every row where it was
# ══════════════════════════════════════════════════════════════════════════

def _long_record(n=60000, seed=11):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(n))
    rec = SeqRecord(Seq(seq), id="LONG", name="LONG",
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    for k, s in enumerate(sorted(rng.sample(range(n - 1000), 40))):
        rec.features.append(SeqFeature(
            FeatureLocation(s, s + rng.randint(60, 900), strand=rng.choice([1, -1])),
            type=rng.choice(["CDS", "misc_feature", "promoter"]),
            qualifiers={"label": [f"f{k}"]}))
    return rec


def test_windowed_seq_text_has_the_full_texts_height_and_rows():
    from textual.content import Content
    rec = _long_record()
    seq = str(rec.seq)
    feats = [{"start": int(f.location.start), "end": int(f.location.end),
              "label": f.qualifiers["label"][0], "type": f.type,
              "color": "green", "strand": f.location.strand, "qualifiers": {}}
             for f in rec.features]
    for lw in (60, 147):
        full = Content.from_rich_text(
            sc._build_seq_text(seq, feats, line_width=lw)).split("\n")
        total = len(full)
        for vp in [(0, 40), (333, 420), (total - 30, total + 40),
                   (total + 5, total + 80)]:
            info: dict = {}
            win = Content.from_rich_text(sc._build_seq_text(
                seq, feats, line_width=lw, viewport_y_range=vp,
                window_info=info)).split("\n")
            assert info["above"] + len(win) + info["below"] == total, (lw, vp)
            for r in range(max(vp[0], info["above"]),
                           min(vp[1], info["above"] + len(win))):
                assert win[r - info["above"]].plain == full[r].plain, (lw, vp, r)


def test_the_panel_paints_the_rows_of_the_full_text():
    """End to end on the mounted widget: every visible row of `#seq-view`,
    padding and all, is the row the old full text had at that offset."""
    import asyncio
    from textual.content import Content
    from textual.containers import ScrollableContainer
    from textual.geometry import Region
    from textual.widgets import Static

    async def run():
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            app._apply_record(_long_record())
            for _ in range(3):
                await pilot.pause()
            panel = app.query_one(sc.SequencePanel)
            scroll = panel.query_one("#seq-scroll", ScrollableContainer)
            view = panel.query_one("#seq-view", Static)
            for target in (0, 400, int(scroll.max_scroll_y)):
                scroll.scroll_to(y=target, animate=False, force=True)
                for _ in range(3):
                    await pilot.pause()
                panel._refresh_view()
                await pilot.pause()
                assert view.styles.padding.top > 0 or target == 0
                disp_seq, disp_feats = panel._get_rotated_state()
                num_w = len(str(len(panel._seq)))
                lw = max(20, panel._seq_render_width() - (num_w + 2))
                full = Content.from_rich_text(sc._build_seq_text(
                    disp_seq, disp_feats, line_width=lw,
                    cursor_pos=panel._abs_to_disp(panel._cursor_pos),
                    show_connectors=panel._show_connectors)).split("\n")
                y0 = int(scroll.scroll_y)
                width = view.outer_size.width
                for y in range(y0, y0 + scroll.size.height):
                    if y >= len(full):
                        break
                    # `render_lines` paints in WIDGET coordinates, padding
                    # included — what reaches the screen at that offset.
                    (strip,) = view.render_lines(Region(0, y, width, 1))
                    assert strip.text.rstrip() == full[y].plain.rstrip(), \
                        (target, y)
    asyncio.run(run())


# ══════════════════════════════════════════════════════════════════════════
# Hardening: the fast paths degrade instead of failing
# ══════════════════════════════════════════════════════════════════════════

def test_blob_reads_fall_back_when_no_thread_can_start(monkeypatch):
    import concurrent.futures as cf
    refs = [sc._blob_write(GB + "y" * i) for i in range(40)]

    class _NoThreads:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            raise RuntimeError("can't start new thread")

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(cf, "ThreadPoolExecutor", _NoThreads)
    got = sc._blob_read_many(refs)
    assert got == {r: GB + "y" * i for i, r in enumerate(refs)}


def test_the_ref_cache_is_bounded(monkeypatch):
    monkeypatch.setattr(sc, "_BLOB_REF_CACHE_MAX", 3)
    texts = [GB + "z" * i for i in range(7)]
    refs = [sc._blob_write(t) for t in texts]
    assert len(sc._BLOB_REF_BY_TEXT) <= 3
    assert [sc._blob_read(r) for r in refs] == texts


def test_refinement_leaves_a_window_it_cannot_align():
    # A character outside the scoring alphabet makes Biopython refuse the
    # window; the rows must come back exactly as given, not raise.
    aq = "ACGTACGTAC*-GTACG--TACGTACGTACGTACGTACGTACGTAC"
    at = "ACGTACGTACGAGTACGTTTACGTACGTACGTACGTACGTACGTAC"
    assert sc._refine_indel_clusters(aq, at) == (aq, at)


def test_global_alignment_survives_odd_inputs():
    rng = random.Random(99)
    shapes = []
    for k in range(600):
        n = rng.randint(1, 400)
        kind = k % 6
        if kind == 0:
            t = "A" * n
            q = "A" * max(1, n - rng.randint(0, 5))
        elif kind == 1:
            t = "".join(rng.choice("ACGT") for _ in range(n))
            q = "N" * max(1, n // 2)
        elif kind == 2:
            t = "".join(rng.choice("ACGT") for _ in range(n))
            q = t[rng.randint(0, n - 1):]
        elif kind == 3:
            t = "".join(rng.choice("ACGTRYN") for _ in range(n))
            q = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 400)))
        elif kind == 4:
            t = "".join(rng.choice("ACGT") for _ in range(n))
            q = "".join(c for c in t if rng.random() > 0.05)
            q = q or "A"
        else:
            t = "".join(rng.choice("AC") for _ in range(n))
            q = "".join(rng.choice("AC") for _ in range(rng.randint(1, 60)))
        shapes.append((q, t))
    for q, t in shapes:
        for free in (False, True):
            r = sc._pairwise_align(q, t, free_read_ends=free)
            aq, at = r["aligned_q"], r["aligned_t"]
            assert len(aq) == len(at)
            assert aq.replace("-", "") == q and at.replace("-", "") == t
            assert (r["n_matches"] + r["n_mismatches"] + r["n_gap_cols"]
                    == len(aq))
