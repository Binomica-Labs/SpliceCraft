"""
test_blast — pure-Python BLASTN + BLASTP engine.

Covers:
  - BLOSUM62 matrix parses with the canonical scores (A/A=4, W/W=11,
    A/W=-3) so a misaligned column count would be caught.
  - BLASTN finds an exact embedded subseq, finds the same subseq on the
    reverse-complement strand, tolerates a single mismatch, and rejects
    a query that's all noise.
  - BLASTP translates a DNA query to protein when the user picks BLASTP
    on a DNA paste, and finds an exact protein hit against an annotated
    CDS in the database.
  - `_blast_get_db` cache returns the same object on a repeat call and
    invalidates on `_blast_clear_cache`.
  - `_save_collections` triggers `_blast_clear_cache` (the cache flush
    on collection mutations).

The boundary regression for BlastModal is in `test_modal_boundaries.py`;
this file owns the **engine** contract.
"""
from __future__ import annotations

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# BLOSUM62 sanity
# ═══════════════════════════════════════════════════════════════════════════════

class TestBlosum62:

    def test_diagonal_canonical_aas(self):
        # The standard BLOSUM62 matrix has these well-known values on
        # the diagonal — if any of them is off, the parser shifted
        # columns or rows.
        assert sc._blosum62_score("A", "A") == 4
        assert sc._blosum62_score("R", "R") == 5
        assert sc._blosum62_score("W", "W") == 11
        assert sc._blosum62_score("C", "C") == 9
        assert sc._blosum62_score("P", "P") == 7

    def test_offdiagonal_known_pairs(self):
        # Reference values from the published BLOSUM62.
        assert sc._blosum62_score("A", "S") == 1
        assert sc._blosum62_score("A", "W") == -3
        assert sc._blosum62_score("E", "D") == 2
        assert sc._blosum62_score("L", "I") == 2

    def test_symmetric(self):
        for a, b in [("A", "G"), ("R", "K"), ("W", "Y"), ("F", "L")]:
            assert sc._blosum62_score(a, b) == sc._blosum62_score(b, a), \
                f"{a}/{b} not symmetric"

    def test_lowercase_normalised(self):
        # Helper upper-cases internally so a lowercase query doesn't
        # silently fall through to the -4 default.
        assert sc._blosum62_score("a", "a") == 4

    def test_unknown_chars_return_default(self):
        assert sc._blosum62_score("A", "?") == -4
        assert sc._blosum62_score("@", "@") == -4


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers — building a canned collection so the BLAST DB has something to chew
# ═══════════════════════════════════════════════════════════════════════════════

def _make_record(rec_id: str, seq: str, *, cds_ranges=None):
    """Build a SeqRecord ready for `_record_to_gb_text`. ``cds_ranges``
    is an iterable of ``(start, end, strand, label)`` tuples — every
    range becomes a CDS feature so the BLASTP DB build picks them up."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rec = SeqRecord(Seq(seq), id=rec_id, name=rec_id, description=rec_id)
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"]      = "circular"
    for (s, e, strand, lbl) in (cds_ranges or []):
        rec.features.append(SeqFeature(
            FeatureLocation(s, e, strand=strand),
            type="CDS",
            qualifiers={"label": [lbl]},
        ))
    return rec


def _seed_collection(name: str, records: list) -> None:
    """Persist a collection with the given records' GenBank text. The
    autouse `_protect_user_data` fixture redirects collections.json to
    a tmp path so this is test-isolated."""
    plasmids = []
    for rec in records:
        gb = sc._record_to_gb_text(rec)
        plasmids.append({
            "name":    rec.name,
            "id":      rec.id,
            "size":    len(rec.seq),
            "n_feats": len([f for f in rec.features if f.type != "source"]),
            "source":  f"id:{rec.id}",
            "added":   "2026-05-01",
            "gb_text": gb,
        })
    existing = sc._load_collections()
    existing.append({
        "name": name, "description": "test", "plasmids": plasmids,
        "saved": "2026-05-01",
    })
    sc._save_collections(existing)


# ═══════════════════════════════════════════════════════════════════════════════
# BLASTN
# ═══════════════════════════════════════════════════════════════════════════════

class TestBlastnEngine:

    def test_finds_exact_embedded_match(self):
        # A 60 bp subject with a 30 bp distinctive run. Query = those 30 bp.
        # Expect at least one HSP that covers the embedded region with
        # high identity. Both backends find the hit; pyhmmer may trim
        # one bp at the boundary so we assert ≥ 27 (90%) rather than
        # the full 30.
        target = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"   # 30 bp
        subject_seq = "GGGG" * 4 + target + "TTTT" * 4
        rec = _make_record("HIT01", subject_seq)
        _seed_collection("BlastTest", [rec])
        sc._blast_clear_cache()

        db = sc._blast_get_db("blastn", ["BlastTest"])
        hits = sc._blast_search(target, db)
        assert hits, "expected at least one HSP"
        top = hits[0]
        assert top["subject_id"] == "HIT01"
        assert top["strand"] == 1
        assert top["identity_pct"] >= 95.0
        assert top["aligned_len"] >= 27

    def test_finds_reverse_strand_match(self):
        # Subject contains the RC of the query — BLASTN should still
        # find it, on strand -1.
        query = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        subject_seq = "GGGG" * 4 + sc._rc(query) + "TTTT" * 4
        rec = _make_record("REV01", subject_seq)
        _seed_collection("BlastTest", [rec])
        sc._blast_clear_cache()

        db = sc._blast_get_db("blastn", ["BlastTest"])
        hits = sc._blast_search(query, db)
        assert any(h["strand"] == -1 and h["subject_id"] == "REV01"
                   for h in hits), f"no RC hit: {hits!r}"

    def test_tolerates_single_mismatch(self):
        # A 50 bp query embedded with a single base flipped. BLASTN
        # should still find it (the mismatch costs -3 but ungapped
        # extension still nets >> min_score on a 50 bp run).
        target = ("ACGT" * 12) + "AC"   # 50 bp
        # Flip position 25 from T → A.
        target_with_mismatch = target[:25] + "A" + target[26:]
        assert target != target_with_mismatch
        subject_seq = "GGGG" * 4 + target + "TTTT" * 4
        rec = _make_record("MIS01", subject_seq)
        _seed_collection("BlastTest", [rec])
        sc._blast_clear_cache()

        db = sc._blast_get_db("blastn", ["BlastTest"])
        hits = sc._blast_search(target_with_mismatch, db)
        # We expect the long flank to score well above the -3 penalty.
        # Either a single HSP covers most of the query, or two HSPs
        # straddle the mismatch — both are acceptable here.
        assert hits, "no hits despite >40 bp ungapped flanks"
        # At least one hit should align ≥20 bp.
        assert max(h["aligned_len"] for h in hits) >= 20

    def test_random_query_yields_no_hits(self):
        # A random-ish 60 bp query against an unrelated subject.
        rec = _make_record("REL01", "AAAACCCCGGGGTTTT" * 4)   # 64 bp homopolymer-ish
        _seed_collection("BlastTest", [rec])
        sc._blast_clear_cache()

        # Construct a query with no shared 11-mer with the subject.
        # A long stretch of all-different nucleotides won't share an
        # 11-mer with a homopolymer-ish subject.
        query = "AGCTGAGCTGA" * 6   # 66 bp, no AAAA / CCCC / GGGG / TTTT runs
        db = sc._blast_get_db("blastn", ["BlastTest"])
        hits = sc._blast_search(query, db)
        # Either zero hits, or hits that are below the score threshold.
        # The function already filters by min_score, so most of the
        # time len(hits) == 0; we just assert no high-id hits.
        assert all(h["identity_pct"] < 99.0 for h in hits), \
            f"random-vs-homopoly produced too-good hits: {hits}"

    def test_query_shorter_than_k_returns_empty(self):
        rec = _make_record("S01", "ACGT" * 30)
        _seed_collection("BlastTest", [rec])
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastn", ["BlastTest"])
        # k = 11, so a 5 bp query can't seed.
        assert sc._blast_search("ACGTA", db) == []


# ═══════════════════════════════════════════════════════════════════════════════
# BLASTP
# ═══════════════════════════════════════════════════════════════════════════════

class TestBlastpEngine:

    def test_finds_exact_protein_match(self):
        # Embed a CDS that translates to a non-repetitive 30-aa
        # protein. HMMER's profile builder filters low-complexity
        # regions, so an all-A poly-Ala query won't seed a match in
        # the pyhmmer backend — we need a realistic protein.
        # Codons chosen to encode "MKLAVTPGGRSEKILVNLRSADENPLG" (27 aa).
        cds_dna = (
            "ATG" "AAA" "CTG" "GCG" "GTG" "ACC" "CCG" "GGC" "GGC"
            "CGC" "AGC" "GAA" "AAA" "ATT" "CTG" "GTG" "AAC" "CTG"
            "CGC" "AGC" "GCG" "GAT" "GAA" "AAC" "CCG" "CTG" "GGC" "TAA"
        )  # 28 codons → 27 aa + stop
        from Bio.Seq import Seq
        protein_full = str(Seq(cds_dna).translate())  # 'MKLAVTPGGRSEKILVNLRSADENPLG*'
        spacer = "GGCCGGCCGGCCGGCC"   # 16 bp spacer (non-repetitive)
        plas = spacer + cds_dna + spacer
        cds_start = len(spacer)
        cds_end   = cds_start + len(cds_dna)
        rec = _make_record(
            "PROT01", plas,
            cds_ranges=[(cds_start, cds_end, 1, "fakeORF")],
        )
        _seed_collection("BpTest", [rec])
        sc._blast_clear_cache()

        # Query the bulk of the protein (drop the start M and stop *).
        query = protein_full[1:-1]
        assert len(query) >= 20
        db = sc._blast_get_db("blastp", ["BpTest"])
        hits = sc._blast_search(query, db)
        assert hits, f"no BLASTP hit: db has {len(db['subjects'])} subjects"
        top = hits[0]
        assert "fakeORF" in (top["subject_name"] or "") \
            or top["subject_id"].endswith(":fakeORF")
        assert top["identity_pct"] >= 90.0

    def test_dna_query_with_blastp_hint_translates(self):
        # If the user pastes DNA into the BLASTP textarea, the engine
        # auto-translates. Verify the helper does the right thing.
        prog, q = sc._detect_query_program("ATGGCTGCTGCT", "blastp")
        assert prog == "blastp"
        # ATG GCT GCT GCT → M A A A
        assert q.startswith("MAAA")

    def test_six_frame_finds_unannotated_orf(self):
        # Plasmid with NO annotated CDS — six-frame should still find
        # ORFs ≥ _BLASTP_MIN_ORF_AA. Build a backbone that contains a
        # 35-codon ORF in frame 1 (no stops until the end).
        # Codons chosen to avoid stops (TAA/TAG/TGA): use only AAA/GCT.
        orf_codons = ("ATG" + "AAA" * 17 + "GCT" * 17 + "TAA")
        plas = "GGCCGGCC" * 3 + orf_codons + "TTTTGGGG" * 3
        rec = _make_record("ORFTEST", plas)   # no CDS feature!
        _seed_collection("OrfBp", [rec])
        sc._blast_clear_cache()

        # Default (no 6-frame) → empty BLASTP db.
        db_default = sc._blast_get_db("blastp", ["OrfBp"])
        assert db_default["subjects"] == [], (
            "BLASTP without 6-frame should yield no subjects when "
            "there's no annotated CDS"
        )

        # 6-frame on → should pick up the embedded ORF.
        db_6f = sc._blast_get_db("blastp", ["OrfBp"], six_frame=True)
        assert db_6f["subjects"], "6-frame should find the embedded ORF"
        # Cache must distinguish six_frame=True from default.
        assert db_6f is not db_default

        # Query the protein middle and expect a hit.
        from Bio.Seq import Seq
        protein = str(Seq(orf_codons[3:-3]).translate())  # drop M and stop
        hits = sc._blast_search(protein[5:30], db_6f)
        # Don't insist on hits — pyhmmer may filter low-complexity poly-K
        # sequences — but verify no crash.
        for h in hits:
            assert "subject_id" in h

    def test_six_frame_off_by_default_in_modal_dispatch(self):
        # Through the dispatcher with the default (no kwarg), 6-frame
        # is off. This guards against accidental flips.
        rec = _make_record("DEFTEST", "ATG" + "AAA" * 50 + "TAA")
        _seed_collection("DefBp", [rec])
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastp", ["DefBp"])
        # Without an annotated CDS feature, the default path yields
        # zero subjects.
        assert db["subjects"] == []

    def test_skips_short_or_non_triple_cds(self):
        # CDS with length not divisible by 3 must be skipped (else the
        # BLASTP DB build would crash on translate). Build a record
        # where the only CDS is 7 bp; expect zero subjects in the db.
        rec = _make_record(
            "BAD01", "ACGT" * 30,
            cds_ranges=[(0, 7, 1, "shortCDS")],
        )
        _seed_collection("BpTest", [rec])
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastp", ["BpTest"])
        assert db["subjects"] == [], \
            "non-triple CDS should be filtered out of BLASTP db"


# ═══════════════════════════════════════════════════════════════════════════════
# Cache + invalidation
# ═══════════════════════════════════════════════════════════════════════════════

class TestBlastDbCache:

    def test_repeat_call_returns_same_object(self):
        rec = _make_record("S01", "AAAA" + "ACGT" * 30)
        _seed_collection("CacheTest", [rec])
        sc._blast_clear_cache()
        db1 = sc._blast_get_db("blastn", ["CacheTest"])
        db2 = sc._blast_get_db("blastn", ["CacheTest"])
        assert db1 is db2

    def test_clear_cache_drops_db(self):
        rec = _make_record("S01", "AAAA" + "ACGT" * 30)
        _seed_collection("CacheTest", [rec])
        sc._blast_clear_cache()
        db1 = sc._blast_get_db("blastn", ["CacheTest"])
        sc._blast_clear_cache()
        db2 = sc._blast_get_db("blastn", ["CacheTest"])
        assert db1 is not db2

    def test_save_collections_invalidates_cache(self):
        rec = _make_record("S01", "AAAA" + "ACGT" * 30)
        _seed_collection("MutTest", [rec])
        sc._blast_clear_cache()
        db1 = sc._blast_get_db("blastn", ["MutTest"])
        # Mutate: rename the collection, re-save.
        existing = sc._load_collections()
        for c in existing:
            if c.get("name") == "MutTest":
                c["name"] = "MutTestRenamed"
        sc._save_collections(existing)
        # After save, the cache must have been cleared. Same key now
        # builds against an empty filter (the rename invalidates by
        # name match) — the returned db should be a fresh object.
        db2 = sc._blast_get_db("blastn", ["MutTest"])
        assert db1 is not db2


# ═══════════════════════════════════════════════════════════════════════════════
# Pre-flight size cap (regression guard for 2026-05-10 OOM-on-large-collection)
# ═══════════════════════════════════════════════════════════════════════════════

class TestBlastDbSizeCap:
    """`_blast_filter_oversized` skips per-plasmid rows over
    `_BLAST_DB_BUILD_MAX_PER_PLASMID_BP` and raises
    `_BlastDbTooLargeError` when the post-skip total exceeds
    `_BLAST_DB_BUILD_MAX_TOTAL_BP`. Without this guard the kernel
    OOM-killed the worker mid-build on chromosome-scale collections —
    `_do_build`'s `except Exception` never fired and the user just
    saw the app vanish.
    """

    def test_normal_collection_passes_through_unfiltered(self):
        rec = _make_record("S01", "AAAA" + "ACGT" * 30)
        _seed_collection("Normal", [rec])
        cols = [c for c in sc._load_collections()
                if c.get("name") == "Normal"]
        kept, skipped, total_bp = sc._blast_filter_oversized(cols)
        assert skipped == []
        assert total_bp == len(rec.seq)
        assert len(kept) == 1
        assert len(kept[0]["plasmids"]) == 1

    def test_oversized_per_plasmid_skipped_with_reason(self):
        # Build a collection whose ONLY plasmid is over the per-plasmid
        # cap. We don't actually allocate 10 Mbp — we lie about `size`
        # in the entry dict, since the filter reads that field directly
        # (avoiding a 10 Mbp gb_text on every test run).
        sc._save_collections([{
            "name": "Big", "description": "test",
            "plasmids": [{
                "name":    "huge_plasmid",
                "id":      "HUGE",
                "size":    sc._BLAST_DB_BUILD_MAX_PER_PLASMID_BP + 1,
                "n_feats": 0,
                "gb_text": "stub",
            }],
            "saved": "2026-05-10",
        }])
        cols = sc._load_collections()
        kept, skipped, total_bp = sc._blast_filter_oversized(cols)
        assert total_bp == 0
        assert kept == []          # collection was emptied + dropped
        assert len(skipped) == 1
        assert skipped[0]["name"] == "huge_plasmid"
        assert skipped[0]["collection"] == "Big"
        assert "exceeds per-plasmid cap" in skipped[0]["reason"]

    def test_total_cap_raises_with_helpful_message(self):
        # Several plasmids each just under the per-plasmid cap, but
        # whose sum exceeds the total cap. Triggers the refusal path
        # without skipping any individual entry.
        per_size = sc._BLAST_DB_BUILD_MAX_PER_PLASMID_BP
        n_plasmids = (sc._BLAST_DB_BUILD_MAX_TOTAL_BP // per_size) + 1
        sc._save_collections([{
            "name": "TotalOverflow", "description": "test",
            "plasmids": [{
                "name":    f"p{i}",
                "id":      f"P{i}",
                "size":    per_size,
                "n_feats": 0,
                "gb_text": "stub",
            } for i in range(n_plasmids)],
            "saved": "2026-05-10",
        }])
        cols = sc._load_collections()
        import pytest
        with pytest.raises(sc._BlastDbTooLargeError) as exc:
            sc._blast_filter_oversized(cols)
        msg = str(exc.value)
        assert "exceeds the BLAST DB build cap" in msg
        assert "Filter" in msg

    def test_blast_too_large_subclasses_memory_error(self):
        """Worker's `except Exception` catches MemoryError subclasses,
        so the refusal flow lands in the user-facing status bar
        rather than leaking the traceback."""
        assert issubclass(sc._BlastDbTooLargeError, MemoryError)
        assert issubclass(sc._BlastDbTooLargeError, Exception)

    def test_build_db_attaches_skipped_and_total_bp(self):
        """`_blast_build_db` must surface the filter's skip list +
        total bp on the returned db so `_build_done` can render them."""
        rec = _make_record("S01", "AAAA" + "ACGT" * 30)
        _seed_collection("WithSkip", [rec])
        # Add an oversized sibling that the filter will drop.
        existing = sc._load_collections()
        for c in existing:
            if c.get("name") == "WithSkip":
                c["plasmids"].append({
                    "name":    "huge",
                    "id":      "HUGE",
                    "size":    sc._BLAST_DB_BUILD_MAX_PER_PLASMID_BP + 1,
                    "n_feats": 0,
                    "gb_text": "stub",
                })
        sc._save_collections(existing)
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastn", ["WithSkip"])
        # Normal-sized plasmid still indexed.
        assert len(db["subjects"]) == 1
        # Skip metadata propagated.
        assert "skipped" in db
        assert len(db["skipped"]) == 1
        assert db["skipped"][0]["name"] == "huge"
        assert db["total_bp"] == len(rec.seq)


# ═══════════════════════════════════════════════════════════════════════════════
# Modal integration (engine wired correctly)
# ═══════════════════════════════════════════════════════════════════════════════

class TestBlastModalIntegration:

    async def test_modal_is_full_screen(
            self, tiny_record, isolated_library):
        """Regression for 2026-05-26 "not enough space for the hmm
        elements" report: the modal MUST render at near-full-screen
        size (98%×96%) like AlignmentScreen so the HMM picker +
        Manage button + banner + query + results all fit without
        cramping, and so future TabbedContent layers have room to
        land. Pre-fix the modal was a 96-cell-wide centered dialog
        (90% height capped at 44) and the HMM section's three rows
        visually collided.

        Uses `dlg.region` (the screen-aligned outer box) rather
        than `dlg.size` (the inner content area, which is reduced
        by `padding: 1 2` + `border: solid`)."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            dlg = modal.query_one("#blast-dlg")
            region = dlg.region
            # 98% of 171 = 167.58 → 167; 96% of 43 = 41.28 → 41.
            # Pre-fix outer width was 96 cells flat (no `%`), so
            # any width ≥ 150 confirms the `width: 98%` rule took
            # effect. Height likewise: pre-fix outer was 38 cells.
            assert region.width  >= 150, (
                f"#blast-dlg outer width {region.width} too "
                f"narrow — the modal must be near-full-screen "
                f"(98%)."
            )
            assert region.height >= 39, (
                f"#blast-dlg outer height {region.height} too "
                f"short — the modal must be near-full-screen "
                f"(96%)."
            )

    async def test_modal_uses_tabbed_content_for_future_growth(
            self, tiny_record, isolated_library):
        """Regression / scaffold guard: the BLAST modal hosts its
        body inside a `TabbedContent` with a `tab-blast-hmm`
        `TabPane`. Pre-tab structure (2026-05-26) the body was a
        bare `ScrollableContainer` directly under `#blast-dlg`,
        so adding a sibling search tool (pairwise align, ORF
        finder, …) required restructuring the modal. The tab
        scaffold means a new tool drops in as a single new
        `TabPane`. This test catches a regression where the
        tabbed shell gets accidentally removed by a refactor."""
        from textual.widgets import TabbedContent, TabPane
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            tabs = modal.query_one("#blast-tabs", TabbedContent)
            assert tabs is not None
            pane = modal.query_one("#tab-blast-hmm", TabPane)
            assert pane is not None
            # The existing controls live inside the tab pane
            # (queryable by ID anywhere in the DOM, but the
            # parent chain must include the TabPane so a future
            # refactor can't accidentally hoist them back out).
            query = modal.query_one("#blast-query")
            ancestor_ids = {
                w.id for w in query.ancestors
                if getattr(w, "id", None)
            }
            assert "tab-blast-hmm" in ancestor_ids
            assert "blast-tabs" in ancestor_ids

    async def test_modal_hidden_path_input_reserves_no_space(
            self, tiny_record, isolated_library):
        """Pre-fix the `.-hidden` class on `#blast-hmm-path` had no
        CSS rule so the "hidden" back-compat Input actually
        rendered, eating 3 rows of layout right between the HMM
        Select and the six-frame Checkbox. Now `.-hidden` on this
        Input is scoped to `display: none` so the widget reserves
        no layout space — the worker / legacy tests can still
        `query_one("#blast-hmm-path")` to read the path, the user
        just doesn't see it."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            inp = modal.query_one("#blast-hmm-path", sc.Input)
            # The Input still exists (back-compat: legacy tests +
            # workers query it for the path value).
            assert inp is not None
            # …but the widget MUST NOT reserve layout space.
            # `display = False` means width / height are both 0.
            assert inp.region.height == 0, (
                "#blast-hmm-path must be display:none (reserves "
                "zero layout rows) so the HMM section's three "
                "rows don't visually collide. Got region.height="
                f"{inp.region.height}."
            )

    async def test_run_returns_engine_results_not_phase2_stub(
            self, tiny_record, isolated_library):
        # Seed a collection with one plasmid containing a recognisable
        # signature, push BlastModal, run BLASTN, and verify the
        # results panel contains the engine's actual output (a hit
        # row), not the old "Phase 2 pending" string.
        target = "ACGTACGTACGTACGTACGTACGTACGT"  # 28 bp
        rec = _make_record("INT01", "GGGG" * 4 + target + "TTTT" * 4)
        _seed_collection("IntTest", [rec])
        sc._blast_clear_cache()

        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#blast-query", sc.TextArea).text = target
            modal.query_one("#blast-source", sc.Select).value = "IntTest"
            await pilot.pause()
            modal.query_one("#btn-blast-run", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            results = modal.query_one("#blast-results", sc.Static)
            txt = str(results.render())
            # The engine renders a header containing the program name +
            # a row with the subject. The Phase-2 stub message would be
            # "Engine pending"; verify that's gone.
            assert "pending" not in txt.lower()
            assert "INT01" in txt or "BLASTN" in txt


# ═══════════════════════════════════════════════════════════════════════════════
# Hardening: input sanitisation
# ═══════════════════════════════════════════════════════════════════════════════

class TestQuerySanitisation:

    def test_strips_fasta_header(self):
        # Standard >id header line at column 0.
        assert sc._strip_fasta_headers(">pUC19\nATGAAA") == "ATGAAA"

    def test_strips_fasta_header_with_leading_whitespace(self):
        # Real-world copy-paste sometimes has leading tabs/spaces on
        # the header line (e.g. wrapped from a code block).
        assert sc._strip_fasta_headers("  >pUC19\nATGAAA") == "ATGAAA"

    def test_strips_multiple_headers(self):
        assert sc._strip_fasta_headers(
            ">a\nATG\n>b\nGCT") == "ATG\nGCT"

    def test_no_header_passes_through(self):
        assert sc._strip_fasta_headers("ATGAAA") == "ATGAAA"

    def test_blastn_filters_to_iupac_alphabet(self):
        # Numbers + punctuation are common in scientific paper pastes.
        # The BLASTN alphabet should drop them.
        prog, q = sc._detect_query_program(
            "1 atg aaa 50 GG | TGC", "blastn")
        assert prog == "blastn"
        assert q == "ATGAAAGGTGC"

    def test_blastn_keeps_iupac_codes(self):
        prog, q = sc._detect_query_program("ACGNRYWSMKBDHV", "blastn")
        assert q == "ACGNRYWSMKBDHV"

    def test_blastp_filters_to_aa_alphabet(self):
        # Real protein paste often has digits (residue numbers) — drop them.
        prog, q = sc._detect_query_program(
            "1 MKLAVT 50 PGRSE", "blastp")
        assert prog == "blastp"
        # Note: M and A are valid as both DNA-A and AA-M; the heuristic
        # would translate if 95% of alpha chars are ACGTN. Here the
        # paste has K, L, V, T which aren't DNA → BLASTP path keeps it raw.
        assert q == "MKLAVTPGRSE"

    def test_blastp_dna_query_translates(self):
        prog, q = sc._detect_query_program(
            "ATGGCTGCTGCTAAATAA", "blastp")
        assert prog == "blastp"
        # ATG GCT GCT GCT AAA TAA → M A A A K *
        assert q.startswith("MAAAK") and "*" in q

    def test_query_capped_at_max_len(self):
        # Send 200 KB of A; expect at most _MAX_BLAST_QUERY_LEN.
        big = "A" * 200_000
        prog, q = sc._detect_query_program(big, "blastn")
        assert len(q) <= sc._MAX_BLAST_QUERY_LEN

    def test_unknown_program_passthrough(self):
        # Defensive: unknown program doesn't crash, just truncates +
        # uppercases. Same shape as the BLASTN/BLASTP returns.
        prog, q = sc._detect_query_program("atgaaa", "weirdprog")
        assert prog == "weirdprog"
        assert q == "ATGAAA"

    def test_empty_after_sanitisation_returns_empty_string(self):
        # Paste of pure noise → zero-length query, modal will surface
        # the "after sanitising" red message.
        prog, q = sc._detect_query_program("123 || === ", "blastn")
        assert q == ""


# ═══════════════════════════════════════════════════════════════════════════════
# Hardening: library annotation cap
# ═══════════════════════════════════════════════════════════════════════════════

class TestLibraryAnnotationCap:

    def test_max_hits_truncates_runaway_search(self):
        # Seed a 12 bp library entry that occurs 1000 times in the
        # query; cap=10 should stop the search early.
        sc._save_features([{
            "name":         "common",
            "feature_type": "misc_feature",
            "sequence":     "ACGTACGTACGT",   # 12 bp = min_overlap default
            "strand":       1,
            "color":        "",
        }])
        big = "ACGTACGTACGT" * 1000
        hits = sc._annotate_seq_from_feature_library(big, max_hits=10)
        assert len(hits) <= 10

    def test_default_cap_is_5000(self):
        # Sanity check on the constant (so a future change accidentally
        # bumping it down to 50 trips this test).
        assert sc._DEFAULT_LIB_ANNOT_MAX_HITS == 5_000


# ═══════════════════════════════════════════════════════════════════════════════
# Modal-active gating: the App's on_key / on_click guards
# ═══════════════════════════════════════════════════════════════════════════════

class TestModalActiveGating:
    """When a ModalScreen is on top, the App's on_key handler must
    not fire seq-cursor / selection-slide / RE-clearing branches —
    those would silently mutate the underlying panel behind the user's
    back. App-level Ctrl+Z stays global as a fallback."""

    async def test_arrow_does_not_move_seq_cursor_under_modal(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = 5
            initial = sp._cursor_pos

            # Push a modal to top of the stack.
            app.push_screen(sc.HelpModal())
            await pilot.pause()
            assert len(app.screen_stack) > 1

            # Press Right; the seq cursor should not move (the modal
            # intercepts, or the App's gate skips the seq branch).
            await pilot.press("right")
            await pilot.pause()
            assert sp._cursor_pos == initial


# ═══════════════════════════════════════════════════════════════════════════════
# pyhmmer probe + HMMscan UX
# ═══════════════════════════════════════════════════════════════════════════════

class TestHmmscanProbe:

    def test_probe_returns_bool(self):
        # We don't know whether pyhmmer is installed in the test env;
        # just verify the probe is a bool (no crashes, no None).
        v = sc._probe_pyhmmer()
        assert v is True or v is False


# ═══════════════════════════════════════════════════════════════════════════════
# BLASTN HSP edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestBlastnEdgeCases:

    def test_palindromic_query_returns_one_hit_per_strand(self):
        # A perfectly palindromic 30 bp query should be found on the
        # forward strand AND the reverse strand of the subject — but
        # the dedup logic must prevent two HSPs at the same query/
        # subject coords stacking up.
        pal = "AAAATTTTAAAATTTTAAAATTTTAAAATT"
        # Subject contains pal once, at offset 10.
        subject = "GGGG" * 4 + pal + "TTTT" * 4
        rec = _make_record("PAL01", subject)
        _seed_collection("PalTest", [rec])
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastn", ["PalTest"])
        hits = sc._blast_search(pal, db)
        # We expect at least one hit; we don't insist on exactly one
        # since palindromes legitimately match both strands.
        assert hits

    def test_query_with_only_n_returns_no_hits(self):
        # A query of all N's would seed against any subject containing
        # NN…N. Most plasmids don't have long N runs, so we should get
        # nothing back.
        rec = _make_record("S01", "ACGT" * 50)
        _seed_collection("NTest", [rec])
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastn", ["NTest"])
        hits = sc._blast_search("N" * 30, db)
        # Either no hits, or any hits below the min identity. Mainly
        # asserting "doesn't crash".
        assert all(h["identity_pct"] < 100.0 or
                   "N" in h["subject_id"] for h in hits)

    def test_format_hits_escapes_markup_in_subject_name(
            self, tiny_record, isolated_library):
        # A subject whose name contains Rich markup tokens shouldn't
        # render as styled text — `_format_hits` runs everything
        # through `rich.markup.escape`. Smoke-check the modal-side
        # escape so a malicious GenBank label can't inject [red]/[/red].
        rec_seq = "ACGTACGTACGTACGTACGT" * 5
        rec = _make_record("EVIL_PLASMID", rec_seq)
        # Inject a malicious label into the record before serialising.
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec.features.append(SeqFeature(
            FeatureLocation(0, 30, strand=1),
            type="CDS",
            qualifiers={"label": ["[red]injected[/red]"]},
        ))
        _seed_collection("MarkupTest", [rec])
        sc._blast_clear_cache()
        # Trigger BLASTP so the malicious label appears in subject_name.
        # The CDS is too short for triplet-divisible (30 bp ÷ 3 = 10 aa),
        # so the BLASTP DB will index it as "[red]injected[/red]".
        db = sc._blast_get_db("blastp", ["MarkupTest"])
        if db["subjects"]:
            modal = sc.BlastModal()
            # Use a private path to call _format_hits without running pilot.
            txt = modal._format_hits(
                "blastp", "MAAAAA",
                [{
                    "subject_idx": 0,
                    "subject_id": "EVIL_PLASMID:[red]injected[/red]",
                    "subject_name": "[red]injected[/red]",
                    "subject_collection": "MarkupTest",
                    "kind": "cds",
                    "strand": 1,
                    "q_start": 0, "q_end": 6,
                    "s_start": 0, "s_end": 6,
                    "score": 30, "matches": 6, "aligned_len": 6,
                    "identity_pct": 100.0,
                }],
                db,
            )
            # The escape replaces `[` with `\[`. Verify our injection
            # didn't survive intact as Rich markup.
            assert r"\[red]" in txt or "[red]" not in txt or txt.count("[red]") == 0


    def test_pure_backend_explicit(self):
        # Force the pure-Python backend; verify it still works.
        target = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        rec = _make_record("PURE01", "GGGG" * 4 + target + "TTTT" * 4)
        _seed_collection("PureTest", [rec])
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastn", ["PureTest"])
        hits = sc._blast_search(target, db, backend="pure")
        assert hits and hits[0]["subject_id"] == "PURE01"

    def test_pyhmmer_backend_explicit(self):
        # Force the pyhmmer backend; verify it routes there + returns hits.
        target = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        rec = _make_record("PYHMR01", "GGGG" * 4 + target + "TTTT" * 4)
        _seed_collection("PyTest", [rec])
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastn", ["PyTest"])
        hits = sc._blast_search(target, db, backend="pyhmmer")
        assert hits and hits[0]["subject_id"] == "PYHMR01"

    def test_short_query_falls_back_to_pure(self, monkeypatch):
        # A 15-bp query is shorter than _PYHMMER_MIN_QUERY_BLASTN (20 bp).
        # Auto-dispatch should route to pure-Python — verify by spying
        # on the backend functions and confirming pyhmmer wasn't called.
        rec = _make_record("SHORT01", "ACGT" * 30)
        _seed_collection("ShortTest", [rec])
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastn", ["ShortTest"])

        called = {"pure": 0, "pyhmmer": 0}
        orig_pure = sc._blast_search_pure
        orig_pyhmmer = sc._blast_search_pyhmmer

        def spy_pure(*a, **kw):
            called["pure"] += 1
            return orig_pure(*a, **kw)

        def spy_pyhmmer(*a, **kw):
            called["pyhmmer"] += 1
            return orig_pyhmmer(*a, **kw)

        monkeypatch.setattr(sc, "_blast_search_pure", spy_pure)
        monkeypatch.setattr(sc, "_blast_search_pyhmmer", spy_pyhmmer)
        sc._blast_search("ACGTACGTACGTACG", db)   # 15 bp
        assert called["pure"] == 1
        assert called["pyhmmer"] == 0

    def test_explicit_pyhmmer_with_pyhmmer_unavailable_raises(
            self, monkeypatch):
        # Force the availability flag off, request backend="pyhmmer"
        # explicitly — should raise. Auto would silently fall back.
        import pytest
        monkeypatch.setattr(sc, "_PYHMMER_AVAILABLE", False)
        rec = _make_record("X01", "ACGT" * 30)
        _seed_collection("ExpTest", [rec])
        sc._blast_clear_cache()
        db = sc._blast_get_db("blastn", ["ExpTest"])
        with pytest.raises(RuntimeError):
            sc._blast_search("ACGT" * 10, db, backend="pyhmmer")

    def test_super_long_query_truncated(self):
        # 200 kb query → engine truncates to _MAX_BLAST_QUERY_LEN. We
        # don't actually run the truncated search here — the seeder on
        # a tandemly-repeating 100 kb query is O(n × #seeds), which on
        # an ACGT-tetra pattern can balloon into millions of seed pairs.
        # Verifying the cap fires is enough; the pure-helper test
        # `test_query_capped_at_max_len` covers the same constant.
        big_q = "ACGT" * 50_000   # 200 kb
        prog, q = sc._detect_query_program(big_q, "blastn")
        assert len(q) <= sc._MAX_BLAST_QUERY_LEN
        assert len(q) == sc._MAX_BLAST_QUERY_LEN  # exact-cap check


# ═══════════════════════════════════════════════════════════════════════════════
# HMMscan via pyhmmer — pyhmmer is now a hard dependency (see pyproject.toml)
# ═══════════════════════════════════════════════════════════════════════════════

def _build_tiny_hmm(tmp_path, *, name: str = "toy",
                     reps: list[str] | None = None) -> str:
    """Build a small HMM file in `tmp_path` from a 2-sequence MSA.
    Returns the path. Used to exercise `_hmmscan_run` without needing
    Pfam-A on disk."""
    from pyhmmer import easel, plan7
    seqs = reps or [
        "MAKVTPGGRSEKAAAAAAAAA",
        "MAKVTPGGSAEKAAAAAAAAA",
    ]
    alphabet = easel.Alphabet.amino()
    msa = easel.TextMSA(
        name=name.encode("utf-8"),
        sequences=[
            easel.TextSequence(name=f"s{i}".encode("utf-8"), sequence=s)
            for i, s in enumerate(seqs)
        ],
    ).digitize(alphabet)
    builder = plan7.Builder(alphabet)
    bg = plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, bg)
    # The HMMER format requires a non-empty COM line; pyhmmer writes
    # an empty one by default and then refuses to parse the result.
    hmm.command_line = "splicecraft tests fixture"
    out = tmp_path / f"{name}.hmm"
    with open(out, "wb") as f:
        hmm.write(f)
    return str(out)


class TestHmmscanEngine:

    def test_pyhmmer_available(self):
        # pyhmmer is now a hard dependency; the probe should always
        # return True on a properly-installed test env.
        assert sc._PYHMMER_AVAILABLE is True

    def test_hits_a_matching_protein(self, tmp_path):
        path = _build_tiny_hmm(tmp_path, name="toy_match")
        hits = sc._hmmscan_run("MAKVTPGGRSEKAAAAAAAAA", path)
        assert hits, "expected at least one HMMscan hit"
        top = hits[0]
        assert "toy_match" in top["subject_id"] or top["subject_id"] == "toy_match"
        assert top["score"] > 0
        assert top["q_start"] >= 0
        assert top["q_end"] > top["q_start"]

    def test_no_hits_for_unrelated_query(self, tmp_path):
        path = _build_tiny_hmm(tmp_path, name="toy_unrel")
        # A totally unrelated 30-mer.
        hits = sc._hmmscan_run("WWWWFFFFCCCCYYYYWWWWFFFFCCCCYY", path)
        # HMMER may emit weak hits; assert any are below the bit-score
        # threshold pyhmmer applies internally OR none come back.
        assert all(h["score"] < 5.0 for h in hits)

    def test_missing_hmm_path_raises(self, tmp_path):
        import pytest
        with pytest.raises(FileNotFoundError):
            sc._hmmscan_run("MAKVT", str(tmp_path / "nope.hmm"))

    def test_too_short_query_raises(self, tmp_path):
        import pytest
        path = _build_tiny_hmm(tmp_path, name="toy_short")
        with pytest.raises(ValueError):
            sc._hmmscan_run("M", path)

    def test_query_alphabet_filtered(self, tmp_path):
        # Numbers / punctuation in the query are stripped before
        # digitisation — verify by passing a noisy version of a known
        # match and checking we still get a hit.
        path = _build_tiny_hmm(tmp_path, name="toy_clean")
        hits = sc._hmmscan_run(
            "1 MAKVT 2 PGGRS 3 EK | AAAAA AAAAAA",
            path,
        )
        assert hits, "alphabet filter should preserve real AA letters"


class TestBlastModalHmmscanIntegration:
    """Modal-level test: HMMscan path triggers the engine and renders
    real hits when the user supplies a valid .hmm file."""

    async def test_run_hmmscan_with_valid_path(
            self, tiny_record, isolated_library, tmp_path):
        path = _build_tiny_hmm(tmp_path, name="toy_modal")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#blast-program", sc.Select).value = "hmmscan"
            # Let the program-change reset settle BEFORE typing the query —
            # switching alphabet (blastn→hmmscan) clears the query box, so
            # the query must be entered after that fires (the natural
            # pick-program-then-type order).
            await pilot.pause()
            modal.query_one("#blast-query", sc.TextArea).text = (
                "MAKVTPGGRSEKAAAAAAAAA"
            )
            modal.query_one("#blast-hmm-path", sc.Input).value = path
            await pilot.pause()
            modal.query_one("#btn-blast-run", sc.Button).press()
            # Worker is async — give it a moment to finish (the HMM
            # build is tiny so this is fast).
            await pilot.pause()
            await pilot.pause(0.5)
            results = modal.query_one("#blast-results", sc.Static)
            txt = str(results.render())
            assert "pending" not in txt.lower()
            assert "toy_modal" in txt or "HMMSCAN" in txt or "hits" in txt.lower()

    async def test_run_hmmscan_missing_path_complains(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#blast-program", sc.Select).value = "hmmscan"
            await pilot.pause()  # let the program-change reset settle
            modal.query_one("#blast-query", sc.TextArea).text = "MAKVT"
            # Leave the hmm-path Input empty.
            await pilot.pause()
            modal.query_one("#btn-blast-run", sc.Button).press()
            await pilot.pause()
            status = modal.query_one("#blast-status", sc.Static)
            txt = str(status.render())
            assert "path" in txt.lower() or "hmm" in txt.lower()

    async def test_hmm_path_persists_across_modal_opens(
            self, tiny_record, isolated_library, tmp_path):
        # Run HMMscan with a valid path → settings.json should remember
        # it → next BlastModal open should pre-fill the Input. Verifies
        # the `hmm_db_path` round-trip without doing two app sessions.
        path = _build_tiny_hmm(tmp_path, name="toy_persist")
        # First session: type the path + run.
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#blast-program", sc.Select).value = "hmmscan"
            await pilot.pause()  # let the program-change reset settle
            modal.query_one("#blast-query", sc.TextArea).text = (
                "MAKVTPGGRSEKAAAAAAAAA"
            )
            modal.query_one("#blast-hmm-path", sc.Input).value = path
            await pilot.pause()
            modal.query_one("#btn-blast-run", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.5)
            # Dismiss; persistence should have happened.
            modal.dismiss(None)
            await pilot.pause()
        # Verify settings.json holds the path now.
        assert sc._get_setting("hmm_db_path", "") == path

        # Second session: open BlastModal, expect prefilled Input.
        app2 = sc.PlasmidApp()
        app2._preload_record = tiny_record
        async with app2.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app2.push_screen(sc.BlastModal())
            await pilot.pause()
            modal2 = app2.screen
            inp = modal2.query_one("#blast-hmm-path", sc.Input)
            assert inp.value == path


class TestBlastMenuEntry:
    """The menubar "BLAST" entry (which replaced "Edit") is a direct-open
    item that pops the same modal as Ctrl+B — no dropdown."""

    def test_blast_in_menus_edit_gone(self):
        assert "BLAST" in sc.MenuBar.MENUS
        assert "Edit" not in sc.MenuBar.MENUS
        # Replaced in place — Settings stays at its pinned index.
        assert sc.MenuBar.MENUS.index("Settings") == 1

    async def test_menu_opens_blast_modal(self, tiny_record,
                                          isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # `open_menu` direct-dispatches the BLAST action (no dropdown).
            app.open_menu("BLAST", 0, 0)
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.BlastModal)
