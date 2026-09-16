"""test_feature_presets — the built-in feature-preset catalogue.

The catalogue (`splicecraft_presets._FEATURE_PRESETS`) is SHIPPED DATA that
users annotate real plasmids with, so a silent edit to one base is a silent
edit to somebody's map. These tests treat the data itself as the thing under
test:

  * **Data integrity** — shape, unique names/sequences, ACGT-only bases, known
    feature types + categories, valid colours, and (for coding entries) a real
    reading frame: divisible by three, recognised start, terminal stop, no
    internal stop. Cross-validated against Biopython's translator, not against
    SpliceCraft's own.
  * **Provenance** — every entry cites a GenBank accession. An entry with no
    source is an entry nobody can check.
  * **Drift guard** — an aggregate digest over the catalogue. Editing a preset
    on purpose means updating one constant; editing one by accident fails.
  * **Presets never reach disk** — the merged library+preset view is read-only.
    `features.json` must only ever hold what the user put there.
  * **Merge semantics** — a user entry shadows a preset of the same
    (name, feature_type) instead of doubling the row.
  * **Annotation scan** — opt-in via `include_presets`, both strands, plus the
    2026-09-15 wrap-tail duplicate regression.
  * **Agent endpoints** — list / get / import, happy and error paths.
  * **UI** — the browser modal and the Synthesis side panel.
"""
from __future__ import annotations

import copy
import hashlib

import pytest

import splicecraft as sc
import splicecraft_presets as P


_TERM = (160, 48)

# Aggregate digest over the whole catalogue: name, feature_type, strand and
# sequence of every entry, in catalogue order. A deliberate data change updates
# this constant in the same commit; an accidental one fails loudly here.
# Descriptions and colours are deliberately OUT of the digest — prose and
# palette tweaks shouldn't need a test edit, but a base never changes silently.
_CATALOGUE_DIGEST = "8aebb5ffd3b84ec364c459474673c64167e32d79652ac912005296b42450568e"
_CATALOGUE_COUNT = 116


def _catalogue_digest() -> str:
    h = hashlib.sha256()
    for e in P._FEATURE_PRESETS:
        h.update(str(e["name"]).encode())
        h.update(b"\x00")
        h.update(str(e["feature_type"]).encode())
        h.update(b"\x00")
        h.update(str(e["strand"]).encode())
        h.update(b"\x00")
        h.update(str(e["sequence"]).encode())
        h.update(b"\x1e")
    return h.hexdigest()


# ═══════════════════════════════════════════════════════════════════════════════
# Data integrity
# ═══════════════════════════════════════════════════════════════════════════════

class TestCatalogueIntegrity:

    def test_catalogue_is_non_empty_and_sized(self):
        assert len(P._FEATURE_PRESETS) == _CATALOGUE_COUNT

    def test_digest_matches(self):
        """Drift guard. If this fails and the change was deliberate, update
        `_CATALOGUE_DIGEST` + `_CATALOGUE_COUNT` in the same commit and say
        in the message which preset changed and why."""
        assert _catalogue_digest() == _CATALOGUE_DIGEST, (
            "feature-preset catalogue changed; update _CATALOGUE_DIGEST "
            "deliberately or revert the edit"
        )

    def test_required_fields_present_and_typed(self):
        for e in P._FEATURE_PRESETS:
            for key in ("name", "feature_type", "category", "strand",
                        "color", "sequence", "description", "source",
                        "aliases"):
                assert key in e, f"{e.get('name')!r} missing {key!r}"
            assert isinstance(e["name"], str) and e["name"].strip()
            assert isinstance(e["description"], str) and e["description"].strip()
            assert isinstance(e["source"], str) and e["source"].strip()
            assert isinstance(e["aliases"], list)
            assert all(isinstance(a, str) and a.strip() for a in e["aliases"])

    def test_names_are_unique(self):
        names = [e["name"] for e in P._FEATURE_PRESETS]
        dupes = {n for n in names if names.count(n) > 1}
        assert not dupes, f"duplicate preset names: {sorted(dupes)}"

    def test_sequences_are_unique(self):
        """Two presets with identical bases would annotate the same span
        twice under two names — noise on every map."""
        seqs: dict[str, str] = {}
        for e in P._FEATURE_PRESETS:
            prev = seqs.get(e["sequence"])
            assert prev is None, (
                f"{e['name']!r} has the same sequence as {prev!r}"
            )
            seqs[e["sequence"]] = e["name"]

    def test_sequences_are_plain_acgt(self):
        """Presets come out of real records, so ambiguity codes would mean a
        transcription error rather than a genuine degenerate element."""
        for e in P._FEATURE_PRESETS:
            bad = sorted(set(e["sequence"]) - set("ACGT"))
            assert not bad, f"{e['name']!r} has non-ACGT bases: {bad}"

    def test_sequence_lengths_are_sane(self):
        for e in P._FEATURE_PRESETS:
            n = len(e["sequence"])
            assert n >= 8, f"{e['name']!r} is {n} bp — too short to be useful"
            assert n <= sc._MAX_FEATURE_SEQ_LEN, (
                f"{e['name']!r} is {n} bp, over the feature-library cap"
            )

    def test_feature_types_are_genbank_vocabulary(self):
        """A type outside `_GENBANK_FEATURE_TYPES` would not round-trip
        through AddFeatureModal's Select or a GenBank export."""
        for e in P._FEATURE_PRESETS:
            assert e["feature_type"] in sc._GENBANK_FEATURE_TYPES, (
                f"{e['name']!r} has unknown feature_type "
                f"{e['feature_type']!r}"
            )

    def test_categories_are_declared(self):
        for e in P._FEATURE_PRESETS:
            assert e["category"] in P._PRESET_CATEGORIES, (
                f"{e['name']!r} has undeclared category {e['category']!r}"
            )

    def test_every_declared_category_is_used(self):
        used = {e["category"] for e in P._FEATURE_PRESETS}
        unused = set(P._PRESET_CATEGORIES) - used
        assert not unused, f"declared but unused categories: {sorted(unused)}"

    def test_colors_are_six_digit_hex(self):
        for e in P._FEATURE_PRESETS:
            col = e["color"]
            assert isinstance(col, str) and sc._HEX6_RE.fullmatch(col), (
                f"{e['name']!r} has non-hex colour {col!r}"
            )

    def test_strand_is_forward(self):
        """Sequences are stored in their own functional orientation (the
        extractor reverse-complemented minus-strand features), so a preset is
        always a forward-reading element."""
        for e in P._FEATURE_PRESETS:
            assert e["strand"] == 1, f"{e['name']!r} strand={e['strand']}"

    def test_every_entry_cites_a_record(self):
        for e in P._FEATURE_PRESETS:
            assert "GenBank" in e["source"], (
                f"{e['name']!r} has no GenBank citation: {e['source']!r}"
            )

    def test_coding_entries_have_a_real_reading_frame(self):
        """Cross-validated against Biopython, deliberately: an in-house
        translator agreeing with itself proves nothing."""
        from Bio.Seq import Seq
        for e in P._FEATURE_PRESETS:
            if e["feature_type"] != "CDS":
                continue
            seq = e["sequence"]
            name = e["name"]
            assert len(seq) % 3 == 0, f"{name!r} is {len(seq)} bp, not /3"
            assert seq[:3] in ("ATG", "GTG", "TTG"), (
                f"{name!r} starts with {seq[:3]!r}"
            )
            prot = str(Seq(seq).translate())
            assert prot.endswith("*"), f"{name!r} has no terminal stop"
            assert prot[:-1].count("*") == 0, (
                f"{name!r} has an internal stop codon"
            )

    def test_tag_entries_carry_no_stop_codon(self):
        """A fusion tag with an in-frame stop silently truncates the protein
        it is fused to. Regression guard for the 2026-09-15 curation pass,
        which caught the GST entry running on through the source vector's
        thrombin site, FLAG tag, polylinker AND stop codon."""
        from Bio.Seq import Seq
        for e in P._FEATURE_PRESETS:
            if e["category"] not in ("Tag", "Protease site", "Linker / 2A",
                                     "Localisation"):
                continue
            seq = e["sequence"]
            if len(seq) % 3:
                continue          # not a whole number of codons; not a fusion
            prot = str(Seq(seq).translate())
            assert "*" not in prot, (
                f"{e['name']!r} contains a stop codon — it would truncate "
                f"any fusion it is placed in front of"
            )

    def test_only_the_intended_presets_nest_inside_each_other(self):
        """One preset's bases sitting inside another's means both annotate
        the same span. Four such pairs are DELIBERATE — the user wants the
        part and the whole (a bare CMV promoter or the enhancer+promoter
        block; one rrnB terminator or the tandem pair). Any OTHER nesting is
        an accident: two names for one element, or a curation slip that
        pasted vector context into an entry. The allow-list is the point —
        it makes a new nesting fail instead of quietly doubling a map."""
        allowed = {
            ("CMV promoter", "CMV enhancer+promoter"),
            ("tet operator (tetO2)", "CMV/TetO2 promoter"),
            ("GAL1 promoter", "GAL1/GAL10 bidirectional promoter"),
            ("rrnB T1 terminator", "rrnB T1+T2 terminator"),
        }
        found = set()
        for a in P._FEATURE_PRESETS:
            for b in P._FEATURE_PRESETS:
                if a is b:
                    continue
                if b["sequence"] in a["sequence"]:
                    found.add((b["name"], a["name"]))
        assert found == allowed, (
            f"unexpected nesting: {sorted(found - allowed)}; "
            f"no longer nesting: {sorted(allowed - found)}"
        )

    def test_only_the_intended_presets_nest_reverse_complemented(self):
        """Same check in the other direction. These three are real biology,
        not duplication: in pUC19 the polylinker and the M13 forward primer
        site genuinely sit inside lacZ-alpha, which reads the other way, and
        the SV40 early and late polyadenylation signals share a region on
        opposite strands. They annotate DIFFERENT spans, so all of them are
        worth drawing — but a NEW reverse-complement nesting is a curation
        slip and should fail here."""
        allowed = {
            ("SV40 poly(A) signal", "SV40 late poly(A) signal"),
            ("M13 forward primer site", "lacZ-alpha"),
            ("pUC19 MCS (EcoRI-HindIII)", "lacZ-alpha"),
        }
        found = set()
        for a in P._FEATURE_PRESETS:
            for b in P._FEATURE_PRESETS:
                if a is b:
                    continue
                if sc._rc(b["sequence"]) in a["sequence"]:
                    found.add((b["name"], a["name"]))
        assert found == allowed, (
            f"unexpected rc-nesting: {sorted(found - allowed)}; "
            f"no longer nesting: {sorted(allowed - found)}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Accessors
# ═══════════════════════════════════════════════════════════════════════════════

class TestPresetAccessors:

    def test_preset_features_returns_a_deep_copy(self):
        a = P._preset_features()
        a[0]["name"] = "clobbered"
        a[0]["aliases"].append("clobbered")
        b = P._preset_features()
        assert b[0]["name"] != "clobbered"
        assert "clobbered" not in b[0]["aliases"]

    def test_categories_are_in_declared_order(self):
        got = P._preset_categories()
        assert list(got) == [c for c in P._PRESET_CATEGORIES if c in got]

    def test_find_preset_is_case_insensitive(self):
        assert P._find_preset("t7 promoter")["name"] == "T7 promoter"
        assert P._find_preset("  T7 Promoter  ")["name"] == "T7 promoter"

    def test_find_preset_returns_a_copy(self):
        got = P._find_preset("loxP")
        got["sequence"] = "AAAA"
        assert P._find_preset("loxP")["sequence"] != "AAAA"

    def test_find_preset_misses_return_none(self):
        assert P._find_preset("no such element") is None
        assert P._find_preset("") is None
        assert P._find_preset(None) is None          # type: ignore[arg-type]
        assert P._find_preset("loxP", "CDS") is None  # wrong feature_type

    def test_preset_matches_searches_aliases_and_category(self):
        his = P._find_preset("6xHis tag")
        assert P._preset_matches(his, "ni-nta")     # alias
        assert P._preset_matches(his, "TAG")        # category, case-folded
        assert P._preset_matches(his, "")           # empty matches everything
        assert not P._preset_matches(his, "zzzzz")

    def test_preset_matches_rejects_non_dicts(self):
        assert P._preset_matches(None, "x") is False   # type: ignore[arg-type]
        assert P._preset_matches("nope", "x") is False  # type: ignore[arg-type]


class TestPresetToLibraryEntry:

    def test_shape_matches_a_library_entry(self):
        entry = P._preset_to_library_entry(P._find_preset("T7 promoter"))
        assert set(entry) == {
            "name", "feature_type", "strand", "color", "sequence",
            "qualifiers", "description",
        }
        assert entry["name"] == "T7 promoter"
        assert entry["sequence"].isupper()

    def test_preset_only_fields_are_dropped(self):
        entry = P._preset_to_library_entry(P._find_preset("loxP"))
        for gone in ("category", "aliases", "source", "preset"):
            assert gone not in entry

    def test_provenance_lands_in_qualifiers(self):
        """The accession has to travel with the annotation, or a plasmid
        exported from SpliceCraft tells the next reader nothing about where
        its features came from."""
        preset = P._find_preset("AmpR (bla)")
        entry = P._preset_to_library_entry(preset)
        note = entry["qualifiers"]["note"][0]
        assert "SpliceCraft preset" in note
        assert "GenBank" in note
        assert entry["qualifiers"]["label"] == ["AmpR (bla)"]

    def test_rejects_non_dict(self):
        with pytest.raises(TypeError):
            P._preset_to_library_entry("nope")   # type: ignore[arg-type]

    def test_round_trips_through_the_real_save_path(self):
        entry = P._preset_to_library_entry(P._find_preset("FRT"))
        sc._save_features([entry])
        back = sc._load_features()
        assert len(back) == 1
        assert back[0]["sequence"] == entry["sequence"]
        assert back[0]["qualifiers"]["note"] == entry["qualifiers"]["note"]


# ═══════════════════════════════════════════════════════════════════════════════
# Merge semantics — and the safety property that presets never reach disk
# ═══════════════════════════════════════════════════════════════════════════════

class TestMergeWithLibrary:

    def test_empty_library_yields_the_whole_catalogue(self):
        sc._save_features([])
        merged = sc._load_features_with_presets()
        assert len(merged) == len(P._FEATURE_PRESETS)
        assert all(e.get("preset") for e in merged)

    def test_user_entries_come_first_and_keep_their_order(self):
        sc._save_features([
            {"name": "zzz", "feature_type": "CDS", "sequence": "ATGAAA",
             "strand": 1, "color": ""},
            {"name": "aaa", "feature_type": "CDS", "sequence": "ATGCCC",
             "strand": 1, "color": ""},
        ])
        merged = sc._load_features_with_presets()
        assert [e["name"] for e in merged[:2]] == ["zzz", "aaa"]
        assert not merged[0].get("preset")
        assert not merged[1].get("preset")

    def test_user_entry_shadows_a_preset_of_the_same_key(self):
        preset = P._find_preset("loxP")
        sc._save_features([{
            "name": "loxP", "feature_type": preset["feature_type"],
            "sequence": "ACGTACGTACGTACGT", "strand": 1, "color": "",
        }])
        merged = sc._load_features_with_presets()
        rows = [e for e in merged if e["name"] == "loxP"]
        assert len(rows) == 1, "shadowed preset should vanish, not double up"
        assert rows[0]["sequence"] == "ACGTACGTACGTACGT"
        assert not rows[0].get("preset")

    def test_same_name_different_type_does_not_shadow(self):
        preset = P._find_preset("loxP")
        assert preset["feature_type"] != "CDS"
        sc._save_features([{
            "name": "loxP", "feature_type": "CDS", "sequence": "ATGAAATAA",
            "strand": 1, "color": "",
        }])
        rows = [e for e in sc._load_features_with_presets()
                if e["name"] == "loxP"]
        assert len(rows) == 2

    def test_merged_view_is_safe_to_mutate(self):
        merged = sc._load_features_with_presets()
        for e in merged:
            e["sequence"] = "MUTATED"
        assert P._FEATURE_PRESETS[0]["sequence"] != "MUTATED"
        assert P._preset_features()[0]["sequence"] != "MUTATED"

    def test_non_dict_user_rows_are_skipped(self):
        merged = P._merge_presets_with_library(["junk", None, 42])
        assert len(merged) == len(P._FEATURE_PRESETS)

    def test_none_library_is_accepted(self):
        assert len(P._merge_presets_with_library(None)) == len(P._FEATURE_PRESETS)

    def test_presets_never_reach_the_user_data_file(self):
        """The safety property. Browsing the catalogue must not be able to
        bake 116 shipped entries into somebody's `features.json` — presets
        live in code, and `_load_features` is the user's own library."""
        sc._save_features([
            {"name": "mine", "feature_type": "CDS", "sequence": "ATGAAATAA",
             "strand": 1, "color": ""},
        ])
        # Build the merged view the way every consumer does...
        merged = sc._load_features_with_presets()
        assert len(merged) > 1
        # ...and the user's own library is untouched by it.
        after = sc._load_features()
        assert [e["name"] for e in after] == ["mine"]
        assert not any(e.get("preset") for e in after)


# ═══════════════════════════════════════════════════════════════════════════════
# Annotation scan
# ═══════════════════════════════════════════════════════════════════════════════

class TestAnnotateWithPresets:

    def _template(self, *names: str) -> tuple[str, list[tuple[str, int, int]]]:
        """Concatenate presets with 40 bp of filler between them; return the
        sequence and the coordinates each one should be found at."""
        filler = "TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTC"
        seq = filler
        spans: list[tuple[str, int, int]] = []
        for n in names:
            p = P._find_preset(n)
            assert p is not None, n
            spans.append((n, len(seq), len(seq) + len(p["sequence"])))
            seq += p["sequence"] + filler
        return seq, spans

    def test_presets_are_opt_in(self):
        sc._save_features([])
        seq, _ = self._template("T7 promoter", "AmpR (bla)")
        assert sc._annotate_seq_from_feature_library(seq) == []

    def test_opt_in_finds_presets_at_exact_coordinates(self):
        sc._save_features([])
        seq, spans = self._template("T7 promoter", "AmpR (bla)", "loxP")
        hits = sc._annotate_seq_from_feature_library(seq, include_presets=True)
        got = {(h["name"], h["start"], h["end"]) for h in hits}
        for name, start, end in spans:
            assert (name, start, end) in got, f"{name} not found at {start}"

    def test_reverse_complemented_preset_is_found_on_the_minus_strand(self):
        sc._save_features([])
        amp = P._find_preset("AmpR (bla)")["sequence"]
        seq = "GGGGGGGGGGGG" + sc._rc(amp) + "CCCCCCCCCCCC"
        hits = sc._annotate_seq_from_feature_library(seq, include_presets=True)
        amp_hits = [h for h in hits if h["name"] == "AmpR (bla)"]
        assert len(amp_hits) == 1
        assert amp_hits[0]["strand"] == -1
        assert amp_hits[0]["start"] == 12

    def test_user_entry_shadows_preset_in_the_scan(self):
        """A user who renamed a preset must see THEIR name on the map."""
        preset = P._find_preset("T7 promoter")
        sc._save_features([{
            "name": "my T7", "feature_type": preset["feature_type"],
            "sequence": preset["sequence"], "strand": 1, "color": "",
        }])
        seq, _ = self._template("T7 promoter")
        hits = sc._annotate_seq_from_feature_library(seq, include_presets=True)
        names = {h["name"] for h in hits}
        assert names == {"my T7"}, (
            "the user's own name must win: the preset carries the same bases "
            "under a different (name, feature_type) key, so it survives the "
            "merge, but the scan must not stack a second band on the span"
        )
        spans = [(h["start"], h["end"], h["strand"]) for h in hits]
        assert len(spans) == len(set(spans))

    def test_circular_scan_does_not_duplicate_hits(self):
        """Regression guard for the 2026-09-15 wrap-tail bug: the circular
        scan appends `seq[:longest-1]`, and a hit found inside that appended
        copy was emitted a SECOND time with an `end` one full turn too large
        (it escaped dedup because only `start` was reduced modulo the
        length). With a 1.6 kb preset in the catalogue the tail reaches most
        of a normal plasmid, so every annotation landed twice."""
        sc._save_features([])
        seq, spans = self._template("T7 promoter", "loxP", "FRT")
        hits = sc._annotate_seq_from_feature_library(
            seq, circular=True, include_presets=True,
        )
        keys = [(h["name"], h["start"], h["strand"]) for h in hits]
        assert len(keys) == len(set(keys)), f"duplicate hits: {keys}"
        for h in hits:
            assert h["end"] <= len(seq), (
                f"{h['name']} ends at {h['end']} past the {len(seq)} bp "
                f"template without wrapping the origin"
            )

    def test_genuine_origin_wrap_still_reports_end_past_total(self):
        """The duplicate guard must not cost us real wrap detection."""
        sc._save_features([])
        lox = P._find_preset("loxP")["sequence"]
        head, tail = lox[:20], lox[20:]
        seq = tail + "ACGT" * 20 + head        # loxP straddles bp 0
        hits = sc._annotate_seq_from_feature_library(
            seq, circular=True, include_presets=True,
        )
        wrapped = [h for h in hits
                   if h["name"] == "loxP" and h["end"] > len(seq)]
        assert wrapped, f"wrap hit missing; got {hits}"
        assert wrapped[0]["start"] == len(seq) - 20

    def test_scan_output_carries_preset_provenance(self):
        sc._save_features([])
        seq, _ = self._template("loxP")
        hits = sc._annotate_seq_from_feature_library(seq, include_presets=True)
        assert hits
        assert hits[0]["description"]
        assert hits[0]["color"]


class TestAnnotateRealPlasmid:
    """End-to-end on a plasmid assembled from real elements — the acceptance
    test that the catalogue actually annotates something a user would load."""

    def _pseudo_puc19(self) -> str:
        parts = ["ori (pUC/pMB1, high-copy)", "AmpR promoter", "AmpR (bla)",
                 "lac promoter", "lacZ-alpha", "pUC19 MCS (EcoRI-HindIII)"]
        spacer = "GATCGATCGATCGATCGATC"
        return spacer.join(P._find_preset(n)["sequence"] for n in parts)

    def test_all_six_elements_are_recovered(self):
        sc._save_features([])
        seq = self._pseudo_puc19()
        hits = sc._annotate_seq_from_feature_library(
            seq, circular=True, include_presets=True,
        )
        names = {h["name"] for h in hits}
        for expected in ("ori (pUC/pMB1, high-copy)", "AmpR (bla)",
                         "lac promoter", "lacZ-alpha",
                         "pUC19 MCS (EcoRI-HindIII)"):
            assert expected in names, f"{expected} missing from {sorted(names)}"

    def test_hits_land_on_the_bases_they_claim(self):
        """Every reported span must actually contain the preset's bases —
        the check that catches an off-by-one in the coordinate mapping."""
        sc._save_features([])
        seq = self._pseudo_puc19()
        hits = sc._annotate_seq_from_feature_library(
            seq, circular=True, include_presets=True,
        )
        assert hits
        for h in hits:
            span = seq[h["start"]:h["end"]]
            if h["end"] > len(seq):
                span = seq[h["start"]:] + seq[:h["end"] - len(seq)]
            want = h["sequence"]
            if h["strand"] == -1:
                want = sc._rc(want)
            assert span == want, (
                f"{h['name']} at {h['start']}..{h['end']} does not match "
                f"its own sequence"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Agent API
# ═══════════════════════════════════════════════════════════════════════════════

def _call(name: str, payload: dict | None = None):
    fn, _write = sc._state._AGENT_HANDLERS[name]
    return fn(None, payload or {})


class TestPresetEndpoints:

    def test_list_returns_every_preset_without_sequences(self):
        r = _call("list-feature-presets")
        assert r["ok"] is True
        assert r["count"] == len(P._FEATURE_PRESETS)
        assert r["categories"] == list(P._preset_categories())
        assert "sequence" not in r["presets"][0]
        assert r["presets"][0]["sequence_length"] > 0

    def test_list_can_include_sequences(self):
        r = _call("list-feature-presets", {"include_sequence": True})
        assert all(p["sequence"] for p in r["presets"])

    def test_list_filters_by_category(self):
        r = _call("list-feature-presets", {"category": "Resistance"})
        assert r["count"] > 0
        assert {p["category"] for p in r["presets"]} == {"Resistance"}

    def test_list_filters_by_search_including_aliases(self):
        r = _call("list-feature-presets", {"search": "ni-nta"})
        assert [p["name"] for p in r["presets"]] == ["6xHis tag"]

    def test_list_rejects_unknown_category(self):
        body, code = _call("list-feature-presets", {"category": "Nope"})
        assert code == 400
        assert "unknown category" in body["error"]

    def test_list_rejects_wrong_types(self):
        assert _call("list-feature-presets", {"category": 7})[1] == 400
        assert _call("list-feature-presets", {"search": []})[1] == 400

    def test_get_returns_the_full_entry(self):
        r = _call("get-feature-preset", {"name": "loxP"})
        assert r["ok"] is True
        assert r["preset"]["sequence"] == P._find_preset("loxP")["sequence"]
        assert "GenBank" in r["preset"]["source"]

    def test_get_is_case_insensitive(self):
        assert _call("get-feature-preset", {"name": "LOXP"})["ok"] is True

    def test_get_missing_is_404(self):
        body, code = _call("get-feature-preset", {"name": "no such thing"})
        assert code == 404
        assert "no preset named" in body["error"]

    def test_get_requires_a_name(self):
        assert _call("get-feature-preset", {})[1] == 400
        assert _call("get-feature-preset", {"name": ""})[1] == 400
        assert _call("get-feature-preset", {"name": 5})[1] == 400

    def test_import_copies_into_the_user_library(self):
        sc._save_features([])
        r = _call("import-feature-preset", {"names": ["loxP", "FRT"]})
        assert r["ok"] is True and r["count"] == 2 and r["replaced"] == 0
        lib = sc._load_features()
        assert [e["name"] for e in lib] == ["loxP", "FRT"]
        assert "GenBank" in lib[0]["qualifiers"]["note"][0]

    def test_import_accepts_a_single_name(self):
        sc._save_features([])
        r = _call("import-feature-preset", {"name": "T7 promoter"})
        assert r["count"] == 1
        assert sc._load_features()[0]["name"] == "T7 promoter"

    def test_reimport_replaces_rather_than_duplicates(self):
        sc._save_features([])
        _call("import-feature-preset", {"name": "loxP"})
        r = _call("import-feature-preset", {"name": "loxP"})
        assert r["replaced"] == 1
        assert len(sc._load_features()) == 1

    def test_import_is_all_or_nothing_on_an_unknown_name(self):
        """A half-applied import is worse than a refused one: the caller
        cannot tell which half landed."""
        sc._save_features([])
        body, code = _call("import-feature-preset",
                           {"names": ["loxP", "definitely not a preset"]})
        assert code == 404
        assert "definitely not a preset" in body["error"]
        assert sc._load_features() == []

    def test_import_validates_its_payload(self):
        assert _call("import-feature-preset", {})[1] == 400
        assert _call("import-feature-preset", {"names": []})[1] == 400
        assert _call("import-feature-preset", {"names": "loxP"})[1] == 400
        assert _call("import-feature-preset", {"names": [1, 2]})[1] == 400
        assert _call("import-feature-preset",
                     {"names": ["loxP"] * 501})[1] == 400

    def test_import_is_registered_as_a_write_endpoint(self):
        assert sc._state._AGENT_HANDLERS["import-feature-preset"][1] is True
        assert sc._state._AGENT_HANDLERS["list-feature-presets"][1] is False
        assert sc._state._AGENT_HANDLERS["get-feature-preset"][1] is False

    def test_list_feature_library_still_shows_only_user_entries(self):
        """The user's library endpoint must not start returning 116 shipped
        rows — an agent diffing it would see a phantom import."""
        sc._save_features([])
        r = _call("list-feature-library")
        assert r["count"] == 0


# ═══════════════════════════════════════════════════════════════════════════════
# UI — preset browser modal
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.asyncio
class TestFeaturePresetsModal:

    async def _open(self, pilot, app):
        app.push_screen(sc.FeaturePresetsModal())
        await pilot.pause()
        await pilot.pause()
        return app.screen

    async def test_opens_and_lists_every_preset(self):
        sc._save_features([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._open(pilot, app)
            assert isinstance(scr, sc.FeaturePresetsModal)
            tbl = scr.query_one("#fpre-table")
            assert tbl.row_count == len(P._FEATURE_PRESETS)

    async def test_category_filter_narrows_the_table(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._open(pilot, app)
            scr.query_one("#fpre-cat").value = "Resistance"
            await pilot.pause()
            await pilot.pause()
            n = len([e for e in P._FEATURE_PRESETS
                     if e["category"] == "Resistance"])
            assert scr.query_one("#fpre-table").row_count == n

    async def test_search_matches_aliases(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._open(pilot, app)
            scr.query_one("#fpre-search").value = "ni-nta"
            await pilot.pause()
            await pilot.pause()
            assert scr.query_one("#fpre-table").row_count == 1

    async def test_toggle_selects_and_add_dismisses_with_entries(self):
        sc._save_features([])
        app = sc.PlasmidApp()
        got: list = []
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            app.push_screen(sc.FeaturePresetsModal(), callback=got.append)
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.action_toggle_selection()
            await pilot.pause()
            assert len(scr._selected) == 1
            scr.query_one("#btn-fpre-add").press()
            await pilot.pause()
            await pilot.pause()
        assert got and got[0]["action"] == "add"
        entry = got[0]["entries"][0]
        assert "qualifiers" in entry and "category" not in entry

    async def test_add_with_no_selection_uses_the_highlighted_row(self):
        app = sc.PlasmidApp()
        got: list = []
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            app.push_screen(sc.FeaturePresetsModal(), callback=got.append)
            await pilot.pause()
            await pilot.pause()
            app.screen.query_one("#btn-fpre-add").press()
            await pilot.pause()
            await pilot.pause()
        assert got and len(got[0]["entries"]) == 1

    async def test_close_dismisses_with_none(self):
        app = sc.PlasmidApp()
        got: list = []
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            app.push_screen(sc.FeaturePresetsModal(), callback=got.append)
            await pilot.pause()
            await pilot.pause()
            app.screen.query_one("#btn-fpre-close").press()
            await pilot.pause()
            await pilot.pause()
        assert got == [None]

    async def test_already_in_library_rows_are_flagged(self):
        preset = P._find_preset("loxP")
        sc._save_features([P._preset_to_library_entry(preset)])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._open(pilot, app)
            assert ("loxP", preset["feature_type"]) in scr._lib_keys

    async def test_modal_opts_out_of_app_undo(self):
        """[INV-41] — the modal hosts an Input, so Ctrl+Z must stay local."""
        assert sc.FeaturePresetsModal._blocks_undo is True


# ═══════════════════════════════════════════════════════════════════════════════
# UI — Feature Library screen + Synthesis side panel
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.asyncio
class TestFeatureLibraryScreenPresets:

    async def test_presets_button_opens_the_browser(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause()
            app.screen.action_presets()
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.FeaturePresetsModal)

    async def test_imported_presets_land_in_the_unsaved_buffer(self):
        """Not on disk: the screen owns one commit point (Ctrl+S), and a
        modal writing behind its back would be overwritten by the next save."""
        sc._save_features([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            entries = [P._preset_to_library_entry(P._find_preset(n))
                       for n in ("loxP", "FRT")]
            for e in entries:
                scr._upsert_entry(e, notice="Added preset", notify=False)
            await pilot.pause()
            assert [e["name"] for e in scr._entries] == ["loxP", "FRT"]
            assert scr._has_pending_changes is True
            assert sc._load_features() == []       # nothing written yet
            assert scr._persist_all() is True
            assert [e["name"] for e in sc._load_features()] == ["loxP", "FRT"]


@pytest.mark.asyncio
class TestSynthesisPanelPresets:

    async def _panel(self, pilot, app):
        app.push_screen(sc.SynthesisScreen())
        await pilot.pause()
        await pilot.pause()
        return app.screen

    async def test_panel_lists_presets_alongside_user_entries(self):
        sc._save_features([{
            "name": "my thing", "feature_type": "CDS",
            "sequence": "ATGAAATAA", "strand": 1, "color": "",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._panel(pilot, app)
            names = [e.get("name") for _k, e in scr._featlib_rows]
            assert "my thing" in names
            assert "T7 promoter" in names
            assert len(scr._featlib_rows) == len(P._FEATURE_PRESETS) + 1

    async def test_preset_rows_are_flagged_as_presets(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._panel(pilot, app)
            row = next(e for _k, e in scr._featlib_rows
                       if e.get("name") == "loxP")
            assert row.get("preset") is True

    async def test_filter_matches_preset_aliases(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._panel(pilot, app)
            scr._refresh_featlib_table(filter_str="ni-nta")
            await pilot.pause()
            assert [e["name"] for _k, e in scr._featlib_rows] == ["6xHis tag"]

    async def test_selected_preset_is_normalised_with_provenance(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._panel(pilot, app)
            scr._refresh_featlib_table(filter_str="loxP")
            await pilot.pause()
            row = next(i for i, (_k, e) in enumerate(scr._featlib_rows)
                       if e.get("name") == "loxP")
            scr.query_one("#syn-featlib-table").move_cursor(row=row)
            await pilot.pause()
            entry = scr._featlib_selected_entry()
            assert entry["name"] == "loxP"
            assert entry["preset"] is True
            assert "GenBank" in entry["qualifiers"]["note"][0]

    async def test_inserting_a_preset_annotates_with_its_provenance(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._panel(pilot, app)
            scr._refresh_featlib_table(filter_str="T7 promoter")
            await pilot.pause()
            row = next(i for i, (_k, e) in enumerate(scr._featlib_rows)
                       if e.get("name") == "T7 promoter")
            scr.query_one("#syn-featlib-table").move_cursor(row=row)
            await pilot.pause()
            scr._featlib_insert_selected(mode="insert")
            await pilot.pause()
            await pilot.pause()
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            assert ed._seq == P._find_preset("T7 promoter")["sequence"]
            assert len(ed._feats) == 1
            feat = ed._feats[0]
            assert feat["label"] == "T7 promoter"
            assert "GenBank" in feat["qualifiers"]["note"][0]

    async def test_editing_a_preset_leaves_the_catalogue_alone(self):
        """Copy-on-write: the Edit button on a preset row opens the modal
        but must not be able to mutate shipped data."""
        before = copy.deepcopy(P._FEATURE_PRESETS)
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._panel(pilot, app)
            scr._featlib_edit_preset(
                dict(P._find_preset("loxP"), preset=True))
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.AddFeatureModal)
        assert P._FEATURE_PRESETS == before

    async def test_editing_a_preset_keeps_its_accession_in_the_prefill(self):
        """The Edit button receives the row ALREADY normalised (that is what
        carries provenance into an insert). Converting it a second time would
        read a `source` key that is no longer there and hand the modal a note
        with the accession stripped out. Regression guard for 2026-09-15."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._panel(pilot, app)
            scr._refresh_featlib_table(filter_str="loxP")
            await pilot.pause()
            row = next(i for i, (_k, e) in enumerate(scr._featlib_rows)
                       if e.get("name") == "loxP")
            scr.query_one("#syn-featlib-table").move_cursor(row=row)
            await pilot.pause()
            entry = scr._featlib_selected_entry()
            scr._featlib_edit_preset(entry)
            await pilot.pause()
            await pilot.pause()
            modal = app.screen
            assert isinstance(modal, sc.AddFeatureModal)
            note = modal._prefill["qualifiers"]["note"][0]
            assert "GenBank" in note
            assert "preset" not in modal._prefill


# ═══════════════════════════════════════════════════════════════════════════════
# Annotating a LOADED plasmid from the library + presets
# ═══════════════════════════════════════════════════════════════════════════════

def _circular_record(seq: str, name: str = "test"):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq(seq), id=name, name=name, description=name)
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    return rec


def _pseudo_vector() -> str:
    parts = ["ori (pUC/pMB1, high-copy)", "AmpR promoter", "AmpR (bla)",
             "T7 promoter", "lacZ-alpha"]
    spacer = "GATCGATCGATCGATCGATC"
    return spacer.join(P._find_preset(n)["sequence"] for n in parts)


class TestPresetAnnotationTransfers:

    def test_finds_every_element_in_the_transfer_shape(self):
        sc._save_features([])
        rec = _circular_record(_pseudo_vector())
        transfers, already, total = sc._preset_annotation_transfers(rec)
        assert already == 0
        assert total == len(transfers)
        labels = {t["label"] for t in transfers}
        for want in ("ori (pUC/pMB1, high-copy)", "AmpR (bla)",
                     "T7 promoter", "lacZ-alpha"):
            assert want in labels
        for t in transfers:
            for key in ("label", "type", "source_start", "source_end",
                        "source_strand", "target_start", "target_end",
                        "target_strand", "length", "qualifiers"):
                assert key in t, f"{t['label']} missing {key}"

    def test_empty_record_is_handled(self):
        assert sc._preset_annotation_transfers(_circular_record("")) == ([], 0, 0)

    def test_library_only_mode_finds_nothing_without_user_entries(self):
        sc._save_features([])
        rec = _circular_record(_pseudo_vector())
        assert sc._preset_annotation_transfers(
            rec, include_presets=False) == ([], 0, 0)

    def test_wrapped_hit_uses_the_transfer_convention(self):
        """The scan reports an origin-spanning hit as `end > len(seq)`; the
        transfer/apply path expects `target_end < target_start`. Regression
        guard for 2026-09-15 — without the conversion a wrapped marker lands
        as a feature running off the end of the plasmid."""
        sc._save_features([])
        lox = P._find_preset("loxP")["sequence"]
        head, tail = lox[:20], lox[20:]
        seq = tail + "ACGT" * 30 + head
        transfers, _already, _total = sc._preset_annotation_transfers(
            _circular_record(seq))
        wrapped = [t for t in transfers if t["label"] == "loxP"]
        assert wrapped, f"loxP not found; got {[t['label'] for t in transfers]}"
        t = wrapped[0]
        assert t["target_end"] < t["target_start"]
        assert t["target_start"] == len(seq) - 20
        assert t["target_end"] == len(lox) - 20

    def test_already_annotated_spans_are_not_offered_again(self):
        """Re-running on an annotated plasmid must be a no-op, not a way to
        stack a second copy of every element."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        sc._save_features([])
        rec = _circular_record(_pseudo_vector())
        transfers, _already, _total = sc._preset_annotation_transfers(rec)
        assert transfers
        for t in transfers:
            if t["target_end"] < t["target_start"]:
                continue
            rec.features.append(SeqFeature(
                FeatureLocation(t["target_start"], t["target_end"],
                                strand=t["target_strand"]),
                type=t["type"], qualifiers={"label": [t["label"]]},
            ))
        again, already, _total = sc._preset_annotation_transfers(rec)
        assert again == []
        assert already == len(transfers)

    def test_a_differently_labelled_feature_does_not_suppress_the_match(self):
        """Only an identical label over an identical span counts as already
        annotated — a user's own unrelated feature at the same coords must
        not hide a preset match."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        sc._save_features([])
        rec = _circular_record(_pseudo_vector())
        transfers, _already, _total = sc._preset_annotation_transfers(rec)
        t = next(x for x in transfers
                 if x["target_end"] > x["target_start"])
        rec.features.append(SeqFeature(
            FeatureLocation(t["target_start"], t["target_end"],
                            strand=t["target_strand"]),
            type="misc_feature", qualifiers={"label": ["something else"]},
        ))
        again, already, _total = sc._preset_annotation_transfers(rec)
        assert already == 0
        assert t["label"] in {x["label"] for x in again}


@pytest.mark.asyncio
class TestAnnotateFromPresetsAction:

    async def test_refuses_with_no_plasmid_loaded(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            app._current_record = None
            app.action_annotate_from_presets()
            await pilot.pause()
            assert not isinstance(app.screen, sc.AnnotationTransferModal)

    async def test_opens_the_preview_and_applies_onto_the_record(self):
        sc._save_features([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            rec = _circular_record(_pseudo_vector(), "pTest")
            app._apply_record(rec)
            await pilot.pause()
            await pilot.pause()
            assert len(app._current_record.features) == 0
            app.action_annotate_from_presets()
            for _ in range(12):
                await pilot.pause()
                if isinstance(app.screen, sc.AnnotationTransferModal):
                    break
            assert isinstance(app.screen, sc.AnnotationTransferModal)
            n_offered = len(app.screen._transfers)
            assert n_offered >= 4
            app.screen.query_one("#btn-annot-apply").press()
            for _ in range(8):
                await pilot.pause()
                if not isinstance(app.screen, sc.AnnotationTransferModal):
                    break
            feats = app._current_record.features
            assert len(feats) == n_offered

    async def test_applied_features_cover_the_bases_they_name(self):
        """The acceptance test: an annotation that does not sit on its own
        sequence is worse than no annotation."""
        sc._save_features([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            rec = _circular_record(_pseudo_vector(), "pTest")
            app._apply_record(rec)
            await pilot.pause()
            await pilot.pause()
            app.action_annotate_from_presets()
            for _ in range(12):
                await pilot.pause()
                if isinstance(app.screen, sc.AnnotationTransferModal):
                    break
            app.screen.query_one("#btn-annot-apply").press()
            for _ in range(8):
                await pilot.pause()
                if not isinstance(app.screen, sc.AnnotationTransferModal):
                    break
            out = app._current_record
            seq = str(out.seq).upper()
            assert out.features
            for feat in out.features:
                label = sc._feat_label(feat)
                preset = P._find_preset(label)
                assert preset is not None, f"unexpected feature {label!r}"
                # `SeqFeature.extract` already reverse-complements a
                # minus-strand feature, so the extraction reads as the
                # element itself on EITHER strand. That is the whole claim
                # an annotation makes.
                got = str(feat.extract(out.seq)).upper()
                assert got == preset["sequence"], (
                    f"{label} was annotated over the wrong bases"
                )
            assert len(seq) == len(_pseudo_vector())

    async def test_undo_restores_the_unannotated_record(self):
        """Bulk annotation is one undo, all-or-nothing."""
        sc._save_features([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app._apply_record(_circular_record(_pseudo_vector(), "pTest"))
            await pilot.pause()
            await pilot.pause()
            app.action_annotate_from_presets()
            for _ in range(12):
                await pilot.pause()
                if isinstance(app.screen, sc.AnnotationTransferModal):
                    break
            app.screen.query_one("#btn-annot-apply").press()
            for _ in range(8):
                await pilot.pause()
                if not isinstance(app.screen, sc.AnnotationTransferModal):
                    break
            assert app._current_record.features
            app.action_undo()
            await pilot.pause()
            await pilot.pause()
            assert len(app._current_record.features) == 0


@pytest.mark.asyncio
class TestAnnotateFromPresetsEndpoint:
    """The headless counterpart. `app.call_from_thread` needs a running app,
    so these drive the handler through a real Pilot."""

    async def _run(self, payload, *, load=True):
        sc._save_features([])
        app = sc.PlasmidApp()
        out = {}
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            if load:
                app._apply_record(_circular_record(_pseudo_vector(), "pTest"))
                await pilot.pause()
                await pilot.pause()
            fn, _write = sc._state._AGENT_HANDLERS["annotate-from-presets"]
            import asyncio
            out["result"] = await asyncio.to_thread(fn, app, payload)
            for _ in range(8):
                await pilot.pause()
            out["features"] = list(app._current_record.features) if load else []
        return out

    async def test_refuses_without_a_loaded_record(self):
        got = await self._run({}, load=False)
        body, code = got["result"]
        assert code == 422
        assert "no plasmid loaded" in body["error"]

    async def test_defaults_to_a_dry_run(self):
        got = await self._run({})
        r = got["result"]
        assert r["applied"] is False
        assert r["count"] >= 4
        assert got["features"] == [], "a dry run must not touch the record"

    async def test_apply_true_writes_the_features(self):
        got = await self._run({"apply": True})
        r = got["result"]
        assert r["applied"] is True
        assert len(got["features"]) == r["count"]

    async def test_include_presets_false_scans_only_the_user_library(self):
        got = await self._run({"include_presets": False})
        assert got["result"]["count"] == 0

    async def test_unrecognised_affirmative_is_refused_not_previewed(self):
        """Mirrors `transfer-annotations`: an agent that thinks it committed
        must not silently receive a preview."""
        got = await self._run({"commit": True})
        body, code = got["result"]
        assert code == 400
        assert "not a recognised switch" in body["error"]

    async def test_rejects_a_non_boolean_include_presets(self):
        got = await self._run({"include_presets": "yes"})
        assert got["result"][1] == 400

    async def test_reapplying_reports_everything_already_annotated(self):
        sc._save_features([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app._apply_record(_circular_record(_pseudo_vector(), "pTest"))
            await pilot.pause()
            await pilot.pause()
            fn, _w = sc._state._AGENT_HANDLERS["annotate-from-presets"]
            import asyncio
            first = await asyncio.to_thread(fn, app, {"apply": True})
            for _ in range(8):
                await pilot.pause()
            assert first["count"] > 0
            second = await asyncio.to_thread(fn, app, {"apply": True})
            for _ in range(8):
                await pilot.pause()
            assert second["count"] == 0
            assert second["already_annotated"] == first["count"]
            assert len(app._current_record.features) == first["count"]


# ═══════════════════════════════════════════════════════════════════════════════
# Hardening — edge cases found by the 2026-09-15 sweep
# ═══════════════════════════════════════════════════════════════════════════════

class TestCorruptLibraryDoesNotBreakTheMerge:
    """`features.json` is user-editable and `_load_features` only checks that a
    row is a dict. One bad row must cost the user that row, not the feature."""

    def test_unhashable_name_does_not_raise(self):
        """Regression guard for 2026-09-15: building the shadow key from raw
        values raised `TypeError: unhashable type` on a list/dict in `name`,
        which took out the Synthesis panel, the preset browser's in-library
        marks and the whole annotate scan at once."""
        bad = [
            {"name": ["a list"], "feature_type": "CDS", "sequence": "ATG"},
            {"name": {"a": 1}, "feature_type": "CDS"},
            {"name": "ok", "feature_type": ["also a list"]},
        ]
        merged = P._merge_presets_with_library(bad)
        assert len(merged) == len(P._FEATURE_PRESETS) + len(bad)

    def test_corrupt_rows_stay_visible_so_the_user_can_fix_them(self):
        bad = [{"name": ["a list"], "feature_type": "CDS"}]
        merged = P._merge_presets_with_library(bad)
        assert merged[0]["name"] == ["a list"]
        assert not merged[0].get("preset")

    def test_a_corrupt_key_shadows_nothing(self):
        merged = P._merge_presets_with_library(
            [{"name": ["loxP"], "feature_type": "misc_recomb"}])
        assert any(e["name"] == "loxP" and e.get("preset") for e in merged)

    def test_shadow_key_is_total(self):
        assert P._shadow_key({}) == ("", "")
        assert P._shadow_key("not a dict") == ("", "")  # type: ignore[arg-type]
        assert P._shadow_key({"name": 5, "feature_type": None}) == ("5", "None")

    def test_scan_survives_a_corrupt_library_row(self):
        sc._save_features([{"name": ["bad"], "feature_type": "CDS",
                            "sequence": "ATGAAATAA"}])
        seq = P._find_preset("T7 promoter")["sequence"] * 2
        hits = sc._annotate_seq_from_feature_library(seq, include_presets=True)
        assert any(h["name"] == "T7 promoter" for h in hits)


class TestPresetStrandCoercionNeverRaises:
    """A bad strand used to surface as a traceback inside a modal callback."""

    @pytest.mark.parametrize("bad,want", [
        ("oops", 1), (None, 1), ([], 1), ({}, 1), (7, 1), ("", 1),
        ("-1", -1), (0, 0), (2, 2), (-1, -1), (1.0, 1),
    ])
    def test_coercion(self, bad, want):
        assert P._coerce_preset_strand(bad) == want

    def test_converter_does_not_raise_on_a_bad_strand(self):
        entry = P._preset_to_library_entry(
            {"name": "x", "sequence": "ACGT", "strand": "oops"})
        assert entry["strand"] == 1


class TestAnnotationPreviewIsCapped:

    def _tandem(self, copies: int) -> str:
        return (P._find_preset("T7 promoter")["sequence"] + "AAAA") * copies

    def test_pathological_repeat_is_capped_and_says_so(self):
        """A tandem-repeat construct produced 4,000 matches (measured). A
        preview that size is unreviewable, and applying it staples 4,000
        features onto the record in one undo step."""
        sc._save_features([])
        rec = _circular_record(self._tandem(900))
        transfers, _already, total = sc._preset_annotation_transfers(rec)
        assert len(transfers) == sc._PRESET_ANNOT_MAX_TRANSFERS
        assert total == 900
        assert total > len(transfers), "truncation must be visible to callers"

    def test_cap_is_configurable_and_floors_at_one(self):
        sc._save_features([])
        rec = _circular_record(self._tandem(20))
        transfers, _a, total = sc._preset_annotation_transfers(
            rec, max_transfers=5)
        assert len(transfers) == 5 and total == 20
        transfers, _a, _t = sc._preset_annotation_transfers(
            rec, max_transfers=0)
        assert len(transfers) == 1, "a zero cap must not silently return none"

    def test_a_normal_plasmid_is_nowhere_near_the_cap(self):
        sc._save_features([])
        transfers, _a, total = sc._preset_annotation_transfers(
            _circular_record(_pseudo_vector()))
        assert total == len(transfers) < 50


class TestScanSurvivesMalformedRecords:

    @pytest.mark.parametrize("rec", [None, object()])
    def test_non_records_return_empty(self, rec):
        assert sc._preset_annotation_transfers(rec) == ([], 0, 0)

    def test_unreadable_feature_does_not_abort_the_scan(self):
        """A hand-edited GenBank can carry a feature whose location won't
        resolve. Skipping it costs one duplicate check, not the whole run."""
        sc._save_features([])
        rec = _circular_record(_pseudo_vector())
        rec.features = ["not a feature", None, 42]
        transfers, _a, _t = sc._preset_annotation_transfers(rec)
        assert transfers

    def test_lowercase_and_missing_topology_are_handled(self):
        sc._save_features([])
        seq = _pseudo_vector()
        lower = _circular_record(seq.lower())
        assert sc._preset_annotation_transfers(lower)[0]
        no_topo = _circular_record(seq)
        no_topo.annotations = {}
        assert sc._preset_annotation_transfers(no_topo)[0]


class TestImportEndpointHardening:

    def test_repeated_names_import_once(self):
        """`count` must equal the number of entries the library gained — the
        pre-fix response said `count: 2, replaced: 1` for one new entry."""
        sc._save_features([])
        r = _call("import-feature-preset", {"names": ["loxP", "loxP", "LOXP"]})
        assert r["count"] == 1 and r["replaced"] == 0
        assert len(sc._load_features()) == 1

    def test_name_and_names_together_are_refused(self):
        sc._save_features([])
        body, code = _call("import-feature-preset",
                           {"name": "FRT", "names": ["loxP"]})
        assert code == 400
        assert "not both" in body["error"]
        assert sc._load_features() == []

    def test_list_rejects_a_non_boolean_include_sequence(self):
        body, code = _call("list-feature-presets", {"include_sequence": "yes"})
        assert code == 400
        assert "must be a boolean" in body["error"]

    def test_list_still_accepts_a_real_boolean(self):
        assert _call("list-feature-presets",
                     {"include_sequence": False})["count"] > 0


@pytest.mark.asyncio
class TestAnnotateEndpointHardening:

    async def _run(self, payload):
        sc._save_features([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app._apply_record(_circular_record(
                (P._find_preset("T7 promoter")["sequence"] + "AAAA") * 900,
                "pRepeat"))
            await pilot.pause()
            await pilot.pause()
            fn, _w = sc._state._AGENT_HANDLERS["annotate-from-presets"]
            import asyncio
            return await asyncio.to_thread(fn, app, payload)

    async def test_truncation_is_reported_not_hidden(self):
        r = await self._run({})
        assert r["truncated"] is True
        assert r["total_found"] == 900
        assert r["count"] == sc._PRESET_ANNOT_MAX_TRANSFERS

    async def test_max_transfers_is_honoured_and_validated(self):
        r = await self._run({"max_transfers": 7})
        assert r["count"] == 7 and r["total_found"] == 900
        assert (await self._run({"max_transfers": 0}))[1] == 400
        assert (await self._run({"max_transfers": "lots"}))[1] == 400


@pytest.mark.asyncio
class TestPresetModalHardening:

    async def _open(self, pilot, app):
        app.push_screen(sc.FeaturePresetsModal())
        await pilot.pause()
        await pilot.pause()
        return app.screen

    async def test_every_action_is_safe_with_no_rows_shown(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._open(pilot, app)
            scr.query_one("#fpre-search").value = "zzzzzzzzzz"
            await pilot.pause()
            await pilot.pause()
            assert scr.query_one("#fpre-table").row_count == 0
            assert scr._current() is None
            scr.action_toggle_selection()
            scr.query_one("#btn-fpre-all").press()
            await pilot.pause()
            assert scr._selected == set()
            scr.query_one("#btn-fpre-add").press()
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.FeaturePresetsModal), (
                "Add with nothing to add must notify, not dismiss"
            )

    async def test_ticks_hidden_by_the_filter_are_counted_out_loud(self):
        """A tick survives a filter change on purpose, but one the filter
        hides would otherwise be added invisibly."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._open(pilot, app)
            scr.query_one("#fpre-cat").value = "Resistance"
            await pilot.pause()
            await pilot.pause()
            scr.query_one("#btn-fpre-all").press()
            await pilot.pause()
            n = len(scr._selected)
            assert n > 0
            scr.query_one("#fpre-cat").value = "Origin"
            await pilot.pause()
            await pilot.pause()
            assert len(scr._selected) == n, "ticks must survive a filter change"
            status = str(scr.query_one("#fpre-status").render())
            assert f"{n} not shown" in status

    async def test_clear_selection_empties_the_hidden_ticks_too(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._open(pilot, app)
            scr.query_one("#fpre-cat").value = "Resistance"
            await pilot.pause()
            await pilot.pause()
            scr.query_one("#btn-fpre-all").press()
            await pilot.pause()
            scr.query_one("#fpre-cat").value = "Origin"
            await pilot.pause()
            scr.query_one("#btn-fpre-none").press()
            await pilot.pause()
            assert scr._selected == set()

    async def test_toggling_keeps_the_cursor_where_it_was(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            scr = await self._open(pilot, app)
            tbl = scr.query_one("#fpre-table")
            tbl.move_cursor(row=5)
            await pilot.pause()
            scr.action_toggle_selection()
            await pilot.pause()
            assert tbl.cursor_row == 5


class TestPresetsCanNeverReachDisk:
    """The sacred property, enforced at the source level rather than by
    hoping every future edit remembers it.

    `_load_features_with_presets` returns a merged browse/scan view. Handing
    that list to `_save_features` would bake 116 shipped entries into the
    user's `features.json` — silently, and permanently, as a side effect of
    opening a panel. Every save must be fed from `_load_features` (the user's
    own library) instead.
    """

    _SOURCES = ("splicecraft.py", "splicecraft_agent.py",
                "splicecraft_modals.py", "splicecraft_dataaccess.py")

    def _functions_calling(self, name: str) -> "set[tuple[str, str]]":
        import ast
        from pathlib import Path
        repo = Path(__file__).resolve().parent.parent
        out: set = set()
        for fname in self._SOURCES:
            tree = ast.parse((repo / fname).read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    continue
                for sub in ast.walk(node):
                    if not isinstance(sub, ast.Call):
                        continue
                    fn = sub.func
                    called = getattr(fn, "id", None) or getattr(fn, "attr", None)
                    if called == name:
                        out.add((fname, node.name))
        return out

    def test_no_function_both_merges_and_saves(self):
        merges = self._functions_calling("_load_features_with_presets")
        saves = self._functions_calling("_save_features")
        overlap = merges & saves
        assert not overlap, (
            "these function(s) call BOTH _load_features_with_presets and "
            f"_save_features, so a merged list can reach disk: {sorted(overlap)}"
        )

    def test_the_merged_view_has_only_read_only_callers(self):
        """Staleguard: if a new caller appears, decide deliberately whether it
        is read-only, then add it here."""
        # The definition itself doesn't call the name, so only real call
        # sites appear here.
        expected = {
            ("splicecraft.py", "_annotate_seq_from_feature_library"),
            ("splicecraft.py", "_refresh_featlib_table"),
        }
        assert self._functions_calling("_load_features_with_presets") == expected

    def test_a_round_trip_through_the_screen_buffer_keeps_presets_out(self):
        """End-to-end version of the same claim: open the library screen on a
        machine with presets available, save, and `features.json` must still
        hold only what the user put there."""
        sc._save_features([{"name": "mine", "feature_type": "CDS",
                            "sequence": "ATGAAATAA", "strand": 1, "color": ""}])
        assert len(sc._load_features_with_presets()) > 1
        entries = sc._load_features()
        sc._save_features(entries)
        after = sc._load_features()
        assert [e["name"] for e in after] == ["mine"]
