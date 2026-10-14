"""`.dna` qualifier values Biopython's reader skips must survive import.

The other editor stores a qualifier value as a `<V>` element whose ATTRIBUTE
says what kind of value it is. Biopython's `.dna` reader keeps only `text`,
`predef` and `int`, so a flag stored as `<V bool="1"/>` came in as
`/pseudo=[]` — and the GenBank writer writes a key with no values as nothing.
Measured on a real library: 3,338 `/pseudo`, 2,116 `/rpt_type` and 101
`/ribosomal_slippage` across 52 imported plasmids, all gone, and GenBank
export then refused every one of those files.
"""
from __future__ import annotations

import pytest

import splicecraft as sc
import splicecraft_fileio as fio


SEQ = ("ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCG"
       "AATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCC"
       "TTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCA"
       "ACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGG")


def _dna(features_xml: str, seq: str = SEQ) -> bytes:
    return (fio._build_commercialsaas_cookie_packet()
            + fio._build_commercialsaas_dna_packet(seq, circular=True)
            + fio._build_commercialsaas_packet(
                0x0A, f"<Features>{features_xml}</Features>".encode())
            + fio._build_commercialsaas_packet(
                0x06, b"<Notes><Type>Synthetic</Type></Notes>"))


def _feature(quals_xml: str, *, ftype: str = "gene", rng: str = "11-100",
             name: str = "geneX") -> str:
    return (f'<Feature name="{name}" type="{ftype}" directionality="1">'
            f'<Segment range="{rng}" color="#ff0000" type="standard"/>'
            f"{quals_xml}</Feature>")


def _load(tmp_path, data: bytes, name: str = "t.dna"):
    p = tmp_path / name
    p.write_bytes(data)
    return sc.load_genbank(str(p))


def _only(rec, ftype: str):
    hits = [f for f in rec.features if f.type == ftype]
    assert len(hits) == 1, [f.type for f in rec.features]
    return hits[0]


class TestEveryValueKindArrives:
    def test_a_flag_is_a_flag(self, tmp_path):
        rec = _load(tmp_path, _dna(_feature(
            '<Q name="pseudo"><V bool="1"/></Q>'
            '<Q name="note"><V text="kept"/></Q>')))
        q = _only(rec, "gene").qualifiers
        assert q["pseudo"] == [""]
        assert q["note"] == ["kept"]

    def test_a_flag_stored_as_unset_is_absent(self, tmp_path):
        rec = _load(tmp_path, _dna(_feature(
            '<Q name="pseudo"><V bool="0"/></Q>')))
        assert "pseudo" not in _only(rec, "gene").qualifiers

    @pytest.mark.parametrize("rng,want", [
        ("12-40", "12..40"),
        ("0-0", ""),                  # the editor's "unset"
    ])
    def test_a_range(self, tmp_path, rng, want):
        rec = _load(tmp_path, _dna(_feature(
            f'<Q name="rpt_unit_range"><V range="{rng}"/></Q>',
            ftype="repeat_region")))
        assert _only(rec, "repeat_region").qualifiers["rpt_unit_range"] == [want]

    @pytest.mark.parametrize("attrs,want", [
        ('year="1952" month="10" day="21"', "21-Oct-1952"),
        ('year="2016" month="3"', "Mar-2016"),
        ('year="2016"', "2016"),
        ('year="0" month="0" day="0"', ""),
    ])
    def test_a_date(self, tmp_path, attrs, want):
        rec = _load(tmp_path, _dna(_feature(
            f'<Q name="collection_date"><V {attrs}/></Q>',
            ftype="misc_feature")))
        assert (_only(rec, "misc_feature").qualifiers["collection_date"]
                == [want])

    def test_an_empty_value_is_present_and_empty(self, tmp_path):
        rec = _load(tmp_path, _dna(_feature(
            '<Q name="anticodon"><V/></Q>', ftype="tRNA")))
        assert _only(rec, "tRNA").qualifiers["anticodon"] == [""]

    def test_mixed_kinds_keep_their_order(self, tmp_path):
        rec = _load(tmp_path, _dna(_feature(
            '<Q name="note"><V text="first"/><V/><V text="third"/></Q>')))
        assert _only(rec, "gene").qualifiers["note"] == ["first", "", "third"]

    def test_the_source_rows_values_arrive_too(self, tmp_path):
        rec = _load(tmp_path, _dna(
            _feature('<Q name="collection_date"><V year="2016"/></Q>'
                     '<Q name="focus"><V bool="1"/></Q>',
                     ftype="source", rng=f"1-{len(SEQ)}", name="source")
            + _feature('<Q name="pseudo"><V bool="1"/></Q>')))
        src = _only(rec, "source").qualifiers
        assert src["collection_date"] == ["2016"]
        assert src["focus"] == [""]
        assert _only(rec, "gene").qualifiers["pseudo"] == [""]


class TestTheyReachTheLibraryAndTheExport:
    def test_a_flag_survives_storage(self, tmp_path):
        rec = _load(tmp_path, _dna(_feature(
            '<Q name="pseudo"><V bool="1"/></Q>'
            '<Q name="ribosomal_slippage"><V bool="1"/></Q>')))
        text = sc._record_to_gb_text(rec)
        back = _only(sc._gb_text_to_record(text), "gene").qualifiers
        assert back["pseudo"] == [""]
        assert back["ribosomal_slippage"] == [""]

    def test_a_flag_is_written_the_way_the_format_spells_it(self, tmp_path):
        rec = _load(tmp_path, _dna(_feature(
            '<Q name="pseudo"><V bool="1"/></Q>')))
        lines = sc._record_to_gb_text(rec).splitlines()
        assert " " * 21 + "/pseudo" in lines
        assert not any('/pseudo="' in ln for ln in lines)

    def test_the_export_accepts_it(self, tmp_path):
        rec = _load(tmp_path, _dna(_feature(
            '<Q name="pseudo"><V bool="1"/></Q>'
            '<Q name="rpt_unit_range"><V range="12-40"/></Q>'
            '<Q name="anticodon"><V/></Q>')))
        out = tmp_path / "out.gb"
        fio._export_genbank_to_path(rec, out)
        back = _only(sc.load_genbank(str(out)), "gene").qualifiers
        assert back["pseudo"] == [""]
        assert back["rpt_unit_range"] == ["12..40"]

    def test_a_key_with_no_values_does_not_block_the_export(self, tmp_path):
        # Whatever put it there, a key holding NO value cannot be written —
        # the export must not refuse the file over it.
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq(SEQ), id="pX", name="pX",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec.features.append(SeqFeature(
            FeatureLocation(10, 100, strand=1), type="gene",
            qualifiers={"gene": ["geneX"], "pseudo": []}))
        out = tmp_path / "out.gb"
        fio._export_genbank_to_path(rec, out)
        assert out.exists()


class TestOlderImportsStillSplice:
    def test_an_entry_the_previous_reader_made_matches_its_original(
            self, tmp_path):
        # 1.3.0-1.3.8 read the file WITHOUT its flags. Today's reading of the
        # same file has them, so an unedited entry those versions imported
        # must be compared with THEIR reading, or its `.dna` export is rebuilt
        # from scratch and drops everything only the original holds.
        data = _dna(_feature('<Q name="pseudo"><V bool="1"/></Q>'))
        p = tmp_path / "old.dna"
        p.write_bytes(data)
        rec = sc.load_genbank(str(p), dna_reader_revision=1)
        assert not _only(rec, "gene").qualifiers.get("pseudo")
        entry = sc._record_to_library_entry(rec, p)
        entry["gb_text"] = entry["gb_text"].replace(
            f"Created by SpliceCraft v{sc.__version__}",
            "Created by SpliceCraft v1.3.8")
        assert fio._dna_sidecar_matches_entry(data, entry)

    def test_an_entry_todays_reader_made_matches_its_original(self, tmp_path):
        data = _dna(_feature('<Q name="pseudo"><V bool="1"/></Q>'))
        rec = _load(tmp_path, data, "new.dna")
        entry = sc._record_to_library_entry(rec, tmp_path / "new.dna")
        assert fio._dna_sidecar_matches_entry(data, entry)

    @pytest.mark.parametrize("stamp,want", [
        (None, 0), ("1.2.71", 0), ("1.3.0", 1), ("1.3.8", 1), ("1.3.9", 2),
        ("1.4.0", 2),
    ])
    def test_which_reading_made_an_entry(self, stamp, want):
        text = "" if stamp is None else f"Created by SpliceCraft v{stamp} on x"
        assert fio._entry_dna_reader_revision(text) == want


class TestRecoverQualifiersFromDna:
    """`recover-qualifiers-from-dna` puts the dropped values back into
    plasmids an older version imported, from the `.dna` originals on disk."""

    def _seed(self, tmp_path, quals_xml='<Q name="pseudo"><V bool="1"/></Q>',
              *, entry_quals=None, sidecars=1, rng="11-100"):
        data = _dna(_feature(quals_xml, rng=rng))
        p = tmp_path / "pOld.dna"
        p.write_bytes(data)
        rec = sc.load_genbank(str(p), dna_reader_revision=1)
        for k, v in (entry_quals or {}).items():
            _only(rec, "gene").qualifiers[k] = v
        entry = sc._record_to_library_entry(rec, p)
        entry["gb_text"] = entry["gb_text"].replace(
            f"Created by SpliceCraft v{sc.__version__}",
            "Created by SpliceCraft v1.3.8")
        for i in range(sidecars):
            side = sc._dna_sidecar_path(f"{entry['id']}{i or ''}")
            side.parent.mkdir(parents=True, exist_ok=True)
            side.write_bytes(data)
        sc._save_collections([{"name": "Work", "plasmids": [entry],
                               "saved": "2026-01-01"}])
        sc._set_active_collection_name("Work")
        sc._restore_library_from_active_collection()
        return entry, data

    def _stored(self):
        c = next(c for c in sc._load_collections() if c["name"] == "Work")
        return sc._rehydrate_entry(c["plasmids"][0])

    def test_dry_run_reports_without_writing(self, tmp_path):
        entry, _ = self._seed(tmp_path)
        r = sc._h_recover_qualifiers_from_dna(None, {})
        assert r["ok"] and r["dry_run"] is True
        assert r["updated_count"] == 1
        assert r["added_total"] == {"pseudo": 1}
        assert self._stored()["gb_text"] == entry["gb_text"]

    def test_apply_puts_the_values_back_and_nothing_else(self, tmp_path):
        entry, _ = self._seed(
            tmp_path, '<Q name="pseudo"><V bool="1"/></Q>'
                      '<Q name="rpt_unit_range"><V range="12-40"/></Q>'
                      '<Q name="note"><V text="kept"/></Q>')
        r = sc._h_recover_qualifiers_from_dna(None, {"dry_run": False})
        assert r["updated_count"] == 1, r
        before = sc._gb_text_to_record(entry["gb_text"], cache=False)
        after = sc._gb_text_to_record(self._stored()["gb_text"], cache=False)
        assert str(after.seq) == str(before.seq)
        assert len(after.features) == len(before.features)
        q = _only(after, "gene").qualifiers
        assert q["pseudo"] == [""]
        assert q["rpt_unit_range"] == ["12..40"]
        for k, v in _only(before, "gene").qualifiers.items():
            assert q[k] == v, k

    def test_a_value_the_entry_already_has_is_left_alone(self, tmp_path):
        self._seed(tmp_path, '<Q name="rpt_unit_range"><V range="12-40"/></Q>',
                   entry_quals={"rpt_unit_range": ["5..9"]})
        r = sc._h_recover_qualifiers_from_dna(None, {"dry_run": False})
        assert r["updated_count"] == 0
        # Left alone by design — not caught late by the proof that nothing
        # but added values changed.
        assert r["skipped"] == []
        after = sc._gb_text_to_record(self._stored()["gb_text"], cache=False)
        assert _only(after, "gene").qualifiers["rpt_unit_range"] == ["5..9"]

    def test_a_feature_that_moved_is_not_touched(self, tmp_path):
        entry, _ = self._seed(tmp_path)
        rec = sc._gb_text_to_record(entry["gb_text"], cache=False)
        from Bio.SeqFeature import FeatureLocation
        _only(rec, "gene").location = FeatureLocation(20, 100, strand=1)
        entry["gb_text"] = sc._record_to_gb_text(rec)
        sc._save_collections([{"name": "Work", "plasmids": [entry],
                               "saved": "2026-01-01"}])
        sc._restore_library_from_active_collection()
        r = sc._h_recover_qualifiers_from_dna(None, {"dry_run": False})
        assert r["updated_count"] == 0

    def test_originals_that_disagree_are_reported_not_guessed(self, tmp_path):
        entry, _ = self._seed(tmp_path)
        other = _dna(_feature('<Q name="ribosomal_slippage"><V bool="1"/></Q>'))
        side = sc._dna_sidecar_path("someOtherName")
        side.write_bytes(other)
        r = sc._h_recover_qualifiers_from_dna(None, {"dry_run": False})
        assert r["updated_count"] == 0
        assert len(r["skipped"]) == 1
        assert "different values" in r["skipped"][0]["reason"]
        assert self._stored()["gb_text"] == entry["gb_text"]

    def test_identical_copies_of_an_original_are_one_answer(self, tmp_path):
        self._seed(tmp_path, sidecars=2)
        r = sc._h_recover_qualifiers_from_dna(None, {"dry_run": False})
        assert r["updated_count"] == 1, r

    def test_running_it_again_changes_nothing(self, tmp_path):
        self._seed(tmp_path)
        sc._h_recover_qualifiers_from_dna(None, {"dry_run": False})
        text = self._stored()["gb_text"]
        r = sc._h_recover_qualifiers_from_dna(None, {"dry_run": False})
        assert r["updated_count"] == 0
        assert self._stored()["gb_text"] == text

    def test_a_repaired_entry_still_matches_its_original(self, tmp_path):
        _, data = self._seed(tmp_path)
        sc._h_recover_qualifiers_from_dna(None, {"dry_run": False})
        assert fio._dna_sidecar_matches_entry(data, self._stored())

    def test_dry_run_must_be_a_boolean(self):
        r = sc._h_recover_qualifiers_from_dna(None, {"dry_run": "no"})
        assert r[1] == 400

    def test_apply_true_applies_and_says_so(self, tmp_path):
        self._seed(tmp_path)
        r = sc._h_recover_qualifiers_from_dna(None, {"apply": True})
        assert r["applied"] is True and r["dry_run"] is False
        after = sc._gb_text_to_record(self._stored()["gb_text"], cache=False)
        assert _only(after, "gene").qualifiers["pseudo"] == [""]

    @pytest.mark.parametrize("syn", ["commit", "write", "persist"])
    def test_an_unknown_apply_word_is_refused_not_previewed(self, tmp_path,
                                                            syn):
        self._seed(tmp_path)
        r = sc._h_recover_qualifiers_from_dna(None, {syn: True})
        assert r[1] == 400
