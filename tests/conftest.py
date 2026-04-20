"""
conftest — shared fixtures for the SpliceCraft test suite.

CRITICAL SAFETY RULE: the _protect_user_data autouse fixture redirects ALL
JSON persistence files (plasmid_library.json, parts_bin.json, primers.json)
to a temporary directory for the ENTIRE test session. No test, fixture, or
ad-hoc import can ever touch the user's real data files.

If you add a new persistence file to splicecraft.py, you MUST add it to the
_DATA_FILES list in _protect_user_data below.
"""
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import pytest


# ═══════════════════════════════════════════════════════════════════════════════
# AUTOUSE: protect user data files from every single test
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture(autouse=True, scope="function")
def _protect_user_data(tmp_path, monkeypatch):
    """Redirect ALL user-data JSON files to tmp_path for EVERY test.

    This is autouse=True so it fires automatically — no test needs to
    explicitly request it. Even if a test forgets `isolated_library` or
    calls `_save_library` directly, the write goes to a throwaway temp
    dir, not the real plasmid_library.json.

    If a future developer adds a new JSON persistence file to
    splicecraft.py, they MUST add the corresponding (attr, cache_attr)
    pair to _DATA_FILES below or the test_no_real_files_touched test
    will catch it.
    """
    import splicecraft as sc

    _DATA_FILES = [
        ("_LIBRARY_FILE",       "_library_cache"),
        ("_PARTS_BIN_FILE",     "_parts_bin_cache"),
        ("_PRIMERS_FILE",       "_primers_cache"),
        ("_CODON_TABLES_FILE",  "_codon_tables_cache"),
        ("_FEATURES_FILE",      "_features_cache"),
        ("_PLANNOTATE_CHECK_CACHE", None),  # not a file, but reset for isolation
    ]

    for file_attr, cache_attr in _DATA_FILES:
        if file_attr == "_PLANNOTATE_CHECK_CACHE":
            monkeypatch.setattr(sc, file_attr, None)
            continue
        # Redirect the file path to a temp location
        real_path = getattr(sc, file_attr)
        tmp_file = tmp_path / real_path.name
        monkeypatch.setattr(sc, file_attr, tmp_file)
        # Clear the in-memory cache so the next load reads from the tmp file
        if cache_attr:
            monkeypatch.setattr(sc, cache_attr, None)

    # Crash-recovery autosave dir: redirect so tests can't leave files in
    # the user's real _DATA_DIR/crash_recovery on disk.
    monkeypatch.setattr(sc, "_CRASH_RECOVERY_DIR", tmp_path / "crash_recovery")


# ═══════════════════════════════════════════════════════════════════════════════
# Session-scoped module import
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture(scope="session")
def splicecraft_module():
    """Import splicecraft once per session and hand it to tests that need
    symbols from it."""
    import splicecraft
    return splicecraft


# ═══════════════════════════════════════════════════════════════════════════════
# Common test fixtures
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def tiny_record():
    """A synthetic ~120 bp SeqRecord with one CDS feature and one misc_feature.

    Small enough to render fast in smoke tests; large enough to have a valid
    6-bp recognition site for EcoRI (GAATTC). Built without touching NCBI.
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    seq_str = (
        "ATGAAAGATCTGGAATTCAAAGGGCCCTAGAAAGCATGCAAAATCGATGTCGACAAAGAATTC"
        "AAATCCTAGGAAAAGGATCCAAAACTCGAGCCCAAAAAATTTGGGCCCAAAATCGA"
        "TAG"
    )
    assert len(seq_str) > 100

    rec = SeqRecord(Seq(seq_str), id="TEST001", name="TEST001",
                    description="Synthetic test plasmid (120 bp)")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"]      = "circular"
    rec.features.append(SeqFeature(
        FeatureLocation(0, 27, strand=1),
        type="CDS",
        qualifiers={"gene": ["testA"], "product": ["test protein A"]},
    ))
    rec.features.append(SeqFeature(
        FeatureLocation(50, 80, strand=-1),
        type="misc_feature",
        qualifiers={"label": ["reverse region"]},
    ))
    return rec


@pytest.fixture
def tiny_gb_path(tmp_path, tiny_record):
    """Write `tiny_record` to a .gb file in a tmp dir and return the path."""
    from Bio import SeqIO
    p = tmp_path / "tiny.gb"
    SeqIO.write(tiny_record, str(p), "genbank")
    return str(p)


@pytest.fixture
def isolated_library(tmp_path, monkeypatch):
    """Redirect `_LIBRARY_FILE` to a tmp path. Kept for backward compat with
    tests that explicitly request it — _protect_user_data already handles
    this automatically, so this fixture is now redundant but harmless."""
    import splicecraft as sc
    tmp_lib = tmp_path / "plasmid_library.json"
    monkeypatch.setattr(sc, "_LIBRARY_FILE", tmp_lib)
    monkeypatch.setattr(sc, "_library_cache", None)
    return tmp_lib
