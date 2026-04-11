"""
conftest — shared fixtures for the SpliceCraft test suite.

The test suite imports `splicecraft` as a module. Because splicecraft.py lives in
the repo root (not packaged), we prepend the repo root to sys.path once here.
"""
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import pytest

# Delay importing splicecraft until fixtures run — it runs _check_deps() at
# import time, which is fine, but keeps the side-effects traceable to a fixture
# call rather than collection.


@pytest.fixture(scope="session")
def splicecraft_module():
    """Import splicecraft once per session and hand it to tests that need
    symbols from it."""
    import splicecraft
    return splicecraft


@pytest.fixture
def tiny_record():
    """A synthetic 120 bp SeqRecord with one CDS feature and one misc_feature.

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
    # Length is ~120 bp — exact value not load-bearing, just long enough
    # to hold several restriction sites and two features.
    assert len(seq_str) > 100

    rec = SeqRecord(Seq(seq_str), id="TEST001", name="TEST001",
                    description="Synthetic test plasmid (120 bp)")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"]      = "circular"
    # Forward-strand CDS 0..24 (ATGAAAGATCTGGAATTCAAAGGG) — starts M, stops at TAG later
    rec.features.append(SeqFeature(
        FeatureLocation(0, 27, strand=1),
        type="CDS",
        qualifiers={"gene": ["testA"], "product": ["test protein A"]},
    ))
    # Reverse-strand misc_feature
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
def isolated_library(tmp_path, monkeypatch, splicecraft_module):
    """Redirect `_LIBRARY_FILE` to a tmp path so tests can't touch the real
    `plasmid_library.json`. Also wipes the in-memory cache so lookups re-read.
    """
    tmp_lib = tmp_path / "plasmid_library.json"
    monkeypatch.setattr(splicecraft_module, "_LIBRARY_FILE", tmp_lib)
    monkeypatch.setattr(splicecraft_module, "_library_cache", None)
    return tmp_lib
