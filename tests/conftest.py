"""
conftest — shared fixtures for the SpliceCraft test suite.

CRITICAL SAFETY RULE: the _protect_user_data autouse fixture redirects
EVERY user-data JSON file enumerated in `_DATA_FILES` (plus the
crash-recovery, .dna-sidecar, and plugin directories) to a temporary
location for the lifetime of each test. No test, fixture, or ad-hoc
import can ever touch the user's real data files.

If you add a new persistence file to splicecraft.py, you MUST add it
to the _DATA_FILES list in _protect_user_data below — the
`test_no_real_files_touched` smoke test will fail otherwise.
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
        ("_LIBRARY_FILE",         "_library_cache"),
        ("_PARTS_BIN_FILE",       "_parts_bin_cache"),
        ("_PARTS_BIN_COLLECTIONS_FILE", "_parts_bin_collections_cache"),
        ("_PRIMERS_FILE",         "_primers_cache"),
        ("_PRIMER_COLLECTIONS_FILE", "_primer_collections_cache"),
        ("_CODON_TABLES_FILE",    "_codon_tables_cache"),
        ("_FEATURES_FILE",        "_features_cache"),
        ("_FEATURE_COLORS_FILE",  "_feature_colors_cache"),
        ("_GRAMMARS_FILE",        "_grammars_cache"),
        ("_ENTRY_VECTORS_FILE",   "_entry_vectors_cache"),
        ("_SETTINGS_FILE",        "_settings_cache"),
        ("_COLLECTIONS_FILE",     "_collections_cache"),
        ("_EXPERIMENTS_FILE",     "_experiments_cache"),
        ("_EXPERIMENT_PROJECTS_FILE", "_experiment_projects_cache"),
        ("_GELS_FILE",            "_gels_cache"),
        ("_PROTEIN_MOTIFS_FILE",  "_protein_motifs_cache"),
        ("_CUSTOM_ENZYMES_FILE",  "_custom_enzymes_cache"),
        ("_ENZYME_COLLECTIONS_FILE", "_enzyme_collections_cache"),
        ("_HMM_DB_CATALOG_FILE",  "_hmm_db_catalog_cache"),   # sweep #28
        ("_AGENT_TOKEN_FILE",     None),   # written when --agent-api is on
    ]

    for file_attr, cache_attr in _DATA_FILES:
        # Redirect the file path to a temp location
        real_path = getattr(sc, file_attr)
        tmp_file = tmp_path / real_path.name
        monkeypatch.setattr(sc, file_attr, tmp_file)
        # Clear the in-memory cache so the next load reads from the tmp file
        if cache_attr:
            monkeypatch.setattr(sc, cache_attr, None)

    # Force the active-collection mirror to write SYNCHRONOUSLY in tests.
    # Its async path uses a persistent background daemon whose deferred
    # writes otherwise bleed across tests — a mirror queued in one test
    # fires during the next, writing a stale snapshot into the shared
    # collections cache → flaky row counts in delete/save tests
    # (TestDeleteFocusRouting, 2026-05-29). Synchronous = deterministic.
    monkeypatch.setattr(sc, "_collection_sync_force_sync", True)

    # Crash-recovery autosave dir: redirect so tests can't leave files in
    # the user's real _DATA_DIR/crash_recovery on disk.
    monkeypatch.setattr(sc, "_CRASH_RECOVERY_DIR", tmp_path / "crash_recovery")
    # CommercialSaaS .dna sidecar storage (Phase 4d): tests that import .dna
    # files would otherwise create copies in the user's real
    # _DATA_DIR/dna_originals on disk. Same isolation pattern.
    monkeypatch.setattr(sc, "_DNA_ORIGINALS_DIR",
                          tmp_path / "dna_originals")
    # Plugin namespace (reserved): redirect so a test that triggers
    # `_check_and_stamp_data_version` doesn't `mkdir` the user's real
    # _DATA_DIR/plugins, and so any snapshot test that includes
    # `_PLUGINS_DIR` reads from the tmp tree.
    monkeypatch.setattr(sc, "_PLUGINS_DIR", tmp_path / "plugins")
    # Data-version stamp file: same reasoning — keeps every test from
    # competing for / racing on the real ~/.local/share/.../.splicecraft-
    # data-version file.
    monkeypatch.setattr(sc, "_DATA_VERSION_FILE",
                          tmp_path / ".splicecraft-data-version")
    # UI snapshot directory (Alt+D in the running app + the
    # `splicecraft logs --bundle` CLI command both write here):
    # redirect so tests of the snapshot system don't litter the
    # user's real ~/.local/share/splicecraft/ui_snapshots/.
    monkeypatch.setattr(sc, "_UI_SNAPSHOTS_DIR", tmp_path / "ui_snapshots")
    # Experiments-notebook attachments directory: image files attached
    # to experiment entries land in `<DATA_DIR>/experiments/<entry_id>/`,
    # which would otherwise litter the user's real data dir during tests
    # that exercise attach. Redirect to tmp_path so the whole experiments
    # surface (entries JSON + per-entry attach dirs) is isolated.
    monkeypatch.setattr(sc, "_EXPERIMENTS_DIR", tmp_path / "experiments")
    # HMM database downloads (sweep #28) — per-DB dir under
    # `<DATA_DIR>/hmm_databases/<id>/`. Redirect so download tests
    # never touch the user's real GBs of Pfam-A.
    monkeypatch.setattr(sc, "_HMM_DATABASES_DIR",
                          tmp_path / "hmm_databases")
    # Content-addressed plasmid blob store (v1.0.23): plasmid sequences
    # live here as immutable sha256-named blobs. Redirect so dehydrate /
    # rehydrate / GC never write to or quarantine the user's real blob dir.
    monkeypatch.setattr(sc, "_PLASMID_BLOB_DIR",
                          tmp_path / "plasmid_blobs")

    # Pre-update snapshot directory (`splicecraft update` data-safety net):
    # tests that exercise the `update` subcommand all the way through to
    # `subprocess.run` go through `_create_pre_update_snapshot`, which by
    # default writes a sibling directory next to `_DATA_DIR` — i.e. in
    # the user's REAL home directory. Redirect the env var so tests
    # always write to tmp_path. Also redirect `_DATA_DIR` to a tmp
    # location for the same reason (so any code that reads it directly
    # — e.g. unsafe-config check inside the update flow — sees a
    # writable, isolated directory).
    monkeypatch.setenv("SPLICECRAFT_UPDATE_BACKUP_DIR",
                         str(tmp_path / "update-backups"))
    monkeypatch.setattr(sc, "_DATA_DIR", tmp_path)

    # L2 chokepoint authorisation (2026-05-22). Every `_save_*` flows
    # through `_safe_save_json` which now refuses unauthorised writes.
    # The fixture has just monkey-patched every `_*_FILE` attr to
    # `tmp_path`, so writes are safe. Use `monkeypatch.setattr` so the
    # flag automatically resets to False on test teardown — function-
    # scoped, so a subsequent test starts fresh and the gate cannot
    # leak into a non-test that imports splicecraft post-pytest.
    monkeypatch.setattr(sc, "_SAVES_AUTHORIZED", True)
    monkeypatch.setattr(sc, "_SAVES_AUTHORIZED_REASON",
                          "pytest _protect_user_data")
    # One-shot backfill flags (2026-05-24 sweep). Default to TRUE in
    # tests so the load-triggered migrations are SKIPPED — most tests
    # seed handcrafted entries with `id != sanitize(name)` and would
    # see their fixture data rewritten by the backfill, which is
    # disruptive and not the test's subject. Tests that EXPLICITLY
    # exercise a backfill (or need the migrated post-state) opt in
    # by setting the relevant flag back to False before triggering
    # the load — see `TestLibraryIdBackfill` for the pattern.
    monkeypatch.setattr(sc, "_id_name_backfill_done", True)
    monkeypatch.setattr(sc, "_parts_bin_sequence_backfill_done", True)
    monkeypatch.setattr(sc, "_collections_backfill_done", True)
    # Origin-history backfill ([INV-93]) — skip by default; tests that
    # exercise it flip the relevant flag to False before loading.
    monkeypatch.setattr(sc, "_origin_history_backfill_done", True)
    monkeypatch.setattr(
        sc, "_collections_origin_history_backfill_done", True,
    )
    monkeypatch.setattr(
        sc, "_parts_bin_collections_backfill_done", True,
    )
    monkeypatch.setattr(sc, "_entry_vectors_name_trim_done", True)
    monkeypatch.setattr(sc, "_primers_name_trim_done", True)
    monkeypatch.setattr(
        sc, "_primer_collections_backfill_done", True,
    )

    # Skip the launch splash for every test by default — the splash modal
    # blocks input until dismissed, which would break every `pilot.click`
    # / `app.action_*` call. Tests that exercise the splash explicitly
    # set `_skip_splash = False` on their app instance before run_test.
    monkeypatch.setattr(sc.PlasmidApp, "_skip_splash", True)


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
def realistic_plasmid():
    """Synthetic ~2700 bp circular plasmid with realistic feature density —
    one resistance-marker-style CDS (~840 bp), one origin-of-replication
    misc_feature, one promoter, two ribosome binding sites, and several
    common restriction-enzyme sites scattered around. Sized to mirror a
    typical lab plasmid (pUC19 = 2686 bp) so render / search / packing
    code is tested under realistic load.

    Use for integration tests where the 120 bp `tiny_record` doesn't
    catch issues that only appear at full-plasmid scale (e.g. seq panel
    chunk caching, BLAST DB build time, feature-packer 2D layout under
    crowding).
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    # Backbone + features sized to look pUC19-ish. Hand-crafted to
    # contain real restriction sites at known positions:
    #   EcoRI  (GAATTC) at  396, 1980
    #   BamHI  (GGATCC) at  608
    #   HindIII (AAGCTT) at 1102
    #   XhoI   (CTCGAG) at 1845
    #   SalI   (GTCGAC) at  220
    # Filler regions are non-repetitive 50-mers seeded with a fixed
    # PRNG so tests stay deterministic.
    import random
    rng = random.Random(0xCAFEBABE)
    bases = "ACGT"
    def filler(n: int) -> str:
        return "".join(rng.choice(bases) for _ in range(n))

    parts = []
    parts.append(filler(220))                # 0..220
    parts.append("GTCGAC")                   # 220 SalI
    parts.append(filler(170))                # 226..396
    parts.append("GAATTC")                   # 396 EcoRI
    parts.append(filler(206))                # 402..608
    parts.append("GGATCC")                   # 608 BamHI
    parts.append(filler(488))                # 614..1102
    parts.append("AAGCTT")                   # 1102 HindIII
    parts.append(filler(737))                # 1108..1845
    parts.append("CTCGAG")                   # 1845 XhoI
    parts.append(filler(129))                # 1851..1980
    parts.append("GAATTC")                   # 1980 EcoRI
    parts.append(filler(700))                # 1986..2686
    seq_str = "".join(parts)
    assert len(seq_str) == 2686

    rec = SeqRecord(
        Seq(seq_str),
        id="SYNREAL", name="SYNREAL",
        description="Synthetic realistic 2.7 kb test plasmid",
    )
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"]      = "circular"

    # Resistance-marker-like CDS, 840 bp on the (+) strand.
    rec.features.append(SeqFeature(
        FeatureLocation(400, 1240, strand=1),
        type="CDS",
        qualifiers={"label": ["AmpR"], "gene": ["bla"],
                     "product": ["beta-lactamase (synthetic)"]},
    ))
    # Promoter upstream of the CDS.
    rec.features.append(SeqFeature(
        FeatureLocation(330, 400, strand=1),
        type="promoter",
        qualifiers={"label": ["AmpR_promoter"]},
    ))
    # RBS just before the CDS.
    rec.features.append(SeqFeature(
        FeatureLocation(390, 400, strand=1),
        type="RBS",
        qualifiers={"label": ["RBS"]},
    ))
    # Ori on the (-) strand, opposite the CDS.
    rec.features.append(SeqFeature(
        FeatureLocation(1700, 2400, strand=-1),
        type="rep_origin",
        qualifiers={"label": ["pMB1_ori"]},
    ))
    # MCS / polylinker spanning the cluster of restriction sites.
    rec.features.append(SeqFeature(
        FeatureLocation(216, 614, strand=0),
        type="misc_feature",
        qualifiers={"label": ["MCS"]},
    ))
    return rec


@pytest.fixture
def realistic_gb_path(tmp_path, realistic_plasmid):
    """Write `realistic_plasmid` to a .gb file and return the path."""
    from Bio import SeqIO
    p = tmp_path / "realistic.gb"
    SeqIO.write(realistic_plasmid, str(p), "genbank")
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


@pytest.fixture
def isolated_parts_bin(tmp_path, monkeypatch):
    """Redirect `_PARTS_BIN_FILE` to a tmp path. Same redundant-but-
    harmless pattern as `isolated_library` — `_protect_user_data` already
    handles redirect; this fixture exists so tests that want to read /
    inspect the redirected path can request it explicitly."""
    import splicecraft as sc
    tmp_bin = tmp_path / "parts_bin.json"
    monkeypatch.setattr(sc, "_PARTS_BIN_FILE", tmp_bin)
    monkeypatch.setattr(sc, "_parts_bin_cache", None)
    return tmp_bin


@pytest.fixture
def isolated_primers(tmp_path, monkeypatch):
    """Redirect `_PRIMERS_FILE` to a tmp path. Same redundant-but-harmless
    pattern as `isolated_library`."""
    import splicecraft as sc
    tmp_p = tmp_path / "primers.json"
    monkeypatch.setattr(sc, "_PRIMERS_FILE", tmp_p)
    monkeypatch.setattr(sc, "_primers_cache", None)
    return tmp_p
