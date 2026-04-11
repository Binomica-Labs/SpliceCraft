"""
test_plannotate — pLannotate integration tests.

pLannotate is GPL-3 and conda-only, so we NEVER import it or require it to be
installed in the test env. Every test here monkeypatches the things that
would normally cross the boundary:
  - `shutil.which` for availability detection
  - `subprocess.run` for the CLI invocation itself

The goal is to cover:
  (a) availability detection for every failure mode
  (b) the size-cap preflight
  (c) subprocess error handling (timeout, missing DB, non-zero exit, missing
      output file, parse failure)
  (d) feature merging — preserves originals, appends pLannotate hits tagged
      with note="pLannotate", skips duplicates
  (e) the cached status helper can be cleared between tests

A single integration-style test (`test_run_plannotate_success_full_path`)
exercises the real `_run_plannotate` happy path with a mocked subprocess that
"writes" a canned GenBank file, proving the parse + tmpdir plumbing works
end-to-end without needing the actual binary.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def clear_plannotate_cache(monkeypatch):
    """Reset the cached availability probe before each test so monkeypatched
    `shutil.which` actually takes effect."""
    monkeypatch.setattr(sc, "_PLANNOTATE_CHECK_CACHE", None)


def _make_which_stub(available: "set[str]"):
    """Return a drop-in replacement for shutil.which that reports the given
    executables as present (with a fake path) and everything else as missing."""
    def _which(cmd, *args, **kwargs):
        return f"/usr/bin/{cmd}" if cmd in available else None
    return _which


@pytest.fixture
def fake_annotated_record(tiny_record):
    """A second SeqRecord that simulates pLannotate's output: same sequence,
    but with two new features added on top of the existing ones."""
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from copy import deepcopy

    out = SeqRecord(
        Seq(str(tiny_record.seq)),
        id=tiny_record.id,
        name=tiny_record.name,
        description=tiny_record.description,
        annotations=dict(tiny_record.annotations),
    )
    for f in tiny_record.features:
        out.features.append(deepcopy(f))
    # Fake pLannotate hits — a promoter and a rep_origin the user didn't have
    out.features.append(SeqFeature(
        FeatureLocation(10, 40, strand=1),
        type="promoter",
        qualifiers={"label": ["lac promoter"]},
    ))
    out.features.append(SeqFeature(
        FeatureLocation(60, 100, strand=-1),
        type="rep_origin",
        qualifiers={"label": ["pMB1 ori"]},
    ))
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Availability detection
# ═══════════════════════════════════════════════════════════════════════════════

class TestAvailabilityDetection:
    def test_all_missing(self, monkeypatch, clear_plannotate_cache):
        import shutil
        monkeypatch.setattr(shutil, "which", _make_which_stub(set()))
        status = sc._plannotate_status()
        assert status["installed"] is False
        assert status["blast"]     is False
        assert status["diamond"]   is False
        assert status["ready"]     is False

    def test_plannotate_only(self, monkeypatch, clear_plannotate_cache):
        import shutil
        monkeypatch.setattr(shutil, "which", _make_which_stub({"plannotate"}))
        status = sc._plannotate_status()
        assert status["installed"] is True
        assert status["blast"]     is False
        assert status["ready"]     is False

    def test_blast_and_diamond_but_no_plannotate(self, monkeypatch,
                                                  clear_plannotate_cache):
        import shutil
        monkeypatch.setattr(
            shutil, "which", _make_which_stub({"blastn", "diamond"})
        )
        status = sc._plannotate_status()
        assert status["installed"] is False
        assert status["blast"]     is True
        assert status["diamond"]   is True
        assert status["ready"]     is False

    def test_all_present(self, monkeypatch, clear_plannotate_cache):
        import shutil
        monkeypatch.setattr(
            shutil, "which",
            _make_which_stub({"plannotate", "blastn", "diamond"}),
        )
        status = sc._plannotate_status()
        assert status["ready"] is True

    def test_status_is_cached(self, monkeypatch, clear_plannotate_cache):
        """Second call must not re-invoke shutil.which."""
        import shutil
        calls = {"n": 0}
        orig_which = _make_which_stub(set())
        def _counting(cmd, *a, **k):
            calls["n"] += 1
            return orig_which(cmd)
        monkeypatch.setattr(shutil, "which", _counting)
        sc._plannotate_status()
        first = calls["n"]
        sc._plannotate_status()
        assert calls["n"] == first   # no new which calls on second access

    def test_install_hint_mentions_conda_and_setupdb(self):
        hint = sc._plannotate_install_hint()
        assert "conda" in hint
        assert "setupdb" in hint


# ═══════════════════════════════════════════════════════════════════════════════
# Size cap preflight
# ═══════════════════════════════════════════════════════════════════════════════

class TestSizeCapPreflight:
    def test_rejects_over_50kb(self, monkeypatch, clear_plannotate_cache):
        import shutil
        monkeypatch.setattr(
            shutil, "which",
            _make_which_stub({"plannotate", "blastn", "diamond"}),
        )
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        too_big = SeqRecord(Seq("A" * 50_001), id="big", name="big",
                            description="")
        too_big.annotations["molecule_type"] = "DNA"
        with pytest.raises(sc.PlannotateTooLarge):
            sc._run_plannotate(too_big)

    def test_max_constant_is_50kb(self):
        assert sc._PLANNOTATE_MAX_BP == 50_000


# ═══════════════════════════════════════════════════════════════════════════════
# Feature merging (pure unit tests — no subprocess)
# ═══════════════════════════════════════════════════════════════════════════════

class TestMergeFeatures:
    def test_preserves_original_features(self, tiny_record, fake_annotated_record):
        merged = sc._merge_plannotate_features(tiny_record, fake_annotated_record)
        original_types = sorted(
            (f.type, int(f.location.start), int(f.location.end))
            for f in tiny_record.features
        )
        merged_types_subset = sorted(
            (f.type, int(f.location.start), int(f.location.end))
            for f in merged.features
            if f.type in {t for t, _, _ in original_types}
        )
        # Every original feature must still appear in the merged record
        for ot in original_types:
            assert ot in merged_types_subset, f"lost original feature {ot}"

    def test_adds_plannotate_features(self, tiny_record, fake_annotated_record):
        merged = sc._merge_plannotate_features(tiny_record, fake_annotated_record)
        assert getattr(merged, "_plannotate_added") == 2
        labels = [
            f.qualifiers.get("label", [""])[0]
            for f in merged.features if f.type in ("promoter", "rep_origin")
        ]
        assert "lac promoter" in labels
        assert "pMB1 ori" in labels

    def test_plannotate_features_tagged_with_note(self, tiny_record,
                                                   fake_annotated_record):
        merged = sc._merge_plannotate_features(tiny_record, fake_annotated_record)
        for f in merged.features:
            if f.type in ("promoter", "rep_origin"):
                notes = f.qualifiers.get("note", [])
                assert any("pLannotate" in n for n in notes), (
                    f"feature {f.type} lacks 'pLannotate' note"
                )

    def test_skips_duplicate_features(self, tiny_record):
        """A pLannotate hit at the exact same type+start+end+strand as an
        existing feature must NOT be duplicated into the merged record."""
        from copy import deepcopy
        # annotated has exactly the same CDS feature the original already has
        annotated = deepcopy(tiny_record)
        merged = sc._merge_plannotate_features(tiny_record, annotated)
        assert merged._plannotate_added == 0
        # Same number of features as the original (ignoring 'source' which
        # neither record has in the fixture)
        assert len(merged.features) == len(tiny_record.features)

    def test_original_features_keep_their_notes(self, tiny_record,
                                                 fake_annotated_record):
        """Original features must not have 'pLannotate' slapped on them —
        only NEW features get that tag."""
        merged = sc._merge_plannotate_features(tiny_record, fake_annotated_record)
        for f in merged.features:
            if f.type in ("CDS", "misc_feature"):
                notes = f.qualifiers.get("note", [])
                notes = notes if isinstance(notes, list) else [notes]
                assert not any("pLannotate" in n for n in notes), (
                    f"original feature {f.type} wrongly tagged as pLannotate"
                )

    def test_sequence_bytes_preserved(self, tiny_record, fake_annotated_record):
        merged = sc._merge_plannotate_features(tiny_record, fake_annotated_record)
        assert str(merged.seq) == str(tiny_record.seq)

    def test_returns_new_record_not_mutated_original(self, tiny_record,
                                                      fake_annotated_record):
        n_original_features_before = len(tiny_record.features)
        sc._merge_plannotate_features(tiny_record, fake_annotated_record)
        assert len(tiny_record.features) == n_original_features_before


# ═══════════════════════════════════════════════════════════════════════════════
# _run_plannotate error paths (subprocess mocked)
# ═══════════════════════════════════════════════════════════════════════════════

class _FakeCompleted:
    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout     = stdout
        self.stderr     = stderr


class TestRunPlannotateErrors:
    def test_not_installed_raises(self, tiny_record, monkeypatch,
                                   clear_plannotate_cache):
        import shutil
        monkeypatch.setattr(shutil, "which", _make_which_stub(set()))
        with pytest.raises(sc.PlannotateNotInstalled) as ei:
            sc._run_plannotate(tiny_record)
        assert "not installed" in ei.value.user_msg.lower()

    def test_missing_blast_raises_not_installed(self, tiny_record, monkeypatch,
                                                 clear_plannotate_cache):
        import shutil
        monkeypatch.setattr(
            shutil, "which", _make_which_stub({"plannotate", "diamond"})
        )
        with pytest.raises(sc.PlannotateNotInstalled) as ei:
            sc._run_plannotate(tiny_record)
        assert "blastn" in ei.value.user_msg

    def test_missing_diamond_raises_not_installed(self, tiny_record, monkeypatch,
                                                   clear_plannotate_cache):
        import shutil
        monkeypatch.setattr(
            shutil, "which", _make_which_stub({"plannotate", "blastn"})
        )
        with pytest.raises(sc.PlannotateNotInstalled) as ei:
            sc._run_plannotate(tiny_record)
        assert "diamond" in ei.value.user_msg

    def test_database_missing_raises_missing_db(self, tiny_record, monkeypatch,
                                                 clear_plannotate_cache):
        """pLannotate signals missing DBs via 'Databases not downloaded' on
        stdout and exit code 0 — the most insidious failure mode."""
        import shutil, subprocess
        monkeypatch.setattr(
            shutil, "which",
            _make_which_stub({"plannotate", "blastn", "diamond"}),
        )
        def _fake_run(*args, **kwargs):
            return _FakeCompleted(
                returncode=0,
                stdout="Databases not downloaded. Run 'plannotate setupdb'...",
            )
        monkeypatch.setattr(subprocess, "run", _fake_run)
        with pytest.raises(sc.PlannotateMissingDb) as ei:
            sc._run_plannotate(tiny_record)
        assert "setupdb" in ei.value.detail or "setupdb" in ei.value.user_msg

    def test_nonzero_exit_raises_failed(self, tiny_record, monkeypatch,
                                         clear_plannotate_cache):
        import shutil, subprocess
        monkeypatch.setattr(
            shutil, "which",
            _make_which_stub({"plannotate", "blastn", "diamond"}),
        )
        def _fake_run(*args, **kwargs):
            return _FakeCompleted(
                returncode=1, stdout="", stderr="diamond: no such file",
            )
        monkeypatch.setattr(subprocess, "run", _fake_run)
        with pytest.raises(sc.PlannotateFailed) as ei:
            sc._run_plannotate(tiny_record)
        assert "failed" in ei.value.user_msg.lower()

    def test_timeout_raises_failed(self, tiny_record, monkeypatch,
                                    clear_plannotate_cache):
        import shutil, subprocess
        monkeypatch.setattr(
            shutil, "which",
            _make_which_stub({"plannotate", "blastn", "diamond"}),
        )
        def _fake_run(*args, **kwargs):
            raise subprocess.TimeoutExpired(cmd="plannotate", timeout=180)
        monkeypatch.setattr(subprocess, "run", _fake_run)
        with pytest.raises(sc.PlannotateFailed) as ei:
            sc._run_plannotate(tiny_record, timeout=180)
        assert "timed out" in ei.value.user_msg.lower()

    def test_file_not_found_at_exec_raises_not_installed(self, tiny_record,
                                                          monkeypatch,
                                                          clear_plannotate_cache):
        """PATH said plannotate existed but it disappeared before exec —
        report as not-installed rather than crash."""
        import shutil, subprocess
        monkeypatch.setattr(
            shutil, "which",
            _make_which_stub({"plannotate", "blastn", "diamond"}),
        )
        def _fake_run(*args, **kwargs):
            raise FileNotFoundError("no plannotate")
        monkeypatch.setattr(subprocess, "run", _fake_run)
        with pytest.raises(sc.PlannotateNotInstalled):
            sc._run_plannotate(tiny_record)

    def test_no_output_file_raises_failed(self, tiny_record, monkeypatch,
                                           clear_plannotate_cache):
        import shutil, subprocess
        monkeypatch.setattr(
            shutil, "which",
            _make_which_stub({"plannotate", "blastn", "diamond"}),
        )
        def _fake_run(*args, **kwargs):
            return _FakeCompleted(returncode=0, stdout="done", stderr="")
        monkeypatch.setattr(subprocess, "run", _fake_run)
        with pytest.raises(sc.PlannotateFailed) as ei:
            sc._run_plannotate(tiny_record)
        assert "no .gbk output" in ei.value.user_msg


# ═══════════════════════════════════════════════════════════════════════════════
# Happy-path end-to-end (subprocess writes a real .gbk to the real tmpdir)
# ═══════════════════════════════════════════════════════════════════════════════

class TestRunPlannotateHappyPath:
    def test_successful_run_returns_parsed_record(
        self, tiny_record, fake_annotated_record, monkeypatch,
        clear_plannotate_cache,
    ):
        """Fake subprocess 'writes' a GenBank file to -o <tmpdir>/annotated.gbk
        and exits 0 — verifies the tmpfile plumbing, the output-file locator,
        and the parse path all work correctly without the real binary."""
        import shutil, subprocess
        monkeypatch.setattr(
            shutil, "which",
            _make_which_stub({"plannotate", "blastn", "diamond"}),
        )

        def _fake_run(cmd, *args, **kwargs):
            # Pull the -o <tmpdir> value out of the cmd list
            out_dir = cmd[cmd.index("-o") + 1]
            file_name = cmd[cmd.index("-f") + 1]
            # Write our fake annotated record as the .gbk output
            out_path = Path(out_dir) / f"{file_name}.gbk"
            out_path.write_text(sc._record_to_gb_text(fake_annotated_record))
            return _FakeCompleted(returncode=0, stdout="ok", stderr="")

        monkeypatch.setattr(subprocess, "run", _fake_run)

        annotated = sc._run_plannotate(tiny_record)
        # The parsed record should include features (at least the ones our
        # fake annotated_record put in)
        assert len(annotated.features) >= len(tiny_record.features)
        # Sequence bytes preserved
        assert str(annotated.seq) == str(tiny_record.seq)
