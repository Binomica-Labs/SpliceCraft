"""Sub-consensus read heterogeneity, and the read-only agent attach.

Two things from the same field report, kept together because the second is
what makes the first reachable while the GUI is open.

`analyse-read-heterogeneity` answers a question neither `read-consensus` nor
`verify-against-reads` can: **is what I sequenced ONE thing?** A consensus is a
single read, so depth is 1 everywhere and anything below it is invisible — a
population of escapers each carrying a different inactivating mutation, none
dominant, consenses back to wild type and reads as a clean clonal culture.
Only an allele FRACTION over raw reads shows it.

The acceptance case is therefore built, not asserted from the implementation:
a synthetic FASTQ of 90% reference plus 10% of a known SNP must report that
position at fraction ~0.10, and a pile of eight different low-frequency
mutations must read as `mixed` rather than `clonal`.
"""

from __future__ import annotations

import random

import pytest

import splicecraft as sc


def _rand(n: int, seed: int) -> str:
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


_TRANSVERSION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def _write_fastq(tmp_path, name, reads, phred=40):
    """A FASTQ with a flat quality of `phred` over every base."""
    ch = chr(33 + phred)
    p = tmp_path / name
    p.write_text("".join(
        f"@read{i}\n{s}\n+\n{ch * len(s)}\n" for i, s in enumerate(reads)))
    return str(p)


def _call(name, payload, app=None):
    return sc._state._AGENT_HANDLERS[name][0](app, payload)


REF = _rand(1200, 7)
SNP_POS = 400
ALT = _TRANSVERSION[REF[SNP_POS]]
MUT = REF[:SNP_POS] + ALT + REF[SNP_POS + 1:]


class TestAlleleFractions:

    def test_a_90_10_mixture_reports_fraction_one_tenth(self, tmp_path):
        """The acceptance case from the field report."""
        reads = [MUT] * 10 + [REF] * 90
        fq = _write_fastq(tmp_path, "panel.fastq", reads)
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": fq, "circular": False})
        assert not isinstance(res, tuple), res
        assert res["verdict"] == "minor_variants"
        assert res["n_reads"] == 100
        hit = [p for p in res["positions"] if p["pos"] == SNP_POS]
        assert hit, f"SNP missed; got {[p['pos'] for p in res['positions']]}"
        p = hit[0]
        assert p["fraction"] == pytest.approx(0.10, abs=0.02)
        assert p["alt"] == ALT
        assert p["ref"] == REF[SNP_POS]
        assert p["alt_reads"] == 10
        assert p["depth"] == 100
        assert p["contradicted"] == 90
        assert p["mean_phred"] == pytest.approx(40.0)

    def test_a_clonal_culture_reads_as_clonal(self, tmp_path):
        fq = _write_fastq(tmp_path, "clonal.fastq", [REF] * 20)
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": fq, "circular": False})
        assert res["verdict"] == "clonal"
        assert res["positions"] == []
        assert res["mean_depth"] == pytest.approx(20.0)

    def test_an_escaper_population_reads_as_mixed(self, tmp_path):
        """THE case the whole endpoint exists for: eight different mutations,
        each at 5%, none dominant. A consensus of these reads is wild type."""
        reads = []
        for i in range(40):
            k = 100 + (i % 8) * 50
            b = _TRANSVERSION[REF[k]]
            reads.append(REF[:k] + b + REF[k + 1:])
        reads += [REF] * 60
        fq = _write_fastq(tmp_path, "escapers.fastq", reads)
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": fq, "circular": False})
        assert res["verdict"] == "mixed"
        assert res["n_positions"] >= 5
        assert max(p["fraction"] for p in res["positions"]) < 0.25

    def test_a_dominant_variant_reads_as_mixed_not_minor(self, tmp_path):
        reads = [MUT] * 50 + [REF] * 50
        fq = _write_fastq(tmp_path, "half.fastq", reads)
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": fq, "circular": False})
        assert res["verdict"] == "mixed"

    def test_min_fraction_filters_below_the_floor(self, tmp_path):
        reads = [MUT] * 2 + [REF] * 98      # 2%
        fq = _write_fastq(tmp_path, "rare.fastq", reads)
        seen = _call("analyse-read-heterogeneity",
                     {"reference": REF, "reads_path": fq, "circular": False,
                      "min_fraction": 0.01})
        assert any(p["pos"] == SNP_POS for p in seen["positions"])
        hidden = _call("analyse-read-heterogeneity",
                       {"reference": REF, "reads_path": fq,
                        "circular": False, "min_fraction": 0.05})
        assert not any(p["pos"] == SNP_POS for p in hidden["positions"])
        assert hidden["verdict"] == "clonal"

    def test_low_quality_basecalls_are_excluded_from_the_fraction(
            self, tmp_path):
        """At 1% the error floor of every instrument looks like a
        sub-population. A variant sitting in unreliable basecalls must not be
        counted as evidence of one."""
        reads = [MUT] * 10 + [REF] * 90
        fq = _write_fastq(tmp_path, "noisy.fastq", reads, phred=5)
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": fq, "circular": False,
                     "min_phred": 20})
        assert res["verdict"] == "clonal"
        assert res["n_filtered_low_quality"] >= 10
        # Lower the bar and the same reads produce the same variant again.
        kept = _call("analyse-read-heterogeneity",
                     {"reference": REF, "reads_path": fq, "circular": False,
                      "min_phred": 0})
        assert any(p["pos"] == SNP_POS for p in kept["positions"])

    def test_a_no_call_is_not_an_allele(self, tmp_path):
        """HARDENING REGRESSION. A basecaller emitting `N` is saying it could
        not read that base. Counting it as a difference turned ONE ordinary
        read carrying a ten-base N run into a ten-position "10%
        sub-population" and the whole sample into `mixed` — and real long-read
        data carries N runs routinely, so this fired on almost every genuine
        dataset."""
        reads = [REF] * 9 + [REF[:150] + "N" * 10 + REF[160:]]
        fq = _write_fastq(tmp_path, "ncalls.fastq", reads)
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": fq, "circular": False})
        assert res["verdict"] == "clonal"
        assert res["positions"] == []
        # Counted, not silently dropped — ambiguous coverage is a fact about
        # the data the caller should be able to see.
        assert res["n_no_call_observations"] >= 10

    def test_a_real_variant_is_still_found_beside_no_calls(self, tmp_path):
        """The other half: over-filtering would hide the thing being looked
        for. The N run and the SNP are in different reads and at different
        positions, and only the SNP is an allele."""
        reads = ([MUT] * 10
                 + [REF[:200] + "N" * 8 + REF[208:]] * 5
                 + [REF] * 85)
        fq = _write_fastq(tmp_path, "mixed_ncalls.fastq", reads)
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": fq, "circular": False})
        hit = [p for p in res["positions"] if p["pos"] == SNP_POS]
        assert hit, f"the real SNP was filtered away: {res['positions']}"
        assert hit[0]["fraction"] == pytest.approx(0.10, abs=0.02)
        assert all(p["alt"] != "N" for p in res["positions"])
        assert res["n_no_call_observations"] >= 8

    def test_the_thresholds_that_produced_the_verdict_are_reported(
            self, tmp_path):
        """They are conventions, not a calibrated model. A caller re-deriving
        its own call from `positions` must be able to see exactly what these
        were."""
        fq = _write_fastq(tmp_path, "t.fastq", [REF] * 5)
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": fq, "circular": False})
        t = res["thresholds"]
        assert set(t) == {"min_fraction", "min_phred", "mixed_fraction",
                          "mixed_min_positions"}

    def test_uncovered_spans_are_part_of_the_answer(self, tmp_path):
        """A plasmid can be clonal over the half that was read."""
        half = REF[:600]
        fq = _write_fastq(tmp_path, "half.fastq", [half] * 10)
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": fq, "circular": False})
        assert res["uncovered_spans"]
        assert res["covered_pct"] < 60

    def test_inline_reads_with_parallel_quality(self):
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads": [MUT] + [REF] * 9,
                     "read_quality": [[40] * len(MUT)] + [None] * 9,
                     "circular": False})
        assert res["reads_source"] == "inline"
        assert res["n_reads"] == 10
        hit = [p for p in res["positions"] if p["pos"] == SNP_POS]
        assert hit and hit[0]["fraction"] == pytest.approx(0.1, abs=0.01)


class TestHeterogeneityInput:

    @pytest.mark.parametrize("body,code", [
        ({"reference": REF}, 400),                            # no reads
        ({"reads": ["ACGT"]}, 400),                           # no reference
        ({"reference": REF, "reads": [], "circular": False}, 400),
        ({"reference": REF, "reads": ["ACGT"], "min_fraction": 0}, 400),
        ({"reference": REF, "reads": ["ACGT"], "min_fraction": 2}, 400),
        ({"reference": REF, "reads": ["ACGT"], "min_phred": 200}, 400),
        ({"reference": REF, "reads": ["ACGT"], "max_reads": 0}, 400),
        ({"reference": REF, "reads": ["ACGT"], "mode": "sideways"}, 400),
        ({"reference": "", "reads": ["ACGT"]}, 400),
    ])
    def test_bad_input_is_refused_with_a_reason(self, body, code):
        res = _call("analyse-read-heterogeneity", body)
        assert isinstance(res, tuple) and res[1] == code, res

    def test_both_reads_and_reads_path_is_refused(self, tmp_path):
        fq = _write_fastq(tmp_path, "x.fastq", [REF])
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads": [REF], "reads_path": fq})
        assert isinstance(res, tuple) and res[1] == 400

    def test_a_missing_fastq_is_404(self):
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": "/nope/missing.fastq"})
        assert isinstance(res, tuple) and res[1] == 404

    def test_a_non_fastq_path_says_why_quality_is_needed(self, tmp_path):
        p = tmp_path / "reads.fasta"
        p.write_text(">r\nACGT\n")
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads_path": str(p)})
        assert isinstance(res, tuple) and res[1] == 400
        assert "quality" in res[0]["error"]

    def test_unknown_reference_name_is_404(self):
        res = _call("analyse-read-heterogeneity",
                    {"reference_name": "no such plasmid", "reads": [REF]})
        assert isinstance(res, tuple) and res[1] == 404

    def test_read_quality_must_be_parallel_to_reads(self):
        res = _call("analyse-read-heterogeneity",
                    {"reference": REF, "reads": [REF, REF],
                     "read_quality": [[40] * len(REF)]})
        assert isinstance(res, tuple) and res[1] == 400
        assert "parallel" in res[0]["error"]


class TestMeanDepth:
    """`mean_depth` lives in `_multi_read_summary`, beside the depth array it
    is derived from, so `read-consensus` reports the same number."""

    def test_mean_depth_is_over_covered_bases_not_the_whole_molecule(self):
        aligns = []
        for _ in range(4):
            aligns.append({"axis": "target", "result": {
                "aligned_q": REF[:300], "aligned_t": REF[:300],
                "t_start": 0, "t_end": 300, "identity_pct": 100.0}})
        out = sc._multi_read_summary(aligns, len(REF))
        # Four reads over the same 300 bp: depth 4 there, 0 elsewhere.
        # Averaging over the whole 1200 bp would report 1.0.
        assert out["mean_depth"] == pytest.approx(4.0, abs=0.01)
        assert out["max_depth"] == 4

    def test_no_reads_reports_zero_not_a_crash(self):
        out = sc._multi_read_summary([], 500)
        assert out["mean_depth"] == 0.0
        assert out["max_depth"] == 0


class TestReadOnlyAttach:
    """`splicecraft --agent --read-only` — attach while the GUI holds the
    data-dir lock, serve every read, refuse every write.

    Two independent layers, and the test exercises both: the 409 at
    `_agent_invoke` makes the refusal legible, and the ABSENCE of
    `_authorize_writes` is what actually protects the data.
    """

    @pytest.fixture
    def read_only(self, monkeypatch):
        monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY", True)
        monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY_HOLDER_PID", 4242)
        yield

    def test_every_write_endpoint_is_refused(self, read_only):
        writes = [n for n, (_fn, w) in sc._state._AGENT_HANDLERS.items() if w]
        assert writes, "no write endpoints registered — test is vacuous"
        for name in writes:
            payload, status = sc._agent_invoke(None, name, {})
            assert status == 409, f"{name} returned {status}"
            assert payload["read_only"] is True
            assert payload["holder_pid"] == 4242

    def test_the_refusal_names_the_holding_process(self, read_only):
        writes = [n for n, (_fn, w) in sc._state._AGENT_HANDLERS.items() if w]
        payload, _ = sc._agent_invoke(None, writes[0], {})
        assert "4242" in payload["error"]
        assert "--read-only" in payload["error"]

    def test_reads_still_answer(self, read_only):
        payload, status = sc._agent_invoke(None, "tools", {})
        assert status == 200
        assert payload

    def test_the_new_endpoints_are_all_reachable_read_only(self, read_only):
        """They are the whole reason the mode was asked for."""
        for name in ("scan-promoters", "scan-terminators",
                     "map-transcription", "analyse-cds",
                     "analyse-read-heterogeneity", "read-consensus",
                     "verify-against-reads"):
            _payload, status = sc._agent_invoke(None, name, {})
            assert status != 409, f"{name} was refused as a write"

    def test_writes_are_unaffected_when_not_read_only(self):
        """The gate must be off by default — a mode that silently stayed on
        would break every ordinary session."""
        assert sc._state._AGENT_READ_ONLY is False
        writes = [n for n, (_fn, w) in sc._state._AGENT_HANDLERS.items() if w]
        payload, status = sc._agent_invoke(None, writes[0], {})
        assert status != 409 or not payload.get("read_only")


class TestReadOnlyStartupWrites:
    """Launch-time writes must be skipped, not merely refused.

    REGRESSION: the first working build attached fine in every unit test and
    crashed on a real launch. `PlasmidApp.compose` runs the collection /
    primer / parts-bin / project migrations before any child mounts, and each
    one either creates a default or REWRITES a live file from its
    active-collection mirror. The L2 chokepoint correctly refused the write —
    and took the whole launch down with it. Only attaching to a genuinely
    locked data dir reached that path.
    """

    _MIGRATIONS = (
        "_ensure_default_collection",
        "_restore_library_from_active_collection",
        "_ensure_default_primer_collection",
        "_restore_primers_from_active_primer_collection",
        "_ensure_default_parts_bin",
        "_restore_parts_bin_from_active_bin",
        "_ensure_default_project",
        "_restore_experiments_from_active_project",
    )

    def test_compose_guards_every_startup_migration(self):
        import inspect
        src = inspect.getsource(sc.PlasmidApp.compose)
        assert "_state._AGENT_READ_ONLY" in src, (
            "compose() no longer guards its startup migrations — a "
            "--read-only attach will crash on the first mirror write")
        guard_at = src.index("_state._AGENT_READ_ONLY")
        for name in self._MIGRATIONS:
            at = src.find(name + "()")
            assert at > guard_at, f"{name} runs before the read-only guard"

    @pytest.mark.parametrize("name", _MIGRATIONS)
    def test_each_migration_really_would_write(self, name, monkeypatch):
        """Keeps the guard honest: if one of these stops writing, the guard
        entry is stale, and if a NEW writing migration is added without being
        guarded this test's list is where the omission shows."""
        assert hasattr(sc, name), f"{name} no longer exists — update the list"
        monkeypatch.setattr(sc._state, "_SAVES_AUTHORIZED", False)
        try:
            getattr(sc, name)()
        except RuntimeError as exc:
            assert "not authorised" in str(exc)
        except Exception:
            # Some migrations short-circuit before writing on an already
            # migrated data dir; that is fine. What must never happen is a
            # SILENT write, which the chokepoint makes impossible.
            pass


class TestLockHolderRead:
    """The lock-holder read must not touch the lockfile — the mode's whole
    promise is that attaching changes nothing on disk."""

    def test_reads_the_pid_without_modifying_the_file(self, tmp_path):
        lf = tmp_path / "splicecraft.lock"
        lf.write_text("9999\n1.2.66\n")
        before = lf.read_bytes()
        assert sc._data_dir_lock_holder_pid(tmp_path) == 9999
        assert lf.read_bytes() == before

    def test_a_missing_lockfile_is_zero_not_an_error(self, tmp_path):
        assert sc._data_dir_lock_holder_pid(tmp_path) == 0

    @pytest.mark.parametrize("body", ["", "garbage\n", "-1\n", "0\n", "\n\n"])
    def test_a_malformed_lockfile_is_zero(self, tmp_path, body):
        (tmp_path / "splicecraft.lock").write_text(body)
        assert sc._data_dir_lock_holder_pid(tmp_path) == 0

    def test_it_does_not_create_a_lockfile(self, tmp_path):
        sc._data_dir_lock_holder_pid(tmp_path)
        assert not (tmp_path / "splicecraft.lock").exists()


class TestReadOnlyTokenFile:
    """A read-only attach is a GUEST on someone else's data dir."""

    def test_the_guest_token_filename_is_separate(self):
        assert sc._AGENT_TOKEN_READONLY_NAME != "agent_token"

    def test_the_cli_prefers_the_host_and_falls_back_to_the_guest(
            self, tmp_path, monkeypatch):
        import splicecraft_cli as cli
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
        monkeypatch.delenv("SPLICECRAFT_READ_ONLY", raising=False)
        # Neither present -> the host path is what the error names.
        assert cli._token_file().name == cli.TOKEN_FILENAME
        # Only the guest -> use it, so `splicecraft-cli status` just works.
        (tmp_path / cli.READONLY_TOKEN_FILENAME).write_text("6702\ntok\n")
        assert cli._token_file().name == cli.READONLY_TOKEN_FILENAME
        # Host appears -> it wins, because it can write.
        (tmp_path / cli.TOKEN_FILENAME).write_text("6701\ntok\n")
        assert cli._token_file().name == cli.TOKEN_FILENAME
        # …unless the guest is pinned explicitly.
        monkeypatch.setenv("SPLICECRAFT_READ_ONLY", "1")
        assert cli._token_file().name == cli.READONLY_TOKEN_FILENAME
