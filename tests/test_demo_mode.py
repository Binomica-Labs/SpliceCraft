"""Demo mode (2026-06-13): the `$SPLICECRAFT_DEMO` sandbox for the web demo.

The load-bearing safety property: a demo session must NEVER resolve its data
dir to a real library (it forces an ephemeral tempdir and ignores
`$SPLICECRAFT_DATA_DIR` / `$XDG`). Plus the web tier's defence-in-depth
lockdown (refuse file-open / NCBI / agent API, with a polite excuse).
"""
from __future__ import annotations

import subprocess
import sys

import splicecraft as sc


class _FakeApp:
    def __init__(self):
        self.notes: list = []

    def notify(self, msg, **kw):
        self.notes.append((msg, kw))


class TestDemoModeDetection:
    def test_resolve_demo_mode(self, monkeypatch):
        for v in ("local", "full", "1", "true", "yes", "on", "LOCAL"):
            monkeypatch.setenv("SPLICECRAFT_DEMO", v)
            assert sc._resolve_demo_mode() == "local", v
        for v in ("web", "public", "restrained", "kiosk", "WEB"):
            monkeypatch.setenv("SPLICECRAFT_DEMO", v)
            assert sc._resolve_demo_mode() == "web", v
        for v in ("", "  ", "nonsense", "0", "off"):
            monkeypatch.setenv("SPLICECRAFT_DEMO", v)
            assert sc._resolve_demo_mode() is None, repr(v)
        monkeypatch.delenv("SPLICECRAFT_DEMO", raising=False)
        assert sc._resolve_demo_mode() is None


class TestDemoDataDirSandbox:
    """THE safety boundary: demo mode forces an ephemeral tempdir and ignores
    every real-path env override, so a demo can't read/clobber a real library.
    Verified in a subprocess because `_DATA_DIR` is computed at import time."""

    def _probe(self, demo_value: str) -> str:
        # Point SPLICECRAFT_DATA_DIR at a (fake) "real" path; demo mode MUST
        # ignore it and resolve `_DATA_DIR` to a tempdir instead.
        code = (
            "import splicecraft as sc, tempfile, pathlib, sys\n"
            "dd = pathlib.Path(sc._state._DATA_DIR).resolve()\n"
            "tmp = pathlib.Path(tempfile.gettempdir()).resolve()\n"
            "assert sc._DEMO_MODE in ('local','web'), sc._DEMO_MODE\n"
            "assert dd.is_relative_to(tmp), ('not under tmp', dd)\n"
            "assert 'fake-real-library' not in str(dd), ('LEAKED real path', dd)\n"
            "assert f'splicecraft-demo-{sc._DEMO_MODE}' in dd.name, dd.name\n"
            "print(dd)\n"
        )
        env = {
            "SPLICECRAFT_DEMO": demo_value,
            "SPLICECRAFT_DATA_DIR": "/some/fake-real-library",
            "PATH": __import__("os").environ.get("PATH", ""),
            "HOME": __import__("os").environ.get("HOME", ""),
        }
        out = subprocess.run([sys.executable, "-c", code], env=env,
                             capture_output=True, text=True, timeout=120)
        assert out.returncode == 0, (
            f"demo={demo_value} sandbox check FAILED:\n{out.stderr}")
        return out.stdout.strip()

    def test_web_demo_ignores_data_dir_override(self):
        dd = self._probe("web")
        assert dd and "fake-real-library" not in dd

    def test_local_demo_ignores_data_dir_override(self):
        dd = self._probe("local")
        assert dd and "fake-real-library" not in dd


class TestWebDemoLockdown:
    def test_refuse_only_fires_in_web(self, monkeypatch):
        app = _FakeApp()
        monkeypatch.setattr(sc, "_DEMO_MODE", "web")
        assert sc._demo_web_refuse(app, "Fetching from NCBI") is True
        assert app.notes and "web demo" in app.notes[-1][0].lower()
        # local + off → never blocks, never notifies
        monkeypatch.setattr(sc, "_DEMO_MODE", "local")
        app2 = _FakeApp()
        assert sc._demo_web_refuse(app2, "Opening files") is False
        assert not app2.notes
        monkeypatch.setattr(sc, "_DEMO_MODE", None)
        assert sc._demo_web_refuse(_FakeApp(), "x") is False

    def test_agent_api_refused_in_web_demo(self, monkeypatch):
        monkeypatch.setattr(sc, "_DEMO_MODE", "web")
        assert sc._start_agent_api(_FakeApp()) is None

    def test_demo_helpers(self, monkeypatch):
        monkeypatch.setattr(sc, "_DEMO_MODE", "web")
        assert sc._demo_active() and sc._demo_is_web()
        monkeypatch.setattr(sc, "_DEMO_MODE", "local")
        assert sc._demo_active() and not sc._demo_is_web()
        monkeypatch.setattr(sc, "_DEMO_MODE", None)
        assert not sc._demo_active() and not sc._demo_is_web()

    def test_egress_blocked_in_web(self, monkeypatch):
        import urllib.error
        monkeypatch.setattr(sc, "_DEMO_MODE", "web")
        for fn, args in [
            (sc.fetch_genbank, ("L09137",)),
            (sc.fetch_protein, ("P12345",)),
            (sc._fetch_latest_pypi_version_ex, ()),
            (sc._ncbi_taxid_search, ("coli",)),
            (sc._codon_fetch_kazusa, ("83333",)),
            (sc._build_hardened_url_opener, ()),
        ]:
            with __import__("pytest").raises(urllib.error.URLError):
                fn(*args)
        # ...and NOT blocked outside web demo (helper is a no-op).
        monkeypatch.setattr(sc, "_DEMO_MODE", None)
        sc._demo_block_network("x")  # must not raise

    def test_master_delete_refused_in_web(self, monkeypatch):
        monkeypatch.setattr(sc, "_DEMO_MODE", "web")
        out = sc._perform_master_delete(
            _FakeApp(), sentinel=sc._MASTER_DELETE_SENTINEL)
        assert out["files_removed"] == 0 and out["dirs_removed"] == 0


class TestWebDemoComputeCap:
    async def test_apply_record_caps_oversize_in_web(self, monkeypatch):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        monkeypatch.setattr(sc, "_DEMO_MODE", "web")
        monkeypatch.setattr(sc, "_DEMO_WEB_MAX_BP", 200)
        big = SeqRecord(Seq("ACGT" * 100), id="big", name="big",   # 400 bp > 200
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 40)) as pilot:
            await pilot.pause()
            before = app._current_record
            app._apply_record(big)        # over cap → refused, no load
            await pilot.pause()
            assert app._current_record is before, "oversize record loaded in web demo"
            # under-cap record loads fine
            small = SeqRecord(Seq("ACGT" * 10), id="small", name="small",
                              annotations={"molecule_type": "DNA",
                                           "topology": "circular"})
            app._apply_record(small)
            await pilot.pause()
            assert app._current_record is not None
            assert "small" in str(getattr(app._current_record, "id", ""))


class TestDemoModeQA:
    """Capstone walkthrough: LOCAL is fully open; WEB locks exactly the
    high-risk surfaces (egress / host-FS / agent / destructive / oversize)
    while the pure-compute tools — Scrub, operon, translate, etc. — STILL
    work on the seed (that's the demo's whole point)."""

    def test_local_mode_blocks_nothing(self, monkeypatch):
        monkeypatch.setattr(sc, "_DEMO_MODE", "local")
        sc._demo_block_network("x")                          # no raise
        assert sc._demo_web_refuse(_FakeApp(), "open") is False
        assert sc._start_agent_api is not None               # not demo-gated
        # local master-delete is NOT short-circuited by the demo guard
        # (it still requires the sentinel; we don't actually run a wipe here).

    def test_web_locks_surfaces(self, monkeypatch):
        import urllib.error
        import pytest as _pt
        monkeypatch.setattr(sc, "_DEMO_MODE", "web")
        with _pt.raises(urllib.error.URLError):
            sc.fetch_genbank("L09137")
        assert sc._demo_web_refuse(_FakeApp(), "Opening files") is True
        assert sc._start_agent_api(_FakeApp()) is None
        assert sc._perform_master_delete(
            _FakeApp(), sentinel=sc._MASTER_DELETE_SENTINEL)["files_removed"] == 0

    async def test_open_file_action_refused_in_web(self, monkeypatch):
        # The host-FS browser never opens in the web demo — every push site
        # (Ctrl+O action, Parts Bin ▸ Open file, entry-vector ▸ Open file) is
        # gated, so the OpenFileModal is never pushed.
        monkeypatch.setattr(sc, "_DEMO_MODE", "web")
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 40)) as pilot:
            await pilot.pause()
            app.action_open_file()
            for _ in range(3):
                await pilot.pause()
            assert not isinstance(app.screen, sc.OpenFileModal), \
                "file dialog opened in the web demo"

    def test_web_tools_still_run_on_seed(self, monkeypatch):
        # The locked surfaces are network/FS/destructive — the SCIENCE tools are
        # pure compute and must keep working so the demo is actually usable.
        monkeypatch.setattr(sc, "_DEMO_MODE", "web")
        scrub = sc._make_demo_scrub_plasmid()
        cds = next(f for f in scrub.features if f.type == "CDS")
        plan = sc._scrub_gb_design(
            str(scrub.seq),
            [{"type": "CDS", "start": int(cds.location.start),
              "end": int(cds.location.end), "strand": 1}],
            ["BsaI", "Esp3I"], circular=True)
        assert plan["ok"] and plan["verified"], "Scrub broke in web demo"
        op = sc._make_demo_operon_plasmid()
        for f in [f for f in op.features if f.type == "CDS"]:
            prot = sc._translate_cds(str(op.seq), int(f.location.start),
                                     int(f.location.end), 1)
            assert prot.count("*") == 1, prot


class TestDemoSeed:
    async def test_seed_touches_widgets_only_from_the_ui_thread(
            self, monkeypatch):
        """The seed runs in a thread worker, and it repopulated the library
        table and scheduled the canvas load straight from that thread —
        DataTable surgery off the event loop, on the public web demo
        (audit 2026-10-23)."""
        import threading
        seen: "dict[str, set]" = {"repopulate": set(), "apply": set()}
        real_repop = sc.LibraryPanel._repopulate
        real_apply = sc.PlasmidApp._apply_record

        def repop(self, *a, **k):
            seen["repopulate"].add(threading.get_ident())
            return real_repop(self, *a, **k)

        def apply(self, *a, **k):
            seen["apply"].add(threading.get_ident())
            return real_apply(self, *a, **k)
        monkeypatch.setattr(sc.LibraryPanel, "_repopulate", repop)
        monkeypatch.setattr(sc.PlasmidApp, "_apply_record", apply)
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 40)) as pilot:
            ui = threading.get_ident()
            await pilot.pause()
            seen["repopulate"].clear()
            seen["apply"].clear()
            app._seed_demo_data()
            try:
                await app.workers.wait_for_complete()
            except Exception:
                pass
            for _ in range(5):
                await pilot.pause()
        assert seen["repopulate"] and seen["apply"], seen
        assert seen["repopulate"] == {ui} and seen["apply"] == {ui}, seen

    async def test_seed_populates_demo_library(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 40)) as pilot:
            await pilot.pause()
            app._seed_demo_data()
            try:
                await app.workers.wait_for_complete()
            except Exception:
                pass
            for _ in range(5):
                await pilot.pause()
        # --- plasmids → "Demo plasmids" collection (the FFE MoClo set) ---
        colls = sc._load_collections()
        names = {c.get("name") for c in colls}
        assert "Demo plasmids" in names, names
        demo = next(c for c in colls if c.get("name") == "Demo plasmids")
        plas = demo.get("plasmids") or []
        assert plas, "demo collection empty"
        assert sc._get_active_collection_name() == "Demo plasmids"
        assert {e.get("source") for e in plas} == {"demo"}
        assert len(plas) == 10, len(plas)
        pnames = {e.get("name") for e in plas}
        for want in ("FFE 6 ENTRY pCambia2300-GREEN", "FFE 1 ENTRY UPD",
                     "FFE 10 CDS cscA", "FFE 13 TU J23100-cscA-T0"):
            assert want in pnames, (want, sorted(pnames))
        # display names stay space-form, never the underscored LOCUS (INV-98)
        assert all("_" not in (e.get("name") or "") for e in plas), pnames
        # --- L0 parts → "Demo parts" bin ---
        bnames = {b.get("name") for b in sc._load_parts_bin_collections()}
        assert "Demo parts" in bnames, bnames
        parts = sc._load_parts_bin()
        assert len(parts) == 7, len(parts)
        assert {"Promoter", "CDS", "Terminator"} <= {p.get("type") for p in parts}
        # --- Golden Braid L0 entry vectors (UPD + Alpha/Omega) ---
        evs = sc._load_entry_vectors()
        assert len(evs) == 5, len(evs)
        assert {"Alpha1", "Alpha2", "Omega1", "Omega2"} <= {
            e.get("role") for e in evs}

    def test_seed_builders_are_valid_and_scrubbable(self):
        # The scrub demo plasmid must actually Golden-Braid-scrub clean, and the
        # operon's two genes must be valid ORFs — else the demo's worked
        # examples wouldn't exercise the tools.
        scrub = sc._make_demo_scrub_plasmid()
        cds = next(f for f in scrub.features if f.type == "CDS")
        s, e = int(cds.location.start), int(cds.location.end)
        plan = sc._scrub_gb_design(str(scrub.seq),
                                   [{"type": "CDS", "start": s, "end": e,
                                     "strand": 1, "label": "demoReporter"}],
                                   ["BsaI", "Esp3I"], circular=True)
        assert plan["ok"] and plan["verified"], plan.get("errors")
        assert not sc._scrub_scan_targets(plan["cured_seq"],
                                          frozenset(["BsaI", "Esp3I"]), True)
        op = sc._make_demo_operon_plasmid()
        cdss = [f for f in op.features if f.type == "CDS"]
        assert len(cdss) == 2
        for f in cdss:
            prot = sc._translate_cds(str(op.seq), int(f.location.start),
                                     int(f.location.end), 1)
            assert prot.startswith("M") and prot.count("*") == 1, prot
