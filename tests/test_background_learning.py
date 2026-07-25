"""Background Learning agent endpoints (learn-start / -status / -results / -list).

The crawl ENGINE lives in the babs repo (bb_learn.py, tested there); these cover the SpliceCraft
side: registration + write flags, the human-armed ``allow_online_lookups`` gate on ``learn-start``
(so an autonomous agent can't self-arm a crawl), slug path-safety, and the read-only session
readers. No real launch — ``learn-start`` is only exercised disarmed (403), with an invalid topic,
or with the launcher monkeypatched, so the suite never spawns a crawl or touches the network."""
import json
import pathlib
import types


import pytest

import splicecraft as sc

H = sc._state._AGENT_HANDLERS
_LEARN = ("learn-start", "learn-status", "learn-results", "learn-list")


def _call(endpoint, payload, app=None):
    return H[endpoint][0](app, payload)


def test_learn_endpoints_registered_with_write_flags():
    for e in _LEARN:
        assert e in H, f"{e} not registered in _AGENT_HANDLERS"
    assert H["learn-start"][1] is True                 # write -> the in-app Babs ask-gates it
    assert not any(H[e][1] for e in _LEARN[1:])         # status/results/list are read-only


def test_slug_is_safe():
    assert sc._learn_slug("T7 expression in E. coli!") == "learn_t7_expression_in_e_coli"
    assert sc._learn_slug("").startswith("learn_")                 # never empty -> always loadable
    s = sc._learn_slug("../../etc/passwd")
    assert s.startswith("learn_") and "/" not in s and "." not in s   # path-traversal defanged


def test_learn_start_refuses_when_disarmed(monkeypatch):
    monkeypatch.setattr(sc, "_get_setting", lambda k, d=None: False)
    body, status = _call("learn-start", {"topic": "widget synthesis"})
    assert status == 403
    assert "disarmed" in body["error"] and "allow_online_lookups" in body["error"]


def test_learn_start_validates_topic_when_armed(monkeypatch):
    monkeypatch.setattr(sc, "_get_setting", lambda k, d=None: True)
    assert _call("learn-start", {})[1] == 400                    # missing topic
    assert _call("learn-start", {"topic": "   "})[1] == 400      # blank
    assert _call("learn-start", {"topic": "x" * 201})[1] == 400  # too long


def test_learn_start_launches_and_caps_when_armed(monkeypatch):
    monkeypatch.setattr(sc, "_get_setting", lambda k, d=None: True)
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: pathlib.Path("/tmp/babs"))
    calls = {}

    def fake_launch(topic, max_docs=None, max_minutes=None):
        calls.update(topic=topic, max_docs=max_docs, max_minutes=max_minutes)
        return {"slug": sc._learn_slug(topic), "session_dir": "/x", "pid": 1234,
                "log": "/x.log", "topic": topic}

    monkeypatch.setattr(sc, "_launch_learning_session", fake_launch)
    body = _call("learn-start", {"topic": "T7 expression", "max_docs": 999, "max_minutes": 999})
    assert body["ok"] and body["slug"] == "learn_t7_expression"
    assert body["max_docs"] == 200 and body["max_minutes"] == 240   # capped at the hub ceilings
    assert calls["max_docs"] == 200 and calls["max_minutes"] == 240  # the cap reaches the launcher


def test_learn_status_slug_validation(monkeypatch):
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: pathlib.Path("/tmp/nope"))
    assert _call("learn-status", {"slug": "../etc/passwd"})[1] == 400
    assert _call("learn-status", {"slug": "not_a_learn_slug"})[1] == 400
    body, status = _call("learn-status", {"slug": "learn_absent"})
    assert status == 404 and body["status"] == "unknown"


def test_learn_status_reads_progress(tmp_path, monkeypatch):
    sess = tmp_path / "knowledge_base" / "packs" / "learn_topic" / "_session"
    sess.mkdir(parents=True)
    (sess / "progress.json").write_text(json.dumps({"kept": 3, "status": "running"}))
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
    body = _call("learn-status", {"slug": "learn_topic"})
    assert body["slug"] == "learn_topic" and body["kept"] == 3 and body["status"] == "running"


def test_learn_results_dedupes_and_sorts_by_score(tmp_path, monkeypatch):
    pack = tmp_path / "knowledge_base" / "packs" / "learn_topic"
    pack.mkdir(parents=True)
    rows = [{"doc_id": "A", "title": "a", "relevance_score": 0.4},
            {"doc_id": "A", "title": "a", "relevance_score": 0.4},   # dup doc -> collapsed
            {"doc_id": "B", "title": "b", "relevance_score": 0.9}]
    (pack / "corpus.jsonl").write_text("\n".join(json.dumps(r) for r in rows))
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
    body = _call("learn-results", {"slug": "learn_topic"})
    assert body["count"] == 2 and [d["doc_id"] for d in body["results"]] == ["B", "A"]


def test_learn_list_scans_only_learn_packs(tmp_path, monkeypatch):
    for slug, kept in (("learn_alpha", 2), ("learn_beta", 5)):
        sess = tmp_path / "knowledge_base" / "packs" / slug / "_session"
        sess.mkdir(parents=True)
        (sess / "anchor.json").write_text(json.dumps({"topic": slug}))
        (sess / "progress.json").write_text(json.dumps({"kept": kept, "status": "done"}))
    (tmp_path / "knowledge_base" / "packs" / "plant_tc").mkdir(parents=True)   # not a learn pack
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
    body = _call("learn-list", {})
    assert {s["slug"] for s in body["sessions"]} == {"learn_alpha", "learn_beta"}


def test_learn_endpoints_error_without_babs(monkeypatch):
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: None)
    assert _call("learn-list", {})[1] == 400
    assert _call("learn-results", {"slug": "learn_x"})[1] == 400
    assert _call("learn-status", {"slug": "learn_x"})[1] == 400


def test_learn_status_reconciles_dead_pid(tmp_path, monkeypatch):
    # a killed / crashed / OOM'd session leaves progress.json frozen at "running"; reconcile against
    # the recorded pid so the reader shows "stopped" instead of a forever-"running" lie
    sess = tmp_path / "knowledge_base" / "packs" / "learn_topic" / "_session"
    sess.mkdir(parents=True)
    (sess / "progress.json").write_text(json.dumps({"status": "running", "kept": 2, "pid": 2147480000}))
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
    assert _call("learn-status", {"slug": "learn_topic"})["status"] == "stopped"


def test_learn_start_refuses_duplicate_running_session(tmp_path, monkeypatch):
    import os
    monkeypatch.setattr(sc, "_get_setting", lambda k, d=None: True)
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
    slug = sc._learn_slug("t7 expression")
    sess = tmp_path / "knowledge_base" / "packs" / slug / "_session"
    sess.mkdir(parents=True)
    (sess / "progress.json").write_text(json.dumps({"status": "running", "pid": os.getpid()}))  # alive
    body, status = _call("learn-start", {"topic": "t7 expression"})
    assert status == 409 and "already running" in body["error"]


def test_learn_status_survives_nondict_progress(tmp_path, monkeypatch):
    # a corrupt/truncated progress.json that parses to a JSON list (not an object) must not crash the
    # reader — it used to reach ``{**prog}`` -> TypeError; now it reads as absent -> 404 unknown
    sess = tmp_path / "knowledge_base" / "packs" / "learn_topic" / "_session"
    sess.mkdir(parents=True)
    (sess / "progress.json").write_text("[1, 2, 3]")
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
    body, status = _call("learn-status", {"slug": "learn_topic"})
    assert status == 404 and body["status"] == "unknown"


def test_learn_results_skips_nondict_and_bad_lines(tmp_path, monkeypatch):
    # a corpus.jsonl line that is valid JSON but not an object (list / scalar) used to AttributeError
    # on rec.get; invalid-JSON lines are skipped too. Only the real records survive.
    pack = tmp_path / "knowledge_base" / "packs" / "learn_topic"
    pack.mkdir(parents=True)
    (pack / "corpus.jsonl").write_text("\n".join([
        '{"doc_id": "A", "title": "a", "relevance_score": 0.5}',
        '["not", "an", "object"]',   # valid JSON, not a dict -> skipped (no crash)
        '42',                         # valid JSON scalar -> skipped
        'not json at all',            # invalid JSON -> skipped
        '{"doc_id": "B", "title": "b", "relevance_score": 0.9}']))
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
    body = _call("learn-results", {"slug": "learn_topic"})
    assert body["count"] == 2 and [d["doc_id"] for d in body["results"]] == ["B", "A"]


def test_learn_list_survives_nondict_session_files(tmp_path, monkeypatch):
    # a truthy non-dict anchor/progress.json ([...]) used to AttributeError on (x or {}).get; now
    # such a session reads as absent and is skipped, not crashed
    sess = tmp_path / "knowledge_base" / "packs" / "learn_bad" / "_session"
    sess.mkdir(parents=True)
    (sess / "anchor.json").write_text('["junk"]')
    (sess / "progress.json").write_text('["junk"]')
    monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
    assert _call("learn-list", {})["sessions"] == []


def test_learn_pid_alive_assumes_alive_off_posix(monkeypatch):
    # on Windows os.kill(pid, 0) routes through TerminateProcess and would KILL the pid, so off-POSIX
    # we must assume-alive WITHOUT ever probing
    monkeypatch.setattr(sc.os, "name", "nt")

    def _boom(*a, **k):
        raise AssertionError("os.kill must not be called off-POSIX")

    monkeypatch.setattr(sc.os, "kill", _boom)
    assert sc._learn_pid_alive(999999) is True


class TestLearningLoopCloses:
    """The learning loop must CLOSE: a session that collects papers has to end up as knowledge
    Babs can actually recall. Pre-fix it did not — `bb_learn` writes only the pack's
    corpus.jsonl, nothing ever ran `bb_index` on it, and no chat path consulted a learn pack.
    Verified on a real box: a completed session's pack existed with kept papers while Chroma
    held no `learn_*` collection at all, so everything it learned was unanswerable.
    """

    def test_launch_chains_the_index_step_scoped_to_the_pack(self, monkeypatch, tmp_path):
        """The launch chain must run bb_index AFTER bb_learn, scoped to the session's own pack
        via BABS_PACK — otherwise the crawl's output is never embedded."""
        home = tmp_path / "babs"
        (home / "knowledge_base" / "jobs").mkdir(parents=True)
        (home / "bb_learn.py").write_text("")
        monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: home)
        monkeypatch.setattr(sc, "_learn_babs_python", lambda h: "/usr/bin/python3")
        seen = {}

        class _P:
            pid = 4242

        def fake_popen(argv, **kwargs):
            seen["argv"] = argv
            return _P()
        import subprocess
        monkeypatch.setattr(subprocess, "Popen", fake_popen)
        info = sc._launch_learning_session("T7 expression", max_docs=10)
        chain = seen["argv"][2]                      # ["bash", "-c", "<chain>"]
        assert "bb_learn.py" in chain
        assert "bb_index.py" in chain, "learned corpus is never indexed — the loop stays open"
        assert chain.index("bb_learn.py") < chain.index("bb_index.py")
        assert f"BABS_PACK={info['slug']}" in chain, "index must be scoped to the session pack"
        # `;` not `&&`: a budget-exhausted or interrupted crawl must still index what it kept.
        assert " ; " in chain and "&&" not in chain

    # — merge into the main corpus —
    def _mk(self, tmp_path, main_lines, pack_lines, slug="learn_x"):
        home = tmp_path / "babs"
        kb = home / "knowledge_base"
        (kb / "packs" / slug).mkdir(parents=True)
        (kb / "corpus.jsonl").write_text("".join(json.dumps(r) + "\n" for r in main_lines))
        (kb / "packs" / slug / "corpus.jsonl").write_text(
            "".join(json.dumps(r) + "\n" for r in pack_lines))
        return home

    def test_merge_appends_new_chunks_and_snapshots_first(self, tmp_path):
        home = self._mk(tmp_path, [{"chunk_id": "a", "text": "old"}],
                        [{"chunk_id": "b", "text": "new"}])
        stats = sc._learn_merge_pack_into_main(home, "learn_x")
        assert stats["added"] == 1 and stats["skipped"] == 0
        body = (home / "knowledge_base" / "corpus.jsonl").read_text()
        assert '"chunk_id": "a"' in body and '"chunk_id": "b"' in body
        # The pre-merge snapshot exists: this appends to hours of crawling.
        assert (home / "knowledge_base" / "corpus.jsonl.premerge-bak").is_file()

    def test_merge_is_idempotent(self, tmp_path):
        """Re-merging must be a no-op, not a duplicate-inflating mistake."""
        home = self._mk(tmp_path, [{"chunk_id": "a", "text": "old"}],
                        [{"chunk_id": "a", "text": "old"}, {"chunk_id": "b", "text": "new"}])
        first = sc._learn_merge_pack_into_main(home, "learn_x")
        assert first == {"added": 1, "skipped": 1, "total": 2}
        second = sc._learn_merge_pack_into_main(home, "learn_x")
        assert second["added"] == 0 and second["skipped"] == 2
        lines = [ln for ln in (home / "knowledge_base" / "corpus.jsonl")
                 .read_text().splitlines() if ln.strip()]
        assert len(lines) == 2

    def test_merge_skips_corrupt_lines(self, tmp_path):
        home = self._mk(tmp_path, [{"chunk_id": "a"}], [{"chunk_id": "b"}])
        pack = home / "knowledge_base" / "packs" / "learn_x" / "corpus.jsonl"
        pack.write_text('not json\n[]\n' + pack.read_text())
        stats = sc._learn_merge_pack_into_main(home, "learn_x")
        assert stats["added"] == 1          # the object line only

    def test_merge_without_a_pack_corpus_is_an_actionable_error(self, tmp_path):
        home = tmp_path / "babs"
        (home / "knowledge_base" / "packs" / "learn_x").mkdir(parents=True)
        with pytest.raises(RuntimeError, match="no corpus"):
            sc._learn_merge_pack_into_main(home, "learn_x")

    # — endpoint —
    def test_merge_endpoint_is_a_write(self):
        assert sc._state._AGENT_HANDLERS["learn-merge"][1] is True

    def test_merge_endpoint_validates_the_slug(self):
        fn = sc._state._AGENT_HANDLERS["learn-merge"][0]
        body, status = fn(None, {"slug": "../../etc"})
        assert status == 400

    def test_merge_endpoint_refuses_a_running_session(self, monkeypatch, tmp_path):
        monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
        monkeypatch.setattr(sc, "_learn_session_running", lambda h, s: True)
        fn = sc._state._AGENT_HANDLERS["learn-merge"][0]
        body, status = fn(None, {"slug": "learn_x"})
        assert status == 409 and "still running" in body["error"]


class TestIngestUrl:
    """`/ingest <url>` — the bridge from "Babs just read this" to "Babs knows this".

    `read-url` puts a page in ONE turn's context and drops it, so the same lookup is paid for
    forever and nothing accumulates. Ingest persists it. The open-licence refusal is inherited
    from the babs engine on purpose: the corpus stays redistributable.
    """

    def _fake(self, monkeypatch, tmp_path, stdout, rc=0):
        import subprocess
        seen = {}

        def run(argv, **kwargs):
            seen["argv"] = argv
            return types.SimpleNamespace(returncode=rc, stdout=stdout, stderr="")
        monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: tmp_path)
        monkeypatch.setattr(sc, "_learn_babs_python", lambda h: "/usr/bin/python3")
        monkeypatch.setattr(sc, "_learn_reindex_main", lambda h: 999)
        monkeypatch.setattr(subprocess, "run", run)
        return seen

    def test_rejects_a_non_http_scheme_without_shelling_out(self, monkeypatch, tmp_path):
        seen = self._fake(monkeypatch, tmp_path, "{}")
        for bad in ("file:///etc/passwd", "ftp://x/y", "javascript:alert(1)", ""):
            res = sc._babs_ingest_url(bad)
            assert res["ok"] is False
        assert "argv" not in seen

    def test_success_indexes_afterwards(self, monkeypatch, tmp_path):
        """Without the chained re-index the ingested text is never embedded — the same open
        loop Background Learning had."""
        seen = self._fake(monkeypatch, tmp_path, json.dumps(
            {"ok": True, "doc_id": "web::abc", "title": "T", "chunks": 4, "license": "CC-BY"}))
        res = sc._babs_ingest_url("https://example.org/paper")
        assert res["ok"] and res["chunks"] == 4 and res["indexed"] is True
        assert "--ingest-url" in seen["argv"]

    def test_parses_json_past_a_diagnostic_line(self, monkeypatch, tmp_path):
        self._fake(monkeypatch, tmp_path,
                   '[warn] noisy\n{"ok": true, "doc_id": "d", "chunks": 1}\n')
        assert sc._babs_ingest_url("https://example.org/x")["ok"] is True

    def test_unparseable_output_is_an_honest_failure(self, monkeypatch, tmp_path):
        self._fake(monkeypatch, tmp_path, "totally not json", rc=1)
        res = sc._babs_ingest_url("https://example.org/x")
        assert res["ok"] is False and "up to date" in res["reason"]

    def test_no_babs_repo_is_actionable(self, monkeypatch):
        monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: None)
        assert "babs-setup" in sc._babs_ingest_url("https://x.org").get("reason", "") or \
            "not found" in sc._babs_ingest_url("https://x.org")["reason"]

    # — endpoint —
    def test_endpoint_is_gated_on_the_human_armed_setting(self, monkeypatch):
        """Ingest reaches the network; an autonomous agent cannot self-arm it."""
        monkeypatch.setattr(sc, "_get_setting",
                            lambda k, d=None: False if k == "allow_online_lookups" else d)
        fn = sc._state._AGENT_HANDLERS["ingest-url"][0]
        body, status = fn(None, {"url": "https://example.org"})
        assert status == 403

    def test_endpoint_is_a_write(self):
        assert sc._state._AGENT_HANDLERS["ingest-url"][1] is True

    def test_endpoint_reports_a_duplicate_as_a_benign_noop(self, monkeypatch):
        monkeypatch.setattr(sc, "_get_setting",
                            lambda k, d=None: True if k == "allow_online_lookups" else d)
        monkeypatch.setattr(sc, "_babs_ingest_url", lambda u, t="": {
            "ok": False, "duplicate": True, "doc_id": "web::abc"})
        fn = sc._state._AGENT_HANDLERS["ingest-url"][0]
        out = fn(None, {"url": "https://example.org"})
        assert not isinstance(out, tuple)          # 200, not an error status
        assert out["duplicate"] is True

    def test_endpoint_maps_a_transient_failure_to_502(self, monkeypatch):
        monkeypatch.setattr(sc, "_get_setting",
                            lambda k, d=None: True if k == "allow_online_lookups" else d)
        monkeypatch.setattr(sc, "_babs_ingest_url", lambda u, t="": {
            "ok": False, "transient": True, "reason": "temporary fetch failure"})
        fn = sc._state._AGENT_HANDLERS["ingest-url"][0]
        _body, status = fn(None, {"url": "https://example.org"})
        assert status == 502


class TestPrivateLibraryPack:
    """Indexing the user's OWN work as recallable knowledge.

    The richest project-specific knowledge on the box is the user's library, and Babs could
    only reach it one endpoint at a time inside a 7B tool loop — "which of my plasmids has a
    KanR and a T7 promoter?" was effectively unanswerable. It is indexed as a Babs pack, but
    it is UNPUBLISHED USER WORK, so the pack is barred from every egress path we control
    rather than merely defaulting to off.
    """

    @pytest.fixture
    def babs(self, tmp_path, monkeypatch):
        home = tmp_path / "babs"
        (home / "packs").mkdir(parents=True)
        (home / "knowledge_base").mkdir(parents=True)
        monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: home)
        return home

    def _seed(self, monkeypatch, gb='     CDS  1..9\n                     /label="KanR"\n'):
        monkeypatch.setattr(sc, "_iter_library_readonly", lambda: [
            {"name": "pTEST", "id": "e1", "size": 2686, "n_feats": 1,
             "source": "manual", "gb_text": "LOCUS x 2686 bp DNA circular\n" + gb}])
        monkeypatch.setattr(sc, "_load_primers", lambda: [
            {"name": "T7-fwd", "sequence": "TAATACGACTCACTATAGGG"}])
        monkeypatch.setattr(sc, "_load_parts_bin", lambda: [])
        monkeypatch.setattr(sc, "_load_experiments", lambda: [])

    def test_every_record_is_flagged_private(self, babs, monkeypatch):
        self._seed(monkeypatch)
        sc._write_library_pack(babs)
        corpus = babs / "knowledge_base" / "packs" / "splicecraft_library" / "corpus.jsonl"
        recs = [json.loads(l) for l in corpus.read_text().splitlines() if l.strip()]
        assert recs and all(r["private"] is True for r in recs)

    def test_pack_dir_carries_the_no_egress_marker(self, babs, monkeypatch):
        """push_corpus.sh keys off this marker to skip the whole directory."""
        self._seed(monkeypatch)
        sc._write_library_pack(babs)
        marker = (babs / "knowledge_base" / "packs" / "splicecraft_library" / ".no-egress")
        assert marker.is_file() and "Do not sync" in marker.read_text()

    def test_feature_labels_are_indexed(self, babs, monkeypatch):
        """The whole point: a feature-level index is what makes "which of my plasmids has a
        KanR?" answerable — a name-and-size index would not be."""
        self._seed(monkeypatch)
        sc._write_library_pack(babs)
        corpus = babs / "knowledge_base" / "packs" / "splicecraft_library" / "corpus.jsonl"
        text = corpus.read_text()
        assert "KanR" in text and "T7-fwd" in text

    def test_generated_pack_module_declares_private(self, babs, monkeypatch):
        self._seed(monkeypatch)
        sc._write_library_pack(babs)
        mod = (babs / "packs" / "splicecraft_library.py").read_text()
        assert "PRIVATE = True" in mod and "no-egress" in mod
        # It must satisfy the babs pack contract or bb_config refuses to load it.
        for field in ("NAME", "DISPLAY_NAME", "COLLECTION", "TOPIC_TERMS", "METHOD_TERMS",
                      "CONTEXT_MARKERS", "TOPIC_MARKERS", "METHOD_MARKERS",
                      "CONSTRUCTION_MARKERS", "MIN_MARKERS", "EXCLUDE_KEYWORDS",
                      "POLICY_TITLE_TERMS", "METHOD_NUMBER_SIGNAL"):
            assert f"{field} =" in mod, f"pack module missing required field {field}"

    def test_rebuild_replaces_rather_than_appends(self, babs, monkeypatch):
        """The pack is rebuilt wholesale, so a deleted plasmid can't linger as a ghost Babs
        still tells the user about."""
        self._seed(monkeypatch)
        sc._write_library_pack(babs)
        first = (babs / "knowledge_base" / "packs" / "splicecraft_library"
                 / "corpus.jsonl").read_text()
        sc._write_library_pack(babs)
        second = (babs / "knowledge_base" / "packs" / "splicecraft_library"
                  / "corpus.jsonl").read_text()
        assert second.count("pTEST") == first.count("pTEST")

    def test_empty_library_is_an_actionable_refusal(self, babs, monkeypatch):
        for fn in ("_load_primers", "_load_parts_bin", "_load_experiments"):
            monkeypatch.setattr(sc, fn, lambda: [])
        monkeypatch.setattr(sc, "_iter_library_readonly", lambda: [])
        with pytest.raises(RuntimeError, match="nothing to index"):
            sc._write_library_pack(babs)

    def test_one_failing_loader_does_not_sink_the_index(self, babs, monkeypatch):
        """A corrupt primer store must not cost the user their plasmid index."""
        self._seed(monkeypatch)
        monkeypatch.setattr(sc, "_load_primers",
                            lambda: (_ for _ in ()).throw(RuntimeError("boom")))
        stats = sc._write_library_pack(babs)
        assert stats["docs"] >= 1

    def test_pack_is_searchable_and_ranks_ahead_of_learn_packs(self, babs, monkeypatch):
        """A question about "my plasmids" must not be crowded out of the recall fan-out by
        learning packs when several exist."""
        self._seed(monkeypatch)
        sc._write_library_pack(babs)
        for slug in ("learn_a", "learn_b"):
            (babs / "packs" / f"{slug}.py").write_text("")
            d = babs / "knowledge_base" / "packs" / slug
            d.mkdir(parents=True)
            (d / "corpus.jsonl").write_text('{"chunk_id": "x"}\n')
        packs = sc._recall_searchable_packs(babs)
        assert packs[0] == ""                       # the default corpus leads
        assert packs[1] == "splicecraft_library"    # then the user's own work

    # — endpoint —
    def test_endpoint_is_a_write(self):
        assert sc._state._AGENT_HANDLERS["index-library"][1] is True

    def test_endpoint_reports_privacy(self, babs, monkeypatch):
        self._seed(monkeypatch)
        monkeypatch.setattr(sc, "_index_library_pack", lambda h: 4242)
        out = sc._state._AGENT_HANDLERS["index-library"][0](None, {})
        assert out["private"] is True and out["pack"] == "splicecraft_library"
        assert out["documents"] >= 1

    def test_endpoint_without_babs_is_actionable(self, monkeypatch):
        monkeypatch.setattr(sc, "_learn_resolve_babs_home", lambda: None)
        body, status = sc._state._AGENT_HANDLERS["index-library"][0](None, {})
        assert status == 400 and "babs" in body["error"].lower()
