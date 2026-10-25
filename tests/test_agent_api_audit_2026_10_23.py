"""
test_agent_api_audit_2026_10_23 — agent-API / CLI findings AG3–AG12 from
the 2026-10-23 audit (re-verified on the uncommitted tree 2026-10-24).

Each test FAILED before the fix:
  * AG3  a snapshot-rebuild handler (replace-sequence, apply-gff3, …)
         installed a record that silently dropped an add-feature it had
         already acknowledged — the in-place edit never moved the load
         counter it checked;
  * AG4  set-origin / reverse-complement / update-primer acted on a record
         loaded meanwhile (update-primer could rename ANOTHER plasmid's
         primer, and reported refused edits as applied);
  * AG5  a `--read-only` guest answered from its first read forever;
  * AG7  two requests with one idempotency key both ran; an invalid key was
         silently ignored;
  * AG8–AG12 CLI tracebacks, torn reads, NaN/Infinity/huge-int 500s,
         `_stale` on editable installs, a double-counted primer collection,
         unbounded connection threads, BrokenPipe tracebacks.
"""
from __future__ import annotations

import http.client
import json
import os
import socket
import threading
import time
from unittest.mock import patch

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import splicecraft as sc
import splicecraft_cli as cli


# ── shared helpers (real PlasmidApp, handler on its own thread) ────────────

def _rec(name, n, seed, feats=()):
    import random
    r = random.Random(seed)
    rec = SeqRecord(Seq("".join(r.choice("ACGT") for _ in range(n))),
                    id=name, name=name,
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    for s, e, t, lab in feats:
        q = {"label": [lab]}
        if t == "primer_bind":
            q["primer_seq"] = [str(rec.seq[s:e])]
        rec.features.append(SeqFeature(FeatureLocation(s, e, strand=1),
                                       type=t, qualifiers=q))
    return rec


def _labels(app):
    return [(f.qualifiers.get("label") or [""])[0]
            for f in app._current_record.features]


async def _run_off_ui(pilot, fn):
    """Run an agent handler on a non-UI thread, as the HTTP server does."""
    out, err = [], []

    def _t():
        try:
            out.append(fn())
        except BaseException as exc:          # surface, never hide
            err.append(exc)
    t = threading.Thread(target=_t, daemon=True)
    t.start()
    for _ in range(400):
        await pilot.pause()
        if not t.is_alive():
            break
    t.join(timeout=5)
    assert not err, err
    assert out, "handler thread did not finish"
    return out[0]


async def _mounted(app, pilot):
    for _ in range(4):
        await pilot.pause()
    while len(app.screen_stack) > 1:
        app.pop_screen()
        await pilot.pause()


# ── AG3: a snapshot-rebuild handler must not drop an acknowledged edit ─────

class TestSnapshotEditsDoNotDropConcurrentEdits:
    async def test_replace_sequence_refuses_rather_than_dropping_an_add_feature(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await _mounted(app, pilot)
            app._apply_record(_rec("RecA", 3000, 1, [(50, 110, "misc_feature", "f1")]))
            await pilot.pause()
            added = {}
            real_rebuild = app._rebuild_record_with_edit

            def rebuild_then_concurrent_add(*a, **k):
                out = real_rebuild(*a, **k)          # the real rebuild
                # a second request lands while the rebuild runs off-UI
                added["r"] = sc._h_add_feature(
                    app, {"start": 1000, "end": 1100, "label": "LATE"})
                return out
            app._rebuild_record_with_edit = rebuild_then_concurrent_add
            r = await _run_off_ui(pilot, lambda: sc._h_replace_sequence(
                app, {"start": 10, "end": 20, "bases": "ACGTACGTAC",
                      "force": True}))
            del app._rebuild_record_with_edit
            assert added["r"].get("ok") is True       # the add was acknowledged…
            assert "LATE" in _labels(app)             # …so it must still be there
            assert isinstance(r, tuple) and r[1] == 409, r

    async def test_apply_gff3_refuses_rather_than_dropping_an_add_feature(
            self, tmp_path, monkeypatch):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await _mounted(app, pilot)
            app._apply_record(_rec("RecB", 3000, 2, [(50, 110, "misc_feature", "f1")]))
            await pilot.pause()
            gff = tmp_path / "x.gff3"
            gff.write_text("##gff-version 3\n##sequence-region RecB 1 3000\n"
                           "RecB\tx\tmisc_feature\t200\t260\t.\t+\t.\tID=g1\n")
            real = sc._gff3_apply_to_loaded_record

            def parse_then_concurrent_add(rec, path):
                n = real(rec, path)
                sc._h_add_feature(app, {"start": 1000, "end": 1100,
                                        "label": "LATE"})
                return n
            monkeypatch.setattr(sc, "_gff3_apply_to_loaded_record",
                                parse_then_concurrent_add)
            r = await _run_off_ui(pilot, lambda: sc._h_apply_gff3(
                app, {"path": str(gff), "force": True}))
            assert "LATE" in _labels(app)
            assert isinstance(r, tuple) and r[1] == 409, r


# ── AG4: whole-record reframes + update-primer re-check the record ─────────

class TestReframesAndPrimerEditsAreStaleGuarded:
    async def test_set_origin_does_not_overwrite_a_record_loaded_meanwhile(
            self, monkeypatch):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await _mounted(app, pilot)
            old = _rec("Old", 3000, 4)
            new = _rec("New", 2000, 5)
            app._apply_record(old)
            await pilot.pause()
            real = sc._agent_resolve_origin_offset

            def resolve_then_load(rec, payload, n):
                out = real(rec, payload, n)
                app.call_from_thread(lambda: app._apply_record(new))
                return out
            monkeypatch.setattr(sc, "_agent_resolve_origin_offset",
                                resolve_then_load)
            r = await _run_off_ui(pilot, lambda: sc._h_set_origin(app, {"bp": 1000}))
            assert str(app._current_record.seq) == str(new.seq)
            assert isinstance(r, tuple) and r[1] == 409, r

    async def test_reverse_complement_does_not_overwrite_a_record_loaded_meanwhile(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await _mounted(app, pilot)
            old = _rec("OldRC", 3000, 7)
            new = _rec("NewRC", 2000, 8)
            app._apply_record(old)
            await pilot.pause()
            real_cft, first = app.call_from_thread, [True]

            def load_first(fn, *a, **k):              # the concurrent load
                if first[0]:
                    first[0] = False
                    real_cft(lambda: app._apply_record(new))
                return real_cft(fn, *a, **k)
            app.call_from_thread = load_first
            r = await _run_off_ui(pilot, lambda: sc._h_reverse_complement(app, {}))
            del app.call_from_thread
            assert str(app._current_record.seq) == str(new.seq)
            assert isinstance(r, tuple) and r[1] == 409, r

    async def test_update_primer_never_edits_another_records_primer(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await _mounted(app, pilot)
            a = _rec("PA", 3000, 11, [(100, 120, "primer_bind", "a_fwd")])
            b = _rec("PB", 3000, 12, [(300, 320, "primer_bind", "b_fwd")])
            app._apply_record(a)
            await pilot.pause()
            real_cft, first = app.call_from_thread, [True]

            def load_first(fn, *a_, **k):
                if first[0]:
                    first[0] = False
                    real_cft(lambda: app._apply_record(b))
                return real_cft(fn, *a_, **k)
            app.call_from_thread = load_first
            r = await _run_off_ui(pilot, lambda: sc._h_update_primer(
                app, {"idx": 0, "label": "RENAMED", "force": True}))
            del app.call_from_thread
            assert _labels(app) == ["b_fwd"]            # B untouched
            assert isinstance(r, tuple) and r[1] == 409, r

    async def test_update_primer_never_reports_a_refused_edit_as_applied(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await _mounted(app, pilot)
            a = _rec("PC", 3000, 13, [(100, 120, "primer_bind", "a_fwd")])
            b = _rec("PD", 3000, 14, [(100, 200, "CDS", "not_a_primer")])
            app._apply_record(a)
            await pilot.pause()
            real_cft, first = app.call_from_thread, [True]

            def load_first(fn, *a_, **k):
                if first[0]:
                    first[0] = False
                    real_cft(lambda: app._apply_record(b))
                return real_cft(fn, *a_, **k)
            app.call_from_thread = load_first
            r = await _run_off_ui(pilot, lambda: sc._h_update_primer(
                app, {"idx": 0, "label": "RENAMED", "force": True}))
            del app.call_from_thread
            assert not (isinstance(r, dict) and r.get("ok")), r


# ── AG5: a --read-only guest reads what the host wrote ─────────────────────

class TestReadOnlyGuestSeesHostWrites:
    def test_list_collections_reflects_a_host_write(self, monkeypatch):
        sc._save_collections([{"name": "Main", "plasmids": []}])
        monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY", True)
        p1, s1 = sc._agent_invoke(None, "list-collections", {})
        assert s1 == 200
        assert [c["name"] for c in p1["collections"]] == ["Main"]
        # The HOST process rewrites the file (atomic replace) — simulated
        # without touching this process's caches, as another process would.
        path = sc._state._COLLECTIONS_FILE
        tmp = path.with_suffix(".host-tmp")
        tmp.write_text(json.dumps({"_schema_version": 1, "entries": [
            {"name": "Main", "plasmids": []},
            {"name": "HostMade", "plasmids": []}]}), encoding="utf-8")
        os.replace(tmp, path)
        p2, s2 = sc._agent_invoke(None, "list-collections", {})
        assert s2 == 200
        assert "HostMade" in [c["name"] for c in p2["collections"]]


# ── AG7: idempotency keys protect overlapping retries; bad keys refused ────

def _free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


def _post(port, token, path, body, headers=None):
    import http.client
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=20)
    h = {"Authorization": f"Bearer {token}", "Content-Type": "application/json"}
    h.update(headers or {})
    c.request("POST", path, body=json.dumps(body).encode(), headers=h)
    r = c.getresponse()
    data = r.read()
    c.close()
    return r.status, json.loads(data or b"{}")


@pytest.fixture
def slow_write_server(tmp_path, monkeypatch):
    """A real agent server with one extra WRITE endpoint that counts its
    executions and parks the FIRST one until released."""
    entered, release = threading.Event(), threading.Event()
    calls = []

    def _h_zz_slow_create(app, payload):
        calls.append(payload.get("n"))
        if len(calls) == 1:
            entered.set()
            release.wait(10)
        return {"ok": True, "created": len(calls)}
    monkeypatch.setitem(sc._state._AGENT_HANDLERS, "zz-slow-create",
                        (_h_zz_slow_create, True))
    monkeypatch.setattr(sc._state, "_AGENT_TOKEN_FILE", tmp_path / "agent_token")
    sc._agent_idempotency_reset()
    srv = sc._start_agent_api(None, port=_free_port())
    assert srv is not None
    port, token = (tmp_path / "agent_token").read_text().split()
    yield int(port), token, calls, entered, release
    release.set()
    sc._stop_agent_api(srv)
    sc._agent_idempotency_reset()


class TestIdempotencyUnderOverlap:
    def test_a_retry_while_the_first_still_runs_does_not_run_twice(
            self, slow_write_server):
        port, token, calls, entered, release = slow_write_server
        key = {"X-Idempotency-Key": "overlap-1"}
        first = {}
        t = threading.Thread(target=lambda: first.setdefault(
            "r", _post(port, token, "/zz-slow-create", {"n": 1}, key)))
        t.start()
        assert entered.wait(5)
        s2, b2 = _post(port, token, "/zz-slow-create", {"n": 1}, key)
        release.set()
        t.join(10)
        assert calls == [1], calls                    # executed ONCE
        assert s2 == 409 and b2.get("in_progress") is True, (s2, b2)
        s3, b3 = _post(port, token, "/zz-slow-create", {"n": 1}, key)
        assert s3 == 200 and b3.get("_idempotent_replay") is True
        assert calls == [1]

    @pytest.mark.parametrize("bad", ["bad/key", "has space", "x" * 129])
    def test_an_invalid_key_is_refused_not_ignored(self, slow_write_server, bad):
        port, token, calls, _entered, release = slow_write_server
        release.set()
        s, b = _post(port, token, "/zz-slow-create", {"n": 2},
                     {"X-Idempotency-Key": bad})
        assert s == 400, (s, b)
        assert calls == []


@pytest.fixture
def server(tmp_path, monkeypatch):
    calls = []

    def _h_zz_slow_read(app, payload):
        calls.append(1)
        time.sleep(0.3)
        return {"ok": True}
    monkeypatch.setitem(sc._state._AGENT_HANDLERS, "zz-slow-read",
                        (_h_zz_slow_read, False))
    monkeypatch.setattr(sc._state, "_AGENT_TOKEN_FILE", tmp_path / "agent_token")
    srv = sc._start_agent_api(None, port=_free_port())
    assert srv is not None
    port, token = (tmp_path / "agent_token").read_text().split()
    yield int(port), token
    sc._stop_agent_api(srv)


def _raw_post(port, token, path, text):
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=20)
    c.request("POST", path, body=text.encode(), headers={
        "Authorization": f"Bearer {token}", "Content-Type": "application/json"})
    r = c.getresponse()
    data = r.read()
    c.close()
    return r.status, json.loads(data or b"{}")


# ── AG10 + NEW: malformed numbers are a 400, never a 500 / dropped socket ──

class TestNonFiniteAndHugeNumbers:
    @pytest.mark.parametrize("body", [
        '{"protein": "MKVLA", "stops": 1e999}',
        '{"protein": "MKVLA", "transl_table": 1e999}',
        '{"protein": "MKVLA", "stops": NaN}',
    ])
    def test_optimize_protein_bad_numbers_are_400(self, server, body):
        port, token = server
        s, b = _raw_post(port, token, "/optimize-protein", body)
        assert s == 400, (s, b)

    def test_an_int_past_the_4300_digit_limit_gets_an_answer(self, server):
        port, token = server
        s, b = _raw_post(port, token, "/get-sequence",
                         '{"start": ' + "1" * 4301 + ', "end": 5}')
        assert s == 400, (s, b)

    @pytest.mark.parametrize("payload, endpoint", [
        ({"name": "F", "sequence": "ATGAAA", "strand": float("nan")}, "create-feature-library"),
        ({"name": "P", "sequence": "ATGAAA", "level": float("inf")}, "create-part"),
        ({"collection": "x", "format": [1], "folder": "/tmp/x"}, "bulk-export-collection"),
        ({"seq": "ATGAAACCCGGGTTTAAA" * 4, "features": [1, "a", None]}, "scrub-plasmid"),
    ])
    def test_in_process_payloads_never_500(self, payload, endpoint):
        # Babs passes dicts straight to `_agent_invoke` -- no JSON parse.
        out, status = sc._agent_invoke(None, endpoint, payload)
        assert status < 500, (endpoint, status, out)


# ── AG9: features / get-feature are one snapshot of one record ─────────────

def _wrap_rec(name, n, seed, note):
    import random
    r = random.Random(seed)
    rec = SeqRecord(Seq("".join(r.choice("ACGT") for _ in range(n))), id=name,
                    name=name, annotations={"molecule_type": "DNA",
                                            "topology": "circular"})
    rec.features.append(SeqFeature(CompoundLocation([
        FeatureLocation(n - 40, n, 1), FeatureLocation(0, 40, 1)]),
        type="misc_feature", qualifiers={"label": [name + "wrap"],
                                         "note": [note]}))
    return rec


class TestFeatureReadsAreCoherent:
    async def _torn_call(self, handler, payload):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            for _ in range(4):
                await pilot.pause()
            while len(app.screen_stack) > 1:
                app.pop_screen()
                await pilot.pause()
            big = _wrap_rec("BIG", 100_000, 1, "BIG note")
            small = _wrap_rec("SMALL", 1000, 2, "SMALL note")
            app._apply_record(big)
            await pilot.pause()
            real_q, fired = app.query_one, [False]

            def q(*a, **k):
                # A swap lands between the handler's reads -- possible only
                # when it reads off the UI thread.
                if (threading.current_thread() is not threading.main_thread()
                        and not fired[0]):
                    fired[0] = True
                    app.call_from_thread(lambda: app._apply_record(small))
                return real_q(*a, **k)
            app.query_one = q
            out = []
            t = threading.Thread(target=lambda: out.append(handler(app, payload)))
            t.start()
            for _ in range(400):
                await pilot.pause()
                if not t.is_alive():
                    break
            t.join(5)
            del app.query_one
            return out[0]

    async def test_features_length_matches_its_own_record(self):
        r = await self._torn_call(sc._h_features, {})
        assert [f["length"] for f in r["features"]] == [80], r

    async def test_get_feature_row_and_qualifiers_from_one_record(self):
        r = await self._torn_call(sc._h_get_feature, {"idx": 0})
        assert (r["label"], r["qualifiers"]["note"]) in (
            ("BIGwrap", ["BIG note"]), ("SMALLwrap", ["SMALL note"])), r


# ── AG11: an editable install is never reported stale ──────────────────────

class TestEditableInstallIsNotStale:
    def test_worker_reports_unknown_for_an_editable_install(self, tmp_path, monkeypatch):
        di = tmp_path / "splicecraft-0.0.1.dist-info"
        di.mkdir()
        (di / "METADATA").write_text("Metadata-Version: 2.1\nName: splicecraft\n"
                                     "Version: 0.0.1\n")
        (di / "direct_url.json").write_text(json.dumps(
            {"dir_info": {"editable": True}, "url": "file:///src/splicecraft"}))
        monkeypatch.syspath_prepend(str(tmp_path))
        # The subprocess is the I/O being stubbed; it reports the frozen
        # dist-info version, exactly as it does on a real editable install.
        monkeypatch.setattr(sc, "_query_installed_version",
                            lambda timeout=8.0: "0.0.1")
        saved = dict(sc._INSTALLED_VERSION_CACHE)
        try:
            sc._installed_version_worker()
            assert sc._INSTALLED_VERSION_CACHE["version"] is None
        finally:
            sc._INSTALLED_VERSION_CACHE.update(saved)


# ── AG12a: no double-counted primer collection ─────────────────────────────

class TestPrimerCollectionsCountOnce:
    def test_the_active_mirror_is_not_listed_as_a_second_collection(self):
        h = lambda n, p: sc._agent_invoke(None, n, p)
        for name, body in [("create-primer-collection", {"name": "PCMain"}),
                           ("set-active-primer-collection", {"name": "PCMain"}),
                           ("create-primer", {"name": "p1",
                                              "sequence": "ACGTACGTACGTACGTACGT"})]:
            out, st = h(name, body)
            assert st == 200, (name, st, out)
        out, st = h("list-primer-collections", {})
        assert st == 200
        assert sum(c["n_primers"] for c in out["primer_collections"]) == 1, out


# ── AG12b/c: connection handling ───────────────────────────────────────────

class TestConnectionHandling:
    def test_a_client_that_hangs_up_leaves_no_traceback(self, server, capfd):
        port, token = server
        for _ in range(2):
            s = socket.create_connection(("127.0.0.1", port))
            body = b"{}"
            s.sendall(b"POST /zz-slow-read HTTP/1.1\r\nHost: 127.0.0.1\r\n"
                      + f"Authorization: Bearer {token}\r\n".encode()
                      + b"Content-Type: application/json\r\n"
                      + f"Content-Length: {len(body)}\r\n\r\n".encode() + body)
            s.setsockopt(socket.SOL_SOCKET, socket.SO_LINGER,
                         b"\x01\x00\x00\x00\x00\x00\x00\x00")   # RST on close
            s.close()
        time.sleep(1.0)
        assert "Traceback" not in capfd.readouterr().err

    def test_idle_connections_do_not_get_a_thread_each(self, server):
        port, token = server
        base = threading.active_count()
        socks = [socket.create_connection(("127.0.0.1", port)) for _ in range(120)]
        time.sleep(0.8)
        grown = threading.active_count() - base
        for s in socks:
            s.close()
        assert grown <= 70, grown


# ── AG8: the CLI never escapes as a traceback; GET never drops a body ──────

def _tok(tmp_path, monkeypatch, content=b"6701\nt\n", as_dir=False):
    monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(tmp_path))
    f = tmp_path / cli.TOKEN_FILENAME
    if as_dir:
        f.mkdir()
    else:
        f.write_bytes(content)


class _Resp:
    def __init__(self, body):
        self._b = body
        self.status = 200

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def read(self, n=-1):
        return self._b


class TestCliFailuresAreClean:
    @pytest.mark.parametrize("exc", [TimeoutError("timed out"),
                                     http.client.RemoteDisconnected("closed")])
    def test_transport_failures_raise_the_cli_error(self, tmp_path, monkeypatch, exc):
        _tok(tmp_path, monkeypatch)
        with patch("urllib.request.urlopen", side_effect=exc):
            with pytest.raises(cli._AgentCliError):
                cli._request("status")

    def test_a_non_utf8_body_is_not_a_traceback(self, tmp_path, monkeypatch):
        _tok(tmp_path, monkeypatch)
        with patch("urllib.request.urlopen", return_value=_Resp(b"\xff\xfe{}")):
            try:
                cli._request("status")
            except cli._AgentCliError:
                pass

    @pytest.mark.parametrize("kind", ["non_utf8", "directory"])
    def test_a_bad_token_file_exits_with_a_message(self, tmp_path, monkeypatch, kind):
        _tok(tmp_path, monkeypatch, content=b"\xff\xfe\n\xff\n",
             as_dir=(kind == "directory"))
        with pytest.raises(SystemExit) as ei:
            cli._read_session()
        assert isinstance(ei.value.code, str)

    def test_get_with_a_json_body_is_refused_not_dropped(self, tmp_path, monkeypatch, capsys):
        _tok(tmp_path, monkeypatch)
        with patch("urllib.request.urlopen",
                   return_value=_Resp(b'{"count": 20}')) as up:
            with pytest.raises(SystemExit) as ei:
                cli.main(["call", "list-primers", "--method", "GET",
                          "--json", '{"collection": "X"}'])
        assert ei.value.code == 1
        assert up.call_count == 0
