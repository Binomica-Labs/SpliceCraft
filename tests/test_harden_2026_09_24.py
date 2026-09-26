"""Regression tests for the HARDENING pass over the 2026-09-22 audit's residual
fixes (the unreleased tail on top of v1.2.71).

Seven read-only reviewers probed the new code against outside references — the
real Opentrons planner and pipette definitions, Textual's own markup parser,
Biopython and an independent `.dna` parser, brute force — and every finding here was
reproduced before it was fixed. Each test fails on the tree as it stood before
this pass.
"""
from __future__ import annotations

import socket

import pytest

import splicecraft as sc
import splicecraft_opentrons as _ot


# ══════════════════════════════════════════════════════════════════════════
# OT-2 / AUTOLAB
# ══════════════════════════════════════════════════════════════════════════

def _resolve_to(monkeypatch, *addrs):
    """Make every name resolve to ``addrs`` (in order), counting lookups."""
    calls = []

    def fake(host, *a, **k):
        calls.append(host)
        return [(socket.AF_INET6 if ":" in ip else socket.AF_INET,
                 socket.SOCK_STREAM, 6, "", (ip, 0)) for ip in addrs]
    monkeypatch.setattr(socket, "getaddrinfo", fake)
    return calls


@pytest.fixture
def _own_net(monkeypatch):
    """Pin this machine's own IPv4 addresses (and empty the cache)."""
    def pin(*ips):
        monkeypatch.setattr(_ot, "_ot2_host_ipv4s",
                            lambda: [(ip, ip.startswith("169.254.")) for ip in ips])
        monkeypatch.setattr(_ot, "_OT2_OWN_NETS_CACHE", [None, ()],
                            raising=False)
    pin()
    return pin


class TestRobotAddressPolicy:
    def test_an_ipv4_address_written_as_ipv6_is_checked_as_ipv4(self, _own_net):
        # `::ffff:169.254.169.254` connects to the metadata service on a
        # dual-stack socket, and `ipaddress` calls the IPv6 form private.
        for host in ("[::ffff:169.254.169.254]", "[::ffff:a9fe:a9fe]",
                     "::ffff:169.254.169.254"):
            with pytest.raises(_ot.OT2Error, match="metadata"):
                _ot._ot2_base_url(host)

    def test_every_metadata_endpoint_is_refused(self, _own_net):
        for host in ("169.254.170.23", "[fd00:ec2::23]", "192.0.0.192"):
            with pytest.raises(_ot.OT2Error, match="metadata"):
                _ot._ot2_base_url(host)

    def test_a_name_is_connected_at_the_address_that_was_checked(
            self, monkeypatch, _own_net):
        # The HTTP client used to resolve the name AGAIN: a name answering
        # the LAN for the check and the metadata service for the connection
        # passed. The URL is built on the checked address now.
        calls = _resolve_to(monkeypatch, "10.0.0.5")
        url = _ot._ot2_base_url("robot.lab")
        assert url == "http://10.0.0.5:31950"
        assert calls == ["robot.lab"]

    def test_a_local_name_is_resolved_and_checked_too(self, monkeypatch, _own_net):
        _resolve_to(monkeypatch, "169.254.169.254")
        with pytest.raises(_ot.OT2Error, match="metadata"):
            _ot._ot2_base_url("jacques.local")

    def test_a_robot_on_this_machines_public_network_is_reachable(
            self, _own_net):
        # Campus networks hand lab devices public addresses; refusing every
        # public address refused those robots, and discovery dropped them.
        _own_net("128.9.33.5")
        assert _ot._ot2_base_url("128.9.33.77") == "http://128.9.33.77:31950"
        with pytest.raises(_ot.OT2Error, match="public"):
            _ot._ot2_base_url("128.9.34.1")

    def test_the_address_the_user_saved_is_trusted(self, monkeypatch, _own_net):
        # `ot2_host` is written only by AUTOLAB — no agent endpoint sets it.
        import splicecraft_state as _state
        monkeypatch.setattr(_state, "_ot2_trusted_host_hook",
                            lambda: "http://93.184.216.34:31950", raising=False)
        assert _ot._ot2_base_url("93.184.216.34") == "http://93.184.216.34:31950"
        with pytest.raises(_ot.OT2Error, match="public"):
            _ot._ot2_base_url("93.184.216.35")

    def test_a_saved_metadata_address_is_still_refused(self, monkeypatch, _own_net):
        import splicecraft_state as _state
        monkeypatch.setattr(_state, "_ot2_trusted_host_hook",
                            lambda: "169.254.169.254", raising=False)
        with pytest.raises(_ot.OT2Error, match="metadata"):
            _ot._ot2_base_url("169.254.169.254")

    def test_an_empty_port_is_no_port(self, _own_net):
        assert _ot._ot2_base_url("192.168.1.56:") == "http://192.168.1.56:31950"
        assert _ot._ot2_base_url("192.168.1.56:8080") == "http://192.168.1.56:8080"

    @pytest.mark.parametrize("host", ["192.168.1.56:0", "192.168.1.56:65536",
                                      "192.168.1.56:80a"])
    def test_a_port_out_of_range_is_refused(self, host, _own_net):
        with pytest.raises(_ot.OT2Error, match="port"):
            _ot._ot2_base_url(host)

    def test_v1_2_71_spellings_still_work(self, monkeypatch, _own_net):
        # A trailing dot and an underscore in a name — both accepted before.
        assert _ot._ot2_base_url("192.168.1.56.") == "http://192.168.1.56:31950"
        _resolve_to(monkeypatch, "192.168.1.40")
        assert _ot._ot2_base_url("ot2_lab.local") == "http://192.168.1.40:31950"


class TestCustomLabwareHeight:
    def test_the_form_requires_a_height(self):
        import asyncio
        from textual.widgets import Input

        async def run():
            app = sc.PlasmidApp()
            async with app.run_test() as pilot:
                app.action_open_autolab()
                await pilot.pause()
                scr = app.screen
                scr.query_one("#autolab-lw-name", Input).value = "Rack"
                scr._create_labware()
                missing = scr._find_custom_labware_def("Rack")
                scr.query_one("#autolab-lw-depth", Input).value = "37.9"
                scr.query_one("#autolab-lw-height", Input).value = "79.85"
                scr._create_labware()
                return missing, scr._find_custom_labware_def("Rack")
        missing, made = asyncio.run(run())
        assert missing is None
        assert {w["z"] for w in made["wells"].values()} == {41.95}


class TestMultichannelRows:
    def _plan(self, labware, wells, *, pipette="p300_multi_gen2", definition=None):
        dst = {"slot": 2}
        dst.update({"definition": definition} if definition
                   else {"labware": labware})
        return {"pipette": pipette,
                "tips": {"labware": "tiprack_300", "slot": 8},
                "labware": {"src": {"labware": "nest_12_reservoir_15ml", "slot": 1},
                            "dst": dst},
                "steps": [{"type": "transfer", "from": "src:A1", "to": f"dst:{w}",
                           "volume": 20} for w in wells]}

    def test_row_b_of_a_384_plate_is_addressable(self):
        # 4.5 mm rows under 9 mm channels: B1 reaches B, D, … P.
        v = _ot._ot2_validate_plan(self._plan("corning_384_wellplate_112ul_flat",
                                              ["A1", "B1", "B24"]))
        assert not [e for e in v["errors"] if "channels" in e], v["errors"]

    def test_row_c_of_a_384_plate_is_not(self):
        v = _ot._ot2_validate_plan(self._plan("corning_384_wellplate_112ul_flat",
                                              ["C1"]))
        assert any("C1" in e and "channels" in e for e in v["errors"])

    def test_a_custom_16_row_plate_is_row_a_only(self):
        # Opentrons decides by the declared FORMAT: a custom definition is
        # `irregular`, so its B3 fails analysis (a transfer) or is silently
        # skipped (a distribute) — the row count alone is not enough
        # (round-2 hardening, 2026-09-25).
        d = _ot._ot2_build_labware_def("Plate 384", 16, 24, spacing=4.5,
                                       depth=11.0, z_dim=14.2)
        v = _ot._ot2_validate_plan(self._plan(None, ["B3"], definition=d))
        assert any("B3" in e and "channels" in e for e in v["errors"])

    def test_a_custom_definition_declared_384_standard_reaches_row_b(self):
        d = _ot._ot2_build_labware_def("Plate 384", 16, 24, spacing=4.5,
                                       depth=11.0, z_dim=14.2)
        d["parameters"]["format"] = "384Standard"
        v = _ot._ot2_validate_plan(self._plan(None, ["B3"], definition=d))
        assert not [e for e in v["errors"] if "channels" in e], v["errors"]

    def test_a_96_plate_is_still_row_a_only(self):
        v = _ot._ot2_validate_plan(self._plan("corning_96_wellplate_360ul_flat",
                                              ["B1"]))
        assert any("B1" in e and "channels" in e for e in v["errors"])

    def test_only_the_compiled_list_is_checked(self):
        # A steps plan ignores a legacy `transfers` list it also carries, so
        # that list's wells cannot make the plan invalid.
        plan = self._plan("corning_96_wellplate_360ul_flat", ["A1"])
        plan["transfers"] = [{"from": "src:A1", "to": "dst:C3", "volume": 20}]
        v = _ot._ot2_validate_plan(plan)
        assert not [e for e in v["errors"] if "channels" in e], v["errors"]


class TestOneTripSteps:
    """Opentrons' planner splits a distribute / consolidate against the
    PIPETTE's maximum but loads trips against the TIP's, and ends the step
    when one piece does not fit (`if not asp_grouped: break`, opentrons 9.1.2):
    every well from there on silently gets nothing. A transfer splits against
    the tip and is fine."""

    def _plan(self, step_type, volume, rack):
        step = {"type": step_type, "volume": volume}
        if step_type == "distribute":
            step.update({"from": "src:A1", "to": ["dst:A1", "dst:A2", "dst:A3"]})
        elif step_type == "consolidate":
            step.update({"from": ["src:A1", "src:A2"], "to": "dst:A1"})
        else:
            step.update({"from": "src:A1", "to": "dst:A1"})
        return {"pipette": "p300_single_gen2",
                "tips": {"labware": rack, "slot": 8},
                "labware": {"src": {"labware": "eppi_24", "slot": 1},
                            "dst": {"labware": "plate_96", "slot": 2}},
                "steps": [step]}

    def _errors(self, *a):
        return [e for e in _ot._ot2_validate_plan(self._plan(*a))["errors"]
                if "will not fit" in e]

    def test_a_distribute_that_cannot_fit_one_tip_is_refused(self):
        # 190 µL + the 20 µL disposal volume > a 200 µL tip.
        assert self._errors("distribute", 190, "opentrons_96_filtertiprack_200ul")
        assert not self._errors("distribute", 180,
                                "opentrons_96_filtertiprack_200ul")

    def test_a_consolidate_that_cannot_fit_one_tip_is_refused(self):
        assert self._errors("consolidate", 250, "opentrons_96_filtertiprack_200ul")
        assert not self._errors("consolidate", 200,
                                "opentrons_96_filtertiprack_200ul")

    def test_full_size_tips_and_transfers_are_unaffected(self):
        assert not self._errors("distribute", 250, "tiprack_300")
        assert not self._errors("transfer", 250, "opentrons_96_filtertiprack_200ul")


def test_the_position_check_writes_wells_as_the_labware_names_them():
    d = _ot._ot2_build_labware_def("Block", 2, 3, depth=10.0, z_dim=20.0)
    d = {**d, "wells": {k.lower(): v for k, v in d["wells"].items()}}
    plan = {"pipette": "p300_single",
            "tips": {"labware": "tiprack_300", "slot": 8},
            "labware": {"blk": {"slot": 1, "definition": d},
                        "p": {"labware": "plate_96", "slot": 2}}}
    proto = _ot._ot2_compile_position_check(plan)
    # The defaults are upper-cased by `_ot2_entry_wells`; the definition's own
    # keys are lower-case, so `lw_blk["A1"]` named a well that does not exist.
    assert 'lw_blk["a1"]' in proto and 'lw_blk["A1"]' not in proto
    proto = _ot._ot2_compile_position_check(plan, wells=["A01"])
    assert 'lw_p["A1"]' in proto and 'lw_blk["a1"]' in proto


# ══════════════════════════════════════════════════════════════════════════
# Agent API
# ══════════════════════════════════════════════════════════════════════════

def _handler(name):
    return sc._state._AGENT_HANDLERS[name][0]


class TestRoutingGuardReadsHelpers:
    """The central guard refuses a routing key a write handler provably never
    reads. It took `f(payload, "const")` to mean "reads the key `const`" for
    EVERY helper — so `_agent_active_collection_only(payload,
    "delete-from-library")`, which reads `collection`, made three endpoints
    refuse the very `collection` they accept."""

    @pytest.mark.parametrize("endpoint", ["delete-from-library",
                                          "set-plasmid-status",
                                          "rename-plasmid"])
    def test_a_collection_the_endpoint_accepts_is_not_refused(self, endpoint):
        assert sc._agent_unread_routing_keys(
            _handler(endpoint), {"name": "p", "collection": "Main"}) == []

    def test_a_key_nothing_reads_is_still_refused(self):
        assert sc._agent_unread_routing_keys(
            _handler("create-parts-bin"), {"name": "B", "collection": "X"}) \
            == ["collection"]

    @pytest.mark.parametrize("value", [None, "", [], {}, False])
    def test_an_empty_routing_value_is_inert(self, value):
        # Typed clients and tool calls send optional fields as null / "".
        assert sc._agent_unread_routing_keys(
            _handler("create-parts-bin"), {"name": "B", "collection": value,
                                           "enzymes": value}) == []

    def test_an_endpoint_reading_keys_in_a_constant_loop_is_analysed(self):
        # `for _syn in ("commit", "write", "persist"): payload.get(_syn)` made
        # the whole endpoint unknown, so its routing keys went unchecked.
        assert sc._agent_payload_keys(_handler("transfer-annotations")) \
            is not None
        assert sc._agent_unread_routing_keys(
            _handler("transfer-annotations"),
            {"source_id": "pSrc", "collection": "Archive"}) == ["collection"]
        assert sc._agent_unread_routing_keys(
            _handler("transfer-annotations"),
            {"source_id": "pSrc", "source_collection": "Archive"}) == []

    def test_an_f_string_key_is_read_as_its_shape(self):
        # `payload.get(f"{which}_collection")` reads vector_/insert_collection,
        # never a bare `collection` — which traditional-clone ignored,
        # saving the product into the ACTIVE collection.
        fn = _handler("traditional-clone")
        assert sc._agent_payload_keys(fn) is not None
        assert sc._agent_unread_routing_keys(
            fn, {"vector_name": "v", "collection": "Archive"}) == ["collection"]

    def test_a_batch_item_gets_the_same_check(self):
        # `domesticate-parts` calls the single-part handler directly, so the
        # dispatcher never saw an item's `enzyme`.
        from tests.test_agent_api import MockApp
        r = sc._h_domesticate_parts(MockApp(), {"parts": [
            {"sequence": "ATGAAACGTTAA", "oh5": "GGAG", "oh3": "AATG",
             "name": "p1", "enzyme": "BsmBI"}]})
        body = r[0] if isinstance(r, tuple) else r
        assert body["results"][0]["code"] == 400
        assert body["results"][0]["unsupported"] == ["enzyme"]


def test_a_comma_list_names_several_enzymes_when_listing_sites():
    rec = sc._make_demo_record() if hasattr(sc, "_make_demo_record") else None
    if rec is None:
        pytest.skip("no demo record")
    import splicecraft_agent as _ag
    class _App:
        _current_record = rec
    fn = _handler("list-restriction-sites")
    both = fn(_App(), {"enzymes": "EcoRI,BamHI", "unique_only": False})
    assert not isinstance(both, tuple), both
    names = {s.get("enzyme") for s in both.get("sites", [])}
    assert names <= {"EcoRI", "BamHI"} and names
    assert _ag is not None


def test_one_request_hashes_the_same_however_its_numbers_are_spelled():
    assert sc._agent_json_canonical(
        {"tm": 60.0, "x": [6e1, 1.5, {"n": 3.0}], "s": "60.0"}) == \
        {"tm": 60, "x": [60, 1.5, {"n": 3}], "s": "60.0"}


class TestAgentHttp:
    """Real HTTP against `_start_agent_api` (the harness `test_agent_api`
    uses)."""

    def _serve(self, tmp_path, monkeypatch):
        from tests.test_agent_api import MockApp, _free_port
        monkeypatch.setattr(sc._state, "_AGENT_TOKEN_FILE",
                            tmp_path / "agent_token")
        srv = sc._start_agent_api(MockApp(), port=_free_port())
        assert srv is not None
        port, token = (tmp_path / "agent_token").read_text().split()
        return srv, int(port), token

    def _post(self, port, token, path, body: bytes, headers=None):
        import http.client
        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
        h = {"Authorization": f"Bearer {token}",
             "Content-Type": "application/json"}
        h.update(headers or {})
        conn.request("POST", path, body=body, headers=h)
        r = conn.getresponse()
        data = r.read()
        conn.close()
        return r.status, data

    def test_a_transfer_coding_other_than_chunked_is_refused(
            self, tmp_path, monkeypatch):
        # `identity` / `gzip` with no Content-Length: the body was dropped and
        # the endpoint ran with `{}`.
        import http.client
        srv, port, token = self._serve(tmp_path, monkeypatch)
        try:
            conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
            conn.putrequest("POST", "/create-primer")
            conn.putheader("Authorization", f"Bearer {token}")
            conn.putheader("Content-Type", "application/json")
            conn.putheader("Transfer-Encoding", "identity")
            conn.endheaders()
            conn.send(b'{"name": "p", "sequence": "ACGTACGTAC"}')
            r = conn.getresponse()
            status = r.status
            r.read()
            conn.close()
        finally:
            sc._stop_agent_api(srv)
        assert status == 411
        assert not any(p.get("name") == "p" for p in sc._load_primers())

    def test_a_retry_that_respells_a_number_is_replayed(
            self, tmp_path, monkeypatch):
        srv, port, token = self._serve(tmp_path, monkeypatch)
        try:
            h = {"X-Idempotency-Key": "k-respell-1"}
            a = self._post(port, token, "/create-primer",
                           b'{"name": "pq", "sequence": "ACGTACGTAC", "tm": 60}', h)
            b = self._post(port, token, "/create-primer",
                           b'{"name": "pq", "sequence": "ACGTACGTAC", "tm": 6e1}', h)
        finally:
            sc._stop_agent_api(srv)
        assert a[0] == 200 and b[0] == 200, (a, b)
        assert sum(p.get("name") == "pq" for p in sc._load_primers()) == 1


def test_a_guest_leaves_another_guests_token_alone(tmp_path, monkeypatch):
    # Every `--read-only` guest publishes `agent_token.readonly`: the second
    # overwrote the first's, and the first to quit deleted the file the other
    # was still serving from.
    from tests.test_agent_api import MockApp, _free_port
    tok = tmp_path / "agent_token.readonly"
    monkeypatch.setattr(sc._state, "_AGENT_TOKEN_FILE", tok)
    monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY", True)
    first = sc._start_agent_api(MockApp(), port=_free_port())
    second = sc._start_agent_api(MockApp(), port=_free_port())
    try:
        assert first is not None and second is not None
        theirs = tok.read_text()
        sc._stop_agent_api(first)
        first = None
        assert tok.exists() and tok.read_text() == theirs
    finally:
        if first is not None:
            sc._stop_agent_api(first)
        sc._stop_agent_api(second)
    assert not tok.exists()


# ══════════════════════════════════════════════════════════════════════════
# Data safety
# ══════════════════════════════════════════════════════════════════════════

def _ent(name, seed=0):
    import random
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    s = "".join(random.Random(seed or abs(hash(name)) % 10**6).choice("ACGT")
                for _ in range(120))
    gb = sc._record_to_gb_text(SeqRecord(
        Seq(s), id=name, name=name, description=name,
        annotations={"molecule_type": "DNA", "topology": "circular"}))
    return {"id": name, "name": name, "size": 120, "n_feats": 0, "gb_text": gb}


def _colls():
    import json
    raw = json.loads(sc._state._COLLECTIONS_FILE.read_text())["entries"]
    return {c["name"]: sorted(p["name"] for p in c.get("plasmids") or [])
            for c in raw}


def _lib():
    import json
    raw = json.loads(sc._state._LIBRARY_FILE.read_text())["entries"]
    return sorted(p["name"] for p in raw)


def _relaunch():
    """The launch sequence `compose()` runs, from cold caches."""
    st = sc._state
    for attr in ("_library_cache", "_collections_cache", "_settings_cache",
                 "_primers_cache", "_primer_collections_cache",
                 "_parts_bin_cache", "_parts_bin_collections_cache",
                 "_experiments_cache", "_experiment_projects_cache"):
        if hasattr(st, attr):
            setattr(st, attr, None)
    sc._ensure_default_collection()
    sc._restore_library_from_active_collection()
    sc._ensure_default_primer_collection()
    sc._restore_primers_from_active_primer_collection()
    sc._ensure_default_parts_bin()
    sc._restore_parts_bin_from_active_bin()


class _Rig:
    _unsaved = False

    def call_from_thread(self, fn, *a, **k):
        pass


def _agent(name, payload):
    return sc._state._AGENT_HANDLERS[name][0](_Rig(), payload)


def _settings_backup(pred):
    import json
    rows = _agent("list-backups", {"label": "settings"})["backups"]
    for r in rows:
        raw = json.loads(open(r["source_path"]).read())
        ents = raw["entries"] if isinstance(raw, dict) else raw
        if pred({e["key"]: e["value"] for e in ents}):
            return r["source_path"]
    raise AssertionError("no matching settings backup")


@pytest.fixture
def _sync_saves(monkeypatch):
    monkeypatch.setattr(sc, "_collection_sync_force_sync", True)
    monkeypatch.setenv("SPLICECRAFT_SKIP_SETTINGS_FLUSH", "1")


class TestSettingsRestoreKeepsEveryContainer:
    """Restoring settings.json moves every mirror's ACTIVE pointer at once."""

    def test_work_a_failed_mirror_left_behind_survives(self, monkeypatch,
                                                         _sync_saves):
        # The library held a plasmid its collection never got (a failed
        # mirror write). Restoring settings re-staged the library from the
        # stale collection and cleared the marker: the plasmid was gone from
        # every live file.
        import splicecraft_dataaccess as DA
        sc._save_collections([{"name": "Main", "plasmids": [_ent("M1")]}])
        sc._activate_collection("Main")
        sc._set_setting("show_restr", True)
        sc._set_setting("show_restr", False)
        orig, hit = DA._safe_save_json, []

        def flaky(path, *a, **k):
            if path == sc._state._COLLECTIONS_FILE and not hit:
                hit.append(1)
                raise OSError(28, "No space left on device")
            return orig(path, *a, **k)
        monkeypatch.setattr(DA, "_safe_save_json", flaky)
        with pytest.raises(OSError):
            sc._save_library(sc._load_library() + [_ent("UNSYNCED")])
        monkeypatch.setattr(DA, "_safe_save_json", orig)
        assert sc._mirror_is_dirty()
        src = _settings_backup(lambda kv: kv.get("active_collection") == "Main")
        r = _agent("restore-backup", {"label": "settings", "source_path": src})
        assert not isinstance(r, tuple), r
        _relaunch()
        sc._save_library(sc._load_library())
        assert "UNSYNCED" in set(_lib()) | {n for v in _colls().values()
                                            for n in v}

    def test_a_restored_empty_pointer_does_not_overwrite_the_first_collection(
            self, _sync_saves):
        sc._save_settings({"active_collection": ""})
        sc._save_collections([{"name": "Main", "plasmids": [_ent("M1")]},
                              {"name": "Other", "plasmids": [_ent("O1")]}])
        sc._state._settings_cache = None
        sc._activate_collection("Other")
        src = _settings_backup(lambda kv: kv.get("active_collection", "?") == "")
        _agent("restore-backup", {"label": "settings", "source_path": src})
        _relaunch()
        sc._save_library(sc._load_library() + [_ent("NEW")])
        c = _colls()
        assert "M1" in c["Main"] and c["Other"] == ["O1"]

    def test_the_primer_pointer_is_restaged_too(self, _sync_saves):
        # The first fix re-staged only the plasmid library: the next primer
        # save wrote the previous collection's primers into the one the
        # restored pointer names.
        import json
        sc._save_primer_collections([
            {"name": "PA", "primers": [{"name": "a1", "sequence": "ACGTACGTAAC"}]},
            {"name": "PB", "primers": [{"name": "b1", "sequence": "TTGGCCAATTG"}]}])
        _agent("set-active-primer-collection", {"name": "PA"})
        _agent("set-active-primer-collection", {"name": "PB"})
        src = _settings_backup(
            lambda kv: kv.get("active_primer_collection") == "PA")
        _agent("restore-backup", {"label": "settings", "source_path": src})
        sc._save_primers(sc._load_primers()
                         + [{"name": "new", "sequence": "GGGCCCAAATT"}])
        raw = json.loads(sc._state._PRIMER_COLLECTIONS_FILE.read_text())
        pc = {c["name"]: sorted(p["name"] for p in c["primers"])
              for c in raw["entries"]}
        assert pc["PA"] == ["a1", "new"] and pc["PB"] == ["b1"]


class TestDeletingTheActiveCollectionFromTheDialog:
    def test_it_is_not_resurrected_as_recovered_plasmids(self, _sync_saves):
        import asyncio
        from textual.widgets import DataTable
        sc._save_collections([{"name": "Keep", "plasmids": [_ent("K1")]},
                              {"name": "Scratch",
                               "plasmids": [_ent("S1"), _ent("S2")]}])
        sc._activate_collection("Scratch")

        async def run():
            app = sc.PlasmidApp()
            async with app.run_test() as pilot:
                modal = sc.CollectionsModal()
                app.push_screen(modal)
                await pilot.pause()
                t = modal.query_one("#coll-table", DataTable)
                keys = [getattr(k, "value", k) for k in t.rows.keys()]
                t.move_cursor(row=keys.index("Scratch"))
                modal._delete(None)
                await pilot.pause()
                app.screen.query_one("#btn-colldel-yes").action_press()
                await pilot.pause()
                app.screen.query_one("#btn-scarydel-yes").action_press()
                await pilot.pause()
                modal.action_cancel()
                await pilot.pause()
        asyncio.run(run())
        assert sc._get_active_collection_name() == "Keep"
        sc._activate_collection("Keep")          # any later switch
        assert _colls() == {"Keep": ["K1"]}


class TestSpillRestoreKeepsSameNamePrimers:
    def test_a_same_name_variant_comes_back(self, _sync_saves):
        # Two primers may share a name. Matching id-less entries on the name
        # alone read the M13F variant as the "twin" of the current M13F and
        # dropped it.
        prim = [{"name": f"p{i}", "sequence": s} for i, s in enumerate(
            ["ACGTACGTAAC", "TTGGCCAATTG", "GGGCCCAAATT", "CATCATCATCA",
             "AGAGAGAGTCT", "TCTCTCTCAGA", "GTGTGTGTACA", "ACACACACTGT"])]
        prim += [{"name": "M13F", "sequence": "GTAAAACGACGGCCAGT"},
                 {"name": "M13F", "sequence": "TGTAAAACGACGGCCAGT"}]
        sc._save_primers(prim)
        sc._save_primers(prim[:3] + [prim[8]])
        before = {(p["name"], p["sequence"]) for p in prim}
        import splicecraft_backup as B
        spills = [r for r in B._list_recoverable_backups(sc._state._PRIMERS_FILE)
                  if "lost_entries" in str(r.get("source_path"))]
        assert spills, "the shrinking save wrote no spill"
        B._restore_from_backup(sc._state._PRIMERS_FILE,
                               __import__("pathlib").Path(
                                   spills[0]["source_path"]), "Primer library")
        import json
        raw = json.loads(sc._state._PRIMERS_FILE.read_text())["entries"]
        after = {(p["name"], p["sequence"]) for p in raw}
        assert before <= after


def test_switching_from_the_default_primer_library_keeps_its_primers(
        _sync_saves):
    import json
    sc._save_primer_collections([{"name": "PA", "primers": []}])
    _agent("set-active-primer-collection", {"name": ""})
    sc._save_primers([{"name": "orph", "sequence": "ACGTTGCAACGT"}])
    _agent("set-active-primer-collection", {"name": "PA"})
    raw = json.loads(sc._state._PRIMER_COLLECTIONS_FILE.read_text())
    held = {p["name"] for c in raw["entries"] for p in c.get("primers") or []}
    assert "orph" in held


def test_a_cancelled_import_is_still_saved(_sync_saves):
    import asyncio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    async def run():
        app = sc.PlasmidApp()
        async with app.run_test() as pilot:
            await pilot.pause()
            app._apply_record(sc._make_demo_record())
            app._unsaved = True                    # edits on the canvas
            rec = SeqRecord(Seq("ACGT" * 60), id="FETCHED1", name="FETCHED1",
                            description="fetched",
                            annotations={"molecule_type": "DNA",
                                         "topology": "circular"})
            app._import_and_persist(rec)
            await pilot.pause()
            assert isinstance(app.screen, sc.UnsavedNavigateModal)
            app.screen.dismiss(None)               # Cancel
            await pilot.pause()
            await pilot.pause()
            return getattr(app._current_record, "id", None)
    shown = asyncio.run(run())
    assert shown != "FETCHED1"                    # the canvas kept the edits
    assert any(e.get("id") == "FETCHED1" for e in sc._load_library())


# ══════════════════════════════════════════════════════════════════════════
# Reads / alignment — judged against the CONSTRUCTION: every read here is the
# reference with a known edit, so the right answer is known exactly.
# ══════════════════════════════════════════════════════════════════════════

def _rnd(n, seed):
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _verify(ref, reads, **kw):
    r = _handler("verify-against-reads")(None, {"reference": ref,
                                                 "reads": reads, **kw})
    assert not isinstance(r, tuple), r
    return r


class TestReadIdentityCountsWhatTheReadSpans:
    """Covered identity skipped EVERY gap column, including a deletion or an
    insertion inside the read — so a read carrying a 300 bp or 1 kb deletion
    was judged on its substitutions alone and verified as a match."""

    REF = _rnd(3000, 7)

    def test_a_full_length_read_with_a_deletion_is_a_mismatch(self):
        r = _verify(self.REF, [self.REF[:1500] + self.REF[1800:]],
                    circular=True)
        assert r["verdict"] == "mismatch" and not r["reads"][0]["passes"]
        # and the consensus reports it, rather than calling the read clean
        assert any(v["type"] == "deletion" and v["length"] == 300
                   for v in r["consensus"]["variants"])

    def test_a_sanger_read_with_an_insertion_is_a_mismatch(self):
        ins = _rnd(100, 8)
        r = _verify(self.REF, [self.REF[1000:1400] + ins + self.REF[1400:1800]],
                    circular=True)
        assert r["verdict"] == "mismatch"

    def test_a_perfect_partial_read_still_matches(self):
        r = _verify(self.REF, [self.REF[1000:1800]], circular=True)
        assert r["verdict"] == "match" and r["reads"][0]["partial"]

    def test_a_read_of_nothing_but_n_is_no_answer(self):
        # It fell back to the BLAST-style identity, which counts an N as a
        # match, and passed.
        r = _verify(self.REF, ["N" * 3000], circular=True)
        assert not r["reads"][0]["passes"] and r["verdict"] != "match"


def test_local_mode_does_not_turn_clipped_ends_into_an_insertion():
    # The read's clipped end was put in the rows against target gaps: a
    # perfect Sanger read across the origin came back with a 250 bp
    # "insertion" at bp 0 and consensus "unconfirmed".
    ref = _rnd(3000, 4)
    r = _verify(ref, [ref[2750:] + ref[:550]], circular=True, mode="local")
    assert r["verdict"] == "match"
    assert not [v for v in r["consensus"]["variants"]
                if v["type"] == "insertion"]


class TestHeterogeneityCountsOnlyWhatEachReadSaw:
    def _hetero(self, **body):
        r = _handler("analyse-read-heterogeneity")(None, body)
        assert not isinstance(r, tuple), r
        return r

    def test_a_no_call_is_not_depth(self):
        # 2 alt, 8 ref, 12 N at bp 500: 2 of the 10 reads that could say is
        # 20% — above the ONT floor. Counting the N reads made it 9%: muted,
        # and the sample "clonal".
        ref = _rnd(1000, 19)
        p = 500
        alt = {"A": "C", "C": "G", "G": "T", "T": "A"}[ref[p]]
        reads = ([ref[:p] + alt + ref[p + 1:]] * 2 + [ref] * 8
                 + [ref[:p] + "N" + ref[p + 1:]] * 12)
        r = self._hetero(reference=ref, reads=reads, circular=False,
                         platform="ont")
        [pos] = [q for q in r["positions"] if q["pos"] == p]
        assert pos["depth"] == 10 and pos["fraction"] == 0.2
        assert r["verdict"] != "clonal"

    def test_an_ambiguity_code_no_read_reached_is_not_a_no_call(self):
        core = _rnd(3000, 12)
        ref = core[:2000] + "NNK" + core[2003:]
        r = self._hetero(reference=ref, reads=[core[:800]] * 10,
                         circular=False)
        assert r["n_no_call_observations"] == 0

    def test_one_reads_unread_end_is_never_another_reads_support(
            self, monkeypatch):
        # The same bases, laid down two ways: read A aligns its first base and
        # deletes TGC inside its extent; read B leaves the first three bases
        # unread (end gaps). Both normalise to (deletion, 0, CTG), and the
        # gate asked "does SOME read's rollup hold this key" — so B's unread
        # start was counted as support. Only the aligner is stubbed, to fix
        # the two placements; the gate under test runs for real.
        x = _rnd(36, 21)
        ref = "CTGC" + x
        rows = [("C---" + x, ref),              # A: a real 3 bp deletion
                ("---C" + x, ref),              # B: bp 0-2 simply unread
                (ref, ref), (ref, ref)]         # two reference reads
        it = iter(rows)

        def fake_pick(read, target, **kw):
            aq, at = next(it)
            return {"aligned_q": aq, "aligned_t": at, "identity_pct": 100.0}
        monkeypatch.setattr(sc._state, "_pick_best_rotation_hook", fake_pick)
        r = self._hetero(reference=ref, reads=["C" + x, "C" + x, ref, ref],
                         circular=False, min_fraction=0.01)
        hit = [q for q in r["positions"] if q["type"] == "deletion"]
        assert hit and hit[0]["alt_reads"] == 1, r["positions"]


def test_a_deletion_across_the_origin_is_one_event_in_every_read():
    # A 1 bp deletion in a homopolymer spanning bp 0 normalised to two
    # positions (the shift stopped at 0), and full-length reads whose own seam
    # fell in the run lost it as "unread".
    import random
    rng = random.Random(9)
    core = ""
    while not core or core[0] == "A" or core[-1] == "A":
        core = "".join(rng.choice("ACGT") for _ in range(1493))
    ref = "AAAA" + core + "AAA"
    mut = "AAA" + core + "AAA"
    reads = [mut[o:] + mut[:o] for o in (0, 2, 700, 1494, 1496)]
    reads += [(mut[o:] + mut[:o])[:500] for o in (1300, 1400)]
    als = [{"result": sc._pick_best_rotation(rd, ref, is_circular=True,
                                             mode="global"),
            "axis": "target", "query_label": f"r{i}"}
           for i, rd in enumerate(reads)]
    summ = sc._multi_read_summary(als, len(ref), circular=True)
    dels = [v for v in summ["variants"] if v["type"] == "deletion"]
    assert len(dels) == 1 and dels[0]["support"] == len(reads)


class TestInversionProbe:
    def test_a_partial_read_is_not_realigned_end_to_end(self, monkeypatch):
        # The unread plasmid scored as "not matching", so the whole alignment
        # was one candidate, realigned 1 kb against 10 kb twice (~0.9 s).
        ref = _rnd(10000, 77)
        r = sc._pairwise_align(ref[4000:5000], ref, mode="local")
        calls = []
        orig = sc._realign_identity_pct
        monkeypatch.setattr(sc, "_realign_identity_pct",
                            lambda *a, **k: calls.append(1) or orig(*a, **k))
        assert sc._alignment_inverted_segments(r["aligned_q"],
                                               r["aligned_t"]) == []
        assert calls == []

    def test_an_origin_inversion_with_a_half_under_min_bp_is_found(self):
        # A 120 bp flip split 24 / 96 across bp 0: the 24 bp half never
        # became a candidate, so nothing was joined and nothing reported.
        import random
        rng = random.Random(1)
        base = "".join(rng.choice("ACGT") for _ in range(4000))
        s, L, k = 1500, 120, 1524
        construct = base[:s] + sc._rc(base[s:s + L]) + base[s + L:]
        res = sc._pick_best_rotation(construct[k:] + construct[:k],
                                     base[k:] + base[:k], is_circular=True,
                                     mode="global")
        truth = {(p - k) % 4000 for p in range(s, s + L)}
        got = set()
        for sg in res["inverted_segments"]:
            got |= set(range(sg["t_start"], sg["t_end"]))
        assert len(truth & got) >= L - 10 and len(got - truth) <= 10


def test_a_fraction_on_a_bin_edge_lands_in_its_own_bin():
    # 0.15 / 0.05 is 2.9999999999999996 in floating point.
    import splicecraft_agent as _ag
    bins = _ag._hetero_shape([0.15] * 4)["bins"]
    assert [b for b in bins if b["n"]] == [{"from": 0.15, "to": 0.2, "n": 4}]


# ══════════════════════════════════════════════════════════════════════════
# Biology — judged against the sigma-70 geometry the scanner itself uses and
# against the codon a /transl_except names.
# ══════════════════════════════════════════════════════════════════════════

_SIGMA70 = "TTGACA" + "GCTAGCTCAGTCCTAGG" + "TATAAT"


class TestArrowlessElementsStillCount:
    """An explicit strand 0 — SpliceCraft's own "no arrow", which `.dna`
    imports now keep — escaped the promoter / terminator default of +1: an
    arrowless promoter drove nothing and an arrowless terminator blocked
    nothing."""

    def _seq(self):
        return (_rnd(200, 15) + _SIGMA70 + _rnd(40, 16) + "ATG"
                + _rnd(297, 17) + "TAA" + _rnd(500, 18))

    def test_an_arrowless_promoter_drives_the_gene_it_points_at(self):
        import splicecraft_seqanalysis as _sa
        seq = self._seq()
        g0 = 200 + len(_SIGMA70) + 40
        r = _sa._map_transcription(seq, [
            {"start": 200, "end": 229, "strand": 0, "type": "promoter",
             "label": "Pdrive"},
            {"start": g0, "end": g0 + 300, "strand": 1, "type": "CDS",
             "label": "gene"}], circular=True, include_predicted=False)
        [gene] = r["genes"]
        assert gene["n_sources"] == 1 and not gene["silent"]
        assert r["promoters"][0]["strand_unknown"] is True
        assert any("no strand" in w for w in r["warnings"])

    def test_an_arrowless_strong_terminator_still_blocks(self):
        import splicecraft_seqanalysis as _sa
        seq = _rnd(3000, 5)
        r = _sa._map_transcription(seq, [
            {"start": 100, "end": 200, "strand": 1, "type": "promoter",
             "label": "Pfwd"},
            {"start": 250, "end": 290, "strand": 0, "type": "terminator",
             "label": "T", "efficiency": "strong"},
            {"start": 300, "end": 900, "strand": 1, "type": "CDS",
             "label": "geneF"}], circular=True, include_predicted=False)
        [src] = r["genes"][0]["reachable_from"]
        assert src["readthrough"] == "blocked"


def test_a_promoter_at_a_linear_3_end_does_not_drive_the_5_end():
    # Its start site fell past the last base and `% n` put it at bp 0.
    import splicecraft_regulatory as _reg
    import splicecraft_seqanalysis as _sa
    s = _rnd(600, 9) + _SIGMA70 + _rnd(3, 10)
    hit = next(h for h in _reg._scan_promoters(s, circular=False)
               if h["strand"] == 1 and h["minus35_start"] == 600)
    assert hit["tss"] is None and hit["end"] <= len(s)
    r = _sa._map_transcription(s, [{"start": 30, "end": 330, "strand": 1,
                                    "type": "CDS", "label": "g"}],
                               circular=False, include_predicted=True)
    assert not [x for x in r["genes"][0]["reachable_from"]
                if x["promoter"].get("start") == 600]


def test_a_recoded_codon_split_across_the_origin_is_honoured():
    # INSDC writes it `pos:join(...)`; only the first range was read, so the
    # declared residue was ignored and still read as a premature stop.
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    n = 300
    rec = SeqRecord(Seq(_rnd(n, 55)), id="t", name="t",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"})
    loc = CompoundLocation([FeatureLocation(250, 300, strand=1),
                            FeatureLocation(0, 49, strand=1)])
    rec.features.append(SeqFeature(loc, type="CDS"))
    cp = sc._cds_coding_positions(n, rec.features[0], circular=True)
    codon = cp[48:51]                   # residue 17: 2 bases at the end, 1 at 0
    assert codon == [298, 299, 0]
    q = "(pos:join(299..300,1..1),aa:Sec)"
    assert sc._transl_except_positions([q], cp) == {17}
    assert sc._transl_except_residues([q], cp) == {17: "U"}


# ══════════════════════════════════════════════════════════════════════════
# File formats — read back by Biopython as well as by SpliceCraft
# ══════════════════════════════════════════════════════════════════════════

_TESTS = __import__("pathlib").Path(__file__).resolve().parent


def _rec_with(feats, seq=None, **ann):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq(seq or _rnd(600, 3)), id="t1", name="t1",
                    description="test",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular", **ann})
    rec.features.extend(feats)
    return rec


def _feat(start, end, strand=1, ftype="misc_feature", **quals):
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    return SeqFeature(FeatureLocation(start, end, strand=strand), type=ftype,
                      qualifiers={k: [v] for k, v in quals.items()})


def test_a_value_ending_in_one_line_break_still_reloads():
    # `"abc\n".splitlines()` is ONE piece, so it was never split and the raw
    # newline made the saved entry unreadable ("Problem with feature").
    for note in ("replication origin\n", "cloned from stock\r\n", "x\r"):
        gb = sc._record_to_gb_text(_rec_with([_feat(10, 30, note=note)]))
        rec = sc._gb_text_to_record(gb)
        assert rec.features[0].qualifiers["note"] == [note.strip()]


def test_a_display_name_with_a_token_too_long_to_wrap_keeps_its_tail():
    name = ("pET28a_His6-TEV-SUMO-MBP-GFP-mCherry-3xFLAG-Avi_construct_"
            "2026_09_24_v3 (final)")
    rec = _rec_with([])
    rec._tui_display_name = name
    back = sc._gb_text_to_record(sc._record_to_gb_text(rec))
    assert back._tui_display_name == name


class TestDnaImportedBeforeTheReaderRevision:
    """Re-reading an old sidecar with today's `.dna` reader (arrowless kept
    arrowless, rich text decoded) can never reproduce an entry v1.2.71
    imported, so every such entry looked edited and exported from scratch —
    dropping the reads, primers and notes splice mode keeps."""

    DNA = _TESTS / "FFE 5 ENTRY O2.dna"
    GB = _TESTS / "FFE 5 ENTRY O2.v1.2.71-entry.gb"

    def test_an_unedited_old_import_still_splices(self):
        import splicecraft_fileio as _fio
        entry = {"id": "FFE_5", "name": "FFE 5 ENTRY O2",
                 "gb_text": self.GB.read_text()}
        assert _fio._dna_sidecar_matches_entry(self.DNA.read_bytes(), entry)

    def test_an_edited_old_import_does_not(self):
        import splicecraft_fileio as _fio
        gb = self.GB.read_text().replace("misc_feature    993..998",
                                         "misc_feature    993..999", 1)
        assert "993..999" in gb
        entry = {"id": "FFE_5", "name": "FFE 5 ENTRY O2", "gb_text": gb}
        assert not _fio._dna_sidecar_matches_entry(self.DNA.read_bytes(),
                                                   entry)


def _dna_bytes(rec):
    import splicecraft_fileio as _fio
    return _fio._write_commercialsaas_dna_bytes(rec)


def _load_dna(tmp_path, data, name="t.dna"):
    p = tmp_path / name
    p.write_bytes(data)
    return sc.load_genbank(str(p))


class TestDnaTextSurvivesEveryReader:
    def test_a_lone_angle_bracket_is_written_plain(self, tmp_path):
        # Wrapped as HTML, Biopython — which never decodes an entity — read
        # "<5 kb" back as "&lt;5 kb".
        from Bio import SeqIO
        rec = _rec_with([_feat(10, 30, note="<5 kb"),
                         _feat(40, 60, note="Tm > 60 C"),
                         _feat(70, 90, note="use < 5 ng and > 2 ng")])
        p = tmp_path / "t.dna"
        p.write_bytes(_dna_bytes(rec))
        bio = [f.qualifiers.get("note") for f in SeqIO.read(str(p), sc._BIOPYTHON_DNA_FMT)
               .features if f.type != "source"]
        assert ["<5 kb"] in bio and ["Tm > 60 C"] in bio
        ours = [f.qualifiers.get("note") for f in sc.load_genbank(str(p))
                .features if f.type != "source"]
        assert ["use < 5 ng and > 2 ng"] in ours

    def test_the_description_keeps_its_angle_brackets(self, tmp_path):
        rec = _rec_with([])
        rec.description = "Cloning vector pTest <v2> & more"
        assert _load_dna(tmp_path, _dna_bytes(rec)).description == \
            "Cloning vector pTest <v2> & more"

    def test_the_segment_colour_beats_a_stale_colour_qualifier(self, tmp_path):
        # Our own export writes ApEinfo_fwdcolor; recoloured in the other
        # editor, the qualifier went stale — and overrode the new colour.
        import splicecraft_fileio as _fio
        data = _dna_bytes(_rec_with([_feat(10, 30, label="f1",
                                           ApEinfo_fwdcolor="#ff0000")]))
        out = b""
        for ptype, _ln, payload in _fio._iter_commercialsaas_packets(data):
            if ptype == 0x0A:
                payload = payload.replace(b'color="#ff0000"',
                                          b'color="#0000ff"')
            out += _fio._build_commercialsaas_packet(ptype, payload)
        [f] = [f for f in _load_dna(tmp_path, out).features
               if f.type != "source"]
        assert f.qualifiers["ApEinfo_fwdcolor"] == ["#0000ff"]

    def test_same_span_features_keep_their_own_qualifiers(self, tmp_path):
        # The strand tier was keyed by the FILE's arrow and looked up by
        # Biopython's, so an arrowless feature was paired with a forward twin
        # at the same span: arrowless "site-A" came back carrying site-B's
        # /gene. Written the way the other editor writes it.
        import splicecraft_fileio as _fio
        xml = ('<?xml version="1.0"?><Features nextValidID="2">'
               '<Feature recentID="0" name="site-A" type="misc_feature">'
               '<Segment range="11-70" color="#ff0000" type="standard"/>'
               '<Q name="note"><V text="arrowless A"/></Q></Feature>'
               '<Feature recentID="1" name="site-B" type="misc_feature" '
               'directionality="1">'
               '<Segment range="11-70" color="#0000ff" type="standard"/>'
               '<Q name="note"><V text="forward B"/></Q>'
               '<Q name="gene"><V text="gB"/></Q></Feature></Features>')
        data = (_fio._build_commercialsaas_cookie_packet()
                + _fio._build_commercialsaas_dna_packet("ATGC" * 150,
                                                        circular=True)
                + _fio._build_commercialsaas_packet(0x0A, xml.encode()))
        got = {f.qualifiers["label"][0]: (f.location.strand,
                                          f.qualifiers.get("note"),
                                          f.qualifiers.get("gene"))
               for f in _load_dna(tmp_path, data).features
               if f.type != "source"}
        assert got == {"site-A": (0, ["arrowless A"], None),
                       "site-B": (1, ["forward B"], ["gB"])}


class TestDnaBlockCap:
    def test_the_writer_refuses_what_the_reader_would_refuse(self, monkeypatch):
        import splicecraft_fileio as _fio
        monkeypatch.setattr(_fio, "_COMMERCIALSAAS_PACKET_MAX_XML", 4096)
        rec = _rec_with([_feat(i, i + 5, note="x" * 40) for i in range(0, 500, 6)])
        with pytest.raises(ValueError, match="read back"):
            _dna_bytes(rec)

    def test_an_annotated_genome_sized_block_is_under_the_cap(self):
        # 7,800 CDS with translations is a ~10 MB features block; the cap
        # was 8 MB and refused SpliceCraft's own export.
        import splicecraft_fileio as _fio
        assert _fio._COMMERCIALSAAS_PACKET_MAX_XML >= 32 * 1024 * 1024


# ══════════════════════════════════════════════════════════════════════════
# Markup — judged by TEXTUAL's own parser, which is what renders it
# ══════════════════════════════════════════════════════════════════════════

def _render(markup: str) -> str:
    from textual.content import Content
    return Content.from_markup(markup).plain


class TestMarkupSafetyNetKnowsTextualStyles:
    """The net classified tags with RICH's style parser: it rejected the
    theme's `$variables` (so `[b $error]⚠ PHYSICAL ROBOT MOTION[/]` showed its
    markup) and lower-cased style words Textual does not (so `[R]` was kept
    as "reverse" and swallowed)."""

    @pytest.mark.parametrize("msg,shown", [
        ("[b $error]Robot moves[/]", "Robot moves"),
        ("[$accent]Selected:[/] gel.png", "Selected: gel.png"),
        ("[green]Imported Primer [R][/green]", "Imported Primer [R]"),
        ("[green]Imported [[x]][/green]", "Imported [[x]]"),
        ("[bold]Name:[/bold] a[b", "Name: a[b"),
        ("Open [@click=app.quit]quit[/]", "Open [@click=app.quit]quit[/]"),
    ])
    def test_notifications_keep_styles_and_show_names(self, msg, shown):
        fixed = sc._markup_escape_unstyled(msg)
        assert fixed is not None and _render(fixed) == shown

    def test_a_bracketed_name_at_the_end_of_a_static_stays(self):
        # An empty trailing span is dropped by Textual, so judging the parsed
        # spans saw nothing wrong and `[v2]` vanished from the dialog.
        from textual import visual
        assert visual.visualize(None, "Delete pUC19 [v2]").plain == \
            "Delete pUC19 [v2]"


def test_the_notebook_header_uses_the_theme_accent():
    # `[accent]` is no style at all in Textual (it needs the `$`).
    header = sc.ExperimentsScreen._format_project_label(None)
    assert "[$accent]" in header and "[accent]" not in header


# ══════════════════════════════════════════════════════════════════════════
# Notebook
# ══════════════════════════════════════════════════════════════════════════

class TestNotebookExportLinks:
    def test_a_link_with_parentheses_is_a_link(self):
        import splicecraft_experiments as E
        html = E._markdown_subset_to_html(
            "see [the paper](https://doi.org/10.1016/0378-1119(85)90120-9)")
        assert 'href="https://doi.org/10.1016/0378-1119(85)90120-9"' in html

    def test_a_long_embedded_image_still_embeds(self):
        # A 2,048-character cap broke every `data:` image: an exported
        # notebook pasted back as a body showed a wall of base64.
        import splicecraft_experiments as E
        uri = "data:image/png;base64," + "A" * 6000
        assert f'<img src="{uri}"' in E._markdown_subset_to_html(
            f"![gel]({uri})")

    def test_a_data_file_in_an_image_tag_is_a_download_link(self):
        import splicecraft_experiments as E
        e = {"id": "e1", "title": "T", "body_md": "trace: ![run.csv](run.csv)",
             "image_paths": ["run.csv"], "tags": [],
             "created_at": "2026-09-24T10:00:00",
             "updated_at": "2026-09-24T10:00:00"}
        html = E._experiment_html_document(
            [e], image_srcs={"run.csv": "data:text/csv;base64,QQ=="})
        assert '<img src="data:text/csv' not in html
        assert 'download="run.csv"' in html


def test_a_bare_string_tag_is_kept():
    # Read as NO tags, the field showed empty and the next save wrote [] —
    # the tag was deleted.
    out = sc._normalise_experiment_entry(
        {"id": "exp-a1", "title": "t", "body_md": "", "tags": "gibson, pcr"},
        fresh=True)
    assert out["tags"] == ["gibson", "pcr"]


def test_a_reference_that_names_no_plasmid_opens_none(_sync_saves):
    # The fallback search matches a SUBSEQUENCE, and its lone hit was taken:
    # `@pY` opened pACYC184.
    def ent(name, eid):
        return {"name": name, "id": eid, "gb_text": "", "size": 10,
                "n_feats": 0}
    sc._save_collections([{"name": "Default", "plasmids": [
        ent("pACYC184", "pACYC184"), ent("35S-GFP", "35S-GFP")]}])
    sc._set_active_collection_name("Default")
    for ref in ("pY", "A184", "GFP"):
        chosen, _ = sc._resolve_plasmid_ref(ref)
        assert chosen is None, ref
    chosen, _ = sc._resolve_plasmid_ref("pacyc184")
    assert chosen == {"collection": "Default", "id": "pACYC184"}


# ══════════════════════════════════════════════════════════════════════════
# Gels
# ══════════════════════════════════════════════════════════════════════════

def test_an_infinite_pcr_size_is_dropped_not_a_crash():
    g = sc._normalise_gel_entry({"name": "g", "lanes": [
        {"name": "p", "source": "pcr", "detail": "", "_pcr_bp": float("inf")},
        {"name": "q", "source": "pcr", "detail": "", "_pcr_bp": 812}]},
        fresh=True)
    assert "_pcr_bp" not in g["lanes"][0] and g["lanes"][1]["_pcr_bp"] == 812


def test_the_gel_library_restores_pcr_lane_sizes(_sync_saves):
    import asyncio
    from textual.app import App
    g = sc._normalise_gel_entry({"id": "gel-aaaa1111", "name": "Two",
                                 "agarose_pct": 1.0, "lanes": [
        {"name": "ladder", "source": "ladder", "detail": ""},
        {"name": "amp A", "source": "pcr", "detail": "", "_pcr_bp": 812},
        {"name": "amp B", "source": "pcr", "detail": "", "_pcr_bp": 2450}]},
        fresh=True)
    sc._save_gels([g])

    class _A(App):
        def on_mount(self):
            self.push_screen(sc.SimulatorScreen(template_seq="ACGT" * 500,
                                                plasmid_name="t"))

    async def run():
        app = _A()
        async with app.run_test(size=(200, 60)) as pilot:
            await pilot.pause()
            scr = app.screen
            got = {}
            orig = app.push_screen
            app.push_screen = lambda screen, callback=None, **k: \
                got.setdefault("cb", callback)
            scr._on_gel_library(None)
            app.push_screen = orig
            got["cb"]("gel-aaaa1111")
            await pilot.pause()
            return [ln.get("_pcr_bp") for ln in scr._lanes]
    assert asyncio.run(run()) == [None, 812, 2450]
