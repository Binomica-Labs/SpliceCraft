"""Regression tests for the SECOND hardening round (2026-09-25), over the
unreleased audit fixes and the first hardening pass (`[INV-207]`, `[INV-208]`).

Seven read-only reviewers attacked that code against outside references: the
real opentrons 9.1.2 planner and shared-data definitions, Textual's own parser,
markdown-it, Biopython and an independent `.dna` parser, brute force. Every finding here was
reproduced before it was fixed, and each test fails on the tree as it stood
before this round (the controls pass on both).
"""
from __future__ import annotations

import asyncio
import threading
import time

import pytest

import splicecraft as sc


# ══════════════════════════════════════════════════════════════════════════
# Data safety
# ══════════════════════════════════════════════════════════════════════════

class TestToastsNeverWaitOnTheUi:
    """A thread holding `_cache_lock` must never wait for the UI thread: the
    UI thread may be waiting for that same lock, and then both wait forever.
    An orphan-rescue toast sent from an agent thread mid-switch did exactly
    that (`call_from_thread` blocks until the callback has run)."""

    def test_a_rescue_toast_sent_under_the_lock_does_not_stall_the_ui(self):
        got: dict = {}

        async def run():
            app = sc.PlasmidApp()
            async with app.run_test() as pilot:
                await pilot.pause()
                holding = threading.Event()

                def worker():
                    with sc._state._cache_lock:
                        holding.set()
                        time.sleep(0.3)       # the UI thread is now waiting
                        sc._announce_recovered_orphans(
                            "kept 1 plasmid in Recovered plasmids")

                t = threading.Thread(target=worker, daemon=True)
                t.start()
                assert holding.wait(5)
                t0 = time.monotonic()
                # The UI thread wanting the lock, as any save or toggle does.
                # A timeout keeps the old code from hanging the suite.
                got["ok"] = sc._state._cache_lock.acquire(timeout=3)
                got["waited"] = time.monotonic() - t0
                if got["ok"]:
                    sc._state._cache_lock.release()
                t.join(5)
                for _ in range(5):
                    await pilot.pause()
                got["toasts"] = [n.message for n in app._notifications]
        asyncio.run(run())
        assert got["ok"], (
            f"the UI thread waited {got['waited']:.1f}s for a lock held by a "
            f"thread that was waiting for the UI")
        assert any("Recovered plasmids" in m for m in got["toasts"])

    def test_a_background_save_failure_toast_is_posted_too(self):
        got: dict = {}

        async def run():
            app = sc.PlasmidApp()
            async with app.run_test() as pilot:
                await pilot.pause()
                holding = threading.Event()

                def worker():
                    with sc._state._cache_lock:
                        holding.set()
                        time.sleep(0.3)
                        sc._bg_notify_save_failure(
                            "Active-parts-bin mirror (Main)",
                            OSError("disk full"))

                t = threading.Thread(target=worker, daemon=True)
                t.start()
                assert holding.wait(5)
                got["ok"] = sc._state._cache_lock.acquire(timeout=3)
                if got["ok"]:
                    sc._state._cache_lock.release()
                t.join(5)
                for _ in range(5):
                    await pilot.pause()
                got["toasts"] = [n.message for n in app._notifications]
        asyncio.run(run())
        assert got["ok"]
        assert any("disk full" in m for m in got["toasts"])


# ══════════════════════════════════════════════════════════════════════════
# OT-2 / AUTOLAB
# ══════════════════════════════════════════════════════════════════════════

def _serve_agent(tmp_path, monkeypatch):
    from tests.test_agent_api import MockApp, _free_port
    monkeypatch.setattr(sc._state, "_AGENT_TOKEN_FILE", tmp_path / "agent_token")
    srv = sc._start_agent_api(MockApp(), port=_free_port())
    assert srv is not None
    port, token = (tmp_path / "agent_token").read_text().split()
    return srv, int(port), token


def _post_json(port, token, path, body, headers=None):
    import http.client
    import json
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=20)
    h = {"Authorization": f"Bearer {token}", "Content-Type": "application/json"}
    h.update(headers or {})
    conn.request("POST", path, body=json.dumps(body).encode(), headers=h)
    r = conn.getresponse()
    data = r.read()
    conn.close()
    return r.status, json.loads(data or b"{}")


class _Robot:
    """A stubbed OT-2: counts runs created, plays and stops. ``mode`` picks how
    the run goes wrong AFTER it exists on the robot."""

    def __init__(self, monkeypatch, mode):
        import splicecraft_opentrons as ot
        self.runs = self.plays = self.stops = 0
        self.mode = mode
        monkeypatch.setattr(ot, "_OT2_RUN_POLL_TIMEOUT", 0.2)
        monkeypatch.setattr(ot, "_OT2_POLL_INTERVAL", 0.02)
        monkeypatch.setattr(ot, "_ot2_analyze", lambda host, text, **k: {
            "result": "ok", "status": "completed", "errors": [], "commands": [],
            "pipettes": [], "labware": [], "protocol_id": "proto-1"})
        monkeypatch.setattr(ot, "_ot2_state", self.state)
        monkeypatch.setattr(ot, "_ot2_request_json", self.request)
        monkeypatch.setattr(ot, "_ot2_set_lights", lambda host, on: {"on": on})
        monkeypatch.setattr(ot, "_ot2_stop_run", self.stop)
        monkeypatch.setattr(ot, "_ot2_pipette_mismatch", lambda a, b: [])

    def state(self, host, run_id=None, **k):
        if run_id is None:                      # the pre-flight check
            return {"reachable": True, "faults": [], "door": {}, "lights": False,
                    "calibration": {"pipettes_calibrated": {"left": True}},
                    "instruments": []}
        if self.mode == "monitor-error":
            raise KeyError("malformed robot response")
        if self.mode == "not-a-dict":
            return ["garbage"]
        return {"run": {"status": "running"}, "faults": []}

    def request(self, host, path, *, method="GET", payload=None, **k):
        if path == "/runs" and method == "POST":
            self.runs += 1
            return {"data": {"id": f"run-{self.runs}"}}
        if path.endswith("/actions"):
            self.plays += 1
            if self.mode == "play-fails":
                import splicecraft_opentrons as ot
                raise ot.OT2Error("connection reset during play")
            return {"data": {}}
        return {"data": {}}

    def stop(self, host, rid):
        self.stops += 1
        return {}


class TestARunThatStartedIsNeverRetried:
    """Once a run exists on the robot, a failure must come back as a result
    the idempotency cache keeps — never a 5xx, which it does not keep and
    every client retries: the retry created and PLAYED a second run over wells
    the first had already filled."""

    @pytest.mark.parametrize("mode", ["timeout", "monitor-error", "not-a-dict",
                                      "play-fails"])
    def test_a_same_key_retry_does_not_move_the_robot_again(
            self, tmp_path, monkeypatch, mode):
        robot = _Robot(monkeypatch, mode)
        srv, port, token = _serve_agent(tmp_path, monkeypatch)
        body = {"host": "192.168.1.56", "protocol": "metadata = {}\n",
                "confirm": True}
        key = {"X-Idempotency-Key": f"run-once-{mode}"}
        try:
            s1, r1 = _post_json(port, token, "/ot2-run", body, key)
            s2, r2 = _post_json(port, token, "/ot2-run", body, key)
        finally:
            sc._stop_agent_api(srv)
        assert robot.runs == 1 and robot.plays == 1, (s1, r1, s2, r2)
        assert s1 == s2 == 422
        assert robot.stops >= 1                 # never left moving unwatched

    def test_the_crash_says_why_and_that_it_was_stopped(self, monkeypatch):
        robot = _Robot(monkeypatch, "timeout")
        out = sc._state._AGENT_HANDLERS["ot2-run"][0](
            None, {"host": "192.168.1.56", "protocol": "metadata = {}\n",
                   "confirm": True})
        body, status = out
        assert status == 422 and body["ran"] is True and body["crashed"] is True
        assert body["run_id"] == "run-1" and body["ok"] is False
        assert "did not finish" in body["error"]
        assert "was stopped" in body["error"] and robot.stops == 1


def _resolve(monkeypatch, *sockaddrs):
    """Make every name resolve to ``sockaddrs`` — (ip, port) or, for IPv6,
    (ip, port, flowinfo, scope_id) exactly as getaddrinfo returns them."""
    import socket
    monkeypatch.setattr(socket, "getaddrinfo", lambda host, *a, **k: [
        (socket.AF_INET6 if ":" in sa[0] else socket.AF_INET,
         socket.SOCK_STREAM, 6, "", sa) for sa in sockaddrs])


@pytest.fixture
def _no_own_net(monkeypatch):
    import splicecraft_opentrons as ot
    monkeypatch.setattr(ot, "_ot2_host_ipv4s", lambda: [])
    monkeypatch.setattr(ot, "_OT2_OWN_NETS_CACHE", [None, ()], raising=False)
    monkeypatch.setattr(sc._state, "_ot2_trusted_host_hook", lambda: "",
                        raising=False)


class TestRobotNamesWithSeveralAddresses:
    """Each address a name resolves to is checked, and the connection is
    pinned to one that passed — so one address that fails the policy is
    skipped, not a reason to refuse the robot."""

    def test_a_dual_stack_robot_is_reached_on_its_lan_address(
            self, monkeypatch, _no_own_net):
        import splicecraft_opentrons as ot
        _resolve(monkeypatch, ("192.168.1.60", 0),
                 ("2a01:4f8:c0c:1::60", 0, 0, 0))
        assert ot._ot2_base_url("dual.local") == "http://192.168.1.60:31950"

    def test_ipv4_is_preferred_over_a_link_local_ipv6_answer(
            self, monkeypatch, _no_own_net):
        # macOS / Windows sort the AAAA of a .local name first.
        import splicecraft_opentrons as ot
        _resolve(monkeypatch, ("fe80::1ff:fe23:4567", 0, 0, 3),
                 ("192.168.1.56", 0))
        assert ot._ot2_base_url("opentrons.local") == \
            "http://192.168.1.56:31950"

    def test_a_link_local_ipv6_address_keeps_its_zone(
            self, monkeypatch, _no_own_net):
        import urllib.parse
        import splicecraft_opentrons as ot
        _resolve(monkeypatch, ("fe80::1ff:fe23:4567", 0, 0, 3))
        url = ot._ot2_base_url("opentrons.local")
        assert url == "http://[fe80::1ff:fe23:4567%3]:31950"
        assert urllib.parse.urlsplit(url).port == 31950

    def test_a_name_with_no_allowed_address_is_still_refused(
            self, monkeypatch, _no_own_net):
        import splicecraft_opentrons as ot
        _resolve(monkeypatch, ("169.254.169.254", 0), ("8.8.8.8", 0))
        with pytest.raises(ot.OT2Error, match="refusing"):
            ot._ot2_base_url("evil.example")

    @pytest.mark.parametrize("host", ["[fd00:ec2::254%1]", "[fd00:ec2::23%lo]"])
    def test_a_zone_does_not_disguise_a_metadata_address(self, host,
                                                        _no_own_net):
        import splicecraft_opentrons as ot
        with pytest.raises(ot.OT2Error, match="metadata"):
            ot._ot2_base_url(host)


def test_this_machines_network_counts_in_the_first_minute_after_boot(tmp_path):
    # `time.monotonic()` counts from boot; a "never computed" stamp of 0.0
    # read as FRESH until uptime passed the 60 s cache life, so the own /24
    # was empty and a public campus-LAN robot was refused right after boot.
    import os
    import subprocess
    import sys
    code = (
        "import time\n"
        "time.monotonic = lambda: 5.0\n"
        "import splicecraft_opentrons as ot\n"
        "ot._ot2_host_ipv4s = lambda: [('128.84.7.20', False)]\n"
        "print(len(ot._ot2_own_networks()))\n")
    env = dict(os.environ, XDG_DATA_HOME=str(tmp_path),
               SPLICECRAFT_SKIP_LOCK="1")
    root = os.path.dirname(os.path.dirname(os.path.abspath(sc.__file__)))
    out = subprocess.run([sys.executable, "-c", code], cwd=os.path.dirname(
        os.path.abspath(sc.__file__)), env=env, capture_output=True,
        text=True, timeout=120)
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.strip().splitlines()[-1] == "1", (root, out.stdout)


class TestRunChecksFailClosed:
    """"Is a run active?" had one answer for "no" and for "could not ask":
    None. A Stop that could not reach the robot said nothing was running and
    sent nothing; Home and Disengage went ahead on the same misreading."""

    @pytest.fixture
    def _unreachable_runs(self, monkeypatch):
        import splicecraft_opentrons as ot
        sent = []

        def try_json(host, path, **k):
            if path == "/runs":
                return None, "timed out"
            return {}, None
        monkeypatch.setattr(ot, "_ot2_try_json", try_json)
        monkeypatch.setattr(ot, "_ot2_run_action",
                            lambda host, rid, action: sent.append(action) or {})
        monkeypatch.setattr(ot, "_ot2_home", lambda host: sent.append("home"))
        monkeypatch.setattr(ot, "_ot2_disengage",
                            lambda host, *a, **k: sent.append("disengage"))
        return sent

    def test_stop_says_it_could_not_ask_not_that_nothing_runs(
            self, _unreachable_runs):
        body, status = sc._state._AGENT_HANDLERS["ot2-run-control"][0](
            None, {"host": "192.168.1.56", "action": "stop"})
        assert status == 502 and "could not ask" in body["error"]
        assert "nothing is running" not in body["error"]

    @pytest.mark.parametrize("ep", ["ot2-home", "ot2-disengage"])
    def test_motion_is_refused_when_the_run_check_cannot_be_made(
            self, _unreachable_runs, ep):
        out = sc._state._AGENT_HANDLERS[ep][0](None, {"host": "192.168.1.56"})
        assert isinstance(out, tuple) and out[1] == 502, out
        assert _unreachable_runs == []            # nothing was sent

    def test_the_library_function_says_so_too(self, _unreachable_runs):
        import splicecraft_opentrons as ot
        with pytest.raises(ot.OT2Error, match="could not ask"):
            ot._ot2_run_control("192.168.1.56", "stop")


class TestNormaliseNeverPipettesIntoAnOccupiedWell:
    """The same-plate check looked only at wells a step READS. A sample the
    normalise skipped, or one left without a destination, is still sitting in
    its well; and a diluent on the destination plate was never checked when
    the plates differed."""

    def _call(self, items, **extra):
        body = {"items": items, "target_ng": 100, "final_volume": 20, **extra}
        return sc._state._AGENT_HANDLERS["ot2-normalize"][0](None, body)

    def test_a_skipped_samples_well_is_not_a_destination(self):
        # s2 has no concentration, so no step reads B1 — it still holds s2.
        items = [{"name": "s1", "well": "A1", "concentration": 50.0},
                 {"name": "s2", "well": "B1", "concentration": None}]
        r = self._call(items, src="plate", dst="plate", dst_wells=["B1"])
        assert isinstance(r, tuple) and r[1] == 400, r
        assert "B1" in r[0]["error"]

    def test_an_unread_samples_well_is_not_a_destination(self):
        items = [{"name": f"s{i}", "well": f"A{i}", "concentration": 50.0}
                 for i in (1, 2, 3)]
        r = self._call(items, src="plate", dst="plate", dst_wells=["A3"])
        assert isinstance(r, tuple) and r[1] == 400, r
        assert "A3" in r[0]["error"]

    def test_a_diluent_on_the_destination_plate_is_not_a_destination(self):
        items = [{"name": "s1", "well": "A1", "concentration": 50.0},
                 {"name": "s2", "well": "A2", "concentration": 80.0}]
        r = self._call(items, src="samples", dst="dstplate",
                       dst_labware="plate_96", diluent_ref="dstplate:A1",
                       target_conc=10, target_ng=None)
        assert isinstance(r, tuple) and r[1] == 400, r
        assert "A1" in r[0]["error"]

    def test_control_a_separate_diluent_reservoir_builds(self):
        items = [{"name": "s1", "well": "A1", "concentration": 50.0},
                 {"name": "s2", "well": "A2", "concentration": 80.0}]
        r = self._call(items, src="samples", dst="dstplate",
                       dst_labware="plate_96", diluent_ref="buf:A1",
                       target_conc=10, target_ng=None)
        assert not isinstance(r, tuple), r
        assert r["steps"]


class TestMultichannelOnCustomLabware:
    """Opentrons lets a multi-channel's first channel sit in row B only on a
    labware whose FORMAT is ``384Standard``. A custom 16-row plate is
    ``irregular`` (the AUTOLAB builder writes that), so the real planner
    skips its row-B wells in a distribute without a word."""

    @staticmethod
    def _plan(step_type, wells, definition):
        step = ({"type": "distribute", "from": "src:A1",
                 "to": [f"dst:{w}" for w in wells], "volume": 20}
                if step_type == "distribute" else
                {"type": "transfer", "from": "src:A1", "to": f"dst:{wells[-1]}",
                 "volume": 20})
        return {"pipette": "p20_multi_gen2",
                "tips": {"labware": "tiprack_20", "slot": 8},
                "labware": {"src": {"labware": "nest_12_reservoir_15ml",
                                    "slot": 1},
                            "dst": {"slot": 2, "definition": definition}},
                "steps": [step]}

    @pytest.mark.parametrize("step_type", ["distribute", "transfer"])
    def test_row_b_of_an_irregular_16_row_plate_is_refused(self, step_type):
        import splicecraft_opentrons as ot
        d = ot._ot2_build_labware_def("My 384", 16, 24, spacing=4.5,
                                      depth=11.0, z_dim=14.4)
        v = ot._ot2_validate_plan(self._plan(step_type, ["A1", "B1"], d))
        assert any("B1" in e and "channels" in e for e in v["errors"]), \
            v["errors"]

    def test_rows_are_counted_once_per_labware_not_per_well(self, monkeypatch):
        # A big plan re-counted a custom definition's rows for every well it
        # named: 30 s to validate a plan the robot analyses in 5.
        import splicecraft_opentrons as ot
        d = ot._ot2_build_labware_def("My 384", 16, 24, spacing=4.5,
                                      depth=11.0, z_dim=14.4)
        calls = []
        real = ot._ot2_labware_rows
        monkeypatch.setattr(ot, "_ot2_labware_rows",
                            lambda entry: calls.append(1) or real(entry))
        wells = [f"A{c}" for c in range(1, 25)] * 10
        ot._ot2_validate_plan(self._plan("distribute", wells, d))
        assert len(calls) <= 2, len(calls)


class TestCustomLabwareFormIsStrict:
    """Every number on the AUTOLAB labware form places wells in space. Height
    was parsed strictly; Depth, Spacing, Rows and Cols fell back to their
    defaults on "38,5" or "38.5 mm" — well z is height − depth, so the tip
    went below the tube bottom, with no message."""

    CASES = [
        ("RackComma", {"#autolab-lw-depth": "38,5"}, "Depth"),
        ("RackUnit", {"#autolab-lw-depth": "38.5 mm"}, "Depth"),
        ("RackPitch", {"#autolab-lw-spacing": "19,5"}, "Spacing"),
        ("RackRows", {"#autolab-lw-rows": "4.5"}, "Rows"),
        ("RackCols", {"#autolab-lw-cols": "six"}, "Cols"),
        ("RackZero", {"#autolab-lw-depth": "0"}, "Depth"),
    ]

    def _run(self, fields):
        from textual.widgets import Input
        got: dict = {}

        async def run():
            app = sc.PlasmidApp()
            notes: list = []
            app.notify = lambda msg, *a, **k: notes.append(str(msg))
            async with app.run_test() as pilot:
                app.action_open_autolab()
                await pilot.pause()
                scr = app.screen
                base = {"#autolab-lw-name": "Rack", "#autolab-lw-height": "80",
                        "#autolab-lw-depth": "38.5", "#autolab-lw-spacing": "19.5",
                        "#autolab-lw-rows": "4", "#autolab-lw-cols": "6"}
                base.update(fields)
                for sel, val in base.items():
                    scr.query_one(sel, Input).value = val
                scr._create_labware()
                got["def"] = scr._find_custom_labware_def(base["#autolab-lw-name"])
                got["notes"] = list(notes)
        asyncio.run(run())
        return got

    @pytest.mark.parametrize("name,fields,label", CASES)
    def test_a_malformed_number_is_refused_by_name(self, name, fields, label):
        got = self._run({"#autolab-lw-name": name, **fields})
        assert got["def"] is None, "saved with a silently defaulted value"
        assert any(label in n for n in got["notes"]), got["notes"]

    def test_control_well_formed_numbers_build_the_labware_as_typed(self):
        got = self._run({"#autolab-lw-name": "RackGood"})
        a1, a2 = got["def"]["wells"]["A1"], got["def"]["wells"]["A2"]
        assert a1["depth"] == 38.5 and a1["z"] == pytest.approx(80 - 38.5)
        assert a2["x"] - a1["x"] == pytest.approx(19.5)


def _opentrons_step_stops_early(step_type, vol, pip_min, pip_max, tip):
    """The opentrons 9.1.2 planner, transcribed: `expand_for_volume_constraints`
    cuts each well's volume against the PIPETTE maximum (less the disposal
    volume for a distribute), then trips are loaded against the working
    volume (the smaller of tip and pipette) and the step ENDS the first time
    a piece does not fit on its own (`if not asp_grouped: break`)."""
    disposal = pip_min if step_type == "distribute" else 0.0
    max_piece = pip_max - disposal

    def expand(volumes):
        for v in volumes:
            while v > max_piece * 2:
                yield max_piece
                v -= max_piece
            if v > max_piece:
                v /= 2
                yield v
            yield v
    work = min(tip, pip_max)
    it = expand([vol] * 3)
    current = next(it)
    done = False
    while not done:
        grouped = []
        try:
            while sum(grouped) + disposal + current <= work:
                grouped.append(current)
                current = next(it)
        except StopIteration:
            done = True
        if not grouped:
            return True
    return False


class TestTipCapacityMatchesTheRealPlanner:
    """The first rule judged a well's WHOLE volume against the tip, and refused
    plans the robot splits and completes — a P300 on 200 µL tips distributes
    300 µL as two 150 µL trips."""

    RACKS = {10.0: "opentrons_96_tiprack_10ul",
             200.0: "opentrons_96_filtertiprack_200ul"}

    def _refused(self, pipette, step_type, vol, rack):
        step = {"type": step_type, "volume": vol}
        if step_type == "distribute":
            step.update({"from": "src:A1",
                         "to": ["dst:A1", "dst:A2", "dst:A3"]})
        else:
            step.update({"from": ["src:A1", "src:A2", "src:A3"],
                         "to": "dst:A1"})
        import splicecraft_opentrons as ot
        plan = {"pipette": pipette, "tips": {"labware": rack, "slot": 8},
                "labware": {"src": {"labware": "eppi_24", "slot": 1},
                            "dst": {"labware": "plate_96", "slot": 2}},
                "steps": [step]}
        return any("will not fit" in e
                   for e in ot._ot2_validate_plan(plan)["errors"])

    @pytest.mark.parametrize("pipette,tip", [("p300_single_gen2", 200.0),
                                             ("p300_single", 200.0),
                                             ("p20_single_gen2", 10.0)])
    @pytest.mark.parametrize("step_type", ["distribute", "consolidate"])
    def test_refused_exactly_when_the_robot_would_stop(self, pipette, tip,
                                                       step_type):
        import splicecraft_opentrons as ot
        pmin, pmax, _ch = ot._OT2_PIPETTES[pipette]
        wrong = []
        for tenth in range(5, int(pmax * 2.5 * 10), 5):
            vol = tenth / 10
            want = _opentrons_step_stops_early(step_type, vol, pmin, pmax, tip)
            if self._refused(pipette, step_type, vol, self.RACKS[tip]) != want:
                wrong.append((vol, want))
        assert not wrong, wrong[:8]

    def test_the_reported_case(self):
        assert not self._refused("p300_single_gen2", "distribute", 300,
                                 self.RACKS[200.0])
        assert not self._refused("p300_single_gen2", "consolidate", 350,
                                 self.RACKS[200.0])
        assert self._refused("p300_single_gen2", "distribute", 250,
                             self.RACKS[200.0])


# ══════════════════════════════════════════════════════════════════════════
# Markup safety
# ══════════════════════════════════════════════════════════════════════════

def _plain_texts(modal):
    """Everything the modal's Statics show, as plain text, after rendering."""
    got: list = []

    async def run():
        from textual.app import App
        from textual.widgets import Static

        class _A(App):
            pass
        app = _A()
        async with app.run_test(size=(140, 50)) as pilot:
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause()
            for st in app.screen.query(Static):
                got.append(st.visual.plain if hasattr(st.visual, "plain")
                           else str(st.visual))
    asyncio.run(run())
    return "\n".join(got)


class TestApprovalDialogsShowWhatTheyApprove:
    """Rich's `escape` leaves `[$var]` alone and Textual reads it as a style:
    a model-supplied ``[$background]`` in an argument painted the rest of the
    request background-on-background in the very dialog asking the user to
    approve it."""

    HIDE = "pA2 [$background]"

    def test_a_single_call_shows_its_arguments(self):
        plain = _plain_texts(sc.BabsToolApprovalModal(
            "delete-from-library", {"names": ["pVeryImportant"],
                                    "note": self.HIDE}))
        assert "[$background]" in plain and "pVeryImportant" in plain

    def test_a_batch_shows_every_call(self):
        plain = _plain_texts(sc.BabsBatchApprovalModal([
            ("rename-plasmid", {"name": "pA", "new_name": self.HIDE}),
            ("delete-collection", {"name": "Main"})]))
        assert "[$background]" in plain
        assert "delete-collection" in plain and "Main" in plain

    def test_the_babs_chat_renderer_escapes_style_variables(self):
        import splicecraft_babs as babs
        from textual.content import Content
        md = babs.md_to_rich("before [$background] hidden after **bold**")
        plain = Content.from_markup(md).plain
        assert "[$background] hidden after" in plain and "bold" in plain


class TestEveryWidgetShowsANameAsWritten:
    """The markup net covered Static/Label. The Select dropdown (OptionList
    prompts) and every label widget (Button, Checkbox, RadioButton, Tab,
    Collapsible) parse markup through other doors, so a plasmid named
    ``pFoo[/b]`` still closed the app from a dropdown. (A name that is itself
    a style — ``[$background]`` — is the documented limit of a net that must
    keep the app's own ``[$accent]``: such text is escaped at its call site.)"""

    NAMES = ["pFoo[/b]", "pUC19 [v2]", "clone [A] [1] [/i] end"]

    def test_a_select_dropdown_with_bracketed_names(self):
        from textual.widgets import Select
        from textual.widgets._option_list import OptionList
        seen: dict = {}

        names = self.NAMES

        async def run():
            from textual.app import App

            class _A(App):
                def compose(self):
                    yield Select([(n, i) for i, n in enumerate(names)],
                                 id="sel")
            app = _A()
            async with app.run_test(size=(120, 40)) as pilot:
                sel = app.query_one("#sel", Select)
                sel.action_show_overlay()
                await pilot.pause()
                await pilot.pause()
                ol = app.query_one(OptionList)
                seen["prompts"] = [
                    str(ol._get_visual(ol.get_option_at_index(i)).plain)
                    if hasattr(ol, "_get_visual") else ""
                    for i in range(ol.option_count)]
                sel.value = 1
                await pilot.pause()
                seen["ok"] = True
        asyncio.run(run())
        assert seen.get("ok")
        shown = " | ".join(seen["prompts"])
        for name in self.NAMES:
            assert name in shown or not shown, (name, shown)

    def test_labels_keep_every_character(self):
        from textual.widgets import Button, Checkbox, RadioButton, Tab
        from textual.widgets import Collapsible
        shown: list = []

        names = self.NAMES

        async def run():
            from textual.app import App

            class _A(App):
                def compose(self):
                    for i, n in enumerate(names):
                        yield Button(n, id=f"b{i}")
                        yield Checkbox(n, id=f"c{i}")
                        yield RadioButton(n, id=f"r{i}")
                        yield Tab(n, id=f"t{i}")
                        yield Collapsible(title=n, id=f"k{i}")
            app = _A()
            async with app.run_test(size=(160, 80)) as pilot:
                await pilot.pause()
                for i, n in enumerate(self.NAMES):
                    shown.append((n, str(app.query_one(f"#b{i}", Button).label)))
                    shown.append((n, str(app.query_one(f"#t{i}", Tab).label)))
        asyncio.run(run())
        for name, label in shown:
            assert name in label, (name, label)


class TestDeletingACollectionFromTheDialog:
    """The Collections dialog deleted a collection — and every plasmid in it —
    on one press with no question, then cleared the canvas even when it held
    unsaved edits."""

    def _open(self, pilot, app):
        from textual.widgets import DataTable
        modal = sc.CollectionsModal()
        app.push_screen(modal)
        return modal

    def test_it_asks_before_deleting(self, _sync_saves):
        from textual.widgets import DataTable
        sc._save_collections([{"name": "Keep", "plasmids": []},
                              {"name": "Drop", "plasmids": []}])
        sc._activate_collection("Keep")
        names: dict = {}

        async def run():
            app = sc.PlasmidApp()
            async with app.run_test() as pilot:
                modal = sc.CollectionsModal()
                app.push_screen(modal)
                await pilot.pause()
                t = modal.query_one("#coll-table", DataTable)
                keys = [getattr(k, "value", k) for k in t.rows.keys()]
                t.move_cursor(row=keys.index("Drop"))
                modal._delete(None)
                await pilot.pause()
                names["asked"] = type(app.screen).__name__
                app.screen.query_one("#btn-colldel-no").action_press()
                await pilot.pause()
        asyncio.run(run())
        assert names["asked"] == "CollectionDeleteConfirmModal"
        assert {c["name"] for c in sc._load_collections()} == {"Keep", "Drop"}

    def test_unsaved_edits_on_the_canvas_are_not_thrown_away(self, _sync_saves):
        from textual.widgets import DataTable
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ACGT" * 30), id="a1", name="a1",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        sc._save_collections([
            {"name": "A", "plasmids": [_entry_for(rec)]},
            {"name": "B", "plasmids": []}])
        sc._activate_collection("A")
        got: dict = {}

        async def run():
            app = sc.PlasmidApp()
            async with app.run_test() as pilot:
                await pilot.pause()
                app._apply_record(rec)
                await pilot.pause()
                app._unsaved = True
                modal = sc.CollectionsModal()
                app.push_screen(modal)
                await pilot.pause()
                t = modal.query_one("#coll-table", DataTable)
                keys = [getattr(k, "value", k) for k in t.rows.keys()]
                t.move_cursor(row=keys.index("A"))
                modal._delete(None)
                await pilot.pause()
                for btn in ("#btn-colldel-yes", "#btn-scarydel-yes"):
                    found = app.screen.query(btn)     # none on the old tree
                    if found:
                        found.first().action_press()
                        await pilot.pause()
                got["record"] = app._current_record
                got["unsaved"] = app._unsaved
        asyncio.run(run())
        assert got["record"] is not None and got["unsaved"] is True
        assert {c["name"] for c in sc._load_collections()} == {"B"}


def _entry_for(rec):
    return {"id": rec.id, "name": rec.name, "size": len(rec.seq),
            "gb_text": sc._record_to_gb_text(rec)}


@pytest.fixture
def _sync_saves(monkeypatch):
    monkeypatch.setattr(sc, "_collection_sync_force_sync", True)
    monkeypatch.setenv("SPLICECRAFT_SKIP_SETTINGS_FLUSH", "1")


# ── restore / switch / rescue ─────────────────────────────────────────────

def _ent(name, seed=None, n=120):
    import random
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rng = random.Random(seed if seed is not None else sum(map(ord, name)))
    s = "".join(rng.choice("ACGT") for _ in range(n))
    gb = sc._record_to_gb_text(SeqRecord(
        Seq(s), id=name, name=name, description=name,
        annotations={"molecule_type": "DNA", "topology": "circular"}))
    return {"id": name, "name": name, "size": n, "n_feats": 0, "gb_text": gb}


def _disk(path):
    import json
    raw = json.loads(path.read_text())
    return raw["entries"] if isinstance(raw, dict) else raw


def _colls():
    return {c["name"]: sorted(p["name"] for p in c.get("plasmids") or [])
            for c in _disk(sc._state._COLLECTIONS_FILE)}


def _pcolls():
    return {c["name"]: sorted(p["name"] for p in c.get("primers") or [])
            for c in _disk(sc._state._PRIMER_COLLECTIONS_FILE)}


class _Rig:
    _unsaved = False

    def call_from_thread(self, fn, *a, **k):
        return fn(*a, **k)


def _agent(name, payload, app=None):
    return sc._state._AGENT_HANDLERS[name][0](app or _Rig(), payload)


def _settings_backup(pred):
    import json
    rows = _agent("list-backups", {"label": "settings"})["backups"]
    for r in rows:
        raw = json.loads(open(r["source_path"]).read())
        ents = raw["entries"] if isinstance(raw, dict) else raw
        if pred({e["key"]: e["value"] for e in ents}):
            return r["source_path"]
    raise AssertionError("no matching settings backup")


class TestSettingsRestoreNeverLandsOnTheWrongContainer:

    def test_a_dirty_library_is_rescued_and_the_restored_collection_kept(
            self, _sync_saves):
        # collections.json restored from before A existed marks the library
        # dirty (it holds A's only copy); restoring settings to B then left
        # the marker, the re-stage refused, and the next save wrote A's
        # plasmids into B over B's own.
        import json
        sc._save_collections([{"name": "B", "plasmids": [_ent("b1"), _ent("b2")]}])
        sc._activate_collection("B")
        sc._save_collections(sc._load_collections()
                             + [{"name": "A", "plasmids": [_ent("a1")]}])
        sc._activate_collection("A")
        rows = _agent("list-backups", {"label": "collections"})["backups"]
        src = next(r["source_path"] for r in rows
                   if [c["name"] for c in (lambda raw: raw["entries"]
                       if isinstance(raw, dict) else raw)(
                           json.loads(open(r["source_path"]).read()))] == ["B"])
        assert _agent("restore-backup", {"label": "collections",
                                         "source_path": src}).get("ok")
        src2 = _settings_backup(lambda kv: kv.get("active_collection") == "B")
        assert _agent("restore-backup", {"label": "settings",
                                         "source_path": src2}).get("ok")
        sc._save_library(sc._load_library() + [_ent("new1")])
        c = _colls()
        assert {"b1", "b2", "new1"} <= set(c["B"]), c
        # ...and A's plasmid is not lost either: it was kept, not overwritten.
        assert any("a1" in v for v in c.values()), c

    def test_a_backup_pointing_at_a_deleted_collection_is_refused(
            self, _sync_saves):
        sc._save_collections([{"name": "A", "plasmids": [_ent("a1")]},
                              {"name": "Gone", "plasmids": []}])
        sc._activate_collection("Gone")
        sc._activate_collection("A")
        src = _settings_backup(lambda kv: kv.get("active_collection") == "Gone")
        _agent("delete-collection", {"name": "Gone"})
        out = _agent("restore-backup", {"label": "settings", "source_path": src})
        assert isinstance(out, tuple) and out[1] == 409, out
        assert "Gone" in out[0]["error"]
        assert sc._get_active_collection_name() == "A"

    def test_no_save_can_land_between_the_restore_steps(self, _sync_saves,
                                                         monkeypatch):
        import splicecraft_agent as agent_mod
        sc._save_collections([{"name": "A", "plasmids": [_ent("a1")]},
                              {"name": "B", "plasmids": [_ent("b1")]}])
        sc._activate_collection("B")
        sc._activate_collection("A")
        src = _settings_backup(lambda kv: kv.get("active_collection") == "B")
        got = {}
        real = agent_mod._restore_from_backup

        def restore_and_probe(*a, **k):
            n = real(*a, **k)
            box = {}
            t = threading.Thread(target=lambda: box.update(
                free=sc._state._cache_lock.acquire(timeout=0.2)))
            t.start()
            t.join(5)
            if box.get("free"):
                sc._state._cache_lock.release()
            got["another_thread_could_save"] = box.get("free")
            return n
        monkeypatch.setattr(agent_mod, "_restore_from_backup", restore_and_probe)
        assert _agent("restore-backup", {"label": "settings",
                                         "source_path": src}).get("ok")
        assert got["another_thread_could_save"] is False


class TestAnEditIsNotHeldByAnOlderCopy:
    def test_an_edit_made_with_no_collection_active_survives_a_switch(
            self, _sync_saves):
        # The pointer names a collection another session deleted. An edit to
        # a1 — whose id ALSO sits, unedited, in C — was called "held" by id
        # and lost at the next switch.
        sc._save_collections([{"name": "A", "plasmids": [_ent("a1")]},
                              {"name": "B", "plasmids": [_ent("b1")]},
                              {"name": "C", "plasmids": [_ent("a1")]}])
        sc._activate_collection("A")
        sc._save_collections([c for c in sc._load_collections()
                              if c["name"] != "A"])   # another session
        edited = _ent("a1", seed=999, n=150)
        sc._save_library([edited])
        sc._activate_collection("B")
        sizes = sorted(p["size"] for c in _disk(sc._state._COLLECTIONS_FILE)
                       for p in c.get("plasmids") or []
                       if p.get("name") == "a1")
        assert 150 in sizes, sizes

    def test_a_second_switch_does_not_recover_it_again(self, _sync_saves):
        sc._save_collections([{"name": "A", "plasmids": [_ent("a1")]},
                              {"name": "B", "plasmids": [_ent("b1")]},
                              {"name": "C", "plasmids": [_ent("a1")]}])
        sc._activate_collection("A")
        sc._save_collections([c for c in sc._load_collections()
                              if c["name"] != "A"])
        sc._save_library([_ent("a1", seed=999, n=150)])
        sc._activate_collection("B")
        n1 = len(_colls().get("Recovered plasmids", []))
        # back to a state with no collection active and the same library
        sc._state._library_cache = [_ent("a1", seed=999, n=150)]
        sc._set_active_collection_name("Deleted-too")
        sc._activate_collection("C")
        assert len(_colls().get("Recovered plasmids", [])) == n1


class TestEveryPrimerSwitchKeepsFreeStandingPrimers:
    def _seed(self):
        sc._save_collections([{"name": "A", "plasmids": [_ent("a1")]}])
        sc._activate_collection("A")
        sc._save_primer_collections([
            {"name": "PA", "primers": [{"name": "pa1", "sequence": "ACGTACGTAAC"}]},
            {"name": "PB", "primers": [{"name": "pb1", "sequence": "TTGGCCAATTG"}]}])
        _agent("set-active-primer-collection", {"name": "PA"})

    def _everywhere(self):
        held = {p for v in _pcolls().values() for p in v}
        return held | {p["name"] for p in _disk(sc._state._PRIMERS_FILE)}

    def test_a_click_on_another_collection_in_the_primers_screen(
            self, _sync_saves):
        from types import SimpleNamespace
        self._seed()

        async def run():
            app = sc.PlasmidApp()
            async with app.run_test() as pilot:
                await pilot.pause()
                _agent("set-active-primer-collection", {"name": ""}, app=app)
                _agent("create-primer", {"name": "orph",
                                         "sequence": "GGGTTTCCCAAA"}, app=app)
                app._apply_record(sc._make_demo_record())
                await pilot.pause()
                app.action_open_primer_design()
                await pilot.pause()
                await pilot.pause()
                scr = app.screen
                if type(scr).__name__ != "PrimerDesignScreen":
                    scr.dismiss(None)
                    await pilot.pause()
                    scr = app.screen
                scr._lib_view_mode = "collections"
                scr._pd_lib_row_selected(
                    SimpleNamespace(row_key=SimpleNamespace(value="PB")))
                await pilot.pause()
        asyncio.run(run())
        assert sc._get_active_primer_collection_name() == "PB"
        assert "orph" in self._everywhere()

    def test_an_agent_switch_away_from_a_deleted_collection(self, _sync_saves):
        self._seed()
        sc._save_primer_collections(sc._load_primer_collections()
                                    + [{"name": "PGone", "primers": []}])
        _agent("set-active-primer-collection", {"name": "PGone"})
        sc._save_primer_collections([c for c in sc._load_primer_collections()
                                     if c["name"] != "PGone"])  # elsewhere
        _agent("create-primer", {"name": "work1", "sequence": "GGGTTTCCCAAA"})
        _agent("set-active-primer-collection", {"name": "PB"})
        assert "work1" in self._everywhere()


class TestASpillRestoresEverythingItHolds:
    """A lost-entries spill is merged back, not written over. Collections
    whose items carry no id (primers, parts, enzyme names) were not
    recognised as collections, and entry vectors were matched by NAME."""

    @staticmethod
    def _spill(path):
        import splicecraft_backup as B
        return [r for r in B._list_recoverable_backups(path)
                if "lost_entries" in str(r.get("source_path"))]

    def test_a_primer_collection_gets_its_primers_back(self, _sync_saves):
        from pathlib import Path
        import splicecraft_backup as B
        P = [{"name": f"p{i}", "sequence": s} for i, s in enumerate(
            ["ACGTACGTAAC", "TTGGCCAATTG", "GGGCCCAAATT", "CATCATCATCA"])]
        full = ([{"name": "PA", "primers": P}]
                + [{"name": f"X{i}", "primers": []} for i in range(5)])
        sc._save_primer_collections(full)
        sc._save_primer_collections([{"name": "PA", "primers": P[:1]}])
        rows = self._spill(sc._state._PRIMER_COLLECTIONS_FILE)
        B._restore_from_backup(sc._state._PRIMER_COLLECTIONS_FILE,
                               Path(rows[0]["source_path"]), "Primer collections")
        got = _pcolls()
        assert got["PA"] == ["p0", "p1", "p2", "p3"], got
        assert {f"X{i}" for i in range(5)} <= set(got)

    def test_an_entry_vector_is_its_grammar_and_role_not_its_name(
            self, _sync_saves):
        from pathlib import Path
        import splicecraft_backup as B
        ev = [{"grammar_id": "gb_l0", "role": r, "name": f"FFE {i}",
               "size": 100, "source": "x", "gb_text": ""}
              for i, r in enumerate(["", "r1", "r2", "r3", "r4"])]
        ev.append({"grammar_id": "my_gb", "role": "r1", "name": "FFE 1",
                   "size": 100, "source": "x", "gb_text": ""})
        sc._save_entry_vectors(ev)
        sc._save_entry_vectors([ev[1]])
        rows = self._spill(sc._state._ENTRY_VECTORS_FILE)
        B._restore_from_backup(sc._state._ENTRY_VECTORS_FILE,
                               Path(rows[0]["source_path"]), "Entry vectors")
        got = {(e["grammar_id"], e["role"])
               for e in _disk(sc._state._ENTRY_VECTORS_FILE)}
        assert got == {(e["grammar_id"], e["role"]) for e in ev}

    def test_an_enzyme_collection_gets_its_enzymes_back(self):
        import splicecraft_backup as B
        merged, n = B._merge_spilled_entries(
            [{"name": "Mine", "enzymes": ["EcoRI"]}],
            [{"name": "Mine", "enzymes": ["EcoRI", "BsaI", "BsmBI"]}])
        assert merged[0]["enzymes"] == ["EcoRI", "BsaI", "BsmBI"] and n == 2

    @pytest.mark.parametrize("bad", [[1, 2], {"a": 1}])
    def test_a_hand_edited_name_does_not_crash_the_merge(self, bad):
        import splicecraft_backup as B
        cur = [{"name": "M13F", "sequence": "GTAAAACGACGGCCAGT"}]
        merged, n = B._merge_spilled_entries(cur, [{"name": bad,
                                                    "sequence": "ACGT"}])
        assert n == 1
        merged, n = B._merge_spilled_entries([{"id": ["x"], "name": "a"}],
                                             [{"id": "y", "name": "b"}])
        assert n == 1

    def test_control_a_same_name_primer_variant_still_comes_back(self):
        import splicecraft_backup as B
        cur = [{"name": "M13F", "sequence": "GTAAAACGACGGCCAGT"}]
        merged, n = B._merge_spilled_entries(
            cur, [{"name": "M13F", "sequence": "TGTAAAACGACGGCCAGT"},
                  {"name": "M13F", "sequence": "gtaaaacgac ggccagt"}])
        assert n == 1 and len(merged) == 2


class TestSmallDataSafetyGaps:
    def test_a_cancelled_import_says_saved_only_once_it_is(self, _sync_saves):
        # The import's name was already in the library, so a name-collision
        # prompt opened — and "Saved pFoo to the library" was shown at once,
        # before the user had chosen anything.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        sc._save_collections([{"name": "A", "plasmids": [_ent("pFoo")]}])
        sc._activate_collection("A")
        got: dict = {}

        async def run():
            app = sc.PlasmidApp()
            notes: list = []
            real = app.notify
            app.notify = lambda msg, *a, **k: (notes.append(str(msg)),
                                               real(msg, *a, **k))[1]
            async with app.run_test() as pilot:
                await pilot.pause()
                app._apply_record(sc._make_demo_record())
                await pilot.pause()
                app._unsaved = True
                rec = SeqRecord(Seq("GATC" * 70), id="NEWID1", name="pFoo",
                                description="imported",
                                annotations={"molecule_type": "DNA",
                                             "topology": "circular"})
                app._import_and_persist(rec)
                await pilot.pause()
                app.screen.dismiss(None)     # keep my edits, not the import
                await pilot.pause()
                await pilot.pause()
                got["prompt"] = type(app.screen).__name__
                got["notes"] = list(notes)
        asyncio.run(run())
        assert got["prompt"] == "NameCollisionModal", got
        assert not any(n.startswith("Saved pFoo") for n in got["notes"]), \
            got["notes"]

    def test_a_token_file_that_is_not_text_does_not_stop_shutdown(
            self, tmp_path, monkeypatch):
        from tests.test_agent_api import MockApp, _free_port
        tok = tmp_path / "agent_token"
        monkeypatch.setattr(sc._state, "_AGENT_TOKEN_FILE", tok)
        srv = sc._start_agent_api(MockApp(), port=_free_port())
        assert srv is not None
        tok.write_bytes(b"\xff\xfe\x00garbage\x80")
        sc._stop_agent_api(srv)                  # must not raise


# ══════════════════════════════════════════════════════════════════════════
# Reads and alignment — every expectation is the planted truth
# ══════════════════════════════════════════════════════════════════════════

def _rnd(n, seed):
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _verify(ref, reads, **kw):
    return sc._state._AGENT_HANDLERS["verify-against-reads"][0](
        None, {"reference": ref, "reads": reads, **kw})


def _hetero(**body):
    return sc._state._AGENT_HANDLERS["analyse-read-heterogeneity"][0](None,
                                                                      body)


class TestEachReadCountsOnce:
    """One indel crossing bp 0 arrives as two pieces at the two ends of the
    rows. Normalised onto one key, the read was counted twice."""

    def test_a_50_50_culture_is_mixed_not_clonal(self):
        core = "C" + _rnd(1492, 5) + "C"
        ref, clone = "AAAA" + core + "AAAA", "AAA" + core + "AAA"
        h = _hetero(reference=ref, reads=[clone] * 5 + [ref] * 5,
                    circular=True, platform="illumina")
        assert h["verdict"] == "mixed", h
        [pos] = [p for p in h["positions"] if p["type"] == "deletion"]
        assert (pos["alt_reads"], pos["depth"]) == (5, 10)

    def test_one_read_never_confirms_its_own_variant(self):
        core = "C" + _rnd(1492, 5) + "C"
        ref = "AAAA" + core + "AAAA"
        n = len(ref)
        res = {"aligned_q": "-" + ref[1:n - 1] + "-", "aligned_t": ref,
               "identity_pct": 98.0}
        summ = sc._multi_read_summary(
            [{"result": res, "axis": "target", "query_label": "r1"}], n,
            circular=True)
        assert summ["verdict"] != "confirmed"
        assert [(v["type"], v["length"], v["support"])
                for v in summ["variants"]] == [("deletion", 2, 1)]

    @pytest.mark.parametrize("clone_of", ["ins", "del"])
    def test_an_indel_at_the_origin_is_one_clonal_event(self, clone_of):
        ref = _rnd(3000, 7)
        clone = ("GATTACAGGT" + ref) if clone_of == "ins" else ref[5:-5]
        m = len(clone)
        reads = [clone[(i * 149) % m:] + clone[:(i * 149) % m]
                 for i in range(20)]
        h = _hetero(reference=ref, reads=reads, circular=True,
                    platform="illumina")
        kind = "insertion" if clone_of == "ins" else "deletion"
        assert [(p["type"], p["alt_reads"], p["depth"])
                for p in h["positions"]] == [(kind, 20, 20)], h["positions"]


class TestAReadIsJudgedOnWhatItRead:
    def test_a_deletion_across_the_origin_fails_verification(self):
        ref = _rnd(3000, 7)
        n = len(ref)
        for L in (40, 120):
            h = L // 2
            clone = ref[h:n - (L - h)]
            for o in (0, 1, 3, 679, len(clone) - 1):
                read = clone[o:] + clone[:o]
                rd = _verify(ref, [read], circular=True)["reads"][0]
                assert rd["passes"] is False, (L, o, rd)

    def test_a_soft_clipped_end_is_not_a_deletion(self):
        ref = _rnd(3000, 7)
        r = _verify(ref, [ref[-1:] + ref[:-1]], circular=True, mode="local",
                    read_quality=[[40] * 3000])
        assert r["read_quality"][0]["verdict"] == "clean", r["read_quality"]
        assert not [v for v in r["consensus"]["variants"]
                    if v["type"] == "deletion"]

    def test_a_linear_reference_is_not_read_as_a_circle(self):
        ref = _rnd(3000, 7)
        r = _verify(ref, [ref[100:3000]], circular=False,
                    read_quality=[[40] * 2900])
        assert r["read_quality"][0]["verdict"] == "clean", r["read_quality"]
        assert r["consensus"]["variants"] == []

    def test_a_partial_read_carrying_a_deletion_is_not_full_length(self):
        ref = _rnd(3000, 7)
        read = ref[100:1500] + ref[1600:2980]
        r = _verify(ref, [read], circular=True,
                    read_quality=[[40] * len(read)])
        dels = [v for v in r["consensus"]["variants"]
                if v["type"] == "deletion"]
        assert [d["length"] for d in dels] == [100], dels
        assert r["read_quality"][0]["n_indel_confident"] == 1

    def test_a_read_too_short_to_say_anything_does_not_pass(self):
        ref = _rnd(3000, 7)
        r = _verify(ref, ["GATTACA"], circular=True)
        assert r["reads"][0]["passes"] is False and r["verdict"] != "match"
        assert r["reads"][0]["too_short"] is True

    def test_a_masked_read_end_is_no_evidence(self):
        ref = _rnd(3000, 7)
        read = "N" * 30 + ref[1030:1740] + "N" * 60
        r = _verify(ref, [read], circular=True)
        assert r["reads"][0]["passes"] is True, r["reads"][0]
        assert r["consensus"]["variants"] == []

    def test_an_infinite_stored_frame_shift_does_not_crash(self):
        import splicecraft_biology as bio
        from splicecraft_util import _phred_in_alignment_frame
        assert bio._read_end_columns("AC--GT", float("inf")) == (0, 5)
        assert _phred_in_alignment_frame(
            [30] * 4, {"query_frame_shift": float("inf")}) == [30] * 4


class TestInversionsAtTheReadEnds:
    """The probe spared the unread columns beside each read end and grew
    candidates into them: a flip at a read's start was reported 35 bp into
    plasmid the read never covered, and an origin-straddling join stitched in
    unread plasmid. A too-short sliver there ended the search before a real
    inversion behind it."""

    @staticmethod
    def _flip(base, s, L):
        n = len(base)
        idx = [(s + k) % n for k in range(L)]
        inv = sc._rc("".join(base[i] for i in idx))
        cons = list(base)
        for k, i in enumerate(idx):
            cons[i] = inv[k]
        return "".join(cons)

    @staticmethod
    def _found(read, base, mode="global"):
        res = sc._pick_best_rotation(read, base, is_circular=True, mode=mode)
        segs = res.get("inverted_segments") or []
        got: set = set()
        for s in segs:
            got |= set(range(s["t_start"], s["t_end"]))
        return got, [(s["t_start"], s["t_end"]) for s in segs]

    def test_a_flip_at_a_read_start_is_reported_where_it_is(self):
        base = _rnd(4000, 101)
        cons = self._flip(base, 1500, 30)
        got, segs = self._found(cons[1500:2300], base)
        assert segs == [(1500, 1530)], segs

    def test_a_flip_at_the_start_of_a_read_at_the_origin_is_found(self):
        base = _rnd(4000, 101)
        got, segs = self._found(self._flip(base, 0, 60)[:800], base)
        truth = set(range(60))
        assert len(got & truth) >= 50 and not (got - truth), segs

    @pytest.mark.parametrize("seed", [101, 102])
    def test_nothing_is_reported_in_plasmid_no_read_covered(self, seed):
        base = _rnd(4000, seed)
        got, segs = self._found(self._flip(base, 3940, 60)[3200:], base)
        truth = set(range(3940, 4000))
        # A flanking base that happens to fit the reverse complement makes
        # the boundary ambiguous by a base or two; the old probe reported a
        # whole second half — 30-53 bp at bp 0 — that no read covered.
        assert len(got & truth) >= 50 and len(got - truth) <= 3, segs


# ══════════════════════════════════════════════════════════════════════════
# Agent API
# ══════════════════════════════════════════════════════════════════════════

def test_a_handler_whose_file_changed_is_not_read_as_another(tmp_path):
    # An update under a running server: `inspect.getsource` read the NEW
    # file at the OLD line numbers and credited the handler with another
    # function's keys.
    import importlib.util
    import splicecraft_agent as ag
    src = tmp_path / "hmod_stale.py"
    src.write_text("def other(app, payload):\n    return 1\n\n\n"
                   "def handler(app, payload):\n"
                   "    return payload.get('name')\n")
    spec = importlib.util.spec_from_file_location("hmod_stale", src)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    src.write_text("def a(app, payload):\n    return 1\n\n\n"
                   "def other(app, payload):\n"
                   "    return payload.get('collection')\n\n\n"
                   "def handler(app, payload):\n"
                   "    return payload.get('name')\n")
    keys = ag._agent_payload_keys(mod.handler)
    assert keys is None or "collection" not in keys, keys


def test_the_server_works_out_every_handlers_keys_when_it_starts(tmp_path,
                                                                  monkeypatch):
    import splicecraft_agent as ag
    monkeypatch.setattr(ag, "_PAYLOAD_KEYS_CACHE", {})
    srv, port, token = _serve_agent(tmp_path, monkeypatch)
    handlers = {e[0] for e in sc._state._AGENT_HANDLERS.values()}
    try:
        deadline = time.monotonic() + 30
        while (not handlers <= set(ag._PAYLOAD_KEYS_CACHE)
               and time.monotonic() < deadline):
            time.sleep(0.05)
    finally:
        sc._stop_agent_api(srv)
    assert handlers <= set(ag._PAYLOAD_KEYS_CACHE)


class TestAgentStatusesSayWhatHappened:
    @pytest.mark.parametrize("v", ["false", "no", 0, "0"])
    def test_a_false_flag_spelled_as_text_is_not_refused(self, v):
        # `carry_annotations: false` was accepted and "false"/"no"/0 refused
        # by the central routing guard (v1.2.71 accepted them all).
        body, status = sc._agent_invoke(_Rig(), "domesticate-part",
                                        {"carry_annotations": v})
        assert "carry_annotations" not in str(body.get("error", "")), body

    def test_a_true_flag_or_a_name_is_still_checked(self):
        body, status = sc._agent_invoke(_Rig(), "domesticate-part",
                                        {"carry_annotations": "true"})
        assert status == 400 and "carry_annotations" in str(body), body
        body, status = sc._agent_invoke(_Rig(), "domesticate-part",
                                        {"parts_bin": "0"})
        assert status == 400 and "parts_bin" in str(body), body

    def test_a_custom_enzyme_named_with_a_comma_is_one_enzyme(self,
                                                              _sync_saves):
        # `create-custom-enzyme` accepts a comma in a name; splitting the
        # string first turned "Cus,1" into an unknown "Cus".
        made = _agent("create-custom-enzyme", {
            "name": "Cus,1", "site": "GCTAGCTTAG", "fwd_cut": 5, "rev_cut": 5})
        assert not isinstance(made, tuple), made
        seq = "A" * 10 + "GCTAGCTTAG" + "A" * 40
        out = _agent("digest", {"sequence": seq, "enzymes": "Cus,1",
                                "circular": False})
        assert not isinstance(out, tuple), out
        assert not out.get("unknown_enzymes")
        assert len(out.get("fragments") or []) == 2, out

    def test_deleting_only_ambiguous_primers_is_a_conflict(self, _sync_saves):
        sc._save_primers([{"name": "M13F", "sequence": "GTAAAACGACGGCCAGT"},
                          {"name": "M13F", "sequence": "TGTAAAACGACGGCCAGT"}])
        out = _agent("delete-primers", {"names": ["M13F"]})
        assert isinstance(out, tuple) and out[1] == 409, out
        assert len(sc._load_primers()) == 2

    def test_a_partial_primer_delete_says_so(self, _sync_saves):
        sc._save_primers([{"name": "a1", "sequence": "ACGTACGTACGTA"},
                          {"name": "b1", "sequence": "TTGGCCAATTGGC"}])
        out = _agent("delete-primers", {"names": ["a1", "nope"]})
        assert not isinstance(out, tuple), out
        assert out["removed"] == 1 and out["partial"] is True
        assert any("nope" in w for w in out["warnings"])

    def test_copying_nothing_is_not_ok(self, _sync_saves):
        sc._save_collections([{"name": "A", "plasmids": [_ent("a1")]},
                              {"name": "B", "plasmids": []}])
        sc._activate_collection("A")
        out = _agent("copy-plasmids", {"names": ["nope"], "to": "B"})
        assert isinstance(out, tuple) and out[1] == 404, out
        out = _agent("copy-plasmids", {"names": ["a1", "nope"], "to": "B"})
        assert not isinstance(out, tuple) and out["partial"] is True, out

    def test_a_save_that_failed_is_not_a_200(self, monkeypatch):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        class _App(_Rig):
            _current_record = SeqRecord(Seq("ACGT" * 30), id="p", name="p")
            _source_path = None
            _last_save_error = "library save failed"
            _last_save_deferred = False

            def call_from_thread(self, fn, *a, **k):
                return False                     # `_do_save` failed

            def _do_save(self):
                return False
        out = sc._state._AGENT_HANDLERS["save"][0](_App(), {"create": True})
        assert isinstance(out, tuple) and out[1] in (409, 500), out
        assert out[0]["ok"] is False

    def test_a_routing_key_the_clone_does_not_read_is_refused(self):
        import splicecraft_agent as ag
        fn = sc._state._AGENT_HANDLERS["traditional-clone"][0]
        assert ag._agent_unread_routing_keys(
            fn, {"source_collection": "X"}) == ["source_collection"]
        assert ag._agent_unread_routing_keys(
            fn, {"vector_collection": "X"}) == []


class TestAgentHttpFraming:
    def test_an_unhashable_idempotent_body_is_refused(self, tmp_path,
                                                      monkeypatch):
        import http.client
        srv, port, token = _serve_agent(tmp_path, monkeypatch)
        deep = "[" * 3000 + "]" * 3000
        body = ('{"name": "pdeep", "sequence": "ACGTACGTACGTAC", "x": '
                + deep + "}").encode()
        try:
            conn = http.client.HTTPConnection("127.0.0.1", port, timeout=20)
            conn.request("POST", "/create-primer", body=body, headers={
                "Authorization": f"Bearer {token}",
                "Content-Type": "application/json",
                "X-Idempotency-Key": "deep-1"})
            r = conn.getresponse()
            status = r.status
            r.read()
            conn.close()
        finally:
            sc._stop_agent_api(srv)
        assert status == 400
        assert not any(p.get("name") == "pdeep" for p in sc._load_primers())

    def test_an_http_1_0_body_without_a_length_is_refused(self, tmp_path,
                                                          monkeypatch):
        import socket as _socket
        srv, port, token = _serve_agent(tmp_path, monkeypatch)
        try:
            s = _socket.create_connection(("127.0.0.1", port), timeout=10)
            s.sendall(("POST /create-primer HTTP/1.0\r\n"
                       f"Authorization: Bearer {token}\r\n"
                       "Content-Type: application/json\r\n\r\n"
                       '{"name": "p10", "sequence": "ACGTACGTACGT"}').encode())
            s.shutdown(_socket.SHUT_WR)
            data = b""
            while chunk := s.recv(65536):
                data += chunk
            s.close()
        finally:
            sc._stop_agent_api(srv)
        assert data.split(b" ", 2)[1] == b"411", data[:200]


# ══════════════════════════════════════════════════════════════════════════
# File formats
# ══════════════════════════════════════════════════════════════════════════

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


def _dna_bytes(rec):
    import splicecraft_fileio as _fio
    return _fio._write_commercialsaas_dna_bytes(rec)


def _load_dna(tmp_path, data, name="t.dna"):
    p = tmp_path / name
    p.write_bytes(data)
    return sc.load_genbank(str(p))


def _v1271_dna_cases():
    import json
    import pathlib
    raw = json.loads((pathlib.Path(__file__).parent
                      / "dna_v1271_entries.json").read_text())
    return raw["cases"]


class TestAnUneditedOldImportStillSplices:
    """An entry v1.2.71 imported, never edited, must still export by splicing
    its original `.dna` — rebuilt from scratch it loses the reads, primers
    and enzymes only the original carries. These entries were written by the
    real v1.2.71 code (`tests/dna_v1271_entries.json`)."""

    @pytest.mark.parametrize("case", ["same_span", "u2028", "cr", "nel"])
    def test_it_matches_its_original(self, case):
        import base64
        import splicecraft_fileio as fio
        c = _v1271_dna_cases()[case]
        entry = {"id": c["id"], "name": c["name"], "gb_text": c["gb_text"]}
        assert fio._dna_sidecar_matches_entry(
            base64.b64decode(c["dna_b64"]), entry), case


def _dna_with_notes(seq, notes_xml, features_xml=None):
    import splicecraft_fileio as fio
    data = (fio._build_commercialsaas_cookie_packet()
            + fio._build_commercialsaas_dna_packet(seq, circular=True))
    if features_xml:
        data += fio._build_commercialsaas_packet(0x0A, features_xml.encode())
    return data + fio._build_commercialsaas_packet(0x06, notes_xml.encode())


class TestDnaIdentityFields:
    def test_a_renamed_entry_exports_under_its_new_name(self, tmp_path,
                                                        _sync_saves):
        # The splice export kept the original Notes packet, whose map label
        # then won over the file name on reopen: the old name came back.
        import splicecraft_fileio as fio
        data = _dna_with_notes(
            _rnd(400, 11),
            "<Notes><Type>Synthetic</Type><CustomMapLabel>pOld</CustomMapLabel>"
            "<UseCustomMapLabel>1</UseCustomMapLabel></Notes>")
        rec = _load_dna(tmp_path, data, "pOld.dna")
        entry = sc._record_to_library_entry(rec, tmp_path / "pOld.dna")
        fio._save_dna_original(entry["id"], data)
        assert fio._dna_sidecar_matches_entry(data, entry)
        entry["name"] = "pNew"
        out = tmp_path / "exported.dna"
        fio._export_commercialsaas_dna(entry, out)
        back = sc.load_genbank(str(out))
        assert getattr(back, "_tui_display_name", None) == "pNew"

    def test_an_unedited_splice_stays_byte_identical(self, tmp_path,
                                                     _sync_saves):
        import splicecraft_fileio as fio
        data = _dna_with_notes(
            _rnd(400, 12),
            "<Notes><Type>Synthetic</Type><CustomMapLabel>pSame</CustomMapLabel>"
            "<UseCustomMapLabel>1</UseCustomMapLabel></Notes>")
        rec = _load_dna(tmp_path, data, "pSame.dna")
        entry = sc._record_to_library_entry(rec, tmp_path / "pSame.dna")
        fio._save_dna_original(entry["id"], data)
        out = tmp_path / "again.dna"
        fio._export_commercialsaas_dna(entry, out)
        assert out.read_bytes() == data

    def test_a_label_is_plain_text(self, tmp_path):
        from xml.sax.saxutils import escape
        data = _dna_with_notes(
            _rnd(300, 13),
            "<Notes><CustomMapLabel>" + escape("<Cas9> donor v2")
            + "</CustomMapLabel><UseCustomMapLabel>1</UseCustomMapLabel>"
              "</Notes>")
        rec = _load_dna(tmp_path, data, "x.dna")
        assert rec._tui_display_name == "<Cas9> donor v2"

    def test_a_label_that_is_nothing_keeps_the_file_name(self, tmp_path):
        data = _dna_with_notes(
            _rnd(300, 14),
            "<Notes><CustomMapLabel>‎</CustomMapLabel>"
            "<UseCustomMapLabel>1</UseCustomMapLabel></Notes>")
        rec = _load_dna(tmp_path, data, "my plasmid.dna")
        assert rec._tui_display_name == "my plasmid"

    def test_a_plain_description_keeps_its_spacing(self, tmp_path):
        data = _dna_with_notes(
            _rnd(300, 15), "<Notes><Comments>a  b   c</Comments></Notes>")
        rec = _load_dna(tmp_path, data, "d.dna")
        assert rec.description == "a  b   c"

    def test_a_distinct_reverse_colour_survives(self, tmp_path):
        rec = _rec_with([_feat(10, 30, ApEinfo_fwdcolor="#ff0000",
                               ApEinfo_revcolor="#0000ff")])
        back = _load_dna(tmp_path, _dna_bytes(rec))
        q = back.features[0].qualifiers
        assert q.get("ApEinfo_revcolor") == ["#0000ff"], q


class TestDnaWritesOnlyWhatItCanReadBack:
    def test_an_oversized_history_is_refused_not_dropped(self, monkeypatch):
        import splicecraft_fileio as fio
        monkeypatch.setattr(fio, "_COMMERCIALSAAS_HISTORY_MAX_XML", 1000)
        with pytest.raises(ValueError, match="construction history"):
            fio._pack_commercialsaas_history_payload("<H>" + "x" * 2000 + "</H>")

    def test_the_agent_answers_422_for_a_refusal(self, monkeypatch, tmp_path):
        import splicecraft_agent as ag
        monkeypatch.setattr(ag, "_export_commercialsaas_dna",
                            lambda entry, path: (_ for _ in ()).throw(
                                ValueError("the .dna features block would be "
                                           "40 MB")))
        monkeypatch.setattr(ag, "_find_library_entry_by_id",
                            lambda *a, **k: {"id": "p", "name": "p",
                                             "gb_text": "x"}, raising=False)
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        class _App(_Rig):
            _current_record = SeqRecord(Seq("ACGT" * 30), id="p", name="p")
        out = sc._state._AGENT_HANDLERS["export-commercialsaas"][0](
            _App(), {"path": str(tmp_path / "o.dna")})
        assert isinstance(out, tuple) and out[1] == 422, out


class TestTheNameMarkerEndsWhereTheNameDoes:
    """The wrap rule could not always tell our own wrapped continuation from a
    COMMENT line another tool appended; the name absorbed that line, and the
    next save deleted it from the COMMENT."""

    @pytest.mark.parametrize("name,appended", [
        ("pGreenII 0229 35S:GFP-NLS construct v2 final",
         "Sequence verified by vendor"),
        ("pUC19 mod", "https://example.org/s/seq-" + "a" * 60),
    ])
    def test_another_tools_line_is_not_part_of_the_name(self, name, appended):
        rec = _rec_with([])
        rec._tui_display_name = name
        gb = sc._record_to_gb_text(rec)
        # another tool appends its own COMMENT line after ours
        lines = gb.split("\n")
        i = max(k for k, ln in enumerate(lines) if ln.startswith("            ")
                or ln.startswith("COMMENT"))
        lines.insert(i + 1, "            " + appended)
        back = sc._gb_text_to_record("\n".join(lines))
        assert getattr(back, "_tui_display_name", None) == name


# ══════════════════════════════════════════════════════════════════════════
# Biology
# ══════════════════════════════════════════════════════════════════════════

def test_a_stale_selenocysteine_declaration_does_not_rename_a_residue():
    # `/transl_except` is written in base coordinates and nothing moves it on
    # an upstream edit, so it named another codon — and residue 15 (a Gln)
    # was reported as the selenocysteine.
    import random
    import types
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rng = random.Random(4)
    sense = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
             if a + b + c not in ("TAA", "TAG", "TGA")]
    cds = ("ATG" + "".join(rng.choice(sense) for _ in range(15)) + "TGA"
           + "".join(rng.choice(sense) for _ in range(40)) + "TAA")
    pre = "".join(rng.choice("ACGT") for _ in range(300))
    post = "".join(rng.choice("ACGT") for _ in range(400))
    seq = pre + cds + post
    s0, e0 = len(pre), len(pre) + len(cds)
    rec = SeqRecord(Seq(seq), id="sel", name="sel",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"})
    rec.features = [SeqFeature(
        FeatureLocation(s0, e0, strand=1), type="CDS",
        qualifiers={"label": ["selP"],
                    "transl_except": [f"(pos:{s0 + 49}..{s0 + 51},aa:Sec)"]})]
    ins, at = "GAATTC", 100
    new = sc.PlasmidApp._rebuild_record_with_edit(
        types.SimpleNamespace(_current_record=None),
        seq[:at] + ins + seq[at:], "insert", at, at, ins, source_record=rec)
    f = [x for x in new.features if x.type == "CDS"][0]
    aa = sc._translate_cds(str(new.seq), int(f.location.start),
                           int(f.location.end), 1)
    plan = sc._protein_edit_plan(new, f, 15, aa[14])
    assert plan["wt_aa"] == aa[14] and plan["silent"] is True, plan


def test_a_linear_terminator_ends_at_the_sequence_end():
    import splicecraft_regulatory as reg
    hairpin = "GCCGCCGCGG" + "AAAA" + "CCGCGGCGGC"
    seq = "ACGT" * 22 + hairpin + "TTTTTTTT"
    hits = reg._scan_terminators(seq, circular=False, both_strands=False)
    assert hits, "the terminator should be found"
    assert all(h["end"] > h["start"] for h in hits), hits
    assert max(h["end"] for h in hits) == len(seq)


def test_a_linear_staggered_cut_is_not_wrapped():
    import splicecraft_crispr as cr
    body = "".join("ACGT"[(i * 7 + i // 3) % 4] for i in range(40))
    seq = body + "TTTA" + "GACCTAGCATCGATTAGCCATGA"   # PAM + 23 nt at the end
    guides = cr._find_guides(seq, variant="lbcas12a", circular=False)
    tail = [g for g in guides if g["strand"] == 1 and g["start"] == 44]
    assert tail, guides
    g = tail[0]
    assert g["end"] == len(seq)
    assert g["cut_site_bottom"] - g["cut_site"] == g["overhang"]


def test_one_golden_gate_circle_is_not_counted_twice():
    # A backbone piece that ligates either way round: the flip support found
    # the reverse-complement reading of one circle from a second seed.
    import random
    import splicecraft_cloning as clo
    comp = str.maketrans("ACGT", "TGCA")

    def rc(s):
        return s.translate(comp)[::-1]
    rng = random.Random(2)

    def rnd(n):
        while True:
            s = "".join(rng.choice("ACGT") for _ in range(n))
            if "GGTCT" not in s and "AGACC" not in s:
                return s
    A, C = "ACGG", "TTAC"
    B = rc(A)
    vec = (rnd(300) + "GGTCTC" + "T" + A + rnd(200) + B + "T" + "GAGACC"
           + rnd(60) + "GGTCTC" + "T" + C)
    part = "GGTCTC" + "A" + B + rnd(120) + C + "T" + "GAGACC"
    r = clo._simulate_golden_gate([part], vec)
    text = " ".join(r.get("warnings") or [])
    import re as _re
    counts = [int(x) for x in _re.findall(r"(\d+) different circles", text)]
    assert not counts or max(counts) <= 2, text


def test_a_both_ways_terminator_blocks_both_strands():
    # ◀▶ arrived as "no strand", was read forward only, and blocked nothing
    # coming the other way — with a warning that it "carries no strand".
    import splicecraft_seqanalysis as sa
    seq = "ACGT" * 500
    feats = [
        {"start": 1500, "end": 1560, "strand": -1, "type": "promoter",
         "label": "Prev"},
        {"start": 1000, "end": 1050, "strand": 2, "type": "terminator",
         "label": "Tboth"},
        {"start": 300, "end": 900, "strand": -1, "type": "CDS",
         "label": "geneR"},
    ]
    out = sa._map_transcription(seq, feats, circular=False,
                                include_predicted=False)
    assert not any("carries no" in w and "Tboth" in w
                   for w in out.get("warnings") or []), out.get("warnings")
    terms = [t for t in out["terminators"] if t.get("label") == "Tboth"]
    assert sorted(t["strand"] for t in terms) == [-1, 1]


def test_a_p300_gen2_is_not_taken_for_a_p50():
    # opentrons-shared-data: a P300 GEN2's `backCompatNames` is the P300 GEN1
    # only — a P50 protocol passed the pre-run check, then failed at
    # `load_instrument` after the run had started.
    import splicecraft_opentrons as ot
    miss = ot._ot2_pipette_mismatch(
        [{"mount": "left", "pipetteName": "p50_single"}],
        [{"mount": "left", "model": "p300_single_v2.1"}])
    assert miss, "a P50 protocol must not pass on a P300 GEN2"
    ok = ot._ot2_pipette_mismatch(
        [{"mount": "left", "pipetteName": "p300_single"}],
        [{"mount": "left", "model": "p300_single_v2.1"}])
    assert not ok, ok


def test_discovery_probes_a_public_campus_lan(monkeypatch):
    import splicecraft_opentrons as ot
    monkeypatch.setattr(ot, "_ot2_host_ipv4s",
                        lambda: [("128.84.7.20", False)])
    monkeypatch.setattr(ot, "_OT2_OWN_NETS_CACHE", [None, ()], raising=False)
    assert ot._ot2_is_lan_ip("128.84.7.50")
    assert not ot._ot2_is_lan_ip("8.8.8.8")


# ══════════════════════════════════════════════════════════════════════════
# Notebook, plasmid refs, gels
# ══════════════════════════════════════════════════════════════════════════

def test_a_markdown_export_links_attachments_a_viewer_will_open(tmp_path):
    # With embedding off, every attachment target became a `file:` URI, which
    # CommonMark renderers refuse: the exported notebook showed nothing.
    import markdown_it
    img = tmp_path / "gel photo.png"
    img.write_bytes(b"\\x89PNG")
    trace = tmp_path / "read1.ab1"
    trace.write_bytes(b"ABIF")
    entry = {"id": "exp-1", "title": "t", "image_paths": [img.name, trace.name],
             "body_md": f"see ![gel]({img.name}) and [the read]({trace.name})"}
    srcs = {img.name: img.absolute().as_uri(),
            trace.name: trace.absolute().as_uri()}
    md = sc._experiment_markdown_document(entry, image_srcs=srcs)
    html = markdown_it.MarkdownIt().render(md)
    assert "file:" not in md
    assert f'src="{img}"' in html.replace("%20", " "), html
    assert f'href="{trace}"' in html.replace("%20", " "), html


class TestNotebookPlasmidRefs:
    @staticmethod
    def _p(name, pid=None):
        return {"id": pid or name, "name": name, "gb_text": "", "size": 10}

    def test_the_same_id_in_two_collections_is_not_opened_blind(self,
                                                                _sync_saves):
        sc._save_collections([{"name": "Other", "plasmids": []},
                              {"name": "Archive", "plasmids": [self._p("pUC19")]},
                              {"name": "Zeta", "plasmids": [self._p("pUC19")]}])
        sc._activate_collection("Other")
        chosen, matches = sc._resolve_plasmid_ref("pUC19")
        assert chosen is None
        assert {m["collection"] for m in matches} == {"Archive", "Zeta"}

    def test_insert_never_writes_a_tag_that_opens_another_copy(self,
                                                               _sync_saves):
        sc._save_collections([{"name": "Main", "plasmids": [self._p("pUC19")]},
                              {"name": "Archive", "plasmids": [self._p("pUC19")]}])
        sc._activate_collection("Main")
        tok = sc._plasmid_ref_token("Archive", "pUC19")
        assert tok is None or sc._resolve_plasmid_ref(tok)[0] == {
            "collection": "Archive", "id": "pUC19"}, tok

    def test_a_case_variant_finds_the_plasmid_among_many_look_alikes(
            self, _sync_saves):
        sc._save_collections([
            {"name": "Main", "plasmids": [self._p("PUC19", "puc19_main")]},
            {"name": "A-Archive", "plasmids": [
                self._p(f"pUC19 clone {i}", f"c{i}") for i in range(12)]}])
        sc._activate_collection("Main")
        chosen, _m = sc._resolve_plasmid_ref("puc19")
        assert chosen == {"collection": "Main", "id": "puc19_main"}, _m

    def test_a_case_variant_prefers_the_active_collection(self, _sync_saves):
        sc._save_collections([
            {"name": "A", "plasmids": [self._p("pUC19", "a1")]},
            {"name": "Main", "plasmids": [self._p("pUC19", "m1")]}])
        sc._activate_collection("Main")
        assert sc._resolve_plasmid_ref("puc19")[0] == {"collection": "Main",
                                                       "id": "m1"}


class TestGelLanesAndPercentages:
    @staticmethod
    def _gel(length):
        out = _agent("simulate-gel", {"lanes": [{"source": "pcr"}],
                                      "pcr_amplicon": {"length": length}})
        return out if isinstance(out, tuple) else (out, 200)

    @pytest.mark.parametrize("length", [1e400, 10**400, True, "²", "9" * 5000],
                             ids=["1e400", "10**400", "true", "superscript",
                                  "5000-digits"])
    def test_an_impossible_pcr_length_draws_no_band(self, length):
        # `1e400` is what JSON `1e400` parses to: `int()` of it raised
        # OverflowError, which the old guard did not catch (a 500), and
        # `true` drew a 1 bp band.
        out, status = self._gel(length)
        assert status == 200, out
        assert out["lanes"][0]["bands"] == []

    @pytest.mark.parametrize("length", [2450, 2450.0, "2450", " 2450 "],
                             ids=["int", "float", "text", "padded-text"])
    def test_a_real_pcr_length_still_draws_its_band(self, length):
        out, status = self._gel(length)
        assert status == 200, out
        assert [b["bp"] for b in out["lanes"][0]["bands"]] == [2450]

    def test_a_gel_at_an_off_menu_percentage_keeps_it_when_loaded(
            self, _sync_saves):
        # The menu showed a saved 0.65% as its nearest choice, and syncing it
        # there fired the change handler, which wrote 0.7 back — so the next
        # save overwrote the gel's own percentage.
        from textual.app import App
        from textual.widgets import Select
        sc._save_gels([sc._normalise_gel_entry({
            "id": "gel-aa11bb22", "name": "run 3", "agarose_pct": 0.65,
            "lanes": [{"name": "L1", "source": "ladder", "detail": "1 kb"}]},
            fresh=True)])
        got: dict = {}

        class _Host(App):
            def on_mount(self) -> None:
                self.push_screen(sc.SimulatorScreen("ATGC" * 100, [], "p",
                                                    "circular"))

        async def run():
            app = _Host()
            async with app.run_test(size=(170, 48)) as pilot:
                await pilot.pause()
                screen = app.screen
                screen.query_one("#sim-tabs").active = "sim-tab-gel"
                await pilot.pause()
                screen._on_gel_library(None)
                await pilot.pause()
                app.screen.dismiss("gel-aa11bb22")
                for _ in range(6):
                    await pilot.pause()
                sel = screen.query_one("#sim-gel-agarose", Select)
                got["loaded"] = (screen._agarose_pct, sel.value)
                sel.value = "0.7"             # the user now picks 0.7
                for _ in range(4):
                    await pilot.pause()
                got["picked"] = screen._agarose_pct
        asyncio.run(run())
        assert got["loaded"] == (0.65, "0.65"), got
        assert got["picked"] == 0.7, got


# ══════════════════════════════════════════════════════════════════════════
# Markup
# ══════════════════════════════════════════════════════════════════════════

def test_a_quoted_tag_in_a_name_is_shown_not_dropped():
    # Textual's style parser reads `b 'x'` as bold, but its markup tokenizer
    # does not take `[b 'x']` for a tag. Kept as one, it swallowed the `[/]`
    # after it, `[/green]` was left unmatched, and the whole toast fell back
    # to raw markup.
    from textual.content import Content
    out = sc._markup_escape_unstyled("[green]Imported [b 'x'][/] ok[/green]")
    assert out is not None
    assert "Imported [b 'x']" in Content.from_markup(out).plain


def test_many_unmatched_closing_tags_are_not_quadratic():
    # Every unmatched closer searched the whole open-tag stack: 20,000 of
    # them after 20,000 opens took seconds (a pasted log in a toast).
    text = "[b]" * 30000 + "[/x]" * 30000
    t0 = time.perf_counter()
    out = sc._markup_escape_unstyled(text, validate=False)
    took = time.perf_counter() - t0
    assert out.count("\\[/x]") == 30000
    assert took < 2.0, f"{took:.2f}s"


# ══════════════════════════════════════════════════════════════════════════
# Exported notebook HTML
# ══════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("body", ["[a](\x02//evil.example/x.png)",
                                  "![a](\x02//evil.example/x.png)",
                                  "[a](\x00\x1f//evil.example/x)"])
def test_a_control_character_does_not_smuggle_an_outside_host(body):
    # A browser strips C0 controls from both ends of a URL, so
    # "\x02//evil.example" is the protocol-relative "//evil.example": opening
    # the export fetched from an outside host.
    html = sc._markdown_subset_to_html(body)
    assert "href=" not in html and "src=" not in html, html


def test_a_backtick_in_a_link_target_stays_in_the_target():
    # Code spans are parked first, so the target carried the placeholder
    # (NUL bytes) into the href instead of the body's own text.
    html = sc._markdown_subset_to_html(
        "[a `c` d](https://x.example/`v2`) and ![p `q`](https://x.example/`i`)")
    assert "\x00" not in html
    assert 'href="https://x.example/`v2`"' in html, html
    assert 'src="https://x.example/`i`" alt="p q"' in html, html
    assert "a <code>c</code> d</a>" in html, html


# ══════════════════════════════════════════════════════════════════════════
# A start codon the annotation reads as Met (C5, "follow the file")
# ══════════════════════════════════════════════════════════════════════════

class TestAnAnnotatedInitiatorIsMet:
    """lacI starts GTG and its protein starts M: an initiator tRNA loads fMet
    on a GTG / TTG start. The AA lane, the protein copy and the residue
    planner read the codon table, so they said V1 where the file, UniProt
    and the commercial editor say M. M is shown on the ANNOTATION's word
    only — a thrombin site annotated as a CDS starts CTG, a start codon of
    the standard code, and no ribosome initiates there."""

    PAD = "ACGTACGTAC" * 3

    @staticmethod
    def _laci():
        import splicecraft_presets as P
        return P._find_preset("lacI", "CDS")["sequence"].upper()

    @classmethod
    def _rec(cls, cds, quals, *, strand=1, wrap_at=None):
        from Bio.Seq import Seq
        from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
        from Bio.SeqRecord import SeqRecord
        body = cds if strand == 1 else sc._rc(cds)
        seq = cls.PAD + body + cls.PAD
        s, e = len(cls.PAD), len(cls.PAD) + len(body)
        loc = FeatureLocation(s, e, strand)
        if wrap_at is not None:                 # re-origin inside the gene
            seq = seq[wrap_at:] + seq[:wrap_at]
            n = len(seq)
            s, e = (s - wrap_at) % n, (e - wrap_at) % n
            loc = CompoundLocation([FeatureLocation(s, n, strand),
                                    FeatureLocation(0, e, strand)])
        rec = SeqRecord(Seq(seq), id="t", name="t")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        rec.features = [SeqFeature(loc, type="CDS",
                                   qualifiers={"label": ["lacI"], **quals})]
        return rec

    @classmethod
    def _ncbi(cls, cds, table=11):
        from Bio.Seq import Seq
        return "M" + str(Seq(cds).translate(table=table)).rstrip("*")[1:]

    @staticmethod
    def _lane(rec):
        f = [d for d in sc.PlasmidMap()._parse(rec) if d.get("type") == "CDS"][0]
        aa, _n, _v = sc._cds_aa_list(str(rec.seq).upper(), f)
        return f, "".join(aa)

    @pytest.mark.parametrize("shape", ["plus", "minus", "across-origin"])
    def test_lacI_reads_M_as_its_file_says(self, shape):
        cds = self._laci()
        assert cds[:3] == "GTG"
        rec = self._rec(cds, {"translation": [self._ncbi(cds)],
                              "transl_table": ["11"]},
                        strand=-1 if shape == "minus" else 1,
                        wrap_at=40 if shape == "across-origin" else None)
        f, lane = self._lane(rec)
        assert lane.startswith("MKPVTLYDV"), lane[:10]
        assert f.get("_init_m") == "translation"
        plan = sc._protein_edit_plan(rec, rec.features[0], 1, "M")
        assert (plan["wt_aa"], plan["initiator"]) == ("M", "translation")
        assert plan["silent"] is True             # GTG -> ATG keeps the Met

    def test_a_declared_met_start_reads_M(self):
        cds = "ACG" + self._laci()[3:]            # not a start codon at all
        pos = f"{len(self.PAD) + 1}..{len(self.PAD) + 3}"
        rec = self._rec(cds, {"transl_except": [f"(pos:{pos},aa:Met)"]})
        f, lane = self._lane(rec)
        assert lane[0] == "M" and f.get("_init_m") == "transl_except"

    def test_a_leu_first_site_is_not_turned_into_a_start(self):
        rec = self._rec("CTGGTGCCGCGCGGCAGC", {})    # thrombin site, as a CDS
        f, lane = self._lane(rec)
        assert lane == "LVPRGS" and not f.get("_init_m")

    def test_a_translation_of_something_else_is_not_followed(self):
        # A /translation whose N-terminus is not these bases' (another frame,
        # another start) says nothing about this codon.
        cds = self._laci()
        tr = self._ncbi(cds)
        other = tr[:4] + ("A" if tr[4] != "A" else "G") + tr[5:]
        _f, lane = self._lane(self._rec(cds, {"translation": [other],
                                              "transl_table": ["11"]}))
        assert lane[0] == "V"

    def test_a_mutation_made_downstream_keeps_the_initiator(self):
        # Nothing rewrites /translation when a residue is edited, so judging
        # the whole protein turned lacI's M1 back into V1 after a K100A-style
        # edit that never touched the start.
        cds = self._laci()
        rec = self._rec(cds, {"translation": [self._ncbi(cds)],
                              "transl_table": ["11"]})
        plan = sc._protein_edit_plan(rec, rec.features[0], 100, "W")
        assert not plan["silent"]
        edited = sc._apply_protein_edit(rec, plan)
        f, lane = self._lane(edited)
        assert lane[0] == "M" and lane[99] == "W"
        assert f.get("_init_m") == "translation"

    def test_the_lane_marks_the_initiator(self):
        cds = self._laci()
        rec = self._rec(cds, {"translation": [self._ncbi(cds)],
                              "transl_table": ["11"]})
        f, _lane = self._lane(rec)
        seq = str(rec.seq).upper()
        arr = [(" ", "")] * len(seq)
        sc._paint_cds_aa(arr, f, 0, len(seq), seq, None)
        first = arr[len(self.PAD) + 1]
        second = arr[len(self.PAD) + 4]
        assert first[0] == "M" and sc._INIT_AA_COLOR in first[1]
        assert second[0] == "K" and sc._INIT_AA_COLOR not in second[1]

    def test_copy_and_the_agent_give_the_files_protein(self, monkeypatch):
        cds = self._laci()
        tr = self._ncbi(cds)
        rec = self._rec(cds, {"translation": [tr], "transl_table": ["11"]})
        copied: list = []
        monkeypatch.setattr(sc, "_copy_to_clipboard_with_fallback",
                            lambda app, text, label="": (copied.append(text)
                                                         or ("osc52", "")))
        got: dict = {}

        async def run():
            app = sc.PlasmidApp()
            async with app.run_test(size=(160, 48)) as pilot:
                app._apply_record(rec)
                for _ in range(6):
                    await pilot.pause()
                out = sc._h_find_feature(app, {"name": "lacI"})
                got["m"] = out["matches"][0]
                sp = app.query_one("#seq-panel", sc.SequencePanel)
                pm = app.query_one("#plasmid-map", sc.PlasmidMap)
                sp._aa_highlight = next(d for d in pm._feats
                                        if d.get("type") == "CDS")
                app._copy_strand(bottom=False)
        asyncio.run(run())
        m = got["m"]
        assert m["translation"].startswith("VKPV")    # the codon table, as before
        # Only residue 1 differs; the stop stays exactly as `translation` has it.
        assert m["protein"] == "M" + m["translation"][1:]
        assert m["protein"].rstrip("*") == tr
        assert m["initiator"] == {"codon": "GTG", "codon_reads": "V",
                                  "declared_by": "translation"}
        assert copied == [tr]
