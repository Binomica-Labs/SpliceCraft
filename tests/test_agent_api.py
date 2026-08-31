# pyright: reportArgumentType=false, reportCallIssue=false, reportAttributeAccessIssue=false
#
# Handlers (`_h_*`) return ``dict | tuple[dict, int]`` (success payload OR
# error tuple). Tests routinely unpack one or the other after asserting a
# status code — pyright can't follow the runtime invariant and tags every
# index op as an arg-type mismatch. Negative-input tests on the sanitizer
# helpers (`_sanitize_label`, `_sanitize_feat_type`, …) deliberately pass
# wrong types to verify rejection. `MockApp` stub methods preserve the
# real `PlasmidApp` signature for duck-typing even when the body ignores
# parameters (e.g. `clear_undo`). All three classes of noise are non-bugs
# and would drown out genuine signal; the project's `pyproject.toml`
# already excludes `tests/**` from pyright analysis for the same reason.
# This file-scope pragma keeps the harness diagnostics quiet too.
"""Tests for the agent-API HTTP server (0.4.6+).

Covers the BYO-AI integration that lets an external CLI agent (Claude
Code, Cursor, aider, …) drive the running SpliceCraft GUI via a
localhost JSON API. Two layers:

  * **Pure handler tests** — call `_h_status` / `_h_features` /
    `_h_add_feature` etc. directly with a fake `app` shim. Fast,
    no socket bind, no Textual mount.

  * **End-to-end HTTP tests** — bind a real `_AgentAPIServer` on a
    free port, send `urllib.request` calls, assert the JSON
    response. Uses a `MockApp` with the `_current_record` /
    `_unsaved` / `_apply_record` / `_do_save` surface the handlers
    actually touch.

We deliberately don't spin up a real `PlasmidApp` here — it would add
seconds per test, and the handler logic is what the tests need to
guard. Smoke-level "real app + real port" coverage lives in
`test_smoke.py`.
"""
from __future__ import annotations

import json
import os
import socket
import threading
import time
import types
import urllib.error
import urllib.request

import pytest

import splicecraft as sc


# ── Helpers ────────────────────────────────────────────────────────────────────


def _free_port() -> int:
    """Bind on port 0 to let the OS pick a free port, then close so
    the test server can rebind. Tiny race window, fine for tests."""
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


class MockApp:
    """Stand-in for `PlasmidApp` that gives the handlers everything
    they touch: `_current_record`, `_unsaved`, `_apply_record`,
    `_do_save`, and `call_from_thread` (which we just call inline,
    since there's no Textual loop here)."""

    def __init__(self, record=None):
        self._current_record  = record
        self._unsaved         = False
        self._source_path     = None
        self._restr_min_len   = 6
        self._restr_unique_only = False
        self._show_restr      = False
        self._applied_records: list = []
        self._saved          = False

    def call_from_thread(self, fn, *args, **kwargs):
        # No real event loop in tests — invoke synchronously. Matches
        # the API's return-value semantics (Textual's call_from_thread
        # returns the callable's result).
        return fn(*args, **kwargs)

    def _apply_record(self, record, *, clear_undo=True):
        self._applied_records.append(record)
        self._current_record = record
        self._unsaved = False

    def _apply_annotation_transfers(self, transfers, target_record):
        """Minimal mirror of `PlasmidApp._apply_annotation_transfers`: append
        one feature per accepted transfer (using its target span) and reseat
        via `_apply_record`, so the transfer-annotations apply path can be
        exercised here. The real method also handles wrap features; the test
        records don't wrap."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from copy import deepcopy
        new_record = deepcopy(target_record)
        n = len(new_record.seq)
        for tr in transfers:
            t_s, t_e = int(tr["target_start"]), int(tr["target_end"])
            loc = (FeatureLocation(t_s, t_e,
                                   strand=int(tr.get("target_strand", 1)))
                   if t_e > t_s else FeatureLocation(t_s, n))
            new_record.features.append(SeqFeature(
                loc, type=tr.get("feature_type") or "misc_feature",
                qualifiers={"label": [tr.get("label", "")]}))
        self._apply_record(new_record)

    def _do_save(self):
        self._saved = True
        self._unsaved = False
        return True

    def _push_undo(self):
        pass

    def _mark_dirty(self):
        self._unsaved = True

    def _mark_clean(self):
        # Mirror of the real controller's `mark_clean` effect on the dirty
        # flag (the part the agent dirty-guard reads); no title/panel/autosave
        # side-effects to mock here.
        self._unsaved = False

    def _discard_changes(self):
        # Minimal mirror: the real method reverts to the library-stored copy,
        # or just marks clean when the record has no entry. The test library
        # is empty, so the no-entry → mark-clean path is what's exercised.
        self._unsaved = False

    def _notify_success(self, msg, **kwargs):
        pass

    def _annotate_with_feature(self, start, end, entry):
        """Stub mirror of `PlasmidApp._annotate_with_feature`. Validates
        the same way (range check, zero-length reject, strand coercion)
        and appends a real BioPython SeqFeature to the record so tests
        that assert on `record.features[-1]` see the same shape they
        would in the running GUI. No panel refresh — there's no Textual
        loop here to refresh against."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        from copy import deepcopy
        if self._current_record is None:
            raise RuntimeError("Load a plasmid first.")
        n = len(self._current_record.seq)
        if not (0 <= start < n):
            raise ValueError(f"start {start} out of range [0, {n})")
        if not (0 <= end <= n):
            raise ValueError(f"end {end} out of range [0, {n}]")
        if end == start:
            raise ValueError("zero-length feature (end == start)")
        try:
            strand = int(entry.get("strand", 1))
        except (TypeError, ValueError):
            strand = 1
        biop_strand = strand if strand in (-1, 1) else None
        if end > start:
            loc = FeatureLocation(start, end, strand=biop_strand)
        else:
            loc = CompoundLocation([
                FeatureLocation(start, n, strand=biop_strand),
                FeatureLocation(0, end, strand=biop_strand),
            ])
        feat_type = entry.get("feature_type") or "misc_feature"
        qualifiers: dict = {
            k: list(v) if isinstance(v, (list, tuple)) else [v]
            for k, v in (entry.get("qualifiers") or {}).items()
        }
        label = (entry.get("name") or "").strip()
        if label and "label" not in qualifiers:
            qualifiers["label"] = [label]
        # Mirror the real `_annotate_with_feature_impl`: strand=2 (double-
        # stranded ◀▶) persists via a SpliceCraft_strand qualifier, and a
        # picked colour writes the ApEinfo_*color pair.
        if strand == 2:
            qualifiers["SpliceCraft_strand"] = ["double"]
        new_color = entry.get("color")
        if isinstance(new_color, str) and new_color.strip():
            qualifiers["ApEinfo_fwdcolor"] = [new_color.strip()]
            qualifiers["ApEinfo_revcolor"] = [new_color.strip()]
        new_feat = SeqFeature(loc, type=feat_type, qualifiers=qualifiers)
        new_rec = deepcopy(self._current_record)
        new_rec.features.append(new_feat)
        self._current_record = new_rec
        self._unsaved = True

    def query_one(self, selector, *args):
        # Handlers that touch `query_one("#plasmid-map", PlasmidMap)`
        # for read-only feature listing — we expose a shim with a
        # `_feats` list pulled straight from the SeqRecord.
        if selector == "#plasmid-map":
            class _PMShim:
                _feats = self._feats_from_record()
                _restr_feats: list = []
                def load_record(self, rec):
                    pass
                def refresh(self):
                    pass
            return _PMShim()
        if selector == "#sidebar":
            class _Sidebar:
                def populate(self, feats):
                    pass
            return _Sidebar()
        if selector == "#seq-panel":
            class _SP:
                def update_seq(self, *a, **k):
                    pass
            return _SP()
        from textual.css.query import NoMatches
        raise NoMatches(selector)

    def _feats_from_record(self):
        rec = self._current_record
        if rec is None:
            return []
        out = []
        for f in rec.features:
            if f.type == "source":
                continue
            out.append({
                "start":  int(f.location.start),
                "end":    int(f.location.end),
                "type":   f.type,
                "label":  (f.qualifiers.get("label") or [""])[0],
                "strand": f.location.strand or 1,
                "color":  None,
            })
        return out


@pytest.fixture
def tiny_app(tiny_record):
    """`MockApp` pre-loaded with the conftest `tiny_record`."""
    return MockApp(record=tiny_record)


# ── Pure handler tests (no socket, no app) ────────────────────────────────────


class TestStatusHandler:
    def test_empty_when_no_record(self):
        app = MockApp(record=None)
        result = sc._h_status(app, {})
        assert result["loaded"] is False
        assert result["length"] == 0
        assert result["dirty"] is False

    def test_reports_loaded_record(self, tiny_app, tiny_record):
        result = sc._h_status(tiny_app, {})
        assert result["loaded"] is True
        assert result["name"]   == tiny_record.name
        assert result["length"] == len(tiny_record.seq)
        assert result["version"] == sc.__version__

    def test_reports_dirty_flag(self, tiny_app):
        tiny_app._unsaved = True
        assert sc._h_status(tiny_app, {})["dirty"] is True

    def test_name_is_display_name_not_locus_slug(self, tiny_app, tiny_record):
        # [INV-98] Agent-API audit (2026-06): status must report the DISPLAY
        # name (spaces / case preserved), not the underscored 16-char GenBank
        # LOCUS slug `rec.name`. The display name is the valid `load-entry`
        # key; the slug is not. Was: `status.name` returned `Zephyr_Reporter_`
        # for a `Zephyr Reporter Construct` entry even after a clean load.
        tiny_record._tui_display_name = "Zephyr Reporter Construct"
        assert sc._h_status(tiny_app, {})["name"] == "Zephyr Reporter Construct"


class TestStatusVersionSkew:
    """Stale-daemon detection (petunia agent-API running-log): `status`
    surfaces the on-disk installed version next to the running one so an
    agent notices when `splicecraft update` upgraded the package under a
    still-running `--agent` daemon (which keeps serving the old code)."""

    @pytest.fixture(autouse=True)
    def _reset_cache(self):
        saved = dict(sc._INSTALLED_VERSION_CACHE)
        sc._INSTALLED_VERSION_CACHE.update(version=None, ts=0.0, inflight=False)
        yield
        sc._INSTALLED_VERSION_CACHE.clear()
        sc._INSTALLED_VERSION_CACHE.update(saved)

    def test_running_version_present_and_aliases_version(self):
        s = sc._h_status(MockApp(), {})
        assert s["running_version"] == sc.__version__ == s["version"]

    def test_no_live_daemon_means_no_subprocess(self, monkeypatch):
        # Pure-handler path: no live app/daemon → installed unknown, never
        # stale, and the background refresh is NOT spawned (keeps the whole
        # test suite from shelling out on every status call). Force the
        # liveness gate to "no app" — an earlier async test may have left a
        # real `_LIVE_APP_REF` set, and that gate is what decides whether the
        # off-thread refresh fires.
        monkeypatch.setattr(sc._state, "_LIVE_APP_REF",
                            types.SimpleNamespace(get=lambda: None))
        s = sc._h_status(MockApp(), {})
        assert s["installed_version"] is None
        assert s["stale"] is False
        assert sc._INSTALLED_VERSION_CACHE["inflight"] is False

    def test_stale_flag_when_installed_differs(self, monkeypatch):
        monkeypatch.setattr(sc, "_installed_version_cached",
                            lambda: "0.0.0-other")
        s = sc._h_status(MockApp(), {})
        assert s["installed_version"] == "0.0.0-other"
        assert s["stale"] is True

    def test_not_stale_when_installed_matches_running(self, monkeypatch):
        monkeypatch.setattr(sc, "_installed_version_cached",
                            lambda: sc.__version__)
        assert sc._h_status(MockApp(), {})["stale"] is False

    def test_blocking_helper_reads_query(self, monkeypatch):
        monkeypatch.setattr(sc, "_query_installed_version",
                            lambda timeout=8.0: "1.2.3")
        assert sc._installed_version_cached(blocking=True) == "1.2.3"

    def test_cache_ttl_serves_without_requery(self, monkeypatch):
        calls = []

        def _q(timeout=8.0):
            calls.append(1)
            return "5.5.5"
        monkeypatch.setattr(sc, "_query_installed_version", _q)
        v1 = sc._installed_version_cached(blocking=True)   # 1 query
        v2 = sc._installed_version_cached()                # within TTL → cache
        assert v1 == v2 == "5.5.5"
        assert len(calls) == 1


class TestToolsHandler:
    def test_lists_registered_endpoints(self):
        result = sc._h_tools(None, {})
        names = {ep["name"] for ep in result["endpoints"]}
        # Spot-check the six starter endpoints.
        for required in ("status", "tools", "features", "fetch",
                          "load-entry", "add-feature", "save"):
            assert required in names, f"missing endpoint {required!r}"

    def test_write_flag_is_correct(self):
        eps = {ep["name"]: ep for ep in sc._h_tools(None, {})["endpoints"]}
        assert eps["status"]["write"]      is False
        assert eps["features"]["write"]    is False
        assert eps["fetch"]["write"]       is True
        assert eps["add-feature"]["write"] is True
        assert eps["save"]["write"]        is True
        # set-setting mutates persisted state; must be token-gated.
        # Regression guard for 2026-05-14 security-audit fix where
        # the @_agent_endpoint decoration was missing `write=True`.
        assert eps["set-setting"]["write"] is True

    def test_tools_handler_returns_bare_shape(self):
        # HANDLER-LAYER check only: `_h_tools` itself must return a bare
        # `{endpoints: [...]}` with no `data` key. This is necessary but
        # NOT sufficient to prove the SERVED response is envelope-free —
        # the `data`/`_stale` fields are added by the dispatcher
        # (`_agent_data_envelope`), a layer this unit call never reaches.
        # The authoritative wire-level guard (that removing the
        # `/tools` + `/healthz` early-returns would let the dispatcher
        # envelope them) lives in `TestEnvelopeContract` below.
        result = sc._h_tools(None, {})
        assert set(result) == {"endpoints"}, (
            "/tools handler must return a bare {endpoints} shape")
        assert isinstance(result["endpoints"], list)

    def test_list_features_alias_resolves(self):
        # `list-features` is a natural-guess alias for `features` (the
        # list-* convention every other collection list follows). Both
        # names must resolve to the SAME handler so a client guessing
        # `list-features` doesn't 404.
        names = {ep["name"] for ep in sc._h_tools(None, {})["endpoints"]}
        assert "list-features" in names
        assert "features" in names
        # Registry values are (fn, write) tuples; the two names are
        # registered by separate decorator applications, so compare the
        # underlying handler FUNCTION (same object) + write flag, not the
        # tuple identity. (Split into separate asserts — a chained
        # `a == b is False` would parse as `(a == b) and (b is False)`,
        # which is not the intended "both are read-only" check.)
        assert (sc._AGENT_HANDLERS["list-features"][0]
                is sc._AGENT_HANDLERS["features"][0])
        assert sc._AGENT_HANDLERS["list-features"][1] is False
        assert sc._AGENT_HANDLERS["features"][1] is False


class TestAddFeatureHandler:
    def test_validates_missing_record(self):
        app = MockApp(record=None)
        result = sc._h_add_feature(app, {"start": 0, "end": 10})
        payload, status = result
        assert status == 422
        assert "no plasmid loaded" in payload["error"]

    def test_validates_missing_start(self, tiny_app):
        result = sc._h_add_feature(tiny_app, {"end": 10})
        payload, status = result
        assert status == 400
        assert "start" in payload["error"]

    def test_validates_zero_length(self, tiny_app):
        result = sc._h_add_feature(tiny_app, {"start": 5, "end": 5})
        payload, status = result
        assert status == 400
        assert "zero-length" in payload["error"]

    def test_validates_out_of_range(self, tiny_app, tiny_record):
        n = len(tiny_record.seq)
        result = sc._h_add_feature(tiny_app, {"start": n + 5, "end": n + 10})
        payload, status = result
        assert status == 400
        assert "out of range" in payload["error"]

    def test_validates_strand(self, tiny_app):
        # 7 is not a valid arrow type (only -1/0/1/2 are).
        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 10, "strand": 7}
        )
        payload, status = result
        assert status == 400
        assert "strand" in payload["error"]

    def test_strand_2_double_stranded_accepted(self, tiny_app):
        # Arrow-type parity with the UI: strand=2 (double-stranded ◀▶) is
        # valid and persists the SpliceCraft_strand qualifier.
        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 9, "strand": 2, "label": "ds"})
        assert isinstance(result, dict) and result["strand"] == 2
        feat = tiny_app._current_record.features[-1]
        assert feat.qualifiers.get("SpliceCraft_strand") == ["double"]

    def test_strand_0_arrowless_accepted(self, tiny_app):
        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 9, "strand": 0, "label": "flat"})
        assert isinstance(result, dict) and result["strand"] == 0

    def test_color_applied(self, tiny_app):
        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 9, "color": "#1f77b4"})
        assert isinstance(result, dict) and result["color"] == "#1f77b4"
        feat = tiny_app._current_record.features[-1]
        assert feat.qualifiers.get("ApEinfo_fwdcolor") == ["#1f77b4"]
        assert feat.qualifiers.get("ApEinfo_revcolor") == ["#1f77b4"]

    def test_invalid_color_rejected(self, tiny_app):
        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 9, "color": "octarine"})
        payload, status = result
        assert status == 400 and "color" in payload["error"]

    def test_qualifiers_applied(self, tiny_app):
        result = sc._h_add_feature(tiny_app, {
            "start": 0, "end": 9,
            "qualifiers": {"gene": "bla", "note": ["a", "b"]}})
        assert isinstance(result, dict)
        feat = tiny_app._current_record.features[-1]
        assert feat.qualifiers.get("gene") == ["bla"]
        assert feat.qualifiers.get("note") == ["a", "b"]

    def test_qualifiers_must_be_dict(self, tiny_app):
        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 9, "qualifiers": ["nope"]})
        payload, status = result
        assert status == 400 and "qualifiers" in payload["error"]

    def test_unknown_keys_echoed(self, tiny_app):
        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 9, "bogus": 1})
        assert isinstance(result, dict) and result["ignored"] == ["bogus"]

    def test_add_on_dirty_canvas_succeeds_without_force(self, tiny_app):
        """P1-10 (agent-API feedback): adding a feature is an additive edit —
        it can't discard unsaved work, so a `*`-dirty canvas accepts the
        add WITHOUT `force`. (The stale-load counter, tested below, is the
        real cross-plasmid protection; the dirty guard only forced
        `force:true` on every call after the first.)"""
        tiny_app._unsaved = True
        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 10, "label": "t"}
        )
        assert isinstance(result, dict), result
        assert result["ok"] is True

    def test_force_is_harmless_noop(self, tiny_app):
        """`force` is no longer required (the dirty guard is gone) but must
        stay ACCEPTED — never 400'd, never echoed under `ignored` — so
        existing scripts that pass it keep working."""
        tiny_app._unsaved = True
        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 10, "label": "t",
                        "force": True}
        )
        assert isinstance(result, dict), result
        assert "force" not in result.get("ignored", [])

    def test_stale_load_counter_rejects(self, tiny_app):
        """Regression guard for 2026-05-17 adversarial audit: a canvas
        swap between handler entry and `_apply` running on the UI thread
        must drop the agent edit with 409. Pre-fix, `_h_add_feature`
        would happily annotate whatever record happened to be on the
        canvas at apply time — wrong molecule corruption."""
        tiny_app._record_load_counter = 0

        # Simulate the race: between handler entry (which captures
        # counter==0) and `_apply` execution, the UI thread loads a new
        # plasmid and bumps `_record_load_counter` to 1.
        orig_call = tiny_app.call_from_thread
        def racy_call(fn, *args, **kwargs):
            tiny_app._record_load_counter = 1
            return orig_call(fn, *args, **kwargs)
        tiny_app.call_from_thread = racy_call

        result = sc._h_add_feature(
            tiny_app, {"start": 0, "end": 10, "label": "t"}
        )
        assert isinstance(result, tuple), result
        payload, status = result
        assert status == 409
        assert "canvas reloaded" in payload["error"]


class TestDeleteUpdateFeatureStaleLoadGuard:
    """Regression guards for 2026-05-17 adversarial audit: the agent
    `delete-feature` and `update-feature` endpoints used to do all
    their work inside `_apply` on the UI thread without first
    capturing `_record_load_counter`. A canvas reload between
    handler entry and the queued `_apply` execution would have the
    handler delete / update feature `idx` of the WRONG molecule —
    silent cross-record data corruption."""

    def _racy_app(self, tiny_app):
        tiny_app._record_load_counter = 0
        orig_call = tiny_app.call_from_thread
        def racy_call(fn, *args, **kwargs):
            tiny_app._record_load_counter = 1
            return orig_call(fn, *args, **kwargs)
        tiny_app.call_from_thread = racy_call
        return tiny_app

    def test_delete_feature_rejects_on_stale_counter(self, tiny_app):
        app = self._racy_app(tiny_app)
        result = sc._h_delete_feature(app, {"idx": 0})
        assert isinstance(result, tuple), result
        payload, status = result
        assert status == 409
        assert "canvas reloaded" in payload["error"]

    def test_update_feature_rejects_on_stale_counter(self, tiny_app):
        app = self._racy_app(tiny_app)
        result = sc._h_update_feature(
            app, {"idx": 0, "label": "should-not-apply"}
        )
        assert isinstance(result, tuple), result
        payload, status = result
        assert status == 409
        assert "canvas reloaded" in payload["error"]


class TestDeleteUpdateFeatureTOCTOUSignature:
    """Sweep #32 (2026-05-26) adversarial audit: delete/update
    feature handlers used to do a bounds check + access inside
    `_apply` on the UI thread, but the agent's `idx` was captured
    on the worker thread. A concurrent `_apply` from another
    request could insert a feature, shifting indices, before
    this request's `_apply` ran — the bounds check then passed
    against the post-insert state but `pm._feats[idx]` referred
    to a DIFFERENT feature than the agent saw. Silent
    cross-feature corruption. Fix captures a signature
    (start/end/type/label) in a pre-flight UI-thread call and
    re-verifies inside `_apply`; mismatch → 409 Conflict."""

    def _make_racy_app(self, tiny_app):
        """Wrap `call_from_thread` so the SECOND call (the
        `_apply` closure, after the pre-flight signature
        capture) sees a record with a NEW feature injected at
        idx 0 — simulating a concurrent agent request that
        inserted a row mid-handler."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        orig_call = tiny_app.call_from_thread
        call_n = [0]
        shadow = SeqFeature(
            FeatureLocation(0, 1, strand=1),
            type="shadow",
            qualifiers={"label": ["injected"]},
        )
        def racy_call(fn, *args, **kwargs):
            call_n[0] += 1
            if call_n[0] == 2:
                # Inject shadow BEFORE the queued `_apply` runs.
                # The mock's `query_one("#plasmid-map", ...)`
                # re-reads `_feats` from the record on every
                # call, so this shift is observable in `_apply`.
                tiny_app._current_record.features = (
                    [shadow]
                    + list(tiny_app._current_record.features)
                )
            return orig_call(fn, *args, **kwargs)
        tiny_app.call_from_thread = racy_call
        return tiny_app

    def test_delete_feature_detects_signature_drift(self, tiny_app):
        app = self._make_racy_app(tiny_app)
        result = sc._h_delete_feature(app, {"idx": 0})
        assert isinstance(result, tuple), result
        payload, status = result
        assert status == 409, (status, payload)
        assert "changed under us" in payload["error"], payload

    def test_update_feature_detects_signature_drift(self, tiny_app):
        app = self._make_racy_app(tiny_app)
        result = sc._h_update_feature(
            app, {"idx": 0, "label": "should-not-apply"},
        )
        assert isinstance(result, tuple), result
        payload, status = result
        assert status == 409
        assert "changed under us" in payload["error"]


class TestUpdateFeatureParity:
    """update-feature reaches the same arrow-type / colour / qualifier
    controls as the edit dialog (the UI-parity gap the user flagged:
    add/update-feature must set TYPE, ARROW TYPE, colour, etc.)."""

    def _app_with_feature(self, tiny_app):
        """Append a feature and return (app, idx, getter). `idx` is the
        appended feature's row in the non-source feature list (tiny_record
        already carries features, so it is NOT 0); `getter()` returns that
        exact SeqFeature after an update."""
        sc._h_add_feature(tiny_app, {"start": 0, "end": 9, "label": "f0",
                                      "force": True})
        tiny_app._unsaved = False
        nonsrc = [f for f in tiny_app._current_record.features
                  if f.type != "source"]
        idx = len(nonsrc) - 1
        # The appended feature is the LAST in record order (source aside).
        getter = lambda: [f for f in tiny_app._current_record.features
                          if f.type != "source"][idx]
        return tiny_app, idx, getter

    def test_update_strand_to_double(self, tiny_app):
        app, idx, feat = self._app_with_feature(tiny_app)
        r = sc._h_update_feature(app, {"idx": idx, "strand": 2, "force": True})
        assert isinstance(r, dict), r
        assert feat().qualifiers.get("SpliceCraft_strand") == ["double"]

    def test_update_strand_away_from_double_clears_qualifier(self, tiny_app):
        app, idx, feat = self._app_with_feature(tiny_app)
        sc._h_update_feature(app, {"idx": idx, "strand": 2, "force": True})
        sc._h_update_feature(app, {"idx": idx, "strand": 1, "force": True})
        assert "SpliceCraft_strand" not in feat().qualifiers

    def test_update_color_then_clear(self, tiny_app):
        app, idx, feat = self._app_with_feature(tiny_app)
        sc._h_update_feature(app, {"idx": idx, "color": "#abc", "force": True})
        assert feat().qualifiers.get("ApEinfo_fwdcolor") == ["#abc"]
        sc._h_update_feature(app, {"idx": idx, "color": "", "force": True})
        assert "ApEinfo_fwdcolor" not in feat().qualifiers

    def test_update_invalid_color_rejected(self, tiny_app):
        app, idx, _ = self._app_with_feature(tiny_app)
        r = sc._h_update_feature(
            app, {"idx": idx, "color": "burgundy", "force": True})
        assert isinstance(r, tuple) and r[1] == 400

    def test_update_qualifiers_merge(self, tiny_app):
        app, idx, feat = self._app_with_feature(tiny_app)
        sc._h_update_feature(
            app, {"idx": idx, "qualifiers": {"gene": "x", "note": "y"},
                  "force": True})
        assert feat().qualifiers.get("gene") == ["x"]
        assert feat().qualifiers.get("note") == ["y"]

    def test_update_invalid_strand_rejected(self, tiny_app):
        app, idx, _ = self._app_with_feature(tiny_app)
        r = sc._h_update_feature(app, {"idx": idx, "strand": 9, "force": True})
        assert isinstance(r, tuple) and r[1] == 400


class TestNewPlasmidHandler:
    """`new-plasmid` — create a plasmid from a raw sequence (the agent
    counterpart to Ctrl+N), loaded onto the canvas but NOT auto-saved."""

    def test_creates_and_loads(self):
        app = MockApp(record=None)
        r = sc._h_new_plasmid(app, {"name": "My Construct v1",
                                     "sequence": "ATGAAA" * 8,
                                     "circular": True})
        assert r["ok"] and r["saved"] is False
        assert r["topology"] == "circular" and r["length"] == 48
        # Spaced display name preserved (not the underscored LOCUS).
        assert sc.PlasmidApp._record_display_name(
            app._current_record) == "My Construct v1"

    def test_features_applied(self):
        app = MockApp(record=None)
        r = sc._h_new_plasmid(app, {
            "name": "p", "sequence": "ATGAAATAG" * 4,
            "features": [{"start": 0, "end": 9, "name": "orf",
                          "feature_type": "CDS"}]})
        assert r["n_features"] == 1

    def test_missing_name_400(self):
        assert sc._h_new_plasmid(
            MockApp(record=None), {"sequence": "ATGC"})[1] == 400

    def test_invalid_sequence_400(self):
        assert sc._h_new_plasmid(
            MockApp(record=None), {"name": "x", "sequence": "ZZZZ"})[1] == 400

    def test_empty_sequence_400(self):
        assert sc._h_new_plasmid(
            MockApp(record=None), {"name": "x", "sequence": ""})[1] == 400

    def test_features_must_be_list(self):
        assert sc._h_new_plasmid(MockApp(record=None), {
            "name": "x", "sequence": "ATGC", "features": "nope"})[1] == 400

    def test_dirty_guard_without_force(self):
        app = MockApp(record=None)
        app._unsaved = True
        assert sc._h_new_plasmid(
            app, {"name": "x", "sequence": "ATGC"})[1] == 409

    def test_unknown_keys_echoed(self):
        r = sc._h_new_plasmid(MockApp(record=None), {
            "name": "x", "sequence": "ATGC", "bogus": 1})
        assert r["ignored"] == ["bogus"]


_SENTINEL = object()


class _UndoMockApp:
    """Minimal undo/redo surface for the wrapper-logic tests. The REAL
    snapshot semantics are covered by the UndoController tests; here we only
    pin the endpoint's 422 / 409 / stack-delta reporting."""

    def __init__(self, undo=2, redo=0, blocked=False, record=_SENTINEL):
        if record is _SENTINEL:
            record = types.SimpleNamespace(seq="ATGCATGCAT")
        self._current_record = record
        self._undo_stack = list(range(undo))
        self._redo_stack = list(range(redo))
        self._blocked = blocked

    def _undo_blocked_by_modal(self):
        return self._blocked

    def _action_undo(self):
        if self._undo_stack:
            self._redo_stack.append(self._undo_stack.pop())

    def _action_redo(self):
        if self._redo_stack:
            self._undo_stack.append(self._redo_stack.pop())

    def call_from_thread(self, fn, *a, **k):
        return fn(*a, **k)


class TestUndoRedoHandlers:
    def test_undo_no_record_422(self):
        app = _UndoMockApp(record=None)
        assert sc._h_undo(app, {})[1] == 422

    def test_redo_no_record_422(self):
        app = _UndoMockApp(record=None)
        assert sc._h_redo(app, {})[1] == 422

    def test_undo_reports_delta(self):
        app = _UndoMockApp(undo=2, redo=0)
        r = sc._h_undo(app, {})
        assert r["ok"] and r["undone"] is True
        assert r["undo_remaining"] == 1 and r["redo_available"] == 1

    def test_undo_empty_stack_is_noop_not_error(self):
        app = _UndoMockApp(undo=0)
        r = sc._h_undo(app, {})
        assert r["ok"] and r["undone"] is False and r["undo_remaining"] == 0

    def test_undo_blocked_409(self):
        app = _UndoMockApp(undo=2, blocked=True)
        assert sc._h_undo(app, {})[1] == 409

    def test_redo_reports_delta(self):
        app = _UndoMockApp(undo=0, redo=2)
        r = sc._h_redo(app, {})
        assert r["ok"] and r["redone"] is True
        assert r["redo_remaining"] == 1 and r["undo_available"] == 1

    def test_redo_blocked_409(self):
        app = _UndoMockApp(redo=1, blocked=True)
        assert sc._h_redo(app, {})[1] == 409


class TestMultiAlignHandler:
    """`multi-align` — batch pairwise alignment (the Alt+A overlay)."""

    def test_explicit_sequences(self):
        r = sc._h_multi_align(MockApp(record=None), {
            "query": "ATGCATGCATGCATGCATGC",
            "targets": [{"sequence": "ATGCATGCATGCATGCATGC", "name": "same"},
                        {"sequence": "TTTTTTTTTTGGGGGGGGGG", "name": "diff"}]})
        assert r["ok"] and len(r["alignments"]) == 2
        assert r["alignments"][0]["identity_pct"] == 100.0
        # The big gapped strings are omitted to keep the batch small.
        assert "aligned_q" not in r["alignments"][0]

    def test_per_target_errors_dont_fail_batch(self):
        r = sc._h_multi_align(MockApp(record=None), {
            "query": "ATGCATGC",
            "targets": [{"sequence": "ATGCATGC"}, {"sequence": "ZZZ"},
                        {"name": "no-such-lib-entry"}]})
        assert r["alignments"][0]["identity_pct"] == 100.0
        assert "error" in r["alignments"][1]   # non-IUPAC
        assert "error" in r["alignments"][2]   # unresolvable

    def test_no_query_no_record_422(self):
        assert sc._h_multi_align(MockApp(record=None),
                                 {"targets": [{"sequence": "AT"}]})[1] == 422

    def test_empty_targets_400(self):
        assert sc._h_multi_align(MockApp(record=None),
                                 {"query": "AT", "targets": []})[1] == 400

    def test_too_many_targets_400(self):
        assert sc._h_multi_align(MockApp(record=None), {
            "query": "AT", "targets": [{"sequence": "AT"}] * 21})[1] == 400

    def test_bad_mode_400(self):
        assert sc._h_multi_align(MockApp(record=None), {
            "query": "AT", "targets": [{"sequence": "AT"}],
            "mode": "diagonal"})[1] == 400

    def test_uses_loaded_record_as_query(self, tiny_app, tiny_record):
        r = sc._h_multi_align(tiny_app, {
            "targets": [{"sequence": str(tiny_record.seq), "name": "self"}]})
        assert r["alignments"][0]["identity_pct"] == 100.0

    def test_origin_shifted_target_is_rotation_aware(self):
        """A circular target whose origin differs from the query's must
        align on its shared backbone via the rotation picker — not smear
        off bp 0 like bare `_pairwise_align` did (2026-07-01 "ori shows
        no alignment" report). Raw-sequence targets default to circular,
        and the row now carries the rotation metadata."""
        import random
        seq = "".join(random.Random(31).choices("ACGT", k=1000))
        shifted = seq[400:] + seq[:400]   # same molecule, origin +400
        r = sc._h_multi_align(MockApp(record=None), {
            "query": seq,
            "targets": [{"sequence": shifted, "name": "origin-shift"}]})
        row = r["alignments"][0]
        assert row["ungapped_identity_pct"] == 100.0
        assert row["picked_rotation"] in ("query", "target")
        assert row["circular"] is True

    def test_circular_false_disables_rotation(self):
        """`circular: false` forces the naive (non-rotation) path — an
        origin-shifted pair then can't reach 100%, matching the old
        behaviour for callers that relied on it."""
        import random
        seq = "".join(random.Random(31).choices("ACGT", k=1000))
        shifted = seq[400:] + seq[:400]
        r = sc._h_multi_align(MockApp(record=None), {
            "query": seq, "circular": False,
            "targets": [{"sequence": shifted, "name": "origin-shift"}]})
        row = r["alignments"][0]
        assert row["picked_rotation"] == "none"
        assert row["ungapped_identity_pct"] < 100.0

    def test_empty_loaded_record_query_422(self):
        """A loaded record with an empty sequence fails once, upfront,
        with a clear 422 — not N identical per-target picker errors."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        empty_rec = SeqRecord(Seq(""), id="empty", name="empty")
        r = sc._h_multi_align(MockApp(record=empty_rec),
                              {"targets": [{"sequence": "ATGC"}]})
        assert r[1] == 422 and "empty" in r[0]["error"].lower()


class TestAttachExperimentImage:
    """`attach-experiment-image` — server-side image attach to a notebook
    entry (the conftest fixture sandboxes + authorises writes)."""

    def _png(self, tmp_path):
        p = tmp_path / "fig.png"
        p.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 64)
        return p

    def test_attach_roundtrip(self, tmp_path):
        ce = sc._h_create_experiment(None, {"title": "E", "body_md": "x"})
        r = sc._h_attach_experiment_image(None, {
            "experiment_id": ce["id"], "path": str(self._png(tmp_path))})
        assert r["ok"] and r["filename"].endswith(".png")
        exp = next(e for e in sc._load_experiments() if e["id"] == ce["id"])
        assert "![" in exp["body_md"]   # markdown ref embedded

    def test_missing_id_400(self, tmp_path):
        assert sc._h_attach_experiment_image(
            None, {"path": str(self._png(tmp_path))})[1] == 400

    def test_unknown_experiment_404(self, tmp_path):
        assert sc._h_attach_experiment_image(None, {
            "experiment_id": "nope", "path": str(self._png(tmp_path))})[1] == 404

    def test_non_image_rejected_400(self, tmp_path):
        ce = sc._h_create_experiment(None, {"title": "E"})
        txt = tmp_path / "x.txt"
        txt.write_text("not an image")
        assert sc._h_attach_experiment_image(None, {
            "experiment_id": ce["id"], "path": str(txt)})[1] == 400


class TestTraditionalCloning:
    """`simulate-traditional-cloning` + `traditional-clone` — a real digest +
    ligation that NEVER picks the insert by fragment size (the catastrophic-
    class never-assume-the-smaller-fragment rule)."""

    # EcoRI · backbone(A×200) · BamHI · stuffer(C×60), as a circular vector.
    VECTOR = "GAATTC" + "A" * 200 + "GGATCC" + "C" * 60
    # pad · EcoRI · payload · BamHI · pad — a PCR insert with primer-added sites.
    INSERT = "GCGC" + "GAATTC" + "ATGAAACCCGGGTTTTAA" + "GGATCC" + "GCGC"
    PAYLOAD = "ATGAAACCCGGGTTTTAA"

    def _base(self):
        return {"vector_seq": self.VECTOR,
                "vector_enzymes": ["EcoRI", "BamHI"],
                "insert_seq": self.INSERT,
                "insert_enzymes": ["EcoRI", "BamHI"]}

    def test_simulate_lists_both_fragments(self):
        r = sc._h_simulate_traditional_cloning(None, self._base())
        assert r["ok"] and r["n_vector_fragments"] == 2
        # The insert ligates into BOTH fragments — both reported, neither
        # silently dropped because it's the smaller piece.
        assert len(r["products"]) == 2
        assert any(p["length"] == 230 for p in r["products"])   # backbone+insert
        # Read response omits the verbose per-feature list.
        assert "features" not in r["products"][0]

    def test_simulate_bad_enzyme_400(self):
        b = self._base()
        b["vector_enzymes"] = ["NotAnEnzyme"]
        assert sc._h_simulate_traditional_cloning(None, b)[1] == 400

    def test_simulate_type_iis_insert_422(self):
        b = self._base()
        b["insert_enzymes"] = ["BsaI"]
        assert sc._h_simulate_traditional_cloning(None, b)[1] == 422

    def test_simulate_no_cuts_422(self):
        b = self._base()
        b["vector_seq"] = "A" * 100   # no EcoRI / BamHI sites
        assert sc._h_simulate_traditional_cloning(None, b)[1] == 422

    def test_clone_ambiguous_409(self):
        # Both fragments accept the insert → refuse, list the options.
        r = sc._h_traditional_clone(None, dict(self._base(), product_name="c"))
        assert isinstance(r, tuple) and r[1] == 409
        assert "options" in r[0] and len(r[0]["options"]) == 2

    def test_clone_disambiguated_saves_circular_golden(self):
        r = sc._h_traditional_clone(
            None, dict(self._base(), product_name="myclone", vector_frag_idx=0))
        assert r["ok"] and r["vector_frag_idx"] == 0
        ent = next(e for e in sc._load_library() if e["name"] == "myclone")
        assert ent["source"] == "agent:traditional"
        rec = sc._gb_text_to_record(ent["gb_text"])
        assert rec.annotations.get("topology") == "circular"
        # GOLDEN: the saved product carries the insert payload AND the
        # backbone — a genuine assembly, not a hand-built final ([INV-127]).
        seq = str(rec.seq).upper()
        assert self.PAYLOAD in seq or self.PAYLOAD in sc._rc(seq)
        assert "A" * 150 in seq or "A" * 150 in sc._rc(seq)

    def test_clone_bad_orientation_400(self):
        r = sc._h_traditional_clone(
            None, dict(self._base(), vector_frag_idx=0, orientation="sideways"))
        assert isinstance(r, tuple) and r[1] == 400

    def test_clone_no_ligation_422(self):
        # Vector cut only by EcoRI, insert only by BamHI → no matching ends.
        r = sc._h_traditional_clone(None, {
            "vector_seq": self.VECTOR, "vector_enzymes": ["EcoRI"],
            "insert_seq": "GCGC" + "GGATCC" + "ATGCATGC" + "GGATCC" + "GCGC",
            "insert_enzymes": ["BamHI"], "product_name": "x",
            "vector_frag_idx": 0})
        assert isinstance(r, tuple) and r[1] == 422

    # ── carry_annotations (agent-API feedback) ─────────────────────────────
    def _seed_vector_entry(self, name="MyVector", label="backbone_marker",
                            seq=None):
        """Save a vector library entry carrying one misc_feature over the
        backbone, so carry_annotations has something to source. ``seq`` lets a
        test store a MISMATCHED sequence to exercise the soft-skip gate."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        s = self.VECTOR if seq is None else seq
        rec = SeqRecord(Seq(s), id="MYVEC", name="MYVEC",
                        annotations={"molecule_type": "DNA",
                                       "topology": "circular"})
        rec.features.append(SeqFeature(
            FeatureLocation(10, 200, strand=1), type="misc_feature",
            qualifiers={"label": [label]}))
        sc._save_library([{
            "id": "MYVEC", "name": name, "size": len(s), "n_feats": 1,
            "added": "2026-07-02", "gb_text": sc._record_to_gb_text(rec)}])
        return label

    def test_carry_annotations_lifts_vector_features(self, isolated_library):
        # vector_frag_idx=0 is the A×200 backbone; its misc_feature must ride
        # onto the ligated product (no follow-up transfer-annotations needed).
        label = self._seed_vector_entry()
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="carried", vector_frag_idx=0,
            carry_annotations=True, vector_name="MyVector"))
        assert r["ok"], r
        assert r["annotations_carried"] is True
        assert not r["carry_warnings"], r["carry_warnings"]
        ent = next(e for e in sc._load_library() if e["name"] == "carried")
        rec = sc._gb_text_to_record(ent["gb_text"])
        labels = [(f.qualifiers.get("label") or [""])[0] for f in rec.features]
        assert label in labels, labels

    def test_carry_annotations_reported_in_simulate(self, isolated_library):
        self._seed_vector_entry()
        r = sc._h_simulate_traditional_cloning(None, dict(
            self._base(), carry_annotations=True, vector_name="MyVector"))
        assert r["ok"] and r["annotations_carried"] is True

    def test_carry_annotations_requires_a_name(self):
        # Asked to carry but gave nothing to source from → fail loud (never a
        # silent ok:true with a feature-bare product).
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="x", vector_frag_idx=0,
            carry_annotations=True))
        assert isinstance(r, tuple) and r[1] == 400

    def test_carry_annotations_unknown_name_404(self, isolated_library):
        sc._save_library([])
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="x", vector_frag_idx=0,
            carry_annotations=True, vector_name="ghost"))
        assert isinstance(r, tuple) and r[1] == 404

    def test_carry_annotations_seq_mismatch_soft_skips(self, isolated_library):
        # The entry exists but its stored sequence differs from the passed
        # vector_seq → carrying rotated/edited coords would mis-place features,
        # so skip that side with a warning (the clone still succeeds).
        self._seed_vector_entry(seq="ACGT" * 80)
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="mismatch", vector_frag_idx=0,
            carry_annotations=True, vector_name="MyVector"))
        assert r["ok"], r
        assert r["annotations_carried"] is False
        assert r["carry_warnings"], "expected a carry warning on seq mismatch"


class TestTraditionalCloningCircularInsert:
    """`insert_circular` + `insert_frag_idx` — cut a cassette OUT of a CIRCULAR
    plasmid (e.g. an Ω multigene → binary vector) and pick which fragment is the
    cassette. The symmetric counterpart to `vector_frag_idx` (verified against a
    real GoldenBraid Ω→pCambia binary build off-suite)."""
    # circular vector: EcoRI · backbone(A×200) · KpnI · stuffer(C×40)
    VECTOR = "GAATTC" + "A" * 200 + "GGTACC" + "C" * 40
    PAYLOAD = "ATGAAACCCGGGTTTAAACCCGGG" * 3            # 72 bp cassette
    # ins-backbone: a DISTINCTIVE non-homopolymer marker (a homopolymer's
    # reverse-complement would alias the vector's A-run and confuse the check).
    INS_BACKBONE = "TCGATCGATCGA" * 7                   # 84 bp, no EcoRI/KpnI
    # circular insert PLASMID: EcoRI · cassette · KpnI · ins-backbone.
    # EcoRI+KpnI cuts it into TWO fragments (the cassette AND the ins-backbone).
    INSERT = "GAATTC" + PAYLOAD + "GGTACC" + INS_BACKBONE

    def _base(self):
        return {"vector_seq": self.VECTOR, "vector_enzymes": ["EcoRI", "KpnI"],
                "insert_seq": self.INSERT, "insert_enzymes": ["EcoRI", "KpnI"],
                "insert_circular": True}

    def test_simulate_enumerates_both_insert_fragments(self):
        r = sc._h_simulate_traditional_cloning(None, self._base())
        assert r["ok"] and r["n_insert_fragments"] == 2
        assert len(r["insert_fragments"]) == 2
        # products carry BOTH fragment indices — neither insert piece dropped.
        assert {p["insert_frag_idx"] for p in r["products"]} == {0, 1}

    def test_clone_ambiguous_lists_insert_frag_idx(self):
        # >1 candidate (2 vec × 2 ins frags) → refuse + list, never guess.
        r = sc._h_traditional_clone(None, dict(self._base(), product_name="amb"))
        assert isinstance(r, tuple) and r[1] == 409
        assert all("insert_frag_idx" in o for o in r[0]["options"])

    def test_clone_picks_cassette_not_ins_backbone(self):
        # Find the product that carries the cassette ON the vector backbone,
        # then build exactly that (vector_frag_idx + insert_frag_idx + orient).
        sim = sc._h_simulate_traditional_cloning(None, self._base())
        tgt = next(p for p in sim["products"]
                   if (self.PAYLOAD in p["sequence"]
                       or self.PAYLOAD in sc._rc(p["sequence"]))
                   and ("A" * 150 in p["sequence"]
                        or "A" * 150 in sc._rc(p["sequence"])))
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="circ_cassette_clone",
            vector_frag_idx=tgt["vector_frag_idx"],
            insert_frag_idx=tgt["insert_frag_idx"],
            orientation=tgt["orientation"]))
        assert isinstance(r, dict) and r["ok"], r
        assert r["insert_frag_idx"] == tgt["insert_frag_idx"]
        ent = next(e for e in sc._load_library() if e["name"] == "circ_cassette_clone")
        seq = str(sc._gb_text_to_record(ent["gb_text"]).seq).upper()
        # cassette IN, ins-backbone OUT — the WRONG insert fragment wasn't used.
        assert self.PAYLOAD in seq or self.PAYLOAD in sc._rc(seq)
        assert (self.INS_BACKBONE not in seq
                and self.INS_BACKBONE not in sc._rc(seq))

    def test_linear_insert_default_single_candidate(self):
        # insert_circular defaults False → the linear PCR insert still yields a
        # single candidate (insert_frag_idx 0); back-compat with old callers.
        r = sc._h_simulate_traditional_cloning(None, {
            "vector_seq": self.VECTOR, "vector_enzymes": ["EcoRI", "KpnI"],
            "insert_seq": "GCGC" + "GAATTC" + self.PAYLOAD + "GGTACC" + "GCGC",
            "insert_enzymes": ["EcoRI", "KpnI"]})
        assert r["n_insert_fragments"] == 1
        assert all(p["insert_frag_idx"] == 0 for p in r["products"])


class TestDiscardChangesHandler:
    """`discard-changes` — agent recovery for a stuck-dirty canvas
    (agent-API feedback P1-10)."""

    def test_no_record_422(self):
        assert sc._h_discard_changes(MockApp(record=None), {})[1] == 422

    def test_marks_clean_when_no_entry(self, tiny_app, isolated_library):
        sc._save_library([])           # tiny_record is NOT in the library
        tiny_app._unsaved = True
        r = sc._h_discard_changes(tiny_app, {})
        assert r["ok"] is True
        assert r["dirty"] is False     # canvas no longer dirty
        assert r["reverted"] is False  # nothing stored to revert to


class TestAddFeaturesBatch:
    """`add-features` — batch insert under one lock + one stale-load check
    (agent-API feedback P1-10)."""

    def test_no_record_422(self):
        assert sc._h_add_features(
            MockApp(record=None), {"features": [{"start": 0, "end": 3}]})[1] == 422

    def test_features_must_be_nonempty_list(self, tiny_app):
        assert sc._h_add_features(tiny_app, {})[1] == 400
        assert sc._h_add_features(tiny_app, {"features": []})[1] == 400
        assert sc._h_add_features(tiny_app, {"features": "x"})[1] == 400

    def test_batch_adds_all_in_one_call(self, tiny_app, tiny_record):
        before = len([f for f in tiny_record.features if f.type != "source"])
        r = sc._h_add_features(tiny_app, {"features": [
            {"start": 0, "end": 9, "label": "a"},
            {"start": 9, "end": 18, "label": "b", "strand": -1},
        ]})
        assert r["ok"] is True and r["added"] == 2
        assert all(f["ok"] for f in r["features"])
        after = len([f for f in tiny_app._current_record.features
                     if f.type != "source"])
        assert after == before + 2

    def test_malformed_item_rejects_whole_call(self, tiny_app):
        # A caller bug (non-int start) fails the call before ANY feature lands.
        r = sc._h_add_features(tiny_app, {"features": [{"start": "x", "end": 3}]})
        assert isinstance(r, tuple) and r[1] == 400

    def test_out_of_range_item_reported_per_item(self, tiny_app):
        # A valid item lands; a biologically-rejected span is reported per-item
        # without aborting the rest.
        r = sc._h_add_features(tiny_app, {"features": [
            {"start": 0, "end": 9, "label": "ok"},
            {"start": 0, "end": 10 ** 9, "label": "bad"},
        ]})
        assert r["ok"] is True and r["added"] == 1
        assert r["features"][0]["ok"] is True
        assert r["features"][1]["ok"] is False and "error" in r["features"][1]


class TestDeletePrimersBatch:
    """`delete-primers` — one-call batch prune (agent-API feedback SC-L)."""

    def test_names_must_be_nonempty_list(self):
        assert sc._h_delete_primers(None, {})[1] == 400
        assert sc._h_delete_primers(None, {"names": []})[1] == 400

    def test_batch_delete_by_name(self):
        sc._save_primers([])
        for nm, sq in (("b1", "ACGTACGTAA"), ("b2", "ACGTACGTCC"),
                       ("b3", "ACGTACGTGG")):
            sc._h_create_primer(None, {"name": nm, "sequence": sq})
        r = sc._h_delete_primers(None, {"names": ["b1", "b3", "nope"]})
        assert r["ok"] is True and r["removed"] == 2
        assert set(r["removed_names"]) == {"b1", "b3"}
        assert r["not_found"] == ["nope"]
        remaining = {p["name"] for p in sc._load_primers()}
        assert remaining == {"b2"}

    def test_ambiguous_name_not_deleted(self):
        sc._save_primers([])
        sc._h_create_primer(None, {"name": "dup", "sequence": "AAAACCCCGG"})
        sc._h_create_primer(None, {"name": "dup", "sequence": "TTTTGGGGAA"})
        r = sc._h_delete_primers(None, {"names": ["dup"]})
        assert r["ambiguous"] == ["dup"] and r["removed"] == 0
        assert len(sc._load_primers()) == 2   # both kept (ambiguous → skipped)

    def test_unknown_collection_404(self):
        r = sc._h_delete_primers(None, {"names": ["x"], "collection": "NoSuch"})
        assert isinstance(r, tuple) and r[1] == 404


class TestDomesticatePartsBatch:
    """`domesticate-parts` — batch domestication (agent-API feedback P1-6)."""

    def test_parts_must_be_nonempty_list(self):
        assert sc._h_domesticate_parts(MockApp(), {})[1] == 400
        assert sc._h_domesticate_parts(MockApp(), {"parts": []})[1] == 400

    def test_per_item_failure_without_entry_vector(self, isolated_library):
        # No entry vector configured → each part fails 422, reported per-item;
        # the batch itself returns ok with built=0 (it does NOT abort).
        sc._save_entry_vectors([])
        app = MockApp()
        parts = [
            {"sequence": "ATGAAACGTTAA", "oh5": "GGAG", "oh3": "AATG", "name": "p1"},
            {"sequence": "ATGGGGCCCTAA", "oh5": "GGAG", "oh3": "AATG", "name": "p2"},
        ]
        r = sc._h_domesticate_parts(app, {"parts": parts})
        assert r["ok"] is True
        assert r["built"] == 0 and r["total"] == 2
        assert all(item["ok"] is False for item in r["results"])
        assert all(item.get("code") == 422 for item in r["results"])


class TestApiTickets:
    """Regression coverage for the ProjectA constructor tickets
    (SPLICECRAFT_API_TICKETS.md): SC-B / SC-D / SC-G / SC-I."""

    def _h(self, name):
        return sc._state._AGENT_HANDLERS[name][0]

    # ── SC-B: delete-primer by name/id ──────────────────────────────────
    def test_sc_b_delete_primer_by_name(self):
        self._h("create-primer")(None, {"name": "oF1",
                                          "sequence": "ACGTACGTAC"})
        assert self._h("delete-primer")(None, {"name": "oF1"})["ok"]
        assert self._h("delete-primer")(None, {"name": "oF1"})[1] == 404

    def test_sc_b_delete_primer_ambiguous_409(self):
        self._h("create-primer")(None, {"name": "dup",
                                         "sequence": "AAAACCCCGG"})
        self._h("create-primer")(None, {"name": "dup",
                                         "sequence": "TTTTGGGGAA"})
        r = self._h("delete-primer")(None, {"name": "dup"})
        assert isinstance(r, tuple) and r[1] == 409
        assert len(r[0]["sequences"]) == 2

    def test_sc_b_delete_primer_missing_all_400(self):
        assert self._h("delete-primer")(None, {})[1] == 400

    # ── SC-D: create-primer honors `collection` (no silent MAIN landing) ──
    def test_sc_d_create_primer_honors_collection(self):
        self._h("create-primer-collection")(None, {"name": "ProjectA"})
        r = self._h("create-primer")(None, {
            "name": "pF", "sequence": "GGGGAAAACC", "collection": "ProjectA"})
        assert r["collection"] == "ProjectA"
        coll = next(c for c in sc._load_primer_collections()
                    if c["name"] == "ProjectA")
        assert any(p["name"] == "pF" for p in coll["primers"])
        # And NOT in the default/active collection.
        assert not any(p.get("name") == "pF" for p in sc._load_primers())

    def test_sc_d_create_primer_unknown_collection_404(self):
        r = self._h("create-primer")(None, {
            "name": "x", "sequence": "CCCCGGGGAA", "collection": "NoSuch"})
        assert isinstance(r, tuple) and r[1] == 404   # loud, not ok:true

    # ── SC-G: one concept, one key (aliases accepted) ───────────────────
    def test_sc_g_sequence_alias_on_rbs_strength(self):
        r = self._h("rbs-strength")(None, {
            "sequence": "AGGAGGACAACAATGAAACGT", "start": 13})
        assert isinstance(r, dict) and r.get("ok")

    def test_sc_g_sequence_alias_on_design_mutagenesis(self):
        r = self._h("design-mutagenesis")(None, {
            "sequence": "ATG" + "GCA" * 40 + "TAA", "mutation": "A5G"})
        # Accepted (not a 400 about a missing cds_dna).
        assert isinstance(r, dict) or r[1] != 400

    # ── SC-I: copy-plasmid into a (non-active) collection ───────────────
    def _seed_collections(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 20), id="FFE1", name="FFE_1_ENTRY",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec._tui_display_name = "FFE 1 ENTRY"
        gb = sc._record_to_gb_text(rec)
        sc._save_collections([
            {"name": "ProjectB", "plasmids": [
                {"id": "FFE1", "name": "FFE 1 ENTRY", "gb_text": gb,
                 "size": 80}]},
            {"name": "ProjectA", "plasmids": []},
        ])

    def test_sc_i_copy_plasmid(self):
        self._seed_collections()
        r = sc._h_copy_plasmid(None, {"name": "FFE 1 ENTRY",
                                       "to": "ProjectA", "from": "ProjectB"})
        assert r["ok"] and r["to"] == "ProjectA"
        coll = next(c for c in sc._load_collections()
                    if c["name"] == "ProjectA")
        copied = next(p for p in coll["plasmids"]
                      if p["name"] == "FFE 1 ENTRY")
        # Name + features (gb_text, incl. the SC-E display stamp) preserved.
        assert "SpliceCraft-name: FFE 1 ENTRY" in copied["gb_text"]
        # Duplicate → 409, missing source → 404, missing 'to' → 400.
        assert sc._h_copy_plasmid(None, {"name": "FFE 1 ENTRY",
                                          "to": "ProjectA"})[1] == 409
        assert sc._h_copy_plasmid(None, {"name": "Nope",
                                          "to": "ProjectA"})[1] == 404
        assert sc._h_copy_plasmid(None, {"name": "FFE 1 ENTRY"})[1] == 400


class TestPrimerCollectionReadsAndMove:
    """SC-J (collection is a real partition for READS) + SC-K (move-primer,
    no create/delete data-loss footgun)."""

    def _h(self, name):
        return sc._state._AGENT_HANDLERS[name][0]

    def _seed(self):
        # 5 in the default library, 3 filed into a named collection.
        # Distinct 2-mer suffixes so the dedup-by-sequence policy keeps all 5.
        for i, suf in enumerate(["AA", "CC", "GG", "TT", "AG"]):
            self._h("create-primer")(None, {
                "name": f"m{i}", "sequence": "ACGT" * 5 + suf})
        self._h("create-primer-collection")(None, {"name": "ProjA"})
        for i, suf in enumerate(["AA", "CC", "GG"]):
            self._h("create-primer")(None, {
                "name": f"p{i}", "sequence": "TTTT" * 5 + suf,
                "collection": "ProjA"})

    # ── SC-J: list-primers / get-primer honor {collection} ──────────────
    def test_list_primers_scoped_to_collection(self):
        self._seed()
        assert self._h("list-primers")(None, {"collection": "ProjA"})["count"] == 3
        assert self._h("list-primers")(None, {})["count"] == 5    # default
        assert self._h("list-primers")(None, {"collection": "Nope"})[1] == 404

    def test_get_primer_by_name_in_collection(self):
        self._seed()
        r = self._h("get-primer")(None, {"name": "p0", "collection": "ProjA"})
        assert r["primer"]["name"] == "p0"
        assert self._h("get-primer")(None, {"name": "p0"})[1] == 404  # not default
        assert self._h("get-primer")(None,
                                      {"name": "m0"})["primer"]["name"] == "m0"

    # ── SC-K: move-primer never loses the primer ────────────────────────
    def test_move_primer_no_data_loss(self):
        self._seed()
        r = self._h("move-primer")(None, {"name": "m0", "to": "ProjA"})
        assert r["ok"] and r["from"] == "" and r["to"] == "ProjA"
        assert self._h("list-primers")(None, {})["count"] == 4          # default −1
        assert self._h("list-primers")(None, {"collection": "ProjA"})["count"] == 4
        # The primer is in ProjA, not lost.
        assert self._h("get-primer")(None,
                                     {"name": "m0", "collection": "ProjA"})["primer"]["name"] == "m0"

    def test_move_primer_errors(self):
        self._seed()
        assert self._h("move-primer")(None, {"name": "m1", "to": "Nope"})[1] == 404
        assert self._h("move-primer")(None, {"to": "ProjA"})[1] == 400  # no selector
        # Already moved → not found in default.
        self._h("move-primer")(None, {"name": "m0", "to": "ProjA"})
        assert self._h("move-primer")(None,
                                      {"name": "m0", "to": "ProjA"})[1] in (404, 409)

    def test_move_to_active_collection_is_noop(self):
        # "" and the ACTIVE collection's name address the SAME store, so moving
        # a primer that's already in the active library "to" that collection is
        # a no-op → 409 (already there). (set-active re-mirrors the live
        # library to the active collection, so `a0` is what "" resolves to.)
        self._h("create-primer-collection")(None, {"name": "ProjA"})
        self._h("set-active-primer-collection")(None, {"name": "ProjA"})
        self._h("create-primer")(None, {"name": "a0",
                                         "sequence": "ACGT" * 5 + "AA"})
        r = self._h("move-primer")(None, {"name": "a0", "to": "ProjA"})
        assert isinstance(r, tuple) and r[1] == 409

    def test_move_out_of_active_collection(self):
        # THE FIX: with a named collection ACTIVE, primers.json holds its live
        # content (the mirror). Moving a primer OUT must route through
        # `_save_primers` so primers.json and the stored mirror stay
        # consistent — previously there was NO safe path ("" was refused while
        # a named collection was active, and passing the active name edited the
        # stored mirror only, which the next `_save_primers` reverted).
        self._h("create-primer-collection")(None, {"name": "GLOW"})
        self._h("create-primer-collection")(None, {"name": "SINK"})
        self._h("set-active-primer-collection")(None, {"name": "GLOW"})
        for i, suf in enumerate(["AA", "CC", "GG"]):
            self._h("create-primer")(None, {"name": f"g{i}",
                                             "sequence": "ACGT" * 5 + suf})
        # Live library (active GLOW) and GLOW's stored mirror agree at 3.
        assert self._h("list-primers")(None, {})["count"] == 3
        assert self._h("list-primers")(None,
                                       {"collection": "GLOW"})["count"] == 3
        # Move g0 OUT of the ACTIVE collection into SINK — previously blocked.
        r = self._h("move-primer")(None,
                                    {"name": "g0", "from": "GLOW", "to": "SINK"})
        assert r["ok"]
        assert self._h("list-primers")(None, {})["count"] == 2        # live −1
        assert self._h("list-primers")(None,
                                       {"collection": "SINK"})["count"] == 1
        # CONSISTENCY: primers.json (live) == GLOW's stored mirror (no desync).
        live = {p["name"] for p in self._h("list-primers")(None, {})["primers"]}
        stored = {p["name"] for p in
                  self._h("list-primers")(None, {"collection": "GLOW"})["primers"]}
        assert live == stored == {"g1", "g2"}

    def test_move_into_active_collection_by_name(self):
        # The inverse, addressing the active collection by its NAME: moving a
        # primer INTO it lands in BOTH primers.json and the stored mirror, and
        # leaves the source — also impossible before the fix.
        self._h("create-primer-collection")(None, {"name": "GLOW"})
        self._h("create-primer-collection")(None, {"name": "SRC"})
        self._h("set-active-primer-collection")(None, {"name": "GLOW"})
        self._h("create-primer")(None,
                                 {"name": "g0", "sequence": "ACGT" * 5 + "AA"})
        self._h("create-primer")(None, {"name": "s0",
                                        "sequence": "TTTT" * 5 + "CC",
                                        "collection": "SRC"})
        # Reference the active collection by its NAME ("GLOW") for `to`.
        r = self._h("move-primer")(None,
                                    {"name": "s0", "from": "SRC", "to": "GLOW"})
        assert r["ok"]
        assert self._h("list-primers")(None,
                                       {"collection": "SRC"})["count"] == 0
        live = {p["name"] for p in self._h("list-primers")(None, {})["primers"]}
        stored = {p["name"] for p in
                  self._h("list-primers")(None, {"collection": "GLOW"})["primers"]}
        assert "s0" in live and live == stored == {"g0", "s0"}

    def test_move_into_active_collection_survives_save_failure(self, monkeypatch):
        # HARDENING: moving INTO the active/live library persists the ADDITION
        # (primers.json + mirror) BEFORE the removal from the named source, so a
        # failed first save NEVER loses the primer — it stays in the source.
        # (The pre-hardening order removed from the source first, so a failed
        # second write dropped the primer entirely.)
        import splicecraft_agent as sca
        self._h("create-primer-collection")(None, {"name": "GLOW"})
        self._h("create-primer-collection")(None, {"name": "SRC"})
        self._h("set-active-primer-collection")(None, {"name": "GLOW"})
        self._h("create-primer")(None, {"name": "s0",
                                         "sequence": "TTTT" * 5 + "CC",
                                         "collection": "SRC"})

        def _boom(*a, **k):
            raise OSError("simulated disk failure")
        monkeypatch.setattr(sca, "_save_primers", _boom)

        r = self._h("move-primer")(None, {"name": "s0",
                                          "from": "SRC", "to": "GLOW"})
        assert isinstance(r, tuple) and r[1] == 500      # surfaced, not silent
        # NOT lost — still in SRC (the removal never ran).
        assert self._h("list-primers")(None,
                                       {"collection": "SRC"})["count"] == 1
        assert self._h("get-primer")(None, {"name": "s0", "collection": "SRC"}
                                     )["primer"]["name"] == "s0"

    def test_move_out_of_active_collection_survives_save_failure(self, monkeypatch):
        # HARDENING (other direction): moving OUT of the active library adds to
        # the named target FIRST, then removes from the live library. A failed
        # removal DUPLICATES the primer (in both stores), never loses it.
        import splicecraft_agent as sca
        self._h("create-primer-collection")(None, {"name": "GLOW"})
        self._h("create-primer-collection")(None, {"name": "SINK"})
        self._h("set-active-primer-collection")(None, {"name": "GLOW"})
        self._h("create-primer")(None, {"name": "g0",
                                         "sequence": "ACGT" * 5 + "AA"})

        def _boom(*a, **k):
            raise OSError("simulated disk failure")
        monkeypatch.setattr(sca, "_save_primers", _boom)

        r = self._h("move-primer")(None, {"name": "g0",
                                          "from": "GLOW", "to": "SINK"})
        assert isinstance(r, tuple) and r[1] == 500
        # DUPLICATED, not lost: the target got it (saved first) AND the live
        # library still has it (the removal save failed).
        assert self._h("get-primer")(None, {"name": "g0", "collection": "SINK"}
                                     )["primer"]["name"] == "g0"
        assert self._h("get-primer")(None, {"name": "g0"})["primer"]["name"] == "g0"

    # ── delete-primer {collection}: prune a collection without data loss ──
    def test_delete_primer_from_collection_keeps_dups_elsewhere(self):
        self._seed()
        # File a copy of a default primer (same sequence) INTO ProjA — the
        # "mirror pollution" shape — then prune it from ProjA only.
        self._h("create-primer")(None, {
            "name": "dup", "sequence": "ACGT" * 5 + "AA", "collection": "ProjA"})
        assert self._h("list-primers")(None, {"collection": "ProjA"})["count"] == 4
        r = self._h("delete-primer")(None, {
            "sequence": "ACGT" * 5 + "AA", "collection": "ProjA"})
        assert r["ok"] and r["removed"] == 1 and r["collection"] == "ProjA"
        # Gone from ProjA, but the original copy survives in the default lib.
        assert self._h("list-primers")(None, {"collection": "ProjA"})["count"] == 3
        assert self._h("get-primer")(None,
                                     {"sequence": "ACGT" * 5 + "AA"})["primer"]["name"] == "m0"

    def test_delete_primer_collection_guards(self):
        self._seed()
        assert self._h("delete-primer")(None,
                                        {"name": "x", "collection": "Nope"})[1] == 404
        assert self._h("delete-primer")(None,
                                        {"name": "absent", "collection": "ProjA"})[1] == 404
        # Active named collection → refuse (mirror safety).
        self._h("set-active-primer-collection")(None, {"name": "ProjA"})
        assert self._h("delete-primer")(None,
                                        {"name": "p0", "collection": "ProjA"})[1] == 409


class TestContainerManagement:
    """rename-/delete- for primer collections + parts bins — the container
    side of full-CRUD parity. Renaming/deleting the ACTIVE container must
    follow the active-pointer and swap the live mirror ([INV-83]); deleting
    the last container is refused."""

    def _h(self, name):
        return sc._state._AGENT_HANDLERS[name][0]

    # ── primer collections ──────────────────────────────────────────────
    def test_rename_primer_collection(self):
        for nm in ("Alpha", "Beta"):
            self._h("create-primer-collection")(None, {"name": nm})
        self._h("set-active-primer-collection")(None, {"name": "Alpha"})
        self._h("create-primer")(None,
                                 {"name": "p1", "sequence": "ACGT" * 5 + "GG"})
        # rename a NON-active bin (case-insensitive lookup).
        r = self._h("rename-primer-collection")(None,
                                                {"name": "beta", "new_name": "Gamma"})
        assert r["name"] == "Gamma" and r["renamed_from"] == "Beta"
        # collision + unknown.
        assert self._h("rename-primer-collection")(None,
                                                   {"name": "Alpha", "new_name": "Gamma"})[1] == 409
        assert self._h("rename-primer-collection")(None,
                                                   {"name": "nope", "new_name": "X"})[1] == 404
        # rename the ACTIVE collection → active-pointer follows, content kept.
        self._h("rename-primer-collection")(None,
                                            {"name": "Alpha", "new_name": "A2"})
        assert sc._get_active_primer_collection_name() == "A2"
        assert self._h("list-primers")(None, {"collection": "A2"})["count"] == 1

    def test_delete_primer_collection(self):
        for nm in ("Alpha", "Beta", "Delta"):
            self._h("create-primer-collection")(None, {"name": nm})
        self._h("set-active-primer-collection")(None, {"name": "Alpha"})
        # delete a NON-active collection — active pointer untouched.
        d = self._h("delete-primer-collection")(None, {"name": "Beta"})
        assert d["deleted"] == "Beta" and d["promoted"] == ""
        assert sc._get_active_primer_collection_name() == "Alpha"
        # delete the ACTIVE collection → promote first remaining + mirror-swap.
        d2 = self._h("delete-primer-collection")(None, {"name": "Alpha"})
        assert d2["deleted"] == "Alpha" and d2["promoted"] != ""
        assert sc._get_active_primer_collection_name() == d2["promoted"]
        # can't delete the last remaining named collection.
        last = sc._load_primer_collections()[0]["name"]
        assert self._h("delete-primer-collection")(None, {"name": last})[1] == 409
        assert self._h("delete-primer-collection")(None, {"name": "ghost"})[1] == 404

    # ── parts bins ──────────────────────────────────────────────────────
    def test_rename_parts_bin(self):
        for nm in ("BinA", "BinB"):
            self._h("create-parts-bin")(None, {"name": nm})
        self._h("set-active-parts-bin")(None, {"name": "BinA"})
        self._h("create-part")(None, {
            "name": "part1", "sequence": "ACGTACGTACGT", "bin": "BinA"})
        r = self._h("rename-parts-bin")(None,
                                        {"name": "binb", "new_name": "BinG"})
        assert r["name"] == "BinG" and r["renamed_from"] == "BinB"
        assert self._h("rename-parts-bin")(None,
                                           {"name": "BinA", "new_name": "BinG"})[1] == 409
        self._h("rename-parts-bin")(None, {"name": "BinA", "new_name": "BinA2"})
        assert sc._get_active_parts_bin_name() == "BinA2"
        assert self._h("list-parts")(None, {"bin": "BinA2"})["count"] == 1

    def test_delete_parts_bin(self):
        for nm in ("BinA", "BinB", "BinD"):
            self._h("create-parts-bin")(None, {"name": nm})
        self._h("set-active-parts-bin")(None, {"name": "BinA"})
        d = self._h("delete-parts-bin")(None, {"name": "BinB"})
        assert d["deleted"] == "BinB" and d["promoted"] == ""
        d2 = self._h("delete-parts-bin")(None, {"name": "BinA"})
        assert d2["promoted"] != ""
        assert sc._get_active_parts_bin_name() == d2["promoted"]
        last = sc._load_parts_bin_collections()[0]["name"]
        assert self._h("delete-parts-bin")(None, {"name": last})[1] == 409
        assert self._h("delete-parts-bin")(None, {"name": "ghost"})[1] == 404


class TestExperimentProjectCRUD:
    """SC-M-class: create/list/delete-experiment honor a ``{project}`` scope
    (a real partition, not silently filed into the active project), and
    move-experiment relocates an entry without a create/delete data-loss
    round-trip. Active-project writes stay mirror-consistent ([INV-83])."""

    def _h(self, name):
        return sc._state._AGENT_HANDLERS[name][0]

    def _seed(self):
        # Two named projects; ProjA active with one live entry.
        self._h("create-experiment-project")(None, {"name": "ProjA"})
        self._h("create-experiment-project")(None, {"name": "ProjB"})
        self._h("set-active-experiment-project")(None, {"name": "ProjA"})
        self._h("create-experiment")(None, {"title": "live A"})  # active ProjA

    # ── create-experiment {project} files into the named project ─────────
    def test_create_experiment_scoped_to_project(self):
        self._seed()
        r = self._h("create-experiment")(None,
                                         {"title": "filed B", "project": "ProjB"})
        assert r["ok"] and r["project"] == "ProjB"
        # Landed in ProjB, NOT the active ProjA.
        assert self._h("list-experiments")(None,
                                           {"project": "ProjB"})["experiments"][0]["title"] == "filed B"
        assert len(self._h("list-experiments")(None,
                                               {"project": "ProjA"})["experiments"]) == 1
        # Unknown project + bad type.
        assert self._h("create-experiment")(None,
                                            {"title": "x", "project": "Nope"})[1] == 404
        assert self._h("create-experiment")(None,
                                            {"title": "x", "project": 5})[1] == 400

    # ── list-experiments {project} is a real read partition ──────────────
    def test_list_experiments_scoped(self):
        self._seed()
        self._h("create-experiment")(None, {"title": "b1", "project": "ProjB"})
        assert len(self._h("list-experiments")(None,
                                              {"project": "ProjB"})["experiments"]) == 1
        assert len(self._h("list-experiments")(None, {})["experiments"]) == 1  # active
        assert self._h("list-experiments")(None, {"project": "Nope"})[1] == 404

    # ── delete-experiment {project} prunes the named project ─────────────
    def test_delete_experiment_scoped_and_active_guard(self):
        self._seed()
        r = self._h("create-experiment")(None,
                                         {"title": "doomed", "project": "ProjB"})
        eid = r["id"]
        d = self._h("delete-experiment")(None, {"id": eid, "project": "ProjB"})
        assert d["ok"] and d["project"] == "ProjB" and d["remaining"] == 0
        # Deleting from the ACTIVE project by name → 409 (switch away first).
        live = self._h("list-experiments")(None, {})["experiments"][0]["id"]
        assert self._h("delete-experiment")(None,
                                           {"id": live, "project": "ProjA"})[1] == 409
        # Unknown project / unknown id.
        assert self._h("delete-experiment")(None,
                                           {"id": eid, "project": "Nope"})[1] == 404
        assert self._h("delete-experiment")(None,
                                           {"id": "exp-ffffffff", "project": "ProjB"})[1] == 404

    # ── move-experiment never loses the entry ────────────────────────────
    def test_move_experiment_no_data_loss(self):
        self._seed()
        live = self._h("list-experiments")(None, {})["experiments"][0]["id"]
        r = self._h("move-experiment")(None,
                                      {"id": live, "to": "ProjB", "from": "ProjA"})
        assert r["ok"] and r["from"] == "ProjA" and r["to"] == "ProjB"
        # Gone from ProjA (active, mirror re-synced), now in ProjB.
        assert len(self._h("list-experiments")(None, {})["experiments"]) == 0
        assert len(self._h("list-experiments")(None,
                                              {"project": "ProjB"})["experiments"]) == 1

    def test_move_experiment_errors(self):
        self._seed()
        live = self._h("list-experiments")(None, {})["experiments"][0]["id"]
        assert self._h("move-experiment")(None,
                                         {"id": live, "to": "Nope"})[1] == 404   # bad dest
        assert self._h("move-experiment")(None, {"to": "ProjB"})[1] == 400        # no id
        assert self._h("move-experiment")(None, {"id": live})[1] == 400           # no dest
        # Auto-find from active, then same-project move → 409.
        assert self._h("move-experiment")(None,
                                         {"id": live, "to": "ProjA"})[1] == 409


class TestMovePlasmid:
    """The last container-CRUD parity gap: move-plasmid relocates a plasmid
    between collections atomically (no copy/delete data-loss round-trip).
    Unlike copy-plasmid it handles the ACTIVE collection on either side,
    re-staging the live plasmid_library.json mirror ([INV-83])."""

    def _ent(self, name, i):
        return {"id": f"id-{name}-{i}", "name": name, "gb_text": "", "size": 80}

    def _names(self, coll):
        c = next(c for c in sc._load_collections() if c["name"] == coll)
        return [p["name"] for p in c["plasmids"]]

    def _lib(self):
        return [p["name"] for p in sc._load_library()]

    def _seed(self):
        # CollC is the (neutral) active collection so non-active moves leave
        # the live mirror untouched; CollA/CollB hold the movable plasmids.
        sc._save_collections([
            {"name": "CollA",
             "plasmids": [self._ent("pA1", 1), self._ent("pA2", 2)],
             "saved": "2026-01-01"},
            {"name": "CollB", "plasmids": [self._ent("pB1", 1)],
             "saved": "2026-01-01"},
            {"name": "CollC", "plasmids": [], "saved": "2026-01-01"},
        ])
        sc._set_active_collection_name("CollC")
        sc._restore_library_from_active_collection()

    def test_move_between_non_active_collections(self):
        self._seed()
        r = sc._h_move_plasmid(None, {"name": "pA1", "to": "CollB"})
        assert r["ok"] and r["from"] == "CollA" and r["to"] == "CollB"
        assert self._names("CollA") == ["pA2"]
        assert "pA1" in self._names("CollB")
        assert self._lib() == []  # active CollC mirror untouched

    def test_move_into_and_out_of_active_remirrors(self):
        self._seed()
        # Into active CollC → live mirror gains the plasmid.
        sc._h_move_plasmid(None, {"name": "pB1", "to": "CollC", "from": "CollB"})
        assert "pB1" in self._names("CollC") and self._lib() == ["pB1"]
        # Out of active CollC → live mirror shrinks (L3-safe, no guard trip).
        sc._h_move_plasmid(None, {"name": "pB1", "to": "CollA"})
        assert self._lib() == [] and "pB1" in self._names("CollA")

    def test_move_plasmid_guards(self):
        self._seed()
        assert sc._h_move_plasmid(None, {"name": "ghost", "to": "CollA"})[1] == 404
        assert sc._h_move_plasmid(None, {"name": "pA2", "to": "Nope"})[1] == 404
        assert sc._h_move_plasmid(None, {"to": "CollA"})[1] == 400
        # already-in-target no-op.
        assert sc._h_move_plasmid(None,
                                  {"name": "pA2", "to": "CollA"})[1] == 409

    def test_move_plasmid_ambiguous_and_collision(self):
        sc._save_collections([
            {"name": "X", "plasmids": [self._ent("dup", 1)], "saved": "2026-01-01"},
            {"name": "Y", "plasmids": [self._ent("dup", 2)], "saved": "2026-01-01"},
        ])
        sc._set_active_collection_name("X")
        sc._restore_library_from_active_collection()
        # Same name in two collections → ambiguous without 'from'.
        assert sc._h_move_plasmid(None, {"name": "dup", "to": "Y"})[1] == 409
        # Name already present in target → refuse (no silent rename).
        assert sc._h_move_plasmid(None,
                                  {"name": "dup", "to": "Y", "from": "X"})[1] == 409


class TestDataEnvelope:
    """SC-C: the dispatcher adds a predictable `data` field to success
    responses (pure helper `_agent_data_envelope`)."""

    def test_single_content_key_unwraps_to_value(self):
        assert sc._agent_data_envelope(
            {"ok": True, "seq": "ACGT"}, 200)["data"] == "ACGT"
        assert sc._agent_data_envelope(
            {"ok": True, "library": [1, 2]}, 200)["data"] == [1, 2]

    def test_multi_content_key_is_dict(self):
        d = sc._agent_data_envelope({"ok": True, "a": 1, "b": 2}, 200)["data"]
        assert d == {"a": 1, "b": 2}

    def test_meta_excluded_originals_preserved(self):
        r = sc._agent_data_envelope(
            {"ok": True, "seq": "AC", "ignored": ["z"], "_stale": True}, 200)
        assert r["data"] == "AC"                       # meta stripped
        assert r["seq"] == "AC" and r["ignored"] == ["z"]   # non-breaking

    def test_errors_and_non2xx_untouched(self):
        assert "data" not in sc._agent_data_envelope({"error": "x"}, 400)
        assert "data" not in sc._agent_data_envelope(
            {"ok": False, "error": "y"}, 422)

    def test_existing_data_not_clobbered(self):
        assert sc._agent_data_envelope(
            {"ok": True, "data": "keep", "x": 1}, 200)["data"] == "keep"

    def test_empty_content_no_data(self):
        assert "data" not in sc._agent_data_envelope({"ok": True}, 200)

    def test_ok_flag_added_when_missing(self):
        # F1/F2: a 2xx success payload without `ok` gets a uniform ok:true, and
        # the data-unwrap shape is unchanged (`ok` is an envelope meta key).
        r = sc._agent_data_envelope({"gels": [1, 2]}, 200)
        assert r["ok"] is True
        assert r["gels"] == [1, 2]        # container key untouched
        assert r["data"] == [1, 2]        # still unwrapped to the bare list

    def test_ok_flag_not_clobbered(self):
        # An explicit `ok` is never overwritten.
        assert sc._agent_data_envelope({"ok": False, "x": 1}, 200)["ok"] is False

    def test_ok_flag_not_added_to_errors_or_non2xx(self):
        assert "ok" not in sc._agent_data_envelope({"error": "x"}, 400)
        assert "ok" not in sc._agent_data_envelope({"matches": []}, 500)


class TestRestartEndpoint:
    """SC-A: `restart` re-execs a HEADLESS daemon to pick up an update;
    refuses in a GUI session (would kill the live TUI)."""

    def test_non_headless_refused_409(self):
        called = []
        orig = sc._agent_schedule_restart
        sc._agent_schedule_restart = lambda: called.append(True)
        try:
            r = sc._h_restart(types.SimpleNamespace(_headless=False), {})
            assert isinstance(r, tuple) and r[1] == 409
            assert not called                # did NOT schedule a re-exec
        finally:
            sc._agent_schedule_restart = orig

    def test_headless_schedules_restart(self):
        called = []
        orig = sc._agent_schedule_restart
        sc._agent_schedule_restart = lambda: called.append(True)
        try:
            r = sc._h_restart(types.SimpleNamespace(_headless=True), {})
            assert r["ok"] and r["restarting"] is True
            assert r["running_version"] == sc.__version__
            assert called == [True]
        finally:
            sc._agent_schedule_restart = orig


class TestShutdownEndpoint:
    """`shutdown` stops a HEADLESS daemon through the normal teardown, so
    the flush-on-exit path runs and the data-dir lock is released — the
    thing a SIGKILL skips. Refuses in a GUI session (it has its own quit
    flow) and on unsaved edits (quitting is not undoable)."""

    class _FakeApp:
        def __init__(self, *, headless=True, unsaved=False, running=True):
            self._headless = headless
            self._unsaved = unsaved
            self.is_running = running
            self.exited = False
            self.marshalled = []

        def exit(self):
            self.exited = True

        def call_from_thread(self, fn, *a, **kw):
            self.marshalled.append(fn)
            return fn(*a, **kw)

    def _stub_schedule(self, sink):
        orig = sc._agent_schedule_shutdown
        sc._agent_schedule_shutdown = lambda app: sink.append(app)
        return orig

    def test_non_headless_refused_409(self):
        called = []
        orig = self._stub_schedule(called)
        try:
            r = sc._h_shutdown(self._FakeApp(headless=False), {})
            assert isinstance(r, tuple) and r[1] == 409
            assert not called            # did NOT schedule a quit
        finally:
            sc._agent_schedule_shutdown = orig

    def test_headless_schedules_shutdown(self):
        called = []
        orig = self._stub_schedule(called)
        try:
            r = sc._h_shutdown(self._FakeApp(), {})
            assert r["ok"] and r["shutting_down"] is True
            # The caller polls this pid to confirm the process is gone.
            assert r["pid"] == os.getpid()
            assert len(called) == 1
        finally:
            sc._agent_schedule_shutdown = orig

    def test_unsaved_edits_refused_409(self):
        called = []
        orig = self._stub_schedule(called)
        try:
            r = sc._h_shutdown(self._FakeApp(unsaved=True), {})
            assert isinstance(r, tuple) and r[1] == 409
            assert r[0]["dirty"] is True
            assert not called            # unsaved work is never dropped silently
        finally:
            sc._agent_schedule_shutdown = orig

    def test_force_overrides_dirty(self):
        called = []
        orig = self._stub_schedule(called)
        try:
            r = sc._h_shutdown(self._FakeApp(unsaved=True), {"force": True})
            assert r["ok"] and r["shutting_down"] is True
            assert len(called) == 1
        finally:
            sc._agent_schedule_shutdown = orig

    def test_exec_marshals_app_exit(self, monkeypatch):
        """Must go through `App.exit()` on the app thread — NOT
        `action_quit()`, which pushes a modal a headless daemon can
        never answer."""
        monkeypatch.setattr(time, "sleep", lambda _s: None)
        app = self._FakeApp()
        sc._agent_shutdown_exec(app)
        assert app.exited is True
        assert app.marshalled == [app.exit]

    def test_exec_noop_when_loop_stopped(self, monkeypatch):
        """A stopped loop must not be marshalled to — `call_from_thread`
        would queue a coroutine nobody awaits."""
        monkeypatch.setattr(time, "sleep", lambda _s: None)
        app = self._FakeApp(running=False)
        sc._agent_shutdown_exec(app)
        assert app.exited is False
        assert app.marshalled == []

    def test_exec_never_force_exits_on_failure(self, monkeypatch):
        """A wedged loop is logged and left alive: `os._exit` here would
        skip the very flushes that make this shutdown graceful."""
        monkeypatch.setattr(time, "sleep", lambda _s: None)

        def _boom(*_a, **_kw):
            raise RuntimeError("loop wedged")

        app = self._FakeApp()
        monkeypatch.setattr(app, "call_from_thread", _boom)
        hard_exit = []
        monkeypatch.setattr(os, "_exit", lambda code: hard_exit.append(code))
        sc._agent_shutdown_exec(app)          # must not raise
        assert hard_exit == []

    def test_registered_as_write_endpoint(self):
        fn, write = sc._state._AGENT_HANDLERS["shutdown"]
        assert fn is sc._h_shutdown
        assert write is True                  # bearer token required


class TestGoldenGate:
    """SC-H: Type IIS (BsaI) Golden Gate / MoClo assembly — overhang-directed,
    order-independent, a real digest + ligation."""

    @staticmethod
    def _cassette(oh5, body, oh3):
        # The canonical L0 part: BsaI sites release a body with (oh5, oh3).
        return "GGTCTC" + "A" + oh5 + body + oh3 + "A" + "GAGACC"

    def _design(self):
        A = self._cassette("GGAG", "ATGAAACCCGGGTTTACGT" * 2, "AATG")
        B = self._cassette("AATG", "TTGCATGCATGCTAGCTAG" * 2, "CGCT")
        V = self._cassette("CGCT", "GGGGCCCCAAAATTTT" * 8, "GGAG")
        return A, B, V

    def test_simulate_assembles_circle(self):
        A, B, V = self._design()
        r = sc._h_simulate_golden_gate(None, {"parts": [A, B], "vector": V})
        assert r["ok"] and r["result"]["ok"]
        assert r["result"]["circular"] and r["result"]["n_parts"] == 2
        assert r["result"]["n_residual_sites"] == 0
        assert r["result"]["warnings"] == []

    def test_order_independent(self):
        A, B, V = self._design()
        r1 = sc._h_simulate_golden_gate(None, {"parts": [A, B], "vector": V})
        r2 = sc._h_simulate_golden_gate(None, {"parts": [B, A], "vector": V})
        assert r1["result"]["length"] == r2["result"]["length"]

    def test_assemble_saves_circular_product(self):
        A, B, V = self._design()
        r = sc._h_golden_gate_assemble(
            None, {"parts": [A, B], "vector": V, "product_name": "gg1"})
        assert r["ok"] and r["saved_name"] == "gg1"
        ent = next(e for e in sc._load_library() if e["name"] == "gg1")
        assert ent["source"] == "agent:golden-gate"
        rec = sc._gb_text_to_record(ent["gb_text"])
        assert rec.annotations.get("topology") == "circular"
        seq = str(rec.seq).upper()
        for body in ("ATGAAACCCGGGTTTACGT", "TTGCATGCATGCTAGCTAG",
                     "GGGGCCCCAAAATTTT"):
            assert body in seq or body in sc._rc(seq), f"{body} missing"

    def test_non_type_iis_rejected(self):
        A, B, V = self._design()
        assert sc._h_golden_gate_assemble(
            None, {"parts": [A, B], "vector": V, "enzyme": "EcoRI",
                   "product_name": "x"})[1] == 422

    def test_empty_parts_400(self):
        _, _, V = self._design()
        assert sc._h_simulate_golden_gate(
            None, {"parts": [], "vector": V})[1] == 400

    def test_part_with_internal_site_fails(self):
        A, B, V = self._design()
        bad = self._cassette("GGAG", "ATG" + "GGTCTC" + "ATTTT", "AATG")
        r = sc._h_simulate_golden_gate(None, {"parts": [bad, B], "vector": V})
        assert r["result"]["ok"] is False and r["result"]["errors"]


class TestGoldenGateMultiCutVector:
    """2026-08-26 report #2: a destination vector carrying BACKGROUND enzyme
    sites is cut into more than two pieces, and the simulator assumed exactly
    one backbone fragment. It answered "the overhangs don't chain… check that
    adjacent parts share a 4 nt overhang" — pointing at parts that were fine,
    which is where the reporter lost the most time in a day's work."""

    HALF = "GGGGCCCCAAAATTTT" * 4

    @staticmethod
    def _cassette(oh5, body, oh3):
        return "GGTCTC" + "A" + oh5 + body + oh3 + "A" + "GAGACC"

    def _parts(self):
        return (self._cassette("GGAG", "ATGAAACCCGGGTTTACGT" * 2, "AATG"),
                self._cassette("AATG", "TTGCATGCATGCTAGCTAG" * 2, "CGCT"))

    def _plain_vector(self):
        return self._cassette("CGCT", "GGGGCCCCAAAATTTT" * 8, "GGAG")

    def _split_vector(self):
        """Same backbone, one extra BsaI site inside it — so the backbone is
        released as TWO pieces joined by a new TTAC overhang."""
        return self._cassette(
            "CGCT", self.HALF + "GGTCTC" + "T" + "TTAC" + self.HALF, "GGAG")

    def test_a_three_piece_vector_still_assembles(self):
        A, B = self._parts()
        r = sc._h_simulate_golden_gate(
            None, {"parts": [A, B], "vector": self._split_vector()})
        assert r["ok"] is True, r["result"]["errors"]
        res = r["result"]
        assert res["n_vector_fragments"] == 3
        assert res["n_vector_fragments_used"] == 2
        assert res["n_vector_sites"] == 3

    def test_the_product_keeps_every_piece_of_the_backbone(self):
        """Not just "it closed" — the extra vector fragment has to BE in the
        product, or the sim quietly built a different plasmid."""
        A, B = self._parts()
        r = sc._h_simulate_golden_gate(
            None, {"parts": [A, B], "vector": self._split_vector()})
        prod = r["result"]["product_seq"]
        doubled = prod + prod              # circular: allow the join to wrap
        for body in ("ATGAAACCCGGGTTTACGT", "TTGCATGCATGCTAGCTAG", self.HALF):
            assert body in doubled, body
        # Both halves of the split backbone, not one of them twice.
        assert doubled.count(self.HALF) >= 2

    def test_it_equals_the_undivided_vector_plus_the_extra_bases(self):
        """The background site adds `GGTCTC`+`T`+`TTAC` (11 bp) to the
        backbone and nothing else — an independent check on the product size
        that doesn't just re-run the assembler."""
        A, B = self._parts()
        plain = sc._h_simulate_golden_gate(
            None, {"parts": [A, B], "vector": self._plain_vector()})
        split = sc._h_simulate_golden_gate(
            None, {"parts": [A, B], "vector": self._split_vector()})
        assert split["result"]["length"] - plain["result"]["length"] == 11

    def test_a_two_piece_vector_reports_one_fragment_used(self):
        A, B = self._parts()
        res = sc._h_simulate_golden_gate(
            None, {"parts": [A, B], "vector": self._plain_vector()})["result"]
        assert res["n_vector_fragments"] == 2
        assert res["n_vector_fragments_used"] == 1      # the stuffer drops out

    def test_failure_shows_the_vector_cuts_and_every_overhang(self):
        """The fix the report asked for: say what was cut and what the ends
        ARE, instead of asserting the parts are at fault."""
        A, _ = self._parts()
        bad = self._cassette("AATG", "TTGCATGCATGCTAGCTAG" * 2, "TTTT")
        r = sc._h_simulate_golden_gate(
            None, {"parts": [A, bad], "vector": self._split_vector()})
        res = r["result"]
        assert r["ok"] is False
        assert res["n_vector_sites"] == 3 and res["n_vector_fragments"] == 3
        assert len(res["vector_cut_bp"]) == 3
        err = res["errors"][0]
        assert "cut 3×" in err
        assert "GGAG→AATG" in err and "AATG→TTTT" in err   # the dangling end
        assert res["part_overhangs"] == ["GGAG→AATG", "AATG→TTTT"]

    def test_ambiguous_overhangs_are_reported_not_silently_picked(self):
        """Two parts that can swap places make TWO different circles. Picking
        one and saying nothing is the wrong-answer shape this batch is
        about."""
        A = self._cassette("GGAG", "AAACCCGGGTTTACGT" * 2, "AATG")
        B = self._cassette("AATG", "TTGCATGCATGCTAGCT" * 2, "AATG")
        C = self._cassette("AATG", "GGTTAACCGGTTAACC" * 2, "AATG")
        D = self._cassette("AATG", "TTAACCGGTTAACCGG" * 2, "CGCT")
        r = sc._h_simulate_golden_gate(
            None, {"parts": [A, B, C, D], "vector": self._plain_vector()})
        assert any("different circles" in w
                   for w in r["result"]["warnings"]), r["result"]["warnings"]

    def test_a_synonym_names_the_assembly_enzyme(self):
        """The engines key on catalog names, so `Esp3I` (or `esp3i`) has to be
        resolved here or the reaction is simulated with an enzyme that cuts
        nothing — a dry run that answers about the wrong chemistry."""
        def cas(a, b, c):
            return "CGTCTC" + "A" + a + b + c + "A" + "GAGACG"   # BsmBI/Esp3I
        A = cas("GGAG", "ATGAAACCCTTTACGT" * 2, "AATG")
        B = cas("AATG", "TTGCATGCATGCTAGC" * 2, "CGCT")
        V = cas("CGCT", "GGGGCCCCAAAATTTT" * 8, "GGAG")
        lengths = set()
        for spelling in ("BsmBI", "Esp3I", "esp3i"):
            r = sc._h_simulate_golden_gate(
                None, {"parts": [A, B], "vector": V, "enzyme": spelling})
            assert r["ok"] is True, (spelling, r["result"]["errors"])
            lengths.add(r["result"]["length"])
        assert len(lengths) == 1

    def test_an_explosive_overhang_set_is_bounded_and_says_unknown(self):
        """Many fragments sharing one overhang make the number of possible
        orders factorial. The search is capped — and an exhausted search must
        NOT be reported as "this design can't be built", because it isn't the
        same answer."""
        parts = [self._cassette("AAAA", "ACGTTGCAACGTTGCA" * 2 + "ACGT" * (i % 4),
                                "AAAA") for i in range(20)]
        parts.append(self._cassette("AAAA", "TTGGCCAATTGGCCAA" * 2, "TTTA"))
        vec = self._cassette("CCCC", "GGGGCCCCAAAATTTT" * 8, "AAAA")
        r = sc._h_simulate_golden_gate(None, {"parts": parts, "vector": vec})
        assert r["ok"] is False
        err = r["result"]["errors"][0]
        assert "step limit" in err
        assert "UNKNOWN" in err

    def test_a_capped_ambiguity_count_says_at_least(self):
        """The enumeration stops at `_GG_MAX_SOLUTIONS` per seed, so the tally
        is a floor. Reporting a cap as an exact number is the quiet kind of
        inaccuracy this batch exists to remove."""
        import splicecraft_cloning as scl
        # A blunt-cutting Type IIS: every released end is blunt, so every
        # fragment ligates to every other and the orderings explode.
        sc._save_custom_enzymes([{"name": "BluntIISx", "site": "GGTCTC",
                                  "fwd_cut": 9, "rev_cut": 9}])
        sc._rebuild_scan_catalog()
        try:
            def cas(body):
                return "GGTCTC" + "AAA" + body + "TTT" + "GAGACC"
            parts = [cas("ACGTTGCAACGTTGCA" * 2 + "ACGT" * (i % 3))
                     for i in range(12)]
            vec = cas("GGGGCCCCAAAATTTT" * 8)
            r = sc._h_simulate_golden_gate(
                None, {"parts": parts, "vector": vec, "enzyme": "BluntIISx"})
            warns = r["result"]["warnings"]
            assert any("At least" in w and "different circles" in w
                       for w in warns), warns
        finally:
            sc._save_custom_enzymes([])
            sc._rebuild_scan_catalog()

    def test_the_assemble_endpoint_saves_the_multi_piece_product(self):
        A, B = self._parts()
        r = sc._h_golden_gate_assemble(
            None, {"parts": [A, B], "vector": self._split_vector(),
                   "product_name": "gg-split"})
        assert r["ok"] and r["saved_name"] == "gg-split"
        ent = next(e for e in sc._load_library() if e["name"] == "gg-split")
        rec = sc._gb_text_to_record(ent["gb_text"])
        assert rec.annotations.get("topology") == "circular"
        assert self.HALF in (str(rec.seq).upper() * 2)


class TestReplaceSequenceSizeCap:
    """Sweep #32 adversarial audit: `_h_replace_sequence` used
    to be capped only on the input `bases` field (1 MB via
    `_sanitize_bases`). A 100 MB pre-loaded record could be
    extended by ~1 MB on every call, growing unbounded. Result
    sequence cap (50 MB) blocks the bloat path."""

    def test_replace_refuses_when_result_too_large(self, tiny_app):
        # Build a record near the 50 MB cap. The mock app's
        # `_seq_len` returns `len(rec.seq)`; rec.seq is a Bio
        # Seq. Patch it to look 50 MB without allocating.
        from unittest.mock import MagicMock
        original_seq = tiny_app._current_record.seq
        big_n = 50 * 1024 * 1024
        tiny_app._current_record = MagicMock(
            seq=MagicMock(__len__=lambda self: big_n),
        )
        # Asking to append 1 KB at the end → final = 50 MB + 1 KB
        result = sc._h_replace_sequence(tiny_app, {
            "start": big_n,
            "end":   big_n,
            "bases": "A" * 1024,
            "force": True,
        })
        # Restore for the rest of the test suite.
        tiny_app._current_record.seq = original_seq
        # Expect 413 Payload Too Large with the cap details.
        assert isinstance(result, tuple), result
        payload, status = result
        assert status == 413, (status, payload)
        assert "result sequence too large" in payload["error"]
        assert payload["limit_bp"] == 50 * 1024 * 1024
        assert payload["final_bp"] > 50 * 1024 * 1024


class TestSaveHandler:
    def test_refuses_when_no_record(self):
        app = MockApp(record=None)
        result = sc._h_save(app, {})
        payload, status = result
        assert status == 422
        assert "nothing to save" in payload["error"]

    def test_calls_do_save(self, tiny_app):
        # The tiny_app record isn't in the library and has no source
        # file, so saving CREATES a new entry — opt in with create:true
        # (the homeless-new-record save gate, snag #16).
        result = sc._h_save(tiny_app, {"create": True})
        assert result["ok"] is True
        assert result["created"] is True
        assert tiny_app._saved is True

    def test_save_gate_blocks_homeless_create(self, tiny_app):
        # Without create:true a record with no source file that isn't in
        # the library refuses to silently create an entry (snag #16 —
        # the fetch->save pollution guard).
        result = sc._h_save(tiny_app, {})
        assert isinstance(result, tuple) and result[1] == 409
        assert result[0].get("would_create") is True
        assert tiny_app._saved is False

    def test_save_echoes_ignored_keys(self, tiny_app):
        # Unknown body keys are surfaced, not silently dropped (snag #14).
        result = sc._h_save(tiny_app, {"create": True, "new_name": "x"})
        assert result["ok"] is True
        assert result["ignored"] == ["new_name"]


class TestAgentApiAuditFixes:
    """Regression coverage for the agent-API audit pass (the agent-driven
    real-world build snags): rename-plasmid, schema discovery via
    `doc_full`, the active-pointer getters, and the ignored-key echo.
    The full HTTP round-trips (cross-collection load-entry name
    preservation, rename, save-gate) are exercised by the verifier."""

    def test_rename_plasmid_requires_old_and_new(self):
        app = MockApp(record=None)
        assert sc._h_rename_plasmid(app, {"new": "x"})[1] == 400
        assert sc._h_rename_plasmid(app, {"old": "x"})[1] == 400
        assert sc._h_rename_plasmid(app, {"old": "x", "new": "  "})[1] == 400

    def test_load_entry_requires_key(self):
        app = MockApp(record=None)
        assert sc._h_load_entry(app, {})[1] == 400

    def test_get_active_parity_endpoints_registered(self):
        H = sc._state._AGENT_HANDLERS
        for ep in ("get-active-collection", "get-active-codon-table",
                   "get-active-primer-collection", "get-active-parts-bin",
                   "get-active-experiment-project",
                   "get-active-hmm-database",
                   "get-active-enzyme-collection", "rename-plasmid"):
            assert ep in H, f"missing endpoint {ep}"

    def test_tools_emits_full_docstring(self):
        eps = sc._h_tools(None, {})["endpoints"]
        assert eps and all("doc_full" in e for e in eps)
        rp = next(e for e in eps if e["name"] == "rename-plasmid")
        # doc_full carries the body schema, not just the one-line summary.
        assert "old" in rp["doc_full"] and "new" in rp["doc_full"]
        assert len(rp["doc_full"]) > len(rp["doc"])

    def test_ignored_keys_helper(self):
        assert sc._agent_ignored_keys({"force": 1, "a": 2}, set()) == ["a"]
        assert sc._agent_ignored_keys({"force": 1}, set()) == []
        assert sc._agent_ignored_keys("not-a-dict", {"x"}) == []


class TestFeaturesHandler:
    def test_empty_when_no_record(self):
        app = MockApp(record=None)
        assert sc._h_features(app, {})["features"] == []

    def test_lists_feature_dicts(self, tiny_app):
        feats = sc._h_features(tiny_app, {})["features"]
        assert len(feats) >= 1
        assert all("idx" in f and "start" in f and "end" in f
                    for f in feats)


class TestFoldRnaHandler:
    """`fold-rna` — pure-Python RNA MFE folding via the agent API. The
    handler is record-independent (folds the payload sequence), so `app`
    is unused."""

    def test_folds_stemloop(self):
        r = sc._h_fold_rna(None, {"sequence": "GGGGAAAACCCC"})
        assert r["ok"] is True
        assert r["structure"] == "((((....))))"
        assert abs(r["dg"] - (-5.40)) < 0.011
        assert r["length"] == 12

    def test_seq_alias_and_dna_t(self):
        r = sc._h_fold_rna(None, {"seq": "GGGGTTTTCCCC"})    # 'seq' alias + T->U
        assert r["structure"] == "((((....))))"

    def test_missing_sequence_400(self):
        assert sc._h_fold_rna(None, {})[1] == 400

    def test_non_string_400(self):
        assert sc._h_fold_rna(None, {"sequence": 123})[1] == 400

    def test_ambiguous_400(self):
        body, code = sc._h_fold_rna(None, {"sequence": "ACGUN"})
        assert code == 400 and "error" in body

    def test_overlength_400(self):
        assert sc._h_fold_rna(None, {"sequence": "A" * 700})[1] == 400

    def test_registered_read_only(self):
        assert "fold-rna" in sc._AGENT_HANDLERS
        _fn, write = sc._AGENT_HANDLERS["fold-rna"]
        assert write is False


class TestCofoldRnaHandler:
    """`cofold-rna` — bound-state heterodimer ΔG (anti-SD : mRNA hybrid)."""

    def test_antisd_duplex(self):
        r = sc._h_cofold_rna(None, {"seq_a": "UAAGGAGGU", "seq_b": "ACCUCCUUA"})
        assert r["ok"] is True and r["dg"] < -10.0

    def test_dna_t_accepted(self):
        assert sc._h_cofold_rna(None, {"seq_a": "GGGGTTTT", "seq_b": "AAAACCCC"})["ok"]

    def test_missing_400(self):
        assert sc._h_cofold_rna(None, {"seq_a": "ACGU"})[1] == 400

    def test_ambiguous_400(self):
        assert sc._h_cofold_rna(None, {"seq_a": "ACGUN", "seq_b": "ACGU"})[1] == 400


class TestRbsStrengthHandler:
    """`rbs-strength` — relative translation-initiation strength."""

    def test_strong_rbs(self):
        r = sc._h_rbs_strength(
            None, {"mrna": "AAUAAAAGGAGGAAUAAAUGAGCAAAGCAACU", "start": 17})
        assert r["ok"] is True and r["rel_strength"] > 1 and "spacing" in r

    def test_missing_mrna_400(self):
        assert sc._h_rbs_strength(None, {"start": 5})[1] == 400

    def test_missing_start_400(self):
        assert sc._h_rbs_strength(None, {"mrna": "AUGAAAAAA"})[1] == 400

    def test_bad_start_400(self):
        assert sc._h_rbs_strength(None, {"mrna": "AUGAAA", "start": 99})[1] == 400

    def test_bool_start_400(self):
        # JSON true must not silently coerce to start=1 (audit #5e)
        assert sc._h_rbs_strength(None, {"mrna": "AUGAAAAAA", "start": True})[1] == 400

    def test_result_is_json_serializable(self):
        # start too close to the 5' end used to return dg_total=Infinity,
        # which is invalid JSON shipped with HTTP 200 (audit #1)
        import json as _json
        r = sc._h_rbs_strength(None, {"mrna": "AAAAAAAUGAAA", "start": 5})
        assert r["ok"] is True and r["dg_total"] is None
        body = _json.dumps(r)
        assert "Infinity" not in body and "NaN" not in body


class TestDesignRbsHandler:
    """`design-rbs` — reverse-design a 5'UTR for a target strength."""

    def test_designs_to_target(self):
        r = sc._h_design_rbs(None, {"cds": "AUGAGCAAAUACUAA", "target": 5.0})
        assert r["ok"] is True and r["full"].endswith("AUGAGCAAAUACUAA")
        assert "on_target" in r and isinstance(r["spacing"], int)

    def test_missing_cds_400(self):
        assert sc._h_design_rbs(None, {"target": 5})[1] == 400

    def test_bad_target_400(self):
        assert sc._h_design_rbs(None, {"cds": "AUGAAAUAA", "target": -1})[1] == 400
        assert sc._h_design_rbs(None, {"cds": "AUGAAAUAA", "target": "x"})[1] == 400

    def test_endpoints_registered_read_only(self):
        for name in ("cofold-rna", "rbs-strength", "design-rbs"):
            assert name in sc._AGENT_HANDLERS
            assert sc._AGENT_HANDLERS[name][1] is False


class TestAssembleOperonHandler:
    """`assemble-operon` — context-aware operon assembly from CDSs."""

    GENES = [{"cds": "AUGAGCAAAUACUAA", "target_strength": 5.0, "name": "A"},
             {"cds": "AUGGCAGAAUGGUAA", "target": 2.0, "name": "B"}]

    def test_assembles(self):
        r = sc._h_assemble_operon(None, {"genes": self.GENES,
                                         "promoter": "TTGACA",
                                         "terminator": "TTTT"})
        assert r["ok"] is True
        assert "U" not in r["sequence"]
        assert [el["kind"] for el in r["layout"]] == \
            ["promoter", "rbs", "cds", "rbs", "cds", "terminator"]
        assert len(r["genes"]) == 2 and "on_target" in r["genes"][0]

    def test_missing_genes_400(self):
        assert sc._h_assemble_operon(None, {})[1] == 400
        assert sc._h_assemble_operon(None, {"genes": []})[1] == 400

    def test_bad_gene_400(self):
        bad = {"genes": [{"cds": "AU", "target": 5}]}
        assert sc._h_assemble_operon(None, bad)[1] == 400

    def test_nonfinite_target_400(self):
        # NaN / Infinity target used to be echoed -> invalid JSON @ 200 (audit #2)
        for bad in (float("inf"), float("nan")):
            r = sc._h_assemble_operon(
                None, {"genes": [{"cds": "AUGAAAUACUAA", "target_strength": bad}]})
            assert r[1] == 400

    def test_registered_read_only(self):
        assert "assemble-operon" in sc._AGENT_HANDLERS
        assert sc._AGENT_HANDLERS["assemble-operon"][1] is False


class TestExportGffHandler:
    """`_h_export_gff` writes the loaded record to disk as GFF3."""

    def test_no_record_returns_422(self):
        app = MockApp(record=None)
        result = sc._h_export_gff(app, {"path": "/tmp/x.gff3"})
        assert result == ({"error": "no plasmid loaded"}, 422)

    def test_missing_path_returns_400(self, tiny_record):
        app = MockApp(record=tiny_record)
        result = sc._h_export_gff(app, {})
        assert isinstance(result, tuple) and result[1] == 400

    def test_writes_file(self, tiny_record, tmp_path):
        app = MockApp(record=tiny_record)
        out = tmp_path / "tiny.gff3"
        result = sc._h_export_gff(app, {"path": str(out)})
        assert isinstance(result, dict)
        assert result["ok"] is True
        assert out.exists()
        assert out.read_text().startswith("##gff-version 3")


class TestTransferAnnotationsHandler:
    """`_h_transfer_annotations` walks a source library entry's
    features and matches them onto the loaded record by sequence
    identity. Defaults to dry-run so an agent can inspect the
    proposed transfers before committing."""

    def test_no_record_loaded_returns_422(self):
        app = MockApp(record=None)
        result = sc._h_transfer_annotations(app, {"source_id": "x"})
        assert result == ({"error": "no plasmid loaded"}, 422)

    def test_missing_source_id_returns_400(self, tiny_record):
        app = MockApp(record=tiny_record)
        result = sc._h_transfer_annotations(app, {})
        assert isinstance(result, tuple) and result[1] == 400

    def test_source_not_in_library_returns_404(self, tiny_record,
                                                  isolated_library):
        sc._save_library([])
        app = MockApp(record=tiny_record)
        result = sc._h_transfer_annotations(
            app, {"source_id": "ghost"},
        )
        assert isinstance(result, tuple) and result[1] == 404

    def test_dry_run_returns_transfers_without_applying(
        self, tiny_record, isolated_library
    ):
        # Source library entry mirrors the loaded record so every
        # feature finds itself.
        sc._save_library([{
            "id":      "src",
            "name":    "src",
            "size":    len(tiny_record.seq),
            "n_feats": len(tiny_record.features),
            "added":   "2026-05-06",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = MockApp(record=tiny_record)
        before = len(tiny_record.features)
        result = sc._h_transfer_annotations(
            app, {"source_id": "src", "dry_run": True},
        )
        assert isinstance(result, dict)
        assert result["applied"] is False
        # Transfer count: only the >= min_len features. tiny_record
        # has a 27-bp CDS and a 30-bp misc_feature; min_len defaults
        # to 30 so only the misc_feature qualifies. Whatever the
        # exact count, the record itself must NOT have been mutated.
        assert len(tiny_record.features) == before

    def test_resolves_source_by_display_name(
        self, tiny_record, isolated_library
    ):
        """P1-11 (agent-API feedback): `source_id` resolves by display NAME,
        not only the LOCUS-safe id. Every other read resolves by name OR
        id; pre-fix a valid display name 404'd here."""
        sc._save_library([{
            "id":      "SRC_LOCUS_ID",
            "name":    "Src Display Name",
            "size":    len(tiny_record.seq),
            "n_feats": len(tiny_record.features),
            "added":   "2026-06-24",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = MockApp(record=tiny_record)
        # By display name (spaces, not the id) — must resolve, not 404.
        by_name = sc._h_transfer_annotations(
            app, {"source_id": "Src Display Name", "dry_run": True},
        )
        assert isinstance(by_name, dict), by_name
        assert by_name["applied"] is False
        # The canonical id still works too.
        by_id = sc._h_transfer_annotations(
            app, {"source_id": "SRC_LOCUS_ID", "dry_run": True},
        )
        assert isinstance(by_id, dict), by_id

    # ── Apply switch (agent-API feedback) ──────────────────────────────
    def _seed_self_source(self, tiny_record):
        """A source library entry that mirrors the loaded record, so every
        feature finds itself (≥ 1 transfer at the default min_len)."""
        sc._save_library([{
            "id": "src", "name": "src", "size": len(tiny_record.seq),
            "n_feats": len(tiny_record.features), "added": "2026-06-24",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])

    def test_apply_true_applies(self, tiny_record, isolated_library):
        # `apply:true` is the discoverable affirmative — pre-fix it was
        # silently ignored and the call stayed a dry run.
        self._seed_self_source(tiny_record)
        app = MockApp(record=tiny_record)
        before = len(app._current_record.features)
        r = sc._h_transfer_annotations(app, {"source_id": "src", "apply": True})
        assert isinstance(r, dict), r
        assert r["applied"] is True and r["count"] >= 1
        assert len(app._current_record.features) > before

    def test_dry_run_false_still_applies(self, tiny_record, isolated_library):
        # Backward compat: the legacy dry_run:false still applies.
        self._seed_self_source(tiny_record)
        app = MockApp(record=tiny_record)
        r = sc._h_transfer_annotations(
            app, {"source_id": "src", "dry_run": False})
        assert isinstance(r, dict) and r["applied"] is True, r

    def test_apply_false_is_dry_run(self, tiny_record, isolated_library):
        self._seed_self_source(tiny_record)
        app = MockApp(record=tiny_record)
        before = len(app._current_record.features)
        r = sc._h_transfer_annotations(
            app, {"source_id": "src", "apply": False})
        assert isinstance(r, dict) and r["applied"] is False, r
        assert len(app._current_record.features) == before

    def test_apply_overrides_dry_run(self, tiny_record, isolated_library):
        # An explicit apply wins over a (conflicting) dry_run.
        self._seed_self_source(tiny_record)
        app = MockApp(record=tiny_record)
        r = sc._h_transfer_annotations(
            app, {"source_id": "src", "apply": True, "dry_run": True})
        assert isinstance(r, dict) and r["applied"] is True, r

    def test_bare_apply_synonym_rejected(self, tiny_record, isolated_library):
        # A `commit`/`write`/`persist` with no apply/dry_run is a loud 400,
        # not a silent dry run (SC-D outcome-changing-ignored-key footgun).
        self._seed_self_source(tiny_record)
        app = MockApp(record=tiny_record)
        for syn in ("commit", "write", "persist"):
            r = sc._h_transfer_annotations(
                app, {"source_id": "src", syn: True})
            assert isinstance(r, tuple) and r[1] == 400, (syn, r)
            assert "apply" in r[0]["error"].lower()
        # A FALSY synonym means "don't apply" — the dry-run default, not a 400.
        ok = sc._h_transfer_annotations(
            app, {"source_id": "src", "commit": False})
        assert isinstance(ok, dict) and ok["applied"] is False, ok

    # ── Cross-collection source resolution (agent-API feedback) ─────────
    def _seed_source_in_other(self, tiny_record):
        """Active collection empty; the source lives in a NON-active one."""
        gb = sc._record_to_gb_text(tiny_record)
        sc._save_collections([
            {"name": "Active", "plasmids": [], "saved": "2026-01-01"},
            {"name": "Other", "plasmids": [{
                "id": "BACKBONE", "name": "Backbone Src",
                "size": len(tiny_record.seq),
                "n_feats": len(tiny_record.features), "gb_text": gb,
            }], "saved": "2026-01-01"},
        ])
        sc._set_active_collection_name("Active")
        sc._restore_library_from_active_collection()

    def test_resolves_source_in_non_active_collection(
            self, tiny_record, isolated_library):
        # The win: a source in a NON-active collection resolves (pre-fix it
        # 404'd, forcing a temp-copy dance).
        self._seed_source_in_other(tiny_record)
        app = MockApp(record=tiny_record)
        r = sc._h_transfer_annotations(
            app, {"source_id": "Backbone Src", "dry_run": True})
        assert isinstance(r, dict), r        # resolved, not 404
        assert r["applied"] is False
        # The canonical id resolves cross-collection too.
        assert isinstance(sc._h_transfer_annotations(
            app, {"source_id": "BACKBONE", "dry_run": True}), dict)

    def test_source_collection_pins(self, tiny_record, isolated_library):
        self._seed_source_in_other(tiny_record)
        app = MockApp(record=tiny_record)
        # Pin to the right collection → resolves.
        assert isinstance(sc._h_transfer_annotations(app, {
            "source_id": "Backbone Src", "source_collection": "Other",
            "dry_run": True}), dict)
        # Pin to a collection that lacks it → 404 (no mirror fallback when
        # the search is explicitly scoped).
        r = sc._h_transfer_annotations(app, {
            "source_id": "Backbone Src", "source_collection": "Active",
            "dry_run": True})
        assert isinstance(r, tuple) and r[1] == 404, r

    def test_cross_collection_ambiguity_409(
            self, tiny_record, isolated_library):
        gb = sc._record_to_gb_text(tiny_record)
        dup = lambda i: {"id": f"dup{i}", "name": "Dup Src",
                         "size": len(tiny_record.seq), "gb_text": gb}
        sc._save_collections([
            {"name": "Active", "plasmids": [], "saved": "2026-01-01"},
            {"name": "X", "plasmids": [dup(1)], "saved": "2026-01-01"},
            {"name": "Y", "plasmids": [dup(2)], "saved": "2026-01-01"},
        ])
        sc._set_active_collection_name("Active")
        sc._restore_library_from_active_collection()
        app = MockApp(record=tiny_record)
        r = sc._h_transfer_annotations(
            app, {"source_id": "Dup Src", "dry_run": True})
        assert isinstance(r, tuple) and r[1] == 409, r
        # source_collection disambiguates.
        assert isinstance(sc._h_transfer_annotations(app, {
            "source_id": "Dup Src", "source_collection": "X",
            "dry_run": True}), dict)


class TestDiffPlasmidHandler:
    """`_h_diff_plasmid` runs `_pairwise_align` between the loaded record
    and a target library entry. Mirrors the GUI diff flow for agent
    consumption — the result dict is the same shape `AlignmentScreen`
    consumes."""

    def test_no_record_loaded_returns_422(self):
        app = MockApp(record=None)
        result = sc._h_diff_plasmid(app, {"target_id": "x"})
        assert result == ({"error": "no plasmid loaded"}, 422)

    def test_missing_target_id_returns_400(self, tiny_record):
        app = MockApp(record=tiny_record)
        result = sc._h_diff_plasmid(app, {})
        assert isinstance(result, tuple) and result[1] == 400

    def test_invalid_mode_rejected(self, tiny_record):
        app = MockApp(record=tiny_record)
        result = sc._h_diff_plasmid(app, {"target_id": "x", "mode": "wat"})
        assert isinstance(result, tuple) and result[1] == 400

    def test_target_not_in_library_returns_404(self, tiny_record,
                                                  isolated_library):
        sc._save_library([])
        app = MockApp(record=tiny_record)
        result = sc._h_diff_plasmid(app, {"target_id": "ghost"})
        assert isinstance(result, tuple) and result[1] == 404

    def test_successful_diff_returns_alignment(self, tiny_record,
                                                  isolated_library):
        # Load a target into the library; diff against current record.
        sc._save_library([{
            "id":      "tgt",
            "name":    "tgt",
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-06",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = MockApp(record=tiny_record)
        result = sc._h_diff_plasmid(app, {"target_id": "tgt"})
        assert isinstance(result, dict)
        assert result["ok"] is True
        assert result["target_id"] == "tgt"
        # Self-vs-self: 100% identity.
        r = result["result"]
        assert r["identity_pct"] == 100.0
        assert r["n_mismatches"] == 0

    def test_circular_rotation_auto_detected_for_circular_target(
            self, tiny_record, isolated_library,
    ):
        """When the target's topology annotation is `circular`, the
        endpoint runs the seed-kmer rotation probe and reports the
        offset alongside the alignment result. Regression for
        2026-05-14 audit finding."""
        # `tiny_record` is annotated `topology=circular` per conftest.
        sc._save_library([{
            "id":      "tgt",
            "name":    "tgt",
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-06",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = MockApp(record=tiny_record)
        result = sc._h_diff_plasmid(app, {"target_id": "tgt"})
        assert result["ok"] is True
        assert result["circular"] is True
        # Self-vs-self at the same origin: no rotation needed.
        assert result["rotation_offset"] == 0

    def test_circular_rotation_can_be_forced(self, tiny_record,
                                                isolated_library):
        """A linear target with `circular: true` in the payload runs
        the rotation probe regardless of annotation."""
        from Bio.SeqRecord import SeqRecord
        # Re-stamp as linear so auto-detect would skip rotation.
        linear_rec = SeqRecord(
            tiny_record.seq, id=tiny_record.id, name=tiny_record.name,
            features=list(tiny_record.features),
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        sc._save_library([{
            "id":      "tgt",
            "name":    "tgt",
            "size":    len(linear_rec.seq),
            "n_feats": 0,
            "added":   "2026-05-06",
            "gb_text": sc._record_to_gb_text(linear_rec),
        }])
        app = MockApp(record=tiny_record)
        result = sc._h_diff_plasmid(
            app, {"target_id": "tgt", "circular": True},
        )
        assert result["ok"] is True
        assert result["circular"] is True

    def test_circular_rotation_can_be_disabled(self, tiny_record,
                                                isolated_library):
        """`circular: false` skips the rotation even for circular
        targets — preserves the pre-0.8.4 behaviour when callers want
        it."""
        sc._save_library([{
            "id":      "tgt",
            "name":    "tgt",
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-06",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = MockApp(record=tiny_record)
        result = sc._h_diff_plasmid(
            app, {"target_id": "tgt", "circular": False},
        )
        assert result["ok"] is True
        assert result["circular"] is False
        assert result["rotation_offset"] == 0


class TestPlasmidsaurusEndpoints:
    """Plasmidsaurus zip alignment endpoints — list-plasmidsaurus-members
    + align-plasmidsaurus-zip.

    Both endpoints take a real path on disk; the tests synthesize a
    minimal zip with one `.gbk` member so the parse + alignment
    pipeline can exercise them without a network round-trip.
    """

    def _make_zip(self, tmp_path, record, member_name: str = "run.gbk"):
        """Build a single-member `.zip` containing the given record as
        GenBank text. Returns the path."""
        import zipfile
        zip_path = tmp_path / "plasmidsaurus.zip"
        gb_text = sc._record_to_gb_text(record)
        with zipfile.ZipFile(str(zip_path), "w") as zf:
            zf.writestr(member_name, gb_text)
        return zip_path

    def test_list_members_returns_gbk_files(self, tiny_record, tmp_path):
        zip_path = self._make_zip(tmp_path, tiny_record)
        result = sc._h_list_plasmidsaurus_members(
            MockApp(), {"path": str(zip_path)},
        )
        assert isinstance(result, dict)
        assert result["ok"] is True
        assert result["count"] == 1
        assert result["members"][0]["name"] == "run.gbk"
        assert result["members"][0]["size"] > 0

    def test_list_members_missing_path_returns_400(self):
        result = sc._h_list_plasmidsaurus_members(MockApp(), {})
        assert isinstance(result, tuple)
        assert result[1] == 400

    def test_list_members_nonexistent_path_returns_400(self, tmp_path):
        # Sweep #25 (2026-05-23): collapsed 422 → uniform 400 (FS
        # oracle reduction — see `_h_list_plasmidsaurus_members`).
        result = sc._h_list_plasmidsaurus_members(
            MockApp(), {"path": str(tmp_path / "does-not-exist.zip")},
        )
        assert isinstance(result, tuple)
        assert result[1] == 400

    def test_list_members_non_zip_rejected(self, tmp_path):
        # Sweep #25: collapsed 422 → 400.
        bogus = tmp_path / "not-a-zip.zip"
        bogus.write_text("hello world")
        result = sc._h_list_plasmidsaurus_members(
            MockApp(), {"path": str(bogus)},
        )
        assert isinstance(result, tuple)
        assert result[1] == 400

    def test_list_members_filters_non_gbk(self, tiny_record, tmp_path):
        """Members with non-`.gbk`/`.gb`/`.genbank` extensions are
        skipped so the agent gets the same picker view the UI uses."""
        import zipfile
        zip_path = tmp_path / "mixed.zip"
        gb_text = sc._record_to_gb_text(tiny_record)
        with zipfile.ZipFile(str(zip_path), "w") as zf:
            zf.writestr("run.gbk",     gb_text)
            zf.writestr("readme.txt",  "ignore me")
            zf.writestr("data.csv",    "a,b,c")
        result = sc._h_list_plasmidsaurus_members(
            MockApp(), {"path": str(zip_path)},
        )
        assert result["ok"] is True
        assert {m["name"] for m in result["members"]} == {"run.gbk"}

    def test_align_self_vs_self_100pct(self, tiny_record, tmp_path,
                                          isolated_library):
        zip_path = self._make_zip(tmp_path, tiny_record)
        sc._save_library([{
            "id":      "tgt",
            "name":    "tgt",
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-06",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        result = sc._h_align_plasmidsaurus_zip(
            MockApp(),
            {
                "path":      str(zip_path),
                "member":    "run.gbk",
                "target_id": "tgt",
            },
        )
        assert isinstance(result, dict)
        assert result["ok"] is True
        assert result["target_id"] == "tgt"
        # Self-vs-self: 100% identity, no rotation needed.
        assert result["result"]["identity_pct"] == 100.0
        assert result["rotation_offset"] == 0
        # `tiny_record` is circular so the endpoint auto-detected it.
        assert result["circular"] is True

    def test_align_resolves_target_by_name(self, tiny_record, tmp_path,
                                             isolated_library):
        """`target_name` is a fallback when the agent doesn't know
        the id. Mirrors `_h_delete_from_library`'s name-based lookup
        contract — the library-entry's display name is the lookup
        key, while the returned `target_name` is the parsed LOCUS
        name from the gb_text (matches `_h_diff_plasmid`)."""
        zip_path = self._make_zip(tmp_path, tiny_record)
        sc._save_library([{
            "id":      "tgt",
            "name":    "Looked Up By Name",
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-06",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        result = sc._h_align_plasmidsaurus_zip(
            MockApp(),
            {
                "path":        str(zip_path),
                "member":      "run.gbk",
                "target_name": "Looked Up By Name",
            },
        )
        assert result["ok"] is True
        # The lookup matched by display name; returned `target_id` is
        # the library entry's id, `target_name` is the parsed LOCUS
        # name from the gb_text (TEST001 here per `tiny_record`).
        assert result["target_id"] == "tgt"
        assert result["target_name"] == "TEST001"

    def test_align_missing_target_returns_404(self, tiny_record,
                                                 tmp_path,
                                                 isolated_library):
        zip_path = self._make_zip(tmp_path, tiny_record)
        sc._save_library([])
        result = sc._h_align_plasmidsaurus_zip(
            MockApp(),
            {
                "path":      str(zip_path),
                "member":    "run.gbk",
                "target_id": "ghost",
            },
        )
        assert isinstance(result, tuple)
        assert result[1] == 404

    def test_align_missing_member_returns_400(self, tiny_record,
                                                 tmp_path):
        zip_path = self._make_zip(tmp_path, tiny_record)
        result = sc._h_align_plasmidsaurus_zip(
            MockApp(),
            {"path": str(zip_path), "target_id": "x"},
        )
        assert isinstance(result, tuple)
        assert result[1] == 400

    def test_align_unknown_zip_member_returns_422(self, tiny_record,
                                                     tmp_path,
                                                     isolated_library):
        zip_path = self._make_zip(tmp_path, tiny_record)
        sc._save_library([{
            "id":      "tgt",
            "name":    "tgt",
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-06",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        result = sc._h_align_plasmidsaurus_zip(
            MockApp(),
            {
                "path":      str(zip_path),
                "member":    "not-in-zip.gbk",
                "target_id": "tgt",
            },
        )
        assert isinstance(result, tuple)
        assert result[1] == 422

    def test_align_invalid_mode_rejected(self, tiny_record, tmp_path):
        zip_path = self._make_zip(tmp_path, tiny_record)
        result = sc._h_align_plasmidsaurus_zip(
            MockApp(),
            {
                "path":      str(zip_path),
                "member":    "run.gbk",
                "target_id": "tgt",
                "mode":      "fast",
            },
        )
        assert isinstance(result, tuple)
        assert result[1] == 400


class TestFindOrfsHandler:
    """`_h_find_orfs` exposes the six-frame ORF scan (added 0.6.0.0).
    Wraps `_find_orfs` — the algorithm itself is covered by
    test_dna_sanity.py::TestFindOrfs; here we just verify the agent
    path returns/normalises shape correctly."""

    def test_empty_when_no_record(self):
        app = MockApp(record=None)
        result = sc._h_find_orfs(app, {})
        assert result == ({"error": "no plasmid loaded"}, 422)

    def test_default_min_aa(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        body = "ATG" + "GCC" * 30 + "TAA"   # 99 bp, 31 aa coding
        rec = SeqRecord(Seq(body + "G" * 21), id="t", name="t")
        rec.annotations["topology"] = "circular"
        app = MockApp(record=rec)
        result = sc._h_find_orfs(app, {})
        assert "orfs" in result and "count" in result
        assert result["count"] >= 1
        # The ATG-stop ORF we built must be present.
        starts = {(o["start"], o["strand"]) for o in result["orfs"]}
        assert (0, 1) in starts

    def test_min_aa_filter(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        body = "ATG" + "GCC" * 19 + "TAA"   # 20 aa coding
        rec = SeqRecord(Seq(body), id="t", name="t")
        app = MockApp(record=rec)
        # 30 aa rejects, 20 aa keeps.
        assert sc._h_find_orfs(app, {"min_aa": 30})["count"] == 0
        assert sc._h_find_orfs(app, {"min_aa": 20})["count"] >= 1

    def test_min_aa_invalid(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ATGAAA"), id="t", name="t")
        app = MockApp(record=rec)
        result = sc._h_find_orfs(app, {"min_aa": "notanint"})
        assert isinstance(result, tuple) and result[1] == 400

    def test_min_aa_below_one_rejected(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ATGAAA"), id="t", name="t")
        app = MockApp(record=rec)
        result = sc._h_find_orfs(app, {"min_aa": 0})
        assert isinstance(result, tuple) and result[1] == 400

    def test_min_length_rejected_not_silently_ignored(self):
        # Agent-API audit (2026-06): `min_length` / `min_bp` are the
        # natural-but-wrong keys — `list-restriction-sites` uses `min_length`
        # in BASE PAIRS, so reusing it here used to be silently dropped (read
        # endpoints don't echo `ignored`) and return a default-30 ORF set that
        # looked authoritative. Now it 400s so the bp-vs-aa mix-up is loud.
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ATG" + "GCC" * 30 + "TAA"), id="t", name="t")
        app = MockApp(record=rec)
        for wrong in ("min_length", "min_bp", "min_len"):
            result = sc._h_find_orfs(app, {wrong: 10})
            assert isinstance(result, tuple) and result[1] == 400, (wrong, result)
            assert "min_aa" in result[0]["error"]
        # The correct key still works (no false trip).
        assert sc._h_find_orfs(app, {"min_aa": 20})["count"] >= 1

    def test_empty_seq_returns_empty_orf_list(self):
        """Regression guard for 2026-05-06: an empty `rec.seq` used to
        traverse `(rec.annotations or {})` — fine — but
        `_find_orfs(seq="")` itself returned `[]`; the agent path
        wrapper is now explicit about the shortcut so a missing /
        empty annotations dict can't surprise."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq(""), id="empty", name="empty")
        # Force a missing annotations attr to mimic a freshly-built
        # record from a partial Biopython parse.
        rec.annotations = None
        app = MockApp(record=rec)
        result = sc._h_find_orfs(app, {})
        assert result == {"orfs": [], "count": 0}


class TestLoadFileSizeCap:
    """Regression guard for 2026-05-06 fix: `_h_load_file` previously
    had NO size cap on disk reads — a malicious or buggy agent script
    could load a 10 GB GenBank file and OOM the worker. Cap is now
    `_BULK_IMPORT_MAX_BYTES` (50 MB) with `force=true` override."""

    def test_oversized_file_rejected_with_400(self, tmp_path, monkeypatch):
        # Sweep #25 (2026-05-23): size-cap response collapsed
        # 413 → uniform 400 (FS oracle reduction — error body no
        # longer carries `size_bytes` / `cap_bytes`; details in logs).
        monkeypatch.setattr(sc, "_BULK_IMPORT_MAX_BYTES", 10)
        big = tmp_path / "huge.gb"
        big.write_bytes(b"X" * 100)
        app = MockApp()
        result = sc._h_load_file(app, {"path": str(big)})
        payload, status = result
        assert status == 400
        assert "log" in payload["error"].lower()

    def test_force_overrides_size_cap(self, tmp_path, monkeypatch, tiny_record):
        """Pass force=true and the cap is bypassed (matches GUI's
        "load anyway" confirmation)."""
        monkeypatch.setattr(sc, "_BULK_IMPORT_MAX_BYTES", 10)
        # Use a real GenBank file so load_genbank succeeds.
        gb = tmp_path / "ok.gb"
        from io import StringIO
        from Bio import SeqIO as _SeqIO
        sio = StringIO()
        _SeqIO.write([tiny_record], sio, "genbank")
        gb.write_text(sio.getvalue())
        app = MockApp()
        # Sweep #25 (2026-05-23): size-cap response collapsed 413 → 400
        # (FS oracle reduction). Without force still rejected; with
        # force still parses.
        result = sc._h_load_file(app, {"path": str(gb)})
        assert isinstance(result, tuple) and result[1] == 400
        # With force: parsed.
        result = sc._h_load_file(app, {"path": str(gb), "force": True})
        assert isinstance(result, dict) and result["ok"] is True

    def test_missing_path_returns_400(self):
        result = sc._h_load_file(MockApp(), {})
        assert result[1] == 400 and "missing" in result[0]["error"]

    def test_nonexistent_path_returns_400(self, tmp_path):
        # Sweep #25: collapsed 404 → 400 (FS oracle reduction).
        result = sc._h_load_file(MockApp(),
                                  {"path": str(tmp_path / "nope.gb")})
        assert result[1] == 400


# ── End-to-end HTTP tests (real socket + JSON wire format) ─────────────────────


@pytest.fixture
def http_server(tiny_app):
    """Bind a real `_AgentAPIServer` on a free port for the test
    duration. Yields `(base_url, token)`."""
    port = _free_port()
    token = "test-token-" + str(port)
    srv = sc._AgentAPIServer(("127.0.0.1", port), tiny_app, token)
    thread = threading.Thread(target=srv.serve_forever, daemon=True)
    thread.start()
    # Tiny settle so the listening socket is ready before the first
    # request — otherwise the very first urlopen() can race the bind.
    time.sleep(0.05)
    try:
        yield f"http://127.0.0.1:{port}", token, tiny_app
    finally:
        srv.shutdown()
        srv.server_close()


def _http(url: str, *, method: str = "GET", body: dict | None = None,
          token: str | None = None,
          timeout: float = 5.0) -> tuple[int, dict]:
    """Tiny urllib helper that returns `(status, json_payload)`."""
    data = json.dumps(body or {}).encode("utf-8") if method == "POST" else None
    req = urllib.request.Request(url, data=data, method=method)
    if token is not None:
        req.add_header("Authorization", f"Bearer {token}")
    if data is not None:
        req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        body_bytes = e.read() if e.fp else b""
        try:
            payload = json.loads(body_bytes.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            payload = {"error": body_bytes.decode("utf-8", errors="replace")}
        return e.code, payload


class TestHTTPRouting:
    def test_status_endpoint(self, http_server):
        base, token, app = http_server
        status, payload = _http(f"{base}/status", token=token)
        assert status == 200
        assert payload["loaded"] is True
        assert payload["version"] == sc.__version__

    def test_tools_endpoint_lists_routes(self, http_server):
        base, token, app = http_server
        status, payload = _http(f"{base}/tools", token=token)
        assert status == 200
        names = [ep["name"] for ep in payload["endpoints"]]
        assert "status" in names

    def test_unknown_endpoint_returns_404(self, http_server):
        base, token, app = http_server
        status, payload = _http(f"{base}/no-such-thing", token=token)
        assert status == 404
        assert "endpoints" in payload   # helpful self-discovery

    def test_root_path_returns_tools(self, http_server):
        base, token, app = http_server
        status, payload = _http(f"{base}/", token=token)
        assert status == 200
        assert "endpoints" in payload


class TestEnvelopeContract:
    """The `data` / `_stale` response envelope is added by the DISPATCHER
    (`_agent_data_envelope` + the stale stamp in `_AgentRequestHandler.
    _handle`) to AUTHENTICATED endpoints only. The two unauthenticated
    endpoints served *ahead* of the wrapper — `/healthz` (readiness) and
    `/tools` (self-describe) — stay bare.

    These drive the REAL socket so they observe the SERVED bytes, not a
    handler return value: the handler runs BELOW the enveloping layer, so
    a handler-only assertion (see `test_tools_handler_returns_bare_shape`)
    cannot catch a dispatcher regression. Concretely, dropping the
    `path_part == "healthz"` or `path_part in ("", "tools")` early-return
    would push the ~150 KB `/tools` corpus (or the `/healthz` probe)
    through the envelope — and the negative tests here would fail, while
    the positive control proves the envelope is genuinely active (so the
    bare-ness is a real signal, not the envelope being globally off)."""

    def test_tools_served_bare_no_data_no_stale(self, http_server):
        base, _token, _app = http_server
        # Unauthenticated on purpose — served ahead of the token check.
        status, payload = _http(f"{base}/tools", token=None)
        assert status == 200
        assert "endpoints" in payload
        assert "data" not in payload, (
            "/tools must stay bare — the dispatcher must NOT envelope it "
            "(mirroring `doc_full` under `data` doubles the largest response)")
        assert "_stale" not in payload

    def test_healthz_served_bare_no_data_no_stale(self, http_server):
        base, _token, _app = http_server
        # `/healthz` is the SECOND envelope-free endpoint (not a lone
        # `/tools` exception): a readiness probe, also served pre-auth.
        status, payload = _http(f"{base}/healthz", token=None)
        assert status == 200
        assert payload.get("ok") is True
        assert {"ok", "status", "version", "headless"} <= set(payload)
        assert "data" not in payload, (
            "/healthz is envelope-free too — the dispatcher must NOT wrap it")
        assert "_stale" not in payload

    def test_authenticated_endpoint_carries_data(self, http_server):
        # POSITIVE control: an authenticated endpoint DOES carry `data`,
        # proving the envelope is active — so the bare-ness of /tools and
        # /healthz above is a meaningful signal, not the wrapper being off.
        base, token, _app = http_server
        status, payload = _http(f"{base}/status", token=token)
        assert status == 200
        assert "data" in payload, (
            "authenticated success responses must carry the `data` envelope")
        # Superset property: the original ad-hoc keys survive alongside `data`.
        assert payload.get("loaded") is True


class TestHTTPAuth:
    def test_write_endpoint_refuses_no_token(self, http_server):
        base, _token, _app = http_server
        status, payload = _http(
            f"{base}/save", method="POST", body={}, token=None,
        )
        assert status == 401
        assert "token" in payload["error"]

    def test_write_endpoint_refuses_wrong_token(self, http_server):
        base, _token, _app = http_server
        status, payload = _http(
            f"{base}/save", method="POST", body={}, token="wrong",
        )
        assert status == 401

    def test_read_endpoint_requires_token(self, http_server):
        """Sweep #25 (2026-05-23): bearer token required on ALL
        endpoints, not just writers. The earlier "reads can't damage
        state" assumption ignored that several read endpoints
        (hmmscan, blast, list-library) leak filesystem state or
        consume CPU/RAM — concrete attack surface for any co-
        resident local process. `tools` stays unauthenticated as the
        self-describe entry point."""
        base, _token, _app = http_server
        # Without token: 401.
        status, _payload = _http(f"{base}/status", token=None)
        assert status == 401
        # `tools` stays open so clients can self-describe.
        status, _payload = _http(f"{base}/tools", token=None)
        assert status == 200
        # With token: 200.
        status, _payload = _http(f"{base}/status", token=_token)
        assert status == 200


class TestHTTPHardening:
    def test_body_size_cap_constant_is_set(self):
        """The handler exposes a `_MAX_BODY_BYTES` cap so a bogus
        `Content-Length: 9999999999` header can't park the handler
        thread on `rfile.read`. We don't drive the cap end-to-end via
        a real TCP request (the localhost half-open race triggers a
        broken-pipe on the client before the server's rejection
        response lands, which is a transport issue not an app one).
        Instead, guard the constant + behavior in `_read_body` via
        a unit-level check below."""
        assert sc._AgentRequestHandler._MAX_BODY_BYTES <= 1 << 20
        assert sc._AgentRequestHandler._MAX_BODY_BYTES >= 1 << 12

    def test_malformed_json_does_not_crash(self, http_server):
        """A non-JSON POST body should be treated as an empty payload
        — never a 500 from a parsing exception leaking out."""
        base, token, _app = http_server
        req = urllib.request.Request(
            f"{base}/add-feature", method="POST",
            data=b"this is not json {{{",
        )
        req.add_header("Authorization", f"Bearer {token}")
        req.add_header("Content-Type", "application/json")
        try:
            with urllib.request.urlopen(req, timeout=5) as resp:
                status = resp.status
        except urllib.error.HTTPError as exc:
            status = exc.code
        # Empty body → handler sees no `start` / `end` → 400.
        assert status == 400

    def test_deeply_nested_json_body_400(self, http_server):
        """A deeply-nested JSON body makes `json.loads` raise RecursionError
        (a RuntimeError, NOT a JSONDecodeError). `_read_body` must catch it and
        the dispatcher must return a clean 400 — pre-fix it escaped the handler
        and dropped the connection with a worker-thread traceback. The body
        stays under `_MAX_BODY_BYTES` so the server reads it fully (no half-open
        race) and the failure is purely in the parse step."""
        base, token, _app = http_server
        depth = 50_000
        data = (b"[" * depth) + (b"]" * depth)  # ~100 KB, well under 1 MiB
        req = urllib.request.Request(
            f"{base}/add-feature", method="POST", data=data,
        )
        req.add_header("Authorization", f"Bearer {token}")
        req.add_header("Content-Type", "application/json")
        try:
            with urllib.request.urlopen(req, timeout=5) as resp:
                status = resp.status
        except urllib.error.HTTPError as exc:
            status = exc.code
        assert status == 400


class TestHTTPAddFeature:
    def test_add_feature_round_trip(self, http_server, tiny_record):
        base, token, app = http_server
        n_before = len(app._current_record.features)
        status, payload = _http(
            f"{base}/add-feature", method="POST",
            body={"start": 30, "end": 40, "label": "agentTest",
                  "type": "misc_feature"},
            token=token,
        )
        assert status == 200, payload
        assert payload["ok"] is True
        # Feature was actually appended to the underlying record.
        assert len(app._current_record.features) == n_before + 1
        new = app._current_record.features[-1]
        assert new.type == "misc_feature"
        assert new.qualifiers["label"] == ["agentTest"]
        assert int(new.location.start) == 30
        assert int(new.location.end)   == 40
        # And the record was marked dirty (so the user sees `*`).
        assert app._unsaved is True

    def test_add_feature_no_dirty_guard(self, http_server):
        """P1-10 (agent-API feedback): an additive feature edit does NOT trip
        the dirty guard. A `*`-dirty canvas — whether from a prior agent
        edit or a user edit in progress — still accepts the add WITHOUT
        `force`. Adding a feature can't discard unsaved work the way a
        load/replace can, so the old guard only forced `force:true` on
        every call after the first (each add re-dirties the canvas)."""
        base, token, app = http_server
        app._unsaved = True
        status, payload = _http(
            f"{base}/add-feature", method="POST",
            body={"start": 30, "end": 40, "label": "additive"},
            token=token,
        )
        assert status == 200, payload
        assert payload["ok"] is True

    def test_add_feature_loop_without_force(self, http_server):
        """The concrete P1-10 scenario: add several features in a row.
        The first dirties the canvas; pre-fix the second 409'd. Now every
        add lands without per-call `force`."""
        base, token, app = http_server
        app._unsaved = False
        for i in range(3):
            status, payload = _http(
                f"{base}/add-feature", method="POST",
                body={"start": 10 + i, "end": 20 + i, "label": f"f{i}"},
                token=token,
            )
            assert status == 200, (i, payload)

    def test_add_feature_force_still_accepted(self, http_server):
        """`force` is now a harmless no-op (the guard is gone) but must
        still be ACCEPTED — never 400'd or echoed under `ignored` — so
        existing agent scripts that pass it keep working."""
        base, token, app = http_server
        app._unsaved = True
        status, payload = _http(
            f"{base}/add-feature", method="POST",
            body={"start": 30, "end": 40, "label": "forced",
                  "force": True},
            token=token,
        )
        assert status == 200, payload
        assert "force" not in payload.get("ignored", [])

    def test_add_wrap_feature(self, http_server, tiny_record):
        """Wrap features (end < start) build a CompoundLocation."""
        base, token, app = http_server
        n = len(tiny_record.seq)
        status, payload = _http(
            f"{base}/add-feature", method="POST",
            body={"start": n - 5, "end": 5, "label": "wrap",
                  "type": "misc_feature"},
            token=token,
        )
        assert status == 200, payload
        from Bio.SeqFeature import CompoundLocation
        new = app._current_record.features[-1]
        assert isinstance(new.location, CompoundLocation)


class TestHTTPRegistration:
    def test_endpoint_decorator_registers(self):
        """Sanity: the decorator populates `_AGENT_HANDLERS` and tags
        write endpoints correctly. Catches a refactor that drops the
        registry by accident."""
        assert "status"      in sc._AGENT_HANDLERS
        assert "add-feature" in sc._AGENT_HANDLERS
        _fn, write = sc._AGENT_HANDLERS["status"]
        assert write is False
        _fn, write = sc._AGENT_HANDLERS["add-feature"]
        assert write is True

    def test_token_file_written_and_cleaned_up(self, tmp_path,
                                                  monkeypatch):
        """`_start_agent_api` writes (port, token) to the token file
        and `_stop_agent_api` removes it."""
        token_path = tmp_path / "agent_token"
        monkeypatch.setattr(sc._state, "_AGENT_TOKEN_FILE", token_path)
        port = _free_port()
        app = MockApp()
        srv = sc._start_agent_api(app, port=port)
        try:
            assert srv is not None
            assert token_path.exists()
            text = token_path.read_text(encoding="utf-8")
            stored_port, stored_token = text.strip().splitlines()
            assert int(stored_port) == port
            assert len(stored_token) >= 16
        finally:
            sc._stop_agent_api(srv)
        # Token file is removed on shutdown so a stale CLI invocation
        # can't accidentally hit a different process that bound the
        # same port later.
        assert not token_path.exists()


# ── Input sanitization (2026-05-01 hardening pass) ────────────────────────────


class TestSanitizeLabel:
    def test_strips_control_chars(self):
        assert sc._sanitize_label("hello\x00\x01world") == "helloworld"

    def test_collapses_newlines(self):
        # CR/LF would corrupt the sidebar's single-row label render.
        assert "\n" not in sc._sanitize_label("a\nb\rc")
        assert "\r" not in sc._sanitize_label("a\nb\rc")

    def test_caps_length(self):
        assert len(sc._sanitize_label("a" * 1000)) == 200
        assert len(sc._sanitize_label("a" * 1000, max_len=10)) == 10

    def test_unicode_survives(self):
        # Emoji + IUPAC-style ASCII labels both legitimate.
        assert sc._sanitize_label("test 🧬 lacZ") == "test 🧬 lacZ"

    def test_empty_returns_empty(self):
        assert sc._sanitize_label(None) == ""
        assert sc._sanitize_label("") == ""
        assert sc._sanitize_label("   ") == ""


class TestSanitizeFeatType:
    def test_default_for_empty(self):
        assert sc._sanitize_feat_type(None) == "misc_feature"
        assert sc._sanitize_feat_type("") == "misc_feature"
        assert sc._sanitize_feat_type("  ") == "misc_feature"

    def test_strips_control_chars(self):
        assert sc._sanitize_feat_type("CDS\x00") == "CDS"

    def test_caps_length(self):
        assert len(sc._sanitize_feat_type("a" * 100)) == 50


class TestSanitizeAccession:
    def test_valid(self):
        assert sc._sanitize_accession("L09137") == "L09137"
        assert sc._sanitize_accession("MW463917.1") == "MW463917.1"
        assert sc._sanitize_accession("NC_001140") == "NC_001140"

    def test_rejects_shell_metacharacters(self):
        # Defends against `accession=L09137; rm -rf /` smuggling.
        assert sc._sanitize_accession("L09137; rm -rf /") is None
        assert sc._sanitize_accession("L09137|cat /etc/passwd") is None
        assert sc._sanitize_accession("../../etc/hosts") is None

    def test_rejects_overlong(self):
        assert sc._sanitize_accession("A" * 33) is None

    def test_empty_returns_none(self):
        assert sc._sanitize_accession(None) is None
        assert sc._sanitize_accession("") is None


class TestSanitizeBases:
    def test_valid_iupac(self):
        s, err = sc._sanitize_bases("acgtnRYWSMKBDHV")
        assert err is None
        assert s == "ACGTNRYWSMKBDHV"

    def test_invalid_char(self):
        s, err = sc._sanitize_bases("ACGZ")
        assert err is not None
        assert "Z" in err

    def test_overlong(self):
        s, err = sc._sanitize_bases("A" * 100, max_len=50)
        assert err is not None
        assert "too long" in err

    def test_missing(self):
        s, err = sc._sanitize_bases(None)
        assert err is not None and "missing" in err


class TestEndpointHardening:
    """Adversarial-input tests: each endpoint must reject malformed
    payloads with a clear 400 error rather than crash or silently
    accept dangerous input."""

    def test_fetch_rejects_shell_meta(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/fetch", method="POST",
            body={"accession": "L09137; rm -rf /"},
            token=token,
        )
        assert status == 400
        assert "accession" in payload.get("error", "")

    def test_add_feature_strips_control_chars_in_label(self, http_server,
                                                        tiny_record):
        base, token, app = http_server
        status, payload = _http(
            f"{base}/add-feature", method="POST",
            body={"start": 30, "end": 40,
                  "label": "evil\x00\nlabel", "type": "misc_feature"},
            token=token,
        )
        assert status == 200, payload
        new = app._current_record.features[-1]
        assert "\x00" not in new.qualifiers["label"][0]
        assert "\n" not in new.qualifiers["label"][0]

    def test_add_feature_invalid_strand(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/add-feature", method="POST",
            body={"start": 30, "end": 40, "strand": 99},
            token=token,
        )
        assert status == 400


class TestTokenHardening:
    """Bearer-token comparison must be timing-safe and the token file
    written atomically with mode 0600 — a local-process attacker
    shouldn't be able to either time-leak the token byte-by-byte or
    race the chmod() to read the token in plaintext."""

    def test_token_compare_is_constant_time(self, http_server):
        # We can't directly time the comparison reliably enough to
        # detect a non-constant-time bug from a unit test, but we can
        # at least confirm `secrets.compare_digest` is in the call
        # path by verifying that two equal-length wrong tokens both
        # 401 (rather than 401-on-prefix-mismatch / 200-on-match).
        base, _token, _ = http_server
        wrong_a = "0" * 32
        wrong_b = "f" * 32
        s1, _ = _http(f"{base}/save", method="POST", body={},
                       token=wrong_a)
        s2, _ = _http(f"{base}/save", method="POST", body={},
                       token=wrong_b)
        assert s1 == s2 == 401

    def test_short_token_rejected_without_crash(self, http_server):
        base, _token, _ = http_server
        # Different length than the real token. compare_digest only
        # returns False here (doesn't raise). Pre-fix, this would
        # have hit a timing oracle; either way it must 401, not 500.
        status, _ = _http(f"{base}/save", method="POST", body={},
                           token="x")
        assert status == 401


class TestNewLibraryEndpoints:
    """Coverage for the parity endpoints added in the hardening pass:
    add-current-to-library, create-collection, delete-collection,
    rename-collection, set-active-collection, bulk-import-folder."""

    def test_create_collection_empty(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/create-collection", method="POST",
            body={"name": "agent-empty"}, token=token,
        )
        assert status == 200, payload
        assert payload["ok"] is True
        assert payload["n_plasmids"] == 0
        names = [c["name"] for c in sc._load_collections()]
        assert "agent-empty" in names

    def test_create_collection_rejects_blank(self, http_server):
        base, token, _ = http_server
        for bad in ("", "   ", "\x00\x00\x00", None):
            status, payload = _http(
                f"{base}/create-collection", method="POST",
                body={"name": bad}, token=token,
            )
            assert status == 400, (bad, payload)

    def test_create_collection_rejects_duplicate(self, http_server):
        sc._save_collections([{"name": "Existing", "plasmids": []}])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/create-collection", method="POST",
            body={"name": "Existing"}, token=token,
        )
        assert status == 409
        assert "already exists" in payload["error"]

    def test_create_collection_with_invalid_folder(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/create-collection", method="POST",
            body={"name": "agent-folder", "folder": "/nope/none/nada"},
            token=token,
        )
        assert status == 400
        assert "not a directory" in payload["error"]

    def test_delete_collection_round_trip(self, http_server):
        sc._save_collections([{"name": "Doomed", "plasmids": []}])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/delete-collection", method="POST",
            body={"name": "Doomed"}, token=token,
        )
        assert status == 200
        assert payload["deleted"] == "Doomed"
        names = [c["name"] for c in sc._load_collections()]
        assert "Doomed" not in names

    def test_delete_collection_404_on_missing(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/delete-collection", method="POST",
            body={"name": "GhostCollection"}, token=token,
        )
        assert status == 404

    def test_rename_collection_updates_active_pointer(self, http_server):
        sc._save_collections([{"name": "Old", "plasmids": []}])
        sc._set_active_collection_name("Old")
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/rename-collection", method="POST",
            body={"old": "Old", "new": "New"}, token=token,
        )
        assert status == 200
        assert sc._get_active_collection_name() == "New"

    def test_rename_collection_rejects_collision(self, http_server):
        sc._save_collections([
            {"name": "A", "plasmids": []},
            {"name": "B", "plasmids": []},
        ])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/rename-collection", method="POST",
            body={"old": "A", "new": "B"}, token=token,
        )
        assert status == 409

    def test_set_active_collection(self, http_server):
        sc._save_collections([
            {"name": "ColA", "plasmids": [
                {"id": "p1", "name": "p1", "size": 10, "gb_text": "X"}
            ]},
            {"name": "ColB", "plasmids": []},
        ])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-active-collection", method="POST",
            body={"name": "ColA"}, token=token,
        )
        assert status == 200
        assert sc._get_active_collection_name() == "ColA"
        assert payload["n_plasmids"] == 1

    def test_set_active_collection_404(self, http_server):
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/set-active-collection", method="POST",
            body={"name": "NotThere"}, token=token,
        )
        assert status == 404

    def test_bulk_import_folder_with_fixtures(self, http_server,
                                                isolated_library):
        from pathlib import Path
        fixtures_dir = Path(__file__).parent
        if not list(fixtures_dir.glob("FFE*.dna")):
            pytest.skip("No FFE .dna fixtures present")
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/bulk-import-folder", method="POST",
            body={"folder": str(fixtures_dir),
                  "collection": "FFE Bulk"},
            token=token,
        )
        assert status == 200, payload
        assert payload["n_imported"] >= 5
        assert payload["n_failed"] == 0
        names = [c["name"] for c in sc._load_collections()]
        assert "FFE Bulk" in names

    def test_bulk_import_folder_refuses_collection_collision(
        self, http_server, isolated_library
    ):
        sc._save_collections([{"name": "Taken", "plasmids": []}])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/bulk-import-folder", method="POST",
            body={"folder": "/tmp", "collection": "Taken"},
            token=token,
        )
        assert status == 409

    def test_bulk_import_folder_validates_folder(self, http_server):
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/bulk-import-folder", method="POST",
            body={"folder": "/no/such/dir/anywhere",
                  "collection": "X"},
            token=token,
        )
        assert status == 400

    def test_search_library_across_collections(self, http_server,
                                                  isolated_library):
        """`search-library` walks every collection on disk and returns
        fuzzy-matching plasmids regardless of which one is active."""
        sc._save_collections([
            {"name": "ColA", "plasmids": [
                {"id": "p1", "name": "pUC19_alpha", "size": 100,
                 "gb_text": "X", "n_feats": 3},
                {"id": "p2", "name": "pET28b", "size": 200,
                 "gb_text": "X", "n_feats": 4},
            ]},
            {"name": "ColB", "plasmids": [
                {"id": "p3", "name": "pUC19_beta", "size": 150,
                 "gb_text": "X", "n_feats": 5},
            ]},
        ])
        sc._set_active_collection_name("ColA")
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/search-library", method="POST",
            body={"query": "puc19"}, token=token,
        )
        assert status == 200, payload
        names = {(m["collection"], m["name"]) for m in payload["matches"]}
        assert ("ColA", "pUC19_alpha") in names
        assert ("ColB", "pUC19_beta") in names
        # pET28b doesn't match `puc19`.
        assert ("ColA", "pET28b") not in names

    def test_search_library_empty_query_lists_everything(
        self, http_server, isolated_library
    ):
        sc._save_collections([
            {"name": "X", "plasmids": [
                {"id": "a", "name": "a", "size": 1, "gb_text": "x"},
                {"id": "b", "name": "b", "size": 1, "gb_text": "x"},
            ]},
        ])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/search-library", method="POST", body={}, token=token,
        )
        assert status == 200
        assert payload["count"] == 2

    def test_search_library_limit_clamped(self, http_server,
                                            isolated_library):
        sc._save_collections([
            {"name": "X", "plasmids": [
                {"id": str(i), "name": f"p{i}", "size": 1, "gb_text": "x"}
                for i in range(50)
            ]},
        ])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/search-library", method="POST",
            body={"limit": 5}, token=token,
        )
        assert status == 200
        assert payload["count"] == 5

    def test_search_library_rejects_non_string_query(self, http_server):
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/search-library", method="POST",
            body={"query": 42}, token=token,
        )
        assert status == 400


class TestNewSearchEndpoints:
    """BLAST + HMMscan parity for agents."""

    def test_blast_returns_hits(self, http_server, isolated_library):
        # Library has the conftest's seeded plasmid; BLASTN against
        # itself should return at least one self-hit. Build the GenBank
        # text via Biopython so LOCUS length and ORIGIN bases agree
        # exactly (no parser warning).
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        bases = "ATGAAATTCCGATTAACCGGTTAAGGGCCATTTGCAAGGACCGGTTTAAA"
        rec = SeqRecord(Seq(bases), id="rec1", name="rec1",
                        annotations={"molecule_type": "DNA",
                                       "topology":      "circular"})
        sc._save_collections([{
            "name": "TestColl",
            "plasmids": [{
                "id":      "rec1",
                "name":    "rec1",
                "size":    len(bases),
                "gb_text": sc._record_to_gb_text(rec),
            }],
        }])
        base, token, _ = http_server
        # Long-enough query for pyhmmer (≥ 20 bp)
        status, payload = _http(
            f"{base}/blast", method="POST",
            body={"query": "ATGAAATTCCGATTAACCGGTTAAGGGCCATTTGC",
                  "program": "blastn", "backend": "pure"},
            token=token,
        )
        assert status == 200, payload
        assert payload["program"] == "blastn"
        assert payload["n_hits"] >= 1

    def test_blast_rejects_empty_query(self, http_server):
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/blast", method="POST",
            body={"query": "", "program": "blastn"},
            token=token,
        )
        assert status == 400

    def test_blast_rejects_invalid_program(self, http_server):
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/blast", method="POST",
            body={"query": "ATGC", "program": "tblastx"},
            token=token,
        )
        assert status == 400

    def test_blast_rejects_oversized_collection_list(self, http_server):
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/blast", method="POST",
            body={"query": "ATGCATGCATGCATGCATGC",
                  "collections": ["x"] * 200},
            token=token,
        )
        assert status == 400

    def test_blast_rejects_invalid_collection_name(self, http_server):
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/blast", method="POST",
            body={"query": "ATGCATGCATGCATGCATGC",
                  "collections": ["valid", "\x00\x00\x00"]},
            token=token,
        )
        assert status == 400

    def test_blast_clamps_max_hits(self, http_server, isolated_library):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/blast", method="POST",
            body={"query": "ATGCATGCATGCATGCATGC",
                  "max_hits": 99999, "backend": "pure"},
            token=token,
        )
        # Clamped to 500 internally; the search itself succeeds.
        assert status == 200

    def test_blast_invalid_backend(self, http_server):
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/blast", method="POST",
            body={"query": "ATGC", "backend": "xyz"},
            token=token,
        )
        assert status == 400

    def test_hmmscan_400_on_missing_path(self, http_server):
        # Sweep #11 (2026-05-20): hmmscan no longer surfaces a 404
        # distinct from 400 — that error differential was a
        # filesystem-state oracle for unauthenticated local
        # processes. All file-not-acceptable responses (not found,
        # symlink, not a regular file, oversize) collapse to a
        # single generic 400 with detail logged for the user.
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/hmmscan", method="POST",
            body={"query": "MAEELFKWILR" * 5,
                  "hmm_path": "/no/such/file.hmm"},
            token=token,
        )
        assert status == 400
        assert "not acceptable" in (payload.get("error") or "").lower()

    def test_hmmscan_400_on_short_query(self, http_server, tmp_path):
        # Build a tiny .hmm so the path-exists check passes; the
        # query length check rejects before we hit pyhmmer.
        fake = tmp_path / "fake.hmm"
        fake.write_text("")  # empty file is enough for the existence check
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/hmmscan", method="POST",
            body={"query": "M",  # below _HMMSCAN_MIN_QUERY_LEN
                  "hmm_path": str(fake)},
            token=token,
        )
        assert status == 400


class TestAdditionalAgentHardening:
    """A grab-bag of attack inputs against the new endpoints — none
    must crash the server or accept dangerous payloads."""

    def test_create_collection_rejects_oversized_name(self, http_server):
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/create-collection", method="POST",
            body={"name": "x" * (sc._MAX_COLLECTION_NAME_LEN + 100)},
            token=token,
        )
        # Long name is *truncated* to the cap, not rejected, so the
        # collection is still created successfully — the cap protects
        # against megabyte-sized JSON, not from semantic validation.
        assert status == 200
        names = [c["name"] for c in sc._load_collections()]
        assert any(len(n) <= sc._MAX_COLLECTION_NAME_LEN for n in names)

    def test_rename_collection_old_equals_new(self, http_server):
        sc._save_collections([{"name": "Same", "plasmids": []}])
        base, token, _ = http_server
        status, _ = _http(
            f"{base}/rename-collection", method="POST",
            body={"old": "Same", "new": "Same"}, token=token,
        )
        assert status == 400

    def test_add_current_to_library_no_record(self, http_server):
        base, token, app = http_server
        # Nuke the current record on the mock app
        app._current_record = None
        status, _ = _http(
            f"{base}/add-current-to-library", method="POST",
            body={}, token=token,
        )
        assert status == 422


class TestTypeStrictSanitisation:
    """Sanitisers must reject non-string inputs (dicts, lists, ints,
    None) rather than silently coerce via ``str()``. A JSON payload
    of ``{"name": {"x": 1}}`` should NOT become a collection literally
    named ``"{'x': 1}"``."""

    def test_sanitize_label_rejects_non_string(self):
        # Each of these used to be accepted via str() coercion; now
        # they must come back as empty.
        assert sc._sanitize_label({"x": 1}) == ""
        assert sc._sanitize_label([1, 2, 3]) == ""
        assert sc._sanitize_label(42) == ""
        assert sc._sanitize_label(None) == ""
        assert sc._sanitize_label(True) == ""

    def test_sanitize_feat_type_rejects_non_string(self):
        assert sc._sanitize_feat_type({"x": 1}) == "misc_feature"
        assert sc._sanitize_feat_type(42)        == "misc_feature"
        assert sc._sanitize_feat_type(None)      == "misc_feature"

    def test_sanitize_accession_rejects_non_string(self):
        assert sc._sanitize_accession({"x": 1}) is None
        assert sc._sanitize_accession([1, 2])   is None
        assert sc._sanitize_accession(42)       is None

    def test_sanitize_path_rejects_non_string(self):
        assert sc._sanitize_path({"x": 1}) is None
        assert sc._sanitize_path([1, 2])   is None
        assert sc._sanitize_path(42)       is None

    def test_create_collection_rejects_non_string_name(self, http_server):
        base, token, _ = http_server
        for bad in ({"x": 1}, [1, 2, 3], 42, None):
            status, payload = _http(
                f"{base}/create-collection", method="POST",
                body={"name": bad}, token=token,
            )
            assert status == 400, (bad, payload)


class TestNumericCoercionHardening:
    """``int(x)`` blows up on ``+/- Infinity`` (OverflowError) and
    ``NaN`` returns silently as 0 in some paths. Both must be caught
    cleanly at every numeric input boundary so a hostile JSON payload
    can't crash the handler thread."""

    def test_blast_max_hits_infinity(self, http_server):
        base, token, _ = http_server
        body_json = '{"query": "ATGCATGCATGCATGCATGC", "max_hits": Infinity}'
        # Send raw JSON (urllib helper auto-encodes a dict, but we
        # want literal JSON Infinity which Python's json.loads accepts).
        import urllib.request, urllib.error
        req = urllib.request.Request(
            f"{base}/blast", data=body_json.encode(),
            headers={"Authorization": f"Bearer {token}",
                       "Content-Type": "application/json"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            status = 200
        except urllib.error.HTTPError as exc:
            status = exc.code
        # Pre-fix this was a 500 (OverflowError trace). Now must 400
        # with a clear "must be a finite number" message.
        assert status == 400

    def test_add_feature_start_infinity(self, http_server, tiny_record):
        base, token, _ = http_server
        body_json = ('{"start": Infinity, "end": 10, "label": "x"}')
        import urllib.request, urllib.error
        req = urllib.request.Request(
            f"{base}/add-feature", data=body_json.encode(),
            headers={"Authorization": f"Bearer {token}",
                       "Content-Type": "application/json"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            status = 200
        except urllib.error.HTTPError as exc:
            status = exc.code
        assert status == 400

    def test_add_feature_strand_infinity(self, http_server, tiny_record):
        """Regression guard for 2026-05-05 retrofit: `add-feature` used
        to call raw `int(payload.get("strand", 1))` which raises
        OverflowError on JSON `Infinity`. Now routes through
        `_coerce_int` and returns a clean 400."""
        base, token, _ = http_server
        body_json = ('{"start": 0, "end": 10, "label": "x", '
                       '"strand": Infinity}')
        import urllib.request, urllib.error
        req = urllib.request.Request(
            f"{base}/add-feature", data=body_json.encode(),
            headers={"Authorization": f"Bearer {token}",
                       "Content-Type": "application/json"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            status = 200
        except urllib.error.HTTPError as exc:
            status = exc.code
        assert status == 400

    def test_update_feature_strand_infinity(self, http_server, tiny_record):
        """Regression guard for 2026-05-05 retrofit: same fix on
        `update-feature`'s optional strand field."""
        base, token, _ = http_server
        body_json = '{"idx": 0, "strand": Infinity}'
        import urllib.request, urllib.error
        req = urllib.request.Request(
            f"{base}/update-feature", data=body_json.encode(),
            headers={"Authorization": f"Bearer {token}",
                       "Content-Type": "application/json"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            status = 200
        except urllib.error.HTTPError as exc:
            status = exc.code
        assert status == 400

    def test_list_restriction_sites_min_length_infinity(self, http_server,
                                                          tiny_record):
        """Regression guard for 2026-05-05 retrofit: `list-restriction-
        sites` now rejects Infinity in `min_length` instead of bubbling
        an OverflowError up to the 500 path."""
        base, token, _ = http_server
        body_json = '{"min_length": Infinity}'
        import urllib.request, urllib.error
        req = urllib.request.Request(
            f"{base}/list-restriction-sites", data=body_json.encode(),
            headers={"Authorization": f"Bearer {token}",
                       "Content-Type": "application/json"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            status = 200
        except urllib.error.HTTPError as exc:
            status = exc.code
        assert status == 400

    def test_list_restriction_sites_rejects_non_string_enzymes(
            self, http_server, tiny_record):
        """Regression guard for 2026-05-17 audit fix: every element of
        `enzymes` must be a string. Pre-fix a mixed-type list like
        ``[1, 2.5, null]`` built a set whose ``not in`` check silently
        filtered every hit to zero — agents got an empty result with
        no signal that their payload was malformed."""
        base, token, _ = http_server
        body_json = '{"enzymes": [1, 2.5, null]}'
        import urllib.request, urllib.error
        req = urllib.request.Request(
            f"{base}/list-restriction-sites", data=body_json.encode(),
            headers={"Authorization": f"Bearer {token}",
                       "Content-Type": "application/json"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            status = 200
        except urllib.error.HTTPError as exc:
            status = exc.code
        assert status == 400

    def test_list_restriction_sites_accepts_all_string_enzymes(
            self, http_server, tiny_record):
        """Positive case for the 2026-05-17 type check: a well-formed
        all-string enzymes list must NOT 400."""
        base, token, _ = http_server
        body_json = '{"enzymes": ["EcoRI", "BamHI"]}'
        import urllib.request, urllib.error
        req = urllib.request.Request(
            f"{base}/list-restriction-sites", data=body_json.encode(),
            headers={"Authorization": f"Bearer {token}",
                       "Content-Type": "application/json"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            status = 200
        except urllib.error.HTTPError as exc:
            status = exc.code
        assert status == 200

    def test_list_restriction_sites_enzyme_singular_alias(self):
        """P1-3 (agent-API feedback): the singular `enzyme` is accepted as an
        alias for `enzymes`, and a bare string works under either key — so
        the natural `{enzyme: "EcoRI"}` guess is HONORED, not silently
        dropped (which would scan the whole catalog and look authoritative
        while answering a different question)."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # One EcoRI (GAATTC) + one BamHI (GGATCC) site, well separated.
        seq = "AAAAAA" + "GAATTC" + "AAAAAA" + "GGATCC" + "AAAAAA"
        rec = SeqRecord(Seq(seq), id="t", name="t")
        rec.annotations["topology"] = "linear"
        app = MockApp(record=rec)
        singular   = sc._h_list_restriction_sites(app, {"enzyme": "EcoRI"})
        plural     = sc._h_list_restriction_sites(app, {"enzymes": ["EcoRI"]})
        str_plural = sc._h_list_restriction_sites(app, {"enzymes": "EcoRI"})
        # Singular alias and bare-string forms == the canonical list form.
        assert singular == plural == str_plural
        # And the scan is actually scoped to EcoRI (BamHI present, excluded).
        assert {s["enzyme"] for s in singular["sites"]} == {"EcoRI"}, singular
        # Routing the OTHER enzyme by its singular key returns a different,
        # correctly-scoped result — proving the param is honored, not ignored.
        bam = sc._h_list_restriction_sites(app, {"enzyme": "BamHI"})
        assert {s["enzyme"] for s in bam["sites"]} == {"BamHI"}, bam


class TestListLibraryCollectionScope:
    """P1-7 (agent-API feedback): `list-library {collection}` scopes the
    listing to ONE collection's plasmids instead of silently returning the
    whole active-collection library (the SC-D "ignored routing param"
    footgun where the result looks scoped but isn't)."""

    def _ent(self, name, size):
        return {"id": f"id-{name}", "name": name, "gb_text": "",
                "size": size, "n_feats": 0}

    def _seed(self):
        sc._save_collections([
            {"name": "Alpha",
             "plasmids": [self._ent("pA1", 10), self._ent("pA2", 20)],
             "saved": "2026-01-01"},
            {"name": "Beta", "plasmids": [self._ent("pB1", 30)],
             "saved": "2026-01-01"},
        ])
        sc._set_active_collection_name("Alpha")
        sc._restore_library_from_active_collection()

    def test_scopes_to_named_non_active_collection(self):
        # Active is Alpha; asking for Beta must return Beta's plasmid, NOT
        # Alpha's — the whole point is reading a non-active collection
        # without a set-active round-trip.
        self._seed()
        r = sc._h_list_library(None, {"collection": "Beta"})
        assert isinstance(r, dict), r
        assert r["collection"] == "Beta"
        assert {e["name"] for e in r["library"]} == {"pB1"}
        assert r["count"] == 1

    def test_scopes_to_active_collection(self):
        self._seed()
        r = sc._h_list_library(None, {"collection": "Alpha"})
        assert {e["name"] for e in r["library"]} == {"pA1", "pA2"}

    def test_unknown_collection_404(self):
        self._seed()
        assert sc._h_list_library(None, {"collection": "Ghost"})[1] == 404

    def test_omitted_collection_returns_full_library_no_scope_key(self):
        self._seed()
        r = sc._h_list_library(None, {})
        # Unscoped → the active-collection library mirror, with NO
        # `collection` key in the response (only present when scoped).
        assert "collection" not in r
        assert {e["name"] for e in r["library"]} == {"pA1", "pA2"}


class TestRateLimitRetryAfter:
    """SC-L (agent-API feedback): a 429 carries a coarse `retry_after`
    (whole seconds, >=1) so a batching agent can self-pace instead of
    busy-retrying."""

    def test_429_carries_retry_after(self, http_server, monkeypatch):
        import time
        base, token, _ = http_server
        # Slow the refill so the bucket stays exhausted across the request
        # round-trip — the default 30 tokens/s would refill a read's worth
        # within the round-trip and flake the assertion.
        monkeypatch.setattr(sc, "_AGENT_RATE_LIMIT_TOKENS_PER_SEC", 1.0)
        sc._agent_rate_limit_reset()
        # Drain this token's bucket (recent timestamp) so the next request
        # 429s rather than refilling from a stale/zero timestamp.
        sc._AGENT_RATE_LIMIT_BUCKETS[token] = (0.0, time.monotonic())
        try:
            status, payload = _http(f"{base}/status", token=token)
        finally:
            sc._agent_rate_limit_reset()
        assert status == 429, payload
        assert isinstance(payload.get("retry_after"), int)
        assert payload["retry_after"] >= 1


class TestDigestHandler:
    """P1-4 (agent-API feedback): `digest` exposes the overhang-aware
    restriction engine on a RAW sequence so an agent can verify junction /
    dropout overhangs without a sandboxed `import splicecraft`."""

    # Two EcoRI (GAATTC) sites in a small circular molecule (52 bp).
    SEQ = "GAATTC" + "A" * 20 + "GAATTC" + "T" * 20

    def test_circular_two_cuts_two_fragments(self):
        r = sc._h_digest(None, {"sequence": self.SEQ, "enzymes": ["EcoRI"]})
        assert isinstance(r, dict), r
        assert r["n_cuts"] == 2
        assert r["n_fragments"] == 2          # circular: N cuts → N fragments
        assert r["circular"] is True
        assert {c["enzyme"] for c in r["cuts"]} == {"EcoRI"}
        # EcoRI leaves a 4-base 5' overhang (AATT) — the QC payload.
        for c in r["cuts"]:
            assert c["overhang_seq"] == "AATT", c
            assert c["kind"] == "5'", c

    def test_matches_real_engine(self):
        # Feedback "done when": the count matches `_enzyme_cuts`.
        r = sc._h_digest(None, {"sequence": self.SEQ, "enzymes": ["EcoRI"],
                                "circular": True})
        assert r["n_cuts"] == len(
            sc._enzyme_cuts(self.SEQ, ["EcoRI"], circular=True))

    def test_linear_n_plus_one_fragments(self):
        r = sc._h_digest(None, {"sequence": self.SEQ, "enzymes": ["EcoRI"],
                                "circular": False})
        # Linear: N cuts → N+1 fragments, free ends marked `linear`.
        assert r["n_fragments"] == r["n_cuts"] + 1
        assert r["fragments"][0]["left"]["kind"] == "linear"
        assert r["fragments"][-1]["right"]["kind"] == "linear"

    def test_singular_enzyme_and_seq_aliases(self):
        a = sc._h_digest(None, {"seq": self.SEQ, "enzyme": "EcoRI"})
        b = sc._h_digest(None, {"sequence": self.SEQ, "enzymes": ["EcoRI"]})
        assert a == b

    def test_unknown_enzyme_reported_not_dropped(self):
        r = sc._h_digest(None, {"sequence": self.SEQ,
                                "enzymes": ["EcoRI", "NotAnEnzyme"]})
        assert r["unknown_enzymes"] == ["NotAnEnzyme"]
        assert r["n_cuts"] == 2                # only EcoRI cuts

    def test_include_fragment_seq_toggle(self):
        off = sc._h_digest(None, {"sequence": self.SEQ, "enzymes": ["EcoRI"]})
        assert "seq" not in off["fragments"][0]
        on = sc._h_digest(None, {"sequence": self.SEQ, "enzymes": ["EcoRI"],
                                 "include_fragment_seq": True})
        assert "seq" in on["fragments"][0]
        # Top strands partition the circle → lengths sum to the molecule.
        assert sum(f["length"] for f in on["fragments"]) == len(self.SEQ)

    def test_validation(self):
        assert sc._h_digest(None, {"enzymes": ["EcoRI"]})[1] == 400  # no seq
        assert sc._h_digest(None, {"sequence": self.SEQ})[1] == 400  # no enz
        assert sc._h_digest(None,
                            {"sequence": self.SEQ, "enzymes": []})[1] == 400
        assert sc._h_digest(None,
                            {"sequence": self.SEQ,
                             "enzymes": [1, 2]})[1] == 400
        # Too many enzymes — DoS guard (max 500).
        assert sc._h_digest(None,
                            {"sequence": self.SEQ,
                             "enzymes": ["EcoRI"] * 501})[1] == 400


class _FakeLibPanel:
    """Records the records handed to `add_entry` (the real LibraryPanel
    isn't mounted under a MockApp)."""
    def __init__(self):
        self.added: list = []

    def add_entry(self, record, **kw):
        self.added.append(record)
        return True


class _AppWithFakeLib(MockApp):
    def __init__(self, record):
        super().__init__(record=record)
        self.fake_lib = _FakeLibPanel()

    def query_one(self, selector, *a, **k):  # type: ignore[override]
        return self.fake_lib


class TestAddCurrentNameIdOverride:
    """P0-8 / P1-9 (agent-API feedback): add-current-to-library accepts
    {id?, name?} so a batch of from-scratch builds (all id 'Cloned') can
    each be filed under a unique LOCUS-safe id + real spaced display name
    instead of silently overwriting each other."""

    def test_override_stamps_copy_not_canvas(self, tiny_record,
                                             isolated_library):
        sc._save_library([])
        app = _AppWithFakeLib(tiny_record)
        original_id = tiny_record.id
        r = sc._h_add_current_to_library(
            app, {"id": "FIRE 10!", "name": "FIRE 10 Module1"})
        assert isinstance(r, dict) and r["ok"] is True, r
        assert r["id"] == "FIRE_10"                 # scrubbed LOCUS-safe
        assert r["name"] == "FIRE 10 Module1"       # spaced display name
        # The record handed to add_entry carries the override identity...
        saved = app.fake_lib.added[-1]
        assert saved.id == "FIRE_10"
        assert getattr(saved, "_tui_display_name") == "FIRE 10 Module1"
        # ...but the LIVE canvas record is untouched (override on a copy).
        assert tiny_record.id == original_id
        assert tiny_record is app._current_record

    def test_distinct_ids_reach_add_entry(self, tiny_record,
                                          isolated_library):
        sc._save_library([])
        app = _AppWithFakeLib(tiny_record)
        sc._h_add_current_to_library(app, {"id": "build-1"})
        sc._h_add_current_to_library(app, {"id": "build-2"})
        assert [x.id for x in app.fake_lib.added] == ["build-1", "build-2"]

    def test_invalid_id_returns_400(self, tiny_record, isolated_library):
        app = _AppWithFakeLib(tiny_record)
        r = sc._h_add_current_to_library(app, {"id": "!!!"})
        assert isinstance(r, tuple) and r[1] == 400

    def test_no_override_preserves_old_behavior(self, tiny_record,
                                                isolated_library):
        sc._save_library([])
        app = _AppWithFakeLib(tiny_record)
        r = sc._h_add_current_to_library(app, {})
        assert r["ok"] is True
        # No copy made — the live record itself is filed.
        assert app.fake_lib.added[-1] is tiny_record

    def test_save_clears_dirty_flag(self, tiny_record, isolated_library):
        # Agent-API audit (2026-06): add-current persists the canvas, so the
        # canvas is no longer dirty. WITHOUT clearing it, a batch loop's next
        # new-plasmid / load-entry / set-active-collection 409s on the
        # dirty-guard and silently strands the build on the previous plasmid.
        sc._save_library([])
        app = _AppWithFakeLib(tiny_record)
        app._unsaved = True
        r = sc._h_add_current_to_library(app, {})
        assert r["ok"] is True
        assert app._unsaved is False

    def test_id_override_copy_leaves_canvas_dirty(self, tiny_record,
                                                  isolated_library):
        # An id/name override files a COPY (the live canvas keeps its own
        # identity), so the canvas is genuinely still unsaved and must stay
        # dirty — only a save of the LIVE record itself clears the flag.
        sc._save_library([])
        app = _AppWithFakeLib(tiny_record)
        app._unsaved = True
        r = sc._h_add_current_to_library(app, {"id": "build-x"})
        assert r["ok"] is True
        assert app._unsaved is True

    def test_long_distinct_ids_dont_truncate_collide(self, tiny_record,
                                                     isolated_library):
        # Two distinct ids sharing a 32-char prefix must stay DISTINCT — a
        # fixed-width truncation let them collide and silently overwrite,
        # the exact batch-loss this override prevents (code-review CONFIRMED).
        sc._save_library([])
        app = _AppWithFakeLib(tiny_record)
        r1 = sc._h_add_current_to_library(
            app, {"id": "cloned_plasmid_batch_run_001_sampleA"})
        r2 = sc._h_add_current_to_library(
            app, {"id": "cloned_plasmid_batch_run_001_sampleB"})
        assert r1["id"] != r2["id"]
        assert {x.id for x in app.fake_lib.added} == {r1["id"], r2["id"]}


class TestLoadFileNameIdOverride:
    """P0-8 / P1-9: load-file accepts {name?, id?} to stamp the loaded
    record's identity, fixing the from-scratch `.dna` 'first-word name'
    truncation + the id-'Cloned' batch collision at the source."""

    def _write_gb(self, tmp_path, locus="Cloned"):
        from io import StringIO
        from Bio import SeqIO as _SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ACGT" * 20), id=locus, name=locus,
                        description="Cloned part FIRE 10")
        rec.annotations["molecule_type"] = "DNA"
        sio = StringIO()
        _SeqIO.write([rec], sio, "genbank")
        gb = tmp_path / "build.gb"
        gb.write_text(sio.getvalue())
        return gb

    def test_override_name_and_id(self, tmp_path):
        gb = self._write_gb(tmp_path)
        app = MockApp()
        r = sc._h_load_file(app, {"path": str(gb),
                                  "id": "FIRE 10!", "name": "FIRE 10 Module1"})
        assert isinstance(r, dict) and r["ok"] is True, r
        assert r["id"] == "FIRE_10"                  # scrubbed LOCUS-safe
        assert r["name"] == "FIRE 10 Module1"        # spaced display name
        # The loaded canvas record carries the stamped identity.
        assert app._current_record.id == "FIRE_10"
        assert sc.PlasmidApp._record_display_name(app._current_record) \
            == "FIRE 10 Module1"

    def test_blank_id_after_scrub_returns_400(self, tmp_path):
        gb = self._write_gb(tmp_path)
        app = MockApp()
        r = sc._h_load_file(app, {"path": str(gb), "id": "!!!"})
        assert isinstance(r, tuple) and r[1] == 400

    def test_no_override_loads_normally(self, tmp_path):
        gb = self._write_gb(tmp_path, locus="MyPlasmid")
        app = MockApp()
        r = sc._h_load_file(app, {"path": str(gb)})
        assert isinstance(r, dict) and r["ok"] is True


class TestAgentAssemblyLineage:
    """P0-2 (agent-API feedback): agent assemblies (gibson / golden-gate /
    traditional) attach REAL parent lineage — each input fragment / part /
    vector becomes a parent node — instead of the flat `createDocument`
    leaf an agent-built construct used to get in the History tab."""

    def test_builder_insert_fragment_with_parents(self):
        xml = sc._build_history_for_agent_assembly(
            name="prod", seq_len=2000, circular=True,
            manipulation="goldenGateAssembly",
            parents=[{"name": "partA", "size": 800},
                     {"name": "partB", "size": 700},
                     {"name": "vector", "size": 500}])
        assert xml
        root = sc._parse_commercialsaas_history(xml)
        d = sc._history_node_to_dict(root)
        assert d["operation"] == "insertFragment"
        assert len(d["parents"]) == 3
        names = " ".join(p["name"] for p in d["parents"])
        assert "partA" in names and "partB" in names and "vector" in names

    def test_builder_empty_parents_still_valid(self):
        xml = sc._build_history_for_agent_assembly(
            name="prod", seq_len=100, circular=True,
            manipulation="gibsonAssembly", parents=[])
        root = sc._parse_commercialsaas_history(xml)
        d = sc._history_node_to_dict(root)
        assert d["operation"] == "insertFragment"
        assert d["parents"] == []

    def test_parent_label_helper(self):
        assert sc._agent_assembly_parent("ACGTACGT", "part 1") == {
            "name": "part 1", "size": 8}
        assert sc._agent_assembly_parent(
            {"sequence": "ACGT", "name": "FFE 28"}, "part 2") == {
            "name": "FFE 28", "size": 4}
        assert sc._agent_assembly_parent({"sequence": "AC"}, "vector") == {
            "name": "vector", "size": 2}

    # Known-good real digest+ligation inputs (mirrors the traditional-clone
    # test): EcoRI · backbone · BamHI · stuffer, circular vector.
    VECTOR = "GAATTC" + "A" * 200 + "GGATCC" + "C" * 60
    INSERT = "GCGC" + "GAATTC" + "ATGAAACCCGGGTTTTAA" + "GGATCC" + "GCGC"

    def test_traditional_clone_attaches_lineage_end_to_end(self):
        r = sc._h_traditional_clone(None, {
            "vector_seq": self.VECTOR, "vector_enzymes": ["EcoRI", "BamHI"],
            "insert_seq": self.INSERT, "insert_enzymes": ["EcoRI", "BamHI"],
            "product_name": "lineageclone", "vector_frag_idx": 0,
            "vector_name": "MyVector", "insert_name": "MyInsert"})
        assert r["ok"], r
        # get-history reads the saved entry's lineage tree.
        h = sc._h_get_history(None, {"name": "lineageclone"})
        assert h["ok"] and h["history"] is not None, h
        tree = h["history"]
        assert tree["operation"] == "insertFragment"   # not createDocument
        assert len(tree["parents"]) == 2               # vector + insert
        pnames = " ".join(p["name"] for p in tree["parents"])
        assert "MyVector" in pnames and "MyInsert" in pnames


class TestDesignGbPartEntryVector:
    """`design-gb-part` ``check_entry_vector``: does the designed part actually
    drop into the CONFIGURED entry vector? Runs the REAL clone sim so the
    verdict matches the clone step (agent-API feedback P1-5 residual)."""

    INSERT = "ATG" + "GGCTCAGGCGGT" * 5 + "TAA"   # 66 bp, no Esp3I/BsaI sites

    def _cds_type(self):
        g = sc._BUILTIN_GRAMMARS["gb_l0"]
        return next((p["type"] for p in g["positions"]
                     if p["type"] in (g.get("coding_types") or [])),
                    g["positions"][0]["type"])

    def _design(self):
        g = sc._BUILTIN_GRAMMARS["gb_l0"]
        d = sc._design_gb_primers(self.INSERT, 0, len(self.INSERT),
                                  self._cds_type(), grammar=g)
        assert not d.get("error"), d.get("error")
        return d

    def _bind_vector(self, oh5, oh3, name="CompatUPD"):
        # An Esp3I acceptor whose dropout releases EXACTLY (oh5, oh3); backbone
        # + dropout are Esp3I/BsaI-free so there are precisely two cut sites.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        backbone = "AAAACCCC" * 100
        dropout = "GGGGTTTT" * 20
        cassette = "CGTCTCA" + oh5 + dropout + oh3 + "AGAGACG"
        seq = backbone[:400] + cassette + backbone[400:]
        rec = SeqRecord(Seq(seq), id="TESTUPD", name="TESTUPD",
                        annotations={"molecule_type": "DNA",
                                       "topology": "circular"})
        sc._set_entry_vector("gb_l0", {
            "name": name, "id": "TESTUPD", "size": len(seq),
            "source": "test", "gb_text": sc._record_to_gb_text(rec)})

    def _call(self, **extra):
        return sc._h_design_gb_part(None, dict(
            {"template": self.INSERT, "start": 0, "end": len(self.INSERT),
             "part_type": self._cds_type(), "grammar": "gb_l0"}, **extra))

    def test_compatible_part_drops_in(self):
        d = self._design()
        self._bind_vector(d["oh5"], d["oh3"])
        chk = self._call(check_entry_vector=True)["result"]["entry_vector_check"]
        assert chk["checkable"] is True
        assert chk["compatible"] is True, chk
        assert chk["vector_name"] == "CompatUPD"
        assert chk["acceptor_oh5"] == d["oh5"] and chk["acceptor_oh3"] == d["oh3"]

    def test_design_auto_conforms_to_configured_vector(self):
        # 2026-07-13 two-tier fix: `design-gb-part` NESTS the grammar's category
        # overhangs inside whatever EXTERNAL overhangs the configured L0 acceptor
        # exposes, so the fragment drops in BY CONSTRUCTION — even for an unusual
        # external pair the part's category overhangs don't themselves match.
        d = self._design()
        assert (d["oh5"], d["oh3"]) != ("AAAA", "TTTT")
        self._bind_vector("AAAA", "TTTT", name="UnusualUPD")
        res = self._call(check_entry_vector=True)["result"]
        # actual Type IIS-cut overhangs now equal the vector's external pair …
        assert res["entry_oh5"] == "AAAA" and res["entry_oh3"] == "TTTT"
        # … while the category identity (what BsaI releases at L1) is unchanged
        assert res["oh5"] == d["oh5"] and res["oh3"] == d["oh3"]
        chk = res["entry_vector_check"]
        assert chk["checkable"] is True and chk["compatible"] is True, chk

    def test_check_still_flags_genuine_mismatch(self):
        # The entry-vector check itself still reports incompatible for a part
        # that genuinely doesn't match the acceptor (unit-level — bypassing the
        # design's auto-conform): a part carrying AATG/GCTT vs a CTCG/TGAG UPD
        # acceptor won't ligate.
        from splicecraft_agent import _agent_entry_vector_check
        self._bind_vector("CTCG", "TGAG", name="UPD")
        g = sc._BUILTIN_GRAMMARS["gb_l0"]
        chk = _agent_entry_vector_check(g, "gb_l0", "",
                                        {"oh5": "AATG", "oh3": "GCTT"})
        assert chk["checkable"] is True and chk["compatible"] is False, chk

    def test_no_vector_configured_reported(self):
        sc._set_entry_vector("gb_l0", None)
        chk = self._call(check_entry_vector=True)["result"]["entry_vector_check"]
        assert chk["configured"] is False

    def test_check_omitted_leaves_no_key(self):
        self._bind_vector("AATG", "GCTT")
        r = self._call()
        assert r["ok"] and "entry_vector_check" not in r["result"]


class TestProteinMotifEndpoints:
    """Agent endpoints list/set/delete-protein-motif — the merge + builtin-
    override + delete-restores-builtin behaviors. The GUI/synthesis side is
    covered in test_synthesis.py; the agent side was untested."""

    def _builtin_name(self):
        # User file is empty under the sandbox → _load_protein_motifs is just
        # the built-ins.
        motifs = sc._load_protein_motifs()
        assert motifs, "no builtin motifs"
        return motifs[0]["name"]

    def _merged(self):
        return {m["name"]: m
                for m in sc._h_list_protein_motifs(None, {})["motifs"]}

    def test_list_returns_builtins(self):
        r = sc._h_list_protein_motifs(None, {})
        assert len(r["motifs"]) >= 1
        assert all("name" in m and "sequence" in m for m in r["motifs"])

    def test_set_creates_user_motif(self):
        r = sc._h_set_protein_motif(None, {
            "name": "MyTag", "sequence": "gsgsgs", "feature_type": "Tag",
            "color": "#ff0000", "description": "a linker"})
        assert r == {"ok": True, "name": "MyTag"}
        assert self._merged()["MyTag"]["sequence"] == "GSGSGS"   # upper-cased

    def test_set_missing_name_400(self):
        assert sc._h_set_protein_motif(None, {"sequence": "GS"})[1] == 400
        assert sc._h_set_protein_motif(
            None, {"name": "   ", "sequence": "GS"})[1] == 400

    def test_set_missing_sequence_400(self):
        assert sc._h_set_protein_motif(None, {"name": "X"})[1] == 400

    def test_set_noncanonical_aa_400(self):
        r = sc._h_set_protein_motif(None, {"name": "Bad", "sequence": "GSBZ1"})
        assert isinstance(r, tuple) and r[1] == 400

    def test_set_overrides_builtin(self):
        bn = self._builtin_name()
        assert sc._h_set_protein_motif(
            None, {"name": bn, "sequence": "WWWWWW"})["ok"]
        assert self._merged()[bn]["sequence"] == "WWWWWW"       # override wins

    def test_delete_override_restores_builtin(self):
        bn = self._builtin_name()
        orig = self._merged()[bn]["sequence"]
        sc._h_set_protein_motif(None, {"name": bn, "sequence": "WWWWWW"})
        assert sc._h_delete_protein_motif(None, {"name": bn}) == {
            "ok": True, "name": bn}
        assert self._merged()[bn]["sequence"] == orig            # restored

    def test_delete_builtin_without_override_404(self):
        r = sc._h_delete_protein_motif(None, {"name": self._builtin_name()})
        assert isinstance(r, tuple) and r[1] == 404

    def test_delete_missing_name_400(self):
        assert sc._h_delete_protein_motif(None, {})[1] == 400

    def test_set_copy_on_write_keeps_others(self):
        sc._h_set_protein_motif(None, {"name": "A", "sequence": "AA"})
        sc._h_set_protein_motif(None, {"name": "B", "sequence": "CC"})
        sc._h_set_protein_motif(None, {"name": "A", "sequence": "DD"})  # override
        merged = self._merged()
        assert merged["A"]["sequence"] == "DD"
        assert merged["B"]["sequence"] == "CC"


class TestDangerousParamGuard:
    """`_agent_reject_dangerous_unknowns`: a routing/selection param passed to
    an endpoint that doesn't accept it fails loud (the SC-D silent-ok:true
    footgun) instead of being quietly dropped. Benign forward-compat unknowns
    stay soft — the endpoint still runs."""

    def test_move_part_rejects_bin_param(self):
        # move-part routes by to/from, not `bin` — a `bin` key is a mistake.
        r = sc._h_move_part(None, {"to": "B", "name": "p", "bin": "wrong"})
        assert isinstance(r, tuple) and r[1] == 400
        assert "bin" in r[0]["error"]

    def test_move_part_benign_unknown_passes_guard(self):
        # A benign unknown doesn't trip the guard; the call proceeds and 404s on
        # the missing bin (NOT the guard's 400).
        r = sc._h_move_part(None, {"to": "Nope", "name": "p", "notes": "hi"})
        assert isinstance(r, tuple) and r[1] == 404

    def test_move_part_normal_call_unaffected(self):
        r = sc._h_move_part(None, {"to": "Nope", "name": "p"})
        assert isinstance(r, tuple) and r[1] == 404

    def test_move_primer_rejects_collection_typo(self):
        r = sc._h_move_primer(None, {"to": "C", "name": "p",
                                     "source_collection": "x"})
        assert isinstance(r, tuple) and r[1] == 400

    def test_move_experiment_rejects_collection(self):
        r = sc._h_move_experiment(None, {"id": "e1", "to": "P",
                                         "collection": "x"})
        assert isinstance(r, tuple) and r[1] == 400


class TestLintSynthesis:
    """`lint-synthesis` — synthesis/assembly readiness pre-flight (internal
    Type IIS sites, GC extremes, homopolymers, tandem repeats, ORF integrity)."""

    def _lint(self, seq, **extra):
        return sc._h_lint_synthesis(None, dict({"sequence": seq}, **extra))

    def test_flags_internal_bsai_site(self):
        r = self._lint("ATCGATCGAT" + "GGTCTC" + "ATCGATCGATTAGCTAGCTAGC")
        assert r["ok"]
        assert any(w["kind"] == "type_iis_site" and w["level"] == "error"
                   for w in r["warnings"])
        assert r["score"] < 100

    def test_clean_mcs_no_type_iis(self):
        # HindIII/BamHI/EcoRI/XbaI/SmaI MCS — none are BsaI/Esp3I/BbsI/SapI.
        r = self._lint("AAGCTTGGATCCGAATTCTCTAGACCCGGG")
        assert not any(w["kind"] == "type_iis_site" for w in r["warnings"])

    def test_flags_homopolymer(self):
        r = self._lint("ATCGATCG" + "A" * 15 + "GCTAGCTA")
        assert any(w["kind"] == "homopolymer" for w in r["warnings"])

    def test_flags_extreme_gc(self):
        r = self._lint("GC" * 60)
        assert any(w["kind"] == "gc_overall" for w in r["warnings"])
        assert r["stats"]["gc"] > 0.9

    def test_expect_cds_clean_orf_ok(self):
        r = self._lint("ATG" + "GCTACG" * 16 + "TAA", expect_cds=True,
                       circular=False)
        assert not any(w["kind"] == "orf_integrity" for w in r["warnings"])

    def test_expect_cds_no_orf_errors(self):
        r = self._lint("TAA" * 20, expect_cds=True, circular=False)
        assert any(w["kind"] == "orf_integrity" and w["level"] == "error"
                   for w in r["warnings"])

    def test_forbidden_enzymes_override(self):
        seq = "ATCGAT" + "GAATTC" + "ATCGATCGTAGCTAGC"      # EcoRI site
        assert not any(w["kind"] == "type_iis_site"
                       for w in self._lint(seq)["warnings"])          # default
        r = self._lint(seq, forbidden_enzymes=["EcoRI"])
        assert any(w["kind"] == "type_iis_site" for w in r["warnings"])

    def test_resolves_by_library_entry(self, isolated_library):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATCGATCG" + "GGTCTC" + "ATCGATCGTAGC"),
                        id="LINTME", name="LINTME",
                        annotations={"molecule_type": "DNA",
                                       "topology": "circular"})
        sc._save_library([{"id": "LINTME", "name": "Lint Me", "size": 0,
                           "n_feats": 0, "added": "2026-07-03",
                           "gb_text": sc._record_to_gb_text(rec)}])
        r = sc._h_lint_synthesis(None, {"name": "Lint Me"})
        assert r["ok"] and r["stats"]["circular"] is True
        assert any(w["kind"] == "type_iis_site" for w in r["warnings"])

    def test_missing_input_400(self):
        r = sc._h_lint_synthesis(None, {})
        assert isinstance(r, tuple) and r[1] == 400

    def test_unknown_entry_404(self, isolated_library):
        sc._save_library([])
        r = sc._h_lint_synthesis(None, {"id": "ghost"})
        assert isinstance(r, tuple) and r[1] == 404

    def test_bad_forbidden_enzymes_type_400(self):
        r = sc._h_lint_synthesis(None, {"sequence": "ATCG",
                                        "forbidden_enzymes": "BsaI"})
        assert isinstance(r, tuple) and r[1] == 400

    def test_length_cap_413(self):
        r = self._lint("A" * 200_001)
        assert isinstance(r, tuple) and r[1] == 413

    def test_empty_or_invalid_sequence_rejected(self):
        # A sequence with no valid IUPAC bases is rejected, not linted as empty.
        r = sc._h_lint_synthesis(None, {"sequence": "123456"})
        assert isinstance(r, tuple) and r[1] in (400, 422)

    def test_too_many_forbidden_enzymes_400(self):
        r = self._lint("ATCGATCG", forbidden_enzymes=["BsaI"] * 33)
        assert isinstance(r, tuple) and r[1] == 400

    def test_unknown_forbidden_enzyme_400(self):
        # A safety check must fail loud on a typo, not silently skip the site.
        r = self._lint("ATCGATCG", forbidden_enzymes=["BsaI", "Bsal"])
        assert isinstance(r, tuple) and r[1] == 400

    def test_lint_synthesis_is_heavy_classified(self):
        assert "lint-synthesis" in sc._AGENT_HEAVY_ENDPOINTS


class TestVerifyAgainstReads:
    """`verify-against-reads` — align sequencing reads to a designed construct
    and report per-read identity + a match/mismatch verdict."""

    REF = "ATGGCAAGCGGTGGTTCTGGTAGCGGTAGCGGTAGCGGTAGCGGTACCGATCGATCGTAA"  # 60

    def _verify(self, reads, **extra):
        return sc._h_verify_against_reads(
            None, dict({"reference": self.REF, "reads": reads,
                        "circular": False}, **extra))

    def _one_mismatch(self):
        i = 30
        return self.REF[:i] + ("T" if self.REF[i] != "T" else "A") + self.REF[i + 1:]

    def test_perfect_read_matches(self):
        r = self._verify([self.REF])
        assert r["ok"] and r["verdict"] == "match"
        assert r["reads"][0]["identity_pct"] == 100.0
        assert r["reads"][0]["passes"] is True

    def test_mismatch_read_flagged(self):
        r = self._verify([self._one_mismatch()])
        assert r["verdict"] == "mismatch"
        assert r["reads"][0]["passes"] is False and r["summary"]["n_fail"] == 1

    def test_rc_read_still_matches(self):
        r = self._verify([sc._rc(self.REF)])
        assert r["reads"][0]["identity_pct"] >= 99.0
        assert r["reads"][0]["rc"] is True

    def test_min_identity_threshold(self):
        # The same single-mismatch read passes under a lenient threshold.
        assert self._verify([self._one_mismatch()],
                            min_identity=90.0)["verdict"] == "match"

    def test_mixed_reads_verdict(self):
        r = self._verify([self.REF, self._one_mismatch()])
        assert r["n_reads"] == 2 and r["verdict"] == "mismatch"
        assert r["summary"]["n_pass"] == 1 and r["summary"]["n_fail"] == 1

    def test_resolve_reference_by_entry(self, isolated_library):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq(self.REF), id="REFME", name="REFME",
                        annotations={"molecule_type": "DNA",
                                       "topology": "linear"})
        sc._save_library([{"id": "REFME", "name": "Ref Me", "size": 0,
                           "n_feats": 0, "added": "2026-07-03",
                           "gb_text": sc._record_to_gb_text(rec)}])
        r = sc._h_verify_against_reads(
            None, {"reference_name": "Ref Me", "reads": [self.REF]})
        assert r["ok"] and r["reads"][0]["identity_pct"] == 100.0

    def test_missing_reference_400(self):
        r = sc._h_verify_against_reads(None, {"reads": ["ACGT"]})
        assert isinstance(r, tuple) and r[1] == 400

    def test_missing_reads_400(self):
        r = sc._h_verify_against_reads(None, {"reference": self.REF})
        assert isinstance(r, tuple) and r[1] == 400

    def test_unknown_reference_404(self, isolated_library):
        sc._save_library([])
        r = sc._h_verify_against_reads(
            None, {"reference_id": "ghost", "reads": [self.REF]})
        assert isinstance(r, tuple) and r[1] == 404

    def test_too_many_reads_413(self):
        r = self._verify([self.REF] * 201)
        assert isinstance(r, tuple) and r[1] == 413


class TestRequestDispatcherHardening:
    """The HTTP dispatcher must hand handlers a real dict (never None,
    never a list) so .get() never crashes."""

    def test_handle_passes_dict_on_empty_body(self, http_server):
        base, token, _ = http_server
        # POST with Content-Length: 0 → handler should still get {}
        import urllib.request, urllib.error
        req = urllib.request.Request(
            f"{base}/save", data=b"",
            headers={"Authorization": f"Bearer {token}"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            status = 200
        except urllib.error.HTTPError as exc:
            status = exc.code
        # 200 (saved), 422 (no record), or 409 (the homeless-new-record
        # create gate, snag #16) are all valid non-crash responses; what
        # we care about is that the handler didn't AttributeError on a
        # None body.
        assert status in (200, 422, 409)

    def test_handle_rejects_non_dict_json(self, http_server):
        # Sweep #25 (2026-05-23): non-dict / malformed JSON body now
        # 400s explicitly (was: silently normalised to {} which
        # masked caller serialisation bugs). The handler still can't
        # 500 from a `.get()` on a list — the dispatcher catches the
        # bad shape before reaching the handler.
        base, token, _ = http_server
        import urllib.request, urllib.error
        req = urllib.request.Request(
            f"{base}/save", data=b'[1,2,3]',
            headers={"Authorization": f"Bearer {token}",
                       "Content-Type": "application/json"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            status = 200
        except urllib.error.HTTPError as exc:
            status = exc.code
        # Either accepted (200) or rejected (400/422) — never a 500.
        assert status in (200, 400, 422)


# ── Plasmid status endpoints (added 2026-05-05 for v1.0) ──────────────────────


class TestPlasmidStatusEndpoints:
    def test_list_plasmid_statuses(self, http_server):
        base, token, _ = http_server
        status, payload = _http(f"{base}/list-plasmid-statuses", token=token)
        assert status == 200
        assert payload["ok"] is True
        # Strict canonical vocabulary — DESIGNING / CLONING /
        # SEQUENCING / VERIFIED / ERROR (the last added in v0.9.24
        # for failed-clone tracking, INV-76).
        assert set(payload["statuses"]) == {
            "DESIGNING", "CLONING", "SEQUENCING", "VERIFIED", "ERROR"
        }
        # Each status carries a hex color; the agent can use it for
        # rendering without re-deriving from the GUI.
        assert all(c.startswith("#") for c in payload["colors"].values())

    def test_set_plasmid_status_round_trip(self, http_server, tiny_record):
        # Seed one library entry the endpoint can target.
        sc._save_library([{"name": "pTest", "id": "pTest",
                            "gb_text": "fake"}])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-plasmid-status", method="POST",
            body={"name": "pTest", "status": "CLONING"}, token=token,
        )
        assert status == 200
        assert payload["status"] == "CLONING"
        # Persisted on disk.
        entry = next(e for e in sc._load_library() if e["name"] == "pTest")
        assert entry["status"] == "CLONING"

    def test_set_plasmid_status_clears_with_empty_string(self, http_server):
        sc._save_library([{"name": "pTest", "id": "pTest",
                            "status": "VERIFIED", "gb_text": "fake"}])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-plasmid-status", method="POST",
            body={"name": "pTest", "status": ""}, token=token,
        )
        assert status == 200
        assert payload["status"] == ""

    def test_set_plasmid_status_invalid_collapses_to_empty(self, http_server):
        """Per `_sanitize_plasmid_status`'s strict-canonical-or-empty
        contract: a non-canonical string (mixed case, garbage)
        silently degrades to "" rather than 400. Documented behaviour
        — the round-trip-exact rule for hand-edited library JSON."""
        sc._save_library([{"name": "pTest", "id": "pTest",
                            "gb_text": "fake"}])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-plasmid-status", method="POST",
            body={"name": "pTest", "status": "Designing"},  # mixed case
            token=token,
        )
        assert status == 200
        assert payload["status"] == ""

    def test_set_plasmid_status_unknown_name_404(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-plasmid-status", method="POST",
            body={"name": "ghost", "status": "DESIGNING"}, token=token,
        )
        assert status == 404
        assert "ghost" in payload["error"]

    def test_set_plasmid_status_rejects_non_string(self, http_server):
        sc._save_library([{"name": "pTest", "id": "pTest",
                            "gb_text": "fake"}])
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-plasmid-status", method="POST",
            body={"name": "pTest", "status": 42}, token=token,
        )
        assert status == 400
        assert "string" in payload["error"]


# ── Entry-vector endpoints (added 2026-05-05 for v1.0) ───────────────────────


def _minimal_gb_text() -> str:
    """Smallest GenBank text that round-trips through SeqIO — used
    so set-entry-vector's parse-validate step has something real to
    chew on without fixture sprawl. Column widths match SeqIO's own
    LOCUS-line formatter so Biopython parses it without warning."""
    return ("LOCUS       test                      10 bp    DNA     "
            "circular SYN 01-JAN-2026\n"
            "FEATURES             Location/Qualifiers\n"
            "ORIGIN      \n"
            "        1 atgcatgcat\n"
            "//\n")


class TestEntryVectorEndpoints:
    def test_list_entry_vectors_empty(self, http_server):
        base, token, _ = http_server
        status, payload = _http(f"{base}/list-entry-vectors", token=token)
        assert status == 200
        assert payload["ok"] is True
        assert payload["entry_vectors"] == []

    def test_set_get_entry_vector_round_trip(self, http_server):
        base, token, _ = http_server
        gb = _minimal_gb_text()
        # SET
        status, payload = _http(
            f"{base}/set-entry-vector", method="POST",
            body={"grammar_id": "gb_l0", "name": "pUPD2",
                   "gb_text": gb, "source": "library:test"},
            token=token,
        )
        assert status == 200, payload
        assert payload["vector"]["name"] == "pUPD2"
        assert payload["vector"]["size"] == 10
        # The set response strips `gb_text` to keep responses small.
        assert "gb_text" not in payload["vector"]
        # GET
        status, payload = _http(
            f"{base}/get-entry-vector", method="POST",
            body={"grammar_id": "gb_l0"}, token=token,
        )
        assert status == 200
        assert payload["vector"]["name"]    == "pUPD2"
        assert payload["vector"]["gb_text"] == gb

    def test_get_entry_vector_returns_null_when_unset(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/get-entry-vector", method="POST",
            body={"grammar_id": "moclo_plant"}, token=token,
        )
        assert status == 200
        assert payload["vector"] is None

    def test_set_entry_vector_clear(self, http_server):
        base, token, _ = http_server
        gb = _minimal_gb_text()
        _http(f"{base}/set-entry-vector", method="POST",
              body={"grammar_id": "gb_l0", "name": "pUPD2", "gb_text": gb},
              token=token)
        status, payload = _http(
            f"{base}/set-entry-vector", method="POST",
            body={"grammar_id": "gb_l0", "clear": True}, token=token,
        )
        assert status == 200
        assert payload["vector"] is None
        assert sc._get_entry_vector("gb_l0") is None

    def test_set_entry_vector_invalid_gb_text(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-entry-vector", method="POST",
            body={"grammar_id": "gb_l0", "name": "x",
                   "gb_text": "not a genbank file"},
            token=token,
        )
        assert status == 400
        assert "parse failed" in payload["error"]

    def test_set_entry_vector_oversized_gb_text(self, http_server):
        base, token, _ = http_server
        # 600 KB of fake bases — over the inner 500 KB cap but under
        # the HTTP transport's 1 MiB body cap, so the inner check is
        # the one that fires.
        big = "A" * (600 * 1024)
        status, payload = _http(
            f"{base}/set-entry-vector", method="POST",
            body={"grammar_id": "gb_l0", "name": "x", "gb_text": big},
            token=token,
        )
        assert status == 400
        assert "too large" in payload["error"]

    def test_set_entry_vector_missing_grammar_id(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-entry-vector", method="POST",
            body={"name": "x", "gb_text": _minimal_gb_text()}, token=token,
        )
        assert status == 400
        assert "grammar_id" in payload["error"]


# ── update-primer endpoint (added 2026-05-05 for v1.0) ───────────────────────


class TestUpdatePrimerEndpoint:
    def test_update_primer_rejects_non_primer_feature(self, http_server,
                                                       tiny_record):
        """The endpoint MUST refuse to mutate a non-primer feature so
        an agent can't smuggle a primer-only field (e.g. `primer_seq`)
        onto a CDS or misc_feature. tiny_record's idx 0 is a CDS."""
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/update-primer", method="POST",
            body={"idx": 0, "label": "x"}, token=token,
        )
        assert status == 400
        assert "primer_bind" in payload["error"]

    def test_update_primer_validates_idx_out_of_range(self, http_server,
                                                       tiny_record):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/update-primer", method="POST",
            body={"idx": 99}, token=token,
        )
        assert status == 400
        assert "out of range" in payload["error"]

    def test_update_primer_rejects_infinity_idx(self, http_server,
                                                  tiny_record):
        base, token, _ = http_server
        body_json = '{"idx": Infinity, "label": "x"}'
        req = urllib.request.Request(
            f"{base}/update-primer", data=body_json.encode(),
            headers={"Authorization": f"Bearer {token}",
                       "Content-Type": "application/json"},
            method="POST",
        )
        try:
            urllib.request.urlopen(req, timeout=5).read()
            code = 200
        except urllib.error.HTTPError as exc:
            code = exc.code
        assert code == 400

    def test_update_primer_rejects_oversized_primer_seq(self, http_server,
                                                         tiny_record):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/update-primer", method="POST",
            body={"idx": 0, "primer_seq": "A" * 600}, token=token,
        )
        assert status == 400
        # Either the non-primer reject (if idx 0 is non-primer) or the
        # length cap. Both are correct rejections.
        assert ("too long" in payload["error"]
                or "primer_bind" in payload["error"])


# ── Settings endpoints (added 2026-05-05 for v1.0) ───────────────────────────


class TestSettingsEndpoints:
    def test_get_settings_returns_allowlisted_keys(self, http_server):
        base, token, _ = http_server
        status, payload = _http(f"{base}/get-settings", token=token)
        assert status == 200
        # Spot-check: every allowlisted key is present, infrastructure
        # keys are not.
        keys = set(payload["settings"].keys())
        for required in ("show_feature_tooltips", "min_primer_binding",
                          "active_grammar"):
            assert required in keys
        for excluded in ("last_known_latest", "last_seen_version",
                          "last_update_check_ts", "hmm_db_path"):
            assert excluded not in keys

    def test_set_setting_round_trip_bool(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-setting", method="POST",
            body={"key": "click_debug", "value": True}, token=token,
        )
        assert status == 200
        assert payload["value"] is True
        assert sc._get_setting("click_debug") is True

    def test_set_setting_round_trip_int_range(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-setting", method="POST",
            body={"key": "min_primer_binding", "value": 18}, token=token,
        )
        assert status == 200
        assert payload["value"] == 18
        assert sc._get_setting("min_primer_binding") == 18

    def test_set_setting_int_range_rejects_out_of_range(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-setting", method="POST",
            body={"key": "min_primer_binding", "value": 100}, token=token,
        )
        assert status == 400
        assert "[1, 60]" in payload["error"]

    def test_set_setting_unknown_key_after_linear_layout_removed(
            self, http_server):
        # `linear_layout` was removed from the allowlist 2026-05-08
        # — flag is the only linear layout. Setting it through the
        # agent now returns an "unknown key" error, NOT the choice-
        # validator's "must be one of …" error.
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-setting", method="POST",
            body={"key": "linear_layout", "value": "flag"}, token=token,
        )
        assert status == 400
        assert "unknown" in payload["error"].lower()

    def test_set_setting_bool_rejects_string(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-setting", method="POST",
            body={"key": "show_restr", "value": "true"}, token=token,
        )
        assert status == 400
        assert "boolean" in payload["error"]

    def test_set_setting_unknown_key(self, http_server):
        base, token, _ = http_server
        status, payload = _http(
            f"{base}/set-setting", method="POST",
            body={"key": "secret_setting", "value": "boom"}, token=token,
        )
        assert status == 400
        assert "unknown setting" in payload["error"]
        # Helpfully lists what the agent CAN write.
        assert "min_primer_binding" in payload["available"]

    def test_set_setting_restr_min_len_only_accepts_4_or_6(self, http_server):
        base, token, _ = http_server
        for good in (4, 6):
            status, _ = _http(
                f"{base}/set-setting", method="POST",
                body={"key": "restr_min_len", "value": good}, token=token,
            )
            assert status == 200
        for bad in (5, 8, 0):
            status, payload = _http(
                f"{base}/set-setting", method="POST",
                body={"key": "restr_min_len", "value": bad}, token=token,
            )
            assert status == 400, (bad, payload)


# ── Simulator agent endpoints (2026-05-17 release) ────────────────────────────
#
# `simulate-pcr` and `simulate-gel` are pure read-only wrappers around
# the SimulatorScreen's underlying functions. Tested at the handler
# layer (no HTTP round-trip) — input validation is the bulk of the
# surface; the underlying physics is covered by `tests/test_simulator.py`.


class TestSimulatePcrHandler:

    _SEQ = ("ATGCGATCGATCGATCGCGT"   # fwd binding site 0..20
            + "A" * 60
            + "GCATCGTAGCTAGCTGATCG") # rev-rc binding site 80..100
    _FWD = "ATGCGATCGATCGATCGCGT"
    _REV = "CGATCAGCTAGCTACGATGC"     # = rc("GCATCGTAGCTAGCTGATCG")

    def test_happy_path_linear(self):
        resp = sc._h_simulate_pcr(None, {
            "template_seq": self._SEQ,
            "fwd_primer":   self._FWD,
            "rev_primer":   self._REV,
            "circular":     False,
        })
        assert isinstance(resp, dict)
        assert resp["ok"] is True
        assert resp["n"] == 1
        assert resp["capped"] is False
        assert resp["amplicons"][0]["length"] == 100
        assert resp["amplicons"][0]["wraps"] is False

    def test_circular_wrap_amplicon(self):
        # Place fwd-binding-site near end, rev-binding-site near start;
        # amplicon must cross the origin.
        # seq pos 20..40: ATGCGATCGATCGATCGCGT  (the rev-binding target)
        # seq pos 50..70: GCATCGTAGCTAGCTGATCG  (the fwd-binding site)
        seq = ("A" * 20 + self._FWD + "A" * 10
                + "GCATCGTAGCTAGCTGATCG" + "A" * 30)
        resp = sc._h_simulate_pcr(None, {
            "template_seq": seq,
            "fwd_primer":   "GCATCGTAGCTAGCTGATCG",
            "rev_primer":   sc._rc(self._FWD),
            "circular":     True,
            "max_amplicon": 200,
        })
        assert isinstance(resp, dict) and resp["ok"] is True
        wrap_amps = [a for a in resp["amplicons"] if a["wraps"]]
        assert wrap_amps, "expected a wrapping amplicon"

    def test_no_match_returns_empty(self):
        resp = sc._h_simulate_pcr(None, {
            "template_seq": "ATGC" * 100,
            "fwd_primer":   "AAAAAAAAAAAAAAAA",
            "rev_primer":   "TTTTTTTTTTTTTTTT",
        })
        assert resp["ok"] is True
        assert resp["n"] == 0
        assert resp["amplicons"] == []

    def test_missing_template_seq_returns_400(self):
        payload, status = sc._h_simulate_pcr(None, {
            "fwd_primer": self._FWD, "rev_primer": self._REV,
        })
        assert status == 400
        assert "template_seq" in payload["error"]

    def test_non_string_template_returns_400(self):
        payload, status = sc._h_simulate_pcr(None, {
            "template_seq": 123,
            "fwd_primer":   self._FWD,
            "rev_primer":   self._REV,
        })
        assert status == 400

    def test_template_over_cap_returns_413(self):
        payload, status = sc._h_simulate_pcr(None, {
            "template_seq": "A" * (sc._PCR_MAX_TEMPLATE_BP + 1),
            "fwd_primer":   self._FWD,
            "rev_primer":   self._REV,
        })
        assert status == 413
        assert "template_seq" in payload["error"]

    def test_missing_primer_returns_400(self):
        payload, status = sc._h_simulate_pcr(None, {
            "template_seq": self._SEQ,
            "fwd_primer":   self._FWD,
        })
        assert status == 400

    def test_short_primer_returns_400(self):
        payload, status = sc._h_simulate_pcr(None, {
            "template_seq": self._SEQ,
            "fwd_primer":   "ATGCG",
            "rev_primer":   self._REV,
        })
        assert status == 400
        assert "at least" in payload["error"]

    def test_long_primer_returns_400(self):
        payload, status = sc._h_simulate_pcr(None, {
            "template_seq": self._SEQ,
            "fwd_primer":   "A" * (sc._PCR_MAX_PRIMER_LEN + 1),
            "rev_primer":   self._REV,
        })
        assert status == 400
        assert "at most" in payload["error"]

    def test_non_acgt_primer_returns_400(self):
        for bad in ("NNNNNNNNNNNNNNN", "ATGCGATCGAT-GATCG", "atgcga"):
            payload, status = sc._h_simulate_pcr(None, {
                "template_seq": self._SEQ,
                "fwd_primer":   bad if not bad.islower() else bad.upper() + "X",
                "rev_primer":   self._REV,
            })
            assert status == 400, (bad, payload)

    def test_max_amplicon_out_of_range_returns_400(self):
        for bad in (0, -1, sc._PCR_AMPLICON_HARD_CAP + 1):
            payload, status = sc._h_simulate_pcr(None, {
                "template_seq": self._SEQ,
                "fwd_primer":   self._FWD,
                "rev_primer":   self._REV,
                "max_amplicon": bad,
            })
            assert status == 400, (bad, payload)

    def test_max_amplicon_non_int_returns_400(self):
        payload, status = sc._h_simulate_pcr(None, {
            "template_seq": self._SEQ,
            "fwd_primer":   self._FWD,
            "rev_primer":   self._REV,
            "max_amplicon": "not-an-int",
        })
        assert status == 400

    def test_empty_primer_strings_return_400(self):
        for fwd, rev in [("", self._REV), (self._FWD, ""), ("   ", "  ")]:
            payload, status = sc._h_simulate_pcr(None, {
                "template_seq": self._SEQ,
                "fwd_primer":   fwd,
                "rev_primer":   rev,
            })
            assert status == 400, (fwd, rev, payload)


class TestSimulateGelHandler:

    def test_happy_path_ladder(self):
        resp = sc._h_simulate_gel(None, {
            "lanes": [{"source": "ladder", "detail": "1 kb"}],
            "agarose_pct": 1.0,
        })
        assert isinstance(resp, dict)
        assert resp["ok"] is True
        assert len(resp["lanes"]) == 1
        assert len(resp["lanes"][0]["bands"]) > 0
        # Every band has bp + form + mobility + row.
        for b in resp["lanes"][0]["bands"]:
            assert {"bp", "form", "mobility", "row"} <= set(b.keys())
            assert 0.0 <= b["mobility"] <= 1.0
            assert 0 <= b["row"] < resp["height"]

    def test_plasmid_lane_with_circular_template(self):
        resp = sc._h_simulate_gel(None, {
            "lanes": [{"source": "plasmid", "detail": ""}],
            "template_seq": "AT" * 1500,
            "template_circular": True,
            "agarose_pct": 1.0,
        })
        assert resp["ok"] is True
        # Circular uncut → SC + nicked = 2 bands.
        assert len(resp["lanes"][0]["bands"]) == 2
        forms = {b["form"] for b in resp["lanes"][0]["bands"]}
        assert "supercoiled" in forms
        assert "nicked" in forms

    def test_digest_lane(self):
        resp = sc._h_simulate_gel(None, {
            "lanes": [{"source": "digest", "detail": "EcoRI"}],
            "template_seq": "GAATTC" + "A" * 100 + "GAATTC" + "A" * 50,
            "template_circular": True,
            "agarose_pct": 1.0,
        })
        assert resp["ok"] is True
        # Two EcoRI sites on a circular template → 2 fragments.
        assert len(resp["lanes"][0]["bands"]) >= 2

    def test_pcr_lane_with_amplicon(self):
        resp = sc._h_simulate_gel(None, {
            "lanes": [{"source": "pcr", "detail": ""}],
            "pcr_amplicon": {"length": 800, "wraps": False,
                              "amplicon_seq": "A" * 800,
                              "start": 0, "end": 800,
                              "fwd_seq": "A" * 20, "rev_seq": "T" * 20,
                              "gc_pct": 0.0, "fwd_tm": None,
                              "rev_tm": None},
        })
        assert resp["ok"] is True
        assert len(resp["lanes"][0]["bands"]) == 1
        assert resp["lanes"][0]["bands"][0]["bp"] == 800

    def test_include_image_returns_text(self):
        resp = sc._h_simulate_gel(None, {
            "lanes": [{"source": "ladder", "detail": "1 kb"}],
            "include_image": True,
            "height": 10, "lane_width": 5,
        })
        assert resp["ok"] is True
        assert "image" in resp
        assert isinstance(resp["image"], str)
        assert "\n" in resp["image"]   # multi-row rendering

    def test_missing_lanes_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {})
        assert status == 400

    def test_empty_lanes_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {"lanes": []})
        assert status == 400

    def test_lanes_not_list_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": {"source": "ladder"},   # dict, not list
        })
        assert status == 400

    def test_too_many_lanes_returns_400(self):
        lanes = [{"source": "empty"}] * (sc._GEL_MAX_LANES + 1)
        payload, status = sc._h_simulate_gel(None, {"lanes": lanes})
        assert status == 400

    def test_lane_missing_source_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": [{"name": "no-source"}],
        })
        assert status == 400
        assert "source" in payload["error"]

    def test_lane_unknown_source_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": [{"source": "gibberish"}],
        })
        assert status == 400

    def test_lane_non_dict_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": ["not-a-dict"],
        })
        assert status == 400

    def test_lane_detail_too_long_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": [{"source": "digest", "detail": "X" * 300}],
        })
        assert status == 400

    def test_lane_detail_wrong_type_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": [{"source": "ladder", "detail": 42}],
        })
        assert status == 400

    def test_agarose_out_of_range_returns_400(self):
        for bad in (0, 0.05, 11.0, -1.0):
            payload, status = sc._h_simulate_gel(None, {
                "lanes": [{"source": "ladder", "detail": "1 kb"}],
                "agarose_pct": bad,
            })
            assert status == 400, (bad, payload)

    def test_agarose_non_numeric_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": [{"source": "ladder", "detail": "1 kb"}],
            "agarose_pct": "high",
        })
        assert status == 400

    def test_height_out_of_range_returns_400(self):
        for bad in (sc._GEL_HEIGHT_MIN - 1, sc._GEL_HEIGHT_MAX + 1, 0, -10):
            payload, status = sc._h_simulate_gel(None, {
                "lanes": [{"source": "ladder", "detail": "1 kb"}],
                "height": bad,
            })
            assert status == 400, (bad, payload)

    def test_lane_width_out_of_range_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": [{"source": "ladder", "detail": "1 kb"}],
            "lane_width": sc._GEL_LANE_WIDTH_MAX + 1,
        })
        assert status == 400

    def test_template_over_cap_returns_413(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": [{"source": "plasmid", "detail": ""}],
            "template_seq": "A" * (sc._PCR_MAX_TEMPLATE_BP + 1),
        })
        assert status == 413

    def test_template_wrong_type_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": [{"source": "plasmid"}],
            "template_seq": 12345,
        })
        assert status == 400

    def test_pcr_amplicon_wrong_type_returns_400(self):
        payload, status = sc._h_simulate_gel(None, {
            "lanes": [{"source": "pcr"}],
            "pcr_amplicon": "not-a-dict",
        })
        assert status == 400

    def test_pcr_lane_without_amplicon_returns_empty_bands(self):
        # Not a validation error — gel renders a lane with no bands
        # and the user sees an empty column. Mirrors UI behaviour.
        resp = sc._h_simulate_gel(None, {
            "lanes": [{"source": "pcr"}],
        })
        assert resp["ok"] is True
        assert resp["lanes"][0]["bands"] == []


class TestSimulatorAgentRegistration:
    """Both new endpoints must be registered as READ-ONLY (write=False)
    so an unauthenticated caller can run simulations without a token —
    matches `simulate-gibson` semantics."""

    def test_simulate_pcr_registered_read_only(self):
        eps = {ep["name"]: ep for ep in sc._h_tools(None, {})["endpoints"]}
        assert "simulate-pcr" in eps
        assert eps["simulate-pcr"]["write"] is False

    def test_simulate_gel_registered_read_only(self):
        eps = {ep["name"]: ep for ep in sc._h_tools(None, {})["endpoints"]}
        assert "simulate-gel" in eps
        assert eps["simulate-gel"]["write"] is False


class TestAddCodonTableGenome:
    """`add-codon-table` source='genome' branch. The NCBI Datasets fetch is
    monkeypatched (no network); validation + save wiring are exercised."""

    def test_missing_accession_and_taxid(self):
        result = sc._h_add_codon_table(None, {"source": "genome"})
        payload, status = result
        assert status == 400
        assert "accession" in payload["error"]

    def test_bad_mode_rejected(self):
        result = sc._h_add_codon_table(
            None, {"source": "genome", "taxid": "1423", "mode": "best"})
        payload, status = result
        assert status == 400
        assert "mode" in payload["error"]

    def test_happy_path_builds_and_saves(self, monkeypatch):
        fake_raw = {"GCT": ("A", 100), "ATG": ("M", 30), "TAA": ("*", 5)}

        def fake_build(query, mode, timeout=60.0):
            return fake_raw, "built ok", {
                "accession": "GCF_TEST.1",
                "taxid": query if str(query).isdigit() else "",
                "organism": "Testus organismus",
                "stats": {"mode": mode, "n_cds_total": 7, "n_codons": 135},
            }
        # _h_add_codon_table moved to splicecraft_agent — patch the sibling
        # namespace it resolves the genome builder in (not the hub re-export).
        monkeypatch.setattr(
            "splicecraft_agent._genome_build_codon_table", fake_build)
        result = sc._h_add_codon_table(
            None, {"source": "genome", "taxid": "1423", "mode": "heg"})
        assert result["ok"] is True
        assert result["entry"]["source"] == "genome"
        assert result["entry"]["taxid"] == "1423"
        got = sc._codon_tables_get("1423")
        assert got is not None and got["source"] == "genome"
        assert got["name"] == "Testus organismus"   # default from organism

    def test_build_failure_returns_502(self, monkeypatch):
        monkeypatch.setattr(
            "splicecraft_agent._genome_build_codon_table",
            lambda q, m, timeout=60.0: (None, "no such assembly", None))
        result = sc._h_add_codon_table(
            None, {"source": "genome", "accession": "GCF_000000000.0"})
        payload, status = result
        assert status == 502
        assert "no such assembly" in payload["error"]


class TestAddCodonTableFile:
    """`add-codon-table` source='file' branch — offline build from a local CDS
    FASTA on disk. Exercises path-safety, validation, and the real save (no
    network); the builder runs for real against an inline CDS file."""

    _CDS = (">lcl|X_cds_1 [gene=rplB] [protein=50S ribosomal protein L2]\n"
            "ATGAAAGCTCGTTGGTGT\n"
            ">lcl|X_cds_2 [gene=dnaA] [protein=replication initiator]\n"
            "ATGTGGGCTAAA\n")

    def test_missing_path_rejected(self):
        payload, status = sc._h_add_codon_table(None, {"source": "file"})
        assert status == 400 and "path" in payload["error"]

    def test_bad_mode_rejected(self, tmp_path):
        p = tmp_path / "cds.fna"
        p.write_text(self._CDS)
        payload, status = sc._h_add_codon_table(
            None, {"source": "file", "path": str(p), "mode": "best"})
        assert status == 400 and "mode" in payload["error"]

    def test_happy_path_builds_and_saves(self, tmp_path):
        p = tmp_path / "ecoli_cds.fna"
        p.write_text(self._CDS)
        result = sc._h_add_codon_table(
            None, {"source": "file", "path": str(p), "mode": "genome"})
        assert result["ok"] is True
        assert result["entry"]["source"] == "file"
        assert result["entry"]["name"] == "ecoli_cds"   # stem -> organism default
        key = result["entry"]["taxid"] or result["entry"]["name"]
        got = sc._codon_tables_get(key)
        assert got is not None and got["source"] == "file"

    def test_nonexistent_file_returns_4xx(self, tmp_path):
        payload, status = sc._h_add_codon_table(
            None, {"source": "file", "path": str(tmp_path / "nope.fna"),
                   "mode": "genome"})
        assert status in (400, 502)   # build returns no data / clean read error


class TestAgentRegistryIntegrity:
    """Guard the agent-endpoint registry as a whole (audit 2026-06-30).

    The endpoint count had silently drifted (docs said ~135; the registry
    held 172) because nothing pinned it. These checks catch an accidental
    mass-unregistration (e.g. a sibling that stops importing its `_h_*`
    handlers) and keep the `/tools` self-describe usable — every handler
    must carry a docstring so an agent can form a request from `doc_full`."""

    def test_registry_has_full_surface(self):
        # A FLOOR, not an exact pin (endpoints are added often). Dropping
        # below this means a whole sibling failed to register — a real bug,
        # not a routine addition.
        assert len(sc._AGENT_HANDLERS) >= 160, (
            f"agent registry shrank to {len(sc._AGENT_HANDLERS)} — did a "
            f"sibling stop importing / registering its _h_* handlers?")

    def test_every_handler_is_wellformed_and_documented(self):
        for name, entry in sc._AGENT_HANDLERS.items():
            assert isinstance(name, str) and name, f"bad endpoint name {name!r}"
            assert isinstance(entry, tuple) and len(entry) == 2, (
                f"{name}: registry value must be (fn, write)")
            fn, write = entry
            assert callable(fn), f"{name}: handler not callable"
            assert isinstance(write, bool), f"{name}: write flag not bool"
            # /tools emits doc_full from the docstring — a missing one makes
            # the endpoint undiscoverable to an agent forming a request.
            assert (fn.__doc__ or "").strip(), (
                f"{name}: handler has no docstring (breaks /tools doc_full)")

    def test_tools_endpoint_describes_every_handler(self):
        tools_fn = sc._AGENT_HANDLERS["tools"][0]
        out = tools_fn(None, {})
        described = {e["name"] for e in out["endpoints"]}
        assert described == set(sc._AGENT_HANDLERS), (
            "tools self-describe is out of sync with the registry")

    def test_heavy_endpoint_set_is_registered(self):
        # A typo in _AGENT_HEAVY_ENDPOINTS would silently fail to gate a heavy
        # op behind the concurrency semaphore (audit 2026-06-30, MED-2a).
        unknown = sc._AGENT_HEAVY_ENDPOINTS - set(sc._AGENT_HANDLERS)
        assert not unknown, (
            f"_AGENT_HEAVY_ENDPOINTS names not in the registry: {unknown}")

    def test_refresh_map_names_are_registered_writes(self):
        # Every live-refresh entry must name a real WRITE endpoint, else a
        # rename would silently stop the open TUI from updating.
        for name in sc._AGENT_REFRESH_MAP:
            entry = sc._AGENT_HANDLERS.get(name)
            assert entry is not None, f"refresh map names unknown endpoint {name!r}"
            assert entry[1] is True, f"refresh map names a READ endpoint {name!r}"


class TestAgentInvokeCore:
    """The shared in-process dispatch core (audit 2026-06-30, Phase 2) that
    gives the HTTP server and the in-process Babs tool-loop identical
    behaviour — same handler, same heavy-op cap, same post-write live refresh."""

    def test_unknown_endpoint_404(self):
        payload, status = sc._agent_invoke(None, "no-such-endpoint", {})
        assert status == 404 and "unknown endpoint" in payload["error"]

    def test_read_passthrough_no_app(self):
        # A read endpoint runs with app=None, returns its payload + 200, and
        # attempts no live-refresh (app is None).
        payload, status = sc._agent_invoke(None, "get-settings", {})
        assert status == 200 and payload.get("ok") is True

    def test_non_dict_body_coerced(self):
        payload, status = sc._agent_invoke(None, "get-settings", "not-a-dict")
        assert status == 200

    def test_write_in_refresh_map_fires_refresh(self):
        # A write endpoint listed in _AGENT_REFRESH_MAP must dispatch a UI-thread
        # repaint via call_from_thread. The stub's `_current_record is None`, so
        # the refresh itself no-ops cleanly — we only assert it was dispatched.
        calls = []

        class _Stub:
            _current_record = None

            def call_from_thread(self, fn, *a, **k):
                calls.append(fn)
                return fn(*a, **k)

        payload, status = sc._agent_invoke(
            _Stub(), "set-feature-color",
            {"feature_type": "CDS", "color": "#123456"})
        assert status == 200
        assert calls, "live-refresh was not dispatched for a mapped write"

    def test_read_does_not_fire_refresh(self):
        calls = []

        class _Stub:
            _current_record = None

            def call_from_thread(self, fn, *a, **k):
                calls.append(fn)
                return fn(*a, **k)

        sc._agent_invoke(_Stub(), "get-settings", {})
        assert not calls, "a read endpoint must not trigger a live-refresh"


class TestFindSequenceEndpoint:
    """`find-sequence` (2026-06-30): locate a DNA subseq / IUPAC motif on the
    loaded record — situational awareness for an autonomous agent."""

    def _app(self, seq, topology="linear"):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq(seq), id="x", name="x",
                        annotations={"topology": topology})
        return types.SimpleNamespace(_current_record=rec)

    def test_no_record_422(self):
        _, st = sc._h_find_sequence(
            types.SimpleNamespace(_current_record=None), {"query": "GAATTC"})
        assert st == 422

    def test_invalid_query_400(self):
        _, st = sc._h_find_sequence(self._app("ACGT" * 10), {"query": "XYZ!"})
        assert st == 400

    def test_forward_hits(self):
        seq = "TT" + "GAATTC" + "AAAA" + "GAATTC" + "TT"
        res = sc._h_find_sequence(self._app(seq),
                                  {"query": "GAATTC", "both_strands": False})
        assert res["n_hits"] == 2
        assert all(h["strand"] == 1 for h in res["hits"])
        assert res["hits"][0]["start"] == 2 and res["hits"][0]["end"] == 8

    def test_reverse_strand_hit(self):
        # query ACGGT (RC = ACCGT); only ACCGT is present → one reverse hit
        seq = "TTTT" + "ACCGT" + "TTTT"
        res = sc._h_find_sequence(self._app(seq), {"query": "ACGGT"})
        assert res["n_hits"] == 1
        assert res["hits"][0]["strand"] == -1 and res["hits"][0]["start"] == 4

    def test_palindrome_reported_once(self):
        seq = "TT" + "GAATTC" + "TT"           # GAATTC is palindromic
        res = sc._h_find_sequence(self._app(seq),
                                  {"query": "GAATTC", "both_strands": True})
        assert res["n_hits"] == 1 and res["hits"][0]["strand"] == 1

    def test_iupac_motif(self):
        seq = "GGACC" + "AA" + "GGTCC"         # GGWCC matches both
        res = sc._h_find_sequence(self._app(seq),
                                  {"query": "GGWCC", "both_strands": False})
        assert res["n_hits"] == 2

    def test_circular_wrap(self):
        seq = "TTC" + "AAAA" + "GAA"           # GAATTC spans the origin
        res = sc._h_find_sequence(self._app(seq, "circular"),
                                  {"query": "GAATTC", "both_strands": False})
        wrap = [h for h in res["hits"] if h["wraps"]]
        assert len(wrap) == 1 and wrap[0]["start"] == 7 and wrap[0]["end"] == 3


class TestExportMigrateArchiveEndpoint:
    """`export-migrate-archive` (2026-06-30): zip the whole workspace. Read-only
    to the data dir; import is deliberately NOT exposed."""

    def test_writes_zip(self, tmp_path):
        out = tmp_path / "ws.zip"
        res = sc._h_export_migrate_archive(None, {"path": str(out)})
        if isinstance(res, tuple):
            payload, status = res
            assert status == 200, payload
        else:
            payload = res
        assert payload.get("ok") is True
        assert out.exists() and out.stat().st_size > 0

    def test_rejects_non_zip(self, tmp_path):
        _, st = sc._h_export_migrate_archive(None, {"path": str(tmp_path / "ws.txt")})
        assert st == 400

    def test_missing_path_400(self):
        _, st = sc._h_export_migrate_archive(None, {})
        assert st == 400


class TestStatusViewState:
    """`status` now carries a `view` block (selection / cursor / map view) so an
    agent can act on what the user has highlighted (2026-06-30)."""

    def test_status_includes_view(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ACGT" * 10), id="x", name="x",
                        annotations={"topology": "circular"})
        app = types.SimpleNamespace(_current_record=rec, _source_path=None,
                                    _unsaved=False)
        res = sc._h_status(app, {})
        assert "view" in res and isinstance(res["view"], dict)
        # no panels mounted on this shim → all view fields null, but the shape
        # is stable so an agent can rely on it.
        assert set(res["view"]) == {"selection", "cursor_bp", "origin_bp",
                                    "map_mode", "selected_feature"}

    def test_status_view_null_when_no_record(self):
        app = types.SimpleNamespace(_current_record=None, _source_path=None,
                                    _unsaved=False)
        res = sc._h_status(app, {})
        assert res["view"] is None


# A 600 bp template (all ACGT, GFP block + a distinct mCherry block) used by
# the primer tests below. Long enough for `design-primers` detection mode
# (≥ 450 bp min product) and non-repeating so check-primer slices stay unique.
_PRIMER_TMPL = (
    "ATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGT"
    "GATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGA"
    "AAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTT"
    "GTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAG"
    "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAG"
    "GTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGC"
    "CGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCC"
    "TTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCAC"
    "CCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGC"
    "GTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGAC"
)


class TestDesignPrimersModes:
    """`design-primers` — detection (default), cloning, and the NEW
    generic (binding-only) mode. The endpoint had no test before this."""

    def test_missing_or_bad_template_400(self):
        assert sc._h_design_primers(None, {})[1] == 400
        assert sc._h_design_primers(None, {"template": "ZZZ", "start": 0,
                                           "end": 3})[1] == 400   # no IUPAC

    def test_bad_mode_400(self):
        r = sc._h_design_primers(None, {"template": _PRIMER_TMPL, "start": 0,
                                        "end": len(_PRIMER_TMPL),
                                        "mode": "wat"})
        assert r[1] == 400

    def test_detection_default(self):
        r = sc._h_design_primers(None, {"template": _PRIMER_TMPL,
                                        "start": 0, "end": len(_PRIMER_TMPL)})
        assert isinstance(r, dict) and r["ok"] and r["mode"] == "detection"
        # detection picks a diagnostic amplicon → fwd_seq/rev_seq + product_size
        assert r["result"].get("fwd_seq") and r["result"].get("rev_seq")
        assert r["result"].get("product_size")

    def test_generic_binding_only(self):
        n = len(_PRIMER_TMPL)
        r = sc._h_design_primers(None, {"template": _PRIMER_TMPL, "start": 0,
                                        "end": n, "mode": "generic"})
        assert isinstance(r, dict) and r["ok"] and r["mode"] == "generic"
        res = r["result"]
        # Binding-only: the forward primer IS a prefix of the template (no
        # 5' tail / RE site / overhang), and positions are reported.
        assert res["fwd_seq"] and _PRIMER_TMPL.startswith(res["fwd_seq"])
        assert "fwd_pos" in res and "rev_pos" in res
        assert res["fwd_tm"] and res["rev_tm"]
        # No tail keys leak in from the other modes.
        assert "fwd_full" not in res

    def test_generic_region_too_short_422(self):
        r = sc._h_design_primers(None, {"template": _PRIMER_TMPL, "start": 0,
                                        "end": 10, "mode": "generic"})
        assert r[1] == 422        # < 18 bp region

    def test_cloning_runs(self):
        # Pre-existing mode — just confirm it dispatches (ok, or a clean 422).
        r = sc._h_design_primers(None, {"template": _PRIMER_TMPL, "start": 20,
                                        "end": 220, "mode": "cloning",
                                        "site_5": "GAATTC", "site_3": "GGATCC"})
        assert isinstance(r, dict) or r[1] == 422


class TestCheckPrimer:
    """`check-primer` — single primer vs a template: Tm, GC%, and every
    3'-anchored binding site (both strands, wrap-aware)."""

    def test_missing_args_400(self):
        assert sc._h_check_primer(None, {})[1] == 400
        assert sc._h_check_primer(None, {"primer": "ACGTACGT"})[1] == 400
        assert sc._h_check_primer(None, {"primer": "ZZZ",
                                         "template": _PRIMER_TMPL})[1] == 400
        assert sc._h_check_primer(None, {"primer": "ACGTACGT",
                                         "template": "ZZZ"})[1] == 400

    def test_forward_site_exact(self):
        fwd = _PRIMER_TMPL[40:62]
        r = sc._h_check_primer(None, {"primer": fwd, "template": _PRIMER_TMPL,
                                      "circular": False})
        assert r["ok"] and r["binds"] and r["length"] == 22
        fwd_sites = [s for s in r["sites"] if s["orientation"] == "forward"]
        assert any(s["foot_start"] == 40 and s["ident_pct"] == 100.0
                   for s in fwd_sites)
        assert r["best_identity"] == 100.0
        assert r["tm"] is not None and 0 <= r["gc_pct"] <= 100
        assert fwd_sites[0]["confidence"] == "✓"   # 100% → check

    def test_reverse_site(self):
        rc = sc._rc(_PRIMER_TMPL[40:62])
        r = sc._h_check_primer(None, {"primer": rc, "template": _PRIMER_TMPL,
                                      "circular": False})
        rev = [s for s in r["sites"] if s["orientation"] == "reverse"]
        assert any(s["foot_start"] == 40 and s["ident_pct"] == 100.0
                   for s in rev)

    def test_no_binding(self):
        r = sc._h_check_primer(None, {"primer": "GGGGGGGGGGGGGGGGGGGG",
                                      "template": "ATATATATATATATATATATAT",
                                      "circular": False})
        assert (r["ok"] and not r["binds"] and r["n_sites"] == 0
                and r["best_identity"] is None)

    def test_internal_mismatch_partial_identity(self):
        # A single mismatch in the 5' half (3' seed intact → site still found)
        # gives identity over the FULL primer, not 100% — exactly what catches
        # a tailed / mutagenic primer. Deterministic: 1 of 22 bases differs.
        core = _PRIMER_TMPL[40:62]
        flip = "T" if core[2] != "T" else "G"
        primer = core[:2] + flip + core[3:]
        r = sc._h_check_primer(None, {"primer": primer,
                                      "template": _PRIMER_TMPL,
                                      "circular": False})
        fwd = [s for s in r["sites"]
               if s["orientation"] == "forward" and s["foot_start"] == 40]
        assert fwd and fwd[0]["mismatches"] == 1
        assert 90.0 < fwd[0]["ident_pct"] < 100.0

    def test_circular_wrap(self):
        # A primer spanning the origin binds only on a circular template.
        n = len(_PRIMER_TMPL)
        wrapper = _PRIMER_TMPL[-10:] + _PRIMER_TMPL[:12]
        circ = sc._h_check_primer(None, {"primer": wrapper,
                                         "template": _PRIMER_TMPL,
                                         "circular": True})
        lin = sc._h_check_primer(None, {"primer": wrapper,
                                        "template": _PRIMER_TMPL,
                                        "circular": False})
        assert any(s["foot_start"] == n - 10 and s["ident_pct"] == 100.0
                   for s in circ["sites"] if s["orientation"] == "forward")
        assert not any(s["foot_start"] == n - 10
                       for s in lin["sites"] if s["orientation"] == "forward")


class TestBatchMoveCopyPlasmids:
    """`copy-plasmids` / `move-plasmids` — bulk move/copy in ONE locked save
    (agent-API feedback P1-6). Per-item results; [INV-83] mirror re-stage on
    move with the active collection on either side."""

    def _ent(self, name, i):
        return {"id": f"id-{name}-{i}", "name": name, "gb_text": "",
                "size": 80}

    def _names(self, coll):
        c = next(c for c in sc._load_collections() if c["name"] == coll)
        return [p["name"] for p in c["plasmids"]]

    def _lib(self):
        return [p["name"] for p in sc._load_library()]

    def _seed(self):
        # CollC is the (neutral) active collection so non-active ops leave the
        # live mirror untouched; CollA/CollB hold the movable plasmids.
        sc._save_collections([
            {"name": "CollA", "plasmids": [self._ent("pA1", 1),
                                           self._ent("pA2", 2),
                                           self._ent("pA3", 3)],
             "saved": "2026-01-01"},
            {"name": "CollB", "plasmids": [self._ent("pB1", 1)],
             "saved": "2026-01-01"},
            {"name": "CollC", "plasmids": [], "saved": "2026-01-01"},
        ])
        sc._set_active_collection_name("CollC")
        sc._restore_library_from_active_collection()

    # ── validation ──────────────────────────────────────────────────────
    def test_bad_inputs(self):
        self._seed()
        assert sc._h_copy_plasmids(None, {})[1] == 400
        assert sc._h_copy_plasmids(None, {"names": [], "to": "CollB"})[1] == 400
        assert sc._h_copy_plasmids(None, {"names": ["pA1"]})[1] == 400
        assert sc._h_copy_plasmids(None,
                                   {"names": [1, 2], "to": "CollB"})[1] == 400
        assert sc._h_copy_plasmids(None,
                                   {"names": ["pA1"], "to": "Nope"})[1] == 404
        assert sc._h_move_plasmids(
            None, {"names": ["pA1"], "to": "CollA", "from": "Nope"})[1] == 404

    def test_copy_refuses_active_target(self):
        self._seed()
        assert sc._h_copy_plasmids(None,
                                   {"names": ["pA1"], "to": "CollC"})[1] == 409

    def test_copy_batch_per_item(self):
        self._seed()
        r = sc._h_copy_plasmids(
            None, {"names": ["pA1", "pA2", "ghost", "pB1"], "to": "CollB"})
        assert r["ok"] and r["copied"] == 2
        assert set(r["copied_names"]) == {"pA1", "pA2"}
        assert r["not_found"] == ["ghost"]
        assert r["conflict"] == ["pB1"]            # already in CollB
        assert set(self._names("CollA")) == {"pA1", "pA2", "pA3"}  # source kept
        assert set(self._names("CollB")) == {"pB1", "pA1", "pA2"}

    def test_copy_dedups_repeats(self):
        self._seed()
        r = sc._h_copy_plasmids(None,
                                {"names": ["pA1", "pA1", "pA1"], "to": "CollB"})
        assert r["copied"] == 1

    def test_copy_ambiguous_needs_from(self):
        sc._save_collections([
            {"name": "X", "plasmids": [self._ent("dup", 1)],
             "saved": "2026-01-01"},
            {"name": "Y", "plasmids": [self._ent("dup", 2)],
             "saved": "2026-01-01"},
            {"name": "W", "plasmids": [], "saved": "2026-01-01"},
            {"name": "Act", "plasmids": [], "saved": "2026-01-01"},
        ])
        sc._set_active_collection_name("Act")
        sc._restore_library_from_active_collection()
        r = sc._h_copy_plasmids(None, {"names": ["dup"], "to": "W"})
        assert r["ok"] and r["copied"] == 0 and r["ambiguous"] == ["dup"]
        r2 = sc._h_copy_plasmids(None,
                                 {"names": ["dup"], "to": "W", "from": "X"})
        assert r2["copied"] == 1

    # ── move ────────────────────────────────────────────────────────────
    def test_move_batch_per_item(self):
        self._seed()
        r = sc._h_move_plasmids(
            None, {"names": ["pA1", "pA3", "ghost"], "to": "CollB"})
        assert r["ok"] and r["moved"] == 2
        assert set(r["moved_names"]) == {"pA1", "pA3"}
        assert r["not_found"] == ["ghost"]
        assert self._names("CollA") == ["pA2"]        # pA1, pA3 removed
        assert set(self._names("CollB")) == {"pB1", "pA1", "pA3"}

    def test_move_into_active_remirrors(self):
        self._seed()
        r = sc._h_move_plasmids(None, {"names": ["pA1", "pA2"], "to": "CollC"})
        assert r["moved"] == 2
        assert set(self._names("CollC")) == {"pA1", "pA2"}
        assert set(self._lib()) == {"pA1", "pA2"}     # active mirror gained

    def test_move_out_of_active_remirrors(self):
        self._seed()
        sc._set_active_collection_name("CollA")
        sc._restore_library_from_active_collection()
        assert set(self._lib()) == {"pA1", "pA2", "pA3"}
        r = sc._h_move_plasmids(None, {"names": ["pA1", "pA2"], "to": "CollB"})
        assert r["moved"] == 2
        assert self._names("CollA") == ["pA3"]
        assert self._lib() == ["pA3"]                 # mirror shrank, L3-safe

    def test_move_name_collision_conflict_no_data_loss(self):
        sc._save_collections([
            {"name": "S", "plasmids": [
                {"id": "s1", "name": "shared", "gb_text": "", "size": 80}],
             "saved": "2026-01-01"},
            {"name": "T", "plasmids": [
                {"id": "t1", "name": "shared", "gb_text": "", "size": 80}],
             "saved": "2026-01-01"},
            {"name": "Act", "plasmids": [], "saved": "2026-01-01"},
        ])
        sc._set_active_collection_name("Act")
        sc._restore_library_from_active_collection()
        r = sc._h_move_plasmids(None,
                                {"names": ["shared"], "to": "T", "from": "S"})
        assert r["moved"] == 0 and r["conflict"] == ["shared"]
        # The source entry must NOT have been removed — no data loss.
        assert self._names("S") == ["shared"]
        assert self._names("T") == ["shared"]


class TestOnlineSearchGate:
    """`blast-online` / `hmmer-web` — the two-key egress gate + validation.
    The remote engines are monkeypatched so no network is touched."""

    def test_disarmed_by_default_403(self):
        sc._set_setting("allow_online_search", False)
        assert sc._h_blast_online(
            None, {"query": "ATGCATGCATGCATGCATGC"})[1] == 403
        assert sc._h_hmmer_web(None, {"query": "MKVLAAGG"})[1] == 403

    def test_setting_is_not_agent_settable(self):
        # The HUMAN half of the two-key gate: the setting exists but an agent
        # cannot flip it via set-setting (absent from the allowlist).
        assert "allow_online_search" not in sc._AGENT_SETTINGS_ALLOWLIST
        r = sc._h_set_setting(
            None, {"key": "allow_online_search", "value": True})
        assert r[1] == 400

    def test_armed_blast_calls_engine(self, monkeypatch):
        sc._set_setting("allow_online_search", True)
        calls = {}

        def fake(query, program, database, max_hits):
            calls.update(query=query, program=program,
                         database=database, max_hits=max_hits)
            return [{"accession": "X1", "description": "hit",
                     "identity_pct": 99.0}]
        monkeypatch.setattr("splicecraft_agent._ncbi_blast_online", fake)
        r = sc._h_blast_online(None,
                               {"query": "ATGCATGCATGCATGCATGCATGC"})
        assert r["ok"] and r["n_hits"] == 1 and r["program"] == "blastn"
        assert r["database"] == "nt" and calls["max_hits"] == 25

    def test_armed_blast_validation(self, monkeypatch):
        sc._set_setting("allow_online_search", True)
        monkeypatch.setattr("splicecraft_agent._ncbi_blast_online",
                            lambda *a: [])
        assert sc._h_blast_online(
            None, {"query": "ATGC", "program": "wat"})[1] == 400
        assert sc._h_blast_online(None, {})[1] == 400        # missing query
        assert sc._h_blast_online(
            None, {"query": "ATGCATGC", "database": "../etc"})[1] == 400

    def test_armed_blast_upstream_error_502(self, monkeypatch):
        sc._set_setting("allow_online_search", True)

        def boom(*a):
            raise RuntimeError("NCBI down")
        monkeypatch.setattr("splicecraft_agent._ncbi_blast_online", boom)
        assert sc._h_blast_online(
            None, {"query": "ATGCATGCATGC"})[1] == 502

    def test_armed_hmmer_calls_engine(self, monkeypatch):
        sc._set_setting("allow_online_search", True)
        seen = {}

        def fake(protein, max_hits):
            seen.update(protein=protein, max_hits=max_hits)
            return [{"acc": "PF00001", "name": "Pkinase"}]
        monkeypatch.setattr("splicecraft_agent._hmmer_web_hmmscan", fake)
        r = sc._h_hmmer_web(None,
                            {"query": "MKVLAAGGAGGSPVT", "max_hits": 5})
        assert r["ok"] and r["n_hits"] == 1 and seen["max_hits"] == 5


# ── blue field report (2026-08-24) ────────────────────────────────────────
# An agent driving the headless API through a plant-binary-vector build hit
# four failures that produced WRONG ANSWERS rather than errors. Each class
# below pins one of them.


def _circ(seq):
    """A circular SeqRecord — Type IIS cuts near an end need the wrap."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq(seq), id="t", name="t")
    rec.annotations["topology"] = "circular"
    rec.annotations["molecule_type"] = "DNA"
    return rec


class TestEnzymeNameResolution:
    """Report #1: `list-restriction-sites` knew BsmBI but not Esp3I, so an
    Esp3I safety check passed VACUOUSLY — `count: 0` reads as "no sites"
    when it actually meant "I never scanned for that name"."""

    # CGTCTC at 4 (BsmBI/Esp3I), GGTCTC at 50 (BsaI/BspTNI).
    SEQ = "AAAA" + "CGTCTC" + "A" * 40 + "GGTCTC" + "T" * 40

    def _app(self):
        return MockApp(record=_circ(self.SEQ))

    def test_isoschizomer_name_finds_its_own_sites(self):
        """The regression itself. Esp3I is in the catalog with BsmBI's exact
        (site, cut, cut) tuple, but the scanner collapses isoschizomers onto
        ONE label for the map — so a post-filter for the losing name came
        back empty."""
        r = sc._h_list_restriction_sites(self._app(), {"enzyme": "Esp3I"})
        assert r["count"] == 1, r
        assert r["sites"][0]["enzyme"] == "Esp3I"
        assert r["sites"][0]["start"] == 4
        # Same physical site the winning name reports.
        bsmbi = sc._h_list_restriction_sites(self._app(), {"enzyme": "BsmBI"})
        assert bsmbi["sites"][0]["start"] == r["sites"][0]["start"]
        assert bsmbi["sites"][0]["cut_bp"] == r["sites"][0]["cut_bp"]

    def test_every_requested_isoschizomer_gets_its_own_row(self):
        """Asking for both names returns both — a script filtering on either
        spelling finds its hits."""
        r = sc._h_list_restriction_sites(
            self._app(), {"enzymes": ["BsmBI", "Esp3I"]})
        assert {s["enzyme"] for s in r["sites"]} == {"BsmBI", "Esp3I"}
        assert {s["start"] for s in r["sites"]} == {4}

    def test_losing_side_of_the_collapse_for_bsai_too(self):
        """Not a CGTCTC special case: BspTNI loses to BsaI the same way."""
        r = sc._h_list_restriction_sites(self._app(), {"enzyme": "BspTNI"})
        assert [s["start"] for s in r["sites"]] == [50], r

    def test_unknown_enzyme_is_400_not_an_empty_result(self):
        """The core of the report: a check that cannot fail is worse than no
        check. `Esp3l` (lowercase L) used to scan nothing and answer clean."""
        out = sc._h_list_restriction_sites(self._app(), {"enzyme": "Esp3l"})
        assert isinstance(out, tuple) and out[1] == 400, out
        assert out[0]["unknown_enzymes"] == ["Esp3l"]
        assert "Esp3I" in out[0]["error"]          # near-miss suggestion

    def test_unknown_name_with_no_near_miss_still_errors(self):
        out = sc._h_list_restriction_sites(
            self._app(), {"enzymes": ["BsaI", "TotallyMadeUpI"]})
        assert isinstance(out, tuple) and out[1] == 400
        assert out[0]["unknown_enzymes"] == ["TotallyMadeUpI"]

    def test_case_insensitive(self):
        """REBASE capitalisation is not something anyone types reliably."""
        for spelling in ("bsai", "BSAI", "BsaI"):
            r = sc._h_list_restriction_sites(self._app(), {"enzyme": spelling})
            assert [s["enzyme"] for s in r["sites"]] == ["BsaI"], spelling

    def test_empty_list_refused(self):
        """`{"enzymes": []}` used to fall through to a full-catalog scan —
        the result looked scoped and wasn't."""
        out = sc._h_list_restriction_sites(self._app(), {"enzymes": []})
        assert isinstance(out, tuple) and out[1] == 400

    def test_equivalent_enzymes_names_what_the_answer_covers(self):
        """A zero-site answer is only trustworthy if you know which spellings
        it accounted for."""
        r = sc._h_list_restriction_sites(self._app(), {"enzyme": "BsmBI"})
        assert set(r["equivalent_enzymes"]) == {"Esp3I", "BsmBI-v2"}
        assert r["enzymes_scanned"] == ["BsmBI"]

    def test_cut_coordinates_and_site_ride_along(self):
        """Report #6: fragment-size arithmetic shouldn't need a private table
        of cut offsets."""
        r = sc._h_list_restriction_sites(self._app(), {"enzyme": "BsaI"})
        hit = r["sites"][0]
        assert hit["site"] == "GGTCTC"
        assert hit["cut_bp"] == 50 + 7            # GGTCTC(1/5): top cut
        assert hit["bottom_cut_bp"] == 50 + 11    # bottom cut

    def test_wrap_continuation_is_not_counted_as_a_site(self):
        """Sacred invariant #6 emits an origin-spanning site as a labeled
        tail piece PLUS an unlabeled head piece. The unlabeled one used to
        surface as a phantom `{"enzyme": ""}` row whenever no enzyme filter
        was passed, so `count` disagreed with itself between a filtered and
        an unfiltered call."""
        # EcoRI (GAATTC) straddling the origin: last 3 bases + first 3.
        seq = "TTC" + "A" * 60 + "GAA"
        r = sc._h_list_restriction_sites(_app_for(seq), {})
        assert all(s["enzyme"] for s in r["sites"]), r
        assert r["count"] == len(r["sites"])
        scoped = sc._h_list_restriction_sites(_app_for(seq), {"enzyme": "EcoRI"})
        assert scoped["count"] == 1
        assert [s["enzyme"] for s in r["sites"] if s["enzyme"] == "EcoRI"] \
            == ["EcoRI"]


def _app_for(seq):
    return MockApp(record=_circ(seq))


class TestCommercialSynonymsOverTheApi:
    """2026-08-26 report #1, one layer out from [INV-187]: the isoschizomers
    already IN the catalog resolved, but the names people actually buy —
    Thermo's Eco31I/LguI/AarI — did not, so `{"enzymes": ["Eco31I"]}` came
    back `400` on an enzyme SpliceCraft knows perfectly well."""

    # GGTCTC at 4 (BsaI/Eco31I), GCTCTTC at 50 (SapI/BspQI/LguI).
    SEQ = "AAAA" + "GGTCTC" + "A" * 40 + "GCTCTTC" + "T" * 40

    def _app(self):
        return MockApp(record=_circ(self.SEQ))

    def test_a_synonym_finds_the_sites(self):
        r = sc._h_list_restriction_sites(self._app(), {"enzyme": "Eco31I"})
        assert [s["start"] for s in r["sites"]] == [4], r

    def test_rows_are_labelled_with_the_name_you_asked_for(self):
        """Reporting `BsaI` to someone who asked for `Eco31I` would leave the
        obvious `s["enzyme"] == "Eco31I"` filter empty — the same silent zero
        one layer down."""
        r = sc._h_list_restriction_sites(self._app(), {"enzyme": "Eco31I"})
        assert [s["enzyme"] for s in r["sites"]] == ["Eco31I"]
        assert r["enzymes_scanned"] == ["Eco31I"]
        # …and the answer still says which catalog enzyme that IS.
        assert "BsaI" in r["equivalent_enzymes"]

    def test_synonym_and_catalog_name_together_both_report(self):
        r = sc._h_list_restriction_sites(
            self._app(), {"enzymes": ["Eco31I", "BsaI"]})
        assert {s["enzyme"] for s in r["sites"]} == {"Eco31I", "BsaI"}
        assert {s["start"] for s in r["sites"]} == {4}

    def test_lgui_reaches_the_sapi_site(self):
        r = sc._h_list_restriction_sites(self._app(), {"enzyme": "LguI"})
        assert [s["start"] for s in r["sites"]] == [50], r
        assert r["sites"][0]["site"] == "GCTCTTC"

    def test_digest_says_what_a_synonym_became(self):
        """`cuts[].enzyme` carries the catalog name, so the mapping has to be
        visible or the caller reads "my enzyme didn't cut"."""
        r = sc._h_digest(None, {"sequence": self.SEQ, "enzymes": ["Eco31I"],
                                "circular": False})
        assert not isinstance(r, tuple), r
        assert r["resolved_enzymes"] == {"Eco31I": "BsaI"}
        assert [c["enzyme"] for c in r["cuts"]] == ["BsaI"]

    def test_digest_omits_the_mapping_when_nothing_was_renamed(self):
        r = sc._h_digest(None, {"sequence": self.SEQ, "enzymes": ["BsaI"],
                                "circular": False})
        assert "resolved_enzymes" not in r


class TestNeoschizomersOverTheApi:
    """Same report: XmaI and SmaI read the same six bases and leave different
    ends, and the endpoint was inheriting the MAP's one-bar-per-span collapse,
    so naming both returned one."""

    SEQ = "AAAA" + "CCCGGG" + "T" * 50 + "GGTACC" + "A" * 50

    def _app(self):
        return MockApp(record=_circ(self.SEQ))

    def test_both_neoschizomers_report_with_their_own_cut(self):
        r = sc._h_list_restriction_sites(
            self._app(), {"enzymes": ["XmaI", "SmaI"]})
        cuts = {s["enzyme"]: s["cut_bp"] for s in r["sites"]}
        assert set(cuts) == {"XmaI", "SmaI"}
        assert cuts["XmaI"] == 4 + 1          # C^CCGGG
        assert cuts["SmaI"] == 4 + 3          # CCC^GGG (blunt)

    def test_acc65i_and_kpni_both_report(self):
        r = sc._h_list_restriction_sites(
            self._app(), {"enzymes": ["Acc65I", "KpnI"]})
        cuts = {s["enzyme"]: s["cut_bp"] for s in r["sites"]}
        assert set(cuts) == {"Acc65I", "KpnI"}
        assert cuts["Acc65I"] == 60 + 1       # G^GTACC
        assert cuts["KpnI"] == 60 + 5         # GGTAC^C

    def test_a_neoschizomer_is_not_offered_as_equivalent(self):
        """`equivalent_enzymes` means "this answer already covers it". SmaI
        does NOT cover XmaI — different overhang."""
        r = sc._h_list_restriction_sites(self._app(), {"enzyme": "XmaI"})
        assert "SmaI" not in (r.get("equivalent_enzymes") or [])


class TestListEnzymes:
    """Report #6: cut offsets weren't reachable, so a script verifying digest
    fragment sizes had to hardcode `SalI G^TCGAC = 1` and `KpnI GGTAC^C = 5`
    from memory."""

    def _get(self, name):
        r = sc._h_list_enzymes(None, {"names": [name]})
        assert not isinstance(r, tuple), r
        return r["enzymes"][0]

    def test_notation_matches_rebase(self):
        """The two the report cited by hand, plus a Type IIS and a blunt."""
        assert self._get("SalI")["notation"] == "G^TCGAC"
        assert self._get("KpnI")["notation"] == "GGTAC^C"
        assert self._get("BsaI")["notation"] == "GGTCTC(1/5)"
        assert self._get("EcoRV")["notation"] == "GAT^ATC"

    def test_overhang_kind_agrees_with_the_digest_engine(self):
        """`overhang_kind` is derived with the same rule `_enzyme_cuts_impl`
        uses, so a size check can't disagree with the fragments it checks."""
        assert (self._get("SalI")["overhang_kind"],
                self._get("SalI")["overhang_len"]) == ("5'", 4)
        assert (self._get("KpnI")["overhang_kind"],
                self._get("KpnI")["overhang_len"]) == ("3'", 4)
        assert (self._get("EcoRV")["overhang_kind"],
                self._get("EcoRV")["overhang_len"]) == ("blunt", 0)

    def test_cut_offsets_reproduce_a_real_digest(self):
        """Independent check: the published offsets must place the cut where
        `_enzyme_cuts` actually puts it."""
        seq = "A" * 20 + "GTCGAC" + "T" * 20
        row = self._get("SalI")
        cuts = sc._enzyme_cuts(seq, ["SalI"], circular=False)
        assert cuts[0]["top"] == 20 + row["fwd_cut"]
        assert cuts[0]["bot"] == 20 + row["rev_cut"]

    def test_type_iis_flag_and_aliases(self):
        assert self._get("BsaI")["type_iis"] is True
        assert self._get("EcoRI")["type_iis"] is False
        # Catalog isoschizomer + the commercial synonym (`_ENZYME_ALIASES`) —
        # a caller who buys "Eco31I" has to be able to see it IS this row.
        assert self._get("BsaI")["aliases"] == ["BspTNI", "Eco31I"]

    def test_unknown_name_400_with_suggestion(self):
        out = sc._h_list_enzymes(None, {"names": ["Bsa1"]})
        assert isinstance(out, tuple) and out[1] == 400
        assert "BsaI" in out[0]["error"]

    def test_search_by_recognition_site(self):
        r = sc._h_list_enzymes(None, {"search": "GGTCTC"})
        assert {e["name"] for e in r["enzymes"]} == {"BsaI", "BspTNI"}

    def test_type_iis_filter(self):
        r = sc._h_list_enzymes(None, {"type_iis_only": True})
        assert r["count"] > 0
        assert all(e["type_iis"] for e in r["enzymes"])


class TestSpanContains:
    """Report #2: with no circular-aware interval helper, "is this feature
    inside the T-DNA" had to be hand-rolled, and the obvious linear form
    reported a false FAIL on a CORRECT construct whose borders wrap the
    origin."""

    L = 10_000
    TDNA = {"start": 9500, "end": 800}          # 1300 bp across the origin

    def _q(self, inner, **kw):
        body = {"outer": self.TDNA, "length": self.L, "inner": inner}
        body.update(kw)
        r = sc._h_span_contains(None, body)
        assert not isinstance(r, tuple), r
        return r

    def test_the_linear_form_this_replaces_is_wrong(self):
        """Not a tautology: assert the naive predicate actually disagrees, so
        this test fails if someone reimplements it linearly."""
        inner = {"start": 100, "end": 700}
        naive = (self.TDNA["start"] <= inner["start"] < self.TDNA["end"])
        assert naive is False                    # what the report hit
        assert self._q(inner)["results"][0]["contains"] is True

    def test_contains_on_both_sides_of_the_origin(self):
        rows = self._q([{"start": 9600, "end": 9900},
                        {"start": 100,  "end": 700},
                        {"start": 9800, "end": 200}])["results"]
        assert [r["contains"] for r in rows] == [True, True, True]
        assert rows[2]["wraps"] is True

    def test_excludes_what_is_outside_and_what_straddles_a_border(self):
        rows = self._q([{"start": 2000, "end": 3000},
                        {"start": 700,  "end": 900}])["results"]
        assert [r["contains"] for r in rows] == [False, False]

    def test_single_bp(self):
        assert self._q({"bp": 9999})["results"][0]["contains"] is True
        assert self._q({"bp": 5000})["results"][0]["contains"] is False

    def test_one_response_shape_for_one_or_many(self):
        """Report #4's complaint applied to a new endpoint: one span in still
        yields a `results` LIST, so callers never branch on the shape."""
        one = self._q({"bp": 9999})
        many = self._q([{"bp": 9999}, {"bp": 5000}])
        assert len(one["results"]) == 1 and len(many["results"]) == 2
        assert one["all_contained"] is True and many["all_contained"] is False
        assert one["outer"]["length"] == 1300 and one["outer"]["wraps"] is True

    def test_bare_start_end_pair_is_one_span_not_two_specs(self):
        r = self._q([100, 700])
        assert len(r["results"]) == 1
        assert (r["results"][0]["start"], r["results"][0]["end"]) == (100, 700)

    def test_length_falls_back_to_the_loaded_record(self):
        app = MockApp(record=_circ("A" * 10_000))
        r = sc._h_span_contains(
            app, {"outer": self.TDNA, "inner": {"bp": 9999}})
        assert r["length"] == 10_000 and r["results"][0]["contains"] is True

    def test_no_record_and_no_length_is_422(self):
        out = sc._h_span_contains(
            MockApp(record=None), {"outer": self.TDNA, "inner": {"bp": 1}})
        assert isinstance(out, tuple) and out[1] == 422

    def test_bad_input_is_refused(self):
        for body in ({"outer": self.TDNA, "length": self.L},
                     {"inner": {"bp": 1}, "length": self.L},
                     {"outer": self.TDNA, "inner": {"bp": 1}, "length": 0},
                     {"outer": self.TDNA, "inner": {"start": "x", "end": 5},
                      "length": self.L}):
            out = sc._h_span_contains(None, body)
            assert isinstance(out, tuple) and out[1] == 400, body


class TestSimulationOkReflectsTheOutcome:
    """Report #4: the dry-run endpoints stamped `ok: true` for a reaction
    that FAILED — `result.ok` said false underneath — so the natural
    `if r["ok"]:` accepted a design that cannot be built."""

    @staticmethod
    def _cassette(oh5, body, oh3):
        return "GGTCTC" + "A" + oh5 + body + oh3 + "A" + "GAGACC"

    def test_failed_golden_gate_reports_ok_false(self):
        # Overhangs that cannot chain into a circle.
        A = self._cassette("GGAG", "ATGAAACCCGGGTTTACGT" * 2, "AATG")
        V = self._cassette("TTTT", "GGGGCCCCAAAATTTT" * 8, "CCCC")
        r = sc._h_simulate_golden_gate(None, {"parts": [A], "vector": V})
        assert r["result"]["ok"] is False
        assert r["ok"] is False, "top-level ok must not outrank result.ok"

    def test_successful_golden_gate_still_reports_ok_true(self):
        A = self._cassette("GGAG", "ATGAAACCCGGGTTTACGT" * 2, "AATG")
        B = self._cassette("AATG", "TTGCATGCATGCTAGCTAG" * 2, "CGCT")
        V = self._cassette("CGCT", "GGGGCCCCAAAATTTT" * 8, "GGAG")
        r = sc._h_simulate_golden_gate(None, {"parts": [A, B], "vector": V})
        assert r["ok"] is True and r["result"]["ok"] is True

    def test_failed_gibson_reports_ok_false(self):
        r = sc._h_simulate_gibson(None, {"fragments": [
            {"name": "a", "sequence": "ACGTACGT" * 20},
            {"name": "b", "sequence": "TTTTGGGG" * 20}]})
        assert r["result"]["success"] is False
        assert r["ok"] is False

    def test_successful_gibson_still_reports_ok_true(self):
        ov = "CGATCGATCGATCGATCGAT"
        r = sc._h_simulate_gibson(None, {"circular": False, "fragments": [
            {"name": "a", "sequence": "A" * 60 + ov},
            {"name": "b", "sequence": ov + "T" * 60}]})
        assert r["ok"] is True and r["result"]["success"] is True


class TestGoldenGateAnnotationParity:
    """Report #7: `golden-gate-assemble` has no `carry_annotations`. It saved
    a feature-BARE product, and passing the flag was silently ignored — a
    bare construct believed to be annotated looks identical to an annotated
    one until you open it."""

    @staticmethod
    def _design():
        cas = lambda o5, b, o3: "GGTCTC" + "A" + o5 + b + o3 + "A" + "GAGACC"
        return (cas("GGAG", "ATGAAACCCGGGTTTACGT" * 2, "AATG"),
                cas("AATG", "TTGCATGCATGCTAGCTAG" * 2, "CGCT"),
                cas("CGCT", "GGGGCCCCAAAATTTT" * 8, "GGAG"))

    def test_carry_annotations_is_refused_not_ignored(self):
        A, B, V = self._design()
        out = sc._h_golden_gate_assemble(
            None, {"parts": [A, B], "vector": V, "product_name": "ggx",
                   "carry_annotations": True})
        assert isinstance(out, tuple) and out[1] == 400, out
        assert "carry_annotations" in out[0]["error"]
        # Nothing was saved.
        assert not any(e["name"] == "ggx" for e in sc._load_library())

    def test_enzyme_is_still_accepted(self):
        """The guard must not over-reject a param the endpoint DOES take."""
        A, B, V = self._design()
        r = sc._h_golden_gate_assemble(
            None, {"parts": [A, B], "vector": V, "product_name": "gge",
                   "enzyme": "BsaI"})
        assert r["saved_name"] == "gge"

    def test_product_is_feature_bare_and_the_docstring_says_so(self):
        """Pins the CURRENT truth so a future `carry_annotations` here has to
        update both the behaviour and the promise."""
        A, B, V = self._design()
        r = sc._h_golden_gate_assemble(
            None, {"parts": [A, B], "vector": V, "product_name": "ggbare"})
        ent = next(e for e in sc._load_library() if e["name"] == "ggbare")
        assert ent["n_feats"] == 0
        rec = sc._gb_text_to_record(ent["gb_text"])
        assert [f for f in rec.features if f.type != "source"] == []
        assert "carries NO features" in sc._h_golden_gate_assemble.__doc__


class TestTraditionalCloneMarkerDefault:
    """Report #5: a directional two-enzyme digest leaves BOTH vector fragments
    with one cut of each kind, so both accept the insert and the caller had to
    read `simulate-traditional-cloning` and hand `vector_frag_idx` back for a
    reaction with exactly one sensible product. The 409 now resolves itself
    when the vector's ORIGIN / RESISTANCE annotations say which half is the
    backbone — never by fragment size."""

    VECTOR = TestTraditionalCloning.VECTOR      # EcoRI · A×200 · BamHI · C×60
    INSERT = TestTraditionalCloning.INSERT

    def _base(self):
        return {"vector_seq": self.VECTOR,
                "vector_enzymes": ["EcoRI", "BamHI"],
                "insert_seq": self.INSERT,
                "insert_enzymes": ["EcoRI", "BamHI"]}

    @staticmethod
    def _seed(feats, name="MarkedVector", seq=None):
        """Seed a vector entry carrying `feats` = [(start, end, type, label)]."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        s = TestTraditionalCloneMarkerDefault.VECTOR if seq is None else seq
        rec = SeqRecord(Seq(s), id="MARKV", name="MARKV",
                        annotations={"molecule_type": "DNA",
                                       "topology": "circular"})
        for start, end, ftype, label in feats:
            rec.features.append(SeqFeature(
                FeatureLocation(start, end, strand=1), type=ftype,
                qualifiers={"label": [label]}))
        sc._save_library([{
            "id": "MARKV", "name": name, "size": len(s), "n_feats": len(feats),
            "added": "2026-08-24", "gb_text": sc._record_to_gb_text(rec)}])

    def test_marker_settles_the_ambiguity(self, isolated_library):
        """AmpR sits in the A×200 half (vector fragment 0). No
        `vector_frag_idx` passed — the clone should just work."""
        self._seed([(20, 190, "CDS", "AmpR")])
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="autoclone", vector_name="MarkedVector"))
        assert not isinstance(r, tuple), r
        assert r["vector_frag_idx"] == 0
        assert r["vector_frag_auto_picked"] is True
        assert "AmpR" in r["vector_frag_pick_reason"]
        assert any(e["name"] == "autoclone" for e in sc._load_library())

    def test_without_a_vector_name_it_still_refuses(self, isolated_library):
        """No annotations to judge from → the 409 stands, unchanged."""
        r = sc._h_traditional_clone(
            None, dict(self._base(), product_name="c"))
        assert isinstance(r, tuple) and r[1] == 409

    def test_markers_on_both_halves_is_a_real_ambiguity(self, isolated_library):
        """If both fragments look like backbone, nothing has been settled —
        refuse rather than pick one."""
        self._seed([(20, 190, "CDS", "AmpR"),
                    (215, 260, "rep_origin", "ColE1 ori")])
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="c", vector_name="MarkedVector"))
        assert isinstance(r, tuple) and r[1] == 409

    def test_it_does_not_pick_by_size(self, isolated_library):
        """The decisive test. Put the marker on the SMALLER fragment (the C×60
        half, vector fragment 1): a size heuristic would take the 206 bp half,
        the marker rule must take the 60 bp one."""
        self._seed([(215, 260, "rep_origin", "ColE1 ori")])
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="smallpick", vector_name="MarkedVector"))
        assert not isinstance(r, tuple), r
        assert r["vector_frag_idx"] == 1, "picked the big half — size crept in"
        assert r["vector_frag_auto_picked"] is True

    def test_an_explicit_index_still_wins(self, isolated_library):
        """The caller's word beats the inference, even against the marker."""
        self._seed([(20, 190, "CDS", "AmpR")])
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="explicit", vector_name="MarkedVector",
            vector_frag_idx=1))
        assert r["vector_frag_idx"] == 1
        assert "vector_frag_auto_picked" not in r

    def test_sequence_mismatch_falls_back_to_the_409(self, isolated_library):
        """Same exact-match gate `carry_annotations` uses: markers read off a
        rotated or edited entry would be attributed to the wrong fragment, so
        an inexact match declines to guess instead of guessing wrong."""
        self._seed([(20, 190, "CDS", "AmpR")], seq=self.VECTOR + "TTTT")
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="c", vector_name="MarkedVector"))
        assert isinstance(r, tuple) and r[1] == 409

    def test_unknown_vector_name_does_not_break_the_clone(self, isolated_library):
        """The inference is ADVISORY — a name that resolves to nothing must
        leave the endpoint exactly as it was, not error."""
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="c", vector_name="ghost",
            vector_frag_idx=0))
        assert not isinstance(r, tuple), r
        assert r["vector_frag_idx"] == 0

    def test_auto_pick_leaves_no_silent_choice(self, isolated_library):
        """An automatic pick you can't see is what the 409 was protecting
        against — the response must always say it happened and why."""
        self._seed([(20, 190, "CDS", "AmpR")])
        r = sc._h_traditional_clone(None, dict(
            self._base(), product_name="loud", vector_name="MarkedVector"))
        assert set(("vector_frag_auto_picked", "vector_frag_pick_reason")) \
            <= set(r)
        assert str(r["vector_frag_idx"]) in r["vector_frag_pick_reason"]


class TestEnzymeInputHardening:
    """Edge-case sweep on the name-resolution surface added for report #1.
    Every case here produced either a wrong answer or an unbounded cost."""

    SEQ = "GAATTC" + "A" * 500

    def _app(self):
        return MockApp(record=_circ(self.SEQ))

    def test_too_many_names_is_refused(self):
        """`_resolve_enzyme_names` is O(catalog) per unrecognised name, so an
        unbounded list is a CPU sink in a read endpoint. 500 matches the cap
        `digest` already had."""
        out = sc._h_list_restriction_sites(
            self._app(), {"enzymes": ["EcoRI"] * 501})
        assert isinstance(out, tuple) and out[1] == 400
        assert "too many" in out[0]["error"]
        ok = sc._h_list_restriction_sites(
            self._app(), {"enzymes": ["EcoRI"] * 500})
        assert not isinstance(ok, tuple)

    def test_a_huge_name_is_refused_on_shape_and_never_echoed(self):
        """A payload can carry a megabyte string; echoing it verbatim turns a
        validation error into a megabyte response. It is now refused on SHAPE
        (no enzyme name is over 64 characters) rather than lowercased, hashed
        and difflib'd on its way to "unknown enzyme" — but the property that
        matters is unchanged: none of it comes back."""
        huge = "X" * 100_000
        out = sc._h_list_restriction_sites(self._app(), {"enzyme": huge})
        assert isinstance(out, tuple) and out[1] == 400
        assert len(out[0]["error"]) < 500
        assert "64" in out[0]["error"]
        assert huge[:200] not in json.dumps(out[0])

    def test_a_64_char_name_still_reaches_the_unknown_path(self):
        """The shape gate must not swallow the near-miss machinery for names
        that are merely wrong rather than absurd."""
        out = sc._h_list_restriction_sites(self._app(), {"enzyme": "X" * 64})
        assert isinstance(out, tuple) and out[1] == 400
        assert out[0]["unknown_enzymes"] == ["X" * 64]

    def test_many_unknowns_report_a_true_total_not_a_truncated_one(self):
        out = sc._h_list_restriction_sites(
            self._app(), {"enzymes": [f"Fake{i}I" for i in range(200)]})
        assert isinstance(out, tuple) and out[1] == 400
        assert out[0]["n_unknown"] == 200
        assert len(out[0]["unknown_enzymes"]) <= 10   # quoted set is bounded
        assert len(out[0]["error"]) < 2000

    def test_suggestion_work_is_bounded(self):
        """Only the first few unknowns get near-misses computed; the rest are
        still reported."""
        names = [f"Nope{i}I" for i in range(100)]
        _resolved, unknown = sc._resolve_enzyme_names(names)
        assert len(unknown) == 100
        assert sum(1 for _w, near in unknown if near) <= \
            sc._RESOLVE_SUGGEST_LIMIT

    def test_min_length_longer_than_any_site_is_refused(self):
        """The vacuous-pass shape again: `min_length: 99` can only ever return
        zero sites, and `count: 0` reads as 'this construct is clean'."""
        out = sc._h_list_restriction_sites(self._app(), {"min_length": 99})
        assert isinstance(out, tuple) and out[1] == 400
        assert "BASE PAIRS" in out[0]["error"]

    def test_the_min_length_bound_comes_from_the_catalog(self):
        """Reading it from the catalog means a longer-site enzyme can't make
        the guard reject a legitimate filter."""
        longest = max(len(v[0]) for v in sc._all_enzymes().values())
        assert not isinstance(
            sc._h_list_restriction_sites(self._app(), {"min_length": longest}),
            tuple)
        assert sc._h_list_restriction_sites(
            self._app(), {"min_length": longest + 1})[1] == 400

    def test_a_negative_min_length_is_treated_as_no_filter(self):
        r = sc._h_list_restriction_sites(self._app(), {"min_length": -5})
        assert not isinstance(r, tuple) and r["count"] >= 1

    def test_list_enzymes_caps_names_too(self):
        out = sc._h_list_enzymes(None, {"names": ["BsaI"] * 501})
        assert isinstance(out, tuple) and out[1] == 400

    def test_digest_resolves_case_but_still_reports_rather_than_errors(self):
        """`digest`'s documented contract is REPORT-don't-error on an unknown
        name. Case-folding must not turn that into a refusal — but a lowercase
        name must actually cut, which it previously did not."""
        seq = "AAAA" + "GAATTC" + "AAAA"
        low = sc._h_digest(None, {"sequence": seq, "enzymes": ["ecori"],
                                  "circular": False})
        assert not isinstance(low, tuple), low
        assert low["n_cuts"] == 1 and low["unknown_enzymes"] == []
        mixed = sc._h_digest(None, {"sequence": seq,
                                    "enzymes": ["EcoRI", "NotAnEnzymeI"],
                                    "circular": False})
        assert not isinstance(mixed, tuple)
        assert mixed["n_cuts"] == 1
        assert mixed["unknown_enzymes"] == ["NotAnEnzymeI"]


class TestSpanContainsHardening:
    def test_an_empty_inner_list_is_refused(self):
        """`all([])` is True, so an empty list answered `all_contained: true`
        having checked nothing — the worst possible reading."""
        out = sc._h_span_contains(
            None, {"outer": {"start": 0, "end": 10}, "length": 100,
                   "inner": []})
        assert isinstance(out, tuple) and out[1] == 400

    def test_out_of_range_coordinates_are_wrapped_but_reported(self):
        """Circular coordinates make `-1` legitimate, so a far-out value is
        normalised rather than refused — but the answer looks equally
        confident either way, so it must say it happened."""
        r = sc._h_span_contains(
            None, {"outer": {"start": 0, "end": 10}, "length": 100,
                   "inner": {"bp": 10 ** 9}})
        assert not isinstance(r, tuple)
        assert r["warnings"] and "outside" in r["warnings"][0]

    def test_in_range_coordinates_carry_no_warning(self):
        r = sc._h_span_contains(
            None, {"outer": {"start": 0, "end": 10}, "length": 100,
                   "inner": {"bp": 5}})
        assert "warnings" not in r

    def test_a_negative_coordinate_still_works(self):
        """-1 means the last base; that has to keep working."""
        r = sc._h_span_contains(
            None, {"outer": {"start": 90, "end": 0}, "length": 100,
                   "inner": {"bp": -1}})
        assert r["results"][0]["start"] == 99
        assert r["results"][0]["contains"] is True
