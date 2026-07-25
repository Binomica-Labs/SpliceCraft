"""splicecraft_opentrons — the Opentrons OT-2 protocol compiler + LAN client (L1).

Pure, app-free core for driving an **Opentrons OT-2** liquid handler from
SpliceCraft. Two halves, both unit-testable in isolation:

1. **The compiler** (offline, no robot, no network). Turns a plain-dict
   *transfer plan* — a pipette, some labware on deck slots, and a list of
   ``from → to`` well transfers with volumes — into a ready-to-run **Opentrons
   Python Protocol API v2** ``.py`` text, after validating it (known pipette /
   labware, in-range wells, sane volumes, unique slots). The design deliberately
   mirrors SpliceCraft's "simulate the steps" ethos: the emitted protocol is a
   real, reviewable recipe the user (or the robot's own analysis) can vet before
   anything moves.

2. **The LAN client** (online, opt-in). A stdlib-``urllib`` wrapper over the
   OT-2 ``robot-server`` HTTP API (port 31950): health / instruments / lights,
   upload-and-analyse a protocol, and a **gated** physical run. The gate is the
   whole point — ``_ot2_run_protocol`` refuses to actuate unless the robot's own
   analysis returns ``ok`` *and* the caller passes ``confirm=True``. Physical
   motion has no ``.bak``; it is treated as strictly more dangerous than a
   data-dir write.

Why a sibling, not hub code: none of this is app-coupled — it's compute +
network. The Textual plate-designer UI and the ``ot2-*`` agent endpoints stay
hub-side / in ``splicecraft_agent``; they call *into* here. This module imports
only L0 siblings (``logging`` / ``state`` / ``util``) and never reaches the hub.

SECURITY / hardening
--------------------
* The OT-2 is a **LAN** device (a private RFC-1918 address / ``.local`` name),
  so its calls deliberately bypass the public-host SSRF assertion in
  ``splicecraft_net`` (which *refuses* private/loopback hosts) — exactly as the
  BABS engine does for a local Ollama server. They are still bounded on every
  axis: explicit timeouts, a hard response-size cap, bounded poll loops, and a
  no-redirect opener (the robot never legitimately 3xx-redirects an API call).
* No SpliceCraft data-dir writes happen here — the engine is compute + network
  only. The compiler returns text; the caller decides where (if anywhere) to
  save it, through the usual chokepoint.
* The physical-run gate (analysis ``ok`` + explicit ``confirm``) means no code
  path here can move the gantry by accident.
"""
from __future__ import annotations

import json
import math
import socket
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

import splicecraft_logging as _logging
import splicecraft_state as _state
import splicecraft_util as _util

_log = _logging._log


# ── Robot-server connection tunables ────────────────────────────────────────────
_OT2_PORT = 31950
# The OT-2 robot-server negotiates an API version via this header; "4" is inside
# the supported range of every OT-2 build this project has tested against and is
# what the upload/analyse/run flow was validated with end-to-end.
_OT2_API_HEADER = "4"

# Timeouts (seconds). Short calls must fail fast (a powered-down robot should not
# hang a worker); the upload + analysis + run budgets are generous because
# server-side protocol analysis and real liquid handling legitimately take time.
_OT2_SHORT_TIMEOUT = 10          # health / instruments / lights / run control
_OT2_UPLOAD_TIMEOUT = 90         # POST /protocols (also kicks off analysis)
_OT2_ANALYSIS_POLL_TIMEOUT = 180 # wait for server-side analysis to complete
_OT2_RUN_POLL_TIMEOUT = 1800     # a physical run: 30 min ceiling
_OT2_POLL_INTERVAL = 1.0

# ── Discovery tunables ──────────────────────────────────────────────────────────
# Find robots by sweeping the local subnet(s) for a live robot-server. A dead host
# hits the full probe timeout (a closed port refuses instantly), so keep it short
# and the pool wide: one /24 (254 hosts) then sweeps in a few seconds.
_OT2_DISCOVER_PROBE_TIMEOUT = 0.8   # per-host /health probe during a sweep
_OT2_DISCOVER_WORKERS = 64          # concurrent probes (bounded pool)
_OT2_DISCOVER_MAX_HOSTS = 1024      # hard cap on addresses probed in one sweep
# A /health reply is well under 4 KB; cap the probe tight so a rogue/large service
# answering on :31950 can't make every worker buffer the 24 MB default across a sweep.
_OT2_DISCOVER_HEALTH_MAX_BYTES = 64 * 1024
# Wall-clock backstop for a whole sweep — a well-behaved /24 finishes in a few
# seconds; this only trips on a host that trickles bytes to stall a worker.
_OT2_DISCOVER_OVERALL_TIMEOUT = 30.0
# One-shot mDNS multicast to catch robots the flat subnet sweep would miss.
_OT2_MDNS_GROUP = "224.0.0.251"
_OT2_MDNS_PORT = 5353
_OT2_MDNS_TIMEOUT = 1.5             # seconds to collect responders

# Response-size cap (bytes). A protocol analysis with hundreds of commands can be
# a few MB; 24 MB is generous vs. real payloads yet defends against a misbehaving
# endpoint streaming without end.
_OT2_JSON_MAX_BYTES = 24 * 1024 * 1024
# Upload cap: a compiled transfer protocol is small; refuse a pathologically
# large one before pushing it over the wire.
_OT2_MAX_PROTOCOL_BYTES = 4 * 1024 * 1024

# Plan sanity caps — defend the compiler + robot from garbage / runaway input.
_OT2_MAX_TRANSFERS = 10000        # far beyond any real plate-prep run
_OT2_MAX_VOLUME_UL = 100000.0     # 100 mL — reject absurd volumes (typos, Infinity)
_OT2_MAX_STEPS = 5000             # a multi-step protocol far beyond any real run
_OT2_MAX_DELAY_S = 86400.0        # 24 h — reject absurd delays
_OT2_MAX_ITEMS = 10000            # normalise: far beyond any real plate-prep run
_OT2_MAX_OFFSET_MM = 5.0          # labware offset is a fine correction — reject fat-fingers
_OT2_MAX_POSCHECK_WELLS = 384     # position-check: cap the per-labware move list (96-well ×4)

# The multi-step "protocol designer" model (a plan may carry `steps` instead of
# the legacy `transfers`). LIQUID steps need a pipette + tip rack + labware;
# control steps (delay/pause/comment) don't.
_OT2_STEP_TYPES = ("transfer", "distribute", "consolidate", "mix",
                   "delay", "pause", "comment")
_OT2_LIQUID_STEPS = ("transfer", "distribute", "consolidate", "mix")

_OT2_DEFAULT_API_LEVEL = "2.13"
_OT2_MULTIPART_BOUNDARY = "----SpliceCraftOT2FormBoundary7b3f"


class OT2Error(Exception):
    """Any OT-2 compiler / client failure with a user-facing, actionable message."""


# ── Deck catalog ────────────────────────────────────────────────────────────────
# A curated slice of the Opentrons Labware Library: the gear a SpliceCraft user
# actually reaches for (tip racks, tube racks, 6/24/48/96-well plates,
# reservoirs). Each entry carries a plate FORMAT (rows × cols) so the compiler
# can validate well references and enumerate wells; an unknown load name is not
# fatal (the robot's own analysis is the final authority) — it just skips
# geometry checks and warns. ``kind`` documents intent and lets the validator
# insist a tip rack is present.
#
# Row-major well names: row letter (A, B, …) + 1-based column. A 24-well plate /
# tube rack is 4×6 (A1…D6); a 96 is 8×12 (A1…H12); a 12-channel reservoir is
# 1×12 (A1…A12); a single-well reservoir/trash is 1×1 (A1).
_OT2_LABWARE: "dict[str, dict[str, Any]]" = {
    # tip racks
    "opentrons_96_tiprack_300ul": {"kind": "tiprack", "rows": 8, "cols": 12},
    "opentrons_96_tiprack_20ul": {"kind": "tiprack", "rows": 8, "cols": 12},
    "opentrons_96_tiprack_1000ul": {"kind": "tiprack", "rows": 8, "cols": 12},
    "opentrons_96_filtertiprack_200ul": {"kind": "tiprack", "rows": 8, "cols": 12},
    # tube racks
    "opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap":
        {"kind": "tuberack", "rows": 4, "cols": 6},
    "opentrons_24_tuberack_eppendorf_2ml_safelock_snapcap":
        {"kind": "tuberack", "rows": 4, "cols": 6},
    "opentrons_24_tuberack_generic_2ml_screwcap":
        {"kind": "tuberack", "rows": 4, "cols": 6},
    "opentrons_24_tuberack_nest_1.5ml_snapcap":
        {"kind": "tuberack", "rows": 4, "cols": 6},
    "opentrons_15_tuberack_15ml_conical": {"kind": "tuberack", "rows": 3, "cols": 5},
    "opentrons_6_tuberack_falcon_50ml_conical": {"kind": "tuberack", "rows": 2, "cols": 3},
    # well plates
    "nest_96_wellplate_100ul_pcr_full_skirt": {"kind": "wellplate", "rows": 8, "cols": 12},
    "nest_96_wellplate_2ml_deep": {"kind": "wellplate", "rows": 8, "cols": 12},
    "corning_96_wellplate_360ul_flat": {"kind": "wellplate", "rows": 8, "cols": 12},
    "corning_48_wellplate_1.6ml_flat": {"kind": "wellplate", "rows": 6, "cols": 8},
    "corning_24_wellplate_3.4ml_flat": {"kind": "wellplate", "rows": 4, "cols": 6},
    "corning_12_wellplate_6.9ml_flat": {"kind": "wellplate", "rows": 3, "cols": 4},
    "corning_6_wellplate_16.8ml_flat": {"kind": "wellplate", "rows": 2, "cols": 3},
    # reservoirs
    "nest_12_reservoir_15ml": {"kind": "reservoir", "rows": 1, "cols": 12},
    "nest_1_reservoir_195ml": {"kind": "reservoir", "rows": 1, "cols": 1},
    "agilent_1_reservoir_290ml": {"kind": "reservoir", "rows": 1, "cols": 1},
}

# Terse, friendly aliases → canonical Opentrons load names, so a plan can say
# "plate_24" instead of "corning_24_wellplate_3.4ml_flat". A raw load name is
# passed through untouched.
_OT2_LABWARE_ALIASES: "dict[str, str]" = {
    "tiprack_300": "opentrons_96_tiprack_300ul",
    "tiprack_20": "opentrons_96_tiprack_20ul",
    "tiprack_1000": "opentrons_96_tiprack_1000ul",
    "eppi_24": "opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap",
    "eppi_24_2ml": "opentrons_24_tuberack_eppendorf_2ml_safelock_snapcap",
    "tuberack_24": "opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap",
    "tuberack_15": "opentrons_15_tuberack_15ml_conical",
    "tuberack_6x50": "opentrons_6_tuberack_falcon_50ml_conical",
    "plate_96": "nest_96_wellplate_100ul_pcr_full_skirt",
    "plate_96_flat": "corning_96_wellplate_360ul_flat",
    "plate_96_deep": "nest_96_wellplate_2ml_deep",
    "plate_48": "corning_48_wellplate_1.6ml_flat",
    "plate_24": "corning_24_wellplate_3.4ml_flat",
    "plate_12": "corning_12_wellplate_6.9ml_flat",
    "plate_6": "corning_6_wellplate_16.8ml_flat",
    "reservoir_12": "nest_12_reservoir_15ml",
    "reservoir_1": "nest_1_reservoir_195ml",
}

# Pipette model → (min µL, max µL, channels). The load name is the Protocol-API
# name; GEN1 pipettes have no "_gen2" suffix (Jacques carries a ``p300_single``).
_OT2_PIPETTES: "dict[str, tuple[float, float, int]]" = {
    "p10_single": (1.0, 10.0, 1),
    "p50_single": (5.0, 50.0, 1),
    "p300_single": (30.0, 300.0, 1),
    "p1000_single": (100.0, 1000.0, 1),
    "p20_single_gen2": (1.0, 20.0, 1),
    "p300_single_gen2": (20.0, 300.0, 1),
    "p1000_single_gen2": (100.0, 1000.0, 1),
    "p20_multi_gen2": (1.0, 20.0, 8),
    "p300_multi_gen2": (20.0, 300.0, 8),
}

_OT2_MOUNTS = ("left", "right")
_OT2_NEW_TIP = ("always", "once", "never")
# OT-2 deck: slots 1-11 hold labware; slot 12 is the fixed trash.
_OT2_SLOTS = tuple(range(1, 12))

_ROW_LETTERS = "ABCDEFGHIJKLMNOP"


# ── Well geometry helpers ───────────────────────────────────────────────────────
def _ot2_resolve_labware(name: str) -> str:
    """Expand a friendly alias to its canonical Opentrons load name (pass-through
    for anything already canonical / unknown)."""
    return _OT2_LABWARE_ALIASES.get(name, name)


def _ot2_wells(rows: int, cols: int) -> "list[str]":
    """Every well name of a rows×cols labware, row-major (A1, A2, …, B1, …)."""
    return [f"{_ROW_LETTERS[r]}{c + 1}" for r in range(rows) for c in range(cols)]


def _ot2_parse_well(well: str) -> "tuple[int, int] | None":
    """Split a well like ``B7`` into (row_index, column_number), or ``None`` if it
    is not a well name at all."""
    well = well.strip().upper()
    # `.isdecimal()` not `.isdigit()`: superscripts (²) are isdigit-true but
    # `int()`-invalid, so isdigit would let a bad well past the gate → crash at int().
    if len(well) < 2 or well[0] not in _ROW_LETTERS or not well[1:].isdecimal():
        return None
    return _ROW_LETTERS.index(well[0]), int(well[1:])


def _ot2_well_ok(load_name: str, well: str) -> "bool | None":
    """Is ``well`` valid for this labware? ``True``/``False`` when the geometry is
    known, ``None`` when the load name is unknown (can't tell — defer to the robot)."""
    info = _OT2_LABWARE.get(load_name)
    if info is None:
        return None
    parsed = _ot2_parse_well(well)
    if parsed is None:
        return False
    row, col = parsed
    return row < info["rows"] and 1 <= col <= info["cols"]


def _ot2_custom_wells(definition: Any) -> "set[str] | None":
    """Valid well names of a custom Opentrons labware definition (from its
    ``wells`` map), or ``None`` if it doesn't declare them."""
    if not isinstance(definition, dict):
        return None
    wells = definition.get("wells")
    if isinstance(wells, dict) and wells:
        return {str(k).upper() for k in wells}
    return None


def _ot2_well_ok_entry(entry: "dict[str, Any]", well: str) -> "bool | None":
    """Well validity for a loaded-labware ENTRY — a built-in load name OR a custom
    definition (``{"definition": {...}}``). ``None`` when it can't be determined."""
    definition = entry.get("definition")
    if isinstance(definition, dict):
        cw = _ot2_custom_wells(definition)
        return None if cw is None else well.strip().upper() in cw
    return _ot2_well_ok(str(entry.get("labware", "")), well)


def _ot2_entry_wells(entry: "dict[str, Any]") -> "list[str]":
    """Ordered well names (row-major: A1, A2, …, B1, …) for a deck labware ENTRY —
    from a custom definition's ``wells`` map or the built-in catalog geometry.
    Empty when neither is known. This is the fill order used to lay a collection's
    plasmids onto a plate/rack (identity-linking)."""
    definition = entry.get("definition")
    if isinstance(definition, dict):
        cw = _ot2_custom_wells(definition)
        if cw:
            def _key(w: str) -> "tuple[int, int]":
                p = _ot2_parse_well(w)
                return p if p is not None else (len(_ROW_LETTERS), 0)
            return sorted(cw, key=_key)
    spec = _OT2_LABWARE.get(_ot2_resolve_labware(str(entry.get("labware", ""))))
    if spec and "rows" in spec and "cols" in spec:
        return _ot2_wells(int(spec["rows"]), int(spec["cols"]))
    return []


# ── Plan normalisation ──────────────────────────────────────────────────────────
def _ot2_safe_var(label: str, prefix: str) -> str:
    """A valid, collision-resistant Python identifier for an emitted variable."""
    cleaned = "".join(ch if (ch.isalnum() or ch == "_") else "_" for ch in str(label))
    if not cleaned:
        cleaned = "x"
    if not prefix and cleaned[0].isdigit():
        cleaned = "_" + cleaned
    return prefix + cleaned


def _ot2_split_ref(ref: str) -> "tuple[str, str] | None":
    """Split a ``labwareId:well`` reference into (labware_id, well)."""
    if not isinstance(ref, str) or ":" not in ref:
        return None
    lid, _, well = ref.partition(":")
    lid, well = lid.strip(), well.strip()
    if not lid or not well:
        return None
    return lid, well


# ── Shared validation helpers (used by both the transfer + step paths) ──────────
def _ot2_ref_errors(ref: "tuple[str, str] | None", labware: "dict[str, Any]",
                    tag: str, end: str) -> "list[str]":
    """Validate a parsed ``(labware_id, well)`` reference against loaded labware."""
    if ref is None:
        return [f"{tag}: '{end}' must look like 'labwareId:well'"]
    lid, well = ref
    if lid not in labware:
        return [f"{tag}: '{end}' references unknown labware id {lid!r}"]
    if _ot2_well_ok_entry(labware[lid], well) is False:
        name = labware[lid].get("labware") or "custom labware"
        return [f"{tag}: well {well!r} is out of range for {name!r}"]
    return []


def _ot2_volume_errors(vol: Any, spec: "tuple[float, float, int] | None",
                       pipette: str, tag: str) -> "tuple[list[str], list[str]]":
    """Validate one volume. Returns ``(errors, warnings)``."""
    if not isinstance(vol, (int, float)) or isinstance(vol, bool):
        return [f"{tag}: volume must be a number, got {vol!r}"], []
    if not math.isfinite(vol):
        return [f"{tag}: volume must be a finite number, got {vol!r}"], []
    if vol <= 0:
        return [f"{tag}: volume must be positive, got {vol}"], []
    if vol > _OT2_MAX_VOLUME_UL:
        return [f"{tag}: volume {vol:g} µL is implausibly large "
                f"(max {_OT2_MAX_VOLUME_UL:g})"], []
    if spec is not None and vol < spec[0]:
        return [], [f"{tag}: {vol} µL is below the {pipette} minimum "
                    f"({spec[0]:g} µL) — dispense will be inaccurate"]
    return [], []


def _ot2_mix_errors(mix: Any, spec: "tuple[float, float, int] | None",
                    pipette: str, tag: str) -> "list[str]":
    """Validate a ``[repetitions, volume]`` mix spec (for transfer mix_before/after)."""
    if not (isinstance(mix, (list, tuple)) and len(mix) == 2):
        return [f"{tag}: must be [repetitions, volume]"]
    reps, vol = mix
    errs: "list[str]" = []
    if not isinstance(reps, int) or isinstance(reps, bool) or reps <= 0:
        errs.append(f"{tag}: repetitions must be a positive integer, got {reps!r}")
    errs += _ot2_volume_errors(vol, spec, pipette, tag)[0]
    return errs


_OT2_DEF_MAX_DEPTH = 32


def _ot2_literal_errors(value: Any, tag: str, _depth: int = 0) -> "list[str]":
    """Reject anything ``repr()`` would NOT render as a valid Python literal.

    A custom labware definition is embedded in the generated protocol via
    ``repr()`` — the right tool, since ``json.dumps`` emits ``true`` /
    ``false`` / ``null``, none of which are Python. ``repr()`` is
    literal-safe for every JSON type EXCEPT non-finite floats:
    ``repr(float('inf'))`` is the bare NAME ``inf``, and ``json.loads``
    accepts ``Infinity`` / ``-Infinity`` / ``NaN`` by DEFAULT. So a
    definition carrying one emitted a protocol that ``NameError``s on the
    robot instead of being refused here, where the message is actionable.
    Depth is capped so a pathologically nested definition can't blow the
    recursion limit inside ``repr()`` during emit.
    """
    if _depth > _OT2_DEF_MAX_DEPTH:
        return [f"{tag}: nested more than {_OT2_DEF_MAX_DEPTH} levels deep"]
    if isinstance(value, bool) or value is None or isinstance(value, (int, str)):
        return []
    if isinstance(value, float):
        return ([] if math.isfinite(value) else
                [f"{tag}: {value!r} is not a finite number — Infinity / NaN "
                 f"cannot be written into a protocol"])
    if isinstance(value, (list, tuple)):
        out: "list[str]" = []
        for i, v in enumerate(value):
            out += _ot2_literal_errors(v, f"{tag}[{i}]", _depth + 1)
        return out
    if isinstance(value, dict):
        out = []
        for k, v in value.items():
            if not isinstance(k, (str, int)) or isinstance(k, bool):
                out.append(f"{tag}: key {k!r} must be a string or integer")
                continue
            out += _ot2_literal_errors(v, f"{tag}[{k!r}]", _depth + 1)
        return out
    return [f"{tag}: value of type {type(value).__name__} cannot be written "
            f"into a protocol"]


def _ot2_normalize_step(step: Any, default_new_tip: str) -> "dict[str, Any]":
    """Normalise one designer step into canonical form (parse refs, coerce types).
    Unknown / malformed steps normalise to ``{"type": <as-given>}`` and are flagged
    by ``_ot2_validate_steps``, never silently dropped."""
    if not isinstance(step, dict):
        return {"type": ""}
    stype = str(step.get("type", "")).lower().strip()
    # distribute / consolidate reuse ONE tip by default ("once"); a transfer
    # follows the plan-level default. An explicit step `new_tip` always wins.
    _nt_default = "once" if stype in ("distribute", "consolidate") else default_new_tip
    nt = str(step.get("new_tip", _nt_default)).lower()
    s: "dict[str, Any]" = {"type": stype,
                           "new_tip": nt if nt in _OT2_NEW_TIP else _nt_default}

    def _reflist(v: Any) -> "list[tuple[str, str] | None]":
        if isinstance(v, str):
            v = [v]
        return [_ot2_split_ref(str(x)) for x in v] if isinstance(v, list) else []

    def _mix(v: Any) -> "list[Any] | None":
        if isinstance(v, (list, tuple)) and len(v) == 2:
            return [v[0], v[1]]
        if isinstance(v, dict):
            return [v.get("reps", v.get("repetitions")), v.get("volume", v.get("vol"))]
        return None

    if stype == "transfer":
        s["from"], s["to"] = step.get("from"), step.get("to")
        s["_src"] = _ot2_split_ref(str(step.get("from", "")))
        s["_dst"] = _ot2_split_ref(str(step.get("to", "")))
        s["volume"] = step.get("volume")
        s["mix_before"] = _mix(step.get("mix_before"))
        s["mix_after"] = _mix(step.get("mix_after"))
        s["blow_out"] = bool(step.get("blow_out"))
        s["touch_tip"] = bool(step.get("touch_tip"))
    elif stype == "distribute":
        s["from"], s["to"] = step.get("from"), step.get("to")
        s["_src"] = _ot2_split_ref(str(step.get("from", "")))
        s["_dsts"] = _reflist(step.get("to"))
        s["volume"] = step.get("volume")
    elif stype == "consolidate":
        s["from"], s["to"] = step.get("from"), step.get("to")
        s["_srcs"] = _reflist(step.get("from"))
        s["_dst"] = _ot2_split_ref(str(step.get("to", "")))
        s["volume"] = step.get("volume")
    elif stype == "mix":
        at = step.get("at", step.get("from"))
        s["at"] = at
        s["_at"] = _ot2_split_ref(str(at or ""))
        s["volume"] = step.get("volume")
        s["repetitions"] = step.get("repetitions", step.get("reps"))
    elif stype == "delay":
        s["seconds"] = step.get("seconds")
        s["message"] = step.get("message") or step.get("msg")
    elif stype == "pause":
        s["message"] = step.get("message") or step.get("msg")
    elif stype == "comment":
        s["text"] = str(step.get("text") or step.get("message") or "")
    return s


def _ot2_normalize_plan(plan: "dict[str, Any]") -> "dict[str, Any]":
    """Coerce a raw plan dict into the canonical form the validator + compiler
    consume. Idempotent. Does NOT validate — that's ``_ot2_validate_plan``.

    Accepted input (all but ``labware`` + ``transfers`` optional)::

        {
          "name": "My prep",
          "pipette": "p300_single",          # or {"model": ..., "mount": ...}
          "mount": "left",                    # default "left"
          "tips": {"labware": "tiprack_300", "slot": 8},   # or a list
          "labware": {"src": {"labware": "eppi_24", "slot": 1},
                      "dst": {"labware": "plate_24", "slot": 2}},
          "transfers": [{"from": "src:A1", "to": "dst:A1", "volume": 50}, ...],
          "new_tip": "always",               # always | once | never
          "api_level": "2.13"
        }
    """
    if not isinstance(plan, dict):
        raise OT2Error("plan must be a dict")
    p: "dict[str, Any]" = {}

    pip = plan.get("pipette", "p300_single")
    if isinstance(pip, dict):
        p["pipette"] = str(pip.get("model", "p300_single"))
        mount = pip.get("mount", plan.get("mount", "left"))
    else:
        p["pipette"] = str(pip)
        mount = plan.get("mount", "left")
    p["mount"] = str(mount).lower()

    tips_in = plan.get("tips", [])
    if isinstance(tips_in, dict):
        tips_in = [tips_in]
    if not isinstance(tips_in, list):     # malformed plan: tips is a scalar / junk
        tips_in = []
    tips: "list[dict[str, Any]]" = []
    for t in tips_in:
        if isinstance(t, dict):
            tips.append({"labware": _ot2_resolve_labware(str(t.get("labware", ""))),
                         "slot": t.get("slot")})
    p["tips"] = tips

    labware: "dict[str, dict[str, Any]]" = {}
    lw_map = plan.get("labware")
    if not isinstance(lw_map, dict):      # malformed plan: labware is a list / junk
        lw_map = {}
    for lid, lw in lw_map.items():
        if isinstance(lw, dict):
            entry: "dict[str, Any]" = {
                "labware": _ot2_resolve_labware(str(lw.get("labware", ""))),
                "slot": lw.get("slot")}
            if isinstance(lw.get("definition"), dict):
                entry["definition"] = lw["definition"]   # custom Opentrons labware def
            # Identity-linking metadata (a plate bound to a library collection).
            # Carried through for provenance + round-trip; the compiler ignores it.
            if isinstance(lw.get("map"), dict):
                entry["map"] = {str(k): v for k, v in lw["map"].items()}
            if lw.get("collection"):
                entry["collection"] = str(lw["collection"])
            # Labware-position offset (x/y/z mm) — applied to a run as a
            # labwareOffset; the compiler ignores it (offsets are a run-time input).
            if isinstance(lw.get("offset"), dict):
                entry["offset"] = {k: lw["offset"].get(k) for k in ("x", "y", "z")}
            labware[str(lid)] = entry
    p["labware"] = labware

    transfers: "list[dict[str, Any]]" = []
    transfers_in = plan.get("transfers")
    if not isinstance(transfers_in, list):   # malformed plan: transfers is a scalar
        transfers_in = []
    for t in transfers_in:
        if not isinstance(t, dict):
            continue
        src = _ot2_split_ref(str(t.get("from", "")))
        dst = _ot2_split_ref(str(t.get("to", "")))
        transfers.append({
            "from": t.get("from"), "to": t.get("to"), "volume": t.get("volume"),
            "_src": src, "_dst": dst,
        })
    p["transfers"] = transfers

    steps_in = plan.get("steps")
    nt_default = str(plan.get("new_tip", "always")).lower()
    nt_default = nt_default if nt_default in _OT2_NEW_TIP else "always"
    p["steps"] = ([_ot2_normalize_step(s, nt_default) for s in steps_in]
                  if isinstance(steps_in, list) else [])
    # The multi-step `steps` designer model takes precedence over the legacy
    # `transfers` list whenever the caller supplies it (even an empty list).
    # IDEMPOTENT: a canonical plan already carries `uses_steps` (validate + compile
    # both re-normalise), so trust it — otherwise a legacy transfers plan, whose
    # normalised `steps` is [], would flip into steps mode on the second pass.
    p["uses_steps"] = (bool(plan.get("uses_steps")) if "uses_steps" in plan
                       else isinstance(steps_in, list))

    new_tip = str(plan.get("new_tip", "always")).lower()
    p["new_tip"] = new_tip if new_tip in _OT2_NEW_TIP else "always"
    p["api_level"] = str(plan.get("api_level", _OT2_DEFAULT_API_LEVEL))
    p["metadata"] = {
        "name": str(plan.get("name") or "SpliceCraft transfer"),
        "author": str(plan.get("author") or "SpliceCraft"),
        "description": str(plan.get("description")
                           or "Generated by SpliceCraft from a transfer plan."),
    }
    return p


# ── Validation ──────────────────────────────────────────────────────────────────
def _ot2_validate_transfers(p: "dict[str, Any]", spec: Any,
                            errors: "list[str]", warnings: "list[str]") -> None:
    """Legacy transfer-list validation (a plan with `transfers`, not `steps`)."""
    n = len(p["transfers"])
    if n == 0:
        errors.append("no transfers: add at least one entry to 'transfers'")
    elif n > _OT2_MAX_TRANSFERS:
        errors.append(f"too many transfers ({n}); the max is {_OT2_MAX_TRANSFERS}")
    # Capped so a pathological count can't make validation itself O(huge).
    for i, t in enumerate(p["transfers"][:_OT2_MAX_TRANSFERS], 1):
        tag = f"transfer #{i}"
        errors += _ot2_ref_errors(t["_src"], p["labware"], tag, "from")
        errors += _ot2_ref_errors(t["_dst"], p["labware"], tag, "to")
        verr, vwarn = _ot2_volume_errors(t.get("volume"), spec, p["pipette"], tag)
        errors += verr
        warnings += vwarn
    if p["new_tip"] == "always":
        tip_capacity = 96 * len([t for t in p["tips"] if t.get("labware")])
        if tip_capacity and n > tip_capacity:
            warnings.append(f"{n} transfers with new_tip='always' need {n} tips "
                            f"but only {tip_capacity} are loaded — add tip racks or use "
                            "new_tip='once'")


def _ot2_validate_steps(p: "dict[str, Any]", spec: Any,
                        errors: "list[str]", warnings: "list[str]") -> None:
    """Multi-step designer validation (a plan with `steps`)."""
    n = len(p["steps"])
    if n == 0:
        errors.append("no steps: add at least one entry to 'steps'")
    elif n > _OT2_MAX_STEPS:
        errors.append(f"too many steps ({n}); the max is {_OT2_MAX_STEPS}")
    lw, pip = p["labware"], p["pipette"]
    for i, s in enumerate(p["steps"][:_OT2_MAX_STEPS], 1):
        st = s["type"]
        tag = f"step #{i} ({st or '?'})"
        if st not in _OT2_STEP_TYPES:
            errors.append(f"step #{i}: unknown step type {st!r} "
                          f"(known: {', '.join(_OT2_STEP_TYPES)})")
            continue
        if st == "transfer":
            errors += _ot2_ref_errors(s["_src"], lw, tag, "from")
            errors += _ot2_ref_errors(s["_dst"], lw, tag, "to")
            verr, vwarn = _ot2_volume_errors(s.get("volume"), spec, pip, tag)
            errors += verr
            warnings += vwarn
            for mk in ("mix_before", "mix_after"):
                if s.get(mk) is not None:
                    errors += _ot2_mix_errors(s[mk], spec, pip, f"{tag} {mk}")
        elif st == "distribute":
            errors += _ot2_ref_errors(s["_src"], lw, tag, "from")
            if not s["_dsts"]:
                errors.append(f"{tag}: 'to' must be a non-empty list of wells")
            for j, r in enumerate(s["_dsts"], 1):
                errors += _ot2_ref_errors(r, lw, f"{tag} to[{j}]", "to")
            verr, vwarn = _ot2_volume_errors(s.get("volume"), spec, pip, tag)
            errors += verr
            warnings += vwarn
        elif st == "consolidate":
            if not s["_srcs"]:
                errors.append(f"{tag}: 'from' must be a non-empty list of wells")
            for j, r in enumerate(s["_srcs"], 1):
                errors += _ot2_ref_errors(r, lw, f"{tag} from[{j}]", "from")
            errors += _ot2_ref_errors(s["_dst"], lw, tag, "to")
            verr, vwarn = _ot2_volume_errors(s.get("volume"), spec, pip, tag)
            errors += verr
            warnings += vwarn
        elif st == "mix":
            errors += _ot2_ref_errors(s["_at"], lw, tag, "at")
            verr, vwarn = _ot2_volume_errors(s.get("volume"), spec, pip, tag)
            errors += verr
            warnings += vwarn
            reps = s.get("repetitions")
            if not isinstance(reps, int) or isinstance(reps, bool) or reps <= 0:
                errors.append(f"{tag}: repetitions must be a positive integer, got {reps!r}")
        elif st == "delay":
            secs = s.get("seconds")
            if (not isinstance(secs, (int, float)) or isinstance(secs, bool)
                    or not math.isfinite(secs) or secs <= 0):
                errors.append(f"{tag}: seconds must be a positive number, got {secs!r}")
            elif secs > _OT2_MAX_DELAY_S:
                errors.append(f"{tag}: delay {secs:g}s is implausibly long "
                              f"(max {_OT2_MAX_DELAY_S:g}s)")
        elif st == "comment":
            if not str(s.get("text", "")).strip():
                warnings.append(f"{tag}: empty comment")
        # pause: an optional message — nothing to validate


def _ot2_validate_plan(plan: "dict[str, Any]") -> "dict[str, list[str]]":
    """Validate a plan (raw or canonical). Returns ``{"errors": [...],
    "warnings": [...]}``. Errors block compilation; warnings don't (e.g. a
    sub-minimum volume the pipette can't dispense accurately, or an unrecognised
    load name the robot will have the final say on). Handles both the legacy
    `transfers` list and the multi-step `steps` designer model."""
    p = _ot2_normalize_plan(plan)
    errors: "list[str]" = []
    warnings: "list[str]" = []

    spec = _OT2_PIPETTES.get(p["pipette"])
    if spec is None:
        errors.append(f"unknown pipette model {p['pipette']!r} "
                      f"(known: {', '.join(sorted(_OT2_PIPETTES))})")
    if p["mount"] not in _OT2_MOUNTS:
        errors.append(f"mount must be 'left' or 'right', got {p['mount']!r}")

    used_slots: "dict[int, str]" = {}

    def _check_slot(slot: Any, what: str) -> None:
        if not isinstance(slot, int) or isinstance(slot, bool):
            errors.append(f"{what}: slot must be an integer 1-11, got {slot!r}")
            return
        if slot not in _OT2_SLOTS:
            errors.append(f"{what}: slot {slot} is not a labware slot "
                          "(use 1-11; slot 12 is the fixed trash)")
            return
        if slot in used_slots:
            errors.append(f"{what}: slot {slot} already holds {used_slots[slot]}")
        else:
            used_slots[slot] = what

    # A protocol of only control steps (delay / pause / comment) needs no tips or
    # labware; require them only when a liquid-handling operation is present.
    needs_liquid = (any(s["type"] in _OT2_LIQUID_STEPS for s in p["steps"])
                    if p["uses_steps"] else bool(p["transfers"]))
    if needs_liquid and not p["tips"]:
        errors.append("no tip rack: add at least one entry to 'tips'")
    for t in p["tips"]:
        _check_slot(t.get("slot"), f"tip rack {t.get('labware')!r}")
        if t.get("labware") and t["labware"] not in _OT2_LABWARE:
            warnings.append(f"tip rack {t['labware']!r} not in the built-in catalog "
                            "— the robot's analysis will validate it")
        elif t.get("labware") and _OT2_LABWARE[t["labware"]]["kind"] != "tiprack":
            errors.append(f"{t['labware']!r} is not a tip rack")

    if needs_liquid and not p["labware"]:
        errors.append("no labware: add at least one entry to 'labware'")
    for lid, lw in p["labware"].items():
        _check_slot(lw.get("slot"), f"labware {lid!r}")
        if isinstance(lw.get("definition"), dict):
            # The definition is `repr()`-embedded into the emitted protocol —
            # refuse anything that wouldn't be a valid Python literal BEFORE
            # it reaches the compiler.
            errors += _ot2_literal_errors(lw["definition"],
                                          f"custom labware (id {lid!r}) definition")
            if _ot2_custom_wells(lw["definition"]) is None:
                warnings.append(f"custom labware (id {lid!r}) declares no 'wells' — well "
                                "checks skipped; the robot will validate it")
        elif lw.get("labware") and lw["labware"] not in _OT2_LABWARE:
            warnings.append(f"labware {lw['labware']!r} (id {lid!r}) not in the built-in "
                            "catalog — well checks skipped; the robot will validate it")

    if p["uses_steps"]:
        _ot2_validate_steps(p, spec, errors, warnings)
    else:
        _ot2_validate_transfers(p, spec, errors, warnings)
    return {"errors": errors, "warnings": warnings}


# ── The compiler ────────────────────────────────────────────────────────────────
def _ot2_emit_step(s: "dict[str, Any]", id_var: "dict[str, str]") -> "list[str]":
    """Emit the Protocol API v2 line(s) for one validated designer step."""
    st = s["type"]

    def well(ref: "tuple[str, str]") -> str:
        return f"{id_var[ref[0]]}[{json.dumps(ref[1])}]"

    def wells(refs: "list[Any]") -> str:
        return "[" + ", ".join(well(r) for r in refs) + "]"

    kw: "list[str]" = []
    if st in ("transfer", "distribute", "consolidate"):
        kw.append(f"new_tip={json.dumps(s['new_tip'])}")
    if st == "transfer":
        if s.get("mix_before"):
            kw.append(f"mix_before=({s['mix_before'][0]!r}, {s['mix_before'][1]!r})")
        if s.get("mix_after"):
            kw.append(f"mix_after=({s['mix_after'][0]!r}, {s['mix_after'][1]!r})")
        if s.get("blow_out"):
            kw.append("blow_out=True")
            kw.append('blowout_location="destination well"')
        if s.get("touch_tip"):
            kw.append("touch_tip=True")
        opts = "".join(f", {k}" for k in kw)
        return [f"    pipette.transfer({s['volume']!r}, {well(s['_src'])}, "
                f"{well(s['_dst'])}{opts})"]
    if st == "distribute":
        opts = "".join(f", {k}" for k in kw)
        return [f"    pipette.distribute({s['volume']!r}, {well(s['_src'])}, "
                f"{wells(s['_dsts'])}{opts})"]
    if st == "consolidate":
        opts = "".join(f", {k}" for k in kw)
        return [f"    pipette.consolidate({s['volume']!r}, {wells(s['_srcs'])}, "
                f"{well(s['_dst'])}{opts})"]
    if st == "mix":
        # A standalone mix needs a tip: pick up, mix, drop (self-contained).
        return ["    pipette.pick_up_tip()",
                f"    pipette.mix({s['repetitions']!r}, {s['volume']!r}, {well(s['_at'])})",
                "    pipette.drop_tip()"]
    if st == "delay":
        msg = f", msg={json.dumps(s['message'])}" if s.get("message") else ""
        return [f"    protocol.delay(seconds={s['seconds']!r}{msg})"]
    if st == "pause":
        return [f"    protocol.pause({json.dumps(s['message']) if s.get('message') else ''})"]
    if st == "comment":
        return [f"    protocol.comment({json.dumps(s['text'])})"]
    return []


def _ot2_compile_protocol(plan: "dict[str, Any]") -> str:
    """Compile a plan into Opentrons Protocol API v2 Python text.

    Handles both the legacy `transfers` list and the multi-step `steps` designer
    model. Raises ``OT2Error`` (listing every problem) if the plan does not
    validate — a caller should surface those rather than ship a broken protocol.
    """
    p = _ot2_normalize_plan(plan)
    report = _ot2_validate_plan(p)
    if report["errors"]:
        raise OT2Error("invalid plan:\n  - " + "\n  - ".join(report["errors"]))

    md = p["metadata"]
    out: "list[str]" = [
        "from opentrons import protocol_api",
        "",
        "metadata = {",
        f"    \"protocolName\": {json.dumps(md['name'])},",
        f"    \"author\": {json.dumps(md['author'])},",
        f"    \"description\": {json.dumps(md['description'])},",
        f"    \"apiLevel\": {json.dumps(p['api_level'])},",
        "}",
        "",
        "# Generated by SpliceCraft. Validated against the built-in deck catalog",
        "# before emit; ALWAYS re-check with the robot's built-in analysis (or the",
        "# Opentrons App's simulate) before running on hardware.",
        "",
        "def run(protocol: protocol_api.ProtocolContext):",
    ]

    tip_vars: "list[str]" = []
    for t in p["tips"]:
        var = _ot2_safe_var(str(t["slot"]), "tips_")
        tip_vars.append(var)
        out.append(f"    {var} = protocol.load_labware("
                   f"{json.dumps(t['labware'])}, {t['slot']})")

    id_var: "dict[str, str]" = {}
    used_vars = set(tip_vars)
    for lid, lw in p["labware"].items():
        var = _ot2_safe_var(lid, "lw_")
        if var in used_vars:          # distinct ids can sanitise to the same name
            k = 2
            while f"{var}_{k}" in used_vars:
                k += 1
            var = f"{var}_{k}"
        used_vars.add(var)
        id_var[lid] = var
        if isinstance(lw.get("definition"), dict):
            # custom labware: embed the Opentrons definition (a plain dict → a
            # valid Python literal) and load it from that.
            out.append(f"    {var} = protocol.load_labware_from_definition("
                       f"{lw['definition']!r}, {lw['slot']})")
        else:
            out.append(f"    {var} = protocol.load_labware("
                       f"{json.dumps(lw['labware'])}, {lw['slot']})")

    out.append(
        f"    pipette = protocol.load_instrument("
        f"{json.dumps(p['pipette'])}, {json.dumps(p['mount'])}, "
        f"tip_racks=[{', '.join(tip_vars)}])"
    )
    out.append("    protocol.home()")

    if p["uses_steps"]:
        for s in p["steps"]:
            out.extend(_ot2_emit_step(s, id_var))
    else:
        volumes = [t["volume"] for t in p["transfers"]]
        srcs = [f"{id_var[t['_src'][0]]}[{json.dumps(t['_src'][1])}]" for t in p["transfers"]]
        dsts = [f"{id_var[t['_dst'][0]]}[{json.dumps(t['_dst'][1])}]" for t in p["transfers"]]
        out.append("    pipette.transfer(")
        out.append(f"        {volumes!r},")
        out.append("        [" + ", ".join(srcs) + "],")
        out.append("        [" + ", ".join(dsts) + "],")
        out.append(f"        new_tip={json.dumps(p['new_tip'])},")
        out.append("    )")
    return "\n".join(out) + "\n"


def _ot2_plan_summary(plan: "dict[str, Any]") -> "dict[str, Any]":
    """A compact, human/agent-friendly digest of a plan (counts, volume, tips)."""
    p = _ot2_normalize_plan(plan)
    report = _ot2_validate_plan(p)

    def _num(v: Any) -> bool:
        return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)

    if p["uses_steps"]:
        liquid = [s for s in p["steps"] if s["type"] in _OT2_LIQUID_STEPS]
        total_vol = sum(s.get("volume") for s in liquid if _num(s.get("volume")))
        n_steps = len(p["steps"])
        n_transfers = len([s for s in p["steps"] if s["type"] == "transfer"])
        tips_needed = len(liquid)   # each liquid step uses ~1 tip
        step_types: "dict[str, int]" = {}
        for s in p["steps"]:
            step_types[s["type"]] = step_types.get(s["type"], 0) + 1
    else:
        total_vol = sum(t["volume"] for t in p["transfers"] if _num(t["volume"]))
        n_steps = len(p["transfers"])
        n_transfers = len(p["transfers"])
        tips_needed = (n_transfers if p["new_tip"] == "always"
                       else (1 if p["new_tip"] == "once" and p["transfers"] else 0))
        step_types = {"transfer": n_transfers} if n_transfers else {}
    return {
        "name": p["metadata"]["name"],
        "pipette": p["pipette"],
        "mount": p["mount"],
        "labware": {lid: (lw.get("labware") or "custom")
                    for lid, lw in p["labware"].items()},
        "steps": n_steps,
        "transfers": n_transfers,
        "step_types": step_types,
        "total_volume_ul": total_vol,
        "new_tip": p["new_tip"],
        "tips_needed": tips_needed,
        # Pre-flight: how much liquid each source well must hold (so you know what
        # to load and whether a well will run dry) + a rough wall-clock estimate.
        "source_volumes": _ot2_source_volumes(p),
        "est_seconds": round(_ot2_estimate_seconds(p), 1),
        "valid": not report["errors"],
        "errors": report["errors"],
        "warnings": report["warnings"],
    }


# ── Resource pre-flight ─────────────────────────────────────────────────────────
# Rough time model: a nominal single-channel flow rate + fixed per-move / per-tip
# overheads. Wall-clock on a real OT-2 varies with speeds + travel, so this is a
# ballpark ("~4 min"), not a promise.
_OT2_TIME_TIP_S = 7.0          # pick up + drop one tip
_OT2_TIME_MOVE_S = 4.0         # one aspirate or dispense (approach + plunger)
_OT2_FLOW_UL_PER_S = 150.0     # nominal single-channel flow


def _ot2_num(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def _ot2_source_volumes(p: "dict[str, Any]") -> "dict[str, float]":
    """Total µL drawn from each source well across the plan — the load budget.
    Keyed by ``labwareId:well``. A mix returns liquid to its own well, so it draws
    nothing net and is skipped."""
    draws: "dict[str, float]" = {}

    def _add(ref: "tuple[str, str] | None", vol: Any) -> None:
        if ref is None or not _ot2_num(vol) or vol <= 0:
            return
        key = f"{ref[0]}:{ref[1]}"
        draws[key] = round(draws.get(key, 0.0) + float(vol), 4)

    if p["uses_steps"]:
        for s in p["steps"]:
            st = s["type"]
            if st == "transfer":
                _add(s.get("_src"), s.get("volume"))
            elif st == "distribute":
                n = len([r for r in (s.get("_dsts") or []) if r])
                _add(s.get("_src"), (s.get("volume") or 0) * n if _ot2_num(s.get("volume")) else None)
            elif st == "consolidate":
                for r in (s.get("_srcs") or []):
                    _add(r, s.get("volume"))
    else:
        for t in p["transfers"]:
            _add(t.get("_src"), t.get("volume"))
    return draws


def _ot2_estimate_seconds(p: "dict[str, Any]") -> float:
    """A rough wall-clock estimate for the plan (see the time model above)."""
    total = 0.0

    def _liquid(vol: Any, moves: int, tips: float) -> float:
        v = float(vol) if _ot2_num(vol) and vol > 0 else 0.0
        return tips * _OT2_TIME_TIP_S + moves * (_OT2_TIME_MOVE_S + v / _OT2_FLOW_UL_PER_S)

    if p["uses_steps"]:
        for s in p["steps"]:
            st = s["type"]
            tips = 0.0 if s.get("new_tip") == "never" else 1.0
            if st == "transfer":
                total += _liquid(s.get("volume"), 2, tips)
            elif st == "distribute":
                total += _liquid(s.get("volume"), 1 + len([r for r in (s.get("_dsts") or []) if r]), tips)
            elif st == "consolidate":
                total += _liquid(s.get("volume"), len([r for r in (s.get("_srcs") or []) if r]) + 1, tips)
            elif st == "mix":
                reps = s.get("repetitions") if isinstance(s.get("repetitions"), int) else 1
                total += _OT2_TIME_TIP_S + max(1, reps) * (_OT2_TIME_MOVE_S)
            elif st == "delay":
                secs = s.get("seconds")
                total += float(secs) if _ot2_num(secs) and secs > 0 else 0.0
    else:
        for i, t in enumerate(p["transfers"]):
            tips = 1.0 if (p["new_tip"] == "always" or (p["new_tip"] == "once" and i == 0)) else 0.0
            total += _liquid(t.get("volume"), 2, tips)
    return total


# ── Concentration normalisation ─────────────────────────────────────────────────
def _ot2_normalize_volumes(items: "list[dict[str, Any]]", *,
                           target_ng: "float | None" = None,
                           target_conc: "float | None" = None,
                           final_volume: "float | None" = None,
                           min_vol: float = 0.0,
                           max_vol: "float | None" = None,
                           resolution: float = 0.1) -> "list[dict[str, Any]]":
    """Per-item sample (+ optional diluent) volumes to normalise concentration.
    Pure — no robot, no plan.

    Each ``item``: ``{"name", "concentration" (ng/µL), "well"?}``.

    Exactly one target mode:
      * ``target_ng`` — deliver a fixed MASS: ``sample = target_ng / conc``. With
        ``final_volume`` set, the rest is diluent (top-up to that volume).
      * ``target_conc`` + ``final_volume`` — deliver ``final_volume`` at
        ``target_conc``: ``sample = target_conc·final_volume / conc``,
        ``diluent = final_volume − sample``.

    Volumes round to ``resolution`` and clamp into ``[min_vol, max_vol]`` (when
    given); an item too dilute to reach target within ``max_vol`` (or so
    concentrated it needs less than ``min_vol``) is FLAGGED + clamped, never
    silently dropped. Returns a per-item list with ``sample_ul`` / ``diluent_ul``
    / ``achieved_ng`` / ``achieved_conc`` / ``ok`` / ``warning``. Raises
    ``OT2Error`` on a contradictory or non-finite request."""
    def _round(v: float) -> float:
        if resolution and resolution > 0:
            return round(round(v / resolution) * resolution, 6)
        return round(v, 6)

    def _pos(v: Any, label: str) -> float:
        if not (isinstance(v, (int, float)) and not isinstance(v, bool)
                and math.isfinite(v) and v > 0):
            raise OT2Error(f"normalise: {label} must be a positive number, got {v!r}")
        return float(v)

    has_ng, has_conc = target_ng is not None, target_conc is not None
    if has_ng == has_conc:
        raise OT2Error("normalise: pass exactly one of 'target_ng' or 'target_conc'")
    tgt_ng = _pos(target_ng, "target_ng") if has_ng else 0.0
    tgt_conc = _pos(target_conc, "target_conc") if has_conc else 0.0
    fv: "float | None" = (_pos(final_volume, "final_volume")
                          if final_volume is not None else None)
    if has_conc and fv is None:
        raise OT2Error("normalise: 'target_conc' mode needs a positive 'final_volume'")
    mv: "float | None" = _pos(max_vol, "max_vol") if max_vol is not None else None

    def _nonneg(v: Any, label: str) -> float:
        if not (isinstance(v, (int, float)) and not isinstance(v, bool)
                and math.isfinite(v) and v >= 0):
            raise OT2Error(f"normalise: {label} must be a non-negative number, got {v!r}")
        return float(v)

    min_vol = _nonneg(min_vol, "min_vol")
    resolution = _nonneg(resolution, "resolution")
    if mv is not None and min_vol > mv:
        raise OT2Error(f"normalise: min_vol ({min_vol:g}) exceeds max_vol ({mv:g})")
    if not isinstance(items, list):
        raise OT2Error("normalise: 'items' must be a list")
    if len(items) > _OT2_MAX_ITEMS:
        raise OT2Error(f"normalise: too many items ({len(items)}); "
                       f"the max is {_OT2_MAX_ITEMS}")

    out: "list[dict[str, Any]]" = []
    for it in items:
        if not isinstance(it, dict):
            continue
        conc = it.get("concentration", it.get("conc"))
        rec: "dict[str, Any]" = {
            "name": str(it.get("name", "") or ""), "well": it.get("well"),
            "concentration": conc, "sample_ul": None, "diluent_ul": None,
            "achieved_ng": None, "achieved_conc": None, "ok": False, "warning": None}
        if not (isinstance(conc, (int, float)) and not isinstance(conc, bool)
                and math.isfinite(conc) and conc > 0):
            rec["warning"] = "no / invalid concentration — skipped"
            out.append(rec)
            continue
        conc = float(conc)
        sample = (tgt_ng / conc) if has_ng else (tgt_conc * (fv or 0.0) / conc)
        warn: "str | None" = None
        if mv is not None and sample > mv:
            warn = (f"needs {sample:.1f} µL to hit target but the max is "
                    f"{mv:g} µL — clamped (under-target)")
            sample = mv
        if sample < min_vol:
            if min_vol > 0:
                warn = (f"only {sample:.2f} µL needed (below the {min_vol:g} µL "
                        "pipette floor) — dilute the stock first")
            sample = max(sample, min_vol)
        sample_r = _round(sample)
        if sample_r <= 0:
            rec["warning"] = warn or "computed volume rounds to zero — skipped"
            out.append(rec)
            continue
        diluent_r = _round(max(0.0, fv - sample_r)) if fv is not None else None
        if fv is not None and sample_r > fv and warn is None:
            warn = (f"sample {sample_r:g} µL exceeds the {fv:g} µL final volume — "
                    "stock too dilute for this target; the well is not diluted")
        mass = sample_r * conc
        total = sample_r + (diluent_r or 0.0)
        rec.update(sample_ul=sample_r, diluent_ul=diluent_r,
                   achieved_ng=round(mass, 3),
                   achieved_conc=(round(mass / total, 4) if total > 0 else None),
                   ok=True, warning=warn)
        out.append(rec)
    return out


def _ot2_normalize_steps(normalized: "list[dict[str, Any]]", *, src_id: str,
                         dst_id: str, dst_wells: "list[str]",
                         diluent_ref: "str | None" = None,
                         new_tip: str = "always") -> "list[dict[str, Any]]":
    """Turn ``_ot2_normalize_volumes`` output into transfer steps: for each OK item,
    diluent first (from ``diluent_ref`` if given) then sample (from ``src_id:well``)
    into the next destination well. Items without a source well, or once the
    destination wells are exhausted, are skipped."""
    steps: "list[dict[str, Any]]" = []
    di = 0
    for rec in normalized:
        if not rec.get("ok"):
            continue
        well = rec.get("well")
        if not well or di >= len(dst_wells):
            continue
        dst = f"{dst_id}:{dst_wells[di]}"
        di += 1
        dil = rec.get("diluent_ul")
        if (diluent_ref and isinstance(dil, (int, float))
                and not isinstance(dil, bool) and dil > 0):
            steps.append({"type": "transfer", "new_tip": new_tip,
                          "from": diluent_ref, "to": dst, "volume": dil})
        steps.append({"type": "transfer", "new_tip": new_tip,
                      "from": f"{src_id}:{well}", "to": dst, "volume": rec["sample_ul"]})
    return steps


def _ot2_cherrypick_steps(picks: "list[dict[str, Any]]", *, src_id: str, dst_id: str,
                          dst_wells: "list[str]", volume: Any,
                          new_tip: str = "always") -> "list[dict[str, Any]]":
    """Transfer ``volume`` from each pick's source well (``src_id:well``) into
    sequential destination wells — the equal-volume cherry-pick / replate."""
    steps: "list[dict[str, Any]]" = []
    di = 0
    for p in picks:
        well = p.get("well") if isinstance(p, dict) else None
        if not well or di >= len(dst_wells):
            continue
        steps.append({"type": "transfer", "new_tip": new_tip,
                      "from": f"{src_id}:{well}", "to": f"{dst_id}:{dst_wells[di]}",
                      "volume": volume})
        di += 1
    return steps


# ── Deck visualizer ─────────────────────────────────────────────────────────────
# The OT-2 deck, physical layout (front row 1-2-3 at the bottom, trash at 12):
#     10  11  12(trash)
#      7   8   9
#      4   5   6
#      1   2   3
_OT2_DECK_LAYOUT = [[10, 11, 12], [7, 8, 9], [4, 5, 6], [1, 2, 3]]

# A single deck slot's physical footprint (mm): the ANSI/SLAS microplate standard,
# long axis left-right. Width:depth ≈ 1.49 — a slot is a landscape rectangle, not a
# square. The Textual deck map (OT2DeckMap) uses this to draw bays at the real
# aspect ratio rather than as squished slivers.
_OT2_SLOT_FOOTPRINT_MM = (127.76, 85.48)


def _ot2_short_labware(name: str) -> str:
    """A compact display label for a labware load name — its friendly alias if
    one exists, else a truncated load name."""
    if not name:
        return ""
    rev = {v: k for k, v in _OT2_LABWARE_ALIASES.items()}
    return (rev.get(name) or name)[:13]


def _ot2_render_deck(plan: "dict[str, Any]") -> str:
    """Render the OT-2 deck as a Unicode grid (slots 1-11 + the fixed trash at 12),
    showing which labware sits in each slot — a text 'deck map' à la the Opentrons
    Protocol Designer."""
    p = _ot2_normalize_plan(plan)
    occ: "dict[int, tuple[str, str]]" = {}
    for t in p["tips"]:
        s = t.get("slot")
        if isinstance(s, int) and not isinstance(s, bool):
            occ[s] = ("tips", _ot2_short_labware(str(t.get("labware", ""))))
    for lid, lw in p["labware"].items():
        s = lw.get("slot")
        if isinstance(s, int) and not isinstance(s, bool):
            label = "custom" if lw.get("definition") else _ot2_short_labware(str(lw.get("labware", "")))
            occ[s] = (lid, label)
    w = 14

    def clip(s: str) -> str:
        return (" " + s)[:w].ljust(w)

    def cell(slot: int) -> "list[str]":
        if slot == 12:
            return [clip("12  TRASH"), clip(""), clip("fixed trash")]
        who, name = occ.get(slot, ("", ""))
        return [clip(str(slot)), clip(who), clip(name)]

    def hbar(left: str, mid: str, right: str) -> str:
        return left + ("─" * w + mid) * 2 + "─" * w + right

    lines = [hbar("┌", "┬", "┐")]
    for ri, row in enumerate(_OT2_DECK_LAYOUT):
        cells = [cell(s) for s in row]
        for li in range(3):
            lines.append("│" + "│".join(cells[c][li] for c in range(3)) + "│")
        lines.append(hbar("├", "┼", "┤") if ri < len(_OT2_DECK_LAYOUT) - 1
                     else hbar("└", "┴", "┘"))
    return "\n".join(lines)


def _ot2_build_labware_def(name: str, rows: int, cols: int, *,
                           category: str = "wellPlate", spacing: float = 9.0,
                           x_off: float = 14.38, y_off: float = 11.24,
                           diameter: float = 6.5, depth: float = 14.0,
                           volume: float = 200.0, x_dim: float = 127.76,
                           y_dim: float = 85.48, z_dim: float = 15.0) -> "dict[str, Any]":
    """Generate a structurally-valid Opentrons labware definition for a REGULAR
    rows×cols grid (well plate / tube rack / reservoir), with SLAS-footprint
    defaults. A1 is back-left; wells march right (columns) and toward the front
    (rows). ALWAYS analyse on the robot before running — the geometry defaults
    are generic, not a substitute for the official Labware Creator's calibration."""
    def _fin(v: Any, default: float, lo: float, hi: float) -> float:
        """Clamp to a finite value in [lo, hi]; reject inf/nan/junk (guards
        int(inf)/int(nan) crashes + inf/nan well coordinates)."""
        try:
            v = float(v)
        except (TypeError, ValueError):
            return default
        return max(lo, min(v, hi)) if math.isfinite(v) else default

    rows = max(1, min(int(_fin(rows, 1, 1, len(_ROW_LETTERS))), len(_ROW_LETTERS)))
    cols = max(1, min(int(_fin(cols, 1, 1, 99)), 99))
    spacing = _fin(spacing, 9.0, 0.1, 100.0)
    x_off, y_off = _fin(x_off, 14.38, 0.0, 200.0), _fin(y_off, 11.24, 0.0, 200.0)
    diameter, depth = _fin(diameter, 6.5, 0.1, 100.0), _fin(depth, 14.0, 0.1, 300.0)
    volume = _fin(volume, 200.0, 0.0, 1e7)
    x_dim, y_dim = _fin(x_dim, 127.76, 1.0, 1000.0), _fin(y_dim, 85.48, 1.0, 1000.0)
    z_dim = _fin(z_dim, 15.0, 1.0, 1000.0)
    ordering: "list[list[str]]" = []
    wells: "dict[str, Any]" = {}
    for c in range(cols):
        col: "list[str]" = []
        for r in range(rows):
            wn = f"{_ROW_LETTERS[r]}{c + 1}"
            col.append(wn)
            wells[wn] = {"depth": depth, "totalLiquidVolume": volume,
                         "shape": "circular", "diameter": diameter,
                         "x": round(x_off + c * spacing, 2),
                         "y": round(y_dim - y_off - r * spacing, 2),
                         "z": round(z_dim - depth, 2)}
        ordering.append(col)
    load_name = ("".join(ch if ch.isalnum() else "_" for ch in name.lower()).strip("_")
                 or "custom_labware")
    return {
        "ordering": ordering,
        "brand": {"brand": "SpliceCraft", "brandId": []},
        "metadata": {"displayName": name, "displayCategory": category,
                     "displayVolumeUnits": "µL", "tags": []},
        "dimensions": {"xDimension": x_dim, "yDimension": y_dim, "zDimension": z_dim},
        "wells": wells,
        "groups": [{"metadata": {"wellBottomShape": "flat"}, "wells": list(wells)}],
        "parameters": {"format": "irregular", "quirks": [],
                       "isMagneticModuleCompatible": False,
                       "loadName": load_name, "isTiprack": False},
        "namespace": "custom_beta", "version": 1, "schemaVersion": 2,
        "cornerOffsetFromSlot": {"x": 0.0, "y": 0.0, "z": 0.0},
    }


# ── LAN client (OT-2 robot-server, port 31950) ──────────────────────────────────
def _ot2_user_agent() -> str:
    return f"SpliceCraft/{_state._sc_version or '?'} (OT-2 client)"


def _ot2_base_url(host: str) -> str:
    """Normalise a user-supplied host into ``http://HOST:31950``.

    Accepts a bare IP / hostname (``192.168.1.56``, ``opentrons.local``), an
    ``host:port`` pair, or a full ``http://…`` URL. Only plain HTTP to the
    robot-server is supported.
    """
    host = (host or "").strip()
    if not host:
        raise OT2Error("no OT-2 host given")
    if "://" in host:
        try:
            parts = urllib.parse.urlsplit(host)
        except ValueError as exc:   # e.g. an unbalanced IPv6 bracket "http://["
            raise OT2Error(f"invalid host {host!r}: {exc}") from exc
        if parts.scheme not in ("http", ""):
            raise OT2Error(f"unsupported scheme {parts.scheme!r} — the OT-2 API is HTTP")
        netloc = parts.netloc or parts.path
    else:
        netloc = host
    if ":" not in netloc:
        netloc = f"{netloc}:{_OT2_PORT}"
    return f"http://{netloc}"


class _OT2NoRedirect(urllib.request.HTTPRedirectHandler):
    """The robot-server never legitimately redirects an API call; treat a 3xx as
    an error rather than silently following it off-box."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[override]
        raise urllib.error.HTTPError(
            newurl, code, f"unexpected redirect to {newurl}", headers, fp)


def _ot2_opener() -> "urllib.request.OpenerDirector":
    # A plain opener (NOT splicecraft_net's public-host-asserting hardened one):
    # the OT-2 lives on a private LAN address the hardened opener refuses. Same
    # deliberate local-service bypass the BABS Ollama client uses.
    return urllib.request.build_opener(_OT2NoRedirect())


def _ot2_read_capped(resp: Any, max_bytes: int) -> bytes:
    data = resp.read(max_bytes + 1)
    if len(data) > max_bytes:
        raise OT2Error(f"OT-2 response exceeded the {max_bytes}-byte cap")
    return data


def _ot2_request_json(host: str, path: str, *, method: str = "GET",
                      payload: "dict[str, Any] | None" = None,
                      timeout: float = _OT2_SHORT_TIMEOUT,
                      max_bytes: int = _OT2_JSON_MAX_BYTES) -> "dict[str, Any]":
    """One JSON call to the robot-server. Raises ``OT2Error`` on any transport /
    HTTP / decode failure with an actionable message."""
    url = _ot2_base_url(host) + path
    headers = {"Opentrons-Version": _OT2_API_HEADER, "User-Agent": _ot2_user_agent()}
    data = None
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with _ot2_opener().open(req, timeout=timeout) as resp:
            body = _ot2_read_capped(resp, max_bytes)
    except urllib.error.HTTPError as exc:
        detail = ""
        try:
            detail = exc.read(64 * 1024).decode("utf-8", "replace")
            exc.close()          # release the socket promptly (up to 1024×/sweep)
        except Exception:
            pass
        raise OT2Error(f"OT-2 {method} {path} → HTTP {exc.code}"
                       + (f": {detail[:400]}" if detail else "")) from exc
    except (urllib.error.URLError, OSError) as exc:
        raise OT2Error(f"cannot reach OT-2 at {host} ({exc}). Is it powered on and "
                       "on this network?") from exc
    try:
        return json.loads(body.decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as exc:
        raise OT2Error(f"OT-2 {method} {path} returned non-JSON") from exc


def _ot2_request_multipart(host: str, path: str, *, filename: str, file_bytes: bytes,
                           timeout: float = _OT2_UPLOAD_TIMEOUT,
                           max_bytes: int = _OT2_JSON_MAX_BYTES) -> "dict[str, Any]":
    """POST a protocol file as ``multipart/form-data`` (field name ``files``)."""
    safe_name = "".join(ch for ch in filename if ch.isalnum() or ch in "._-") or "protocol.py"
    boundary = _OT2_MULTIPART_BOUNDARY
    body = b"".join([
        f"--{boundary}\r\n".encode(),
        f'Content-Disposition: form-data; name="files"; filename="{safe_name}"\r\n'.encode(),
        b"Content-Type: text/x-python\r\n\r\n",
        file_bytes, b"\r\n",
        f"--{boundary}--\r\n".encode(),
    ])
    headers = {
        "Opentrons-Version": _OT2_API_HEADER,
        "User-Agent": _ot2_user_agent(),
        "Content-Type": f"multipart/form-data; boundary={boundary}",
    }
    req = urllib.request.Request(_ot2_base_url(host) + path, data=body,
                                 headers=headers, method="POST")
    try:
        with _ot2_opener().open(req, timeout=timeout) as resp:
            raw = _ot2_read_capped(resp, max_bytes)
    except urllib.error.HTTPError as exc:
        detail = ""
        try:
            detail = exc.read(64 * 1024).decode("utf-8", "replace")
            exc.close()          # release the socket promptly (up to 1024×/sweep)
        except Exception:
            pass
        raise OT2Error(f"OT-2 protocol upload → HTTP {exc.code}"
                       + (f": {detail[:400]}" if detail else "")) from exc
    except (urllib.error.URLError, OSError) as exc:
        raise OT2Error(f"cannot reach OT-2 at {host} ({exc})") from exc
    try:
        return json.loads(raw.decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as exc:
        raise OT2Error("OT-2 protocol upload returned non-JSON") from exc


def _ot2_health(host: str) -> "dict[str, Any]":
    """Robot identity + software versions (also the reachability probe)."""
    return _ot2_request_json(host, "/health")


def _ot2_instruments(host: str) -> "dict[str, Any]":
    """Attached pipettes / instruments as the robot reports them."""
    return _ot2_request_json(host, "/instruments")


def _ot2_set_lights(host: str, on: bool) -> "dict[str, Any]":
    """Toggle the status light — a zero-motion control check."""
    return _ot2_request_json(host, "/robot/lights", method="POST", payload={"on": bool(on)})


def _ot2_home(host: str) -> "dict[str, Any]":
    """Home the gantry — all axes to their reference position. Safe motion: no
    plunger actuation, no descent into labware."""
    return _ot2_request_json(host, "/robot/home", method="POST",
                             payload={"target": "robot"})


# The six OT-2 gantry axes: x/y (deck), z/a (left/right mount carriage), b/c
# (left/right plunger). Disengaging all six releases every motor's holding torque.
_OT2_ALL_AXES = ("x", "y", "z", "a", "b", "c")


def _ot2_disengage(host: str, *, axes: "list[str] | None" = None) -> "dict[str, Any]":
    """De-energise the gantry motors (``POST /motors/disengage``) so the carriage can
    be pushed by hand and no motor holds torque — a zero-descent power-down. Homes
    nothing and never actuates a plunger; safe when the robot is idle. Defaults to
    all six axes. NB: disengaging DURING a run would drop the moving gantry, so the
    caller must refuse while a run is active (the endpoint + UI both enforce that)."""
    ax = [str(a).strip().lower() for a in (axes or _OT2_ALL_AXES) if str(a).strip()]
    if not ax:
        ax = list(_OT2_ALL_AXES)
    return _ot2_request_json(host, "/motors/disengage", method="POST",
                             payload={"axes": ax})


# ── Robot discovery (network sweep + mDNS + USB link-local) ──────────────────────
# Find OT-2 robots reachable from this machine WITHOUT any new dependency. Three
# best-effort, pure-stdlib sources are merged: a concurrent ``/health`` sweep of the
# local subnet(s), a one-shot mDNS multicast (catches robots the flat sweep misses),
# and the USB link-local interface (169.254/16) an OT-2 USB-ethernet gadget brings
# up. Every candidate is CONFIRMED by a real ``GET /health`` carrying an Opentrons
# signature, so only genuine robots are ever returned.

def _ot2_probe_robot(host: str, *, timeout: float = _OT2_DISCOVER_PROBE_TIMEOUT,
                     source: str = "network") -> "dict[str, Any] | None":
    """Probe one host's ``/health``. Returns a robot-summary dict for a genuine OT-2
    (health JSON carrying an Opentrons version signature), else ``None``. NEVER
    raises — a dead or non-robot host is just a miss, so the sweep can fan out over
    a whole subnet safely."""
    try:
        h = _ot2_request_json(host, "/health", timeout=timeout,
                              max_bytes=_OT2_DISCOVER_HEALTH_MAX_BYTES)
    except Exception:
        return None
    if not isinstance(h, dict):
        return None
    # A random service answering on :31950 that ISN'T a robot-server won't carry
    # these Opentrons-only version fields; require at least one so discovery can
    # never list a false robot.
    sig = h.get("api_version") or h.get("fw_version") or h.get("system_version")
    if not sig:
        return None
    name = h.get("name")
    return {
        "host": host,
        "name": str(name) if name else host,
        "model": h.get("robot_model") or "OT-2",
        "fw_version": h.get("fw_version"),
        "api_version": h.get("api_version"),
        "system_version": h.get("system_version"),
        "serial": h.get("serial_number") or h.get("serialNumber"),
        "source": source,
    }


def _ot2_host_ipv4s() -> "list[tuple[str, bool]]":
    """This machine's own IPv4 address(es) as ``(addr, is_link_local)`` pairs,
    stdlib-only (no ``netifaces``/``ifaddr``). The UDP-connect trick reads the
    source address the OS would use to reach a target WITHOUT sending a packet, one
    per probed route — a normal LAN address plus, if a USB-ethernet gadget is up, a
    ``169.254`` link-local one. ``getaddrinfo(hostname)`` supplements. Deduped;
    loopback dropped."""
    found: "dict[str, bool]" = {}
    for tgt in ("10.255.255.255", "172.31.255.255", "192.168.255.255",
                "169.254.255.255"):
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            try:
                s.connect((tgt, 9))          # no traffic — only fixes the local addr
                ip = s.getsockname()[0]
            finally:
                s.close()
        except OSError:
            continue
        if ip and not ip.startswith("127."):
            found[ip] = ip.startswith("169.254.")
    try:
        for info in socket.getaddrinfo(socket.gethostname(), None, socket.AF_INET):
            ip = str(info[4][0])          # AF_INET sockaddr is (addr, port)
            if ip and not ip.startswith("127."):
                found.setdefault(ip, ip.startswith("169.254."))
    except OSError:
        pass
    return list(found.items())


def _ot2_sweep_candidates(*, max_hosts: int = _OT2_DISCOVER_MAX_HOSTS
                          ) -> "list[tuple[str, str]]":
    """Every host address in the local ``/24``(s) to probe, as ``(host, source)``
    pairs. Each of this machine's interfaces contributes its ``/24`` (254 hosts) —
    enough for a home/lab LAN, small enough to sweep in a few seconds; a wider
    prefix is clamped to ``/24`` so a ``/16`` can't explode into 65k probes. A
    link-local (``169.254``) interface is tagged ``"usb"`` (that's the address
    family an OT-2 USB gadget uses); everything else ``"network"``. Deduped, own
    address skipped, hard-capped."""
    import ipaddress
    out: "list[tuple[str, str]]" = []
    seen: "set[str]" = set()
    for ip, is_ll in _ot2_host_ipv4s():
        try:
            net = ipaddress.ip_network(f"{ip}/24", strict=False)
        except ValueError:
            continue
        src = "usb" if is_ll else "network"
        for host in net.hosts():
            h = str(host)
            if h == ip or h in seen:
                continue
            seen.add(h)
            out.append((h, src))
            if len(out) >= max_hosts:
                return out
    return out


def _ot2_build_mdns_query(service: str = "_http._tcp.local") -> bytes:
    """A minimal one-question mDNS packet: a PTR query for ``service``. Header (id 0,
    flags 0, 1 question, 0 answers) + length-prefixed QNAME + QTYPE 12 (PTR) +
    QCLASS ``0x8001`` — IN (1) with the top **QU (unicast-response) bit** set, so
    responders reply directly to our ephemeral socket instead of only to the
    multicast group (which, un-joined + off port 5353, we'd never receive)."""
    import struct
    header = struct.pack(">HHHHHH", 0, 0, 1, 0, 0, 0)
    qname = b"".join(bytes([len(lbl)]) + lbl.encode("ascii")
                     for lbl in service.split(".") if lbl) + b"\x00"
    return header + qname + struct.pack(">HH", 12, 0x8001)


def _ot2_mdns_responders(*, timeout: float = _OT2_MDNS_TIMEOUT) -> "list[str]":
    """Best-effort mDNS: multicast an ``_http._tcp`` PTR query and collect the
    SOURCE address of every responder. We deliberately do NOT parse the DNS answer
    (name compression is a footgun) — we harvest who replied and let ``/health``
    decide which are OT-2s. Pure stdlib, swallow-all (mDNS is a bonus on top of the
    subnet sweep). Returns responder IPs (loopback dropped)."""
    ips: "set[str]" = set()
    sock = None
    try:
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        try:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        except OSError:
            pass
        sock.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, 2)
        sock.settimeout(0.5)
        sock.sendto(_ot2_build_mdns_query(), (_OT2_MDNS_GROUP, _OT2_MDNS_PORT))
        deadline = _util._monotonic() + max(0.1, float(timeout))
        while _util._monotonic() < deadline:
            try:
                _data, addr = sock.recvfrom(2048)
            except socket.timeout:
                continue
            except OSError:
                break
            if addr and addr[0] and not addr[0].startswith("127."):
                ips.add(addr[0])
    except OSError:
        pass
    finally:
        if sock is not None:
            try:
                sock.close()
            except OSError:
                pass
    return sorted(ips)


def _ot2_usb_serial_hint() -> "str | None":
    """Best-effort, non-fatal: is an OT-2 cabled over USB as a serial device? Linux
    only (reads ``/dev/serial/by-id``). Returns a short human note or ``None``. It
    does NOT imply serial control — the OT-2 API is HTTP; this only flags a cabled
    robot that may need USB networking (or Wi-Fi/Ethernet) to expose an address we
    can actually reach."""
    import glob
    import os
    try:
        for path in sorted(glob.glob("/dev/serial/by-id/*")):
            low = os.path.basename(path).lower()
            if "opentrons" in low or "ot-2" in low or "ot2" in low:
                return (f"USB serial device present ({os.path.basename(path)}) — "
                        "enable USB networking (or Wi-Fi/Ethernet) so it gets an "
                        "address SpliceCraft can reach.")
    except OSError:
        pass
    return None


def _ot2_is_lan_ip(ip: str) -> bool:
    """True iff ``ip`` is a private (RFC-1918) or link-local IPv4 address — the only
    ranges an OT-2 legitimately lives on. Discovery filters every AUTO-discovered
    candidate (subnet sweep + mDNS responder) through this so a spoofed mDNS reply
    can't point the probe at an arbitrary off-LAN address. The user's own persisted
    host is exempt (it's config, not discovery, and may be a ``.local`` name)."""
    import ipaddress
    try:
        addr = ipaddress.ip_address(ip)
    except ValueError:
        return False
    return bool(addr.is_private or addr.is_link_local)


def _ot2_discover(*, probe_timeout: float = _OT2_DISCOVER_PROBE_TIMEOUT,
                  max_hosts: int = _OT2_DISCOVER_MAX_HOSTS,
                  workers: int = _OT2_DISCOVER_WORKERS, use_mdns: bool = True,
                  extra_hosts: "list[str] | None" = None,
                  on_progress: "Any" = None, cancel: "Any" = None
                  ) -> "list[dict[str, Any]]":
    """Discover OT-2 robots reachable from this machine. Merges three best-effort,
    dependency-free candidate sources — mDNS responders, the local ``/24`` subnet
    sweep (incl. any USB link-local interface), and caller ``extra_hosts`` (e.g. the
    persisted host) — then CONFIRMS each with a real ``GET /health`` carrying an
    Opentrons signature, concurrently over a bounded pool. Returns a list of robot
    dicts, sorted known/network-first then by name. ``on_progress(done, total)``
    fires as probes finish; a truthy ``cancel()`` aborts the sweep early. Egress is
    gated by the fail-closed demo hook, so a web-demo session never scans a LAN."""
    import concurrent.futures
    _state._demo_block_network_hook("OT-2 discovery")

    cands: "list[tuple[str, str]]" = []      # (host, source), deduped
    seen: "set[str]" = set()

    def _add(host: str, source: str) -> None:
        host = (host or "").strip()
        if not host or host in seen:
            return
        # An auto-discovered candidate (sweep / mDNS responder) MUST be a LAN IP so a
        # forged mDNS source address can't redirect the probe off-box; the user's own
        # persisted host ("known") is exempt.
        if source != "known" and not _ot2_is_lan_ip(host):
            return
        seen.add(host)
        cands.append((host, source))

    # Known/persisted hosts first (so a reachable one sorts to the top), then mDNS,
    # then the flat sweep.
    for h in (extra_hosts or []):
        _add(str(h), "known")
    if use_mdns:
        try:
            for ip in _ot2_mdns_responders():
                _add(ip, "usb" if ip.startswith("169.254.") else "network")
        except Exception:
            _log.debug("[ot2] mDNS discovery failed (ignored)", exc_info=True)
    try:
        for host, source in _ot2_sweep_candidates(max_hosts=max_hosts):
            _add(host, source)
    except Exception:
        _log.debug("[ot2] subnet enumeration failed (ignored)", exc_info=True)

    total = len(cands)
    found: "list[dict[str, Any]]" = []
    if not total:
        return found

    def _cancelled() -> bool:
        try:
            return bool(cancel and cancel())
        except Exception:
            return False

    done = 0
    deadline = _util._monotonic() + _OT2_DISCOVER_OVERALL_TIMEOUT
    ex = concurrent.futures.ThreadPoolExecutor(max_workers=max(1, int(workers)))
    try:
        futs = {ex.submit(_ot2_probe_robot, host, timeout=probe_timeout, source=src): host
                for host, src in cands}
        # Drain in short waits (not a bare `as_completed`) so `cancel()` and the
        # overall deadline stay responsive even while every in-flight probe is
        # blocked on a slow host — a bare as_completed would only re-check between
        # completions and a wedged batch would freeze the whole loop.
        pending = set(futs)
        while pending and not _cancelled() and _util._monotonic() < deadline:
            batch, pending = concurrent.futures.wait(
                pending, timeout=0.25,
                return_when=concurrent.futures.FIRST_COMPLETED)
            for fut in batch:
                done += 1
                if on_progress:
                    try:
                        on_progress(done, total)
                    except Exception:
                        pass
                try:
                    robot = fut.result()
                except Exception:
                    robot = None
                if robot:
                    found.append(robot)
    finally:
        # Don't block on the remaining probes when cancelled or erroring — a sweep
        # of dead hosts would otherwise hang on their timeouts.
        ex.shutdown(wait=False, cancel_futures=True)

    rank = {"known": 0, "network": 1, "usb": 2}
    found.sort(key=lambda r: (rank.get(str(r.get("source")), 3),
                              str(r.get("name") or "").lower()))
    return found


# ── Telemetry / fault detection ─────────────────────────────────────────────────
# The OT-2 REST API does not stream live gantry XYZ, so "position" is tracked at
# the PROTOCOL level: which command the robot is currently executing and — on a
# crash — which command FAILED and why. Combined with the peripheral state
# (instrument OK flags + specs, motor engagement, deck-calibration health, module
# sensor readings, robot settings/variables, reachability) this gives direct,
# immediate crash detection: a pipette that hits labware fails its move / pick-up
# command and the run goes to ``failed`` with an error the monitor surfaces the
# instant it polls. All of it is read-only — nothing here actuates.

def _ot2_try_json(host: str, path: str, *,
                  timeout: float = _OT2_SHORT_TIMEOUT) -> "tuple[dict[str, Any] | None, str | None]":
    """Best-effort GET: ``(data, None)`` or ``(None, error)`` — never raises, so
    one dead endpoint can't sink a whole state snapshot."""
    try:
        return _ot2_request_json(host, path, timeout=timeout), None
    except OT2Error as exc:
        return None, str(exc)


def _ot2_run_commands(host: str, rid: str, *, page_length: int = 200
                      ) -> "tuple[list[dict[str, Any]], int]":
    """Per-command status for a run (page-capped) PLUS the run's true command total.

    Returns ``(commands, total)``. Each command carries id / type / status / error
    detail (or ``None``) and a brief position hint (labware / well / mount) — the
    closest thing to 'where the pipette is'. ``total`` is the run's full command
    count so far (``meta.totalLength``): the OT-2 creates run commands as the
    protocol executes, so it climbs from 0 toward the analysed length and makes a
    real progress numerator (the page itself holds at most ``page_length`` items)."""
    data, _ = _ot2_try_json(host, f"/runs/{rid}/commands?pageLength={int(page_length)}")
    if not isinstance(data, dict):   # a non-dict (e.g. a top-level list) must not raise
        data = {}
    out: "list[dict[str, Any]]" = []
    for c in (data.get("data") or []):
        params = c.get("params") or {}
        at = {k: params[k] for k in ("labwareId", "wellName", "pipetteId", "mount")
              if k in params}
        err = c.get("error")
        out.append({
            "id": c.get("id"),
            "commandType": c.get("commandType"),
            "status": c.get("status"),
            "error": (err.get("detail") if isinstance(err, dict) else err),
            "at": at or None,
        })
    meta = data.get("meta")
    total = meta.get("totalLength") if isinstance(meta, dict) else None
    # Reject bool (True is an int), non-int, too-small, or absurdly-large values —
    # any of which would poison the progress numerator — falling back to the page.
    if (not isinstance(total, int) or isinstance(total, bool)
            or total < len(out) or total > 1_000_000):
        total = len(out)
    return out, total


def _ot2_active_run(host: str) -> "str | None":
    """The id of the current run, if any — so the monitor can watch a run started
    from the Opentrons App too, not just one SpliceCraft launched."""
    data, _ = _ot2_try_json(host, "/runs")
    runs = (data or {}).get("data") or []
    for r in runs:
        if r.get("current"):
            return r.get("id")
    for r in reversed(runs):
        if r.get("status") in ("running", "idle", "paused", "finishing",
                               "blocked-by-open-door"):
            return r.get("id")
    return None


def _ot2_run_state(host: str, rid: str) -> "dict[str, Any]":
    """Live run state: overall status, the current command (position), any failed
    commands, run-level errors, and a progress numerator (``command_total`` = run
    commands created so far vs ``command_count`` = the page actually seen)."""
    data, err = _ot2_try_json(host, f"/runs/{rid}")
    if data is None:
        return {"id": rid, "status": "unknown", "error": err, "errors": [],
                "current_command": None, "failed_commands": [], "command_count": 0,
                "command_total": 0}
    d = data.get("data") or {}
    cmds, total = _ot2_run_commands(host, rid)
    running = [c for c in cmds if c["status"] == "running"]
    progressed = [c for c in cmds if c["status"] in ("running", "succeeded")]
    current = running[0] if running else (progressed[-1] if progressed else None)
    return {
        "id": rid,
        "status": d.get("status"),
        "errors": d.get("errors") or [],
        "current_command": current,
        "failed_commands": [c for c in cmds if c["status"] == "failed"],
        "command_count": len(cmds),
        "command_total": total,
    }


def _ot2_detect_faults(state: "dict[str, Any]") -> "list[str]":
    """Pure crash/fault extraction from a state snapshot — the direct signals a
    physical mishap produces: unreachable robot, a pipette subsystem reporting
    not-ok, bad/singular deck calibration, a module in an error state, a failed
    run, or a specific failed command (with where it failed)."""
    faults: "list[str]" = list(state.get("faults", []))
    if not state.get("reachable", False):
        return faults or ["unreachable"]
    for inst in state.get("instruments", []):
        if inst.get("ok") is False:
            faults.append(f"instrument fault: {inst.get('mount')} "
                          f"{inst.get('model')} reports not-ok")
    cal = state.get("calibration") or {}
    if cal.get("deck_status") not in (None, "OK", "IDENTITY"):
        faults.append(f"deck calibration status {cal.get('deck_status')}")
    if cal.get("marked_bad"):
        faults.append("deck calibration marked bad")
    for m in state.get("modules", []):
        s = str(m.get("status") or "").lower()
        if "error" in s or "fault" in s:
            faults.append(f"module {m.get('id') or m.get('moduleType')}: {m.get('status')}")
    run = state.get("run") or {}
    if run.get("status") == "failed":
        errs = run.get("errors") or []
        detail = errs[0].get("detail") if errs and isinstance(errs[0], dict) else None
        faults.append(f"run failed: {detail or 'see run errors'}")
    for c in run.get("failed_commands", []):
        where = f" at {c.get('at')}" if c.get("at") else ""
        faults.append(f"command failed: {c.get('commandType')} — {c.get('error')}{where}")
    return faults


def _ot2_pipette_offsets(host: str) -> "list[dict[str, Any]]":
    """Per-pipette offset calibrations the robot has stored (empty when a pipette
    was never calibrated in the Opentrons App). Best-effort — never raises."""
    data, _ = _ot2_try_json(host, "/calibration/pipette_offset")
    return (data or {}).get("data") or []


def _ot2_tip_lengths(host: str) -> "list[dict[str, Any]]":
    """Stored tip-length calibrations (per pipette + tip-rack pairing)."""
    data, _ = _ot2_try_json(host, "/calibration/tip_length")
    return (data or {}).get("data") or []


def _ot2_state(host: str, *, run_id: "str | None" = None) -> "dict[str, Any]":
    """A full snapshot of everything the robot exposes — reachability + versions,
    pipette OK flags + volume specs, motor engagement, deck / instrument
    calibration health, attached modules and their sensor readings, the status
    light, robot settings (variables), and (for the active run or ``run_id``) live
    run / command state — with detected ``faults`` and an ``ok`` verdict."""
    state: "dict[str, Any]" = {"host": host, "reachable": False, "faults": []}
    health, herr = _ot2_try_json(host, "/health")
    if health is None:
        state["faults"].append(f"unreachable: {herr}")
        state["ok"] = False
        return state
    state["reachable"] = True
    state["health"] = {k: health.get(k) for k in
                       ("name", "api_version", "fw_version", "system_version", "robot_model")}

    instruments: "list[dict[str, Any]]" = []
    idata, _ = _ot2_try_json(host, "/instruments")
    for it in ((idata or {}).get("data") or []):
        d = it.get("data") or {}
        instruments.append({
            "mount": it.get("mount"),
            "model": it.get("instrumentModel") or it.get("instrumentName"),
            "ok": it.get("ok"),
            "min_volume": d.get("min_volume"),
            "max_volume": d.get("max_volume"),
            "channels": d.get("channels"),
        })
    state["instruments"] = instruments

    motors, _ = _ot2_try_json(host, "/motors/engaged")
    state["motors"] = motors or {}

    lights, _ = _ot2_try_json(host, "/robot/lights")
    state["lights"] = (lights or {}).get("on")

    # Door switch: status + whether the robot enforces closed-for-protocol. Present
    # even mid-run (an open door pauses the run) — cheap and safety-relevant, so it
    # is NOT skipped by the lean-poll path. Best-effort: a robot/firmware without the
    # endpoint just leaves ``door`` absent.
    door, _ = _ot2_try_json(host, "/robot/door/status")
    if isinstance(door, dict):
        inner = door.get("data")
        dd = inner if isinstance(inner, dict) else door
        state["door"] = {
            "status": dd.get("status"),
            "required_closed": bool(dd.get("doorRequiredClosedForProtocol")),
        }

    calibration: "dict[str, Any]" = {}
    cal, _ = _ot2_try_json(host, "/calibration/status")
    if cal:
        dc = cal.get("deckCalibration") or {}
        marked = ((dc.get("data") or {}).get("status") or {}).get("markedBad")
        calibration = {"deck_status": dc.get("status"), "marked_bad": bool(marked)}
    # Pipette-offset + tip-length calibration presence. Kept OUT of the hard-fault
    # list (below) so a position check can still move an uncalibrated gantry — the
    # liquid-run gate is what insists on a calibrated pipette. Skipped while polling
    # a live run (run_id set): calibration can't change mid-run, so the monitor stays
    # lean (2 fewer HTTP calls per tick).
    if run_id is None:
        poff = _ot2_pipette_offsets(host)
        cal_mounts = {str(p.get("mount")).lower() for p in poff if p.get("mount")}
        calibration["pipette_offsets"] = [{"mount": p.get("mount"), "pipette": p.get("pipette")}
                                          for p in poff]
        calibration["tip_lengths"] = len(_ot2_tip_lengths(host))
        calibration["pipettes_calibrated"] = {
            str(i.get("mount")).lower(): (str(i.get("mount")).lower() in cal_mounts)
            for i in instruments if i.get("mount")}
    if calibration:
        state["calibration"] = calibration

    modules: "list[dict[str, Any]]" = []
    mdata, _ = _ot2_try_json(host, "/modules")
    for m in ((mdata or {}).get("data") or []):
        modules.append({"id": m.get("id"), "moduleType": m.get("moduleType"),
                        "status": m.get("status"), "data": m.get("data")})
    state["modules"] = modules

    sdata, _ = _ot2_try_json(host, "/settings")
    if sdata:
        state["settings"] = {s.get("id"): s.get("value")
                             for s in (sdata.get("settings") or []) if s.get("id")}

    if run_id is None:
        run_id = _ot2_active_run(host)
    if run_id:
        state["run"] = _ot2_run_state(host, run_id)

    state["faults"] = _ot2_detect_faults(state)
    state["ok"] = not state["faults"]
    return state


# ── Run control (in-flight actuation: pause / resume / stop) ─────────────────────
# The only three commands you can send a live run. "resume" is the same "play"
# action that starts a run — the robot picks up where a pause left off. These are
# the manual counterpart to the automatic stop-on-fault the run monitor performs.
_OT2_RUN_ACTIONS: "dict[str, str]" = {
    "pause": "pause", "resume": "play", "stop": "stop",
    "play": "play", "cancel": "stop",   # aliases
}
# A run id is interpolated into the request path — keep it to the UUID-shaped
# charset so a crafted value (slashes / '..') can't manipulate the URL path.
_OT2_RUN_ID_CHARS = frozenset(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.")


def _ot2_valid_run_id(rid: Any) -> str:
    rid = str(rid)
    if not rid or len(rid) > 128 or any(c not in _OT2_RUN_ID_CHARS for c in rid):
        raise OT2Error(f"invalid run id {rid!r}")
    return rid


def _ot2_run_action(host: str, rid: str, action: str) -> "dict[str, Any]":
    """Send a run-control action to an in-flight run — ``pause`` / ``resume`` /
    ``stop`` (``cancel`` is an alias for ``stop``). Raises ``OT2Error`` on an
    unknown action or a path-unsafe run id."""
    act = _OT2_RUN_ACTIONS.get(str(action).strip().lower())
    if act is None:
        raise OT2Error(f"unknown run action {action!r} (use pause / resume / stop)")
    rid = _ot2_valid_run_id(rid)
    return _ot2_request_json(host, f"/runs/{rid}/actions", method="POST",
                             payload={"data": {"actionType": act}})


def _ot2_stop_run(host: str, rid: str) -> "dict[str, Any]":
    """Halt a run (stop action) — used to abort the instant a fault is detected."""
    return _ot2_run_action(host, rid, "stop")


def _ot2_run_control(host: str, action: str, *,
                     run_id: "str | None" = None) -> "dict[str, Any]":
    """Resolve the current run (or use ``run_id``) and send it a control action.
    Returns ``{"ok", "action", "run_id"}``; raises ``OT2Error`` when there is no
    active run, the action is unknown, or the transport fails."""
    rid = run_id or _ot2_active_run(host)
    if not rid:
        raise OT2Error("no active run to control (nothing is running on the robot)")
    _ot2_run_action(host, rid, action)
    _log.info("[ot2] run-control %s on %s (run %s)", str(action).strip().lower(), host, rid)
    return {"ok": True, "action": str(action).strip().lower(), "run_id": rid}


def _ot2_monitor(host: str, *, run_id: "str | None" = None,
                 on_state: "Any" = None, interval: float = _OT2_POLL_INTERVAL,
                 max_seconds: float = _OT2_RUN_POLL_TIMEOUT,
                 stop: "Any" = None) -> "dict[str, Any]":
    """Poll robot state until the (active) run reaches a terminal status, a fault
    is detected, ``stop()`` returns truthy, or ``max_seconds`` elapses — calling
    ``on_state(snapshot)`` each tick. Returns the final snapshot (with a
    ``stopped_reason``). Read-only: it observes, it never actuates."""
    deadline = _util._monotonic() + max_seconds
    while True:
        snap = _ot2_state(host, run_id=run_id)
        if on_state:
            try:
                on_state(snap)
            except Exception:
                pass
        if snap.get("faults"):
            snap["stopped_reason"] = "fault"
            return snap
        run = snap.get("run") or {}
        if run.get("status") in ("succeeded", "failed", "stopped"):
            snap["stopped_reason"] = f"run-{run.get('status')}"
            return snap
        if stop and stop():
            snap["stopped_reason"] = "cancelled"
            return snap
        if _util._monotonic() > deadline:
            snap["stopped_reason"] = "timeout"
            return snap
        time.sleep(interval)


def _ot2_analyze(host: str, protocol_text: str, *,
                 filename: str = "splicecraft_protocol.py") -> "dict[str, Any]":
    """Upload a protocol and wait for the robot's server-side analysis to finish.

    Returns ``{"protocol_id", "analysis_id", "status", "result", "errors",
    "commands", "pipettes", "labware"}``. ``result`` is ``"ok"`` / ``"not-ok"``;
    on ``not-ok`` the ``errors`` list explains why (this is the pre-flight the
    physical-run gate insists on).
    """
    if len(protocol_text.encode("utf-8")) > _OT2_MAX_PROTOCOL_BYTES:
        raise OT2Error("protocol is too large to upload "
                       f"({_OT2_MAX_PROTOCOL_BYTES // (1024 * 1024)} MB max)")
    up = _ot2_request_multipart(host, "/protocols", filename=filename,
                                file_bytes=protocol_text.encode("utf-8"))
    prot = up.get("data", {})
    pid = prot.get("id")
    summaries = prot.get("analysisSummaries") or []
    aid = summaries[-1]["id"] if summaries else None
    if not pid or not aid:
        raise OT2Error("OT-2 upload did not return a protocol/analysis id")

    deadline = _util._monotonic() + _OT2_ANALYSIS_POLL_TIMEOUT
    adata: "dict[str, Any]" = {}
    while True:
        adata = _ot2_request_json(host, f"/protocols/{pid}/analyses/{aid}").get("data", {})
        if adata.get("status") == "completed":
            break
        if _util._monotonic() > deadline:
            raise OT2Error("OT-2 analysis did not complete within "
                           f"{_OT2_ANALYSIS_POLL_TIMEOUT}s")
        time.sleep(_OT2_POLL_INTERVAL)

    return {
        "protocol_id": pid,
        "analysis_id": aid,
        "status": adata.get("status"),
        "result": adata.get("result"),
        "errors": adata.get("errors", []),
        "commands": adata.get("commands", []),
        "pipettes": adata.get("pipettes", []),
        "labware": adata.get("labware", []),
    }


def _ot2_labware_offsets(analysis_labware: "list[dict[str, Any]]",
                         plan: "dict[str, Any]") -> "list[dict[str, Any]]":
    """Build a run's ``labwareOffsets`` from a plan's per-labware ``offset`` (x/y/z
    mm) matched to each analysed labware's ``definitionUri`` + slot. An offset with
    no matching analysed labware — or a zero vector (a no-op) — is skipped."""
    p = _ot2_normalize_plan(plan)
    slot_off: "dict[str, dict[str, Any]]" = {}
    for lw in p["labware"].values():
        off = lw.get("offset")
        if isinstance(off, dict) and lw.get("slot") is not None:
            slot_off[str(lw["slot"])] = off

    def _f(off: "dict[str, Any]", k: str) -> float:
        v = off.get(k)
        return (float(v) if isinstance(v, (int, float)) and not isinstance(v, bool)
                and math.isfinite(v) else 0.0)

    out: "list[dict[str, Any]]" = []
    for al in (analysis_labware or []):
        if not isinstance(al, dict):
            continue        # untrusted robot-analysis JSON: a non-dict item is skipped
        loc = al.get("location")
        if not isinstance(loc, dict):
            loc = {}
        slot = str(loc.get("slotName") or loc.get("slot") or "")
        uri = al.get("definitionUri")
        off = slot_off.get(slot)
        if off is None or not uri:
            continue
        vec = {"x": _f(off, "x"), "y": _f(off, "y"), "z": _f(off, "z")}
        # HARD safety cap: offsets are applied at run creation AFTER analysis, so the
        # robot's own gate never sees them. A large / negative value would drive the
        # pipette below the well top (a descent) — reject rather than move unsafely.
        for axis, val in vec.items():
            if abs(val) > _OT2_MAX_OFFSET_MM:
                raise OT2Error(
                    f"labware offset {axis}={val:g} mm on slot {slot} exceeds the "
                    f"±{_OT2_MAX_OFFSET_MM:g} mm safety limit — an offset is a fine "
                    "correction; re-seat the labware if it is off by more than that")
        if vec == {"x": 0.0, "y": 0.0, "z": 0.0}:
            continue
        out.append({"definitionUri": uri, "location": {"slotName": slot}, "vector": vec})
    return out


def _ot2_compile_position_check(plan: "dict[str, Any]", *,
                                wells: "list[str] | None" = None) -> str:
    """Compile a MOTION-ONLY 'position check' protocol: load the deck, home, move
    the pipette to the TOP of each loaded labware's reference wells (A1 + the far
    corner by default), then home again. Emits NO aspirate / dispense / tip
    pick-up — the plunger is never actuated and the pipette never descends into a
    well — so it is safe on occupied slots (alignment / offset verification).
    Raises ``OT2Error`` when the deck has no loadable labware."""
    p = _ot2_normalize_plan(plan)
    loaded: "list[tuple[str, dict[str, Any]]]" = []
    for t in p["tips"]:
        if t.get("labware") and t.get("slot") is not None:
            loaded.append((f"tips_{t['slot']}", t))
    for lid, lw in p["labware"].items():
        if lw.get("slot") is not None:
            loaded.append((lid, lw))
    if not loaded:
        raise OT2Error("position check: the deck has no labware to move to")
    # Validate slots BEFORE interpolating them into the emitted Python (the main
    # compiler does this via _ot2_validate_plan; this path must not skip it, or a
    # non-int slot could emit malformed code / a non-.top() move).
    for lid, lw in loaded:
        slot = lw.get("slot")
        if not (isinstance(slot, int) and not isinstance(slot, bool) and 1 <= slot <= 11):
            raise OT2Error(f"position check: labware {lid!r} has an invalid slot "
                           f"{slot!r} (slots are integers 1-11)")
        # Same reason for the `repr()`-embedded custom definition below.
        if isinstance(lw.get("definition"), dict):
            lit = _ot2_literal_errors(lw["definition"],
                                      f"position check: labware {lid!r} definition")
            if lit:
                raise OT2Error("\n  - ".join(["invalid custom labware:"] + lit))

    md = p["metadata"]
    out: "list[str]" = [
        "from opentrons import protocol_api",
        "",
        "metadata = {",
        f"    \"protocolName\": {json.dumps(md['name'] + ' — position check')},",
        f"    \"author\": {json.dumps(md['author'])},",
        "    \"description\": \"SpliceCraft position check — gantry moves to well "
        "tops only (no liquid handling).\",",
        f"    \"apiLevel\": {json.dumps(p['api_level'])},",
        "}",
        "",
        "def run(protocol: protocol_api.ProtocolContext):",
    ]
    id_var: "dict[str, str]" = {}
    used: "set[str]" = set()
    for lid, lw in loaded:
        var = _ot2_safe_var(lid, "lw_")
        while var in used:
            var += "_"
        used.add(var)
        id_var[lid] = var
        if isinstance(lw.get("definition"), dict):
            out.append(f"    {var} = protocol.load_labware_from_definition("
                       f"{lw['definition']!r}, {lw['slot']})")
        else:
            out.append(f"    {var} = protocol.load_labware("
                       f"{json.dumps(_ot2_resolve_labware(str(lw.get('labware', ''))))}, "
                       f"{lw['slot']})")
    out.append(f"    pipette = protocol.load_instrument({json.dumps(p['pipette'])}, "
               f"{json.dumps(p['mount'])})")
    out.append("    protocol.home()")
    for lid, lw in loaded:
        ew = _ot2_entry_wells(lw)
        targets = wells or ([ew[0], ew[-1]] if ew else ["A1"])
        seen_set: "set[str]" = set()
        seen: "list[str]" = []
        for w in targets[:_OT2_MAX_POSCHECK_WELLS]:   # cap + O(n) dedup (was O(n²))
            if w not in seen_set:
                seen_set.add(w)
                seen.append(w)
        for w in seen:
            out.append(f"    pipette.move_to({id_var[lid]}[{json.dumps(w)}].top())")
    out.append("    protocol.home()")
    return "\n".join(out) + "\n"


def _ot2_pipette_base(name: "Any") -> str:
    """A pipette identifier without its version/generation suffix so an analysed
    ``pipetteName`` (``p300_single``) matches an attached ``instrumentModel``
    (``p300_single_v1.5`` / ``p300_single_gen2``). Lower-cased; empty on junk."""
    s = str(name or "").strip().lower()
    for sep in ("_v", "_gen"):
        i = s.find(sep)
        if i != -1 and i + len(sep) < len(s) and s[i + len(sep)].isdigit():
            s = s[:i]
    return s


def _ot2_pipette_mismatch(analysis_pipettes: "list[dict[str, Any]]",
                          instruments: "list[dict[str, Any]]") -> "list[str]":
    """Which pipettes the analysed protocol loads are NOT satisfied by an attached
    instrument on the same mount + (version-insensitive) model. Returns
    human-readable mismatch strings; empty when every required pipette is present.
    The robot's *analysis* simulates with the requested pipette regardless of what
    is physically attached, so a wrong / absent pipette passes analysis and only
    fails at ``load_instrument`` once the run starts — this catches it BEFORE any
    motion so the gate can refuse."""
    attached: "dict[str, str]" = {}
    for it in (instruments or []):
        if not isinstance(it, dict):
            continue
        mount = str(it.get("mount") or "").strip().lower()
        if mount:
            attached[mount] = _ot2_pipette_base(it.get("model"))
    problems: "list[str]" = []
    for p in (analysis_pipettes or []):
        if not isinstance(p, dict):
            continue
        want_mount = str(p.get("mount") or "").strip().lower()
        want = _ot2_pipette_base(p.get("pipetteName") or p.get("pipetteModel")
                                 or p.get("model"))
        if not want:
            continue
        have = attached.get(want_mount)
        if have is None:
            problems.append(f"protocol needs {want} on the {want_mount or '?'} mount "
                            "but nothing is attached there")
        elif have != want:
            problems.append(f"protocol needs {want} on the {want_mount} mount "
                            f"but a {have} is attached")
    return problems


def _ot2_run_protocol(host: str, protocol_text: str, *, confirm: bool = False,
                      filename: str = "splicecraft_protocol.py", poll: bool = True,
                      on_state: "Any" = None, stop_on_fault: bool = True,
                      offset_plan: "dict[str, Any] | None" = None,
                      require_pipette_cal: bool = True,
                      indicator_lights: bool = True,
                      on_analysis: "Any" = None) -> "dict[str, Any]":
    """Analyse, then (only if it passed AND the caller confirmed) run a protocol
    on real hardware, monitoring state throughout for a crash.

    The gate is deliberate and non-negotiable: this refuses to move the gantry
    unless ``_ot2_analyze`` returns ``result == "ok"`` *and* ``confirm=True`` was
    passed. It also refuses to start on an already-faulted robot (pre-flight
    ``_ot2_state`` check — bad calibration, a pipette subsystem error, etc.), when
    the door-safety switch is enabled but the door is open (``reason:
    "door-open"``), and when the attached pipette does not match what the protocol
    loads (``reason: "pipette-mismatch"`` — analysis can't see the hardware, so a
    wrong pipette would otherwise only fail mid-run).

    While the run executes it polls a full ``_ot2_state`` snapshot each tick,
    calls ``on_state(snapshot)`` (for a live UI / log), turns the rail lights on as
    a "robot is moving" indicator (restored afterwards), and — the instant a fault
    is detected (a failed command, an instrument fault, …) — records it and, if
    ``stop_on_fault``, halts the run. If the run overruns ``_OT2_RUN_POLL_TIMEOUT``
    it is STOPPED on the robot before the timeout is raised (never left moving
    unattended). The result carries ``crashed`` plus the ``faults`` /
    ``failed_commands`` that explain what went wrong and where. ``indicator_lights``
    (default on) governs the rail-light signalling.
    """
    analysis = _ot2_analyze(host, protocol_text, filename=filename)
    if analysis["result"] != "ok":
        return {"ran": False, "reason": "analysis-failed", **analysis}
    # Hand the analysed command list to the caller (for a determinate progress bar)
    # the moment it's known — before the confirm gate, so a dry run reports it too.
    if on_analysis:
        try:
            on_analysis(analysis)
        except Exception:
            _log.debug("[ot2] on_analysis callback raised (ignored)", exc_info=True)
    if not confirm:
        return {"ran": False, "reason": "confirm-required",
                "detail": "physical run needs confirm=True (analysis passed)",
                **analysis}

    pre = _ot2_state(host)
    if pre.get("faults"):
        return {"ran": False, "reason": "robot-unhealthy",
                "faults": pre["faults"], "state": pre, **analysis}
    # Door interlock: only when the robot itself enforces closed-for-protocol (its
    # door-safety switch is on). If enforced AND open, the robot would refuse/pause
    # anyway — surface it as a clear reason rather than a mid-run stall.
    door = pre.get("door") or {}
    if door.get("required_closed") and str(door.get("status") or "").lower() == "open":
        return {"ran": False, "reason": "door-open",
                "detail": "the robot's door-safety switch is enabled and the door is "
                          "open — close the door before running",
                "state": pre, **analysis}
    if require_pipette_cal:
        calib = (pre.get("calibration") or {}).get("pipettes_calibrated") or {}
        if not calib:
            # Fail CLOSED: an empty map means we couldn't confirm a mounted pipette's
            # calibration (a transient /instruments read failure would otherwise let
            # an uncalibrated liquid run slip through).
            return {"ran": False, "reason": "calibration-unknown",
                    "detail": "could not confirm pipette calibration (robot busy, or the "
                              "instruments read failed) — retry, or check the robot",
                    "state": pre, **analysis}
        uncal = sorted(m for m, ok in calib.items() if not ok)
        if uncal:
            return {"ran": False, "reason": "pipette-not-calibrated",
                    "detail": f"pipette(s) {', '.join(uncal)} have no offset calibration "
                              "— calibrate in the Opentrons App first",
                    "state": pre, **analysis}

    # Attached-pipette-vs-protocol match: calibration confirms a pipette is
    # calibrated, NOT that it is the RIGHT one. Analysis simulates with whatever the
    # protocol requests, so a wrong / absent pipette only fails at load_instrument
    # once the run starts — refuse here, before any motion.
    mismatch = _ot2_pipette_mismatch(analysis.get("pipettes") or [],
                                     pre.get("instruments") or [])
    if mismatch:
        return {"ran": False, "reason": "pipette-mismatch",
                "detail": "; ".join(mismatch), "state": pre, **analysis}

    run_data: "dict[str, Any]" = {"protocolId": analysis["protocol_id"]}
    if offset_plan is not None:
        offsets = _ot2_labware_offsets(analysis.get("labware") or [], offset_plan)
        if offsets:
            run_data["labwareOffsets"] = offsets
            _log.info("[ot2] applying %d labware offset(s)", len(offsets))
    created = _ot2_request_json(host, "/runs", method="POST", payload={"data": run_data})
    rid = created.get("data", {}).get("id")
    if not rid:
        raise OT2Error("OT-2 did not return a run id")
    _log.info("[ot2] starting physical run %s on %s", rid, host)
    _ot2_request_json(host, f"/runs/{rid}/actions", method="POST",
                      payload={"data": {"actionType": "play"}})

    if not poll:
        return {"ran": True, "run_id": rid, "run_status": "running", **analysis}

    # Rail lights ON for the duration as a physical "robot is moving" signal, then
    # restored to their prior state when we stop watching. The finally makes the
    # restore fire on a clean finish, a fault, OR the timeout. Best-effort — a
    # lights hiccup never affects the run itself.
    prev_lights = pre.get("lights")
    if indicator_lights:
        try:
            _ot2_set_lights(host, True)
        except OT2Error:
            pass
    try:
        deadline = _util._monotonic() + _OT2_RUN_POLL_TIMEOUT
        while True:
            try:
                snap = _ot2_state(host, run_id=rid)
            except Exception:
                # SAFETY: an unexpected monitor error (e.g. a malformed robot
                # response) must never leave the gantry moving unwatched — stop the
                # run, then re-raise (the finally still restores the lights).
                _log.exception("[ot2] monitor error during run %s — stopping it", rid)
                try:
                    _ot2_stop_run(host, rid)
                except OT2Error:
                    pass
                raise
            if on_state:
                try:
                    on_state(snap)
                except Exception:
                    _log.debug("[ot2] on_state callback raised (ignored)", exc_info=True)
            run = snap.get("run") or {}
            status = run.get("status", "unknown")
            faults = snap.get("faults") or []

            if faults and status != "succeeded":
                # A crash mid-flight. The robot already halts on a hard move error
                # (status 'failed'); for any other detected fault, stop it ourselves.
                if stop_on_fault and status not in ("failed", "stopped"):
                    try:
                        _ot2_stop_run(host, rid)
                    except OT2Error:
                        pass
                _log.warning("[ot2] fault during run %s: %s", rid, faults)
                return {"ran": True, "run_id": rid, "run_status": status, "crashed": True,
                        "faults": faults, "failed_commands": run.get("failed_commands") or [],
                        "run_errors": run.get("errors") or [], "state": snap, **analysis}

            if status in ("succeeded", "failed", "stopped"):
                _log.info("[ot2] run %s finished: %s", rid, status)
                return {"ran": True, "run_id": rid, "run_status": status, "crashed": False,
                        "faults": [], "failed_commands": run.get("failed_commands") or [],
                        "run_errors": run.get("errors") or [], "state": snap, **analysis}

            if _util._monotonic() > deadline:
                # SAFETY: never leave the gantry moving unattended. Stop the run on
                # the robot BEFORE surfacing the timeout (this used to just stop
                # watching while the robot kept running).
                _log.warning("[ot2] run %s exceeded %ss — stopping it",
                             rid, _OT2_RUN_POLL_TIMEOUT)
                try:
                    _ot2_stop_run(host, rid)
                except OT2Error:
                    pass
                raise OT2Error(f"OT-2 run {rid} did not finish within "
                               f"{_OT2_RUN_POLL_TIMEOUT}s — the run was stopped")
            time.sleep(_OT2_POLL_INTERVAL)
    finally:
        if indicator_lights:
            try:
                _ot2_set_lights(host, bool(prev_lights))
            except OT2Error:
                pass


def _ot2_run_position_check(host: str, plan: "dict[str, Any]", *,
                            wells: "list[str] | None" = None, confirm: bool = False,
                            poll: bool = True, on_state: "Any" = None,
                            stop_on_fault: bool = True) -> "dict[str, Any]":
    """Compile + run a MOTION-ONLY position check (a move-to-top tour). Gated like a
    run (needs ``confirm=True``) but RELAXES the pipette-calibration requirement —
    you position-check to VERIFY alignment, often before calibrating. Reachability +
    deck calibration are still required (they gate any safe motion). Applies the
    plan's labware offsets so you can see whether an offset lands the pipette dead
    centre. No plunger, no descent — safe on occupied slots."""
    proto = _ot2_compile_position_check(plan, wells=wells)
    res = _ot2_run_protocol(host, proto, confirm=confirm, poll=poll, on_state=on_state,
                            stop_on_fault=stop_on_fault, require_pipette_cal=False,
                            offset_plan=plan)
    res["position_check"] = True
    return res
