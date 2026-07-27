"""Layered-import guardrail for the flat-sibling refactor.

The monolith is being carved into top-level ``splicecraft_*.py`` sibling
modules (continuing the ``splicecraft_biology.py`` / ``splicecraft_cli.py``
precedent). ``splicecraft.py`` stays the hub: it imports the siblings and
re-exports their public surface. The one architectural invariant is:

    A sibling may import only same-or-lower layers. Never a higher layer,
    and NEVER the ``splicecraft`` hub itself (that is a guaranteed cycle).

This is the static guard. It catches three failure modes before they become
runtime ``ImportError``s (often surfacing as confusing "partially initialized
module" errors):

  1. a sibling importing the ``splicecraft`` hub (cycle: hub imports sibling),
  2. an upward import (e.g. a widget importing a screen),
  3. any import cycle among siblings.

It also fails if a ``splicecraft_*`` module is not classified into a layer --
so a new sibling cannot silently appear without an explicit layer decision.
"""
from __future__ import annotations

import ast
import tomllib
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent

# Layer assignment by module-name prefix, lowest first. Longest matching prefix
# wins. Add a rule here when you introduce a new sibling family -- an
# unclassified ``splicecraft_*`` module is a hard failure (see test below).
_LAYER_RULES = [
    ("splicecraft_errors", 0),
    ("splicecraft_logging", 0),
    ("splicecraft_constants", 0),
    ("splicecraft_paths", 0),
    ("splicecraft_persistence", 0),
    ("splicecraft_util", 0),        # pure cross-cutting helpers (natural sort, sanitise, ...)
    ("splicecraft_net", 0),         # shared SSRF-hardened network-fetch primitives (state + stdlib only)
    ("splicecraft_record", 1),      # GenBank <-> SeqRecord serialization (+ LRU parse cache)
    ("splicecraft_dataaccess", 1),  # _typed_clone + (Phase D) the domain _load_X/_save_X accessors
    ("splicecraft_gels", 1),        # [SUB-gels] agarose-gel sim/render: imports biology L0 + util L0
    ("splicecraft_experiments", 1), # [SUB-experiments] entry normalise/id/cross-ref/migrate: imports util L0
    ("splicecraft_backup", 1),      # user-data backup/restore/migrate engine: imports persistence/util/logging/state L0
    ("splicecraft_babs", 1),        # BABS chat engine: Ollama HTTP + HF search + chat protocol; imports net/util/logging/state L0
    ("splicecraft_opentrons", 1),   # OT-2 protocol compiler + robot-server LAN client; imports logging/util/state L0
    ("splicecraft_mapimage", 1),    # circular plasmid-map image export (SVG + PNG); imports biology/util/logging L0
    ("splicecraft_splice", 0),      # plant splice-site PWM scoring; imports logging L0 only
    ("splicecraft_splice_model", 0),  # generated PWM data, no imports (longest-prefix wins over _splice)
    ("splicecraft_cassette", 3),    # eukaryotic expression-cassette assembler; imports biology L0 + dataaccess L1
    ("splicecraft_cli_parser", 0),
    ("splicecraft_biology", 0),
    ("splicecraft_cli", 0),        # stdlib-only sidecar (standalone client)
    ("splicecraft_demo_plasmids", 0),  # pure seed-data module (demo mode)
    ("splicecraft_state", 0),          # shared mutable process state (flags/caches)
    ("splicecraft_logging", 0),        # logging primitives (_log, _log_event, filters)
    ("splicecraft_render", 1),
    ("splicecraft_history", 2),
    ("splicecraft_fileio", 2),      # single-file seq-format I/O: imports record L1 + persistence/logging L0
    ("splicecraft_codon", 2),       # codon optimizer: imports dataaccess L1 + biology/util/logging L0
    ("splicecraft_search", 1),      # online NCBI search ops + HMM-DB downloader: net/util/logging L0 + dataaccess/persistence L1
    ("splicecraft_primer", 2),      # primer/mutagenesis design: imports codon L2 (mutagenesis-inner) + biology/util/logging L0
    ("splicecraft_widgets", 3),
    ("splicecraft_modals", 4),
    ("splicecraft_cloning", 3),  # construction sim: imports history L2/dataaccess L1, used by modals L4
    ("splicecraft_seqanalysis", 3),  # ORF finder + part classifier: imports cloning L3 (same-layer, no cycle) + codon L2/record L1
    ("splicecraft_screens", 6),
    ("splicecraft_agent", 7),
    ("splicecraft_live_ref", 8),
    ("splicecraft_app", 8),
]


def _sibling_modules() -> "list[str]":
    """Every top-level ``splicecraft_*.py`` module name in the repo root."""
    return sorted(
        p.stem
        for p in _REPO.glob("splicecraft_*.py")
        if not p.stem.endswith("_baseline")
    )


def _layer_of(module: str) -> "int | None":
    best = None
    for prefix, layer in _LAYER_RULES:
        if (module == prefix or module.startswith(prefix)) and (
            best is None or len(prefix) > best[0]
        ):
            best = (len(prefix), layer)
    return None if best is None else best[1]


def _sibling_imports(module: str) -> "set[str]":
    """All ``splicecraft*`` modules imported by ``module`` (top-level or nested,
    including the bare ``splicecraft`` hub)."""
    src = (_REPO / f"{module}.py").read_text(encoding="utf-8")
    tree = ast.parse(src, filename=f"{module}.py")
    found: "set[str]" = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split(".")[0]
                if root == "splicecraft" or root.startswith("splicecraft_"):
                    found.add(root)
        elif isinstance(node, ast.ImportFrom):
            if node.level:  # relative import -- siblings must use absolute
                raise AssertionError(
                    f"{module}.py uses a relative import (level={node.level}); "
                    "top-level siblings must import absolutely."
                )
            if node.module:
                root = node.module.split(".")[0]
                if root == "splicecraft" or root.startswith("splicecraft_"):
                    found.add(root)
    return found


def test_every_sibling_is_classified():
    unclassified = [m for m in _sibling_modules() if _layer_of(m) is None]
    assert not unclassified, (
        f"Unclassified splicecraft_* sibling(s): {unclassified}. Add a layer "
        "rule in tests/test_import_layers.py::_LAYER_RULES -- every sibling "
        "needs an explicit layer."
    )


def test_no_sibling_imports_the_hub():
    """Siblings must never import the ``splicecraft`` hub: the hub imports them,
    so this is always a cycle."""
    offenders = {m: i for m in _sibling_modules() if "splicecraft" in (i := _sibling_imports(m))}
    assert not offenders, (
        "Sibling(s) import the `splicecraft` hub, creating an import cycle: "
        f"{sorted(offenders)}. Move the needed symbol down into a lower-layer "
        "sibling, or pass it in -- never import the hub from a sibling."
    )


def test_no_upward_imports():
    violations = []
    for m in _sibling_modules():
        lm = _layer_of(m)
        for dep in _sibling_imports(m):
            if dep == "splicecraft":
                continue  # reported by test_no_sibling_imports_the_hub
            ld = _layer_of(dep)
            if ld is not None and ld > lm:
                violations.append(f"{m} (L{lm}) imports {dep} (L{ld})")
    assert not violations, (
        "Upward import(s) violate the layer rule (a module may import only "
        f"same-or-lower layers):\n  " + "\n  ".join(violations)
    )


def test_no_import_cycles_among_siblings():
    graph = {m: {d for d in _sibling_imports(m) if d != "splicecraft"}
             for m in _sibling_modules()}
    visiting, done = set(), set()
    stack: "list[str]" = []

    def visit(node: str):
        if node in done:
            return
        if node in visiting:
            cycle = stack[stack.index(node):] + [node]
            raise AssertionError("Import cycle among siblings: " + " -> ".join(cycle))
        visiting.add(node)
        stack.append(node)
        for nxt in graph.get(node, ()):  # ignore deps that aren't local siblings
            if nxt in graph:
                visit(nxt)
        stack.pop()
        visiting.discard(node)
        done.add(node)

    for m in list(graph):
        visit(m)


def test_all_siblings_are_packaged():
    """Every splicecraft_*.py sibling the hub re-imports MUST appear in
    pyproject.toml's wheel `only-include` AND sdist `include`, or the released
    package breaks: `import splicecraft` fails at the re-import step for end
    users (it works in the source tree, so the rest of the suite wouldn't catch
    it). Guards the broken-wheel bug the flat-sibling refactor could introduce."""
    pyproject = tomllib.loads((_REPO / "pyproject.toml").read_text(encoding="utf-8"))
    siblings = {p.name for p in _REPO.glob("splicecraft_*.py")}
    targets = pyproject["tool"]["hatch"]["build"]["targets"]
    wheel = set(targets["wheel"]["only-include"])
    sdist = set(targets["sdist"]["include"])
    assert not (siblings - wheel), (
        "siblings missing from pyproject wheel only-include (would ship a broken "
        f"wheel): {sorted(siblings - wheel)}"
    )
    assert not (siblings - sdist), (
        f"siblings missing from pyproject sdist include: {sorted(siblings - sdist)}"
    )
