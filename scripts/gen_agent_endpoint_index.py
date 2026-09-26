#!/usr/bin/env python3
"""Regenerate `docs/agent-endpoints.md` from the live handler registry.

`docs/agent-api.md` is the prose guide and is hand-written; it documented
roughly half the endpoints, and the half it missed was invisible to anyone
reading the docs (audit 2026-09-22). This generates the COMPLETE list from
`_state._AGENT_HANDLERS`, so a new endpoint cannot be absent from the docs
without `tests/test_audit_2026_09_22b.py` failing.

Sandboxes `XDG_DATA_HOME` before importing splicecraft — see CLAUDE.md's
data-dir rules; this script only reads, but the rule is the rule.
"""
from __future__ import annotations

import inspect
import os
import re
import sys
import tempfile
from pathlib import Path

os.environ["XDG_DATA_HOME"] = tempfile.mkdtemp(prefix="sc-docgen-")
os.environ.setdefault("SPLICECRAFT_SKIP_LOCK", "1")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import splicecraft as sc  # noqa: E402

_SUMMARY_MAX = 150

HEADER = """# Agent API — complete endpoint index

Every registered endpoint, generated from the live handler registry
(`_state._AGENT_HANDLERS`) so it cannot drift from the code. `R` reads,
`W` writes (a write also takes the canvas dirty-guard). The prose guide,
with request/response shapes and worked examples, is `agent-api.md`;
this file exists so no endpoint is invisible.

Regenerate with `python3 scripts/gen_agent_endpoint_index.py`.

"""


def _summary(handler) -> str:
    doc = inspect.getdoc(handler) or ""
    if not doc:
        return ""
    first = " ".join(doc.split("\n\n")[0].split())
    first = re.sub(r"``([^`]+)``", r"`\1`", first)
    if len(first) > _SUMMARY_MAX:
        first = first[:_SUMMARY_MAX - 3].rsplit(" ", 1)[0] + "…"
    return first.replace("|", "\\|")


def build() -> str:
    rows = []
    for name in sorted(sc._state._AGENT_HANDLERS):
        entry = sc._state._AGENT_HANDLERS[name]
        handler = entry[0] if isinstance(entry, tuple) else entry
        handler = handler if callable(handler) else getattr(
            handler, "handler", None)
        write = "W" if (isinstance(entry, tuple) and len(entry) > 1
                        and entry[1]) else "R"
        rows.append((name, write, _summary(handler)))
    lines = ["| Endpoint | R/W | What it does |", "|---|---|---|"]
    lines += [f"| `{n}` | {w} | {d} |" for n, w, d in rows]
    return HEADER + "\n".join(lines) + "\n"


if __name__ == "__main__":
    out = Path(__file__).resolve().parent.parent / "docs" / "agent-endpoints.md"
    out.write_text(build())
    n = len(sc._state._AGENT_HANDLERS)
    print(f"wrote {out} ({n} endpoints)")
