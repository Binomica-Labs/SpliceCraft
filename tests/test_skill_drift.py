"""
test_skill_drift — the Claude Skill may not drift from the code.

The skill at `plugin/skills/splicecraft/` is how an agent learns this API, and
it is shipped to the Anthropic plugin directory, so a wrong statement in it is
a wrong statement in front of every user of the plugin. A 2026-10-16 audit
found six: `apply-gff3` documented as dry-run when it writes on the first
call, `recover-history-from-dna` documented as taking `apply` when it reads
`dry_run`, two endpoint names that do not exist
(`protocol-from-history`, `import-plasmidsaurus-*`), a three-mode
`design-primers` that has four, and stale endpoint counts in five places.

Every one of those is what this file catches. `docs/agent-endpoints.md` has
had a two-way drift test since `[INV-207]`; the skill had nothing.

The rule is deliberately asymmetric:
  * every endpoint-shaped name the skill mentions MUST exist (it is telling an
    agent to call it), but
  * the skill need NOT mention every endpoint — it is a guide, and SKILL.md
    says in as many words that the live `tools` list outranks it.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

import splicecraft as sc

_SKILL_DIR = Path(__file__).resolve().parent.parent / "plugin" / "skills" / "splicecraft"
_FILES = ["SKILL.md", "references/endpoints.md", "references/workflows.md"]


def _registry() -> "set[str]":
    return set(sc._state._AGENT_HANDLERS)


def _skill_text() -> "dict[str, str]":
    out = {}
    for rel in _FILES:
        p = _SKILL_DIR / rel
        if p.exists():
            out[rel] = p.read_text(encoding="utf-8")
    return out


# Backticked lower-kebab tokens. An endpoint name is `[a-z0-9]+(-[a-z0-9]+)+`,
# which also matches ordinary hyphenated prose in backticks (`read-only`,
# `end-exclusive`), so the comparison below is against the registry and a
# curated allowlist of non-endpoint tokens rather than a bare "looks like an
# endpoint" guess — a heuristic that fired on prose would make this test noise
# and noise gets deleted.
_TOKEN_RE = re.compile(r"`([a-z0-9]+(?:-[a-z0-9]+)+)`")

# Hyphenated prose, CLI words, file/format names and settings keys that are
# legitimately backticked but are not endpoints. Anything NOT here and NOT in
# the registry fails, so a genuinely bogus endpoint name cannot hide.
_NOT_ENDPOINTS = frozenset({
    # CLI + shell
    "splicecraft-cli", "read-only", "end-exclusive", "pipe-separated",
    "comma-separated", "self-describe", "round-trip", "round-trips",
    "dry-run", "drop-in", "built-in", "built-ins", "follow-up", "add-on",
    "well-formed", "case-insensitive", "cross-collection", "in-process",
    "in-memory", "on-disk", "per-base", "per-read", "per-cluster",
    "per-endpoint", "sub-population", "whole-codon", "single-read",
    "no-reads", "minor-variants", "name-hit", "doc-hit", "0-based",
    "1-based", "double-stranded", "reverse-complement-of",
    "lower-kebab", "two-way", "four-letter", "high-GC", "low-quality",
    "trace-quality", "read-modify-write", "all-or-nothing",
    "exact-match", "wild-type", "amino-acid", "order-cart",
    "primer-only", "map-feature", "non-numeric", "server-side",
    "client-side", "loop-back", "write-twin",
})


def test_skill_files_exist():
    """The skill must be where CLAUDE.md says it is — `plugin/skills/`, not
    `.claude/skills/` (the path is fixed by the plugin format plus the
    directory's 5 MiB per-file ceiling)."""
    missing = [rel for rel in _FILES if not (_SKILL_DIR / rel).exists()]
    assert not missing, f"skill file(s) missing: {missing}"


def test_every_endpoint_name_in_the_skill_exists():
    """THE test. A name the skill tells an agent to call must resolve in the
    live handler registry, or the agent gets a 404 and reports the feature as
    broken."""
    live = _registry()
    bogus: "dict[str, set[str]]" = {}
    for rel, text in _skill_text().items():
        for tok in set(_TOKEN_RE.findall(text)):
            if tok in live or tok in _NOT_ENDPOINTS:
                continue
            bogus.setdefault(rel, set()).add(tok)
    assert not bogus, (
        "The skill names tokens that are neither live endpoints nor known "
        f"prose: {({k: sorted(v) for k, v in bogus.items()})}\n"
        "If one is a renamed/removed endpoint, fix the skill. If it is "
        "ordinary hyphenated prose, add it to _NOT_ENDPOINTS."
    )


def test_stated_endpoint_counts_match_the_registry():
    """The skill states a total in several places. A stale count is how an
    agent concludes the surface is smaller than it is — the `[INV-202]`
    failure one layer out."""
    n = len(_registry())
    stale: "list[str]" = []
    for rel, text in _skill_text().items():
        for m in re.finditer(r"(\d{3})\s+endpoints", text):
            if int(m.group(1)) != n:
                stale.append(f"{rel}: '{m.group(0)}' (registry has {n})")
        for m in re.finditer(r"~(\d{3})\s+endpoints", text):
            if int(m.group(1)) != n:
                stale.append(f"{rel}: '~{m.group(1)}' (registry has {n})")
        for m in re.finditer(r"×\s*(\d{3})", text):
            if int(m.group(1)) != n:
                stale.append(f"{rel}: '× {m.group(1)}' (registry has {n})")
    assert not stale, "Stale endpoint counts in the skill: " + "; ".join(stale)


def test_dry_run_table_matches_the_handlers():
    """The commit-key table in SKILL.md §3 is the one an agent follows to make
    a write actually happen. Each listed endpoint must really read the key the
    table names — sending the wrong one is a silent no-op reported as done.

    `apply-gff3` is the case that went wrong: documented as dry-run, it writes
    on the first call and has no `applied` field at all."""
    import inspect
    live = sc._state._AGENT_HANDLERS
    # (endpoint, key it must read)
    expect = [
        ("transfer-annotations", "apply"),
        ("annotate-from-presets", "apply"),
        ("recover-qualifiers-from-dna", "apply"),
        ("recover-history-from-dna", "dry_run"),
    ]
    problems = []
    for name, key in expect:
        entry = live.get(name)
        if entry is None:
            problems.append(f"{name}: endpoint missing")
            continue
        src = inspect.getsource(entry[0])
        if f'"{key}"' not in src:
            problems.append(f"{name}: does not read {key!r}")
    # apply-gff3 must NOT read either key — it is not dry-run.
    gff = live.get("apply-gff3")
    if gff is not None:
        src = inspect.getsource(gff[0])
        for key in ("payload.get(\"apply\")", "payload.get(\"dry_run\")"):
            if key in src:
                problems.append(
                    f"apply-gff3 now reads {key} — it became dry-run, so the "
                    f"commit-key table in SKILL.md §3 must be updated")
    assert not problems, "; ".join(problems)


def test_egress_table_names_both_settings():
    """There are TWO independent egress settings and the skill used to merge
    them, so a relayed refusal named the wrong toggle half the time."""
    text = (_SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    for key in ("allow_online_search", "allow_online_lookups"):
        assert key in text, f"SKILL.md must name {key} in its egress section"


@pytest.mark.parametrize("mode", ["detection", "cloning", "generic",
                                   "allele_specific"])
def test_design_primers_modes_are_documented(mode):
    """All four modes. `allele_specific` was invisible in the skill, so the
    SNP/genotyping mode could not be found by an agent reading it."""
    text = (_SKILL_DIR / "references" / "endpoints.md").read_text(
        encoding="utf-8")
    assert mode in text, f"design-primers mode {mode!r} is undocumented"
