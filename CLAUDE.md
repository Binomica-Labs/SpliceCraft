# CLAUDE.md — AI Agent Context for SpliceCraft

Agent handoff. Read before touching the codebase.

## ⚠ SACRED: data-dir safety (READ FIRST)

The user's plasmid library + collections + primers + parts live in `~/.local/share/splicecraft/` (or `$XDG_DATA_HOME/splicecraft/`). **The data is the product.** A wrong write here destroys hours-to-years of user work. Three hard rules:

1. **Never `import splicecraft` from an ad-hoc script (`/tmp/*.py`, REPL, probe) without first sandboxing `XDG_DATA_HOME`.** `_DATA_DIR` is computed at import time and won't budge afterwards. Sandbox by:
   ```python
   import os, tempfile; os.environ["XDG_DATA_HOME"] = tempfile.mkdtemp(prefix="sc-")
   os.environ.setdefault("SPLICECRAFT_SKIP_LOCK", "1")
   import splicecraft as sc
   assert "sc-" in str(sc._DATA_DIR), f"unsandboxed: {sc._DATA_DIR}"
   sc._authorize_writes_for_sandbox(sc._DATA_DIR)   # L2 chokepoint opt-in
   ```
2. **`_save_*` helpers are nuclear-coded.** Calling `_save_collections`, `_save_library`, `_save_primers`, `_save_parts_bin`, `_save_features`, `_save_grammars`, `_save_entry_vectors`, `_save_codon_tables`, `_save_protein_motifs`, `_save_experiments`, `_save_experiment_projects`, `_save_gels`, `_save_custom_enzymes`, `_save_enzyme_collections`, or `_safe_save_json` directly from outside the four sanctioned callers (`PlasmidApp.main()`, pytest `_protect_user_data` fixture, agent HTTP server, sandboxed verifier harness) raises `RuntimeError` since the L2 chokepoint landed — sandbox first or use the GUI.
3. **Verifier scripts always go through `.claude/skills/verifier-splicecraft.md`.** It enforces the sandbox + authorization at the top. Don't roll your own.

**Caught failure (2026-05-22):** an unsandboxed `/tmp/sc_probe.py` ran `_save_collections([{"name": "Default", "plasmids": []}])` for test setup. It wrote directly to the user's real 160 MB `collections.json`, rotating the previous good state to `.bak`. The four-layer safety net + lost-entries spillover recovered the data, but the lesson stands: there is NO "I'll be careful this once" version of writing to the data dir. Sandbox or refuse.

---

Bioinformatician + Claude. **Near-single-file architecture** — `splicecraft.py` (~105k lines) + extracted biology module `splicecraft_biology.py` + stdlib-only sidecar `splicecraft_cli.py`. Single-file constraint is intentional (greppable); biology extraction is the first deliberate exception (pure functions/constants, no `PlasmidApp` coupling). See `CONTRIBUTING.md` three-test rule.

## What is SpliceCraft?

Terminal-based circular plasmid map viewer, sequence editor, cloning/mutagenesis workbench. Python 3.10+ / Textual / Biopython. Unicode braille-dot maps, per-base sequence panel, restriction overlays, collection-driven library, Golden Braid L0 + MoClo grammars, Primer3, SOE-PCR mutagenesis, in-process BLASTN/BLASTP/HMMscan via pyhmmer.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` · **PyPI:** `splicecraft` · `__version__` in `splicecraft.py` and `pyproject.toml`.

## How to run

```bash
python3 splicecraft.py                       # empty canvas (or auto-loads first library entry)
python3 splicecraft.py L09137                # fetch pUC19 from NCBI
python3 splicecraft.py myplasmid.gb          # local GenBank (.gb/.gbk/.dna)
python3 -m pytest -n auto -q                 # full suite (~5–6 min on 8 cores)
python3 -m pytest tests/test_dna_sanity.py   # biology only (<2 s, fast inner loop)
./release.py X.Y.Z                           # bump, test, build, tag, push (PyPI via OIDC)
```

End users: `pipx install splicecraft && splicecraft`.

No-arg launch shows empty canvas (or first library entry). Demo plasmid (`_make_demo_record` / `_DEMO_PLASMID_SEQ`) kept in source for tests but `main()` no longer pre-sets `_preload_demo_record`. First-run NCBI seed (`_seed_default_library` → MW463917.1) suppressed via `_skip_seed = True`. Dev builds flip `_skip_seed = False` for auto-seed.

Logs: `~/.local/share/splicecraft/logs/splicecraft.log` (override `$SPLICECRAFT_LOG`). 8-char session ID prefix per line.

## Where to find more (grep first, ahead of dispatch)

The long-form rules and subsystem deep-dives live in split files. **Each entry has a tag — grep before editing the matching subsystem.** Each file starts with a tag→topic table.

| File | Holds | Grep when touching |
|---|---|---|
| `docs/invariants.md` | `[PIT-01]`…`[PIT-35]` known pitfalls; `[INV-36]`…`[INV-86]` sweep history; `[PREFS]` settings; `[ARCH]` pointers; `[CONV]` conventions; `[SISTER]` ScriptoScope; `[RECIPE]` new-feature playbook | bare-except, wrap features, cache contracts, agent endpoints, master delete, synthesis/protein, collision modals, settings persistence, NEW FEATURE checklists |
| `docs/subsystems.md` | `[SUB-plasmidsaurus]`, `[SUB-experiments]`, `[SUB-gels]` | sequencing zip ingestion, lab notebook, gel snapshots |
| `docs/architecture.md` | Single-file rationale, test pyramid, concurrency model, observability | high-level structural decisions |
| `docs/PLATFORMS.md` | Supported OS / terminal matrix | cross-platform behaviour, terminal capability checks |
| `docs/agent-api.md`, `docs/features.md` | User-facing reference docs | adding/renaming endpoints, user-visible features |

**Rule:** before any non-trivial edit, run `grep -ni '<keyword>' docs/invariants.md docs/subsystems.md`. If a relevant `[PIT-NN]` / `[INV-NN]` / `[SUB-xxx]` exists, read it first. Dispatching a sub-agent? Quote the matching tag in the prompt so it knows where to look.

## Sacred invariants (DO NOT BREAK)

Each has at least one test in `tests/`. Touching `_scan_restriction_sites`, `_rc`, `_iupac_pattern`, `_translate_cds`, `_bp_in`, `_feat_len`, the wrap-midpoint formula, or `_rebuild_record_with_edit` trips tests immediately.

1. **Palindromic enzymes scanned forward only.** Bottom-strand hit emitted as `recut`. Scanning both strands double-counts.
2. **Reverse-strand resite positions use forward coordinate.** Reverse hit at `p` (after RC) stored as `p`, not `n - p - site_len`. Cut column maps via `rev_cut_col = site_len - fwd_cut`.
3. **`_rc()` handles full IUPAC** via `_IUPAC_COMP`, not just ACGT.
4. **IUPAC regex patterns cached** in `_PATTERN_CACHE`.
5. **Circular wrap midpoint:** `arc_len = (end - start) % total; mid = (start + arc_len // 2) % total`. Naive form puts label opposite actual arc.
6. **Circular wrap RE scan** scans `seq + seq[:max_site_len-1]`. Each wrap hit emits **two resite pieces** (labeled tail `[p, n)` + unlabeled head `[0, (p+site_len) - n)`) and **one recut** at `(p + fwd_cut) % n`. Resite-counting code must count only labeled pieces.
7. **Data-file saves always back up.** Always go through `_safe_save_json` (`.bak` + `tempfile.mkstemp` + `os.fsync` + `os.replace`). Schema envelope `{"_schema_version": 1, "entries": [...]}`; `_extract_entries` accepts legacy bare-list (pre-0.3.1). **`_safe_save_json` re-raises on failure** so callers can notify — silent swallow used to desync UI from disk.
8. **Wrap-aware feature length.** Use `_feat_len(start, end, total)` — returns `(total - start) + end` when `end < start`. All sort keys, length displays, biology checks route through it.
9. **Wrap-feature integrity in record edits.** `int(CompoundLocation.start)` returns `min(parts.start)` and silently flattens. `_rebuild_record_with_edit` per-part shifts wrap features and only collapses to FeatureLocation when 1 part survives.
10. **Undo snapshots deepcopied.** `_push_undo`, `_action_undo`, `_action_redo` all `deepcopy(self._current_record)`.

> Known pitfalls #1–#35 live in `docs/invariants.md` as `[PIT-01]`…`[PIT-35]`. Persistent-settings recipe → `[PREFS]`. Architecture pointers → `[ARCH]`. Conventions → `[CONV]`. Sister-project crib → `[SISTER]`. New-feature playbook ("Borrow before respinning") → `[RECIPE]`.

## For future agents

1. Read this file first, then **grep `docs/invariants.md` + `docs/subsystems.md` for the area you're touching** before any edit, then `git log --oneline` for recent context.
2. `python3 -m pytest -n auto -q` before and after any change. `tests/test_dna_sanity.py` (<2 s) is the fast inner loop.
3. Don't break sacred invariants. Don't bypass `_safe_save_json`. Don't add bare `except`.
4. Eyeball real-world behaviour on pUC19 (`L09137`) and pACYC184 (`MW463917.1`).
5. Past fix history is in git — `git show <hash>` beats stale prose.
6. **Dispatching a sub-agent?** Quote the relevant `[PIT-NN]` / `[INV-NN]` / `[SUB-xxx]` / `[RECIPE]` tag in its prompt so it greps the right file before working.
