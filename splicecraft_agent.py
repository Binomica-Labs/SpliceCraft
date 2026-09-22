"""splicecraft_agent — the data-only agent-API HTTP endpoint handlers (Phase D, layer L7).

The pure / data-only `_h_*` endpoint handlers (read queries + data mutations that
go through the sibling accessors, NOT the running app's widgets) extracted from
the hub, plus the agent-API helper layer they share: the `_agent_endpoint`
registration decorator, the response dict-builders (`_agent_*_dict` /
`_custom_enzyme_meta` / `_parts_bin_entry_summary` / `_history_node_to_dict`),
payload coercers / validators (`_coerce_int` / `_sanitize_bases` /
`_settings_validator_*` / `_agent_validate_custom_enzyme_payload`), the agent
path-safety checks, and the write-guard (`_agent_dirty_guard` / `_agent_save_or_500`).

The handlers register into `_state._AGENT_HANDLERS` (the shared registry) via the
decorator at import; the hub's HTTP server dispatches from that same dict. The
APP-COUPLED handlers (those that `query_one(PlasmidMap/SequencePanel/LibraryPanel)`
or touch `app._current_record` / `app.call_from_thread`) STAY hub-side and register
into the SAME dict through the re-exported decorator.

Reaches hub-pinned bits via _state: `_AGENT_HANDLERS` (registry), `_LIVE_APP_REF`
(live-app singleton, for save-failure notify), `_resolve_data_attr_hook`,
`_all_enzymes_hook`, `_cache_lock`, plus 12 deferred-handler hooks for the deep
engines that stay hub-side (in-process BLAST/HMM via pyhmmer, the alignment-
rotation picker, the PCR sim, settings-flush, the master-delete + assembly-
fragment cache busts, the experiment blob delete, the entry-vector binder). The
"deferred chase" relocated the movable blockers — active-name getters/finders →
dataaccess, the restore-from-backup engine → backup, the PCR caps → cloning — so
ALL 109 data-only endpoints live here; only the 26 app-coupled handlers stay hub.
`_sanitize_path` / `_ENZYME_CUT_RANGE` were relocated to L0 (util / biology) as
prerequisites. Top layer among the siblings (imports the whole domain stack ≤L3);
nothing imports it but the hub.
"""
from __future__ import annotations

import inspect
import re
from datetime import date as _date, datetime as _datetime
from pathlib import Path

import splicecraft_state as _state
import splicecraft_opentrons as _ot2
from splicecraft_dataaccess import (_load_custom_labware, _load_protocol_collections,
                                    _save_custom_labware, _save_protocol_collections)
from splicecraft_backup import (_list_pre_update_snapshots)
from splicecraft_biology import (_ENZYME_CUT_RANGE, _assemble_operon, _feat_len, _iupac_pattern, _rbs_design, _rbs_strength, _rc, _rna_cofold, _rna_fold, _seq_len, _span_in_span)
from splicecraft_cloning import (_GIBSON_MAX_OVERLAP_BP, _GIBSON_MIN_OVERLAP_BP, _excise_fragment_pair, _excise_pcr_insert, _scrub_gb_design, _simulate_gibson_assembly, _simulate_golden_gate, _simulate_traditional_cloning_multi, _enzyme_is_type_iis)
from splicecraft_codon import (_CODON_RAMP_CODONS, _CODON_RARE_W, _codon_relative_adaptiveness, _CODON_GC3_HIGH, _CODON_GC3_LOW, _CODON_GC3_MIN_CODONS, _CODON_GC_WINDOW_DEFAULT, _CODON_GENETIC_CODE, _CODON_MODES, _CODON_REPEAT_RUN_DEFAULT, _CODON_SCRUB_MAX_CODONS, _codon_cai, _codon_diversify, _codon_fetch_kazusa, _codon_fix_gc_window, _codon_fix_sites, _codon_gc, _codon_gc3, _codon_gc_window_range, _codon_hazard_motifs, _codon_kmer_set, _codon_optimize, _codon_shared_runs, _codon_tables_add, _file_build_codon_table, _genome_build_codon_table)
from splicecraft_dataaccess import (_BUILTIN_GRAMMARS, _all_grammars, _clear_entry_vectors_for_grammar, _codon_tables_get, _codon_tables_load, _codon_tables_save, _find_gel, _find_hmm_db_entry, _find_library_entry_by_id, _get_active_collection_name, _get_active_primer_collection_name, _get_entry_vector, _get_setting, _hmm_db_name_taken, _iter_all_experiments, _iter_collections_readonly, _iter_library_readonly, _iter_parts_bin_readonly, _load_custom_enzymes, _load_custom_grammars, _load_entry_vectors, _load_enzyme_collections, _load_experiment_projects, _load_experiments, _load_feature_colors, _load_features, _load_gels, _load_hmm_db_catalog, _load_library, _load_parts_bin, _load_primer_collections, _load_primers, _load_protein_motifs, _normalise_hmm_db_entry, _sanitize_hmm_db_id, _sanitize_hmm_db_url, _save_custom_enzymes, _save_custom_grammars, _save_enzyme_collections, _save_experiment_projects, _save_experiments, _save_feature_colors, _save_features, _save_gels, _save_hmm_db_catalog, _save_library, _save_parts_bin, _save_parts_bin_collections, _save_primer_collections, _save_primers, _save_protein_motifs, _search_collections_library, _set_active_primer_collection_name, _set_entry_vector, _set_setting, _typed_clone)
from splicecraft_experiments import (_EXPERIMENT_ACTIONS, _experiment_duplicate, _experiment_html_document, _experiment_markdown_document, _experiment_project_markdown, _experiment_ref_token_ok, _experiment_search, _experiment_snippet, _experiment_search_terms, _experiment_template_markdown, _experiments_referencing, _new_experiment_id, _normalise_experiment_entry, _protocol_steps_markdown, _sanitize_experiment_id)
from splicecraft_fileio import (_PLASMIDSAURUS_ZIP_MAX_BYTES, _export_commercialsaas_dna, _export_embl_to_path, _list_gbk_members_in_zip, _parse_commercialsaas_history, _plasmidsaurus_zip_to_entries)
from splicecraft_gels import (_new_gel_id, _normalise_gel_entry)
from splicecraft_history import (_HISTORY_NODE_MAX_DEPTH, _HISTORY_NODE_MAX_NODES, _history_build_steps, _history_node_warnings)
from splicecraft_crispr import (  # noqa: E402
    _CAS_VARIANTS, _cas_variant, _design_guides, _guide_cloning_oligos,
    _guide_offtargets, _score_guide,
)
from splicecraft_presets import (
    _find_preset, _preset_categories, _preset_features, _preset_matches,
    _preset_to_library_entry,
)
from splicecraft_logging import (_log, _log_event)
from splicecraft_net import (_sanitize_accession)
from splicecraft_persistence import (_safe_file_size_check, _safe_load_json)
from splicecraft_primer import (_mut_design_inner, _mut_design_outer, _scrub_design, _scrub_qc_primers, _scrub_qc_verify)
from splicecraft_record import (_gb_text_to_record, _normalize_primer_seq,
                                _topology_from_gb_text)
from splicecraft_search import (_ONLINE_LOOKUP_MAX_HITS, _ONLINE_LOOKUP_QUERY_MAX, _PLASMIDSAURUS_ITEMS_LIMIT, _PLASMIDSAURUS_ITEMS_TRUNCATED_HINT, _PLASMIDSAURUS_RESULT_KINDS, _delete_hmm_db_files, _europepmc_search, _fpbase_search, _hmm_db_acquire_download_slot, _hmm_db_perform_download, _hmm_db_pressed, _hmm_db_release_download_slot, _hmmer_web_hmmscan, _ncbi_blast_db_for, _ncbi_blast_online, _ncbi_db_search, _online_clean_query, _online_max_query_len, _patent_search, _plasmidsaurus_credentials, _plasmidsaurus_fetch_item_zip, _plasmidsaurus_item_has_results, _plasmidsaurus_list_items, _plasmidsaurus_oauth_token, _read_url, _sanitize_plasmidsaurus_item_code, _uniprot_search, _web_search, _wikipedia_search)
from splicecraft_seqanalysis import (_classify_part_from_plasmid, _ev_frag_input_features, _find_orfs, _fragment_has_backbone_marker, _synthesis_lint, _fragment_backbone_marker_labels, _predict_transcript, _map_transcription)
from splicecraft_util import (_PLASMID_STATUS_VALUES, _check_export_extension, _feat_bounds, _feat_label, _normalize_collection_name, _notify_save_failure, _primer_tm_safe, _record_is_circular, _safe_color_for_write, _sanitize_feat_type, _sanitize_gel_id, _sanitize_label, _sanitize_note, _sanitize_path, _scrub_path)
from splicecraft_widgets import (_PLASMID_STATUS_COLORS)
from splicecraft_backup import (_AGENT_BACKUP_LABELS, _PRE_UPDATE_NAME_RE, _export_migrate_archive, _list_recoverable_backups, _resolve_backup_label, _restore_from_backup, _restore_pre_update_snapshot)
from splicecraft_biology import (_digest_with_enzymes, _enzyme_aliases, _enzyme_cuts, _enzyme_resolve_one, _enzyme_signature, _resolve_enzyme_names, _scan_restriction_sites)
from splicecraft_cloning import (_PCR_AMPLICON_HARD_CAP, _PCR_DEFAULT_MAX_AMPLICON, _PCR_MAX_AMPLICONS, _PCR_MAX_PRIMER_LEN, _PCR_MAX_TEMPLATE_BP, _PCR_MIN_PRIMER_LEN, _build_synthesis_l0_fragment, _design_gb_primers, _entry_vector_acceptor_overhangs, _grammar_position_by_type, _l0_part_from_syn_fragment)
from splicecraft_dataaccess import (_active_enzyme_allowed_set, _agent_scan_library_for_key, _find_enzyme_collection, _find_library_entry_by_name, _find_parts_bin, _find_project, _get_active_enzyme_collection_name, _get_active_parts_bin_name, _get_active_project_name, _load_parts_bin_collections, _set_active_enzyme_collection_name, _set_active_parts_bin_name, _set_active_project_name)
from splicecraft_fileio import (_export_fasta_to_path, _export_genbank_to_path, _export_gff_to_path, _extract_gbk_member, _fastq_path_to_records)
from splicecraft_gels import (_AGAROSE_CHOICES, _GEL_HEIGHT_MAX, _GEL_HEIGHT_MIN, _GEL_LANE_WIDTH_MAX, _GEL_LANE_WIDTH_MIN, _GEL_MAX_LANES, _agarose_mobility, _gel_bands_for_lane, _render_gel_image)
from splicecraft_persistence import (_safe_save_json_mirror)
from splicecraft_regulatory import (_PROMOTER_MIN_SCORE, _TERM_MAX_STEM, _TERM_MIN_STEM, _TERM_MIN_U_TRACT, _TERM_U_TRACT_WINDOW, _scan_promoters, _scan_terminators)
from splicecraft_primer import (_AS_MIN_LEN, _design_allele_specific_primer, _design_cloning_primers_raw, _design_detection_primers, _design_generic_primers, _primer_binding_sites, _primer_check_confidence, _primer_tm)


def _custom_enzyme_meta(name: str) -> "dict | None":
    """Return the full custom-enzyme dict (incl. type, supplier) for
    ``name`` if present, else None. Used by the modal to show the
    extra columns."""
    if not isinstance(name, str) or not name:
        return None
    for entry in _load_custom_enzymes():
        if isinstance(entry, dict) and entry.get("name") == name:
            return entry
    return None


def _check_agent_read_path_ancestors(path: Path) -> "str | None":
    """Tighten read-endpoint path validation by walking every parent
    component for symlinks. Mirrors the sweep-#4 hardening in
    ``_check_agent_write_path`` but for read endpoints (Plasmidsaurus
    zip ingestion, etc.).

    Returns an error message string when an ancestor symlink would
    redirect the read; None when safe. Sweep #26 (2026-05-23) —
    closes the gap audit M8 flagged: read endpoints rejected the
    target's `~user` expansion but a parent like `Documents/zips/`
    being a symlink to `/etc` would still be followed silently by
    the downstream ``os.open``.

    The path itself does NOT need to exist — symlink refusal applies
    only to existing components.
    """
    try:
        parent = path.parent
        if not parent.exists():
            # Nothing to redirect through if it doesn't exist yet.
            return None
        resolved_parent = parent.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        return f"could not resolve parent directory: {exc}"
    try:
        lexical_parent = parent.absolute()
    except OSError as exc:
        return f"could not normalise parent directory: {exc}"
    if str(resolved_parent) != str(lexical_parent):
        return (
            f"parent path resolves through a symlink: "
            f"{lexical_parent!s} → {resolved_parent!s}"
        )
    return None


def _check_agent_write_path(path: Path) -> "str | None":
    """Tighter validation for agent write endpoints (`export-*`,
    `save`, etc.). Returns an error message string when the path is
    rejected; None when safe to write.

    Rejects:
    * Symlinks at the destination — an agent shouldn't get to write
      through a pre-placed symlink to `/etc/passwd` or similar.
    * Existing symlinks in ANY parent component up to root (TOCTOU
      defense — a racing process can't redirect the write via a
      grandparent symlink either).
    * Paths whose parent doesn't exist (forces the user to mkdir
      first rather than us auto-creating arbitrary directories).

    Audit hardening 2026-05-14: previously the agent's export
    endpoints used `_sanitize_path` only, which expands `~` and does
    nothing else — symlink-as-destination was unprotected.

    Audit sweep #4 2026-05-15: previously this only checked the
    immediate parent. A symlink at any deeper ancestor (e.g.
    `/home/<user>/Documents` → `/etc`) could redirect every write
    under it. Walk the full chain via `resolve()` so the canonical
    absolute path is what we compare; an ancestor symlink would
    yield a `resolve()` result differing from the lexical path,
    flagging the redirect.
    """
    if path.is_symlink():
        return f"refusing to write through symlink at {path!s}"
    parent = path.parent
    if not parent.exists():
        return f"parent directory does not exist: {parent!s}"
    # Resolve the parent through every symlink hop. If the result
    # differs from the lexical absolute path (modulo `..` / `.`
    # collapsing), some ancestor is a symlink.
    try:
        resolved_parent = parent.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        return f"could not resolve parent directory: {exc}"
    # Walk every ancestor segment under the parent and lstat
    # each. `parent.resolve()` follows symlinks, so if any
    # intermediate component IS a symlink, `lexical_parent` and
    # `resolved_parent` diverge. Refuse the divergence outright —
    # the user can re-target through the resolved path if they
    # actually meant to write there.
    try:
        lexical_parent = parent.absolute()
    except OSError as exc:
        return f"could not normalise parent directory: {exc}"
    if str(resolved_parent) != str(lexical_parent):
        return (
            f"parent path resolves through a symlink: "
            f"{lexical_parent!s} → {resolved_parent!s}"
        )
    # Per-segment lstat as defense in depth — `resolve()` already
    # catches symlinks via the divergence check above, but an
    # attacker swapping a regular dir for a symlink between the
    # resolve() and the open() (TOCTOU race) would slip through.
    # Walking each ancestor with lstat at least narrows the race
    # window to the open() itself.
    cur = parent
    seen: set = set()
    while True:
        try:
            if cur.is_symlink():
                return (
                    f"ancestor directory is a symlink: {cur!s}"
                )
        except OSError:
            # A permission error mid-walk is a refusal — we can't
            # tell whether an ancestor is safe.
            return f"could not stat ancestor: {cur!s}"
        if cur.parent == cur or str(cur) in seen:
            break
        seen.add(str(cur))
        cur = cur.parent
    return None


def _coerce_bool(value, *, name: str = "value") -> "bool | str":
    """Type-safe bool coercion for agent-API JSON payloads.

    Returns the bool on success, or a human-readable error message (string) on
    failure — the same `value | str` shape as `_coerce_int`, so an
    `isinstance(result, str)` guard narrows it.

    A bare `bool(value)` is the trap this exists to close: `bool("false")` is
    True, so a caller passing the STRING "false" — which shell-built JSON and
    hand-edited config do constantly — silently gets the opposite of what they
    asked for. That is merely untidy for a cosmetic flag and dangerous for a
    flag like `homopolymer_filter`, whose whole purpose is to let someone turn
    OFF a filter that could otherwise hide the variant they are hunting.

    Accepts real booleans, the usual textual spellings, and 0/1. Refuses
    anything else rather than guessing.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        v = value.strip().lower()
        if v in ("true", "yes", "on", "1"):
            return True
        if v in ("false", "no", "off", "0"):
            return False
    return (f"'{name}' must be true or false "
            f"(got {type(value).__name__} {value!r})")


def _coerce_int(value, *, name: str = "value") -> "int | str":
    """Type-safe int coercion for agent-API JSON payloads.

    Returns the integer on success, or a human-readable error message
    (string) on failure. Accepts ``bool`` / ``int`` / finite ``float``
    / digit ``str``. Rejects ``None`` / dict / list / NaN / +-Inf —
    all of which would either AttributeError on a downstream `.get`
    or raise ``OverflowError`` on a naked ``int(value)`` (the case
    that bit us when an agent sent ``{"max_hits": Infinity}`` and a
    downstream ``int(...)`` blew up before our range check).

    The return shape is ``int | str`` rather than the older
    ``tuple[int | None, str | None]`` so that a caller-side
    ``isinstance(result, str)`` guard narrows the value to ``int``
    automatically — no separate ``assert value is not None`` needed
    and no tuple unpacking that loses pyright's discriminated-union
    narrowing.
    """
    import math as _m
    if isinstance(value, bool):
        # `True` is not 1 for an API that takes coordinates and indices.
        # Coercing it silently answered a different question than the caller
        # asked, and `rbs-strength` already refused it — one policy, not two.
        return f"{name!r} must be an integer, not a boolean"
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not _m.isfinite(value):
            return f"{name!r} must be a finite number"
        if value != int(value):
            # Truncating 1.9 to 1 loses 0.9 bp and returns 200 OK, so a model
            # that emits a fractional coordinate gets a wrong-by-one answer
            # with no indication anything was adjusted.
            return f"{name!r} must be a whole number, not {value!r}"
        return int(value)
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            return f"{name!r} must be an integer"
    return f"{name!r} must be an integer"


def _agent_record_is_circular(app) -> bool:
    """Topology of the record currently on the canvas.

    Only an EXPLICIT `linear` annotation means linear, matching
    `splicecraft_util._record_is_circular` and the app's own reading."""
    rec = getattr(app, "_current_record", None)
    return str(
        (getattr(rec, "annotations", {}) or {}).get("topology", "")
    ).strip().lower() != "linear"


def _sanitize_bases(s: "str | None", *,
                     max_len: int = 1_000_000) -> "tuple[str, str | None]":
    """Validate IUPAC DNA. Returns ``(cleaned, error_msg)``;
    ``error_msg`` is None on success. Caps at 1 MB to keep
    `_rebuild_record_with_edit` from hanging on adversarial input —
    matches the agent-API body cap so the whole request would have
    been rejected upstream anyway, but the cap-check here lets non-
    HTTP callers (modal `_gather`) reuse the same bound."""
    if s is None:
        return "", "missing 'bases'"
    s = str(s).strip().upper()
    if len(s) > max_len:
        return s, f"'bases' too long ({len(s)} > {max_len})"
    invalid = [c for c in s if c not in "ACGTNRYWSMKBDHV"]
    if invalid:
        return s, f"non-IUPAC characters: {''.join(sorted(set(invalid)))!r}"
    return s, None


def _agent_endpoint(name: str, *, write: bool = False):
    """Decorator: register a handler at `/<name>`.
    Handlers take `(app, payload)` and return either a `dict` (200) or
    `(dict, status_code)`. `write=True` flags state-mutating endpoints —
    these require the bearer token AND refuse if `app._unsaved` is True
    (unless the payload has `{"force": true}`)."""
    def deco(fn):
        _state._AGENT_HANDLERS[name] = (fn, write)
        return fn
    return deco


def _agent_save_or_500(save_fn, label: str = "library"):
    """Run `save_fn()`; on OSError/RuntimeError, return
    ``({"error": "save failed: ..."}, 500)`` so the handler can
    propagate. Returns None on success — caller pattern is
    ``err = _agent_save_or_500(...); if err: return err``.

    Per sacred invariant #7, `_safe_save_json` re-raises on disk
    failure so callers can surface; agent endpoints used to let those
    exceptions bubble up as generic 500s with no actionable detail.
    Routes through `_notify_save_failure` so the in-process UI user
    also sees the failure toast (the agent path and the GUI share an
    `app` object).
    """
    try:
        save_fn()
    except (OSError, RuntimeError) as exc:
        _notify_save_failure(_state._LIVE_APP_REF.get(), label, exc)
        return ({"error": f"save failed for {label}: {_scrub_path(str(exc))}"}, 500)
    return None


def _agent_dirty_guard(app, payload):
    """Return None if writes may proceed, else (error_dict, 409). The
    ``force`` field in the POST JSON body overrides the dirty check.
    (Writes are POST-only and the dispatcher parses only the JSON body
    for POST, so a ``?force=1`` query param is NOT honoured — pass
    ``{"force": true}`` in the body.)"""
    if getattr(app, "_unsaved", False) and not bool(payload.get("force")):
        return ({"error":
                  "unsaved changes — pass {\"force\": true} to override",
                  "dirty": True}, 409)
    return None


def _agent_ignored_keys(payload, known):
    """Body keys the handler does NOT act on — for the response's
    ``ignored`` echo (snag #14). Silent param-dropping hid real caller
    bugs: a ``save {name: …}`` that renamed nothing still returned
    ``ok: true``, so the caller believed a rename happened. Surfacing the
    unconsumed keys makes the no-op visible WITHOUT hard-rejecting
    forward-compatible callers (a 400 on an unknown key would break an
    agent written against a newer schema). ``force`` is always allowed
    (the universal dirty-guard override, applied by the dispatcher)."""
    if not isinstance(payload, dict):
        return []
    allowed = set(known) | {"force"}
    return sorted(str(k) for k in payload if k not in allowed)


# Routing / selection params that change WHERE or HOW an operation applies.
# Silently dropping one is the SC-D footgun (a call that "succeeds" while
# routing somewhere else). This is a CLOSED vocabulary of current API params —
# NOT a catch-all for unknown keys, so a newer-schema caller passing a brand-new
# key still gets the forward-compat leniency `_agent_ignored_keys` documents.
# Exact-match (case-sensitive) on the CURRENT param names — add any NEW
# routing/selection param here as it's introduced so the guard keeps covering it.
_AGENT_DANGEROUS_PARAMS = frozenset({
    "collection", "source_collection", "bin", "parts_bin",
    "enzyme", "enzymes", "orientation", "rename",
    # `carry_annotations` decides whether the saved product keeps its parents'
    # features. Only `traditional-clone` implements it; silently ignoring it
    # elsewhere hands back a feature-BARE construct to a caller who asked for
    # an annotated one and has no way to tell (blue field report #7).
    "carry_annotations",
})


def _agent_reject_dangerous_unknowns(payload, allowed):
    """Fail loud (400) when ``payload`` carries a routing/selection param from
    `_AGENT_DANGEROUS_PARAMS` that this endpoint does NOT accept — the "silent
    ok:true while doing something else" footgun (agent-API feedback SC-D). Any
    OTHER unknown key stays soft (forward-compat — same policy as
    `_agent_ignored_keys`). ``allowed`` is the endpoint's recognised key-set
    (the same one it passes to `_agent_ignored_keys`). Returns ``(err, 400)``
    or None; call it EARLY, before any state change."""
    if not isinstance(payload, dict):
        return None
    ok = set(allowed) | {"force"}
    bad = sorted(str(k) for k in payload
                 if k in _AGENT_DANGEROUS_PARAMS and k not in ok)
    if bad:
        keys = ", ".join(repr(k) for k in bad)
        return ({"error":
                  f"unrecognized routing/selection parameter(s) {keys} for "
                  "this endpoint — they would change WHERE or HOW the operation "
                  "applies, but this endpoint doesn't accept them (check the "
                  "spelling / the endpoint). Nothing was changed."}, 400)
    return None


def _agent_sanitize_qualifiers(raw):
    """Validate + sanitise a feature ``qualifiers`` object from an agent
    payload. Returns ``(cleaned_dict, error_msg)`` — `error_msg` is None on
    success, in which case `cleaned_dict` is safe to merge into a SeqFeature.

    ``None`` input → ``({}, None)`` (qualifiers are optional). Keys are
    sanitised labels; each value is coerced to a LIST of sanitised strings (a
    bare string becomes a 1-element list, matching GenBank's repeated-
    qualifier model). Caps at 64 keys × 32 values/key so a hostile or
    accidental blob can't bloat the record / `.gb` export. Lets an agent set
    arbitrary GenBank qualifiers (``gene``, ``product``, ``note``, …) the way
    the AddFeatureModal's fields + colour picker do."""
    if raw is None:
        return {}, None
    if not isinstance(raw, dict):
        return None, "'qualifiers' must be an object (key -> string/list)"
    if len(raw) > 64:
        return None, "too many qualifiers (max 64)"
    out: dict = {}
    for k, v in raw.items():
        key = _sanitize_label(str(k), max_len=50)
        if not key:
            continue
        if isinstance(v, (list, tuple)):
            vals = [_sanitize_note(str(x), max_len=2000) for x in list(v)[:32]]
        else:
            vals = [_sanitize_note(str(v), max_len=2000)]
        out[key] = [x for x in vals if x]
    return out, None


@_agent_endpoint("get-sequence")
def _h_get_sequence(app, payload):
    """Return DNA from `[start, end)`. Body: ``{start, end, bottom?}``.
    `bottom: true` returns the reverse-complement (5'→3' on the
    bottom strand). Wrap-aware: `end < start` is interpreted as a
    span that wraps the origin."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    seq = str(rec.seq).upper()
    n   = len(seq)
    try:
        if "start" not in payload or "end" not in payload:
            raise KeyError("start/end")
        start = _coerce_int(payload["start"], name="start")
        end = _coerce_int(payload["end"], name="end")
    except (KeyError, ValueError, TypeError, OverflowError):
        return ({"error": "missing or invalid 'start'/'end'"}, 400)
    if isinstance(start, str):
        return ({"error": start}, 400)
    if isinstance(end, str):
        return ({"error": end}, 400)
    if not (0 <= start <= n) or not (0 <= end <= n):
        return ({"error": f"start/end out of range [0, {n}]"}, 400)
    if end >= start:
        sub = seq[start:end]
    elif not _agent_record_is_circular(app):
        # `end < start` means "wrap the origin", and a LINEAR molecule has no
        # origin to wrap. Answering it anyway returned bases that do not exist
        # as one contiguous span on that DNA.
        return ({"error": "end < start means an origin wrap, but this record "
                          "is linear — pass end >= start"}, 400)
    else:
        sub = seq[start:] + seq[:end]
    if bool(payload.get("bottom")):
        sub = _rc(sub)
    return {
        "ok":     True,
        "start":  start,
        "end":    end,
        "bottom": bool(payload.get("bottom")),
        "length": len(sub),
        "seq":    sub,
    }


@_agent_endpoint("fold-rna")
def _h_fold_rna(app, payload):
    """Fold a sequence to its minimum-free-energy RNA secondary structure.
    Body: ``{sequence}`` (alias ``{seq}``). DNA ``T`` is read as RNA
    ``U``. Returns ``{ok, structure, dg, length}`` — `structure` is the
    dot-bracket MFE secondary structure and `dg` the free energy in
    kcal/mol. Pure-Python Turner-2004 nearest-neighbor model (no external
    dependency); inputs must be unambiguous A/C/G/U(T) and within the
    folding length cap, else 400 with the reason."""
    seq = payload.get("sequence")
    if seq is None:
        seq = payload.get("seq")
    if not isinstance(seq, str):
        return ({"error": "missing or non-string 'sequence'"}, 400)
    try:
        structure, dg = _rna_fold(seq)
    except ValueError as exc:
        return ({"error": str(exc)}, 400)
    return {
        "ok":        True,
        "structure": structure,
        "dg":        round(dg, 2),
        "length":    len(structure),
    }


@_agent_endpoint("cofold-rna")
def _h_cofold_rna(app, payload):
    """Bound-state heterodimer ΔG of two strands. Body: ``{seq_a, seq_b}``.
    Returns ``{ok, dg}`` — the free energy (kcal/mol) of strand B bound to
    strand A (e.g. a 16S anti-SD : mRNA hybrid; DNA ``T`` read as ``U``).
    Pure-Python Turner-2004 cofold (no external dependency); ambiguous /
    over-length input → 400 with the reason."""
    a = payload.get("seq_a")
    b = payload.get("seq_b")
    if not isinstance(a, str) or not isinstance(b, str):
        return ({"error": "missing or non-string 'seq_a' / 'seq_b'"}, 400)
    try:
        dg = _rna_cofold(a, b)
    except ValueError as exc:
        return ({"error": str(exc)}, 400)
    return {"ok": True, "dg": round(dg, 2)}


@_agent_endpoint("rbs-strength")
def _h_rbs_strength(app, payload):
    """Relative E. coli translation-initiation strength of a ribosome
    binding site. Body: ``{mrna, start}`` — `start` is the 0-based index of
    the start codon (DNA ``T`` read as ``U``). Returns ``{ok, dg_total,
    dg_mrna, dg_hybrid, spacing, rel_strength}``. `rel_strength` is RELATIVE
    (only ratios between RBSs are meaningful — a biophysically-grounded
    ranking score, not an absolute rate). Bad input → 400."""
    mrna = payload.get("mrna") or payload.get("sequence")   # SC-G alias
    if not isinstance(mrna, str):
        return ({"error": "missing or non-string 'mrna'"}, 400)
    if isinstance(payload.get("start"), bool):    # JSON true/false ≠ an index
        return ({"error": "'start' must be an integer, not a boolean"}, 400)
    try:
        start = int(payload["start"])
    except (KeyError, ValueError, TypeError, OverflowError):
        return ({"error": "missing or invalid 'start' (0-based int)"}, 400)
    try:
        result = _rbs_strength(mrna, start)
    except ValueError as exc:
        return ({"error": str(exc)}, 400)
    return {"ok": True, **result}


@_agent_endpoint("design-rbs")
def _h_design_rbs(app, payload):
    """Reverse-design a 5'UTR (Shine-Dalgarno + spacer) for a target
    relative RBS strength. Body: ``{cds, target}`` — `cds` begins with the
    start codon, `target` is a non-negative relative strength (on the
    `rbs-strength` scale); optional ``{upstream}`` (5' context). Returns
    ``{ok, utr, full, sd, spacing, rel_strength, dg_total,
    achievable_min, achievable_max, on_target}`` — the design closest to
    the target (nearest achievable + `on_target=false` if out of range).
    Runs ~80 fold/cofold evaluations, so it blocks for a few seconds. Bad
    input → 400."""
    cds = payload.get("cds") or payload.get("sequence")   # SC-G alias
    if not isinstance(cds, str):
        return ({"error": "missing or non-string 'cds'"}, 400)
    target = payload.get("target")
    if not isinstance(target, (int, float)) or isinstance(target, bool):
        return ({"error": "missing or invalid 'target' (non-negative number)"}, 400)
    kw = {}
    up = payload.get("upstream")
    if isinstance(up, str):
        kw["upstream"] = up
    try:
        result = _rbs_design(cds, target, **kw)
    except ValueError as exc:
        return ({"error": str(exc)}, 400)
    return {"ok": True, **result}


# Real bacterial operons carry a handful of genes; each gene triggers a
# multi-eval context-aware RBS design (~seconds of pure-Python folding), so an
# unbounded `genes` list is an authenticated-local CPU DoS (~20k minimal genes
# fit inside the 1 MiB body cap). Cap it like the sibling assembly endpoints
# (simulate-gibson ≤64, simulate-golden-gate ≤30); the endpoint is also
# heavy-classified (`_AGENT_HEAVY_ENDPOINTS`) so a burst can't pin every worker.
_ASSEMBLE_OPERON_MAX_GENES = 32


@_agent_endpoint("assemble-operon")
def _h_assemble_operon(app, payload):
    """Assemble a contiguous bacterial operon from CDSs, each RBS
    reverse-designed (context-aware) to its target relative strength.
    Body: ``{genes: [{cds, target_strength, name?}, ...]}`` + optional
    ``{promoter, terminator}`` (DNA). Returns ``{ok, sequence, layout,
    genes}`` — `sequence` is DNA; `layout` tiles it (promoter / rbs / cds /
    terminator); `genes` reports target / achieved rel_strength /
    on_target. A target unreachable in its assembled context yields the
    nearest achievable + on_target=false. Bad input → 400."""
    genes = payload.get("genes")
    if not isinstance(genes, list) or not genes:
        return ({"error": "missing or empty 'genes' list"}, 400)
    if len(genes) > _ASSEMBLE_OPERON_MAX_GENES:
        return ({"error": f"too many genes ({len(genes)}); cap is "
                 f"{_ASSEMBLE_OPERON_MAX_GENES}. Each gene runs a multi-eval "
                 f"RBS design — split the operon or assemble in stages."}, 400)
    kw = {}
    for key in ("promoter", "terminator"):
        val = payload.get(key)
        if isinstance(val, str):
            kw[key] = val
    try:
        result = _assemble_operon(genes, **kw)
    except (ValueError, TypeError) as exc:
        return ({"error": str(exc)}, 400)
    return {"ok": True, **result}


_EXPORT_DNA_EXTS     = (".dna",)


_EXPORT_EMBL_EXTS    = (".embl",)


@_agent_endpoint("export-commercialsaas", write=True)
def _h_export_commercialsaas(app, payload):
    """Write the current record to `path` as CommercialSaaS `.dna`.
    Body: ``{path}``. Requires the loaded record to be in the active
    library (`.dna` writer keys off the library entry's id / sidecar);
    returns 422 otherwise. Splice mode (preserves all original packets
    + updated history) when a sidecar exists; from-scratch mode
    (sequence + features + notes + history only — cosmetic packets
    omitted) otherwise. Mirrors `action_export_commercialsaas`."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    path = _sanitize_path(payload.get("path"))
    if path is None:
        return ({"error": "missing 'path'"}, 400)
    if (ext_err := _check_export_extension(
            path, _EXPORT_DNA_EXTS, "CommercialSaaS .dna")) is not None:
        return ({"error": ext_err}, 400)
    err = _check_agent_write_path(path)
    if err is not None:
        return ({"error": err}, 403)
    rec_id = getattr(rec, "id", "") or ""
    # Sweep #26: targeted clone via `_find_library_entry_by_id`.
    entry = _find_library_entry_by_id(rec_id) if rec_id else None
    if entry is None:
        return ({"error": (
            "loaded plasmid is not in the library; .dna export keys "
            "off the library entry's id / sidecar. Add to library "
            "first, or use export-genbank."
        )}, 422)
    try:
        out_path = _export_commercialsaas_dna(entry, path)
    except (OSError, ValueError) as exc:
        return ({"error": f"export failed: {_scrub_path(str(exc))}"}, 500)
    try:
        size = Path(out_path).stat().st_size
    except OSError:
        size = 0
    return {"ok": True, "path": str(out_path), "bytes": size}


@_agent_endpoint("export-embl", write=True)
def _h_export_embl(app, payload):
    """Write the current record to `path` as EMBL flatfile. Body: ``{path}``.
    EMBL is a near-equivalent of GenBank — same feature model, different
    text layout. Useful for tools that consume EMBL directly (ENA
    submissions, some annotation pipelines). Returns ``{ok, path, bp,
    features}``."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    path = _sanitize_path(payload.get("path"))
    if path is None:
        return ({"error": "missing 'path'"}, 400)
    if (ext_err := _check_export_extension(
            path, _EXPORT_EMBL_EXTS, "EMBL")) is not None:
        return ({"error": ext_err}, 400)
    err = _check_agent_write_path(path)
    if err is not None:
        return ({"error": err}, 403)
    try:
        result = _export_embl_to_path(rec, path)
    except (OSError, ValueError) as exc:
        return ({"error": f"export failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, **result}


@_agent_endpoint("export-map-image", write=True)
def _h_export_map_image(app, payload):
    """Render the LOADED plasmid's circular map to `path` as PNG or SVG. Body:
    ``{path, format?, size?, transparent?, labels?, sites?}``. ``format`` ∈
    {"png","svg"} (inferred from the path extension when omitted). Fills the
    single-plasmid gap that `bulk-export-collection` (whole collection) left.
    Path is sanitised + symlink-refused like the other exporters; the write is
    atomic. Returns ``{ok, path, fmt, bp, features, bytes}``."""
    import splicecraft_mapimage as _mi
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    path = _sanitize_path(payload.get("path"))
    if path is None:
        return ({"error": "missing 'path'"}, 400)
    fmt = str(payload.get("format") or path.suffix[1:] or "png").lower()
    if fmt not in _mi._MAP_EXPORT_FORMATS:
        return ({"error": f"unknown format {fmt!r}; "
                          f"choose one of {sorted(_mi._MAP_EXPORT_FORMATS)}"}, 400)
    if path.suffix.lower() != f".{fmt}":
        return ({"error": f"path extension must be .{fmt} for format {fmt!r}"},
                400)
    err = _check_agent_write_path(path)
    if err is not None:
        return ({"error": err}, 403)
    size = _coerce_int(payload.get("size", 1400), name="size")
    if isinstance(size, str):
        return ({"error": size}, 400)
    show_labels = bool(payload.get("labels", True))
    show_sites = bool(payload.get("sites", True))
    transparent = bool(payload.get("transparent", False))
    title = str(getattr(rec, "_tui_display_name", None)
                or getattr(rec, "name", "") or "")
    try:
        sites = _mi._map_sites_from_record(rec) if show_sites else []
        result = _mi.export_plasmid_map(
            path, record=rec, fmt=fmt, title=title, sites=sites,
            size=int(size), transparent=transparent,
            show_labels=show_labels, show_sites=show_sites)
    except (OSError, ValueError) as exc:
        return ({"error": f"map export failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, **result}


@_agent_endpoint("list-library")
def _h_list_library(app, payload):
    """Plasmid library entries. Returns name, id, length (bp),
    n_features, topology, source for each. Subset of the on-disk JSON
    tuned for an agent's "what plasmids do I have" question.

    Body: ``{collection?}``. Pass ``collection`` to scope the listing to
    ONE collection's plasmids (``404`` if no such collection, ``400`` if
    blank); omit it for the full library across all collections. P1-7
    (agent-API feedback): pre-fix a ``collection`` arg was silently ignored
    and the full library returned — the SC-D footgun where the result
    looks scoped but isn't. The response echoes the resolved
    ``collection`` when scoped.

    Sweep #26 (2026-05-25): fixed projection bug. Pre-fix the handler
    read non-existent JSON fields (`sequence`, `features`, `topology`,
    `source_path`) and returned `{length: 0, n_features: 0, topology:
    "", source_path: ""}` for every entry — the on-disk schema uses
    `size`, `n_feats`, `source` and topology lives in the GenBank
    LOCUS line. Now projects the real fields and pulls topology out
    of the LOCUS line when available (defaults to "circular" — the
    overwhelming majority of library entries).
    """
    # P1-7: scope to a single collection when asked. Resolve the name the
    # same way `load-entry` does (`_normalize_collection_name` + exact
    # match on the readonly cache view, mirroring `_find_collection`) so
    # `{collection: "RUBYFIRE"}` behaves identically across endpoints.
    coll_in = payload.get("collection")
    scope: "str | None" = None
    if coll_in not in (None, ""):
        coll_name = _normalize_collection_name(coll_in)
        if coll_name is None:
            return ({"error": "invalid 'collection'"}, 400)
        col = next((c for c in _iter_collections_readonly()
                    if isinstance(c, dict) and c.get("name") == coll_name),
                   None)
        if col is None:
            return ({"error": f"no collection named {coll_name!r}"}, 404)
        # Read-only projection of the collection's own plasmid entries —
        # never mutate them (pitfall #17).
        source_entries: "list[dict]" = [
            e for e in (col.get("plasmids") or []) if isinstance(e, dict)
        ]
        scope = coll_name
    else:
        # Sweep #25 (2026-05-23): read-only iterator — pre-fix paid full
        # library deepcopy just to project six fields. Build the output
        # list from immutable reads.
        source_entries = _iter_library_readonly()

    out: list[dict] = []
    for e in source_entries:
        # Topology from the LOCUS line's topology FIELD (whole token, that
        # line only). Library entries omit explicit topology in the JSON, so
        # it comes from the embedded gb_text; default "circular" because
        # ~99% of stored entries are plasmids. This used to be
        # `"linear" in gb_text[:200].lower()`, which also read the entry's
        # own NAME and its DEFINITION — an entry called `pLinear2` reported
        # itself linear to every agent caller.
        topology = _topology_from_gb_text(e.get("gb_text") or "")
        out.append({
            "name":        e.get("name", ""),
            "id":          e.get("id", ""),
            "length":      int(e.get("size") or 0),
            "n_features":  int(e.get("n_feats") or 0),
            "topology":    topology,
            "source":      e.get("source", "") or "",
        })
    result = {"library": out, "count": len(out)}
    if scope is not None:
        result["collection"] = scope
    return result


@_agent_endpoint("list-collections")
def _h_list_collections(app, payload):
    """Plasmid collection buckets. Returns name + plasmid count for
    each collection, and which one is currently active."""
    # Sweep #25 (2026-05-23): read-only iterator — pre-fix paid full
    # ~160 MB collections.json deepcopy just to project name +
    # n_plasmids per collection.
    cols = _iter_collections_readonly()
    active = _get_active_collection_name()
    return {
        "active": active,
        "collections": [
            {"name": c.get("name", ""),
             "n_plasmids": len(c.get("plasmids", []) or [])}
            for c in cols
        ],
    }


@_agent_endpoint("search-library")
def _h_search_library(app, payload):
    """Cross-collection plasmid name search. Body: ``{query?: str,
    limit?: int}``. Empty `query` returns the first `limit` plasmids
    across all collections (default `limit` = 200, capped at 1000).
    Each match is ``{collection, name, id, size, status, n_feats}``.

    Use `set-active-collection` + `load-entry` (or the GUI
    equivalent) to actually load a hit. The fuzzy matcher is the
    same one the library panel's in-collection search uses, so a
    query of `puc` matches `pUC19`."""
    query = payload.get("query")
    if query is not None and not isinstance(query, str):
        return ({"error": "'query' must be a string"}, 400)
    # Cap query length to 200 chars so an agent's runaway prompt
    # can't pin the fuzzy matcher (audit hardening 2026-05-14).
    if query is not None and len(query) > 200:
        return ({"error": "'query' exceeds 200-char cap"}, 400)
    raw_limit = payload.get("limit", 200)
    limit = _coerce_int(raw_limit, name="limit")
    if isinstance(limit, str):
        return ({"error": limit}, 400)
    limit = max(1, min(1000, limit))
    matches = _search_collections_library(query or "", limit=limit)
    return {"matches": matches, "count": len(matches)}


@_agent_endpoint("list-plasmidsaurus-members")
def _h_list_plasmidsaurus_members(app, payload):
    """List the GenBank-format members of a Plasmidsaurus result
    zip. Body: ``{path}``. Returns ``{members: [{name, size}, ...]}``.

    Cap-protected — zips above `_PLASMIDSAURUS_ZIP_MAX_BYTES`
    (500 MB), members above `_PLASMIDSAURUS_MEMBER_MAX_BYTES` and
    listings beyond `_PLASMIDSAURUS_MAX_MEMBERS` are refused or
    skipped to keep the picker snappy and resistant to malformed
    archives. Symlinks at the path are rejected via
    `_safe_file_size_check` upstream.

    Read-only — does not extract or align anything; the agent uses
    the returned member name with `align-plasmidsaurus-zip` to run
    the actual alignment.
    """
    raw_path = payload.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        return ({"error": "missing 'path'"}, 400)
    path = _sanitize_path(raw_path)
    if path is None:
        return ({"error": "could not sanitize 'path'"}, 400)
    # Sweep #26 (2026-05-23): defense-in-depth ancestor symlink walk.
    # Pre-sweep an attacker-placed symlink at any parent (e.g.
    # `~/Documents` → `/etc`) was silently followed by the
    # downstream `os.open`.
    anc_err = _check_agent_read_path_ancestors(path)
    if anc_err is not None:
        _log.warning("agent list-plasmidsaurus-members: %s", anc_err)
        return ({"error": "zip rejected (see splicecraft log)"}, 400)
    # Sweep #25 (2026-05-23): collapse path-shape errors to a uniform
    # 400 so the differentiated error responses don't act as a
    # filesystem-state oracle (pre-fix an unauthenticated caller
    # could probe arbitrary paths via "not found" / "could not open
    # zip" / "zip too large" responses). Logs still carry detail.
    ok, reason = _safe_file_size_check(
        path, _PLASMIDSAURUS_ZIP_MAX_BYTES, "Plasmidsaurus zip",
    )
    if not ok:
        _log.warning("agent list-plasmidsaurus-members: rejected (%s): %s",
                     reason or "unsafe", path)
        return ({"error": "zip rejected (see splicecraft log)"}, 400)
    try:
        members = _list_gbk_members_in_zip(path)
    except (ValueError, OSError) as exc:
        _log.warning(
            "agent list-plasmidsaurus-members: walk failed (%s)", exc,
        )
        return ({"error": "zip rejected (see splicecraft log)"}, 400)
    return {"ok": True, "path": str(path),
            "members": members, "count": len(members)}


@_agent_endpoint("lint-synthesis")
def _h_lint_synthesis(app, payload):
    """Synthesis / assembly-readiness pre-flight for a construct — the "is it
    safe to order + assemble?" check. Body:
    ``{sequence | id | name, circular?=true, expect_cds?=false,
        forbidden_enzymes?: [str]}``.

    Resolves the sequence from a bare ``sequence`` OR a saved library entry
    (``id``/``name`` — its stored topology sets ``circular`` unless you pass
    ``circular`` explicitly). Scans for internal Type IIS sites (Golden Gate /
    MoClo assembly killers → errors), extreme overall + windowed GC, long
    homopolymer runs, tandem repeats, degenerate bases, and — with
    ``expect_cds`` — a clean full-length ORF. Returns ``{ok, score, warnings,
    stats}`` where ``score`` is 0-100 (100 = clean) and each ``warnings`` entry
    is ``{level, kind, message, start?, end?}``. Read-only."""
    circular = bool(payload.get("circular", True))
    raw_seq = payload.get("sequence")
    if isinstance(raw_seq, str) and raw_seq.strip():
        seq, e = _sanitize_bases(raw_seq)
        if e:
            return ({"error": f"'sequence' rejected: {e}"}, 400)
    else:
        key = _sanitize_label(payload.get("id") or payload.get("name"),
                              max_len=200)
        if not key:
            return ({"error": "provide 'sequence', or 'id'/'name' of a saved "
                              "library entry to lint"}, 400)
        entry = (_find_library_entry_by_id(key)
                 or _find_library_entry_by_name(key))
        if entry is None:
            return ({"error": f"no library entry named {key!r}"}, 404)
        try:
            rec = _gb_text_to_record(entry.get("gb_text") or "")
        except Exception as exc:
            return ({"error":
                      f"couldn't parse {key!r}: {_scrub_path(str(exc))}"}, 500)
        seq = str(getattr(rec, "seq", "") or "")
        if "circular" not in payload:
            circular = (rec.annotations.get("topology") or "").lower() != "linear"
    if not seq:
        return ({"error": "sequence is empty"}, 422)
    if len(seq) > _PAIRWISE_MAX_LEN:
        return ({"error": f"sequence exceeds {_PAIRWISE_MAX_LEN:,} bp cap"}, 413)
    kwargs: dict = {}
    fenz = payload.get("forbidden_enzymes")
    if fenz is not None:
        if not isinstance(fenz, list) or not all(isinstance(x, str) and x
                                                  for x in fenz):
            return ({"error":
                      "'forbidden_enzymes' must be a list of enzyme names"}, 400)
        if len(fenz) > 32:
            return ({"error": "'forbidden_enzymes' too long (max 32)"}, 400)
        # Fail loud on an unknown/typo'd enzyme: a synthesis-SAFETY check must
        # never silently skip a site it was asked to scan for (a clean score
        # would be false confidence).
        catalog = _state._all_enzymes_hook()
        unknown = [x for x in fenz if x not in catalog]
        if unknown:
            return ({"error":
                      f"unknown enzyme(s) in 'forbidden_enzymes': {unknown[:8]}"},
                    400)
        kwargs["forbidden_enzymes"] = tuple(fenz)
    try:
        result = _synthesis_lint(seq, circular=circular,
                                 expect_cds=bool(payload.get("expect_cds")),
                                 **kwargs)
    except Exception as exc:
        _log.exception("agent lint-synthesis: lint failed")
        return ({"error": f"lint failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, **result,
            "ignored": _agent_ignored_keys(payload, {
                "sequence", "id", "name", "circular", "expect_cds",
                "forbidden_enzymes"})}


@_agent_endpoint("find-orfs")
def _h_find_orfs(app, payload):
    """Six-frame ORF scan over the loaded record. Body:
    ``{min_aa?: int, include_alt_starts?: bool}``. Defaults: 30 aa,
    ATG only. The length cutoff is in AMINO ACIDS (`min_aa`); `min_length` /
    `min_bp` are rejected (400) so a bp-vs-aa unit mix-up can't silently return
    a default-length result. Returns
    ``{orfs: [{start, end, strand, length_aa, aa_seq}, ...]}``
    sorted by length descending. ORFs that cross the origin on a
    circular plasmid are reported with `end < start` — the same
    wrap-feature convention used elsewhere."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    # find-orfs filters by AMINO-ACID length (`min_aa`). `min_length` / `min_bp`
    # are the natural-but-wrong guesses — `list-restriction-sites` uses
    # `min_length` in BASE PAIRS, so an agent reusing that key here would have
    # it silently dropped (read endpoints don't echo `ignored`) and get a
    # default-30 ORF set that LOOKS authoritative — the SC-D / SC-G footgun the
    # agent-API feedback flags as most dangerous. Fail loud so the caller fixes
    # both the key AND the aa-vs-bp unit, rather than trusting a wrong answer.
    for _wrong in ("min_length", "min_bp", "min_len"):
        if _wrong in payload:
            return ({"error":
                      f"unknown parameter {_wrong!r}; find-orfs filters by "
                      f"amino-acid length — use 'min_aa' (amino acids, "
                      f"default 30). Note: ORF length here is in aa, not bp."},
                    400)
    min_aa = _coerce_int(payload.get("min_aa", 30), name="min_aa")
    if isinstance(min_aa, str):
        return ({"error": min_aa}, 400)
    if min_aa < 1:
        return ({"error": "'min_aa' must be ≥ 1"}, 400)
    alt = bool(payload.get("include_alt_starts", False))
    seq = str(rec.seq) if rec.seq is not None else ""
    if not seq:
        return {"orfs": [], "count": 0}
    annotations = getattr(rec, "annotations", None) or {}
    is_circular = annotations.get("topology") == "circular"
    orfs = _find_orfs(
        seq,
        min_aa=min_aa,
        include_alt_starts=alt,
        circular=is_circular,
    )
    return {"orfs": orfs, "count": len(orfs)}


@_agent_endpoint("list-cas-variants")
def _h_list_cas_variants(app, payload):
    """The Cas nucleases ``design-guides`` supports. Returns
    ``{variants: [{key, label, pam, pam_side, guide_len, cut_offsets, notes},
    ...]}``.

    Three distinct PAM GEOMETRIES rather than a long catalogue of near
    duplicates — a short 3' PAM, a long degenerate 3' PAM, and a 5' PAM with a
    staggered cut. Use `key` as the ``variant`` argument elsewhere."""
    return {"variants": [
        {"key": k, "label": v["label"], "pam": v["pam"],
         "pam_side": v["pam_side"], "guide_len": v["guide_len"],
         "cut_offsets": list(v["cut_offsets"]), "notes": v["notes"]}
        for k, v in _CAS_VARIANTS.items()
    ]}


def _crispr_target_sequence(app, payload):
    """Resolve the sequence to scan: an explicit ``sequence``, or the loaded
    plasmid when none is given.

    Returns ``(seq, is_circular, err_tuple)``. A raw ``sequence`` is treated as
    LINEAR unless ``circular`` says otherwise, because a pasted exon is the
    common case and silently wrapping it would invent guides across a junction
    that does not exist.
    """
    raw = payload.get("sequence")
    if raw is not None:
        seq, err = _sanitize_bases(raw)
        if err:
            return ("", False, ({"error": err}, 400))
        if not seq:
            return ("", False, ({"error": "'sequence' is empty"}, 400))
        return (seq, bool(payload.get("circular", False)), None)
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ("", False, ({"error": "no plasmid loaded and no 'sequence' "
                                      "given"}, 422))
    seq = str(rec.seq) if rec.seq is not None else ""
    if not seq:
        return ("", False, ({"error": "loaded plasmid has no sequence"}, 422))
    anns = getattr(rec, "annotations", None) or {}
    circ = (anns.get("topology") == "circular")
    return (seq, bool(payload.get("circular", circ)), None)


@_agent_endpoint("design-guides")
def _h_design_guides(app, payload):
    """CRISPR guide design over the loaded plasmid or a supplied sequence.

    Body: ``{sequence?, circular?, variant?: "spcas9"|"sacas9"|"lbcas12a",
    region?: [start, end], u6_driven?: bool, offtarget?: bool|str,
    max_mismatch?: int, limit?: int}``.

    Returns ``{variant, label, circular, n_found, guides: [...], truncated,
    offtarget_searched_bp}``. Each guide carries its forward-frame coordinates,
    strand, cut site, a `score` block of NAMED flags (not a fabricated
    efficiency number) and, when off-target search ran, `offtargets`.

    **`offtarget_searched_bp` is the honest scope of any off-target claim.**
    `offtarget: true` searches the target sequence itself; pass a SEQUENCE
    string to search something else (a host genome you hold). Omit it and no
    off-target claim is made at all — `offtarget_searched_bp` is then 0, so a
    caller cannot read silence as "clean". There is no genome-wide prediction
    here and none is implied; see `splicecraft_crispr`'s module docstring.
    """
    seq, circ, err = _crispr_target_sequence(app, payload)
    if err:
        return err
    variant = payload.get("variant", "spcas9")
    try:
        _cas_variant(variant)
    except ValueError as exc:
        return ({"error": str(exc)}, 400)

    region = payload.get("region")
    if region is not None:
        if (not isinstance(region, (list, tuple))) or len(region) != 2:
            return ({"error": "'region' must be [start, end] (0-based, "
                              "half-open)"}, 400)
        raw_lo = _coerce_int(region[0], name="region[0]")
        if isinstance(raw_lo, str):
            return ({"error": raw_lo}, 400)
        raw_hi = _coerce_int(region[1], name="region[1]")
        if isinstance(raw_hi, str):
            return ({"error": raw_hi}, 400)
        lo, hi = int(raw_lo), int(raw_hi)
        if lo < 0 or hi < 0 or lo > len(seq) or hi > len(seq):
            return ({"error": f"'region' out of range for a {len(seq)} bp "
                              f"sequence"}, 400)
        if lo == hi:
            return ({"error": "'region' is empty (start == end)"}, 400)
        region = (lo, hi)

    limit = _coerce_int(payload.get("limit", 50), name="limit")
    if isinstance(limit, str):
        return ({"error": limit}, 400)
    if limit < 1:
        return ({"error": "'limit' must be >= 1"}, 400)
    mm = _coerce_int(payload.get("max_mismatch", 4), name="max_mismatch")
    if isinstance(mm, str):
        return ({"error": mm}, 400)
    if not (0 <= mm <= 10):
        return ({"error": "'max_mismatch' must be 0-10"}, 400)

    ot = payload.get("offtarget")
    ot_seq = None
    if isinstance(ot, str) and ot.strip():
        cleaned, ot_err = _sanitize_bases(ot)
        if ot_err:
            return ({"error": f"'offtarget': {ot_err}"}, 400)
        ot_seq = cleaned
    elif ot is True:
        ot_seq = seq
    try:
        return _design_guides(
            seq, variant=variant, circular=circ,
            u6_driven=bool(payload.get("u6_driven", True)),
            region=region, offtarget_in=ot_seq, max_mismatch=mm,
            limit=limit,
        )
    except ValueError as exc:
        return ({"error": str(exc)}, 400)


@_agent_endpoint("score-guide")
def _h_score_guide(app, payload):
    """Triage a spacer you already have. Body:
    ``{guide: str, u6_driven?: bool}``.

    Returns ``{guide, gc_pct, homopolymer, self_complement, starts_with_g,
    flags, tier}``. `flags` names the mechanism behind each concern (a Pol III
    terminator truncates the transcript; a homopolymer slips in synthesis) —
    deliberately NOT an efficiency score, because a number that looks like one
    and is not would be worse than no number."""
    guide, err = _sanitize_bases(payload.get("guide"), max_len=200)
    if err:
        return ({"error": f"'guide': {err}"}, 400)
    if not guide:
        return ({"error": "'guide' is required"}, 400)
    out = _score_guide(guide, u6_driven=bool(payload.get("u6_driven", True)))
    out["guide"] = guide
    return out


@_agent_endpoint("guide-offtargets")
def _h_guide_offtargets(app, payload):
    """Mismatch-tolerant search for a spacer in a sequence. Body:
    ``{guide: str, sequence?, circular?, variant?, max_mismatch?: int,
    require_pam?: bool}``; `sequence` defaults to the loaded plasmid.

    Returns ``{searched_bp, max_mismatch, require_pam, hits, truncated}``.
    `searched_bp` bounds the claim — this searches WHAT YOU GAVE IT, never a
    genome. `seed_mismatches` counts mismatches in the PAM-proximal bases,
    where Cas9 tolerates them least, and hits are ordered most-dangerous
    first."""
    guide, err = _sanitize_bases(payload.get("guide"), max_len=200)
    if err:
        return ({"error": f"'guide': {err}"}, 400)
    if not guide:
        return ({"error": "'guide' is required"}, 400)
    seq, circ, serr = _crispr_target_sequence(app, payload)
    if serr:
        return serr
    variant = payload.get("variant", "spcas9")
    try:
        _cas_variant(variant)
    except ValueError as exc:
        return ({"error": str(exc)}, 400)
    mm = _coerce_int(payload.get("max_mismatch", 4), name="max_mismatch")
    if isinstance(mm, str):
        return ({"error": mm}, 400)
    if not (0 <= mm <= 10):
        return ({"error": "'max_mismatch' must be 0-10"}, 400)
    return _guide_offtargets(
        guide, seq, variant=variant, circular=circ, max_mismatch=mm,
        require_pam=bool(payload.get("require_pam", True)),
    )


@_agent_endpoint("guide-cloning-oligos")
def _h_guide_cloning_oligos(app, payload):
    """The annealed oligo pair that clones a spacer into a Type IIS guide
    vector. Body: ``{guide: str, top_prefix?, bottom_prefix?, add_g?: bool}``.

    Returns ``{top, bottom, duplex, guide, added_g, overhang_top,
    overhang_bottom}``. Body may override ``overhang_top`` / ``overhang_bottom``
    (defaults are the BbsI pair from the pX330 / pSpCas9(BB) family most lab
    guide vectors descend from; a BsmBI lentiCRISPR backbone wants different
    ones).

    The bottom oligo carries the complement of the transcription-start G as
    well as of the spacer, so both oligos come out the same length and anneal
    with no gap — writing it as `AAAC + rc(guide)` is the classic error.
    `added_g` reports whether that G was added (a guide already starting with G
    supplies its own)."""
    guide, err = _sanitize_bases(payload.get("guide"), max_len=200)
    if err:
        return ({"error": f"'guide': {err}"}, 400)
    if not guide:
        return ({"error": "'guide' is required"}, 400)
    kw = {}
    for key in ("overhang_top", "overhang_bottom"):
        if key in payload:
            val, verr = _sanitize_bases(payload.get(key), max_len=40)
            if verr:
                return ({"error": f"{key!r}: {verr}"}, 400)
            kw[key] = val
    try:
        return _guide_cloning_oligos(
            guide, add_g=bool(payload.get("add_g", True)), **kw)
    except ValueError as exc:
        return ({"error": str(exc)}, 400)


@_agent_endpoint("analyse-cds")
def _h_analyse_cds(app, payload):
    """Read the codon usage of an EXISTING CDS — the read half of
    ``optimize-protein``, which until now only wrote.

    Body: ``{sequence | id | name, feature?, taxid?, transl_table?,
    rare_w?, ramp_codons?}``. Pass a bare coding ``sequence``, or name a saved
    plasmid with ``id`` / ``name`` plus the CDS ``feature`` label (omit
    ``feature`` and the record's only CDS is used). ``taxid`` picks the codon
    table; the active table is used otherwise, falling back to E. coli K12.

    Returns ``{ok, source, cds_label, length_bp, n_codons, gc, gc3, cai,
    rare_pct, n_rare, tandem_rare_runs, worst_window, ramp, table, warnings}``.

      * ``cai`` is computed by the SAME function ``optimize-protein`` reports,
        against the same table, so the two can be compared directly.
      * ``rare_pct`` counts codons whose relative adaptiveness is below
        ``rare_w`` (default 0.2, echoed in the response — it is a convention,
        not a derived constant).
      * ``tandem_rare_runs`` lists consecutive rare codons, which stall a
        ribosome far more than the same codons scattered; a run is where the
        cost actually lands, and a percentage alone hides them.
      * ``worst_window`` is the lowest-CAI stretch, so a locally terrible
        region inside a healthy average is visible.
      * ``ramp`` reports the 5' end separately — a deliberately slow
        translational ramp there is normal, and averaging it into the whole
        gene flags healthy designs and hides unhealthy bodies.

    Read-only. Nothing about the sequence is changed or saved."""
    warnings: "list[str]" = []
    raw_seq = payload.get("sequence")
    cds_label = None
    if isinstance(raw_seq, str) and raw_seq.strip():
        dna, err = _sanitize_bases(raw_seq)
        if err:
            return ({"error": f"'sequence' rejected: {err}"}, 400)
        source = "sequence"
    else:
        key = _sanitize_label(payload.get("id") or payload.get("name"),
                              max_len=200)
        if not key:
            return ({"error": "provide 'sequence' (a CDS) or 'id' / 'name' "
                              "naming a saved plasmid"}, 400)
        entry = (_find_library_entry_by_id(key)
                 or _find_library_entry_by_name(key))
        if entry is None:
            return ({"error": f"no library entry named {key!r}"}, 404)
        try:
            rec = _gb_text_to_record(entry.get("gb_text") or "")
        except Exception as exc:
            return ({"error": f"couldn't parse {key!r}: "
                              f"{_scrub_path(str(exc))}"}, 500)
        seq = str(getattr(rec, "seq", "") or "")
        total = len(seq)
        circular = _record_is_circular(rec)
        want = _sanitize_label(payload.get("feature"), max_len=200)
        cands = []
        for feat in (getattr(rec, "features", []) or []):
            if str(getattr(feat, "type", "")).upper() != "CDS":
                continue
            lbl = _feat_label(feat) or ""
            if want and want.lower() not in lbl.lower():
                continue
            b = _feat_bounds(feat, total, circular=circular)
            if b is None:
                continue
            cands.append((lbl, b, feat))
        if not cands:
            return ({"error": (f"no CDS feature matching {want!r} on {key!r}"
                               if want else
                               f"no CDS feature on {key!r}")}, 404)
        if len(cands) > 1:
            return ({"error": f"{len(cands)} CDS features match — pass "
                              f"'feature' to pick one: "
                              f"{sorted(c[0] or '(unlabelled)' for c in cands)}"},
                    400)
        cds_label, (fs, fe, fstrand), cds_feat = cands[0]
        # ── Exons, not the outer span ────────────────────────────────
        # A SPLICED CDS read as one contiguous block runs straight through
        # its introns: every codon after the first intron is out of frame and
        # the GC3 / CAI / rare-codon numbers are then confidently wrong for a
        # gene the caller believes was measured. Reconstruct the coding
        # sequence from the exons instead.
        #
        # The wrap-vs-splice test is the SAME one `_agent_transcript_feature_dicts`
        # uses — an origin-WRAPPING feature is also a two-part compound
        # location (sacred invariant #9), and reading its parts as exons would
        # invent an intron spanning the backbone. `e >= s` means `_feat_bounds`
        # did not already resolve the location as a wrap.
        exons: "list[tuple[int, int]]" = []
        if fe >= fs:
            for part in (getattr(getattr(cds_feat, "location", None),
                                 "parts", None) or []):
                try:
                    exons.append((int(part.start), int(part.end)))
                except (TypeError, ValueError, AttributeError):
                    exons = []
                    break
        if len(exons) > 1:
            # Plus-strand order, then reverse-complement as a whole for a
            # minus-strand gene — which is what `complement(join(...))` means.
            exons.sort()
            dna = "".join(seq[a:b] for a, b in exons)
            if len(dna) != (fe - fs):
                warnings.append(
                    f"spliced CDS — {len(exons)} exons joined "
                    f"({len(dna)} bp); the {(fe - fs) - len(dna)} bp of "
                    f"intron between them are not counted")
        else:
            dna = (seq[fs:fe] if fe >= fs else seq[fs:] + seq[:fe])
        if (fstrand or 1) < 0:
            dna = _rc(dna)
        # ── /codon_start ─────────────────────────────────────────────
        # A CDS whose coding frame begins 1 or 2 bases into the feature (a
        # partial 5' end) reads in the WRONG FRAME otherwise — every codon,
        # and so every number in the answer, shifted by one base.
        cs_quals = getattr(cds_feat, "qualifiers", None) or {}
        cs_raw = (cs_quals.get("codon_start") or ["1"])[0]
        try:
            codon_start = int(str(cs_raw).strip())
        except (TypeError, ValueError):
            codon_start = 1
        if codon_start in (2, 3):
            dna = dna[codon_start - 1:]
            warnings.append(
                f"/codon_start={codon_start} — the first {codon_start - 1} "
                f"base(s) are outside the reading frame and were trimmed")
        source = f"{key}:{cds_label or 'CDS'}"

    if len(dna) < 3:
        return ({"error": "sequence is shorter than one codon"}, 422)
    if len(dna) % 3:
        warnings.append(f"length {len(dna)} is not a multiple of 3 — the "
                        f"trailing {len(dna) % 3} base(s) are ignored")

    # ── codon table ──────────────────────────────────────────────────
    taxid = _sanitize_label(payload.get("taxid"), max_len=40)
    if not taxid:
        taxid = str(_get_setting("active_codon_table", "") or "").strip()
    table = _codon_tables_get(taxid) if taxid else None
    if table is None:
        if taxid:
            return ({"error": f"no codon table with taxid {taxid!r}; call "
                              f"list-codon-tables"}, 404)
        tables = _codon_tables_load()
        table = tables[0] if tables else None
        if table is None:
            return ({"error": "no codon tables installed; call "
                              "add-codon-table"}, 422)
        warnings.append(f"no active codon table set — used "
                        f"{table.get('name') or table.get('taxid')}")
    raw = table.get("codons") or table.get("raw") or {}
    if not raw:
        return ({"error": f"codon table {table.get('name')!r} has no usage "
                          f"data"}, 422)
    tt = payload.get("transl_table")
    if tt is not None:
        tt = _coerce_int(tt, name="transl_table")
        if isinstance(tt, str):
            return ({"error": tt}, 400)

    try:
        rare_w = float(payload.get("rare_w", _CODON_RARE_W))
    except (TypeError, ValueError):
        return ({"error": "'rare_w' must be a number"}, 400)
    if not (0.0 < rare_w <= 1.0):
        return ({"error": "'rare_w' must be in (0, 1]"}, 400)
    ramp_n = _coerce_int(payload.get("ramp_codons", _CODON_RAMP_CODONS),
                         name="ramp_codons")
    if isinstance(ramp_n, str):
        return ({"error": ramp_n}, 400)
    if ramp_n < 0:
        return ({"error": "'ramp_codons' must be >= 0"}, 400)

    try:
        weights = _codon_relative_adaptiveness(dna, raw, transl_table=tt)
        cai = _codon_cai(dna, raw, transl_table=tt)
    except Exception as exc:
        _log.exception("analyse-cds: codon scoring failed")
        return ({"error": f"codon analysis failed: "
                          f"{_scrub_path(str(exc))}"}, 500)

    rare = [w for w in weights if w["w"] < rare_w]
    # Tandem runs — consecutive CONTRIBUTING codons that are all rare. Indices
    # come from `_codon_relative_adaptiveness`, which skips stops and
    # single-codon families, so adjacency is judged on that list rather than on
    # raw codon index: a Met between two rare codons does not break the run in
    # any way a ribosome experiences.
    runs: "list[dict]" = []
    run_start = None
    for pos, w in enumerate(weights):
        is_rare = w["w"] < rare_w
        if is_rare and run_start is None:
            run_start = pos
        elif not is_rare and run_start is not None:
            if pos - run_start >= 2:
                runs.append({"start_codon": weights[run_start]["index"],
                             "length": pos - run_start,
                             "codons": [x["codon"]
                                        for x in weights[run_start:pos]]})
            run_start = None
    if run_start is not None and len(weights) - run_start >= 2:
        runs.append({"start_codon": weights[run_start]["index"],
                     "length": len(weights) - run_start,
                     "codons": [x["codon"] for x in weights[run_start:]]})
    runs.sort(key=lambda r: -r["length"])

    # Worst CAI window — a local collapse hidden inside a healthy average.
    import math as _math
    win = min(30, len(weights))
    worst = None
    if win >= 5:
        logs = [_math.log(max(x["w"], 1e-10)) for x in weights]
        acc = 0.0
        for i in range(len(logs)):
            acc += logs[i]
            if i >= win:
                acc -= logs[i - win]
            if i >= win - 1:
                val = _math.exp(acc / win)
                if worst is None or val < worst[0]:
                    worst = (val, weights[i - win + 1]["index"],
                             weights[i]["index"])

    ramp = None
    if ramp_n and len(weights) > ramp_n:
        head = [_math.log(max(x["w"], 1e-10)) for x in weights[:ramp_n]]
        tail = [_math.log(max(x["w"], 1e-10)) for x in weights[ramp_n:]]
        ramp = {
            "codons": ramp_n,
            "cai_first": round(_math.exp(sum(head) / len(head)), 4),
            "cai_rest": round(_math.exp(sum(tail) / len(tail)), 4),
        }
        ramp["slower_than_body"] = ramp["cai_first"] < ramp["cai_rest"]

    n_codons = len(dna) // 3
    _log_event("analyse.cds", n_codons=n_codons, cai=round(cai, 3),
               taxid=table.get("taxid"), via="agent")
    return {
        "ok": True,
        "source": source,
        "cds_label": cds_label,
        "length_bp": len(dna),
        "n_codons": n_codons,
        "n_scored_codons": len(weights),
        "gc": round(_codon_gc(dna), 2),
        "gc3": round(_codon_gc3(dna), 2),
        "cai": round(cai, 4),
        "rare_w": rare_w,
        "n_rare": len(rare),
        "rare_pct": (round(len(rare) * 100.0 / len(weights), 2)
                     if weights else 0.0),
        "rare_codons": [{"index": r["index"], "codon": r["codon"],
                         "aa": r["aa"], "w": round(r["w"], 4)}
                        for r in rare[:200]],
        "tandem_rare_runs": runs[:50],
        "worst_window": ({"cai": round(worst[0], 4),
                          "start_codon": worst[1],
                          "end_codon": worst[2],
                          "window_codons": win} if worst else None),
        "ramp": ramp,
        "table": {"name": table.get("name"), "taxid": table.get("taxid")},
        "warnings": warnings,
        "ignored": _agent_ignored_keys(payload, {
            "sequence", "id", "name", "feature", "taxid", "transl_table",
            "rare_w", "ramp_codons"}),
    }


@_agent_endpoint("list-codon-tables")
def _h_list_codon_tables(app, payload):
    """Available codon usage tables. Returns name, taxid, source
    (builtin/kazusa/user) and date added for each entry, plus the
    `active_taxid` field showing which is currently the persisted
    default (read by SynthesisScreen's `_init_codon_table`). Use
    the taxid as the `table` arg to ``optimize-protein``."""
    return {
        "tables": [
            {"name":   e.get("name", "?"),
             "taxid":  e.get("taxid", ""),
             "source": e.get("source", "user"),
             "added":  e.get("added", "")}
            for e in _codon_tables_load()
        ],
        "active_taxid": str(
            _get_setting("active_codon_table", "") or "",
        ),
    }


@_agent_endpoint("set-active-codon-table", write=True)
def _h_set_active_codon_table(app, payload):
    """Persist the active codon-table preference. Body:
    ``{taxid: str}`` — pass empty string to clear (Synthesis then
    falls back to the first registry entry, K12 by convention).

    Sweep #20: the Synthesis protein tab honors this setting on
    open via `_init_codon_table`, so an agent pre-setting the
    user's preferred organism (e.g. "83333" for E. coli K12,
    "559292" for S. cerevisiae) means the next time the user
    opens Synthesis they're already in the right table. Unknown
    taxids are rejected 404; empty string is accepted (clears)."""
    taxid_raw = payload.get("taxid")
    if taxid_raw is None:
        return ({"error": "missing 'taxid'"}, 400)
    if not isinstance(taxid_raw, str):
        return ({"error": "'taxid' must be a string"}, 400)
    taxid = taxid_raw.strip()
    if taxid:
        # Verify the taxid resolves to a real registry entry —
        # otherwise the setting is dead weight.
        if _codon_tables_get(taxid) is None:
            return ({"error":
                      f"no codon table with taxid {taxid!r}; "
                      f"see list-codon-tables"}, 404)
    _set_setting("active_codon_table", taxid)
    _log_event(
        "settings.changed",
        key="active_codon_table",
        value=taxid, via="agent",
    )
    return {"ok": True, "active_taxid": taxid}


@_agent_endpoint("optimize-protein")
def _h_optimize_protein(app, payload):
    """Codon-optimize a 1-letter AA sequence to DNA using a codon
    table from the registry. Body: ``{protein, table?, stops?,
    transl_table?, mode?, hazard_hosts?, forbidden_motifs?, min_gc?,
    max_gc?, gc_window?, avoid_repeats_with?, max_repeat?}`` where

    `hazard_hosts` pulls built-in host-expression hazards to scrub out of the
    CDS — ``"plant"`` (polyA signal AATAAA, cryptic splice donor GTRAG),
    ``"mammalian"``, ``"bacterial"`` (internal Shine-Dalgarno) — and
    `forbidden_motifs` adds your own, either ``["AATAAA", ...]`` or
    ``{"name": "MOTIF"}``. Any IUPAC pattern is accepted; the reverse
    complement is scrubbed too. `min_gc` / `max_gc` (PERCENT) then pull every
    `gc_window`-base window (default 50) into that band by further synonymous
    swaps, guarded so they cannot re-create a motif the scrub just removed.

    `avoid_repeats_with` takes DNA sequences you have ALREADY built and breaks
    any perfect-identity run of `max_repeat`+ bases (default 25) that the new
    design shares with them. Same-protein variants optimized one at a time all
    converge on the same table and end up sharing hundreds of bases, which is
    a recombination substrate when they meet in one cell; pass the panel built
    so far and each new member comes out diversified.

    The GC-window and diversify passes are heavier (roughly quadratic), so
    they are capped at 1500 codons — a protein longer than that with `min_gc`/
    `max_gc` or `avoid_repeats_with` set returns 400 (optimize it without those
    passes, or scrub in segments). The plain optimize and the motif scrub have
    no such limit.

    All of these are opt-in — omit them and the output is exactly what it was.
    Because a constraint can be physically unreachable (a Lys/Ile/Phe/Asn
    stretch has no synonym richer than one G/C per codon, so it tops out near
    33% and a 35% floor is impossible there), the response REPORTS what was
    achieved rather than implying success: `fixes` lists every swap made,
    `remaining_motifs` names any motif still present, `remaining_repeats`
    counts shared windows left, and `gc_window_min` / `gc_window_max` give the
    achieved window range. Check those, don't assume.

    The response also reports `gc3` (third-position "wobble" GC) alongside
    overall `gc`, and a `warnings` list. GC3 is where synonymous choice lives,
    so it swings far more than overall GC — and `max_cai` on an AT- or GC-biased
    host stacks that host's single most-frequent codon everywhere, collapsing
    GC3 (a silencing / mRNA-instability risk) while overall GC still looks
    healthy. When GC3 is extreme a warning says so; for such a host `frequency`
    mode (which matches the native GC3) is usually the safer choice than
    `max_cai`, or set `min_gc` to pull the wobble position back up.

    `mode` picks the strategy — ``"frequency"`` (default) matches the host's
    codon-usage DISTRIBUTION, deliberately using rare codons at their natural
    rate, so the result can score a lower CAI than the wild-type gene and an
    individual residue can get a worse synonym than a native gene used;
    ``"max_cai"`` gives every position its most-frequent synonym, which
    maximises CAI and can never downgrade a residue, at the cost of a much
    more repetitive sequence. Only one table (E. coli K12) ships — use
    ``add-codon-table`` to fetch your host from Kazusa by taxid, build one
    from an NCBI genome, or import a local CDS FASTA. `table`
    is a taxid (defaults to E. coli K12 = 83333) and `stops`
    (0–3, default 1) is how many stop codons to append when `protein`
    has no trailing '*'. A trailing '*' run in `protein` (1–3) is honored
    verbatim and overrides `stops`; a single stop is TAA, 2–3 are
    frequency-matched to the table's stop usage. `transl_table` is an
    optional NCBI genetic-code id (default 1 = standard); set it (e.g. 4
    Mycoplasma, 6 ciliate, a mito code) so a host with a reassigned codon
    optimises and terminates against the code it actually reads. Read-only —
    doesn't touch the loaded record."""
    protein = _sanitize_label(payload.get("protein"),
                                max_len=10_000).upper()
    if not protein:
        return ({"error": "missing 'protein'"}, 400)
    invalid_aa = [c for c in protein if c not in "ACDEFGHIKLMNPQRSTVWY*"]
    if invalid_aa:
        return ({"error":
                  f"non-canonical amino acids in 'protein': "
                  f"{''.join(sorted(set(invalid_aa)))!r}"}, 400)
    stops_raw = payload.get("stops", 1)
    try:
        stops = int(stops_raw)
    except (TypeError, ValueError):
        return ({"error": "'stops' must be an integer between 0 and 3"}, 400)
    if not 0 <= stops <= 3:
        return ({"error": "'stops' must be between 0 and 3"}, 400)
    # Optional NCBI genetic-code id (transl_table): default None ⇒ standard
    # table 1. Lets a host with a reassigned codon (Mycoplasma 4, ciliate 6,
    # a mito code) optimise + terminate correctly; an unknown id degrades to
    # the standard code (with a log warning) rather than erroring.
    transl_raw = payload.get("transl_table")
    transl_table = None
    if transl_raw is not None:
        try:
            transl_table = int(transl_raw)
        except (TypeError, ValueError):
            return ({"error": "'transl_table' must be an NCBI genetic-code "
                      "integer (e.g. 1, 4, 6, 11)"}, 400)
    mode = payload.get("mode", "frequency")
    if mode not in _CODON_MODES:
        return ({"error": f"'mode' must be one of "
                          f"{', '.join(repr(m) for m in _CODON_MODES)}"}, 400)
    taxid = _sanitize_accession(payload.get("table")) or "83333"
    entry = _codon_tables_get(taxid)
    if entry is None:
        # Point at the ADD path too: only K12 ships, so "no table for my host"
        # is the expected first result, and a caller who only hears
        # "see list-codon-tables" concludes the host is unsupported and
        # hand-builds a table instead of fetching one.
        return ({"error": f"no codon table with taxid {taxid!r}; "
                          f"see list-codon-tables, or add one with "
                          f"add-codon-table (source 'kazusa' fetches by "
                          f"taxid, 'genome' builds from an NCBI assembly, "
                          f"'file' from a local CDS FASTA)"}, 404)
    # ── optional scrub passes (all no-ops unless asked for) ──────────────
    scrub: "dict[str, str]" = {}
    hosts = payload.get("hazard_hosts")
    if hosts is not None:
        if not isinstance(hosts, (list, str)):
            return ({"error": "'hazard_hosts' must be a list of host names"},
                    400)
        try:
            scrub.update(_codon_hazard_motifs(hosts))
        except ValueError as exc:
            return ({"error": str(exc)}, 400)
    motifs = payload.get("forbidden_motifs")
    if motifs is not None:
        if isinstance(motifs, dict):
            pairs = [(str(k), v) for k, v in motifs.items()]
        elif isinstance(motifs, list):
            pairs = [(str(m), m) for m in motifs]
        else:
            return ({"error": "'forbidden_motifs' must be a list of motifs "
                              "or a {name: motif} object"}, 400)
        for name, motif in pairs:
            text = str(motif or "").upper()
            if not text:
                return ({"error": "'forbidden_motifs' contains an empty "
                                  "motif"}, 400)
            try:
                _iupac_pattern(text)
            except ValueError:
                return ({"error": f"motif {text!r} is not a valid IUPAC "
                                  f"pattern"}, 400)
            scrub[name] = text
    gc_window = _CODON_GC_WINDOW_DEFAULT
    if payload.get("gc_window") is not None:
        gc_window = _coerce_int(payload["gc_window"], name="gc_window")
        if isinstance(gc_window, str):
            return ({"error": gc_window}, 400)
    if gc_window <= 0:
        return ({"error": "'gc_window' must be a positive integer"}, 400)
    bounds: "dict[str, float]" = {}
    for key in ("min_gc", "max_gc"):
        if payload.get(key) is None:
            continue
        try:
            bounds[key] = float(payload[key])
        except (TypeError, ValueError):
            return ({"error": f"'{key}' must be a number (percent, 0-100)"},
                    400)
        if not 0.0 <= bounds[key] <= 100.0:
            return ({"error": f"'{key}' must be between 0 and 100 (percent)"},
                    400)
    refs_in = payload.get("avoid_repeats_with")
    refs: list = []
    if refs_in is not None:
        if isinstance(refs_in, str):
            refs_in = [refs_in]
        if not isinstance(refs_in, list):
            return ({"error": "'avoid_repeats_with' must be a DNA string or "
                              "a list of DNA strings"}, 400)
        for ref in refs_in:
            if not isinstance(ref, str):
                return ({"error": "'avoid_repeats_with' entries must be DNA "
                                  "strings"}, 400)
            text, err = _sanitize_bases(ref)
            if err or not text:
                return ({"error": "'avoid_repeats_with' contains a value that "
                                  "is not a DNA sequence"
                                  + (f": {err}" if err else "")}, 400)
            refs.append(text)
    max_repeat = _CODON_REPEAT_RUN_DEFAULT
    if payload.get("max_repeat") is not None:
        max_repeat = _coerce_int(payload["max_repeat"], name="max_repeat")
        if isinstance(max_repeat, str):
            return ({"error": max_repeat}, 400)
    if max_repeat <= 0:
        return ({"error": "'max_repeat' must be a positive integer"}, 400)

    # The GC-window and diversify passes are O(n²); refuse them past a
    # generous CDS length rather than pegging a worker thread (the plain
    # optimize + the linear motif scrub stay uncapped).
    n_codons = len(protein.rstrip("*"))
    if (bounds or refs) and n_codons > _CODON_SCRUB_MAX_CODONS:
        which = " and ".join(
            w for w, on in (("GC-window (min_gc/max_gc)", bool(bounds)),
                            ("repeat diversification (avoid_repeats_with)",
                             bool(refs))) if on)
        return ({"error":
                  f"{which} is limited to {_CODON_SCRUB_MAX_CODONS} codons; "
                  f"this protein is {n_codons}. Optimize without that pass, or "
                  f"scrub the CDS in segments."}, 400)

    try:
        dna = _codon_optimize(protein, entry["raw"], stops=stops,
                              transl_table=transl_table, mode=mode)
        body = protein.rstrip("*")
        # Does `dna` actually END in a stop codon we appended? The scrub
        # passes skip their last codon when it does, so that they never
        # substitute a synthetic stop. Leaving that skip on for a STOP-FREE
        # request (`stops=0`, e.g. a domain destined for a C-terminal fusion)
        # made the last real codon unreachable: a forbidden site overlapping
        # the C-terminus came back as "cannot be removed without changing the
        # protein" when a plain synonymous swap was available. A residual Type
        # IIS site at the C-terminus is exactly what kills a Golden Gate
        # reaction, so the claim mattered.
        n_stops = (len(protein) - len(body)) or max(0, int(stops))
        has_stop = n_stops > 0
        fixes: list = []
        if scrub:
            dna, site_fixes = _codon_fix_sites(
                dna, body, entry["raw"], sites=scrub,
                has_appended_stop=has_stop,
                transl_table=transl_table)
            fixes.extend(site_fixes)
        if refs:
            dna, rep_fixes = _codon_diversify(
                dna, body, entry["raw"], refs, min_run=max_repeat,
                sites=scrub or None, has_appended_stop=has_stop,
                transl_table=transl_table)
            fixes.extend(rep_fixes)
        if bounds:
            # LAST, and guarded by BOTH the motif set and the repeat k-mers,
            # so a GC swap can't re-create a motif the scrub removed nor a
            # shared run the diversification broke. `remaining_*` below are
            # still measured on the FINAL sequence, not trusted from here.
            dna, gc_fixes = _codon_fix_gc_window(
                dna, body, entry["raw"], window=gc_window,
                min_gc=bounds.get("min_gc"), max_gc=bounds.get("max_gc"),
                sites=scrub or None, has_appended_stop=has_stop,
                transl_table=transl_table,
                avoid_kmers=_codon_kmer_set(refs, max_repeat) if refs else None,
                kmer_len=max_repeat if refs else 0)
            fixes.extend(gc_fixes)
    except ValueError as exc:
        return ({"error": str(exc)}, 400)

    # Re-MEASURE rather than assume: a constraint can be unreachable (a
    # Lys/Ile/Phe/Asn stretch tops out near 33% GC, so a 35% floor is
    # impossible there), and reporting the achieved numbers beats implying
    # success. `remaining_motifs` does the same for the scrub.
    win_lo, win_hi = _codon_gc_window_range(dna, gc_window)
    remaining_repeats = len(_codon_shared_runs(dna, refs, max_repeat)) if refs \
        else 0
    remaining = sorted({
        name for name, motif in scrub.items()
        if _iupac_pattern(motif).search(dna)
        or _iupac_pattern(_rc(motif)).search(dna)
    })
    # GC3 (third-position GC) is where synonymous choice lives, so it swings
    # far more than overall GC — a healthy overall number can hide a collapsed
    # wobble position, which is a real silencing / mRNA-instability risk. Warn
    # when it is extreme so the caller isn't lulled by the top-line GC.
    gc3 = round(_codon_gc3(dna), 2)
    warnings: list = []
    if ((len(dna) // 3) >= _CODON_GC3_MIN_CODONS
            and (gc3 < _CODON_GC3_LOW or gc3 > _CODON_GC3_HIGH)):
        side = "well below" if gc3 < _CODON_GC3_LOW else "well above"
        msg = (f"third-position GC (GC3) is {gc3}% — {side} the typical range, "
               f"so the wobble position is nearly monotonous. That is a "
               f"silencing / mRNA-instability / synthesis-complexity risk the "
               f"overall GC ({round(_codon_gc(dna), 2)}%) hides.")
        if mode == "max_cai":
            msg += (" max_cai stacks each residue's single most-frequent codon, "
                    "which skews GC3 hardest on AT- or GC-biased hosts; "
                    "'frequency' mode matches the host's native GC3, or set "
                    "min_gc to pull it back.")
        warnings.append(msg)
    return {
        "ok":     True,
        "protein":      protein,
        "table":        taxid,
        "table_name":   entry.get("name", "?"),
        "transl_table": transl_table or 1,
        "mode":         mode,
        "dna":          dna,
        "length":       len(dna),
        "n_codons":     len(dna) // 3,
        "fixes":             fixes,
        "remaining_motifs":  remaining,
        "max_repeat":        max_repeat,
        "remaining_repeats": remaining_repeats,
        "gc_window":         gc_window,
        "gc_window_min":    round(win_lo, 2),
        "gc_window_max":    round(win_hi, 2),
        "warnings":          warnings,
        # Display/report-only metrics, matching the Synthesis + Mutato toasts.
        # Scored against the SAME genetic code the sequence was optimised under
        # so a reassigned-codon host isn't mis-scored. Additive (V1_GATE-safe).
        "cai":          round(_codon_cai(dna, entry["raw"],
                                         transl_table=transl_table), 4),
        "gc":           round(_codon_gc(dna), 2),
        "gc3":          gc3,
    }


@_agent_endpoint("list-plasmid-statuses")
def _h_list_plasmid_statuses(app, payload):
    """Return the canonical workflow-status vocabulary. Discoverability
    helper so agents don't have to hard-code `DESIGNING`/`CLONING`/
    `SEQUENCING`/`VERIFIED` strings."""
    return {
        "ok":       True,
        "statuses": list(_PLASMID_STATUS_VALUES),
        "colors":   dict(_PLASMID_STATUS_COLORS),
    }


@_agent_endpoint("list-entry-vectors")
def _h_list_entry_vectors(app, payload):
    """Return every grammar/role that currently has an assigned entry
    vector. Each item carries ``{grammar_id, role, name, size, source}``
    — `gb_text` is omitted to keep the response small; fetch it via
    `get-entry-vector` for a specific grammar.

    A single grammar can legitimately hold SEVERAL entry vectors under
    different ``role`` slots (Golden Braid's Alpha1 / Alpha2 / Omega1 /
    Omega2 acceptors, plus the empty-role L0 vector), so the listing now
    surfaces ``role`` — what looked like "the same grammar listed twice"
    two distinct roles. Genuine ``(grammar_id, role)`` duplicates are
    collapsed."""
    out = []
    seen = set()
    for e in _load_entry_vectors():
        gid = e.get("grammar_id", "")
        role = e.get("role", "") or ""
        if (gid, role) in seen:
            continue
        seen.add((gid, role))
        out.append({
            "grammar_id": gid,
            "role":       role,
            "name":       e.get("name", ""),
            "size":       e.get("size", 0),
            "source":     e.get("source", ""),
        })
    return {"ok": True, "entry_vectors": out}


@_agent_endpoint("get-entry-vector")
def _h_get_entry_vector(app, payload):
    """Return the entry vector for a grammar. Body: ``{grammar_id}``.
    Includes the full `gb_text` in the response so the agent can parse
    + render it without a follow-up call. Returns `null` for `vector`
    if the grammar has no assigned vector."""
    gid = payload.get("grammar_id")
    if not isinstance(gid, str) or not gid:
        return ({"error": "missing or non-string 'grammar_id'"}, 400)
    vec = _get_entry_vector(gid)
    return {"ok": True, "grammar_id": gid, "vector": vec}


def _agent_grammar_dict(payload: dict) -> "dict | str":
    """Coerce + validate a payload describing a cloning grammar.
    Returns the cleaned grammar dict or a string error message.

    Shape-checks each field type without auto-defaulting missing
    ones; the agent must pass everything it wants persisted. Tighter
    than the UI because we don't have an editing modal to fix up
    on save."""
    gid = payload.get("id")
    if not isinstance(gid, str) or not gid.strip():
        return "missing or non-string 'id'"
    if gid in _BUILTIN_GRAMMARS:
        return f"refusing to overwrite built-in grammar {gid!r}"
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return "missing or non-string 'name'"
    enzyme = payload.get("enzyme")
    if not isinstance(enzyme, str) or not enzyme.strip():
        return "missing or non-string 'enzyme'"
    site = payload.get("site")
    if not isinstance(site, str) or not site.strip():
        return "missing or non-string 'site'"
    spacer = payload.get("spacer", "")
    pad = payload.get("pad", "")
    if not isinstance(spacer, str) or not isinstance(pad, str):
        return "'spacer' and 'pad' must be strings"
    forb = payload.get("forbidden_sites", {})
    if not isinstance(forb, dict) or not all(
            isinstance(k, str) and isinstance(v, str)
            for k, v in forb.items()):
        return "'forbidden_sites' must be {enzyme: site} string→string"
    positions = payload.get("positions")
    if not isinstance(positions, list) or not positions:
        return "missing non-empty 'positions' list"
    cleaned_positions: list[dict] = []
    for i, pos in enumerate(positions):
        if not isinstance(pos, dict):
            return f"positions[{i}] is not a dict"
        for k in ("name", "type", "oh5", "oh3"):
            if not isinstance(pos.get(k), str):
                return f"positions[{i}] missing string field {k!r}"
        cleaned_positions.append({
            "name": _sanitize_label(pos["name"]),
            "type": _sanitize_label(pos["type"]),
            "oh5":  pos["oh5"].upper(), "oh3": pos["oh3"].upper(),
            "color": _sanitize_label(pos.get("color")) or "white",
        })
    coding_types = payload.get("coding_types") or []
    if not isinstance(coding_types, list) or not all(
            isinstance(x, str) for x in coding_types):
        return "'coding_types' must be a list of strings"
    type_to_insdc = payload.get("type_to_insdc") or {}
    if not isinstance(type_to_insdc, dict) or not all(
            isinstance(k, str) and isinstance(v, str)
            for k, v in type_to_insdc.items()):
        return "'type_to_insdc' must be {part_type: insdc_type} strings"
    return {
        "id":              gid.strip(),
        "name":            _sanitize_label(name),
        "enzyme":          _sanitize_label(enzyme),
        "level_up_enzyme": (payload.get("level_up_enzyme")
                              if isinstance(payload.get("level_up_enzyme"), str)
                              else ""),
        "site":            site.strip().upper(),
        "spacer":          spacer,
        "pad":             pad,
        "forbidden_sites": dict(forb),
        "positions":       cleaned_positions,
        "coding_types":    list(coding_types),
        "type_to_insdc":   dict(type_to_insdc),
        "catalog":         [],
        "editable":        True,
    }


@_agent_endpoint("list-grammars")
def _h_list_grammars(app, payload):
    """Return every grammar (built-in + custom). Each item carries
    ``{id, name, enzyme, level_up_enzyme, editable, n_positions,
    catalog_size}``. Fetch the full grammar via `get-grammar`."""
    out = []
    for gid, g in _all_grammars().items():
        out.append({
            "id":              gid,
            "name":            g.get("name", gid),
            "enzyme":          g.get("enzyme", ""),
            "level_up_enzyme": g.get("level_up_enzyme", ""),
            "editable":        bool(g.get("editable")),
            "n_positions":     len(g.get("positions") or []),
            "catalog_size":    len(g.get("catalog") or []),
        })
    return {"ok": True, "grammars": out}


@_agent_endpoint("get-grammar")
def _h_get_grammar(app, payload):
    """Return the full grammar dict for `grammar_id`. Body:
    ``{grammar_id}``."""
    gid = payload.get("grammar_id")
    if not isinstance(gid, str) or not gid:
        return ({"error": "missing or non-string 'grammar_id'"}, 400)
    grammars = _all_grammars()
    g = grammars.get(gid)
    if g is None:
        return ({"error": f"unknown grammar id {gid!r}"}, 404)
    return {"ok": True, "grammar": g}


@_agent_endpoint("create-grammar", write=True)
def _h_create_grammar(app, payload):
    """Create a new custom grammar. Body: every field in the schema
    listed above. Returns 400 on validation failure, 409 if a grammar
    with the same id already exists. Built-in ids are reserved."""
    g = _agent_grammar_dict(payload)
    if isinstance(g, str):
        return ({"error": g}, 400)
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        custom = _load_custom_grammars()
        if any(e.get("id") == g["id"] for e in custom):
            return ({"error": (
                f"grammar id {g['id']!r} already exists; use update-grammar"
            )}, 409)
        custom.append(g)
        if (err := _agent_save_or_500(
                lambda: _save_custom_grammars(custom),
                "cloning_grammars")) is not None:
            return err
    return {"ok": True, "grammar_id": g["id"]}


@_agent_endpoint("update-grammar", write=True)
def _h_update_grammar(app, payload):
    """Replace an existing custom grammar by id. Built-in ids are
    refused (mirrors the UI's `editable=False` gate)."""
    g = _agent_grammar_dict(payload)
    if isinstance(g, str):
        return ({"error": g}, 400)
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        custom = _load_custom_grammars()
        for i, e in enumerate(custom):
            if e.get("id") == g["id"]:
                # Preserve the existing catalog list — `_agent_grammar_dict`
                # always sets `catalog: []`, but the user might have built
                # one up via Parts Bin. Keep it unless the payload
                # explicitly clears it.
                if "catalog" not in payload:
                    g["catalog"] = list(e.get("catalog") or [])
                custom[i] = g
                if (err := _agent_save_or_500(
                        lambda: _save_custom_grammars(custom),
                        "cloning_grammars")) is not None:
                    return err
                return {"ok": True, "grammar_id": g["id"]}
    return ({"error": f"unknown grammar id {g['id']!r}"}, 404)


@_agent_endpoint("delete-grammar", write=True)
def _h_delete_grammar(app, payload):
    """Delete a custom grammar by id. Built-ins are refused. Also
    clears any entry vector bound to the grammar so a future
    create-grammar with the same id starts fresh."""
    gid = payload.get("grammar_id")
    if not isinstance(gid, str) or not gid:
        return ({"error": "missing or non-string 'grammar_id'"}, 400)
    if gid in _BUILTIN_GRAMMARS:
        return ({"error": (
            f"refusing to delete built-in grammar {gid!r}"
        )}, 400)
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        custom = _load_custom_grammars()
        new_list = [e for e in custom if e.get("id") != gid]
        if len(new_list) == len(custom):
            return ({"error": f"unknown grammar id {gid!r}"}, 404)
        if (err := _agent_save_or_500(
                lambda: _save_custom_grammars(new_list),
                "cloning_grammars")) is not None:
            return err
    # Also clear the bound entry vector if any. Failure here is
    # non-fatal — the orphan record gets pruned next save.
    try:
        _set_entry_vector(gid, None)
    except (OSError, RuntimeError):
        pass
    return {"ok": True, "grammar_id": gid}


def _agent_primer_collection_list(coll):
    """SC-J: the primers belonging to collection `coll`. Returns
    ``(list, None)`` or ``(None, (error_dict, status))``. ``coll`` empty / None
    → the default (live) library via `_load_primers()`; a named collection →
    ITS OWN stored ``primers`` (so a collection is a real partition for reads,
    not just writes), 404 if it doesn't exist."""
    if coll in (None, ""):
        return _load_primers(), None
    if not isinstance(coll, str):
        return None, ({"error": "'collection' must be a string"}, 400)
    target = next((c for c in _load_primer_collections()
                   if (c.get("name") or "") == coll.strip()), None)
    if target is None:
        return None, ({"error":
                        f"no primer collection named {coll.strip()!r}"}, 404)
    return list(target.get("primers") or []), None


@_agent_endpoint("list-primers")
def _h_list_primers(app, payload):
    """Return the primer library. Each item: ``{name, sequence,
    source, tm, status, type, date, notes}``. Sequence + notes can
    be long; clients on small terminals should paginate via the
    `limit` + `offset` body params.

    Pass ``{collection: "X"}`` to list only collection X's OWN primers (SC-J);
    omit it for the default/live library. The response echoes ``collection``."""
    coll_in = payload.get("collection")
    entries, cerr = _agent_primer_collection_list(coll_in)
    if cerr is not None:
        return cerr
    assert entries is not None          # cerr is None ⇒ entries populated
    raw_limit = payload.get("limit")
    raw_offset = payload.get("offset")
    # Sweep #35 (2026-05-26): clamp pagination to sane bounds and
    # reject JSON floats / booleans. Python's list slicing already
    # tolerates out-of-range slice indices, so pre-fix `int(1e308)`
    # didn't actually DoS, but it's tidier to validate at the
    # boundary than rely on slice clamping. Booleans pass
    # `isinstance(x, int)` in Python, so they're filtered out
    # explicitly.
    _LIST_PRIMERS_LIMIT_MAX = 10_000
    def _coerce_pagination(v, default: int, cap: int) -> int:
        if isinstance(v, bool) or not isinstance(v, int):
            return default
        return max(0, min(v, cap))
    limit = _coerce_pagination(
        raw_limit, default=1000, cap=_LIST_PRIMERS_LIMIT_MAX,
    )
    offset = _coerce_pagination(
        raw_offset, default=0, cap=_LIST_PRIMERS_LIMIT_MAX,
    )
    sliced = entries[offset:offset + limit]
    out = []
    for e in sliced:
        out.append({
            "name":     e.get("name", ""),
            "sequence": e.get("sequence", ""),
            "source":   e.get("source", ""),
            "tm":       e.get("tm"),
            "status":   e.get("status", ""),
            "type":     e.get("type", ""),
            "date":     e.get("date", ""),
            "notes":    e.get("notes", ""),
        })
    return {"ok": True, "primers": out, "count": len(out),
            "total": len(entries),
            "collection": (coll_in.strip()
                           if isinstance(coll_in, str) and coll_in.strip()
                           else "")}


@_agent_endpoint("get-primer")
def _h_get_primer(app, payload):
    """Return a primer-library entry by ``{name|id}`` or exact
    (case-insensitive) ``{sequence}``. Optional ``{collection}`` scopes the
    lookup to that collection's own primers (SC-J); omit it for the default
    library. An ambiguous name returns 409 with the candidate sequences; no
    match returns 404."""
    pool, cerr = _agent_primer_collection_list(payload.get("collection"))
    if cerr is not None:
        return cerr
    assert pool is not None             # cerr is None ⇒ pool populated
    seq = payload.get("sequence")
    name = payload.get("name") or payload.get("id")
    if isinstance(seq, str) and seq.strip():
        target = seq.strip().upper()
        for e in pool:
            if (e.get("sequence") or "").upper() == target:
                return {"ok": True, "primer": e}
        return ({"error": "no primer with that sequence"}, 404)
    if isinstance(name, str) and name.strip():
        key = name.strip()
        matches = [e for e in pool if (e.get("name") or "") == key]
        if not matches:
            return ({"error": f"no primer named {key!r}"}, 404)
        if len(matches) > 1:
            return ({"error":
                      f"{key!r} matches {len(matches)} primers — pass "
                      f"{{\"sequence\": …}} to pick one",
                      "sequences": [e.get("sequence", "")
                                     for e in matches]}, 409)
        return {"ok": True, "primer": matches[0]}
    return ({"error": "missing 'name', 'id', or 'sequence'"}, 400)


_PRIMER_STATUS_VALUES: tuple = ("Designed", "Ordered", "Validated")


def _agent_primer_dict(payload: dict) -> "dict | str":
    """Coerce + validate a primer payload. Returns the cleaned dict
    or a string error. Used by both create and update paths.

    Sequence is the canonical identity field — `_dedupe_primers_by_sequence`
    treats two primers with the same uppercase sequence as duplicates
    and keeps the first.
    """
    name = _sanitize_label(payload.get("name"))
    if not name:
        return "missing or non-string 'name'"
    raw_seq = payload.get("sequence")
    if not isinstance(raw_seq, str) or not raw_seq.strip():
        return "missing or non-string 'sequence'"
    seq, seq_err = _sanitize_bases(raw_seq.upper())
    if seq_err:
        return f"'sequence' rejected: {seq_err}"
    if not seq:
        return "'sequence' empty after sanitisation"
    if len(seq) > 500:
        return "'sequence' too long (max 500 bp)"
    source = payload.get("source", "")
    if source is not None and not isinstance(source, str):
        return "'source' must be a string"
    status = payload.get("status", "Designed")
    if status and status not in _PRIMER_STATUS_VALUES:
        return (f"'status' must be one of {_PRIMER_STATUS_VALUES}")
    ptype = payload.get("type", "")
    if ptype is not None and not isinstance(ptype, str):
        return "'type' must be a string"
    notes = _sanitize_note(payload.get("notes", ""))
    tm = payload.get("tm")
    if tm is not None and not isinstance(tm, (int, float)):
        return "'tm' must be a number or null"
    return {
        "name":     name,
        "sequence": seq,
        "source":   _sanitize_label(source, max_len=300),
        "status":   status or "Designed",
        "type":     _sanitize_label(ptype),
        "tm":       tm,
        "notes":    notes or "",
        "date":     _datetime.now().strftime("%Y-%m-%d"),
    }


@_agent_endpoint("create-primer", write=True)
def _h_create_primer(app, payload):
    """Add a primer to the library. Body: ``{name, sequence, source?,
    status?, type?, tm?, notes?, collection?}``.

    Without ``collection`` the primer lands in the ACTIVE primer collection.
    With ``collection`` it's filed into that NAMED collection — and if that
    collection doesn't exist the call FAILS with 404 (create it first with
    `create-primer-collection`). SC-D: ``collection`` used to be silently
    ignored, so a primer the caller believed was filed into a project
    collection actually landed in the default — an outcome-changing param
    must never be dropped with ``ok:true``. Duplicate sequences
    (case-insensitive, within the target collection) return 409."""
    p = _agent_primer_dict(payload)
    if isinstance(p, str):
        return ({"error": p}, 400)
    target = p["sequence"]
    coll = payload.get("collection")
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        if coll not in (None, ""):
            if not isinstance(coll, str):
                return ({"error": "'collection' must be a string"}, 400)
            coll = coll.strip()
            active = _get_active_primer_collection_name() or ""
            if coll != active:
                # File into a SPECIFIC non-active named collection: write its
                # stored primers list directly (the live primers.json mirror
                # belongs to the ACTIVE collection and must stay untouched).
                colls = _load_primer_collections()
                tgt = next((c for c in colls
                            if (c.get("name") or "") == coll), None)
                if tgt is None:
                    return ({"error":
                              f"no primer collection named {coll!r}; create "
                              f"it with create-primer-collection"}, 404)
                bucket = tgt.setdefault("primers", [])
                for e in bucket:
                    if (e.get("sequence") or "").upper() == target:
                        return ({"error": (
                            f"primer with sequence {target!r} already exists "
                            f"in collection {coll!r} (name "
                            f"{e.get('name', '?')!r})"),
                            "existing_name": e.get("name", "")}, 409)
                bucket.insert(0, p)
                if (err := _agent_save_or_500(
                        lambda: _save_primer_collections(colls),
                        "primer_collections")) is not None:
                    return err
                return {"ok": True, "name": p["name"],
                        "sequence": p["sequence"], "collection": coll}
            # coll == active → fall through to the mirrored active path.
        entries = _load_primers()
        for e in entries:
            if (e.get("sequence") or "").upper() == target:
                return ({"error": (
                    f"primer with sequence {target!r} already exists "
                    f"(name {e.get('name', '?')!r}); use update-primer "
                    f"or delete-primer first."
                ), "existing_name": e.get("name", "")}, 409)
        entries.insert(0, p)  # MRU at index 0
        if (err := _agent_save_or_500(
                lambda: _save_primers(entries),
                "primers")) is not None:
            return err
    return {"ok": True, "name": p["name"], "sequence": p["sequence"],
            "collection": _get_active_primer_collection_name() or ""}


@_agent_endpoint("delete-primer", write=True)
def _h_delete_primer(app, payload):
    """Remove a primer by ``{name}`` (or ``{id}``) or exact ``{sequence}``,
    optionally scoped to a ``{collection}``. Body: one of
    ``{name}``/``{id}``/``{sequence}`` (+ optional ``collection``).

    Without ``collection`` it deletes from the default / live library; with
    ``collection`` it removes from THAT named collection's OWN primers — so
    you can prune a collection (e.g. strip stale mirror-pollution) without
    touching the default. Removing a primer that also lives in another
    collection is safe: the other copy stays. By name → the one named primer
    (409 if the name is ambiguous, listing sequences); by sequence → every
    exact match. 404 if no match. Refuses (409) to delete from a named
    collection while it's ACTIVE — switch to the default first so the live
    mirror stays consistent. Each write triggers the standard ``.bak``
    rotation (Settings → Restore from backup)."""
    seq = payload.get("sequence")
    name = payload.get("name") or payload.get("id")
    has_seq = isinstance(seq, str) and bool(seq.strip())
    has_name = isinstance(name, str) and bool(name.strip())
    if not (has_seq or has_name):
        return ({"error": "missing 'name', 'id', or 'sequence'"}, 400)
    coll = payload.get("collection")
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        named = None
        colls = None
        if coll not in (None, ""):
            if not isinstance(coll, str):
                return ({"error": "'collection' must be a string"}, 400)
            coll = coll.strip()
            if coll == (_get_active_primer_collection_name() or ""):
                return ({"error":
                          f"{coll!r} is the active primer collection — "
                          f"set-active-primer-collection to \"\" (default) "
                          f"before deleting from it, so the live mirror stays "
                          f"consistent"}, 409)
            colls = _load_primer_collections()
            named = next((c for c in colls
                          if (c.get("name") or "") == coll), None)
            if named is None:
                return ({"error":
                          f"no primer collection named {coll!r}"}, 404)
            entries = named.get("primers") or []
        else:
            entries = _load_primers()
        where = f" in collection {coll!r}" if named is not None else ""
        if has_seq:
            target = seq.strip().upper()
            new_list = [e for e in entries
                        if (e.get("sequence") or "").upper() != target]
            if len(new_list) == len(entries):
                return ({"error":
                          f"no primer with sequence {target!r}{where}"}, 404)
        else:
            key = name.strip()
            matches = [e for e in entries if (e.get("name") or "") == key]
            if not matches:
                return ({"error": f"no primer named {key!r}{where}"}, 404)
            if len(matches) > 1:
                return ({"error":
                          f"{key!r} matches {len(matches)} primers{where} — "
                          f"pass {{\"sequence\": …}} to pick one",
                          "sequences": [e.get("sequence", "")
                                         for e in matches]}, 409)
            # Remove exactly the matched primer (by identity), so a different
            # primer that happens to share its sequence is left intact.
            new_list = [e for e in entries if e is not matches[0]]
        if named is not None:
            assert colls is not None     # `named` set ⇒ colls was loaded
            named["primers"] = new_list
            if (err := _agent_save_or_500(
                    lambda: _save_primer_collections(colls),
                    "primer_collections")) is not None:
                return err
        else:
            if (err := _agent_save_or_500(
                    lambda: _save_primers(new_list), "primers")) is not None:
                return err
    return {"ok": True, "removed": len(entries) - len(new_list),
            "collection": coll if named is not None else ""}


@_agent_endpoint("delete-primers", write=True)
def _h_delete_primers(app, payload):
    """Batch-delete primers by name in ONE locked save — the machine-friendly
    path for a large cleanup (agent-API feedback SC-L), so a 600-item prune
    completes in a single call instead of N rate-limited `delete-primer`s.
    Body: ``{names: [str, …]}`` (``ids`` accepted as an alias), optional
    ``{collection}`` to scope to that named collection's OWN primers (refused
    409 if it's ACTIVE — same rule as `delete-primer`, so the live mirror
    stays consistent).

    Each name is matched against the target list: exactly one match → removed;
    zero → reported under ``not_found``; more than one → reported under
    ``ambiguous`` and NOT removed (pass the single `delete-primer {sequence}`
    to disambiguate). Every removal commits under a single ``_cache_lock`` +
    one ``.bak``-rotating save.

    Returns ``{ok, removed, removed_names, not_found, ambiguous, collection}``."""
    names = payload.get("names")
    if names is None:
        names = payload.get("ids")
    if not isinstance(names, list) or not names:
        return ({"error": "'names' (or 'ids') must be a non-empty list"}, 400)
    if len(names) > 5000:
        return ({"error": "too many names (max 5000 per batch)"}, 400)
    cleaned = []
    for nm in names:
        if not isinstance(nm, str) or not nm.strip():
            return ({"error":
                      "'names' must contain only non-empty strings"}, 400)
        cleaned.append(nm.strip())
    coll = payload.get("collection")
    with _state._cache_lock:
        named = None
        colls = None
        if coll not in (None, ""):
            if not isinstance(coll, str):
                return ({"error": "'collection' must be a string"}, 400)
            coll = coll.strip()
            if coll == (_get_active_primer_collection_name() or ""):
                return ({"error":
                          f"{coll!r} is the active primer collection — "
                          f"set-active-primer-collection to \"\" (default) "
                          f"before deleting from it, so the live mirror stays "
                          f"consistent"}, 409)
            colls = _load_primer_collections()
            named = next((c for c in colls
                          if (c.get("name") or "") == coll), None)
            if named is None:
                return ({"error":
                          f"no primer collection named {coll!r}"}, 404)
            entries = named.get("primers") or []
        else:
            entries = _load_primers()
        # Index by name so multiplicity (→ ambiguous) is detected in one pass.
        by_name: "dict[str, list]" = {}
        for e in entries:
            by_name.setdefault(e.get("name") or "", []).append(e)
        remove_ids: "set[int]" = set()
        removed_names, not_found, ambiguous = [], [], []
        for key in cleaned:
            matches = by_name.get(key, [])
            if not matches:
                not_found.append(key)
            elif len(matches) > 1:
                ambiguous.append(key)
            else:
                # Mark by identity so a different primer sharing this one's
                # sequence is left intact (mirrors `delete-primer`).
                remove_ids.add(id(matches[0]))
                removed_names.append(key)
        if remove_ids:
            new_list = [e for e in entries if id(e) not in remove_ids]
            if named is not None:
                assert colls is not None     # `named` set ⇒ colls was loaded
                named["primers"] = new_list
                if (err := _agent_save_or_500(
                        lambda: _save_primer_collections(colls),
                        "primer_collections")) is not None:
                    return err
            else:
                if (err := _agent_save_or_500(
                        lambda: _save_primers(new_list), "primers")) is not None:
                    return err
    return {"ok": True, "removed": len(removed_names),
            "removed_names": removed_names, "not_found": not_found,
            "ambiguous": ambiguous,
            "collection": coll if named is not None else ""}


@_agent_endpoint("list-primer-collections")
def _h_list_primer_collections(app, payload):
    """Every named primer collection. Each item: ``{name, n_primers}``.
    The unnamed "default" collection (= top-level primers.json) is
    surfaced as ``{name: "", n_primers: <count>}`` if it has content."""
    colls = _load_primer_collections()
    out = []
    # Top-level primers — surfaced as the "default" collection.
    top = _load_primers()
    if top:
        out.append({"name": "", "n_primers": len(top)})
    for c in colls:
        out.append({
            "name":      c.get("name", "?"),
            "n_primers": len(c.get("primers", []) or []),
        })
    return {"ok": True, "primer_collections": out,
            "active": _get_active_primer_collection_name() or ""}


@_agent_endpoint("set-active-primer-collection", write=True)
def _h_set_active_primer_collection(app, payload):
    """Switch the active primer collection AND re-point the live
    ``primers.json`` mirror at it. Body: ``{name}`` where `name` is one of
    the collections from `list-primer-collections` (empty string = the
    free-standing default library).

    Setting the active pointer alone is NOT enough: ``primers.json`` is the
    live mirror of the active collection (`_save_primers` keeps them in
    lock-step), so if the switch left it holding the PREVIOUS collection's
    primers, the very next `create-primer` / `_save_primers` would mirror that
    stale content back over the newly-active collection and destroy its
    primers. So — like the plasmid `set-active-collection` handler and the
    startup restore ([INV-83]) — a switch to a NAMED collection rewrites
    ``primers.json`` from that collection's stored primers (via
    `_safe_save_json_mirror`, which permits the legitimate shrink) and busts
    the cache. Switching to the ``""`` default leaves ``primers.json`` as-is
    (it becomes the free-standing default), matching
    `_restore_primers_from_active_primer_collection`, which also no-ops for
    ``""``. If the mirror write fails the active pointer is reverted, so the
    live library and the active-name setting never disagree."""
    name = payload.get("name")
    if name is None:
        return ({"error": "missing 'name' (use \"\" for default)"}, 400)
    if not isinstance(name, str):
        return ({"error": "'name' must be a string"}, 400)
    name = name.strip()
    # Validate + swap under ONE lock so a concurrent create-primer can't
    # interleave between the pointer move and the mirror re-point.
    with _state._cache_lock:
        colls = _load_primer_collections()
        valid_names = {(c.get("name") or "") for c in colls}
        valid_names.add("")  # default sentinel
        if name not in valid_names:
            return ({"error": (
                f"unknown primer collection {name!r}; "
                f"valid: {sorted(n for n in valid_names if n)}"
            )}, 404)
        prev = _get_active_primer_collection_name() or ""
        if (err := _agent_save_or_500(
                lambda: _set_active_primer_collection_name(name or None),
                "primer_collections")) is not None:
            return err
        if name:
            # Named collection → swap the live mirror to its stored primers.
            coll = next((c for c in colls
                         if (c.get("name") or "") == name), None)
            primers = [dict(p) for p in ((coll or {}).get("primers") or [])
                       if isinstance(p, dict)]
            if (err := _agent_save_or_500(
                    lambda: _safe_save_json_mirror(
                        _state._PRIMERS_FILE, primers, "Primer library"),
                    "primer library")) is not None:
                # Mirror failed AFTER the pointer moved — revert so the live
                # library and the active-name setting stay consistent.
                _set_active_primer_collection_name(prev or None)
                return err
            _state._primers_cache = None
    return {"ok": True, "active": name}


@_agent_endpoint("create-primer-collection", write=True)
def _h_create_primer_collection(app, payload):
    """Create a new (empty) named primer collection. Body:
    ``{name, description?}``.

    The plasmid-side parallel to ``create-collection``. Primer
    collections were previously creatable only through the GUI, so an
    agent could ``set-active-primer-collection`` to *switch* but had no
    way to *make* one (set-active 404s on an unknown name, and
    ``create-primer`` files into whatever collection is already active).
    To put a project's primers in their own collection: create it here,
    switch with ``set-active-primer-collection {name}``, then
    ``create-primer``.

    A name already in use (case-insensitive) returns 409; the reserved
    default collection is the empty string and can't be re-created. The
    new collection is written as ``{name, description, primers: [],
    saved}`` — the same envelope the GUI writes."""
    name = _normalize_collection_name(payload.get("name"))
    if name is None:
        return ({"error": "missing or invalid 'name'"}, 400)
    desc = _sanitize_note(payload.get("description"), max_len=2000)
    # RMW under the cache lock: a concurrent GUI/agent create of the
    # same name can't clobber via last-writer-wins (mirrors
    # create-collection sweep #11). The case-insensitive check also
    # guards the reserved "" default, which `_normalize_collection_name`
    # already rejects as blank.
    with _state._cache_lock:
        colls = _load_primer_collections()
        existing = {(c.get("name") or "").strip().casefold() for c in colls}
        if name.casefold() in existing:
            return ({"error":
                      f"primer collection {name!r} already exists"}, 409)
        colls.append({
            "name":        name,
            "description": desc,
            "primers":     [],
            "saved":       _date.today().isoformat(),
        })
        if (err := _agent_save_or_500(
                lambda: _save_primer_collections(colls),
                "primer_collections")) is not None:
            return err
    return {"ok": True, "name": name, "n_primers": 0,
            "ignored": _agent_ignored_keys(payload, {"name", "description"})}


@_agent_endpoint("rename-primer-collection", write=True)
def _h_rename_primer_collection(app, payload):
    """Rename a primer collection. Body: ``{name, new_name}`` (both
    named collections; the reserved default is the empty string and
    can't be renamed). Matching is case-insensitive; 404 if ``name``
    is unknown, 409 if ``new_name`` collides with another collection.
    If the renamed collection is the active one the active-pointer
    follows, so later ``create-primer`` mirrors land in the right slot.
    The GUI parallel is ``_rename_primer_collection`` (active-only); the
    agent can rename any collection by name."""
    old_in = payload.get("name")
    if not isinstance(old_in, str) or not old_in.strip():
        return ({"error": "missing 'name'"}, 400)
    old = old_in.strip()
    new = _normalize_collection_name(payload.get("new_name"))
    if new is None:
        return ({"error": "missing or invalid 'new_name'"}, 400)
    # RMW under the cache lock — a concurrent GUI/agent rename of the same
    # collection can't clobber via last-writer-wins (mirrors create-).
    with _state._cache_lock:
        colls = _load_primer_collections()
        matched = next(
            (c.get("name") or "" for c in colls
             if (c.get("name") or "").strip().casefold() == old.casefold()),
            None)
        if matched is None:
            valid = sorted(c.get("name") or "" for c in colls if c.get("name"))
            return ({"error":
                      f"unknown primer collection {old!r}; valid: {valid}"},
                    404)
        if (new.casefold() != matched.casefold()
                and any((c.get("name") or "").strip().casefold()
                        == new.casefold() for c in colls)):
            return ({"error":
                      f"primer collection {new!r} already exists"}, 409)
        for c in colls:
            if (c.get("name") or "") == matched:
                c["name"] = new
                break
        if (err := _agent_save_or_500(
                lambda: _save_primer_collections(colls),
                "primer_collections")) is not None:
            return err
        if (_get_active_primer_collection_name() or "") == matched:
            _set_active_primer_collection_name(new)
            _state._settings_flush_sync_hook()
    _log_event("primer_collection.rename", old=matched, new=new, via="agent")
    return {"ok": True, "name": new, "renamed_from": matched}


@_agent_endpoint("delete-primer-collection", write=True)
def _h_delete_primer_collection(app, payload):
    """Delete a named primer collection + its primers. Body: ``{name}``.
    Refuses to delete the last remaining named collection (matches the
    GUI last-collection guard); the reserved default ("") can't be
    deleted. Matching is case-insensitive. If the deleted collection was
    active, the first remaining collection is promoted to active and the
    live ``primers.json`` mirror is swapped to its primers ([INV-83]).
    The GUI parallel is ``_delete_primer_collection`` (active-only); the
    agent can delete any non-active collection without touching the
    mirror."""
    name_in = payload.get("name")
    if not isinstance(name_in, str) or not name_in.strip():
        return ({"error": "missing 'name'"}, 400)
    name = name_in.strip()
    with _state._cache_lock:
        colls = _load_primer_collections()
        matched = next(
            (c.get("name") or "" for c in colls
             if (c.get("name") or "").strip().casefold() == name.casefold()),
            None)
        if matched is None:
            valid = sorted(c.get("name") or "" for c in colls if c.get("name"))
            return ({"error":
                      f"unknown primer collection {name!r}; valid: {valid}"},
                    404)
        if len(colls) <= 1:
            return ({"error":
                      "cannot delete the last remaining primer collection"},
                    409)
        n = len(next((c.get("primers") or [] for c in colls
                      if (c.get("name") or "") == matched), []))
        remaining = [c for c in colls if (c.get("name") or "") != matched]
        if (err := _agent_save_or_500(
                lambda: _save_primer_collections(remaining),
                "primer_collections")) is not None:
            return err
        # Auto-promote + mirror-swap stays under the same lock hold (RLock,
        # so the inner load/save re-enter safely) — mirrors the
        # delete-experiment-project auto-promote.
        promoted = ""
        if (_get_active_primer_collection_name() or "") == matched:
            promoted = remaining[0].get("name") or ""
            _set_active_primer_collection_name(promoted or None)
            _state._settings_flush_sync_hook()
            primers = [dict(p) for p in (remaining[0].get("primers") or [])
                       if isinstance(p, dict)]
            try:
                # [INV-83] Mirror-swap to the promoted collection's primers.
                _safe_save_json_mirror(
                    _state._PRIMERS_FILE, primers, "Primer library")
            except (OSError, RuntimeError) as exc:
                # Source-of-truth write already succeeded; mirror self-heals
                # on next reload. Log rather than swallow silently.
                _log.warning(
                    "agent primer-collection-delete: mirror-swap to promoted "
                    "collection failed (self-heals on next reload): %s", exc)
            _state._primers_cache = None
    _log_event("primer_collection.delete", name=matched,
                n_primers=n, via="agent", promoted=promoted)
    return {"ok": True, "deleted": matched, "n_primers": n,
            "promoted": promoted,
            "remaining": [c.get("name") or "" for c in remaining]}


@_agent_endpoint("move-primer", write=True)
def _h_move_primer(app, payload):
    """Move a primer between collections WITHOUT a create/delete round-trip
    (SC-K — that round-trip lost primers: `create-primer` rejected the
    duplicate, then the follow-up `delete-primer` removed the only copy).
    Body: ``{name|id|sequence, to, from?}``. ``to``/``from`` are collection
    names (``""`` = the default / live library).

    The live library (``primers.json``) IS the ACTIVE collection's content —
    `_save_primers` mirrors one into the other — so ``""`` and the active
    collection's own NAME address the SAME store. Both are canonicalised to
    ``""`` and the move routes through `_save_primers` (the only writer that
    keeps ``primers.json`` and the stored mirror consistent). That makes
    moving a primer OUT of (or INTO) the active collection safe: before this
    there was no such path — ``""`` was refused while a named collection was
    active, and passing the active collection's name edited its stored mirror
    directly, a change the next `_save_primers` silently reverted (leaving the
    primer duplicated).

    Resolves the primer in ``from`` (or searches every DISTINCT store when
    ``from`` is omitted; 409 if it lives in more than one), removes it there
    and adds it to ``to`` — all under ONE lock, so the primer is never
    momentarily absent. 404 if the primer or a named collection is missing;
    409 on an ambiguous name, a sequence already present in ``to``, or when
    ``to`` is the store the primer already lives in."""
    if (derr := _agent_reject_dangerous_unknowns(
            payload, {"name", "id", "sequence", "to", "from"})) is not None:
        return derr
    to_in = payload.get("to")
    if not isinstance(to_in, str):
        return ({"error":
                  "missing 'to' (collection name; \"\" = default)"}, 400)
    to_raw = to_in.strip()
    from_in = payload.get("from")
    seq = payload.get("sequence")
    name = payload.get("name") or payload.get("id")
    has_seq = isinstance(seq, str) and bool(seq.strip())
    has_name = isinstance(name, str) and bool(name.strip())
    if not (has_seq or has_name):
        return ({"error": "missing 'name', 'id', or 'sequence'"}, 400)

    def _match(p):
        if has_seq:
            return (p.get("sequence") or "").upper() == seq.strip().upper()
        return (p.get("name") or "") == name.strip()

    with _state._cache_lock:
        colls = _load_primer_collections()
        by_name = {(c.get("name") or ""): c for c in colls}
        default_list = _load_primers()
        active = _get_active_primer_collection_name() or ""

        # The live library (primers.json) IS the active collection's content
        # (`_save_primers` mirrors one into the other, [INV-50] / pitfall #10),
        # so "" and the active collection's own NAME are the SAME store.
        # Canonicalise any reference to the active collection down to "" so the
        # move routes through `_save_primers` — the ONLY writer that keeps
        # primers.json and the stored mirror consistent. Editing the active
        # collection's stored `primers` list directly (the pre-fix bug) desynced
        # it from primers.json, and the next `_save_primers` silently reverted.
        def _canon(cn):
            return "" if (cn == "" or (active and cn == active)) else cn

        def _store(cn):        # cn is a CANONICAL id ("" or a non-active name)
            if cn == "":
                return default_list
            c = by_name.get(cn)
            return c.get("primers") if c is not None else None

        # Validate + canonicalise the target.
        if to_raw and to_raw not in by_name:
            return ({"error":
                      f"no primer collection named {to_raw!r}; create it "
                      f"with create-primer-collection"}, 404)
        to_name = _canon(to_raw)
        # Resolve + canonicalise the source.
        if isinstance(from_in, str) and from_in.strip():
            from_raw = from_in.strip()
            if from_raw and from_raw not in by_name:
                return ({"error":
                          f"no primer collection named {from_raw!r}"}, 404)
            from_name = _canon(from_raw)
        elif from_in in (None, ""):
            # DISTINCT stores only: the active collection is already covered by
            # "" (same data), so exclude its name to avoid a spurious
            # "in multiple collections" 409.
            canon_stores = [""] + [n for n in by_name if n and n != active]
            cands = [cn for cn in canon_stores
                     if any(_match(p) for p in (_store(cn) or []))]
            if not cands:
                return ({"error": "primer not found in any collection"}, 404)
            if len(cands) > 1:
                return ({"error":
                          "primer is in multiple collections — pass 'from'",
                          "collections": cands}, 409)
            from_name = cands[0]
        else:
            return ({"error": "'from' must be a string"}, 400)
        if from_name == to_name:
            return ({"error":
                      f"primer already in {to_name or 'default'!r}"}, 409)
        src = _store(from_name)
        assert src is not None             # from_name validated above
        idxs = [i for i, p in enumerate(src) if _match(p)]
        if not idxs:
            return ({"error":
                      f"primer not found in {from_name or 'default'!r}"}, 404)
        if len(idxs) > 1:
            return ({"error":
                      f"ambiguous — {len(idxs)} primers match in "
                      f"{from_name or 'default'!r}; pass 'sequence'"}, 409)
        primer = dict(src[idxs[0]])
        tgt = _store(to_name)
        assert tgt is not None             # to_name validated above
        tseq = (primer.get("sequence") or "").upper()
        if any((p.get("sequence") or "").upper() == tseq for p in tgt):
            return ({"error":
                      f"a primer with that sequence already exists in "
                      f"{to_name or 'default'!r}"}, 409)
        # The move — under the single lock, so no momentary gap. Persist the
        # store the primer is ADDED to (`to`) BEFORE the store it's removed
        # from, so a crash between the two atomic writes DUPLICATES the primer
        # (recoverable) rather than LOSING it (the data-safety contract is
        # "never lose, may duplicate"). Mirror correctness holds because the
        # ordering keeps `_save_primers`' active re-mirror as the last word on
        # the active entry.
        del src[idxs[0]]
        tgt.insert(0, primer)
        if to_name == "":
            # Added to the live/active library: make primers.json + the active
            # stored mirror durable FIRST (both via `_save_primers`), then the
            # removal from the named source.
            if (e := _agent_save_or_500(
                    lambda: _save_primers(default_list), "primers")) is not None:
                return e
            if from_name != "":
                # `_save_primers` just re-mirrored the active collection into the
                # collections file from the (deduped) saved library. Refresh
                # `colls`' active entry to match, so the removal-save below keeps
                # it rather than reverting to the stale pre-move snapshot.
                if active and active in by_name:
                    by_name[active]["primers"] = _load_primers()
                if (e := _agent_save_or_500(
                        lambda: _save_primer_collections(colls),
                        "primer_collections")) is not None:
                    return e
        else:
            # Added to a NAMED collection: make that store durable FIRST, then
            # (only when the source is the live library) persist its removal.
            # `_save_primers` runs LAST here, so its re-mirror leaves the active
            # entry correct even though `colls` carried a stale copy of it.
            if (e := _agent_save_or_500(
                    lambda: _save_primer_collections(colls),
                    "primer_collections")) is not None:
                return e
            if from_name == "":
                if (e := _agent_save_or_500(
                        lambda: _save_primers(default_list), "primers")) is not None:
                    return e
    _log_event("primer.moved", name=primer.get("name", ""),
               source=from_name, dest=to_name, via="agent")
    return {"ok": True, "name": primer.get("name"),
            "sequence": primer.get("sequence"),
            "from": from_name, "to": to_name,
            "ignored": _agent_ignored_keys(
                payload, {"name", "id", "sequence", "to", "from"})}


@_agent_endpoint("list-hmm-databases")
def _h_list_hmm_databases(app, payload):
    """Every registered HMM database. Each item: ``{id, name, ready,
    builtin}`` (`ready` = pressed/usable on disk). `active` is the
    picked id used by hmmscan."""
    out = []
    for e in _load_hmm_db_catalog():
        eid = e.get("id", "")
        out.append({
            "id":      eid,
            "name":    e.get("name", eid),
            "ready":   _hmm_db_pressed(eid),
            "builtin": bool(e.get("builtin")),
        })
    return {"ok": True, "hmm_databases": out,
            "active": (_get_setting("hmm_db_active_id", "") or "")}


@_agent_endpoint("get-hmm-database")
def _h_get_hmm_database(app, payload):
    """Detail for one HMM database by `id`."""
    eid = payload.get("id")
    if not isinstance(eid, str) or not eid.strip():
        return ({"error": "missing 'id'"}, 400)
    eid = eid.strip()
    for e in _load_hmm_db_catalog():
        if e.get("id") == eid:
            return {"ok": True, "id": eid, "name": e.get("name", eid),
                    "ready": _hmm_db_pressed(eid),
                    "builtin": bool(e.get("builtin")),
                    "url": e.get("url", ""),
                    "version": e.get("version", "")}
    return ({"error": "hmm database not found"}, 404)


@_agent_endpoint("delete-hmm-database", write=True)
def _h_delete_hmm_database(app, payload):
    """Delete the DOWNLOADED files for an HMM database by `id`
    (un-downloads it; the catalog entry stays so it can be re-fetched).
    Files removal is L2-gated via `_delete_hmm_db_files`."""
    eid = payload.get("id")
    if not isinstance(eid, str) or not eid.strip():
        return ({"error": "missing 'id'"}, 400)
    eid = eid.strip()
    if eid not in {e.get("id") for e in _load_hmm_db_catalog()}:
        return ({"error": "hmm database not found"}, 404)
    try:
        removed = _delete_hmm_db_files(eid)
    except (OSError, RuntimeError) as exc:
        _log.exception("agent delete-hmm-database failed")
        return ({"error": f"delete failed: {_scrub_path(str(exc))}"}, 500)
    _log_event("hmm_db.delete", id=eid, files_removed=removed, via="agent")
    return {"ok": True, "id": eid, "files_removed": removed}


@_agent_endpoint("add-hmm-database", write=True)
def _h_add_hmm_database(app, payload):
    """Register a CUSTOM HMM database in the catalog — the headless
    equivalent of the GUI catalog modal's "Add" form
    (`HmmDbAddEditModal`). Body: ``{name, url, version_url?,
    description?}``.

    ``url`` must be a well-formed ``http(s)://`` link (no whitespace)
    pointing at a gzipped HMMER3 ``.hmm.gz``; ``http`` is downloaded only
    if the `hmm_db_allow_http` setting is on (enforced later, at download
    time). The id is derived from the name exactly as the form does
    (``name.replace(" ", "_").lower()`` → `_sanitize_hmm_db_id`); a
    colliding id (incl. a builtin) or a taken display name is refused
    (409). This REGISTERS the entry only — call `download-hmm-database`
    to fetch + index the files. Catalog write is L2-gated via
    `_save_hmm_db_catalog`."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing or empty 'name'"}, 400)
    name = name.strip()
    url = payload.get("url")
    if not isinstance(url, str):
        return ({"error": "missing 'url'"}, 400)
    url = url.strip()
    if _sanitize_hmm_db_url(url) is None:
        return ({"error": "invalid 'url' — must be a well-formed "
                 "http(s):// link with no whitespace"}, 400)
    vurl = payload.get("version_url", "")
    if not isinstance(vurl, str):
        return ({"error": "'version_url' must be a string"}, 400)
    vurl = vurl.strip()
    if vurl and _sanitize_hmm_db_url(vurl) is None:
        return ({"error": "invalid 'version_url' — must be http(s):// "
                 "or omitted"}, 400)
    desc = payload.get("description", "")
    if not isinstance(desc, str):
        return ({"error": "'description' must be a string"}, 400)
    entry_id = _sanitize_hmm_db_id(name.replace(" ", "_").lower())
    if entry_id is None:
        return ({"error": "name can't be reduced to a valid id "
                 "(needs letters / digits / `_` / `-`)"}, 400)
    if _find_hmm_db_entry(entry_id) is not None:
        return ({"error": f"an entry with id '{entry_id}' already exists"},
                409)
    if _hmm_db_name_taken(name):
        return ({"error": f"display name '{name}' is taken"}, 409)
    normalised = _normalise_hmm_db_entry({
        "id":          entry_id,
        "name":        name[:200],
        "url":         url,
        "version_url": vurl,
        "format":      "hmm-gz",
        "builtin":     False,
        "description": desc.strip(),
    })
    if normalised is None:
        return ({"error": "entry failed validation"}, 400)
    catalog = _load_hmm_db_catalog()
    catalog.append(normalised)
    err = _agent_save_or_500(
        lambda: _save_hmm_db_catalog(catalog), "HMM database catalog")
    if err:
        return err
    _log_event("hmm_db.add", id=entry_id, via="agent")
    return {"ok": True, "id": entry_id,
            "name": normalised.get("name", entry_id)}


@_agent_endpoint("download-hmm-database", write=True)
def _h_download_hmm_database(app, payload):
    """Download + index an HMM database by `id` so hmmscan can use it —
    the headless equivalent of the GUI catalog modal's Download button.
    Body: ``{id}`` (a builtin like `pfam-a` / `ncbifam`, or a custom
    entry from `add-hmm-database`; see `list-hmm-databases`).

    Runs SYNCHRONOUSLY: the request blocks through download → decompress
    → hmmpress, so a builtin (Pfam-A ~300 MB, NCBIfam ~600 MB) can take
    minutes — size your client timeout accordingly. The cross-process
    `_HMM_DB_DOWNLOAD_INFLIGHT` slot returns 409 if the same DB is
    already downloading (e.g. a GUI download in flight). Files land in
    `<DATA_DIR>/hmm_databases/<id>/` exactly as a GUI download leaves
    them — on success the DB reports `ready: true` in
    `list-hmm-databases` and becomes selectable via
    `set-active-hmm-database`. Everything below the slot guard runs in
    the shared, L2-gated `_hmm_db_perform_download`."""
    eid = payload.get("id")
    if not isinstance(eid, str) or not eid.strip():
        return ({"error": "missing 'id'"}, 400)
    eid = eid.strip()
    entry = _find_hmm_db_entry(eid)
    if entry is None:
        return ({"error": "hmm database not found"}, 404)
    if not entry.get("url"):
        return ({"error": "catalog entry has no download URL"}, 400)
    if not _hmm_db_acquire_download_slot(eid):
        return ({"error": "a download for this database is already "
                 "running"}, 409)
    try:
        result = _hmm_db_perform_download(entry)
        _log_event("hmm_db.downloaded", id=eid,
                   n_profiles=result["n_profiles"],
                   pressed=result["pressed"], bytes=result["bytes"],
                   via="agent")
        return {"ok": True, "id": eid, "name": entry.get("name", eid),
                "n_profiles": result["n_profiles"],
                "pressed": result["pressed"], "version": result["version"],
                "sha256": result["sha256"], "bytes": result["bytes"],
                "ready": _hmm_db_pressed(eid)}
    except (OSError, ValueError, RuntimeError) as exc:
        _log.exception("agent download-hmm-database failed for %r", eid)
        _log_event("hmm_db.download.failed", id=eid,
                   error=str(exc)[:200], via="agent")
        return ({"error": f"download failed: {_scrub_path(str(exc))}"}, 500)
    except Exception as exc:
        _log.exception("agent download-hmm-database crashed for %r", eid)
        _log_event("hmm_db.download.crashed", id=eid,
                   exc_type=type(exc).__name__, via="agent")
        return ({"error": "download failed: unexpected error "
                 "(see splicecraft log)"}, 500)
    finally:
        _hmm_db_release_download_slot(eid)


# `list-plasmidsaurus-items` is the canonical name (every other list endpoint
# uses the `list-` prefix); `plasmidsaurus-items` stays as an alias for
# V1_GATE back-compat. Same handler, registered under both names.
@_agent_endpoint("list-plasmidsaurus-items")
@_agent_endpoint("plasmidsaurus-items")
def _h_plasmidsaurus_items(app, payload):
    """List your Plasmidsaurus items (sequencing orders), most-recent first,
    over the official OAuth2 REST API. Read-only; no body.

    Credentials resolve env-first (``PLASMIDSAURUS_CLIENT_ID`` /
    ``PLASMIDSAURUS_CLIENT_SECRET``) then the Settings values — returns 400
    if neither is set. Each item carries ``{code, status, order_name,
    product_name, quantity, done_date, gross, has_results}``.

    **Pick by ``has_results``, not by ``status``.** Roughly 40% of a real
    listing is shipping-label orders (``product_name:
    "ups_shipping_label"``), which are marked ``status: "complete"`` with a
    null ``done_date`` and sort to the TOP — but have no results zip, so
    feeding one to ``download-plasmidsaurus`` 404s. ``status`` also takes the
    value ``"canceled"``, likewise with nothing to download. ``has_results``
    folds all of that into one flag; ``downloadable`` counts them.

    The server truncates at 1000 items with no pagination available — if
    ``truncated`` is true, older orders exist that the API cannot reach.
    Network/credential failures surface as 502 (the secret is never echoed);
    a 429 means the account's 10-requests-per-minute budget is spent."""
    cid, sec = _plasmidsaurus_credentials()
    if not cid or not sec:
        return ({"error": "no Plasmidsaurus API credentials — set "
                 "PLASMIDSAURUS_CLIENT_ID / PLASMIDSAURUS_CLIENT_SECRET or "
                 "add them in Settings"}, 400)
    try:
        token = _plasmidsaurus_oauth_token(cid, sec)
        items = _plasmidsaurus_list_items(token)
    except (OSError, ValueError) as exc:
        _log.warning("agent plasmidsaurus-items failed: %s", exc)
        return ({"error": f"Plasmidsaurus API error: {_scrub_path(str(exc))}"},
                502)
    # `order_name` is the human-meaningful label the user typed when placing
    # the order; without it a caller only sees `plasmid_low_copy` + a date and
    # can't tell which run is which.
    fields = ("code", "status", "order_name", "product_name", "quantity",
              "done_date", "gross")
    out = []
    for it in items:
        if not isinstance(it, dict):
            continue
        row = {k: it.get(k) for k in fields}
        row["has_results"] = _plasmidsaurus_item_has_results(it)
        out.append(row)
    n_dl = sum(1 for r in out if r["has_results"])
    resp = {"ok": True, "count": len(out), "downloadable": n_dl, "items": out}
    # Truncation is about what the SERVER returned, so measure the raw list —
    # `out` drops any non-dict row, and a single junk entry in a full page
    # would otherwise put the count under the ceiling and silently retract
    # the warning on exactly the listing that needed it.
    if len(items) >= _PLASMIDSAURUS_ITEMS_LIMIT:
        resp["truncated"] = True
        resp["warning"] = _PLASMIDSAURUS_ITEMS_TRUNCATED_HINT
    _log_event("plasmidsaurus.items", count=len(out), downloadable=n_dl,
               via="agent")
    return resp


@_agent_endpoint("download-plasmidsaurus", write=True)
def _h_download_plasmidsaurus(app, payload):
    """Download a Plasmidsaurus item's ``results`` zip by item code and import
    every .gbk sample into the library. Body: ``{item_code, [kind]}``.

    ``item_code`` is the 6-char order code from ``plasmidsaurus-items`` /
    your dashboard (validated against ``^[A-Z0-9]{6}$`` before any request).
    ``kind`` defaults to ``"results"`` — only that archive carries the .gbk
    assemblies; ``reads`` / ``pod5`` are raw data with nothing to import, so
    they're refused here. Samples are ADDED as new library entries (never
    overwriting existing ones), tagged ``source: plasmidsaurus:<code>:<sample>``.

    Runs SYNCHRONOUSLY (download → parse → save) so a big run can take a
    while — size your client timeout accordingly. The zip lands in a temp
    dir and is deleted after import. Credentials resolve env-first then
    Settings (400 if unset); download/parse failures are 502, an empty/
    sampleless archive is 422."""
    if not isinstance(payload, dict):
        return ({"error": "expected a JSON object body"}, 400)
    code = _sanitize_plasmidsaurus_item_code(payload.get("item_code"))
    if code is None:
        return ({"error": "invalid or missing 'item_code' (expected 6 "
                 "characters, A-Z / 0-9)"}, 400)
    kind = payload.get("kind", "results")
    if kind not in _PLASMIDSAURUS_RESULT_KINDS:
        return ({"error": "'kind' must be one of "
                 f"{sorted(_PLASMIDSAURUS_RESULT_KINDS)}"}, 400)
    if kind != "results":
        return ({"error": "only kind='results' carries .gbk assemblies to "
                 "import; 'reads'/'pod5' are raw data"}, 400)
    cid, sec = _plasmidsaurus_credentials()
    if not cid or not sec:
        return ({"error": "no Plasmidsaurus API credentials — set "
                 "PLASMIDSAURUS_CLIENT_ID / PLASMIDSAURUS_CLIENT_SECRET or "
                 "add them in Settings"}, 400)
    import shutil
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix="splicecraft-ps-")
    try:
        try:
            zip_path = _plasmidsaurus_fetch_item_zip(
                code, tmpdir, kind="results",
                client_id=cid, client_secret=sec)
        except (OSError, ValueError) as exc:
            _log.warning("agent download-plasmidsaurus fetch failed "
                         "for %s: %s", code, exc)
            return ({"error": f"download failed: {_scrub_path(str(exc))}"},
                    502)
        try:
            entries, warnings = _plasmidsaurus_zip_to_entries(
                zip_path, run_id=code)
        except (OSError, ValueError) as exc:
            return ({"error": "could not read results zip: "
                     f"{_scrub_path(str(exc))}"}, 502)
        if not entries:
            return ({"error": "no .gbk samples found in the results zip",
                     "warnings": warnings}, 422)
        lib = _load_library()
        # Uniquify the imported ids against the existing library before
        # appending (mirrors the GUI fetch worker). Library `id` is the
        # canonical key: `_find_library_entry_by_id` returns the first match
        # and delete-by-id filters on `id != entry_id`, so a duplicate id
        # would make a later single-entry delete remove BOTH copies. A
        # re-import of the same run therefore suffix-bumps (`<id>_2`, …).
        seen = {e.get("id") for e in lib if isinstance(e, dict)}
        for e in entries:
            base = e["id"]
            n = 1
            while e["id"] in seen:
                n += 1
                e["id"] = f"{base}_{n}"
            seen.add(e["id"])
            lib.append(e)
        # Sync mirror (NOT async_sync): the agent is headless, so deferring the
        # active-collection mirror to the background worker buys nothing and
        # opens a crash window where the freshly-imported plasmids are in
        # `plasmid_library.json` but not yet in the active collection — a
        # crash + restart would then rebuild the live file from the active
        # collection (`_restore_library_from_active_collection`) and LOSE the
        # imports. Blocking the import response on the mirror is fine here.
        if (err := _agent_save_or_500(
                lambda: _save_library(lib), "library")):
            return err
        added = [{"name": e["name"], "size": e["size"],
                  "n_feats": e["n_feats"]} for e in entries]
        _log_event("plasmidsaurus.import", item_code=code,
                   n_added=len(added), n_warn=len(warnings), via="agent")
        return {"ok": True, "item_code": code, "n_added": len(added),
                "added": added, "warnings": warnings}
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _settings_validator_bool(value):
    if isinstance(value, bool):
        return value, None
    return None, "must be a boolean"


def _settings_validator_int_range(lo: int, hi: int):
    def _v(value):
        # Sweep #9 (2026-05-19): reject bool BEFORE coercing.
        # `_coerce_int(True)` returns 1 (Python truthiness), but
        # invariant #41 documents strict bool-vs-int separation
        # for settings — `_validate_settings` enforces it on the
        # disk-load path, this validator now does the same for
        # the agent-API path. Without this, `set-setting
        # min_primer_binding=true` silently succeeded with value 1.
        if isinstance(value, bool):
            return None, "must be int (got bool)"
        result = _coerce_int(value, name="value")
        if isinstance(result, str):
            return None, result
        if not (lo <= result <= hi):
            return None, f"must be in [{lo}, {hi}]"
        return result, None
    return _v


def _settings_validator_min_len_4_or_6(value):
    result = _coerce_int(value, name="value")
    if isinstance(result, str):
        return None, result
    if result not in (4, 6):
        return None, "must be 4 or 6"
    return result, None


def _settings_validator_custom_enzymes_csv(value):
    """Comma-separated list of enzyme names, each in the combined
    catalog (built-in NEB ∪ user-added custom enzymes). Whitespace
    tolerated; duplicates silently collapsed.

    Defensive parse — we'd rather drop an unknown name than reject the
    whole list, so a typo or a HF-variant rename doesn't strand the
    user. Returns the canonicalised CSV form (sorted, deduped, valid
    names only)."""
    if value is None or value == "":
        return "", None
    if not isinstance(value, str):
        return None, "must be a comma-separated string of enzyme names"
    raw_names = [
        s.strip() for s in value.replace(";", ",").split(",")
        if s and s.strip()
    ]
    catalog = _state._all_enzymes_hook()
    valid: list[str] = []
    for nm in raw_names:
        if nm in catalog and nm not in valid:
            valid.append(nm)
    return ",".join(sorted(valid)), None


def _settings_validator_collection_name(value):
    """Empty string / null is allowed to clear the active collection.
    Otherwise route through `_normalize_collection_name`."""
    if value is None or value == "":
        return "", None
    norm = _normalize_collection_name(value)
    if norm is None:
        return None, "invalid collection name"
    return norm, None


def _settings_validator_grammar_id(value):
    if not isinstance(value, str) or not value:
        return None, "must be a non-empty string"
    # Cap defensively; grammar IDs are short identifiers.
    if len(value) > 100:
        return None, "too long (max 100 chars)"
    return value, None


# Allowlist of user-facing toggle settings the agent may read / write.
# Infrastructure caches (`last_known_latest`, `last_seen_version`,
# `last_update_check_ts`, `hmm_db_path`) are deliberately excluded —
# those are session bookkeeping, not user preferences.
_AGENT_SETTINGS_ALLOWLIST: "dict[str, tuple]" = {
    # key: (validator, default)
    "show_feature_tooltips": (_settings_validator_bool,                  False),
    "click_debug":           (_settings_validator_bool,                  False),
    "check_updates":         (_settings_validator_bool,                  True),
    "show_restr":            (_settings_validator_bool,                  False),
    "restr_unique_only":     (_settings_validator_bool,                  True),
    "show_connectors":       (_settings_validator_bool,                  False),
    "restr_min_len":         (_settings_validator_min_len_4_or_6,        6),
    "restr_custom_enzymes":  (_settings_validator_custom_enzymes_csv,    ""),
    "restr_use_custom_list": (_settings_validator_bool,                  False),
    "min_primer_binding":    (_settings_validator_int_range(1, 60),      15),
    "sanger_min_phred":      (_settings_validator_int_range(0, 93),      20),
    "active_collection":     (_settings_validator_collection_name,       ""),
    "active_enzyme_collection": (_settings_validator_collection_name,     ""),
    "active_grammar":        (_settings_validator_grammar_id,            "gb_l0"),
    "constructor_filter_by_grammar": (_settings_validator_bool,           True),
}


@_agent_endpoint("get-settings")
def _h_get_settings(app, payload):
    """Return every allowlisted user-toggle setting with its current
    value, default, and validator hint. Useful for an agent that
    wants to inspect which toggles are exposed before writing."""
    out = {}
    for key, (_validator, default) in _AGENT_SETTINGS_ALLOWLIST.items():
        out[key] = {
            "value":   _get_setting(key, default),
            "default": default,
        }
    return {"ok": True, "settings": out}


@_agent_endpoint("set-setting", write=True)
def _h_set_setting(app, payload):
    """Persist a single user-toggle setting. Body: ``{key, value}``.
    `key` must be in the allowlist (see `get-settings`); `value` is
    type-checked + range-checked per the key's validator. Persists
    via `_set_setting`. Live in-memory app state mirrors are NOT
    refreshed here — the change takes effect on the next session
    restart for most toggles. (The GUI re-applies certain toggles
    immediately via dedicated action methods; mirroring that
    semantically across the agent surface would require dispatching
    each key to a UI-thread handler. Out of scope for v1.)

    Marked ``write=True`` so it requires the agent-API bearer token
    same as the other mutating endpoints (closes inconsistency where
    settings.json could be mutated token-free)."""
    key = payload.get("key")
    if not isinstance(key, str) or not key:
        return ({"error": "missing or non-string 'key'"}, 400)
    if key not in _AGENT_SETTINGS_ALLOWLIST:
        return ({"error": f"unknown setting {key!r}",
                  "available": list(_AGENT_SETTINGS_ALLOWLIST)}, 400)
    validator, _default = _AGENT_SETTINGS_ALLOWLIST[key]
    cleaned, err = validator(payload.get("value"))
    if err is not None:
        return ({"error": f"{key!r}: {err}"}, 400)
    _set_setting(key, cleaned)
    return {"ok": True, "key": key, "value": cleaned}


@_agent_endpoint("get-feature-colors")
def _h_get_feature_colors(app, payload):
    """Return the user's custom feature-type → colour overrides. Body: `{}`.

    Each entry maps a feature `type` (e.g. ``"CDS"``, ``"promoter"``) to a
    hex colour (``#RGB`` / ``#RRGGBB``). Types NOT listed fall back to the
    built-in palette — this endpoint reports only the user's overrides (the
    same set ``set-feature-color`` writes), so an empty map means "all
    types use their defaults"."""
    return {"ok": True, "feature_colors": _load_feature_colors()}


@_agent_endpoint("set-feature-color", write=True)
def _h_set_feature_color(app, payload):
    """Set or clear the render colour for a feature type. Body:
    ``{feature_type, color}`` — `color` is a hex string (``#1f77b4`` or the
    short ``#abc`` form). Pass ``{feature_type, clear: true}`` (or
    ``color: ""``/``null``) to REMOVE the override so the type falls back to
    the built-in palette. Persisted to ``feature_colors.json``; changes how
    that type renders on every map. Clearing a type that has no override
    returns 404."""
    # Validate the RAW value first: `_sanitize_feat_type` defaults None/""
    # to "misc_feature", which would silently set a colour for the wrong
    # type on a missing key.
    ft_raw = payload.get("feature_type")
    if not isinstance(ft_raw, str) or not ft_raw.strip():
        return ({"error": "missing or invalid 'feature_type'"}, 400)
    ft = _sanitize_feat_type(ft_raw)
    clear = bool(payload.get("clear")) or payload.get("color") in (None, "")
    color = None
    if not clear:
        color = _safe_color_for_write(payload.get("color"))
        if color is None:
            return ({"error":
                      "invalid 'color' (expected hex #RGB / #RRGGBB, "
                      "or a terminal palette ref color(0)-color(255))"}, 400)
    # RMW under the cache lock so a concurrent colour edit can't clobber.
    with _state._cache_lock:
        mapping = _load_feature_colors()
        if clear:
            if ft not in mapping:
                return ({"error":
                          f"no custom colour set for type {ft!r}"}, 404)
            del mapping[ft]
        else:
            assert color is not None       # validated above (else returned 400)
            mapping[ft] = color
        if (err := _agent_save_or_500(
                lambda: _save_feature_colors(mapping),
                "Feature colors")) is not None:
            return err
    return {"ok": True, "feature_type": ft, "color": color, "cleared": clear,
            "ignored": _agent_ignored_keys(
                payload, {"feature_type", "color", "clear"})}


def _parts_bin_entry_summary(p: dict) -> dict:
    """Compact view of a parts-bin row for list endpoints. Drops the
    `gb_text` blob so a 500-entry bin doesn't return 50+ MB; agents
    that need full text call `get-part`."""
    return {
        "name":     p.get("name", ""),
        "type":     p.get("type", ""),
        "level":    p.get("level", 0),
        "position": p.get("position", ""),
        "grammar":  p.get("grammar", ""),
        "oh5":      p.get("oh5", ""),
        "oh3":      p.get("oh3", ""),
        "size":     int(p.get("size") or len(str(p.get("sequence") or ""))),
    }


def _agent_parts_bin_list(bin_name):
    """SC-M: the parts belonging to bin `bin_name`. Returns ``(list, None)``
    or ``(None, (error_dict, status))``. ``bin_name`` empty / None → the ACTIVE
    bin (live); a named bin → ITS OWN stored ``parts`` (so a bin is a real
    partition for reads), 404 if it doesn't exist. Parallels
    `_agent_primer_collection_list`."""
    if bin_name in (None, ""):
        return list(_iter_parts_bin_readonly()), None
    if not isinstance(bin_name, str):
        return None, ({"error": "'bin' must be a string"}, 400)
    b = _find_parts_bin(bin_name.strip())
    if b is None:
        return None, ({"error":
                        f"no parts bin named {bin_name.strip()!r}"}, 404)
    return list(b.get("parts") or []), None


@_agent_endpoint("list-parts")
def _h_list_parts(app, payload):
    """List parts in a parts bin. Optional body:
    ``{bin?, grammar?: str, level?: int, position?: str}``. Pass ``{bin: "X"}``
    to list only bin X's OWN parts (SC-M); omit it for the active bin. Returns
    compact rows (no `gb_text`); call `get-part` for full content."""
    parts_src, berr = _agent_parts_bin_list(payload.get("bin"))
    if berr is not None:
        return berr
    assert parts_src is not None         # berr is None ⇒ parts_src is populated
    grammar = payload.get("grammar")
    level   = payload.get("level")
    position = payload.get("position")
    if level is not None:
        lvl = _coerce_int(level, name="level")
        if isinstance(lvl, str):
            return ({"error": lvl}, 400)
    else:
        lvl = None
    if grammar is not None and not isinstance(grammar, str):
        return ({"error": "'grammar' must be string"}, 400)
    if position is not None and not isinstance(position, str):
        return ({"error": "'position' must be string"}, 400)
    rows = []
    # `_parts_bin_entry_summary` builds a fresh dict from `p.get(...)` reads.
    for p in parts_src:
        if grammar and (p.get("grammar") or "") != grammar:
            continue
        if lvl is not None and int(p.get("level") or 0) != lvl:
            continue
        if position and (p.get("position") or "") != position:
            continue
        rows.append(_parts_bin_entry_summary(p))
    return {"ok": True, "parts": rows, "count": len(rows),
            "bin": (payload.get("bin") or "").strip()
                   if isinstance(payload.get("bin"), str) else ""}


@_agent_endpoint("get-part")
def _h_get_part(app, payload):
    """Fetch a single part by `name` (and optional `grammar` to
    disambiguate when two grammars carry a part with the same name).
    Returns the full entry including `gb_text` + `sequence` so the
    agent can parse it locally."""
    name = _sanitize_label(payload.get("name"), max_len=200)
    if not name:
        return ({"error": "missing 'name'"}, 400)
    grammar = payload.get("grammar")
    if grammar is not None and not isinstance(grammar, str):
        return ({"error": "'grammar' must be string"}, 400)
    for p in _load_parts_bin():
        if p.get("name") != name:
            continue
        if grammar and (p.get("grammar") or "") != grammar:
            continue
        return {"ok": True, "part": p}
    return ({"error": f"no part named {name!r}"
              + (f" in grammar {grammar!r}" if grammar else "")}, 404)


def _agent_part_dict(payload: dict) -> "dict | str":
    """Coerce + validate a parts-bin payload. Returns the cleaned
    dict or a string error message. Full parity with `PartEditModal`'s
    Edit-mode form: ``{name, sequence|gb_text, grammar?, type?, level?,
    position?, oh5?, oh3?, backbone?, marker?, fwd_primer?, rev_primer?,
    notes?}``. One of ``sequence``/``gb_text`` must be present so the part
    has content; everything else is optional. ``fwd_primer``/``rev_primer``
    are IUPAC DNA (≤1 kb) and their Tms are DERIVED on save (the modal
    never lets a user type a Tm directly)."""
    name = _sanitize_label(payload.get("name"), max_len=200)
    if not name:
        return "missing or non-string 'name'"
    grammar = payload.get("grammar", "")
    if grammar is not None and not isinstance(grammar, str):
        return "'grammar' must be a string"
    sequence = payload.get("sequence")
    gb_text  = payload.get("gb_text")
    if not isinstance(sequence, str) and not isinstance(gb_text, str):
        return "must provide 'sequence' or 'gb_text'"
    if sequence is not None and not isinstance(sequence, str):
        return "'sequence' must be a string"
    if gb_text is not None and not isinstance(gb_text, str):
        return "'gb_text' must be a string"
    clean_seq = ""
    if isinstance(sequence, str):
        clean_seq, seq_err = _sanitize_bases(sequence.upper())
        if seq_err:
            return f"'sequence' rejected: {seq_err}"
    elif isinstance(gb_text, str):
        try:
            rec = _gb_text_to_record(gb_text)
            clean_seq = str(rec.seq).upper()
        except Exception as exc:
            return f"gb_text parse failed: {exc}"
    if not clean_seq:
        return "no usable sequence after sanitisation"
    ptype = payload.get("type", "")
    if ptype is not None and not isinstance(ptype, str):
        return "'type' must be a string"
    level = payload.get("level", 0)
    if not isinstance(level, (int, float)) or isinstance(level, bool):
        return "'level' must be a number"
    position = payload.get("position", "")
    if position is not None and not isinstance(position, str):
        return "'position' must be a string"
    oh5 = payload.get("oh5", "")
    oh3 = payload.get("oh3", "")
    for label, val in (("oh5", oh5), ("oh3", oh3)):
        if val is not None and not isinstance(val, str):
            return f"{label!r} must be a string"
    # Vector + domestication-primer fields — the PartEditModal exposes these
    # in its Edit mode (backbone / selection marker / forward + reverse
    # domestication primers), so the agent must be able to set them too. The
    # Tms are DERIVED from the primer sequences (matching the modal, which
    # never lets the user type a Tm directly), so they're not request params.
    backbone = payload.get("backbone", "")
    marker   = payload.get("marker", "")
    for label, val in (("backbone", backbone), ("marker", marker)):
        if val is not None and not isinstance(val, str):
            return f"{label!r} must be a string"
    clean_fwd = clean_rev = ""
    for label, raw in (("fwd_primer", payload.get("fwd_primer", "")),
                        ("rev_primer", payload.get("rev_primer", ""))):
        if raw in (None, ""):
            continue
        if not isinstance(raw, str):
            return f"{label!r} must be a string"
        cleaned, perr = _sanitize_bases(raw.upper(), max_len=1000)
        if perr:
            return f"{label!r} rejected: {perr}"
        if label == "fwd_primer":
            clean_fwd = cleaned
        else:
            clean_rev = cleaned

    def _tm(seq: str) -> float:
        if not seq:
            return 0.0
        tm = _primer_tm_safe(seq)
        return round(float(tm), 1) if tm is not None else 0.0

    notes = _sanitize_note(payload.get("notes", ""))
    return {
        "name":       name,
        "grammar":    _sanitize_label(grammar),
        "type":       _sanitize_label(ptype),
        "level":      int(level),
        "position":   _sanitize_label(position),
        "oh5":        (oh5 or "").upper(),
        "oh3":        (oh3 or "").upper(),
        "backbone":   _sanitize_label(backbone, max_len=120),
        "marker":     _sanitize_label(marker, max_len=120),
        "fwd_primer": clean_fwd,
        "rev_primer": clean_rev,
        "fwd_tm":     _tm(clean_fwd),
        "rev_tm":     _tm(clean_rev),
        "sequence":   clean_seq,
        "size":       len(clean_seq),
        "gb_text":    gb_text if isinstance(gb_text, str) else "",
        "notes":      notes or "",
        "date":       _datetime.now().strftime("%Y-%m-%d"),
    }


@_agent_endpoint("create-part", write=True)
def _h_create_part(app, payload):
    """Add a part to a parts bin. Body: see `_agent_part_dict` for the part
    schema, plus optional ``{bin}``.

    Without ``bin`` the part lands in the ACTIVE bin. With ``bin`` it's filed
    into that NAMED bin — and if the bin doesn't exist the call FAILS with 404
    (create it first with `create-parts-bin`). SC-M: ``bin`` used to be
    silently ignored, so a part the caller believed was filed into a project
    bin actually landed in the active one. Returns 409 if a part with the same
    ``(name, grammar)`` already exists in the target bin."""
    p = _agent_part_dict(payload)
    if isinstance(p, str):
        return ({"error": p}, 400)
    bin_in = payload.get("bin")
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        if bin_in not in (None, ""):
            if not isinstance(bin_in, str):
                return ({"error": "'bin' must be a string"}, 400)
            bin_name = bin_in.strip()
            if bin_name != (_get_active_parts_bin_name() or ""):
                # File into a SPECIFIC non-active named bin: write its stored
                # parts list directly (the live parts_bin.json mirror belongs
                # to the ACTIVE bin and must stay untouched).
                bins = _load_parts_bin_collections()
                tgt = next((b for b in bins
                            if (b.get("name") or "") == bin_name), None)
                if tgt is None:
                    return ({"error":
                              f"no parts bin named {bin_name!r}; create it "
                              f"with create-parts-bin"}, 404)
                bucket = tgt.setdefault("parts", [])
                for e in bucket:
                    if (e.get("name") == p["name"]
                            and (e.get("grammar") or "") == p["grammar"]):
                        return ({"error":
                                  f"part {p['name']!r} already exists in bin "
                                  f"{bin_name!r} (grammar {p['grammar']!r}); "
                                  f"use update-part or rename."}, 409)
                bucket.append(p)
                if (err := _agent_save_or_500(
                        lambda: _save_parts_bin_collections(bins),
                        "parts_bin_collections")) is not None:
                    return err
                return {"ok": True, "name": p["name"],
                        "grammar": p["grammar"], "bin": bin_name}
            # bin == active → fall through to the mirrored active path.
        entries = _load_parts_bin()
        for e in entries:
            if (e.get("name") == p["name"]
                    and (e.get("grammar") or "") == p["grammar"]):
                return ({"error": (
                    f"part {p['name']!r} already exists in grammar "
                    f"{p['grammar']!r}; use update-part or rename."
                )}, 409)
        entries.append(p)
        if (err := _agent_save_or_500(
                lambda: _save_parts_bin(entries),
                "parts_bin")) is not None:
            return err
    return {"ok": True, "name": p["name"], "grammar": p["grammar"],
            "bin": _get_active_parts_bin_name() or ""}


@_agent_endpoint("make-l0-part-from-fragment", write=True)
def _h_make_l0_part_from_fragment(app, payload):
    """Turn a saved synthetic fragment into a Level-0 PART **and the L0
    plasmid it lives on**, by cloning it into the grammar's entry vector. Body:
    ``{fragment, part_type, [grammar], [name], [bin], [collection],
    [product_name], [save_plasmid]}``.

    ``fragment`` is a library entry id or name — a FRAG built by Synthesis →
    L0 Fragment (or ``design-synthesis-fragment``), i.e. one already wrapped
    in the grammar's Type IIS ends. ``part_type`` is the position the part
    will occupy (``CDS``, ``Promoter``, …); its fusion overhangs come from
    THAT grammar's position table, never from guessing at the sequence — which
    is what makes the two-tier nested layout come out right. ``grammar``
    defaults to the active one; ``bin`` routes exactly as in `create-part`.

    The clone is SIMULATED, not assembled by hand ([INV-127]): both the
    fragment and the entry vector are digested with the grammar's enzyme, the
    released insert is ligated into the surviving backbone, and the circle is
    closed. Anything that can't happen is a 400 naming the cause — a plain
    (unwrapped) fragment, a part type this grammar doesn't define (the error
    lists the ones it does), a part type whose overhangs aren't on this
    fragment, or a grammar with no entry vector assigned. 409 if a part of
    that ``(name, grammar)`` already exists.

    **Both artifacts are kept**, because the point of cloning a fragment is to
    end up holding the plasmid: the part goes to the bin, and the circular L0
    plasmid is SAVED to the library (active collection, or the one named by
    ``collection``) with its construction lineage attached, so History shows
    the insert and the entry vector it came from. That save is delegated to
    `domesticate-part` — the same engine the GUI Domesticator uses — so there
    is one definition of what an L0 clone is. The plasmid is named
    ``<part> (L0)`` so it can't be mistaken for the FRAGMENT it came from
    (which usually sits in the library under the part's own name);
    ``product_name`` overrides that, and a clash gets a ``-2`` suffix rather
    than overwriting anything. Pass ``save_plasmid: false`` to file only the
    part.

    The part's stored ``sequence`` is the body BETWEEN its overhangs, matching
    what `create-part` and the Domesticator file — `oh5`/`oh3` are held apart
    and re-added by the assembly, so the part re-clones through
    `_clone_part_into_entry_vector` into the very plasmid saved here.

    ``via`` is an internal label for the log line only (the Constructor button
    routes through this same handler and passes ``"ui"``); it defaults to
    ``"agent"`` and changes nothing about what is built."""
    if not isinstance(payload, dict):
        return ({"error": "expected a JSON object body"}, 400)
    known = {"fragment", "part_type", "grammar", "name", "bin", "collection",
             "product_name", "save_plasmid", "via"}
    # A mistyped `bin`/`collection` would silently route the part or the
    # plasmid somewhere else and still answer ok (SC-D).
    if (derr := _agent_reject_dangerous_unknowns(payload, known)) is not None:
        return derr
    # Guard ONCE, up front (`domesticate-parts`' pattern): the delegated save
    # guards too, and finding out the canvas is dirty AFTER the part is filed
    # would leave the two artifacts out of step.
    guard = _agent_dirty_guard(app, payload)
    if guard is not None:
        return guard
    frag_ref = payload.get("fragment")
    if not isinstance(frag_ref, str) or not frag_ref.strip():
        return ({"error": "missing or non-string 'fragment' "
                 "(library entry id or name)"}, 400)
    part_type = payload.get("part_type")
    if not isinstance(part_type, str) or not part_type.strip():
        return ({"error": "missing or non-string 'part_type'"}, 400)
    entry = (_find_library_entry_by_id(frag_ref)
             or _find_library_entry_by_name(frag_ref))
    if entry is None:
        return ({"error": f"no library entry {frag_ref!r}"}, 404)
    gid = payload.get("grammar")
    if gid in (None, ""):
        gid = _get_setting("active_grammar", "gb_l0") or "gb_l0"
    if not isinstance(gid, str):
        return ({"error": "'grammar' must be a string"}, 400)
    # `_all_grammars()` is keyed BY id — iterating it yields keys, not dicts.
    grammar = _all_grammars().get(gid)
    if not isinstance(grammar, dict):
        return ({"error": f"no grammar {gid!r}"}, 404)
    vec = _get_entry_vector(gid)
    if not isinstance(vec, dict) or not str(vec.get("gb_text") or ""):
        return ({"error": f"grammar {gid!r} has no entry vector assigned — "
                 f"set one first (Settings → Entry Vectors, or "
                 f"set-entry-vector)."}, 400)
    save_plasmid = payload.get("save_plasmid", True)
    if not isinstance(save_plasmid, bool):
        return ({"error": "'save_plasmid' must be true or false"}, 400)
    # Confirm the destination EXISTS before anything is filed. `domesticate-
    # part`'s `_agent_resolve_save_dest` is still the resolver that routes the
    # save — this only moves the "no such collection" failure ahead of the
    # part, so a typo can't leave a part filed with no plasmid beside it.
    coll_in = payload.get("collection")
    if coll_in not in (None, ""):
        # Same resolution as every other endpoint here: normalise, then exact
        # match on the readonly cache view (mirrors `_find_collection`).
        coll_name = _normalize_collection_name(coll_in)
        if coll_name is None:
            return ({"error": "invalid 'collection'"}, 400)
        if not any(isinstance(c, dict) and c.get("name") == coll_name
                   for c in _iter_collections_readonly()):
            return ({"error": f"no collection named {coll_name!r}"}, 404)
    try:
        frag_rec = _gb_text_to_record(str(entry.get("gb_text") or ""))
        vec_rec = _gb_text_to_record(str(vec.get("gb_text") or ""))
    except Exception as exc:                       # noqa: BLE001 — reported
        return ({"error": f"could not parse a record: {exc}"}, 400)
    # A CIRCULAR entry is not a synthesis fragment, and digesting one as if it
    # were linear is not an error — it happily "releases" whatever sits
    # between two sites, so pointing this at a plasmid (pUPD2 itself, say)
    # would file its stuffer as a part. Refuse by name instead.
    if (frag_rec.annotations.get("topology") or "").lower() == "circular":
        return ({"error":
                  f"{entry.get('name') or frag_ref!r} is a circular plasmid, "
                  f"not a linear synthesis fragment — this builds a part from "
                  f"a FRAG wrapped in {grammar.get('enzyme') or 'Type IIS'} "
                  f"ends (Synthesis → L0 Fragment)."}, 400)
    try:
        built = _l0_part_from_syn_fragment(
            str(frag_rec.seq), str(vec_rec.seq), grammar=grammar,
            part_type=part_type.strip(),
            name=str(payload.get("name") or entry.get("name") or "part"),
        )
    except ValueError as exc:
        return ({"error": str(exc)}, 400)
    # Delegate the actual filing so bin routing, the (name, grammar) dedup
    # and the save-failure path stay in ONE place (`create-part`).
    res = _h_create_part(app, {
        "name":     built["name"],
        "sequence": built["sequence"],
        "grammar":  built["grammar"],
        "type":     built["type"],
        "position": built["position"],
        "oh5":      built["oh5"],
        "oh3":      built["oh3"],
        "level":    0,
        "bin":      payload.get("bin"),
    })
    if isinstance(res, tuple):          # (error_dict, status) — pass through
        return res
    # ── The plasmid the part lives on ────────────────────────────────────
    # Delegated to `domesticate-part` (via the shared registry — this sibling
    # can't import the hub) rather than saving our own simulated record: that
    # handler already owns the two-tier nesting, the fail-loud "no stub"
    # policy, name/id dedup against the destination, the L0 lineage history
    # and collection routing. `_clone_part_into_entry_vector` under it
    # rebuilds the identical circle we simulated — pinned by
    # `test_the_filed_part_reclones_into_the_same_plasmid` — so delegating
    # costs a second digest and buys one definition of an L0 clone.
    saved: "dict | None" = None
    plasmid_error = ""
    if save_plasmid:
        dom = _state._AGENT_HANDLERS.get("domesticate-part")
        if dom is None:                 # registry gutted — should be unreachable
            plasmid_error = "the domesticate-part engine is unavailable"
        else:
            dres = dom[0](app, {
                "name":     built["name"],
                "sequence": built["sequence"],
                "oh5":      built["oh5"],
                "oh3":      built["oh3"],
                "grammar":  built["grammar"],
                "type":     built["type"],
                # Synthesised, not amplified: there are no domestication
                # primers to attach, and inventing some would put a PCR that
                # never happened into the construct's History.
                "collection":   payload.get("collection"),
                # Distinct from the FRAGMENT it was built from, which is
                # usually sitting in the library under exactly this name —
                # letting them collide would leave the user with "MyFrag" and
                # a dedup-suffixed "MyFrag-2" and no way to tell which is the
                # plasmid. Mirrors the Domesticator's "(domesticated)" tag.
                "product_name": (payload.get("product_name")
                                 or f"{built['name']} (L0)"),
                "force":        True,   # guarded once, at the top
            })
            if isinstance(dres, tuple):
                plasmid_error = str((dres[0] or {}).get(
                    "error", "the L0 plasmid could not be saved"))
                _log.warning("part.from_fragment: part %r filed but the "
                             "plasmid save failed: %s",
                             built["name"], plasmid_error)
            else:
                saved = dres
    via = payload.get("via")
    _log_event("part.from_fragment", grammar=gid, part_type=built["type"],
               nested=built["nested"], plasmid_saved=bool(saved),
               via="ui" if via == "ui" else "agent")
    out = dict(res)
    out.update({
        "part_type":   built["type"],
        "position":    built["position"],
        "oh5":         built["oh5"],
        "oh3":         built["oh3"],
        "insert_len":  built["insert_len"],
        "part_len":    len(built["sequence"]),
        "nested":      built["nested"],
        # The saved artifact's own length when there is one — the same circle
        # either way (`test_the_filed_part_reclones_into_the_same_plasmid`),
        # but a number that came off a real entry beats one that didn't.
        "plasmid_len": ((saved or {}).get("length")
                        or len(built["cloned_seq"])),
        "entry_vector": vec.get("name") or "",
        # The plasmid half of the result. `saved_id` is None when the part
        # filed but the plasmid didn't — reported rather than swallowed,
        # since the caller is now holding one artifact instead of two.
        "saved_name":  (saved or {}).get("saved_name"),
        "saved_id":    (saved or {}).get("saved_id"),
        "collection":  (saved or {}).get("collection"),
        "ignored":     _agent_ignored_keys(payload, known),
    })
    if plasmid_error:
        out["plasmid_error"] = plasmid_error
    return out


@_agent_endpoint("update-part", write=True)
def _h_update_part(app, payload):
    """Update an existing part. Body: `{name, grammar?, ...new fields}`.
    Lookup is by `(name, grammar)`; if grammar is omitted matches the
    first part by name across all grammars."""
    p_new = _agent_part_dict(payload)
    if isinstance(p_new, str):
        return ({"error": p_new}, 400)
    target_grammar = p_new["grammar"]
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        entries = _load_parts_bin()
        for i, e in enumerate(entries):
            if e.get("name") != p_new["name"]:
                continue
            if target_grammar and (e.get("grammar") or "") != target_grammar:
                continue
            # Preserve original date so update doesn't masquerade as add.
            p_new["date"] = e.get("date") or p_new["date"]
            entries[i] = p_new
            if (err := _agent_save_or_500(
                    lambda: _save_parts_bin(entries),
                    "parts_bin")) is not None:
                return err
            return {"ok": True, "name": p_new["name"],
                    "grammar": p_new["grammar"]}
    return ({"error": (
        f"no part {p_new['name']!r}"
        + (f" in grammar {target_grammar!r}" if target_grammar else "")
    )}, 404)


@_agent_endpoint("create-parts-bin", write=True)
def _h_create_parts_bin(app, payload):
    """Create a new (empty) named parts bin (SC-M) — the parts-side parallel
    to `create-collection` / `create-primer-collection`. Body:
    ``{name, description?}``.

    Parts bins were creatable only in the GUI, so an agent could
    `set-active-parts-bin` to *switch* but not *make* one. To put a project's
    parts in their own bin: create it here, then file with
    ``create-part {bin}`` (or switch to it with `set-active-parts-bin`). A name
    already in use (case-insensitive) returns 409."""
    name = _normalize_collection_name(payload.get("name"))
    if name is None:
        return ({"error": "missing or invalid 'name'"}, 400)
    desc = _sanitize_note(payload.get("description"), max_len=2000)
    with _state._cache_lock:
        bins = _load_parts_bin_collections()
        if any((b.get("name") or "").strip().casefold() == name.casefold()
               for b in bins):
            return ({"error": f"parts bin {name!r} already exists"}, 409)
        bins.append({
            "name":        name,
            "description": desc,
            "parts":       [],
            "saved":       _date.today().isoformat(),
        })
        if (err := _agent_save_or_500(
                lambda: _save_parts_bin_collections(bins),
                "parts_bin_collections")) is not None:
            return err
    return {"ok": True, "name": name, "n_parts": 0,
            "ignored": _agent_ignored_keys(payload, {"name", "description"})}


@_agent_endpoint("rename-parts-bin", write=True)
def _h_rename_parts_bin(app, payload):
    """Rename a parts bin. Body: ``{name, new_name}``. Matching is
    case-insensitive; 404 if ``name`` is unknown, 409 if ``new_name``
    collides with another bin. If the renamed bin is the active one the
    active-pointer follows so later ``create-part`` mirrors land in the
    right slot. The parts parallel to ``rename-primer-collection``."""
    old_in = payload.get("name")
    if not isinstance(old_in, str) or not old_in.strip():
        return ({"error": "missing 'name'"}, 400)
    old = old_in.strip()
    new = _normalize_collection_name(payload.get("new_name"))
    if new is None:
        return ({"error": "missing or invalid 'new_name'"}, 400)
    with _state._cache_lock:
        bins = _load_parts_bin_collections()
        matched = next(
            (b.get("name") or "" for b in bins
             if (b.get("name") or "").strip().casefold() == old.casefold()),
            None)
        if matched is None:
            valid = sorted(b.get("name") or "" for b in bins if b.get("name"))
            return ({"error": f"unknown parts bin {old!r}; valid: {valid}"},
                    404)
        if (new.casefold() != matched.casefold()
                and any((b.get("name") or "").strip().casefold()
                        == new.casefold() for b in bins)):
            return ({"error": f"parts bin {new!r} already exists"}, 409)
        for b in bins:
            if (b.get("name") or "") == matched:
                b["name"] = new
                break
        if (err := _agent_save_or_500(
                lambda: _save_parts_bin_collections(bins),
                "parts_bin_collections")) is not None:
            return err
        if _get_active_parts_bin_name() == matched:
            _set_active_parts_bin_name(new)
            _state._settings_flush_sync_hook()
    _log_event("parts_bin.rename", old=matched, new=new, via="agent")
    return {"ok": True, "name": new, "renamed_from": matched}


@_agent_endpoint("delete-parts-bin", write=True)
def _h_delete_parts_bin(app, payload):
    """Delete a named parts bin + its parts. Body: ``{name}``. Refuses to
    delete the last remaining bin (one must always exist). Matching is
    case-insensitive. If the deleted bin was active, the first remaining
    bin is promoted to active and the live ``parts_bin.json`` mirror is
    swapped to its parts ([INV-83]). The parts parallel to
    ``delete-primer-collection`` / ``delete-experiment-project``."""
    name_in = payload.get("name")
    if not isinstance(name_in, str) or not name_in.strip():
        return ({"error": "missing 'name'"}, 400)
    name = name_in.strip()
    with _state._cache_lock:
        bins = _load_parts_bin_collections()
        matched = next(
            (b.get("name") or "" for b in bins
             if (b.get("name") or "").strip().casefold() == name.casefold()),
            None)
        if matched is None:
            valid = sorted(b.get("name") or "" for b in bins if b.get("name"))
            return ({"error": f"unknown parts bin {name!r}; valid: {valid}"},
                    404)
        if len(bins) <= 1:
            return ({"error":
                      "cannot delete the last remaining parts bin"}, 409)
        n = len(next((b.get("parts") or [] for b in bins
                      if (b.get("name") or "") == matched), []))
        remaining = [b for b in bins if (b.get("name") or "") != matched]
        if (err := _agent_save_or_500(
                lambda: _save_parts_bin_collections(remaining),
                "parts_bin_collections")) is not None:
            return err
        promoted = ""
        if _get_active_parts_bin_name() == matched:
            promoted = remaining[0].get("name") or ""
            _set_active_parts_bin_name(promoted or None)
            _state._settings_flush_sync_hook()
            parts = [p for p in (remaining[0].get("parts") or [])
                     if isinstance(p, dict)]
            try:
                # [INV-83] Mirror-swap to the promoted bin's parts.
                _safe_save_json_mirror(
                    _state._PARTS_BIN_FILE, parts, "Parts bin")
            except (OSError, RuntimeError) as exc:
                _log.warning(
                    "agent parts-bin-delete: mirror-swap to promoted bin "
                    "failed (self-heals on next reload): %s", exc)
            _state._parts_bin_cache = None
            _state._clear_assembly_fragment_cache_hook()
    _log_event("parts_bin.delete", name=matched,
                n_parts=n, via="agent", promoted=promoted)
    return {"ok": True, "deleted": matched, "n_parts": n,
            "promoted": promoted,
            "remaining": [b.get("name") or "" for b in remaining]}


@_agent_endpoint("move-part", write=True)
def _h_move_part(app, payload):
    """Move a part between bins without a create/delete round-trip (SC-M, the
    parts analog of `move-primer`). Body: ``{name|id, to, from?, grammar?}``.

    Identifies the part by ``name`` (+ optional ``grammar`` to disambiguate),
    removes it from ``from`` (auto-found if omitted; 409 if it's in more than
    one bin) and adds it to ``to`` — atomically, under one lock. The active
    bin's live mirror is re-synced when touched, so a move into/out of the
    active bin stays consistent. 404 if the part or a bin is missing; 409 on an
    ambiguous name or a ``(name, grammar)`` collision in ``to``."""
    if (derr := _agent_reject_dangerous_unknowns(
            payload, {"name", "id", "to", "from", "grammar"})) is not None:
        return derr
    to_in = payload.get("to")
    if not isinstance(to_in, str) or not to_in.strip():
        return ({"error": "missing 'to' (a parts bin name)"}, 400)
    to_name = to_in.strip()
    name = _sanitize_label(payload.get("name") or payload.get("id"),
                            max_len=200)
    if not name:
        return ({"error": "missing 'name' or 'id'"}, 400)
    grammar = payload.get("grammar")
    if grammar is not None and not isinstance(grammar, str):
        return ({"error": "'grammar' must be a string"}, 400)
    from_in = payload.get("from")

    def _match(pp):
        if pp.get("name") != name:
            return False
        if grammar and (pp.get("grammar") or "") != grammar:
            return False
        return True

    with _state._cache_lock:
        bins = _load_parts_bin_collections()
        by_name = {(b.get("name") or ""): b for b in bins}
        active = _get_active_parts_bin_name() or ""
        if to_name not in by_name:
            return ({"error":
                      f"no parts bin named {to_name!r}; create it with "
                      f"create-parts-bin"}, 404)
        if isinstance(from_in, str) and from_in.strip():
            from_name = from_in.strip()
            if from_name not in by_name:
                return ({"error": f"no parts bin named {from_name!r}"}, 404)
        elif from_in in (None, ""):
            cands = [bn for bn, b in by_name.items()
                     if any(_match(pp) for pp in (b.get("parts") or []))]
            if not cands:
                return ({"error": "part not found in any bin"}, 404)
            if len(cands) > 1:
                return ({"error": "part is in multiple bins — pass 'from'",
                          "bins": cands}, 409)
            from_name = cands[0]
        else:
            return ({"error": "'from' must be a string"}, 400)
        if from_name == to_name:
            return ({"error": f"part already in {to_name!r}"}, 409)
        src = by_name[from_name].get("parts") or []
        idxs = [i for i, pp in enumerate(src) if _match(pp)]
        if not idxs:
            return ({"error": f"part {name!r} not found in {from_name!r}"}, 404)
        if len(idxs) > 1:
            return ({"error":
                      f"{name!r} matches {len(idxs)} parts in {from_name!r} — "
                      f"pass 'grammar' to pick one"}, 409)
        part = dict(src[idxs[0]])
        tgt = by_name[to_name].setdefault("parts", [])
        if any(pp.get("name") == part.get("name")
               and (pp.get("grammar") or "") == (part.get("grammar") or "")
               for pp in tgt):
            return ({"error":
                      f"a part {part.get('name')!r} (grammar "
                      f"{part.get('grammar')!r}) already exists in "
                      f"{to_name!r}"}, 409)
        del by_name[from_name]["parts"][idxs[0]]
        tgt.append(part)
        if (err := _agent_save_or_500(
                lambda: _save_parts_bin_collections(bins),
                "parts_bin_collections")) is not None:
            return err
        # Re-sync the live mirror if the ACTIVE bin's parts changed. Idempotent
        # on the canonical record (writes the same parts back), and keeps
        # parts_bin.json from drifting from the bin it mirrors.
        if active in (from_name, to_name) and active in by_name:
            ab = by_name[active].get("parts") or []
            if (err := _agent_save_or_500(
                    lambda: _save_parts_bin(ab), "parts_bin")) is not None:
                return err
    return {"ok": True, "name": name, "from": from_name, "to": to_name,
            "grammar": part.get("grammar"),
            "ignored": _agent_ignored_keys(
                payload, {"name", "id", "to", "from", "grammar"})}


def _agent_feature_dict(payload: dict) -> "dict | str":
    """Coerce + validate a feature-library payload. Schema:
    `{name, sequence, feature_type?, strand?, color?, notes?}`.
    Sequence is the canonical identity field for dedupe purposes.
    """
    name = _sanitize_label(payload.get("name"), max_len=200)
    if not name:
        return "missing or non-string 'name'"
    raw_seq = payload.get("sequence")
    if not isinstance(raw_seq, str) or not raw_seq.strip():
        return "missing or non-string 'sequence'"
    seq, seq_err = _sanitize_bases(raw_seq.upper())
    if seq_err:
        return f"'sequence' rejected: {seq_err}"
    if not seq:
        return "'sequence' empty after sanitisation"
    if len(seq) > 100_000:
        return "'sequence' too long (max 100 kb)"
    ftype = payload.get("feature_type", "misc_feature")
    if ftype is not None and not isinstance(ftype, str):
        return "'feature_type' must be a string"
    strand = payload.get("strand", 1)
    if isinstance(strand, bool) or not isinstance(strand, (int, float)):
        return "'strand' must be -1, 0, or 1"
    if int(strand) not in (-1, 0, 1):
        return "'strand' must be -1, 0, or 1"
    color = payload.get("color")
    if color is not None and not isinstance(color, str):
        return "'color' must be a string"
    notes = _sanitize_note(payload.get("notes", ""))
    return {
        "name":         name,
        "sequence":     seq,
        "feature_type": _sanitize_feat_type(ftype),
        "strand":       int(strand),
        "color":        _sanitize_label(color),
        "notes":        notes or "",
    }


@_agent_endpoint("list-feature-library")
def _h_list_feature_library(app, payload):
    """Return the feature library. Each item: ``{name, feature_type,
    strand, sequence_length, color}``. Sequences themselves can be
    long; agents needing the raw bases call `get-feature-library`.

    Named `list-feature-library` (not `list-features`) to avoid
    collision with the record-level feature endpoints (`add-feature`,
    `delete-feature`, `update-feature`, `get-feature`) which all
    operate on the LOADED record, not the persistent snippet library."""
    rows = []
    for e in _load_features():
        rows.append({
            "name":            e.get("name", ""),
            "feature_type":    e.get("feature_type", ""),
            "strand":          int(e.get("strand", 1)),
            "sequence_length": len(e.get("sequence", "") or ""),
            "color":           e.get("color", ""),
        })
    return {"ok": True, "features": rows, "count": len(rows)}


@_agent_endpoint("get-feature-library")
def _h_get_feature_library(app, payload):
    """Fetch a feature-library entry by (name, feature_type). Body:
    `{name, feature_type?}`. If feature_type is omitted, matches the
    first entry by name. Returns the full entry including `sequence`
    + `notes` so the agent can use it for downstream substring/BLAST
    matches.

    Named `get-feature-library` to disambiguate from `get-feature`
    which fetches a feature on the loaded record by idx."""
    name = _sanitize_label(payload.get("name"), max_len=200)
    if not name:
        return ({"error": "missing or non-string 'name'"}, 400)
    ftype = payload.get("feature_type")
    if ftype is not None and not isinstance(ftype, str):
        return ({"error": "'feature_type' must be a string"}, 400)
    for e in _load_features():
        if e.get("name") != name:
            continue
        if ftype and (e.get("feature_type") or "") != ftype:
            continue
        return {"ok": True, "feature": e}
    return ({"error": (
        f"no feature named {name!r}"
        + (f" of type {ftype!r}" if ftype else "")
    )}, 404)


@_agent_endpoint("create-feature-library", write=True)
def _h_create_feature_library(app, payload):
    """Add a feature snippet to the persistent feature library. Body:
    see `_agent_feature_dict` for the schema. Returns 409 if a
    feature with the same (name, feature_type) already exists — pick
    a different name or use `update-feature-library`.

    Named `create-feature-library` to disambiguate from `add-feature`
    which adds a feature to the loaded record."""
    f = _agent_feature_dict(payload)
    if isinstance(f, str):
        return ({"error": f}, 400)
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        entries = _load_features()
        for e in entries:
            if (e.get("name") == f["name"]
                    and (e.get("feature_type") or "") == f["feature_type"]):
                return ({"error": (
                    f"feature {f['name']!r} (type {f['feature_type']!r}) "
                    f"already exists; use update-feature or pick a "
                    f"different name."
                )}, 409)
        entries.append(f)
        if (err := _agent_save_or_500(
                lambda: _save_features(entries),
                "features")) is not None:
            return err
    return {"ok": True, "name": f["name"],
            "feature_type": f["feature_type"]}


@_agent_endpoint("list-feature-presets")
def _h_list_feature_presets(app, payload):
    """List the built-in feature presets — the shipped, GenBank-verified
    catalogue of common plasmid elements (markers, origins, promoters,
    terminators, polyA signals, reporters, tags, recombination sites).

    Body (all optional): ``{category, search, include_sequence}``.
    ``category`` must be one of the names in `categories`; ``search``
    matches name / type / category / description / aliases. Sequences are
    omitted unless ``include_sequence`` is true, because the full catalogue
    is ~35 kb of bases — fetch one with `get-feature-preset` instead.

    Presets are READ-ONLY and are NOT part of the user's feature library;
    `list-feature-library` still returns only the user's own entries. Copy
    one across with `import-feature-preset`."""
    category = payload.get("category")
    if category is not None and not isinstance(category, str):
        return ({"error": "'category' must be a string"}, 400)
    known = _preset_categories()
    if category and category not in known:
        return ({"error": (
            f"unknown category {category!r}; known: {list(known)}"
        )}, 400)
    search = payload.get("search")
    if search is not None and not isinstance(search, str):
        return ({"error": "'search' must be a string"}, 400)
    want_seq = payload.get("include_sequence", False)
    if not isinstance(want_seq, bool):
        # Strict, because it changes the response SHAPE. The sibling
        # `annotate-from-presets` validates its boolean the same way; two
        # endpoints in one family disagreeing about what a flag accepts is
        # exactly what an agent trips over.
        return ({"error": "'include_sequence' must be a boolean"}, 400)
    rows = []
    for p in _preset_features():
        if category and p.get("category") != category:
            continue
        if search and not _preset_matches(p, search):
            continue
        row = {
            "name":            p.get("name", ""),
            "feature_type":    p.get("feature_type", ""),
            "category":        p.get("category", ""),
            "strand":          int(p.get("strand", 1)),
            "color":           p.get("color", ""),
            "sequence_length": len(p.get("sequence", "") or ""),
            "description":     p.get("description", ""),
            "source":          p.get("source", ""),
            "aliases":         list(p.get("aliases") or []),
        }
        if want_seq:
            row["sequence"] = p.get("sequence", "")
        rows.append(row)
    return {"ok": True, "presets": rows, "count": len(rows),
            "categories": list(known)}


@_agent_endpoint("get-feature-preset")
def _h_get_feature_preset(app, payload):
    """Fetch one built-in preset by name, including its bases. Body:
    ``{name, feature_type?}``. Name match is exact but case-insensitive.
    404 when there is no such preset."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing or non-string 'name'"}, 400)
    ftype = payload.get("feature_type")
    if ftype is not None and not isinstance(ftype, str):
        return ({"error": "'feature_type' must be a string"}, 400)
    preset = _find_preset(name, ftype)
    if preset is None:
        return ({"error": f"no preset named {name!r}"}, 404)
    return {"ok": True, "preset": preset}


@_agent_endpoint("import-feature-preset", write=True)
def _h_import_feature_preset(app, payload):
    """Copy one or more built-in presets into the user's feature library.

    Body: ``{name}`` OR ``{names: [...]}`` — passing both is a 400 rather
    than a silent choice. Repeats within ``names`` collapse to one import, so
    ``count`` always equals the number of entries the library gained.
    Each preset becomes an ordinary
    library entry (provenance folded into ``qualifiers.note``), replacing
    any existing entry with the same (name, feature_type) — the same
    latest-write-wins rule `create-feature-library` uses on a rename.
    404 if any requested name is unknown, and nothing is written in that
    case: the import is all-or-nothing rather than half-applied."""
    names = payload.get("names")
    if names is not None and payload.get("name") is not None:
        # Fail loud rather than pick one. A caller that sent both does not
        # agree with itself about what to import, and silently honouring
        # `names` while dropping `name` is the ignored-outcome-key footgun
        # the agent API has been bitten by before.
        return ({"error": "pass 'name' OR 'names', not both"}, 400)
    if names is None:
        one = payload.get("name")
        if not isinstance(one, str) or not one.strip():
            return ({"error": "missing 'name' or 'names'"}, 400)
        names = [one]
    if not isinstance(names, list) or not names:
        return ({"error": "'names' must be a non-empty list of strings"}, 400)
    if len(names) > 500:
        return ({"error": "'names' capped at 500 entries per call"}, 400)
    resolved: list[dict] = []
    missing: list[str] = []
    seen_keys: set = set()
    for n in names:
        if not isinstance(n, str) or not n.strip():
            return ({"error": "every entry in 'names' must be a string"}, 400)
        preset = _find_preset(n)
        if preset is None:
            missing.append(n)
            continue
        # De-duplicate on the resolved identity, not the spelling: "loxP"
        # twice, or "loxP" and "LOXP", are one import. Without this the
        # response reported `count: 2, replaced: 1` for a library that
        # gained exactly one entry, which is a lie an agent would act on.
        key = (preset.get("name"), preset.get("feature_type"))
        if key in seen_keys:
            continue
        seen_keys.add(key)
        resolved.append(preset)
    if missing:
        return ({"error": f"unknown preset(s): {missing}"}, 404)
    if not resolved:
        return ({"error": "no presets to import"}, 400)
    try:
        new_entries = [_preset_to_library_entry(p) for p in resolved]
    except (TypeError, ValueError) as exc:
        return ({"error": f"could not convert preset: {exc}"}, 400)
    # RMW under the cache lock so a concurrent feature write can't be lost.
    with _state._cache_lock:
        entries = _load_features()
        replaced = 0
        for entry in new_entries:
            key = (entry.get("name"), entry.get("feature_type"))
            before = len(entries)
            entries = [e for e in entries
                       if (e.get("name"), e.get("feature_type")) != key]
            replaced += before - len(entries)
            entries.append(entry)
        if (err := _agent_save_or_500(
                lambda: _save_features(entries),
                "features")) is not None:
            return err
    return {"ok": True, "imported": [e["name"] for e in new_entries],
            "count": len(new_entries), "replaced": replaced}


@_agent_endpoint("update-feature-library", write=True)
def _h_update_feature_library(app, payload):
    """Update a feature-library entry, looked up by (name,
    feature_type). Body: every field in `_agent_feature_dict`'s
    schema; missing fields preserve existing values.

    Named `update-feature-library` (not `update-feature`) because
    `update-feature` is already taken by the record-level handler
    that edits features on the loaded plasmid map."""
    name = _sanitize_label(payload.get("name"), max_len=200)
    if not name:
        return ({"error": "missing or non-string 'name'"}, 400)
    ftype = payload.get("feature_type")
    if ftype is not None and not isinstance(ftype, str):
        return ({"error": "'feature_type' must be a string"}, 400)
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        entries = _load_features()
        target_idx = None
        for i, e in enumerate(entries):
            if e.get("name") != name:
                continue
            if ftype and (e.get("feature_type") or "") != ftype:
                continue
            target_idx = i
            break
        if target_idx is None:
            return ({"error": (
                f"no feature named {name!r}"
                + (f" of type {ftype!r}" if ftype else "")
            )}, 404)
        # Merge: take existing entry as base, overlay user-provided
        # fields. `_agent_feature_dict` requires sequence — patch it
        # to use the existing sequence if the payload omits it.
        base = entries[target_idx]
        merged = {**base, **payload}
        # Default sequence/feature_type from existing if missing.
        if "sequence" not in payload:
            merged["sequence"] = base.get("sequence", "")
        if "feature_type" not in payload:
            merged["feature_type"] = base.get("feature_type", "misc_feature")
        f = _agent_feature_dict(merged)
        if isinstance(f, str):
            return ({"error": f}, 400)
        entries[target_idx] = f
        if (err := _agent_save_or_500(
                lambda: _save_features(entries),
                "features")) is not None:
            return err
    return {"ok": True, "name": f["name"],
            "feature_type": f["feature_type"]}


@_agent_endpoint("delete-feature-library", write=True)
def _h_delete_feature_library(app, payload):
    """Remove a feature-library entry by (name, feature_type). Body:
    `{name, feature_type?}`. Named `delete-feature-library` so it
    doesn't collide with `delete-feature` (which removes a feature
    from the loaded record)."""
    name = _sanitize_label(payload.get("name"), max_len=200)
    if not name:
        return ({"error": "missing or non-string 'name'"}, 400)
    ftype = payload.get("feature_type")
    if ftype is not None and not isinstance(ftype, str):
        return ({"error": "'feature_type' must be a string"}, 400)
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        entries = _load_features()
        new_list = [
            e for e in entries
            if not (e.get("name") == name
                    and (not ftype or (e.get("feature_type") or "") == ftype))
        ]
        if len(new_list) == len(entries):
            return ({"error": (
                f"no feature named {name!r}"
                + (f" of type {ftype!r}" if ftype else "")
            )}, 404)
        if (err := _agent_save_or_500(
                lambda: _save_features(new_list),
                "features")) is not None:
            return err
    return {"ok": True, "removed": len(entries) - len(new_list)}


@_agent_endpoint("classify-part")
def _h_classify_part(app, payload):
    """Classify a candidate part by digest-overhang matching, without
    persisting. Body: ``{sequence, circular?: bool = true,
    features?: list}`` — same interface `_classify_part_from_plasmid`
    uses internally. Returns the (grammar, position, type, level,
    oh5, oh3) tuple if a match is found, or `match=null` if not.
    """
    seq = payload.get("sequence")
    if not isinstance(seq, str) or not seq.strip():
        return ({"error": "missing or non-string 'sequence'"}, 400)
    seq_clean = "".join(ch for ch in seq.upper()
                          if ch in "ACGTRYWSMKBDHVN")
    if not seq_clean:
        return ({"error": "no IUPAC bases in 'sequence'"}, 400)
    # Cap input size so a multi-MB paste can't pin the classifier on
    # multiple-grammar digest loops.
    if len(seq_clean) > 1_000_000:
        return ({"error": "sequence exceeds 1 Mbp cap for classifier"},
                413)
    circular = bool(payload.get("circular", True))
    raw_feats = payload.get("features") or []
    if not isinstance(raw_feats, list):
        return ({"error": "'features' must be a list"}, 400)
    # Sweep #26: defense-in-depth length cap on the features list.
    # Pre-fix the only upstream bound was the 1 MiB body cap — a
    # payload of ~50k micro-dicts would saturate the classifier's
    # inner loops without exceeding the body cap.
    if len(raw_feats) > 10_000:
        return ({"error":
                  "'features' too long (max 10,000 entries)"}, 413)
    feats = [f for f in raw_feats if isinstance(f, dict)]
    try:
        result = _classify_part_from_plasmid(
            seq_clean, circular=circular, features=feats,
        )
    except Exception as exc:
        _log.exception("agent classify-part: classifier failed")
        return ({"error": f"classification failed: {_scrub_path(str(exc))}"},
                500)
    return {"ok": True, "match": result}


@_agent_endpoint("add-codon-table", write=True)
def _h_add_codon_table(app, payload):
    """Add or replace a codon-usage table. Body either fetches from
    Kazusa or stamps a raw dict directly:

      * Fetch:   ``{taxid: str, name?: str, source: "kazusa"}``
      * Genome:  ``{source: "genome", accession|taxid: str,
                   mode?: "heg"|"genome", name?: str}``
      * File:    ``{source: "file", path: str, mode?: "heg"|"genome",
                   name?: str}`` — build from a LOCAL CDS FASTA (offline)
      * Raw:     ``{name: str, taxid?: str, raw: {<codon>: count, ...}}``

    The Kazusa fetch is size-capped at `_KAZUSA_MAX_RESPONSE_BYTES`
    and timeout-bounded; the genome build downloads the CDS from the
    NCBI Datasets API (size-capped at `_NCBI_CDS_ZIP_MAX_BYTES`,
    default ``mode="heg"`` → ribosomal-protein bias); the raw path
    caps at 64 codons, validates each value is a non-negative int, and
    derives each codon's amino acid from the standard genetic code so
    the stored table is the canonical ``{codon: (aa, count)}`` shape the
    optimizer consumes (snag #17).
    """
    source = payload.get("source", "user")
    if not isinstance(source, str):
        return ({"error": "'source' must be string"}, 400)
    if source == "kazusa":
        taxid_raw = payload.get("taxid")
        if not isinstance(taxid_raw, str) or not taxid_raw.strip():
            return ({"error": "'taxid' required for kazusa source"}, 400)
        taxid = _sanitize_label(taxid_raw, max_len=32)
        if not taxid or not taxid.isdigit():
            return ({"error": "'taxid' must be a digit string"}, 400)
        name_in = _sanitize_label(payload.get("name"), max_len=200)
        try:
            raw, msg = _codon_fetch_kazusa(taxid)
        except Exception as exc:
            _log.exception("agent add-codon-table: Kazusa fetch failed")
            return ({"error": f"Kazusa fetch failed: {_scrub_path(str(exc))}"}, 502)
        if raw is None:
            return ({"error": msg or "Kazusa returned no data"}, 502)
        display = name_in or f"Species (taxid {taxid})"
        try:
            entry = _codon_tables_add(display, taxid, raw, source="kazusa")
        except (OSError, RuntimeError) as exc:
            _log.exception("agent add-codon-table: save failed")
            _notify_save_failure(_state._LIVE_APP_REF.get(),
                                  "Codon tables", exc)
            return ({"error": f"save failed: {_scrub_path(str(exc))}"}, 500)
        return {"ok": True, "entry": {
            "name":   entry["name"],
            "taxid":  entry["taxid"],
            "source": entry["source"],
        }}
    if source == "genome":
        # Build from an NCBI genome's CDS — accession (GCF_…/GCA_…) or taxid.
        query_raw = payload.get("accession") or payload.get("taxid")
        if not isinstance(query_raw, str) or not query_raw.strip():
            return ({"error": "'accession' or 'taxid' required for genome "
                              "source"}, 400)
        query = _sanitize_accession(query_raw) or ""
        if not query:
            return ({"error": "invalid 'accession' or 'taxid'"}, 400)
        mode = payload.get("mode", "heg")
        if mode not in ("heg", "genome"):
            return ({"error": "'mode' must be 'heg' or 'genome'"}, 400)
        name_in = _sanitize_label(payload.get("name"), max_len=200)
        try:
            raw, msg, meta = _genome_build_codon_table(query, mode)
        except Exception as exc:
            _log.exception("agent add-codon-table: genome build failed")
            return ({"error": f"genome build failed: {_scrub_path(str(exc))}"}, 502)
        if raw is None or meta is None:
            return ({"error": msg or "genome build returned no data"}, 502)
        display = (name_in or meta.get("organism")
                   or meta.get("accession") or "Genome codon table")
        try:
            entry = _codon_tables_add(display, meta.get("taxid", ""), raw,
                                      source="genome")
        except (OSError, RuntimeError) as exc:
            _log.exception("agent add-codon-table: save failed")
            _notify_save_failure(_state._LIVE_APP_REF.get(),
                                  "Codon tables", exc)
            return ({"error": f"save failed: {_scrub_path(str(exc))}"}, 500)
        return {"ok": True, "message": msg, "entry": {
            "name":   entry["name"],
            "taxid":  entry["taxid"],
            "source": entry["source"],
        }}
    if source == "file":
        # Build from a LOCAL CDS FASTA on disk — offline, no NCBI.
        path_raw = payload.get("path")
        if not isinstance(path_raw, str) or not path_raw.strip():
            return ({"error": "'path' required for file source"}, 400)
        mode = payload.get("mode", "heg")
        if mode not in ("heg", "genome"):
            return ({"error": "'mode' must be 'heg' or 'genome'"}, 400)
        p = _sanitize_path(path_raw)
        if p is None:
            return ({"error": "invalid 'path'"}, 400)
        # Agent read-endpoint hardening: refuse an ancestor symlink redirect.
        ancestor_err = _check_agent_read_path_ancestors(p)
        if ancestor_err:
            return ({"error": ancestor_err}, 400)
        name_in = _sanitize_label(payload.get("name"), max_len=200)
        try:
            raw, msg, meta = _file_build_codon_table(p, mode)
        except Exception as exc:
            _log.exception("agent add-codon-table: file build failed")
            return ({"error": f"file build failed: {_scrub_path(str(exc))}"}, 502)
        if raw is None or meta is None:
            return ({"error": msg or "file build returned no data"}, 400)
        display = name_in or meta.get("organism") or "Codon table from file"
        try:
            entry = _codon_tables_add(display, meta.get("taxid", ""), raw,
                                      source="file")
        except (OSError, RuntimeError) as exc:
            _log.exception("agent add-codon-table: save failed")
            _notify_save_failure(_state._LIVE_APP_REF.get(),
                                  "Codon tables", exc)
            return ({"error": f"save failed: {_scrub_path(str(exc))}"}, 500)
        return {"ok": True, "message": msg, "entry": {
            "name":   entry["name"],
            "taxid":  entry["taxid"],
            "source": entry["source"],
        }}
    # Raw path.
    name_in = _sanitize_label(payload.get("name"), max_len=200)
    if not name_in:
        return ({"error": "missing 'name'"}, 400)
    raw = payload.get("raw")
    if not isinstance(raw, dict):
        return ({"error": "'raw' must be a dict of {codon: count}"}, 400)
    if len(raw) > 64:
        return ({"error": "'raw' has more than 64 codons"}, 400)
    valid = set("ACGTU")
    cleaned: dict = {}
    for codon, count in raw.items():
        if not isinstance(codon, str) or len(codon) != 3:
            return ({"error": f"bad codon {codon!r} — must be 3-char str"},
                    400)
        if set(codon.upper()) - valid:
            return ({"error": f"non-IUPAC codon {codon!r}"}, 400)
        n = _coerce_int(count, name=f"raw[{codon!r}]")
        if isinstance(n, str):
            return ({"error": n}, 400)
        if n < 0:
            return ({"error": f"negative count for codon {codon!r}"}, 400)
        # Snag #17: the registry shape is ``{codon: (aa, count)}`` (see
        # `_parse_codon_tsv`), NOT ``{codon: count}``. Storing a bare int
        # used to "succeed" here, then crash the MISSION-CRITICAL codon
        # optimizer later with "cannot unpack non-iterable int object"
        # when it did ``for codon, (aa, n) in raw.items()``. Derive the
        # amino acid from the standard genetic code so the table is valid
        # and immediately usable. `.get` guards the (already excluded by
        # the IUPAC check above) non-ACGT-after-fold case.
        folded = codon.upper().replace("U", "T")
        aa = _CODON_GENETIC_CODE.get(folded)
        if aa is None:
            return ({"error": f"{codon!r} is not a valid codon"}, 400)
        cleaned[folded] = (aa, n)
    taxid_raw = payload.get("taxid", "")
    taxid = (_sanitize_label(taxid_raw, max_len=32) if taxid_raw else "")
    try:
        entry = _codon_tables_add(name_in, taxid, cleaned, source="user")
    except (OSError, RuntimeError) as exc:
        _log.exception("agent add-codon-table: save failed")
        return ({"error": f"save failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, "entry": {
        "name":   entry["name"],
        "taxid":  entry["taxid"],
        "source": entry["source"],
    }}


@_agent_endpoint("delete-codon-table", write=True)
def _h_delete_codon_table(app, payload):
    """Remove a codon-usage table by `taxid` or `name`. Built-in
    tables (source='builtin') cannot be removed."""
    key_raw = payload.get("taxid") or payload.get("name")
    if not isinstance(key_raw, str) or not key_raw.strip():
        return ({"error": "missing 'taxid' or 'name'"}, 400)
    key = key_raw.strip().lower()
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        entries = _codon_tables_load()
        target = None
        for e in entries:
            if (str(e.get("taxid") or "").lower() == key
                    or str(e.get("name") or "").lower() == key):
                target = e
                break
        if target is None:
            return ({"error": f"no codon table for {key_raw!r}"}, 404)
        if (target.get("source") or "") == "builtin":
            return ({"error": "built-in codon tables cannot be removed"}, 409)
        kept = [e for e in entries
                if (e.get("taxid") or e.get("name")) !=
                   (target.get("taxid") or target.get("name"))]
        try:
            _codon_tables_save(kept)
        except (OSError, RuntimeError) as exc:
            _log.exception("agent delete-codon-table: save failed")
            _notify_save_failure(_state._LIVE_APP_REF.get(),
                                  "Codon tables", exc)
            return ({"error": f"save failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, "removed": {
        "name":   target.get("name", ""),
        "taxid":  target.get("taxid", ""),
    }}


def _history_node_to_dict(node) -> "dict | None":
    """Convert a `_CommercialSaaSHistoryNode` tree to a plain-dict JSON
    shape. Mirrors the viewer's per-node fields so an agent gets the
    same surface the UI sees without driving a TUI.

    Sweep #10 (2026-05-20): iterative DFS with explicit depth + node
    caps. Pre-fix this recursed through `node.parents`, blowing the
    Python recursion limit (~1000) on hostile `.dna` imports carrying
    deeply-nested `<HistoryTree>` blocks (32 MB cap can yield 3M+
    nodes). Sibling helpers `walk`, `_history_node_count`,
    `HistoryScreen.populate` were already iterative for this exact
    reason; this one was overlooked.
    """
    if node is None:
        return None

    def _shell(n) -> dict:
        # Every modelled field, so an agent reading `get-history` sees
        # exactly what the History viewer shows. The pre-2026-08-04 shape
        # carried 6 keys and dropped the date, the node identity, the PCR
        # primers, the run conditions, the end chemistry and the primer
        # BINDING record — i.e. most of what a `.dna` actually documents.
        return {
            "name":       getattr(n, "name", "") or "",
            "operation":  getattr(n, "operation", "") or "",
            "seq_len":    int(getattr(n, "seq_len", 0) or 0),
            "circular":   bool(getattr(n, "circular", False)),
            "regenerated_sites": list(
                getattr(n, "regenerated_sites", []) or []
            ),
            "input_summaries":   list(
                getattr(n, "input_summaries", []) or []
            ),
            "node_id":       int(getattr(n, "node_id", 0) or 0),
            "date":          getattr(n, "date", "") or "",
            "resurrectable": bool(getattr(n, "resurrectable", False)),
            "oligos":        list(getattr(n, "oligos", []) or []),
            "primers":       list(getattr(n, "primer_details", []) or []),
            "hybridization_params": dict(
                getattr(n, "hybridization_params", {}) or {}),
            "parameters":    list(getattr(n, "parameters", []) or []),
            "end_modifications": dict(
                getattr(n, "end_modifications", {}) or {}),
            "sticky_ends":   dict(getattr(n, "sticky_ends", {}) or {}),
            "custom_map_label": getattr(n, "custom_map_label", "") or "",
            "feature_snapshot_count": int(
                getattr(n, "feature_snapshot_count", 0) or 0),
            "source_collapsed": bool(getattr(n, "source_collapsed", False)),
            # Backstop, same as the viewer's "Other recorded fields": a
            # field this codebase doesn't model still reaches the caller.
            "extra_attributes": dict(
                getattr(n, "extra_attributes", {}) or {}),
            "other_children": [
                {"tag": t, "count": c}
                for t, c in (getattr(n, "other_children", []) or [])
            ],
            "parents":    [],
        }

    root_dict = _shell(node)
    # Stack frame: (source_node, target_dict_under_construction, depth).
    stack: list = [(node, root_dict, 0)]
    n_seen = 1
    truncated = False
    while stack:
        src, dst, depth = stack.pop()
        if depth >= _HISTORY_NODE_MAX_DEPTH:
            # Bail out the parents-walk at this branch — root_dict
            # still carries the upstream shape. Mark on dst so the
            # caller surfaces a "history truncated" hint.
            dst["_truncated"] = "depth_cap"
            truncated = True
            continue
        for p in (getattr(src, "parents", []) or []):
            if n_seen >= _HISTORY_NODE_MAX_NODES:
                dst["_truncated"] = "node_cap"
                truncated = True
                break
            child_dict = _shell(p)
            dst["parents"].append(child_dict)
            n_seen += 1
            stack.append((p, child_dict, depth + 1))
    if truncated:
        root_dict["_truncated_at_node"] = n_seen
    return root_dict


@_agent_endpoint("get-history")
def _h_get_history(app, payload):
    """Return the construction-history tree for a library entry as
    nested JSON. Body: ``{name}`` (or ``{id}``). Returns ``history=null``
    when the entry has no recorded history (NOT a 404)."""
    name = _sanitize_label(payload.get("name"), max_len=200)
    eid_raw = payload.get("id")
    eid = (_sanitize_label(eid_raw, max_len=200)
           if isinstance(eid_raw, str) else "")
    if not name and not eid:
        return ({"error": "missing 'name' or 'id'"}, 400)
    # Targeted finder clones only the matched entry instead of `_load_library()`
    # deep-cloning the ENTIRE active-collection cache just to read one entry's
    # history_xml — the per-call full-clone anti-pattern that sweep #11's
    # `_find_library_entry_by_*` helpers exist to kill (this read-only handler
    # was missed). `id` wins over `name`, matching the diff-plasmid precedence.
    entry = _find_library_entry_by_id(eid) if eid else None
    if entry is None and name:
        entry = _find_library_entry_by_name(name)
    if entry is None:
        return ({"error": f"no entry matching name={name!r} id={eid!r}"},
                404)
    xml = entry.get("history_xml") or ""
    if not xml:
        return {"ok": True, "history": None,
                "name": entry.get("name"), "id": entry.get("id")}
    try:
        root = _parse_commercialsaas_history(xml)
    except ValueError as exc:
        _log.warning("agent get-history: malformed history_xml: %s", exc)
        return ({"error": f"malformed history XML: {exc}"}, 422)
    # Check the ROOT node's claims against the entry's actual bases — the
    # one node we can identify with certainty (an ancestor is only a name
    # and a length in the record, and resolving that to a library entry can
    # land on the wrong molecule). STRICTLY READ-ONLY: this reports a
    # disagreement, it never edits the sequence or the history to resolve
    # one. A parse failure just means no check, never a failed request.
    warnings: "list[str]" = []
    try:
        seq = str(getattr(
            _gb_text_to_record(entry.get("gb_text") or ""), "seq", "") or "")
        # `_parse_commercialsaas_history` returns None for a `<HistoryTree>`
        # with no `<Node>` — nothing to check against.
        if seq and root is not None:
            warnings = _history_node_warnings(
                root, seq, enzymes=_state._all_enzymes_hook())
    except Exception:
        _log.debug("agent get-history: consistency check skipped for %r",
                    entry.get("name"), exc_info=True)
    return {
        "ok": True,
        "name": entry.get("name"),
        "id":   entry.get("id"),
        "history": _history_node_to_dict(root),
        # Claims the record makes that the sequence doesn't support.
        # Empty list = everything checks out.
        "warnings": warnings,
    }


@_agent_endpoint("check-primer")
def _h_check_primer(app, payload):
    """Check ONE primer against a template: melting temp, GC%, and every
    3'-anchored binding site (both strands, wrap-aware on a circular
    template). Body: ``{primer, template?, circular?: bool = true,
    min_identity?: float = 0, max_sites?: int = 50}``.

    `template` defaults to the LOADED plasmid (with its own topology). With
    no template and nothing loaded, the melting temperature and GC% are
    still returned — `binds` is then null and a `note` says binding was not
    checked — so "what is this primer's Tm?" never needs a sequence.

    The read-only analog of the GUI "Primer Check" tab for a single
    oligo: `simulate-pcr` needs a PAIR; this answers "does this one
    primer anneal, where, and how well." A site requires an EXACT match
    over the primer's 3' seed; identity is then scored over the FULL
    primer, so a primer carrying a 5' tail (cloning / mutagenesis)
    still reports its site with the tail counted as mismatches. Use it
    to confirm a designed primer binds exactly once (unique) before
    saving, or to locate where a mystery oligo lands.

    Read-only — never touches the loaded record or any saved file.
    Returns ``{ok, primer, length, tm, gc_pct, circular, n_sites,
    binds, best_identity, sites:[{orientation, strand, foot_start,
    length, ident_pct, mismatches, confidence}]}`` — `foot_start` is the
    0-based footprint start on the cleaned template; sites are
    best-first (highest identity)."""
    primer_raw = payload.get("primer")
    if not isinstance(primer_raw, str) or not primer_raw.strip():
        return ({"error": "missing or non-string 'primer'"}, 400)
    template_raw = payload.get("template")
    if template_raw is None:
        template_raw = payload.get("sequence")
    if template_raw is not None and not isinstance(template_raw, str):
        return ({"error": "'template' must be a string"}, 400)
    _iupac = "ACGTRYWSMKBDHVN"
    primer = "".join(ch for ch in primer_raw.upper() if ch in _iupac)
    if not primer:
        return ({"error": "no IUPAC bases in 'primer'"}, 400)
    if len(primer) > 1000:
        return ({"error": "'primer' too long (max 1000 bp)"}, 413)
    loaded_default = False
    if template_raw is None or not template_raw.strip():
        rec = getattr(app, "_current_record", None)
        loaded = str(getattr(rec, "seq", "") or "") if rec is not None else ""
        if not loaded:
            gc = sum(1 for c in primer if c in "GCS")
            return {
                "ok": True, "primer": primer, "length": len(primer),
                "tm": _primer_tm(primer),
                "gc_pct": round(100.0 * gc / len(primer), 1),
                "binds": None, "n_sites": 0, "sites": [],
                "note": "no template given and no plasmid loaded — melting "
                        "temperature and GC% only; binding was not checked",
            }
        template_raw = loaded
        loaded_default = True
        if payload.get("circular") is None:
            payload = {**payload, "circular": _record_is_circular(rec)}
    template = "".join(ch for ch in template_raw.upper() if ch in _iupac)
    if not template:
        return ({"error": "no IUPAC bases in 'template'"}, 400)
    if len(template) > _PAIRWISE_MAX_LEN:
        return ({"error":
                  f"'template' exceeds {_PAIRWISE_MAX_LEN:,} bp cap"}, 413)
    circular = bool(payload.get("circular", True))
    min_identity = payload.get("min_identity", 0.0)
    try:
        min_identity = float(min_identity)
    except (TypeError, ValueError):
        return ({"error": "'min_identity' must be a number"}, 400)
    if not (0.0 <= min_identity <= 100.0):
        return ({"error": "'min_identity' must be in [0, 100]"}, 400)
    max_sites = _coerce_int(payload.get("max_sites", 50), name="max_sites")
    if isinstance(max_sites, str):
        return ({"error": max_sites}, 400)
    max_sites = max(1, min(max_sites, 500))
    try:
        sites = _primer_binding_sites(
            primer, template, len(template), circular=circular,
            min_identity_pct=min_identity, max_sites=max_sites,
        )
    except ValueError as exc:
        return ({"error": f"alignment failed: {_scrub_path(str(exc))}"}, 400)
    out_sites = []
    for s in sites:
        glyph, _color = _primer_check_confidence(s.get("ident_pct"))
        out_sites.append({
            "orientation": "forward" if s.get("strand") == 1 else "reverse",
            "strand":      s.get("strand"),
            "foot_start":  s.get("foot_start"),
            "length":      s.get("length"),
            "ident_pct":   round(float(s.get("ident_pct") or 0.0), 1),
            "mismatches":  s.get("mismatches"),
            "confidence":  glyph,
        })
    # S = G|C, so count it toward GC for the ambiguous-base case.
    gc = sum(1 for c in primer if c in "GCS")
    gc_pct = round(100.0 * gc / len(primer), 1)
    best = max((s["ident_pct"] for s in out_sites), default=None)
    return {
        "ok":            True,
        "primer":        primer,
        "length":        len(primer),
        "tm":            _primer_tm(primer),
        "gc_pct":        gc_pct,
        "circular":      circular,
        "n_sites":       len(out_sites),
        "binds":         bool(out_sites),
        "best_identity": best,
        "sites":         out_sites,
        **({"template": "the loaded plasmid"} if loaded_default else {}),
    }


@_agent_endpoint("check-primer-duplicates")
def _h_check_primer_duplicates(app, payload):
    """Scan the primer library for duplicate sequences. Returns a
    list of groups where multiple entries share the same canonical
    primer sequence — useful for an agent to flag before saving a
    designed primer. Empty list when the library is clean.

    Read-only (no `write=True`) — running it via the agent doesn't
    mutate state. The on-launch `PrimerDuplicatesModal` (sweep #3)
    does the dedupe, but it's splash-time only; this endpoint lets
    an agent inspect anytime.
    """
    entries = _load_primers()
    by_seq: dict = {}
    for e in entries:
        if not isinstance(e, dict):
            continue
        seq = _normalize_primer_seq(e.get("sequence") or "")
        if not seq:
            continue
        by_seq.setdefault(seq, []).append({
            "name":     e.get("name", ""),
            "sequence": seq,
            "tm":       e.get("tm"),
            "primer_type": e.get("primer_type", ""),
            "source":   e.get("source", ""),
        })
    duplicates = [grp for grp in by_seq.values() if len(grp) > 1]
    return {
        "ok": True,
        "duplicates": duplicates,
        "n_groups": len(duplicates),
        "n_total_entries": len(entries),
    }


@_agent_endpoint("list-pre-update-snapshots")
def _h_list_pre_update_snapshots(app, payload):
    """List pre-update snapshots created by `splicecraft update`.
    Returns each snapshot's id, from_version, mtime, and counts —
    same data the `--restore-pre-update` CLI shows."""
    snaps = _list_pre_update_snapshots()
    rows = []
    for s in snaps:
        rows.append({
            "id":                  s.get("id", ""),
            "path":                str(s.get("path", "")),
            "mtime":               s.get("mtime"),
            "from_version":        s.get("from_version", "?"),
            "from_python_version": s.get("from_python_version", "?"),
            "from_platform":       s.get("from_platform", "?"),
            "schema_version":      s.get("schema_version", 1),
            "n_files":             s.get("n_files", 0),
            "n_dirs":              s.get("n_dirs", 0),
            "total_size":          s.get("total_size", 0),
        })
    return {"ok": True, "snapshots": rows, "count": len(rows)}


@_agent_endpoint("simulate-gibson")
def _h_simulate_gibson(app, payload):
    """Dry-run a Gibson assembly without saving. Body:
    ``{fragments: [{name, sequence, features?}, ...],
        min_overlap?: int = 15, circular?: bool = true}``.

    Returns the full simulator result dict (overlaps, errors,
    warnings, product features, product sequence). Read-only — pair
    with ``gibson-assemble`` to commit.
    """
    fragments = payload.get("fragments")
    if not isinstance(fragments, list):
        return ({"error": "'fragments' must be a list"}, 400)
    if len(fragments) > 64:
        return ({"error": "too many fragments (max 64)"}, 400)
    cleaned: list[dict] = []
    for i, f in enumerate(fragments):
        if not isinstance(f, dict):
            return ({"error": f"fragment[{i}] is not a dict"}, 400)
        seq = f.get("sequence")
        if not isinstance(seq, str):
            return ({"error": f"fragment[{i}].sequence missing or "
                      "non-string"}, 400)
        if len(seq) > 1_000_000:
            return ({"error": f"fragment[{i}] exceeds 1 Mbp"}, 413)
        name = _sanitize_label(f.get("name"), max_len=80) or f"F{i+1}"
        feats_raw = f.get("features") or []
        if not isinstance(feats_raw, list):
            return ({"error": f"fragment[{i}].features must be a list"},
                    400)
        cleaned.append({
            "name":     name,
            "sequence": seq,
            "features": [x for x in feats_raw if isinstance(x, dict)],
        })
    min_overlap = _coerce_int(
        payload.get("min_overlap", _GIBSON_MIN_OVERLAP_BP),
        name="min_overlap",
    )
    if isinstance(min_overlap, str):
        return ({"error": min_overlap}, 400)
    # 2026-05-27 (audit-5 GB M3): bump min_overlap floor from 1 to
    # 10 bp at the agent endpoint. Pre-fix `min_overlap=1` was
    # accepted — a chain of four fragments all ending in `A` and
    # starting with `A` would "assemble" into junk with a 1-bp
    # overlap. Real Gibson assemblies use 20-40 bp homology; 10 bp
    # is already permissive while refusing the obvious abuse.
    _GIBSON_AGENT_MIN_FLOOR = 10
    if not (_GIBSON_AGENT_MIN_FLOOR <= min_overlap <= _GIBSON_MAX_OVERLAP_BP):
        return ({"error": f"min_overlap out of range "
                  f"[{_GIBSON_AGENT_MIN_FLOOR}, "
                  f"{_GIBSON_MAX_OVERLAP_BP}] (1-9 bp overlaps would "
                  f"accept biologically meaningless assemblies)"}, 400)
    circular = bool(payload.get("circular", True))
    try:
        result = _simulate_gibson_assembly(
            cleaned, min_overlap=min_overlap, circular=circular,
        )
    except Exception as exc:
        _log.exception("agent simulate-gibson: simulator failed")
        return ({"error": f"simulator failed: {_scrub_path(str(exc))}"}, 500)
    # See `simulate-golden-gate`: `ok` is the assembly's verdict, not the
    # simulator's exit status. A Gibson dry run that found no usable overlap
    # reports `ok: false` with the reason in `result`.
    return {"ok": bool(result.get("success")), "result": result}


_SC_STRAND_QUAL_CARRY = "SpliceCraft_strand"


def _agent_carry_feature_dicts(rec) -> "list[dict]":
    """Convert a SeqRecord's annotations into the simple
    ``{start, end, strand, type, label, color}`` dicts the traditional-cloning
    engine carries through a digest→ligate (the same shape
    `_gibson_record_from_result` reads back). Skips the ``source`` feature;
    wrap features are emitted with ``end < start`` so the engine's
    `_split_features_at_cuts` halves them at the origin. Feeds
    ``carry_annotations`` on `traditional-clone`."""
    out: "list[dict]" = []
    total = _seq_len(rec)
    _circular = _record_is_circular(rec)
    for f in getattr(rec, "features", None) or []:
        if getattr(f, "type", "") == "source":
            continue
        bounds = _feat_bounds(f, total, circular=_circular)
        if bounds is None:
            continue
        s, e, strand = bounds
        quals = getattr(f, "qualifiers", {}) or {}
        label = (quals.get("label") or quals.get("product") or [f.type])[0]
        # `int(strand or 1)` turned every ARROWLESS parent feature into a
        # forward-pointing one on the clone: 0 is a real strand value, not
        # a falsy default. The double-stranded convention rides a
        # qualifier rather than the location, so read it back to the dict
        # model's `2` or it degrades to arrowless here (2026-09-12).
        sc_q = quals.get(_SC_STRAND_QUAL_CARRY)
        dict_strand = 0 if strand is None else int(strand)
        if (isinstance(sc_q, list) and sc_q
                and str(sc_q[0]).strip().lower() == "double"):
            dict_strand = 2
        d = {"start": int(s), "end": int(e), "strand": dict_strand,
             "type": f.type, "label": str(label)[:200]}
        color = (quals.get("ApEinfo_fwdcolor") or [""])[0] or ""
        if color:
            d["color"] = str(color)
        note = (quals.get("note") or [None])[0]
        if note:
            d["note"] = str(note)
        out.append(d)
    return out


def _agent_carry_clean_bases(seq: str) -> str:
    """Uppercase + strip to the IUPAC alphabet (matching `_sanitize_bases`'s
    cleaning) so a passed ``vector_seq`` and a saved entry's stored sequence
    compare on equal footing for the carry_annotations exact-match gate."""
    return "".join(c for c in (seq or "").upper() if c in "ACGTRYWSMKBDHVN")


def _agent_resolve_carry_entry(name, coll_in, which):
    """Resolve ``name`` to ONE library entry for a ``carry_annotations``
    source, searching EVERY collection (active first) the way `load-entry`
    and `transfer-annotations` do. Returns ``(entry, None)`` or
    ``(None, (err, status))``.

    Was active-collection-only, via the `_find_library_entry_by_*`
    mirror finders. A vector in one collection and an insert in another
    therefore could not be named in the same `traditional-clone` call, so
    the caller had to create a staging collection, copy both parents in,
    clone, move the product out and delete the staging collection (agent
    field report 2026-09-12). `transfer-annotations` already resolved
    across collections and its docstring argued the case."""
    coll: "str | None" = None
    if coll_in not in (None, ""):
        coll = _normalize_collection_name(coll_in)
        if coll is None:
            return None, ({"error": f"invalid {which}_collection"}, 400)
        if not any(isinstance(c, dict) and c.get("name") == coll
                    for c in _iter_collections_readonly()):
            return None, ({"error":
                            f"no collection named {coll!r}"}, 404)
    hits = _agent_scan_library_for_key(name, coll)
    if not hits and coll is None:
        # Live active-mirror fallback: a fresh/minimal state can hold the
        # entry in the mirror without a materialised collection for the
        # scan to walk (same fallback `transfer-annotations` keeps).
        fb = (_find_library_entry_by_id(name)
              or _find_library_entry_by_name(name))
        if fb is not None:
            hits = [(_get_active_collection_name() or "", fb)]
    if not hits:
        where = (f" in collection {coll!r}" if coll else " in any collection")
        return None, ({"error": f"no library entry named {name!r}{where} to "
                                f"carry {which} annotations from"}, 404)
    if len(hits) > 1:
        active = _get_active_collection_name()
        active_hits = [h for h in hits if h[0] == active]
        if len(active_hits) == 1:
            hits = active_hits
        else:
            holders = sorted({h[0] for h in hits})
            return None, ({"error":
                            f"{name!r} matches entries in multiple "
                            f"collections ({', '.join(holders)}) — pass "
                            f"{{\"{which}_collection\": \"…\"}} to "
                            f"disambiguate",
                            "collections": holders}, 409)
    return hits[0][1], None


def _agent_carry_source_features(payload, vec_seq, ins_seq):
    """Source feature dicts for ``carry_annotations`` on `traditional-clone`
    from the ``vector_name`` / ``insert_name`` library entries the caller
    already passes for lineage. Returns
    ``(vec_feats, ins_feats, warnings, carried, err)`` where ``carried`` is
    ``{"vector": bool, "insert": bool, "any": bool}`` and ``err`` is
    ``(dict, status)`` or None. A no-op unless ``carry_annotations`` is set.

    Names resolve across EVERY collection (pass ``vector_collection`` /
    ``insert_collection`` to pin one; 409 on a genuine ambiguity), so a
    backbone and an insert filed apart need no staging copy.

    Fail-loud where it matters (agent-API feedback — never a silent ok:true):
      * ``carry_annotations:true`` with NEITHER name → 400 (nothing to source
        from).
      * A named entry that doesn't resolve → 404 (the caller named a
        non-existent entry).
    Soft-skip (reported in ``warnings`` AND in the per-side ``carried``
    flags, never corrupts coordinates):
      * A named entry whose stored sequence doesn't EXACTLY match the passed
        seq — carrying rotated/edited coordinates would mis-place every feature
        ([INV-127] catastrophic-class), so that side is skipped with a warning.
        Rotation-matching is intentionally NOT attempted.

    Pass ``carry_strict: true`` to turn that skip into a 422 instead. The
    soft skip is right by default (half an annotation set beats none), but
    it made the natural failure a call that answers ``ok: true`` with a
    product missing half its features, signalled only inside a nested
    ``warnings`` list — so a caller checking ``ok`` sails past it. The
    per-side ``carried`` flags are the cheap check; ``carry_strict`` is the
    loud one."""
    if not bool(payload.get("carry_annotations")):
        return [], [], [], {"vector": False, "insert": False, "any": False}, None
    vname = _sanitize_label(payload.get("vector_name"), max_len=200)
    iname = _sanitize_label(payload.get("insert_name"), max_len=200)
    if not vname and not iname:
        return None, None, None, None, (
            {"error": "carry_annotations:true needs 'vector_name' and/or "
                      "'insert_name' to source features from a saved library "
                      "entry"}, 400)
    strict = bool(payload.get("carry_strict"))
    warnings: "list[str]" = []
    carried = {"vector": False, "insert": False, "any": False}

    def _one(name, target_clean, which):
        """-> (feats, err). Records the per-side outcome in `carried`."""
        if not name:
            return [], None
        entry, err = _agent_resolve_carry_entry(
            name, payload.get(f"{which}_collection"), which)
        if err is not None:
            return None, err
        # `err is None` ⇒ the entry is populated; assert so pyright narrows
        # it off the (dict | None) union the tuple return widens to.
        assert entry is not None

        def _skip(msg):
            warnings.append(msg)
            if strict:
                return None, ({"error": msg, "carried": dict(carried)}, 422)
            return [], None

        gb = entry.get("gb_text") or ""
        if not gb:
            return _skip(f"{which} entry {name!r} has no stored sequence; "
                          "annotations not carried")
        try:
            rec = _gb_text_to_record(gb)
        except Exception:
            return _skip(f"{which} entry {name!r} could not be parsed; "
                          "annotations not carried")
        if _agent_carry_clean_bases(str(getattr(rec, "seq", "") or "")) \
                != target_clean:
            return _skip(
                f"{which}_seq doesn't match saved entry {name!r} "
                "(rotation-matching not supported) — annotations not carried; "
                "pass the entry's sequence verbatim")
        feats = _agent_carry_feature_dicts(rec)
        if feats:
            carried[which] = True
            carried["any"] = True
        return feats, None

    vec_feats, verr = _one(vname, _agent_carry_clean_bases(vec_seq), "vector")
    if verr is not None:
        return None, None, None, None, verr
    ins_feats, ierr = _one(iname, _agent_carry_clean_bases(ins_seq), "insert")
    if ierr is not None:
        return None, None, None, None, ierr
    return vec_feats, ins_feats, warnings, carried, None


def _agent_traditional_cloning_candidates(payload):
    """Shared digest + ligate core for `simulate-traditional-cloning` /
    `traditional-clone`. Returns ``(info, None)`` on success or
    ``(None, (error_dict, status))`` on a validation / biology failure.

    `info` = ``{insert, n_vector_fragments, products, incompatible}`` where
    ``products`` lists every COMPATIBLE ligation product (each with its full
    `features`, so the write path can build a record) and ``incompatible``
    lists the ``(vector_frag_idx, orientation)`` pairs whose sticky ends
    didn't match — surfaced for diagnostics.

    A LINEAR insert (default) is excised as the purified middle fragment
    (`_excise_pcr_insert` — the bench "cut + gel-purify the off-cuts" step,
    so primer pad / outside-the-cut bases never reach the product). A CIRCULAR
    insert (``insert_circular=true`` — cutting a cassette OUT of a plasmid) is
    digested with `_excise_fragment_pair` into its candidate fragments, EACH a
    candidate cassette (`insert_frag_idx`). The vector is likewise digested
    with `_excise_fragment_pair`, which REFUSES an ambiguous ≥3-cut circular
    digest. EVERY (vector fragment × insert fragment × orientation) is tried —
    nothing is picked by size (the never-assume-the-smaller-fragment rule); the
    caller disambiguates with ``vector_frag_idx`` / ``insert_frag_idx``."""
    vec_seq, e = _sanitize_bases(payload.get("vector_seq"))
    if e:
        return None, ({"error": f"'vector_seq' rejected: {e}"}, 400)
    ins_seq, e = _sanitize_bases(payload.get("insert_seq"))
    if e:
        return None, ({"error": f"'insert_seq' rejected: {e}"}, 400)
    vec_enz = payload.get("vector_enzymes")
    ins_enz = payload.get("insert_enzymes")
    if (not isinstance(vec_enz, list) or not (1 <= len(vec_enz) <= 6)
            or not all(isinstance(x, str) and x for x in vec_enz)):
        return None, ({"error":
                        "'vector_enzymes' must be a list of 1-6 enzyme "
                        "names"}, 400)
    if (not isinstance(ins_enz, list) or not (1 <= len(ins_enz) <= 2)
            or not all(isinstance(x, str) and x for x in ins_enz)):
        return None, ({"error":
                        "'insert_enzymes' must be a list of 1 or 2 enzyme "
                        "names"}, 400)
    # Resolve names the same way `list-restriction-sites` does — exact first,
    # then case-insensitively — so `bsai` and `BsaI` behave identically here
    # and the error suggests the real spelling instead of just refusing. The
    # RESOLVED names replace the caller's, because the digest engine looks
    # enzymes up by exact catalog key.
    resolved_by_field: "dict[str, list[str]]" = {}
    for grp, field in ((vec_enz, "vector_enzymes"), (ins_enz, "insert_enzymes")):
        resolved, unknown = _resolve_enzyme_names(grp)
        if unknown:
            return None, _agent_unknown_enzyme_error(unknown, field=field)
        resolved_by_field[field] = resolved
    vec_enz = resolved_by_field["vector_enzymes"]
    ins_enz = resolved_by_field["insert_enzymes"]
    # carry_annotations (agent-API feedback): lift the vector/insert entries'
    # own features onto the ligated product so a raw cut/ligate clone isn't
    # feature-bare (no follow-up transfer-annotations pass needed). Sourced
    # from the `vector_name`/`insert_name` entries; fails loud when it can't be
    # honored, soft-skips (reported) on a seq mismatch. No-op (empty) unless
    # `carry_annotations` is set — the engine already splits + shifts these onto
    # the product coordinates.
    (carry_vec_feats, carry_ins_feats, carry_warnings,
     carried_sides, carry_err) = _agent_carry_source_features(
        payload, vec_seq, ins_seq)
    if carry_err is not None:
        return None, carry_err
    vector_circular = bool(payload.get("vector_circular", True))
    insert_circular = bool(payload.get("insert_circular", False))
    enz_l = ins_enz[0]
    enz_r = ins_enz[1] if len(ins_enz) > 1 else ins_enz[0]
    # 1) Build the candidate insert fragment(s).
    #   * Linear insert (default — a PCR product): the unambiguous purified
    #     MIDDLE fragment (`_excise_pcr_insert`), so primer pad / off-cut bases
    #     never reach the product. ONE candidate (insert_frag_idx 0).
    #   * Circular insert (`insert_circular=true` — excising a cassette OUT of a
    #     plasmid, e.g. an Ω multigene): a 2-cut digest leaves TWO fragments and
    #     EITHER could be the cassette, so enumerate both (insert_frag_idx 0/1)
    #     and let the caller pick — the same never-assume-the-smaller-fragment
    #     rule the vector side honours (a ≥3-cut circular insert is refused).
    if insert_circular:
        ins_frags, ierr = _excise_fragment_pair(
            ins_seq, list(ins_enz), circular=True,
            features=carry_ins_feats, source_label="insert")
        if ierr is not None:
            return None, ({"error": f"insert: {ierr['error']}",
                            "cuts": ierr.get("cuts", [])}, 422)
    else:
        single, ierr = _excise_pcr_insert(ins_seq, enz_l, enz_r,
                                            features=carry_ins_feats,
                                            source_label="insert")
        if ierr is not None:
            return None, ({"error": f"insert: {ierr}"}, 422)
        assert single is not None          # ierr is None ⇒ frag is populated
        ins_frags = [single]
    # 2) Digest the vector — refuses an ambiguous ≥3-cut circular digest.
    vec_frags, verr = _excise_fragment_pair(
        vec_seq, list(vec_enz), circular=vector_circular,
        features=carry_vec_feats, source_label="vector")
    if verr is not None:
        return None, ({"error": f"vector: {verr['error']}",
                        "cuts": verr.get("cuts", [])}, 422)
    # 3) Try every (vector fragment × insert fragment × orientation); keep the
    #    compatible ones, each tagged with BOTH fragment indices.
    products: list = []
    incompatible: list = []
    for vi, vfrag in enumerate(vec_frags):
        for ii, ifrag in enumerate(ins_frags):
            try:
                result = _simulate_traditional_cloning_multi([ifrag], vfrag)
            except Exception as exc:
                _log.exception("agent traditional-cloning: simulate failed")
                return None, ({"error":
                                f"simulator failed: {_scrub_path(str(exc))}"},
                              500)
            for orient in ("forward", "reverse"):
                prod = result.get(orient) or {}
                tag = {"vector_frag_idx": vi, "insert_frag_idx": ii,
                       "orientation": orient}
                if prod.get("compatible"):
                    top = prod.get("top_seq", "")
                    products.append({
                        **tag,
                        "length":     len(top),
                        "n_features": len(prod.get("features") or []),
                        "sequence":   top,
                        "features":   prod.get("features") or [],
                        # Off-target closures for THIS vector fragment
                        # (empty vector, double insert, …) — the same
                        # report the Constructor shows [INV-192]. Per
                        # product because it depends on which digest
                        # fragment plays the backbone.
                        "self_ligation": result.get("self_ligation"),
                    })
                else:
                    incompatible.append(tag)
    # Per-candidate insert summary (insert_frag_idx → length + cut ends).
    insert_fragments = [
        {"insert_frag_idx": ii,
         "length":      len(f.get("top_seq", "")),
         "left_enzyme":  (f.get("left")  or {}).get("enzyme", "") or enz_l,
         "right_enzyme": (f.get("right") or {}).get("enzyme", "") or enz_r}
        for ii, f in enumerate(ins_frags)]
    return ({
        # `insert` keeps the single-fragment summary for back-compat (the
        # linear/default case has exactly one candidate); `insert_fragments`
        # lists every candidate — the disambiguation set when the insert is
        # circular.
        "insert": {"length": len(ins_frags[0].get("top_seq", "")),
                   "left_enzyme": enz_l, "right_enzyme": enz_r},
        "insert_fragments":   insert_fragments,
        "n_vector_fragments": len(vec_frags),
        "n_insert_fragments": len(ins_frags),
        # Canonical catalog spellings of what was actually used — the write
        # path re-digests the vector for its backbone-marker check and must
        # not re-resolve the caller's raw names a second time.
        "vector_enzymes":     list(vec_enz),
        "insert_enzymes":     list(ins_enz),
        "products":           products,
        "incompatible":       incompatible,
        "carry_warnings":     carry_warnings,
        # Per-side outcome, not just an any-of bool: a caller that checks
        # `ok` alone must still be able to see that ONE parent's features
        # were skipped (agent field report 2026-09-12).
        "carried":             carried_sides,
        "annotations_carried": bool((carried_sides or {}).get("any")),
    }, None)


def _agent_vector_backbone_flags(payload, vec_enzymes, circular):
    """Which digested vector fragments carry a bacterial-backbone marker (an
    origin of replication or a resistance gene)?

    Returns ``(flags, why)`` — ``flags[i]`` is True when vector fragment `i`
    carries one and ``why[i]`` names the features that said so — or
    ``(None, {})`` when the vector's annotations aren't available to judge
    from. **Every failure path here is SOFT.** This feeds a convenience
    default on `traditional-clone`; if it can't tell, the caller just gets the
    409 it would have got anyway. It must never turn a workable clone into an
    error.

    Deliberately independent of ``carry_annotations``: knowing WHICH half of a
    digest is the backbone and putting the parent's features onto the product
    are different questions, and a caller shouldn't have to ask for the second
    to get the first. It reuses `carry_annotations`' exact-sequence gate for
    the same reason that gate exists — markers read off a rotated or edited
    entry would be attributed to the wrong fragment ([INV-127]).

    It does NOT look at fragment size. A stacked-TU insert can outgrow its
    carrier vector, which is exactly why the size heuristic is banned here."""
    name = _sanitize_label(payload.get("vector_name"), max_len=200)
    if not name:
        return None, {}
    entry = (_find_library_entry_by_id(name)
             or _find_library_entry_by_name(name))
    if entry is None or not (entry.get("gb_text") or ""):
        return None, {}
    try:
        rec = _gb_text_to_record(entry["gb_text"])
    except Exception:
        return None, {}
    vec_seq, err = _sanitize_bases(payload.get("vector_seq"))
    if err:
        return None, {}
    if (_agent_carry_clean_bases(str(getattr(rec, "seq", "") or ""))
            != _agent_carry_clean_bases(vec_seq)):
        return None, {}
    feats = _agent_carry_feature_dicts(rec)
    if not feats:
        return None, {}
    try:
        frags, ferr = _excise_fragment_pair(
            vec_seq, list(vec_enzymes), circular=circular,
            features=feats, source_label="vector")
    except Exception:                       # pragma: no cover - defensive
        return None, {}
    if ferr is not None or not frags:
        return None, {}
    flags: "list[bool]" = []
    why: "dict[int, str]" = {}
    for i, f in enumerate(frags):
        hit = _fragment_has_backbone_marker(f)
        flags.append(hit)
        if hit:
            why[i] = ", ".join(_fragment_backbone_marker_labels(f)[:4])
    return flags, why


@_agent_endpoint("simulate-traditional-cloning")
def _h_simulate_traditional_cloning(app, payload):
    """Dry-run a traditional restriction-cloning reaction: digest a vector
    and an insert with the chosen enzymes, then ligate. Body:
    ``{vector_seq, vector_enzymes, insert_seq, insert_enzymes,
    vector_circular?=true, insert_circular?=false}``.

    `insert_enzymes` is 1 or 2 enzyme names (the sites flanking the cassette).
    A LINEAR insert (default — e.g. a PCR product whose primers add the sites)
    is excised as the purified middle fragment. Set ``insert_circular=true`` to
    cut a cassette OUT of a PLASMID (e.g. an Ω multigene): the insert is
    digested like the vector, leaving TWO candidate fragments — each surfaced
    in ``insert_fragments`` with an ``insert_frag_idx`` — and NEITHER is picked
    by size. `vector_enzymes` digests the backbone — a circular vector (or
    insert) must cut EXACTLY twice (an ambiguous ≥3-cut digest is REFUSED with
    the cut list). Type IIS enzymes (BsaI / BbsI / BsmBI …) are refused for the
    insert — use the Golden-Braid path (`design-gb-part` /
    `scrub-plasmid method:golden_braid`).

    EVERY (vector fragment × insert fragment × orientation) is tried; the
    response lists every COMPATIBLE product as ``{vector_frag_idx,
    insert_frag_idx, orientation, length, n_features, sequence,
    self_ligation}`` (nothing is picked by size). ``incompatible`` lists the
    combos whose sticky ends didn't match.

    ``self_ligation`` is the off-target report for that product's backbone
    fragment: ``vector_self_closes`` (the empty vector re-circularises —
    single-enzyme, blunt, or a compatible-cohesive pair such as SalI/XhoI;
    with ``empty_vector_bp`` and the re-closed joint's ``vector_self_junction``
    so you know whether a diagnostic digest linearises empty colonies),
    ``double_insert`` (tandem / inverted two-copy inserts also close),
    ``orientation_ambiguous``, and a ``risks`` list of ``{kind, severity,
    message, advice}``; ``clean: true`` with a ``clean_note`` when nothing
    off-target closes. Read-only — pair with `traditional-clone` to save a
    chosen product."""
    info, err = _agent_traditional_cloning_candidates(payload)
    if err is not None:
        return err
    assert info is not None             # err is None ⇒ info is populated
    # Trim the verbose per-feature lists from the read response (the write
    # path recomputes them); keep the product sequence for inspection.
    products = [{k: v for k, v in p.items() if k != "features"}
                for p in info["products"]]
    return {"ok": True,
            "insert":              info["insert"],
            "insert_fragments":    info["insert_fragments"],
            "n_vector_fragments":  info["n_vector_fragments"],
            "n_insert_fragments":  info["n_insert_fragments"],
            "products":            products,
            "incompatible":        info["incompatible"],
            "annotations_carried": info["annotations_carried"],
            "carried":             info.get("carried", {}),
            "carry_warnings":      info["carry_warnings"],
            "ignored": _agent_ignored_keys(payload, {
                "vector_seq", "vector_enzymes", "insert_seq",
                "insert_enzymes", "vector_circular", "insert_circular",
                "carry_annotations", "carry_strict", "vector_name",
                "insert_name", "vector_collection", "insert_collection"})}


def _agent_gg_carry_features(payload, parts, vseq):
    """Source per-part + vector feature dicts for ``carry_annotations`` on
    `golden-gate-assemble`. Returns
    ``(part_features, vector_features, warnings, carried, err)``.

    Same contract as `_agent_carry_source_features` on `traditional-clone`,
    one layer wider: each part and the vector may carry a ``name`` (the
    spec shape ``{sequence, name?}`` already used for lineage) naming the
    library entry to lift features from, resolved across EVERY collection
    and gated on an EXACT sequence match so coordinates can never be
    mis-placed. ``collection`` inside a part spec (or ``vector_collection``)
    pins one.

    Before this there was no supported way to get an annotated one-pot
    product at all: the endpoint 400'd on ``carry_annotations`` and pointed
    at `transfer-annotations` (sequence-similarity — wrong instrument, and
    its 30 bp floor drops every promoter/RBS/overhang) or
    `assemble-into-entry-vector` (needs a grammar-registered entry vector,
    which a custom UPD acceptor is not). 26 saved one-pot products landed
    with zero features (agent field report 2026-09-12)."""
    if not bool(payload.get("carry_annotations")):
        return None, None, [], {"vector": False, "parts": [], "any": False}, None
    raw_parts = payload.get("parts") or []
    strict = bool(payload.get("carry_strict"))
    warnings: "list[str]" = []
    carried = {"vector": False, "parts": [False] * len(parts), "any": False}
    # Asking to carry annotations with nothing to carry them FROM is a
    # caller error, not an empty result — same fail-loud rule
    # `traditional-clone` applies to a missing vector_name/insert_name.
    _named = [r for r in ([payload.get("vector")] + list(raw_parts))
               if isinstance(r, dict) and str(r.get("name") or "").strip()]
    if not _named and not str(payload.get("vector_name") or "").strip():
        return None, None, None, None, (
            {"error": "carry_annotations:true needs a 'name' on the vector "
                      "and/or at least one part (pass each input as "
                      '{"sequence": …, "name": "<library entry>"}) to source '
                      "features from a saved library entry"}, 400)

    def _spec(raw):
        """-> (name, collection) from a part/vector spec."""
        if isinstance(raw, dict):
            nm = _sanitize_label(raw.get("name"), max_len=200)
            return nm, raw.get("collection")
        return "", None

    def _one(name, coll_in, clean_seq, which):
        if not name:
            warnings.append(
                f"{which} has no 'name' — pass {{\"sequence\": …, "
                f"\"name\": \"<library entry>\"}} to carry its features")
            return [], None
        entry, err = _agent_resolve_carry_entry(name, coll_in, which)
        if err is not None:
            return None, err
        # `err is None` ⇒ the entry is populated; assert so pyright narrows
        # it off the (dict | None) union the tuple return widens to.
        assert entry is not None

        def _skip(msg):
            warnings.append(msg)
            if strict:
                return None, ({"error": msg, "carried": dict(carried)}, 422)
            return [], None

        gb = entry.get("gb_text") or ""
        if not gb:
            return _skip(f"{which} entry {name!r} has no stored sequence; "
                          "annotations not carried")
        try:
            rec = _gb_text_to_record(gb)
        except Exception:
            return _skip(f"{which} entry {name!r} could not be parsed; "
                          "annotations not carried")
        if _agent_carry_clean_bases(str(getattr(rec, "seq", "") or "")) \
                != clean_seq:
            return _skip(
                f"{which} sequence doesn't match saved entry {name!r} "
                "(rotation-matching not supported) — annotations not "
                "carried; pass the entry's sequence verbatim")
        return _agent_carry_feature_dicts(rec), None

    vname, vcoll = _spec(payload.get("vector"))
    if not vname:
        vname = _sanitize_label(payload.get("vector_name"), max_len=200)
        vcoll = payload.get("vector_collection")
    vec_feats, err = _one(vname, vcoll, _agent_carry_clean_bases(vseq),
                           "vector")
    if err is not None:
        return None, None, None, None, err
    if vec_feats:
        carried["vector"] = True
        carried["any"] = True

    part_feats: "list[list[dict]]" = []
    for i, pseq in enumerate(parts):
        raw = raw_parts[i] if i < len(raw_parts) else None
        pname, pcoll = _spec(raw)
        feats, err = _one(pname, pcoll,
                           _agent_carry_clean_bases(pseq), f"part {i + 1}")
        if err is not None:
            return None, None, None, None, err
        part_feats.append(feats or [])
        if feats:
            carried["parts"][i] = True
            carried["any"] = True
    return part_feats, vec_feats, warnings, carried, None


def _agent_golden_gate_inputs(payload):
    """Parse + validate a Golden Gate payload ``{parts, vector, enzyme?}``.
    Returns ``(parts, vector_seq, enzyme, None)`` or
    ``(None, None, None, (error_dict, status))``. Each part / the vector may be
    a bare sequence string or ``{sequence}``."""
    parts = payload.get("parts")
    if not isinstance(parts, list) or not parts:
        return None, None, None, (
            {"error": "'parts' must be a non-empty list of sequences"}, 400)
    if len(parts) > 30:
        return None, None, None, ({"error": "too many parts (max 30)"}, 400)
    clean: list = []
    for i, p in enumerate(parts):
        raw = (p if isinstance(p, str)
               else p.get("sequence") if isinstance(p, dict) else None)
        seq, err = _sanitize_bases(raw)
        if err:
            return None, None, None, ({"error": f"part[{i}]: {err}"}, 400)
        clean.append(seq)
    vraw = payload.get("vector")
    vraw = (vraw if isinstance(vraw, str)
            else vraw.get("sequence") if isinstance(vraw, dict) else None)
    vseq, verr = _sanitize_bases(vraw)
    if verr:
        return None, None, None, ({"error": f"'vector': {verr}"}, 400)
    enzyme = _sanitize_label(payload.get("enzyme") or "BsaI", max_len=40)
    # Resolve case + commercial synonyms the same way every other enzyme
    # input does — the digest engines key on the CATALOG spelling, so
    # `Eco31I` (Thermo's BsaI) has to become `BsaI` here or the reaction is
    # simulated with an enzyme that cuts nothing.
    match = _enzyme_resolve_one(enzyme)
    if match is None:
        return None, None, None, ({"error": f"unknown enzyme {enzyme!r}"}, 400)
    return clean, vseq, match[0], None


@_agent_endpoint("simulate-golden-gate")
def _h_simulate_golden_gate(app, payload):
    """Dry-run a Golden Gate / MoClo (Type IIS) one-pot assembly. Body:
    ``{parts: [seq, …], vector, enzyme?="BsaI"}``.

    Digests each part + the vector with the Type IIS enzyme (BsaI / BsmBI /
    BbsI / SapI / Esp3I) and chains the released fragments by their 4 nt
    overhangs into a circular product — parts may be in ANY order, the
    overhangs set the order. Each part must release exactly one fragment
    (flank it with two inward enzyme sites).

    **The vector may be cut more than twice.** A destination plasmid carrying
    background enzyme sites is released as several pieces, and the real
    one-pot reaction puts ALL of them back around the parts; the design stays
    determinate as long as the overhangs are unique. Every released vector
    piece goes into the chain, and ``n_vector_fragments`` /
    ``n_vector_fragments_used`` say how many there were and how many the
    product kept.

    The ``result`` carries ``{ok, product_seq, length, n_parts, order,
    junctions:[{overhang}], n_residual_sites, n_vector_sites,
    n_vector_fragments, n_vector_fragments_used, vector_cut_bp, warnings,
    errors}`` — `warnings` flags a non-unique junction overhang, more than
    one circle that can form from the same fragments, or a residual enzyme
    site (it would be re-cut). When nothing closes, the error lists the
    vector's cut positions and EVERY fragment's two overhangs (plus
    ``vector_fragment_overhangs`` / ``part_overhangs``) so the dangling one
    is visible, instead of asserting the parts are at fault. Read-only —
    pair with `golden-gate-assemble` to save.

    ``carry_annotations: true`` previews the product's FEATURE map too,
    sourced from the ``name`` on each ``{sequence, name?}`` part / vector
    spec, exactly as `golden-gate-assemble` saves it (``result.features``,
    plus ``carried`` / ``carry_warnings``). The pair takes the same payload
    on purpose: a preview that silently dropped the key would disagree with
    the save that honoured it."""
    parts, vseq, enzyme, err = _agent_golden_gate_inputs(payload)
    if err is not None:
        return err
    assert parts is not None and vseq is not None and enzyme is not None
    (part_feats, vec_feats, carry_warnings,
     carried, carry_err) = _agent_gg_carry_features(payload, parts, vseq)
    if carry_err is not None:
        return carry_err
    try:
        result = _simulate_golden_gate(
            parts, vseq, enzyme=enzyme,
            part_features=part_feats, vector_features=vec_feats)
    except Exception as exc:
        _log.exception("agent simulate-golden-gate: assembler failed")
        return ({"error": f"assembler failed: {_scrub_path(str(exc))}"}, 500)
    # `ok` reports whether the ASSEMBLY worked, not whether the simulator
    # ran. Hardcoding it True meant a dry run that failed to close a circle
    # still answered `ok: true` at the top level while `result.ok` said
    # false — so the natural `if r["ok"]:` accepted a design that cannot be
    # built, and callers had to write `r.get("ok") or
    # r.get("result", {}).get("ok")` to get a straight answer (blue field
    # report #4). Still HTTP 200: the request was served, the reaction it
    # modelled is what failed. `errors`/`warnings` explain why.
    return {"ok": bool(result.get("ok")), "result": result,
            "carried": carried, "carry_warnings": carry_warnings,
            "annotations_carried": bool((carried or {}).get("any")),
            "ignored": _agent_ignored_keys(
                payload, {"parts", "vector", "enzyme", "carry_annotations",
                           "carry_strict", "vector_name",
                           "vector_collection"})}


_AGENT_MUT_RE = re.compile(r"^([A-Z])(\d{1,5})([A-Z\*])$")


@_agent_endpoint("design-mutagenesis")
def _h_design_mutagenesis(app, payload):
    """Design SOE-PCR primers for a single-site mutation. Body:
    ``{cds_dna, mutation, codon_taxid?}``. `mutation` is a string
    like ``"W140F"`` (WT-aa, 1-based position, mutant-aa).

    Returns ``{outer, inner}`` — the outer + inner primer pairs
    `_mut_design_outer` / `_mut_design_inner` produce. The agent is
    free to save the result via `update-primer` or by setting up a
    Constructor save flow.
    """
    cds_dna = payload.get("cds_dna") or payload.get("sequence")  # SC-G alias
    if not isinstance(cds_dna, str) or not cds_dna.strip():
        return ({"error": "missing or non-string 'cds_dna'"}, 400)
    cds_clean = "".join(ch for ch in cds_dna.upper()
                          if ch in "ACGTRYWSMKBDHVN")
    if not cds_clean:
        return ({"error": "no IUPAC bases in 'cds_dna'"}, 400)
    if len(cds_clean) > 30_000:
        return ({"error": "'cds_dna' exceeds 30 kbp cap"}, 413)
    if len(cds_clean) % 3 != 0:
        return ({"error": f"'cds_dna' length {len(cds_clean)} is not a "
                  "multiple of 3 (CDS must be whole codons)"}, 400)
    mut_raw = payload.get("mutation")
    if not isinstance(mut_raw, str):
        return ({"error": "missing or non-string 'mutation'"}, 400)
    m = _AGENT_MUT_RE.match(mut_raw.strip())
    if m is None:
        return ({"error": f"mutation {mut_raw!r} doesn't match the "
                  "pattern '<WT-aa><1-based-pos><mut-aa>' "
                  "(e.g. 'W140F')"}, 400)
    wt_aa, pos_str, mut_aa = m.group(1), m.group(2), m.group(3)
    pos = int(pos_str)
    if pos < 1:
        return ({"error": "mutation position must be >= 1"}, 400)
    if (pos - 1) * 3 + 3 > len(cds_clean):
        return ({"error": f"mutation position {pos} is past the end of "
                  f"the {len(cds_clean) // 3}-aa CDS"}, 400)
    codon_taxid = payload.get("codon_taxid")
    codon_raw = None
    if codon_taxid is not None:
        if not isinstance(codon_taxid, str):
            return ({"error": "'codon_taxid' must be string"}, 400)
        entry = _codon_tables_get(codon_taxid)
        if entry is None:
            return ({"error": f"unknown codon_taxid {codon_taxid!r}"},
                    404)
        codon_raw = entry.get("raw")
    try:
        outer = _mut_design_outer(cds_clean)
        inner = _mut_design_inner(cds_clean, pos, mut_aa, wt_aa,
                                    codon_table=codon_raw)
    except (ValueError, RuntimeError) as exc:
        return ({"error": f"mutagenesis design failed: {exc}"}, 422)
    except Exception as exc:
        _log.exception("agent design-mutagenesis: unexpected failure")
        return ({"error": f"unexpected failure: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, "mutation": mut_raw, "outer": outer, "inner": inner}


def _record_to_scrub_feats(record) -> "list[dict]":
    """Build scrub-compatible feature dicts (type/start/end/strand, plus
    codon_start/transl_table/_exons for CDS) from a SeqRecord for the
    headless agent path — the subset of `PlasmidMap._parse` that
    `_scrub_design` consumes, using the hardened `_feat_bounds` for
    wrap-aware coordinates. CDS features carry their reading frame so the
    scrub stays protein-preserving without a mounted UI."""
    total = _seq_len(record)
    circular = _record_is_circular(record)
    out: "list[dict]" = []
    for feat in getattr(record, "features", []) or []:
        if feat.type == "source":
            continue
        b = _feat_bounds(feat, total, circular=circular)
        if b is None:
            continue
        start, end, strand = b
        d: dict = {"type": feat.type, "start": start, "end": end,
                   "strand": strand or 1, "label": _feat_label(feat)}
        loc = getattr(feat, "location", None)
        try:
            from Bio.SeqFeature import CompoundLocation
            # Spliced CDS: pass exon parts so `_translate_cds` splices introns
            # out before checking that a scrub edit is synonymous.
            #
            # The test MUST mirror `PlasmidMap._parse` — a TWO-part join is the
            # commonest spliced CDS there is (one intron), and the old
            # `len(parts) > 2` gate skipped every one of them. The scrub then
            # read exon 2 through the intron in a shifted frame and cheerfully
            # accepted an edit that changes the real protein. Exclude only the
            # shapes where the parts are not exons: an origin WRAP (the two
            # halves are already encoded in start/end) and a CONTIGUOUS split
            # (abutting parts — no intron, nothing spliced out).
            if isinstance(loc, CompoundLocation) and len(loc.parts) >= 2:
                parts = sorted(loc.parts, key=lambda p: int(p.start))
                is_contiguous = all(
                    int(parts[i].end) == int(parts[i + 1].start)
                    for i in range(len(parts) - 1)
                )
                is_wrap = (
                    circular and total > 0 and len(parts) == 2
                    and int(parts[0].start) == 0
                    and int(parts[-1].end) == total
                    and int(parts[0].end) < int(parts[-1].start)
                )
                if not is_wrap and not is_contiguous:
                    d["_exons"] = [
                        (int(p.start), int(p.end)) for p in parts
                    ]
        except (ImportError, TypeError, ValueError):
            pass
        if feat.type.upper() == "CDS":
            try:
                cs = int((feat.qualifiers.get("codon_start", ["1"]) or ["1"])[0])
                if cs in (2, 3):
                    d["codon_start"] = cs
            except (TypeError, ValueError, IndexError):
                pass
            try:
                tt = int((feat.qualifiers.get("transl_table", ["1"])
                          or ["1"])[0])
                if tt and tt != 1:
                    d["transl_table"] = tt
            except (TypeError, ValueError, IndexError):
                pass
        out.append(d)
    return out


@_agent_endpoint("scrub-plasmid")
def _h_scrub_plasmid(app, payload):
    """Plan a clone-free restriction-site scrub. Body:
    ``{seq?, features?, enzymes?, overlap?, method?, circular?, codon_taxid?}``.

    With no ``seq``, scrubs the plasmid currently on the canvas (using its
    CDS features so the cure stays synonymous / protein-preserving). When
    ``seq`` is given, pass ``features`` (a list of ``{type,start,end,strand,
    codon_start?,transl_table?}`` dicts) to protect overlapping CDSes —
    WITHOUT them every site is treated as non-coding and a coding site could
    change a protein. ``enzymes`` is a list of names (default BsaI/Esp3I/BbsI)
    — matched case-insensitively, commercial synonyms accepted (``Eco31I`` is
    ``BsaI``), and an unrecognised name is a **400**, never a scrub that
    quietly removes nothing and reports `ok: true` with an empty
    ``sites_removed``; ``overlap`` is ``"improved"`` (default) or
    ``"classic"``; ``codon_taxid`` (a registered codon-usage table id) makes
    coding cures prefer that host's frequent synonymous codons.

    ``method`` picks the re-circularization route: ``"quikchange"`` (default —
    one whole-plasmid amplicon that self-circularises) or ``"golden_braid"``
    (split into BsaI-tailed PCR fragments that a Golden Gate reaction
    reassembles seamlessly; BsaI is force-cured as the assembly enzyme).

    QuikChange returns ``{ok, method, enzymes, cured_seq, edits,
    sites_removed, sites_skipped, n_rounds, rounds, warnings}``. Golden Braid
    returns ``{ok, method, enzyme, cured_seq, edits, sites_removed,
    sites_skipped, n_fragments, fragments, verified, warnings, errors}`` —
    each fragment carrying its BsaI-tailed primer pair + junction overhangs.
    Design-only: never mutates the canvas or any file.
    """
    seq = payload.get("seq")
    feats: "list[dict]" = []
    if seq is not None:
        if not isinstance(seq, str) or not seq.strip():
            return ({"error": "'seq' must be a non-empty string"}, 400)
        feats_in = payload.get("features", [])
        if feats_in is None:
            feats_in = []
        if not isinstance(feats_in, list):
            return ({"error": "'features' must be a list of feature dicts"}, 400)
        feats = feats_in
    else:
        rec = getattr(app, "_current_record", None)
        if rec is None:
            return ({"error": "no 'seq' provided and no plasmid is loaded"}, 400)
        seq = str(rec.seq)
        try:
            feats = _record_to_scrub_feats(rec)
        except Exception:
            _log.exception("scrub endpoint: feature extraction failed")
            feats = []
    if len(seq) > 1_000_000:
        return ({"error": "'seq' exceeds the 1 Mbp cap"}, 413)
    enzymes = payload.get("enzymes")
    if enzymes is not None:
        # `_scrub_resolve_sites` SKIPS a name the catalog doesn't know, so an
        # unrecognised (or merely mis-cased, or supplier-named) enzyme used to
        # scrub nothing and answer `ok: true, sites_removed: []` — which reads
        # as "your plasmid is already clean". Same vacuous pass
        # `list-restriction-sites` refuses, in the endpoint where being wrong
        # costs the most. [INV-189]
        if (shape_err := _agent_check_enzyme_list(enzymes)) is not None:
            return shape_err
        resolved_enz, unknown_enz = _resolve_enzyme_names(enzymes)
        if unknown_enz:
            return _agent_unknown_enzyme_error(
                unknown_enz, hint=" — call 'list-enzymes' for the catalog")
        if not resolved_enz:
            return ({"error":
                      "'enzymes' is empty — omit the key to scrub the "
                      "default set (BsaI / Esp3I / BbsI)"}, 400)
        enzymes = resolved_enz
    overlap = payload.get("overlap", "improved")
    if overlap not in ("improved", "classic"):
        return ({"error": "'overlap' must be 'improved' or 'classic'"}, 400)
    method = payload.get("method", "quikchange")
    if method not in ("quikchange", "golden_braid"):
        return ({"error": "'method' must be 'quikchange' or 'golden_braid'"}, 400)
    circular = payload.get("circular", True)
    if not isinstance(circular, bool):
        return ({"error": "'circular' must be a boolean"}, 400)
    codon_taxid = payload.get("codon_taxid")
    codon_raw = None
    if codon_taxid is not None:
        if not isinstance(codon_taxid, str):
            return ({"error": "'codon_taxid' must be a string"}, 400)
        entry = _codon_tables_get(codon_taxid)
        if entry is None:
            return ({"error": f"unknown codon_taxid {codon_taxid!r}"}, 404)
        codon_raw = entry.get("raw")
    rounds: list = []
    try:
        if method == "golden_braid":
            plan = _scrub_gb_design(seq, feats, enzymes, circular=circular,
                                    codon_raw=codon_raw)
        else:
            plan = _scrub_design(seq, feats, enzymes, circular=circular,
                                 codon_raw=codon_raw)
            if plan.get("ok"):
                for i, cl in enumerate(plan.get("clusters", []), 1):
                    rounds.append(_scrub_qc_primers(
                        plan["cured_seq"], cl["positions"],
                        circular=circular, overlap=overlap, round_no=i))
                if any(not r.get("error") for r in rounds):
                    qv_ok, _qv = _scrub_qc_verify(
                        plan.get("orig_seq", ""), plan["cured_seq"],
                        rounds, len(plan["cured_seq"]))
                    plan["verified"] = qv_ok
    except Exception as exc:
        _log.exception("agent scrub-plasmid: unexpected failure")
        return ({"error": f"unexpected failure: {_scrub_path(str(exc))}"}, 500)
    if method == "golden_braid":
        return {
            "ok": plan.get("ok", False),
            "method": "golden_braid",
            "enzyme": plan.get("enzyme", "BsaI"),
            "cured_seq": plan.get("cured_seq", ""),
            "edits": plan.get("edits", []),
            "sites_removed": plan.get("sites_removed", []),
            "sites_skipped": plan.get("sites_skipped", []),
            "n_fragments": plan.get("n_fragments", 0),
            "fragments": plan.get("fragments", []),
            "verified": plan.get("verified", False),
            # Longest homopolymer run in the CURED sequence. The cure scorer
            # avoids lengthening one, but where every silent option does, the
            # caller has to be able to see it rather than assert it by hand.
            "max_homopolymer_run": plan.get("max_homopolymer_run", 0),
            "warnings": plan.get("warnings", []),
            "errors": plan.get("errors", []),
        }
    return {
        "ok": plan.get("ok", False),
        "method": "quikchange",
        "verified": plan.get("verified", False),
        "enzymes": plan.get("enzymes", []),
        "cured_seq": plan.get("cured_seq", ""),
        "edits": plan.get("edits", []),
        "sites_removed": plan.get("sites_removed", []),
        "sites_skipped": plan.get("sites_skipped", []),
        "n_rounds": plan.get("n_rounds", 0),
        "rounds": rounds,
        "max_homopolymer_run": plan.get("max_homopolymer_run", 0),
        "warnings": plan.get("warnings", []),
    }


@_agent_endpoint("get-experiment")
def _h_get_experiment(app, payload):
    """Fetch one notebook entry (full body + metadata) by id.
    Body: ``{id: str}``."""
    eid = _sanitize_experiment_id(payload.get("id"))
    if eid is None:
        return ({"error": "missing or invalid 'id'"}, 400)
    for e in _load_experiments():
        if e.get("id") == eid:
            return {"experiment": _typed_clone(e)}
    return ({"error": f"no experiment with id {eid!r}"}, 404)


@_agent_endpoint("create-experiment", write=True)
def _h_create_experiment(app, payload):
    """Create a new notebook entry in the active project. Body:
    ``{title?: str = "Untitled entry",
        body_md?: str = "",
        tags?: list[str] = []}``.

    Routes through `_normalise_experiment_entry` so size caps + tag
    dedup + plasmid/action/gel xref extraction apply automatically.
    Returns the freshly-stamped entry id."""
    title = payload.get("title") or "Untitled entry"
    body  = payload.get("body_md") or ""
    tags  = payload.get("tags") or []
    if not isinstance(tags, list):
        return ({"error": "'tags' must be a list of strings"}, 400)
    proj_in = payload.get("project")
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        # Determine the target store: a SPECIFIC named project (SC-M-class —
        # `project` used to be silently ignored, filing into the active one),
        # or the active project's live store.
        named_proj = None
        projs = None
        if proj_in not in (None, ""):
            if not isinstance(proj_in, str):
                return ({"error": "'project' must be a string"}, 400)
            pname = proj_in.strip()
            if pname != (_get_active_project_name() or ""):
                projs = _load_experiment_projects()
                named_proj = next((p for p in projs
                                   if (p.get("name") or "") == pname), None)
                if named_proj is None:
                    return ({"error":
                              f"no experiment project named {pname!r}; create "
                              f"it with create-experiment-project"}, 404)
                entries = named_proj.setdefault("experiments", [])
            else:
                entries = _load_experiments()
        else:
            entries = _load_experiments()
        existing_ids: set[str] = {
            e.get("id") for e in entries if e.get("id")  # type: ignore[misc]
        }
        new_id = _new_experiment_id(existing_ids)
        entry = _normalise_experiment_entry({
            "id":      new_id,
            "title":   title if isinstance(title, str) else "",
            "body_md": body  if isinstance(body, str)  else "",
            "tags":    tags,
        }, fresh=True)
        entries.append(entry)
        if named_proj is not None:
            assert projs is not None     # `named_proj` set ⇒ projs was loaded
            err = _agent_save_or_500(
                lambda: _save_experiment_projects(projs),
                "experiment_projects")
        else:
            err = _agent_save_or_500(
                lambda: _save_experiments(entries), "Experiments")
        if err:
            return err
    _log_event("experiments.new", eid=new_id, via="agent")
    return {"ok": True, "id": new_id, "experiment": _typed_clone(entry),
            "project": (proj_in.strip()
                        if isinstance(proj_in, str) and proj_in.strip()
                        else _get_active_project_name() or "")}


@_agent_endpoint("update-experiment", write=True)
def _h_update_experiment(app, payload):
    """Update an existing notebook entry. Body:
    ``{id: str, title?: str, body_md?: str, tags?: list[str]}``.
    Fields not supplied keep their prior value. Re-normalises so
    body-byte cap / tag dedup / xref extraction stay consistent."""
    eid = _sanitize_experiment_id(payload.get("id"))
    if eid is None:
        return ({"error": "missing or invalid 'id'"}, 400)
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        entries = _load_experiments()
        for i, e in enumerate(entries):
            if e.get("id") != eid:
                continue
            merged = dict(e)
            if "title" in payload:
                t = payload["title"]
                merged["title"] = t if isinstance(t, str) else ""
            if "body_md" in payload:
                b = payload["body_md"]
                merged["body_md"] = b if isinstance(b, str) else ""
            if "tags" in payload:
                tg = payload["tags"]
                if not isinstance(tg, list):
                    return ({"error": "'tags' must be a list of strings"}, 400)
                merged["tags"] = tg
            entries[i] = _normalise_experiment_entry(merged, fresh=False)
            err = _agent_save_or_500(
                lambda: _save_experiments(entries), "Experiments",
            )
            if err:
                return err
            _log_event("experiments.save", eid=eid, via="agent")
            return {"ok": True, "experiment": _typed_clone(entries[i])}
    return ({"error": f"no experiment with id {eid!r}"}, 404)


@_agent_endpoint("create-experiment-project", write=True)
def _h_create_experiment_project(app, payload):
    """Create a new (empty) experiment project. Body:
    ``{name: str, description?: str = ""}``. Active-project pointer
    is NOT modified by this endpoint — call
    ``set-active-experiment-project`` afterwards if you want to
    work in the new project."""
    name = _sanitize_label(payload.get("name"), max_len=200)
    if not name:
        return ({"error": "missing 'name'"}, 400)
    desc = _sanitize_note(payload.get("description"))
    # Sweep #26: RMW under `_state._cache_lock`. Pre-fix the name-collision
    # check ran outside the lock so two concurrent creates with the
    # same name could both pass and both append.
    with _state._cache_lock:
        projs = _load_experiment_projects()
        if any((p.get("name") or "").strip() == name for p in projs):
            return ({"error":
                      f"project named {name!r} already exists"}, 409)
        projs.append({
            "name":        name,
            "description": desc,
            "experiments": [],
            "saved":       _date.today().isoformat(),
        })
        err = _agent_save_or_500(
            lambda: _save_experiment_projects(projs),
            "Experiment projects",
        )
        if err:
            return err
    _log_event("project.created", name=name, via="agent")
    return {"ok": True, "name": name}


@_agent_endpoint("list-gels")
def _h_list_gels(app, payload):
    """List every saved gel snapshot. Returns id, name, lane count,
    agarose %, and timestamps per gel (omits per-lane detail — fetch
    via ``get-gel``)."""
    out: "list[dict]" = []
    for g in _load_gels():
        lanes = g.get("lanes") or []
        n_lanes = len(lanes) if isinstance(lanes, list) else 0
        out.append({
            "id":          g.get("id", ""),
            "name":        g.get("name", ""),
            "n_lanes":     n_lanes,
            "agarose_pct": g.get("agarose_pct", 1.0),
            "created_at":  g.get("created_at", ""),
            "updated_at":  g.get("updated_at", ""),
        })
    return {"gels": out}


@_agent_endpoint("get-gel")
def _h_get_gel(app, payload):
    """Fetch one saved gel snapshot (full lane payload + notes)."""
    gid = _sanitize_gel_id(payload.get("id"))
    if gid is None:
        return ({"error": "missing or invalid 'id'"}, 400)
    entry = _find_gel(gid)
    if entry is None:
        return ({"error": f"no gel with id {gid!r}"}, 404)
    return {"gel": _typed_clone(entry)}


@_agent_endpoint("create-gel", write=True)
def _h_create_gel(app, payload):
    """Save a new gel snapshot. Body:
    ``{name: str, lanes: list[dict], agarose_pct?: float = 1.0,
        notes?: str = ""}``.
    Routes through `_normalise_gel_entry` (caps + clamps). Returns
    the freshly-stamped id."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing 'name'"}, 400)
    name = name.strip()
    lanes = payload.get("lanes") or []
    if not isinstance(lanes, list):
        return ({"error": "'lanes' must be a list of dicts"}, 400)
    # Sweep #26 (2026-05-25): hold `_state._cache_lock` across the load +
    # mutate + save so concurrent agent calls can't both pass the
    # name-collision check and both append. Pre-fix two simultaneous
    # `create-gel` calls each loaded the same snapshot, each appended,
    # each saved — the second save's `_typed_clone(entries)` overwrote
    # the first's gel from the cache. RLock allows the nested re-entry
    # from `_save_gels`'s own locking.
    with _state._cache_lock:
        gels = _load_gels()
        if any((g.get("name") or "").strip() == name for g in gels):
            return ({"error": f"gel named {name!r} already exists"}, 409)
        existing_gel_ids: set[str] = {
            g.get("id") for g in gels if g.get("id")  # type: ignore[misc]
        }
        new_id = _new_gel_id(existing_gel_ids)
        entry = _normalise_gel_entry({
            "id":          new_id,
            "name":        name,
            "lanes":       lanes,
            "agarose_pct": payload.get("agarose_pct", 1.0),
            "notes":       payload.get("notes", ""),
        }, fresh=True)
        gels.append(entry)
        err = _agent_save_or_500(lambda: _save_gels(gels), "Gels")
    if err:
        return err
    _log_event("gel.created", gid=new_id, name=name, via="agent")
    return {"ok": True, "id": new_id, "gel": _typed_clone(entry)}


@_agent_endpoint("update-gel", write=True)
def _h_update_gel(app, payload):
    """Replace a saved gel's fields by id. Body:
    ``{id: str, name?: str, lanes?: list[dict],
        agarose_pct?: float, notes?: str}``. Fields not supplied
    keep their prior value. Re-normalises through `_normalise_gel_entry`."""
    gid = _sanitize_gel_id(payload.get("id"))
    if gid is None:
        return ({"error": "missing or invalid 'id'"}, 400)
    # Sweep #26: RMW under `_state._cache_lock` so concurrent renames /
    # field-merges don't drop each other's updates.
    with _state._cache_lock:
        gels = _load_gels()
        for i, g in enumerate(gels):
            if g.get("id") != gid:
                continue
            merged = dict(g)
            if "name" in payload:
                nm = payload["name"]
                if not isinstance(nm, str) or not nm.strip():
                    return ({"error":
                              "'name' must be a non-empty string"}, 400)
                new_name = nm.strip()
                if new_name != merged.get("name") and any(
                    (other.get("name") or "").strip() == new_name
                    for j, other in enumerate(gels) if j != i
                ):
                    return ({"error":
                              f"gel named {new_name!r} already exists"}, 409)
                merged["name"] = new_name
            if "lanes" in payload:
                ln = payload["lanes"]
                if not isinstance(ln, list):
                    return ({"error": "'lanes' must be a list"}, 400)
                merged["lanes"] = ln
            if "agarose_pct" in payload:
                merged["agarose_pct"] = payload["agarose_pct"]
            if "notes" in payload:
                nt = payload["notes"]
                merged["notes"] = nt if isinstance(nt, str) else ""
            gels[i] = _normalise_gel_entry(merged, fresh=False)
            err = _agent_save_or_500(lambda: _save_gels(gels), "Gels")
            if err:
                return err
            _log_event("gel.renamed", gid=gid, via="agent") \
                if "name" in payload else None
            return {"ok": True, "gel": _typed_clone(gels[i])}
    return ({"error": f"no gel with id {gid!r}"}, 404)


@_agent_endpoint("delete-gel", write=True)
def _h_delete_gel(app, payload):
    """Delete a saved gel snapshot by id."""
    gid = _sanitize_gel_id(payload.get("id"))
    if gid is None:
        return ({"error": "missing or invalid 'id'"}, 400)
    # Sweep #26: RMW under `_state._cache_lock` so a concurrent create-gel
    # can't insert into the snapshot we read but didn't persist.
    with _state._cache_lock:
        gels = _load_gels()
        new_gels = [g for g in gels if g.get("id") != gid]
        if len(new_gels) == len(gels):
            return ({"error": f"no gel with id {gid!r}"}, 404)
        err = _agent_save_or_500(
            lambda: _save_gels(new_gels), "Gels",
        )
    if err:
        return err
    _log_event("gel.deleted", gid=gid, via="agent")
    return {"ok": True, "id": gid, "remaining": len(new_gels)}


@_agent_endpoint("list-protein-motifs")
def _h_list_protein_motifs(app, payload):
    """List the merged protein-motif library (built-ins + user
    overrides). Each entry carries name, feature_type, sequence,
    color, description. Sourced from `_load_protein_motifs`."""
    return {"motifs": [
        {
            "name":         m.get("name", ""),
            "feature_type": m.get("feature_type", "Motif"),
            "sequence":     m.get("sequence", ""),
            "color":        m.get("color", ""),
            "description":  m.get("description", ""),
        }
        for m in _load_protein_motifs()
    ]}


@_agent_endpoint("set-protein-motif", write=True)
def _h_set_protein_motif(app, payload):
    """Create or override a protein motif. Body:
    ``{name: str, sequence: str, feature_type?: str = "Motif",
        color?: str = "", description?: str = ""}``.

    Persists only the user-modified entries (built-ins stay in code);
    `_load_protein_motifs` merges on read. If `name` matches a
    built-in, the user override replaces it in the merged list."""
    raw_name = payload.get("name")
    if not isinstance(raw_name, str) or not raw_name.strip():
        return ({"error": "missing 'name'"}, 400)
    # Sweep #26: defense-in-depth — route name + feature_type +
    # description + color through `_sanitize_label` so an agent
    # cannot smuggle terminal escape codes into the on-disk JSON.
    name = _sanitize_label(raw_name.strip(), max_len=200)
    if not name:
        return ({"error": "missing 'name'"}, 400)
    seq = payload.get("sequence")
    if not isinstance(seq, str) or not seq.strip():
        return ({"error": "missing 'sequence'"}, 400)
    seq = seq.strip().upper()
    # Same single-letter AA validation as `_h_optimize_protein`.
    invalid = [c for c in seq if c not in "ACDEFGHIKLMNPQRSTVWY*"]
    if invalid:
        return ({"error":
                  "non-canonical amino acids in 'sequence': "
                  f"{''.join(sorted(set(invalid)))!r}"}, 400)
    raw_ftype = payload.get("feature_type") or "Motif"
    ftype = (
        _sanitize_label(raw_ftype, max_len=200) or "Motif"
        if isinstance(raw_ftype, str) else "Motif"
    )
    raw_color = payload.get("color") or ""
    color = _sanitize_label(raw_color, max_len=64) if isinstance(raw_color, str) else ""
    raw_desc = payload.get("description") or ""
    desc = _sanitize_label(raw_desc, max_len=2000) if isinstance(raw_desc, str) else ""
    # Sweep #26: RMW under `_state._cache_lock` so concurrent set-protein-motif
    # calls can't both load the user file, both append/override, and
    # both save — second save wins.
    with _state._cache_lock:
        # Load existing user file directly (NOT the merged list, which
        # includes built-ins). `_load_protein_motifs` merges on read, so
        # we re-read the raw user file to know which entries to persist.
        try:
            user_entries, _ = _safe_load_json(
                _state._PROTEIN_MOTIFS_FILE, "Protein motifs",
            )
        except Exception as exc:
            _log.exception("agent set-protein-motif: load failed")
            return ({"error": f"load failed: {_scrub_path(str(exc))}"}, 500)
        user_entries = [
            e for e in user_entries if isinstance(e, dict)
        ]
        # Drop any existing entry with the same name (copy-on-write).
        user_entries = [
            e for e in user_entries if e.get("name") != name
        ]
        user_entries.append({
            "name":         name,
            "feature_type": ftype,
            "sequence":     seq,
            "color":        color,
            "description":  desc,
        })
        err = _agent_save_or_500(
            lambda: _save_protein_motifs(user_entries),
            "Protein motifs",
        )
        if err:
            return err
    _log_event("synthesis.protein.motif_edit", name=name, via="agent")
    return {"ok": True, "name": name}


@_agent_endpoint("delete-protein-motif", write=True)
def _h_delete_protein_motif(app, payload):
    """Delete a user-stored protein motif override. Built-in
    motifs cannot be deleted — only user overrides. If `name`
    matches a built-in but the user has not edited it, returns 404
    (nothing to delete). If the user has overridden a built-in,
    deleting the override restores the original built-in entry."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing 'name'"}, 400)
    name = name.strip()
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        try:
            user_entries, _ = _safe_load_json(
                _state._PROTEIN_MOTIFS_FILE, "Protein motifs",
            )
        except Exception as exc:
            _log.exception("agent delete-protein-motif: load failed")
            return ({"error": f"load failed: {_scrub_path(str(exc))}"}, 500)
        user_entries = [
            e for e in user_entries if isinstance(e, dict)
        ]
        new_entries = [e for e in user_entries if e.get("name") != name]
        if len(new_entries) == len(user_entries):
            return ({"error":
                      f"no user-stored override for {name!r} "
                      "(built-ins cannot be deleted)"}, 404)
        err = _agent_save_or_500(
            lambda: _save_protein_motifs(new_entries),
            "Protein motifs",
        )
        if err:
            return err
    _log_event(
        "synthesis.protein.motif_edit",
        name=name, action="delete", via="agent",
    )
    return {"ok": True, "name": name}


_ENZYME_IUPAC_ALPHABET = frozenset("ACGTRYSWKMBDHVN")


def _agent_validate_custom_enzyme_payload(payload: dict) -> "dict | str":
    """Validate + canonicalise an agent-supplied custom-enzyme dict.
    Returns the persistable dict on success or an error string on
    failure — mirrors `AddCustomEnzymeModal._validate` so the agent
    surface accepts the same shape the UI does."""
    name = _sanitize_label(payload.get("name"), max_len=200)
    if not name:
        return "missing or non-string 'name'"
    if not (1 <= len(name) <= 64):
        return "'name' must be 1-64 characters"
    site = payload.get("site")
    if not isinstance(site, str) or not site.strip():
        return "missing or non-string 'site'"
    site = site.strip().upper()
    if not (4 <= len(site) <= 30):
        return "'site' must be 4-30 characters"
    bad = set(site) - _ENZYME_IUPAC_ALPHABET
    if bad:
        return f"'site' has non-IUPAC characters: {''.join(sorted(bad))!r}"
    fwd_raw = payload.get("fwd_cut")
    rev_raw = payload.get("rev_cut")
    if fwd_raw is None or rev_raw is None:
        return "missing 'fwd_cut' or 'rev_cut'"
    try:
        fwd_cut = int(fwd_raw)
        rev_cut = int(rev_raw)
    except (TypeError, ValueError):
        return "'fwd_cut' and 'rev_cut' must be integers"
    # Sweep #26: shared `_ENZYME_CUT_RANGE` constant — was hardcoded
    # ±30 in two places.
    lo, hi = -_ENZYME_CUT_RANGE, len(site) + _ENZYME_CUT_RANGE
    if not (lo <= fwd_cut <= hi) or not (lo <= rev_cut <= hi):
        return f"cut positions must be in {lo}..{hi}"
    ftype = _sanitize_label(payload.get("type"), max_len=64) or "other"
    supplier = _sanitize_label(payload.get("supplier"), max_len=64)
    return {
        "name":     name,
        "site":     site,
        "fwd_cut":  fwd_cut,
        "rev_cut":  rev_cut,
        "type":     ftype,
        "supplier": supplier,
    }


@_agent_endpoint("list-custom-enzymes")
def _h_list_custom_enzymes(app, payload):
    """List every user-added custom enzyme. Built-in NEB enzymes are
    NOT included — fetch the combined view via list-restriction-sites
    or `_state._all_enzymes_hook()` introspection. Each item carries name, site,
    fwd_cut, rev_cut, type, supplier."""
    return {"ok": True, "enzymes": _load_custom_enzymes()}


@_agent_endpoint("get-custom-enzyme")
def _h_get_custom_enzyme(app, payload):
    """Return a single custom enzyme by `name`. Body: ``{name}``."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing or non-string 'name'"}, 400)
    name = name.strip()
    meta = _custom_enzyme_meta(name)
    if meta is None:
        return ({"error": f"unknown custom enzyme {name!r}"}, 404)
    return {"ok": True, "enzyme": meta}


@_agent_endpoint("create-custom-enzyme", write=True)
def _h_create_custom_enzyme(app, payload):
    """Add a new custom enzyme. Body:
    ``{name, site, fwd_cut, rev_cut, type?, supplier?}``. Returns 409
    if `name` collides with an existing built-in OR custom enzyme.
    Mirrors `AddCustomEnzymeModal` validation rules.

    Sweep #25 (2026-05-23) — Built-in NEB enzyme names are
    intentionally reserved: ``create-custom-enzyme`` refuses any
    name in ``_state._all_enzymes_hook()``, and ``update-custom-enzyme`` only
    matches existing CUSTOM entries (not built-ins). The earlier
    "add a custom enzyme with the same name to override" note
    in this docstring was aspirational — the actual implementation
    treats built-in names as read-only. Users wanting an alternate
    cut convention for a known enzyme name should append a suffix
    (e.g. ``BsaI-isoschiz``) instead.

    Sweep #25 (2026-05-23): collision-check + load + append + save
    wrapped in ``_state._cache_lock`` so two concurrent agent calls can't
    both pass the collision check on an identical cache snapshot and
    then both append (silently dropping the first writer's entry on
    the second save). RLock allows re-entry from ``_save_custom_enzymes``.
    Matches the worker-side RMW pattern from INV-51.
    """
    payload_or_err = _agent_validate_custom_enzyme_payload(payload)
    if isinstance(payload_or_err, str):
        return ({"error": payload_or_err}, 400)
    with _state._cache_lock:
        if payload_or_err["name"] in _state._all_enzymes_hook():
            return ({"error":
                      f"enzyme {payload_or_err['name']!r} already exists; "
                      "use update-custom-enzyme to modify"}, 409)
        entries = _load_custom_enzymes()
        entries.append(payload_or_err)
        if (err := _agent_save_or_500(
                lambda: _save_custom_enzymes(entries),
                "Custom enzymes")) is not None:
            return err
    _log_event("custom_enzyme.added",
                name=payload_or_err["name"], via="agent")
    return {"ok": True, "name": payload_or_err["name"]}


@_agent_endpoint("update-custom-enzyme", write=True)
def _h_update_custom_enzyme(app, payload):
    """Replace an existing custom enzyme by `name`. Built-in NEB
    enzymes are not editable via this endpoint (refuse with 400) —
    add a custom enzyme with the same name to override.

    Sweep #25: full RMW under ``_state._cache_lock`` (see
    ``_h_create_custom_enzyme`` docstring for rationale).
    """
    payload_or_err = _agent_validate_custom_enzyme_payload(payload)
    if isinstance(payload_or_err, str):
        return ({"error": payload_or_err}, 400)
    name = payload_or_err["name"]
    with _state._cache_lock:
        entries = _load_custom_enzymes()
        for i, e in enumerate(entries):
            if e.get("name") == name:
                entries[i] = payload_or_err
                if (err := _agent_save_or_500(
                        lambda: _save_custom_enzymes(entries),
                        "Custom enzymes")) is not None:
                    return err
                return {"ok": True, "name": name}
    return ({"error": f"unknown custom enzyme {name!r}"}, 404)


@_agent_endpoint("delete-custom-enzyme", write=True)
def _h_delete_custom_enzyme(app, payload):
    """Delete a custom enzyme by `name`. Built-in NEB enzymes are
    refused. Any enzyme collection that referenced the deleted name
    keeps the stale row — `_active_enzyme_allowed_set` filters
    against `_state._all_enzymes_hook()` so the stale name is dropped at
    scan time without separate housekeeping.

    Sweep #25: full RMW under ``_state._cache_lock`` (see
    ``_h_create_custom_enzyme`` docstring for rationale).
    """
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing or non-string 'name'"}, 400)
    name = name.strip()
    with _state._cache_lock:
        entries = _load_custom_enzymes()
        new_entries = [e for e in entries if e.get("name") != name]
        if len(new_entries) == len(entries):
            return ({"error": f"unknown custom enzyme {name!r}"}, 404)
        if (err := _agent_save_or_500(
                lambda: _save_custom_enzymes(new_entries),
                "Custom enzymes")) is not None:
            return err
    _log_event("custom_enzyme.deleted", name=name, via="agent")
    return {"ok": True, "name": name}


@_agent_endpoint("list-enzyme-collections")
def _h_list_enzyme_collections(app, payload):
    """List every enzyme collection (named subset of the master
    catalog). Each item carries ``{name, enzymes: [str, ...]}``."""
    return {"ok": True, "collections": _load_enzyme_collections()}


@_agent_endpoint("create-enzyme-collection", write=True)
def _h_create_enzyme_collection(app, payload):
    """Create a new enzyme collection. Body:
    ``{name, enzymes?: [str, ...] = []}``. Returns 409 if a collection
    with the same name already exists. Unknown enzyme names are
    accepted at create time (they're filtered at scan time by
    `_active_enzyme_allowed_set` so the user can add custom enzymes
    later and have them participate retroactively)."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing or non-string 'name'"}, 400)
    name = name.strip()
    enzymes = payload.get("enzymes") or []
    if not isinstance(enzymes, list) or not all(
            isinstance(e, str) for e in enzymes):
        return ({"error": "'enzymes' must be a list of strings"}, 400)
    # Sweep #25: full RMW under `_state._cache_lock` (see
    # `_h_create_custom_enzyme` docstring for rationale).
    with _state._cache_lock:
        entries = _load_enzyme_collections()
        if any(e.get("name") == name for e in entries):
            return ({"error":
                      f"enzyme collection {name!r} already exists; "
                      "use update-enzyme-collection to modify"}, 409)
        entries.append({"name": name, "enzymes": sorted(set(enzymes))})
        if (err := _agent_save_or_500(
                lambda: _save_enzyme_collections(entries),
                "Enzyme collections")) is not None:
            return err
    return {"ok": True, "name": name}


@_agent_endpoint("clear-entry-vectors-for-grammar", write=True)
def _h_clear_entry_vectors_for_grammar(app, payload):
    """Drop every entry-vector binding for one grammar id. Body:
    ``{grammar_id: str}``. Used by grammar-delete flows; agents
    rarely need this directly but it round-trips for symmetry."""
    gid = payload.get("grammar_id")
    if not isinstance(gid, str) or not gid:
        return ({"error": "missing 'grammar_id'"}, 400)
    try:
        n = _clear_entry_vectors_for_grammar(gid)
    except (OSError, RuntimeError) as exc:
        return ({"error": f"clear failed: {_scrub_path(str(exc))}"}, 500)
    _log_event(
        "entry_vector.cleared",
        grammar=gid, n_cleared=n, via="agent",
    )
    return {"ok": True, "grammar_id": gid, "n_cleared": n}


# Deferred data handlers (Phase D: relocated-blocker + hooked-engine endpoints)
# Hard cap on per-side sequence length to keep DP memory bounded.
# 200 kb on each side = ~80 GB if naive O(NM); Biopython's aligner is
# smarter than that but we still don't want the user pasting a
# chromosome.
_PAIRWISE_MAX_LEN = 200_000


# Output formats accepted by the bulk-export helper. Each maps to a
# writer function + canonical extension. Keep in sync with the format
# Select in `BulkExportCollectionModal` and the agent endpoint's
# `format` validation.
_BULK_EXPORT_FORMATS: dict = {
    "genbank": ("gb",   "GenBank"),
    "embl":    ("embl", "EMBL"),
    "fasta":   ("fa",   "FASTA"),
    "dna":     ("dna",  "CommercialSaaS .dna"),
    "png":     ("png",  "Map image (PNG)"),
    "svg":     ("svg",  "Map image (SVG)"),
}


# Module-level cache globals that hold cleared persisted state. Every
# entry MUST exist as `_<name>_cache: ... | None = None` somewhere in
# this module — `_perform_master_delete` looks them up by name and
# resets to None so the next read re-loads from (now empty) disk.
# Adding a new cache? Append it here AND to `_protect_user_data` in
# tests/conftest.py.
_MASTER_DELETE_CACHE_ATTRS: tuple = (
    "_library_cache",
    "_collections_cache",
    "_parts_bin_cache",
    "_parts_bin_collections_cache",
    "_primers_cache",
    "_primer_collections_cache",
    "_primer_usage_cache",            # sweep #20: was missed, primer use-count map
    "_features_cache",
    "_feature_colors_cache",
    "_feature_library_index_cache",   # sweep #20: was missed, feature-library lookup index
    "_grammars_cache",
    "_entry_vectors_cache",
    "_codon_tables_cache",
    "_settings_cache",
    "_experiments_cache",
    "_experiment_projects_cache",
    "_gels_cache",
    "_protein_motifs_cache",   # sweep #15: protein-motif user overrides
    "_custom_enzymes_cache",      # user-added restriction enzymes
    "_enzyme_collections_cache",  # enzyme-set collections
    "_hmm_db_catalog_cache",      # sweep #28: HMM database registry
    "_protein_collections_cache",  # named protein-sequence collections (operon workbench)
    "_model_collections_cache",   # BABS model-picker collections (INV-139)
    "_protocol_collections_cache",  # OT-2 saved protocol designs
    "_custom_labware_cache",      # OT-2 custom labware definitions
)


# Per-lane source kinds. The UI surfaces these in a Select dropdown.
_LANE_SOURCES: tuple[tuple[str, str], ...] = (
    ("Empty",          "empty"),
    ("Ladder",         "ladder"),
    ("Plasmid (uncut)", "plasmid"),
    ("Digest",         "digest"),
    ("PCR amplicon",   "pcr"),
)


@_agent_endpoint("tools")
def _h_tools(app, payload):
    """Self-describe: every endpoint with its read/write mode and full
    request-body documentation.

    Each entry is ``{name, method, write, doc, doc_full}`` — ``doc`` is
    the one-line summary; ``doc_full`` (snag #19) is the COMPLETE
    docstring, which documents the request body (required / optional
    keys, aliases, enums, size caps) so a caller forms a correct request
    in one round-trip instead of N trial-and-error 400s.

    Global error-shape contract every endpoint follows:
      * Success     → ``200`` with handler-specific JSON.
      * Bad input   → ``400`` with ``{"error": "..."}`` describing the
                       field that failed validation.
      * Auth        → ``401`` for write endpoints called without the
                       bearer token.
      * Missing     → ``404`` when the requested object (library entry,
                       collection, hmm file) doesn't exist.
      * Stale ref   → ``409`` when the request was valid at receive time
                       but the canvas record / library state moved on
                       before apply (audit sweep #4 added this to
                       `transfer-annotations`; mirrors
                       `replace-sequence`). Retry against current state.
      * Loaded?     → ``422`` for endpoints that need a record loaded
                       but `_current_record is None`.
      * Save fail   → ``500`` with ``{"error": "save failed for X: ..."}``
                       when the disk write raised (disk-full, RO mount,
                       EACCES). The in-process UI user ALSO sees a
                       failure toast via `_notify_save_failure`. Agent
                       endpoints route through `_agent_save_or_500` so
                       the shape is uniform across 7 write paths
                       (delete-from-library, create/delete/rename-
                       collection, set-active-collection,
                       bulk-import-folder, set-plasmid-status).
      * Crash       → ``500`` with ``{"error": "...", "type": "..."}``
                       (the dispatcher's catch-all).
    """
    # `/tools` is served unauthenticated and AHEAD of the dispatcher's
    # envelope wrapper (which adds the universal `data`/`_stale` fields),
    # so it deliberately returns a bare ``{endpoints: [...]}`` — no `data`
    # mirror. It is ONE of the two pre-envelope bare responses (the other is
    # `/healthz`, the readiness probe); the "every response carries `data`"
    # contract holds only for AUTHENTICATED endpoints (docs/agent-api.md).
    # Duplicating this response's `doc_full` corpus under `data` would double
    # the single largest payload in the API (~150 KB → ~300 KB) for zero
    # ergonomic gain — a client calling the discovery endpoint reads
    # `endpoints` by name.
    return {"endpoints": [
        {
            "name":   name,
            "method": "POST" if write else "GET",
            "write":  write,
            "doc":    (fn.__doc__ or "").strip().split("\n")[0],
            # Snag #19: the FULL docstring, not just its first line. Each
            # handler's body shape (required/optional keys, aliases,
            # enums, caps) is documented in prose here — emitting it lets
            # an agent get the request schema in ONE self-describe call
            # instead of N trial-and-error 400s. `inspect.getdoc` strips
            # the uniform leading indentation so it reads cleanly.
            "doc_full": inspect.getdoc(fn) or "",
        }
        for name, (fn, write) in sorted(_state._AGENT_HANDLERS.items())
    ]}


_EXPORT_GENBANK_EXTS = (".gb", ".gbk", ".genbank")


_EXPORT_GFF_EXTS     = (".gff", ".gff3")


_EXPORT_FASTA_EXTS   = (".fa", ".fasta", ".fna")


@_agent_endpoint("export-genbank", write=True)
def _h_export_genbank(app, payload):
    """Write the current record to `path` as GenBank. Body: ``{path}``.
    Path may be relative (resolved against $HOME). Doesn't change the
    record's `_source_path`, so subsequent `save` calls still target
    the original location."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    path = _sanitize_path(payload.get("path"))
    if path is None:
        return ({"error": "missing 'path'"}, 400)
    if (ext_err := _check_export_extension(
            path, _EXPORT_GENBANK_EXTS, "GenBank")) is not None:
        return ({"error": ext_err}, 400)
    err = _check_agent_write_path(path)
    if err is not None:
        return ({"error": err}, 403)
    try:
        result = _export_genbank_to_path(rec, path)
    except (OSError, ValueError) as exc:
        return ({"error": f"export failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, **result}


@_agent_endpoint("export-gff", write=True)
def _h_export_gff(app, payload):
    """Write the current record to `path` as GFF3. Body: ``{path}``.
    Wrap features become two GFF3 rows joined by a shared `ID=...`
    attribute. Returns ``{ok, path, bp, features}``."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    path = _sanitize_path(payload.get("path"))
    if path is None:
        return ({"error": "missing 'path'"}, 400)
    if (ext_err := _check_export_extension(
            path, _EXPORT_GFF_EXTS, "GFF3")) is not None:
        return ({"error": ext_err}, 400)
    err = _check_agent_write_path(path)
    if err is not None:
        return ({"error": err}, 403)
    try:
        result = _export_gff_to_path(rec, path)
    except OSError as exc:
        return ({"error": f"export failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, **result}


@_agent_endpoint("export-fasta", write=True)
def _h_export_fasta(app, payload):
    """Write the current sequence to `path` as FASTA. Body: ``{path}``.
    Single-record FASTA with the plasmid name as the header."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    path = _sanitize_path(payload.get("path"))
    if path is None:
        return ({"error": "missing 'path'"}, 400)
    if (ext_err := _check_export_extension(
            path, _EXPORT_FASTA_EXTS, "FASTA")) is not None:
        return ({"error": ext_err}, 400)
    err = _check_agent_write_path(path)
    if err is not None:
        return ({"error": err}, 403)
    try:
        # The DISPLAY name, not the underscored LOCUS slug — a FASTA header
        # has neither the 16-character limit nor the no-space rule. [INV-98]
        result = _export_fasta_to_path(
            (getattr(rec, "_tui_display_name", None)
             or rec.name or rec.id or "plasmid"),
            str(rec.seq),
            path,
        )
    except (OSError, ValueError) as exc:
        return ({"error": f"export failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, **result}


@_agent_endpoint("bulk-export-collection", write=True)
def _h_bulk_export_collection(app, payload):
    """Write every plasmid in `collection` to `dir` as files of
    `format`. Body: ``{collection, format, dir}``. ``format`` ∈
    {"genbank", "embl", "fasta", "dna"}. Returns ``{ok, total,
    written, failures, dir}`` mirroring `_bulk_export_collection`.
    Per-entry failures don't abort the batch — each is reported in
    the `failures` list. Path is sanitised + symlink-refused via the
    same checks `export-genbank` uses on the target dir."""
    coll = payload.get("collection")
    if not coll or not isinstance(coll, str):
        return ({"error": "missing 'collection'"}, 400)
    fmt = payload.get("format")
    if fmt not in _BULK_EXPORT_FORMATS:
        return ({"error": (
            f"unknown format {fmt!r}; "
            f"choose one of {sorted(_BULK_EXPORT_FORMATS)}"
        )}, 400)
    target = _sanitize_path(payload.get("dir"))
    if target is None:
        return ({"error": "missing 'dir'"}, 400)
    if target.exists() and not target.is_dir():
        return ({"error": f"target is not a directory: {target}"}, 400)
    err = _check_agent_write_path(target)
    if err is not None:
        return ({"error": err}, 403)
    try:
        result = _state._bulk_export_collection_hook(coll, fmt, target)
    except (OSError, ValueError) as exc:
        return ({"error": f"bulk export failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, **result}


@_agent_endpoint("find-sequence")
def _h_find_sequence(app, payload):
    """Locate a DNA subsequence / IUPAC motif on the loaded record. Body:
    ``{query, both_strands?, max_hits?}``.

    * ``query`` — ACGT/U or IUPAC bases (e.g. ``GAATTC``, ``GGWCC``,
      ``NNNGATCNNN``); case-insensitive, U treated as T.
    * ``both_strands`` — also search the reverse-complement strand (default
      ``true``); reverse hits report the FORWARD coordinate with ``strand: -1``
      (sacred invariant #2). A palindromic match is reported ONCE, on the
      forward strand (mirrors invariant #1).
    * ``max_hits`` — cap on returned hits (default 100, max 10000).

    Wrap-aware on a circular plasmid (a motif spanning the origin is found).
    Matches are non-overlapping (``re.finditer``). Returns ``hits:
    [{start, end, strand, match, wraps}]`` — 0-based half-open ``[start, end)``
    in FORWARD coordinates; ``end <= start`` with ``wraps: true`` when the hit
    crosses the origin. Read-only — does not change the record."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    raw = payload.get("query")
    if not isinstance(raw, str) or not raw.strip():
        return ({"error": "missing 'query' (ACGT / IUPAC bases)"}, 400)
    q = raw.strip().upper().replace("U", "T")
    if not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", q):
        return ({"error": "'query' must be ACGT / IUPAC bases"}, 400)
    seq = str(getattr(rec, "seq", "") or "").upper()
    n = len(seq)
    if n == 0 or len(q) > n:
        return {"ok": True, "query": q, "n_hits": 0, "hits": []}
    both = payload.get("both_strands", True)
    if not isinstance(both, bool):
        both = str(both).strip().lower() not in ("false", "0", "no", "off", "")
    max_hits = _coerce_int(payload.get("max_hits", 100), name="max_hits")
    if isinstance(max_hits, str):
        return ({"error": max_hits}, 400)
    max_hits = max(1, min(int(max_hits), 10000))
    circular = (rec.annotations.get("topology") or "linear") == "circular"
    extra = (len(q) - 1) if circular else 0
    search = (seq + seq[:extra]) if extra else seq
    hits: "list[dict]" = []
    seen: "set[tuple[int, int]]" = set()

    def _scan(pattern, strand) -> bool:
        for m in pattern.finditer(search):
            s = m.start()
            if s >= n:                       # duplicate from the wrap region
                break
            e = s + len(q)
            wraps = e > n
            end = (e - n) if wraps else e
            key = (s, end)
            if key in seen:                  # palindrome already found forward
                continue
            seen.add(key)
            hits.append({"start": s, "end": end, "strand": strand,
                         "match": m.group(), "wraps": wraps})
            if len(hits) >= max_hits:
                return True
        return False

    full = _scan(_iupac_pattern(q), 1)
    if both and not full:
        _scan(_iupac_pattern(_rc(q)), -1)
    hits.sort(key=lambda h: (h["start"], -h["strand"]))
    return {"ok": True, "query": q, "circular": circular,
            "n_hits": len(hits), "hits": hits[:max_hits]}


@_agent_endpoint("export-migrate-archive", write=True)
def _h_export_migrate_archive(app, payload):
    """Package ALL user data into one portable, compressed ``.zip`` (backup, or
    moving to another machine). Body: ``{path, include_hmm?}``. ``path`` must end
    in ``.zip``; ``include_hmm`` (default false) bundles the large HMM databases
    too (re-derivable from the always-included catalog otherwise).

    READ-ONLY with respect to the data dir — it only copies user data OUT, never
    modifies it. The inverse (IMPORT) is deliberately NOT exposed over the agent
    API: it replaces ALL data, the same risk class as Master Delete. Returns
    ``{path, bytes, n_files, n_dirs, included_hmm}``."""
    path = _sanitize_path(payload.get("path"))
    if path is None:
        return ({"error": "missing 'path'"}, 400)
    if (ext_err := _check_export_extension(
            path, (".zip",), "migrate archive")) is not None:
        return ({"error": ext_err}, 400)
    err = _check_agent_write_path(path)
    if err is not None:
        return ({"error": err}, 403)
    include_hmm = bool(payload.get("include_hmm"))
    try:
        result = _export_migrate_archive(path, include_hmm=include_hmm)
    except OSError as exc:
        return ({"error": f"export failed: {_scrub_path(str(exc))}"}, 500)
    return {"ok": True, **result,
            "ignored": _agent_ignored_keys(payload, {"path", "include_hmm"})}


@_agent_endpoint("list-restriction-sites")
def _h_list_restriction_sites(app, payload):
    """Scan the loaded record for restriction sites. Body:
    ``{enzymes?: [str, ...], min_length?: int, unique_only?: bool,
       respect_active_collection?: bool = true}``.

    Default scans the combined catalog (built-in NEB ∪ user-added
    custom enzymes) with `min_length=4` and `unique_only=False`.

    When ``enzymes`` is omitted AND ``respect_active_collection`` is
    true (the default), the agent endpoint mirrors the UI overlay: if
    the user has set an active enzyme collection, only those enzymes
    surface in the result. Pass ``respect_active_collection=false`` to
    force a full-catalog scan regardless. Explicit ``enzymes`` always
    wins.

    Returns each hit as ``{enzyme, start, end, strand, cut_bp,
    bottom_cut_bp, site, wraps}`` — `cut_bp`/`bottom_cut_bp` are absolute
    top-strand coordinates of the two strand cuts, `site` the recognition
    sequence, so fragment sizes can be checked without a private table of
    cut offsets.

    **An unknown enzyme name is a 400, never an empty result.** Names are
    matched case-insensitively against the combined catalog; anything that
    doesn't match comes back as `unknown_enzymes` with near-miss
    suggestions. This matters because the failure it replaces was
    wrong-answer-shaped: a typo used to return `count: 0`, which reads as
    "no sites here" rather than "I never looked".

    **Isoschizomers resolve by the name you asked for.** The scan collapses
    enzymes that share a recognition site and cut (BsmBI / Esp3I / BsmBI-v2,
    BsaI / BspTNI, EcoRI / EcoRI-HF) onto ONE label so the map stays
    readable; asking for a name on the losing side of that collapse used to
    return nothing. Named enzymes now go into the scan as the hand-picked
    set, and every hit is reported under each requested spelling.
    `equivalent_enzymes` lists the names that cut identically but weren't
    asked for.

    **Commercial synonyms resolve too.** The same enzyme is sold under
    different names by different suppliers — Thermo's `Eco31I` is NEB's
    `BsaI`, `LguI` is `SapI`, `AarI` is `PaqCI` — and the scan now accepts
    either. Hits come back labelled with the name YOU asked for.

    **NEOschizomers stay apart.** XmaI (`C^CCGGG`) and SmaI (`CCC^GGG`) read
    the same six bases and leave different ends, so they are NOT
    interchangeable and this endpoint reports both — even though the plasmid
    map deliberately draws one bar for the pair. Same for Acc65I/KpnI,
    ApaI/PspOMI, NheI/BmtI, KasI/NarI, AatII/ZraI, SacI/Eco53kI.

    Note that `min_length` does not apply to an explicit `enzymes` list (a
    hand-picked short cutter is kept) — passing both returns a `warnings`
    entry saying so."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    enzymes = payload.get("enzymes")
    # Accept the singular `enzyme` as an alias and a bare string under
    # either key (agent-API feedback P1-3): the natural first guess is
    # `{enzyme: "BsaI"}`, and silently dropping it (then scanning the whole
    # catalog) is the SC-D "ignored routing param" footgun — the result
    # looks authoritative while answering a different question.
    if enzymes is None and payload.get("enzyme") not in (None, ""):
        enzymes = payload.get("enzyme")
    if isinstance(enzymes, str):
        enzymes = [enzymes]
    resolved: "list[str]" = []
    if enzymes is not None:
        # Shape, size and element type up front. A mixed-type list (e.g.
        # `[1, 2.5, null]`) used to build a set whose `not in` check silently
        # filtered every hit to zero; the size cap bounds the name resolution
        # below.
        if (shape_err := _agent_check_enzyme_list(enzymes)) is not None:
            return shape_err
        # Every requested name must BE an enzyme. Pre-fix an unrecognised
        # name — a typo, or an isoschizomer spelling the scan files under a
        # different label — filtered every hit away and returned
        # `{"sites": [], "count": 0}`, so "check this construct carries no
        # Esp3I site" answered "clean" without ever having scanned for one.
        # A safety check that cannot fail is worse than no check at all, so
        # an unknown name is a 400 now. An empty list is refused for the
        # same reason: it reads as "no enzymes matched", not "scan
        # everything".
        resolved, unknown = _resolve_enzyme_names(enzymes)
        if unknown:
            return _agent_unknown_enzyme_error(
                unknown, hint=" — call 'list-enzymes' for the catalog, or "
                               "'create-custom-enzyme' to add one")
        if not resolved:
            return ({"error":
                      "'enzymes' is empty — omit the key to scan the "
                      "catalog"}, 400)
    min_len = _coerce_int(payload.get("min_length", 4),
                            name="min_length")
    if isinstance(min_len, str):
        return ({"error": min_len}, 400)
    # A `min_length` longer than the longest recognition site in the catalog
    # can only ever return zero sites — and `count: 0` reads as "this
    # construct is clean", not as "that filter excluded every enzyme there
    # is". Same vacuous-pass shape as an unrecognised enzyme name, so it gets
    # the same refusal. The bound is read from the catalog so it can't go
    # stale when a longer-site enzyme is added.
    try:
        longest_site = max((len(str(v[0])) for v in
                            (_state._all_enzymes_hook() or {}).values()),
                           default=0)
    except Exception:                       # pragma: no cover - hook absent
        longest_site = 0
    if longest_site and min_len > longest_site:
        return ({"error":
                  f"min_length {min_len} is longer than the longest "
                  f"recognition site in the catalog ({longest_site} bp), so "
                  "this would report zero sites for every enzyme — check the "
                  "units (min_length is in BASE PAIRS)"}, 400)
    min_len = max(0, min_len)
    unique = bool(payload.get("unique_only", False))
    respect_active = bool(payload.get("respect_active_collection", True))
    seq = str(rec.seq)
    is_circular = rec.annotations.get("topology") == "circular"
    warnings: "list[str]" = []

    # A hand-picked name list goes INTO the scan as `allowed_enzymes` rather
    # than being post-filtered off a full-catalog scan. `_scan_restriction_
    # sites_impl` collapses isoschizomers and HF/v2 variants that land on the
    # same span (its `placed` set) so the map doesn't stack three bars on one
    # CGTCTC — which means a post-filter for "Esp3I" loses to the catalog-
    # order winner "BsmBI" and comes back empty. `allowed_enzymes` runs
    # BEFORE that collapse, and is also exactly what the UI overlay passes,
    # so this endpoint now genuinely mirrors it rather than claiming to.
    if resolved:
        allowed: "frozenset[str] | None" = frozenset(resolved)
        if "min_length" in payload:
            warnings.append(
                "'min_length' does not apply when 'enzymes' is given — a "
                "hand-picked list keeps short recognition sites")
    elif respect_active:
        allowed = _active_enzyme_allowed_set()
    else:
        allowed = None

    # `distinct_cuts=True`: this is a QUERY, so a NEOschizomer must not be
    # swallowed by the render's one-bar-per-span collapse. XmaI (`C^CCGGG`)
    # and SmaI (`CCC^GGG`) recognise the same six bases and leave DIFFERENT
    # ends, so answering `{"enzymes": ["XmaI", "SmaI"]}` with only whichever
    # the catalog lists first reads as "no XmaI site" on a plasmid that has
    # one. Same family as the isoschizomer loss [INV-187] fixed; the map keeps
    # collapsing (default False) because there one bar per span is the point.
    sites = _scan_restriction_sites(
        seq, min_recognition_len=min_len,
        unique_only=unique, circular=is_circular,
        allowed_enzymes=allowed, distinct_cuts=True,
    )

    # Re-expand the collapse for the names the caller actually asked about:
    # the scan emits ONE label per (span, recognition site), so a request
    # naming both BsmBI and Esp3I comes back labelled BsmBI only. Report the
    # hit under every requested name that cuts identically — the names
    # differ, the biology does not. `display` is what the CALLER typed (its
    # canonical spelling): a commercial synonym scans under its catalog
    # target, and reporting `BsaI` rows to someone who asked for `Eco31I`
    # would leave their `s["enzyme"] == "Eco31I"` filter empty — the same
    # silent zero one layer down.
    want_by_sig: "dict[tuple | None, list[str]]" = {}
    display_names: "list[str]" = []
    for raw in (enzymes or []):
        match = _enzyme_resolve_one(raw)
        if match is None:                   # already 400'd above
            continue
        catalog_name, display = match
        names = want_by_sig.setdefault(_enzyme_signature(catalog_name), [])
        if display not in names:
            names.append(display)
        if display not in display_names:
            display_names.append(display)

    out = []
    for s in sites:
        if s.get("type") != "resite":
            continue
        label = s.get("label", "")
        # Skip the unlabeled wrap continuation — sacred invariant #6 emits an
        # origin-spanning site as a labeled tail piece PLUS an unlabeled head
        # piece, and counting code must count only the labeled one. Pre-fix an
        # unfiltered scan emitted the head as a phantom `{"enzyme": ""}` row,
        # so `count` disagreed with itself depending on whether the caller
        # passed an enzyme filter.
        if not label:
            continue
        for nm in (want_by_sig.get(_enzyme_signature(label), [label])
                   if resolved else [label]):
            out.append({
                "enzyme":        nm,
                "start":         s.get("start"),
                "end":           s.get("end"),
                "strand":        s.get("strand", 1),
                "cut_bp":        s.get("top_cut_bp", -1),
                # Bottom-strand cut + recognition sequence ride along so a
                # script can check fragment-size arithmetic without keeping a
                # private table of cut offsets (blue field report #6).
                "bottom_cut_bp": s.get("bottom_cut_bp", -1),
                "site":          (_enzyme_signature(nm) or ("", 0, 0))[0],
                # An origin-spanning site has end <= start in plasmid
                # coordinates; say so rather than making the caller infer it.
                "wraps":         bool(s.get("rec_start") is not None
                                      and s.get("end") == len(seq)
                                      and s.get("rec_end", 0) > 0),
            })

    # `allowed_enzymes` deliberately overrides `unique_only` inside the
    # scanner ("unique cutters of MY hand-picked list" hides the multi-cutter
    # the user just typed in), so re-apply it here for callers that asked.
    if unique and allowed is not None:
        counts: "dict[str, int]" = {}
        for r in out:
            counts[r["enzyme"]] = counts.get(r["enzyme"], 0) + 1
        out = [r for r in out if counts[r["enzyme"]] == 1]

    resp: dict = {"sites": out, "count": len(out)}
    if resolved:
        # Report the spellings the ROWS are labelled with — for a commercial
        # synonym that is the caller's own name, and the catalog's name for it
        # then shows up under `equivalent_enzymes` (so `Eco31I` → scanned
        # Eco31I, equivalent BsaI/BspTNI: the answer says what it looked for
        # AND what that enzyme is otherwise called).
        resp["enzymes_scanned"] = display_names or resolved
        # Names that cut identically to a requested one. A zero-site answer is
        # only trustworthy if the caller knows which spellings it covered.
        eq = sorted({a for nm in resolved for a in _enzyme_aliases(nm)}
                    - set(display_names or resolved))
        if eq:
            resp["equivalent_enzymes"] = eq
    if warnings:
        resp["warnings"] = warnings
    return resp


def _agent_enzyme_dict(name: str, site: str, fwd_cut: int,
                       rev_cut: int, custom_names: "set[str]") -> dict:
    """One catalog row for `list-enzymes`.

    `fwd_cut` / `rev_cut` are the catalog's own offsets from the START of the
    recognition site to the top- and bottom-strand cut. Everything else is
    derived from them with the SAME rule `_enzyme_cuts_impl` uses (blunt when
    the two coincide, 5' when the top cut comes first, 3' when it doesn't), so
    a script checking fragment arithmetic against this row can't disagree with
    the digest engine that produced the fragments."""
    site_u   = str(site).upper()
    site_len = len(site_u)
    overhang = rev_cut - fwd_cut
    kind = "blunt" if overhang == 0 else ("5'" if overhang > 0 else "3'")
    # Human notation: a within-site cutter carries a caret at the top-strand
    # cut (SalI G^TCGAC); one that reaches outside gets REBASE's offset form
    # (BsaI GGTCTC(1/5), BaeI (10/15)ACNNNNGTAYC for an upstream cutter).
    if 0 < fwd_cut < site_len:
        notation = site_u[:fwd_cut] + "^" + site_u[fwd_cut:]
    elif fwd_cut < 0:
        notation = "({}/{}){}".format(-fwd_cut, -rev_cut, site_u)
    else:
        notation = "{}({}/{})".format(
            site_u, fwd_cut - site_len, rev_cut - site_len)
    return {
        "name":            name,
        "site":            site_u,
        "recognition_len": site_len,
        "fwd_cut":         fwd_cut,
        "rev_cut":         rev_cut,
        "notation":        notation,
        "overhang_len":    abs(overhang),
        "overhang_kind":   kind,
        "type_iis":        _enzyme_is_type_iis(name),
        "custom":          name in custom_names,
        # Catalog names that cut identically. The scan collapses these onto
        # one label, so a caller needs to know which spellings a result
        # already covers.
        "aliases":         [a for a in _enzyme_aliases(name) if a != name],
    }


def _agent_transcript_feature_dicts(rec) -> "list[dict]":
    """`_agent_carry_feature_dicts`' shape PLUS the location parts.

    A spliced CDS or mRNA is a compound location, and the gaps between its
    parts ARE the introns — the standard GenBank way of saying so. The
    carry-annotations shape drops them (a digest only needs outer bounds),
    and dropping them here would silently turn a spliced transcript into an
    unspliced one, which is the exact failure `predict-transcript` exists to
    prevent.

    An ORIGIN-WRAPPING feature is ALSO a two-part compound location (sacred
    invariant #9), so `parts` is emitted only when `_feat_bounds` did not
    already resolve the location as a wrap — otherwise every feature crossing
    bp 0 would be read as carrying an intron it doesn't have."""
    out: "list[dict]" = []
    total = _seq_len(rec)
    # Topology matters here: on a LINEAR record the two-part shape that looks
    # like an origin wrap is a spliced feature, and resolving it as a wrap
    # suppressed the `parts` emit below — silently turning a spliced
    # transcript into an unspliced one.
    _circular = _record_is_circular(rec)
    for f in getattr(rec, "features", None) or []:
        if getattr(f, "type", "") == "source":
            continue
        bounds = _feat_bounds(f, total, circular=_circular)
        if bounds is None:
            continue
        s, e, strand = bounds
        quals = getattr(f, "qualifiers", {}) or {}
        label = (quals.get("label") or quals.get("product") or [f.type])[0]
        # 0 is a strand, not a falsy default — see `_agent_carry_feature_dicts`.
        d = {"start": int(s), "end": int(e),
             "strand": 0 if strand is None else int(strand),
             "type": f.type, "label": str(label)[:200]}
        if e >= s:                       # not a wrap — parts may be exons
            parts: "list[list[int]]" = []
            for part in (getattr(getattr(f, "location", None), "parts", None)
                         or []):
                try:
                    parts.append([int(part.start), int(part.end)])
                except (TypeError, ValueError, AttributeError):
                    parts = []
                    break
            if len(parts) > 1:
                d["parts"] = parts
        out.append(d)
    return out


# Caps shared by every endpoint that takes a list of enzyme names. 500
# matches `digest`'s existing limit; the echo cap keeps a validation error
# from becoming a megabyte response when the payload carries a huge string.
_AGENT_MAX_ENZYME_NAMES = 500
_AGENT_NAME_ECHO_MAX = 64
# Longest an enzyme NAME may be. The catalog's longest is 10 characters and
# `create-custom-enzyme` caps at 40; this is the shape gate that stops a
# multi-megabyte string being lowercased and hashed on its way to "unknown".
_AGENT_MAX_ENZYME_NAME_LEN = 64


def _agent_unknown_enzyme_error(unknown, *, field: str = "enzymes",
                                hint: str = ""):
    """The 400 body for unrecognised enzyme name(s), in ONE place.

    `unknown` is `_resolve_enzyme_names`' second return value. Names are
    CLAMPED before they reach the message and only the first few are quoted:
    a caller can send thousands of bad names, each arbitrarily long, and
    echoing them verbatim turns a validation error into a denial of service
    against whoever reads the response. `n_unknown` reports the real total so
    the truncation can't be mistaken for the whole story."""
    shown = unknown[:10]
    bits = []
    for want, near in shown:
        clipped = want[:_AGENT_NAME_ECHO_MAX]
        if len(want) > _AGENT_NAME_ECHO_MAX:
            clipped += "…"
        bits.append(repr(clipped) + (
            " (did you mean " + ", ".join(near) + "?)" if near else ""))
    more = len(unknown) - len(shown)
    body = {
        "error": (f"unknown enzyme(s) in '{field}': " + "; ".join(bits)
                  + (f" (+{more} more)" if more else "") + hint),
        "unknown_enzymes": [w[:_AGENT_NAME_ECHO_MAX] for w, _ in shown],
        "n_unknown": len(unknown),
    }
    return (body, 400)


def _agent_check_enzyme_list(names, *, field: str = "enzymes"):
    """Shape + size check for an enzyme-name list. Returns an error tuple or
    None. Bounds the list BEFORE resolution so the per-name near-miss search
    can't be turned into a CPU sink."""
    if not isinstance(names, list):
        return ({"error": f"'{field}' must be a list (or a single string)"},
                400)
    if len(names) > _AGENT_MAX_ENZYME_NAMES:
        return ({"error": f"too many enzyme names in '{field}' "
                          f"(max {_AGENT_MAX_ENZYME_NAMES})"}, 400)
    if not all(isinstance(e, str) for e in names):
        return ({"error": f"'{field}' must contain only strings"}, 400)
    # No enzyme name is anywhere near this long (the longest in the catalog is
    # 10 characters, and a custom one is capped at 40), so a longer string is
    # a payload error — refuse it on SHAPE rather than lowercasing and hashing
    # a megabyte on the way to "unknown enzyme".
    if any(len(e) > _AGENT_MAX_ENZYME_NAME_LEN for e in names):
        return ({"error": f"an enzyme name in '{field}' exceeds "
                          f"{_AGENT_MAX_ENZYME_NAME_LEN} characters"}, 400)
    return None


def _agent_resolve_scan_target(app, payload, *, want_record: bool = False):
    """Resolve the molecule a regulatory scan runs over.

    Returns ``(seq, circular, record_or_None, label, None)`` or
    ``("", True, None, "", (error_dict, status))``.

    As with `_agent_resolve_read_reference`, the error path carries sentinels
    rather than ``None`` so the success values keep concrete types; the fifth
    element is the one every caller checks.

    Accepts a bare ``sequence``, a saved plasmid by ``id`` / ``name``, or —
    when neither is given — the currently loaded record. `want_record` makes a
    bare sequence an error, because a scan that needs FEATURES cannot get them
    from loose bases and silently scanning with an empty feature list would
    report every gene as absent rather than saying the input was wrong.
    """
    raw = payload.get("sequence")
    if isinstance(raw, str) and raw.strip():
        if want_record:
            return ("", True, None, "",
                    ({"error": "this endpoint needs annotations — pass 'id' "
                               "or 'name' for a saved plasmid, or load one, "
                               "rather than a bare 'sequence'"}, 400))
        seq, err = _sanitize_bases(raw)
        if err:
            return ("", True, None, "",
                    ({"error": f"'sequence' rejected: {err}"}, 400))
        circ = payload.get("circular")
        return (seq, (True if circ is None else bool(circ)), None,
                "sequence", None)
    key = _sanitize_label(payload.get("id") or payload.get("name"),
                          max_len=200)
    if key:
        entry = (_find_library_entry_by_id(key)
                 or _find_library_entry_by_name(key))
        if entry is None:
            return ("", True, None, "",
                    ({"error": f"no library entry named {key!r}"}, 404))
        try:
            rec = _gb_text_to_record(entry.get("gb_text") or "")
        except Exception as exc:
            return ("", True, None, "",
                    ({"error": f"couldn't parse {key!r}: "
                               f"{_scrub_path(str(exc))}"}, 500))
        label = entry.get("name") or key
    else:
        rec = getattr(app, "_current_record", None)
        if rec is None:
            return ("", True, None, "",
                    ({"error": "no plasmid loaded and no 'id' / 'name' / "
                               "'sequence' given"}, 422))
        # Same resolution order as `PlasmidApp._record_display_name`, which
        # is hub-side and so out of reach here: the TYPED display name first,
        # LOCUS only as a fallback. Never `record.name` alone — the LOCUS is
        # a sanitised, underscore-bearing 28-column field, not the name the
        # user gave the plasmid ([INV-98]).
        label = (_sanitize_label(getattr(rec, "_tui_display_name", None),
                                 max_len=200)
                 or _sanitize_label(getattr(rec, "name", None), max_len=200)
                 or _sanitize_label(getattr(rec, "id", None), max_len=200)
                 or "?")
    seq = str(getattr(rec, "seq", "") or "")
    if not seq:
        return ("", True, None, "",
                ({"error": f"{label!r} has no sequence"}, 422))
    if len(seq) > _PAIRWISE_MAX_LEN:
        return ("", True, None, "",
                ({"error": f"{label!r} exceeds {_PAIRWISE_MAX_LEN:,} bp"},
                 413))
    circ_override = payload.get("circular")
    circular = (_record_is_circular(rec) if circ_override is None
                else bool(circ_override))
    return (seq, circular, rec, label, None)


@_agent_endpoint("scan-promoters")
def _h_scan_promoters(app, payload):
    """Scan for sigma-70 (-35 / spacer / -10) promoters.

    Body: ``{sequence | id | name, circular?, both_strands?=true,
    min_score?=0.75, limit?}``. Coordinates are FORWARD-strand and 0-based for
    hits on either strand; on a circular molecule a promoter spanning the
    origin is found at its real start.

    Each hit: ``{start, end, strand, minus35, minus35_start, spacer, minus10,
    minus10_start, tss, score, strength, components, model}``.

    **`score` ranks, it does not predict.** It is a composite of named
    components — how well each hexamer matches the consensus weighted by which
    positions are conserved, the spacer geometry, and an extended -10 bonus —
    and `components` breaks it out so a caller can show its working. There is
    no trained model behind it and it is not a transcription rate. What it IS
    good for is the question that gets asked: "does this insert carry a
    cryptic promoter, and how does it compare with the others?"

    `min_score` is calibrated rather than arbitrary: at the 0.75 default,
    random DNA yields about one hit per kb per strand, while real promoters
    from the built-in preset catalogue score 0.69-1.00 (median 0.77). Raise it
    to 0.90 for unambiguous promoters only (~1 chance hit per 50 kb); lower it
    to hunt weak ones and expect the noise. Read-only."""
    seq, circular, _rec, label, err = _agent_resolve_scan_target(app, payload)
    if err is not None:
        return err
    both = payload.get("both_strands")
    if both is None:
        both = True
    else:
        both = _coerce_bool(both, name="both_strands")
        if isinstance(both, str):
            return ({"error": both}, 400)
    min_score = payload.get("min_score")
    if min_score is not None:
        try:
            min_score = float(min_score)
        except (TypeError, ValueError):
            return ({"error": "'min_score' must be a number"}, 400)
        if not (0.0 <= min_score <= 1.0):
            return ({"error": "'min_score' must be in [0, 1]"}, 400)
    limit = _coerce_int(payload.get("limit", 500), name="limit")
    if isinstance(limit, str):
        return ({"error": limit}, 400)
    if limit < 1:
        return ({"error": "'limit' must be >= 1"}, 400)
    try:
        hits = _scan_promoters(seq, circular=circular, both_strands=both,
                               min_score=min_score)
    except Exception as exc:
        _log.exception("scan-promoters failed")
        return ({"error": f"scan failed: {_scrub_path(str(exc))}"}, 500)
    _log_event("scan.promoters", length=len(seq), n_hits=len(hits),
               via="agent")
    return {"ok": True, "plasmid": label, "length": len(seq),
            "circular": circular, "both_strands": both,
            "min_score": (_PROMOTER_MIN_SCORE if min_score is None
                          else min_score),
            "n_hits": len(hits), "truncated": len(hits) > limit,
            "hits": hits[:limit],
            "ignored": _agent_ignored_keys(payload, {
                "sequence", "id", "name", "circular", "both_strands",
                "min_score", "limit"})}


@_agent_endpoint("scan-terminators")
def _h_scan_terminators(app, payload):
    """Scan for intrinsic (Rho-independent) terminators: a stem-loop followed
    by a U-tract.

    Body: ``{sequence | id | name, circular?, both_strands?=true, min_stem?=5,
    min_u_tract?=4, limit?}``. Coordinates are FORWARD-strand and 0-based for
    hits on either strand.

    Each hit: ``{start, end, strand, hairpin_start, hairpin_end, stem_len,
    loop_len, u_tract, u_tract_start, dg, efficiency}``.

    **ONE RECORD PER HAIRPIN.** Overlapping stem/loop registers of the same
    structure are merged to the best one — a scan that reports every register
    separately turns a single hairpin into twenty-odd "terminators", which is
    worse than not scanning.

    ``dg`` is a real Turner-2004 fold of the hairpin in kcal/mol (``null`` if
    it could not be folded), not a heuristic. ``efficiency`` is a coarse named
    tier — ``strong`` / ``moderate`` / ``weak`` — calibrated so that `strong`
    covers every real intrinsic terminator in the built-in preset catalogue
    (rrnB T1, T7 Tphi, lambda t0, T7Te, BBa_B1008: dg -12.1 to -31.3) while
    random DNA reaches it about once per 25 kb. Never a percentage: how well a
    terminator works in vivo depends on the polymerase and the strain.

    Rho-DEPENDENT termination is not modelled and is not attempted — it has no
    hairpin signature to scan for. Read-only."""
    seq, circular, _rec, label, err = _agent_resolve_scan_target(app, payload)
    if err is not None:
        return err
    both = payload.get("both_strands")
    if both is None:
        both = True
    else:
        both = _coerce_bool(both, name="both_strands")
        if isinstance(both, str):
            return ({"error": both}, 400)
    min_stem = _coerce_int(payload.get("min_stem", _TERM_MIN_STEM),
                           name="min_stem")
    if isinstance(min_stem, str):
        return ({"error": min_stem}, 400)
    # The ceilings are the scanner's own search bounds, not round numbers. A
    # `min_stem` above `_TERM_MAX_STEM` (or a `min_u_tract` past the window the
    # tract is measured in) can only ever return zero hits — and "0 terminators"
    # reads as "this construct has no terminator", not as "that filter excluded
    # every candidate there is". Same vacuous-pass refusal as an out-of-range
    # `min_length` on `list-restriction-sites`.
    if not (2 <= min_stem <= _TERM_MAX_STEM):
        return ({"error": f"'min_stem' must be 2-{_TERM_MAX_STEM} — the "
                          f"scanner searches stems up to {_TERM_MAX_STEM} bp, "
                          f"so a higher floor would report zero terminators "
                          f"for every sequence"}, 400)
    min_u = _coerce_int(payload.get("min_u_tract", _TERM_MIN_U_TRACT),
                        name="min_u_tract")
    if isinstance(min_u, str):
        return ({"error": min_u}, 400)
    if not (1 <= min_u <= _TERM_U_TRACT_WINDOW):
        return ({"error": f"'min_u_tract' must be 1-{_TERM_U_TRACT_WINDOW} — "
                          f"the U-tract is measured over a "
                          f"{_TERM_U_TRACT_WINDOW}-base window, so a higher "
                          f"floor would report zero terminators for every "
                          f"sequence"}, 400)
    limit = _coerce_int(payload.get("limit", 500), name="limit")
    if isinstance(limit, str):
        return ({"error": limit}, 400)
    if limit < 1:
        return ({"error": "'limit' must be >= 1"}, 400)
    try:
        hits = _scan_terminators(seq, circular=circular, both_strands=both,
                                 min_stem=min_stem, min_u_tract=min_u)
    except Exception as exc:
        _log.exception("scan-terminators failed")
        return ({"error": f"scan failed: {_scrub_path(str(exc))}"}, 500)
    _log_event("scan.terminators", length=len(seq), n_hits=len(hits),
               via="agent")
    return {"ok": True, "plasmid": label, "length": len(seq),
            "circular": circular, "both_strands": both,
            "min_stem": min_stem, "min_u_tract": min_u,
            "n_hits": len(hits), "truncated": len(hits) > limit,
            "hits": hits[:limit],
            "ignored": _agent_ignored_keys(payload, {
                "sequence", "id", "name", "circular", "both_strands",
                "min_stem", "min_u_tract", "limit"})}


@_agent_endpoint("map-transcription")
def _h_map_transcription(app, payload):
    """Walk the circle from every promoter to every CDS — "what ELSE
    transcribes this gene?"

    Body: ``{id? | name?, circular?, include_predicted?=true,
    min_promoter_score?, max_distance?}``. Omit `id`/`name` to use the loaded
    plasmid. Needs an annotated record (CDS features), so a bare sequence is
    refused rather than silently answered with no genes.

    **Why this is not `predict-transcript`.** That endpoint reconstructs the
    mature message of ONE unit, and to do it the caller must already name the
    promoter driving the CDS. It therefore cannot answer the question that
    costs bench time: a plasmid is a circle, every promoter on it transcribes
    until something stops it, and on a construct with no designed insulators
    that can be the whole way round. Checking terminators DOWNSTREAM of a
    cassette — the natural thing to do — cannot see a read-through, because
    the read-through arrives from behind.

    Returns ``{ok, length, circular, n_promoters, n_terminators, n_cds,
    promoters, terminators, genes, warnings}``. Per gene: ``reachable_from``
    (each ``{promoter, distance_bp, direction, terminators_between,
    readthrough}``), ``primary``, ``silent``, plus ``host_driven_by`` /
    ``silent_in_host`` / ``host_primary``.

    ``readthrough`` is ``clear`` / ``attenuated`` / ``blocked`` — a coarse
    named call over the terminators in the way, never a percentage.

    **`silent` and `silent_in_host` are different claims and the gap between
    them is the point.** `silent` means nothing at all is aimed at the gene.
    `silent_in_host` means nothing the HOST polymerase can read reaches it —
    so a T7-only cassette is `silent_in_host` in a strain with no T7 RNAP,
    UNLESS something else on the circle reads through into it, which is
    exactly the case that looks like a mystery at the bench. Each promoter
    carries ``host_recognised`` with the ``host_basis`` that decided it.

    With `include_predicted` (default true), scanned sigma-70 promoters and
    intrinsic terminators are folded in beside the annotated ones, each tagged
    ``source``. Read-only, heavy."""
    seq, circular, rec, label, err = _agent_resolve_scan_target(
        app, payload, want_record=True)
    if err is not None:
        return err
    inc = payload.get("include_predicted")
    if inc is None:
        inc = True
    else:
        inc = _coerce_bool(inc, name="include_predicted")
        if isinstance(inc, str):
            return ({"error": inc}, 400)
    min_score = payload.get("min_promoter_score")
    if min_score is not None:
        try:
            min_score = float(min_score)
        except (TypeError, ValueError):
            return ({"error": "'min_promoter_score' must be a number"}, 400)
        if not (0.0 <= min_score <= 1.0):
            return ({"error": "'min_promoter_score' must be in [0, 1]"}, 400)
    max_distance = payload.get("max_distance")
    if max_distance is not None:
        max_distance = _coerce_int(max_distance, name="max_distance")
        if isinstance(max_distance, str):
            return ({"error": max_distance}, 400)
        if max_distance < 1:
            return ({"error": "'max_distance' must be >= 1"}, 400)
    feats = _agent_transcript_feature_dicts(rec)
    try:
        out = _map_transcription(seq, feats, circular=circular,
                                 include_predicted=inc,
                                 min_promoter_score=min_score,
                                 max_distance=max_distance)
    except Exception as exc:
        _log.exception("map-transcription failed")
        return ({"error": f"map failed: {_scrub_path(str(exc))}"}, 500)
    out["plasmid"] = label
    out["ignored"] = _agent_ignored_keys(payload, {
        "id", "name", "circular", "include_predicted", "min_promoter_score",
        "max_distance", "host"})
    return out


@_agent_endpoint("predict-transcript")
def _h_predict_transcript(app, payload):
    """Reconstruct the MATURE mRNA of an annotated transcription unit on the
    loaded record, and score translation initiation on it.

    Body: ``{promoter?, cds?, terminator?, tx_start?, tx_end?,
    check_splice?=true, min_uorf_aa?=1, include_sequences?=true}``. Set
    ``include_sequences: false`` to drop `pre_mrna` / `mature_mrna` / the UTR
    sequences and keep only the coordinates and verdicts — a long unit
    otherwise puts three copies of itself in the response, which a screening
    loop over a whole collection doesn't need. Each of `promoter` / `cds` /
    `terminator` takes a feature index (as `list-features` numbers them) or a
    label; omit one to pick by feature type — the longest CDS, and the nearest
    promoter/terminator on the same strand, measured AROUND the circle rather
    than by raw coordinate order. `tx_start` / `tx_end` bound the unit
    directly (plus-strand coordinates) when the record annotates no promoter
    or terminator.

    **Why this exists.** When a construct carries a genomic 5'UTR intron, the
    biology depends on the SPLICED leader while every check written against
    the plasmid sequence sees the unspliced one — so an upstream-ATG screen
    reports intronic ATGs that splicing removes. They look exactly like real
    uORF hazards. The response separates the two: ``uorfs`` are the ones a
    ribosome actually meets on the mature message, and
    ``removed_by_splicing`` names each ATG that is in the DNA and is NOT in
    the message, with the intron that removes it.

    Returns ``{ok, unit, pre_mrna, mature_mrna, spliced, exons, introns,
    five_utr, cds, three_utr, kozak, uorfs, removed_by_splicing, splice,
    warnings}``. Transcript coordinates are offsets into `mature_mrna`
    (`removed_by_splicing` is in `pre_mrna` offsets and also carries the
    genomic position). The unit is taken to run from the END of the promoter
    (a promoter is bound, not transcribed) to the end of the terminator, and
    it is rotation- and strand-invariant — the same construct answers
    identically whichever way it is drawn.

    **Introns come from ANNOTATION only** — `intron` features, or the gaps
    between a spliced location's parts. Nothing here predicts where a
    spliceosome would cut; `splice` carries a separate ADVISORY scan of the
    pre-mRNA with the calibrated plant PWM, marking sites that coincide with
    an annotated boundary so only genuinely cryptic ones are reported. If the
    model is unavailable the scan reports `skipped` with a reason — never
    an empty list that would read as clean."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    for key in ("promoter", "cds", "terminator"):
        v = payload.get(key)
        if v is not None and not isinstance(v, (int, str)) or isinstance(v, bool):
            return ({"error": f"'{key}' must be a feature index or a label"},
                    400)
    bounds: "dict[str, int]" = {}
    for key in ("tx_start", "tx_end"):
        if payload.get(key) is None:
            continue
        v = _coerce_int(payload.get(key), name=key)
        if isinstance(v, str):
            return ({"error": v}, 400)
        bounds[key] = int(v)
    if len(bounds) == 1:
        return ({"error": "pass BOTH 'tx_start' and 'tx_end', or neither"},
                400)
    min_uorf_aa = _coerce_int(payload.get("min_uorf_aa", 1),
                              name="min_uorf_aa")
    if isinstance(min_uorf_aa, str):
        return ({"error": min_uorf_aa}, 400)
    try:
        result = _predict_transcript(
            str(rec.seq),
            _agent_transcript_feature_dicts(rec),
            circular=(rec.annotations.get("topology") == "circular"),
            promoter=payload.get("promoter"),
            cds=payload.get("cds"),
            terminator=payload.get("terminator"),
            tx_start=bounds.get("tx_start"),
            tx_end=bounds.get("tx_end"),
            check_splice=bool(payload.get("check_splice", True)),
            min_uorf_aa=max(0, int(min_uorf_aa)),
        )
    except Exception as exc:
        _log.exception("agent predict-transcript failed")
        return ({"error":
                  f"transcript prediction failed: {_scrub_path(str(exc))}"},
                500)
    if isinstance(result, dict) and not bool(
            payload.get("include_sequences", True)):
        # The message itself is the point of this endpoint, so it ships by
        # default — but a screening loop over a whole collection only wants
        # the verdicts, and a long transcription unit puts three copies of
        # itself in the response. Mirrors `digest`'s `include_fragment_seq`.
        for key in ("pre_mrna", "mature_mrna"):
            result.pop(key, None)
        for block in ("five_utr", "three_utr"):
            if isinstance(result.get(block), dict):
                result[block].pop("seq", None)
    if not result.get("ok"):
        # The unit couldn't be read — 422, with the reason. Never a partial
        # transcript: half a message answered confidently is the failure mode
        # this endpoint was written to remove.
        return ({"error": result.get("error", "could not read a transcription "
                                              "unit from this record"),
                  "warnings": result.get("warnings", [])}, 422)
    return result


@_agent_endpoint("list-enzymes")
def _h_list_enzymes(app, payload):
    """Look up the restriction-enzyme catalog: recognition site, cut offsets,
    overhang and isoschizomers, with no sequence involved.

    Body: ``{names?: [str, ...] | str, search?: str, type_iis_only?: bool,
    custom_only?: bool, limit?: int = 500}``.

    Covers the COMBINED catalog — built-in NEB plus the user's custom
    enzymes (``list-custom-enzymes`` returns only the latter). Exists so a
    script can verify a digest's fragment sizes against the enzyme's real
    offsets instead of hardcoding a private table of them (blue field report
    #6), and so `list-restriction-sites`'s "unknown enzyme" 400 has somewhere
    to point.

    * ``names`` — exact rows, matched case-insensitively. An unrecognised
      name is a 400 with near-miss suggestions, matching
      `list-restriction-sites`; use ``search`` for open-ended lookup.
    * ``search`` — case-insensitive substring over the name AND the
      recognition site, so ``"GGTCTC"`` finds BsaI and its isoschizomers.

    Each row is ``{name, site, recognition_len, fwd_cut, rev_cut, notation,
    overhang_len, overhang_kind, type_iis, custom, aliases}``. `fwd_cut` /
    `rev_cut` are offsets from the start of the recognition site to the top-
    and bottom-strand cut; `notation` renders them the way REBASE does
    (``G^TCGAC``, ``GGTCTC(1/5)``)."""
    catalog = _state._all_enzymes_hook() or {}
    names = payload.get("names")
    if isinstance(names, str):
        names = [names]
    if names is not None:
        if (shape_err := _agent_check_enzyme_list(
                names, field="names")) is not None:
            return shape_err
        resolved, unknown = _resolve_enzyme_names(names)
        if unknown:
            return _agent_unknown_enzyme_error(unknown, field="names")
        wanted = resolved
    else:
        wanted = sorted(catalog)

    search = payload.get("search")
    if search is not None and not isinstance(search, str):
        return ({"error": "'search' must be a string"}, 400)
    limit = _coerce_int(payload.get("limit", 500), name="limit")
    if isinstance(limit, str):
        return ({"error": limit}, 400)
    limit = max(1, min(int(limit), 2000))
    type_iis_only = bool(payload.get("type_iis_only", False))
    custom_only   = bool(payload.get("custom_only", False))
    custom_names = {str(e.get("name")) for e in _load_custom_enzymes()
                    if isinstance(e, dict) and e.get("name")}

    rows = []
    needle = (search or "").strip().upper()
    for name in wanted:
        entry = catalog.get(name)
        if not entry:
            continue
        try:
            site, fwd_cut, rev_cut = entry
            row = _agent_enzyme_dict(name, str(site), int(fwd_cut),
                                     int(rev_cut), custom_names)
        except (TypeError, ValueError):
            # A corrupt custom enzyme — same tolerance as the scan-catalog
            # rebuild, which skips it rather than failing the whole read.
            continue
        if needle and needle not in name.upper() and needle not in row["site"]:
            continue
        if type_iis_only and not row["type_iis"]:
            continue
        if custom_only and not row["custom"]:
            continue
        rows.append(row)

    truncated = len(rows) > limit
    return {"enzymes": rows[:limit], "count": min(len(rows), limit),
            "total": len(rows), "truncated": truncated}


def _agent_parse_span(spec, label: str):
    """Coerce one span spec into ``(start, end)`` or return an error string.

    Accepts ``{start, end}``, ``{bp}`` (a single base — the half-open span
    ``[bp, bp+1)``), or a two-element ``[start, end]`` list."""
    if isinstance(spec, (list, tuple)):
        if len(spec) != 2:
            return None, f"'{label}' as a list must be [start, end]"
        spec = {"start": spec[0], "end": spec[1]}
    if not isinstance(spec, dict):
        return None, f"'{label}' must be an object with start/end (or bp)"
    if "bp" in spec and "start" not in spec:
        bp = _coerce_int(spec.get("bp"), name=f"{label}.bp")
        if isinstance(bp, str):
            return None, bp
        return (int(bp), int(bp) + 1), None
    if "start" not in spec or "end" not in spec:
        return None, f"'{label}' needs both 'start' and 'end' (or 'bp')"
    start = _coerce_int(spec.get("start"), name=f"{label}.start")
    if isinstance(start, str):
        return None, start
    end = _coerce_int(spec.get("end"), name=f"{label}.end")
    if isinstance(end, str):
        return None, end
    return (int(start), int(end)), None


@_agent_endpoint("span-contains")
def _h_span_contains(app, payload):
    """Circular-aware containment: does a span contain a base or another
    span, on a molecule where either may cross the origin?

    Body: ``{outer: {start, end}, inner: {start, end} | {bp} | [...],
    length?: int}``.

    `length` defaults to the loaded record's; pass it to work on
    coordinates from somewhere else. Spans are half-open ``[start, end)``,
    the same convention `list-features` and `list-restriction-sites` report,
    and a span with ``end <= start`` crosses the origin.

    `inner` may be one span or a list of them, and the response is always
    the same shape: ``{results: [{start, end, length, wraps, contains}],
    all_contained, outer: {start, end, length, wraps}, length}``. One span
    in, one row out — there is no separate single-span response to branch on.

    This exists because ``(x - a) % L < (b - a) % L`` is the kind of
    arithmetic that looks right and isn't: the obvious linear form
    ``a <= x < b`` reports every base of an origin-spanning T-DNA, marker or
    operon as OUTSIDE it, and a construct that is actually correct fails the
    check (blue field report #2). It shares `_feat_len` with the rest of the
    app, so containment and reported length cannot disagree."""
    length = payload.get("length")
    if length is None:
        rec = getattr(app, "_current_record", None)
        if rec is None:
            return ({"error":
                      "no plasmid loaded — pass 'length' to test coordinates "
                      "without one"}, 422)
        total = _seq_len(rec)
    else:
        coerced = _coerce_int(length, name="length")
        if isinstance(coerced, str):
            return ({"error": coerced}, 400)
        total = int(coerced)
    if total <= 0:
        return ({"error": "'length' must be > 0"}, 400)

    if "outer" not in payload:
        return ({"error": "missing 'outer'"}, 400)
    outer, err = _agent_parse_span(payload.get("outer"), "outer")
    if err is not None:
        return ({"error": err}, 400)
    assert outer is not None
    o_start, o_end = outer[0] % total, outer[1] % total

    if "inner" not in payload:
        return ({"error": "missing 'inner'"}, 400)
    raw_inner = payload.get("inner")
    # A bare [start, end] pair is one span, not a list of two specs — only a
    # list whose members are themselves dicts/lists is a batch.
    many = (isinstance(raw_inner, list)
            and all(isinstance(x, (dict, list, tuple)) for x in raw_inner))
    specs = list(raw_inner) if many else [raw_inner]
    if len(specs) > 1000:
        return ({"error": "too many 'inner' spans (max 1000)"}, 400)
    # An empty list would answer `all_contained: true` having checked
    # nothing — "everything fits" is exactly the wrong reading of "I was
    # given nothing to check".
    if not specs:
        return ({"error": "'inner' is empty — nothing to test"}, 400)

    # Coordinates are circular, so a negative or past-the-end value is
    # normalised rather than refused (`-1` meaning the last base is genuinely
    # useful). But a value far outside the molecule is far more often a units
    # or off-by-a-plasmid mistake than an intentional lap, and the answer that
    # comes back looks equally confident either way — so say when it happened.
    out_of_range = sum(1 for v in (outer[0], outer[1]) if not 0 <= v < total)
    results = []
    for k, spec in enumerate(specs):
        span, err = _agent_parse_span(spec, f"inner[{k}]" if many else "inner")
        if err is not None:
            return ({"error": err}, 400)
        assert span is not None
        out_of_range += sum(1 for v in span if not 0 <= v < total)
        i_start, i_end = span[0] % total, span[1] % total
        results.append({
            "start":    i_start,
            "end":      i_end,
            "length":   _feat_len(i_start, i_end, total),
            "wraps":    bool(i_end < i_start),
            "contains": _span_in_span(i_start, i_end, o_start, o_end, total),
        })

    resp = {
        "results":       results,
        "all_contained": all(r["contains"] for r in results),
        "outer": {
            "start":  o_start,
            "end":    o_end,
            "length": _feat_len(o_start, o_end, total),
            "wraps":  bool(o_end < o_start),
        },
        "length": total,
    }
    if out_of_range:
        resp["warnings"] = [
            f"{out_of_range} coordinate(s) were outside [0, {total}) and were "
            "wrapped around the molecule — check the units if that wasn't "
            "intended"]
    return resp


@_agent_endpoint("digest")
def _h_digest(app, payload):
    """Digest a RAW sequence with restriction enzymes and report the cuts
    + resulting fragments with their overhangs — the overhang-aware QC a
    Golden-Braid / restriction build needs, on ANY sequence, without a
    sandboxed `import splicecraft` (agent-API feedback P1-4).

    Operates on the sequence in the BODY, not the loaded canvas, so an
    agent can verify a junction or a dropout overhang it just designed.

    Body: ``{sequence|seq, enzymes[]|enzyme, circular?=true,
             include_fragment_seq?=false}``.

    * ``sequence`` — IUPAC DNA, capped at 5 Mb (``seq`` accepted as alias).
    * ``enzymes`` — list of names from the combined catalog (built-in NEB ∪
      your custom enzymes); a singular ``enzyme`` string is accepted. Names
      the catalog doesn't know are reported under ``unknown_enzymes`` (and
      contribute no cuts) rather than silently dropped.
    * ``circular`` (default ``true``) — origin-spanning cuts + fragment
      count (a circular molecule with N cuts → N fragments; linear → N+1).
    * ``include_fragment_seq`` (default ``false``) — include each fragment's
      top-strand sequence; off by default so a megabase digest doesn't echo
      its whole self back.

    Returns:
      * ``cuts`` — ``[{enzyme, top, bottom, kind, overhang_seq}]`` in cut
        order; ``top``/``bottom`` are absolute 0-based top/bottom-strand cut
        coordinates, ``kind`` ∈ ``5'`` / ``3'`` / ``blunt``.
      * ``fragments`` — ``[{index, length, left, right}]`` in cut order
        around the molecule (5'→3' for linear); ``left``/``right`` are
        ``{overhang_seq, kind, enzyme}`` (the empty-enzyme ``linear`` kind
        marks a free end of a linear input). ``seq`` is added per fragment
        only when ``include_fragment_seq`` is set.
      * ``n_cuts``, ``n_fragments``, ``circular``, ``unknown_enzymes``.
      * ``resolved_enzymes`` — present only when a requested name is a
        commercial synonym (Thermo's ``Eco31I`` is NEB's ``BsaI``, ``LguI``
        is ``SapI``, ``AarI`` is ``PaqCI``): ``{asked: catalog_name}``, so
        the name in ``cuts[].enzyme`` can be traced back to the one you sent.
    """
    raw = (payload.get("sequence")
           if payload.get("sequence") not in (None, "")
           else payload.get("seq"))
    # `_sanitize_bases` enforces the cap (and IUPAC charset) itself — pass
    # the 5 Mb digest ceiling so it bounds the work (an HTTP caller is also
    # limited by the ~1 MiB JSON body cap upstream).
    bases, berr = _sanitize_bases(raw, max_len=_PCR_MAX_TEMPLATE_BP)
    if berr is not None:
        # `_sanitize_bases` speaks "'bases'"; this endpoint's field is
        # 'sequence' — translate so the error names a key the agent sent.
        return ({"error": berr.replace("'bases'", "'sequence'")}, 400)
    if not bases:
        return ({"error": "missing/empty 'sequence' (IUPAC DNA)"}, 400)
    enzymes = payload.get("enzymes")
    if enzymes is None and payload.get("enzyme") not in (None, ""):
        enzymes = payload.get("enzyme")
    if isinstance(enzymes, str):
        enzymes = [enzymes]
    if not isinstance(enzymes, list) or not enzymes:
        return ({"error":
                  "'enzymes' must be a non-empty list (or 'enzyme' a string)"},
                400)
    if len(enzymes) > 500:
        return ({"error": "too many enzymes (max 500)"}, 400)
    if not all(isinstance(e, str) and e for e in enzymes):
        return ({"error": "'enzymes' must contain only non-empty strings"}, 400)
    circular = bool(payload.get("circular", True))
    include_seq = bool(payload.get("include_fragment_seq", False))

    # Report names the catalog doesn't know (they contribute no cuts) so a
    # typo'd enzyme is visible instead of silently scanning nothing — the
    # SC-D "looks authoritative while answering a different question" trap.
    # Resolve the same way `list-restriction-sites` does — exact first, then
    # case-insensitively — so `bsai` cuts here too. This endpoint keeps its
    # REPORT-don't-error contract (a caller digesting with a long list wants
    # the cuts it can get plus a list of what was skipped), so only genuinely
    # unrecognised names land in `unknown_enzymes`; the resolved spellings are
    # what the engines actually receive, since they look enzymes up by exact
    # catalog key.
    resolved, unresolved = _resolve_enzyme_names(enzymes)
    unknown = sorted({w[:_AGENT_NAME_ECHO_MAX] for w, _ in unresolved})

    try:
        cuts = _enzyme_cuts(bases, resolved, circular=circular)
        frags = _digest_with_enzymes(bases, resolved, circular=circular)
    except Exception as exc:
        _log.exception("agent digest failed")
        return ({"error": f"digest failed: {_scrub_path(str(exc))}"}, 500)

    cuts_out = [{
        "enzyme":       c.get("enzyme", ""),
        "top":          c.get("top"),
        "bottom":       c.get("bot"),
        "kind":         c.get("kind", ""),
        "overhang_seq": c.get("overhang_seq", ""),
    } for c in cuts]

    frags_out = []
    for i, f in enumerate(frags):
        left  = f.get("left",  {}) or {}
        right = f.get("right", {}) or {}
        top_seq = f.get("top_seq", "") or ""
        entry = {
            "index":  i,
            "length": len(top_seq),
            "left":  {"overhang_seq": left.get("overhang_seq", ""),
                       "kind": left.get("kind", ""),
                       "enzyme": left.get("enzyme", "")},
            "right": {"overhang_seq": right.get("overhang_seq", ""),
                       "kind": right.get("kind", ""),
                       "enzyme": right.get("enzyme", "")},
        }
        if include_seq:
            entry["seq"] = top_seq
        frags_out.append(entry)

    out: dict = {
        "cuts":            cuts_out,
        "fragments":       frags_out,
        "n_cuts":          len(cuts_out),
        "n_fragments":     len(frags_out),
        "circular":        circular,
        "unknown_enzymes": unknown,
    }
    # `cuts[].enzyme` carries the CATALOG name (the engines key on it), so a
    # caller who digested with a commercial synonym would find their own
    # spelling nowhere in the answer and read that as "it didn't cut". Say
    # which name each one became.
    renamed = {}
    for raw in enzymes:
        match = _enzyme_resolve_one(raw)
        if match is not None and match[0] != match[1]:
            renamed[match[1]] = match[0]
    if renamed:
        out["resolved_enzymes"] = renamed
    return out


@_agent_endpoint("diff-plasmid")
def _h_diff_plasmid(app, payload):
    """Pairwise alignment of the loaded record against another plasmid
    in the library. Body: ``{target_id, mode?, circular?}``.

    `circular` (default ``auto`` — detected from the target's topology
    annotation) gates the rotation picker. Pass ``false`` to skip
    rotation; pass ``true`` to force it even when the target's
    topology isn't annotated.

    Runs the same `_pick_best_rotation` pipeline `_diff_align_worker`
    uses (canvas_axis="query"): tries plain forward, plain RC,
    query-rotation (forward + RC), target-rotation (forward + RC).
    Picks whichever crosses ``_ROTATION_EARLY_STOP_PCT`` first or has
    the highest overall identity. Returns the full result plus
    rotation metadata so an agent can replicate the alignment:

      * ``picked_rotation``  — ``"none"`` / ``"query"`` / ``"target"``.
      * ``query_rotation``   — bp offset applied to the query (0 unless
                               ``picked_rotation == "query"``).
      * ``target_rotation``  — bp offset applied to the target (0 unless
                               ``picked_rotation == "target"``).
      * ``inverted_segments`` — spans where the query matches the target's
                               REVERSE COMPLEMENT (an insert cloned the
                               wrong way round), each
                               ``{t_start, t_end, q_start, q_end, length,
                               identity_pct, forward_identity_pct}``.
                               Forward alignment scores such a region like
                               unrelated sequence, so identity alone reads
                               it as plain divergence ([INV-193]).
      * ``query_rc``         — bool: was the query reverse-complemented
                               to find this alignment?
      * ``rotation_offset``  — back-compat field: equals
                               ``target_rotation`` (the legacy single-
                               rotation API only rotated the target).

    Skips the UI and feeds an agent the same numbers
    `AlignmentScreen` surfaces — a "how similar are pUC19 and my new
    construct" question can be answered in one round-trip."""
    rec = getattr(app, "_current_record", None)
    if rec is None:
        return ({"error": "no plasmid loaded"}, 422)
    target_id = _sanitize_label(   # SC-G: accept source_id / id aliases
        payload.get("target_id") or payload.get("source_id")
        or payload.get("id"), max_len=200)
    if not target_id:
        return ({"error": "missing 'target_id'"}, 400)
    mode = payload.get("mode", "global")
    if mode not in ("global", "local"):
        return ({"error": "'mode' must be 'global' or 'local'"}, 400)
    # Sweep #25: helper avoids the per-call full-library deepcopy.
    target_entry = _find_library_entry_by_id(target_id)
    if target_entry is None:
        return ({"error": f"no library entry id={target_id!r}"}, 404)
    gb_text = target_entry.get("gb_text", "")
    if not gb_text:
        return ({"error": "target entry has no gb_text"}, 422)
    try:
        target_record = _gb_text_to_record(gb_text)
    except Exception as exc:
        # Sweep #27: scrub exception text.
        _log.warning("agent diff-plasmid: target parse failed (%s)", exc)
        return ({"error": f"target parse failed: {_scrub_path(str(exc))}"},
                500)
    # Auto-detect circular from topology annotation; agents can override
    # with an explicit `circular` boolean.
    circ_raw = payload.get("circular")
    if circ_raw is None:
        target_annotations = getattr(target_record, "annotations",
                                       None) or {}
        is_circular = (target_annotations.get("topology") == "circular")
    else:
        is_circular = bool(circ_raw)
    query_seq = str(rec.seq)
    target_seq = str(target_record.seq)
    # Pre-cap both seqs at `_PAIRWISE_MAX_LEN` BEFORE the picker runs.
    # `_pick_best_rotation` may double the target inside
    # `_find_circular_alignment_offset` (t + t) — a 50 MB library entry
    # would otherwise allocate 100 MB before reaching the alignment cap.
    if (len(query_seq) > _PAIRWISE_MAX_LEN
            or len(target_seq) > _PAIRWISE_MAX_LEN):
        return ({
            "error": (
                f"sequence too long for diff "
                f"(query={len(query_seq):,}, target={len(target_seq):,}, "
                f"cap={_PAIRWISE_MAX_LEN:,})"
            ),
        }, 413)
    # Use the same picker the UI workers use — INV-72 (2026-05-25)
    # closes the agent/UI divergence: pre-sweep this endpoint ran
    # bare `_pairwise_align` after a single `_find_circular_alignment_
    # offset` call, so agents missed RC-orientation detection +
    # the multi-rotation best-of-N pick that `_diff_align_worker`
    # gained in `[INV-71]`. canvas_axis="query" matches the UI flow:
    # the loaded record is the query, the picked library entry is
    # the target.
    try:
        result = _state._pick_best_rotation_hook(
            query_seq, target_seq,
            is_circular=is_circular,
            mode=mode, canvas_axis="query",
        )
    except ValueError as exc:
        return ({"error": f"alignment rejected: {exc}"}, 400)
    except Exception as exc:
        _log.exception("diff-plasmid: pick_best_rotation raised")
        return ({"error": f"alignment failed: {_scrub_path(str(exc))}"}, 500)
    picked = result.get("picked_rotation", "none")
    target_rotation = int(result.get("target_rotation") or 0)
    return {
        "ok":              True,
        "target_id":       target_id,
        "target_name":     target_record.name or target_record.id or "",
        "circular":        is_circular,
        # Back-compat: the legacy single-rotation API exposed
        # `rotation_offset` as the bp the target was rotated by.
        # Picker may instead rotate the query or pick "none" — in
        # those cases `rotation_offset` is 0 and the new
        # `query_rotation` / `picked_rotation` fields tell the full
        # story.
        "rotation_offset": target_rotation,
        "picked_rotation": picked,
        "query_rotation":  int(result.get("query_rotation") or 0),
        "target_rotation": target_rotation,
        "query_rc":        bool(result.get("query_rc", False)),
        "result":          result,
    }


_VERIFY_READS_MAX = 200   # one alignment per read — cap the batch size


def _agent_resolve_read_reference(payload):
    """Resolve the reference a read-comparison endpoint works against.

    Returns ``(ref_seq, circular, None)`` on success, or
    ``("", True, (error_dict, status))`` to return straight to the caller.

    The error path carries harmless SENTINELS rather than ``None`` so the
    success values keep their concrete types. A caller that checks the third
    element (as every caller does) never reads them, and typing them
    ``str | None`` would force an assert at each callsite to say something the
    control flow already guarantees.

    Accepts a bare ``reference`` sequence, or ``reference_id`` /
    ``reference_name`` naming a saved plasmid (whose stored topology then sets
    ``circular`` unless the body overrides it). Shared by
    `verify-against-reads` and `analyse-read-heterogeneity` rather than copied:
    two resolvers that disagree about what ``circular`` defaults to would make
    the same reads answer differently depending on which endpoint asked, and a
    wrong topology silently re-frames every coordinate in the answer.
    """
    circular_override = payload.get("circular")
    circular = bool(circular_override) if circular_override is not None else True
    ref_raw = payload.get("reference")
    if isinstance(ref_raw, str) and ref_raw.strip():
        ref_seq, e = _sanitize_bases(ref_raw)
        if e:
            return ("", True, ({"error": f"'reference' rejected: {e}"}, 400))
    else:
        key = _sanitize_label(payload.get("reference_id")
                              or payload.get("reference_name"), max_len=200)
        if not key:
            return ("", True,
                    ({"error": "provide 'reference' (a sequence) or "
                               "'reference_id' / 'reference_name'"}, 400))
        entry = (_find_library_entry_by_id(key)
                 or _find_library_entry_by_name(key))
        if entry is None:
            return ("", True,
                    ({"error": f"no library entry named {key!r}"}, 404))
        try:
            rec = _gb_text_to_record(entry.get("gb_text") or "")
        except Exception as exc:
            return ("", True,
                    ({"error": f"couldn't parse {key!r}: "
                               f"{_scrub_path(str(exc))}"}, 500))
        ref_seq = str(getattr(rec, "seq", "") or "")
        if circular_override is None:
            circular = (rec.annotations.get("topology") or "").lower() != "linear"
    if not ref_seq:
        return ("", True, ({"error": "reference sequence is empty"}, 422))
    if len(ref_seq) > _PAIRWISE_MAX_LEN:
        return ("", True,
                ({"error": f"reference exceeds {_PAIRWISE_MAX_LEN:,} bp"}, 413))
    return (ref_seq, circular, None)


@_agent_endpoint("verify-against-reads")
def _h_verify_against_reads(app, payload):
    """Verify a designed construct against sequencing reads — the "I built X, I
    got reads back, do they match?" check. Body:
    ``{reference | reference_id | reference_name, reads: [str, ...],
        circular?, mode?="global", min_identity?=99.0}``.

    ``read_quality`` (optional) is a per-read list of per-base Phred arrays,
    PARALLEL to ``reads`` — pass ``null`` for a read that has none. When given,
    each read also gets a ``quality`` verdict splitting its differences into the
    ones the read actually SUPPORTS and the ones sitting in unreliable
    basecalls, and ``read_quality`` in the response carries them. Both ends of a
    Sanger read are ragged, so without this a good clone can score as divergent
    on instrument noise. ``min_phred`` (default 20) is the cutoff.

    ``reference`` is a bare sequence, OR resolve a saved plasmid via
    ``reference_id`` / ``reference_name`` (its stored topology sets
    ``circular`` unless you pass ``circular``). Each read in ``reads`` (bare
    sequences — Nanopore / Sanger / a Plasmidsaurus consensus) is aligned to
    the reference with the SAME rotation + RC-aware picker the Plasmidsaurus
    aligner and the UI use, and its identity% reported. Returns
    ``{ok, reference_len, n_reads, min_identity, verdict, circular,
    reads:[{index, length, identity_pct, rc, inversions, inverted_bp, passes}],
    summary:{mean, min, max, n_pass, n_fail, n_inverted}}``.

    ``verdict`` is "match" when EVERY read ≥ ``min_identity``; "inverted" when
    every read that falls short does so because part of it matches the
    reference BACKWARDS; else "mismatch" (failing reads flagged
    ``passes:false``). ``inversions`` lists each such span as
    ``{t_start, t_end, q_start, q_end, length, identity_pct,
    forward_identity_pct}`` in reference coordinates — an insert cloned the
    wrong way round otherwise scores like an unrelated plasmid, so identity
    alone can't separate "wrong construct" from "right construct, wrong
    orientation" ([INV-193]). ``rc`` means the WHOLE read aligned as the
    reverse complement. Read-only, heavy — one alignment per read, capped at
    200 reads."""
    mode = payload.get("mode", "global")
    if mode not in ("global", "local"):
        return ({"error": "'mode' must be 'global' or 'local'"}, 400)
    try:
        min_identity = float(payload.get("min_identity", 99.0))
    except (TypeError, ValueError):
        return ({"error": "'min_identity' must be a number"}, 400)
    if not (0.0 <= min_identity <= 100.0):
        return ({"error": "'min_identity' must be in [0, 100]"}, 400)
    ref_seq, circular, ref_err = _agent_resolve_read_reference(payload)
    if ref_err is not None:
        return ref_err
    reads = payload.get("reads")
    if not isinstance(reads, list) or not reads:
        return ({"error": "'reads' must be a non-empty list of sequences"}, 400)
    if len(reads) > _VERIFY_READS_MAX:
        return ({"error": f"too many reads (max {_VERIFY_READS_MAX})"}, 413)
    clean_reads = []
    for i, r in enumerate(reads):
        rs, e = _sanitize_bases(r if isinstance(r, str) else "")
        if e or not rs:
            return ({"error": f"reads[{i}] rejected: {e or 'empty'}"}, 400)
        if len(rs) > _PAIRWISE_MAX_LEN:
            return ({"error": f"reads[{i}] exceeds {_PAIRWISE_MAX_LEN:,} bp"},
                    413)
        clean_reads.append(rs)
    # Per-read Phred arrays, parallel to `reads`. Validated up front so a
    # malformed entry is a 400 naming the read rather than a traceback halfway
    # through the alignment loop.
    read_quals: list = []
    rq = payload.get("read_quality")
    if rq is not None:
        if not isinstance(rq, list):
            return ({"error": "'read_quality' must be a list parallel to "
                              "'reads' (null for a read with no quality)"},
                    400)
        if len(rq) != len(reads):
            return ({"error": f"'read_quality' has {len(rq)} entries for "
                              f"{len(reads)} reads — they must be parallel"},
                    400)
        for qi, one in enumerate(rq):
            if one is None:
                read_quals.append(None)
                continue
            if not isinstance(one, list) or not all(
                    isinstance(v, (int, float)) and not isinstance(v, bool)
                    for v in one):
                return ({"error": f"read_quality[{qi}] must be a list of "
                                  f"numbers or null"}, 400)
            read_quals.append([int(v) for v in one])
    min_phred = _coerce_int(payload.get("min_phred", 20), name="min_phred")
    if isinstance(min_phred, str):
        return ({"error": min_phred}, 400)
    if not (0 <= min_phred <= 93):
        return ({"error": "'min_phred' must be 0-93"}, 400)
    out_reads = []
    _consensus_inputs: list = []
    for i, rs in enumerate(clean_reads):
        try:
            res = _state._pick_best_rotation_hook(
                rs, ref_seq, is_circular=circular, mode=mode,
                canvas_axis="target")
            # Coerce inside the try: a non-numeric identity_pct would otherwise
            # escape the scrubbed-500 guard as a raw float() traceback.
            ident = round(float(res.get("identity_pct") or 0.0), 2)
        except ValueError as exc:
            return ({"error": f"reads[{i}] alignment rejected: {exc}"}, 400)
        except Exception as exc:
            _log.exception("verify-against-reads: alignment raised")
            return ({"error": f"alignment failed: {_scrub_path(str(exc))}"}, 500)
        # Inverted spans [INV-193] — a read whose insert went in backwards
        # scores like an unrelated plasmid, so identity alone reports it as
        # a plain mismatch and the caller cannot tell "wrong construct" from
        # "right construct, wrong way round".
        inv = [s for s in (res.get("inverted_segments") or [])
               if isinstance(s, dict)]
        # Keep the alignment itself so the cross-read rollup below needs no
        # second pass: one read agreeing with itself is not evidence, and the
        # caller cannot combine per-read identities into per-variant support.
        _entry = {"result": res, "axis": "target",
                  "query_label": f"read {i + 1}"}
        # Basecall quality, when the caller supplied it for this read. Computed
        # HERE because the per-base array cannot be persisted anywhere — a
        # library entry is `gb_text` — so it is consumed while in hand and
        # reduced to a verdict the consensus can reuse.
        one_q = read_quals[i] if i < len(read_quals) else None
        if one_q:
            try:
                _entry["quality"] = _state._trace_verification_summary_hook(
                    res.get("aligned_q") or "", res.get("aligned_t") or "",
                    one_q, min_phred=min_phred)
            except Exception:
                _log.exception("verify-against-reads: quality summary failed")
        _consensus_inputs.append(_entry)
        out_reads.append({"index": i, "length": len(rs), "identity_pct": ident,
                          "rc": bool(res.get("query_rc", False)),
                          "inversions": inv,
                          "inverted_bp": sum(
                              max(0, int(s.get("length") or 0)) for s in inv),
                          "passes": ident >= min_identity})
    idents = [r["identity_pct"] for r in out_reads]
    n_pass = sum(1 for r in out_reads if r["passes"])
    n_inverted = sum(1 for r in out_reads if r["inversions"])
    # "inverted" is a THIRD verdict, not a flavour of mismatch: the bases are
    # right and the orientation isn't, which is a different next step at the
    # bench (re-screen colonies) from a failed build (start over). Only
    # reported when every read that falls short does so for THAT reason.
    if n_pass == len(out_reads):
        verdict = "match"
    elif n_inverted and all(r["passes"] or r["inversions"] for r in out_reads):
        verdict = "inverted"
    else:
        verdict = "mismatch"
    _log_event("verify.reads.agent", n_reads=len(out_reads), verdict=verdict,
               min_identity=min_identity, n_inverted=n_inverted, via="agent")
    return {"ok": True, "reference_len": len(ref_seq),
            "n_reads": len(out_reads), "min_identity": min_identity,
            "verdict": verdict, "circular": circular, "reads": out_reads,
            "summary": {"mean": round(sum(idents) / len(idents), 2),
                        "min": min(idents), "max": max(idents),
                        "n_pass": n_pass, "n_fail": len(out_reads) - n_pass,
                        "n_inverted": n_inverted},
            # Cross-read rollup: which differences the reads AGREE on, how much
            # of the reference they actually covered, and what nobody read.
            # Per-read identity cannot answer any of those — two reads at 99.9%
            # might disagree about entirely different bases.
            "consensus": _state._multi_read_summary_hook(
                _consensus_inputs, len(ref_seq)),
            "read_quality": [e.get("quality") for e in _consensus_inputs],
            "ignored": _agent_ignored_keys(payload, {
                "reference", "reference_id", "reference_name", "reads",
                "circular", "mode", "min_identity", "read_quality",
                "min_phred"})}


# ── Sub-consensus read heterogeneity ──────────────────────────────────
# `read-consensus` answers "do my reads agree with each other", and
# `verify-against-reads` answers "do my reads match the design". Neither
# answers "is what I sequenced ONE thing" — a consensus is a single read, so a
# population of escapers each carrying a different inactivating mutation,
# none of them dominant, consenses back to wild type and reads as clean.
_HETERO_READS_MAX = 400          # one full alignment per read — bound the batch
# A read count alone does NOT bound the work: the rotation-aware aligner costs
# roughly 0.18 s per read against a 3 kb reference and 4.5 s against 30 kb
# (measured), so 400 reads on a 30 kb molecule is half an hour in ONE request
# — long enough to pin a worker thread past any client timeout. The real
# budget is reads x reference length, and it is enforced rather than hoped for.
#
# 1,000,000 bp-reads works out to ~330 reads on a 3 kb plasmid, ~100 on 10 kb
# and ~33 on 30 kb, each landing in the 60-150 s range. `max_reads` can only
# LOWER this, never raise it; a caller who needs more depth on a large
# molecule subsamples and calls again. Whenever the budget bites, the response
# says so in `capped_by` and in a warning — a quietly truncated read set would
# change every allele fraction in the answer without saying why.
_HETERO_WORK_BUDGET_BP = 1_000_000
_HETERO_MIN_FRACTION_DEFAULT = 0.01     # 1% — below this is basecall noise
# A position this much of the population disagrees at is not a "minor" variant
# any more; the sample is substantially not one thing.
_HETERO_MIXED_FRACTION = 0.25

# ── Platform noise, and why the verdict cannot be a COUNT ─────────────
# The first cut called `mixed` when 5+ positions cleared 1%. Measured against a
# simulated but realistic CLONAL nanopore run (30 reads, 3 kb, homopolymer-
# biased indels), that rule produced: verdict `mixed`, 1,347 positions, 78% of
# them inside a homopolymer, 72% indels, top fraction 0.304 — i.e. it cleared
# the count rule AND the 0.25 dominant-fraction rule, on a sample that was one
# clone. A heterogeneity call that says `mixed` for every long-read dataset is
# worth exactly as much as one that always says `clean`.
#
# Two things fix it, and both are about the SHAPE of the evidence:
#
#   1. HOMOPOLYMER CONTEXT. Nanopore miscalls run LENGTH, so indels pile into
#      homopolymers. Those observations are still REPORTED (a real frameshift
#      can live in a homopolymer, and hiding it would be worse), flagged
#      `context: "homopolymer"` — they just don't get a vote, because the
#      platform cannot distinguish them from its own error mode.
#   2. DISTRIBUTION SHAPE, not a count. A genuine sub-population is a
#      HAPLOTYPE: every read from the escaper carries ALL of that clone's
#      differences, so those positions share one fraction and appear as a MODE.
#      Platform noise is independent per position and decays monotonically from
#      the floor. So the test is "is there a bump", not "are there many".
_HETERO_PLATFORMS = ("ont", "illumina", "sanger", "unknown")
# Per-platform defaults: (noise_floor, homopolymer_filter).
# `noise_floor` is the fraction below which a call does not get a VOTE in the
# verdict. It is not `min_fraction` — that still controls what is REPORTED, so
# lowering the reporting floor never silently changes the verdict.
_HETERO_PLATFORM_PROFILE = {
    # Long reads: ~1-2% substitution but heavy homopolymer indel noise, which
    # on a 30x pile reaches 20-30% at individual positions.
    "ont":      (0.10, True),
    # Short reads: substitutions dominate, indels are rare and homopolymer
    # slippage is far less severe.
    "illumina": (0.02, False),
    # Few reads, ragged ends already handled by `min_phred`; slippage exists
    # but a Sanger pile is too shallow for a fraction to mean much below 5%.
    "sanger":   (0.05, False),
    # Caller did not say. A homopolymer INDEL is a suspect call on every
    # platform — all three slip in runs, long reads just far more — so it is
    # muted here too rather than letting an undeclared nanopore run read
    # `mixed` on its own error mode. Nothing is hidden: those calls are still
    # in `positions` with `context: "homopolymer"`, counted in
    # `n_muted_homopolymer`, and a warning fires when muting removed anything.
    # Screening for a frameshift IN a homopolymer tract is the one case that
    # wants `homopolymer_filter: false`, and the warning says so.
    "unknown":  (0.05, True),
}
# Homopolymer run length at which a platform starts miscalling it.
_HETERO_HOMOPOLYMER_MIN_RUN = 4
# Shape test. Fractions above the noise floor are binned; noise decays
# monotonically, a haplotype makes a bump. A bin must beat the one below it by
# this many observations to count as a mode rather than binning jitter.
_HETERO_SHAPE_BIN = 0.05
_HETERO_MODE_MIN_EXCESS = 2
# A mode of one position is not a mode.
_HETERO_MODE_MIN_POSITIONS = 2


def _hetero_homopolymer_run(seq: str, pos: int, *, circular: bool) -> int:
    """Length of the homopolymer run the reference carries AT `pos`.

    1 for an isolated base, 0 for an out-of-range or ambiguous one. Walks both
    ways from `pos`, wrapping on a circular molecule (a run straddling the
    origin is still one run) and bounded by the sequence length so a
    single-base molecule cannot loop forever.
    """
    n = len(seq)
    if n == 0 or not (0 <= pos < n):
        return 0
    b = seq[pos]
    if b not in "ACGT":
        return 0
    run = 1
    i = pos
    while run < n:
        j = i - 1
        if j < 0:
            if not circular:
                break
            j = n - 1
        if seq[j] != b:
            break
        run += 1
        i = j
    i = pos
    while run < n:
        j = i + 1
        if j >= n:
            if not circular:
                break
            j = 0
        if seq[j] != b:
            break
        run += 1
        i = j
    return run


def _hetero_homopolymer_context(seq: str, pos: int, *, circular: bool,
                                min_run: int) -> int:
    """Longest homopolymer run touching `pos` or either neighbour, else 0.

    The neighbours matter: a run-length miscall reports its indel at the run's
    EDGE as often as inside it, so testing only `pos` misses half of them.
    """
    n = len(seq)
    if n == 0:
        return 0
    best = 0
    for d in (-1, 0, 1):
        p = pos + d
        if circular:
            p %= n
        elif not (0 <= p < n):
            continue
        best = max(best, _hetero_homopolymer_run(seq, p, circular=circular))
    return best if best >= min_run else 0


def _hetero_shape(fractions: "list[float]", floor: float = 0.0,
                  max_fraction: "float | None" = None) -> dict:
    """Is this distribution a decaying noise tail, or does it carry a MODE?

    Returns ``{n, max_fraction, bins, monotonic_decay, mode_fraction}``.

    Independent per-position error decays monotonically away from the floor.
    A real sub-population is linked — every read of the escaper carries all of
    its differences — so its positions share one fraction and stack into a bin
    that stands ABOVE the bin below it. `mode_fraction` is the centre of the
    lowest such bin, or None when the distribution just decays.

    **`fractions` must include the calls BELOW the noise floor**, because they
    are the baseline the mode is judged against. Binning only the above-floor
    calls makes the first bin empty by construction, so anything just above
    the floor looks like a mode standing over nothing — which made a clonal
    nanopore pile read `mixed` on its own substitution noise. Conversely,
    discarding the sub-floor tail loses the decay that tells eight escaper
    mutations sharing one fraction apart from a noise skirt.

    A mode is only honoured at or above `floor`: below it, a bump is noise
    structure, not a sub-population worth reporting. `max_fraction` reports
    the VOTING maximum when one is supplied, since that is what the dominant
    test keys on.
    """
    if not fractions:
        return {"n": 0, "max_fraction": round(float(max_fraction or 0.0), 4),
                "max_observed_fraction": 0.0,
                "bins": [], "monotonic_decay": True, "mode_fraction": None}
    top = max(fractions)
    n_bins = max(1, int(top / _HETERO_SHAPE_BIN) + 1)
    counts = [0] * n_bins
    for f in fractions:
        counts[min(n_bins - 1, max(0, int(f / _HETERO_SHAPE_BIN)))] += 1
    mode = None
    for i in range(1, n_bins):
        # A bin's centre can sit past 1.0 for the topmost bin (a fraction of
        # exactly 1.0 lands in [1.00, 1.05)), and a reported "mode_fraction"
        # above 1 is not a fraction. Clamp to the bin's own upper edge.
        centre = min(1.0, (i + 0.5) * _HETERO_SHAPE_BIN)
        if centre < floor:
            continue
        if (counts[i] >= _HETERO_MODE_MIN_POSITIONS
                and counts[i] - counts[i - 1] >= _HETERO_MODE_MIN_EXCESS):
            mode = round(centre, 4)
            break
    bins = [{"from": round(i * _HETERO_SHAPE_BIN, 4),
             "to": round(min(1.0, (i + 1) * _HETERO_SHAPE_BIN), 4),
             "n": counts[i]} for i in range(n_bins)]
    return {"n": len(fractions),
            # The VOTING maximum — what the dominant test keys on — alongside
            # the maximum actually present in the binned baseline. Reporting
            # only the former left `max_fraction: 0.0` sitting next to
            # populated bins, which a caller cannot reconcile.
            "max_fraction": round(float(top if max_fraction is None
                                        else max_fraction), 4),
            "max_observed_fraction": round(float(top), 4),
            "bins": bins,
            "monotonic_decay": mode is None, "mode_fraction": mode}


@_agent_endpoint("analyse-read-heterogeneity")
def _h_analyse_read_heterogeneity(app, payload):
    """Per-base ALLELE FRACTIONS from raw reads — "is this culture one thing?"

    Body: ``{reference | reference_id | reference_name,
    reads_path | reads: [str, ...], read_quality?, min_fraction?=0.01,
    min_phred?=20, max_reads?, circular?, mode?="global"}``.

    `reads_path` points at a FASTQ (`.fastq` / `.fq`); its per-base Phred
    arrays are read straight off the file, which is the normal way to call
    this. `reads` passes sequences inline instead, with optional parallel
    `read_quality` arrays — the same shape `verify-against-reads` takes.

    **Why this is not `read-consensus`.** That endpoint reads the alignments
    STORED on a plasmid, and a Plasmidsaurus consensus is one read: depth 1
    everywhere, so anything below the consensus is invisible. A heterogeneous
    population — many different mutations, none dominant — consenses to wild
    type and reads as a clean, clonal culture. Only the raw reads carry that
    signal, and only as a FRACTION: `alt_reads / depth`.

    Returns ``{ok, reference_len, n_reads, reads_source, circular, mean_depth,
    covered_bp, covered_pct, uncovered_spans, verdict, positions, thresholds,
    n_positions, n_filtered_low_quality, dropped_reads}`` where each entry of
    `positions` is ``{pos, type, ref, alt, alt_reads, depth, fraction,
    mean_phred, contradicted}``, sorted by descending fraction.

    `verdict` is:
      * ``clonal`` — nothing above `min_fraction`; the reads are one sequence.
      * ``minor_variants`` — a sub-population exists but is small and focal.
      * ``mixed`` — a position at/above 25% of the reads, OR 5+ distinct
        positions carrying real sub-populations (the escaper signature).
      * ``no_reads`` — nothing aligned.

    **`fraction`'s denominator is every read that COVERED the base**, including
    the ones that disagreed, so a variant inside a region only two reads reached
    cannot read as 100%. `mean_phred` is the mean basecall quality of the reads
    SUPPORTING the call; observations whose backing basecall is below
    `min_phred` are excluded from `alt_reads` and counted in
    `n_filtered_low_quality`, because at 1% every instrument's error floor
    otherwise looks like a sub-population. An observation with no quality at
    all is counted (unknown is not the same as bad).

    The thresholds that produced `verdict` are echoed in `thresholds` — they
    are conventions, not a calibrated model, and a caller re-deriving its own
    call from `positions` should be able to see exactly what these were.
    Read-only, heavy — one alignment per read."""
    ref_seq, circular, ref_err = _agent_resolve_read_reference(payload)
    if ref_err is not None:
        return ref_err
    mode = payload.get("mode", "global")
    if mode not in ("global", "local"):
        return ({"error": "'mode' must be 'global' or 'local'"}, 400)
    try:
        min_fraction = float(payload.get("min_fraction",
                                         _HETERO_MIN_FRACTION_DEFAULT))
    except (TypeError, ValueError):
        return ({"error": "'min_fraction' must be a number"}, 400)
    if not (0.0 < min_fraction <= 1.0):
        return ({"error": "'min_fraction' must be in (0, 1]"}, 400)
    min_phred = _coerce_int(payload.get("min_phred", 20), name="min_phred")
    if isinstance(min_phred, str):
        return ({"error": min_phred}, 400)
    if not (0 <= min_phred <= 93):
        return ({"error": "'min_phred' must be 0-93"}, 400)
    platform = payload.get("platform", "unknown")
    if not isinstance(platform, str) or \
            platform.strip().lower() not in _HETERO_PLATFORMS:
        return ({"error": f"'platform' must be one of "
                          f"{', '.join(_HETERO_PLATFORMS)}"}, 400)
    platform = platform.strip().lower()
    _default_floor, _default_hp = _HETERO_PLATFORM_PROFILE[platform]
    noise_floor = payload.get("noise_floor")
    if noise_floor is None:
        noise_floor = _default_floor
    else:
        try:
            noise_floor = float(noise_floor)
        except (TypeError, ValueError):
            return ({"error": "'noise_floor' must be a number"}, 400)
        if not (0.0 <= noise_floor <= 1.0):
            return ({"error": "'noise_floor' must be in [0, 1]"}, 400)
    hp_filter = payload.get("homopolymer_filter")
    if hp_filter is None:
        hp_filter = _default_hp
    else:
        hp_filter = _coerce_bool(hp_filter, name="homopolymer_filter")
        if isinstance(hp_filter, str):
            return ({"error": hp_filter}, 400)
    hp_min_run = _coerce_int(
        payload.get("homopolymer_min_run", _HETERO_HOMOPOLYMER_MIN_RUN),
        name="homopolymer_min_run")
    if isinstance(hp_min_run, str):
        return ({"error": hp_min_run}, 400)
    if not (2 <= hp_min_run <= 50):
        return ({"error": "'homopolymer_min_run' must be 2-50"}, 400)
    max_reads = _coerce_int(payload.get("max_reads", _HETERO_READS_MAX),
                            name="max_reads")
    if isinstance(max_reads, str):
        return ({"error": max_reads}, 400)
    if not (1 <= max_reads <= _HETERO_READS_MAX):
        return ({"error": f"'max_reads' must be 1-{_HETERO_READS_MAX}"}, 400)
    # Work budget — see `_HETERO_WORK_BUDGET_BP`. Applied here so it bounds
    # BOTH read sources, and only ever downwards.
    budget_reads = max(1, _HETERO_WORK_BUDGET_BP // max(1, len(ref_seq)))
    capped_by = None
    if budget_reads < max_reads:
        max_reads = budget_reads
        capped_by = "work_budget"
    elif payload.get("max_reads") is not None:
        capped_by = "max_reads"

    # ── Reads: a FASTQ on disk, or inline sequences ──────────────────
    reads_path = payload.get("reads_path")
    inline = payload.get("reads")
    if reads_path is not None and inline is not None:
        return ({"error": "pass 'reads_path' OR 'reads', not both"}, 400)
    seqs: "list[str]" = []
    quals: "list" = []
    dropped: "list[dict]" = []
    if reads_path is not None:
        p = _sanitize_path(reads_path)
        if p is None:
            return ({"error": "'reads_path' is not a usable path"}, 400)
        if not p.exists():
            return ({"error": f"no such file: {_scrub_path(str(p))}"}, 404)
        if p.suffix.lower() not in (".fastq", ".fq"):
            return ({"error": f"'reads_path' must be .fastq / .fq "
                              f"(got {p.suffix or 'no extension'!r}). Per-base "
                              f"quality is what separates a real 1% variant "
                              f"from the instrument's error floor, and only "
                              f"FASTQ carries it."}, 400)
        try:
            recs = _fastq_path_to_records(str(p))
        except ValueError as exc:
            return ({"error": f"could not read FASTQ: "
                              f"{_scrub_path(str(exc))}"}, 400)
        except OSError as exc:
            return ({"error": f"could not read FASTQ: "
                              f"{_scrub_path(str(exc))}"}, 500)
        source = "fastq"
        for i, rec in enumerate(recs[:max_reads]):
            rs, e = _sanitize_bases(str(getattr(rec, "seq", "") or ""))
            if e or not rs:
                dropped.append({"index": i, "reason": e or "empty"})
                continue
            if len(rs) > _PAIRWISE_MAX_LEN:
                dropped.append({"index": i, "reason": "over length cap"})
                continue
            seqs.append(rs)
            q = (getattr(rec, "letter_annotations", {}) or {}).get(
                "phred_quality")
            quals.append([int(v) for v in q] if q else None)
        n_available = len(recs)
    else:
        if not isinstance(inline, list) or not inline:
            return ({"error": "provide 'reads_path' (a FASTQ) or 'reads' "
                              "(a non-empty list of sequences)"}, 400)
        source = "inline"
        rq = payload.get("read_quality")
        if rq is not None:
            if not isinstance(rq, list):
                return ({"error": "'read_quality' must be a list parallel to "
                                  "'reads' (null for a read with no quality)"},
                        400)
            if len(rq) != len(inline):
                return ({"error": f"'read_quality' has {len(rq)} entries for "
                                  f"{len(inline)} reads — they must be "
                                  f"parallel"}, 400)
        for i, r in enumerate(inline[:max_reads]):
            rs, e = _sanitize_bases(r if isinstance(r, str) else "")
            if e or not rs:
                return ({"error": f"reads[{i}] rejected: {e or 'empty'}"}, 400)
            if len(rs) > _PAIRWISE_MAX_LEN:
                return ({"error": f"reads[{i}] exceeds "
                                  f"{_PAIRWISE_MAX_LEN:,} bp"}, 413)
            seqs.append(rs)
            one = rq[i] if isinstance(rq, list) and i < len(rq) else None
            if one is None:
                quals.append(None)
                continue
            if not isinstance(one, list) or not all(
                    isinstance(v, (int, float)) and not isinstance(v, bool)
                    for v in one):
                return ({"error": f"read_quality[{i}] must be a list of "
                                  f"numbers or null"}, 400)
            quals.append([int(v) for v in one])
        n_available = len(inline)

    thresholds = {"min_fraction": min_fraction, "min_phred": min_phred,
                  "mixed_fraction": _HETERO_MIXED_FRACTION,
                  "platform": platform, "noise_floor": noise_floor,
                  "homopolymer_filter": hp_filter,
                  "homopolymer_min_run": hp_min_run,
                  "shape_bin": _HETERO_SHAPE_BIN}
    hetero_warnings: "list[str]" = []
    if n_available > max_reads:
        hetero_warnings.append(
            f"{n_available} reads available but only {max_reads} were used"
            + (f" — the {_HETERO_WORK_BUDGET_BP:,} bp-read work budget allows "
               f"{max_reads} against a {len(ref_seq):,} bp reference. Every "
               f"fraction below is computed from those {max_reads} reads."
               if capped_by == "work_budget"
               else f" ('max_reads'). Every fraction below is computed from "
                    f"those {max_reads} reads."))
    if dropped:
        hetero_warnings.append(f"{len(dropped)} read(s) were unusable and "
                               f"excluded; see 'dropped_reads'")
    base = {"ok": True, "reference_len": len(ref_seq), "circular": circular,
            "reads_source": source, "reads_available": n_available,
            "max_reads": max_reads, "capped_by": capped_by,
            "work_budget_bp": _HETERO_WORK_BUDGET_BP,
            "warnings": hetero_warnings,
            "thresholds": thresholds,
            "dropped_reads": dropped,
            "ignored": _agent_ignored_keys(payload, {
                "reference", "reference_id", "reference_name", "reads",
                "reads_path", "read_quality", "min_fraction", "min_phred",
                "max_reads", "circular", "mode", "platform", "noise_floor",
                "homopolymer_filter", "homopolymer_min_run"})}
    if not seqs:
        return {**base, "n_reads": 0, "mean_depth": 0.0, "covered_bp": 0,
                "covered_pct": 0.0,
                "uncovered_spans": [[0, len(ref_seq)]] if ref_seq else [],
                "verdict": "no_reads", "positions": [], "n_positions": 0,
                "n_filtered_low_quality": 0, "n_no_call_observations": 0,
                "n_counted_in_verdict": 0, "n_muted_homopolymer": 0,
                "n_below_noise_floor": 0, "platform": platform,
                "distribution": _hetero_shape([])}

    # ── Align every read, then reuse the cross-read rollup ───────────
    aligned: "list[dict]" = []
    for i, rs in enumerate(seqs):
        try:
            res = _state._pick_best_rotation_hook(
                rs, ref_seq, is_circular=circular, mode=mode,
                canvas_axis="target")
        except ValueError as exc:
            return ({"error": f"reads[{i}] alignment rejected: {exc}"}, 400)
        except Exception as exc:
            _log.exception("analyse-read-heterogeneity: alignment raised")
            return ({"error": f"alignment failed: "
                              f"{_scrub_path(str(exc))}"}, 500)
        aligned.append({"result": res, "axis": "target",
                        "query_label": f"read {i + 1}"})
    # `_multi_read_summary` owns the per-bp depth sweep, the coverage and the
    # uncovered spans — the denominator of every fraction below. Reused rather
    # than re-derived: a second depth sweep that disagreed with the one
    # `read-consensus` reports would make the two endpoints contradict each
    # other about the same reads.
    roll = _state._multi_read_summary_hook(aligned, len(ref_seq))
    depth_of = {(int(v["target_pos"]), v["type"], v.get("alt", "")): v
                for v in (roll.get("variants") or [])}

    # ── Quality-filtered per-position tally ──────────────────────────
    # Support is recounted HERE rather than taken from the rollup because a
    # fraction is only as meaningful as the basecalls behind it: at 1% the
    # error floor of every instrument looks like a sub-population.
    tally: dict = {}
    n_filtered = 0
    n_no_call = 0
    for i, entry in enumerate(aligned):
        res = entry["result"]
        aq = res.get("aligned_q") or ""
        at = res.get("aligned_t") or ""
        q = quals[i] if i < len(quals) else None
        for v in _state._alignment_variants_in_axis_hook(entry):
            if v.get("type") == "truncated":
                continue
            key = (int(v.get("target_pos", 0) or 0), v.get("type"),
                   v.get("alt", ""))
            if key not in depth_of:
                # Outside the read's observed extent — the rollup already
                # excluded it (a global alignment reports the unread remainder
                # as one long deletion). Honour that same decision.
                continue
            # A NO-CALL IS NOT AN ALLELE. A basecaller that emits `N` is saying
            # it could not read that base, and counting it as a difference
            # turns one ordinary read with a ten-base N run into a ten-position
            # "10% sub-population" and the whole sample into `mixed`. Real
            # long-read data carries N runs routinely, so this fired on almost
            # every genuine dataset. An ambiguity code in the REFERENCE is the
            # same problem from the other side: every read then looks variant
            # against a base the reference does not actually specify.
            # Counted, not silently dropped — ambiguous coverage is a fact
            # about the data the caller should see.
            _alt = str(v.get("alt", "") or "")
            _ref = str(v.get("ref", "") or "")
            if (_alt and any(c not in "ACGT" for c in _alt)) or \
                    (_ref and any(c not in "ACGT" for c in _ref)):
                n_no_call += 1
                continue
            ph = None
            if q:
                try:
                    ph = _state._variant_phred_hook(v, aq, at, q)
                except Exception:
                    _log.exception("analyse-read-heterogeneity: phred lookup "
                                   "failed")
                    ph = None
            if ph is not None and ph < min_phred:
                n_filtered += 1
                continue
            slot = tally.setdefault(key, {"n": 0, "phred_sum": 0,
                                          "phred_n": 0})
            slot["n"] += 1
            if ph is not None:
                slot["phred_sum"] += ph
                slot["phred_n"] += 1

    positions: "list[dict]" = []
    for key, slot in tally.items():
        pos, vtype, alt = key
        ref_v = depth_of[key]
        d = int(ref_v.get("depth") or 0)
        if d <= 0:
            continue
        frac = slot["n"] / d
        if frac < min_fraction:
            continue
        # Homopolymer context. Reported for EVERY call so a caller can see it
        # regardless of platform; it only removes a vote when the filter is on
        # AND the call is an indel, which is the error mode it describes. A
        # SUBSTITUTION inside a homopolymer is not a run-length miscall and
        # keeps its vote.
        hp_len = _hetero_homopolymer_context(ref_seq, pos, circular=circular,
                                             min_run=hp_min_run)
        is_indel = vtype != "snp"
        muted = bool(hp_filter and hp_len and is_indel)
        positions.append({
            "pos":          pos,
            "type":         vtype,
            "ref":          ref_v.get("ref", ""),
            "alt":          alt,
            "alt_reads":    slot["n"],
            "depth":        d,
            "fraction":     round(frac, 4),
            "contradicted": max(0, d - slot["n"]),
            "mean_phred":   (round(slot["phred_sum"] / slot["phred_n"], 1)
                             if slot["phred_n"] else None),
            "context":      ("homopolymer" if hp_len else None),
            "homopolymer_len": hp_len or None,
            # Whether this call got a VOTE in the verdict, and if not, WHY.
            # Reported per position so "why is this mixed / not mixed" is
            # answerable from the response alone rather than by re-deriving
            # the thresholds. The homopolymer filter is checked first because
            # it is the reason a caller can turn off.
            "counted_in_verdict": bool(not muted and frac >= noise_floor),
            "muted_by": ("homopolymer" if muted
                         else "noise_floor" if frac < noise_floor
                         else None),
        })
    positions.sort(key=lambda p: (-p["fraction"], p["pos"]))

    n_pos = len(positions)
    # ── Verdict: SHAPE, not count ────────────────────────────────────
    # Only calls that earned a vote are judged. `mixed` needs either one
    # substantial sub-population (a position at/above `mixed_fraction`) or a
    # MODE in the distribution — a bin standing above the one below it, which
    # is what a linked haplotype looks like and what independent per-position
    # error does not. A long decaying tail of low-fraction calls is platform
    # noise no matter how many entries it has.
    voting = [p for p in positions if p["counted_in_verdict"]]
    # The shape baseline is every call the platform filter did NOT mute —
    # including the sub-floor tail, which is what makes a decay recognisable
    # as a decay. Homopolymer-muted calls are left out: they are the
    # instrument's error mode, and letting them pad the low bins would mask a
    # real cluster sitting just above them.
    shape = _hetero_shape(
        [p["fraction"] for p in positions if p["muted_by"] != "homopolymer"],
        noise_floor,
        max_fraction=max((p["fraction"] for p in voting), default=0.0))
    n_muted_hp = sum(1 for p in positions if p["muted_by"] == "homopolymer")
    n_below_floor = sum(1 for p in positions
                        if p["muted_by"] == "noise_floor")
    if not voting:
        verdict = "clonal"
    elif (shape["max_fraction"] >= _HETERO_MIXED_FRACTION
            or shape["mode_fraction"] is not None):
        verdict = "mixed"
    else:
        verdict = "minor_variants"

    # Muting is never silent: it is the one thing that can turn a genuinely
    # mixed sample into `clonal`, and a frameshift escaper in a polyA tract is
    # exactly the case that lives in a homopolymer.
    if n_muted_hp:
        hetero_warnings.append(
            f"{n_muted_hp} indel call(s) inside a homopolymer run were "
            f"reported but NOT counted towards the verdict — every platform "
            f"miscalls run length, long reads most of all. They are in "
            f"`positions` with context 'homopolymer'. If you are screening "
            f"for a frameshift INSIDE a homopolymer tract, re-run with "
            f"homopolymer_filter: false.")
    # The noise floor is the other way a real sub-population can vanish, and
    # it fails in the DANGEROUS direction: a false `mixed` costs a re-screen,
    # a false `clonal` lets an escaper through. So whenever the floor removed
    # something and the answer was not `mixed`, say what the biggest excluded
    # call was — that is the number that tells a caller whether anything
    # worth looking at was dropped.
    if n_below_floor and verdict != "mixed":
        _top_cut = max((p["fraction"] for p in positions
                        if p["muted_by"] == "noise_floor"), default=0.0)
        hetero_warnings.append(
            f"{n_below_floor} call(s) were reported but fell below the "
            f"{noise_floor:g} noise floor for platform '{platform}' and did "
            f"not count towards the verdict — the highest was {_top_cut:g}. "
            f"If this is short-read data, pass platform: \"illumina\" (floor "
            f"{_HETERO_PLATFORM_PROFILE['illumina'][0]:g}) or set "
            f"noise_floor explicitly.")
    # The nanopore signature, measured rather than assumed: if most calls are
    # indels sitting in homopolymers, this is long-read data and the caller
    # should say so — the platform profile changes the noise floor, and an
    # undeclared run gets judged against a floor built for short reads.
    if platform == "unknown" and n_pos >= 20:
        _sig = sum(1 for p in positions
                   if p["type"] != "snp" and p["context"] == "homopolymer")
        if _sig >= 0.5 * n_pos:
            hetero_warnings.append(
                f"{_sig * 100 // n_pos}% of the calls are indels inside "
                f"homopolymers — the signature of long-read data. Pass "
                f"platform: \"ont\" so the noise floor matches the "
                f"instrument; judged as 'unknown' this verdict may overstate "
                f"heterogeneity.")

    total_bp = int(roll.get("total_bp") or len(ref_seq))
    covered_bp = int(roll.get("covered_bp") or 0)
    _log_event("analyse.heterogeneity", n_reads=len(seqs), verdict=verdict,
               n_positions=n_pos, source=source, via="agent")
    return {**base, "n_reads": len(seqs),
            "mean_depth": roll.get("mean_depth", 0.0),
            "max_depth": roll.get("max_depth", 0),
            "covered_bp": covered_bp,
            "covered_pct": roll.get("covered_pct", 0.0),
            "uncovered_spans": roll.get("uncovered_spans") or [],
            "total_bp": total_bp,
            "verdict": verdict, "positions": positions, "n_positions": n_pos,
            "n_filtered_low_quality": n_filtered,
            "n_no_call_observations": n_no_call,
            # The verdict's working, so it can be checked rather than trusted.
            "n_counted_in_verdict": len(voting),
            "n_muted_homopolymer": n_muted_hp,
            "n_below_noise_floor": n_below_floor,
            "distribution": shape,
            "platform": platform}


@_agent_endpoint("align-plasmidsaurus-zip")
def _h_align_plasmidsaurus_zip(app, payload):
    """Align a Plasmidsaurus zip member against a library plasmid.
    Body: ``{path, member, target_id?, target_name?, mode?, circular?}``.

    Either `target_id` (library entry id) or `target_name` is
    required. `mode` is ``global`` (default) or ``local``.
    `circular` defaults to the target's topology annotation;
    passing it explicitly overrides the auto-detect.

    Runs the same pipeline `_align_worker` uses in the UI:
    extract the `.gbk` member from the zip, parse it as a
    SeqRecord, find the circular alignment offset against the
    target via `_find_circular_alignment_offset`, run
    `_pairwise_align`. Returns the alignment result plus the
    rotation offset so the agent can map matches back to the
    target's original coordinates.

    Read-only — does not register an overlay in the UI or persist
    anything. The library is consulted to resolve the target only.
    """
    raw_path = payload.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        return ({"error": "missing 'path'"}, 400)
    member = payload.get("member")
    if not isinstance(member, str) or not member:
        return ({"error": "missing 'member' (zip member filename)"}, 400)
    target_id_raw = payload.get("target_id")
    target_name_raw = payload.get("target_name")
    target_id = (_sanitize_label(target_id_raw, max_len=200)
                 if isinstance(target_id_raw, str) else "")
    target_name = (_sanitize_label(target_name_raw, max_len=200)
                   if isinstance(target_name_raw, str) else "")
    if not target_id and not target_name:
        return ({"error": "missing 'target_id' or 'target_name'"}, 400)
    mode = payload.get("mode", "global")
    if mode not in ("global", "local"):
        return ({"error": "'mode' must be 'global' or 'local'"}, 400)
    path = _sanitize_path(raw_path)
    if path is None:
        return ({"error": "could not sanitize 'path'"}, 400)
    # Sweep #26 (2026-05-23): defense-in-depth ancestor symlink walk
    # (see `_h_list_plasmidsaurus_members` rationale).
    anc_err = _check_agent_read_path_ancestors(path)
    if anc_err is not None:
        _log.warning("agent align-plasmidsaurus-zip: %s", anc_err)
        return ({"error": "zip rejected (see splicecraft log)"}, 400)
    # Sweep #25 (2026-05-23): collapse path errors to a uniform 400
    # (see `_h_list_plasmidsaurus_members` rationale).
    ok, reason = _safe_file_size_check(
        path, _PLASMIDSAURUS_ZIP_MAX_BYTES, "Plasmidsaurus zip",
    )
    if not ok:
        _log.warning("agent align-plasmidsaurus-zip: rejected (%s): %s",
                     reason or "unsafe", path)
        return ({"error": "zip rejected (see splicecraft log)"}, 400)
    # Resolve the target in the active library. Cross-collection
    # search isn't this endpoint's job — agents should
    # `set-active-collection` first when the target lives elsewhere.
    # Sweep #25: id-then-name helpers avoid the per-call full-library
    # deepcopy. id takes precedence (matches pre-sweep semantics —
    # `target_id` was checked first on each iteration).
    target_entry: "dict | None" = None
    if target_id:
        target_entry = _find_library_entry_by_id(target_id)
    if target_entry is None and target_name:
        target_entry = _find_library_entry_by_name(target_name)
    if target_entry is None:
        descriptor = (f"id={target_id!r}" if target_id
                      else f"name={target_name!r}")
        return ({"error": f"no library entry matching {descriptor}"},
                404)
    # Sweep #26 (2026-05-23): capture the resolved id so the post-
    # alignment re-check can detect a concurrent library mutation
    # (target deleted/renamed mid-flight). The alignment itself runs
    # against the snapshot we extracted above, but the result's
    # `target_name` would otherwise reflect the pre-rename name with
    # no way for the agent to detect the drift.
    resolved_target_id = target_entry.get("id") or ""
    gb_text_target = target_entry.get("gb_text") or ""
    if not gb_text_target:
        return ({"error": "target entry has no gb_text"}, 422)
    try:
        target_record = _gb_text_to_record(gb_text_target)
    except Exception as exc:
        # Sweep #27: scrub exception text on all three parse paths.
        _log.warning("agent align-plasmidsaurus: target parse failed (%s)",
                      exc)
        return ({"error": f"target parse failed: {_scrub_path(str(exc))}"},
                500)
    # Extract + parse the zip member.
    try:
        gb_text_query = _extract_gbk_member(path, member)
    except ValueError as exc:
        _log.warning("agent align-plasmidsaurus: extract failed (%s)", exc)
        return ({"error": "could not extract member (see splicecraft log)"},
                422)
    except Exception as exc:
        _log.exception("agent align-plasmidsaurus-zip: extract failed")
        return ({"error": "extract failed (see splicecraft log)",
                 "type":  type(exc).__name__}, 500)
    try:
        query_record = _gb_text_to_record(gb_text_query)
    except Exception as exc:
        _log.warning("agent align-plasmidsaurus: query parse failed (%s)",
                      exc)
        return ({"error": f"query parse failed: {_scrub_path(str(exc))}"},
                422)
    query_seq = str(query_record.seq)
    target_seq = str(target_record.seq)
    if not query_seq or not target_seq:
        return ({"error": "query or target sequence is empty"}, 422)
    # Length-cap before kicking off the C-loop alignment. Even though
    # `_pairwise_align` enforces the same cap, surfacing it as 413
    # here gives a clearer error than a generic "rejected" downstream.
    if len(query_seq) > _PAIRWISE_MAX_LEN:
        return ({"error": f"query exceeds {_PAIRWISE_MAX_LEN:,} bp"},
                413)
    if len(target_seq) > _PAIRWISE_MAX_LEN:
        return ({"error": f"target exceeds {_PAIRWISE_MAX_LEN:,} bp"},
                413)
    # Circular rotation: same auto-detect as `diff-plasmid`.
    circ_raw = payload.get("circular")
    if circ_raw is None:
        target_annotations = getattr(target_record, "annotations",
                                       None) or {}
        is_circular = (target_annotations.get("topology") == "circular")
    else:
        is_circular = bool(circ_raw)
    # INV-72 (2026-05-25): use the same picker the UI uses. Pre-sweep
    # this endpoint ran bare `_pairwise_align` after a single
    # `_find_circular_alignment_offset` call — agents missed the
    # RC-orientation + multi-rotation best-of-N pick that `_align_worker`
    # gained in `[INV-71]`. canvas_axis defaults to "target" (matches
    # the UI: the read is the query, the library entry is the target).
    try:
        result = _state._pick_best_rotation_hook(
            query_seq, target_seq,
            is_circular=is_circular,
            mode=mode, canvas_axis="target",
        )
    except ValueError as exc:
        return ({"error": f"alignment rejected: {exc}"}, 400)
    except Exception as exc:
        _log.exception("align-plasmidsaurus-zip: pick_best_rotation raised")
        return ({"error": f"alignment failed: {_scrub_path(str(exc))}"}, 500)
    # Sweep #26 (2026-05-23): post-alignment drift check. If the
    # target entry was deleted between resolve and alignment
    # completion (multi-second CPU-bound `_pairwise_align`), surface
    # 410 Gone — the alignment result is technically valid but the
    # named target no longer exists, so any agent follow-up against
    # it (set-active, load-entry) would 404.
    if resolved_target_id:
        current = _find_library_entry_by_id(resolved_target_id)
        if current is None:
            return ({"error": (
                "target deleted mid-flight; alignment ran against "
                "the resolved snapshot but the library entry is gone"
            ), "target_id": resolved_target_id}, 410)
        # If rename happened, the result still ships, but flag the
        # drift in the payload so the agent's follow-ups can use the
        # current name.
        current_name = current.get("name") or ""
        if current_name and current_name != (target_record.name
                                              or target_record.id or ""):
            result = dict(result)
            result["_target_renamed_to"] = current_name
    picked = result.get("picked_rotation", "none")
    target_rotation = int(result.get("target_rotation") or 0)
    query_rotation = int(result.get("query_rotation") or 0)
    _log_event(
        "alignment.agent",
        path=str(path), member=member,
        target=target_record.name or target_record.id or "",
        identity_pct=round(float(result.get("identity_pct") or 0), 1),
        picked_rotation=picked,
        query_rotation=query_rotation,
        target_rotation=target_rotation,
        query_rc=bool(result.get("query_rc", False)),
        mode=mode,
    )
    return {
        "ok":              True,
        "path":            str(path),
        "member":          member,
        "query_name":      query_record.name or query_record.id or "",
        "target_id":       target_entry.get("id") or "",
        "target_name":     target_record.name or target_record.id or "",
        "circular":        is_circular,
        # Back-compat: `rotation_offset` mirrors `target_rotation`
        # (the legacy single-rotation API only ever rotated the target).
        # New `picked_rotation` / `query_rotation` / `query_rc` fields
        # carry the full picker decision so agents can replicate the
        # alignment.
        "rotation_offset": target_rotation,
        "picked_rotation": picked,
        "query_rotation":  query_rotation,
        "target_rotation": target_rotation,
        "query_rc":        bool(result.get("query_rc", False)),
        "result":          result,
    }


@_agent_endpoint("blast")
def _h_blast(app, payload):
    """In-process BLASTN / BLASTP against the user's plasmid
    collections. Body: ``{query, program?, collections?, max_hits?,
    six_frame?, backend?}``. ``program`` is auto-detected from the
    query alphabet when omitted; ``collections`` defaults to all
    collections; ``max_hits`` defaults to 25 (capped at 500);
    ``six_frame`` enables BLASTP six-frame ORF indexing (off by
    default); ``backend`` is ``auto`` / ``pyhmmer`` / ``pure``.

    Snag #18: the in-process DB is built from the chosen collections, so
    a large WHOLE library can exceed the DB build cap and 400. Pass
    ``collections`` to scope the DB to just the collection(s) you're
    searching — that both avoids the cap and speeds the build.

    Read-only — never touches the loaded record or any saved file.
    Hits are returned as a list of dicts (subject name + collection,
    coords, score, identity %, alignment fragments)."""
    raw_query = payload.get("query")
    if not raw_query or not isinstance(raw_query, str):
        return ({"error": "missing or non-string 'query'"}, 400)
    program_hint = str(payload.get("program") or "").lower()
    if program_hint and program_hint not in ("blastn", "blastp"):
        return ({"error": "'program' must be 'blastn' or 'blastp'"}, 400)
    program, cleaned = _state._detect_query_program_hook(raw_query, program_hint)
    if not cleaned:
        return ({"error": "query is empty after sanitisation"}, 400)

    # Validate collections list — must be a JSON array of strings,
    # capped at a sane size so a malicious payload can't make us scan
    # 10,000 names.
    coll_raw = payload.get("collections")
    coll_names: "list[str] | None" = None
    if coll_raw is not None:
        if not isinstance(coll_raw, list):
            return ({"error": "'collections' must be a list of names"}, 400)
        if len(coll_raw) > 100:
            return ({"error": "'collections' too long (max 100)"}, 400)
        coll_names = []
        for n in coll_raw:
            norm = _normalize_collection_name(n)
            if norm is None:
                return ({"error":
                          f"invalid collection name in list: {n!r}"}, 400)
            coll_names.append(norm)

    max_hits = _coerce_int(payload.get("max_hits", 25),
                             name="max_hits")
    if isinstance(max_hits, str):
        return ({"error": max_hits}, 400)
    max_hits = max(1, min(max_hits, 500))

    six_frame = bool(payload.get("six_frame", False))
    backend = str(payload.get("backend") or "auto").lower()
    if backend not in ("auto", "pyhmmer", "pure"):
        return ({"error": "'backend' must be 'auto' / 'pyhmmer' / 'pure'"},
                400)

    try:
        db = _state._blast_get_db_hook(program, coll_names, six_frame=six_frame)
        hits = _state._blast_search_hook(cleaned, db, max_hits=max_hits, backend=backend)
    except Exception as exc:
        _log.exception("agent-api blast failed")
        return ({"error": f"blast failed: {exc}",
                 "type": type(exc).__name__}, 500)

    return {
        "ok":           True,
        "program":      program,
        "backend":      backend,
        "query_length": len(cleaned),
        "n_hits":       len(hits),
        "hits":         hits,
    }


@_agent_endpoint("hmmscan")
def _h_hmmscan(app, payload):
    """Scan a query AA sequence against a HMMER 3 profile DB. Body:
    ``{query, hmm_path, max_hits?}``. ``hmm_path`` must point at a
    readable .hmm / .h3m / .h3p file on the server filesystem; the
    file is read lazily so Pfam-scale DBs don't pre-fetch into RAM.
    Returns hits in the same shape as ``blast``.

    Path safety (audit sweep #4 2026-05-15): `hmm_path` is routed
    through `_safe_file_size_check` (2 GB cap, symlink rejection,
    `S_ISREG` check) before opening — without this an agent caller
    could pass `/dev/zero` or a symlink to it and DoS the worker
    via `pyhmmer.HMMFile`'s streaming read."""
    raw_query = payload.get("query")
    if not raw_query or not isinstance(raw_query, str):
        return ({"error": "missing or non-string 'query'"}, 400)
    hmm_path = _sanitize_path(payload.get("hmm_path"))
    if hmm_path is None:
        return ({"error": "missing 'hmm_path'"}, 400)
    # Sweep #11 (2026-05-20): collapse "file-not-acceptable" responses
    # to a generic 400 so an unauthenticated caller can't probe
    # filesystem state via the error differential (pre-fix the
    # responses for "not found" / "not a regular file" / "too large"
    # were distinguishable, letting a local attacker enumerate which
    # paths exist + what kind). The detail still lands in the log
    # for the SpliceCraft user's own diagnostic bundle.
    _HMMSCAN_HMM_MAX_BYTES = 2 * 1024 * 1024 * 1024  # 2 GB
    _generic_rejection = (
        {"error": "hmm file not acceptable (see splicecraft log)"}, 400,
    )
    # Ancestor-symlink parity (audit 2026-06-30): the size/symlink check below
    # covers the target itself, but a PARENT being a symlink (Documents → /etc)
    # would still redirect the open. Same generic rejection — no probe leak.
    if _check_agent_read_path_ancestors(hmm_path) is not None:
        _log.info("hmmscan: refused path (ancestor symlink): %s", hmm_path)
        return _generic_rejection
    if not hmm_path.exists():
        _log.info("hmmscan: refused path (not found): %s", hmm_path)
        return _generic_rejection
    # Symlink + size + regular-file check. Asymmetric vs
    # `_h_load_file` which already routes through `_safe_file_size_check`;
    # without this, `pyhmmer.HMMFile(/dev/zero)` streams forever and an
    # agent caller can DoS the worker. Cap at 2 GB — real Pfam-A.hmm
    # bundles weigh in around 1.3 GB so this is generous but bounded.
    ok, reason = _safe_file_size_check(
        hmm_path, _HMMSCAN_HMM_MAX_BYTES, "hmm file",
    )
    if not ok:
        _log.info(
            "hmmscan: refused path %s — %s", hmm_path, reason,
        )
        return _generic_rejection
    max_hits = _coerce_int(payload.get("max_hits", 25),
                             name="max_hits")
    if isinstance(max_hits, str):
        return ({"error": max_hits}, 400)
    max_hits = max(1, min(max_hits, 500))
    try:
        hits = _state._hmmscan_run_hook(raw_query, str(hmm_path), max_hits=max_hits)
    except FileNotFoundError as exc:
        # Sweep #26: collapse the TOCTOU race-window 404 to a uniform
        # log-only message. The sibling `_safe_file_size_check` /
        # path-shape checks earlier in the handler already returned
        # the same opaque string; matching the wording here closes
        # the path-existence oracle that an attacker could probe
        # via "file removed between exists() and open()" races.
        _log.warning("agent hmmscan: %s", _scrub_path(str(exc)))
        return ({"error": "hmm file not acceptable (see splicecraft log)"},
                404)
    except ValueError as exc:
        return ({"error": _scrub_path(str(exc))}, 400)
    except Exception as exc:
        _log.exception("agent-api hmmscan failed")
        return ({"error": f"hmmscan failed: {_scrub_path(str(exc))}",
                 "type": type(exc).__name__}, 500)
    return {"ok": True, "n_hits": len(hits), "hits": hits}


# ── Online search (NCBI BLAST URL-API + EBI HMMER web) ─────────────────────
# Off by default. The TWO-KEY egress gate: (1) the human arms
# `allow_online_search` (Settings → "Allow agent online BLAST/HMMER" — the
# setting is DELIBERATELY absent from `_AGENT_SETTINGS_ALLOWLIST`, so an
# autonomous agent cannot flip it on itself), AND (2) an agent calls the
# endpoint by name. Even when armed, the engines route through the hardened
# opener (`_online_http` → `_build_hardened_url_opener`) whose egress is
# fail-closed-gated by `_state._demo_block_network_hook`, so the web demo
# (DEMO_MODE) stays offline regardless. Both are HEAVY (a remote job can
# block several minutes) and so sit in `_AGENT_HEAVY_ENDPOINTS` (hub) to
# bound concurrency. The in-process `blast` / `hmmscan` never consult the
# gate — they keep everything local.
_ONLINE_SEARCH_MAX_HITS = 100          # polite cap for a public service
_VALID_ONLINE_BLAST_PROGRAMS = frozenset(
    {"blastn", "blastp", "blastx", "tblastn", "tblastx"})


def _online_search_armed() -> bool:
    """The human half of the two-key gate: True only when the user has armed
    `allow_online_search` (Settings → "Allow agent online BLAST/HMMER", or a
    hand-edited settings.json). Read live at call time so a fresh toggle takes
    effect without a restart."""
    return bool(_get_setting("allow_online_search", False))


def _online_search_refusal():
    """Uniform 403 for a disarmed online-search call, pointing at the human-
    only arming path so the caller knows it can't self-enable."""
    return ({"error":
              "online search is disarmed — an agent may not send sequences "
              "off-box. The SpliceCraft user must enable Settings → \"Allow "
              "agent online BLAST/HMMER\" (or set allow_online_search in "
              "settings.json) first. The in-process `blast` / `hmmscan` "
              "endpoints stay available and never leave the machine."}, 403)


# ── Online reference-database lookups (Babs / agent-API) ───────────────────
# A SEPARATE egress gate from `allow_online_search` above: these lookups send
# only a short QUERY STRING to public databases (FPbase, UniProt, Europe PMC,
# NCBI, Wikipedia, patents, the web) — never the user's sequence. Armed by the
# human-only `allow_online_lookups` setting (NOT in `_AGENT_SETTINGS_ALLOWLIST`,
# so an autonomous agent can't self-arm). Read-only, but egress → arm-required.
def _online_lookups_armed() -> bool:
    """True only when the user has armed `allow_online_lookups` (Settings →
    "Allow Babs online database lookups", or a hand-edited settings.json).
    Read live so a fresh toggle takes effect without a restart."""
    return bool(_get_setting("allow_online_lookups", False))


def _online_lookups_refusal():
    """Uniform 403 for a disarmed online-lookup call, telling the caller the
    human-only path to arm it (so Babs can relay it to the user)."""
    return ({"error":
              "online lookups are disarmed — Babs may not query external "
              "databases. The SpliceCraft user must enable Settings → \"Allow "
              "Babs online database lookups\" (or set allow_online_lookups in "
              "settings.json) first. Only a short query string is ever sent; "
              "your sequences never leave the machine."}, 403)


def _validate_lookup_payload(payload, cap):
    """Validate the shared ``{query, max_hits?}`` body. Returns
    ``(query, max_hits, None)`` on success or ``(None, 0, (err, status))`` so a
    caller can ``if err: return err``. `cap` bounds max_hits per service."""
    q = payload.get("query")
    if not isinstance(q, str) or not q.strip():
        return None, 0, ({"error": "missing or non-string 'query'"}, 400)
    q = q.strip()
    if len(q) > _ONLINE_LOOKUP_QUERY_MAX:
        return None, 0, ({"error":
                          f"'query' too long (max "
                          f"{_ONLINE_LOOKUP_QUERY_MAX} chars)"}, 400)
    n = _coerce_int(payload.get("max_hits", 10), name="max_hits")
    if isinstance(n, str):
        return None, 0, ({"error": n}, 400)
    return q, max(1, min(n, cap)), None


def _online_lookup_error(exc, label):
    """Map an engine exception to an agent-API error tuple. Network / upstream
    failures (RuntimeError, incl. the rate-limit hints) → 502; anything else is
    logged and 500'd. Paths scrubbed either way."""
    if isinstance(exc, RuntimeError):
        return ({"error": f"{label} lookup failed: {_scrub_path(str(exc))}"}, 502)
    _log.exception("agent-api %s failed", label)
    return ({"error": f"{label} lookup failed: {_scrub_path(str(exc))}",
             "type": type(exc).__name__}, 500)


@_agent_endpoint("blast-online")
def _h_blast_online(app, payload):
    """Remote NCBI BLAST via the public URL API. Body:
    ``{query, program?, database?, max_hits?}``.

    EGRESS WARNING: this SHIPS your query sequence to NCBI's servers
    (unlike the in-process `blast`, which never leaves the box). Refused
    403 unless the SpliceCraft user has armed Settings → "Allow agent
    online BLAST/HMMER" — an autonomous agent cannot enable it itself
    (the setting is not in the agent allowlist).

    ``program`` is one of blastn / blastp / blastx / tblastn / tblastx
    (auto-detected as blastn or blastp from the query alphabet when
    omitted; the translated programs need an explicit ask); ``database``
    defaults to NCBI `nt` (nucleotide programs) or `nr` (protein) and
    may be any short NCBI db name; ``max_hits`` defaults to 25, capped
    at 100. Blocks until NCBI returns (up to a few minutes) or errors —
    502 on a network / server failure. Read-only; never touches the
    loaded record or any saved file. Hits: accession, description,
    identity %, e-value, bit score, query/subject coords."""
    if not _online_search_armed():
        return _online_search_refusal()
    raw_query = payload.get("query")
    if not raw_query or not isinstance(raw_query, str):
        return ({"error": "missing or non-string 'query'"}, 400)
    program_hint = str(payload.get("program") or "").lower().strip()
    if program_hint:
        if program_hint not in _VALID_ONLINE_BLAST_PROGRAMS:
            return ({"error":
                      f"'program' must be one of "
                      f"{sorted(_VALID_ONLINE_BLAST_PROGRAMS)}"}, 400)
        program = program_hint
    else:
        # Auto-detect blastn/blastp from the alphabet (the in-process
        # detector); the three translated programs need an explicit ask.
        program, _unused = _state._detect_query_program_hook(raw_query, "")
    cleaned = _online_clean_query(raw_query, program)
    if not cleaned:
        return ({"error": "query is empty after sanitisation"}, 400)
    cap = _online_max_query_len(program)
    if len(cleaned) > cap:
        return ({"error":
                  f"query is {len(cleaned):,} chars; NCBI's {program} limit "
                  f"is {cap:,} — search a shorter region"}, 413)
    database = payload.get("database")
    if database in (None, ""):
        database = _ncbi_blast_db_for(program)
    elif not isinstance(database, str):
        return ({"error": "'database' must be a string"}, 400)
    else:
        database = database.strip()
        # NCBI db NAMES (not paths): letters/digits/_.- only. No slash and no
        # ".." so a hostile `database` can't smuggle a path-ish token into the
        # outbound DATABASE param.
        if (not database or len(database) > 64 or ".." in database
                or not re.fullmatch(r"[A-Za-z0-9_.-]+", database)):
            return ({"error":
                      "'database' must be a short NCBI db name "
                      "(e.g. nt, nr, refseq_rna)"}, 400)
    max_hits = _coerce_int(payload.get("max_hits", 25), name="max_hits")
    if isinstance(max_hits, str):
        return ({"error": max_hits}, 400)
    max_hits = max(1, min(max_hits, _ONLINE_SEARCH_MAX_HITS))
    try:
        hits = _ncbi_blast_online(cleaned, program, database, max_hits)
    except RuntimeError as exc:
        # Network / NCBI-server error (timeout, rejected job, 5xx). Upstream
        # failure, not a caller bug → 502.
        return ({"error": f"NCBI BLAST failed: {_scrub_path(str(exc))}"}, 502)
    except Exception as exc:
        _log.exception("agent-api blast-online failed")
        return ({"error": f"blast-online failed: {_scrub_path(str(exc))}",
                 "type": type(exc).__name__}, 500)
    return {"ok": True, "program": program, "database": database,
            "query_length": len(cleaned), "n_hits": len(hits), "hits": hits}


@_agent_endpoint("hmmer-web")
def _h_hmmer_web(app, payload):
    """Remote hmmscan vs Pfam via the EBI HMMER web API. Body:
    ``{query, max_hits?}`` — `query` is a PROTEIN sequence.

    EGRESS WARNING: this SHIPS your protein query to EBI's servers
    (unlike the in-process `hmmscan` against a locally-downloaded Pfam
    database). Refused 403 unless the SpliceCraft user has armed
    Settings → "Allow agent online BLAST/HMMER"; an autonomous agent
    cannot enable it itself.

    ``max_hits`` defaults to 25, capped at 100. Blocks until EBI
    returns (up to a few minutes) or errors — 502 on a network / server
    failure. Read-only. Hits: Pfam accession + family name +
    description, e-value, bit score, # domains, clan."""
    if not _online_search_armed():
        return _online_search_refusal()
    raw_query = payload.get("query")
    if not raw_query or not isinstance(raw_query, str):
        return ({"error": "missing or non-string 'query'"}, 400)
    protein = _online_clean_query(raw_query, "hmmscan")
    if not protein:
        return ({"error": "query is empty after sanitisation"}, 400)
    cap = _online_max_query_len("hmmscan")
    if len(protein) > cap:
        return ({"error":
                  f"query is {len(protein):,} aa; the cap is {cap:,}"}, 413)
    max_hits = _coerce_int(payload.get("max_hits", 25), name="max_hits")
    if isinstance(max_hits, str):
        return ({"error": max_hits}, 400)
    max_hits = max(1, min(max_hits, _ONLINE_SEARCH_MAX_HITS))
    try:
        hits = _hmmer_web_hmmscan(protein, max_hits)
    except RuntimeError as exc:
        return ({"error": f"EBI HMMER failed: {_scrub_path(str(exc))}"}, 502)
    except Exception as exc:
        _log.exception("agent-api hmmer-web failed")
        return ({"error": f"hmmer-web failed: {_scrub_path(str(exc))}",
                 "type": type(exc).__name__}, 500)
    return {"ok": True, "query_length": len(protein),
            "n_hits": len(hits), "hits": hits}


@_agent_endpoint("fpbase-search")
def _h_fpbase_search(app, payload):
    """Search FPbase (fpbase.org) for FLUORESCENT PROTEINS by name. Body:
    ``{query, max_hits?}`` (`query` is a substring, e.g. "mCherry", "GFP").

    Requires the user to arm Settings → "Allow Babs online database lookups"
    (403 otherwise); only the query string is sent — never your sequence.
    Read-only. Each hit: name, spectra (ex/em maxima, quantum yield, extinction
    coefficient, brightness), oligomerization, switch type, GenBank/UniProt
    xrefs, the amino-acid sequence, and the FPbase URL."""
    if not _online_lookups_armed():
        return _online_lookups_refusal()
    query, n, err = _validate_lookup_payload(payload, _ONLINE_LOOKUP_MAX_HITS)
    if err:
        return err
    assert query is not None  # validated non-None above; narrows for the engine
    try:
        results = _fpbase_search(query, n)
    except Exception as exc:
        return _online_lookup_error(exc, "fpbase-search")
    return {"ok": True, "source": "FPbase", "query": query,
            "count": len(results), "results": results}


@_agent_endpoint("uniprot-search")
def _h_uniprot_search(app, payload):
    """Search UniProtKB for PROTEINS. Body: ``{query, max_hits?}`` (`query` is
    a UniProt query, e.g. "Cas9", "NPTII", "gene:BAR AND organism_id:3702").

    Requires arming Settings → "Allow Babs online database lookups" (403
    otherwise); only the query string leaves the box. Read-only. Each hit:
    accession, protein name, organism, the curated FUNCTION annotation,
    keywords, reviewed (Swiss-Prot) flag, and the UniProt URL."""
    if not _online_lookups_armed():
        return _online_lookups_refusal()
    query, n, err = _validate_lookup_payload(payload, _ONLINE_LOOKUP_MAX_HITS)
    if err:
        return err
    assert query is not None  # validated non-None above; narrows for the engine
    try:
        results = _uniprot_search(query, n)
    except Exception as exc:
        return _online_lookup_error(exc, "uniprot-search")
    return {"ok": True, "source": "UniProt", "query": query,
            "count": len(results), "results": results}


@_agent_endpoint("literature-search")
def _h_literature_search(app, payload):
    """Search the scientific literature via Europe PMC (PubMed + PMC +
    preprints). Body: ``{query, max_hits?}``.

    Requires arming Settings → "Allow Babs online database lookups" (403
    otherwise); only the query string leaves the box. Read-only. Each hit:
    title, authors, journal, year, DOI/PMID/PMCID, a truncated abstract,
    open-access flag, and a resolvable URL."""
    if not _online_lookups_armed():
        return _online_lookups_refusal()
    query, n, err = _validate_lookup_payload(payload, _ONLINE_LOOKUP_MAX_HITS)
    if err:
        return err
    assert query is not None  # validated non-None above; narrows for the engine
    try:
        results = _europepmc_search(query, n)
    except Exception as exc:
        return _online_lookup_error(exc, "literature-search")
    return {"ok": True, "source": "Europe PMC", "query": query,
            "count": len(results), "results": results}


@_agent_endpoint("genbank-search")
def _h_genbank_search(app, payload):
    """Search NCBI for SEQUENCE RECORDS by term. Body:
    ``{query, db?, max_hits?}`` — `db` is "nucleotide" (default) or "protein".

    Requires arming Settings → "Allow Babs online database lookups" (403
    otherwise); only the query string leaves the box. Read-only. Returns
    `total_matches` plus records (accession, title, organism, length, URL) —
    feed a hit's accession to the `fetch` endpoint to load the full record."""
    if not _online_lookups_armed():
        return _online_lookups_refusal()
    query, n, err = _validate_lookup_payload(payload, _ONLINE_LOOKUP_MAX_HITS)
    if err:
        return err
    assert query is not None  # validated non-None above; narrows for the engine
    db = str(payload.get("db") or "nucleotide").lower().strip()
    if db == "nuccore":
        db = "nucleotide"
    if db not in ("nucleotide", "protein"):
        return ({"error": "'db' must be 'nucleotide' or 'protein'"}, 400)
    try:
        results, total = _ncbi_db_search(query, db, n)
    except Exception as exc:
        return _online_lookup_error(exc, "genbank-search")
    return {"ok": True, "source": "NCBI", "db": db, "query": query,
            "total_matches": total, "count": len(results), "results": results}


@_agent_endpoint("wikipedia-search")
def _h_wikipedia_search(app, payload):
    """Search Wikipedia for encyclopedic background. Body:
    ``{query, max_hits?}``.

    Requires arming Settings → "Allow Babs online database lookups" (403
    otherwise); only the query string leaves the box. Read-only. Each hit:
    title, a lead-section summary (or search snippet), and the article URL."""
    if not _online_lookups_armed():
        return _online_lookups_refusal()
    query, n, err = _validate_lookup_payload(payload, _ONLINE_LOOKUP_MAX_HITS)
    if err:
        return err
    assert query is not None  # validated non-None above; narrows for the engine
    try:
        results = _wikipedia_search(query, n)
    except Exception as exc:
        return _online_lookup_error(exc, "wikipedia-search")
    return {"ok": True, "source": "Wikipedia", "query": query,
            "count": len(results), "results": results}


@_agent_endpoint("web-search")
def _h_web_search(app, payload):
    """General web search. Body: ``{query, max_hits?}``.

    Requires arming Settings → "Allow Babs online database lookups" (403
    otherwise); only the query string leaves the box. Read-only. Uses the
    official Brave Search API when a `brave_search_api_key` is set in Settings
    (robust); otherwise falls back to DuckDuckGo's keyless HTML endpoint, which
    is BEST-EFFORT and rate-limited (a 502 with a "rate-limited" hint if
    blocked). Each hit: title, url, snippet."""
    if not _online_lookups_armed():
        return _online_lookups_refusal()
    query, n, err = _validate_lookup_payload(payload, _ONLINE_LOOKUP_MAX_HITS)
    if err:
        return err
    assert query is not None  # validated non-None above; narrows for the engine
    brave_key = (_get_setting("brave_search_api_key", "") or "").strip()
    try:
        results = _web_search(query, n, brave_key)
    except Exception as exc:
        return _online_lookup_error(exc, "web-search")
    return {"ok": True, "source": "web", "query": query,
            "provider": "Brave" if brave_key else "DuckDuckGo (best-effort)",
            "count": len(results), "results": results}


def _url_host(url) -> str:
    """Just the hostname of a URL, for privacy-preserving egress logging — no
    path or query string (a query can carry sensitive search terms)."""
    try:
        from urllib.parse import urlsplit
        return urlsplit(str(url)).hostname or "?"
    except (ValueError, TypeError):
        return "?"


@_agent_endpoint("read-url")
def _h_read_url(app, payload):
    """Open ONE web page and return its readable text (HTML → text; NO
    JavaScript). Body: ``{url, max_chars?}``.

    Pairs with web-search: open a result URL and actually READ the page, not
    just its snippet. Requires arming Settings → "Allow Babs online database
    lookups" (403 otherwise); only the URL is sent. Read-only, egress.
    SSRF-hardened — public http(s) hosts only (private / loopback / link-local
    refused) on the initial URL AND every redirect hop, the response size is
    capped, and page scripts are never run. Binary pages (PDF, images,
    downloads) are refused with a 502. ``max_chars`` bounds the returned text
    (default 20000, hard max 100000). Returns ``{url (final, post-redirect),
    title, content_type, text, chars, truncated}``."""
    if not _online_lookups_armed():
        return _online_lookups_refusal()
    url = payload.get("url")
    if not isinstance(url, str) or not url.strip():
        return {"error": "missing or non-string 'url'"}, 400
    mc = payload.get("max_chars")
    if mc is not None:
        mc = _coerce_int(mc, name="max_chars")
        if isinstance(mc, str):
            return {"error": mc}, 400
    try:
        data = _read_url(url) if mc is None else _read_url(url, mc)
    except Exception as exc:
        # Audit arbitrary-URL egress (host only — no query string) so a failed
        # fetch is diagnosable from the log the user attaches to a bug report.
        _log_event("agent.read_url.fail", host=_url_host(url),
                   reason=_scrub_path(str(exc))[:200], via="agent")
        return _online_lookup_error(exc, "read-url")
    _log_event("agent.read_url.ok", host=_url_host(data.get("url") or url),
               content_type=data.get("content_type"), chars=data.get("chars"),
               truncated=data.get("truncated"), via="agent")
    return {"ok": True, "source": "web", **data}


@_agent_endpoint("patent-search")
def _h_patent_search(app, payload):
    """Search patents. Body: ``{query, max_hits?}``.

    Requires arming Settings → "Allow Babs online database lookups" (403
    otherwise); only the query string leaves the box. Read-only. Uses the
    official PatentsView API when a `patentsview_api_key` is set in Settings
    (robust, US patents); otherwise falls back to Google Patents' keyless XHR
    endpoint, which is BEST-EFFORT and rate-limited (a 502 with a "rate-limited"
    hint if blocked). Returns `total_matches` plus records (number, title,
    snippet, assignee, dates, Google Patents URL)."""
    if not _online_lookups_armed():
        return _online_lookups_refusal()
    query, n, err = _validate_lookup_payload(payload, _ONLINE_LOOKUP_MAX_HITS)
    if err:
        return err
    assert query is not None  # validated non-None above; narrows for the engine
    pv_key = (_get_setting("patentsview_api_key", "") or "").strip()
    try:
        results, total = _patent_search(query, n, pv_key)
    except Exception as exc:
        return _online_lookup_error(exc, "patent-search")
    return {"ok": True, "source": "patents", "query": query,
            "provider": "PatentsView" if pv_key else "Google Patents (best-effort)",
            "total_matches": total, "count": len(results), "results": results}


@_agent_endpoint("list-parts-bins")
def _h_list_parts_bins(app, payload):
    """Every named parts bin. Each item: ``{name, n_parts, description}``;
    `active` is the currently-selected bin name (empty if none)."""
    out = []
    for b in _load_parts_bin_collections():
        out.append({
            "name":        b.get("name", "?"),
            "n_parts":     len(b.get("parts", []) or []),
            "description": str(b.get("description") or "")[:200],
        })
    return {"ok": True, "parts_bins": out,
            "active": _get_active_parts_bin_name() or ""}


@_agent_endpoint("set-active-parts-bin", write=True)
def _h_set_active_parts_bin(app, payload):
    """Switch the active parts bin. Body: ``{name}`` (one of
    `list-parts-bins`). The named bin in `parts_bin_collections.json` is
    the source of truth; its parts are mirrored into the live
    `parts_bin.json` via `_safe_save_json_mirror` so the switch may
    legitimately shrink the mirror without tripping the L3 guard
    ([INV-83]). RMW under `_state._cache_lock`, matching `PartsBinPickerModal`."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing 'name'"}, 400)
    name = name.strip()
    with _state._cache_lock:
        target = _find_parts_bin(name)
        if target is None:
            valid = [b.get("name") for b in _load_parts_bin_collections()]
            return ({"error":
                      f"unknown parts bin {name!r}; "
                      f"valid: {sorted(v for v in valid if v)}"}, 404)
        prev = _get_active_parts_bin_name()
        _set_active_parts_bin_name(name)
        _state._settings_flush_sync_hook()
        raw_parts = target.get("parts", [])
        if not isinstance(raw_parts, list):
            raw_parts = []
        target_parts = [p for p in raw_parts if isinstance(p, dict)]
        err = _agent_save_or_500(
            lambda: _safe_save_json_mirror(
                _state._PARTS_BIN_FILE, target_parts, "Parts bin"),
            "Parts bin")
        if err:
            # Mirror failed AFTER the pointer moved — revert so the active name
            # and the live parts_bin.json can't disagree (the desync would
            # clobber the target bin on the next save). [INV-161]
            _set_active_parts_bin_name(prev)
            _state._settings_flush_sync_hook()
            return err
        _state._parts_bin_cache = None
        _state._clear_assembly_fragment_cache_hook()
    _log_event("parts_bin.switch", name=name, via="agent")
    return {"ok": True, "active": name, "n_parts": len(target_parts)}


@_agent_endpoint("set-active-hmm-database", write=True)
def _h_set_active_hmm_database(app, payload):
    """Pick the active HMM database. Body: ``{id}`` (one of
    `list-hmm-databases`). Persists the `hmm_db_active_id` setting."""
    eid = payload.get("id")
    if not isinstance(eid, str) or not eid.strip():
        return ({"error": "missing 'id'"}, 400)
    eid = eid.strip()
    if eid not in {e.get("id") for e in _load_hmm_db_catalog()}:
        return ({"error": "hmm database not found"}, 404)
    _set_setting("hmm_db_active_id", eid)
    _state._settings_flush_sync_hook()
    _log_event("hmm_db.set_active", id=eid, via="agent")
    return {"ok": True, "active": eid}


@_agent_endpoint("list-backups")
def _h_list_backups(app, payload):
    """List recoverable backups for user-data files. Body: ``{label?}``.

    WITH ``label`` (one of plasmid_library / collections / parts_bin /
    parts_bin_collections / primers / primer_collections / features /
    feature_colors / grammars / entry_vectors / codon_tables / settings
    / experiments / experiment_projects / gels / protein_motifs /
    custom_enzymes / enzyme_collections / hmm_databases) returns the
    individual backup rows from the four-layer JSON safety net:
    `legacy_bak` (`<file>.bak`), `rotating_bak` (timestamped),
    `snapshot` (daily), `lost_entries` (shrink-guard spillover). Newest
    first.

    WITHOUT ``label`` (snag #6) returns a per-label SUMMARY across every
    user-data file — ``{labels: [{label, target_path, count, newest},
    ...], total}`` — so this "list" verb no longer requires an argument.
    Drill into a specific label for the full rows.
    """
    raw_label = payload.get("label")
    if raw_label in (None, ""):
        summary = []
        for lbl in sorted(_AGENT_BACKUP_LABELS):
            path, _msg = _resolve_backup_label(lbl)
            if path is None:
                continue
            rows = _list_recoverable_backups(path)
            summary.append({
                "label":       lbl,
                "target_path": str(path),
                "count":       len(rows),
                # Rows are newest-first, so [0] is the most recent.
                "newest":      (str(rows[0].get("mtime_str") or "")
                                if rows else ""),
            })
        return {"ok": True, "labels": summary,
                "total": sum(row["count"] for row in summary)}
    path, msg = _resolve_backup_label(raw_label)
    if path is None:
        return ({"error": msg}, 400)
    rows = []
    for b in _list_recoverable_backups(path):
        rows.append({
            "kind":        b.get("kind", ""),
            "source_path": str(b.get("source_path", "")),
            "n_entries":   int(b.get("n_entries") or 0),
            "mtime_str":   str(b.get("mtime_str") or ""),
        })
    return {"ok": True, "label": msg, "target_path": str(path),
            "backups": rows, "count": len(rows)}


@_agent_endpoint("restore-backup", write=True)
def _h_restore_backup(app, payload):
    """Restore one user-data file from a chosen backup. Body:
    ``{label, source_path}`` where ``source_path`` must come from a
    prior `list-backups` response (verified to belong to one of the
    four backup tiers for the same target).

    Triggers a fresh rotating backup of the current state as part of
    `_safe_save_json`'s atomic write, so the pre-restore data is
    NOT lost. Caches for the restored label are busted afterwards.
    """
    guard = _agent_dirty_guard(app, payload)
    if guard is not None:
        return guard
    path, msg = _resolve_backup_label(payload.get("label"))
    if path is None:
        return ({"error": msg}, 400)
    src_raw = payload.get("source_path")
    if not isinstance(src_raw, str) or not src_raw:
        return ({"error": "missing 'source_path'"}, 400)
    # Verify the source_path matches one of the registered backup
    # rows — otherwise an agent could read or write arbitrary files
    # under the user's uid via this endpoint.
    src = Path(src_raw).expanduser()
    legitimate = {Path(b["source_path"])
                  for b in _list_recoverable_backups(path)}
    if src not in legitimate:
        return ({"error": "source_path is not a registered backup for "
                  f"label {msg!r}"}, 403)
    try:
        n = _restore_from_backup(path, src, label=msg)
    except ValueError as exc:
        return ({"error": f"unreadable backup: {exc}"}, 422)
    except OSError as exc:
        _log.exception("agent restore-backup: write failed")
        return ({"error": f"restore write failed: {_scrub_path(str(exc))}"}, 500)
    # Bust the in-memory cache for the restored label so the next
    # read picks up the restored content.
    cache_attr = {
        "plasmid_library":       "_library_cache",
        "collections":           "_collections_cache",
        "parts_bin":             "_parts_bin_cache",
        "parts_bin_collections": "_parts_bin_collections_cache",
        "primers":               "_primers_cache",
        "primer_collections":    "_primer_collections_cache",
        "features":              "_features_cache",
        "feature_colors":        "_feature_colors_cache",
        "grammars":              "_grammars_cache",
        "entry_vectors":         "_entry_vectors_cache",
        "codon_tables":          "_codon_tables_cache",
        "settings":              "_settings_cache",
        "experiments":           "_experiments_cache",
        "experiment_projects":   "_experiment_projects_cache",
        "gels":                  "_gels_cache",
        "protein_motifs":        "_protein_motifs_cache",
        "custom_enzymes":        "_custom_enzymes_cache",
        "enzyme_collections":    "_enzyme_collections_cache",
        "hmm_db_catalog":        "_hmm_db_catalog_cache",
        "protein_collections":   "_protein_collections_cache",
        "model_collections":     "_model_collections_cache",
        "protocol_collections":  "_protocol_collections_cache",
        "custom_labware":        "_custom_labware_cache",
    }.get(msg)
    if cache_attr is not None:
        _state._reset_master_delete_cache_hook(cache_attr)
    # Mirror files (plasmid_library / primers / parts_bin) revert on the next
    # launch unless the restore is written through into the owning active
    # collection — reconcile via the hub hook (guarded for the import window).
    reconcile = getattr(_state, "_reconcile_mirror_after_restore_hook", None)
    if reconcile is not None:
        reconcile(path)
    return {"ok": True, "label": msg, "entries_restored": n,
            "source_path": str(src)}


@_agent_endpoint("restore-pre-update-snapshot", write=True)
def _h_restore_pre_update_snapshot(app, payload):
    """Restore a pre-update snapshot by id (or "latest"). Body:
    ``{id?: str}``. When ``id`` is missing or "latest", restores the
    newest snapshot.

    Same sacred four checks as the CLI path (schema version ≤
    current, attr in `_USER_DATA_FILE_ATTRS`, name regex
    `_PRE_UPDATE_NAME_RE`, SHA-256 re-verify before `os.replace`).
    """
    guard = _agent_dirty_guard(app, payload)
    if guard is not None:
        return guard
    raw_id = payload.get("id") or "latest"
    if not isinstance(raw_id, str):
        return ({"error": "'id' must be string"}, 400)
    raw_id = raw_id.strip() or "latest"
    if raw_id != "latest":
        # Enforce the regex on the wire boundary so a malformed id
        # never reaches `_restore_pre_update_snapshot`'s validator.
        if not _PRE_UPDATE_NAME_RE.match(raw_id):
            return ({"error": f"id {raw_id!r} fails the pre-update "
                      "snapshot name regex"}, 400)
    try:
        if raw_id == "latest":
            snaps = _list_pre_update_snapshots()
            if not snaps:
                return ({"error": "no pre-update snapshots available"},
                        404)
            target = snaps[0]["path"]
        else:
            target = raw_id
        report = _restore_pre_update_snapshot(target)
    except FileNotFoundError as exc:
        return ({"error": f"snapshot not found: {exc}"}, 404)
    except ValueError as exc:
        return ({"error": f"snapshot rejected: {exc}"}, 422)
    except OSError as exc:
        _log.exception("agent restore-pre-update-snapshot: write failed")
        return ({"error": f"restore write failed: {_scrub_path(str(exc))}"}, 500)
    # Bust every cached user-data so post-restore reads see the
    # restored state. Sweep #25 (2026-05-23) — drives the enumeration
    # from `_MASTER_DELETE_CACHE_ATTRS` so new caches enroll
    # automatically. Pre-sweep this list was hand-maintained and
    # drifted: sweep #20 added `_state._primer_usage_cache` +
    # `_state._feature_library_index_cache` (derived caches that must reset
    # when their source data changes); sweep #20 also added
    # `_state._primer_collections_cache`. None propagated here. Settings
    # cache IS reset here (unlike the UI restore path) because the
    # `splicecraft update --restore` flow is invoked outside the
    # running app session.
    for cache_attr in _MASTER_DELETE_CACHE_ATTRS:
        _state._reset_master_delete_cache_hook(cache_attr)
    return {"ok": True, "report": report}


@_agent_endpoint("simulate-pcr")
def _h_simulate_pcr(app, payload):
    """Run an in-silico PCR. Body:
    ``{template_seq, fwd_primer, rev_primer,
        circular?: bool = true,
        max_amplicon?: int = 20000}``.
    Natural-guess aliases accepted: ``template`` / ``sequence`` →
    ``template_seq``; ``forward`` / ``fwd`` → ``fwd_primer``;
    ``reverse`` / ``rev`` → ``rev_primer``.

    Binding model is exact-match (no mismatch tolerance, no Tm-aware
    annealing — see `_simulate_pcr` for the contract). Primers must be
    10–80 bp, ACGT only. Templates above ``_PCR_MAX_TEMPLATE_BP``
    (5 Mb) are refused rather than risking a chromosome-scale find.
    Returns up to ``_PCR_MAX_AMPLICONS`` (50) amplicons sorted by
    length descending; the ``capped`` field flags mispriming runaway.

    Read-only. To save an amplicon as a linear library entry, use the
    Simulator screen (which constructs a SeqRecord with primer_bind
    features at the ends) — the agent flow is meant for analysis.
    """
    # Natural-guess aliases (ergonomics nit) — canonical name wins when
    # present, else fall back. `or` is safe: an empty/missing canonical
    # is invalid anyway, so deferring to an alias only ever helps.
    template = (payload.get("template_seq")
                or payload.get("template") or payload.get("sequence"))
    if not isinstance(template, str):
        return ({"error": "missing or non-string 'template_seq' "
                  "(aliases: template / sequence)"}, 400)
    if len(template) > _PCR_MAX_TEMPLATE_BP:
        return ({"error": f"'template_seq' exceeds "
                  f"{_PCR_MAX_TEMPLATE_BP:,} bp PCR-sim cap"}, 413)
    fwd = (payload.get("fwd_primer")
           or payload.get("forward") or payload.get("fwd"))
    rev = (payload.get("rev_primer")
           or payload.get("reverse") or payload.get("rev"))
    if not isinstance(fwd, str) or not isinstance(rev, str):
        return ({"error": "missing or non-string 'fwd_primer' / "
                  "'rev_primer' (aliases: forward/fwd, reverse/rev)"}, 400)
    fwd_clean = fwd.upper().strip()
    rev_clean = rev.upper().strip()
    if not fwd_clean or not rev_clean:
        return ({"error": "'fwd_primer' and 'rev_primer' must be "
                  "non-empty"}, 400)
    if len(fwd_clean) < _PCR_MIN_PRIMER_LEN or \
            len(rev_clean) < _PCR_MIN_PRIMER_LEN:
        return ({"error": f"primers must be at least "
                  f"{_PCR_MIN_PRIMER_LEN} bp"}, 400)
    if len(fwd_clean) > _PCR_MAX_PRIMER_LEN or \
            len(rev_clean) > _PCR_MAX_PRIMER_LEN:
        return ({"error": f"primers must be at most "
                  f"{_PCR_MAX_PRIMER_LEN} bp"}, 400)
    if any(c not in "ACGT" for c in fwd_clean):
        return ({"error": "'fwd_primer' must be ACGT only "
                  "(IUPAC ambiguity not supported)"}, 400)
    if any(c not in "ACGT" for c in rev_clean):
        return ({"error": "'rev_primer' must be ACGT only "
                  "(IUPAC ambiguity not supported)"}, 400)
    max_amp = _coerce_int(payload.get("max_amplicon",
                                        _PCR_DEFAULT_MAX_AMPLICON),
                           name="max_amplicon")
    if isinstance(max_amp, str):
        return ({"error": max_amp}, 400)
    if not (1 <= max_amp <= _PCR_AMPLICON_HARD_CAP):
        return ({"error": f"'max_amplicon' must be in "
                  f"[1, {_PCR_AMPLICON_HARD_CAP}]"}, 400)
    circular = bool(payload.get("circular", True))
    try:
        amps = _state._simulate_pcr_hook(template, fwd_clean, rev_clean,
                              circular=circular, max_amplicon=max_amp)
    except Exception as exc:
        _log.exception("agent simulate-pcr: simulator failed")
        return ({"error": f"simulator failed: {_scrub_path(str(exc))}"}, 500)
    _log_event("simulator.pcr.agent",
                template_bp=len(template),
                circular=circular,
                fwd_len=len(fwd_clean), rev_len=len(rev_clean),
                max_amplicon=max_amp,
                n_amplicons=len(amps),
                capped=(len(amps) >= _PCR_MAX_AMPLICONS))
    return {
        "ok":         True,
        "n":          len(amps),
        "capped":     len(amps) >= _PCR_MAX_AMPLICONS,
        "amplicons":  amps,
    }


@_agent_endpoint("simulate-gel")
def _h_simulate_gel(app, payload):
    """Render an in-silico agarose gel. Body:
    ``{lanes:           [{name?, source, detail?}, ...],
        agarose_pct?:    float = 1.0,
        template_seq?:   str = "",
        template_circular?: bool = true,
        pcr_amplicon?:   dict | null = null,
        height?:         int = 22,
        lane_width?:     int = 7,
        include_image?:  bool = false}``.

    `source` is one of ``"empty" | "ladder" | "plasmid" | "digest" |
    "pcr"``. `detail` semantics depend on source: ladder → ladder
    name (``"1 kb Plus"``, ``"1 kb"``, ``"100 bp"``,
    ``"Lambda/HindIII"``); digest → comma-separated enzyme names;
    plasmid/pcr/empty → ignored. `pcr_amplicon` is the output of
    ``simulate-pcr`` (when source = ``"pcr"``).

    Returns structured band data per lane (bp / form / mobility /
    display row). When ``include_image=true``, also returns the
    rendered gel as a plain-text string (one line per row).

    Read-only.
    """
    lanes = payload.get("lanes")
    if not isinstance(lanes, list):
        return ({"error": "'lanes' must be a list"}, 400)
    if not lanes:
        return ({"error": "'lanes' must contain at least one lane"}, 400)
    if len(lanes) > _GEL_MAX_LANES:
        return ({"error": f"too many lanes (max {_GEL_MAX_LANES})"}, 400)
    cleaned_lanes: list[dict] = []
    for i, lane in enumerate(lanes):
        if not isinstance(lane, dict):
            return ({"error": f"lanes[{i}] must be a dict"}, 400)
        src = lane.get("source")
        if not isinstance(src, str) or not src:
            return ({"error": f"lanes[{i}].source missing or "
                      "non-string"}, 400)
        valid_sources = {kind for _, kind in _LANE_SOURCES}
        if src.lower() not in valid_sources:
            return ({"error": f"lanes[{i}].source {src!r} not in "
                      f"{sorted(valid_sources)}"}, 400)
        name = _sanitize_label(lane.get("name"), max_len=80) or f"L{i+1}"
        detail = lane.get("detail", "")
        if not isinstance(detail, str):
            return ({"error": f"lanes[{i}].detail must be a string"},
                    400)
        if len(detail) > 200:
            return ({"error": f"lanes[{i}].detail exceeds 200 chars"},
                    400)
        cleaned_lanes.append({
            "name":   name,
            "source": src.lower(),
            "detail": detail,
        })
    template = payload.get("template_seq", "")
    if not isinstance(template, str):
        return ({"error": "'template_seq' must be a string"}, 400)
    if len(template) > _PCR_MAX_TEMPLATE_BP:
        return ({"error": f"'template_seq' exceeds "
                  f"{_PCR_MAX_TEMPLATE_BP:,} bp cap"}, 413)
    template_circular = bool(payload.get("template_circular", True))
    pcr_amplicon = payload.get("pcr_amplicon")
    if pcr_amplicon is not None and not isinstance(pcr_amplicon, dict):
        return ({"error": "'pcr_amplicon' must be a dict or null"}, 400)
    agarose = payload.get("agarose_pct", 1.0)
    try:
        agarose = float(agarose)
    except (TypeError, ValueError):
        return ({"error": "'agarose_pct' must be a number"}, 400)
    # _agarose_mobility snaps to nearest configured %, but reject
    # absurdly out-of-range values up front so the agent gets a
    # clear error instead of silent snapping.
    if not (0.1 <= agarose <= 10.0):
        return ({"error": "'agarose_pct' must be in [0.1, 10.0]"}, 400)
    height = _coerce_int(payload.get("height", 22), name="height")
    if isinstance(height, str):
        return ({"error": height}, 400)
    if not (_GEL_HEIGHT_MIN <= height <= _GEL_HEIGHT_MAX):
        return ({"error": f"'height' must be in "
                  f"[{_GEL_HEIGHT_MIN}, {_GEL_HEIGHT_MAX}]"}, 400)
    lane_width = _coerce_int(payload.get("lane_width", 7),
                              name="lane_width")
    if isinstance(lane_width, str):
        return ({"error": lane_width}, 400)
    if not (_GEL_LANE_WIDTH_MIN <= lane_width <= _GEL_LANE_WIDTH_MAX):
        return ({"error": f"'lane_width' must be in "
                  f"[{_GEL_LANE_WIDTH_MIN}, "
                  f"{_GEL_LANE_WIDTH_MAX}]"}, 400)
    include_image = bool(payload.get("include_image", False))
    # Compute per-lane bands + mobility.
    lane_results: list[dict] = []
    try:
        for li, lane in enumerate(cleaned_lanes):
            bands = _gel_bands_for_lane(
                lane,
                template_seq=template,
                template_circular=template_circular,
                pcr_amplicon=pcr_amplicon,
            )
            band_info: list[dict] = []
            for bp, form in bands:
                mob = _agarose_mobility(bp, agarose, dna_form=form)
                row = max(0, min(height - 1, int(round(mob * (height - 1)))))
                band_info.append({
                    "bp":       bp,
                    "form":     form,
                    "mobility": mob,
                    "row":      row,
                })
            lane_results.append({
                "index":  li,
                "name":   lane["name"],
                "source": lane["source"],
                "detail": lane["detail"],
                "bands":  band_info,
            })
    except Exception as exc:
        _log.exception("agent simulate-gel: band computation failed")
        return ({"error": f"band computation failed: {_scrub_path(str(exc))}"}, 500)
    response: dict = {
        "ok":          True,
        "agarose_pct": min(_AGAROSE_CHOICES,
                            key=lambda g: abs(g - agarose)),
        "height":      height,
        "lane_width":  lane_width,
        "lanes":       lane_results,
    }
    if include_image:
        try:
            rt = _render_gel_image(
                cleaned_lanes,
                template_seq=template,
                template_circular=template_circular,
                pcr_amplicon=pcr_amplicon,
                agarose_pct=agarose,
                height=height,
                lane_width=lane_width,
            )
            # `Text.plain` strips Rich styles — agent gets a plain
            # ASCII rendering it can paste into a terminal.
            response["image"] = rt.plain
        except Exception as exc:
            _log.exception("agent simulate-gel: image render failed")
            return ({"error": f"image render failed: {_scrub_path(str(exc))}"}, 500)
    source_mix: dict[str, int] = {}
    for lane in cleaned_lanes:
        source_mix[lane["source"]] = source_mix.get(lane["source"], 0) + 1
    _log_event("simulator.gel.agent",
                n_lanes=len(cleaned_lanes),
                agarose_pct=agarose,
                template_bp=len(template),
                sources=source_mix,
                include_image=include_image)
    return response


def _agent_entry_vector_check(grammar, gid, role, part):
    """For `design-gb-part`'s ``check_entry_vector``: would the designed
    ``part`` (its 4-nt ``oh5``/``oh3``) actually drop into the CONFIGURED entry
    vector for ``(gid, role)``? Digests that vector with the grammar's Type IIS
    enzyme and checks whether either released fragment's acceptor overhangs
    match the part's — in either orientation — using the exact `overhang_seq`
    equality that real ligation uses (the STRICT compatibility rule of
    `_clone_part_into_entry_vector`'s primer path).

    NB we deliberately do NOT run the full clone engine: it has a synthetic-
    insert FALLBACK that force-fits a closed product even when the overhangs
    don't match, so "it produced a plasmid" can't answer "is this compatible?".
    The overhang match is the actual MoClo/Golden-Braid compatibility criterion.

    Advisory — NEVER fails the design; the verdict rides in the response so a
    mismatch surfaces at DESIGN time not only at the clone step (agent-API
    feedback P1-5 residual)."""
    role = role or ""
    oh5 = (part.get("oh5") or "").upper()
    oh3 = (part.get("oh3") or "").upper()
    ev = _get_entry_vector(gid, role)
    if ev is None:
        return {"configured": False, "role": role,
                "part_oh5": oh5, "part_oh3": oh3,
                "note": f"no entry vector configured for grammar {gid!r}"
                        + (f" role {role!r}" if role else "")
                        + " — nothing to check the part against"}
    ev_name = ev.get("name") or "?"
    enzyme = grammar.get("enzyme") if isinstance(grammar, dict) else None
    if not enzyme:
        return {"configured": True, "checkable": False, "vector_name": ev_name,
                "role": role, "part_oh5": oh5, "part_oh3": oh3,
                "note": f"grammar {gid!r} declares no Type IIS enzyme to "
                        "digest the entry vector with"}
    try:
        ev_rec = _gb_text_to_record(ev.get("gb_text") or "")
        ev_seq = str(getattr(ev_rec, "seq", "") or "").upper()
        frags, err = (_excise_fragment_pair(
                          ev_seq, [enzyme], circular=True,
                          features=_ev_frag_input_features(ev_rec))
                      if ev_seq else ([], {"error": "no sequence"}))
    except Exception as exc:
        _log.exception("agent design-gb-part: entry-vector digest failed")
        return {"configured": True, "checkable": False, "vector_name": ev_name,
                "role": role, "part_oh5": oh5, "part_oh3": oh3,
                "note": f"couldn't digest entry vector: {_scrub_path(str(exc))}"}
    if err is not None or len(frags) != 2:
        return {"configured": True, "checkable": False, "vector_name": ev_name,
                "role": role, "enzyme": enzyme, "part_oh5": oh5, "part_oh3": oh3,
                "note": f"digesting {ev_name!r} with {enzyme} didn't yield a "
                        "clean 2-fragment acceptor (misconfigured / extra "
                        "sites?)"}
    # The stuffer (smaller fragment) carries the acceptor's boundary overhangs;
    # report them either way. Compatibility checks BOTH fragments (fwd + RC),
    # matching how the real clone locates the dropout regardless of size.
    # Report the DROPOUT's overhangs (its boundary == the acceptor's): identify
    # it by marker-ABSENCE, not size — size mislabels an acceptor whose dropout
    # outgrew its backbone. Fall back to the smaller fragment when the marker
    # signal is ambiguous. If a fragment actually MATCHES the part below, report
    # ITS overhangs so `compatible` and the shown overhangs never disagree.
    non_marker = [f for f in frags if not _fragment_has_backbone_marker(f)]
    stuffer = (non_marker[0] if len(non_marker) == 1
               else min(frags, key=lambda f: len(f.get("top_seq", ""))))
    acc5 = ((stuffer.get("left") or {}).get("overhang_seq") or "").upper()
    acc3 = ((stuffer.get("right") or {}).get("overhang_seq") or "").upper()
    orientation = None
    for f in frags:
        fl = ((f.get("left") or {}).get("overhang_seq") or "").upper()
        fr = ((f.get("right") or {}).get("overhang_seq") or "").upper()
        if not (oh5 and oh3 and fl and fr):
            continue
        if oh5 == fl and oh3 == fr:
            orientation = "forward"; acc5, acc3 = fl, fr; break
        if oh5 == _rc(fr) and oh3 == _rc(fl):
            orientation = "reverse"; acc5, acc3 = _rc(fr), _rc(fl); break
    compatible = orientation is not None
    return {"configured": True, "checkable": True, "compatible": compatible,
            "orientation": orientation, "vector_name": ev_name, "role": role,
            "enzyme": enzyme, "acceptor_oh5": acc5, "acceptor_oh3": acc3,
            "part_oh5": oh5, "part_oh3": oh3,
            "note": ("part overhangs match the acceptor — it will drop in"
                     if compatible else
                     f"part overhangs {oh5}/{oh3} don't match the {ev_name!r} "
                     f"acceptor {acc5}/{acc3} — it won't assemble into that "
                     "vector")}


@_agent_endpoint("design-gb-part")
def _h_design_gb_part(app, payload):
    """Design Golden Braid / MoClo domestication primers. Body:
    ``{template, start, end, part_type, grammar?, target_tm?,
        codon_taxid?, check_entry_vector?, entry_vector_role?}``.

    `template` is the source DNA sequence; `start`/`end` are 0-based
    half-open coords (wrap-aware: end < start means origin-spanning
    on a circular plasmid). `part_type` must be a position type in
    the named grammar (default `gb_l0`). When `part_type` is a
    coding type and `codon_taxid` resolves to a codon table,
    internal forbidden sites are silently repaired via synonymous
    substitution.

    If a singleton L0 entry vector is CONFIGURED for the grammar and its
    external acceptor overhangs differ from the grammar's category overhangs
    (a two-tier Golden-Braid UPD vector — external CTCG/TGAG vs category
    AATG/GCTT), the primers are automatically NESTED: the Type IIS cut presents
    the external pair so the fragment enters the vector, with the category pair
    one layer in (released later by the vector's own nested sites). The
    fragment's actual cut overhangs come back as ``result.entry_oh5`` /
    ``entry_oh3``. With no entry vector configured (or one whose overhangs
    already equal the category pair) the classic one-tier design is unchanged.

    Pass ``check_entry_vector: true`` to validate the designed part's 4-nt
    overhangs against the CONFIGURED entry vector for the grammar (agent-API
    feedback) — the design digests that vector, reads its acceptor overhangs,
    and reports whether the part will actually drop in, so a mismatch surfaces
    HERE instead of only at the clone step. ``entry_vector_role`` picks a
    Golden-Braid role (``Alpha1`` / ``Alpha2`` / ``Omega1`` / ``Omega2``;
    default ``""`` = the singleton L0 acceptor). The verdict rides under
    ``result.entry_vector_check`` (``compatible`` / ``acceptor_oh5`` /
    ``acceptor_oh3`` / ``part_oh5`` / ``part_oh3`` / ``note``); it's advisory
    and never fails the design.

    Returns the `_design_gb_primers` dict (insert_seq, primer pairs,
    mutations list, overhang labels).
    """
    template = payload.get("template") or payload.get("sequence")  # SC-G alias
    if not isinstance(template, str) or not template.strip():
        return ({"error": "missing or non-string 'template'"}, 400)
    template_clean = "".join(ch for ch in template.upper()
                              if ch in "ACGTRYWSMKBDHVN")
    if not template_clean:
        return ({"error": "no IUPAC bases in 'template'"}, 400)
    if len(template_clean) > _PAIRWISE_MAX_LEN:
        return ({"error": f"'template' exceeds "
                  f"{_PAIRWISE_MAX_LEN:,} bp cap"}, 413)
    start = _coerce_int(payload.get("start"), name="start")
    if isinstance(start, str):
        return ({"error": start}, 400)
    end = _coerce_int(payload.get("end"), name="end")
    if isinstance(end, str):
        return ({"error": end}, 400)
    n = len(template_clean)
    if not (0 <= start < n) or not (0 < end <= n):
        return ({"error": f"start/end out of range for "
                  f"{n}-bp template"}, 400)
    # 2026-05-27 (audit-5 domesticator M2): accept optional `strand`
    # parameter. The GUI's `_resolve_source` already respects the
    # feature strand (H1 fix); programmatic callers can now do the
    # same. Strand = -1 means the CDS is on the bottom strand; we
    # RC the slice before passing to `_design_gb_primers` so codon
    # repair runs in the right frame.
    strand = payload.get("strand", 1)
    try:
        strand = int(strand)
    except (TypeError, ValueError):
        return ({"error": "'strand' must be 1 or -1"}, 400)
    if strand not in (1, -1):
        return ({"error": "'strand' must be 1 or -1"}, 400)
    if strand == -1:
        # Wrap-aware slice + RC. Mirrors `_resolve_source`'s logic.
        if end < start:
            fwd_slice = (template_clean[start:n] + template_clean[:end])
        else:
            fwd_slice = template_clean[start:end]
        template_clean = _rc(fwd_slice)
        start = 0
        end = len(template_clean)
        n = len(template_clean)
    part_type = payload.get("part_type")
    if not isinstance(part_type, str) or not part_type:
        return ({"error": "missing 'part_type'"}, 400)
    gid = payload.get("grammar", "gb_l0")
    if not isinstance(gid, str):
        return ({"error": "'grammar' must be string"}, 400)
    grammar = _all_grammars().get(gid)
    if grammar is None:
        return ({"error": f"unknown grammar id {gid!r}"}, 404)
    target_tm = payload.get("target_tm", 60.0)
    try:
        target_tm = float(target_tm)
    except (TypeError, ValueError):
        return ({"error": "'target_tm' must be a number"}, 400)
    if not (40.0 <= target_tm <= 90.0):
        return ({"error": "'target_tm' must be in [40, 90] °C"}, 400)
    codon_raw = None
    codon_taxid = payload.get("codon_taxid")
    if codon_taxid is not None:
        if not isinstance(codon_taxid, str):
            return ({"error": "'codon_taxid' must be string"}, 400)
        entry = _codon_tables_get(codon_taxid)
        if entry is None:
            return ({"error": f"unknown codon_taxid {codon_taxid!r}"},
                    404)
        codon_raw = entry.get("raw")
    # Two-tier Golden-Braid nesting: derive the destination L0 acceptor's
    # EXTERNAL overhangs from the configured entry vector (role="" = the
    # singleton domestication vector) and nest the grammar's category overhangs
    # inside them, so the designed fragment actually drops into that vector.
    # None (no vector configured / not cleanly digestible) → classic one-tier.
    enzyme = grammar.get("enzyme") if isinstance(grammar, dict) else None
    entry_overhangs = (_entry_vector_acceptor_overhangs(gid, "", enzyme)
                       if enzyme else None)
    try:
        result = _design_gb_primers(
            template_clean, start, end, part_type,
            target_tm=target_tm, codon_raw=codon_raw, grammar=grammar,
            entry_overhangs=entry_overhangs,
        )
    except Exception as exc:
        _log.exception("agent design-gb-part: design failed")
        return ({"error": f"design failed: {_scrub_path(str(exc))}"}, 500)
    if isinstance(result, dict) and result.get("error"):
        return ({"error": result["error"]}, 422)
    # Optional entry-vector compatibility check (agent-API feedback): does the
    # part's (oh5, oh3) actually match the configured acceptor? Advisory —
    # attached to the result, never fails the design.
    if isinstance(result, dict) and payload.get("check_entry_vector"):
        role = payload.get("entry_vector_role", "") or ""
        if not isinstance(role, str):
            return ({"error": "'entry_vector_role' must be a string"}, 400)
        try:
            result["entry_vector_check"] = _agent_entry_vector_check(
                grammar, gid, role,
                {"oh5": result.get("entry_oh5") or result.get("oh5") or "",
                 "oh3": result.get("entry_oh3") or result.get("oh3") or ""})
        except Exception as exc:
            _log.exception("agent design-gb-part: entry-vector check failed")
            result["entry_vector_check"] = {
                "checkable": False,
                "note": f"entry-vector check failed: {_scrub_path(str(exc))}"}
    return {"ok": True, "result": result}


@_agent_endpoint("design-synthesis-fragment")
def _h_design_synthesis_fragment(app, payload):
    """Wrap a sequence in the correct nested overhangs so the directly-
    SYNTHESISED fragment (order it as a gBlock — no PCR) becomes a Level-0 part
    once cloned into the configured UPD/pUPD entry vector with the grammar's
    Type IIS enzyme. The direct-synthesis counterpart to ``design-gb-part``
    (which designs PCR domestication primers for an existing template).

    Body: ``{sequence, part_type?, oh5?, oh3?, grammar?="gb_l0"}``.
      * ``sequence`` — the insert to wrap (IUPAC DNA).
      * ``part_type`` — a grammar POSITION type (``CDS`` / ``Promoter`` /
        ``Terminator`` …) whose category overhangs flank the insert; OR
      * ``oh5`` / ``oh3`` — CUSTOM 4-nt category overhangs (override part_type).
      * ``grammar`` — cloning grammar id (default ``gb_l0``).

    Auto-nests exactly like ``design-gb-part``: if the grammar's configured L0
    entry vector exposes EXTERNAL overhangs that differ from the category pair
    (a two-tier UPD vector — e.g. external CTCG/TGAG vs category AATG/GCTT), the
    category pair nests inside so the fragment cuts to the vector's external
    overhangs. With no entry vector configured, the category overhangs are used
    directly (one-tier).

    Returns ``{ok, result}`` — result = ``{fragment, length, oh5, oh3,
    entry_oh5, entry_oh3, part_type, enzyme, nested, entry_vector,
    internal_sites}``. ``entry_oh5``/``entry_oh3`` are the fragment's ACTUAL
    Type IIS-cut overhangs (what enters the vector); ``internal_sites`` (when
    non-empty) flags grammar-forbidden recognition sites in the insert that make
    the fragment self-cut — scrub them (codon-optimise) before ordering."""
    seq, serr = _sanitize_bases(payload.get("sequence") or payload.get("template"))
    if serr is not None:
        return ({"error": f"'sequence': {serr}"}, 400)
    if not seq:
        return ({"error": "missing/empty 'sequence'"}, 400)
    if len(seq) > _PAIRWISE_MAX_LEN:
        return ({"error": f"'sequence' exceeds {_PAIRWISE_MAX_LEN:,} bp cap"}, 413)
    gid = payload.get("grammar", "gb_l0")
    if not isinstance(gid, str):
        return ({"error": "'grammar' must be a string"}, 400)
    grammar = _all_grammars().get(gid)
    if grammar is None:
        return ({"error": f"unknown grammar id {gid!r}"}, 404)
    # Category overhangs: explicit oh5/oh3 override a part_type lookup.
    part_type = payload.get("part_type") or ""
    if payload.get("oh5") is not None or payload.get("oh3") is not None:
        # `_sanitize_bases` permits full IUPAC; a Type IIS sticky end must be
        # ACGT-only (ambiguous bases won't ligate) — match the GUI modal's
        # ACGT filter so an agent/Babs caller can't design a non-cloneable frag.
        c5, e5 = _sanitize_bases(payload.get("oh5"), max_len=8)
        c3, e3 = _sanitize_bases(payload.get("oh3"), max_len=8)
        _acgt = set("ACGT")
        if (e5 is not None or e3 is not None or len(c5) != 4 or len(c3) != 4
                or not set(c5) <= _acgt or not set(c3) <= _acgt):
            return ({"error": "custom 'oh5'/'oh3' must each be exactly "
                     "4 ACGT bases"}, 400)
        cat5, cat3 = c5, c3
    elif part_type:
        if not isinstance(part_type, str):
            return ({"error": "'part_type' must be a string"}, 400)
        pos = _grammar_position_by_type(grammar, part_type)
        if pos is None:
            avail = ", ".join(p.get("type", "?")
                              for p in grammar.get("positions", []))
            return ({"error": f"part_type {part_type!r} not in grammar "
                     f"{gid!r}. Available: {avail}"}, 400)
        cat5, cat3 = pos.get("oh5", ""), pos.get("oh3", "")
        # A position with blank overhangs (malformed custom grammar) would yield
        # an un-cloneable "L0 part" with no fusion sites — reject like the modal
        # (which filters positions lacking oh5/oh3).
        if len(cat5) != 4 or len(cat3) != 4:
            return ({"error": f"grammar position {part_type!r} has no valid "
                     "4-nt fusion overhangs"}, 400)
    else:
        return ({"error": "supply either 'part_type' or custom "
                 "'oh5'+'oh3'"}, 400)
    enzyme = grammar.get("enzyme") or ""
    ext = _entry_vector_acceptor_overhangs(gid, "", enzyme) if enzyme else None
    try:
        result = _build_synthesis_l0_fragment(
            seq, cat5, cat3, grammar=grammar,
            part_type=part_type, entry_overhangs=ext)
    except Exception as exc:
        _log.exception("agent design-synthesis-fragment: build failed")
        return ({"error": f"build failed: {_scrub_path(str(exc))}"}, 500)
    ev = _get_entry_vector(gid, "")
    result["entry_vector"] = (ev.get("name") if isinstance(ev, dict) else "") or ""
    return {"ok": True, "result": result}


_AS_OFF_TARGET_MAX = 20          # named off-target templates per call
_AS_SEED_LEN = 12                # 3'-anchor length for a binding call


def _agent_design_allele_specific(template: str, payload):
    """`design-primers` with ``mode: "allele_specific"``.

    Body: ``{template, variant_pos, alt_base?, orientation?="auto",
    partner_primer?, off_targets?, target_tm?=60, destabilize?=false,
    circular?=false}``.

    Designs the DISCRIMINATING primer — the one whose 3'-terminal base sits on
    the variant — and reports how specific it actually is. ``off_targets`` is
    ``{name: sequence}`` (or a list of ``{name, sequence}``); every template
    given, plus the main one, is searched for 3'-seeded binding sites, because
    "this primer is allele-specific" is a claim about where it does NOT bind
    and cannot be made from the primer alone.

    ``partner_primer`` is the common (non-discriminating) primer of the pair;
    when given, its binding sites are reported too and the expected product
    size is computed.
    """
    variant_pos = _coerce_int(payload.get("variant_pos"), name="variant_pos")
    if isinstance(variant_pos, str):
        return ({"error": variant_pos}, 400)
    n = len(template)
    if not (0 <= variant_pos < n):
        return ({"error": f"'variant_pos' {variant_pos} out of range for a "
                          f"{n}-bp template"}, 400)
    alt_raw = payload.get("alt_base")
    alt_base = None
    if alt_raw is not None:
        if not isinstance(alt_raw, str) or len(alt_raw.strip()) != 1 \
                or alt_raw.strip().upper() not in "ACGT":
            return ({"error": "'alt_base' must be a single A/C/G/T base"}, 400)
        alt_base = alt_raw.strip().upper()
    orientation = payload.get("orientation", "auto")
    if orientation not in ("auto", "fwd", "rev"):
        return ({"error": "'orientation' must be 'auto', 'fwd' or 'rev'"}, 400)
    try:
        target_tm = float(payload.get("target_tm", 60.0))
    except (TypeError, ValueError):
        return ({"error": "'target_tm' must be a number"}, 400)
    if not (40.0 <= target_tm <= 90.0):
        return ({"error": "'target_tm' must be in [40, 90] °C"}, 400)
    circular = bool(payload.get("circular", False))
    destabilize = bool(payload.get("destabilize", False))

    best = _design_allele_specific_primer(
        template, variant_pos, alt_base=alt_base, orientation=orientation,
        target_tm=target_tm, destabilize=destabilize, circular=circular)
    if best is None:
        return ({"error": f"could not place an allele-specific primer at "
                          f"{variant_pos} — the template has fewer than "
                          f"{_AS_MIN_LEN} bases on the needed side (pass "
                          f"circular: true if it wraps)"}, 422)

    # ── Specificity: where ELSE does this 3' end anchor? ─────────────
    off_raw = payload.get("off_targets")
    templates: "list[tuple[str, str]]" = [("template", template)]
    if off_raw is not None:
        items = []
        if isinstance(off_raw, dict):
            items = list(off_raw.items())
        elif isinstance(off_raw, list):
            for e in off_raw:
                if not isinstance(e, dict):
                    return ({"error": "'off_targets' list entries must be "
                                      "{name, sequence} objects"}, 400)
                items.append((e.get("name"), e.get("sequence")))
        else:
            return ({"error": "'off_targets' must be {name: sequence} or a "
                              "list of {name, sequence}"}, 400)
        if len(items) > _AS_OFF_TARGET_MAX:
            return ({"error": f"too many off_targets (max "
                              f"{_AS_OFF_TARGET_MAX})"}, 413)
        for nm, sq in items:
            nm = _sanitize_label(nm, max_len=120) or "?"
            clean, err = _sanitize_bases(sq if isinstance(sq, str) else "")
            if err or not clean:
                return ({"error": f"off_target {nm!r} rejected: "
                                  f"{err or 'empty'}"}, 400)
            if len(clean) > _PAIRWISE_MAX_LEN:
                return ({"error": f"off_target {nm!r} exceeds "
                                  f"{_PAIRWISE_MAX_LEN:,} bp"}, 413)
            templates.append((nm, clean))

    def _covers(site, pos, total):
        """Does this footprint span `pos`, wrapping if circular?"""
        if total <= 0:
            return False
        off = (pos - site["foot_start"]) % total
        return off < site["length"]

    specificity: "list[dict]" = []
    n_intended = 0
    n_off = 0
    for name, seq in templates:
        try:
            sites = _primer_binding_sites(best["seq"], seq, len(seq),
                                          circular=circular,
                                          seed_len=_AS_SEED_LEN)
        except ValueError as exc:
            return ({"error": f"binding check on {name!r} failed: {exc}"}, 400)
        rows = []
        for s in sites:
            # The INTENDED site is the one on the main template whose
            # footprint covers the variant. Everything else with a full 3'
            # seed is a specificity risk, because the seed is precisely what
            # allele-specific PCR stakes its discrimination on.
            intended = (name == "template"
                        and _covers(s, variant_pos, len(seq)))
            if intended:
                n_intended += 1
            else:
                n_off += 1
            rows.append({"strand": s["strand"], "start": s["foot_start"],
                         "ident_pct": s["ident_pct"],
                         "mismatches": s["mismatches"],
                         "intended": intended})
        specificity.append({"name": name, "n_sites": len(sites),
                            "sites": rows[:10]})
    # A primer carrying the NON-reference allele should NOT exact-seed on the
    # reference template — that failure to bind IS the discrimination, so it
    # is reported as a property rather than as a missing site.
    discriminates = (n_intended == 0) and not best["matches_template"]
    n_extra = n_off

    partner = None
    partner_raw = payload.get("partner_primer")
    if partner_raw is not None:
        pseq, err = _sanitize_bases(partner_raw
                                    if isinstance(partner_raw, str) else "")
        if err or not pseq:
            return ({"error": f"'partner_primer' rejected: "
                              f"{err or 'empty'}"}, 400)
        try:
            psites = _primer_binding_sites(pseq, template, n,
                                           circular=circular,
                                           seed_len=_AS_SEED_LEN)
        except ValueError as exc:
            return ({"error": f"'partner_primer' check failed: {exc}"}, 400)
        product = None
        if psites:
            p0 = psites[0]
            if best["orientation"] == "fwd" and p0["strand"] == -1:
                product = (p0["foot_start"] + p0["length"]
                           - best["binding_start"]) % (n or 1)
            elif best["orientation"] == "rev" and p0["strand"] == 1:
                product = (best["binding_end"] - p0["foot_start"]) % (n or 1)
        partner = {
            "seq": pseq, "tm": _primer_tm_safe(pseq),
            "n_sites": len(psites),
            "strand": (psites[0]["strand"] if psites else None),
            "product_size": product,
        }

    warnings: "list[str]" = []
    if best["clamp_note"]:
        warnings.append(best["clamp_note"])
    if n_extra > 0:
        warnings.append(f"{n_extra} additional 3'-seeded binding site(s) "
                        f"found — allele specificity depends on the 3' end "
                        f"being unique, so check the `specificity` block "
                        f"before ordering")
    if best["matches_template"] and alt_base is not None:
        warnings.append(f"'alt_base' {best['allele']} is the base the "
                        f"template already carries at {variant_pos} — this "
                        f"primer is specific for the REFERENCE allele")

    _log_event("primers.allele_specific", variant_pos=variant_pos,
               orientation=best["orientation"], gc_clamp=best["gc_clamp"],
               n_off_targets=len(templates) - 1, via="agent")
    return {"ok": True, "mode": "allele_specific", "result": best,
            "specificity": specificity, "partner": partner,
            "gc_clamp": best["gc_clamp"], "clamp_note": best["clamp_note"],
            # `discriminates` True means the primer's 3' seed does NOT match
            # the reference template — the intended behaviour for a primer
            # designed against a non-reference allele, and the thing a caller
            # would otherwise have to infer from an empty site list.
            "discriminates": discriminates,
            "n_intended_sites": n_intended,
            "n_offtarget_sites": n_off,
            "warnings": warnings,
            "ignored": _agent_ignored_keys(payload, {
                "template", "sequence", "mode", "variant_pos", "alt_base",
                "orientation", "partner_primer", "off_targets", "target_tm",
                "destabilize", "circular"})}


@_agent_endpoint("design-primers")
def _h_design_primers(app, payload):
    """Primer-pair design over a target region. Body:
    ``{template, start, end, mode?: "detection"|"cloning"|"generic",
        target_tm?: float, site_5?: str, site_3?: str}``.

    Wraps `_design_detection_primers`, `_design_cloning_primers_raw`,
    and `_design_generic_primers`:

      * **detection** (default) — picks a diagnostic amplicon WITHIN
        the selected region (no overhangs). Returns `fwd_seq`/`rev_seq`
        + `product_size`. Needs a region ≥ the minimum product size.
      * **cloning** — requires `site_5` + `site_3` recognition-site
        strings and appends them to the primers as RE-cloning tails.
        Returns `fwd_full`/`rev_full` (tail + binding) + `fwd_binding`/
        `rev_binding`. For cloning-grammar overhangs (Golden Braid /
        MoClo Type IIS) use `design-gb-part` instead.
      * **generic** — binding-only primers (NO tails / RE sites /
        overhangs): the forward anneals at the region start, the
        reverse (reverse-complement) at the end. The whole amplicon is
        the selected region. Returns `fwd_seq`/`rev_seq` +
        `fwd_pos`/`rev_pos`. Needs ≥ 18 bp.

    All modes return `fwd_tm`/`rev_tm`. The primer dict is under the
    `result` key, with `mode` echoed alongside.
    """
    template = payload.get("template") or payload.get("sequence")  # SC-G alias
    if not isinstance(template, str) or not template.strip():
        return ({"error": "missing or non-string 'template'"}, 400)
    template_clean = "".join(ch for ch in template.upper()
                              if ch in "ACGTRYWSMKBDHVN")
    if not template_clean:
        return ({"error": "no IUPAC bases in 'template'"}, 400)
    if len(template_clean) > _PAIRWISE_MAX_LEN:
        return ({"error": f"'template' exceeds "
                  f"{_PAIRWISE_MAX_LEN:,} bp cap"}, 413)
    # Allele-specific is dispatched BEFORE start/end are parsed: it is anchored
    # on a single base, not a region, so requiring a region it never uses would
    # be a mandatory parameter with no meaning.
    if payload.get("mode") == "allele_specific":
        return _agent_design_allele_specific(template_clean, payload)
    start = _coerce_int(payload.get("start"), name="start")
    if isinstance(start, str):
        return ({"error": start}, 400)
    end = _coerce_int(payload.get("end"), name="end")
    if isinstance(end, str):
        return ({"error": end}, 400)
    n = len(template_clean)
    if not (0 <= start < n) or not (0 < end <= n):
        return ({"error": f"start/end out of range for "
                  f"{n}-bp template"}, 400)
    mode = payload.get("mode", "detection")
    if mode not in ("cloning", "detection", "generic", "allele_specific"):
        return ({"error":
                  "'mode' must be 'cloning', 'detection', 'generic', or "
                  "'allele_specific'"},
                400)
    target_tm = payload.get("target_tm", 60.0)
    try:
        target_tm = float(target_tm)
    except (TypeError, ValueError):
        return ({"error": "'target_tm' must be a number"}, 400)
    if not (40.0 <= target_tm <= 90.0):
        return ({"error": "'target_tm' must be in [40, 90] °C"}, 400)
    try:
        if mode == "detection":
            result = _design_detection_primers(
                template_clean, start, end, target_tm=target_tm,
            )
        elif mode == "generic":
            # Binding-only primers (no tails / RE sites / overhangs): the
            # forward anneals at the region start, the reverse at the end.
            result = _design_generic_primers(
                template_clean, start, end, target_tm=target_tm,
            )
        else:
            site_5 = payload.get("site_5", "")
            site_3 = payload.get("site_3", "")
            if not isinstance(site_5, str) or not isinstance(site_3, str):
                return ({"error": "cloning mode requires 'site_5' + "
                          "'site_3' recognition-site strings"}, 400)
            result = _design_cloning_primers_raw(
                template_clean, start, end,
                site_5=site_5.upper(), site_3=site_3.upper(),
                target_tm=target_tm,
            )
    except Exception as exc:
        _log.exception("agent design-primers: design failed")
        return ({"error": f"design failed: {_scrub_path(str(exc))}"}, 500)
    if isinstance(result, dict) and result.get("error"):
        return ({"error": result["error"]}, 422)
    return {"ok": True, "mode": mode, "result": result}


def _agent_project_experiments_list(project):
    """SC-M-class: the experiments belonging to project `project`. Empty/None
    → the ACTIVE project (live); a named project → its OWN stored experiments
    (a real partition for reads), 404 if missing. Parallels
    `_agent_primer_collection_list`/`_agent_parts_bin_list`."""
    if project in (None, ""):
        return list(_load_experiments()), None
    if not isinstance(project, str):
        return None, ({"error": "'project' must be a string"}, 400)
    p = _find_project(project.strip())
    if p is None:
        return None, ({"error":
                        f"no experiment project named {project.strip()!r}"}, 404)
    return list(p.get("experiments") or []), None


@_agent_endpoint("list-experiments")
def _h_list_experiments(app, payload):
    """List entries in an experiment project's notebook. Pass ``{project: "X"}``
    to list only project X's own entries (SC-M); omit it for the active
    project. Returns id, title, tag list, timestamps + body byte length per
    entry (omits `body_md` itself to keep responses small — fetch via
    ``get-experiment``)."""
    src, perr = _agent_project_experiments_list(payload.get("project"))
    if perr is not None:
        return perr
    assert src is not None               # perr is None ⇒ src is populated
    out: "list[dict]" = []
    for e in src:
        body = e.get("body_md") or ""
        body_bytes = len(body.encode("utf-8", errors="replace")) \
            if isinstance(body, str) else 0
        out.append({
            "id":          e.get("id", ""),
            "title":       e.get("title", ""),
            "tags":        list(e.get("tags") or []),
            "created_at":  e.get("created_at", ""),
            "updated_at":  e.get("updated_at", ""),
            "body_bytes":  body_bytes,
            "attached_plasmid_ids": list(
                e.get("attached_plasmid_ids") or [],
            ),
            "attached_gel_ids":     list(e.get("attached_gel_ids") or []),
            "attached_actions":     list(e.get("attached_actions") or []),
            "image_paths":          list(e.get("image_paths") or []),
        })
    pj = payload.get("project")
    return {"experiments": out,
            "active_project": _get_active_project_name() or "",
            "project": (pj.strip() if isinstance(pj, str) and pj.strip()
                        else "")}


@_agent_endpoint("delete-experiment", write=True)
def _h_delete_experiment(app, payload):
    """Delete a notebook entry by id, optionally scoped to a ``{project}``.
    Without ``project`` it deletes from the active project; with ``project`` it
    removes from THAT named project (SC-M). Refuses (409) to delete from a named
    project while it's ACTIVE — switch away first. Removes the attachment dir
    (`_EXPERIMENTS_DIR/<id>/`) too. Body: ``{id, project?}``."""
    eid = _sanitize_experiment_id(payload.get("id"))
    if eid is None:
        return ({"error": "missing or invalid 'id'"}, 400)
    proj_in = payload.get("project")
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        if proj_in not in (None, ""):
            if not isinstance(proj_in, str):
                return ({"error": "'project' must be a string"}, 400)
            pname = proj_in.strip()
            if pname == (_get_active_project_name() or ""):
                return ({"error":
                          f"{pname!r} is the active experiment project — "
                          f"set-active-experiment-project to another before "
                          f"deleting from it, so the live mirror stays "
                          f"consistent"}, 409)
            projs = _load_experiment_projects()
            tgt = next((p for p in projs
                        if (p.get("name") or "") == pname), None)
            if tgt is None:
                return ({"error":
                          f"no experiment project named {pname!r}"}, 404)
            entries = tgt.get("experiments") or []
            new_entries = [e for e in entries if e.get("id") != eid]
            if len(new_entries) == len(entries):
                return ({"error":
                          f"no experiment with id {eid!r} in project "
                          f"{pname!r}"}, 404)
            tgt["experiments"] = new_entries
            err = _agent_save_or_500(
                lambda: _save_experiment_projects(projs),
                "experiment_projects")
            if err:
                return err
            try:
                _state._delete_experiment_attach_dir_hook(eid)
            except OSError as exc:
                _log.warning("agent delete-experiment: attach cleanup "
                             "failed: %s", exc)
            return {"ok": True, "id": eid, "remaining": len(new_entries),
                    "project": pname}
        entries = _load_experiments()
        new_entries = [e for e in entries if e.get("id") != eid]
        if len(new_entries) == len(entries):
            return ({"error": f"no experiment with id {eid!r}"}, 404)
        err = _agent_save_or_500(
            lambda: _save_experiments(new_entries), "Experiments",
        )
        if err:
            return err
    try:
        _state._delete_experiment_attach_dir_hook(eid)
    except OSError as exc:
        _log.warning(
            "agent delete-experiment: attach dir cleanup failed: %s",
            exc,
        )
    _log_event("experiments.delete", eid=eid, via="agent")
    return {"ok": True, "id": eid,
            "remaining": len(new_entries)}


@_agent_endpoint("move-experiment", write=True)
def _h_move_experiment(app, payload):
    """Move a notebook entry between projects without a create/delete round-trip
    (SC-M, the experiment analog of `move-primer`/`move-part`). Body:
    ``{id, to, from?}``. Removes the entry from ``from`` (auto-found if omitted;
    409 if it's in more than one project) and adds it to ``to`` — atomically.
    The active project's live mirror is re-synced when touched. 404 if the entry
    or a project is missing; 409 on an ambiguous id or an id already in ``to``."""
    if (derr := _agent_reject_dangerous_unknowns(
            payload, {"id", "to", "from"})) is not None:
        return derr
    eid = _sanitize_experiment_id(payload.get("id"))
    if eid is None:
        return ({"error": "missing or invalid 'id'"}, 400)
    to_in = payload.get("to")
    if not isinstance(to_in, str) or not to_in.strip():
        return ({"error": "missing 'to' (an experiment project name)"}, 400)
    to_name = to_in.strip()
    from_in = payload.get("from")
    with _state._cache_lock:
        projs = _load_experiment_projects()
        by_name = {(p.get("name") or ""): p for p in projs}
        active = _get_active_project_name() or ""
        if to_name not in by_name:
            return ({"error":
                      f"no experiment project named {to_name!r}; create it "
                      f"with create-experiment-project"}, 404)
        if isinstance(from_in, str) and from_in.strip():
            from_name = from_in.strip()
            if from_name not in by_name:
                return ({"error":
                          f"no experiment project named {from_name!r}"}, 404)
        elif from_in in (None, ""):
            cands = [pn for pn, p in by_name.items()
                     if any(e.get("id") == eid
                            for e in (p.get("experiments") or []))]
            if not cands:
                return ({"error": "experiment not found in any project"}, 404)
            if len(cands) > 1:
                return ({"error":
                          "experiment is in multiple projects — pass 'from'",
                          "projects": cands}, 409)
            from_name = cands[0]
        else:
            return ({"error": "'from' must be a string"}, 400)
        if from_name == to_name:
            return ({"error": f"experiment already in {to_name!r}"}, 409)
        src = by_name[from_name].get("experiments") or []
        idxs = [i for i, e in enumerate(src) if e.get("id") == eid]
        if not idxs:
            return ({"error":
                      f"experiment {eid!r} not found in {from_name!r}"}, 404)
        entry = dict(src[idxs[0]])
        tgt = by_name[to_name].setdefault("experiments", [])
        if any(e.get("id") == eid for e in tgt):
            return ({"error":
                      f"an experiment with id {eid!r} already exists in "
                      f"{to_name!r}"}, 409)
        del by_name[from_name]["experiments"][idxs[0]]
        tgt.append(entry)
        if (err := _agent_save_or_500(
                lambda: _save_experiment_projects(projs),
                "experiment_projects")) is not None:
            return err
        # Re-sync the live mirror if the ACTIVE project's entries changed.
        if active in (from_name, to_name) and active in by_name:
            ab = by_name[active].get("experiments") or []
            if (err := _agent_save_or_500(
                    lambda: _save_experiments(ab), "Experiments")) is not None:
                return err
    return {"ok": True, "id": eid, "from": from_name, "to": to_name,
            "ignored": _agent_ignored_keys(payload, {"id", "to", "from"})}


@_agent_endpoint("list-experiment-projects")
def _h_list_experiment_projects(app, payload):
    """List experiment projects + the currently active one.
    Returns project name, description, save-stamp, and entry count
    per project (no body content — fetch via ``list-experiments``
    after switching active)."""
    projs = _load_experiment_projects()
    out: "list[dict]" = []
    for p in projs:
        raw_entries = p.get("experiments") or []
        n = len(raw_entries) if isinstance(raw_entries, list) else 0
        out.append({
            "name":        p.get("name", ""),
            "description": p.get("description", ""),
            "saved":       p.get("saved", ""),
            "n_entries":   n,
        })
    return {"projects": out,
            "active": _get_active_project_name() or ""}


@_agent_endpoint("set-active-experiment-project", write=True)
def _h_set_active_experiment_project(app, payload):
    """Switch the active experiment project. Body: ``{name: str}``.
    Mirrors the GUI flow (invariant #43 + sweep #9): updates the
    active-pointer, forces a sync settings flush, then writes the
    target project's entries into `experiments.json`. Power-loss
    safe across the two writes."""
    name = payload.get("name")
    if not isinstance(name, str) or not name:
        return ({"error": "missing 'name'"}, 400)
    # Sweep #26 (2026-05-25): hold `_state._cache_lock` across the
    # (find-target → set-pointer → save-disk → null-cache) flow so
    # a concurrent `_save_experiments` from another thread can't
    # interleave: pre-fix the OTHER thread's save could complete +
    # reseat `_state._experiments_cache` between our `_safe_save_json` and
    # our `_state._experiments_cache = None`, leaving the
    # cache in the OTHER thread's freshly-saved state — which the
    # next `_load_experiments` would then return, masking our
    # project switch. RLock allows the nested re-entry from
    # `_settings_flush_sync` (which may touch `_set_setting`'s
    # locked path).
    with _state._cache_lock:
        target = _find_project(name)
        if target is None:
            return ({"error": f"no project named {name!r}"}, 404)
        raw_entries = target.get("experiments") or []
        if not isinstance(raw_entries, list):
            raw_entries = []
        target_entries = [e for e in raw_entries if isinstance(e, dict)]
        prev = _get_active_project_name()
        _set_active_project_name(name)
        _state._settings_flush_sync_hook()
        try:
            # [INV-83, sweep #27] Agent project-switch mirror-swap;
            # source-of-truth is experiment_projects.json.
            _safe_save_json_mirror(
                _state._EXPERIMENTS_FILE, target_entries, "Experiments",
            )
        except (OSError, RuntimeError) as exc:
            # Mirror failed AFTER the pointer moved — revert so the active name
            # and the live experiments.json can't disagree (the desync would
            # clobber the target project on the next save). [INV-161]
            _set_active_project_name(prev)
            _state._settings_flush_sync_hook()
            _notify_save_failure(
                _state._LIVE_APP_REF.get(), "Experiments", exc,
            )
            return ({"error": f"save failed: {_scrub_path(str(exc))}"}, 500)
        _state._experiments_cache = None
    return {"ok": True, "active": name,
            "n_entries": len(target_entries)}


@_agent_endpoint("rename-experiment-project", write=True)
def _h_rename_experiment_project(app, payload):
    """Rename a project. Body: ``{name: str, new_name: str}``. If
    the renamed project is the active one, the active-pointer
    follows so subsequent ``_save_experiments`` mirrors land in the
    correct slot."""
    old = payload.get("name")
    if not isinstance(old, str) or not old.strip():
        return ({"error": "missing 'name'"}, 400)
    old = old.strip()
    new = _sanitize_label(payload.get("new_name"), max_len=200)
    if not new:
        return ({"error": "missing 'new_name'"}, 400)
    if old == new:
        return ({"error": "old and new names are identical"}, 400)
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        projs = _load_experiment_projects()
        if not any((p.get("name") or "").strip() == old for p in projs):
            return ({"error": f"no project named {old!r}"}, 404)
        if any((p.get("name") or "").strip() == new for p in projs):
            return ({"error": f"name {new!r} already taken"}, 409)
        for p in projs:
            if p.get("name") == old:
                p["name"]  = new
                p["saved"] = _date.today().isoformat()
                break
        err = _agent_save_or_500(
            lambda: _save_experiment_projects(projs),
            "Experiment projects",
        )
        if err:
            return err
        if _get_active_project_name() == old:
            _set_active_project_name(new)
            _state._settings_flush_sync_hook()
    _log_event("project.renamed", old=old, new=new, via="agent")
    return {"ok": True, "name": new}


@_agent_endpoint("delete-experiment-project", write=True)
def _h_delete_experiment_project(app, payload):
    """Delete a project + its entries. Body: ``{name: str}``.
    Refuses to delete the last remaining project (matches the GUI
    last-project guard). If the deleted project was active,
    promotes the first remaining project to active."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing 'name'"}, 400)
    name = name.strip()
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        projs = _load_experiment_projects()
        if not any((p.get("name") or "").strip() == name for p in projs):
            return ({"error": f"no project named {name!r}"}, 404)
        if len(projs) <= 1:
            return ({"error":
                      "cannot delete the last remaining project"}, 409)
        was_active = _get_active_project_name() == name
        new_projs = [p for p in projs if p.get("name") != name]
        err = _agent_save_or_500(
            lambda: _save_experiment_projects(new_projs),
            "Experiment projects",
        )
        if err:
            return err
        # The auto-promote + mirror-swap must stay under the same
        # _state._cache_lock hold (RLock, so the inner load/save re-enter
        # safely) — otherwise a concurrent _save_experiments can
        # reseat the cache between our save and the cache-null below,
        # desyncing cache from disk. Matches the locked UI sibling
        # ExperimentProjectsPickerModal._do_delete.
        promoted = ""
        if was_active and new_projs:
            promoted = new_projs[0].get("name") or ""
            _set_active_project_name(promoted)
            _state._settings_flush_sync_hook()
            raw_entries = new_projs[0].get("experiments") or []
            if not isinstance(raw_entries, list):
                raw_entries = []
            target_entries = [
                e for e in raw_entries if isinstance(e, dict)
            ]
            try:
                # [INV-83, sweep #27] Agent project-delete auto-promotes;
                # mirror-swap to the promoted project's entries.
                _safe_save_json_mirror(
                    _state._EXPERIMENTS_FILE, target_entries, "Experiments",
                )
            except (OSError, RuntimeError) as exc:
                # The primary source-of-truth write (the projects file) already
                # succeeded above; this mirror self-heals on the next reload.
                # Log it rather than swallowing silently so a recurring failure
                # is at least visible in the diagnostic trail.
                _log.warning(
                    "agent project-delete: mirror-swap to promoted project "
                    "failed (self-heals on next reload): %s", exc)
            _state._experiments_cache = None
    _log_event("project.deleted", name=name, via="agent",
                promoted=promoted)
    return {"ok": True, "name": name, "promoted": promoted}


# ── Notebook search / backlinks / protocol / export (2026-09-17) ────────────
# The read side of the lab notebook. Until these landed, the `attached_*`
# xrefs were rebuilt on every save and read by nothing, and the only way to
# find an entry by its contents was to `get-experiment` each id in turn.

def _agent_experiment_rows(payload) -> "tuple[list[tuple[str, dict]] | None, tuple | None]":
    """`[(project, entry), …]` for the scope named in `payload`.

    ``{project: "X"}`` scopes to one project, ``{all_projects: true}``
    spans every project, and the default is the active project only —
    the same scoping `list-experiments` uses, plus the cross-project
    option the GUI's Ctrl+F search has.
    """
    if payload.get("all_projects"):
        if payload.get("project"):
            return None, ({"error": "pass 'project' or 'all_projects', "
                                     "not both"}, 400)
        return _iter_all_experiments(), None
    src, perr = _agent_project_experiments_list(payload.get("project"))
    if perr is not None:
        return None, perr
    assert src is not None
    pname = payload.get("project") or (_get_active_project_name() or "")
    return [(pname, e) for e in src], None


def _agent_experiment_hit(project: str, entry: dict, *,
                           fields=(), snippet: str = "") -> dict:
    """Search-result row. Carries no `body_md` (a hit list of 1 MB bodies
    is not a response) — `get-experiment` fetches the text."""
    body = entry.get("body_md") or ""
    return {
        "id":         entry.get("id", ""),
        "title":      entry.get("title", ""),
        "project":    project,
        "tags":       list(entry.get("tags") or []),
        "created_at": entry.get("created_at", ""),
        "updated_at": entry.get("updated_at", ""),
        "body_bytes": len(body.encode("utf-8", errors="replace"))
                      if isinstance(body, str) else 0,
        "matched_in": list(fields),
        "snippet":    snippet,
    }


@_agent_endpoint("search-experiments")
def _h_search_experiments(app, payload):
    """Full-text search across notebook entries. Body:
    ``{query?: str, tags?: list[str], project?: str,
       all_projects?: bool, limit?: int}``.

    `query` terms must ALL appear, in the title, tags or body (AND
    across terms, OR across fields) — so ``"gibson failed"`` finds the
    Gibson entry that mentions failure without matching every entry
    that says "gibson". `tags` narrows to entries carrying ALL the named
    tags and works on its own, which is the "show me every Gibson entry"
    browse. Returns `matched_in` (which fields hit) and a `snippet` of
    body context per row; `body_md` itself comes from
    ``get-experiment``."""
    rows, err = _agent_experiment_rows(payload)
    if err is not None:
        return err
    assert rows is not None
    query = payload.get("query") or ""
    if not isinstance(query, str):
        return ({"error": "'query' must be a string"}, 400)
    tags = payload.get("tags") or []
    if not isinstance(tags, list):
        return ({"error": "'tags' must be a list of strings"}, 400)
    if not query.strip() and not tags:
        return ({"error": "pass 'query' and/or 'tags'"}, 400)
    limit = _coerce_int(payload.get("limit", 200), name="limit")
    if isinstance(limit, str):
        return ({"error": limit}, 400)
    limit = max(1, min(2000, limit))
    # Search per project so each hit keeps naming the project it came
    # from, then re-rank globally (score, then recency).
    grouped: "dict[str, list[dict]]" = {}
    for proj, e in rows:
        grouped.setdefault(proj, []).append(e)
    scored: "list[tuple[int, str, str, dict, list]]" = []
    for proj, entries in grouped.items():
        for h in _experiment_search(entries, query, tags=tags):
            scored.append((h["score"],
                           h["entry"].get("updated_at") or "",
                           proj, h["entry"], h["fields"]))
    scored.sort(key=lambda s: (s[0], s[1]), reverse=True)
    terms = _experiment_search_terms(query)
    out = [
        _agent_experiment_hit(
            proj, e, fields=fields,
            snippet=_experiment_snippet(e.get("body_md"), terms),
        )
        for _sc, _up, proj, e, fields in scored[:limit]
    ]
    return {"experiments": out, "count": len(out),
            "truncated": len(scored) > limit,
            "scope": ("all" if payload.get("all_projects")
                      else (payload.get("project")
                            or _get_active_project_name() or ""))}


@_agent_endpoint("experiment-backlinks")
def _h_experiment_backlinks(app, payload):
    """Which notebook entries reference an object — the reverse of the
    `@plasmid` / `!action` / `&gel` body tags. Body:
    ``{ref: str | list[str], kind?: "plasmid"|"action"|"gel",
       project?: str, all_projects?: bool}``.

    `ref` accepts a list so a plasmid can be asked for by BOTH its
    library entry id and its display name: the GUI's `Plasmid ref`
    button writes the id, while a hand-typed `@ref` is usually the name.
    Matching is case-insensitive.

    Refs are re-extracted from each body rather than read from the
    stored `attached_*` xref — the xref is rebuilt on every save, but a
    hand-edited `experiments.json` can carry a stale one, and answering
    "nothing references that" when something does is the failure this
    endpoint exists to prevent."""
    ref = payload.get("ref")
    refs = [ref] if isinstance(ref, str) else ref
    if not isinstance(refs, list) or not any(
            isinstance(r, str) and r.strip() for r in refs):
        return ({"error": "missing 'ref' (string or list of strings)"}, 400)
    kind = payload.get("kind") or "plasmid"
    if kind not in ("plasmid", "action", "gel"):
        return ({"error": "'kind' must be one of "
                           "plasmid, action, gel"}, 400)
    rows, err = _agent_experiment_rows(payload)
    if err is not None:
        return err
    assert rows is not None
    grouped: "dict[str, list[dict]]" = {}
    for proj, e in rows:
        grouped.setdefault(proj, []).append(e)
    out: "list[dict]" = []
    for proj, entries in grouped.items():
        for e in _experiments_referencing(entries, refs, kind=kind):
            out.append(_agent_experiment_hit(proj, e))
    out.sort(key=lambda h: h.get("updated_at") or "", reverse=True)
    return {"experiments": out, "count": len(out), "kind": kind,
            "ref": [r for r in refs if isinstance(r, str) and r.strip()]}


@_agent_endpoint("get-plasmid-protocol")
def _h_get_plasmid_protocol(app, payload):
    """A plasmid's recorded construction history as a bench protocol.
    Body: ``{id: str, collection?: str}``.

    Reads the same `_history_build_steps` the History viewer renders, so
    the answer is what the program actually recorded doing — not a
    reconstruction. Returns the structured `steps` AND a `markdown`
    rendering ready to paste into a notebook entry.

    Only names that resolve in the library become `@refs` in the
    markdown; the rest stay inline code, so the protocol never plants a
    dangling ref. Steps the program never saw — transformation,
    miniprep, sequencing — are absent by construction; use
    ``experiment-step-template`` for those."""
    key = payload.get("id") or payload.get("name")
    if not isinstance(key, str) or not key.strip():
        return ({"error": "missing 'id'"}, 400)
    key = key.strip()
    coll = payload.get("collection")
    if coll is not None and not isinstance(coll, str):
        return ({"error": "'collection' must be a string"}, 400)
    hits = _agent_scan_library_for_key(key, coll)
    if not hits:
        return ({"error": f"no library plasmid named {key!r}"
                           + (f" in collection {coll!r}" if coll else "")},
                404)
    _coll_name, entry = hits[0]
    name = entry.get("name") or key
    history_xml = entry.get("history_xml")
    if not history_xml:
        return ({"error": f"no construction history recorded for {name!r}"},
                404)
    try:
        root = _parse_commercialsaas_history(history_xml)
    except ValueError as exc:
        return ({"error": f"malformed history: {exc}"}, 422)
    if root is None:
        return ({"error": f"history for {name!r} is empty"}, 404)
    steps = _history_build_steps(root)
    known: "set[str]" = set()
    for s in steps:
        cands = [s.get("product"), s.get("backbone")]
        cands += list(s.get("inputs") or [])
        for c in cands:
            if (isinstance(c, str) and c and c not in known
                    and _experiment_ref_token_ok(c)
                    and _agent_scan_library_for_key(c)):
                known.add(c)
    return {
        "plasmid": name,
        "collection": _coll_name,
        "n_steps": len(steps),
        "steps": [
            {"op": s.get("op", ""), "product": s.get("product", ""),
             "inputs": list(s.get("inputs") or []),
             "backbone": s.get("backbone", ""),
             "enzymes": list(s.get("enzymes") or []),
             "where": s.get("where", ""),
             "seq_len": s.get("seq_len", 0),
             "circular": bool(s.get("circular"))}
            for s in steps
        ],
        "markdown": _protocol_steps_markdown(
            steps, product=name, known_ids=known),
    }


@_agent_endpoint("experiment-step-template")
def _h_experiment_step_template(app, payload):
    """A heading-per-step skeleton for the named bench actions. Body:
    ``{actions: list[str], heading?: str}``.

    The counterpart to ``get-plasmid-protocol``: that one emits what the
    program DID, this one emits the steps you are ABOUT to do, as empty
    headings carrying their `!action` tags in catalog order. Unknown
    action ids are accepted (the catalog is curated, not enforced) and
    titled by the id. Returns `markdown` to prepend/append to a body via
    ``update-experiment``; ``list-experiment-actions`` names the
    catalog."""
    actions = payload.get("actions")
    if not isinstance(actions, list) or not any(
            isinstance(a, str) and a.strip() for a in actions):
        return ({"error": "missing 'actions' (list of action ids)"}, 400)
    heading = payload.get("heading") or ""
    if not isinstance(heading, str):
        return ({"error": "'heading' must be a string"}, 400)
    picked = [a.strip() for a in actions
              if isinstance(a, str) and a.strip()]
    known = {row[1] for row in _EXPERIMENT_ACTIONS}
    md = _experiment_template_markdown(
        picked, _EXPERIMENT_ACTIONS, heading=heading,
    )
    return {"markdown": md, "actions": picked,
            "unknown_actions": [a for a in picked if a not in known]}


@_agent_endpoint("list-experiment-actions")
def _h_list_experiment_actions(app, payload):
    """The curated bench-action catalog behind the `!<id>` body tags —
    `[{group, id, description}, …]`. Curated, not enforced: a body may
    carry any `!my-step` tag that matches the ref pattern."""
    return {"actions": [
        {"group": g, "id": a, "description": d}
        for g, a, d in _EXPERIMENT_ACTIONS
    ], "count": len(_EXPERIMENT_ACTIONS)}


@_agent_endpoint("export-experiment")
def _h_export_experiment(app, payload):
    """Render a notebook entry — or a whole project — as a document.
    Body: ``{id?: str, project?: str, format?: "markdown"|"html"}``.

    Returns the document TEXT rather than writing a file, so the caller
    chooses the destination (and no agent path-safety check is needed).

    Markdown copies each body verbatim — cross-ref sigils intact, so an
    exported entry pastes back into a new one with working refs — and
    adds a metadata header, a resolved References block and the
    attachment list around it. HTML renders the documented markdown
    subset into a standalone page with the stylesheet inlined.

    Image attachments are referenced by their stored filename; the
    bytes are not embedded (that is the GUI export's job, where the
    attachment directory is at hand). Use ``get-experiment`` +
    ``list-experiments`` for the raw fields."""
    # Type-check BEFORE `.lower()`: `(payload.get("format") or "markdown")`
    # passes a truthy non-string straight through, so `{"format": 1}` used
    # to raise AttributeError — a 500 where a 400 belongs (found by
    # junk-payload fuzzing, 2026-09-17).
    fmt_raw = payload.get("format")
    if fmt_raw is None or fmt_raw == "":
        fmt_raw = "markdown"
    if not isinstance(fmt_raw, str):
        return ({"error": "'format' must be a string "
                           "('markdown' or 'html')"}, 400)
    fmt = fmt_raw.lower()
    if fmt not in ("markdown", "md", "html"):
        return ({"error": "'format' must be 'markdown' or 'html'"}, 400)
    eid_in = payload.get("id")
    src, perr = _agent_project_experiments_list(payload.get("project"))
    if perr is not None:
        return perr
    assert src is not None
    project = payload.get("project") or (_get_active_project_name() or "")
    if eid_in is not None:
        eid = _sanitize_experiment_id(eid_in)
        if eid is None:
            return ({"error": "invalid 'id'"}, 400)
        entry = next((e for e in src if e.get("id") == eid), None)
        if entry is None:
            return ({"error": f"no experiment with id {eid!r}"
                               + (f" in project {project!r}" if project
                                  else "")}, 404)
        entries = [entry]
    else:
        entries = sorted(src, key=lambda e: (e.get("updated_at") or ""),
                          reverse=True)
        if not entries:
            return ({"error": "no entries to export"}, 404)
    if fmt == "html":
        text = _experiment_html_document(
            entries,
            title=(entries[0].get("title") if len(entries) == 1
                    else project),
            project=project,
        )
    elif len(entries) == 1:
        text = _experiment_markdown_document(entries[0], project=project)
    else:
        text = _experiment_project_markdown(entries, project=project)
    return {"format": "html" if fmt == "html" else "markdown",
            "project": project, "n_entries": len(entries),
            "bytes": len(text.encode("utf-8", errors="replace")),
            "document": text}


@_agent_endpoint("duplicate-experiment", write=True)
def _h_duplicate_experiment(app, payload):
    """Copy a notebook entry under a fresh id — "start from the write-up
    I already have". Body: ``{id: str, title?: str}``.

    Operates on the ACTIVE project only, and refuses a ``project`` key
    rather than quietly duplicating whatever shares that id in the
    active project: unlike `create-experiment`, this both READS and
    WRITES, and doing that against a non-active project's stored entries
    while the live mirror holds the active one is the desync `[INV-161]`
    guards. Switch with ``set-active-experiment-project``, or duplicate
    here and ``move-experiment`` the copy.

    Image attachments are deliberately NOT copied: they live in a
    per-entry directory keyed by entry id, so carrying the paths over
    would point the copy into the ORIGINAL's directory — deleting either
    entry would then break the other's images. `attachments_copied` is
    always false and says so rather than leaving the caller to find
    out."""
    eid = _sanitize_experiment_id(payload.get("id"))
    if eid is None:
        return ({"error": "missing or invalid 'id'"}, 400)
    if payload.get("project"):
        return ({"error":
                  "duplicate-experiment works on the ACTIVE project only — "
                  "set-active-experiment-project first, or duplicate here "
                  "and move-experiment the copy"}, 400)
    title = payload.get("title")
    if title is not None and not isinstance(title, str):
        return ({"error": "'title' must be a string"}, 400)
    # Sweep #26: RMW under `_state._cache_lock`.
    with _state._cache_lock:
        entries = _load_experiments()
        src = next((e for e in entries if e.get("id") == eid), None)
        if src is None:
            return ({"error": f"no experiment with id {eid!r}"}, 404)
        new_id = _new_experiment_id(
            {i for i in (e.get("id") for e in entries)
             if isinstance(i, str)})
        copy = _experiment_duplicate(src, new_id=new_id, title=title)
        entries.append(copy)
        err = _agent_save_or_500(lambda: _save_experiments(entries),
                                  "experiments")
        if err:
            return err
    _log_event("experiments.duplicated", src=eid, eid=new_id, via="agent")
    return {"ok": True, "id": new_id, "source_id": eid,
            "title": copy.get("title", ""),
            "attachments_copied": False,
            "n_entries": len(entries)}


@_agent_endpoint("get-enzyme-collection")
def _h_get_enzyme_collection(app, payload):
    """Return one enzyme collection by `name`. Body: ``{name}``."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing or non-string 'name'"}, 400)
    name = name.strip()
    coll = _find_enzyme_collection(name)
    if coll is None:
        return ({"error": f"unknown enzyme collection {name!r}"}, 404)
    return {"ok": True, "collection": coll}


@_agent_endpoint("update-enzyme-collection", write=True)
def _h_update_enzyme_collection(app, payload):
    """Replace an existing enzyme collection by `name`. Body:
    ``{name, enzymes?: [str, ...], new_name?: str}``. ``new_name``
    renames the collection (and updates the active pointer if it
    matched). Omitting ``enzymes`` keeps the existing list."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing or non-string 'name'"}, 400)
    name = name.strip()
    new_name = payload.get("new_name")
    if new_name is not None:
        new_name = _sanitize_label(new_name, max_len=200)
        if not new_name:
            return ({"error": "'new_name' must be a non-empty string"}, 400)
    enzymes_arg = payload.get("enzymes")
    if enzymes_arg is not None:
        if not isinstance(enzymes_arg, list) or not all(
                isinstance(e, str) for e in enzymes_arg):
            return ({"error":
                      "'enzymes' must be a list of strings"}, 400)
    # Sweep #25: full RMW under `_state._cache_lock` (see
    # `_h_create_custom_enzyme` docstring for rationale).
    with _state._cache_lock:
        entries = _load_enzyme_collections()
        for i, e in enumerate(entries):
            if e.get("name") == name:
                updated = dict(e)
                if new_name is not None:
                    # Reject rename onto an already-taken name.
                    if new_name != name and any(
                            x.get("name") == new_name for x in entries):
                        return ({"error":
                                  f"name {new_name!r} already exists"}, 409)
                    updated["name"] = new_name
                if enzymes_arg is not None:
                    updated["enzymes"] = sorted(set(enzymes_arg))
                entries[i] = updated
                if (err := _agent_save_or_500(
                        lambda: _save_enzyme_collections(entries),
                        "Enzyme collections")) is not None:
                    return err
                # Update active pointer on rename.
                if new_name is not None and new_name != name and \
                        _get_active_enzyme_collection_name() == name:
                    _set_active_enzyme_collection_name(new_name)
                return {"ok": True, "name": updated["name"]}
    return ({"error": f"unknown enzyme collection {name!r}"}, 404)


@_agent_endpoint("delete-enzyme-collection", write=True)
def _h_delete_enzyme_collection(app, payload):
    """Delete an enzyme collection by `name`. If the deleted
    collection was the active one, clears the active pointer."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing or non-string 'name'"}, 400)
    name = name.strip()
    # Sweep #25: full RMW under `_state._cache_lock` (see
    # `_h_create_custom_enzyme` docstring for rationale).
    with _state._cache_lock:
        entries = _load_enzyme_collections()
        new_entries = [e for e in entries if e.get("name") != name]
        if len(new_entries) == len(entries):
            return ({"error": f"unknown enzyme collection {name!r}"}, 404)
        if (err := _agent_save_or_500(
                lambda: _save_enzyme_collections(new_entries),
                "Enzyme collections")) is not None:
            return err
    if _get_active_enzyme_collection_name() == name:
        _set_active_enzyme_collection_name(None)
    return {"ok": True, "name": name}


@_agent_endpoint("get-active-enzyme-collection")
def _h_get_active_enzyme_collection(app, payload):
    """Return the currently active enzyme collection name (or
    ``None`` if no collection is active — the scanner then uses the
    full catalog)."""
    return {"ok": True, "name": _get_active_enzyme_collection_name()}


# Read-only `get-active-*` counterparts for the `set-active-*` family.
# Pre-2026-06-22 only enzyme-collection had a getter, so an agent had to
# infer the active item from list flags to no-op safely. Each getter
# mirrors its setter's response key.


@_agent_endpoint("get-active-collection")
def _h_get_active_collection(app, payload):
    """Return the active plasmid collection name (or ``None`` if none is
    active). Read-only counterpart to ``set-active-collection``."""
    return {"ok": True, "name": _get_active_collection_name()}


@_agent_endpoint("get-active-codon-table")
def _h_get_active_codon_table(app, payload):
    """Return the active codon-table taxid preference (empty string if
    unset → Synthesis falls back to the first registry entry, K12 by
    convention). Counterpart to ``set-active-codon-table``."""
    return {"ok": True, "taxid": _get_setting("active_codon_table", "")}


@_agent_endpoint("get-active-primer-collection")
def _h_get_active_primer_collection(app, payload):
    """Return the active primer collection name (or ``None``).
    Counterpart to ``set-active-primer-collection``."""
    return {"ok": True, "name": _get_active_primer_collection_name()}


@_agent_endpoint("get-active-parts-bin")
def _h_get_active_parts_bin(app, payload):
    """Return the active parts-bin name (or ``None``). Counterpart to
    ``set-active-parts-bin``."""
    return {"ok": True, "name": _get_active_parts_bin_name()}


@_agent_endpoint("get-active-experiment-project")
def _h_get_active_experiment_project(app, payload):
    """Return the active experiment-project name (or ``None``).
    Counterpart to ``set-active-experiment-project``."""
    return {"ok": True, "name": _get_active_project_name()}


@_agent_endpoint("get-active-hmm-database")
def _h_get_active_hmm_database(app, payload):
    """Return the active HMM database id (empty string if unset).
    Counterpart to ``set-active-hmm-database``."""
    return {"ok": True, "active": _get_setting("hmm_db_active_id", "")}


@_agent_endpoint("set-active-enzyme-collection", write=True)
def _h_set_active_enzyme_collection(app, payload):
    """Set or clear the active enzyme collection pointer. Body:
    ``{name: str | null}`` — ``null`` clears the pointer. Refuses
    names that don't resolve to an existing collection."""
    name = payload.get("name")
    if name is not None and (not isinstance(name, str) or not name.strip()):
        return ({"error":
                  "'name' must be a non-empty string or null"}, 400)
    if isinstance(name, str):
        name = name.strip()
        if _find_enzyme_collection(name) is None:
            return ({"error":
                      f"unknown enzyme collection {name!r}"}, 404)
    _set_active_enzyme_collection_name(name)
    return {"ok": True, "name": name}


@_agent_endpoint("auto-detect-entry-vectors", write=True)
def _h_auto_detect_entry_vectors(app, payload):
    """Run the Type-IIS-digest auto-detection pass across every
    library entry and bind newly-detected acceptors to their grammar
    roles. Body: ``{grammar_ids?: list[str] = <all>}``.

    Sacred — existing user bindings are NEVER clobbered (invariant
    #47). Returns the consolidated summary string + per-grammar
    bound-role count. Use ``list-entry-vectors`` afterwards to see
    the new bindings."""
    grammar_ids = payload.get("grammar_ids")
    if grammar_ids is not None and not isinstance(grammar_ids, list):
        return ({"error":
                  "'grammar_ids' must be a list of grammar id strings"
                }, 400)
    entries = _load_library()
    summary = _state._auto_bind_entry_vectors_from_entries_hook(
        entries,
        grammar_ids=grammar_ids,
    )
    return {"ok": True, "summary": summary}


# ── OT-2 (Opentrons) liquid-handler endpoints ───────────────────────────────────
# Compile a plate-transfer plan to an Opentrons protocol, analyse it on the robot
# (server-side simulate — no motion), watch the robot's full state / sensors for a
# crash, and — gated behind an explicit confirm + a clean analysis + a pre-flight
# health check — run it on real hardware. The compiler + client live in the
# app-free `splicecraft_opentrons` sibling; these endpoints are its data / network
# surface. See docs/agent-api.md.

def _agent_ot2_host(payload):
    """Resolve the OT-2 host: explicit ``host`` in the body, else the persisted
    ``ot2_host`` setting. Returns "" if neither is set."""
    return str(payload.get("host") or _get_setting("ot2_host", "") or "").strip()


def _agent_ot2_plan(payload):
    """The transfer plan: the ``plan`` object if present, else the payload itself
    (so a caller can send the plan fields at the top level)."""
    plan = payload.get("plan")
    return plan if isinstance(plan, dict) else payload


@_agent_endpoint("ot2-compile")
def _h_ot2_compile(app, payload):
    """Compile a plate-transfer plan into an Opentrons OT-2 protocol (offline —
    no robot needed). Validates the plan against the built-in deck catalog first.

    Body: a transfer plan, either under ``plan`` or at the top level::

        {"pipette": "p300_single", "mount": "left",
         "tips": {"labware": "tiprack_300", "slot": 8},
         "labware": {"src": {"labware": "eppi_24", "slot": 1},
                     "dst": {"labware": "plate_24", "slot": 2}},
         "transfers": [{"from": "src:A1", "to": "dst:A1", "volume": 50}],
         "new_tip": "always"}

    Or a multi-step protocol via ``steps`` (which takes precedence over
    ``transfers``) — an ordered sequence of typed operations::

        "steps": [
          {"type": "transfer", "from": "src:A1", "to": "dst:A1", "volume": 50,
           "mix_after": [3, 30], "blow_out": true, "touch_tip": true},
          {"type": "distribute", "from": "src:A1", "to": ["dst:A1", "dst:A2"], "volume": 40},
          {"type": "consolidate", "from": ["src:A1", "src:A2"], "to": "dst:B1", "volume": 60},
          {"type": "mix", "at": "dst:A1", "volume": 100, "repetitions": 5},
          {"type": "delay", "seconds": 30, "message": "incubate"},
          {"type": "pause", "message": "swap plates"},
          {"type": "comment", "text": "done"}]

    A control-only protocol (delay / pause / comment) needs no tips or labware.

    Friendly labware aliases (``tiprack_300``, ``eppi_24``, ``plate_24``,
    ``plate_96``, ``reservoir_12`` …) expand to canonical Opentrons load names.

    Returns ``{summary, valid, errors, warnings, protocol?}`` — ``protocol`` (the
    ``.py`` text) is present only when the plan validates. Always re-check with
    ``ot2-analyze`` before running on hardware.
    """
    plan = _agent_ot2_plan(payload)
    if not isinstance(plan, dict):
        return ({"error": "missing transfer plan (send plan fields or a 'plan' object)"}, 400)
    try:
        report = _ot2._ot2_validate_plan(plan)
        summary = _ot2._ot2_plan_summary(plan)
        out = {"summary": summary, "valid": not report["errors"],
               "errors": report["errors"], "warnings": report["warnings"]}
        if not report["errors"]:
            out["protocol"] = _ot2._ot2_compile_protocol(plan)
    except _ot2.OT2Error as exc:
        return ({"error": str(exc)}, 400)
    return out


@_agent_endpoint("ot2-status")
def _h_ot2_status(app, payload):
    """Full OT-2 state snapshot — the monitoring surface. Reports reachability +
    versions, pipette OK flags + volume specs, motor engagement, deck / instrument
    calibration health, attached modules and their sensor readings, robot settings
    (variables), the status light, and — for the active run or a given
    ``run_id`` — live run / command state, plus any detected ``faults`` and an
    ``ok`` verdict (so a crash is visible directly).

    Body: ``{"host": "192.168.1.56", "run_id"?: "..."}`` (``host`` falls back to
    the persisted ``ot2_host`` setting).
    """
    host = _agent_ot2_host(payload)
    if not host:
        return ({"error": "no OT-2 host (pass 'host' or set the 'ot2_host' setting)"}, 400)
    run_id = payload.get("run_id")
    return _ot2._ot2_state(host, run_id=str(run_id) if run_id else None)


@_agent_endpoint("ot2-analyze")
def _h_ot2_analyze(app, payload):
    """Compile (if given a plan) and analyse a protocol on the robot itself — the
    robot's built-in simulate. Catches deck / labware / tip errors a run would
    hit, WITHOUT any motion. The pre-flight before ``ot2-run``.

    Body: ``{"host": ..., "plan"|<plan fields>}`` or
    ``{"host": ..., "protocol": "<py text>"}``. Returns ``{result, status, errors,
    commands, pipettes, labware, protocol_id}``; ``result`` is ``"ok"`` /
    ``"not-ok"``.
    """
    host = _agent_ot2_host(payload)
    if not host:
        return ({"error": "no OT-2 host (pass 'host' or set the 'ot2_host' setting)"}, 400)
    protocol = payload.get("protocol")
    if not (isinstance(protocol, str) and protocol.strip()):
        try:
            protocol = _ot2._ot2_compile_protocol(_agent_ot2_plan(payload))
        except _ot2.OT2Error as exc:
            return ({"error": f"cannot compile plan: {exc}"}, 400)
    try:
        return _ot2._ot2_analyze(host, protocol)
    except _ot2.OT2Error as exc:
        return ({"error": str(exc)}, 502)


@_agent_endpoint("ot2-run", write=True)
def _h_ot2_run(app, payload):
    """Run a plate-transfer protocol on the OT-2 — GATED PHYSICAL ACTUATION.

    Refuses to move the gantry unless BOTH the robot's own analysis passes AND the
    body carries ``{"confirm": true}`` — and it also refuses on an already-faulted
    robot (pre-flight state check). Throughout the run it monitors full robot
    state and, the instant a fault is detected (a failed command, an instrument
    fault, …), halts the run and reports it.

    Body: ``{"host": ..., "plan"|<plan fields> | "protocol": "<py>",
    "confirm": true, "wait"?: true, "stop_on_fault"?: true}``.

    * ``confirm`` (default false) — without it you get a dry, analysis-only result
      and NO motion.
    * ``wait`` (default true) — block until the run finishes and return the final
      result (``crashed``, ``faults``, ``run_status``, ``failed_commands``). Pass
      ``wait: false`` to start it and return immediately with ``run_id``, then
      poll ``ot2-status`` with that ``run_id`` to watch it live.

    A ``write`` endpoint: needs the agent bearer token; a dirty library blocks it
    unless you also pass ``{"force": true}``.
    """
    host = _agent_ot2_host(payload)
    if not host:
        return ({"error": "no OT-2 host (pass 'host' or set the 'ot2_host' setting)"}, 400)
    protocol = payload.get("protocol")
    offset_plan = None
    if not (isinstance(protocol, str) and protocol.strip()):
        plan = _agent_ot2_plan(payload)
        try:
            protocol = _ot2._ot2_compile_protocol(plan)
        except _ot2.OT2Error as exc:
            return ({"error": f"cannot compile plan: {exc}"}, 400)
        offset_plan = plan   # so any per-labware offsets in the plan are applied
    try:
        return _ot2._ot2_run_protocol(
            host, protocol,
            confirm=bool(payload.get("confirm")),
            poll=bool(payload.get("wait", True)),
            stop_on_fault=bool(payload.get("stop_on_fault", True)),
            offset_plan=offset_plan,
        )
    except _ot2.OT2Error as exc:
        return ({"error": str(exc)}, 502)


@_agent_endpoint("ot2-run-control", write=True)
def _h_ot2_run_control(app, payload):
    """Pause / resume / stop the OT-2's current run — MANUAL run control (the
    counterpart to ``ot2-run``'s automatic stop-on-fault). Lets you hold a run to
    intervene, resume it, or abort it without waiting for a fault.

    Body: ``{"host": ..., "action": "pause"|"resume"|"stop", "run_id"?: "..."}``.
    The active run is resolved automatically when ``run_id`` is omitted (so this
    also controls a run started from the Opentrons App).

    A ``write`` endpoint: needs the agent bearer token; a dirty library blocks it
    unless you also pass ``{"force": true}``.
    """
    host = _agent_ot2_host(payload)
    if not host:
        return ({"error": "no OT-2 host (pass 'host' or set the 'ot2_host' setting)"}, 400)
    action = str(payload.get("action", "") or "").strip().lower()
    if action not in ("pause", "resume", "stop", "cancel"):
        return ({"error": "missing/invalid 'action' (use pause | resume | stop)"}, 400)
    rid = payload.get("run_id")
    rid = str(rid) if rid else _ot2._ot2_active_run(host)
    if not rid:
        return ({"error": "no active run to control (nothing is running on the robot)"}, 409)
    try:
        return _ot2._ot2_run_control(host, action, run_id=rid)
    except _ot2.OT2Error as exc:
        return ({"error": str(exc)}, 502)


def _ot2_opt_float(payload, key, default=None):
    """Read an optional finite positive-or-any float from a payload (None if absent
    / unparseable)."""
    v = payload.get(key, default)
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    import math as _m
    return f if _m.isfinite(f) else None


@_agent_endpoint("ot2-normalize")
def _h_ot2_normalize(app, payload):
    """Compute per-sample volumes to NORMALISE concentration across a set of wells —
    the everyday plate-prep task (equalise DNA mass, or hit a target ng/µL). Pure
    compute; no robot needed.

    Body::

        {"items": [{"name": "pGFP-1", "well": "A1", "concentration": 82.0}, ...],
         "target_ng": 200,          # OR "target_conc": 20 + "final_volume": 20
         "final_volume"?: 20,       # µL, required for target_conc / diluent top-up
         "pipette"?: "p300_single", # its min/max default the volume floor/ceiling
         "min_vol"?, "max_vol"?, "resolution"?: 0.1,
         # optional — also emit transfer steps ready for ot2-compile:
         "src"?: "src", "dst"?: "dst", "dst_labware"?: "plate_96",
         "dst_wells"?: ["A1", ...], "diluent_ref"?: "buf:A1", "new_tip"?: "always"}

    Returns ``{normalized: [{name, well, sample_ul, diluent_ul, achieved_ng,
    achieved_conc, ok, warning}], warnings, n_ok}`` — plus ``steps`` when a
    ``src`` + ``dst`` (+ wells) target is given. Concentrations too low/high to hit
    target within the pipette range are flagged, never silently dropped.
    """
    items = payload.get("items")
    if not isinstance(items, list) or not items:
        return ({"error": "missing 'items' (a list of {name, well, concentration})"}, 400)
    spec = _ot2._OT2_PIPETTES.get(str(payload.get("pipette", "p300_single")))
    min_vol = _ot2_opt_float(payload, "min_vol", spec[0] if spec else 0.0) or 0.0
    max_vol = _ot2_opt_float(payload, "max_vol", spec[1] if spec else None)
    try:
        normalized = _ot2._ot2_normalize_volumes(
            items,
            target_ng=_ot2_opt_float(payload, "target_ng"),
            target_conc=_ot2_opt_float(payload, "target_conc"),
            final_volume=_ot2_opt_float(payload, "final_volume"),
            min_vol=min_vol, max_vol=max_vol,
            resolution=_ot2_opt_float(payload, "resolution", 0.1) or 0.1)
    except _ot2.OT2Error as exc:
        return ({"error": str(exc)}, 400)
    out = {"ok": True, "normalized": normalized,
           "n_ok": sum(1 for r in normalized if r.get("ok")),
           "warnings": [f"{r['name'] or r.get('well')}: {r['warning']}"
                        for r in normalized if r.get("warning")]}
    src, dst = payload.get("src"), payload.get("dst")
    if src and dst:
        dst_wells = payload.get("dst_wells")
        if not (isinstance(dst_wells, list) and dst_wells):
            dst_wells = _ot2._ot2_entry_wells({"labware": str(payload.get("dst_labware", "") or "")})
        diluent_ref = str(payload["diluent_ref"]) if payload.get("diluent_ref") else None
        # A dilution target with no diluent well would build undiluted samples while
        # the `normalized` preview reports the diluted concentration — warn loudly.
        if not diluent_ref and any((r.get("diluent_ul") or 0) > 0
                                   for r in normalized if r.get("ok")):
            out["warnings"].append(
                "samples need diluent to reach target but no 'diluent_ref' was given — "
                "the steps deliver undiluted sample only; set 'diluent_ref' (e.g. 'buf:A1')")
        out["steps"] = _ot2._ot2_normalize_steps(
            normalized, src_id=str(src), dst_id=str(dst),
            dst_wells=[str(w) for w in (dst_wells or [])],
            diluent_ref=diluent_ref, new_tip=str(payload.get("new_tip", "always")))
    return out


@_agent_endpoint("ot2-plate-map")
def _h_ot2_plate_map(app, payload):
    """Map a plasmid collection onto a labware's wells (row-major) — the "a plate IS
    a collection" link. Body: ``{"collection": "MyPlasmids", "labware": "eppi_24"}``.
    Returns ``{map: {well: {id, name}}, wells, n, overflow}`` (``overflow`` = how
    many plasmids didn't fit). Feed the wells into ``ot2-normalize`` / ``ot2-compile``
    to cherry-pick or replate the collection by identity, not bare coordinates.
    """
    coll = _normalize_collection_name(payload.get("collection"))
    if coll is None:
        return ({"error": "missing or invalid 'collection'"}, 400)
    col = next((c for c in _iter_collections_readonly()
                if isinstance(c, dict) and c.get("name") == coll), None)
    if col is None:
        return ({"error": f"no collection named {coll!r}"}, 404)
    labware = str(payload.get("labware", "") or "")
    if not labware:
        return ({"error": "missing 'labware' (a load name or alias)"}, 400)
    wells = _ot2._ot2_entry_wells({"labware": labware})
    if not wells:
        return ({"error": f"unknown labware {labware!r} — no known geometry"}, 400)
    plasmids = [p for p in (col.get("plasmids") or []) if isinstance(p, dict)]
    mapping = {w: {"id": p.get("id", ""), "name": p.get("name", "")}
               for w, p in zip(wells, plasmids)}
    return {"ok": True, "collection": coll,
            "labware": _ot2._ot2_resolve_labware(labware), "map": mapping,
            "wells": wells, "n": len(mapping),
            "overflow": max(0, len(plasmids) - len(wells))}


@_agent_endpoint("ot2-calibration")
def _h_ot2_calibration(app, payload):
    """OT-2 calibration status — deck calibration + per-pipette offset calibration +
    tip-length calibrations. Body: ``{host}``. Returns the ``calibration`` block plus
    a ``ready`` verdict (deck OK AND every mounted pipette calibrated) and, when not
    ready, ``needs_calibration`` (the mounts to calibrate in the Opentrons App).
    A liquid run (``ot2-run``) is blocked until the pipette is calibrated; a
    ``ot2-position-check`` can still run to verify alignment first.
    """
    host = _agent_ot2_host(payload)
    if not host:
        return ({"error": "no OT-2 host (pass 'host' or set the 'ot2_host' setting)"}, 400)
    st = _ot2._ot2_state(host)
    if not st.get("reachable"):
        return ({"error": "robot unreachable", "faults": st.get("faults", [])}, 502)
    cal = st.get("calibration") or {}
    needs = sorted(m for m, ok in (cal.get("pipettes_calibrated") or {}).items() if not ok)
    # Strict: only an explicit "OK" deck status counts as ready. "IDENTITY" is the
    # identity transform = deck NEVER calibrated; None = unknown/unavailable — neither
    # is "ready", so a caller can't trust `ready` and run on an uncalibrated deck.
    deck_ok = cal.get("deck_status") == "OK" and not cal.get("marked_bad")
    return {"ok": True, "calibration": cal, "ready": bool(deck_ok and not needs),
            "deck_ok": bool(deck_ok), "needs_calibration": needs,
            "instruments": st.get("instruments", [])}


@_agent_endpoint("ot2-home", write=True)
def _h_ot2_home(app, payload):
    """Home the OT-2 gantry (all axes to their reference). Safe motion — no plunger,
    no labware descent. Body: ``{host}``. A ``write`` endpoint (bearer token; a dirty
    library blocks it unless ``{"force": true}``)."""
    host = _agent_ot2_host(payload)
    if not host:
        return ({"error": "no OT-2 host (pass 'host' or set the 'ot2_host' setting)"}, 400)
    if _ot2._ot2_active_run(host):
        return ({"error": "a run is active — stop it before homing"}, 409)
    try:
        _ot2._ot2_home(host)
    except _ot2.OT2Error as exc:
        return ({"error": str(exc)}, 502)
    return {"ok": True, "homed": True}


@_agent_endpoint("ot2-lights", write=True)
def _h_ot2_lights(app, payload):
    """Toggle the OT-2 rail (status) lights — a zero-motion control, handy to
    blink / identify the robot, confirm you're talking to the right one, or turn the
    lights off when the bench is done for the night.

    Body: ``{host, "on": true|false}`` (``on`` defaults to true). A ``write``
    endpoint (bearer token; a dirty library blocks it unless ``{"force": true}``)."""
    host = _agent_ot2_host(payload)
    if not host:
        return ({"error": "no OT-2 host (pass 'host' or set the 'ot2_host' setting)"}, 400)
    on = payload.get("on", True)
    if isinstance(on, str):
        on = on.strip().lower() in ("1", "true", "yes", "on")
    try:
        res = _ot2._ot2_set_lights(host, bool(on))
    except _ot2.OT2Error as exc:
        return ({"error": str(exc)}, 502)
    state = res.get("on") if isinstance(res, dict) else None
    return {"ok": True, "on": bool(state) if state is not None else bool(on)}


@_agent_endpoint("ot2-disengage", write=True)
def _h_ot2_disengage(app, payload):
    """De-energise the OT-2 gantry motors so the carriage can be moved by hand and no
    motor holds torque — a zero-descent power-down (homes nothing; no plunger, no
    labware descent). Refused while a run is active (disengaging mid-run would drop
    the moving gantry).

    Body: ``{host, "axes"?: ["x","y","z","a","b","c"]}`` (defaults to all six axes).
    A ``write`` endpoint (bearer token; a dirty library blocks it unless
    ``{"force": true}``)."""
    host = _agent_ot2_host(payload)
    if not host:
        return ({"error": "no OT-2 host (pass 'host' or set the 'ot2_host' setting)"}, 400)
    if _ot2._ot2_active_run(host):
        return ({"error": "a run is active — stop it before disengaging the motors"}, 409)
    # Validate axes against the six-axis set (drop unknowns, cap the list); an empty
    # or absent list means the engine default (all six).
    raw = payload.get("axes")
    axes = None
    if isinstance(raw, list) and raw:
        valid = set(_ot2._OT2_ALL_AXES)
        axes = [a for a in (str(x).strip().lower() for x in raw[:6]) if a in valid] or None
    try:
        _ot2._ot2_disengage(host, axes=axes)
    except _ot2.OT2Error as exc:
        return ({"error": str(exc)}, 502)
    return {"ok": True, "disengaged": axes or list(_ot2._OT2_ALL_AXES)}


@_agent_endpoint("ot2-position-check", write=True)
def _h_ot2_position_check(app, payload):
    """Run a MOTION-ONLY position check — the gantry tours the deck (moves to the top
    of each loaded labware's wells) so you can verify alignment + labware offsets.
    GATED like a run: needs ``{"confirm": true}`` and the robot's analysis to pass,
    but NO aspirate/dispense/tip pick-up is ever emitted — the plunger is never
    actuated and the pipette never descends into a well, so it is safe on occupied
    slots. The pipette need NOT be calibrated (that is what you're checking); deck
    calibration + reachability are still required. Any plan labware ``offset`` is
    applied so you can see whether it lands the pipette dead centre.

    Body: ``{host, plan|<plan fields>, wells?: ["A1", ...], confirm?, wait?}``.
    Without ``confirm`` you get a dry, analysis-only result and NO motion.
    """
    host = _agent_ot2_host(payload)
    if not host:
        return ({"error": "no OT-2 host (pass 'host' or set the 'ot2_host' setting)"}, 400)
    plan = _agent_ot2_plan(payload)
    if not isinstance(plan, dict):
        return ({"error": "missing plan (send plan fields or a 'plan' object)"}, 400)
    p = _ot2._ot2_normalize_plan(plan)
    if not any(lw.get("slot") is not None for lw in p["labware"].values()) and \
       not any(t.get("slot") is not None for t in p["tips"]):
        return ({"error": "position check needs at least one labware / tip rack with a "
                          "slot on the deck"}, 400)
    wells = payload.get("wells")
    if isinstance(wells, list) and wells:
        wells = [str(w) for w in wells[:_ot2._OT2_MAX_POSCHECK_WELLS]]
    else:
        wells = None
    try:
        return _ot2._ot2_run_position_check(
            host, plan, wells=wells, confirm=bool(payload.get("confirm")),
            poll=bool(payload.get("wait", True)),
            stop_on_fault=bool(payload.get("stop_on_fault", True)))
    except _ot2.OT2Error as exc:
        return ({"error": str(exc)}, 502)


# ── OT-2 protocol library + custom-labware library endpoints ────────────────────
# Two single-file collection stores (collections embed items): protocols carry a
# "plan" (a compile-able multi-step design), custom labware carry a "definition"
# (Opentrons labware JSON). Same CRUD shape, so a small generic layer drives both.

def _cstore_items(load_fn, item_key: str) -> "list[dict]":
    out: "list[dict]" = []
    for c in load_fn():
        cn = c.get("name", "")
        for it in c.get(item_key, []) or []:
            out.append({"name": it.get("name"), "collection": cn})
    return out


def _cstore_find(load_fn, item_key: str, name):
    key = str(name or "").strip().casefold()
    for c in load_fn():
        for it in c.get(item_key, []) or []:
            if str(it.get("name", "")).strip().casefold() == key:
                return c.get("name", ""), it
    return None, None


def _cstore_collections(load_fn, item_key: str) -> "list[dict]":
    return [{"name": c.get("name", ""), "count": len(c.get(item_key, []) or [])}
            for c in load_fn()]


_OT2_MAX_STORE_BYTES = 512 * 1024   # per-item cap for the protocol / custom-labware stores


def _cstore_too_big(obj) -> bool:
    """True if `obj` serialises past the per-item store cap (defence-in-depth on
    top of the 1 MiB HTTP body cap — bounds a single pathological plan/def)."""
    import json
    try:
        return len(json.dumps(obj, default=str).encode("utf-8", "replace")) > _OT2_MAX_STORE_BYTES
    except (TypeError, ValueError):
        return False


def _cstore_save_item(load_fn, save_fn, item_key, coll, item, label):
    with _state._cache_lock:
        colls = load_fn()
        target = next((c for c in colls
                       if (c.get("name") or "").strip().casefold() == coll.casefold()), None)
        if target is None:
            target = {"name": coll, item_key: []}
            colls.append(target)
        target.setdefault(item_key, [])
        target[item_key] = [x for x in target[item_key] if x.get("name") != item.get("name")]
        target[item_key].append(item)
        return _agent_save_or_500(lambda: save_fn(colls), label)


def _cstore_delete_item(load_fn, save_fn, item_key, name, label, collection=None):
    # Delete ONE item — the first case-insensitive name match, optionally scoped to
    # a collection. Was: strip the name from EVERY collection, which silently nuked
    # same-named items elsewhere and used case-sensitive matching (find was not).
    # (Audit 2026-07-04, findings #2 + #5.)
    key = str(name or "").strip().casefold()
    coll_key = str(collection).strip().casefold() if collection else None
    with _state._cache_lock:
        colls = load_fn()
        for c in colls:
            if coll_key is not None and (c.get("name") or "").strip().casefold() != coll_key:
                continue
            items = c.get(item_key, []) or []
            for i, x in enumerate(items):
                if str(x.get("name", "")).strip().casefold() == key:
                    c[item_key] = items[:i] + items[i + 1:]
                    return _agent_save_or_500(lambda: save_fn(colls), label)
        return ({"error": f"{name!r} not found"}, 404)


def _cstore_create_collection(load_fn, save_fn, item_key, name, label):
    with _state._cache_lock:
        colls = load_fn()
        if any((c.get("name") or "").strip().casefold() == name.casefold() for c in colls):
            return ({"error": f"collection {name!r} already exists"}, 409)
        colls.append({"name": name, item_key: []})
        return _agent_save_or_500(lambda: save_fn(colls), label)


def _cstore_delete_collection(load_fn, save_fn, item_key, name, label):
    with _state._cache_lock:
        colls = load_fn()
        kept = [c for c in colls if (c.get("name") or "").strip().casefold() != name.casefold()]
        if len(kept) == len(colls):
            return ({"error": f"collection {name!r} not found"}, 404)
        return _agent_save_or_500(lambda: save_fn(kept), label)


@_agent_endpoint("list-protocols")
def _h_list_protocols(app, payload):
    """Every saved OT-2 protocol as ``[{name, collection}]``. Load a plan with
    ``get-protocol``; save one with ``save-protocol``."""
    return {"ok": True, "protocols": _cstore_items(_load_protocol_collections, "protocols")}


@_agent_endpoint("get-protocol")
def _h_get_protocol(app, payload):
    """The saved ``plan`` for a protocol. Body: ``{name}``. The plan is compile-able
    by ``ot2-compile`` / ``ot2-run``."""
    coll, it = _cstore_find(_load_protocol_collections, "protocols", payload.get("name"))
    if it is None:
        return ({"error": "protocol not found (pass 'name')"}, 404)
    return {"ok": True, "name": it.get("name"), "collection": coll, "plan": it.get("plan")}


@_agent_endpoint("save-protocol", write=True)
def _h_save_protocol(app, payload):
    """Save (or overwrite) a protocol. Body: ``{name, plan, collection?}``. ``plan`` is
    a compile-able plan (see ``ot2-compile``) — it is validated before saving. A new
    ``collection`` is created on demand (default ``Default``)."""
    name = _normalize_collection_name(payload.get("name"))
    if name is None:
        return ({"error": "missing or invalid 'name'"}, 400)
    plan = payload.get("plan")
    if not isinstance(plan, dict):
        return ({"error": "'plan' must be a compile-able plan object"}, 400)
    report = _ot2._ot2_validate_plan(plan)
    if report["errors"]:
        return ({"error": "invalid plan", "errors": report["errors"]}, 400)
    if _cstore_too_big(plan):
        return ({"error": f"plan too large (> {_OT2_MAX_STORE_BYTES // 1024} KB serialized)"}, 413)
    coll = _normalize_collection_name(payload.get("collection")) or "Default"
    err = _cstore_save_item(_load_protocol_collections, _save_protocol_collections, "protocols",
                            coll, {"name": name, "plan": plan}, "protocol_collections")
    return err if err else {"ok": True, "name": name, "collection": coll}


@_agent_endpoint("delete-protocol", write=True)
def _h_delete_protocol(app, payload):
    """Delete a saved protocol — the FIRST match for ``name`` (case-insensitive),
    optionally scoped to ``collection``. Body: ``{name, collection?}``."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing 'name'"}, 400)
    err = _cstore_delete_item(_load_protocol_collections, _save_protocol_collections, "protocols",
                              name, "protocol_collections", payload.get("collection"))
    return err if err else {"ok": True, "deleted": name}


@_agent_endpoint("list-protocol-collections")
def _h_list_protocol_collections(app, payload):
    """Named protocol collections as ``[{name, count}]``."""
    return {"ok": True,
            "protocol_collections": _cstore_collections(_load_protocol_collections, "protocols")}


@_agent_endpoint("create-protocol-collection", write=True)
def _h_create_protocol_collection(app, payload):
    """Create an empty named protocol collection. Body: ``{name}``."""
    name = _normalize_collection_name(payload.get("name"))
    if name is None:
        return ({"error": "missing or invalid 'name'"}, 400)
    err = _cstore_create_collection(_load_protocol_collections, _save_protocol_collections,
                                    "protocols", name, "protocol_collections")
    return err if err else {"ok": True, "name": name}


@_agent_endpoint("delete-protocol-collection", write=True)
def _h_delete_protocol_collection(app, payload):
    """Delete a protocol collection AND every protocol in it. Body: ``{name}``."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing 'name'"}, 400)
    err = _cstore_delete_collection(_load_protocol_collections, _save_protocol_collections,
                                    "protocols", name, "protocol_collections")
    return err if err else {"ok": True, "deleted": name}


@_agent_endpoint("list-custom-labware")
def _h_list_custom_labware(app, payload):
    """Every custom labware definition as ``[{name, collection}]``."""
    return {"ok": True, "custom_labware": _cstore_items(_load_custom_labware, "labware")}


@_agent_endpoint("get-custom-labware")
def _h_get_custom_labware(app, payload):
    """The Opentrons ``definition`` for a custom labware. Body: ``{name}``. Reference it
    from a plan's labware entry as ``{"labware": name, "slot": N, "definition": <def>}``."""
    coll, it = _cstore_find(_load_custom_labware, "labware", payload.get("name"))
    if it is None:
        return ({"error": "custom labware not found (pass 'name')"}, 404)
    return {"ok": True, "name": it.get("name"), "collection": coll,
            "definition": it.get("definition")}


@_agent_endpoint("save-custom-labware", write=True)
def _h_save_custom_labware(app, payload):
    """Save (or overwrite) a custom labware. Body: ``{name, definition, collection?}`` —
    ``definition`` is an Opentrons labware-definition object with a ``wells`` map."""
    name = _normalize_collection_name(payload.get("name"))
    if name is None:
        return ({"error": "missing or invalid 'name'"}, 400)
    definition = payload.get("definition")
    if not isinstance(definition, dict) or not isinstance(definition.get("wells"), dict):
        return ({"error": "'definition' must be an Opentrons labware object with a 'wells' map"}, 400)
    if _cstore_too_big(definition):
        return ({"error": f"definition too large (> {_OT2_MAX_STORE_BYTES // 1024} KB)"}, 413)
    coll = _normalize_collection_name(payload.get("collection")) or "Default"
    err = _cstore_save_item(_load_custom_labware, _save_custom_labware, "labware", coll,
                            {"name": name, "definition": definition}, "custom_labware")
    return err if err else {"ok": True, "name": name, "collection": coll}


@_agent_endpoint("delete-custom-labware", write=True)
def _h_delete_custom_labware(app, payload):
    """Delete a custom labware — the FIRST match for ``name`` (case-insensitive),
    optionally scoped to ``collection``. Body: ``{name, collection?}``."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing 'name'"}, 400)
    err = _cstore_delete_item(_load_custom_labware, _save_custom_labware, "labware", name,
                              "custom_labware", payload.get("collection"))
    return err if err else {"ok": True, "deleted": name}


@_agent_endpoint("list-labware-collections")
def _h_list_labware_collections(app, payload):
    """Named custom-labware collections as ``[{name, count}]``."""
    return {"ok": True,
            "labware_collections": _cstore_collections(_load_custom_labware, "labware")}


@_agent_endpoint("create-labware-collection", write=True)
def _h_create_labware_collection(app, payload):
    """Create an empty named custom-labware collection. Body: ``{name}``."""
    name = _normalize_collection_name(payload.get("name"))
    if name is None:
        return ({"error": "missing or invalid 'name'"}, 400)
    err = _cstore_create_collection(_load_custom_labware, _save_custom_labware, "labware", name,
                                    "custom_labware")
    return err if err else {"ok": True, "name": name}


@_agent_endpoint("delete-labware-collection", write=True)
def _h_delete_labware_collection(app, payload):
    """Delete a custom-labware collection AND every labware in it. Body: ``{name}``."""
    name = payload.get("name")
    if not isinstance(name, str) or not name.strip():
        return ({"error": "missing 'name'"}, 400)
    err = _cstore_delete_collection(_load_custom_labware, _save_custom_labware, "labware", name,
                                    "custom_labware")
    return err if err else {"ok": True, "deleted": name}
