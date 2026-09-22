"""splicecraft_experiments — experiment-entry processing (Phase D, layer L1).

The [SUB-experiments] lab-notebook entry pipeline, extracted from the hub:
entry normalisation + id minting/sanitising, the `@plasmid` / `!action` / `&gel`
cross-reference extractors, and the legacy-tag-format migration. Pure (app-free)
— operates on plain experiment dicts / body markdown. The data-safety pieces
(`_save_experiment_image` blob write, `_delete_experiment_attach_dir` filesystem
delete) stay hub-side; `_migrate_legacy_tag_format` is re-exported because the
hub-side body-readers + the `_state._migrate_experiment_body_hook` registration
also call it. Depends only on util (L0). Re-exported by the hub so `sc.<name>` +
every existing call site (modals, agent endpoints, notes rendering) resolves
unchanged.
"""
from __future__ import annotations

import re
import uuid as _uuid

from splicecraft_util import (_NOTE_CTRL_RE, _is_image_path, _now_iso,
                              _sanitize_label)


# Per-entry body cap. 1 MB of markdown is ~250 k words — far past any
# realistic single-entry use. Larger entries stutter the live preview
# anyway (Markdown re-renders on every debounce).
_EXPERIMENT_BODY_MAX_BYTES = 1_000_000

_EXPERIMENT_TITLE_MAX_LEN = 200

_EXPERIMENT_TAG_MAX_LEN   = 60

_EXPERIMENT_TAGS_MAX      = 20

# Plasmid cross-reference token: `@<id>` inline anywhere in the body.
# Single-sigil format (refactor 2026-05-18) so the editor displays
# just the tag id without the noisy `@plasmid:` prefix. The negative
# lookbehind rejects matches where `@` is preceded by a word char
# (avoids email-like patterns: `user@example.com` doesn't tag
# `example.com`). The first id char must be a letter so numeric
# prose like "rev 2 @ 5pm" doesn't trigger.
#
# Sweep #9 (2026-05-19) atomic-group pattern (`(?=(...))\1`):
# the id captures inside a lookahead THEN is consumed via the
# `\1` backreference, which prevents backtracking — without this
# trick the trailing `(?![;=])` reject would just shorten the
# match (e.g. `&amp;` would match `&am`). With it, any id
# followed by `;` or `=` is rejected ENTIRELY (HTML entities
# `&amp;`/`&nbsp;`/`&copy;`, URL params `?foo=bar`). Python 3.10
# lacks possessive quantifiers (`{0,63}+`) and atomic groups
# (`(?>...)`) so the lookahead+backref idiom is the portable
# stand-in. The captured id is still `m.group(1)` and the full
# match (sigil + id) is still `m.group(0)`.
_PLASMID_REF_RE = re.compile(
    r"(?<![\w@])@(?=([A-Za-z](?:[\w.\-]{0,62}\w)?))\1(?![;=])"
)

# Action cross-reference token: `!<id>` inline anywhere in the body.
# Same single-sigil rationale as `@<id>`. `!` doesn't conflict with
# markdown image syntax `![alt](url)` because our regex requires the
# next char to be a letter, while images require `[`. Same atomic-
# group + trailing-reject hardening as the plasmid pattern (sweep #9).
_ACTIONS_REF_RE = re.compile(
    r"(?<![\w!])!(?=([A-Za-z](?:[\w.\-]{0,62}\w)?))\1(?![;=])"
)

# Gel cross-reference token: `&<id>` inline anywhere in the body
# (2026-05-19). Distinct sigil from plasmid + action so the three
# object kinds stay visually separable in the editor.
#
# The sweep #9 atomic-group + trailing-reject hardening was a real
# bug fix here (not "cosmetic at worst" as the original comment
# said): the pre-fix regex matched the entity name inside any
# pasted HTML or markdown export (`&amp;`, `&nbsp;`, `&copy;`...),
# polluting `attached_gel_ids` on save, false-highlighting in the
# editor, and surfacing a misleading "no such gel" notify on
# Ctrl+G click-through.
_GEL_REF_RE = re.compile(
    r"(?<![\w&])&(?=([A-Za-z](?:[\w.\-]{0,62}\w)?))\1(?![;=])"
)

# Filesystem-id constraint. Entry ids are mechanically generated as
# `exp-<8 hex>` (see `_new_experiment_id`), but accept the wider
# `[A-Za-z0-9][A-Za-z0-9._-]{0,63}` form so a hand-edited JSON with a
# sensible custom id still loads. Rejects empty, separators, `..`,
# NUL, shell metacharacters so the id can be path-joined safely.
# Mirrors the spirit of `_dna_sidecar_path`'s sanitisation.
_EXPERIMENT_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._\-]{0,63}$")


def _sanitize_experiment_id(raw: object) -> "str | None":
    """Return `raw` (str) if it passes `_EXPERIMENT_ID_RE`, else `None`.

    Rejects non-strings, empty strings, NUL embeds, `..` traversal, and
    any path separator (forward OR back slash). Used at every callsite
    that joins an entry id under `_EXPERIMENTS_DIR`.
    """
    if not isinstance(raw, str) or not raw:
        return None
    if "\x00" in raw or ".." in raw or "/" in raw or "\\" in raw:
        return None
    if not _EXPERIMENT_ID_RE.match(raw):
        return None
    return raw


def _as_str_list(value: object) -> "list[str]":
    """The strings in `value`, or `[]` — never a string's CHARACTERS.

    Hardening 2026-09-17. `_normalise_experiment_entry` guards `tags` and
    `image_paths` with an `isinstance(list)` check, so anything the app
    saved is well-shaped — but a hand-edited `experiments.json` (or a
    direct call) can put a bare string where a list belongs, and a plain
    `[x for x in value if isinstance(x, str)]` then iterates it as
    CHARACTERS. That produced real nonsense rather than an error:
    `"tags": "gibson"` matched a filter for the tag `g`,
    `"image_paths": "img.png"` exported seven bogus attachment links
    (`[i](i)`, `[m](m)`, …), and a protocol step with `inputs: "abc"`
    read `` `a` + `b` + `c` ``. A bare string is one item or nothing, and
    treating it as one item would guess at intent — so it yields [].
    """
    if isinstance(value, str) or not isinstance(value, (list, tuple)):
        return []
    return [v for v in value if isinstance(v, str)]


def _new_experiment_id(existing: "set[str] | None" = None) -> str:
    """Generate a fresh `exp-<8 hex>` id. `existing` (the current
    entries' id set) is consulted to avoid collision; bounded retries
    keep the loop deterministic for tests that monkeypatch `_uuid`."""
    seen = existing or set()
    for _ in range(64):
        eid = f"exp-{_uuid.uuid4().hex[:8]}"
        if eid not in seen:
            return eid
    return f"exp-{_uuid.uuid4().hex}"


def _migrate_legacy_tag_format(body: str) -> str:
    """Rewrite legacy `@plasmid:<id>` / `@actions:<id>` tokens to the
    single-sigil format `@<id>` / `!<id>` (refactor 2026-05-18 — the
    editor now shows the bare tag id without a noisy prefix). One-way
    migration applied on load; once the migrated body lands back on
    disk through `_save_experiments`, the old format is gone."""
    if not body:
        return body
    if "@plasmid:" in body:
        body = body.replace("@plasmid:", "@")
    if "@actions:" in body:
        body = body.replace("@actions:", "!")
    return body


def _extract_plasmid_refs(body_md: str) -> "list[str]":
    """Return the unique plasmid ids referenced via `@<id>` in
    `body_md`, preserving first-appearance order. Used to maintain
    the denormalised `attached_plasmid_ids` xref on save."""
    if not body_md or "@" not in body_md:
        return []
    seen: "list[str]" = []
    seen_set: "set[str]" = set()
    for m in _PLASMID_REF_RE.finditer(body_md):
        ref = m.group(1)
        if ref and ref not in seen_set:
            seen.append(ref)
            seen_set.add(ref)
    return seen


def _extract_action_refs(body_md: str) -> "list[str]":
    """Return the unique action ids referenced via `!<id>` in
    `body_md`, preserving first-appearance order. Mirrors
    `_extract_plasmid_refs` — used to maintain the denormalised
    `attached_actions` xref on save (2026-05-18)."""
    if not body_md or "!" not in body_md:
        return []
    seen: "list[str]" = []
    seen_set: "set[str]" = set()
    for m in _ACTIONS_REF_RE.finditer(body_md):
        ref = m.group(1)
        if ref and ref not in seen_set:
            seen.append(ref)
            seen_set.add(ref)
    return seen


def _extract_gel_refs(body_md: str) -> "list[str]":
    """Return the unique gel ids referenced via `&<id>` in
    `body_md`, preserving first-appearance order. Mirrors
    `_extract_plasmid_refs` — used to maintain the denormalised
    `attached_gel_ids` xref on save (2026-05-19)."""
    if not body_md or "&" not in body_md:
        return []
    seen: "list[str]" = []
    seen_set: "set[str]" = set()
    for m in _GEL_REF_RE.finditer(body_md):
        ref = m.group(1)
        if ref and ref not in seen_set:
            seen.append(ref)
            seen_set.add(ref)
    return seen


def _normalise_experiment_entry(entry: dict, *, fresh: bool = False
                                  ) -> dict:
    """Normalise an entry dict in place-style (returns a new dict).

    Caps title length, drops empty tags, deduplicates the
    `attached_plasmid_ids` xref from the live body, bumps `updated_at`.
    `fresh=True` also stamps `created_at` (used on new-entry create).

    Truncates `body_md` to `_EXPERIMENT_BODY_MAX_BYTES`; over-cap input
    is rare (1 MB markdown) but a deterministic truncate beats a save
    refusal that loses the user's work outright.
    """
    out = dict(entry) if isinstance(entry, dict) else {}
    eid = _sanitize_experiment_id(out.get("id"))
    if eid is None:
        # Caller bug: never persist an entry without a valid id.
        # Fall back to a fresh id rather than crash a save batch.
        eid = _new_experiment_id()
    out["id"] = eid
    title = out.get("title")
    if not isinstance(title, str):
        title = ""
    # Strip control bytes (terminal-escape defence) + length cap — the
    # title renders in the experiments-list DataTable.
    out["title"] = _sanitize_label(title, max_len=_EXPERIMENT_TITLE_MAX_LEN)
    body = out.get("body_md")
    if not isinstance(body, str):
        body = ""
    # Sweep #30 (2026-05-28): strip terminal-control bytes (preserving
    # \t / \n so multi-paragraph Markdown survives) — body_md renders in
    # the Compose TextArea + Markdown view; pre-fix an agent create/update
    # could persist an ESC/OSC escape that fired when the note opened.
    # [INV-85]
    body = _NOTE_CTRL_RE.sub("", body)
    # Sweep #9 (2026-05-19): re-apply legacy `@plasmid:` / `@actions:`
    # tag migration on every save, not only on load. Without this,
    # a body that arrived into in-memory state via paste / import
    # AFTER the initial load (when `_migrate_legacy_tag_format`
    # already ran) would persist back to disk with the old format
    # and remain unhighlighted / unclickable until the next launch.
    if "@plasmid:" in body or "@actions:" in body:
        body = _migrate_legacy_tag_format(body)
    encoded = body.encode("utf-8", errors="replace")
    if len(encoded) > _EXPERIMENT_BODY_MAX_BYTES:
        # Sweep #9 (2026-05-19): byte-cap truncation in one pass.
        # Pre-fix iterated 1024-char shrinks re-encoding the whole
        # body each pass — quadratic on multi-MB non-ASCII bodies
        # (e.g. 3 MB-encoded Chinese / emoji-heavy markdown
        # triggered seconds of UI freeze on save). New approach:
        # slice the encoded bytes to the cap, decode with
        # `errors="ignore"` so a truncation mid-multibyte-sequence
        # drops the partial sequence cleanly. Single encode +
        # single decode regardless of body size.
        body = encoded[:_EXPERIMENT_BODY_MAX_BYTES].decode(
            "utf-8", errors="ignore",
        )
    out["body_md"] = body
    raw_tags = out.get("tags") or []
    tags: list[str] = []
    if isinstance(raw_tags, list):
        for t in raw_tags:
            if not isinstance(t, str):
                continue
            # Sweep #30 (2026-05-28): strip control bytes + cap via
            # _sanitize_label — tags render in the #exp-tags-input Input;
            # the prior strip()+slice let an agent smuggle an escape. [INV-85]
            t = _sanitize_label(t, max_len=_EXPERIMENT_TAG_MAX_LEN)
            if not t:
                continue
            tags.append(t)
            if len(tags) >= _EXPERIMENT_TAGS_MAX:
                break
    out["tags"] = tags
    out["attached_plasmid_ids"] = _extract_plasmid_refs(body)
    out["attached_actions"] = _extract_action_refs(body)
    out["attached_gel_ids"] = _extract_gel_refs(body)
    image_paths = out.get("image_paths") or []
    if not isinstance(image_paths, list):
        image_paths = []
    out["image_paths"] = [
        p for p in image_paths if isinstance(p, str) and p
    ]
    now = _now_iso()
    if fresh or not isinstance(out.get("created_at"), str):
        out["created_at"] = now
    out["updated_at"] = now
    return out


# ── Curated bench-action catalog ─────────────────────────────────────────
#
# `(group, id, description)` rows offered by `ActionsPickerModal` as
# `!<id>` tags. Curated for typical molecular-cloning workflows but NOT
# enforced — a hand-typed `!my-step` is just as valid, because
# `_ACTIONS_REF_RE` accepts any matching id. Lives here (L1) rather than
# hub-side because it is pure data with three readers.
_EXPERIMENT_ACTIONS: tuple = (
    ("Design",       "design",         "Design primers / vector / insert"),
    ("PCR",          "pcr",            "Standard PCR amplification"),
    ("PCR",          "soe-pcr",        "Splicing by overlap extension"),
    ("PCR",          "colony-pcr",     "Colony screen by PCR"),
    ("Modification", "mutagenesis",    "Site-directed mutagenesis"),
    ("Restriction",  "digest",         "Restriction digestion"),
    ("Restriction",  "ligate",         "Ligation"),
    ("Restriction",  "dephosphorylate", "Dephosphorylation (rSAP/CIP)"),
    ("Assembly",     "gibson",         "Gibson / NEBuilder assembly"),
    ("Assembly",     "golden-gate",    "Golden Gate assembly"),
    ("Assembly",     "moclo",          "MoClo assembly"),
    ("Assembly",     "in-fusion",      "In-Fusion assembly"),
    ("Purification", "miniprep",       "Plasmid miniprep"),
    ("Purification", "gel-extract",    "Gel extraction"),
    ("Purification", "pcr-cleanup",    "PCR cleanup"),
    ("Biological",   "transform",      "Bacterial transformation"),
    ("Biological",   "colony-pick",    "Colony picking"),
    ("Validation",   "gel-check",      "Analytical gel electrophoresis"),
    ("Validation",   "sanger-seq",     "Sanger sequencing"),
)


# ── Notebook search ──────────────────────────────────────────────────────
#
# A query is split into terms that must ALL appear (AND). Capped so a
# pathological paste ("a b c d …" × 500) can't turn one keystroke into a
# 500-pass scan of every body in the project.
_EXPERIMENT_SEARCH_MAX_TERMS  = 12
_EXPERIMENT_SEARCH_SNIPPET    = 120


def _experiment_search_terms(query: object) -> "list[str]":
    """Split a search query into deduplicated lowercase terms.

    Notebook search is substring-AND, deliberately NOT the fuzzy
    subsequence match `LibrarySearchModal` uses for plasmid names. A
    plasmid name is short and half-remembered, so a subsequence match
    helps ("puc19" → "pUC19-lacZ-K12"). An entry body is up to 1 MB of
    prose, where subsequence matching finds the letters of "gibson"
    scattered across three paragraphs and calls it a hit. Two remembered
    words that both genuinely appear is the query people mean.
    """
    if not isinstance(query, str):
        return []
    seen: "list[str]" = []
    seen_set: "set[str]" = set()
    for raw in query.split():
        t = raw.strip().lower()
        if not t or t in seen_set:
            continue
        seen.append(t)
        seen_set.add(t)
        if len(seen) >= _EXPERIMENT_SEARCH_MAX_TERMS:
            break
    return seen


def _experiment_match_fields(entry: dict, terms: "list[str]") -> "list[str]":
    """Fields of `entry` that the query matched, or `[]` for no match.

    Every term must appear in at least one of title / tags / body
    (AND across terms, OR across fields) — so "gibson failed" finds the
    entry titled "Gibson attempt 3" whose body says the ligation failed,
    without matching every entry that merely says "gibson".

    The returned field list is what drives the "where it hit" column:
    a title hit and a body hit are worth showing differently.
    """
    if not terms or not isinstance(entry, dict):
        return []
    title = (entry.get("title") or "")
    title_l = title.lower() if isinstance(title, str) else ""
    body = entry.get("body_md") or ""
    body_l = body.lower() if isinstance(body, str) else ""
    raw_tags = entry.get("tags") or []
    tags_l = " ".join(
        t.lower() for t in raw_tags if isinstance(t, str)
    ) if isinstance(raw_tags, list) else ""
    fields: "list[str]" = []
    for t in terms:
        hit = False
        if t in title_l:
            hit = True
            if "title" not in fields:
                fields.append("title")
        if t in tags_l:
            hit = True
            if "tags" not in fields:
                fields.append("tags")
        if t in body_l:
            hit = True
            if "body" not in fields:
                fields.append("body")
        if not hit:
            # AND semantics: one missing term disqualifies the entry.
            return []
    return fields


def _experiment_snippet(body: object, terms: "list[str]",
                        width: int = _EXPERIMENT_SEARCH_SNIPPET) -> str:
    """One-line context window around the earliest term hit in `body`.

    Newlines/tabs collapse to spaces so the snippet occupies exactly one
    DataTable row; ellipses mark a window that doesn't start/end at a
    body boundary. Returns "" when nothing matched in the body (the hit
    was title- or tag-only, and the caller says so instead).
    """
    if not isinstance(body, str) or not body or not terms:
        return ""
    low = body.lower()
    at = -1
    for t in terms:
        if not t:
            # `"".find` returns 0, which would silently snippet the head
            # of the body for a query that matched nothing.
            continue
        i = low.find(t)
        if i >= 0 and (at < 0 or i < at):
            at = i
    if at < 0:
        return ""
    half = max(8, width // 2)
    start = max(0, at - half)
    end = min(len(body), at + half)
    chunk = body[start:end]
    # Collapse every run of whitespace — a markdown body is full of
    # newlines and a raw slice would break the single-row layout.
    chunk = " ".join(chunk.split())
    if start > 0:
        chunk = "…" + chunk
    if end < len(body):
        chunk = chunk + "…"
    return chunk


def _experiment_search(entries: "list[dict]", query: object, *,
                       tags: "list[str] | None" = None,
                       limit: "int | None" = None) -> "list[dict]":
    """Rank notebook entries against `query` (and an optional tag filter).

    Returns `[{entry, fields, snippet, score}, …]`, best first. `tags`
    narrows to entries carrying ALL the named tags (case-insensitive) and
    works on its own — an empty query plus `tags=["gibson"]` is the
    "show me every Gibson entry" browse, which is why the tag filter is
    not folded into the query terms.

    NO criteria (empty query, no tags) returns every entry: this is a
    FILTER, and an empty filter matches everything — the same thing an
    empty filter box means in the entries pane. Callers that want "no
    criteria" to mean "no answer" say so themselves; the
    `search-experiments` endpoint refuses it with a 400 rather than
    quietly dumping the notebook.

    Scoring is deliberately coarse: a title hit outranks a tag hit
    outranks a body hit, ties broken by recency at the callsite's sort.
    There is no relevance model here worth pretending to.
    """
    terms = _experiment_search_terms(query)
    want_tags = [t.strip().lower() for t in (tags or [])
                 if isinstance(t, str) and t.strip()]
    out: "list[dict]" = []
    for e in entries or []:
        if not isinstance(e, dict):
            continue
        if want_tags:
            have = {t.lower() for t in _as_str_list(e.get("tags"))}
            if not all(w in have for w in want_tags):
                continue
        if terms:
            fields = _experiment_match_fields(e, terms)
            if not fields:
                continue
        else:
            # Tag-only browse: every tag-matching entry qualifies.
            fields = []
        score = 0
        if "title" in fields:
            score += 3
        if "tags" in fields:
            score += 2
        if "body" in fields:
            score += 1
        out.append({
            "entry":   e,
            "fields":  fields,
            "snippet": _experiment_snippet(e.get("body_md"), terms),
            "score":   score,
        })
    out.sort(
        key=lambda h: (h["score"], (h["entry"].get("updated_at") or "")),
        reverse=True,
    )
    if isinstance(limit, int) and limit > 0:
        return out[:limit]
    return out


# ── Backlinks (which entries mention this object) ────────────────────────

# `kind` → the body extractor that owns that sigil. Keyed by the same
# names the agent endpoints take so a typo'd kind fails loudly at the
# lookup rather than silently returning "no references".
_EXPERIMENT_REF_EXTRACTORS = {
    "plasmid": _extract_plasmid_refs,
    "action":  _extract_action_refs,
    "gel":     _extract_gel_refs,
}


def _experiments_referencing(entries: "list[dict]", ref_id: object, *,
                             kind: str = "plasmid") -> "list[dict]":
    """Entries whose body references `ref_id` — the reverse of `@<id>`.

    `ref_id` may be one id or a list of alternatives, matching ANY. The
    `Plasmid ref` button writes the library ENTRY ID, but a hand-typed
    ref is usually the display NAME, so the plasmid-side backlink asks
    for both — one of them missing notes the user knows they wrote is
    the whole failure mode this view exists to fix.

    Re-extracts the refs from each body rather than trusting the stored
    `attached_*` xref. The xref is rebuilt by
    `_normalise_experiment_entry` on every save, so for anything written
    through the app the two agree — but a hand-edited `experiments.json`
    can carry a stale index, and a backlink query that silently answers
    "none" is worse than one that costs a linear scan. The extractors
    short-circuit on a body with no sigil at all, so the common case is
    a substring check.

    Matching is case-insensitive: `@puc19` finds the entries you wrote
    about `pUC19`. Plasmid names differing only in case are not a real
    library, and the alternative is a backlink view that misses notes.
    """
    extract = _EXPERIMENT_REF_EXTRACTORS.get(kind)
    if extract is None:
        raise ValueError(
            f"unknown reference kind {kind!r} — expected one of "
            f"{sorted(_EXPERIMENT_REF_EXTRACTORS)}"
        )
    if ref_id is None:
        raw: "list" = []
    elif isinstance(ref_id, str):
        raw = [ref_id]
    elif isinstance(ref_id, (list, tuple, set, frozenset)):
        raw = list(ref_id)
    else:
        # Fail loud, like the unknown-`kind` raise above. A dict is
        # iterable (yielding its KEYS) and a generator is consumable
        # once, so silently coercing either would answer "nothing
        # references that" for a query that was simply malformed.
        raise TypeError(
            f"ref_id must be a string, a list/tuple/set of strings, or "
            f"None — got {type(ref_id).__name__}"
        )
    want = {r.strip().lower() for r in raw
            if isinstance(r, str) and r.strip()}
    if not want:
        return []
    out: "list[dict]" = []
    for e in entries or []:
        if not isinstance(e, dict):
            continue
        body = e.get("body_md")
        if not isinstance(body, str) or not body:
            continue
        refs = extract(body)
        if any(r.lower() in want for r in refs):
            out.append(e)
    out.sort(key=lambda e: (e.get("updated_at") or ""), reverse=True)
    return out


# ── Protocol from construction history ──────────────────────────────────
#
# A history step records what the PROGRAM did. The notebook's `!action`
# catalog names what happens at the BENCH. The two overlap only where the
# mapping is unambiguous, and this table is that overlap — nothing more.
# An op that isn't here keeps its verb as plain text and gets no ref,
# because a protocol line claiming `!ligate` for a step that might have
# been a Golden Gate reaction is worse than a line that just says what
# the history recorded. Wet-lab-only steps (transformation, miniprep,
# sequencing) are never inferred: the simulation never saw them, so the
# template flow (`_experiment_template_markdown`) is where those come from.
_PROTOCOL_ACTION_BY_OP: "dict[str, str]" = {
    "Golden Gate":   "golden-gate",
    "Gibson":        "gibson",
    "HiFi assembly": "gibson",     # NEBuilder HiFi — what `!gibson` covers
    "PCR":           "pcr",
    "mutagenesis":   "mutagenesis",
}

# Inputs listed per step before "+N more". Mirrors the History viewer's
# `_HISTORY_PROTOCOL_INPUT_MAX` (kept as its own constant so L1 needn't
# import history L2, and drift-guarded by a test).
_PROTOCOL_INPUT_MAX = 8

# A name is only emitted as a `@<id>` ref when it can BE one. The id group
# of `_PLASMID_REF_RE` is mirrored here as a full match: a display name
# with a space ("pUC19 lacZ cassette") would otherwise ref just `@pUC19`
# and silently point at the wrong plasmid — or nothing.
_EXPERIMENT_REF_ID_RE = re.compile(r"^[A-Za-z](?:[\w.\-]{0,62}\w)?$")


def _experiment_ref_token_ok(name: object) -> bool:
    """True when `name` can be written as a `@<id>` / `&<id>` token and
    round-trip through the body extractors unchanged."""
    return bool(isinstance(name, str) and name
                and _EXPERIMENT_REF_ID_RE.match(name))


def _protocol_action_for_op(op: object) -> str:
    """The `!action` id for a history op label, or "" when the mapping
    would be a guess. See `_PROTOCOL_ACTION_BY_OP`."""
    if not isinstance(op, str):
        return ""
    return _PROTOCOL_ACTION_BY_OP.get(op.strip(), "")


def _protocol_steps_markdown(steps: "list[dict]", *,
                             product: str = "",
                             known_ids: "set[str] | None" = None,
                             action_refs: bool = True) -> str:
    """Render `_history_build_steps` output as notebook markdown.

    Takes the step dicts rather than a history root so this stays pure
    (L1 cannot import history at L2) — the caller parses the `.dna`
    history and hands over the steps, the same "callers supply the
    parts" split `splicecraft_cassette` uses.

    `known_ids` is the set of names that really resolve in the library;
    only those become clickable `@<id>` refs. Everything else is emitted
    as inline code, so a protocol never plants a dangling ref that
    pollutes `attached_plasmid_ids` and then reports "no such plasmid"
    on Ctrl+G.
    """
    known = known_ids or set()

    def _name(raw: object) -> str:
        nm = raw if isinstance(raw, str) else ""
        if not nm:
            return "`?`"
        if nm in known and _experiment_ref_token_ok(nm):
            return f"@{nm}"
        return f"`{nm}`"

    head = "## Protocol"
    if isinstance(product, str) and product.strip():
        head += f" — {product.strip()}"
    lines = [head, ""]
    if not steps:
        lines.append(
            "*No construction steps recorded for this plasmid — the "
            "history holds a single record.*"
        )
        return "\n".join(lines) + "\n"
    lines.append(
        "*Reconstructed from the recorded construction history. Steps the "
        "program never saw — transformation, miniprep, sequencing — are "
        "not listed.*"
    )
    lines.append("")
    for i, s in enumerate(steps, 1):
        if not isinstance(s, dict):
            continue
        verb = (s.get("op") or "").strip()
        prod = s.get("product") or ""
        ins = _as_str_list(s.get("inputs"))
        backbone = s.get("backbone") or ""
        enzymes = _as_str_list(s.get("enzymes"))
        where = (s.get("where") or "").strip()
        bits: "list[str]" = []
        if verb:
            bits.append(f"**{verb}**")
        # In-place edit — a single input that IS the product reads as
        # "X → X"; name the product once and say what was done to it.
        if len(ins) == 1 and not backbone and ins[0] == prod:
            bits.append(f"{_name(prod)}")
        else:
            ingredients: "list[str]" = []
            if ins:
                # Same cap as the History viewer's protocol pane: a real
                # MoClo/GB assembly can list dozens of parts, and an
                # uncapped join put 200 names on one 2,500-character
                # line. Shown count + "+N more", so the line stays
                # readable and never claims fewer inputs than there were.
                shown = ins[:_PROTOCOL_INPUT_MAX]
                joined = " + ".join(_name(x) for x in shown)
                extra = len(ins) - len(shown)
                if extra > 0:
                    joined += f" +{extra} more"
                ingredients.append(joined)
            if backbone:
                lead = "into " if ingredients else ""
                ingredients.append(f"{lead}{_name(backbone)}")
            if ingredients:
                bits.append(f"{' '.join(ingredients)} → {_name(prod)}")
            else:
                bits.append(_name(prod))
        if enzymes:
            bits.append("· " + " + ".join(f"`{e}`" for e in enzymes))
        if where:
            bits.append(f"({where})")
        if action_refs:
            act = _protocol_action_for_op(verb)
            if act:
                bits.append(f"!{act}")
        lines.append(f"{i}. " + " ".join(bits))
    lines.append("")
    return "\n".join(lines) + "\n"


# ── Entry templates + duplication ───────────────────────────────────────

def _experiment_template_markdown(action_ids: "list[str]",
                                  catalog: "tuple | list" = (),
                                  *, heading: str = "") -> str:
    """A skeleton body: one section per picked `!action`, in catalog order.

    `catalog` is the hub's `_EXPERIMENT_ACTIONS` rows
    (`(group, id, description)`) — passed in rather than imported so this
    stays pure. Unknown ids still get a section (the catalog is curated,
    not enforced — a hand-typed `!my-step` is as valid as `!miniprep`),
    titled by the id itself.

    Sections are headings plus blank space on purpose. A template that
    pre-writes "Cells: — Plates: — Colonies:" imposes one lab's form on
    every entry; the heading and the ref are what you would have retyped.
    """
    desc_by_id = {}
    order = {}
    for i, row in enumerate(catalog or ()):
        try:
            _grp, aid, desc = row[0], row[1], row[2]
        except (IndexError, TypeError):
            continue
        if isinstance(aid, str) and aid:
            desc_by_id.setdefault(aid, desc if isinstance(desc, str) else aid)
            order.setdefault(aid, i)
    wanted = [a for a in _as_str_list(action_ids) if a.strip()]
    # Catalog order, not click order — a protocol reads design → PCR →
    # digest → ligate → transform, which is how the catalog is arranged.
    seen: "set[str]" = set()
    picked: "list[str]" = []
    for a in sorted(wanted, key=lambda x: order.get(x, 10_000)):
        if a in seen:
            continue
        seen.add(a)
        picked.append(a)
    lines: "list[str]" = []
    if heading.strip():
        lines += [f"# {heading.strip()}", ""]
    for aid in picked:
        heading_text = desc_by_id.get(aid, aid)
        # Only write a `!tag` the body extractor will actually match.
        # `_ACTIONS_REF_RE` caps an id at 64 chars and restricts its
        # charset, so a longer or oddly-shaped id would render a tag that
        # can never be highlighted, clicked or denormalised into
        # `attached_actions` — the same "never emit a ref that doesn't
        # round-trip" rule `_protocol_steps_markdown` follows.
        if _experiment_ref_token_ok(aid):
            lines.append(f"## {heading_text}  !{aid}")
        else:
            lines.append(f"## {heading_text}")
        lines.append("")
        lines.append("")
    return "\n".join(lines)


def _experiment_duplicate(entry: dict, *, new_id: str,
                          title: "str | None" = None) -> dict:
    """Copy `entry` under `new_id` — the "use this write-up as a template"
    move.

    `image_paths` is deliberately NOT copied. Attachments live in a
    per-entry directory keyed by entry id, so carrying the paths over
    would leave the duplicate pointing into the ORIGINAL's directory:
    deleting either entry would then break the other's images, and the
    per-entry 100 MB cap would be computed against a directory the copy
    doesn't own. The copy starts with no attachments and says so.
    """
    src = dict(entry) if isinstance(entry, dict) else {}
    base = src.get("title") if isinstance(src.get("title"), str) else ""
    out = {
        "id":      new_id,
        "title":   title if isinstance(title, str) and title.strip()
                   else (f"{base} (copy)" if base else "(untitled) (copy)"),
        "body_md": src.get("body_md") if isinstance(src.get("body_md"), str)
                   else "",
        "tags":    [t for t in (src.get("tags") or []) if isinstance(t, str)],
        "image_paths": [],
    }
    return _normalise_experiment_entry(out, fresh=True)


# ── Export ──────────────────────────────────────────────────────────────

def _experiment_reference_block(entry: dict) -> str:
    """A "References" section listing the objects the body cites.

    An exported entry leaves the program behind, so `@pUC19-lacZ` stops
    being clickable and becomes a bare token in someone else's PDF. The
    block restates the refs as a list under named headings, which is what
    makes the export readable as a standalone record.
    """
    body = entry.get("body_md") if isinstance(entry, dict) else ""
    if not isinstance(body, str) or not body:
        return ""
    groups = (
        ("Plasmids", _extract_plasmid_refs(body), "@"),
        ("Actions",  _extract_action_refs(body),  "!"),
        ("Gels",     _extract_gel_refs(body),     "&"),
    )
    parts: "list[str]" = []
    for label, refs, sigil in groups:
        if not refs:
            continue
        parts.append(f"**{label}:** "
                     + ", ".join(f"`{sigil}{r}`" for r in refs))
    if not parts:
        return ""
    return "## References\n\n" + "\n\n".join(parts) + "\n"


def _experiment_markdown_document(entry: dict, *,
                                  project: "str | None" = None,
                                  image_srcs: "dict[str, str] | None" = None
                                  ) -> str:
    """One entry as a standalone markdown document.

    The body is copied VERBATIM — cross-ref sigils and all — so an
    exported entry can be pasted back into a new entry and still carry
    working `@`/`!`/`&` refs. Everything the export adds (metadata
    header, references, attachment list) goes around the body, never
    through it.
    """
    e = entry if isinstance(entry, dict) else {}
    title = (e.get("title") or "").strip() or "(untitled)"
    lines = [f"# {title}", ""]
    meta: "list[str]" = []
    if isinstance(project, str) and project.strip():
        meta.append(f"- **Project:** {project.strip()}")
    if e.get("id"):
        meta.append(f"- **Entry:** `{e.get('id')}`")
    if e.get("created_at"):
        meta.append(f"- **Created:** {e.get('created_at')}")
    if e.get("updated_at"):
        meta.append(f"- **Updated:** {e.get('updated_at')}")
    tags = [t for t in (e.get("tags") or []) if isinstance(t, str)]
    if tags:
        meta.append("- **Tags:** " + ", ".join(tags))
    if meta:
        lines += meta + [""]
    body = e.get("body_md")
    if isinstance(body, str) and body.strip():
        lines += [body.rstrip(), ""]
    refs = _experiment_reference_block(e)
    if refs:
        lines += [refs]
    images = _as_str_list(e.get("image_paths"))
    if images:
        lines += ["## Attachments", ""]
        srcs = image_srcs or {}
        for p in images:
            src = srcs.get(p, p)
            # `![...]` only for a raster image. A trace / CSV / datasheet
            # gets a plain link — an image tag around a `.ab1` renders as
            # a broken-image box in every viewer.
            bang = "!" if _is_image_path(p) else ""
            lines.append(f"{bang}[{p}]({src})")
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _experiment_project_markdown(entries: "list[dict]", *,
                                 project: "str | None" = None,
                                 image_srcs: "dict[str, str] | None" = None,
                                 generated_at: "str | None" = None) -> str:
    """A whole project as one markdown document, newest entry first."""
    name = (project or "").strip() or "Notebook"
    # Count the entries actually RENDERED, not `len(entries)`: a
    # non-list argument (or a list with junk in it) otherwise produced a
    # header claiming "3 entries" above a document containing none.
    rows = [e for e in (entries or []) if isinstance(e, dict)] \
        if not isinstance(entries, str) else []
    out = [f"# {name}", ""]
    out.append(f"*{len(rows)} "
               f"{'entry' if len(rows) == 1 else 'entries'} — "
               f"exported {generated_at or _now_iso()}*")
    out.append("")
    for e in rows:
        doc = _experiment_markdown_document(
            e, project=None, image_srcs=image_srcs,
        )
        # Demote the per-entry H1 to H2 so the project name stays the
        # document's only top-level heading.
        out.append("---")
        out.append("")
        out.append(doc.replace("# ", "## ", 1) if doc.startswith("# ")
                   else doc)
    return "\n".join(out).rstrip() + "\n"


# The markdown subset the HTML exporter understands. Named and bounded on
# purpose: a half-implemented CommonMark that mangles a table or an
# indented code block is worse than one that documents what it covers.
# Anything outside the subset is emitted as escaped literal text, so the
# worst case is "this reads like the source", never lost content.
_MARKDOWN_HTML_SUBSET = (
    "headings, paragraphs, hard line breaks, bold, italic, inline code, "
    "fenced code blocks, links, images, bullet and numbered lists, "
    "blockquotes and horizontal rules"
)

_MD_HEADING_RE  = re.compile(r"^(#{1,6})\s+(.*)$")
_MD_HR_RE       = re.compile(r"^\s*(?:-{3,}|\*{3,}|_{3,})\s*$")
_MD_UL_RE       = re.compile(r"^\s*[-*+]\s+(.*)$")
_MD_OL_RE       = re.compile(r"^\s*\d+[.)]\s+(.*)$")
_MD_QUOTE_RE    = re.compile(r"^\s*>\s?(.*)$")
_MD_FENCE_RE    = re.compile(r"^\s*```+\s*([A-Za-z0-9_+-]*)\s*$")
_MD_CODE_RE     = re.compile(r"`([^`\n]+)`")
# Alt / label text may not contain a bracket of either kind: allowing '['
# let a label swallow a whole run of them, so a body of 20,000 '[' took
# quadratic time to export, and it is exactly the shape of a link nested in
# an image's alt text — which used to be rewritten INSIDE the <img> tag.
_MD_IMAGE_RE    = re.compile(r"!\[([^\[\]]*)\]\(([^)\s]+)\)")
_MD_LINK_RE     = re.compile(r"\[([^\[\]]+)\]\(([^)\s]+)\)")
_MD_BOLD_RE     = re.compile(r"\*\*(\S(?:[^*]*\S)?)\*\*")
# Only `*italic*`, never `_italic_`: a lab notebook is full of
# snake_case identifiers and file names, and turning `pUC19_v2` into
# italics mid-word is a mangling the user can't undo.
_MD_ITALIC_RE   = re.compile(r"(?<![\w*])\*(\S(?:[^*]*\S)?)\*(?![\w*])")


def _html_safe_url(raw: object) -> str:
    """Return `raw` if it is a safe href/src, else "".

    An exported notebook is an HTML file someone opens in a browser, so a
    `javascript:` href pasted into a body must not survive the export.
    Allowlist: http(s), mailto, in-page anchors, `data:image/` (the
    inlined attachments), and scheme-less relative paths. Anything with
    an unrecognised scheme is dropped.
    """
    if not isinstance(raw, str):
        return ""
    # A browser's URL parser deletes tab / CR / LF anywhere in a URL, so
    # "/\t/host" IS "//host": judge the string the browser will use.
    s = re.sub(r"[\t\r\n]", "", raw).strip()
    if not s:
        return ""
    low = s.lower()
    if low.startswith(("http://", "https://", "mailto:", "#")):
        return s
    if low.startswith("data:image/"):
        return s
    # `//host/x` is a protocol-RELATIVE URL: in the exported file a
    # browser resolves it to `https://host/x`, so a link pasted into a
    # body becomes an outbound request the moment the export is opened.
    # `\\host\share` is the Windows UNC equivalent. Explicit http(s) is
    # already allowed above, so refusing these costs no capability and
    # makes "scheme-less means relative path" actually true. (Hardening
    # 2026-09-17.)
    # Browsers also read '\' as '/' in an http(s) URL, so "/\host" and
    # "\/host" are protocol-relative too: any two leading slashes of either
    # kind.
    if len(s) >= 2 and s[0] in "/\\" and s[1] in "/\\":
        return ""
    # Scheme-less → relative path. A colon before the first slash means
    # some other scheme (javascript:, vbscript:, file:, splicecraft:).
    if ":" in s.split("/", 1)[0]:
        return ""
    return s


def _md_inline_to_html(text: str) -> str:
    """Inline markdown → HTML for one already-block-classified line.

    Every tag this emits is built from the RAW text, escaped exactly once,
    and parked behind a placeholder until the end — so no later pass can
    match inside a tag an earlier one produced. (Running the passes over the
    already-escaped, already-rewritten line let a link nested in an image's
    alt text be rewritten INSIDE the <img> tag, breaking out of the
    attribute, and escaped every '&' in a URL twice.)"""
    from html import escape as _esc
    stash: "list[str]" = []

    def _park(html: str) -> str:
        stash.append(html)
        return f"\x00{len(stash) - 1}\x00"
    out = str(text).replace("\x00", "")
    # Code spans first, so `**literal**` inside backticks stays literal.
    out = _MD_CODE_RE.sub(
        lambda m: _park(f"<code>{_esc(m.group(1), quote=False)}</code>"), out)

    def _img(m) -> str:
        src = _html_safe_url(m.group(2))
        alt = m.group(1)
        if not src:
            return _park(_esc(f"[image: {alt}]", quote=False))
        return _park(f'<img src="{_esc(src, quote=True)}" '
                     f'alt="{_esc(alt, quote=True)}">')
    out = _MD_IMAGE_RE.sub(_img, out)

    def _link(m) -> str:
        href = _html_safe_url(m.group(2))
        label = m.group(1)
        if not href:
            return label
        # Only the tags are parked; the label stays in the stream so it is
        # escaped and can still carry bold / italic.
        return (_park(f'<a href="{_esc(href, quote=True)}" rel="noopener '
                      f'noreferrer">') + label + _park("</a>"))
    out = _MD_LINK_RE.sub(_link, out)
    out = _esc(out, quote=False)
    out = _MD_BOLD_RE.sub(r"<strong>\1</strong>", out)
    out = _MD_ITALIC_RE.sub(r"<em>\1</em>", out)
    return re.sub(r"\x00(\d+)\x00", lambda m: stash[int(m.group(1))], out)


def _markdown_subset_to_html(md: object) -> str:
    """Convert the documented markdown subset to an HTML fragment.

    See `_MARKDOWN_HTML_SUBSET` for exactly what is covered. Input is
    HTML-escaped first, so raw HTML in a body renders as visible source
    rather than executing — an exported notebook is a file people open in
    a browser, and a body is free-form text from anywhere (pasted
    protocols, supplier notes).

    A single newline becomes a hard break rather than a space. Markdown's
    paragraph-folding is the wrong default for a lab notebook, where the
    line breaks in a pasted gel table or a list of colony counts are the
    formatting.
    """
    from html import escape as _esc
    text = md if isinstance(md, str) else ""
    # Strip the placeholder sentinel so a hand-crafted body can't smuggle
    # a fake code-span marker through `_md_inline_to_html`.
    text = text.replace("\x00", "")
    lines = text.splitlines()
    out: "list[str]" = []
    para: "list[str]" = []
    list_kind = ""        # "ul" | "ol" | ""
    quote: "list[str]" = []
    fence_lang = None
    fence: "list[str]" = []

    def _flush_para() -> None:
        if para:
            out.append("<p>" + "<br>".join(para) + "</p>")
            para.clear()

    def _flush_list() -> None:
        nonlocal list_kind
        if list_kind:
            out.append(f"</{list_kind}>")
            list_kind = ""

    def _flush_quote() -> None:
        if quote:
            out.append("<blockquote>" + "<br>".join(quote)
                        + "</blockquote>")
            quote.clear()

    def _flush_all() -> None:
        _flush_para()
        _flush_list()
        _flush_quote()

    for raw in lines:
        fence_m = _MD_FENCE_RE.match(raw)
        if fence_lang is not None:
            if fence_m:
                cls = (f' class="language-{_esc(fence_lang, quote=True)}"'
                       if fence_lang else "")
                out.append(f"<pre><code{cls}>"
                           + _esc("\n".join(fence), quote=False)
                           + "</code></pre>")
                fence.clear()
                fence_lang = None
            else:
                fence.append(raw)
            continue
        if fence_m:
            _flush_all()
            fence_lang = fence_m.group(1) or ""
            continue
        if not raw.strip():
            _flush_all()
            continue
        h = _MD_HEADING_RE.match(raw)
        if h:
            _flush_all()
            lvl = len(h.group(1))
            out.append(f"<h{lvl}>{_md_inline_to_html(h.group(2))}</h{lvl}>")
            continue
        if _MD_HR_RE.match(raw):
            _flush_all()
            out.append("<hr>")
            continue
        q = _MD_QUOTE_RE.match(raw)
        if q:
            _flush_para()
            _flush_list()
            quote.append(_md_inline_to_html(q.group(1)))
            continue
        _flush_quote()
        ul = _MD_UL_RE.match(raw)
        ol = _MD_OL_RE.match(raw) if ul is None else None
        li = ul if ul is not None else ol
        if li is not None:
            _flush_para()
            want = "ul" if ul is not None else "ol"
            if list_kind != want:
                _flush_list()
                out.append(f"<{want}>")
                list_kind = want
            out.append(f"<li>{_md_inline_to_html(li.group(1))}</li>")
            continue
        _flush_list()
        para.append(_md_inline_to_html(raw))
    if fence_lang is not None:
        # Unclosed fence — emit what we collected rather than dropping it.
        cls = (f' class="language-{_esc(fence_lang, quote=True)}"'
               if fence_lang else "")
        out.append(f"<pre><code{cls}>" + _esc("\n".join(fence), quote=False)
                   + "</code></pre>")
    _flush_all()
    return "\n".join(out)


_EXPERIMENT_HTML_CSS = """
:root {
  --bg: #ffffff; --fg: #1d2021; --muted: #6b7280; --rule: #e5e7eb;
  --code-bg: #f4f4f5; --accent: #2563eb; --quote: #9ca3af;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #17181c; --fg: #e6e6e6; --muted: #9ca3af; --rule: #33353b;
    --code-bg: #22242a; --accent: #7aa2f7; --quote: #4b5563;
  }
}
* { box-sizing: border-box; }
body {
  background: var(--bg); color: var(--fg); margin: 0;
  font: 16px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
        "Helvetica Neue", Arial, sans-serif;
}
main { max-width: 46rem; margin: 0 auto; padding: 2.5rem 16px 4rem; }
h1, h2, h3, h4, h5, h6 { line-height: 1.25; margin: 2rem 0 .6rem; }
h1 { font-size: 1.9rem; margin-top: 0; }
h2 { font-size: 1.4rem; border-bottom: 1px solid var(--rule);
     padding-bottom: .3rem; }
h3 { font-size: 1.15rem; }
p, ul, ol, blockquote, pre { margin: 0 0 1rem; }
ul, ol { padding-left: 1.4rem; }
li { margin: .2rem 0; }
a { color: var(--accent); }
code {
  background: var(--code-bg); border-radius: 4px; padding: .1em .35em;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: .88em;
}
pre {
  background: var(--code-bg); border-radius: 6px; padding: .8rem 1rem;
  overflow-x: auto;
}
pre code { background: none; padding: 0; font-size: .85em; }
blockquote {
  border-left: 3px solid var(--quote); margin-left: 0;
  padding: .1rem 0 .1rem 1rem; color: var(--muted);
}
hr { border: 0; border-top: 1px solid var(--rule); margin: 2rem 0; }
img { max-width: 100%; height: auto; border-radius: 4px; }
table { border-collapse: collapse; }
.sc-meta {
  color: var(--muted); font-size: .88rem; margin: 0 0 2rem;
  padding: .7rem 0 0; border-top: 1px solid var(--rule); list-style: none;
}
.sc-meta li { margin: .15rem 0; }
.sc-entry { margin: 0 0 3rem; }
.sc-foot {
  color: var(--muted); font-size: .8rem; margin-top: 3rem;
  border-top: 1px solid var(--rule); padding-top: .8rem;
}
@media print {
  :root { --bg: #fff; --fg: #000; }
  main { max-width: none; padding: 0; }
  a { color: inherit; text-decoration: underline; }
}
"""


def _experiment_html_document(entries: "list[dict]", *,
                              title: "str | None" = None,
                              project: "str | None" = None,
                              image_srcs: "dict[str, str] | None" = None,
                              generated_at: "str | None" = None) -> str:
    """One standalone HTML page holding `entries`, newest first.

    Self-contained: the stylesheet is inlined, and `image_srcs` maps a
    stored image path to the `src` to emit (the caller passes `data:`
    URIs to make the page survive being emailed on its own).
    """
    from html import escape as _esc
    rows = [e for e in (entries or []) if isinstance(e, dict)]
    doc_title = (title or project or "").strip() or "Lab notebook"
    srcs = image_srcs or {}
    body: "list[str]" = [f"<h1>{_esc(doc_title, quote=False)}</h1>"]
    if project and len(rows) != 1:
        body.append(
            f'<p class="sc-meta">{len(rows)} '
            f'{"entry" if len(rows) == 1 else "entries"}</p>'
        )
    for e in rows:
        body.append('<article class="sc-entry">')
        e_title = (e.get("title") or "").strip() or "(untitled)"
        body.append(f"<h2>{_esc(e_title, quote=False)}</h2>")
        meta: "list[str]" = []
        if project:
            meta.append(f"<li><strong>Project:</strong> "
                        f"{_esc(str(project), quote=False)}</li>")
        for label, key in (("Entry", "id"), ("Created", "created_at"),
                           ("Updated", "updated_at")):
            val = e.get(key)
            if isinstance(val, str) and val:
                meta.append(f"<li><strong>{label}:</strong> "
                            f"{_esc(val, quote=False)}</li>")
        tags = [t for t in (e.get("tags") or []) if isinstance(t, str)]
        if tags:
            meta.append("<li><strong>Tags:</strong> "
                        + _esc(", ".join(tags), quote=False) + "</li>")
        if meta:
            body.append('<ul class="sc-meta">' + "".join(meta) + "</ul>")
        body.append(_markdown_subset_to_html(e.get("body_md") or ""))
        refs = _experiment_reference_block(e)
        if refs:
            body.append(_markdown_subset_to_html(refs))
        images = _as_str_list(e.get("image_paths"))
        if images:
            body.append("<h2>Attachments</h2>")
            for p in images:
                src = _html_safe_url(srcs.get(p, p))
                if not src:
                    continue
                if _is_image_path(p):
                    body.append(f'<p><img src="{_esc(src, quote=True)}" '
                                f'alt="{_esc(p, quote=True)}"></p>')
                else:
                    # A non-image attachment is a link, not a broken
                    # <img>. `download` names the file as stored so a
                    # `data:` URI saves under something meaningful.
                    body.append(
                        f'<p><a href="{_esc(src, quote=True)}" '
                        f'download="{_esc(p, quote=True)}" rel="noopener '
                        f'noreferrer">{_esc(p, quote=False)}</a></p>')
        body.append("</article>")
    body.append(f'<p class="sc-foot">Exported from SpliceCraft — '
                f'{_esc(generated_at or _now_iso(), quote=False)}</p>')
    return (
        "<!DOCTYPE html>\n"
        '<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, '
        'initial-scale=1">\n'
        f"<title>{_esc(doc_title, quote=False)}</title>\n"
        f"<style>{_EXPERIMENT_HTML_CSS}</style>\n"
        "</head>\n<body>\n<main>\n"
        + "\n".join(body)
        + "\n</main>\n</body>\n</html>\n"
    )
