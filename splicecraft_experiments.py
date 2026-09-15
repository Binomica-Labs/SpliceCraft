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

from splicecraft_util import _NOTE_CTRL_RE, _now_iso, _sanitize_label


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
