---
name: splicecraft
description: Drive SpliceCraft — a terminal plasmid-map viewer, sequence editor and molecular-cloning workbench — through its localhost agent API (~267 endpoints via `splicecraft-cli`). Use this whenever the work involves plasmids, vectors, constructs or a `.gb`/`.gbk`/`.dna` file: annotating a sequence, designing or checking primers, simulating a Golden Gate / Golden Braid / MoClo / Gibson / restriction clone, codon-optimizing a protein, scanning restriction sites or regulatory elements, designing CRISPR guides, verifying a build against sequencing reads, driving an Opentrons OT-2, or managing a plasmid library. Use it even when the user doesn't name SpliceCraft — if `splicecraft-cli` is on PATH it is almost certainly the right tool, and it beats hand-rolling Biopython because it simulates the actual bench steps and keeps the user's real library consistent.
---

# Driving SpliceCraft

SpliceCraft is the user's plasmid workbench: their real library, collections,
primers, parts and lab notebook live inside it. A localhost JSON API covers
essentially every action the GUI can perform, so you can do genuine
molecular-biology work rather than describing it.

**It simulates the steps, not the answer.** A cloning endpoint makes real
amplicons, runs a real digest, ligates real ends, and returns a product whose
provenance you can inspect. The result therefore carries the same failure modes
as the bench — an incompatible overhang, a regenerated site, background from a
self-closing vector. Never hand-build a construct by string concatenation when an
endpoint will simulate it; the string is what you *hoped* would happen.

**The data is the product.** These files are hours-to-years of the user's work.
The API is a *sanctioned writer* — it validates, backs up and rotates every save.
Direct file editing is not. See [Data safety](#6-data-safety).

## 1. Connect

```bash
splicecraft-cli status
```

It exits `0` when connected and `1` on every failure, and the failures are
four genuinely different situations — read the message, don't just retry:

| What it prints | What it means | Do this |
|---|---|---|
| JSON with `"ok": true` | Connected. | Start working. |
| `No SpliceCraft session found.` + `Expected token file: <path>` | No session for **that data dir**. | Either nothing is running, or a session *is* running against a different data dir. Check the path it names before starting a second session. |
| `Could not reach SpliceCraft at 127.0.0.1:<port> ([Errno 111] Connection refused).` | A token file exists but nothing is listening — a dead session's leftover. | Start a session; the stale file is harmless and gets overwritten. |
| `Malformed token file at <path> (expected `port\ntoken`).` | The file is truncated or corrupt. | Start a session to rewrite it. Don't hand-edit it. |
| `command not found` | Not installed. | Tell the user: `pipx install splicecraft`. Don't improvise a substitute. |

```bash
splicecraft --headless </dev/null &   # API only, no pty needed; default port 6701
until splicecraft-cli status >/dev/null 2>&1; do sleep 0.3; done
```

Keep the explicit `--headless`. Bare `--agent` only auto-headlesses when stdin
genuinely isn't a TTY, so `splicecraft --agent &` from an interactive shell starts
a **real TUI in the background** and stops on `SIGTTIN`.

Don't hardcode the port: if 6701 is busy the daemon walks to the next free one
(only an explicit `--agent-port` is never moved). Let `splicecraft-cli`
discover the session rather than polling a URL you guessed — it locates the
running session and handles authentication itself, so you never need to touch
the user's credentials. Use `splicecraft-cli status` as the readiness probe and
for every call after it.

**If the user already has the GUI open**, that process holds a single-instance
lock and your launch will refuse, naming the PID. Attach alongside it:

```bash
splicecraft --agent --read-only    # serves every read, 409s every write
```

That's the right move when the user is working in the app — don't ask them to
quit, and never set `SPLICECRAFT_SKIP_LOCK=1` to force in, which can corrupt the
library cache.

Use the `--read-only` **flag** and leave `SPLICECRAFT_READ_ONLY` unset. A
read-only session publishes its port to `agent_token.readonly`, and the CLI falls
back to that automatically when no host token exists — so discovery just works.
The env var is only a *token-file selector* for the CLI: exporting it pins the
CLI to the guest file, so when a normal host session is what's running you get
"No SpliceCraft session found" and nothing works. When both files exist, the host
token wins.

One diagnostic worth knowing: **`/healthz` answering `200` while `status` returns
`401`** means some unrelated process has taken the port your token file names —
not that the token needs re-reading. A clean shutdown deletes the token file, so
a leftover one usually means the previous session was killed.

Shut down a session **you** started with `splicecraft-cli call shutdown --json '{}'`
so pending writes flush. Leave the user's own session running.

## 2. Read the endpoint's docs before calling it

```bash
splicecraft-cli tools --json      # {name, method, write, doc, doc_full} × 267
splicecraft-cli call <endpoint> --json '{...}'    # POST with --json, else GET
```

`doc_full` is the complete docstring: request body, response shape, and the
specific mistake that endpoint was hardened against. It will tell you the footgun
before you hit it. Grep for the *verb* you want (`tools --json | grep -i primer`)
rather than guessing a name — the vocabulary is specific and the surface is wider
than you expect. Before concluding SpliceCraft can't do something, search.

Aliases that save a 404: `features` == `list-features`,
`plasmidsaurus-items` == `list-plasmidsaurus-items`.

**The live `tools` list beats anything written down, including this skill.** The
endpoint surface grows every release, so if an endpoint named here 404s, the
user's SpliceCraft is older than this skill — check `status.version` and work
with what `tools` actually reports rather than telling them the feature is
broken. Equally, `tools` may list endpoints this skill never mentions; that is
expected, and `doc_full` is how you learn them.

## 3. The mental model

### Canvas vs. library

One record is "loaded" (the canvas) — what `status`, `features`, `get-sequence`
and most analysis endpoints act on. The library is the persistent store of
collections. **Getting a record onto the canvas is separate from persisting it**,
and this is the most common way to lose work:

```
fetch / load-file / new-plasmid  →  canvas only, saved:false, NOT in the library
add-current-to-library           →  files it into the ACTIVE collection
save                             →  persists an EXISTING entry;
                                    creating a new one needs {"create": true}
```

A bare `save` on a freshly-fetched record returns `409` on purpose — silently
filing an inspection under its GenBank LOCUS slug would pollute the active
collection.

`plasmid-overview` is the best first look at whatever is loaded: name, length,
topology, unsaved state, features, and every enzyme that cuts exactly once
together with the feature that cut lands in.

### Coordinates: 0-based, end-exclusive, always forward-strand

The API states this in its own responses:

> `"0-based, end-exclusive (start 0 is the first base; end < start means the
> feature crosses the origin). GenBank-style 1-based inclusive = start+1..end."`

- A − strand feature still reports **forward** coordinates, so `start` is the
  lower number, *not* the biological 5′ end. On `amplify-feature` for a − feature
  it's the **reverse** primer that sits at the feature's start — the response's
  `primer_at_feature_start` tells you, and this is the concrete place strand goes
  wrong.
- `end < start` means the feature crosses the origin. Length is
  `(total - start) + end` — the response gives you `length` and `wraps`, so use
  them. **Never compute `end - start`.**
- **Never hand-roll `a <= x < b`** for "is this base inside that feature" — use
  `span-contains`. The linear form reports every base of an origin-spanning
  T-DNA, marker or operon as *outside* it: a false negative exactly where being
  wrong costs most.
- `find-orfs` takes `min_aa` in **amino acids** (`min_length`/`min_bp` are
  rejected 400); take ORF length from `nt_len`/`length_aa`.
- **The CLI's human-readable output is already 1-based**, while the JSON is not.
  The same feature prints as `2,487..2,591` from `splicecraft-cli features` and
  comes back as `start: 2486, end: 2591` over the API. Don't mix the two in one
  calculation; when in doubt take the JSON and convert once, at the end.
- When showing a position to a **human**, use 1-based inclusive
  (`start+1..end`) — what GenBank and SnapGene show — and say which you're using.

### Preview, then commit

`simulate-X` is read-only and returns what *would* happen; its write twin takes
the same payload and saves:

| Preview | Commit |
|---|---|
| `simulate-golden-gate` | `golden-gate-assemble` |
| `simulate-gibson` | `gibson-assemble` |
| `simulate-traditional-cloning` | `traditional-clone` |
| `plan-residue-edit` | `edit-residue` |

Some *write* endpoints are also dry-run by default and need `{"apply": true}`:
`transfer-annotations`, `annotate-from-presets`, `apply-gff3`,
`recover-history-from-dna`. (`commit`/`write`/`persist` are rejected 400 rather
than silently previewing.) Check the `applied` field — reporting a preview as a
completed action is a real failure mode.

## 4. Four ways a call can look successful and not be

This section is the difference between driving this tool well and being
confidently wrong with it.

### (a) HTTP 200 with `ok: false`

On `simulate-golden-gate` and `simulate-gibson`, **`ok` reports the reaction, not
the request.** A design that cannot close a circle is served correctly — HTTP
`200` — and returns `ok: false` with the reason in `result.errors`:

```
HTTP: 200   top-level ok: False
result.errors: ["part 1 released 0 fragment(s) when cut with BsaI
                (expected exactly 1 between two BsaI sites — check the part's
                 flanking sites / orientation)."]
```

An agent checking only the status code reports an unbuildable design as a
success. **Check HTTP status *and* `ok` *and*, on simulate/assembly endpoints,
the nested `result`.** The same shape appears on `bulk-import-folder`
(`partial: true` + `failures`) and the bulk `copy-plasmids`/`move-plasmids`/
`delete-primers`.

### (b) `ok: true` with a headline number that misleads

Several endpoints return a figure that is wrong on its own, with an adjacent
field that says what it really means:

| Endpoint | Headline | Sibling that qualifies it |
|---|---|---|
| `verify-against-reads` | `identity_pct` | Coverage. BLAST-style over the whole alignment, so a **perfect** 800 bp Sanger read of a 5 kb plasmid scores ~16%. A `warnings` entry says the verdict covers the sequenced windows only. |
| `read-consensus` | identity | `uncovered_spans` — "100% over the half that was read" is how a partial read passes as a full one. |
| `golden-gate-assemble`, `traditional-clone` | `ok: true` | `carried` / `carry_warnings` — omit `carry_annotations` and the product saves with **zero features**. |
| `transfer-annotations` | `count` | `skipped` — the default `min_len` of 30 drops every promoter, operator, RBS, overhang and restriction site *before* the search. |
| `optimize-protein` | `dna`, `cai` | `remaining_motifs`, `remaining_repeats`, `fixes`, `gc_window_*` — a constraint can be physically unreachable, so the response reports what was *achieved*. |
| `design-guides` | guide list | `offtarget_searched_bp` — `0` means **no off-target claim was made**. Silence is not "clean". |
| `read-consensus` vs `analyse-read-heterogeneity` | verdict | A Plasmidsaurus consensus is depth 1, so a heterogeneous population consenses to wild type and reads as a clean clonal culture. |

### (c) An empty read that means "nothing is loaded", not "nothing is there"

With **no record on the canvas**, `get-sequence` and `list-restriction-sites`
return `422 {"error": "no plasmid loaded"}` — unmistakable. But `features`
returns `200 {"ok": true, "features": []}`, which is indistinguishable from a
loaded plasmid that genuinely has no annotations. That is exactly the state a
freshly-fetched bare NCBI record is in, so the two are easy to confuse.

**Check `status.loaded` before believing an empty result.** `status` also gives
you `length` (`0` when nothing is loaded) and `dirty` in the same call, so make
it your first call in any sequence rather than inferring state from a read.

### (d) `ok: true` but it acted on the wrong thing

- **`ignored`** lists body keys the handler didn't consume. If a parameter seems
  to have had no effect, check it — you probably misspelled it. (Write endpoints
  only; reads don't echo it.)
- **These act on the ACTIVE collection only**, whatever name you pass:
  `rename-plasmid`, `set-plasmid-status`, `delete-from-library`, and
  `diff-plasmid`'s lookup. `set-active-collection` first, or use `search-library`
  to find the plasmid's home. Silent wrong-target is the hazard.
- **`?force=1` in a query string is silently ignored on writes** (POST-only, JSON
  body only). Pass `{"force": true}` in the body.
- **Batch loops:** a from-scratch `.dna` carries no name packet, so every such
  record loads as id `Cloned`; the library dedups by id, so a loop of
  `add-current-to-library` calls silently **replaces** each prior entry down to
  one survivor. Pass a unique `{id}` (scrubbed to `[A-Za-z0-9_-]`, 32 chars) and
  a spaced `{name}` per build. Note the override applies to a *copy*, so the
  canvas stays **dirty** afterwards — `discard-changes` or the next load 409s.

## 5. Making writes safe to retry

The hazard on an unattended batch: a call succeeds, the connection drops before
you see the response, and you re-send. Named creates refuse duplicates rather
than making a second copy, so what you get back is `409 "… already exists"` — an
*error* for an operation that actually worked. A batch loop reading that as
failure will abort or mis-report a build that is, in fact, fine.

Pass `--idempotency-key` on writes you may need to re-send:

```bash
splicecraft-cli call create-collection --json '{"name":"Build 07"}' \
    --idempotency-key build-07-collection
```

Within 60 s, the same key **with the same body** replays the first response —
you get the original success back, stamped `_idempotent_replay: true`, instead
of the 409. The same key with a *different* body is refused (`409`, nothing
done), which is the point: it can't silently hand you another request's answer.
Use a **fresh key per distinct request**. Keys are `[A-Za-z0-9_-]`, 128 chars
max. `409`s and `5xx`s are never cached, so genuine failures still retry
normally.

Better still, avoid the situation: the bulk endpoints — `add-features`,
`copy-plasmids`, `move-plasmids`, `delete-primers` — do N items in one locked
save, so there's one outcome to check instead of N. They also sidestep the rate
limiter, which matters because the API is throttled per session (~30 ops/s,
burst 60, **a write costs double a read**). On `429`, back off by `retry_after`.

## 6. Data safety

The library's default location is platform-dependent — don't assume the Linux
path:

| Platform | Data dir |
|---|---|
| Linux | `~/.local/share/splicecraft/` |
| macOS | `~/Library/Application Support/splicecraft/` |
| Windows | `%LOCALAPPDATA%\splicecraft\` (Local, **not** Roaming `%APPDATA%`) |

1. **Write through the API, never the filesystem.** Hand-editing
   `collections.json` or a blob bypasses validation, atomic replace and backup
   rotation. This project has lost 160 MB of real library data to exactly one
   ad-hoc script that "just needed to set up a test".
2. **To experiment, sandbox with `SPLICECRAFT_DATA_DIR`** and run a whole
   throwaway daemon — cheap, and completely isolated from the real library:
   ```bash
   mkdir -p /tmp/sc-sandbox
   export SPLICECRAFT_DATA_DIR=/tmp/sc-sandbox              # app AND cli, one value
   export SPLICECRAFT_LOG=/tmp/sc-sandbox/sandbox.log       # spare the real log
   splicecraft --headless --agent-port 6899 </dev/null &    # own data dir AND lock
   until splicecraft-cli status >/dev/null 2>&1; do sleep 0.3; done
   test -f /tmp/sc-sandbox/collections.json || echo "NOT sandboxed — stop"
   ```
   (Pick a fresh directory per run if you want a clean slate each time.)
   Use `SPLICECRAFT_DATA_DIR`, **not `XDG_DATA_HOME`**: it is honoured first, on
   every platform, sets the app and the CLI from one value, and points *directly*
   at the dir (nothing is appended). `XDG_DATA_HOME` works on Linux and macOS but
   is a **no-op on Windows**, where an "XDG sandbox" silently resolves to the
   user's real library — invisible until you have already written to it. Pass an
   absolute path with no surrounding whitespace (the app and the CLI disagree on
   trimming, and would then use different directories).

   **Verify the sandbox on the filesystem, not over the API.** No endpoint
   reports which data dir the session is using — `status` has no `data_dir`
   field — so the surest check is that the daemon created its own files in the
   directory you named, as in the snippet above (`collections.json` is written
   at launch).
3. **Never `import splicecraft` in a scratch script** against the real data dir:
   the path is fixed at import time, so it is too late to redirect afterwards.
   Prefer the API. If you genuinely need Python-level access, set the env var
   before the import and assert it took effect.
4. **Confirm before destroying.** `delete-from-library` is undoable by the user
   in the running app for that session only; a deleted *collection* is not.
   `restore-backup` overwrites user data — only on explicit request, naming which
   backup. The whole-data master wipe has no endpoint by design.

### Environment variables

| Variable | Use it? |
|---|---|
| `SPLICECRAFT_DATA_DIR` | **Yes** — the sandbox lever. Absolute path, no whitespace. |
| `SPLICECRAFT_LOG` | Yes, in a sandbox, so a probe doesn't rotate the user's real log history away. |
| `SPLICECRAFT_HEADLESS=1` | Fine — same as `--headless`. |
| `SPLICECRAFT_SKIP_LOCK=1` | **Never.** Defeats the single-instance lock; concurrent instances can corrupt the library cache. |
| `SPLICECRAFT_DEMO` | **Never.** Forces an ephemeral data dir that *ignores* your sandbox variable, and `web` disables the agent API outright. |
| `SPLICECRAFT_READ_ONLY` | Prefer the `--read-only` flag (above). |

All of these are read **once, at import**, so setting one after the process
starts changes nothing — export it before you launch.

## 7. Actions that reach outside the app

- **Physical actuation.** `ot2-run` and `ot2-position-check` move a real robot.
  Both need `{"confirm": true}` *and* a passing analysis. Never pass `confirm`
  unless the user asked for the run; walk the chain first (`ot2-compile` →
  `ot2-analyze` → `ot2-position-check` → `ot2-run`).
- **Egress.** `blast-online`, `hmmer-web`, `fpbase-search`, `read-url`,
  `learn-start` ship data off the machine and return `403` unless the user has
  armed the setting. **You cannot enable it yourself** — it's deliberately outside
  the agent allowlist. Report the refusal and name the toggle; prefer in-process
  `blast`, which never leaves the box.

## 8. Report what was actually established

The tool is scrupulous about the line between measured and guessed, and that has
to survive into your summary. Promoter and RBS scores **rank candidates, they
don't predict rates** — `scan-promoters` says so explicitly, `rbs-strength` is
relative-only, `score-guide` returns named mechanistic flags rather than a
fabricated efficiency number, and `simulate-pcr`'s binding model is exact-match,
not thermodynamic. `guide-offtargets` reports `searched_bp` so its scope can't be
overstated; there is no genome-wide search. Exact-match annotation failing means
"not identical", not "not present".

Name what was simulated, over what, and what wasn't checked. *"Three guides pass
the hygiene filters; off-targets were searched only across the 8.2 kb you
supplied"* is useful. *"Three high-efficiency guides with no off-targets"* invents
a claim the tool refused to make.

For a named reagent, strain, enzyme or antibiotic working concentration, prefer
`reference-lookup` (curated local registries) over your own recall.

## 9. When the API says no

Errors are self-explanatory and name the fix — read the message.

| Code | Meaning | Response |
|---|---|---|
| `400` `"malformed JSON body"` | **Four causes share this one message**: invalid JSON, a body over the ~1 MiB cap, a UTF-8 decode failure, or the socket breaking mid-body. Check the payload's *size* before debugging your serializer — a 1.2 MB body reports as malformed, not `413`. This is why `load-file` takes a server-side *path*, not bytes. | Shrink or chunk the body; validate the JSON separately. |
| `400` + `unsupported: [...]` | You passed a routing key (`collection`, `bin`, `enzyme`, `orientation`, `carry_annotations`, …) this endpoint ignores — it would have acted on its **default** target. Nothing ran. | Fix the key, or use an endpoint that accepts it. Never re-send without it and assume it worked. |
| `401` | Authentication failed; fires before handler lookup, so any path 401s. Credentials rotate on every launch. | Re-run through `splicecraft-cli`, which re-discovers the current session. If it persists, the session likely restarted — check `splicecraft-cli status`. |
| `403` | Egress gate unarmed, forbidden Host, or a refused write path (symlink). | You can't arm egress — tell the user which setting. |
| `404` | Unknown endpoint (body lists them all) or entity not found. | Check `tools`; use `list-library`/`search-library` for entities. |
| `405` | `GET` on a write endpoint. Writes are POST-only so a pasted URL can't fire one. | Use POST. Reads accept GET with a query string as the payload. |
| `409` dirty canvas | Unsaved edits. Note the asymmetry: `add-feature`/`add-features` are **exempt** (they're additive edits to the loaded record, so guarding them would demand `force` on every call after the first), while `load-entry`, `fetch`, `new-plasmid`, `replace-sequence` and friends do trip it. | `save`, or **`discard-changes`** (the clean recovery path — reverts and clears the undo stack, returning `{reverted}`). `{"force": true}` *in the body* discards the user's edits — ask first. |
| `409` stale canvas | The record swapped between your read and the apply step. | Re-read `status`, recompute coordinates, retry. |
| `409` ambiguous name | Exists in several collections. | Body lists holders; pass `collection`. |
| `409` "switch away first" | Deleting from a collection/bin/project while it's **active**, or copying *into* the active collection. | `set-active-*` elsewhere first. |
| `409` name collision on save | **A modal is open in front of the human**; nothing is saved until they answer. | Ask the user to answer it, or rename first. |
| `409` read-only | Attached `--read-only`; body names the lock holder's PID. | Do the read half; tell the user what needs the app closed. |
| `409` idempotency | Same key, different body. Nothing done. | Use a new key. |
| `411` | A body arrived without `Content-Length` (chunked/streaming). The body was **not read**. | Send `Content-Length`, no `Transfer-Encoding`. `splicecraft-cli` never trips this; hand-rolled HTTP does. |
| `413` | A per-endpoint input cap (sequence length, list length). | Chunk it. |
| `422` | Needs a record loaded, or the endpoint **refuses to guess**. | Information, not an obstacle. `domesticate-part` / `assemble-into-entry-vector` 422 rather than emit a stub construct you couldn't tell from a real clone. Fix the input; don't route around it. |
| `429` + `retry_after` | Rate limit (writes cost double). | Back off; prefer bulk endpoints. |
| `502` | Upstream (NCBI, EBI, Plasmidsaurus) failed. | Retry with backoff. |
| `503` + `retry_after` | Heavy-endpoint concurrency cap (half the cores). **Nothing ran.** | Wait and retry — safe. Don't fan heavy calls out in parallel. The heavy set includes several you wouldn't guess: `find-orfs`, `list-restriction-sites`, `diff-plasmid`, `fold-rna`, `lint-synthesis`, `design-rbs`. |
| `500` | Handler crash, or a save failed. | Surface it — on a save failure the data may not have persisted. |
| `_stale` on a response | The daemon is older than the version installed on disk. | Mention it; a **headless** session can `POST /restart` (a GUI `--agent` session refuses — the user restarts it). |

## Reference material

Load as needed — detail, not prerequisites:

- **`references/workflows.md`** — worked, verified end-to-end recipes: annotate a
  bare NCBI record, design and verify primers, restriction cloning, Golden Gate,
  codon optimization, read verification, notebook. Read the one matching the job.
- **`references/endpoints.md`** — all 267 endpoints organised by what you're
  trying to do, with the per-cluster gotchas worth knowing before you call.

SpliceCraft keeps its own reference docs — `docs/agent-api.md` (the contract in
prose), `docs/cli.md`, `docs/features.md` — in the project's source tree. Those
paths exist in a checkout; from a `pipx` or `pip` install they don't, so look
them up in the project's repository rather than telling the user a file is
missing.

## A note on shells

The snippets here are POSIX shell. On Windows outside WSL the shell is usually
PowerShell, which sets environment variables through its own `env:` drive rather
than `export`, has no `mktemp -d` (make a directory under the temp path
instead), and has neither `until … done` nor `</dev/null`. The *endpoint*
behaviour is identical everywhere — only the shell syntax changes, so translate
rather than reporting that something is unsupported. `splicecraft-cli` itself is
cross-platform.
