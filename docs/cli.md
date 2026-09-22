# CLI sidecar (`splicecraft-cli`)

A stdlib-only Python script that connects to a running
`splicecraft --agent` session via its localhost JSON API. Imports
complete in **~50 ms** (vs ~1.5 s for the full GUI module), so an AI
agent firing dozens of commands per session doesn't pay startup cost
on every call.

Installed alongside the main `splicecraft` entry point — pipx /
pip / uv / conda all ship it.

## Usage

```bash
splicecraft --agent                  # in one terminal
splicecraft-cli status               # in another
splicecraft-cli features
splicecraft-cli fetch L09137
splicecraft-cli load-entry pUC19
splicecraft-cli add-feature 100 200 --label lacZ --type CDS --strand 1
splicecraft-cli save
```

`splicecraft-cli tools` lists every endpoint the running session
exposes (handy when wiring a new tool surface for the agent).

## Connection discovery

The CLI auto-discovers the running session via the token file at
`<DATA_DIR>/agent_token`:

| Platform | Path                                          |
|----------|-----------------------------------------------|
| Linux    | `~/.local/share/splicecraft/agent_token`      |
| macOS    | `~/Library/Application Support/splicecraft/agent_token` |
| Windows  | `%APPDATA%\splicecraft\agent_token`           |

Override the data dir with `SPLICECRAFT_DATA_DIR=/path/to/dir`.

The token file is two lines: `port\ntoken`. The CLI:

- Refuses to read a token file larger than **1 KB** (so a hostile or
  runaway process can't fill it and DoS the CLI).
- Refuses a malformed file (single line, non-integer port, empty)
  with a useful message.
- Caps the response size at **50 MB** (mirrors the server-side cap).

## Common subcommands

| Subcommand                  | Purpose                                                   |
|-----------------------------|-----------------------------------------------------------|
| `status`                    | Current record snapshot                                   |
| `tools`                     | List every endpoint + one-line doc                        |
| `features`                  | Features on the loaded record                             |
| `fetch <accession>`         | Fetch a record from NCBI by accession                     |
| `load-entry <name>`         | Load a library entry by name                              |
| `load-file <path>`          | Load a `.gb` / `.gbk` / `.dna` from a server-side path    |
| `add-feature <s> <e>`       | Add a feature to the loaded record                        |
| `save`                      | Save the loaded record (file + library)                   |
| `get-sequence <s> <e>`      | Extract DNA from a bp range (use `--bottom` for RC)       |
| `replace-sequence <s> <e> <bases>` | Replace bp range with new bases (mutagenesis)      |
| `delete-feature <idx>`      | Delete the feature at index `idx`                         |
| `update-feature <idx>`      | Update label / type / strand                              |
| `get-feature <idx>`         | Feature detail (qualifiers, etc.)                         |
| `export-genbank <path>`     | Write loaded record to PATH as GenBank                    |
| `export-fasta <path>`       | Write loaded record's seq to PATH as FASTA                |
| `list-library`              | List saved plasmid library entries                        |
| `list-collections`          | List collections + active one                             |
| `delete-from-library <name>` | Remove a library entry by name                           |
| `list-restriction-sites`    | Scan the record for restriction sites                     |
| `list-codon-tables`         | List available codon usage tables                         |
| `optimize-protein <aa>`     | Codon-optimize AA sequence to DNA                         |
| `call <endpoint>`           | Generic passthrough — call ANY endpoint (`--json '{...}'`, `--method`) |

Most subcommands support `--json` for machine-readable output and
`--force` to override unsaved-changes guards.

The named subcommands above cover the common path, but the running
session exposes ~267 endpoints (the authoritative live list is whatever the `tools` endpoint returns). `call` reaches **any** of them with the
same auth / host / port plumbing, so there's no need to hand-roll HTTP:

```bash
splicecraft-cli call rename-plasmid --json '{"old":"My Plasmid","new":"My Plasmid v2"}'
splicecraft-cli call list-entry-vectors          # GET (no body)
splicecraft-cli call blast --json '{"query":"ATGC...","collections":["MyCollection"]}'
```

Run `splicecraft-cli tools --json` first; each endpoint's `doc_full`
documents its request body.

## Authentication

The CLI attaches `Authorization: Bearer <token>` to every request.
The token is the same one written to the token file by
`splicecraft --agent` at startup — rotated each launch.

## Error mode

Errors from a named subcommand return as `Error: <message> (HTTP
<code>)` on stderr and the CLI exits non-zero. `call` instead emits the
server's structured error as JSON (with an `http_code` field) and exits
non-zero, so a calling script can parse the failure. Either way an HTTP
4xx/5xx is now a *catchable* error inside `_request` (it raises rather
than hard-`sys.exit`-ing), so importing the helper no longer risks a
4xx killing a batch mid-run. Connection refusal (server not running)
surfaces with a hint to start the GUI with `--agent`. Bad input (e.g.
`--strand 2`) fails at argparse before reaching the network.

## Why a separate sidecar?

A typical agent session calls 10–50 CLI subcommands. Importing the
full SpliceCraft module on each call pays the Textual + Biopython +
pyhmmer import cost (~1.5 s each) — turning a 5-second session into
~30 seconds of imports. The sidecar's stdlib-only design keeps each
call below ~100 ms total, dominated by HTTP roundtrip rather than
Python startup.
