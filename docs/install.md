# Install

Requires **Python 3.10+** and a terminal at least **100×30**.

## pipx (recommended)

```bash
pipx install splicecraft
splicecraft           # empty canvas
```

`pipx` installs SpliceCraft and its deps (Textual, Biopython,
primer3-py, platformdirs, pyhmmer) into an isolated virtual environment
and places the `splicecraft` command on your `PATH`. This is the right
approach on modern Debian, Ubuntu, Fedora, and WSL2, where `pip install`
into the system Python is blocked by
[PEP 668](https://peps.python.org/pep-0668/).

If you don't already have pipx:

```bash
sudo apt install pipx                # Debian / Ubuntu / WSL2
brew install pipx                    # macOS
python -m pip install --user pipx    # everywhere else
pipx ensurepath                      # one-time; adds ~/.local/bin to PATH
```

## pip inside a venv

```bash
python3 -m venv ~/.venvs/splicecraft
~/.venvs/splicecraft/bin/pip install splicecraft
~/.venvs/splicecraft/bin/splicecraft
```

Plain `pip install splicecraft` into system Python works on older
distros and inside conda envs, but is rejected by PEP 668 on any recent
Debian-family system — use `pipx` or a venv instead.

## uv tool

```bash
uv tool install splicecraft
```

Equivalent to `pipx install splicecraft` but with `uv`'s faster
resolver. Comes with the same `splicecraft update` recovery story
because the updater detects the install method at runtime.

## conda / bioconda

The bioconda recipe lives at `conda-recipe/` in the source tree. Once
the recipe lands at
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes):

```bash
conda install -c bioconda splicecraft
```

Bioconda's autotick bot handles version bumps after the first
acceptance; see [`conda-recipe/README.md`](
https://github.com/Binomica-Labs/SpliceCraft/blob/master/conda-recipe/README.md)
for the submission workflow.

## From source

```bash
git clone https://github.com/Binomica-Labs/SpliceCraft.git
cd SpliceCraft
pip install -e ".[dev]"   # inside a venv
splicecraft
```

The `[dev]` extra pulls in pytest, hypothesis, pytest-xdist, ruff,
pyright, coverage, and build — the same toolchain CI runs.

## Windows

SpliceCraft installs and runs on **native Windows** (Windows Terminal
recommended) via `pipx install splicecraft` — the editor, maps, primer
design, and BLASTN/BLASTP all work. Two Windows-specific notes:

* **Requires Python 3.10+.** If `pipx install splicecraft` fails with
  `No matching distribution found for splicecraft (from versions:
  none)`, your interpreter is older than 3.10 and pip has filtered out
  every release. Check with `python --version` (and `pipx --version`,
  which shows the interpreter pipx itself runs on), install Python
  3.10+ (`winget install Python.Python.3.12`), then retry — or pin it
  for this install with `pipx install --python 3.12 splicecraft`.
* **Local HMMscan needs WSL2.** HMMscan is powered by `pyhmmer`
  (HMMER 3), whose C core has no Windows build, so it is omitted from
  the native-Windows install and the HMMscan button explains this.
  Everything else works natively; for HMMscan, run SpliceCraft inside
  [WSL2](https://learn.microsoft.com/windows/wsl/install): `wsl
  --install`, then `pipx install splicecraft` in the Ubuntu shell.

## User-data location

User data (collections, library, parts, primers, features, codon
tables, settings) lives in the platform-appropriate data directory:

| Platform | Path                                          |
|----------|-----------------------------------------------|
| Linux    | `~/.local/share/splicecraft/`                 |
| macOS    | `~/Library/Application Support/splicecraft/`  |
| Windows  | `%APPDATA%\splicecraft\`                      |

Override with `SPLICECRAFT_DATA_DIR=/path/to/dir splicecraft`.

## Updating + recovery

```bash
splicecraft update          # detects pipx / pip / uv / pixi / conda
splicecraft update --check  # see what's available without installing
splicecraft update 0.8.10   # downgrade to a known-good version
splicecraft update --list-snapshots
splicecraft update --restore-pre-update latest
```

The updater **always snapshots your user data first** — pre-update
snapshots live in a sibling directory of the data dir (override
`SPLICECRAFT_UPDATE_BACKUP_DIR`) so a hypothetical recursive-wipe bug
in a new version cannot kill recovery. See
[Data safety](data-safety.md) for the full backup story.
