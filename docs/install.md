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

## ARM64 Linux + Apple Silicon — one-time toolchain

x86-64 Linux, Intel macOS, and Windows install entirely from prebuilt
wheels — nothing to compile. On **ARM64 Linux** (64-bit Raspberry Pi,
Graviton, ARM cloud VMs) and **Apple Silicon with Python 3.10+**, one
dependency — `primer3-py`, the primer-design engine — publishes no ARM
wheel for those targets, so it compiles from source at install. Install
a C toolchain once, then install normally:

* Linux: `sudo apt install build-essential python3-dev`
* macOS: `xcode-select --install`

Then `pipx install splicecraft`. (`edlib`, the optional turbo aligner,
also lacks an ARM64-Linux wheel but transparently falls back to the
built-in pure-Python Myers aligner — no build needed; ~12× the old
Biopython fallback, identical results.)
See [`PLATFORMS.md`](PLATFORMS.md) for the full per-platform matrix.

## Troubleshooting `pipx install splicecraft`

Almost every install resolves to prebuilt wheels and just works. When it
doesn't, it's nearly always one of the cases below — each with a
one-command fix. (A release gate, `scripts/check_dep_wheels.py`, verifies
every required dependency has a wheel on every supported platform ×
Python version *before* each release, so these keep shrinking.)

### `failed to build <package>` / `Python.h: No such file or directory`

A dependency had no prebuilt wheel for your exact Python + CPU, so pip
tried to compile it from source — and the machine has no C compiler /
Python headers. Install them once, then re-run the install:

| OS | Command |
|---|---|
| Debian / Ubuntu / Raspberry Pi OS | `sudo apt install build-essential python3-dev` |
| Fedora / RHEL / Rocky | `sudo dnf install gcc gcc-c++ python3-devel` |
| Arch | `sudo pacman -S base-devel` |
| macOS | `xcode-select --install` |

Then `pipx install splicecraft` again. As of this release the only
dependency that hits this is **`primer3-py`** (the primer-design engine),
and only on **ARM64 Linux** and **Apple Silicon with Python ≥3.10**,
where it ships no wheel — it's a small, fast compile. The optional turbo
aligner **`edlib` is never the culprit**: SpliceCraft transparently falls
back to its built-in pure-Python Myers aligner wherever edlib has no wheel
(ARM, and brand-new Pythons like 3.14), so it never needs a compiler.

### `No matching distribution found for splicecraft (from versions: none)`

Your interpreter is older than **Python 3.10**. Check both
`python3 --version` **and** `pipx --version` (pipx runs on its *own*
interpreter, which may differ). Fix by pointing pipx at a 3.10+ Python:

```bash
pipx install --python python3.12 splicecraft
```

…or upgrade the system Python.

### Brand-new Python / fresh distro (e.g. Ubuntu 26.04 ships Python 3.14)

Fully supported. Most dependencies publish wheels for a new CPython
within days, and SpliceCraft's release gate checks through the newest
release. If one dependency lags, you'll see the *failed to build* case
above (install a toolchain) — or just pin an interpreter that already has
wheels for everything:

```bash
pipx install --python python3.13 splicecraft
```

### Still stuck?

pipx prints the path to the full pip log in its error (e.g.
`~/.local/state/pipx/log/cmd_*_pip_errors.log`). Open an issue with that
log plus your platform — `python3 -VV`, `uname -m` (CPU), and OS version
— and we'll add a wheel/marker fix and a gate entry so it can't recur.
See [`PLATFORMS.md`](PLATFORMS.md) for the full per-platform matrix.

## User-data location

User data (collections, library, parts, primers, features, codon
tables, settings) lives in the platform-appropriate data directory:

| Platform | Path                                          |
|----------|-----------------------------------------------|
| Linux    | `~/.local/share/splicecraft/`                 |
| macOS    | `~/Library/Application Support/splicecraft/`  |
| Windows  | `%LOCALAPPDATA%\splicecraft\`                 |

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
