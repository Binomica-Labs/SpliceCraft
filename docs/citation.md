# Citing SpliceCraft

If SpliceCraft did part of the work in something you publish — a construct, a
primer set, a codon-optimised gene, a figure — please cite it. Software
citations are how tools like this stay fundable and maintained, and they let a
reader reproduce what you did.

Every tagged release is archived on [Zenodo](https://zenodo.org/) straight from
GitHub and gets a permanent DOI, so a methods section can point at the exact
version that produced the design.

## The one-liner

```bash
splicecraft --citation
```

That prints an APA reference and a BibTeX entry **pinned to the version you are
actually running** — no need to remember which release made the plasmid:

```
Cocioba, S. (2026). SpliceCraft: a terminal-native plasmid workbench for
molecular cloning (Version 1.2.53) [Computer software]. Zenodo.
https://doi.org/10.5281/zenodo.21960400

@software{splicecraft,
  author    = {Cocioba, Sebastian},
  title     = {SpliceCraft: a terminal-native plasmid workbench for molecular cloning},
  version   = {1.2.53},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21960400},
  url       = {https://doi.org/10.5281/zenodo.21960400}
}
```

The same text is reachable from inside the app at the bottom of the `?` help
screen.

## From GitHub

The repository ships a [`CITATION.cff`](https://github.com/Binomica-Labs/SpliceCraft/blob/master/CITATION.cff)
file, so the **Cite this repository** button in the right-hand sidebar of the
[repo page](https://github.com/Binomica-Labs/SpliceCraft) exports APA or BibTeX
in one click.

## Which DOI to use

Zenodo mints two kinds of DOI, and they answer different questions:

| DOI | Resolves to | Use it when |
|---|---|---|
| **Concept DOI** | Always the newest archived release | You're citing the tool in general — a review, a methods mention, a software list |
| **Version DOI** | One specific release, forever | You're reporting work someone should be able to reproduce exactly |

`splicecraft --citation` prints the concept DOI. To cite a specific version,
open the Zenodo record and pick that version from the **Versions** list in the
sidebar — each row carries its own DOI.

## For maintainers

The citation metadata lives in three places, kept in step automatically:

- **`.zenodo.json`** — the authoritative record metadata (title, creators,
  ORCID, licence, abstract, keywords). Zenodo reads this when archiving a
  GitHub release. When both this file and `CITATION.cff` exist, Zenodo uses
  **only** `.zenodo.json`.
- **`CITATION.cff`** — what GitHub's *Cite this repository* button reads.
- **`splicecraft_util.py`** — `_ZENODO_CONCEPT_DOI` plus the `_citation_*`
  formatters behind `--citation` and the help screen.

`release.py` stamps the new version and release date into the first two on
every release, and `tests/test_citation.py` fails the suite if any of them
drift from `splicecraft.__version__` — or if the DOI string is present in some
files but not others.

### Where the DOI lives

The concept DOI is `10.5281/zenodo.21960400`, minted with v1.2.53 on
2026-08-15 — the first release archived under the Zenodo integration. It is
written down in four places, and `TestDoiIsAllOrNothing` fails the suite if any
one of them disagrees with the others:

- `splicecraft_util.py` — `_ZENODO_CONCEPT_DOI`, the single source of truth
- `CITATION.cff` — the `identifiers:` entry
- `docs/citation.md` — the worked example above
- `README.md` — the citation section

None of these need touching again. The concept DOI is stable across releases,
per-version DOIs are minted automatically by the webhook on each GitHub
Release, and the README badge resolves through the repository id rather than
any DOI string. The one thing that does move every release is the version and
date stamped into `CITATION.cff` / `.zenodo.json`, which `release.py` handles.
