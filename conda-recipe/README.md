# Bioconda recipe for SpliceCraft

This directory holds the **reference copy** of the bioconda recipe. The
canonical copy lives in the bioconda-recipes repo on GitHub. Keeping a
copy here makes it easy to keep the two in sync and to review the recipe
alongside source changes.

## First-time submission

Bioconda review is lightweight for pure-Python, PyPI-hosted tools.
One-time steps:

```bash
# 1. Fork https://github.com/bioconda/bioconda-recipes on GitHub
# 2. Clone your fork and create a branch
git clone https://github.com/<your-user>/bioconda-recipes
cd bioconda-recipes
git checkout -b add-splicecraft

# 3. Copy the recipe into place
mkdir -p recipes/splicecraft
cp /home/seb/SpliceCraft/conda-recipe/meta.yaml recipes/splicecraft/

# 4. Regenerate source.sha256 for the current PyPI tarball
VERSION=0.2.2  # or whatever the current version is
curl -sL "https://pypi.io/packages/source/s/splicecraft/splicecraft-${VERSION}.tar.gz" \
  | sha256sum
# paste the 64-char hash into the sha256 field of meta.yaml

# 5. Lint locally (optional but fast)
# See https://bioconda.github.io/contributor/building-locally.html
bioconda-utils lint --packages splicecraft

# 6. Commit + push + open PR against bioconda/bioconda-recipes
git add recipes/splicecraft/
git commit -m "Add splicecraft recipe"
git push origin add-splicecraft
# then open a PR on GitHub
```

Bioconda's bot will build the recipe on Linux and macOS, run the
smoke tests declared in `meta.yaml` (`splicecraft --version`), and
report any linter issues.

## Subsequent releases

Once accepted, bioconda's **regro-cf-autotick-bot** detects each new
PyPI release automatically and opens a PR bumping the version +
sha256. You merely approve and merge it. Manual intervention is only
needed when dependencies or test commands change.

## Why bioconda?

Bench scientists running BLAST, Biopython, pLannotate, or any of the
thousands of bioconda-packaged tools can then do:

```bash
conda install -c bioconda splicecraft
```

…which sidesteps PEP 668 on Debian/Ubuntu/WSL2, keeps versions pinned,
and composes cleanly with the rest of a conda-based bioinformatics
workflow. pLannotate itself is on bioconda — a user can install both in
the same env.
