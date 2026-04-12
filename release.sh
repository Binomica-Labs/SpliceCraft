#!/usr/bin/env bash
# release.sh — cut a new SpliceCraft release
#
# Usage:  ./release.sh 0.2.1
#
# Does:
#   1. Bumps the version in BOTH pyproject.toml and splicecraft.py
#   2. Runs the test suite (aborts on any failure)
#   3. Builds the package locally (aborts if twine check fails)
#   4. Commits the bump, tags v<version>, and pushes
#   5. GitHub Actions (publish.yml) then builds + uploads to PyPI
#
# Prereqs (one-time, see README):
#   - Trusted Publishing configured at pypi.org for this project
#   - git clean working tree (commit everything else first)
#
# The whole thing takes ~2 minutes from run to "live on PyPI".

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <new-version>   (e.g. $0 0.2.1)"
    exit 1
fi

NEW_VERSION="$1"

# Basic X.Y.Z sanity check
if [[ ! "$NEW_VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([a-z0-9.-]*)?$ ]]; then
    echo "Error: version must look like X.Y.Z (e.g. 0.2.1)"
    exit 1
fi

# Refuse if working tree is dirty
if [[ -n "$(git status --porcelain)" ]]; then
    echo "Error: working tree is dirty. Commit or stash first."
    git status --short
    exit 1
fi

# Refuse if tag already exists
if git rev-parse "v$NEW_VERSION" >/dev/null 2>&1; then
    echo "Error: tag v$NEW_VERSION already exists."
    exit 1
fi

echo "── Bumping version to $NEW_VERSION ──"
# pyproject.toml: line 'version         = "X.Y.Z"'
sed -i -E "s|^(version[[:space:]]*=[[:space:]]*)\"[^\"]+\"|\1\"$NEW_VERSION\"|" pyproject.toml
# splicecraft.py: '__version__ = "X.Y.Z"'
sed -i -E "s|^(__version__[[:space:]]*=[[:space:]]*)\"[^\"]+\"|\1\"$NEW_VERSION\"|" splicecraft.py

# Verify both got updated
if ! grep -q "^version\s*=\s*\"$NEW_VERSION\"" pyproject.toml; then
    echo "Error: failed to update pyproject.toml"
    exit 1
fi
if ! grep -q "^__version__\s*=\s*\"$NEW_VERSION\"" splicecraft.py; then
    echo "Error: failed to update splicecraft.py"
    exit 1
fi

echo "── Running test suite ──"
python3 -m pytest -q --tb=short

echo "── Building sdist + wheel ──"
rm -rf dist/ build/ *.egg-info
python3 -m build

echo "── Verifying package metadata ──"
python3 -m twine check dist/*

echo "── Committing + tagging + pushing ──"
git add pyproject.toml splicecraft.py
git commit -m "Release v$NEW_VERSION"
git tag "v$NEW_VERSION"
git push origin master
git push origin "v$NEW_VERSION"

echo ""
echo "═════════════════════════════════════════════════════════════"
echo " Release v$NEW_VERSION pushed."
echo " GitHub Actions will publish to PyPI in ~2 minutes."
echo " Watch:  https://github.com/Binomica-Labs/SpliceCraft/actions"
echo " Verify: https://pypi.org/project/splicecraft/"
echo "═════════════════════════════════════════════════════════════"
