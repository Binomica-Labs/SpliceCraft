"""release.py's confidential-name gate must read every file the release ships.

The release commit is a `git add -A`, so a brand-new, untracked file ships
with it. The curated `.private-names` gate listed files with a bare
`git ls-files` — tracked files only — so it never read a new one: a release
adding six new modules and test files scanned none of them (2026-10-13). And
`.split()` on that listing broke a file name with a space in it into pieces
that were then skipped as "not a file".
"""
from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture
def release_mod():
    spec = importlib.util.spec_from_file_location(
        "release_under_test", _ROOT / "release.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _git(cwd: Path, *args: str) -> None:
    subprocess.run(["git", "-c", "user.email=t@example.invalid",
                    "-c", "user.name=t", *args],
                   cwd=cwd, check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path, monkeypatch, release_mod) -> Path:
    root = tmp_path / "repo"
    root.mkdir()
    _git(root, "init", "-q")
    (root / ".gitignore").write_text(".private-names\nscratch/\n")
    (root / ".private-names").write_text("# comment lines are skipped\n"
                                         "democlientname\n")
    (root / "tracked.py").write_text("x = 1\n")
    _git(root, "add", ".gitignore", "tracked.py")
    _git(root, "commit", "-q", "-m", "init")
    monkeypatch.setattr(release_mod, "REPO_ROOT", root)
    monkeypatch.setattr(release_mod, "_PRIVATE_NAMES_FILE",
                        root / ".private-names")
    monkeypatch.delenv("SPLICECRAFT_SKIP_PRIVATE_CHECK", raising=False)
    return root


def test_a_new_untracked_file_is_read(release_mod, repo, capsys):
    (repo / "new_module.py").write_text("NAME = 'DemoClientName'\n")
    with pytest.raises(SystemExit):
        release_mod._check_private_tokens()
    assert "new_module.py:1" in capsys.readouterr().err


def test_a_file_name_with_a_space_is_read(release_mod, repo, capsys):
    (repo / "a note.md").write_text("about democlientname\n")
    with pytest.raises(SystemExit):
        release_mod._check_private_tokens()
    assert "a note.md:1" in capsys.readouterr().err


def test_ignored_files_are_not_read(release_mod, repo, capsys):
    (repo / "scratch").mkdir()
    (repo / "scratch" / "notes.txt").write_text("democlientname\n")
    release_mod._check_private_tokens()                 # no SystemExit
    assert "Clean" in capsys.readouterr().out
