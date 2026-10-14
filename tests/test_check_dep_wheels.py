"""scripts/check_dep_wheels.py resolves "the newest allowed version" the way a
real pip does: per Python version, honouring each release's Requires-Python.

Asked once for the interpreter running the gate, a release that dropped an old
Python (platformdirs 4.13.0, Requires-Python >=3.11) became the "newest
allowed" on 3.10 as well, had no 3.10 wheel, and blocked the v1.3.9 release —
although every 3.10 install resolves to 4.12.4, which has one (2026-10-13).
No network: pip is stubbed at the subprocess boundary.
"""
from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

import pytest
from packaging.requirements import Requirement

_ROOT = Path(__file__).resolve().parent.parent

# What PyPI says: 4.13.0 is Python >=3.11 only; 4.12.4 supports 3.10.
_VERSIONS = {"3.10": "4.12.4, 4.12.3, 4.11.11", "3.11": "4.13.0, 4.12.4, 4.11.11"}
_WHEELS = {("4.12.4", "3.10"), ("4.13.0", "3.11"), ("4.12.4", "3.11")}


@pytest.fixture
def gate(monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "check_dep_wheels_under_test", _ROOT / "scripts" / "check_dep_wheels.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    def fake_run(cmd, *a, **k):
        if cmd[3:5] == ["index", "versions"]:
            py = cmd[cmd.index("--python-version") + 1] if "--python-version" in cmd else "3.11"
            out = f"platformdirs ({_VERSIONS[py].split(',')[0]})\n" \
                  f"Available versions: {_VERSIONS[py]}\n"
            return subprocess.CompletedProcess(cmd, 0, out, "")
        if cmd[3] == "download":
            want, py = cmd[4], cmd[cmd.index("--python-version") + 1]
            if "==" in want:
                ok = (want.split("==")[1], py) in _WHEELS
            else:                           # the range: any wheel for this py
                ok = any(p == py for _v, p in _WHEELS)
            return subprocess.CompletedProcess(cmd, 0 if ok else 1, "", "")
        raise AssertionError(f"unexpected command {cmd}")

    monkeypatch.setattr(mod.subprocess, "run", fake_run)
    monkeypatch.setattr(mod, "_required_deps",
                        lambda: [Requirement("platformdirs>=4.11.11")])
    monkeypatch.setattr(mod, "PLATFORMS", mod.PLATFORMS[:1])
    monkeypatch.setattr(mod, "PYTHONS", ["3.10", "3.11"])
    return mod


def test_newest_allowed_is_asked_per_python(gate, tmp_path):
    req = Requirement("platformdirs>=4.11.11")
    assert gate._newest_allowed(req, str(tmp_path), "3.10") == "4.12.4"
    assert gate._newest_allowed(req, str(tmp_path), "3.11") == "4.13.0"


def test_a_release_that_drops_an_old_python_does_not_block(gate, capsys):
    assert gate.check() == 0
    assert "every required dependency resolves to a wheel" in capsys.readouterr().out


def test_a_release_that_drops_a_platform_wheel_still_blocks(gate, capsys):
    """The shape the pin exists for (Dependabot PR #26): the newest allowed
    release has no wheel for a platform, an older one does."""
    _WHEELS.discard(("4.13.0", "3.11"))
    try:
        assert gate.check() == 1
        assert "resolves only by backtracking" in capsys.readouterr().out
    finally:
        _WHEELS.add(("4.13.0", "3.11"))
