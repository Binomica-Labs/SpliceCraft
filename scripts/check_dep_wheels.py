#!/usr/bin/env python3
"""Release gate — every REQUIRED dependency must install from a WHEEL on
every supported platform × Python where its environment marker includes it.

Why this exists (2026-06-05): edlib sat in the required `dependencies`
with only a `sys_platform != 'win32'` marker, but edlib ships **no
Linux-aarch64 wheel**. So `pipx install splicecraft` on a 64-bit ARM
machine (Raspberry Pi, Graviton, ARM VM) had no wheel to use, fell back
to compiling edlib's C++ sdist, and failed on any box without a
toolchain — while `docs/PLATFORMS.md` promised a "clean pipx install".
A required dependency that can't install from a wheel anywhere in the
support matrix is an install-blocker; this gate makes that condition
fail the release instead of a user's machine.

Policy enforced: a dependency in `[project].dependencies` must resolve to
a wheel for every (platform, Python) its marker selects. A compiled
dependency that lacks universal wheels must either (a) carry a marker
that excludes the wheel-less platforms (and the code must degrade
gracefully there), or (b) move to an optional extra. Extras in
`[project.optional-dependencies]` are EXEMPT — they're opt-in and may
compile.

Run standalone: `python scripts/check_dep_wheels.py` (exit 1 on a gap).
Wired into `release.py` pre-flight. Needs network (queries PyPI via pip).
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

from packaging.requirements import Requirement
from packaging.version import InvalidVersion, Version

ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = ROOT / "pyproject.toml"

# Supported Python versions — from `requires-python` (the floor) up to the
# latest CPython users can already be on. ADD a new release here the moment
# it ships (distros pick it up fast: Ubuntu 26.04 shipped 3.14, which is
# exactly how the edlib-no-cp314-wheel install failure reached a user). An
# unmaintained C-ext dep silently lacks wheels for each new Python, so the
# top of this list is the early-warning tripwire.
PYTHONS = ["3.10", "3.11", "3.12", "3.13", "3.14"]

# Support matrix rows — the "Fully supported" platforms in
# docs/PLATFORMS.md. Each row: (label, sys_platform, platform_machine,
# [pip --platform tags]). The tag lists are generous (multiple
# manylinux floors + macOS OS versions + universal2) so a wheel built
# against any reasonable floor still matches. 32-bit ARM is "Limited"
# (source-compiles by design) and intentionally NOT gated here.
_MANYLINUX = lambda arch: [
    f"manylinux2014_{arch}", f"manylinux_2_17_{arch}",
    f"manylinux_2_28_{arch}", f"manylinux_2_34_{arch}",
]
_MACOS = lambda arch: [
    f"macosx_{v}_{a}"
    for v in ("10_9", "11_0", "12_0", "13_0", "14_0")
    for a in (arch, "universal2")
]
PLATFORMS = [
    ("linux-x86_64",  "linux",  "x86_64",  _MANYLINUX("x86_64")),
    ("linux-aarch64", "linux",  "aarch64", _MANYLINUX("aarch64")),
    ("macos-x86_64",  "darwin", "x86_64",  _MACOS("x86_64")),
    ("macos-arm64",   "darwin", "arm64",   _MACOS("arm64")),
    ("win-amd64",     "win32",  "AMD64",   ["win_amd64"]),
]

# Documented, ACCEPTED wheel gaps — a required dep that genuinely has no
# wheel on a platform, where the project's deliberate policy is "that
# platform compiles it from source (install a C toolchain)". Each entry
# is (dependency-name, platform-label) and MUST be mirrored by a note in
# docs/PLATFORMS.md. The gate WARNS on these but does NOT fail; it FAILS
# on any gap NOT listed here — so a dep newly dropping wheels (the
# biopython-1.87 / edlib-aarch64 class) still blocks the release, while a
# known upstream limitation we've chosen to document doesn't wedge it.
# Keep this list SHORT and justified; prefer a marker + in-code fallback
# over an accepted gap whenever the feature can degrade — but ONLY when the
# fallback returns the SAME results and only the speed differs. edlib was
# such a marker (Python <= 3.13): its alignments differed from the built-in
# engine's, the release machine could not install it, and wrong variant
# calls shipped for 16 days with only CI seeing them ([INV-212]).
ACCEPTED_GAPS = {
    # primer3-py ships no Linux-aarch64 wheel and no Apple-Silicon wheel
    # for Python >=3.10 (only cp39 arm64); no release fixes this. Primer
    # design is core with no fallback, so primer3-py stays required and
    # ARM users compile it (build-essential / Xcode CLT). See the ARM
    # rows in docs/PLATFORMS.md.
    ("primer3-py", "linux-aarch64"),
    ("primer3-py", "macos-arm64"),
}


def _marker_env(sys_platform: str, machine: str, py: str) -> dict:
    """A PEP 508 marker environment for one (platform, Python) cell."""
    return {
        "sys_platform": sys_platform,
        "platform_machine": machine,
        "python_version": py,
        "python_full_version": f"{py}.0",
        "os_name": "nt" if sys_platform == "win32" else "posix",
        "platform_system": {"linux": "Linux", "darwin": "Darwin",
                            "win32": "Windows"}[sys_platform],
        "platform_release": "",
        "implementation_name": "cpython",
        "implementation_version": f"{py}.0",
        "platform_python_implementation": "CPython",
        "extra": "",
    }


def _required_deps() -> list[Requirement]:
    data = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    return [Requirement(d) for d in data["project"]["dependencies"]]


def _wheel_exists(spec: str, py: str, plat_tags: list[str], cache: str) -> bool:
    """True iff pip can resolve `spec` (full ``name>=ver`` requirement —
    the version floor MUST be included so pip can't backtrack to an
    ancient release that happens to have a wheel) to a wheel for the
    target (python, platform). I.e. a real user's pip wouldn't have to
    build."""
    abi = "cp" + py.replace(".", "")
    cmd = [sys.executable, "-m", "pip", "download", spec,
           "--no-deps", "--only-binary=:all:",
           "--python-version", py, "--implementation", "cp", "--abi", abi,
           "-d", cache]
    for tag in plat_tags:
        cmd += ["--platform", tag]
    return subprocess.run(cmd, capture_output=True, text=True).returncode == 0


def _newest_allowed(req: Requirement, cache: str) -> "str | None":
    """The version a real ``pip install`` would actually pick for *req* — the
    highest release its specifier allows. None if that can't be determined.

    This exists because `_wheel_exists` answers a subtly different question
    than the one that matters. It runs ``pip download --only-binary=:all:``,
    which FORCES pip to find a wheel, so pip is free to backtrack to an older
    release inside the range that still has one. A real user's pip is not so
    restricted: it takes the NEWEST allowed version and, finding only an sdist
    for their platform, compiles it.

    So a cap raised past a release that DROPPED a platform's wheels reported
    clean — the exact shape of Dependabot PR #26 (biopython `<1.87` → `<1.89`,
    where 1.87 and 1.88 ship no Intel-macOS wheels but 1.86 does). Range
    coverage was green; a real Intel-Mac install would have compiled biopython
    from source. Pinning the newest allowed version is what closes that gap.

    Implemented with pip rather than a PyPI JSON call so it honours whatever
    index the caller is configured against.
    """
    r = subprocess.run(
        [sys.executable, "-m", "pip", "index", "versions", req.name,
         "--disable-pip-version-check"],
        capture_output=True, text=True, cwd=cache,
    )
    if r.returncode != 0:
        return None
    # "name (1.2.3)\n  Available versions: 1.2.3, 1.2.2, ..."
    for line in r.stdout.splitlines():
        if "Available versions:" not in line:
            continue
        raw = line.split("Available versions:", 1)[1]
        versions = []
        for chunk in raw.split(","):
            try:
                versions.append(Version(chunk.strip()))
            except InvalidVersion:
                continue
        allowed = [v for v in versions
                   if req.specifier.contains(v, prereleases=False)]
        if allowed:
            return str(max(allowed))
        return None
    return None


def check() -> int:
    deps = _required_deps()
    # entries are (name, requirement-str, platform-label, py)
    missing: list[tuple[str, str, str, str]] = []
    checked = 0
    with tempfile.TemporaryDirectory(prefix="wheelcheck-") as cache:
        for req in deps:
            # The version a real resolver picks. Checking the RANGE alone lets
            # pip backtrack to an older release that still has a wheel, which
            # hides a cap raised past a release that dropped one — see
            # `_newest_allowed`. Pin it and ask again.
            newest = _newest_allowed(req, cache)
            if newest is None:
                print(f"  note: couldn't determine the newest allowed version "
                      f"for {str(req)!r}; checked the range only.")
            for label, sysplat, machine, tags in PLATFORMS:
                for py in PYTHONS:
                    env = _marker_env(sysplat, machine, py)
                    if req.marker is not None and not req.marker.evaluate(env):
                        continue  # marker excludes this dep here — not required
                    checked += 1
                    # Ask about the PINNED newest first: it is what a real
                    # install resolves to, and if it has a wheel the range
                    # trivially does too — so the common all-green path costs
                    # one probe per cell, not two. Only when the pin fails is
                    # the range worth checking, and then only to tell a hard
                    # miss (nothing in range has a wheel) apart from the blind
                    # spot (range resolves by backtracking, newest does not).
                    if newest is not None:
                        if _wheel_exists(f"{req.name}=={newest}", py, tags,
                                         cache):
                            continue
                        spec = f"{req.name}{req.specifier}"
                        if _wheel_exists(spec, py, tags, cache):
                            missing.append(
                                (req.name,
                                 f"{req} (newest allowed: {newest} — the range "
                                 f"resolves only by backtracking to an older "
                                 f"release)", label, py))
                        else:
                            missing.append((req.name, str(req), label, py))
                        continue
                    spec = f"{req.name}{req.specifier}"  # keep the version floor
                    if not _wheel_exists(spec, py, tags, cache):
                        missing.append((req.name, str(req), label, py))

    accepted = [m for m in missing if (m[0], m[2]) in ACCEPTED_GAPS]
    blocking = [m for m in missing if (m[0], m[2]) not in ACCEPTED_GAPS]
    print(f"checked {checked} required (dependency × platform × python) cells "
          f"across {len(deps)} required deps")

    if accepted:
        # Collapse the per-Python rows to one line per (dep, platform).
        seen = sorted({(m[1], m[2]) for m in accepted})
        print("\n⚠ accepted (documented) source-build platforms — these "
              "compile from source by policy (see docs/PLATFORMS.md), NOT a "
              "regression:")
        for dep, label in seen:
            print(f"    {dep!r}  →  {label}")

    if blocking:
        print("\n✗ REQUIRED dependencies with NO wheel and NO accepted-gap "
              "entry (a clean install would be forced to compile / fail):")
        for _name, dep, label, py in blocking:
            print(f"    {dep!r}  →  {label} / cp{py.replace('.', '')}")
        print("\nFix one of:\n"
              "  • narrow the dep's marker to the platforms that HAVE wheels "
              "(and rely on an in-code fallback that gives the SAME results "
              "elsewhere), or\n"
              "  • move it to an optional extra in "
              "[project.optional-dependencies], or\n"
              "  • if the platform must compile it by deliberate policy, add "
              "(name, platform) to ACCEPTED_GAPS here AND document it in "
              "docs/PLATFORMS.md.")
        return 1

    print("\n✓ every required dependency resolves to a wheel across the whole "
          "support matrix (modulo the documented accepted-gap platforms) — no "
          "UNEXPECTED clean install can be forced to compile.")
    return 0


if __name__ == "__main__":
    raise SystemExit(check())
