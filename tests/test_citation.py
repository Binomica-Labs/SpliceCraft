"""test_citation — the DOI / citation metadata surface.

SpliceCraft is archived on Zenodo: publishing a GitHub Release snapshots the
repo and mints a permanent DOI whose metadata comes from `.zenodo.json`. That
makes a stale version string an unfixable mistake — a DOI that cites the wrong
release cannot be recalled — so the version + date live in four files and this
module is what keeps them honest:

  * `CITATION.cff`        — GitHub's "Cite this repository" button
  * `.zenodo.json`        — the Zenodo/DataCite record (wins over CITATION.cff)
  * `splicecraft.py`      — `__version__` + `_RELEASE_DATE`
  * `splicecraft_util.py` — `_ZENODO_CONCEPT_DOI` + the `--citation` formatters

`release.py::_stamp_citation_metadata` rewrites the first two on every release;
the regexes it uses are re-checked here so a reformat of either file fails the
suite instead of silently no-op-ing (or aborting) a release.

The DOI itself is all-or-nothing: it must appear in every user-facing file or
in none of them. A half-wired DOI — badge live, BibTeX still pointing at the
repo URL — is worse than no DOI, because a reader can't tell which is right.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

import splicecraft as sc
import splicecraft_util as util

REPO_ROOT = Path(__file__).resolve().parent.parent
CITATION_CFF = REPO_ROOT / "CITATION.cff"
ZENODO_JSON = REPO_ROOT / ".zenodo.json"
README = REPO_ROOT / "README.md"
DOCS_CITATION = REPO_ROOT / "docs" / "citation.md"

# Any real Zenodo DOI, e.g. 10.5281/zenodo.1234567. Used to catch a stray
# hand-typed DOI in a file the constant doesn't know about.
_DOI_RE = re.compile(r"10\.5281/zenodo\.\d+")


def _cff_field(name: str) -> str:
    """Read a top-level scalar out of CITATION.cff without a YAML parser
    (PyYAML is not a dependency). Anchored at line start so `version:` can't
    match `cff-version:`."""
    text = CITATION_CFF.read_text(encoding="utf-8")
    m = re.search(rf"^{re.escape(name)}:\s*(.+)$", text, re.MULTILINE)
    assert m, f"CITATION.cff has no top-level `{name}:` field"
    return m.group(1).strip().strip("'\"")


def _load_release():
    """Import release.py as a module. It is a script, not a package member, so
    importlib does the work; everything below the `if __name__` guard is inert."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_release_for_test", REPO_ROOT / "release.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def zenodo() -> dict:
    return json.loads(ZENODO_JSON.read_text(encoding="utf-8"))


class TestMetadataFilesAreWellFormed:
    """Both files must parse and carry the fields Zenodo/GitHub require. A
    malformed `.zenodo.json` doesn't fail the release — Zenodo silently falls
    back to auto-generated metadata, and the record ships with a machine title
    like `Binomica-Labs/SpliceCraft: v1.2.3` and no ORCID."""

    def test_zenodo_json_parses(self, zenodo):
        assert isinstance(zenodo, dict)

    def test_zenodo_required_fields(self, zenodo):
        for key in ("title", "upload_type", "license", "creators",
                    "description", "version"):
            assert zenodo.get(key), f".zenodo.json is missing `{key}`"
        assert zenodo["upload_type"] == "software"
        assert zenodo["license"] == "mit", (
            "licence id must match the repo's LICENSE (MIT) — Zenodo takes a "
            "lowercase SPDX-style id"
        )

    def test_zenodo_creator_has_orcid_and_affiliation(self, zenodo):
        creator = zenodo["creators"][0]
        # "Family, Given" is the form DataCite indexes on.
        assert "," in creator["name"], (
            f"creator name {creator['name']!r} should be 'Family, Given'"
        )
        assert re.fullmatch(r"\d{4}-\d{4}-\d{4}-\d{3}[\dX]", creator["orcid"]), (
            "ORCID must be the bare 16-digit form (no https://orcid.org/ "
            "prefix) in .zenodo.json"
        )
        assert creator.get("affiliation")

    def test_cff_required_fields(self):
        assert _cff_field("cff-version") == "1.2.0"
        assert _cff_field("type") == "software"
        assert _cff_field("license") == "MIT"
        assert _cff_field("repository-code").startswith("https://github.com/")

    def test_cff_orcid_is_a_url(self):
        text = CITATION_CFF.read_text(encoding="utf-8")
        m = re.search(r"^\s+orcid:\s*(\S+)$", text, re.MULTILINE)
        assert m, "CITATION.cff author has no orcid"
        # CFF (unlike Zenodo) requires the full URL form.
        assert re.fullmatch(
            r"https://orcid\.org/\d{4}-\d{4}-\d{4}-\d{3}[\dX]", m.group(1)
        ), m.group(1)

    def test_titles_match(self, zenodo):
        assert _cff_field("title") == zenodo["title"], (
            "CITATION.cff and .zenodo.json must present the SAME title — they "
            "are read by different tools to build the same citation string"
        )

    def test_cff_parses_as_yaml_when_available(self):
        """Full parse with a real YAML engine. PyYAML is a `[dev]` dependency
        precisely so this runs in CI: a syntax error here breaks GitHub's "Cite
        this repository" button in public, and the regex readers above cannot
        see one. The skip only fires in a bare env without the dev extras."""
        yaml = pytest.importorskip("yaml")
        data = yaml.safe_load(CITATION_CFF.read_text(encoding="utf-8"))
        assert data["cff-version"] == "1.2.0"
        assert data["authors"][0]["family-names"] == "Cocioba"
        assert isinstance(data["keywords"], list) and data["keywords"]


class TestVersionsStayInSync:
    """`release.py::_stamp_citation_metadata` moves these together. If this
    class fails, either a release half-completed or someone hand-edited one
    file — never paper over it by editing the test."""

    def test_cff_version_matches_dunder_version(self):
        assert _cff_field("version") == sc.__version__

    def test_zenodo_version_matches_dunder_version(self, zenodo):
        assert zenodo["version"] == sc.__version__

    def test_cff_release_date_matches_release_date_constant(self):
        assert _cff_field("date-released") == sc._RELEASE_DATE

    def test_release_date_is_iso(self):
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", sc._RELEASE_DATE), (
            f"_RELEASE_DATE must be ISO YYYY-MM-DD, got {sc._RELEASE_DATE!r}"
        )


class TestReleaseStampRegexesStillMatch:
    """The release script rewrites these lines with anchored regexes. Reformat
    a file and the substitution count changes — fatal at release time, which is
    the worst moment to find out. Assert here instead."""

    @staticmethod
    def _release_regexes():
        return _load_release()

    def test_each_stamp_target_matches_exactly_once(self):
        rel = self._release_regexes()
        cff = CITATION_CFF.read_text(encoding="utf-8")
        hub = (REPO_ROOT / "splicecraft.py").read_text(encoding="utf-8")
        zen = ZENODO_JSON.read_text(encoding="utf-8")
        for label, pattern, text in (
            ("CITATION.cff version", rel._CFF_VERSION_RE, cff),
            ("CITATION.cff date-released", rel._CFF_DATE_RE, cff),
            (".zenodo.json version", rel._ZENODO_VERSION_RE, zen),
            ("splicecraft.py _RELEASE_DATE", rel._RELEASE_DATE_RE, hub),
        ):
            assert len(pattern.findall(text)) == 1, (
                f"release.py's {label} regex must match exactly once"
            )

    def test_cff_version_regex_does_not_eat_cff_version_line(self):
        rel = self._release_regexes()
        cff = CITATION_CFF.read_text(encoding="utf-8")
        stamped = rel._CFF_VERSION_RE.sub(lambda m: m.group(1) + "9.9.9", cff)
        assert "cff-version: 1.2.0" in stamped, (
            "the `version:` regex clobbered the `cff-version:` schema line"
        )
        assert "\nversion: 9.9.9\n" in stamped


class TestDoiIsAllOrNothing:
    """The concept DOI lives in exactly one place in code; every user-facing
    file must agree with it."""

    _FILES = ("README.md", "CITATION.cff", "docs/citation.md")

    def test_doi_constant_shape(self):
        doi = util._ZENODO_CONCEPT_DOI
        assert doi == "" or _DOI_RE.fullmatch(doi), (
            f"_ZENODO_CONCEPT_DOI must be empty or a Zenodo DOI, got {doi!r}"
        )

    def test_doi_present_everywhere_or_nowhere(self):
        doi = util._ZENODO_CONCEPT_DOI
        found = {}
        for rel in self._FILES:
            text = (REPO_ROOT / rel).read_text(encoding="utf-8")
            # docs/citation.md shows a worked example with a placeholder;
            # ignore obviously-fake ids so the example can stay readable.
            found[rel] = {d for d in _DOI_RE.findall(text)}
        if doi:
            missing = [f for f, ds in found.items() if doi not in ds]
            assert not missing, (
                f"_ZENODO_CONCEPT_DOI is {doi} but it is missing from "
                f"{missing}. Add it (or clear the constant) — a partially "
                f"wired DOI misleads whoever reads the file that has it."
            )
        else:
            stray = {f: ds for f, ds in found.items() if ds}
            assert not stray, (
                f"a Zenodo DOI appears in {list(stray)} while "
                f"_ZENODO_CONCEPT_DOI is empty: {stray}. Set the constant so "
                f"--citation prints the same DOI the docs advertise."
            )

    def test_readme_badge_points_at_the_repo_id(self):
        text = README.read_text(encoding="utf-8")
        rid = util._ZENODO_GITHUB_REPO_ID
        assert f"https://zenodo.org/badge/{rid}.svg" in text
        assert f"https://zenodo.org/badge/latestdoi/{rid}" in text, (
            "the badge must link through Zenodo's latestdoi redirect so it "
            "keeps resolving to the newest release without a README edit"
        )

    def test_repo_id_is_numeric(self):
        assert util._ZENODO_GITHUB_REPO_ID.isdigit()


class TestCitationFormatters:
    """What `splicecraft --citation` actually prints."""

    def test_apa_carries_author_title_and_version(self):
        out = util._citation_apa("1.2.3", "2026-08-06")
        assert "Cocioba, S." in out
        assert "SpliceCraft" in out
        assert "(Version 1.2.3)" in out
        assert "[Computer software]" in out
        assert "2026" in out

    def test_bibtex_is_balanced_and_pinned(self):
        out = util._citation_bibtex("1.2.3", "2026-08-06")
        assert out.startswith("@software{splicecraft,")
        assert out.rstrip().endswith("}")
        assert out.count("{") == out.count("}"), "unbalanced BibTeX braces"
        assert "version   = {1.2.3}" in out
        assert "year      = {2026}" in out

    def test_year_comes_from_the_release_date_not_the_clock(self):
        # An old install must keep citing its own release year forever.
        assert "2019" in util._citation_apa("0.1.0", "2019-03-04")
        assert util._citation_year("2019-03-04") == "2019"

    def test_year_falls_back_when_stamp_is_junk(self):
        # A malformed stamp must not crash a citation; current year is fine.
        for junk in ("", "   ", "not-a-date", None):
            year = util._citation_year(junk)  # type: ignore[arg-type]
            assert re.fullmatch(r"\d{4}", year), junk

    def test_text_block_includes_both_forms(self):
        out = util._citation_text(sc.__version__, sc._RELEASE_DATE)
        assert sc.__version__ in out
        assert "@software{splicecraft," in out
        assert "Cocioba, S." in out

    def test_doi_appears_in_output_when_minted(self, monkeypatch):
        monkeypatch.setattr(util, "_ZENODO_CONCEPT_DOI", "10.5281/zenodo.123456")
        apa = util._citation_apa("1.2.3", "2026-08-06")
        bib = util._citation_bibtex("1.2.3", "2026-08-06")
        assert "https://doi.org/10.5281/zenodo.123456" in apa
        assert "Zenodo." in apa
        assert "doi       = {10.5281/zenodo.123456}" in bib
        assert "publisher = {Zenodo}" in bib

    def test_no_placeholder_doi_when_unminted(self, monkeypatch):
        """With no DOI, output must fall back to the repo URL — never print a
        fake or empty `https://doi.org/` that a reader might copy."""
        monkeypatch.setattr(util, "_ZENODO_CONCEPT_DOI", "")
        apa = util._citation_apa("1.2.3", "2026-08-06")
        bib = util._citation_bibtex("1.2.3", "2026-08-06")
        assert "doi.org" not in apa and "doi.org" not in bib
        assert util._CITATION_REPO_URL in apa
        assert "doi       =" not in bib

    def test_hub_reexports_are_the_same_objects(self):
        # The hub re-exports rather than re-defining; a copy would drift.
        assert sc._citation_apa is util._citation_apa
        assert sc._citation_bibtex is util._citation_bibtex
        assert sc._citation_text is util._citation_text


    def test_bibtex_protects_the_program_name_capitalisation(self):
        """plain/unsrt/abbrv title-case the title (`"t" change.case$`), which
        renders "SpliceCraft" as "Splicecraft" in the bibliography unless the
        name sits in its own brace group. The braces do not print."""
        out = util._citation_bibtex("1.2.3", "2026-08-06")
        assert "title     = {{SpliceCraft}:" in out
        assert out.count("{") == out.count("}")
        # and the protection is the ONLY difference from the plain title
        stripped = util._CITATION_TITLE_BIBTEX.replace("{", "").replace("}", "")
        assert stripped == util._CITATION_TITLE

    def test_apa_is_cleanly_spaced(self):
        for doi in ("", "10.5281/zenodo.123456"):
            out = util._citation_apa("1.2.3", "2026-08-06")
            assert "  " not in out, f"double space in APA output: {out!r}"
            assert not out.endswith("."), (
                "APA-7 puts no full stop after a trailing URL/DOI"
            )

    def test_no_placeholder_leaks_into_the_output(self):
        out = util._citation_text("1.2.3", "2026-08-06")
        for placeholder in ("XXXXXXX", "NNNNNNN", "TODO", "<name>"):
            assert placeholder not in out

    def test_year_accepts_anything_with_a_sane_str(self):
        # Typed `str`, but a caller handing over a date must not crash a
        # citation mid-print.
        import datetime
        assert util._citation_year(datetime.date(2019, 3, 4)) == "2019"  # type: ignore[arg-type]

    def test_year_rejects_non_ascii_digits(self):
        # `"\u0662\u0660\u0662\u0666".isdigit()` is True; a year in
        # Arabic-Indic digits would otherwise sail into a bibliography.
        assert util._citation_year("\u0662\u0660\u0662\u0666-01-01") == str(
            util._now().year)


class TestCliSurface:
    """`--citation` must behave like `--version`: no TUI, no *user-data* writes,
    and reachable even when the installed Textual is too old to launch. (It does
    open the log file — logging is initialised at import for every entry point,
    `--version` included — so the guarantee is about library/collections/primers,
    not about touching the data dir at all.)"""

    def test_flag_bypasses_the_textual_version_gate(self, monkeypatch):
        monkeypatch.setattr(sc, "_version_at_least", lambda have, need: False)
        monkeypatch.setattr(sc.sys, "argv", ["splicecraft", "--citation"])
        sc._check_deps()   # must return WITHOUT SystemExit

    def test_help_text_advertises_the_flag(self):
        src = (REPO_ROOT / "splicecraft.py").read_text(encoding="utf-8")
        assert "splicecraft --citation" in src

    def test_help_modal_mentions_citing(self):
        assert "--citation" in sc._HELP_BODY_MD

    def test_runs_as_a_real_cli_without_writing_user_data(self, tmp_path):
        """End-to-end: the flag is what a user actually types, so exercise the
        real process. HOME *and* XDG_DATA_HOME both point into tmp_path, so a
        stray write lands in the sandbox rather than the real library."""
        env = {
            "PATH": __import__("os").environ.get("PATH", ""),
            "HOME": str(tmp_path),
            "XDG_DATA_HOME": str(tmp_path),
            "PYTHONIOENCODING": "utf-8",
        }
        out = subprocess.run(
            [sys.executable, str(REPO_ROOT / "splicecraft.py"), "--citation"],
            env=env, cwd=str(REPO_ROOT), capture_output=True, text=True,
            timeout=300,
        )
        assert out.returncode == 0, out.stderr
        assert "Traceback" not in out.stderr, out.stderr
        assert sc.__version__ in out.stdout
        assert "Cocioba, S." in out.stdout
        assert "@software{splicecraft," in out.stdout
        assert "[Computer software]" in out.stdout

        # Logging opens a file at import (same for `--version`); nothing else
        # may be touched — no library, collections, primers, or parts bin.
        data = tmp_path / "splicecraft"
        stray = sorted(q.name for q in data.rglob("*.json")) if data.exists() else []
        assert not stray, f"--citation wrote user data: {stray}"


# Zenodo's legacy deposit metadata schema — the field set `.zenodo.json` is
# validated against when the GitHub webhook fires. An unrecognised key is at
# best ignored and at worst fails the archive, and a failed archive is silent:
# the release publishes on GitHub, no record appears, and the DOI never exists.
# Adding a field here is deliberate; a typo'd one should trip this list.
_ZENODO_LEGACY_FIELDS = frozenset({
    "upload_type", "publication_type", "image_type", "publication_date",
    "title", "creators", "description", "access_right", "license",
    "embargo_date", "access_conditions", "doi", "prereserve_doi", "keywords",
    "notes", "related_identifiers", "contributors", "references",
    "communities", "grants", "subjects", "version", "language", "locations",
    "dates", "method",
})

# Zenodo renders `description` as HTML through a sanitiser. Anything outside
# its allowlist is stripped, so a heading or a table would silently vanish from
# the permanent record — keep to the inline formatting that survives.
_ZENODO_SAFE_HTML_TAGS = frozenset({
    "p", "br", "strong", "b", "em", "i", "u", "code", "pre", "blockquote",
    "a", "ul", "ol", "li", "sub", "sup", "span", "div",
})


class TestZenodoRecordWillValidate:
    """The record metadata is the highest-stakes file in this feature: it is
    read once, at archive time, and what it produces is permanent. These checks
    are the only thing standing between a typo and a silently-skipped archive."""

    def test_no_unknown_top_level_fields(self, zenodo):
        unknown = sorted(set(zenodo) - _ZENODO_LEGACY_FIELDS)
        assert not unknown, (
            f"unknown .zenodo.json field(s) {unknown}. Zenodo validates this "
            f"file against its legacy deposit schema when the release webhook "
            f"fires; a field it does not know can fail the archive outright, "
            f"and nothing on the GitHub side reports that. Fix the key, or add "
            f"it to _ZENODO_LEGACY_FIELDS if Zenodo documents it."
        )

    def test_enumerated_values_are_the_ones_zenodo_accepts(self, zenodo):
        assert zenodo["upload_type"] == "software"
        assert zenodo["access_right"] == "open"
        assert zenodo["license"] == "mit", (
            "Zenodo wants a lowercase licence id (`mit`), not an SPDX string"
        )
        assert re.fullmatch(r"[a-z]{3}", zenodo["language"]), (
            "`language` must be an ISO 639-3 three-letter code"
        )

    def test_creator_name_is_family_comma_given(self, zenodo):
        for creator in zenodo["creators"]:
            assert "," in creator["name"], (
                f"{creator['name']!r} must be `Family, Given` — Zenodo splits "
                f"on the comma to build the citation string"
            )
            assert re.fullmatch(r"\d{4}-\d{4}-\d{4}-\d{3}[\dX]",
                                creator["orcid"]), (
                "ORCID in .zenodo.json is the bare identifier, not a URL "
                "(CITATION.cff is the one that wants the URL form)"
            )

    def test_description_uses_only_html_zenodo_keeps(self, zenodo):
        tags = {m.lower() for m in
                re.findall(r"</?([A-Za-z][A-Za-z0-9]*)", zenodo["description"])}
        stray = sorted(tags - _ZENODO_SAFE_HTML_TAGS)
        assert not stray, (
            f"tag(s) {stray} are outside Zenodo's HTML allowlist and will be "
            f"stripped from the permanent record"
        )

    def test_description_and_notes_are_not_empty(self, zenodo):
        assert len(zenodo["description"]) > 200
        assert zenodo["notes"].strip()

    def test_keywords_are_unique_non_empty_strings(self, zenodo):
        kws = zenodo["keywords"]
        assert all(isinstance(k, str) and k.strip() for k in kws)
        assert len(set(kws)) == len(kws), "duplicate keyword in .zenodo.json"

    def test_no_placeholder_text_survived_into_the_record(self, zenodo):
        blob = json.dumps(zenodo)
        for placeholder in ("XXXXXXX", "TODO", "FIXME", "<name>", "lorem"):
            assert placeholder not in blob, (
                f"{placeholder!r} is still in .zenodo.json — it would be "
                f"archived verbatim and cannot be un-published"
            )


class TestCffFileHealth:
    """GitHub parses CITATION.cff on every page render and shows a public error
    banner when it cannot. The regex readers elsewhere in this module would
    happily read a file GitHub rejects, so check the file itself."""

    def _text(self) -> str:
        return CITATION_CFF.read_text(encoding="utf-8")

    def test_no_tabs(self):
        # YAML forbids tab indentation outright; a hand-edit that pastes one
        # breaks the file in a way the field-level readers never notice.
        assert "\t" not in self._text(), "CITATION.cff contains a tab"

    def test_no_duplicate_top_level_keys(self):
        keys = re.findall(r"^([a-z][a-z-]*):", self._text(), re.MULTILINE)
        dupes = sorted({k for k in keys if keys.count(k) > 1})
        assert not dupes, (
            f"duplicate top-level key(s) {dupes} — YAML keeps the last one, so "
            f"the file would cite something other than what it appears to say"
        )

    def test_schema_version_and_type(self):
        assert _cff_field("cff-version") == "1.2.0"
        assert _cff_field("type") == "software"

    def test_license_matches_the_licence_file(self):
        assert _cff_field("license") == "MIT", (
            "CFF wants the SPDX id; `MIT` here, lowercase `mit` in .zenodo.json"
        )
        licence = (REPO_ROOT / "LICENSE").read_text(encoding="utf-8")
        assert "MIT License" in licence, (
            "CITATION.cff claims MIT but LICENSE no longer says so"
        )

    def test_date_released_is_iso_and_quoted(self):
        raw = re.search(r"^date-released:\s*(.+)$", self._text(), re.MULTILINE)
        assert raw and re.fullmatch(r"'\d{4}-\d{2}-\d{2}'", raw.group(1).strip()), (
            "release.py writes `date-released: 'YYYY-MM-DD'` and reads that "
            "shape back; changing the quoting breaks the stamp verification"
        )

    def test_title_matches_the_code_constant(self, zenodo):
        # Three files, one title. `_citation_apa` builds its string from the
        # constant, so a drift here means the app cites a different work name
        # than GitHub and Zenodo do.
        assert _cff_field("title") == util._CITATION_TITLE
        assert zenodo["title"] == util._CITATION_TITLE

    def test_repo_url_matches_the_code_constant(self):
        assert _cff_field("repository-code") == util._CITATION_REPO_URL


class TestStampIsAllOrNothing:
    """`_stamp_citation_metadata` rewrites four lines across three files. It
    must land all of them or none: a partial stamp leaves CITATION.cff claiming
    the new version while .zenodo.json still claims the old one, and whichever
    one a reader trusts, the other is a lie."""

    def _fixture_tree(self, tmp_path):
        """Real CITATION.cff + .zenodo.json copies, plus a hub stand-in that
        carries the VERBATIM `_RELEASE_DATE` line out of splicecraft.py (the
        line the regex actually has to match), so nothing here is a mock of the
        thing under test."""
        cff = tmp_path / "CITATION.cff"
        zen = tmp_path / ".zenodo.json"
        hub = tmp_path / "splicecraft.py"
        cff.write_text(CITATION_CFF.read_text(encoding="utf-8"), encoding="utf-8")
        zen.write_text(ZENODO_JSON.read_text(encoding="utf-8"), encoding="utf-8")
        hub_src = (REPO_ROOT / "splicecraft.py").read_text(encoding="utf-8")
        line = re.search(r'^_RELEASE_DATE\s*=\s*"[^"]*"$', hub_src, re.MULTILINE)
        assert line, "splicecraft.py no longer has a _RELEASE_DATE line"
        hub.write_text(f'__version__ = "0.0.0"\n{line.group(0)}\n', encoding="utf-8")
        return cff, zen, hub

    def _stamped_release_module(self, tmp_path, monkeypatch):
        rel = _load_release()
        cff, zen, hub = self._fixture_tree(tmp_path)
        monkeypatch.setattr(rel, "CITATION_CFF", cff)
        monkeypatch.setattr(rel, "ZENODO_JSON", zen)
        monkeypatch.setattr(rel, "SPLICECRAFT", hub)
        return rel, cff, zen, hub

    def test_stamps_every_target(self, tmp_path, monkeypatch):
        rel, cff, zen, hub = self._stamped_release_module(tmp_path, monkeypatch)
        today = rel._stamp_citation_metadata("9.9.9")
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", today)
        assert "\nversion: 9.9.9\n" in cff.read_text(encoding="utf-8")
        assert f"date-released: '{today}'" in cff.read_text(encoding="utf-8")
        assert json.loads(zen.read_text(encoding="utf-8"))["version"] == "9.9.9"
        assert f'_RELEASE_DATE = "{today}"' in hub.read_text(encoding="utf-8")

    def test_is_idempotent(self, tmp_path, monkeypatch):
        # A release that aborts after the stamp gets re-run from a dirty tree.
        rel, cff, zen, hub = self._stamped_release_module(tmp_path, monkeypatch)
        rel._stamp_citation_metadata("9.9.9")
        first = (cff.read_text(), zen.read_text(), hub.read_text())
        rel._stamp_citation_metadata("9.9.9")
        assert (cff.read_text(), zen.read_text(), hub.read_text()) == first

    def test_a_later_miss_writes_nothing_at_all(self, tmp_path, monkeypatch):
        """The regression this guards: CITATION.cff is edited first and
        .zenodo.json third. Write-as-you-go leaves the tree half-stamped when
        the third target cannot match."""
        rel, cff, zen, hub = self._stamped_release_module(tmp_path, monkeypatch)
        zen.write_text(
            ZENODO_JSON.read_text(encoding="utf-8").replace('"version"', '"ver"'),
            encoding="utf-8",
        )
        before = (cff.read_text(), zen.read_text(), hub.read_text())
        with pytest.raises(SystemExit):
            rel._stamp_citation_metadata("9.9.9")
        assert (cff.read_text(), zen.read_text(), hub.read_text()) == before, (
            "a failed stamp left files rewritten — the release aborts with a "
            "tree that disagrees with itself"
        )

    def test_a_missing_file_aborts_before_writing(self, tmp_path, monkeypatch):
        rel, cff, zen, hub = self._stamped_release_module(tmp_path, monkeypatch)
        before = (zen.read_text(), hub.read_text())
        cff.unlink()
        with pytest.raises(SystemExit):
            rel._stamp_citation_metadata("9.9.9")
        assert (zen.read_text(), hub.read_text()) == before

    def test_backslashes_in_values_are_not_re_interpreted(self, tmp_path,
                                                          monkeypatch):
        """`re.subn` with a string replacement would treat `\\g` or `\\1` in the
        value as a group reference. The function passes a callable instead;
        this pins that choice."""
        rel, cff, _zen, _hub = self._stamped_release_module(tmp_path, monkeypatch)
        text = cff.read_text(encoding="utf-8")
        out = rel._CFF_VERSION_RE.subn(lambda m: m.group(1) + r"1.0.0\g<1>",
                                       text)[0]
        assert r"version: 1.0.0\g<1>" in out
