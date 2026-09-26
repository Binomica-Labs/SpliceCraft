"""
test_notebook_tools — the lab notebook's READ side (2026-09-17).

The notebook was write-only: `attached_plasmid_ids` / `attached_actions` /
`attached_gel_ids` were rebuilt from the body on every save and read by
nothing, tags were normalised and never surfaced, and there was no way to
find an entry by its contents short of fetching every one in turn. These
tests cover the readers:

  * Search — AND-across-terms, OR-across-fields, the `#tag` filter, the
    tag-only browse, ranking, snippets, the term cap.
  * Backlinks — the reverse of `@<id>`, re-extracted from the body (so a
    stale on-disk xref can't answer "nothing references that"), matching
    a LIST of spellings (entry id + display name), case-insensitively.
  * `_iter_all_experiments` — cross-project iteration that reads the
    ACTIVE project from the live mirror, not its lagging stored copy.
  * Protocol-from-history — history step dicts → notebook markdown, and
    the two things it must never do: map an op it isn't sure about onto
    an `!action`, or emit a `@ref` for a name that wouldn't round-trip.
  * Step templates + duplication (attachments deliberately not copied).
  * Export — markdown that round-trips, and an HTML subset that escapes
    what it doesn't understand and drops `javascript:` URLs.
  * The seven agent endpoints.
"""
from __future__ import annotations

import pytest

import splicecraft as sc
import splicecraft_agent as _agent
import splicecraft_dataaccess as _da


H = sc._state._AGENT_HANDLERS


def _entry(eid, title="", body="", tags=(), updated="2026-09-01T00:00:00+00:00",
           **extra):
    e = {"id": eid, "title": title, "body_md": body, "tags": list(tags),
         "created_at": updated, "updated_at": updated}
    e.update(extra)
    return e


@pytest.fixture
def entries():
    return [
        _entry("exp-1", "Gibson attempt 3",
               "The ligation failed.\nUsed @pUC19 with !gibson.\nSee &gel-aa.",
               ["cloning", "gibson"], "2026-09-10T00:00:00+00:00",
               image_paths=["img-1.png"]),
        _entry("exp-2", "Miniprep", "Clean prep of @pUC19 today.",
               ["prep"], "2026-09-11T00:00:00+00:00"),
        _entry("exp-3", "Unrelated notes", "Nothing to see here.",
               [], "2026-09-12T00:00:00+00:00"),
    ]


# ═════════════════════════════════════════════════════════════════════════════
# Search
# ═════════════════════════════════════════════════════════════════════════════

class TestSearchTerms:
    def test_splits_and_lowercases(self):
        assert sc._experiment_search_terms("Gibson  FAILED") == \
            ["gibson", "failed"]

    def test_dedupes(self):
        assert sc._experiment_search_terms("a A a") == ["a"]

    def test_caps_terms(self):
        many = " ".join(f"t{i}" for i in range(50))
        assert len(sc._experiment_search_terms(many)) == \
            sc._EXPERIMENT_SEARCH_MAX_TERMS

    def test_non_string_is_empty(self):
        assert sc._experiment_search_terms(None) == []
        assert sc._experiment_search_terms(7) == []


class TestSearch:
    def test_and_across_terms(self, entries):
        """Both words must really appear — this is the whole reason
        search isn't the fuzzy subsequence match the plasmid picker
        uses."""
        hits = sc._experiment_search(entries, "gibson failed")
        assert [h["entry"]["id"] for h in hits] == ["exp-1"]
        assert sc._experiment_search(entries, "gibson zebrafish") == []

    def test_or_across_fields(self, entries):
        """One term in the title, one in the body still matches."""
        hits = sc._experiment_search(entries, "gibson ligation")
        assert [h["entry"]["id"] for h in hits] == ["exp-1"]
        assert set(hits[0]["fields"]) >= {"title", "body"}

    def test_tag_match_counts(self, entries):
        hits = sc._experiment_search(entries, "cloning")
        assert [h["entry"]["id"] for h in hits] == ["exp-1"]
        assert "tags" in hits[0]["fields"]

    def test_title_outranks_body(self, entries):
        """exp-2 has pUC19 in the title-adjacent sense (body) — the one
        whose TITLE matches sorts first."""
        hits = sc._experiment_search(entries, "miniprep")
        assert [h["entry"]["id"] for h in hits] == ["exp-2"]
        ranked = sc._experiment_search(entries, "puc19")
        assert {h["entry"]["id"] for h in ranked} == {"exp-1", "exp-2"}

    def test_scores_order_title_tag_body(self):
        e_title = _entry("t", "gibson", "x")
        e_tag = _entry("g", "x", "y", ["gibson"])
        e_body = _entry("b", "x", "gibson")
        hits = sc._experiment_search([e_body, e_tag, e_title], "gibson")
        assert [h["entry"]["id"] for h in hits] == ["t", "g", "b"]

    def test_tag_only_browse(self, entries):
        """Empty query + tags is a browse, not a no-op."""
        hits = sc._experiment_search(entries, "", tags=["gibson"])
        assert [h["entry"]["id"] for h in hits] == ["exp-1"]

    def test_tag_filter_is_and(self, entries):
        assert sc._experiment_search(entries, "", tags=["gibson", "prep"]) == []

    def test_tag_filter_case_insensitive(self, entries):
        assert sc._experiment_search(entries, "", tags=["GIBSON"])

    def test_tag_filter_composes_with_query(self, entries):
        assert [h["entry"]["id"] for h in
                sc._experiment_search(entries, "failed", tags=["gibson"])] \
            == ["exp-1"]
        assert sc._experiment_search(entries, "miniprep",
                                     tags=["gibson"]) == []

    def test_limit(self, entries):
        e = [_entry(f"x{i}", "hit") for i in range(5)]
        assert len(sc._experiment_search(e, "hit", limit=2)) == 2

    def test_no_criteria_matches_everything(self, entries):
        """It's a FILTER — an empty filter matches all, the same thing an
        empty filter box means in the entries pane. The endpoint is the
        layer that refuses "no criteria" with a 400."""
        assert len(sc._experiment_search(entries, "")) == len(entries)
        assert len(sc._experiment_search(entries, "", tags=[])) == len(entries)

    def test_ignores_non_dict_rows(self):
        assert sc._experiment_search(["nope", None, _entry("a", "hit")],
                                     "hit")[0]["entry"]["id"] == "a"


class TestSnippet:
    def test_context_around_first_hit(self):
        body = "x" * 200 + " the ligation failed " + "y" * 200
        s = sc._experiment_snippet(body, ["ligation"])
        assert "ligation failed" in s
        assert s.startswith("…") and s.endswith("…")

    def test_collapses_newlines(self):
        """The snippet occupies exactly one DataTable row."""
        s = sc._experiment_snippet("a\n\nb\tc ligation", ["ligation"])
        assert "\n" not in s and "\t" not in s

    def test_empty_when_body_has_no_hit(self):
        assert sc._experiment_snippet("nothing", ["gibson"]) == ""
        assert sc._experiment_snippet("", ["x"]) == ""
        assert sc._experiment_snippet(None, ["x"]) == ""

    def test_picks_earliest_term(self):
        s = sc._experiment_snippet("alpha beta", ["beta", "alpha"])
        assert s.startswith("alpha")


# ═════════════════════════════════════════════════════════════════════════════
# Backlinks
# ═════════════════════════════════════════════════════════════════════════════

class TestBacklinks:
    def test_finds_plasmid_refs(self, entries):
        hits = sc._experiments_referencing(entries, "pUC19")
        assert [e["id"] for e in hits] == ["exp-2", "exp-1"]   # newest first

    def test_case_insensitive(self, entries):
        assert sc._experiments_referencing(entries, "puc19")

    def test_kinds(self, entries):
        assert [e["id"] for e in
                sc._experiments_referencing(entries, "gel-aa", kind="gel")] \
            == ["exp-1"]
        assert [e["id"] for e in
                sc._experiments_referencing(entries, "gibson",
                                            kind="action")] == ["exp-1"]

    def test_accepts_a_list_of_spellings(self):
        """`Plasmid ref` writes the entry ID; a hand-typed ref is the
        display NAME. The backlink must find both."""
        rows = [_entry("a", "", "about @lib-0001"),
                _entry("b", "", "about @pMyPlasmid")]
        hits = sc._experiments_referencing(
            rows, ["lib-0001", "pMyPlasmid"])
        assert {e["id"] for e in hits} == {"a", "b"}

    def test_no_match_is_empty(self, entries):
        assert sc._experiments_referencing(entries, "nope") == []
        assert sc._experiments_referencing(entries, "") == []
        assert sc._experiments_referencing(entries, None) == []

    def test_unknown_kind_raises(self, entries):
        """A typo'd kind fails loudly — returning [] would read as
        "nothing references that"."""
        with pytest.raises(ValueError, match="bogus"):
            sc._experiments_referencing(entries, "x", kind="bogus")

    def test_reextracts_rather_than_trusting_stored_xref(self):
        """A hand-edited experiments.json can carry a stale index. The
        BODY is the source of truth."""
        stale = _entry("s", "", "mentions @pReal",
                       attached_plasmid_ids=["pWrong"])
        assert [e["id"] for e in
                sc._experiments_referencing([stale], "pReal")] == ["s"]
        assert sc._experiments_referencing([stale], "pWrong") == []

    def test_partial_name_does_not_match(self, entries):
        """`@pUC19` must not answer a query for `pUC`."""
        assert sc._experiments_referencing(entries, "pUC") == []


class TestIterAllExperiments:
    def test_reads_active_project_from_the_live_mirror(self, monkeypatch):
        """`_save_experiments` writes the live file first and mirrors into
        the project record after, so within a session the stored copy can
        lag. Reading the mirror for the active project is what makes a
        just-edited entry findable."""
        live = [_entry("fresh", "just saved")]
        stored_stale = [_entry("stale", "one save behind")]
        monkeypatch.setattr(_da, "_load_experiment_projects",
                            lambda: [{"name": "P1",
                                      "experiments": stored_stale},
                                     {"name": "P2",
                                      "experiments": [_entry("other")]}])
        monkeypatch.setattr(_da, "_get_active_project_name",
                            lambda: "P1")
        monkeypatch.setattr(_da, "_load_experiments",
                            lambda: live)
        rows = sc._iter_all_experiments()
        assert ("P1", live[0]) in rows
        assert not any(e["id"] == "stale" for _p, e in rows)
        assert ("P2", {"id": "other", "title": "", "body_md": "",
                       "tags": [],
                       "created_at": "2026-09-01T00:00:00+00:00",
                       "updated_at": "2026-09-01T00:00:00+00:00"}) in rows

    def test_active_pointer_with_no_project_record_still_yields_live(
            self, monkeypatch):
        """First run before `_ensure_default_project`, or a hand-edited
        settings file — the live entries must still be searchable."""
        live = [_entry("orphan")]
        monkeypatch.setattr(_da, "_load_experiment_projects",
                            lambda: [])
        monkeypatch.setattr(_da, "_get_active_project_name",
                            lambda: "Ghost")
        monkeypatch.setattr(_da, "_load_experiments",
                            lambda: live)
        assert sc._iter_all_experiments() == [("Ghost", live[0])]

    def test_skips_non_dict_entries(self, monkeypatch):
        monkeypatch.setattr(_da, "_load_experiment_projects",
                            lambda: [{"name": "P", "experiments":
                                      ["junk", None, _entry("ok")]}])
        monkeypatch.setattr(_da, "_get_active_project_name",
                            lambda: None)
        monkeypatch.setattr(_da, "_load_experiments", lambda: [])
        rows = sc._iter_all_experiments()
        assert [e["id"] for _p, e in rows] == ["ok"]


# ═════════════════════════════════════════════════════════════════════════════
# Protocol from construction history
# ═════════════════════════════════════════════════════════════════════════════

def _step(op, product, inputs=(), backbone="", enzymes=(), where=""):
    return {"op": op, "product": product, "inputs": list(inputs),
            "backbone": backbone, "enzymes": list(enzymes), "where": where}


class TestProtocolMarkdown:
    def test_numbered_steps_with_arrow(self):
        md = sc._protocol_steps_markdown(
            [_step("PCR", "amp", ["tmpl"], where="region 1-99")],
            product="pFinal")
        assert md.startswith("## Protocol — pFinal")
        assert "1. **PCR** `tmpl` → `amp` (region 1-99)" in md

    def test_maps_only_unambiguous_ops_to_action_refs(self):
        """An op the table isn't sure about keeps its verb and gets NO
        `!action` — a protocol line claiming `!ligate` for what might
        have been a Golden Gate reaction is worse than one that just
        says what was recorded."""
        md = sc._protocol_steps_markdown([
            _step("Golden Gate", "p1", ["a"], backbone="v", enzymes=["BsaI"]),
            _step("Gibson", "p2", ["b"]),
            _step("PCR", "p3", ["c"]),
            _step("mutagenesis", "p4", ["d"]),
            _step("set origin", "p5", ["p5"]),
            _step("insert site", "p6", ["e"]),
        ])
        for ref in ("!golden-gate", "!gibson", "!pcr", "!mutagenesis"):
            assert ref in md, ref
        assert "!set origin" not in md and "!set-origin" not in md
        assert "!insert site" not in md and "!insert-site" not in md

    def test_hifi_maps_to_gibson(self):
        assert "!gibson" in sc._protocol_steps_markdown(
            [_step("HiFi assembly", "p", ["a"])])

    def test_action_refs_can_be_disabled(self):
        md = sc._protocol_steps_markdown([_step("PCR", "p", ["a"])],
                                         action_refs=False)
        assert "!pcr" not in md and "**PCR**" in md

    def test_only_known_names_become_refs(self):
        """A `@ref` to something not in the library would dangle, then
        pollute `attached_plasmid_ids` and report "no such plasmid"."""
        md = sc._protocol_steps_markdown(
            [_step("", "prod", ["known", "unknown"], backbone="vec")],
            known_ids={"known", "vec"})
        assert "@known" in md and "@vec" in md
        assert "`unknown`" in md and "@unknown" not in md
        assert "`prod`" in md and "@prod" not in md

    def test_never_refs_a_name_that_would_not_round_trip(self):
        """A display name with a space would ref only its first word."""
        md = sc._protocol_steps_markdown(
            [_step("", "my plasmid", ["a b"])],
            known_ids={"my plasmid", "a b"})
        assert "@" not in md

    def test_in_place_edit_does_not_read_as_x_to_x(self):
        md = sc._protocol_steps_markdown([_step("edit", "p", ["p"])])
        assert "→" not in md and "**edit**" in md

    def test_enzymes_all_listed(self):
        md = sc._protocol_steps_markdown(
            [_step("", "p", ["a"], enzymes=["KpnI", "XbaI"])])
        assert "`KpnI` + `XbaI`" in md

    def test_no_steps_says_single_record(self):
        md = sc._protocol_steps_markdown([], product="pX")
        assert "single record" in md and "pX" in md

    def test_says_what_it_does_not_cover(self):
        """The wet-lab steps the simulation never saw must not be
        implied."""
        md = sc._protocol_steps_markdown([_step("PCR", "p", ["a"])])
        assert "transformation" in md and "not listed" in md
        for never in ("!transform", "!miniprep", "!sanger-seq"):
            assert never not in md

    def test_skips_non_dict_steps(self):
        md = sc._protocol_steps_markdown(["junk", _step("PCR", "p", ["a"])])
        assert "**PCR**" in md


class TestRefTokenOk:
    def test_accepts_plain_ids(self):
        assert sc._experiment_ref_token_ok("pUC19")
        assert sc._experiment_ref_token_ok("p.a-b_c9")

    def test_rejects_unroundtrippable(self):
        for bad in ("", None, "my plasmid", "9lead", "-lead", "a" * 100,
                    "trail-"):
            assert not sc._experiment_ref_token_ok(bad), bad

    def test_matches_the_body_extractor(self):
        """Whatever this accepts must survive a round-trip through the
        `@<id>` regex, or the protocol plants a broken ref."""
        for name in ("pUC19", "p.a-b_c9", "aB"):
            assert sc._extract_plasmid_refs(f"see @{name} here") == [name]


# ═════════════════════════════════════════════════════════════════════════════
# Step templates + duplication
# ═════════════════════════════════════════════════════════════════════════════

class TestStepTemplate:
    def test_catalog_order_not_click_order(self):
        md = sc._experiment_template_markdown(
            ["transform", "design", "pcr"], sc._EXPERIMENT_ACTIONS)
        assert md.index("Design primers") < md.index("Standard PCR") \
            < md.index("Bacterial transformation")

    def test_emits_action_refs_and_headings(self):
        md = sc._experiment_template_markdown(["miniprep"],
                                              sc._EXPERIMENT_ACTIONS)
        assert "## Plasmid miniprep  !miniprep" in md

    def test_heading_optional(self):
        assert sc._experiment_template_markdown(
            ["pcr"], sc._EXPERIMENT_ACTIONS, heading="Run 4"
        ).startswith("# Run 4")
        assert not sc._experiment_template_markdown(
            ["pcr"], sc._EXPERIMENT_ACTIONS).startswith("# ")

    def test_unknown_id_still_gets_a_section(self):
        """The catalog is curated, not enforced."""
        assert "## my-step  !my-step" in sc._experiment_template_markdown(
            ["my-step"], sc._EXPERIMENT_ACTIONS)

    def test_dedupes(self):
        md = sc._experiment_template_markdown(["pcr", "pcr"],
                                              sc._EXPERIMENT_ACTIONS)
        assert md.count("!pcr") == 1

    def test_empty_is_empty(self):
        assert sc._experiment_template_markdown([],
                                                sc._EXPERIMENT_ACTIONS) == ""

    def test_tolerates_a_malformed_catalog(self):
        md = sc._experiment_template_markdown(["x"], [("a",), None, 7])
        assert "## x  !x" in md


class TestDuplicate:
    def test_new_id_and_copy_suffix(self, entries):
        d = sc._experiment_duplicate(entries[0], new_id="exp-9")
        assert d["id"] == "exp-9"
        assert d["title"] == "Gibson attempt 3 (copy)"

    def test_explicit_title_wins(self, entries):
        assert sc._experiment_duplicate(
            entries[0], new_id="x", title="Run 4")["title"] == "Run 4"

    def test_body_and_tags_carry(self, entries):
        d = sc._experiment_duplicate(entries[0], new_id="x")
        assert d["body_md"] == entries[0]["body_md"]
        assert d["tags"] == ["cloning", "gibson"]
        assert d["attached_plasmid_ids"] == ["pUC19"]

    def test_attachments_are_not_copied(self, entries):
        """Attachment dirs are keyed by entry id — carrying the paths
        over would point the copy into the ORIGINAL's directory, so
        deleting either entry would break the other's images."""
        assert entries[0]["image_paths"] == ["img-1.png"]
        assert sc._experiment_duplicate(entries[0],
                                        new_id="x")["image_paths"] == []

    def test_timestamps_are_fresh(self, entries):
        d = sc._experiment_duplicate(entries[0], new_id="x")
        assert d["created_at"] > entries[0]["created_at"]
        assert d["created_at"] == d["updated_at"]

    def test_untitled_source(self):
        assert sc._experiment_duplicate(_entry("a"), new_id="b")["title"] \
            == "(untitled) (copy)"


# ═════════════════════════════════════════════════════════════════════════════
# Export
# ═════════════════════════════════════════════════════════════════════════════

class TestMarkdownExport:
    def test_body_is_verbatim_so_it_round_trips(self, entries):
        """An exported entry must paste back into a new one with its
        refs still working."""
        doc = sc._experiment_markdown_document(entries[0], project="Main")
        assert "@pUC19" in doc and "!gibson" in doc and "&gel-aa" in doc
        assert sc._extract_plasmid_refs(doc) == ["pUC19"]

    def test_metadata_header(self, entries):
        doc = sc._experiment_markdown_document(entries[0], project="Main")
        assert doc.startswith("# Gibson attempt 3")
        for frag in ("**Project:** Main", "`exp-1`", "**Tags:** cloning, "
                     "gibson", "**Updated:**"):
            assert frag in doc, frag

    def test_reference_block_groups_by_kind(self, entries):
        doc = sc._experiment_markdown_document(entries[0])
        assert "**Plasmids:** `@pUC19`" in doc
        assert "**Actions:** `!gibson`" in doc
        assert "**Gels:** `&gel-aa`" in doc

    def test_no_reference_block_when_no_refs(self, entries):
        assert "## References" not in \
            sc._experiment_markdown_document(entries[2])

    def test_attachments_use_supplied_srcs(self, entries):
        doc = sc._experiment_markdown_document(
            entries[0], image_srcs={"img-1.png": "data:image/png;base64,AA"})
        assert "## Attachments" in doc and "data:image/png;base64,AA" in doc

    def test_untitled_entry(self):
        assert sc._experiment_markdown_document(_entry("a")) \
            .startswith("# (untitled)")

    def test_project_document_demotes_entry_headings(self, entries):
        doc = sc._experiment_project_markdown(
            entries, project="Main", generated_at="2026-09-17")
        assert doc.startswith("# Main")
        assert "## Gibson attempt 3" in doc
        assert "3 entries" in doc
        assert doc.count("\n# ") == 0    # only ONE top-level heading

    def test_project_document_singular_noun(self, entries):
        assert "1 entry" in sc._experiment_project_markdown(entries[:1])


class TestHtmlSubset:
    def test_covers_the_documented_subset(self):
        html = sc._markdown_subset_to_html(
            "# H1\n\n## H2\n\n**b** *i* `c`\n\n- a\n- b\n\n1. x\n2. y\n\n"
            "> quoted\n\n---\n\n```py\nx = 1\n```")
        for frag in ("<h1>H1</h1>", "<h2>H2</h2>", "<strong>b</strong>",
                     "<em>i</em>", "<code>c</code>", "<ul>", "<li>a</li>",
                     "<ol>", "<blockquote>", "<hr>",
                     '<pre><code class="language-py">'):
            assert frag in html, frag

    def test_escapes_raw_html(self):
        """A body is free-form text from anywhere; an export is a file
        someone opens in a browser."""
        html = sc._markdown_subset_to_html("<script>alert(1)</script>")
        assert "<script>" not in html and "&lt;script&gt;" in html

    def test_drops_javascript_urls(self):
        # A target that parses is dropped: the label survives, the URL does not.
        html = sc._markdown_subset_to_html("[click](javascript:void0)")
        assert "javascript" not in html and "<a" not in html
        html = sc._markdown_subset_to_html("![x](javascript:void0)")
        assert "javascript" not in html and "<img" not in html
        # A target holding '(' does not parse as a link at all (the target
        # stops at '(' so a run of them can't go quadratic — audit
        # 2026-09-22): it is shown as escaped source, never as an href/src.
        for md in ("[click](javascript:alert(1))", "![x](javascript:alert(1))"):
            html = sc._markdown_subset_to_html(md)
            assert "<a" not in html and "<img" not in html
            assert 'href="javascript' not in html
            assert 'src="javascript' not in html

    def test_code_spans_protect_their_contents(self):
        assert "<strong>" not in sc._markdown_subset_to_html("`**x**`")

    def test_underscores_are_not_italics(self):
        """Lab notebooks are full of snake_case identifiers."""
        assert "<em>" not in sc._markdown_subset_to_html("pUC19_v2_final")

    def test_single_newline_is_a_hard_break(self):
        """The line breaks in a pasted colony count ARE the formatting."""
        assert "<br>" in sc._markdown_subset_to_html("a\nb")

    def test_unclosed_fence_keeps_its_content(self):
        assert "x = 1" in sc._markdown_subset_to_html("```\nx = 1")

    def test_placeholder_sentinel_cannot_be_smuggled(self):
        assert "<code>" not in sc._markdown_subset_to_html("a\x000\x00b")

    def test_switching_list_kind_closes_the_previous(self):
        html = sc._markdown_subset_to_html("- a\n1. b")
        assert html.index("</ul>") < html.index("<ol>")

    def test_non_string_input(self):
        assert sc._markdown_subset_to_html(None) == ""


class TestHtmlSafeUrl:
    @pytest.mark.parametrize("url", [
        "https://example.org/x", "http://a.b", "mailto:a@b.c", "#anchor",
        "data:image/png;base64,AA", "img-1.png", "./sub/img.png",
    ])
    def test_allows(self, url):
        assert sc._html_safe_url(url) == url

    @pytest.mark.parametrize("url", [
        "javascript:alert(1)", "JavaScript:alert(1)", "vbscript:x",
        "file:///etc/passwd", "splicecraft://plasmid/x", "data:text/html,x",
        "", "   ", None, 7,
    ])
    def test_refuses(self, url):
        assert sc._html_safe_url(url) == ""


class TestHtmlDocument:
    def test_standalone_and_print_ready(self, entries):
        html = sc._experiment_html_document(
            entries, project="Main", title="Main",
            generated_at="SEP 17 2026")
        assert html.startswith("<!DOCTYPE html>") and "</html>" in html
        assert "<title>Main</title>" in html
        assert "<style>" in html          # stylesheet inlined
        assert "prefers-color-scheme" in html and "@media print" in html
        assert "SEP 17 2026" in html

    def test_renders_each_entry(self, entries):
        html = sc._experiment_html_document(entries, project="Main")
        for e in entries:
            assert e["title"] in html

    def test_escapes_titles(self):
        html = sc._experiment_html_document(
            [_entry("a", "<script>x</script>")])
        assert "<script>x" not in html

    def test_skips_unsafe_image_srcs(self, entries):
        html = sc._experiment_html_document(
            entries, image_srcs={"img-1.png": "javascript:alert(1)"})
        assert "javascript" not in html


# ═════════════════════════════════════════════════════════════════════════════
# Agent endpoints
# ═════════════════════════════════════════════════════════════════════════════

class TestEndpointsRegistered:
    @pytest.mark.parametrize("name", [
        "search-experiments", "experiment-backlinks", "get-plasmid-protocol",
        "experiment-step-template", "list-experiment-actions",
        "export-experiment", "duplicate-experiment",
    ])
    def test_registered(self, name):
        assert name in H, f"{name} not in _AGENT_HANDLERS"

    def test_only_duplicate_is_a_write(self):
        """Everything else is a read — a read endpoint marked `write`
        would take the dirty guard for no reason."""
        assert H["duplicate-experiment"][1] is True
        for n in ("search-experiments", "experiment-backlinks",
                  "get-plasmid-protocol", "experiment-step-template",
                  "list-experiment-actions", "export-experiment"):
            assert H[n][1] is not True, n


@pytest.fixture
def notebook(monkeypatch, entries):
    """Three entries in the active project 'Main', one in 'Side'."""
    side = [_entry("exp-9", "Side note", "About @pACYC.", ["other"],
                   "2026-09-05T00:00:00+00:00")]
    store = {"live": list(entries)}

    def _save(new):
        store["live"] = list(new)
    for mod in (sc, _da, _agent):
        monkeypatch.setattr(mod, "_load_experiments",
                            lambda: list(store["live"]), raising=False)
        monkeypatch.setattr(mod, "_save_experiments", _save, raising=False)
        monkeypatch.setattr(mod, "_get_active_project_name",
                            lambda: "Main", raising=False)
        monkeypatch.setattr(
            mod, "_load_experiment_projects",
            lambda: [{"name": "Main", "experiments": list(store["live"])},
                     {"name": "Side", "experiments": side}],
            raising=False)
    return store


class TestSearchEndpoint:
    def test_active_project_by_default(self, notebook):
        body = H["search-experiments"][0](None, {"query": "gibson failed"})
        assert [h["id"] for h in body["experiments"]] == ["exp-1"]
        assert body["count"] == 1 and body["scope"] == "Main"

    def test_hit_carries_context_but_not_the_body(self, notebook):
        h = H["search-experiments"][0](
            None, {"query": "ligation"})["experiments"][0]
        assert "body_md" not in h
        assert h["body_bytes"] > 0
        assert "body" in h["matched_in"] and "ligation" in h["snippet"]
        assert h["project"] == "Main"

    def test_all_projects(self, notebook):
        body = H["search-experiments"][0](
            None, {"query": "side", "all_projects": True})
        assert [(h["project"], h["id"]) for h in body["experiments"]] == \
            [("Side", "exp-9")]
        assert body["scope"] == "all"

    def test_named_project(self, notebook):
        body = H["search-experiments"][0](
            None, {"query": "side", "project": "Side"})
        assert [h["id"] for h in body["experiments"]] == ["exp-9"]

    def test_tags_only(self, notebook):
        body = H["search-experiments"][0](None, {"tags": ["gibson"]})
        assert [h["id"] for h in body["experiments"]] == ["exp-1"]

    def test_requires_a_query_or_tags(self, notebook):
        body, status = H["search-experiments"][0](None, {})
        assert status == 400 and "query" in body["error"]

    def test_rejects_project_plus_all_projects(self, notebook):
        body, status = H["search-experiments"][0](
            None, {"query": "x", "project": "Side", "all_projects": True})
        assert status == 400 and "not both" in body["error"]

    def test_rejects_bad_types(self, notebook):
        assert H["search-experiments"][0](None, {"query": 5})[1] == 400
        assert H["search-experiments"][0](
            None, {"query": "x", "tags": "nope"})[1] == 400

    def test_unknown_project_404s(self, notebook):
        assert H["search-experiments"][0](
            None, {"query": "x", "project": "Ghost"})[1] == 404

    def test_limit_and_truncated_flag(self, notebook):
        body = H["search-experiments"][0](
            None, {"query": "e", "limit": 1})
        assert len(body["experiments"]) <= 1
        assert isinstance(body["truncated"], bool)


class TestBacklinksEndpoint:
    def test_finds_referencing_entries(self, notebook):
        body = H["experiment-backlinks"][0](None, {"ref": "pUC19"})
        assert [h["id"] for h in body["experiments"]] == ["exp-2", "exp-1"]
        assert body["kind"] == "plasmid" and body["ref"] == ["pUC19"]

    def test_accepts_a_list_of_spellings(self, notebook):
        body = H["experiment-backlinks"][0](
            None, {"ref": ["pUC19", "nope"]})
        assert body["count"] == 2

    def test_kinds(self, notebook):
        assert H["experiment-backlinks"][0](
            None, {"ref": "gel-aa", "kind": "gel"})["count"] == 1
        assert H["experiment-backlinks"][0](
            None, {"ref": "gibson", "kind": "action"})["count"] == 1

    def test_bad_kind_400s(self, notebook):
        body, status = H["experiment-backlinks"][0](
            None, {"ref": "x", "kind": "bogus"})
        assert status == 400 and "kind" in body["error"]

    def test_missing_ref_400s(self, notebook):
        for payload in ({}, {"ref": ""}, {"ref": []}, {"ref": [""]},
                        {"ref": 7}):
            assert H["experiment-backlinks"][0](None, payload)[1] == 400

    def test_no_match_is_an_empty_list_not_an_error(self, notebook):
        body = H["experiment-backlinks"][0](None, {"ref": "pNothing"})
        assert body["experiments"] == [] and body["count"] == 0

    def test_all_projects(self, notebook):
        body = H["experiment-backlinks"][0](
            None, {"ref": "pACYC", "all_projects": True})
        assert [(h["project"], h["id"]) for h in body["experiments"]] == \
            [("Side", "exp-9")]


class TestProtocolEndpoint:
    def test_requires_an_id(self):
        assert H["get-plasmid-protocol"][0](None, {})[1] == 400
        assert H["get-plasmid-protocol"][0](None, {"id": "  "})[1] == 400

    def test_unknown_plasmid_404s(self, monkeypatch):
        monkeypatch.setattr(_agent, "_agent_scan_library_for_key",
                            lambda k, c=None: [])
        body, status = H["get-plasmid-protocol"][0](None, {"id": "pNope"})
        assert status == 404 and "pNope" in body["error"]

    def test_no_history_404s(self, monkeypatch):
        monkeypatch.setattr(
            _agent, "_agent_scan_library_for_key",
            lambda k, c=None: [("C", {"name": "pX"})])
        body, status = H["get-plasmid-protocol"][0](None, {"id": "pX"})
        assert status == 404 and "history" in body["error"]

    def test_malformed_history_422s(self, monkeypatch):
        monkeypatch.setattr(
            _agent, "_agent_scan_library_for_key",
            lambda k, c=None: [("C", {"name": "pX", "history_xml": "<x>"})])

        def _boom(_xml):
            raise ValueError("bad xml")
        monkeypatch.setattr(_agent, "_parse_commercialsaas_history", _boom)
        body, status = H["get-plasmid-protocol"][0](None, {"id": "pX"})
        assert status == 422 and "bad xml" in body["error"]

    def test_returns_steps_and_markdown(self, monkeypatch):
        monkeypatch.setattr(
            _agent, "_agent_scan_library_for_key",
            lambda k, c=None: ([("C", {"name": "pX",
                                       "history_xml": "<x>"})]
                               if k in ("pX", "tmpl") else []))
        monkeypatch.setattr(_agent, "_parse_commercialsaas_history",
                            lambda _x: object())
        monkeypatch.setattr(
            _agent, "_history_build_steps",
            lambda _r: [_step("PCR", "amp", ["tmpl"], where="region 1-9")])
        body = H["get-plasmid-protocol"][0](None, {"id": "pX"})
        assert body["plasmid"] == "pX" and body["n_steps"] == 1
        assert body["steps"][0]["op"] == "PCR"
        assert body["steps"][0]["inputs"] == ["tmpl"]
        assert "**PCR**" in body["markdown"] and "!pcr" in body["markdown"]
        # `tmpl` resolves in the library, so it may be a ref; `amp` does
        # not, so it must not be.
        assert "@amp" not in body["markdown"]

    def test_empty_history_404s(self, monkeypatch):
        monkeypatch.setattr(
            _agent, "_agent_scan_library_for_key",
            lambda k, c=None: [("C", {"name": "pX", "history_xml": "<x>"})])
        monkeypatch.setattr(_agent, "_parse_commercialsaas_history",
                            lambda _x: None)
        assert H["get-plasmid-protocol"][0](None, {"id": "pX"})[1] == 404

    def test_rejects_bad_collection_type(self):
        assert H["get-plasmid-protocol"][0](
            None, {"id": "pX", "collection": 7})[1] == 400


class TestStepTemplateEndpoint:
    def test_returns_markdown_in_catalog_order(self):
        body = H["experiment-step-template"][0](
            None, {"actions": ["transform", "design"]})
        assert body["markdown"].index("Design primers") < \
            body["markdown"].index("Bacterial transformation")
        assert body["unknown_actions"] == []

    def test_names_unknown_actions_without_refusing(self):
        body = H["experiment-step-template"][0](
            None, {"actions": ["pcr", "my-step"]})
        assert body["unknown_actions"] == ["my-step"]
        assert "!my-step" in body["markdown"]

    def test_heading(self):
        body = H["experiment-step-template"][0](
            None, {"actions": ["pcr"], "heading": "Run 4"})
        assert body["markdown"].startswith("# Run 4")

    def test_validates(self):
        for payload in ({}, {"actions": []}, {"actions": "pcr"},
                        {"actions": [""]}):
            assert H["experiment-step-template"][0](None, payload)[1] == 400
        assert H["experiment-step-template"][0](
            None, {"actions": ["pcr"], "heading": 7})[1] == 400


class TestActionsCatalogEndpoint:
    def test_lists_the_catalog(self):
        body = H["list-experiment-actions"][0](None, {})
        assert body["count"] == len(sc._EXPERIMENT_ACTIONS) == 19
        first = body["actions"][0]
        assert set(first) == {"group", "id", "description"}
        assert {a["id"] for a in body["actions"]} >= {
            "pcr", "golden-gate", "gibson", "miniprep", "transform"}


class TestExportEndpoint:
    def test_single_entry_markdown(self, notebook):
        body = H["export-experiment"][0](None, {"id": "exp-1"})
        assert body["format"] == "markdown" and body["n_entries"] == 1
        assert body["document"].startswith("# Gibson attempt 3")
        assert body["bytes"] == len(body["document"].encode("utf-8"))

    def test_whole_project_markdown(self, notebook):
        body = H["export-experiment"][0](None, {})
        assert body["n_entries"] == 3
        assert body["document"].startswith("# Main")

    def test_html(self, notebook):
        body = H["export-experiment"][0](
            None, {"id": "exp-1", "format": "html"})
        assert body["format"] == "html"
        assert body["document"].startswith("<!DOCTYPE html>")

    def test_md_alias(self, notebook):
        assert H["export-experiment"][0](
            None, {"id": "exp-1", "format": "md"})["format"] == "markdown"

    def test_bad_format_400s(self, notebook):
        body, status = H["export-experiment"][0](
            None, {"format": "pdf"})
        assert status == 400 and "format" in body["error"]

    def test_unknown_id_404s(self, notebook):
        assert H["export-experiment"][0](None, {"id": "exp-nope"})[1] == 404

    def test_invalid_id_400s(self, notebook):
        assert H["export-experiment"][0](None, {"id": "../x"})[1] == 400

    def test_named_project(self, notebook):
        body = H["export-experiment"][0](None, {"project": "Side"})
        assert body["n_entries"] == 1 and body["project"] == "Side"


class TestDuplicateEndpoint:
    def test_duplicates_into_the_live_project(self, notebook):
        body = H["duplicate-experiment"][0](None, {"id": "exp-1"})
        assert body["ok"] is True and body["source_id"] == "exp-1"
        assert body["title"] == "Gibson attempt 3 (copy)"
        assert body["n_entries"] == 4
        assert any(e["id"] == body["id"] for e in notebook["live"])

    def test_reports_that_attachments_did_not_come_along(self, notebook):
        body = H["duplicate-experiment"][0](None, {"id": "exp-1"})
        assert body["attachments_copied"] is False
        copy = next(e for e in notebook["live"] if e["id"] == body["id"])
        assert copy["image_paths"] == []

    def test_explicit_title(self, notebook):
        body = H["duplicate-experiment"][0](
            None, {"id": "exp-1", "title": "Run 4"})
        assert body["title"] == "Run 4"

    def test_validates(self, notebook):
        assert H["duplicate-experiment"][0](None, {})[1] == 400
        assert H["duplicate-experiment"][0](None, {"id": "../x"})[1] == 400
        assert H["duplicate-experiment"][0](
            None, {"id": "exp-1", "title": 7})[1] == 400

    def test_unknown_id_404s(self, notebook):
        assert H["duplicate-experiment"][0](None, {"id": "exp-nope"})[1] == 404

    def test_refuses_a_project_key_rather_than_guessing(self, notebook):
        """It both reads AND writes, so pointing it at a non-active
        project would desync the live mirror ([INV-161]). Refuse loudly
        instead of duplicating whatever shares that id in the active
        project."""
        body, status = H["duplicate-experiment"][0](
            None, {"id": "exp-1", "project": "Side"})
        assert status == 400 and "ACTIVE project only" in body["error"]
        assert len(notebook["live"]) == 3      # nothing written


# ═════════════════════════════════════════════════════════════════════════════
# Live mount — the widgets really compose, and the flows really run
# ═════════════════════════════════════════════════════════════════════════════

_TERM = (160, 48)


class TestNotebookUI:
    """Mounts the real screen. `#exp-compose-btns` became a
    `HorizontalScroll` and the entries pane gained a filter row + a
    fourth button row, so a CSS or compose slip here is a hard crash
    rather than a wrong number."""

    async def test_new_widgets_compose(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            assert isinstance(scr, sc.ExperimentsScreen)
            for sel in ("#exp-filter-input", "#exp-filter-count",
                        "#btn-exp-find", "#btn-exp-duplicate",
                        "#btn-exp-steps", "#btn-exp-protocol",
                        "#btn-exp-export"):
                assert scr.query_one(sel) is not None, sel

    async def test_filter_narrows_and_counts(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.action_new_entry()
            await pilot.pause()
            scr.query_one("#exp-title-input").value = "Gibson run"
            scr.query_one("#exp-tags-input").value = "cloning"
            scr.query_one("#exp-body").text = "The ligation failed."
            assert scr._persist_current()
            await pilot.pause()
            table = scr.query_one("#exp-entries-table")
            count = scr.query_one("#exp-filter-count")

            scr.query_one("#exp-filter-input").value = "gibson"
            scr._refresh_entries_table()
            assert table.row_count == 1
            assert "1 of 1" in str(count.render())

            scr.query_one("#exp-filter-input").value = "ligation gibson"
            scr._refresh_entries_table()
            assert table.row_count == 1      # AND across title + body

            scr.query_one("#exp-filter-input").value = "zebrafish"
            scr._refresh_entries_table()
            assert table.row_count == 0
            # An over-narrow filter must read as a filter, not an empty
            # notebook.
            assert "no match" in str(count.render())

            scr.query_one("#exp-filter-input").value = "#cloning"
            scr._refresh_entries_table()
            assert table.row_count == 1

            scr.query_one("#exp-filter-input").value = "#nosuchtag"
            scr._refresh_entries_table()
            assert table.row_count == 0

    async def test_row_keys_survive_filtering(self):
        """Selection resolves the row KEY, so a filtered table still
        opens the right entry (sacred #33 sort/lookup symmetry)."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([
                _entry("exp-aa", "Gibson run", "ligation failed"),
                _entry("exp-bb", "Miniprep", "clean"),
            ])
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#exp-filter-input").value = "miniprep"
            scr._refresh_entries_table()
            await pilot.pause()
            table = scr.query_one("#exp-entries-table")
            assert table.row_count == 1
            table.move_cursor(row=0)
            scr._on_btn_open(None)
            await pilot.pause()
            assert scr._current_entry is not None
            assert scr._current_entry["id"] == "exp-bb"

    async def test_duplicate_flow(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([_entry("exp-aa", "Gibson run", "notes",
                                         image_paths=["img-1.png"])])
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#exp-entries-table").move_cursor(row=0)
            scr._duplicate_current()
            await pilot.pause()
            await pilot.pause()
            rows = sc._load_experiments()
            assert len(rows) == 2
            copy = next(e for e in rows if e["title"].endswith("(copy)"))
            assert copy["image_paths"] == []      # not carried over
            assert scr._current_entry["id"] == copy["id"]

    async def test_steps_template_inserts_in_catalog_order(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.action_new_entry()
            await pilot.pause()
            scr.action_insert_steps()
            await pilot.pause()
            await pilot.pause()
            picker = app.screen
            assert isinstance(picker, sc.ActionsPickerModal)
            picker._picked = {"miniprep", "pcr"}
            picker._ok(None)
            await pilot.pause()
            await pilot.pause()
            body = scr.query_one("#exp-body").text
            assert "!pcr" in body and "!miniprep" in body
            assert body.index("Standard PCR") < body.index("Plasmid miniprep")

    async def test_actions_picker_multi_mode(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.ActionsPickerModal(multi=True))
            await pilot.pause()
            await pilot.pause()
            m = app.screen
            table = m.query_one("#actpick-table")
            assert len(table.columns) == 4       # pick column added
            ok = m.query_one("#btn-actpick-ok")
            assert ok.disabled is True           # nothing picked yet
            table.move_cursor(row=0)
            m.action_toggle_pick()
            await pilot.pause()
            assert ok.disabled is False and "(1)" in str(ok.label)
            # The label must also FIT: a 10-column button clipped
            # "Insert (1)" back to exactly "Insert", so the count was
            # set and invisible.
            assert ok.size.width >= len(str(ok.label)), \
                (ok.size.width, str(ok.label))
            m.action_toggle_pick()               # toggles off
            await pilot.pause()
            assert m._picked == set() and ok.disabled is True

    async def test_single_mode_picker_unchanged(self):
        """The existing single-select contract still returns a str."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            got = {}
            app.push_screen(sc.ActionsPickerModal(),
                             callback=lambda p: got.update(v=p))
            await pilot.pause()
            await pilot.pause()
            m = app.screen
            m.query_one("#actpick-table").move_cursor(row=0)
            m._ok(None)
            await pilot.pause()
            await pilot.pause()
            assert got.get("v") == "design"

    async def test_search_modal_finds_and_backlinks(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([
                _entry("exp-aa", "Gibson run",
                       "ligation failed on @pUC19", ["cloning"]),
                _entry("exp-bb", "Miniprep", "clean"),
            ])
            app.push_screen(sc.ExperimentSearchModal())
            await pilot.pause()
            await pilot.pause()
            m = app.screen
            assert isinstance(m, sc.ExperimentSearchModal)
            m.query_one("#expsearch-input").value = "ligation"
            m._refresh()
            await pilot.pause()
            assert [e["id"] for _p, e in m._matches] == ["exp-aa"]
            assert m.query_one("#expsearch-table").row_count == 1
            m.query_one("#expsearch-input").value = "#cloning"
            m._refresh()
            await pilot.pause()
            assert [e["id"] for _p, e in m._matches] == ["exp-aa"]

    async def test_search_modal_backlink_mode(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([
                _entry("exp-aa", "Gibson run", "about @pUC19"),
                _entry("exp-bb", "Other", "nothing"),
            ])
            app.push_screen(sc.ExperimentSearchModal(
                backlink=("plasmid", ["pUC19", "lib-1"])))
            await pilot.pause()
            await pilot.pause()
            m = app.screen
            assert [e["id"] for _p, e in m._matches] == ["exp-aa"]
            assert "pUC19" in m._title_text()

    async def test_export_modal_writes_both_formats(self, tmp_path):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            entry = _entry("exp-aa", "Gibson run",
                           "# Notes\n\nligation failed on @pUC19",
                           ["cloning"])
            sc._save_experiments([entry])
            for fname, fmt, head in (("out.md", "md", "# Gibson run"),
                                      ("out.html", "html", "<!DOCTYPE html>")):
                m = sc.ExperimentExportModal(entry=entry, project="Main")
                app.push_screen(m)
                await pilot.pause()
                await pilot.pause()
                m._selected_dir = str(tmp_path)
                m.query_one("#expexp-fmt").value = fmt
                m.query_one("#expexp-filename").value = fname
                out = tmp_path / fname
                m._do_export()
                # The write runs on a worker thread, and the modal
                # dismisses itself when it finishes — which cancels the
                # worker, so `wait_for_complete()` would raise. Poll the
                # file instead.
                for _ in range(200):
                    if out.is_file():
                        break
                    await pilot.pause()
                assert out.is_file(), fname
                text = out.read_text(encoding="utf-8")
                assert text.startswith(head), (fname, text[:40])
                assert "@pUC19" in text        # refs survive both formats
                if app.screen is m:
                    app.pop_screen()
                    await pilot.pause()

    async def test_export_modal_refuses_a_mismatched_extension(self, tmp_path):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            entry = _entry("exp-aa", "Gibson run", "notes")
            sc._save_experiments([entry])
            m = sc.ExperimentExportModal(entry=entry, project="Main")
            app.push_screen(m)
            await pilot.pause()
            await pilot.pause()
            m._selected_dir = str(tmp_path)
            m.query_one("#expexp-fmt").value = "html"
            m.query_one("#expexp-filename").value = "out.md"
            m._do_export()
            await pilot.pause()
            assert "must end in .html" in str(
                m.query_one("#expexp-status").render())
            assert not (tmp_path / "out.md").exists()

    async def test_library_row_n_opens_backlinks(self):
        """`n` on a library row posts the message the app turns into a
        backlink view. The App also binds `n` (find-next) but without
        priority, so the panel wins while the library has focus."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            panel = app.query_one(sc.LibraryPanel)
            assert any(b.key == "n" and b.action == "request_notebook_refs"
                       for b in panel.BINDINGS)
            app._library_notebook_refs_requested(
                sc.LibraryPanel.NotebookRefsRequested(None))
            await pilot.pause()
            # No row highlighted → a warning, not a crash or a modal.
            assert not isinstance(app.screen, sc.ExperimentSearchModal)

    async def test_filter_debounce_is_cancelled_on_unmount(self):
        """A timer that outlives the screen refreshes unmounted widgets
        — the v1.2.48 post-unmount re-arm class."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#exp-filter-input").value = "x"
            scr._on_filter_changed(None)
            assert scr._filter_timer is not None
            app.pop_screen()
            await pilot.pause()
            assert scr._filter_timer is None


# ═════════════════════════════════════════════════════════════════════════════
# Attachments past images
# ═════════════════════════════════════════════════════════════════════════════

_TINY_PNG = bytes.fromhex(
    "89504e470d0a1a0a0000000d494844520000000100000001080600000"
    "01f15c4890000000a49444154789c6300010000050001"
    "0d0a2db40000000049454e44ae426082"
)


class TestAttachableExtensions:
    """Attachments were image-only, so the `.ab1` the vendor sent and the
    plate-reader CSV had nowhere to live — even though the app parses
    `.ab1` natively elsewhere."""

    def test_image_list_still_image_only(self):
        """`_IMAGE_EXTS` must not grow — it drives the inline preview,
        the clipboard path and the `![...]` markdown form."""
        assert set(sc._IMAGE_EXTS) == {
            ".png", ".jpg", ".jpeg", ".gif", ".bmp", ".webp", ".tiff", ".tif"}

    def test_attach_set_is_images_plus_data(self):
        assert set(sc._EXPERIMENT_ATTACH_EXTS) == (
            set(sc._IMAGE_EXTS) | set(sc._EXPERIMENT_DATA_EXTS))
        for want in (".ab1", ".csv", ".tsv", ".pdf", ".gb", ".fasta",
                     ".fastq", ".xlsx", ".txt"):
            assert want in sc._EXPERIMENT_ATTACH_EXTS, want

    def test_data_and_image_lists_are_disjoint(self):
        assert not (set(sc._IMAGE_EXTS) & set(sc._EXPERIMENT_DATA_EXTS))

    def test_is_image_path(self):
        assert sc._is_image_path("a.PNG") and sc._is_image_path("x/y.jpg")
        from pathlib import Path as _P
        assert sc._is_image_path(_P("b.tif"))
        for no in ("a.ab1", "a.pdf", "a.csv", "noext", "", None, 7):
            assert not sc._is_image_path(no), no

    def test_saves_a_data_file_under_the_att_prefix(self):
        out = sc._save_experiment_image("exp-att1234", b"ABIF\x00trace",
                                        "read1.ab1")
        assert out is not None
        assert out.suffix == ".ab1"
        assert out.name.startswith("att-"), out.name

    def test_images_keep_the_img_prefix(self):
        out = sc._save_experiment_image("exp-att1234", _TINY_PNG, "snap.png")
        assert out is not None
        assert out.name.startswith("img-") and out.suffix == ".png"

    def test_unknown_extension_still_falls_back_to_png(self):
        """Existing contract: an unrecognised suffix becomes .png (the
        clipboard-paste path has no filename at all)."""
        out = sc._save_experiment_image("exp-att1234", _TINY_PNG, "x.weird")
        assert out is not None
        assert out.suffix == ".png" and out.name.startswith("img-")

    def test_caps_still_apply_to_data_files(self):
        big = b"x" * (sc._EXPERIMENT_IMAGE_MAX_BYTES + 1)
        assert sc._save_experiment_image("exp-att1234", big, "big.csv") is None

    def test_picker_filter_offers_data_files(self):
        assert ".ab1" in sc._EXPERIMENT_ATTACH_FILE_FILTER
        assert ".ab1" not in sc._IMAGE_FILE_FILTER   # preview stays images

    def test_markdown_export_links_non_images(self):
        """An `![...]` around a `.ab1` renders as a broken-image box."""
        e = _entry("exp-1", "T", "body", image_paths=["img-1.png",
                                                      "att-1.ab1"])
        doc = sc._experiment_markdown_document(e)
        assert "![img-1.png](img-1.png)" in doc
        assert "[att-1.ab1](att-1.ab1)" in doc
        assert "![att-1.ab1]" not in doc

    def test_html_export_links_non_images(self):
        e = _entry("exp-1", "T", "body", image_paths=["img-1.png",
                                                      "att-1.ab1"])
        html = sc._experiment_html_document([e])
        assert '<img src="img-1.png"' in html
        assert 'href="att-1.ab1"' in html and 'download="att-1.ab1"' in html
        assert '<img src="att-1.ab1"' not in html

    def test_agent_endpoint_accepts_a_data_file(self, tmp_path):
        sc._save_experiments([_entry("exp-att1234", "T", "")])
        p = tmp_path / "read1.ab1"
        p.write_bytes(b"ABIF\x00trace")
        body = sc._state._AGENT_HANDLERS["attach-experiment-image"][0](
            None, {"experiment_id": "exp-att1234", "path": str(p)})
        assert body["ok"] is True
        assert body["is_image"] is False
        assert body["filename"].endswith(".ab1")
        # Body reference must be a link, not an image tag.
        entry = next(e for e in sc._load_experiments()
                     if e["id"] == "exp-att1234")
        assert f"[{body['filename']}]({body['filename']})" in entry["body_md"]
        assert f"![{body['filename']}]" not in entry["body_md"]

    def test_agent_endpoint_still_flags_an_image(self, tmp_path):
        sc._save_experiments([_entry("exp-att1234", "T", "")])
        p = tmp_path / "snap.png"
        p.write_bytes(_TINY_PNG)
        body = sc._state._AGENT_HANDLERS["attach-experiment-image"][0](
            None, {"experiment_id": "exp-att1234", "path": str(p)})
        assert body["is_image"] is True
        entry = next(e for e in sc._load_experiments()
                     if e["id"] == "exp-att1234")
        assert f"![{body['filename']}]" in entry["body_md"]

    def test_agent_endpoint_still_refuses_an_unlisted_extension(self, tmp_path):
        sc._save_experiments([_entry("exp-att1234", "T", "")])
        p = tmp_path / "payload.exe"
        p.write_bytes(b"MZ")
        out = sc._state._AGENT_HANDLERS["attach-experiment-image"][0](
            None, {"experiment_id": "exp-att1234", "path": str(p)})
        body, status = out
        assert status == 400 and "not an attachable file" in body["error"]


# ═════════════════════════════════════════════════════════════════════════════
# Hardening — malformed input, junk payloads, empty states
# ═════════════════════════════════════════════════════════════════════════════

class TestStringWhereListBelongs:
    """A hand-edited `experiments.json` can put a bare string where a
    list belongs, and `[x for x in value if isinstance(x, str)]` then
    iterates it as CHARACTERS. Every one of these produced confidently
    wrong output rather than an error."""

    def test_as_str_list_refuses_a_bare_string(self):
        assert sc._as_str_list("gibson") == []
        assert sc._as_str_list(["a", None, 7, "b"]) == ["a", "b"]
        assert sc._as_str_list(("a",)) == ["a"]
        for junk in (None, 7, {"a": 1}, {"a"}, 1.5):
            assert sc._as_str_list(junk) == [], junk

    def test_tag_filter_does_not_match_a_character(self):
        e = _entry("a", "t", "b")
        e["tags"] = "gibson"          # not a list
        assert sc._experiment_search([e], "", tags=["g"]) == []
        # ...but it IS the tag `gibson` — the comma-separated text the tags
        # field takes. Reading it as NO tags made the next save delete it
        # (hardening 2026-09-24, `_tag_values`).
        hits = sc._experiment_search([e], "", tags=["gibson"])
        assert [h["entry"]["id"] for h in hits] == ["a"]

    def test_protocol_inputs_are_not_split_into_characters(self):
        md = sc._protocol_steps_markdown(
            [{"op": "PCR", "product": "p", "inputs": "abc",
              "backbone": "", "enzymes": "BsaI", "where": ""}])
        assert "`a` + `b` + `c`" not in md
        assert "`B` + `s`" not in md

    def test_template_actions_are_not_split_into_characters(self):
        assert sc._experiment_template_markdown(
            "pcr", sc._EXPERIMENT_ACTIONS) == ""

    def test_export_attachments_are_not_split_into_characters(self):
        e = _entry("a", "t", "b")
        e["image_paths"] = "img.png"      # not a list
        doc = sc._experiment_markdown_document(e)
        assert "[i](i)" not in doc and "## Attachments" not in doc
        assert "[i](i)" not in sc._experiment_html_document([e])


class TestProtocolInputCap:
    def test_caps_inputs_like_the_history_viewer_does(self):
        step = {"op": "Golden Gate", "product": "p",
                "inputs": [f"part{i}" for i in range(200)],
                "backbone": "v", "enzymes": ["BsaI"], "where": ""}
        md = sc._protocol_steps_markdown([step])
        assert f"+{200 - sc._PROTOCOL_INPUT_MAX} more" in md
        assert "part9" not in md          # past the cap
        assert len(md) < 600              # was one 2,500-char line

    def test_cap_matches_the_history_viewer_constant(self):
        """Drift guard — the notebook protocol mirrors the History pane's
        cap without importing history L2."""
        assert sc._PROTOCOL_INPUT_MAX == sc._HISTORY_PROTOCOL_INPUT_MAX

    def test_no_more_marker_under_the_cap(self):
        md = sc._protocol_steps_markdown(
            [{"op": "", "product": "p", "inputs": ["a", "b"],
              "backbone": "", "enzymes": [], "where": ""}])
        assert "more" not in md


class TestTemplateTagRoundTrips:
    def test_an_unusable_id_gets_a_heading_but_no_tag(self):
        """`_ACTIONS_REF_RE` caps an id at 64 chars, so a longer one would
        render a tag that can never be highlighted or extracted."""
        long_id = "z" * 5000
        md = sc._experiment_template_markdown([long_id], ())
        assert f"## {long_id}" in md
        assert "!" not in md
        assert sc._extract_action_refs(md) == []

    def test_a_usable_id_keeps_its_tag_and_round_trips(self):
        md = sc._experiment_template_markdown(["my-step"], ())
        assert "!my-step" in md
        assert sc._extract_action_refs(md) == ["my-step"]


class TestBacklinkRefTypes:
    def test_none_and_empty_are_empty(self, entries):
        assert sc._experiments_referencing(entries, None) == []
        assert sc._experiments_referencing(entries, "") == []
        assert sc._experiments_referencing(entries, []) == []

    def test_wrong_type_raises_rather_than_answering_none(self, entries):
        """A dict is iterable (yielding KEYS) and a generator is
        consumable once, so coercing either would answer "nothing
        references that" for a query that was simply malformed."""
        for junk in (7, 1.5, {"a": 1}, (r for r in ["x"])):
            with pytest.raises(TypeError, match="ref_id must be"):
                sc._experiments_referencing(entries, junk)


class TestSnippetEmptyTerm:
    def test_an_empty_term_does_not_snippet_the_head(self):
        """`"".find` returns 0, which would show the start of the body
        for a query that matched nothing."""
        assert sc._experiment_snippet("hello world", [""]) == ""
        assert sc._experiment_snippet("hello world", ["", "world"]) != ""


class TestProjectDocumentCount:
    def test_counts_what_it_rendered_not_the_argument(self):
        """A non-list argument produced a header claiming "3 entries"
        above a document containing none."""
        doc = sc._experiment_project_markdown("abc", project="P")
        assert "0 entries" in doc
        doc2 = sc._experiment_project_markdown(
            ["junk", None, _entry("a", "T")], project="P")
        assert "1 entry" in doc2


class TestExportUrlSafety:
    @pytest.mark.parametrize("url", [
        "//evil.com/pixel.png",        # protocol-relative → https://evil.com
        "\\\\host\\share\\x.png",      # Windows UNC
    ])
    def test_refuses_scheme_inheriting_urls(self, url):
        """An export is an HTML file opened in a browser; a pasted
        protocol-relative link would fire an outbound request. Explicit
        http(s) is already allowed, so refusing these costs nothing."""
        assert sc._html_safe_url(url) == ""

    def test_still_allows_a_plain_relative_path(self):
        assert sc._html_safe_url("/abs/path.png") == "/abs/path.png"
        assert sc._html_safe_url("sub/img.png") == "sub/img.png"

    def test_html_export_drops_them(self):
        e = _entry("a", "T", "b", image_paths=["x.png"])
        html = sc._experiment_html_document(
            [e], image_srcs={"x.png": "//evil.com/p.png"})
        assert "evil.com" not in html


class TestAgentPayloadFuzz:
    """Junk payloads must produce a 4xx with an `error` key, never an
    unhandled exception (which the HTTP layer turns into a 500). This
    found `{"format": 1}` raising AttributeError in export-experiment."""

    _ENDPOINTS = ["search-experiments", "experiment-backlinks",
                  "get-plasmid-protocol", "experiment-step-template",
                  "list-experiment-actions", "export-experiment",
                  "duplicate-experiment", "attach-experiment-image"]

    _JUNK = [None, "", 0, -1, 1.5, True, False, "x" * 5000, [], {}, [None],
             {"a": 1}, ["ok", None, 7], float("inf"), float("nan"),
             "../../etc/passwd", "\x00bad", "@#$%^&*", [[]], {"k": []}]

    _KEYS = ["id", "name", "ref", "kind", "query", "tags", "project",
             "all_projects", "limit", "format", "actions", "heading",
             "title", "collection", "experiment_id", "path"]

    def test_no_unhandled_exception_and_shape_holds(self):
        sc._save_experiments([_entry("exp-aa", "T", "@pUC19", ["a"])])
        bad = []
        for name in self._ENDPOINTS:
            fn = sc._state._AGENT_HANDLERS[name][0]
            payloads = [{}] + [{k: v} for k in self._KEYS
                               for v in self._JUNK]
            for pl in payloads:
                try:
                    out = fn(None, pl)
                except Exception as exc:          # noqa: BLE001
                    bad.append((name, pl, type(exc).__name__, str(exc)[:90]))
                    continue
                if isinstance(out, tuple):
                    body, status = out
                    if not isinstance(body, dict):
                        bad.append((name, pl, "BAD_BODY", repr(out)[:90]))
                    elif not isinstance(status, int) or not (
                            400 <= status < 600):
                        bad.append((name, pl, "BAD_STATUS", status))
                    elif "error" not in body:
                        bad.append((name, pl, "NO_ERROR_KEY",
                                    repr(body)[:90]))
                elif not isinstance(out, dict):
                    bad.append((name, pl, "BAD_SHAPE", repr(out)[:90]))
        assert not bad, f"{len(bad)} bad responses: {bad[:10]}"


class TestNotebookUIEmptyStates:
    async def test_copy_with_nothing_selected_warns(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([])
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            assert scr.query_one("#exp-entries-table").row_count == 0
            scr._duplicate_current()          # must not raise
            await pilot.pause()
            assert sc._load_experiments() == []

    async def test_export_with_an_empty_notebook_warns(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([])
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.action_export_entry()
            await pilot.pause()
            # No modal — nothing to export.
            assert not isinstance(app.screen, sc.ExperimentExportModal)

    async def test_protocol_and_steps_need_an_entry(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([])
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            assert scr._current_entry is None
            scr.action_insert_protocol()
            await pilot.pause()
            assert not isinstance(app.screen, sc.LibrarySearchModal)
            scr.action_insert_steps()
            await pilot.pause()
            assert not isinstance(app.screen, sc.ActionsPickerModal)

    async def test_jump_to_a_vanished_entry_clears_instead_of_crashing(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([_entry("exp-aa", "T", "b")])
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr._jump_to_entry("", "exp-gone")
            await pilot.pause()
            assert scr._current_entry is None

    async def test_jump_to_a_vanished_project_reports_it(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([_entry("exp-aa", "T", "b")])
            app.push_screen(sc.ExperimentsScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            before = [e["id"] for e in sc._load_experiments()]
            scr._jump_to_entry("No Such Project", "exp-aa")
            await pilot.pause()
            # Refused the switch; the live mirror is untouched.
            assert [e["id"] for e in sc._load_experiments()] == before

    async def test_search_modal_with_an_empty_notebook(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([])
            app.push_screen(sc.ExperimentSearchModal())
            await pilot.pause()
            await pilot.pause()
            m = app.screen
            assert m._matches == []
            m._on_query_submitted(None)       # must not raise
            await pilot.pause()
            assert isinstance(app.screen, sc.ExperimentSearchModal)

    async def test_search_modal_with_an_unknown_backlink_kind(self):
        """`_experiments_referencing` raises on a bad kind; the modal
        shows nothing rather than everything."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([_entry("exp-aa", "T", "@pUC19")])
            app.push_screen(sc.ExperimentSearchModal(
                backlink=("bogus", "pUC19")))
            await pilot.pause()
            await pilot.pause()
            assert app.screen._matches == []

    async def test_export_modal_scope_omits_entry_when_there_is_none(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            sc._save_experiments([_entry("exp-aa", "T", "b")])
            m = sc.ExperimentExportModal(entry=None, project="Main")
            app.push_screen(m)
            await pilot.pause()
            await pilot.pause()
            assert m.query_one("#expexp-scope").value == "project"
            assert m._entries_for_scope()[0]["id"] == "exp-aa"

    async def test_export_modal_reads_the_entry_back_from_disk(self):
        """A save after the modal opened must not export stale text."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            entry = _entry("exp-aa", "T", "old body")
            sc._save_experiments([entry])
            m = sc.ExperimentExportModal(entry=entry, project="Main")
            app.push_screen(m)
            await pilot.pause()
            await pilot.pause()
            fresh = dict(entry)
            fresh["body_md"] = "new body"
            sc._save_experiments([fresh])
            assert m._entries_for_scope()[0]["body_md"] == "new body"


class TestExportWriteIsAtomic:
    async def test_no_partial_file_on_a_failed_render(self, tmp_path,
                                                       monkeypatch):
        """The export worker is `exclusive=True` and gets cancelled when
        the modal closes, so a plain `write_text` could leave a truncated
        document that still looks finished."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            entry = _entry("exp-aa", "T", "body")
            sc._save_experiments([entry])
            m = sc.ExperimentExportModal(entry=entry, project="Main")
            app.push_screen(m)
            await pilot.pause()
            await pilot.pause()

            def _boom(*_a, **_k):
                raise OSError("disk full")
            monkeypatch.setattr(sc, "_atomic_write_bytes", _boom)
            # A subdirectory, because `tmp_path` itself is the sandboxed
            # data dir the _protect_user_data fixture points at.
            dest = tmp_path / "exports"
            dest.mkdir()
            m._selected_dir = str(dest)
            m.query_one("#expexp-filename").value = "out.md"
            m._do_export()
            for _ in range(60):
                await pilot.pause()
            # Nothing written, modal still up with the reason.
            assert not (dest / "out.md").exists()
            assert list(dest.iterdir()) == []          # no stray .tmp
            assert "disk full" in str(m.query_one("#expexp-status").render())

    async def test_unreadable_attachment_clears_the_embedded_flag(
            self, tmp_path, monkeypatch):
        """One unreadable image means the document is NOT self-contained;
        reporting a clean embed would hide a broken image."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            saved = sc._save_experiment_image("exp-aa", _TINY_PNG, "a.png")
            assert saved is not None
            entry = _entry("exp-aa", "T", "body",
                           image_paths=[saved.name])
            sc._save_experiments([entry])
            m = sc.ExperimentExportModal(entry=entry, project="Main")
            app.push_screen(m)
            await pilot.pause()
            await pilot.pause()
            srcs, embedded = m._image_srcs([entry], embed=True)
            assert embedded is True
            assert srcs[saved.name].startswith("data:image/")

            def _no_read(self, *_a, **_k):
                raise OSError("gone")
            monkeypatch.setattr(sc.Path, "read_bytes", _no_read)
            srcs2, embedded2 = m._image_srcs([entry], embed=True)
            assert embedded2 is False
            # local fallback — as a file: URI, which survives a Windows drive
            # letter and a space in the path (audit 2026-09-22)
            assert srcs2[saved.name] == saved.absolute().as_uri()
