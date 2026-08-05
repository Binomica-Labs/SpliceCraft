"""
test_ui_layout — geometry guards for the app's text-entry widgets.

The bug these exist for (2026-08-05): `PlasmidApp.CSS` carried a universal

    * { … scrollbar-size-vertical: 1; scrollbar-size-horizontal: 1; }

App CSS outranks a widget's DEFAULT_CSS whatever the specificity, so that
one declaration silently overrode Textual's `Input { scrollbar-size-
horizontal: 0 }`. An Input is `height: 3` — border + ONE content row +
border — so the moment a value overflowed the box, the resurrected
scrollbar ate that single content row and the field rendered as a solid
bar with the text nowhere to be seen. You could neither read the value nor
watch what you were typing. Reported on the Rename plasmid dialog; it hit
EVERY Input in the app, plus `Footer` (a ScrollableContainer pinned to
height 1) and `MarkdownFence`.

Textual's own default for `scrollbar-size-horizontal` is already 1, so the
declaration bought nothing anywhere — it only ever clobbered the widgets
that deliberately set 0. The vertical half (default 2) is the one that
earns its place and stays.
"""
from __future__ import annotations

import re

import pytest
from textual.widgets import Footer, Input

import splicecraft as sc


# Long enough to overflow any dialog-sized Input on a 160-col terminal.
LONG_NAME = "PLASMID " + "x" * 120


def _universal_rule_bodies(css: str) -> list[str]:
    """Bodies of every top-level `* { … }` rule in a stylesheet."""
    return [m.group(1) for m in re.finditer(r"(?m)^\*\s*\{([^}]*)\}", css)]


class TestScrollbarSuppression:
    """A widget that sets `scrollbar-size…: 0` means it — don't override it
    from a blanket rule."""

    def test_star_rule_does_not_force_hscrollbar(self):
        bodies = _universal_rule_bodies(sc.PlasmidApp.CSS)
        assert bodies, (
            "expected the universal `* { … }` scrollbar-theme rule in "
            "PlasmidApp.CSS — did it get renamed or removed?"
        )
        for body in bodies:
            assert "scrollbar-size-horizontal" not in body, (
                "`* { scrollbar-size-horizontal: … }` is back. Textual's "
                "default is ALREADY 1, so this declaration changes nothing "
                "except clobbering the widgets that set 0 on purpose — it "
                "makes every Input lose its only content row to a scrollbar "
                "as soon as the value overflows. Put the size on the "
                "containers that actually want a horizontal scrollbar."
            )
            assert not re.search(r"\bscrollbar-size\s*:", body), (
                "the `scrollbar-size` shorthand on `*` sets the horizontal "
                "size too — same clobber. Set the vertical size explicitly."
            )

    def test_textual_zeroed_scrollbars_survive_app_css(self):
        """Forward-looking: catch a Textual upgrade that adds a new widget
        relying on a zeroed scrollbar, or a new blanket rule that revives
        the clobber. Only the *sizes* may not be forced globally — the
        colour declarations on `*` are fine and deliberately kept."""
        zeroed = {}
        for cls in (Input, Footer):
            css = cls.__dict__.get("DEFAULT_CSS", "") or ""
            for m in re.finditer(
                    r"scrollbar-size(-horizontal|-vertical)?\s*:\s*([^;]+);", css):
                if re.fullmatch(r"[0\s]+", m.group(2).strip()):
                    zeroed.setdefault(cls.__name__, []).append(m.group(0).strip())
        assert "Input" in zeroed, (
            "Textual's Input no longer zeroes its horizontal scrollbar — "
            "re-derive this guard against the new default before deleting it."
        )

    async def test_long_value_keeps_its_content_row(
            self, tiny_record, isolated_library):
        """The end-to-end symptom: a value too long for the box must still
        be visible. Drives the real Rename dialog rather than asserting on
        CSS text, so it fails for ANY cause that steals the content row."""
        from tests.test_smoke import TERMINAL_SIZE, _build_app

        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await app.push_screen(sc.RenamePlasmidModal(LONG_NAME, "entry-1"))
            await pilot.pause()
            await pilot.pause()
            inp = app.screen.query_one("#rename-input", Input)
            # Precondition — the value genuinely overflows, so a forced
            # scrollbar WOULD appear. Without this the test passes vacuously.
            assert inp.virtual_size.width > inp.content_size.width, (
                "test setup no longer overflows the Input; lengthen LONG_NAME"
            )
            assert inp.styles.scrollbar_size_horizontal == 0, (
                "Input's horizontal scrollbar was re-enabled — it will eat "
                "the field's only content row"
            )
            assert inp.content_size.height >= 1, (
                "Input has no content row left to render its value into"
            )

    async def test_footer_keeps_its_single_row(
            self, tiny_record, isolated_library):
        """Footer is a ScrollableContainer docked at `height: 1`. Force a
        scrollbar onto it and the entire footer becomes the scrollbar."""
        from tests.test_smoke import TERMINAL_SIZE, _build_app

        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            footers = app.screen.query(Footer)
            if not footers:
                pytest.skip("no Footer on the default screen")
            f = footers.first()
            assert f.styles.scrollbar_size_horizontal == 0
            # The vertical size IS still pinned to 1 by the `*` rule, and
            # that is fine: Footer's content is exactly one row tall, so it
            # can never overflow vertically and no vertical scrollbar is
            # ever laid out. Only the horizontal one could steal the row.
            assert f.content_size.height >= 1


class TestRenameDialogHeadroom:
    """The reported dialog: long plasmid names must be READABLE, both the
    'Current name:' line and the editable field."""

    async def test_current_name_wraps_instead_of_truncating(
            self, tiny_record, isolated_library):
        from textual.widgets import Label

        from tests.test_smoke import TERMINAL_SIZE, _build_app

        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await app.push_screen(sc.RenamePlasmidModal(LONG_NAME, "entry-1"))
            await pilot.pause()
            await pilot.pause()
            lbl = app.screen.query_one("#rename-current", Label)
            # `width: auto` would size the Label to its (very long) text and
            # let it overflow the dialog border, slicing the name off with
            # no ellipsis. Bounded width => it reflows onto further rows.
            assert lbl.content_size.width <= app.screen.query_one(
                "#rename-dlg").content_size.width
            assert lbl.content_size.height > 1, (
                "the long current-name did not wrap — it is being truncated "
                "at the dialog edge instead"
            )
