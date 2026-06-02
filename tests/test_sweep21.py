"""
test_sweep21 — modal-staleness audit (2026-05-21).

Staleness window: a modal opens at T=0 capturing some piece of canvas
state (feature idx, selection range, record reference). Between T=0
and the modal's dismiss at T=2, an agent-API endpoint (or background
worker) can mutate `_current_record` via `_h_load_entry` / `_h_fetch`
/ `_h_load_file`. The modal's dismiss callback then applies user
edits to the WRONG molecule at the WRONG idx.

PlasmidApp already had `_guard_callback(callback, label)` that wraps
a callback to refuse if `_record_load_counter` shifted while the
modal was up. Pre-sweep, it was used by exactly ONE path
(`action_edit_seq`). Sweep #21 extends coverage to the other
captured-state callbacks:

  * `_open_feature_editor` — FeatureEditModal / PrimerEditModal
    capture `idx` at T=0; agent-swap of `_current_record` would land
    the user's edit on the new record's feature[idx], which is a
    different feature entirely.

  * `action_add_feature` — AddFeatureModal "annotate" path captures
    selection `(start, end)` at T=0; same risk.

Both now wrap their callback via `_guard_callback`.

`action_transfer_annotations` already had an inline counter check
(lines 73066-73082). `action_edit_seq` was already guarded.

White-box test: monkey-patch the dismiss callback path to verify
the guard short-circuits when the counter has shifted.
"""
from __future__ import annotations

from unittest.mock import MagicMock

import splicecraft as sc


class TestGuardCallbackBehavior:
    """`_guard_callback` itself is the reusable primitive; verify its
    contract first so callsites can rely on it.
    """

    def _make_app(self):
        # Minimal stub — only `_record_load_counter` and `notify` are
        # exercised by `_guard_callback`.
        app = MagicMock()
        app._record_load_counter = 0
        app.notify = MagicMock()
        # Bind the real method to this stub via method binding so
        # `self.X` semantics work.
        app._guard_callback = sc.PlasmidApp._guard_callback.__get__(
            app, sc.PlasmidApp,
        )
        return app

    def test_callback_fires_when_counter_unchanged(self):
        app = self._make_app()
        inner = MagicMock()
        wrapped = app._guard_callback(inner, "Test")
        wrapped("result")
        inner.assert_called_once_with("result")
        app.notify.assert_not_called()

    def test_callback_dropped_when_counter_shifted(self):
        app = self._make_app()
        inner = MagicMock()
        wrapped = app._guard_callback(inner, "Test")
        # Simulate canvas swap between modal open and dismiss.
        app._record_load_counter = 5
        wrapped("result")
        inner.assert_not_called()
        app.notify.assert_called_once()
        # The notification should mention the label so the user
        # understands what was dropped.
        first_arg = app.notify.call_args[0][0]
        assert "Test" in first_arg
        assert "dropped" in first_arg.lower() \
            or "changed" in first_arg.lower()

    def test_label_default(self):
        app = self._make_app()
        wrapped = app._guard_callback(lambda _: None)
        app._record_load_counter = 1
        wrapped("payload")
        # Default label is "Action" — used when callsite doesn't
        # specify (currently no callsites omit, but the API allows).
        first_arg = app.notify.call_args[0][0]
        assert "Action" in first_arg


class TestGuardCallbackCallsites:
    """Verify the new sweep #21 callsites actually wrap their
    callback through `_guard_callback`. White-box check by scanning
    the source for the guarded-dispatch line.
    """

    def test_open_feature_editor_uses_guard(self):
        # Grep the method body via inspect.
        import inspect
        src = inspect.getsource(sc.PlasmidApp._open_feature_editor)
        assert "_guard_callback" in src, (
            "_open_feature_editor should wrap its dismiss callback "
            "via _guard_callback so a canvas reload mid-edit drops "
            "the edit instead of corrupting the new record."
        )

    def test_action_add_feature_uses_guard(self):
        import inspect
        src = inspect.getsource(sc.PlasmidApp.action_add_feature)
        assert "_guard_callback" in src, (
            "action_add_feature should wrap its callback via "
            "_guard_callback so an agent-API record swap between "
            "modal open and Annotate-button click doesn't land the "
            "(start, end) coords on the wrong molecule."
        )

    def test_action_edit_seq_uses_guard(self):
        # Existing usage — regression-lock so a future refactor
        # doesn't accidentally remove it. The guard lives in the shared
        # opener `_open_seq_edit_dialog` since 2026-06-01 (Ctrl+E and the
        # Delete key both route through it); `action_edit_seq` delegates.
        import inspect
        assert "_open_seq_edit_dialog" in inspect.getsource(
            sc.PlasmidApp.action_edit_seq)
        assert "_guard_callback" in inspect.getsource(
            sc.PlasmidApp._open_seq_edit_dialog)

    def test_action_transfer_annotations_uses_inline_counter(self):
        # Different shape — inline `_record_load_counter` check
        # rather than `_guard_callback` because the action splits
        # across two modal callbacks (source picker + apply modal)
        # with state captured in between. Both layers need the same
        # counter; `_guard_callback` would only cover one.
        import inspect
        src = inspect.getsource(
            sc.PlasmidApp.action_transfer_annotations,
        )
        assert "_record_load_counter" in src
        assert "entry_counter" in src
