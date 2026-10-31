"""Attack-surface / input-sterilization sweep (2026-06-14).

Locks in the hardening from the fresh attack-surface review:

  * control-char / ANSI / C1 / lone-surrogate / line-separator stripping on
    user-supplied names BEFORE they are persisted and reflected to the
    terminal (`_sanitize_label`, `_sanitize_note`, the logger scrub filter,
    `_record_display_name`, the custom-enzyme + experiment-project agent
    handlers, the primer-CSV import);
  * the spellcheck markdown-link/image ReDoS fix (now O(n), not O(n^2));
  * zip member-name tightening (control bytes + lone surrogates rejected);
  * the foreign-file ingest cap (`_GB_INGEST_MAX_BYTES`, 256 MB);
  * the empty-GenBank-text guard;
  * the agent server's DNS-rebinding Host check.

All pure / handler-level (fast) bar one E2E for the Host header. The autouse
`_protect_user_data` fixture (conftest) sandboxes every data-file write.
"""
from __future__ import annotations

import socket
import threading
import time
import types
import urllib.error
import urllib.request

import pytest

import splicecraft as sc

ESC = "\x1b"        # C0 escape — terminal-injection prefix
CSI8 = "\x9b"       # 8-bit CSI (C1 block)
SURR = "\udce9"     # lone surrogate (surrogateescape of a non-UTF-8 byte)
LSEP = "\u2028"    # Unicode line separator


# ── Control-char sterilization helpers ──────────────────────────────────────

def test_sanitize_label_strips_control_c1_surrogate_separator():
    clean = sc._sanitize_label(f"p{ESC}[2J{CSI8}{SURR}{LSEP}UC19")
    assert not any(ch in clean for ch in (ESC, CSI8, SURR, LSEP))
    assert clean.startswith("p") and clean.endswith("UC19")


def test_sanitize_label_preserves_spaces_and_plus():
    # [INV-98] no-underscores rule: internal spaces + '+' survive untouched.
    assert sc._sanitize_label("FFE 6") == "FFE 6"
    assert sc._sanitize_label("pED + insert") == "pED + insert"


def test_sanitize_label_control_only_becomes_empty():
    assert sc._sanitize_label(SURR + ESC + CSI8) == ""


def test_sanitize_note_keeps_tab_newline_strips_esc():
    out = sc._sanitize_note(f"para 1\n\tindented{ESC}[31m\npara 2")
    assert "\n" in out and "\t" in out
    assert ESC not in out and CSI8 not in out


def test_log_scrub_escapes_control_bytes():
    # A terminal-escape smuggled via any logged value is rendered as a visible
    # \xNN token, never a live ESC byte that fires when the log is cat'd.
    scrubbed = sc._SurrogateScrubFilter._scrub(f"name{ESC}[2Jx")
    assert ESC not in scrubbed
    assert "\\x1b" in scrubbed


def test_log_scrub_makes_surrogate_utf8_safe():
    scrubbed = sc._SurrogateScrubFilter._scrub("a" + SURR + "b")
    scrubbed.encode("utf-8")  # must not raise (the raw value would)


def test_record_display_name_sanitizes_terminal_escape():
    rec = types.SimpleNamespace(name="locus", id="id1")
    rec._tui_display_name = f"My{ESC}[2J Plasmid"
    name = sc.PlasmidApp._record_display_name(rec)
    assert ESC not in name
    assert " " in name  # display name keeps spaces


def test_record_display_name_falls_back_through_name():
    rec = types.SimpleNamespace(name=f"loc{ESC}us", id="id1")
    name = sc.PlasmidApp._record_display_name(rec)
    assert ESC not in name and name


# ── Agent custom-enzyme + experiment-project name sterilization ─────────────

def test_custom_enzyme_payload_sanitizes_text_fields():
    out = sc._agent_validate_custom_enzyme_payload({
        "name":     f"Eco{ESC}[31mRX",
        "site":     "GAATTC",
        "fwd_cut":  1, "rev_cut": 5,
        "type":     f"II{CSI8}5",
        "supplier": f"acme{ESC}co",
    })
    assert isinstance(out, dict), out
    for field in ("name", "type", "supplier"):
        assert ESC not in out[field] and CSI8 not in out[field]
    assert out["site"] == "GAATTC"  # unchanged


def test_experiment_project_create_sanitizes_name():
    res = sc._h_create_experiment_project(None, {"name": f"Proj{ESC}[2JX"})
    assert isinstance(res, dict) and res.get("ok") is True, res
    assert ESC not in res["name"]


# ── Primer CSV import name sterilization ────────────────────────────────────

def test_primer_csv_import_sanitizes_name(tmp_path):
    csv_path = tmp_path / "primers.csv"
    csv_path.write_text(
        "Name,Sequence\n" f"fwd{ESC}[2J,ACGTACGTACGTACGTACGT\n",
        encoding="utf-8",
    )
    res = sc._import_primers_from_csv(str(csv_path))
    assert res["primers"], res["skipped"]
    assert ESC not in res["primers"][0]["name"]


# ── Spellcheck markdown link/image ReDoS ────────────────────────────────────

def test_spellcheck_strip_is_linear_on_bracket_run():
    bomb = "[" * 60000  # O(n^2) form took seconds here; linear form ~1 ms
    t0 = time.process_time()
    out = sc._spellcheck_strip_code(bomb)
    spent = time.process_time() - t0
    # CPU time: the claim is "linear", and a wall-clock budget inside
    # `-n auto` fails on scheduler starvation rather than on a regression.
    assert spent < 1.0, f"strip took {spent:.3f} s CPU (budget 1.0)"
    assert len(out) == len(bomb)  # offset alignment preserved


def test_spellcheck_strip_still_masks_links_and_images():
    masked = sc._spellcheck_strip_code("see [text](http://x/y) and ![a](u) ok")
    assert "see " in masked and " ok" in masked
    assert "text" not in masked and "http" not in masked


# ── Zip member-name tightening ──────────────────────────────────────────────

@pytest.mark.parametrize("bad", [
    "sample\x07.gbk",    # bell (C0)
    "sample\rfoo.gbk",   # carriage return
    "sample\x9bfoo.gbk",  # C1 8-bit CSI
    "sample\udce9.gbk",  # lone surrogate
])
def test_is_safe_zip_member_name_rejects_control_and_surrogate(bad):
    assert sc._is_safe_zip_member_name(bad) is False


def test_is_safe_zip_member_name_accepts_normal():
    assert sc._is_safe_zip_member_name("run/34XK5N_sample.gbk") is True


# ── Foreign-file ingest cap + empty-GenBank guard ───────────────────────────

def test_gb_ingest_cap_is_256mb_and_below_json_cap():
    assert sc._GB_INGEST_MAX_BYTES == 256 * 1024 * 1024
    # MUST stay below the data-dir JSON cap so the two never get conflated.
    assert sc._GB_INGEST_MAX_BYTES < sc._state._SAFE_LOAD_JSON_MAX_BYTES


def test_gb_text_to_record_empty_raises_clean():
    with pytest.raises(ValueError):
        sc._gb_text_to_record("")


# ── Agent server DNS-rebinding Host check (E2E) ─────────────────────────────

class _TinyApp:
    _current_record = None
    _unsaved = False


def _free_port() -> int:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


@pytest.fixture
def host_server():
    port = _free_port()
    srv = sc._AgentAPIServer(("127.0.0.1", port), _TinyApp(), f"tok-{port}")
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    time.sleep(0.05)
    try:
        yield port
    finally:
        srv.shutdown()
        srv.server_close()


def _get_status(port, path, host_header=None) -> int:
    req = urllib.request.Request(
        f"http://127.0.0.1:{port}{path}", method="GET")
    if host_header is not None:
        req.add_header("Host", host_header)
    try:
        with urllib.request.urlopen(req, timeout=5) as resp:
            return resp.status
    except urllib.error.HTTPError as e:
        return e.code


def test_agent_rejects_nonloopback_host(host_server):
    # DNS-rebinding: connect to the loopback IP, send an attacker Host header.
    assert _get_status(host_server, "/tools", "evil.example.com") == 403
    assert _get_status(host_server, "/", "attacker.test") == 403


def test_agent_allows_loopback_host(host_server):
    assert _get_status(host_server, "/tools") == 200  # default 127.0.0.1 Host
    assert _get_status(host_server, "/tools", f"localhost:{host_server}") == 200


# ── audit 2026-10-23: escape sequences in a FILE's feature text ─────────────

_ESC = "\x1b"
_OSC52 = _ESC + "]52;c;UFdORUQ=" + _ESC + "\\"      # OSC 52 clipboard write
_OSC0 = _ESC + "]0;PWNED" + _ESC + "\\"              # OSC 0 window title


def _hostile_gb() -> str:
    seq = "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCG" * 8
    lines = [f"LOCUS       ESCTEST               {len(seq)} bp    DNA     "
             f"circular SYN 01-JAN-2026",
             "FEATURES             Location/Qualifiers",
             "     CDS             1..300",
             f'                     /label="lab{_OSC52}el"',
             f'                     /note="note{_OSC52}{_OSC0}"',
             "     misc_feature    310..400",
             f'                     /SpliceCraft_type="t{_OSC52}"',
             "ORIGIN"]
    for i in range(0, len(seq), 60):
        c = seq[i:i + 60]
        lines.append(f"{i + 1:>9} " + " ".join(c[j:j + 10]
                                              for j in range(0, len(c), 10)))
    return "\n".join(lines + ["//", ""])


def test_feature_labels_and_types_carry_no_terminal_escapes():
    """A `/label` is painted straight to the terminal (map, sidebar, sequence
    lanes): its ESC was a live `ESC ] 52` clipboard write on OSC-52
    terminals. Restored long types are file text too."""
    rec = sc._gb_text_to_record(_hostile_gb(), cache=False)
    for f in rec.features:
        assert _ESC not in sc._feat_label_full(f)
        assert _ESC not in f.type


async def test_hostile_labels_never_reach_the_terminal_stream():
    rec = sc._gb_text_to_record(_hostile_gb(), cache=False)
    app = sc.PlasmidApp()
    async with app.run_test(size=(200, 60)) as pilot:
        for _ in range(4):
            await pilot.pause()
        app._apply_record(rec)
        for _ in range(6):
            await pilot.pause()
        out = app.screen._compositor.render_update(full=True) \
            .render_segments(app.console)
    assert (_ESC + "]") not in out


def _hostile_migrate_archive(tmp_path, settings_entries, plugin_files):
    import hashlib
    import json
    import zipfile
    settings_bytes = json.dumps(
        {"_schema_version": 1, "entries": settings_entries}).encode()
    manifest = {"schema_version": 1,
                "files": [{"attr": "_SETTINGS_FILE", "name": "settings.json",
                           "sha256": hashlib.sha256(settings_bytes)
                           .hexdigest()}],
                "directories": [{"attr": "_PLUGINS_DIR", "name": "plugins"}]}
    zpath = tmp_path / "colleague-library.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("splicecraft-migrate.json", json.dumps(
            {"format": "splicecraft-migrate", "format_version": 1,
             "snapshot_subdir": "data"}))
        zf.writestr("data/manifest.json", json.dumps(manifest))
        zf.writestr("data/settings.json", settings_bytes)
        for name in plugin_files:
            zf.writestr(f"data/plugins/{name}", "")
    return zpath


def test_a_migrate_archive_cannot_choose_where_babs_runs_code(
        tmp_path, monkeypatch):
    """An archive is a file someone can send you. Its settings.json pointed
    `babs_home` at its own `plugins/` (holding `bb_learn.py` + `rag_bot.py`),
    and every BABS chat turn then ran the archive's script (audit
    2026-10-23)."""
    for k in ("SPLICECRAFT_BABS_HOME", "BABS_HOME"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    sc._set_setting("allow_online_search", False)
    sc._settings_flush_sync()
    plugins = sc._state._PLUGINS_DIR
    zpath = _hostile_migrate_archive(
        tmp_path,
        [{"key": "babs_home", "value": str(plugins)},
         {"key": "allow_online_search", "value": True},
         {"key": "show_restr", "value": True}],
        ["bb_learn.py", "rag_bot.py", "bb_snowball.py"])
    summary = sc._import_migrate_archive(zpath)
    sc._bust_all_data_caches()
    assert not (plugins / "rag_bot.py").exists(), "archive plugins installed"
    assert sc._learn_resolve_babs_home() is None
    assert sc._get_setting("babs_home", "") in ("", None)
    assert sc._get_setting("allow_online_search", None) is False
    # Ordinary settings still travel with the archive.
    assert sc._get_setting("show_restr", None) is True
    assert "babs_home" in summary["kept_local_settings"]


def test_babs_home_inside_the_data_dir_is_refused(tmp_path):
    fake = sc._state._DATA_DIR / "somewhere"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "bb_learn.py").write_text("")
    assert sc._babs_home_allowed(fake) is False
    elsewhere = tmp_path.parent / (tmp_path.name + "-babs")
    elsewhere.mkdir(exist_ok=True)
    assert sc._babs_home_allowed(elsewhere) is True


def test_version_probe_never_imports_from_the_launch_directory(
        tmp_path, monkeypatch):
    """`python -c` puts the CWD first on sys.path; the metadata read imports
    stdlib modules, so a planted `operator.py` in the launch directory ran on
    `splicecraft update` and every authenticated agent response."""
    planted = tmp_path / "launch-dir"
    planted.mkdir()
    marker = tmp_path / "ran"
    for mod in ("operator", "json", "email", "zipfile", "pathlib",
                "contextlib", "typing", "warnings", "csv", "textwrap"):
        (planted / f"{mod}.py").write_text(
            f"open({str(marker)!r}, 'a').write({mod!r})\n"
            f"raise ImportError('planted')\n")
    monkeypatch.chdir(planted)
    sc._query_installed_version(timeout=30)
    assert not marker.exists(), marker.read_text()


def test_cds_gz_is_never_fully_inflated_before_the_size_check(
        tmp_path, monkeypatch):
    """The size check ran AFTER `gzip.decompress` inflated the whole file —
    a 1.8 MB `.gz` became 1.75 GiB in RAM first (OOM from the GUI picker or
    the agent's `add-codon-table`)."""
    import gzip
    import splicecraft_codon as scc
    monkeypatch.setattr(scc, "_CDS_FILE_MAX_BYTES", 1_000_000)
    bomb = tmp_path / "cds_from_genomic.fna.gz"
    bomb.write_bytes(gzip.compress(b">x\n" + b"A" * 5_000_000))
    whole = []
    real = gzip.decompress
    monkeypatch.setattr(gzip, "decompress",
                        lambda *a, **k: whole.append(1) or real(*a, **k))
    raw, msg, meta = scc._file_build_codon_table(bomb)
    assert raw is None and "too large" in msg
    assert not whole, "the whole bomb was inflated before the cap applied"


def test_notebook_attachments_are_confined_to_the_entry_folder(tmp_path):
    """An absolute / `..` image path in experiments.json (hand-edited, or
    from an imported archive) made the HTML export read and embed an
    arbitrary local file into a document the user then shares."""
    import base64
    secret = tmp_path / "id_rsa_DECOY"
    secret.write_text("DECOY-PRIVATE-KEY-MATERIAL\n")
    entry = sc._normalise_experiment_entry(
        {"id": "exp1", "title": "t", "body_md": "x",
         "image_paths": [str(secret), "../../x", "sub/a.png",
                         "img-1-abc.png"]})
    assert entry["image_paths"] == ["img-1-abc.png"]
    # The exporter refuses an escape even if a raw entry reaches it.
    sc._experiment_attach_dir("exp1", create=True)
    raw = dict(entry, image_paths=[str(secret)])
    srcs, _ = sc.ExperimentExportModal._image_srcs(None, [raw], embed=True)
    assert not srcs
    html = sc._experiment_html_document([raw], image_srcs=srcs,
                                        generated_at="t")
    assert base64.b64encode(secret.read_bytes()).decode() not in html
