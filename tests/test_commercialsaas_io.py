# pyright: reportArgumentType=false, reportCallIssue=false, reportAttributeAccessIssue=false, reportOptionalMemberAccess=false, reportPossiblyUnboundVariable=false
#
# Tests pass `None`-returning helpers through downstream calls
# (`_parse_commercialsaas_history` returning None for malformed XML;
# `_CommercialSaaSHistoryNode.walk` on optional nodes) to verify the
# defensive paths. The project's `pyproject.toml` already excludes
# `tests/**` from pyright; this file-scope pragma keeps the editor /
# harness diagnostics aligned with that policy.
"""
test_commercialsaas_io — low-level .dna packet I/O.

Covers `_iter_commercialsaas_packets`, `_extract_commercialsaas_history_xml`,
`_pack_commercialsaas_history_payload`, `_inject_commercialsaas_history`,
and `_build_commercialsaas_packet`. These are the foundation for full
.dna round-trip including the construction-history packet (0x07),
which neither BioPython nor any other open-source library handles.

Tests use both synthetic byte streams (so the suite stays self-
contained) AND a tiny check against real .dna sample files when
they're available — point ``SPLICECRAFT_DNA_FIXTURES_DIR`` at a
local directory of fixtures to enable the integration tests; they
skip on machines that don't have them.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

import splicecraft as sc
import splicecraft_fileio as _fileio  # .dna codec lives here (Phase D L2); cap reads resolve in this namespace


_SAMPLE_DIR = Path(os.environ.get(
    "SPLICECRAFT_DNA_FIXTURES_DIR", "/nonexistent-default"
))
_SAMPLE_WITH_HISTORY    = _SAMPLE_DIR / "ECXA.dna"
_SAMPLE_WITHOUT_HISTORY = _SAMPLE_DIR / "EcoRI-T5-NcoI-PhiLys-SacI.dna"


def _make_minimal_dna(*packets: tuple[int, bytes]) -> bytes:
    """Build a minimal .dna byte stream from ``(type_byte, payload)``
    pairs. Always starts with a cookie packet so downstream parsers
    don't choke."""
    # 8-byte format magic + 3 unsigned shorts (seqType, exp, imp)
    cookie_payload = sc._COMMERCIALSAAS_COOKIE_MAGIC + bytes([0, 1, 0, 15, 0, 19])
    out = sc._build_commercialsaas_packet(0x09, cookie_payload)
    for type_byte, payload in packets:
        out += sc._build_commercialsaas_packet(type_byte, payload)
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Packet iterator
# ──────────────────────────────────────────────────────────────────────────────


class TestIterCommercialSaaSPackets:
    def test_iterates_simple_two_packet_stream(self):
        data = _make_minimal_dna(
            (0x06, b"<Notes><Type>Synthetic</Type></Notes>"),
        )
        packets = list(sc._iter_commercialsaas_packets(data))
        assert [p[0] for p in packets] == [0x09, 0x06]
        # Length matches payload length.
        assert packets[0][1] == 14   # 8s + 6 bytes
        assert packets[1][1] == len(b"<Notes><Type>Synthetic</Type></Notes>")

    def test_returns_payload_bytes_unchanged(self):
        payload = b"the quick brown fox jumps over the lazy dog"
        data = _make_minimal_dna((0x42, payload))
        ((_t1, _l1, _p1), (t2, l2, p2)) = list(sc._iter_commercialsaas_packets(data))
        assert t2 == 0x42
        assert l2 == len(payload)
        assert p2 == payload

    def test_truncated_header_stops_cleanly(self):
        # A complete cookie packet plus 2 stray bytes (incomplete header).
        data = _make_minimal_dna() + b"\x05\x00"
        packets = list(sc._iter_commercialsaas_packets(data))
        # Only the cookie comes through; the truncated bit is silently
        # dropped (no exception).
        assert len(packets) == 1

    def test_truncated_payload_raises(self):
        # 2026-05-27 (audit-3 H2): a length-overrun used to silently
        # stop the iterator, and `_inject_commercialsaas_history`
        # would then re-emit only the packets seen — irreversibly
        # truncating the file on save. Now raises so the load
        # refuses cleanly instead of producing data loss on round-
        # trip. Header claims 100 bytes but only 5 follow.
        bad = bytes([0x07]) + (100).to_bytes(4, "big") + b"hello"
        data = _make_minimal_dna() + bad
        with pytest.raises(ValueError, match="length overrun"):
            list(sc._iter_commercialsaas_packets(data))

    def test_empty_input_yields_nothing(self):
        assert list(sc._iter_commercialsaas_packets(b"")) == []


# ──────────────────────────────────────────────────────────────────────────────
# Build packet
# ──────────────────────────────────────────────────────────────────────────────


class TestBuildCommercialSaaSPacket:
    def test_round_trips_via_iterator(self):
        payload = b"\x00\x01\x02hello world\xff"
        # 2026-05-27 (audit-3 M9): the iterator now requires a valid
        # cookie packet at file offset 0 so a malformed file lacking
        # the format magic can't load past the iterator. Wrap the
        # arbitrary packet behind a cookie so the iterator accepts
        # the stream; we want to verify packet round-trip, not the
        # cookie guard itself (that's covered separately).
        cookie = sc._build_commercialsaas_cookie_packet()
        pkt = sc._build_commercialsaas_packet(0x42, payload)
        packets = list(sc._iter_commercialsaas_packets(cookie + pkt))
        # 2 packets: cookie + our arbitrary one.
        assert len(packets) == 2
        t, length, p = packets[1]
        assert t == 0x42
        assert length == len(payload)
        assert p == payload

    def test_rejects_oversized_payload(self):
        # 32-bit big-endian length cap is 4 GB - 1.
        with pytest.raises(ValueError):
            sc._build_commercialsaas_packet(0x07, b"x" * (0xFFFFFFFF + 1))

    def test_rejects_out_of_range_type(self):
        with pytest.raises(ValueError):
            sc._build_commercialsaas_packet(-1, b"")
        with pytest.raises(ValueError):
            sc._build_commercialsaas_packet(256, b"")


# ──────────────────────────────────────────────────────────────────────────────
# History XML extract / pack / inject (round-trip core)
# ──────────────────────────────────────────────────────────────────────────────


class TestCommercialSaaSHistoryExtract:
    def test_returns_none_when_no_history_packet(self):
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x06, b"<Notes/>"),
        )
        assert sc._extract_commercialsaas_history_xml(data) is None

    def test_extracts_xml_from_xz_compressed_payload(self):
        """`_extract_commercialsaas_history_xml` decompresses the 0x07
        payload as xz / LZMA — matches the format real CommercialSaaS
        files use. Build a synthetic file that pairs the packer
        with the extractor and verify the round-trip preserves
        the XML byte-for-byte."""
        xml = ('<?xml version="1.0"?>\n'
                '<HistoryTree>'
                '<Node name="test.dna" type="DNA" seqLen="100" '
                'circular="1" operation="insertFragment"/>'
                '</HistoryTree>')
        payload = sc._pack_commercialsaas_history_payload(xml)
        data = _make_minimal_dna((0x07, payload))
        recovered = sc._extract_commercialsaas_history_xml(data)
        assert recovered == xml

    def test_rejects_invalid_xz(self):
        data = _make_minimal_dna((0x07, b"not actually xz"))
        with pytest.raises(ValueError) as exc:
            sc._extract_commercialsaas_history_xml(data)
        assert "xz" in str(exc.value).lower()

    def test_rejects_decompression_bomb(self, monkeypatch):
        """A crafted .dna with a tiny 0x07 packet that decompresses to
        gigabytes must NOT silently allocate that much. The cap lives
        at `_COMMERCIALSAAS_HISTORY_MAX_XML`; lower it for the test so we
        don't actually have to build a multi-GB file."""
        monkeypatch.setattr(_fileio, "_COMMERCIALSAAS_HISTORY_MAX_XML", 100)
        big_xml = "<HistoryTree>" + ("X" * 500) + "</HistoryTree>"
        payload = sc._pack_commercialsaas_history_payload(big_xml)
        data = _make_minimal_dna((0x07, payload))
        with pytest.raises(ValueError) as exc:
            sc._extract_commercialsaas_history_xml(data)
        assert "too large" in str(exc.value).lower()

    def test_rejects_streaming_bomb_aborts_early(self, monkeypatch):
        """Regression guard for 2026-05-06 fix: the previous implementation
        called `_lzma.decompress(payload)` to completion BEFORE checking the
        size — a 10 MB compressed payload that expands to 50 GB would OOM
        before the cap even ran. The fix uses `LZMADecompressor` with
        `max_length`, so the decoder aborts at the cap and never
        materialises the full plaintext."""
        # 10 KB compressed cap with a payload that decompresses to ~100 KB.
        # The old code would allocate the full 100 KB before checking
        # — verify behaviour is now bounded by the cap.
        monkeypatch.setattr(_fileio, "_COMMERCIALSAAS_HISTORY_MAX_XML", 1_000)
        big_xml = "<HistoryTree>" + ("Z" * 100_000) + "</HistoryTree>"
        payload = sc._pack_commercialsaas_history_payload(big_xml)
        data = _make_minimal_dna((0x07, payload))
        with pytest.raises(ValueError, match="too large"):
            sc._extract_commercialsaas_history_xml(data)


class TestCommercialSaaSHistoryInject:
    def test_replaces_existing_history(self):
        old_xml = "<HistoryTree><Node name='old.dna'/></HistoryTree>"
        new_xml = "<HistoryTree><Node name='new.dna'/></HistoryTree>"
        old_payload = sc._pack_commercialsaas_history_payload(old_xml)
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x07, old_payload),
            (0x06, b"<Notes/>"),
        )
        out = sc._inject_commercialsaas_history(data, new_xml)
        # Recovered XML matches the new one.
        assert sc._extract_commercialsaas_history_xml(out) == new_xml
        # Other packets preserved verbatim.
        types = [p[0] for p in sc._iter_commercialsaas_packets(out)]
        assert types == [0x09, 0x00, 0x07, 0x06]

    def test_inserts_when_no_existing_history(self):
        new_xml = "<HistoryTree><Node name='fresh.dna'/></HistoryTree>"
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x06, b"<Notes/>"),
        )
        out = sc._inject_commercialsaas_history(data, new_xml)
        assert sc._extract_commercialsaas_history_xml(out) == new_xml
        # New packet inserted right after the cookie.
        types = [p[0] for p in sc._iter_commercialsaas_packets(out)]
        assert types == [0x09, 0x07, 0x00, 0x06]

    def test_strips_history_when_passed_none(self):
        old_xml = "<HistoryTree><Node name='old.dna'/></HistoryTree>"
        old_payload = sc._pack_commercialsaas_history_payload(old_xml)
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x07, old_payload),
            (0x06, b"<Notes/>"),
        )
        out = sc._inject_commercialsaas_history(data, None)
        assert sc._extract_commercialsaas_history_xml(out) is None
        types = [p[0] for p in sc._iter_commercialsaas_packets(out)]
        assert 0x07 not in types

    def test_preserves_unknown_packets_byte_exact(self):
        """`_inject_commercialsaas_history` must NOT mangle packets it
        doesn't understand — the whole point of the splice approach
        is that we can round-trip files with packets we haven't yet
        decoded (alignments, custom enzymes, etc.) without losing
        them or corrupting their bytes."""
        weird_payload = b"\x00\x01\x02\x03binary garbage\xff\xfe"
        data = _make_minimal_dna(
            (0x42, weird_payload),
            (0x99, b"another unknown"),
        )
        new_xml = "<HistoryTree/>"
        out = sc._inject_commercialsaas_history(data, new_xml)
        # Re-extract the unknown packets and check byte-equality.
        recovered = {t: p for t, _l, p in sc._iter_commercialsaas_packets(out)}
        assert recovered[0x42] == weird_payload
        assert recovered[0x99] == b"another unknown"


# ──────────────────────────────────────────────────────────────────────────────
# Real-file smoke checks — gated on sample files being present
# ──────────────────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────────────────
# History tree parser + serializer
# ──────────────────────────────────────────────────────────────────────────────


_SAMPLE_HISTORY_XML = """<?xml version="1.0" encoding="UTF-8"?>
<HistoryTree>
  <Node name="result.dna" type="DNA" seqLen="2603" strandedness="double"
        ID="7" circular="1" operation="insertFragment">
    <RegeneratedSite name="KpnI" pos="869" siteCount="1"/>
    <RegeneratedSite name="HindIII" pos="1961" siteCount="1"/>
    <InputSummary manipulation="replace" name1="HindIII" name2="KpnI"
                  val1="1826" val2="869" siteCount1="1" siteCount2="1"/>
    <Node name="parent_a.dna" type="DNA" seqLen="2468" strandedness="double"
          ID="2" circular="1" resurrectable="1" operation="insertFragment">
      <Node name="grandparent.dna" type="DNA" seqLen="1500"
            strandedness="double" ID="0" circular="0"
            operation="insertFragment"/>
    </Node>
    <Node name="parent_b.dna" type="DNA" seqLen="500" strandedness="double"
          ID="3" circular="0" resurrectable="1" operation="insert"/>
  </Node>
</HistoryTree>"""


class TestParseCommercialSaaSHistory:
    def test_parses_top_level_node(self):
        root = sc._parse_commercialsaas_history(_SAMPLE_HISTORY_XML)
        assert root is not None
        assert root.name == "result.dna"
        assert root.operation == "insertFragment"
        assert root.seq_len == 2603
        assert root.circular is True
        assert root.node_id == 7

    def test_typed_getters_for_attributes(self):
        root = sc._parse_commercialsaas_history(_SAMPLE_HISTORY_XML)
        assert root.regenerated_sites == [
            {"name": "KpnI",    "pos": 869,  "siteCount": 1},
            {"name": "HindIII", "pos": 1961, "siteCount": 1},
        ]
        summaries = root.input_summaries
        assert len(summaries) == 1
        assert summaries[0]["manipulation"] == "replace"
        assert summaries[0]["val1"] == 1826

    def test_parents_returns_direct_children(self):
        root = sc._parse_commercialsaas_history(_SAMPLE_HISTORY_XML)
        parents = root.parents
        assert [p.name for p in parents] == ["parent_a.dna", "parent_b.dna"]
        # Resurrectable flag captured on parent nodes.
        assert all(p.resurrectable for p in parents)
        # The top-level result node is NOT resurrectable.
        assert root.resurrectable is False

    def test_walk_pre_order_depth_first(self):
        root = sc._parse_commercialsaas_history(_SAMPLE_HISTORY_XML)
        names = [n.name for n in root.walk()]
        assert names == [
            "result.dna",
            "parent_a.dna",
            "grandparent.dna",
            "parent_b.dna",
        ]

    def test_returns_none_on_empty_xml(self):
        assert sc._parse_commercialsaas_history("") is None
        assert sc._parse_commercialsaas_history("   \n  ") is None

    def test_returns_none_when_history_tree_has_no_nodes(self):
        assert sc._parse_commercialsaas_history(
            '<?xml version="1.0"?><HistoryTree></HistoryTree>'
        ) is None

    def test_rejects_wrong_root_element(self):
        with pytest.raises(ValueError) as exc:
            sc._parse_commercialsaas_history(
                '<?xml version="1.0"?><WrongRoot/>'
            )
        assert "HistoryTree" in str(exc.value)

    def test_rejects_malformed_xml(self):
        with pytest.raises(ValueError):
            sc._parse_commercialsaas_history("<HistoryTree><Node oops")

    def test_malformed_int_attrs_default_to_zero(self):
        """A node with unparseable seqLen must NOT raise — return 0
        and let downstream code decide what to do. Defends against
        future CommercialSaaS versions that might change attribute formats."""
        xml = ('<?xml version="1.0"?><HistoryTree>'
                '<Node name="x" seqLen="not-a-number" '
                'circular="1" operation="insertFragment"/>'
                '</HistoryTree>')
        root = sc._parse_commercialsaas_history(xml)
        assert root is not None
        assert root.seq_len == 0
        assert root.node_id == 0


class TestSerializeCommercialSaaSHistory:
    def test_round_trip_via_parse_then_serialize(self):
        root = sc._parse_commercialsaas_history(_SAMPLE_HISTORY_XML)
        serialized = sc._serialize_commercialsaas_history(root)
        # Re-parse: must give back an equivalent tree.
        root2 = sc._parse_commercialsaas_history(serialized)
        assert root2 is not None
        assert root2.name           == root.name
        assert root2.operation      == root.operation
        assert root2.seq_len        == root.seq_len
        assert [p.name for p in root2.parents] == \
               [p.name for p in root.parents]

    def test_emits_xml_declaration(self):
        """CommercialSaaS's history XML always opens with the standard
        `<?xml ...?>` declaration. Make sure our serializer keeps it
        — file parsers downstream sometimes refuse declaration-less
        XML."""
        root = sc._parse_commercialsaas_history(_SAMPLE_HISTORY_XML)
        out = sc._serialize_commercialsaas_history(root)
        assert out.startswith('<?xml version="1.0" encoding="UTF-8"?>')

    def test_serializes_empty_tree(self):
        out = sc._serialize_commercialsaas_history(None)
        # Just the HistoryTree shell, valid + parseable.
        assert "<HistoryTree" in out
        assert sc._parse_commercialsaas_history(out) is None


class TestCommercialSaaSHistoryNodeMutation:
    def test_new_creates_canonical_attributes(self):
        node = sc._CommercialSaaSHistoryNode.new(
            name="fresh.dna", seq_len=1000, circular=True,
            operation="insertFragment", node_id=42,
        )
        assert node.name        == "fresh.dna"
        assert node.seq_len     == 1000
        assert node.circular    is True
        assert node.node_id     == 42
        assert node.operation   == "insertFragment"
        # Type defaulted to DNA.
        assert node.element.get("type") == "DNA"
        # Strandedness defaulted to double.
        assert node.element.get("strandedness") == "double"

    def test_add_parent_appends_to_tree(self):
        root = sc._CommercialSaaSHistoryNode.new(
            name="result.dna", seq_len=2000, circular=True,
            operation="insertFragment",
        )
        parent = sc._CommercialSaaSHistoryNode.new(
            name="parent.dna", seq_len=1500, circular=True,
            operation="insertFragment",
        )
        root.add_parent(parent)
        assert [p.name for p in root.parents] == ["parent.dna"]

    def test_add_regenerated_site(self):
        root = sc._CommercialSaaSHistoryNode.new(
            name="r", seq_len=100, circular=True,
            operation="insertFragment",
        )
        root.add_regenerated_site("EcoRI", 50, site_count=1)
        sites = root.regenerated_sites
        assert sites == [{"name": "EcoRI", "pos": 50, "siteCount": 1}]

    def test_add_input_summary(self):
        root = sc._CommercialSaaSHistoryNode.new(
            name="r", seq_len=100, circular=True,
            operation="insertFragment",
        )
        root.add_input_summary(manipulation="replace", name1="EcoRI",
                                  name2="BamHI", val1=10, val2=50)
        summaries = root.input_summaries
        assert len(summaries) == 1
        assert summaries[0]["manipulation"] == "replace"
        assert summaries[0]["name1"]        == "EcoRI"
        assert summaries[0]["val1"]         == 10

    def test_built_tree_round_trips_through_serialize(self):
        """Build a fresh tree from scratch, serialise it, parse it
        back. The parser must see every parent + regenerated site +
        input summary we attached."""
        root = sc._CommercialSaaSHistoryNode.new(
            name="construct.dna", seq_len=5000, circular=True,
            operation="insertFragment", node_id=1,
        )
        root.add_regenerated_site("EcoRI", 100)
        root.add_input_summary(manipulation="replace",
                                  name1="EcoRI", name2="BamHI",
                                  val1=100, val2=600)
        parent_a = sc._CommercialSaaSHistoryNode.new(
            name="vec.dna", seq_len=3500, circular=True,
            operation="insertFragment", node_id=2,
        )
        parent_b = sc._CommercialSaaSHistoryNode.new(
            name="ins.dna", seq_len=1500, circular=False,
            operation="insertFragment", node_id=3,
        )
        root.add_parent(parent_a)
        root.add_parent(parent_b)
        out = sc._serialize_commercialsaas_history(root)
        re_parsed = sc._parse_commercialsaas_history(out)
        assert re_parsed is not None
        assert re_parsed.name == "construct.dna"
        parent_names = [p.name for p in re_parsed.parents]
        assert parent_names == ["vec.dna", "ins.dna"]
        # Sites + input summaries also survived.
        assert re_parsed.regenerated_sites == [
            {"name": "EcoRI", "pos": 100, "siteCount": 1},
        ]
        assert re_parsed.input_summaries[0]["manipulation"] == "replace"

    def test_full_chain_dna_to_dna_round_trip(self):
        """End-to-end: take an in-memory .dna byte stream, extract
        history, parse to typed tree, append a new step, re-serialise,
        re-inject, re-extract — must give the modified tree back."""
        # Start with a synthetic .dna that has a history.
        original_xml = _SAMPLE_HISTORY_XML
        original_payload = sc._pack_commercialsaas_history_payload(original_xml)
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x07, original_payload),
        )
        # Parse → modify → serialise.
        xml = sc._extract_commercialsaas_history_xml(data)
        root = sc._parse_commercialsaas_history(xml)
        # Append a new step via add_input_summary.
        root.add_input_summary(manipulation="ligate-via-splicecraft",
                                  name1="EcoRI", name2="BamHI",
                                  val1=10, val2=20)
        new_xml = sc._serialize_commercialsaas_history(root)
        out = sc._inject_commercialsaas_history(data, new_xml)
        # Round-trip read.
        recovered_xml = sc._extract_commercialsaas_history_xml(out)
        recovered_root = sc._parse_commercialsaas_history(recovered_xml)
        # New step is in.
        all_summaries = recovered_root.input_summaries
        labels = [s["manipulation"] for s in all_summaries]
        assert "ligate-via-splicecraft" in labels


class TestHistoryDateStamp:
    """Every history step SpliceCraft builds going forward carries a `date`
    attribute, rendered in the History tab in the universal slash-free
    "JUN 9 2026 14:30" form. Reconstructed lineage / starting material stays
    undated (no fabricated date). Imported CommercialSaaS "YYYY.MM.DD" dates
    render too."""

    def test_human_dt_iso_with_time(self):
        assert sc._history_human_dt("2026-06-09T14:30") == "JUN 9 2026 14:30"
        assert sc._history_human_dt("2026-06-09T09:05") == "JUN 9 2026 09:05"

    def test_human_dt_date_only_and_commercialsaas_dotted(self):
        # Date-only (no time) and the CommercialSaaS Notes "YYYY.MM.DD" form.
        assert sc._history_human_dt("2026-12-25") == "DEC 25 2026"
        assert sc._history_human_dt("2026.01.03") == "JAN 3 2026"

    def test_human_dt_rejects_garbage_and_empty(self):
        for bad in ("", None, "not a date", "2026", "abcd-ef-gh",
                    "2026-13-01", "2026-06-45"):
            assert sc._history_human_dt(bad) == "", repr(bad)

    def test_human_dt_drops_bad_time_keeps_date(self):
        # Hour/minute out of range → keep the date, drop the time.
        assert sc._history_human_dt("2026-06-09T99:99") == "JUN 9 2026"

    def test_human_dt_tolerates_non_string(self):
        # Defensive str() coercion — never raises on an odd caller value.
        assert sc._history_human_dt(20260609) == ""        # int → no match
        assert sc._history_human_dt(["2026-06-09"]) == ""  # list → no match

    def test_now_str_round_trips_through_formatter(self):
        now = sc._history_now_str()
        assert sc._history_human_dt(now), f"now-stamp unrenderable: {now!r}"

    def test_new_stores_and_omits_date(self):
        dated = sc._CommercialSaaSHistoryNode.new(
            name="x.dna", seq_len=10, circular=True,
            operation="insertFragment", date="2026-06-09T14:30")
        assert dated.date == "2026-06-09T14:30"
        # Default: no date attr (reconstructed lineage / starting material).
        undated = sc._CommercialSaaSHistoryNode.new(
            name="y.dna", seq_len=10, circular=True,
            operation="insertFragment")
        assert undated.date == ""
        assert undated.element.get("date") is None

    def test_date_survives_serialize_round_trip(self):
        root = sc._CommercialSaaSHistoryNode.new(
            name="r.dna", seq_len=10, circular=True,
            operation="insertFragment", date="2026-06-09T14:30")
        reparsed = sc._parse_commercialsaas_history(
            sc._serialize_commercialsaas_history(root))
        assert reparsed.date == "2026-06-09T14:30"

    def test_fresh_builders_stamp_a_date(self):
        # The origin builder stamps "now"; the parsed root renders a date.
        xml = sc._build_origin_history_xml(
            name="p", seq_len=100, circular=True, source="library:abc")
        root = sc._parse_commercialsaas_history(xml)
        assert root is not None and sc._history_human_dt(root.date)

    def test_parent_leaf_is_undated(self):
        # `_parent_node_for_entry` reconstructs PRE-EXISTING material → no
        # fabricated date.
        leaf = sc.TraditionalCloningPane._parent_node_for_entry(
            {"name": "old-vec", "id": "v1", "size": 3000})
        assert leaf is not None and leaf.date == ""

    def test_tree_label_and_detail_show_date(self):
        dated = sc._CommercialSaaSHistoryNode.new(
            name="x.dna", seq_len=10, circular=True,
            operation="insertFragment", date="2026-06-09T14:30")
        assert "JUN 9 2026 14:30" in sc._history_tree_label(dated)
        assert any("JUN 9 2026 14:30" in ln
                   for ln in sc._history_detail_lines(dated))
        undated = sc._CommercialSaaSHistoryNode.new(
            name="y.dna", seq_len=10, circular=True,
            operation="insertFragment")
        assert "JUN" not in sc._history_tree_label(undated)

    # ── .dna file-date import ────────────────────────────────────────────

    def test_extract_file_date_from_notes_created(self):
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x06, b"<Notes><Created>2026.06.09</Created></Notes>"))
        assert sc._extract_commercialsaas_file_date(data) == "2026-06-09"

    def test_extract_file_date_none_without_notes(self):
        data = _make_minimal_dna((0x00, bytes([0x01]) + b"ATGC"))
        assert sc._extract_commercialsaas_file_date(data) is None

    def test_stamp_root_date_only_when_dateless(self):
        xml = "<HistoryTree><Node name='x.dna' ID='0'/></HistoryTree>"
        out = sc._stamp_history_root_date(xml, "2026-06-09")
        assert sc._parse_commercialsaas_history(out).date == "2026-06-09"
        # An already-dated root is left untouched (verbatim string back).
        dated = sc._serialize_commercialsaas_history(
            sc._CommercialSaaSHistoryNode.new(
                name="x.dna", seq_len=1, circular=True,
                operation="insertFragment", date="2020-01-01"))
        assert sc._stamp_history_root_date(dated, "2026-06-09") == dated
        # No date to add → unchanged.
        assert sc._stamp_history_root_date(xml, None) == xml

    def test_bulk_import_stamps_file_date_on_root(self, tmp_path):
        # Dateless history root + a Notes file date → the bulk-import
        # extractor stamps the file date onto the root, rendered "JUN 6 2025".
        hist = "<HistoryTree><Node name='child.dna' ID='0'/></HistoryTree>"
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x06, b"<Notes><LastModified>2025.06.06</LastModified></Notes>"),
            (0x07, sc._pack_commercialsaas_history_payload(hist)),
        )
        path = tmp_path / "dated.dna"
        path.write_bytes(data)
        root = sc._parse_commercialsaas_history(
            sc._try_extract_history_xml_from_dna_path(path))
        assert root.date == "2025-06-06"
        assert sc._history_human_dt(root.date) == "JUN 6 2025"


# ──────────────────────────────────────────────────────────────────────────────
# Library-entry import wiring (Phase 4a)
# ──────────────────────────────────────────────────────────────────────────────


class TestLibraryEntryHistory:
    def test_record_to_library_entry_omits_history_when_none(
            self, tiny_record):
        from pathlib import Path as _Path
        entry = sc._record_to_library_entry(tiny_record,
                                              _Path("foo.gb"),
                                              history_xml=None)
        assert "history_xml" not in entry

    def test_record_to_library_entry_attaches_history(self, tiny_record):
        from pathlib import Path as _Path
        xml = "<HistoryTree><Node name='x'/></HistoryTree>"
        entry = sc._record_to_library_entry(tiny_record,
                                              _Path("foo.dna"),
                                              history_xml=xml)
        assert entry["history_xml"] == xml

    def test_extract_history_returns_none_for_non_dna(self, tmp_path):
        gb = tmp_path / "x.gb"
        gb.write_bytes(b"LOCUS dummy\n//\n")
        assert sc._try_extract_history_xml_from_dna_path(gb) is None

    def test_extract_history_returns_none_for_history_less_dna(
            self, tmp_path):
        # Build a minimal .dna without a 0x07 packet.
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
        )
        path = tmp_path / "no_history.dna"
        path.write_bytes(data)
        assert sc._try_extract_history_xml_from_dna_path(path) is None

    def test_extract_history_returns_xml_for_dna_with_history(
            self, tmp_path):
        xml = "<HistoryTree><Node name='child.dna'/></HistoryTree>"
        payload = sc._pack_commercialsaas_history_payload(xml)
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x07, payload),
        )
        path = tmp_path / "with_history.dna"
        path.write_bytes(data)
        assert sc._try_extract_history_xml_from_dna_path(path) == xml

    def test_extract_history_swallows_malformed_payload(self, tmp_path):
        """A .dna with a malformed history packet must NOT crash the
        bulk-import path — log + return None so the rest of the
        plasmid still imports."""
        data = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x07, b"definitely not xz"),
        )
        path = tmp_path / "bad_history.dna"
        path.write_bytes(data)
        assert sc._try_extract_history_xml_from_dna_path(path) is None

    def test_extract_history_swallows_unreadable_file(self, tmp_path):
        # Path to a file that doesn't exist.
        ghost = tmp_path / "ghost.dna"
        assert sc._try_extract_history_xml_from_dna_path(ghost) is None


# ──────────────────────────────────────────────────────────────────────────────
# History viewer modal (Phase 4c)
# ──────────────────────────────────────────────────────────────────────────────


class TestHistoryViewerModal:
    async def test_renders_root_and_parents_in_tree(self, tiny_record):
        """Build a 3-deep history tree, push the modal, verify the
        Textual `Tree` widget contains a node per history step."""
        from textual.widgets import Tree as _TreeWidget
        # Construct: root → parent → grandparent.
        gp = sc._CommercialSaaSHistoryNode.new(
            name="grandparent.dna", seq_len=1000, circular=True,
            operation="insertFragment", node_id=2,
        )
        parent = sc._CommercialSaaSHistoryNode.new(
            name="parent.dna", seq_len=2500, circular=True,
            operation="insertFragment", node_id=1,
        )
        parent.add_parent(gp)
        root = sc._CommercialSaaSHistoryNode.new(
            name="result.dna", seq_len=5000, circular=True,
            operation="insertFragment", node_id=0,
        )
        root.add_parent(parent)
        # Drive the modal headlessly via PlasmidApp.run_test.
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        # The modal's compose doesn't depend on the loaded record;
        # using `tiny_record` just so the app mounts cleanly.
        app = _build_app(tiny_record, isolated_library=None)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.HistoryViewerModal("test-plasmid", root)
            await app.push_screen(modal)
            await pilot.pause()
            tree = modal.query_one("#hist-tree", _TreeWidget)
            # Walk the Textual tree, collect labels.
            labels: list[str] = []
            def _walk(n):
                # Tree.root has no useful label of its own (root is hidden).
                if n is not tree.root:
                    labels.append(str(n.label))
                for c in n.children:
                    _walk(c)
            _walk(tree.root)
            # Three nodes total: result + parent + grandparent.
            assert len(labels) == 3, labels
            # Labels strip the cosmetic `.dna` suffix (de-noised render,
            # 2026-05-30) — assert on the bare names.
            assert any("result" in lab for lab in labels)
            assert any("grandparent" in lab for lab in labels)
            assert any(lab.split()[0] == "parent"
                       or "parent" in lab for lab in labels)

    async def test_h_key_with_history_opens_modal(
            self, tiny_record, isolated_library):
        """Library panel binding `h` → app handler → push
        HistoryViewerModal when the entry has `history_xml`."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        # Seed library with one entry that has history.
        history_xml = sc._serialize_commercialsaas_history(
            sc._CommercialSaaSHistoryNode.new(
                name="seeded.dna", seq_len=1000, circular=True,
                operation="insertFragment", node_id=0,
            )
        )
        sc._save_library([{
            "id":   "seeded",
            "name": "seeded",
            "size": 1000,
            "history_xml": history_xml,
        }])
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            # Switch to plasmids view + cursor on the seeded row.
            from textual.widgets import DataTable
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "plasmids"
            await pilot.pause()
            try:
                tbl = lib.query_one("#lib-table", DataTable)
            except Exception:
                pytest.skip("library panel layout differs in this build")
            # Use cursor moves to land on the seeded row (other rows may
            # have been auto-loaded by the app).
            tbl.move_cursor(row=0)
            await pilot.pause()
            # Trigger the action directly (bypasses keyboard simulation
            # which is fragile in Textual run_test).
            lib.action_request_history()
            await pilot.pause()
            # Modal should now be on the screen stack.
            assert any(isinstance(s, sc.HistoryViewerModal)
                        for s in app.screen_stack)

    async def test_h_key_without_history_notifies(
            self, tiny_record, isolated_library):
        """Library entry without `history_xml` → user-friendly
        notify, NOT the modal."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        sc._save_library([{
            "id":   "no_hist",
            "name": "no_hist",
            "size": 1000,
            # No history_xml field.
        }])
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            from textual.widgets import DataTable
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "plasmids"
            await pilot.pause()
            try:
                tbl = lib.query_one("#lib-table", DataTable)
            except Exception:
                pytest.skip("library panel layout differs in this build")
            tbl.move_cursor(row=0)
            await pilot.pause()
            lib.action_request_history()
            await pilot.pause()
            # No modal pushed.
            assert not any(isinstance(s, sc.HistoryViewerModal)
                            for s in app.screen_stack)


# ──────────────────────────────────────────────────────────────────────────────
# Full-screen HistoryScreen + menu / key wiring
# ──────────────────────────────────────────────────────────────────────────────


class TestHistoryScreen:
    """Cover the full-screen viewer (`HistoryScreen`) and the
    `action_show_history` path that drives it from the `History` menu
    tab, F5, and Ctrl+H. Regression guard for 2026-05-11: the panel
    promotion (modal → fullscreen Screen) must not lose either the
    tree-rendering behaviour from `HistoryViewerModal` or the
    "current-record / library-entry lookup" semantics."""

    async def test_renders_root_and_parents_in_tree(self, tiny_record):
        """Same data-shape assertion as the modal: a 3-deep history
        tree round-trips through `HistoryScreen.on_mount` into the
        Textual Tree widget."""
        from textual.widgets import Tree as _TreeWidget
        gp = sc._CommercialSaaSHistoryNode.new(
            name="grandparent.dna", seq_len=1000, circular=True,
            operation="insertFragment", node_id=2,
        )
        parent = sc._CommercialSaaSHistoryNode.new(
            name="parent.dna", seq_len=2500, circular=True,
            operation="insertFragment", node_id=1,
        )
        parent.add_parent(gp)
        root = sc._CommercialSaaSHistoryNode.new(
            name="result.dna", seq_len=5000, circular=True,
            operation="insertFragment", node_id=0,
        )
        root.add_parent(parent)
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library=None)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            screen = sc.HistoryScreen("test-plasmid", root)
            await app.push_screen(screen)
            await pilot.pause()
            tree = screen.query_one("#hist-scr-tree", _TreeWidget)
            labels: list[str] = []
            def _walk(n):
                if n is not tree.root:
                    labels.append(str(n.label))
                for c in n.children:
                    _walk(c)
            _walk(tree.root)
            assert len(labels) == 3, labels
            # Labels strip the cosmetic `.dna` suffix (de-noised render,
            # 2026-05-30) — assert on the bare names.
            assert any("result" in lab for lab in labels)
            assert any("grandparent" in lab for lab in labels)
            assert any(lab.split()[0] == "parent"
                       or "parent" in lab for lab in labels)

    async def test_action_show_history_pushes_screen_when_loaded(
            self, tiny_record, isolated_library):
        """`PlasmidApp.action_show_history` looks the loaded record up
        in the active library by id and, if it has `history_xml`,
        pushes `HistoryScreen` (not the legacy modal). Wired to F6 /
        Ctrl+H / the `History` top-bar menu tab.

        Seeds the library AFTER the app mounts so the auto-persist
        of the preloaded record doesn't overwrite our seed — the
        app's load path saves the record without `history_xml` if it
        wasn't already on disk."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        history_xml = sc._serialize_commercialsaas_history(
            sc._CommercialSaaSHistoryNode.new(
                name=tiny_record.id, seq_len=len(tiny_record.seq),
                circular=True, operation="insertFragment", node_id=0,
            )
        )
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            # Stamp history_xml onto the auto-persisted entry.
            entries = sc._load_library()
            for e in entries:
                if e.get("id") == tiny_record.id:
                    e["history_xml"] = history_xml
            sc._save_library(entries)
            app.action_show_history()
            await pilot.pause()
            assert any(isinstance(s, sc.HistoryScreen)
                        for s in app.screen_stack)

    async def test_action_show_history_notifies_when_no_history(
            self, tiny_record, isolated_library):
        """A loaded record whose library entry has NO `history_xml` — a
        legacy entry pre-dating origin-history — makes
        `action_show_history` notify rather than push a blank screen.

        As of the origin-history change the auto-persist path now
        BACKFILLS a minimal root, so we strip it here to exercise the
        genuine no-history notify branch (legacy / history-less entries
        still hit it)."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            # Simulate a legacy entry with no recorded history.
            lib = sc._load_library()
            for e in lib:
                e.pop("history_xml", None)
            sc._save_library(lib)
            app.action_show_history()
            await pilot.pause()
            assert not any(isinstance(s, sc.HistoryScreen)
                            for s in app.screen_stack)

    def test_history_in_menu_bar(self):
        """The `History` menu tab is registered alongside `Features`
        as a direct-open (no-dropdown) entry."""
        assert "History" in sc.MenuBar.MENUS

    def test_app_has_history_and_restore_bindings(self):
        """F6 + Ctrl+H route to `show_history`; F5 restores the
        multi-panel view via `focus_panel_all`. Reverted from the
        0.7.11.0 swap that broke a user's muscle memory (GH #15) —
        F5 returns to "all panels" again as the F1-F4 inverse. Surface
        assertions only — the click and modal-stack tests above cover
        behavior."""
        keys_to_action = {
            (b.key, b.action) for b in sc.PlasmidApp.BINDINGS
        }
        assert ("f5", "focus_panel_all") in keys_to_action
        assert ("f6", "show_history") in keys_to_action
        assert ("ctrl+h", "show_history") in keys_to_action


class TestHistoryScreenHardening:
    """Edge-case coverage for `HistoryScreen` — the visible-from-the-
    outside-world history viewer. Catches recursion-limit blowouts on
    deep trees, Rich-markup injection from XML-controlled fields,
    overflow on long names / many regenerated sites / many parents,
    and the `Tree` iteration helpers (`_build_tree`, `expand_all`,
    `collapse_all`). Regression guard for 2026-05-11 hardening."""

    def _make_chain(self, depth: int) -> "sc._CommercialSaaSHistoryNode":
        """Build a single-chain history tree N levels deep, returning
        the root. Each node hangs off its child's `parents` list."""
        root = sc._CommercialSaaSHistoryNode.new(
            name="leaf_0.dna", seq_len=100, circular=True,
            operation="insertFragment", node_id=0,
        )
        cur = root
        for i in range(1, depth):
            up = sc._CommercialSaaSHistoryNode.new(
                name=f"gen_{i}.dna", seq_len=100 + i, circular=True,
                operation="insertFragment", node_id=i,
            )
            cur.add_parent(up)
            cur = up
        return root

    def test_node_count_helper_counts_all_nodes(self):
        """`_history_node_count` must walk iteratively (mirrors
        `_CommercialSaaSHistoryNode.walk`) — a deep chain shouldn't
        require a deep recursion."""
        root = self._make_chain(200)
        assert sc._history_node_count(root) == 200

    def test_node_count_root_only(self):
        node = sc._CommercialSaaSHistoryNode.new(
            name="solo.dna", seq_len=10, circular=True,
            operation="insertFragment", node_id=0,
        )
        assert sc._history_node_count(node) == 1

    def test_tree_label_truncates_long_name(self):
        """Names longer than `_HISTORY_LABEL_NAME_MAX` are truncated
        with an ellipsis so the tree column stays usable on hostile
        / pathologically-long names."""
        long_name = "x" * (sc._HISTORY_LABEL_NAME_MAX + 50)
        node = sc._CommercialSaaSHistoryNode.new(
            name=long_name, seq_len=100, circular=True,
            operation="insertFragment", node_id=0,
        )
        label = sc._history_tree_label(node)
        # No bare `[` in the rendered string (markup injection would
        # show up as a literal Rich tag mid-label).
        assert long_name not in label
        # Display name is bounded.
        assert "…" in label

    def test_tree_label_escapes_markup_in_name(self):
        """Rich markup characters in the node name must be escaped
        — a node named ``pUC[red]boom[/red]`` should NOT paint red
        text in the tree. `rich.markup.escape` rewrites the leading
        `[` as `\\[`, so the output carries `\\[red]` (literal) rather
        than the active tag form. We assert the escape sentinel was
        applied; Rich's rendering pipeline treats `\\[red]` as plain
        text and will not turn on the `red` style."""
        node = sc._CommercialSaaSHistoryNode.new(
            name="pUC[red]boom[/red]", seq_len=100, circular=True,
            operation="insertFragment", node_id=0,
        )
        label = sc._history_tree_label(node)
        assert "\\[red]" in label, label

    def test_tree_label_empty_fields_render_placeholders(self):
        """An XML node with no operation / no name produces a readable
        row rather than a whitespace-only label. De-noised render
        (2026-05-30): an empty name still shows ``(unnamed)``, but an
        empty operation / parent-less node no longer prints a
        ``(no operation)`` placeholder on every row — the tag is simply
        omitted (it was per-row noise; the detail pane still records it)."""
        import xml.etree.ElementTree as ET
        # Build a node with empty operation + name attributes — bypassing
        # `.new()` because it always sets them.
        el = ET.Element("Node")
        el.set("name", "")
        el.set("seqLen", "0")
        el.set("circular", "0")
        el.set("operation", "")
        node = sc._CommercialSaaSHistoryNode(el)
        label = sc._history_tree_label(node)
        assert "(unnamed)" in label
        assert label.strip()                 # never whitespace-only
        assert "no operation" not in label   # placeholder dropped

    def test_title_truncates_long_plasmid_name(self):
        """A library plasmid with a 200-char name shouldn't blow the
        title bar — HistoryScreen ellipsises the title up front."""
        root = sc._CommercialSaaSHistoryNode.new(
            name="x.dna", seq_len=10, circular=True,
            operation="insertFragment", node_id=0,
        )
        long_title = "Z" * 200
        screen = sc.HistoryScreen(long_title, root)
        assert len(screen._title) <= sc._HISTORY_TITLE_NAME_MAX
        assert screen._title.endswith("…")

    async def test_screen_handles_deep_history_without_recursion(
            self, tiny_record):
        """Pushing a HistoryScreen with a 200-deep chain must mount
        without tripping CPython's recursion limit. Smoke for
        `_build_tree`'s iterative shape."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from textual.widgets import Tree as _TreeWidget
        root = self._make_chain(200)
        app = _build_app(tiny_record, isolated_library=None)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            screen = sc.HistoryScreen("deep", root)
            await app.push_screen(screen)
            await pilot.pause()
            tree = screen.query_one("#hist-scr-tree", _TreeWidget)
            # Count labels — must equal node count.
            cnt = 0
            stack = [tree.root]
            while stack:
                n = stack.pop()
                if n is not tree.root:
                    cnt += 1
                stack.extend(n.children)
            assert cnt == 200

    async def test_detail_pane_truncates_many_sites_and_parents(
            self, tiny_record):
        """A node with 50 regenerated sites + 50 parents renders only
        the cap (`_HISTORY_DETAIL_LIST_MAX`) inline, with a
        `(+N more)` suffix so the detail column doesn't blow up."""
        root = sc._CommercialSaaSHistoryNode.new(
            name="result.dna", seq_len=5000, circular=True,
            operation="insertFragment", node_id=0,
        )
        for i in range(50):
            root.add_regenerated_site(f"Enz{i}", pos=i, site_count=1)
        for i in range(50):
            root.add_parent(sc._CommercialSaaSHistoryNode.new(
                name=f"parent_{i}.dna", seq_len=100, circular=True,
                operation="insertFragment", node_id=i + 1,
            ))
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library=None)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            screen = sc.HistoryScreen("packed", root)
            await app.push_screen(screen)
            await pilot.pause()
            # Force a detail render by selecting the root.
            from textual.widgets import Tree as _TreeWidget
            tree = screen.query_one("#hist-scr-tree", _TreeWidget)
            top = tree.root.children[0]
            tree.select_node(top)
            await pilot.pause(); await pilot.pause(0.05)
            from textual.widgets import Static
            detail = screen.query_one("#hist-scr-detail-text", Static)
            txt = str(detail.content)
            assert "more" in txt, txt
            # And the per-list cap is respected — Enz25 should not
            # appear (only Enz0..Enz11 shown for cap=12).
            assert "Enz25" not in txt

    async def test_expand_and_collapse_actions_no_error(self, tiny_record):
        """Pressing `e` / `c` while the tree has focus drives the
        expand-all / collapse-all helpers without raising."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        root = self._make_chain(5)
        app = _build_app(tiny_record, isolated_library=None)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            screen = sc.HistoryScreen("toggle", root)
            await app.push_screen(screen)
            await pilot.pause()
            screen.action_collapse_all()
            await pilot.pause()
            screen.action_expand_all()
            await pilot.pause()
            # No assertion needed — the test fails if either call raises.

    async def test_viewport_fits_key_widgets_at_160x48(self, tiny_record):
        """Smoke check that every load-bearing widget — title bar,
        subtitle, tree, detail pane, button row, footer — has a
        non-zero rendered region inside the 160×48 standard terminal
        size. Regression guard against accidentally cropping critical
        UI when restyling."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        root = self._make_chain(3)
        app = _build_app(tiny_record, isolated_library=None)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            screen = sc.HistoryScreen("vp", root)
            await app.push_screen(screen)
            await pilot.pause(); await pilot.pause(0.05)
            for sel in (
                "#hist-scr-title", "#hist-scr-subtitle",
                "#hist-scr-tree", "#hist-scr-detail",
                "#hist-scr-btns", "#btn-hist-scr-close",
            ):
                w = screen.query_one(sel)
                r = w.region
                assert r.width > 0 and r.height > 0, (sel, r)
                # Inside the terminal bounds (top-left at 0,0).
                assert r.x >= 0 and r.y >= 0, (sel, r)
                assert r.x + r.width <= TERMINAL_SIZE[0], (sel, r)
                assert r.y + r.height <= TERMINAL_SIZE[1], (sel, r)


# ──────────────────────────────────────────────────────────────────────────────
# Constructor → history wiring
# ──────────────────────────────────────────────────────────────────────────────


class TestConstructorHistory:
    """Constructor assemblies attach a `history_xml` to their library
    entry, parents (the entry vector + every L0 part) appear as
    nested nodes, and previously-saved parents inherit their full
    subtree so the lineage chains through multi-step builds."""

    def test_persist_assembly_attaches_history_xml(
            self, isolated_library, isolated_parts_bin):
        """`_persist_assembly` writes `history_xml` onto the library
        entry. The root node carries the assembly's name + size; the
        vector and every part hang off the root as parents."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("AAAA" * 100), id="MyTU", name="MyTU")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=0,
            entry_vector={"name": "alpha1_vec",
                            "gb_text": "LOCUS x 1 bp DNA\n//\n"},
            parts=[
                {"name": "P", "oh5": "GGAG", "oh3": "TGAC", "level": 0},
                {"name": "T", "oh5": "TGAC", "oh3": "CGCT", "level": 0},
            ],
            backbone_role="Alpha1",
        )
        entry = sc._load_library()[0]
        assert "history_xml" in entry, entry
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root is not None
        assert root.name == "MyTU.dna"
        assert root.seq_len == len(rec.seq)
        # 3 parents: backbone + 2 parts.
        parents = root.parents
        names = sorted(p.name for p in parents)
        assert names == ["P.dna", "T.dna", "alpha1_vec.dna"], names
        # Input summary records the grammar id so a downstream reader
        # can identify which assembly style produced the plasmid.
        sums = root.input_summaries
        assert sums and "gb_l0" in sums[0]["manipulation"], sums

    def test_persist_assembly_nests_parent_history(
            self, isolated_library, isolated_parts_bin):
        """When a part used in the assembly already has its OWN
        `history_xml` in the library (because it came from an earlier
        Save), the parent node in the new tree carries that full
        subtree — so the lineage chains through L0 → TU → MOD."""
        # Seed library with a "part_with_lineage" that itself has a
        # one-step history (representing an earlier assembly).
        parent_history = sc._serialize_commercialsaas_history(
            sc._CommercialSaaSHistoryNode.new(
                name="part_with_lineage.dna", seq_len=500,
                circular=True, operation="insertFragment", node_id=0,
            )
        )
        sc._save_library([{
            "id":          "lineage_part",
            "name":        "part_with_lineage",
            "size":        500,
            "n_feats":     0,
            "source":      "test",
            "added":       "2026-05-11",
            "gb_text":     "LOCUS x 500 bp DNA\n//\n",
            "history_xml": parent_history,
        }])
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("AAAA" * 100), id="MyTU2", name="MyTU2")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=0,
            entry_vector={"name": "alpha1_vec",
                            "gb_text": "LOCUS x 1 bp DNA\n//\n"},
            parts=[
                {"name": "part_with_lineage", "oh5": "GGAG",
                  "oh3": "CGCT", "level": 0},
            ],
            backbone_role="Alpha1",
        )
        entry = next(e for e in sc._load_library()
                       if e.get("id") == "MyTU2")
        assert "history_xml" in entry
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root is not None
        # Locate the nested parent — it must be the same node that
        # was serialised above, not a freshly synthesised leaf.
        nested = [p for p in root.parents
                    if p.name == "part_with_lineage.dna"]
        assert nested, [p.name for p in root.parents]
        # `_parent_node_for_entry` re-uses parent's existing element,
        # so node attributes round-trip through the parse.
        nested_node = nested[0]
        assert nested_node.seq_len == 500
        assert nested_node.operation == "insertFragment"


# ──────────────────────────────────────────────────────────────────────────────
# Sidecar storage + .dna export (Phase 4d)
# ──────────────────────────────────────────────────────────────────────────────


class TestDnaSidecarStorage:
    def test_save_and_load_round_trip(self):
        data = b"\xfdfake commercialsaas bytes"
        assert sc._save_dna_original("foo", data) is True
        recovered = sc._load_dna_original("foo")
        assert recovered == data

    def test_save_refuses_empty_bytes(self):
        # Don't silently nuke an existing sidecar with empty bytes.
        assert sc._save_dna_original("foo", b"") is False

    def test_save_refuses_oversized(self, monkeypatch):
        monkeypatch.setattr(_fileio, "_DNA_SIDECAR_MAX_BYTES", 8)
        # Cap is 8 bytes; 16-byte payload rejected.
        assert sc._save_dna_original("foo", b"x" * 16) is False
        # Sidecar file should NOT exist.
        assert sc._load_dna_original("foo") is None

    def test_load_returns_none_for_unknown_id(self):
        assert sc._load_dna_original("never_saved") is None

    def test_delete_removes_existing_sidecar(self):
        sc._save_dna_original("victim", b"some bytes")
        assert sc._load_dna_original("victim") is not None
        assert sc._delete_dna_original("victim") is True
        assert sc._load_dna_original("victim") is None
        # Idempotent — second delete is a no-op.
        assert sc._delete_dna_original("victim") is False

    def test_save_handles_path_traversal_attempts(self):
        """A user-controlled `entry_id` containing slashes must NOT
        let a malicious entry write outside the sidecar dir. The
        helper sanitises by replacing path separators."""
        evil = "../../etc/passwd"
        sc._save_dna_original(evil, b"trying to escape")
        # The actual file lives under the sidecar dir, not at the
        # escaped path.
        actual = sc._dna_sidecar_path(evil)
        assert actual.parent == sc._state._DNA_ORIGINALS_DIR

    def test_sidecar_path_rejects_pure_dot_segments(self):
        """Regression guard for 2026-05-06 fix: an `entry_id` of just
        ``..`` or ``.`` is normalised to a sentinel so the resulting
        sidecar can't shadow / collide with the parent dir entry on
        weird filesystems."""
        for evil in ("..", ".", "...", "./..", "/", "\\", "//"):
            p = sc._dna_sidecar_path(evil)
            assert p.parent == sc._state._DNA_ORIGINALS_DIR
            assert p.name not in (".dna", "..dna")  # never a dotfile
            assert "/" not in p.name
            assert "\\" not in p.name

    def test_sidecar_path_rejects_nul_bytes(self):
        """NUL bytes in a path raise on POSIX; normalise to underscore
        instead so the sentinel rule applies cleanly."""
        p = sc._dna_sidecar_path("foo\x00bar")
        assert p.parent == sc._state._DNA_ORIGINALS_DIR
        assert "\x00" not in p.name

    def test_sidecar_path_rejects_absolute_path_id(self):
        """An entry_id that's a fully-qualified path must still produce
        a sidecar inside the originals dir (basename only)."""
        p = sc._dna_sidecar_path("/etc/passwd")
        assert p.parent == sc._state._DNA_ORIGINALS_DIR
        # Underscores replace slashes before basename extraction, so the
        # full sanitised id ends up as the filename.
        assert p.name.endswith(".dna")


class TestExportCommercialSaaSDna:
    def test_export_round_trips_through_history_splice(self, tmp_path):
        # Build a synthetic .dna with one history step.
        original_xml = "<HistoryTree><Node name='orig.dna'/></HistoryTree>"
        original_payload = sc._pack_commercialsaas_history_payload(original_xml)
        sidecar = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x07, original_payload),
            (0x06, b"<Notes/>"),
        )
        sc._save_dna_original("test_entry", sidecar)
        # Build a NEW history XML to splice in.
        new_xml = ('<?xml version="1.0"?>'
                    '<HistoryTree><Node name="modified.dna" type="DNA" '
                    'seqLen="100" circular="1" '
                    'operation="insertFragment"/></HistoryTree>')
        entry = {"id": "test_entry", "name": "test",
                  "history_xml": new_xml}
        out_path = tmp_path / "exported.dna"
        result = sc._export_commercialsaas_dna(entry, out_path)
        assert result == str(out_path.resolve())
        assert out_path.exists()
        # Re-extract and confirm the new history is in.
        recovered = sc._extract_commercialsaas_history_xml(out_path.read_bytes())
        assert recovered is not None
        assert "modified.dna" in recovered
        assert "orig.dna" not in recovered

    def test_export_strips_history_when_field_absent(self, tmp_path):
        original_xml = "<HistoryTree><Node name='orig.dna'/></HistoryTree>"
        sidecar = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x07, sc._pack_commercialsaas_history_payload(original_xml)),
        )
        sc._save_dna_original("entry_no_hist", sidecar)
        # Entry has NO history_xml field — export should remove the
        # 0x07 packet entirely.
        entry = {"id": "entry_no_hist", "name": "x"}
        out_path = tmp_path / "stripped.dna"
        sc._export_commercialsaas_dna(entry, out_path)
        assert sc._extract_commercialsaas_history_xml(out_path.read_bytes()) is None

    def test_export_preserves_unknown_packets(self, tmp_path):
        """Splice-only export must keep every non-history packet
        byte-equal — that's the whole reason we sidecar the
        original instead of regenerating from scratch."""
        weird = b"\x00\xff\xaa\x55weird payload"
        sidecar = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
            (0x42, weird),
        )
        sc._save_dna_original("weirdo", sidecar)
        new_xml = ("<HistoryTree><Node name='x' type='DNA' seqLen='4' "
                    "circular='1' operation='insertFragment'/></HistoryTree>")
        entry = {"id": "weirdo", "name": "weirdo", "history_xml": new_xml}
        out_path = tmp_path / "weirdo.dna"
        sc._export_commercialsaas_dna(entry, out_path)
        out_packets = {t: p for t, _l, p in
                        sc._iter_commercialsaas_packets(out_path.read_bytes())}
        assert out_packets[0x42] == weird

    def test_export_falls_back_to_from_scratch_when_no_sidecar(
            self, tmp_path, tiny_record):
        """Phase 3 wiring: when there's no sidecar, the export
        rebuilds a `.dna` from the entry's GenBank text instead of
        raising. Result must round-trip back through BioPython
        (sequence + features preserved)."""
        # Build an entry from a real SeqRecord; no sidecar saved.
        from pathlib import Path as _P
        entry = sc._record_to_library_entry(tiny_record, _P("tiny.gb"))
        out_path = tmp_path / "from_scratch.dna"
        result = sc._export_commercialsaas_dna(entry, out_path)
        assert result == str(out_path.resolve())
        assert out_path.exists()
        # BioPython must accept it.
        from Bio import SeqIO
        rec2 = SeqIO.read(str(out_path), sc._BIOPYTHON_DNA_FMT)
        assert len(rec2.seq) == len(tiny_record.seq)

    def test_export_raises_on_id_less_entry(self, tmp_path):
        with pytest.raises(ValueError):
            sc._export_commercialsaas_dna({"name": "no-id"},
                                       tmp_path / "out.dna")

    def test_export_raises_when_no_gb_text_and_no_sidecar(self, tmp_path):
        entry = {"id": "ghost", "name": "ghost"}
        with pytest.raises(ValueError) as exc:
            sc._export_commercialsaas_dna(entry, tmp_path / "out.dna")
        assert "gb_text" in str(exc.value).lower()


class TestWriteCommercialSaaSDnaBytes:
    """Exercise the from-scratch writer (Phase 3) directly: build a
    SeqRecord, write `.dna` bytes, round-trip back through both our
    own packet iterator AND BioPython's CommercialSaaS reader. Catches
    schema regressions (a feature attribute we accidentally drop,
    a bad XML entity escape, etc.) without needing CommercialSaaS Viewer
    in the loop."""

    def test_writes_minimum_viable_packets(self, tiny_record):
        data = sc._write_commercialsaas_dna_bytes(tiny_record)
        types = [t for t, _l, _p in sc._iter_commercialsaas_packets(data)]
        # Cookie first, then DNA, features, notes (no primers / no
        # cosmetic packets in v1).
        assert types[:4] == [0x09, 0x00, 0x0A, 0x06]

    def test_writer_round_trips_via_biopython(self, tiny_record, tmp_path):
        """A BioPython re-read of our writer's output must recover
        the sequence (length-equal, content-equal up to case),
        topology, and the same number of features."""
        data = sc._write_commercialsaas_dna_bytes(tiny_record)
        out = tmp_path / "rt.dna"
        out.write_bytes(data)
        from Bio import SeqIO
        rec2 = SeqIO.read(str(out), sc._BIOPYTHON_DNA_FMT)
        assert len(rec2.seq) == len(tiny_record.seq)
        assert str(rec2.seq).upper() == str(tiny_record.seq).upper()
        assert (rec2.annotations.get("topology") ==
                tiny_record.annotations.get("topology"))
        # Same feature count (excluding source).
        rt_feats = [f for f in rec2.features if f.type != "source"]
        src_feats = [f for f in tiny_record.features
                       if f.type != "source"]
        assert len(rt_feats) == len(src_feats)

    def test_writer_includes_history_when_provided(self, tiny_record,
                                                       tmp_path):
        history_xml = ('<?xml version="1.0"?>'
                        '<HistoryTree>'
                        '<Node name="x.dna" type="DNA" seqLen="100" '
                        'circular="1" operation="insertFragment"/>'
                        '</HistoryTree>')
        data = sc._write_commercialsaas_dna_bytes(
            tiny_record, history_xml=history_xml)
        # 0x07 history packet present.
        types = [t for t, _l, _p in sc._iter_commercialsaas_packets(data)]
        assert 0x07 in types
        # Round-trip the history XML.
        recovered = sc._extract_commercialsaas_history_xml(data)
        assert recovered == history_xml

    def test_writer_emits_directionality_for_strand(self, tiny_record):
        """Forward / reverse strand maps to CommercialSaaS's
        `directionality="1"` / `"2"` attribute on `<Feature>`."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("A" * 100), id="P", name="P",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        rec.features = [
            SeqFeature(FeatureLocation(10, 20, strand=1),
                        type="CDS", qualifiers={"label": ["fwd"]}),
            SeqFeature(FeatureLocation(30, 40, strand=-1),
                        type="CDS", qualifiers={"label": ["rev"]}),
            SeqFeature(FeatureLocation(50, 60, strand=0),
                        type="misc_feature",
                        qualifiers={"label": ["unstranded"]}),
        ]
        data = sc._write_commercialsaas_dna_bytes(rec)
        # Pull out the features XML payload.
        for t, _l, payload in sc._iter_commercialsaas_packets(data):
            if t == 0x0A:
                xml_text = payload.decode("utf-8")
                break
        else:
            pytest.fail("features packet not emitted")
        assert 'directionality="1"' in xml_text
        assert 'directionality="2"' in xml_text
        # Unstranded must NOT have a directionality attribute (matches
        # CommercialSaaS's convention of omitting it for non-directional
        # features). Match exactly the `<Feature ...>` open tag that
        # carries `name="unstranded"`; checking the surrounding
        # context too broadly picks up the prior feature's attrs.
        import re as _re
        m = _re.search(r'<Feature [^>]*?name="unstranded"[^>]*?>',
                         xml_text)
        assert m is not None, "couldn't locate unstranded <Feature>"
        assert "directionality" not in m.group(0), m.group(0)

    def test_writer_emits_compoundlocation_as_multiple_segments(self):
        """CompoundLocation features (wrap-around, spliced CDSes)
        get one `<Segment>` per part — matches CommercialSaaS's convention
        for storing multi-part features."""
        from Bio.SeqFeature import (SeqFeature, FeatureLocation,
                                       CompoundLocation)
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("A" * 200), id="P", name="P",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        # Wrap feature 180..200 + 0..20.
        rec.features = [
            SeqFeature(CompoundLocation([
                FeatureLocation(180, 200, strand=1),
                FeatureLocation(0, 20, strand=1),
            ]), type="misc_feature", qualifiers={"label": ["wrap"]}),
        ]
        data = sc._write_commercialsaas_dna_bytes(rec)
        for t, _l, p in sc._iter_commercialsaas_packets(data):
            if t == 0x0A:
                xml_text = p.decode("utf-8")
                break
        # Two Segment elements for the wrap feature.
        assert xml_text.count("<Segment ") == 2
        # Ranges are 1-based.
        assert 'range="181-200"' in xml_text
        assert 'range="1-20"' in xml_text

    def test_writer_rejects_empty_sequence(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq(""), id="empty", name="empty")
        with pytest.raises(ValueError) as exc:
            sc._write_commercialsaas_dna_bytes(rec)
        assert "empty" in str(exc.value).lower()

    def test_writer_rejects_none_record(self):
        with pytest.raises(ValueError):
            sc._write_commercialsaas_dna_bytes(None)

    def test_writer_emits_default_primers_packet(self):
        """The from-scratch writer now emits a default 0x05 Primers
        packet (just ``HybridizationParams`` defaults; no ``<Primer>``
        entries — primer features still ride on the 0x0A features
        packet) so the output's packet inventory matches what the
        commercial editor produces. Pinned byte-for-byte against the
        FFE_* fixtures' default 217-byte payload."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 100), id="t",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        data = sc._write_commercialsaas_dna_bytes(rec)
        primers_payloads = [
            payload for type_byte, _length, payload
            in sc._iter_commercialsaas_packets(data)
            if type_byte == 0x05
        ]
        assert len(primers_payloads) == 1, (
            "writer must emit exactly one 0x05 Primers packet"
        )
        text = primers_payloads[0].decode("utf-8")
        assert text.startswith('<?xml version="1.0"?>')
        assert "<Primers nextValidID=\"0\">" in text
        # HybridizationParams defaults must match what real CommercialSaaS
        # files carry — these aren't arbitrary: they're the editor's
        # save-time defaults and Viewer reads them for primer search.
        assert 'minContinuousMatchLen="10"' in text
        assert 'allowMismatch="1"' in text
        assert 'minMeltingTemperature="40"' in text
        assert 'showAdditionalFivePrimeMatches="1"' in text
        assert 'minimumFivePrimeAnnealing="15"' in text

    def test_writer_emits_default_addprops_packet(self):
        """The from-scratch writer also emits a default 0x08
        AdditionalSequenceProperties packet matching the FFE_* fixtures'
        289-byte default — Upstream/DownstreamStickiness=0 (blunt) and
        FivePrimePhosphorylated end modifications. Real CommercialSaaS
        files emit this even on circular plasmids, so the editor's
        Sequence Properties inspector renders our output without
        falling back to (empty) defaults."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 100), id="t",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        data = sc._write_commercialsaas_dna_bytes(rec)
        addprops_payloads = [
            payload for type_byte, _length, payload
            in sc._iter_commercialsaas_packets(data)
            if type_byte == 0x08
        ]
        assert len(addprops_payloads) == 1
        text = addprops_payloads[0].decode("utf-8")
        assert text.startswith("<AdditionalSequenceProperties>")
        assert "<UpstreamStickiness>0</UpstreamStickiness>" in text
        assert "<DownstreamStickiness>0</DownstreamStickiness>" in text
        assert "FivePrimePhosphorylated" in text


class TestWriterHardening:
    """Phase 5 — adversarial inputs the writer must handle without
    blowing up. CommercialSaaS Viewer compatibility is gated separately
    on real-machine validation; these checks catch bugs that would
    crash the writer or produce malformed XML."""

    def _record(self, seq="A" * 100, features=None, *, circular=True):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq(seq), id="P", name="P",
            annotations={"molecule_type": "DNA",
                          "topology": "circular" if circular else "linear"},
        )
        if features:
            rec.features = features
        return rec

    def test_idempotent_round_trip(self, tiny_record):
        """Write → re-read → re-write → re-read must converge: the
        second cycle's output should match the first's. Catches
        schema drift where some attribute is added or dropped on
        every pass."""
        from Bio import SeqIO
        import io as _io
        first = sc._write_commercialsaas_dna_bytes(tiny_record)
        rec1 = SeqIO.read(_io.BytesIO(first), sc._BIOPYTHON_DNA_FMT)
        second = sc._write_commercialsaas_dna_bytes(rec1)
        rec2 = SeqIO.read(_io.BytesIO(second), sc._BIOPYTHON_DNA_FMT)
        # Sequence + topology + feature count stable across cycles.
        assert str(rec2.seq).upper() == str(rec1.seq).upper()
        assert (rec2.annotations.get("topology")
                == rec1.annotations.get("topology"))
        assert len([f for f in rec2.features if f.type != "source"]) == \
                len([f for f in rec1.features if f.type != "source"])

    def test_xml_escapes_special_chars_in_qualifier_values(self):
        """A qualifier containing `<>&"'` characters must be XML-
        escaped on write — otherwise the resulting XML wouldn't
        parse. ElementTree does this automatically; the test locks
        that behaviour in so a future "manual XML build" refactor
        doesn't regress it."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._record(features=[
            SeqFeature(FeatureLocation(0, 50, strand=1),
                        type="CDS",
                        qualifiers={
                            "label": ["test"],
                            "note":  ['Has "quotes" & <tags>'],
                        }),
        ])
        data = sc._write_commercialsaas_dna_bytes(rec)
        # Pull the features XML and re-parse it via xml.etree to make
        # sure the escaping survives a full read.
        for t, _l, p in sc._iter_commercialsaas_packets(data):
            if t == 0x0A:
                xml_text = p.decode("utf-8")
                break
        import xml.etree.ElementTree as _ET
        root = _ET.fromstring(xml_text)
        # Find the V element with our note.
        for v in root.iter("V"):
            text = v.get("text") or ""
            if "quotes" in text:
                assert text == 'Has "quotes" & <tags>'
                break
        else:
            pytest.fail("note V element with our text not found")

    def test_writer_handles_no_features(self):
        """Empty features list emits a valid `<Features
        nextValidID="0">` — must round-trip."""
        rec = self._record(features=[])
        data = sc._write_commercialsaas_dna_bytes(rec)
        from Bio import SeqIO
        import io as _io
        rt = SeqIO.read(_io.BytesIO(data), sc._BIOPYTHON_DNA_FMT)
        assert len([f for f in rt.features if f.type != "source"]) == 0

    def test_writer_handles_iupac_bases(self):
        """Sequences with IUPAC ambiguity codes (N, R, Y, etc.) pass
        through unchanged — CommercialSaaS tolerates them in the DNA packet."""
        rec = self._record(seq="ATGCNRYWSKMBDHV" * 5)
        data = sc._write_commercialsaas_dna_bytes(rec)
        from Bio import SeqIO
        import io as _io
        rt = SeqIO.read(_io.BytesIO(data), sc._BIOPYTHON_DNA_FMT)
        # BioPython lowercases on read; compare uppercase.
        assert str(rt.seq).upper() == str(rec.seq).upper()

    def test_writer_strips_control_chars_from_feature_name(self):
        """A feature label with control bytes (NUL, DEL) must NOT
        propagate into the XML — those would break the reader and
        could be used for terminal-escape injection."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._record(features=[
            SeqFeature(FeatureLocation(0, 50, strand=1),
                        type="CDS",
                        qualifiers={"label": ["nasty\x00\x07\x1bname"]}),
        ])
        data = sc._write_commercialsaas_dna_bytes(rec)
        # Find the XML, confirm no control bytes.
        for t, _l, p in sc._iter_commercialsaas_packets(data):
            if t == 0x0A:
                xml_text = p.decode("utf-8")
                break
        # `\x00` etc. should NOT be in the XML.
        for ch in "\x00\x07\x1b":
            assert ch not in xml_text, f"control byte {ch!r} leaked"
        # The clean part of the label survives.
        assert "nastyname" in xml_text

    def test_writer_handles_very_long_feature_label(self):
        """Labels are capped at 200 chars to avoid pathological XML
        from a hand-crafted SeqRecord with a multi-MB qualifier."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        long_label = "X" * 5000
        rec = self._record(features=[
            SeqFeature(FeatureLocation(0, 50, strand=1),
                        type="CDS",
                        qualifiers={"label": [long_label]}),
        ])
        data = sc._write_commercialsaas_dna_bytes(rec)
        for t, _l, p in sc._iter_commercialsaas_packets(data):
            if t == 0x0A:
                xml_text = p.decode("utf-8")
                break
        # The cap is 200 chars (per `_commercialsaas_feat_name`).
        # Find the name attribute in the rendered XML and check length.
        import re as _re
        m = _re.search(r'name="(X+)"', xml_text)
        assert m is not None
        assert len(m.group(1)) <= 200

    def test_writer_works_on_linear_topology(self):
        """Linear topology = flag byte 0x00 (vs 0x01 for circular).
        BioPython must round-trip the topology field."""
        rec = self._record(circular=False)
        data = sc._write_commercialsaas_dna_bytes(rec)
        from Bio import SeqIO
        import io as _io
        rt = SeqIO.read(_io.BytesIO(data), sc._BIOPYTHON_DNA_FMT)
        assert rt.annotations.get("topology") == "linear"

    def test_writer_skips_source_features(self):
        """The synthetic `source` feature BioPython adds for GenBank
        records is NOT a real annotation — CommercialSaaS's features
        packet shouldn't include it."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._record(features=[
            SeqFeature(FeatureLocation(0, 100, strand=1), type="source",
                        qualifiers={"organism": ["synthetic"]}),
            SeqFeature(FeatureLocation(10, 50, strand=1), type="CDS",
                        qualifiers={"label": ["real"]}),
        ])
        data = sc._write_commercialsaas_dna_bytes(rec)
        for t, _l, p in sc._iter_commercialsaas_packets(data):
            if t == 0x0A:
                xml_text = p.decode("utf-8")
                break
        assert "source" not in xml_text
        assert "real" in xml_text


class TestExportCommercialSaaSAction:
    """The File → Export as CommercialSaaS (.dna)… action surfaces a clear
    notify when the loaded record has no CommercialSaaS sidecar to splice
    into. End-to-end happy path is exercised via direct
    `_export_commercialsaas_dna` calls in `TestExportCommercialSaaSDna`."""

    async def test_no_record_loaded_notifies(self, isolated_library):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        # Explicitly bypass the auto-preload path so `_current_record`
        # stays None at the moment the action fires.
        app = _build_app(None, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Force the no-record state regardless of whatever the
            # bootstrap may have auto-loaded.
            app._current_record = None
            n_screens_before = len(app.screen_stack)
            app.action_export_commercialsaas()
            await pilot.pause()
            assert len(app.screen_stack) == n_screens_before, (
                "action shouldn't push a modal when no plasmid loaded")

    async def test_no_sidecar_falls_back_to_from_scratch(
            self, tiny_record, isolated_library):
        """Phase 3 wiring (was Phase 4d notify): when an entry has
        no sidecar, the export action now still pushes the modal —
        the export will use `_write_commercialsaas_dna_bytes` to build
        from scratch instead of erroring out."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from pathlib import Path as _P
        # tiny_record is in the library but has NO sidecar.
        entry = sc._record_to_library_entry(tiny_record, _P("t.gb"))
        sc._save_library([entry])
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.action_export_commercialsaas()
            await pilot.pause()
            # Modal IS pushed; from-scratch export will run on submit.
            assert any(isinstance(s, sc.ExportCommercialSaaSModal)
                        for s in app.screen_stack), (
                f"expected modal on stack, got "
                f"{[type(s).__name__ for s in app.screen_stack]}")

    async def test_sidecar_present_pushes_modal(
            self, tiny_record, isolated_library):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        # Seed library + sidecar.
        sidecar_bytes = _make_minimal_dna(
            (0x00, bytes([0x01]) + b"ATGC"),
        )
        sc._save_library([{
            "id":   tiny_record.id,
            "name": "tiny",
            "size": len(tiny_record.seq),
        }])
        sc._save_dna_original(tiny_record.id, sidecar_bytes)
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.action_export_commercialsaas()
            await pilot.pause()
            assert any(isinstance(s, sc.ExportCommercialSaaSModal)
                        for s in app.screen_stack), (
                f"expected ExportCommercialSaaSModal on stack, got "
                f"{[type(s).__name__ for s in app.screen_stack]}")


# ──────────────────────────────────────────────────────────────────────────────
# Phase 5 hardening — corpus round-trip
# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(not _SAMPLE_DIR.exists(),
                    reason="CommercialSaaS sample corpus not available")
class TestWriterCorpusRoundTrip:
    """Confirm the from-scratch writer's output is BioPython-readable
    across a representative sample of real `.dna` files. This catches
    schema regressions (a feature attribute we accidentally drop, an
    XML escape that breaks parsing, etc.) without needing CommercialSaaS
    Viewer in the loop. CommercialSaaS-Viewer compatibility is the next
    validation step; this is the BioPython-level smoke test."""

    @pytest.fixture(scope="class")
    def corpus_files(self):
        samples = sorted(_SAMPLE_DIR.glob("*.dna"))
        # Cap each file at 500 KB to keep the test fast (the full
        # corpus has chromosome dumps; round-trip on those is real
        # work but not what this regression guard targets).
        return [s for s in samples if s.stat().st_size < 500_000]

    def test_corpus_round_trip_via_biopython(self, corpus_files):
        if not corpus_files:
            pytest.skip("no <500KB sample files found")
        from Bio import SeqIO
        import io as _io
        import warnings as _warnings
        _warnings.simplefilter("ignore")
        failures = []
        for path in corpus_files:
            try:
                original = SeqIO.read(str(path), sc._BIOPYTHON_DNA_FMT)
                rewritten = sc._write_commercialsaas_dna_bytes(original)
                recovered = SeqIO.read(_io.BytesIO(rewritten),
                                          sc._BIOPYTHON_DNA_FMT)
                if str(recovered.seq).upper() != str(original.seq).upper():
                    failures.append((path.name, "seq mismatch"))
                elif (recovered.annotations.get("topology")
                       != original.annotations.get("topology")):
                    failures.append(
                        (path.name,
                         f"topology mismatch: "
                         f"{original.annotations.get('topology')} "
                         f"→ {recovered.annotations.get('topology')}"))
            except Exception as exc:
                failures.append((path.name,
                                  f"{type(exc).__name__}: "
                                  f"{str(exc)[:120]}"))
        assert not failures, (
            f"\n{len(failures)} round-trip failure(s):\n" +
            "\n".join(f"  {n}: {e}" for n, e in failures[:10])
        )


@pytest.mark.skipif(not _SAMPLE_WITH_HISTORY.exists(),
                    reason="CommercialSaaS sample corpus not available")
class TestRealCommercialSaaSFiles:
    def test_extracts_history_from_real_file(self):
        """ECXA.dna in the test corpus carries a populated history.
        The XML root must be ``<HistoryTree>`` and contain at least
        one ``<Node>`` element."""
        data = _SAMPLE_WITH_HISTORY.read_bytes()
        xml = sc._extract_commercialsaas_history_xml(data)
        assert xml is not None
        assert "<HistoryTree>" in xml
        assert "<Node" in xml
        # The top-level node names the file (this is CommercialSaaS's
        # convention — re-saving updates it).
        assert "ECXA.dna" in xml or "ECXA" in xml

    def test_returns_none_for_history_less_file(self):
        if not _SAMPLE_WITHOUT_HISTORY.exists():
            pytest.skip("history-less sample file unavailable")
        data = _SAMPLE_WITHOUT_HISTORY.read_bytes()
        # The PCR-product fixture has no construction history.
        assert sc._extract_commercialsaas_history_xml(data) is None

    def test_round_trip_history_on_real_file(self):
        """Extract → re-inject → re-extract must give the same XML
        back. Note: the `.dna` byte stream itself is NOT byte-equal
        after round-trip because xz compression is non-deterministic
        (different presets / dictionaries can encode the same input
        differently). The semantic invariant is the decoded XML."""
        data = _SAMPLE_WITH_HISTORY.read_bytes()
        xml = sc._extract_commercialsaas_history_xml(data)
        assert xml is not None
        out = sc._inject_commercialsaas_history(data, xml)
        recovered = sc._extract_commercialsaas_history_xml(out)
        assert recovered == xml
        # Every non-history packet should also survive the round-trip
        # byte-exact (the splicer only touches 0x07).
        before = [(t, p) for t, _l, p in sc._iter_commercialsaas_packets(data)
                  if t != 0x07]
        after  = [(t, p) for t, _l, p in sc._iter_commercialsaas_packets(out)
                  if t != 0x07]
        assert before == after


# ══════════════════════════════════════════════════════════════════════════════
# `.dna` import augmentation (regression guard for 2026-05-10 user report).
# ══════════════════════════════════════════════════════════════════════════════
#
# `_augment_dna_record_from_packets` recovers info BioPython's
# commercial-SaaS-format parser drops:
#   1. per-feature colours from the 0x0A features packet — without
#      this, every imported .dna feature gets a deterministic-but-
#      unrelated colour from `_FEATURE_PALETTE` rotation.
#   2. primer_seq stamps on every `primer_bind` feature — without
#      this, primers visualise as plain bars instead of through the
#      seq-panel primer machinery (flap, weak-primer arrow, tooltip).
#   3. primer entries for `primers.json` — both standalone 0x05
#      `<Primer>` entries AND primer_bind features get added to the
#      persistent primer library so the user's primer DB mirrors what
#      they had in the source editor.

_FFE_FIXTURE = Path(__file__).resolve().parent / "FFE 1 ENTRY UPD.dna"


def _fixture_available() -> bool:
    return _FFE_FIXTURE.exists()


@pytest.mark.skipif(not _fixture_available(), reason="FFE fixture not present")
class TestDnaImportAugmentation:
    """End-to-end via `load_genbank` against the real FFE_1 fixture.
    Verifies the augmentation actually fires on the commercial SaaS
    parse path and stamps the qualifiers / primer entries we expect."""

    def test_per_feature_colors_recovered(self):
        rec = sc.load_genbank(str(_FFE_FIXTURE))
        non_source = [f for f in rec.features if f.type != "source"]
        # Every non-source feature in the FFE fixture has a Segment
        # colour in the 0x0A XML, so every one should now have the
        # qualifier stamped.
        for f in non_source:
            assert "ApEinfo_revcolor" in f.qualifiers, (
                f"feature {f.qualifiers.get('label', ['?'])[0]} "
                f"({f.type}) missing colour qualifier"
            )
            assert "ApEinfo_fwdcolor" in f.qualifiers
            c = f.qualifiers["ApEinfo_revcolor"][0]
            # Hex CSS colour shape — 3 or 6 nibbles after #.
            assert c.startswith("#")
            assert len(c) in (4, 7)

    def test_primer_bind_features_carry_primer_seq(self):
        rec = sc.load_genbank(str(_FFE_FIXTURE))
        primers = [f for f in rec.features if f.type == "primer_bind"]
        assert len(primers) >= 1
        for p in primers:
            assert "primer_seq" in p.qualifiers, (
                f"primer_bind {p.qualifiers.get('label', ['?'])[0]} "
                f"missing primer_seq stamp"
            )
            ps = p.qualifiers["primer_seq"][0]
            # Sanity: pure ACGT, length matches the bound region.
            assert set(ps) <= set("ACGTN")
            start = int(p.location.start)
            end   = int(p.location.end)
            assert len(ps) == end - start

    def test_reverse_strand_primer_gets_rc_of_bound_region(self):
        """For a `(-)`-strand primer, the stamped `primer_seq` should
        read 5'→3' on the bottom strand — i.e., the RC of the top
        strand bases at the bound region."""
        rec = sc.load_genbank(str(_FFE_FIXTURE))
        seq_str = str(rec.seq).upper()
        for p in rec.features:
            if p.type != "primer_bind":
                continue
            strand = p.location.strand
            if strand != -1:
                continue
            start = int(p.location.start)
            end   = int(p.location.end)
            ps = p.qualifiers["primer_seq"][0]
            assert ps == sc._rc(seq_str[start:end]), (
                f"reverse primer at {start}..{end} primer_seq doesn't "
                f"match RC of bound region"
            )

    def test_primer_db_entries_stashed_on_record(self):
        """The augment helper stashes a list of primer DB dicts on
        `_dna_primer_entries`. `_apply_record` (App-level) is what
        actually flushes them to `primers.json`; here we just verify
        the stash."""
        rec = sc.load_genbank(str(_FFE_FIXTURE))
        entries = getattr(rec, "_dna_primer_entries", None)
        assert isinstance(entries, list)
        assert len(entries) >= 1
        for e in entries:
            assert e["primer_type"] == "imported"
            assert e["source"] == ".dna import"
            assert e["status"] == "Imported"
            assert e["sequence"]
            assert e["name"]

    def test_palette_color_overridden_by_qualifier(self):
        """Regression: pre-fix, `PlasmidMap._parse` unconditionally
        assigned a palette colour, throwing away the user's colours.
        After the fix, the qualifier value wins."""
        rec = sc.load_genbank(str(_FFE_FIXTURE))
        # Pick a feature whose colour we KNOW from inspecting the
        # 0x0A packet: M13 fwd is purple (#a020f0).
        m13 = next(
            (f for f in rec.features
             if f.qualifiers.get("label", [""])[0] == "M13 fwd"),
            None,
        )
        assert m13 is not None
        assert m13.qualifiers["ApEinfo_revcolor"][0].lower() == "#a020f0"

    def test_pre_stamped_primer_seq_still_appends_db_entry(self):
        """Regression for 2026-05-10: pre-fix, `_augment_dna_record_from_packets`
        early-continued when a `primer_bind` feature already carried a
        `primer_seq` qualifier. That `continue` skipped both the
        sequence-derivation AND the primer-DB entry append, so any
        ``.dna`` file round-tripped through splicecraft (or exported
        from ApE / a tool that stamps `primer_seq`) lost its primers
        from the imported DB.

        AB303066.dna is the round-trip case in our fixture set: 2
        `primer_bind` features, both already carrying `primer_seq` from
        the splicecraft writer. Without the fix, `_dna_primer_entries`
        would be empty; with the fix, both primers land in the queue."""
        ab_fixture = Path(__file__).resolve().parent.parent / "AB303066.dna"
        if not ab_fixture.exists():
            pytest.skip("AB303066.dna fixture not present")
        rec = sc.load_genbank(str(ab_fixture))
        n_pb = len([f for f in rec.features if f.type == "primer_bind"])
        entries = getattr(rec, "_dna_primer_entries", [])
        assert n_pb >= 1, "AB303066 should have at least one primer_bind"
        assert len(entries) == n_pb, (
            f"Expected {n_pb} primer DB entries (one per primer_bind), "
            f"got {len(entries)} — the pre-stamped-primer_seq case "
            f"used to skip the append"
        )

    def test_bulk_import_folder_flushes_primers_to_db(
            self, tmp_path, monkeypatch,
    ):
        """Regression for 2026-05-10: pre-fix, `_bulk_import_folder`
        called `load_genbank` on every .dna in the folder (which DID
        run the augment helper and stash `_dna_primer_entries` on each
        SeqRecord), then converted each record to a library entry —
        but the rec was discarded after that and the primer entries
        with it. Bulk-imported primers never reached `primers.json`,
        and the user didn't see them in the Primer Library.

        After the fix, `_bulk_import_folder` itself accumulates the
        primer entries across the batch and writes them to the primer
        DB at the end (dedupe by sequence)."""
        import shutil
        # Build a folder containing both an FFE-style .dna (primer_bind
        # WITHOUT primer_seq qualifier) and the AB303066 case
        # (primer_bind WITH primer_seq) so we cover both code paths
        # the augment helper takes.
        src_folder = tmp_path / "src"
        src_folder.mkdir()
        ffe = Path(__file__).resolve().parent / "FFE 1 ENTRY UPD.dna"
        ab  = Path(__file__).resolve().parent.parent / "AB303066.dna"
        if not (ffe.exists() and ab.exists()):
            pytest.skip(".dna fixtures not present")
        shutil.copy(ffe, src_folder / ffe.name)
        shutil.copy(ab,  src_folder / ab.name)
        # `_protect_user_data` (autouse) already redirected
        # `_PRIMERS_FILE`. Start from an empty primer library so the
        # assertions check the bulk-import path's contribution alone.
        sc._save_primers([])
        entries, failures = sc._bulk_import_folder(src_folder)
        assert failures == []
        assert len(entries) == 2
        primers = sc._load_primers()
        # FFE_1 contributes M13 fwd + M13 rev. AB303066 contributes
        # KanR Promoter-DET-F + KanR Promoter-DET-R. Total 4 unique
        # sequences, none colliding → all 4 land in primers.json.
        seqs = {p["sequence"] for p in primers}
        assert "GTAAAACGACGGCCAGT" in seqs   # M13 fwd
        assert "CAGGAAACAGCTATGAC" in seqs   # M13 rev
        assert "CGTTGTGTCTCAAAATCTCTGATGT" in seqs  # KanR-F (pre-stamped)
        assert "ACACCCCTTGTATTACTGTTTATGT" in seqs  # KanR-R (pre-stamped)
        # Numeric Tm on all four — no `tm=None` legacy crashes.
        for p in primers:
            assert isinstance(p["tm"], (int, float))

    def test_bulk_import_dedupes_primers_across_files(
            self, tmp_path, monkeypatch,
    ):
        """Multiple .dna files in the same folder often share primers
        (e.g. M13 fwd/rev appearing in every pUC-derived plasmid).
        The bulk-import flush must dedupe by sequence so the user
        doesn't end up with 5× M13 fwd entries from a 5-plasmid folder."""
        import shutil
        src_folder = tmp_path / "src"
        src_folder.mkdir()
        # Copy all five FFE fixtures — each has M13 fwd + M13 rev.
        test_dir = Path(__file__).resolve().parent
        names = [
            "FFE 1 ENTRY UPD.dna",
            "FFE 2 ENTRY A1.dna",
            "FFE 3 ENTRY A2.dna",
            "FFE 4 ENTRY O1.dna",
            "FFE 5 ENTRY O2.dna",
        ]
        present = [n for n in names if (test_dir / n).exists()]
        if len(present) < 2:
            pytest.skip("Need ≥2 FFE fixtures for dedupe test")
        for n in present:
            shutil.copy(test_dir / n, src_folder / n)
        sc._save_primers([])
        entries, failures = sc._bulk_import_folder(src_folder)
        assert failures == []
        primers = sc._load_primers()
        # M13 fwd + M13 rev are the only two unique sequences across
        # all FFE fixtures. The dedupe must collapse them to exactly 2.
        seqs = {p["sequence"] for p in primers}
        assert seqs == {"GTAAAACGACGGCCAGT", "CAGGAAACAGCTATGAC"}, (
            f"Expected dedupe to 2 unique primers, got {len(seqs)}: {seqs}"
        )

    def test_pre_stamped_primer_seq_preserved_verbatim(self):
        """When a primer_bind has a pre-existing `primer_seq` qualifier
        (e.g. with a 5' flap, longer than the bound region), the augment
        helper must use THAT sequence verbatim in the DB entry — NOT
        re-derive from the bound region and drop the flap."""
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        seq = "ACGTACGT" * 10  # 80 bp
        # A primer with a 5' flap: bound region is seq[0:20] (20 bp);
        # the full primer (with flap) is 5 bp longer.
        bound = seq[0:20]
        full_primer_with_flap = "GGGGG" + bound  # 25 bp
        prim = SeqFeature(
            location=FeatureLocation(0, 20, strand=1),
            type="primer_bind",
            qualifiers={
                "label": ["test_with_flap"],
                "primer_seq": [full_primer_with_flap],
            },
        )
        rec = SeqRecord(Seq(seq), id="syn", name="syn", features=[prim])
        # Empty .dna byte stream (no 0x0A / 0x05) — we're testing the
        # primer_bind branch only.
        from tests.test_commercialsaas_io import _make_minimal_dna
        data = _make_minimal_dna()
        entries = sc._augment_dna_record_from_packets(rec, data)
        assert len(entries) == 1
        assert entries[0]["sequence"] == full_primer_with_flap, (
            "should preserve the flap-bearing primer_seq verbatim"
        )


class TestAugmentHelperUnit:
    """Direct unit tests on `_augment_dna_record_from_packets` using
    a synthesised .dna byte stream — runs even without fixture files."""

    def _build_dna_with_one_primer(self):
        """Synthesize a minimal .dna with: cookie + DNA seq + 0x0A
        features (one CDS, one primer_bind) + 0x05 with one standalone
        Primer entry."""
        seq = "ATGAAACGCGGGAAATAACCC" * 5  # 105 bp
        # 0x00 DNA packet: 1-byte topology + seq bytes (lowercase is
        # the editor's convention; the BioPython parser tolerates both).
        dna_payload = b"\x01" + seq.encode("ascii")
        # 0x0A features XML — two segments with distinct colours.
        features_xml = (
            '<?xml version="1.0"?>'
            '<Features nextValidID="2">'
            '<Feature recentID="0" name="test_cds" type="CDS" '
            'directionality="1" allowSegmentOverlaps="0" '
            'consecutiveTranslationNumbering="1">'
            '<Segment range="1-21" color="#33ccff" type="standard"/>'
            '</Feature>'
            '<Feature recentID="1" name="test_primer" type="primer_bind" '
            'directionality="1" allowSegmentOverlaps="0" '
            'consecutiveTranslationNumbering="1">'
            '<Segment range="22-42" color="#ff9966" type="standard"/>'
            '</Feature>'
            '</Features>'
        ).encode("utf-8")
        # 0x05 primers XML with one standalone <Primer>.
        primers_xml = (
            '<?xml version="1.0"?>'
            '<Primers nextValidID="1">'
            '<HybridizationParams minContinuousMatchLen="10" '
            'allowMismatch="1" minMeltingTemperature="40" '
            'showAdditionalFivePrimeMatches="1" '
            'minimumFivePrimeAnnealing="15"/>'
            '<Primer name="standalone_X" sequence="GATTACAGATTACA"/>'
            '</Primers>'
        ).encode("utf-8")
        return _make_minimal_dna(
            (0x00, dna_payload),
            (0x0A, features_xml),
            (0x05, primers_xml),
        )

    def _make_synthetic_rec(self, seq_str: str):
        """Build a minimal SeqRecord that mirrors what BioPython would
        produce after parsing our synthetic .dna bytes — one CDS + one
        primer_bind, both forward strand, no colour qualifiers yet."""
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        cds = SeqFeature(
            location=FeatureLocation(0, 21, strand=1),
            type="CDS",
            qualifiers={"label": ["test_cds"]},
        )
        prim = SeqFeature(
            location=FeatureLocation(21, 42, strand=1),
            type="primer_bind",
            qualifiers={"label": ["test_primer"]},
        )
        return SeqRecord(Seq(seq_str), id="syn", name="syn",
                         features=[cds, prim])

    def test_color_stamped_from_synthesized_packet(self):
        data = self._build_dna_with_one_primer()
        rec = self._make_synthetic_rec("ATGAAACGCGGGAAATAACCC" * 5)
        sc._augment_dna_record_from_packets(rec, data)
        # CDS picked up the first <Segment> colour.
        cds = rec.features[0]
        assert cds.qualifiers["ApEinfo_revcolor"][0].lower() == "#33ccff"
        # primer_bind picked up the second.
        prim = rec.features[1]
        assert prim.qualifiers["ApEinfo_revcolor"][0].lower() == "#ff9966"
        # primer_seq stamped from the bound region.
        assert "primer_seq" in prim.qualifiers

    def test_standalone_primer_entry_returned(self):
        """A 0x05 `<Primer name="standalone_X" sequence="GATTACA..."/>`
        must surface as a primer DB entry in the returned list."""
        data = self._build_dna_with_one_primer()
        rec = self._make_synthetic_rec("ATGAAACGCGGGAAATAACCC" * 5)
        extras = sc._augment_dna_record_from_packets(rec, data)
        names = {e["name"] for e in extras}
        seqs  = {e["sequence"] for e in extras}
        assert "standalone_X" in names
        assert "GATTACAGATTACA" in seqs

    def test_dedupe_by_sequence(self):
        """If the same sequence appears in both 0x05 and as a
        primer_bind, the merge step must dedupe."""
        # Forge a 0x05 entry whose sequence happens to match the
        # primer_bind's bound region.
        seq = "ATGAAACGCGGGAAATAACCC" * 5
        # primer_bind is at [21, 42); the bound bases on the top
        # strand are seq[21:42].
        bound = seq[21:42]
        primers_xml = (
            '<?xml version="1.0"?>'
            '<Primers>'
            f'<Primer name="dup" sequence="{bound}"/>'
            '</Primers>'
        ).encode("utf-8")
        features_xml = (
            '<?xml version="1.0"?>'
            '<Features>'
            '<Feature name="test_primer" type="primer_bind" '
            'directionality="1">'
            '<Segment range="22-42" color="#ff9966" type="standard"/>'
            '</Feature>'
            '</Features>'
        ).encode("utf-8")
        data = _make_minimal_dna(
            (0x00, b"\x01" + seq.encode("ascii")),
            (0x0A, features_xml),
            (0x05, primers_xml),
        )
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        prim = SeqFeature(
            location=FeatureLocation(21, 42, strand=1),
            type="primer_bind",
            qualifiers={"label": ["test_primer"]},
        )
        rec = SeqRecord(Seq(seq), id="syn", name="syn", features=[prim])
        extras = sc._augment_dna_record_from_packets(rec, data)
        # Same sequence → one entry only. Standalone (`name="dup"`)
        # comes first in the merge order so it wins the name.
        assert len(extras) == 1
        assert extras[0]["name"] == "dup"

    def test_malformed_color_rejected(self):
        """Defensive: a malformed colour string (no #, wrong length)
        must NOT land in the qualifier — we don't want arbitrary
        strings leaking into our colour-rendering paths."""
        seq = "ATGCCC" * 10
        bad_features_xml = (
            '<?xml version="1.0"?>'
            '<Features>'
            '<Feature name="bad" type="misc_feature">'
            '<Segment range="1-6" color="javascript:alert(1)" '
            'type="standard"/>'
            '</Feature>'
            '</Features>'
        ).encode("utf-8")
        data = _make_minimal_dna(
            (0x00, b"\x01" + seq.encode("ascii")),
            (0x0A, bad_features_xml),
        )
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        f = SeqFeature(
            location=FeatureLocation(0, 6, strand=1),
            type="misc_feature",
            qualifiers={"label": ["bad"]},
        )
        rec = SeqRecord(Seq(seq), id="syn", name="syn", features=[f])
        sc._augment_dna_record_from_packets(rec, data)
        # The malformed colour should NOT be in the qualifier.
        assert "ApEinfo_revcolor" not in rec.features[0].qualifiers


class TestColorQualifierReadInPlasmidMap:
    """`PlasmidMap._parse` must read the ApEinfo colour qualifiers
    before falling back to `_FEATURE_PALETTE`. Applies to all imports
    (not just .dna) — ApE / Geneious .gb files use the same convention."""

    def test_qualifier_overrides_palette(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        f1 = SeqFeature(
            location=FeatureLocation(0, 100, strand=1),
            type="CDS",
            qualifiers={
                "label": ["test"],
                "ApEinfo_revcolor": ["#abcdef"],
                "ApEinfo_fwdcolor": ["#abcdef"],
            },
        )
        rec = SeqRecord(Seq("A" * 500), id="syn", name="syn",
                        features=[f1])
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(rec)
        assert len(feats) == 1
        assert feats[0]["color"].lower() == "#abcdef"

    def test_palette_fallback_when_no_qualifier(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        f1 = SeqFeature(
            location=FeatureLocation(0, 100, strand=1),
            type="CDS",
            qualifiers={"label": ["test"]},
        )
        rec = SeqRecord(Seq("A" * 500), id="syn", name="syn",
                        features=[f1])
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(rec)
        assert len(feats) == 1
        # Must come from the palette — non-empty, not the dummy
        # malformed string.
        assert feats[0]["color"]
        assert feats[0]["color"] in sc._FEATURE_PALETTE


class TestGH17LabelOverride:
    """Regression guard for GH #17 (a user, 2026-05-13): feature
    names containing whitespace landed as `lac\\operator` instead of
    `lac operator` after .dna import. Root cause was BioPython's
    commercial SaaS format parser mangling whitespace on some payloads;
    the fix in v0.8.0 (`_augment_dna_record_from_packets`) re-pins
    `qualifiers["label"]` to the verbatim `<Feature name=...>` XML
    attribute. This test simulates the mangle by handing the augment
    helper a SeqRecord whose label is already corrupted and asserts
    the override restores it from the 0x0A packet.
    """

    def test_label_override_restores_whitespace(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        seq = "A" * 100
        # Simulate the BioPython-mangled output: the label that lands
        # on the parsed SeqRecord has backslashes where the source
        # XML had spaces.
        bad = SeqFeature(
            location=FeatureLocation(0, 12, strand=1),
            type="misc_feature",
            qualifiers={"label": ["lac\\operator"]},
        )
        rec = SeqRecord(Seq(seq), id="syn", name="syn", features=[bad])
        # 0x0A packet carries the clean XML attribute "lac operator".
        features_xml = (
            '<?xml version="1.0"?>'
            '<Features nextValidID="1">'
            '<Feature recentID="0" name="lac operator" type="misc_feature" '
            'directionality="1" allowSegmentOverlaps="0">'
            '<Segment range="1-12" color="#ff0000" type="standard"/>'
            '</Feature>'
            '</Features>'
        ).encode("utf-8")
        data = _make_minimal_dna((0x0A, features_xml))
        sc._augment_dna_record_from_packets(rec, data)
        # Override fires: the XML name wins over the mangled value.
        assert rec.features[0].qualifiers["label"] == ["lac operator"], (
            f"Label override didn't restore whitespace; got "
            f"{rec.features[0].qualifiers.get('label')!r}"
        )

    def test_label_override_preserves_other_printables(self):
        """Confirm the override doesn't accidentally strip spaces,
        slashes, dots, or hyphens — only the control-char set is
        scrubbed."""
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        for label in (
                "Integration Seq",
                "M13 fwd",
                "Lambda T0 Terminator",
                "5'/3' UTR",
                "promoter.minimal",
                "tac-promoter",
        ):
            seq = "A" * 100
            bad = SeqFeature(
                location=FeatureLocation(0, 20, strand=1),
                type="misc_feature",
                qualifiers={"label": ["MANGLED"]},
            )
            rec = SeqRecord(Seq(seq), id="syn", name="syn",
                            features=[bad])
            # Escape XML metachars in the test data so the parser
            # accepts the name attribute verbatim.
            xml_label = (label
                         .replace("&", "&amp;")
                         .replace("<", "&lt;")
                         .replace(">", "&gt;")
                         .replace("'", "&apos;")
                         .replace('"', "&quot;"))
            features_xml = (
                '<?xml version="1.0"?>'
                f'<Features nextValidID="1">'
                f'<Feature recentID="0" name="{xml_label}" type="misc_feature" '
                f'directionality="1" allowSegmentOverlaps="0">'
                f'<Segment range="1-20" color="#ff0000" type="standard"/>'
                f'</Feature>'
                f'</Features>'
            ).encode("utf-8")
            data = _make_minimal_dna((0x0A, features_xml))
            sc._augment_dna_record_from_packets(rec, data)
            assert rec.features[0].qualifiers["label"] == [label], (
                f"label {label!r} not preserved verbatim; got "
                f"{rec.features[0].qualifiers.get('label')!r}"
            )

    def test_label_override_skipped_when_xml_name_empty(self):
        """Some third-party .dna writers omit the `name` attribute on
        Feature elements. In that case the override is skipped and
        whatever BioPython parsed survives — empty XML name must NOT
        clobber a non-empty BioPython label with `[""]`."""
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        seq = "A" * 100
        good = SeqFeature(
            location=FeatureLocation(0, 12, strand=1),
            type="misc_feature",
            qualifiers={"label": ["from_biopython"]},
        )
        rec = SeqRecord(Seq(seq), id="syn", name="syn", features=[good])
        features_xml = (
            '<?xml version="1.0"?>'
            '<Features nextValidID="1">'
            '<Feature recentID="0" type="misc_feature" '
            'directionality="1" allowSegmentOverlaps="0">'
            '<Segment range="1-12" color="#ff0000" type="standard"/>'
            '</Feature>'
            '</Features>'
        ).encode("utf-8")
        data = _make_minimal_dna((0x0A, features_xml))
        sc._augment_dna_record_from_packets(rec, data)
        # BioPython's label survives — override only fires on non-empty xml_name.
        assert rec.features[0].qualifiers["label"] == ["from_biopython"]

    def test_label_override_strips_control_chars(self):
        """Embedded NUL / CR / LF in the XML name would break a
        single-row sidebar render. The override scrubs those before
        pinning the label."""
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        seq = "A" * 100
        bad = SeqFeature(
            location=FeatureLocation(0, 12, strand=1),
            type="misc_feature",
            qualifiers={"label": ["x"]},
        )
        rec = SeqRecord(Seq(seq), id="syn", name="syn", features=[bad])
        # XML `name` attribute carries embedded control chars (encoded
        # as numeric character references because raw bytes are
        # forbidden in XML 1.0). After the parser decodes them,
        # `_CONTROL_CHARS_RE` should strip them.
        features_xml = (
            '<?xml version="1.0"?>'
            '<Features nextValidID="1">'
            '<Feature recentID="0" name="dirty&#10;name" type="misc_feature" '
            'directionality="1" allowSegmentOverlaps="0">'
            '<Segment range="1-12" color="#ff0000" type="standard"/>'
            '</Feature>'
            '</Features>'
        ).encode("utf-8")
        data = _make_minimal_dna((0x0A, features_xml))
        sc._augment_dna_record_from_packets(rec, data)
        # Newline scrubbed; the space and other text survives.
        assert rec.features[0].qualifiers["label"] == ["dirtyname"]


# ══════════════════════════════════════════════════════════════════════════════
# De-noised history viewer — render helpers, protocol summary, dedup tree, and
# the single-file `.dna` open → save → export history round-trip (2026-05-30).
# ══════════════════════════════════════════════════════════════════════════════


def _hist_leaf(name, bp, circ=True):
    return sc._CommercialSaaSHistoryNode.new(
        name=name, seq_len=bp, circular=circ, operation="insertFragment")


def _hist_asm(name, bp, backbone, parts, enzyme):
    r = sc._CommercialSaaSHistoryNode.new(
        name=name, seq_len=bp, circular=True, operation="insertFragment")
    r.add_regenerated_site(enzyme, 0, 1)
    r.add_input_summary(manipulation="gbAssembly",
                          name1=backbone.name, name2="parts")
    r.add_parent(backbone)
    for p in parts:
        r.add_parent(p)
    return r


class TestHistoryRenderHelpers:
    """Pure-unit coverage of the de-noised render toolkit. No event loop —
    these are the helpers both viewers delegate to."""

    def test_clean_name_strips_dna_suffix(self):
        assert sc._history_clean_name("pUC19.dna") == "pUC19"
        assert sc._history_clean_name("pUC19.DNA") == "pUC19"
        assert sc._history_clean_name("plain") == "plain"
        assert sc._history_clean_name("  spaced.dna  ") == "spaced"
        assert sc._history_clean_name("") == "(unnamed)"
        assert sc._history_clean_name(None) == "(unnamed)"  # type: ignore[arg-type]

    def test_size_label_units(self):
        assert sc._history_size_label(712) == "712 bp"
        assert sc._history_size_label(13_900) == "13.9 kb"
        assert sc._history_size_label(2_400_000) == "2.40 Mb"
        assert sc._history_size_label(-5) == "0 bp"      # clamps negatives
        assert sc._history_size_label("nope") == "0 bp"  # type: ignore[arg-type]

    def test_op_label_friendly_and_passthrough(self):
        assert sc._history_op_label("insertFragment") == "assemble"
        assert sc._history_op_label("gibsonAssembly") == "Gibson"
        # Unknown ops pass through verbatim — never assume a closed set.
        assert sc._history_op_label("someFutureOp") == "someFutureOp"
        assert sc._history_op_label("") == ""

    def test_signature_matches_same_node_not_different(self):
        a = _hist_leaf("pENTR.dna", 2604)
        b = _hist_leaf("pENTR.dna", 2604)
        c = _hist_leaf("pENTR.dna", 2999)
        assert sc._history_node_signature(a) == sc._history_node_signature(b)
        assert sc._history_node_signature(a) != sc._history_node_signature(c)
        # `.dna` suffix is normalised away in the signature.
        d = _hist_leaf("pENTR", 2604)
        assert sc._history_node_signature(a) == sc._history_node_signature(d)

    def test_tree_label_omits_op_tag_on_starting_material(self):
        """Raw starting material carries no op tag — but that's because
        real files stamp it `operation="invalid"` (2,116 such nodes in
        the reference library), NOT because it lacks parents.

        The old rule ("tag only nodes with parents") ALSO suppressed
        genuine root actions: an `importFile` / `createDocument` leaf is
        exactly how a molecule entered the project, and it rendered as a
        bare name with no hint of its origin."""
        material = sc._CommercialSaaSHistoryNode.new(
            name="pProm.dna", seq_len=712, circular=True,
            operation="invalid")                     # what real files use
        lbl = sc._history_tree_label(material)
        assert "pProm" in lbl and "712 bp" in lbl
        assert "invalid" not in lbl                  # sentinel never echoed
        # A parent-less node WITH a real root action IS tagged.
        imported = sc._CommercialSaaSHistoryNode.new(
            name="pUC19.dna", seq_len=2686, circular=True,
            operation="importFile")
        assert "import" in sc._history_tree_label(imported)
        asm = _hist_asm("TU.dna", 4000, _hist_leaf("v.dna", 2000),
                         [_hist_leaf("p.dna", 500)], "BsaI")
        assert "assemble" in sc._history_tree_label(asm)  # has parents → tag

    def test_tree_label_flags_linear_only(self):
        assert "linear" not in sc._history_tree_label(_hist_leaf("c", 100, True))
        assert "linear" in sc._history_tree_label(_hist_leaf("l", 100, False))

    def test_tree_label_escapes_markup(self):
        node = _hist_leaf("ev[red]il.dna", 100)
        assert "\\[red]" in sc._history_tree_label(node)  # escaped, not raw

    def test_build_steps_orders_earliest_first_and_dedups(self):
        # Two identical TU subtrees reused in a MOD must collapse to ONE
        # step; steps order earliest (deepest) → final product last.
        def tu():
            return _hist_asm("TU.dna", 4000, _hist_leaf("pENTR.dna", 2000),
                              [_hist_leaf("prom.dna", 500),
                               _hist_leaf("cds.dna", 800)], "Esp3I")
        mod = _hist_asm("MOD.dna", 9000, _hist_leaf("pENTR2.dna", 2500),
                         [tu(), tu()], "BsaI")
        steps = sc._history_build_steps(mod)
        assert [s["product"] for s in steps] == ["TU", "MOD"]
        tu_step = steps[0]
        assert tu_step["enzyme"] == "Esp3I"
        assert tu_step["backbone"] == "pENTR"
        assert set(tu_step["inputs"]) == {"prom", "cds"}

    def test_build_steps_preserves_sibling_order(self):
        mod = _hist_asm("MOD.dna", 9000, _hist_leaf("bb.dna", 2500), [
            _hist_asm("TU_A.dna", 4000, _hist_leaf("v.dna", 2000),
                       [_hist_leaf("a.dna", 500)], "Esp3I"),
            _hist_asm("TU_B.dna", 4000, _hist_leaf("v.dna", 2000),
                       [_hist_leaf("b.dna", 500)], "Esp3I"),
        ], "BsaI")
        assert [s["product"] for s in sc._history_build_steps(mod)] == [
            "TU_A", "TU_B", "MOD"]

    def test_protocol_lines_single_record_placeholder(self):
        """A record with no step nodes still says how it got here.

        A stub history is not "nothing happened" — its root records
        whether the molecule was imported, created de novo, or extracted
        from a selection. The placeholder names that action instead of
        flatly reporting no steps; only a root with NO recorded operation
        (CommercialSaaS's `invalid` sentinel) falls back to the bare form.
        """
        no_op = sc._CommercialSaaSHistoryNode.new(
            name="pUC19.dna", seq_len=2686, circular=True,
            operation="invalid")
        lines = sc._history_protocol_lines(no_op)
        assert len(lines) == 1 and "no construction steps" in lines[0]
        imported = sc._CommercialSaaSHistoryNode.new(
            name="pUC19.dna", seq_len=2686, circular=True,
            operation="importFile")
        lines = sc._history_protocol_lines(imported)
        assert len(lines) == 1 and "import" in lines[0]

    def test_protocol_lines_render_step(self):
        tu = _hist_asm("TU.dna", 4000, _hist_leaf("pENTR.dna", 2000),
                        [_hist_leaf("prom.dna", 500)], "Esp3I")
        joined = sc._history_protocol_lines(tu)[0]
        for token in ("TU", "prom", "pENTR", "Esp3I"):
            assert token in joined

    def test_protocol_lists_all_digest_enzymes(self):
        """A restriction digest uses ≥2 enzymes — the protocol must list
        BOTH, not just the first (user report 2026-06-01: a KpnI + XbaI
        double digest showed only "✂ KpnI")."""
        root = sc._CommercialSaaSHistoryNode.new(
            name="BATCH 55-C.dna", seq_len=5470, circular=True,
            operation="insertFragment", node_id=0)
        root.add_regenerated_site("KpnI", 0, 1)
        root.add_regenerated_site("XbaI", 0, 1)
        root.add_parent(sc._CommercialSaaSHistoryNode.new(
            name="backbone.dna", seq_len=4000, circular=True,
            operation="insertFragment", node_id=0))
        step = sc._history_step_from_node(root)
        assert step["enzymes"] == ["KpnI", "XbaI"]
        assert step["enzyme"] == "KpnI"                  # back-compat
        line = sc._history_protocol_lines(root)[0]
        assert "KpnI" in line and "XbaI" in line         # BOTH shown
        assert line.count("✂") == 1                       # one scissor prefix

    def test_protocol_dedups_and_drops_reverse_sentinel(self):
        root = sc._CommercialSaaSHistoryNode.new(
            name="P.dna", seq_len=100, circular=True,
            operation="insertFragment", node_id=0)
        root.add_regenerated_site("EcoRI", 0, 1)
        root.add_regenerated_site("EcoRI", 0, 1)              # duplicate
        root.add_regenerated_site("(reverse insert)", 0, 1)  # sentinel
        root.add_parent(sc._CommercialSaaSHistoryNode.new(
            name="bb.dna", seq_len=80, circular=True,
            operation="insertFragment", node_id=0))
        step = sc._history_step_from_node(root)
        assert step["enzymes"] == ["EcoRI"]    # deduped, sentinel dropped

    def test_detail_lines_has_properties_and_escapes(self):
        node = sc._CommercialSaaSHistoryNode.new(
            name="ev[red]il.dna", seq_len=1234, circular=False,
            operation="insertFragment")
        text = "\n".join(sc._history_detail_lines(node))
        assert "Properties" in text and "1,234 bp" in text and "linear" in text
        assert "\\[red]" in text   # name markup escaped

    def test_oligos_property_parses_primers(self):
        import xml.etree.ElementTree as ET
        node = sc._CommercialSaaSHistoryNode.new(
            name="PCR.dna", seq_len=520, circular=False,
            operation="amplifyFragment")
        for nm, seq in (("Fwd-1", "ACGT" * 5), ("Rev-1", "TTTT")):
            o = ET.SubElement(node.element, "Oligo")
            o.set("name", nm)
            o.set("sequence", seq)
        olis = node.oligos
        assert [x["name"] for x in olis] == ["Fwd-1", "Rev-1"]
        assert olis[0]["sequence"] == "ACGT" * 5

    def test_detail_amplify_shows_region_and_primers_not_question_marks(self):
        """A SnapGene amplify (PCR) step stores no name1/name2 — the
        detail must show the region + <Oligo> primers (extracted from the
        .dna) instead of a bare "amplify ? ↔ ?" (user report 2026-06-01)."""
        import xml.etree.ElementTree as ET
        from rich.text import Text
        node = sc._CommercialSaaSHistoryNode.new(
            name="PCR-X.dna", seq_len=500, circular=False,
            operation="amplifyFragment")
        node.add_input_summary(manipulation="amplify", val1=100, val2=600)
        o = ET.SubElement(node.element, "Oligo")
        o.set("name", "X-F")
        o.set("sequence", "GCGCGGTACCgccatg")
        node.add_parent(sc._CommercialSaaSHistoryNode.new(
            name="template.dna", seq_len=4000, circular=True,
            operation="insertFragment"))
        plain = "\n".join(Text.from_markup(ln).plain
                          for ln in sc._history_detail_lines(node))
        assert "?" not in plain                       # no bare question marks
        assert "region 100–600" in plain              # amplified region
        assert "Primers" in plain and "X-F" in plain  # primer surfaced
        assert "GCGCGGTACCgccatg" in plain            # primer sequence

    def test_add_oligo_round_trips(self):
        node = sc._CommercialSaaSHistoryNode.new(
            name="x.dna", seq_len=10, circular=False,
            operation="amplifyFragment")
        node.add_oligo(name="F", sequence="ACGT", description="fwd")
        node.add_oligo(name="R", sequence="TTTT")
        assert [(o["name"], o["sequence"]) for o in node.oligos] == [
            ("F", "ACGT"), ("R", "TTTT")]
        assert node.oligos[0]["description"] == "fwd"

    def test_build_amplicon_history_matches_dna_element_set(self):
        """A de-novo amplicon's generated history carries the same
        elements a `.dna` import does — amplifyFragment op, amplify
        region InputSummary, and primer <Oligo>s — so the History detail
        is uniform regardless of origin."""
        from rich.text import Text
        xml = sc._build_amplicon_history_xml(
            name="PCR-X", seq_len=520, fwd_seq="GCGCGGTACCxxx",
            rev_seq="GCGCTCTAGAyyy", start_1based=100, end_1based=620)
        root = sc._parse_commercialsaas_history(xml)
        assert root.operation == "amplifyFragment"
        assert [o["name"] for o in root.oligos] == [
            "forward primer", "reverse primer"]
        assert root.oligos[0]["sequence"] == "GCGCGGTACCxxx"
        sm = root.input_summaries[0]
        assert sm["manipulation"] == "amplify"
        assert sm["val1"] == 100 and sm["val2"] == 620
        plain = "\n".join(Text.from_markup(ln).plain
                          for ln in sc._history_detail_lines(root))
        assert "Primers" in plain and "forward primer" in plain
        assert "region 100–620" in plain

    def test_op_label_maps_sentinels_to_empty(self):
        for s in ("invalid", "Invalid", "unknown", "none", "unspecified", ""):
            assert sc._history_op_label(s) == ""
        assert sc._history_op_label("amplifyFragment") == "PCR"
        # Every op real .dna files carry now has a bench verb — the map
        # was derived from a full scan of the reference library, so
        # `goldenGateAssembly` (the 4th most common manipulation there)
        # no longer reads as raw camelCase mid-protocol.
        assert sc._history_op_label("goldenGateAssembly") == "Golden Gate"
        # A genuinely unknown op still passes through verbatim, never dropped.
        assert sc._history_op_label("someFutureOp2027") == "someFutureOp2027"

    def test_detail_invalid_operation_shows_no_op_not_literal(self):
        """A CommercialSaaS node with operation='invalid' (its sentinel
        for a base/starting sequence with no recorded operation) must
        render '(no operation recorded)', not a literal 'invalid' (user
        report 2026-06-01 on pIT6-NEO)."""
        from rich.text import Text
        node = sc._CommercialSaaSHistoryNode.new(
            name="pIT6.dna", seq_len=8732, circular=True,
            operation="invalid")
        plain = "\n".join(Text.from_markup(ln).plain
                          for ln in sc._history_detail_lines(node))
        assert "invalid" not in plain.lower()
        assert "no operation recorded" in plain

    def test_reopen_nudge_only_for_dna_sources(self):
        # Fires for .dna-sourced entries (both import paths carry .dna)…
        assert sc._history_reopen_nudge("file:DEMO 27.dna")
        assert sc._history_reopen_nudge("/home/x/p.DNA")     # path, any case
        # …never for hand-built constructs or non-.dna imports.
        assert sc._history_reopen_nudge("constructor:gb_l0:backbone") == ""
        assert sc._history_reopen_nudge("file:plasmid.gb") == ""
        assert sc._history_reopen_nudge("") == ""


class TestHistoryViewerDeNoise:
    """Drive the real `HistoryScreen` through Pilot and confirm the
    de-noising fires end-to-end: a protocol pane, and repeated ancestors
    rendered once with later occurrences collapsed to references."""

    async def test_protocol_pane_and_reference_dedup(self, tiny_record):
        from textual.widgets import Tree as _TreeWidget, Static as _Static
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        # Two TUs share the SAME backbone (pENTR, 2000 bp). The shared
        # backbone must render fully once, then as a "shown above" ref.
        mod = _hist_asm("MOD.dna", 9000, _hist_leaf("pENTR2.dna", 2500), [
            _hist_asm("TU1.dna", 4000, _hist_leaf("pENTR.dna", 2000),
                       [_hist_leaf("a.dna", 500)], "Esp3I"),
            _hist_asm("TU2.dna", 4000, _hist_leaf("pENTR.dna", 2000),
                       [_hist_leaf("b.dna", 500)], "Esp3I"),
        ], "BsaI")
        app = _build_app(tiny_record, isolated_library=None)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            screen = sc.HistoryScreen("MOD", mod)
            await app.push_screen(screen)
            await pilot.pause()
            # Protocol pane: 3 build steps (TU1, TU2, MOD). The pane now
            # renders a Rich Group(legend, hanging-indent table); render
            # that SAME renderable to text (the pane shows exactly this).
            assert screen.query_one("#hist-scr-proto-text", _Static) is not None
            from rich.console import Console as _Console
            from io import StringIO as _SIO
            _buf = _SIO()
            _Console(file=_buf, width=120).print(
                sc._history_protocol_renderable(mod))
            proto = _buf.getvalue()
            assert "TU1" in proto and "TU2" in proto and "MOD" in proto
            assert "ingredients" in proto   # symbol legend present
            # One ✂ in the legend key + one per build step (3 steps).
            assert proto.count("✂") == 4
            tree = screen.query_one("#hist-scr-tree", _TreeWidget)
            refs, full_shared = [], []
            def _walk(n):
                if n is not tree.root:
                    lab = str(n.label)
                    if "shown above" in lab:
                        refs.append(n)
                    elif lab.split() and lab.split()[0] == "pENTR":
                        full_shared.append(n)
                for c in n.children:
                    _walk(c)
            _walk(tree.root)
            # Shared backbone drawn fully exactly once; the repeat is a ref.
            assert len(full_shared) == 1, [str(x.label) for x in full_shared]
            assert len(refs) >= 1
            # Reference nodes never redraw their subtree.
            for rn in refs:
                assert len(rn.children) == 0


@pytest.mark.skipif(not _SAMPLE_WITH_HISTORY.exists(),
                     reason="CommercialSaaS history fixture not available")
class TestDnaSingleOpenStash:
    """`load_genbank` on a single `.dna` must stash the construction
    history + original bytes on the record so a later `add_entry`
    round-trips them (pre-fix only the bulk path did this)."""

    def test_load_genbank_stashes_history_and_bytes(self):
        rec = sc.load_genbank(str(_SAMPLE_WITH_HISTORY))
        assert getattr(rec, "_dna_history_xml", None)
        assert "<HistoryTree>" in rec._dna_history_xml
        assert getattr(rec, "_dna_original_bytes", None)


class TestDnaSingleOpenRoundTrip:
    """The fix for the user's report: open one `.dna` → save to library →
    export must preserve the history (and, via the sidecar, every
    CommercialSaaS packet). Exercised at the `add_entry` seam with a
    fabricated record carrying the load-time stash."""

    async def test_single_open_save_export_preserves_history(
            self, tiny_record, isolated_library, tmp_path):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        xml = ("<HistoryTree><Node name='kid.dna' seqLen='80' circular='1' "
               "operation='insertFragment'/></HistoryTree>")
        original = _make_minimal_dna((0x00, bytes([0x01]) + b"ATGCATGCATGC"))
        rec = SeqRecord(Seq("ATGC" * 25), id="rt_open", name="rt_open")
        rec.annotations["topology"] = "circular"
        rec.annotations["molecule_type"] = "DNA"
        rec._dna_history_xml = xml          # stamped by load_genbank's .dna branch
        rec._dna_original_bytes = original
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            assert lib.add_entry(rec) is True
            await pilot.pause()
            # 1. Library entry carries the history.
            entry = next((e for e in sc._load_library()
                          if e.get("id") == "rt_open"), None)
            assert entry is not None and entry.get("history_xml") == xml
            # 2. The original bytes were stashed as a sidecar.
            assert sc._load_dna_original("rt_open") == original
            # 3. Export re-injects the history (splice mode off the sidecar).
            out = tmp_path / "rt_open.dna"
            sc._export_commercialsaas_dna(entry, out)
            recovered = sc._extract_commercialsaas_history_xml(out.read_bytes())
            assert recovered == xml

    async def test_resave_after_edit_keeps_history(
            self, tiny_record, isolated_library):
        """Re-saving an EDITED plasmid (same id, rebuilt record with no
        `_dna_history_xml` attr) must keep the history already on the
        entry — same preserve-on-resave contract as `status`."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        xml = ("<HistoryTree><Node name='h.dna' seqLen='40' circular='1' "
               "operation='insertFragment'/></HistoryTree>")
        sc._save_library([{"id": "edited", "name": "edited", "size": 40,
                           "history_xml": xml}])
        rec = SeqRecord(Seq("ATGC" * 10), id="edited", name="edited")
        rec.annotations["topology"] = "circular"
        rec.annotations["molecule_type"] = "DNA"
        # NOTE: no _dna_history_xml on the record (an edit rebuilt it).
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            assert lib.add_entry(rec) is True
            await pilot.pause()
            entry = next((e for e in sc._load_library()
                          if e.get("id") == "edited"), None)
            assert entry is not None and entry.get("history_xml") == xml


class TestHistoryManipulationDetail:
    """The `<InputSummary>` val1/val2 coordinates — parsed since the
    beginning, rendered nowhere until 2026-08-04. The protocol said WHAT
    combined but never WHERE on the molecule it landed."""

    @staticmethod
    def _node(manip, v1=0, v2=0, name1="", name2=""):
        n = sc._CommercialSaaSHistoryNode.new(
            name="x.dna", seq_len=100, circular=True, operation="replace")
        n.add_input_summary(manipulation=manip, val1=v1, val2=v2,
                             name1=name1, name2=name2)
        return n

    def test_range_manipulations(self):
        assert sc._history_manipulation_detail(
            self._node("amplify", 469, 969)) == "region 469–969"
        assert sc._history_manipulation_detail(
            self._node("remove", 450, 458)) == "bases 450–458"
        assert sc._history_manipulation_detail(
            self._node("editSequence", 7821, 7830)) == "7,821 → 7,830 bp"

    def test_point_manipulations_are_not_rendered_as_ranges(self):
        """`insertAt`/`renumberBase` carry ONE position and leave val2 at
        0 — rendering them as a range produced "renumberBase 610–0".

        The point-vs-range split keys off the manipulation NAME, not off
        which value happens to be zero, because a genuine range can
        legitimately start at base 0 (`amplify 0–542` is a real region)."""
        assert sc._history_manipulation_detail(
            self._node("renumberBase", 610, 0)) == "origin → 610"
        assert sc._history_manipulation_detail(
            self._node("insertAt", 4707, 0)) == "at 4,707"
        # …and a range starting at 0 stays a range.
        assert sc._history_manipulation_detail(
            self._node("amplify", 0, 542)) == "region 0–542"

    def test_named_summary_yields_no_coordinates(self):
        """An assembly's val1/val2 are cut positions on two DIFFERENT
        molecules; a range across them would be meaningless. Those
        summaries name their enzymes instead, and the ✂ tag covers it."""
        assert sc._history_manipulation_detail(
            self._node("insert", 68, 253, name1="BsaI", name2="BsaI")) == ""

    def test_no_coordinates_at_all(self):
        assert sc._history_manipulation_detail(self._node("assembly")) == ""

    def test_unknown_manipulation_still_shows_its_numbers(self):
        """A manipulation this codebase has never seen must still surface
        its coordinates rather than swallow them."""
        out = sc._history_manipulation_detail(self._node("warpDrive", 12, 34))
        assert "12" in out and "34" in out and "warpDrive" in out

    def test_in_place_step_names_what_was_done(self):
        """A single-input step whose input IS the product used to collapse
        to a flat "· edited", so consecutive set-origin / delete / insert
        steps all read identically."""
        from rich.text import Text
        node = sc._CommercialSaaSHistoryNode.new(
            name="P.dna", seq_len=6600, circular=True, operation="setOrigin")
        node.add_input_summary(manipulation="renumberBase", val1=6491)
        node.add_parent(sc._CommercialSaaSHistoryNode.new(
            name="P.dna", seq_len=6600, circular=True, operation="invalid"))
        cells = sc._history_protocol_step_cells(node)
        text = " ".join(Text.from_markup(c).plain for _n, c in cells)
        assert "set origin" in text and "origin → 6,491" in text
        assert "edited" not in text

    def test_inputs_block_and_tag_share_one_formatter(self):
        """The detail pane's Inputs block used to hardcode
        "(region v1–v2)", so a single-position manipulation rendered as
        "renumberBase (region 6,491–0)" there while the tree row said
        "origin → 6,491". Both now go through `_history_coord_phrase`."""
        from rich.text import Text
        node = sc._CommercialSaaSHistoryNode.new(
            name="P.dna", seq_len=6600, circular=True, operation="setOrigin")
        node.add_input_summary(manipulation="renumberBase", val1=6491)
        plain = "\n".join(Text.from_markup(ln).plain
                          for ln in sc._history_detail_lines(node))
        assert "origin → 6,491" in plain
        assert "6,491–0" not in plain

    def test_former_name_ignores_locus_vs_display_form(self):
        """`Demo_12_TU_x` (LOCUS-safe id) and `Demo 12 TU x` (display
        name) are the SAME plasmid — a caller passing the id must not
        make the viewer claim the plasmid was renamed."""
        root = sc._CommercialSaaSHistoryNode.new(
            name="Demo 12 TU Pdemo-GeneX-Tdemo.dna", seq_len=6560,
            circular=True, operation="setOrigin")
        assert sc._history_former_name(
            root, "Demo_12_TU_Pdemo-GeneX-Tdemo") == ""
        assert sc._history_former_name(
            root, "Demo 12 TU Pdemo-GeneX-Tdemo") == ""
        # A genuine rename still reports.
        assert sc._history_former_name(root, "pDemo99") == \
            "Demo 12 TU Pdemo-GeneX-Tdemo"


class TestHistoryConsistencyCheck:
    """`_history_node_warnings` reports claims the DNA contradicts.

    STRICTLY READ-ONLY by contract: it never edits a sequence or a history
    record to make a warning go away. Both are the user's data — the
    history is a faithful account of what a tool did, and silently
    "repairing" either destroys the evidence that they diverged.

    Calibrated against a 1,184-plasmid reference library, where it must
    report ZERO — every apparent error there turned out to be an artefact
    of checking the wrong node (see `test_ancestor_claims_are_never_checked`).
    """

    ENZ = {"KpnI": ("GGTACC", 1, 5), "SacI": ("GAGCTC", 5, 1)}

    def _node(self, **kw):
        kw.setdefault("name", "P.dna")
        kw.setdefault("seq_len", 40)
        kw.setdefault("circular", True)
        kw.setdefault("operation", "insertRestrictionSite")
        return sc._CommercialSaaSHistoryNode.new(**kw)

    def test_absent_enzyme_site_is_reported(self):
        n = self._node()
        n.add_regenerated_site("KpnI", 20, 1)
        w = sc._history_node_warnings(n, "A" * 40, enzymes=self.ENZ)
        assert len(w) == 1 and "KpnI" in w[0] and "does not occur" in w[0]

    def test_present_enzyme_site_is_silent(self):
        n = self._node()
        n.add_regenerated_site("KpnI", 20, 1)
        assert sc._history_node_warnings(
            n, "A" * 17 + "GGTACC" + "A" * 17, enzymes=self.ENZ) == []

    def test_our_own_enzyme_marker_never_warns(self):
        """SpliceCraft records the enzyme an assembly USED as a
        `RegeneratedSite` at pos 0 ([INV-93]) — not a claim that the site
        survives. Type IIS assembly consumes its sites, so checking those
        would fire on 126 correct records in the reference library."""
        n = self._node(operation="insertFragment")
        n.add_regenerated_site("KpnI", 0, 1)
        assert sc._history_node_warnings(n, "A" * 40, enzymes=self.ENZ) == []

    def test_ancestor_claims_are_never_checked(self):
        """THE invariant that keeps false positives at zero.

        A real case from the reference library: a plasmid's lineage
        inserted a KpnI site, then two same-length `replace` steps swapped
        it for SacI. The final molecule has no KpnI — correctly. Checking
        the ANCESTOR's KpnI claim against the CURRENT sequence reports an
        error that does not exist, because an ancestor node describes an
        EARLIER molecule.

        The checker only ever receives a sequence for a node the caller
        can identify with certainty; callers pass "" for everything else,
        and no sequence means no sequence-based check.
        """
        ancestor = self._node()
        ancestor.add_regenerated_site("KpnI", 20, 1)
        final_molecule = "A" * 34 + "GAGCTC"          # SacI, no KpnI
        # With a sequence it would (wrongly) fire …
        assert sc._history_node_warnings(
            ancestor, final_molecule, enzymes=self.ENZ)
        # … which is exactly why an ancestor is passed no sequence.
        assert sc._history_node_warnings(ancestor, "", enzymes=self.ENZ) == []

    def test_primer_that_moved_is_distinguished_from_one_that_vanished(self):
        n = self._node(operation="amplifyFragment")
        import xml.etree.ElementTree as ET
        block = ET.SubElement(n.element, "Primers")
        pr = ET.SubElement(block, "Primer")
        pr.set("name", "F1")
        bs = ET.SubElement(pr, "BindingSite")
        bs.set("location", "5-14")
        bs.set("boundStrand", "0")
        bs.set("annealedBases", "GGGGCCCCTT")
        # present, but at 20 rather than 5 -> stale coordinate, still binds
        moved = "A" * 20 + "GGGGCCCCTT" + "A" * 10
        w = sc._history_node_warnings(moved and n, moved, enzymes=self.ENZ)
        assert len(w) == 1 and "stale" in w[0] and "still binds" in w[0]
        # absent entirely -> the sequence changed
        w = sc._history_node_warnings(n, "T" * 40, enzymes=self.ENZ)
        assert len(w) == 1 and "does not bind" in w[0]
        # exactly where claimed -> silent
        ok = "A" * 5 + "GGGGCCCCTT" + "A" * 25
        assert sc._history_node_warnings(ok and n, ok, enzymes=self.ENZ) == []

    def test_length_arithmetic_needs_no_sequence(self):
        """These checks are safe on EVERY node — they compare only the
        lengths the record itself states."""
        n = self._node(seq_len=100, operation="setOrigin")
        n.add_parent(self._node(name="par.dna", seq_len=120,
                                 operation="invalid"))
        w = sc._history_node_warnings(n, "", enzymes=self.ENZ)
        assert len(w) == 1 and "should not change length" in w[0]

    def test_deletion_arithmetic_uses_inclusive_spans(self):
        """0-based INCLUSIVE — calibrated against the reference library,
        where 193/193 deletions balance on this rule and none do on the
        exclusive reading."""
        n = self._node(seq_len=94, operation="remove")
        n.add_input_summary(manipulation="remove", val1=10, val2=15)  # 6 bases
        n.add_parent(self._node(name="par.dna", seq_len=100,
                                 operation="invalid"))
        assert sc._history_node_warnings(n, "", enzymes=self.ENZ) == []
        bad = self._node(seq_len=90, operation="remove")
        bad.add_input_summary(manipulation="remove", val1=10, val2=15)
        bad.add_parent(self._node(name="par.dna", seq_len=100,
                                   operation="invalid"))
        assert sc._history_node_warnings(bad, "", enzymes=self.ENZ)

    def test_summary_is_empty_when_clean(self):
        assert sc._history_consistency_summary([]) == ""
        assert "1" in sc._history_consistency_summary(["x"])
        assert "2" in sc._history_consistency_summary(["x", "y"])


class TestHistoryHardening:
    """Hostile / degenerate `.dna` input must never crash a viewer or hang
    a render. Every case here was found by an edge-case sweep of the
    consistency checker before release."""

    ENZ = {"KpnI": ("GGTACC", 1, 5)}

    @staticmethod
    def _with_site(loc, annealed, strand=0, seq_len=40):
        import xml.etree.ElementTree as ET
        nd = sc._CommercialSaaSHistoryNode.new(
            name="p.dna", seq_len=seq_len, circular=True,
            operation="amplifyFragment")
        block = ET.SubElement(nd.element, "Primers")
        pr = ET.SubElement(block, "Primer")
        pr.set("name", "F")
        bs = ET.SubElement(pr, "BindingSite")
        bs.set("location", loc)
        bs.set("boundStrand", str(strand))
        bs.set("annealedBases", annealed)
        return nd

    def test_manipulation_name_cannot_inject_a_format_field(self):
        """`manipulation` comes straight off attacker-controllable XML. A
        value like `{evil}` folded into a format string raises KeyError —
        which would take out the tree label AND the protocol pane."""
        n = sc._CommercialSaaSHistoryNode.new(
            name="x.dna", seq_len=10, circular=True, operation="remove")
        n.add_input_summary(manipulation="{evil}", val1=1, val2=5)
        assert sc._history_manipulation_detail(n) == "{evil} 1–5"
        n2 = sc._CommercialSaaSHistoryNode.new(
            name="y.dna", seq_len=10, circular=True, operation="remove")
        n2.add_input_summary(manipulation="{0}{1}%s%%", val1=7, val2=9)
        assert "7" in sc._history_manipulation_detail(n2)

    def test_overlong_manipulation_name_is_capped(self):
        n = sc._CommercialSaaSHistoryNode.new(
            name="x.dna", seq_len=10, circular=True, operation="remove")
        n.add_input_summary(manipulation="Z" * 5000, val1=1, val2=5)
        assert len(sc._history_manipulation_detail(n)) < 80

    def test_negative_binding_offset_never_indexes_from_the_end(self):
        """A negative start would slice from the END of the sequence and
        compare the WRONG bases — a false pass or a false warning, both
        worse than skipping."""
        seq = "ACGT" * 10
        assert sc._history_node_warnings(
            self._with_site("-5-6", "ACGT"), seq, enzymes=self.ENZ) == []

    def test_absurd_or_unparseable_coordinates_are_skipped(self):
        seq = "ACGT" * 10
        for loc in ("999999-1000000", "x-y", "--", "1e9-2e9"):
            assert sc._history_node_warnings(
                self._with_site(loc, "ACGT"), seq, enzymes=self.ENZ) == []

    def test_annealed_run_longer_than_the_molecule_is_skipped(self):
        seq = "ACGT" * 10
        assert sc._history_node_warnings(
            self._with_site("0-9", "A" * 500), seq, enzymes=self.ENZ) == []

    def test_empty_recognition_site_is_not_treated_as_always_present(self):
        """An empty site compiles to a pattern matching at every offset —
        i.e. the check would silently stop checking."""
        n = sc._CommercialSaaSHistoryNode.new(
            name="q.dna", seq_len=40, circular=True,
            operation="insertRestrictionSite")
        n.add_regenerated_site("BOGUS", 5, 1)
        assert sc._history_node_warnings(
            n, "ACGT" * 10, enzymes={"BOGUS": ("", 0, 0)}) == []

    def test_output_and_work_are_both_bounded(self):
        """The warning cap bounds OUTPUT; a node whose sites all match
        yields no warnings and still costs a scan each, so the number of
        claims inspected is capped separately."""
        n = sc._CommercialSaaSHistoryNode.new(
            name="b.dna", seq_len=40, circular=True,
            operation="insertRestrictionSite")
        for i in range(500):
            n.add_regenerated_site("KpnI", i + 1, 1)
        w = sc._history_node_warnings(n, "T" * 40, enzymes=self.ENZ)
        assert len(w) == sc._HISTORY_WARN_MAX

    def test_chromosome_scale_falls_back_to_arithmetic_only(self):
        """Each scan is O(len(seq)); on a 4.5 Mb chromosome that stalls
        `compose()`. Above the cap only the recorded-length checks run."""
        big = "ACGTTGCA" * (sc._HISTORY_CHECK_MAX_BP // 8 + 100)
        n = self._with_site("10-29", "TTTTTTTTTT", seq_len=len(big))
        assert sc._history_node_warnings(n, big, enzymes=self.ENZ) == []
        # …and the same claim IS checked on a normal-sized molecule.
        small = "ACGTTGCA" * 10
        n2 = self._with_site("10-29", "TTTTTTTTTT", seq_len=len(small))
        assert sc._history_node_warnings(n2, small, enzymes=self.ENZ)

    def test_never_raises_on_a_malformed_node(self):
        """It runs inside `compose()` — a bad value must degrade to "no
        warnings", never to a History screen that won't render."""
        class _Hostile:
            name = "x"
            operation = "remove"
            seq_len = "not-an-int"
            circular = True
            parents = property(lambda self: (_ for _ in ()).throw(
                RuntimeError("boom")))
            regenerated_sites = []
            primer_details = []
            input_summaries = []
        assert sc._history_node_warnings(
            _Hostile(), "ACGT", enzymes=self.ENZ) == []


# ══════════════════════════════════════════════════════════════════════════════
# `.dna` label/colour stamping is keyed on COORDINATES, not list position.
# Regression guard for the 2026-08-05 user report. [INV-179]
# ══════════════════════════════════════════════════════════════════════════════
#
# `_augment_dna_record_from_packets` used to zip the 0x0A `<Feature>` list
# against `rec.features` by enumeration index. Any difference in length or
# order silently slid every label onto its neighbour. The trigger in the wild:
# a `.dna` authored by importing a GenBank record carries an explicit
# `<Feature type="source">` as entry #0, but the stamping loop skips the
# record's `source` row — so the XML source entry displaced every subsequent
# name by one. A user's 4-arm homology-arm reference came out with each arm
# wearing its neighbour's name (HomArm2's DNA labelled "HomArm1", …) while
# the sequence itself was byte-perfect: wrong in the way that reads as right.
#
# These tests drive the REAL matching code — the record passed in is the
# input, never a stand-in for the mechanism being verified.


class TestSegmentSpanDecoding:
    """`_commercialsaas_segment_span` — the 1-based-inclusive → BioPython-frame
    decode that the whole coordinate match rests on."""

    def _segs(self, *ranges):
        import xml.etree.ElementTree as ET
        return [ET.fromstring(f'<Segment range="{r}"/>') for r in ranges]

    def test_single_range_decodes_to_biopython_frame(self):
        # The writer emits f"{int(part.start) + 1}-{int(part.end)}", so a
        # feature at BioPython (13, 2239) round-trips through "14-2239".
        assert _fileio._commercialsaas_segment_span(
            self._segs("14-2239")) == (13, 2239)

    def test_multi_segment_join_reduces_to_min_start_max_end(self):
        """A join / origin-spanning feature has several <Segment>s; BioPython
        reports the CompoundLocation's min-start / max-end, so we must too."""
        assert _fileio._commercialsaas_segment_span(
            self._segs("7000-7825", "1-500")) == (0, 7825)

    def test_malformed_and_inverted_ranges_are_skipped(self):
        assert _fileio._commercialsaas_segment_span(self._segs("not-a-range")) is None
        assert _fileio._commercialsaas_segment_span(self._segs("500-100")) is None
        # A 0 or negative start is malformed, not "a feature at coordinate 0".
        assert _fileio._commercialsaas_segment_span(self._segs("0-100")) is None
        # …and one good range among bad ones still lands.
        assert _fileio._commercialsaas_segment_span(
            self._segs("500-100", "14-2239")) == (13, 2239)

    def test_no_segments_is_unmatchable_not_a_crash(self):
        assert _fileio._commercialsaas_segment_span([]) is None


class TestDnaLabelsFollowCoordinates:
    """The stamping loop must hand every feature ITS OWN name."""

    @staticmethod
    def _features_packet(*entries: tuple[str, str, str]) -> bytes:
        """entries: (name, type, "start-end") 1-based inclusive."""
        body = "".join(
            f'<Feature recentID="{i}" name="{name}" type="{ftype}">'
            f'<Segment range="{rng}" color="#33ccff" type="standard"/>'
            f'</Feature>'
            for i, (name, ftype, rng) in enumerate(entries)
        )
        return (f'<?xml version="1.0"?><Features nextValidID="{len(entries)}">'
                f'{body}</Features>').encode("utf-8")

    @staticmethod
    def _rec(*feats: tuple[str, str, int, int]):
        """feats: (label, type, start0, end) in BioPython frame."""
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        return SeqRecord(
            Seq("ACGT" * 2500), id="syn", name="syn",
            features=[SeqFeature(location=FeatureLocation(s, e, strand=1),
                                 type=t, qualifiers={"label": [lab]})
                      for lab, t, s, e in feats])

    @staticmethod
    def _labels(rec):
        return [(f.type, int(f.location.start), int(f.location.end),
                 (f.qualifiers.get("label") or [""])[0])
                for f in rec.features]

    def test_source_entry_in_packet_does_not_shift_every_label(self):
        """THE bug. The 0x0A packet leads with a `source` <Feature>, which the
        old index-zip let displace every following name by one — the tell-tale
        symptom being the first real feature coming out labelled "source"."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("source",    "source",       "1-10000"),
            ("HomArm1",  "misc_feature", "3096-3595"),
            ("HomArm2",  "misc_feature", "3596-4095"),
            ("HomArm3",  "misc_feature", "4096-4595"),
            ("HomArm4",  "misc_feature", "4596-5095"),
        )))
        rec = self._rec(
            ("parsed_source", "source",       0,    10000),
            ("x",             "misc_feature", 3095, 3595),
            ("x",             "misc_feature", 3595, 4095),
            ("x",             "misc_feature", 4095, 4595),
            ("x",             "misc_feature", 4595, 5095),
        )
        sc._augment_dna_record_from_packets(rec, data)
        assert self._labels(rec) == [
            ("source",       0,    10000, "parsed_source"),   # untouched
            ("misc_feature", 3095, 3595,  "HomArm1"),
            ("misc_feature", 3595, 4095,  "HomArm2"),
            ("misc_feature", 4095, 4595,  "HomArm3"),
            ("misc_feature", 4595, 5095,  "HomArm4"),
        ]

    def test_gene_and_cds_at_identical_coords_keep_distinct_names(self):
        """NCBI-derived records carry a gene/CDS pair at the same span; the
        type disambiguates which XML entry belongs to which."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("QFC41_RS02110",      "gene", "2257-4140"),
            ("Beta-galactosidase", "CDS",  "2257-4140"),
        )))
        rec = self._rec(("x", "CDS",  2256, 4140),
                        ("x", "gene", 2256, 4140))
        sc._augment_dna_record_from_packets(rec, data)
        assert self._labels(rec) == [
            ("CDS",  2256, 4140, "Beta-galactosidase"),
            ("gene", 2256, 4140, "QFC41_RS02110"),
        ]

    def test_feature_with_no_packet_entry_keeps_its_parsed_label(self):
        """`primer_bind` features come from the 0x05 packet and have no 0x0A
        entry. They must keep the name BioPython parsed rather than be handed
        whatever entry happens to sit at their coordinates."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("RAND-2A", "misc_feature", "1-25"),
        )))
        rec = self._rec(("x",          "misc_feature", 0, 25),
                        ("RAND-2A-F",  "primer_bind",  0, 25))
        sc._augment_dna_record_from_packets(rec, data)
        assert self._labels(rec) == [
            ("misc_feature", 0, 25, "RAND-2A"),
            ("primer_bind",  0, 25, "RAND-2A-F"),   # NOT renamed to RAND-2A
        ]

    def test_packet_order_independent_of_record_order(self):
        """Nothing pins a third-party writer's `<Feature>` order to
        BioPython's emission order, so the match must not depend on it."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("last",   "misc_feature", "900-1000"),
            ("first",  "misc_feature", "1-100"),
            ("middle", "misc_feature", "400-500"),
        )))
        rec = self._rec(("x", "misc_feature", 0,   100),
                        ("x", "misc_feature", 399, 500),
                        ("x", "misc_feature", 899, 1000))
        sc._augment_dna_record_from_packets(rec, data)
        assert [t[3] for t in self._labels(rec)] == ["first", "middle", "last"]

    def test_colour_follows_the_same_match_as_the_label(self):
        """Colours were stamped off the same index, so they shifted too."""
        packet = (
            '<?xml version="1.0"?><Features nextValidID="3">'
            '<Feature recentID="0" name="source" type="source">'
            '<Segment range="1-10000" color="#ffffff" type="standard"/></Feature>'
            '<Feature recentID="1" name="A" type="misc_feature">'
            '<Segment range="101-200" color="#112233" type="standard"/></Feature>'
            '<Feature recentID="2" name="B" type="misc_feature">'
            '<Segment range="301-400" color="#445566" type="standard"/></Feature>'
            '</Features>').encode("utf-8")
        rec = self._rec(("s", "source",       0,   10000),
                        ("x", "misc_feature", 100, 200),
                        ("x", "misc_feature", 300, 400))
        sc._augment_dna_record_from_packets(rec, _make_minimal_dna((0x0A, packet)))
        assert rec.features[1].qualifiers["ApEinfo_revcolor"] == ["#112233"]
        assert rec.features[2].qualifiers["ApEinfo_revcolor"] == ["#445566"]
        # The source row is skipped by design, so it never takes a colour.
        assert "ApEinfo_revcolor" not in rec.features[0].qualifiers

    def test_empty_xml_name_leaves_the_parsed_label_alone(self):
        """Some third-party writers omit the `name` attribute."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("", "misc_feature", "1-100"),
        )))
        rec = self._rec(("biopython_parsed", "misc_feature", 0, 100))
        sc._augment_dna_record_from_packets(rec, data)
        assert self._labels(rec) == [("misc_feature", 0, 100, "biopython_parsed")]

    def test_unmatched_packet_entry_is_logged(self, caplog):
        """An entry nobody claims means the file's coordinates and BioPython's
        disagree — exactly what the index-zip papered over. It must surface in
        the log instead of silently mislabelling."""
        import logging
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("ghost", "misc_feature", "9000-9500"),
        )))
        rec = self._rec(("kept", "misc_feature", 0, 100))
        # `sc._log` sets propagate=False, so caplog's root handler never sees
        # the record — attach it to that logger for the duration instead.
        logger = logging.getLogger(sc._log.name)
        prev_level = logger.level
        logger.addHandler(caplog.handler)
        logger.setLevel(logging.WARNING)
        try:
            sc._augment_dna_record_from_packets(rec, data)
        finally:
            logger.removeHandler(caplog.handler)
            logger.setLevel(prev_level)
        assert "did not match a parsed feature by coordinate" in caplog.text
        assert "ghost" in caplog.text
        # …and the unmatched record feature kept its own label.
        assert self._labels(rec) == [("misc_feature", 0, 100, "kept")]

    def test_same_span_and_type_are_told_apart_by_strand(self):
        """Two features can share coordinates AND type and differ only in
        direction — the last way one could still be handed the other's name."""
        packet = ('<?xml version="1.0"?><Features nextValidID="2">'
                  '<Feature recentID="0" name="fwd-one" type="promoter" '
                  'directionality="1">'
                  '<Segment range="101-200" color="#111111" type="standard"/>'
                  '</Feature>'
                  '<Feature recentID="1" name="rev-one" type="promoter" '
                  'directionality="2">'
                  '<Segment range="101-200" color="#222222" type="standard"/>'
                  '</Feature></Features>').encode("utf-8")
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        rec = SeqRecord(
            Seq("ACGT" * 100), id="syn", name="syn",
            features=[  # reverse first — the OPPOSITE of the packet order
                SeqFeature(location=FeatureLocation(100, 200, strand=-1),
                           type="promoter", qualifiers={"label": ["x"]}),
                SeqFeature(location=FeatureLocation(100, 200, strand=1),
                           type="promoter", qualifiers={"label": ["x"]}),
            ])
        sc._augment_dna_record_from_packets(rec, _make_minimal_dna((0x0A, packet)))
        assert [f.qualifiers["label"][0] for f in rec.features] == \
            ["rev-one", "fwd-one"]

    def test_a_directionless_packet_entry_still_matches(self):
        """The writer omits `directionality` when it has no opinion; the strand
        tier must degrade to the parts/type match rather than drop the name."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("no-direction", "misc_feature", "101-200"),
        )))
        rec = self._rec(("x", "misc_feature", 100, 200))
        sc._augment_dna_record_from_packets(rec, data)
        assert self._labels(rec) == [("misc_feature", 100, 200, "no-direction")]

    def test_two_wrap_features_are_told_apart_by_their_parts(self):
        """Both origin-spanning features reduce to min-start/max-end = the
        whole molecule, so a span-only key hands them each other's name. The
        exact segment set keeps them straight."""
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        from Bio.Seq import Seq
        packet = ('<?xml version="1.0"?><Features nextValidID="2">'
                  '<Feature recentID="0" name="wrapA" type="misc_feature">'
                  '<Segment range="9900-10000" color="#111111" type="standard"/>'
                  '<Segment range="1-50" color="#111111" type="standard"/>'
                  '</Feature>'
                  '<Feature recentID="1" name="wrapB" type="misc_feature">'
                  '<Segment range="9500-10000" color="#222222" type="standard"/>'
                  '<Segment range="1-200" color="#222222" type="standard"/>'
                  '</Feature></Features>').encode("utf-8")
        def wrap(a1, b1, a2, b2):
            return CompoundLocation([FeatureLocation(a1, b1, strand=1),
                                     FeatureLocation(a2, b2, strand=1)])
        rec = SeqRecord(
            Seq("ACGT" * 2500), id="syn", name="syn",
            features=[
                # deliberately the OPPOSITE order to the packet
                SeqFeature(location=wrap(9499, 10000, 0, 200),
                           type="misc_feature", qualifiers={"label": ["x"]}),
                SeqFeature(location=wrap(9899, 10000, 0, 50),
                           type="misc_feature", qualifiers={"label": ["x"]}),
            ])
        sc._augment_dna_record_from_packets(rec, _make_minimal_dna((0x0A, packet)))
        assert [(f.qualifiers["label"][0]) for f in rec.features] == \
            ["wrapB", "wrapA"]

    def test_a_source_entry_never_claims_a_real_feature(self):
        """A `source` entry has no legitimate partner — the record's own
        `source` row is skipped by design. Left matchable, the span-only
        fallback hands it to whatever genuine feature spans the whole
        molecule, which is the "labelled `source`" symptom all over again."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("source",        "source",     "1-14098"),
            ("gRNA scaffold", "terminator", "1-14098"),
        )))
        # BioPython gave the whole-molecule terminator a label already; the
        # packet's own terminator entry is the only thing allowed to rename it.
        rec = self._rec(("parsed", "terminator", 0, 14098))
        sc._augment_dna_record_from_packets(rec, data)
        assert self._labels(rec) == [("terminator", 0, 14098, "gRNA scaffold")]

    def test_a_source_entry_alone_leaves_the_feature_untouched(self):
        """…and with no real entry to claim it, the feature keeps its own
        parsed name rather than being renamed "source"."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("source", "source", "1-14098"),
        )))
        rec = self._rec(("gRNA scaffold", "terminator", 0, 14098))
        sc._augment_dna_record_from_packets(rec, data)
        assert self._labels(rec) == [("terminator", 0, 14098, "gRNA scaffold")]

    def test_two_entries_at_one_span_are_handed_out_in_packet_order(self):
        """Same coordinates AND same type is legal (two overlapping notes).
        Each record feature must get a DIFFERENT entry, never the same one
        twice."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("first-note",  "misc_feature", "1-100"),
            ("second-note", "misc_feature", "1-100"),
        )))
        rec = self._rec(("x", "misc_feature", 0, 100),
                        ("x", "misc_feature", 0, 100))
        sc._augment_dna_record_from_packets(rec, data)
        assert [t[3] for t in self._labels(rec)] == ["first-note", "second-note"]

    def test_span_match_wins_when_the_type_string_differs(self):
        """BioPython does not always reproduce the packet's `type` verbatim.
        A lone entry at exactly these coordinates is still this feature's
        name, so the span-only fallback claims it — deliberately."""
        data = _make_minimal_dna((0x0A, self._features_packet(
            ("Ori*", "rep_origin", "1-600"),
        )))
        rec = self._rec(("x", "misc_feature", 0, 600))
        sc._augment_dna_record_from_packets(rec, data)
        assert self._labels(rec) == [("misc_feature", 0, 600, "Ori*")]

    def test_multi_segment_feature_matches_its_compound_location(self):
        """An origin-spanning feature is one <Feature> with two <Segment>s and
        one BioPython CompoundLocation."""
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        from Bio.Seq import Seq
        # Two <Segment>s make the span (0, 10000), matching what BioPython
        # reports for the CompoundLocation below.
        packet = ('<?xml version="1.0"?><Features nextValidID="1">'
                  '<Feature recentID="0" name="wrapped" type="misc_feature">'
                  '<Segment range="9500-10000" color="#33ccff" type="standard"/>'
                  '<Segment range="1-200" color="#33ccff" type="standard"/>'
                  '</Feature></Features>').encode("utf-8")
        loc = CompoundLocation([FeatureLocation(9499, 10000, strand=1),
                                FeatureLocation(0, 200, strand=1)])
        rec = SeqRecord(Seq("ACGT" * 2500), id="syn", name="syn",
                        features=[SeqFeature(location=loc, type="misc_feature",
                                             qualifiers={"label": ["x"]})])
        sc._augment_dna_record_from_packets(rec, _make_minimal_dna((0x0A, packet)))
        assert rec.features[0].qualifiers["label"] == ["wrapped"]


# ══════════════════════════════════════════════════════════════════════════════
# Re-saving a library entry must not drop its entry-level state. [INV-179]
# ══════════════════════════════════════════════════════════════════════════════
#
# `map_mode` (the per-plasmid linear/circular view) and `alignments` (saved
# BLAST / Plasmidsaurus results) live on the LIBRARY ENTRY, not the record — so
# a rebuilt record carried neither, and `add_entry` overwrote them with a fresh
# dict on every save. Surfaced by the 2026-08-05 label repair: re-saving 155
# entries cost 3 of them their view mode and one its alignments.


@pytest.mark.asyncio
class TestReSavePreservesEntryState:
    """Drives the REAL `LibraryPanel.add_entry` through a live app — the
    preservation logic is exercised, not stood in for."""

    @staticmethod
    def _rec(seq="ACGTACGTACGTACGTACGT", rid="P"):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        r = SeqRecord(Seq(seq), id=rid, name=rid)
        r.annotations = {"molecule_type": "DNA"}
        return r

    async def _resave_with(self, tiny_record, isolated_library, mutate,
                           second=None):
        """Save once, stamp entry-level state, then re-save under the same id."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one(sc.LibraryPanel)
            lib.add_entry(self._rec())
            await pilot.pause(0.1)
            # stamp the entry-level state the UI would have written
            entries = sc._load_library()
            for e in entries:
                if e.get("id") == "P":
                    mutate(e)
            sc._save_library(entries)
            await pilot.pause(0.1)
            lib.add_entry(second() if second else self._rec("ACGT" * 8))
            await pilot.pause(0.1)
            out = sc._find_library_entry_by_id("P")
            app.exit()
            return out

    async def test_map_mode_survives_a_resave(self, tiny_record, isolated_library):
        e = await self._resave_with(tiny_record, isolated_library,
                                    lambda x: x.__setitem__("map_mode", "linear"))
        assert e and e.get("map_mode") == "linear"

    async def test_alignments_survive_a_resave(self, tiny_record, isolated_library):
        al = [{"id": "a1", "label": "1 1A", "target_id": "X"}]
        e = await self._resave_with(tiny_record, isolated_library,
                                    lambda x: x.__setitem__("alignments", al))
        assert e and e.get("alignments") == al

    async def test_a_live_view_toggle_wins_over_the_stored_value(
            self, tiny_record, isolated_library):
        def _toggled():
            r = self._rec("ACGT" * 8)
            r._tui_map_mode = "circular"       # user just pressed 'v'
            return r
        e = await self._resave_with(tiny_record, isolated_library,
                                    lambda x: x.__setitem__("map_mode", "linear"),
                                    second=_toggled)
        assert e and e.get("map_mode") == "circular"

    async def test_a_brand_new_entry_gets_neither_key(
            self, tiny_record, isolated_library):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.query_one(sc.LibraryPanel).add_entry(self._rec(rid="FRESH"))
            await pilot.pause(0.1)
            e = sc._find_library_entry_by_id("FRESH")
            app.exit()
        assert e is not None
        assert "map_mode" not in e
        assert "alignments" not in e
