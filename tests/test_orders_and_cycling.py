"""Oligo order sheets (plate layout) and thermocycler programs.

Both close the same kind of gap: the app designs something and the user then
retypes it into a vendor form or a machine by hand, which is where the
transcription errors happen.
"""

from __future__ import annotations

import csv

import pytest

import splicecraft as sc


class TestPlateWells:

    def test_row_major_fill(self):
        """Row-major (A1..A12, B1..) is what every vendor plate template
        assumes; the other order silently transposes the whole set."""
        wells = sc._plate_well_positions(14)
        assert [w for _p, w in wells[:13]] == [
            "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10",
            "A11", "A12", "B1"]

    def test_overflow_rolls_onto_a_second_plate(self):
        """Dropping oligo 97 silently would be worse than a second plate."""
        wells = sc._plate_well_positions(100)
        assert wells[95] == (1, "H12")
        assert wells[96] == (2, "A1")
        assert len(wells) == 100

    def test_zero_and_negative(self):
        assert sc._plate_well_positions(0) == []
        assert sc._plate_well_positions(-5) == []

    def test_plate_is_96_wells(self):
        assert sc._PLATE_SIZE == 96
        assert len({w for _p, w in sc._plate_well_positions(96)}) == 96


class TestPlateExport:

    def _primers(self, n):
        return [{"name": f"oligo{i}", "sequence": "ACGT" * 5}
                for i in range(n)]

    def test_plate_columns_and_rows(self, tmp_path):
        out = tmp_path / "plate.csv"
        res = sc._export_primers_to_csv(self._primers(3), out,
                                        order_format="plate")
        assert res["count"] == 3
        rows = list(csv.reader(out.read_text().splitlines()))
        assert rows[0] == ["Plate", "Well Position", "Name", "Sequence"]
        assert rows[1] == ["1", "A1", "oligo0", "ACGTACGTACGTACGTACGT"]
        assert rows[3][1] == "A3"

    def test_plate_export_spans_two_plates(self, tmp_path):
        out = tmp_path / "big.csv"
        sc._export_primers_to_csv(self._primers(97), out,
                                  order_format="plate")
        rows = list(csv.reader(out.read_text().splitlines()))[1:]
        assert rows[95][:2] == ["1", "H12"]
        assert rows[96][:2] == ["2", "A1"]

    def test_the_other_formats_still_work(self, tmp_path):
        for fmt, header in (("generic", ["Name", "Sequence", "Length", "Tm"]),
                            ("idt", ["Name", "Sequence", "Scale",
                                     "Purification"])):
            out = tmp_path / f"{fmt}.csv"
            sc._export_primers_to_csv(self._primers(2), out,
                                      order_format=fmt)
            rows = list(csv.reader(out.read_text().splitlines()))
            assert rows[0] == header

    def test_unknown_format_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="unknown order_format"):
            sc._export_primers_to_csv(self._primers(1), tmp_path / "x.csv",
                                      order_format="nope")


class TestPcrProgram:

    def test_extension_scales_with_product_and_polymerase(self):
        """Taq is ~1 kb/min, Q5 ~30 s/kb — a 6 kb product cannot share a
        20-second extension with a 300 bp one."""
        taq = sc._pcr_program(6000, [55.0, 56.0], polymerase="taq")
        q5 = sc._pcr_program(6000, [65.0, 66.0], polymerase="q5")
        assert taq["extend_s"] == 360
        assert q5["extend_s"] == 180
        assert taq["extend_s"] > q5["extend_s"]

    def test_extension_rounds_up(self):
        """4 seconds short of full extension is a smeared band, and no
        thermocycler cares about the difference."""
        p = sc._pcr_program(1010, [60.0, 60.0], polymerase="q5")
        assert p["extend_s"] % 5 == 0
        assert p["extend_s"] >= 1010 * 30 / 1000

    def test_minimum_extension(self):
        assert sc._pcr_program(10, [60.0], polymerase="q5")["extend_s"] >= 10

    def test_anneal_from_the_lowest_tm(self):
        p = sc._pcr_program(1000, [58.0, 66.0], polymerase="taq")
        assert p["anneal_c"] == 53          # 58 - 5, the Taq rule
        assert "Tm" in p["anneal_rule"]

    def test_hifi_uses_its_own_rule(self):
        p = sc._pcr_program(1000, [60.0, 64.0], polymerase="q5")
        assert p["anneal_c"] == 63          # 60 + 3

    def test_high_tm_hifi_anneals_at_the_extension_temp(self):
        p = sc._pcr_program(1000, [70.0, 71.0], polymerase="q5")
        assert p["anneal_c"] == 72
        assert p["two_step"] is True
        assert any("Two-step" in n for n in p["notes"])
        assert not any(st["phase"] == "Anneal" for st in p["steps"])

    def test_no_tm_says_so_rather_than_inventing_a_derivation(self):
        p = sc._pcr_program(1000, [])
        assert p["anneal_c"] == 55
        assert p["anneal_rule"] == "conventional default"
        assert any("No primer Tm" in n for n in p["notes"])

    def test_anneal_is_clamped_and_reported(self):
        p = sc._pcr_program(1000, [20.0], polymerase="taq")
        assert p["anneal_c"] == 45
        assert any("clamped" in n for n in p["notes"])

    def test_cycles_are_bounded(self):
        assert sc._pcr_program(1000, [60.0], cycles=0)["cycles"] == 1
        assert sc._pcr_program(1000, [60.0], cycles=999)["cycles"] == 45

    def test_steps_cover_the_whole_run(self):
        p = sc._pcr_program(2000, [62.0, 63.0], polymerase="phusion")
        phases = [st["phase"] for st in p["steps"]]
        assert phases[0] == "Initial denature"
        assert "Denature" in phases and "Extend" in phases
        assert phases[-1] == "Hold"
        assert p["total_minutes"] > 0

    def test_hold_can_be_omitted(self):
        p = sc._pcr_program(1000, [60.0], final_hold_c=None)
        assert all(st["phase"] != "Hold" for st in p["steps"])

    def test_unknown_polymerase_names_the_options(self):
        with pytest.raises(ValueError, match="q5"):
            sc._pcr_program(1000, [60.0], polymerase="pfu-turbo-9000")

    def test_every_step_carries_a_temperature_and_duration(self):
        p = sc._pcr_program(1500, [61.0, 62.0])
        for st in p["steps"]:
            assert isinstance(st["temp_c"], int)
            assert isinstance(st["seconds"], int)
            assert int(st.get("cycles", 1)) >= 1

    def test_text_is_typable(self):
        p = sc._pcr_program(3200, [62.1, 64.5], polymerase="q5")
        txt = sc._pcr_program_text(p)
        assert "Initial denature" in txt and "Extend" in txt
        assert "1:40" in txt                # 100 s rendered as mm:ss
        assert "x30" in txt
        assert "min total" in txt

    def test_text_tolerates_junk(self):
        assert sc._pcr_program_text({}) == ""
        assert sc._pcr_program_text(None) == ""   # type: ignore[arg-type]


class TestAgentEndpoints:
    """Agent parity for both: the GUI could order oligos and the agent
    could not."""

    class _App:
        _current_record = None

    def _call(self, name, payload):
        import splicecraft_state as st
        fn, _w = st._AGENT_HANDLERS[name]
        return fn(self._App(), payload)

    def test_endpoints_registered(self):
        import splicecraft_state as st
        assert st._AGENT_HANDLERS["export-primers"][1] is True   # write
        assert "pcr-program" in st._AGENT_HANDLERS
        assert "list-polymerases" in st._AGENT_HANDLERS

    def test_list_polymerases_shows_the_rules(self):
        out = self._call("list-polymerases", {})
        keys = {p["key"] for p in out["polymerases"]}
        assert {"taq", "phusion", "q5"} <= keys
        assert all(p["anneal_rule"] for p in out["polymerases"])

    def test_pcr_program_endpoint(self):
        out = self._call("pcr-program",
                         {"product_bp": 3200, "primer_tms": [62.1, 64.5]})
        assert out["anneal_c"] == 65
        assert "Initial denature" in out["text"]

    @pytest.mark.parametrize("payload,needle", [
        ({}, "product_bp"),
        ({"product_bp": 0}, ">= 1"),
        ({"product_bp": 100, "polymerase": "nope"}, "unknown polymerase"),
        ({"product_bp": 100, "primer_tms": ["hot"]}, "must be numbers"),
    ])
    def test_pcr_program_bad_payloads(self, payload, needle):
        body, code = self._call("pcr-program", payload)
        assert code == 400 and needle in body["error"]

    def test_export_primers_rejects_a_non_csv_path(self):
        body, code = self._call("export-primers", {"path": "/tmp/x.txt"})
        assert code == 400 and ".csv" in body["error"]

    def test_export_primers_needs_a_path(self):
        body, code = self._call("export-primers", {})
        assert code == 400 and "path" in body["error"]

    def test_export_primers_rejects_an_unknown_format(self):
        body, code = self._call("export-primers",
                                {"path": "/tmp/x.csv", "format": "nope"})
        assert code == 400 and "plate" in body["error"]
