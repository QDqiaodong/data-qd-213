"""Behavioral tests for the SpectrumArchive count-conserving delivery API."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from wavelength_calibration import SpectrumArchive, process

SRC = str(Path(__file__).resolve().parents[1] / "src")

# Worked example from the contract: pixels 0/10/20/30 -> 400/421/444/469,
# exact fit W(p) = 400 + 2p + p^2/100 with zero residuals.
A_LINES = [
    {"id": "l0", "pixel": "0", "wavelength_nm": "400"},
    {"id": "l1", "pixel": "10", "wavelength_nm": "421"},
    {"id": "l2", "pixel": "20", "wavelength_nm": "444"},
    {"id": "l3", "pixel": "30", "wavelength_nm": "469"},
]
A_DATA = {"max_residual_nm": "0", "reference_lines": A_LINES}

B_LINES = [
    {"id": "b0", "pixel": "0", "wavelength_nm": "500"},
    {"id": "b1", "pixel": "10", "wavelength_nm": "521"},
    {"id": "b2", "pixel": "20", "wavelength_nm": "544"},
    {"id": "b3", "pixel": "30", "wavelength_nm": "569"},
]
B_DATA = {"max_residual_nm": "0", "reference_lines": B_LINES}

EA_PIXELS = ["0", "10", "20"]
EA_COUNTS = ["21", "23"]
EA_WAVES = ["400", "410", "421", "430", "444"]


def shifted(lines, delta):
    return {
        "max_residual_nm": "0",
        "reference_lines": [
            {"id": l["id"], "pixel": l["pixel"], "wavelength_nm": str(int(l["wavelength_nm"]) + delta)}
            for l in lines
        ],
    }


class ArchiveTestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = self._tmp.name
        self.arc = SpectrumArchive(self.root)

    def tearDown(self):
        self.arc.close()
        self._tmp.cleanup()

    def reopen(self):
        self.arc.close()
        self.arc = SpectrumArchive(self.root)
        return self.arc

    def make_a(self):
        return self.arc.set_calibration("A", 0, A_DATA)

    def make_ea(self):
        self.make_a()
        self.arc.add_exposure("EA", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)
        return self.arc.build("EA")


class CalibrationApiTests(ArchiveTestCase):
    def test_create_returns_full_evidence(self):
        rec = self.make_a()
        self.assertEqual(rec["slot_id"], "A")
        self.assertEqual(rec["revision"], 1)
        self.assertEqual(rec["input"], A_DATA)
        self.assertEqual(rec["coefficients"], {"c0": "400", "c1": "2", "c2": "1/100"})
        self.assertEqual(rec["max_abs_residual_nm"], "0")
        self.assertEqual(rec["status"], "accepted")
        self.assertEqual([r["id"] for r in rec["reference_residuals"]], ["l0", "l1", "l2", "l3"])
        # endpoints were used as the query pixels
        self.assertEqual(
            rec["calibrated"],
            [
                {"pixel": "0", "wavelength_nm": "400"},
                {"pixel": "30", "wavelength_nm": "469"},
            ],
        )

    def test_fit_delegates_to_existing_process(self):
        rec = self.make_a()
        direct = process({**A_DATA, "query_pixels": ["0", "30"]})
        for key in ("coefficients", "reference_residuals", "max_abs_residual_nm", "status", "calibrated"):
            self.assertEqual(rec[key], direct[key])

    def test_revision_sequence_and_history(self):
        self.make_a()
        corrected = shifted(A_LINES, 0)
        corrected["reference_lines"][2]["id"] = "l2-renamed"
        rec2 = self.arc.set_calibration("A", 1, corrected)
        self.assertEqual(rec2["revision"], 2)
        self.assertEqual(rec2["input"]["reference_lines"][2]["id"], "l2-renamed")
        # history is preserved verbatim and current defaults to the newest
        self.assertEqual(self.arc.get_calibration("A")["revision"], 2)
        rev1 = self.arc.get_calibration("A", 1)
        self.assertEqual(rev1["input"]["reference_lines"][2]["id"], "l2")
        self.assertEqual(rev1["revision"], 1)
        with self.assertRaises(ValueError):
            self.arc.get_calibration("A", 3)
        with self.assertRaises(ValueError):
            self.arc.get_calibration("missing")

    def test_version_conflict_rejected(self):
        self.make_a()
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", 0, A_DATA)  # already exists
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", 2, A_DATA)  # not current
        with self.assertRaises(ValueError):
            self.arc.set_calibration("B", 1, B_DATA)  # new lineage needs 0
        self.assertEqual(self.arc.get_calibration("A")["revision"], 1)

    def test_residual_rejection_consumes_nothing(self):
        bad = shifted(A_LINES, 0)
        bad["reference_lines"][2]["wavelength_nm"] = "445"
        bad["max_residual_nm"] = "0"
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", 0, bad)
        with self.assertRaises(ValueError):
            self.arc.get_calibration("A")
        rec = self.arc.set_calibration("A", 0, A_DATA)  # revision 1 still available
        self.assertEqual(rec["revision"], 1)

    def test_zero_derivative_candidate_rejected(self):
        # Increasing reference lines, old fit accepts with zero residuals, but
        # W(p) = 400 + p^2 has derivative 0 at pixel 0.
        data = {
            "max_residual_nm": "0",
            "reference_lines": [
                {"id": "a", "pixel": "0", "wavelength_nm": "400"},
                {"id": "b", "pixel": "1", "wavelength_nm": "401"},
                {"id": "c", "pixel": "2", "wavelength_nm": "404"},
            ],
        }
        probe = process({**data, "query_pixels": ["0", "2"]})
        self.assertEqual(probe["status"], "accepted")  # old contract unchanged
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", 0, data)
        with self.assertRaises(ValueError):
            self.arc.get_calibration("A")

    def test_nonpositive_mapped_wavelength_rejected(self):
        # Derivative stays positive on [0,30] but the fit maps pixel 0 to a
        # nonpositive wavelength; residuals pass the loose tolerance.
        data = {
            "max_residual_nm": "4",
            "reference_lines": [
                {"id": "a", "pixel": "0", "wavelength_nm": "1"},
                {"id": "b", "pixel": "10", "wavelength_nm": "2"},
                {"id": "c", "pixel": "20", "wavelength_nm": "14"},
                {"id": "d", "pixel": "30", "wavelength_nm": "15"},
            ],
        }
        probe = process({**data, "query_pixels": ["0", "30"]})
        self.assertEqual(probe["status"], "accepted")
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", 0, data)

    def test_invalid_numeric_forms_rejected(self):
        for bad in ("1/2", "1e3", "nan", "inf", "-inf", "", 5, 1.5, None, "0.1.2", "x"):
            data = {"max_residual_nm": "0", "reference_lines": [dict(l) for l in A_LINES]}
            data["reference_lines"][0]["pixel"] = bad
            with self.assertRaises(ValueError, msg=f"pixel={bad!r}"):
                self.arc.set_calibration("A", 0, data)
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", 0, {"max_residual_nm": "1/2", "reference_lines": A_LINES})
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", 0, {"max_residual_nm": "0", "reference_lines": "nope"})
        with self.assertRaises(ValueError):
            self.arc.set_calibration("", 0, A_DATA)
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", -1, A_DATA)
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", "0", A_DATA)

    def test_rejected_candidates_keep_ready_products(self):
        product = self.make_ea()
        bad_residual = shifted(A_LINES, 0)
        bad_residual["reference_lines"][1]["wavelength_nm"] = "500"
        for mutate in (
            bad_residual,  # residual rejected
            {  # zero derivative
                "max_residual_nm": "0",
                "reference_lines": [
                    {"id": "a", "pixel": "0", "wavelength_nm": "400"},
                    {"id": "b", "pixel": "1", "wavelength_nm": "401"},
                    {"id": "c", "pixel": "2", "wavelength_nm": "404"},
                ],
            },
            {"max_residual_nm": "0", "reference_lines": [A_LINES[0]]},  # too few lines
        ):
            with self.assertRaises(ValueError):
                self.arc.set_calibration("A", 1, mutate)
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", 5, A_DATA)  # version conflict
        info = self.arc.inspect_exposure("EA")
        self.assertEqual(info["state"], "ready")
        self.assertEqual(info["last_product"], product)
        self.assertEqual(self.arc.get_calibration("A")["revision"], 1)


class ExposureTests(ArchiveTestCase):
    def test_register_and_reopen(self):
        self.make_a()
        rec = self.arc.add_exposure("EA", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)
        self.assertEqual(
            rec,
            {
                "exposure_id": "EA",
                "slot_id": "A",
                "pixel_edges": EA_PIXELS,
                "counts": EA_COUNTS,
                "wavelength_edges": EA_WAVES,
            },
        )
        self.reopen()
        info = self.arc.inspect_exposure("EA")
        self.assertEqual(info["state"], "unbuilt")
        self.assertIsNone(info["last_product"])
        self.assertEqual(info["counts"], EA_COUNTS)

    def test_duplicate_id_fails_even_with_same_content(self):
        self.make_a()
        self.arc.add_exposure("EA", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)
        with self.assertRaises(ValueError):
            self.arc.add_exposure("EA", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)

    def test_unknown_slot_fails(self):
        with self.assertRaises(ValueError):
            self.arc.add_exposure("EA", "ghost", EA_PIXELS, EA_COUNTS, EA_WAVES)

    def test_registration_allows_uncovered_edges(self):
        self.make_a()  # calibrated pixels [0,30] -> [400,469]
        self.arc.add_exposure("WIDE", "A", ["0", "31"], ["5"], ["400", "470"])
        info = self.arc.inspect_exposure("WIDE")
        self.assertEqual(info["state"], "unbuilt")

    def test_validation(self):
        self.make_a()
        bad_calls = [
            ("", "A", EA_PIXELS, EA_COUNTS, EA_WAVES),  # empty id
            ("E1", "A", ["0"], ["1"], EA_WAVES),  # too few pixel edges
            ("E2", "A", EA_PIXELS, ["1"], EA_WAVES),  # counts length mismatch
            ("E3", "A", EA_PIXELS, EA_COUNTS, ["400"]),  # too few wave edges
            ("E4", "A", ["0", "5", "5"], ["1", "1"], EA_WAVES),  # not increasing
            ("E5", "A", ["-1", "10", "20"], EA_COUNTS, EA_WAVES),  # negative pixel
            ("E6", "A", EA_PIXELS, EA_COUNTS, ["0", "410", "421", "430", "444"]),  # zero wave
            ("E7", "A", EA_PIXELS, ["1", "-1"], EA_WAVES),  # negative count
            ("E8", "A", ["0", "10", "20"], ["1", "1"], ["400", "420", "420", "430", "444"]),
            ("E9", "A", ["0", "1/2"], ["1"], ["400", "410"]),  # fraction literal
            ("E10", "A", ["0", 10], ["1"], ["400", "410"]),  # python number
        ]
        for args in bad_calls:
            with self.assertRaises(ValueError, msg=args[0]):
                self.arc.add_exposure(*args)
        # nothing was registered
        for args in bad_calls:
            if args[0]:
                with self.assertRaises(ValueError):
                    self.arc.inspect_exposure(args[0])

    def test_bin_count_limit(self):
        self.make_a()
        with self.assertRaises(ValueError):
            self.arc.add_exposure("BIG", "A", [str(i) for i in range(2002)], ["1"] * 2001, ["400", "469"])
        with self.assertRaises(ValueError):
            self.arc.add_exposure("BIGW", "A", ["0", "10"], ["1"], [str(i + 400) for i in range(2002)])


class BuildTests(ArchiveTestCase):
    def test_worked_example_exact(self):
        product = self.make_ea()
        self.assertEqual(product["product_id"], "P000001")
        self.assertEqual(product["exposure_id"], "EA")
        self.assertEqual(product["slot_id"], "A")
        self.assertEqual(product["calibration_revision"], 1)
        self.assertEqual(product["pixel_edges"], EA_PIXELS)
        self.assertEqual(product["counts"], EA_COUNTS)
        self.assertEqual(product["wavelength_edges"], EA_WAVES)
        self.assertEqual(product["source_wavelength_edges"], ["400", "421", "444"])
        self.assertEqual(product["target_counts"], ["10", "11", "9", "14"])
        self.assertEqual(product["dropped_counts"], "0")
        self.assertEqual(product["calibration"]["coefficients"]["c1"], "2")
        self.assertEqual(product["calibration"]["revision"], 1)

    def test_narrowed_target_drops_rest(self):
        self.make_a()
        self.arc.add_exposure("EN", "A", EA_PIXELS, EA_COUNTS, ["410", "430"])
        product = self.arc.build("EN")
        self.assertEqual(product["target_counts"], ["20"])
        self.assertEqual(product["dropped_counts"], "24")

    def test_non_integer_shares_are_exact(self):
        self.make_a()
        self.arc.add_exposure("EF", "A", ["0", "10"], ["5"], ["400", "410.5"])
        product = self.arc.build("EF")
        # density 5/21 per nm over [400,421]; 10.5 nm overlap -> 5/2
        self.assertEqual(product["target_counts"], ["5/2"])
        self.assertEqual(product["dropped_counts"], "5/2")

    def test_zero_counts(self):
        self.make_a()
        self.arc.add_exposure("EZ", "A", EA_PIXELS, ["0", "0"], EA_WAVES)
        product = self.arc.build("EZ")
        self.assertEqual(product["target_counts"], ["0", "0", "0", "0"])
        self.assertEqual(product["dropped_counts"], "0")

    def test_target_edges_coinciding_with_source_edges(self):
        self.make_a()
        self.arc.add_exposure("EE", "A", EA_PIXELS, EA_COUNTS, ["400", "421", "444"])
        product = self.arc.build("EE")
        self.assertEqual(product["target_counts"], ["21", "23"])
        self.assertEqual(product["dropped_counts"], "0")

    def test_conservation_invariant(self):
        self.make_a()
        from fractions import Fraction

        self.arc.add_exposure("EC", "A", ["0", "10", "20"], ["7", "11"], ["403.5", "417", "440.25"])
        product = self.arc.build("EC")
        total = sum(Fraction(x) for x in product["target_counts"]) + Fraction(product["dropped_counts"])
        self.assertEqual(total, Fraction(18))

    def test_coverage_failures_leave_nothing(self):
        self.make_a()
        self.arc.add_exposure("WIDE", "A", ["0", "31"], ["5"], ["400", "470"])
        self.arc.add_exposure("FAR", "A", EA_PIXELS, EA_COUNTS, ["400", "445"])
        with self.assertRaises(ValueError):
            self.arc.build("WIDE")  # pixel 31 outside [0,30]
        with self.assertRaises(ValueError):
            self.arc.build("FAR")  # 445 beyond mapped 444
        self.assertEqual(self.arc.inspect_exposure("WIDE")["state"], "unbuilt")
        self.assertEqual(self.arc.inspect_exposure("FAR")["state"], "unbuilt")
        # failed builds consumed no product identity
        self.arc.add_exposure("OK", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)
        self.assertEqual(self.arc.build("OK")["product_id"], "P000001")

    def test_rebuild_same_revision_returns_same_product(self):
        product = self.make_ea()
        again = self.arc.build("EA")
        self.assertEqual(again, product)
        self.reopen()
        self.assertEqual(self.arc.build("EA"), product)

    def test_correction_forces_new_product_identity(self):
        product = self.make_ea()
        renamed = shifted(A_LINES, 0)
        renamed["reference_lines"][0]["id"] = "l0-renamed"
        self.arc.set_calibration("A", 1, renamed)
        rebuilt = self.arc.build("EA")
        self.assertNotEqual(rebuilt["product_id"], product["product_id"])
        self.assertEqual(rebuilt["calibration_revision"], 2)
        # numerically identical content, but a new evidence version
        self.assertEqual(rebuilt["target_counts"], product["target_counts"])
        self.assertEqual(rebuilt["calibration"]["input"]["reference_lines"][0]["id"], "l0-renamed")
        # the old product is still frozen in the archive
        self.assertEqual(self.arc.inspect_exposure("EA")["last_product"], rebuilt)

    def test_unknown_exposure(self):
        with self.assertRaises(ValueError):
            self.arc.build("ghost")
        with self.assertRaises(ValueError):
            self.arc.inspect_exposure("ghost")


class StaleAndIsolationTests(ArchiveTestCase):
    def setUp(self):
        super().setUp()
        self.make_a()
        self.arc.set_calibration("B", 0, B_DATA)
        self.arc.add_exposure("EA", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)
        self.arc.add_exposure("EA2", "A", EA_PIXELS, ["1", "2"], EA_WAVES)
        self.arc.add_exposure("EB", "B", EA_PIXELS, ["10", "14"], ["500", "521", "544"])
        self.arc.build("EA")
        self.arc.build("EA2")
        self.arc.build("EB")

    def correct_a_by_rename(self):
        renamed = shifted(A_LINES, 0)
        renamed["reference_lines"][2]["id"] = "l2-renamed"
        return self.arc.set_calibration("A", 1, renamed)

    def test_correction_stales_only_own_slot(self):
        self.correct_a_by_rename()
        self.assertEqual(self.arc.inspect_exposure("EA")["state"], "stale")
        self.assertEqual(self.arc.inspect_exposure("EA2")["state"], "stale")
        self.assertEqual(self.arc.inspect_exposure("EB")["state"], "ready")

    def test_batch_fails_whole_and_single_slot_still_delivers(self):
        self.correct_a_by_rename()
        with self.assertRaises(ValueError):
            self.arc.release("req-batch", ["EA", "EB"])
        # request_id was not consumed by the failed batch
        solo = self.arc.release("req-eb", ["EB"])
        self.assertEqual([p["exposure_id"] for p in solo["products"]], ["EB"])

    def test_rebuild_restores_only_itself(self):
        self.correct_a_by_rename()
        old_product = self.arc.inspect_exposure("EA")["last_product"]
        rebuilt = self.arc.build("EA")
        self.assertNotEqual(rebuilt["product_id"], old_product["product_id"])
        self.assertEqual(rebuilt["calibration_revision"], 2)
        self.assertEqual(self.arc.inspect_exposure("EA")["state"], "ready")
        # the sibling exposure was not revived along with EA
        self.assertEqual(self.arc.inspect_exposure("EA2")["state"], "stale")
        # the previously failed batch request now succeeds
        batch = self.arc.release("req-batch", ["EA", "EB"])
        self.assertEqual([p["exposure_id"] for p in batch["products"]], ["EA", "EB"])

    def test_wavelength_shift_breaks_coverage_but_keeps_evidence(self):
        self.correct_a_by_rename()
        rebuilt = self.arc.build("EA")
        self.arc.set_calibration("A", 2, shifted(A_LINES, 1))  # all wavelengths +1
        # mapped source range now starts at 401; EA's fixed target starts at 400
        with self.assertRaises(ValueError):
            self.arc.build("EA")
        info = self.arc.inspect_exposure("EA")
        self.assertEqual(info["state"], "stale")
        self.assertEqual(info["last_product"], rebuilt)  # old product still auditable
        with self.assertRaises(ValueError):
            self.arc.release("req-late", ["EA"])
        # history of both calibration revisions is intact
        self.assertEqual(self.arc.get_calibration("A", 1)["input"]["reference_lines"][0]["wavelength_nm"], "400")
        self.assertEqual(self.arc.get_calibration("A", 3)["input"]["reference_lines"][0]["wavelength_nm"], "401")


class ReleaseTests(ArchiveTestCase):
    def setUp(self):
        super().setUp()
        self.make_a()
        self.arc.set_calibration("B", 0, B_DATA)
        self.arc.add_exposure("EA", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)
        self.arc.add_exposure("EB", "B", EA_PIXELS, ["10", "14"], ["500", "521", "544"])
        self.arc.build("EA")
        self.arc.build("EB")

    def test_release_snapshot_sorted_and_numbered(self):
        rel = self.arc.release("req-1", ["EB", "EA"])
        self.assertEqual(rel["release_id"], "R000001")
        self.assertEqual(rel["release_number"], 1)
        self.assertEqual(rel["request_id"], "req-1")
        self.assertEqual([p["exposure_id"] for p in rel["products"]], ["EA", "EB"])
        rel2 = self.arc.release("req-2", ["EB"])
        self.assertEqual(rel2["release_number"], 2)
        self.assertEqual(rel2["release_id"], "R000002")

    def test_replay_returns_original_verbatim(self):
        rel = self.arc.release("req-1", ["EA"])
        replay = self.arc.release("req-1", ["EA"])
        self.assertEqual(replay, rel)
        # replay does not consume a number
        self.assertEqual(self.arc.release("req-2", ["EB"])["release_number"], 2)

    def test_replay_survives_correction_and_rebuild(self):
        rel = self.arc.release("req-1", ["EA"])
        renamed = shifted(A_LINES, 0)
        renamed["reference_lines"][0]["id"] = "l0-x"
        self.arc.set_calibration("A", 1, renamed)
        # stale EA cannot start a new delivery, but the old request replays
        replay = self.arc.release("req-1", ["EA"])
        self.assertEqual(replay, rel)
        self.assertEqual(replay["products"][0]["calibration_revision"], 1)
        self.arc.build("EA")  # rebuild at revision 2
        replay2 = self.arc.release("req-1", ["EA"])
        self.assertEqual(replay2, rel)
        # a genuinely new request delivers the new product with the next number
        fresh = self.arc.release("req-3", ["EA"])
        self.assertEqual(fresh["release_number"], 2)
        self.assertEqual(fresh["products"][0]["calibration_revision"], 2)
        self.assertNotEqual(fresh["products"][0]["product_id"], rel["products"][0]["product_id"])
        # the historical delivery is frozen
        self.assertEqual(self.arc.get_release(rel["release_id"]), rel)

    def test_request_id_bound_to_exposure_set(self):
        self.arc.release("req-1", ["EA"])
        with self.assertRaises(ValueError):
            self.arc.release("req-1", ["EB"])
        with self.assertRaises(ValueError):
            self.arc.release("req-1", ["EA", "EB"])
        # order does not matter for the set
        self.arc.release("req-2", ["EA", "EB"])
        replay = self.arc.release("req-2", ["EB", "EA"])
        self.assertEqual(replay["release_number"], 2)

    def test_legality_checks_precede_replay(self):
        self.arc.release("req-1", ["EA"])
        with self.assertRaises(ValueError):
            self.arc.release("req-1", ["EA", "EA"])  # duplicates still rejected
        with self.assertRaises(ValueError):
            self.arc.release("req-1", [])
        with self.assertRaises(ValueError):
            self.arc.release("", ["EA"])
        with self.assertRaises(ValueError):
            self.arc.release("req-1", ["ghost"])

    def test_failed_request_consumes_nothing(self):
        self.arc.add_exposure("EU", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)  # unbuilt
        with self.assertRaises(ValueError):
            self.arc.release("req-x", ["EU"])
        with self.assertRaises(ValueError):
            self.arc.release("req-y", ["EA", "EU"])
        self.arc.build("EU")
        rel = self.arc.release("req-x", ["EU"])  # same request_id now succeeds
        self.assertEqual(rel["release_number"], 1)  # no gap from the failures
        rel2 = self.arc.release("req-y", ["EA", "EU"])
        self.assertEqual(rel2["release_number"], 2)

    def test_get_release_isolated_from_caller_mutation(self):
        rel = self.arc.release("req-1", ["EA", "EB"])
        rel["products"][0]["target_counts"][0] = "999"
        rel["products"].clear()
        rel["request_id"] = "forged"
        again = self.arc.get_release("R000001")
        self.assertEqual(again["request_id"], "req-1")
        self.assertEqual(len(again["products"]), 2)
        self.assertEqual(again["products"][0]["target_counts"][0], "10")
        with self.assertRaises(ValueError):
            self.arc.get_release("R999999")

    def test_release_frozen_across_reopen(self):
        rel = self.arc.release("req-1", ["EA"])
        self.reopen()
        self.assertEqual(self.arc.release("req-1", ["EA"]), rel)
        self.assertEqual(self.arc.get_release("R000001"), rel)


class FaultHookTests(ArchiveTestCase):
    def test_hook_called_once_per_committed_operation(self):
        calls = []
        self.arc.close()
        self.arc = SpectrumArchive(self.root, fault_hook=lambda ph, op: calls.append((ph, op)))
        self.arc.set_calibration("A", 0, A_DATA)
        self.arc.add_exposure("EA", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)
        self.arc.build("EA")
        self.arc.release("req-1", ["EA"])
        self.assertEqual(
            calls,
            [
                ("before_commit", "set_calibration"),
                ("after_commit", "set_calibration"),
                ("before_commit", "add_exposure"),
                ("after_commit", "add_exposure"),
                ("before_commit", "build"),
                ("after_commit", "build"),
                ("before_commit", "release"),
                ("after_commit", "release"),
            ],
        )
        # reads, build cache hits and replays do not fire the hook
        calls.clear()
        self.arc.get_calibration("A")
        self.arc.inspect_exposure("EA")
        self.arc.build("EA")
        self.arc.release("req-1", ["EA"])
        self.arc.get_release("R000001")
        self.assertEqual(calls, [])

    def test_hook_exception_propagates_unwrapped(self):
        def boom(phase, op):
            raise RuntimeError(f"boom at {phase}")

        self.arc.close()
        self.arc = SpectrumArchive(self.root, fault_hook=boom)
        with self.assertRaises(RuntimeError):
            self.arc.set_calibration("A", 0, A_DATA)
        with self.assertRaises(ValueError):
            self.arc.get_calibration("A")  # before_commit: nothing persisted

    def test_after_commit_exception_still_persists(self):
        def boom(phase, op):
            if phase == "after_commit":
                raise RuntimeError("boom")

        self.arc.close()
        self.arc = SpectrumArchive(self.root, fault_hook=boom)
        with self.assertRaises(RuntimeError):
            self.arc.set_calibration("A", 0, A_DATA)
        self.reopen()
        self.assertEqual(self.arc.get_calibration("A")["revision"], 1)

    def test_invalid_input_never_calls_hook(self):
        calls = []
        self.arc.close()
        self.arc = SpectrumArchive(self.root, fault_hook=lambda ph, op: calls.append(ph))
        with self.assertRaises(ValueError):
            self.arc.set_calibration("A", 1, A_DATA)
        with self.assertRaises(ValueError):
            self.arc.add_exposure("EA", "ghost", EA_PIXELS, EA_COUNTS, EA_WAVES)
        with self.assertRaises(ValueError):
            self.arc.release("req", [])
        self.assertEqual(calls, [])


CHILD_COMMON = """
import os, sys
sys.path.insert(0, os.environ["WC_SRC"])
from wavelength_calibration import SpectrumArchive

A_DATA = {a_data!r}
EA_PIXELS = {ea_pixels!r}
EA_COUNTS = {ea_counts!r}
EA_WAVES = {ea_waves!r}
root, phase, op = sys.argv[1], sys.argv[2], sys.argv[3]

def hook(ph, operation):
    if ph == phase and operation == op:
        os._exit(42)

arc = SpectrumArchive(root, fault_hook=hook)
"""

CHILD_SET_CALIBRATION = CHILD_COMMON + """
arc.set_calibration("A", 0, A_DATA)
"""

CHILD_RELEASE = CHILD_COMMON + """
arc.set_calibration("A", 0, A_DATA)
arc.add_exposure("EA", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)
arc.build("EA")
arc.release("req-1", ["EA"])
"""


class CrashRecoveryTests(unittest.TestCase):
    """Kill a child process at the commit boundaries and verify all-or-nothing."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def run_child(self, script, phase, op):
        env = dict(os.environ, WC_SRC=SRC)
        code = script.format(a_data=A_DATA, ea_pixels=EA_PIXELS, ea_counts=EA_COUNTS, ea_waves=EA_WAVES)
        result = subprocess.run(
            [sys.executable, "-c", code, self.root, phase, op],
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertEqual(result.returncode, 42, msg=result.stderr)

    def test_set_calibration_before_commit_is_invisible(self):
        self.run_child(CHILD_SET_CALIBRATION, "before_commit", "set_calibration")
        arc = SpectrumArchive(self.root)
        with self.assertRaises(ValueError):
            arc.get_calibration("A")
        # the lineage can still be created at revision 1 afterwards
        self.assertEqual(arc.set_calibration("A", 0, A_DATA)["revision"], 1)
        arc.close()

    def test_set_calibration_after_commit_is_complete(self):
        self.run_child(CHILD_SET_CALIBRATION, "after_commit", "set_calibration")
        arc = SpectrumArchive(self.root)
        rec = arc.get_calibration("A")
        self.assertEqual(rec["revision"], 1)
        self.assertEqual(rec["coefficients"]["c1"], "2")
        arc.close()

    def test_release_before_commit_leaves_no_trace(self):
        self.run_child(CHILD_RELEASE, "before_commit", "release")
        arc = SpectrumArchive(self.root)
        self.assertEqual(arc.inspect_exposure("EA")["state"], "ready")
        with self.assertRaises(ValueError):
            arc.get_release("R000001")
        # the request_id and the first release number are both still free
        rel = arc.release("req-1", ["EA"])
        self.assertEqual(rel["release_number"], 1)
        self.assertEqual(rel["release_id"], "R000001")
        arc.close()

    def test_release_after_commit_is_fully_visible(self):
        # The child commits the release and dies at after_commit, simulating a
        # lost response.  Reopening must see the complete delivery.
        self.run_child(CHILD_RELEASE, "after_commit", "release")
        arc = SpectrumArchive(self.root)
        # cross-slot isolation: an unrelated lineage stays ready throughout
        arc.set_calibration("B", 0, B_DATA)
        arc.add_exposure("EB", "B", EA_PIXELS, ["10", "14"], ["500", "521", "544"])
        arc.build("EB")
        # correct A after reopen; EA becomes stale, EB stays ready
        renamed = json.loads(json.dumps(A_DATA))
        renamed["reference_lines"][0]["id"] = "l0-x"
        arc.set_calibration("A", 1, renamed)
        self.assertEqual(arc.inspect_exposure("EA")["state"], "stale")
        self.assertEqual(arc.inspect_exposure("EB")["state"], "ready")
        # retrying the crashed request replays the original delivery verbatim:
        # same release_id, release_number, product identity and evidence
        rel = arc.release("req-1", ["EA"])
        self.assertEqual(rel["release_id"], "R000001")
        self.assertEqual(rel["release_number"], 1)
        self.assertEqual(rel["products"][0]["product_id"], "P000001")
        self.assertEqual(rel["products"][0]["calibration_revision"], 1)
        self.assertEqual(arc.get_release("R000001"), rel)
        # a new request for stale EA fails; req-1 with a different set fails
        with self.assertRaises(ValueError):
            arc.release("req-2", ["EA"])
        with self.assertRaises(ValueError):
            arc.release("req-1", ["EB"])
        with self.assertRaises(ValueError):
            arc.release("req-9", ["EA", "EB"])
        # EB alone still delivers, taking the next consecutive number
        self.assertEqual(arc.release("req-eb", ["EB"])["release_number"], 2)
        # rebuilding EA does not alter the historical delivery
        arc.build("EA")
        self.assertEqual(arc.release("req-1", ["EA"]), rel)
        fresh = arc.release("req-3", ["EA"])
        self.assertEqual(fresh["release_number"], 3)
        self.assertEqual(fresh["products"][0]["calibration_revision"], 2)
        self.assertEqual(arc.get_release("R000001"), rel)
        arc.close()

    def test_add_exposure_and_build_survive_reopen(self):
        arc = SpectrumArchive(self.root)
        arc.set_calibration("A", 0, A_DATA)
        arc.add_exposure("EA", "A", EA_PIXELS, EA_COUNTS, EA_WAVES)
        arc.close()
        arc = SpectrumArchive(self.root)
        self.assertEqual(arc.inspect_exposure("EA")["state"], "unbuilt")
        product = arc.build("EA")
        arc.close()
        arc = SpectrumArchive(self.root)
        info = arc.inspect_exposure("EA")
        self.assertEqual(info["state"], "ready")
        self.assertEqual(info["last_product"], product)
        arc.close()


if __name__ == "__main__":
    unittest.main()
