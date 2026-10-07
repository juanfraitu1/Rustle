#!/usr/bin/env python3
"""Tests of table.py's parsing helpers (2026-10-07 verification: the wall-time regex matched nothing, so every time column was empty).
Run: python3 -m unittest bench/seed_pool/test_table.py"""
import json
import os
import sys
import tempfile
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import table as T  # noqa: E402

TIME_V = """\tCommand being timed: "bash rlock.sh heavy bash rustle_pipeline.sh assemble"
\tUser time (seconds): 71.08
\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:24.19
\tMaximum resident set size (kbytes): 746572
\tExit status: 0
"""


class WallTests(unittest.TestCase):
    def write(self, text):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "x.stderr")
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def test_minutes_and_seconds(self):
        secs, gb = T.wall(self.write(TIME_V))
        self.assertAlmostEqual(secs, 24.19)
        self.assertAlmostEqual(gb, 0.746572)

    def test_hours_minutes_seconds(self):
        secs, _ = T.wall(self.write(TIME_V.replace("0:24.19", "1:02:03.50")))
        self.assertAlmostEqual(secs, 3723.5)

    def test_a_missing_file_or_a_log_without_the_line_gives_none(self):
        self.assertEqual(T.wall("/nonexistent/x.stderr"), (None, None))
        self.assertEqual(T.wall(self.write("nothing here\n")), (None, None))


class TargetRowTests(unittest.TestCase):
    FS = {"truths": {"soto": {"per_family": [
        {"family_id": "ID_154", "n_truth": "14", "hit": "13", "sens": "0.928571", "prec": "0.812500", "f": "0.866667"},
        {"family_id": "ID_2", "n_truth": "3", "hit": "0", "sens": "0.000000", "prec": "0.000000", "f": "0.000000"}]}}}

    def test_the_named_families_are_picked_in_the_requested_order(self):
        rows = T.target_rows(self.FS, {"soto": ["ID_2", "ID_154"]})
        self.assertEqual([r["family_id"] for r in rows], ["ID_2", "ID_154"])
        self.assertEqual((rows[1]["hit"], rows[1]["n_truth"]), (13, 14))
        self.assertAlmostEqual(rows[1]["f"], 0.866667)

    def test_a_family_that_is_absent_is_reported_not_skipped(self):
        rows = T.target_rows(self.FS, {"soto": ["ID_999"]})
        self.assertEqual(rows[0]["family_id"], "ID_999")
        self.assertIsNone(rows[0]["f"])

    def test_the_target_option_is_parsed(self):
        self.assertEqual(T.parse_targets("compara:CF153;u2:ID_154,ID_149;soto:ID_154"),
                         {"compara": ["CF153"], "u2": ["ID_154", "ID_149"], "soto": ["ID_154"]})
        self.assertEqual(T.parse_targets(""), {})


class BannerTests(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        with open(f"{self.d}/comp.json", "w") as fh:
            fh.write("{}")

    def gates(self, obj, newer=True):
        p = f"{self.d}/gates.json"
        with open(p, "w") as fh:
            json.dump(obj, fh)
        t = os.path.getmtime(f"{self.d}/comp.json") + (10 if newer else -10)
        os.utime(p, (t, t))

    def test_missing_gates_are_announced(self):
        self.assertIn("not run", T.gate_banner(self.d))

    def test_a_failed_gate_is_announced(self):
        self.gates({"G0": {"ok": True, "note": ""}, "G2": {"ok": False, "note": "x"}})
        self.assertIn("FAILED: G2", T.gate_banner(self.d))

    def test_gates_older_than_the_last_scoring_are_announced(self):
        self.gates({"G0": {"ok": True, "note": ""}}, newer=False)
        self.assertIn("older than", T.gate_banner(self.d))

    def test_passing_current_gates_say_so(self):
        self.gates({"G0": {"ok": True, "note": ""}, "G1": {"ok": None, "note": "n/a"}})
        self.assertIn("pass", T.gate_banner(self.d))


if __name__ == "__main__":
    unittest.main()
