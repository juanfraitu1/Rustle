#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_artifact_data.py"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import artifact_data as A  # noqa: E402


class Compact(unittest.TestCase):
    def test_reads_lose_names_and_positions(self):
        self.assertEqual(A.compact_reads([["SRR1.5", "TIED", "CM1", 100, 0.01]]), [[0, "TIED", None, None, 0.01]])

    def test_pairs_lose_names_only(self):
        self.assertEqual(A.compact_pairs([["SRR1.5", "TIED", 0.05, 60, 0.001, 39]]), [[0, "TIED", 0.05, 60, 0.001, 39]])

    def test_unmapped_summary_is_parsed(self):
        txt = ("R_unm on mat: 959 reads, 132 with a primary at query coverage >= 0.8\n  region CM054594.2:95958044-96028773: 125 reads, median MAPQ 60\n"
               "  region CM054589.2:5-9: 1 reads, median MAPQ 60\nR_unm on pat: 959 reads, 3 with a primary at query coverage >= 0.8\n"
               "  region CM1:1-2: 1 reads, median MAPQ 1")
        self.assertEqual(A.parse_unm(txt), {"total": 959, "mat": 132, "pat": 3, "mat_region": {"acc": "CM054594.2", "start": 95958044, "end": 96028773, "n": 125}})

    def test_unmapped_summary_of_nothing_is_none(self):
        self.assertIsNone(A.parse_unm(""))

    def test_shared_control_drops_its_read_list(self):
        out = A.light({"n": 3, "fates": {}, "reads": [["a", "TIED", None, None, None]] * 3})
        self.assertNotIn("reads", out)
        self.assertEqual(out["n"], 3)


if __name__ == "__main__":
    unittest.main()
