#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_fate.py"""
import collections
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fate as F  # noqa: E402
from common import Rec  # noqa: E402


def rec(primary=True, mapq=60, score=1000, qcov=1.0, ref="A", start=0, end=100, de=0.01):
    return Rec(primary, mapq, score, qcov, ref, start, end, de)


class Verdict(unittest.TestCase):
    def fr(self, a, u, t):
        return {"absorbed": a, "unmapped": u, "tied": t}

    def test_refuted(self):
        self.assertEqual(F.verdict(self.fr(0.9, 0.05, 0.05)), "REFUTED")

    def test_partly(self):
        self.assertEqual(F.verdict(self.fr(0.7, 0.2, 0.1)), "PARTLY")

    def test_supported_when_lost_reads_reach_half(self):
        self.assertEqual(F.verdict(self.fr(0.4, 0.3, 0.3)), "SUPPORTED")

    def test_between_bars(self):
        self.assertEqual(F.verdict(self.fr(0.85, 0.0, 0.15)), "BETWEEN_BARS")

    def test_zero_reads_has_no_verdict(self):
        self.assertEqual(F.verdict(F.fractions(collections.Counter(), 0)), "NO_READS")


class Bar(unittest.TestCase):
    def test_bar_applies_to_large_absent_loci_only(self):
        self.assertTrue(F.bar_applies("catalog", 20))
        self.assertFalse(F.bar_applies("catalog", 19))
        self.assertFalse(F.bar_applies("sex", 500))
        self.assertFalse(F.bar_applies("lrpap1_desc", 500))
        self.assertTrue(F.bar_applies("lrpap1", 83))


class Regions(unittest.TestCase):
    def test_hits_within_the_gap_form_one_region(self):
        hits = [("A", 100, 200), ("A", 250, 400), ("A", 9000, 9100), ("B", 5, 50), ("A", 150, 260)]
        self.assertEqual(F.cluster_regions(hits, gap=100), [("A", 100, 400, 3), ("A", 9000, 9100, 1), ("B", 5, 50, 1)])

    def test_largest_first_in_the_summary(self):
        self.assertEqual(F.cluster_regions([("A", 1, 2), ("B", 1, 2), ("B", 3, 4)], gap=10, largest_first=True)[0][0], "B")


class Rows(unittest.TestCase):
    def test_counts_and_de(self):
        recs = {"r1": [rec(ref="P", start=10, end=90, de=0.002)],           # absorbed on the nearest paralog
                "r2": [rec(ref="Q", de=0.05)],                                # absorbed elsewhere
                "r3": [rec(mapq=0)],                                           # tied
                "r4": [],                                                      # unmapped
                "r5": [rec(qcov=0.5)],                                         # partial
                "r6": [rec()]}                                                 # shared, not counted for the locus
        labels = {"r1": "L", "r2": "L", "r3": "L", "r4": "L", "r5": "L", "r6": "shared", "r7": "ambiguous"}
        out = F.fate_rows(recs, labels, {"L": ("P", 0, 100)})
        self.assertEqual(out["L"]["n"], 5)
        self.assertEqual(dict(out["L"]["fates"]), {"ABSORBED_NEAREST": 1, "ABSORBED_OTHER": 1, "TIED": 1, "UNMAPPED": 1, "PARTIAL": 1})
        self.assertEqual(sorted(out["L"]["de"]), [0.002, 0.05])
        self.assertEqual(out["shared"]["n"], 1)
        self.assertNotIn("ambiguous", out)

    def test_read_missing_from_the_bam_is_skipped(self):
        self.assertEqual(F.fate_rows({}, {"x": "L"}, {}), {})

    def test_fractions_merge_partial_into_unmapped(self):
        fr = F.fractions(collections.Counter({"UNMAPPED": 1, "PARTIAL": 1, "TIED": 2, "ABSORBED_OTHER": 6}), 10)
        self.assertEqual(fr, {"absorbed": 0.6, "unmapped": 0.2, "tied": 0.2})


if __name__ == "__main__":
    unittest.main()
