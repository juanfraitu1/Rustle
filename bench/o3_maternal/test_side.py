#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_side.py"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import side as S  # noqa: E402
import testutil  # noqa: E402
from common import Rec  # noqa: E402


class Track(unittest.TestCase):
    def run_one(self, cigar, start=0, lo=0, hi=100, n=2):
        cov, mis = [0] * n, [0] * n
        S.add_alignment(cov, mis, start, cigar, lo, hi)
        return cov, mis

    def test_mismatch_lands_in_its_bin(self):
        self.assertEqual(self.run_one([(7, 10), (8, 1), (7, 89)]), ([50, 50], [1, 0]))

    def test_insertion_adds_a_mismatch_without_coverage(self):
        self.assertEqual(self.run_one([(7, 60), (1, 3), (7, 40)]), ([50, 50], [0, 1]))

    def test_deletion_is_covered_and_mismatched(self):
        self.assertEqual(self.run_one([(7, 40), (2, 20), (7, 40)]), ([50, 50], [10, 10]))

    def test_intron_covers_nothing(self):
        self.assertEqual(self.run_one([(7, 20), (3, 500), (7, 20)], hi=1000), ([20, 20], [0, 0]))

    def test_outside_the_interval_is_ignored(self):
        self.assertEqual(self.run_one([(7, 100)], start=500), ([0, 0], [0, 0]))


class Pairs(unittest.TestCase):
    def rec(self, de, mapq, primary=True):
        return Rec(primary, mapq, 100, 1.0, "A", 0, 10, de)

    def test_pairs_need_a_primary_on_both(self):
        ref = {"a": [self.rec(0.05, 60)], "b": [self.rec(0.05, 60)], "c": []}
        oth = {"a": [self.rec(0.001, 39)], "b": [], "c": [self.rec(0.0, 60)]}
        rows = S.pair_rows(["a", "b", "c"], ref, oth, {"a": "ABSORBED_OTHER", "b": "TIED", "c": "UNMAPPED"})
        self.assertEqual(rows, [["a", "ABSORBED_OTHER", 0.05, 60, 0.001, 39]])


class BamTrack(unittest.TestCase):
    def test_track_counts_only_the_named_primaries(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.bam")
            testutil.write_bam(p, {"chr1_mat_hsa1": 5000}, [
                dict(name="a", ref="chr1_mat_hsa1", start=100, cigar="50=1X49="),
                dict(name="b", ref="chr1_mat_hsa1", start=100, cigar="100="),                    # not named
                dict(name="a", flag=256, ref="chr1_mat_hsa1", start=100, cigar="100=", mapq=0),   # secondary: ignored
            ])
            t = S.track(p, {"a"}, ("CM1", 100, 200), {("mat", "1"): "CM1"}, nbins=2)
            self.assertEqual(t["cov"], [50, 50])
            self.assertEqual(t["mis"], [0, 1])
            self.assertEqual(t["target"], ["CM1", 100, 200])


if __name__ == "__main__":
    unittest.main()
