#!/usr/bin/env python3
"""Tests of paralog_support.py's placement rule. Run: python3 -m unittest bench/seed_pool/test_paralog_support.py"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import paralog_support as P  # noqa: E402

SPANS = {
    "c1": dict(chrom="chr1", lo=1000, hi=2000),
    "c2": dict(chrom="chr1", lo=5000, hi=6000),
    "c3": dict(chrom="chr2", lo=1000, hi=2000),
}


class WhereTests(unittest.TestCase):
    def test_a_primary_in_another_copys_territory_is_other_copy(self):
        self.assertEqual(P.where(("chr1", 5100, 5400), "c1", SPANS), "other_copy")

    def test_a_primary_in_the_same_copys_territory_is_same_copy(self):
        self.assertEqual(P.where(("chr1", 1100, 1400), "c1", SPANS), "same_copy")

    def test_a_primary_outside_every_territory_or_unseen_is_elsewhere(self):
        self.assertEqual(P.where(("chr1", 9000, 9500), "c1", SPANS), "elsewhere")
        self.assertEqual(P.where(None, "c1", SPANS), "elsewhere")

    def test_the_same_coordinates_on_another_contig_are_another_territory(self):
        self.assertEqual(P.where(("chr2", 1100, 1400), "c1", SPANS), "other_copy")
        self.assertEqual(P.where(("chr3", 1100, 1400), "c1", SPANS), "elsewhere")

    def test_a_primary_straddling_two_territories_counts_as_other_copy(self):
        spans = dict(SPANS, c4=dict(chrom="chr1", lo=1900, hi=2500))
        self.assertEqual(P.where(("chr1", 1500, 2200), "c1", spans), "other_copy")

    def test_touching_territories_do_not_overlap(self):
        self.assertEqual(P.where(("chr1", 2000, 2100), "c1", SPANS), "elsewhere")


if __name__ == "__main__":
    unittest.main()
