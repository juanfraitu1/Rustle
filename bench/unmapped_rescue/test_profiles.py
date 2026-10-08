#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_profiles.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import profiles as P  # noqa: E402

LINES = [
    "# target name accession query name accession hmmfrom hmm to alifrom ali to envfrom env to sq len strand E-value score bias description",
    "cl1                   -          GWFAM5               -                1    500        10       400      8      402   2000    +   1e-30  120.5   0.1  -",
    "cl1                   -          GWFAM5               -                1    300      1500      1200   1190     1510   2000    -   1e-10   40.0   0.0  -",
    "cl1                   -          GWFAM7               -                1    200       100       250     95      255   2000    +   0.0001  25.0   0.0  -",
    "# Program: nhmmer",
]


class Tblout(unittest.TestCase):
    def test_parse_skips_comments_and_reads_spans_either_strand(self):
        rows = P.parse_tblout(LINES)
        self.assertEqual(rows[0], ("cl1", "GWFAM5", 10, 400, 1e-30, 120.5))
        self.assertEqual(rows[1], ("cl1", "GWFAM5", 1200, 1500, 1e-10, 40.0))      # minus strand: from > to, normalised
        self.assertEqual(len(rows), 3)

    def test_evalue_cut(self):
        self.assertEqual(len(P.parse_tblout(LINES, max_evalue=1e-5)), 2)

    def test_cover_scores_are_the_union_of_spans_per_family(self):
        sc = P.cover_scores_from_rows(P.parse_tblout(LINES))
        self.assertEqual(sc["cl1"]["GWFAM5"], (400 - 10 + 1) + (1500 - 1200 + 1))
        self.assertEqual(sc["cl1"]["GWFAM7"], 151)

    def test_bit_scores_keep_the_best_hit_per_family(self):
        self.assertEqual(P.bit_scores_from_rows(P.parse_tblout(LINES))["cl1"], {"GWFAM5": 120.5, "GWFAM7": 25.0})


if __name__ == "__main__":
    unittest.main()
