#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_pool.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pool as P  # noqa: E402

PANEL = [dict(fam="F1", mask_gene="F1:1", mask=["chrA", 1000, 2000], keep_gene="F1:0", keep=["chrB", 500, 900]),
         dict(fam="F2", mask_gene="F2:1", mask=["chrA", 5000, 6000], keep_gene="F2:0", keep=["chrC", 0, 100])]


class Label(unittest.TestCase):
    def test_read_starting_in_an_erased_interval_is_D_of_that_family(self):
        rows = {"r1": ("chrA", 1500), "r2": ("chrA", 5999), "r3": ("chrB", 600), "r4": ("chrA", 2000), "r5": ("chrX", 5)}
        lab = P.label_reads(rows, PANEL)
        self.assertEqual(lab["r1"], ("F1", "F1:1", "D"))
        self.assertEqual(lab["r2"], ("F2", "F2:1", "D"))
        self.assertEqual(lab["r3"], ("F1", "F1:0", "S"))
        self.assertEqual(lab["r4"], (None, None, "other"))           # end is exclusive
        self.assertEqual(lab["r5"], (None, None, "other"))

    def test_unknown_read_is_other(self):
        self.assertEqual(P.label_reads({}, PANEL), {})


if __name__ == "__main__":
    unittest.main()
