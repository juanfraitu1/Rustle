#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_lrpap1.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lrpap1 as L  # noqa: E402

ROWS = [
    dict(cid="c01", pat_acc="P1", pat_s="100", pat_e="200", mat_acc="M1", mat_s="1000", mat_e="1100"),
    dict(cid="c03", pat_acc="P2", pat_s="500", pat_e="600", mat_acc="M1", mat_s="1000", mat_e="1100"),   # lifts to the same mat site as c01
    dict(cid="c07", pat_acc="P3", pat_s="900", pat_e="950", mat_acc="", mat_s="", mat_e=""),             # chrY: no mat site
    dict(cid="p12", pat_acc="P1", pat_s="300", pat_e="400", mat_acc="M2", mat_s="2000", mat_e="2100"),
]


class Loci(unittest.TestCase):
    def test_pat_intervals_one_per_locus(self):
        iv = L.intervals(ROWS, "pat")
        self.assertEqual([x[0] for x in iv], ["c01", "c03", "c07", "p12"])

    def test_mat_sites_merge_loci_that_lift_to_the_same_interval_and_skip_absent(self):
        iv = L.intervals(ROWS, "mat")
        self.assertEqual(sorted(x[0] for x in iv), ["c01|c03", "p12"])

    def test_locus_at_needs_half_the_span_inside(self):
        iv = L.intervals(ROWS, "pat")
        self.assertEqual(L.locus_at(iv, "P1", 120, 190), "c01")
        self.assertEqual(L.locus_at(iv, "P1", 190, 290), None)       # 10 of 100 bp inside
        self.assertEqual(L.locus_at(iv, "P9", 120, 190), None)       # wrong chromosome
        self.assertEqual(L.locus_at(iv, "P1", 150, 350), None)       # 50 of 200 bp inside c01, 50 inside p12: neither reaches half


class Groups(unittest.TestCase):
    def test_p12_and_p14_are_one_group_others_are_themselves(self):
        self.assertEqual(L.group("p12"), L.group("p14"))
        self.assertNotEqual(L.group("c01"), L.group("c03"))
        self.assertIsNone(L.group(None))


class Placement(unittest.TestCase):
    def test_best_identity_times_coverage_per_locus(self):
        iv = L.intervals(ROWS, "pat")
        recs = [
            dict(ref="P1", start=110, end=190, de=0.002, qcov=1.0),
            dict(ref="P1", start=120, end=180, de=0.020, qcov=0.9),   # second record on c01, worse
            dict(ref="P1", start=310, end=390, de=0.010, qcov=1.0),   # p12
            dict(ref="P5", start=0, end=100, de=0.0, qcov=1.0),       # elsewhere
        ]
        pl = L.placements(recs, iv)
        self.assertAlmostEqual(pl["c01"], 0.998)
        self.assertAlmostEqual(pl["p12"], 0.990)
        self.assertEqual(set(pl), {"c01", "p12"})

    def test_o3_flag_needs_the_reference_below_the_bar_and_the_other_haplotype_at_or_above_it(self):
        self.assertTrue(L.o3_flag(ref_best=0.987, other_best=0.999))
        self.assertFalse(L.o3_flag(ref_best=0.9995, other_best=1.0))   # the reference has it
        self.assertFalse(L.o3_flag(ref_best=0.5, other_best=0.99))     # the other does not have it either
        self.assertTrue(L.o3_flag(ref_best=None, other_best=0.999))     # nothing on the reference at all


if __name__ == "__main__":
    unittest.main()
