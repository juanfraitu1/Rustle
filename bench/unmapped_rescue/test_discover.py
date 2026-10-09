#!/usr/bin/env python3
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discover as D


class DenseControl(unittest.TestCase):
    """Amendment 33 (control): reads from loci that hold >= min_reads reads in one bin, so the control can form clusters"""

    def test_keeps_only_dense_bins(self):
        recs = [("r%d" % i, "chr1", 10_000 + 100 * i) for i in range(6)] + [("s0", "chr1", 90_000), ("t0", "chr2", 10_100)]
        got = D.dense_reads(recs, bin_size=5000, min_reads=5)
        self.assertEqual(sorted(got), ["r%d" % i for i in range(6)])

    def test_bins_do_not_mix_contigs(self):
        recs = [("a%d" % i, "chr1", 1000 + i) for i in range(3)] + [("b%d" % i, "chr2", 1000 + i) for i in range(3)]
        self.assertEqual(D.dense_reads(recs, bin_size=5000, min_reads=5), [])

    def test_sample_is_seeded_and_capped(self):
        names = ["r%d" % i for i in range(100)]
        a, b = D.sample_names(names, 10, seed=5), D.sample_names(names, 10, seed=5)
        self.assertEqual(a, b)
        self.assertEqual(len(a), 10)
        self.assertNotEqual(a, D.sample_names(names, 10, seed=6))
        self.assertEqual(D.sample_names(names[:4], 10, seed=5), sorted(names[:4]))


class Counts(unittest.TestCase):
    def test_class_counts(self):
        rows = [dict(cls="PRESENT", reads=5), dict(cls="PRESENT", reads=3), dict(cls="NOVEL", reads=7)]
        c = D.class_counts(rows)
        self.assertEqual(c["PRESENT"], [2, 8])
        self.assertEqual(c["NOVEL"], [1, 7])
        self.assertEqual(c["ELSEWHERE"], [0, 0])

    def test_specific_bar(self):
        rows = [dict(cls="PRESENT", reads=90), dict(cls="ALLELE-LIKE", reads=6), dict(cls="NOVEL", reads=4)]
        self.assertAlmostEqual(D.specific_fraction(rows), 0.96)
        self.assertEqual(D.specific_fraction([]), None)


class Config(unittest.TestCase):
    def test_every_animal_has_a_reference_and_a_bam(self):
        for name, a in D.ANIMALS.items():
            self.assertIn("bam", a, name)
            self.assertIn("R", a, name)
            self.assertTrue(a["others"], name)
            self.assertNotIn("R", a["others"], name)


if __name__ == "__main__":
    unittest.main()
