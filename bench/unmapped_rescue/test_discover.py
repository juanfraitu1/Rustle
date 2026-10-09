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


class Rescore(unittest.TestCase):
    """Amendment 35: the saved scores are re-classified with the corrected rule; UNSUPPORTED is kept"""

    def test_rescore(self):
        rows = [dict(k="a", reads=10, cls="ELSEWHERE", scores=dict(R=0.988, GGO=0.965)),
                dict(k="b", reads=30, cls="ELSEWHERE", scores=dict(R=None, GGO=0.99)),
                dict(k="c", reads=5, cls="UNSUPPORTED", scores=dict(R=0.5, GGO=0.99)),
                dict(k="d", reads=4, cls="NOVEL", scores=dict(R=0.3, GGO=0.4))]
        new = D.rescore_rows(rows)
        self.assertEqual([r["cls"] for r in new], ["DIVERGED", "ELSEWHERE", "UNSUPPORTED", "NOVEL"])
        self.assertEqual(rows[0]["cls"], "ELSEWHERE")                              # the input is not modified
        self.assertEqual([r["k"] for r in D.changed(rows, new)], ["a"])

    def test_rescore_cross_species_scores_count_only_when_the_reference_lacks_the_sequence(self):
        rows = [dict(k="a", reads=5, cls="ELSEWHERE", scores=dict(R=0.97, HSA=0.994, GGO=0.994)),
                dict(k="b", reads=5, cls="ELSEWHERE", scores=dict(R=None, HSA=0.5, GGO=0.99))]
        self.assertEqual([r["cls"] for r in D.rescore_rows(rows, cross=("HSA", "GGO"))], ["DIVERGED", "ELSEWHERE"])
        self.assertEqual([r["cls"] for r in D.rescore_rows(rows)], ["ELSEWHERE", "ELSEWHERE"])      # as same-species assemblies the 35b rule keeps the first

    def test_every_animal_has_a_cross_entry_naming_real_assemblies(self):
        for name, a in D.ANIMALS.items():
            self.assertTrue(set(D.CROSS[name]) <= set(a["others"]), name)

    def test_rescore_legacy_rows_with_haplotype_assemblies_are_not_touched(self):
        rows = [dict(k="a", reads=3, cls="CONFIRMED", scores=dict(R=0.5, Tm=1.0, Tp=0.2))]
        self.assertEqual(D.rescore_rows(rows, animal=None)[0]["cls"], "CONFIRMED")


class Config(unittest.TestCase):
    def test_every_animal_has_a_reference_and_a_bam(self):
        for name, a in D.ANIMALS.items():
            self.assertIn("bam", a, name)
            self.assertIn("R", a, name)
            self.assertTrue(a["others"], name)
            self.assertNotIn("R", a["others"], name)


if __name__ == "__main__":
    unittest.main()
