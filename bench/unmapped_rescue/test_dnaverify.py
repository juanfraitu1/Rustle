#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_dnaverify.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dna_verify as D  # noqa: E402


class Covered(unittest.TestCase):
    def test_a_base_is_covered_if_any_kmer_containing_it_is_present(self):
        # a 30-base sequence has 10 21-mers (starts 0..9); only the k-mer at start 5 is present: it covers bases 5..25
        counts = [0] * 10
        counts[5] = 3
        self.assertEqual(D.covered_bases(counts, 21), 21)

    def test_count_below_the_threshold_does_not_cover(self):
        counts = [1] * 10
        self.assertEqual(D.covered_bases(counts, 21, min_count=2), 0)
        self.assertEqual(D.covered_bases(counts, 21, min_count=1), 30)

    def test_a_junction_gap_leaves_the_flanks_covered(self):
        # 100-base sequence, 80 k-mers; the 20 junction k-mers (starts 30..49) are absent, the others present
        counts = [5] * 80
        for i in range(30, 50):
            counts[i] = 0
        self.assertEqual(D.covered_bases(counts, 21), 100)          # every base still lies in some present k-mer

    def test_empty(self):
        self.assertEqual(D.covered_bases([], 21), 0)


class Summary(unittest.TestCase):
    def test_fractions_and_median(self):
        counts = [0, 0, 4, 6, 8, 10]
        s = D.summarize(counts, 26)
        self.assertEqual(s["kmers"], 6)
        self.assertAlmostEqual(s["frac_ge2"], 4 / 6)
        self.assertEqual(s["median_ge2"], 7)


if __name__ == "__main__":
    unittest.main()
