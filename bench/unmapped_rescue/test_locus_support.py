#!/usr/bin/env python3
"""Amendment 44: read support of a flag among all reads of its locus."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import locus_support as L


def rec(q, t, matches, tp="P"):
    return "\t".join([q, "1000", "0", "1000", "+", t, "5000", "0", "1000", str(matches), "1000", "60", f"tp:A:{tp}"])


class Assign(unittest.TestCase):
    def test_primary_record_decides(self):
        lines = [rec("r1", "cons", 990), rec("r2", "ref", 995), rec("r3", "cons", 900, "S"), rec("r3", "ref", 980), rec("r4", "other", 999)]
        self.assertEqual(L.assign(lines, "cons", "ref"), (1, 2))

    def test_empty(self):
        self.assertEqual(L.assign([], "cons", "ref"), (0, 0))


class Stats(unittest.TestCase):
    def test_share(self):
        self.assertAlmostEqual(L.share(3, 9), 0.25)
        self.assertIsNone(L.share(0, 0))

    def test_minor_is_a_one_sided_binomial_against_half(self):
        self.assertTrue(L.minor(2, 100))             # 2 of 100: far below an allele's half
        self.assertFalse(L.minor(45, 100))
        self.assertFalse(L.minor(3, 5))
        self.assertFalse(L.minor(0, 0))              # no reads: no evidence
        self.assertFalse(L.minor(1, 6))              # P(X <= 1 | 6, 1/2) = 0.109 > 0.01
        self.assertTrue(L.minor(0, 8))               # P = 0.0039

    def test_sample_every(self):
        self.assertEqual(L.sample_every(list(range(10)), 5), [0, 2, 4, 6, 8])
        self.assertEqual(L.sample_every(list(range(3)), 5), [0, 1, 2])


if __name__ == "__main__":
    unittest.main()
