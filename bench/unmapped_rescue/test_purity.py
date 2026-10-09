#!/usr/bin/env python3
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import control_purity as P


class Loci(unittest.TestCase):
    def test_is_multi(self):
        self.assertFalse(P.is_multi([("c1", 10), ("c1", 5000)]))
        self.assertTrue(P.is_multi([("c1", 10), ("c1", 500_000)]))
        self.assertTrue(P.is_multi([("c1", 1), ("c2", 1)]))
        self.assertFalse(P.is_multi([("c1", 7)]))
        self.assertFalse(P.is_multi([("c1", 0), ("c1", 100_000)]))        # span exactly 100 kb is still one locus
        self.assertTrue(P.is_multi([("c1", 0), ("c1", 100_001)]))

    def test_locus_count_is_single_linkage(self):
        self.assertEqual(P.locus_count([("c1", 0), ("c1", 90_000), ("c1", 180_000)]), 1)   # a chain of 90 kb gaps is one locus
        self.assertEqual(P.locus_count([("c1", 0), ("c1", 500_000), ("c2", 5), ("c1", 500_100)]), 3)
        self.assertEqual(P.locus_count([]), 0)

    def test_share(self):
        rows = [dict(cls="NOVEL", reads=4, multi=True), dict(cls="DIVERGED", reads=6, multi=False), dict(cls="PRESENT", reads=10, multi=True)]
        self.assertAlmostEqual(P.multi_share(rows, ("NOVEL", "DIVERGED")), 0.4)
        self.assertAlmostEqual(P.multi_share(rows, ("PRESENT",)), 1.0)
        self.assertIsNone(P.multi_share(rows, ("UNSUPPORTED",)))


if __name__ == "__main__":
    unittest.main()
