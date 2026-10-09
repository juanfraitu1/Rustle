#!/usr/bin/env python3
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import species_check as S


def rec(q, qlen, matches, tp="P"):
    return "\t".join([q, str(qlen), "0", str(qlen), "+", "chr1", "1000", "0", "1000", str(matches), str(qlen), "60", f"tp:A:{tp}", "de:f:0.01"])


class PafScores(unittest.TestCase):
    def test_sums_primary_records_over_query_length(self):
        s = S.paf_scores([rec("a", 1000, 400), rec("a", 1000, 500), rec("b", 200, 100)])
        self.assertAlmostEqual(s["a"], 0.9)
        self.assertAlmostEqual(s["b"], 0.5)

    def test_ignores_secondary_records_and_caps_at_one(self):
        s = S.paf_scores([rec("a", 1000, 900), rec("a", 1000, 800, tp="S"), rec("c", 100, 90), rec("c", 100, 90)])
        self.assertAlmostEqual(s["a"], 0.9)
        self.assertEqual(s["c"], 1.0)

    def test_empty(self):
        self.assertEqual(S.paf_scores([]), {})


class Best(unittest.TestCase):
    def test_best_species_is_the_maximum_and_a_missing_score_is_zero(self):
        self.assertEqual(S.best_species({"HSA": None, "GGO": 0.97, "PTR": 0.995}), ("PTR", 0.995))
        self.assertEqual(S.best_species({"HSA": None, "GGO": None}), (None, 0.0))

    def test_same_species(self):
        self.assertTrue(S.same_species(("GGO", 0.991), "GGO"))
        self.assertFalse(S.same_species(("GGO", 0.98), "GGO"))
        self.assertFalse(S.same_species(("PTR", 0.999), "GGO"))
        self.assertFalse(S.same_species((None, 0.0), "GGO"))

    def test_read_share(self):
        rows = [dict(reads=30, own=True), dict(reads=10, own=False), dict(reads=10, own=True)]
        self.assertAlmostEqual(S.read_share(rows, lambda r: r["own"]), 0.8)
        self.assertIsNone(S.read_share([], lambda r: True))


if __name__ == "__main__":
    unittest.main()
