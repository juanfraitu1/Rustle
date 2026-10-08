#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_hybrid.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hybrid as H  # noqa: E402


class Combine(unittest.TestCase):
    def test_pooled_attribution_wins_and_the_fallback_fills_the_rest(self):
        clusters = {"c1": ["a", "b", "c"], "c2": ["d", "e", "f"]}
        att = H.combine(clusters, {"c1": "F1", "c2": None}, {"d": "F2", "e": None, "f": "F9", "x": "F3"}, reads=["a", "b", "c", "d", "e", "f", "x", "y"])
        self.assertEqual(att, {"a": "F1", "b": "F1", "c": "F1", "d": "F2", "e": None, "f": "F9", "x": "F3", "y": None})

    def test_a_pooled_read_ignores_the_fallback(self):
        att = H.combine({"c1": ["a", "b", "c"]}, {"c1": "F1"}, {"a": "F7"}, reads=["a", "b", "c"])
        self.assertEqual(att["a"], "F1")

    def test_residual_reads_are_the_ones_without_a_pooled_family(self):
        clusters = {"c1": ["a", "b", "c"], "c2": ["d", "e", "f"]}
        self.assertEqual(H.residual(clusters, {"c1": "F1", "c2": None}, ["a", "b", "c", "d", "e", "f", "x"]), ["d", "e", "f", "x"])


class Metrics(unittest.TestCase):
    def test_read_level_counts(self):
        truth = {"a": "F1", "b": "F1", "c": "F2", "bg1": "bg", "bg2": "bg"}
        att = {"a": "F1", "b": "F9", "c": None, "bg1": "F1", "bg2": None}
        m = H.read_metrics(att, truth, d_reads={"a", "b", "c"})
        self.assertEqual((m["rescued_correct"], m["wrong"], m["unattributed_d"]), (1, 2, 1))     # b wrong, bg1 wrong
        self.assertAlmostEqual(m["wrong_fraction_of_joined"], 2 / 3)


if __name__ == "__main__":
    unittest.main()
