#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_score.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as S  # noqa: E402

TRUTH = {"a1": "A", "a2": "A", "a3": "A", "a4": "A", "b1": "B", "b2": "B", "b3": "B", "x1": "A", "bg1": "bg", "u1": None}


class Cluster(unittest.TestCase):
    def test_coverage_purity_fragmentation_and_background(self):
        clusters = {1: ["a1", "a2", "a3", "x1"], 2: ["b1", "b2", "b3", "a4"], 3: ["bg1", "u1", "a4x"]}
        m = S.cluster_metrics(clusters, TRUTH, d_reads={"a1", "a2", "a3", "a4", "b1", "b2", "b3", "x1"})
        self.assertEqual(m["n_d"], 8)
        self.assertEqual(m["d_clustered"], 8)                   # a4x is unknown to the truth and ignored
        self.assertAlmostEqual(m["coverage"], 1.0)
        self.assertAlmostEqual(m["purity"], 7 / 8)              # a4 sits in cluster 2 whose majority is B
        self.assertEqual(m["clusters_per_family"], {"A": 1, "B": 1})
        self.assertEqual(m["background_in_clusters"], 1)

    def test_no_clusters(self):
        m = S.cluster_metrics({}, TRUTH, d_reads={"a1"})
        self.assertEqual((m["coverage"], m["purity"]), (0.0, None))


class Rescue(unittest.TestCase):
    def test_rescued_correct_wrong_joins_and_abstentions(self):
        truth = {"a1": "A", "a2": "A", "a3": "A", "b1": "B", "b2": "B", "b3": "B", "bg1": "bg", "a4": "A", "a5": "A", "a6": "A", "u1": None}
        clusters = {1: ["a1", "a2", "a3"], 2: ["b1", "b2", "b3", "bg1", "u1"], 3: ["a4", "a5", "a6"]}
        d = {r for r, f in truth.items() if f in ("A", "B")}
        m = S.rescue_metrics(clusters, {1: "A", 2: "A", 3: None}, truth, d)
        self.assertEqual(m["rescued_correct"], 3)
        self.assertEqual(m["wrong_joins"], 4)                    # b1-b3 joined to A, plus the background read
        self.assertEqual(m["unknown_joined"], 1)                 # u1 has no truth label: reported, not judged
        self.assertEqual(m["abstained_reads"], 3)
        self.assertEqual((m["clusters_attributed"], m["clusters_abstained"]), (2, 1))
        self.assertAlmostEqual(m["cluster_accuracy"], 0.5)
        self.assertEqual(m["copies_reached"], 1)
        self.assertEqual(m["copies_with_unmapped_reads"], 2)

    def test_nothing_attributed(self):
        m = S.rescue_metrics({1: ["a1", "a2", "a3"]}, {1: None}, {"a1": "A", "a2": "A", "a3": "A"}, {"a1", "a2", "a3"})
        self.assertEqual((m["rescued_correct"], m["cluster_accuracy"], m["copies_reached"]), (0, None, 0))


if __name__ == "__main__":
    unittest.main()
