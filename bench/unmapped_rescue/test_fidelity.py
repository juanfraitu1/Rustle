#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_fidelity.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fidelity as F  # noqa: E402

ERASED = {"F1": ("chrA", 1000, 2000), "F2": ("chrA", 5000, 6000)}
TRUTH = {"a": "F1", "b": "F1", "c": "F1", "d": "F2", "e": "F2", "f": "F2"}


class Registered(unittest.TestCase):
    def test_registered_metric_is_identity_times_coverage(self):
        clusters = {"c1": ["a", "b", "c"], "c2": ["d", "e", "f"]}
        hits = {"c1": (0.9995, 0.9996, "chrA", 1200, 1900), "c2": (0.9995, 0.90, "chrA", 5100, 5900)}
        m = F.fidelity(clusters, TRUTH, hits, ERASED)
        self.assertEqual(m["on_erased_copy"], 2)
        self.assertEqual(m["identity_ge_0_999"], 2)
        self.assertEqual(m["identity_x_coverage_ge_0_999"], 1)          # c2 has coverage .90

    def test_restricted_to_attributed_clusters(self):
        clusters = {"c1": ["a", "b", "c"], "c2": ["d", "e", "f"]}
        hits = {"c1": (1.0, 1.0, "chrA", 1200, 1900), "c2": (1.0, 1.0, "chrA", 5100, 5900)}
        m = F.fidelity(clusters, TRUTH, hits, ERASED, only={"c2"})
        self.assertEqual((m["clusters_with_family"], m["identity_x_coverage_ge_0_999"]), (1, 1))

    def test_a_hit_elsewhere_is_not_on_the_erased_copy(self):
        m = F.fidelity({"c1": ["a", "b", "c"]}, TRUTH, {"c1": (1.0, 1.0, "chrB", 1200, 1900)}, ERASED)
        self.assertEqual((m["on_erased_copy"], m["identity_x_coverage_ge_0_999"]), (0, 0))


if __name__ == "__main__":
    unittest.main()
