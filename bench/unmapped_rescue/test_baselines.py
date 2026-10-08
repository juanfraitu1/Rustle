#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_baselines.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import baselines as B  # noqa: E402


def paf(q, ql, qs, qe, t, matches, de):
    return "\t".join([q, str(ql), str(qs), str(qe), "+", t, "9999", "0", "999", str(matches), str(qe - qs), "60", f"de:f:{de}"])


class B0(unittest.TestCase):
    def test_best_hit_with_coverage_and_divergence_floors(self):
        lines = [paf("r1", 1000, 0, 800, "F1:0", 600, 0.15),       # coverage .8, de .15: usable
                 paf("r1", 1000, 0, 900, "F2:0", 700, 0.25),       # too divergent
                 paf("r2", 1000, 0, 400, "F1:0", 300, 0.10),       # coverage .4: too short
                 paf("r3", 1000, 0, 600, "F3:0", 400, 0.10),
                 paf("r3", 1000, 0, 600, "F4:0", 400, 0.10)]       # two families tie
        out = B.single_read_nucleotide(lines)
        self.assertEqual(out, {"r1": "F1", "r3": None})

    def test_unmapped_read_is_absent(self):
        self.assertEqual(B.single_read_nucleotide([]), {})


class Sample(unittest.TestCase):
    def test_seeded_sample_is_stable_and_capped(self):
        reads = [f"r{i}" for i in range(100)]
        self.assertEqual(B.sample(reads, 10, seed=1), B.sample(reads, 10, seed=1))
        self.assertEqual(len(B.sample(reads, 10, seed=1)), 10)
        self.assertEqual(sorted(B.sample(["a", "b"], 10, seed=1)), ["a", "b"])


if __name__ == "__main__":
    unittest.main()
