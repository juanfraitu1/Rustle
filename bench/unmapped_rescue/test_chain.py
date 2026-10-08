#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_chain.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chain as C  # noqa: E402


def key(a, b):
    return (a, b) if a < b else (b, a)


class Star(unittest.TestCase):
    def test_a_contained_bridge_goes_to_the_longest_representative_only(self):
        lens = {"A": 1000, "A1": 950, "A2": 940, "B": 300, "C": 900, "C1": 880, "C2": 870}
        compat = {key("A", x) for x in ("A1", "A2", "B")} | {key("A1", "A2"), key("A1", "B"), key("A2", "B")}
        compat |= {key("C", x) for x in ("C1", "C2", "B")} | {key("C1", "C2"), key("C1", "B"), key("C2", "B")}
        out = C.star_clusters(lens, compat, min_size=3)
        self.assertEqual(sorted(sorted(c) for c in out), [["A", "A1", "A2", "B"], ["C", "C1", "C2"]])

    def test_a_chain_of_three_reads_is_dropped(self):
        lens = {"A": 1000, "B": 300, "C": 900}
        out = C.star_clusters(lens, {key("A", "B"), key("B", "C")}, min_size=3)
        self.assertEqual(out, [])

    def test_a_compatible_set_stays_one_cluster(self):
        lens = {r: 1000 - i for i, r in enumerate("abcde")}
        compat = {key(x, y) for x in "abcde" for y in "abcde" if x < y}
        self.assertEqual([sorted(c) for c in C.star_clusters(lens, compat, 3)], [list("abcde")])

    def test_ties_in_length_break_by_name(self):
        lens = {"b": 500, "a": 500, "c": 500}
        out = C.star_clusters(lens, {key("a", "b"), key("a", "c"), key("b", "c")}, 3)
        self.assertEqual(out[0][0], "a")


class Refine(unittest.TestCase):
    def test_small_components_are_refined_big_ones_kept(self):
        calls = []

        def allvsall(reads):
            calls.append(sorted(reads))
            return []                       # no compatible pair at all
        clusters = {"x": ["r1", "r2", "r3"], "y": [f"s{i}" for i in range(70)]}
        seqs = {r: "A" * 100 for rs in clusters.values() for r in rs}
        out = C.refine(clusters, seqs, allvsall, delta=0.00958, min_size=3, max_component=60)
        self.assertEqual(sorted(out), ["y"])            # x had no compatible pair and is dropped; y kept as it was
        self.assertEqual(calls, [["r1", "r2", "r3"]])

    def test_a_refined_cluster_is_named_after_its_representative(self):
        lines = []

        def paf(q, t):
            return "\t".join([q, "1000", "0", "1000", "+", t, "1000", "0", "1000", "995", "1000", "60", "NM:i:5", "de:f:0.001"])
        for a, b in (("r1", "r2"), ("r1", "r3"), ("r2", "r3")):
            lines += [paf(a, b), paf(b, a)]
        seqs = {"r1": "A" * 1000, "r2": "A" * 990, "r3": "A" * 980}
        out = C.refine({"x": ["r1", "r2", "r3"]}, seqs, lambda reads: lines, 0.00958, 3, 60)
        self.assertEqual(out, {"r1": ["r1", "r2", "r3"]})


if __name__ == "__main__":
    unittest.main()
