#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_graph.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import graph as G  # noqa: E402


def paf(q, t, ql, tl, alnlen, de, qs=0, qe=None):
    qe = qe if qe is not None else ql
    return "\t".join([q, str(ql), str(qs), str(qe), "+", t, str(tl), "0", str(tl), str(alnlen - 5), str(alnlen), "60", f"de:f:{de}"])


class Edges(unittest.TestCase):
    def test_edge_needs_low_divergence_and_long_enough_overlap(self):
        lines = [paf("a", "b", 1000, 1200, 900, 0.004),      # ok: block .90 of the shorter read
                 paf("a", "c", 1000, 1200, 900, 0.02),       # too divergent
                 paf("a", "d", 1000, 1200, 400, 0.001),      # block .40 of the shorter read
                 paf("a", "e", 1000, 1200, 500, 0.001)]      # block .50: boundary is an edge
        self.assertEqual(sorted(G.edges(lines, delta=0.00958, min_frac=0.5)), [("a", "b"), ("a", "e")])

    def test_self_hits_and_missing_de_are_dropped(self):
        lines = [paf("a", "a", 1000, 1000, 900, 0.0), "\t".join(["a", "1000", "0", "900", "+", "b", "1000", "0", "900", "890", "900", "60"])]
        self.assertEqual(list(G.edges(lines, 0.00958, 0.5)), [])

    def test_edges_are_unordered_and_unique(self):
        lines = [paf("b", "a", 1000, 1000, 900, 0.001), paf("a", "b", 1000, 1000, 900, 0.001)]
        self.assertEqual(sorted(set(G.edges(lines, 0.00958, 0.5))), [("a", "b")])


class Components(unittest.TestCase):
    def test_components_and_the_size_floor(self):
        comp = G.components([("a", "b"), ("b", "c"), ("x", "y")], nodes=["a", "b", "c", "x", "y", "z"])
        self.assertEqual(comp["a"], comp["c"])
        self.assertNotEqual(comp["a"], comp["x"])
        self.assertEqual(sorted(len(v) for v in G.clusters(comp, min_size=3).values()), [3])

    def test_singletons_get_their_own_component(self):
        comp = G.components([], nodes=["a", "b"])
        self.assertNotEqual(comp["a"], comp["b"])


if __name__ == "__main__":
    unittest.main()
