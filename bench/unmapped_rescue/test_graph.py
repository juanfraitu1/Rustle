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


def paf2(q, t, ql, tl, qs, qe, ts, te, de=0.001, strand="+"):
    blk = max(qe - qs, te - ts)
    return "\t".join([q, str(ql), str(qs), str(qe), strand, t, str(tl), str(ts), str(te), str(blk - 5), str(blk), "60", f"de:f:{de}"])


class ProperOverlap(unittest.TestCase):
    def kept(self, line):
        return len(list(G.edges([line], 0.00958, 0.5, proper=True))) == 1

    def test_containment_and_suffix_prefix_pass(self):
        self.assertTrue(self.kept(paf2("a", "b", 600, 3000, 0, 600, 1000, 1600)))          # a inside b
        self.assertTrue(self.kept(paf2("a", "b", 3000, 3000, 0, 3000, 0, 3000)))           # identical
        self.assertTrue(self.kept(paf2("a", "b", 3000, 3000, 1000, 3000, 0, 2000)))        # a's tail on b's head (dovetail)

    def test_a_shared_middle_is_not_an_edge(self):
        # both reads continue differently on both sides: 1500 of 3000 each
        self.assertFalse(self.kept(paf2("a", "b", 3000, 3000, 800, 2300, 700, 2200)))

    def test_ragged_ends_within_the_tolerance_pass(self):
        # remainders of 20 bases on both reads at the right end, block 2900 -> tolerance 27.8
        self.assertTrue(self.kept(paf2("a", "b", 3000, 3000, 0, 2980, 0, 2980)))
        self.assertTrue(self.kept(paf2("a", "b", 3000, 3000, 0, 2980, 25, 3005 - 5)))      # shifted start by 25: a reaches its start

    def test_one_end_interior_on_both_reads_fails(self):
        self.assertFalse(self.kept(paf2("a", "b", 3000, 3000, 0, 2000, 0, 2000)))          # right end: 1000 left on each

    def test_reverse_strand_pairs_the_query_start_with_the_target_end(self):
        # reverse strand: the query start pairs with the target END and the query end with the target START
        self.assertTrue(self.kept(paf2("a", "b", 3000, 3000, 0, 2000, 0, 2000, strand="-")))     # a's head at b's tail (b reaches 0 on one end), proper
        self.assertFalse(self.kept(paf2("a", "b", 3000, 3000, 500, 2500, 500, 2500, strand="-")))  # both reads continue at both ends

    def test_without_the_flag_the_frozen_rule_is_unchanged(self):
        line = paf2("a", "b", 3000, 3000, 800, 2300, 700, 2200)
        self.assertEqual(len(list(G.edges([line], 0.00958, 0.5))), 1)


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
