#!/usr/bin/env python3
"""Amendment 42: locus-binned clustering."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import binclust as B


class Bins(unittest.TestCase):
    def test_overlapping_spans_same_strand_share_a_bin(self):
        recs = [("a", "c1", 100, 500, "+"), ("b", "c1", 400, 900, "+"), ("c", "c1", 950, 1200, "+"),
                ("d", "c1", 450, 800, "-"), ("e", "c2", 100, 500, "+")]
        self.assertEqual(sorted(sorted(b) for b in B.bin_reads(recs)), [["a", "b"], ["c"], ["d"], ["e"]])

    def test_touching_spans_do_not_join(self):
        self.assertEqual(len(B.bin_reads([("a", "c", 0, 100, "+"), ("b", "c", 100, 200, "+")])), 2)

    def test_chained_overlaps_join(self):
        recs = [("a", "c", 0, 100, "+"), ("b", "c", 90, 300, "+"), ("c", "c", 250, 400, "+")]
        self.assertEqual(B.bin_reads(recs), [["a", "b", "c"]])


def pairs(*ps):
    return {tuple(sorted(p)) for p in ps}


class ClusterBin(unittest.TestCase):
    def test_small_bin_is_the_star_step(self):
        lens = {"a": 10, "b": 9, "c": 8, "d": 7, "e": 6}
        compat = pairs(("a", "b"), ("a", "c"), ("d", "e"))
        cl = B.cluster_bin(list(lens), lens, lambda head: compat, lambda rest, centers: {}, cap=500)
        self.assertEqual(cl, {"a": ["a", "b", "c"]})               # d-e is a cluster of 2: dropped

    def test_reads_beyond_the_cap_join_a_centre(self):
        lens = {f"r{i}": 100 - i for i in range(8)}
        compat = pairs(("r0", "r1"), ("r0", "r2"))
        assign = lambda rest, centers: {r: "r0" for r in rest if r in ("r5", "r6")}
        cl = B.cluster_bin(list(lens), lens, lambda head: {p for p in compat if p[0] in head and p[1] in head}, assign, cap=3)
        self.assertEqual(sorted(cl["r0"]), ["r0", "r1", "r2", "r5", "r6"])

    def test_leftovers_get_their_own_pass(self):
        lens = {f"r{i}": 100 - i for i in range(9)}
        compat = pairs(("r0", "r1"), ("r0", "r2"), ("r3", "r4"), ("r3", "r5"))
        calls = []

        def cf(head):
            calls.append(list(head))
            return {p for p in compat if p[0] in head and p[1] in head}
        cl = B.cluster_bin(list(lens), lens, cf, lambda rest, centers: {}, cap=3)
        self.assertEqual(sorted(cl), ["r0", "r3"])
        self.assertEqual(sorted(cl["r3"]), ["r3", "r4", "r5"])
        self.assertEqual(calls[0], ["r0", "r1", "r2"])

    def test_no_progress_stops(self):
        lens = {f"r{i}": 100 - i for i in range(6)}
        n = []
        cl = B.cluster_bin(list(lens), lens, lambda head: n.append(1) or set(), lambda rest, centers: {}, cap=3)
        self.assertEqual(cl, {})
        self.assertEqual(len(n), 1)

    def test_interrupted_run_resumes_to_the_same_result(self):
        lens = {f"r{i}": 100 - i for i in range(9)}
        compat = pairs(("r0", "r1"), ("r0", "r2"), ("r3", "r4"), ("r3", "r5"), ("r6", "r7"), ("r6", "r8"))
        cf = lambda head: {p for p in compat if p[0] in head and p[1] in head}
        full = B.cluster_bin(list(lens), lens, cf, lambda rest, centers: {}, cap=3)
        state = {}
        part = B.cluster_bin(list(lens), lens, cf, lambda rest, centers: {}, cap=3, state=state, stop=lambda: True)
        self.assertIsNone(part)                      # stopped after one pass: no result yet
        self.assertEqual(sorted(state["clusters"]), ["r0"])
        resumed = B.cluster_bin(list(lens), lens, cf, lambda rest, centers: {}, cap=3, state=state)
        self.assertEqual(resumed, full)

    def test_first_centre_in_order_wins(self):
        self.assertEqual(B.pick_centres({"x": {"c2", "c1"}, "y": {"c2"}}, ["c1", "c2"]), {"x": "c1", "y": "c2"})


class Minimap(unittest.TestCase):
    """integration: the minimap2-backed helpers on reads with known structure"""

    def test_assign_and_compat(self):
        import random
        import shutil
        if not shutil.which("minimap2"):
            self.skipTest("minimap2")
        rnd = random.Random(3)
        t = "".join(rnd.choice("ACGT") for _ in range(2000))
        other = "".join(rnd.choice("ACGT") for _ in range(2000))
        seqs = {"c": t, "a": t[100:1900], "b": t[300:2000], "x": other[0:1800]}
        a = B.assign_fn(seqs, "/tmp/binclust_test")(["a", "b", "x"], ["c"])
        self.assertEqual(a, {"a": "c", "b": "c"})
        comp = B.compat_fn(seqs, "/tmp/binclust_test")(["c", "a", "b", "x"])
        self.assertIn(("a", "c"), comp)
        self.assertFalse(any("x" in p for p in comp))


if __name__ == "__main__":
    unittest.main()
