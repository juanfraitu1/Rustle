#!/usr/bin/env python3
"""Amendment 41: helpers of the maternal-assembly truth."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mat_truth as T


class CsDivergence(unittest.TestCase):
    def test_positions_on_the_target(self):
        # 10 matches, a substitution, 5 matches, a 2 bp deletion from the target, 3 matches, a 2 bp insertion, 2 matches
        ev, end = T.cs_divergence(100, ":10*ag:5-ac:3+tt:2")
        self.assertEqual(ev, [(110, 1), (116, 2), (121, 2)])
        self.assertEqual(end, 123)

    def test_long_runs(self):
        ev, end = T.cs_divergence(0, ":1000")
        self.assertEqual((ev, end), ([], 1000))


class Mask(unittest.TestCase):
    def test_uncovered_and_divergent_windows(self):
        covered = [(0, 1000)]
        events = [(10, 1), (20, 1), (30, 1), (600, 1)]          # 3 bases in window 0 (> 2.5), 1 in window 1
        self.assertEqual(T.mask(1500, covered, events, window=500, max_div=0.005), [(0, 500), (1000, 1500)])

    def test_merging_and_clipping(self):
        self.assertEqual(T.mask(1000, [(0, 400), (450, 1000)], [], window=500, max_div=0.005), [(400, 450)])
        self.assertEqual(T.mask(100, [], [], window=500, max_div=0.005), [(0, 100)])
        self.assertEqual(T.mask(100, [(0, 100)], [], window=500, max_div=0.005), [])

    def test_a_long_deletion_masks_its_window(self):
        self.assertEqual(T.mask(1000, [(0, 1000)], [(700, 40)], window=500, max_div=0.005), [(500, 1000)])


class AssemblyVersion(unittest.TestCase):
    def test_blocks_skip_introns_and_insertions(self):
        g = "AAAAACCCCCGGGGGTTTTT" * 5
        fetch = lambda a, b: g[a:b]
        # 5 aligned, 2 inserted read bases, 3 aligned, 1 deleted, 10 intron, 4 aligned, soft clip
        self.assertEqual(T.assembly_version(fetch, 0, "5M2I3M1D10N4M3S"), g[0:9] + g[19:23])
        self.assertEqual(T.assembly_version(fetch, 2, "4=1X"), g[2:7])


class Loci(unittest.TestCase):
    def test_grouping(self):
        rs = [("c1", 100, 500, "a"), ("c1", 1400, 1600, "b"), ("c1", 2700, 2800, "c"), ("c2", 100, 200, "d")]
        loci = T.group_loci(rs, gap=1000)
        self.assertEqual([(l["contig"], l["start"], l["end"], sorted(l["names"])) for l in loci],
                         [("c1", 100, 1600, ["a", "b"]), ("c1", 2700, 2800, ["c"]), ("c2", 100, 200, ["d"])])


class Testable(unittest.TestCase):
    def test_identical_chromosomes_are_never_tested(self):
        t = T.testable_contigs("mat")
        self.assertNotIn("CM054587.2", t)        # maternal chr5 = the primary's chr5
        self.assertIn("CM054594.2", t)           # maternal chr12, partner of the paternal primary chr12
        self.assertEqual(len([c for c in t if c.startswith("CM")]), 15)


def rec(de):
    return ["r", "1000", "0", "1000", "+", "c", "5000", "0", "1000", "990", "1000", "60", "tp:A:P", f"de:f:{de}"]


class Status(unittest.TestCase):
    def test_lower_divergence_decides(self):
        self.assertEqual(T.read_status("r", {"r": rec(0.001)}, {"r": True}, {"r": rec(0.004)}, {"r": False}), "mat")
        self.assertEqual(T.read_status("r", {"r": rec(0.004)}, {"r": True}, {"r": rec(0.001)}, {"r": False}), "shared")
        self.assertEqual(T.read_status("r", {}, {}, {"r": rec(0.002)}, {"r": True}), "pat")

    def test_tie_is_not_specific_and_unclean_is_unplaced(self):
        self.assertEqual(T.read_status("r", {"r": rec(0.002)}, {"r": True}, {"r": rec(0.002)}, {"r": False}), "shared")
        self.assertEqual(T.read_status("r", {"r": rec(0.02)}, {"r": True}, {"r": rec(0.03)}, {}), "unplaced")
        self.assertEqual(T.read_status("r", {}, {}, {}, {}), "unplaced")


class Verdict(unittest.TestCase):
    def test_majority_of_placed_reads(self):
        self.assertEqual(T.verdict(["mat", "mat", "shared"]), "TRUE")
        self.assertEqual(T.verdict(["pat", "shared", "unplaced", "unplaced"]), "WRONG")      # placed: 1 specific of 2, not most
        self.assertEqual(T.verdict(["pat", "mat", "shared", "unplaced"]), "TRUE")
        self.assertEqual(T.verdict(["unplaced", "unplaced"]), "WRONG")
        self.assertEqual(T.verdict(["shared"]), "WRONG")


if __name__ == "__main__":
    unittest.main()
