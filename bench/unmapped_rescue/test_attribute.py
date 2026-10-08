#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_attribute.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import attribute as A  # noqa: E402
import consensus as C  # noqa: E402


def row(q, t, bits, ev=1e-20):
    return (q, t, float(bits), float(ev))


class Family(unittest.TestCase):
    def test_family_of_a_target_is_the_prefix_before_the_colon(self):
        self.assertEqual(A.family_of("GWFAM12:3"), "GWFAM12")

    def test_best_family_wins_when_strictly_above_every_other(self):
        rows = [row("c1", "F1:0", 100), row("c1", "F1:1", 90), row("c1", "F2:0", 80)]
        self.assertEqual(A.attribute(rows)["c1"], ("F1", 100.0, 80.0))

    def test_tie_between_families_abstains(self):
        rows = [row("c1", "F1:0", 100), row("c1", "F2:0", 100)]
        self.assertEqual(A.attribute(rows)["c1"], (None, 100.0, 100.0))

    def test_evalue_cut_and_a_query_with_no_hit(self):
        rows = [row("c1", "F1:0", 100, 1e-2), row("c2", "F1:0", 50, 1e-9)]
        out = A.attribute(rows, max_evalue=1e-5)
        self.assertNotIn("c1", out)
        self.assertEqual(out["c2"][0], "F1")

    def test_single_family_hit_has_zero_runner_up(self):
        self.assertEqual(A.attribute([row("c", "F1:0", 70)])["c"], ("F1", 70.0, 0.0))


class Cover(unittest.TestCase):
    def test_union_length_merges_overlaps_and_counts_inclusive_bases(self):
        self.assertEqual(A.union_length([(1, 10), (5, 20), (30, 39)]), 30)
        self.assertEqual(A.union_length([]), 0)

    def test_family_score_is_the_covered_bases_of_its_best_target(self):
        hsps = [("c1", "F1:0", 1, 100), ("c1", "F1:0", 50, 150), ("c1", "F1:1", 1, 80), ("c1", "F2:0", 200, 260)]
        self.assertEqual(A.cover_scores(hsps), {"c1": {"F1": 150, "F2": 61}})

    def test_strand_flipped_hsp_coordinates_are_normalised(self):
        self.assertEqual(A.cover_scores([("c1", "F1:0", 100, 1)]), {"c1": {"F1": 100}})

    def test_margin_rule(self):
        sc = {"win": {"F1": 110, "F2": 100}, "tie": {"F1": 100, "F2": 100}, "close": {"F1": 109, "F2": 100}, "alone": {"F1": 5}, "none": {}}
        out = A.attribute_cover(sc, margin=1.10)
        self.assertEqual(out["win"], ("F1", 110, 100))
        self.assertEqual(out["tie"][0], None)
        self.assertEqual(out["close"][0], None)
        self.assertEqual(out["alone"], ("F1", 5, 0))
        self.assertEqual(out["none"], (None, 0, 0))


class Pick(unittest.TestCase):
    def test_longest_first_ties_by_name_and_capped(self):
        lens = {"b": 10, "a": 10, "c": 30, "d": 5}
        self.assertEqual(C.pick_reads(["a", "b", "c", "d"], lens, 3), ["c", "a", "b"])


if __name__ == "__main__":
    unittest.main()
