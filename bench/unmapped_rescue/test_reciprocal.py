#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_reciprocal.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reciprocal as RC  # noqa: E402


class Segments(unittest.TestCase):
    def test_exons_joined_introns_dropped_deletions_kept_insertions_ignored(self):
        # target starts at 100: 50 matched, intron of 1000, 30 matched, 5 deleted from the query (kept: they are reference bases), 3 inserted in the query, 20 matched
        segs = RC.target_segments("50=1000N30=5D3I20=", 100)
        self.assertEqual(segs, [(100, 150), (1150, 1205)])

    def test_adjacent_blocks_are_merged(self):
        self.assertEqual(RC.target_segments("10M2X8M", 0), [(0, 20)])

    def test_mismatch_and_match_ops_both_advance_the_target(self):
        self.assertEqual(RC.target_segments("5=1X5=", 7), [(7, 18)])


class Overlap(unittest.TestCase):
    def test_same_chromosome_and_half_of_the_shorter(self):
        self.assertTrue(RC.same_locus(("c1", 100, 200), ("c1", 150, 400)))      # 50 of the shorter 100
        self.assertFalse(RC.same_locus(("c1", 100, 200), ("c1", 160, 400)))     # 40 of 100
        self.assertFalse(RC.same_locus(("c1", 100, 200), ("c2", 100, 200)))
        self.assertFalse(RC.same_locus(None, ("c1", 100, 200)))


class Transcript(unittest.TestCase):
    def test_minus_strand_transcripts_are_reverse_complemented(self):
        self.assertEqual(RC.orient("AACG", "+"), "AACG")
        self.assertEqual(RC.orient("AACG", "-"), "CGTT")


if __name__ == "__main__":
    unittest.main()
