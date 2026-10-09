#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_rescue.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flagmetric as F  # noqa: E402
import rescue as R  # noqa: E402


def paf(q, qlen, qs, qe, matches, blk):
    return "\t".join([q, str(qlen), str(qs), str(qe), "+", "w", "9999", "100", "200", str(matches), str(blk), "60"])


class EndsAndUnion(unittest.TestCase):
    def test_union_length(self):
        self.assertEqual(R.union_length([(0, 10), (5, 20), (30, 40)]), 30)
        self.assertEqual(R.union_length([]), 0)

    def test_unaligned_end_segments_need_the_minimum_length(self):
        self.assertEqual(R.end_segments(137, 5694, 5694), [(0, 137)])
        self.assertEqual(R.end_segments(10, 5000, 5100), [(5000, 5100)])           # a 10-base 5' end is too short, the 100-base 3' end counts
        self.assertEqual(R.end_segments(10, 5000, 5100, min_len=200), [])
        self.assertEqual(R.end_segments(0, 100, 100), [])


class Pieces(unittest.TestCase):
    def test_pieces_with_enough_identity_become_consensus_intervals(self):
        segs = [(0, 137)]
        lines = [paf("seg0", 137, 47, 137, 90, 90), paf("seg0", 137, 2, 48, 46, 46), paf("seg0", 137, 60, 100, 30, 40)]   # the third is 75% identical: dropped
        iv, m, b = R.rescued(segs, lines)
        self.assertEqual(sorted(iv), [(2, 48), (47, 137)])
        self.assertEqual((m, b), (136, 136))

    def test_segment_offsets_are_added(self):
        iv, m, b = R.rescued([(5000, 5100)], [paf("seg0", 100, 10, 90, 80, 80)])
        self.assertEqual(iv, [(5010, 5090)])

    def test_a_piece_overlapping_the_main_alignment_does_not_double_count(self):
        total, ident = R.combine(qlen=1000, main=(100, 1000, 899, 900), pieces=([(40, 130)], 88, 90))
        self.assertEqual(total, 960)                       # 40..1000, the 30 bases 100-130 counted once
        self.assertAlmostEqual(ident, (899 + 88) / (900 + 90))


class GateAwareTotal(unittest.TestCase):
    def test_leading_pure_g_excluded_when_the_gate_is_open(self):
        self.assertAlmostEqual(F.gate_aware_total(1.0, 998, 1000, 2, "GG", True), 1.0)
        self.assertAlmostEqual(F.gate_aware_total(1.0, 998, 1000, 2, "GG", False), 0.998)
        self.assertAlmostEqual(F.gate_aware_total(1.0, 960, 1000, 40, "G" * 40, True), 0.96)     # longer than 3: untouched
        self.assertAlmostEqual(F.gate_aware_total(0.99, 1000, 1000, 0, "", True), 0.99)


if __name__ == "__main__":
    unittest.main()
