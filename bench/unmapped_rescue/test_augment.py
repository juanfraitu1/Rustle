#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_augment.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import augment as G  # noqa: E402


def rec(score, de=0.01, qcov=1.0, ref="X", mapq=60):
    return dict(score=score, de=de, qcov=qcov, ref=ref, mapq=mapq)


class Primary(unittest.TestCase):
    def test_consensus_wins_only_with_a_strictly_higher_score(self):
        self.assertEqual(G.new_primary(rec(100), rec(101))[0], "consensus")
        self.assertEqual(G.new_primary(rec(100), rec(100))[0], "genome")
        self.assertEqual(G.new_primary(rec(100), rec(99))[0], "genome")

    def test_unmapped_read_takes_a_consensus_record_with_enough_coverage(self):
        self.assertEqual(G.new_primary(None, rec(50, qcov=0.8))[0], "consensus")
        self.assertEqual(G.new_primary(None, rec(50, qcov=0.79))[0], "none")
        self.assertEqual(G.new_primary(None, None)[0], "none")

    def test_genome_only(self):
        self.assertEqual(G.new_primary(rec(100), None)[0], "genome")


class Moves(unittest.TestCase):
    def test_counts_moves_correct_moves_and_de_change(self):
        cons_fam = {"cA": "F1", "cB": "F2"}
        rows = [  # (class, true family, genome rec or None, consensus rec or None, consensus name)
            ("D_abs", "F1", rec(100, 0.03), rec(120, 0.001, ref="cA"), "cA"),     # moves, correct
            ("D_abs", "F1", rec(100, 0.03), rec(150, 0.001, ref="cB"), "cB"),     # moves, wrong family
            ("D_abs", "F1", rec(100, 0.03), None, None),                           # stays
            ("S", "F9", rec(100, 0.001), rec(130, 0.0005, ref="cA"), "cA"),       # false move
            ("S", "F9", rec(100, 0.001), None, None),
        ]
        m = G.move_metrics(rows, cons_fam)
        self.assertEqual(m["D_abs"]["n"], 3)
        self.assertEqual((m["D_abs"]["moved"], m["D_abs"]["moved_to_own_family"]), (2, 1))
        self.assertAlmostEqual(m["D_abs"]["median_de_before_moved"], 0.03)
        self.assertEqual((m["S"]["n"], m["S"]["moved"]), (2, 1))
        self.assertAlmostEqual(m["S"]["move_fraction"], 0.5)

    def test_unmapped_reads_rescued_to_own_family_need_identity(self):
        cons_fam = {"cA": "F1", "cB": "F2"}
        rows = [
            ("D_unm", "F1", None, rec(100, 0.005, ref="cA"), "cA"),   # rescued, own family
            ("D_unm", "F1", None, rec(100, 0.05, ref="cA"), "cA"),    # identity below 0.98: not a move
            ("D_unm", "F1", None, rec(100, 0.005, ref="cB"), "cB"),   # rescued onto another family
            ("D_unm", "F1", None, None, None),
        ]
        m = G.move_metrics(rows, cons_fam, min_identity=0.98)
        self.assertEqual((m["D_unm"]["n"], m["D_unm"]["moved"], m["D_unm"]["moved_to_own_family"]), (4, 2, 1))


if __name__ == "__main__":
    unittest.main()
