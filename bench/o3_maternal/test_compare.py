#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_compare.py"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare as K  # noqa: E402

LOCI = [dict(locus="A", chrom="CM1", start=100, end=200), dict(locus="B", chrom="CM1", start=500, end=600)]


class Outputs(unittest.TestCase):
    def test_outputs_on_the_locus_at_the_floor(self):
        best = {"o1": (1.0, "CM1", 120, 180, 1.0), "o2": (0.998, "CM1", 120, 180, 0.998), "o3": (1.0, "CM2", 120, 180, 1.0),
                "o4": (0.9995, "CM1", 550, 590, 0.9995)}
        self.assertEqual(K.outputs_on(best, LOCI), {"A": ["o1"], "B": ["o4"]})

    def test_new_outputs_are_the_unlinked_ones(self):
        rows = [dict(output="F|x one", linked="0"), dict(output="F|y", linked="1")]
        self.assertEqual(K.new_outputs(rows), {"F|x"})


class Inhouse(unittest.TestCase):
    def test_recovering(self):
        rows = [dict(candidate="k1", n_transcripts=3, recovers=["A"]), dict(candidate="k2", n_transcripts=2, recovers=["B"])]
        self.assertEqual(K.recovering(rows, "A"), [dict(candidate="k1", n_transcripts=3)])

    def test_near_uses_accessions_and_overlap(self):
        al = {("pat", "5"): "CM1"}
        rows = [dict(candidate="c1", n_clusters="4", d="0.0003", flagged="0", nearest_locus="chr5_pat_hsa5:150-400"),
                dict(candidate="c2", n_clusters="2", d="0.2", flagged="1", nearest_locus="chr5_pat_hsa5:1000-2000")]
        out = K.inhouse_near(rows, LOCI, al)
        self.assertEqual(out["A"], [dict(candidate="c1", n_clusters=4, d=0.0003, flagged=0)])
        self.assertEqual(out["B"], [])


if __name__ == "__main__":
    unittest.main()
