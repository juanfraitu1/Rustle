#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_common.py"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import testutil  # noqa: E402
from common import Rec  # noqa: E402


def rec(primary=True, mapq=60, score=1000, qcov=1.0, ref="chr1", start=100, end=200, de=0.01):
    return Rec(primary, mapq, score, qcov, ref, start, end, de)


class Fate(unittest.TestCase):
    def test_no_record_is_unmapped(self):
        self.assertEqual(C.classify_fate([], None), "UNMAPPED")

    def test_only_secondary_is_unmapped(self):
        self.assertEqual(C.classify_fate([rec(primary=False)], None), "UNMAPPED")

    def test_low_coverage_primary_is_partial(self):
        self.assertEqual(C.classify_fate([rec(qcov=0.79)], None), "PARTIAL")

    def test_coverage_boundary_is_mapped(self):
        self.assertEqual(C.classify_fate([rec(qcov=0.80)], None), "ABSORBED_OTHER")

    def test_mapq0_is_tied(self):
        self.assertEqual(C.classify_fate([rec(mapq=0)], None), "TIED")

    def test_close_secondary_is_tied(self):
        self.assertEqual(C.classify_fate([rec(score=1000), rec(primary=False, score=981)], None), "TIED")

    def test_far_secondary_is_not_tied(self):
        self.assertEqual(C.classify_fate([rec(score=1000), rec(primary=False, score=979)], None), "ABSORBED_OTHER")

    def test_absorbed_on_nearest_paralog(self):
        self.assertEqual(C.classify_fate([rec(ref="A", start=100, end=200)], ("A", 150, 400)), "ABSORBED_NEAREST")

    def test_other_chromosome_is_other(self):
        self.assertEqual(C.classify_fate([rec(ref="B", start=100, end=200)], ("A", 150, 400)), "ABSORBED_OTHER")


class Place(unittest.TestCase):
    def test_untied_best(self):
        p = C.place([rec(score=1000, start=1), rec(primary=False, score=900, start=2)])
        self.assertEqual(p.start, 1)

    def test_tied_is_none(self):
        self.assertIsNone(C.place([rec(score=1000), rec(primary=False, score=985)]))

    def test_empty_is_none(self):
        self.assertIsNone(C.place([]))


class Names(unittest.TestCase):
    def test_index_name_maps_to_accession(self):
        self.assertEqual(C.accession("chr3_mat_hsa4", {("mat", "3"): "CM1"}), "CM1")

    def test_unknown_name_passes_through(self):
        self.assertEqual(C.accession("scaffold_12", {}), "scaffold_12")
        self.assertEqual(C.accession("chr9_mat_hsa9", {}), "chr9_mat_hsa9")


class Bam(unittest.TestCase):
    def test_read_records(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.bam")
            testutil.write_bam(p, {"chr1_mat_hsa1": 5000}, [
                dict(name="r1", ref="chr1_mat_hsa1", start=100, AS=190, de=0.01),
                dict(name="r1", flag=256, ref="chr1_mat_hsa1", start=900, mapq=0, AS=185, de=0.02),
                dict(name="r1", flag=2048, ref="chr1_mat_hsa1", start=900, cigar="50S50M", mapq=0, AS=90),
                dict(name="r2", flag=4),
                dict(name="r3", ref="chr1_mat_hsa1", start=300, cigar="20S80M", AS=70),
            ])
            recs = C.read_records(p, {("mat", "1"): "CM1"})
            self.assertEqual(len(recs["r1"]), 2)
            self.assertTrue(recs["r1"][0].primary)
            self.assertEqual(recs["r1"][0].score, 190)
            self.assertEqual(recs["r1"][0].ref, "CM1")
            self.assertAlmostEqual(recs["r1"][0].qcov, 1.0)
            self.assertFalse(recs["r1"][1].primary)
            self.assertEqual(recs["r2"], [])
            self.assertAlmostEqual(recs["r3"][0].qcov, 0.8)


class Paf(unittest.TestCase):
    def test_hits_and_best(self):
        line = lambda q, t, ts, te, m, aln, qs, qe, ql: "\t".join(
            [q, str(ql), str(qs), str(qe), "+", t, "9999", str(ts), str(te), str(m), str(aln), "60"])
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.paf")
            open(p, "w").write("\n".join([
                line("q1", "chr1_mat_hsa1", 10, 1010, 990, 1000, 0, 1000, 1000),      # ident .99, cov 1.0
                line("q1", "chr2_mat_hsa2", 50, 550, 480, 500, 0, 500, 1000),         # cov .5: not a hit
                line("q2", "chr1_mat_hsa1", 70, 1070, 800, 1000, 0, 1000, 1000),      # ident .80: not a hit
            ]) + "\n")
            al = {("mat", "1"): "CM1", ("mat", "2"): "CM2"}
            h = C.paf_hits(p, al)
            self.assertEqual(h["q1"], [("CM1", 10, 1010, 0.99, 1.0)])
            self.assertNotIn("q2", h)
            b = C.best_hits(p, al)
            self.assertEqual(b["q1"][1:4], ("CM1", 10, 1010))
            self.assertAlmostEqual(b["q1"][0], 0.99)


if __name__ == "__main__":
    unittest.main()
