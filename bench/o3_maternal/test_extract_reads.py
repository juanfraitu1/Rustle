#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_extract_reads.py"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import extract_reads as E  # noqa: E402
import testutil  # noqa: E402


class Net(unittest.TestCase):
    def bam(self, d):
        p = os.path.join(d, "t.bam")
        testutil.write_bam(p, {"chrA": 100000}, [
            dict(name="a", ref="chrA", start=1000),
            dict(name="b", ref="chrA", start=1500),
            dict(name="b", flag=256, ref="chrA", start=40000, mapq=0),       # secondary on locus 2 only
            dict(name="c", ref="chrA", start=90000),                          # outside every locus
            dict(name="s", flag=2048, ref="chrA", start=1200, cigar="50S50M"),  # supplementary: ignored
        ])
        return p

    def test_reads_on_a_locus_primary_or_secondary(self):
        with tempfile.TemporaryDirectory() as d:
            n = E.net_names(self.bam(d), [("c1", "g1", "chrA", 500, 3000), ("c2", "g2", "chrA", 39000, 42000)])
            self.assertEqual(n, {"a": ["c1"], "b": ["c1", "c2"]})

    def test_cap_is_per_locus_and_seeded(self):
        with tempfile.TemporaryDirectory() as d:
            n1 = E.net_names(self.bam(d), [("c1", "g1", "chrA", 500, 3000)], cap=1)
            n2 = E.net_names(self.bam(d), [("c1", "g1", "chrA", 500, 3000)], cap=1)
            self.assertEqual(len(n1), 1)
            self.assertEqual(n1, n2)


class Merge(unittest.TestCase):
    def test_overlap_read_keeps_family_label(self):
        r34 = {"x": "GWFAM9"}
        rows, overlap = E.merge_labels(r34, {"x": ["p12"], "y": ["p12", "p14"]})
        self.assertEqual(overlap, ["x"])
        self.assertEqual(rows, [("y", "LRPAP1", "p12,p14")])


if __name__ == "__main__":
    unittest.main()
