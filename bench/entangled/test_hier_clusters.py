#!/usr/bin/env python3
"""Tests of hier_clusters.py (Amendment 2 of docs/PREREG_locus_units_2026-10-06.md). Run: python3 -m unittest test_hier_clusters"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hier_clusters as H

LOCI = "\n".join([
    "##gff-version 3",
    "c1\t.\tgene\t101\t500\t.\t+\t.\tID=gene-a;Name=a",
    "c1\t.\texon\t101\t200\t.\t+\t.\tParent=gene-a;gene=a",
    "c1\t.\tgene\t1001\t1500\t.\t+\t.\tID=gene-b;Name=b",
    "c1\t.\texon\t1001\t1100\t.\t+\t.\tParent=gene-b;gene=b",
    "c2\t.\tgene\t101\t500\t.\t-\t.\tID=gene-c;Name=c",
    "c2\t.\texon\t101\t200\t.\t-\t.\tParent=gene-c;gene=c",
    "c2\t.\tgene\t2001\t2500\t.\t-\t.\tID=gene-d;Name=d",
    "c2\t.\texon\t2001\t2100\t.\t-\t.\tParent=gene-d;gene=d",
    "c3\t.\tgene\t101\t500\t.\t+\t.\tID=gene-e;Name=e",
    "c3\t.\texon\t101\t200\t.\t+\t.\tParent=gene-e;gene=e",
]) + "\n"
GRAPH = "c1:101-500\tc1:1001-1500\t0.9\nc1:1001-1500\tc2:101-500\t0.4\nc2:2001-2500\tc9:1-2\t0.7\n"


class HierTests(unittest.TestCase):
    def test_components_of_the_graph_become_clusters_over_locus_keys_only(self):
        with tempfile.TemporaryDirectory() as d:
            open(d + "/x.fam.loci.gff3", "w").write(LOCI)
            open(d + "/g.tsv", "w").write(GRAPH)
            rows = H.components(d + "/x.fam.loci.gff3", d + "/g.tsv")
        by = {}
        for cid, chrom, s, e in rows:
            by.setdefault(cid, []).append((chrom, s, e))
        sizes = sorted(len(v) for v in by.values())
        self.assertEqual(sizes, [3])                       # a (c1), b (c1), c (c2) are connected; d joins only an unknown node c9:1-2 (not a locus), e is a singleton
        members = next(iter(by.values()))
        self.assertEqual(sorted(members), [("c1", 101, 500), ("c1", 1001, 1500), ("c2", 101, 500)])

    def test_the_cluster_file_has_the_registered_header_and_one_row_per_member(self):
        with tempfile.TemporaryDirectory() as d:
            open(d + "/x.fam.loci.gff3", "w").write(LOCI)
            open(d + "/g.tsv", "w").write(GRAPH)
            rows = H.components(d + "/x.fam.loci.gff3", d + "/g.tsv")
            H.write_clusters(rows, d + "/out.fam.clusters.tsv")
            lines = open(d + "/out.fam.clusters.tsv").read().splitlines()
        self.assertEqual(lines[0].split("\t"), ["cluster_id", "size", "density", "frac_in", "corroborated", "chrom", "start", "end"])
        self.assertEqual(len(lines) - 1, 3)
        self.assertTrue(all(l.split("\t")[1] == "3" for l in lines[1:]))


if __name__ == "__main__":
    unittest.main()
