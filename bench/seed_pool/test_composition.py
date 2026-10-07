#!/usr/bin/env python3
"""Tests of composition.py (docs/PREREG_seed_pool_real_reads_2026-10-07.md section 3). Run: python3 -m unittest bench/seed_pool/test_composition.py"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import composition as C  # noqa: E402

NODES = os.path.join(HERE, "..", "default_rescore", "nodes.py")


def load(path):
    with open(path) as fh:
        return json.load(fh)


def locus(c, s1, e, st, exons):
    """a locus as read_loci returns it: gene row span (1-based s1..e) and exons (0-based half-open)."""
    return dict(c=c, s1=s1, e=e, st=st, ex=exons)


def copy(chrom, strand, exons, span):
    return dict(chrom=chrom, strand=strand, exons=exons, span=span)


COPIES = {
    "c1": copy("chr1", "+", [[1000, 1100], [2000, 2100]], (1000, 2100)),
    "c2": copy("chr1", "+", [[5000, 5100], [6000, 6100]], (5000, 6100)),
    "c3": copy("chr1", "-", [[9000, 9100]], (9000, 9100)),
}


class JoinTests(unittest.TestCase):
    def test_a_locus_takes_the_cluster_ids_of_its_span(self):
        loci = {"A": locus("chr1", 1001, 2100, "+", [(1000, 1100)])}
        ids, shared = C.cluster_ids(loci, [("chr1", 1001, 2100, "MCL0")])
        self.assertEqual(ids, {"A": {"MCL0"}})
        self.assertEqual(shared, 0)

    def test_loci_on_one_span_share_the_ids_of_that_span_and_are_counted(self):
        loci = {"A": locus("chr1", 1001, 2100, "+", [(1000, 1100)]), "B": locus("chr1", 1001, 2100, "+", [(1000, 1100), (2000, 2100)])}
        ids, shared = C.cluster_ids(loci, [("chr1", 1001, 2100, "MCL0"), ("chr1", 1001, 2100, "MCL7")])
        self.assertEqual(ids["A"], {"MCL0", "MCL7"})
        self.assertEqual(ids["B"], {"MCL0", "MCL7"})
        self.assertEqual(shared, 1)

    def test_a_locus_without_a_cluster_row_is_unclustered(self):
        loci = {"A": locus("chr1", 1001, 2100, "+", [(1000, 1100)])}
        ids, _ = C.cluster_ids(loci, [])
        self.assertEqual(ids, {"A": set()})

    def test_a_cluster_row_on_no_locus_stops_the_run(self):
        with self.assertRaises(SystemExit):
            C.cluster_ids({"A": locus("chr1", 1001, 2100, "+", [(1000, 1100)])}, [("chr1", 5, 6, "MCL0")])


class AnalyseTests(unittest.TestCase):
    def setUp(self):
        self.loci = {
            "L1": locus("chr1", 1001, 2100, "+", [(1000, 1100), (2000, 2100)]),   # on c1
            "L2": locus("chr1", 5001, 6100, "+", [(5000, 5100)]),                   # on c2
            "L3": locus("chr1", 9001, 9100, "-", [(9000, 9100)]),                   # on c3, own cluster
            "L4": locus("chr1", 1301, 1400, "+", [(1300, 1400)]),                   # inside c1's span, off its exons (in_span)
            "L5": locus("chr1", 5001, 5100, "-", [(5000, 5100)]),                   # antisense to c2
            "L6": locus("chr1", 20001, 20100, "+", [(20000, 20100)]),               # elsewhere
            "L7": locus("chr1", 5601, 5700, "+", [(5600, 5700)]),                   # in c2's span, unclustered
        }
        self.rows = [("chr1", 1001, 2100, "MCL0"), ("chr1", 5001, 6100, "MCL0"), ("chr1", 9001, 9100, "MCL1"),
                     ("chr1", 1301, 1400, "MCL0"), ("chr1", 5001, 5100, "MCL0"), ("chr1", 20001, 20100, "MCL0")]

    def test_the_family_clusters_are_those_that_hold_a_locus_on_a_copy(self):
        r = C.analyse(self.loci, self.rows, COPIES)
        self.assertEqual(r["family_clusters"], ["MCL0", "MCL1"])
        self.assertEqual(r["own"], ["c1", "c2", "c3"])

    def test_node_classes_and_precision(self):
        r = C.analyse(self.loci, self.rows, COPIES)
        self.assertEqual(r["classes"], {"on_copy": 3, "in_span": 1, "antisense": 1, "elsewhere": 1})
        self.assertEqual(r["nodes"], 6)           # L7 is unclustered: not a node of a family cluster
        self.assertAlmostEqual(r["np"], 3 / 6)

    def test_a_copy_whose_only_locus_is_unclustered_has_no_own_node(self):
        rows = [x for x in self.rows if x[3] != "MCL1"]            # L3 loses its row
        r = C.analyse(self.loci, rows, COPIES)
        self.assertEqual(r["own"], ["c1", "c2"])

    def test_a_locus_on_the_wrong_strand_is_not_on_the_copy(self):
        loci = {"A": locus("chr1", 1001, 2100, "-", [(1000, 1100)])}
        r = C.analyse(loci, [("chr1", 1001, 2100, "MCL0")], COPIES)
        self.assertEqual(r["own"], [])
        self.assertEqual(r["family_clusters"], [])
        self.assertEqual(r["nodes"], 0)


class CohesionTests(unittest.TestCase):
    """M6: the copies that sit in ONE family cluster (the cluster with the most distinct copies, ties to the smaller id)."""

    def test_kstar_is_the_family_cluster_holding_most_copies(self):
        loci = {"L1": locus("chr1", 1001, 2100, "+", [(1000, 1100)]), "L2": locus("chr1", 5001, 6100, "+", [(5000, 5100)]),
                "L3": locus("chr1", 9001, 9100, "-", [(9000, 9100)])}
        rows = [("chr1", 1001, 2100, "MCL0"), ("chr1", 5001, 6100, "MCL0"), ("chr1", 9001, 9100, "MCL1")]
        r = C.analyse(loci, rows, COPIES)
        self.assertEqual(r["kstar"], "MCL0")
        self.assertEqual(r["copies_in_kstar"], 2)
        self.assertEqual(r["copies_with_node"], 3)

    def test_a_tie_goes_to_the_smaller_cluster_id_and_no_family_cluster_means_zero(self):
        loci = {"L1": locus("chr1", 1001, 2100, "+", [(1000, 1100)]), "L3": locus("chr1", 9001, 9100, "-", [(9000, 9100)])}
        r = C.analyse(loci, [("chr1", 1001, 2100, "MCL7"), ("chr1", 9001, 9100, "MCL3")], COPIES)
        self.assertEqual((r["kstar"], r["copies_in_kstar"]), ("MCL3", 1))
        r = C.analyse(loci, [], COPIES)
        self.assertEqual((r["kstar"], r["copies_in_kstar"]), (None, 0))

    def test_a_copy_on_a_locus_of_two_clusters_counts_in_each(self):
        loci = {"A": locus("chr1", 1001, 2100, "+", [(1000, 1100)]), "B": locus("chr1", 1001, 2100, "+", [(1000, 1100)]),
                "L2": locus("chr1", 5001, 6100, "+", [(5000, 5100)])}
        r = C.analyse(loci, [("chr1", 1001, 2100, "MCL0"), ("chr1", 1001, 2100, "MCL1"), ("chr1", 5001, 6100, "MCL1")], COPIES)
        self.assertEqual((r["kstar"], r["copies_in_kstar"]), ("MCL1", 2))


class AgreementWithNodesPy(unittest.TestCase):
    """composition.py's own-node flags equal nodes.py's wherever nodes.py is defined (no shared span)."""

    def test_same_flags_on_a_small_case(self):
        with tempfile.TemporaryDirectory() as d:
            copies = [("c1", "chr1", 1000, 2100, "+"), ("c2", "chr1", 5000, 6100, "+"), ("c3", "chr1", 9000, 9100, "-")]
            with open(f"{d}/copies.tsv", "w") as fh:
                fh.write("cid\tfamily\tname\tchrom\tterr_lo0\tterr_hi\tstrand\tisoform_gene\n")
                for cid, ch, lo, hi, st in copies:
                    fh.write(f"{cid}\tFAM\t{cid}\t{ch}\t{lo}\t{hi}\t{st}\t{cid}\n")
            with open(f"{d}/truth.gtf", "w") as fh:
                for cid, exs, st in (("c1", [(1000, 1100), (2000, 2100)], "+"), ("c2", [(5000, 5100), (6000, 6100)], "+"), ("c3", [(9000, 9100)], "-")):
                    for s0, e in exs:
                        fh.write(f'chr1\tx\texon\t{s0 + 1}\t{e}\t.\t{st}\t.\tgene_id "{cid}"; transcript_id "{cid}.1";\n')
            loci = {"L1": ("chr1", 1001, 2100, "+", [(1000, 1100), (2000, 2100)]), "L2": ("chr1", 5001, 5100, "+", [(5000, 5100)]),
                    "L3": ("chr1", 9001, 9100, "-", [(9000, 9100)]), "L6": ("chr1", 20001, 20100, "+", [(20000, 20100)])}
            with open(f"{d}/loci.gff3", "w") as fh:
                fh.write("##gff-version 3\n")
                for n, (ch, s1, e, st, exs) in loci.items():
                    fh.write(f"{ch}\t.\tgene\t{s1}\t{e}\t.\t{st}\t.\tID=gene-{n};Name={n}\n")
                    for s0, ee in exs:
                        fh.write(f"{ch}\t.\texon\t{s0 + 1}\t{ee}\t.\t{st}\t.\tParent=gene-{n};gene={n}\n")
            with open(f"{d}/clusters.tsv", "w") as fh:
                fh.write("cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n")
                for cl, n in (("MCL0", "L1"), ("MCL0", "L2"), ("MCL1", "L3"), ("MCL0", "L6")):
                    ch, s1, e, _, _ = loci[n]
                    fh.write(f"{cl}\t2\t1\t0\tNA\t{ch}\t{s1}\t{e}\n")
            arm = f"A={d}/loci.gff3,{d}/clusters.tsv"
            base = ["--copies", f"{d}/copies.tsv", "--truth", f"{d}/truth.gtf", "--family", "FAM", "--arm", arm]
            subprocess.run([sys.executable, NODES, *base, "--out", f"{d}/nodes_old.json"], check=True, capture_output=True)
            subprocess.run([sys.executable, os.path.join(HERE, "composition.py"), *base, "--out", f"{d}/nodes_new.json", "--report", f"{d}/comp.json"],
                           check=True, capture_output=True)
            old = {r["cid"]: r["node"] for r in load(f"{d}/nodes_old.json")["rows"]}
            new = {r["cid"]: r["node"] for r in load(f"{d}/nodes_new.json")["rows"]}
            self.assertEqual(old, new)
            self.assertEqual(old, {"c1": {"A": True}, "c2": {"A": True}, "c3": {"A": True}})
            rep = load(f"{d}/comp.json")["arms"]["A"]
            self.assertEqual(rep["classes"]["elsewhere"], 1)
            self.assertEqual(rep["copies_with_node"], 3)


if __name__ == "__main__":
    unittest.main()
