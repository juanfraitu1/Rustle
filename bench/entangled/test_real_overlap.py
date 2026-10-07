#!/usr/bin/env python3
"""Tests of real_truth.py and real_score.py (docs/PREREG_gorilla_overlap_2026-10-07.md). Run: python3 -m unittest bench/entangled/test_real_overlap.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import real_truth as T
import real_score as S

GFF = "\n".join([
    "##gff-version 3",
    "c1\tx\tgene\t101\t900\t.\t+\t.\tID=gene-A;Name=A;gene_biotype=protein_coding",
    "c1\tx\tmRNA\t101\t900\t.\t+\t.\tID=rna-A1;Parent=gene-A;gbkey=mRNA",
    "c1\tx\texon\t101\t200\t.\t+\t.\tID=exon-A1-1;Parent=rna-A1",
    "c1\tx\texon\t301\t400\t.\t+\t.\tID=exon-A1-2;Parent=rna-A1",
    "c1\tx\texon\t801\t900\t.\t+\t.\tID=exon-A1-3;Parent=rna-A1",
    "c1\tx\tmRNA\t101\t500\t.\t+\t.\tID=rna-A2;Parent=gene-A",
    "c1\tx\texon\t101\t200\t.\t+\t.\tID=exon-A2-1;Parent=rna-A2",
    "c1\tx\texon\t301\t400\t.\t+\t.\tID=exon-A2-2;Parent=rna-A2",
    "c1\tx\tpseudogene\t1001\t1500\t.\t-\t.\tID=gene-P;Name=P;gene_biotype=pseudogene",
    "c1\tx\texon\t1001\t1100\t.\t-\t.\tID=exon-P-1;Parent=gene-P",
    "c1\tx\texon\t1401\t1500\t.\t-\t.\tID=exon-P-2;Parent=gene-P",
    "c1\tx\tgene\t301\t700\t.\t+\t.\tID=gene-B;Name=B;gene_biotype=lncRNA",
    "c1\tx\tlnc_RNA\t301\t700\t.\t+\t.\tID=rna-B1;Parent=gene-B",
    "c1\tx\texon\t301\t400\t.\t+\t.\tID=exon-B1-1;Parent=rna-B1",
    "c1\tx\texon\t601\t700\t.\t+\t.\tID=exon-B1-2;Parent=rna-B1",
]) + "\n"


class ParseTests(unittest.TestCase):
    def test_genes_transcripts_and_exons_are_zero_based_half_open(self):
        genes = T.parse_gff(GFF)
        self.assertEqual(set(genes), {"gene-A", "gene-P", "gene-B"})
        a = genes["gene-A"]
        self.assertEqual((a["chrom"], a["strand"], a["name"], a["biotype"]), ("c1", "+", "A", "protein_coding"))
        self.assertEqual(a["transcripts"]["rna-A1"], [(100, 200), (300, 400), (800, 900)])
        self.assertEqual(a["transcripts"]["rna-A2"], [(100, 200), (300, 400)])

    def test_a_pseudogene_with_exons_directly_under_the_gene_gets_one_transcript(self):
        p = T.parse_gff(GFF)["gene-P"]
        self.assertEqual(list(p["transcripts"].values()), [[(1000, 1100), (1400, 1500)]])
        self.assertEqual(p["strand"], "-")


def motif_ok(chrom, strand, i0, i1):
    return ("GT", "AG")


def motif_bad_second(chrom, strand, i0, i1):
    return ("GT", "AG") if i0 < 300 else ("AA", "CC")


class ChainTests(unittest.TestCase):
    def test_valid_chain_is_the_ordered_intron_list_of_canonical_long_enough_introns(self):
        self.assertEqual(T.valid_chain("c1", "+", [(100, 200), (300, 400), (800, 900)], motif_ok), ((200, 300), (400, 800)))

    def test_a_short_intron_makes_the_chain_invalid(self):
        self.assertIsNone(T.valid_chain("c1", "+", [(100, 200), (230, 400)], motif_ok))

    def test_a_non_canonical_intron_makes_the_chain_invalid(self):
        self.assertIsNone(T.valid_chain("c1", "+", [(100, 200), (300, 400), (800, 900)], motif_bad_second))

    def test_single_exon_transcripts_have_no_chain(self):
        self.assertIsNone(T.valid_chain("c1", "+", [(100, 200)], motif_ok))


class StrataTests(unittest.TestCase):
    def setUp(self):
        genes = T.parse_gff(GFF)
        for g in genes.values():
            g["chains"] = {}
            for tid, ex in g["transcripts"].items():
                ch = T.valid_chain(g["chrom"], g["strand"], ex, motif_ok)
                if ch:
                    g["chains"][tid] = ch
        self.genes = genes

    def test_same_strand_exon_overlap_of_100_bp_is_entangled_and_junction_sharing_is_detected(self):
        # A (exons 100-200, 300-400, 800-900) and B (300-400, 600-700) share 100 exonic bp; they share no junction
        expressed = {"gene-A": True, "gene-B": True, "gene-P": False}
        st = T.strata(self.genes, expressed)
        self.assertEqual(st["gene-A"]["stratum"], "E_both")
        self.assertEqual(st["gene-B"]["stratum"], "E_both")
        self.assertEqual(st["gene-A"]["junction_sharing"], 0)
        self.assertEqual(st["gene-P"]["stratum"], "N")

    def test_an_unexpressed_partner_gives_E_one(self):
        st = T.strata(self.genes, {"gene-A": True, "gene-B": False, "gene-P": False})
        self.assertEqual(st["gene-A"]["stratum"], "E_one")

    def test_opposite_strand_overlap_only_is_A(self):
        genes = self.genes
        genes["gene-B"]["strand"] = "-"
        st = T.strata(genes, {"gene-A": True, "gene-B": True, "gene-P": False})
        self.assertEqual(st["gene-A"]["stratum"], "A")


class SupportTests(unittest.TestCase):
    def test_pools_count_distinct_reads_and_apply_the_secondary_as_rule(self):
        chains = {("c1", "+", ((200, 300), (400, 800))): "gene-A"}
        spans = {"gene-A": ("c1", "+", 100, 900)}
        best_as = {"r1": 1000, "r2": 1000, "r3": 1000, "r4": 1000, "r5": 1000}
        aln = [  # (name, flag, AS, contig, junction list)
            ("r1", 0, 1000, "c1", [(200, 300), (400, 800)]),
            ("r2", 16, 1000, "c1", [(200, 300), (400, 800)]),
            ("r3", 256, 985, "c1", [(200, 300), (400, 800)]),      # secondary, AS >= .98 x best: P2 only
            ("r4", 256, 900, "c1", [(200, 300), (400, 800)]),      # secondary below .98: neither pool
            ("r5", 2048, 1000, "c1", [(200, 300), (400, 800)]),    # supplementary: ignored
            ("r1", 256, 995, "c1", [(200, 300), (400, 800)]),      # the same read again: counted once
            ("r6", 0, 1000, "c1", [(200, 300)]),                   # a different chain
        ]
        p1, p2 = T.chain_support(aln, chains, spans, best_as)
        key = ("c1", "+", ((200, 300), (400, 800)))
        self.assertEqual(len(p1[key]), 2)
        self.assertEqual(len(p2[key]), 3)

    def test_only_junctions_inside_the_gene_span_count(self):
        chains = {("c1", "+", ((200, 300),)): "gene-A"}
        spans = {"gene-A": ("c1", "+", 100, 400)}
        aln = [("r1", 0, 1000, "c1", [(200, 300), (5000, 5200)])]   # the second junction lies outside the span
        p1, _ = T.chain_support(aln, chains, spans, {"r1": 1000})
        self.assertEqual(len(p1[("c1", "+", ((200, 300),))]), 1)


class BamTests(unittest.TestCase):
    def test_junctions_follow_the_eqx_cigar_convention(self):
        import pysam
        rec = pysam.AlignedSegment()
        rec.reference_start = 1000
        rec.cigarstring = "50=200N60=10D30=100N20X5S"
        self.assertEqual(T.junctions_of(rec), [(1050, 1250), (1350, 1450)])

    def test_introns_shorter_than_50_are_not_junctions_but_still_advance_the_reference(self):
        import pysam
        rec = pysam.AlignedSegment()
        rec.reference_start = 0
        rec.cigarstring = "30=40N30=100N30="
        self.assertEqual(T.junctions_of(rec), [(100, 200)])      # 30 + 40 (short intron, advances) + 30 = 100


class ScoreTests(unittest.TestCase):
    def test_recall_match_share_and_resolution_on_a_toy_arm(self):
        truth = dict(
            genes={"gA": dict(gene="gA", name="A", chrom="c1", strand="+", union=[(100, 900)], expressed={((200, 300),), ((200, 300), (400, 800))},
                              stratum="E_both", junction_sharing=0),
                   "gB": dict(gene="gB", name="B", chrom="c1", strand="+", union=[(300, 700)], expressed={((400, 600),)}, stratum="E_both", junction_sharing=0)},
            all_chains={("c1", "+", ((200, 300),)): {"gA"}, ("c1", "+", ((200, 300), (400, 800))): {"gA"}, ("c1", "+", ((400, 600),)): {"gB"}})
        arm = {"t1": dict(chrom="c1", strand="+", gene="x1", exons=[(100, 200), (300, 400), (800, 900)]),     # chain ((200,300),(400,800)): gA exact
               "t2": dict(chrom="c1", strand="+", gene="x1", exons=[(100, 200), (300, 500)]),                 # chain ((200,300),): gA exact
               "t3": dict(chrom="c1", strand="+", gene="x2", exons=[(300, 400), (600, 700)]),                 # chain ((400,600),): gB exact
               "t4": dict(chrom="c1", strand="+", gene="x2", exons=[(300, 400), (650, 700)])}                 # unannotated chain overlapping gB
        res = S.score_arm(arm, truth)
        e = res["strata"]["E_both"]
        self.assertEqual((e["chains"], e["recovered"]), (3, 3))
        self.assertEqual(e["genes"], 2)
        self.assertEqual(e["complete"], 2)
        self.assertEqual((e["match_n"], e["match_d"]), (3, 4))
        self.assertEqual(e["resolved"], 2)      # x1 carries only gA chains, x2 only gB chains
        arm["t3"]["gene"] = "x1"                # now one gene_id carries exact chains of both genes
        res2 = S.score_arm(arm, truth)
        self.assertEqual(res2["strata"]["E_both"]["resolved"], 0)
        self.assertEqual(res2["merged_ids"], 1)
