#!/usr/bin/env python3
"""Tests of locus_units.py (Rule 1 and Rule 2 of docs/PREREG_locus_units_2026-10-06.md). Run: python3 -m unittest bench/entangled/test_locus_units.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import locus_units as L


def J(*names):
    return tuple((n, n + 1000) for n in names)   # a chain of junctions named by one integer each


class SeparatorTests(unittest.TestCase):
    def test_isoforms_sharing_two_junctions_have_no_separator(self):
        sep, comp = L.separators_components({"t1": J(1, 2, 3), "t2": J(1, 2, 4)})
        self.assertEqual(sep, set())
        self.assertEqual(comp["t1"], comp["t2"])

    def test_readthrough_containing_a_transcript_of_each_group_is_a_separator(self):
        chains = {"A1": J(1, 2), "A2": J(1, 3), "B1": J(11, 12), "B2": J(11, 13), "R": J(1, 2, 99, 11, 12)}
        sep, comp = L.separators_components(chains)
        self.assertEqual(sep, {"R"})
        self.assertEqual(comp["A1"], comp["A2"])
        self.assertEqual(comp["B1"], comp["B2"])
        self.assertNotEqual(comp["A1"], comp["B1"])
        self.assertNotIn("R", comp)

    def test_two_genes_sharing_one_junction_stay_together(self):
        chains = {"A1": J(1, 50), "A2": J(2, 50), "B1": J(11, 50), "B2": J(12, 50)}
        sep, comp = L.separators_components(chains)
        self.assertEqual(sep, set())            # the shared junction is not a transcript: not cut by this rule
        self.assertEqual(len(set(comp.values())), 1)

    def test_transcripts_sharing_nothing_are_their_own_components_and_never_separators(self):
        sep, comp = L.separators_components({"t1": J(1, 2), "t2": J(3, 4)})
        self.assertEqual(sep, set())
        self.assertNotEqual(comp["t1"], comp["t2"])

    def test_the_concatenation_of_two_single_transcript_groups_is_a_separator(self):
        sep, comp = L.separators_components({"T1": J(1, 2), "T": J(1, 2, 5, 8, 9), "T2": J(8, 9)})
        self.assertEqual(sep, {"T"})
        self.assertNotEqual(comp["T1"], comp["T2"])

    def test_a_private_junction_between_junctions_of_two_groups_makes_a_separator_even_without_containing_either_chain(self):
        sep, comp = L.separators_components({"T1": J(1, 2), "T": J(1, 5, 9), "T2": J(9, 8)})
        self.assertEqual(sep, {"T"})
        self.assertNotEqual(comp["T1"], comp["T2"])

    def test_a_link_without_a_private_junction_between_the_groups_is_not_a_separator(self):
        sep, comp = L.separators_components({"T1": J(1, 2), "T": J(1, 9), "T2": J(9, 8)})
        self.assertEqual(sep, set())
        self.assertEqual(len({comp["T1"], comp["T"], comp["T2"]}), 1)

    def test_an_isoform_that_is_the_only_link_to_a_one_transcript_attachment_is_not_a_separator(self):
        chains = {"A1": J(1, 2), "A2": J(1, 2, 3), "leaf": J(3, 77)}      # no junction of A2 lies between the two groups
        sep, comp = L.separators_components(chains)
        self.assertEqual(sep, set())
        self.assertEqual(len({comp["A1"], comp["A2"], comp["leaf"]}), 1)

    def test_the_private_junction_must_lie_between_the_two_groups_in_chain_order(self):
        # the private junction 99 comes after both group junctions: a tail, not a bridge
        sep, comp = L.separators_components({"T1": J(1, 2), "T": J(1, 9, 99), "T2": J(9, 8)})
        self.assertEqual(sep, set())


GTF = "\n".join([
    # locus g1: A group (a1, a2), B group (b1), readthrough R bridging them; plus a single-exon transcript s overlapping A
    'c1\tx\ttranscript\t1\t900\t.\t+\t.\tgene_id "g1"; transcript_id "A1"; reads "10";',
    'c1\tx\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "A1"; exon_number "1";',
    'c1\tx\texon\t200\t300\t.\t+\t.\tgene_id "g1"; transcript_id "A1"; exon_number "2";',
    'c1\tx\texon\t400\t500\t.\t+\t.\tgene_id "g1"; transcript_id "A1"; exon_number "3";',
    'c1\tx\ttranscript\t1\t900\t.\t+\t.\tgene_id "g1"; transcript_id "A2"; reads "6";',
    'c1\tx\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "A2"; exon_number "1";',
    'c1\tx\texon\t200\t300\t.\t+\t.\tgene_id "g1"; transcript_id "A2"; exon_number "2";',
    'c1\tx\texon\t600\t700\t.\t+\t.\tgene_id "g1"; transcript_id "A2"; exon_number "3";',
    'c1\tx\ttranscript\t1000\t1900\t.\t+\t.\tgene_id "g1"; transcript_id "B1"; reads "8";',
    'c1\tx\texon\t1000\t1100\t.\t+\t.\tgene_id "g1"; transcript_id "B1"; exon_number "1";',
    'c1\tx\texon\t1200\t1300\t.\t+\t.\tgene_id "g1"; transcript_id "B1"; exon_number "2";',
    'c1\tx\texon\t1400\t1500\t.\t+\t.\tgene_id "g1"; transcript_id "B1"; exon_number "3";',
    'c1\tx\ttranscript\t1\t1900\t.\t+\t.\tgene_id "g1"; transcript_id "R"; reads "9";',
    'c1\tx\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "R"; exon_number "1";',
    'c1\tx\texon\t200\t300\t.\t+\t.\tgene_id "g1"; transcript_id "R"; exon_number "2";',
    'c1\tx\texon\t400\t500\t.\t+\t.\tgene_id "g1"; transcript_id "R"; exon_number "3";',
    'c1\tx\texon\t1200\t1300\t.\t+\t.\tgene_id "g1"; transcript_id "R"; exon_number "4";',
    'c1\tx\texon\t1400\t1500\t.\t+\t.\tgene_id "g1"; transcript_id "R"; exon_number "5";',
    'c1\tx\ttranscript\t10\t60\t.\t+\t.\tgene_id "g1"; transcript_id "S"; reads "3";',
    'c1\tx\texon\t10\t60\t.\t+\t.\tgene_id "g1"; transcript_id "S"; exon_number "1";',
    # locus g2: a plain two-isoform gene, untouched by Rule 2
    'c1\tx\ttranscript\t5000\t5900\t.\t-\t.\tgene_id "g2"; transcript_id "G1"; reads "5";',
    'c1\tx\texon\t5000\t5100\t.\t-\t.\tgene_id "g2"; transcript_id "G1"; exon_number "1";',
    'c1\tx\texon\t5200\t5300\t.\t-\t.\tgene_id "g2"; transcript_id "G1"; exon_number "2";',
    'c1\tx\ttranscript\t5000\t5900\t.\t-\t.\tgene_id "g2"; transcript_id "G2"; reads "2";',
    'c1\tx\texon\t5000\t5100\t.\t-\t.\tgene_id "g2"; transcript_id "G2"; exon_number "1";',
    'c1\tx\texon\t5200\t5300\t.\t-\t.\tgene_id "g2"; transcript_id "G2"; exon_number "2";',
    'c1\tx\texon\t5400\t5500\t.\t-\t.\tgene_id "g2"; transcript_id "G2"; exon_number "3";',
]) + "\n"


class GtfTests(unittest.TestCase):
    def test_rule1_adds_a_large_constant_to_reads_of_primary_supported_chains_only(self):
        recs = L.parse_gtf(GTF)
        out = L.apply_rule1(recs, {("c1", "-", ((5100, 5199),))})     # G1's single junction (0-based end 5100, next start 5199)
        reads = {r.tid: r.attrs.get("reads") for r in out if r.feature == "transcript"}
        self.assertEqual(reads["G1"], str(5 + L.PRIMARY_BONUS))
        self.assertEqual(reads["G2"], "2")
        self.assertEqual(reads["A1"], "10")

    def test_rule2_drops_the_separator_and_splits_the_locus(self):
        recs = L.parse_gtf(GTF)
        out, side = L.apply_rule2(recs, pool=None)
        tids = {r.tid for r in out if r.feature == "transcript"}
        self.assertNotIn("R", tids)
        gid = {r.tid: r.gid for r in out if r.feature == "transcript"}
        self.assertEqual(gid["A1"], gid["A2"])
        self.assertNotEqual(gid["A1"], gid["B1"])
        self.assertEqual(gid["S"], gid["A1"])          # the single-exon transcript attaches to the group it overlaps
        self.assertEqual(gid["G1"], "g2")              # a locus without separator keeps its gene_id
        self.assertEqual(gid["G2"], "g2")
        self.assertEqual([s["tid"] for s in side], ["R"])

    def test_rule2_best_representative_group_keeps_the_gene_id(self):
        recs = L.parse_gtf(GTF)
        out, _ = L.apply_rule2(recs, pool=None)
        gid = {r.tid: r.gid for r in out if r.feature == "transcript"}
        self.assertEqual(gid["A1"], "g1")              # A1 carries the most reads of the locus
        self.assertTrue(gid["B1"].startswith("g1.s"))


if __name__ == "__main__":
    unittest.main()
