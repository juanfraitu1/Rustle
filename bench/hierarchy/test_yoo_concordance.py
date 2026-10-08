#!/usr/bin/env python3
"""Tests of the Y block (concordance with Yoo et al. 2025) of docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md, Amendment 3.

Run: python3 -B -m unittest bench/hierarchy/test_yoo_concordance.py

Written BEFORE the implementation. Synthetic tables only; the real inputs are never read here.
"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import yoo_concordance as yc  # noqa: E402


class MatchRule(unittest.TestCase):
    def test_overlap_is_closed_interval_bp(self):
        self.assertEqual(yc.overlap((1, 10), (10, 20)), 1)
        self.assertEqual(yc.overlap((1, 10), (11, 20)), 0)
        self.assertEqual(yc.overlap((5, 15), (1, 100)), 11)

    def test_exactly_half_of_the_shorter_matches(self):
        a = ("c", 1, 100)
        self.assertTrue(yc.matches(a, ("c", 51, 150)))   # overlap 50 = 0.5 * 100
        self.assertFalse(yc.matches(a, ("c", 52, 151)))  # overlap 49

    def test_shorter_locus_sets_the_denominator(self):
        self.assertTrue(yc.matches(("c", 1, 1000), ("c", 10, 19)))
        self.assertFalse(yc.matches(("c", 1, 1000), ("c", 995, 1100)))  # 6 of 106

    def test_other_chromosome_never_matches(self):
        self.assertFalse(yc.matches(("c1", 1, 100), ("c2", 1, 100)))

    def test_adjacent_loci_do_not_match(self):
        self.assertFalse(yc.matches(("c", 1, 100), ("c", 101, 200)))


class OneToOne(unittest.TestCase):
    def test_sizes(self):
        self.assertEqual(yc.max_one_to_one([(0, "a"), (0, "b"), (1, "a")]), 2)
        self.assertEqual(yc.max_one_to_one([(0, "a"), (1, "a")]), 1)
        self.assertEqual(yc.max_one_to_one([]), 0)

    def test_augmenting_path_is_needed(self):
        # greedy 0->a blocks 1; a maximum matching reroutes 0->b
        self.assertEqual(yc.max_one_to_one([(0, "a"), (0, "b"), (1, "a"), (2, "b")]), 2)
        self.assertEqual(yc.max_one_to_one([(0, "a"), (1, "a"), (1, "b")]), 2)


class Readers(unittest.TestCase):
    def test_unit_loci_get_refseq_names(self):
        with tempfile.TemporaryDirectory() as d:
            chrmap = os.path.join(d, "map.tsv")
            with open(chrmap, "w") as f:
                f.write("chr\tgenbank\trefseq\tlength\n1\tCM1\tNC_1.2\t100\n16\tCM16\tNC_16.2\t90\n")
            unit = os.path.join(d, "unit.tsv")
            with open(unit, "w") as f:
                f.write("chr\trefseq\tstart\tend\tgene\tflag\nchr1\tNC_1.2\t10\t20\tMAPKBP1\t\nchr16\tNC_16.2\t5\t50\tSPTBN5\t***\n")
            self.assertEqual(yc.read_chr_map(chrmap), {"chr1": "NC_1.2", "chr16": "NC_16.2"})
            loci = yc.read_unit_loci(unit)
            self.assertEqual(loci, [("MAPKBP1", "NC_1.2", 10, 20), ("SPTBN5", "NC_16.2", 5, 50)])

    def test_clusters_skip_the_header(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "clusters.tsv")
            with open(p, "w") as f:
                f.write("cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\nMCL0\t2\t1\t1\tNA\tNC_1.2\t11\t30\nMCL0\t2\t1\t1\tNA\tNC_1.2\t51\t90\n")
            self.assertEqual(yc.read_clusters(p), [("MCL0", "NC_1.2", 11, 30), ("MCL0", "NC_1.2", 51, 90)])


class CopiesReader(unittest.TestCase):
    def test_copy_rows_become_one_based_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "copies.tsv")
            hdr = "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\tsource\tgene_id\tcore_hull\tsd_depth\tcore_bp\trep_frac\tmember_status\tlocus_start\tlocus_end\n"
            with open(p, "w") as f:
                f.write(hdr)
                f.write("MCL0\t0\tDN_x\tNC_1.2\t48946409\t48952180\t2\t+\t2\t1-2\t1.0\tlocus_rep\tg\tNA\tNA\tNA\tNA\tungated\t48946000\t48952500\n")
            self.assertEqual(yc.read_copies(p), [("MCL0", "NC_1.2", 48946410, 48952180)])
            self.assertEqual(yc.read_copies(p, use_locus=True), [("MCL0", "NC_1.2", 48946001, 48952500)])


class PerGene(unittest.TestCase):
    def setUp(self):
        # two genes: A has two Yoo loci, B one; family F1 holds three members, F2 one
        self.yoo = [("A", "c", 100, 199), ("A", "c", 300, 399), ("B", "c", 1000, 1099)]
        self.clusters = [
            ("F1", "c", 100, 199),    # matches A#0
            ("F1", "c", 1000, 1099),  # matches B#0
            ("F1", "c", 5000, 5099),  # matches nothing
            ("F2", "c", 320, 380),    # matches A#1
        ]

    def test_counts(self):
        r = yc.per_gene(self.yoo, self.clusters)
        a, b = r["A"], r["B"]
        self.assertEqual((a["n_yoo"], a["yoo_matched"], a["one_to_one"]), (2, 2, 2))
        self.assertEqual(a["families"], ["F1", "F2"])
        self.assertEqual(a["members_matching"], 2)
        self.assertEqual((b["n_yoo"], b["yoo_matched"], b["one_to_one"]), (1, 1, 1))
        self.assertEqual(b["families"], ["F1"])

    def test_unmatched_members_of_touched_families(self):
        r = yc.per_gene(self.yoo, self.clusters)
        # F1 has 3 members and 2 match some Yoo locus (of any gene): 1 matches none
        self.assertEqual(r["A"]["touched_members_without_any_yoo_locus"], 1)
        self.assertEqual(r["B"]["touched_members_without_any_yoo_locus"], 1)

    def test_one_member_covering_two_loci_counts_two_loci_but_matching_one(self):
        yoo = [("A", "c", 100, 199), ("A", "c", 200, 299)]
        clusters = [("F", "c", 100, 299)]  # one long member covers both
        r = yc.per_gene(yoo, clusters)["A"]
        self.assertEqual((r["yoo_matched"], r["members_matching"], r["one_to_one"]), (2, 1, 1))

    def test_same_family_table(self):
        r = yc.per_gene(self.yoo, self.clusters)
        shared = yc.shared_families(r)
        self.assertEqual(shared, {("A", "B"): ["F1"]})

    def test_gene_with_no_match(self):
        r = yc.per_gene([("Z", "c", 1, 10)], [("F", "c", 100, 200)])["Z"]
        self.assertEqual((r["yoo_matched"], r["one_to_one"], r["families"], r["members_matching"]), (0, 0, [], 0))


if __name__ == "__main__":
    unittest.main()
