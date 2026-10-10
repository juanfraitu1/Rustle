#!/usr/bin/env python3
"""Unplaced O3 candidates as FASTA + GTF (+ BED of the nearest reference locus), for loading beside the assembler's GTF."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o3_export as X

CAND = dict(id="O3cand_000007", seq="ACGTACGTAA", cls="DIVERGED", reads=12, R=0.97312, median_read_divergence=0.0011,
            nearest=("NC_073236.2", 95963504, 96028675, "+"), read_share=0.4812)


class Gtf(unittest.TestCase):
    def test_transcript_and_exon_on_their_own_sequence(self):
        lines = X.gtf_lines(CAND)
        self.assertEqual(len(lines), 2)
        t, e = (l.split("\t") for l in lines)
        self.assertEqual(t[:8], ["O3cand_000007", "rustle_o3", "transcript", "1", "10", "12", "+", "."])
        self.assertEqual(e[2], "exon")
        self.assertIn('gene_id "O3cand_000007"', t[8])
        self.assertIn('transcript_id "O3cand_000007.t1"', t[8])
        self.assertIn('placement "requires_WGS"', t[8])
        self.assertIn('nearest_ref_locus "NC_073236.2:95963505-96028675(+)"', t[8])
        self.assertIn('ref_identity_x_coverage "0.9731"', t[8])
        self.assertIn('locus_read_share "0.4812"', t[8])
        self.assertIn('exon_number "1"', e[8])

    def test_missing_values_are_written_as_na(self):
        c = dict(CAND, R=None, nearest=None, read_share=None, cls="NOVEL")
        t = X.gtf_lines(c)[0].split("\t")[8]
        self.assertIn('ref_identity_x_coverage "NA"', t)
        self.assertIn('nearest_ref_locus "NA"', t)
        self.assertIn('locus_read_share "NA"', t)
        self.assertIn('o3_class "NOVEL"', t)

    def test_score_column_is_the_read_count_capped(self):
        self.assertEqual(X.gtf_lines(dict(CAND, reads=5000))[0].split("\t")[5], "1000")


class Bed(unittest.TestCase):
    def test_nearest_locus_bed(self):
        self.assertEqual(X.bed_line(CAND), "NC_073236.2\t95963504\t96028675\tO3cand_000007\t12\t+")
        self.assertIsNone(X.bed_line(dict(CAND, nearest=None)))


class Fasta(unittest.TestCase):
    def test_fasta_header_carries_the_summary(self):
        h, s = X.fasta_record(CAND).split("\n")[:2]
        self.assertEqual(h, ">O3cand_000007 class=DIVERGED reads=12 ref_idcov=0.9731 nearest=NC_073236.2:95963505-96028675(+) placement=requires_WGS")
        self.assertEqual(s, "ACGTACGTAA")


class Select(unittest.TestCase):
    def test_flagged_classes_only_and_stable_ids(self):
        rows = [dict(k="cl3", cls="PRESENT"), dict(k="cl1", cls="DIVERGED"), dict(k="cl2", cls="NOVEL"), dict(k="cl9", cls="UNSUPPORTED"),
                dict(k="cl4", cls="ALLELE-LIKE")]
        self.assertEqual([r["k"] for r in X.flagged(rows)], ["cl1", "cl2"])
        self.assertEqual(X.candidate_id(7), "O3cand_000007")


if __name__ == "__main__":
    unittest.main()
