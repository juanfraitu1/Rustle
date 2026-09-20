#!/usr/bin/env python3
"""Small policy regressions for the typed O1 gene-family audit."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import o1_gene_family_audit as audit


class GeneFamilyAuditTest(unittest.TestCase):
    def test_forward_paf_edge_is_typed_separately_from_reverse(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            paf = Path(directory) / "x.paf"
            paf.write_text(
                "a\t100\t0\t80\t-\tb\t100\t0\t80\t75\t80\t60\tde:f:0.05\n"
                "a\t100\t0\t80\t+\tc\t100\t0\t80\t75\t80\t60\tde:f:0.05\n"
            )
            all_edges, forward = audit.paf_edges(paf)
            self.assertEqual(all_edges, {("a", "b"), ("a", "c")})
            self.assertEqual(forward, {("a", "c")})

    def test_sd_match_does_not_automatically_become_gene_family(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidates = root / "candidates.tsv"
            fields = [
                "locus", "near_family", "n_genomic_paralog_loci", "n_chroms",
                "best_soto_id", "best_soto_cov", "best_soto_family", "class",
            ]
            rows = [
                ["chr1:100-200", "ID_x", "2", "1", "0.95", "0.90", "ID_400", "Soto-family copy"],
                ["chr1:300-400", "ID_x", "2", "1", "0.95", "0.90", "ID_212", "Soto-family copy"],
                ["chr1:500-600", "ID_x", "2", "1", "0.95", "0.90", "ID_179", "Soto-family copy"],
                ["chr1:700-800", "ID_x", "1", "1", "1.00", "1.00", "ID_245", "single-copy"],
            ]
            with candidates.open("w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
                writer.writerow(fields)
                writer.writerows(rows)
            gff = root / "genes.gff"
            gff.write_text(
                "##gff-version 3\n"
                "chr1\ttest\tgene\t101\t200\t.\t+\t.\tID=gene-NBPF8;Name=NBPF8;"
                "description=NBPF member 8;gene_biotype=protein_coding\n"
                "chr1\ttest\tgene\t301\t400\t.\t+\t.\tID=gene-X;Name=X;"
                "description=unrelated;gene_biotype=protein_coding\n"
                "chr1\ttest\tgene\t501\t600\t.\t+\t.\tID=gene-CHRFAM7A;Name=CHRFAM7A;"
                "description=CHRNA7 fusion;gene_biotype=protein_coding\n"
            )
            out = root / "out"
            out.mkdir()
            got = audit.candidate_audit(candidates, gff, out)
            self.assertEqual(got[0]["disposition"], "SUPPORTED_SOTO_MISSING_GENE_COPY")
            self.assertEqual(got[1]["disposition"], "SD_EXTENSION_ONLY")
            self.assertEqual(got[2]["disposition"], "REJECT_GENE_FAMILY_ASSIGNMENT")
            self.assertEqual(got[3]["disposition"], "REJECT_NOT_MULTICOPY")

    def test_named_other_gene_is_review_not_automatic_family_member(self) -> None:
        row = {
            "seed_genes": ".",
            "overlapping_refseq_genes": "LOC100505915;PDXDC1",
        }
        self.assertEqual(audit.conflicting_gene_symbols(row, "NPIP"), {"PDXDC1"})

        same_family = {
            "seed_genes": "NPIPB3",
            "overlapping_refseq_genes": "LOC100190986;NPIPB3",
        }
        self.assertEqual(audit.conflicting_gene_symbols(same_family, "NPIP"), set())


if __name__ == "__main__":
    unittest.main()
