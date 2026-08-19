#!/usr/bin/env python3
"""Unit tests for coordinate matching in the fresh O1 emission check."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import o1_fresh_emission_validation as audit


class FreshEmissionValidationTest(unittest.TestCase):
    def test_best_match_uses_overlap_of_shorter_then_jaccard(self) -> None:
        fresh = [
            {"chrom": "chr1", "start": "90", "end": "210", "family_id": "F0", "copy_idx": "0"},
            {"chrom": "chr1", "start": "100", "end": "200", "family_id": "F1", "copy_idx": "0"},
        ]
        row, containment, jaccard = audit.best_match(("chr1", 100, 200), fresh)
        self.assertEqual(row["family_id"], "F1")
        self.assertEqual(containment, 1.0)
        self.assertEqual(jaccard, 1.0)

    def test_disjoint_locus_does_not_match(self) -> None:
        fresh = [
            {"chrom": "chr1", "start": "300", "end": "400", "family_id": "F0", "copy_idx": "0"},
        ]
        row, containment, jaccard = audit.best_match(("chr1", 100, 200), fresh)
        self.assertIsNone(row)
        self.assertEqual((containment, jaccard), (0.0, 0.0))

    def test_conflict_in_target_has_highest_visual_priority(self) -> None:
        rows = [
            {
                "genes": "NBPF4",
                "case": "NBPF",
                "disposition": "KEEP_KNOWN",
                "in_fresh_target_family": 1,
            },
            {
                "genes": "TTC6",
                "case": "NBPF",
                "disposition": "REVIEW_CONFLICTING_GENE",
                "in_fresh_target_family": 1,
            },
        ]
        status, genes, case = audit.fresh_node_status(rows)
        self.assertEqual(status, "CONFLICT_IN_TARGET")
        self.assertEqual(genes, "NBPF4,TTC6")
        self.assertEqual(case, "NBPF")

    def test_reemitted_conflict_outside_target_is_separated(self) -> None:
        rows = [
            {
                "genes": "DNAH14",
                "case": "NBPF",
                "disposition": "REVIEW_CONFLICTING_GENE",
                "in_fresh_target_family": 0,
            }
        ]
        status, _, _ = audit.fresh_node_status(rows)
        self.assertEqual(status, "CONFLICT_SEPARATED")

    def test_related_outgroup_gets_its_own_visual_class(self) -> None:
        rows = [
            {
                "genes": "GOLGA2,SWI5",
                "case": "GOLGA6_8",
                "disposition": "REVIEW_RELATED_OUTGROUP",
                "in_fresh_target_family": 1,
            }
        ]
        status, _, _ = audit.fresh_node_status(rows)
        self.assertEqual(status, "RELATED_OUTGROUP")


if __name__ == "__main__":
    unittest.main()
