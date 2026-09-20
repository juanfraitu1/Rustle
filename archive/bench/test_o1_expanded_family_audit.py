#!/usr/bin/env python3
"""Policy-level tests for the expanded O1 family audit."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import o1_expanded_family_audit as audit


class ExpandedFamilyAuditTest(unittest.TestCase):
    def test_family_prefix_and_loc_placeholder_typing(self) -> None:
        row = {"all_pc_genes": "NBPF10,NOTCH2NLB,LOC124905564"}
        self.assertTrue(audit.matches_family("NBPF10", ("NBPF",)))
        self.assertEqual(
            audit.informative_other_genes(row, ("NBPF",)),
            {"NOTCH2NLB"},
        )

    def test_rgpd_accepts_ranbp2_but_not_plasminogen(self) -> None:
        row = {"all_pc_genes": "RANBP2,RGPD5,PLGLB1"}
        self.assertEqual(
            audit.informative_other_genes(row, ("RANBP2", "RGPD")),
            {"PLGLB1"},
        )

    def test_conflicting_node_cannot_bridge_an_untyped_candidate(self) -> None:
        edges = {audit.pair("known", "conflict"), audit.pair("conflict", "unknown")}
        conflict = {"conflict"}
        clean = {edge for edge in edges if not (set(edge) & conflict)}
        self.assertEqual(audit.reachable({"known"}, clean), {"known"})

    def test_golga2_is_a_related_outgroup_not_swi5_overlap(self) -> None:
        row = {"best_pc_gene": "GOLGA2", "all_pc_genes": "GOLGA2,SWI5"}
        self.assertEqual(audit.best_gene(row), "GOLGA2")
        self.assertIn(
            audit.best_gene(row),
            audit.RELATED_OUTGROUPS[("HSA", "GOLGA6_8")],
        )


if __name__ == "__main__":
    unittest.main()
