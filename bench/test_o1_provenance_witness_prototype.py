#!/usr/bin/env python3
"""Unit tests for the typed O1 provenance-witness prototype."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import o1_provenance_witness_prototype as provenance


class ProvenanceWitnessPrototypeTest(unittest.TestCase):
    def test_paf_identity_prefers_de_tag_and_coverage_uses_shorter_sequence(self) -> None:
        record = provenance.parse_paf_line(
            "0\t100\t10\t90\t+\t1\t200\t20\t100\t70\t80\t60\tde:f:0.1000"
        )
        self.assertAlmostEqual(record.identity, 0.9)
        self.assertAlmostEqual(record.coverage, 0.8)

    def test_paf_coverage_uses_target_when_target_is_shorter(self) -> None:
        record = provenance.parse_paf_line(
            "0\t200\t20\t140\t+\t1\t100\t10\t70\t50\t60\t60"
        )
        self.assertAlmostEqual(record.coverage, 0.6)

    def test_witness_selection_applies_forward_and_threshold_rules(self) -> None:
        reverse = provenance.parse_paf_line(
            "0\t100\t0\t100\t-\t1\t100\t0\t100\t95\t100\t60"
        )
        forward = provenance.parse_paf_line(
            "0\t100\t0\t80\t+\t1\t100\t0\t80\t72\t80\t60"
        )
        witness = provenance.passing_witness([reverse, forward], 0.6, 0.5, True)
        self.assertEqual(witness, forward)
        self.assertIsNone(provenance.passing_witness([reverse], 0.6, 0.5, True))

    def test_components_include_isolated_loci(self) -> None:
        sizes = provenance.component_sizes({"a", "b", "c", "d"}, {("a", "b"), ("b", "c")})
        self.assertEqual(sizes, [3, 1])

    def test_false_positive_roles_are_not_family_members(self) -> None:
        self.assertEqual(
            provenance.membership_role("REVIEW_RELATED_OUTGROUP"), "CONTEXT_OUTGROUP"
        )
        self.assertEqual(
            provenance.membership_role("REVIEW_CONFLICTING_GENE"), "EXCLUDE_CONFOUND"
        )


if __name__ == "__main__":
    unittest.main()
