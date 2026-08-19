#!/usr/bin/env python3
"""Unit tests for the single-outgroup rooting proof of concept."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import o1_outgroup_rooting_poc as rooting


def paf(
    query: str,
    query_len: int,
    query_start: int,
    query_end: int,
    strand: str,
    target: str,
    target_len: int,
    target_start: int,
    target_end: int,
    identity: float = 0.98,
    mapq: int = 60,
) -> rooting.PafRecord:
    span = query_end - query_start
    return rooting.PafRecord(
        query, query_len, query_start, query_end, strand, target, target_len,
        target_start, target_end, int(span * identity), span, mapq, identity,
    )


class OutgroupRootingPocTest(unittest.TestCase):
    def test_parse_paf_uses_de_tag(self) -> None:
        record = rooting.parse_paf_line(
            "q\t100\t0\t90\t+\tt\t200\t10\t100\t80\t90\t60\tde:f:0.025"
        )
        self.assertAlmostEqual(record.identity, 0.975)
        self.assertAlmostEqual(record.query_coverage, 0.9)

    def test_collinear_reverse_split_anchor_is_chained(self) -> None:
        records = [
            paf("RIGHT_SRC", 25000, 0, 13294, "-", "chr2", 1000000, 50000, 63294),
            paf("RIGHT_SRC", 25000, 13494, 25000, "-", "chr2", 1000000, 38300, 49800),
        ]
        anchors = rooting.chain_anchor_records(records)
        self.assertEqual(len(anchors), 1)
        self.assertEqual(anchors[0].record_count, 2)
        self.assertAlmostEqual(anchors[0].query_coverage, 24800 / 25000)
        self.assertTrue(rooting.passing_flank(anchors[0]))

    def test_distant_split_records_are_not_chained(self) -> None:
        records = [
            paf("RIGHT_SRC", 25000, 0, 10000, "+", "chr2", 1000000, 10000, 20000),
            paf("RIGHT_SRC", 25000, 12000, 25000, "+", "chr2", 1000000, 22000, 35000),
        ]
        self.assertEqual(len(rooting.chain_anchor_records(records)), 2)

    def test_two_sided_reverse_synteny_returns_inner_interval(self) -> None:
        left = rooting.chain_anchor_records([
            paf("LEFT_SRC", 25000, 0, 25000, "-", "chr2", 1000000, 200000, 225000)
        ])
        right = rooting.chain_anchor_records([
            paf("RIGHT_SRC", 25000, 0, 25000, "-", "chr2", 1000000, 100000, 125000)
        ])
        pairs = rooting.synteny_pairs(left, right)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0][2:], (125000, 200000))

    def test_index_and_fasta_aliases_match_by_sequence_length(self) -> None:
        locus = paf("SRC", 20000, 0, 20000, "+", "chr13_mat_hsa9", 137287480, 100, 20100)
        anchor = rooting.AnchorHit(
            "LEFT_SRC", 25000, "+", "CM054595.2", 137287480,
            50, 100, 0.98, 1.0, 60, 1,
        )
        self.assertTrue(rooting.overlaps(locus, anchor, 90, 20200))


if __name__ == "__main__":
    unittest.main()
