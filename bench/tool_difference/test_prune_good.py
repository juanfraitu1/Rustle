"""Tests for prune_good.py (docs/PREREG_o2_pruning_and_eichler_2026-10-12.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import prune_good  # noqa: E402


def test_drops_only_good_pool_only_transcripts_touching_a_contradicted_gene():
    # transcripts: id -> (chrom, exons)
    good = {"g1": ("c", [(100, 200), (300, 400)]),      # good-only chain, touches the gene -> dropped
            "g2": ("c", [(100, 200), (500, 600)]),      # chain also in the primary pool -> kept
            "g3": ("c", [(5000, 5100), (5200, 5300)]),  # good-only but far from the gene -> kept
            "g4": ("c", [(250, 380)]),                  # mono-exon, touches the gene, no primary exon overlaps -> dropped
            "g5": ("c", [(520, 580)])}                  # mono-exon touching the gene but overlapped by a primary exon -> kept
    primary = {"p1": ("c", [(90, 200), (500, 610)])}
    contradicted_exons = {"c": [(120, 180), (250, 300), (530, 540)]}
    assert prune_good.to_drop(good, primary, contradicted_exons) == {"g1", "g4"}


def test_no_contradicted_genes_drops_nothing():
    good = {"g1": ("c", [(100, 200), (300, 400)])}
    assert prune_good.to_drop(good, {}, {}) == set()
