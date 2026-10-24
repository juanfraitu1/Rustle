"""Tests for o2_moved_loci.py (docs/archive/2026-10/merged/PREREG_o2_population_2026-10-12.md, amendment 1)."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import o2_moved_loci as ml  # noqa: E402


def test_count_moved_reads_per_gene_and_unannotated():
    index = ml.build_index({("c", "g1"): [(100, 200), (300, 400)], ("c", "g2"): [(1000, 1100)]})
    spans = [("c", 150, 350), ("c", 120, 130), ("c", 1050, 1060), ("c", 5000, 5100), ("d", 100, 200)]
    per, unannotated = ml.count_per_gene(spans, index)
    assert per == {("c", "g1"): 2, ("c", "g2"): 1} and unannotated == 2


def test_a_span_between_exons_does_not_hit_the_gene():
    index = ml.build_index({("c", "g1"): [(100, 200), (300, 400)]})
    per, unannotated = ml.count_per_gene([("c", 210, 290)], index)
    assert per == {} and unannotated == 1
