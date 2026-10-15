import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import genomic_explain as GE  # noqa: E402

REF = "ACGTACGTAC"
SITES = {2: "T", 5: "G"}                      # alt bases at positions 2 and 5
CONS = {2: "alt", 5: "alt", 7: "ref"}


def test_a_locus_carrying_every_alt_and_the_reference_elsewhere_explains():
    assert GE.explains(CONS, SITES, REF, {2: "T", 5: "G", 7: "T"})


def test_no_sequence_at_an_alt_is_no_explanation():
    # what a representative lacking the exon looks like: the locus has no base at position 5
    assert not GE.explains(CONS, SITES, REF, {2: "T", 7: "T"})


def test_a_different_base_at_an_alt_or_a_conflict_at_a_reference_site_is_no_explanation():
    assert not GE.explains(CONS, SITES, REF, {2: "T", 5: "C", 7: "T"})
    assert not GE.explains(CONS, SITES, REF, {2: "T", 5: "G", 7: "A"})


def test_merge_spans_joins_nearby_alignments_per_chromosome():
    spans = [("c1", 100, 200), ("c1", 4_000, 5_000), ("c1", 20_000, 21_000), ("c2", 0, 10)]
    assert GE.merge_spans(spans, gap=5_000) == [("c1", 100, 5_000), ("c1", 20_000, 21_000), ("c2", 0, 10)]


def test_cluster_consensus_is_the_majority_base_where_enough_reads_cover():
    covs = [{0: "A", 1: "C", 2: "T"}, {0: "A", 1: "C", 2: "T"}, {0: "A", 1: "G", 2: "T"}, {0: "A", 2: "G"}]
    assert GE.consensus(covs, min_reads=3) == {0: "A", 1: "C", 2: "T"}


def test_divergence_counts_mismatches_over_the_bases_both_cover():
    cons = {0: "A", 1: "C", 2: "T", 3: "G"}
    d, cover = GE.divergence(cons, {0: "A", 1: "G", 2: "T", 9: "C"})
    assert (d, cover) == (1 / 3, 0.75)
    assert GE.divergence(cons, {}) == (1.0, 0.0)
    assert GE.divergence(cons, {0: "-", 1: "C"}) == (0.0, 0.25)   # a deleted base is no comparison
