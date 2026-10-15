import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import family_best as FB  # noqa: E402

SPANS = [("chr1", 1_000, 2_000), ("chr2", 5_000, 6_000)]


def test_best_loci_keep_every_record_at_the_top_score():
    recs = [(True, 900, "chr9", [(100, 400)]), (False, 900, "chr1", [(1_900, 2_300)]), (False, 700, "chr2", [(5_100, 5_200)])]
    assert FB.best_loci(recs) == [("chr9", [(100, 400)]), ("chr1", [(1_900, 2_300)])]
    assert FB.in_family(FB.best_loci(recs), SPANS)


def test_a_lower_scoring_family_alignment_does_not_count():
    recs = [(False, 700, "chr1", [(1_100, 1_400)]), (True, 900, "chr9", [(100, 400)])]
    assert not FB.in_family(FB.best_loci(recs), SPANS)


def test_without_scores_the_primary_decides_and_an_intron_is_no_aligned_base():
    assert FB.best_loci([(False, None, "chr1", [(1_100, 1_400)]), (True, None, "chr9", [(1, 2)])]) == [("chr9", [(1, 2)])]
    assert not FB.in_family([("chr1", [(500, 900), (2_100, 2_400)])], SPANS)
    assert FB.best_loci([]) == []
