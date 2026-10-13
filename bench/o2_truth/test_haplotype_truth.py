import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import haplotype_truth as H


def test_sam_hit_aligned_bases_and_end():
    line = "r1\t0\tchr1\t101\t60\t5S10=2X3I20=100N10=\t*\t0\t0\tACGT\t*\tNM:i:5\tAS:i:50"
    (h,) = H.parse_sam([line])["r1"]
    assert h == ("chr1", 100, 100 + 10 + 2 + 20 + 100 + 10, 5, 10 + 2 + 20 + 10)


def test_near_hits_need_a_near_perfect_best_and_keep_everything_within_two_edits():
    a = ("mat", "c1", 0, 100, 1, 1000)
    allele = ("pat", "c1x", 0, 100, 2, 1000)    # the same copy's other allele: kept, resolved at copy level
    far = ("pat", "c2", 0, 100, 4, 1000)        # 3 edits worse: dropped
    assert H.near_hits([a, allele, far]) == [a, allele]
    bad = ("mat", "c1", 0, 100, 9, 1000)        # 9 edits > max(2, 3) -> no truth
    assert H.near_hits([bad]) is None
    assert H.near_hits([]) is None


def test_allelic_hits_of_one_copy_give_one_truth_copy():
    loci = {("F", "0"): {("mat", "c1", 0, 100), ("pat", "c1x", 0, 100)}, ("F", "1"): {("mat", "c2", 0, 100)}}
    near = [("mat", "c1", 10, 90, 1, 80), ("pat", "c1x", 10, 90, 2, 80)]
    assert H.truth_copies(near, loci) == {("F", "0")}


def test_copy_loci_keep_hits_within_two_edits_of_the_copy_best_per_haplotype():
    hits = [("mat", "c1", 0, 10, 0, 500), ("mat", "c9", 0, 10, 2, 500), ("mat", "c8", 0, 10, 3, 500),
            ("pat", "c3", 0, 10, 15, 500)]
    assert H.copy_loci(hits) == {("mat", "c1", 0, 10), ("mat", "c9", 0, 10), ("pat", "c3", 0, 10)}


def test_truth_is_the_copies_whose_loci_overlap_the_best_hit():
    loci = {("F", "0"): {("pat", "c1", 100, 200)}, ("F", "1"): {("pat", "c1", 500, 600)}, ("G", "0"): {("mat", "c2", 0, 9)}}
    assert H.truth_copies([("pat", "c1", 150, 180, 0, 30)], loci) == {("F", "0")}
    assert H.truth_copies([("mat", "c3", 0, 10, 0, 10)], loci) == set()
