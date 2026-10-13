import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import anchor_trim as A


def hit(t, qs, qe, nm):
    return {"tname": t, "qs": qs, "qe": qe, "nmatch": nm}


def test_shared_interval_is_the_overlap_of_the_truth_hit_and_the_best_sibling_hit():
    hits = [hit("0", 0, 2000, 1990), hit("1", 400, 1800, 1390), hit("2", 300, 900, 590)]
    assert A.shared_interval(hits, "0", min_shared=300) == (400, 1800, "1")


def test_no_sibling_hit_means_no_test_case():
    assert A.shared_interval([hit("0", 0, 2000, 1990)], "0", min_shared=300) is None


def test_a_short_overlap_is_not_a_test_case():
    hits = [hit("0", 0, 2000, 1990), hit("1", 1800, 2000, 199)]
    assert A.shared_interval(hits, "0", min_shared=300) is None


def test_the_truth_copy_must_have_a_hit():
    assert A.shared_interval([hit("1", 0, 2000, 1990)], "0", min_shared=300) is None


def test_truth_copy_is_the_copy_with_the_largest_overlap():
    copies = [("0", "c1", 100, 1000), ("1", "c1", 900, 3000), ("2", "c2", 0, 5000)]
    assert A.truth_copy(copies, "c1", 950, 2500) == "1"
    assert A.truth_copy(copies, "c1", 5000, 6000) is None
    assert A.truth_copy(copies, "c3", 0, 10) is None


def test_original_orientation_reverses_minus_strand_records():
    assert A.original_seq("AACG", 16) == "CGTT"
    assert A.original_seq("AACG", 0) == "AACG"


def test_aligned_ref_span_from_cigar():
    assert A.ref_end(100, "10S50M200N30M5I20M") == 100 + 50 + 200 + 30 + 20
