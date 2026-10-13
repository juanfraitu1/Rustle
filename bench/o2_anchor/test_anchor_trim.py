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


import score as S


def test_distinguishable_needs_two_differing_covered_columns_against_every_rival():
    cols = [(1, "A", "AC"), (2, "G", "GT"), (3, ".", "CA")]
    assert S.distinguishable([5, 7], cols, {5})
    assert not S.distinguishable([5, 7], cols[:1] + cols[2:], {5}), "one covered column only"
    assert not S.distinguishable([5, 7], cols, {9}), "truth not a candidate"


def test_a_nested_truth_set_does_not_count_its_members_as_rivals():
    cols = [(1, "A", "AAC"), (2, "G", "GGT")]
    assert S.distinguishable([5, 6, 7], cols, {5, 6})


def test_overlap_is_same_chrom_and_half_open():
    assert S.overlaps(("c", 0, 10), ("c", 9, 20))
    assert not S.overlaps(("c", 0, 10), ("c", 10, 20))
    assert not S.overlaps(("c", 0, 10), ("d", 0, 10))


def test_upper95_rule_of_three_at_zero_wrong():
    assert S.upper95(0, 100) == round(1 - 0.05 ** (1 / 100), 4)


def test_amended_a1_uses_the_truth_noise_as_the_bar():
    import collections
    c = collections.Counter(t0_n=1000, t0_right=985, wrong=2, wrong_own_sibling=2)
    r = S.amended_a1(c, 100)
    assert r["d"] == 0.015 and r["bar"] == 0.025 and r["verdict"] == "PASS"
    c["wrong"] = 3
    assert S.amended_a1(c, 100)["verdict"] == "FAIL"
    assert S.amended_a1(collections.Counter(t0_n=100, t0_right=90), 100)["verdict"].startswith("STOP")
    assert S.amended_a1(c, 29)["verdict"] == "UNDERPOWERED"
