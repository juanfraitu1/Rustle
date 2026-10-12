"""Tests for s1_o2.py (docs/PREREG_s1_o2_certificate_2026-10-12.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import s1_o2  # noqa: E402

EQ, X, I, D, N, S = 7, 8, 1, 2, 3, 4  # pysam cigar op codes


def test_columns_forward_walk_skips_insertions_deletions_and_introns():
    # 2S 3= 1I 2X 1D 2= 5N 1=  (query length 2+3+1+2+2+1 = 11)
    cols = s1_o2.aligned_columns(100, [(S, 2), (EQ, 3), (I, 1), (X, 2), (D, 1), (EQ, 2), (N, 5), (EQ, 1)], reverse=False)
    assert cols == {2: (100, True), 3: (101, True), 4: (102, True), 6: (103, False), 7: (104, False),
                    8: (106, True), 9: (107, True), 10: (113, True)}


def test_columns_reverse_index_from_the_read_start():
    # read length 4, aligned reverse: aligned query position q -> original position 3 - q
    cols = s1_o2.aligned_columns(50, [(EQ, 2), (X, 1), (S, 1)], reverse=True)
    assert cols == {3: (50, True), 2: (51, True), 1: (52, False)}


def test_certify_counts_only_columns_where_the_placements_differ():
    # o -> (ref pos, matched): A matches at 0..4, B mismatches at 0..4 where the bases differ -> 5 columns for A
    a = {o: (o, True) for o in range(10)}
    b = {o: (o, o >= 5) for o in range(10)}
    base_a = lambda o: "A"
    base_b = lambda o: "C" if o < 5 else "A"
    verdict, ka, kb, m, p = s1_o2.certify_pair(a, b, base_a, base_b)
    assert (verdict, ka, kb, m) == ("A", 5, 0, 5) and p < 1e-3


def test_certify_identical_copies_abstain():
    a = {o: (o, True) for o in range(20)}
    verdict, ka, kb, m, p = s1_o2.certify_pair(a, dict(a), lambda o: "G", lambda o: "G")
    assert (verdict, ka, kb, m, p) == (None, 0, 0, 0, 1.0)


def test_certify_one_column_never_assigns():
    a = {0: (0, True)}
    b = {0: (0, False)}
    verdict, ka, kb, m, p = s1_o2.certify_pair(a, b, lambda o: "A", lambda o: "T")
    assert verdict is None and (ka, kb, m) == (1, 0, 1)


def test_certify_mismatch_in_both_counts_as_a_column_for_neither():
    a = {0: (0, True), 1: (1, True), 2: (2, False)}
    b = {0: (0, False), 1: (1, False), 2: (2, False)}
    verdict, ka, kb, m, p = s1_o2.certify_pair(a, b, lambda o: "A", lambda o: "T")
    assert (ka, kb, m) == (2, 0, 3)


def test_gene_class():
    assert s1_o2.gene_class(2, 0) == "supported"
    assert s1_o2.gene_class(2, 9) == "supported"
    assert s1_o2.gene_class(0, 2) == "contradicted"
    assert s1_o2.gene_class(1, 5) == "undecidable"
    assert s1_o2.gene_class(0, 1) == "undecidable"


def test_eichler_assigns_the_unique_best_by_ten_and_discards_closer_rivals():
    # (AS of copy B, primary's AS, molecule best, molecule second)
    assert s1_o2.eichler_verdict(2000, 1980, 2000, 1980) == "B"
    assert s1_o2.eichler_verdict(1980, 2000, 2000, 1980) == "A"
    assert s1_o2.eichler_verdict(1990, 1985, 2000, 1990) == "other"       # margin exactly 10 assigns, to a third placement here
    assert s1_o2.eichler_verdict(1990, 1985, 2000, 1991) == "discard"     # margin 9
    assert s1_o2.eichler_verdict(2000, 2000, 2000, 2000) == "discard"     # a true tie


def test_eichler_best_elsewhere_is_other():
    assert s1_o2.eichler_verdict(1950, 1960, 2000, 1960) == "other"
