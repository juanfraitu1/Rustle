"""Unit tests for the population table (run: python3 -m pytest -q bench/popfix)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import table as T


def test_alignment_verdict_overrides_kmer_verdict():
    assert T.combine({"status": "FIXED"}, {"status": "POLYMORPHIC"}) == ("POLYMORPHIC", "A")
    assert T.combine({"status": "POLYMORPHIC"}, {"status": "FIXED"}) == ("FIXED", "A")


def test_kmer_verdict_stands_without_alignment():
    assert T.combine({"status": "FIXED"}, {"status": "NO_ALIGN"}) == ("FIXED", "K")
    assert T.combine({"status": "FIXED"}, None) == ("FIXED", "K")
    assert T.combine({"status": "NOT_ASSESSED"}, None) == ("NOT_ASSESSED", "-")


def test_alignment_uncovered_falls_back_to_kmer():
    assert T.combine({"status": "FIXED"}, {"status": "UNCOVERED"}) == ("FIXED", "K")


# ---- merge_sites: one (chrom, pos) site that belongs to two columns ----

def _km(status, n_carry="40"):
    return {"status": status, "n_men": "48", "n_cover": "40", "n_carry_all": n_carry}


def _two_column_inputs(first_status, second_status):
    """Columns F1/0 and F1/1 both contain the site chrY:100. Returns (cols, km) with the given column order."""
    cols = {("F1", 0): [(0, "chrY", 100, "A", "+")],
            ("F1", 1): [(0, "chrY", 100, "G", "+")]}
    km = {("F1", 0): _km(first_status, n_carry="41"),
          ("F1", 1): _km(second_status, n_carry="37")}
    return cols, km


def _both_orders(first_status, second_status):
    cols, km = _two_column_inputs(first_status, second_status)
    forward = dict(cols)
    reverse = dict(reversed(list(cols.items())))
    return [(forward, km), (reverse, km)]


def test_polymorphic_column_wins_over_fixed_column_in_either_order():
    for cols, km in _both_orders("FIXED", "POLYMORPHIC"):
        merged = T.merge_sites(cols, km, {})
        status, method, column, *_ = merged[("chrY", 100)]
        assert status == "POLYMORPHIC"
        assert column == ("F1", 1)


def test_not_assessed_column_wins_over_fixed_column_in_either_order():
    for cols, km in _both_orders("FIXED", "NOT_ASSESSED"):
        merged = T.merge_sites(cols, km, {})
        status, method, column, *_ = merged[("chrY", 100)]
        assert status == "NOT_ASSESSED"
        assert column == ("F1", 1)


def test_merged_row_carries_kmer_counts_and_alignment_mismatch_of_winning_column():
    cols, km = _two_column_inputs("FIXED", "POLYMORPHIC")
    al = {("F1", 1): {"status": "POLYMORPHIC", "n_mismatch": "3"}}
    merged = T.merge_sites(cols, km, al)
    status, method, column, n_men, n_cover, n_carry, n_mismatch = merged[("chrY", 100)]
    assert (status, method, column) == ("POLYMORPHIC", "A", ("F1", 1))
    assert (n_men, n_cover, n_carry, n_mismatch) == ("48", "40", "37", "3")
