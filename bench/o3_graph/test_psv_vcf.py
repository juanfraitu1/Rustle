"""Tests for psv_vcf.py's pure parts: O2's variable columns as VCF rows, and the comparison with another caller's calls."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import psv_vcf as V  # noqa: E402


def obs_two_alleles():
    # column 100: 10 A (ref) / 8 G at Q30; column 200: 12 C only; column 300: 11 T / 1 A (one carrier: not an allele)
    obs = {100: [(f"a{i}", "A", 30) for i in range(10)] + [(f"g{i}", "G", 30) for i in range(8)],
           200: [(f"c{i}", "C", 30) for i in range(12)],
           300: [(f"t{i}", "T", 30) for i in range(11)] + [("x", "A", 30)]}
    return obs


def test_vcf_rows_only_variable_columns_with_poisson_qual():
    ref = {100: "A", 200: "C", 300: "T"}
    rows = V.vcf_rows("chr1", "F", "c0", obs_two_alleles(), lambda col: ref[col])
    assert [r["pos"] for r in rows] == [101]                 # 1-based, columns 200 / 300 are not variable
    r = rows[0]
    assert r["id"] == "F:c0:101" and r["ref"] == "A" and r["alt"] == ["G"]
    assert r["dp"] == 18 and r["ad"] == {"A": 10, "G": 8} and r["m"] == 1
    assert 0 < r["p"] < 1e-6 and r["qual"] > 60            # 8 G against lambda = 18 * 1e-3 / 3


def test_vcf_rows_ref_not_among_alleles_lists_both_as_alt():
    obs = {5: [(f"a{i}", "A", 30) for i in range(6)] + [(f"g{i}", "G", 30) for i in range(6)]}
    rows = V.vcf_rows("chr1", "F", "c0", obs, lambda col: "T")
    assert rows[0]["alt"] == ["A", "G"]


def test_compare_counts_shared_and_caller_only_with_reasons():
    o2 = {("chr1", 101): {"G"}, ("chr1", 501): {"C"}}
    calls = {("chr1", 101): ("A", ["G"], "snv"), ("chr1", 201): ("C", ["T"], "snv"), ("chr1", 301): ("CA", ["C"], "indel")}
    obs = {200: [(f"c{i}", "C", 30) for i in range(12)] + [("y", "T", 30)]}
    out = V.compare(o2, calls, lambda chrom, col: obs.get(col, []))
    assert out["o2_cols"] == 2 and out["o2_called"] == 1 and out["o2_missed"] == [("chr1", 501)]
    assert out["caller_only"] == [(("chr1", 201), "one_allele")]
    assert out["caller_indels"] == 1


def test_reason_names_the_failed_clause():
    assert V.reason([("r", "A", 30)] * 10 + [("s", "G", 30)]) == "one_allele"
    assert V.reason([("r", "A", 30)] * 10 + [("s", "G", 5), ("t", "G", 5)]) == "one_allele"      # Q < 10 bases are missing
    assert V.reason([("r", "A", 30)] * 1000 + [("s", "G", 10), ("t", "G", 10)]).startswith("p=")   # error explains 2 of 1,000
    assert V.reason([]) == "no_pile"
    xs = [("r", "A", 30)] * 30 + [("s", "G", 30), ("t", "G", 30), ("u", "G", 30)]   # p ~ 1e-6: variable alone, not among 10,000 columns
    assert V.reason(xs) == "variable" and V.reason(xs, m=10 ** 5) == "bonferroni"
