import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import advisor_table as A  # noqa: E402


def test_row_reads_a_hier_match_json(tmp_path):
    d = {"coverage": {"labels": 10, "labels_evaluable": 6, "labels_one_member": 1, "labels_no_member": 3,
                      "records_on_members": 20, "records_on_unclustered_locus": 4, "records_without_locus": 6},
         "hierarchy": {"exact": 3, "exact_rate": 0.5, "mean_jaccard": 0.8, "macro_recall": 0.9, "macro_precision": 0.85,
                       "macro_gene_recall": 0.7, "by_level": {"family": 4, "expansion": 2}},
         "flat": {"exact": 2, "exact_rate": 1 / 3, "mean_jaccard": 0.7, "macro_recall": 0.75, "macro_precision": 0.9},
         "exact_any": {"expansion": 3}, "null": {"global": {"p_exact": 0.000999}, "within_family": {"p_exact": 0.02}}}
    r = A.row("s", d)
    assert (r["evaluable"], r["hier_exact"], r["hier_exact_of_all"], r["lands_expansion"], r["lands_superfamily"]) == (6, 3, 0.3, 2, 0)
    assert (r["flat_exact_rate"], r["p_global"], r["p_within_family"]) == (0.333, 0.001, 0.02)
