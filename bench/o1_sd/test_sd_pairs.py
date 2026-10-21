"""PREREG_o1_families_sd_pairs_2026-10-21: copies of one family paired by an independent SD row."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o1_sd_pairs", os.path.join(HERE, "sd_pairs.py"))
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


def _c(fam, i, chrom, exons):
    return dict(key=(fam, i), family=fam, chrom=chrom, exons=list(exons))


def test_two_copies_of_one_family_on_the_two_sides_of_an_sd_are_paired_but_not_across_families():
    copies = [_c("F", 0, "c1", [(100, 200)]), _c("F", 1, "c2", [(1000, 1100)]), _c("G", 0, "c3", [(5000, 5100)])]
    sd = [("c1", 150, 400, "c2", 900, 1050), ("c1", 150, 400, "c3", 5000, 5050)]
    assert S.sd_paired(copies, sd) == {("F", 0), ("F", 1)}


def test_an_sd_inside_an_intron_pairs_nothing():
    copies = [_c("F", 0, "c1", [(100, 200), (900, 1000)]), _c("F", 1, "c1", [(5000, 5100)])]
    assert S.sd_paired(copies, [("c1", 300, 800, "c1", 5000, 5100)]) == set()


def test_both_sides_on_one_copy_do_not_pair_it_with_itself():
    copies = [_c("F", 0, "c1", [(100, 200), (900, 1000)]), _c("F", 1, "c9", [(1, 2)])]
    assert S.sd_paired(copies, [("c1", 100, 200, "c1", 900, 1000)]) == set()


def test_shuffling_family_labels_keeps_family_sizes():
    copies = [_c(f, i, "c", [(i, i + 1)]) for f, n in (("A", 3), ("B", 2)) for i in range(n)]
    sh = S.shuffle_families(copies, seed=1)
    assert sorted(c["family"] for c in sh) == ["A", "A", "A", "B", "B"] and [c["key"] for c in sh] == [c["key"] for c in copies]


def test_identity_bins():
    assert [S.identity_bin(x) for x in (0.99, 0.98, 0.97, 0.95, 0.92, 0.90, 0.89, None)] == [
        ">= 0.98", ">= 0.98", "0.95-0.98", "0.95-0.98", "0.90-0.95", "0.90-0.95", "< 0.90", "n/a"]
