import importlib.util
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2vc_rule", os.path.join(HERE, "o2vc_rule.py"))
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)


def col(base, qual=30, **g):
    """The read's base and quality, each copy's genotype in the read's frame (F_0 -> F.0; None = unknown)."""
    return dict(base=base, qual=qual, g={k.replace("_", "."): None if v is None else frozenset(v) for k, v in g.items()})


def test_a_copy_s_genotype_is_its_record_s_gt_alleles_whatever_the_filter():
    C = R.Call
    assert R.genotype(None, 6, "A") == frozenset("A")
    assert R.genotype(C("A", ("G",), (0, 1), "LowQual", 9), 9, "A") == frozenset("AG")
    assert R.genotype(C("A", ("G",), (1, 0), "PASS", 9), 9, "A") == frozenset("AG")
    assert R.genotype(C("A", ("G",), (1, 1), "PASS", 9), 9, "A") == frozenset("G")
    assert R.genotype(C("A", ("G", "T"), (1, 2), "Multiallelic", 9), 9, "A") == frozenset("GT")
    assert R.genotype(C("A", ("G",), (0, 0), "HomRef", 9), 9, "A") == frozenset("A")


def test_a_copy_is_unknown_below_depth_six_or_at_an_editing_or_dense_record():
    C = R.Call
    assert R.genotype(None, 5, "A") is None
    assert R.genotype(C("A", ("G",), (0, 1), "RnaEdit", 30), 30, "A") is None
    assert R.genotype(C("A", ("G",), (0, 1), "Dense", 30), 30, "A") is None


def test_alleles_are_complemented_into_the_read_s_frame_for_a_copy_aligned_opposite_to_the_read():
    assert R.to_read_frame(frozenset("AG"), True) == frozenset("TC")
    assert R.to_read_frame(frozenset("AG"), False) == frozenset("AG")
    assert R.to_read_frame(None, True) is None


def test_a_copy_is_placed_when_a_separating_column_votes_for_it_and_none_for_any_rival():
    assert R.assign(["F.0", "F.1"], [col("C", F_0="C", F_1="T"), col("A", F_0="A", F_1="A")]) == ("F.0", frozenset({"F.0"}), 1)
    three = [col("C", F_0="C", F_1="T", F_2="C"), col("G", F_0="G", F_1="G", F_2="A")]
    assert R.assign(["F.0", "F.1", "F.2"], three) == ("F.0", frozenset({"F.0"}), 2)


def test_contradicting_columns_and_shared_alleles_leave_the_copy_set():
    assert R.assign(["F.0", "F.1"], [col("C", F_0="C", F_1="T"), col("G", F_0="A", F_1="G")]) == (None, frozenset({"F.0", "F.1"}), 0)
    assert R.assign(["F.0", "F.1"], [col("C", F_0="CT", F_1="T")]) == (None, frozenset({"F.0", "F.1"}), 0)


def test_an_unknown_copy_is_never_beaten_and_a_low_quality_base_never_votes():
    assert R.assign(["F.0", "F.1", "F.2"], [col("C", F_0="C", F_1="T", F_2=None)]) == (None, frozenset({"F.0", "F.2"}), 0)
    assert R.assign(["F.0", "F.1"], [col("C", qual=9, F_0="C", F_1="T")]) == (None, frozenset({"F.0", "F.1"}), 0)


def test_a_foreign_allele_needs_every_copy_known_and_none_holding_the_read_s_base():
    assert R.foreign(col("G", F_0="C", F_1="T")) == "foreign"
    assert R.foreign(col("G", F_0="C", F_1=None)) == ""
    assert R.foreign(col("C", F_0="C", F_1="T")) == ""
    assert R.foreign(col("G", qual=5, F_0="C", F_1="T")) == ""
    assert R.foreign(col("N", F_0="C", F_1="T")) == ""
    assert R.foreign(col("G", F_0="A", F_1="AT")) == "editing"          # A>G in the read's (transcript) frame


def _k(pos, b, copy="F.0", chrom="chr1"):
    return (copy, chrom, pos, b)


def test_a_flag_is_two_reads_sharing_two_linked_foreign_alleles_at_different_columns():
    by_read = {"r1": [(10, frozenset({_k(100, "G"), _k(500, "G", "F.1", "chr2")})), (20, frozenset({_k(110, "T")}))],
               "r2": [(11, frozenset({_k(100, "G")})), (21, frozenset({_k(110, "T")}))],
               "r3": [(12, frozenset({_k(100, "G")}))]}
    f = R.flags(by_read)
    assert [x["reads"] for x in f] == [["r1", "r2"]]
    assert f[0]["alleles"] == [_k(100, "G"), _k(110, "T")] and f[0]["id"] == 0


def test_one_shared_foreign_column_seen_on_two_copies_makes_no_flag():
    same = frozenset({_k(100, "G"), _k(500, "G", "F.1", "chr2")})
    assert R.flags({"r1": [(10, same)], "r2": [(11, same)]}) == []


def test_flags_do_not_depend_on_the_order_reads_arrive():
    by_read = {f"r{i}": [(1, frozenset({_k(100 + i % 2, "G")})), (2, frozenset({_k(200, "T")}))] for i in range(6)}
    keys = list(by_read)
    random.Random(7).shuffle(keys)
    assert R.flags({k: by_read[k] for k in keys}) == R.flags(by_read)


def test_a_record_whose_gt_calls_no_allele_leaves_the_copy_unknown():
    assert R.genotype(R.Call("A", ("G",), (), "PASS", 9), 9, "A") is None
