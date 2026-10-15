import collections
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import subset_lib as L  # noqa: E402

# Our blocks: B1 = {a, b, c, d, i}, B2 = {e, f, j}, singletons g and h.
BLOCK_OF = {"a": "B1", "b": "B1", "c": "B1", "d": "B1", "i": "B1", "e": "B2", "f": "B2", "j": "B2",
            "g": "Lg", "h": "Lh"}
# Soto: S1 and S2 contained; S3 needs our singleton g (coverage); S4 shares d with S2 (cover); S5 is a conflict.
SOTO_SETS = {"S1": frozenset("ab"), "S2": frozenset("cd"), "S3": frozenset("eg"), "S4": frozenset("df"),
             "S5": frozenset("ij")}
MULTI = frozenset("d")
SOTO_OF = {"a": "S1", "b": "S1", "c": "S2", "e": "S3", "g": "S3", "f": "S4", "i": "S5", "j": "S5", "h": "U_h"}


# ------------------------------------------------------------------------------------------------------------ T1

def test_containment_count():
    assert L.containment_count(SOTO_SETS, BLOCK_OF) == 2


def test_classify_violation_each_class():
    size = collections.Counter(BLOCK_OF.values())
    assert L.classify_violation(sorted(SOTO_SETS["S3"]), BLOCK_OF, size, MULTI) == "coverage"
    assert L.classify_violation(sorted(SOTO_SETS["S4"]), BLOCK_OF, size, MULTI) == "cover"
    assert L.classify_violation(sorted(SOTO_SETS["S5"]), BLOCK_OF, size, MULTI) == "conflict"


def test_conditional_entropies_by_hand():
    members = sorted(SOTO_OF)
    # H(Soto | ours): B1 holds S1, S1, S2, S5 (1.5 bits, weight 4/9); B2 holds S3, S4, S5 (log2 3, weight 3/9).
    assert abs(L.cond_entropy(SOTO_OF, BLOCK_OF, members) - (2 / 3 + math.log2(3) / 3)) < 1e-12
    # H(ours | Soto): S3 and S5 each split over two blocks (1 bit, weight 2/9 each).
    assert abs(L.cond_entropy(BLOCK_OF, SOTO_OF, members) - 4 / 9) < 1e-12
    assert abs(L.direction(SOTO_OF, BLOCK_OF) - (2 / 3 + math.log2(3) / 3 - 4 / 9)) < 1e-12


def test_p_upper():
    assert L.p_upper([0, 1, 2, 3], 3) == 2 / 5
    assert L.p_upper([], 5) == 1.0
    assert math.isnan(L.p_upper([0.1, 0.2], float("nan")))


def test_size_null_is_seeded_and_order_free():
    a = L.size_null(SOTO_SETS, BLOCK_OF, 50, 7)
    b = L.size_null(SOTO_SETS, dict(reversed(list(BLOCK_OF.items()))), 50, 7)
    assert a == b and len(a) == 50


def test_size_null_probability():
    # Blocks of sizes 2 and 1 over x, y, z: the family {x, y} is contained with probability 1/3.
    null = L.size_null({"S": frozenset("xy")}, {"x": "B", "y": "B", "z": "C"}, 3000, 7)
    assert abs(sum(null) / 3000 - 1 / 3) < 0.035


def test_one_block_never_beats_its_null():
    one = {m: "B" for m in BLOCK_OF}
    obs = L.containment_count(SOTO_SETS, one)
    assert obs == 5 and L.p_upper(L.size_null(SOTO_SETS, one, 20, 7), obs) == 1.0


def test_units():
    assert L.units(SOTO_OF, BLOCK_OF) == [list("abcefgij"), ["h"]]


def test_bootstrap_identical_partitions_is_zero():
    same = {m: "B" + SOTO_OF[m] for m in SOTO_OF}
    assert L.direction_bootstrap(SOTO_OF, same, 50, 7) == (0.0, 0.0)


def test_bootstrap_one_block_three_soto_pairs():
    soto = {"a": "S1", "b": "S1", "c": "S2", "d": "S2", "e": "S3", "f": "S3"}
    lo, hi = L.direction_bootstrap(soto, {m: "B" for m in soto}, 50, 7)
    assert abs(lo - math.log2(3)) < 1e-12 and abs(hi - math.log2(3)) < 1e-12


def test_bootstrap_is_seeded_and_order_free():
    a = L.direction_bootstrap(SOTO_OF, BLOCK_OF, 200, 7)
    b = L.direction_bootstrap(dict(reversed(list(SOTO_OF.items()))), BLOCK_OF, 200, 7)
    assert a == b and a[0] <= L.direction(SOTO_OF, BLOCK_OF) <= a[1]


def test_t1_on_the_toy():
    clustered = frozenset(m for m, b in BLOCK_OF.items() if b in ("B1", "B2"))
    r = L.t1(SOTO_SETS, SOTO_OF, BLOCK_OF, MULTI, clustered, 20, 20, 7)
    assert r["contained"] == 2 and r["families_ge2"] == 5 and r["rate_all"] == 2 / 5
    assert r["violation_classes"] == {"coverage": 1, "cover": 1, "conflict": 1}
    assert r["by_size"] == {"le3": {"families": 5, "contained": 2}, "gt3": {"families": 0, "contained": 0}}
    assert r["plan1"] == {"evaluable": 4, "contained": 2}
    assert [v["label"] for v in r["violations"]] == ["S3", "S4", "S5"]
    assert r["t1a_pass"] == (r["p"] < 0.001) and r["t1b_pass"] == (r["D_ci95"][0] > 0)


def test_t1_without_any_family():
    r = L.t1({"S": frozenset("a")}, {"a": "S"}, {"a": "B", "b": "C"}, frozenset(), frozenset("ab"), 10, 10, 7)
    assert r["contained"] == 0 and r["families_ge2"] == 0 and r["p"] == 1.0
    assert not r["t1a_pass"] and not r["t1b_pass"]


# ------------------------------------------------------------------------------------------------------------ T2

def test_cn_components_runs():
    block_of = {"a": "B1", "b": "B1", "c": "B1", "d": "B1", "e": "B2", "f": "B2", "g": "B2"}
    cn = {"a": 2.0, "b": 3.5, "c": 5.0, "d": 9.0, "e": 2.0, "f": 3.0}
    comp = L.cn_components(block_of, cn)
    assert comp["a"] == comp["b"] == comp["c"] != comp["d"]
    assert comp["e"] == comp["f"] and comp["e"] != comp["a"] and "g" not in comp


def test_cn_gap_of_exactly_two_cuts():
    # Soto's rule keeps a pair when MAD = |a - b| / 2 < 1: a difference of exactly 2 cuts, 1.999 joins.
    comp = L.cn_components({"x": "B", "y": "B", "z": "B"}, {"x": 2.0, "y": 4.0, "z": 5.999})
    assert comp["x"] != comp["y"] and comp["y"] == comp["z"]


def test_cn_shuffle_keeps_each_block_multiset():
    block_of = {"a": "B1", "b": "B1", "c": "B1", "d": "B2", "e": "B2"}
    cn = {"a": 1.0, "b": 2.0, "c": 3.0, "d": 7.0, "e": 9.0}
    s1 = L.cn_shuffle(block_of, cn, random.Random(7))
    s2 = L.cn_shuffle(block_of, dict(reversed(list(cn.items()))), random.Random(7))
    assert s1 == s2 and set(s1) == set(cn)
    assert sorted(s1[m] for m in "abc") == [1.0, 2.0, 3.0] and sorted(s1[m] for m in "de") == [7.0, 9.0]


def test_ari_by_hand():
    t = {"a": "S1", "b": "S1", "c": "S2", "d": "S2"}
    assert L.ari(t, t, "abcd") == 1.0
    assert L.ari(t, {m: "one" for m in "abcd"}, "abcd") == 0.0
    assert math.isnan(L.ari(t, t, "a"))


def test_exact_counts_distinct_sets():
    soto = {"S1": frozenset("abc"), "S2": frozenset("abc"), "S3": frozenset("dx")}
    comp = {"a": "c1", "b": "c1", "c": "c1", "d": "c2"}
    assert L.exact_count(soto, comp, "abcd") == 1  # {a, b, c} counted once; S3 keeps one member inside


def test_hungarian_jaccard_by_hand():
    soto = {"S1": frozenset("abc"), "S2": frozenset("de")}
    comp = {"a": "C1", "b": "C1", "c": "C2", "d": "C2", "e": "C2"}
    assert abs(L.hungarian_jaccard(soto, comp, "abcde") - 2 / 3) < 1e-12
    one = {m: "C" for m in "abcd"}
    assert abs(L.hungarian_jaccard({"S1": frozenset("ab"), "S2": frozenset("cd")}, one, "abcd") - 0.25) < 1e-12
    assert math.isnan(L.hungarian_jaccard({"S1": frozenset("a")}, {"a": "C"}, "a"))


def test_meet():
    assert L.meet({"a": "S1", "b": "S1", "c": "S2"}, {"a": "B1", "b": "B2", "c": "B1"}) == {
        "a": "S1|B1", "b": "S1|B2", "c": "S2|B1"}


def test_scores_drops_fam90a_only_for_its_ari():
    soto_of = {"a": "S1", "b": "S1", "c": L.FAM90A, "d": L.FAM90A}
    soto = {"S1": frozenset("ab"), L.FAM90A: frozenset("cd")}
    comp = {"a": "x", "b": "x", "c": "y", "d": "z"}
    s = L.scores(soto, soto_of, comp, "abcd")
    assert s["exact"] == 1 and s["single"] == 4 and s["ari_no_fam90a"] == 1.0 and s["ari"] < 1.0


def test_t2_reproduces_soto_and_null_rate():
    soto = {"S1": frozenset("ab"), "S2": frozenset("cd")}
    soto_of = {"a": "S1", "b": "S1", "c": "S2", "d": "S2"}
    cn = {"a": 2.0, "b": 2.0, "c": 10.0, "d": 10.0}
    r = L.t2(soto, soto_of, {m: "B" for m in "abcd"}, cn, 3000, 7)
    assert r["observed"]["exact"] == 2 and r["observed"]["ari"] == 1.0
    # 2 of the 6 distinct arrangements of {2, 2, 10, 10} rebuild Soto, so p is near 1/3 and the bar fails.
    assert abs(r["p_exact"] - 1 / 3) < 0.035 and not r["passes"]


def test_t2_strong_signal_has_small_p():
    labels = "abcdefghijkl"
    soto_of = {m: f"S{i // 2}" for i, m in enumerate(labels)}
    soto = {f"S{k}": frozenset(labels[2 * k:2 * k + 2]) for k in range(6)}
    cn = {m: 10.0 * (i // 2) for i, m in enumerate(labels)}
    r = L.t2(soto, soto_of, {m: "B" for m in labels}, cn, 50, 7)
    # 720 of 12! / 2^6 = 7,484,400 arrangements rebuild Soto (about 1e-4 per draw).
    assert r["observed"]["exact"] == 6 and r["p_exact"] < 0.05 and r["p_ari"] < 0.05
    assert r["passes"] == (r["p_exact"] < 0.001 and r["p_ari"] < 0.001)


def test_t2_with_one_member_is_not_a_pass():
    r = L.t2({"S1": frozenset("a")}, {"a": "S1"}, {"a": "B"}, {"a": 2.0}, 10, 7)
    assert math.isnan(r["p_ari"]) and not r["passes"]
