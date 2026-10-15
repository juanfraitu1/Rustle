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
