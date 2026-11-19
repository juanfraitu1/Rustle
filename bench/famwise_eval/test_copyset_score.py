"""tests of copyset_score.py (docs/GOAL_multicopy_families_2026-11-19.md)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R
import copyset_score as C

LOCI = [("c1", 100, 1000, "+", 1000), ("c1", 5000, 6000, "+", 1000)]
KA, KB = ("c1", "+", ((200, 300),)), ("c1", "+", ((5200, 5300),))

def t(kind, origin):
    return R.Truth("S2", "x", kind, frozenset(origin), LOCI)

def rec(s, e, ch, ms, prim):
    return R.Rec("c1", s, e, "+", ch, ms, ms, prim)

def claim(*spans_keys):
    return dict(keys=frozenset(k for k, _ in spans_keys), spans=[sp for _, sp in spans_keys])

CA = claim((KA, ("c1", "+", 150, 900)))
CB = claim((KB, ("c1", "+", 5150, 5900)))
CAB = claim((KA, ("c1", "+", 150, 900)), (KB, ("c1", "+", 5150, 5900)))

def test_verdicts():
    assert C.read_verdict({0}, {0}) == "EQUAL" and C.read_verdict({0}, {0, 1}) == "NARROW"
    assert C.read_verdict({0, 1}, {0}) == "WIDE" and C.read_verdict({1}, {0}) == "OFF" and C.read_verdict({0, -1}, {1, 0, 2} - {2} | {1}) == "MIXED"

def test_claimed_local_and_regime():
    assert C.claimed_local(CAB, t("several", {0, 1})) == {0, 1}
    far = claim((KA, ("c2", "+", 1, 50)))
    assert C.claimed_local(far, t("one", {0})) == {-1}
    assert C.n_places(CAB["spans"]) == 2 and C.n_places([("c1", "+", 1, 10), ("c1", "+", 5, 20)]) == 1

def test_score_and_recall():
    truth = {"a": t("one", {0}), "b": t("one", {0}), "s1": t("several", {0, 1}), "s2": t("several", {0, 1})}
    recs = {"a": [rec(100, 900, KA[2], 500, True)], "b": [rec(100, 900, KA[2], 500, True)],
            "s1": [rec(100, 900, KA[2], 500, True), rec(5100, 5900, KB[2], 500, False)],
            "s2": [rec(5100, 5900, KB[2], 500, True), rec(100, 900, KA[2], 500, False)]}
    sup, us = C.supporters(truth, recs), C.units(truth, recs, 2)
    assert {u: v[:2] for u, v in us.items()} == {frozenset({KA}): (2, "one"), frozenset({KA, KB}): (2, "set")}
    per_copy = C.score_arm([CA, CB], sup, truth)
    assert per_copy[("one", "MIXED")] == 1 and per_copy[("one", "NARROW")] == 1      # A: 2 EQUAL vs 2 NARROW -> tie -> MIXED; B: NARROW
    fam = C.score_arm([CAB], sup, truth)
    assert fam[("set", "MIXED")] == 1                                                  # 2 EQUAL (s1, s2) vs 2 WIDE (a, b): tie -> MIXED
    assert C.plurality(["EQUAL", "EQUAL", "WIDE"]) == "EQUAL" and fam[("set", "w_judged")] == 4
    r1, r2 = C.recall_arm([CA, CB], us), C.recall_arm([CAB], us)
    assert (r1[("set", "structure")], r1[("set", "resolution")]) == (1, 0)
    assert (r2[("set", "structure")], r2[("set", "resolution")]) == (1, 1) and r2[("one", "resolution")] == 0

def test_containment():
    k = ("c1", "+", ((1, 2), (3, 4), (5, 6)))
    cont = C.Containment({k})
    assert cont.holders("c1", "+", ((3, 4), (5, 6))) == [k] and cont.holders("c1", "+", ((1, 2), (5, 6))) == []
    us = {frozenset({("c1", "+", ((3, 4),))}): (2, "one", [("c1", 0, 10, "+", 1)])}
    cl = dict(keys=frozenset({k}), spans=[("c1", "+", 0, 10)])
    assert C.recall_arm([cl], us)[("one", "structure")] == 0
    r = C.recall_arm([cl], us, contained=True)
    assert r[("one", "structure")] == 1 and r[("one", "resolution")] == 1
