"""Tests for score.py (docs/PREREG_asj_o2_2026-10-10.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as S  # noqa: E402


def _row(gene, testsA=(), testsR=(), both=0, disc=0, anchor=None):
    return dict(gene=gene, chrom="c", testsA=list(testsA), testsR=list(testsR), T1=dict(both=both, discordant=disc), anchor=anchor,
                armA=dict(S=1, B=1, abstain=2), n_reads=4)


def test_calls_apply_bh_over_the_arm_and_both_bars():
    t = lambda p, d: dict(intron=[1, 2], p=p, dpsi=d, u0=0, s0=1, u1=1, s1=1)  # noqa: E731
    rows = [_row("g1", testsA=[t(0.001, 0.5), t(0.001, 0.1)]), _row("g2", testsA=[t(0.9, 0.9)])]
    c = S.calls(rows, "testsA")
    assert [(x["gene"], x["dpsi"]) for x in c] == [("g1", 0.5)]


def test_t1_discordance():
    assert S.t1([_row("a", both=100, disc=1), _row("b", both=100, disc=0)]) == (1, 200)


def test_t2_medians_and_power():
    al = [dict(**{"class": "ALLELIC", "conc": 0.99, "n_conc": 12})] * 10
    na = [dict(**{"class": "NON-ALLELIC", "conc": 0.6, "n_conc": 30})] * 9 + [dict(**{"class": "NON-ALLELIC", "conc": 0.9, "n_conc": 5})]
    rows = [_row(f"g{i}", anchor=a) for i, a in enumerate(al + na)]
    r = S.t2(rows)
    assert r["ALLELIC"] == (10, 0.99) and r["NON-ALLELIC"] == (9, 0.6)
