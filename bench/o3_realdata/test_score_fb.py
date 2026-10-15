import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score_fb as F  # noqa: E402


def cl(fam, host, k, status, n, est=0.02):
    return dict(family=fam, host=str(host), cluster=str(k), status=status, n_reads=str(n), est_div=str(est))


def test_family_counts_copy_and_explained():
    clusters = [cl("A", 0, 0, "explained", 10), cl("A", 0, 1, "unexplained", 5, 0.02),     # A: copy + explained
                cl("B", 0, 0, "unexplained", 5, 0.005),                                     # B: allele_or_copy only
                cl("C", 0, 0, "explained", 8), cl("C", 1, 0, "unlinked", 4, 0.05)]           # C: explained, unlinked
    assert F.family_counts(["A", "B", "C", "D"], clusters) == dict(families=4, copy=1, explained=2)


def world():
    # masked copies of A (10 reads, 8 placed) and B (4 placed: not eligible)
    absent = {"A": ("A:2", {f"a{i}" for i in range(10)}), "B": ("B:1", {"b0", "b1", "b2", "b3"})}
    placed = {"A": {f"a{i}" for i in range(8)} | {"s0"}, "B": {"b0", "b1", "b2", "b3"}}
    return absent, placed


def test_s2_counts_recoveries_on_the_same_eligible_copies_with_the_placed_denominator():
    absent, placed = world()
    old = ([cl("A", 0, 0, "unexplained", 8)], {("A", "0", "0"): {f"a{i}" for i in range(8)}})
    # the new rule keeps only 6 of the 8 placed reads in the cluster: 6 / 8 = 0.75 < 0.8, not recovered
    new = ([cl("A", 0, 0, "unexplained", 6)], {("A", "0", "0"): {f"a{i}" for i in range(6)}})
    r = F.removed_side(absent, placed, placed, old, new)
    assert (r["eligible"], r["recovered_old"], r["recovered_new"]) == (1, 1, 0)


def test_s2_refuses_placed_sets_that_differ_between_rules():
    absent, placed = world()
    other = {k: set(v) for k, v in placed.items()}
    other["A"].discard("a0")
    try:
        F.removed_side(absent, placed, other, ([], {}), ([], {}))
    except ValueError as e:
        assert "placed" in str(e)
    else:
        raise AssertionError("different placed sets must be refused")


def test_bars():
    v = F.bars(fams=40, copy_old=10, copy_new=4, expl_old=30, expl_new=28, eligible=12, rec_old=4, rec_new=3)
    assert v == dict(S1=0.1, S1_pass=True, S1b_pass=True, S2_measurable=True, S2_pass=True, G_pass=True)
    v = F.bars(fams=40, copy_old=10, copy_new=6, expl_old=30, expl_new=26, eligible=4, rec_old=4, rec_new=1)
    assert v == dict(S1=0.15, S1_pass=False, S1b_pass=False, S2_measurable=False, S2_pass=None, G_pass=False)
    assert F.bars(fams=40, copy_old=0, copy_new=0, expl_old=30, expl_new=30, eligible=6, rec_old=2, rec_new=0)["S1b_pass"] is None
