import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rd_lib as L  # noqa: E402


def cp(fam, idx, chrom, s, e, n_reads=50):
    return dict(family=fam, idx=idx, chrom=chrom, start=s, end=e, exons=[(s, e)], n=n_reads)


def test_pick_masks_the_last_copy_and_applies_every_rule():
    copies = [cp("A", 0, "c1", 100, 200), cp("A", 1, "c1", 1000, 1100), cp("A", 2, "c2", 50, 150),          # kept
              cp("B", 0, "c1", 5000, 5100), cp("B", 1, "c1", 6000, 6100),                                  # 2 copies
              cp("C", 0, "c3", 0, 100), cp("C", 1, "c3", 500, 600), cp("C", 2, "c3", 900, 1000, n_reads=5),  # low reads
              cp("D", 0, "c4", 0, 100), cp("D", 1, "c4", 500, 600), cp("D", 2, "c4", 900, 400_000),         # span
              cp("E", 0, "c5", 0, 100), cp("E", 1, "c5", 500, 600), cp("E", 2, "c5", 900, 1000),            # excluded
              cp("F", 0, "c6", 0, 100), cp("F", 1, "c6", 500, 600), cp("F", 2, "c6", 900, 1000)]            # overlap
    other = [dict(family="Z", idx=0, chrom="c6", start=950, end=1200, exons=[(950, 1200)], n=50)]
    picks = L.pick_families(copies, catalog=copies + other, exclude=[("c5", 0, 2000)], cap=10, seed=1)
    assert [(p["family"], p["masked"]["idx"]) for p in picks] == [("A", 2)]


def test_pick_caps_with_a_seeded_sample():
    copies = [cp(f"F{i}", k, f"c{i}", 1000 * k, 1000 * k + 100) for i in range(6) for k in range(3)]
    a = L.pick_families(copies, catalog=copies, exclude=[], cap=3, seed=1)
    b = L.pick_families(copies, catalog=copies, exclude=[], cap=3, seed=1)
    assert len(a) == 3 and [p["family"] for p in a] == [p["family"] for p in b]


def test_reads_touching_one_copy_only():
    copies = {("A", 0): ("c1", [(100, 200)]), ("A", 1): ("c1", [(150, 300)]), ("B", 0): ("c2", [(0, 100)])}
    alns = [("r1", "c1", [(110, 120)]), ("r2", "c1", [(160, 170)]), ("r3", "c2", [(50, 60)]), ("r4", "c2", [(500, 600)]),
            ("r5", "c1", [(250, 260)])]
    assert L.reads_by_copy(alns, copies) == {"r1": ("A", 0), "r3": ("B", 0), "r5": ("A", 1)}


def test_sample_per_copy_is_seeded_and_capped():
    rb = {f"r{i}": ("A", i % 2) for i in range(10)}
    a, b = L.sample_per_copy(rb, n=3, seed=1), L.sample_per_copy(rb, n=3, seed=1)
    assert a == b and sum(1 for v in a.values() if v == ("A", 0)) == 3


def clusters(*rows):
    return [dict(family=f, host=h, cluster=c, n_reads=str(n), status=st, est_div=str(d)) for f, h, c, n, st, d in rows]


def test_recovered_needs_recall_and_purity_and_labels_by_divergence():
    d = {f"d{i}" for i in range(10)}
    placed = d | {"s1", "s2"}
    members = {("A", "0", "1"): set(list(d)[:9]), ("A", "0", "2"): {"s1", "s2"}}
    cl = clusters(("A", "0", "1", 9, "unexplained", 0.02), ("A", "0", "2", 2, "explained", 0.0))
    r = L.score_absent("A", d, placed, cl, members)
    assert (r["n_placed"], r["recovered"], r["recall"], r["purity"], r["label"]) == (10, 1, 0.9, 1.0, "copy")
    cl2 = clusters(("A", "0", "1", 9, "unexplained", 0.005), ("A", "0", "2", 2, "explained", 0.0))
    assert L.score_absent("A", d, placed, cl2, members)["label"] == "allele_or_copy"
    impure = {("A", "0", "1"): set(list(d)[:9]) | {f"x{i}" for i in range(3)}}
    cl3 = clusters(("A", "0", "1", 12, "unexplained", 0.02))
    assert L.score_absent("A", d, placed, cl3, impure)["recovered"] == 0          # purity 0.75
    few = {("A", "0", "1"): set(list(d)[:7])}
    cl4 = clusters(("A", "0", "1", 7, "unexplained", 0.02))
    assert L.score_absent("A", d, placed, cl4, few)["recovered"] == 0             # recall 0.7
    cl5 = clusters(("A", "0", "1", 9, "explained", 0.02))
    assert L.score_absent("A", d, placed, cl5, members)["recovered"] == 0         # not unexplained


def test_only_the_family_s_placed_reads_count():
    d = {f"d{i}" for i in range(10)}
    placed = {f"d{i}" for i in range(5)}                                        # 5 absorbed elsewhere / unmapped
    members = {("A", "0", "1"): placed}
    r = L.score_absent("A", d, placed, clusters(("A", "0", "1", 5, "unexplained", 0.02)), members)
    assert (r["n_absent"], r["n_placed"], r["recovered"]) == (10, 5, 1)


def test_cleared_after_add_back():
    d = {f"d{i}" for i in range(10)}
    members = {("A", "3", "1"): {f"d{i}" for i in range(6)}}
    assert L.cleared("A", d, d, clusters(("A", "3", "1", 6, "explained", 0.0)), members) is True
    assert L.cleared("A", d, d, clusters(("A", "3", "1", 6, "unexplained", 0.02)), members) is False
    members2 = {("A", "3", "1"): {f"d{i}" for i in range(4)}}
    assert L.cleared("A", d, d, clusters(("A", "3", "1", 4, "unexplained", 0.02)), members2) is True   # < 50%


def test_control_flags():
    cl = clusters(("G1", "0", "1", 9, "unexplained", 0.02), ("G2", "0", "1", 9, "unexplained", 0.004),
                  ("G3", "0", "1", 9, "explained", 0.0), ("G4", "0", "1", 9, "unlinked", 0.03))
    f = L.control_flags(["G1", "G2", "G3", "G4", "G5"], cl)
    assert f == {"G1": "copy", "G2": "allele_or_copy", "G3": "-", "G4": "-", "G5": "-"}


def test_heterozygosity_bins():
    assert [L.het_bin(x) for x in (0.0, 0.0005, 0.003, 0.009, 0.02)] == ["0", "<=0.1%", "<=0.5%", "<=delta", ">delta"]


def test_pick_drops_a_family_whose_masked_copy_overlaps_its_own_sibling():
    copies = [cp("G", 0, "c7", 0, 100), cp("G", 1, "c7", 500, 1000), cp("G", 2, "c7", 600, 700)]   # G2 is last, inside G1
    assert L.pick_families(copies, catalog=copies, exclude=[], cap=10, seed=1) == []


def test_pick_admits_two_copy_families_when_asked():
    copies = [cp("B", 0, "c1", 5000, 5100), cp("B", 1, "c1", 6000, 6100), cp("Q", 0, "c9", 0, 100)]
    picks = L.pick_families(copies, catalog=copies, exclude=[], min_copies=2, cap=10, seed=1)
    assert [(p["family"], p["masked"]["idx"], len(p["survivors"])) for p in picks] == [("B", 1, 1)]


def test_pick_with_min_expressed_masks_the_last_expressed_copy_and_keeps_silent_survivors():
    copies = [cp("A", 0, "c1", 100, 200), cp("A", 1, "c1", 1000, 1100), cp("A", 2, "c2", 50, 150, n_reads=3),  # A2 silent
              cp("B", 0, "c3", 0, 100), cp("B", 1, "c3", 500, 600, n_reads=3), cp("B", 2, "c3", 900, 1000, n_reads=3)]
    picks = L.pick_families(copies, catalog=copies, exclude=[], min_expressed=2, cap=10, seed=1)
    # A: two expressed copies, the last expressed one (A1) is masked and the silent A2 stays a survivor; B: one expressed
    assert [(p["family"], p["masked"]["idx"], sorted(s["idx"] for s in p["survivors"])) for p in picks] == [("A", 1, [0, 2])]
    # the default (every copy expressed) still drops A
    assert L.pick_families(copies, catalog=copies, exclude=[], cap=10, seed=1) == []


def test_pick_skips_excluded_families():
    copies = [cp(f, k, f"c{f}", 1000 * k, 1000 * k + 100) for f in ("A", "B") for k in range(3)]
    picks = L.pick_families(copies, catalog=copies, exclude=[], exclude_families={"A"}, cap=10, seed=1)
    assert [p["family"] for p in picks] == ["B"]
