"""PREREG_o1_self_hits_2026-10-21: self-only copies and the bars."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o1_self_hits_eval", os.path.join(HERE, "self_hits_eval.py"))
E = importlib.util.module_from_spec(spec)
spec.loader.exec_module(E)


def test_a_record_on_the_main_diagonal_is_a_self_hit_and_a_tandem_shift_is_not():
    assert E.is_self_hit(("c", 1001, 9000), 2000, 8000, "+", ("c", 3001, 12000), 0, 6000)
    assert not E.is_self_hit(("c", 1001, 9000), 0, 3000, "+", ("c", 3001, 12000), 3000, 6000)
    assert not E.is_self_hit(("c", 1001, 9000), 2000, 8000, "-", ("c", 3001, 12000), 0, 6000)
    assert not E.is_self_hit(("c", 1001, 9000), 2000, 8000, "+", ("d", 3001, 12000), 0, 6000)


def test_a_copy_is_self_only_when_its_best_same_family_evidence_at_090_is_a_self_hit():
    best = {"a": {"self": 1.0, "non": 0.85}, "b": {"self": 1.0, "non": 0.95}, "c": {"non": 0.99}, "d": {"self": 0.8}}
    assert E.self_only(best) == {"a"}


def test_b1_counts_self_only_copies_fixed_when_out_of_a_family_or_with_a_non_self_partner():
    fixed_best = {"a": {"non": 0.92}, "b": {"self": 1.0}}
    in_family = {"a", "b"}
    n_fixed, n = E.b1({"a", "b", "c"}, fixed_best, in_family)
    assert (n_fixed, n) == (2, 3)          # a: non-self partner; b: still self-only -> not fixed; c: left the families -> fixed


def test_b2_counts_control_sd_paired_copies_still_paired():
    assert E.b2({"a", "b", "c", "d"}, {"a", "b", "c"}) == (3, 4)


def _cp(key, family, chrom, strand, exons):
    return dict(key=key, family=family, chrom=chrom, strand=strand, start=exons[0][0], end=exons[-1][1], exons=exons)


def test_b2_by_locus_keeps_a_copy_whose_locus_is_sd_paired_through_a_representative_switch():
    ctl = {"a": _cp("a", "F1", "c", "+", [(100, 200)]), "b": _cp("b", "F1", "c", "+", [(5000, 5100)])}
    fix = [_cp("a2", "G1", "c", "-", [(150, 250)]), _cp("z", "G1", "d", "+", [(1, 50)])]
    assert E.b2_locus({"a", "b"}, ctl, {"a2"}, fix) == (1, 2)   # a: its exon mate a2 is SD-paired; b: no fix copy at its locus


def test_cross_family_exon_pairs_are_split_by_strand_and_ignore_same_family_pairs():
    cps = [_cp("x", "F1", "c", "+", [(100, 200)]), _cp("y", "F2", "c", "+", [(150, 160)]), _cp("z", "F3", "c", "-", [(190, 300)]),
           _cp("w", "F1", "c", "+", [(120, 130)]), _cp("v", "F4", "c", "+", [(200, 210)])]
    assert E.cross_family_exon_pairs(cps) == (1, 2)             # same strand: x-y; opposite: x-z, z-v; x-w same family; x-v abut


def test_g1_tolerates_the_larger_of_five_pairs_and_five_percent():
    assert E.g1(100, 105) and not E.g1(100, 106)
    assert E.g1(0, 5) and not E.g1(0, 6)
    assert E.g1(1368, 1437) and not E.g1(1368, 1438)


def test_span_disjoint_loci():
    assert E.span_disjoint(("c", 1, 100), ("c", 101, 200)) and E.span_disjoint(("c", 1, 100), ("d", 1, 100))
    assert not E.span_disjoint(("c", 1, 100), ("c", 100, 200))


def test_g2_keeps_a_distant_pair_still_sharing_a_family_through_exon_mates():
    ctl = {"a": _cp("a", "F1", "c", "+", [(100, 200)]), "b": _cp("b", "F1", "d", "+", [(100, 200)]),
           "p": _cp("p", "F2", "c", "+", [(9000, 9100)]), "q": _cp("q", "F2", "e", "+", [(1, 90)])}
    fix = [_cp("a2", "G1", "c", "+", [(150, 160)]), _cp("b", "G1", "d", "+", [(100, 200)]),
           _cp("p", "G2", "c", "+", [(9000, 9100)]), _cp("q", "G3", "e", "+", [(1, 90)])]
    assert E.g2({("a", "b"), ("p", "q")}, ctl, fix) == (2, 4)   # a-b kept (a2 ~ a); p-q split


def test_g3_allows_the_larger_of_two_copies_and_ten_percent_of_the_control():
    assert E.g3(2, 0) and not E.g3(3, 0)
    assert E.g3(5, 50) and not E.g3(6, 50)


def test_verdict_fails_on_any_guard_and_is_underpowered_below_ten_self_only_copies():
    assert E.verdict((9, 9), True, True, True, True) == "UNDERPOWERED"
    assert E.verdict((9, 9), True, False, True, True) == "FAIL"
    assert E.verdict((90, 100), True, True, True, True) == "PASS"
    assert E.verdict((89, 100), True, True, True, True) == "FAIL"


def test_best_by_family_collects_distant_non_self_pairs_at_090(tmp_path):
    rows = [["c:1-1000", 1000, 0, 1000, "+", "d:1-1000", 1000, 0, 1000, 950, 1000, 60],      # distant, non-self 0.95: a pair
            ["c:1-1000", 1000, 100, 900, "+", "c:101-2000", 1900, 0, 800, 800, 800, 60],     # diagonal self-hit: not a pair
            ["c:1-1000", 1000, 0, 1000, "+", "e:1-1000", 1000, 0, 1000, 850, 1000, 60]]      # distant but 0.85: not a pair
    paf = tmp_path / "x.paf"
    paf.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")
    key2gene = {"c:1-1000": "a", "d:1-1000": "b", "c:101-2000": "c", "e:1-1000": "e"}
    fam = {"a": "F", "b": "F", "c": "F", "e": "F"}
    best, pairs = E.best_by_family(str(paf), key2gene, {"control": fam, "fix": fam})
    assert pairs == {("a", "b")}
    assert best["control"]["a"] == {"non": 0.95, "self": 1.0}
