import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pairs_lib as P  # noqa: E402

LF = P.LF


def locus(source, cls, contig, s, e, short=False):
    """A Liftoff locus as `_liftoff.read_loci` returns it (start0 / end stay strings, as in the file)."""
    return {"source_id": source, "cls": cls, "contig": contig, "start0": str(s), "end": str(e), "iv": [(s, e)],
            "short": short, "seqid_f": 1.0}


def fam(family, contig, s, e):
    return {"contig": contig, "family": family, "iv": [(s, e)]}


def test_dist_bin_edges():
    a = locus("S1", "in_place", "chr1", 1_000, 2_000)
    at = lambda s: locus("S2", "in_place", "chr1", s, s + 100)  # noqa: E731
    assert P.dist_bin(a, at(1_500)) == "lt100k"            # overlapping spans: gap 0
    assert P.dist_bin(a, at(2_000 + 99_999)) == "lt100k"
    assert P.dist_bin(a, at(2_000 + 100_000)) == "100k_1M"
    assert P.dist_bin(a, at(2_000 + 1_000_000)) == "1M_10M"
    assert P.dist_bin(a, at(2_000 + 10_000_000)) == "ge10M"
    assert P.dist_bin(a, locus("S2", "in_place", "chr2", 1_000, 2_000)) == "other_contig"


def test_cover_families_matches_figure8_rule():
    a = locus("S1", "in_place", "chr1", 1_000, 2_000)
    c = locus("S1", "extra_copy", "chr1", 50_000, 51_000)
    half = locus("S2", "in_place", "chr1", 90_000, 91_000)
    fams = [fam("F1", "chr1", 1_000, 2_000), fam("F1", "chr1", 50_400, 51_000), fam("F2", "chr1", 90_000, 90_499)]
    cov = P.cover_families([a, c, half], fams)
    assert cov[LF.locus_key(a)] == {"F1"} and cov[LF.locus_key(c)] == {"F1"}  # 600 of 1,000 bases: covered
    assert cov[LF.locus_key(half)] == frozenset()                            # 499 of 1,000 bases: not covered
    (src, both, shared), = LF.pair_families([(a, c)], fams)
    assert both and shared == cov[LF.locus_key(a)] & cov[LF.locus_key(c)]


def toy_pool():
    a = locus("S1", "in_place", "chr1", 1_000, 2_000)
    c = locus("S1", "extra_copy", "chr1", 50_000, 51_000)        # gap 48,000: lt100k
    x1 = locus("S2", "in_place", "chr1", 60_000, 61_000)         # lt100k, covered: a level-0 decoy
    x2 = locus("S3", "in_place", "chr1", 500_000, 501_000)       # 100k_1M, covered
    x3 = locus("S2", "in_place", "chr1", 1_500, 1_800)           # overlaps a: never a decoy
    x4 = locus("S4", "in_place", "chr1", 70_000, 71_000)         # not covered
    x5 = locus("S1", "extra_copy", "chr1", 80_000, 81_000)       # a's own source: never a decoy
    x6 = locus("S5", "in_place", "chr2", 1_000, 2_000)           # other contig, covered
    fams = {LF.locus_key(r): frozenset(f) for r, f in ((a, {"F1"}), (c, {"F1"}), (x1, {"F2"}), (x2, {"F3"}),
                                                       (x3, {"F1"}), (x4, set()), (x5, {"F1"}), (x6, {"F4"}))}
    return a, c, [x1, x2, x3, x4, x5, x6], fams


def names(pool, idx):
    return [pool.rows[i]["start0"] for i in idx]


def test_candidates_match_and_widen():
    a, c, xs, fams = toy_pool()
    pool = P.Pool([a, c] + xs, fams)
    idx, level = pool.candidates(a, c, True)
    assert level == 0 and names(pool, idx) == ["60000"]
    pool = P.Pool([a, c] + [x for x in xs if x["start0"] != "60000"], fams)
    idx, level = pool.candidates(a, c, True)
    assert level == 1 and names(pool, idx) == ["500000"]
    pool = P.Pool([a, c] + [x for x in xs if x["start0"] not in ("60000", "500000")], fams)
    idx, level = pool.candidates(a, c, True)
    assert level == 2 and names(pool, idx) == ["1000"] and pool.rows[idx[0]]["contig"] == "chr2"


def test_candidates_for_another_contig_and_uncovered_copy():
    a, _, xs, fams = toy_pool()
    far = locus("S1", "extra_copy", "chr2", 9_000, 9_500)
    fams[LF.locus_key(far)] = frozenset()
    pool = P.Pool([a, far] + xs, fams)
    idx, level = pool.candidates(a, far, False)   # uncovered copy on another contig: none matches, so level 2
    assert level == 2 and names(pool, idx) == ["70000"]
    idx, level = pool.candidates(a, far, True)    # covered copy on another contig: chr2's covered locus
    assert level == 0 and names(pool, idx) == ["1000"] and pool.rows[idx[0]]["contig"] == "chr2"


def test_no_candidate_leaves_the_test():
    a = locus("S1", "in_place", "chr1", 1_000, 2_000)
    c = locus("S1", "extra_copy", "chr1", 50_000, 51_000)
    fams = {LF.locus_key(a): frozenset({"F1"}), LF.locus_key(c): frozenset({"F1"})}
    r = P.run_test([(a, c)], P.Pool([a, c], fams), fams, 10, 7)
    assert r["levels"] == {"0": 0, "1": 0, "2": 0, "none": 1} and r["n"] == 0 and r["p"] != r["p"]


def test_run_test_decoys_never_shared():
    a, c, xs, fams = toy_pool()
    r = P.run_test([(a, c)], P.Pool([a, c] + xs, fams), fams, 50, 7)
    assert (r["n"], r["k"], r["null_max"]) == (1, 1, 0) and r["p"] == 1 / 51 and r["levels"]["0"] == 1


def test_run_test_decoys_always_shared():
    a, c, xs, fams = toy_pool()
    fams[LF.locus_key(xs[0])] = frozenset({"F1", "F9"})   # the only level-0 decoy shares a's family
    r = P.run_test([(a, c)], P.Pool([a, c] + xs, fams), fams, 50, 7)
    assert r["k"] == 1 and r["null_mean"] == 1.0 and r["p"] == 1.0


def test_run_test_is_seeded():
    a, c, xs, fams = toy_pool()
    pool = P.Pool([a, c] + xs, fams)
    far = locus("S1", "extra_copy", "chr2", 9_000, 9_500)
    fams[LF.locus_key(far)] = frozenset({"F1"})
    pairs = [(a, c), (a, far)]
    assert P.run_test(pairs, pool, fams, 200, 7) == P.run_test(pairs, pool, fams, 200, 7)
