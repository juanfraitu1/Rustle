import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import v1_compare as V  # noqa: E402


def sedef_row(c1, s1, e1, c2, s2, e2, ident):
    f = [c1, str(s1), str(e1), c2, str(s2), str(e2)] + ["x"] * 14 + [str(ident)] + ["y"] * 13
    assert len(f) == 34
    return "\t".join(f) + "\n"


def rs_row(c1, s1, e1, c2, s2, e2, est, status):
    return "\t".join([c1, str(s1), str(e1), c2, str(s2), str(e2), "+", "100", str(est), str(est), status]) + "\n"


def test_sedef_rows_are_kept_by_identity_and_the_shorter_side():
    lines = ["#header\n", sedef_row("chr1", 0, 2000, "chr2", 0, 2000, 0.95), sedef_row("chr1", 0, 2000, "chr2", 0, 2000, 0.89),
             sedef_row("chr1", 0, 2000, "chr2", 0, 999, 0.99)]
    assert V.read_pairs(lines, "sedef", 0.90, 1000) == [(("chr1", 0, 2000), ("chr2", 0, 2000), 0.95)]


def test_sedefrs_rows_keep_ok_only_with_percent_identity():
    lines = ["#h\n", rs_row("chr1", 0, 2000, "chr2", 0, 2000, 97.6, "ok"), rs_row("chr1", 0, 2000, "chr1", 0, 2000, 100.0, "self"),
             rs_row("chr3", 0, 2000, "chr4", 0, 2000, 95.0, "weak"), rs_row("chr5", 0, 2000, "chr6", 0, 2000, 89.9, "ok")]
    p = V.read_pairs(lines, "sedefrs", 0.90, 1000)
    assert [(a, b) for a, b, _ in p] == [(("chr1", 0, 2000), ("chr2", 0, 2000))]
    assert abs(p[0][2] - 0.976) < 1e-12


def test_recall_and_precision_count_only_universe_bases():
    ref = [(("chr1", 0, 1000), ("chr1", 5000, 7000), 0.95)]
    test = [(("chr1", 500, 6000), ("chr2", 0, 1000), 0.95)]
    res, _ = V.compare(ref, test, {"chr1": [(0, 10_000)]})   # chr2 is outside the universe
    assert (res["ref_bp"], res["ref_covered"]) == (3000, 1500)
    assert (res["test_bp"], res["test_covered"]) == (5500, 1500)
    assert res["recall"] == 0.5
    assert abs(res["precision"] - 1500 / 5500) < 1e-12


def test_universe_from_chain_is_the_v2_side_of_the_aligned_blocks():
    chain = ["chain 100 chr1 1000 + 0 1000 chr1 1010 + 0 1010 1\n", "300 0 10\n", "700\n"]
    # v1 [0,300) -> v2 [0,300); 10 bases inserted on v2; v1 [300,1000) -> v2 [310,1010)
    assert V.universe_from_chain(chain) == {"chr1": [(0, 300), (310, 1010)]}


def test_recall_bins_split_by_identity_and_length():
    ref = [(("chr1", 0, 2000), ("chr1", 10_000, 12_000), 0.96), (("chr1", 20_000, 26_000), ("chr1", 30_000, 36_000), 0.995)]
    b = V.recall_bins(ref, {"chr1": [(0, 2000), (20_000, 23_000)]}, {"chr1": [(0, 100_000)]})
    assert b["identity"]["0.95-0.98"] == {"bp": 4000, "recall": 0.5}
    assert b["identity"]["0.99-1.0000001"] == {"bp": 12000, "recall": 0.25}
    assert b["length"]["1000-5000"]["recall"] == 0.5
    assert b["length"]["5000-20000"]["recall"] == 0.25
