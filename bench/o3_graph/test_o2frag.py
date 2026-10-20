import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2frag", os.path.join(HERE, "o2frag.py"))
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)


def test_fragment_windows_are_seeded_inside_the_read_and_within_the_length_bounds():
    w = F.windows("r1", 2000, 3, 300, 1200, seed=7)
    assert len(w) == 3 and w == F.windows("r1", 2000, 3, 300, 1200, seed=7) != F.windows("r1", 2000, 3, 300, 1200, seed=8)
    assert all(300 <= n <= 1200 and 0 <= s and s + n <= 2000 for s, n in w)
    assert F.windows("r2", 250, 3, 300, 1200, seed=7) == []                        # shorter than a fragment
    assert all(s + n <= 500 for s, n in F.windows("r3", 500, 3, 300, 1200, seed=1))  # capped at the read


def test_a_source_read_s_copy_is_its_truth_only_when_its_full_alignment_is_clear_by_the_margin_and_on_its_copy():
    T = F.source_truth
    al = [dict(AS=2000, copy="F.0", prim=True), dict(AS=1950, copy="F.1", prim=False), dict(AS=1990, copy=None, prim=False)]
    assert T(al[:2], "F.0", margin=30) == "F.0"                                    # 50 AS clear of its paralog
    assert T(al, "F.0", margin=30) is None                                         # an off-copy locus 10 AS away
    assert T(al[:2], "F.1", margin=30) is None                                     # realigned onto another copy
    assert T(al[:1], "F.0", margin=30) == "F.0"                                    # no rival at all
    assert T([dict(AS=2000, copy="F.0", prim=True), dict(AS=2000, copy="F.0", prim=False)], "F.0", margin=30) == "F.0"


def test_a_fragment_s_name_carries_its_source_read():
    assert F.source_of("m64/1/ccs|f2|130|700") == "m64/1/ccs" and F.source_of("m64/1/ccs") is None


def test_the_sample_keeps_the_reads_of_lowest_seeded_rank_whatever_the_order_they_come_in():
    ranks = {f"r{i}": F._rank(1, f"r{i}") for i in range(50)}
    want = sorted(ranks, key=ranks.get)[:5]
    for order in (sorted(ranks), sorted(ranks, reverse=True)):
        h = []
        for r in order:
            F.keep_lowest(h, 5, ranks[r], r)
        assert sorted(x for _, x in h) == sorted(want)
