import collections
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2perfamily", os.path.join(HERE, "o2perfamily.py"))
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)


def test_pair_identities_keep_same_family_pairs_aligned_over_half_the_shorter_copy():
    paf = ["A.0\t1000\t0\t900\tx\tA.1\t1000\t0\t900\t855\t900\n",       # 0.95 over 90%
           "A.0\t1000\t0\t300\tx\tA.2\t1000\t0\t300\t300\t300\n",       # 30% of the shorter: left out
           "A.0\t1000\t0\t900\tx\tB.0\t1000\t0\t900\t900\t900\n"]       # another family
    assert F.pair_identities(paf) == {("A.0", "A.1"): 0.95}


def test_a_family_qualifies_only_with_every_pair_aligned_at_90_to_under_100_percent():
    assert F.qualifies(2, [0.95]) and not F.qualifies(2, [1.0]) and not F.qualifies(2, [0.89])
    assert not F.qualifies(3, [0.95, 0.97])                  # one pair of three unaligned
    assert F.qualifies(3, [0.95, 0.97, 0.99])


def test_a_family_works_with_enough_correct_placements_true_sets_and_no_worse_precision_than_minimap2():
    n = collections.Counter({"right": 19, "wrong": 1, "judged sets": 40, "set holds truth": 39, "minimap2 right": 30,
                             "minimap2 wrong": 10})
    assert F.works(n)
    assert not F.works({**n, "right": 3, "wrong": 0})                               # too few placements
    assert not F.works({**n, "wrong": 2})                                           # 2 of 21 wrong
    assert not F.works({**n, "set holds truth": 37})                                # sets miss the truth too often
    assert not F.works({**n, "minimap2 right": 40, "minimap2 wrong": 0})            # minimap2 is better


def test_fragments_are_judged_by_their_source_read_for_o2_and_for_minimap2_at_mapq_1():
    rows = [dict(read="s1|f0|0|500", family="A", call="A.0", candidates="A.0,A.1", primary="A.0"),
            dict(read="s1|f1|600|500", family="A", call="-", candidates="A.0,A.1", primary="A.1"),
            dict(read="s2|f0|0|500", family="A", call="A.1", candidates="A.0,A.1", primary="A.1"),
            dict(read="s3|f0|0|500", family="B", call="B.0", candidates="B.0,B.1", primary="B.0"),
            dict(read="s4|f0|0|500", family="A", call="A.0", candidates="A.1,A.2", primary="A.1")]   # truth not a candidate
    truth = {"s1|f0|0|500": "A.0", "s1|f1|600|500": "A.0", "s2|f0|0|500": "A.0", "s3|f0|0|500": "B.0", "s4|f0|0|500": "A.0"}
    mapq = {"s1|f0|0|500": 3, "s1|f1|600|500": 1, "s2|f0|0|500": 0, "s3|f0|0|500": 5, "s4|f0|0|500": 9}
    got = F.fair_counts(rows, truth, mapq, families={"A"})
    assert got["A"] == {"o2 right": 1, "o2 wrong": 1, "mm2 right": 0, "mm2 wrong": 1}   # s1: O2 right; s2: O2 wrong (MAPQ 0)
    assert "B" not in got                                                               # not a qualifying family
