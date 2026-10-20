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
