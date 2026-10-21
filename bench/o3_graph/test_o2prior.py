import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2prior", os.path.join(HERE, "o2prior.py"))
PR = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PR)


def test_rescue_gives_a_tied_read_to_the_copy_with_the_most_uniquely_placed_reads():
    assert PR.rescue({"a", "b"}, {"a": 12, "b": 3}) == "a"
    assert PR.rescue({"a", "b"}, {"a": 3, "b": 3}) is None                    # equal support: abstain
    assert PR.rescue({"a", "b"}, {}) is None
    assert PR.rescue({"a", "b", "c"}, {"a": 1, "c": 5}) == "c"


def test_parsimony_assigns_only_when_exactly_one_tied_copy_is_open():
    assert PR.parsimony({"a", "b"}, {"a": 1}) == "a"                         # b has no uniquely placed read: closed
    assert PR.parsimony({"a", "b"}, {"a": 40, "b": 1}) is None                # both open: abstain, however unequal
    assert PR.parsimony({"a", "b"}, {}) is None                               # none open
    assert PR.parsimony({"a", "b", "c"}, {"b": 2, "z": 9}) == "b"


def test_world_keyed_truth_is_unkeyed_only_for_its_own_world_and_pacbio_names_survive():
    t = {"A_s1_CROSS/s001": ("determined", "a"), "B_s1_CROSS/s001": ("determined", "b")}
    assert PR.own_truth(t, "A_s1_CROSS") == {"s001": ("determined", "a")}
    pb = {"m84036_1/261690178/ccs/10900_12112": ("determined", "c")}
    assert PR.own_truth(pb, "hg002") == pb


def test_an_equal_score_tie_needs_every_kept_copy_at_the_same_best_as():
    assert PR.equal_as({"a", "b"}, {"a": 900, "b": 900, "c": 950}) is True
    assert PR.equal_as({"a", "b"}, {"a": 900, "b": 899}) is False
    assert PR.equal_as({"a", "b"}, {"a": 900}) is False                       # a kept copy without an alignment


def test_a_fragment_s_own_source_is_left_out_of_the_counts():
    per_source = {"s1": {"a": 2}, "s2": {"a": 1, "b": 1}}
    total = PR.total(per_source)
    assert total == {"a": 3, "b": 1}
    assert PR.without(total, per_source.get("s2", {})) == {"a": 2, "b": 0}
