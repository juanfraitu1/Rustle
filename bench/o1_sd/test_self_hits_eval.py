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
