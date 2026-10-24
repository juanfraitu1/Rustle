"""Tests for seed_reach.py (docs/archive/2026-10/merged/PREREG_lrpap1_seed_reach_2026-10-11.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import seed_reach as R  # noqa: E402

COPIES = {"a": ("c1", 100, 200), "b": ("c1", 500, 600), "c": ("c2", 100, 200), "d": ("c3", 0, 50)}


def test_finds_by_overlap_on_the_same_chromosome():
    hits = [("a", "c1", 150, 160), ("a", "c2", 190, 300), ("a", "c3", 60, 90), ("b", "c1", 590, 700), ("b", "c9", 0, 10)]
    f, extra = R.finds(hits, COPIES)
    assert f == {"a": {"a", "c"}, "b": {"b"}, "c": set(), "d": set()}
    assert extra == {"a": 1, "b": 1, "c": 0, "d": 0}


def test_one_pass_and_closure():
    F = {"a": {"b"}, "b": {"c"}, "c": set(), "d": set()}
    assert R.one_pass({"a"}, F) == {"a", "b"}
    assert R.closure({"a"}, F) == {"a", "b", "c"}
    assert R.closure({"a", "d"}, F) == {"a", "b", "c", "d"}


def test_min_seeds_and_curve():
    F = {"a": {"b"}, "b": {"c"}, "c": set(), "d": set()}
    ids = sorted(F)
    assert R.min_seeds(F, ids, R.closure) == 2            # a reaches b, c; d needs its own seed
    assert R.min_seeds(F, ids, R.one_pass) == 3           # {a, b, d} or {a, c, d}
    cur = R.curve(F, ids, R.closure)
    assert cur[0] == dict(k=1, n_sets=4, min=1, median=1.5, max=3, all_found=0)
    assert cur[-1] == dict(k=4, n_sets=1, min=4, median=4, max=4, all_found=1)


def test_source_components_count_seeds_needed_with_iteration():
    F = {"a": {"b"}, "b": {"a"}, "c": {"a"}, "d": set()}
    assert R.source_components(F) == [["c"], ["d"]]
    assert R.min_seeds(F, sorted(F), R.closure) == len(R.source_components(F))


def test_threshold_sweep_raises_seeds_needed_as_identity_rises():
    copies = {"a": ("c", 0, 10), "b": ("c", 20, 30), "c": ("c", 40, 50)}
    hits = [("a", "c", 0, 10, 1.0), ("b", "c", 20, 30, 1.0), ("c", "c", 40, 50, 1.0),
            ("a", "c", 20, 30, 0.99), ("b", "c", 0, 10, 0.99), ("a", "c", 40, 50, 0.95), ("c", "c", 0, 10, 0.95)]
    rows = R.sweep(hits, copies, ["a", "b", "c"], [0.90, 0.97, 0.995])
    assert [(r["threshold"], r["min_seeds_one_pass"], r["min_seeds_iterate"], r["components"]) for r in rows] == \
        [(0.90, 1, 1, [["a", "b", "c"]]), (0.97, 2, 2, [["a", "b"], ["c"]]), (0.995, 3, 3, [["a"], ["b"], ["c"]])]
