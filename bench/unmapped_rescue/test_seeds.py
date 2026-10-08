#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_seeds.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import seeds as S  # noqa: E402


class Choose(unittest.TestCase):
    def test_even_spread_over_the_length_order(self):
        reads = {f"r{i}": 1000 - i for i in range(10)}        # r0 longest
        self.assertEqual(S.choose_seeds(reads, 5), ["r0", "r2", "r4", "r6", "r8"])

    def test_fewer_reads_than_seeds_takes_all(self):
        self.assertEqual(sorted(S.choose_seeds({"a": 5, "b": 9}, 10)), ["a", "b"])


def fake_map(groups):
    """edges between a read and every seed of its own group (what a perfect mapper returns)"""
    g = {r: k for k, rs in groups.items() for r in rs}

    def f(queries, seed_names):
        return [(q, s) for q in queries for s in seed_names if q != s and g[q] == g[s]]
    return f


class Rounds(unittest.TestCase):
    def test_every_group_of_three_is_found_and_unseeded_groups_use_later_rounds(self):
        groups = {"A": [f"a{i}" for i in range(20)], "B": [f"b{i}" for i in range(4)], "C": ["c0", "c1"]}
        lens = {r: 1000 - j for j, r in enumerate(sum(groups.values(), []))}
        comp = S.run_rounds(lens, fake_map(groups), n_seeds=3, min_size=3, max_rounds=5)
        cl = S.G.clusters(comp, 3)
        self.assertEqual(sorted(len(v) for v in cl.values()), [4, 20])           # B needed a later round; C is too small
        self.assertEqual(len({comp[r] for r in groups["A"]}), 1)

    def test_stops_when_a_round_makes_no_progress(self):
        lens = {"x": 10, "y": 9, "z": 8}
        comp = S.run_rounds(lens, lambda q, s: [], n_seeds=2, min_size=3, max_rounds=9)
        self.assertEqual(len(set(comp.values())), 3)


class Cache(unittest.TestCase):
    def test_a_finished_round_is_read_back_without_running_minimap2(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            line = "\t".join(["q", "1000", "0", "900", "+", "s", "1000", "0", "900", "890", "900", "60", "de:f:0.001"])
            open(f"{d}/round1.paf", "w").write(line + "\n")
            open(f"{d}/round1.paf.done", "w").write("ok")
            old = os.environ.get("PATH")
            os.environ["PATH"] = ""
            try:
                f = S.minimap_map_fn({"q": "A", "s": "A"}, d, 0.00958, 0.5)
                self.assertEqual(f(["q"], ["s"]), [("q", "s")])
            finally:
                os.environ["PATH"] = old

    def test_a_fresh_round_pauses_after_the_budget(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            f = S.minimap_map_fn({"q": "ACGT" * 300, "s": "ACGT" * 300}, d, 0.00958, 0.5, fresh_budget=0)
            with self.assertRaises(S.Pause):
                f(["q"], ["s"])


if __name__ == "__main__":
    unittest.main()
