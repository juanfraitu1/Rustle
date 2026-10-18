import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_pool16", os.path.join(HERE, "pool16.py"))
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)


def test_round_16_takes_every_unused_pool_family_with_three_copies():
    pool = [dict(family=n, copies=[{}] * k) for n, k in (("B", 3), ("A", 5), ("C", 2), ("D", 57), ("E", 4))]
    assert [f["family"] for f in P.remaining(pool, {"D"})] == ["A", "B", "E"]       # by name; D used, C too small
