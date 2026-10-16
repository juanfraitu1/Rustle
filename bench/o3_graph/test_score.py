import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    """By file path under a unique name: bench has a `{name}.py` in many directories (score.py in 13), and a plain import would
    shadow another directory's module in a combined pytest run."""
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load("score")

E = {f"FE|iso0|{i}" for i in range(10)}
A = {f"FA|iso0|{i}" for i in range(10)}


def test_a_call_recovers_the_absent_copy_with_recall_and_purity():
    calls = [E - {"FE|iso0|0"}, set(list(A)[:3])]
    det, false = S.family_outcome(calls, source=E, placed_source=E)
    assert det and false                       # 9 / 10 at purity 1.0; the A-only call is a false call


def test_without_an_absent_copy_every_call_is_false_and_low_recall_misses():
    assert S.family_outcome([set(list(A)[:5])], source=set(), placed_source=set()) == (None, True)
    assert S.family_outcome([], source=set(), placed_source=set()) == (None, False)
    det, false = S.family_outcome([set(list(E)[:5])], source=E, placed_source=E)
    assert det is False and false is False     # 5 / 10 of the source: no recovery, but not a false call either
