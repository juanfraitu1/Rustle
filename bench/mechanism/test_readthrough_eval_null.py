import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load():
    spec = importlib.util.spec_from_file_location("readthrough_eval_null", os.path.join(HERE, "readthrough_eval.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


RE = _load()


def test_null_candidates_leave_out_the_flagged_the_protected_and_the_rare():
    # review I3: a junction the fusion protection kept is no NULL candidate either (the unprotected list: arm would remove
    # the fusion reads arm R kept)
    S = {(100, 200, "+"): 5, (300, 400, "+"): 7, (500, 600, "+"): 9, (700, 800, "+"): 1, (900, 950, "-"): 3}
    flagged = {(100, 200, "+")}
    protected = {(300, 400, "+")}
    cand = RE.null_candidates(S, flagged, protected, lambda s, e, st: True)
    assert cand == {(500, 600, "+"): 9, (900, 950, "-"): 3}
    assert RE.null_candidates(S, flagged, set(), lambda s, e, st: st == "+") == {(300, 400, "+"): 7, (500, 600, "+"): 9}
