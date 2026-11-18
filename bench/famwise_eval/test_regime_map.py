"""Tests of regime_map.py (bench/famwise_eval): the cell specs and the several-fold rule."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import regime_map as M

def test_cells():
    assert len(M.CELLS) == 38 and "k0_1to1_full" in M.CELLS and "k20_silent_frag" in M.CELLS and "array_frag" in M.CELLS

def test_spec_for():
    s = M.spec_for("k0_1to1_full", 3, "/t/template.json")
    assert s["seed"] == 3 and s["template"] == {"source": "file", "path": "/t/template.json"}
    a, b = s["copies"]
    assert (a["id"], a["expression"], a["ops"]) == ("A", 200, []) and (b["expression"], b["ops"]) == (200, [])     # k = 0: no op
    assert a["isoforms"] == b["isoforms"] == [{"skip": [3], "weight": 0.3}]
    assert "trunc5_frac" not in s["reads"]
    s = M.spec_for("k20_10to1_frag", 1, "/t/template.json")
    assert s["copies"][1]["ops"] == [{"op": "snp", "rate": 0.01, "region": "exons"}] and s["copies"][1]["expression"] == 20
    assert (s["reads"]["trunc5_frac"], s["reads"]["trunc5_max"]) == (1.0, 0.7)
    assert M.spec_for("k4_silent_full", 1, "/t/template.json")["copies"][1]["expression"] == 0
    s = M.spec_for("array_full", 1, "/t/template.json")
    assert len(s["copies"]) == 8 and all(c["ops"] == [] and c["expression"] == 25 for c in s["copies"])
    assert M.spec_for("k0_1to1_full", 1, None)["template"] == {"source": "synthetic", "exons": [300, 200, 400, 250, 150, 700], "introns": [300] * 5}

def test_several_fold():
    assert M.several_fold([(0, 2), (0, 2), (0, 1)]) == (0, 2, "yes")        # P recovers none, O some
    assert M.several_fold([(2, 3), (2, 3), (2, 4)]) == (2, 3, "no")         # medians 3 / 2
    assert M.several_fold([(1, 2), (2, 4), (1, 3)]) == (1, 3, "yes")
    assert M.several_fold([(0, 0), (0, 0)]) == (0, 0, "no")
