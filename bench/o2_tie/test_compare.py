import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import compare as C


def r(status, copy, margin):
    return {"status": status, "assigned_copy": copy, "margin": margin, "p_value": "1e-9"}


def test_a_zero_margin_assigned_row_turning_ambiguous_is_the_tie_fix():
    assert C.classify(r("assigned", "1", "0.000"), r("ambiguous", "1", "0.000")) == "status_tie_fix"


def test_a_real_margin_status_change_is_not_allowed():
    assert C.classify(r("assigned", "1", "6.905"), r("ambiguous", "1", "6.905")) == "status_OTHER"


def test_a_copy_change_without_status_change_is_reported_separately():
    assert C.classify(r("tied", "0", "0.000"), r("tied", "1", "0.000")) == "copy_only"


def test_identical_rows():
    assert C.classify(r("tied", "0", "0.000"), r("tied", "0", "0.000")) is None
