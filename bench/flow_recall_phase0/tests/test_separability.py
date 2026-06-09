import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import separability as sep

def test_best_threshold_separates():
    rows = [{"label": "real", "longcov": v} for v in (10, 12, 8, 15)] + \
           [{"label": "spurious", "longcov": v} for v in (1, 2, 1, 3)]
    auc = sep.auc_for_feature(rows, "longcov")
    assert auc > 0.9

def test_no_separation_gives_half_auc():
    rows = [{"label": "real", "longcov": v} for v in (5, 5, 5, 5)] + \
           [{"label": "spurious", "longcov": v} for v in (5, 5, 5, 5)]
    auc = sep.auc_for_feature(rows, "longcov")
    assert abs(auc - 0.5) < 1e-9
