"""Tests for o2_population.py (docs/PREREG_o2_population_2026-10-12.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import o2_population as op  # noqa: E402


def test_molecule_class_priority():
    # each pair = (verdict, columns)
    assert op.molecule_class([("A", 5), ("B", 4)]) == "moved"
    assert op.molecule_class([("A", 5), ("A", 9)]) == "confirmed"
    assert op.molecule_class([(None, 0), (None, 0)]) == "identical"
    assert op.molecule_class([("A", 5), (None, 0)]) == "undecided"
    assert op.molecule_class([(None, 3)]) == "undecided"


def test_wilson_interval_contains_the_estimate():
    lo, hi = op.wilson(30, 1000)
    assert lo < 0.03 < hi and 0.0 <= lo and hi <= 1.0
    assert op.wilson(0, 100)[0] == 0.0
