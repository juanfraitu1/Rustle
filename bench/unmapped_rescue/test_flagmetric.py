#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_flagmetric.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flagmetric as F  # noqa: E402


class GateAware(unittest.TestCase):
    def test_registered_metric(self):
        self.assertAlmostEqual(F.registered(1.0, 2, 1002, 1002), 1000 / 1002)

    def test_pure_g_prefix_with_open_gate_leaves_the_denominator(self):
        self.assertAlmostEqual(F.gate_aware(1.0, 2, 1002, 1002, "GG", True), 1.0)
        self.assertAlmostEqual(F.gate_aware(0.9995, 3, 1003, 1003, "GGG", True), 0.9995)

    def test_closed_gate_is_the_registered_metric(self):
        self.assertAlmostEqual(F.gate_aware(1.0, 2, 1002, 1002, "GG", False), F.registered(1.0, 2, 1002, 1002))

    def test_prefix_not_pure_g_or_longer_than_three_is_scored_as_before(self):
        self.assertAlmostEqual(F.gate_aware(1.0, 3, 1003, 1003, "ACT", True), F.registered(1.0, 3, 1003, 1003))
        self.assertAlmostEqual(F.gate_aware(1.0, 3, 1003, 1003, "GAG", True), F.registered(1.0, 3, 1003, 1003))
        self.assertAlmostEqual(F.gate_aware(1.0, 4, 1004, 1004, "GGGG", True), F.registered(1.0, 4, 1004, 1004))

    def test_no_prefix_is_unchanged_and_the_3_prime_end_is_never_touched(self):
        self.assertAlmostEqual(F.gate_aware(1.0, 0, 1000, 1002, "", True), F.registered(1.0, 0, 1000, 1002))

    def test_cannot_pass_by_dropping_a_real_mismatch(self):
        # identity still multiplies: a consensus at 99.5% identity does not pass because its G prefix left the denominator
        self.assertLess(F.gate_aware(0.995, 2, 1002, 1002, "GG", True), 0.999)


if __name__ == "__main__":
    unittest.main()
