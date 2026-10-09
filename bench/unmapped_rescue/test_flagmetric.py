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


class O3Class(unittest.TestCase):
    def test_copy_allele_present_and_none(self):
        self.assertEqual(F.o3_class(0.980, 1.0), "COPY")          # more than the allele cutoff from every reference locus
        self.assertEqual(F.o3_class(0.995, 1.0), "ALLELE")        # within the allele cutoff of a reference locus but not identical
        self.assertEqual(F.o3_class(0.9995, 1.0), "PRESENT")
        self.assertEqual(F.o3_class(None, 1.0), "COPY")           # nothing on the reference at all
        self.assertIsNone(F.o3_class(0.5, 0.9))                   # not reconstructed on the truth genome: no call

    def test_boundaries(self):
        self.assertEqual(F.o3_class(1 - 0.00958, 0.999), "ALLELE")        # exactly the cutoff is still an allele
        self.assertEqual(F.o3_class(1 - 0.00958 - 1e-6, 0.999), "COPY")
        self.assertEqual(F.o3_class(0.999, 0.999), "PRESENT")

    def test_the_registered_flag_is_copy_or_allele(self):
        self.assertTrue(F.registered_flag(0.995, 1.0))
        self.assertTrue(F.registered_flag(0.980, 1.0))
        self.assertFalse(F.registered_flag(0.9995, 1.0))
        self.assertFalse(F.registered_flag(0.5, 0.9))


class Discovery(unittest.TestCase):
    def test_classes(self):
        self.assertEqual(F.discovery_class(0.95, 0.5, 1.0), "CONFIRMED")          # absent from the primary, present in the paternal assembly
        self.assertEqual(F.discovery_class(0.95, 1.0, 0.2), "CONFIRMED")
        self.assertEqual(F.discovery_class(0.95, 0.9, 0.9), "DIVERGED")          # a hit exists, nowhere at the bar
        self.assertEqual(F.discovery_class(0.85, 0.89, 0.5), "NOVEL")            # no assembly holds 90% of it
        self.assertEqual(F.discovery_class(None, None, None), "NOVEL")
        self.assertEqual(F.discovery_class(0.995, 1.0, 1.0), "ALLELE-LIKE")
        self.assertEqual(F.discovery_class(1.0, 1.0, 1.0), "PRESENT")

    def test_boundaries(self):
        self.assertEqual(F.discovery_class(1 - 0.00958, 0.0, 0.0), "ALLELE-LIKE")
        self.assertEqual(F.discovery_class(1 - 0.00958 - 1e-6, 0.999, 0.0), "CONFIRMED")
        self.assertEqual(F.discovery_class(0.95, 0.9989, 0.0), "DIVERGED")
        self.assertEqual(F.discovery_class(0.95, 0.0, 0.0), "DIVERGED")          # the primary itself holds 95%: not absent
        self.assertEqual(F.discovery_class(0.5, 0.899, 0.0), "NOVEL")
        self.assertEqual(F.discovery_class(0.999, 0.0, 0.0), "PRESENT")


class DiscoveryOtherAnimals(unittest.TestCase):
    """Amendment 33: R is the animal's own reference; E is the best score over every other assembly"""

    def test_classes(self):
        self.assertEqual(F.elsewhere_class(1.0, []), "PRESENT")
        self.assertEqual(F.elsewhere_class(0.995, [0.2]), "ALLELE-LIKE")
        self.assertEqual(F.elsewhere_class(0.4, [0.99, 0.5]), "ELSEWHERE")        # absent from its own reference, an assembly holds it
        self.assertEqual(F.elsewhere_class(None, [None, 0.93]), "ELSEWHERE")      # E is the maximum, a missing score is 0
        self.assertEqual(F.elsewhere_class(0.85, [0.89, 0.5]), "NOVEL")           # nothing holds 90%
        self.assertEqual(F.elsewhere_class(None, []), "NOVEL")
        self.assertEqual(F.elsewhere_class(None, [None, None]), "NOVEL")
        self.assertEqual(F.elsewhere_class(0.95, [0.5]), "DIVERGED")              # its own reference holds a distant hit, no other assembly does

    def test_boundaries(self):
        self.assertEqual(F.elsewhere_class(1 - 0.00958, [0.0]), "ALLELE-LIKE")
        self.assertEqual(F.elsewhere_class(0.5, [0.9]), "ELSEWHERE")
        self.assertEqual(F.elsewhere_class(0.5, [0.8999]), "NOVEL")
        self.assertEqual(F.elsewhere_class(0.999, [0.0]), "PRESENT")

    def test_amendment_35_other_assembly_must_beat_the_reference(self):
        self.assertEqual(F.elsewhere_class(0.988, [0.965]), "DIVERGED")           # the own reference explains it better
        self.assertEqual(F.elsewhere_class(0.974, [0.974]), "DIVERGED")           # a tie is not "elsewhere"
        self.assertEqual(F.elsewhere_class(0.95, [0.9]), "DIVERGED")
        self.assertEqual(F.elsewhere_class(1 - 0.00958 - 1e-6, [0.9]), "DIVERGED")
        self.assertEqual(F.elsewhere_class(0.95, [0.97]), "ELSEWHERE")            # another assembly holds it better
        self.assertEqual(F.elsewhere_class(1 - 0.00958 - 1e-6, [0.999]), "ELSEWHERE")     # the other assembly holds it completely
        self.assertEqual(F.elsewhere_class(0.9, [0.95]), "ELSEWHERE")
        self.assertEqual(F.elsewhere_class(0.91, [0.8999]), "DIVERGED")
        self.assertEqual(F.elsewhere_class(None, [0.9]), "ELSEWHERE")

    def test_amendment_35b_tolerance(self):
        self.assertEqual(F.elsewhere_class(0.982, [0.983]), "DIVERGED")           # the gorilla-control near tie
        self.assertEqual(F.elsewhere_class(0.95, [0.9595]), "DIVERGED")           # under the allele tolerance above R
        self.assertEqual(F.elsewhere_class(0.95, [0.9596]), "ELSEWHERE")
        self.assertEqual(F.elsewhere_class(0.9, [0.9095]), "DIVERGED")            # 0.9 + 0.00958 = 0.90958
        self.assertEqual(F.elsewhere_class(0.9, [0.9097]), "ELSEWHERE")
        self.assertEqual(F.elsewhere_class(0.99, [0.9989]), "DIVERGED")           # under 0.999 and within the tolerance of R
        self.assertEqual(F.elsewhere_class(0.99, [0.999]), "ELSEWHERE")
        self.assertEqual(F.elsewhere_class(None, [0.89]), "NOVEL")

    def test_amendment_35c_cross_species_assemblies_need_the_reference_to_lack_it(self):
        # the chimp seed-6 controls: the own reference holds 0.97, a human / gorilla homolog scores higher
        self.assertEqual(F.elsewhere_class(0.971, [], cross=[0.968, 0.981]), "DIVERGED")
        self.assertEqual(F.elsewhere_class(0.970, [], cross=[0.994, 0.994]), "DIVERGED")
        self.assertEqual(F.elsewhere_class(0.062, [], cross=[0.989]), "ELSEWHERE")       # the own reference lacks it
        self.assertEqual(F.elsewhere_class(None, [], cross=[0.9]), "ELSEWHERE")
        self.assertEqual(F.elsewhere_class(0.8999, [], cross=[0.95]), "ELSEWHERE")
        self.assertEqual(F.elsewhere_class(0.9, [], cross=[0.999]), "DIVERGED")
        self.assertEqual(F.elsewhere_class(None, [], cross=[0.8999]), "NOVEL")
        self.assertEqual(F.elsewhere_class(0.5, [], cross=[0.5]), "NOVEL")

    def test_amendment_35c_same_species_keeps_the_tolerance_rule(self):
        self.assertEqual(F.elsewhere_class(0.97, [0.999], cross=[0.5]), "ELSEWHERE")     # copy absent from the primary, divergent paralog on it, whole sequence on the other haplotype
        self.assertEqual(F.elsewhere_class(0.982, [0.983], cross=[0.99]), "DIVERGED")
        self.assertEqual(F.elsewhere_class(None, [0.95], cross=[]), "ELSEWHERE")
        self.assertEqual(F.elsewhere_class(0.97, [0.8], cross=[0.999]), "DIVERGED")       # a cross-species homolog does not count when R holds it


if __name__ == "__main__":
    unittest.main()
