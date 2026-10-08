#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_trim5g.py"""
import collections
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import synth_world as W  # noqa: E402
import trim5g as T  # noqa: E402


class Prefix(unittest.TestCase):
    def test_prefix_of_the_best_covering_copy_of_the_family(self):
        hsps = [("c1", "f01:A", 3, 400), ("c1", "f01:A", 410, 900),      # copy A covers 3-400 and 410-900
                ("c1", "f01:B", 6, 200),                                    # copy B covers less
                ("c1", "f02:A", 1, 900)]                                    # another family: ignored
        self.assertEqual(T.best_copy_prefix(hsps, "f01"), 2)               # first HSP of A starts at query base 3: two bases before it

    def test_no_hsp_for_the_family_gives_none(self):
        self.assertIsNone(T.best_copy_prefix([("c1", "f02:A", 1, 900)], "f01"))
        self.assertIsNone(T.best_copy_prefix([], "f01"))

    def test_hsp_starting_at_base_one_has_no_prefix(self):
        self.assertEqual(T.best_copy_prefix([("c1", "f01:A", 1, 500)], "f01"), 0)


class Decision(unittest.TestCase):
    def test_trims_a_short_pure_g_prefix(self):
        self.assertEqual(T.decide("GGACGTACGT", 2), ("ACGTACGT", 2, "trimmed"))

    def test_leaves_longer_non_g_or_missing_prefixes(self):
        self.assertEqual(T.decide("GGGGACGT", 4), ("GGGGACGT", 0, "too_long"))
        self.assertEqual(T.decide("GACACGT", 2), ("GACACGT", 0, "not_g"))
        self.assertEqual(T.decide("ACGT", 0), ("ACGT", 0, "no_prefix"))
        self.assertEqual(T.decide("ACGT", None), ("ACGT", 0, "no_hsp"))

    def test_boundary_three_is_trimmed(self):
        self.assertEqual(T.decide("GGGACGT", 3), ("ACGT", 3, "trimmed"))


class PafPrefix(unittest.TestCase):
    def paf(self, q, t, qstart, matches):
        return "\t".join([q, "2500", str(qstart), "2500", "+", t, "9000", "100", "5000", str(matches), "2500", "60"])

    def test_best_record_by_matches_among_the_family_gives_the_prefix(self):
        lines = [self.paf("c1", "f01:A", 2, 2400), self.paf("c1", "f01:B", 0, 2300), self.paf("c1", "f02:A", 5, 2499)]
        self.assertEqual(T.prefixes_from_paf(lines, {"c1": "f01"}), {"c1": 2})

    def test_cluster_names_drop_the_n_suffix_and_unaligned_clusters_are_none(self):
        lines = [self.paf("c1|n=80", "f01:A", 1, 2400)]
        self.assertEqual(T.prefixes_from_paf(lines, {"c1": "f01", "c2": "f01"}), {"c1": 1, "c2": None})

    def test_clusters_without_a_family_are_skipped(self):
        self.assertEqual(T.prefixes_from_paf([self.paf("c1", "f01:A", 1, 2400)], {"c1": None}), {})


class G5Reads(unittest.TestCase):
    def fq(self, n):
        out = []
        for i in range(n):
            out += [f"@r{i}", "ACGTACGTAC", "+", "IIIIIIIIII"]
        return out

    def test_most_reads_get_a_g_run_with_matching_quality(self):
        out = W.add_g5(self.fq(4000), seed=11)
        carrying, lens = 0, collections.Counter()
        for i in range(0, len(out), 4):
            seq, qual = out[i + 1], out[i + 3]
            self.assertEqual(len(seq), len(qual))
            n = 0
            while seq[n] == "G":
                n += 1
            extra = len(seq) - 10
            self.assertEqual(seq[extra:], "ACGTACGTAC")
            if extra:
                carrying += 1
                lens[extra] += 1
                self.assertEqual(seq[:extra], "G" * extra)
        self.assertTrue(0.93 <= carrying / 4000 <= 0.97)
        tot = sum(lens.values())
        self.assertTrue(0.18 <= lens[1] / tot <= 0.26 and 0.38 <= lens[2] / tot <= 0.48 and 0.20 <= lens[3] / tot <= 0.28)
        self.assertTrue(max(lens) <= 6)

    def test_deterministic(self):
        self.assertEqual(W.add_g5(self.fq(50), 3), W.add_g5(self.fq(50), 3))


if __name__ == "__main__":
    unittest.main()
