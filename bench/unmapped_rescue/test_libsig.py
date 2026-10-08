#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_libsig.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import libsig as L  # noqa: E402
import trim5g as T  # noqa: E402


def sam(name, flag, cigar, seq, mapq=60, de=0.001):
    return "\t".join([name, str(flag), "chr1", "100", str(mapq), cigar, "*", "0", "0", seq, "*", f"de:f:{de}"])


class Clip(unittest.TestCase):
    def test_forward_read_clip_is_the_left_soft_clip(self):
        self.assertEqual(L.five_prime_clip(0, "3S7=", "GGGACGTACG"), (3, "GGG"))

    def test_reverse_read_clip_is_the_right_soft_clip_reverse_complemented(self):
        # SEQ is stored on the reference strand: the original read began with GG, which appears as CC at the right end
        self.assertEqual(L.five_prime_clip(16, "8=2S", "ACGTACGTCC"), (2, "GG"))

    def test_no_clip(self):
        self.assertEqual(L.five_prime_clip(0, "10=", "ACGTACGTAC"), (0, ""))
        self.assertEqual(L.five_prime_clip(16, "10=", "ACGTACGTAC"), (0, ""))

    def test_hard_clips_are_ignored(self):
        self.assertEqual(L.five_prime_clip(0, "2H8=", "ACGTACGT"), (0, ""))


class Signature(unittest.TestCase):
    def test_counts_pure_single_base_clips_of_1_to_3_among_clean_primary_reads(self):
        lines = [
            sam("a", 0, "2S8=", "GGACGTACGT"),                 # pure G, 2
            sam("b", 16, "8=1S", "ACGTACGTC"),                  # reverse: original starts with G, 1
            sam("c", 0, "10=", "ACGTACGTAC"),                   # no clip
            sam("d", 0, "2S8=", "ATACGTACGT"),                  # mixed clip: other
            sam("e", 0, "4S6=", "GGGGACGTAC"),                  # pure G but 4 long: other
            sam("f", 0, "1S9=", "TACGTACGTA"),                  # pure T, 1
            sam("g", 0, "2S8=", "GGACGTACGT", mapq=5),          # low MAPQ: skipped
            sam("h", 0, "2S8=", "GGACGTACGT", de=0.05),         # divergent: skipped
            sam("i", 256, "2S8=", "GGACGTACGT"),                # secondary: skipped
            "@SQ\tSN:chr1\tLN:1000",
        ]
        s = L.signature(lines)
        self.assertEqual(s["reads"], 6)
        self.assertEqual(s["no_clip"], 1)
        self.assertEqual((s["pure"]["G"], s["pure"]["T"], s["pure"]["A"], s["pure"]["C"]), (2, 1, 0, 0))
        self.assertEqual(s["g_lengths"], {1: 1, 2: 1})
        self.assertEqual(s["other_clip"], 2)

    def test_keep_filter(self):
        lines = [sam("a", 0, "2S8=", "GGACGTACGT"), sam("b", 0, "2S8=", "GGACGTACGT")]
        self.assertEqual(L.signature(lines, keep=lambda n: n == "a")["reads"], 1)


class Gate(unittest.TestCase):
    def test_g_dominated_clips_pass(self):
        ok, p = L.gate({"pure": {"G": 100, "A": 1, "C": 0, "T": 1}})
        self.assertTrue(ok)
        self.assertLess(p, 1e-6)

    def test_balanced_clips_fail(self):
        self.assertFalse(L.gate({"pure": {"G": 30, "A": 30, "C": 30, "T": 30}})[0])

    def test_too_few_clips_fail(self):
        self.assertFalse(L.gate({"pure": {"G": 15, "A": 0, "C": 0, "T": 0}})[0])


class TrimLeadingG(unittest.TestCase):
    def test_trims_a_run_of_1_to_3(self):
        self.assertEqual(T.trim_leading_g("GACGT"), ("ACGT", 1, "trimmed"))
        self.assertEqual(T.trim_leading_g("GGGACGT"), ("ACGT", 3, "trimmed"))

    def test_leaves_longer_runs_and_none(self):
        self.assertEqual(T.trim_leading_g("GGGGACGT"), ("GGGGACGT", 0, "too_long"))
        self.assertEqual(T.trim_leading_g("ACGT"), ("ACGT", 0, "no_run"))


if __name__ == "__main__":
    unittest.main()
