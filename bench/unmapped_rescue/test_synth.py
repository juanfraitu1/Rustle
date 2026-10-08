#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_synth.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import synth_world as W  # noqa: E402


class Specs(unittest.TestCase):
    def setUp(self):
        self.specs = W.family_specs(seed=7, n_per_class=2, classes=(0.005, 0.02, 0.08))

    def test_counts_and_classes(self):
        self.assertEqual(len(self.specs), 6)
        self.assertEqual([s["divergence_class"] for s in self.specs], [0.005, 0.005, 0.02, 0.02, 0.08, 0.08])
        self.assertEqual(len({s["name"] for s in self.specs}), 6)

    def test_three_copies_one_erased_with_the_class_rate_and_unique_ids(self):
        for s in self.specs:
            ids = [c["id"] for c in s["copies"]]
            self.assertEqual(len(set(ids)), 3)
            self.assertTrue(all(i.startswith(s["name"]) for i in ids))
            a, b, e = s["copies"]
            self.assertEqual(a.get("ops", []), [])
            self.assertEqual(b["ops"][0]["rate"], 1.5 * s["divergence_class"])
            self.assertEqual(e["ops"][0]["rate"], s["divergence_class"])
            self.assertIs(e["in_reference"], False)
            self.assertNotIn("in_reference", a)

    def test_gene_structure_within_the_registered_ranges(self):
        for s in self.specs:
            ex, intr = s["template"]["exons"], s["template"]["introns"]
            self.assertTrue(5 <= len(ex) <= 10)
            self.assertEqual(len(intr), len(ex) - 1)
            self.assertTrue(100 <= ex[0] <= 250 and 250 <= ex[-1] <= 800)
            self.assertTrue(all(90 <= x <= 350 for x in ex[1:-1]))
            self.assertTrue(all(300 <= x <= 2500 for x in intr))

    def test_deterministic_and_seed_dependent(self):
        self.assertEqual(W.family_specs(7, 2, (0.005, 0.02, 0.08)), self.specs)
        self.assertNotEqual(W.family_specs(8, 2, (0.005, 0.02, 0.08)), self.specs)

    def test_reads_block_matches_the_registered_model(self):
        r = self.specs[0]["reads"]
        self.assertEqual((r["per_copy"], r["err"], r["indel"], r["jitter"]), (80, 0.001, 0.0003, 30))


class PolyA(unittest.TestCase):
    def test_every_read_gets_20_to_30_a_and_matching_quality(self):
        fq = ["@r1", "ACGT", "+", "IIII", "@r2", "GGGG", "+", "IIII"]
        out = W.add_polya(fq, seed=3)
        self.assertEqual(len(out), 8)
        for i in (0, 4):
            seq, qual = out[i + 1], out[i + 3]
            tail = len(seq) - 4
            self.assertTrue(20 <= tail <= 30)
            self.assertEqual(seq[4:], "A" * tail)
            self.assertEqual(len(qual), len(seq))

    def test_deterministic(self):
        fq = ["@r1", "ACGT", "+", "IIII"]
        self.assertEqual(W.add_polya(fq, 3), W.add_polya(fq, 3))


class Ends(unittest.TestCase):
    def test_exact_alignment_has_zero_offsets(self):
        import ends as E
        self.assertEqual(E.end_offsets([("=", 100)]), (0, 0))

    def test_consensus_extra_bases_at_the_5_prime_end_are_positive(self):
        import ends as E
        # 6 consensus bases that the transcript lacks, then 100 matches, then 4 consensus bases the transcript lacks
        self.assertEqual(E.end_offsets([("I", 6), ("=", 100), ("I", 4)]), (6, 4))

    def test_consensus_missing_transcript_bases_are_negative(self):
        import ends as E
        self.assertEqual(E.end_offsets([("D", 12), ("=", 100), ("D", 3)]), (-12, -3))

    def test_a_short_mismatch_run_inside_the_end_does_not_move_the_anchor(self):
        import ends as E
        # leading: 2 matches, 1 mismatch, 2 matches (no run of 8 yet), 50 matches: anchor after 5 columns consuming 5 on both sides -> offset 0
        self.assertEqual(E.end_offsets([("=", 2), ("X", 1), ("=", 2), ("=", 50)]), (0, 0))

    def test_net_offset_counts_both_sides_of_the_anchor(self):
        import ends as E
        # 3 inserted, 1 mismatch, 2 deleted before the anchor: consensus consumed 3+1, transcript 1+2 -> net +1
        self.assertEqual(E.end_offsets([("I", 3), ("X", 1), ("D", 2), ("=", 100)])[0], 1)

    def test_cigar_string_parse(self):
        import ends as E
        self.assertEqual(E.parse_edlib_cigar("3=1X2I10="), [("=", 3), ("X", 1), ("I", 2), ("=", 10)])

    def test_core_identity_ignores_the_end_gaps_before_the_anchors(self):
        import ends as E
        ops = [("D", 20), ("=", 100), ("X", 1), ("=", 99), ("I", 5)]
        self.assertAlmostEqual(E.core_identity(ops), 1 - 1 / 200)

    def test_core_identity_counts_internal_indels(self):
        import ends as E
        ops = [("=", 50), ("I", 2), ("=", 50), ("D", 3), ("=", 45)]
        self.assertAlmostEqual(E.core_identity(ops), 1 - 5 / 150)


if __name__ == "__main__":
    unittest.main()
