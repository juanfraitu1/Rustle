#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_t4.py"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import libsig as L  # noqa: E402
import synth_world as W  # noqa: E402
import trim5g as T  # noqa: E402

ART = {0: 0.20, 1: 0.32, 2: 0.38, 3: 0.07, 4: 0.03}


def sam(name, flag, cigar, seq, mapq=60, de=0.001):
    return "\t".join([name, str(flag), "chr1", "100", str(mapq), cigar, "*", "0", "0", seq, "*", f"de:f:{de}"])


def reads_with(j, n, rng, tail="ACGTTACGGATC"):
    """n read sequences: j templated G's then an artifact run drawn from ART, then a non-G start"""
    out = []
    for _ in range(n):
        u, acc, a = rng.random(), 0.0, 0
        for l, w in ART.items():
            acc += w
            a = l
            if u <= acc:
                break
        out.append("G" * (j + a) + tail)
    return out


class Distribution(unittest.TestCase):
    def test_signature_keeps_pure_g_lengths_up_to_eight(self):
        lines = [sam("a", 0, "5S8=", "GGGGGACGTACGT"), sam("b", 0, "2S8=", "GGACGTACGT"), sam("c", 0, "10=", "ACGTACGTAC"), sam("d", 0, "9S8=", "GGGGGGGGGACGTACGT")]
        s = L.signature(lines)
        self.assertEqual(s["g_lengths_all"], {5: 1, 2: 1})
        self.assertEqual(s["g_lengths"], {2: 1})          # the registered <= 3 keys are unchanged
        self.assertEqual(s["pure"]["G"], 1)

    def test_artifact_distribution_is_clip_length_over_clean_reads(self):
        sig = {"reads": 100, "no_clip": 20, "g_lengths_all": {1: 30, 2: 40, 5: 4}}
        a = L.artifact_distribution(sig)
        self.assertAlmostEqual(a[0], 0.20)
        self.assertAlmostEqual(a[2], 0.40)
        self.assertAlmostEqual(a[5], 0.04)
        self.assertEqual(a[3], 1e-4)              # never seen: the floor
        self.assertEqual(sorted(a)[:9], list(range(9)))


class Estimate(unittest.TestCase):
    def test_run_length(self):
        self.assertEqual(T.leading_g_run("GGAC"), 2)
        self.assertEqual(T.leading_g_run("ACGG"), 0)

    def test_recovers_the_templated_count_from_the_reads(self):
        a = L.artifact_distribution({"reads": 100, "no_clip": 20, "g_lengths_all": {1: 32, 2: 38, 3: 7, 4: 3}})
        rng = random.Random(1)
        for j in (0, 1, 2, 3):
            runs = [T.leading_g_run(r) for r in reads_with(j, 120, rng)]
            self.assertEqual(T.estimate_templated(runs, a), j, f"j={j}")

    def test_ties_go_to_the_larger_templated_count(self):
        a = {l: 1e-4 for l in range(9)}
        a.update({0: 0.5, 1: 0.5})
        self.assertEqual(T.estimate_templated([1, 1, 1, 1], a), 1)    # j=0 (artifact 1 each time) and j=1 (artifact 0) are equally likely: pick 1


class Correct(unittest.TestCase):
    A = L.artifact_distribution({"reads": 100, "no_clip": 20, "g_lengths_all": {1: 32, 2: 38, 3: 7, 4: 3}})

    def test_removes_only_the_artifact_part_of_the_run(self):
        rng = random.Random(2)
        reads = reads_with(2, 100, rng)
        cons = "GGGGACGTTACGGATC"          # 2 templated + 2 artifact
        new, t, j = T.correct_cluster(cons, reads, self.A)
        self.assertEqual((j, t, new), (2, 2, "GGACGTTACGGATC"))

    def test_no_templated_run_removes_all(self):
        rng = random.Random(3)
        new, t, j = T.correct_cluster("GGACGTTACGGATC", reads_with(0, 100, rng), self.A)
        self.assertEqual((j, t, new), (0, 2, "ACGTTACGGATC"))

    def test_run_not_longer_than_the_templated_count_is_left(self):
        rng = random.Random(4)
        new, t, j = T.correct_cluster("GGACGTTACGGATC", reads_with(3, 100, rng), self.A)
        self.assertEqual((t, new), (0, "GGACGTTACGGATC"))

    def test_small_clusters_are_left_alone(self):
        new, t, j = T.correct_cluster("GGACGT", ["GGACGT"] * 4, self.A)
        self.assertEqual((t, j), (0, None))


class LeadGWorld(unittest.TestCase):
    def test_families_get_0_to_3_forced_leading_g_in_turn(self):
        specs = W.family_specs(seed=5, n_per_class=4, classes=(0.01, 0.02), lead_g=True)
        self.assertEqual([s["lead_g"] for s in specs], [0, 1, 2, 3, 0, 1, 2, 3])
        self.assertTrue(all("lead_g" not in s for s in W.family_specs(seed=5, n_per_class=2, classes=(0.01,))))


if __name__ == "__main__":
    unittest.main()
