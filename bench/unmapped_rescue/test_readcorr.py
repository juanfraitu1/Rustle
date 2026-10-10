#!/usr/bin/env python3
"""Amendment 45: IsoCon-style majority correction of a read by its neighbours."""
import os
import random
import shutil
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import readcorr as RC


class Pileup(unittest.TestCase):
    def test_counts_from_long_cs(self):
        cnt, ins, ngap = RC.pileup_cs(8, [(0, "=ACG*tc=ACGT")])
        self.assertEqual(cnt[3], [0, 1, 0, 0, 0])           # C over a T
        self.assertEqual(cnt[0], [1, 0, 0, 0, 0])
        cnt, ins, ngap = RC.pileup_cs(8, [(2, "=GT-ac=GT")])
        self.assertEqual(cnt[4][4], 1)                       # deletion of the target's A at 4
        self.assertEqual(cnt[5][4], 1)
        self.assertEqual(sum(cnt[0]), 0)                     # not covered
        cnt, ins, ngap = RC.pileup_cs(8, [(0, "=ACGT+gg=ACGT")])
        self.assertEqual(ins[4], {"GG": 1})
        self.assertEqual(ngap[4], 1)


class Majority(unittest.TestCase):
    def test_substitution_needs_two_neighbours_and_a_strict_majority(self):
        read = "ACGTACGT"
        one = RC.pileup_cs(8, [(0, "=ACG*tc=ACGT")])
        self.assertEqual(RC.majority(read, one), [])
        two = RC.pileup_cs(8, [(0, "=ACG*tc=ACGT"), (0, "=ACG*tc=ACGT")])
        self.assertEqual(RC.majority(read, two), [("sub", 3, "C")])
        split = RC.pileup_cs(8, [(0, "=ACG*tc=ACGT"), (0, "=ACGTACGT")])
        self.assertEqual(RC.majority(read, split), [])        # 1 of 2 is not more than half

    def test_deletion_and_insertion(self):
        read = "ACGTACGT"
        d = RC.pileup_cs(8, [(0, "=ACGT-a=CGT")] * 3)
        self.assertEqual(RC.majority(read, d), [("del", 4, None)])
        i = RC.pileup_cs(8, [(0, "=ACGT+gg=ACGT")] * 2 + [(0, "=ACGTACGT")])
        self.assertEqual(RC.majority(read, i), [("ins", 4, "GG")])

    def test_correct_read_end_to_end(self):
        self.assertEqual(RC.correct("ACGTACGT", [(0, "=ACG*tc=ACGT")] * 3), "ACGCACGT")


class Minimap(unittest.TestCase):
    def test_noisy_reads_become_identical(self):
        if not shutil.which("minimap2"):
            self.skipTest("minimap2")
        rnd = random.Random(4)
        t = "".join(rnd.choice("ACGT") for _ in range(1500))
        reads = {}
        for i in range(6):
            s = list(t)
            for _ in range(15):                             # 1% random substitutions per read
                p = rnd.randrange(len(s))
                s[p] = rnd.choice([b for b in "ACGT" if b != s[p]])
            reads[f"r{i}"] = "".join(s)
        out = RC.correct_set(reads, "/tmp/readcorr_test")
        before = sum(a != b for a, b in zip(reads["r0"], t))
        after = sum(a != b for a, b in zip(out["r0"], t))
        self.assertEqual(before, 15)
        self.assertLessEqual(after, 1)


if __name__ == "__main__":
    unittest.main()
