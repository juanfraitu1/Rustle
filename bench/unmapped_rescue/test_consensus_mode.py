#!/usr/bin/env python3
"""Amendment 37 S1: fragments of unequal extent. Run with /home/juanfra/miniforge3/bin/python (pyabpoa, edlib); skipped without them."""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import partition as PT

try:
    import edlib
    import pyabpoa  # noqa: F401
    HAVE = True
except ImportError:
    HAVE = False


def fragments(seed, n=12, length=2000, max_start=500, max_cut=200, err=0.005):
    rnd = random.Random(seed)
    truth = "".join(rnd.choice("ACGT") for _ in range(length))
    reads = []
    for _ in range(n):
        s, e = rnd.randint(0, max_start), length - rnd.randint(0, max_cut)
        out = []
        for ch in truth[s:e]:
            if rnd.random() < err:
                kind = rnd.random()
                if kind < 0.6:
                    out.append(rnd.choice([b for b in "ACGT" if b != ch]))
                elif kind < 0.8:
                    out.append(ch + rnd.choice("ACGT"))
                continue
            out.append(ch)
        reads.append("".join(out))
    return truth, reads


def identity(cons, truth):
    r = edlib.align(cons, truth, mode="HW", task="distance")
    return 1 - r["editDistance"] / len(cons)


@unittest.skipUnless(HAVE, "needs edlib and pyabpoa")
class Modes(unittest.TestCase):
    def test_default_is_global_and_the_modes_are_selectable(self):
        truth, reads = fragments(1)
        self.assertEqual(PT.abpoa_consensus(reads), PT.abpoa_consensus(reads, mode="g"))
        for m in ("g", "l", "e"):
            self.assertTrue(PT.abpoa_consensus(reads, mode=m))

    def test_fragments_of_unequal_extent(self):
        ids = {m: [] for m in ("g", "l", "e")}
        for seed in range(1, 9):
            truth, reads = fragments(seed)
            for m in ids:
                ids[m].append(identity(PT.abpoa_consensus(reads, mode=m), truth))
        med = {m: sorted(v)[len(v) // 2] for m, v in ids.items()}
        print("median identity to the truth by mode:", {m: round(v, 4) for m, v in med.items()}, "worst:", {m: round(min(v), 4) for m, v in ids.items()})
        self.assertLess(med["g"], 0.99)                           # S1: global mode is poor on fragments
        self.assertGreaterEqual(max(med["l"], med["e"]), 0.995)   # and a local / extension mode is not


if __name__ == "__main__":
    unittest.main()
