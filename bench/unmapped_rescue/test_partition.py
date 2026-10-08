#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/unmapped_rescue/test_partition.py"""
import collections
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import partition as PT  # noqa: E402

L = 300
_RNG = random.Random(7)
BASE = "".join(_RNG.choice("ACGT") for _ in range(L))   # one generator: a fresh Random(7) per character made BASE a homopolymer (found by the independent review)
SIB_COLS = [40, 70, 100, 150, 200, 230]


def mutate(seq, cols):
    out = list(seq)
    for c in cols:
        out[c] = {"A": "C", "C": "G", "G": "T", "T": "A"}[out[c]]
    return "".join(out)


def noisy(seq, rng, rate=0.001):
    return "".join(({"A": "C", "C": "G", "G": "T", "T": "A"}[b] if rng.random() < rate else b) for b in seq)


def make_reads(groups, seed=1):
    """groups = [(count, columns mutated)] -> {name: seq}, truth {name: group index}"""
    rng = random.Random(seed)
    reads, truth = {}, {}
    for gi, (n, cols) in enumerate(groups):
        ref = mutate(BASE, cols)
        for i in range(n):
            nm = f"g{gi}r{i}"
            reads[nm] = noisy(ref, rng)
            truth[nm] = gi
    return reads, truth


def align_fn(cons, reads):
    lines = []
    for nm, s in reads.items():
        cig = "".join("1=" if a == b else "1X" for a, b in zip(cons, s))
        lines.append("\t".join([nm, "0", "cons", "1", "60", cig, "*", "0", "0", s, "*"]))
    return lines


def consensus_fn(seqs):
    return "".join(collections.Counter(col).most_common(1)[0][0] for col in zip(*seqs))


class Partition(unittest.TestCase):
    def leaves(self, groups, seed=1):
        reads, truth = make_reads(groups, seed)
        return PT.partition(reads, align_fn, consensus_fn), truth

    def purity(self, leaves, truth):
        return all(len({truth[r] for r in lf["reads"]}) == 1 for lf in leaves)

    def test_two_siblings_split_into_two_pure_leaves(self):
        leaves, truth = self.leaves([(30, []), (30, SIB_COLS)])
        self.assertEqual(len(leaves), 2)
        self.assertTrue(self.purity(leaves, truth))
        self.assertEqual(sorted(len(lf["reads"]) for lf in leaves), [30, 30])
        self.assertEqual(sorted(lf["cons"] for lf in leaves), sorted([BASE, mutate(BASE, SIB_COLS)]))

    def test_unbalanced_members_still_split(self):
        leaves, truth = self.leaves([(50, []), (10, SIB_COLS)])
        self.assertEqual(len(leaves), 2)
        self.assertTrue(self.purity(leaves, truth))

    def test_one_transcript_with_noise_stays_one_leaf(self):
        for seed in (1, 2, 3):
            leaves, _ = self.leaves([(60, [])], seed)
            self.assertEqual(len(leaves), 1, seed)

    def test_a_lone_hotspot_column_is_not_a_block(self):
        leaves, _ = self.leaves([(52, []), (8, [120])])
        self.assertEqual(len(leaves), 1)

    def test_too_few_reads_are_left_alone(self):
        leaves, _ = self.leaves([(3, []), (2, SIB_COLS)])
        self.assertEqual(len(leaves), 1)

    def test_three_members_give_three_leaves(self):
        leaves, truth = self.leaves([(30, []), (30, SIB_COLS), (30, [20, 60, 90, 130, 180, 260])])
        self.assertEqual(len(leaves), 3)
        self.assertTrue(self.purity(leaves, truth))

    def test_variants_near_the_ends_are_not_used(self):
        # the only difference sits in the first 15 columns: not a block
        leaves, _ = self.leaves([(30, []), (30, [2, 5, 9])])
        self.assertEqual(len(leaves), 1)


class Blocks(unittest.TestCase):
    def test_linked_variants_form_a_block_and_noise_does_not(self):
        import numpy as np
        n = 60
        C = np.zeros((6, n), dtype=bool)
        C[0:3, :20] = True          # three variants carried by the same 20 reads
        C[3, [1, 30]] = True        # singletons
        C[4, [5, 40]] = True
        C[5, [50, 51]] = True
        K = np.ones((6, n), dtype=bool)
        cols = [10, 20, 30, 40, 50, 60]
        blocks = PT.find_blocks(C, K, cols)
        self.assertEqual([sorted(b) for b in blocks], [[0, 1, 2]])

    def test_same_column_alleles_are_not_linked(self):
        import numpy as np
        C = np.zeros((2, 40), dtype=bool)
        C[0, :20] = True
        C[1, 20:] = True
        self.assertEqual(PT.find_blocks(C, np.ones((2, 40), dtype=bool), [5, 5]), [])


@unittest.skipUnless(__import__("importlib").util.find_spec("edlib"), "edlib needed (miniforge python)")
class EdlibAlign(unittest.TestCase):
    def test_edlib_unit_cost_prefers_a_messy_alignment_to_a_long_deletion(self):
        # documents the limitation that made Amendment 19's aligner poor: a 100-base skip costs 100 edits, a messy alignment of the read elsewhere costs less, so edlib
        # does not report the deletion run (the spliced aligner of Amendment 20 does). Found by the independent review, which also found this test asserted the opposite
        # on a homopolymer fixture.
        cons = BASE
        read = BASE[20:100] + BASE[200:290]
        f = PT.edlib_align_fn()(cons, {"r": read})[0].split("\t")
        ops = [(int(n), o) for n, o in PT.P.CIG.findall(f[5])]
        self.assertNotIn((100, "D"), ops)

    def test_the_skipped_segment_is_found_as_a_block_of_deletions(self):
        full, skip = BASE, BASE[:120] + BASE[220:]
        reads = {f"a{i}": full for i in range(30)}
        reads.update({f"b{i}": skip for i in range(15)})
        leaves = PT.partition(reads, PT.edlib_align_fn(), consensus_fn)
        self.assertEqual(len(leaves), 2)
        self.assertEqual(sorted(len(lf["reads"]) for lf in leaves), [15, 30])


class SplicedCarry(unittest.TestCase):
    def test_an_n_gap_makes_a_read_carry_the_deletion_only_with_the_flag(self):
        cons = BASE
        line = "\t".join(["r", "0", "cons", "1", "60", "100=50N150=", "*", "0", "0", BASE[:100] + BASE[150:], "*"])
        V = [("del", 120, None)]
        _n, C, K = PT.read_carry([line], cons, V, n_as_del=True)
        self.assertTrue(C[0, 0] and K[0, 0])
        _n, C, K = PT.read_carry([line], cons, V)
        self.assertFalse(C[0, 0])


class Cap(unittest.TestCase):
    def test_keeps_the_variants_carried_by_most_reads(self):
        import numpy as np
        C = np.zeros((4, 10), dtype=bool)
        C[0, :2] = True
        C[1, :6] = True
        C[2, :4] = True
        C[3, :1] = True
        K = np.ones_like(C)
        V = [("sub", 20, "A"), ("sub", 21, "A"), ("sub", 22, "A"), ("sub", 23, "A")]
        V2, C2, K2 = PT.cap_variants(V, C, K, 2)
        self.assertEqual(V2, [("sub", 21, "A"), ("sub", 22, "A")])
        self.assertEqual(C2.shape, (2, 10))
        self.assertEqual(PT.cap_variants(V, C, K, 10)[0], V)


if __name__ == "__main__":
    unittest.main()
