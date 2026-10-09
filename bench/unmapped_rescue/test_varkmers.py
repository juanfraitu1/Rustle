#!/usr/bin/env python3
"""Amendment 43b: variant k-mers and their DNA label."""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import varkmers as V

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def slow(seq, k=21):
    out = []
    for i in range(len(seq) - k + 1):
        w = seq[i:i + k]
        if any(c not in CODE for c in w):
            continue
        fw = 0
        rv = 0
        for c in w:
            fw = (fw << 2) | CODE[c]
        for c in reversed(w):
            rv = (rv << 2) | (3 - CODE[c])
        out.append((i, min(fw, rv)))
    return out


class Codes(unittest.TestCase):
    def test_matches_the_c_tool_encoding(self):
        rnd = random.Random(1)
        s = "".join(rnd.choice("ACGT") for _ in range(300)) + "N" + "".join(rnd.choice("ACGT") for _ in range(50))
        pos, codes = V.kmer_codes(s)
        self.assertEqual(list(zip(pos.tolist(), codes.tolist())), slow(s))

    def test_reverse_complement_is_the_same_canonical_kmer(self):
        s = "ACGTTGCAAGGCTTACGATCGAT"
        rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        self.assertEqual(sorted(V.kmer_codes(s)[1].tolist()), sorted(V.kmer_codes(rc)[1].tolist()))

    def test_short_sequence(self):
        self.assertEqual(len(V.kmer_codes("ACGT")[1]), 0)


class Junctions(unittest.TestCase):
    def test_query_junctions_from_cigar(self):
        self.assertEqual(V.query_junctions("100M500N50M2I30M1000N20M"), [100, 182])
        self.assertEqual(V.query_junctions("10S100M"), [])

    def test_kmers_across_a_junction_are_dropped(self):
        rnd = random.Random(2)
        s = "".join(rnd.choice("ACGT") for _ in range(100))
        pos, codes = V.kmer_codes(s)
        keep = V.not_across(pos, [50])
        self.assertTrue(all(p + 21 <= 50 or p >= 50 for p in pos[keep].tolist()))
        self.assertEqual(int(keep.sum()), len(pos) - 20)


class Label(unittest.TestCase):
    def test_rule(self):
        self.assertEqual(V.label([2, 3, 5, 0, 0, 7, 9, 4, 2, 6]), "DNA-SUPPORTED")      # 8 of 10 seen >= 2
        self.assertEqual(V.label([0] * 9 + [5]), "RNA-ONLY")                             # 1 of 10
        self.assertEqual(V.label([0] * 5 + [5] * 5), "DNA-SUPPORTED")                    # exactly half
        self.assertEqual(V.label([0] * 6 + [5] * 4), "UNDECIDED")
        self.assertEqual(V.label([5] * 9), "UNDECIDED")                                  # fewer than 10 variant k-mers
        self.assertEqual(V.label([1] * 20), "RNA-ONLY")                                  # seen once is not seen


if __name__ == "__main__":
    unittest.main()
