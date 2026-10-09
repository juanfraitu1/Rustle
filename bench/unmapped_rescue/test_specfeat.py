#!/usr/bin/env python3
"""Amendment 43: differences of a consensus to the primary, from minimap2's long cs string."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import specfeat as F


class Cs(unittest.TestCase):
    def test_counts(self):
        f = F.cs_features("=ACGT*ag=TTTT+t=CA-g=CC~gt100ag=AAAA")
        self.assertEqual((f["sub"], f["ins"], f["dele"], f["ins_bp"], f["del_bp"]), (1, 1, 1, 1, 1))
        self.assertEqual(f["ag"], 1)                                  # a>g is the editing signature

    def test_homopolymer_indels(self):
        f = F.cs_features("=ACTTT+t=GCA")                           # one more T after TTT
        self.assertEqual(f["hp_indel"], 1)
        f = F.cs_features("=ACG+tt=GCA")                            # TT inserted after G, followed by G: not a homopolymer change
        self.assertEqual(f["hp_indel"], 0)
        f = F.cs_features("=AC-ggg=GGCA")                           # GGG deleted before GG: homopolymer
        self.assertEqual(f["hp_indel"], 1)
        f = F.cs_features("=AAC-g=TAA")                             # single G deleted between C and T: no
        self.assertEqual(f["hp_indel"], 0)

    def test_editing_signature_both_strands(self):
        self.assertEqual(F.cs_features("=AA*ag=CC*tc=GG*ac=TT")["ag"], 2)   # a>g and t>c; a>c is not
        self.assertEqual(F.cs_features("=AA*ag=CC*tc=GG*ac=TT")["sub"], 3)

    def test_matches_and_intron(self):
        f = F.cs_features("=ACGT~gt50ag=AC")
        self.assertEqual((f["match"], f["sub"], f["ins"], f["dele"]), (6, 0, 0, 0))


if __name__ == "__main__":
    unittest.main()
