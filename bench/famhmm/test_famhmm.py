#!/usr/bin/env python3
"""Run: pixi run -e bench python3 -B -m unittest bench/famhmm/test_famhmm.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famhmm as F  # noqa: E402

ORFS = [
    ("orf1 source=c1 coords=5..124 length=40 frame=1 desc=x", "M" * 40),
    ("orf2 source=c1 coords=10..189 length=60 frame=2 desc=x", "M" * 60),
]


class Helpers(unittest.TestCase):
    def test_copies_by_family_groups_on_first_field(self):
        recs = [("MCL0|0|chr1:1-9|+|nexon=1", "A"), ("MCL1|0|chr2:1-9|+|nexon=1", "C"), ("MCL0|1|chr1:20-29|-|nexon=1", "G")]
        fams = F.copies_by_family_records(recs)
        self.assertEqual(sorted(fams), ["MCL0", "MCL1"])
        self.assertEqual(len(fams["MCL0"]), 2)

    def test_longest_orf_picks_longest_and_parses_coords(self):
        hdr, frame, off, aa = F.longest_orf(ORFS)
        self.assertTrue(hdr.startswith("orf2"))
        self.assertEqual((frame, off, aa), (2, 9, "M" * 60))

    def test_longest_orf_none_when_empty(self):
        self.assertIsNone(F.longest_orf([]))


if __name__ == "__main__":
    unittest.main()
