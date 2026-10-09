#!/usr/bin/env python3
"""Amendment 38: trim the consensus to the columns covered by >= 2 reads. Run with /home/juanfra/miniforge3/bin/python."""
import os
import random
import shutil
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import endtrim as E
import partition as PT
from test_consensus_mode import fragments, identity

try:
    import edlib  # noqa: F401
    import pyabpoa  # noqa: F401
    HAVE = shutil.which("minimap2") is not None
except ImportError:
    HAVE = False


class Window(unittest.TestCase):
    def test_support_window(self):
        self.assertEqual(E.support_window(100, [(0, 100), (10, 90), (20, 80)]), (10, 90))
        self.assertEqual(E.support_window(100, [(0, 100), (0, 100)]), (0, 100))
        self.assertEqual(E.support_window(100, [(0, 100), (50, 100)]), (50, 100))        # one read alone covers 0-50: trimmed

    def test_empty_window_keeps_the_consensus(self):
        self.assertEqual(E.support_window(100, [(0, 100)]), (0, 100))                     # a single aligned read: nothing to compare with
        self.assertEqual(E.support_window(100, []), (0, 100))
        self.assertEqual(E.support_window(100, [(0, 30), (60, 100)]), (0, 100))            # never two reads on a column

    def test_interior_hole_does_not_split(self):
        self.assertEqual(E.support_window(100, [(0, 100), (5, 95), (5, 40), (60, 95)]), (5, 95))

    def test_ref_span_from_cigar(self):
        self.assertEqual(E.ref_span(10, "5S20M2D10M3I5N15M4S"), (10, 10 + 20 + 2 + 10 + 5 + 15))
        self.assertEqual(E.ref_span(0, "100="), (0, 100))


class Flags(unittest.TestCase):
    """Amendment 39: an end is single-read when the extreme read sticks out past the next one by more than the spread of all the others"""

    def test_five_prime(self):
        f = E.end_flags([0, 300, 310, 320], [1000, 1000, 1000, 1000])
        self.assertTrue(f["flag5"])
        self.assertEqual(f["gap5"], 300)
        self.assertFalse(f["flag3"])

    def test_three_prime(self):
        f = E.end_flags([0, 0, 0, 0], [700, 710, 720, 1100])
        self.assertTrue(f["flag3"])
        self.assertEqual(f["gap3"], 380)
        self.assertFalse(f["flag5"])

    def test_spread_ends_are_not_single_read(self):
        f = E.end_flags([0, 40, 90, 160, 250], [900, 960, 1000, 1010, 1100])
        self.assertFalse(f["flag5"])
        self.assertFalse(f["flag3"])

    def test_equal_ends_never_flag(self):
        f = E.end_flags([5, 5, 5], [100, 100, 100])
        self.assertEqual((f["flag5"], f["flag3"], f["gap5"], f["gap3"]), (False, False, 0, 0))

    def test_boundary_is_strict(self):
        self.assertFalse(E.end_flags([0, 100, 150, 200], [10, 10, 10, 10])["flag5"])   # gap 100 vs spread 100
        self.assertTrue(E.end_flags([0, 101, 150, 200], [10, 10, 10, 10])["flag5"])

    def test_fewer_than_three_reads_is_undefined(self):
        f = E.end_flags([0, 300], [1000, 1000])
        self.assertIsNone(f["flag5"])
        self.assertIsNone(f["flag3"])
        self.assertEqual(E.end_flags([], [])["n"], 0)

    def test_order_does_not_matter(self):
        self.assertEqual(E.end_flags([320, 0, 310, 300], [1000] * 4), E.end_flags([0, 300, 310, 320], [1000] * 4))


@unittest.skipUnless(HAVE, "needs minimap2, edlib and pyabpoa")
class Trim(unittest.TestCase):
    def test_unsupported_tail_is_removed(self):
        rnd = random.Random(5)
        truth, reads = fragments(301)
        longest = max(range(len(reads)), key=lambda i: len(reads[i]))
        reads[longest] = reads[longest] + "".join(rnd.choice("ACGT") for _ in range(300))
        cons = PT.abpoa_consensus(reads, mode="l")
        out = E.trim_consensus(cons, reads, "/tmp/endtrim_test")
        self.assertGreater(len(cons), len(out))
        self.assertGreaterEqual(identity(out, truth), 0.99)

    def test_clean_cluster_is_not_harmed(self):
        truth, reads = fragments(302)
        cons = PT.abpoa_consensus(reads, mode="l")
        out = E.trim_consensus(cons, reads, "/tmp/endtrim_test")
        self.assertGreaterEqual(identity(out, truth), 0.99)
        self.assertGreater(len(out), 0.8 * len(cons))


if __name__ == "__main__":
    unittest.main()
