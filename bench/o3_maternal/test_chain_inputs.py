#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_chain_inputs.py"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chain_inputs as I  # noqa: E402


class Loci(unittest.TestCase):
    def test_merge_within_gap(self):
        self.assertEqual(I.merge_loci([("c1", 100, 200), ("c1", 4000, 4500), ("c1", 90000, 91000), ("c2", 10, 20)]),
                         [("c1", 100, 4500), ("c1", 90000, 91000), ("c2", 10, 20)])

    def test_merge_is_order_independent(self):
        a = [("c1", 4000, 4500), ("c1", 100, 200)]
        self.assertEqual(I.merge_loci(a), I.merge_loci(list(reversed(a))))

    def test_panel_mask_first_then_keep(self):
        p = I.build_panel({"F1": [("c1", 100, 200), ("c2", 5, 9)], "F2": [("c1", 1, 2)]}, ["F1", "F2", "F3"])
        self.assertEqual([x["fam"] for x in p], ["F1", "F2"])          # F3 has no mat locus: dropped
        self.assertEqual(p[0]["mask"], ["c1", 100, 200, "F1:0"])
        self.assertEqual(p[0]["keep"], [["c2", 5, 9, "F1:1"]])
        self.assertEqual(p[1]["keep"], [])


class Rename(unittest.TestCase):
    def test_accessions_become_index_names_and_index_names_stay(self):
        m = {"CM1": "chr1_mat_hsa1"}
        self.assertEqual(I.rename_hits([("CM1", 1, 2), ("chr2_mat_hsa3", 3, 4)], m), [("chr1_mat_hsa1", 1, 2), ("chr2_mat_hsa3", 3, 4)])


class Names(unittest.TestCase):
    def test_name_map_pairs_by_order_and_length(self):
        self.assertEqual(I.name_map([("CM1", 10), ("CM2", 20)], [("chr1_mat_hsa1", 10), ("chr2_mat_hsa2", 20)]),
                         {"CM1": "chr1_mat_hsa1", "CM2": "chr2_mat_hsa2"})

    def test_name_map_pairs_by_length_when_the_orders_differ(self):
        self.assertEqual(I.name_map([("CM2", 20), ("CM1", 10)], [("chr1_mat_hsa1", 10), ("chr2_mat_hsa3", 20)]),
                         {"CM1": "chr1_mat_hsa1", "CM2": "chr2_mat_hsa3"})

    def test_name_map_pairs_a_duplicated_length_in_relative_order(self):
        self.assertEqual(I.name_map([("A", 5), ("B", 9), ("C", 5)], [("x", 5), ("y", 5), ("z", 9)]), {"A": "x", "C": "y", "B": "z"})

    def test_name_map_refuses_a_length_mismatch(self):
        with self.assertRaises(AssertionError):
            I.name_map([("CM1", 10)], [("chr1_mat_hsa1", 11)])


if __name__ == "__main__":
    unittest.main()
