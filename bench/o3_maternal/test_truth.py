#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_truth.py"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import truth as T  # noqa: E402
from common import Rec  # noqa: E402


class Catalog(unittest.TestCase):
    def test_only_pat_rows(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "bonly.tsv")
            open(p, "w").write("family\tlocus\thap\tchrom\tstart\tend\n"
                               "F1\tF1_B0\tpat\tCM1\t10\t50\nF1\tF1_B1\tmat\tCM2\t5\t9\n")
            rows = T.catalog_loci(p, other="pat")
            self.assertEqual([r["locus"] for r in rows], ["F1_B0"])
            self.assertEqual([r["locus"] for r in T.catalog_loci(p, other="mat")], ["F1_B1"])
            self.assertEqual((rows[0]["chrom"], rows[0]["start"], rows[0]["end"], rows[0]["kind"]), ("CM1", 10, 50, "catalog"))


class Lrpap1(unittest.TestCase):
    loci = [("c0", "L0", "NC1", 100, 200), ("c1", "L1", "NC1", 300, 400), ("c2", "L2", "NC2", 100, 200),
            ("c3", "L3", "NCY", 100, 200), ("c4", "L4", "NC3", 100, 200)]
    chrmap = {"NC1": dict(same_hap="pat", same_name="CMP1"), "NC2": dict(same_hap="mat", same_name="CMM2"),
              "NC3": dict(same_hap="pat", same_name="CMP3")}
    lift = {"c0": dict(lift_frac="1.0", B_chrom="CMM1", B_start="1000", B_end="1100"),
            "c1": dict(lift_frac="1.0", B_chrom="CMM1", B_start="2000", B_end="2100"),
            "c2": dict(lift_frac="1.0", B_chrom="CMP2", B_start="50", B_end="150"),
            "c4": dict(lift_frac="0.2", B_chrom="CMM3", B_start="1", B_end="2")}
    hits = {"mat": {"q0": [("CMM1", 1010, 1090, 0.99, 1.0)], "q1": [("CMM9", 5, 50, 0.99, 1.0)], "q4": [("CMM3", 1, 2, 0.99, 1.0)]},
            "pat": {"q2": [("CMP2", 60, 140, 0.99, 1.0)]}}
    qof = {"c0": "q0", "c1": "q1", "c2": "q2", "c3": "q3", "c4": "q4"}

    def iv(self):
        return T.hap_intervals(self.loci, self.chrmap, self.lift, self.hits, self.qof, {("pat", "Y"): "CMY"})

    def test_pat_chromosome_with_a_mat_ortholog(self):
        d = self.iv()["c0"]
        self.assertEqual((d["pat"], d["mat"]), (("CMP1", 100, 200), ("CMM1", 1000, 1100)))

    def test_pat_chromosome_without_a_mat_ortholog_is_mat_absent(self):
        d = self.iv()["c1"]
        self.assertEqual((d["pat"], d["mat"]), (("CMP1", 300, 400), None))

    def test_mat_chromosome_with_a_pat_ortholog(self):
        d = self.iv()["c2"]
        self.assertEqual((d["mat"], d["pat"]), (("CMM2", 100, 200), ("CMP2", 50, 150)))

    def test_chry_is_pat_only_and_sex(self):
        d = self.iv()["c3"]
        self.assertEqual((d["pat"], d["mat"], d["sex"]), (("CMY", 100, 200), None, True))

    def test_chrY_listed_in_chrmap_without_a_partner_is_sex(self):
        chrmap = dict(self.chrmap, NCY=dict(same_hap="pat", same_name="CMPY", B_name=""))
        d = T.hap_intervals([("c3", "L3", "NCY", 100, 200)], chrmap, {}, self.hits, {}, {("pat", "Y"): "CMPY"})["c3"]
        self.assertEqual((d["pat"], d["mat"], d["sex"]), (("CMPY", 100, 200), None, True))

    def test_poor_lift_means_no_counterpart(self):
        self.assertIsNone(self.iv()["c4"]["mat"])

    def test_query_of_picks_the_overlapping_query(self):
        qs = ["NC1:101-200", "NC1:5000-6000", "NC2:101-200"]
        out = T.query_of([("c0", "L0", "NC1", 100, 200), ("c2", "L2", "NC2", 100, 200)], qs)
        self.assertEqual(out, {"c0": "NC1:101-200", "c2": "NC2:101-200"})


class Lrpap1Rows(unittest.TestCase):
    loci = [("c0", "L0", "NC1", 100, 200), ("c1", "L1", "NC1", 300, 400), ("c2", "L2", "NCY", 100, 200),
            ("c3", "L3", "NC1", 500, 600), ("c4", "L4", "NC2", 100, 200)]
    chrmap = {"NC1": dict(same_hap="pat"), "NC2": dict(same_hap="mat")}
    iv = {"c0": dict(sex=False, pat=("P", 100, 200), mat=("M", 1, 2)), "c1": dict(sex=False, pat=("P", 300, 400), mat=None),
          "c2": dict(sex=True, pat=("PY", 100, 200), mat=None), "c3": dict(sex=False, pat=("P", 500, 600), mat=("M", 5, 6)),
          "c4": dict(sex=False, pat=("P2", 1, 2), mat=("M2", 100, 200))}
    lift = {"c0": dict(**{"class": "T2d"}), "c1": dict(**{"class": "T2d"}), "c3": dict(**{"class": "T?"}), "c4": dict(**{"class": "T?"})}

    def kinds(self, ref, other):
        rows = T.lrpap1_rows(self.loci, self.iv, self.lift, self.chrmap, ref, other)
        return {r["locus"]: r["kind"] for r in rows}

    def test_mat_reference(self):
        self.assertEqual(self.kinds("mat", "pat"), {"LRPAP1_c1": "lrpap1", "LRPAP1_c2": "sex", "LRPAP1_c3": "lrpap1_desc"})

    def test_t_question_on_the_reference_sourced_chromosome_is_not_descriptive(self):
        self.assertNotIn("LRPAP1_c4", self.kinds("mat", "pat"))

    def test_row_carries_the_truth_haplotype_interval(self):
        r = [x for x in T.lrpap1_rows(self.loci, self.iv, self.lift, self.chrmap, "mat", "pat") if x["locus"] == "LRPAP1_c3"][0]
        self.assertEqual((r["chrom"], r["start"], r["end"], r["family"]), ("P", 500, 600, "LRPAP1"))


class Labels(unittest.TestCase):
    loci = [dict(locus="F1_B0", kind="catalog", family="F1", chrom="CM1", start=100, end=200),
            dict(locus="LRPAP1_c9", kind="lrpap1", family="LRPAP1", chrom="CM2", start=100, end=200)]

    def rec(self, ref, s, e):
        return Rec(True, 60, 100, 1.0, ref, s, e, 0.0)

    def test_catalog_needs_family_match(self):
        lab = T.label_reads({"a": self.rec("CM1", 120, 180), "b": self.rec("CM1", 120, 180)}, self.loci, {"a": "F1", "b": "F2"}, {})
        self.assertEqual(lab, {"a": "F1_B0", "b": "shared"})

    def test_lrpap1_needs_net_membership(self):
        lab = T.label_reads({"a": self.rec("CM2", 120, 180), "b": self.rec("CM2", 120, 180)}, self.loci, {}, {"a": ["c9"]})
        self.assertEqual(lab, {"a": "LRPAP1_c9", "b": "shared"})

    def test_descriptive_locus_needs_net_membership(self):
        loci = self.loci + [dict(locus="LRPAP1_c3", kind="lrpap1_desc", family="LRPAP1", chrom="CM3", start=100, end=200)]
        lab = T.label_reads({"a": self.rec("CM3", 120, 180), "b": self.rec("CM3", 120, 180)}, loci, {}, {"a": ["c3"]})
        self.assertEqual(lab, {"a": "LRPAP1_c3", "b": "shared"})

    def test_descriptive_locus_is_labelled_by_the_primary_even_when_tied(self):
        loci = self.loci + [dict(locus="LRPAP1_c3", kind="lrpap1_desc", family="LRPAP1", chrom="CM3", start=100, end=200)]
        prim = {"a": self.rec("CM3", 120, 180), "b": self.rec("CM3", 120, 180), "c": self.rec("CM3", 120, 180)}
        lab = T.label_reads({"a": None, "b": None, "c": self.rec("CM7", 1, 9)}, loci, {}, {"a": ["c3"], "c": ["c3"]}, prim)
        self.assertEqual(lab, {"a": "LRPAP1_c3", "b": "ambiguous", "c": "shared"})      # b: not in the LRPAP1 net; c: untied elsewhere

    def test_catalog_locus_never_uses_the_tied_primary(self):
        prim = {"a": self.rec("CM1", 120, 180)}
        lab = T.label_reads({"a": None}, self.loci, {"a": "F1"}, {}, prim)
        self.assertEqual(lab, {"a": "ambiguous"})

    def test_tied_is_ambiguous_and_elsewhere_is_shared(self):
        lab = T.label_reads({"a": None, "b": self.rec("CM7", 1, 9)}, self.loci, {}, {})
        self.assertEqual(lab, {"a": "ambiguous", "b": "shared"})


if __name__ == "__main__":
    unittest.main()
