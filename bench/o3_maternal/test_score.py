#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_score.py"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score as S  # noqa: E402

DELTA = 0.00958


def L(locus, kind="catalog", family="F1", chrom="CM1", start=100, end=200, n=50, ident=0.91):
    return dict(locus=locus, kind=kind, family=family, chrom=chrom, start=start, end=end, n=n, ident=ident)


def cand(name, fam, n, contigs=None):
    return dict(candidate=name, family=fam, n_transcripts=n, contigs=contigs or [name])


def hit(score, acc, s, e):
    return (score, acc, s, e, score)


class Recovery(unittest.TestCase):
    def test_recovers_by_overlap_at_the_floor(self):
        pat = {"c1": hit(0.9995, "CM1", 120, 180)}
        self.assertEqual(S.recovered_loci(["c1"], pat, [L("A")]), {"A"})

    def test_below_the_floor_or_elsewhere_does_not(self):
        self.assertEqual(S.recovered_loci(["c1"], {"c1": hit(0.998, "CM1", 120, 180)}, [L("A")]), set())
        self.assertEqual(S.recovered_loci(["c1"], {"c1": hit(0.9995, "CM2", 120, 180)}, [L("A")]), set())

    def test_class(self):
        self.assertEqual(S.candidate_class(["c"], {}, {}, {"A"}), "a_recovered")
        self.assertEqual(S.candidate_class(["c"], {"c": hit(1.0, "M", 1, 2)}, {"c": hit(1.0, "P", 1, 2)}, set()), "ref")
        self.assertEqual(S.candidate_class(["c"], {"c": hit(0.99, "M", 1, 2)}, {"c": hit(1.0, "P", 1, 2)}, set()), "b_other")
        self.assertEqual(S.candidate_class(["c"], {"c": hit(0.95, "M", 1, 2)}, {}, set()), "c_unmatched")
        self.assertEqual(S.candidate_class(["c"], {}, {}, set()), "c_unmatched")


class Matching(unittest.TestCase):
    def test_kuhn(self):
        self.assertEqual(S.max_matching([(0, 0), (1, 0), (1, 1)], 2), 2)
        self.assertEqual(S.max_matching([(0, 0), (1, 0)], 2), 1)
        self.assertEqual(S.max_matching([], 0), 0)


class Evaluate(unittest.TestCase):
    loci = [L("BIG", n=281, ident=0.91), L("SMALL", n=13, ident=0.92), L("NEAR", n=77, ident=0.995),
            L("SEX", kind="sex", family="LRPAP1", chrom="CMY", n=30, ident=0.98),
            L("LRPAP1_p12", kind="lrpap1", family="LRPAP1", chrom="CM12", n=83, ident=0.9984)]
    fams = ["F1", "F2", "F3", "LRPAP1"]

    def run_eval(self, cands, pat):
        return S.evaluate(cands, self.loci, self.fams, {}, pat)

    def test_r1_passes_when_the_big_beyond_delta_locus_is_recovered(self):
        out = self.run_eval([cand("k1", "F1", 3)], {"k1": hit(1.0, "CM1", 120, 180)})
        self.assertEqual(out["R1"]["verdict"], "PASS")
        self.assertEqual(out["R1"]["targets"], ["BIG"])

    def test_empty_candidates_fail_r1(self):
        out = self.run_eval([], {})
        self.assertEqual(out["R1"]["verdict"], "FAIL")
        self.assertEqual(out["R4"]["flagged"], 0)

    def test_one_transcript_is_not_a_flag(self):
        out = self.run_eval([cand("k1", "F1", 1)], {"k1": hit(1.0, "CM1", 120, 180)})
        self.assertEqual(out["R1"]["verdict"], "FAIL")

    def test_r2_flags_a_recovered_within_delta_locus(self):
        pat = {"k1": hit(1.0, "CM1", 120, 180)}
        loci = [L("NEAR", ident=0.995, n=77)]
        out = S.evaluate([cand("k1", "F1", 2)], loci, self.fams, {}, pat)
        self.assertEqual(out["R2"]["verdict"], "FAIL")

    def test_sex_locus_excluded(self):
        out = self.run_eval([cand("k1", "LRPAP1", 2)], {"k1": hit(1.0, "CMY", 120, 180)})
        self.assertNotIn("SEX", out["R1"]["targets"])
        self.assertEqual(out["R4"]["expressed"], 4)            # BIG, SMALL, NEAR, p12; SEX not counted

    def test_p12_not_in_r1_or_r2(self):
        out = self.run_eval([cand("k1", "LRPAP1", 2)], {"k1": hit(1.0, "CM12", 120, 180)})
        self.assertNotIn("LRPAP1_p12", out["R1"]["targets"])
        self.assertEqual(out["R2"]["verdict"], "PASS")

    def test_r3_counts_false_flags_in_families_without_a_truth_locus(self):
        out = self.run_eval([cand("k1", "F2", 2)], {"k1": hit(1.0, "CM9", 1, 9)})
        self.assertEqual(out["R3"]["denominator"], 2)          # F2, F3 (F1 and LRPAP1 hold expressed loci)
        self.assertEqual(out["R3"]["false_families"], ["F2"])
        self.assertEqual(out["R3"]["fraction"], 0.5)


class Descriptive(unittest.TestCase):
    def test_descriptive_locus_is_not_an_absent_locus(self):
        loci = [L("BIG", n=281, ident=0.91), L("DESC", kind="lrpap1_desc", family="LRPAP1", chrom="CM12", n=83, ident=0.9984)]
        out = S.evaluate([], loci, ["F1", "LRPAP1"], {}, {})
        self.assertEqual(out["R4"]["expressed"], 1)                       # BIG only
        self.assertEqual(out["R3"]["denominator"], 1)                     # LRPAP1 holds no absent expressed locus, so it is in the R3 denominator

    def test_candidate_recovering_a_descriptive_locus_is_not_a_false_flag(self):
        loci = [L("DESC", kind="lrpap1_desc", family="LRPAP1", chrom="CM12", n=83, ident=0.9984)]
        out = S.evaluate([cand("k1", "LRPAP1", 2)], loci, ["LRPAP1"], {}, {"k1": hit(1.0, "CM12", 120, 180)})
        self.assertEqual(out["R3"]["false_families"], [])


class Types(unittest.TestCase):
    def test_distinct_targets(self):
        mat = {"s1": hit(1.0, "M14", 10, 90)}
        pat = {"s2": hit(1.0, "P12", 10, 90), "s3": hit(0.9, "P14", 10, 90)}
        t = {"p12@pat": ("pat", "P12", 0, 100), "p14@pat": ("pat", "P14", 0, 100), "p14@mat": ("mat", "M14", 0, 100)}
        self.assertEqual(S.type_hits(mat, pat, t), {"p12@pat": ["s2"], "p14@pat": [], "p14@mat": ["s1"]})


if __name__ == "__main__":
    unittest.main()
