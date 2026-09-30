#!/usr/bin/env python3
"""Unit tests for figures/_lrc.py (%LRC; docs/PREREG_lrc_metric_2026-09-29.md).

    python3 figures/test_lrc.py        (stdlib unittest; the BAM and GTF are written to a temporary directory)

The fixture (SAM POS and GTF coordinates 1-based; the expected values below in 0-based half-open [s, e)):

  reads on c1 (contigs c1 = 1000 bp, c2 = 500 bp; nothing on c2)
    r1  flag 0     POS 101  50M100N50M     blocks [100,150) [250,300)      primary
    r2  flag 256   POS 151  100M           [150,250)                       SECONDARY: not counted
    r5  flag 4     POS 151  *              unmapped (placed)               not counted
    r3  flag 2048  POS 301  50M            [300,350)                       SUPPLEMENTARY: not counted
    r4  flag 0     POS 321  10M5D10M       [320,330) [335,345)             MAPQ 0 primary counts; the deletion does not
    r6  flag 16    POS 401  5S20M3I20M     [400,420) [420,440)             reverse strand counts; clip not, insertion
    r7  flag 0     POS 461  10=1X9=        [460,480)                       = and X count

  union of c1 = [100,150) [250,300) [320,330) [335,345) [400,440) [460,480); 4 primary alignments

  models (GTF exons)                        length   covered (by hand)                          %LRC
    T1 c1 101-150, 251-350                  50+100   50 (r1) + 50 (r1) + 10 + 10 (r4) = 120     0.8
    T2 c1 151-250                           100      0 (r1 skips it with N; r2, r5 excluded)    0
    T3 c1 401-420, 421-440, 431-440         40       40 (r6); overlapping records merged        1.0  (3 exon records)
    T4 c2 1-100                             100      0 (no reads on c2)                         0
    T5 c1 461-480                           20       20 (r7)                                    1.0

  summary of T1-T5: mean 2.8 / 5 = 0.56, median 0.8, > 0.98: 2/5, 0.75-0.98: 1/5, < 0.75: 2/5, = 1: 2/5, = 0: 2/5
  (had r2 counted T2 would be 1.0; had r3, T1 = 150/150; had the deletion, T1 = 125/150)
"""
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _lrc as L  # noqa: E402

SAM = """@HD	VN:1.6	SO:coordinate
@SQ	SN:c1	LN:1000
@SQ	SN:c2	LN:500
r1	0	c1	101	60	50M100N50M	*	0	0	{s100}	*
r2	256	c1	151	0	100M	*	0	0	*	*
r5	4	c1	151	0	*	*	0	0	{s10}	*
r3	2048	c1	301	60	50M	*	0	0	{s50}	*
r4	0	c1	321	0	10M5D10M	*	0	0	{s20}	*
r6	16	c1	401	60	5S20M3I20M	*	0	0	{s48}	*
r7	0	c1	461	60	10=1X9=	*	0	0	{s20}	*
"""

GTF = """# comment line
c1	t	transcript	101	350	.	+	.	gene_id "G1"; transcript_id "T1";
c1	t	exon	101	150	.	+	.	gene_id "G1"; transcript_id "T1";
c1	t	exon	251	350	.	+	.	gene_id "G1"; transcript_id "T1";
c1	t	exon	151	250	.	-	.	gene_id "G2"; transcript_id "T2";
c1	t	exon	401	420	.	+	.	gene_id "G3"; transcript_id "T3";
c1	t	exon	421	440	.	+	.	gene_id "G3"; transcript_id "T3";
c1	t	exon	431	440	.	+	.	gene_id "G3"; transcript_id "T3";
c2	t	exon	1	100	.	+	.	gene_id "G4"; transcript_id "T4";
c1	t	exon	461	480	.	+	.	gene_id "G5"; transcript_id "T5";
"""

EXPECTED = {"T1": (2, 150, 120), "T2": (1, 100, 0), "T3": (3, 40, 40), "T4": (1, 100, 0), "T5": (1, 20, 20)}


class PureTest(unittest.TestCase):
    def test_merge_intervals(self):
        s, e = L.merge_intervals([30, 0, 5, 12, 20], [40, 10, 12, 15, 30])
        self.assertEqual(list(zip(s.tolist(), e.tolist())), [(0, 15), (20, 40)])

    def test_covered_bases(self):
        us, ue = np.array([10, 30]), np.array([20, 40])
        # [0,100) holds 20 union bases; [15,35) holds 5 + 5; [20,30) none; [12,13) one
        self.assertEqual(L.covered_bases([0, 15, 20, 12], [100, 35, 30, 13], us, ue).tolist(), [20, 10, 0, 1])
        self.assertEqual(L.covered_bases([0], [5], np.empty(0, np.int64), np.empty(0, np.int64)).tolist(), [0])

    def test_class_boundaries(self):
        # 49/50 = 0.98 and 3/4 = 0.75 are in the middle class (both ends included); 99/100 is above; 74/100 below
        c = L.classes([50, 4, 100, 100], [49, 3, 99, 74])
        self.assertEqual(c["frac_gt98"].tolist(), [False, False, True, False])
        self.assertEqual(c["frac_75to98"].tolist(), [True, True, False, False])
        self.assertEqual(c["frac_lt75"].tolist(), [False, False, False, True])


class BamTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import pysam
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.tmp.name
        sam = os.path.join(d, "reads.sam")
        with open(sam, "w") as fh:
            fh.write(SAM.format(**{f"s{n}": "A" * n for n in (10, 20, 48, 50, 100)}))
        cls.bam = os.path.join(d, "reads.bam")
        pysam.sort("-o", cls.bam, sam)
        pysam.index(cls.bam)
        cls.gtf = os.path.join(d, "models.gtf")
        with open(cls.gtf, "w") as fh:
            fh.write(GTF)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def unions(self):
        return {c: L.contig_union(self.bam, c, threads=1) for c in ("c1", "c2")}

    def test_union(self):
        u = self.unions()
        us, ue, n = u["c1"]
        self.assertEqual(n, 4)
        self.assertEqual(list(zip(us.tolist(), ue.tolist())),
                         [(100, 150), (250, 300), (320, 330), (335, 345), (400, 440), (460, 480)])
        self.assertEqual((u["c2"][0].size, u["c2"][2]), (0, 0))

    def test_models_by_hand(self):
        u = self.unions()
        scored = L.score_models(L.read_models(self.gtf), lambda c: u[c][:2])
        got = {}
        for c, (ids, nx, length, cov) in scored.items():
            for i, t in enumerate(ids):
                got[t] = (int(nx[i]), int(length[i]), int(cov[i]))
        self.assertEqual(got, EXPECTED)
        lrc = {t: v[2] / v[1] for t, v in got.items()}
        self.assertAlmostEqual(lrc["T1"], 0.8)
        self.assertEqual((lrc["T2"], lrc["T3"], lrc["T4"], lrc["T5"]), (0.0, 1.0, 0.0, 1.0))

    def test_contig_filter(self):
        self.assertEqual(sorted(L.read_models(self.gtf, {"c2"})), ["c2"])

    def test_summary_by_hand(self):
        length = [EXPECTED[t][1] for t in sorted(EXPECTED)]
        covered = [EXPECTED[t][2] for t in sorted(EXPECTED)]
        n, mean, median, gt98, mid, lt75, full, zero = L.summarize(length, covered)
        self.assertEqual(n, 5)
        self.assertAlmostEqual(mean, 0.56)
        self.assertAlmostEqual(median, 0.8)
        self.assertEqual((gt98, mid, lt75, full, zero), (0.4, 0.2, 0.4, 0.4, 0.4))


if __name__ == "__main__":
    unittest.main()
