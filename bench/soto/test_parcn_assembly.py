#!/usr/bin/env python3
"""Unit tests for bench/soto/parcn_assembly.py (docs/archive/2026-09/PREREG_soto_parcn_assembly_2026-09-29.md §2, C0(a)).

    python3 bench/soto/test_parcn_assembly.py        (stdlib unittest; numpy; meryl / samtools tests skip if absent)

Every k-mer function is checked against a plain-Python brute force on a random soft-masked sequence with N runs
(the C0(a) design of the prereg, at toy scale): canonical encoding, exact counts (all-A k-mer skipped), the 90-
neighbour edit depth, the region-level parCN / famCN rule on hand-made count arrays, the meryl `count` route on a
tiny FASTA, and the region k-mer extraction through samtools faidx.
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parcn_assembly as pa  # noqa: E402

K = pa.K
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rand_seq(rng, n, lower=0.2, n_runs=2):
    s = [rng.choice("ACGT") for _ in range(n)]
    for i in range(n):
        if rng.random() < lower:
            s[i] = s[i].lower()
    for _ in range(n_runs):
        p = rng.randrange(0, n - 5)
        s[p:p + 3] = ["N", "n", "N"]
    return "".join(s)


def py_code(kmer):
    """Forward 2-bit code (A0 C1 G2 T3) of an upper-case k-mer string."""
    v = 0
    for ch in kmer:
        v = (v << 2) | "ACGT".index(ch)
    return v


def py_canonical_list(seq):
    """Brute force: (canonical code, touches lower-case) of every ACGT-only window."""
    out = []
    for i in range(len(seq) - K + 1):
        w = seq[i:i + K]
        u = w.upper()
        if any(c not in "ACGT" for c in u):
            continue
        rc = "".join(COMP[c] for c in reversed(u))
        out.append((min(py_code(u), py_code(rc)), any(c.islower() for c in w)))
    return out


def py_counts(seq):
    d = {}
    for code, _ in py_canonical_list(seq):
        d[code] = d.get(code, 0) + 1
    return d


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


class KmerCodes(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(20260929)
        self.seq = rand_seq(self.rng, 3000)

    def test_canonical_kmers_matches_brute_force(self):
        can, msk = pa.canonical_kmers(self.seq, with_mask=True)
        ref = py_canonical_list(self.seq)
        self.assertEqual(len(can), len(ref))
        self.assertEqual([int(x) for x in can], [c for c, _ in ref])
        self.assertEqual([bool(x) for x in msk], [m for _, m in ref])

    def test_revcomp_and_decode_roundtrip(self):
        codes = np.array([c for c, _ in py_canonical_list(self.seq)][:200], dtype=np.uint64)
        chars = pa.decode_codes(codes)
        self.assertTrue(np.array_equal(pa.encode_strings(chars), codes))
        rc = pa.revcomp_codes(codes)
        self.assertTrue(np.array_equal(pa.revcomp_codes(rc), codes))
        for code, row in zip(codes[:20], chars[:20]):
            s = row.tobytes().decode()
            self.assertEqual(py_code(s), int(code))

    def test_stream_counts_skip_zero(self):
        d = tempfile.mkdtemp()
        try:
            recs = [("c1", self.seq), ("c2", rand_seq(self.rng, 1500) + "A" * 40), ("skip", self.seq[:500])]
            fa = os.path.join(d, "g.fa")
            write_fasta(fa, recs)
            ref = {}
            for name, s in recs[:2]:
                for k, v in py_counts(s).items():
                    ref[k] = ref.get(k, 0) + v
            present = sorted(ref)[: 300]
            absent = [int(x) for x in self.rng.sample(range(1, 1 << 60), 100) if x not in ref]
            targets = np.unique(np.array(present + absent + [0], dtype=np.uint64))
            got = pa.stream_kmer_counts(fa, targets, exclude={"skip"}, chunk=700)
            for t, g in zip(targets.tolist(), got.tolist()):
                self.assertEqual(g, 0 if t == 0 else ref.get(t, 0), t)
            self.assertGreater(ref.get(0, 0), 0)   # the all-A k-mer occurs in c2 ...
            self.assertEqual(int(got[targets == 0][0]), 0)   # ... and is skipped, as in QuicK-mer2
        finally:
            shutil.rmtree(d)

    def test_edit_depth_matches_brute_force(self):
        d = tempfile.mkdtemp()
        try:
            genome = rand_seq(self.rng, 4000, lower=0.0, n_runs=0)
            fa = os.path.join(d, "g.fa")
            write_fasta(fa, [("c1", genome)])
            cnt = py_counts(genome)
            cands = np.array(sorted(self.rng.sample(sorted(cnt), 40)), dtype=np.uint64)
            # plant near-duplicates: a genome that contains two 1-substitution neighbours of the first candidate
            got = pa.edit_depth_batch(cands, fa)
            for x, g in zip(cands.tolist(), got.tolist()):
                dsum = 0
                for p in range(K):
                    b = (x >> (2 * p)) & 3
                    for e in (1, 2, 3):
                        y = (x & ~(3 << (2 * p))) | (((b + e) & 3) << (2 * p))
                        yc = min(y, int(pa.revcomp_codes(np.array([y], dtype=np.uint64))[0]))
                        dsum += min(cnt.get(yc, 0), 255)
                self.assertEqual(g, dsum)
            var = pa.variants_1sub(cands[:1])
            self.assertEqual(var.shape, (1, 90))
        finally:
            shutil.rmtree(d)


class RegionRule(unittest.TestCase):
    def test_region_stats_rule(self):
        """120 positions: 110 SPEC (100 with HG002 count 2, 10 with 0), 10 non-SPEC; positions 0-4 soft-masked, so
        FAM = positions 5-119 (115). HG002 (s = 1): parCN = median 2; chimp (s = 2): counts 1 on positions 0-99, 0
        elsewhere -> parCN 2, famCN 2 x median over the present FAM positions (5-99) = 2, p = 95/115."""
        Q = np.arange(1, 200, dtype=np.uint64)
        ii = np.arange(120)
        mm = np.zeros(120, dtype=bool)
        mm[:5] = True                     # 5 soft-masked positions -> FAM = 115
        spec = np.zeros(len(Q), dtype=bool)
        spec[:110] = True                 # 110 SPEC k-mers among the region's positions
        hg = np.zeros(len(Q), dtype=np.int64)
        hg[:100] = 2
        hg[110:120] = 3
        ptr = np.zeros(len(Q), dtype=np.int64)
        ptr[:100] = 1
        chm = np.ones(len(Q), dtype=np.int64)
        C = {"CHM13": chm, "HG002": hg, "PTR": ptr}
        S = {"CHM13": 2, "HG002": 1, "PTR": 2}
        row = pa.region_stats(ii, mm, spec, spec, C, S, ["PTR"])
        self.assertEqual((row["n_spec"], row["n_fam"]), (110, 115))
        self.assertEqual(row["par_HG002"], 2.0)
        self.assertEqual(row["par_CHM13"], 2.0)
        self.assertEqual(row["par_PTR"], 2.0)
        self.assertAlmostEqual(row["parmean_HG002"], 200 / 110)
        self.assertAlmostEqual(row["specpres_HG002"], 100 / 110)
        self.assertAlmostEqual(row["p_PTR"], 95 / 115)
        self.assertEqual(row["fam_PTR"], 2.0)
        self.assertEqual(row["fam_HG002"], 2.0)
        self.assertAlmostEqual(row["p_HG002"], 105 / 115)   # HG002 present on 5-99 and 110-119
        self.assertAlmostEqual(row["r_PTR"], (100 / 110) / (95 / 115))
        # CHM13-private share: SPEC k-mers absent from every non-reference genome = the 10 with HG002 0 and PTR 0
        self.assertAlmostEqual(row["spec_chm13_private"], 10 / 110)
        # below the resolution floor everything is unresolved; below the presence floor famCN is 0 (absent)
        row2 = pa.region_stats(ii[:50], mm[:50], spec, spec, C, S, ["PTR"], min_k=100)
        self.assertTrue(pa.isnan(row2["par_HG002"]) and pa.isnan(row2["fam_PTR"]))
        far = {"CHM13": chm, "HG002": hg, "PTR": np.zeros(len(Q), dtype=np.int64)}
        row3 = pa.region_stats(ii, mm, spec, spec, far, S, ["PTR"])
        self.assertEqual(row3["fam_PTR"], 0.0)
        self.assertEqual(row3["par_PTR"], 0.0)

    def test_margin_bins(self):
        self.assertEqual(pa.margin_bin(-3), "(-1000000000.0, 0.0]")
        self.assertEqual(pa.margin_bin(0.5), "(0.0, 1.0]")
        self.assertEqual(pa.margin_bin(2.0), "(1.0, 2.0]")
        self.assertEqual(pa.margin_bin(50), "(2.0, 1000000000.0]")
        self.assertEqual(pa.margin_bin(float("nan")), "nan")


@unittest.skipUnless(shutil.which("samtools"), "samtools not on PATH")
class RegionKmers(unittest.TestCase):
    def test_region_kmers_through_faidx(self):
        d = tempfile.mkdtemp()
        try:
            rng = random.Random(7)
            g = {"chr1": rand_seq(rng, 2000), "chr2": rand_seq(rng, 1500)}
            fa = os.path.join(d, "g.fa")
            write_fasta(fa, sorted(g.items()))
            subprocess.run(["samtools", "faidx", fa], check=True)
            regions = [dict(rid="S0", kind="S1E", s1e_row=0, gene="A", gene_id="G1", chrom="chr1", start=100, end=700),
                       dict(rid="C0", kind="CTRL", s1e_row=-1, gene="B", gene_id="", chrom="chr2", start=0, end=400)]
            Q, idx, mask, offs = pa.region_kmers(regions, fa)
            ref = py_canonical_list(g["chr1"][100:700]) + py_canonical_list(g["chr2"][0:400])
            self.assertEqual(list(offs), [0, len(py_canonical_list(g["chr1"][100:700])), len(ref)])
            self.assertEqual([int(Q[i]) for i in idx], [c for c, _ in ref])
            self.assertEqual([bool(m) for m in mask], [m for _, m in ref])
            self.assertTrue(np.all(np.diff(Q.astype(np.int64)) > 0))
        finally:
            shutil.rmtree(d)


@unittest.skipUnless(os.path.exists(pa.MERYL) or shutil.which("meryl"), "meryl not installed")
class MerylCount(unittest.TestCase):
    def test_meryl_count_matches_brute_force(self):
        d = tempfile.mkdtemp()
        try:
            rng = random.Random(3)
            recs = [("c1", rand_seq(rng, 2500)), ("c2", rand_seq(rng, 1200) + "A" * 35), ("chrY", rand_seq(rng, 800))]
            fa = os.path.join(d, "g.fa")
            write_fasta(fa, recs)
            subprocess.run(["samtools", "faidx", fa], check=True)
            ref = {}
            for name, s in recs[:2]:
                for k, v in py_counts(s).items():
                    ref[k] = ref.get(k, 0) + v
            ychrs = py_counts(recs[2][1])
            present = sorted(ref)[:250]
            y_only = [k for k in ychrs if k not in ref][:20]
            absent = [int(x) for x in rng.sample(range(1, 1 << 60), 50) if x not in ref and x not in ychrs]
            Q = np.unique(np.array(present + y_only + absent, dtype=np.uint64))
            out = os.path.join(d, "c.i32")
            meryl = pa.MERYL if os.path.exists(pa.MERYL) else shutil.which("meryl")
            got = pa.count_with_meryl(Q, fa, out, os.path.join(d, "w"), exclude={"chrY"}, meryl=meryl, threads=1,
                                      memory=1)
            self.assertTrue(np.array_equal(got, np.fromfile(out, dtype=np.int32)))
            for q, c in zip(Q.tolist(), got.tolist()):
                self.assertEqual(c, ref.get(q, 0), q)
            self.assertGreater(sum(1 for q in y_only if ref.get(q, 0) == 0), 0)  # chrY k-mers really excluded
        finally:
            shutil.rmtree(d)


if __name__ == "__main__":
    unittest.main()
