#!/usr/bin/env python3
"""Allele-specific junctions: quantify splice junctions whose USAGE depends on the allele a molecule
carries. Long HiFi reads observe BOTH the heterozygous-SNP allele AND the read's junctions directly,
so allele->junction linkage is per-molecule (no statistical phasing). The het loci that CONFOUNDED
copy-detection (task c) are exactly the substrate here: heterozygosity = the phasing signal.

Per gene locus:
  1. call a balanced het ANCHOR SNP from the pileup (biallelic, freq in [0.35,0.65], depth>=D,
     >=C reads per allele) -> partitions reads into 2 haplotypes by the allele they carry.
  2. per read (primary, SEQ): allele at the anchor (CIGAR walk) + its introns (CIGAR N).
  3. per junction j seen >= MIN_J times: among reads carrying each anchor allele AND spanning j's
     locus, PSI = (reads using j)/(reads spanning j). Fisher exact on the 2x2 (allele x uses-j);
     effect = |PSI_alleleX - PSI_alleleY|.
ASJ = significant (BH-FDR q<0.05, in the aggregator) AND |dPSI| >= DPSI.

Run with /home/juanfra/miniforge3/bin/python (scipy). --chrom to parallelize. Deterministic.
"""
import argparse
import os
from collections import defaultdict

import numpy as np
import pysam
from scipy.stats import fisher_exact

META = "/tmp/gene_reps_gw.meta.tsv"
BASES = set("ACGT")
D_MIN = 12          # anchor SNP min depth
C_MIN = 5           # anchor SNP min reads per allele
BAL_LO, BAL_HI = 0.35, 0.65   # anchor SNP balanced-frequency window (het, not error/edit-skew)
MIN_J = 3          # junction observed >= this many times to be tested
MIN_SPAN = 5       # >= this many spanning reads PER ALLELE to test a junction
DPSI = 0.30        # report |dPSI| >= this (effect size); FDR applied in the aggregator


def allele_at(aln, ref_pos):
    qpos, rpos = 0, aln.reference_start
    seq = aln.query_sequence
    for op, length in aln.cigartuples:
        if op in (0, 7, 8):           # M/=/X consume both
            if rpos <= ref_pos < rpos + length:
                return seq[qpos + (ref_pos - rpos)]
            qpos += length; rpos += length
        elif op in (2, 3):            # D/N consume ref
            if rpos <= ref_pos < rpos + length:
                return None
            rpos += length
        elif op in (1, 4):            # I/S consume query
            qpos += length
    return None


def introns(aln):
    out = []
    rpos = aln.reference_start
    for op, length in aln.cigartuples:
        if op == 3:
            out.append((rpos, rpos + length)); rpos += length
        elif op in (0, 2, 7, 8):
            rpos += length
    return out


def call_anchor(bam, chrom, start, end):
    """Return (pos, alleleX, alleleY) of the most-balanced, best-supported het SNP, or None."""
    cov = bam.count_coverage(chrom, start, end, quality_threshold=0,
                             read_callback=lambda r: not r.is_unmapped)
    stack = np.vstack([np.asarray(cov[i], dtype=np.int64) for i in range(4)])  # A,C,G,T
    depth = stack.sum(axis=0)
    srt = np.sort(stack, axis=0)
    major, minor = srt[3], srt[2]
    ok = (depth >= D_MIN) & (minor >= C_MIN) & (minor >= BAL_LO * depth) & (major <= BAL_HI * depth)
    cols = np.nonzero(ok)[0]
    best = None
    AC = "ACGT"
    for c in cols:
        col = stack[:, c]
        order = np.argsort(col)[::-1]
        x, y = AC[order[0]], AC[order[1]]
        bal = int(col[order[1]])      # minor count = phasing power
        if best is None or bal > best[0]:
            best = (bal, start + int(c), x, y)
    if best is None:
        return None
    return best[1], best[2], best[3]


def scan_gene(bam, chrom, start, end):
    anc = call_anchor(bam, chrom, start, end)
    if not anc:
        return []
    pos, ax, ay = anc
    # per read: (allele, introns set, ref_start, ref_end)
    reads = []
    for aln in bam.fetch(chrom, start, end):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.query_sequence is None:
            continue
        if aln.reference_start > pos or aln.reference_end <= pos:
            continue
        b = allele_at(aln, pos)
        if b not in (ax, ay):
            continue
        reads.append((b, set(introns(aln)), aln.reference_start, aln.reference_end))
    if sum(1 for r in reads if r[0] == ax) < C_MIN or sum(1 for r in reads if r[0] == ay) < C_MIN:
        return []
    # candidate junctions: introns seen >= MIN_J times
    jc = defaultdict(int)
    for _, js, _, _ in reads:
        for j in js:
            jc[j] += 1
    out = []
    for j, n in jc.items():
        if n < MIN_J:
            continue
        d, a = j
        cnt = {ax: [0, 0], ay: [0, 0]}   # allele -> [used, spanning]
        for al, js, rs, re in reads:
            if rs <= d and re >= a:       # read spans the junction locus
                cnt[al][1] += 1
                if j in js:
                    cnt[al][0] += 1
        ux, sx = cnt[ax]; uy, sy = cnt[ay]
        if sx < MIN_SPAN or sy < MIN_SPAN:
            continue
        psi_x, psi_y = ux / sx, uy / sy
        # skip junctions constitutive in BOTH alleles (can't be allele-specific) — focuses the
        # multiple-testing burden on alternatively-spliced junctions where ASJ can occur.
        if (psi_x >= 0.98 and psi_y >= 0.98) or (psi_x <= 0.02 and psi_y <= 0.02):
            continue
        eff = abs(psi_x - psi_y)
        _, p = fisher_exact([[ux, sx - ux], [uy, sy - uy]])
        out.append((chrom, pos, ax, ay, d, a, ux, sx, uy, sy,
                    round(psi_x, 3), round(psi_y, 3), round(eff, 3), p))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", default="/home/juanfra/winloci_scratch/GGO.bam")
    ap.add_argument("--chrom", default=None)
    ap.add_argument("--region", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")

    if args.region:
        c, rng = args.region.split(":")
        s, e = (int(x.replace(",", "")) for x in rng.split("-"))
        for r in scan_gene(bam, c, s, e):
            print(r)
        return

    genes = []
    with open(META) as fh:
        fh.readline()
        for line in fh:
            g, c, s, e, strand, ln = line.rstrip("\n").split("\t")
            if args.chrom and c != args.chrom:
                continue
            genes.append((g, c, int(s), int(e)))

    out = args.out or f"/tmp/gw/asj_{args.chrom or 'all'}.tsv"
    rows = []
    n_anchor = 0
    for g, c, s, e in genes:
        r = scan_gene(bam, c, s, e)
        if r:
            n_anchor += 1
        for x in r:
            rows.append((g,) + x)
    with open(out, "w") as fh:
        fh.write("gene\tchrom\tanchor\talleleX\talleleY\tdonor\tacceptor\t"
                 "usedX\tspanX\tusedY\tspanY\tpsiX\tpsiY\tdPSI\tpval\n")
        for x in rows:
            fh.write("\t".join(str(v) for v in x) + "\n")
    sig = sum(1 for x in rows if x[-2] >= DPSI and x[-1] < 0.05)
    print(f"[{args.chrom or 'ALL'}] genes={len(genes)} phaseable={n_anchor} junctions_tested={len(rows)} "
          f"|dPSI|>={DPSI}&p<0.05={sig} -> {out}")


if __name__ == "__main__":
    main()
