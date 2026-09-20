#!/usr/bin/env python3
"""Allele-specific junctions with MULTI-SNP HAPLOTYPE phasing (extends the single-anchor version).

Single-anchor phasing requires a read to cover ONE specific het SNP. Multi-SNP phasing assigns each
read to one of the two diploid haplotypes using ALL het SNPs in the gene (read-backed 2-haplotype
clustering), so reads covering any subset of SNPs get phased -> more reach (more junctions reach the
spanning-read minimum) and the phasing QUALITY (how cleanly reads bipartition) is itself a confound
check: a clean 2-way split = diploid het; a messy/>2 split = paralog/artifact.

Per gene: call het SNPs -> read x SNP allele matrix -> deterministic 2-means haplotype phasing ->
per junction, PSI per HAPLOTYPE (Fisher exact). BH-FDR + compare-to-single in asjm_aggregate.py.

Run with /home/juanfra/miniforge3/bin/python (scipy). --chrom to parallelize. Deterministic.
"""
import argparse
import os
from collections import defaultdict

import numpy as np
import pysam
from scipy.stats import fisher_exact

META = "/tmp/gene_reps_gw.meta.tsv"
AC = "ACGT"
D_MIN, C_MIN = 12, 5
BAL_LO, BAL_HI = 0.30, 0.70
MAX_SNP = 40            # cap het SNPs per gene (highest-coverage)
MIN_SNP = 1            # >=1 het SNP; with 1 it degenerates to single-anchor (so this is a SUPERSET)
MIN_J = 3
MIN_SPAN = 5
DPSI = 0.30
AGREE = 0.80           # a read is confidently phased if it agrees with one haplotype at >=this frac
PHASE_QUAL_MIN = 0.70  # >= this fraction of SNP-covering reads confidently phased
TRANSITION = ({"A", "G"}, {"C", "T"})


def het_snps(bam, chrom, start, end):
    cov = bam.count_coverage(chrom, start, end, quality_threshold=0,
                             read_callback=lambda r: not r.is_unmapped)
    stack = np.vstack([np.asarray(cov[i], dtype=np.int64) for i in range(4)])
    depth = stack.sum(axis=0)
    srt = np.sort(stack, axis=0)
    major, minor = srt[3], srt[2]
    ok = (depth >= D_MIN) & (minor >= C_MIN) & (minor >= BAL_LO * depth) & (major <= BAL_HI * depth)
    snps = []
    for c in np.nonzero(ok)[0]:
        col = stack[:, c]
        order = np.argsort(col)[::-1]
        snps.append((int(col[order[1]]), start + int(c), AC[order[0]], AC[order[1]]))  # (minorcount,pos,a,b)
    snps.sort(reverse=True)  # by coverage (minor count) desc
    return [(p, a, b) for _, p, a, b in snps[:MAX_SNP]]


def read_alleles_introns(aln, snp_pos):
    """One CIGAR walk: base at each het SNP position + the read's introns."""
    qpos, rpos = 0, aln.reference_start
    seq = aln.query_sequence
    alleles = {}
    introns = []
    si = 0
    ns = len(snp_pos)
    for op, length in aln.cigartuples:
        if op in (0, 7, 8):                 # M/=/X
            end = rpos + length
            while si < ns and snp_pos[si] < rpos:
                si += 1
            j = si
            while j < ns and snp_pos[j] < end:
                alleles[snp_pos[j]] = seq[qpos + (snp_pos[j] - rpos)]
                j += 1
            qpos += length; rpos += length
        elif op == 3:                        # N
            introns.append((rpos, rpos + length)); rpos += length
        elif op == 2:                        # D
            rpos += length
        elif op in (1, 4):                   # I/S
            qpos += length
    return alleles, introns


def phase(reads, snps):
    """Deterministic 2-means over read allele-vectors. reads: list of dict{pos:base}. snps:[(pos,a,b)].
    Returns hap[read_index] in {0,1,None} (None = not confidently phased) and n_snps_diff."""
    pos_idx = {p: i for i, (p, _, _) in enumerate(snps)}
    k = len(snps)
    # seed: SNP 0 (highest coverage). h0/h1 = its two alleles; others None.
    h0 = [None] * k; h1 = [None] * k
    h0[0], h1[0] = snps[0][1], snps[0][2]
    assign = [None] * len(reads)
    for _ in range(6):
        # assign reads
        for ri, al in enumerate(reads):
            a0 = a1 = cov = 0
            for p, b in al.items():
                i = pos_idx.get(p)
                if i is None:
                    continue
                if h0[i] is not None:
                    cov += 1
                    if b == h0[i]: a0 += 1
                    if h1[i] is not None and b == h1[i]: a1 += 1
            if cov == 0:
                assign[ri] = None
            else:
                assign[ri] = 0 if a0 >= a1 else 1
        # re-estimate consensus haplotypes per SNP per group
        for i in range(k):
            for grp, hh in ((0, h0), (1, h1)):
                tally = defaultdict(int)
                for ri, al in enumerate(reads):
                    if assign[ri] == grp:
                        b = al.get(snps[i][0])
                        if b in (snps[i][1], snps[i][2]):
                            tally[b] += 1
                hh[i] = max(sorted(tally), key=tally.get) if tally else None
    # final confident assignment: agree fraction with chosen hap over covered informative SNPs
    n_diff = sum(1 for i in range(k) if h0[i] is not None and h1[i] is not None and h0[i] != h1[i])
    hap = [None] * len(reads)
    for ri, al in enumerate(reads):
        a0 = a1 = cov = 0
        for p, b in al.items():
            i = pos_idx.get(p)
            if i is None or h0[i] is None or h1[i] is None or h0[i] == h1[i]:
                continue
            cov += 1
            if b == h0[i]: a0 += 1
            elif b == h1[i]: a1 += 1
        if cov == 0:
            continue
        best, frac = (0, a0 / cov) if a0 >= a1 else (1, a1 / cov)
        if frac >= AGREE:
            hap[ri] = best
    return hap, n_diff


def scan_gene(bam, chrom, start, end):
    snps = het_snps(bam, chrom, start, end)
    if len(snps) < MIN_SNP:
        return []
    snp_pos = sorted(p for p, _, _ in snps)   # sorted positions for the single-pass CIGAR walk
    reads_al = []
    reads_meta = []  # (introns set, rstart, rend)
    for aln in bam.fetch(chrom, start, end):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.query_sequence is None:
            continue
        al, introns = read_alleles_introns(aln, snp_pos)
        al = {p: b for p, b in al.items() if b in "ACGT"}
        if not al:
            continue
        reads_al.append(al)
        reads_meta.append((set(introns), aln.reference_start, aln.reference_end))
    if len(reads_al) < 2 * C_MIN:
        return []
    hap, n_diff = phase(reads_al, snps)
    if n_diff < MIN_SNP:
        return []
    phased = [i for i in range(len(hap)) if hap[i] is not None]
    snp_covering = sum(1 for al in reads_al if any(p in {s[0] for s in snps} for p in al))
    if snp_covering == 0:
        return []
    phase_qual = len(phased) / snp_covering
    n0 = sum(1 for i in phased if hap[i] == 0)
    n1 = sum(1 for i in phased if hap[i] == 1)
    if n0 < C_MIN or n1 < C_MIN:
        return []
    n_tv = sum(1 for p, a, b in snps if {a, b} not in TRANSITION)  # transversion het SNPs = genetic
    # candidate junctions
    jc = defaultdict(int)
    for i in phased:
        for j in reads_meta[i][0]:
            jc[j] += 1
    out = []
    for j, nj in jc.items():
        if nj < MIN_J:
            continue
        d, a = j
        cnt = [[0, 0], [0, 0]]  # [hap][used,spanning]
        for i in phased:
            js, rs, re = reads_meta[i]
            if rs <= d and re >= a:
                cnt[hap[i]][1] += 1
                if j in js:
                    cnt[hap[i]][0] += 1
        (u0, s0), (u1, s1) = cnt[0], cnt[1]
        if s0 < MIN_SPAN or s1 < MIN_SPAN:
            continue
        p0, p1 = u0 / s0, u1 / s1
        if (p0 >= 0.98 and p1 >= 0.98) or (p0 <= 0.02 and p1 <= 0.02):
            continue
        _, pv = fisher_exact([[u0, s0 - u0], [u1, s1 - u1]])
        out.append((chrom, len(snps), n_tv, round(phase_qual, 2), d, a,
                    u0, s0, u1, s1, round(p0, 3), round(p1, 3), round(abs(p0 - p1), 3), pv))
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
    out = args.out or f"/tmp/gw/asjm_{args.chrom or 'all'}.tsv"
    rows = []
    nph = 0
    for g, c, s, e in genes:
        r = scan_gene(bam, c, s, e)
        if r:
            nph += 1
        for x in r:
            rows.append((g,) + x)
    with open(out, "w") as fh:
        fh.write("gene\tchrom\tn_snp\tn_tv\tphase_qual\tdonor\tacceptor\t"
                 "used0\tspan0\tused1\tspan1\tpsi0\tpsi1\tdPSI\tpval\n")
        for x in rows:
            fh.write("\t".join(str(v) for v in x) + "\n")
    sig = sum(1 for x in rows if x[-2] >= DPSI and x[-1] < 0.05)
    print(f"[{args.chrom or 'ALL'}] genes={len(genes)} multiSNP_phaseable={nph} "
          f"junctions_tested={len(rows)} |dPSI|>={DPSI}&p<0.05={sig} -> {out}")


if __name__ == "__main__":
    main()
