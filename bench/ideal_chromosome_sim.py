#!/usr/bin/env python3
"""Ideal-scenario chromosome simulation, per `docs/PREREG_ideal_chromosome_sim_2026-09-21.md`
(md5 `ff226f41`).

Establishes the CEILING of node construction under ideal input, with single-copy genes as the negative
control. §6n0 simulated ONE family (26 NPIP transcripts) and scored transcript completeness; this
simulates a whole chromosome and scores the node/family endpoint, with a single-copy stratum.

Arms:
  ideal  full-length reads, ends jittered, no readthrough      -> the ceiling
  trunc  + 5' degradation                                      -> cost of truncation
  rt     + readthrough molecules at the MEASURED 7.51%         -> cost of readthrough

⚠Jitter is mandatory: identical (chrom,pos,CIGAR) reads collapse under dedup (§6n0).
⚠The readthrough rate is measured, not invented: 33,058 of 439,985 real primary MAPQ-60 chr16 reads
 (7.51%) hit the exons of >= 2 distinct annotated genes at >= 25 bp.
"""
import argparse
import collections
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_reads import simulate_reads  # noqa: E402

TX = {'mRNA', 'transcript', 'ncRNA', 'lncRNA', 'pseudogenic_transcript', 'primary_transcript',
      'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'miRNA', 'misc_RNA',
      'V_gene_segment', 'C_gene_segment', 'J_gene_segment'}
RC = str.maketrans('ACGTacgtN', 'TGCAtgcaN')


def load_transcripts(gff, chrom):
    """transcript_id -> (gene_name, strand, [exons]) for one chromosome."""
    gene_of, t2g, ex, strand = {}, {}, collections.defaultdict(list), {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] in ('gene', 'pseudogene'):
            n = re.search(r'Name=([^;]+)', f[8]); i = re.search(r'ID=([^;]+)', f[8])
            if n and i:
                gene_of[i.group(1)] = n.group(1)
        elif f[2] in TX:
            i = re.search(r'ID=([^;]+)', f[8]); p = re.search(r'Parent=([^;,]+)', f[8])
            if i and p and p.group(1) in gene_of:
                t2g[i.group(1)] = gene_of[p.group(1)]; strand[i.group(1)] = f[6]
        elif f[2] == 'exon':
            p = re.search(r'Parent=([^;,]+)', f[8])
            if p and p.group(1) in t2g:
                ex[p.group(1)].append((int(f[3]), int(f[4])))
    out = {}
    for t, e in ex.items():
        e.sort()
        if e:
            out[t] = (t2g[t], strand[t], e)
    return out


def spliced(fa, chrom, exons, strand):
    s = ''.join(fa.fetch(chrom, a - 1, b) for a, b in exons).upper()
    return s.translate(RC)[::-1] if strand == '-' else s


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--genome', '--chrom', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--arm', required=True, choices=['ideal', 'trunc', 'rt'])
    ap.add_argument('--reads-per-tx', type=int, default=10)
    ap.add_argument('--err', type=float, default=0.001)
    ap.add_argument('--jitter', type=int, default=30, help='max bp trimmed from EACH end (mandatory)')
    ap.add_argument('--rt-rate', type=float, default=0.0751, help='MEASURED on real chr16 reads')
    ap.add_argument('--seed', type=int, default=17)
    a = ap.parse_args()

    import pysam
    fa = pysam.FastaFile(a.genome)
    tx = load_transcripts(a.gff, a.chrom)
    print(f'{a.chrom}: {len(tx)} transcripts over '
          f'{len({v[0] for v in tx.values()})} genes', file=sys.stderr)

    # neighbouring same-strand transcript, for readthrough molecules
    order = sorted(tx.items(), key=lambda kv: kv[1][2][0][0])
    nxt = {}
    for i, (t, (g, s, e)) in enumerate(order):
        for j in range(i + 1, min(i + 25, len(order))):
            t2, (g2, s2, e2) = order[j]
            if g2 != g and s2 == s and e2[0][0] > e[-1][1]:
                nxt[t] = t2; break

    n_rt = 0
    seqs = {}
    with open(a.out, 'w') as fh:
        for idx, (t, (g, s, e)) in enumerate(order):
            if t not in seqs:
                seqs[t] = spliced(fa, a.chrom, e, s)
            base = seqs[t]
            if len(base) < 120:
                continue
            for k in range(a.reads_per_tx):
                rng = random.Random((a.seed * 7919) ^ (idx * 104729) ^ (k * 1299709))
                body = base
                tag = 'fl'
                if a.arm == 'rt' and nxt.get(t) and rng.random() < a.rt_rate:
                    t2 = nxt[t]
                    if t2 not in seqs:
                        seqs[t2] = spliced(fa, a.chrom, tx[t2][2], tx[t2][1])
                    if len(seqs[t2]) >= 120:
                        body = base + seqs[t2]; tag = 'rt'; n_rt += 1
                # MANDATORY end jitter -- identical reads collapse under dedup (§6n0)
                lo = rng.randint(0, a.jitter); hi = rng.randint(0, a.jitter)
                body = body[lo:len(body) - hi] if len(body) - hi > lo + 100 else body
                tf = 0.30 if a.arm == 'trunc' else 0.0
                for ri, (rd, _q) in enumerate(simulate_reads(body, 1, err=a.err, indel=a.err / 3,
                                                             seed=(idx * 131 + k), trunc_frac=tf)):
                    fh.write(f'>{t}|{g}|{tag}|{k}\n{rd}\n')
    print(f'arm={a.arm}  wrote {a.out}  (readthrough molecules: {n_rt})', file=sys.stderr)


if __name__ == '__main__':
    main()
