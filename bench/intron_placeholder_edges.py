#!/usr/bin/env python3
"""Intron placeholders: does preserving exon SPACING separate real containments from repeat ones?
Per `docs/PREREG_intron_placeholder_2026-09-20.md`.

The containment problem (§6t7) survives every pairwise signal: 19 TRUE vs 57 FALSE among the pairs the
shipped rule rejects at containment >= 0.90, with identity 0.987 vs 0.989 and a flat precision curve.
§6t8 showed a variation graph's multiplicity channel is only AUC 0.681. Both pointed at STRUCTURE.

This sits between two measured failures:
  * the GENOMIC span drowns the exon signal in intron sequence, and its length asymmetry is the problem;
  * bare exon CONCATENATION destroys the spacing -- r505: "the single-record rule penalises
    CONCATENATED EXONS specifically".

So render each gene four ways and align each all-vs-all with the same minimap2 invocation:

    G      genomic gene body (shipped substrate)
    S      exons concatenated, introns deleted
    P50    exons joined by a fixed 50 bp N spacer
    Pprop  exons joined by an N spacer of min(true intron length, 500)

Usage: intron_placeholder_edges.py --gff G --genome FA --chrom chrN --out PREFIX [--threads 4]
"""
import argparse
import collections
import os
import re
import subprocess

COMP = str.maketrans('ACGTacgtN', 'TGCAtgcaN')


def exon_structure(gff, chrom):
    """span-name -> (strand, [(s,e)] exons of the transcript with most exonic bases)."""
    by_tx = collections.defaultdict(list)
    tx_gene, tx_strand = {}, {}
    span_of = {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] in ('gene', 'pseudogene'):
            n = re.search(r'Name=([^;]+)', f[8])
            if n:
                span_of[n.group(1)] = f'{f[0]}:{f[3]}-{f[4]}'
        elif f[2] == 'exon':
            p = re.search(r'Parent=([^;]+)', f[8]); g = re.search(r'gene=([^;]+)', f[8])
            if p and g:
                by_tx[p.group(1)].append((int(f[3]), int(f[4])))
                tx_gene[p.group(1)] = g.group(1); tx_strand[p.group(1)] = f[6]
    best = {}
    for tx, segs in by_tx.items():
        g = tx_gene[tx]; n = sum(e - s + 1 for s, e in segs)
        if g not in best or n > best[g][0]:
            best[g] = (n, tx_strand[tx], sorted(segs))
    out = {}
    for g, (n, st, segs) in best.items():
        if g in span_of:
            out[span_of[g]] = (st, segs)
    return out


def render(fa, chrom, st, segs, mode):
    parts = []
    prev_end = None
    for (s, e) in segs:
        if prev_end is not None:
            if mode == 'P50':
                parts.append('N' * 50)
            elif mode == 'Pprop':
                parts.append('N' * min(max(s - prev_end - 1, 0), 500))
        parts.append(fa.fetch(chrom, s - 1, e).upper())
        prev_end = e
    seq = ''.join(parts)
    if st == '-':
        seq = seq.translate(COMP)[::-1]
    return seq


def align(fa_path, paf, threads):
    if os.path.exists(paf):
        return
    with open(paf, 'w') as fh:
        subprocess.run(['minimap2', '-x', 'asm20', '-c', '-X', '-N', '50', '-p', '0.1',
                        '-t', str(threads), fa_path, fa_path], stdout=fh,
                       stderr=subprocess.DEVNULL, check=True)


def pair_stats(paf):
    """(a,b) -> (best identity, best containment)."""
    out = {}
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12 or f[0] == f[5]:
            continue
        la, lb, m, bl = int(f[1]), int(f[6]), int(f[9]), int(f[10])
        k = tuple(sorted((f[0], f[5])))
        ident = m / bl if bl else 0.0
        cont = m / min(la, lb) if min(la, lb) else 0.0
        cur = out.get(k, (0.0, 0.0))
        out[k] = (max(cur[0], ident), max(cur[1], cont))
    return out


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--genome', '--chrom', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--threads', default='4')
    a = ap.parse_args()
    import pysam
    fa = pysam.FastaFile(a.genome)
    struct = exon_structure(a.gff, a.chrom)
    for mode in ('S', 'P50', 'Pprop'):
        path = f'{a.out}.{mode}.fa'
        if not os.path.exists(path):
            with open(path, 'w') as fh:
                for name, (st, segs) in struct.items():
                    seq = render(fa, a.chrom, st, segs, mode)
                    if len(seq) >= 200:
                        fh.write(f'>{name}\n{seq}\n')
        align(path, f'{a.out}.{mode}.paf', a.threads)
        print(f'  {a.chrom} {mode}: {sum(1 for _ in open(path)) // 2} seqs, '
              f'{sum(1 for _ in open(f"{a.out}.{mode}.paf"))} paf records')


if __name__ == '__main__':
    main()
