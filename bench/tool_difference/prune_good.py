#!/usr/bin/env python3
"""Prune the good pool's shadow transcripts (docs/PREREG_o2_pruning_and_eichler_2026-10-12.md section 3).

usage: prune_good.py --good GOOD.gtf --primary PRIMARY.gtf --ref REF.gtf[.gz] --genes A.genes.tsv[,B.genes.tsv] --rule o2|eichler --out OUT.gtf
Drops every GOOD transcript that touches an exon of a gene the rule calls contradicted AND is good-pool-only: its intron chain is
absent from PRIMARY on that chromosome (mono-exon: no PRIMARY exon overlaps it). Genes come from s1_o2.py's genes.tsv (`class` for
O2, `eichler_class` for Eichler)."""
import argparse
import collections
import gzip
import re
import sys


def touches(exons, ivs):
    return any(s < e2 and e > s2 for s, e in exons for s2, e2 in ivs)


def to_drop(good, primary, contradicted_exons):
    """good / primary: {tid: (chrom, sorted exons)}; contradicted_exons: {chrom: [(start, end)]} -> tids to drop"""
    chains = collections.defaultdict(set); pex = collections.defaultdict(list)
    for chrom, ex in primary.values():
        chains[chrom].add(tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)))
        pex[chrom].extend(ex)
    out = set()
    for t, (chrom, ex) in good.items():
        civs = contradicted_exons.get(chrom)
        if not civs or not touches(ex, civs):
            continue
        if len(ex) > 1:
            if tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)) not in chains[chrom]:
                out.add(t)
        elif not touches(ex, pex[chrom]):
            out.add(t)
    return out


def read_gtf(p, chroms=None):
    tx = collections.defaultdict(lambda: [None, []])
    for ln in (gzip.open(p, 'rt') if p.endswith('.gz') else open(p)):
        f = ln.split('\t')
        if len(f) < 9 or f[2] != 'exon' or (chroms is not None and f[0] not in chroms):
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        tx[t][0] = f[0]; tx[t][1].append((int(f[3]) - 1, int(f[4])))
    return {t: (c, sorted(ex)) for t, (c, ex) in tx.items()}


def main():
    ap = argparse.ArgumentParser()
    for a in ('--good', '--primary', '--ref', '--genes', '--rule', '--out'):
        ap.add_argument(a, required=True)
    args = ap.parse_args()
    col = 'class' if args.rule == 'o2' else 'eichler_class'
    bad = set()
    for path in args.genes.split(','):
        with open(path) as fh:
            hdr = fh.readline().rstrip('\n').split('\t'); ci = hdr.index(col)
            for ln in fh:
                f = ln.rstrip('\n').split('\t')
                if f[ci] == 'contradicted':
                    bad.add((f[0], f[1]))
    cex = collections.defaultdict(list)
    for ln in (gzip.open(args.ref, 'rt') if args.ref.endswith('.gz') else open(args.ref)):
        f = ln.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        g = re.search(r'gene_id "([^"]*)"', f[8]).group(1) or t
        if (f[0], g) in bad:
            cex[f[0]].append((int(f[3]) - 1, int(f[4])))
    chroms = set(cex)
    drop = to_drop(read_gtf(args.good, chroms), read_gtf(args.primary, chroms), cex)
    kept = dropped = 0
    with open(args.good) as i, open(args.out, 'w') as o:
        for ln in i:
            m = re.search(r'transcript_id "([^"]+)"', ln)
            if m and m.group(1) in drop:
                dropped += ln.split('\t')[2] == 'transcript'
                continue
            kept += len(ln.split('\t')) > 2 and ln.split('\t')[2] == 'transcript'
            o.write(ln)
    print(f"[prune_good] rule {args.rule}: {len(bad):,} contradicted genes; dropped {dropped:,} transcripts, kept {kept:,}", file=sys.stderr)


if __name__ == '__main__':
    main()
