#!/usr/bin/env python3
"""Where O2's moved reads land (docs/PREREG_o2_population_2026-10-12.md, amendment 1).

usage: o2_moved_loci.py --molecules PREFIX.molecules.tsv --ref REF.gtf[.gz] --gtf NAME=PATH [--gtf ...]
Each moved molecule's winning placement (aligned span) is mapped to the reference genes whose exons it overlaps (none = unannotated).
Reports genes receiving >= 1 / 2 / 5 / 10 moved reads and, for genes with >= 2, which GTFs emit any transcript touching the gene's exons."""
import argparse
import collections
import gzip
import re

BIN = 10_000


def build_index(gene_exons):
    """{gene key (chrom, id): [(start, end)]} -> {(chrom, bin): [(start, end, key)]} over the merged exons"""
    index = collections.defaultdict(list)
    for gk, ex in gene_exons.items():
        merged = []
        for s, e in sorted(ex):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        for s, e in merged:
            for b in range(s // BIN, (e - 1) // BIN + 1):
                index[(gk[0], b)].append((s, e, gk))
    return index


def hits(index, chrom, s, e):
    out = set()
    for b in range(s // BIN, (e - 1) // BIN + 1):
        for xs, xe, gk in index.get((chrom, b), ()):
            if xs < e and xe > s:
                out.add(gk)
    return out


def count_per_gene(spans, index):
    per = collections.Counter(); unannotated = 0
    for chrom, s, e in spans:
        h = hits(index, chrom, s, e)
        if not h:
            unannotated += 1
        for gk in h:
            per[gk] += 1
    return dict(per), unannotated


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--molecules', required=True); ap.add_argument('--ref', required=True)
    ap.add_argument('--gtf', action='append', default=[])
    args = ap.parse_args()
    gene_ex = collections.defaultdict(list)
    for ln in opener(args.ref):
        f = ln.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        g = re.search(r'gene_id "([^"]*)"', f[8]).group(1) or t
        gene_ex[(f[0], g)].append((int(f[3]) - 1, int(f[4])))
    index = build_index(gene_ex)
    spans = []
    n = collections.Counter()
    for ln in open(args.molecules):
        f = ln.rstrip('\n').split('\t')
        if f[0] == 'read':
            continue
        n[f[1]] += 1
        if f[1] == 'moved' and f[8] != '-':
            spans.append((f[8], int(f[9]), int(f[10])))
    per, unannotated = count_per_gene(spans, index)
    print(f"molecules {sum(n.values()):,}: " + ', '.join(f"{k} {v:,}" for k, v in n.most_common()))
    print(f"moved reads {len(spans):,}: on annotated genes {len(spans) - unannotated:,}, unannotated {unannotated:,}")
    for k in (1, 2, 5, 10):
        g = [x for x, c in per.items() if c >= k]
        print(f"  genes with >= {k:2d} moved reads: {len(g):6,}  (holding {sum(per[x] for x in g):,} moved reads)")
    ge2 = {x for x, c in per.items() if c >= 2}
    have = {}
    for spec in args.gtf:
        name, path = spec.split('=', 1)
        hit = set()
        for ln in open(path):
            f = ln.split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            hit |= hits(index, f[0], int(f[3]) - 1, int(f[4])) & ge2
        have[name] = hit
        print(f"  of the {len(ge2):,} genes with >= 2 moved reads, {name} has a transcript touching {len(hit):,}")
    if {'StringTie', 'FLAIR'} <= set(have):
        miss = ge2 - have['StringTie'] - have['FLAIR']
        print(f"  genes with >= 2 moved reads and NO StringTie or FLAIR transcript: {len(miss):,} (holding {sum(per[x] for x in miss):,} moved reads)"
              + (f"; ours (OP_good) touches {len(miss & have['OP_good']):,} of them" if 'OP_good' in have else ''))
        top = sorted(miss, key=lambda x: -per[x])[:15]
        print("    top: " + ', '.join(f"{x[1]} ({per[x]})" for x in top))


if __name__ == '__main__':
    main()
