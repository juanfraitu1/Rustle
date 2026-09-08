#!/usr/bin/env python3
"""Annotation-degradation arms for the O1 ablation (PREREG_annotation_ablation_2026-09-07, md5 65efc5b2).
Only gene/pseudogene records and their exons on the given contigs are touched; every other line is dropped
(the three-contig substrate). Deterministic: numpy default_rng(1337).
usage: o1_degrade_gff.py <in.gff> <arm> <out.gff> <contig,...>
arms: A0 (identity) | A1 (+-5kb jitter) | A2 (+-20kb) | A3 (50% dropout) | A4 (span-only exons) | A5 (middle 50%)
"""
import sys, re, collections
import numpy as np
inp, arm, out, contigs = sys.argv[1], sys.argv[2], sys.argv[3], set(sys.argv[4].split(','))
rng = np.random.default_rng(1337)
MINLEN = 500
genes, exons = [], collections.defaultdict(list)   # gene: (chrom, start, end, strand, name, raw_fields)
for l in open(inp):
    if l[0] == '#':
        continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[0] not in contigs:
        continue
    if f[2] in ('gene', 'pseudogene'):
        m = re.search(r'ID=([^;]*)', f[8]); n = re.search(r'gene=([^;]*)', f[8])
        genes.append([f[0], int(f[3]), int(f[4]), f[6], (n.group(1) if n else (m.group(1) if m else '?')), f])
    elif f[2] == 'exon':
        n = re.search(r'gene=([^;]*)', f[8])
        if n:
            exons[n.group(1)].append((int(f[3]), int(f[4])))
# one exon set per gene name: the union of its exon intervals (mcl_families joins exons by gene=)
kept, n_drop = [], 0
for g in genes:
    c, s, e, st, name, f = g
    ex = sorted(set(exons.get(name, [])))
    if arm == 'A0':
        pass
    elif arm in ('A1', 'A2'):
        w = 5000 if arm == 'A1' else 20000
        s = max(1, s + int(rng.integers(-w, w + 1))); e = e + int(rng.integers(-w, w + 1))
    elif arm == 'A3':
        if rng.random() < 0.5:
            n_drop += 1; continue
    elif arm == 'A4':
        ex = [(s, e)]
    elif arm == 'A5':
        L = e - s + 1; q = L // 4; s, e = s + q, e - q
    else:
        sys.exit(f'unknown arm {arm}')
    if e - s + 1 < MINLEN:
        n_drop += 1; continue
    ex = [(max(a, s), min(b, e)) for a, b in ex if b >= s and a <= e]
    if arm != 'A4' and not ex:
        ex = [(s, e)]
    kept.append((c, s, e, st, name, f, ex))
kept.sort(key=lambda r: (r[0], r[1], r[2]))
with open(out, 'w') as o:
    o.write('##gff-version 3\n')
    for c, s, e, st, name, f, ex in kept:
        o.write('\t'.join([c, f[1], f[2], str(s), str(e), f[5], st, f[7], f[8]]) + '\n')
        for i, (a, b) in enumerate(ex):
            o.write('\t'.join([c, f[1], 'exon', str(a), str(b), '.', st, '.', f'ID={name}-e{i};gene={name}']) + '\n')
print(f'{arm}: {len(kept)} records written, {n_drop} dropped, {sum(len(r[6]) for r in kept)} exons')
