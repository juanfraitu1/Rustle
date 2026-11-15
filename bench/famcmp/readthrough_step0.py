#!/usr/bin/env python3
"""PREREG r1623 step 0 (docs/PREREG_cc_readthrough_allele_2026-11-10.md): confined vs crossing reads at fused-locus boundaries
and at every intron of control genes. usage: readthrough_step0.py FUSED_TSV REFSEQ_GTF BAM OUT_TSV"""
import sys, re, gzip, hashlib, collections, pysam
A = lambda k, a: (re.search(k + r' "([^"]+)"', a) or [None, ''])[1]
gene = {}
for l in gzip.open(sys.argv[2], 'rt'):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon': continue
    n = A('gene_name', f[8]) or A('gene_id', f[8])
    gene.setdefault(n, (f[0], f[6], {}))[2].setdefault(A('transcript_id', f[8]), []).append((int(f[3]) - 1, int(f[4])))
bam = pysam.AlignmentFile(sys.argv[3])
def reads(c, s, a, b):
    return [(r.reference_start, r.reference_end) for r in bam.fetch(c, a, b)
            if not (r.is_secondary or r.is_supplementary) and r.is_reverse == (s == '-')]
cut = lambda L, R, X: L >= 2 and R >= 2 and X < min(L, R)
span = lambda n: (min(e[0][0] for e in map(sorted, gene[n][2].values())), max(e[-1][1] for e in map(sorted, gene[n][2].values())))
out = open(sys.argv[4], 'w'); fused = set(); n = k = nb = 0
for l in open(sys.argv[1]):
    f = l.split()
    if f[0] == 'class' or f[4] == 'mono_blob': continue
    g = f[-1].split(',')[:2]; fused.update(f[-1].split(','))
    if any(x not in gene for x in g): out.write(f'fused\t{f[1]}\t{f[-1]}\tno_gene\n'); n += 1; continue
    (a1, b1), (a2, b2) = sorted(span(x) for x in g); c, s, _ = gene[g[0]]; n += 1
    if a2 < b1: nb += 1; out.write(f'fused\t{f[1]}\t{f[-1]}\tno_boundary\n'); continue
    R = reads(c, s, a1, b2)
    L = sum(e <= a2 and e > a1 for st, e in R); Rr = sum(st >= b1 and st < b2 for st, e in R); X = sum(st < b1 and e > a2 for st, e in R)
    k += cut(L, Rr, X); out.write(f'fused\t{f[1]}\t{f[-1]}\t{L}\t{Rr}\t{X}\t{int(cut(L, Rr, X))}\t{f[5]}\n')
print(f'fused spliced loci {n}: no boundary (gene spans overlap) {nb}; cut statistic holds at {k}')
m = fc = 0
for g in sorted((x for x in gene if x not in fused and re.fullmatch(r'chr[0-9XY]+', gene[x][0])), key=lambda x: hashlib.md5(x.encode()).hexdigest()):
    c, s, T = gene[g]; intr = {(e[i][1], e[i + 1][0]) for e in map(sorted, T.values()) for i in range(len(e) - 1)}
    if not intr: continue
    a, b = span(g); R = reads(c, s, a, b)
    if len(R) < 20: continue
    m += 1; hit = [(i, j) for i, j in sorted(intr) if cut(sum(e <= i for st, e in R), sum(st >= j for st, e in R), sum(st < i and e > j for st, e in R))]
    fc += bool(hit); out.write(f'control\t{g}\t{len(R)}\t{len(intr)}\t{len(hit)}\n')
    if m == 1000: break
print(f'control genes {m}: false cut at some intron {fc} = {fc / max(1, m):.3f}')
