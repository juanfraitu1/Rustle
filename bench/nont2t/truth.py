#!/usr/bin/env python3
"""PREREG r1504 truth: copies.fa aligned to the non-T2T genome (PAF with cg/NM via minimap2 -c). A copy is PRESENT when a hit has
identity >= 0.99 over >= 0.95 of the copy and its locus (target, strand, span) is not claimed by a better-identity copy; else
MISSING. EXPRESSED = >= 3 primary reads in the T2T BAM over its exons (counts given as gene<TAB>reads).
  truth.py COPIES.paf GENES.tsv EXPR.tsv > TRUTH.tsv   (gene name family status expressed best_identity best_cov locus reason)"""
import collections
import sys

paf, genes, expr = sys.argv[1:4]
meta = {}
for l in open(genes):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'gene':
        meta[f[0]] = f[1:3]
reads = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(expr)}
hits = collections.defaultdict(list)
for l in open(paf):
    f = l.rstrip('\n').split('\t')
    q, ql, qs, qe, st, t, ts, te, m, blk = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
    hits[q].append((m / blk, (qe - qs) / ql, t, st, ts, te))
cands = sorted(((h[0], q, h) for q, hs in hits.items() for h in hs if h[0] >= 0.99 and h[1] >= 0.95), reverse=True)
claimed, status = [], {}
for ident, q, (i, c, t, st, ts, te) in cands:
    if q in status:
        continue
    if any(t == t2 and st == s2 and ts < e2 and s1 < te for t2, s2, s1, e2 in claimed):
        continue
    claimed.append((t, st, ts, te))
    status[q] = ('present', i, c, f'{t}:{ts}-{te}{st}')
# why a copy is missing (report only): collapsed = it has a >= 99% / 95% hit, but on a locus a better copy claimed (the older
# assembly merged them); absent = no such hit anywhere
good = {q for _, q, _ in cands}
print('gene\tname\tfamily\tstatus\texpressed\tbest_identity\tbest_cov\tlocus\treason')
for g, (name, fam) in meta.items():
    best = max(hits.get(g, [(0, 0, '-', '.', 0, 0)]))
    s = status.get(g, ('missing', best[0], best[1], '-'))
    why = '-' if s[0] == 'present' else ('collapsed' if g in good else 'absent')
    print(f"{g}\t{name}\t{fam}\t{s[0]}\t{int(reads.get(g, 0) >= 3)}\t{s[1]:.4f}\t{s[2]:.3f}\t{s[3]}\t{why}")
