#!/usr/bin/env python3
"""exon_rules.py LOCI_GFF3 GATED_PAF GENES_TSV [T=0.5]  (2026-10-29 per-exon shared-exon measurement; GENES_TSV = chrom,start,end,strand,name,biotype)

Per locus pair of the gated all-vs-all: shared-exon evidence as (a) column fraction of the smaller locus's exon bases (the
--shared-exon-columns rule, best record), (b) EXON fraction: exons of the smaller locus with >= T of their bases joined to exon bases of
the other locus, over the union of the pair's records at identity >= 0.70, divided by the smaller locus's exon count; and (c) the same
with the union of column stretches. Prints one row per pair with gene names, to be joined against the missed true pairs and the admitted edges."""
import os, sys, re, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from family_container import exon_columns, mask_ranges
gff, paf, genes_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
T = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
blocks = collections.defaultdict(list); span_of = {}
for line in open(gff):
    if line.startswith('#'): continue
    r = line.rstrip('\n').split('\t')
    if r[2] == 'exon': blocks[re.search(r'gene=([^;]+)', r[8]).group(1)].append((int(r[3]) - 1, int(r[4])))
    elif r[2] == 'gene': span_of[f'{r[0]}:{r[3]}-{r[4]}'] = re.search(r'ID=gene-([^;]+)', r[8]).group(1)
cache = {}
def L(name):
    if name in cache: return cache[name]
    c, s, e = re.match(r'(.+):(\d+)-(\d+)$', name).groups()
    m = []
    for s0, e0 in sorted(blocks[span_of[name]]):
        if m and s0 <= m[-1][1]: m[-1][1] = max(m[-1][1], e0)
        else: m.append([s0, e0])
    cache[name] = (c, int(s) - 1, [x[0] for x in m], [x[1] for x in m], sum(e0 - s0 for s0, e0 in m))
    return cache[name]
# per pair: best-record column count, per-exon joined bases (union over records), on each side
acc = {}
for line in open(paf):
    f = line.rstrip('\n').split('\t'); q, t = f[0], f[5]
    if q == t: continue
    ident = int(f[9]) / max(1, int(f[10]))
    qc, qoff, qst, qen, qlen = L(q); tc, toff, tst, ten, tlen = L(t)
    qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
    cg = next(x[5:] for x in f[12:] if x.startswith('cg:Z:'))
    cols = list(exon_columns(cg, f[4], qs, qe, ts, te, qoff, toff, (qst, qen), (tst, ten)))
    key = (q, t) if q < t else (t, q); flip = q > t
    a = acc.setdefault(key, dict(best=0, bestid=0.0, maxid=0.0, qiv=collections.defaultdict(list), tiv=collections.defaultdict(list)))
    n = sum(qh - ql for ql, qh, *_ in cols)
    a['maxid'] = max(a['maxid'], ident)
    if n > a['best']: a['best'], a['bestid'] = n, ident
    if ident >= 0.70:
        for ql, qh, qb, tl, th, tb in cols:
            (a['tiv'] if flip else a['qiv'])[qb].append((ql, qh))
            (a['qiv'] if flip else a['tiv'])[tb].append((tl, th))
def ulen(iv):
    n, cur = 0, None
    for s, e in sorted(iv):
        if cur is None or s > cur[1]:
            if cur: n += cur[1] - cur[0]
            cur = [s, e]
        else: cur[1] = max(cur[1], e)
    if cur: n += cur[1] - cur[0]
    return n
genes = collections.defaultdict(list)
for l in open(genes_tsv):
    c, s, e, st, n, bt = l.rstrip('\n').split('\t'); genes[c].append((int(s), int(e), n, bt))
def names(name):
    c, s, e = re.match(r'(.+):(\d+)-(\d+)$', name).groups(); s, e = int(s), int(e)
    return ','.join(n for gs, ge, n, bt in genes.get(c, ()) if gs <= e and s <= ge and bt == 'protein_coding') or \
           ','.join(n for gs, ge, n, bt in genes.get(c, ()) if gs <= e and s <= ge) or '-'
print('a\tb\tsmaller_bp\tn_exons_smaller\tbest_col_frac\tbest_id\tmax_id\tunion_col_frac\texon_frac\texons_joined\tnames_a\tnames_b')
for (a_, b_), a in acc.items():
    A, B = L(a_), L(b_)
    small_is_a = A[4] <= B[4]
    S, iv = (A, a['qiv']) if small_is_a else (B, a['tiv'])
    nex = len(S[2]); joined = 0; ucols = 0
    for bi in range(nex):
        ex_len = S[3][bi] - S[2][bi]
        j = ulen(iv.get(bi, []))
        ucols += j
        if ex_len > 0 and j / ex_len >= T: joined += 1
    print(f"{a_}\t{b_}\t{S[4]}\t{nex}\t{a['best']/max(1,S[4]):.3f}\t{a['bestid']:.3f}\t{a['maxid']:.3f}\t{ucols/max(1,S[4]):.3f}\t{joined/max(1,nex):.3f}\t{joined}\t{names(a_)}\t{names(b_)}")
