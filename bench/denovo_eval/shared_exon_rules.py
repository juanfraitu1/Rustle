#!/usr/bin/env python3
"""Usage: shared_exon_rules.py LOCI.gff3 LOCI.gated.paf GENES.tsv FRAC  (GENES.tsv: contig start end strand name biotype, from the RefSeq GFF)

Per locus pair of the gated all-vs-all: the shared-exon evidence under the interval rule (min of exon bases inside the
record's query span and target span) and under the column rule (CIGAR columns joining exon to exon), each maxed over the
pair's records, as fractions of the smaller locus's exon length. Lists the pairs the column rule rejects at 0.60."""
import os, sys, re, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from family_container import aligned_runs, exon_columns, mask_ranges

gff, paf, genes_tsv, frac = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
blocks = collections.defaultdict(list)
for line in open(gff):
    if line.startswith('#'):
        continue
    r = line.rstrip('\n').split('\t')
    if r[2] != 'exon':
        continue
    g = re.search(r'gene=([^;]+)', r[8]).group(1)
    blocks[g].append((int(r[3]) - 1, int(r[4])))
# locus name in the PAF is CONTIG:START-END (1-based start); the gene id in the GFF is DN_CONTIG_START0_N -- map by span
span_of = {}
for line in open(gff):
    if line.startswith('#'):
        continue
    r = line.rstrip('\n').split('\t')
    if r[2] == 'gene':
        g = re.search(r'ID=gene-([^;]+)', r[8]).group(1)
        span_of[f'{r[0]}:{r[3]}-{r[4]}'] = g
def locus(name):
    g = span_of[name]
    c, s, e = re.match(r'(.+):(\d+)-(\d+)$', name).groups()
    b = sorted(blocks[g])
    m = []
    for s0, e0 in b:
        if m and s0 <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e0)
        else:
            m.append([s0, e0])
    st, en = [x[0] for x in m], [x[1] for x in m]
    return c, int(s) - 1, st, en, sum(e0 - s0 for s0, e0 in m)
cache = {}
def L(name):
    if name not in cache:
        cache[name] = locus(name)
    return cache[name]
def exonic_in(st, en, a, b):
    return sum(hi - lo for lo, hi, _ in mask_ranges(st, en, a, b))
best = {}
for line in open(paf):
    f = line.rstrip('\n').split('\t')
    q, t = f[0], f[5]
    if q == t:
        continue
    qc, qoff, qst, qen, qlen = L(q)
    tc, toff, tst, ten, tlen = L(t)
    if qc == tc and qoff <= int(f[8]) + toff and toff <= int(f[3]) + qoff and False:
        pass
    qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
    cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
    iv = min(exonic_in(qst, qen, qoff + qs, qoff + qe), exonic_in(tst, ten, toff + ts, toff + te))
    col = sum(qh - ql for ql, qh, _, _, _, _ in exon_columns(cg, f[4], qs, qe, ts, te, qoff, toff, (qst, qen), (tst, ten)))
    key = (q, t) if q < t else (t, q)
    b = best.setdefault(key, [0, 0, min(qlen, tlen), 0.0])
    b[0] = max(b[0], iv); b[1] = max(b[1], col)
    b[3] = max(b[3], int(f[9]) / max(1, int(f[10])))
genes = collections.defaultdict(list)
for l in open(genes_tsv):
    c, s, e, st, n, bt = l.rstrip('\n').split('\t'); genes[c].append((int(s), int(e), n, bt))
def names(name):
    c, s, e = re.match(r'(.+):(\d+)-(\d+)$', name).groups(); s, e = int(s), int(e)
    return ','.join(f'{n}' for gs, ge, n, bt in genes.get(c, ()) if gs <= e and s <= ge) or '-'
n_pairs = len(best)
both = sum(1 for v in best.values() if v[0] / max(1, v[2]) >= frac and v[1] / max(1, v[2]) >= frac)
iv_only = [(k, v) for k, v in best.items() if v[0] / max(1, v[2]) >= frac and v[1] / max(1, v[2]) < frac]
print(f'pairs with records: {n_pairs}; pass both rules at {frac}: {both}; pass interval only: {len(iv_only)}; '
      f'pass column only: {sum(1 for v in best.values() if v[0] / max(1, v[2]) < frac and v[1] / max(1, v[2]) >= frac)}')
hist = collections.Counter()
for k, v in iv_only:
    hist[min(int(v[1] / max(1, v[2]) * 10), 5)] += 1
print('column fraction of the interval-only pairs (deciles 0.0-0.5):', sorted(hist.items()))
print('pair\tinterval_frac\tcolumn_frac\tsmaller_exon_bp\tbest_identity\tgenes_q\tgenes_t')
for k, v in sorted(iv_only, key=lambda kv: -kv[1][3]):
    print(f'{k[0]} {k[1]}\t{v[0]/max(1,v[2]):.2f}\t{v[1]/max(1,v[2]):.2f}\t{v[2]}\t{v[3]:.3f}\t{names(k[0])}\t{names(k[1])}')
