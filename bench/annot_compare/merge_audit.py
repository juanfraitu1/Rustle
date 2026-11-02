#!/usr/bin/env python3
"""Merge audit (2026-11-02, descriptive): which loci (gene_id) of an assembled GTF span >= 2 coding RefSeq genes at their exons, and why:
  shared-exon cluster = those RefSeq genes share exonic bases with each other (PCDH/UGT1A-like: one RNA unit, separate genes)
  fused               = they do not (two separate genes joined: readthrough or a merge artefact)
usage: merge_audit.py ASSEMBLY.gtf[.gz] REFSEQ.gtf.gz OUT.tsv  -> summary on stdout"""
import collections, gzip, re, sys
sys.path.insert(0, __import__('os').path.dirname(__file__))
from cat_vs_refseq_lib import Idx, merge, share  # noqa
asm, rs, out = sys.argv[1:4]
op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
A = lambda k, a: (re.search(k + r' "([^"]+)"', a) or [None, ''])[1]


def genes(p, coding_only=False):
    ex, meta, tx = collections.defaultdict(list), {}, collections.defaultdict(set)
    for l in op(p):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        g = A('gene_id', f[8])
        ex[g].append((int(f[3]) - 1, int(f[4])))
        meta[g] = (f[0], f[6], A('gene_name', f[8]) or g)
        tx[g].add(A('transcript_id', f[8]))
    keep = {g for g in ex if not coding_only or any(t.startswith(('rna-NM_', 'rna-XM_')) for t in tx[g])}
    return {g: (meta[g][0], meta[g][1], merge(ex[g])) for g in keep}, meta, tx


R, rmeta, _ = genes(rs, True)
G, _, gtx = genes(asm)
ri = Idx(R)
n = collections.Counter()
rows = []
for g, (c, s, ex) in G.items():
    h = ri.hits(c, s, ex)
    if len(h) < 2:
        continue
    pairs = [(a, b) for i, a in enumerate(h) for b in h[i + 1:]]
    kind = 'shared-exon cluster' if all(share(R[a][2], R[b][2]) for a, b in pairs) else (
        'fused' if not any(share(R[a][2], R[b][2]) for a, b in pairs) else 'mixed')
    n[kind] += 1
    rows.append((g, c, ex[0][0], ex[-1][1], s, len(h), kind, len(gtx[g]), ','.join(sorted(rmeta[x][2] for x in h))))
with open(out, 'w') as fh:
    fh.write('locus\tchrom\tstart\tend\tstrand\trefseq_coding_genes\tkind\ttranscripts\tgenes\n')
    for r in sorted(rows, key=lambda r: -r[5]):
        fh.write('\t'.join(map(str, r)) + '\n')
print(f'loci {len(G)}; spanning >= 2 coding RefSeq genes: {len(rows)} -> ' + ', '.join(f'{k} {v}' for k, v in n.most_common()))
for r in sorted(rows, key=lambda r: -r[5])[:12]:
    print(f'  {r[0]} {r[1]}:{r[2]}-{r[3]} {r[6]} {r[5]} genes: {r[8][:120]}')
