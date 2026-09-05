#!/usr/bin/env python3
"""Adapt an `mcl_families` catalog to the old `cat.copies.tsv` layout so `bench/soto_adjudicate.py` scores
both O1 definitions with one scorer.  usage: mcl_to_cat_copies.py <prefix> <paf> <out_dir_units> <out_dir_all>
  <out_dir_units>/cat.copies.tsv — one row per emitted UNIT (read-supported members; "RNA disposes")
  <out_dir_all>/cat.copies.tsv   — one row per cluster MEMBER (annotation only; "DNA proposes")
`max_family_identity` = the member's best PAF identity (matches / block length) to a cluster mate."""
import sys, os, csv, collections
prefix, paf, out_u, out_a = sys.argv[1:5]
rows = [r for r in csv.DictReader(open(prefix + '.clusters.tsv'), delimiter='\t')]
key = lambda r: f"{r['chrom']}:{r['start']}-{r['end']}"
fam_of = {key(r): r['cluster_id'] for r in rows}
best = collections.defaultdict(float)
for l in open(paf):
    f = l.split('\t')
    if f[0] == f[5] or f[0] not in fam_of or fam_of.get(f[5]) != fam_of[f[0]]: continue
    ident = int(f[9]) / max(1, int(f[10]))
    best[f[0]] = max(best[f[0]], ident); best[f[5]] = max(best[f[5]], ident)
hdr = 'family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\n'
os.makedirs(out_a, exist_ok=True); os.makedirs(out_u, exist_ok=True)
idx = collections.Counter()
with open(f'{out_a}/cat.copies.tsv', 'w') as o:
    o.write(hdr)
    for r in rows:
        k = key(r); i = idx[r['cluster_id']]; idx[r['cluster_id']] += 1
        o.write(f"{r['cluster_id']}\t{i}\tANN_{r['chrom']}_{r['start']}\t{r['chrom']}\t{int(r['start'])-1}\t{r['end']}\tNA\t.\tNA\tNA\t{best.get(k,0):.6f}\n")
n = 0
with open(f'{out_u}/cat.copies.tsv', 'w') as o:
    o.write(hdr)
    for r in csv.DictReader(open(prefix + '.units.tsv'), delimiter='\t'):
        o.write('\t'.join([r['family_id'], r['copy_idx'], r['tid'], r['chrom'], r['start'], r['end'], r['n_exon'], r['strand'], r['n_reads'], r['exons'], r['nearest_ident'] if r['nearest_ident'] != 'NA' else '0']) + '\n'); n += 1
print(f'members {len(rows)} -> {out_a}; units {n} -> {out_u}')
