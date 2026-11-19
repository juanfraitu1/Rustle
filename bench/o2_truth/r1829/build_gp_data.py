#!/usr/bin/env python3
"""r1829 page data for the gene / pseudogene section: gene_pseudo.json (gap bins per sample, pair and score) and
gene_pseudo_summary.json (per sample and pair: reads and shares), each {rows, meta} with this code.
  build_gp_data.py GENE_PSEUDO.tsv OUTDIR"""
import collections, json, os, sys

src, out = sys.argv[1:3]
CODE = open(os.path.abspath(__file__)).read()
GP_CODE = open(os.path.join(os.path.dirname(os.path.abspath(src)), 'gene_pseudo.py')).read()
meta = [{'label': 'build_gp_data.py', 'value': CODE, 'type': 'code', 'lang': 'python'},
        {'label': 'gene_pseudo.py', 'value': GP_CODE, 'type': 'code', 'lang': 'python'},
        {'label': 'Input', 'value': src + ' (gene_pseudo.py over the as_tables of A119b and HG002, CHM13 v2.0, CAT gene biotypes)'}]
C = collections.defaultdict(collections.Counter)
rows = []
for l in open(src):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'sample':
        continue
    C[(f[0], f[1])][(f[2], f[3])] = int(f[4])
    if f[2] in ('ms', 'AS', 'implied_divergence'):
        rows.append({'sample': f[0], 'pair': f[1], 'metric': f[2], 'bin': f[3], 'reads': int(f[4])})
json.dump({'rows': rows, 'meta': meta}, open(os.path.join(out, 'gene_pseudo.json'), 'w'))
summ = []
for (s, p), c in sorted(C.items()):
    n = c[('reads', '-')]
    near = lambda m: sum(c[(m, b)] for b in ('-10..-1', '0', '1..5', '6..10'))
    summ.append({'key': f'{s}|{p}', 'sample': s, 'pair': p, 'reads': n,
                 'other_is_primary': round(c[('other_is_primary', '-')] / n, 4),
                 'coding_spliced_other_not': round(c[('coding_spliced_other_not', '-')] / n, 4),
                 'AS_flips_to_other': round(c[('AS_prefers_other_while_ms_prefers_coding', '-')] / n, 4),
                 'ms_within_10': round(near('ms') / n, 4), 'AS_within_10': round(near('AS') / n, 4)})
json.dump({'rows': summ, 'meta': meta}, open(os.path.join(out, 'gene_pseudo_summary.json'), 'w'))
print(len(rows), len(summ))
