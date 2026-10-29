#!/usr/bin/env python3
"""Judge sketchfam edges (tid pairs) with the famtruth-based Judge: verdict x shared/containment/masked-share bins.
  sketch_judge.py LOCI_GFF3 EDGES_TSV [--species gorilla]"""
import sys, os, re, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exon_rule_judge import Judge
gff, edges = sys.argv[1], sys.argv[2]
species = sys.argv[sys.argv.index('--species') + 1] if '--species' in sys.argv else 'gorilla'
key = {}
for line in open(gff):
    if line[0] != '#':
        f = line.rstrip('\n').split('\t')
        if f[2] == 'gene': key[re.search(r'ID=gene-([^;]+)', f[8]).group(1)] = f"{f[0]}:{f[3]}-{f[4]}"
J = Judge(species, gff)
rows = []
for line in open(edges):
    x, y, n, c, m = line.rstrip('\n').split('\t')
    rows.append((J.verdict(key[x], key[y]), int(n), float(c), float(m), key[x].split(':')[0] == key[y].split(':')[0]))
def b_n(n): return '5-9' if n < 10 else '10-29' if n < 30 else '30-99' if n < 100 else '>=100'
def b_c(c): return '<0.3' if c < .3 else '<0.5' if c < .5 else '<0.8' if c < .8 else '>=0.8'
def b_m(m): return '0' if m == 0 else '<0.5' if m < .5 else '<1' if m < 1 else '1'
for name, fn in (('shared', lambda r: b_n(r[1])), ('containment', lambda r: b_c(r[2])), ('masked share', lambda r: b_m(r[3])), ('same chrom', lambda r: str(r[4]))):
    t = collections.defaultdict(collections.Counter)
    for r in rows: t[r[0]][fn(r)] += 1
    print(f"--- {name}")
    for v in sorted(t): print(f"  {v:18s}", dict(sorted(t[v].items())))
