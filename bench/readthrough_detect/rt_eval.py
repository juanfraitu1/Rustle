#!/usr/bin/env python3
"""register 1492: structural readthrough calls vs the annotation's own labels (gene level, protein_coding / lncRNA genes).
usage: rt_eval.py ANNOT.gff[.gz] [CONTIGS]  (CONTIGS comma list; default all) -> recall, precision, misses / extras on stdout"""
import os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import units_from_annotation as U

gff = sys.argv[1]
contigs = set(sys.argv[2].split(',')) if len(sys.argv) > 2 and sys.argv[2] else None
bt = {}
for line in U.open_text(gff):
    f = line.split('\t')
    if len(f) < 9 or f[2] not in U.GENE_TYPES or (contigs and f[0] not in contigs):
        continue
    m, n = U._BIOTYPE.search(f[8]), U._ATTR['Name'].search(f[8])
    i = U._ATTR['ID'].search(f[8]); bt[f'{n.group(1) if n else (i.group(1) if i else "")}|{f[3]}'] = m.group(1).lower() if m else ''
keep = lambda k: bt.get(k) in ('protein_coding', 'lncrna')
lab, st = {}, {}
U.read_genes(gff, contigs, 'labels', None, lab)
U.read_genes(gff, contigs, 'structural', None, st)
L, S = {k for k in lab if keep(k)}, {k for k in st if keep(k)}
tp = len(L & S)
print(f'{os.path.basename(gff)}\tlabelled {len(L)}\tstructural {len(S)}\tboth {tp}\trecall {tp / max(len(L), 1):.3f}\tprecision {tp / max(len(S), 1):.3f}')
for tag, xs in (('MISS', sorted(L - S)), ('EXTRA', sorted(S - L))):
    print(f'{tag}\t{len(xs)}\t' + ' '.join(xs[:40]))
