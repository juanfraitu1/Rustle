#!/usr/bin/env python3
"""Diagnosis for seqscore (2026-10-30): why 98%-identity pairs are missed (gene in no family / placed apart) and what the predicted
pairs below tier 0.85 are.  seqdiag.py SET RUNDIR ARM.fam.copies.tsv  (run inside ~/runs/grid_gorilla; SET.seqtruth.pkl from seqscore)"""
import sys, os, csv, collections, itertools, pickle, subprocess, tempfile
sys.path.insert(0, os.path.expanduser('~/Rustle/bench/denovo_eval'))
import famtruth as F
p, d, arm = sys.argv[1:4]
G = F.load_ggo(F.load_cat()); idx = F.exon_index(G)
expressed = F.expressed_from_gtf(os.path.expanduser(f'~/runs/{d}/{p}.families.gtf'), idx)
best = pickle.load(open(f'{p}.seqtruth.pkl', 'rb'))
# all alignments incl. cov < 0.5 for the false-pair diagnosis: recompute quickly from cache? cache only has cov>=0.5; redo a light pass
genes = sorted(g for g in expressed if G[g]['exons'])
fams, _, _ = F.load_catalog(arm, G, idx)
gene_fam = collections.defaultdict(set)
for f, v in fams.items():
    for g in v: gene_fam[g].add(f)
pred = {(x, y) for v in fams.values() for x, y in itertools.combinations(sorted(v), 2) if x in expressed and y in expressed}
# loci -> gene, to know whether a gene has a de novo locus at all and whether it is in any family
loc_gene = {}
for l in open(os.path.expanduser(f'~/runs/{d}/{p}.fam.loci.gff3')):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if f[2] != 'gene': continue
has_locus = expressed  # expressed = has a transcript in the GTF
edges = {}
for l in open(arm.replace('.fam.copies.tsv', '.edges.tsv')):
    a_, b_, n, c, m = l.split('\t')[:5]; edges[(a_, b_)] = float(c)
truth98 = {k for k, i in best.items() if i >= 0.98}
miss = truth98 - pred
cat = collections.Counter(); ex = collections.defaultdict(list)
for x, y in miss:
    bx, by = G[x]['biotype'], G[y]['biotype']
    if not gene_fam[x] and not gene_fam[y]: k = 'neither gene in any family'
    elif not gene_fam[x] or not gene_fam[y]: k = 'one gene in no family'
    else: k = 'both placed, different families'
    cat[(k, 'coding-coding' if bx == by == 'protein_coding' else 'with pseudogene' if 'pseudo' in bx or 'pseudo' in by else 'other')] += 1
    if len(ex[k]) < 6: ex[k].append(f"{G[x]['name']}/{G[y]['name']} {G[x]['chrom'] == G[y]['chrom'] and 'same-chrom' or 'diff-chrom'}")
print(f"{p} {os.path.basename(arm)}: truth@0.98 {len(truth98)}, missed {len(miss)}")
for k, n in sorted(cat.items(), key=lambda x: -x[1]): print(f"   {n:4d}  {k[0]:35s} {k[1]}")
for k, v in ex.items(): print(f"   e.g. {k}: {v}")
# false pairs at 0.85: identity of best alignment at cov>=0.5 (in cache) or none
fp = [k for k in pred if best.get(k, 0) < 0.85]
c2 = collections.Counter('aligned cov>=0.5 but identity < 0.85' if k in best else 'no alignment covering half of the shorter' for k in fp)
print(f"   predicted {len(pred)}, below tier 0.85: {len(fp)} -> {dict(c2)}")
bt = collections.Counter()
for x, y in fp:
    if (x, y) not in best:
        bx, by = sorted([G[x]['biotype'], G[y]['biotype']]); bt[f"{bx}+{by}"] += 1
print('   no-alignment pairs by biotype:', dict(bt.most_common(6)))
