#!/usr/bin/env python3
"""r1528: catalog coverage of the family panel's expressed members (T >= 2 placements): how many overlap a copy of the given
copy table on the gene's strand, split by whether base assembled a model there.  usage: catalog_cov.py GENES.tsv NAME=COPIES.tsv ...
(GENES.tsv = family_panel.py's *.genes.tsv with a 'base' column)"""
import csv, sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F
G = F.load_cat()
mem = [r for r in csv.DictReader(open(sys.argv[1]), delimiter='\t') if int(r['tied_molecules']) >= 2]
for arg in sys.argv[2:]:
    name, path = arg.split('=', 1)
    copies = collections.defaultdict(list); fam_n = collections.Counter()
    for r in csv.DictReader(open(path), delimiter='\t'):
        copies[(r['chrom'], r['strand'])].append((int(r['locus_start']), int(r['locus_end']), r['family_id'])); fam_n[r['family_id']] += 1
    tab = collections.Counter()
    for m in mem:
        v = G[m['gene']]; lo, hi = v['exons'][0][0], v['exons'][-1][1]
        cov = any(a < hi and b > lo for a, b, _ in copies[(v['chrom'], v['strand'])])
        tab[('covered' if cov else 'uncovered', 'hit' if m['base'] != '0' else 'miss')] += 1
    c = tab[('covered', 'hit')] + tab[('covered', 'miss')]
    print(f"[cov] {name}: expressed members {len(mem)}, covered by a catalog copy {c} (hit {tab[('covered','hit')]}, miss {tab[('covered','miss')]}); "
          f"uncovered hit {tab[('uncovered','hit')]} miss {tab[('uncovered','miss')]}; families {len(fam_n)}, copies {sum(fam_n.values())}, largest {max(fam_n.values())}")
