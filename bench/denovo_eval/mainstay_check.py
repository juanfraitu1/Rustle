#!/usr/bin/env python3
"""Mainstay / Y-ampliconic family check (human only): per named family, expressed members, families holding them, largest, foreign
protein-coding genes in those families.  mainstay_check.py BASE.families.gtf COPIES_TSV...   (env as famtruth)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
G = F.load_cat(); idx = F.exon_index(G)
expressed = F.expressed_from_gtf(sys.argv[1], idx)
_, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS['simiiformes'])
pats = dict(F.MAINSTAYS); pats.update(F.YAGS)
for p in sys.argv[2:]:
    fams, _, _ = F.load_catalog(p, G, idx)
    rep = F.named(G, expressed, fams, universe, pats)
    print('==', os.path.basename(p))
    for n, v in rep.items():
        if v['expressed'] == 0: continue
        print(f"{n:8s} exp {v['expressed']:2d}/{v['members']:2d} fams {v['families']} largest {v['largest']:2d} foreign {v['foreign']:2d} {v['foreign_names'][:6]}")
