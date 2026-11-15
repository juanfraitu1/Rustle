#!/usr/bin/env python3
"""PREREG r1638 (docs/PREREG_flair_profile_2026-11-11.md): score the support rules on subset_merge.py's all-chains table.

  flair_profile.py ALL_TSV N_DEFAULT N_DEFAULT_CORRECT N_FLAIR
ALL_TSV: contig, introns, full, merged unique, merged any, in default, refseq, stringtie, flair (read chains only; default
chains that are no read's chain are counted from N_DEFAULT and taken as outside FLAIR, a lower bound for the union arms)."""
import sys
rows = [tuple(map(int, l.split('\t')[2:])) for l in open(sys.argv[1])]
nd, ndc, nfl = map(int, sys.argv[2:5]); seen_d = sum(r[3] for r in rows); seen_dc = sum(r[3] and r[4] for r in rows)
print(f'read chains {len(rows)}; FLAIR chains that are a read chain {sum(r[6] for r in rows)} of {nfl}; default chains that are a read chain {seen_d} of {nd}')
print('arm\tchains\tcorrect\tprecision\tshared\treached\tinside\tjaccard')
def line(name, n, ok, sh): print(f'{name}\t{n}\t{ok}\t{ok / max(n, 1):.3f}\t{sh}\t{sh / nfl:.3f}\t{sh / max(n, 1):.3f}\t{sh / (n + nfl - sh):.3f}')
best = None
for m, col in (('U', 1), ('A', 2)):
    for S in range(2, 7):
        for union in (0, 1):
            keep = [r for r in rows if r[0] + r[col] >= S or (union and r[3])]
            n = len(keep) + union * (nd - seen_d); ok = sum(r[4] for r in keep) + union * (ndc - seen_dc); sh = sum(r[6] for r in keep)
            line(f'{m}{S}' + ('+default' if union else ''), n, ok, sh)
