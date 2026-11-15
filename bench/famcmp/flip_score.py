#!/usr/bin/env python3
"""PREREG r1613 (docs/PREREG_secondary_as_primary_2026-11-10.md): family members modelled, loosely and exactly, per GTF.

  flip_score.py PANEL_GENES_TSV CAT_TSV REF_GTF NAME=GTF ...
"""
import sys, re, gzip, collections
def models(path):  # (chrom, strand) -> [exons] per transcript
    ex = collections.defaultdict(list)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
    by = collections.defaultdict(list)
    for (c, s, _), e in ex.items(): e.sort(); by[(c, s)].append(e)
    return by
chain = lambda e: tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))
with open(sys.argv[1]) as fh:
    h = fh.readline().rstrip('\n').split('\t'); P = [dict(zip(h, l.rstrip('\n').split('\t'))) for l in fh]
cls = lambda r: 'no-primary' if float(r['primary_molecules']) == 0 else 'ambiguous' if float(r['mapq0_share']) >= 0.5 else 'ordinary'
want = {r['gene'] for r in P}; G = {}
for l in open(sys.argv[2]):
    f = l.rstrip('\n').split('\t')
    if f[0] in want: G[f[0]] = (f[4], f[7], [tuple(map(int, b.split('-'))) for b in f[9].split(',')])
ref = {(c, s, chain(e)) for (c, s), ts in models(sys.argv[3]).items() for e in ts if len(e) > 1}
n = collections.Counter(cls(r) for r in P); print('classes', dict(n))
for a in sys.argv[4:]:
    name, path = a.split('=', 1); M = models(path); m1 = collections.Counter(); m2 = collections.Counter(); fam = set()
    for ts in M.values(): ts.sort(key=lambda e: e[0][0])
    for r in P:
        c, s, gx = G[r['gene']]; lo, hi = gx[0][0], gx[-1][1]; hit = exact = False
        for e in M.get((c, s), []):
            if e[0][0] >= hi: break
            if e[-1][1] <= lo or not any(x < b and a_ < y for a_, b in e for x, y in gx): continue
            hit = True
            if len(e) > 1:
                k = (c, s, chain(e)); fam.add(k); exact |= k in ref
        m1[cls(r)] += hit; m2[cls(r)] += exact
    tot = sum(len(v) for v in M.values())
    print(name, f'transcripts {tot}', *(f'{k}: M1 {m1[k]}/{n[k]} M2 {m2[k]}/{n[k]}' for k in ('no-primary', 'ambiguous', 'ordinary')),
          f'P family {len(fam & ref)}/{len(fam)} = {len(fam & ref) / max(len(fam), 1):.3f}', sep=' | ')
