#!/usr/bin/env python3
"""mode_table.py REF NAME=GTF[,KEY | ,@TSV] ...   (register 1641) de novo block: multi-exon exact intron chains vs REF.
KEY = transcript attribute to cut on, higher is better (reads, cov). @TSV = bambu mcols table (TXNAME, NDR), lower is better;
a TSV whose header starts with feature_id is an IsoQuant counts table (higher is better)."""
import sys, re, collections
from curve import load
CUT = {'hi': (2, 3, 4, 6, 11, 21), 'ndr': (0.9, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05)}
ref = set(load(sys.argv[1])); sets = {}
rep = lambda n, s: print(f'{n}\t{len(s)}\t{len(s & ref)}\t{len(s & ref) / max(1, len(s)):.3f}')
print('set\tchains\tcorrect\tprecision')
for a in sys.argv[2:]:
    name, spec = a.split('=', 1); path, _, key = spec.partition(',')
    if key.startswith('@'):  # bambu: NDR by transcript id; a chain keeps its best (lowest) NDR
        rows = [l.rstrip('\n').split('\t') for l in open(key[1:]) if not l.startswith('__')]; h = rows[0]; cnt = h[0] == 'feature_id'  # IsoQuant counts table: higher is better, stored negated
        ndr = {r[0]: -float(r[1]) for r in rows[1:]} if cnt else {r[h.index('TXNAME')]: float(r[h.index('NDR')]) for r in rows[1:] if r[h.index('NDR')] != 'NA'}
        ex = collections.defaultdict(list)
        for l in open(path):
            f = l.split('\t')
            if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
        d = {}
        for (c, s, t), e in ex.items():
            e.sort()
            if len(e) > 1:
                ch = (c, s, tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))); d[ch] = min(d.get(ch, 9), ndr.get(t, 0.0))
        cuts = [(f'{name} count >= {t}', {c for c, v in d.items() if -v >= t}) for t in CUT['hi']] if cnt else [(f'{name} NDR <= {t}', {c for c, v in d.items() if v <= t}) for t in CUT['ndr']]
    else:
        d = load(path, key or None); cuts = [(f'{name} {key} >= {t}', {c for c, v in d.items() if v[0] >= t}) for t in CUT['hi']] if key else []
    sets[name] = set(d); rep(name, set(d))
    if any(c[1] not in '+-' for c in d): print(f'# {name}: {sum(c[1] not in "+-" for c in d)} multi-exon chains without a strand (never match)')
    for n, s in cuts: rep(n, s)
names = list(sets)
for i, x in enumerate(names):
    for y in names[i + 1:]: rep(f'{x} & {y}', sets[x] & sets[y])
for x in names: rep(f'{x} only', sets[x] - set().union(*(sets[y] for y in names if y != x)))
