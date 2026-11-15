#!/usr/bin/env python3
"""Rows 1719 / 1723: why is the minor chain of a row 1714 pair (minor >= 5 reads, only the major emitted) missing?

  lost_why.py PAIRS_TSV DEFAULT_GTF [NAME=GTF ...]   class of each lost chain, and the share present in each further GTF"""
import sys, csv, collections, ast
from alt_events import load
pairs, gtf = sys.argv[1:3]; X = {a.split('=')[0]: load(a.split('=')[1]) for a in sys.argv[3:]}; csv.field_size_limit(1 << 30)
E = collections.defaultdict(dict)
for (c, s), d in load(gtf, 'reads').items():
    for ch, v in d.items(): E[c, s][ch] = (max(t[2] or 0 for t in v), ch[0][0], ch[-1][1])
T = collections.Counter(); seen = set()
for r in csv.DictReader(open(pairs), delimiter='\t'):
    a, b = int(r['reads_a']), int(r['reads_b'])
    if min(a, b) < 5 or r['default'] != 'major' or r['class'].startswith('NEST'): continue
    g = 'SITE' if r['class'].startswith('SITE') else r['class']
    m = ast.literal_eval(r['chain_b'] if b <= a else r['chain_a']); n = min(a, b); k = (r['contig'], r['strand'], m)
    if (g, k) in seen: continue
    seen.add((g, k)); em = E[r['contig'], r['strand']]; lo, hi = m[0][0], m[-1][1]
    ov = [(ch, v) for ch, v in em.items() if v[1] < hi and v[2] > lo]
    sub = [v[0] for ch, v in ov if len(ch) > len(m) and any(ch[i:i + len(m)] == m for i in range(len(ch) - len(m) + 1))]
    mx = max((v[0] for ch, v in ov), default=0)
    fz = any(len(ch) == len(m) and ch != m and all(abs(x[0] - y[0]) <= 10 and abs(x[1] - y[1]) <= 10 for x, y in zip(ch, m)) for ch, v in ov)
    why = 'sub-chain of an emitted chain, reads < 0.7 x it' if sub and n < 0.7 * max(sub) else 'sub-chain, reads >= 0.7 x' if sub else 'under 0.02 of the locus top chain' if n < 0.02 * mx else 'within 10 bp of an emitted chain' if fz else 'other'
    T[g, why] += 1; T[g] += 1
    for nm, x in X.items(): T[g, why, nm] += m in x.get((r['contig'], r['strand']), ())
for g in ('SITE', 'CASSETTE', 'AFE', 'ALE'):
    print(g, T[g], '; '.join(f'{w}: {v}' + ''.join(f' [{nm} {T[g, w, nm]}]' for nm in X) for (gg, w), v in sorted(((k, v) for k, v in T.items() if isinstance(k, tuple) and len(k) == 2), key=lambda x: -x[1]) if gg == g))
