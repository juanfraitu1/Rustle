#!/usr/bin/env python3
"""PREREG r1626 bars. usage: mono_bars.py REF_GTF MONO_READS_TSV DEFAULT_GTF FINAL_GTF NAME=GTF ...
Per GTF: single-exon models, E3, P, guard set S found, multi-exon chains equal to FINAL's, single-exon models >= 20 kb."""
import sys, re, gzip, collections
op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
def load(path):
    T = collections.defaultdict(list)
    for l in op(path):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': T[(f[0], f[6], re.search(r'transcript_id "([^"]*)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    for e in T.values(): e.sort()
    return T
def index(items):
    b = collections.defaultdict(list)
    for c, a, e, v in items:
        for k in range(a // 10000, e // 10000 + 1): b[(c, k)].append((a, e, v))
    return b
hits = lambda b, c, a, e: {x for k in range(a // 10000, e // 10000 + 1) for x in b.get((c, k), ()) if x[0] <= e and a <= x[1]}
rec = lambda a, b, x, y: min(b, y) - max(a, x) + 1 >= 0.5 * max(b - a + 1, y - x + 1)
ref = load(sys.argv[1]); refmono = index((c, e[0][0], e[0][1], s) for (c, s, t), e in ref.items() if len(e) == 1)
R = {}
for l in open(sys.argv[2]):
    c, a, b, s, t, same, anys, spl = l.rstrip('\n').split('\t'); R[t] = (c, int(a) + 1, int(b), s, int(same))
def view(p):
    T = load(p); mono = [(c, e[0][0], e[0][1], s) for (c, s, t), e in T.items() if len(e) == 1]
    chains = {(c, s, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))) for (c, s, t), e in T.items() if len(e) > 1}
    idx = index(mono)
    found = {t for t, (c, a, b, s, n) in R.items() if any(v == s and rec(a, b, x, y) for x, y, v in hits(idx, c, a, b))}
    match = sum(any(v == s and rec(a, b, x, y) for x, y, v in hits(refmono, c, a, b)) for c, a, b, s in mono)
    return mono, chains, found, match
dm, dc, dfound, _ = view(sys.argv[3]); S = {t for t in dfound if R[t][4] < 2}
fm, fc, ffound, fmatch = view(sys.argv[4]); E = {t for t, v in R.items() if v[4] >= 3}
print(f'expressed at 3: {len(E)}; guard set S: {len(S)}; bars: E3 >= {1.2 * len(ffound & E):.0f}, P >= {fmatch / max(1, len(fm)) - 0.05:.3f}, S >= {0.9 * len(S):.0f}')
print('arm\tsingle-exon models\tE3\tP\tS found\tmulti-exon chains = final\tmodels >= 20 kb\tverdict')
for arg in ['default=' + sys.argv[3], 'final=' + sys.argv[4]] + sys.argv[5:]:
    n, p = arg.split('=', 1); m, ch, found, match = view(p); P = match / max(1, len(m)); big = sum(b - a >= 20000 for c, a, b, s in m)
    ok = len(found & E) >= 1.2 * len(ffound & E) and P >= fmatch / max(1, len(fm)) - 0.05 and len(found & S) >= 0.9 * len(S) and ch == fc and big <= sum(b - a >= 20000 for c, a, b, s in fm)
    print(n, len(m), len(found & E), f'{P:.3f}', len(found & S), ch == fc, big, 'PASS' if ok else 'FAIL', sep='\t')
