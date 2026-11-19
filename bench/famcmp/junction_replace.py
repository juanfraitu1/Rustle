#!/usr/bin/env python3
"""Row 1819 (docs/PREREG_junction_replace_2026-11-19.md): junction_replace.py NAME TRUTH.gtf GTF [CONTIG_PREFIX]
Analysis only, writes no GTF: replaces an intron of a chain of GTF by the best-supported intron within 10 bp at both ends when
that one has at least R times its support (support = summed `reads` of the transcripts carrying the intron), for R in 2, 5, 10,
20, 50, and prints one tab-separated row per R: NAME, R, chains, correct, precision, chains after, correct after, precision
after, correct ratio, share of wrong chains fewer, wrong->correct, of them gained, wrong->wrong, correct->wrong, correct->correct."""
import sys, bisect, collections
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from curve import load
name, truth, gtf = sys.argv[1:4]; pre = sys.argv[4] if len(sys.argv) > 4 else ''
sel = lambda d: {k: v for k, v in d.items() if k[0].startswith(pre)}
T = set(sel(load(truth))); A = sel(load(gtf, 'reads'))
sup = collections.Counter()
for (c, s, ch), v in A.items():
    for i in ch: sup[c, s, i] += v[0]
by = collections.defaultdict(list)
for (c, s, i) in sup: by[c, s].append(i)
for v in by.values(): v.sort()
def best(c, s, i, R):
    v = by[c, s]; b = None
    for j in v[bisect.bisect_left(v, (i[0] - 10, 0)):bisect.bisect_right(v, (i[0] + 10, 1 << 60))]:
        if j != i and abs(j[1] - i[1]) <= 10 and sup[c, s, j] >= R * sup[c, s, i] and (b is None or sup[c, s, j] > sup[c, s, b]): b = j
    return b
def ok(ch): return all(a < b for a, b in ch) and all(ch[k][1] < ch[k + 1][0] for k in range(len(ch) - 1))
assert ok(((1, 5), (8, 9))) and not ok(((1, 5), (4, 9)))
c0 = sum(k in T for k in A)
for R in (2, 5, 10, 20, 50):
    new = set(); n = collections.Counter(); gained = 0
    for k in A:
        c, s, ch = k; ch2 = tuple(best(c, s, i, R) or i for i in ch)
        if ch2 == ch or not ok(ch2): new.add(k); continue
        k2 = (c, s, ch2); new.add(k2); n[('correct' if k in T else 'wrong') + ' -> ' + ('correct' if k2 in T else 'wrong')] += 1
        gained += k not in T and k2 in T and k2 not in A
    c1 = sum(k in T for k in new)
    print(name, R, len(A), c0, f'{c0 / len(A):.4f}', len(new), c1, f'{c1 / len(new):.4f}', f'{c1 / c0:.4f}', f'{1 - (len(new) - c1) / (len(A) - c0):.4f}',
          n['wrong -> correct'], gained, n['wrong -> wrong'], n['correct -> wrong'], n['correct -> correct'], sep='\t')
