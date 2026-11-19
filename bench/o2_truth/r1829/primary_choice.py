#!/usr/bin/env python3
"""r1829 side table (descriptive, not pre-registered; user 2026-11-19): when a multimapper's primary and its best rival
place (as_table's best-AS non-primary record, at another place) are within a gap, which one did minimap2 make primary?
Per sample and per score (AS, ms): gap bands of |primary - rival|; within each band, primary higher / equal / rival higher.
  primary_choice.py SAMPLE AS_TABLE > SAMPLE.choice.tsv"""
import collections, sys

BANDS = [(0, 0), (1, 5), (6, 10), (11, 25), (26, 50), (51, 100), (101, 10 ** 9)]
def band(g):
    return next(f"{a}" if a == b else f">{a - 1}" if b == 10 ** 9 else f"{a}-{b}" for a, b in BANDS if a <= g <= b)

s, t = sys.argv[1:3]
C = collections.Counter()
for l in open(t):
    if l[0] == '#':
        continue
    f = l.rstrip('\n').split('\t')
    if f[4] in ('NA', '*'):
        continue
    C['molecules'] += 1
    if f[11] == 'NA':
        continue
    if f[14] == f[4] and abs(int(f[16]) - int(f[15])) < 5000:
        C['same_place'] += 1; continue
    C['multimappers'] += 1
    # whether the primary carries the read's best score over ALL its records
    C['primary_is_best_AS'] += f[5] == f[1]
    C['primary_is_best_ms'] += f[8] == f[6]
    for m, p, r in (('AS', int(f[5]), int(f[11])), ('ms', int(f[8]), int(f[12]))):
        g = abs(p - r); w = 'primary higher' if p > r else 'equal' if p == r else 'rival higher'
        C[(m, band(g), w)] += 1
print('sample\tmetric\tband\twho\tcount')
for k in ('molecules', 'same_place', 'multimappers', 'primary_is_best_AS', 'primary_is_best_ms'):
    print(f"{s}\t-\t-\t{k}\t{C[k]}")
for k, v in sorted((k, v) for k, v in C.items() if isinstance(k, tuple)):
    print(f"{s}\t{k[0]}\t{k[1]}\t{k[2]}\t{v}")
