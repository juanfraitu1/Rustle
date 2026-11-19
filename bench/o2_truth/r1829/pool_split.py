#!/usr/bin/env python3
"""r1829 side table (descriptive; user 2026-11-19): the seeding pool (rival AS >= 0.98 x best AS) split by the ms gap between
the primary and the rival place, read as the number of sequence differences (PSVs) the read covers (about 5 ms each): gap 0 =
no PSV (exact tie), 1-10 = 1-2, 11-25 = about 3-5, 26-50 = about 6-10, > 50 = more (or a longer aligned part at one place).
Also the near ties (ms gap <= 10) that are NOT in the pool (AS also charges introns, which ms does not).
Same place rule as near_tie_census.py. Shares are of all molecules with a primary.
  pool_split.py SAMPLE AS_TABLE > SAMPLE.pool.tsv"""
import collections, sys

s, t = sys.argv[1:3]
C = collections.Counter()
for l in open(t):
    if l[0] == '#':
        continue
    f = l.rstrip('\n').split('\t')
    if f[4] in ('NA', '*'):
        continue
    C['molecules'] += 1
    if f[11] == 'NA' or (f[14] == f[4] and abs(int(f[16]) - int(f[15])) < 5000):
        continue
    g = abs(int(f[8]) - int(f[12]))
    pool = int(f[11]) >= 0.98 * int(f[1])
    if pool:
        C['pool'] += 1
        C['pool_gap_' + ('0' if g == 0 else '1-10' if g <= 10 else '11-25' if g <= 25 else '26-50' if g <= 50 else 'over50')] += 1
    elif g <= 10:
        C['near_not_pool'] += 1
print('sample\tkey\tcount\tshare')
for k in sorted(C):
    print(f"{s}\t{k}\t{C[k]}\t{C[k] / C['molecules']:.6f}")
