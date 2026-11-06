#!/usr/bin/env python3
"""Families with no copy overlapping a `gene` feature of a GFF3: one family id per line (row 1551).
usage: annotated_families.py COPIES_TSV GENES_GFF3 > unannotated.txt"""
import bisect, collections, sys

copies, gff = sys.argv[1], sys.argv[2]
genes = collections.defaultdict(list)
for line in open(gff):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5 and f[2] == 'gene':
        genes[f[0]].append((int(f[3]), int(f[4])))
starts, maxend = {}, {}
for c, v in genes.items():
    v.sort()
    starts[c] = [a for a, _ in v]
    m, cur = [], 0
    for _, b in v:
        cur = max(cur, b)
        m.append(cur)
    maxend[c] = m

def hit(c, a, b):
    if c not in genes:
        return False
    i = bisect.bisect_right(starts[c], b)
    return i > 0 and maxend[c][i - 1] >= a

fams, annotated = [], set()
with open(copies) as fh:
    header = fh.readline().rstrip('\n').split('\t')
    col = {k: i for i, k in enumerate(header)}
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) < len(header):
            continue
        fid = f[col['family_id']]
        if fid not in annotated and fid not in fams:
            fams.append(fid)
        if hit(f[col['chrom']], int(f[col['start']]), int(f[col['end']])):
            annotated.add(fid)
seen = set()
for fid in fams:
    if fid not in annotated and fid not in seen:
        seen.add(fid)
        print(fid)
print(f"[annotated_families] {len(set(fams))} families, {len(annotated)} annotated, {len(seen)} unannotated", file=sys.stderr)
