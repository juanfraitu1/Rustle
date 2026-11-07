#!/usr/bin/env python3
"""Families with no copy overlapping a `gene` feature of a GFF3: one family id per line (row 1551).
usage: annotated_families.py COPIES_TSV GENES_GFF3 > unannotated.txt
With --majority-in (row 1559): families whose copies MOSTLY overlap an element of a BED (0-based, cols 1-3),
e.g. T2T censat satellites: usage: annotated_families.py COPIES_TSV ELEMENTS_BED --majority-in > satellite.txt"""
import bisect, collections, sys

copies, gff = sys.argv[1], sys.argv[2]
majority = '--majority-in' in sys.argv[3:]
genes = collections.defaultdict(list)
for line in open(gff):
    if line.startswith('#') or line.startswith('track'):
        continue
    f = line.rstrip('\n').split('\t')
    # a BED is told by its columns (2-3 integers, not a GFF3 type in column 3), not by the file suffix
    if len(f) >= 3 and f[1].isdigit() and f[2].isdigit():
        genes[f[0]].append((int(f[1]) + 1, int(f[2])))
    elif len(f) >= 5 and f[2] == 'gene':
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

def die(msg):
    print(f"[annotated_families] ERROR: {msg}", file=sys.stderr)
    sys.exit(2)

def hit(c, a, b):
    if c not in genes:
        return False
    i = bisect.bisect_right(starts[c], b)
    return i > 0 and maxend[c][i - 1] >= a

fams, annotated, chroms = {}, set(), set()
n_in, n_all = collections.Counter(), collections.Counter()
with open(copies) as fh:
    header = fh.readline().rstrip('\n').split('\t')
    col = {k: i for i, k in enumerate(header)}
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) < len(header):
            continue
        fid = f[col['family_id']]
        fams.setdefault(fid, None)
        chroms.add(f[col['chrom']])
        # copies.tsv is 0-based half-open: its first base is start + 1 in the 1-based closed intervals above
        h = hit(f[col['chrom']], int(f[col['start']]) + 1, int(f[col['end']]))
        n_all[fid] += 1
        n_in[fid] += h
        if h:
            annotated.add(fid)
# a different contig naming, or an annotation of another genome, would list EVERY family as unannotated
if chroms and not chroms & genes.keys():
    die(f"no copy contig ({sorted(chroms)[:3]}) is in {gff} ({sorted(genes)[:3]}): wrong genome or naming")
if fams and not annotated:  # a small panel can hold no satellite / no gene: say so, do not stop
    print(f"[annotated_families] WARNING: none of {len(fams)} families overlaps an element of {gff}", file=sys.stderr)
seen = set()
for fid in fams:
    out = 2 * n_in[fid] > n_all[fid] if majority else fid not in annotated
    if out and fid not in seen:
        seen.add(fid)
        print(fid)
print(f"[annotated_families] {len(fams)} families, {len(annotated)} with a hit, {len(seen)} listed ({'majority-in' if majority else 'no hit'})", file=sys.stderr)
