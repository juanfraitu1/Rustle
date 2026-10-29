#!/usr/bin/env python3
"""Family-level overlap merge (2026-10-29, mainstay diagnosis: TBC1D3 and NPIPB4 carry several overlapping de novo loci per copy, each
locus type forming its own family). Two families are merged when the member loci overlap pairwise on the same strand in >= --frac of
the smaller family's members (and >= 2 members); transitive. Loci stay as they are (no sketch change).
  famoverlap.py IN.fam.copies.tsv OUT.fam.copies.tsv [--frac 0.5]"""
import sys, csv, collections, bisect
inp, out = sys.argv[1:3]; FRAC = float(sys.argv[sys.argv.index('--frac') + 1]) if '--frac' in sys.argv else 0.5
rows = [r for r in csv.DictReader(open(inp), delimiter='\t')]
fams = collections.defaultdict(list)
for r in rows: fams[r['family_id']].append(r)
by = collections.defaultdict(list)
for r in rows: by[(r['chrom'], r['strand'])].append((int(r['start']), int(r['end']), r['family_id'], r['tid']))
for v in by.values(): v.sort()
pair_hits = collections.defaultdict(set)  # (famA, famB) -> set of (locusA) overlapping some locus of B
for k, v in by.items():
    starts = [x[0] for x in v]; maxlen = max(e - s for s, e, _, _ in v)
    for s, e, f, t in v:
        i = bisect.bisect_left(starts, s - maxlen)
        while i < len(v) and v[i][0] <= e:
            s2, e2, f2, t2 = v[i]; i += 1
            if f2 != f and s2 <= e and s <= e2:
                pair_hits[(f, f2)].add(t)
parent = {f: f for f in fams}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
merged = 0
for (a, b), ta in pair_hits.items():
    if a < b:
        tb = pair_hits.get((b, a), set())
        need = max(2, FRAC * min(len(fams[a]), len(fams[b])))
        if len(ta) >= need and len(tb) >= need:
            if find(a) != find(b): parent[find(a)] = find(b); merged += 1
with open(out, 'w') as w:
    w.write('\t'.join(rows[0].keys()) + '\n')
    groups = collections.defaultdict(list)
    for f, v in fams.items(): groups[find(f)] += v
    for i, (g, v) in enumerate(sorted(groups.items(), key=lambda x: -len(x[1]))):
        for j, r in enumerate(v):
            r = dict(r); r['family_id'] = f"FO{i}"; r['copy_idx'] = str(j); w.write('\t'.join(r.values()) + '\n')
print(f"families {len(fams)} -> {len(set(find(f) for f in fams))} ({merged} merges)", file=sys.stderr)
