#!/usr/bin/env python3
"""Row 1945 (docs/STEP0_vg_mpmap_2026-11-25.md): one copy per read from a GAF on a family graph whose paths are named
FAMILY|COPY|KIND|ID. Best alignment per read (aligned share, then identity); >= 0.90 of the read aligned; compatible with
every path that holds all its nodes; assigned iff those paths are one copy index.
    graph_pick.py GRAPH.gfa ALN.gaf > picks.tsv   (read, copy | abstain reason; counts on stderr)"""
import sys, re, collections
gfa, gaf = sys.argv[1:3]
paths = {}
for l in open(gfa):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'P': paths[f[1]] = {x[:-1] for x in f[2].split(',')}
    elif f[0] == 'W': paths[f[3] if f[1] == '*' else f'{f[1]}#{f[2]}#{f[3]}'] = set(re.findall(r'[<>]([^<>]+)', f[6]))
best = {}
for l in open(gaf):
    f = l.rstrip('\n').split('\t')
    if len(f) < 12 or f[5] == '*': continue
    q = f[0].split()[0]; q = q[2:] if q.startswith('R|') else q
    cov = (int(f[3]) - int(f[2])) / int(f[1]); idn = next((float(t[5:]) for t in f[12:] if t.startswith('id:f:')), int(f[9]) / max(1, int(f[10])))
    if q not in best or (cov, idn) > best[q][:2]: best[q] = (cov, idn, set(re.findall(r'[<>]([^<>]+)', f[5])))
n = collections.Counter(aligned=len(best))
for q, (cov, _, nodes) in sorted(best.items()):
    if cov < 0.90: r = 'abstain_cover'
    else:
        ks = {p.split('|')[1] for p, ns in paths.items() if nodes <= ns and '|' in p}
        r = 'abstain_no_path' if not ks else ('abstain_several_copies' if len(ks) > 1 else ks.pop())
    n[r if r.startswith('abstain') else 'assigned'] += 1; print(q, r, sep='\t')
print(dict(n), file=sys.stderr)
