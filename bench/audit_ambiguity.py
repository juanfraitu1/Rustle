#!/usr/bin/env python3
"""Audit whether O2's `tied` / `ambiguous` molecules are ABSOLUTELY ambiguous, per CANDIDATE.

From `--dump-star`: per molecule, every column where the family's copies differ, the read's base and each
candidate's base (`.` = not covered). For each molecule and each candidate: mismatches (covered columns where
the candidate's base != the read's). The CONSISTENT set S = candidates with 0 mismatches.

  |S| >= 2  -> absolutely ambiguous among S: every member matches the read at every column the read
               covers, so no substitution PSV in the read's span separates them (a `tied` verdict is right)
  |S| == 1  -> a unique candidate explains every column: NOT ambiguous by PSVs (a `tied` verdict is wrong,
               an `ambiguous` verdict means the margin/certificate, not the evidence, blocked it)
  |S| == 0  -> the read matches no candidate perfectly: best-vs-second mismatch counts say how close

usage: audit_ambiguity.py <star_reads.tsv> <assignments.tsv>
"""
import sys, csv, collections
star_p, asg_p = sys.argv[1:3]
A = {r['read_name']: r for r in csv.DictReader(open(asg_p), delimiter='\t')}
by = collections.defaultdict(list)
for r in csv.DictReader(open(star_p), delimiter='\t'):
    a = A.get(r['read_name'])
    if not a or a['status'] not in ('tied', 'ambiguous') or int(r['n_candidates']) < 2:
        continue
    cands = r['candidates'].split(',')
    mm = [0] * len(cands); cov = [0] * len(cands)
    for col in r['columns'].split(','):
        if not col:
            continue
        pos, rb, cb = col.split(':', 2)
        for k, b in enumerate(cb[:len(cands)]):
            if b == '.':
                continue
            cov[k] += 1
            if b != rb:
                mm[k] += 1
    S = [k for k in range(len(cands)) if cov[k] > 0 and mm[k] == 0]
    covered = [k for k in range(len(cands)) if cov[k] > 0]
    srt = sorted(mm[k] for k in covered) if covered else []
    best_mm = srt[0] if srt else None
    second_mm = srt[1] if len(srt) > 1 else None
    key = (a['status'],
           'origin_rejected' if a['origin_rejected'] == '1' else 'explained',
           'tie_outside' if a.get('tie_outside_catalog') == '1' else 'inside')
    by[key].append((len(S), best_mm, second_mm, len(covered), int(r['n_cols'])))

for key in sorted(by):
    v = by[key]; n = len(v)
    s2 = sum(1 for x in v if x[0] >= 2); s1 = sum(1 for x in v if x[0] == 1); s0 = sum(1 for x in v if x[0] == 0)
    print(f"{' / '.join(key):45s} n={n:5d}   |S|>=2 (absolutely ambiguous) {s2:5d} ({s2/n:5.1%})   "
          f"|S|==1 (a unique perfect candidate) {s1:5d} ({s1/n:5.1%})   |S|==0 {s0:5d}")
    if s1:
        g = sorted(x[2] - x[1] for x in v if x[0] == 1 and x[2] is not None)
        if g: print(f"{'':45s}   when |S|==1: runner-up's extra mismatches  min {g[0]}  median {g[len(g)//2]}  max {g[-1]}")
    if s0:
        b = sorted(x[1] for x in v if x[0] == 0 and x[1] is not None); g = sorted(x[2] - x[1] for x in v if x[0] == 0 and x[2] is not None and x[1] is not None)
        if b: print(f"{'':45s}   when |S|==0: best candidate's mismatches median {b[len(b)//2]} max {b[-1]};  runner-up gap median {g[len(g)//2] if g else '-'}")
    if key[0] == 'tied' and s2:
        sz = collections.Counter(x[0] for x in v if x[0] >= 2)
        print(f"{'':45s}   |S| sizes among the absolutely-ambiguous: {dict(sorted(sz.items()))}")
