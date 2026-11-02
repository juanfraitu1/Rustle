#!/usr/bin/env python3
"""5' / 3' end distance of exact intron-chain matches (gffcompare `=`) to their reference transcript, per assembly (report only,
2026-11-02: does coverage-based end trimming have anything left to fix?). Signed in transcript orientation: + = the model extends
past the reference end, - = it stops short. Also on the reference transcripts matched by EVERY assembly (same denominator).
  end_accuracy.py REF.gtf NAME=QUERY.gtf:QUERY.tmap ..."""
import collections
import csv
import re
import statistics
import sys


def spans(path):
    s = {}
    for l in open(path):
        f = l.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        a, b = int(f[3]), int(f[4])
        c = s.get(t)
        s[t] = (f[6], min(a, c[1]) if c else a, max(b, c[2]) if c else b)
    return s


def dist(q, r):
    st, qa, qb = q
    _, ra, rb = r
    return (ra - qa, qb - rb) if st != '-' else (qb - rb, ra - qa)  # (5', 3')


def summary(d):
    a = [abs(x) for x in d]
    return (f"n {len(d)}  median |d| {statistics.median(a):.0f}  <=10 {sum(x <= 10 for x in a) / len(a):.3f}  <=50 "
            f"{sum(x <= 50 for x in a) / len(a):.3f}  >200 {sum(x > 200 for x in a) / len(a):.3f}  longer {sum(x > 0 for x in d) / len(d):.3f}")


ref = spans(sys.argv[1])
runs = {}
for arg in sys.argv[2:]:
    name, rest = arg.split('=', 1)
    g, tm = rest.split(':')
    q = spans(g)
    best = {}
    for r in csv.DictReader(open(tm), delimiter='\t'):
        if r['class_code'] == '=' and r['ref_id'] in ref and r['qry_id'] in q:
            d = dist(q[r['qry_id']], ref[r['ref_id']])
            # one model per reference transcript: the closest by summed |d|
            if r['ref_id'] not in best or sum(map(abs, d)) < sum(map(abs, best[r['ref_id']])):
                best[r['ref_id']] = d
    runs[name] = best
common = set.intersection(*(set(v) for v in runs.values()))
for name, best in runs.items():
    for label, keys in (('all', best), ('common', common)):
        v = [best[k] for k in keys]
        print(f"{name}\t{label}\t5'\t{summary([x[0] for x in v])}")
        print(f"{name}\t{label}\t3'\t{summary([x[1] for x in v])}")
