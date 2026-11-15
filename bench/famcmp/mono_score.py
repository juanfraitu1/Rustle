#!/usr/bin/env python3
"""Single-exon models against RefSeq. usage: mono_score.py REF_GTF NAME=GTF ...
A reference single-exon transcript is FOUND by a single-exon model that overlaps it by >= 50% of both lengths (strand must
agree unless the model has none). Every single-exon model gets one class: single-exon ref match / inside one exon of a
multi-exon ref transcript (>= 80% of the model) / overlaps a ref gene otherwise / intergenic."""
import sys, re, gzip, collections, bisect
op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
def load(path):
    T = collections.defaultdict(list)
    for l in op(path):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': T[(f[0], f[6], re.search(r'transcript_id "([^"]*)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    return T
class Idx:  # intervals per chrom in 10 kb bins
    def __init__(s): s.b = collections.defaultdict(list)
    def add(s, c, a, b, v):
        for k in range(a // 10000, b // 10000 + 1): s.b[(c, k)].append((a, b, v))
    def hits(s, c, a, b):
        return {x for k in range(a // 10000, b // 10000 + 1) for x in s.b.get((c, k), ()) if x[0] <= b and a <= x[1]}
ref = load(sys.argv[1]); mono = Idx(); exon = Idx(); span = Idx(); nref = 0
for (c, s, t), e in ref.items():
    e.sort(); span.add(c, e[0][0], e[-1][1], s)
    if len(e) == 1: mono.add(c, e[0][0], e[0][1], (s, t)); nref += 1
    else:
        for a, b in e: exon.add(c, a, b, s)
ok = lambda s, m: m not in '+-' or s == m
print(f'reference single-exon transcripts {nref}')
print('tool\tsingle-exon models\tref single-exon found\tmodels: ref match\tinside a multi-exon exon\tother genic\tintergenic')
for arg in sys.argv[2:]:
    name, p = arg.split('=', 1); found = set(); C = collections.Counter()
    for (c, s, t), e in load(p).items():
        if len(e) != 1: continue
        a, b = e[0]; L = b - a + 1
        m = [(x, y, v) for x, y, v in mono.hits(c, a, b) if ok(v[0], s) and min(b, y) - max(a, x) + 1 >= 0.5 * max(L, y - x + 1)]
        if m: found.update(v[1] for x, y, v in m); C['match'] += 1
        elif any(ok(v, s) and min(b, y) - max(a, x) + 1 >= 0.8 * L for x, y, v in exon.hits(c, a, b)): C['exon'] += 1
        elif any(ok(v, s) for x, y, v in span.hits(c, a, b)): C['genic'] += 1
        else: C['inter'] += 1
    print(name, sum(C.values()), len(found), C['match'], C['exon'], C['genic'], C['inter'], sep='\t')
