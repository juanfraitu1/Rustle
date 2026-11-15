#!/usr/bin/env python3
"""Step 0 for the nascent row (2026-11-11): size of StringTie's `mark_if_nascent` class in an assembly (rlink.cpp, -N).

  nascent_class.py GTF REFSEQ_GTF [TOOL_GTF ...]      (selftest: nascent_class.py selftest)

p is a nascent of n when both are on one strand, n has more exons and higher cov, p starts no later than the end of n's first
exon (5' side), p's introns are n's first introns, and p runs past the end of n's next exon boundary into the following intron."""
import sys, re, gzip, collections

def models(lines):
    ex = collections.defaultdict(list); cov = {}
    for l in lines:
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        k = (f[0], f[6], re.search(r'transcript_id "([^"]*)"', f[8]).group(1))
        if f[2] == 'exon': ex[k].append((int(f[3]) - 1, int(f[4])))
        else:
            m = re.search(r'cov "([^"]*)"', f[8]); cov[k] = float(m.group(1)) if m else 0.0
    return [(k[0], k[1], sorted(e), cov.get(k, 0.0)) for k, e in ex.items()]
def flip(e): return [(-b, -a) for a, b in reversed(e)]  # '-' strand: mirror, so 5' is always first
def is_nascent(p, n):
    k = len(p) - 1
    if len(n) <= len(p) or p[0][0] >= n[0][1]: return False
    if any(p[i][1] != n[i][1] or p[i + 1][0] != n[i + 1][0] for i in range(k)): return False
    return p[-1][1] > n[k][1]
def mark(ms):
    by = collections.defaultdict(list); out = []
    for c, s, e, cov in ms:
        if s in '+-': by[(c, s)].append((e if s == '+' else flip(e), cov, e))
    for (c, s), L in by.items():
        first = collections.defaultdict(list)  # 5' donor (or contig-wide for single-exon p) -> candidates
        for e, cov, _ in L:
            if len(e) > 1: first[e[0][1]].append((e, cov))
        donors = sorted(first)
        for e, cov, raw in L:
            if len(e) > 1: cand = first[e[0][1]]
            else:
                import bisect
                i = bisect.bisect_right(donors, e[0][0]); cand = [x for d in donors[i:] if d < e[0][1] for x in first[d]]
            if any(nc > cov and is_nascent(e, ne) for ne, nc in cand): out.append((c, s, raw))
    return out
chain = lambda c, s, e: (c, s, tuple((a[1], b[0]) for a, b in zip(e, e[1:])))
def read(p): return models((gzip.open if p.endswith('.gz') else open)(p, 'rt'))

if sys.argv[1] == 'selftest':
    g = lambda t, s, cov, exons: [f'c\tx\ttranscript\t1\t9\t.\t{s}\t.\ttranscript_id "{t}"; cov "{cov}";'] + [f'c\tx\texon\t{a + 1}\t{b}\t.\t{s}\t.\ttranscript_id "{t}";' for a, b in exons]
    N = [(0, 10), (20, 30), (40, 50)]
    ms = models(g('n', '+', 9, N) + g('p1', '+', 1, [(0, 10), (20, 35)]) + g('p0', '+', 1, [(2, 15)]) + g('frag', '+', 1, [(0, 10), (20, 28)]) + g('hi', '+', 99, [(0, 10), (20, 36)])
                + g('m', '-', 9, N) + g('q', '-', 1, [(15, 30), (40, 50)]) + g('notq', '-', 1, [(0, 10), (20, 35)]))
    got = sorted(tuple(e) for _, _, e in mark(ms)); assert got == sorted([((0, 10), (20, 35)), ((2, 15),), ((15, 30), (40, 50))]), got
    print('ok'); sys.exit()
ms = read(sys.argv[1]); ref = {chain(c, s, e) for c, s, e, _ in read(sys.argv[2]) if len(e) > 1}; nas = mark(ms)
multi = {chain(c, s, e) for c, s, e in nas if len(e) > 1}; allm = {chain(c, s, e) for c, s, e, _ in ms if len(e) > 1}
print(f'models {len(ms)}; nascent-marked {len(nas)} (single-exon {sum(len(e) == 1 for _, _, e in nas)}, multi-exon chains {len(multi)})')
print(f'multi-exon marked: RefSeq-correct {len(multi & ref)} = {len(multi & ref) / max(len(multi), 1):.3f}; all chains {len(allm)} correct {len(allm & ref)} = {len(allm & ref) / max(len(allm), 1):.3f}')
for p in sys.argv[3:]:
    t = {chain(c, s, e) for c, s, e, _ in read(p) if len(e) > 1}; print(f'  marked chains also in {p.split("/")[-1]}: {len(multi & t)} of {len(multi)}')
