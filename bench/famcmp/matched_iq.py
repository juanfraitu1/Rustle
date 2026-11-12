# matched_iq.py REF OURS ST IQ_GTF IQ_COUNTS BAMBU_GTF BAMBU_TSV : precision of each curve at IsoQuant's correct counts (first cut reaching the count)
import sys, re, collections  # run with PYTHONPATH or from bench/famcmp (imports curve)
from curve import load
ref = set(load(sys.argv[1]))
def by_id(gtf, score, worst):
    ex = collections.defaultdict(list); d = {}
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    for (c, s, t), e in ex.items():
        e.sort()
        if len(e) > 1:
            ch = (c, s, tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))); d[ch] = max(d.get(ch, worst), score.get(t, worst))
    return d
rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[5]) if not l.startswith('__')]
iq = by_id(sys.argv[4], {r[0]: float(r[1]) for r in rows[1:]}, 0.0)
rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[7])]; h = rows[0]
bb = by_id(sys.argv[6], {r[h.index('TXNAME')]: -float(r[h.index('NDR')]) for r in rows[1:] if r[h.index('NDR')] != 'NA'}, -9.0)
arms = {'isoquant pt': iq, 'ours': {c: v[0] for c, v in load(sys.argv[2], 'reads').items()}, 'stringtie': {c: v[0] for c, v in load(sys.argv[3], 'cov').items()}, 'bambu': bb}
def at(d, k):  # strictest cut (score >= t) with at least k correct; ties enter together
    n = ok = 0; out = None; items = sorted(d.items(), key=lambda x: -x[1])
    for i, (c, v) in enumerate(items):
        n += 1; ok += c in ref
        if ok >= k and (i + 1 == len(items) or items[i + 1][1] != v): return f'{ok / n:.3f} ({ok})'
    return '-'
assert at({('a',): 2, ('b',): 2, ('c',): 1}, 0) == '0.000 (0)'  # ties enter together
print('K\t' + '\t'.join(arms))
for t in (0, 3, 6, 11, 21):
    k = sum(c in ref for c, v in iq.items() if v >= t)
    print(f'{k} (count>={t})\t' + '\t'.join(at(d, k) for d in arms.values()))
