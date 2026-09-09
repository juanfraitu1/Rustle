#!/usr/bin/env python3
"""Genome-wide subfamily scan: run the identity-gap detector on every family.

Builds per-family within-family pairwise identities by streaming the all-vs-all PAF and resolving
each PAF key (contig:start-end of a gene span) to a catalog member by interval overlap, then applies
`identity_gap`'s statistic + nulls to each family with enough pairs.

The point is CALIBRATION: a p-value on one family is uninterpretable, but the distribution of p over
hundreds of families says whether a given family's gap is extreme.

usage: gw_subfamily_scan.py <clusters.tsv> <paf> <out.tsv> [--min-pairs 10]
"""
import sys, re, bisect, collections, math, random, statistics as st, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import modality as M

clusters_p, paf_p, out_p = sys.argv[1:4]
MINP = int(sys.argv[sys.argv.index('--min-pairs') + 1]) if '--min-pairs' in sys.argv else 10

# ---------------------------------------------------------------- members
import csv
byc = collections.defaultdict(list)       # contig -> [(start, end, family)]
fam_n = collections.Counter()
for r in csv.DictReader(open(clusters_p), delimiter='\t'):
    byc[r['chrom']].append((int(r['start']), int(r['end']), r['cluster_id']))
    fam_n[r['cluster_id']] += 1
for c in byc:
    byc[c].sort()
starts = {c: [x[0] for x in v] for c, v in byc.items()}
print(f'{sum(len(v) for v in byc.values())} members in {len(fam_n)} families', file=sys.stderr)

KEY = re.compile(r'^(.+):(\d+)-(\d+)$')
cache = {}
def resolve(tok):
    """Which catalog member does this PAF key overlap by >50% of the shorter interval?"""
    if tok in cache:
        return cache[tok]
    m = KEY.match(tok)
    r = None
    if m:
        c, s, e = m.group(1), int(m.group(2)), int(m.group(3))
        v = byc.get(c)
        if v:
            i = bisect.bisect_right(starts[c], e)
            for ms, me, fam in v[max(0, i - 40):i]:
                ov = min(e, me) - max(s, ms)
                if ov > 0.5 * min(e - s, me - ms):
                    r = (fam, ms, me)
                    break
    cache[tok] = r
    return r

# ---------------------------------------------------------------- pairs
best = collections.defaultdict(dict)      # family -> {(a,b): identity}
seen = 0
for line in open(paf_p):
    f = line.split('\t', 12)
    if len(f) < 12:
        continue
    a, b = resolve(f[0]), resolve(f[5])
    if not a or not b or a[0] != b[0] or a[1:] == b[1:]:
        continue
    bl = int(f[10])
    if bl < 1000:
        continue
    idn = int(f[9]) / bl
    k = tuple(sorted((a[1:], b[1:])))
    d = best[a[0]]
    if idn > d.get(k, 0):
        d[k] = idn
    seen += 1
print(f'{seen} in-family alignments over {len(best)} families', file=sys.stderr)

# ---------------------------------------------------------------- statistic
def largest_gap(xs):
    lo, hi = max(1, int(0.1 * len(xs))), min(len(xs) - 1, int(0.9 * len(xs)))
    b, at = 0.0, None
    for i in range(lo, hi):
        d = xs[i] - xs[i - 1]
        if d > b:
            b, at = d, (xs[i - 1] + xs[i]) / 2
    return b, at

def pvals(v, reps=2000, seed=0):
    n = len(v)
    obs, cut = largest_gap(v)
    mu, sd = st.mean(v), st.pstdev(v)
    if sd <= 0:
        return obs, cut, {}
    rnd = random.Random(seed)
    out = {}
    t = mu * (1 - mu) / (sd * sd) - 1 if 0 < mu < 1 else -1
    if t > 0:
        al, be = mu * t, (1 - mu) * t
        out['beta'] = sum(1 for _ in range(reps)
                          if largest_gap(sorted(rnd.betavariate(al, be) for _ in range(n)))[0] >= obs) / reps
    q = st.quantiles(v, n=4)
    h = 0.9 * min(sd, (q[2] - q[0]) / 1.34) * n ** -0.2
    if h > 0:
        out['smooth'] = sum(1 for _ in range(reps)
                            if largest_gap(sorted(min(1.0, max(0.0, rnd.choice(v) + rnd.gauss(0, h)))
                                                  for _ in range(n)))[0] >= obs) / reps
    return obs, cut, out

# dump every in-family pair once, so any later statistic reads this instead of the 1 GB PAF
with open(out_p + '.pairs.tsv', 'w') as ph:
    ph.write('family\ta\tb\tidentity\n')
    for fam, d in sorted(best.items()):
        for (a, b), idn in d.items():
            ph.write(f'{fam}\t{a[0]}-{a[1]}\t{b[0]}-{b[1]}\t{idn:.6f}\n')

with open(out_p, 'w') as fh:
    fh.write('family\tn_members\tn_pairs\tmin_id\tmedian_id\tmax_id\tgap\tcut\tp_beta\tp_smooth\tp_worst\tsilverman_p\th_crit\tdbic\n')
    done = 0
    for fam, d in sorted(best.items()):
        v = sorted(d.values())
        if len(v) < MINP:
            continue
        obs, cut, pv = pvals(v)
        if not pv:
            continue
        pw = max(pv.values())
        sp, hc = M.silverman_p(v, B=200)
        dbic = M.mixture_dbic(v)
        fh.write(f'{fam}\t{fam_n[fam]}\t{len(v)}\t{v[0]:.4f}\t{st.median(v):.4f}\t{v[-1]:.4f}\t'
                 f'{obs:.4f}\t{cut:.4f}\t{pv.get("beta", float("nan")):.4f}\t'
                 f'{pv.get("smooth", float("nan")):.4f}\t{pw:.4f}\t{sp:.4f}\t{hc:.4f}\t{dbic:.1f}\n')
        done += 1
        if done % 50 == 0:
            print(f'  {done} families scored', file=sys.stderr)
print(f'wrote {out_p}: {done} families', file=sys.stderr)
