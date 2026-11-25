#!/usr/bin/env python3
"""PREREG r1950 (docs/PREREG_read_kind_2026-11-25.md): how many of a near-tied read's places to keep beside the primary,
decided by how often the two loci tie over the whole sample.

  read_kind.py --ties TIES.tsv.gz --gtf LOCI.gtf --fit TRUTH,SPLIT,SEL.bam --judge TRUTH,SPLIT,SEL.bam --out OUT.tsv"""
import argparse, collections, gzip, itertools, os, re, sys, bisect
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R
ap = argparse.ArgumentParser(); [ap.add_argument(x, required=True) for x in ('--ties', '--gtf', '--fit', '--judge', '--out')]
a = ap.parse_args()
GRID = (0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)

# loci
span = {}
for l in open(a.gtf):
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'exon':
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); s, e = int(f[3]) - 1, int(f[4])
        c = span.setdefault(g, [f[0], f[6], s, e]); c[2] = min(c[2], s); c[3] = max(c[3], e)
by = collections.defaultdict(list)
for g, (c, st, s, e) in span.items(): by[(c, st)].append((s, e, g))
idx = {}
for k, v in by.items(): v.sort(); idx[k] = ([x[0] for x in v], v, max(e - s for s, e, _ in v))
def locus(rec):
    k = idx.get((rec.chrom, rec.strand))
    if not k: return None
    starts, v, longest = k; best = None; i = bisect.bisect_left(starts, rec.start - longest)
    while i < len(v) and v[i][0] < rec.end:
        o = min(rec.end, v[i][1]) - max(rec.start, v[i][0])
        if o > 0 and (best is None or o > best[0]): best = (o, v[i][2])
        i += 1
    return best and best[1]

# how often two loci tie
placed = collections.Counter()
with open(a.ties + '.placed.tsv') as fh:
    fh.readline()
    for l in fh: g, n = l.split('\t'); placed[g] = int(n)
tied = collections.Counter(); of = collections.defaultdict(set)   # a read's rows are not adjacent in the table (BAM order)
with gzip.open(a.ties, 'rt') as fh:
    fh.readline(); fh.readline()
    for l in fh:
        q, g, gap, _ = l.split('\t')
        if g != 'outside': of[q].add(g)
for loci in of.values():
    for x, y in itertools.combinations(sorted(loci), 2): tied[x, y] += 1
del of
def rate(x, y):
    if x is None or y is None: return None
    if x == y: return 1.0
    t = tied[(x, y) if x < y else (y, x)]; d = t + placed[x] + placed[y]
    return t / d if d else 0.0

def molecules(spec):
    tr, sp, bam = spec.split(','); truth, recs = R.load_truth(tr, sp), R.load_records(bam); out = []
    for q in sorted(set(truth) & set(recs)):
        t, rs = truth[q], recs[q]
        if not (R.near_tied(rs, 'ms') and t.kind in R.JUDGED): continue
        best = max(r.ms for r in rs); prim = [r for r in rs if r.primary]; pl = locus(prim[0]) if prim else None
        rows = []   # (truth place or None, primary, gap, rate to the primary's locus)
        for r in rs:
            if best - r.ms > 10 and not r.primary: continue
            rows.append((R.locus_index(r, t), r.primary, best - r.ms, 1.0 if r.primary else rate(pl, locus(r))))
        out.append((t.kind, t.origin, rows))
    return out

def stmt(rows, rule, tau=None):
    if rule == 'primary': keep = [x for x in rows if x[1]]
    elif rule == 'exact band': keep = [x for x in rows if x[1] or x[2] <= 0]
    elif rule == 'within 10': keep = rows
    elif rule == 'T': keep = [x for x in rows if x[1] or x[2] <= 0 or (x[3] is not None and x[3] >= tau)]
    elif rule == 'X': keep = [x for x in rows if x[1] or (x[3] is not None and x[3] >= tau) or (x[3] is None and x[2] <= 0)]
    return {x[0] for x in keep}

def score(mols, rule, tau=None):
    n = collections.Counter()
    for kind, origin, rows in mols:
        if rule == 'oracle': s = stmt(rows, 'primary' if kind == 'one' else 'within 10')
        else: s = stmt(rows, rule, tau)
        if None in s: n[kind, 'nov'] += 1; continue
        n[kind, 'n'] += 1; n[kind, 'eq'] += s == origin
    return n
def pooled(n): return sum(n[k, 'eq'] for k in R.JUDGED), sum(n[k, 'n'] for k in R.JUDGED)

fit, judge = molecules(a.fit), molecules(a.judge)
L = []; tau = {}
for rule in ('T', 'X'):
    res = [(pooled(score(fit, rule, t))[0], -i, t) for i, t in enumerate(GRID)]
    tau[rule] = max(res)[2]
    L.append(f'# fitting quarter, rule {rule}: pooled equals by tau ' + ', '.join(f'{t}: {k}' for k, _, t in res) + f' -> tau {tau[rule]}')
def line(name, n):
    cells = []
    for ks in (('one',), ('several',), R.JUDGED):
        k, m = sum(n[x, 'eq'] for x in ks), sum(n[x, 'n'] for x in ks); lo, hi = R.wilson(k, m)
        cells.append(f'{k}\t{m}\t{k / m if m else 0:.3f}\t{lo:.3f}\t{hi:.3f}')
    return name + '\t' + '\t'.join(cells) + f"\t{n['one', 'nov'] + n['several', 'nov']}"
for label, mols in (('fitting quarter 0-24', fit), ('JUDGING quarter 25-49', judge)):
    L.append(f'## {label}'); L.append('statement\tone_eq\tone_n\tone\tlo\thi\tsev_eq\tsev_n\tseveral\tlo\thi\tpool_eq\tpool_n\tpooled\tlo\thi\tno_verdict')
    for rule in ('primary', 'exact band', 'within 10'): L.append(line(rule, score(mols, rule)))
    for rule in ('T', 'X'): L.append(line(f'{rule}(tau {tau[rule]})', score(mols, rule, tau[rule])))
    L.append(line('oracle (kind known)', score(mols, 'oracle')))
    # raw separating power: records beside the primary, truly in the origin set, by rate bin
    bins = collections.Counter()
    for kind, origin, rows in mols:
        for i, prim, gap, r in rows:
            if prim or i is None: continue
            b = 'no locus' if r is None else next(f'<{t}' for t in (*GRID, 1.01) if r < t)
            bins['exact' if gap <= 0 else 'gap 1-10', b, 'n'] += 1; bins['exact' if gap <= 0 else 'gap 1-10', b, 'in'] += i in origin
    L.append('records beside the primary: class\trate bin\tin the origin set\trecords\tshare')
    for c in ('exact', 'gap 1-10'):
        for b in [f'<{t}' for t in (*GRID, 1.01)] + ['no locus']:
            if bins[c, b, 'n']: L.append(f"  {c}\t{b}\t{bins[c, b, 'in']}\t{bins[c, b, 'n']}\t{bins[c, b, 'in'] / bins[c, b, 'n']:.3f}")
open(a.out, 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
