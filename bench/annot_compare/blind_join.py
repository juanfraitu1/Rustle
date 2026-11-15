#!/usr/bin/env python3
"""Row 1712 (docs/PREREG_blind_join_tag_2026-11-14.md): score a blind gene-join tag (x.gtf only) with row 1707's bars.
usage: blind_join.py X.gtf REFSEQ.gtf[.gz] JOINS.tsv [TANDEM_FUSED.tsv]   (JOINS.tsv from join_truth.py)"""
import bisect, collections, csv, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fusion_audit as FA
from join_truth import tx_table
from join_support import crossing
from cat_vs_refseq_lib import Idx, share

FORMS = ('F1', 'F2', 'F3', 'F4')
ov = lambda x, y: any(a[0] < b[1] and b[0] < a[1] for a in x for b in y)


def subchain(a, t):
    """a's intron chain is a contiguous slice of t's."""
    if not a or a[0] not in t: return not a
    i = t.index(a[0]); return t[i:i + len(a)] == a


def blind(X):
    """transcript -> set of forms that tag it."""
    grp = collections.defaultdict(list); out = collections.defaultdict(set)
    for t, v in X.items():
        if len(v['ex']) > 1: grp[(v['c'], v['s'])].append(t)
    for key, ts in sorted(grp.items()):
        t0 = time.time(); byend = sorted(ts, key=lambda t: X[t]['ex'][-1][1]); ends = [X[t]['ex'][-1][1] for t in byend]
        bystart = sorted(ts, key=lambda t: X[t]['ex'][0][0]); starts = [X[t]['ex'][0][0] for t in bystart]; cache = {}
        for t in ts:
            v = X[t]; ch = v['chain'][2]
            for i in range(len(v['ex']) - 1):
                L, R = v['ex'][i], v['ex'][i + 1]; s, e = L[1], R[0]
                if (L, R) not in cache:
                    A = [a for a in byend[bisect.bisect_right(ends, L[0]):bisect.bisect_left(ends, e)] if ov(X[a]['ex'], [L])]
                    B = [b for b in bystart[bisect.bisect_right(starts, s):bisect.bisect_left(starts, R[1])] if ov(X[b]['ex'], [R])]
                    cache[(L, R)] = (A, B)
                A, B = ([x for x in S if x != t and not subchain(X[x]['chain'][2], ch)] for S in cache[(L, R)])
                if not A or not B: continue
                pair = lambda A, B: any(not ov(X[a]['ex'], X[b]['ex']) for a in A for b in B)
                if not pair(A, B): continue
                out[t].add('F1'); rd = lambda x: X[x]['reads']
                if pair([a for a in A if rd(a) >= 5], [b for b in B if rd(b) >= 5]): out[t].add('F2')
                if pair([a for a in A if rd(a) >= v['reads']], [b for b in B if rd(b) >= v['reads']]): out[t].add('F3')
                if pair([a for a in A if X[a]['ex'][-1][0] >= s], B) or pair(A, [b for b in B if X[b]['ex'][0][1] <= e]): out[t].add('F4')
        print(f'[blind] {key[0]}{key[1]} {len(ts)} transcripts {time.time() - t0:.0f} s', file=sys.stderr, flush=True)
    return out


def main(xgtf, rs, joins, fused=None, tagger=None):
    """tagger: X -> {name: set of transcripts}, default the model forms of row 1712; returns the tag sets."""
    X = tx_table(xgtf); R, rmeta = FA.load(rs, 'gene_id'); name = {g: rmeta[g][2].removeprefix('gene-') for g in R}; ri = Idx(R)
    kind = {}; own = collections.defaultdict(set)  # gene -> non-joining transcripts
    for t, v in X.items():
        hits = ri.hits(v['c'], v['s'], [list(x) for x in v['ex']]); cr = crossing(v['ex'], hits, R)
        dis = any(not share(R[a][2], R[b][2]) for i, a in enumerate(hits) for b in hits[i + 1:])
        on = [{g for g in hits if share([list(e)], R[g][2])} for e in v['ex']]
        kind[t] = 'join' if cr else 'inside one gene' if on and set.intersection(*on) else 'touches two genes, no crossing intron' if dis else 'exon outside any gene'
        if not dis:
            for g in hits: own[name[g]].add(t)
    J = list(csv.DictReader(open(joins), delimiter='\t')); ab = [r for r in J if r['class'] in ('a', 'b1')]
    has = lambda gs: any(own[g] for g in gs.split(',') if g)
    ceil = sum(has(r['left_genes']) and has(r['right_genes']) for r in ab)
    print(f'{xgtf}\nstep 0 ceiling: a + b1 joining transcripts with a non-joining model on both flanking genes {ceil} of {len(ab)} = {ceil / max(len(ab), 1):.3f}')
    tags = {'current fusion_of': {t for t, v in X.items() if v['tag']}}
    if tagger: tags.update(tagger(X))
    else:
        bl = blind(X)
        for f in FORMS: tags[f] = {t for t, fs in bl.items() if f in fs}
    fj = []
    if fused:
        want = set()
        for l in list(open(fused))[1:]:
            f = l.split()
            if f[6] != 'NA': c, se = f[6].split(':'); want.add((c, *map(int, se.split('-'))))
        have = collections.defaultdict(set)
        for t, v in X.items():
            for a, b in v['chain'][2]:
                for w in ((v['c'], a + da, b + db) for da in (-1, 0, 1) for db in (-1, 0, 1)):
                    if w in want: have[w].add(t)
        fj = list(have.values()); print(f'fused tandem junctions listed {len(want)}, present in this x.gtf {len(fj)}')
    print('tag\ttagged\tinside_one_gene\tshare(bar<=0.10)\tinside/(inside+join)\tab_untagged\tshare(bar<=0.50)\tfused_junctions_tagged\tkinds\ttag_rate_by_class\tgenes_without_untagged_model')
    for n, T in tags.items():
        k = collections.Counter(kind[t] for t in T); un = sum(r['transcript'] not in T for r in ab)
        rate = {c: f"{sum(r['transcript'] in T for r in J if r['class'] == c)}/{sum(r['class'] == c for r in J)}" for c in ('a', 'b1', 'b2', 'b0', 'w', 'c')}
        lost = sum(1 for g, m in own.items() if m and m <= T)
        print(n, len(T), k['inside one gene'], f"{k['inside one gene'] / max(len(T), 1):.3f}", f"{k['inside one gene'] / max(k['inside one gene'] + k['join'], 1):.3f}", f'{un}/{len(ab)}', f'{un / max(len(ab), 1):.3f}', f'{sum(bool(h & T) for h in fj)}/{len(fj)}' if fused else 'NA', dict(k), rate, lost, sep='\t', flush=True)
    return tags


def _selfcheck():
    assert subchain(((5, 9),), ((1, 3), (5, 9))) and not subchain(((5, 8),), ((1, 3), (5, 9))) and subchain((), ((1, 2),))
    m = lambda ex, r=9: dict(c='c', s='+', reads=r, ex=ex, chain=('c', '+', tuple((a[1], b[0]) for a, b in zip(ex, ex[1:]))))
    X = {'T': m([(0, 10), (20, 30), (100, 110), (120, 130)], 3), 'a': m([(0, 10), (20, 30), (40, 50)]), 'b': m([(80, 90), (100, 110), (120, 130)]),
         'frag': m([(100, 110), (120, 130)]), 'U': m([(200, 210), (220, 230), (300, 310)]), 'u1': m([(200, 210), (220, 230)])}
    e, sys.stderr = sys.stderr, open(os.devnull, 'w')
    try: o = blind(X)
    finally: sys.stderr = e
    assert o['T'] == set(FORMS) and 'U' not in o and 'a' not in o, dict(o)  # U has only truncated copies of itself


if __name__ == '__main__':
    _selfcheck(); main(*sys.argv[1:5])
