#!/usr/bin/env python3
"""Row 1681: family topology as the witness for a refused pair (docs/PREREG_family_topology_2026-11-13.md).

  edgetopo.py TRAIN.tsv[,TRAIN.tsv] TEST.tsv[,TEST.tsv]      edgetopo.py --selftest
"""
import sys, collections
ARMS = [(k, f) for k in (2, 3) for f in (0.5, 0.67, 1.0)]


def pairs(rows):
    """rows (gene_a, gene_b, truth, linked, stage) -> [(truth, {direction: (support, frac, family size)})] per refused pair"""
    fam, nb = collections.defaultdict(set), collections.defaultdict(set)
    for a, b, t, l, st in rows:
        nb[a].add(b); nb[b].add(a)
        if l == '1': fam[a].add(b); fam[b].add(a)
    out = []
    for a, b, t, l, st in rows:
        if l == '1' or st == 'edge': continue
        d = []
        for x, y in ((a, b), (b, a)):
            F = fam[y] | {y}; s = len(F & nb[x]); d.append((s, s / len(F), len(F)))
        out.append((t == '1', d))
    return out


def count(ps, k, f):
    t = n = 0
    for tr, d in ps:
        if any(sz >= 2 and s >= k and fr >= f for s, fr, sz in d): n += 1; t += tr
    return t, n - t


def load(p):
    rows = [l.rstrip('\n').split('\t') for l in open(p)]; h = rows[0]; c = [h.index(x) for x in ('gene_a', 'gene_b', 'truth', 'linked', 'stage')]
    return pairs([[r[i] for i in c] for r in rows[1:]])


def selftest():
    rows = [('b', 'c', '1', '1', 'edge'), ('b', 'd', '1', '1', 'edge'), ('a', 'b', '1', '0', 'low_share'), ('a', 'c', '1', '0', 'low_share'), ('x', 'b', '0', '0', 'low_share')]
    ps = pairs(rows); assert count(ps, 2, 0.5) == (2, 0) and count(ps, 3, 0.5) == (0, 0) and count(ps, 2, 1.0) == (1, 0), ps   # family(c) = {c, b}: direct links only
    print('selftest ok')


def main(a):
    tr = [p for f in a[0].split(',') for p in load(f)]; best = None
    print(f'train refused pairs {len(tr)}, truth {sum(t for t, _ in tr)}')
    for k, f in ARMS:
        t, w = count(tr, k, f); pr = t / max(1, t + w); print(f'  K {k} F {f}: true {t} false {w} precision {pr:.3f}')
        if pr >= 0.75 and t >= 20 and (best is None or t > best[0]): best = (t, k, f)
    if best is None: print('STOP: no arm reaches precision 0.75 with 20 true'); return
    T = W = 0
    for f in a[1].split(','):
        t, w = count(load(f), best[1], best[2]); T += t; W += w; print(f'  test {f}: true {t} false {w}')
    print(f'arm K {best[1]} F {best[2]}; pooled test true {T} false {W} precision {T / max(1, T + W):.3f} -> {"PASS" if T >= 20 and T / max(1, T + W) >= 0.60 else "FAIL"}')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(sys.argv[1:])
