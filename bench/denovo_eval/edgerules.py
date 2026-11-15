#!/usr/bin/env python3
"""Row 1674: search a conjunction of <= 3 threshold conditions over edgefeat.py features on one sample, score it frozen on others.

  edgerules.py TRAIN.tsv TEST.tsv [TEST2.tsv ...]      edgerules.py --apply 'f>=x&g<=y' TSV...      edgerules.py --selftest
Population: gene pairs with a record, not linked, not at stage `edge`. Thresholds = deciles of the TRAIN true pairs.
Selected rule = most TRAIN true pairs among conjunctions with TRAIN precision >= 0.80 and >= 20 true (ties: fewer conditions).
"""
import sys, itertools
import numpy as np
MINP, MINT = 0.80, 20
HOMOLOG = False


def load(p):
    rows = [l.rstrip('\n').split('\t') for l in open(p)]
    h = rows[0]; rows = [r for r in rows[1:] if r[h.index('linked')] == '0' and r[h.index('stage')] != 'edge']
    end = h.index('locus_a'); feats = h[7:end]
    if HOMOLOG:   # row 1676: positives = truth or related; undecidable pairs leave the population
        rows = [r for r in rows if r[4] == '1' or r[h.index('related')] != 'NA']
    X = np.array([[float(v) for v in r[7:end]] for r in rows]).reshape(len(rows), len(feats))
    y = np.array([r[4] == '1' or (HOMOLOG and r[h.index('related')] == '1') for r in rows])
    return feats, X, y, rows


def conds(feats, X, y):
    out = []
    for j, f in enumerate(feats):
        v = X[y, j]
        if len(v) == 0: continue
        for t in sorted(set(np.percentile(v, range(10, 100, 10)).round(4))):
            out.append((f, '>=', t, X[:, j] >= t)); out.append((f, '<=', t, X[:, j] <= t))
    return out


def score(m, y):
    t = int((m & y).sum()); n = int(m.sum())
    return t, n - t, (t / n if n else 0.0)


def search(feats, X, y):
    C = conds(feats, X, y); cand = []
    for c in C: cand.append(((c,), c[3]))
    for a, b in itertools.combinations(C, 2):
        if a[0] != b[0] or a[1] != b[1]: cand.append(((a, b), a[3] & b[3]))
    two = sorted((x for x in cand if len(x[0]) == 2 and score(x[1], y)[0] >= MINT), key=lambda x: -score(x[1], y)[2])[:200]
    for cs, m in two:
        for c in C:
            if all(c[0] != d[0] or c[1] != d[1] for d in cs): cand.append((cs + (c,), m & c[3]))
    ok = [(score(m, y), cs) for cs, m in cand]
    ok = [x for x in ok if x[0][2] >= MINP and x[0][0] >= MINT]
    ok.sort(key=lambda x: (-x[0][0], len(x[1]), -x[0][2]))
    return ok


def text(cs): return '&'.join(f'{f}{op}{t:g}' for f, op, t, *_ in cs)


def apply(rule, feats, X):
    m = np.ones(len(X), bool)
    for c in rule.split('&'):
        op = '>=' if '>=' in c else '<='; f, t = c.split(op)
        m &= (X[:, feats.index(f)] >= float(t)) if op == '>=' else (X[:, feats.index(f)] <= float(t))
    return m


def selftest():
    X = np.array([[1, 0], [1, 1], [0, 1], [0, 0]] * 30, float); y = np.array([True, False, False, False] * 30)
    r = search(['a', 'b'], X, y)[0]
    assert r[0][:2] == (30, 0) and set(text(r[1]).split('&')) == {'a>=1', 'b<=0'}, r
    assert apply('a>=1&b<=0', ['a', 'b'], X).sum() == 30
    print('selftest ok')


def main(a):
    if a[0] == '--apply':
        for p in a[2:]:
            feats, X, y, rows = load(p); m = apply(a[1], feats, X)
            print(p, 'population', len(y), 'true', int(y.sum()), 'rule: true %d false %d precision %.3f' % score(m, y))
        return
    feats, X, y, _ = load(a[0])
    print(f'train {a[0]}: population {len(y)}, true {int(y.sum())}')
    ok = search(feats, X, y)
    if not ok:
        print('no conjunction reaches the training floor'); return
    for (t, f, p), cs in ok[:5]: print(f'  train true {t} false {f} precision {p:.3f}  {text(cs)}')
    rule = text(ok[0][1]); print('SELECTED', rule)
    for p in a[1:]:
        feats2, X2, y2, _ = load(p); m = apply(rule, feats2, X2)
        print(f'  test {p}: population {len(y2)}, true {int(y2.sum())}; rule true %d false %d precision %.3f' % score(m, y2))


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(sys.argv[1:])
