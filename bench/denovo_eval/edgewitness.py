#!/usr/bin/env python3
"""Row 1680: a second sample as witness for a refused family pair (docs/PREREG_cross_sample_witness_2026-11-13.md).

  edgewitness.py TRAIN_A.tsv TRAIN_B.tsv TEST_A.tsv TEST_B.tsv      edgewitness.py --selftest
Tables are edgefeat.py outputs. Row-1675 model; witness score = min of a gene-name pair's scores in the two samples.
"""
import sys, zlib
import numpy as np
import edgerules as R, edgescore as S

key = lambda r: '|'.join(sorted((r[2], r[3])))


def witness(ka, sa, ya, kb, sb, yb):
    """pairs present in both samples -> (keys, min score, truth in either)"""
    a, b = {}, {}
    for d, ks, s, y in ((a, ka, sa, ya), (b, kb, sb, yb)):
        for k, v, t in zip(ks, s, y):
            if k not in d or v > d[k][0]: d[k] = (v, bool(t))   # several locus pairs of one name pair: the best
    ks = sorted(a.keys() & b.keys())
    return ks, np.array([min(a[k][0], b[k][0]) for k in ks]), np.array([a[k][1] or b[k][1] for k in ks], bool)


def selftest():
    ks, s, y = witness(['x', 'y', 'x'], [1., 5., 3.], [1, 0, 1], ['x', 'z'], [2., 9.], [0, 1])
    assert ks == ['x'] and s.tolist() == [2.] and y.tolist() == [True]
    print('selftest ok')


def main(a):
    T = [R.load(p) for p in a]; feats = T[0][0]
    K = [[key(r) for r in t[3]] for t in T]; X = [S.prep(t[0], t[1]) for t in T]; Y = [t[2] for t in T]
    both = set(K[0]) & set(K[1]); n0 = len({k for i in (0, 1) for k, y in zip(K[i], Y[i]) if y and k in both})
    print(f'step 0: pairs refused in both train samples {len(both)}; truth pairs among them {n0} -> {"GO" if n0 >= S.MINT else "STOP"}')
    if n0 < S.MINT: return
    Xp, yp, kp = np.vstack(X[:2]), np.concatenate(Y[:2]), K[0] + K[1]
    fold = np.array([zlib.crc32(k.encode()) % S.FOLDS for k in kp]); oof = np.zeros(len(yp))
    for f in range(S.FOLDS):
        oof[fold == f] = S.predict(S.fit(Xp[fold != f], yp[fold != f].astype(float)), Xp[fold == f])
    n = len(K[0]); ks, s, y = witness(K[0], oof[:n], Y[0], K[1], oof[n:], Y[1]); th = S.threshold(s, y)
    print(f'out-of-fold witness population {len(ks)}, truth {int(y.sum())}; threshold (score, true, false): {th}')
    if th is None:
        o = np.argsort(-s); t = np.cumsum(y[o])
        for m in (20, 40):
            if t[-1] >= m: i = np.searchsorted(t, m); print(f'  at {m} true: false {i + 1 - m} precision {m / (i + 1):.3f}')
        print('STOP'); return
    m = S.fit(Xp, yp.astype(float))
    ks, s, y = witness(K[2], S.predict(m, X[2]), Y[2], K[3], S.predict(m, X[3]), Y[3]); k = s >= th[0]
    t, f = int((k & y).sum()), int((k & ~y).sum())
    print(f'test: witness population {len(ks)}, truth {int(y.sum())}; picked true {t} false {f} precision {t / max(1, t + f):.3f} -> {"PASS" if t >= 20 and t / max(1, t + f) >= 0.60 else "FAIL"}')
    for kk, sc, yy in sorted(zip(ks, s, y), key=lambda v: -v[1]):
        if sc >= th[0]: print(f'  {kk}\t{sc:.3f}\t{int(yy)}')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(sys.argv[1:])
