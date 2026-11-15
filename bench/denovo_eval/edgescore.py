#!/usr/bin/env python3
"""Row 1675: one weighted score (L2 logistic, IRLS) over edgefeat.py features (docs/PREREG_edge_score_2026-11-13.md).

  edgescore.py --train A.tsv,H.tsv --test OR.tsv,KB.tsv [--out-prefix P]      edgescore.py --selftest
Population as edgerules.py. Count-like features enter as log1p. Threshold = the lowest score whose 5-fold out-of-fold
precision on the training pool is >= 0.75 with the most true pairs; frozen before the test tables are read.
"""
import sys, zlib
import numpy as np
import edgerules as R
LOGF = ('aln_bp', 'n_records', 'exon_bp', 'rank_max', 'degree_min', 'degree_max', 'neutral_bp')
L2, FOLDS, MINP, MINT = 1.0, 5, 0.75, 20
TRUE = set()


def prep(feats, X):
    X = X.copy()
    for j, f in enumerate(feats):
        if f in LOGF: X[:, j] = np.log1p(X[:, j])
    return X


def fit(X, y):
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Z = np.hstack([np.ones((len(X), 1)), (X - mu) / sd]); w = np.zeros(Z.shape[1]); reg = L2 * np.eye(Z.shape[1]); reg[0, 0] = 0
    for _ in range(30):
        p = 1 / (1 + np.exp(-np.clip(Z @ w, -30, 30)))
        g = Z.T @ (p - y) + reg @ w; H = (Z * (p * (1 - p))[:, None]).T @ Z + reg
        w -= np.linalg.solve(H, g)
    return mu, sd, w


def predict(m, X):
    mu, sd, w = m
    return np.hstack([np.ones((len(X), 1)), (X - mu) / sd]) @ w


def threshold(s, y):
    """lowest score with precision >= MINP among the top-ranked; (thr, true, false) or None"""
    o = np.argsort(-s); t = np.cumsum(y[o]); n = np.arange(1, len(o) + 1)
    ok = np.where((t / n >= MINP) & (t >= MINT))[0]
    if len(ok) == 0: return None
    i = ok[-1]
    return float(s[o][i]), int(t[i]), int(n[i] - t[i])


def selftest():
    rng = np.random.default_rng(1); X = rng.normal(size=(4000, 3)); y = (X[:, 0] + 0.3 * rng.normal(size=4000) > 2).astype(float)
    m = fit(X, y); s = predict(m, X); th = threshold(s, y.astype(bool))
    assert th and th[1] >= 20 and th[1] / (th[1] + th[2]) >= 0.75 and abs(m[2][1]) > abs(m[2][2]) * 3, (th, m[2])
    assert threshold(np.array([3., 2., 1.]), np.array([False, True, True])) is None
    print('selftest ok')


def pool(paths):
    Xs, ys, ks = [], [], []; feats = None
    for p in paths:
        feats, X, y, rows = R.load(p); Xs.append(prep(feats, X)); ys.append(y); ks += [f'{p}|{r[0]}|{r[1]}' for r in rows]
        TRUE.update(f'{p}|{r[0]}|{r[1]}' for r in rows if r[4] == '1')
    return feats, np.vstack(Xs), np.concatenate(ys), ks


def main(a):
    global MINP
    opt = lambda k: a[a.index(k) + 1]
    R.HOMOLOG = '--homolog' in a
    if '--minp' in a: MINP = float(opt('--minp'))
    bar = (0.85, 50) if R.HOMOLOG else (0.60, 20)
    feats, X, y, ks = pool(opt('--train').split(','))
    fold = np.array([zlib.crc32(k.encode()) % FOLDS for k in ks]); oof = np.zeros(len(y))
    for f in range(FOLDS):
        oof[fold == f] = predict(fit(X[fold != f], y[fold != f].astype(float)), X[fold == f])
    th = threshold(oof, y)
    print(f'train pool {len(y)} pairs, {int(y.sum())} positives; out-of-fold threshold (score, correct, wrong):', th)
    if th is not None and R.HOMOLOG:
        young = sum(1 for k, s_ in zip(ks, oof) if s_ >= th[0] and k in TRUE)
        print(f'  out-of-fold picks: truth pairs {young} -> go/no-go {"GO" if th[1] >= 100 and young >= 30 else "STOP"}')
        if th[1] < 100 or young < 30: return
    if th is None:
        print('STOP: no out-of-fold threshold reaches precision 0.75 with 20 true'); 
        o = np.argsort(-oof); t = np.cumsum(y[o])
        for n in (20, 40, 80): 
            i = np.searchsorted(t, n); print(f'  out-of-fold at {n} true: false {i + 1 - n} precision {n / (i + 1):.3f}')
        return
    m = fit(X, y.astype(float)); thr = th[0]
    print('weights (standardized):', ' '.join(f'{f} {w:+.2f}' for f, w in zip(['intercept'] + feats, m[2])))
    tt = tf = ty = 0
    for p in opt('--test').split(','):
        f2, X2, y2, rows = R.load(p); s = predict(m, prep(f2, X2)); k = s >= thr
        t, f = int((k & y2).sum()), int((k & ~y2).sum()); tt += t; tf += f
        young = sum(1 for r, kk in zip(rows, k) if kk and r[4] == '1'); ty += young
        print(f'  test {p}: population {len(y2)} positives {int(y2.sum())}; score >= {thr:.3f}: correct {t} (truth pairs {young}) wrong {f} precision {t / max(1, t + f):.3f}')
        if '--out-prefix' in a:
            with open(opt('--out-prefix') + '.' + p.split('/')[-1], 'w') as w:
                h = open(p).readline().rstrip('\n').split('\t'); c = [h.index(x) for x in ('locus_a', 'locus_b', 'name_a', 'name_b', 'truth')] + [h.index(x) for x in ('related', 'rel_kind') if x in h]
                w.write('locus_a\tlocus_b\tname_a\tname_b\ttruth' + ('\trelated\trel_kind' if len(c) > 5 else '') + '\tscore\n')
                w.writelines('\t'.join(r[i] for i in c) + f'\t{sc:.4f}\n' for r, kk, sc in zip(rows, k, s) if kk)
    print(f'pooled test: correct {tt} (truth pairs {ty}) wrong {tf} precision {tt / max(1, tt + tf):.3f} -> {"PASS" if tt >= bar[1] and tt / max(1, tt + tf) >= bar[0] and (not R.HOMOLOG or ty >= 20) else "FAIL"}')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(sys.argv[1:])
