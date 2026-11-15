#!/usr/bin/env python3
"""Row 1678: second stage on the relative tier's picks (docs/PREREG_young_links_2026-11-13.md).

  edgeyoung.py --train A.tsv,H.tsv --test OR.tsv,KB.tsv [--out-prefix P]     (tables from edgerelated.py)
Stage 1 = edgescore.py's homolog model fitted on the training tables, threshold THR1 (row 1676, frozen).
Stage 2 = the same logistic on the stage-1 picks, positives = truth pairs."""
import sys, zlib
import numpy as np
import edgerules as R, edgescore as S
THR1 = 0.05105613186294716


def tables(paths):
    out = []
    for p in paths:
        R.HOMOLOG = True; feats, X, y, rows = R.load(p)
        out.append((p, feats, S.prep(feats, X), y, np.array([r[4] == '1' for r in rows]), rows))
    return out


def main(a):
    opt = lambda k: a[a.index(k) + 1]
    tr, te = tables(opt('--train').split(',')), tables(opt('--test').split(','))
    m1 = S.fit(np.vstack([t[2] for t in tr]), np.concatenate([t[3] for t in tr]).astype(float))
    pick = lambda t: S.predict(m1, t[2]) >= THR1
    X = np.vstack([t[2][pick(t)] for t in tr]); y = np.concatenate([t[4][pick(t)] for t in tr])
    ks = [f'{t[0]}|{r[0]}|{r[1]}' for t in tr for r, k in zip(t[5], pick(t)) if k]
    fold = np.array([zlib.crc32(k.encode()) % S.FOLDS for k in ks]); oof = np.zeros(len(y))
    for f in range(S.FOLDS): oof[fold == f] = S.predict(S.fit(X[fold != f], y[fold != f].astype(float)), X[fold == f])
    th = S.threshold(oof, y)
    print(f'stage-1 picks on the training tables {len(y)}, truth pairs {int(y.sum())}; out-of-fold threshold (score, true, false): {th}')
    o = np.argsort(-oof); t = np.cumsum(y[o])
    for n in (20, 40, 80, 120):
        if n <= t[-1]: i = int(np.searchsorted(t, n)); print(f'  out-of-fold at {n} true: false {i + 1 - n} precision {n / (i + 1):.3f}')
    if th is None:
        print('STOP: no out-of-fold threshold reaches precision 0.75 with 20 true'); return
    m2 = S.fit(X, y.astype(float)); tt = tf = 0
    print('stage-2 weights:', ' '.join(f'{f} {w:+.2f}' for f, w in zip(['intercept'] + tr[0][1], m2[2])))
    for p, feats, X2, _, y2, rows in te:
        k1 = S.predict(m1, X2) >= THR1; k = k1 & (S.predict(m2, X2) >= th[0])
        t_, f_ = int((k & y2).sum()), int((k & ~y2).sum()); tt += t_; tf += f_
        print(f'  test {p}: stage-1 picks {int(k1.sum())} (truth {int((k1 & y2).sum())}); shortlist true {t_} false {f_} precision {t_ / max(1, t_ + f_):.3f}')
        if '--out-prefix' in a:
            h = open(p).readline().rstrip('\n').split('\t'); c = [h.index(x) for x in ('locus_a', 'locus_b', 'name_a', 'name_b', 'truth', 'related', 'rel_kind')]
            with open(opt('--out-prefix') + '.' + p.split('/')[-1], 'w') as w:
                w.write('locus_a\tlocus_b\tname_a\tname_b\ttruth\trelated\trel_kind\n'); w.writelines('\t'.join(r[i] for i in c) + '\n' for r, kk in zip(rows, k) if kk)
    print(f'pooled test: true {tt} false {tf} precision {tt / max(1, tt + tf):.3f} -> {"PASS" if tt >= 20 and tt / max(1, tt + tf) >= 0.60 else "FAIL"}')


if __name__ == '__main__':
    main(sys.argv[1:])
