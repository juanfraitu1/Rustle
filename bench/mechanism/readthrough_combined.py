#!/usr/bin/env python3
"""Combine per-junction BAM-tag medians (+/- the R, Q locus ratios) in one logistic model; train chr16, test
chr20 + gorilla. Needs <tag>.junctions.tsv and <tag>.junction_tags.tsv. usage: readthrough_combined.py OUTDIR"""
import sys, csv, math
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
O = sys.argv[1]
TAGS = ['de', 'NM_per_bp', 'AS_per_bp', 's2_over_s1', 'cm', 'mapq', 'softclip', 'min_anchor', 'local_err']


def load(t):
    J = {(r['start'], r['end']): r for r in csv.DictReader(open(f'{O}/{t}.junctions.tsv'), delimiter='\t')}
    X, y = [], []
    for r in csv.DictReader(open(f'{O}/{t}.junction_tags.tsv'), delimiter='\t'):
        j = J.get((r['start'], r['end']))
        if j is None or j['canonical'] != 'True':
            continue
        v = [float(r[k]) for k in TAGS]
        if any(math.isnan(x) for x in v):
            continue
        X.append(v + [math.log(int(j['reads'])), float(j['R']), float(j['Q'])]); y.append(r['label'] == 'RT')
    return np.array(X), np.array(y)


sets = {'T': list(range(10)), 'RQ': [10, 11], 'T+RQ': list(range(12))}
data = {t: load(t) for t in ('hsa16', 'hsa20', 'ggo44')}
for t, (X, y) in data.items():
    print(f'{t}: {y.sum()} RT junctions, {(~y).sum()} control junctions')
for name, cols in sets.items():
    X, y = data['hsa16']
    sc = StandardScaler().fit(X[:, cols])
    m = LogisticRegression(class_weight='balanced', max_iter=5000).fit(sc.transform(X[:, cols]), y)
    res = []
    for t in ('hsa16', 'hsa20', 'ggo44'):
        Xt, yt = data[t]
        p = m.predict_proba(sc.transform(Xt[:, cols]))[:, 1]
        a = roc_auc_score(yt, p)
        th = np.sort(p[yt])[::-1][int(0.5 * yt.sum()) - 1]          # FPR at recall .5
        res.append(f'{t} AUC {a:.3f} FPR@rec.5 {np.mean(p[~yt] >= th):.4f}')
    print(f'{name:5s} | ' + ' | '.join(res))
    if name == 'T':
        print('      coefs:', ', '.join(f'{k} {c:+.2f}' for k, c in zip(TAGS + ['log_reads'], m.coef_[0])))

# ---- depth control: is it the tags or just "readthroughs are rarely-read junctions"? ----
print('\n-- depth control --')
extra = {'count only': [9], 'tags, no count': list(range(9))}
for name, cols in extra.items():
    X, y = data['hsa16']
    sc = StandardScaler().fit(X[:, cols])
    m = LogisticRegression(class_weight='balanced', max_iter=5000).fit(sc.transform(X[:, cols]), y)
    print(f'{name:15s} | ' + ' | '.join(f'{t} AUC {roc_auc_score(data[t][1], m.predict_proba(sc.transform(data[t][0][:, cols]))[:, 1]):.3f}'
                                       for t in ('hsa16', 'hsa20', 'ggo44')))
# depth-matched: controls re-sampled to the RT read-count distribution (log2 bins), then tags-only AUC
rng = np.random.default_rng(1)
def matched(X, y):
    b = np.floor(X[:, 9] / math.log(2)).astype(int)
    keep = list(np.where(y)[0])
    for k in np.unique(b[y]):
        pool = np.where((~y) & (b == k))[0]
        need = 5 * int(((b == k) & y).sum())
        if len(pool):
            keep += list(rng.choice(pool, size=min(need, len(pool)), replace=False))
    keep = np.array(keep); return X[keep], y[keep]
cols = list(range(9))
Xd, yd = matched(*data['hsa16'])
sc = StandardScaler().fit(Xd[:, cols]); m = LogisticRegression(class_weight='balanced', max_iter=5000).fit(sc.transform(Xd[:, cols]), yd)
out = []
for t in ('hsa16', 'hsa20', 'ggo44'):
    Xm, ym = matched(*data[t])
    out.append(f'{t} AUC {roc_auc_score(ym, m.predict_proba(sc.transform(Xm[:, cols]))[:, 1]):.3f} (RT {ym.sum()}, ctrl {(~ym).sum()})')
print('tags, depth-matched | ' + ' | '.join(out))
cols = [10, 11]
sc = StandardScaler().fit(Xd[:, cols]); m = LogisticRegression(class_weight='balanced', max_iter=5000).fit(sc.transform(Xd[:, cols]), yd)
print('R+Q,  depth-matched | ' + ' | '.join(f"{t} AUC {roc_auc_score(ym, m.predict_proba(sc.transform(Xm[:, cols]))[:, 1]):.3f}"
                                             for t, (Xm, ym) in ((t, matched(*data[t])) for t in ('hsa16', 'hsa20', 'ggo44'))))
