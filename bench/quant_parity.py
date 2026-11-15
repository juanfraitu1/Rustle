#!/usr/bin/env python3
"""PREREG row 1697: per-chain quantification parity with StringTie.  quant_parity.py OURS.gtf STRINGTIE.gtf[.gz]"""
import sys, re, gzip, collections, numpy as np
from scipy.stats import spearmanr, pearsonr

def load(path, keys):
    ex = collections.defaultdict(list); val = {}
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        t = (f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))
        if f[2] == 'exon': ex[t].append((int(f[3]), int(f[4])))
        else: val[t] = [float(m.group(1)) if (m := re.search(r'(?<![A-Za-z_])' + k + r' "([^"]+)"', f[8])) else 0.0 for k in keys]
    out = collections.defaultdict(lambda: np.zeros(len(keys)))
    for t, e in ex.items():
        if len(e) > 1 and t in val:
            e.sort(); out[(t[0], t[1], tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))] += val[t]
    return out

if __name__ == '__main__':
    a = load(sys.argv[1], ['reads', 'TPM']); b = load(sys.argv[2], ['cov', 'TPM']); sh = sorted(set(a) & set(b))
    for name, i, j in (('our reads vs StringTie cov', 0, 0), ('our TPM vs StringTie TPM', 1, 1)):
        x = np.array([a[k][i] for k in sh]); y = np.array([b[k][j] for k in sh])
        print(f'{name}\tshared chains {len(sh)} (ours {len(a)}, StringTie {len(b)})\tSpearman {spearmanr(x, y)[0]:.3f}\tPearson log1p {pearsonr(np.log1p(x), np.log1p(y))[0]:.3f}\tmedian ratio {np.median(x[y > 0] / y[y > 0]):.2f}')
