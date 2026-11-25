#!/usr/bin/env python3
"""PREREG r1677: rules table (Part 1) and cross-sample score (Part 2) from chain_rules_feat.py tables.

  chain_rules.py NAME=feat.tsv NAME=feat.tsv [test:NAME=feat.tsv]

The first two samples are fitted on each other; a `test:` sample is scored once by the pooled fit of the first two."""
import sys, csv, math, bisect
import numpy as np

def load(path):
    R = list(csv.DictReader(open(path), delimiter='\t')); d = {k: np.array([r[k] for r in R]) for k in ('chrom', 'key') if k in R[0]}
    for k in R[0]:
        if k not in ('chrom', 'strand', 'key'): d[k] = np.array([float(r[k]) for r in R])
    return d

LOGS = 'full partial cont5 cont3 min_isup locus_chains n_introns length'.split()
OLD = ['l_reads', 'l_n_introns', 'l_length', 'l_locus_frac', 'l_min_isup', 'l_min_isup_rel', 'wobble', 'retention', 'a_rich']  # r1619's set
def design(d, reads_key, keep, cols=None, raw=False):
    f = {'l_reads': np.log(np.maximum(d[reads_key], 1))}
    for k in LOGS: f['l_' + k] = np.log1p(d[k])
    f['l_locus_frac'] = np.log(d['locus_frac'] + 1e-3); f['l_min_isup_rel'] = np.log(d['min_isup_rel'] + 1e-3)
    f['r_cont5'] = np.log((d['cont5'] + 1) / (d['full'] + 1)); f['r_cont3'] = np.log((d['cont3'] + 1) / (d['full'] + 1))
    f['r_partial'] = np.log((d['partial'] + 1) / (d['full'] + 1)); f['full0'] = (d['full'] == 0).astype(float)
    for k in 'wobble retention end3_conc end5_conc pas a_rich mq0_share de'.split(): f[k] = d[k]
    names = cols or sorted(f); X = np.array([f[k][keep] for k in names]).T
    if raw: return X, names  # register 1813: the features before the within-sample standardisation
    sd = X.std(0); sd[sd == 0] = 1
    return np.hstack([(X - X.mean(0)) / sd, np.ones((len(X), 1))]), names  # standardised within the sample

def design_decile(X, noncanon, min_intron):
    """Register 1769 (C1): each standardised column of design() as nine indicators (above its within-sample deciles) plus the
    column itself, then the count of non-GT-AG introns (capped at 2) and the standardised log of the shortest intron."""
    cols = []
    for j in range(X.shape[1] - 1):
        for t in np.quantile(X[:, j], np.arange(0.1, 1, 0.1)): cols.append((X[:, j] > t).astype(float))
        cols.append(X[:, j])
    mi = np.log(min_intron)
    cols += [np.minimum(noncanon, 2), (mi - mi.mean()) / (mi.std() or 1), np.ones(len(X))]
    return np.array(cols).T

def design_decile_ties(X, noncanon, min_intron):
    """Register 1813: design_decile with tied values split: a chain whose value v is shared gets, for decile k, the share of the
    tie's rank interval [lo, hi) (as fractions of the sample) above k/10; untied chains get 0 or 1 as before. Continuous in the
    chain set, where `X > quantile` flips a whole tie at once."""
    cols = []; n = len(X); q = np.arange(0.1, 1, 0.1)
    for j in range(X.shape[1] - 1):
        v = X[:, j]; u, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
        hi = np.cumsum(cnt) / n; lo = hi - cnt / n
        for t in q: cols.append(np.clip((hi - t) / (hi - lo), 0, 1)[inv])
        cols.append(v)
    mi = np.log(min_intron)
    cols += [np.minimum(noncanon, 2), (mi - mi.mean()) / (mi.std() or 1), np.ones(n)]
    return np.array(cols).T

def design_frozen(Xraw, noncanon, min_intron, m):
    """Register 1813, second attempt: design_decile with nothing taken from the sample. Columns are standardised by the fit
    samples' pooled mean and sd (m['mu'], m['sd']) and cut at their pooled deciles (m['thr']); the shortest intron likewise
    (m['mi']). A chain's score then depends on its own features only, so it cannot move with the rest of the chain set."""
    X = (Xraw - np.array(m['mu'])) / np.array(m['sd']); cols = []
    for j in range(X.shape[1]):
        for t in m['thr'][j]: cols.append((X[:, j] > t).astype(float))
        cols.append(X[:, j])
    cols += [np.minimum(noncanon, 2), (np.log(min_intron) - m['mi'][0]) / m['mi'][1], np.ones(len(X))]
    return np.array(cols).T

def fit(X, y):
    w = np.zeros(X.shape[1])
    for _ in range(40):
        p = 1 / (1 + np.exp(-np.clip(X @ w, -30, 30))); H = (X * (p * (1 - p))[:, None]).T @ X + 1e-3 * np.eye(len(w))
        w += np.linalg.solve(H, X.T @ (y - p) - 1e-3 * w)
    return w

def curve(score, y):  # ties in random order (seeded), so a tied block reads as its mean
    o = np.random.default_rng(0).permutation(len(y)); o = o[np.argsort(-score[o], kind='stable')]
    c = np.cumsum(y[o]); return c, c / np.arange(1, len(y) + 1), o
def prec_at(score, y, C):
    c, p, _ = curve(score, y); return float(p[min(len(y) - 1, bisect.bisect_left(c, C))])

def table(S):
    strata = [('2', 2, 2), ('3-5', 3, 5), ('6-20', 6, 20), ('21+', 21, 1e18)]
    cls = [('cont5 >= full', lambda d: d['cont5'] >= d['full']), ('cont3 >= full', lambda d: d['cont3'] >= d['full']),
           ('partial = 0', lambda d: d['partial'] == 0), ('partial >= 5 x full', lambda d: d['partial'] >= 5 * d['full']),
           ('min_isup_rel < 0.05', lambda d: d['min_isup_rel'] < 0.05), ('wobble', lambda d: d['wobble'] == 1),
           ('retention', lambda d: d['retention'] == 1), ('locus_frac < 0.05', lambda d: d['locus_frac'] < 0.05),
           ('locus_chains >= 10', lambda d: d['locus_chains'] >= 10), ('end3_conc < 0.5', lambda d: d['end3_conc'] < 0.5),
           ('pas = 0', lambda d: d['pas'] == 0), ('a_rich >= 12', lambda d: d['a_rich'] >= 12),
           ('mq0_share >= 0.5', lambda d: d['mq0_share'] >= 0.5), ('de > median', lambda d: d['de'] > np.median(d['de'][d['in_default'] == 1])),
           ('n_introns = 1', lambda d: d['n_introns'] == 1)]
    print('## Part 1: share correct, class yes / no, within read strata (default chains)')
    print('class\t' + '\t'.join(f'{n} {s[0]}' for n in S for s in strata) + '\t' + '\t'.join(f'{n} n_yes\t{n} pooled ratio' for n in S) + '\tstable')
    for name, fn in cls:
        cells = []; tails = []; sides = set(); ok = True
        for n, d in S.items():
            D = d['in_default'] == 1; yes = fn(d) & D; a = b = 0.0
            for _, lo, hi in strata:
                st = D & (d['default_reads'] >= lo) & (d['default_reads'] <= hi); y1 = st & yes; y0 = st & ~yes
                if y1.sum() < 50 or y0.sum() < 50: cells.append('-'); continue
                p1 = d['in_ref'][y1].mean(); p0 = d['in_ref'][y0].mean(); cells.append(f'{p1:.3f}/{p0:.3f}'); sides.add(p1 > p0)
                a += p1 * y1.sum(); b += p0 * y1.sum()  # observed vs expected correct among yes, at the no-class rate per stratum
            r = a / b if b else float('nan'); tails.append(f'{int(yes.sum())}\t{r:.2f}'); ok &= (r <= 0.67 or r >= 1.5)
        print(name + '\t' + '\t'.join(cells) + '\t' + '\t'.join(tails) + '\t' + ('STABLE' if ok and len(sides) == 1 else 'no'))

def report_A(tag, score, d, keep):
    y = d['in_ref'][keep]; rd = d['default_reads'][keep]; C0 = y.sum(); out = []
    for q in (0.95, 0.90, 0.85, 0.75, 0.50):
        out.append(f'{q:.2f}: {prec_at(score, y, q * C0):.3f} vs {prec_at(rd, y, q * C0):.3f} ({prec_at(score, y, q * C0) - prec_at(rd, y, q * C0):+.3f})')
    print(f'A  {tag}\t' + '\t'.join(out))
def report_B(tag, score, d):
    D = d['in_default'] == 1; y = d['in_ref']; n = int(D.sum()); top = np.zeros(len(y), bool); top[curve(score, y)[2][:n]] = True
    c0 = y[D].sum(); c1 = y[top].sum(); new = top & ~D; lost = D & ~top
    print(f'B  {tag}\tkept {n}\tcorrect {c1:.0f} vs default {c0:.0f} ({c1 / c0 - 1:+.4f})\tnew {new.sum()} / correct {y[new].sum():.0f} ({y[new].mean() if new.sum() else 0:.3f})\tlost {lost.sum()} / correct {y[lost].sum():.0f} ({y[lost].mean() if lost.sum() else 0:.3f})')

if __name__ == '__main__':
    args = [a.split('=') for a in sys.argv[1:]]; S = {n: load(p) for n, p in args if not n.startswith('test:')}
    T = {n[5:]: load(p) for n, p in args if n.startswith('test:')}
    for n, d in {**S, **T}.items():
        D = d['in_default'] == 1; print(f"# {n}: superset {len(D)} / correct {d['in_ref'].sum():.0f}; default {D.sum()} / correct {d['in_ref'][D].sum():.0f} ({d['in_ref'][D].mean():.3f})")
    table(S); a, b = list(S)
    print('\n## Part 2: score curve vs read-count curve at shares of the default correct chains (A), and equal-count swap on the superset (B)')
    for src, dst in ((a, b), (b, a)):
        for label, cols in (('all features', None), ('r1619 features', OLD)):
            Ds = S[src]['in_default'] == 1; Dd = S[dst]['in_default'] == 1
            Xs, names = design(S[src], 'default_reads', Ds, cols); w = fit(Xs, S[src]['in_ref'][Ds])
            report_A(f'fit {src} -> {dst}, {label}', design(S[dst], 'default_reads', Dd, cols)[0] @ w, S[dst], Dd)
            if cols is None: print('   coef', ' '.join(f'{k}={v:+.2f}' for k, v in sorted(zip(names, w), key=lambda x: -abs(x[1]))))
        allr = np.ones(len(S[src]['in_ref']), bool); Xs, names = design(S[src], 'reads', allr); w = fit(Xs, S[src]['in_ref'])
        allr = np.ones(len(S[dst]['in_ref']), bool); report_B(f'fit {src} -> {dst}', design(S[dst], 'reads', allr)[0] @ w, S[dst])
        print('   coef', ' '.join(f'{k}={v:+.2f}' for k, v in sorted(zip(names, w), key=lambda x: -abs(x[1]))))
    for n, d in S.items():  # within-sample ceiling, 5 folds by chromosome
        D = d['in_default'] == 1; X, _ = design(d, 'default_reads', D); y = d['in_ref'][D]
        fold = np.unique(d['chrom'][D], return_inverse=True)[1] % 5; s = np.zeros(len(y))
        for k in range(5): s[fold == k] = X[fold == k] @ fit(X[fold != k], y[fold != k])
        report_A(f'{n} 5-fold by chromosome', s, d, D)
    for n, d in T.items():  # frozen test: pooled fit of the two seen samples
        for kind, key in (('A', 'default_reads'), ('B', 'reads')):
            Xs = []; ys = []
            for m in S.values():
                k = (m['in_default'] == 1) if kind == 'A' else np.ones(len(m['in_ref']), bool); Xs.append(design(m, key, k)[0]); ys.append(m['in_ref'][k])
            w = fit(np.vstack(Xs), np.concatenate(ys)); k = (d['in_default'] == 1) if kind == 'A' else np.ones(len(d['in_ref']), bool)
            (report_A(f'FROZEN pooled -> {n}', design(d, key, k)[0] @ w, d, k) if kind == 'A' else report_B(f'FROZEN pooled -> {n}', design(d, key, k)[0] @ w, d))

COMP = str.maketrans('ACGT', 'TGCA')
def intron_extras(d, keep, fa_p):
    """per chain: introns that are not GT-AG on the model's strand, length of the shortest intron (genome only)"""
    import pysam
    fa = pysam.FastaFile(fa_p); nc = []; mi = []
    for c, key, s in zip(d['chrom'][keep], d['key'][keep], d['strand'][keep]):
        k = 0; m = 1 << 60
        for iv in key.split(','):
            a, b = map(int, iv.split('-')); m = min(m, b - a); dn = (fa.fetch(c, a, a + 2) + fa.fetch(c, b - 2, b)).upper()
            if s == '-': dn = dn.translate(COMP)[::-1]
            k += dn != 'GTAG'
        nc.append(k); mi.append(m)
    return np.array(nc, float), np.array(mi, float)
