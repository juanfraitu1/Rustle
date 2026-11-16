#!/usr/bin/env python3
"""Addendum 3 of PREREG_ptr_ppa_frozen_rules_2026-11-15: gate and rescoring with sister-species confirmed chains.

  sister_score.py NAME=feat.tsv:sister.tsv:refseq.gtf[:dropped_keys.txt] ... [fit:NAME=feat.tsv ...]   |   --selftest
feat.tsv = chain_rules_feat.py table with the key column and the real 1-read-floor GTF; sister.tsv = xspecies_chains.py;
dropped_keys.txt = chrom<TAB>key of the chains a stage removed (shiftdrop). fit: samples = the frozen score's fit set."""
import sys, re, csv, collections
import numpy as np
from chain_rules import load, design, fit, prec_at

NB = ((1, 1), (2, 3), (4, 7), (8, 10**9)); RB = ((1, 1), (2, 2), (3, 5), (6, 20), (21, 10**9))
bin_ = lambda x, B: next(i for i, (a, b) in enumerate(B) if a <= x <= b)

def std_ratio(cls, yard, conf, cell, floor=20):
    """class rate, yardstick rate standardised to the class's cells, share of the class in cells the yardstick fills."""
    num = den = cov = 0
    for c in np.unique(cell[cls]):
        w = (cls & (cell == c)).sum(); y = yard & (cell == c)
        if y.sum() >= floor: num += w * conf[y].mean(); den += w
    cr = conf[cls].mean() if cls.sum() else 0.0; yr = num / den if den else 0.0
    return cr, yr, den / max(1, cls.sum())

def subchain(keys, extra):
    """bool per key: a contiguous sub-chain of a longer chain among keys + extra (same contig)."""
    U = {}
    for c, k in list(keys) + list(extra): U[(c, k)] = tuple(k.split(','))
    by = collections.defaultdict(list)
    for (c, k), t in U.items():
        for i in set(t): by[(c, i)].append(t)
    out = []
    for c, k in keys:
        t = U[(c, k)]; n = len(t); hit = False
        for b in by[(c, t[0])]:
            if len(b) > n and any(b[i:i + n] == t for i in range(len(b) - n + 1)): hit = True; break
        out.append(hit)
    return np.array(out)

def ref_keys(p):
    ex = collections.defaultdict(list)
    for l in open(p):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    out = set()
    for (c, _), e in ex.items():
        if len(e) > 1: e.sort(); out.add((c, ','.join(f'{a[1]}-{b[0] - 1}' for a, b in zip(e, e[1:]))))
    return out

def sample(name, spec):
    fp, sp, rp, *dp = spec.split(':'); d = load(fp); S = {(r['chrom'], r['key']): (int(r['lifted']), int(r['sister_reads'])) for r in csv.DictReader(open(sp), delimiter='\t')}
    K = list(zip(d['chrom'], d['key'])); miss = sum(k not in S for k in K); assert miss == 0, f'{name}: {miss} feature rows have no sister row (key convention?)'
    lift = np.array([S[k][0] for k in K]) == 1; conf = np.array([S[k][1] for k in K]) >= 1; y = d['in_ref'] == 1; D = d['in_default'] == 1; rd = d['default_reads']
    wide = y | conf; sub = subchain(K, ref_keys(rp)); ns = ~sub
    cell = np.array([bin_(n, NB) * 10 + bin_(max(r, 1), RB) for n, r in zip(d['n_introns'], np.where(D, rd, d['reads']))])
    c5 = D & (d['cont5'] >= d['full']) & (d['full'] >= 1)
    print(f"\n# {name}: default chains {D.sum()}, lift whole {lift[D].mean():.3f}, RefSeq precision {y[D].mean():.3f}, RefSeq or sister {wide[D].mean():.3f}; not a sub-chain {(D & ns).sum()}: {y[D & ns].mean():.3f} / {wide[D & ns].mean():.3f}")
    print(f"  sister-confirmed among RefSeq-exact default chains that lift {conf[D & y & lift].mean():.3f}; among RefSeq-absent {conf[D & ~y & lift].mean():.3f}")
    gate = {}
    for pop, m in (('all', np.ones(len(y), bool)), ('not a sub-chain', ns)):
        cr, yr, cov = std_ratio(c5 & lift & m, D & y & lift & m, conf, cell); gate[pop] = cr / yr if yr else 0
        print(f"  1 GATE [{pop}] cont5 class lifted {(c5 & lift & m).sum()} (lift {lift[c5 & m].mean():.3f}): confirmed {cr:.3f} vs RefSeq-exact standardised {yr:.3f} = ratio {gate[pop]:.2f} (cells covered {cov:.2f}) -> {'INVALID for end questions' if gate[pop] >= 0.5 else 'valid'}")
    for pop, m in (('all', np.ones(len(y), bool)), ('not a sub-chain', ns)):
        Dm = D & m; drop = c5 & m; k = Dm & ~drop; out = []
        for lab, t in (('RefSeq', y), ('wide', wide)):
            c = t[k].sum(); p = c / k.sum(); ref = prec_at(rd[Dm], t[Dm].astype(float), c); out.append(f"{lab}: dropped {t[drop].mean():.3f} correct, kept {p:.3f} vs reads curve {ref:.3f} ({p - ref:+.3f})")
        print(f"  4 trunc5 [{pop}] drops {drop.sum()} of {Dm.sum()}; " + '; '.join(out))
    one = ~D & (d['reads'] == 1); two = D & (rd == 2)
    for pop, m in (('all', np.ones(len(y), bool)), ('not a sub-chain', ns)):
        nb = np.array([bin_(n, NB) for n in d['n_introns']]); o = []
        for lab, t in (('RefSeq', y), ('wide', wide)):
            cr, yr, cov = std_ratio(one & m, two & m, t, nb); o.append(f"{lab} 1-read {cr:.3f} vs 2-read standardised {yr:.3f} (x{cr / yr if yr else 0:.2f})")
        print(f"  5 1-read class [{pop}] {(one & m).sum()} chains vs {(two & m).sum()} default 2-read chains; " + '; '.join(o))
    if dp:
        X = {tuple(l.rstrip('\n').split('\t')) for l in open(dp[0])}; dr = np.array([k in X for k in K]) & D
        print(f"  6 shiftdrop dropped {dr.sum()}: in RefSeq {y[dr].mean():.3f}, sister-confirmed {conf[dr & lift].mean():.3f} of {(dr & lift).sum()} lifted (default chains overall {conf[D & lift].mean():.3f})")
    return d, D, y, wide, rd

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        k = [('c', '1-2,3-4,5-6'), ('c', '3-4,5-6'), ('c', '1-2,5-6'), ('c', '3-4')]
        assert subchain(k, [('c', '0-1,1-2,3-4,5-6')]).tolist() == [True, True, False, True]
        cls = np.array([1, 1, 0, 0, 0, 0], bool); yard = ~cls; conf = np.array([1, 0, 1, 1, 0, 0], bool); cell = np.array([0, 1, 0, 0, 1, 1])
        assert std_ratio(cls, yard, conf, cell, floor=1) == (0.5, 0.5, 1.0)
        print('selftest ok'); sys.exit()
    fits = [load(a.split('=')[1]) for a in sys.argv[1:] if a.startswith('fit:')]; R = {}
    for a in sys.argv[1:]:
        if not a.startswith('fit:'): n, spec = a.split('='); R[n] = sample(n, spec)
    if fits:
        Xs = np.vstack([design(m, 'default_reads', m['in_default'] == 1)[0] for m in fits]); ys = np.concatenate([m['in_ref'][m['in_default'] == 1] for m in fits]); w = fit(Xs, ys)
        print('\n# 6 frozen logistic score (fit set pooled), score curve minus read-count curve at 0.90 of the correct default chains')
        for n, (d, D, y, wide, rd) in R.items():
            s = design(d, 'default_reads', D)[0] @ w; o = []
            for lab, t in (('RefSeq', y), ('wide', wide)):
                t = t[D].astype(float); C = 0.9 * t.sum(); a = prec_at(s, t, C); b = prec_at(rd[D], t, C); o.append(f"{lab} {a:.3f} vs {b:.3f} ({a - b:+.3f})")
            print(f"  {n}: " + '; '.join(o))
