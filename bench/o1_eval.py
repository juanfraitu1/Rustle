#!/usr/bin/env python3
"""O1 evaluation (§6ft): sensitivity, specificity and bipartite truth coverage of a catalog against a truth set of
loci. Matching = optimal 1:1 (Hungarian) on reciprocal overlap ov / max(len_pred, len_true), a pair counted when
either interval is ≥ 50 % covered (bench/soto/bipartite_size_match.py's rule). Per matched pair: truth coverage
= ov / len_true (how much of the ground-truth member the rediscovered unit covers), size ratio, in-band (0.5–2×).
  sensitivity = matched truth / truth ;  specificity = matched predictions / predictions (1 − spurious fraction)
usage: o1_eval.py <units.tsv or copies.tsv> <truth.bed: chrom start end name> [--family F] [--extent] [--out prefix]
"""
import sys, csv, numpy as np
from scipy.optimize import linear_sum_assignment
args = sys.argv[1:]; units_p, truth_p = args[0], args[1]
fam = args[args.index('--family') + 1] if '--family' in args else None
use_extent = '--extent' in args; out = args[args.index('--out') + 1] if '--out' in args else None
pred = []
for r in csv.DictReader(open(units_p), delimiter='\t'):
    if fam and r['family_id'] != fam: continue
    s, e = (int(r['locus_start']), int(r['locus_end'])) if use_extent and r.get('locus_start', 'NA') not in ('NA', '') else (int(r['start']), int(r['end']))
    pred.append((r['family_id'], r.get('copy_idx', '?'), r['chrom'], s, e, r.get('member_status', 'NA'), r.get('n_reads', '?')))
truth = []
for l in open(truth_p):
    f = l.rstrip('\n').split('\t')
    if len(f) < 3 or f[0].startswith('#'): continue
    truth.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else f'{f[0]}:{f[1]}'))
P, T = len(pred), len(truth)
W = np.zeros((P, T)); cov = np.zeros((P, T)); ok = np.zeros((P, T), dtype=bool)
for i, (_, _, c, s, e, _, _) in enumerate(pred):
    for j, (tc, ts, te, _) in enumerate(truth):
        if tc != c: continue
        ov = max(0, min(e, te) - max(s, ts))
        if ov <= 0: continue
        W[i, j] = ov / max(e - s, te - ts); cov[i, j] = ov / (te - ts)
        ok[i, j] = ov >= 0.5 * (e - s) or ov >= 0.5 * (te - ts)
rows, cols = linear_sum_assignment(-W)
pairs = [(i, j) for i, j in zip(rows, cols) if W[i, j] > 0 and ok[i, j]]
mt = {j for _, j in pairs}; mp = {i for i, _ in pairs}
NONMEMBER = ('dropped', 'readthrough', 'partner', 'noncoding')  # candidates and derived objects, never in the specificity denominator
members = [i for i in range(P) if pred[i][5] not in NONMEMBER]; cands = [i for i in range(P) if pred[i][5] in NONMEMBER]
sens = len(mt) / max(T, 1); spec = len([i for i in members if i in mp]) / max(len(members), 1)
ratios = [(pred[i][4] - pred[i][3]) / (truth[j][2] - truth[j][1]) for i, j in pairs]
covs = [cov[i, j] for i, j in pairs]
inband = sum(1 for r in ratios if 0.5 <= r <= 2)
print(f"O1 eval: predictions {P} | truth {T} | matched 1:1 {len(pairs)}")
print(f"  sensitivity (truth rediscovered) {len(mt)}/{T} = {sens:.3f} | specificity (family MEMBERS that are truth) {len([i for i in members if i in mp])}/{len(members)} = {spec:.3f} | candidates (dropped / readthrough / noncoding) {len(cands)}, of which truth {len([i for i in cands if i in mp])}")
print(f"  truth coverage by the matched unit: median {np.median(covs):.2f}, ≥ 0.5: {sum(1 for c in covs if c >= 0.5)}/{len(covs)}, ≥ 0.9: {sum(1 for c in covs if c >= 0.9)}/{len(covs)}")
print(f"  size ratio pred/true: median {np.median(ratios):.2f}, in-band 0.5–2×: {inband}/{len(ratios)} = {inband/max(1,len(ratios)):.2f}, truncated ≤ 0.5×: {sum(1 for r in ratios if r <= 0.5)}, over-extended ≥ 2×: {sum(1 for r in ratios if r >= 2)}")
miss = [truth[j] for j in range(T) if j not in mt]; extra = [pred[i] for i in range(P) if i not in mp]
print(f"  missed truth ({len(miss)}): " + '; '.join(f"{m[3]} {m[0]}:{m[1]}-{m[2]}" for m in miss[:8]))
print(f"  unmatched predictions ({len(extra)}): " + '; '.join(f"{p[0]}:{p[1]} {p[5]} reads={p[6]}" for p in extra[:8]))
if out:
    with open(out + '.pairs.tsv', 'w') as o:
        o.write('family\tcopy_idx\tmember_status\tn_reads\tpred\ttruth_name\ttruth\ttruth_coverage\tsize_ratio\tin_band\n')
        for (i, j), r, c in zip(pairs, ratios, covs):
            p, t = pred[i], truth[j]; o.write(f"{p[0]}\t{p[1]}\t{p[5]}\t{p[6]}\t{p[2]}:{p[3]}-{p[4]}\t{t[3]}\t{t[0]}:{t[1]}-{t[2]}\t{c:.3f}\t{r:.3f}\t{int(0.5 <= r <= 2)}\n")
