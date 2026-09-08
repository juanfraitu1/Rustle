#!/usr/bin/env python3
"""Score one ablation arm (PREREG_annotation_ablation_2026-09-07, md5 65efc5b2).
M1 fragmentation (selection-free): clusters holding >=1 truth locus.
M2 sensitivity: truth loci recovered by the dominant cluster / by any cluster.
M3 precision: dominant cluster's members that are truth / its members.
M4 bipartite (Hungarian 1:1 on the CORE HULL): median truth coverage, in-band 0.5-2x.
A truth locus is "held" by a unit when they overlap by >= 50% of either interval (the o1_eval rule).
usage: o1_ablation_score.py <arm.units.tsv> <truth.bed> [--baseline <A0.units.tsv>] [--json out.json]
"""
import sys, csv, json
import numpy as np
from scipy.optimize import linear_sum_assignment
a = sys.argv[1:]
units_p, truth_p = a[0], a[1]
base_p = a[a.index('--baseline') + 1] if '--baseline' in a else None
out_j = a[a.index('--json') + 1] if '--json' in a else None
truth = []
for l in open(truth_p):
    f = l.rstrip('\n').split('\t')
    if len(f) >= 3 and not f[0].startswith('#'):
        truth.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else f'{f[0]}:{f[1]}'))
rows = list(csv.DictReader(open(units_p), delimiter='\t'))
def ov(a1, b1, a2, b2): return max(0, min(b1, b2) - max(a1, a2))
def hit(r, t):
    if r['chrom'] != t[0]: return False
    s, e = int(r['start']), int(r['end'])
    o = ov(s, e, t[1], t[2])
    return o > 0 and (o >= 0.5 * (e - s) or o >= 0.5 * (t[2] - t[1]))
members = [r for r in rows if r['member_status'] != 'dropped']
by_fam = {}
for r in members: by_fam.setdefault(r['family_id'], []).append(r)
fam_truth = {}
for f, rs in by_fam.items():
    n = sum(1 for r in rs if any(hit(r, t) for t in truth))
    if n: fam_truth[f] = n
M1 = len(fam_truth)
dom = max(fam_truth, key=lambda f: (fam_truth[f], len(by_fam[f]))) if fam_truth else None
dom_rows = by_fam.get(dom, [])
held_any = {t[3] for t in truth if any(hit(r, t) for r in members)}
held_dom = {t[3] for t in truth if any(hit(r, t) for r in dom_rows)}
M3 = fam_truth.get(dom, 0) / max(len(dom_rows), 1)
# M4 on the core hull of the dominant cluster's members
pred = []
for r in dom_rows:
    h = r.get('core_hull', 'NA')
    if h in ('NA', ''): continue
    s, e = h.split('-'); pred.append((r['chrom'], int(s), int(e)))
W = np.zeros((len(pred), len(truth))); cov = np.zeros_like(W); ok = np.zeros(W.shape, dtype=bool)
for i, (c, s, e) in enumerate(pred):
    for j, t in enumerate(truth):
        if t[0] != c: continue
        o = ov(s, e, t[1], t[2])
        if o <= 0: continue
        W[i, j] = o / max(e - s, t[2] - t[1]); cov[i, j] = o / (t[2] - t[1])
        ok[i, j] = o >= 0.5 * (e - s) or o >= 0.5 * (t[2] - t[1])
pairs = []
if len(pred):
    ri, ci = linear_sum_assignment(-W)
    pairs = [(i, j) for i, j in zip(ri, ci) if W[i, j] > 0 and ok[i, j]]
ratios = [(pred[i][2] - pred[i][1]) / (truth[j][2] - truth[j][1]) for i, j in pairs]
covs = [cov[i, j] for i, j in pairs]
res = dict(arm=units_p, M1_clusters_with_truth=M1, dominant=dom,
           dominant_members=len(dom_rows), dominant_truth_members=fam_truth.get(dom, 0),
           M2_sens_dominant=f"{len(held_dom)}/{len(truth)}", M2_sens_any=f"{len(held_any)}/{len(truth)}",
           M3_precision=round(M3, 4),
           M4_cov_median=round(float(np.median(covs)), 3) if covs else None,
           M4_inband=f"{sum(1 for r in ratios if 0.5 <= r <= 2)}/{len(ratios)}",
           missed=sorted({t[3] for t in truth} - held_dom))
if base_p:
    b = {(r['chrom'], r['start'], r['end']) for r in csv.DictReader(open(base_p), delimiter='\t')
         if r['member_status'] != 'dropped'}
    # baseline dominant set, by the same rule
    brows = list(csv.DictReader(open(base_p), delimiter='\t')); bm = [r for r in brows if r['member_status'] != 'dropped']
    bfam = {}
    for r in bm: bfam.setdefault(r['family_id'], []).append(r)
    bft = {f: sum(1 for r in rs if any(hit(r, t) for t in truth)) for f, rs in bfam.items()}
    bft = {f: n for f, n in bft.items() if n}
    bdom = max(bft, key=lambda f: (bft[f], len(bfam[f])))
    bset = {(r['chrom'], int(r['start']), int(r['end'])) for r in bfam[bdom]}
    aset = {(r['chrom'], int(r['start']), int(r['end'])) for r in dom_rows}
    # compare by truth-locus identity, not coordinates (degraded arms move boundaries)
    def labels(rs):
        out = set()
        for r in rs:
            for t in truth:
                if hit(r, t): out.add(t[3])
        return out
    res['baseline_members'] = len(bfam[bdom]); res['baseline_truth_labels'] = len(labels(bfam[bdom]))
    res['labels_lost_vs_baseline'] = sorted(labels(bfam[bdom]) - labels(dom_rows))
    res['labels_gained_vs_baseline'] = sorted(labels(dom_rows) - labels(bfam[bdom]))
    res['exact_coord_match'] = (aset == bset)
for k, v in res.items(): print(f"  {k}: {v}")
if out_j: json.dump(res, open(out_j, 'w'), indent=1)
