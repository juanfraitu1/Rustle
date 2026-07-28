#!/usr/bin/env python3
"""Bipartite (optimal 1:1) matching of PREDICTED copies to GROUND-TRUTH members, then size comparison.

Why bipartite matching rather than "does it overlap?": overlap-based scoring is many-to-many — one predicted
copy can overlap two truth members (an OVER-MERGE) and one truth member can be covered by two predicted
copies (a SPLIT). Both are invisible to a per-member recall number. Forcing an optimal **1:1** assignment
(Hungarian / `scipy.optimize.linear_sum_assignment`, maximising total reciprocal overlap) makes the leftovers
explicit and lets predicted vs true SIZES be compared pair-by-pair:

  matched pair with size ratio >> 1  -> predicted copy is OVER-EXTENDED (over-merge signature)
  matched pair with size ratio << 1  -> predicted copy is TRUNCATED    (fragmentation signature)
  unmatched prediction              -> spurious / extra copy
  unmatched truth                   -> missed member

Usage: bipartite_size_match.py <copies.tsv> [truth.bed] [--label NAME]
"""
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment

copies_tsv = sys.argv[1]
truth_bed  = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else "bench/soto/80_fams.chr.bed"
label      = sys.argv[sys.argv.index("--label")+1] if "--label" in sys.argv else copies_tsv.split("/")[-1]

pred = []
for i, ln in enumerate(open(copies_tsv)):
    if i == 0: continue
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 6:
        pred.append((f[0], f[3], int(f[4]), int(f[5])))          # family, chrom, start, end
truth = []
for ln in open(truth_bed):
    c, s, e, name, *_ = ln.rstrip("\n").split("\t")
    truth.append((name.split("|")[1], name.split("|")[0], c, int(s), int(e)))

# reciprocal overlap = overlap / max(len_pred, len_truth): 1.0 only when the intervals coincide, so it
# rewards getting the SIZE right, not merely touching the right place.
P, T = len(pred), len(truth)
gain = np.zeros((P, T))
for i, (_pf, pc, ps, pe) in enumerate(pred):
    for j, (_tf, _tg, tc, ts, te) in enumerate(truth):
        if pc != tc: continue
        ov = min(pe, te) - max(ps, ts)
        if ov > 0:
            gain[i, j] = ov / max(pe - ps, te - ts, 1)

ri, cj = linear_sum_assignment(-gain)                              # maximise total reciprocal overlap
pairs = [(i, j) for i, j in zip(ri, cj) if gain[i, j] > 0]         # drop zero-overlap forced matches

matched_pred = {i for i, _ in pairs}
matched_truth = {j for _, j in pairs}
ratios, over, under, exact = [], [], [], 0
for i, j in pairs:
    plen = pred[i][3] - pred[i][2]
    tlen = truth[j][4] - truth[j][3]
    r = plen / max(tlen, 1)
    ratios.append(r)
    if r >= 2.0:   over.append((truth[j][1], plen, tlen, r))
    elif r <= 0.5: under.append((truth[j][1], plen, tlen, r))
    else:          exact += 1

ratios = np.array(ratios)
log2r = np.log2(np.clip(ratios, 1e-9, None))
print(f"=== bipartite size match: {label} ===")
print(f"  predicted copies : {P}")
print(f"  truth members    : {T}")
print(f"  MATCHED 1:1      : {len(pairs)}  ({100*len(pairs)/T:.1f}% of truth)")
print(f"  unmatched truth  : {T - len(matched_truth)}   (missed members)")
print(f"  unmatched pred   : {P - len(matched_pred)}   (spurious / extra copies)")
print()
print(f"  size ratio (predicted / true), over matched pairs:")
print(f"    median            {np.median(ratios):.2f}")
print(f"    median log2 bias  {np.median(log2r):+.2f}   ({'over' if np.median(log2r)>0 else 'under'}-prediction)")
print(f"    IQR               {np.percentile(ratios,25):.2f} - {np.percentile(ratios,75):.2f}")
print(f"    within 2x (0.5-2) {exact}  ({100*exact/max(len(pairs),1):.0f}% of matched)")
print(f"    OVER-EXTENDED >=2x  {len(over)}   (over-merge signature)")
print(f"    TRUNCATED    <=0.5x {len(under)}  (fragmentation signature)")
if over:
    print("\n  worst over-extended (gene, pred_bp, true_bp, ratio):")
    for g, p_, t_, r in sorted(over, key=lambda x: -x[3])[:8]:
        print(f"    {g:14s} {p_:>9} {t_:>9} {r:>8.1f}x")
if under:
    print("\n  worst truncated:")
    for g, p_, t_, r in sorted(under, key=lambda x: x[3])[:8]:
        print(f"    {g:14s} {p_:>9} {t_:>9} {r:>8.2f}x")
