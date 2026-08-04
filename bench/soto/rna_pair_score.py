#!/usr/bin/env python3
"""Pair-level partition score for an RNA `copies.tsv` against Soto's family assignment.

Usage: rna_pair_score.py <copies.tsv> [chr1,chr15]

Two properties this scorer is built around, both learned the hard way (see
bench/soto/EDGE_RULE_DETECTION_COUPLING.txt):

1. THE UNIVERSE IS FIXED. Every truth member is scored, whether or not the run placed it in a
   family. A member the run never emitted is given its own unique singleton label, so it
   contributes zero true pairs instead of vanishing from the denominator. Scoring only the
   members both sides placed lets a change that DELETES members look like a precision win.

2. DETECTION AND COPY COUNT ARE REPORTED ALONGSIDE PRECISION. Because a family must span >= 2
   distinct loci to be emitted at all, tightening an edge rule removes whole loci rather than
   only wrong edges. A precision number quoted without these two is not interpretable.
"""
import sys
from collections import defaultdict
from itertools import combinations

COP = sys.argv[1]
CHROMS = set(sys.argv[2].split(",")) if len(sys.argv) > 2 else None
BED = "bench/soto/80_fams.gene_preferred.bed"

# ---- truth: one entry per gene, labelled by Soto family id -------------------------------------
truth = []  # (gene, famid, chrom, start, end)
for ln in open(BED):
    c, s, e, name, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and c not in CHROMS:
        continue
    gene, fam = name.split("|")[0], name.split("|")[1]
    truth.append((gene, fam, c, int(s), int(e)))

# ---- prediction: copies grouped by emitted family -----------------------------------------------
by_chrom = defaultdict(list)
n_copies = 0
fams = set()
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 6:
        continue
    fid, c, s, e = f[0], f[3], int(f[4]), int(f[5])
    if CHROMS and c not in CHROMS:
        continue
    by_chrom[c].append((s, e, fid))
    fams.add(fid)
    n_copies += 1

# ---- label each truth member by the predicted family covering it --------------------------------
pred = {}
detected = 0
for idx, (gene, fam, c, s, e) in enumerate(truth):
    # Best-overlap, not first-hit: several emitted copies can overlap one truth member (isoforms,
    # nested loci), and "whichever came first in the file" would make the score depend on row order.
    ov = [(min(e, pe) - max(s, ps), fid) for (ps, pe, fid) in by_chrom.get(c, []) if s < pe and e > ps]
    hit = max(ov)[1] if ov else None
    if hit is None:
        pred[gene] = f"__unplaced_{idx}"  # unique singleton: keeps the universe fixed
    else:
        pred[gene] = hit
        detected += 1

# ---- pair-level precision / recall ---------------------------------------------------------------
genes = [t[0] for t in truth]
tfam = {t[0]: t[1] for t in truth}
tp = fp = fn = 0
for a, b in combinations(genes, 2):
    same_t = tfam[a] == tfam[b]
    same_p = pred[a] == pred[b]
    if same_t and same_p:
        tp += 1
    elif same_p:
        fp += 1
    elif same_t:
        fn += 1
P = tp / (tp + fp) if tp + fp else 0.0
R = tp / (tp + fn) if tp + fn else 0.0
F1 = 2 * P * R / (P + R) if P + R else 0.0

print(f"members={len(truth)}  true_pairs={tp+fn}  detected={detected}  copies={n_copies}  "
      f"families={len(fams)}")
print(f"P={P:.3f}  R={R:.3f}  F1={F1:.3f}   (TP={tp} FP={fp} FN={fn})")
