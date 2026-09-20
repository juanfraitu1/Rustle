#!/usr/bin/env python3
"""How much of DSFAM42's 'tied' (n_decisive=0, unassignable) pile is the TRUE information-theoretic
floor vs improvable coverage-face? Discriminator = read completeness. An ASSIGNED read (n_decisive>=1)
is, by definition, long enough to reach a distinguishing PSV. A TIED read reaches none — either
because (a) it is too SHORT/partial (a full-length capture of that molecule would likely hit a PSV =
coverage-face, IMPROVABLE by longer reads), or (b) it is FULL-LENGTH yet still spans no PSV (the copies
are identical over the whole molecule = TRUE FLOOR, no read length or EM/ML helps).

So compare the exonic-length / exon-count distributions of TIED vs ASSIGNED reads."""
import sys, pysam
import numpy as np

BAM = "/home/juanfra/winloci_scratch/GGO.bam"
REGION = ("NC_073247.2", 59690000, 59790000)
ASN = "/tmp/dsfam42.assignments.tsv"

status = {}
for line in open(ASN):
    c = line.rstrip("\n").split("\t")
    if c[0] == "read_name":
        continue
    status[c[0]] = c[3]   # assigned / ambiguous / tied

bam = pysam.AlignmentFile(BAM, "rb")
feat = {}  # read -> (exonic_len, n_exons, genomic_span)
for a in bam.fetch(*REGION):
    if a.is_secondary or a.is_supplementary or a.is_unmapped:
        continue
    if a.query_name not in status or a.query_name in feat:
        continue
    exonic = sum(l for op, l in (a.cigartuples or []) if op in (0, 7, 8))   # M/=/X
    n_introns = sum(1 for op, l in (a.cigartuples or []) if op == 3)
    feat[a.query_name] = (exonic, n_introns + 1, a.reference_end - a.reference_start)

def dist(label):
    vals = [feat[r] for r, s in status.items() if s == label and r in feat]
    if not vals:
        return None
    ex = np.array([v[0] for v in vals]); nex = np.array([v[1] for v in vals])
    return dict(n=len(vals), med_exonic=int(np.median(ex)), mean_exonic=int(ex.mean()),
               med_exons=int(np.median(nex)), single_exon=int((nex == 1).sum()),
               ex=ex, nex=nex)

A = dist("assigned"); T = dist("tied"); AM = dist("ambiguous")
print("DSFAM42 — completeness of TIED vs ASSIGNED reads (n_decisive=0 vs >=1):")
print(f"{'status':10s} {'n':>5} {'med_exonic_bp':>13} {'mean_bp':>8} {'med_exons':>9} {'single_exon%':>12}")
for lab, D in (("assigned", A), ("ambiguous", AM), ("tied", T)):
    if D:
        print(f"{lab:10s} {D['n']:>5} {D['med_exonic']:>13} {D['mean_exonic']:>8} "
              f"{D['med_exons']:>9} {100*D['single_exon']/D['n']:>11.1f}%")

# classify TIED reads: full-length (floor) vs short (coverage-face), using the ASSIGNED length as the
# 'informative length' reference. A tied read >= 70% of the median assigned exonic length is treated as
# full-length-yet-tied (floor); shorter = coverage-face (a full-length molecule would likely reach a PSV).
ref = A["med_exonic"]
thr = 0.70 * ref
floor = int((T["ex"] >= thr).sum())
cover = int((T["ex"] < thr).sum())
single = T["single_exon"]
print(f"\nReference 'informative length' = median ASSIGNED exonic = {ref} bp; full-length threshold = {int(thr)} bp")
print(f"TIED reads (n={T['n']}):")
print(f"  TRUE FLOOR (>= {int(thr)} bp, full-length yet span no PSV; copies identical over the molecule;"
      f" no read length / EM / ML helps): {floor}  ({100*floor/T['n']:.0f}%)")
print(f"  COVERAGE-FACE (< {int(thr)} bp, short/partial; a full-length capture would likely reach a PSV;"
      f" IMPROVABLE by longer reads): {cover}  ({100*cover/T['n']:.0f}%)")
print(f"  (of the tied, single-exon = {single} ({100*single/T['n']:.0f}%) — no junction axis either)")
print(f"\nINTERPRETATION: tied median {T['med_exonic']}bp vs assigned median {ref}bp — "
      f"{'tied are SHORTER -> coverage-face significant -> longer reads would recover some' if T['med_exonic'] < 0.8*ref else 'tied ~ same length -> dominated by TRUE FLOOR (identical copies), not recoverable'}.")
