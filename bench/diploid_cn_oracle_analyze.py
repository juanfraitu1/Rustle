#!/usr/bin/env python3
"""Corroboration analysis for the diploid CN oracle (non-circular DNA vs RNA)."""
import os
from collections import defaultdict
import numpy as np
from scipy.stats import spearmanr, pearsonr

SCRATCH = "/home/juanfra/winloci_scratch"
VALIDATED = os.path.join(SCRATCH, "validated_families.tsv")
ORACLE = "bench/diploid_cn_oracle.tsv"
CHRX, CHRY = "NC_073247.2", "NC_073248.2"

# family -> set of chroms
fam_chroms = defaultdict(set)
with open(VALIDATED) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        fam_chroms[f[0]].add(f[2])

rows = []
with open(ORACLE) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        d = dict(zip(hdr, line.rstrip("\n").split("\t")))
        rows.append(d)

def num(x):
    try:
        return float(x)
    except Exception:
        return None

def sexlink(fam):
    ch = fam_chroms[fam]
    if ch == {CHRX}:
        return "X"
    if ch == {CHRY}:
        return "Y"
    if CHRX in ch or CHRY in ch:
        return "mixed_sex"
    return "auto"

# ---- assembly haploid CN: sex-aware ----
for d in rows:
    cm, cp = int(d["hap_CN_mat"]), int(d["hap_CN_pat"])
    sl = sexlink(d["family"])
    d["sex"] = sl
    if sl == "X":
        d["asm_hapCN"] = cm            # X present only on maternal
    elif sl == "Y":
        d["asm_hapCN"] = cp            # Y present only on paternal
    else:
        d["asm_hapCN"] = max(cm, cp)   # autosome: robust to a 1-hap assembly gap
    # sex-aware class
    if sl in ("X", "Y"):
        h = d["asm_hapCN"]
        d["class_sexaware"] = ("multi_copy" if h >= 2 else
                               "single_locus_allele" if h == 1 else "absent")
    else:
        d["class_sexaware"] = d["class"]

# ---- distributions ----
def dist(key):
    c = defaultdict(int)
    for d in rows:
        c[d[key]] += 1
    return dict(sorted(c.items(), key=lambda kv: (-kv[1], str(kv[0]))))

print("=== raw class (literal rule, TSV) ===", dist("class"))
print("=== sex-aware class ===", dist("class_sexaware"))
print("=== sex linkage ===", dist("sex"))

# how many cnv are just sex-linkage
cnv = [d for d in rows if d["class"] == "cnv"]
cnv_sex = [d for d in cnv if d["sex"] in ("X", "Y", "mixed_sex")]
print("cnv total=%d of which sex-linked=%d (autosomal-CNV=%d)" %
      (len(cnv), len(cnv_sex), len(cnv) - len(cnv_sex)))
print("  autosomal CNV families:",
      [(d["family"], d["gene"], d["hap_CN_mat"], d["hap_CN_pat"]) for d in cnv
       if d["sex"] == "auto"])

# ---- corroboration: assembly CN vs reference n_loci (all 196) ----
def corr(xs, ys, label):
    xs, ys = np.array(xs, float), np.array(ys, float)
    rho, p = spearmanr(xs, ys)
    r, pp = pearsonr(xs, ys)
    # permutation null
    rng = np.random.default_rng(0)
    perm = np.array([spearmanr(xs, rng.permutation(ys))[0] for _ in range(10000)])
    pnull = (np.sum(np.abs(perm) >= abs(rho)) + 1) / (len(perm) + 1)
    print("%-42s n=%3d  spearman rho=%.3f (p=%.2e, perm_p=%.4f)  pearson r=%.3f" %
          (label, len(xs), rho, p, pnull, r))
    return rho

print("\n=== CORROBORATION (Spearman + 10k permutation null) ===")
allf = rows
corr([d["asm_hapCN"] for d in allf], [float(d["n_loci_ref"]) for d in allf],
     "asm_hapCN vs n_loci_ref (all 196)")

# a1 subset with K_read / chi_H
a1 = [d for d in rows if num(d["K_read"]) is not None]
corr([d["asm_hapCN"] for d in a1], [num(d["K_read"]) for d in a1],
     "asm_hapCN vs K_read (a1, NON-circular)")
corr([d["asm_hapCN"] for d in a1], [num(d["chi_H"]) for d in a1],
     "asm_hapCN vs chi_H (a1)")
corr([d["asm_hapCN"] for d in a1], [float(d["n_loci_ref"]) for d in a1],
     "asm_hapCN vs n_loci_ref (a1 subset)")
sed = [d for d in a1 if num(d["sedef_partners"]) is not None]
corr([d["asm_hapCN"] for d in sed], [num(d["sedef_partners"]) for d in sed],
     "asm_hapCN vs sedef_partners (a1)")

# how many a1 families have asm_hapCN>=2 (assembly-confirmed multi-copy)
a1_conf = sum(1 for d in a1 if d["asm_hapCN"] >= 2)
print("\na1 families assembly-confirmed multi-copy (asm_hapCN>=2): %d/%d" %
      (a1_conf, len(a1)))

# ---- disagreements ----
print("\n=== STRONG DISAGREEMENTS ===")
print("-- asm_hapCN vs n_loci_ref (ratio >=3x or <=1/3, min diff 3) --")
dis = []
for d in rows:
    a, n = d["asm_hapCN"], float(d["n_loci_ref"])
    if a == 0:
        continue
    if abs(a - n) >= 3 and (a >= 3 * n or n >= 3 * a):
        dis.append(d)
for d in sorted(dis, key=lambda x: -abs(x["asm_hapCN"] - float(x["n_loci_ref"]))):
    print("  fam%-4s %-16s n_loci_ref=%-3s asm_hapCN=%-3d (mat=%s pat=%s) sex=%s K_read=%s" %
          (d["family"], d["gene"], d["n_loci_ref"], d["asm_hapCN"],
           d["hap_CN_mat"], d["hap_CN_pat"], d["sex"], d["K_read"]))

print("-- asm_hapCN vs K_read (a1; |diff|>=3 and ratio>=3x) --")
for d in a1:
    a, k = d["asm_hapCN"], num(d["K_read"])
    if k == 0:
        note = "K_read=0 (unresolved by reads)"
        if a >= 3:
            print("  fam%-4s %-16s K_read=0 asm_hapCN=%-3d n_loci_ref=%s sex=%s [%s]" %
                  (d["family"], d["gene"], a, d["n_loci_ref"], d["sex"], note))
        continue
    if abs(a - k) >= 3 and (a >= 3 * k or k >= 3 * a):
        print("  fam%-4s %-16s K_read=%-3g asm_hapCN=%-3d n_loci_ref=%s sex=%s" %
              (d["family"], d["gene"], k, a, d["n_loci_ref"], d["sex"]))

# ---- A1 het-risk survivors resolution ----
print("\n=== A1 het-risk survivors resolution ===")
for fam in ("94", "145"):
    d = next(x for x in rows if x["family"] == fam)
    print("  fam%s %-14s mat=%s pat=%s sex=%s asm_hapCN=%d -> %s" %
          (fam, d["gene"], d["hap_CN_mat"], d["hap_CN_pat"], d["sex"],
           d["asm_hapCN"], d["class_sexaware"]))

# all het_risk families (from a1) resolution
HETRISK = ["114","58","65","71","87","91","120","94","145","151","69"]
print("\n=== all a1 het_risk=1 families -> oracle verdict ===")
for fam in HETRISK:
    d = next((x for x in rows if x["family"] == fam), None)
    if d:
        print("  fam%-4s %-16s mat=%s pat=%s sex=%-5s asm_hapCN=%d -> %s" %
              (fam, d["gene"], d["hap_CN_mat"], d["hap_CN_pat"], d["sex"],
               d["asm_hapCN"], d["class_sexaware"]))
