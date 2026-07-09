#!/usr/bin/env python3
"""Cross-modal validation of the RNA E_r family catalog (companion spec Axes B + C).

Axis B  RNA-98 vs SEDEF SD98 : what fraction of the RNA >=98% families are DNA-confirmed (>=1 copy inside a
                               SEDEF SD98 >=98% segdup interval) -> a circularity-free RNA-vs-DNA check.
Axis C  Soto famCN table     : annotate each RNA family's copies with RefSeq genes, join to the phased-assembly
                               oracle (diploid_cn_oracle asm_hapCN), and compare famCN (genome projection) to
                               the assembly copy number.

Inputs (paths hardcoded to the run dir; override via argv):
  <cat>.copies.tsv  <cat>.famcn.tsv   (from gw_family_catalog --homology-primary [--min-identity 0.98] --enumerate-copies)
  sd98_intervals.bed   genes.bed      (prepped)   diploid_cn_oracle.tsv
Uses bedtools for interval overlaps.
"""
import subprocess, sys, os, collections

SC = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/fae6260a-25c9-415b-a8a3-d36c0bc4ca96/scratchpad/xval"
CAT = sys.argv[1] if len(sys.argv) > 1 else f"{SC}/rna98"
ORACLE = "/mnt/c/Users/jfris/Desktop/Rustle/bench/diploid_cn_oracle.tsv"
BT = "/home/juanfra/miniforge3/bin/bedtools"

def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout

# --- load catalog copies -> per-family loci + a BED ---
copies = collections.defaultdict(list)   # fam -> [(chrom,start,end)]
with open(f"{CAT}.copies.tsv") as f:
    next(f)
    for l in f:
        c = l.rstrip("\n").split("\t")
        if len(c) < 6: continue
        copies[c[0]].append((c[3], int(c[4]), int(c[5])))
n_fam = len(copies)
copies_bed = f"{SC}/cat_copies.bed"
with open(copies_bed, "w") as o:
    for fam, loci in copies.items():
        for ch, s, e in loci:
            o.write(f"{ch}\t{s}\t{e}\t{fam}\n")

# --- famCN ---
famcn = {}
if os.path.exists(f"{CAT}.famcn.tsv"):
    with open(f"{CAT}.famcn.tsv") as f:
        next(f)
        for l in f:
            c = l.rstrip("\n").split("\t")
            if len(c) >= 3: famcn[c[0]] = (int(c[1]), int(c[2]))  # (n_rna, famCN)

# === Axis B: DNA-confirmation via SEDEF SD98 ===
conf = sh(f"{BT} intersect -a {copies_bed} -b {SC}/sd98_intervals.bed -u")
fam_dna_confirmed = set(l.split("\t")[3] for l in conf.strip().split("\n") if l)
b_conf = len(fam_dna_confirmed)
print("="*70)
print("AXIS B — RNA-98 families DNA-confirmed by SEDEF SD98 (>=98% segdup)")
print(f"  RNA families: {n_fam}")
print(f"  DNA-confirmed (>=1 copy in a SD98 interval): {b_conf}  ({100*b_conf/max(n_fam,1):.1f}%)")
print(f"  RNA-only (no SD98 overlap = below the DNA >=98% regime, i.e. more divergent): {n_fam-b_conf}")

# === Axis C: annotate families with genes, join oracle asm_hapCN, compare famCN ===
ann = sh(f"{BT} intersect -a {copies_bed} -b {SC}/genes.bed -wa -wb")
fam_genes = collections.defaultdict(set)
for l in ann.strip().split("\n"):
    if not l: continue
    c = l.split("\t")
    fam_genes[c[3]].add(c[7])  # gene Name is col8 (0-based 7) of the -wb append (genes.bed col4)

# oracle: gene -> (asm_hapCN, class)
oracle = {}
with open(ORACLE) as f:
    hdr = next(f).rstrip("\n").split("\t")
    gi, ai, ci = hdr.index("gene"), hdr.index("asm_hapCN"), hdr.index("class")
    for l in f:
        c = l.rstrip("\n").split("\t")
        if len(c) > max(gi, ai, ci) and c[gi] != "NA":
            oracle[c[gi]] = (c[ai], c[ci])

print("="*70)
print("AXIS C — Soto famCN table: RNA family -> genes -> assembly asm_hapCN vs famCN(projection)")
print(f"{'family':10} {'famCN':>6} {'n_rna':>6} {'asm_hapCN':>9}  genes(oracle-matched)")
rows = 0
for fam, genes in sorted(fam_genes.items()):
    matched = [g for g in genes if g in oracle]
    if not matched: continue
    hapcns = sorted(set(oracle[g][0] for g in matched))
    nr, fc = famcn.get(fam, (len(copies[fam]), len(copies[fam])))
    print(f"{fam:10} {fc:>6} {nr:>6} {','.join(hapcns):>9}  {','.join(sorted(matched)[:6])}")
    rows += 1
    if rows >= 40: print("  ... (truncated)"); break
print(f"  ({rows} RNA families matched to an oracle-annotated multi-copy gene)")
