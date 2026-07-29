#!/usr/bin/env python3
"""ANNOTATION-FREE test for "this unspliced cluster is intronic pre-mRNA, not a transcript".

Section 14 found unspliced clusters sitting inside introns (SRGAP2's two copies lie in one ~93kb intron) that
satisfy overlap-based member recall while covering zero exons. Excluding them via a GFF would break the
method's annotation-independence, which is the point of the approach.

The reads already carry the answer. A spliced read's `N` CIGAR gap IS an intron, observed directly. So:

    a cluster is INTRONIC  iff  spliced reads at that locus have an intron that CONTAINS it

That is a witness, not a threshold: the same reads that build the catalog testify that the interval is spliced
OUT of the mature transcript. Annotation is used here ONLY to grade the rule, never to compute it.

Label (annotation, evaluation only): intronic = the cluster overlaps 0 annotated exonic bases.
Prediction (reads only)           : n_reads whose intron contains >= CONTAIN_FRAC of the cluster.
"""
import subprocess, gzip, re, sys
from collections import defaultdict
import numpy as np

SAM = "/home/juanfra/miniforge3/bin/samtools"
BAM = "/home/juanfra/winloci_scratch/soto_cache/soto_regions.bam"
GFF = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
CAT = "/home/juanfra/winloci_scratch/soto_cache/definitive.copies.tsv"
BED = "bench/soto/80_fams.gene_preferred.bed"
CONTAIN_FRAC = 0.90

members = []
for ln in open(BED):
    c, s, e, name, *_ = ln.rstrip("\n").split("\t")
    members.append((name.split("|")[0], c, int(s), int(e)))

exons = defaultdict(list)
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        exons[f[0]].append((int(f[3]) - 1, int(f[4])))
for c in exons:
    iv = sorted(exons[c]); out = []
    for a, b in iv:
        if out and a <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], b))
        else: out.append((a, b))
    exons[c] = out

def exon_overlap_bp(c, s, e):
    iv = exons.get(c, []); lo = 0
    import bisect
    i = bisect.bisect_left([x[1] for x in iv], s)
    tot = 0
    while i < len(iv) and iv[i][0] < e:
        tot += max(0, min(e, iv[i][1]) - max(s, iv[i][0])); i += 1
    return tot

# unspliced predicted copies inside a Soto member
clusters = []
for i, ln in enumerate(open(CAT)):
    if i == 0: continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or int(f[6]) != 1: continue
    c, s, e = f[3], int(f[4]), int(f[5])
    if any(tc == c and ts < e and s < te for _g, tc, ts, te in members):
        clusters.append((c, s, e, int(f[8])))

def intron_witnesses(c, s, e):
    """n spliced primary reads whose N-gap contains >= CONTAIN_FRAC of [s,e); reads only."""
    out = subprocess.run([SAM, "view", "-F", "2308", BAM, f"{c}:{s+1}-{e}"],
                         capture_output=True, text=True).stdout
    n = 0
    L = e - s
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 6 or "N" not in f[5]: continue
        ref = int(f[3]) - 1
        for num, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
            num = int(num)
            if op == "N":
                a, b = ref, ref + num
                if min(e, b) - max(s, a) >= CONTAIN_FRAC * L: n += 1; break
                ref = b
            elif op in "M=XD": ref += num
    return n

rows = []
for c, s, e, nr in clusters:
    lab_intronic = exon_overlap_bp(c, s, e) == 0
    rows.append((lab_intronic, intron_witnesses(c, s, e), c, s, e, nr))

P = sum(1 for r in rows if r[0]); N = len(rows) - P
print(f"unspliced clusters inside a Soto member: {len(rows)}")
print(f"  annotation label: INTRONIC (0 exonic bp) {P}   exon-overlapping {N}   base rate {P/max(len(rows),1):.3f}")
print(f"\nrule: >= k spliced reads whose intron CONTAINS the cluster (annotation-free)")
print(f"  {'k':>3} {'TP':>5} {'FP':>5} {'FN':>5} {'precision':>10} {'recall':>7} {'F1':>6}")
for k in (1, 2, 3, 5, 10):
    tp = sum(1 for l, w, *_ in rows if l and w >= k)
    fp = sum(1 for l, w, *_ in rows if not l and w >= k)
    fn = sum(1 for l, w, *_ in rows if l and w < k)
    pr = tp / max(tp + fp, 1); rc = tp / max(tp + fn, 1)
    print(f"  {k:>3} {tp:>5} {fp:>5} {fn:>5} {pr:>10.3f} {rc:>7.3f} {2*pr*rc/max(pr+rc,1e-9):>6.3f}")
