#!/usr/bin/env python3
"""Does the method find the FULL locus of each ground-truth member -- at least for complete genes?

Distinct from the size question (section 8/12), which asked whether ONE predicted copy has the right span.
A locus covered by three fragments is still FOUND, just fragmented. So here: take the UNION of every
predicted copy overlapping a truth member, and ask what fraction of that gene's ANNOTATED EXONIC bases it
covers. Stratified by gene_biotype, because "complete gene" is the question -- protein_coding members are
real gene models, whereas a pseudogene's annotated span may not correspond to anything transcribed.

Reports the exon-coverage distribution, not a single number: "found the locus" is a matter of degree and the
distribution is what makes it honest.
"""
import gzip, re, sys
from collections import defaultdict
import numpy as np

GFF = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
BED = "bench/soto/80_fams.gene_preferred.bed"
CAT = sys.argv[1] if len(sys.argv) > 1 else "/home/juanfra/winloci_scratch/soto_cache/definitive.copies.tsv"

members = []
for ln in open(BED):
    c, s, e, name, *_ = ln.rstrip("\n").split("\t")
    members.append((name.split("|")[0], c, int(s), int(e)))

biotype, exons = {}, defaultdict(list)
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9: continue
        if f[2] in ("gene", "pseudogene"):
            m = re.search(r'(?:^|;)Name=([^;]+)', f[8]); b = re.search(r'gene_biotype=([^;]+)', f[8])
            if m: biotype[(f[0], m.group(1).upper())] = b.group(1) if b else f[2]
        elif f[2] == "exon":
            m = re.search(r'(?:^|;)gene=([^;]+)', f[8])
            if m: exons[(f[0], m.group(1).upper())].append((int(f[3]) - 1, int(f[4])))

pred = defaultdict(list)
for i, ln in enumerate(open(CAT)):
    if i == 0: continue
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 7: pred[f[3]].append((int(f[4]), int(f[5]), int(f[6])))

def merge(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], b))
        else: out.append((a, b))
    return out

def cov_frac(target, cover):
    """fraction of `target` intervals' bases covered by `cover` intervals"""
    tot = sum(b - a for a, b in target)
    if tot == 0: return None
    hit = 0
    for a, b in target:
        for x, y in cover:
            lo, hi = max(a, x), min(b, y)
            if hi > lo: hit += hi - lo
    return hit / tot

rows = []
for gene, c, s, e in members:
    key = (c, gene.upper())
    ex = merge([(a, b) for a, b in exons.get(key, []) if a < e and s < b])
    if not ex: continue
    bt = biotype.get(key, "unannotated")
    ov = [(a, b, n) for a, b, n in pred.get(c, []) if a < e and s < b]
    spans   = merge([(a, b) for a, b, _n in ov])
    spl     = merge([(a, b) for a, b, n in ov if n >= 2])
    rows.append((gene, bt, len(ov), cov_frac(ex, spans) or 0.0, cov_frac(ex, spl) or 0.0))

def show(label, sel):
    r = [x for x in rows if sel(x)]
    if not r: return
    det = [x for x in r if x[2] > 0]
    a = np.array([x[3] for x in det]) if det else np.array([0.0])
    print(f"\n{label}  (n={len(r)})")
    print(f"  detected at all (>=1 overlapping copy) : {len(det)}/{len(r)} = {100*len(det)/len(r):.0f}%")
    print(f"  exon coverage over DETECTED members    : median {np.median(a):.2f}   mean {a.mean():.2f}")
    for t in (0.90, 0.75, 0.50, 0.25):
        n = int((a >= t).sum())
        print(f"    >= {int(t*100):>3}% of annotated exons covered : {n:>3}/{len(det)}  ({100*n/max(len(det),1):.0f}% of detected, "
              f"{100*n/len(r):.0f}% of all)")

show("ALL members with annotated exons", lambda x: True)
show("COMPLETE GENES (protein_coding)", lambda x: x[1] == "protein_coding")
show("transcribed_pseudogene", lambda x: x[1] == "transcribed_pseudogene")
show("pseudogene", lambda x: x[1] == "pseudogene")

pc = [x for x in rows if x[1] == "protein_coding" and x[2] > 0]
print("\nprotein_coding members with the WORST exon coverage:")
for gene, bt, n, cv, cs in sorted(pc, key=lambda x: x[3])[:12]:
    print(f"  {gene:14s} copies={n:<3} exon_cov={cv:.2f}  (spliced-only {cs:.2f})")
print("\nprotein_coding: exon coverage using SPLICED copies only:")
sp = np.array([x[4] for x in pc])
print(f"  median {np.median(sp):.2f}   >=90%: {int((sp>=0.9).sum())}/{len(pc)}   >=50%: {int((sp>=0.5).sum())}/{len(pc)}")
