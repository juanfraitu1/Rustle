#!/usr/bin/env python3
"""How many of Soto's 83 families does RNA mode actually find -- completely or partially?

Reports at TWO detection standards, because section 14 showed they differ a lot:
  TOUCHED    >=1 predicted copy overlaps the member's interval   (what member-recall has always meant)
  RECOVERED  >=50% of the member's ANNOTATED EXONIC bases are covered by the union of overlapping copies
             (annotation used for SCORING only; members with no annotated exons cannot be judged this way
             and are counted separately rather than silently dropped)

Family verdicts (a family needs >= 2 copies to BE a family, so 1 member is not a family):
  COMPLETE  every member of the family detected
  PARTIAL   >= 2 members detected, but not all
  SINGLETON exactly 1 member detected -- the locus is seen but no family can be formed
  MISSED    0 members detected
"""
import gzip, re, bisect, sys
from collections import defaultdict
import numpy as np

GFF="/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
BED="bench/soto/80_fams.gene_preferred.bed"
CATS=[("definitive (1 genome-wide --cross-chrom pass)","/home/juanfra/winloci_scratch/soto_cache/definitive.copies.tsv"),
      ("per-chrom corrected (21 chrom + xchrom)","/home/juanfra/winloci_scratch/soto_cache/perchrom_catalog.copies.tsv")]

fam=defaultdict(list)
for ln in open(BED):
    c,s,e,name,*_=ln.rstrip("\n").split("\t")
    g,fid=name.split("|")[0],name.split("|")[1]
    fam[fid].append((g,c,int(s),int(e)))

exons=defaultdict(list)
with gzip.open(GFF,"rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="exon": continue
        m=re.search(r'(?:^|;)gene=([^;]+)',f[8])
        if m: exons[(f[0],m.group(1).upper())].append((int(f[3])-1,int(f[4])))
def merge(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def frac(target,cover):
    tot=sum(b-a for a,b in target)
    if not tot: return None
    h=0
    for a,b in target:
        for x,y in cover:
            lo,hi=max(a,x),min(b,y)
            if hi>lo: h+=hi-lo
    return h/tot

print(f"Soto benchmark: {len(fam)} families, {sum(len(v) for v in fam.values())} members\n")
for label,CAT in CATS:
    pred=defaultdict(list)
    ncop=0
    for i,ln in enumerate(open(CAT)):
        if i==0: continue
        f=ln.rstrip("\n").split("\t")
        if len(f)>=7: pred[f[3]].append((int(f[4]),int(f[5]))); ncop+=1
    res={}
    for fid,mem in fam.items():
        touched=recovered=judgeable=0
        for g,c,s,e in mem:
            ov=merge([(a,b) for a,b in pred.get(c,[]) if a<e and s<b])
            if ov: touched+=1
            ex=merge([(a,b) for a,b in exons.get((c,g.upper()),[]) if a<e and s<b])
            if ex:
                judgeable+=1
                if ov and (frac(ex,ov) or 0)>=0.50: recovered+=1
        res[fid]=(len(mem),touched,recovered,judgeable)
    def verdicts(idx, denom_idx=0):
        C=P=S=M=0
        for fid,(n,t,r,j) in res.items():
            d = n if denom_idx==0 else res[fid][denom_idx]
            k = (t,r)[idx-1]
            if d==0: continue
            if k==0: M+=1
            elif k==1: S+=1
            elif k>=d: C+=1
            else: P+=1
        return C,P,S,M
    print(f"=== {label}  ({ncop} copies) ===")
    C,P,S,M=verdicts(1)
    print(f"  TOUCHED (>=1 overlapping copy), denominator = all members")
    print(f"    COMPLETE (all members)  {C:>3}/83 = {100*C/83:.0f}%")
    print(f"    PARTIAL  (>=2, not all) {P:>3}/83 = {100*P/83:.0f}%")
    print(f"    -> family formable (COMPLETE+PARTIAL) {C+P:>3}/83 = {100*(C+P)/83:.0f}%")
    print(f"    SINGLETON (1 member only) {S:>3}   MISSED (0) {M:>3}")
    C,P,S,M=verdicts(2,3)
    tot=C+P+S+M
    print(f"  RECOVERED (>=50% of annotated exons), denominator = members with annotated exons")
    print(f"    COMPLETE {C:>3}/{tot}   PARTIAL {P:>3}/{tot}   -> formable {C+P:>3}/{tot} = {100*(C+P)/max(tot,1):.0f}%")
    print(f"    SINGLETON {S:>3}   MISSED {M:>3}")
    print()

# ---------------------------------------------------------------------------------------------------------
# ISOFORM REQUIREMENT. A family should not count as found when it is represented only by unspliced clusters:
# section 14 showed those can be purely INTRONIC (SRGAP2's two copies lie inside one ~93kb intron, covering
# zero exons) yet still satisfy an overlap test. Requiring at least one grouped member to carry a real
# ISOFORM -- a spliced, multi-exon copy, i.e. a transcript model with observed junctions -- closes that
# loophole using only the catalog (no annotation).
print("=" * 100)
print("ISOFORM-REQUIRED family counts (>=2 members in ONE predicted family AND >=1 of them has a spliced copy)")
for label, CAT in CATS:
    pred2 = defaultdict(list)
    for i, ln in enumerate(open(CAT)):
        if i == 0: continue
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 7: pred2[f[3]].append((int(f[4]), int(f[5]), f[0], int(f[6])))
    found = iso = all_iso = comp = part = 0
    dropped = []
    for fid, mem in fam.items():
        byfam, isom = defaultdict(set), defaultdict(set)
        for g, c, s, e in mem:
            for a, b, pf, nx in pred2.get(c, []):
                if a < e and s < b:
                    byfam[pf].add(g)
                    if nx >= 2: isom[pf].add(g)
        best = None
        for pf, ms in byfam.items():
            if len(ms) >= 2 and (best is None or len(ms) > len(byfam[best])): best = pf
        if best is None: continue
        found += 1
        if len(isom.get(best, ())) >= 1:
            iso += 1
            comp += len(byfam[best]) >= len(mem); part += len(byfam[best]) < len(mem)
            all_iso += len(isom[best]) >= len(byfam[best])
        else:
            dropped.append(fid)
    print(f"\n  {label}")
    print(f"    grouped only (>=2 members in one family)        {found:>3}/83")
    print(f"    + >=1 grouped member carries an isoform         {iso:>3}/83 = {100*iso/83:.0f}%   "
          f"(COMPLETE {comp}, PARTIAL {part})")
    print(f"    + ALL grouped members carry an isoform          {all_iso:>3}/83 = {100*all_iso/83:.0f}%")
    print(f"    dropped as unspliced-only: {len(dropped)}  {sorted(dropped)[:10]}")
