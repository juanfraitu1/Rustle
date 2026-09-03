#!/usr/bin/env python3
"""Tier every seeded candidate locus by the READ evidence at it (ledger §6be step 2).

The seeded projection says WHERE to look; this says what the reads support there.
Emits three tiers so the seeded mode can never be mistaken for a definition:
  O1 member    - the family has >=2 loci with enough PRIMARY support to build a node
  recoverable  - >=2 loci expressed once AS-tied secondaries are counted (the starved stratum)
  O3 candidate - fewer than 2 expressed loci: duplicated in DNA, not seen in this tissue
"""
import collections, subprocess, sys
BASE="/mnt/linuxdisk/home/juanfraitu"
BAM=f"{BASE}/npip_cat/npip3.bam"
MINCOV=float(sys.argv[1]) if len(sys.argv)>1 else 0.30
PRIM_FLOOR, SEC_FLOOR = 3, 50          # node floor (NODE_MIN_READS=2 +1 margin); starved bar

def pr(s):
    c,r=s.split(":"); a,b=r.split("-"); return c,int(a)-1,int(b)

pairs=[]
for l in open(f"{BASE}/seedmode/multiseed.paf"):
    p,q,t,ts,te,st,nres,blk=l.rstrip("\n").split("\t")
    qc,qa,qb=pr(q); ts,te,blk=int(ts),int(te),int(blk); qlen=qb-qa
    if blk<MINCOV*qlen: continue
    if t==qc and (min(qb,te)-max(qa,ts))>0.5*min(qlen,te-ts): continue
    pairs.append(((qc,qa,qb),(t,ts,te)))
ivs=sorted({x for s,t in pairs for x in (s,t)})
loci=[]
for c,a,b in ivs:
    if loci and loci[-1][0]==c and a<=loci[-1][2]: loci[-1][2]=max(loci[-1][2],b)
    else: loci.append([c,a,b])
byc=collections.defaultdict(list)
for i,L in enumerate(loci): byc[L[0]].append(i)
def lof(iv):
    c,a,b=iv
    for i in byc[c]:
        L=loci[i]
        if min(b,L[2])-max(a,L[1])>0: return i
par=list(range(len(loci)))
def f(x):
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
for s,t in pairs:
    i,j=lof(s),lof(t)
    if i is not None and j is not None:
        a,b=f(i),f(j)
        if a!=b: par[a]=b
fam=collections.defaultdict(list)
for i in range(len(loci)): fam[f(i)].append(i)
multi=[v for v in fam.values() if len(v)>=2]

# ⛔ FIXED 2026-09-02 (ledger §6cn, meta-ledger row 12). This used `samtools view -c REGION`, which
# counts every read whose SPAN overlaps the interval -- including reads that splice straight OVER it
# contributing no aligned base. On this substrate that inflated the "expressed" tier from 21 to 49 and
# the headline from 40% to 50%, because seeded candidate intervals routinely sit inside a larger gene's
# intron. A read spliced OVER a locus is no evidence for it; count ALIGNED BLOCKS.
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locus_reads import aligned_blocks, reads_with_block_in, SAMTOOLS

def _block_count(L, flag_args):
    """Reads matching `flag_args` with >= 1 aligned base inside L."""
    out = subprocess.run([SAMTOOLS, "view", *flag_args, BAM, f"{L[0]}:{L[1]+1}-{L[2]}"],
                         capture_output=True, text=True).stdout
    n = 0
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) > 5 and f[5] != "*" and any(bs < L[2] and L[1] < be
                                              for bs, be in aligned_blocks(int(f[3]) - 1, f[5])):
            n += 1
    return n

def depth(L):
    p = reads_with_block_in(BAM, L[0], L[1], L[2])      # primaries, -F 2308
    s = _block_count(L, ["-f", "256"])                  # secondaries, same block rule
    return p, s
dep={}
for v in multi:
    for i in v:
        if i not in dep: dep[i]=depth(loci[i])

cop=[]
with open(f"{BASE}/npip_cat/arm_f2/cat.copies.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        fl=l.rstrip("\n").split("\t")
        cop.append((fl[ci["chrom"]],int(fl[ci["start"]]),int(fl[ci["end"]]),fl[ci["family_id"]]))
def incat(L):
    return {fid for c,a,b,fid in cop if c==L[0] and min(L[2],b)-max(L[1],a)>0}

tiers=collections.Counter(); tier_loci=collections.Counter()
newfam=[]; newcopies=0
for v in multi:
    nprim=sum(1 for i in v if dep[i][0]>=PRIM_FLOOR)
    nexpr=sum(1 for i in v if dep[i][0]>=PRIM_FLOOR or dep[i][1]>=SEC_FLOOR)
    known=set().union(*[incat(loci[i]) for i in v]) if v else set()
    t = "O1 member" if nprim>=2 else ("recoverable" if nexpr>=2 else "O3 candidate")
    tiers[t]+=1; tier_loci[t]+=len(v)
    if t in ("O1 member","recoverable") and not known:
        newfam.append((t,len(v),nexpr)); newcopies+=nexpr
print(f"seeded candidate families (>=2 loci): {len(multi)}   loci: {len(dep)}")
for t in ("O1 member","recoverable","O3 candidate"):
    print(f"  {t:<14} families {tiers[t]:>3}   loci {tier_loci[t]:>4}")
print(f"\nfamilies in an expressed tier with NO overlap of the read catalog: {len(newfam)}"
      f"  ({newcopies} expressed loci)")
for t,n,e in sorted(newfam, key=lambda x:-x[2])[:10]:
    print(f"   {t:<12} {n:>3} loci, {e:>3} expressed")

# ── §6be's "loci in 0-CORROBORATED families" breakdown, made REPRODUCIBLE (added 2026-09-02, §6cn).
# The ledger published 49 / 103 / 153 over 305 loci from an ad-hoc slice this script never printed, so
# the figure could not be re-derived when `depth()` was found to be span-based (meta-ledger row 12).
# It is emitted here so it is never again a number without a command that produces it.
zero = [v for v in multi if not set().union(*[incat(loci[i]) for i in v])]
zl = [i for v in zero for i in v]
expressed = sum(1 for i in zl if dep[i][0] >= PRIM_FLOOR)
starved = sum(1 for i in zl if dep[i][0] < PRIM_FLOOR and dep[i][1] >= SEC_FLOOR)
silent = len(zl) - expressed - starved
print(f"\n0-CORROBORATED families: {len(zero)}   loci: {len(zl)}")
print(f"  >=3 primaries WITH AN ALIGNED BLOCK — expressed, node never built : {expressed}")
print(f"  <3 primaries but >=50 secondaries — STARVED                       : {starved}")
print(f"  effectively silent                                                : {silent}")
if zl:
    print(f"  => {expressed + starved}/{len(zl)} = {(expressed + starved) / len(zl):.1%} EXPRESSED but absent from the catalog")
