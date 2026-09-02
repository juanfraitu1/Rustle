#!/usr/bin/env python3
"""Evaluate the all-annotation multi-seed projection against the read-only catalog."""
import sys, collections
BASE="/mnt/linuxdisk/home/juanfraitu"
MINCOV=float(sys.argv[1]) if len(sys.argv)>1 else 0.30

def parse_region(s):
    c,r=s.split(":"); a,b=r.split("-"); return c,int(a)-1,int(b)

truth=[parse_region(l.strip()) for l in open(f"{BASE}/o1_oracle/npip31.regions") if l.strip()]
cop=[]
with open(f"{BASE}/npip_cat/arm_f2/cat.copies.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t"); cop.append((f[ci["chrom"]],int(f[ci["start"]]),int(f[ci["end"]])))

pairs=[]      # (seed_interval, landing_interval)
seeds=set()
for l in open(f"{BASE}/seedmode/multiseed.paf"):
    p,q,t,ts,te,strand,nres,blk=l.rstrip("\n").split("\t")
    qc,qa,qb=parse_region(q); ts,te,blk,nres=int(ts),int(te),int(blk),int(nres)
    qlen=qb-qa
    seeds.add((qc,qa,qb))
    if blk < MINCOV*qlen: continue
    ov = (min(qb,te)-max(qa,ts)) if t==qc else 0
    if ov > 0.5*min(qlen, te-ts): continue        # self-hit
    pairs.append(((qc,qa,qb),(t,ts,te),blk/qlen,nres/max(blk,1)))

firing={s for s,_,_,_ in pairs}
print(f"multi-seed projection  (min coverage of seed = {MINCOV:.2f})")
print(f"  seeds supplied          : {len(seeds)} annotated genes/pseudogenes")
print(f"  seeds that FIRE (>=1 non-self copy): {len(firing)} = {len(firing)/max(len(seeds),1):.1%}")
print(f"  non-self copy relations : {len(pairs)}")

# ---- merge every interval that participates (seed side or landing side) into loci
ivs=[]
for s,t,_,_ in pairs: ivs.append(s); ivs.append(t)
ivs=sorted(set(ivs))
loci=[]
for c,a,b in ivs:
    if loci and loci[-1][0]==c and a<=loci[-1][2]:
        loci[-1][2]=max(loci[-1][2],b)
    else: loci.append([c,a,b])
def locus_of(iv):
    c,a,b=iv
    lo,hi=0,len(loci)-1
    while lo<=hi:
        m=(lo+hi)//2
        L=loci[m]
        if L[0]<c or (L[0]==c and L[2]<=a): lo=m+1
        elif L[0]>c or (L[0]==c and L[1]>=b): hi=m-1
        else: return m
    return None
print(f"  distinct loci involved  : {len(loci)}")

# ---- union-find into families
par=list(range(len(loci)))
def find(x):
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
def uni(x,y):
    x,y=find(x),find(y)
    if x!=y: par[x]=y
for s,t,_,_ in pairs:
    i,j=locus_of(s),locus_of(t)
    if i is not None and j is not None: uni(i,j)
fam=collections.defaultdict(list)
for i in range(len(loci)): fam[find(i)].append(i)
multi={k:v for k,v in fam.items() if len(v)>=2}
sizes=sorted((len(v) for v in multi.values()), reverse=True)
print(f"  families (>=2 loci)     : {len(multi)}")
print(f"  copies in them          : {sum(sizes)}")
print(f"  size distribution       : max {sizes[0] if sizes else 0}, "
      f"median {sizes[len(sizes)//2] if sizes else 0}, top10 {sizes[:10]}")

# ---- NPIP recall
covered=[]
memb=[loci[i] for v in multi.values() for i in v]
for t in truth:
    hit=any(m[0]==t[0] and min(t[2],m[2])-max(t[1],m[1])>0 for m in memb)
    covered.append(hit)
print(f"\n  NPIP recall (member of a >=2-locus family): {sum(covered)}/31")
readhit=[any(c==t[0] and min(t[2],b)-max(t[1],a)>0 for c,a,b in cop) for t in truth]
print(f"  read-only catalog                        : {sum(readhit)}/31")
gain=[t for t,h,r in zip(truth,covered,readhit) if h and not r]
lost=[t for t,h,r in zip(truth,covered,readhit) if r and not h]
print(f"  recovered by seeds, missed by reads      : {len(gain)}")
for t in gain: print(f"     {t[0]}:{t[1]}-{t[2]}  len={t[2]-t[1]}")
print(f"  found by reads, missed by seeds          : {len(lost)}")
for t in lost: print(f"     {t[0]}:{t[1]}-{t[2]}  len={t[2]-t[1]}")
