#!/usr/bin/env python3
"""Which GROUPING is better, on IDENTICAL loci? (the fair form of "better as a catalog")

Both earlier framings were biased: seeded loci ARE annotated gene intervals, so "points at a
gene" is near-entailed for them and incidental for read loci. Here the locus set is held fixed
-- only loci present in BOTH catalogs are scored -- and each locus gets the SAME CDS. The only
variable left is which loci each catalog puts in one family.

Independent relation: a CDS-CDS edge under the project's own E_r rule (id >= 0.60, cov >= 0.50
of the shorter). Scored over PAIRS of shared loci, so neither catalog's family count enters.
"""
import collections, os, subprocess, sys, itertools
BASE="/mnt/linuxdisk/home/juanfraitu"
WORK=f"{BASE}/cdscoh"; os.makedirs(WORK, exist_ok=True)
REF=f"{BASE}/npip_cat/npip3_contigs.fa"; GFF=f"{BASE}/winloci_data/GGO_genomic.gff"
CONTIGS={"NC_073241.2","NC_073242.2","NC_073244.2"}
MINCOV=0.30; ER_ID=0.60; ER_COV=0.50

exec(open("bench/catalog_cds_coherence.py").read().split("# ------------------------------------------------------------------ catalogs")[0]
     .split("import collections, os, random, subprocess, sys")[1].replace('random.seed(11)',''), globals())

def seeded_families():
    def pr(s):
        c,r=s.split(":"); a,b=r.split("-"); return c,int(a)-1,int(b)
    pairs=[]
    for l in open(f"{BASE}/seedmode/multiseed.paf"):
        p,q,t,ts,te,st,nres,blk=l.rstrip("\n").split("\t")
        qc,qa,qb=pr(q); ts,te,blk=int(ts),int(te),int(blk); qlen=qb-qa
        if blk<MINCOV*qlen: continue
        if t==qc and (min(qb,te)-max(qa,ts))>0.5*min(qlen,te-ts): continue
        pairs.append(((qc,qa,qb),(t,ts,te)))
    ivs=sorted({x for s,t in pairs for x in (s,t)}); loci=[]
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
    return loci,{i:f"SEED{f(i)}" for i in range(len(loci))}

loci,seedfam = seeded_families()
readcop=[]
with open(f"{BASE}/npip_cat/arm_f2/cat.copies.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t")
        readcop.append((f[ci["chrom"]],int(f[ci["start"]]),int(f[ci["end"]]),f[ci["family_id"]]))

# ---- shared loci: a seeded locus that overlaps exactly one read copy (and vice versa) ----
shared=[]
for i,L in enumerate(loci):
    hits=[r for r in readcop if r[0]==L[0] and min(L[2],r[2])-max(L[1],r[1])>0]
    if not hits: continue
    # max-overlap, not 1:1: requiring a unique partner excluded 98% of the shared set and left
    # 4 positive pairs, which is not a measurement.
    r=max(hits, key=lambda x: min(L[2],x[2])-max(L[1],x[1]))
    g=best_cds_in(L[0],L[1],L[2])
    if not g: continue
    shared.append({"iv":tuple(L),"seed":seedfam[i],"read":r[3],"cds":g})
sys.stderr.write(f"shared loci with a CDS: {len(shared)}\n")

# ---- one CDS per shared locus, all-vs-all, E_r ----
regions=[]
for n,s in enumerate(shared):
    for ea,eb in s["cds"][6]: regions.append(f"{s['cds'][0]}:{ea+1}-{eb}")
open(f"{WORK}/shared.regions","w").write("\n".join(regions)+"\n")
seq=subprocess.run(["samtools","faidx","-r",f"{WORK}/shared.regions",REF],
                   capture_output=True,text=True).stdout
blocks=[];cur=[]
for line in seq.split("\n"):
    if line.startswith(">"):
        if cur: blocks.append("".join(cur))
        cur=[]
    elif line: cur.append(line)
if cur: blocks.append("".join(cur))
comp=str.maketrans("ACGTacgtN","TGCAtgcaN")
with open(f"{WORK}/shared.fa","w") as o:
    k=0
    for n,s in enumerate(shared):
        m=len(s["cds"][6]); q="".join(blocks[k:k+m]); k+=m
        if s["cds"][3]=="-": q=q.translate(comp)[::-1]
        o.write(f">L{n}\n{q}\n")
paf=subprocess.run(["minimap2","-x","asm20","-c","-N","200","-p","0.02","-t","4",
                    f"{WORK}/shared.fa",f"{WORK}/shared.fa"],capture_output=True,text=True).stdout
cdsedge=set()
for l in paf.split("\n"):
    if not l: continue
    f=l.split("\t"); q,ql,t,tl=f[0],int(f[1]),f[5],int(f[6])
    nres,blk=int(f[9]),int(f[10])
    if q==t: continue
    if nres/max(blk,1)<ER_ID or blk<ER_COV*min(ql,tl): continue
    a,b=int(q[1:]),int(t[1:]); cdsedge.add((min(a,b),max(a,b)))

pos=len(cdsedge); n=len(shared); tot=n*(n-1)//2
print(f"shared loci scored          : {n}      pairs: {tot}")
print(f"CDS-homologous pairs (E_r)  : {pos}   base rate {pos/max(tot,1):.4f}")
print()
print(f"{'grouping':<10}{'same-family pairs':>19}{'recall':>10}{'precision':>11}{'F1':>8}")
for key in ("seed","read"):
    same={(i,j) for i,j in itertools.combinations(range(n),2)
          if shared[i][key]==shared[j][key]}
    tp=len(same & cdsedge)
    rec=tp/pos if pos else 0
    prec=tp/len(same) if same else 0
    f1=2*rec*prec/(rec+prec) if rec+prec else 0
    print(f"{key:<10}{len(same):>19}{rec:>10.4f}{prec:>11.4f}{f1:>8.4f}")
