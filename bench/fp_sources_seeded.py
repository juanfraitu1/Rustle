#!/usr/bin/env python3
"""Where do the seeded catalog's false merges come from, and can any feature filter them?

The §6bf comparison was confined to the 118 loci both catalogs see. FP-source analysis is not,
so this uses EVERY gamma-refined seeded locus that carries a CDS (~10x the pairs).

Label: a same-family pair is CONCORDANT if the two loci's CDS align under the project's own E_r
rule, DISCORDANT otherwise. DISCORDANT is a CANDIDATE false positive, not a proven one -- real
paralogues can diverge past 0.60 in coding sequence. Features are scored TP-vs-FP so a filter can
be judged on the trade it makes, never on FP removal alone.
"""
import bisect, collections, os, subprocess, sys
BASE="/mnt/linuxdisk/home/juanfraitu"
REF=f"{BASE}/_from_wsl/winloci_scratch/GGO.fasta"; GFF=f"{BASE}/winloci_data/GGO_genomic.gff"
WORK=f"{BASE}/fpsrc"; os.makedirs(WORK,exist_ok=True)
MINCOV=0.30; ER_ID=0.60; ER_COV=0.50

per=collections.defaultdict(lambda:{"chrom":None,"strand":None,"ex":[]})
with open(GFF) as fh:
    for l in fh:
        if l[0]=="#": continue
        f=l.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="CDS": continue
        pid=None
        for kv in f[8].split(";"):
            if kv.startswith("protein_id="): pid=kv[11:]; break
        if pid is None: continue
        d=per[pid]; d["chrom"]=f[0]; d["strand"]=f[6]; d["ex"].append((int(f[3])-1,int(f[4])))
bychrom=collections.defaultdict(list)
for pid,d in per.items():
    d["ex"].sort(); ln=sum(b-a for a,b in d["ex"])
    bychrom[d["chrom"]].append(((d["ex"][0][0]+d["ex"][-1][1])//2, ln, d["strand"], d["ex"]))
for k in bychrom: bychrom[k].sort()
mids={k:[x[0] for x in v] for k,v in bychrom.items()}
def best_cds_in(c,a,b):
    v=bychrom.get(c)
    if not v: return None
    lo=bisect.bisect_left(mids[c],a); hi=bisect.bisect_right(mids[c],b)
    best=None
    for i in range(lo,hi):
        if best is None or v[i][1]>best[1]: best=v[i]
    return (c,)+best if best else None

def pr(s):
    c,r=s.rsplit(":",1); a,b=r.split("-"); return c,int(a)-1,int(b)
rel=[]
for l in open(f"{BASE}/seedgw/multiseed_gw.paf"):
    p,q,t,ts,te,st,nres,blk=l.rstrip("\n").split("\t")
    qc,qa,qb=pr(q); ts,te,blk,nres=int(ts),int(te),int(blk),int(nres); qlen=qb-qa
    if blk<MINCOV*qlen: continue
    if t==qc and (min(qb,te)-max(qa,ts))>0.5*min(qlen,te-ts): continue
    rel.append(((qc,qa,qb),(t,ts,te),blk/qlen,nres/max(blk,1),blk))
ivs=sorted({x for a,b,_,_,_ in rel for x in (a,b)}); loci=[]
for c,a,b in ivs:
    if loci and loci[-1][0]==c and a<=loci[-1][2]: loci[-1][2]=max(loci[-1][2],b)
    else: loci.append([c,a,b])
byc=collections.defaultdict(list)
for i,L in enumerate(loci): byc[L[0]].append(i)
starts={k:[loci[i][1] for i in v] for k,v in byc.items()}
def lof(iv):
    c,a,b=iv; v=byc.get(c)
    if not v: return None
    j=bisect.bisect_right(starts[c],b)-1
    while j>=0:
        i=v[j]; L=loci[i]
        if min(b,L[2])-max(a,L[1])>0: return i
        if L[1]<a-5_000_000: break
        j-=1
    return None
direct={}                                  # (i,j) -> best (coverage, identity, block)
edges=[]
for s,t,cov,idn,blk in rel:
    i,j=lof(s),lof(t)
    if i is None or j is None or i==j: continue
    k=(min(i,j),max(i,j)); edges.append(k)
    if k not in direct or blk>direct[k][2]: direct[k]=(cov,idn,blk)
r=subprocess.run([f"{BASE}/rustle_target/release/gamma_refine"],
                 input="\n".join(f"{u}\t{v}" for u,v in sorted(set(edges)))+"\n",
                 capture_output=True,text=True)
sys.stderr.write(r.stderr)
famof={}
for line in r.stdout.split("\n"):
    if line:
        b,m=line.split("\t"); famof[int(m)]=b
fam=collections.defaultdict(list)
for m,b in famof.items(): fam[b].append(m)
fam={k:v for k,v in fam.items() if len(v)>=2}
sys.stderr.write(f"gamma-refined families: {len(fam)}  loci: {sum(len(v) for v in fam.values())}\n")

# ---- CDS per locus ----
cds={}
for b,mem in fam.items():
    for i in mem:
        L=loci[i]; g=best_cds_in(L[0],L[1],L[2])
        if g: cds[i]=g
sys.stderr.write(f"loci with a CDS: {len(cds)}\n")
order=sorted(cds)
regions=[f"{cds[i][0]}:{ea+1}-{eb}" for i in order for ea,eb in cds[i][4]]
open(f"{WORK}/c.regions","w").write("\n".join(regions)+"\n")
seq=subprocess.run(["samtools","faidx","-r",f"{WORK}/c.regions",REF],capture_output=True,text=True).stdout
blocks=[];cur=[]
for line in seq.split("\n"):
    if line.startswith(">"):
        if cur: blocks.append("".join(cur))
        cur=[]
    elif line: cur.append(line)
if cur: blocks.append("".join(cur))
comp=str.maketrans("ACGTacgtN","TGCAtgcaN")
with open(f"{WORK}/c.fa","w") as o:
    k=0
    for i in order:
        m=len(cds[i][4]); q="".join(blocks[k:k+m]); k+=m
        if cds[i][3]=="-": q=q.translate(comp)[::-1]
        o.write(f">L{i}\n{q}\n")
paf=subprocess.run(["minimap2","-x","asm20","-c","-N","200","-p","0.02","-t","4",
                    f"{WORK}/c.fa",f"{WORK}/c.fa"],capture_output=True,text=True).stdout
cdsedge=set()
for l in paf.split("\n"):
    if not l: continue
    g=l.split("\t"); q,ql,t,tl=g[0],int(g[1]),g[5],int(g[6])
    nres,blk=int(g[9]),int(g[10])
    if q==t or nres/max(blk,1)<ER_ID or blk<ER_COV*min(ql,tl): continue
    a,b=int(q[1:]),int(t[1:]); cdsedge.add((min(a,b),max(a,b)))

# ---- label same-family pairs and attach features ----
import itertools, statistics
rows=[]
for b,mem in fam.items():
    m=[i for i in mem if i in cds]
    if len(m)<2: continue
    n=len(mem)
    dens_e=sum(1 for u,v in set(edges) if u in set(mem) and v in set(mem))
    dens=2*dens_e/(n*(n-1)) if n>1 else 1.0
    for i,j in itertools.combinations(sorted(m),2):
        k=(i,j)
        d=direct.get(k)
        rows.append({"pair":k,"tp":k in cdsedge,"famsize":n,"dens":dens,
                     "direct":d is not None,"cov":d[0] if d else 0.0,"idn":d[1] if d else 0.0,
                     "blk":d[2] if d else 0,
                     "samechrom":loci[i][0]==loci[j][0],
                     "spanratio":min(loci[i][2]-loci[i][1],loci[j][2]-loci[j][1])/
                                 max(loci[i][2]-loci[i][1],loci[j][2]-loci[j][1]),
                     "minspan":min(loci[i][2]-loci[i][1],loci[j][2]-loci[j][1])})
tp=[r for r in rows if r["tp"]]; fp=[r for r in rows if not r["tp"]]
print(f"same-family pairs with CDS on both sides: {len(rows)}   concordant {len(tp)}   discordant {len(fp)}")
print(f"baseline precision: {len(tp)/len(rows):.4f}\n")
def med(v,k): 
    x=[r[k] for r in v]
    return statistics.median(x) if x else float('nan')
print(f"{'feature':<12}{'concordant':>13}{'discordant':>13}")
for k in ("famsize","dens","cov","idn","blk","spanratio","minspan"):
    print(f"{k:<12}{med(tp,k):>13.4f}{med(fp,k):>13.4f}")
for k in ("direct","samechrom"):
    print(f"{k:<12}{sum(r[k] for r in tp)/len(tp):>13.4f}{sum(r[k] for r in fp)/len(fp):>13.4f}")
