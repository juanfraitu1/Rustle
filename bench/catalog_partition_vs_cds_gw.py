#!/usr/bin/env python3
"""Genome-wide form of bench/catalog_partition_vs_cds.py.

Three contigs gave 40 shared loci / 116 CDS-homologous pairs, half of them inside ONE seeded
component -- underpowered and blob-driven. This repeats the identical-loci test on the whole
assembly: 41,193 annotated seeds vs the genome-wide fibroblast catalog (1,070 copies).
Only the GROUPING varies; the locus set and its CDS are shared.
"""
import bisect, collections, itertools, os, subprocess, sys
BASE="/mnt/linuxdisk/home/juanfraitu"
REF=f"{BASE}/_from_wsl/winloci_scratch/GGO.fasta"
GFF=f"{BASE}/winloci_data/GGO_genomic.gff"
WORK=f"{BASE}/cdscoh_gw"; os.makedirs(WORK,exist_ok=True)
MINCOV=0.30; ER_ID=0.60; ER_COV=0.50

# ------------------------------------------------------------------- CDS
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
    mid=(d["ex"][0][0]+d["ex"][-1][1])//2
    bychrom[d["chrom"]].append((mid,ln,d["strand"],d["ex"]))
for k in bychrom: bychrom[k].sort()
mids={k:[x[0] for x in v] for k,v in bychrom.items()}
sys.stderr.write(f"CDS transcripts genome-wide: {len(per)}\n")
def best_cds_in(chrom,a,b):
    v=bychrom.get(chrom)
    if not v: return None
    lo=bisect.bisect_left(mids[chrom],a); hi=bisect.bisect_right(mids[chrom],b)
    best=None
    for i in range(lo,hi):
        if best is None or v[i][1]>best[1]: best=v[i]
    return (chrom,)+best if best else None

# --------------------------------------------------------------- seeded
def pr(s):
    c,r=s.rsplit(":",1); a,b=r.split("-"); return c,int(a)-1,int(b)
pairs=[]
for l in open(f"{BASE}/seedgw/multiseed_gw.paf"):
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
starts={k:[loci[i][1] for i in v] for k,v in byc.items()}
def lof(iv):
    c,a,b=iv; v=byc.get(c)
    if not v: return None
    j=bisect.bisect_right(starts[c],b)-1
    while j>=0:
        i=v[j]; L=loci[i]
        if L[2]<=a and loci[v[max(j-1,0)]][2]<=a and j<len(v)-1: pass
        if min(b,L[2])-max(a,L[1])>0: return i
        if L[1]<a-5_000_000: break
        j-=1
    return None
# GAMMA: the seeded relation is refined by the SHIPPED gamma-quasi-clique rule, not by raw
# connected components. Without this the comparison varies the PARTITIONER as well as the edges
# -- the read catalog is gamma-refined and the seeded one was not, which is why it carried a
# 191-member blob. `gamma_refine` calls family_definition::refine_component directly.
GB=f"{BASE}/rustle_target/release/gamma_refine"
el=[]
for s_,t_ in pairs:
    i,j=lof(s_),lof(t_)
    if i is not None and j is not None and i!=j: el.append(f"{i}\t{j}")
GAMMA=os.environ.get("RUSTLE_GAMMA","")
env=dict(os.environ)
r=subprocess.run([GB],input="\n".join(el)+"\n",capture_output=True,text=True,env=env)
sys.stderr.write(r.stderr)
block={}
for line in r.stdout.split("\n"):
    if not line: continue
    b,m=line.split("\t"); block[int(m)]=b
def f(x): return block.get(x, f"S{x}")
sys.stderr.write(f"seeded loci: {len(loci)}  gamma-refined blocks: {len(set(block.values()))}\n")

# ----------------------------------------------------------- read catalog
readcop=collections.defaultdict(list)
with open(f"{BASE}/o1_replicate/fibro_gwcat.copies.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        g=l.rstrip("\n").split("\t")
        readcop[g[ci["chrom"]]].append((int(g[ci["start"]]),int(g[ci["end"]]),g[ci["family_id"]]))
for k in readcop: readcop[k].sort()
nread=sum(len(v) for v in readcop.values())
sys.stderr.write(f"read-catalog copies: {nread}\n")

shared=[]
for i,L in enumerate(loci):
    hits=[r for r in readcop.get(L[0],()) if min(L[2],r[1])-max(L[1],r[0])>0]
    if not hits: continue
    r=max(hits,key=lambda x:min(L[2],x[1])-max(L[1],x[0]))
    g=best_cds_in(L[0],L[1],L[2])
    if not g: continue
    shared.append({"iv":tuple(L),"seed":f"SEED{f(i)}","read":r[2],"cds":g})
print(f"shared loci with a CDS: {len(shared)}")
if len(shared)<2: sys.exit(0)

regions=[]
for s in shared:
    for ea,eb in s["cds"][4]: regions.append(f"{s['cds'][0]}:{ea+1}-{eb}")
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
        m=len(s["cds"][4]); q="".join(blocks[k:k+m]); k+=m
        if s["cds"][3]=="-": q=q.translate(comp)[::-1]
        o.write(f">L{n}\n{q}\n")
sys.stderr.write(f"blocks consumed {k}/{len(blocks)}\n")
paf=subprocess.run(["minimap2","-x","asm20","-c","-N","200","-p","0.02","-t","4",
                    f"{WORK}/shared.fa",f"{WORK}/shared.fa"],capture_output=True,text=True).stdout
cdsedge=set()
for l in paf.split("\n"):
    if not l: continue
    g=l.split("\t"); q,ql,t,tl=g[0],int(g[1]),g[5],int(g[6])
    nres,blk=int(g[9]),int(g[10])
    if q==t or nres/max(blk,1)<ER_ID or blk<ER_COV*min(ql,tl): continue
    a,b=int(q[1:]),int(t[1:]); cdsedge.add((min(a,b),max(a,b)))

def score(idx,key):
    same={(i,j) for i,j in itertools.combinations(sorted(idx),2) if shared[i][key]==shared[j][key]}
    pos={(i,j) for i,j in cdsedge if i in idx and j in idx}
    tp=len(same&pos)
    rec=tp/len(pos) if pos else float("nan")
    prec=tp/len(same) if same else float("nan")
    f1=2*rec*prec/(rec+prec) if (rec and prec and rec+prec) else 0
    return len(pos),len(same),rec,prec,f1
allidx=set(range(len(shared)))
big=collections.Counter(s["seed"] for s in shared).most_common(1)[0]
noblob={i for i in allidx if shared[i]["seed"]!=big[0]}
for label,idx in (("ALL shared loci",allidx),(f"largest seeded block ({big[1]} loci) EXCLUDED",noblob)):
    pos,_,_,_,_=score(idx,"seed")
    print(f"\n-- {label}: {len(idx)} loci, {len(idx)*(len(idx)-1)//2} pairs, "
          f"{pos} CDS-homologous")
    print(f"     {'grouping':<8}{'same-family':>12}{'recall':>9}{'precision':>11}{'F1':>8}")
    for key in ("seed","read"):
        pos,same,rec,prec,f1=score(idx,key)
        print(f"     {key:<8}{same:>12}{rec:>9.4f}{prec:>11.4f}{f1:>8.4f}")
