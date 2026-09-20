#!/usr/bin/env python3
"""Is the annotation-seeded catalog BETTER AS A CATALOG than the read-only one?

Independent criterion: neither catalog was built from coding sequence, so ask whether a
family's members share a CDS. The test is the project's OWN edge rule (identity >= 0.60,
coverage >= 0.50 of the shorter) applied to the longest annotated CDS at each member locus.
A real gene family should form ONE component under it; a repeat clique should not.

Reports the CANDIDATE COUNT (families with >=2 CDS-bearing members) before any rate, and a
SIZE-MATCHED null, because an edge-count-matched null proves nothing.
"""
import collections, os, random, subprocess, sys
BASE="/mnt/linuxdisk/home/juanfraitu"
WORK=f"{BASE}/cdscoh"; os.makedirs(WORK, exist_ok=True)
REF=f"{BASE}/npip_cat/npip3_contigs.fa"
GFF=f"{BASE}/winloci_data/GGO_genomic.gff"
CONTIGS={"NC_073241.2","NC_073242.2","NC_073244.2"}
MINCOV=0.30; ER_ID=0.60; ER_COV=0.50
random.seed(11)

# ---------------------------------------------------------------- CDS by gene
def load_cds():
    per=collections.defaultdict(lambda:{"gene":None,"chrom":None,"strand":None,"ex":[]})
    with open(GFF) as fh:
        for l in fh:
            if l[0]=="#": continue
            f=l.rstrip("\n").split("\t")
            if len(f)<9 or f[2]!="CDS" or f[0] not in CONTIGS: continue
            a=f[8]
            pid=None; gene=None
            for kv in a.split(";"):
                if kv.startswith("protein_id="): pid=kv[11:]
                elif kv.startswith("gene="): gene=kv[5:]
            if pid is None: continue
            d=per[pid]; d["gene"]=gene; d["chrom"]=f[0]; d["strand"]=f[6]
            d["ex"].append((int(f[3])-1,int(f[4])))
    out=[]
    for pid,d in per.items():
        d["ex"].sort(); ln=sum(b-a for a,b in d["ex"])
        out.append((d["chrom"],d["ex"][0][0],d["ex"][-1][1],d["strand"],ln,pid,d["ex"]))
    return out
CDS=load_cds()
sys.stderr.write(f"CDS transcripts on the 3 contigs: {len(CDS)}\n")
bychrom=collections.defaultdict(list)
for c in CDS: bychrom[c[0]].append(c)
for k in bychrom: bychrom[k].sort(key=lambda x:x[1])

def best_cds_in(chrom,a,b):
    """The longest CDS this locus POINTS AT (midpoint inside), taken WHOLE.

    Containment ("CDS fully inside the locus") was the first rule and it is BIASED: seeded
    loci are annotated gene intervals and contain a complete CDS by construction, while
    read loci are transcript spans (54.6% single-exon) that clip it. Midpoint-overlap asks
    only which gene a locus identifies, which both catalogs can answer on equal terms."""
    best=None
    for c in bychrom.get(chrom,()):
        mid=(c[1]+c[2])//2
        if a<=mid<b and (best is None or c[4]>best[4]): best=c
    return best

# ------------------------------------------------------------------ catalogs
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
    fam=collections.defaultdict(list)
    for i in range(len(loci)): fam[f(i)].append(i)
    return {f"SEED{k}":[tuple(loci[i]) for i in v] for k,v in fam.items() if len(v)>=2}

def read_families():
    fam=collections.defaultdict(list)
    with open(f"{BASE}/npip_cat/arm_f2/cat.copies.tsv") as fh:
        h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
        for l in fh:
            f=l.rstrip("\n").split("\t")
            fam[f[ci["family_id"]]].append((f[ci["chrom"]],int(f[ci["start"]]),int(f[ci["end"]])))
    return {k:v for k,v in fam.items() if len(v)>=2}

def null_families(template, pool_chrom_len):
    """size-matched null: same family SIZE distribution, spans drawn from the real ones"""
    out={}
    spans=[b-a for v in template.values() for c,a,b in v]
    for n,(k,v) in enumerate(template.items()):
        mem=[]
        for _ in v:
            c=random.choice(sorted(pool_chrom_len)); L=random.choice(spans)
            s=random.randint(0, max(pool_chrom_len[c]-L-1,1))
            mem.append((c,s,s+L))
        out[f"NULL{n}"]=mem
    return out

clen={}
for l in open(f"{REF}.fai"):
    f=l.split("\t"); clen[f[0]]=int(f[1])

cats={"seeded":seeded_families(), "read":read_families()}
import statistics as _st
for _c,_f in cats.items():
    _sp=[b-a for v in _f.values() for _,a,b in v]
    sys.stderr.write(f"  {_c}: {len(_f)} families, {len(_sp)} member loci, "
                     f"median span {_st.median(_sp):,.0f} bp\n")
cats["null"]=null_families(cats["seeded"], clen)

# ----------------------------------------------- one CDS per member, one FASTA
regions=[]; index=[]
for cat,fams in cats.items():
    for fid,mem in fams.items():
        for i,(c,a,b) in enumerate(mem):
            g=best_cds_in(c,a,b)
            if not g: continue
            name=f"{cat}|{fid}|{i}"
            index.append((name,g))
            for ea,eb in g[6]: regions.append(f"{c}:{ea+1}-{eb}")
with open(f"{WORK}/cds.regions","w") as o: o.write("\n".join(regions)+"\n")
seq=subprocess.run(["samtools","faidx","-r",f"{WORK}/cds.regions",REF],
                   capture_output=True,text=True).stdout
blocks=[]; cur=[]
for line in seq.split("\n"):
    if line.startswith(">"):
        if cur: blocks.append("".join(cur))
        cur=[]
    elif line: cur.append(line)
if cur: blocks.append("".join(cur))
comp=str.maketrans("ACGTacgtN","TGCAtgcaN")
with open(f"{WORK}/cds.fa","w") as o:
    k=0
    for name,g in index:
        n=len(g[6]); s="".join(blocks[k:k+n]); k+=n
        if g[3]=="-": s=s.translate(comp)[::-1]
        o.write(f">{name}\n{s}\n")
sys.stderr.write(f"CDS sequences written: {len(index)}  (blocks consumed {k}/{len(blocks)})\n")

# --------------------------------------------------- all-vs-all, then the E_r rule
paf=subprocess.run(["minimap2","-x","asm20","-c","-N","200","-p","0.02","-t","4",
                    f"{WORK}/cds.fa",f"{WORK}/cds.fa"],capture_output=True,text=True).stdout
qlen={}; edges=collections.defaultdict(set)
for l in paf.split("\n"):
    if not l: continue
    f=l.split("\t")
    q,ql,t,tl=f[0],int(f[1]),f[5],int(f[6])
    nres,blk=int(f[9]),int(f[10])
    qlen[q]=ql; qlen[t]=tl
    if q==t: continue
    if nres/max(blk,1) < ER_ID: continue
    if blk < ER_COV*min(ql,tl): continue
    edges[q].add(t); edges[t].add(q)

# --------------------------------------------------------------- score families
print(f"{'catalog':<9}{'families':>10}{'testable':>10}{'median largest':>16}{'fully':>9}{'singletons':>12}")
print(f"{'':<9}{'(>=2 mem)':>10}{'(>=2 CDS)':>10}{'CDS component':>16}{'coherent':>9}{'(no CDS edge)':>12}")
for cat,fams in cats.items():
    have=collections.defaultdict(list)
    for name,_ in index:
        c,fid,i=name.split("|")
        if c==cat: have[fid].append(name)
    testable={k:v for k,v in have.items() if len(v)>=2}
    nmem=sum(len(v) for v in fams.values()); ncds=sum(len(v) for v in have.values())
    fracs=[]; full=0; single=0
    for fid,names in testable.items():
        ns=set(names); seen=set(); best=0
        for n in names:
            if n in seen: continue
            stack=[n]; comp=set()
            while stack:
                x=stack.pop()
                if x in comp: continue
                comp.add(x)
                for y in edges.get(x,()):
                    if y in ns and y not in comp: stack.append(y)
            seen|=comp; best=max(best,len(comp))
        fracs.append(best/len(names))
        if best==len(names): full+=1
        if best==1: single+=1
    fracs.sort()
    med=fracs[len(fracs)//2] if fracs else 0
    print(f"{cat:<9}{len(fams):>10}{len(testable):>10}{med:>16.3f}"
          f"{f'{full}/{len(testable)}':>9}{f'{single}/{len(testable)}':>12}"
          f"   CDS-yield {ncds}/{nmem} = {ncds/max(nmem,1):.1%}")
