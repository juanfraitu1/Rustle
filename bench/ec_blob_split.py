#!/usr/bin/env python3
"""Can SHARED MULTI-MAPPING READS (E_c) split the blobs? (the advisor's edge idea, scoped)

E_c was refuted as a DEFINITION (§6bb: 18% coverage, blind below identity 0.95) and as a
conjunction term. This asks something different and much narrower: inside an over-merged family,
does E_c mark the genuinely-together members? Its high-identity specificity, fatal for a
definition, is what a blob splitter wants.

Substrate: arm_f2 (3 contigs) -- the only catalog with E_c computed on the same rep set. Its
largest family GWFAM79 (54 members) supplies 34.3% of all same-family pairs, the same blob shape
as GWFAM2 genome-wide.
Truth proxy: CDS concordance under E_r (see §6bg for why it is RELATIVE, not a rate).
"""
import bisect, collections, itertools, os, subprocess, sys
BASE="/mnt/linuxdisk/home/juanfraitu"
REF=f"{BASE}/npip_cat/npip3_contigs.fa"; GFF=f"{BASE}/winloci_data/GGO_genomic.gff"
DUMP=f"{BASE}/npip_cat/arm_f2/dump"; WORK=f"{BASE}/ecblob"; os.makedirs(WORK,exist_ok=True)
CONTIGS={"NC_073241.2","NC_073242.2","NC_073244.2"}
ER_ID=0.60; ER_COV=0.50

per=collections.defaultdict(lambda:{"chrom":None,"strand":None,"ex":[]})
for l in open(GFF):
    if l[0]=="#": continue
    f=l.rstrip("\n").split("\t")
    if len(f)<9 or f[2]!="CDS" or f[0] not in CONTIGS: continue
    pid=None
    for kv in f[8].split(";"):
        if kv.startswith("protein_id="): pid=kv[11:]; break
    if pid is None: continue
    d=per[pid]; d["chrom"]=f[0]; d["strand"]=f[6]; d["ex"].append((int(f[3])-1,int(f[4])))
bych=collections.defaultdict(list)
for pid,d in per.items():
    d["ex"].sort(); ln=sum(b-a for a,b in d["ex"])
    bych[d["chrom"]].append(((d["ex"][0][0]+d["ex"][-1][1])//2, ln, d["strand"], d["ex"]))
for k in bych: bych[k].sort()
mids={k:[x[0] for x in v] for k,v in bych.items()}
def best_cds_in(c,a,b):
    v=bych.get(c)
    if not v: return None
    lo=bisect.bisect_left(mids[c],a); hi=bisect.bisect_right(mids[c],b); r=None
    for i in range(lo,hi):
        if r is None or v[i][1]>r[1]: r=v[i]
    return (c,)+r if r else None

copies=[]
with open(f"{BASE}/npip_cat/arm_f2/cat.copies.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t")
        copies.append({"fam":f[ci["family_id"]],"chrom":f[ci["chrom"]],
                       "start":int(f[ci["start"]]),"end":int(f[ci["end"]])})
byc=collections.defaultdict(list)
for n,c in enumerate(copies): byc[c["chrom"]].append(n)
def copy_at(c,a,b):
    best=None; bov=0
    for n in byc.get(c,()):
        d=copies[n]; ov=min(b,d["end"])-max(a,d["start"])
        if ov>bov: bov,best=ov,n
    return best
direct=set(); nodemap={}
with open(f"{DUMP}/e.edges.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t")
        i=copy_at(f[ci["chrom_i"]],int(f[ci["start_i"]]),int(f[ci["end_i"]]))
        j=copy_at(f[ci["chrom_j"]],int(f[ci["start_j"]]),int(f[ci["end_j"]]))
        nodemap[f[ci["node_key_i"]]]=i; nodemap[f[ci["node_key_j"]]]=j
        if i is not None and j is not None and i!=j: direct.add((min(i,j),max(i,j)))
with open(f"{DUMP}/e.nodes.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t")
        nodemap.setdefault(f[ci["node_key"]],
                           copy_at(f[ci["chrom"]],int(f[ci["start"]]),int(f[ci["end"]])))
ec=collections.defaultdict(int)
p=f"{BASE}/ecband/allpairs_ec.tsv"
if os.path.exists(p):
    with open(p) as fh:
        fh.readline()
        for l in fh:
            a,b,n0,n50=l.rstrip("\n").split("\t")
            i,j=nodemap.get(a),nodemap.get(b)
            if i is None or j is None or i==j: continue
            k=(min(i,j),max(i,j)); ec[k]=max(ec[k],int(n0))
sys.stderr.write(f"copies {len(copies)}  direct E_r pairs {len(direct)}  E_c copy-pairs {len(ec)}\n")

cds={}
for n,c in enumerate(copies):
    g=best_cds_in(c["chrom"],c["start"],c["end"])
    if g: cds[n]=g
order=sorted(cds)
open(f"{WORK}/c.regions","w").write("\n".join(
    f"{cds[i][0]}:{ea+1}-{eb}" for i in order for ea,eb in cds[i][4])+"\n")
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
sys.stderr.write(f"copies with a CDS {len(cds)}   CDS-homologous pairs {len(cdsedge)}\n")
import pickle
pickle.dump({"copies":copies,"direct":direct,"ec":dict(ec),"cds":cds,"cdsedge":cdsedge},
            open(f"{WORK}/state.pkl","wb"))
print("state written")
