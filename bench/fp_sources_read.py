#!/usr/bin/env python3
"""The false-merge analysis of bench/fp_sources_seeded.py, run on the SHIPPED read catalog.

If the dominant FP source (membership carried by transitive closure rather than a direct E_r
edge) is shared, it is a property of defining families as connected components -- not of the
seeded mode, and it applies to what ships.
"""
import bisect, collections, itertools, os, subprocess, sys
BASE="/mnt/linuxdisk/home/juanfraitu"
REF=f"{BASE}/_from_wsl/winloci_scratch/GGO.fasta"; GFF=f"{BASE}/winloci_data/GGO_genomic.gff"
WORK=f"{BASE}/fpsrc"; os.makedirs(WORK,exist_ok=True)
ER_ID=0.60; ER_COV=0.50

per=collections.defaultdict(lambda:{"chrom":None,"strand":None,"ex":[]})
for l in open(GFF):
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
    lo=bisect.bisect_left(mids[c],a); hi=bisect.bisect_right(mids[c],b); best=None
    for i in range(lo,hi):
        if best is None or v[i][1]>best[1]: best=v[i]
    return (c,)+best if best else None

copies=[]
with open(f"{BASE}/o1_replicate/fibro_gwcat.copies.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t")
        copies.append({"fam":f[ci["family_id"]],"chrom":f[ci["chrom"]],
                       "start":int(f[ci["start"]]),"end":int(f[ci["end"]])})
sys.stderr.write(f"copies {len(copies)}  families {len({c['fam'] for c in copies})}\n")
ered=[]
with open(f"{BASE}/o1_replicate/dump/fibro.edges.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t")
        ered.append(((f[ci["chrom_i"]],int(f[ci["start_i"]]),int(f[ci["end_i"]])),
                     (f[ci["chrom_j"]],int(f[ci["start_j"]]),int(f[ci["end_j"]]))))
sys.stderr.write(f"E_r edges {len(ered)}\n")
bych=collections.defaultdict(list)
for n,c in enumerate(copies): bych[c["chrom"]].append(n)
def copy_at(iv):
    c,a,b=iv; best=None; bov=0
    for n in bych.get(c,()):
        d=copies[n]; ov=min(b,d["end"])-max(a,d["start"])
        if ov>bov: bov,best=ov,n
    return best
direct=set()
for u,v in ered:
    i,j=copy_at(u),copy_at(v)
    if i is not None and j is not None and i!=j: direct.add((min(i,j),max(i,j)))
sys.stderr.write(f"E_r edges mapped onto copy pairs: {len(direct)}\n")
cds={}
for n,c in enumerate(copies):
    g=best_cds_in(c["chrom"],c["start"],c["end"])
    if g: cds[n]=g
sys.stderr.write(f"copies with a CDS: {len(cds)}\n")
order=sorted(cds)
open(f"{WORK}/r.regions","w").write("\n".join(
    f"{cds[i][0]}:{ea+1}-{eb}" for i in order for ea,eb in cds[i][4])+"\n")
seq=subprocess.run(["samtools","faidx","-r",f"{WORK}/r.regions",REF],capture_output=True,text=True).stdout
blocks=[];cur=[]
for line in seq.split("\n"):
    if line.startswith(">"):
        if cur: blocks.append("".join(cur))
        cur=[]
    elif line: cur.append(line)
if cur: blocks.append("".join(cur))
comp=str.maketrans("ACGTacgtN","TGCAtgcaN")
with open(f"{WORK}/r.fa","w") as o:
    k=0
    for i in order:
        m=len(cds[i][4]); q="".join(blocks[k:k+m]); k+=m
        if cds[i][3]=="-": q=q.translate(comp)[::-1]
        o.write(f">L{i}\n{q}\n")
paf=subprocess.run(["minimap2","-x","asm20","-c","-N","200","-p","0.02","-t","4",
                    f"{WORK}/r.fa",f"{WORK}/r.fa"],capture_output=True,text=True).stdout
cdsedge=set()
for l in paf.split("\n"):
    if not l: continue
    g=l.split("\t"); q,ql,t,tl=g[0],int(g[1]),g[5],int(g[6])
    nres,blk=int(g[9]),int(g[10])
    if q==t or nres/max(blk,1)<ER_ID or blk<ER_COV*min(ql,tl): continue
    a,b=int(q[1:]),int(t[1:]); cdsedge.add((min(a,b),max(a,b)))
fam=collections.defaultdict(list)
for n,c in enumerate(copies): fam[c["fam"]].append(n)
rows=[]
for f_,mem in fam.items():
    if len(mem)<2: continue
    m=[i for i in mem if i in cds]
    for i,j in itertools.combinations(sorted(m),2):
        rows.append({"tp":(i,j) in cdsedge,"direct":(i,j) in direct})
tp=[r for r in rows if r["tp"]]; fp=[r for r in rows if not r["tp"]]
print(f"READ CATALOG: same-family pairs with CDS both sides {len(rows)}"
      f"  concordant {len(tp)}  discordant {len(fp)}  precision {len(tp)/max(len(rows),1):.4f}")
print(f"  directly E_r-linked -- concordant {sum(r['direct'] for r in tp)/max(len(tp),1):.4f}"
      f"  discordant {sum(r['direct'] for r in fp)/max(len(fp),1):.4f}")
print(f"\n{'family size':<14}{'pairs':>8}{'prec|direct':>13}{'prec|transitive':>17}{'lift':>8}")
for lo,hi in ((2,4),(5,10),(11,30),(31,10**9)):
    a=b=c=d=0
    for f_,mem in fam.items():
        if not (lo<=len(mem)<=hi): continue
        m=[i for i in mem if i in cds]
        for i,j in itertools.combinations(sorted(m),2):
            t=(i,j) in cdsedge
            if (i,j) in direct: a+=t; b+=1
            else: c+=t; d+=1
    lab=f"{lo}-{hi}" if hi<10**9 else f"{lo}+"
    print(f"{lab:<14}{b+d:>8}{(a/b if b else float('nan')):>13.4f}"
          f"{(c/d if d else float('nan')):>17.4f}{((a/b-c/d) if b and d else float('nan')):>8.4f}")
kept=[r for r in rows if r["direct"]]
if kept and cdsedge:
    print(f"\nfilter 'direct E_r edge required':  precision {len(tp)/len(rows):.4f} -> "
          f"{sum(r['tp'] for r in kept)/len(kept):.4f}   pairs {len(rows)} -> {len(kept)}")
