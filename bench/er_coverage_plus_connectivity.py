#!/usr/bin/env python3
"""Does connectivity add anything ON TOP of the two-sided coverage clause?

§6bh found connectivity subsumed by the direct-edge rule. This is a different baseline: the
coverage clause (§6bk/§6bl) changes WHICH EDGES EXIST, so the remaining graph is a different
object and the features may behave differently -- or may now be redundant with it.

Arms are edge-set filters, then the shipped gamma binary, then scored on the non-ZNF column.
Run on the genome-wide fibroblast catalog (167 non-ZNF positive pairs).
⚠ many combinations are being scored against ONE proxy, so anything that wins here must also be
checked on the second substrate before it means anything.
"""
import collections, itertools, os, re, subprocess, statistics, sys
BASE="/mnt/linuxdisk/home/juanfraitu"; DUMP=f"{BASE}/o1_replicate/dump"
ID_FLOOR=0.60; COV_FLOOR=0.50
src=open("bench/fp_sources_read.py").read().split("fam=collections.defaultdict(list)\nfor n,c in enumerate(copies)")[0]
g={"__name__":"__main__"}; exec(compile(src,"r","exec"),g)
copies,cds,cdsedge = g["copies"],g["cds"],g["cdsedge"]
genes=collections.defaultdict(list)
for l in open(f"{BASE}/winloci_data/GGO_genomic.gff"):
    if l[0]=="#": continue
    f=l.rstrip("\n").split("\t")
    if len(f)<9 or f[2] not in ("gene","pseudogene"): continue
    d=re.search(r"(?:^|;)description=([^;]+)",f[8])
    genes[f[0]].append((int(f[3])-1,int(f[4]),(d.group(1).lower() if d else "")))
znf={n:any(min(c["end"],gb)-max(c["start"],ga)>0 and "zinc finger" in dd
           for ga,gb,dd in genes.get(c["chrom"],())) for n,c in enumerate(copies)}
byc=collections.defaultdict(list)
for n,c in enumerate(copies): byc[c["chrom"]].append(n)
def copy_at(c,a,b):
    best=None;bov=0
    for n in byc.get(c,()):
        d=copies[n]; ov=min(b,d["end"])-max(a,d["start"])
        if ov>bov: bov,best=ov,n
    return best
node2copy={}
with open(f"{DUMP}/fibro.nodes.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t")
        node2copy[int(f[ci["idx"]])]=copy_at(f[ci["chrom"]],int(f[ci["start"]]),int(f[ci["end"]]))
best={}
for l in open(f"{DUMP}/fibro.er._k11_w5.0.paf"):
    f=l.rstrip("\n").split("\t")
    if len(f)<12 or f[4]!="+": continue
    q,ql,qs,qe=int(f[0]),float(f[1]),float(f[2]),float(f[3])
    t,tl,ts,te=int(f[5]),float(f[6]),float(f[7]),float(f[8])
    nm,bl=float(f[9]),float(f[10])
    if bl<=0 or ql<=0 or tl<=0 or q==t: continue
    de=None
    for x in f[12:]:
        if x.startswith("de:f:"):
            try: de=float(x[5:])
            except ValueError: pass
    ident=(1.0-de) if de is not None else nm/bl
    cs,cl = ((qe-qs)/ql,(te-ts)/tl) if ql<=tl else ((te-ts)/tl,(qe-qs)/ql)
    k=(min(q,t),max(q,t))
    if k not in best or cs>best[k][0]: best[k]=(cs,cl,ident)
SHIPPED={k for k,v in best.items() if v[2]>=ID_FLOOR and v[0]>=COV_FLOOR}
ALLPOS=len(cdsedge); NZPOS=sum(1 for p in cdsedge if not znf[p[0]] and not znf[p[1]])

def triangles(el):
    adj=collections.defaultdict(set)
    for u,v in el: adj[u].add(v); adj[v].add(u)
    return {(u,v) for u,v in el if adj[u]&adj[v]}
def kcore(el,k):
    adj=collections.defaultdict(set)
    for u,v in el: adj[u].add(v); adj[v].add(u)
    alive={n for n in adj}
    while True:
        drop={n for n in alive if len(adj[n]&alive)<k}
        if not drop: break
        alive-=drop
    return {(u,v) for u,v in el if u in alive and v in alive}
def build(el):
    r=subprocess.run([f"{BASE}/rustle_target/release/gamma_refine"],
        input="\n".join(f"{u}\t{v}" for u,v in sorted(el))+"\n",capture_output=True,text=True)
    fam=collections.defaultdict(set)
    for line in r.stdout.split("\n"):
        if line:
            b,m=line.split("\t"); c=node2copy[int(m)]
            if c is not None: fam[b].add(c)
    return {k:sorted(v) for k,v in fam.items() if len(v)>=2}
def score(fam):
    out={}
    for tag,sel in (("all",lambda p:True),("nz",lambda p: not znf[p[0]] and not znf[p[1]])):
        tp=n=0
        for f_,mem in fam.items():
            m=[i for i in mem if i in cds]
            for p in itertools.combinations(sorted(m),2):
                if not sel(p): continue
                n+=1; tp += p in cdsedge
        out[tag]=(n,tp)
    return out
COV30={k for k in SHIPPED if best[k][1]>=0.30}
arms=[("shipped", SHIPPED),
      ("cov_long >= 0.30", COV30),
      ("triangle >= 1 only", triangles(SHIPPED)),
      ("k-core 2 only", kcore(SHIPPED,2)),
      ("cov 0.30 + triangle", triangles(COV30)),
      ("cov 0.30 + k-core 2", kcore(COV30,2)),
      ("cov 0.30 + tri + kcore2", kcore(triangles(COV30),2))]
print(f"non-ZNF positives {NZPOS}   all positives {ALLPOS}\n")
print(f"{'arm':<26}{'edges':>7}{'fams':>6}{'cop':>6}{'nZ prec':>9}{'nZ rec':>8}{'nZ F1':>8}")
for lab,el in arms:
    fam=build(el); sc=score(fam)
    nz,tpz=sc["nz"]
    p=tpz/nz if nz else 0; r=tpz/NZPOS; f1=2*p*r/(p+r) if p+r else 0
    print(f"{lab:<26}{len(el):>7}{len(fam):>6}{sum(len(v) for v in fam.values()):>6}"
          f"{p:>9.4f}{r:>8.4f}{f1:>8.4f}")
