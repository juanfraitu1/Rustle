#!/usr/bin/env python3
"""§6bk's two-sided coverage clause, re-run on a SECOND substrate.

arm_f2 is 3 contigs and its non-ZNF positive set was only 51 pairs. This repeats the arms on the
genome-wide fibroblast catalog (13,196 reps / 2,075 edges / 1,070 copies / 356 families), same
pipeline and same parameters (identity 0.60, coverage 0.50, -k11 -w5).

NPIP recall is NOT available here -- §6bh measured this catalog covering the truth set at 1/31 --
so the safety axis is catalog size and CDS recall instead.
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
node2copy={}; nodeinfo={}
with open(f"{DUMP}/fibro.nodes.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t"); i=int(f[ci["idx"]])
        node2copy[i]=copy_at(f[ci["chrom"]],int(f[ci["start"]]),int(f[ci["end"]]))
        nodeinfo[i]=int(f[ci["n_exon"]])
sys.stderr.write(f"nodes {len(node2copy)}  mapped to a copy {sum(1 for v in node2copy.values() if v is not None)}\n")
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
    if ql<=tl: cs,cl,sl,ll,al=(qe-qs)/ql,(te-ts)/tl,ql,tl,(qe-qs)
    else:      cs,cl,sl,ll,al=(te-ts)/tl,(qe-qs)/ql,tl,ql,(te-ts)
    k=(min(q,t),max(q,t))
    if k not in best or cs>best[k][0]: best[k]=(cs,cl,ident,sl,ll,al)
sh={k:v for k,v in best.items() if v[2]>=ID_FLOOR and v[0]>=COV_FLOOR}
print(f"shipped E_r edges reproduced: {len(sh)}   (params.tsv says 2075)")
print(f"median cov_longer over them: {statistics.median(v[1] for v in sh.values()):.4f}")
drop={k:v for k,v in sh.items() if v[1]<0.30}; keep={k:v for k,v in sh.items() if v[1]>=0.30}
print(f"\nMECHANISM CHECK (does the arm_f2 signature reproduce?)  dropped {len(drop)}  kept {len(keep)}")
print(f"{'':<28}{'DROPPED':>12}{'KEPT':>12}")
for lab,i in (("aligned bp on the shorter",5),("shorter length bp",3),
              ("longer length bp",4),("identity",2)):
    print(f"{lab:<28}{statistics.median(v[i] for v in drop.values()):>12,.3f}"
          f"{statistics.median(v[i] for v in keep.values()):>12,.3f}")
print(f"{'length ratio short/long':<28}"
      f"{statistics.median(v[3]/v[4] for v in drop.values()):>12.3f}"
      f"{statistics.median(v[3]/v[4] for v in keep.values()):>12.3f}")
for lab,fn in (("either end single-exon", lambda k: nodeinfo[k[0]]==1 or nodeinfo[k[1]]==1),):
    print(f"{lab:<28}{sum(1 for k in drop if fn(k))/len(drop):>12.3f}"
          f"{sum(1 for k in keep if fn(k))/len(keep):>12.3f}")
ALLPOS=len(cdsedge); NZPOS=sum(1 for p in cdsedge if not znf[p[0]] and not znf[p[1]])
print(f"\nCDS-homologous pairs: {ALLPOS} total, {NZPOS} non-ZNF")
def build(lf):
    el={(i,j) for (i,j),v in sh.items() if lf is None or v[1]>=lf}
    r=subprocess.run([f"{BASE}/rustle_target/release/gamma_refine"],
        input="\n".join(f"{u}\t{v}" for u,v in sorted(el))+"\n",capture_output=True,text=True)
    fam=collections.defaultdict(set)
    for line in r.stdout.split("\n"):
        if line:
            b,m=line.split("\t"); c=node2copy[int(m)]
            if c is not None: fam[b].add(c)
    return len(el), {k:sorted(v) for k,v in fam.items() if len(v)>=2}
print(f"\n{'longer-side floor':<22}{'edges':>7}{'fams':>6}{'cop':>6}"
      f"{'nZ prec':>9}{'nZ rec':>8}{'nZ F1':>8}{'all prec':>10}{'all rec':>9}")
for lab,lf in (("none (shipped)",None),("0.20",0.20),("0.30",0.30),("0.40",0.40),("0.50",0.50)):
    ne,fam=build(lf)
    res={}
    for tag,sel in (("all",lambda p:True),("nz",lambda p: not znf[p[0]] and not znf[p[1]])):
        tp=n=0
        for f_,mem in fam.items():
            m=[i for i in mem if i in cds]
            for p in itertools.combinations(sorted(m),2):
                if not sel(p): continue
                n+=1; tp += p in cdsedge
        res[tag]=(n,tp)
    na,tpa=res["all"]; nz,tpz=res["nz"]
    pz=tpz/nz if nz else 0; rz=tpz/NZPOS if NZPOS else 0
    f1z=2*pz*rz/(pz+rz) if pz+rz else 0
    print(f"{lab:<22}{ne:>7}{len(fam):>6}{sum(len(v) for v in fam.values()):>6}"
          f"{pz:>9.4f}{rz:>8.4f}{f1z:>8.4f}{(tpa/na if na else 0):>10.4f}{tpa/ALLPOS:>9.4f}")
