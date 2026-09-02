#!/usr/bin/env python3
"""Should E_r charge coverage on BOTH sequences (BLAST qcovs/scovs), not just the shorter?

The shipped clause is `cov = aligned span on the SHORTER / shorter length >= 0.50`
(`denovo_pipeline.rs:3110`), computed here VERBATIM including the `de:f:`-then-`nm/bl` identity
fallback. `cov_longer` is the same span measured on the LONGER sequence -- already DISCLOSED by
§6ba but never gated on.

Precedent for the idea inside the codebase: the divergent-PROTEIN tier already requires
`qcov/tcov >= min_coverage`, and its doc says the guard exists "to reject lone-shared-domain
merges" -- exactly the partial-homology failure mode.

Each arm rebuilds families from its own edge set (components + the shipped gamma binary) and is
scored on NPIP recall and CDS concordance. ⚠ offline reconstruction, so arms are comparable to
EACH OTHER, not to the shipped catalog.
"""
import collections, itertools, os, pickle, re, subprocess, sys
BASE="/mnt/linuxdisk/home/juanfraitu"; DUMP=f"{BASE}/npip_cat/arm_f2/dump"
S=pickle.load(open(f"{BASE}/ecblob/state.pkl","rb"))
copies,cds,cdsedge = S["copies"],S["cds"],S["cdsedge"]
ID_FLOOR=0.60; COV_FLOOR=0.50

nodes=[]
with open(f"{DUMP}/e.nodes.tsv") as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for l in fh:
        f=l.rstrip("\n").split("\t")
        nodes.append((f[ci["chrom"]],int(f[ci["start"]]),int(f[ci["end"]])))
byc=collections.defaultdict(list)
for n,c in enumerate(copies): byc[c["chrom"]].append(n)
def copy_at(c,a,b):
    best=None;bov=0
    for n in byc.get(c,()):
        d=copies[n]; ov=min(b,d["end"])-max(a,d["start"])
        if ov>bov: bov,best=ov,n
    return best
node2copy=[copy_at(*nd) for nd in nodes]

best=collections.defaultdict(lambda: None)      # (i,j) -> (ident, cov_short, cov_long)
with open(f"{DUMP}/e.er._k11_w5.0.paf") as fh:
    for l in fh:
        f=l.rstrip("\n").split("\t")
        if len(f)<12 or f[4]!="+": continue
        q,ql,qs,qe = int(f[0]),float(f[1]),float(f[2]),float(f[3])
        t,tl,ts,te = int(f[5]),float(f[6]),float(f[7]),float(f[8])
        nm,bl = float(f[9]),float(f[10])
        if bl<=0 or ql<=0 or tl<=0 or q==t: continue
        de=None
        for x in f[12:]:
            if x.startswith("de:f:"):
                try: de=float(x[5:])
                except ValueError: pass
        ident = (1.0-de) if de is not None else nm/bl
        if ql<=tl: cs,cln = (qe-qs)/ql, (te-ts)/tl
        else:      cs,cln = (te-ts)/tl, (qe-qs)/ql
        k=(min(q,t),max(q,t))
        cur=best[k]
        if cur is None or cs>cur[1]: best[k]=(ident,cs,cln)
sys.stderr.write(f"forward rep-pair records reduced to {len(best)} pairs\n")

genes=collections.defaultdict(list)
for l in open(f"{BASE}/winloci_data/GGO_genomic.gff"):
    if l[0]=="#": continue
    f=l.rstrip("\n").split("\t")
    if len(f)<9 or f[2] not in ("gene","pseudogene"): continue
    if f[0] not in ("NC_073241.2","NC_073242.2","NC_073244.2"): continue
    d=re.search(r"(?:^|;)description=([^;]+)",f[8])
    genes[f[0]].append((int(f[3])-1,int(f[4]),(d.group(1).lower() if d else "")))
znf={n:any(min(c["end"],gb)-max(c["start"],ga)>0 and "zinc finger" in dd
           for ga,gb,dd in genes.get(c["chrom"],())) for n,c in enumerate(copies)}
truth=[]
for l in open(f"{BASE}/o1_oracle/npip31.regions"):
    l=l.strip()
    if l:
        c,r=l.split(":"); a,b=r.split("-"); truth.append((c,int(a)-1,int(b)))

def build(long_floor):
    el=set()
    for (i,j),(ident,cs,cln) in best.items():
        if ident>=ID_FLOOR and cs>=COV_FLOOR and (long_floor is None or cln>=long_floor):
            el.add((i,j))
    r=subprocess.run([f"{BASE}/rustle_target/release/gamma_refine"],
        input="\n".join(f"{u}\t{v}" for u,v in sorted(el))+"\n",capture_output=True,text=True)
    fam=collections.defaultdict(set)
    for line in r.stdout.split("\n"):
        if line:
            b,m=line.split("\t")
            c=node2copy[int(m)]
            if c is not None: fam[b].add(c)
    return len(el), {k:sorted(v) for k,v in fam.items() if len(v)>=2}
def score(fam):
    npip=sum(1 for t in truth if any(copies[i]["chrom"]==t[0] and
             min(t[2],copies[i]["end"])-max(t[1],copies[i]["start"])>0
             for mem in fam.values() for i in mem))
    out={}
    for lab,sel in (("all",lambda p:True),
                    ("ZNF",lambda p: znf[p[0]] and znf[p[1]]),
                    ("nonZNF",lambda p: not znf[p[0]] and not znf[p[1]])):
        tp=n=0
        for f_,mem in fam.items():
            m=[i for i in mem if i in cds]
            for p in itertools.combinations(sorted(m),2):
                if not sel(p): continue
                n+=1; tp += p in cdsedge
        out[lab]=(n, tp/n if n else float("nan"), tp)
    return npip,out
ALLPOS=len(cdsedge)
NZPOS=sum(1 for p in cdsedge if not znf[p[0]] and not znf[p[1]])
print(f"CDS-homologous pairs: {ALLPOS} total, {NZPOS} non-ZNF\n")
print(f"{'coverage rule (short stays >=0.50)':<34}{'edges':>7}{'fams':>6}{'cop':>6}{'NPIP':>7}"
      f"{'nZ prec':>9}{'nZ rec':>8}{'nZ F1':>8}{'all F1':>8}")
for lab,lf in (("shipped: shorter only",None),
               ("both >= 0.30",0.30),("both >= 0.40",0.40),
               ("both >= 0.50 (symmetric)",0.50),("both >= 0.60",0.60)):
    ne,fam=build(lf); npip,sc=score(fam)
    npz,ppz,tpz=sc["nonZNF"]; na,pa,tpa=sc["all"]
    rz=tpz/NZPOS; f1z=2*ppz*rz/(ppz+rz) if ppz+rz else 0
    ra=tpa/ALLPOS; f1a=2*pa*ra/(pa+ra) if pa+ra else 0
    print(f"{lab:<34}{ne:>7}{len(fam):>6}{sum(len(v) for v in fam.values()):>6}"
          f"{f'{npip}/31':>7}{ppz:>9.4f}{rz:>8.4f}{f1z:>8.4f}{f1a:>8.4f}")
