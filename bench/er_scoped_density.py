#!/usr/bin/env python3
"""Cliques / denser clusters, applied INSIDE families (ledger §6bm follow-up).

§6bm showed a global triangle or k-core requirement is a "families must have >=3 members" rule --
it cannot keep a 2-member family, and 72% of the catalog is 2-member. The scoped form is the only
viable one, so this sweeps the DENSITY demanded inside a family instead of an edge predicate:

  build families at the shipped gamma=0.20, then RE-REFINE each family of >= MINSZ members at a
  higher gamma'; families below MINSZ are left intact. gamma' = 1.0 demands a CLIQUE.

Both substrates are run in one pass, because §6bm's composite gained on one and vanished on the
other, and a threshold picked on a single proxy is not a result.
"""
import collections, itertools, os, pickle, re, subprocess, sys
BASE="/mnt/linuxdisk/home/juanfraitu"
GB=f"{BASE}/rustle_target/release/gamma_refine"
ID_FLOOR=0.60; COV_FLOOR=0.50; COVLONG=0.30; MINSZ=3

def load_paf(paf):
    best={}
    for l in open(paf):
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
        cs,cl=((qe-qs)/ql,(te-ts)/tl) if ql<=tl else ((te-ts)/tl,(qe-qs)/ql)
        k=(min(q,t),max(q,t))
        if k not in best or cs>best[k][0]: best[k]=(cs,cl,ident)
    return best
def gamma(el, g=None):
    env=dict(os.environ)
    if g is not None: env["RUSTLE_GAMMA"]=str(g)
    r=subprocess.run([GB],input="\n".join(f"{u}\t{v}" for u,v in sorted(el))+"\n",
                     capture_output=True,text=True,env=env)
    fam=collections.defaultdict(set)
    for line in r.stdout.split("\n"):
        if line:
            b,m=line.split("\t"); fam[b].add(int(m))
    return fam
def scoped(el, gprime, minsz=MINSZ):
    """gamma=0.20 families; each with >= minsz members re-refined at gprime, others untouched."""
    out={}; k=0
    for b,mem in gamma(el).items():
        sub={(u,v) for u,v in el if u in mem and v in mem}
        if len(mem)<minsz or not sub:
            out[f"F{k}"]=mem; k+=1; continue
        for b2,mem2 in gamma(sub, gprime).items():
            out[f"F{k}"]=mem2; k+=1
    return out
def run(name, dumpdir, prefix, copies, cds, cdsedge, znf, truth):
    best=load_paf(f"{dumpdir}/{prefix}.er._k11_w5.0.paf")
    SHIPPED={k for k,v in best.items() if v[2]>=ID_FLOOR and v[0]>=COV_FLOOR}
    COV={k for k in SHIPPED if best[k][1]>=COVLONG}
    byc=collections.defaultdict(list)
    for n,c in enumerate(copies): byc[c["chrom"]].append(n)
    def copy_at(c,a,b):
        bb=None;bov=0
        for n in byc.get(c,()):
            d=copies[n]; ov=min(b,d["end"])-max(a,d["start"])
            if ov>bov: bov,bb=ov,n
        return bb
    n2c={}
    with open(f"{dumpdir}/{prefix}.nodes.tsv") as fh:
        h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
        for l in fh:
            f=l.rstrip("\n").split("\t")
            n2c[int(f[ci["idx"]])]=copy_at(f[ci["chrom"]],int(f[ci["start"]]),int(f[ci["end"]]))
    NZPOS=sum(1 for p in cdsedge if not znf[p[0]] and not znf[p[1]])
    print(f"\n=== {name}  (E_r {len(SHIPPED)}, cov{COVLONG} {len(COV)}, non-ZNF positives {NZPOS}) ===")
    hdr=f"{'arm':<30}{'fams':>6}{'size2':>7}{'cop':>6}"
    if truth: hdr+=f"{'NPIP':>7}"
    print(hdr+f"{'nZ prec':>9}{'nZ rec':>8}{'nZ F1':>8}")
    arms=[("shipped", lambda: gamma(SHIPPED)), (f"cov_long >= {COVLONG}", lambda: gamma(COV))]
    for gp in (0.30,0.50,0.70,1.00):
        arms.append((f"cov + gamma'={gp:.2f} in fams>={MINSZ}", (lambda g: (lambda: scoped(COV,g)))(gp)))
    for lab,fn in arms:
        fam=fn(); cf={}
        for b,mem in fam.items():
            c={n2c[m] for m in mem if n2c.get(m) is not None}
            if len(c)>=2: cf[b]=sorted(c)
        tp=n=0
        for f_,mem in cf.items():
            m=[i for i in mem if i in cds]
            for p in itertools.combinations(sorted(m),2):
                if znf[p[0]] or znf[p[1]]: continue
                n+=1; tp += p in cdsedge
        pr=tp/n if n else 0; rc=tp/NZPOS if NZPOS else 0
        f1=2*pr*rc/(pr+rc) if pr+rc else 0
        row=(f"{lab:<30}{len(cf):>6}{sum(1 for v in cf.values() if len(v)==2):>7}"
             f"{sum(len(v) for v in cf.values()):>6}")
        if truth:
            npip=sum(1 for t in truth if any(copies[i]["chrom"]==t[0] and
                     min(t[2],copies[i]["end"])-max(t[1],copies[i]["start"])>0
                     for mem in cf.values() for i in mem))
            row+=f"{f'{npip}/31':>7}"
        print(row+f"{pr:>9.4f}{rc:>8.4f}{f1:>8.4f}")

def znf_map(copies, contigs=None):
    genes=collections.defaultdict(list)
    for l in open(f"{BASE}/winloci_data/GGO_genomic.gff"):
        if l[0]=="#": continue
        f=l.rstrip("\n").split("\t")
        if len(f)<9 or f[2] not in ("gene","pseudogene"): continue
        if contigs and f[0] not in contigs: continue
        d=re.search(r"(?:^|;)description=([^;]+)",f[8])
        genes[f[0]].append((int(f[3])-1,int(f[4]),(d.group(1).lower() if d else "")))
    return {n:any(min(c["end"],gb)-max(c["start"],ga)>0 and "zinc finger" in dd
                  for ga,gb,dd in genes.get(c["chrom"],())) for n,c in enumerate(copies)}

S=pickle.load(open(f"{BASE}/ecblob/state.pkl","rb"))
truth=[]
for l in open(f"{BASE}/o1_oracle/npip31.regions"):
    l=l.strip()
    if l:
        c,r=l.split(":"); a,b=r.split("-"); truth.append((c,int(a)-1,int(b)))
run("arm_f2 (3 contigs)", f"{BASE}/npip_cat/arm_f2/dump", "e", S["copies"], S["cds"], S["cdsedge"],
    znf_map(S["copies"], {"NC_073241.2","NC_073242.2","NC_073244.2"}), truth)
src=open("bench/fp_sources_read.py").read().split("fam=collections.defaultdict(list)\nfor n,c in enumerate(copies)")[0]
g={"__name__":"__main__"}; exec(compile(src,"r","exec"),g)
run("genome-wide fibroblast", f"{BASE}/o1_replicate/dump", "fibro",
    g["copies"], g["cds"], g["cdsedge"], znf_map(g["copies"]), None)
