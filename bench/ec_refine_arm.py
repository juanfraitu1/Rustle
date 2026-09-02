#!/usr/bin/env python3
"""E_c within-family refinement, evaluated as a catalog change (ledger §6bi follow-up).

E_c reaches only ~13% of same-family pairs, so splitting EVERY family by its E_c components
would shatter the ones E_c cannot see. The rule therefore carries a GUARD: a family is refined
only when it is large enough AND carries enough E_c evidence; otherwise it is left intact.

Reports, for each guard setting: catalog size, NPIP recall (the safety check), and CDS
concordance split by ZNF / non-ZNF -- because §6bi showed the instrument fails on ZNF, so the
non-ZNF column is where generalisation is actually tested.
"""
import collections, itertools, pickle, re, sys
BASE="/mnt/linuxdisk/home/juanfraitu"
S=pickle.load(open(f"{BASE}/ecblob/state.pkl","rb"))
copies,direct,ec,cds,cdsedge = S["copies"],S["direct"],S["ec"],S["cds"],S["cdsedge"]
genes=collections.defaultdict(list)
for l in open(f"{BASE}/winloci_data/GGO_genomic.gff"):
    if l[0]=="#": continue
    f=l.rstrip("\n").split("\t")
    if len(f)<9 or f[2] not in ("gene","pseudogene"): continue
    if f[0] not in ("NC_073241.2","NC_073242.2","NC_073244.2"): continue
    d=re.search(r"(?:^|;)description=([^;]+)",f[8])
    genes[f[0]].append((int(f[3])-1,int(f[4]),(d.group(1).lower() if d else "")))
znf={}
for n,c in enumerate(copies):
    znf[n]=any(min(c["end"],gb)-max(c["start"],ga)>0 and "zinc finger" in dd
               for ga,gb,dd in genes.get(c["chrom"],()))
fam0=collections.defaultdict(list)
for n,c in enumerate(copies): fam0[c["fam"]].append(n)
truth=[]
for l in open(f"{BASE}/o1_oracle/npip31.regions"):
    l=l.strip()
    if l:
        c,r=l.split(":"); a,b=r.split("-"); truth.append((c,int(a)-1,int(b)))
def npip(fams):
    hit=0
    for t in truth:
        if any(copies[i]["chrom"]==t[0] and min(t[2],copies[i]["end"])-max(t[1],copies[i]["start"])>0
               for mem in fams.values() for i in mem): hit+=1
    return hit
def ec_components(mem, minreads):
    adj=collections.defaultdict(set)
    for i,j in itertools.combinations(sorted(mem),2):
        if ec.get((i,j),0)>=minreads: adj[i].add(j); adj[j].add(i)
    seen=set(); out=[]
    for s in mem:
        if s in seen: continue
        st=[s]; c=set()
        while st:
            x=st.pop()
            if x in seen: continue
            seen.add(x); c.add(x); st.extend(y for y in adj[x] if y not in seen)
        out.append(c)
    return out
def refine(min_size, min_frac, minreads=1):
    out={}; nsplit=0
    for f_,mem in fam0.items():
        n=len(mem)
        if n<2: continue
        npairs=n*(n-1)//2
        nec=sum(1 for p in itertools.combinations(sorted(mem),2) if ec.get(p,0)>=minreads)
        if n>=min_size and npairs and nec/npairs>=min_frac:
            comps=[c for c in ec_components(mem,minreads) if len(c)>=2]
            if len(comps)>1: nsplit+=1
            for k,c in enumerate(comps): out[f"{f_}.{k}"]=sorted(c)
        else:
            out[f_]=sorted(mem)
    return out,nsplit
def score(fams, sel):
    tp=n=0
    for f_,mem in fams.items():
        m=[i for i in mem if i in cds]
        for p in itertools.combinations(sorted(m),2):
            if not sel(p): continue
            n+=1; tp += p in cdsedge
    return n, (tp/n if n else float("nan"))
ALL=lambda p: True
NZ=lambda p: not znf[p[0]] and not znf[p[1]]
ZZ=lambda p: znf[p[0]] and znf[p[1]]
print(f"arm_f2 baseline: {len(fam0)} families, {sum(len(v) for v in fam0.values())} copies")
print(f"{'guard (min size, min E_c frac)':<32}{'fams':>6}{'copies':>8}{'split':>7}{'NPIP':>7}"
      f"{'prec all':>10}{'prec ZNF':>10}{'prec nonZNF':>13}")
base,_=refine(10**9,2.0)
for lab,args in (("none (shipped)",(10**9,2.0)),
                 ("size>=10, E_c frac>=0.02",(10,0.02)),
                 ("size>=10, E_c frac>=0.05",(10,0.05)),
                 ("size>=20, E_c frac>=0.05",(20,0.05)),
                 ("size>=5,  E_c frac>=0.05",(5,0.05)),
                 ("size>=2,  E_c frac>=0.00",(2,0.0))):
    f,ns=refine(*args)
    na,pa=score(f,ALL); nz,pz=score(f,ZZ); nn,pn=score(f,NZ)
    print(f"{lab:<32}{len(f):>6}{sum(len(v) for v in f.values()):>8}{ns:>7}"
          f"{f'{npip(f)}/31':>7}{pa:>10.4f}{pz:>10.4f}{pn:>13.4f}")
