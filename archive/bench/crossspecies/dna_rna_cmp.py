#!/usr/bin/env python3
"""DNA vs RNA E_r on identical node sets, shipped two-tier rule, aggregated coverage."""
from collections import defaultdict
W='/home/juanfra/winloci_scratch/seedfam/dnapr/meeting'
def ulen(iv):
    iv=sorted(iv); tot=0; cs=ce=None
    for s,e in iv:
        if cs is None: cs,ce=s,e
        elif s>ce: tot+=ce-cs; cs,ce=s,e
        else: ce=max(ce,e)
    return tot+(ce-cs if cs is not None else 0)
def cov(p,floor,ns):
    recs=defaultdict(list)
    for line in open(p):
        f=line.rstrip('\n').split('\t')
        if len(f)<12: continue
        q,ql,qs,qe,t,tl=f[0],int(f[1]),int(f[2]),int(f[3]),f[5],int(f[6])
        nm,bl=int(f[9]),int(f[10])
        if q==t or bl==0 or q not in ns or t not in ns or nm/bl<floor: continue
        recs[(q,t)].append((qs,qe,ql,tl))
    m={}
    for (q,t),v in recs.items():
        c=ulen([(x[0],x[1]) for x in v])/max(min(v[0][2],v[0][3]),1)
        k=frozenset((q,t)); m[k]=max(m.get(k,0),c)
    return m
def E(pre,ns):
    a=cov(f'{pre}_a.paf',0.80,ns); b=cov(f'{pre}_b.paf',0.60,ns)
    m={}
    for d in (a,b):
        for k,c in d.items(): m[k]=max(m.get(k,0),c)
    return {k for k,c in m.items() if c>=0.50}
def comps(ns,ed):
    p={x:x for x in ns}
    def fi(x):
        while p[x]!=x: p[x]=p[p[x]]; x=p[x]
        return x
    for e in ed:
        a,b=tuple(e); ra,rb=fi(a),fi(b)
        if ra!=rb: p[ra]=rb
    g=defaultdict(list)
    for x in ns: g[fi(x)].append(x)
    return sorted(g.values(),key=len,reverse=True)
print(f"\n{'family':<9}{'|V|':>5}{'rna|V|':>8}{'DNA E':>7}{'RNA E':>7}{'shared':>8}{'DNAonly':>9}{'RNAonly':>9}{'Jaccard':>9}{'RNA comps':>11}")
for F,lab in (('npip','NPIP'),('tbc','TBC1D3')):
    ns=[l.strip() for l in open(f'{W}/{F}_nodes.txt')]
    rna=[l[1:].strip() for l in open(f'{W}/{F}_rna.fa') if l[0]=='>']
    Ed=E(f'{W}/{F}_dna',set(rna)); Er=E(f'{W}/{F}_rna',set(rna))
    inter=Ed&Er
    cs=comps(rna,Er)
    print(f"{lab:<9}{len(ns):>5}{len(rna):>8}{len(Ed):>7}{len(Er):>7}{len(inter):>8}"
          f"{len(Ed-Er):>9}{len(Er-Ed):>9}{len(inter)/max(len(Ed|Er),1):>9.3f}"
          f"{len(cs):>11}")
    n=len(rna); poss=n*(n-1)//2
    if n<2:
        print(f"{'':>9}   ⚠ only {n} node(s) have read support -- NOT a result, check the BAM path/mount")
        continue
    print(f"{'':>9}   possible {poss}; DNA density {2*len(Ed)/(n*(n-1)):.3f}  RNA density {2*len(Er)/(n*(n-1)):.3f}")
