#!/usr/bin/env python3
"""Evaluate E_c as a blob splitter (companion to bench/ec_blob_split.py)."""
import collections, itertools, pickle, sys
BASE="/mnt/linuxdisk/home/juanfraitu"
S=pickle.load(open(f"{BASE}/ecblob/state.pkl","rb"))
copies,direct,ec,cds,cdsedge = S["copies"],S["direct"],S["ec"],S["cds"],S["cdsedge"]
fam=collections.defaultdict(list)
for n,c in enumerate(copies): fam[c["fam"]].append(n)
fam={k:v for k,v in fam.items() if len(v)>=2}
ALLPOS=len(cdsedge)
rows=[]
for f_,mem in fam.items():
    m=[i for i in mem if i in cds]
    for p in itertools.combinations(sorted(m),2):
        rows.append({"fam":f_,"p":p,"tp":p in cdsedge,"direct":p in direct,"ec":ec.get(p,0)})
n=len(rows)
print(f"arm_f2: {len(fam)} families, same-family pairs with CDS both sides {n}, "
      f"baseline precision {sum(r['tp'] for r in rows)/n:.4f}")
def rep(name, sel, pool=rows, ref=None):
    s=[r for r in pool if sel(r)]
    if not s: print(f"  {name:<32} (empty)"); return
    tp=sum(r["tp"] for r in s); prec=tp/len(s); rec=tp/ALLPOS
    f1=2*prec*rec/(prec+rec) if prec+rec else 0
    print(f"  {name:<32}{len(s):>7}{prec:>11.4f}{rec:>9.4f}{f1:>8.4f}")
print(f"\n{'MARGINAL':<34}{'pairs':>7}{'precision':>11}{'recall':>9}{'F1':>8}")
rep("none", lambda r: True)
rep("direct E_r edge", lambda r: r["direct"])
rep("E_c >= 1 shared tied read", lambda r: r["ec"]>=1)
rep("E_c >= 2", lambda r: r["ec"]>=2)
rep("E_c >= 10", lambda r: r["ec"]>=10)
rep("direct OR E_c >= 1", lambda r: r["direct"] or r["ec"]>=1)
d=[r for r in rows if r["direct"]]
t=[r for r in rows if not r["direct"]]
print(f"\n{'IS E_c ORTHOGONAL TO direct?':<34}{'pairs':>7}{'precision':>11}{'recall':>9}{'F1':>8}")
rep("direct only (reference)", lambda r: True, d)
rep("direct + E_c >= 1", lambda r: r["ec"]>=1, d)
rep("direct + E_c == 0", lambda r: r["ec"]==0, d)
print(f"\n{'WITHIN TRANSITIVE pairs':<34}{'pairs':>7}{'precision':>11}")
for nm,sel in (("all transitive",lambda r: True),
               ("transitive + E_c >= 1",lambda r: r["ec"]>=1),
               ("transitive + E_c == 0",lambda r: r["ec"]==0)):
    s=[r for r in t if sel(r)]
    print(f"  {nm:<32}{len(s):>7}{(sum(r['tp'] for r in s)/len(s) if s else 0):>11.4f}")
# ---- blob split ----
big=max(fam, key=lambda k: len(fam[k]))
mem=fam[big]; ms=set(mem)
print(f"\nBLOB SPLIT — {big}, {len(mem)} members "
      f"({len(mem)*(len(mem)-1)//2} pairs = "
      f"{len(mem)*(len(mem)-1)//2/sum(len(v)*(len(v)-1)//2 for v in fam.values()):.1%} of all)")
m=[i for i in mem if i in cds]
base=[(i,j) for i,j in itertools.combinations(sorted(m),2)]
bt=sum(1 for p in base if p in cdsedge)
print(f"  as one family: pairs {len(base)}  CDS-concordant {bt} = {bt/max(len(base),1):.4f}")
for lab,adjsel in (("E_c >= 1", lambda p: ec.get(p,0)>=1),
                   ("E_c >= 2", lambda p: ec.get(p,0)>=2),
                   ("E_r direct", lambda p: p in direct)):
    adj=collections.defaultdict(set)
    for i,j in itertools.combinations(sorted(mem),2):
        if adjsel((i,j)): adj[i].add(j); adj[j].add(i)
    seen=set(); comps=[]
    for s0 in mem:
        if s0 in seen: continue
        st=[s0]; c=set()
        while st:
            x=st.pop()
            if x in seen: continue
            seen.add(x); c.add(x); st.extend(y for y in adj[x] if y not in seen)
        comps.append(c)
    keep=[c for c in comps if len(c)>=2]
    tp=tot=0
    for c in keep:
        mm=[i for i in c if i in cds]
        for p in itertools.combinations(sorted(mm),2):
            tot+=1; tp += p in cdsedge
    lost=sum(1 for p in base if p in cdsedge)-tp
    print(f"  split by {lab:<12} -> {len(keep):>2} sub-families (sizes "
          f"{sorted((len(c) for c in keep),reverse=True)[:6]}), pairs {tot:>4}, "
          f"concordant {tp:>3} = {(tp/tot if tot else 0):.4f}, true pairs lost {lost}")
