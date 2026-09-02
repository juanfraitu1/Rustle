#!/usr/bin/env python3
"""Are there connectivity levers ORTHOGONAL to the direct-edge rule? (ledger §6bg follow-up)

Candidate levers, all computed on the family's induced E_r/DNA subgraph:
  triangle support - common neighbours of the two endpoints (a chain link has 0)
  min endpoint degree - k-core-style; also an upper bound on lambda
  articulation point - does the pair hang on a single cut vertex (Tarjan)
  family min-degree  - the family-level k-core floor (lambda <= min degree)
`lambda` itself already ships per family (`family_split::edge_connectivity`, emitted as the
`cut_certified` certificate) but is never used as a filter; min-degree bounds it from above.

A lever is only interesting if it adds signal ON TOP OF `direct`, so everything is reported
BOTH marginally and conditioned on direct/transitive.
"""
import collections, itertools, pickle, sys
S=pickle.load(open("/mnt/linuxdisk/home/juanfraitu/fpsrc/state.pkl","rb"))
loci,edges,direct,cds,cdsedge = S["loci"],S["edges"],S["direct"],S["cds"],S["cdsedge"]
import subprocess, os
BASE="/mnt/linuxdisk/home/juanfraitu"
E=sorted(set(edges))
r=subprocess.run([f"{BASE}/rustle_target/release/gamma_refine"],
    input="\n".join(f"{u}\t{v}" for u,v in E)+"\n",capture_output=True,text=True)
fam=collections.defaultdict(list)
for line in r.stdout.split("\n"):
    if line:
        b,m=line.split("\t"); fam[b].append(int(m))
fam={k:v for k,v in fam.items() if len(v)>=2}
ES=set(E)

def artic(nodes, adj):
    """Tarjan articulation points of the induced subgraph."""
    disc={}; low={}; ap=set(); t=[0]
    def dfs(u,parent):
        disc[u]=low[u]=t[0]; t[0]+=1; ch=0
        for v in adj[u]:
            if v not in disc:
                ch+=1; dfs(v,u); low[u]=min(low[u],low[v])
                if parent is not None and low[v]>=disc[u]: ap.add(u)
            elif v!=parent: low[u]=min(low[u],disc[v])
        if parent is None and ch>1: ap.add(u)
    for n in nodes:
        if n not in disc: dfs(n,None)
    return ap

rows=[]
for b,mem in fam.items():
    ms=set(mem)
    adj={i:set() for i in mem}
    for u,v in ES:
        if u in ms and v in ms: adj[u].add(v); adj[v].add(u)
    ap=artic(mem,adj)
    fam_mindeg=min(len(adj[i]) for i in mem)
    m=[i for i in mem if i in cds]
    for i,j in itertools.combinations(sorted(m),2):
        k=(i,j)
        rows.append({"tp":k in cdsedge,"direct":k in ES,
                     "tri":len(adj[i]&adj[j]),
                     "mindeg":min(len(adj[i]),len(adj[j])),
                     "ap":(i in ap) or (j in ap),
                     "fam_mindeg":fam_mindeg,"famsize":len(mem)})
n=len(rows); base=sum(r["tp"] for r in rows)/n
ALLPOS=len(cdsedge)
print(f"pairs {n}   baseline precision {base:.4f}")
def rep(name, sel, rs=rows):
    s=[r for r in rs if sel(r)]
    if not s: print(f"  {name:<34} (empty)"); return
    tp=sum(r["tp"] for r in s)
    prec=tp/len(s); rec=tp/ALLPOS; f1=2*prec*rec/(prec+rec) if prec+rec else 0
    print(f"  {name:<34}{len(s):>7}{prec:>11.4f}{rec:>9.4f}{f1:>8.4f}")
print(f"\n{'MARGINAL':<36}{'pairs':>7}{'precision':>11}{'recall':>9}{'F1':>8}")
rep("none", lambda r: True)
rep("direct edge", lambda r: r["direct"])
rep("triangle support >= 1", lambda r: r["tri"]>=1)
rep("triangle support >= 2", lambda r: r["tri"]>=2)
rep("min endpoint degree >= 2", lambda r: r["mindeg"]>=2)
rep("min endpoint degree >= 3", lambda r: r["mindeg"]>=3)
rep("neither endpoint articulation", lambda r: not r["ap"])
rep("family min-degree >= 2", lambda r: r["fam_mindeg"]>=2)
print(f"\n{'CONDITIONED ON direct (is it orthogonal?)':<36}{'pairs':>7}{'precision':>11}{'recall':>9}{'F1':>8}")
d=[r for r in rows if r["direct"]]
print(f"  {'direct only (reference)':<34}{len(d):>7}"
      f"{sum(r['tp'] for r in d)/len(d):>11.4f}"
      f"{sum(r['tp'] for r in d)/ALLPOS:>9.4f}"
      f"{2*(sum(r['tp'] for r in d)/len(d))*(sum(r['tp'] for r in d)/ALLPOS)/((sum(r['tp'] for r in d)/len(d))+(sum(r['tp'] for r in d)/ALLPOS)):>8.4f}")
rep("direct + triangle >= 1", lambda r: r["direct"] and r["tri"]>=1, d)
rep("direct + triangle >= 2", lambda r: r["direct"] and r["tri"]>=2, d)
rep("direct + not articulation", lambda r: r["direct"] and not r["ap"], d)
rep("direct + min endpoint deg >= 3", lambda r: r["direct"] and r["mindeg"]>=3, d)
print(f"\n{'WITHIN TRANSITIVE pairs (can any lever rescue them?)':<36}{'pairs':>7}{'precision':>11}")
t=[r for r in rows if not r["direct"]]
for name,sel in (("all transitive",lambda r: True),
                 ("transitive + triangle >= 1",lambda r: r["tri"]>=1),
                 ("transitive + triangle >= 2",lambda r: r["tri"]>=2),
                 ("transitive + not articulation",lambda r: not r["ap"])):
    s=[r for r in t if sel(r)]
    print(f"  {name:<34}{len(s):>7}{(sum(r['tp'] for r in s)/len(s) if s else 0):>11.4f}")
