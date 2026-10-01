# post hoc (after the prereg's numbers were seen): span-level fusion, absorbed copies, where the off-copy NPIP-family members sit,
# and whether ALL's extra loci survive into ANY multi-member cluster.
import json, collections, csv
D="/mnt/linuxdisk/tmp/readpool_npip"
d=json.load(open(f"{D}/npip_read_pool.json"))
cp=d["copies"]
def gff(arm):
    out={}
    for ln in open(f"{D}/{arm}.gff3"):
        f=ln.rstrip("\n").split("\t")
        if len(f)<9: continue
        at=dict(x.split("=",1) for x in f[8].split(";") if "=" in x)
        if f[2]=="gene": out[at["Name"]]=dict(c=f[0],s0=int(f[3])-1,e=int(f[4]),st=f[6])
    return out
def clusters(arm, L):
    by={(v["c"],v["s0"]+1,v["e"]):n for n,v in L.items()}; out={}
    for r in csv.DictReader(open(f"{D}/{arm}.fam.clusters.tsv"),delimiter="\t"): out[by[(r["chrom"],int(r["start"]),int(r["end"]))]]=r["cluster_id"]
    return out
def exov(s0,e,ex): return sum(max(0,min(e,y)-max(s0,x)) for x,y in ex)
res={}
for arm in ("P","GOOD","ALL"):
    L=gff(arm); cl=clusters(arm,L)
    span_on={n:[c["name"] for c in cp if c["strand"]==v["st"] and exov(v["s0"],v["e"],c["exons"])>0] for n,v in L.items()}
    rows=d["loci"][arm]
    rep_on={n:r["copies"] for n,r in rows.items()}
    span_fused=[n for n in L if len(span_on[n])>=2]
    absorbed=[c["name"] for c in cp if not any(c["cid"] in rep_on.get(n,[]) for n in L) and any(c["name"] in span_on[n] for n in L)]
    fam=[n for n,r in rows.items() if r["in_fam"]]
    off=[n for n in fam if rows[n]["klass"]!="npip_copy"]
    # where the off-copy family members sit: inside a copy's span on the other strand / inside an NPIP copy's intron or flank / elsewhere
    def where(n):
        r=rows[n]; mid=(r["s0"]+r["e"])//2
        for c in cp:
            if c["s0"]<r["e"] and r["s0"]<c["e"]:
                return "overlaps a copy, opposite strand" if c["strand"]!=r["strand"] else "inside a copy's span, off its exons"
        dist=min(min(abs(c["s0"]-r["e"]),abs(r["s0"]-c["e"])) for c in cp)
        return "within 100 kb of a copy" if dist<=100_000 else "farther than 100 kb"
    res[arm]=dict(span_fused=len(span_fused), span_fused_ex=[(n,span_on[n]) for n in span_fused][:6], absorbed=absorbed,
                  fam=len(fam), off=len(off), off_where=dict(collections.Counter(where(n) for n in off)),
                  multi_clustered=sum(1 for n in L if n in cl), loci=len(L))
    print(arm, json.dumps(res[arm])[:900])
# ALL-added loci: how many end up in any multi-member cluster
G=gff("GOOD"); A=gff("ALL"); clA=clusters("ALL",A)
gs=collections.defaultdict(list)
for v in G.values(): gs[(v["c"],v["st"])].append(v)
print("note: ALL-added set = summary all_added =", d["summary"]["all_added"])
def gffx(arm):
    out={}
    for ln in open(f"{D}/{arm}.gff3"):
        f=ln.rstrip("\n").split("\t")
        if len(f)<9: continue
        at=dict(x.split("=",1) for x in f[8].split(";") if "=" in x)
        if f[2]=="gene": out[at["Name"]]=dict(c=f[0],s0=int(f[3])-1,e=int(f[4]),st=f[6],ex=[])
        elif f[2]=="exon": out[at["gene"]]["ex"].append((int(f[3])-1,int(f[4])))
    return out
GX=gffx("GOOD"); AX=gffx("ALL")
gidx=collections.defaultdict(list)
for v in GX.values(): gidx[v["st"]].append(v)
def hit(L):
    return any(g["s0"]<L["e"] and L["s0"]<g["e"] and sum(max(0,min(b,y)-max(a,x)) for a,b in L["ex"] for x,y in g["ex"])>0 for g in gidx[L["st"]])
added=[n for n,L in AX.items() if not hit(L)]
rowsA=d["loci"]["ALL"]
fam_id=collections.Counter(clA[n] for n,r in rowsA.items() if r["in_fam"] and n in clA).most_common(1)[0][0]
in_any=sum(1 for n in added if n in clA); in_npip=sum(1 for n in added if clA.get(n)==fam_id)
print(f"ALL-added loci {len(added)}: in any multi-member cluster {in_any} ({in_any/len(added):.1%}), singletons {len(added)-in_any}, in the NPIP family {in_npip}")
json.dump(dict(posthoc=res, added=len(added), added_clustered=in_any, added_in_npip=in_npip), open(f"{D}/posthoc.json","w"))
print("--- NPIP family as NODES (clusters.tsv members) vs LOCI (incl. records folded into a member)")
nodes_out={}
for arm in ("P","GOOD","ALL"):
    L=gff(arm); cl=clusters(arm,L); rows=d["loci"][arm]
    fid=collections.Counter(cl[n] for n,r in rows.items() if r["in_fam"] and n in cl).most_common(1)[0][0]
    nodes=[n for n in L if cl.get(n)==fid]
    folded=[n for n,r in rows.items() if r["in_fam"] and n not in cl]
    k=collections.Counter(rows[n]["klass"] for n in nodes)
    cov=set(c for n in nodes for c in rows[n]["copies"])
    nodes_out[arm]=dict(nodes=len(nodes), by_class=dict(k), folded=len(folded), copies_with_node=len(cov),
                        folded_by_class=dict(collections.Counter(rows[n]["klass"] for n in folded)))
    print(arm, nodes_out[arm])
r=json.load(open(f"{D}/posthoc.json")); r["family_nodes"]=nodes_out; json.dump(r,open(f"{D}/posthoc.json","w"))
print("--- every cluster that holds >= 1 locus on an NPIP copy (nodes; folded loci counted separately)")
allc={}
for arm in ("P","GOOD","ALL"):
    L=gffx(arm); cl=clusters(arm,L)
    on={n:[c["cid"] for c in cp if c["strand"]==v["st"] and sum(max(0,min(b,y)-max(a,x)) for a,b in v["ex"] for x,y in c["exons"])>0] for n,v in L.items()}
    ncl={cl[n] for n in L if n in cl and on[n]}
    nodes=[n for n in L if cl.get(n) in ncl]
    def klass(n):
        v=L[n]
        if on[n]: return "copy"
        if any(c["strand"]!=v["st"] and c["s0"]<v["e"] and v["s0"]<c["e"] for c in cp): return "antisense"
        if any(c["strand"]==v["st"] and c["s0"]<v["e"] and v["s0"]<c["e"] for c in cp): return "offexon"
        return "elsewhere"
    k=collections.Counter(klass(n) for n in nodes)
    covered=set(c for n in nodes for c in on[n])
    allc[arm]=dict(clusters=len(ncl), nodes=len(nodes), by_class=dict(k), copies_with_node=len(covered))
    print(arm, allc[arm])
r=json.load(open(f"{D}/posthoc.json")); r["npip_clusters"]=allc; json.dump(r,open(f"{D}/posthoc.json","w"))
