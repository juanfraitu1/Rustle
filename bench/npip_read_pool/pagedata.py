import json, csv, collections
D="/mnt/linuxdisk/tmp/readpool_npip"
S=json.load(open(f"{D}/npip_read_pool.json")); P=json.load(open(f"{D}/posthoc.json")); C=json.load(open(f"{D}/cointoss.json")); F=json.load(open(f"{D}/figdata.json"))
cp=S["copies"]
def gffx(arm):
    out={}
    for ln in open(f"{D}/{arm}.gff3"):
        f=ln.rstrip("\n").split("\t")
        if len(f)<9: continue
        at=dict(x.split("=",1) for x in f[8].split(";") if "=" in x)
        if f[2]=="gene": out[at["Name"]]=dict(c=f[0],s0=int(f[3])-1,e=int(f[4]),st=f[6],ex=[])
        elif f[2]=="exon": out[at["gene"]]["ex"].append((int(f[3])-1,int(f[4])))
    return out
def exov(a,b): return sum(max(0,min(y,q)-max(x,p)) for x,y in a for p,q in b)
node_copies={}
for arm in ("P","GOOD","ALL"):
    L=gffx(arm); by={(v["c"],v["s0"]+1,v["e"]):n for n,v in L.items()}
    cl={by[(r["chrom"],int(r["start"]),int(r["end"]))]:r["cluster_id"] for r in csv.DictReader(open(f"{D}/{arm}.fam.clusters.tsv"),delimiter="\t")}
    on={n:[c["cid"] for c in cp if c["strand"]==v["st"] and exov(v["ex"],c["exons"])>0] for n,v in L.items()}
    fam={cl[n] for n in L if n in cl and on[n]}
    node_copies[arm]=sorted({c for n in L if cl.get(n) in fam for c in on[n]})
tied={r["cid"]:r for r in C}
rows=[]
for c in sorted(cp,key=lambda c:c["s0"]):
    t=tied[c["cid"]]
    rows.append(dict(cid=c["cid"],name=c["name"],s0=c["s0"],e=c["e"],st=c["strand"],prim=t["primaries"],tied=t["tied"],
        loci={a:S["summary"]["arms"][a]["loci_per_copy"].get(c["cid"],0) for a in ("P","GOOD","ALL")},
        node={a:c["cid"] in node_copies[a] for a in ("P","GOOD","ALL")}))
arms=S["summary"]["arms"]
out=dict(rows=rows, arms={a:dict(transcripts={"P":8673,"GOOD":9473,"ALL":14183}[a], loci=arms[a]["loci_chr16"], paf=arms[a]["paf_records"],
        mm2={"P":"39 s","GOOD":"46 s","ALL":"≈ 20 min"}[a], clusters=arms[a]["clusters"], loci_on_copies=arms[a]["loci_on_copies"],
        copies_covered=arms[a]["copies_covered"], span_fused=P["posthoc"][a]["span_fused"], npip=P["npip_clusters"][a],
        copies_with_node=len(node_copies[a])) for a in ("P","GOOD","ALL")},
    prereg=dict(f=S["summary"]["f"], fp=S["summary"]["all_added_fp"], fp_in=S["summary"]["all_added_fp_in_npip_family"], added=P["added"],
        added_clustered=P["added_clustered"], verdict=S["summary"]["verdict"]),
    windows=F["windows"])
print({a:(out["arms"][a]["copies_with_node"], out["arms"][a]["npip"]["copies_with_node"]) for a in out["arms"]})
s=json.dumps(out,separators=(",",":")); assert "</" not in s
open(f"{D}/pagedata.json","w").write(s); print(len(s),"bytes")
