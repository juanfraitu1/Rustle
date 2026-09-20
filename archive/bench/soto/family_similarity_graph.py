#!/usr/bin/env python3
"""Pairwise similarity graph for a Soto family: nodes = members, edges = exon-sum alignment identity.

Uses the SAME substrate the family edge is computed on (the spliced exon-sum, asm20 + the sensitive tier),
so the graph shows what the method actually sees rather than a separate notion of similarity.
"""
import gzip, re, subprocess, sys, json
from collections import defaultdict
HFA="/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
HGFF="/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
SAM="/home/juanfra/miniforge3/bin/samtools"; MM2="/home/juanfra/miniforge3/bin/minimap2"
OUT="/home/juanfra/winloci_scratch/famgraph"; subprocess.run(["mkdir","-p",OUT],check=False)

members=defaultdict(list)
for ln in open("bench/soto/80_fams.chr.bed"):
    c,s,e,name,*_=ln.rstrip("\n").split("\t")
    members[name.split("|")[1]].append((name.split("|")[0],c,int(s),int(e)))
exons=defaultdict(list)
with gzip.open(HGFF,"rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="exon": continue
        m=re.search(r'(?:^|;)gene=([^;]+)',f[8])
        if m: exons[(f[0],m.group(1).upper())].append((int(f[3])-1,int(f[4])))
def merge(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o

out={}
for fid in sys.argv[1:]:
    seqs={}
    for g,c,s,e in members[fid]:
        ex=merge([(a,b) for a,b in exons.get((c,g.upper()),[]) if a<e and s<b])
        if not ex:   # no annotation -> fall back to the genomic span so the member still appears
            r=subprocess.run([SAM,"faidx",HFA,f"{c}:{s+1}-{e}"],capture_output=True,text=True).stdout
            seqs[g]="".join(r.split("\n")[1:]); continue
        parts=[]
        for a,b in ex:
            r=subprocess.run([SAM,"faidx",HFA,f"{c}:{a+1}-{b}"],capture_output=True,text=True).stdout
            parts.append("".join(r.split("\n")[1:]))
        seqs[g]="".join(parts)
    fa=f"{OUT}/{fid}.fa"
    with open(fa,"w") as fh:
        for k,v in seqs.items(): fh.write(f">{k}\n{v}\n")
    edges={}; selfrep=set()
    for label,args in (("asm20",["-x","asm20"]),("sensitive",["-k","11","-w","5"])):
        paf=subprocess.run([MM2,"-c","--eqx","-X",*args,"-N","50","-p","0.01",fa,fa],
                           capture_output=True,text=True).stdout
        for ln in paf.splitlines():
            f=ln.split("\t")
            if len(f)<12: continue
            a,b=sorted([f[0],f[5]]); al=int(f[10])
            # a==b is a SELF-alignment at a different offset -- an internal repeat inside the gene, not a
            # relationship between two members. NPIP has 13 of them, which is itself informative, but they
            # must not become graph edges.
            if not al: continue
            if a==b: selfrep.add(a); continue
            ident=int(f[9])/al
            shorter=min(int(f[1]),int(f[6]))
            cov=(int(f[3])-int(f[2]))/max(shorter,1)
            k=(a,b)
            if k not in edges or ident>edges[k]["identity"]:
                edges[k]={"identity":round(ident,4),"coverage":round(min(cov,1.0),3),"tier":label}
    out[fid]={"members":{k:len(v) for k,v in seqs.items()},
              "self_repeats":sorted(selfrep),
              "edges":[{"a":a,"b":b,**v} for (a,b),v in sorted(edges.items())]}
    n=len(seqs); npair=n*(n-1)//2
    print(f"{fid}: {n} members, {len(edges)}/{npair} pairs with an alignment, "
          f"{sum(1 for v in edges.values() if v['identity']>=0.80)} at >=0.80 identity")
json.dump(out, open(f"{OUT}/graphs.json","w"), indent=1)
print(f"wrote {OUT}/graphs.json")
