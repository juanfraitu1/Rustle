#!/usr/bin/env python3
"""Export the Rustle E_r family-definition graph for the Soto NPIP/FAM72 loci.

This mirrors ``homology_edges_all_reps``' nucleotide edge rules: an undirected
edge is admitted by the asm20 tier (identity >= .80) OR the sensitive tier
(-k11 -w5, identity >= .60), provided its alignment covers >= .50 of the
shorter member.  This is the graph on which gamma-quasi-clique partitioning
operates; gamma is .20.
"""
from __future__ import annotations
import csv, itertools, subprocess
from collections import defaultdict, deque
from pathlib import Path

BED = Path("bench/soto/80_fams.chr.bed")
GENOME = Path("../Reference/chm13v2.0.fa")
OUT = Path("bench/soto/family_definition_graphs")
GAMMA, MIN_COV = .20, .50
TIERS = (("asm20", ("-x", "asm20"), .80), ("sensitive", ("-k", "11", "-w", "5"), .60))
COLORS = {"NPIPA":"#4285f4", "NPIPB":"#ea4335", "NPIPP":"#fbbc04", "FAM72":"#34a853"}

def family_members(fid):
    out=[]
    for r in csv.reader(BED.open(), delimiter="\t"):
        if len(r)>=4 and r[3].endswith("|"+fid):
            name,_=r[3].rsplit("|",1); out.append((name,r[0],int(r[1]),int(r[2])))
    return out

def cls(name):
    for k in COLORS:
        if name.startswith(k): return k
    return "other"

def sequences(members):
    seqs=[]
    for name,chrom,start,end in members:
        x=subprocess.check_output(["samtools","faidx",str(GENOME),f"{chrom}:{start+1}-{end}"],text=True)
        seqs.append("".join(z for z in x.splitlines() if not z.startswith(">")))
    return seqs

def edges(seqs, args, min_ident):
    """Pairwise form of nucleotide_edges.

    The production routine batches all copies; this bounded exporter aligns one
    pair at a time to avoid minimap2's all-vs-all repeat-memory blow-up on the
    49-kb NPIPB2 locus.  A pair passes if either endpoint's aligned span covers
    >= MIN_COV of the shorter member, equivalent to the production PAF test.
    """
    ref, qry = OUT/".definition.ref.fa", OUT/".definition.qry.fa"
    ans={}
    for i,j in itertools.combinations(range(len(seqs)),2):
        ref.write_text(f">{i}\n{seqs[i]}\n"); qry.write_text(f">{j}\n{seqs[j]}\n")
        out=subprocess.check_output(["minimap2","-c","--secondary=no",*args,str(ref),str(qry)], text=True, stderr=subprocess.DEVNULL)
        for line in out.splitlines():
            f=line.split("\t")
            if len(f)<12: continue
            ql,tl=float(f[1]),float(f[6]); shorter=min(ql,tl)
            cov=max(float(f[3])-float(f[2]), float(f[8])-float(f[7]))/shorter
            de=next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")),None)
            ident=1-de if de is not None else float(f[9])/max(float(f[10]),1)
            if ident>=min_ident and cov>=MIN_COV:
                old=ans.get((i,j))
                if old is None or ident*cov>old[0]: ans[(i,j)]=(ident*cov,ident,cov)
    ref.unlink(missing_ok=True); qry.unlink(missing_ok=True)
    return ans

def components(n, es):
    a=defaultdict(set)
    for i,j in es: a[i].add(j);a[j].add(i)
    seen=set();out=[]
    for i in range(n):
        if i not in seen:
            q=[i];seen.add(i);c=[]
            while q:
                x=q.pop();c.append(x)
                for y in a[x]:
                    if y not in seen: seen.add(y);q.append(y)
            out.append(sorted(c))
    return out

def export(stem, fid):
    m=family_members(fid); s=sequences(m); OUT.mkdir(parents=True,exist_ok=True)
    bypair={}
    for tier,args,mi in TIERS:
        for pair,val in edges(s,args,mi).items():
            if pair not in bypair or val[0]>bypair[pair][1][0]: bypair[pair]=(tier,val)
    n=len(m); possible=n*(n-1)//2; es=set(bypair); density=len(es)/possible
    comps=components(n,es)
    if len(comps)!=1 or density<GAMMA:
        raise RuntimeError(f"{stem} is not a definition-family certificate: components={len(comps)}, density={density:.3f}")
    with (OUT/(stem+".gfa")).open("w") as f:
        f.write(f"H\tVN:Z:1.0\tTS:Z:rustle_E_r_definition\tGA:f:{GAMMA:.2f}\n")
        for name,chrom,start,end in m:
            f.write(f"S\t{name}\tN\tLN:i:{end-start}\tCL:Z:{cls(name)}\tLO:Z:{chrom}:{start}-{end}\n")
        for (i,j),(tier,(score,ident,cov)) in sorted(bypair.items()):
            f.write(f"L\t{m[i][0]}\t+\t{m[j][0]}\t+\t0M\tET:Z:{tier}\tID:f:{ident:.4f}\tCV:f:{cov:.4f}\tSC:f:{score:.4f}\n")
        for name,*_ in m: f.write(f"P\t{name}\t{name}+\t*\n")
    with (OUT/(stem+".colours.csv")).open("w") as f:
        f.write("Node,Colour\n"); [f.write(f"{x[0]},{COLORS.get(cls(x[0]),'#9aa0a6')}\n") for x in m]
    with (OUT/(stem+".certificate.tsv")).open("w") as f:
        f.write("family\tn_members\tn_edges\tpossible_edges\tdensity\tgamma\tconnected_components\tcertificate\n")
        f.write(f"{stem}\t{n}\t{len(es)}\t{possible}\t{density:.4f}\t{GAMMA:.2f}\t{len(comps)}\tPASS\n")
    with (OUT/(stem+".edges.tsv")).open("w") as f:
        f.write("member_1\tmember_2\ttier\tidentity\tcoverage_of_shorter\tscore\n")
        for (i,j),(tier,(score,ident,cov)) in sorted(bypair.items()): f.write(f"{m[i][0]}\t{m[j][0]}\t{tier}\t{ident:.6f}\t{cov:.6f}\t{score:.6f}\n")
    print(f"{stem}: {len(es)}/{possible} edges, density={density:.3f}, components={len(comps)}; PASS")

if __name__ == '__main__':
    export("ID_154_NPIP.definition", "ID_154")
    export("ID_354_FAM72.definition", "ID_354")
