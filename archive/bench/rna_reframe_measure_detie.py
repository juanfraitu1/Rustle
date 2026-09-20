#!/usr/bin/env python3
"""FAITHFUL E_c edge: reproduce read_conflict.rs de_tied exactly.
Read r de-ties loci i,j iff |de_i - de_j| <= delta AND max(de_i,de_j) <= de_max, using r's BEST de per locus
(primary+secondary alignments carry the de tag). Edge iff a locus-pair has >= min_reads de-tie reads.
Family E_c-visible iff >=1 such edge; else E_c drops it (fragments to singletons)."""
import csv, sys, os, pysam
from collections import defaultdict, Counter
BAM="/home/juanfra/winloci_scratch/GGO_mm.bam"
# Output TSV: repo-relative bench/ by default; override with RNA_REFRAME_OUT.
OUT=os.environ.get("RNA_REFRAME_OUT","bench/rna_reframe_detie_rows.tsv")
META="/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
DELTA=0.005; DE_MAX=0.05; MIN_READS=3
coord={}
for r in csv.DictReader(open(META),delimiter='\t'):
    coord[r['id']]=(r['chrom'],int(r['start']),int(r['end']))
bam=pysam.AlignmentFile(BAM,"rb"); refs=set(bam.references)
MAXWIN=400_000; QCAP=8000
def per_read_de(loci):
    """loci: list of (chrom,s,e). Return dict qname -> {locus_idx: best(min) de}."""
    rd=defaultdict(dict)
    for li,(c,s,e) in enumerate(loci):
        if c not in refs: continue
        n=0
        for a in bam.fetch(c,max(0,s),min(e,s+MAXWIN)):
            if a.is_unmapped: continue
            try: de=a.get_tag('de')
            except KeyError: continue
            cur=rd[a.query_name].get(li)
            if cur is None or de<cur: rd[a.query_name][li]=de
            n+=1
            if n>=QCAP: break
    return rd
def detie_edges(rd, nloci):
    pair=Counter()
    for qn,dd in rd.items():
        if len(dd)<2: continue
        items=list(dd.items())
        for x in range(len(items)):
            for y in range(x+1,len(items)):
                li,de_i=items[x]; lj,de_j=items[y]
                if abs(de_i-de_j)<=DELTA and max(de_i,de_j)<=DE_MAX:
                    pair[(min(li,lj),max(li,lj))]+=1
    return {k:v for k,v in pair.items() if v>=MIN_READS}
size_cap=int(sys.argv[1]) if len(sys.argv)>1 else 40
fams=[]
for r in csv.DictReader(open('bench/denovo_families.tsv'),delimiter='\t'):
    fams.append((r['family_id'],int(r['n_copies']),r['members'].split(',')))
rows=[]
for fid,ncop,mem in fams:
    if ncop>size_cap:
        rows.append((fid,ncop,"SKIP_BLOB",0,0,0)); continue
    loci=[coord[dn] for dn in mem if dn in coord]
    rd=per_read_de(loci)
    covered=len(set(li for dd in rd.values() for li in dd))
    edges=detie_edges(rd,len(loci))
    # nodes touched by an edge
    ncov_loci=sum(1 for dn in mem if dn in coord and coord[dn][0] in refs)
    n_detie_reads=sum(edges.values())
    if ncov_loci<2: cls="UNCOVERED"
    elif edges: cls="EC_VISIBLE"
    else: cls="EC_DROPPED"
    rows.append((fid,ncop,cls,ncov_loci,len(edges),n_detie_reads))
cls_fam=Counter(r[2] for r in rows); cls_cop=Counter()
for fid,ncop,cls,cov,ne,nr in rows: cls_cop[cls]+=ncop
print("=== FAITHFUL de-tie E_c edge (delta=%.3f de_max=%.2f min_reads=%d), size_cap=%d ==="%(DELTA,DE_MAX,MIN_READS,size_cap))
print("families:",dict(cls_fam)); print("copies  :",dict(cls_cop))
ev=cls_fam['EC_VISIBLE']+cls_fam['EC_DROPPED']
print("EC_DROPPED/evaluable families: %d/%d = %.1f%%"%(cls_fam['EC_DROPPED'],ev,100*cls_fam['EC_DROPPED']/ev))
evc=cls_cop['EC_VISIBLE']+cls_cop['EC_DROPPED']
print("EC_DROPPED/evaluable copies  : %d/%d = %.1f%%"%(cls_cop['EC_DROPPED'],evc,100*cls_cop['EC_DROPPED']/evc))
# size-stratified drop rate
strat=defaultdict(lambda:[0,0])
for fid,ncop,cls,cov,ne,nr in rows:
    if cls in ("EC_VISIBLE","EC_DROPPED"):
        strat[ncop][0]+= (cls=="EC_DROPPED"); strat[ncop][1]+=1
print("\ndrop-rate by family size:")
for sz in sorted(strat):
    d,t=strat[sz]; print("  size %2d: dropped %d/%d = %.0f%%"%(sz,d,t,100*d/t))
with open(OUT,"w") as fh:
    fh.write("family_id\tn_copies\tclass\tcovered_loci\tn_edges\tn_detie_reads\n")
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
