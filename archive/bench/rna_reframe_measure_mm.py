#!/usr/bin/env python3
"""Exact E_c edge measurement: a read multimapping across >=2 loci of an E_r family = a de-tie (E_c edge).
Family is E_c-visible iff >=1 such intra-family multimapper exists; else E_c drops it (fragments)."""
import csv, sys, os, pysam
from collections import defaultdict, Counter
BAM="/home/juanfra/winloci_scratch/GGO_mm.bam"
# Output TSV: repo-relative bench/ by default; override with RNA_REFRAME_OUT.
OUT=os.environ.get("RNA_REFRAME_OUT","bench/rna_reframe_mm_rows.tsv")
META="/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
coord={}
for r in csv.DictReader(open(META),delimiter='\t'):
    coord[r['id']]=(r['chrom'],int(r['start']),int(r['end']))
bam=pysam.AlignmentFile(BAM,"rb"); refs=set(bam.references)
MAXWIN=400_000; QCAP=6000
def qnames(chrom,s,e):
    if chrom not in refs: return set()
    out=set(); n=0
    for a in bam.fetch(chrom,max(0,s),min(e,s+MAXWIN)):
        if a.is_unmapped: continue
        out.add(a.query_name); n+=1
        if n>=QCAP: break
    return out
size_cap=int(sys.argv[1]) if len(sys.argv)>1 else 40
fams=[]
for r in csv.DictReader(open('bench/denovo_families.tsv'),delimiter='\t'):
    fams.append((r['family_id'],int(r['n_copies']),r['members'].split(',')))
rows=[]
for fid,ncop,mem in fams:
    if ncop>size_cap:
        rows.append((fid,ncop,"SKIP_BLOB",0,0)); continue
    locus_q=[]
    for dn in mem:
        if dn in coord:
            c,s,e=coord[dn]; locus_q.append(qnames(c,s,e))
    covered=sum(1 for q in locus_q if q)
    # count reads appearing at >=2 loci
    cnt=Counter()
    for q in locus_q:
        for nm in q: cnt[nm]+=1
    mm=sum(1 for nm,c2 in cnt.items() if c2>=2)
    if covered<2: cls="UNCOVERED"
    elif mm>=1: cls="EC_VISIBLE"
    else: cls="EC_DROPPED"
    rows.append((fid,ncop,cls,covered,mm))
cls_fam=Counter(r[2] for r in rows); cls_cop=Counter()
for fid,ncop,cls,cov,mm in rows: cls_cop[cls]+=ncop
print("=== EXACT intra-family multimapper (E_c de-tie) classification, size_cap=%d ==="%size_cap)
print("families:",dict(cls_fam)); print("copies  :",dict(cls_cop))
ev=cls_fam['EC_VISIBLE']+cls_fam['EC_DROPPED']
print("EC_DROPPED/evaluable families: %d/%d = %.1f%%"%(cls_fam['EC_DROPPED'],ev,100*cls_fam['EC_DROPPED']/ev))
evc=cls_cop['EC_VISIBLE']+cls_cop['EC_DROPPED']
print("EC_DROPPED/evaluable copies  : %d/%d = %.1f%%"%(cls_cop['EC_DROPPED'],evc,100*cls_cop['EC_DROPPED']/evc))
# show the EC_VISIBLE families
print("\nEC_VISIBLE families (fid,ncop,covered,n_multimapper_reads):")
for r in rows:
    if r[2]=="EC_VISIBLE": print("  ",r[0],r[1],r[3],r[4])
with open(OUT,"w") as fh:
    fh.write("family_id\tn_copies\tclass\tcovered\tn_mm\n")
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
