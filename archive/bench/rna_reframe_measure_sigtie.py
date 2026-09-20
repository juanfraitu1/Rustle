#!/usr/bin/env python3
"""Significance-tie variant (RUSTLE_CONFLICT_SIG=1): sig_tied uses m_x=de_x*aln_len_x; tied iff
eps^|m_i-m_j| >= alpha  (eps=0.001, alpha=1e-3  => |m_i-m_j| <= 1 mismatch) AND max de <= de_max."""
import csv, sys, pysam
from collections import defaultdict, Counter
BAM="/home/juanfra/winloci_scratch/GGO_mm.bam"
META="/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
EPS=0.001; ALPHA=1e-3; DE_MAX=0.05; MIN_READS=3
import math
coord={}
for r in csv.DictReader(open(META),delimiter='\t'):
    coord[r['id']]=(r['chrom'],int(r['start']),int(r['end']))
bam=pysam.AlignmentFile(BAM,"rb"); refs=set(bam.references)
MAXWIN=400_000; QCAP=8000
def per_read(loci):
    rd=defaultdict(dict)  # qn -> {li:(de,alnlen)} keep min-de
    for li,(c,s,e) in enumerate(loci):
        if c not in refs: continue
        n=0
        for a in bam.fetch(c,max(0,s),min(e,s+MAXWIN)):
            if a.is_unmapped: continue
            try: de=a.get_tag('de')
            except KeyError: continue
            ln=a.reference_length or a.query_alignment_length or 0
            cur=rd[a.query_name].get(li)
            if cur is None or de<cur[0]: rd[a.query_name][li]=(de,ln)
            n+=1
            if n>=QCAP: break
    return rd
def sig_edges(rd):
    pair=Counter()
    for qn,dd in rd.items():
        if len(dd)<2: continue
        items=list(dd.items())
        for x in range(len(items)):
            for y in range(x+1,len(items)):
                li,(di,ln_i)=items[x]; lj,(dj,ln_j)=items[y]
                if max(di,dj)>DE_MAX: continue
                delta_cols=abs(di*ln_i - dj*ln_j)
                if EPS**delta_cols >= ALPHA:
                    pair[(min(li,lj),max(li,lj))]+=1
    return {k:v for k,v in pair.items() if v>=MIN_READS}
size_cap=int(sys.argv[1]) if len(sys.argv)>1 else 40
fams=[]
for r in csv.DictReader(open('bench/denovo_families.tsv'),delimiter='\t'):
    fams.append((r['family_id'],int(r['n_copies']),r['members'].split(',')))
rows=[]
for fid,ncop,mem in fams:
    if ncop>size_cap:
        rows.append((fid,ncop,"SKIP_BLOB")); continue
    loci=[coord[dn] for dn in mem if dn in coord]
    ncov=sum(1 for c,s,e in loci if c in refs)
    rd=per_read(loci); edges=sig_edges(rd)
    cls="UNCOVERED" if ncov<2 else ("EC_VISIBLE" if edges else "EC_DROPPED")
    rows.append((fid,ncop,cls))
cls_fam=Counter(r[2] for r in rows); cls_cop=Counter()
for fid,ncop,cls in rows: cls_cop[cls]+=ncop
print("=== SIG-tie E_c edge (eps=%.3f alpha=%.0e => |dm|<=1), size_cap=%d ==="%(EPS,ALPHA,size_cap))
print("families:",dict(cls_fam)); print("copies  :",dict(cls_cop))
ev=cls_fam['EC_VISIBLE']+cls_fam['EC_DROPPED']
print("EC_DROPPED/evaluable families: %d/%d = %.1f%%"%(cls_fam['EC_DROPPED'],ev,100*cls_fam['EC_DROPPED']/ev))
evc=cls_cop['EC_VISIBLE']+cls_cop['EC_DROPPED']
print("EC_DROPPED/evaluable copies  : %d/%d = %.1f%%"%(cls_cop['EC_DROPPED'],evc,100*cls_cop['EC_DROPPED']/evc))
