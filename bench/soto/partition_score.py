#!/usr/bin/env python3
"""Partition-quality of a --from-genome copies.tsv vs Soto's 83 families (over/under-merge, homogeneity,
completeness), NOT just member recall. Usage: partition_score.py <copies.tsv>"""
import sys
from collections import defaultdict
COP=sys.argv[1]
BED=sys.argv[2] if len(sys.argv)>2 else "bench/soto/80_fams.chr.bed"
copies=[]
for i,ln in enumerate(open(COP)):
    if i==0: continue
    f=ln.rstrip("\n").split("\t")
    if len(f)>=6: copies.append((f[0],f[3],int(f[4]),int(f[5])))
famcnt=defaultdict(int)
for (gw,c,s,e) in copies: famcnt[gw]+=1
members=[]
for ln in open(BED):
    c,s,e,name,*_=ln.rstrip("\n").split("\t")
    members.append((name.split("|")[1],c,int(s),int(e)))
def member_at(c,s,e):
    for (fid,mc,ms,me) in members:
        if mc==c and s<me and e>ms: return fid
    return None
gw2soto=defaultdict(set); soto2gw=defaultdict(set)
recovered=set()
for (gw,c,s,e) in copies:
    if famcnt[gw]<2: continue
    sf=member_at(c,s,e)
    if sf: gw2soto[gw].add(sf); soto2gw[sf].add(gw); recovered.add((sf,c,s))
nmem=len(members)
# member recall
rec=sum(1 for (fid,c,s,e) in members if any(True for (gw,cc,ss,ee) in copies if famcnt[gw]>=2 and cc==c and s<ee and e>ss))
nfam=len([g for g in famcnt if famcnt[g]>=2])
clean=[sf for sf,gws in soto2gw.items() if len(gws)==1 and len(gw2soto[list(gws)[0]])==1]
merges=[gw for gw,sfs in gw2soto.items() if len(sfs)>=2]
splits=[sf for sf,gws in soto2gw.items() if len(gws)>=2]
pure=sum(1 for gw,sfs in gw2soto.items() if len(sfs)==1)
print(f"member_recall={rec}/{nmem}={100*rec/nmem:.1f}%  dna_families={nfam}  "
      f"clean_1:1={len(clean)}  over_merge_fams={len(merges)}(fusing {sum(len(gw2soto[g]) for g in merges)})  "
      f"split_soto_fams={len(splits)}  homogeneity={100*pure/max(len(gw2soto),1):.0f}%  "
      f"completeness={100*len(clean)/max(len(soto2gw),1):.0f}%")
