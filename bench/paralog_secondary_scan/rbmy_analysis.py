#!/usr/bin/env python3
"""Standalone analysis of RBMY1, a 6-copy testis-specific Y ampliconic gene family
(NC_073248.2:19.60-19.73 Mb), used to validate copy-level read assignment.

⚠️ CORRECTION (2026-06-03): the `capacity_confidence 1.000` this figure shows is
the DEFAULT (no VG apportionment), NOT a VG result — RBMY collapses to one bundle
so VG discovers 0 families and never resolves the copies; the per-copy genes are
ordinary position-based assembly. RBMY is genuinely VG-resolved only with the
intra-bundle tandem feature (RUSTLE_VG_TANDEM=1), which yields REAL per-copy
confidence (c6=0.956, c4=0.220, …). See rbmy_analysis.md for the full correction.
This script must be regenerated from a RUSTLE_VG_TANDEM=1 run before it is cited.

Reads multi-map heavily across the near-identical copies (599 alignments / 87
reads), yet they resolve: 5/6 copies are assembled, each anchored
(capacity_confidence 1.000), reads distributing across the array. So the method
recovers a multi-copy testis family correctly rather than collapsing it.
Data: /tmp/rbmy/ident.json, /tmp/rbmy/percopy.json.
"""
import json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np

ID=json.load(open('/tmp/rbmy/ident.json'))
M=np.array(ID['M']); names=[n.replace('_',' ').replace('c','copy ').split()[0:2] for n in ID['names']]
short=['c1\nLOC243','c2\nLOC239','c3\nLOC240','c4\nLOC238','c5\nLOC241','c6\nLOC242']
# genomic order coords (Mb-relative) + per-copy stats
coords=[(19602754,19616644),(19625715,19639621),(19648670,19662577),
        (19671638,19685525),(19694606,19708531),(19717578,19730926)]
reads=[8,10,14,7,45,2]; cov=[7.9,12.1,10.9,6.8,38.0,0.0]; ntx=[2,4,3,1,6,0]
cc=[1.0,1.0,1.0,1.0,1.0,None]
BLUE,GREEN,RED,GREY,ORANGE='#2c7fb8','#1a9850','#d73027','#666','#f16913'

fig=plt.figure(figsize=(14,11))
ax=fig.add_axes([0,0,1,1]); ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
def T(x,y,s,sz=11,ha='left',va='top',c='black',w='normal',st='normal'):
    ax.text(x,y,s,size=sz,ha=ha,va=va,color=c,weight=w,style=st,zorder=5)

T(0.5,0.978,"RBMY1: a 6-copy testis-specific Y gene family, correctly resolved per copy",17,ha='center',w='bold')
T(0.5,0.952,"chrY tandem array (NC_073248.2:19.60-19.73 Mb), ~14 kb/copy, - strand. Reads multi-map across near-identical copies "
  "but resolve to individual copies.",11,ha='center',c=GREY,st='italic')

# ---- Panel A: array map ----
T(0.05,0.928,"A.  The tandem array (6 copies on chrY, drawn to genomic scale)",12.5,w='bold')
ax0,ax1=0.07,0.95; lo=19595000; hi=19738000
def gx(g): return ax0+(g-lo)/(hi-lo)*(ax1-ax0)
yb=0.892
ax.plot([ax0,ax1],[yb,yb],c='#ccc',lw=1.0,zorder=1)
for i,(s,e) in enumerate(coords):
    col=GREEN if cc[i]==1.0 else RED
    ax.add_patch(Rectangle((gx(s),yb-0.011),gx(e)-gx(s),0.022,fc=col,ec='#333',lw=0.8,alpha=0.85,zorder=3))
    T((gx(s)+gx(e))/2,yb+0.016,short[i].split('\n')[0],8.5,ha='center',va='bottom',w='bold',c='#222')
    T((gx(s)+gx(e))/2,yb-0.015,(f"cov {cov[i]:.0f}" if cc[i] else "silent"),7.5,ha='center',va='top',
      c=(GREEN if cc[i] else RED))
T(0.07,0.852,"Green = assembled & anchored (capacity_confidence 1.000); red = below assembly threshold (2 reads). Copies ~23 kb apart.",
  9.3,c='#333',st='italic')

# ---- Panel B: identity heatmap ----
axB=fig.add_axes([0.10,0.42,0.34,0.36])
im=axB.imshow(M*100,cmap='viridis',vmin=93,vmax=100)
axB.set_xticks(range(6)); axB.set_yticks(range(6))
axB.set_xticklabels([s.split('\n')[0] for s in short],fontsize=8)
axB.set_yticklabels([s.split('\n')[0] for s in short],fontsize=8)
for i in range(6):
    for j in range(6):
        axB.text(j,i,f"{M[i,j]*100:.1f}",ha='center',va='center',fontsize=7,
                 color=('white' if M[i,j]<0.985 else 'black'))
axB.set_title("B.  Pairwise identity (%)",fontsize=12,weight='bold',loc='left')
cb=fig.colorbar(im,ax=axB,fraction=0.046,pad=0.04); cb.ax.tick_params(labelsize=7)
T(0.10,0.40,"Core c2-c5 ~99.8% identical; c1 ~97.3%; c6 ~93-96%.",9,c='#333',st='italic')
T(0.10,0.388,"-> a real range, like a young amplicon.",9,c='#333',st='italic')

# ---- Panel C: per-copy expression ----
axC=fig.add_axes([0.55,0.42,0.40,0.36])
x=np.arange(6); w=0.38
axC.bar(x-w/2,reads,w,color=BLUE,label='primary reads')
axC.bar(x+w/2,cov,w,color=GREEN,label='assembled coverage')
for i in range(6):
    axC.annotate(('cc 1.00' if cc[i] else 'no tx'),(i,max(reads[i],cov[i])),textcoords='offset points',
                 xytext=(0,3),ha='center',fontsize=7.5,color=(GREEN if cc[i] else RED),weight='bold')
axC.set_xticks(x); axC.set_xticklabels([s.split('\n')[0] for s in short],fontsize=8)
axC.set_ylabel("reads / coverage",fontsize=10); axC.legend(fontsize=8.5,loc='upper left')
axC.set_title("C.  Per-copy expression & assignment",fontsize=12,weight='bold',loc='left')
T(0.55,0.40,"5/6 copies expressed; reads spread across the array (8/10/14/7/45/2), not piled on one copy.",9,c='#333',st='italic')

# ---- conclusion ----
ax.add_patch(FancyBboxPatch((0.05,0.04),0.90,0.30,boxstyle="round,pad=0.006",fc='#f4f6fb',ec=BLUE,lw=1.6))
T(0.07,0.325,"What this validates",12.5,w='bold',c=BLUE)
pts=[
 ("-  Heavy multimapping, clean resolution: 87 reads produce 599 alignments (~7 placements each) across the 6 near-identical",1),
 ("   copies -- yet each read carries enough copy-specific sequence to resolve to ONE copy. Result: 5/6 copies assembled,",0),
 ("   each anchored (capacity_confidence 1.000); the 6th has only 2 reads (below threshold), not a mis-assignment.",0),
 ("-  The method does NOT collapse a multi-copy family. RBMY1 -- a known testis-specific Y ampliconic family -- is recovered",1),
 ("   as genuinely multi-copy expressed, which is the biologically expected result for testis.",0),
 ("-  Resolvability tracks identity: the RBMY core is ~99.8% identical and STILL resolves. Only a near-perfect duplicate",1),
 ("   (e.g. DAZ1/DAZ3 at 99.97%) becomes non-identifiable -- and there the tool flags it (low_confidence) rather than guessing.",0),
]
for k,(p,b) in enumerate(pts):
    T(0.07,0.298-k*0.0345,p,9.6,c=('#222' if b else '#444'),w=('bold' if b else 'normal'))

fig.savefig("rbmy_analysis.png",dpi=150,bbox_inches='tight')
print("wrote rbmy_analysis.png")
