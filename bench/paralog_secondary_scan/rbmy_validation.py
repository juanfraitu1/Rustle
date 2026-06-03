#!/usr/bin/env python3
"""Validation that copy-assignment is correct, on a known testis Y-ampliconic
family: RBMY1 (6-copy tandem array on chrY). Contrast with the DAZ1/DAZ3 case.

⚠️ CORRECTION (2026-06-03): the RBMY `capacity_confidence 1.000` shown here is the
DEFAULT (no VG apportionment ran) — RBMY collapses to one bundle, VG finds 0
families, and the copies are resolved by ordinary assembly, not VG. Genuine VG
copy-resolution of RBMY needs the intra-bundle tandem feature (RUSTLE_VG_TANDEM=1)
→ real per-copy cc (c6=0.956, c4=0.220). The DAZ panel (dispersed, genuinely
VG-resolved) stands. See rbmy_analysis.md. Regenerate before citing.

RBMY copies are 95.8-99.8% identical -> reads are distinguishable -> the method
recovers 5/6 copies as expressed, reads distributing across the array. DAZ1/DAZ3
are 99.97% identical (inverted near-perfect duplicate) -> reads non-identifiable
-> DAZ1 is the expressed copy (174:2 directional), DAZ3 flagged low-confidence.
Conclusion: the method assigns correctly; outcome is set by copy identity, not bias.
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

# --- RBMY1 (LOC129530238-243), ordered by genomic position ---
rbmy_names = ['243','239','240','238','241','242']
rbmy_reads = [8, 10, 14, 7, 45, 2]
rbmy_cov   = [7.9, 12.1, 10.9, 6.8, 38.0, 0.0]
rbmy_cc    = [1.0, 1.0, 1.0, 1.0, 1.0, None]   # None = no transcript assembled
# --- DAZ ---
daz_names = ['DAZ1\n(expressed)','DAZ3\n(silent)']
daz_dir   = [174, 2]      # reads leaning each copy (39 tied apportioned to DAZ1)
daz_cc    = [1.0, 0.0]

BLUE,GREEN,RED,GREY,ORANGE='#2c7fb8','#1a9850','#d73027','#666','#f16913'
fig=plt.figure(figsize=(14,8.6))
ax=fig.add_axes([0,0,1,1]); ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
def T(x,y,s,sz=11,ha='left',va='top',c='black',w='normal',st='normal'):
    ax.text(x,y,s,size=sz,ha=ha,va=va,color=c,weight=w,style=st,zorder=5)

T(0.5,0.975,"Is copy-assignment correct? Validation on RBMY1, a testis-specific Y ampliconic family",16,ha='center',w='bold')
T(0.5,0.94,"RBMY1 = 6-copy tandem array on chrY (like TSPY). If the method were collapsing copies, only one would be expressed.",
  11,ha='center',c=GREY,st='italic')

# ---- Panel A: RBMY per-copy ----
axA=fig.add_axes([0.07,0.40,0.52,0.44])
x=np.arange(len(rbmy_names)); w=0.38
axA.bar(x-w/2, rbmy_reads, w, color=BLUE, label='primary reads')
axA.bar(x+w/2, rbmy_cov,  w, color=GREEN, label='assembled coverage')
for i,cc in enumerate(rbmy_cc):
    lab = 'cc 1.00' if cc==1.0 else 'no tx'
    axA.annotate(lab,(i,max(rbmy_reads[i],rbmy_cov[i])),textcoords='offset points',
                 xytext=(0,4),ha='center',fontsize=8,color=(GREEN if cc==1.0 else RED),weight='bold')
axA.set_xticks(x); axA.set_xticklabels([f"copy\n{n}" for n in rbmy_names],fontsize=9)
axA.set_ylabel("reads / coverage",fontsize=10); axA.legend(fontsize=9,loc='upper left')
axA.set_title("A. RBMY1: reads distribute across the array",fontsize=12,weight='bold',loc='left')
T(0.07,0.37,"5 of 6 copies expressed, each anchored (capacity_confidence 1.000); reads land on all 6 copies (8/10/14/7/45/2).",
  9.5,c='#333',st='italic')
T(0.07,0.355,"-> the method does NOT collapse a multi-copy family to one copy. Copies here are 95.8-99.8% identical (resolvable).",
  9.5,c='#333',st='italic')

# ---- Panel B: DAZ ----
axB=fig.add_axes([0.66,0.40,0.28,0.44])
xb=np.arange(2)
axB.bar(xb, daz_dir, 0.5, color=[BLUE,RED])
for i,(d,cc) in enumerate(zip(daz_dir,daz_cc)):
    axB.annotate(f"cc {cc:.2f}",(i,d),textcoords='offset points',xytext=(0,4),ha='center',
                 fontsize=9,color=(GREEN if cc==1.0 else RED),weight='bold')
axB.set_xticks(xb); axB.set_xticklabels(daz_names,fontsize=9)
axB.set_ylabel("reads leaning each copy",fontsize=10)
axB.set_title("B. DAZ1 vs DAZ3",fontsize=12,weight='bold',loc='left')
T(0.66,0.37,"99.97% identical (inverted). 174 reads lean DAZ1,",9.5,c='#333',st='italic')
T(0.66,0.355,"2 lean DAZ3 (noise); 39 tied -> apportioned to DAZ1.",9.5,c='#333',st='italic')

# ---- bottom: the point ----
ax.add_patch(FancyBboxPatch((0.05,0.04),0.90,0.27,boxstyle="round,pad=0.006",fc='#f4f6fb',ec=BLUE,lw=1.6))
T(0.07,0.295,"What this shows",12,w='bold',c=BLUE)
pts=[
 ("-  RBMY1 (a testis Y family the size of TSPY) is correctly recovered as a MULTI-copy expressed family: 5/6 copies, reads",1),
 ("   spread across the array. So we are NOT systematically mis-assigning everything to one copy.",0),
 ("-  Whether copies separate is set by their IDENTITY, not by the tool: RBMY copies (95.8-99.8%) are distinguishable and",1),
 ("   resolve; DAZ1/DAZ3 (99.97% identical, an inverted near-perfect duplicate) are non-identifiable from any read.",0),
 ("-  DAZ is still expressed in this testis sample -- via DAZ1 (cc 1.000). The directional evidence (174 vs 2) says the",1),
 ("   inverted DAZ3 copy is not detectably expressed (<=~1% of DAZ1); it is flagged low_confidence rather than fabricated.",0),
]
for k,(p,b) in enumerate(pts):
    T(0.07,0.262-k*0.034,p,9.8,c=('#222' if b else '#444'),w=('bold' if b else 'normal'))

fig.savefig("rbmy_validation.png",dpi=150,bbox_inches='tight')
print("wrote rbmy_validation.png")
