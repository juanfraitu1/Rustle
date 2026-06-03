#!/usr/bin/env python3
"""Read-intrinsic copy evidence for DAZ: everything from the READ's own BAM tags
(AS / de / NM / minimap2 s1,s2), NOT from comparing the two reference loci.

This is what rustle's VG EM actually consumes. The 7 'diagnostic SNVs' from the
locus-vs-locus comparison are mostly intronic and uncovered by spliced IsoSeq
reads -> the fingerprint term is neutral; the working signals are per-read
alignment identity (NM/de), junction-chain compatibility, and the anchored prior.

Data: /tmp/daz_tags.json (built by /tmp/daz_tags.py from /tmp/daz_aln.sam).
"""
import json, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

rows=json.load(open('/tmp/daz_tags.json'))
EX={'65405150':'R1 resolved','50660416':'R2 parked/tied','24709554':'R3 tied'}
def l1p(x): return math.log10(1+x)

fig=plt.figure(figsize=(13.5,11))
gs=fig.add_gridspec(2,2,height_ratios=[1.35,1],hspace=0.28,wspace=0.22)
BLUE,RED,GREY,GREEN,ORANGE='#2c7fb8','#d73027','#666','#1a9850','#f16913'

# ---------- Panel A: NM(DAZ1) vs NM(DAZ3), colored by minimap2 s2/s1 ----------
axA=fig.add_subplot(gs[0,:])
xs=[l1p(r['NM1']) for r in rows]; ys=[l1p(r['NM3']) for r in rows]
amb=[r['ambig'] for r in rows]
sc=axA.scatter(xs,ys,c=amb,cmap='RdYlBu_r',vmin=0.85,vmax=1.0,s=34,ec='#333',lw=0.3,zorder=3)
mx=max(max(xs),max(ys))*1.05
axA.plot([0,mx],[0,mx],'--',c=GREY,lw=1.2,zorder=1)
axA.text(mx*0.62,mx*0.69,"equal fit to both copies\n= TIED (apportion)",c=GREY,fontsize=10,rotation=37,va='center')
# tied band
axA.fill_between([0,mx],[ -0.07,mx-0.07],[0.07,mx+0.07],color=GREY,alpha=0.07,zorder=0)
cb=fig.colorbar(sc,ax=axA,pad=0.01); cb.set_label("minimap2 ambiguity  s2/s1  (1.0 = 2nd copy fits just as well)",fontsize=9)
for r in rows:
    for k,lab in EX.items():
        if k in r['q']:
            axA.annotate(f"{lab}\nNM {r['NM1']} vs {r['NM3']}",(l1p(r['NM1']),l1p(r['NM3'])),
                         textcoords='offset points',xytext=(9,6),fontsize=8.5,weight='bold',
                         color=(GREEN if 'resolved' in lab else RED),
                         arrowprops=dict(arrowstyle='-',color='#999',lw=0.7))
axA.set_xlabel("edit distance to DAZ1  (NM, log scale)",fontsize=10.5)
axA.set_ylabel("edit distance to DAZ3  (NM)",fontsize=10.5)
axA.set_title("A.  All copy evidence is ON THE READ: per-placement BAM tags, no locus comparison\n"
              "Off-diagonal = the read fits one copy better (assign).  On-diagonal = equal fit (apportion).",
              fontsize=11.5,weight='bold',loc='left')
ticks=[0,1,2,3]; axA.set_xticks([l1p(10**t-1) for t in ticks]); axA.set_xticklabels([f"{10**t-1:.0f}" for t in ticks])
axA.set_yticks([l1p(10**t-1) for t in ticks]); axA.set_yticklabels([f"{10**t-1:.0f}" for t in ticks])
axA.set_xlim(-0.05,mx); axA.set_ylim(-0.05,mx)

# ---------- Panel B: minimap2 s2/s1 separates tied vs resolved by itself ----------
axB=fig.add_subplot(gs[1,0])
tied=[r['ambig'] for r in rows if r['dAS']<=5]; res=[r['ambig'] for r in rows if r['dAS']>5]
axB.hist(res,bins=np.linspace(0.80,1.001,22),color=BLUE,alpha=0.75,label=f"resolved (n={len(res)})")
axB.hist(tied,bins=np.linspace(0.80,1.001,22),color=RED,alpha=0.8,label=f"tied (n={len(tied)})")
axB.axvline(0.98,color=GREY,ls='--',lw=1)
axB.set_xlabel("minimap2 s2/s1 on the primary  (read-intrinsic ambiguity)",fontsize=10)
axB.set_ylabel("reads",fontsize=10)
axB.set_title("B.  minimap2 flags the tied reads itself:\n40/40 tied have s2/s1 >= 0.98",fontsize=10.5,weight='bold',loc='left')
axB.legend(fontsize=9,loc='upper left')

# ---------- Panel C: how rustle's VG EM uses these read tags ----------
axC=fig.add_subplot(gs[1,1]); axC.axis('off'); axC.set_xlim(0,1); axC.set_ylim(0,1)
axC.set_title("C.  What rustle's VG EM consumes (all read-intrinsic)",fontsize=10.5,weight='bold',loc='left')
def box(x,y,w,h,fc,ec):
    axC.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.006",fc=fc,ec=ec,lw=1.3))
def tt(x,y,s,sz=9,c='black',w='normal',st='normal'): axC.text(x,y,s,fontsize=sz,color=c,weight=w,style=st,va='center')
rowsC=[
 (BLUE,'#eaf2f8',"alignment identity",  "from NM / de  -> read_alignment_identity_score",True),
 (GREEN,'#e9f5ee',"junction compatibility","read's splice-junction chain vs each copy",True),
 (ORANGE,'#fdf0e3',"anchored prior",     "unique + decisively-assigned reads' mass",True),
 (GREY,'#f0f0f0',"sequence fingerprint", "the 7 SNVs -> NEUTRAL (intronic, ~0 covered)",False),
]
y=0.82
for c,fc,name,desc,on in rowsC:
    box(0.02,y-0.085,0.96,0.115,fc,c if on else '#bbb')
    tt(0.05,y-0.005,name,9.5,c=(c if on else '#999'),w='bold')
    tt(0.05,y-0.05,desc,8.2,c=('#333' if on else '#999'),st=('normal' if on else 'italic'))
    if not on: tt(0.80,y-0.027,"OFF",9,c='#999',w='bold')
    y-=0.135
box(0.02,0.10,0.96,0.10,'#f4f6fb',BLUE)
tt(0.06,0.15,"=>  per-read weight gamma (sums to 1 across the family's copies)",9,c=BLUE,w='bold')
tt(0.06,0.115,"resolved -> ~1 on the better copy;  tied -> apportioned, never hard-called",8.3,c='#333')

fig.savefig("daz_readtags.png",dpi=145,bbox_inches='tight')
print("wrote daz_readtags.png")
