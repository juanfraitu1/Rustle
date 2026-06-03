#!/usr/bin/env python3
"""A genuinely non-identifiable (tied) DAZ read still carries the SHARED isoform
structure of both copies: it feeds both copies' splice graph identically, so the
transcript STRUCTURE is recovered regardless of assignment -- only the ABUNDANCE
is ambiguous and must be apportioned. Numbers from the adversarial audit
(workflow wf_6c506f4b-cda), each channel independently verified. mathtext-safe.
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

fig=plt.figure(figsize=(14.5,10)); ax=fig.add_axes([0,0,1,1]); ax.axis('off')
ax.set_xlim(0,1); ax.set_ylim(0,1)
BLUE,GREEN,RED,ORANGE,PURP,GREY='#2c7fb8','#1a9850','#d73027','#f16913','#7570b3','#555'
def txt(x,y,s,sz=12,ha='left',va='top',c='black',w='normal',mono=False,st='normal'):
    ax.text(x,y,s,size=sz,ha=ha,va=va,color=c,weight=w,style=st,
            family='monospace' if mono else 'sans-serif',zorder=5)
def box(x,y,w,h,fc,ec='#888',lw=1.2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.004",fc=fc,ec=ec,lw=lw))
def arrow(x1,y1,x2,y2,c='#444',lw=1.6,ls='-'):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=14,
                 lw=lw,color=c,zorder=3,linestyle=ls))

# 7-exon shared pattern (fractions of the model span)
PAT=[(0.00,0.05),(0.12,0.17),(0.26,0.30),(0.40,0.52),(0.61,0.66),(0.74,0.81),(0.91,1.00)]
def gene_model(x0,x1,y,color,fc,exfracs,read=None,h=0.020):
    ax.plot([x0,x1],[y,y],color=color,lw=1.1,zorder=2)            # intron line
    for s,e in exfracs:
        ax.add_patch(Rectangle((x0+s*(x1-x0),y-h/2),(e-s)*(x1-x0),h,fc=fc,ec=color,lw=1.1,zorder=3))
    if read is not None:                                          # read exon blocks
        for s,e in read:
            ax.add_patch(Rectangle((x0+s*(x1-x0),y-h/2),(e-s)*(x1-x0),h,fc=color,ec=color,lw=0,zorder=4))

# ---------- title ----------
txt(0.5,0.982,"A tied read carries the shared structure of both copies: structure is recovered, only abundance is apportioned",
    16.5,ha='center',w='bold')
txt(0.5,0.958,
    "40 of 215 DAZ multimappers are score-tied (|dAS|<=5). 39/40 are byte-identical across both copies on every channel; "
    "0/39 separable by any of 6 verified channels.",
    11.5,ha='center',c=GREY,st='italic')

# ================= PANEL 1: the read feeds both graphs identically =================
box(0.03,0.555,0.94,0.375,'#f7f9fb')
txt(0.05,0.918,"1.  The same tied read maps to -- and supports -- both copies' splice graph (the copies are inverted)",13,w='bold')
X0,X1=0.20,0.90
# DAZ1 model (top)
txt(0.18,0.862,"DAZ1",10.5,ha='right',va='center',w='bold',c=BLUE); txt(0.18,0.846,"(- strand)",8,ha='right',va='center',c=GREY)
gene_model(X0,X1,0.860,BLUE,'#cfe2f0',PAT)
# READ (middle) -- covers exons 2..7 (a 5'-spliced read)
RD=[(0.12,0.17),(0.26,0.30),(0.40,0.52),(0.61,0.66),(0.74,0.81),(0.91,1.00)]
txt(0.18,0.760,"tied READ",10.5,ha='right',va='center',w='bold',c='#222')
txt(0.18,0.744,"16 ex / 15 jx",8,ha='right',va='center',c=GREY)
ax.plot([X0+0.12*(X1-X0),X1],[0.760,0.760],color='#888',lw=1.0,zorder=2)
for s,e in RD:
    ax.add_patch(Rectangle((X0+s*(X1-X0),0.760-0.010),(e-s)*(X1-X0),0.020,fc='#444',ec='#444',zorder=4))
# DAZ3 model (bottom) -- inverted = mirror the exon fractions
PAT3=[(1-e,1-s) for s,e in PAT]
txt(0.18,0.662,"DAZ3",10.5,ha='right',va='center',w='bold',c=GREEN); txt(0.18,0.646,"(+ strand, inverted)",8,ha='right',va='center',c=GREY)
gene_model(X0,X1,0.660,GREEN,'#cde6d7',PAT3)
# "contributes to both" arrows (a few representative exons; down-arrows cross = inversion)
for s,e in [(0.26,0.30),(0.40,0.52),(0.74,0.81)]:
    xm=X0+((s+e)/2)*(X1-X0)
    arrow(xm,0.770,xm,0.850,c=BLUE,lw=1.4,ls=(0,(2,1)))                       # up to DAZ1 (same x)
    xm3=X0+(1-(s+e)/2)*(X1-X0)
    arrow(xm,0.750,xm3,0.670,c=GREEN,lw=1.4,ls=(0,(2,1)))                     # down to DAZ3 (mirrored x)
txt(0.05,0.578,
    "Every exon/junction the read carries exists in BOTH copies (inverted). It adds identical coverage to each graph -- "
    "junction count identical in both for 40/40 reads; 96.6% of junctions and 95.7% of aligned bp on shared edges; 0 genuine copy-specific bp.",
    10,c='#333')

# ================= PANEL 2: structure vs abundance =================
box(0.03,0.305,0.455,0.225,'#eef7f0',ec=GREEN,lw=1.6)
txt(0.05,0.518,"2a.  Isoform STRUCTURE",12.5,w='bold',c=GREEN)
txt(0.05,0.498,"= the transcript model (exon-junction chain)",9.5,c='#333')
# same model drawn small twice
gene_model(0.07,0.45,0.455,BLUE,'#cfe2f0',PAT,h=0.016)
txt(0.07,0.430,"assign -> DAZ1",8.5,c=BLUE,w='bold')
gene_model(0.07,0.45,0.385,GREEN,'#cde6d7',PAT3,h=0.016)
txt(0.07,0.360,"assign -> DAZ3",8.5,c=GREEN,w='bold')
txt(0.06,0.335,"IDENTICAL model either way  ->  recovered deterministically",10,c=GREEN,w='bold')

box(0.515,0.305,0.455,0.225,'#fdf3ec',ec=ORANGE,lw=1.6)
txt(0.535,0.518,"2b.  Isoform ABUNDANCE",12.5,w='bold',c=ORANGE)
txt(0.535,0.498,"= coverage / TPM split between the copies",9.5,c='#333')
# a +1 read going to a fork
txt(0.545,0.470,"+1 read",10,mono=True,w='bold')
arrow(0.60,0.463,0.72,0.445,c=BLUE); txt(0.725,0.448,"DAZ1 ?",10,va='center',c=BLUE,w='bold')
arrow(0.60,0.460,0.72,0.405,c=GREEN); txt(0.725,0.408,"DAZ3 ?",10,va='center',c=GREEN,w='bold')
txt(0.535,0.360,"This is the ONLY ambiguous quantity.",10,c='#333')
txt(0.535,0.335,"Apportion fractionally (EM / expression) -- never hard-assign.",10,c=ORANGE,w='bold')

# ================= banner =================
box(0.03,0.025,0.94,0.255,'#f4f6fb',ec=BLUE,lw=1.8)
txt(0.05,0.267,"What the adversarial audit established (6 channels, each independently re-verified)",13,w='bold',c=BLUE)
pts=[
 "-  Genuinely-tied reads separable by ANY channel: 0 / 39.  Channels: SNV, copy-specific indel, splice junction, transcript strand, position/soft-clip, splice-graph.",
 "   The 40th read is at the |dAS|=5 boundary (not truly tied; DAZ3 already its primary) and its sequence evidence is self-contradictory -- threshold noise.",
 "-  Strand cannot help the inverted pair: each spliced read is canonical in one orientation consistent with BOTH copies; ts=+ on 80/80 placements; 0/40 poly-A.",
 "-  Splice-graph contribution of the tied reads: junctions 258/267 = 96.6% shared (98.9% inversion-invariant); aligned exon bp 87,219/91,098 = 95.7% on shared",
 "   backbone; junction count identical in both copies 40/40; 0/40 cover any of the 7 SNVs; 30/40 spliced, 10/40 single-exon.",
 "-  Therefore: the tied reads ARE the shared isoform structure of both copies. Assemble that structure deterministically; apportion only the abundance.",
]
for k,p in enumerate(pts):
    txt(0.05,0.243-k*0.0245,p,10,c=('#222' if not p.startswith('   ') else '#444'),
        w=('bold' if (not p.startswith('   ') and k in (0,5)) else 'normal'))

fig.savefig("splicegraph_tied_contribution.png",dpi=150,bbox_inches='tight')
print("wrote splicegraph_tied_contribution.png")
