#!/usr/bin/env python3
"""Corrected DAZ summary (v3): real gene models (intron/exon) + an honest
CONFIDENCE SPECTRUM (resolvable is a gradient, not yes/no).

Read-intrinsic: per-read AS/de/NM + minimap2 s2/s1. 6/7 SNVs intronic.
Data: /tmp/daz_tags.json, /tmp/igv_daz/daz1_ex.txt, /tmp/igv_daz/daz3_ex.txt
"""
import json, math, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

rows=json.load(open('/tmp/daz_tags.json'))
DAZ1=(42783133,42859657); DAZ3=(42879918,42945552); BLK=(42828911,42839709)
SNVS=[42812422,42815771,42823876,42840676,42842573,42849269,42855142]; EXONIC=42840676
def load_ex(p): return [tuple(map(int,l.split())) for l in open(p) if l.strip()]
DAZ1_EX=load_ex('/tmp/igv_daz/daz1_ex.txt'); DAZ3_EX=load_ex('/tmp/igv_daz/daz3_ex.txt')
EX={'65405150':('R1 resolved','#1a9850'),'50660416':('R2 parked','#d73027'),'24709554':('R3 tied','#d73027')}
BLUE,GREEN,RED,ORANGE,GREY,GOLD='#2c7fb8','#1a9850','#d73027','#f16913','#666','#e8a200'

# confidence tiers by the AS-ratio of the two placements (always present on every
# record; matches minimap2 s2/s1 to ~0.004 but does not depend on the primary-only s2 tag)
def tier(s): return 0 if s>=0.999 else 1 if s>=0.95 else 2
TC=collections.Counter(tier(r['ambig']) for r in rows)
N=len(rows); n_tied,n_weak,n_conf=TC[0],TC[1],TC[2]

fig=plt.figure(figsize=(15.5,17.6))
ax=fig.add_axes([0,0,1,1]); ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
def T(x,y,s,sz=11,ha='left',va='top',c='black',w='normal',mono=False,st='normal'):
    ax.text(x,y,s,size=sz,ha=ha,va=va,color=c,weight=w,style=st,
            family='monospace' if mono else 'sans-serif',zorder=5)
def panel(x,y,w,h,fc,ec='#bbb',lw=1.2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.004",fc=fc,ec=ec,lw=lw))

T(0.5,0.990,"Assigning DAZ multimapping reads to copies — the decision is on the READS (real BAM)",18,ha='center',w='bold')
T(0.5,0.973,f"{N} reads place on both DAZ1 and DAZ3 (99.97% identical). Each read is assigned by its OWN alignment to each copy",12,ha='center',c=GREY,st='italic')
T(0.5,0.962,"(AS / de / NM / ms — tags present on every placement) — and 'resolvable' is a confidence gradient, not yes/no.",12,ha='center',c=GREY,st='italic')

# ===================== Panel 1: real gene models =====================
panel(0.03,0.770,0.94,0.180,'#f6f8fb')
T(0.05,0.943,"1.  The two genes (real NCBI models): the few SNVs sit mostly in INTRONS, between exons",13,w='bold')
X0,X1=0.145,0.865
def gxc(g,lo,hi): return X0+(g-lo)/(hi-lo)*(X1-X0)
def gene_track(ymid,exons,lo,hi,col,fillc,name,strand):
    ax.plot([X0,X1],[ymid,ymid],c=col,lw=1.0,zorder=2)
    h=0.015
    for s,e in exons:
        xs=gxc(s,lo,hi); xe=gxc(e,lo,hi)
        ax.add_patch(Rectangle((xs,ymid-h/2),max(xe-xs,0.0009),h,fc=fillc,ec=col,lw=0.6,zorder=3))
    T(X0-0.008,ymid,name,10.5,ha='right',va='center',w='bold',c=col)
    T(X1+0.008,ymid,strand,9.5,ha='left',va='center',c=GREY)
bcx=(gxc(BLK[0],*DAZ1)+gxc(BLK[1],*DAZ1))/2
# DAZ1 (top)
y1=0.898
ax.add_patch(Rectangle((gxc(BLK[0],*DAZ1),y1-0.013),gxc(BLK[1],*DAZ1)-gxc(BLK[0],*DAZ1),0.026,
             fc=ORANGE,ec='none',alpha=0.20,zorder=1))
T(bcx,0.921,"10.8 kb DAZ1-only block",8.5,ha='center',va='bottom',c=ORANGE,w='bold')
gene_track(y1,DAZ1_EX,*DAZ1,BLUE,'#bcd8ee',"DAZ1","(- str)")
for d in SNVS:
    exonic=(d==EXONIC); dx=gxc(d,*DAZ1)
    if exonic:
        ax.plot([dx,dx],[y1-0.009,y1+0.016],c=RED,lw=2.2,zorder=6)        # spans the track
        ax.plot(dx,y1+0.016,marker='v',ms=7,c=RED,zorder=6)
        T(dx,0.879,"exonic SNV",8,ha='center',va='top',c=RED,w='bold')   # label BELOW the track
    else:
        ax.plot([dx,dx],[y1+0.008,y1+0.016],c='#e26',lw=1.1,zorder=6)
        ax.plot(dx,y1+0.016,marker='v',ms=3.5,c='#e26',zorder=6)
# DAZ3 (bottom)
y3=0.851
gene_track(y3,DAZ3_EX,*DAZ3,GREEN,'#c4e3d2',"DAZ3","(+ str)")
T(bcx,y3-0.011,"(no block here)",7.5,ha='center',va='top',c=GREY,st='italic')
T(0.05,0.812,"Red ticks = the 7 copy-distinguishing SNVs. 6 sit in introns (in the gaps between exon boxes) so a spliced read never reads",10,c='#333',st='italic')
T(0.05,0.800,f"them; only 1 is exonic. DAZ1 = {len(DAZ1_EX)} exons, DAZ3 = {len(DAZ3_EX)} exons, and the exons are nearly identical -> almost no per-base signal on a read.",10,c='#333',st='italic')

# ===================== Panel 2: read-intrinsic scatter =====================
panel(0.03,0.448,0.94,0.318,'#ffffff')
T(0.05,0.758,"2.  What separates a read = its OWN alignment to each copy (no locus comparison)",13,w='bold')
axS=fig.add_axes([0.085,0.490,0.40,0.235])
def l1p(x): return math.log10(1+x)
xs=[l1p(r['NM1']) for r in rows]; ys=[l1p(r['NM3']) for r in rows]; amb=[r['ambig'] for r in rows]
sc=axS.scatter(xs,ys,c=amb,cmap='RdYlBu_r',vmin=0.85,vmax=1.0,s=30,ec='#333',lw=0.3,zorder=3)
mx=max(max(xs),max(ys))*1.05
axS.plot([0,mx],[0,mx],'--',c=GREY,lw=1.1,zorder=1)
axS.fill_between([0,mx],[-0.09,mx-0.09],[0.09,mx+0.09],color=GREY,alpha=0.08,zorder=0)
axS.text(mx*0.46,mx*0.6,"equal fit = TIED",c=GREY,fontsize=10,rotation=40,va='center')
for r in rows:
    for k,(lab,col) in EX.items():
        if k in r['q']:
            axS.annotate(lab,(l1p(r['NM1']),l1p(r['NM3'])),textcoords='offset points',
                         xytext=(8,6),fontsize=9,weight='bold',color=col,
                         arrowprops=dict(arrowstyle='-',color='#aaa',lw=0.7))
axS.set_xlabel("edit distance to DAZ1  (NM)",fontsize=11); axS.set_ylabel("edit distance to DAZ3  (NM)",fontsize=11)
tk=[0,1,2,3]; axS.set_xticks([l1p(10**t-1) for t in tk]); axS.set_xticklabels([f"{10**t-1:.0f}" for t in tk],fontsize=9)
axS.set_yticks([l1p(10**t-1) for t in tk]); axS.set_yticklabels([f"{10**t-1:.0f}" for t in tk],fontsize=9)
axS.set_xlim(-0.05,mx); axS.set_ylim(-0.05,mx)
cax=fig.add_axes([0.497,0.490,0.011,0.235]); cb=fig.colorbar(sc,cax=cax)
cb.set_label("AS ratio  (DAZ1 vs DAZ3)",fontsize=9.5); cb.ax.tick_params(labelsize=8.5)
ex0=0.585
T(ex0,0.730,"Each read aligns to BOTH copies; the BAM",11,c='#222',w='bold')
T(ex0,0.716,"scores every placement (always present):",11,c='#222',w='bold')
for i,lab in enumerate(["AS   alignment score   (higher = better)",
                        "ms   best-segment score(higher = better)",
                        "de   divergence        (lower  = better)",
                        "NM   edit distance     (lower  = better)"]):
    T(ex0,0.696-i*0.016,lab,9.6,mono=True,c='#222')
T(ex0,0.630,"Ambiguity = AS ratio (worse/better) of the two placements",8.8,c='#555',st='italic')
T(ex0,0.620,"= minimap2 s2/s1, but s2 is primary-only so we use AS.",8.8,c='#555',st='italic')
T(ex0,0.600,"OFF the diagonal  ->  fits one copy better",11,c=GREEN,w='bold')
T(ex0,0.585,"->  ASSIGN (confidently if the gap is big).",11,c=GREEN,w='bold')
T(ex0,0.560,"NEAR the diagonal  ->  fits both almost",11,c=RED,w='bold')
T(ex0,0.545,"equally  ->  low / no confidence.",11,c=RED,w='bold')
T(ex0,0.522,"R1: NM 5 vs 626 -> DAZ1.  R3: NM 0=0 -> tied.",9.6,c='#333',st='italic')
T(0.05,0.470,"Most points hug the diagonal: the copies are so similar that the per-read fit barely differs. Only a tail (e.g. R1, a "
  "609 bp insertion on DAZ3) is decisive.",10,c='#333',st='italic')

# ===================== Panel 3: confidence spectrum =====================
panel(0.03,0.232,0.94,0.205,'#ffffff')
T(0.05,0.430,"3.  How resolvable are the reads, really?  A CONFIDENCE SPECTRUM (AS ratio of the two placements — present on every read)",13,w='bold')
# stacked spectrum bar
bx0,bx1=0.06,0.93; by=0.388; bh=0.026
segs=[(n_conf,GREEN,'confident\n(AS ratio < 0.95)'),(n_weak,GOLD,'low confidence\n(0.95-0.999)'),(n_tied,RED,'fully tied\n(AS ratio ~ 1)')]
cur=bx0
for cnt,col,lab in segs:
    w=(cnt/N)*(bx1-bx0)
    ax.add_patch(Rectangle((cur,by),w,bh,fc=col,ec='white',lw=1.2,zorder=3))
    T(cur+w/2,by+bh/2,f"{cnt}",11,ha='center',va='center',c='white',w='bold')
    T(cur+w/2,by-0.006,lab,8.6,ha='center',va='top',c=col,w='bold')
    cur+=w
T(0.06,0.345,f"Only {n_conf}/{N} reads are confidently assignable. {n_weak} are low-confidence (DAZ3 scores >=95% as well as DAZ1) and "
  f"{n_tied} are fully tied.",10.3,c='#333',w='bold')
T(0.06,0.330,"So yes — far more than the fully-tied set is uncertain. That is expected for 99.97%-identical copies.",10.3,c='#333',st='italic')
# what the EM does
panel(0.055,0.247,0.89,0.073,'#f0f6ff',ec=BLUE,lw=1.2)
T(0.07,0.312,"The EM never hard-assigns. Each read gets a SOFT fractional weight =",10.2,c=BLUE,w='bold')
T(0.07,0.298,"   (per-read fit likelihood from NM/de + junction compatibility)  x  anchored prior,   normalized to sum to 1 across copies.",9.6,mono=True,c='#222')
T(0.07,0.283,"Confident reads -> ~1 on the better copy. Low-confidence & tied reads -> apportioned by the ANCHORED PRIOR, which holds them",9.6,c='#333')
T(0.07,0.270,"toward DAZ1 (it carries the unique-read evidence). That is why DAZ3 settles at ~4, not an inflated 50/50 share.",9.6,c='#333')

# ===================== Summary =====================
panel(0.03,0.012,0.94,0.205,'#f4f6fb',ec=BLUE,lw=1.8)
T(0.05,0.205,"Summary (corrected)",13,w='bold',c=BLUE)
pts=[
 ("-  Evidence is ON THE READ: every read is scored against BOTH copies (AS / de / NM / ms, present on every placement). We assign",1),
 ("   from how the read fits each copy, not from where the reference loci differ.",0),
 (f"-  Resolvability is a GRADIENT: {n_conf} confident, {n_weak} low-confidence, {n_tied} fully tied (of {N}). The EM apportions ALL of them",1),
 ("   softly (fractional weights) — nothing is hard-assigned; the anchored prior keeps low-confidence reads off DAZ3 (cov 163 -> ~4).",0),
 ("-  The 7 locus SNVs do NOT drive this: 6/7 are intronic (a spliced read never covers them); only 1 is exonic.",1),
 ("   rustle's VG EM uses alignment identity (NM/de) + junction-chain compatibility + the anchored prior; the SNV term is neutral.",0),
 ("-  The 10.8 kb block is a LOCUS property; on a read it just appears as a large NM/de gap (e.g. R1's 609 bp insertion on DAZ3).",1),
]
for k,(p,bold) in enumerate(pts):
    T(0.05,0.183-k*0.0215,p,10.4,c=('#222' if bold else '#444'),w=('bold' if bold else 'normal'))

fig.savefig("summary_daz_tied_corrected.png",dpi=150,bbox_inches='tight')
print(f"wrote summary_daz_tied_corrected.png  (tiers: confident={n_conf} weak={n_weak} tied={n_tied})")
