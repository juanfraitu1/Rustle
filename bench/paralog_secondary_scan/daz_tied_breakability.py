#!/usr/bin/env python3
"""Can the tied DAZ reads be broken with the BAM tags? Evidence figure.

The 39 fully-tied reads are byte-identical on all 9 alignment tags (provably
non-identifiable: they cover sequence identical in both copies). The 84
low-confidence reads are NOT truly tied (all 84 differ on >=1 tag -> sharpenable
by combining tags in the EM). The anchored prior contains the unbreakable 39.
Data: /tmp/daz_tied.json (built by /tmp/daz_tied_data.py).
"""
import json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

J=json.load(open('/tmp/daz_tied.json'))
N=J['N']; n_tied=J['n_tied']; n_weak=J['n_weak']; n_conf=J['n_conf']
TAGS=J['TAGS']; tagdiff=J['tagdiff']
exq,exc=J['examples'][0]; exq_short=exq.split('/')[-2] if '/' in exq else exq
BLUE,GREEN,RED,ORANGE,GREY,GOLD='#2c7fb8','#1a9850','#d73027','#f16913','#666','#e8a200'

fig=plt.figure(figsize=(14.5,15))
ax=fig.add_axes([0,0,1,1]); ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
def T(x,y,s,sz=11,ha='left',va='top',c='black',w='normal',mono=False,st='normal'):
    ax.text(x,y,s,size=sz,ha=ha,va=va,color=c,weight=w,style=st,
            family='monospace' if mono else 'sans-serif',zorder=5)
def panel(x,y,w,h,fc,ec='#bbb',lw=1.2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.004",fc=fc,ec=ec,lw=lw))

T(0.5,0.988,"Can the tied DAZ reads be broken with the alignment tags?  No — and that is provable",18,ha='center',w='bold')
T(0.5,0.968,f"Of {N} reads on both copies, {n_tied} are 'fully tied'. On those, every BAM tag is identical between DAZ1 and DAZ3 — "
  "there is no signal on the molecule to break.",11.5,ha='center',c=GREY,st='italic')

# ===================== Panel A: 9-tag table with a worked example =====================
panel(0.03,0.640,0.94,0.305,'#ffffff')
T(0.05,0.937,f"A.  The {n_tied} fully-tied reads are byte-identical on ALL 9 alignment tags (worked example: read {exq_short})",13,w='bold')
cx=[0.06,0.165,0.55,0.70,0.835]
T(cx[0],0.905,"tag",11,w='bold',c='#222',mono=True)
T(cx[1],0.905,"what it is",11,w='bold',c='#222')
T(cx[2],0.905,"DAZ1",11,w='bold',c=BLUE); T(cx[3],0.905,"DAZ3",11,w='bold',c=GREEN)
T(cx[4],0.905,f"differ / {n_tied}",11,w='bold',c='#222')
def val(rec,t):
    v=rec[t]
    return f"{v:.4f}" if t=='de' else str(v)
for i,(t,desc) in enumerate(TAGS):
    y=0.882-i*0.0245
    same=(tagdiff[t]==0)
    T(cx[0],y,t,10.5,mono=True,w='bold',c='#222')
    T(cx[1],y,desc,10,c='#333')
    T(cx[2],y,val(exc['DAZ1'],t),10.5,mono=True,c='#222')
    T(cx[3],y,val(exc['DAZ3'],t),10.5,mono=True,c='#222')
    T(cx[4],y,f"{tagdiff[t]} / {n_tied}",10.5,mono=True,w='bold',c=(GREEN if same else RED))
T(0.05,0.660,f"Every value in the DAZ1 and DAZ3 columns is identical, for all {n_tied} reads (the 'differ' column is 0 for every tag). "
  "MAPQ is 0 on both — minimap2 also calls them ambiguous.",9.6,c='#333',st='italic')

# ===================== Panel B: why =====================
panel(0.03,0.380,0.94,0.250,'#f6f8fb')
T(0.05,0.622,"B.  Why no tag can EVER break them: a tied read covers sequence that is identical in both copies",13,w='bold')
X0,X1=0.10,0.66
def bar(y,col,fc,lab,blk):
    ax.add_patch(Rectangle((X0,y),X1-X0,0.020,fc=fc,ec=col,lw=1.1))
    if blk: ax.add_patch(Rectangle((0.50,y-0.002),0.085,0.024,fc=ORANGE,ec='#a33',lw=1.0))
    T(X0-0.008,y+0.010,lab,10,ha='right',va='center',w='bold',c=col)
bar(0.588,BLUE,'#cfe2f0',"DAZ1",True)
bar(0.560,GREEN,'#cde6d7',"DAZ3",False)
T(0.5425,0.610,"copy-specific (DAZ1-only block)",7.5,ha='center',va='bottom',c=ORANGE,w='bold')
T(0.27,0.547,"shared backbone — identical sequence in both copies",8.5,ha='center',va='top',c=GREY,st='italic')
# tied read on shared
ax.add_patch(Rectangle((0.16,0.524),0.18,0.012,fc=RED,ec='none',alpha=0.8))
T(0.25,0.513,"a TIED read sits only on shared sequence",8.6,ha='center',va='top',c=RED,w='bold')
ax.add_patch(FancyArrowPatch((0.34,0.530),(0.67,0.562),arrowstyle='-|>',mutation_scale=12,color='#999',lw=1.1))
T(0.685,0.575,"same bases in both copies  ->  minimap2 builds the SAME",9.5,c='#222',w='bold')
T(0.685,0.561,"alignment to each  ->  AS, ms, NM, de, ... all identical",9.5,c='#222',w='bold')
T(0.685,0.547,"->  AS-ratio = 1.00  ->  non-identifiable.",9.5,c=RED,w='bold')
# resolvable read on copy-specific
ax.add_patch(Rectangle((0.47,0.470),0.14,0.012,fc=GREEN,ec='none',alpha=0.85))
T(0.54,0.458,"a RESOLVABLE read overlaps copy-specific sequence",8.6,ha='center',va='top',c=GREEN,w='bold')
T(0.685,0.470,"different bases  ->  the wrong copy costs mismatches /",9.5,c='#222')
T(0.685,0.457,"an insertion  ->  tags differ  ->  AS-ratio < 1  ->  ASSIGN.",9.5,c=GREEN,w='bold')
T(0.05,0.398,"Phasing/linkage cannot rescue the tied reads either: they cover ZERO copy-distinguishing positions, so there is nothing to link on.",9.6,c='#333',st='italic')

# ===================== Panel C: what this means / what reduces ambiguity =====================
panel(0.03,0.020,0.94,0.350,'#ffffff')
T(0.05,0.362,"C.  So how much is truly unbreakable, and what CAN reduce the ambiguity?",13,w='bold')
# spectrum
bx0,bx1=0.06,0.93; by=0.318; bh=0.028
segs=[(n_conf,GREEN,'confident'),(n_weak,GOLD,'low confidence'),(n_tied,RED,'fully tied')]
cur=bx0
for cnt,col,lab in segs:
    w=(cnt/N)*(bx1-bx0)
    ax.add_patch(Rectangle((cur,by),w,bh,fc=col,ec='white',lw=1.2))
    T(cur+w/2,by+bh/2,f"{cnt}",12,ha='center',va='center',c='white',w='bold')
    T(cur+w/2,by-0.007,lab,9,ha='center',va='top',c=col,w='bold')
    cur+=w
T(0.06,0.272,f"Only {n_tied} of {N} reads are truly unbreakable. The decision is not 'tied vs resolved' — it is a spectrum.",10.5,c='#222',w='bold')
def card(x,w,y,h,title,tc,lines):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.005",fc='#f7f9fc',ec=tc,lw=1.4))
    T(x+0.012,y+h-0.012,title,10.3,w='bold',c=tc)
    for i,l in enumerate(lines):
        T(x+0.012,y+h-0.034-i*0.0185,l,8.9,c='#222')
card(0.05,0.285,0.045,0.205,f"{n_tied} fully tied  ->  apportion",RED,
     ["Identical on every tag; no per-read","signal exists. The EM does NOT","guess: the ANCHORED PRIOR sends","their mass to the copy that other,","resolvable reads anchor (DAZ1).","-> no DAZ3 inflation, no fabrication."])
card(0.355,0.285,0.045,0.205,f"{n_weak} low-confidence  ->  sharpen",GOLD,
     [f"NOT truly tied: all {n_weak}/{n_weak} differ on",">=1 tag (small AS/NM/de gaps).","No single tag decides them, but","COMBINING tags into one likelihood","(AS+de+ms+junction compat) firms","them up -- exactly the EM E-step."])
card(0.66,0.285,0.045,0.205,f"{n_conf} confident  +  other families",GREEN,
     ["Clear winner from the tags.","To break ties in OTHER paralog","families, a read covering even one","PSV + cross-read LINKAGE resolves","it. DAZ's 39 cover none -> they are","the worst-case floor, not a tooling gap."])
T(0.05,0.038,"Bottom line: the 39 are provably non-identifiable (a property of the sequence, not a missing feature); everything near them is sharpenable "
  "by using the tags jointly; the prior keeps the unbreakable core from doing harm.",9.4,c='#333',st='italic')

fig.savefig("daz_tied_breakability.png",dpi=150,bbox_inches='tight')
print("wrote daz_tied_breakability.png")
