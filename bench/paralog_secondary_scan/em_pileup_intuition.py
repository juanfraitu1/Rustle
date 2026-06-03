#!/usr/bin/env python3
"""Pileup-EM intuition, hardened against pushback (panel wf_077aa42b).
Core idea: each read carries ONE vote; EM SPLITS that vote between copies in the
same proportion the DECIDED reads voted (177:1 for DAZ1:DAZ3) -- never copies it
onto both (double-count), and when nothing distinguishes the copies it abstains.
Real DAZ numbers throughout; full-length (FLNC) read geometry; mass-conservation
made visible. mathtext-safe. Output: em_pileup_intuition.png
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

# ---- real DAZ numbers (em_abundance_comparison.py / final_tally.py) ----
N=215; dec1=174; dec3=1; tied=40; uniq1=3; uniq3=0
anch1=dec1+uniq1; anch3=dec3+uniq3          # 177 : 1
pi3=anch3/(anch1+anch3)                       # 0.0056
A_1n=0.5*N; B_incl=1.0*N; C_em=dec3*1.0+tied*pi3   # 107.5 / 215 / 1.22

fig=plt.figure(figsize=(16,13.4)); ax=fig.add_axes([0,0,1,1]); ax.axis('off')
ax.set_xlim(0,1); ax.set_ylim(0,1)
BLUE,GREEN,RED,ORANGE,GREY,SLATE,PURP='#2c7fb8','#1a9850','#d73027','#f16913','#666','#5b6b7a','#7570b3'
def txt(x,y,s,sz=12,ha='left',va='top',c='black',w='normal',mono=False,st='normal'):
    ax.text(x,y,s,size=sz,ha=ha,va=va,color=c,weight=w,style=st,
            family='monospace' if mono else 'sans-serif',zorder=6)
def box(x,y,w,h,fc,ec='#888',lw=1.2,z=1):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.004",fc=fc,ec=ec,lw=lw,zorder=z))
def rect(x,y,w,h,fc,ec='none',lw=0,z=2,hatch=None,alpha=1):
    ax.add_patch(Rectangle((x,y),w,h,fc=fc,ec=ec,lw=lw,zorder=z,hatch=hatch,alpha=alpha))
PAT=[(0.0,0.08),(0.17,0.24),(0.35,0.43),(0.55,0.66),(0.78,0.86),(0.94,1.0)]
def flnc(x0,x1,y,color,exf,h=0.0095,alpha=0.95,ec=None):
    ax.plot([x0,x1],[y,y],color=(ec or color),lw=0.6,alpha=alpha*0.8,zorder=3)
    for s,e in exf:
        rect(x0+s*(x1-x0),y-h/2,(e-s)*(x1-x0),h,color,ec=(ec or color),lw=0.4,z=4,alpha=alpha)

# =================== ZONE 0: header ===================
txt(0.5,0.987,"How the EM splits a multi-mapping read between copies — shown on the WORST case",17,ha='center',w='bold')
txt(0.5,0.966,"Each read = one full-length IsoSeq (FLNC) molecule spanning the whole transcript. DAZ1 & DAZ3 are 99.97% identical, so most reads cover NO position that tells them apart.",
    11,ha='center',c=GREY,st='italic')
txt(0.5,0.949,"This is spectrum row E (zero information). Even here, EM is provably NOT 1/n and CANNOT inflate a copy.",
    11,ha='center',c=PURP if False else SLATE,w='bold')

# =================== ZONE 1: the pileup ===================
box(0.03,0.665,0.94,0.265,'#f7f9fb',ec=BLUE,lw=1.6)
txt(0.05,0.922,"1.  The pileup — real DAZ data (215 reads place on both copies)",13,w='bold',c=BLUE)
LX0,LX1=0.085,0.46
# DAZ1 (minus) + decided blue reads
ax.plot([LX0,LX1],[0.898,0.898],color=BLUE,lw=1.2,zorder=2)
for s,e in PAT: rect(LX0+s*(LX1-LX0),0.898-0.011,(e-s)*(LX1-LX0),0.022,'#cfe2f0',ec=BLUE,lw=1.0,z=3)
txt(LX0-0.006,0.898,"DAZ1 (−)",8.5,ha='right',va='center',w='bold',c=BLUE)
for i in range(5): flnc(LX0,LX1,0.880-0.0115*i,BLUE,PAT)
txt(LX0,0.815,"174 reads DECIDED → DAZ1",8.5,c=BLUE,w='bold')
# tied reads straddling (grey)
for i in range(3): flnc(LX0,LX1,0.792-0.0115*i,GREY,PAT,alpha=0.6)
txt(LX0,0.749,"40 reads TIED — cover no distinguishing position",8.5,c=GREY,w='bold')
# DAZ3 (plus, inverted) + its single decided green read
PAT3=[(1-e,1-s) for s,e in PAT]
ax.plot([LX0,LX1],[0.726,0.726],color=GREEN,lw=1.2,zorder=2)
for s,e in PAT3: rect(LX0+s*(LX1-LX0),0.726-0.011,(e-s)*(LX1-LX0),0.022,'#cde6d7',ec=GREEN,lw=1.0,z=3)
txt(LX0-0.006,0.726,"DAZ3 (+, inv)",8.5,ha='right',va='center',w='bold',c=GREEN)
flnc(LX0,LX1,0.708,GREEN,PAT3)
txt(LX0,0.692,"1 read DECIDED → DAZ3",8.5,c=GREEN,w='bold')
# right column: the 177:1 argument + preempts
RX=0.50
txt(RX,0.900,"DECIDED  (dNM ≥ 2):",11,w='bold',c='#222')
txt(RX+0.02,0.880,"174 → DAZ1      1 → DAZ3      (+ unique 3 / 0)",10.5,mono=True)
txt(RX+0.02,0.861,"⇒  decided ratio  =  177 : 1",11.5,mono=True,w='bold',c=PURP)
txt(RX,0.835,"TIED  (|dNM| < 2):  40 reads",11,w='bold',c=GREY)
txt(RX+0.02,0.816,"cover 0 of the 7 SNVs; AS, de, NM all identical across copies",9.5,c='#333')
box(RX,0.690,0.45,0.110,'#ffffff',ec='#bbb',lw=1.0)
txt(RX+0.012,0.788,"Why this isn't the \"one SNP is trivial\" case:",10,w='bold',c=RED)
txt(RX+0.012,0.769,"the 174 decided reads are NOT one-SNP calls — each accrues ≥ 2",9.3,c='#333')
txt(RX+0.012,0.753,"distinguishing positions (dNM = count, not a rate) over its FULL",9.3,c='#333')
txt(RX+0.012,0.737,"length. The hard part is the 40 residual tied reads + not letting",9.3,c='#333')
txt(RX+0.012,0.721,"them inflate mass — which 1/n and naive inclusion both botch.",9.3,c='#333')
txt(RX+0.012,0.702,"Inverted pair, ts=+ on 80/80 → strand can't disambiguate either.",9.0,c=GREY,st='italic')

# =================== ZONE 2: one-vote triptych ===================
box(0.03,0.405,0.94,0.245,'#fbfbfc',ec=SLATE,lw=1.6)
txt(0.05,0.642,"2.  One TIED read = one vote.  Three ways to handle that vote:",13,w='bold')
txt(0.05,0.623,"\"Undecided reads split in the same proportion the decided reads voted (177:1).\"",10.5,c=SLATE,st='italic')
def token(cx,cy,label,col,sz=10):
    rect(cx-0.018,cy-0.013,0.036,0.026,col,ec='#333',lw=1.0,z=4)
    txt(cx,cy,label,sz,ha='center',va='center',c='white',w='bold')
# (A) double-count
box(0.05,0.420,0.285,0.165,'#fdecea',ec=RED,lw=1.4)
txt(0.06,0.575,"(A) include 2ndary, no correction",10.5,w='bold',c=RED)
token(0.12,0.535,"1.0",BLUE); txt(0.12,0.512,"on DAZ1",8,ha='center',c=BLUE)
token(0.25,0.535,"1.0",GREEN); txt(0.25,0.512,"on DAZ3",8,ha='center',c=GREEN)
txt(0.06,0.487,"the SAME read counted FULLY on both  →  total 2.0",9,c='#333')
txt(0.06,0.470,"mass per read > 1  →  fabricates the copy.",9,c=RED,w='bold')
txt(0.06,0.450,"rustle --vg sits here for DAZ today  →  DAZ3 = 215  ✗",9,c=RED,w='bold')
# (B) EM mass-conserving
box(0.355,0.420,0.285,0.165,'#eef7f0',ec=GREEN,lw=1.4)
txt(0.365,0.575,"(B) EM — split the one vote",10.5,w='bold',c=GREEN)
# a 1.0 bar split 0.994 / 0.006
bx,bw=0.375,0.235
rect(bx,0.533,bw*0.994,0.018,BLUE,ec='#333',lw=0.8,z=4)
rect(bx+bw*0.994,0.533,max(bw*0.006,0.004),0.018,GREEN,ec='#333',lw=0.8,z=4)
txt(bx,0.560,"DAZ1  0.994",8.5,c=BLUE,w='bold')
ax.add_patch(FancyArrowPatch((bx+bw,0.527),(bx+bw+0.012,0.500),arrowstyle='-|>',mutation_scale=9,lw=0.9,color=GREEN,zorder=5))
txt(bx+bw+0.014,0.498,"DAZ3 0.006",8.5,c=GREEN,w='bold')
txt(0.365,0.487,"one vote SPLIT 177:1, never duplicated  →  sum = 1.0",9,c='#333')
txt(0.365,0.470,"mass conserved  →  cannot inflate.",9,c=GREEN,w='bold')
txt(0.365,0.450,"DAZ3 = 1 + 40×0.006 ≈ 1.2   (truth 2–6)  ✓",9,c=GREEN,w='bold')
# (C) abstain
box(0.66,0.420,0.31,0.165,'#eef0f3',ec=SLATE,lw=1.4)
txt(0.67,0.575,"(C) zero-anchor → ABSTAIN",10.5,w='bold',c=SLATE)
token(0.715,0.535,"?",SLATE,12)
txt(0.67,0.508,"no copy-specific read  →  decided ratio = 0/0",9,c='#333')
txt(0.67,0.491,"(undefined).  Nothing to split by.",9,c='#333')
txt(0.67,0.470,"The method ABSTAINS — reports aggregate +",9,c=SLATE,w='bold')
txt(0.67,0.453,"flags unresolvable.  It does NOT invent a split.",9,c=SLATE,w='bold')
txt(0.67,0.432,"(spectrum row E, fully non-identifiable)",8.3,c=GREY,st='italic')

# =================== ZONE 3: payoff bar ===================
box(0.03,0.215,0.94,0.180,'#fbf7ee',ec='#caa',lw=1.4)
txt(0.05,0.387,"3.  DAZ3 read-mass under each scheme  (the whole argument, on real data)",13,w='bold')
axb=fig.add_axes([0.085,0.245,0.55,0.115])
bars=['1/n (uniform)','include 2ndary\n(no correction)','EM\n(split + conserve)']
vals=[A_1n,B_incl,C_em]; cols=[ORANGE,RED,GREEN]
axb.barh(range(3),vals,color=cols,edgecolor='#333',height=0.62,zorder=3)
axb.set_xscale('log'); axb.set_xlim(0.7,400); axb.set_yticks(range(3)); axb.set_yticklabels(bars,fontsize=9)
axb.invert_yaxis(); axb.axvspan(2,6,color='#888',alpha=0.22,zorder=1)
for i,v in enumerate(vals):
    axb.text(v*1.12,i,(f"{v:.1f}" if v>=10 else f"{v:.2f}"),va='center',fontsize=10,weight='bold',color=cols[i])
axb.text(3.4,2.62,'truth ~2–6',fontsize=8,color='#333',style='italic',ha='center')
axb.set_xlabel("DAZ3 read-mass (log scale)",fontsize=9); axb.tick_params(labelsize=8)
for sp in ['top','right']: axb.spines[sp].set_visible(False)
txt(0.66,0.360,"1/n = 107.5  (~50×): ignores the 177:1 evidence.",9.5,c=ORANGE,w='bold')
txt(0.66,0.341,"include 2ndary = 215  (~100×): double-count,",9.5,c=RED,w='bold')
txt(0.675,0.324,"total mass across copies = 430 > N=215.",9.0,c=RED)
txt(0.66,0.302,"EM = 1.2: follows evidence AND conserves mass",9.5,c=GREEN,w='bold')
txt(0.675,0.285,"(DAZ1 213.8 + DAZ3 1.2 = 215).",9.0,c=GREEN)
txt(0.66,0.260,"DAZ3 ≈ 1 of 215 is prior-driven, high-variance",9.0,c='#333',st='italic')
txt(0.66,0.244,"(1 decisive read): \"essentially unexpressed,\"",9.0,c='#333',st='italic')
txt(0.66,0.228,"consistent with 2–6 — NOT a recovered copy.",9.0,c='#333',st='italic')

# =================== ZONE 4: two jobs (replaces self-defeating footer) ===================
box(0.03,0.02,0.455,0.180,'#eef3f7',ec=BLUE,lw=1.5)
txt(0.05,0.192,"Structure — deterministic (EM does NOT create the copy)",11,w='bold',c=BLUE)
txt(0.05,0.168,"Secondary alignments are retained at ingest (bundle.rs:1422).",9.5,c='#333')
txt(0.05,0.151,"The assembler builds DAZ3's exon-junction structure from them",9.5,c='#333')
txt(0.05,0.134,"with or WITHOUT the EM. The reads carry both copies' shared",9.5,c='#333')
txt(0.05,0.117,"structure (96.6% of junctions shared).",9.5,c='#333')
txt(0.05,0.090,"So the EM is not what recovers the copy —",9.5,c=BLUE,w='bold')
txt(0.05,0.073,"it owns the next job:",9.5,c=BLUE,w='bold')
box(0.515,0.02,0.455,0.180,'#fdf3ec',ec=ORANGE,lw=1.5)
txt(0.535,0.192,"Abundance — apportioned (the job EM owns)",11,w='bold',c=ORANGE)
txt(0.535,0.168,"Getting the split wrong (107.5 or 215 vs 1.2) fabricates the",9.5,c='#333')
txt(0.535,0.151,"copy at the QUANTIFICATION level — a phantom isoform with",9.5,c='#333')
txt(0.535,0.134,"50–100× its real coverage.",9.5,c='#333')
txt(0.535,0.107,"EM's contribution: follow the decided ratio, conserve mass,",9.5,c=ORANGE,w='bold')
txt(0.535,0.090,"and abstain when there is no ratio — so a copy is never",9.5,c=ORANGE,w='bold')
txt(0.535,0.073,"inflated or invented.",9.5,c=ORANGE,w='bold')
txt(0.535,0.045,"(On informative families — rows B–D, reads cover diagnostic",8.3,c=GREY,st='italic')
txt(0.535,0.031,"sites — the per-read fit moves the split across passes; DAZ is the row-E worst case.)",8.3,c=GREY,st='italic')

fig.savefig("em_pileup_intuition.png",dpi=150,bbox_inches='tight')
print(f"wrote em_pileup_intuition.png | 177:1 pi3={pi3:.4f} A={A_1n} B={B_incl} C={C_em:.2f}")
