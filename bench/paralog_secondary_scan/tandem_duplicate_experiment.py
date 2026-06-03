#!/usr/bin/env python3
"""Are near-identical tandem duplicates merged / losing transcripts? Experiment summary.

Synthetic tandem: two gene copies (3 exons each) at a 2 kb gap, identity and
per-copy expression varied; reads simulated (full-length, HiFi error), aligned,
assembled with rustle (de-novo and --vg). Data: /tmp/tandem_synth_results.json.
"""
import json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

R=json.load(open('/tmp/tandem_synth_results.json'))
BLUE,ORANGE,GREEN,RED,GREY='#2c7fb8','#f16913','#1a9850','#d73027','#666'

fig=plt.figure(figsize=(14.5,11))
ax=fig.add_axes([0,0,1,1]); ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
def T(x,y,s,sz=11,ha='left',va='top',c='black',w='normal',mono=False,st='normal'):
    ax.text(x,y,s,size=sz,ha=ha,va=va,color=c,weight=w,style=st,family='monospace' if mono else 'sans-serif',zorder=5)
def panel(x,y,w,h,fc,ec='#bbb',lw=1.2): ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.004",fc=fc,ec=ec,lw=lw))

T(0.5,0.985,"Do near-identical tandem duplicates get merged / lose transcripts? (synthetic test)",17,ha='center',w='bold')
T(0.5,0.963,"Two tandem gene copies (3 exons, 2 kb apart); identity and per-copy expression varied; full-length reads simulated, "
  "aligned, assembled with rustle.",11,ha='center',c=GREY,st='italic')

# ---- Panel A: bar chart of per-copy coverage (true vs rustle) ----
panel(0.03,0.485,0.94,0.45,'#ffffff')
T(0.05,0.925,"A.  Recovered coverage per copy: TRUE (hollow) vs rustle (solid).  copyA = expressed; copyB = the test copy.",12.5,w='bold')
axB=fig.add_axes([0.085,0.55,0.87,0.32])
labels=[f"{r['label']}\n({r['mode']})" for r in R]
x=np.arange(len(R)); w=0.19
trueA=[r['trueA'] for r in R]; covA=[r['covA'] for r in R]
trueB=[r['trueB'] for r in R]; covB=[r['covB'] for r in R]
axB.bar(x-1.5*w,trueA,w,facecolor='none',edgecolor=BLUE,lw=1.4,label='copyA true')
axB.bar(x-0.5*w,covA,w,color=BLUE,label='copyA rustle')
axB.bar(x+0.5*w,trueB,w,facecolor='none',edgecolor=ORANGE,lw=1.4,label='copyB true')
axB.bar(x+1.5*w,covB,w,color=ORANGE,label='copyB rustle')
for i,r in enumerate(R):
    # flag fabrication / drop / split
    if r['covB']>r['trueB']*2+3:
        axB.annotate('fabricated',(i+1.5*w,r['covB']),textcoords='offset points',xytext=(0,4),ha='center',fontsize=8,color=RED,weight='bold')
    if r['mode']=='VG' and r['flag']!='-':
        axB.annotate(f"flag:{r['flag']}",(i,max(r['covA'],r['covB'])),textcoords='offset points',xytext=(0,16),ha='center',fontsize=8,color=GREEN,weight='bold')
    if r['covB']<1 and r['trueB']>1:
        axB.annotate('dropped',(i+1.5*w,1),textcoords='offset points',xytext=(0,4),ha='center',fontsize=8,color=RED,weight='bold')
axB.set_xticks(x); axB.set_xticklabels(labels,fontsize=9)
axB.set_ylabel("coverage",fontsize=10.5); axB.legend(fontsize=9,ncol=4,loc='upper center',bbox_to_anchor=(0.5,1.13))
axB.set_ylim(0,52)
T(0.05,0.505,"Both copies are always assembled as SEPARATE loci (never merged into one, no chimeras at this 2 kb gap). The problem is "
  "ABUNDANCE, not structure.",9.8,c='#333',st='italic')

# ---- Panel B: findings ----
panel(0.03,0.255,0.46,0.215,'#f6f8fb')
T(0.05,0.460,"B.  What happens (structure & abundance)",12,w='bold')
lines=[
 ("STRUCTURE: copies stay separate.",GREEN,True),
 ("Both copies recovered at gap>=200 bp, even at 100%",None,False),
 ("identity; no chimeric transcripts. (gap=0, i.e. copies",None,False),
 ("directly abutting -> 0 transcripts: pathological edge.)",None,False),
 ("ABUNDANCE fabrication only at ~100% identity:",RED,True),
 ("equal expr -> ~correct; UNEQUAL (40 vs 3) de-novo ->",None,False),
 ("copyB fabricated to 26 (true 3). Resolves fully by 99%",None,False),
 ("identity (copyB -> 3.0): ~15 diffs/read is enough.",None,False),
]
y=0.438
for txt,col,bold in lines:
    T(0.055,y,txt,9.2,c=(col or '#222'),w=('bold' if bold else 'normal')); y-=0.0205

# ---- Panel C: VG behaviour + the bug ----
panel(0.51,0.255,0.46,0.215,'#f6f8fb')
T(0.53,0.460,"C.  VG mode, and a real bug it exposed",12,w='bold')
lines2=[
 ("VG at 100% id: flags family_identifiability=\"none\",",GREEN,True),
 ("copy_confidence 0.5, splits 50/50 -> honest, but cannot",None,False),
 ("recover truth (no unique reads to anchor = the floor).",None,False),
 ("VG at 99% id: over-suppresses the low copy (3 -> 0).",ORANGE,True),
 ("BUG: capacity_confidence is stuck at 1.000.",RED,True),
 ("was_unique=(nh<=1) but minimap2 emits NO NH tag ->",None,False),
 ("every read looks unique -> channel never flags. Real",None,False),
 ("DAZ run had the same 1.000. Fix: derive multimap from",None,False),
 ("tp:A:S / family, not the absent NH tag.",None,False),
]
y=0.438
for txt,col,bold in lines2:
    T(0.535,y,txt,9.2,c=(col or '#222'),w=('bold' if bold else 'normal')); y-=0.0205

# ---- Panel D: real-data context ----
panel(0.03,0.02,0.94,0.225,'#eef2f8')
T(0.05,0.235,"D.  Real GGO data: the risk loci exist, but separation mostly holds",12.5,w='bold')
T(0.05,0.212,"-  Genome scan: 397 adjacent near-identical paralog pairs (identity >=90%). The 100%-identical ones are mostly UNEXPRESSED",10,c='#222')
T(0.07,0.193,"tandem arrays (0-2 reads). Expressed near-identical tandems include the protocadherin-beta (PCDHB) cluster.",10,c='#222')
T(0.05,0.168,"-  PCDHB cluster (90-97% identity, single-exon tandem genes): rustle produces ONE transcript per expressed copy at the",10,c='#222')
T(0.07,0.149,"right positions -- it does NOT collapse the cluster -- with a few low-coverage CHIMERIC transcripts spanning adjacent copies.",10,c='#222')
T(0.05,0.122,"Bottom line for the advisor:",10.5,w='bold',c=BLUE)
T(0.05,0.103,"-  We do NOT merge tandem copies into one locus or lose a copy's transcripts (they assemble separately).",10,c='#222')
T(0.05,0.084,"-  The residual risk is ABUNDANCE: at ~100% identity with uneven expression, reads are non-identifiable and the silent copy",10,c='#222')
T(0.07,0.065,"can be over-counted (same wall as DAZ). VG flags this honestly via family_identifiability; the capacity_confidence channel",10,c='#222')
T(0.07,0.046,"needs the NH-tag fix to also flag it. >=1% divergence removes the problem entirely.",10,c='#222')

fig.savefig("tandem_duplicate_experiment.png",dpi=150,bbox_inches='tight')
print("wrote tandem_duplicate_experiment.png")
