#!/usr/bin/env python3
"""DAZ3 abundance under three multimapper schemes, from the real BAM:
(A) StringTie uniform 1/n, (B) include secondaries with no correction (double-count),
(C) mass-conserving EM (anchored-prior apportionment). Shows 1/n over-counts ~50x and
naive inclusion ~100x vs the ~2-6 sequence-supported truth; EM lands at it. mathtext-safe.
"""
import re, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

DAZ1=(42783133,42859657); DAZ3=(42879918,42945552)
def ref_len(c): return sum(int(n) for n,o in re.findall(r'(\d+)([MIDNSH=X])',c) if o in 'MDN=X')
def copy_of(p,e):
    m=(p+e)//2
    return 'DAZ1' if DAZ1[0]<=m<=DAZ1[1] else 'DAZ3' if DAZ3[0]<=m<=DAZ3[1] else None
reads=collections.defaultdict(dict)
for line in open('/tmp/daz_aln.sam'):
    f=line.split('\t'); c=f[5]
    if c=='*': continue
    p=int(f[3]); cp=copy_of(p,p+ref_len(c))
    if cp is None: continue
    AS=next((int(t[5:]) for t in f[11:] if t.startswith('AS:i:')),-1)
    if cp not in reads[f[0]] or AS>reads[f[0]][cp]: reads[f[0]][cp]=AS
uniq1=sum(1 for c in reads.values() if 'DAZ1' in c and 'DAZ3' not in c)
uniq3=sum(1 for c in reads.values() if 'DAZ3' in c and 'DAZ1' not in c)
both={q:c for q,c in reads.items() if 'DAZ1' in c and 'DAZ3' in c}
N=len(both); fav1=fav3=tied=0
for c in both.values():
    d=c['DAZ1']-c['DAZ3']
    if d>5: fav1+=1
    elif d<-5: fav3+=1
    else: tied+=1
anch1=uniq1+fav1; anch3=uniq3+fav3; pi3=anch3/max(1,anch1+anch3)
A=0.5*N; B=1.0*N; C=fav3*1.0+tied*pi3

BLUE,GREEN,RED,ORANGE,GREY='#2c7fb8','#1a9850','#d73027','#f16913','#555'
fig=plt.figure(figsize=(14,8))
axb=fig.add_axes([0.07,0.12,0.52,0.74])
bars=['(A) StringTie\nuniform 1/n','(B) include 2ndary,\nno correction','(C) mass-conserving\nEM (apportion)']
vals=[A,B,C]; cols=[ORANGE,RED,GREEN]
x=range(3)
axb.bar(x,vals,color=cols,width=0.6,zorder=3,edgecolor='#333',linewidth=1)
axb.set_yscale('log'); axb.set_ylim(0.5,400)
axb.axhspan(2,6,color='#888',alpha=0.25,zorder=1)
axb.text(2.45,4,'ground truth\n~2-6 reads',ha='left',va='center',fontsize=10,color='#333',style='italic')
for i,v in enumerate(vals):
    axb.text(i,v*1.12,f"{v:.1f}" if v>=10 else f"{v:.2f}",ha='center',va='bottom',
             fontsize=13,weight='bold',color=cols[i])
axb.set_xticks(list(x)); axb.set_xticklabels(bars,fontsize=10.5)
axb.set_ylabel("DAZ3 apportioned read-mass  (log scale)",fontsize=11.5)
axb.set_title("DAZ3 abundance under three multimapper schemes (215 reads on both copies)",
              fontsize=12.5,weight='bold')
axb.spines[['top','right']].set_visible(False)
axb.annotate('~50x too high',xy=(0,A),xytext=(0,250),ha='center',fontsize=9.5,color=ORANGE,
             arrowprops=dict(arrowstyle='->',color=ORANGE))
axb.annotate('~100x too high\n(= double-count;\nrustle default --vg sits here\nfor DAZ today)',
             xy=(1,B),xytext=(1.0,300),ha='center',fontsize=9,color=RED,
             arrowprops=dict(arrowstyle='->',color=RED))
axb.annotate('matches truth',xy=(2,C),xytext=(2,30),ha='center',fontsize=9.5,color=GREEN,
             arrowprops=dict(arrowstyle='->',color=GREEN))

# right-side text
ax=fig.add_axes([0,0,1,1]); ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
def txt(x,y,s,sz=11,c='black',w='normal',st='normal',mono=False):
    ax.text(x,y,s,size=sz,ha='left',va='top',color=c,weight=w,style=st,
            family='monospace' if mono else 'sans-serif')
def box(x,y,w,h,fc,ec): ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.004",fc=fc,ec=ec,lw=1.4))
txt(0.5,0.965,"Why uniform 1/n and uncorrected inclusion both fail; EM fixes it",16,w='bold')
box(0.63,0.55,0.345,0.30,'#fbfbfb','#bbb')
txt(0.645,0.835,"The three schemes",12,w='bold')
txt(0.645,0.805,"A  1/n: split every read evenly. Ignores",10,c=ORANGE,w='bold')
txt(0.645,0.787,"     the evidence AND the copies' expression.",9.5,c='#333')
txt(0.645,0.762,"B  include 2ndary at full weight: each read",10,c=RED,w='bold')
txt(0.645,0.744,"     counts fully on BOTH copies -> mass>1,",9.5,c='#333')
txt(0.645,0.727,"     double-count, fabricates the copy.",9.5,c='#333')
txt(0.645,0.700,"C  EM: each read distributes fractional",10,c=GREEN,w='bold')
txt(0.645,0.682,"     weight summing to 1 across copies;",9.5,c='#333')
txt(0.645,0.665,"     gamma = likelihood x relative-abundance",9.5,c='#333')
txt(0.645,0.648,"     prior. Anchored reads set the prior;",9.5,c='#333')
txt(0.645,0.631,"     tied reads follow it.",9.5,c='#333')
txt(0.645,0.598,f"anchored: DAZ1={anch1}  DAZ3={anch3}  ->  prior pi_DAZ3={pi3:.3f}",9.5,c='#222',mono=True,w='bold')

box(0.63,0.12,0.345,0.40,'#eef0f7','#2c7fb8')
txt(0.645,0.505,"My read on your advisor's point",12,w='bold',c=BLUE)
notes=[
 ("He's right that 1/n is wrong.",True),
 ("1/n is the zero-information estimator; it",False),
 ("ignores both the per-read evidence and the",False),
 ("relative-abundance prior. The field (RSEM,",False),
 ("Salmon, Cufflinks) uses EM for exactly this.",False),
 ("",False),
 ("His qualm about uncorrected inclusion is",True),
 ("also right -- and it's WORSE than 1/n:",False),
 ("full inclusion double-counts (215 vs truth).",False),
 ("",False),
 ("But: include secondaries FOR STRUCTURE",True),
 ("(you can't assemble a starved copy without",False),
 ("the reads); apportion their weight (EM,",False),
 ("mass-conserving) FOR ABUNDANCE; feed the",False),
 ("apportioned coverage to the flow.",False),
]
yy=0.475
for s,b in notes:
    txt(0.645,yy,s,9.5,c=('#222' if b else '#444'),w=('bold' if b else 'normal'))
    yy-=0.0205

fig.savefig("em_abundance_comparison.png",dpi=150,bbox_inches='tight')
print(f"wrote em_abundance_comparison.png | A={A} B={B} C={C:.2f} pi3={pi3:.4f} N={N} fav1={fav1} fav3={fav3} tied={tied}")
