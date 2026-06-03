#!/usr/bin/env python3
"""Presentation figure: the operational definition of a multi-copy gene family
(rustle classify_family), with an EXAMPLE that satisfies it (NBPF) and a
COUNTEREXAMPLE that fails it (SORD + a spillover paralog). Reads are drawn as
FULL-LENGTH IsoSeq (FLNC) alignments — each spans the whole transcript, so it
covers many distinguishing positions (the long-read advantage for anchoring).
Grounded in docs/superpowers/specs/2026-06-01-multicopy-family-definition-design.md.
mathtext-safe. Output: family_definition_figure.png
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

fig=plt.figure(figsize=(15.5,10.5)); ax=fig.add_axes([0,0,1,1]); ax.axis('off')
ax.set_xlim(0,1); ax.set_ylim(0,1)
BLUE,GREEN,RED,ORANGE,PURP,GREY='#2c7fb8','#1a9850','#d73027','#f16913','#7570b3','#555'
def txt(x,y,s,sz=12,ha='left',va='top',c='black',w='normal',mono=False,st='normal'):
    ax.text(x,y,s,size=sz,ha=ha,va=va,color=c,weight=w,style=st,
            family='monospace' if mono else 'sans-serif',zorder=5)
def box(x,y,w,h,fc,ec='#888',lw=1.2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.004",fc=fc,ec=ec,lw=lw))

PAT=[(0.0,0.07),(0.16,0.22),(0.33,0.40),(0.52,0.62),(0.74,0.81),(0.92,1.0)]

def model(x0,x1,y,color,fc,h=0.020):
    """An annotated copy (gene model): thin intron line + filled exon boxes."""
    ax.plot([x0,x1],[y,y],color=color,lw=1.0,zorder=2)
    for s,e in PAT:
        ax.add_patch(Rectangle((x0+s*(x1-x0),y-h/2),(e-s)*(x1-x0),h,fc=fc,ec=color,lw=1.0,zorder=3))

def flnc(x0,x1,y,color,h=0.011,hatch=None,alpha=0.95,ec=None):
    """A FULL-LENGTH (FLNC) read: spans the whole transcript — intron line +
    exon blocks across the entire model span (5'->3'), not a fragment."""
    ax.plot([x0,x1],[y,y],color=(ec or color),lw=0.6,alpha=alpha*0.8,zorder=3)
    for s,e in PAT:
        ax.add_patch(Rectangle((x0+s*(x1-x0),y-h/2),(e-s)*(x1-x0),h,
                     fc=color,ec=(ec or color),lw=0.5,alpha=alpha,hatch=hatch,zorder=4))

# ---------------- title ----------------
txt(0.5,0.984,"What is a multi-copy gene family?  (the operational definition)",19,ha='center',w='bold')
txt(0.5,0.961,"A connected component of the read-multi-mapping graph that classify_family accepts. Reads are full-length IsoSeq (FLNC): each spans the whole transcript.",
    11,ha='center',c=GREY,st='italic')

# ---------------- DEFINITION box ----------------
box(0.03,0.645,0.94,0.295,'#f4f6fb',ec=BLUE,lw=1.8)
txt(0.05,0.930,"Definition  — a family must satisfy all of:",13.5,w='bold',c=BLUE)
defs=[
 ("(M) multiplicity","≥ 2 copies (member loci / path-bundles)."),
 ("(H) homology / connectivity","copies share ≥ 1 exon-class, i.e. linked by multi-mapping reads.  [primary structural gate]"),
 ("(X) expression","≥ 2 copies INDEPENDENTLY expressed: each owns ≥ k=3 full-length reads ANCHORED to it (below)."),
 ("(scope) spliced","each copy has ≥ 1 intron (single-exon paralogs out of scope)."),
]
yy=0.906
for tag,desc in defs:
    txt(0.06,yy,tag,11,w='bold',c='#222'); txt(0.275,yy,desc,10.3,c='#333'); yy-=0.024

# ---- worked anchoring mini-diagram (what NM / dNM mean) ----
txt(0.05,0.800,"How the (X) anchor decides a read's copy — edit distance over the FULL read, not a single SNP:",
    10.5,w='bold',c=PURP)
# 12-column toy alignment; distinguishing columns + 1 error column
seqC  = list("ACGTACGTACTG")
dist  = [2,5,8]                 # copy-distinguishing columns (fixed paralog differences)
errc  = 4                       # a sequencing error in the read
seqCp = seqC[:];
for d,b in zip(dist,"TGT"): seqCp[d]=b
read  = seqC[:]; read[errc]='C'  # error: differs from BOTH copies here
gx0,gx1=0.135,0.46; n=12; dx=(gx1-gx0)/n
def cx(i): return gx0+dx*(i+0.5)
rC,rCp,rR=0.778,0.760,0.742
# column shading
for i in range(n):
    if i in dist: ax.add_patch(Rectangle((cx(i)-dx*0.5,rR-0.010),dx,0.046,fc='#ece7f5',ec='none',zorder=1))
    if i==errc:   ax.add_patch(Rectangle((cx(i)-dx*0.5,rR-0.010),dx,0.046,fc='#fdecea',ec='none',zorder=1))
txt(0.13,rC,"copy c",8.5,ha='right',va='center',w='bold',c=BLUE)
txt(0.13,rCp,"copy c'",8.5,ha='right',va='center',w='bold',c=GREEN)
txt(0.13,rR,"read",8.5,ha='right',va='center',w='bold',c='#222')
for i in range(n):
    txt(cx(i),rC,seqC[i],8.5,ha='center',va='center',mono=True,c=BLUE)
    txt(cx(i),rCp,seqCp[i],8.5,ha='center',va='center',mono=True,c=GREEN)
    rc = RED if i==errc else (BLUE if i in dist else '#222')
    txt(cx(i),rR,read[i],8.5,ha='center',va='center',mono=True,c=rc,w=('bold' if (i in dist or i==errc) else 'normal'))
txt(0.135,0.716,"purple = copy-distinguishing positions   ·   red = a sequencing error",8,c=GREY,st='italic')
txt(0.135,0.700,"the read carries copy c's allele at all 3 distinguishing positions",8,c=BLUE)
# NM tally box
box(0.50,0.700,0.21,0.092,'#ffffff',ec='#bbb',lw=1.0)
txt(0.512,0.785,"NM = edit distance",9.5,w='bold')
txt(0.512,0.766,"NM(c)  = 1   (just the error)",9,mono=True,c=BLUE)
txt(0.512,0.749,"NM(c') = 4   (error + 3 diffs)",9,mono=True,c=GREEN)
txt(0.512,0.730,"dNM = 4-1 = 3  >= T(2)",9,mono=True,w='bold',c=PURP)
txt(0.512,0.712,"=> read ANCHORED to c",9,mono=True,w='bold',c=BLUE)
txt(0.50,0.690,"The error hits BOTH copies -> cancels in dNM; what remains = the COUNT of distinguishing",8.3,c='#333')
txt(0.50,0.675,"positions (raw count, not per-base rate). A full-length read covers many -> decisive.",8.3,c='#333')
# reported-per-family compact box (right)
box(0.725,0.700,0.235,0.092,'#ffffff',ec='#bbb',lw=1.0)
txt(0.737,0.785,"Also reported (not gates):",9.3,w='bold',c='#222')
txt(0.737,0.767,"- identifiability:",8.6,c='#333')
txt(0.737,0.752,"  full / partial / none",8.6,c='#333')
txt(0.737,0.735,"  (non-identifiable but expressed",8.2,c=GREY,st='italic')
txt(0.737,0.722,"   = still FAMILY, flagged)",8.2,c=GREY,st='italic')
txt(0.737,0.705,"- arrangement: tandem/segmental",8.6,c='#333')

# ================= EXAMPLE (left) =================
box(0.03,0.05,0.46,0.585,'#eef7f0',ec=GREEN,lw=1.8)
txt(0.05,0.622,"EXAMPLE  -  satisfies the definition",12.5,w='bold',c=GREEN)
txt(0.05,0.601,"NBPF (neuroblastoma breakpoint family): M=23 copies, segmental.",9.8,c='#333')
xa,xb=0.10,0.45
model(xa,xb,0.572,BLUE,'#cfe2f0'); txt(0.095,0.572,"copy A",8.5,ha='right',va='center',w='bold',c=BLUE)
for i in range(3): flnc(xa,xb,0.553-0.013*i,BLUE)
txt(0.10,0.510,"3 full-length FLNC reads anchored to A (dNM>=2 toward A)",8,c=BLUE)
model(xa,xb,0.486,GREEN,'#cde6d7'); txt(0.095,0.486,"copy B",8.5,ha='right',va='center',w='bold',c=GREEN)
for i in range(3): flnc(xa,xb,0.467-0.013*i,GREEN)
txt(0.10,0.424,"3 full-length FLNC reads anchored to B (dNM>=2 toward B)",8,c=GREEN)
for i in range(2): flnc(xa,xb,0.405-0.012*i,GREY,alpha=0.65)
txt(0.10,0.378,"shared reads (|dNM|<T): apportioned, don't decide membership",8,c=GREY,st='italic')
txt(0.06,0.352,"Criteria:",10.5,w='bold')
chk=[("(M)","23 copies  >= 2"),("(H)","shared exon-classes, multi-mapped"),
     ("(X)","A and B each own >= 3 anchored reads => 2 expressed classes"),
     ("(scope)","copies spliced")]
yy=0.330
for tag,d in chk:
    txt(0.07,yy,"✓ "+tag,10,w='bold',c=GREEN); txt(0.17,yy,d,9.0,c='#333'); yy-=0.0235
txt(0.06,0.232,"identifiability: FULL     arrangement: segmental",9.3,mono=True,c='#222')
box(0.06,0.066,0.40,0.085,'#ffffff',ec=GREEN,lw=1.6)
txt(0.075,0.141,"VERDICT:  FAMILY  ✓",13,w='bold',c=GREEN)
txt(0.075,0.114,"2+ independently-expressed copies, each with its own full-",9.0,c='#333')
txt(0.075,0.097,"length sequence-decisive reads. A genuine multi-copy family",9.0,c='#333')
txt(0.075,0.080,"StringTie's primary-only assembly cannot reconstruct.",9.0,c='#333')

# ================= COUNTEREXAMPLE (right) =================
box(0.51,0.05,0.46,0.585,'#fdecea',ec=RED,lw=1.8)
txt(0.53,0.622,"COUNTEREXAMPLE  -  fails the definition",12.5,w='bold',c=RED)
txt(0.53,0.601,"SORD + predicted paralog LOC: 2 annotated loci, tandem.",9.8,c='#333')
xc,xd=0.58,0.93
model(xc,xd,0.572,BLUE,'#cfe2f0'); txt(0.575,0.572,"SORD",8.5,ha='right',va='center',w='bold',c=BLUE)
for i in range(4): flnc(xc,xd,0.553-0.012*i,BLUE)
txt(0.58,0.500,"full-length reads anchored to SORD (the genuinely expressed gene)",8,c=BLUE)
model(xc,xd,0.476,GREY,'#e6e6e6'); txt(0.575,0.476,"LOC",8.5,ha='right',va='center',w='bold',c=GREY)
# the SAME full-length reads multi-map onto LOC, but are spillover (fit SORD better)
for i in range(4):
    flnc(xc,xd,0.457-0.012*i,'#f3c9c2',hatch='///',ec='#c66')
    ax.add_patch(FancyArrowPatch((0.70+0.03*i,0.460),(0.70+0.03*i,0.548),
                 arrowstyle='-|>',mutation_scale=8,lw=0.8,color=RED,alpha=0.45,zorder=4))
txt(0.58,0.405,"same FLNC reads multi-map to LOC but fit SORD better over their full",8,c=RED)
txt(0.58,0.390,"length (dNM favors SORD)  ->  0 reads own LOC  (spillover / echo)",8,c=RED,w='bold')
txt(0.54,0.360,"Criteria:",10.5,w='bold')
chk2=[("✓","(M)","2 annotated loci  >= 2",GREEN),
      ("✓","(H)","connected by multi-mapping reads",GREEN),
      ("✗","(X)","only SORD expressed; LOC has 0 anchored reads",RED),
      ("","","=> 1 expressed copy, not >= 2",RED)]
yy=0.338
for mark,tag,d,col in chk2:
    txt(0.55,yy,(mark+" "+tag).strip(),10,w='bold',c=col); txt(0.66,yy,d,9.0,c='#333'); yy-=0.0235
txt(0.54,0.244,"naive secondary use would FABRICATE the LOC copy;",9.0,c=RED)
txt(0.54,0.228,"the anchor test (X) rejects it.",9.0,c=RED,w='bold')
box(0.54,0.066,0.40,0.085,'#ffffff',ec=RED,lw=1.6)
txt(0.555,0.141,"VERDICT:  spillover  ✗",13,w='bold',c=RED)
txt(0.555,0.114,"M and H pass, but X fails: a single expressed gene whose full-",9.0,c='#333')
txt(0.555,0.097,"length reads echo onto a homologous locus. NOT a multi-copy",9.0,c='#333')
txt(0.555,0.080,"family - emitting the second copy would be a false positive.",9.0,c='#333')

fig.savefig("family_definition_figure.png",dpi=150,bbox_inches='tight')
print("wrote family_definition_figure.png")
