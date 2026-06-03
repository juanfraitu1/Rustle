#!/usr/bin/env python3
"""Advisor deep-dive on pane 3 of summary_daz_tied.png.

Shows, IGV-style with REAL reads from /tmp/daz_aln.sam, what actually separates
a DAZ multimapping read between the two copies. Honest mechanism: the separator
is the read's *sequence cost* against each copy reference (minimap2 AS / de / NM),
NOT a coordinate overlap with the 10.8 kb block. Three traceable exemplars:
  - 65405150 : carries DAZ1-only sequence -> DAZ3 placement pays 626 mismatches +
               a 609 bp insertion -> RESOLVED to DAZ1.
  - 50660416 : coordinate-overlaps the block, but its in-block bases re-align to
               DAZ3 backbone equally well (NM 4 vs 3) -> minimap2 PARKED -> TIED.
  - 24709554 : confined to shared backbone, NM 0 vs 0 -> non-identifiable -> apportion.
Data: /tmp/daz_igv.json (built by /tmp/daz_igv_data.py).
"""
import json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon

J=json.load(open('/tmp/daz_igv.json'))
recs=J['recs']; DIAG=J['DIAG']; BLK=tuple(J['BLK'])
DAZ1=tuple(J['DAZ1']); DAZ3=tuple(J['DAZ3'])
D1OFF=DAZ1[0]; D3OFF=DAZ3[0]
BLUE,GREEN,RED,PURPLE,ORANGE,GREY='#2c7fb8','#1a9850','#d73027','#7b3294','#f16913','#666'
SPAN=77000.0  # copy-relative width to normalize both lanes (DAZ1 is ~76.5kb)

fig=plt.figure(figsize=(15,16)); ax=fig.add_axes([0,0,1,1]); ax.axis('off')
ax.set_xlim(0,1); ax.set_ylim(0,1)
def T(x,y,s,sz=11,ha='left',va='top',c='black',w='normal',mono=False,st='normal',rot=0):
    ax.text(x,y,s,size=sz,ha=ha,va=va,color=c,weight=w,style=st,rotation=rot,
            family='monospace' if mono else 'sans-serif',zorder=6)
def panel(x,y,w,h,fc,ec='#bbb'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.004",fc=fc,ec=ec,lw=1.2,zorder=0.5))

# coordinate helpers: copy-relative genomic -> figure-x within a lane [x0,x1]
def lane_x(g,off,x0,x1): return x0+(g-off)/SPAN*(x1-x0)

T(0.5,0.986,"Pane 3, decomposed: what really separates a DAZ multimapping read (real reads, IGV view)",
  17,ha='center',w='bold')
T(0.5,0.967,"Separator = the read's sequence cost against each copy reference (minimap2 AS / de / NM). "
  "Coordinate-overlap with the 10.8 kb block is only a proxy.",12,ha='center',c=GREY,st='italic')

# ====================== SECTION A: reference / SNV provenance ======================
ax0,ax1=0.07,0.93
panel(0.03,0.793,0.94,0.158,'#f6f8fb')
T(0.05,0.945,"A.  Where the 'diagnostic SNVs' live: they are GENOME differences between the two assembled copies",12.5,w='bold')
T(0.05,0.928,"Align the DAZ1 reference sequence to the DAZ3 reference (copy-vs-copy). The 7 mismatches are fixed positions where "
  "DAZ1's genome base != DAZ3's.",10,c='#333')
T(0.05,0.913,"They are NOT base-calls from individual reads -> they cannot be sequencing errors. A read merely *reads out* which "
  "copy's allele it carries at those 7 spots.",10,c='#333')
# DAZ1 ref bar
yb=0.862
T(0.055,yb+0.008,"DAZ1 ref",9,ha='right',va='center',w='bold',c=BLUE)
ax.add_patch(Rectangle((ax0,yb),ax1-ax0,0.016,fc='#cfe2f0',ec=BLUE,lw=1.1))
# block
bx0=lane_x(BLK[0],D1OFF,ax0,ax1); bx1=lane_x(BLK[1],D1OFF,ax0,ax1)
ax.add_patch(Rectangle((bx0,yb-0.002),bx1-bx0,0.020,fc=ORANGE,ec='#a33',lw=1.0))
T((bx0+bx1)/2,yb+0.030,"10.8 kb block\n(in DAZ1, absent in DAZ3)",8,ha='center',va='bottom',c=ORANGE,w='bold')
for d in DIAG:
    dx=lane_x(d,D1OFF,ax0,ax1); ax.plot([dx,dx],[yb,yb+0.016],c=RED,lw=1.4,zorder=5)
T(ax1+0.004,yb+0.008,"7 SNVs",8.5,va='center',c=RED,w='bold')
# DAZ3 ref bar (block absent -> shown as gap)
yc=0.822
T(0.055,yc+0.008,"DAZ3 ref",9,ha='right',va='center',w='bold',c=GREEN)
ax.add_patch(Rectangle((ax0,yc),bx0-ax0,0.016,fc='#cde6d7',ec=GREEN,lw=1.1))
ax.add_patch(Rectangle((bx1,yc),ax1-bx1,0.016,fc='#cde6d7',ec=GREEN,lw=1.1))
T((bx0+bx1)/2,yc+0.008,"block absent",7.5,ha='center',va='center',c=GREY,st='italic')
for d in DIAG:
    dx=lane_x(d,D1OFF,ax0,ax1)
    if not (bx0<dx<bx1): ax.plot([dx,dx],[yc,yc+0.016],c=RED,lw=1.4,zorder=5)
T(0.05,0.806,"=> Only 7 SNVs in 76 kb (~1 every 11 kb): too sparse to label most reads. The bulk signal is the dense per-base "
  "match/mismatch the aligner scores everywhere.",9.5,c='#333',st='italic')

# ====================== SECTION B: 3 IGV reads ======================
T(0.05,0.778,"B.  Three real reads, each shown at BOTH placements (what IGV draws: primary + secondary alignment of the same read)",12.5,w='bold')
T(0.05,0.766,"The two lanes per read are two SEPARATE loci at their own genomic coordinates (not homology-aligned). Each lane spans its copy.",
  9,c=GREY,st='italic')

def draw_lane(x0,x1,ylane,rec,off,base_c,label,show_feats):
    """Draw one IGV alignment row for `rec` on copy-relative axis."""
    h=0.012
    # backbone line across span (faint)
    ax.plot([x0,x1],[ylane+h/2,ylane+h/2],c='#e8e8e8',lw=0.8,zorder=1)
    # exons + introns
    ex=rec['exons']; introns=rec['introns']
    # intron lines
    for s,e in introns:
        xs=lane_x(s,off,x0,x1); xe=lane_x(e,off,x0,x1)
        ax.plot([xs,xe],[ylane+h/2,ylane+h/2],c=base_c,lw=0.7,zorder=2)
    for s,e in ex:
        xs=lane_x(s,off,x0,x1); xe=lane_x(e,off,x0,x1)
        ax.add_patch(Rectangle((xs,ylane),max(xe-xs,0.0008),h,fc=base_c,ec='none',zorder=3))
    # mismatches (only where SEQ present)
    for m in rec['mm']:
        mx=lane_x(m,off,x0,x1); ax.plot([mx,mx],[ylane-0.002,ylane+h+0.002],c=RED,lw=1.2,zorder=5)
    # insertions (big = DAZ1-only content forced into the wrong copy)
    for ipos,ilen in rec['ins']:
        if ilen>=30:
            ix=lane_x(ipos,off,x0,x1)
            ax.add_patch(Polygon([[ix-0.006,ylane+h+0.004],[ix+0.006,ylane+h+0.004],[ix,ylane+h-0.002]],
                         closed=True,fc=PURPLE,ec='none',zorder=6))
            T(ix,ylane+h+0.007,f"{ilen} bp\ninsertion",7,ha='center',va='bottom',c=PURPLE,w='bold')
    # SNV coverage markers
    if show_feats:
        for d in rec['snv_cov']:
            dx=lane_x(d,off,x0,x1)
            ax.plot(dx,ylane+h/2,marker='D',ms=5,c=RED,zorder=7,mec='white',mew=0.5)
    T(x0-0.004,ylane+h/2,label,8.5,ha='right',va='center',w='bold',c=base_c)
    return h

def read_group(ytop,key,title,verdict,vc,note):
    panel(0.03,ytop-0.205,0.94,0.205,'#ffffff')
    T(0.05,ytop-0.012,title,11.5,w='bold')
    r1=recs[key].get('DAZ1'); r3=recs[key].get('DAZ3')
    lx0,lx1=0.16,0.74
    # DAZ1 lane
    T(0.05,ytop-0.034,"DAZ1 locus",9,w='bold',c=BLUE)
    draw_lane(lx0,lx1,ytop-0.058,r1,D1OFF,BLUE,"DAZ1",True)
    # block shading on DAZ1 lane
    bx0=lane_x(BLK[0],D1OFF,lx0,lx1); bx1=lane_x(BLK[1],D1OFF,lx0,lx1)
    ax.add_patch(Rectangle((bx0,ytop-0.064),bx1-bx0,0.024,fc=ORANGE,ec='none',alpha=0.18,zorder=1.2))
    # SNV gridlines on DAZ1 lane
    for d in DIAG:
        dx=lane_x(d,D1OFF,lx0,lx1); ax.plot([dx,dx],[ytop-0.066,ytop-0.040],c=RED,lw=0.5,alpha=0.35,zorder=1.1)
    T(0.755,ytop-0.052,f"AS {r1['AS']}   NM {r1['NM']}   de {r1['de']}",9.5,va='center',mono=True,c='#222',w='bold')
    T(lx0,ytop-0.072,f"{r1['exons'][0][0]:,}",7,ha='center',c=GREY)
    T(lx1,ytop-0.072,f"..{r1['exons'][-1][1]:,}",7,ha='center',c=GREY)
    # DAZ3 lane
    T(0.05,ytop-0.094,"DAZ3 locus",9,w='bold',c=GREEN)
    draw_lane(lx0,lx1,ytop-0.118,r3,D3OFF,GREEN,"DAZ3",True)
    T(0.755,ytop-0.112,f"AS {r3['AS']}   NM {r3['NM']}   de {r3['de']}",9.5,va='center',mono=True,c='#222',w='bold')
    T(lx0,ytop-0.132,f"{r3['exons'][0][0]:,}",7,ha='center',c=GREY)
    T(lx1,ytop-0.132,f"..{r3['exons'][-1][1]:,}",7,ha='center',c=GREY)
    # verdict box
    ax.add_patch(FancyBboxPatch((0.05,ytop-0.198),0.88,0.052,boxstyle="round,pad=0.004",
                 fc='#fbfbf0',ec=vc,lw=1.4,zorder=2))
    T(0.065,ytop-0.157,"verdict:",10,w='bold',c=vc)
    T(0.145,ytop-0.157,verdict,10,w='bold',c=vc)
    T(0.065,ytop-0.180,note,9.0,c='#222')

read_group(0.760,'resolved',
  "Read 1  m64076_221110_210557/65405150/ccs  -- carries real DAZ1-only sequence",
  "RESOLVED -> DAZ1",GREEN,
  "DAZ1 fits cleanly (NM 5). Forced onto DAZ3 the same bases cost a 609 bp insertion + 626 mismatches (NM 626). "
  "It also reads the DAZ1 allele at SNV 42,840,676 (red diamond). dAS = 562.")
read_group(0.534,'parked',
  "Read 2  .../50660416/ccs  -- coordinate-overlaps the block but is NOT copy-specific",
  "TIED -> apportion (block overlap is a minimap2 'parking' artifact)",RED,
  "Its 'in-block' bases (52+92+147 bp) re-align to shared DAZ3 backbone with NM 3 ~ NM 4. Equal fit on both copies, "
  "covers no SNV. Overlapping the block's COORDINATES is not the same as carrying DAZ1-only SEQUENCE. dAS = 5.")
read_group(0.308,'tied',
  "Read 3  .../24709554/ccs  -- confined to the 99.97%-identical shared backbone",
  "TIED -> apportion (non-identifiable from sequence)",RED,
  "Both placements are a perfect tie: AS 2518=2518, NM 0=0, de 0=0, no block, no SNV. No evidence exists to assign it; "
  "hard-calling it inflates one copy. Correct handling = fractional apportionment. dAS = 0.")

# legend
ly=0.092
panel(0.03,0.012,0.94,0.083,'#eef2f8')
T(0.05,0.088,"Legend / how to read each row",10.5,w='bold')
def chip(x,y,c,shape='rect'):
    if shape=='rect': ax.add_patch(Rectangle((x,y),0.018,0.011,fc=c,ec='none'))
    elif shape=='tick': ax.plot([x+0.009,x+0.009],[y-0.002,y+0.013],c=c,lw=1.4)
    elif shape=='tri': ax.add_patch(Polygon([[x,y+0.013],[x+0.018,y+0.013],[x+0.009,y]],closed=True,fc=c))
    elif shape=='dia': ax.plot(x+0.009,y+0.0055,marker='D',ms=5,c=c,mec='white',mew=0.5)
items=[(BLUE,'rect','exon block (aligned, M)'),(RED,'tick','per-base mismatch (where SEQ stored)'),
       (PURPLE,'tri','>=30 bp insertion = DAZ1-only bases forced into DAZ3'),
       (RED,'dia','read covers one of the 7 diagnostic SNVs'),
       (ORANGE,'rect','10.8 kb DAZ1-only block (shaded)')]
xx=0.05
for c,sh,lab in items:
    chip(xx,0.058,c,sh); T(xx+0.024,0.063,lab,8.3,va='center'); xx+=0.185
T(0.05,0.030,"Manual trace: grep the read name in /tmp/daz_aln.sam (or load GGO.bam at NC_073248.2:42,783,133-42,945,552 in IGV). "
  "Compare the two records' NM/de/AS. Large gap + an insertion on one copy = real copy-specific sequence; equal NM on both = tied.",
  8.6,c='#333',st='italic')

fig.savefig('daz_pane3_igv.png',dpi=145,bbox_inches='tight')
print("wrote daz_pane3_igv.png")
