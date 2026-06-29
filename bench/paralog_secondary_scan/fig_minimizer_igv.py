#!/usr/bin/env python3
"""IGV evidence for the minimizer predictor, on REAL reads (GGO.bam), ONE FAMILY PER FIGURE.

  fig2a_rabl2_igv.png : RABL2A + RABL2B  — the IDENTIFIABLE family
  fig2b_daz_igv.png   : DAZ1   + DAZ3    — the SEQUENCE-FLOOR family

Minimizers use minimap2's OWN scheme (mm_sketch: canonical k-mer + invertible hash64,
k=15, w=5 — the splice:hq seeding that aligned these very reads). Canonical k-mers are
strand-symmetric, so NO manual reverse-complement is needed. Each panel: real IsoSeq
pileup over one copy (blue=primary, red=secondary); track below = minimap2 minimizers vs
the OTHER copy of the family (GREEN=discriminative/copy-unique, GREY=shared → reads tie).

Computed live from GGO.fasta / GGO.bam.
"""
import numpy as np, subprocess, pysam
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

BAM    = "/mnt/c/Users/jfris/Desktop/GGO.bam"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
PRIM, SEC, GREEN, GREY, INK = "#1f6fb2", "#c0392b", "#2e9e5b", "#c4cbd1", "#222"
K, W = 15, 5     # minimap2 splice:hq defaults

RABL2A = ("NC_073235.2", 15_131_653, 15_147_533, "RABL2A")
RABL2B = ("NC_086018.1", 48_818_440, 48_831_973, "RABL2B")
DAZ1   = ("NC_073248.2", 42_783_133, 42_859_657, "DAZ1")
DAZ3   = ("NC_073248.2", 42_879_918, 42_945_552, "DAZ3")

# ---- minimap2's exact minimizer scheme (sketch.c) ----
def hash64(key, mask):
    key = (~key + (key << 21)) & mask
    key = key ^ (key >> 24)
    key = ((key + (key << 3)) + (key << 8)) & mask
    key = key ^ (key >> 14)
    key = ((key + (key << 2)) + (key << 4)) & mask
    key = key ^ (key >> 28)
    key = (key + (key << 31)) & mask
    return key
NT4 = {65:0,67:1,71:2,84:3, 97:0,99:1,103:2,116:3}
def mm_minpos(seq, k=K, w=W):
    """minimap2 minimizers: list of (kmer_start_offset, hash), canonical + hash64."""
    mask=(1<<(2*k))-1; shift1=2*(k-1)
    kmer=[0,0]; l=0; stream=[]   # (hash, offset) per position; hash None if invalid
    for i,ch in enumerate(seq):
        c=NT4.get(ch,4)
        if c<4:
            kmer[0]=((kmer[0]<<2)|c)&mask
            kmer[1]=(kmer[1]>>2)|((3-c)<<shift1)
            l+=1
            if l>=k and kmer[0]!=kmer[1]:
                z=0 if kmer[0]<kmer[1] else 1
                stream.append((hash64(kmer[z],mask), i-k+1))
            else:
                stream.append((None, i-k+1))
        else:
            l=0; stream.append((None, i))
    chosen={}
    for s in range(0, max(len(stream)-w+1, 1)):
        win=[x for x in stream[s:s+w] if x[0] is not None]
        if win:
            h,off=min(win, key=lambda x:x[0])
            chosen[off]=h
    return [(off,h) for off,h in chosen.items()]

def fetch(c,a,b):
    s=subprocess.run(["samtools","faidx",GENOME,f"{c}:{a+1}-{b}"],capture_output=True,text=True).stdout
    return "".join(l for l in s.splitlines() if not l.startswith(">")).encode()
def reads_in(c,a,b):
    out=[]
    with pysam.AlignmentFile(BAM,"rb") as bam:
        for r in bam.fetch(c,a,b):
            if r.is_unmapped or r.is_supplementary: continue
            blk=r.get_blocks()
            if blk: out.append((r.reference_start,r.reference_end,blk,r.is_secondary))
    return out
def pack(reads,pad=300):
    reads=sorted(reads,key=lambda r:r[0]); row_end=[]; placed=[]
    for s,e,blk,sec in reads:
        row=None
        for i,en in enumerate(row_end):
            if s>en+pad: row=i; row_end[i]=e; break
        if row is None: row=len(row_end); row_end.append(e)
        placed.append((row,s,e,blk,sec))
    return placed,len(row_end)
def exon_blocks(reads):
    """Merged read-covered intervals = the EXONS the IsoSeq reads actually sample."""
    cov=set()
    for s,e,blk,sec in reads:
        for bs,be in blk: cov.add((bs,be))
    merged=[]
    for s,e in sorted(cov):
        if merged and s<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],e))
        else: merged.append((s,e))
    return merged
def exonic_minpos(c, blocks):
    """minimap2 minimizers computed ONLY inside exon blocks (genomic positions)."""
    out=[]
    for s,e in blocks:
        if e-s>=15:
            for off,h in mm_minpos(fetch(c,s,e)): out.append((s+off,h))
    return out

def draw_panel(ax, copy, sib_set):
    c,a,b,label=copy
    reads=reads_in(c,a,b)
    blocks=exon_blocks(reads)
    own=exonic_minpos(c,blocks)                       # minimizers ON EXONS ONLY
    disc=sum(1 for _,h in own if h not in sib_set)/len(own) if own else 0
    placed,nrow=pack(reads)
    n_sec=sum(1 for _,_,_,_,s in placed if s); n_pri=len(placed)-n_sec
    track_lo,track_hi=0.05,1.15
    for row,s,e,blk,sec in placed:
        y=track_hi+0.55+row*0.92; col=SEC if sec else PRIM
        ax.add_line(Line2D([s,e],[y,y],color=col,lw=0.7,alpha=0.5,zorder=2))
        for bs,be in blk:
            ax.add_patch(Rectangle((bs,y-0.32),be-bs,0.64,fc=col,ec="none",alpha=0.9,zorder=3))
    # track band; shade ONLY the exon blocks so introns are visibly empty
    ax.add_patch(Rectangle((a,track_lo),b-a,track_hi-track_lo,fc="white",ec="#d7dde2",lw=0.8,zorder=1))
    for s,e in blocks:
        ax.add_patch(Rectangle((s,track_lo),e-s,track_hi-track_lo,fc="#eef1f3",ec="none",zorder=1))
    for off,h in own:
        d=h not in sib_set
        ax.add_line(Line2D([off,off],[track_lo+0.07,track_hi-0.07],
                    color=GREEN if d else GREY,lw=0.7,alpha=0.95 if d else 0.55,zorder=2))
    ax.text(a-(b-a)*0.012,(track_lo+track_hi)/2,"exon\nminimizers",ha="right",va="center",fontsize=8.2,color="#555")
    ax.set_xlim(a-(b-a)*0.06,b+(b-a)*0.02)
    ax.set_ylim(-0.25,track_hi+0.55+max(nrow,1)*0.92+0.5)
    ax.set_yticks([])
    xt=np.linspace(a,b,5); ax.set_xticks(xt); ax.set_xticklabels([f"{x/1e6:.2f}" for x in xt],fontsize=8.5)
    ax.set_xlabel(f"{c}  (Mb)",fontsize=9)
    for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
    ax.set_title(f"{label}   ·   {n_pri} primary / {n_sec} secondary reads   ·   discriminative EXON minimizers = {disc:.0%}",
                 fontsize=11.5,fontweight="bold",loc="left",color=INK,pad=6)
    return disc

def make_family(top, bottom, fam_title, verdict, accent, outfile):
    fig,axes=plt.subplots(2,1,figsize=(13.6,9.6),
                          gridspec_kw=dict(hspace=0.36,left=0.075,right=0.975,top=0.85,bottom=0.10))
    fig.suptitle(fam_title, fontsize=15.5, fontweight="bold", y=0.965)
    fig.text(0.075,0.905,
             "Real IsoSeq reads (GGO.bam). Minimizers are computed ONLY on EXONS (shaded) — the spliced reads carry no intron sequence, so introns are excluded.",
             fontsize=10,color="#444")
    setT=set(h for _,h in exonic_minpos(top[0], exon_blocks(reads_in(*top[:3]))))
    setB=set(h for _,h in exonic_minpos(bottom[0], exon_blocks(reads_in(*bottom[:3]))))
    dT=draw_panel(axes[0], top,    setB)   # top's sibling = bottom (exon minimizers)
    dB=draw_panel(axes[1], bottom, setT)   # bottom's sibling = top
    leg=[Line2D([0],[0],color=PRIM,lw=6,label="primary read"),
         Line2D([0],[0],color=SEC,lw=6,label="secondary read (multimaps)"),
         Line2D([0],[0],color=GREEN,lw=2.6,label="discriminative minimizer (tie-breaker present)"),
         Line2D([0],[0],color=GREY,lw=2.6,label="shared minimizer (no tie-breaker)")]
    fig.legend(handles=leg,loc="lower center",ncol=4,fontsize=9.2,frameon=False,bbox_to_anchor=(0.5,0.005))
    fig.text(0.5,0.05, verdict, ha="center", fontsize=11.5, fontweight="bold", color=accent)
    fig.savefig(outfile,dpi=160,facecolor="white")
    print(f"wrote {outfile} | {top[3]} disc={dT:.3f}  {bottom[3]} disc={dB:.3f}")

make_family(RABL2A, RABL2B,
            "RABL2 family — both copies carry tie-breakers → the multimapping is RESOLVABLE",
            "→ IDENTIFIABLE: VG recovers both copies (real win)", GREEN,
            "/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/fig2a_rabl2_igv.png")
make_family(DAZ1, DAZ3,
            "DAZ family — the copies are near-identical → almost no tie-breaker",
            "→ SEQUENCE FLOOR: reads tie; emitting a copy here would be fabrication", SEC,
            "/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/fig2b_daz_igv.png")
