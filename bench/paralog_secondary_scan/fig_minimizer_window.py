#!/usr/bin/env python3
"""Slide: how to READ a k=15, w=5 minimizer against the sequence.
A window = 5 CONSECUTIVE 15-mers (so it spans 15+4 = 19 bp); slide it one base at a
time; the minimizer is the smallest-hash 15-mer in the window. Real minimap2 hashes.
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch

K, W = 15, 5
SEQ = "CATTTTTCTCAGTTTAGAATTAA"          # real DAZ1 slice (23 bp → 9 fifteen-mers)
BLUE, GREEN, GREY, INK, AMB = "#1f6fb2", "#2e9e5b", "#aeb6bd", "#222", "#e08a1e"

# minimap2 canonical hash64
def hash64(key, mask):
    key=(~key+(key<<21))&mask; key=key^(key>>24)
    key=((key+(key<<3))+(key<<8))&mask; key=key^(key>>14)
    key=((key+(key<<2))+(key<<4))&mask; key=key^(key>>28)
    key=(key+(key<<31))&mask; return key
NT={"A":0,"C":1,"G":2,"T":3}
def mm_hash(kmer):
    f=r=0
    for ch in kmer: f=(f<<2)|NT[ch]
    for ch in reversed(kmer): r=(r<<2)|(3-NT[ch])
    return hash64(min(f,r),(1<<(2*len(kmer)))-1)

nk = len(SEQ)-K+1                                  # number of 15-mers
kmers = [(p, SEQ[p:p+K], mm_hash(SEQ[p:p+K]) % 100) for p in range(nk)]

fig, ax = plt.subplots(figsize=(13.6, 8.6)); ax.axis("off")
ax.set_xlim(-9, 26); ax.set_ylim(-1.5, nk+4)

ax.text(-9, nk+3.2, "How to read a  k = 15,  w = 5  minimizer", fontsize=15.5, fontweight="bold")
ax.text(-9, nk+2.4, "A window is 5 CONSECUTIVE 15-mers (it spans 15 + 4 = 19 bp). Slide it 1 base at a time; keep the smallest-hash 15-mer in each window.",
        fontsize=10.2, color="#444")

# ---- the sequence ruler ----
seq_y = nk+1.1
for i,c in enumerate(SEQ):
    ax.text(i+0.5, seq_y, c, ha="center", va="center", fontsize=12, family="monospace", color=INK)
ax.text(-0.5, seq_y, "seq", ha="right", va="center", fontsize=9, color="#555")
for i in range(0, len(SEQ)+1, 5):
    ax.text(i, seq_y+0.55, str(i), ha="center", fontsize=7, color="#aaa")

# ---- cascade of 15-mers (one per row), each a bar + hash ----
def row_y(r): return nk-1-r           # row 0 near the top, just below the sequence
for p, km, h in kmers:
    y = row_y(p)
    ax.add_patch(Rectangle((p, y-0.32), K, 0.64, fc="#eef2f5", ec="#c7ced4", lw=1.0, zorder=2))
    ax.text(p+K/2, y, f"15-mer @{p}", ha="center", va="center", fontsize=8.0, family="monospace", color="#667")
    ax.text(p+K+0.3, y, f"hash {h:02d}", ha="left", va="center", fontsize=8.6, color=INK)

# ---- three rolling-window positions: draw brackets + dotted line to the chosen min ----
windows = [0,1,2]
wx = [-2.0, -4.3, -6.6]            # staggered left brackets so they don't overlap
wcol = [BLUE, AMB, GREEN]
picked = {}                        # min-kmer position -> (windows that picked it, color)
for s,xx,col in zip(windows,wx,wcol):
    rows=list(range(s,s+W))
    ytop,ybot=row_y(s)+0.45,row_y(s+W-1)-0.45
    ax.add_patch(FancyBboxPatch((xx-0.5,ybot),0.55,ytop-ybot,
                 boxstyle="round,pad=0.02,rounding_size=0.15",fc="none",ec=col,lw=2.0,zorder=4))
    ax.text(xx-0.25,(ytop+ybot)/2,f"window\n@{s}",ha="center",va="center",fontsize=8.2,
            color=col,fontweight="bold",rotation=90)
    mp=min((kmers[r] for r in rows),key=lambda x:x[2])[0]
    ax.add_line(plt.Line2D([xx,mp],[(ytop+ybot)/2,row_y(mp)],color=col,lw=1.0,ls=":",zorder=3))
    picked.setdefault(mp,[[],col])[0].append(s)

# highlight each distinct minimizer once, label which window(s) chose it
for mp,(ws,col) in picked.items():
    y=row_y(mp); h=kmers[mp][2]
    ax.add_patch(Rectangle((mp,y-0.36),K,0.72,fc="none",ec=col,lw=2.6,zorder=5))
    wtxt = f"window @{ws[0]}" if len(ws)==1 else "windows @"+" & @".join(map(str,ws))
    ax.text(mp+K+3.0,y,f"◄ smallest in {wtxt}  →  MINIMIZER (15-mer @{mp}, hash {h:02d})",
            ha="left",va="center",fontsize=8.6,color=col,fontweight="bold")

ax.text(-9,-0.4,
        "Window @0 picks 15-mer @1.  Windows @1 and @2 both pick 15-mer @5 — once it enters, its hash (09) is smallest until it drops out, so @5 is emitted ONCE.",
        fontsize=9.4,color=INK,fontweight="bold")
ax.text(-9,-1.05,
        "That persistence is the point: a 15-mer wins every window it sits in until a smaller-hash one enters — so minimizers are SPARSE (≈1 per w bases), not one per base.",
        fontsize=9.3,color="#555")

out="/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/fig1c_window_reading.png"
fig.savefig(out, dpi=160, facecolor="white", bbox_inches="tight"); print("wrote", out, "| picked:", {p:w for p,(w,_) in picked.items()})
