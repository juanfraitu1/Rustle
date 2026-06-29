#!/usr/bin/env python3
"""Figure 1 — WHAT a minimizer is, and why discriminative minimizers predict
paralog identifiability. Didactic: real FNV-1a hashing, small k/w for legibility,
then the real RABL2-vs-DAZ anchor numbers. No claims beyond the algorithm + the
two measured anchors.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
import numpy as np

# ---- minimap2's EXACT minimizer hash (sketch.c: hash64 over the CANONICAL k-mer) ----
# small k/w below so the boxes are legible; minimap2 (splice:hq) uses k=15, w=5.
def hash64(key, mask):
    key = (~key + (key << 21)) & mask
    key = key ^ (key >> 24)
    key = ((key + (key << 3)) + (key << 8)) & mask
    key = key ^ (key >> 14)
    key = ((key + (key << 2)) + (key << 4)) & mask
    key = key ^ (key >> 28)
    key = (key + (key << 31)) & mask
    return key
_NT = {"A":0, "C":1, "G":2, "T":3}
def _canon(kmer):               # 2-bit pack; canonical = min(forward, reverse-complement)
    f = r = 0
    for ch in kmer: f = (f << 2) | _NT[ch]
    for ch in reversed(kmer): r = (r << 2) | (3 - _NT[ch])
    return min(f, r), (1 << (2 * len(kmer))) - 1
def mm_hash(kmer):
    v, mask = _canon(kmer); return hash64(v, mask)
def h2(s):                      # display form (real minimap2 hash, mod 100) used in panels
    return mm_hash(s) % 100

K_demo, W_demo = 5, 5           # teaching scale; minimap2 runs k=15, w=5
# A real ~40 bp DAZ1 slice + a one-base paralog variant: one substitution makes only the
# windows that SPAN it discriminative; the rest stay shared (this is why a near-identical
# copy has FEW discriminative minimizers — the visual point of the DAZ floor).
SEQA = "AGACAGTATCAGCAATAGGTAGTAGGAACATTTTTCTCAG"
SEQB = "AGACAGTATCAGCAATAGGTACTAGGAACATTTTTCTCAG"   # single substitution (G->C) mid-sequence

BLUE, GREEN, GREY, RED, INK = "#1f6fb2", "#2e9e5b", "#9aa3ab", "#c0392b", "#222"

def minimizers_with_pos(seq, k=K_demo, w=W_demo):
    """Return list of (window_start, minimizer_kmer, minimizer_pos, hashval)."""
    n = len(seq) - k + 1
    out = []
    for ws in range(0, max(n - w, 1)):
        we = min(ws + w, n)
        best = None
        for i in range(ws, we):
            hv = h2(seq[i:i+k])
            if best is None or hv < best[0]:
                best = (hv, i, seq[i:i+k])
        out.append((ws, best[2], best[1], best[0]))
    return out

fig = plt.figure(figsize=(13.5, 12.6))
gs = fig.add_gridspec(4, 1, height_ratios=[1.15, 1.0, 1.25, 0.85], hspace=0.5,
                      left=0.05, right=0.975, top=0.925, bottom=0.045)

fig.suptitle("What a minimizer is — and why it predicts whether two gene copies are tellable apart",
             fontsize=16, fontweight="bold", y=0.978)
fig.text(0.05, 0.948,
         "A minimizer is the smallest-HASH k-mer in each sliding window. minimap2 SEEDS on these, so they decide where a read can map.",
         fontsize=10.5, style="italic", color="#444")

# ================= Panel 0: what the hash is =================
ax0 = fig.add_subplot(gs[0]); ax0.axis("off")
ax0.set_xlim(0, 10); ax0.set_ylim(0, 10)
ax0.text(0, 9.4, "⓪ The hash (minimap2's own): turn a k-mer into a number — so “smallest” is well-defined, spread out, and identical in every copy",
         fontsize=12.3, fontweight="bold")
# the recipe + a fully worked example
ax0.text(0.0, 8.25, "Step 1  2 bits/base: A=00 C=01 G=10 T=11, glued into one binary number.   Step 2  canonical = min(it, reverse-complement).   Step 3  hash64.",
         fontsize=9.0, family="monospace", color=INK)
ax0.text(0.0, 7.35, "worked:   TATCA → 11 00 11 01 00 = 820     canonical = min(820, revcomp TGATA=908) = 820     hash64(820) = 223",
         fontsize=9.4, family="monospace", color=BLUE, fontweight="bold")
# the table of demo k-mers
demo_kmers = ["TATCA", "GTATC", "AGTAT"]
ax0.text(0.0, 6.3, "k-mer", fontsize=9.5, fontweight="bold", color="#555")
ax0.text(1.7, 6.3, "hash64  (k=5 here → 10-bit, 0–1023; real minimap2 k=15 → ~10⁹)", fontsize=9.5, fontweight="bold", color="#555")
for i, km in enumerate(demo_kmers):
    y = 5.4 - i * 0.95
    hv = mm_hash(km)
    ax0.text(0.0, y, km, fontsize=10.5, family="monospace", color=INK)
    ax0.text(2.0, y, f"{hv:>4}", fontsize=9.8, family="monospace", color=BLUE)
    ax0.text(3.6, y, f"→ {hv % 100:02d}  (mod 100, the form used below)", fontsize=8.5, color="#888")
ax0.text(0.0, 1.95,
         "(1) DETERMINISTIC — the same k-mer ALWAYS gives the same number, so where two copies share sequence they pick the SAME minimizer",
         fontsize=9.5, color="#444")
ax0.text(0.55, 1.18, "(a read landing there cannot tell the copies apart → it ties).   ", fontsize=9.5, color=RED)
ax0.text(3.7, 1.18, "(2) CANONICAL — the k-mer and its reverse complement hash the same, so strand never matters.",
         fontsize=9.5, color="#444")
ax0.text(0.0, 0.3,
         "(3) SPREAD-OUT — scatters k-mers evenly, so “smallest” is an arbitrary-but-stable pick.   This is the SAME hash minimap2 used to align these reads.",
         fontsize=9.5, color="#444")

# ================= Panel A: the sliding-window mechanic =================
axA = fig.add_subplot(gs[1]); axA.axis("off")
axA.set_xlim(0, len(SEQA)); axA.set_ylim(-3.2, 3.6)
axA.text(0, 3.3, "① One window → one minimizer", fontsize=12.5, fontweight="bold")
axA.text(0, 2.75, f"(shown at k={K_demo}, w={W_demo} for legibility; minimap2 uses k=15, w=5 — same rule)",
         fontsize=9.5, color="#555")

# draw the sequence letters
for i, c in enumerate(SEQA):
    axA.text(i + 0.5, 1.7, c, ha="center", va="center", fontsize=12.5, family="monospace", color=INK)

# show ONE window (window starting at index 3): its w k-mers and their hashes, circle the min
ws = 3
kmers = [(i, SEQA[i:i+K_demo], h2(SEQA[i:i+K_demo])) for i in range(ws, ws + W_demo)]
minh = min(k[2] for k in kmers)
# bracket over the window's span
x0, x1 = ws, ws + W_demo + K_demo - 1
axA.annotate("", xy=(x1, 2.25), xytext=(x0, 2.25),
             arrowprops=dict(arrowstyle="-", color="#888", lw=1.2))
axA.text((x0 + x1) / 2, 2.45, f"window = {W_demo} consecutive {K_demo}-mers", ha="center", fontsize=9.5, color="#555")
for r, (i, km, hv) in enumerate(kmers):
    yy = 0.7 - r * 0.85
    is_min = (hv == minh)
    fc = "#eaf6ef" if is_min else "#f4f5f6"
    ec = GREEN if is_min else "#ccc"
    axA.add_patch(Rectangle((i, yy - 0.32), K_demo, 0.64, fc=fc, ec=ec, lw=1.6 if is_min else 1.0, zorder=2))
    axA.text(i + K_demo/2, yy, km, ha="center", va="center", fontsize=11, family="monospace",
             color=INK, zorder=3)
    axA.text(i + K_demo + 0.25, yy, f"hash={hv:02d}", ha="left", va="center", fontsize=9.5,
             color=GREEN if is_min else "#888", fontweight="bold" if is_min else "normal")
    if is_min:
        axA.text(i + K_demo + 2.7, yy, "◄ smallest = MINIMIZER", ha="left", va="center",
                 fontsize=10, color=GREEN, fontweight="bold")

# ================= Panel B: shared vs discriminative across a SNP =================
axB = fig.add_subplot(gs[2]); axB.axis("off")
axB.set_xlim(0, len(SEQA)); axB.set_ylim(-1.2, 6.2)
axB.text(0, 5.9, "② Two near-identical gene copies — which minimizers separate them?",
         fontsize=12.5, fontweight="bold")

snp = next(i for i in range(len(SEQA)) if SEQA[i] != SEQB[i])  # the differing base
def draw_copy(ax, seq, y, label, color):
    ax.text(-0.3, y, label, ha="right", va="center", fontsize=11, fontweight="bold", color=color)
    for i, c in enumerate(seq):
        diff = (seq[i] != (SEQB[i] if seq is SEQA else SEQA[i]))
        ax.text(i + 0.5, y, c, ha="center", va="center", fontsize=12, family="monospace",
                color=RED if diff else INK, fontweight="bold" if diff else "normal")
        if diff:
            ax.add_patch(Rectangle((i, y - 0.35), 1, 0.7, fc="none", ec=RED, lw=1.4, zorder=1))

draw_copy(axB, SEQA, 4.7, "copy A", BLUE)
draw_copy(axB, SEQB, 3.9, "copy B", "#8e44ad")
axB.text(snp + 0.5, 5.35, "1 substitution", ha="center", fontsize=9, color=RED)
axB.annotate("", xy=(snp + 0.5, 5.05), xytext=(snp + 0.5, 5.25),
             arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))

# compute minimizers per window for both copies; tick each window shared/discriminative
mA = minimizers_with_pos(SEQA)
mB = minimizers_with_pos(SEQB)
n_disc = 0
for (wsA, kmA, posA, hA), (wsB, kmB, posB, hB) in zip(mA, mB):
    shared = (hA == hB and kmA == kmB)
    if not shared:
        n_disc += 1
    cx = wsA + (W_demo + K_demo - 1) / 2
    col = GREY if shared else GREEN
    axB.add_patch(Rectangle((wsA, 2.55), W_demo + K_demo - 1, 0.5,
                            fc=col, ec="white", alpha=0.85 if not shared else 0.4, lw=0.5))
axB.text(-0.3, 2.8, "windows", ha="right", va="center", fontsize=9.5, color="#555")

# legend / interpretation
axB.add_patch(Rectangle((0.3, 1.15), 0.9, 0.45, fc=GREY, alpha=0.4, ec="white"))
axB.text(1.4, 1.37, "SHARED minimizer — identical in both copies → a read here SEEDS on both → MAPQ 0, copy is a coin-flip (TIE)",
         va="center", fontsize=10, color="#555")
axB.add_patch(Rectangle((0.3, 0.45), 0.9, 0.45, fc=GREEN, alpha=0.85, ec="white"))
axB.text(1.4, 0.67, "DISCRIMINATIVE minimizer — window spans the substitution → unique to ONE copy → the read is PLACEABLE",
         va="center", fontsize=10, color=GREEN, fontweight="bold")
axB.text(0.3, -0.55,
         f"ONE substitution makes only {n_disc} of {len(mA)} windows discriminative — the rest are shared.   "
         "disc_frac = discriminative / all windows, summed over the copy:",
         fontsize=10.5, color=INK, fontweight="bold")
axB.text(0.3, -1.0,
         "a copy peppered with substitutions (RABL2) → high disc_frac → assignable;   a near-identical copy (DAZ core) → almost none → reads tie.",
         fontsize=10, color="#555", style="italic")

# ================= Panel C: the real anchor numbers =================
axC = fig.add_subplot(gs[3]); axC.axis("off")
axC.set_xlim(0, 10); axC.set_ylim(0, 4)
axC.add_patch(FancyBboxPatch((0.05, 0.15), 9.9, 3.7, boxstyle="round,pad=0.02,rounding_size=0.12",
                             fc="#f7f9fb", ec="#cbd5dd", lw=1.2))
axC.text(0.4, 3.45, "③ The score, computed live on this genome (GGO.fasta) — see the IGV figure for the read-level consequence",
         fontsize=12, fontweight="bold")

# numbers are the live pairwise disc_frac reported by fig_minimizer_igv.py (minimap2 minimizers)
rows = [
    ("RABL2A vs RABL2B", 0.30, "IDENTIFIABLE", GREEN,
     "tie-breakers along the gene → reads resolvable → genuine VG copy win"),
    ("DAZ3 vs DAZ1", 0.001, "SEQUENCE FLOOR", RED,
     "~0 tie-breakers (near-identical copies) → reads tie → emitting = fabrication risk"),
]
axC.text(0.4, 2.78, "copy pair", fontsize=10, fontweight="bold", color="#555")
axC.text(3.05, 2.78, "discriminative\nminimizers", fontsize=10, fontweight="bold", color="#555", ha="center")
axC.text(5.4, 2.78, "verdict", fontsize=10, fontweight="bold", color="#555", ha="center")
for r, (name, disc, verdict, col, why) in enumerate(rows):
    y = 2.05 - r * 0.95
    axC.text(0.4, y, name, fontsize=11, fontweight="bold", color=INK)
    axC.add_patch(Rectangle((2.55, y - 0.12), 1.0, 0.24, fc="#e7ebef", ec="#cbd5dd"))
    axC.add_patch(Rectangle((2.55, y - 0.12), 1.0 * max(disc, 0.004), 0.24, fc=col, ec="none"))
    axC.text(3.05, y + 0.27, f"{disc:.0%}" if disc >= 0.01 else "~0%", ha="center", fontsize=10, color=col, fontweight="bold")
    axC.text(5.4, y, verdict, ha="center", fontsize=10.5, color=col, fontweight="bold")
    axC.text(0.4, y - 0.34, why, fontsize=9, color="#666", style="italic")

axC.text(0.4, 0.32,
         "Predictor (NO reads needed): disc_frac ≥ 0.30 → IDENTIFIABLE   ·   0.12–0.30 → BORDERLINE   ·   < 0.12 → SEQUENCE FLOOR.  "
         "% identity is NOT the axis — discriminative-minimizer fraction is.",
         fontsize=9.5, color="#444")

out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/fig1_minimizer_explained.png"
fig.savefig(out, dpi=160, facecolor="white")
print("wrote", out)
