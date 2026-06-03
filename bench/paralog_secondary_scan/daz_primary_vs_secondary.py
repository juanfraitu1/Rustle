#!/usr/bin/env python3
"""DAZ two-locus figure: DAZ1 carries the PRIMARY alignments; DAZ3 is covered
almost entirely by SECONDARY alignments.

Pure factual claim from GGO.bam (gorilla IsoSeq). No expression/recovery claim.
Tracks:
  samtools depth -a -r CHR:SPAN -G 0x900          GGO.bam  -> primary  (excl. 0x100+0x800)
  samtools depth -a -r CHR:SPAN -g 0x100 -G 0x800 GGO.bam  -> secondary (req 0x100, excl 0x800)
Counts (samtools view -c):
  DAZ1: 200 primary / 20 secondary     DAZ3: 20 primary / 197 secondary
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CHR = "NC_073248.2"
DAZ1 = (42_783_133, 42_859_657)   # - strand
DAZ3 = (42_879_918, 42_945_552)   # + strand  (LOC129530216)
COUNTS = {"DAZ1": (200, 20), "DAZ3": (20, 197)}  # (primary, secondary)

PRIM = "#1f6fb2"   # blue
SEC  = "#c0392b"   # red
mb = 1e6

def load(path):
    a = np.loadtxt(path, usecols=(1, 2), dtype=np.int64)
    return a[:, 0], a[:, 1]

pp, pdp = load("/tmp/daz_depth_primary.txt")
sp, sdp = load("/tmp/daz_depth_secondary.txt")

def binned(pos, dep, step=200):
    n = (len(pos) // step) * step
    return pos[:n].reshape(-1, step).mean(1), dep[:n].reshape(-1, step).mean(1)

bpp, bpd = binned(pp, pdp)
bsp, bsd = binned(sp, sdp)

fig = plt.figure(figsize=(11, 7))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.25], hspace=0.32)

# ---------- Panel A: grouped read-count bars ----------
axA = fig.add_subplot(gs[0])
loci = ["DAZ1\n(−strand)", "DAZ3 / LOC129530216\n(+strand)"]
x = np.arange(2)
w = 0.38
prim_counts = [COUNTS["DAZ1"][0], COUNTS["DAZ3"][0]]
sec_counts  = [COUNTS["DAZ1"][1], COUNTS["DAZ3"][1]]
b1 = axA.bar(x - w/2, prim_counts, w, label="Primary alignments", color=PRIM)
b2 = axA.bar(x + w/2, sec_counts,  w, label="Secondary alignments", color=SEC)
axA.bar_label(b1, fontsize=11, fontweight="bold", padding=2)
axA.bar_label(b2, fontsize=11, fontweight="bold", padding=2)
axA.set_xticks(x); axA.set_xticklabels(loci, fontsize=11)
axA.set_ylabel("read alignments\nin locus", fontsize=11)
axA.set_ylim(0, 230)
axA.legend(fontsize=10, frameon=False, loc="upper center", ncol=2)
for s in ("top", "right"):
    axA.spines[s].set_visible(False)
axA.set_title("Read-alignment census per locus", fontsize=11.5, loc="left", color="0.3")

# ---------- Panel B: overlaid depth tracks ----------
axB = fig.add_subplot(gs[1])
axB.axvspan(DAZ1[0]/mb, DAZ1[1]/mb, color=PRIM, alpha=0.05)
axB.axvspan(DAZ3[0]/mb, DAZ3[1]/mb, color=SEC,  alpha=0.05)
axB.fill_between(bpp/mb, bpd, color=PRIM, alpha=0.80, linewidth=0, label="Primary depth")
axB.fill_between(bsp/mb, bsd, color=SEC,  alpha=0.55, linewidth=0, label="Secondary depth")
axB.set_xlabel("%s position (Mb)" % CHR, fontsize=11)
axB.set_ylabel("read depth", fontsize=11)
axB.set_ylim(0, 215)
axB.legend(fontsize=10, frameon=False, loc="upper right")
for s in ("top", "right"):
    axB.spines[s].set_visible(False)

c1 = (DAZ1[0]+DAZ1[1])/2/mb
c3 = (DAZ3[0]+DAZ3[1])/2/mb
axB.annotate("DAZ1", xy=(c1, 1), xycoords=("data", "axes fraction"),
             xytext=(0, 2), textcoords="offset points", ha="center", va="bottom",
             fontsize=11, fontweight="bold", color=PRIM)
axB.annotate("DAZ3", xy=(c3, 1), xycoords=("data", "axes fraction"),
             xytext=(0, 2), textcoords="offset points", ha="center", va="bottom",
             fontsize=11, fontweight="bold", color=SEC)
# call out the key fact: primary depth ~0 on DAZ3
axB.annotate("primary depth ≈ 0\n(invisible to a\nprimary-only assembler)",
             xy=(c3, 1.5), xytext=(c3, 70), ha="center", fontsize=9.5, color=PRIM,
             arrowprops=dict(arrowstyle="->", color=PRIM, lw=1.2))
axB.set_title("Per-base depth across the DAZ region", fontsize=11.5, loc="left", color="0.3")

fig.suptitle("DAZ: two paralog copies, opposite alignment classes  (gorilla IsoSeq, GGO.bam)",
             fontsize=13.5, fontweight="bold", y=0.97)
fig.text(0.5, 0.01,
         "DAZ1 carries the primary alignments; DAZ3 is covered almost entirely by secondary alignments "
         "(primary depth ≈0). A primary-only assembler never instantiates the DAZ3 locus.",
         ha="center", fontsize=9.5, style="italic", color="0.35")

fig.subplots_adjust(top=0.91, bottom=0.07)
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/daz_primary_vs_secondary.png"
fig.savefig(out, dpi=160)
print("wrote", out)
