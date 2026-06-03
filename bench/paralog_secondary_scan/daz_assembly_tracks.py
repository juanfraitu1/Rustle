#!/usr/bin/env python3
"""IGV-style transcript-model tracks for the DAZ region.

Compares assembly output:
  - StringTie -L            (primary-only)  -> /tmp/st_daz.gtf
  - rustle (no --vg)        (primary-only)  -> /tmp/ru_prim_daz.gtf
  - rustle --vg             (secondary-aware)-> /tmp/ru_vg_daz.gtf

Each transcript is drawn IGV-style: thick boxes = exons, thin line = introns.
Factual claim: primary-only assemblers emit transcripts only at DAZ1; including
secondary alignments yields transcript models at the DAZ3 copy as well.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import defaultdict

CHR = "NC_073248.2"
DAZ1 = (42_783_133, 42_859_657)
DAZ3 = (42_879_918, 42_945_552)
SPAN = (42_780_000, 42_948_000)
mb = 1e6
PRIM = "#1f6fb2"
SEC  = "#c0392b"
GREY = "#444444"

def parse_gtf(path):
    """Return {tx_id: (strand, [(exon_start,exon_end),...])} for transcripts on CHR in SPAN."""
    tx = defaultdict(list)
    strand = {}
    try:
        fh = open(path)
    except FileNotFoundError:
        return {}
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != CHR or f[2] != "exon":
            continue
        s, e = int(f[3]), int(f[4])
        if e < SPAN[0] or s > SPAN[1]:
            continue
        attrs = f[8]
        tid = None
        for kv in attrs.split(";"):
            kv = kv.strip()
            if kv.startswith("transcript_id"):
                tid = kv.split('"')[1]
                break
        if tid is None:
            continue
        tx[tid].append((s, e))
        strand[tid] = f[6]
    fh.close()
    return {t: (strand[t], sorted(ex)) for t, ex in tx.items()}

def locus_of(exons):
    mid = (exons[0][0] + exons[-1][1]) / 2
    if DAZ1[0] <= mid <= DAZ1[1]:
        return "DAZ1"
    if DAZ3[0] <= mid <= DAZ3[1]:
        return "DAZ3"
    return "other"

TRACKS = [
    ("StringTie  (primary-only)",      "/tmp/st_daz.gtf",      PRIM),
    ("rustle, no --vg  (primary-only)", "/tmp/ru_prim_daz.gtf", PRIM),
    ("rustle --vg  (includes secondary alignments)", "/tmp/ru_vg_daz.gtf", SEC),
]

parsed = [(name, parse_gtf(path), col) for name, path, col in TRACKS]
for name, txs, _ in parsed:
    by = defaultdict(int)
    for _, (_, ex) in txs.items():
        by[locus_of(ex)] += 1
    print("%-45s DAZ1=%d DAZ3=%d other=%d" % (name, by["DAZ1"], by["DAZ3"], by["other"]))

# ----- plot -----
n = len(parsed)
fig, axes = plt.subplots(n, 1, figsize=(13, 2.0 * n + 1.4), sharex=True)
if n == 1:
    axes = [axes]

def draw_tx(ax, y, exons, color):
    x0 = exons[0][0] / mb
    x1 = exons[-1][1] / mb
    ax.plot([x0, x1], [y, y], color=color, lw=0.9, zorder=1)   # intron line
    for s, e in exons:
        ax.add_patch(Rectangle((s/mb, y - 0.32), (e - s)/mb, 0.64,
                               facecolor=color, edgecolor="none", zorder=2))

for ax, (name, txs, col) in zip(axes, parsed):
    ax.axvspan(DAZ1[0]/mb, DAZ1[1]/mb, color=PRIM, alpha=0.05, zorder=0)
    ax.axvspan(DAZ3[0]/mb, DAZ3[1]/mb, color=SEC,  alpha=0.05, zorder=0)
    # order transcripts by start, stack them
    items = sorted(txs.values(), key=lambda v: v[1][0][0])
    for i, (strand, exons) in enumerate(items):
        draw_tx(ax, i, exons, col)
    nd1 = sum(1 for _, (_, ex) in txs.items() if locus_of(ex) == "DAZ1")
    nd3 = sum(1 for _, (_, ex) in txs.items() if locus_of(ex) == "DAZ3")
    ax.set_ylim(-1, max(3, len(items)))
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_ylabel(name, fontsize=10.5, rotation=0, ha="right", va="center",
                  color=col, fontweight="bold")
    # per-locus count badges
    ax.text(0.985, 0.9, "DAZ1: %d tx    DAZ3: %d tx" % (nd1, nd3),
            transform=ax.transAxes, ha="right", va="top", fontsize=11,
            fontweight="bold",
            color=(GREY if nd3 == 0 else SEC),
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=col, alpha=0.9))

# locus headers on top axis
c1 = (DAZ1[0]+DAZ1[1])/2/mb
c3 = (DAZ3[0]+DAZ3[1])/2/mb
axes[0].annotate("DAZ1  (−strand)", xy=(c1, 1), xycoords=("data", "axes fraction"),
                 xytext=(0, 6), textcoords="offset points", ha="center", va="bottom",
                 fontsize=12, fontweight="bold", color=PRIM)
axes[0].annotate("DAZ3 / LOC129530216  (+strand)", xy=(c3, 1), xycoords=("data", "axes fraction"),
                 xytext=(0, 6), textcoords="offset points", ha="center", va="bottom",
                 fontsize=12, fontweight="bold", color=SEC)

axes[-1].set_xlim(SPAN[0]/mb, SPAN[1]/mb)
axes[-1].set_xlabel("%s position (Mb)" % CHR, fontsize=12)

fig.suptitle("Assembled transcripts in the DAZ region — primary-only vs secondary-aware",
             fontsize=14.5, fontweight="bold", y=0.985)
fig.text(0.5, 0.008,
         "Thick = exon, thin = intron. Primary-only tools already assemble DAZ3 (1 truncated 16-exon model). "
         "rustle --vg yields 5 longer DAZ3 isoforms but redistributes weight off DAZ1 (10→2 tx); the added DAZ3 "
         "evidence is cross-mapped DAZ1 reads, not validated DAZ3 expression.",
         ha="center", fontsize=9.0, style="italic", color="0.35")

fig.subplots_adjust(left=0.27, right=0.985, top=0.92, bottom=0.10, hspace=0.25)
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/daz_assembly_tracks.png"
fig.savefig(out, dpi=160)
print("wrote", out)
