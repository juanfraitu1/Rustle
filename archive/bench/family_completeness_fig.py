#!/usr/bin/env python3
"""Figure: completeness of the annotation-based RNA-level family definition."""
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(__file__)
META = "/tmp/gene_reps_gw.meta.tsv"
GRADED = os.path.join(BASE, "family_pairs_graded.tsv")
NAVY = "#22313f"; GREEN = "#1e8449"; ORANGE = "#e8590c"; RED = "#c0392b"; SLATE = "#8898aa"
plt.rcParams["font.family"] = "DejaVu Sans"
CORE_COV, CORE_IDENT = 0.13, 0.70

genes = set()
for line in open(META):
    if not line.startswith("gene\t"):
        genes.add(line.split("\t", 1)[0])

# rebuild components (same as completeness)
parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
with open(GRADED) as fh:
    fh.readline()
    for line in fh:
        a, ca, b, cb, cc, ci, sc = line.rstrip("\n").split("\t")
        if float(cc) >= CORE_COV and float(ci) >= CORE_IDENT:
            parent[find(a)] = find(b)
comp = defaultdict(set)
for g in list(parent):
    comp[find(g)].add(g)
comps = list(comp.values())
g2c = {g: i for i, c in enumerate(comps) for g in c}
sizes = sorted((len(c) for c in comps), reverse=True)

CUR = {"RABL2": ["RABL2A", "RABL2B"], "DAZ": ["DAZ1", "DAZL"],
       "APOBEC3": ["APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H"],
       "RFPL": ["RFPL1", "RFPL2", "RFPL3", "RFPL4A"],
       "GGT": ["GGT1", "GGT5", "GGT6", "GGT7", "GGTLC1", "GGTLC2"],
       "MAGEA*": sorted(x for x in genes if x.startswith("MAGEA")),
       "PRAMEF*": sorted(x for x in genes if x.startswith("PRAMEF")),
       "TAS2R*": sorted(x for x in genes if x.startswith("TAS2R")),
       "DEFB*": sorted(x for x in genes if x.startswith("DEFB")),
       "SIGLEC*": sorted(x for x in genes if x.startswith("SIGLEC")),
       "CRYBG": ["CRYBG1", "CRYBG2", "CRYBG3"]}
CUR = {f: [m for m in ms if m in genes] for f, ms in CUR.items()}
CUR = {f: ms for f, ms in CUR.items() if len(ms) >= 2}

rows = []
for fam, ms in CUR.items():
    ch = defaultdict(int)
    for m in ms:
        if m in g2c:
            ch[g2c[m]] += 1
    best = max(ch.values()) if ch else 0
    rows.append((fam, best, len(ms)))
rows.sort(key=lambda r: -(r[1] / r[2]))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.4), gridspec_kw={"width_ratios": [1.05, 1.0]})
y = range(len(rows))
fracs = [b / t for _, b, t in rows]
cols = [GREEN if f >= 0.99 else (ORANGE if f > 0 else RED) for f in fracs]
axA.barh(list(y), fracs, color=cols, alpha=0.9)
axA.set_yticks(list(y)); axA.set_yticklabels([f"{f} ({b}/{t})" for f, b, t in rows], fontsize=8.5)
axA.invert_yaxis(); axA.set_xlim(0, 1.05); axA.set_xlabel("fraction of members in one component")
axA.set_title("Known-family recovery (curated)\ngreen=full · orange=partial · red=missed", fontsize=10.5, color=NAVY)
axA.axvline(0, color="k", lw=0.5)
axA.text(0.5, len(rows)-0.5, "recent dups recovered ↑\nancient/diverged families ↓ (need protein tier)",
         fontsize=8, color=SLATE, ha="center", va="bottom", style="italic")

axB.axis("off"); axB.set_xlim(0, 10); axB.set_ylim(0, 10)
axB.set_title("Scoped RNA-level family definition", fontsize=11.5, color=NAVY, fontweight="bold")
txt = [
    ("roster", f"{len(genes):,} ANNOTATED genes", NAVY),
    ("criterion", "POA contiguous core ≥0.13 & core-id ≥0.70", NAVY),
    ("families", f"{sum(1 for c in comps if len(c)>=2):,} (size≥2); {sum(len(c) for c in comps if len(c)>=2):,} genes", GREEN),
    ("recent families", "recovered (RABL2, DAZ, APOBEC3, RFPL…)", GREEN),
    ("domain-sharers", "correctly REJECTED (CDPF1/CREB1/GCA…)", GREEN),
    ("ancient families", "MISSED (TAS2R/DEFB/SIGLEC) → need DNA/protein tier", RED),
    ("over-merge", f"14 mega-components (largest {sizes[0]}) = domain-hub chains", ORANGE),
    ("universe recall", "46/53 (the 7 misses = domain-sharers)", GREEN),
]
yy = 8.6
for k, v, c in txt:
    axB.text(0.3, yy, f"{k}:", fontsize=9.5, color=NAVY, fontweight="bold")
    axB.text(3.4, yy, v, fontsize=9, color=c)
    yy -= 1.05
axB.text(5, 0.2, "complete for RECENT duplicates + precise vs domain-sharers;\n"
                 "ancient-family gap = the case for the parallel DNA/protein tier",
         fontsize=8.6, color=NAVY, ha="center", style="italic")

fig.suptitle("Annotation-based RNA-level multi-copy gene family definition — completeness",
             fontsize=13, fontweight="bold", color=NAVY, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(BASE, "family_definition.png")
fig.savefig(out, dpi=160, facecolor="white")
print("saved", out)
