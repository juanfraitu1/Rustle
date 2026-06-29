#!/usr/bin/env python3
"""Flagship topological figure: a multi-copy family is a CLIQUE (backbone graph) whose copies are
PATHS through a PSV variation graph; identifiability = min-path-cover = χ(H). RABL2 (2 copies -> 2
paths, fully resolvable) vs RFPL4A (5 copies -> 2 paths, the K-frontier: copies 2-5 are PSV-identical
and collapse to ONE path). Data from bench/psv_graph_demo.json."""
import json, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

d = {f["label"]: f for f in json.load(open("bench/psv_graph_demo.json"))}
BY = {"A": 0, "C": 1, "G": 2, "T": 3}
fig, axes = plt.subplots(2, 2, figsize=(15, 8), gridspec_kw={"height_ratios": [1, 1.3]})

def panel(col, key, title, sub):
    f = d[key]
    names = list(f["haplotypes"].keys())
    haps = f["haplotypes"]
    # group copies by identical PSV-haplotype (= distinct paths = χ(H))
    groups = {}
    for nm in names:
        groups.setdefault(haps[nm], []).append(nm)
    chi = len(groups)
    palette = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]
    grp_color = {h: palette[i % len(palette)] for i, h in enumerate(groups)}

    # ---- TOP: the CLIQUE (family = dense reciprocal-backbone graph) ----
    axc = axes[0, col]
    n = len(names)
    ang = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, n, endpoint=False)
    pos = {nm: (np.cos(a), np.sin(a)) for nm, a in zip(names, ang)}
    for a, b in itertools.combinations(names, 2):   # clique: every pair has a ~B edge
        axc.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], "-", color="0.7", lw=1, zorder=1)
    for nm in names:
        c = grp_color[haps[nm]]
        axc.add_patch(Circle(pos[nm], 0.14, color=c, zorder=3))
        axc.text(pos[nm][0] * 1.32, pos[nm][1] * 1.32, nm, ha="center", va="center", fontsize=8)
    axc.set_xlim(-1.7, 1.7); axc.set_ylim(-1.6, 1.6); axc.set_aspect("equal"); axc.axis("off")
    axc.set_title(f"{title}\nfamily = clique (~B backbone): {n} copies", fontsize=10)

    # ---- BOTTOM: the VARIATION GRAPH (copies = paths through PSV bubbles) ----
    axg = axes[1, col]
    L = len(next(iter(haps.values())))
    # draw each copy's path; identical copies overlap (collapse), alpha stacks -> bolder
    for nm in names:
        h = haps[nm]
        xs = list(range(L)); ys = [BY.get(b, 1.5) for b in h]
        axg.plot(xs, ys, "-", color=grp_color[h], lw=2.2, alpha=0.45, zorder=2)
    # bubble markers: columns where >=2 alleles among copies
    for x in range(L):
        col_bases = {haps[nm][x] for nm in names}
        if len(col_bases) >= 2:
            axg.axvspan(x - 0.3, x + 0.3, color="0.92", zorder=0)
    axg.set_yticks([0, 1, 2, 3]); axg.set_yticklabels(["A", "C", "G", "T"])
    axg.set_xlabel(f"PSV column (bubble)   —   {L} PSV sites")
    axg.set_xlim(-1, L); axg.set_ylim(-0.6, 3.6)
    # legend: one entry per distinct path (group)
    for h, members in groups.items():
        lab = members[0] if len(members) == 1 else f"{{{members[0]}..{members[-1]}}} ×{len(members)}"
        axg.plot([], [], "-", color=grp_color[h], lw=3, label=lab)
    axg.legend(title=f"χ(H) = {chi} distinct path(s)", fontsize=8, loc="upper right", framealpha=0.9)
    axg.set_title(sub, fontsize=10, color="#b00" if chi < n else "#070")

panel(0, "RABL2 (2 copies)", "RABL2",
      "2 copies → 2 paths  (χ(H)=2 = #copies)\nFULLY RESOLVABLE · 58/58 reads validated")
panel(1, "RFPL4A array (5 copies)", "RFPL4A",
      "5 copies → 2 paths  (χ(H)=2 < 5)\nK-FRONTIER: copies 2–5 PSV-identical → collapse to ONE path")

fig.suptitle("A multi-copy family is a clique whose copies are paths through a PSV variation graph;\n"
             "identifiability = minimum path-cover = χ(H).  The K-frontier is a graph property.",
             fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig("bench/flagship_topology.png", dpi=140)
print("wrote bench/flagship_topology.png")
print(f"RABL2:  {len(d['RABL2 (2 copies)']['haplotypes'])} copies -> χ(H)=2")
print(f"RFPL4A: {len(d['RFPL4A array (5 copies)']['haplotypes'])} copies -> χ(H)=2 (K-frontier)")
