#!/usr/bin/env python3
"""Figure for psv_graph_demo: the PSV-aware variation graph (copies = paths braiding
through PSV bubbles) + the read assignment, for a resolvable family (RABL2) and a
K-frontier family (RFPL4A). Reads psv_graph_demo.json.
Run: /home/juanfra/miniforge3/bin/python bench/psv_graph_figure.py
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NAVY, TEAL, ORANGE, GREY = "#1F2D5A", "#127A6E", "#D96B27", "#9aa0a6"
BASE_Y = {"A": 0.9, "C": 0.3, "G": -0.3, "T": -0.9, "N": 0.0, "-": 0.0}
d = {f["label"]: f for f in json.load(open(os.path.join(HERE, "psv_graph_demo.json")))}


def draw(ax, fam, path_styles, nbub=9):
    haps = fam["haplotypes"]; names = fam["copy_names"]
    cols = min(nbub, len(next(iter(haps.values()))))
    xs = [1 + i * 1.15 for i in range(cols)]
    # bubble guides + allele letters
    for i, x in enumerate(xs):
        ax.plot([x, x], [-1.25, 1.25], color="#eee", lw=8, zorder=0)
        for a in sorted({haps[nm][i] for nm in names}):
            ax.text(x, BASE_Y.get(a, 0), a, ha="center", va="center", fontsize=8.5,
                    color="#555", zorder=2,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="#ccc", lw=0.6))
    ax.text(xs[0] - 0.85, 0, "copies →", ha="right", va="center", fontsize=9, color="#666")
    # paths: each drawn group once
    drawn = set()
    for nm in names:
        grp = next(tuple(g) for g in fam["groups"] if nm in g)
        if grp in drawn:
            continue
        drawn.add(grp)
        ys = [BASE_Y.get(haps[nm][i], 0) for i in range(cols)]
        st = path_styles(grp)
        ax.plot(xs, ys, color=st["c"], lw=st["lw"], alpha=0.9, zorder=3,
                solid_capstyle="round")
        ax.text(xs[-1] + 0.35, ys[-1], st["label"], va="center", fontsize=9.5,
                color=st["c"], weight="bold")
    ax.set_xlim(0, xs[-1] + 3.0); ax.set_ylim(-1.7, 1.6); ax.axis("off")


fig, ax = plt.subplots(2, 1, figsize=(13, 8.6))

# ---- RABL2: resolvable ----
fam = d["RABL2 (2 copies)"]
def rabl2_style(grp):
    return {"RABL2A": dict(c=TEAL, lw=2.6, label="RABL2A"),
            "RABL2B": dict(c=ORANGE, lw=2.6, label="RABL2B")}[grp[0]]
draw(ax[0], fam, rabl2_style)
pc = fam["per_copy"]
ax[0].set_title(f"RESOLVABLE — RABL2 (2 copies):  the paths separate at the PSV bubbles, so reads get a copy",
                fontsize=12.5, weight="bold", color=NAVY, loc="left")
ax[0].text(0, -1.55,
           f"{fam['psvs']} PSV bubbles.  {sum(pc.values())} reads threaded to a single copy "
           f"(RABL2A {pc.get('RABL2A',0)}, RABL2B {pc.get('RABL2B',0)}) — {fam['validation']} agree with best-mapping copy (100%).   "
           f"{fam['assignment'].get('unassignable (no PSV covered)',0)} reads stay on the shared backbone → unassignable (K-frontier).",
           fontsize=9.2, color="#333")

# ---- RFPL4A: K-frontier ----
fam = d["RFPL4A array (5 copies)"]
def rfpl_style(grp):
    if grp == ("RFPL4A",):
        return dict(c=NAVY, lw=2.6, label="RFPL4A")
    return dict(c=GREY, lw=6.0, label="copy2–5  (×4, identical)")
draw(ax[1], fam, rfpl_style)
a = fam["assignment"]
ax[1].set_title("K-FRONTIER — RFPL4A array (5 copies):  4 copies share every path → the graph shows they are indistinguishable",
                fontsize=12.5, weight="bold", color=NAVY, loc="left")
ax[1].text(0, -1.55,
           f"{fam['psvs']} PSV bubbles, but only {fam['resolvable_groups']} resolvable haplotypes: RFPL4A vs the collapsed "
           f"{{copy2–5}}.   {a.get('resolved to ONE copy',0)} reads → RFPL4A;  {a.get('resolved to a COLLAPSED group',0)} reads → "
           f"the cluster (can't pick which of 4 — they take the same allele at every bubble);  "
           f"{a.get('unassignable (no PSV covered)',0)} unassignable.",
           fontsize=9.2, color="#333")

fig.suptitle("A PSV-aware variation graph assigns reads to copies — and exposes where it is impossible (identifiability)",
             fontsize=13.5, weight="bold", color=NAVY, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(os.path.join(HERE, "psv_graph.png"), dpi=150, bbox_inches="tight")
print("[+] wrote psv_graph.png")
