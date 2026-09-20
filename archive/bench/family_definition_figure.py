#!/usr/bin/env python3
"""Figure for advisor interest #1: the read-conflict graph defines multi-copy gene families from IsoSeq reads.

Two panels:
  (A) The conflict graph. Vertices = candidate loci (named gene families + decoys). SOLID green edges = de-tie
      conflict edges (our definition) -> connected components are the families found. DASHED red edges = the
      FALSE-POSITIVE edges AS-tie would add (retrocopies, diverged paralog) that our definition correctly omits.
      Decoy loci that stay isolated (domain-sharers, RFPL, resolvable annotated MAGEA genes) = correct negatives.
  (B) The FP/FN ledger: per candidate, de-conflict reads (ours) vs AS-conflict reads (legacy), coloured by
      verdict (TP/TN), with the 3 AS false positives our definition avoids marked.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

J = json.load(open("/mnt/c/Users/jfris/Desktop/Rustle/bench/family_definition_demo.json"))
MIN_READS = J["operating_point"]["min_reads"]
pc = J["pair_counts"]            # "a~b" -> {shared, n_de, n_as, both_mapq0}
rows = J["rows"]

# build node set + family membership (for colouring) from the demo rows
fam_of = {}
truth_of = {}
for r in rows:
    pass
# reconstruct member lists from the panel embedded in the demo (re-read tsv-style from json rows is not enough);
# instead derive nodes/edges directly from pair_counts.
edges_de = [(k.split("~")[0], k.split("~")[1], v["n_de"]) for k, v in pc.items() if v["n_de"] >= MIN_READS]
edges_as_only = [(k.split("~")[0], k.split("~")[1], v["n_as"]) for k, v in pc.items()
                 if v["n_as"] >= MIN_READS and v["n_de"] < MIN_READS]
nodes = sorted({n for k in pc for n in k.split("~")})

G = nx.Graph()
G.add_nodes_from(nodes)
for a, b, w in edges_de:
    G.add_edge(a, b, w=w)

# components of the de-tie graph = families
comps = [c for c in nx.connected_components(G) if len(c) >= 2]
fam_idx = {}
for i, c in enumerate(comps):
    for n in c:
        fam_idx[n] = i
PAL = plt.cm.tab10.colors

fig = plt.figure(figsize=(15, 7.5))

# ---- Panel A: the conflict graph ----
axA = fig.add_subplot(1, 2, 1)
H = nx.Graph()
H.add_nodes_from(nodes)
H.add_edges_from([(a, b) for a, b, _ in edges_de] + [(a, b) for a, b, _ in edges_as_only])
pos = nx.spring_layout(H, seed=7, k=0.9, iterations=200)
node_colors = [PAL[fam_idx[n] % 10] if n in fam_idx else "#dddddd" for n in nodes]
nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=node_colors, node_size=520,
                       edgecolors="#333333", linewidths=0.8, ax=axA)
nx.draw_networkx_edges(G, pos, edgelist=[(a, b) for a, b, _ in edges_de], edge_color="#1a9850",
                       width=2.4, ax=axA)
nx.draw_networkx_edges(H, pos, edgelist=[(a, b) for a, b, _ in edges_as_only], edge_color="#d73027",
                       width=2.0, style="dashed", ax=axA)
short = {n: (n.replace("LOC", "L").replace("MAGd", "M") if len(n) > 8 else n) for n in nodes}
nx.draw_networkx_labels(G, pos, labels=short, font_size=7.0, ax=axA)
axA.set_title("(A) Read-conflict graph  →  families = connected components\n"
              "green = de-tie conflict edge (ours) · red dashed = AS-tie false-positive edge (avoided)",
              fontsize=10)
axA.axis("off")
axA.legend(handles=[
    Line2D([0], [0], color="#1a9850", lw=2.4, label="de-tie edge (our definition)"),
    Line2D([0], [0], color="#d73027", lw=2.0, ls="--", label="AS-tie edge we correctly OMIT"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#dddddd", markersize=11,
           markeredgecolor="#333", label="isolated locus = correct negative"),
], loc="lower left", fontsize=8, framealpha=0.9)

# ---- Panel B: FP/FN ledger ----
axB = fig.add_subplot(1, 2, 2)
rows_sorted = sorted(rows, key=lambda r: (-r["de_reads"], r["name"]))
names = [r["name"] for r in rows_sorted]
de_r = [max(r["de_reads"], 0.3) for r in rows_sorted]
as_r = [max(r["as_reads"], 0.3) for r in rows_sorted]
y = range(len(names))
axB.barh([i + 0.2 for i in y], de_r, height=0.38, color="#1a9850", label="de-conflict reads (ours)")
axB.barh([i - 0.2 for i in y], as_r, height=0.38, color="#fc8d59", label="AS-conflict reads (legacy)")
axB.axvline(MIN_READS, color="#333", ls=":", lw=1.2)
axB.text(MIN_READS, len(names) - 0.5, f" edge threshold ={MIN_READS}", fontsize=7.5, va="top")
axB.set_xscale("symlog")
axB.set_yticks(list(y))
ylabels = []
for r in rows_sorted:
    tag = {"TP": "✓ family", "TN": "· not", "FP": "✗FP", "FN": "✗FN"}[r["verdict"]]
    avoided = (not r["truth"]) and r["as_reads"] >= MIN_READS
    ylabels.append(f"{r['name']}  [{tag}]" + ("  ⟵AS-FP avoided" if avoided else ""))
axB.set_yticklabels(ylabels, fontsize=7.2)
axB.set_xlabel("conflict-supporting reads (symlog)")
axB.set_title("(B) per-candidate ledger — de-tie: 7 TP / 10 TN / 0 FP / 0 FN\n"
              "AS-tie over-links APOBEC3, EEF1A1, CNN2 (orange ≫ threshold where green is 0)", fontsize=10)
axB.legend(loc="lower right", fontsize=8)
axB.invert_yaxis()

fig.suptitle("Defining a multi-copy gene family from IsoSeq reads: the read-conflict graph (de-tie criterion)",
             fontsize=12, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/family_definition_figure.png"
fig.savefig(out, dpi=160, bbox_inches="tight")
print(f"wrote {out}")
print(f"families (de-tie components): {[sorted(c) for c in comps]}")
print(f"AS-only false-positive edges: {[(a, b) for a, b, _ in edges_as_only]}")
