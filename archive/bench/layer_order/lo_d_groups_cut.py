#!/usr/bin/env python3
"""Are the member D groups exact components of the D graph? For each member D group: distinct genes outside it with a
D edge to it, and gene-level edges leaving it (distinct (inside gene, outside gene) pairs), from D.edges.corrected rows
(rows touching member-group keys are all in D.edges, layer_dna.py rule). Also the raw catalog edges of the TBC1D26 record
(audit 2026-09-16: TBC1D26 is labelled MCL24 through the ZNF286A-TBC1D26 readthrough fold, not through a TBC1D3 edge).
Output: integrate_slim/d_groups_cut.out"""
import collections, csv
L = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light"
rows = list(csv.DictReader(open(f"{L}/D.edges.corrected.tsv"), delimiter="\t"))
grp = {r["gene_id"]: r["group_id"] for r in csv.DictReader(open(f"{L}/D.groups.corrected.tsv"), delimiter="\t")
       if r["status"] == "clustered"}
out = collections.defaultdict(dict)
gene_edges = collections.defaultdict(set)
inside = collections.Counter()
for r in rows:
    gu, gv = grp.get(r["u_gene_id"]), grp.get(r["v_gene_id"])
    if gu and gu == gv:
        inside[gu] += 1
        continue
    for g, iid, oid, oname, ogrp in ((gu, r["u_gene_id"], r["v_gene_id"], r["v_name"], r["v_group"]),
                                     (gv, r["v_gene_id"], r["u_gene_id"], r["u_name"], r["u_group"])):
        if g and oid != "?":
            k = f"{oname}({ogrp or 'unclustered'})"
            out[g][k] = max(out[g].get(k, 0.0), float(r["weight"]))
            gene_edges[g].add((iid, oid))
for g in sorted(set(grp.values())):
    o = sorted(out[g].items(), key=lambda kv: -kv[1])
    print(f"{g}: within-group edge rows {inside[g]}; distinct outside genes with an edge {len(o)}; gene-level edges leaving "
          f"the group {len(gene_edges[g])}; strongest: " + ", ".join(f"{k} {v:.3f}" for k, v in o[:8]))
genes = {r["gene_id"]: r for r in csv.DictReader(open(f"{L}/work/refseq/genes.tsv"), delimiter="\t")}
by_key = collections.defaultdict(list)
for g in genes.values():
    by_key[f"{g['chrom']}:{int(g['start0']) + 1}-{g['end']}"].append(g["name"])
t = genes["gene-TBC1D26"]
key = f"{t['chrom']}:{int(t['start0']) + 1}-{t['end']}"
print(f"TBC1D26 record {key}: raw catalog edges (light/work/D/c15_17_22.e1.graph.tsv):")
for line in open(f"{L}/work/D/c15_17_22.e1.graph.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) == 3 and key in (f[0], f[1]):
        other = f[1] if f[0] == key else f[0]
        print(f"   {other} ({','.join(by_key.get(other, ['?']))}) weight {f[2]}")
loci = {r["annotation"]: r["representative"] for r in csv.DictReader(open(f"{L}/work/D/c15_17_22.e1.loci.tsv"), delimiter="\t")}
rep = loci.get(key, key)
print(f"TBC1D26 record folded into {rep} ({','.join(by_key.get(rep, ['?']))}); records folded into that locus: "
      + ", ".join(f"{a} ({','.join(by_key.get(a, ['?']))})" for a, b in loci.items() if b == rep))
