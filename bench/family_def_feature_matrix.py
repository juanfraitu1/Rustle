#!/usr/bin/env python3
"""family_def_feature_matrix.py — rich per-edge feature matrix for UNBIASED differentiator discovery.

Population: both-spliced cDNA-homology edges (retrocopies already removed; both-intronless excluded
since junction co-threading is undefined there). For each edge (a<b) we compute every cheap structural
feature we can, INCLUDING ones not yet tried (graph topology: edge betweenness = the canonical
bridge signal, common-neighbours, degree, neighbourhood Jaccard), and merge the co-threading features.

Labels are INDEPENDENT of the structural features and used ONLY to orient / sanity-check (never to
define families here -- avoids the circularity that sank earlier attempts):
  compara_pos  : both genes in one Compara checkable paralog family (external positive)
  panel_real / panel_counter : the hand-validated airtight panel
  group        : connected-component id (post-retrocopy graph) -> for GROUPED CV (hold a whole
                 component out, so a classifier cannot memorise a family)

Run: /home/juanfra/miniforge3/bin/python bench/family_def_feature_matrix.py
Out: bench/family_def_features.tsv
"""
import collections
import json
import math
import os
import sys

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_retrocopy_filter import build_family_graph, load_intron_index
from family_def_strand_probe import edge_strand

GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
COTHREAD_TSV = os.path.join(HERE, "family_def_cothread_genomewide.tsv")
COMPARA = os.path.join(HERE, "compara_paralog_relation.json")
PANEL_TSV = os.path.join(HERE, "family_def_airtight_panel.tsv")
NA = ""   # missing marker


def load_cothread():
    w = {}
    if not os.path.exists(COTHREAD_TSV):
        return w
    with open(COTHREAD_TSV) as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for line in f:
            t = line.rstrip("\n").split("\t")
            e = tuple(sorted((t[ix["ref"]], t[ix["mem"]])))
            w[e] = dict(frac_ref=float(t[ix["frac_ref"]]), frac_mem=float(t[ix["frac_mem"]]),
                        aln_cov_small=float(t[ix["aln_cov_small"]]), contig=int(t[ix["contig"]]),
                        co=int(t[ix["co"]]), ref_j=int(t[ix["ref_j"]]), mem_j=int(t[ix["mem_j"]]))
    return w


def compara_pos_pairs():
    pos = set()
    try:
        cj = json.load(open(COMPARA))
    except (OSError, ValueError):
        return pos
    for fam, info in cj.get("families", {}).items():
        if not info.get("checkable_pair"):
            continue
        gs = info.get("mapped_genes", [])
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                pos.add(tuple(sorted((gs[i], gs[j]))))
    return pos


def panel_pairs():
    real, counter = set(), set()
    if not os.path.exists(PANEL_TSV):
        return real, counter
    with open(PANEL_TSV) as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for line in f:
            t = line.rstrip("\n").split("\t")
            mem = t[ix["members"]].split(",")
            dst = real if t[ix["kind"]] == "real" else counter
            for i in range(len(mem)):
                for j in range(i + 1, len(mem)):
                    dst.add(tuple(sorted((mem[i], mem[j]))))
    return real, counter


def main():
    Hd, _ = dna_homology()
    idx = load_intron_index()
    strand = edge_strand()
    ct = load_cothread()
    compos = compara_pos_pairs()
    preal, pcounter = panel_pairs()
    coord = {}
    for line in open(GENES_BED):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            coord[p[3]] = (p[0], int(p[1]), int(p[2]))

    G, prov, part = build_family_graph(Hd, idx, filter_retrocopy=True)
    both_spliced = [tuple(sorted(e)) for e in part["both_spliced_kept"]]

    # component id per node + per-component edge betweenness
    comp_id = {}
    betw = {}
    for ci, c in enumerate(nx.connected_components(G)):
        for g in c:
            comp_id[g] = ci
        if len(c) >= 3:
            sub = G.subgraph(c)
            eb = nx.edge_betweenness_centrality(sub, normalized=True)
            for (x, y), v in eb.items():
                betw[tuple(sorted((x, y)))] = v
        elif len(c) == 2:
            x, y = list(c)
            betw[tuple(sorted((x, y)))] = 1.0

    def clen(g):
        v = idx.get(g)
        return sum(v["exon_lengths"]) if v and v.get("exon_lengths") else None

    def njunc(g):
        v = idx.get(g)
        return v.get("n_intron") if v else None

    cols = ["a", "b", "seq_id", "cov_min", "cov_max", "len_min", "len_max", "len_ratio",
            "njunc_a", "njunc_b", "njunc_min", "njunc_diff", "njunc_ratio",
            "frac_ref", "frac_mem", "aln_cov_small", "contig", "co_junc",
            "antisense", "same_chrom", "log_dist", "overlap", "disjoint",
            "edge_betw", "jaccard_nbr", "common_nbr", "deg_min", "deg_max",
            "compara_pos", "panel_real", "panel_counter", "group"]
    rows = []
    for a, b in both_spliced:
        r = Hd[(a, b)] if (a, b) in Hd else Hd[(b, a)]
        la, lb = clen(a), clen(b)
        na, nb = njunc(a), njunc(b)
        ca, cb = coord.get(a), coord.get(b)
        same_chrom = int(bool(ca and cb and ca[0] == cb[0]))
        overlap = int(bool(ca and cb and ca[0] == cb[0] and min(ca[2], cb[2]) - max(ca[1], cb[1]) > 0))
        if ca and cb and ca[0] == cb[0]:
            dist = max(0, max(ca[1], cb[1]) - min(ca[2], cb[2]))
            log_dist = math.log10(dist + 1)
        else:
            log_dist = NA
        e = (a, b)
        cof = ct.get(e, {})
        # topology
        na_nbr = set(G.neighbors(a)) - {b}
        nb_nbr = set(G.neighbors(b)) - {a}
        inter = na_nbr & nb_nbr
        union = na_nbr | nb_nbr
        jac = len(inter) / len(union) if union else 0.0
        row = {
            "a": a, "b": b,
            "seq_id": round(r["id"], 4),
            "cov_min": round(min(r["cov_a"], r["cov_b"]), 4),
            "cov_max": round(max(r["cov_a"], r["cov_b"]), 4),
            "len_min": min(la, lb) if la and lb else NA,
            "len_max": max(la, lb) if la and lb else NA,
            "len_ratio": round(min(la, lb) / max(la, lb), 4) if la and lb else NA,
            "njunc_a": na if na is not None else NA,
            "njunc_b": nb if nb is not None else NA,
            "njunc_min": min(na, nb) if na is not None and nb is not None else NA,
            "njunc_diff": abs(na - nb) if na is not None and nb is not None else NA,
            "njunc_ratio": round(min(na, nb) / max(na, nb), 4) if na and nb else NA,
            "frac_ref": cof.get("frac_ref", NA),
            "frac_mem": cof.get("frac_mem", NA),
            "aln_cov_small": cof.get("aln_cov_small", NA),
            "contig": cof.get("contig", NA),
            "co_junc": cof.get("co", NA),
            "antisense": int(strand.get(e) == "-"),
            "same_chrom": same_chrom,
            "log_dist": round(log_dist, 3) if log_dist != NA else NA,
            "overlap": overlap,
            "disjoint": int(not overlap),
            "edge_betw": round(betw.get(e, NA), 6) if betw.get(e) is not None else NA,
            "jaccard_nbr": round(jac, 4),
            "common_nbr": len(inter),
            "deg_min": min(G.degree(a), G.degree(b)),
            "deg_max": max(G.degree(a), G.degree(b)),
            "compara_pos": int(e in compos),
            "panel_real": int(e in preal),
            "panel_counter": int(e in pcounter),
            "group": comp_id.get(a, -1),
        }
        rows.append(row)

    with open(os.path.join(HERE, "family_def_features.tsv"), "w") as out:
        out.write("\t".join(cols) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in cols) + "\n")
    npos = sum(r["compara_pos"] for r in rows)
    print(f"[+] {len(rows)} both-spliced edges x {len(cols)-3} features -> bench/family_def_features.tsv")
    print(f"    independent labels: compara_pos={npos}, panel_real={sum(r['panel_real'] for r in rows)}, "
          f"panel_counter={sum(r['panel_counter'] for r in rows)}")
    print(f"    groups (components): {len(set(r['group'] for r in rows))}")
    print(f"    edges missing co-threading: {sum(1 for r in rows if r['frac_ref']==NA)}")


if __name__ == "__main__":
    main()
