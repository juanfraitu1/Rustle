#!/usr/bin/env python3
"""Cross-check of the headline containment numbers from the ORIGINAL light tables (P.groups.tsv, D.groups.tsv,
C.groups.tsv) with the verification fixes re-applied by hand here (not via universe.corrected.tsv), plus the TBC1D3 tie
sensitivity (drop one disputed gene at a time). Independent code path from lo_analysis.py.
"""
import csv
import itertools

L = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light"
H = "/mnt/linuxdisk/home/juanfraitu/o1_falsemerge"


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


P = {r["name"]: r["group_id"] for r in tsv(f"{L}/P.groups.tsv") if r["p_status"] == "in_P"}   # fix: no P singletons
D = {r["name"]: r["group_id"] for r in tsv(f"{L}/D.groups.tsv") if r["status"] == "clustered"}  # fix: no D|none rows
# coding genes in the D universe that D.groups.tsv does not list (not in a member cluster): take their raw catalog cluster
idx = {r["name"] for r in tsv(f"{L}/work/P/proteins.index.tsv")}
raw = {(r["chrom"], int(r["start"]), int(r["end"])): r["cluster_id"] for r in tsv(f"{H}/human2/guided.clusters.tsv")}
e1 = {(r["chrom"], int(r["start"]), int(r["end"])): r["cluster_id"] for r in tsv(f"{H}/lit/aj_dev/refseq_e1.clusters.tsv")}
loci = {r["annotation"]: r["representative"] for r in tsv(f"{H}/lit/aj_dev/refseq_e1.loci.tsv")}
genes = {r["name"]: r for r in tsv(f"{L}/work/refseq/genes.tsv") if r["chrom"] in ("chr17", "chr16")}
for n in ("TBC1D26",):
    g = genes[n]
    key = f"{g['chrom']}:{int(g['start0']) + 1}-{g['end']}"
    rep = loci.get(key, key)
    c, r = rep.rsplit(":", 1)
    a, b = map(int, r.split("-"))
    D[n] = f"D|c15_17_22|{e1.get((c, a, b), 'unclustered')}"
    print(f"TBC1D26 record {key} -> representative {rep} -> E1/D cluster {D[n]}")
C = {r["name"]: r for r in tsv(f"{L}/C.groups.tsv") if r["in_reference_trees"] == "yes"}
fam = {r["name"]: r["family"] for r in tsv(f"{L}/members.corrected.tsv")}
side = dict(fam)
for n in ("TBC1D26", "USP6NL", "DHX40"):
    side[n] = "TBC1D3"
side["PKD1"] = "NPIP"


def pairs(lab, S):
    return {(a, b) for a, b in itertools.combinations(sorted(S), 2) if lab[a] == lab[b]}


def c(X, Y, S):
    pY, pX = pairs(Y, S), pairs(X, S)
    return len(pY & pX), len(pY)


coding_in_D = {n for n in D if n in idx}
for n in ("PKD1", "DHX40"):
    P[n] = f"P|other:{n}"
PU = set(P)
DU = set(D) | {"NPIPB10P", "NPIPB14P"}  # non-coding D genes are in D.groups already; set used only with PU below
for f in ("NPIP", "TBC1D3"):
    S = {n for n in PU & set(D) if side.get(n) == f}
    a, b = c(P, D, S), c(D, P, S)
    print(f"{f}: genes in P∩D {len(S)}; c(P⊇D) {a[0]}/{a[1]}; c(D⊇P) {b[0]}/{b[1]}")
    if f == "TBC1D3":
        for drop in ("DHX40", "TBC1D26"):
            S2 = S - {drop}
            a, b = c(P, D, S2), c(D, P, S2)
            print(f"   drop {drop}: c(P⊇D) {a[0]}/{a[1]} = {a[0] / a[1]:.3f}; c(D⊇P) {b[0]}/{b[1]} = {b[0] / b[1]:.3f}")
CL1 = {n: r["clade_L1"] for n, r in C.items()}
for n, r in C.items():
    if r["family"] == "TBC1D3":
        CL1[n] = f"singleton:{n}"  # fix: positional clusters are not clades
S = set(CL1) & PU
print(f"C_L1 vs P: genes {len(S)}; c(P⊇C_L1) {c(P, CL1, S)}; c(C_L1⊇P) {c(CL1, P, S)}")
S = set(CL1) & set(D)
print(f"C_L1 vs D: genes {len(S)}; c(D⊇C_L1) {c(D, CL1, S)}; c(C_L1⊇D) {c(CL1, D, S)}")
