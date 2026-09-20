#!/usr/bin/env python3
"""Nested edge-test lattice, step 3: levels G_0..G_3 (primary tests in lattice_common.tests), their 3-truss (triangle)
variant, T1/T2 sanity checks, chaining evidence, operationalisation variants.

Inputs: lattice/nodes.tsv, lattice/edges.tsv (lattice_edges.py), lattice/expr_counts.tsv (lattice_expr.py).
Outputs (lattice/):
  groups.tsv          per gene: component label (gene_id of the component representative) per level x variant
  levels.tsv          per level x variant: edges, components, largest, member-holding groups
  member_groups.tsv   every member-holding group per level x variant: size, members, pulled-in non-members
  sanity.tsv          T1 (nesting across levels) and T2 (expression views) violation counts (must be 0), with the
                      non-vacuity counts: coarse blocks with >= 2 genes, how many the finer partition splits, member-holding
  expr_views.tsv      member-holding expressed components per expression set x variant x level, including both 3-truss
                      views comp(truss(G_k[X])) and comp(truss(G_k)[X])
  chaining.tsv        genes of interest (PKD1 / readthrough / DHX40 / RNFT1 / TBC-domain neighbours): group per level and
                      the shortest connecting path to the family anchor with its edge attributes
  triangle_drops.tsv  member-holding groups of G_k that lose genes in the 3-truss variant (2-copy groups included)
  levels.out          log
"""
import collections
import csv
import sys
import os
# §6r9: repo root from THIS file, so the tool runs from any clone (it used to hardcode
# /mnt/c/Users/jfris/Desktop/Rustle, which only ever worked on one machine).
_RUSTLE_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import time

sys.path.insert(0, os.path.join(_RUSTLE_REPO, 'bench', 'layer_order'))
from lattice_common import (LEVELS, OUT, components, fnum, groups, refines, split_counts, tests, truss3, tsv,  # noqa: E402
                            write)

T0 = time.time()
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


nodes = {r["gene_id"]: r for r in tsv(f"{OUT}/nodes.tsv")}
NAME = {g: r["name"] for g, r in nodes.items()}
MEM = {g for g, r in nodes.items() if r["is_member"] == "yes"}
FAM = {g: r["member_family"] for g, r in nodes.items() if r["is_member"] == "yes"}
_ec = tsv(f"{OUT}/expr_counts.tsv")
reads = {r["gene_id"]: int(r["n_reads_any"]) for r in _ec}
reads_u = {r["gene_id"]: int(r["n_reads_unique"]) for r in _ec}
V = set(nodes)

# ---------------------------------------------------------------------------------------------- edges (compact)
NEED = ["gene_a", "gene_b", "same_locus", "p_qualifies_6ko", "p_aa_identity", "d_c2_approx", "d_e1_edge", "d_e1_identity",
        "d_e1_cov_longer", "d_shared_exon_frac", "d_e1_identity_gapexcl", "s2_max_identity", "d_c2_genebody", "d_c2_exon",
        "p_cov_longer", "d_c2_gb_chain_identity", "d_c2_gb_best_frac", "d_c2_exon_best_frac", "d_c2x_approx",
        "d_c2nostrand_approx", "d_c2exontgt_approx", "d_c2loose_approx", "d_w98_gapexcl", "d_w98_gapincl"]
ROWS = []
with open(f"{OUT}/edges.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t")
    hdr = next(rd)
    ix = {c: hdr.index(c) for c in NEED}
    for f in rd:
        ROWS.append({c: f[i] for c, i in ix.items()})
say(f"[load] nodes {len(V)}; edge rows {len(ROWS)}; {time.time() - T0:.0f}s")
EDGE = {(r["gene_a"], r["gene_b"]): r for r in ROWS}

VARIANTS = {  # primary: L1 = clause-2 approx with v-exon overlap + strand check; L3 = single-record w_98, gap-excluded
    "primary": dict(),
    "with_same_locus": dict(with_same_locus=True),
    "L0=(P_and_aa>=0.50)_or_D": dict(p_aa_min=0.50),
    "L1=c2_no_strand_check": dict(l1="c2_nostrand"),
    "L1=c2_vexon_on_exon_proxy_only": dict(l1="c2_exontgt"),
    "L1=c2_no_vexon_no_strand": dict(l1="c2_loose"),
    "L1=c2x_extrapolated": dict(l1="c2x"),
    "L1=E1_as_built(D graph)": dict(l1="e1"),
    "L1=E1_at_0.80/0.50": dict(l1="e1c2"),
    "L3=w98_gap-inclusive": dict(id_which="w98_gapincl"),
    "L3=pooled_gap-excluded": dict(id_which="pooled_gapexcl"),
    "L3=pooled_gap-inclusive": dict(id_which="pooled_gapincl"),
    "L3=S2_SD98_mapback": dict(id_which="s2"),
    "17:03_tests_exact": dict(l1="c2_loose", id_which="pooled_gapexcl"),  # the 17:12 report's tests, unrounded thresholds
}


def level_edges(**kw):
    E = {k: [] for k in LEVELS}
    for r in ROWS:
        t = tests(r, **kw)
        for k, ok in zip(LEVELS, t):
            if ok:
                E[k].append((r["gene_a"], r["gene_b"]))
    return E


def summarize(tag, lab, E, nodeset, rows_levels, rows_groups):
    G = groups(lab)
    nonsingle = [S for S in G.values() if len(S) >= 2]
    mem_groups = [S for S in G.values() if S & MEM]
    rows_levels.append({"variant": tag[0], "level": tag[1], "nodes": len(nodeset), "edges": len(E),
                        "components_ge2": len(nonsingle), "singletons": sum(1 for S in G.values() if len(S) == 1),
                        "largest": max(len(S) for S in G.values()),
                        "member_groups": len(mem_groups),
                        "member_groups_sizes": ";".join(f"{len(S)}(NPIP {sum(1 for g in S & MEM if FAM[g] == 'NPIP')},TBC1D3 "
                                                        f"{sum(1 for g in S & MEM if FAM[g] == 'TBC1D3')})"
                                                        for S in sorted(mem_groups, key=lambda s: (-len(s), sorted(s))))})
    for S in sorted(mem_groups, key=lambda s: (-len(s), sorted(s))):
        mems = sorted(NAME[g] for g in S & MEM)
        non = sorted(S - MEM, key=lambda g: (int(nodes[g]["l0_hops_from_U"]), NAME[g]))
        bt = collections.Counter(nodes[g]["biotype"] for g in S)
        ch = collections.Counter(nodes[g]["chrom"] for g in S)
        rows_groups.append({"variant": tag[0], "level": tag[1], "group_rep": NAME[min(S)], "size": len(S),
                            "n_NPIP_members": sum(1 for g in S & MEM if FAM[g] == "NPIP"),
                            "n_TBC1D3_members": sum(1 for g in S & MEM if FAM[g] == "TBC1D3"),
                            "members": ",".join(mems), "n_nonmembers": len(non),
                            "nonmembers_first80": ",".join(NAME[g] for g in non[:80]),
                            "n_outside_E1_catalogs": sum(1 for g in S if nodes[g]["catalog"] == "none"),
                            "n_protein_never_searched": sum(1 for g in S if nodes[g]["protein_searched"] == "no"),
                            "biotypes": ";".join(f"{k}:{v}" for k, v in bt.most_common()),
                            "chroms": ";".join(f"{k}:{v}" for k, v in ch.most_common())})


rows_levels, rows_groups, gl_rows = [], [], {g: {"gene_id": g, "name": NAME[g], "is_member": nodes[g]["is_member"],
                                                   "member_family": nodes[g]["member_family"],
                                                   "n_reads_any": reads.get(g, 0)} for g in V}
LAB = {}
EDGES = {}
for vname, kw in VARIANTS.items():
    E = level_edges(**kw)
    EDGES[vname] = E
    for k in LEVELS:
        lab = components(V, E[k])
        LAB[(vname, k)] = lab
        summarize((vname, k), lab, E[k], V, rows_levels, rows_groups)
        for g in V:
            gl_rows[g][f"{vname}|{k}"] = lab[g]
    say(f"[levels] {vname}: edges " + ", ".join(f"{k} {len(E[k])}" for k in LEVELS) + f"; {time.time() - T0:.0f}s")

# ---- triangle (3-truss) variant of the primary levels (and of the 17:12 report's tests, unrounded)
TRI = {}
for base, tag in (("primary", "triangle"), ("17:03_tests_exact", "17:03_tests_exact_triangle")):
    for k in LEVELS:
        kept, dropped, depth = truss3(EDGES[base][k])
        TRI[(tag, k)] = kept
        lab = components(V, kept)
        LAB[(tag, k)] = lab
        summarize((tag if tag != "triangle" else "triangle(3-truss)", k), lab, kept, V, rows_levels, rows_groups)
        for g in V:
            gl_rows[g][f"{tag}|{k}"] = lab[g]
        say(f"[{tag}] {k}: edges {len(EDGES[base][k])} -> {len(kept)} (dropped {dropped}, peel depth {depth}); "
            f"{time.time() - T0:.0f}s")

# ---- node variant: readthrough records removed (clause 1: readthrough spans must not be nodes)
RT = {g for g in V if nodes[g]["readthrough"] == "yes"}
V_nrt = V - RT
for k in LEVELS:
    E = [e for e in EDGES["primary"][k] if e[0] in V_nrt and e[1] in V_nrt]
    lab = components(V_nrt, E)
    LAB[("no_readthrough_nodes", k)] = lab
    summarize(("no_readthrough_nodes", k), lab, E, V_nrt, rows_levels, rows_groups)
    for g in V:
        gl_rows[g][f"no_readthrough_nodes|{k}"] = lab[g] if g in lab else "removed"
say(f"[no-readthrough] readthrough records in V removed: {len(RT)} ({sorted(NAME[g] for g in RT & MEM)} are members)")

write(f"{OUT}/levels.tsv", rows_levels)
write(f"{OUT}/member_groups.tsv", rows_groups)
gcols = ["gene_id", "name", "is_member", "member_family", "n_reads_any"] + [c for c in next(iter(gl_rows.values())) if "|" in c]
write(f"{OUT}/groups.tsv", [gl_rows[g] for g in sorted(V, key=lambda g: NAME[g])], gcols)

# ---------------------------------------------------------------------------------------------- sanity T1 / T2
srows = []


def add(check, variant, detail, fine, coarse, nodes_, viol):
    """one sanity row; groups_checked = all fine groups (singletons included, the 17:12 count); non-vacuity columns count
    coarse blocks with >= 2 genes that the fine partition actually splits."""
    nb, ns, nm = split_counts(fine, coarse, nodes_, MEM)
    srows.append({"check": check, "variant": variant, "detail": detail,
                  "groups_checked": len(groups({n: fine[n] for n in nodes_})), "violations": len(viol),
                  "coarse_blocks_ge2": nb, "coarse_blocks_split": ns, "member_holding_blocks_split": nm,
                  "example": ";".join(",".join(sorted(NAME[g] for g in S))[:200] for S in viol[:2])})


T1_VARIANTS = list(VARIANTS) + ["triangle", "17:03_tests_exact_triangle", "no_readthrough_nodes"]
for vname in T1_VARIANTS:
    for i in range(len(LEVELS)):
        for j in range(i + 1, len(LEVELS)):
            fine, coarse = LAB[(vname, LEVELS[j])], LAB[(vname, LEVELS[i])]
            viol = refines(fine, coarse)
            add("T1 nesting: components of G_j refine G_i", vname, f"{LEVELS[j]} in {LEVELS[i]}", fine, coarse, set(fine), viol)
# triangle inside plain components (truss(E_k) subset of E_k)
for tag, base in (("triangle", "primary"), ("17:03_tests_exact_triangle", "17:03_tests_exact")):
    for k in LEVELS:
        viol = refines(LAB[(tag, k)], LAB[(base, k)])
        add("triangle components refine plain components", f"{tag} vs {base}", k, LAB[(tag, k)], LAB[(base, k)],
            set(LAB[(tag, k)]), viol)

XSETS = (("any-overlap reads>=3", {g for g in V if reads.get(g, 0) >= 3}),
         ("any-overlap reads>=1", {g for g in V if reads.get(g, 0) >= 1}),
         ("unique reads>=3", {g for g in V if reads_u.get(g, 0) >= 3}))
erows = []


def expr_row(xname, vname, view, k, lab):
    G = groups(lab)
    mg = sorted([S for S in G.values() if S & MEM and len(S) >= 2], key=lambda S: (-len(S), sorted(S)))
    erows.append({"expression_set": xname, "variant": vname, "view": view, "level": k,
                  "member_components": ";".join(f"{len(S)}({len(S & MEM)})" for S in mg),
                  "members": " | ".join(",".join(sorted(NAME[g] for g in S & MEM)) for S in mg),
                  "nonmembers": " | ".join(",".join(sorted(NAME[g] for g in S - MEM)) for S in mg),
                  "member_singletons": ",".join(sorted(NAME[g] for S in G.values() if len(S) == 1 for g in S & MEM))})
    say(f"[EXPR {xname}] {vname} {view} {k}: member-holding expressed components (>= 2 genes): "
        + " | ".join(f"{len(S)}: members {sorted(NAME[g] for g in S & MEM)} + {len(S - MEM)} non-members" for S in mg))


for xname, X in XSETS:
    say(f"[T2] expression set X: {xname}: {len(X)} of {len(V)} nodes (members {len(X & MEM)})")
    for vname in ("primary", "with_same_locus", "L1=E1_as_built(D graph)", "17:03_tests_exact"):
        labX = {}
        for k in LEVELS:
            EX = [e for e in EDGES[vname][k] if e[0] in X and e[1] in X]
            labX[k] = components(X, EX)
            viol = refines(labX[k], LAB[(vname, k)], X)
            add(f"T2b: G_k[X] components refine G_k restricted to X ({xname})", vname, k, labX[k], LAB[(vname, k)], X, viol)
            if vname in ("primary", "17:03_tests_exact"):
                expr_row(xname, vname, "comp(G_k[X])", k, labX[k])
        for i in range(len(LEVELS) - 1):
            viol = refines(labX[LEVELS[i + 1]], labX[LEVELS[i]])
            add(f"T2a: G_(k+1)[X] components refine G_k[X] ({xname})", vname, f"{LEVELS[i + 1]} in {LEVELS[i]}",
                labX[LEVELS[i + 1]], labX[LEVELS[i]], X, viol)
    # triangle: two expression views, both nest (truss(G_k[X]) subset of truss(G_k)[X] subset of G_k[X])
    for base, tag in (("primary", "triangle"), ("17:03_tests_exact", "17:03_tests_exact_triangle")):
        labTX, labTX2 = {}, {}
        for k in LEVELS:
            EX = [e for e in EDGES[base][k] if e[0] in X and e[1] in X]
            keptX, _, _ = truss3(EX)
            labTX[k] = components(X, keptX)                                               # comp(truss(G_k[X]))
            labTX2[k] = components(X, [e for e in TRI[(tag, k)] if e[0] in X and e[1] in X])  # comp(truss(G_k)[X])
            viol = refines(labTX[k], LAB[(tag, k)], X)
            add(f"T2b: truss(G_k[X]) components refine truss(G_k) restricted to X ({xname})", tag, k, labTX[k],
                LAB[(tag, k)], X, viol)
            viol = refines(labTX[k], labTX2[k], X)
            add(f"truss views: comp(truss(G_k[X])) refines comp(truss(G_k)[X]) ({xname})", tag, k, labTX[k], labTX2[k], X,
                viol)
            viol = refines(labTX2[k], LAB[(tag, k)], X)
            add(f"T2b: comp(truss(G_k)[X]) refines truss(G_k) restricted to X ({xname})", tag, k, labTX2[k],
                LAB[(tag, k)], X, viol)
            expr_row(xname, base, "comp(truss(G_k[X]))", k, labTX[k])
            expr_row(xname, base, "comp(truss(G_k)[X])", k, labTX2[k])
        for i in range(len(LEVELS) - 1):
            viol = refines(labTX[LEVELS[i + 1]], labTX[LEVELS[i]])
            add(f"T2a: truss(G_(k+1)[X]) refine truss(G_k[X]) ({xname})", tag, f"{LEVELS[i + 1]} in {LEVELS[i]}",
                labTX[LEVELS[i + 1]], labTX[LEVELS[i]], X, viol)
            viol = refines(labTX2[LEVELS[i + 1]], labTX2[LEVELS[i]])
            add(f"T2a: truss(G_(k+1))[X] refine truss(G_k)[X] ({xname})", tag, f"{LEVELS[i + 1]} in {LEVELS[i]}",
                labTX2[LEVELS[i + 1]], labTX2[LEVELS[i]], X, viol)
    say(f"[T2] {xname} done; {time.time() - T0:.0f}s")
write(f"{OUT}/sanity.tsv", srows)
write(f"{OUT}/expr_views.tsv", erows)
say(f"[sanity] checks {len(srows)}; total violations {sum(r['violations'] for r in srows)}")

# ---------------------------------------------------------------------------------------------- chaining
INTEREST = ["PKD1", "PKD1P1", "PKD1P2", "PKD1P3", "PKD1P6", "PKD1P3-NPIPA1", "LOC131696449", "PKD1P4-NPIPA8",
            "PKD1P5-LOC105376752", "PKD1P6-NPIPP1", "PDXDC2P-NPIPB14P", "NPIPB1P", "NPIPB14P", "LOC100505915",
            "DHX40", "DHX40P1", "RNFT1", "RNFT1-DT", "RNFT1P3", "TBC1D3P1-DHX40P1", "TBC1D3P1", "USP6", "USP6NL", "TBC1D26",
            "TBC1D29P", "TBC1D28", "LOC100420408", "TBC1D3P5", "TBC1D3P7", "LOC124905656", "TBC1D3P6", "LOC100420289"]
by_name = collections.defaultdict(list)
for g in V:
    by_name[NAME[g]].append(g)
ANCHOR = {"NPIP": by_name["NPIPB2"][0], "TBC1D3": by_name["TBC1D3"][0]}


def side_of(name):
    return "NPIP" if any(x in name for x in ("PKD1", "NPIP", "PDXDC2P", "LOC131696449", "LOC100505915")) else "TBC1D3"


ADJ = {}
for _k in LEVELS:
    ADJ[_k] = collections.defaultdict(set)
    for _a, _b in EDGES["primary"][_k]:
        ADJ[_k][_a].add(_b)
        ADJ[_k][_b].add(_a)


def shortest_path(adj, src, dst):
    prev = {src: None}
    q = collections.deque([src])
    while q:
        x = q.popleft()
        if x == dst:
            break
        for y in sorted(adj[x]):
            if y not in prev:
                prev[y] = x
                q.append(y)
    if dst not in prev:
        return None
    path = [dst]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    return path[::-1]


def edge_desc(a, b):
    r = EDGE.get((a, b)) or EDGE.get((b, a))
    return (f"{NAME[a]}-{NAME[b]}[p {r['p_qualifies_6ko']} aa {fmt3(r['p_aa_identity'])}; c2 {r['d_c2_approx']} "
            f"(gb {r['d_c2_genebody']}, ex {r['d_c2_exon']}); sef {fmt3(r['d_shared_exon_frac'])}; w98 {fmt3(r['d_w98_gapexcl'])}; "
            f"pooled id {fmt3(r['d_e1_identity_gapexcl'])}; same_locus {r['same_locus']}]")


def fmt3(s):
    v = fnum(s)
    return "NA" if v is None else f"{v:.3f}"


crow = []
for nm in INTEREST:
    for g in by_name.get(nm, []):
        fam = side_of(nm)
        anc = ANCHOR[fam]
        row = {"gene": nm, "is_member": nodes[g]["is_member"], "catalog": nodes[g]["catalog"], "anchor": NAME[anc],
               "readthrough": nodes[g]["readthrough"], "n_reads_any": reads.get(g, 0)}
        for vname in ("primary", "triangle", "no_readthrough_nodes"):
            for k in LEVELS:
                lab = LAB[(vname, k)]
                if g not in lab:
                    row[f"{vname}|{k}"] = "removed"
                    continue
                row[f"{vname}|{k}"] = "with anchor" if lab[g] == lab[anc] else f"apart ({len(groups(lab)[lab[g]])})"
        for k in LEVELS:
            p = shortest_path(ADJ[k], anc, g) if LAB[("primary", k)][g] == LAB[("primary", k)][anc] else None
            row[f"path|{k}"] = " > ".join(edge_desc(p[i], p[i + 1]) for i in range(len(p) - 1)) if p and len(p) <= 8 else (
                f"path of {len(p) - 1} edges" if p else "")
        crow.append(row)
write(f"{OUT}/chaining.tsv", crow)

# ---------------------------------------------------------------------------------------------- triangle drops
trows = []
for k in LEVELS:
    G = groups(LAB[("primary", k)])
    for S in G.values():
        if not (S & MEM) or len(S) < 2:
            continue
        labT = LAB[("triangle", k)]
        parts = groups({g: labT[g] for g in S})
        big = max(parts.values(), key=lambda P: (len(P & MEM), len(P), sorted(P)))  # the part holding most members
        lost = S - big
        trows.append({"level": k, "group_rep": NAME[min(S)], "size": len(S), "members": ",".join(sorted(NAME[g] for g in S & MEM)),
                      "n_members": len(S & MEM), "triangle_parts": len(parts),
                      "triangle_singletons": sum(1 for P in parts.values() if len(P) == 1),
                      "main_part_size": len(big), "main_part_n_members": len(big & MEM),
                      "largest_part_size": max(len(P) for P in parts.values()),
                      "genes_outside_main_part": len(lost),
                      "members_outside_main_part": ",".join(sorted(NAME[g] for g in lost & MEM)),
                      "nonmembers_outside_main_part_first40": ",".join(sorted(NAME[g] for g in lost - MEM)[:40]),
                      "is_2copy_group": len(S) == 2})
write(f"{OUT}/triangle_drops.tsv", trows)
with open(f"{OUT}/levels.out", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
say(f"[done] {time.time() - T0:.0f}s")
