#!/usr/bin/env python3
"""NPIP/TBC1D3 layer order (integration, revised after the 2026-09-16 audit): containment, tournament, enforcement cost,
truth agreement, EXPR operator, disagreement lists — layers P (protein), D (= the RefSeq E1 guided catalog, §6kl; NOT
§0★★ clause 2/4), C (subfamily clades) on the CORRECTED tables (bench/layer_order/lo_corrected_tables.py).
Spec: docs/superpowers/specs/2026-09-16-family-layer-order-design.md (SCOPE AMENDMENT + Definitions, incl. the
2026-09-16 audit amendments). S1/S2/S3 are deferred (user scope 2026-09-16 15:03) and not read.

Conventions (all stated in bench/LAYER_ORDER_NPIP_TBC1D3.md):
  layer universe  P: U gene with a §6ko protein (work/P/proteins.index.tsv), including coding genes P does not place with a
                  member ('P|other:<gene>' = its own group). D: U gene that is a node of one of the two E1 catalogs.
                  C layers: U gene that is a leaf of the §6js reference trees (literature C: 31 leaves; C_tree: the 30 leaves
                  common to the exon and intron trees). Outside its universe a gene is 'not in layer', never a singleton.
  C layers        literature-anchored (CIRCULAR reference): C_L1, C_mid (:= as-built C_mid ∨ C_fine), C_mid_as_built,
                  C_fine. Clause-5 (§0★★) split-system layers: all SH-aLRT > 75 splits of either tree restricted to the
                  common leaves, smaller side = cluster, kept iff compatible with every other kept split (Buneman).
                  Ctree_top = maximal clusters (the partition whose pairs are 'share >= 1 subfamily'); Ctree_min = minimal.
  c(X ⊇ Y)        |pairs(Y) ∩ pairs(X)| / |pairs(Y)| over genes in U_X ∩ U_Y (∩ family); 0 pairs -> vacuous (verdict NA).
  nesting         fraction of Y groups (>= 2 genes after restriction to U_X ∩ U_Y) inside one X group.
  JOIN(P over D)  components of the graph whose vertices are D groups (+ P-universe genes outside D as own vertices), two
                  vertices joined when two of their genes share a P group. Two variants: 'U' (D groups cut to U, as first
                  reported) and 'whole' (whole catalog groups; P groups for genes outside U from the k = 2 MCL on N_2,
                  integrate_slim/P_N2_clusters.tsv; genes outside N_2 are their own P group), U re-closed afterwards.
  REFINE(X in Y)  X's grouping recomputed inside each Y group: P -> mcl_port.mcl (inflation 2.8) on P.edges induced on
                  X group ∩ Y group; C -> clade co-membership ∩ Y group. A gene outside Y's universe is a FREE CHOICE,
                  variants reported: 'attach' (P: to the refined part with the largest summed X-edge weight; C: kept
                  with the largest part of its clade), 'single' (P: its own group; C: isolated), and for C 'together'
                  (the clade's outside genes form one part of their own).
  EXPR(L)         within each L group, connected components (>= 2 genes) of the subgraph induced on expressed genes.
                  Expressed = reads >= t, t = 3 in the main tables (sweep t = 1..5). Read modes: any (PRIMARY, spec rule),
                  unique (the read hits exons of one RefSeq record genome-wide), unique_mr (same, ignoring the 6 readthrough
                  records that overlap a member on the same strand). Edges: P -> P.edges; D -> D.edges.corrected; C ->
                  clade co-membership; P_join -> P.edges ∪ D.edges; P_ref -> P.edges. EXPR(L) ⊆ L holds by construction;
                  with co-membership edges no L group can split, so T2 holds trivially there (reported as guaranteed).
  truth scores    (A) on U (recall CONDITIONED ON THE PREDICTION: U is the closure of the scored layers);
                  (B) member-anchored, layer-independent: pairs with >= 1 member endpoint over every gene of a truth group
                  that holds a member (HGNC gene group; Soto family), restricted to the layer's universe genome-wide
                  (P: genes with a protein; D: catalog nodes). Bipartite F only when both sides have >= 1 pair.
"""
import collections
import csv
import itertools
import sys
from fractions import Fraction

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench")
sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench/layer_order")
sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light/scripts")
import guided_pipeline as gp  # noqa: E402
import mcl_port  # noqa: E402

LIGHT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light"
INT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/integrate_slim"
H = "/mnt/linuxdisk/home/juanfraitu/o1_falsemerge"
HGNC = "/mnt/linuxdisk/home/juanfraitu/winloci_data/hgnc/hgnc_complete_set.txt"
PAF = {"c15_17_22": f"{H}/human2/genes.asm20.paf", "c16_19_20": f"{H}/lit/aj_ho/refseq/all.paf"}
NODES_C16 = f"{H}/lit/aj_ho/refseq/nodes.tsv"
FAMS = ("NPIP", "TBC1D3", "pooled")
READTHROUGH_OVER_MEMBERS = ("PKD1P3-NPIPA1", "LOC131696449", "PKD1P4-NPIPA8", "PKD1P5-LOC105376752", "PDXDC2P-NPIPB14P",
                            "TBC1D3P1-DHX40P1")


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def write(path, rows, cols=None):
    cols = cols or (list(rows[0].keys()) if rows else ["empty"])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def fmt(x):
    if x is None:
        return "NA"
    if isinstance(x, Fraction):
        return f"{float(x):.3f}"
    if isinstance(x, float):
        return "NA" if x != x else f"{x:.3f}"
    return str(x)


# ---------------------------------------------------------------------------------------------------------- data
U = {r["gene_id"]: r for r in tsv(f"{LIGHT}/universe.corrected.tsv")}
NAME = {g: r["name"] for g, r in U.items()}
MEMB = {g for g, r in U.items() if r["is_member"] == "yes"}
SIDE = {g: r["family_side"] for g, r in U.items()}


def fam_genes(fam, side=None):
    side = side or SIDE
    return set(side) if fam == "pooled" else {g for g in side if side[g] == fam}


def layer_from(col, uni_col):
    return {g: r[col] for g, r in U.items() if r[uni_col] == "yes"}


def _join_labels(a, b):
    """partition join of two labelings on the same genes (union-find over shared labels)."""
    parent = {g: g for g in a}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for lab in (a, b):
        first = {}
        for g in sorted(a):
            if lab[g] in first:
                parent[find(g)] = find(first[lab[g]])
            else:
                first[lab[g]] = g
    comp = collections.defaultdict(set)
    for g in a:
        comp[find(g)].add(g)
    out = {}
    for G in comp.values():
        names = sorted({lab[g] for lab in (a, b) for g in G if "singleton" not in lab[g] and "|" in lab[g]})
        if len(G) == 1:
            out[next(iter(G))] = a[next(iter(G))]
            continue
        for g in G:
            out[g] = "Cmid|" + ("+".join(x.split("|", 1)[1] for x in names) if names else "+".join(sorted(G)))
    return out


LAYERS = {"P": layer_from("P_group", "in_P_universe"), "D": layer_from("D_group", "in_D_universe"),
          "C_L1": layer_from("C_L1", "in_C_universe"), "C_mid": layer_from("C_mid", "in_C_universe"),
          "C_fine": layer_from("C_fine", "in_C_universe")}
C_MID_AS_BUILT = dict(LAYERS["C_mid"])
# C_mid as built holds ONE recovered group (the named NPIPB subfamily) and singletons elsewhere, so it is not a level
# between C_L1 and C_fine (A6-9 and B6-9 are singletons in it). Hierarchical mid level = C_mid ∨ C_fine (by construction
# c(C_mid ⊇ C_fine) = 1). The as-built layer is kept as C_mid_ab in every table.
LAYERS["C_mid"] = _join_labels(C_MID_AS_BUILT, LAYERS["C_fine"])
LAYERS["C_mid_ab"] = C_MID_AS_BUILT

P_EDGES = {}
for r in tsv(f"{LIGHT}/P.edges.tsv"):
    P_EDGES[frozenset((r["u_gene_id"], r["v_gene_id"]))] = (float(r["weight"]), float(r["identity"]),
                                                             float(r["coverage_longer"]))
P_W = {k: v[0] for k, v in P_EDGES.items()}
D_EDGE_ROWS = tsv(f"{LIGHT}/D.edges.corrected.tsv")
D_EDGES = {}
for r in D_EDGE_ROWS:
    if "?" in (r["u_gene_id"], r["v_gene_id"]) or r["u_gene_id"] == r["v_gene_id"]:
        continue
    k = frozenset((r["u_gene_id"], r["v_gene_id"]))
    D_EDGES[k] = max(D_EDGES.get(k, 0.0), float(r["weight"]))
_rec = {r["gene_id"]: r for r in tsv(f"{INT}/expr_recount.tsv")}
READS = {g: {"unique": int(r["n_reads_unique"]), "any": int(r["n_reads_any"]),
             "unique_mr": int(_rec[NAME[g]]["n_reads_unique_mr"])} for g, r in U.items()}
MODES = ("any", "unique", "unique_mr")


def groups_of(lab, genes):
    out = collections.defaultdict(set)
    for g in genes:
        out[lab[g]].add(g)
    return out


def pairs_of(lab, genes):
    s = set()
    for G in groups_of(lab, genes).values():
        for a, b in itertools.combinations(sorted(G), 2):
            s.add((a, b))
    return s


# ---------------------------------------------------------------------------------------------------------- catalogs
def catalog_all():
    """gene_id -> (catalog, cluster_id or '', record key, folded representative or '') for every gene of both E1 catalogs
    (lo_corrected_tables.catalog_keys/membership), plus gene rows."""
    import lo_corrected_tables as CT
    genes = {r["gene_id"]: r for r in tsv(f"{LIGHT}/work/refseq/genes.tsv")}
    by_coord = collections.defaultdict(list)
    for g in genes.values():
        by_coord[(g["chrom"], int(g["start0"]) + 1, int(g["end"]))].append(g)
    out = {}
    for tag, c in CT.CATALOGS.items():
        k2g = CT.catalog_keys(by_coord, c["kind"], c["path"])
        for g, (cl, k, f) in CT.membership(k2g, c["D"]).items():
            out[g] = (tag, cl, k, f)
    return out, genes


def soto_ok(db, exons, genes, g):
    """Soto families of RefSeq gene g under the U table's flag rule (lo_corrected_tables): 'ok' iff matched, match quality
    not weak, not ambiguous. Returns the family string or None."""
    import soto_map
    r = genes[g]
    m = soto_map.map_gene(db, r["name"], r["chrom"], r["strand"], exons.get(g) or [(int(r["start0"]), int(r["end"]))])
    ok = m["soto_gene_id"] and m["soto_match_quality"] != "weak" and m["soto_ambiguous"] != "yes"
    return m["soto_families"] if ok and m["soto_families"] else None


def hgnc_all(genes):
    by_id, by_sym = {}, {}
    for r in tsv(HGNC):
        by_id[r["hgnc_id"]] = r
        by_sym[r["symbol"]] = r
    dbx = dict(line.rstrip("\n").split("\t") for line in open(f"{LIGHT}/work/refseq/gene_dbxref.tsv"))
    out = {}
    for g, r in genes.items():
        ids = [x[5:] for x in dbx.get(g, "").split(",") if x.startswith("HGNC:")]
        h = by_id[ids[0]] if ids and ids[0] in by_id else by_sym.get(r["name"])  # = lo_corrected_tables (dbxref, else symbol)
        if h and h["gene_group_id"]:
            out[g] = h["gene_group_id"]
    return out


# ---------------------------------------------------------------------------------------------------------- C_tree
def c_tree():
    """§0★★ clause 5 estimator on the §6js reference trees: supported (SH-aLRT > 75, either tree), pairwise-compatible
    split system on the leaves common to both trees; cluster = smaller side of each split."""
    rows = tsv(f"{LIGHT}/C.supported_clades.tsv")
    trees = collections.defaultdict(list)
    leaves = collections.defaultdict(set)
    for r in rows:
        s, c = set(r["side"].split(",")), set(r["complement"].split(","))
        trees[(r["family"], r["tree"])].append((s, c, r["sh_alrt"], r["ufboot"]))
        leaves[(r["family"], r["tree"])] |= s | c
    name2gid = {NAME[g]: g for g in U}
    lab_min, lab_top, lab_root, cl_rows, lit_rows, clusters = {}, {}, {}, [], [], {}
    lit = {r["name"]: r for r in tsv(f"{LIGHT}/truth_literature_subfamilies.corrected.tsv")}
    for fam in ("NPIP", "TBC1D3"):
        L = leaves[(fam, "exon")] & leaves[(fam, "intron")]
        splits = {}
        for tr in ("exon", "intron"):
            for s, c, sh, bs in trees[(fam, tr)]:
                s2, c2 = frozenset(s & L), frozenset(c & L)
                if min(len(s2), len(c2)) < 2:
                    continue
                small = min((s2, c2), key=lambda x: (len(x), sorted(x)))
                splits.setdefault(small, []).append(f"{tr} {sh}/{bs}")

        def compat(a, b):  # unrooted splits a|L-a and b|L-b: compatible iff one of the four intersections is empty
            A2, B2 = L - a, L - b
            return not (a & b) or not (a & B2) or not (A2 & b) or not (A2 & B2)
        K = [s for s in splits if all(compat(s, t) for t in splits if t != s)]
        conflicts = {s: sorted(",".join(sorted(t)) for t in splits if t != s and not compat(s, t)) for s in splits}
        minimal = [s for s in K if not any(t < s for t in K)]
        maximal = [s for s in K if not any(s < t for t in K)]
        clusters[fam] = (L, K)
        for s in minimal:
            for n in s:
                lab_min[name2gid[n]] = f"Ctree_min|{fam}|{'+'.join(sorted(s))}"
        for s in maximal:
            for n in s:
                lab_top[name2gid[n]] = f"Ctree_top|{fam}|{'+'.join(sorted(s))}"
        # rooted variant: root at the most balanced compatible split; its two sides are the top clusters
        root = max(K, key=lambda s: (min(len(s), len(L) - len(s)), sorted(s)))
        for side in (root, L - root):
            for n in side:
                lab_root[name2gid[n]] = f"Ctree_root|{fam}|{'+'.join(sorted(side))}"
        for n in L:
            lab_min.setdefault(name2gid[n], f"Ctree|singleton:{name2gid[n]}")
            lab_top.setdefault(name2gid[n], f"Ctree|singleton:{name2gid[n]}")
        for s in sorted(splits, key=lambda x: (len(x), sorted(x))):
            cl_rows.append({"family": fam, "cluster_smaller_side": ",".join(sorted(s)), "size": len(s),
                            "support": "; ".join(splits[s]), "compatible_with_all": "yes" if s in K else "no",
                            "minimal": "yes" if s in minimal else "no", "maximal": "yes" if s in maximal else "no",
                            "conflicts_with": " | ".join(conflicts[s])})
        for lvl in ("level1", "level2", "npipb_named_subfamily"):
            grp = collections.defaultdict(set)
            for n in L:
                if lit.get(n) and lit[n][lvl] and lit[n][lvl] != "no":
                    grp["named NPIPB {B3,B4,B5,B11,B12,B13}" if lvl == "npipb_named_subfamily" else lit[n][lvl]].add(n)
            for gname, G in sorted(grp.items()):
                if len(G) < 2:
                    continue
                Gf, comp = frozenset(G), frozenset(L - G)
                key = Gf if Gf in splits else comp if comp in splits else None
                lit_rows.append({"family": fam, "level": lvl, "literature_group": gname, "genes_on_common_leaves": len(G),
                                 "supported_split_any_tree": "yes" if key else "no",
                                 "support": "; ".join(splits.get(Gf, []) + splits.get(comp, [])),
                                 "in_compatible_system": "yes" if (Gf in K or comp in K) else "no",
                                 "conflicts_with": " | ".join(conflicts[key]) if key else ""})
    return {"Ctree_min": lab_min, "Ctree_top": lab_top, "Ctree_root": lab_root}, cl_rows, lit_rows, clusters


CT_LAB, CT_ROWS, CT_LIT, CT_CLUSTERS = c_tree()
LAYERS.update({"Ctree_top": CT_LAB["Ctree_top"], "Ctree_min": CT_LAB["Ctree_min"]})
LAYERS_VARIANT_C = {"Ctree_root": CT_LAB["Ctree_root"]}


# ---------------------------------------------------------------------------------------------------------- containment
def containment(X, Y, fam, layers, side=None):
    LX, LY = layers[X], layers[Y]
    S = set(LX) & set(LY) & fam_genes(fam, side)
    pY, pX = pairs_of(LY, S), pairs_of(LX, S)
    c = Fraction(len(pY & pX), len(pY)) if pY else None
    gY = [G for G in groups_of(LY, S).values() if len(G) >= 2]
    inside = sum(1 for G in gY if len({LX[g] for g in G}) == 1)
    return {"family": fam, "X": X, "Y": Y, "n_genes_in_both": len(S), "pairs_Y": len(pY), "pairs_Y_in_X": len(pY & pX),
            "c_X_contains_Y": "vacuous (0 pairs)" if c is None else fmt(c), "_c": c,
            "groups_Y_ge2": len(gY), "groups_Y_inside_one_X": inside,
            "nesting": "NA" if not gY else fmt(Fraction(inside, len(gY)))}


def verdict(a, b):
    """X above Y iff c(X ⊇ Y) > c(Y ⊇ X); NA when either containment is over 0 pairs."""
    if a["_c"] is None or b["_c"] is None:
        return "NA (0 pairs on one side)"
    if a["_c"] > b["_c"]:
        return f"{a['X']} above {a['Y']}"
    if b["_c"] > a["_c"]:
        return f"{a['Y']} above {a['X']}"
    return "tie"


def tournament(names, fam, layers, side=None):
    rows = []
    for X, Y in itertools.combinations(names, 2):
        a, b = containment(X, Y, fam, layers, side), containment(Y, X, fam, layers, side)
        rows.append({"family": fam, "X": X, "Y": Y, "c(X>=Y)": a["c_X_contains_Y"], "n_pairs_Y": a["pairs_Y"],
                     "groups_Y_inside_X": f"{a['groups_Y_inside_one_X']}/{a['groups_Y_ge2']}",
                     "c(Y>=X)": b["c_X_contains_Y"], "n_pairs_X": b["pairs_Y"],
                     "groups_X_inside_Y": f"{b['groups_Y_inside_one_X']}/{b['groups_Y_ge2']}",
                     "n_genes": a["n_genes_in_both"], "verdict": verdict(a, b)})
    return rows


# ---------------------------------------------------------------------------------------------------------- operators
def uf_components(nodes, links):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in links:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = collections.defaultdict(set)
    for n in nodes:
        comp[find(n)].add(n)
    return list(comp.values())


def join_P_over_D(LP, LD):
    genes = set(LP) | set(LD)
    vert = {g: ("D", LD[g]) if g in LD else ("Ponly", g) for g in genes}
    links = []
    for G in groups_of(LP, set(LP)).values():
        G = sorted(G)
        links += [(vert[G[0]], vert[h]) for h in G[1:]]
    comps = uf_components(set(vert.values()), links)
    lab = {}
    for comp in comps:
        members = sorted(g for g in genes if vert[g] in comp)
        name = "Pjoin|" + ("+".join(sorted({v[1] for v in comp if v[0] == "D"})) or members[0])
        for g in members:
            lab[g] = name
    return lab


def join_P_over_D_whole(cat, n2, pidx_genes):
    """JOIN on whole catalog groups. cat: gene -> (catalog, cluster, key, folded) for all catalog genes; n2: gene -> k=2 MCL
    cluster for N_2 genes (member clusters = PC1/PC2); pidx_genes: genes with a protein. Returns labels for every gene of
    a member-containing component, and the P label used per gene."""
    Pl = {}
    for g in pidx_genes:
        Pl[g] = n2.get(g, f"P|own:{g}")
    vert = {}
    for g in set(cat) | set(Pl):
        if g in cat and cat[g][1]:
            vert[g] = ("D", f"D|{cat[g][0]}|{cat[g][1]}")
        elif g in cat:
            vert[g] = ("Dsingle", g)
        else:
            vert[g] = ("Ponly", g)
    links = []
    for G in groups_of(Pl, set(Pl)).values():
        if len(G) < 2:
            continue
        G = sorted(G)
        links += [(vert[G[0]], vert[h]) for h in G[1:]]
    comps = uf_components(set(vert.values()), links)
    by_vert = collections.defaultdict(list)
    for g, v in vert.items():
        by_vert[v].append(g)
    lab = {}
    for comp in comps:
        genes_c = sorted(g for v in comp for g in by_vert[v])
        if not set(genes_c) & MEMB:
            continue
        name = "PjoinW|" + "+".join(sorted(v[1] for v in comp if v[0] == "D"))
        for g in genes_c:
            lab[g] = name
    return lab, Pl


def refine(LX, LY, edge_w, use_mcl, outside="attach"):
    """X's grouping recomputed inside each Y group (see module docstring). outside: 'attach' | 'single' | 'together'
    (the X group's genes outside Y's universe form one part of their own)."""
    lab = {}
    for xg, G in groups_of(LX, set(LX)).items():
        inside = collections.defaultdict(set)
        out_genes = []
        for g in G:
            if g in LY:
                inside[LY[g]].add(g)
            else:
                out_genes.append(g)
        parts = []
        for yg, part in inside.items():
            if len(part) == 1 or not use_mcl:
                parts.append(set(part))
                continue
            E = {tuple(sorted(k)): w for k, w in edge_w.items() if k <= part}
            cl = mcl_port.mcl(E) if E else []
            got = set().union(*map(set, cl)) if cl else set()
            parts += [set(c) for c in cl] + [{g} for g in part - got]
        if outside == "together" and out_genes:
            parts.append(set(out_genes))
            out_genes = []
        for g in sorted(out_genes):
            if outside == "single":
                parts.append({g})
            elif use_mcl:
                best = max(parts, key=lambda p: (sum(edge_w.get(frozenset((g, h)), 0.0) for h in p), len(p)), default=None)
                if best is not None and sum(edge_w.get(frozenset((g, h)), 0.0) for h in best) > 0:
                    best.add(g)
                else:
                    parts.append({g})
            else:
                best = max(parts, key=len, default=None)
                if best is not None:
                    best.add(g)
                else:
                    parts.append({g})
        for i, p in enumerate(sorted(parts, key=lambda p: sorted(p))):
            for g in p:
                lab[g] = f"{xg}#r{i}" if len(parts) > 1 else xg
    return lab


def expr_layer(lab, expressed, edges):
    """EXPR(L): {L group: [components >= 2]}, plus dropped expressed singletons per group."""
    out, dropped = {}, {}
    for grp, G in groups_of(lab, set(lab)).items():
        Ex = G & expressed
        links = [tuple(k) for k in edges(Ex)]
        comps = [c for c in uf_components(Ex, links)]
        out[grp] = [c for c in comps if len(c) >= 2]
        dropped[grp] = sorted(g for c in comps if len(c) == 1 for g in c)
    return out, dropped


def edge_set(kind):
    return {"P": [P_EDGES], "P_ref": [P_EDGES], "D": [D_EDGES], "P_join": [P_EDGES, D_EDGES]}.get(kind)


def edge_fn(kind):
    pools = edge_set(kind)
    if kind == "comembership" or pools is None:  # C layers: clade co-membership
        return lambda S: [frozenset(p) for p in itertools.combinations(sorted(S), 2)]
    return lambda S: [k for pool in pools for k in pool if k <= S]


def expr_labels(comps_by_group):
    lab = {}
    for grp, comps in comps_by_group.items():
        for i, c in enumerate(sorted(comps, key=lambda c: sorted(c))):
            for g in c:
                lab[g] = f"{grp}#e{i}"
    return lab


# ---------------------------------------------------------------------------------------------------------- truth scores
def bip_jaccard(pred, true):
    P, T = sorted(set(pred), key=str), sorted(set(true), key=str)
    M = np.zeros((len(T), len(P)), dtype=int)
    for p, t in zip(pred, true):
        M[T.index(t), P.index(p)] += 1
    ts, ps = M.sum(axis=1), M.sum(axis=0)
    J = np.zeros(M.shape)
    for i in range(len(T)):
        for j in range(len(P)):
            if M[i, j]:
                J[i, j] = M[i, j] / (ts[i] + ps[j] - M[i, j])
    r, c = linear_sum_assignment(-J)
    matched = sum(M[i, j] for i, j in zip(r, c) if M[i, j] > 0)
    msize = sum(ps[j] for i, j in zip(r, c) if M[i, j] > 0)
    return matched / len(pred), (matched / msize if msize else float("nan"))


def f1(r, p):
    return 2 * r * p / (r + p) if r + p and r == r and p == p else float("nan")


def score(lab, truth, genes):
    gs = sorted(g for g in genes if g in lab and g in truth)
    if not gs:
        return {"n_genes": 0}
    pred = [lab[g] for g in gs]
    true = [truth[g] for g in gs]
    tp = sum(1 for i, j in itertools.combinations(range(len(gs)), 2) if pred[i] == pred[j] and true[i] == true[j])
    npp = sum(1 for i, j in itertools.combinations(range(len(gs)), 2) if pred[i] == pred[j])
    ntp = sum(1 for i, j in itertools.combinations(range(len(gs)), 2) if true[i] == true[j])
    out = {"n_genes": len(gs), "truth_pairs": ntp, "pred_pairs": npp, "tp_pairs": tp,
           "pair_precision": fmt(tp / npp) if npp else "NA", "pair_recall": fmt(tp / ntp) if ntp else "NA"}
    if ntp == 0 or npp == 0:  # audit: no bipartite F when either side has no pairs (singleton-to-singleton matches)
        out.update({"bip_F_count(§6ks)": "NA (a side has 0 pairs)", "bip_R": "NA", "bip_P": "NA", "bip_F_jaccard": "NA"})
        return out
    br, bp = gp.bipartite(pred, true)
    jr, jp = bip_jaccard(pred, true)
    out.update({"bip_F_count(§6ks)": fmt(f1(br, bp)), "bip_R": fmt(br), "bip_P": fmt(bp), "bip_F_jaccard": fmt(f1(jr, jp))})
    return out


def score_member_anchored(lab, truth_parts, genes, members):
    """pairs (a, b), a != b, both in genes, >= 1 in members. pred: same lab (a gene absent from lab is its own group);
    truth: share >= 1 truth part (truth_parts: gene -> frozenset of group ids)."""
    gs = sorted(g for g in genes if truth_parts.get(g))
    tp = npp = ntp = 0
    for a, b in itertools.combinations(gs, 2):
        if a not in members and b not in members:
            continue
        pr = a in lab and b in lab and lab[a] == lab[b]
        tr = bool(truth_parts[a] & truth_parts[b])
        tp += pr and tr
        npp += pr
        ntp += tr
    return {"n_genes": len(gs), "n_members": len(set(gs) & members), "truth_pairs": ntp, "pred_pairs": npp,
            "tp_pairs": tp, "pair_precision": fmt(tp / npp) if npp else "NA", "pair_recall": fmt(tp / ntp) if ntp else "NA"}


# ---------------------------------------------------------------------------------------------------------- DNA pair stats
def dna_pair_stats(pairs):
    """identity / coverage for D edges (catalog, key_u, key_v) from the catalog PAF, mcl_families graph rule
    (records >= 300 bp, identity >= 0.70, pooled identity, union of aligned intervals on the longer gene / its exon-union
    length). Exon-union lengths: c16_19_20 nodes.tsv; c15_17_22 light/work/refseq/exons.tsv (checked via weight)."""
    exlen = {}
    for r in tsv(NODES_C16):
        exlen[("c16_19_20", f"{r['chrom']}:{int(r['start']) + 1}-{r['end']}")] = sum(
            int(b) - int(a) for a, b in (x.split("-") for x in r["exons"].split(",")))
    for r in tsv(f"{LIGHT}/work/refseq/exons.tsv"):
        exlen.setdefault(("c15_17_22", f"{r['chrom']}:{int(r['start0']) + 1}-{r['end']}"),
                         sum(int(b) - int(a) for a, b in (x.split("-") for x in r["exons"].split(","))))
    want = collections.defaultdict(set)
    for cat, a, b in pairs:
        want[cat].add(frozenset((a, b)))
    acc = {}
    for cat, S in want.items():
        keys = set().union(*S)
        with open(PAF[cat]) as fh:
            for line in fh:
                f = line.split("\t", 12)
                q, t = f[0], f[5]
                if q == t or q not in keys or t not in keys or frozenset((q, t)) not in S:
                    continue
                nm, bl = int(f[9]), int(f[10])
                if bl < 300 or nm / max(bl, 1) < 0.70:
                    continue
                k = (cat, frozenset((q, t)))
                e = acc.setdefault(k, {"nm": 0, "bl": 0, "iv": collections.defaultdict(list)})
                e["nm"] += nm
                e["bl"] += bl
                e["iv"][q].append((int(f[2]), int(f[3])))
                e["iv"][t].append((int(f[7]), int(f[8])))
    out = {}
    for (cat, pr), e in acc.items():
        a, b = sorted(pr)
        la, lb = exlen.get((cat, a)), exlen.get((cat, b))
        if la is None or lb is None:
            continue
        longer = a if la >= lb else b
        merged = gp.merge(sorted(e["iv"][longer]))
        cov = min(1.0, sum(y - x for x, y in merged) / max(la if longer == a else lb, 1))
        out[(cat, pr)] = (e["nm"] / e["bl"], cov)
    return out


# ---------------------------------------------------------------------------------------------------------- main
def main():
    log = []

    def say(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        log.append(s)

    LP, LD = LAYERS["P"], LAYERS["D"]
    say(f"[universe] U {len(U)}; members {len(MEMB)}; universes P {len(LP)} D {len(LD)} C_lit {len(LAYERS['C_L1'])} "
        f"C_tree {len(LAYERS['Ctree_top'])}")
    tab_p = {r["gene_id"]: r["group_id"] for r in tsv(f"{LIGHT}/P.groups.corrected.tsv")}
    tab_d = {r["gene_id"]: r["group_id"] for r in tsv(f"{LIGHT}/D.groups.corrected.tsv") if r["group_id"] != "D|NA"}
    say(f"[check] P.groups.corrected == universe P labels: {tab_p == LP}; D.groups.corrected == universe D labels: "
        f"{tab_d == LD}")

    # ------------------------------------------------ containment + tournament
    main_layers = ["P", "D", "Ctree_top", "Ctree_min"]
    ref_layers = ["C_L1", "C_mid", "C_mid_ab", "C_fine"]
    allL = main_layers + ref_layers
    crow, trow = [], []
    for fam in FAMS:
        for X in allL:
            for Y in allL:
                if X != Y:
                    r = containment(X, Y, fam, LAYERS)
                    r.pop("_c")
                    crow.append(r)
        trow += tournament(allL, fam, LAYERS)
    write(f"{INT}/containment.tsv", crow)
    write(f"{INT}/tournament.tsv", trow)
    for r in trow:
        say(f"[tournament] {r['family']:6s} {r['X']:>9s} vs {r['Y']:<9s}: c(X⊇Y) {r['c(X>=Y)']} [{r['n_pairs_Y']}] "
            f"{r['groups_Y_inside_X']} | c(Y⊇X) {r['c(Y>=X)']} [{r['n_pairs_X']}] {r['groups_X_inside_Y']} | genes "
            f"{r['n_genes']} -> {r['verdict']}")
    # one common gene set for every main layer
    common = set.intersection(*(set(LAYERS[x]) for x in main_layers + ["C_L1"]))
    lay_c = {x: {g: LAYERS[x][g] for g in common} for x in allL}
    crow2 = []
    for fam in FAMS:
        crow2 += tournament(allL, fam, lay_c)
    write(f"{INT}/tournament_common_genes.tsv", crow2)
    say(f"[common] genes in every layer universe (P ∩ D ∩ C_tree ∩ C_lit): {len(common)} "
        f"(NPIP {len(common & fam_genes('NPIP'))}, TBC1D3 {len(common & fam_genes('TBC1D3'))})")
    for r in crow2:
        if r["family"] == "pooled":
            say(f"[common] pooled {r['X']:>9s} vs {r['Y']:<9s}: {r['c(X>=Y)']} [{r['n_pairs_Y']}] {r['groups_Y_inside_X']} | "
                f"{r['c(Y>=X)']} [{r['n_pairs_X']}] {r['groups_X_inside_Y']} -> {r['verdict']}")
    # all compatible clusters (not only the partition levels) inside one P / D group
    for fam, (L, K) in CT_CLUSTERS.items():
        for X in ("P", "D"):
            gid = {NAME[g]: g for g in LAYERS[X]}
            ok = [s for s in K if all(n in gid for n in s)]
            inside = sum(1 for s in ok if len({LAYERS[X][gid[n]] for n in s}) == 1)
            say(f"[C_tree hierarchy] {fam}: compatible supported clusters {len(K)}; inside one {X} group {inside}/{len(ok)} "
                f"(clusters with every leaf in {X}'s universe)")
    write(f"{INT}/ctree_clusters.tsv", CT_ROWS)
    write(f"{INT}/ctree_literature_groups.tsv", CT_LIT)

    # ------------------------------------------------ enforcement
    cat, genes = catalog_all()
    pidx_genes = {r["gene_id"] for r in tsv(f"{LIGHT}/work/P/proteins.index.tsv")}
    n2 = {}
    for r in tsv(f"{INT}/P_N2_clusters.tsv"):
        n2[r["gene_id"]] = r["k2_cluster"]
    # the k = 2 MCL member clusters must equal P's groups
    for grp, G in groups_of(LP, {g for g in LP if not LP[g].startswith("P|other")}).items():
        assert len({n2[g] for g in G}) == 1 and sum(1 for x in n2.values() if x == n2[next(iter(G))]) == len(G), grp
    Pj = join_P_over_D(LP, LD)
    Pjw, Pl_w = join_P_over_D_whole(cat, n2, pidx_genes)
    Pr = refine(LP, LD, P_W, use_mcl=True, outside="attach")
    Prs = refine(LP, LD, P_W, use_mcl=True, outside="single")
    enf = dict(LAYERS)
    enf.update({"P_join": Pj, "P_join_whole": Pjw, "P_ref": Pr, "P_ref_single": Prs})
    cl_names = ("C_L1", "C_mid", "C_fine", "Ctree_top", "Ctree_min")
    for cl in cl_names:
        for Yn, LY in (("D", LD), ("P", LP)):
            enf[f"{cl}_ref{Yn}"] = refine(LAYERS[cl], LY, {}, use_mcl=False, outside="attach")
            enf[f"{cl}_ref{Yn}_single"] = refine(LAYERS[cl], LY, {}, use_mcl=False, outside="single")
            enf[f"{cl}_ref{Yn}_together"] = refine(LAYERS[cl], LY, {}, use_mcl=False, outside="together")
    extra_w = sorted(genes[g]["name"] for g in set(Pjw) - set(U))
    say(f"[JOIN whole] member-containing components: " + " | ".join(
        f"{k}: {len(v)} genes" for k, v in sorted(groups_of(Pjw, set(Pjw)).items())) + f"; genes outside U {len(extra_w)}: "
        f"{extra_w}")
    for g in sorted(set(Pjw) - set(U), key=lambda x: genes[x]["name"]):
        say(f"   {genes[g]['name']} biotype {genes[g]['biotype']} catalog group {cat[g][1] if g in cat else '-'} P label "
            f"{Pl_w.get(g, 'no protein')}")
    changes = []
    plan = [("P_join", Pj, LP), ("P_join_whole", Pjw, LP), ("P_ref", Pr, LP), ("P_ref_single", Prs, LP)] + [
        (f"{c}_{s}", enf[f"{c}_{s}"], LAYERS[c]) for c in cl_names for s in ("refD", "refD_single", "refD_together", "refP", "refP_single",
                                                                "refP_together")]
    for name, lab, ref in plan:
        S = set(ref)
        before, after = pairs_of(ref, S), pairs_of(lab, S & set(lab))
        added = sorted((NAME[a], NAME[b]) for a, b in after - before)
        removed = sorted((NAME[a], NAME[b]) for a, b in before - after)
        extra = sorted(genes[g]["name"] for g in set(lab) - S)
        changes.append({"enforced": name, "pairs_on_original_universe_before": len(before), "after": len(after),
                        "pairs_added": len(added), "pairs_removed": len(removed),
                        "genes_added_to_universe": len(extra),
                        "added_examples": "; ".join(f"{a}-{b}" for a, b in added[:12]),
                        "removed_examples": "; ".join(f"{a}-{b}" for a, b in removed[:12])})
        say(f"[enforce] {name}: pairs on original universe {len(before)} -> {len(after)} (+{len(added)} / -{len(removed)}); "
            f"universe +{len(extra)} genes")
    write(f"{INT}/enforcement_changes.tsv", changes)
    for X, Y in (("P_join", "D"), ("P_join_whole", "D"), ("D", "P_ref"), ("D", "P_ref_single"), ("D", "C_L1_refD"),
                 ("P", "C_L1_refP"), ("P", "C_L1_refP_single"), ("D", "Ctree_top_refD"), ("P", "Ctree_top_refP")):
        r = containment(X, Y, "pooled", enf, side={**SIDE, **{g: "x" for g in set(enf[X]) | set(enf[Y]) if g not in SIDE}})
        say(f"[nesting] c({X} ⊇ {Y}) pooled = {r['c_X_contains_Y']} (pairs {r['pairs_Y']}); groups inside one "
            f"{r['groups_Y_inside_one_X']}/{r['groups_Y_ge2']}")
    erow = []
    for g in sorted(set(U) | set(Pjw), key=lambda g: (SIDE.get(g, "~"), genes[g]["name"])):
        erow.append({"gene_id": g, "name": genes[g]["name"], "family_side": SIDE.get(g, "outside U (whole-group JOIN)"),
                     "is_member": "yes" if g in MEMB else "no",
                     **{k: enf[k].get(g, "NA") for k in ("P", "D", "P_join", "P_join_whole", "P_ref", "P_ref_single",
                                                         "C_L1", "C_mid", "C_fine", "Ctree_top", "Ctree_min")}})
    write(f"{INT}/enforced_partitions.tsv", erow)

    # ------------------------------------------------ truths (A): on U
    hg = {g: r["hgnc_gene_group_id"] for g, r in U.items() if r["hgnc_gene_group_id"]}
    so = {g: r["soto_families"] for g, r in U.items() if r["soto_flag"] == "ok" and r["soto_families"]}
    e0 = {g: r["E0_group"] for g, r in U.items() if r["E0_group"]}
    lit = {r["gene_id"]: r for r in tsv(f"{LIGHT}/truth_literature_subfamilies.corrected.tsv")}
    litL1 = {g: r["level1"] for g, r in lit.items() if r["level1"] and r["family"] == "NPIP"}
    litL2 = {g: r["level2"] for g, r in lit.items() if r["level2"]}
    named = {g: ("named" if r["npipb_named_subfamily"] == "yes" else f"other:{g}") for g, r in lit.items()
             if r["family"] == "NPIP" and r["in_literature_truth"] == "yes"}
    all_hg = hgnc_all(genes)
    import soto_map
    db = soto_map.load()
    exons = {r["gene_id"]: soto_map.parse_blocks(r["exons"]) for r in tsv(f"{LIGHT}/work/refseq/exons.tsv")}
    for g in set(Pjw) - set(U):
        if g in all_hg:
            hg[g] = all_hg[g]
        fs = soto_ok(db, exons, genes, g)
        if fs:
            so[g] = fs
    say(f"[truth A] whole-group JOIN genes outside U: HGNC {sorted((genes[g]['name'], hg[g]) for g in set(Pjw) - set(U) if g in hg)}"
        f"; Soto ok {sorted((genes[g]['name'], so[g]) for g in set(Pjw) - set(U) if g in so)}")
    PU = set(LP)
    PU_join_w = {g for g in Pjw if g in pidx_genes} | PU
    trows = []
    plan = [("P", "as built", LP, PU), ("P", "JOIN over D, U-restricted (P universe)", Pj, PU),
            ("P", "JOIN over D, U-restricted (own universe)", Pj, set(Pj)),
            ("P", "JOIN over D, whole groups (P universe, re-closed)", Pjw, PU_join_w),
            ("P", "JOIN over D, whole groups (own universe, re-closed)", Pjw, set(Pjw)),
            ("P", "REFINE in D, USP6NL attached", Pr, PU), ("P", "REFINE in D, USP6NL single", Prs, PU),
            ("D", "as built = after (D unchanged)", LD, set(LD))]
    for layer, variant, lab, uni in plan:
        for tname, truth in (("HGNC gene_group_id (superfamily-level for TBC1D3)", hg), ("Soto family (flag ok)", so)):
            for fam in FAMS:
                side = {g: SIDE.get(g, "TBC1D3") for g in uni}  # whole-group JOIN genes are all on the TBC1D3 side
                r = score(lab, truth, fam_genes(fam, side) & uni)
                trows.append({"layer": layer, "variant": variant, "truth": tname, "family": fam, **r})
    for layer, variant, lab, tname, truth in (
            ("D", "as built", LD, "E0 guided catalog (construction sensitivity, not truth)", e0),
            ("C_L1", "as built", LAYERS["C_L1"], "literature L1 NPIPA|NPIPB (CIRCULAR)", litL1),
            ("C_mid_ab", "as built", C_MID_AS_BUILT, "literature named NPIPB subfamily (CIRCULAR)", named),
            ("C_fine", "as built", LAYERS["C_fine"], "literature L2 paralog groups (CIRCULAR)", litL2),
            ("Ctree_top (clause 5)", "as built", LAYERS["Ctree_top"], "literature L1 (NPIP only)", litL1),
            ("Ctree_top (clause 5)", "as built", LAYERS["Ctree_top"], "literature L2 paralog groups", litL2),
            ("Ctree_min (clause 5)", "as built", LAYERS["Ctree_min"], "literature L2 paralog groups", litL2),
            ("Ctree_root (clause 5, rooted variant)", "as built", CT_LAB["Ctree_root"], "literature L1 (NPIP only)", litL1)):
        for fam in FAMS:
            r = score(lab, truth, fam_genes(fam) & set(lab))
            trows.append({"layer": layer, "variant": variant, "truth": tname, "family": fam, **r})
    trows.append({"layer": "C_mid (= C_mid_ab ∨ C_fine)", "variant": "not scored",
                  "truth": "named NPIPB ∨ L2: identical by construction (both sides are the same join)", "family": "-"})
    cols = ["layer", "variant", "truth", "family", "n_genes", "truth_pairs", "pred_pairs", "tp_pairs", "pair_precision",
            "pair_recall", "bip_F_count(§6ks)", "bip_R", "bip_P", "bip_F_jaccard"]
    write(f"{INT}/truth_agreement.tsv", trows, cols)
    for r in trows:
        if r["family"] in ("TBC1D3", "pooled", "NPIP") and r.get("n_genes"):
            say(f"[truth A] {r['layer']} | {r['variant']} | {r['truth'][:30]} | {r['family']} n={r['n_genes']} "
                f"P {r['pair_precision']} R {r['pair_recall']} bipF {r['bip_F_count(§6ks)']} (pairs pred {r['pred_pairs']} "
                f"truth {r['truth_pairs']})")

    # ------------------------------------------------ truths (B): member-anchored, layer-independent gene sets
    hg_parts = {g: frozenset(v.split("|")) for g, v in all_hg.items()}
    mem_hg = set().union(*(hg_parts.get(m, frozenset()) for m in MEMB))
    hg_group_genes = {g for g, p in hg_parts.items() if p & mem_hg}
    n_2227_sym = sum(1 for r in tsv(HGNC) if "2227" in r["gene_group_id"].split("|"))
    say(f"[truth B] HGNC groups holding a member: {sorted(mem_hg)}; HGNC symbols in those groups {n_2227_sym}; RefSeq "
        f"genes in them {len(hg_group_genes)}, with a §6ko protein {len(hg_group_genes & pidx_genes)}, catalog nodes "
        f"{len(hg_group_genes & set(cat))}")
    so_parts = {g: frozenset(v.split(";")) for g, v in so.items() if g in U}
    mem_so = set().union(*(so_parts.get(m, frozenset()) for m in MEMB))
    cand = {r["best_refseq_gene_id"] for r in tsv(f"{LIGHT}/truth_soto_families.tsv")
            if r["family_id"] in mem_so and r["best_refseq_gene_id"]}
    n_soto_genes = sum(1 for r in tsv(f"{LIGHT}/truth_soto_families.tsv") if r["family_id"] in mem_so)
    for g in cand - set(U):
        fs = soto_ok(db, exons, genes, g)
        if fs and set(fs.split(";")) & mem_so:
            so_parts[g] = frozenset(fs.split(";"))
    so_group_genes = {g for g, p in so_parts.items() if p & mem_so}
    say(f"[truth B] Soto families holding a member (flag ok): {sorted(mem_so)}; Soto genes in them {n_soto_genes}; RefSeq "
        f"genes forward-mapped into them (flag ok) {len(so_group_genes)} (outside U "
        f"{sorted(genes[g]['name'] for g in so_group_genes - set(U))})")
    brows = []
    for tname, parts, grp_genes in (("HGNC gene group", hg_parts, hg_group_genes), ("Soto family", so_parts, so_group_genes)):
        for layer, variant, lab, uni in (
                ("P", "as built", LP, pidx_genes), ("P", "JOIN U-restricted", Pj, pidx_genes),
                ("P", "JOIN whole groups", Pjw, pidx_genes), ("P", "REFINE attach", Pr, pidx_genes),
                ("P", "REFINE single", Prs, pidx_genes), ("D", "as built", LD, set(cat))):
            for fam in ("NPIP", "TBC1D3", "pooled"):
                mem_f = {m for m in MEMB if fam == "pooled" or SIDE[m] == fam}
                gset = ((set(U) | grp_genes) & uni)
                # a gene outside U is outside every member group of every layer here, except whole-group JOIN genes
                labx = dict(lab)
                if layer == "D":
                    for g in gset - set(labx):
                        if g in cat:
                            labx[g] = f"D|{cat[g][0]}|{cat[g][1]}" if cat[g][1] else f"D|single:{g}"
                r = score_member_anchored(labx, parts, gset, mem_f & gset)
                brows.append({"truth": tname, "layer": layer, "variant": variant, "family": fam,
                              "gene_set": "U ∪ all genes of member-holding truth groups, ∩ layer universe genome-wide", **r})
                say(f"[truth B] {tname} | {layer} {variant} | {fam}: genes {r['n_genes']} (members {r['n_members']}) "
                    f"P {r['pair_precision']} R {r['pair_recall']} (pred {r['pred_pairs']} truth {r['truth_pairs']} tp "
                    f"{r['tp_pairs']})")
    write(f"{INT}/truth_member_anchored.tsv", brows)

    # ------------------------------------------------ EXPR
    xrows, t2rows, grows, sweep = [], [], [], []
    layer_set = {"P": LP, "D": LD, "P_join": Pj, "P_ref": Pr, "C_L1": LAYERS["C_L1"], "C_mid": LAYERS["C_mid"],
                 "C_fine": LAYERS["C_fine"], "Ctree_top": LAYERS["Ctree_top"], "Ctree_min": LAYERS["Ctree_min"]}
    # T2 precondition E_M ⊆ E_L inside each L group (layer edges)
    pre = {}
    for L, M in itertools.permutations(layer_set, 2):
        S = set(layer_set[L]) & set(layer_set[M])
        gM = [G for G in groups_of(layer_set[M], S).values() if len(G) >= 2]
        if not gM or any(len({layer_set[L][g] for g in G}) != 1 for G in gM):
            continue
        fM, fL = edge_fn(M), edge_fn(L)
        miss = 0
        tot = 0
        for G in groups_of(layer_set[M], S).values():
            if len(G) < 2:
                continue
            eM, eL = set(fM(G)), set(fL(G))
            tot += len(eM)
            miss += len(eM - eL)
        pre[(L, M)] = (tot, miss)
    for (L, M), (tot, miss) in sorted(pre.items()):
        say(f"[T2 precondition] {L} ⊇ {M} on data: M-edges inside M groups {tot}; not in E_L {miss} -> "
            f"{'holds' if miss == 0 else 'FAILS'}")
    for mode in MODES:
        for t in (1, 2, 3, 4, 5):
            expressed = {g for g in U if READS[g][mode] >= t}
            mem_ex = {f: sum(1 for m in MEMB if SIDE[m] == f and m in expressed) for f in ("NPIP", "TBC1D3")}
            for emode in ("layer_edges", "comembership"):
                EX, splits, distinct = {}, [], set()
                for L, lab in layer_set.items():
                    kind = "comembership" if emode == "comembership" else L
                    comps, dropped = expr_layer(lab, expressed, edge_fn(kind))
                    EX[L] = comps
                    for grp, cs in comps.items():
                        for c in cs:
                            distinct.add(frozenset(c))
                        if len(cs) > 1:
                            splits.append(f"{L}:{grp}: " + " | ".join(",".join(sorted(NAME[g] for g in c)) for c in cs))
                    if t == 3:
                        for grp, G in groups_of(lab, set(lab)).items():
                            ex = sorted(G & expressed)
                            if len(G) < 2 or not (G & MEMB):
                                continue
                            xrows.append({"expr_mode": mode, "edges": emode, "layer": L, "group": grp, "n_genes": len(G),
                                          "n_expressed": len(ex), "n_EXPR_groups": len(comps[grp]),
                                          "EXPR_groups": " | ".join(",".join(sorted(NAME[g] for g in c))
                                                                    for c in sorted(comps[grp], key=lambda c: -len(c))),
                                          "expressed_dropped_singletons": ",".join(NAME[g] for g in dropped[grp]),
                                          "split": "yes" if len(comps[grp]) > 1 else "no"})
                        bad = sum(1 for grp, cs in comps.items() for c in cs if len({lab[g] for g in c}) != 1)
                        grows.append({"expr_mode": mode, "edges": emode, "layer": L,
                                      "EXPR_groups": sum(len(cs) for cs in comps.values()),
                                      "not_inside_one_L_group (0 by construction)": bad})
                nchk = nviol = 0
                for L, M in itertools.permutations(layer_set, 2):
                    if (L, M) not in pre:
                        continue
                    S = set(layer_set[L]) & set(layer_set[M])
                    eL = expr_labels(EX[L])
                    viol, checked = [], 0
                    for c in (c for cs in EX[M].values() for c in cs):
                        cS = c & S
                        if len(cS) < 2:
                            continue
                        checked += 1
                        if len({eL.get(g, f"none:{g}") for g in cS}) != 1:
                            viol.append(",".join(sorted(NAME[g] for g in cS)))
                    nchk += checked
                    nviol += len(viol)
                    if t == 3:
                        t2rows.append({"expr_mode": mode, "edges": emode, "coarse_L": L, "fine_M": M,
                                       "precondition_E_M_in_E_L": ("guaranteed (co-membership)" if emode == "comembership"
                                                                   else "holds" if pre[(L, M)][1] == 0 else
                                                                   f"fails ({pre[(L, M)][1]} of {pre[(L, M)][0]} M edges)"),
                                       "EXPR_M_groups_checked": checked, "violations": len(viol),
                                       "detail": " || ".join(viol)})
                sweep.append({"expr_mode": mode, "t": t, "edges": emode, "NPIP_members_expressed": mem_ex["NPIP"],
                              "TBC1D3_members_expressed": mem_ex["TBC1D3"], "U_genes_expressed": len(expressed),
                              "distinct_EXPR_groups": len(distinct), "L_groups_split": len(splits),
                              "T2_checks": nchk, "T2_violations": nviol,
                              "splits": " || ".join(splits)})
                say(f"[EXPR sweep] {mode:9s} t>={t} {emode:12s}: members NPIP {mem_ex['NPIP']}/27 TBC1D3 "
                    f"{mem_ex['TBC1D3']}/19; U expressed {len(expressed)}; distinct EXPR groups {len(distinct)}; L groups "
                    f"split {len(splits)}; T2 checks {nchk} violations {nviol}" + (f"; splits: {splits}" if splits else ""))
    write(f"{INT}/expr_groups.tsv", xrows)
    write(f"{INT}/expr_nesting.tsv", grows)
    write(f"{INT}/expr_T2.tsv", t2rows)
    write(f"{INT}/expr_sweep.tsv", sweep)
    um = []
    for g in sorted(MEMB, key=lambda g: (SIDE[g], NAME[g])):
        rd = READS[g]
        um.append({"gene_id": g, "name": NAME[g], "family": U[g]["member_family"], "biotype": U[g]["biotype"],
                   "n_reads_any": rd["any"], "n_reads_unique": rd["unique"], "n_reads_unique_mr": rd["unique_mr"],
                   "expressed_any_ge3": "yes" if rd["any"] >= 3 else "no",
                   "expressed_unique_ge3": "yes" if rd["unique"] >= 3 else "no",
                   "expressed_unique_mr_ge3": "yes" if rd["unique_mr"] >= 3 else "no",
                   "layers": ",".join(L for L in ("P", "D", "C_L1", "Ctree_top") if g in LAYERS[L])})
    write(f"{INT}/member_expression.tsv", um)
    for fam in ("NPIP", "TBC1D3"):
        sub = [r for r in um if r["family"] == fam]
        say(f"[EXPR] {fam} members {len(sub)}: >=3 any {sum(r['expressed_any_ge3'] == 'yes' for r in sub)}, unique "
            f"{sum(r['expressed_unique_ge3'] == 'yes' for r in sub)}, unique_mr "
            f"{sum(r['expressed_unique_mr_ge3'] == 'yes' for r in sub)}; any<3: "
            f"{[r['name'] for r in sub if r['expressed_any_ge3'] == 'no']}")

    # ------------------------------------------------ member reconciliation vs §6jg truth (22 NPIP + 9 TBC1D3 records)
    mrec = []
    mrows = {r["gene_id"]: r for r in tsv(f"{LIGHT}/members.corrected.tsv")}
    for g in sorted(MEMB, key=lambda g: (SIDE[g], NAME[g])):
        r = mrows[g]
        inlit = r["in_lit_truth_31"] == "yes"
        if inlit:
            why = "in §6jg truth"
        elif "readthrough" in r["member_basis"]:
            why = "readthrough record whose family part has no gene record of its own (member rule)"
        elif r["biotype"] == "protein_coding":
            why = f"coding RefSeq LOC record ('{r['description']}'); §6jg's 22 NPIP records predate it"
        else:
            why = f"{r['biotype']} ('{r['description']}'); §6jg's TBC1D3 truth is the 9 protein-coding copies" \
                if r["family"] == "TBC1D3" else f"{r['biotype']} ('{r['description']}'); not among §6jg's 22 records"
        if r["span_inside_member_record"]:
            why += f"; span lies inside member record {r['span_inside_member_record']} (same strand)"
        mrec.append({"member": r["name"], "family": r["family"], "biotype": r["biotype"], "chrom": r["chrom"],
                     "in_6jg_truth": "yes" if inlit else "no", "reason": why})
    write(f"{INT}/member_reconciliation.tsv", mrec)
    for fam in ("NPIP", "TBC1D3"):
        say(f"[members] {fam}: {sum(1 for x in mrec if x['family'] == fam)} members; in §6jg truth "
            f"{sum(1 for x in mrec if x['family'] == fam and x['in_6jg_truth'] == 'yes')}; extra: "
            f"{[x['member'] for x in mrec if x['family'] == fam and x['in_6jg_truth'] == 'no']}")

    # ------------------------------------------------ disagreement lists
    dP, dD, dPD, dC = [], [], [], []
    for g in sorted(U, key=lambda g: (SIDE[g], NAME[g])):
        if g in LP:
            mates = [m for m in MEMB if m != g and LP.get(m) == LP[g]]
            if mates:
                not_d = [m for m in mates if not (g in LD and m in LD and LD[g] == LD[m])]
                if not_d and len(not_d) == len(mates):
                    best = max(mates, key=lambda m: P_EDGES.get(frozenset((g, m)), (0, 0, 0))[0])
                    w, i, c = P_EDGES.get(frozenset((g, best)), (float("nan"),) * 3)
                    ws = [P_EDGES[frozenset((g, m))] for m in mates if frozenset((g, m)) in P_EDGES]
                    dP.append({"gene": NAME[g], "biotype": U[g]["biotype"], "chrom": U[g]["chrom"],
                               "is_member": U[g]["is_member"], "P_group": LP[g],
                               "D_status": LD.get(g, "not in D universe (outside both catalogs)"),
                               "D_folded_into": U[g]["D_folded_into"],
                               "member_D_group": ";".join(sorted({LD[m] for m in mates if m in LD})),
                               "n_member_mates": len(mates), "n_member_P_edges": len(ws),
                               "best_member": NAME[best], "blastp_identity": fmt(i), "coverage_longer": fmt(c),
                               "weight": fmt(w),
                               "identity_range": f"{min(x[1] for x in ws):.3f}-{max(x[1] for x in ws):.3f}" if ws else "",
                               "coverage_range": f"{min(x[2] for x in ws):.3f}-{max(x[2] for x in ws):.3f}" if ws else ""})
        if g in LD:
            mates = [m for m in MEMB if m != g and LD.get(m) == LD[g]]
            if mates and not all(g in LP and m in LP and LP[g] == LP[m] for m in mates):
                pm = [m for m in mates if g in LP and m in LP and LP[g] == LP[m]]
                if pm:
                    continue  # P also groups g with some member: not a D-only co-membership
                es = [(D_EDGES[frozenset((g, m))], m) for m in mates if frozenset((g, m)) in D_EDGES]
                best = max(es)[1] if es else None
                es_any = [(D_EDGES[frozenset((g, m))], m) for m in MEMB if m != g and frozenset((g, m)) in D_EDGES]
                best_any = max(es_any) if es_any else None
                grp_genes = [h for h in LD if h != g and LD[h] == LD[g]]
                ea = [(D_EDGES[frozenset((g, h))], h) for h in grp_genes if frozenset((g, h)) in D_EDGES]
                bany = max(ea)[1] if ea else None
                rt = lambda h: bool(U[h]["D_folded_into"]) or "readthrough" in U[h]["description"]  # noqa: E731
                flags = []
                if U[g]["D_folded_into"]:
                    flags.append("gene folded into another locus")
                if "readthrough" in U[g]["description"]:
                    flags.append("gene is a readthrough record")
                if all(rt(m) for m in mates):
                    flags.append("its only member mates are readthrough/folded records")
                if not es and bany is not None and rt(bany):
                    flags.append(f"no direct member edge; strongest group edge is to readthrough/folded {NAME[bany]}")
                dD.append({"gene": NAME[g], "gene_id": g, "biotype": U[g]["biotype"], "chrom": U[g]["chrom"],
                           "is_member": U[g]["is_member"], "D_group": LD[g], "folded_into": U[g]["D_folded_into"],
                           "P_status": LP.get(g, "not in P universe (non-coding / r2-excluded)"),
                           "n_member_mates": len(mates), "n_member_D_edges_same_group": len(es),
                           "best_member_same_D_group": NAME[best] if best else "", "best_member_id": best or "",
                           "D_weight_same_group": fmt(max(es)[0]) if es else "",
                           "best_member_any_D_group": f"{NAME[best_any[1]]} ({best_any[0]:.3f}, {LD.get(best_any[1])})"
                           if best_any else "",
                           "best_group_partner_if_no_member_edge": "" if es or bany is None else
                           f"{NAME[bany]} ({max(ea)[0]:.3f})",
                           "readthrough_fold_flags": "; ".join(flags)})
    keyrow = {}
    for r in D_EDGE_ROWS:
        keyrow[frozenset((r["u_gene_id"], r["v_gene_id"]))] = (r["catalog"], r["u_key"], r["v_key"], r["u_gene_id"])
    want = []
    for r in dD:
        k = frozenset((r["gene_id"], r["best_member_id"]))
        if r["best_member_id"] and k in keyrow:
            c_, uk, vk, _ = keyrow[k]
            want.append((c_, uk, vk))
    st = dna_pair_stats(want)
    n_ok = n_chk = 0
    for r in dD:
        k = frozenset((r["gene_id"], r["best_member_id"]))
        r["DNA_identity"], r["DNA_cov_longer_exonic"], r["id_x_cov_equals_weight"] = "", "", ""
        if r["best_member_id"] and k in keyrow:
            c_, uk, vk, _ = keyrow[k]
            s = st.get((c_, frozenset((uk, vk))))
            if s:
                r["DNA_identity"], r["DNA_cov_longer_exonic"] = fmt(s[0]), fmt(s[1])
                eq = abs(s[0] * s[1] - float(r["D_weight_same_group"])) < 5e-4
                r["id_x_cov_equals_weight"] = "yes" if eq else f"no ({s[0] * s[1]:.3f})"
                n_chk += 1
                n_ok += eq
    say(f"[disagree] D-only genes {len(dD)}; DNA identity x coverage (best member edge in the same D group) reproduces the "
        f"D weight for {n_ok}/{n_chk}")
    for g in sorted(MEMB, key=lambda g: (SIDE[g], NAME[g])):
        if g in LP and g in LD:
            mp = {NAME[m] for m in U if m != g and m in LP and m in LD and LP[m] == LP[g]}
            md = {NAME[m] for m in U if m != g and m in LP and m in LD and LD[m] == LD[g]}
            if mp != md:
                dPD.append({"member": NAME[g], "P_group": LP[g], "D_group": LD[g],
                            "P_only_comates": ",".join(sorted(mp - md)), "D_only_comates": ",".join(sorted(md - mp)),
                            "member_comates_P_only": ",".join(sorted(x for x in mp - md if "gene-" + x in MEMB)),
                            "member_comates_D_only": ",".join(sorted(x for x in md - mp if "gene-" + x in MEMB))})
    for fam in ("NPIP", "TBC1D3"):
        cnt = collections.Counter(LD[m] for m in MEMB if SIDE[m] == fam and m in LD)
        main_g = cnt.most_common(1)[0][0]
        for m in sorted(MEMB, key=lambda x: NAME[x]):
            if SIDE[m] == fam and m in LD and LD[m] != main_g:
                dPD.append({"member": NAME[m], "P_group": LP.get(m, "not in P universe"), "D_group": LD[m],
                            "P_only_comates": "", "D_only_comates": f"outside the family's main D group {main_g} "
                                                                    f"(folded into {U[m]['D_folded_into'] or '-'})",
                            "member_comates_P_only": "", "member_comates_D_only": ""})
    for cl in ("C_L1", "C_mid", "C_fine", "Ctree_top", "Ctree_min"):
        LC = LAYERS[cl]
        for other_name, LO in (("P", LP), ("D", LD)):
            S = set(LC) & set(LO)
            bad = sorted(pairs_of(LC, S) - pairs_of(LO, S))
            dC.append({"clade_level": cl, "vs": other_name, "clade_pairs_on_shared_genes": len(pairs_of(LC, S)),
                       "clade_pairs_split_by_other": len(bad),
                       "examples": "; ".join(f"{NAME[a]}-{NAME[b]}" for a, b in bad[:10]),
                       "clade_genes_outside_other_universe": ",".join(sorted(NAME[g] for g in set(LC) - set(LO)))})
    write(f"{INT}/disagree_P_not_D.tsv", dP)
    write(f"{INT}/disagree_D_not_P.tsv", [{k: v for k, v in r.items() if k not in ("gene_id", "best_member_id")}
                                          for r in dD])
    write(f"{INT}/disagree_members_P_vs_D.tsv", dPD)
    write(f"{INT}/disagree_clades.tsv", dC)
    say(f"[disagree] P-with-member-not-D {len(dP)}: {[r['gene'] for r in dP]}")
    say(f"[disagree] members with different P vs D co-mates / off-main D group {len(dPD)}")
    for r in dC:
        say(f"[disagree] {r['clade_level']} vs {r['vs']}: clade pairs {r['clade_pairs_on_shared_genes']}, split "
            f"{r['clade_pairs_split_by_other']}; clade genes outside {r['vs']}: {r['clade_genes_outside_other_universe']}")
    with open(f"{INT}/analysis.out", "w") as fh:
        fh.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
