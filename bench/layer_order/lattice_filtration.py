#!/usr/bin/env python3
"""Nested edge-test lattice, step 5: identity filtration inside L1 (a threshold filtration = single linkage on FIXED
evidence; it is not a 'more information' view: T3(b) is not tested here).

For each member-holding L1 group (NPIP, TBC1D3) and each threshold t in (L1 itself, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99,
1.00): components of the L1 edges inside the group whose identity field w >= t, with the shared-exon clause OFF
(L1 AND w >= t) and ON (L2 AND w >= t). An edge whose field is missing passes no threshold (unsatisfiable, not imputed).
Fields (all compared unrounded):
  w98_gapexcl     §0★★★.1 single-record w_98, gap-excluded (the primary L3 field: ON at 0.98 = L3 exactly). Monotone.
  w98_gapincl     single-record w_98, gap-inclusive.
  pooled_gapexcl  identity pooled over the E1 records, gap-excluded (the 17:12 report's field). NOT monotone.
  pooled_gapincl  pooled, gap-inclusive. NOT monotone.
w_98 needs a record witnessing shared-exon >= 0.30 of the smaller exon union, so its OFF sweep is not shared-exon-free;
the nested listing therefore shows OFF with pooled_gapexcl and ON with w98_gapexcl (and ON with pooled_gapexcl).

Single-linkage view: for genes x, y of the group, b(x, y) = max over paths of the minimum edge identity (bottleneck; from a
maximum spanning forest). Components at threshold t are exactly the classes of b >= t, so the thresholds form a dendrogram
(nested by construction). A reference group G (literature subfamily, or a clause-5 split) inside its reference set R
APPEARS at t iff it is whole (min over pairs in G of b >= t) and separated from R \\ G (max over x in G, y in R \\ G of
b < t); the appearance interval is (sep, whole], empty when sep >= whole (G fragments before it separates).
Reference sets: NPIP = literature records (Dishuck 2025 truth) that are E1 catalog nodes (NPIPB1P excluded: outside both
catalogs); TBC1D3 = the 9 clause-5 C_tree leaves.

Outputs: lattice/filtration.txt (nested listing + appearance table), lattice/filtration_groups.tsv,
lattice/filtration_appearance.tsv.
usage: lattice_filtration.py                 primary L1 (clause-2 approx with v-exon overlap + strand check)
       lattice_filtration.py --l1 c2_loose   the 17:12 report's L1 (no v-exon/strand requirements), unrounded fields;
                                             outputs get the suffix .17_03_tests_exact
"""
import collections
import csv
import sys

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench/layer_order")
from lattice_common import ID_COL, OUT, SEF_MIN, UF, components, fnum, groups, tests, tsv, write  # noqa: E402

FIELDS = ("w98_gapexcl", "w98_gapincl", "pooled_gapexcl", "pooled_gapincl")

L1_OPT = sys.argv[sys.argv.index("--l1") + 1] if "--l1" in sys.argv else "c2"
SUF = "" if L1_OPT == "c2" else ".17_03_tests_exact"
GRID = [None, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99, 1.00]
nodes = {r["gene_id"]: r for r in tsv(f"{OUT}/nodes.tsv")}
grp = {r["gene_id"]: r for r in tsv(f"{OUT}/groups.tsv")}
GCOL = {"c2": "primary|L1", "c2_loose": "L1=c2_no_vexon_no_strand|L1"}[L1_OPT]
NAME = {g: r["name"] for g, r in nodes.items()}
BYNAME = collections.defaultdict(list)
for g, n in NAME.items():
    BYNAME[n].append(g)
MEM = {g for g, r in nodes.items() if r["is_member"] == "yes"}
FAM = {g: nodes[g]["member_family"] for g in MEM}
OUTL = []


def out(s=""):
    print(s)
    OUTL.append(s)


# L1 edges with identity fields / shared-exon attribute: (a, b, {field: value}, sef)
E1 = []
with open(f"{OUT}/edges.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        t = tests(r, l1=L1_OPT)
        if t[1]:
            E1.append((r["gene_a"], r["gene_b"], {f: fnum(r[ID_COL[f]]) for f in FIELDS}, fnum(r["d_shared_exon_frac"])))
            assert (t[3] == (t[2] and E1[-1][2]["w98_gapexcl"] is not None and E1[-1][2]["w98_gapexcl"] >= 0.98))
L1lab = {g: grp[g][GCOL] for g in nodes}
anchor = {"NPIP": BYNAME["NPIPB2"][0], "TBC1D3": BYNAME["TBC1D3"][0]}
REF = {
    "NPIP": {g for g in nodes if nodes[g]["lit_in_truth"] == "yes" and FAM.get(g) == "NPIP" and nodes[g]["catalog"] != "none"},
    "TBC1D3": {g for g in nodes if nodes[g]["Ctree_top"] and FAM.get(g) == "TBC1D3"},
}


def names(S):
    return sorted(NAME[g] for g in S)


def short(n):
    return n.replace("NPIP", "").replace("TBC1D3", "T3") if n.startswith(("NPIP", "TBC1D3")) else n


def edges_at(S, t, sef_on, which):
    for a, b, w, sef in E1:
        if a in S and b in S:
            if sef_on and (sef is None or sef < SEF_MIN):
                continue
            if t is not None:
                v = w[which]
                if v is None or v < t:
                    continue
            yield a, b


def bottleneck(S, sef_on, which):
    """b(x, y) for all x, y in S via Kruskal on the field (descending); returns dict of dicts restricted to ref ∪ members."""
    es = sorted(((e[2][which], e[0], e[1]) for e in E1 if e[0] in S and e[1] in S and e[2][which] is not None
                 and (not sef_on or (e[3] is not None and e[3] >= SEF_MIN))), reverse=True)
    uf = UF(S)
    members = {g: {g} for g in S}
    b = collections.defaultdict(dict)
    keep = set(REF["NPIP"]) | set(REF["TBC1D3"]) | MEM
    for w, x, y in es:
        rx, ry = uf.find(x), uf.find(y)
        if rx == ry:
            continue
        A = [g for g in members[rx] if g in keep]
        B = [g for g in members[ry] if g in keep]
        for p in A:
            for q in B:
                b[p][q] = w
                b[q][p] = w
        uf.union(rx, ry)
        r = uf.find(rx)
        members[r] = members[rx] | members[ry]
    return b


def appear(G, R, b):
    G = set(G) & R
    rest = R - G
    whole = min((b[x].get(y, float("-inf")) for x in G for y in G if x < y), default=float("inf"))
    sep = max((b[x].get(y, float("-inf")) for x in G for y in rest), default=float("-inf"))
    return sep, whole


def fmt(v):
    if v == float("inf"):
        return "+inf"
    if v == float("-inf"):
        return "-inf"
    return f"{v:.6f}"


rows = []
for fam in ("NPIP", "TBC1D3"):
    S = {g for g in nodes if L1lab[g] == L1lab[anchor[fam]]}
    fam_members = S & MEM
    out(f"=== {fam}: member-holding L1 group ({GCOL}), {len(S)} genes ({len(fam_members)} members); reference set for appearance: "
        f"{len(REF[fam])} genes {names(REF[fam])}")
    for sef_on, field in ((False, "pooled_gapexcl"), (True, "w98_gapexcl"), (True, "pooled_gapexcl")):
        out(f"--- {fam}, shared-exon clause {'ON (L2 AND w >= t)' if sef_on else 'OFF (L1 AND w >= t)'}; field = {field}")
        for t in GRID:
            lab = components(S, list(edges_at(S, t, sef_on, field)))
            G = groups(lab)
            memg = sorted((C for C in G.values() if C & MEM), key=lambda C: (-len(C & MEM), -len(C), names(C)))
            nonmem_only = [C for C in G.values() if not (C & MEM) and len(C) >= 2]
            tag = "L1" if t is None and not sef_on else ("L2" if t is None else f"w>={t:.2f}")
            if sef_on and t == 0.98 and field == "w98_gapexcl":
                tag += " (= L3)"
            venn = " ".join("{" + " ".join(short(n) for n in names(C & MEM)) + (f" +{len(C - MEM)}" if C - MEM else "") + "}"
                            for C in memg)
            out(f"  [{tag}] groups holding members: {len(memg)}; non-member groups (>=2): {len(nonmem_only)}")
            out(f"     {venn}")
            for C in memg:
                if C - MEM and len(C - MEM) <= 25:
                    out(f"       non-members with {{{' '.join(short(n) for n in names(C & MEM))[:60]}...}}: {', '.join(names(C - MEM))}")
                elif C - MEM:
                    out(f"       non-members with {{{' '.join(short(n) for n in names(C & MEM))[:60]}...}}: {len(C - MEM)} genes, "
                        f"e.g. {', '.join(names(C - MEM)[:12])}")
            for C in memg:
                rows.append({"family": fam, "shared_exon": "on" if sef_on else "off", "field": field, "threshold": tag,
                             "group_members": ",".join(names(C & MEM)), "n_members": len(C & MEM), "size": len(C),
                             "nonmembers": ",".join(names(C - MEM)) if len(C - MEM) <= 200 else f"{len(C - MEM)} genes"})
        out()

write(f"{OUT}/filtration_groups{SUF}.tsv", rows)

# ---------------------------------------------------------------------------------------------- appearance thresholds
lit = {g: nodes[g] for g in REF["NPIP"]}
REFGROUPS = [
    ("NPIP", "NPIPA (lit L1)", {g for g in REF["NPIP"] if lit[g]["lit_level1"] == "NPIPA"}),
    ("NPIP", "NPIPB (lit L1)", {g for g in REF["NPIP"] if lit[g]["lit_level1"] == "NPIPB"}),
    ("NPIP", "A6-9 (lit L2)", {g for g in REF["NPIP"] if lit[g]["lit_level2"] == "A6-9"}),
    ("NPIP", "B3-5 (lit L2)", {g for g in REF["NPIP"] if lit[g]["lit_level2"] == "B3-5"}),
    ("NPIP", "B6-9 (lit L2)", {g for g in REF["NPIP"] if lit[g]["lit_level2"] == "B6-9"}),
    ("NPIP", "B12/13 (lit L2)", {g for g in REF["NPIP"] if lit[g]["lit_level2"] == "B12/13"}),
    ("NPIP", "named NPIPB {B3,B4,B5,B11,B12,B13}", {g for g in REF["NPIP"] if lit[g]["lit_named_npipb"] == "yes"}),
    ("TBC1D3", "clause-5 {B,F,G,H}", {BYNAME[n][0] for n in ("TBC1D3B", "TBC1D3F", "TBC1D3G", "TBC1D3H")}),
    ("TBC1D3", "clause-5 {TBC1D3,D,E,K}", {BYNAME[n][0] for n in ("TBC1D3", "TBC1D3D", "TBC1D3E", "TBC1D3K")}),
]
arows = []
out("=== Appearance thresholds (G appears at t iff whole(G) >= t > sep(G); grid first hit and exact interval)")
for which, sef_on in (("pooled_gapexcl", False), ("pooled_gapexcl", True), ("pooled_gapincl", False),
                     ("w98_gapexcl", True), ("w98_gapincl", True)):
    if True:
        B = {}
        for fam in ("NPIP", "TBC1D3"):
            S = {g for g in nodes if L1lab[g] == L1lab[anchor[fam]]}
            B[fam] = bottleneck(S, sef_on, which)
        for fam, gname, G in REFGROUPS:
            sep, whole = appear(G, REF[fam], B[fam])
            first = next((t for t in GRID[1:] if whole >= t > sep), None)
            arows.append({"identity": which, "shared_exon": "on" if sef_on else "off", "family": fam, "group": gname,
                          "genes": ",".join(names(G)), "sep_max_bottleneck_to_rest": fmt(sep),
                          "whole_min_bottleneck_inside": fmt(whole),
                          "appears_interval": f"({fmt(sep)}, {fmt(whole)}]" if whole > sep else "never (fragments before it separates)",
                          "first_grid_threshold": "none" if first is None else f"{first:.2f}"})
        # A vs B as a split (both sides separated from each other)
        A = next(G for f, n, G in REFGROUPS if n.startswith("NPIPA"))
        Bb = next(G for f, n, G in REFGROUPS if n.startswith("NPIPB"))
        ab = max((B["NPIP"][x].get(y, float("-inf")) for x in A for y in Bb), default=float("-inf"))
        arows.append({"identity": which, "shared_exon": "on" if sef_on else "off", "family": "NPIP",
                      "group": "NPIPA | NPIPB split (no A-B pair together)", "genes": "",
                      "sep_max_bottleneck_to_rest": fmt(ab), "whole_min_bottleneck_inside": "",
                      "appears_interval": f"t > {fmt(ab)}",
                      "first_grid_threshold": next((f"{t:.2f}" for t in GRID[1:] if t > ab), "none")})
        tb = [G for f, n, G in REFGROUPS if f == "TBC1D3"]
        x = max((B["TBC1D3"][p].get(q, float("-inf")) for p in tb[0] for q in tb[1]), default=float("-inf"))
        arows.append({"identity": which, "shared_exon": "on" if sef_on else "off", "family": "TBC1D3",
                      "group": "{B,F,G,H} | {TBC1D3,D,E,K} split (no cross pair together)", "genes": "",
                      "sep_max_bottleneck_to_rest": fmt(x), "whole_min_bottleneck_inside": "", "appears_interval": f"t > {fmt(x)}",
                      "first_grid_threshold": next((f"{t:.2f}" for t in GRID[1:] if t > x), "none")})
for r in arows:
    out(f"  [{r['identity']}, shared-exon {r['shared_exon']}] {r['family']} {r['group']}: sep {r['sep_max_bottleneck_to_rest']}, "
        f"whole {r['whole_min_bottleneck_inside']} -> appears {r['appears_interval']}; first grid threshold {r['first_grid_threshold']}")
write(f"{OUT}/filtration_appearance{SUF}.tsv", arows)

# ---------------------------------------------------------------------------------------------- member dendrogram (merge heights)
out()
out("=== Single-linkage merge heights among members (pooled gap-excluded identity with shared-exon OFF; w98 gap-excluded with "
    "shared-exon ON): each line = a merge of two member sets at bottleneck h (paths may pass through non-members of the L1 group)")
for fam in ("NPIP", "TBC1D3"):
    S = {g for g in nodes if L1lab[g] == L1lab[anchor[fam]]}
    M = sorted(S & MEM, key=lambda g: NAME[g])
    for sef_on, field in ((False, "pooled_gapexcl"), (True, "w98_gapexcl")):
        b = bottleneck(S, sef_on, field)
        pairs = sorted(((b[x].get(y, float("-inf")), x, y) for i, x in enumerate(M) for y in M[i + 1:]), reverse=True)
        uf = UF(M)
        sets = {g: {g} for g in M}
        out(f"--- {fam} members ({len(M)}), shared-exon {'ON' if sef_on else 'OFF'}, field {field}")
        for h, x, y in pairs:
            rx, ry = uf.find(x), uf.find(y)
            if rx == ry:
                continue
            A, Bs = sets[rx], sets[ry]
            if h == float("-inf"):
                roots = {uf.find(g) for g in M}
                out("   never joined by any identity-bearing path: " + " | ".join(
                    "{" + " ".join(short(NAME[g]) for g in sorted(sets[r0], key=lambda g: NAME[g])) + "}" for r0 in sorted(roots)))
                break
            out(f"   h={h:.6f}: {{{' '.join(short(NAME[g]) for g in sorted(A, key=lambda g: NAME[g]))}}} + "
                f"{{{' '.join(short(NAME[g]) for g in sorted(Bs, key=lambda g: NAME[g]))}}}")
            uf.union(rx, ry)
            r = uf.find(rx)
            sets[r] = A | Bs
with open(f"{OUT}/filtration{SUF}.txt", "w") as fh:
    fh.write("\n".join(OUTL) + "\n")
