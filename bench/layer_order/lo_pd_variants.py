#!/usr/bin/env python3
"""NPIP/TBC1D3 layer order — how robust is the P ~ D tie? (audit 2026-09-16). Imports lo_analysis (module-level code only
reads files). Every variant re-closes U (members ∪ genes the variant layers place with a member) and scores
c(P ⊇ D), c(D ⊇ P) with the same containment() as the main analysis.

Variants
  P universe rule   spec: every U gene with a §6ko protein (34; PKD1 and DHX40 as 'P|other');
                    verifier: only the 32 genes P clusters with a member; pre-audit tables: P 32 and D 62 (no TBC1D26 row)
  P edge rule       plain §6ko (as built); aa identity >= 0.50, k = 2 MCL member clusters (integrate_slim/P_variant_labels.tsv,
                    stable k = 2..5); aa >= 0.50 exact connected components (threshold-free closure: adds USP6, USP32)
  D membership      as built (a record folded into another locus inherits that locus' cluster); fold-inherited records
                    removed from D's universe (they have no node of their own)
  D annotation      RefSeq E1 (as built); GENCODE/CAT E1 built by the same construction (§6kl: lit/annot_gencode/e1 for
                    chr15/17/22, lit/aj_ho/gencode/e1 for chr16/19/20). RefSeq gene -> GENCODE node with the largest shared
                    exonic bp on the same strand; the node's (folded) cluster is the gene's GENCODE-E1 group.
  leave-one-out / leave-two-out over the 68 U genes (spec reading).
Outputs: integrate_slim/pd_variants.tsv, pd_leave_out.tsv, gencode_map.tsv; stdout -> pd_variants.out
"""
import collections
import csv
import itertools
import sys

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench/layer_order")
import lo_analysis as LA  # noqa: E402

LIGHT, INT, H = LA.LIGHT, LA.INT, LA.H
GENCODE = {"c15_17_22": f"{H}/lit/annot_gencode", "c16_19_20": f"{H}/lit/aj_ho/gencode"}
CHROMS = {"c15_17_22": {"chr15", "chr17", "chr22"}, "c16_19_20": {"chr16", "chr19", "chr20"}}


def tsv(p):
    return LA.tsv(p)


def key_of(s):
    c, r = s.rsplit(":", 1)
    a, b = r.split("-")
    return (c, int(a), int(b))


def ov_bp(a, b):
    i = j = t = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            t += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return t


def blocks(s):
    return sorted(tuple(map(int, x.split("-"))) for x in s.split(","))


def gencode_labels(genes):
    """RefSeq gene_id -> (GENCODE-E1 label, node name, shared exonic bp, shared / min exon-union length)."""
    rex = {r["gene_id"]: blocks(r["exons"]) for r in tsv(f"{LIGHT}/work/refseq/exons.tsv")}
    out = {}
    for tag, d in GENCODE.items():
        names = {r["idx"]: r["name"] for r in tsv(f"{d}/nodes.tsv.names.tsv")}
        rep = {key_of(r["annotation"]): key_of(r["representative"]) for r in tsv(f"{d}/e1.loci.tsv")}
        cl = {(r["chrom"], int(r["start"]), int(r["end"])): r["cluster_id"] for r in tsv(f"{d}/e1.clusters.tsv")}
        by_chrom = collections.defaultdict(list)
        for r in tsv(f"{d}/nodes.tsv"):
            k = (r["chrom"], int(r["start"]) + 1, int(r["end"]))
            by_chrom[r["chrom"]].append((int(r["start"]), int(r["end"]), r["strand"], k, names[r["idx"]],
                                         blocks(r["exons"])))
        for g, r in genes.items():
            if r["chrom"] not in CHROMS[tag] or g not in rex:
                continue
            s0, e0 = int(r["start0"]), int(r["end"])
            ex = rex[g]
            L = sum(b - a for a, b in ex)
            best = None
            for a0, a1, st, k, nm, nex in by_chrom[r["chrom"]]:
                if st != r["strand"] or a1 <= s0 or a0 >= e0:
                    continue
                o = ov_bp(ex, nex)
                if o <= 0:
                    continue
                NL = sum(b - a for a, b in nex)
                cand = (o, o / (L + NL - o), k, nm, o / min(L, NL))
                if best is None or cand[:2] > best[:2]:
                    best = cand
            if best is None:
                continue
            o, jac, k, nm, frac = best
            rk = rep.get(k, k)
            c = cl.get(rk, "")
            out[g] = (f"G|{tag}|{c}" if c else f"G|{tag}|single:{rk[0]}:{rk[1]}-{rk[2]}", nm, o, frac)
    return out


def close_universe(P, D, genes_side_seed):
    """U_v = members ∪ genes sharing a P or D group with a member; side = family of those members."""
    memb = LA.MEMB
    side = {m: LA.SIDE[m] for m in memb}
    for lab in (P, D):
        grp_f = collections.defaultdict(set)
        for m in memb:
            if m in lab:
                grp_f[lab[m]].add(LA.SIDE[m])
        for g, x in lab.items():
            if x in grp_f and "|other" not in x and "single" not in x:
                fs = grp_f[x]
                assert len(fs) == 1, (g, x, fs)
                side.setdefault(g, next(iter(fs)))
    return side


def score_pd(tag, P, D, rows, note=""):
    side = close_universe(P, D, None)
    lay = {"P": {g: x for g, x in P.items() if g in side}, "D": {g: x for g, x in D.items() if g in side}}
    out = []
    for fam in LA.FAMS:
        a = LA.containment("P", "D", fam, lay, side)
        b = LA.containment("D", "P", fam, lay, side)
        v = LA.verdict(a, b)
        rows.append({"variant": tag, "family": fam, "U_size": len(side), "P_universe": len(lay["P"]),
                     "D_universe": len(lay["D"]), "genes_in_both": a["n_genes_in_both"],
                     "c(P>=D)": f"{a['pairs_Y_in_X']}/{a['pairs_Y']}", "c(P>=D)_value": a["c_X_contains_Y"],
                     "D_groups_inside_one_P": f"{a['groups_Y_inside_one_X']}/{a['groups_Y_ge2']}",
                     "c(D>=P)": f"{b['pairs_Y_in_X']}/{b['pairs_Y']}", "c(D>=P)_value": b["c_X_contains_Y"],
                     "P_groups_inside_one_D": f"{b['groups_Y_inside_one_X']}/{b['groups_Y_ge2']}",
                     "verdict": v.replace("P above D", "P > D").replace("D above P", "D > P"), "note": note})
        out.append(rows[-1])
        print(f"[{tag}] {fam:6s} U {len(side)} P {len(lay['P'])} D {len(lay['D'])} genes {a['n_genes_in_both']}: "
              f"c(P⊇D) {a['pairs_Y_in_X']}/{a['pairs_Y']} {a['groups_Y_inside_one_X']}/{a['groups_Y_ge2']} | c(D⊇P) "
              f"{b['pairs_Y_in_X']}/{b['pairs_Y']} {b['groups_Y_inside_one_X']}/{b['groups_Y_ge2']} -> {rows[-1]['verdict']}",
              flush=True)
    return out


def main():
    cat, genes = LA.catalog_all()
    pidx = {r["gene_id"] for r in tsv(f"{LIGHT}/work/P/proteins.index.tsv")}
    P0, D0 = dict(LA.LAYERS["P"]), dict(LA.LAYERS["D"])

    def d_all_label(g):
        t, c, k, f = cat[g]
        return f"D|{t}|{c}" if c else f"D|{t}|singleton:{g}"

    D_full = {g: d_all_label(g) for g in cat}          # every catalog gene (closure is recomputed per variant)
    assert all(D_full[g] == x for g, x in D0.items())
    P_full = dict(P0)
    for g in pidx:                                      # coding genes outside U: not with a member
        P_full.setdefault(g, f"P|other:{g}")
    rows = []
    score_pd("as built (spec P-universe rule: every U gene with a protein)", P_full, D_full, rows)
    P32 = {g: x for g, x in P_full.items() if not x.startswith("P|other")}
    score_pd("verifier P-universe rule (only genes P clusters with a member)", P32, D_full, rows)
    D62 = {g: x for g, x in D_full.items() if LA.NAME.get(g) != "TBC1D26"}
    score_pd("pre-audit group tables on disk (P 32 rows; D without the TBC1D26 row)", P32, D62, rows)

    lab = {r["gene_id"]: r for r in tsv(f"{INT}/P_variant_labels.tsv")}
    P50m = {g: (f"P|{lab[g]['P_aa50_mcl']}" if g in lab and lab[g]["P_aa50_mcl"] else f"P|other:{g}") for g in pidx}
    P50c = {g: (f"P|{lab[g]['P_aa50_comp']}" if g in lab and lab[g]["P_aa50_comp"] else f"P|other:{g}") for g in pidx}
    score_pd("P aa>=0.50, k=2 MCL member clusters", P50m, D_full, rows)
    score_pd("P aa>=0.50, exact components (threshold-free closure)", P50c, D_full, rows)

    Dnf = {g: x for g, x in D_full.items() if not cat[g][3]}
    folded_U = sorted(LA.NAME[g] for g in LA.U if g in cat and cat[g][3])
    print(f"[fold] U genes whose catalog record is folded into another locus (removed from D): {folded_U}")
    score_pd("D without fold-inherited membership", P_full, Dnf, rows, note="removed: " + ",".join(folded_U))
    score_pd("P aa>=0.50 MCL + D without folds", P50m, Dnf, rows)

    gl = gencode_labels(genes)
    DG = {g: v[0] for g, v in gl.items()}
    mrows = []
    for g in sorted(set(LA.U) | {g for g, v in gl.items() if any(DG.get(m) == v[0] for m in LA.MEMB)},
                    key=lambda x: (genes[x]["chrom"], int(genes[x]["start0"]))):
        v = gl.get(g)
        mrows.append({"gene_id": g, "name": genes[g]["name"], "chrom": genes[g]["chrom"], "in_U": "yes" if g in LA.U else "no",
                      "refseq_E1_group": D_full.get(g, "not in catalogs"), "gencode_node": v[1] if v else "",
                      "shared_exonic_bp": v[2] if v else "", "shared_frac_of_shorter": f"{v[3]:.3f}" if v else "",
                      "gencode_E1_group": v[0] if v else "no overlapping GENCODE node"})
    LA.write(f"{INT}/gencode_map.tsv", mrows)
    ug = [r for r in mrows if r["in_U"] == "yes"]
    print(f"[gencode] U genes in a RefSeq catalog {sum(1 for r in ug if r['refseq_E1_group'] != 'not in catalogs')}; "
          f"mapped to a GENCODE node {sum(1 for r in ug if r['gencode_node'])}; genes outside U in a member GENCODE group "
          f"{[r['name'] for r in mrows if r['in_U'] == 'no']}")
    for m in ("NPIPA1", "PKD1", "PKD1P6-NPIPP1", "TBC1D3", "TBC1D26", "DHX40", "TBC1D3P1-DHX40P1", "RNFT1-DT"):
        r = next(x for x in mrows if x["name"] == m)
        print(f"   {m}: RefSeq {r['refseq_E1_group']} | GENCODE node {r['gencode_node']} ({r['shared_frac_of_shorter']}) "
              f"-> {r['gencode_E1_group']}")
    score_pd("D = GENCODE E1 (same construction)", P_full, DG, rows)
    score_pd("P aa>=0.50 MCL + D = GENCODE E1", P50m, DG, rows)
    DG5 = {g: v[0] for g, v in gl.items() if v[3] >= 0.5}
    print(f"[gencode] mapping sensitivity: RefSeq genes kept at shared >= 0.5 of the shorter exon union "
          f"{len(DG5)} of {len(DG)}; U genes dropped: {sorted(LA.NAME[g] for g in LA.U if g in DG and g not in DG5)}")
    score_pd("D = GENCODE E1, mapping shared >= 0.5 of the shorter exon union", P_full, DG5, rows)
    LA.write(f"{INT}/pd_variants.tsv", rows)

    # ---------------------------------------------------------------- leave-one-out / leave-two-out (spec reading)
    lay = {"P": LA.LAYERS["P"], "D": LA.LAYERS["D"]}
    base = {f: LA.verdict(LA.containment("P", "D", f, lay), LA.containment("D", "P", f, lay)) for f in LA.FAMS}
    lrows = []
    for g in sorted(LA.U, key=lambda x: LA.NAME[x]):
        l2 = {k: {h: x for h, x in v.items() if h != g} for k, v in lay.items()}
        side = {h: s for h, s in LA.SIDE.items() if h != g}
        v = {f: LA.verdict(LA.containment("P", "D", f, l2, side), LA.containment("D", "P", f, l2, side)) for f in LA.FAMS}
        if v != base:
            lrows.append({"dropped": LA.NAME[g], **{f: v[f] for f in LA.FAMS}})
            print(f"[leave-one-out] drop {LA.NAME[g]}: {v}")
    print(f"[leave-one-out] genes whose removal changes a P-D verdict: {len(lrows)} of {len(LA.U)}")
    dist = collections.Counter()
    for g, h in itertools.combinations(sorted(LA.U), 2):
        l2 = {k: {x: y for x, y in v.items() if x not in (g, h)} for k, v in lay.items()}
        side = {x: s for x, s in LA.SIDE.items() if x not in (g, h)}
        dist[LA.verdict(LA.containment("P", "D", "pooled", l2, side), LA.containment("D", "P", "pooled", l2, side))] += 1
    print(f"[leave-two-out] pooled P-D verdicts over {sum(dist.values())} gene pairs: {dict(dist)}")
    lrows.append({"dropped": "leave-two-out pooled distribution", "pooled": str(dict(dist))})
    LA.write(f"{INT}/pd_leave_out.tsv", lrows, ["dropped", "NPIP", "TBC1D3", "pooled"])

    # ---------------------------------------------------------------- rooted C_tree variant vs P and D
    layc = {"P": LA.LAYERS["P"], "D": LA.LAYERS["D"], "Ctree_root": LA.LAYERS_VARIANT_C["Ctree_root"],
            "Ctree_top": LA.LAYERS["Ctree_top"]}
    for fam in LA.FAMS:
        for r in LA.tournament(list(layc), fam, layc):
            if "Ctree_root" in (r["X"], r["Y"]):
                print(f"[C_tree rooted] {fam} {r['X']} vs {r['Y']}: {r['c(X>=Y)']} [{r['n_pairs_Y']}] {r['groups_Y_inside_X']}"
                      f" | {r['c(Y>=X)']} [{r['n_pairs_X']}] {r['groups_X_inside_Y']} -> {r['verdict']}")


if __name__ == "__main__":
    main()
