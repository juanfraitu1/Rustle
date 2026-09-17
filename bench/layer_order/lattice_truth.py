#!/usr/bin/env python3
"""Nested edge-test lattice, step 4: agreement of every level with the truths, side by side with the prior study's
report-only groupings (P MCL, D = E1 MCL, literature C, clause-5 C_tree) on the SAME genes.

Metrics (same definitions as bench/layer_order/lo_analysis.score, re-implemented with group counts for speed and checked
against integrate_slim/truth_agreement.tsv): genes = those with a truth label (and in the prediction's universe); pairwise
precision / recall over same-group pairs; bipartite F with ONE-TO-ONE JACCARD matching (lo_analysis.bip_jaccard: Hungarian
assignment maximising summed Jaccard; R = matched genes / genes, P = matched genes / genes of matched predicted groups);
F is NA when the truth or the prediction has 0 same-group pairs.

Truths (per gene, lattice/nodes.tsv): HGNC gene_group_id (superfamily-level for TBC1D3); Soto family (flag ok: matched,
not weak, not ambiguous — exon-overlap mapping of light/scripts/soto_map.py); literature L1 NPIPA|NPIPB (NPIP only);
literature L2 paralog groups (NPIP A/B groups; TBC1D3 positional AE/CDKL as in the prior study).

Gene sets:
  U     the 68-gene universe of the prior study, per family side (NPIP / TBC1D3 / pooled).
        'U all'   : every U gene with a truth label (a lattice node without edges is its own group)
        'U ∩ X'   : the report-only layer X's universe (P: coding; D: E1 catalog nodes; C: tree leaves)
  V     the closure V (8,070 genes), pooled; lattice levels only (recall conditioned on V).
  group the member-holding group of each family anchor (NPIPB2, TBC1D3) at each level: pair precision of ALL its genes
        with a truth label (members and pulled-in non-members), with the truth-family composition. U-scope scores only see
        members; this is the in-group precision that must be read next to them.
Circularity: Soto families are SD98 (>= 98% identity) duplications with a shared-exon map-back, the conventions L2
(shared-exon >= 0.30) and L3 (identity >= 0.98) test, so Soto agreement at L2/L3 is partly by construction. Clause 5
(C_tree) was developed on NPIP (§6jp-§6js) and is not an independent comparator on these families.
Outputs: lattice/truth.tsv, lattice/truth_ingroup.tsv, lattice/truth.out
"""
import collections
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench/layer_order")
from lattice_common import INT, OUT, tsv, write  # noqa: E402

LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def c2(n):
    return n * (n - 1) // 2


def score(pred, truth, genes):
    gs = [g for g in genes if g in pred and g in truth]
    out = {"n_genes": len(gs)}
    if not gs:
        return out
    cell = collections.Counter((pred[g], truth[g]) for g in gs)
    pc = collections.Counter(pred[g] for g in gs)
    tc = collections.Counter(truth[g] for g in gs)
    tp = sum(c2(n) for n in cell.values())
    npp = sum(c2(n) for n in pc.values())
    ntp = sum(c2(n) for n in tc.values())
    out.update({"truth_pairs": ntp, "pred_pairs": npp, "tp_pairs": tp,
                "pair_precision": tp / npp if npp else "NA", "pair_recall": tp / ntp if ntp else "NA"})
    if not npp or not ntp:
        out.update({"bip_R_jaccard": "NA", "bip_P_jaccard": "NA", "bip_F_jaccard": "NA (a side has 0 pairs)"})
        return out
    T, P = sorted(tc), sorted(pc)
    ti, pi = {t: i for i, t in enumerate(T)}, {p: j for j, p in enumerate(P)}
    J = np.zeros((len(T), len(P)))
    M = {}
    for (p, t), n in cell.items():
        J[ti[t], pi[p]] = n / (tc[t] + pc[p] - n)
        M[(ti[t], pi[p])] = n
    r, c = linear_sum_assignment(-J)
    matched = sum(M.get((i, j), 0) for i, j in zip(r, c))
    msize = sum(pc[P[j]] for i, j in zip(r, c) if M.get((i, j), 0) > 0)
    R = matched / len(gs)
    Pp = matched / msize if msize else float("nan")
    F = 2 * R * Pp / (R + Pp) if R + Pp and Pp == Pp else float("nan")
    out.update({"bip_R_jaccard": R, "bip_P_jaccard": Pp, "bip_F_jaccard": F})
    return out


nodes = {r["gene_id"]: r for r in tsv(f"{OUT}/nodes.tsv")}
grp = {r["gene_id"]: r for r in tsv(f"{OUT}/groups.tsv")}
V = set(nodes)
U = {g for g, r in nodes.items() if r["in_U"] == "yes"}
SIDE = {g: nodes[g]["family_side"] for g in U}

LATTICE = {  # label -> groups.tsv column
    "L0": "primary|L0", "L1": "primary|L1", "L2": "primary|L2", "L3": "primary|L3",
    "L0Δ": "triangle|L0", "L1Δ": "triangle|L1", "L2Δ": "triangle|L2", "L3Δ": "triangle|L3",
    "L1[E1@.80/.50]": "L1=E1_at_0.80/0.50|L1", "L2[E1@.80/.50]": "L1=E1_at_0.80/0.50|L2", "L3[E1@.80/.50]": "L1=E1_at_0.80/0.50|L3",
    "L1[no v-exon/strand]": "L1=c2_no_vexon_no_strand|L1", "L2[no v-exon/strand]": "L1=c2_no_vexon_no_strand|L2",
    "L3[no v-exon/strand]": "L1=c2_no_vexon_no_strand|L3",
    "L3[w98 gap-incl]": "L3=w98_gap-inclusive|L3", "L3[pooled gap-excl]": "L3=pooled_gap-excluded|L3",
    "L3[pooled gap-incl]": "L3=pooled_gap-inclusive|L3", "L3[id S2]": "L3=S2_SD98_mapback|L3",
    "17:03 L0": "17:03_tests_exact|L0", "17:03 L1": "17:03_tests_exact|L1", "17:03 L2": "17:03_tests_exact|L2",
    "17:03 L3": "17:03_tests_exact|L3", "17:03 L3Δ": "17:03_tests_exact_triangle|L3",
}
LAT = {k: {g: grp[g][c] for g in V} for k, c in LATTICE.items()}
REPORT = {
    "P (§6ko MCL)": ({g: nodes[g]["P_group"] for g in U if nodes[g]["in_P_universe"] == "yes"}),
    "D (E1 MCL)": ({g: nodes[g]["D_group"] for g in U if nodes[g]["in_D_universe"] == "yes"}),
    "C_L1 (lit, circular)": ({g: nodes[g]["C_L1"] for g in U if nodes[g]["in_C_universe"] == "yes"}),
    "C_fine (lit, circular)": ({g: nodes[g]["C_fine"] for g in U if nodes[g]["in_C_universe"] == "yes"}),
    "Ctree_top (clause 5)": ({g: nodes[g]["Ctree_top"] for g in U if nodes[g]["Ctree_top"]}),
    "Ctree_min (clause 5)": ({g: nodes[g]["Ctree_min"] for g in U if nodes[g]["Ctree_min"]}),
}
TRUTH = {
    "HGNC gene group": {g: r["hgnc_gene_group_id"] for g, r in nodes.items() if r["hgnc_gene_group_id"]},
    "Soto family (flag ok)": {g: r["soto_families"] for g, r in nodes.items() if r["soto_flag"] == "ok" and r["soto_families"]},
    "literature L1 (NPIPA|NPIPB)": {g: r["lit_level1"] for g, r in nodes.items() if r["lit_level1"] and r["member_family"] == "NPIP"},
    "literature L2": {g: r["lit_level2"] for g, r in nodes.items() if r["lit_level2"]},
}
for t, d in TRUTH.items():
    say(f"[truth] {t}: genes labelled in V {len(d)}, in U {len(set(d) & U)}")

rows = []


def fam_genes(fam):
    return set(U) if fam == "pooled" else {g for g in U if SIDE[g] == fam}


# ---- check the scorer against the prior study (integrate_slim/truth_agreement.tsv, P as built / D as built, Soto)
prior = {(r["layer"], r["variant"], r["truth"], r["family"]): r for r in tsv(f"{INT}/truth_agreement.tsv")}
for lay, var, name in (("P", "as built", "P (§6ko MCL)"), ("D", "as built = after (D unchanged)", "D (E1 MCL)")):
    for tname, ptname in (("Soto family (flag ok)", "Soto family (flag ok)"),
                          ("HGNC gene group", "HGNC gene_group_id (superfamily-level for TBC1D3)")):
        for fam in ("NPIP", "TBC1D3", "pooled"):
            pr = prior[(lay, var, ptname, fam)]
            lab = REPORT[name]
            s = score(lab, TRUTH[tname], fam_genes(fam) & set(lab))
            fj = s.get("bip_F_jaccard")
            fj = fj if isinstance(fj, str) else ("NA" if fj is None else f"{fj:.3f}")
            mine = (s["n_genes"], s.get("pred_pairs"), s.get("tp_pairs"), fj)
            ok = (str(s["n_genes"]) == pr["n_genes"] and str(s.get("pred_pairs", "")) == pr.get("pred_pairs", "")
                  and str(s.get("tp_pairs", "")) == pr.get("tp_pairs", "")
                  and (fj == pr["bip_F_jaccard"] or (fj.startswith("NA") and pr["bip_F_jaccard"] == "NA")))
            say(f"[check vs prior] {name} {tname} {fam}: mine n {mine[0]} pred {mine[1]} tp {mine[2]} F_jac {mine[3]} | prior n "
                f"{pr['n_genes']} pred {pr.get('pred_pairs')} tp {pr.get('tp_pairs')} F_jac {pr['bip_F_jaccard']} -> {'same' if ok else 'DIFFERENT'}")

# ---- (A) on U, same genes
for tname, truth in TRUTH.items():
    fams = ("NPIP",) if tname.startswith("literature L1") else ("NPIP", "TBC1D3", "pooled")
    for fam in fams:
        F = fam_genes(fam)
        sets = {"U all": F}
        for rname, rlab in REPORT.items():
            sets[f"U ∩ {rname}"] = F & set(rlab)
        for sname, S in sets.items():
            for lname, lab in list(LAT.items()) + ([(sname[4:], REPORT[sname[4:]])] if sname != "U all" else []):
                s = score(lab, truth, S)
                rows.append({"scope": "U", "truth": tname, "family": fam, "gene_set": sname, "layer": lname,
                             "kind": "report-only" if lname in REPORT else "lattice", **s})
# ---- (B) on V (closure), pooled, lattice only
for tname in ("HGNC gene group", "Soto family (flag ok)"):
    for lname, lab in LAT.items():
        s = score(lab, TRUTH[tname], V)
        rows.append({"scope": "V (L0 closure)", "truth": tname, "family": "pooled", "gene_set": "V", "layer": lname,
                     "kind": "lattice", **s})
cols = ["scope", "truth", "family", "gene_set", "layer", "kind", "n_genes", "truth_pairs", "pred_pairs", "tp_pairs",
        "pair_precision", "pair_recall", "bip_R_jaccard", "bip_P_jaccard", "bip_F_jaccard"]
write(f"{OUT}/truth.tsv", rows, cols)
say(f"[write] truth.tsv rows {len(rows)}")

# ---- (C) in-group precision: every labelled gene of the anchor's group (members and pulled-in non-members)
NAME = {g: r["name"] for g, r in nodes.items()}
ANCH = {"NPIP": next(g for g in V if NAME[g] == "NPIPB2"), "TBC1D3": next(g for g in V if NAME[g] == "TBC1D3")}
MEMS = {g for g, r in nodes.items() if r["is_member"] == "yes"}
grows = []
for lname, lab in LAT.items():
    for fam, anc in ANCH.items():
        S = {g for g in V if lab[g] == lab[anc]}
        for tname in ("Soto family (flag ok)", "HGNC gene group"):
            truth = TRUTH[tname]
            lg = sorted(g for g in S if g in truth)
            fams = collections.Counter(truth[g] for g in lg)
            pairs = c2(len(lg))
            same = sum(c2(n) for n in fams.values())
            comp = []
            for f, n in fams.most_common():
                gs = sorted(NAME[g] for g in lg if truth[g] == f)
                comp.append(f"{f}:{n}[{','.join(gs[:12])}{',...' if len(gs) > 12 else ''}]")
            grows.append({"layer": lname, "anchor_family": fam, "group_size": len(S), "group_members": len(S & MEMS),
                          "truth": tname, "labelled_genes": len(lg), "labelled_members": len(set(lg) & MEMS),
                          "labelled_pairs": pairs, "same_family_pairs": same,
                          "in_group_pair_precision": same / pairs if pairs else "NA", "n_truth_families": len(fams),
                          "truth_family_composition": " ".join(comp)})
write(f"{OUT}/truth_ingroup.tsv", grows)
say(f"[write] truth_ingroup.tsv rows {len(grows)}")
with open(f"{OUT}/truth.out", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
