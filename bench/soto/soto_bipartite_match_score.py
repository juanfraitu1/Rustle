#!/usr/bin/env python3
"""Score a predicted gene->family_id assignment against Soto's own published truth via OPTIMAL BIPARTITE
MATCHING between predicted and true families (Hungarian algorithm, scipy.optimize.linear_sum_assignment),
maximizing total gene overlap over a 1:1 assignment of predicted families to true families.

SCOPE NOTE (standing project rule): this is a bipartite match used ONLY as an evaluation/scoring
technique to compare two ALREADY-COMPUTED partitions after the fact. It plays no role in how families
are defined or how genes are assigned to them anywhere in this pipeline -- the standing rule ("no
bipartite matching or facility-location step" in the family-DEFINITION method) is about modeling choices,
not about how a finished result gets graded against an external benchmark, which is a different, common,
and unrelated use of the same algorithm.

Reports, per matched (predicted, true) family pair: PRECISION (overlap / |predicted family|) and RECALL
(overlap / |true family|, i.e. "sensitivity" -- what fraction of the true family this predicted family
correctly recovered). Aggregated three ways, since they answer different questions and none alone is the
"right" one (this project's own standing METRIC TRAPS list: never pick one summary without naming which):
  - MICRO (gene-weighted): sum(overlap) / sum(|predicted|) and sum(overlap) / sum(|true|) over ALL matched
    pairs -- dominated by large families.
  - MACRO (family-weighted): mean of each matched pair's own precision/recall -- every family counts
    equally regardless of size.
  - Also reports true families that matched NOTHING (0 overlap with their assigned predicted partner,
    i.e. structurally undetected) and predicted families matched to a true family with 0 real overlap
    (spurious clusters the matching still had to pair with something, since scipy's assignment is total).

Reuses soto_score_against_truth.py's load_truth() (same ambiguous-gene exclusion, same singleton
handling) so this script's universe is directly comparable to that script's ARI/exact-match/pair-P-R-F1
numbers -- not a different, incompatible truth definition.

Usage: soto_bipartite_match_score.py --predicted <gene_id,family_id TSV> --truth soto_famCN_S1C.tsv
                                      [--eligible-only-universe soto_1793_geneset.tsv]
"""
import argparse, csv, os, sys
from collections import defaultdict

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from soto_score_against_truth import load_truth  # noqa: E402


def build_families(label_map, universe):
    fam = defaultdict(set)
    for g in universe:
        lbl = label_map.get(g, "")
        if lbl:
            fam[lbl].add(g)
    return fam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predicted", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--eligible-only-universe")
    ap.add_argument("--show-worst", type=int, default=5,
                     help="print this many worst-matched true families (by recall) for inspection")
    a = ap.parse_args()

    clean_truth, ambiguous = load_truth(a.truth)
    predicted = {}
    with open(a.predicted) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("family_id"):
                predicted[r["gene_id"]] = r["family_id"]

    universe = set(clean_truth)
    if a.eligible_only_universe:
        with open(a.eligible_only_universe) as fh:
            restrict = {r["gene_id"] for r in csv.DictReader(fh, delimiter="\t")}
        universe &= restrict

    # "Unassigned_*" truth labels are Soto's own true singletons -- exclude them as TRUE families (a
    # singleton has no meaningful bipartite partner), but leave the genes in the universe for the
    # predicted side's precision accounting (a predicted family claiming one of these genes still pays a
    # precision cost if that gene isn't real overlap with any other true family it's matched against).
    true_fam = defaultdict(set)
    for g in universe:
        t = clean_truth.get(g, "")
        if t and not t.startswith("Unassigned"):
            true_fam[t].add(g)
    pred_fam = build_families(predicted, universe)

    true_ids = sorted(true_fam)
    pred_ids = sorted(pred_fam)
    n_true, n_pred = len(true_ids), len(pred_ids)
    print(f"scored universe: {len(universe)} genes ({len(ambiguous)} ambiguous multi-family genes excluded)")
    print(f"true families: {n_true}   predicted families: {n_pred}")

    # overlap[i][j] = |pred_fam[pred_ids[i]] & true_fam[true_ids[j]]|
    overlap = np.zeros((n_pred, n_true), dtype=int)
    for i, p in enumerate(pred_ids):
        pf = pred_fam[p]
        for j, t in enumerate(true_ids):
            overlap[i, j] = len(pf & true_fam[t])

    # linear_sum_assignment requires a square-ish cost matrix; pad with zero-overlap dummies so every
    # true family gets a (possibly dummy, zero-overlap) predicted partner and vice versa, then match on
    # -overlap (maximize overlap == minimize its negative).
    n = max(n_pred, n_true)
    padded = np.zeros((n, n), dtype=int)
    padded[:n_pred, :n_true] = overlap
    row_ind, col_ind = linear_sum_assignment(-padded)

    matches = []  # (pred_id_or_None, true_id_or_None, overlap)
    for i, j in zip(row_ind, col_ind):
        p = pred_ids[i] if i < n_pred else None
        t = true_ids[j] if j < n_true else None
        matches.append((p, t, int(padded[i, j])))

    micro_num_p = micro_den_p = micro_num_r = micro_den_r = 0
    macro_p_list, macro_r_list = [], []
    zero_overlap_true = []
    for p, t, ov in matches:
        if t is None:
            continue  # a predicted family matched to a dummy (more predicted than true families)
        tsize = len(true_fam[t])
        psize = len(pred_fam[p]) if p is not None else 0
        micro_num_r += ov
        micro_den_r += tsize
        if p is not None:
            micro_num_p += ov
            micro_den_p += psize
            macro_p_list.append(ov / psize if psize else 0.0)
        else:
            macro_p_list.append(0.0)
        macro_r_list.append(ov / tsize if tsize else 0.0)
        if ov == 0:
            zero_overlap_true.append((t, tsize, p))

    micro_prec = micro_num_p / micro_den_p if micro_den_p else 0.0
    micro_rec = micro_num_r / micro_den_r if micro_den_r else 0.0
    macro_prec = sum(macro_p_list) / len(macro_p_list) if macro_p_list else 0.0
    macro_rec = sum(macro_r_list) / len(macro_r_list) if macro_r_list else 0.0

    print()
    print("=== bipartite-matched (Hungarian algorithm, max total overlap) ===")
    print(f"MICRO (gene-weighted)  precision={micro_prec:.3f}  recall/sensitivity={micro_rec:.3f}")
    print(f"MACRO (family-weighted) precision={macro_prec:.3f}  recall/sensitivity={macro_rec:.3f}")
    print(f"true families with ZERO overlap in their matched predicted partner "
          f"(structurally undetected): {len(zero_overlap_true)}/{n_true} "
          f"({100*len(zero_overlap_true)/n_true:.1f}%)")

    if a.show_worst and zero_overlap_true:
        by_size = sorted(zero_overlap_true, key=lambda x: -x[1])[: a.show_worst]
        print(f"\nlargest undetected true families (up to {a.show_worst}):")
        for t, size, p in by_size:
            print(f"  {t}: {size} members, matched to predicted family {p!r} (0 shared genes)")


if __name__ == "__main__":
    main()
