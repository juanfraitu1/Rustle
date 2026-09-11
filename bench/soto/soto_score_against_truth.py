#!/usr/bin/env python3
"""Score a predicted gene->family_id assignment against Soto's own published truth (soto_famCN_S1C.tsv),
using the FIXED-UNIVERSE methodology this project settled on after retracting an earlier shrinking-set
version of this exact metric (docs/o1_ledger.md, project_soto_full_replication.md 2026-08-02 section):
every gene in the truth universe is scored, including ones the prediction never placed (given a unique
singleton label rather than being silently dropped from the comparison).

Genes with MORE THAN ONE distinct Family ID across their S1C rows (their own table has 149 such genes,
mostly among the ~541 non-eligible-biotype "extra" genes outside their 1,793-gene family-eligible set --
docs/o1_ledger.md §6ih) are EXCLUDED from scoring: a partition-comparison metric needs one ground-truth
label per gene, and picking one of several real, simultaneously-true Family IDs for such a gene would be
arbitrary. This is a disclosed exclusion (reported in the output), not a silent one.

By default scores over ALL 2,185 genes with unambiguous ground truth (Soto's real family-membership
universe, including non-eligible-biotype members) -- the fair, complete comparison. Pass
--eligible-only-universe to instead reproduce the narrower 1,793-gene comparison this project's earlier
sections used (valid as a DIFFERENT, more limited question -- "how well do we cluster the family-eligible
genes" -- not a substitute for the complete one).

Usage: soto_score_against_truth.py --predicted <gene_id,family_id TSV> --truth soto_famCN_S1C.tsv
                                    [--eligible-only-universe soto_1793_geneset.tsv]
"""
import argparse, csv, sys
from collections import defaultdict
from sklearn.metrics import adjusted_rand_score


def load_truth(path):
    """Returns (gene_to_family: clean single-family genes only, ambiguous: set of excluded gene ids)."""
    families = defaultdict(set)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fid = r["Family ID"]
            if fid:
                families[r["Gene ID"]].add(fid)
    ambiguous = {g for g, f in families.items() if len(f) > 1}
    clean = {g: next(iter(f)) for g, f in families.items() if len(f) == 1}
    return clean, ambiguous


def pairs_of(label_map, genes):
    by_label = defaultdict(list)
    for g in genes:
        by_label[label_map[g]].append(g)
    pairs = set()
    for lbl, members in by_label.items():
        if lbl.startswith("__"):
            continue
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                pairs.add(frozenset((members[i], members[j])))
    return pairs


def score(truth, predicted, universe_genes):
    scored = sorted(universe_genes)
    truth_labels, pred_labels = [], []
    for i, g in enumerate(scored):
        t = truth.get(g, "")
        truth_labels.append(t if (t and not t.startswith("Unassigned")) else f"__s{i}")
        p = predicted.get(g, "")
        pred_labels.append(p if p else f"__o{i}")
    ari = adjusted_rand_score(truth_labels, pred_labels)

    truth_map, pred_map = dict(zip(scored, truth_labels)), dict(zip(scored, pred_labels))
    truth_fam = defaultdict(set)
    for g, t in truth_map.items():
        if not t.startswith("__s"):
            truth_fam[t].add(g)
    pred_fam = defaultdict(set)
    for g, p in pred_map.items():
        if not p.startswith("__o"):
            pred_fam[p].add(g)
    pred_sets = {frozenset(v) for v in pred_fam.values()}
    n_exact = sum(1 for v in truth_fam.values() if frozenset(v) in pred_sets)

    truth_pairs, pred_pairs = pairs_of(truth_map, scored), pairs_of(pred_map, scored)
    tp = len(truth_pairs & pred_pairs)
    prec = tp / len(pred_pairs) if pred_pairs else 0.0
    rec = tp / len(truth_pairs) if truth_pairs else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return dict(n_genes=len(scored), ari=ari, n_exact=n_exact, n_truth_fam=len(truth_fam),
                n_pred_fam=len(pred_fam), precision=prec, recall=rec, f1=f1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predicted", required=True, help="TSV with gene_id, family_id columns (empty family_id = unplaced)")
    ap.add_argument("--truth", required=True, help="soto_famCN_S1C.tsv")
    ap.add_argument("--eligible-only-universe", help="restrict scoring to this geneset's gene_id column "
                     "(reproduces the narrower 1,793-gene comparison; omit for the full, honest universe)")
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

    r = score(clean_truth, predicted, universe)
    print(f"scored universe: {r['n_genes']} genes ({len(ambiguous)} ambiguous multi-family genes "
          f"excluded from all scoring, {'restricted to --eligible-only-universe' if a.eligible_only_universe else 'full truth universe'})")
    print(f"ARI: {r['ari']:.4f}")
    print(f"exact family match: {r['n_exact']}/{r['n_truth_fam']} = {100*r['n_exact']/r['n_truth_fam']:.1f}%")
    print(f"pair precision/recall/F1: {r['precision']:.3f}/{r['recall']:.3f}/{r['f1']:.3f}")
    print(f"predicted families: {r['n_pred_fam']}")


if __name__ == "__main__":
    main()
