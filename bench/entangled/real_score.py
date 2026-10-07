#!/usr/bin/env python3
"""Score transcript sets on the truth tables of real_truth.py (docs/PREREG_gorilla_overlap_2026-10-07.md).

    real_score.py --truth TRUTHPREFIX --denominator P1|P2 --arm NAME=GTF [--arm NAME=GTF ...] --out OUTPREFIX

Strata of the evaluated genes (those with >= 1 expressed chain under the chosen denominator): E_both (and its subsets E_both_j / E_both_x by junction sharing), E_one, A, N.
Per arm and stratum: genes, chains (expressed, distinct), recovered (exact chain of >= 1 arm transcript, same chromosome and strand), complete genes, match_n / match_d (the annotation-match share: arm multi-exon
transcripts that overlap an evaluated gene of the stratum on its strand, and those among them whose chain equals ANY valid annotated chain), resolved genes (an arm gene_id carries an exact chain of the gene and
an exact valid annotated chain of no other gene) and merged_ids (arm gene_ids carrying exact valid annotated chains of >= 2 genes). Writes OUTPREFIX.<arm>.json lines into OUTPREFIX.json.
"""
import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arms_score as A

LABELS = ["E_both", "E_both_j", "E_both_x", "E_one", "A", "N"]


def load_truth(prefix, denominator="P1"):
    k = 1 if denominator == "P1" else 2
    genes, all_chains = {}, collections.defaultdict(set)
    gstrat = {}
    with open(prefix + ".genes.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            u = [tuple(map(int, x.split("-"))) for x in r["union"].split(",")] if r["union"] else []
            genes[r["gene"]] = dict(gene=r["gene"], name=r["name"], chrom=r["chrom"], strand=r["strand"], union=u, expressed=set(),
                                    stratum=r[f"stratum_P{k}"], junction_sharing=int(r[f"jsharing_P{k}"]))
    with open(prefix + ".chains.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ch = tuple(tuple(map(int, x.split("-"))) for x in r["chain"].split(","))
            all_chains[(r["chrom"], r["strand"], ch)].add(r["gene"])
            if int(r[f"n_P{k}"]) >= 3:
                genes[r["gene"]]["expressed"].add(ch)
    return dict(genes=genes, all_chains=all_chains)


def score_arm(tx, truth):
    genes, all_chains = truth["genes"], truth["all_chains"]
    ev = {gid: g for gid, g in genes.items() if g["expressed"] and g["union"]}
    labels_of = collections.defaultdict(list)
    for gid, g in ev.items():
        labels_of[gid].append(g["stratum"])
        if g["stratum"] == "E_both":
            labels_of[gid].append("E_both_j" if g["junction_sharing"] else "E_both_x")
    idx = A.Index({gid: g for gid, g in ev.items()})
    arm_chains, multi = set(), {}
    for t, d in tx.items():
        if len(d["exons"]) >= 2:
            ch = A.chain_of(d["exons"])
            multi[t] = ch
            arm_chains.add((d["chrom"], d["strand"], ch))
    id_exact = collections.defaultdict(set)
    for t, ch in multi.items():
        owners = all_chains.get((tx[t]["chrom"], tx[t]["strand"], ch))
        if owners:
            id_exact[tx[t]["gene"]] |= owners
    ids_for_gene = collections.defaultdict(set)
    for i, owners in id_exact.items():
        for o in owners:
            ids_for_gene[o].add(i)
    strata = {lab: collections.Counter() for lab in LABELS}
    for gid, g in ev.items():
        got = [c for c in g["expressed"] if (g["chrom"], g["strand"], c) in arm_chains]
        resolved = any(id_exact[i] == {gid} for i in ids_for_gene.get(gid, ()))
        for lab in labels_of[gid]:
            s = strata[lab]
            s["genes"] += 1
            s["chains"] += len(g["expressed"]); s["recovered"] += len(got)
            s["complete"] += int(len(got) == len(g["expressed"])); s["resolved"] += int(resolved)
    seen = collections.defaultdict(set)
    for t, ch in multi.items():
        d = tx[t]
        hit = idx.genes_hit(d["chrom"], d["strand"], d["exons"])
        labs = {lab for gid in hit for lab in labels_of[gid]}
        for lab in labs:
            strata[lab]["match_d"] += 1
            strata[lab]["match_n"] += int((d["chrom"], d["strand"], ch) in all_chains)
    merged = sum(1 for s in id_exact.values() if len(s) >= 2)
    return dict(transcripts=len(tx), multi_exon=len(multi), merged_ids=merged, strata={k: dict(v) for k, v in strata.items()})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--denominator", choices=["P1", "P2"], default="P1")
    ap.add_argument("--arm", action="append", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    truth = load_truth(a.truth, a.denominator)
    out = dict(denominator=a.denominator, evaluated_genes=sum(1 for g in truth["genes"].values() if g["expressed"]), arms={})
    for spec in a.arm:
        name, path = spec.split("=", 1)
        res = score_arm(A.load_arm(path), truth)
        out["arms"][name] = res
        row = " | ".join(f"{lab}: {res['strata'][lab].get('recovered', 0)}/{res['strata'][lab].get('chains', 0)} chains, {res['strata'][lab].get('resolved', 0)}/{res['strata'][lab].get('genes', 0)} resolved, match {res['strata'][lab].get('match_n', 0)}/{res['strata'][lab].get('match_d', 0)}" for lab in ("E_both", "E_both_x", "N"))
        print(f"{name:6} tx {res['transcripts']:7} | {row}")
    json.dump(out, open(a.out + ".json", "w"), indent=1)


if __name__ == "__main__":
    main()
