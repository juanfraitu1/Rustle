#!/usr/bin/env python3
"""Benchmark IsoCon-style iterative copy pruning on the sim5x PSV ladder.

Runs copy_assign with/without --iterative-prune for K=0..8, reports detected copy count and
assigned-read accuracy (ground-truth copy is encoded in each read name). Expected effect:
at K=1 the 5 input copies contain a collision (copies 0 and 4 share the only allele), so
pruning should reduce the count toward the identifiable number without hurting accuracy.
"""
import os
import subprocess
import sys

S = "/home/juanfra/winloci_scratch/sim5x"
BIN = "/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign"
KS = [0, 1, 2, 4, 8]


def run_k(k, prune):
    bam = f"{S}/sim5x_K{k}.bam"
    fa = f"{S}/sim5x_K{k}.ref.fa"
    length = open(f"{fa}.fai").read().split("\t")[1].strip()
    out = f"/tmp/ipeval_K{k}_{'prune' if prune else 'base'}"
    cmd = (
        f"bash -c 'ulimit -v 12000000; {BIN} --bam {bam} --fasta {fa} "
        f"--region SIM5X_K{k}:0-{length} --out {out} --min-copies 2"
    )
    if prune:
        cmd += " --iterative-prune"
    cmd += "'"
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # copy count from families.tsv (one family expected)
    families = [l.rstrip("\n").split("\t") for l in open(f"{out}.families.tsv")][1:]
    n_copies = int(families[0][2]) if families else 0

    # assignment accuracy vs read-name truth
    rows = [l.rstrip("\n").split("\t") for l in open(f"{out}.assignments.tsv")][1:]
    total = assigned = strict_correct = 0
    from collections import Counter, defaultdict

    pred_to_true = defaultdict(Counter)
    for r in rows:
        name = r[0]
        true_c = int(name.split("_c")[1].split("_")[0])
        pred_c = int(r[2])
        status = r[3]
        total += 1
        if status == "assigned":
            assigned += 1
            if pred_c == true_c:
                strict_correct += 1
            pred_to_true[pred_c][true_c] += 1

    strict_acc = strict_correct / assigned if assigned else 0.0
    pct_assigned = 100.0 * assigned / total if total else 0.0

    # Cluster accuracy: each predicted copy is mapped to its plurality true copy.
    # This is fair when pruning merges indistinguishable copies and renames the survivors.
    cluster_correct = 0
    for pred_c, counter in pred_to_true.items():
        majority_true = counter.most_common(1)[0][0]
        cluster_correct += counter[majority_true]
    cluster_acc = cluster_correct / assigned if assigned else 0.0

    return n_copies, assigned, strict_correct, strict_acc, pct_assigned, cluster_correct, cluster_acc


def main():
    print(
        "K  | base copies | prune copies | base assigned% | prune assigned% | "
        "base strict | prune strict | base cluster | prune cluster"
    )
    print("-" * 100)
    rows = []
    for k in KS:
        base = run_k(k, False)
        prune = run_k(k, True)
        rows.append((k, base, prune))
        print(
            f"{k:<3}| {base[0]:<12}| {prune[0]:<13}| "
            f"{base[4]:<15.1f}| {prune[4]:<16.1f}| "
            f"{base[3]:<12.3f}| {prune[3]:<13.3f}| {base[6]:<13.3f}| {prune[6]:.3f}"
        )

    # Save TSV for figures.
    with open("iterative_prune_eval.tsv", "w") as fh:
        fh.write(
            "K\tbase_copies\tprune_copies\tbase_assigned_pct\tprune_assigned_pct\t"
            "base_acc\tprune_acc\tbase_cluster_acc\tprune_cluster_acc\n"
        )
        for k, base, prune in rows:
            fh.write(
                f"{k}\t{base[0]}\t{prune[0]}\t{base[4]:.2f}\t{prune[4]:.2f}\t"
                f"{base[3]:.4f}\t{prune[3]:.4f}\t{base[6]:.4f}\t{prune[6]:.4f}\n"
            )
    print("[wrote iterative_prune_eval.tsv]")


if __name__ == "__main__":
    main()
