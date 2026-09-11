#!/usr/bin/env python3
"""The ACTUAL Dennis-lab family-clustering algorithm, reverse-engineered from their own released code
(github.com/mydennislab/HSD_brain_evolution, section_I&II/B_SD98_families.ipynb -- fetched and read
directly, 2026-09-11), not from the paper's own natural-language METHODS description. Every prior
reimplementation in this project (soto_replicate_clustering.py, soto_cluster_from_shared.py's
step5_step6) modeled their step 5 as "build one shared-exon graph, take connected components, split any
component whose famCN MAD is too high into smaller coherent sub-groups" -- a literal reading of "groupings
where the MAD of CN was less than one were selected". Their real code does something structurally
different, confirmed by reading it directly (verbatim quotes below), not by re-guessing from prose:

    get_mad(elements) = stats.median_abs_deviation(median(wssd rows for elements))   # MEDIAN, not MEAN
    low_dispersion_clusters  = [c for c in raw_clusters if get_mad(c) <  1]
    high_dispersion_clusters = [c for c in raw_clusters if get_mad(c) >= 1]           # DISCARDED outright,
                                                                                       # never split further
    genes = [protein_coding/unprocessed_pseudogene elements of any low_dispersion_cluster]
    for gene in genes:
        gene_cluster = [gene]; i = 0
        while True:
            if gene_cluster[i] is protein_coding/unprocessed_pseudogene:
                for cluster in low_dispersion_clusters:
                    if gene_cluster[i] in cluster:
                        gene_cluster = list(set(gene_cluster) | cluster)   # merge the WHOLE cluster in
            i += 1
            if i == len(gene_cluster):
                families.append(gene_cluster); break
    families = dedup(sorted(families))   # the same underlying group is found once per coding seed it holds

So: (1) MAD filtering happens FIRST, as a pass/fail gate on each RAW shared-exon component, not as a
post-hoc split of an over-large one -- a component that fails is dropped from family formation entirely,
not partitioned into smaller MAD-coherent pieces (the "greedy sorted agglomeration" every prior
reimplementation invented for "how do you split it" was answering a question their code never asks).
(2) Two low-dispersion clusters that share ONLY a non-coding gene (lncRNA / processed pseudogene) are
NEVER merged -- expansion propagates only through protein-coding/unprocessed-pseudogene bridge genes.
This is a real, structural difference from plain connected-components over the whole shared-exon graph
(which would merge on ANY shared gene, coding or not), independent of the MAD statistic question
(mean vs median, see mad_mean/mad_median in soto_cluster_from_shared.py) and worth testing separately.

Usage: soto_cluster_dennislab_algorithm.py --shared shared_exons_1793_final.tsv --geneset soto_1793_geneset.tsv
                                            --famcn soto_famCN_S1C.tsv --mad-statistic median --out <path>
"""
import argparse, csv, sys
from collections import defaultdict

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from soto_cluster_from_shared import mad_mean, mad_median, ELIGIBLE  # noqa: E402


def raw_components(shared):
    """Connected components of the shared-exon graph -- the SAME initial object their pipeline calls
    `clusters` (loaded from their data/SD98_exon_clusters.txt), before any MAD filtering."""
    seen, comps = set(), []
    for g in shared:
        if g in seen:
            continue
        stack, comp = [g], []
        seen.add(g)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in shared[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(frozenset(comp))
    return comps


def dennislab_families(comps, biotype, famcn, mad_threshold, mad_fn):
    """Faithful port of the algorithm quoted in this file's own header. Returns a list of frozensets
    (deduplicated families) plus the set of genes that appear in some LOW-dispersion cluster but never
    became part of any final family (i.e. non-coding-only clusters with no coding bridge -- these keep
    membership information but never get a Family ID, matching "SD98 genes associated with other gene
    features... were also assigned a gene family ID" ONLY when they actually join one via a coding bridge)
    and the set of genes whose ONLY cluster was high-dispersion (discarded outright).
    """
    low, high_genes = [], set()
    for comp in comps:
        vals = [famcn[g] for g in comp if g in famcn]
        is_low = len(vals) < 2 or mad_fn(vals) < mad_threshold
        if is_low:
            low.append(comp)
        else:
            high_genes.update(comp)

    seeds = [g for c in low for g in c if biotype.get(g) in ELIGIBLE]

    families_raw = []
    for seed in seeds:
        cluster_list = [seed]
        seen_local = {seed}
        i = 0
        while i < len(cluster_list):
            if biotype.get(cluster_list[i]) in ELIGIBLE:
                for c in low:
                    if cluster_list[i] in c:
                        for g in c:
                            if g not in seen_local:
                                seen_local.add(g)
                                cluster_list.append(g)
            i += 1
        families_raw.append(frozenset(cluster_list))

    families = sorted(set(families_raw), key=lambda f: sorted(f))
    in_family = set().union(*families) if families else set()
    low_genes = set().union(*low) if low else set()
    orphaned_low = low_genes - in_family  # non-coding-only low-dispersion clusters: never bridged in
    return families, orphaned_low, high_genes - in_family


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shared", required=True)
    ap.add_argument("--geneset", required=True)
    ap.add_argument("--famcn", required=True)
    ap.add_argument("--mad", type=float, default=1.0)
    ap.add_argument("--mad-statistic", choices=["mean", "median"], default="median",
                     help="default 'median' here (unlike soto_cluster_from_shared.py's default 'mean'), "
                          "since this script specifically exists to test their REAL algorithm faithfully")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    mad_fn = mad_median if a.mad_statistic == "median" else mad_mean

    genes, biotype = set(), {}
    with open(a.geneset) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes.add(r["gene_id"])
            biotype[r["gene_id"]] = r.get("biotype", "")

    shared = defaultdict(set)
    with open(a.shared) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            shared[r["gene_a"]].add(r["gene_b"])
            shared[r["gene_b"]].add(r["gene_a"])

    famcn = {}
    with open(a.famcn) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = r.get("gene_id") or r.get("Gene ID")
            v = r.get("famCN") or r.get("famCN_median") or r.get("Median famCN")
            try:
                famcn[gid] = float(v)
            except (TypeError, ValueError):
                pass

    comps = raw_components(shared)
    families, orphaned_low, orphaned_high = dennislab_families(comps, biotype, famcn, a.mad, mad_fn)

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "biotype", "family_id", "n_members", "famCN", "status"])
        assigned = set()
        for i, fam in enumerate(families):
            fid = f"DENNISFAM{i}"
            for g in sorted(fam):
                w.writerow([g, biotype.get(g, ""), fid, len(fam),
                            f"{famcn[g]:.2f}" if g in famcn else "", "family"])
                assigned.add(g)
        for g in sorted(genes - assigned):
            status = "low_dispersion_no_coding_bridge" if g in orphaned_low else \
                     "high_dispersion_discarded" if g in orphaned_high else "singleton"
            w.writerow([g, biotype.get(g, ""), "", 1, f"{famcn[g]:.2f}" if g in famcn else "", status])

    n_low = sum(1 for c in comps
                if len([famcn[g] for g in c if g in famcn]) < 2
                or mad_fn([famcn[g] for g in c if g in famcn]) < a.mad)
    print(f"[step4-input] {len(genes)} SD98 genes, {len(shared)} genes with >=1 shared exon, "
          f"{len(comps)} raw components ({n_low} low-dispersion, {len(comps) - n_low} high-dispersion)",
          file=sys.stderr)
    print(f"[done] {len(families)} families (mad-statistic={a.mad_statistic}), "
          f"{len(orphaned_low)} low-dispersion-but-no-coding-bridge, "
          f"{len(orphaned_high)} high-dispersion-discarded -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
