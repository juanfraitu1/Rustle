#!/usr/bin/env python3
"""vg_family_prototype_fp_characterize.py — characterize FP families in the VG-native catalog.

Rebuilds the default mmseqs-pairs VG, loads the per-family eval TSV, and computes
structural / graph features for each family.  Reports distributions of FP vs clean
families and tests a few empirical gate ideas.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype_fp_characterize.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
import statistics
from collections import defaultdict, Counter

import pysam

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)

from vg_family_prototype import (
    load_meta, load_skeletons, extract_exon_paths, build_pair_vg,
    SCRATCH, FASTA,
)
from family_level_pr_current import build_ctx

CATALOG_TSV = os.path.join(BENCH, "vg_family_prototype.tsv")
EVAL_TSV = os.path.join(BENCH, "vg_family_prototype_eval.tsv")
OUT_TSV = os.path.join(BENCH, "vg_family_prototype_fp_features.tsv")
OUT_JSON = os.path.join(BENCH, "vg_family_prototype_fp_characterize.json")

REPEAT_THRESH = 30


def recip_overlap(a, b):
    ca, sa, ea = a; cb, sb, eb = b
    if ca != cb:
        return 0.0
    ov = min(ea, eb) - max(sa, sb)
    if ov <= 0:
        return 0.0
    la = ea - sa; lb = eb - sb
    if la <= 0 or lb <= 0:
        return 0.0
    return min(ov / la, ov / lb)


def repeat_frac(fa, chrom, start, end):
    seq = fa.fetch(chrom, start, end)
    if not seq:
        return 0.0
    return sum(1 for c in seq if c.islower()) / len(seq)


def repeat_frac_seq(seq):
    if not seq:
        return 0.0
    return sum(1 for c in seq if c.islower()) / len(seq)


def characterize():
    print("[*] building ctx (annotation / protein families) ...", flush=True)
    ctx, *_ = build_ctx()
    gene_of_dn = ctx["gene_of_dn"]
    g2f = ctx["g2f"]
    mega = ctx["mega"]

    print("[*] loading catalog + eval ...", flush=True)
    fam_members = defaultdict(list)
    with open(CATALOG_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            fam_members[int(row["fam_id"])].append(row["member"])

    eval_rows = {}
    with open(EVAL_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            eval_rows[int(row["family_id"])] = row

    print("[*] loading meta + skeletons + genome ...", flush=True)
    meta = load_meta()
    loci = load_skeletons(meta)
    fa = pysam.FastaFile(FASTA)

    print("[*] extracting exon paths ...", flush=True)
    exons, locus_paths = extract_exon_paths(fa, loci)

    print("[*] building default mmseqs-pairs VG (identity=0.90) ...", flush=True)
    node_seq, node_loci, locus_path_nodes = build_pair_vg(
        exons, locus_paths, identity=0.90, threads=4)

    # precompute repeat fraction per exon occurrence from genomic interval
    print("[*] computing per-exon repeat fractions ...", flush=True)
    exon_repeat_frac = {}
    for ex in exons:
        seq = fa.fetch(ex["chrom"], ex["start"], ex["end"])
        exon_repeat_frac[ex["oid"]] = repeat_frac_seq(seq) if seq else 0.0

    # precompute per-node repeat fraction (mean over member exons)
    print("[*] computing per-node repeat fractions ...", flush=True)
    node_repeat_frac = {}
    # mmseqs-pairs mode: node_seq is concatenation of two exon canonical seqs.
    # Map node -> set of (lid, eidx) occurrences via locus_path_nodes and exons.
    # Build oid -> occurrence list.
    oid_occurrences = defaultdict(list)
    for ex in exons:
        oid_occurrences[ex["oid"]].append(ex)

    for node, loci_set in node_loci.items():
        fracs = []
        for lid in loci_set:
            if lid not in locus_path_nodes:
                continue
            nodes = locus_path_nodes[lid]
            # find positions of this node in the path
            for pos, n in enumerate(nodes):
                if n != node:
                    continue
                # pair node spans exon pos and pos+1
                oids = locus_paths[lid]
                if pos < len(oids) - 1:
                    for epos in (pos, pos + 1):
                        oid = oids[epos]
                        fracs.append(exon_repeat_frac.get(oid, 0.0))
        node_repeat_frac[node] = round(statistics.mean(fracs), 4) if fracs else 0.0

    # RNA-only exact-exon-node architecture (no clustering, no protein truth)
    print("[*] building exact exon-node architecture map ...", flush=True)
    canon_to_exon_node = {}
    oid_to_exon_node = {}
    for ex in exons:
        c = ex["canon"]
        if c not in canon_to_exon_node:
            canon_to_exon_node[c] = ex["oid"]
        oid_to_exon_node[ex["oid"]] = canon_to_exon_node[c]

    locus_exon_nodes = {}
    locus_exon_path = {}
    locus_intron_adj = {}
    for lid, oids in locus_paths.items():
        nodes = [oid_to_exon_node[oid] for oid in oids]
        locus_exon_nodes[lid] = set(nodes)
        locus_exon_path[lid] = tuple(nodes)
        locus_intron_adj[lid] = set(zip(nodes[:-1], nodes[1:]))

    print("[*] computing per-family features ...", flush=True)
    features = []
    for fid, members in fam_members.items():
        eval_row = eval_rows.get(fid, {})
        n_distinct = int(eval_row.get("n_loci_distinct", 0))
        n_genes = int(eval_row.get("n_genes", 0))
        ep_impure = int(eval_row.get("ep_impure", 0))
        n_prot_fam = int(eval_row.get("n_protein_families", 0))
        fp_multi = int(eval_row.get("fp_multifam", 0))
        fp_over = int(eval_row.get("fp_oversize", 0))
        fp_allele = int(eval_row.get("fp_allele", 0))
        any_fp = int(fp_multi or fp_over or fp_allele)

        # coordinate / strand / exon stats
        chroms = [meta[m]["chrom"] for m in members if m in meta]
        starts = [meta[m]["start"] for m in members if m in meta]
        ends = [meta[m]["end"] for m in members if m in meta]
        strands = [meta[m]["strand"] for m in members if m in meta]
        nexons = [meta[m]["n_exon"] for m in members if m in meta]

        n_chrom = len(set(chroms))
        same_chrom_pairs = 0
        recip_overlaps = []
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a = members[i]; b = members[j]
                if a in meta and b in meta:
                    if meta[a]["chrom"] == meta[b]["chrom"]:
                        same_chrom_pairs += 1
                        recip_overlaps.append(recip_overlap(
                            (meta[a]["chrom"], meta[a]["start"], meta[a]["end"]),
                            (meta[b]["chrom"], meta[b]["start"], meta[b]["end"])))
        mean_recip_ov = round(statistics.mean(recip_overlaps), 4) if recip_overlaps else 0.0
        max_recip_ov = round(max(recip_overlaps), 4) if recip_overlaps else 0.0
        strand_major = max(Counter(strands).values()) / len(strands) if strands else 0.0

        # pair-node multiplicity stats within family
        node_mults = []
        node_repeats = []
        hub_frac = 0.0
        repeat_hub_frac = 0.0
        for m in members:
            if m in locus_path_nodes:
                for node in locus_path_nodes[m]:
                    mult = len(node_loci.get(node, []))
                    node_mults.append(mult)
                    node_repeats.append(node_repeat_frac.get(node, 0.0))
                    if mult >= REPEAT_THRESH and node_repeat_frac.get(node, 0.0) >= 0.5:
                        repeat_hub_frac += 1
        if node_mults:
            node_mults.sort()
            max_mult = node_mults[-1]
            mean_mult = round(statistics.mean(node_mults), 2)
            median_mult = node_mults[len(node_mults) // 2]
            hub_frac = round(sum(1 for x in node_mults if x >= REPEAT_THRESH) / len(node_mults), 4)
            repeat_hub_frac = round(repeat_hub_frac / len(node_mults), 4)
        else:
            max_mult = mean_mult = median_mult = 0

        # shared-node coverage: fraction of each member's nodes that are shared (mult>=2)
        shared_node_fracs = []
        for m in members:
            if m in locus_path_nodes:
                nodes = locus_path_nodes[m]
                if nodes:
                    shared = sum(1 for node in nodes if len(node_loci.get(node, [])) >= 2)
                    shared_node_fracs.append(shared / len(nodes))
        mean_shared_node_frac = round(statistics.mean(shared_node_fracs), 4) if shared_node_fracs else 0.0
        min_shared_node_frac = round(min(shared_node_fracs), 4) if shared_node_fracs else 0.0
        max_shared_node_frac = round(max(shared_node_fracs), 4) if shared_node_fracs else 0.0

        # within-family pair-node coverage and Jaccard (does the family co-thread, or just share one hub?)
        fam_pair_lids = [m for m in members if m in locus_path_nodes]
        family_pair_counts = Counter()
        for lid in fam_pair_lids:
            family_pair_counts.update(locus_path_nodes[lid])
        family_shared_node_fracs = []
        pair_jaccs = []
        if fam_pair_lids:
            for lid in fam_pair_lids:
                nodes = locus_path_nodes[lid]
                if nodes:
                    shared = sum(1 for node in nodes if family_pair_counts[node] >= 2)
                    family_shared_node_fracs.append(shared / len(nodes))
            for i, a in enumerate(fam_pair_lids):
                sa = set(locus_path_nodes[a])
                for b in fam_pair_lids[i + 1:]:
                    sb = set(locus_path_nodes[b])
                    inter = len(sa & sb)
                    uni = len(sa | sb)
                    pair_jaccs.append(inter / uni if uni else 0.0)
        mean_family_shared_node_frac = round(statistics.mean(family_shared_node_fracs), 4) if family_shared_node_fracs else 0.0
        mean_pair_jaccard = round(statistics.mean(pair_jaccs), 4) if pair_jaccs else 0.0

        # RNA-only exon architecture signals for domain-share prevention
        fam_lids = [m for m in members if m in locus_exon_nodes]
        exon_jaccs = []
        intron_jaccs = []
        max_unshared_runs = []
        if len(fam_lids) >= 2:
            family_exon_counts = Counter()
            for lid in fam_lids:
                family_exon_counts.update(locus_exon_nodes[lid])
            for i, a in enumerate(fam_lids):
                sa = locus_exon_nodes[a]
                ia = locus_intron_adj[a]
                path_a = locus_exon_path[a]
                cur_run = 0
                max_run = 0
                for n in path_a:
                    if family_exon_counts[n] == 1:
                        cur_run += 1
                        max_run = max(max_run, cur_run)
                    else:
                        cur_run = 0
                max_unshared_runs.append(max_run)
                for b in fam_lids[i + 1:]:
                    sb = locus_exon_nodes[b]
                    ib = locus_intron_adj[b]
                    inter = len(sa & sb)
                    uni = len(sa | sb)
                    exon_jaccs.append(inter / uni if uni else 0.0)
                    inter_i = len(ia & ib)
                    uni_i = len(ia | ib)
                    intron_jaccs.append(inter_i / uni_i if uni_i else 0.0)
        mean_exon_jaccard = round(statistics.mean(exon_jaccs), 4) if exon_jaccs else 0.0
        mean_intron_jaccard = round(statistics.mean(intron_jaccs), 4) if intron_jaccs else 0.0
        mean_max_unshared_run = round(statistics.mean(max_unshared_runs), 4) if max_unshared_runs else 0.0

        # repeat fraction (soft-masked) averaged over member spans
        rfracs = [repeat_frac(fa, meta[m]["chrom"], meta[m]["start"], meta[m]["end"]) for m in members if m in meta]
        mean_rfrac = round(statistics.mean(rfracs), 4) if rfracs else 0.0

        # protein family composition
        prot_fams = set()
        for m in members:
            g = gene_of_dn.get(m)
            if g and g in g2f and g2f[g] not in mega:
                prot_fams.add(g2f[g])

        features.append(dict(
            family_id=fid,
            n_members=len(members),
            n_distinct_loci=n_distinct,
            n_genes=n_genes,
            n_chrom=n_chrom,
            same_chrom_pairs=same_chrom_pairs,
            mean_recip_overlap=mean_recip_ov,
            max_recip_overlap=max_recip_ov,
            strand_majority=round(strand_major, 4),
            mean_n_exons=round(statistics.mean(nexons), 2) if nexons else 0,
            mean_pair_mult=mean_mult,
            median_pair_mult=median_mult,
            max_pair_mult=max_mult,
            pair_hub_frac=hub_frac,
            mean_repeat_frac=mean_rfrac,
            mean_node_repeat_frac=round(statistics.mean(node_repeats), 4) if node_repeats else 0.0,
            max_node_repeat_frac=round(max(node_repeats), 4) if node_repeats else 0.0,
            repeat_hub_frac=repeat_hub_frac,
            mean_shared_node_frac=mean_shared_node_frac,
            min_shared_node_frac=min_shared_node_frac,
            max_shared_node_frac=max_shared_node_frac,
            mean_family_shared_node_frac=mean_family_shared_node_frac,
            mean_pair_jaccard=mean_pair_jaccard,
            mean_exon_jaccard=mean_exon_jaccard,
            mean_intron_jaccard=mean_intron_jaccard,
            mean_max_unshared_run=mean_max_unshared_run,
            n_protein_families=len(prot_fams),
            ep_impure=ep_impure,
            fp_multifam=fp_multi,
            fp_oversize=fp_over,
            fp_allele=fp_allele,
            any_dna_fp=any_fp,
        ))

    # write feature TSV
    cols = list(features[0].keys())
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in features:
            w.writerow(row)
    print(f"[write] {OUT_TSV}", flush=True)

    # summarize distributions
    def summarize(group, name):
        if not group:
            return {}
        return {
            "n": len(group),
            "mean_distinct_loci": round(statistics.mean(r["n_distinct_loci"] for r in group), 2),
            "median_distinct_loci": statistics.median(r["n_distinct_loci"] for r in group),
            "mean_genes": round(statistics.mean(r["n_genes"] for r in group), 2),
            "mean_chrom": round(statistics.mean(r["n_chrom"] for r in group), 2),
            "mean_recip_overlap": round(statistics.mean(r["mean_recip_overlap"] for r in group), 4),
            "max_recip_overlap": round(statistics.mean(r["max_recip_overlap"] for r in group), 4),
            "mean_pair_mult": round(statistics.mean(r["mean_pair_mult"] for r in group), 2),
            "max_pair_mult": round(statistics.mean(r["max_pair_mult"] for r in group), 2),
            "pair_hub_frac": round(statistics.mean(r["pair_hub_frac"] for r in group), 4),
            "mean_repeat_frac": round(statistics.mean(r["mean_repeat_frac"] for r in group), 4),
            "mean_node_repeat_frac": round(statistics.mean(r["mean_node_repeat_frac"] for r in group), 4),
            "max_node_repeat_frac": round(statistics.mean(r["max_node_repeat_frac"] for r in group), 4),
            "repeat_hub_frac": round(statistics.mean(r["repeat_hub_frac"] for r in group), 4),
            "mean_shared_node_frac": round(statistics.mean(r["mean_shared_node_frac"] for r in group), 4),
            "min_shared_node_frac": round(statistics.mean(r["min_shared_node_frac"] for r in group), 4),
            "max_shared_node_frac": round(statistics.mean(r["max_shared_node_frac"] for r in group), 4),
            "mean_family_shared_node_frac": round(statistics.mean(r["mean_family_shared_node_frac"] for r in group), 4),
            "mean_pair_jaccard": round(statistics.mean(r["mean_pair_jaccard"] for r in group), 4),
            "mean_exon_jaccard": round(statistics.mean(r["mean_exon_jaccard"] for r in group), 4),
            "mean_intron_jaccard": round(statistics.mean(r["mean_intron_jaccard"] for r in group), 4),
            "mean_max_unshared_run": round(statistics.mean(r["mean_max_unshared_run"] for r in group), 4),
            "mean_n_exons": round(statistics.mean(r["mean_n_exons"] for r in group), 2),
            "n_protein_families": round(statistics.mean(r["n_protein_families"] for r in group), 2),
        }

    ep_impure_rows = [r for r in features if r["ep_impure"]]
    ep_pure_rows = [r for r in features if not r["ep_impure"]]
    dna_fp_rows = [r for r in features if r["any_dna_fp"]]
    clean_rows = [r for r in features if not r["any_dna_fp"]]

    summary = {
        "total_families": len(features),
        "ep_impure": summarize(ep_impure_rows, "ep_impure"),
        "ep_pure": summarize(ep_pure_rows, "ep_pure"),
        "dna_fp": summarize(dna_fp_rows, "dna_fp"),
        "dna_clean": summarize(clean_rows, "dna_clean"),
    }

    # simple rule sweeps on held-out flag sets (using the whole catalog for discovery)
    def apply_rule(pred):
        flagged = [r for r in features if pred(r)]
        removed_fp = sum(1 for r in flagged if r["any_dna_fp"])
        removed_ep = sum(1 for r in flagged if r["ep_impure"])
        removed_real = len(flagged) - removed_fp - removed_ep
        return dict(flagged=len(flagged), removed_dna_fp=removed_fp,
                    removed_ep_impure=removed_ep, collateral_real=removed_real)

    rules = {
        "max_pair_mult>=200": apply_rule(lambda r: r["max_pair_mult"] >= 200),
        "mean_pair_mult>=80": apply_rule(lambda r: r["mean_pair_mult"] >= 80),
        "pair_hub_frac>=0.5": apply_rule(lambda r: r["pair_hub_frac"] >= 0.5),
        "n_protein_families>=2": apply_rule(lambda r: r["n_protein_families"] >= 2),
        "n_genes>n_distinct_loci": apply_rule(lambda r: r["n_genes"] > r["n_distinct_loci"]),
        "max_recip_overlap>=0.5": apply_rule(lambda r: r["max_recip_overlap"] >= 0.5),
        "mean_recip_overlap>=0.3": apply_rule(lambda r: r["mean_recip_overlap"] >= 0.3),
        "strand_majority<1.0 AND n_chrom==1": apply_rule(
            lambda r: r["strand_majority"] < 1.0 and r["n_chrom"] == 1),
        "mean_exon_jaccard<=0.5": apply_rule(lambda r: r["mean_exon_jaccard"] <= 0.5),
        "mean_intron_jaccard<=0.3": apply_rule(lambda r: r["mean_intron_jaccard"] <= 0.3),
        "mean_max_unshared_run>=5": apply_rule(lambda r: r["mean_max_unshared_run"] >= 5),
        "mean_exon_jaccard<=0.6 AND mean_max_unshared_run>=3": apply_rule(
            lambda r: r["mean_exon_jaccard"] <= 0.6 and r["mean_max_unshared_run"] >= 3),
        "mean_family_shared_node_frac<=0.5": apply_rule(
            lambda r: r["mean_family_shared_node_frac"] <= 0.5),
        "mean_pair_jaccard<=0.5": apply_rule(lambda r: r["mean_pair_jaccard"] <= 0.5),
        "mean_pair_jaccard<=0.7 AND mean_max_unshared_run>=3": apply_rule(
            lambda r: r["mean_pair_jaccard"] <= 0.7 and r["mean_max_unshared_run"] >= 3),
    }
    summary["rule_sweeps"] = rules

    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    print(f"[write] {OUT_JSON}", flush=True)

    print("\n=== SUMMARY ===")
    print(f"Total families: {len(features)}")
    print(f"E_p-impure: {len(ep_impure_rows)} ; E_p-pure: {len(ep_pure_rows)}")
    print(f"DNA-confirmed FP: {len(dna_fp_rows)}")
    print("\nE_p-impure vs pure means:")
    for k in summary["ep_pure"]:
        if k == "n":
            continue
        print(f"  {k:22s} pure={summary['ep_pure'][k]:>8}  impure={summary['ep_impure'][k]:>8}")
    print("\nRule sweeps:")
    for name, res in rules.items():
        print(f"  {name:35s} flagged={res['flagged']:3d}  FP={res['removed_dna_fp']}  "
              f"EP={res['removed_ep_impure']:2d}  collateral_real={res['collateral_real']:3d}")


if __name__ == "__main__":
    characterize()
