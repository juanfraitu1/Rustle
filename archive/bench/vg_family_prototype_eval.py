#!/usr/bin/env python3
"""vg_family_prototype_eval.py — evaluate the VG-native family catalog.

Loads bench/vg_family_prototype.tsv and runs the same family-level P/R
machinery used for the shipped family_rna_refine.tsv:

  * protein E_p block purity
  * DNA-loose cDNA component recovery
  * diploid-DNA oracle precision/recall
  * residual FP roster (multifam / oversize / allele-as-copy)

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype_eval.py
"""
import argparse
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
from collections import defaultdict
from itertools import combinations

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)

import graph_def_refine_sweep as SW
import genome_family_def as G
import rna_only_edge_oracle as RN
from family_level_pr_current import (
    build_ctx,
    oracle_multicopy_genes,
    dna_block_metrics,
    impure_block_roster,
    _distinct_fp_blocks,
    _block_oracle_flags,
    relabel_recovered_genes,
)

PROTOTYPE_TSV = os.path.join(BENCH, "vg_family_prototype.tsv")
OUT_JSON = os.path.join(BENCH, "vg_family_prototype_eval.json")
OUT_TSV = os.path.join(BENCH, "vg_family_prototype_eval.tsv")
CURRENT_JSON = os.path.join(BENCH, "family_level_pr_current.json")


def load_prototype_blocks():
    """Load prototype catalog as list of sets (one set per family)."""
    blocks = defaultdict(set)
    with open(PROTOTYPE_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            blocks[row["fam_id"]].add(row["member"])
    order = sorted(blocks.keys(), key=lambda k: int(k))
    return [blocks[k] for k in order]


def oracle_pr(res, mc_oracle, extra_credit=frozenset()):
    """Scalar P/R from RN.oracle_residuals, with optional relabel credit."""
    orac = res["orac_blocks"]
    n_allele = len(res["allele"])
    n_over = len(res["oversize"])
    n_multi = len(res["multifam"])
    flags = n_allele + n_over + n_multi
    distinct_fp = _distinct_fp_blocks(res)
    base_recovered = set(res["recovered"])
    credit_new = (set(extra_credit) & mc_oracle) - base_recovered
    recovered_set = base_recovered | credit_new
    recov_mc = len(recovered_set & mc_oracle)
    return dict(
        oracle_mapped_families=orac,
        n_allele=n_allele,
        n_oversize=n_over,
        n_multifam=n_multi,
        flag_instances=flags,
        distinct_fp_blocks=distinct_fp,
        precision_oracle_taskformula=round(1 - flags / orac, 4) if orac else None,
        precision_oracle_dedup=round(1 - distinct_fp / orac, 4) if orac else None,
        oracle_genes_recovered_multicopy=recov_mc,
        oracle_multicopy_genes=len(mc_oracle),
        recall_oracle=round(recov_mc / len(mc_oracle), 4) if mc_oracle else None,
        relabel_credited_genes=sorted(credit_new),
    )


def main():
    global PROTOTYPE_TSV, OUT_JSON, OUT_TSV

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=PROTOTYPE_TSV,
                        help="prototype catalog TSV to evaluate (default: vg_family_prototype.tsv)")
    parser.add_argument("--output-prefix", default=None,
                        help="prefix for output JSON/TSV (default: same base as input)")
    args = parser.parse_args()

    PROTOTYPE_TSV = args.input
    if args.output_prefix:
        OUT_JSON = args.output_prefix + ".json"
        OUT_TSV = args.output_prefix + ".tsv"
    else:
        base, _ = os.path.splitext(PROTOTYPE_TSV)
        OUT_JSON = base + "_eval.json"
        OUT_TSV = base + "_eval.tsv"

    if not os.path.exists(PROTOTYPE_TSV):
        print(f"[!] prototype catalog not found: {PROTOTYPE_TSV}", flush=True)
        print("    run bench/vg_family_prototype.py first", flush=True)
        return 1

    print("[load] ctx (graph / projection / truth / oracle) ...", flush=True)
    ctx, *_ = build_ctx()
    genes = ctx["genes"]
    gene_of_dn = ctx["gene_of_dn"]
    g2rows = ctx["g2rows"]

    print("[load] prototype catalog ...", flush=True)
    proto_blocks = load_prototype_blocks()
    print(f"       {len(proto_blocks)} raw families", flush=True)

    # restrict to the loci that the shipped eval ctx knows about (same gene/interval projection)
    eval_nodes = set(ctx["genes"])
    proto_blocks = [b & eval_nodes for b in proto_blocks]
    proto_blocks = [b for b in proto_blocks if b]
    print(f"       {len(proto_blocks)} raw families after intersecting eval ctx", flush=True)

    proto_fams = SW.filter_multicopy(proto_blocks, genes)
    print(f"       {len(proto_fams)} multi-copy families (>=2 distinct loci)", flush=True)

    # ---- same metrics as family_level_pr_current ----
    proto_eval = SW.eval_partition("vg_family_prototype", proto_blocks, ctx)
    proto_res = RN.oracle_residuals(proto_fams, gene_of_dn, genes, g2rows, G, margin=SW.DIPLOID_MARGIN)
    mc_oracle = oracle_multicopy_genes(ctx)

    relabel_credit, relabel_note, relabel_detail = relabel_recovered_genes(proto_blocks, genes)
    orc = oracle_pr(proto_res, mc_oracle, extra_credit=relabel_credit)

    # E_p purity
    ep = dict(
        total_blocks=proto_eval["n_families"],
        impure_blocks=proto_eval["impure_blocks"],
        pure_blocks=proto_eval["n_families"] - proto_eval["impure_blocks"],
        precision_Ep=round(1 - proto_eval["impure_blocks"] / proto_eval["n_families"], 4)
        if proto_eval["n_families"] else None,
    )

    # DNA-loose block metrics
    dnab = dna_block_metrics(proto_fams, ctx)
    dna = dict(
        pair_projection_recall_real_cdna=proto_eval["tp_retention"],
        pair_projection_precision_nonovermerge=proto_eval["genuine_precision"],
        block_overmerge_precision=dnab["block_overmerge_precision"],
        component_recovery_recall=dnab["component_recovery_recall"],
        n_multigene_blocks=dnab["n_multigene_blocks"],
        n_blocks_with_genuine_pair=dnab["n_blocks_with_genuine_pair"],
        n_real_dna_components=dnab["n_real_dna_components"],
        n_components_recovered=dnab["n_components_recovered"],
    )

    # residual FP roster
    fp_multifam = [dict(genes=e[0], distinct_loci=e[1]) for e in
                   sorted(proto_res["multifam"], key=lambda x: (-x[1], x[0][0]))]
    fp_oversize = [dict(genes=e[0], distinct_loci=e[1], diploid_CN=e[2], ratio=e[3]) for e in
                   sorted(proto_res["oversize"], key=lambda x: -x[3])]
    fp_allele = [dict(genes=e[0], distinct_loci=e[1], asm_hapCN=e[2]) for e in proto_res["allele"]]

    dna_named_genes = set()
    for e in proto_res["multifam"] + proto_res["oversize"] + proto_res["allele"]:
        dna_named_genes.update(e[0])
    impure_roster = impure_block_roster(proto_fams, ctx)
    impure_not_dna = [b for b in impure_roster if not (set(b["genes"]) & dna_named_genes)]

    # per-family TSV
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "family_id", "n_loci_distinct", "n_genes", "genes",
            "ep_impure", "n_protein_families", "oracle_genes", "diploid_CN",
            "fp_multifam", "fp_oversize", "fp_allele",
        ])
        for bi, b in enumerate(proto_fams):
            gs = sorted({gene_of_dn[dn] for dn in b if dn in gene_of_dn})
            fam_of = defaultdict(list)
            for g in gs:
                if g in ctx["g2f"] and ctx["g2f"][g] not in ctx["mega"]:
                    fam_of[ctx["g2f"][g]].append(g)
            fl = _block_oracle_flags(b, ctx)
            w.writerow([
                bi,
                G.distinct_loci(b, genes),
                len(gs),
                ",".join(gs),
                int(bi in {x["block_index"] for x in impure_roster}),
                len(fam_of),
                ",".join(fl["oracle"]),
                fl["dip"],
                fl["multifam"],
                fl["oversize"],
                fl["allele"],
            ])
    print(f"[write] {OUT_TSV}", flush=True)

    # summary JSON
    out = dict(
        prototype_catalog=PROTOTYPE_TSV,
        n_raw_families=len(proto_blocks),
        n_multicopy_families=len(proto_fams),
        truth1_protein_Ep=ep,
        truth2_dna_loose=dna,
        truth3_diploid_oracle=orc,
        graph_level=dict(
            tp_retention=proto_eval["tp_retention"],
            overmerge_cut=proto_eval["overmerge_cut"],
            truthbar_retention=proto_eval["truthbar_retention"],
            genuine_precision=proto_eval["genuine_precision"],
            ge_genuine_cut=proto_eval["ge_genuine_cut"],
            ge_real_kept=proto_eval["ge_real_kept"],
        ),
        residual_fp_roster=dict(
            multifam=fp_multifam,
            oversize=fp_oversize,
            allele=fp_allele,
            ep_impure_not_dna_named=impure_not_dna,
            ep_impure_total=len(impure_roster),
        ),
        gene_projection_relabel=dict(
            note=relabel_note,
            credited_genes=sorted(relabel_credit & mc_oracle),
            n_newly_credited=len(orc["relabel_credited_genes"]),
            recall_without_relabel=round(
                (orc["oracle_genes_recovered_multicopy"] - len(orc["relabel_credited_genes"]))
                / len(mc_oracle), 4) if mc_oracle else None,
            recall_with_relabel=orc["recall_oracle"],
        ),
        determinism=dict(seed=SW.SEED, gamma=G.GAMMA,
                         pythonhashseed=os.environ.get("PYTHONHASHSEED")),
    )
    json.dump(out, open(OUT_JSON, "w"), indent=2, sort_keys=True)
    print(f"[write] {OUT_JSON}", flush=True)

    # print summary
    P = print
    P("\n================ VG-NATIVE PROTOTYPE :: FAMILY-LEVEL P/R ================")
    P(f"Prototype: {len(proto_blocks)} raw families -> {len(proto_fams)} multi-copy families")
    P("\n---- TRUTH 1 :: PROTEIN E_p PURITY ----")
    P(f"  pure {ep['pure_blocks']}/{ep['total_blocks']}  impure {ep['impure_blocks']}  ->  P_Ep = {ep['precision_Ep']}")
    P("\n---- TRUTH 2 :: DNA-loose cDNA ----")
    P(f"  pair-recall(real cDNA)={dna['pair_projection_recall_real_cdna']}  "
      f"block-overmerge-precision={dna['block_overmerge_precision']}  "
      f"component-recovery-recall={dna['component_recovery_recall']} "
      f"({dna['n_components_recovered']}/{dna['n_real_dna_components']})")
    P("\n---- TRUTH 3 :: DIPLOID DNA ORACLE ----")
    P(f"  oracle-mapped {orc['oracle_mapped_families']}  "
      f"FP(allele/oversize/multifam)={orc['n_allele']}/{orc['n_oversize']}/{orc['n_multifam']} "
      f"(flags {orc['flag_instances']}, distinct blocks {orc['distinct_fp_blocks']})")
    P(f"  P_oracle(task)={orc['precision_oracle_taskformula']}  "
      f"P_oracle(dedup)={orc['precision_oracle_dedup']}")
    P(f"  R_oracle = {orc['oracle_genes_recovered_multicopy']}/{orc['oracle_multicopy_genes']} = {orc['recall_oracle']}  "
      f"(+{len(orc['relabel_credited_genes'])} relabel credit)")
    P("\n---- RESIDUAL FP ROSTER ----")
    P(f"  multifam: {len(fp_multifam)}")
    for e in fp_multifam[:10]:
        P(f"      {'+'.join(e['genes'])} (dl={e['distinct_loci']})")
    P(f"  oversize: {len(fp_oversize)}")
    for e in fp_oversize[:10]:
        P(f"      {'+'.join(e['genes'])} (dl={e['distinct_loci']}, dipCN={e['diploid_CN']}, {e['ratio']}x)")
    P(f"  allele-as-copy: {len(fp_allele)}")
    for e in fp_allele[:10]:
        P(f"      {'+'.join(e['genes'])} (dl={e['distinct_loci']}, hapCN={e['asm_hapCN']})")
    P(f"  E_p-impure not DNA-named: {len(impure_not_dna)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
