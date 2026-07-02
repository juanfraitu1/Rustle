#!/usr/bin/env python
"""FAMILY-LEVEL (block-level) precision/recall of the CURRENT DEFAULT family
definition = the RNA-only recall-preserving refined oracle.

CURRENT DEFAULT catalog  : bench/family_rna_refine.tsv (607 multi-copy families,
                           rule KEEP iff core_recip>=0.19 AND aln_frac>=0.24 ->
                           gamma-quasi-clique(0.20) -> allele-demote; committed cc472a4).
We LOAD that catalog directly (member_dn per family) -- we DO NOT re-derive it.

We reuse the SHIPPED eval machinery so the numbers match the committed catalogs:
    graph_def_refine_sweep.eval_partition   (Ep-impurity + diploid-oracle tallies)
    rna_only_edge_oracle.oracle_residuals   (named residual-FP roster)

THREE clearly-labelled family/BLOCK-level truths:
  (1) PROTEIN E_p purity     : block PURE iff its member genes span <=1 non-mega
                               protein family E_p.  P_Ep = pure_blocks/total_blocks.
                               (conservative: E_p counts divergent-paralog "impurity"
                                that may be real.)
  (2) DNA-loose (cDNA)       : gene-pair truth in_dna_loose (family_er_pr.tsv, same
                               values as family_def_dna_pr_edges.tsv). Pair-projection
                               tp_retention = recall of real cDNA families; block over-
                               merge precision = blocks carrying no GENUINE (over-merge)
                               pair; plus block-level component recovery.
  (3) DIPLOID DNA ORACLE     : the only assembly-independent scalar (asm_hapCN, mat+pat).
                               P_oracle = 1-(multifam+oversize+allele)/oracle_mapped;
                               R_oracle = oracle_genes_recovered/oracle_multicopy_genes.

Legacy comparison = shipped gamma-quasi-clique(core_recip>=0.13, gamma 0.20), recomputed
here AND cross-checked vs the committed graph_def_refine_sweep.json baseline.

Deterministic (PYTHONHASHSEED=0).
Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_level_pr_current.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
from collections import defaultdict
from itertools import combinations

import networkx as nx
from networkx.algorithms.community import louvain_communities

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import family_er_pr as FP
import genome_family_def as G
import graph_def_refine_sweep as SW
import rna_only_edge_oracle as RN

CATALOG_TSV = os.path.join(BENCH, "family_rna_refine.tsv")
OUT_TSV = os.path.join(BENCH, "family_level_pr_current.tsv")
OUT_JSON = os.path.join(BENCH, "family_level_pr_current.json")
SWEEP_JSON = os.path.join(BENCH, "graph_def_refine_sweep.json")
CAND_TSV = os.path.join(BENCH, "candidate_generation_recall.tsv")


# --------------------------------------------------------------------------- ctx
def build_ctx():
    """Reconstruct the exact ctx graph_def_refine_sweep.main() builds."""
    meta = FP.load_meta(); annot = FP.load_annot(); gene_of = FP.gene_of_factory(annot)
    raw_fams = FP.load_raw_families(); edge_pairs = FP.load_edges()
    genes, gene_of_dn, *_ = FP.build_genes_dict(raw_fams, meta, gene_of)
    ew = SW.load_weighted_edges()
    lab = SW.load_pair_labels()
    g2f, mega = FP.load_prfam()
    g2rows = SW.load_oracle()

    all_nodes = set()
    for f in raw_fams:
        all_nodes.update(f)
    Gr = nx.Graph(); Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b, weight=ew.get(frozenset((a, b)), 0.13))

    REAL = {k for k, v in lab.items() if v == "REAL_cdna"}
    GEN = {k for k, v in lab.items() if v == "GENUINE"}
    TB = {k for k, v in lab.items() if v == "TRUTHBAR"}

    ctx = dict(genes=genes, gene_of_dn=gene_of_dn, lab=lab, g2f=g2f, mega=mega,
               g2rows=g2rows, REAL=REAL, GEN=GEN, TB=TB, G=Gr, ew=ew)
    return ctx, raw_fams, edge_pairs


def legacy_gamma_blocks(raw_fams, edge_pairs):
    """Shipped gamma-quasi-clique(0.20) on the raw core_recip>=0.13 components."""
    adj = defaultdict(set); fam_of = {}
    for k, m in enumerate(raw_fams):
        for i in m:
            fam_of[i] = k
    for u, v in edge_pairs:
        if u in fam_of and fam_of.get(u) == fam_of.get(v):
            adj[u].add(v); adj[v].add(u)
    out = []
    for m in raw_fams:
        for block in G._refine_component(set(m), adj, G.GAMMA, louvain_communities, SW.SEED):
            out.append(set(block))
    return out


def load_catalog_blocks():
    """CURRENT DEFAULT catalog: family_id -> set(member_dn).  Also keep per-family
    metadata (dominant_gene, member_gene per dn) for FP naming."""
    blocks = defaultdict(set)
    fam_dom = {}
    member_gene = {}
    with open(CATALOG_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            fid = row["family_id"]
            blocks[fid].add(row["member_dn"])
            fam_dom[fid] = row["dominant_gene"]
            member_gene[row["member_dn"]] = row["member_gene"]
    order = sorted(blocks, key=lambda f: int(f))
    block_list = [blocks[f] for f in order]
    return block_list, order, fam_dom, member_gene


# --------------------------------------------------------------------------- gene-projection relabel
def relabel_recovered_genes(cat_blocks, genes):
    """GENE-PROJECTION RELABEL (recall-metric matching; MEASUREMENT-ONLY, monotone).

    The diploid-oracle recall credits an oracle gene as 'recovered' iff the max-overlap gene_of
    projection of some de-novo locus in a MULTI-COPY catalog family equals the oracle gene NAME.
    In a segdup cluster a single de-novo copy overlaps several paralogous genes, and the max-overlap
    projection labels it after a NEIGHBOUR gene (e.g. the ARL17A / LRRC37A cluster) -- so an oracle
    gene whose OWN copy IS grouped into a recovered multi-copy family is not credited by NAME. The
    catalog GROUPING is correct; only the gene-name projection/measurement misses it.

    This restores the missing credit by GENOMIC / copy linkage instead of gene NAME. An oracle gene is
    credited iff its OWN de-novo copy participates in a POA-passing cross-copy homology edge
    (candidate_generation_recall.tsv: touches_gene_copy==1 & passes_scoring==1) at least one endpoint
    of which is a member of a MULTI-COPY catalog family -- i.e. the gene's own copy is demonstrably
    linked into a recovered multi-copy family, only mislabelled. This is exactly the "the catalog
    grouping is correct; the projection mislabels" case, and it excludes the genuine misses (a gene
    whose own copy FAILS the POA gate, or whose linked loci were NOT retained in the catalog).

    Deliberately NARROW (only the diagnosed candidate rows can ever be credited) and MEASUREMENT-ONLY:
    it augments the recall NUMERATOR, never gene_of / gene_of_dn, the block projection, the FP flags or
    the family grouping -- so precision and the family_rna_refine catalog stay byte-identical. If the
    diagnosis TSV is absent it credits nothing (falls back to the pure name-projection recall).
    Deterministic (sorted / set membership only).  Returns (credited_genes, note, detail)."""
    mc_loci = set()
    for b in cat_blocks:
        if G.distinct_loci(b, genes) >= 2:
            mc_loci |= set(b)
    if not os.path.exists(CAND_TSV):
        return set(), "candidate_generation_recall.tsv absent -> name-projection recall only (no relabel)", {}
    credited = set()
    detail = {}
    with open(CAND_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row.get("touches_gene_copy") == "1" and row.get("passes_scoring") == "1":
                ci, cj = row["ci_member"], row["cj_member"]
                if ci in mc_loci or cj in mc_loci:
                    g = row["gene"]
                    credited.add(g)
                    detail.setdefault(g, []).append(dict(
                        edge=[ci, cj], core_recip=row.get("core_recip"),
                        endpoint_in_multicopy_family=[e for e in (ci, cj) if e in mc_loci]))
    return credited, "ok", {g: detail[g] for g in sorted(detail)}


# --------------------------------------------------------------------------- oracle multicopy denominator
def oracle_multicopy_genes(ctx):
    """distinct oracle REP genes whose diploid-DNA class is multi_copy (recall denom)."""
    mc = set()
    for g, rows in ctx["g2rows"].items():
        if any(r["cls"] == "multi_copy" for r in rows):
            mc.add(g)
    return mc


# --------------------------------------------------------------------------- block-level DNA-loose
def dna_block_metrics(fams, ctx):
    """BLOCK-level DNA-loose over-merge purity + component recovery recall.
       real DNA family = connected component (>=2 genes) of the in_dna_loose gene graph
       restricted to genes present in the catalog's multi-copy blocks."""
    gene_of_dn = ctx["gene_of_dn"]; lab = ctx["lab"]
    REAL = ctx["REAL"]; GEN = ctx["GEN"]; TB = ctx["TB"]

    # gene -> set(block index); block -> gene set
    block_genes = []
    gene_blocks = defaultdict(set)
    for bi, b in enumerate(fams):
        gs = sorted({gene_of_dn[dn] for dn in b if dn in gene_of_dn})
        block_genes.append(gs)
        for g in gs:
            gene_blocks[g].add(bi)

    catalog_genes = set(gene_blocks)

    # --- block over-merge precision (block carries a GENUINE cross-gene pair?) ---
    n_multigene = 0
    n_with_labeled = 0
    blocks_with_genuine = []      # (block_index, dominant, genuine_pairs)
    for bi, gs in enumerate(block_genes):
        if len(gs) < 2:
            continue
        n_multigene += 1
        pairs = [frozenset(p) for p in combinations(gs, 2)]
        labeled = [p for p in pairs if p in lab]
        if labeled:
            n_with_labeled += 1
        gen = [tuple(sorted(p)) for p in labeled if lab[p] == "GENUINE"]
        if gen:
            blocks_with_genuine.append((bi, gen))
    block_precision = 1 - len(blocks_with_genuine) / n_with_labeled if n_with_labeled else 1.0

    # --- component recovery recall (real DNA families fully co-membered) ---
    H = nx.Graph()
    H.add_nodes_from(catalog_genes)
    for k in REAL:
        a, b = tuple(k)
        if a in catalog_genes and b in catalog_genes:
            H.add_edge(a, b)
    comps = [c for c in nx.connected_components(H) if len(c) >= 2]
    recovered_comp = 0
    for c in comps:
        blks = set()
        for g in c:
            blks |= gene_blocks[g]
        # recovered iff every gene of the component sits in ONE common block
        common = None
        for g in c:
            common = gene_blocks[g] if common is None else (common & gene_blocks[g])
        if common:
            recovered_comp += 1
    comp_recall = recovered_comp / len(comps) if comps else 0.0

    return dict(n_multigene_blocks=n_multigene, n_blocks_with_labeled_crossgene=n_with_labeled,
                n_blocks_with_genuine_pair=len(blocks_with_genuine),
                block_overmerge_precision=round(block_precision, 4),
                n_real_dna_components=len(comps), n_components_recovered=recovered_comp,
                component_recovery_recall=round(comp_recall, 4),
                _blocks_with_genuine=blocks_with_genuine, _block_genes=block_genes)


# --------------------------------------------------------------------------- Ep-impure enumeration
def impure_block_roster(fams, ctx):
    """Enumerate the Ep-impure blocks (>=2 distinct NON-mega protein families) with
       their member genes -- exactly the (b) predicate in eval_partition."""
    gene_of_dn = ctx["gene_of_dn"]; g2f = ctx["g2f"]; mega = ctx["mega"]
    roster = []
    for bi, b in enumerate(fams):
        gs = sorted({gene_of_dn[dn] for dn in b if dn in gene_of_dn})
        fam_of = defaultdict(list)
        for g in gs:
            if g in g2f and g2f[g] not in mega:
                fam_of[g2f[g]].append(g)
        if len(fam_of) >= 2:
            roster.append(dict(block_index=bi, n_genes=len(gs),
                               distinct_loci=G.distinct_loci(b, ctx["genes"]),
                               n_protein_families=len(fam_of),
                               genes=gs,
                               protein_families={fid: sorted(v) for fid, v in fam_of.items()}))
    return roster


# --------------------------------------------------------------------------- main
def main():
    print("[load] ctx (graph / projection / truth / oracle) ...", flush=True)
    ctx, raw_fams, edge_pairs = build_ctx()
    genes = ctx["genes"]; gene_of_dn = ctx["gene_of_dn"]; g2rows = ctx["g2rows"]

    cat_blocks, fam_ids, fam_dom, member_gene = load_catalog_blocks()
    print(f"       CURRENT catalog family_rna_refine.tsv: {len(cat_blocks)} families")

    # ----- run the SHIPPED eval on the CURRENT catalog and the LEGACY gamma -----
    cur = SW.eval_partition("family_rna_refine_DEFAULT", cat_blocks, ctx)
    leg_blocks = legacy_gamma_blocks(raw_fams, edge_pairs)
    leg = SW.eval_partition("gamma_shipped_0.13->0.20_LEGACY", leg_blocks, ctx)

    cur_fams = SW.filter_multicopy(cat_blocks, genes)
    leg_fams = SW.filter_multicopy(leg_blocks, genes)
    cur_res = RN.oracle_residuals(cur_fams, gene_of_dn, genes, g2rows, G, margin=SW.DIPLOID_MARGIN)
    leg_res = RN.oracle_residuals(leg_fams, gene_of_dn, genes, g2rows, G, margin=SW.DIPLOID_MARGIN)

    # cross-check legacy vs committed baseline
    committed = json.load(open(SWEEP_JSON))["baseline_gamma_0_20"]
    xcheck = dict(
        legacy_n_families=(leg["n_families"], committed["n_families"]),
        legacy_impure=(leg["impure_blocks"], committed["impure_blocks"]),
        legacy_oracle_blocks=(leg["oracle_blocks"], committed["oracle_blocks"]),
        legacy_oracle_residual_total=(
            leg["oracle_multifam"] + leg["oracle_oversize_diploid"] + leg["oracle_allele_as_copy"],
            committed["oracle_multifam"] + committed["oracle_oversize_diploid"] + committed["oracle_allele_as_copy"]),
        legacy_recovered=(leg["oracle_recovered_genes"], committed["oracle_recovered_genes"]),
    )
    xcheck_ok = all(a == b for a, b in xcheck.values())

    mc_oracle = oracle_multicopy_genes(ctx)
    n_mc_oracle = len(mc_oracle)

    # GENE-PROJECTION RELABEL (recall numerator only; see relabel_recovered_genes). Credits oracle
    # genes whose own copy IS grouped into a recovered multi-copy family but is mislabelled by the
    # max-overlap projection (segdup neighbour). NEVER touches gene_of / block projection / FP flags.
    relabel_credit, relabel_note, relabel_detail = relabel_recovered_genes(cat_blocks, genes)
    relabel_credit_mc = relabel_credit & mc_oracle

    # ================= TRUTH 1 : PROTEIN E_p purity =================
    ep_cur = dict(total_blocks=cur["n_families"], impure_blocks=cur["impure_blocks"],
                  pure_blocks=cur["n_families"] - cur["impure_blocks"],
                  precision_Ep=round(1 - cur["impure_blocks"] / cur["n_families"], 4))
    ep_leg = dict(total_blocks=leg["n_families"], impure_blocks=leg["impure_blocks"],
                  pure_blocks=leg["n_families"] - leg["impure_blocks"],
                  precision_Ep=round(1 - leg["impure_blocks"] / leg["n_families"], 4))

    # ================= TRUTH 2 : DNA-loose (cDNA) =================
    dnab_cur = dna_block_metrics(cur_fams, ctx)
    dnab_leg = dna_block_metrics(leg_fams, ctx)
    dna_cur = dict(
        pair_projection_recall_real_cdna=cur["tp_retention"],       # shipped
        pair_projection_precision_nonovermerge=cur["genuine_precision"],  # shipped
        block_overmerge_precision=dnab_cur["block_overmerge_precision"],
        component_recovery_recall=dnab_cur["component_recovery_recall"],
        n_multigene_blocks=dnab_cur["n_multigene_blocks"],
        n_blocks_with_genuine_pair=dnab_cur["n_blocks_with_genuine_pair"],
        n_real_dna_components=dnab_cur["n_real_dna_components"],
        n_components_recovered=dnab_cur["n_components_recovered"])
    dna_leg = dict(
        pair_projection_recall_real_cdna=leg["tp_retention"],
        pair_projection_precision_nonovermerge=leg["genuine_precision"],
        block_overmerge_precision=dnab_leg["block_overmerge_precision"],
        component_recovery_recall=dnab_leg["component_recovery_recall"],
        n_multigene_blocks=dnab_leg["n_multigene_blocks"],
        n_blocks_with_genuine_pair=dnab_leg["n_blocks_with_genuine_pair"],
        n_real_dna_components=dnab_leg["n_real_dna_components"],
        n_components_recovered=dnab_leg["n_components_recovered"])

    # ================= TRUTH 3 : DIPLOID DNA ORACLE =================
    def oracle_pr(res, recovered_field, extra_credit=frozenset()):
        orac = res["orac_blocks"]
        n_allele = len(res["allele"]); n_over = len(res["oversize"]); n_multi = len(res["multifam"])
        flags = n_allele + n_over + n_multi
        # de-duplicated distinct offending blocks (a block can be both oversize+multifam)
        distinct_fp = _distinct_fp_blocks(res)
        # RECALL NUMERATOR ONLY: augment the name-projection recovered set with the gene-projection
        # relabel credit (extra_credit). This changes neither orac_blocks, the FP flags (allele/
        # oversize/multifam), nor distinct_fp -> precision (P_oracle) is byte-identical; only the
        # recall numerator rises. extra_credit is empty for the legacy comparison.
        base_recovered = set(res["recovered"])
        credit_new = (set(extra_credit) & mc_oracle) - base_recovered
        recovered_set = base_recovered | credit_new
        recov = len(recovered_set)
        recov_mc = len(recovered_set & mc_oracle)
        return dict(oracle_mapped_families=orac,
                    n_allele=n_allele, n_oversize=n_over, n_multifam=n_multi,
                    flag_instances=flags, distinct_fp_blocks=distinct_fp,
                    precision_oracle_taskformula=round(1 - flags / orac, 4) if orac else None,
                    precision_oracle_dedup=round(1 - distinct_fp / orac, 4) if orac else None,
                    oracle_genes_recovered=recov,
                    oracle_genes_recovered_multicopy=recov_mc,
                    oracle_multicopy_genes=n_mc_oracle,
                    relabel_credited_genes=sorted(credit_new),
                    recall_oracle=round(recov_mc / n_mc_oracle, 4) if n_mc_oracle else None)
    orc_cur = oracle_pr(cur_res, "cur", extra_credit=relabel_credit_mc)
    orc_leg = oracle_pr(leg_res, "leg")

    # ================= ENUMERATE RESIDUAL FP SET (current catalog) =================
    fp_multifam = [dict(genes=e[0], distinct_loci=e[1]) for e in
                   sorted(cur_res["multifam"], key=lambda x: (-x[1], x[0][0]))]
    fp_oversize = [dict(genes=e[0], distinct_loci=e[1], diploid_CN=e[2], ratio=e[3]) for e in
                   sorted(cur_res["oversize"], key=lambda x: -x[3])]
    fp_allele = [dict(genes=e[0], distinct_loci=e[1], asm_hapCN=e[2]) for e in cur_res["allele"]]

    # E_p-impure blocks NOT already named by the DNA oracle
    dna_named_genes = set()
    for e in cur_res["multifam"] + cur_res["oversize"] + cur_res["allele"]:
        dna_named_genes.update(e[0])
    impure_roster = impure_block_roster(cur_fams, ctx)
    impure_not_dna = [b for b in impure_roster
                      if not (set(b["genes"]) & dna_named_genes)]

    # ---- per-family purity flags TSV ----
    write_family_tsv(cat_blocks, fam_ids, fam_dom, ctx, cur_res, impure_roster)

    # ================= REPORT =================
    P = print
    P("\n================ FAMILY-LEVEL (block-level) P/R :: CURRENT DEFAULT ================")
    P(f"CURRENT  = family_rna_refine.tsv  ({cur['n_families']} multi-copy families)")
    P(f"LEGACY   = shipped gamma-quasi-clique(core_recip>=0.13, gamma 0.20)  ({leg['n_families']} families)")
    P(f"legacy vs committed graph_def_refine_sweep.json cross-check: {'OK' if xcheck_ok else 'MISMATCH'}  {xcheck}")

    P("\n---- TRUTH 1 :: PROTEIN E_p PURITY (block precision; conservative) ----")
    P(f"  CURRENT : pure {ep_cur['pure_blocks']}/{ep_cur['total_blocks']}  "
      f"impure {ep_cur['impure_blocks']}  ->  P_Ep = {ep_cur['precision_Ep']}")
    P(f"  LEGACY  : pure {ep_leg['pure_blocks']}/{ep_leg['total_blocks']}  "
      f"impure {ep_leg['impure_blocks']}  ->  P_Ep = {ep_leg['precision_Ep']}")

    P("\n---- TRUTH 2 :: DNA-loose cDNA (in_dna_loose) ----")
    P(f"  CURRENT : pair-recall(real cDNA)={dna_cur['pair_projection_recall_real_cdna']}  "
      f"block-overmerge-precision={dna_cur['block_overmerge_precision']} "
      f"(blocks w/ genuine pair {dna_cur['n_blocks_with_genuine_pair']}/{dnab_cur['n_blocks_with_labeled_crossgene']})  "
      f"component-recovery-recall={dna_cur['component_recovery_recall']} "
      f"({dna_cur['n_components_recovered']}/{dna_cur['n_real_dna_components']})")
    P(f"  LEGACY  : pair-recall(real cDNA)={dna_leg['pair_projection_recall_real_cdna']}  "
      f"block-overmerge-precision={dna_leg['block_overmerge_precision']} "
      f"(blocks w/ genuine pair {dna_leg['n_blocks_with_genuine_pair']}/{dnab_leg['n_blocks_with_labeled_crossgene']})  "
      f"component-recovery-recall={dna_leg['component_recovery_recall']} "
      f"({dna_leg['n_components_recovered']}/{dna_leg['n_real_dna_components']})")

    P("\n---- TRUTH 3 :: INDEPENDENT DIPLOID DNA ORACLE (assembly-independent) ----")
    P(f"  CURRENT : oracle-mapped {orc_cur['oracle_mapped_families']}  "
      f"FP(allele/oversize/multifam)={orc_cur['n_allele']}/{orc_cur['n_oversize']}/{orc_cur['n_multifam']} "
      f"(flags {orc_cur['flag_instances']}, distinct blocks {orc_cur['distinct_fp_blocks']})")
    P(f"            P_oracle(task) = {orc_cur['precision_oracle_taskformula']}   "
      f"P_oracle(dedup) = {orc_cur['precision_oracle_dedup']}")
    P(f"            R_oracle = {orc_cur['oracle_genes_recovered_multicopy']}/{orc_cur['oracle_multicopy_genes']} "
      f"= {orc_cur['recall_oracle']}  (n_recovered all-classes={orc_cur['oracle_genes_recovered']})")
    _rc = orc_cur["relabel_credited_genes"]
    P(f"            gene-projection RELABEL credited +{len(_rc)} {_rc}  "
      f"(own copy grouped into a recovered multi-copy family but max-overlap-mislabelled; "
      f"recall NUMERATOR only -> precision/grouping byte-identical)  [{relabel_note}]")
    P(f"  LEGACY  : oracle-mapped {orc_leg['oracle_mapped_families']}  "
      f"FP(allele/oversize/multifam)={orc_leg['n_allele']}/{orc_leg['n_oversize']}/{orc_leg['n_multifam']} "
      f"(flags {orc_leg['flag_instances']}, distinct blocks {orc_leg['distinct_fp_blocks']})")
    P(f"            P_oracle(task) = {orc_leg['precision_oracle_taskformula']}   "
      f"P_oracle(dedup) = {orc_leg['precision_oracle_dedup']}")
    P(f"            R_oracle = {orc_leg['oracle_genes_recovered_multicopy']}/{orc_leg['oracle_multicopy_genes']} "
      f"= {orc_leg['recall_oracle']}  (n_recovered all-classes={orc_leg['oracle_genes_recovered']})")

    P("\n================ RESIDUAL FALSE-POSITIVE ROSTER (CURRENT DEFAULT) ================")
    P(f"(a) DNA-CONFIRMED OVER-MERGES (multifam: block spans >=2 oracle genes):")
    for e in fp_multifam:
        P(f"      {'+'.join(e['genes'])}  (distinct_loci={e['distinct_loci']})")
    P(f"    OVER-SIZE (distinct_loci > 1.5x diploid DNA CN):")
    for e in fp_oversize:
        P(f"      {'+'.join(e['genes'])}  (dl={e['distinct_loci']}, dipCN={e['diploid_CN']}, {e['ratio']}x)")
    P(f"(b) RESIDUAL ALLELE-AS-COPY : {len(fp_allele)}"
      + ("" if not fp_allele else "  " + ", ".join('+'.join(e['genes']) for e in fp_allele)))
    P(f"(c) E_p-IMPURE blocks NOT named by DNA oracle : {len(impure_not_dna)} of {len(impure_roster)} impure")
    for b in sorted(impure_not_dna, key=lambda x: -x["n_genes"])[:40]:
        fams_str = "; ".join(f"{fid}:{','.join(v)}" for fid, v in b["protein_families"].items())
        P(f"      block#{b['block_index']} dl={b['distinct_loci']} nProtFam={b['n_protein_families']} :: {fams_str}")
    if len(impure_not_dna) > 40:
        P(f"      ... ({len(impure_not_dna)-40} more; full list in JSON)")

    # determinism echo
    det = dict(seed=SW.SEED, gamma=G.GAMMA, pythonhashseed=os.environ.get("PYTHONHASHSEED"))

    out = dict(
        current_catalog=CATALOG_TSV, n_families_current=cur["n_families"],
        n_families_legacy=leg["n_families"],
        legacy_crosscheck=dict(ok=xcheck_ok, pairs={k: list(v) for k, v in xcheck.items()}),
        truth1_protein_Ep=dict(current=ep_cur, legacy=ep_leg),
        truth2_dna_loose=dict(current=dna_cur, legacy=dna_leg),
        truth3_diploid_oracle=dict(current=orc_cur, legacy=orc_leg),
        gene_projection_relabel=dict(
            note=relabel_note,
            mechanism=("recall-numerator-only credit: an oracle multi-copy gene is credited if its OWN "
                       "de-novo copy sits in a POA-passing cross-copy edge (candidate_generation_recall: "
                       "touches_gene_copy & passes_scoring) with an endpoint in a multi-copy catalog "
                       "family -- i.e. the gene's copy IS grouped into a recovered family but is "
                       "max-overlap-mislabelled after a segdup neighbour. Never touches gene_of / block "
                       "projection / FP flags / grouping."),
            credited_genes=sorted(relabel_credit_mc),
            newly_credited=orc_cur["relabel_credited_genes"],
            n_newly_credited=len(orc_cur["relabel_credited_genes"]),
            detail={g: relabel_detail[g] for g in sorted(relabel_detail) if g in mc_oracle},
            recall_without_relabel=(round(
                (orc_cur["oracle_genes_recovered_multicopy"] - len(orc_cur["relabel_credited_genes"]))
                / n_mc_oracle, 4) if n_mc_oracle else None),
            recall_with_relabel=orc_cur["recall_oracle"]),
        residual_fp_roster=dict(
            multifam=fp_multifam, oversize=fp_oversize, allele=fp_allele,
            ep_impure_not_dna_named=impure_not_dna,
            ep_impure_total=len(impure_roster)),
        determinism=det,
    )
    json.dump(out, open(OUT_JSON, "w"), indent=1, sort_keys=True)
    P(f"\n[write] {OUT_TSV}")
    P(f"[write] {OUT_JSON}")
    return out


def _distinct_fp_blocks(res):
    """count DISTINCT offending blocks (a block may be flagged both oversize & multifam)."""
    seen = set()
    for e in res["multifam"] + res["oversize"] + res["allele"]:
        seen.add(frozenset(e[0]))
    return len(seen)


def _block_oracle_flags(b, ctx):
    """Per-BLOCK diploid-oracle FP flags, mirroring oracle_residuals EXACTLY (keyed by
    the block itself, not by its oracle-gene set -- two blocks can share an oracle gene)."""
    gene_of_dn = ctx["gene_of_dn"]; g2rows = ctx["g2rows"]; genes = ctx["genes"]
    gs = {gene_of_dn[dn] for dn in b if dn in gene_of_dn}
    og = sorted(g for g in gs if g in g2rows)
    if not og:
        return dict(oracle=[], dip="", multifam=0, oversize=0, allele=0)
    rows = [r for g in og for r in g2rows[g]]
    classes = {r["cls"] for r in rows}
    dip = max(r["dip"] for r in rows)
    dl = G.distinct_loci(b, genes)
    multifam = int(len(og) >= 2)
    allele = int(dl >= 2 and "single_locus_allele" in classes and "multi_copy" not in classes)
    oversize = 0
    if "single_locus_allele" not in classes and dip > 0 and dl > SW.DIPLOID_MARGIN * dip:
        oversize = 1
    return dict(oracle=og, dip=dip, multifam=multifam, oversize=oversize, allele=allele)


def write_family_tsv(cat_blocks, fam_ids, fam_dom, ctx, res, impure_roster):
    gene_of_dn = ctx["gene_of_dn"]; g2f = ctx["g2f"]; mega = ctx["mega"]; genes = ctx["genes"]
    impure_idx = {b["block_index"] for b in impure_roster}
    with open(OUT_TSV, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family_id", "dominant_gene", "n_loci_distinct", "n_genes", "genes",
                    "ep_impure", "n_protein_families", "oracle_genes", "diploid_CN",
                    "fp_multifam", "fp_oversize", "fp_allele"])
        for bi, (fid, b) in enumerate(zip(fam_ids, cat_blocks)):
            gs = sorted({gene_of_dn[dn] for dn in b if dn in gene_of_dn})
            fam_of = defaultdict(list)
            for g in gs:
                if g in g2f and g2f[g] not in mega:
                    fam_of[g2f[g]].append(g)
            fl = _block_oracle_flags(b, ctx)
            w.writerow([fid, fam_dom.get(fid, ""), G.distinct_loci(b, genes), len(gs),
                        ",".join(gs), int(bi in impure_idx), len(fam_of),
                        ",".join(fl["oracle"]), fl["dip"],
                        fl["multifam"], fl["oversize"], fl["allele"]])


if __name__ == "__main__":
    main()
