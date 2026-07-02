#!/usr/bin/env python
"""RNA-ONLY family refinement -- PRODUCTION opt-in stage.

WHAT THIS IS
------------
The shipped RNA multi-copy family definition (bench/denovo_families.py) confirms a
family edge iff `core_recip >= 0.13` (edge-CREATION threshold), then families = the
gamma-quasi-clique refinement (genome_family_def.refine_families, gamma=0.20, seed=0)
of the connected components.  bench/rna_only_edge_oracle.py + bench/RNA_ONLY_EDGE_ORACLE.md
LEARNED + VALIDATED an RNA-only refinement gate on top of that.  THIS FILE productionizes
the RECALL-PRESERVING operating point of that oracle as a clean, opt-in, default-OFF stage.

It does NOT touch bench/denovo_families.py, so the shipped byte-identical path is
untouched by construction.  It re-uses the oracle's exact loaders/thresholds and the
shipped refiner -- nothing here is re-derived.

THE RULE (recall-preserving deploy point, RNA_ONLY_EDGE_ORACLE.md sec 2)
-----------------------------------------------------------------------
  1. KEEP a family edge iff  core_recip >= 0.19  AND  aln_frac >= 0.24   (else CUT).
     - core_recip : max whole-transcript reciprocal homology weight over the DN edges
       of the gene pair (bench/denovo_family_edges.tsv).  Absent => 0.0 (transitive-
       closure / non-arbitration pair -> CUT).  Matches rna_only_edge_oracle.decide_recall.
     - aln_frac   : leakage-free UNIVERSAL longest shared spliced-exon-body fraction
       (bench/ri_sharedlen_universal.tsv).  Absent => 0.0 -> CUT.  The universal cache
       column `in_ep` (protein label) is NEVER read (leakage-free by construction).
     - within-gene / unannotated DN edges (ga is None / gb is None / ga == gb) are ALWAYS
       kept -- they are never a cross-gene over-merge (matches oracle build_kept).
  2. gamma-quasi-clique refinement: genome_family_def.refine_families(gamma=0.20, seed=0)
     (unchanged shipped operator; includes the >=2-distinct-loci multi-copy predicate).
  3. ALLELE-DEMOTE: a same-gene multi-locus family whose dominant gene is a balanced
     diploid het (balanced_frac >= 0.90 AND copy_like <= 0.10, read-consensus O1 signal
     bench/a1_read_consensus_o1.tsv) is ALLELIC, not multi-copy -> split to singletons
     (dropped from the catalog).  Exact thresholds/logic reused from
     rna_only_edge_oracle.apply_demote / demote_gene.

RNA-ONLY GUARD (thesis-critical)
--------------------------------
The INFERENCE feature set is exactly {core_recip, aln_frac} (edge decision) +
{balanced_frac, copy_like} (demote) and is hard-asserted DISJOINT from any
DNA/protein/genome column (in_dna_loose, in_ep, ep_tier, sedef, asm_hapCN, bridge_mask,
...).  DNA/protein/genome enter ONLY the VALIDATION report (residual-FP scoring), never
a decision.

OPT-IN, DEFAULT-OFF
-------------------
Runs only with env RUSTLE_RNA_ORACLE=1 OR flag --rna-oracle.  Otherwise prints one line
and exits 0 WITHOUT writing anything (the shipped path stays untouched).

DETERMINISM
-----------
PYTHONHASHSEED=0 (re-exec), fixed gamma=0.20 seed=0, sorted writes.  Re-runs are
byte-identical (see bench/test_family_rna_refine.py).

Writes: bench/family_rna_refine.tsv (family_id -> member loci/genes) + bench/family_rna_refine.json
Run:    RUSTLE_RNA_ORACLE=1 /home/juanfra/miniforge3/bin/python bench/family_rna_refine.py
   or:  /home/juanfra/miniforge3/bin/python bench/family_rna_refine.py --rna-oracle
"""
import os
import sys

# --- determinism: pin the hash seed BEFORE anything imports (re-exec preserves argv) ---
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import argparse
import json
from collections import Counter

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)

# re-used modules: shipped loaders + shipped refiner + the oracle's EXACT feature loaders,
# demote thresholds, residual roster and validation eval.  Nothing re-derived here.
import family_er_pr as FP
import genome_family_def as G
import graph_def_refine_sweep as SW
import rna_only_edge_oracle as RO

# --------------------------------------------------------------------------- constants
# RECALL-PRESERVING deploy point (RNA_ONLY_EDGE_ORACLE.md sec 2; do NOT re-derive):
CORE_MIN = 0.19          # core_recip threshold
ALN_MIN = 0.24           # aln_frac  threshold
GAMMA = G.GAMMA          # 0.20 (shipped gamma-quasi-clique cohesion)
SEED = SW.SEED           # 0    (shipped splitter witness seed)
# allele DEMOTE thresholds (reused from rna_only_edge_oracle.demote_gene):
DEMOTE_BAL_MIN = 0.90    # balanced_frac >= 0.90  (~0.5 minor-allele = diploid het)
DEMOTE_COPY_MAX = 0.10   # copy_like    <= 0.10  (not ~1/K = a real copy)

# RNA-only inference feature contract (hard-asserted):
EDGE_DECISION_FEATURES = ("core_recip", "aln_frac")
DEMOTE_FEATURES = ("balanced_frac", "copy_like")
DNA_FORBIDDEN = {
    "in_dna_loose", "in_dna", "in_ep", "ep_tier", "class", "cls", "cls_auth",
    "sedef", "sedef_identity", "sedef_corr", "asm_hapCN", "hap_CN_mat", "hap_CN_pat",
    "dip", "hap", "bridge_mask", "abl_bridge_mask", "mask_a", "mask_b",
}

OUT_TSV = os.path.join(BENCH, "family_rna_refine.tsv")
OUT_JSON = os.path.join(BENCH, "family_rna_refine.json")

assert abs(GAMMA - 0.20) < 1e-9 and SEED == 0, "gamma/seed drifted from the shipped constants"


# --------------------------------------------------------------------------- guards
def rna_only_guard():
    """Hard-assert the inference feature set is exactly the RNA contract and disjoint
    from every DNA/protein/genome column.  Fails LOUD if any label leaks into a decision."""
    infer = set(EDGE_DECISION_FEATURES) | set(DEMOTE_FEATURES)
    assert infer == {"core_recip", "aln_frac", "balanced_frac", "copy_like"}, \
        f"inference feature set drifted: {sorted(infer)}"
    leak = infer & DNA_FORBIDDEN
    assert not leak, f"DNA/protein/genome column in the inference path: {sorted(leak)}"


# --------------------------------------------------------------------------- build (RNA-only)
def build_catalog():
    """Apply the recall-preserving RNA-only gate + shipped gamma refinement + allele demote.
    Returns dict with the multi-copy catalog and the RNA-only bookkeeping.  No DNA read here."""
    rna_only_guard()

    # ---- RNA features (exact oracle loaders) ----
    gene_of_dn = RO.load_gene_of_dn()            # DN locus -> gene symbol (floored projection)
    pair_core = RO.load_pair_core(gene_of_dn)    # gene-pair -> max core_recip (denovo_family_edges.tsv)
    univ_aln = RO.load_universal_aln()           # gene-pair -> aln_frac (ri_sharedlen_universal.tsv; in_ep IGNORED)
    allele = RO.load_allele()                    # gene -> balanced_frac/copy_like/... (a1_read_consensus_o1.tsv)

    # ---- shipped graph context ----
    meta = FP.load_meta(); annot = FP.load_annot(); gene_of = FP.gene_of_factory(annot)
    raw_fams = FP.load_raw_families(); edge_pairs = FP.load_edges()
    genes, gene_of_dn2, *_ = FP.build_genes_dict(raw_fams, meta, gene_of)
    all_nodes = set()
    for f in raw_fams:
        all_nodes.update(f)

    # ---- RNA-only KEEP/CUT decision on cross-gene pairs ----
    def decide(k):
        c = pair_core.get(k)
        c = c if c is not None else 0.0
        a = univ_aln.get(k)
        a = a if a is not None else 0.0
        return (c >= CORE_MIN) and (a >= ALN_MIN)

    import networkx as nx
    Gr = nx.Graph(); Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b)

    kept = set()
    kept_pairs, cut_pairs = set(), set()
    n_dn_within, n_dn_cross_kept, n_dn_cross_cut = 0, 0, 0
    for u, v in Gr.edges():
        ga, gb = gene_of_dn2.get(u), gene_of_dn2.get(v)
        if ga is None or gb is None or ga == gb:
            kept.add(frozenset((u, v)))          # within-gene / unannotated: never an over-merge
            n_dn_within += 1
            continue
        k = frozenset((ga, gb))
        if decide(k):
            kept.add(frozenset((u, v))); n_dn_cross_kept += 1; kept_pairs.add(k)
        else:
            n_dn_cross_cut += 1; cut_pairs.add(k)

    # ---- shipped gamma-quasi-clique refinement (unchanged operator, gamma=0.20 seed=0) ----
    comps = SW.components_from_edges(all_nodes, kept)
    refined = G.refine_families(comps, [tuple(e) for e in kept], genes, GAMMA, SEED)

    # ---- allele DEMOTE (RNA read signal only; exact oracle logic) ----
    def demote_gene(g):
        a = allele.get(g)
        return a is not None and a["balanced_frac"] >= DEMOTE_BAL_MIN and a["copy_like"] <= DEMOTE_COPY_MAX

    catalog, demotions = [], []
    for b in refined:
        gs = [gene_of_dn2[dn] for dn in b if dn in gene_of_dn2]
        if gs:
            dom, cnt = Counter(gs).most_common(1)[0]
            homog = cnt / len(gs)
            if demote_gene(dom) and homog >= 0.5 and G.distinct_loci(b, genes) >= 2:
                demotions.append(dict(
                    gene=dom, n_loci=G.distinct_loci(b, genes),
                    balanced_frac=allele[dom]["balanced_frac"],
                    copy_like=allele[dom]["copy_like"],
                    dna_confirmed=(dom in RO.DNA_RESIDUAL_FP["allele_as_copy"])))
                continue                          # alleles -> singletons, dropped from the catalog
        catalog.append(sorted(b))

    return dict(
        catalog=catalog, demotions=demotions,
        gene_of_dn=gene_of_dn2, genes=genes, raw_fams=raw_fams, edge_pairs=edge_pairs,
        n_dn_edges_total=Gr.number_of_edges(),
        n_dn_within=n_dn_within, n_dn_cross_kept=n_dn_cross_kept, n_dn_cross_cut=n_dn_cross_cut,
        n_cross_pairs_kept=len(kept_pairs), n_cross_pairs_cut=len(cut_pairs - kept_pairs),
    )


# --------------------------------------------------------------------------- validate (DNA = scoring only)
def validate(built):
    """VALIDATION ONLY (never gates an edge): score the residual DNA-confirmed FP roster
    from GRAPH_DEF_REFINE_SWEEP with the oracle's exact assembly-independent residual eval,
    and count how many of the 12 named FP the RNA-only definition removes vs shipped gamma."""
    catalog = built["catalog"]; gene_of_dn = built["gene_of_dn"]; genes = built["genes"]
    g2rows = SW.load_oracle()                    # diploid CN oracle (DNA -- scoring only)

    res = RO.oracle_residuals(catalog, gene_of_dn, genes, g2rows, G, margin=SW.DIPLOID_MARGIN)
    # shipped gamma baseline catalog (same refiner, no RNA gate) for the removed-count
    shipped = G.refine_families(built["raw_fams"], built["edge_pairs"], genes, GAMMA, SEED)
    res_ship = RO.oracle_residuals(shipped, gene_of_dn, genes, g2rows, G, margin=SW.DIPLOID_MARGIN)

    def genes_in(examples):
        return {g for e in examples for g in e[0]}
    allele_mine, allele_ship = genes_in(res["allele"]), genes_in(res_ship["allele"])
    over_mine, over_ship = genes_in(res["oversize"]), genes_in(res_ship["oversize"])

    # named-FP transitions (present in shipped -> removed in RNA-only), mirrors the oracle's tracking:
    rm_allele = sum(1 for g in RO.DNA_RESIDUAL_FP["allele_as_copy"]
                    if g in allele_ship and g not in allele_mine)
    rm_oversize = sum(1 for g in RO.DNA_RESIDUAL_FP["oversize_diploid"]
                      if g in over_ship and g not in over_mine)
    # multifam: GSTM2 hub instance-count delta + explicit spanned-oracle-gene pairs
    def hub_count(examples, hub):
        return sum(1 for e in examples if hub in e[0])
    def comembered(examples, a, b):
        return any({a, b} <= set(e[0]) for e in examples)
    rm_multifam = max(0, hub_count(res_ship["multifam"], "GSTM2") - hub_count(res["multifam"], "GSTM2"))
    for a, b in [("FOXO1", "LOC115933254"), ("LOC101142904", "LOC129526550")]:
        if comembered(res_ship["multifam"], a, b) and not comembered(res["multifam"], a, b):
            rm_multifam += 1
    named_removed = rm_allele + rm_oversize + rm_multifam

    remaining = dict(allele=len(res["allele"]), oversize=len(res["oversize"]),
                     multifam=len(res["multifam"]))
    shipped_counts = dict(allele=len(res_ship["allele"]), oversize=len(res_ship["oversize"]),
                          multifam=len(res_ship["multifam"]))
    return dict(
        residual_remaining=remaining,
        residual_remaining_total=sum(remaining.values()),
        shipped_residual=shipped_counts,
        shipped_residual_total=sum(shipped_counts.values()),
        named_removed=named_removed,
        named_removed_breakdown=dict(allele=rm_allele, oversize=rm_oversize, multifam=rm_multifam),
        oracle_genes_recovered=res["n_recovered"],
        oracle_genes_recovered_shipped=res_ship["n_recovered"],
        residual_examples=dict(
            allele=RO._fmt_examples(res["allele"], "allele"),
            oversize=RO._fmt_examples(res["oversize"], "oversize"),
            multifam=RO._fmt_examples(res["multifam"], "multifam")),
    )


# --------------------------------------------------------------------------- write
def write_outputs(built, val):
    catalog = built["catalog"]; gene_of_dn = built["gene_of_dn"]; genes = built["genes"]
    # deterministic family_id: sort families by their sorted member tuple, then long-format rows
    fams_sorted = sorted(catalog, key=lambda b: tuple(sorted(b)))
    with open(OUT_TSV, "w") as out:
        out.write("family_id\tn_loci\tdominant_gene\tmember_dn\tmember_gene\tchrom\tstart\tend\n")
        for fid, b in enumerate(fams_sorted):
            gs = [gene_of_dn.get(dn, "NA") for dn in b]
            dom = Counter(g for g in gs if g != "NA").most_common(1)
            dom = dom[0][0] if dom else "NA"
            nl = G.distinct_loci(b, genes)
            for dn in sorted(b):
                g = genes[dn]
                out.write(f"{fid}\t{nl}\t{dom}\t{dn}\t{gene_of_dn.get(dn, 'NA')}\t"
                          f"{g['chrom']}\t{g['start']}\t{g['end']}\n")

    summary = dict(
        stage="family_rna_refine (RNA-only recall-preserving refinement; opt-in)",
        rule=dict(edge="KEEP iff core_recip>=%.2f AND aln_frac>=%.2f" % (CORE_MIN, ALN_MIN),
                  core_recip_min=CORE_MIN, aln_frac_min=ALN_MIN,
                  gamma=GAMMA, seed=SEED,
                  demote="balanced_frac>=%.2f AND copy_like<=%.2f" % (DEMOTE_BAL_MIN, DEMOTE_COPY_MAX),
                  demote_balanced_frac_min=DEMOTE_BAL_MIN, demote_copy_like_max=DEMOTE_COPY_MAX),
        n_families=len(catalog),
        edges=dict(
            n_dn_edges_total=built["n_dn_edges_total"],
            n_dn_within_gene_kept=built["n_dn_within"],
            n_dn_cross_gene_kept=built["n_dn_cross_kept"],
            n_dn_cross_gene_cut=built["n_dn_cross_cut"],
            n_cross_gene_pairs_kept=built["n_cross_pairs_kept"],
            n_cross_gene_pairs_cut=built["n_cross_pairs_cut"]),
        n_alleles_demoted=len(built["demotions"]),
        alleles_demoted=sorted(built["demotions"], key=lambda d: d["gene"]),
        residual_fp=dict(
            note=("12 DNA-confirmed residual FP in shipped gamma "
                  "(2 allele + 4 oversize + 6 multifam = %d); "
                  "RNA-only recall-preserving+demote removes %d named FP "
                  "(RNA_ONLY_EDGE_ORACLE.md recall-preserving row: oracle_allele %d, "
                  "oracle_oversize %d, oracle_multifam %d remaining)"
                  % (val["shipped_residual_total"], val["named_removed"],
                     val["residual_remaining"]["allele"], val["residual_remaining"]["oversize"],
                     val["residual_remaining"]["multifam"])),
            **val),
        guards=dict(
            edge_decision_features=list(EDGE_DECISION_FEATURES),
            demote_features=list(DEMOTE_FEATURES),
            no_dna_in_inference=True,
            gamma=GAMMA, seed=SEED),
        inputs=dict(
            edges="bench/denovo_family_edges.tsv",
            aln_frac="bench/ri_sharedlen_universal.tsv",
            allele="bench/a1_read_consensus_o1.tsv"),
        outputs=dict(catalog_tsv="bench/family_rna_refine.tsv",
                     summary_json="bench/family_rna_refine.json"),
    )
    with open(OUT_JSON, "w") as out:
        json.dump(summary, out, sort_keys=True, indent=1,
                  default=lambda x: None if (isinstance(x, float) and x != x) else x)
    return summary


# --------------------------------------------------------------------------- driver
def run(write=True):
    built = build_catalog()
    val = validate(built)
    summary = write_outputs(built, val) if write else None
    return built, val, summary


def _report(built, val, summary):
    P = print
    P("\n==================== RNA-ONLY FAMILY REFINEMENT (opt-in) ====================")
    P(f"rule : KEEP iff core_recip>={CORE_MIN:.2f} AND aln_frac>={ALN_MIN:.2f}  ->  "
      f"gamma-refine (gamma={GAMMA}, seed={SEED})  ->  allele-demote "
      f"(balanced_frac>={DEMOTE_BAL_MIN:.2f} AND copy_like<={DEMOTE_COPY_MAX:.2f})")
    P(f"DN edges         : total={built['n_dn_edges_total']}  within-gene kept={built['n_dn_within']}  "
      f"cross-gene kept={built['n_dn_cross_kept']}  cross-gene CUT={built['n_dn_cross_cut']}")
    P(f"cross-gene pairs : kept={built['n_cross_pairs_kept']}  cut={built['n_cross_pairs_cut']}")
    P(f"n_families       : {len(built['catalog'])}")
    P(f"alleles demoted  : {len(built['demotions'])}  "
      + ", ".join(f"{d['gene']}(dl={d['n_loci']},bal={d['balanced_frac']:.2f},"
                  f"copy_like={d['copy_like']:.2f}{',DNA-confirmed' if d['dna_confirmed'] else ',novel'})"
                  for d in sorted(built['demotions'], key=lambda d: d['gene'])))
    P(f"residual FP      : shipped total={val['shipped_residual_total']} "
      f"(allele {val['shipped_residual']['allele']}/oversize {val['shipped_residual']['oversize']}/"
      f"multifam {val['shipped_residual']['multifam']})  ->  remaining={val['residual_remaining_total']} "
      f"(allele {val['residual_remaining']['allele']}/oversize {val['residual_remaining']['oversize']}/"
      f"multifam {val['residual_remaining']['multifam']})")
    P(f"named FP removed : {val['named_removed']}/12  "
      f"(allele {val['named_removed_breakdown']['allele']}, "
      f"oversize {val['named_removed_breakdown']['oversize']}, "
      f"multifam {val['named_removed_breakdown']['multifam']})")
    P(f"oracle recovery  : shipped {val['oracle_genes_recovered_shipped']} -> "
      f"RNA-only {val['oracle_genes_recovered']}")
    if summary is not None:
        P(f"wrote {OUT_TSV}\nwrote {OUT_JSON}")
    P("============================================================================")


def main(argv=None):
    ap = argparse.ArgumentParser(description="RNA-only family refinement (opt-in, default-OFF).")
    ap.add_argument("--rna-oracle", action="store_true",
                    help="enable the opt-in RNA-only refinement stage")
    ap.add_argument("--no-write", action="store_true",
                    help="run + report but do NOT write outputs (used by the self-check)")
    args = ap.parse_args(argv)
    enabled = args.rna_oracle or os.environ.get("RUSTLE_RNA_ORACLE") == "1"
    if not enabled:
        print("opt-in: set RUSTLE_RNA_ORACLE=1 or --rna-oracle")
        return 0
    built, val, summary = run(write=not args.no_write)
    _report(built, val, summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
