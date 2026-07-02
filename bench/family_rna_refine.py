#!/usr/bin/env python
"""RNA-ONLY family refinement -- PRODUCTION opt-in stage.

WHAT THIS IS
------------
The shipped RNA multi-copy family definition (bench/denovo_families.py) confirms a
family edge iff `core_recip >= 0.13` (edge-CREATION threshold), then families = the
gamma-quasi-clique refinement (genome_family_def.refine_families, gamma=0.20, seed=0)
of the connected components.  bench/rna_only_edge_oracle.py + bench/RNA_ONLY_EDGE_ORACLE.md
LEARNED + VALIDATED an RNA-only refinement gate on top of that.  THIS FILE productionizes
the RECALL-PRESERVING operating point of that oracle as the DEFAULT family definition
(default-ON; opt out with --legacy).

It does NOT touch bench/denovo_families.py, so the LEGACY core_recip>=0.13 path stays
bit-for-bit reproducible (--legacy / RUSTLE_RNA_ORACLE=0).  It re-uses the oracle's exact
loaders/thresholds and the shipped refiner -- nothing here is re-derived.

THE RULE (recall-preserving deploy point, RNA_ONLY_EDGE_ORACLE.md sec 2)
-----------------------------------------------------------------------
  1. KEEP a family edge iff  core_recip >= 0.19  AND  aln_frac >= 0.24
     AND NOT( min_shared_mult >= 20 )   (else CUT).
     - core_recip : max whole-transcript reciprocal homology weight over the DN edges
       of the gene pair (bench/denovo_family_edges.tsv).  Absent => 0.0 (transitive-
       closure / non-arbitration pair -> CUT).  Matches rna_only_edge_oracle.decide_recall.
     - aln_frac   : leakage-free UNIVERSAL longest shared spliced-exon-body fraction
       (bench/ri_sharedlen_universal.tsv).  Absent => 0.0 -> CUT.  The universal cache
       column `in_ep` (protein label) is NEVER read (leakage-free by construction).
     - min_shared_mult : REPEAT-HUB GATE (DEFAULT-ON; --no-repeat-gate / RUSTLE_NO_REPEAT_GATE=1
       ablates it).  The LOWEST canonical-minimizer multiplicity among the VG nodes the two
       genes SHARE (bench/vg_repeat_catalog.py / .tsv per-edge rows; multiplicity = # distinct
       genes traversing a node, LIBRARY-FREE = read/structure-derived, NOT soft-mask).  A cross-
       gene edge whose shared sequence is ONLY an EXTREME repeat (min_shared_mult >= REPEAT_MULT_MIN)
       is CUT even if core_recip/aln_frac pass -- this is the ONE residual over-merge class that
       alignment cannot cut (extreme repeat-bridge hubs, e.g. fam17: 27 genes / 16 protein families
       joined by a shared Alu/poly-A bridge node, mult up to 503).  REPEAT_MULT_MIN = 20 is chosen
       at the extreme tail (VG_REPEAT_CATALOG.md sec 4: mult>=20 -> 92.7% RepeatMasker-concordant)
       and cleanly SEPARATES the true fam17 hub (per-edge min_shared_mult up to 38, 43/301 edges
       >=20) from the negative controls GSTM2 (protein domain, per-edge max 9) and MAGE (cardinality,
       per-edge max 8) -- both have ZERO edges >= 20, so the gate NEVER touches them.  Absent => no
       repeat cut (falls through to core+aln).  min_shared_mult is loaded from the VG catalog output;
       the minimizers are NOT re-derived here.
     - within-gene / unannotated DN edges (ga is None / gb is None / ga == gb) are ALWAYS
       kept -- they are never a cross-gene over-merge (matches oracle build_kept), and are NEVER
       subject to the repeat-hub gate.
  2. gamma-quasi-clique refinement: genome_family_def.refine_families(gamma=0.20, seed=0)
     (unchanged shipped operator; includes the >=2-distinct-loci multi-copy predicate).
  3. ALLELE-DEMOTE: a same-gene multi-locus family whose dominant gene is a balanced
     diploid het (balanced_frac >= 0.90 AND copy_like <= 0.10, read-consensus O1 signal
     bench/a1_read_consensus_o1.tsv) is ALLELIC, not multi-copy -> split to singletons
     (dropped from the catalog).  Exact thresholds/logic reused from
     rna_only_edge_oracle.apply_demote / demote_gene.

RNA-ONLY / LIBRARY-FREE GUARD (thesis-critical)
-----------------------------------------------
The INFERENCE feature set is exactly {core_recip, aln_frac} (alignment edge decision) +
{min_shared_mult, cyclic} (repeat-hub gate) + {balanced_frac, copy_like} (demote), and is
hard-asserted DISJOINT from any DNA/protein/genome column (in_dna_loose, in_ep, ep_tier,
sedef, asm_hapCN, bridge_mask, ...) AND from any soft-mask / RepeatMasker / RepBase / Dfam
column.  The repeat-hub feature {min_shared_mult, cyclic} is VG canonical-minimizer
MULTIPLICITY (# distinct genes traversing a node) -- read/structure-derived and LIBRARY-FREE;
soft-mask is used NOWHERE in the gate (it is only VALIDATION in vg_repeat_catalog.py sec 4).
DNA/protein/genome/soft-mask enter ONLY the VALIDATION report, never a decision.

DEFAULT-ON (opt out with --legacy; ablate the repeat gate with --no-repeat-gate)
--------------------------------------------------------------------------------
The RNA-only refinement (core+aln + repeat-hub gate + gamma + demote) is now the DEFAULT
family definition -- runs by default.  --no-repeat-gate (or env RUSTLE_NO_REPEAT_GATE=1)
ABLATES just the repeat-hub gate (keeps core+aln+gamma+demote) and recovers the
pre-repeat-gate catalog.  The legacy core_recip>=0.13 catalog is recovered with --legacy OR
env RUSTLE_RNA_ORACLE=0 (prints one line, exits 0 without writing; run bench/denovo_families.py
for the legacy catalog).  Because this stage never edits denovo_families.py, the legacy path
remains bit-for-bit reproducible.

HIGH-PRECISION (opt in with --high-precision / RUSTLE_HIGH_PRECISION=1)
----------------------------------------------------------------------
--high-precision swaps ONLY the gamma-quasi-clique cohesion from the recall-preserving default
GAMMA=0.20 to HIGH_PRECISION_GAMMA=0.40 (bench/PRECISION_RECALL_FRONTIER.md recommended point);
core/aln thresholds, the repeat-hub gate and the allele demote are UNCHANGED.  It removes the two
collapsed-array OVERSIZE blobs (MPHOSPH8, LOC134758618) -> distinct FP 6->4, P_fixed48 0.917,
recall held 48/57 (nFam 606 -> 623).  HONEST costs, carried in the summary JSON + report (not
dropped): off-oracle KRAB-ZNF over-split (gamma>=0.27) and the MAGE X-array DNA-only floor that
survives at every RNA point.  Default (no flag) stays byte-identical; --high-precision writes the
gamma=0.40 catalog and records the active gamma + caveats in the summary.

DETERMINISM
-----------
PYTHONHASHSEED=0 (re-exec), fixed gamma=0.20 seed=0, sorted writes.  Re-runs are
byte-identical (see bench/test_family_rna_refine.py).

Writes: bench/family_rna_refine.tsv (family_id -> member loci/genes) + bench/family_rna_refine.json
Run:    /home/juanfra/miniforge3/bin/python bench/family_rna_refine.py                   (default: refined catalog)
hi-prec:/home/juanfra/miniforge3/bin/python bench/family_rna_refine.py --high-precision  (gamma=0.40 catalog)
legacy: /home/juanfra/miniforge3/bin/python bench/family_rna_refine.py --legacy          (opt out -> nothing written)
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
# REPEAT-HUB GATE threshold (VG_REPEAT_CATALOG.md; do NOT re-derive -- picked from the data):
#   min_shared_mult >= REPEAT_MULT_MIN => the pair's shared sequence is ONLY an extreme
#   repeat node (# distinct genes traversing >= 20).  REPEAT_MULT_MIN = 20 is the extreme
#   tail where RepeatMasker concordance = 92.7% (VG_REPEAT_CATALOG.md sec 4) AND the value
#   that SEPARATES the fam17 hub (per-edge min_shared_mult up to 38; 43/301 edges >= 20)
#   from the negative controls GSTM2 (per-edge max 9) and MAGE (per-edge max 8), both of
#   which have ZERO edges >= 20 => the gate cannot touch them.
REPEAT_MULT_MIN = 20     # VG min_shared_mult cut (library-free minimizer multiplicity)
# allele DEMOTE thresholds (reused from rna_only_edge_oracle.demote_gene):
DEMOTE_BAL_MIN = 0.90    # balanced_frac >= 0.90  (~0.5 minor-allele = diploid het)
DEMOTE_COPY_MAX = 0.10   # copy_like    <= 0.10  (not ~1/K = a real copy)

# HIGH-PRECISION operating point (bench/PRECISION_RECALL_FRONTIER.md recommended point; do NOT
# re-derive).  --high-precision (or env RUSTLE_HIGH_PRECISION=1) swaps ONLY the gamma-quasi-clique
# cohesion from the recall-preserving default GAMMA=0.20 to HIGH_PRECISION_GAMMA=0.40.  Everything
# else (core/aln thresholds, repeat-hub gate, allele demote) is UNCHANGED.  From the frontier
# (gamma=0.40 row):
#   - removes the two collapsed-array OVERSIZE blobs MPHOSPH8 + LOC134758618 (the fam17 repeat-hubs
#     are ALREADY removed by the default-on repeat gate in BOTH modes): distinct over-merge FP
#     blocks 6 -> 4, fixed-denominator precision P_fixed48 0.875 -> 0.917 (moving P_dedup 0.920);
#   - recall HELD at 48/57 recovered, ZERO on-oracle genes lost (nFam 606 -> 623).
# HONEST costs (carried in the summary JSON + report; do NOT drop them):
#   - OFF-ORACLE KRAB-ZNF cost: gamma>=0.27 trips the ZNF716/KRAB-ZNF knife-edge (family density
#     0.261) and over-splits divergent KRAB-ZNF paralog families the sparse high-CN diploid oracle
#     cannot see (NOT captured by the recovered-somewhere recall metric); default gamma=0.20 keeps them;
#   - over-split surrogates: undersize 33 -> 37 (+4), 15 divergent-paralog (TRUTHBAR) pairs cut;
#   - MAGE floor: the dense-uniform X-array LOC129529978+LOC129529986 (Q=-0.001) SURVIVES at every
#     RNA operating point -- the DNA-only cardinality floor; no gamma removes it.
# best-2-cut is NOT wired (DOMINATED: net-zero distinct-FP, cosmetic GSTM2 relabel at over-split cost).
HIGH_PRECISION_GAMMA = 0.40     # frontier high-precision gamma-quasi-clique cohesion
HIGH_PRECISION_NOTE = dict(
    source="bench/PRECISION_RECALL_FRONTIER.md (recommended high-precision operating point)",
    default_gamma=GAMMA,
    high_precision_gamma=HIGH_PRECISION_GAMMA,
    frontier_row_gamma040=("nFam~623, distinct over-merge FP blocks 4 (down from 6), "
                           "P_fixed48 0.917 (vs 0.875 default), moving P_dedup 0.920, "
                           "R 0.842 (48/57), zero on-oracle genes lost"),
    precision_impact=("gamma=0.40 removes the two collapsed-array OVERSIZE blobs MPHOSPH8 + "
                      "LOC134758618 (fam17 repeat-hubs already removed by the default-on repeat "
                      "gate in both modes) -> distinct FP 6->4 / P_fixed48 0.917"),
    offoracle_krabznf_cost=("HONEST off-oracle cost: gamma>=0.27 trips the ZNF716/KRAB-ZNF "
                            "knife-edge (family density 0.261) and over-splits divergent KRAB-ZNF "
                            "paralog families the sparse high-CN diploid oracle cannot see; NOT "
                            "captured by the recovered-somewhere recall metric. Default gamma=0.20 "
                            "preserves them."),
    over_split_cost="undersize 33->37 (+4); 15 divergent-paralog (TRUTHBAR) co-membership pairs cut",
    mage_floor=("MAGE-class dense-uniform X-array LOC129529978+LOC129529986 (Q=-0.001) SURVIVES at "
                "every RNA operating point -- the DNA-only cardinality floor; no gamma removes it"),
    best2cut="NOT wired -- DOMINATED (net-zero distinct-FP, cosmetic GSTM2 relabel at over-split cost)",
)

# RNA-only / library-free inference feature contract (hard-asserted):
EDGE_DECISION_FEATURES = ("core_recip", "aln_frac")            # alignment edge decision
REPEAT_GATE_FEATURES = ("min_shared_mult", "cyclic")          # VG minimizer multiplicity (library-free)
DEMOTE_FEATURES = ("balanced_frac", "copy_like")
DNA_FORBIDDEN = {
    "in_dna_loose", "in_dna", "in_ep", "ep_tier", "class", "cls", "cls_auth",
    "sedef", "sedef_identity", "sedef_corr", "asm_hapCN", "hap_CN_mat", "hap_CN_pat",
    "dip", "hap", "bridge_mask", "abl_bridge_mask", "mask_a", "mask_b",
}
# soft-mask / RepeatMasker / library columns are FORBIDDEN in the repeat-hub gate -- the gate
# is pure VG minimizer MULTIPLICITY (library-free); soft-mask is only external VALIDATION.
LIBRARY_FORBIDDEN = {
    "softmask", "soft_mask", "softmask_frac", "node_softmask", "mean_softmask",
    "mean_softmask_hi", "mean_softmask_all", "frac_softmasked", "repeatmasker",
    "repbase", "dfam", "rmsk", "te_class", "te_family", "mask",
}

# LIBRARY-FREE repeat-hub multiplicity source: VG catalog per-edge rows (NOT re-derived here).
VG_REPEAT_TSV = os.path.join(BENCH, "vg_repeat_catalog.tsv")

OUT_TSV = os.path.join(BENCH, "family_rna_refine.tsv")
OUT_JSON = os.path.join(BENCH, "family_rna_refine.json")

assert abs(GAMMA - 0.20) < 1e-9 and SEED == 0, "gamma/seed drifted from the shipped constants"
assert REPEAT_MULT_MIN == 20, "REPEAT_MULT_MIN drifted from the VG_REPEAT_CATALOG.md tail (20)"
assert abs(HIGH_PRECISION_GAMMA - 0.40) < 1e-9, "HIGH_PRECISION_GAMMA drifted from the frontier point (0.40)"


# --------------------------------------------------------------------------- guards
def rna_only_guard():
    """Hard-assert the inference feature set is exactly the RNA/library-free contract and disjoint
    from every DNA/protein/genome column AND every soft-mask/RepeatMasker/library column.  Fails
    LOUD if any external label leaks into a decision (edge, repeat-hub gate, or demote)."""
    infer = set(EDGE_DECISION_FEATURES) | set(DEMOTE_FEATURES)
    assert infer == {"core_recip", "aln_frac", "balanced_frac", "copy_like"}, \
        f"inference feature set drifted: {sorted(infer)}"
    leak = infer & DNA_FORBIDDEN
    assert not leak, f"DNA/protein/genome column in the inference path: {sorted(leak)}"
    # repeat-hub gate must be VG minimizer multiplicity ONLY -- library-free, no soft-mask/DNA.
    rep = set(REPEAT_GATE_FEATURES)
    assert rep == {"min_shared_mult", "cyclic"}, f"repeat-gate feature set drifted: {sorted(rep)}"
    leak_rep = rep & (DNA_FORBIDDEN | LIBRARY_FORBIDDEN)
    assert not leak_rep, f"soft-mask/RepeatMasker/DNA column in the repeat-hub gate: {sorted(leak_rep)}"


# --------------------------------------------------------------------------- repeat-hub multiplicity (library-free)
def load_repeat_mult():
    """Load the per-edge VG min_shared_mult keyed by gene-pair from bench/vg_repeat_catalog.tsv
    (the LIBRARY-FREE canonical-minimizer multiplicity catalog).  REUSES the VG catalog's exact
    per-edge computation -- the minimizers are NOT re-derived here.  Reads only the `min_shared_mult`
    column of the per-edge section; the soft-mask column is NEVER read.  Absent/blank => omitted
    (no repeat cut for that pair).  Deterministic (single ordered pass, dict keyed by frozenset)."""
    if not os.path.exists(VG_REPEAT_TSV):
        raise FileNotFoundError(
            f"repeat-hub gate is DEFAULT-ON but {VG_REPEAT_TSV} is missing; "
            f"run bench/vg_repeat_catalog.py, or ablate with --no-repeat-gate / RUSTLE_NO_REPEAT_GATE=1")
    out = {}
    with open(VG_REPEAT_TSV) as fh:
        in_edges, ix = False, None
        for ln in fh:
            if ln.startswith("# SECTION edges"):
                in_edges, ix = True, None
                continue
            if not in_edges:
                continue
            if ln.startswith("gene_a\t"):
                hdr = ln.rstrip("\n").split("\t")
                ix = {h: i for i, h in enumerate(hdr)}
                # LIBRARY-FREE guard: we consult ONLY min_shared_mult, never a soft-mask column.
                assert "min_shared_mult" in ix, "vg_repeat_catalog.tsv missing min_shared_mult column"
                continue
            if ix is None:
                continue
            f = ln.rstrip("\n").split("\t")
            msm = f[ix["min_shared_mult"]]
            if msm == "":
                continue
            out[frozenset((f[ix["gene_a"]], f[ix["gene_b"]]))] = int(msm)
    return out


# --------------------------------------------------------------------------- build (RNA-only)
def build_catalog(repeat_gate=True, gamma=GAMMA):
    """Apply the recall-preserving RNA-only gate + repeat-hub gate + shipped gamma refinement +
    allele demote.  Returns dict with the multi-copy catalog and the RNA-only bookkeeping.  No DNA
    read here.  repeat_gate=False ablates ONLY the repeat-hub gate (keeps core+aln+gamma+demote).
    gamma selects the gamma-quasi-clique cohesion: default GAMMA=0.20 (recall-preserving) or
    HIGH_PRECISION_GAMMA=0.40 (--high-precision; PRECISION_RECALL_FRONTIER.md).  Nothing else changes."""
    rna_only_guard()

    # ---- RNA features (exact oracle loaders) ----
    gene_of_dn = RO.load_gene_of_dn()            # DN locus -> gene symbol (floored projection)
    pair_core = RO.load_pair_core(gene_of_dn)    # gene-pair -> max core_recip (denovo_family_edges.tsv)
    univ_aln = RO.load_universal_aln()           # gene-pair -> aln_frac (ri_sharedlen_universal.tsv; in_ep IGNORED)
    allele = RO.load_allele()                    # gene -> balanced_frac/copy_like/... (a1_read_consensus_o1.tsv)
    # LIBRARY-FREE repeat-hub multiplicity (VG canonical-minimizer catalog; not re-derived):
    pair_repeat_mult = load_repeat_mult() if repeat_gate else {}

    # ---- shipped graph context ----
    meta = FP.load_meta(); annot = FP.load_annot(); gene_of = FP.gene_of_factory(annot)
    raw_fams = FP.load_raw_families(); edge_pairs = FP.load_edges()
    genes, gene_of_dn2, *_ = FP.build_genes_dict(raw_fams, meta, gene_of)
    all_nodes = set()
    for f in raw_fams:
        all_nodes.update(f)

    # ---- alignment KEEP/CUT decision on cross-gene pairs ----
    def core_aln_keep(k):
        c = pair_core.get(k)
        c = c if c is not None else 0.0
        a = univ_aln.get(k)
        a = a if a is not None else 0.0
        return (c >= CORE_MIN) and (a >= ALN_MIN)

    # ---- repeat-hub gate: shared sequence is ONLY an extreme repeat (library-free) ----
    def repeat_hub(k):
        m = pair_repeat_mult.get(k)              # absent => no repeat cut (fall through to core+aln)
        return (m is not None) and (m >= REPEAT_MULT_MIN)

    import networkx as nx
    Gr = nx.Graph(); Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b)

    kept = set()
    kept_pairs, cut_pairs, repeat_cut_pairs = set(), set(), set()
    n_dn_within, n_dn_cross_kept, n_dn_cross_cut, n_dn_cross_cut_repeat = 0, 0, 0, 0
    for u, v in Gr.edges():
        ga, gb = gene_of_dn2.get(u), gene_of_dn2.get(v)
        if ga is None or gb is None or ga == gb:
            kept.add(frozenset((u, v)))          # within-gene / unannotated: never an over-merge, never repeat-gated
            n_dn_within += 1
            continue
        k = frozenset((ga, gb))
        keep_ca = core_aln_keep(k)
        repeat_only = repeat_gate and repeat_hub(k)   # passes core+aln but shares ONLY extreme repeat
        if keep_ca and not repeat_only:
            kept.add(frozenset((u, v))); n_dn_cross_kept += 1; kept_pairs.add(k)
        else:
            n_dn_cross_cut += 1; cut_pairs.add(k)
            if keep_ca and repeat_only:               # cut SPECIFICALLY by the repeat-hub gate
                n_dn_cross_cut_repeat += 1; repeat_cut_pairs.add(k)

    # ---- shipped gamma-quasi-clique refinement (unchanged operator; gamma threaded, seed=0) ----
    comps = SW.components_from_edges(all_nodes, kept)
    refined = G.refine_families(comps, [tuple(e) for e in kept], genes, gamma, SEED)

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
        repeat_gate=repeat_gate, gamma=gamma,
        n_dn_edges_total=Gr.number_of_edges(),
        n_dn_within=n_dn_within, n_dn_cross_kept=n_dn_cross_kept, n_dn_cross_cut=n_dn_cross_cut,
        n_dn_cross_cut_repeat=n_dn_cross_cut_repeat,
        n_cross_pairs_kept=len(kept_pairs), n_cross_pairs_cut=len(cut_pairs - kept_pairs),
        n_cross_pairs_cut_repeat=len(repeat_cut_pairs - kept_pairs),
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
    gamma = built["gamma"]
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

    repeat_gate = built["repeat_gate"]
    summary = dict(
        stage="family_rna_refine (RNA-only recall-preserving refinement + repeat-hub gate; "
              "DEFAULT-ON, opt out --legacy / ablate gate --no-repeat-gate)",
        rule=dict(edge="KEEP iff core_recip>=%.2f AND aln_frac>=%.2f AND NOT(min_shared_mult>=%d)"
                       % (CORE_MIN, ALN_MIN, REPEAT_MULT_MIN),
                  core_recip_min=CORE_MIN, aln_frac_min=ALN_MIN,
                  repeat_gate_enabled=repeat_gate, repeat_mult_min=REPEAT_MULT_MIN,
                  gamma=gamma, seed=SEED,
                  demote="balanced_frac>=%.2f AND copy_like<=%.2f" % (DEMOTE_BAL_MIN, DEMOTE_COPY_MAX),
                  demote_balanced_frac_min=DEMOTE_BAL_MIN, demote_copy_like_max=DEMOTE_COPY_MAX),
        n_families=len(catalog),
        edges=dict(
            n_dn_edges_total=built["n_dn_edges_total"],
            n_dn_within_gene_kept=built["n_dn_within"],
            n_dn_cross_gene_kept=built["n_dn_cross_kept"],
            n_dn_cross_gene_cut=built["n_dn_cross_cut"],
            n_dn_cross_gene_cut_by_repeat_gate=built["n_dn_cross_cut_repeat"],
            n_cross_gene_pairs_kept=built["n_cross_pairs_kept"],
            n_cross_gene_pairs_cut=built["n_cross_pairs_cut"],
            n_cross_gene_pairs_cut_by_repeat_gate=built["n_cross_pairs_cut_repeat"]),
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
            repeat_gate_features=list(REPEAT_GATE_FEATURES),
            demote_features=list(DEMOTE_FEATURES),
            no_dna_in_inference=True,
            repeat_gate_library_free=True,
            no_softmask_in_repeat_gate=True,
            gamma=gamma, seed=SEED),
        inputs=dict(
            edges="bench/denovo_family_edges.tsv",
            aln_frac="bench/ri_sharedlen_universal.tsv",
            repeat_mult="bench/vg_repeat_catalog.tsv (min_shared_mult; library-free VG multiplicity)",
            allele="bench/a1_read_consensus_o1.tsv"),
        outputs=dict(catalog_tsv="bench/family_rna_refine.tsv",
                     summary_json="bench/family_rna_refine.json"),
    )
    # HIGH-PRECISION disclosure: only present when the flag/env selects gamma=0.40 (keeps the
    # default catalog byte-identical).  Carries the frontier's precision impact AND its HONEST
    # off-oracle KRAB-ZNF + MAGE-floor caveats -- do NOT drop these.
    if abs(gamma - GAMMA) > 1e-9:
        summary["high_precision"] = dict(
            active=True, active_gamma=gamma, n_families=len(catalog),
            live_precision_signal=(
                "residual oversize %d (default gamma=%.2f: 3 = MPHOSPH8 + LOC134758618 + MAGE "
                "X-array); gamma=%.2f removes the two collapsed-array OVERSIZE blobs "
                "MPHOSPH8 + LOC134758618, leaving only the MAGE X-array DNA-only floor"
                % (val["residual_remaining"]["oversize"], GAMMA, gamma)),
            oracle_genes_recovered=val["oracle_genes_recovered"],
            **HIGH_PRECISION_NOTE,
        )
    with open(OUT_JSON, "w") as out:
        json.dump(summary, out, sort_keys=True, indent=1,
                  default=lambda x: None if (isinstance(x, float) and x != x) else x)
    return summary


# --------------------------------------------------------------------------- driver
def run(write=True, repeat_gate=True, gamma=GAMMA):
    built = build_catalog(repeat_gate=repeat_gate, gamma=gamma)
    val = validate(built)
    summary = write_outputs(built, val) if write else None
    return built, val, summary


def _report(built, val, summary):
    P = print
    rg = built["repeat_gate"]
    gm = built["gamma"]
    hp = abs(gm - GAMMA) > 1e-9
    P("\n==================== RNA-ONLY FAMILY REFINEMENT (%s) ===================="
      % ("HIGH-PRECISION gamma=%.2f" % gm if hp else "default"))
    P(f"rule : KEEP iff core_recip>={CORE_MIN:.2f} AND aln_frac>={ALN_MIN:.2f} AND "
      f"NOT(min_shared_mult>={REPEAT_MULT_MIN}) [repeat-hub gate {'ON' if rg else 'OFF (ablated)'}]  ->  "
      f"gamma-refine (gamma={gm}{' [HIGH-PRECISION]' if hp else ''}, seed={SEED})  ->  allele-demote "
      f"(balanced_frac>={DEMOTE_BAL_MIN:.2f} AND copy_like<={DEMOTE_COPY_MAX:.2f})")
    P(f"DN edges         : total={built['n_dn_edges_total']}  within-gene kept={built['n_dn_within']}  "
      f"cross-gene kept={built['n_dn_cross_kept']}  cross-gene CUT={built['n_dn_cross_cut']}  "
      f"(of which by repeat-hub gate={built['n_dn_cross_cut_repeat']})")
    P(f"cross-gene pairs : kept={built['n_cross_pairs_kept']}  cut={built['n_cross_pairs_cut']}  "
      f"(repeat-hub-gate cut pairs={built['n_cross_pairs_cut_repeat']})")
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
    if hp:
        P(f"HIGH-PRECISION   : gamma {GAMMA} -> {gm} (PRECISION_RECALL_FRONTIER.md); "
          f"n_families -> {len(built['catalog'])} (frontier gamma=0.40 row: ~623)")
        P(f"  precision      : {HIGH_PRECISION_NOTE['precision_impact']}")
        P(f"                   live residual oversize -> {val['residual_remaining']['oversize']} "
          f"(default 3: MPHOSPH8 + LOC134758618 + MAGE X-array)")
        P(f"  CAVEAT off-orc : {HIGH_PRECISION_NOTE['offoracle_krabznf_cost']}")
        P(f"  CAVEAT oversplt: {HIGH_PRECISION_NOTE['over_split_cost']}")
        P(f"  CAVEAT MAGE flr: {HIGH_PRECISION_NOTE['mage_floor']}")
    if summary is not None:
        P(f"wrote {OUT_TSV}\nwrote {OUT_JSON}")
    P("============================================================================")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="RNA-only family refinement (DEFAULT-ON; opt out with --legacy).")
    ap.add_argument("--legacy", action="store_true",
                    help="opt OUT of the RNA-only refinement: write nothing and recover the "
                         "legacy core_recip>=0.13 shipped catalog (via bench/denovo_families.py)")
    ap.add_argument("--rna-oracle", action="store_true",
                    help="(deprecated no-op; the RNA-only refinement is now the default)")
    ap.add_argument("--no-repeat-gate", action="store_true",
                    help="ablation: DISABLE just the repeat-hub gate (min_shared_mult>=%d cut); "
                         "keeps core+aln+gamma+demote and recovers the pre-repeat-gate catalog "
                         "(also RUSTLE_NO_REPEAT_GATE=1)" % REPEAT_MULT_MIN)
    ap.add_argument("--high-precision", action="store_true",
                    help="HIGH-PRECISION operating point: swap ONLY the gamma-quasi-clique cohesion "
                         "GAMMA=%.2f -> %.2f (PRECISION_RECALL_FRONTIER.md); everything else "
                         "(core/aln thresholds, repeat gate, demote) UNCHANGED. Removes the two "
                         "collapsed-array OVERSIZE blobs (MPHOSPH8, LOC134758618): distinct FP 6->4, "
                         "P_fixed48 0.917, recall held 48/57 (nFam 606 -> 623). HONEST costs: "
                         "off-oracle KRAB-ZNF over-split (gamma>=0.27) + MAGE X-array DNA-only floor "
                         "survives. Also RUSTLE_HIGH_PRECISION=1." % (GAMMA, HIGH_PRECISION_GAMMA))
    ap.add_argument("--no-write", action="store_true",
                    help="run + report but do NOT write outputs (used by the self-check)")
    args = ap.parse_args(argv)
    # DEFAULT-ON: the RNA-only refinement IS the family definition unless the legacy opt-out.
    disabled = args.legacy or os.environ.get("RUSTLE_RNA_ORACLE") == "0"
    if disabled:
        print("legacy: RNA-only refinement DISABLED "
              "(core_recip>=0.13 shipped catalog; run bench/denovo_families.py for the legacy path)")
        return 0
    # DEFAULT-ON repeat-hub gate; --no-repeat-gate / RUSTLE_NO_REPEAT_GATE=1 ablates ONLY the gate.
    repeat_gate = not (args.no_repeat_gate or os.environ.get("RUSTLE_NO_REPEAT_GATE") == "1")
    # HIGH-PRECISION: --high-precision / RUSTLE_HIGH_PRECISION=1 swaps ONLY gamma (0.20 -> 0.40).
    high_precision = args.high_precision or os.environ.get("RUSTLE_HIGH_PRECISION") == "1"
    gamma = HIGH_PRECISION_GAMMA if high_precision else GAMMA
    built, val, summary = run(write=not args.no_write, repeat_gate=repeat_gate, gamma=gamma)
    _report(built, val, summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
