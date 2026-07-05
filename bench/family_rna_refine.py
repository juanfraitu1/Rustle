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

DEFAULT-ON (opt out with --legacy; ablate the gates with --no-*)
---------------------------------------------------------------
The RNA-only refinement (core+aln + repeat-hub gate + antisense/reciprocal-overlap gate + gamma +
recombinant-split gate + multi-repeat-bridge gate + demote) is now the DEFAULT family definition --
runs by default.  Each gate has a matching ablation flag/env (all DEFAULT-ON):
  --no-repeat-gate        / RUSTLE_NO_REPEAT_GATE=1         (single-extreme-edge repeat-hub gate)
  --no-antisense-gate     / RUSTLE_NO_ANTISENSE_GATE=1      (4th gate: edge-level antisense/reciprocal-
                                                            overlap; same-locus opposite-strand nested
                                                            gene; FP_EXCLUSION_DISCRIMINATORS.md)
  --no-split-recombinants / RUSTLE_NO_SPLIT_RECOMBINANTS=1  (recombinant/mosaic path-colinearity split)
  --no-repeat-bridge-gate / RUSTLE_NO_REPEAT_BRIDGE_GATE=1  (3rd gate: family-level multi-repeat-bridge
                                                            conjunction; MULTI_REPEAT_BRIDGE_GATE.md)
--no-repeat-bridge-gate recovers the pre-repeat-bridge catalog (recombinant-split ON, repeat-bridge
OFF) BYTE-IDENTICAL (md5 5e58378a).  The legacy core_recip>=0.13 catalog is recovered with --legacy OR
env RUSTLE_RNA_ORACLE=0 (prints one line, exits 0 without writing; run bench/denovo_families.py for the
legacy catalog).  Because this stage never edits denovo_families.py, the legacy path remains
bit-for-bit reproducible.

MULTI-REPEAT-BRIDGE GATE -- DEPLOYED SCOPE + SENSITIVITY (honest disclosure)
---------------------------------------------------------------------------
The gate is a PREDICATE ("CUT a family iff ...", MULTI_REPEAT_BRIDGE_GATE.md), so it is deployed
catalog-wide on EVERY multi-copy family -- the standing rule that composes with --high-precision (no
fid roster).  At the conservative T=8/C=2 point this cuts 61 DISCONNECTED families (default 605 -> 573):
the 15 pre-identified roster DISCONNECTED-FP over-merges MULTI_REPEAT_BRIDGE_GATE.md scored ("~13-15/35")
PLUS 46 additional DISCONNECTED families the SAME predicate flags blindly (the committed
gate_sweep_catalog T8C2 point -- NOT a 15-family roster lookup; reported honestly).  Every cut is
DISCONNECTED (no full shared exon) by construction.
  PRECISION (improves; the headline): R_oracle 50/57 = 0.8772 HELD; E_p block purity 0.8694 -> 0.8918
  (+0.0224); distinct-FP blocks 6 -> 4; collapsed-array OVERSIZE blobs MPHOSPH8 + LOC134758618 removed
  (oversize 3 -> 1; the MAGE X-array DNA-only floor survives).  GSTM2/MAGE controls + the 23 real-paralog
  controls: 0 cut (share a full exon => single full-exon component => the gate structurally cannot fire).
  SENSITIVITY COST (catalog-scope; R_oracle and E_p are BLIND to it by construction -- the 57-gene diploid
  oracle and the mega-family-excluded E_p do not see the single-locus glued passengers that the <2-loci
  filter drops): the catalog shrinks, so the conservative cDNA-pair truth moves -- pair-recall
  0.9196 -> 0.9087 (== kept_REAL -5) and DNA component-recovery 182/187 -> 177/182
  (bench/family_level_pr_current.py truth2_dna_loose = the operative sensitivity probe).  Of the 5 lost
  REAL cDNA pairs at the DEPLOYED catalog-scope T8C2, 4 are repeat-glue LABEL-NOISE (Alu-glued no-shared-
  exon pairs, mult 231-503, that the cDNA-loose oracle mislabelled REAL -- arguably CORRECT to cut:
  ASB6+NTMT1, CCDC152+SELENOP, ETFDH+PPID, GATC+SRSF9) and only 1 is near-threshold-divergent
  (LOC109024259+ZNF224, exon-id 0.686, nicked at the strict 0.70 cutoff).  ZNF224 itself KEEPS all 10 of
  its KRAB-ZNF paralogs (ZNF155/221/222/225/230/234/284 + 2 LOCs; only the single divergent partner
  LOC109024259 drops) and RFPL3 is UNTOUCHED (its bridge is mult 5 < T=8).  The FULL ZNF224+RFPL3 nick
  (both near-threshold pairs) occurs only at the NON-deployed aggressive T=5/C=1 max-removal (kept_REAL -10).

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
byte-identical (see bench/test_family_rna_refine.py).  The DEFAULT catalog (all four gates ON) is
md5 548029ad across runs (566 families); --no-antisense-gate recovers the pre-antisense md5 dca64cbd
(573 families) BYTE-IDENTICAL.  An outlier "default" md5 means a gate-ablation flag / env leaked in
(e.g. RUSTLE_HIGH_PRECISION=1 -> gamma=0.40; RUSTLE_NO_ANTISENSE_GATE=1 -> dca64cbd), never a
nondeterministic default.

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
# MULTI-REPEAT-BRIDGE GATE thresholds (bench/MULTI_REPEAT_BRIDGE_GATE.md; do NOT re-derive).  This is
# the 3rd VG-native gate = the FAMILY-LEVEL multi-repeat extension of the single-extreme-EDGE repeat-hub
# gate above.  It CUTS a family iff  (NO full shared exon: >=2 full-exon components AND cross-component
# best-exon-id < REPEAT_BRIDGE_EXON_ID)  AND  (REPEAT-bridged: >= REPEAT_BRIDGE_COUNT_MIN distinct
# cross-component shared VG minimizer nodes each with global multiplicity >= REPEAT_BRIDGE_MULT_MIN)
# AND the SAME-GENE guard (never separate same-gene loci) -> replace the family by its full-exon
# components (drop <2-loci).  This removes the DISCONNECTED repeat-bridge FP class (41% of over-merges;
# sub-20 / distributed family-level bridges) that SURVIVES the single-extreme-edge repeat-hub gate.
# CONSERVATIVE T=8/C=2 (MULTI_REPEAT_BRIDGE_GATE.md: ~13-15/35 removed, P_Ep +0.022, R_oracle 50/57
# held, GSTM2/MAGE spared structurally, 0 genuine paralogs lost).  The aggressive T=5/C=1 removes 24/35
# but nicks ~2 near-0.70 divergent-ZNF/RFPL cross-gene paralogs -- AVOIDED.  The no-full-shared-exon
# conjunct is what lets T drop from the repeat-hub 20 to 8 WITHOUT touching GSTM2 (internal Alu mult 426,
# but a single full-exon component -> gate structurally cannot fire) / MAGE (share a full exon).
REPEAT_BRIDGE_MULT_MIN = 8      # T: per-node canonical-minimizer multiplicity (library-free)
REPEAT_BRIDGE_COUNT_MIN = 2     # C: >= C distinct cross-component shared VG nodes at mult >= T
REPEAT_BRIDGE_EXON_ID = 0.70    # NO full shared exon iff cross-component best-exon-id < this (== ID_THRESH)
# allele DEMOTE thresholds (reused from rna_only_edge_oracle.demote_gene):
DEMOTE_BAL_MIN = 0.90    # balanced_frac >= 0.90  (~0.5 minor-allele = diploid het)
DEMOTE_COPY_MAX = 0.10   # copy_like    <= 0.10  (not ~1/K = a real copy)
# ANTISENSE / RECIPROCAL-OVERLAP GATE thresholds (4th VG-external gate; bench/FP_EXCLUSION_DISCRIMINATORS.md;
# do NOT re-derive).  This is the ONE clean NEW rule surviving the four-axis FP-exclusion investigation: a
# pure GENOME-ARCHITECTURE axis (annotation COORDINATES + STRAND, ZERO sequence similarity), ORTHOGONAL to
# the four exhausted similarity axes (nucleotide core/aln, protein E_p, TE/repeat multiplicity, VG topology).
# A within-family CROSS-GENE edge is a FALSE MERGE when the two GENES occupy the SAME genomic region on
# OPPOSITE strands (sense/antisense or nested gene) -- they cannot be two COPIES of one gene.  CUT the edge iff
#   ga != gb  AND  same_contig(ga,gb)  AND  OPPOSITE genomic strand  AND
#   reciprocal_overlap = overlap_bp / min(span_ga, span_gb) >= ANTISENSE_RECIP_OVERLAP_MIN  AND
#   NEITHER gene is a MEGA-SPAN array locus (span >= MEGA_SPAN_MAX).
# 0.50 is a PRINCIPLED floor: the REAL antisense pair MPDU1/MPDU1-AS1 sits at reciprocal-overlap 0.4855 (just
# below) and is SPARED, while the 9 validated genuine over-merges (RASA1+CCNH, RNASEH2C+KAT5, ARHGEF39+CCDC107,
# HDGFL3+TM6SF1, TRMT10B+EXOSC3, +4) sit at 0.55-1.00 and are CUT -- at 0 confirmed-paralog / 0 diploid-oracle
# collateral (FP_EXCLUSION_DISCRIMINATORS.md).  The MEGA-SPAN guard (500 kb) protects CARDINALITY_ARRAY
# truth-artifacts whose (mis-annotated) gene span spuriously "contains" nested LOCs on the opposite strand --
# e.g. GSTM2's 1.18 Mb array span -- which are ALREADY truth-artifacts, not FPs; the guard leaves GSTM2 (and
# MAGE) byte-identical.  Real protein-coding genes in the FP class are all <100 kb, so 500 kb excludes array
# loci without touching a real gene; and the guard can only ever PREVENT a cut (conservative direction).
ANTISENSE_RECIP_OVERLAP_MIN = 0.50   # reciprocal span-overlap floor (MPDU1/MPDU1-AS1 0.4855 -> spared)
MEGA_SPAN_MAX = 500_000              # gene span >= this = CARDINALITY_ARRAY array locus, excluded from the
                                     # antisense cut (GSTM2 1.18 Mb / MAGE spared; only ever prevents a cut)

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
# ANTISENSE gate = pure GENOME-ARCHITECTURE (annotation coordinates + strand); NOT a similarity/DNA-copy-
# number/homology feature -- orthogonal to the four exhausted similarity axes.  Asserted disjoint from the
# DNA/protein/genome and soft-mask/library forbidden columns below (coordinates+strand are already used by
# the gene_of projection; strand is the one extra annotation column).
ANTISENSE_GATE_FEATURES = ("gene_strand", "gene_span")        # per-gene annotation span + strand (architecture)
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

# ANTISENSE gate GENOME-ARCHITECTURE source: per-gene annotation SPAN + STRAND, extracted from the
# GGO_genomic.gff 'gene' features (col 7 = strand) into gene_meta_strand.tsv (chrom/start/end/strand/
# gene/biotype).  Coordinates + strand ONLY -- no sequence, no DNA copy number, no homology label.  The
# coordinates are byte-consistent with the annotation loader (family_er_pr.load_annot / annot_intervals.tsv,
# same GFF); strand is the column the annotation loader drops.
GENE_STRAND_TSV = next((p for p in (os.path.join(BENCH, "gene_meta_strand.tsv"),
                                    "/home/juanfra/winloci_scratch/gene_meta_strand.tsv")
                        if os.path.exists(p)), os.path.join(BENCH, "gene_meta_strand.tsv"))

OUT_TSV = os.path.join(BENCH, "family_rna_refine.tsv")
OUT_JSON = os.path.join(BENCH, "family_rna_refine.json")

assert abs(GAMMA - 0.20) < 1e-9 and SEED == 0, "gamma/seed drifted from the shipped constants"
assert REPEAT_MULT_MIN == 20, "REPEAT_MULT_MIN drifted from the VG_REPEAT_CATALOG.md tail (20)"
assert abs(HIGH_PRECISION_GAMMA - 0.40) < 1e-9, "HIGH_PRECISION_GAMMA drifted from the frontier point (0.40)"
assert REPEAT_BRIDGE_MULT_MIN == 8 and REPEAT_BRIDGE_COUNT_MIN == 2, \
    "repeat-bridge T/C drifted from the MULTI_REPEAT_BRIDGE_GATE.md conservative point (T=8, C=2)"
assert abs(REPEAT_BRIDGE_EXON_ID - 0.70) < 1e-9, "REPEAT_BRIDGE_EXON_ID drifted from ID_THRESH (0.70)"
assert abs(ANTISENSE_RECIP_OVERLAP_MIN - 0.50) < 1e-9, \
    "ANTISENSE_RECIP_OVERLAP_MIN drifted from the principled 0.50 floor (MPDU1/MPDU1-AS1 0.4855 spared)"
assert MEGA_SPAN_MAX == 500_000, "MEGA_SPAN_MAX drifted from the 500 kb array-locus guard (GSTM2/MAGE spared)"


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
    # antisense gate = pure GENOME-ARCHITECTURE (annotation coordinates + strand); NOT a similarity/DNA
    # copy-number/homology column -- disjoint from the DNA/protein/genome + soft-mask/library forbidden sets.
    anti = set(ANTISENSE_GATE_FEATURES)
    assert anti == {"gene_strand", "gene_span"}, f"antisense-gate feature set drifted: {sorted(anti)}"
    leak_anti = anti & (DNA_FORBIDDEN | LIBRARY_FORBIDDEN)
    assert not leak_anti, f"DNA/soft-mask column in the antisense gate: {sorted(leak_anti)}"


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


# --------------------------------------------------------------------------- antisense gate (genome architecture)
def load_gene_strand():
    """gene -> (chrom, start, end, strand) for the antisense/reciprocal-overlap gate, from the
    external gene_meta_strand.tsv (per-gene annotation SPAN + STRAND, extracted from the GGO_genomic.gff
    'gene' features; col 7 = strand; format chrom<TAB>start<TAB>end<TAB>strand<TAB>gene<TAB>biotype).
    GENOME-ARCHITECTURE ONLY: coordinates + strand -- reads NO sequence / DNA-copy-number / homology
    column.  Deterministic: a single ordered pass keeping the FIRST occurrence per gene symbol (the file
    is coordinate-sorted; the few PAR-region duplicate symbols keep their first, lowest-coordinate locus).
    A gene absent from the file => omitted (no antisense cut for a pair touching it)."""
    if not os.path.exists(GENE_STRAND_TSV):
        raise FileNotFoundError(
            f"antisense gate is DEFAULT-ON but {GENE_STRAND_TSV} is missing; extract it from the GFF "
            f"'gene' features, or ablate with --no-antisense-gate / RUSTLE_NO_ANTISENSE_GATE=1")
    out = {}
    with open(GENE_STRAND_TSV) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, start, end, strand, gene = f[0], f[1], f[2], f[3], f[4]
            if gene in out:
                continue                     # first (lowest-coordinate) occurrence -- deterministic
            try:
                out[gene] = (chrom, int(start), int(end), strand)
            except ValueError:
                continue
    return out


# --------------------------------------------------------------------------- build (RNA-only)
def build_catalog(repeat_gate=True, gamma=GAMMA, split_recombinants=True, repeat_bridge_gate=True,
                  antisense_gate=True):
    """Apply the recall-preserving RNA-only gate + repeat-hub gate + antisense/reciprocal-overlap gate
    + shipped gamma refinement + recombinant-split gate + multi-repeat-bridge gate + allele demote.
    Returns dict with the multi-copy catalog and the RNA-only bookkeeping.  No DNA read here.
    repeat_gate=False ablates ONLY the repeat-hub gate.  antisense_gate=False ablates ONLY the
    antisense/reciprocal-overlap gate (edge-level; when OFF the catalog is BYTE-IDENTICAL to the
    pre-antisense golden).  split_recombinants=False ablates ONLY the recombinant-split gate.
    repeat_bridge_gate=False ablates ONLY the multi-repeat-bridge gate (keeps core+aln+repeat+gamma+
    recombinant-split+demote -> recovers the pre-repeat-bridge catalog).  gamma selects the
    gamma-quasi-clique cohesion: default GAMMA=0.20 (recall-preserving) or HIGH_PRECISION_GAMMA=0.40
    (--high-precision; PRECISION_RECALL_FRONTIER.md).  Nothing else changes."""
    rna_only_guard()

    # ---- RNA features (exact oracle loaders) ----
    gene_of_dn = RO.load_gene_of_dn()            # DN locus -> gene symbol (floored projection)
    pair_core = RO.load_pair_core(gene_of_dn)    # gene-pair -> max core_recip (denovo_family_edges.tsv)
    univ_aln = RO.load_universal_aln()           # gene-pair -> aln_frac (ri_sharedlen_universal.tsv; in_ep IGNORED)
    allele = RO.load_allele()                    # gene -> balanced_frac/copy_like/... (a1_read_consensus_o1.tsv)
    # LIBRARY-FREE repeat-hub multiplicity (VG canonical-minimizer catalog; not re-derived):
    pair_repeat_mult = load_repeat_mult() if repeat_gate else {}
    # GENOME-ARCHITECTURE per-gene span + strand (annotation coordinates only; not re-derived):
    gene_strand = load_gene_strand() if antisense_gate else {}

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

    # ---- antisense / reciprocal-overlap gate: the two GENES occupy the SAME genomic region on
    #      OPPOSITE strands (sense/antisense or nested gene) -> cannot be two copies of one gene
    #      (pure GENOME-ARCHITECTURE: coordinates + strand only; NO sequence).  Cut iff same contig
    #      AND opposite strand AND reciprocal span-overlap >= ANTISENSE_RECIP_OVERLAP_MIN AND neither
    #      gene is a MEGA-SPAN array locus (span >= MEGA_SPAN_MAX; GSTM2/MAGE truth-artifact guard). ----
    def antisense_overlap(k):
        a, b = tuple(k)
        ia = gene_strand.get(a); ib = gene_strand.get(b)
        if ia is None or ib is None:             # a gene span/strand unknown => no antisense cut
            return False
        ca, sa, ea, ta = ia
        cb, sb, eb, tb = ib
        if ca != cb or ta == tb:                 # require SAME contig AND OPPOSITE strand
            return False
        spa, spb = ea - sa, eb - sb
        if spa <= 0 or spb <= 0:
            return False
        if spa >= MEGA_SPAN_MAX or spb >= MEGA_SPAN_MAX:   # MEGA-SPAN array-locus guard (GSTM2/MAGE)
            return False
        ov = min(ea, eb) - max(sa, sb)
        if ov <= 0:
            return False
        return (ov / min(spa, spb)) >= ANTISENSE_RECIP_OVERLAP_MIN

    import networkx as nx
    Gr = nx.Graph(); Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b)

    kept = set()
    kept_pairs, cut_pairs, repeat_cut_pairs, antisense_cut_pairs = set(), set(), set(), set()
    n_dn_within, n_dn_cross_kept, n_dn_cross_cut = 0, 0, 0
    n_dn_cross_cut_repeat, n_dn_cross_cut_antisense = 0, 0
    for u, v in Gr.edges():
        ga, gb = gene_of_dn2.get(u), gene_of_dn2.get(v)
        if ga is None or gb is None or ga == gb:
            kept.add(frozenset((u, v)))          # within-gene / unannotated: never an over-merge, never gated
            n_dn_within += 1
            continue
        k = frozenset((ga, gb))
        keep_ca = core_aln_keep(k)
        repeat_only = repeat_gate and repeat_hub(k)       # passes core+aln but shares ONLY extreme repeat
        antisense_only = antisense_gate and antisense_overlap(k)  # same-locus opposite-strand nested gene
        if keep_ca and not repeat_only and not antisense_only:
            kept.add(frozenset((u, v))); n_dn_cross_kept += 1; kept_pairs.add(k)
        else:
            n_dn_cross_cut += 1; cut_pairs.add(k)
            if keep_ca and repeat_only:               # cut SPECIFICALLY by the repeat-hub gate
                n_dn_cross_cut_repeat += 1; repeat_cut_pairs.add(k)
            if keep_ca and antisense_only:            # cut SPECIFICALLY by the antisense/overlap gate
                n_dn_cross_cut_antisense += 1; antisense_cut_pairs.add(k)

    # ---- shipped gamma-quasi-clique refinement (unchanged operator; gamma threaded, seed=0) ----
    comps = SW.components_from_edges(all_nodes, kept)
    refined = G.refine_families(comps, [tuple(e) for e in kept], genes, gamma, SEED)

    # ---- recombinant-split gate (VG path-colinearity; DEFAULT-ON; --no-split-recombinants ablates) ----
    # Split a family held together by a RECOMBINANT/mosaic bridge (a locus whose exon-chain splits
    # colinearly into >=2 DIFFERENT sub-families, articulation point, best single-neighbour colinear
    # cover < 0.60) into its colinear sub-families.  This is the ONE over-merge class the pairwise-edge
    # graph transitively merges (density-2/3 barbell survives gamma) that only path-level colinearity
    # catches (RECOMBINATION_BRIDGE_DETECTOR.md; the K-frontier recombination obstruction).  Sub-families
    # that fall below the >=2-distinct-loci multi-copy predicate after splitting are dropped (not multi-copy).
    split_info = []
    n_families_split = 0
    if split_recombinants:
        import recombinant_split as RS
        pre_split_n = len(refined)
        refined, split_info = RS.split_families(refined, genes, gene_of_dn2)
        n_families_split = len(split_info)
        refined = [b for b in refined if G.distinct_loci(b, genes) >= 2]

    # ---- multi-repeat-bridge gate (3rd VG-native gate; DEFAULT-ON; --no-repeat-bridge-gate ablates) ----
    # Runs AFTER the recombinant-split stage (on the split catalog, matching how MULTI_REPEAT_BRIDGE_GATE.md
    # measured it on the post-recombinant-split catalog) and BEFORE allele-demote.  CUTS a family iff
    # (>=2 full-exon components AND cross-component best-exon-id < REPEAT_BRIDGE_EXON_ID = NO full shared
    # exon) AND (>= REPEAT_BRIDGE_COUNT_MIN distinct cross-component shared VG nodes each with multiplicity
    # >= REPEAT_BRIDGE_MULT_MIN) AND the same-gene guard -> replace by full-exon components (drop <2-loci).
    # This removes the DISCONNECTED repeat-bridge FP class (sub-20 distributed bridges) that survives the
    # single-extreme-edge repeat-hub gate.  Every cut family is DISCONNECTED (no full shared exon) by
    # construction.  MULTI_REPEAT_BRIDGE_GATE.md: GSTM2/MAGE share a full exon -> single component -> the
    # gate STRUCTURALLY cannot fire (spared at every threshold).
    repeat_bridge_info = []
    n_families_repeat_bridge_split = 0
    if repeat_bridge_gate:
        import multi_repeat_bridge_gate as MRB
        assert abs(MRB.ID_THRESH - REPEAT_BRIDGE_EXON_ID) < 1e-9, \
            "multi_repeat_bridge_gate.ID_THRESH drifted from REPEAT_BRIDGE_EXON_ID (0.70)"
        refined, repeat_bridge_info = MRB.split_families_repeat_bridge(
            refined, genes, gene_of_dn2, REPEAT_BRIDGE_MULT_MIN, REPEAT_BRIDGE_COUNT_MIN)
        n_families_repeat_bridge_split = len(repeat_bridge_info)
        refined = [b for b in refined if G.distinct_loci(b, genes) >= 2]

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
        split_recombinants=split_recombinants, n_families_split=n_families_split, split_info=split_info,
        repeat_bridge_gate=repeat_bridge_gate,
        n_families_repeat_bridge_split=n_families_repeat_bridge_split,
        repeat_bridge_info=repeat_bridge_info,
        antisense_gate=antisense_gate,
        n_dn_edges_total=Gr.number_of_edges(),
        n_dn_within=n_dn_within, n_dn_cross_kept=n_dn_cross_kept, n_dn_cross_cut=n_dn_cross_cut,
        n_dn_cross_cut_repeat=n_dn_cross_cut_repeat,
        n_dn_cross_cut_antisense=n_dn_cross_cut_antisense,
        n_cross_pairs_kept=len(kept_pairs), n_cross_pairs_cut=len(cut_pairs - kept_pairs),
        n_cross_pairs_cut_repeat=len(repeat_cut_pairs - kept_pairs),
        n_cross_pairs_cut_antisense=len(antisense_cut_pairs - kept_pairs),
        antisense_cut_pairs=sorted("+".join(sorted(k)) for k in (antisense_cut_pairs - kept_pairs)),
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
        stage="family_rna_refine (RNA-only recall-preserving refinement + repeat-hub gate + "
              "antisense/reciprocal-overlap gate + recombinant-split gate + multi-repeat-bridge gate; "
              "DEFAULT-ON, opt out --legacy / ablate gates --no-repeat-gate / --no-antisense-gate / "
              "--no-split-recombinants / --no-repeat-bridge-gate)",
        rule=dict(edge="KEEP iff core_recip>=%.2f AND aln_frac>=%.2f AND NOT(min_shared_mult>=%d) "
                       "AND NOT(antisense_recip_overlap>=%.2f)"
                       % (CORE_MIN, ALN_MIN, REPEAT_MULT_MIN, ANTISENSE_RECIP_OVERLAP_MIN),
                  core_recip_min=CORE_MIN, aln_frac_min=ALN_MIN,
                  repeat_gate_enabled=repeat_gate, repeat_mult_min=REPEAT_MULT_MIN,
                  antisense_gate_enabled=built["antisense_gate"],
                  antisense_recip_overlap_min=ANTISENSE_RECIP_OVERLAP_MIN, mega_span_max=MEGA_SPAN_MAX,
                  split_recombinants_enabled=built["split_recombinants"],
                  repeat_bridge_gate_enabled=built["repeat_bridge_gate"],
                  repeat_bridge_mult_min=REPEAT_BRIDGE_MULT_MIN,
                  repeat_bridge_count_min=REPEAT_BRIDGE_COUNT_MIN,
                  repeat_bridge_exon_id=REPEAT_BRIDGE_EXON_ID,
                  gamma=gamma, seed=SEED,
                  demote="balanced_frac>=%.2f AND copy_like<=%.2f" % (DEMOTE_BAL_MIN, DEMOTE_COPY_MAX),
                  demote_balanced_frac_min=DEMOTE_BAL_MIN, demote_copy_like_max=DEMOTE_COPY_MAX),
        n_families=len(catalog),
        recombinant_split=dict(
            enabled=built["split_recombinants"],
            n_families_split=built["n_families_split"],
            split_families=built["split_info"],
            note=("VG path-colinearity split of recombinant/mosaic-bridge over-merges (colinear-cover<0.60, "
                  "articulation point, HIGH-confidence distributed mosaic); RECALL-SAFE (never separates "
                  "same-gene loci; HIGH-confidence only) -> 0 oracle-recall cost. RECOMBINATION_BRIDGE_DETECTOR.md. "
                  "--no-split-recombinants / RUSTLE_NO_SPLIT_RECOMBINANTS=1 ablates.")),
        multi_repeat_bridge=dict(
            enabled=built["repeat_bridge_gate"],
            mult_min=REPEAT_BRIDGE_MULT_MIN, count_min=REPEAT_BRIDGE_COUNT_MIN,
            exon_id=REPEAT_BRIDGE_EXON_ID,
            n_families_cut=built["n_families_repeat_bridge_split"],
            # every cut family is DISCONNECTED (no full shared exon: >=2 full-exon components AND
            # cross-component best-exon-id < exon_id) by construction of the gate predicate
            n_disconnected_removed=built["n_families_repeat_bridge_split"],
            families_cut=sorted(built["repeat_bridge_info"], key=lambda d: d["name"]),
            families_cut_named=[d["name"] for d in
                                sorted(built["repeat_bridge_info"], key=lambda d: d["name"])],
            scope=("catalog-wide PREDICATE (applied to every multi-copy family, not a 15-family roster "
                   "lookup): the T=8/C=2 cut count is the committed gate_sweep_catalog point = the "
                   "~13-15 roster DISCONNECTED-FP over-merges PLUS the additional DISCONNECTED families "
                   "the same predicate flags blindly. R_oracle/E_p precision headline HOLDS/IMPROVES; "
                   "the catalog-scope SENSITIVITY cost (single-locus glued passengers dropped by the "
                   "<2-loci filter) is measured by bench/family_level_pr_current.py truth2_dna_loose "
                   "(pair-recall 0.9196->0.9087, DNA component-recovery 182/187->177/182), to which "
                   "R_oracle (57-gene diploid oracle) and E_p (mega-family-excluded) are blind."),
            note=("3rd VG-native gate = family-level multi-repeat extension of the single-extreme-edge "
                  "repeat-hub gate: CUT iff (>=2 full-exon components AND cross-component best-exon-id < "
                  "%.2f) AND (>= %d distinct cross-component shared VG nodes at multiplicity >= %d) AND "
                  "same-gene guard -> replace by full-exon components (drop <2-loci). Removes the "
                  "DISCONNECTED repeat-bridge FP class (sub-20 distributed bridges) that survives the "
                  "single-extreme-edge repeat-hub gate; GSTM2/MAGE share a full exon -> single component "
                  "-> the gate STRUCTURALLY cannot fire (spared). Conservative T=%d/C=%d point "
                  "(MULTI_REPEAT_BRIDGE_GATE.md: R_oracle 50/57 held, P_Ep +, 0 genuine paralogs lost). "
                  "RECALL-SAFE (never separates same-gene loci). "
                  "--no-repeat-bridge-gate / RUSTLE_NO_REPEAT_BRIDGE_GATE=1 ablates."
                  % (REPEAT_BRIDGE_EXON_ID, REPEAT_BRIDGE_COUNT_MIN, REPEAT_BRIDGE_MULT_MIN,
                     REPEAT_BRIDGE_MULT_MIN, REPEAT_BRIDGE_COUNT_MIN))),
        antisense_overlap_gate=dict(
            enabled=built["antisense_gate"],
            recip_overlap_min=ANTISENSE_RECIP_OVERLAP_MIN, mega_span_max=MEGA_SPAN_MAX,
            n_cross_gene_pairs_cut=built["n_cross_pairs_cut_antisense"],
            pairs_cut=built["antisense_cut_pairs"],
            axis="genome_architecture (annotation coordinates + strand; zero sequence similarity)",
            scope=("catalog-wide EDGE-LEVEL predicate (same placement as the repeat-hub gate): a within-"
                   "family cross-gene DN edge is CUT iff its projected gene pair is same-contig, "
                   "OPPOSITE-strand, reciprocal span-overlap >= %.2f, and NEITHER gene is a >= %d bp "
                   "MEGA-SPAN array locus. Two genes at one locus on opposite strands (sense/antisense or "
                   "nested gene) cannot be two copies of one gene." % (ANTISENSE_RECIP_OVERLAP_MIN,
                                                                       MEGA_SPAN_MAX)),
            note=("4th default-on FP gate = the ONE clean NEW rule from FP_EXCLUSION_DISCRIMINATORS.md, "
                  "orthogonal to the 4 exhausted similarity axes (nucleotide/protein/TE/VG-topology). "
                  "Cuts the 9 validated genuine over-merges (RASA1+CCNH, RNASEH2C+KAT5, ARHGEF39+CCDC107, "
                  "HDGFL3+TM6SF1, TRMT10B+EXOSC3, +4) at 0 confirmed-paralog / 0 diploid-oracle collateral; "
                  "the real antisense pair MPDU1/MPDU1-AS1 sits at reciprocal-overlap 0.4855 (< %.2f) and "
                  "is SPARED; the MEGA-SPAN guard leaves the GSTM2 (1.18 Mb array span) + MAGE "
                  "CARDINALITY_ARRAY truth-artifacts byte-identical. RECALL-NEUTRAL (never a copy-of-one-"
                  "gene edge). --no-antisense-gate / RUSTLE_NO_ANTISENSE_GATE=1 ablates (byte-identical "
                  "to the pre-antisense golden)." % ANTISENSE_RECIP_OVERLAP_MIN)),
        edges=dict(
            n_dn_edges_total=built["n_dn_edges_total"],
            n_dn_within_gene_kept=built["n_dn_within"],
            n_dn_cross_gene_kept=built["n_dn_cross_kept"],
            n_dn_cross_gene_cut=built["n_dn_cross_cut"],
            n_dn_cross_gene_cut_by_repeat_gate=built["n_dn_cross_cut_repeat"],
            n_dn_cross_gene_cut_by_antisense_gate=built["n_dn_cross_cut_antisense"],
            n_cross_gene_pairs_kept=built["n_cross_pairs_kept"],
            n_cross_gene_pairs_cut=built["n_cross_pairs_cut"],
            n_cross_gene_pairs_cut_by_repeat_gate=built["n_cross_pairs_cut_repeat"],
            n_cross_gene_pairs_cut_by_antisense_gate=built["n_cross_pairs_cut_antisense"]),
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
            antisense_gate_features=list(ANTISENSE_GATE_FEATURES),
            demote_features=list(DEMOTE_FEATURES),
            no_dna_in_inference=True,
            repeat_gate_library_free=True,
            no_softmask_in_repeat_gate=True,
            antisense_gate_genome_architecture=True,
            no_dna_in_antisense_gate=True,
            gamma=gamma, seed=SEED),
        inputs=dict(
            edges="bench/denovo_family_edges.tsv",
            aln_frac="bench/ri_sharedlen_universal.tsv",
            repeat_mult="bench/vg_repeat_catalog.tsv (min_shared_mult; library-free VG multiplicity)",
            gene_strand=("%s (per-gene annotation span + strand from the GGO_genomic.gff 'gene' "
                         "features; genome architecture only)" % GENE_STRAND_TSV),
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
def run(write=True, repeat_gate=True, gamma=GAMMA, split_recombinants=True, repeat_bridge_gate=True,
        antisense_gate=True):
    built = build_catalog(repeat_gate=repeat_gate, gamma=gamma, split_recombinants=split_recombinants,
                          repeat_bridge_gate=repeat_bridge_gate, antisense_gate=antisense_gate)
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
    P(f"antisense gate   : {'ON' if built['antisense_gate'] else 'OFF (ablated)'}  "
      f"cross-gene pairs cut={built['n_cross_pairs_cut_antisense']} "
      f"(same-contig opposite-strand recip-overlap>={ANTISENSE_RECIP_OVERLAP_MIN:.2f}, "
      f"mega-span guard {MEGA_SPAN_MAX} bp -> GSTM2/MAGE spared)"
      + ("  " + "; ".join(built['antisense_cut_pairs']) if built['antisense_cut_pairs'] else ""))
    P(f"recombinant split: {'ON' if built['split_recombinants'] else 'OFF (ablated)'}  "
      f"families split={built['n_families_split']}"
      + ("  " + "; ".join("|".join(g for g in si['subfam_genes'] if g != 'NA')
                          for si in built['split_info']) if built['split_info'] else ""))
    P(f"repeat-bridge gate: {'ON' if built['repeat_bridge_gate'] else 'OFF (ablated)'}  "
      f"families cut={built['n_families_repeat_bridge_split']} (all DISCONNECTED, no full shared exon; "
      f"catalog-wide predicate = roster ~13-15 DISCONNECTED-FP + additional DISCONNECTED)"
      + ("  " + "; ".join(si['name'] for si in
                          sorted(built['repeat_bridge_info'], key=lambda d: d['name']))
         if built['repeat_bridge_info'] else ""))
    if built['repeat_bridge_gate']:
        P("  precision HOLDS/IMPROVES (R_oracle 50/57 HELD, E_p 0.8694->0.8918, distinct-FP 6->4, "
          "GSTM2/MAGE/23-controls 0 cut); SENSITIVITY cost (catalog-scope, R_oracle/E_p blind): "
          "pair-recall 0.9196->0.9087, DNA component-recovery 182/187->177/182 "
          "(family_level_pr_current.py truth2_dna_loose)")
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
    ap.add_argument("--no-split-recombinants", action="store_true",
                    help="ablation: DISABLE just the recombinant-split gate (VG path-colinearity split "
                         "of recombinant/mosaic-bridge over-merges, e.g. fid 210 GALNT17|LOC101126070); "
                         "keeps core+aln+repeat+gamma+demote and recovers the pre-split catalog "
                         "(also RUSTLE_NO_SPLIT_RECOMBINANTS=1)")
    ap.add_argument("--no-repeat-bridge-gate", action="store_true",
                    help="ablation: DISABLE just the multi-repeat-bridge gate (3rd VG-native gate; "
                         "family-level multi-repeat extension of the repeat-hub gate -- CUT a family "
                         "iff >=2 full-exon components AND cross-component best-exon-id < %.2f AND >= %d "
                         "cross-component shared VG nodes at multiplicity >= %d AND same-gene guard); "
                         "keeps core+aln+repeat+gamma+split+demote and recovers the pre-repeat-bridge "
                         "catalog (also RUSTLE_NO_REPEAT_BRIDGE_GATE=1)"
                         % (REPEAT_BRIDGE_EXON_ID, REPEAT_BRIDGE_COUNT_MIN, REPEAT_BRIDGE_MULT_MIN))
    ap.add_argument("--no-antisense-gate", action="store_true",
                    help="ablation: DISABLE just the antisense/reciprocal-overlap gate (4th default-on FP "
                         "gate; edge-level GENOME-ARCHITECTURE cut -- CUT a cross-gene edge iff its gene "
                         "pair is same-contig, OPPOSITE-strand, reciprocal span-overlap >= %.2f, and "
                         "neither gene span >= %d bp); keeps core+aln+repeat+gamma+split+bridge+demote and "
                         "recovers the pre-antisense catalog BYTE-IDENTICAL "
                         "(also RUSTLE_NO_ANTISENSE_GATE=1)" % (ANTISENSE_RECIP_OVERLAP_MIN, MEGA_SPAN_MAX))
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
    # DEFAULT-ON recombinant-split gate; --no-split-recombinants / RUSTLE_NO_SPLIT_RECOMBINANTS=1 ablates it.
    split_recombinants = not (args.no_split_recombinants
                              or os.environ.get("RUSTLE_NO_SPLIT_RECOMBINANTS") == "1")
    # DEFAULT-ON multi-repeat-bridge gate; --no-repeat-bridge-gate / RUSTLE_NO_REPEAT_BRIDGE_GATE=1 ablates it.
    repeat_bridge_gate = not (args.no_repeat_bridge_gate
                              or os.environ.get("RUSTLE_NO_REPEAT_BRIDGE_GATE") == "1")
    # DEFAULT-ON antisense/reciprocal-overlap gate; --no-antisense-gate / RUSTLE_NO_ANTISENSE_GATE=1 ablates it.
    antisense_gate = not (args.no_antisense_gate
                          or os.environ.get("RUSTLE_NO_ANTISENSE_GATE") == "1")
    # HIGH-PRECISION: --high-precision / RUSTLE_HIGH_PRECISION=1 swaps ONLY gamma (0.20 -> 0.40).
    high_precision = args.high_precision or os.environ.get("RUSTLE_HIGH_PRECISION") == "1"
    gamma = HIGH_PRECISION_GAMMA if high_precision else GAMMA
    built, val, summary = run(write=not args.no_write, repeat_gate=repeat_gate, gamma=gamma,
                              split_recombinants=split_recombinants, repeat_bridge_gate=repeat_bridge_gate,
                              antisense_gate=antisense_gate)
    _report(built, val, summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
