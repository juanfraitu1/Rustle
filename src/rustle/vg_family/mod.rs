//! Variation-graph family analysis for novel gene-family copies.
//!
//! Builds a per-exon family graph (`family_graph`) over paralog copies and
//! drives structural detectors (mosaic, segdup, hidden_copy, positional) plus
//! k-mer-based novel-copy rescue. See `docs/UNMAPPED_FAMILY_RESCUE.md` for why
//! the aligner misses the reads this module rescues.

pub mod collapse_gate; // O2: admit a COLLAPSED single-rep locus as a multi-copy family (ambiguity test, then chi(H)).
pub mod minimizers; // O1 over-merge-gate FOUNDATION: canonical (k,w)-minimizers (Rust port of vg_repeat_catalog.py `minimizers`; byte-parity tested).
pub mod repeat_catalog; // O1 over-merge-gate: node-multiplicity REPEAT CATALOG (load_skeletons/dn_exons/NodeCatalog; Rust port of vg_repeat_catalog.py catalog-production core; byte-parity tested).
pub mod family_definition; // O1 multi-copy predicate `distinct_loci` (Rust port of genome_family_def.py; byte-parity tested).
pub mod family_loaders; // O1 family-definition LOADERS (meta/annot/edges/families + gene_of projection; Rust port of family_er_pr.py loaders; byte-parity tested).
pub mod edge_oracle; // O1 RNA-only EDGE-ORACLE loaders (pair_core/universal_aln/allele) + allele-demote decision (Rust port of rna_only_edge_oracle.py; byte-parity tested).
pub mod bridge_detector; // O1 over-merge-gate SHARED LIBRARY: exon graph (family_exons/exon_match_tensor/build_graph/colinear_cov/subfams + aln_id HW-DP; Rust port of recombination_bridge_detector.py library core; byte-parity tested).
pub mod recombinant_split; // O1 over-merge-gate: recombinant-bridge SPLIT gate (articulation_points/_recombinant_bridges/split_block; Rust port of recombinant_split.py; byte-parity tested).
pub mod multi_repeat_bridge; // O1 over-merge-gate: multi-repeat-bridge gate (characterize/gate_cut/split_families_repeat_bridge + locus_node_set/load_node_mult; Rust port of multi_repeat_bridge_gate.py WIRED path; byte-parity tested).
pub mod driver; // O1 family-definition DRIVER: build_catalog/write_outputs/load_repeat_mult/apply_demote orchestration (Rust port of family_rna_refine.py; reproduces the shipped catalog md5 dca64cbd).
pub mod family_graph;
pub mod rescue;
pub mod diagnostic;
pub mod positional;
pub mod tandem;
pub mod consensus; // cross-copy consensus error-correction (subtractive precision lever)
pub mod mosaic;
pub mod segdup;
pub mod hidden_copy;
pub mod phasing;
pub mod secondary_index;
pub mod layer2;
pub mod psv_linkage; // Layer-2 "C": within-molecule PSV->junction linkage (PSV-column extraction).
pub mod allele_specific_junctions; // ASJ: junctions whose usage depends on a molecule's het-SNP allele.
pub mod asj_strand_bias; // O3 ASJ analysis layer: StrandOddsRatio (SOR) strand-bias filter over asj_calls.tsv (Rust port of asj_strand_bias.py; reuses O2 noodles indexed-BAM fetch + CIGAR walk + allele_specific_junctions::fisher_exact_2x2; byte-parity tested vs GGO_mm.bam).
pub mod asj_verify; // O3 ASJ analysis layer: confound control (frac_mq0 MAPQ-0 fraction at the anchor + anchor->junction dist + high_confidence) over asj_calls.tsv -> asj_calls_verified.tsv (Rust port of asj_verify.py; reuses the O2 noodles indexed-BAM fetch w/ 600-cap; byte-parity tested vs GGO.bam).
pub mod asj_genetic_core; // O3 ASJ DELIVERABLE: the reproducible 54-call genetic core. Pure-TSV row-aligned join of asj_calls_verified.tsv (high_confidence) + asj_calls_strandbias.tsv (sor/sor_pass) -> 3-filter funnel (transversion 475->120, non-LOC 120->76, SOR-clean 76->54) -> asj_genetic_core.tsv (Rust port of asj_genetic_core.py; byte-parity tested, \r\n csv.DictWriter bytes).
pub mod copy_split; // Joint read-coherence + PSV decomposition into (copy, isoform) units.
pub mod absent_copy; // Admission gate for reference-ABSENT (collapsed) copy candidates.
pub mod copy_assign; // Copy ASSIGNMENT: resolve a read to a known copy via PSV + junction likelihood.
pub mod o2_margin_gate; // O2 MARGIN gate: PSV+junction logL assign-with-margin (Rust port of copy_assign.py::assign_read, the gate o2_vg_visualization.materialize_family calls; byte-parity tested). Distinct from copy_assign.rs (the min_p significance gate).
pub mod o2_columns; // O2 VG-materialization step 2 (crux): per-column paralogous-allele extraction from a BAM on a backbone contig (Rust port of psv_graph_genomewide.py::column_alleles, the pysam-pileup parser; byte-parity tested vs real + adversarial BAMs).
pub mod o2_materialize; // O2 VG-materialization step 3: the pre-alignment LOADERS + copy-dedup materialize_family runs first (dedup_copies/load_families/load_gene_labels/FLAGSHIP; Rust port of psv_graph_genomewide.py::dedup_copies + o2_vg_visualization.py loaders; byte-parity tested).
pub mod recombinant_abstain; // O2 recombinant-read ABSTAIN gate: belongs-to-no-copy bi-copy switch test over the family VG (Rust port of recombinant_abstain.py; byte-parity tested).
pub mod family_rescue; // Family-aware copy RESCUE: borrow-strength POA confirm of under-assembled copies.
pub mod family_detect; // Strand-aware de-novo family DETECTION: loci collapse + kmer prefilter + POA edges.
pub mod family_split; // De-novo family DECOMPOSITION: connected components + weighted-modularity Louvain.
pub mod read_conflict; // OPERATIONAL family criterion: read cross-mapping conflict graph (mutual-mappability).
pub mod denovo_assemble; // Integration: Pass-1 read-coherence skeletons + general-purpose assemble gate.
pub mod denovo_pipeline; // Integration: de-novo family DETECTION driver (pass1->gate->collapse->detect->split).
pub mod copy_assign_pipeline; // Integration: per-read COPY ASSIGNMENT driver (PSV + junction, discover+assign).
pub mod rescue_pipeline; // Integration: family-aware RESCUE thin-locus scan (borrow-strength copy recovery).
pub mod genome_projection; // Liftoff-style famCN copy enumeration (spec §7): project a family consensus onto the genome via in-engine minimap2 to enumerate near-identical genomic copies, recovering K=0 collapses.
pub mod em_copy_assign; // EM copy-assignment CORE: pure e_step/m_step/loglik (max-likelihood soft relaxation of SDA PSV correlation-clustering, Vollger 2019, on the PSV-aware VG); logl/pi wiring to Task 1's ReadEvidence arrives separately.
pub mod readonly_copy_number; // Reference-free per-family copy number (Task R1): chi_h (PSV conflict-structure lower bound, Rust port of family_copy_number.py's copyonly_K) + depth_cn (read-depth E_fam/lambda_global leg, recovers Tier-3 collapsed copies chi_h misses).
pub mod vg_realign; // VG re-align supplement (Task 1): candidate-read selection (is_candidate + RealignParams) for poor-fit/unmapped reads to be re-aligned to O1's family copy-paths.

pub use family_graph::{ExonClass, FamilyGraph, JunctionEdge};
pub use diagnostic::{RescueClass, classify_internal, classify_external, cigar_has_long_indel};  // Task 6.1
pub use secondary_index::{SecondaryAlignment, SecondaryIndex};
