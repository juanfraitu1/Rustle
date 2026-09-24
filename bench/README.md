# bench/ — the analysis scripts and per-topic reports

Regenerated 2026-09-24 after waves 5–6. **38 scripts**, each part of a current objective result, the Soto replication,
or the shared infrastructure; every retired analysis script is at git tags `notebook-2026-09-23b` / `-23c` and in
`~/Desktop/Rustle_attic/`. Reports (`*.md`) are the per-topic records the ledger and register cite; not pruned.

## Shared infrastructure

| script | what it does |
|---|---|
| `soto/rustlib.py` | Canonical primitives for the bench scripts — ONE implementation of each rule that has been got wrong. |
| `sim_reads.py` | Shared full-length HiFi transcript-read simulator (deterministic). IsoSeq reads are full-length, |

## O1 — RNA-level definition machinery (nested edge-test lattice) and the guided mode with its truths

| script | what it does |
|---|---|
| `layer_order/lattice_check_c2.py` | Independent check of lattice_edges.py's clause-2 GENE-BODY re-implementation: for every PAF pair touching the NPIP and |
| `layer_order/lattice_common.py` | Nested edge-test lattice (NPIP/TBC1D3, human CHM13) — shared paths, loaders and THE LEVEL TESTS. |
| `layer_order/lattice_edges.py` | Nested edge-test lattice, step 1: the unified edge table and the L0 closure. |
| `layer_order/lattice_expr.py` | Nested edge-test lattice, step 2: testis read counts for every node of V (lattice/nodes.tsv). |
| `layer_order/lattice_filtration.py` | Nested edge-test lattice, step 5: identity filtration inside L1 (a threshold filtration = single linkage on FIXED |
| `layer_order/lattice_levels.py` | Nested edge-test lattice, step 3: levels G_0..G_3 (primary tests in lattice_common.tests), their 3-truss (triangle) |
| `layer_order/lattice_report_tables.py` | Nested edge-test lattice, step 6: every table quoted in bench/NESTED_LATTICE_NPIP_TBC1D3.md, regenerated from the |
| `layer_order/lattice_truth.py` | Nested edge-test lattice, step 4: agreement of every level with the truths, side by side with the prior study's |
| `layer_order/lo_analysis.py` | NPIP/TBC1D3 layer order (integration, revised after the 2026-09-16 audit): containment, tournament, enforcement cost, |
| `layer_order/lo_corrected_tables.py` | NPIP/TBC1D3 layer order (integration) — apply the verification fixes to the light/heavy builds. |
| `layer_order/lo_expr_recount.py` | NPIP/TBC1D3 layer order (integration) — testis read counts for EXPR on the corrected universe. |
| `layer_order/soto_map.py` | Map RefSeq CHM13 genes to the Soto et al. 2025 gene-ID convention (CAT CHM13_G* / Liftoff LOFF_G* ids). |
| `guided_pipeline.py` | Guided O1 pipeline with the Addendum T fixes, on a leave-out of an annotated truth table. |
| `adjudicated_truth.py` | Prereg Addendum AK: adjudicated two-annotation ground truth for multi-copy gene families, and its scorer. |
| `annotation_nodes.py` | Prereg Addenda AI/AJ/AK: gene-level node tables from an annotation, for `bench/node_graph_mcl.py prep`. |
| `protein_families.py` | Prereg Addendum AN: protein-space multi-copy gene families, and cross-annotation scoring by CDS overlap. |
| `locus_reads.py` | Counting reads at a locus — the ONE correct way, and the wrong way named so it cannot be reached |
| `heldout_family_score.py` | Score `mcl_families` clusters against symbol-root truth families, per |
| `soto_vs_us_referee.py` | Us vs Soto, scored against a NEUTRAL referee. |
| `rna_truth_from_protein.py` | Build a NON-CIRCULAR RNA-level truth from protein families, and measure the ceiling it implies. |
| `protein_edge_gap.py` | Does a protein-level edge close §6o8's no-edge gap? Per `docs/PREREG_protein_edges_2026-09-20.md`. |
| `mcl_port.py` | Python MCL comparator — now a thin shim over the bit-faithful Rust bin `mcl_port` (§6z3, r1047). |
| `ideal_chromosome_sim.py` | Ideal-scenario chromosome simulation, per `docs/PREREG_ideal_chromosome_sim_2026-09-21.md` |

## O2 — read-level truth (sim + score), excision robustness, hard-locus tool bakeoff (calls + compare), Eichler comparator

| script | what it does |
|---|---|
| `o2_read_truth.py` | O2 read-level truth (docs/PREREG_o2_read_truth_2026-09-23.md), one script, two modes. |
| `o2_excision.py` | PREREG adj/excise: remove copy X from a family, rerun copy_assign (genomic read-star default), follow X's |
| `o2_tool_bakeoff.py` | O2 hard-locus tool bakeoff (docs/PREREG_tool_bakeoff_2026-09-08.md, PREREG hard_locus_bakeoff 5ca5c7e4), one script, |
| `eichler_compare.py` | Eichler-style AS-margin assignment, computed alongside ours and compared. |

## O3 — simulations behind the RNA-only chain (transcript / genomic / shuffled)

| script | what it does |
|---|---|
| `o3_sim_copies.py` | O3 simulations with truth (docs/PREREG_o3_reference_bias_2026-09-23.md arm A, docs/PREREG_o3_rna_only_2026-09-23.md |

## Soto 2025 replication chain and scorer

| script | what it does |
|---|---|
| `soto/famcn_from_wssd.py` | Compute famCN (WSSD read-depth copy number) at ARBITRARY coordinates, replicating Soto's CN leg. |
| `soto/soto_attach_noncoding_members.py` | Attach non-eligible-biotype genes (lncRNA, processed_pseudogene, and similar) to an ALREADY-FORMED |
| `soto/soto_bipartite_match_score.py` | Score a predicted gene->family_id assignment against Soto's own published truth via OPTIMAL BIPARTITE |
| `soto/soto_cluster_dennislab_algorithm.py` | The ACTUAL Dennis-lab family-clustering algorithm, reverse-engineered from their own released code |
| `soto/soto_cluster_from_shared.py` | Steps 5-6 of Soto's clustering (connected components -> famCN MAD split -> family/singleton call), |
| `soto/soto_replicate_clustering.py` | Soto 2025 gene-family clustering, steps 4-6, reimplemented from their STAR Methods verbatim. |
| `soto/soto_replicate_from_sedef.py` | Redo Soto's SD98/shared-exon/famCN replication using a fresh, unmerged, native CHM13 v2.0 SEDEF |
| `soto/soto_score_against_truth.py` | Score a predicted gene->family_id assignment against Soto's own published truth (soto_famCN_S1C.tsv), |

## Reports (80)

| report | title | last change |
|---|---|---|
| `ASSEMBLER_WIDENING.md` | Assembler read-isoform widening + single-exon strand — vs `docs/PREREG_assembler_widening_2026-09-18 | 2026-09-18 |
| `ASSEMBLE_ONLY_MODE.md` | `--assemble-only`: the assembler product, with none of the all-vs-all work | 2026-09-19 |
| `ASSEMBLY_POLISH.md` | Assembly polish: matching StringTie in `--assemble-only` mode (§6p8, 2026-09-19) | 2026-09-23 |
| `CHIMERA_POLICY.md` | `--chimera-policy` — measured against `docs/PREREG_chimera_policy_2026-09-18.md` (md5 57f578ff0c4e2e | 2026-09-18 |
| `CHR16_JUNCTION_MAJORITY_ARM.md` | The chr16 arm `build_spliced_seq_with` demands before `RUSTLE_JUNCTION_MAJORITY` can be a default | 2026-09-18 |
| `CHR20_ASSEMBLER_COMPARISON.md` | Chr20 assembler comparison: ours vs StringTie vs FLAIR (gffcompare + SQANTI3) | 2026-09-16 |
| `CLUSTERING_OPERATOR_BAKEOFF.md` | Clustering-operator bakeoff on the L2 copy graph — DESCRIPTIVE, development families only | 2026-09-18 |
| `COPY_ASSIGNMENT_AND_GATE.md` | Copy Assignment And Gate (consolidated) | 2026-09-01 |
| `CROSS_SPECIES_NPIP_CONJUNCT.md` | Agent 3 of 4 — the certificates | 2026-09-17 |
| `DEFINITIONS_FORMAL.md` | Five Concepts, Four Baselines: Paralog, Segmental Duplication, Multi-Copy Gene Family, Expansion, an | 2026-08-14 |
| `DENOVO_PIPELINE.md` | Denovo Pipeline (consolidated) | 2026-09-01 |
| `FALSE_NEGATIVES.md` | False negatives: what the pipeline misses, and why | 2026-08-26 |
| `FAMILY_CERTIFICATES_NPIP_TBC1D3.md` | Family certificates on NPIP and TBC1D3: where each family is an exact connected component, and with  | 2026-09-17 |
| `FAMILY_DEF.md` | Family Definition (consolidated) | 2026-09-01 |
| `FAMILY_LEVELS_AND_RELATED.md` | Family Levels (RNA/DNA/Protein) & Related Methods (consolidated) | 2026-09-01 |
| `FEX_SWEEP_LORO.md` | Raising `--min-shared-exon-frac` from 0.30 to 0.60 — leave-one-region-out validated | 2026-09-18 |
| `GAP_CLOSED_FRACTION.md` | Gap-closed fraction: read-isoform widening and evidence-based admission floors, measured relative to | 2026-09-17 |
| `GATE_CENSUS_NPIP.md` | Where pass-1 skeletons die, and why the "ceiling" was not a depth limit | 2026-09-18 |
| `GENOME_WIDE_BAKEOFF_2026-09-22.md` | Genome-wide `--assemble-only` bakeoff: ours vs StringTie vs FLAIR | 2026-09-23 |
| `GFFCOMPARE_CHR20_2026_09_19.md` | gffcompare on chr20: ours vs StringTie vs FLAIR, with today's assembly settings | 2026-09-19 |
| `GTF_SECONDARY_POOL.md` | `RUSTLE_GTF_SECONDARY` — admitting secondary alignments into the `--gtf` read pool | 2026-09-18 |
| `IDEAL_WIDENING_K5.md` | Read-isoform widening at k = 5 on the idealized simulated substrate (agent 1 of 2) | 2026-09-17 |
| `ISOSEQ_FLAIR_MECHANISMS.md` | isoseq collapse and FLAIR: as comparison arms, and as mechanisms (§6q6, 2026-09-19) | 2026-09-19 |
| `JUNCTION_AND_READTHROUGH_RULES.md` | Adopting the two proposed definition changes: what each form costs, and why none of them is adoptabl | 2026-09-17 |
| `JUNCTION_MAJORITY_CHR16.md` | Splitting the `build_spliced_seq` bucket — and the chr16 arm the code asked for | 2026-09-18 |
| `LAB_DATASET_BAKEOFF.md` | Against the lab's own isoseq / StringTie / FLAIR runs (§6q7, 2026-09-19) | 2026-09-19 |
| `LATTICE_RULE_STRENGTHENERS.md` | Strengthening the DNA levels of the nested edge-test lattice: five strengtheners, swept, on NPIP and | 2026-09-17 |
| `LAYER_ORDER_NPIP_TBC1D3.md` | Layer order for NPIP and TBC1D3 (human T2T-CHM13): protein P, DNA catalog E1 ("D"), subfamily clades | 2026-09-16 |
| `LOCUS_ASSEMBLY_NPIP.md` | Locus assembly at NPIP — measured against `docs/PREREG_locus_assembly_2026-09-18.md` (md5 02237f9d6d | 2026-09-18 |
| `LOCUS_WIDTH_GAP.md` | Locus formation (node width and composition) vs annotated gene records — agent 1 of 2 | 2026-09-17 |
| `MERGED_LOCI_LAYER.md` | The MERGED-LOCUS layer: real fusions recorded as dual membership, outside the partition | 2026-09-19 |
| `NESTED_LATTICE_NPIP_TBC1D3.md` | Nested edge-test lattice on NPIP and TBC1D3 (human T2T-CHM13) | 2026-09-16 |
| `NODE_CUT_RULE.md` | The NODE CUT rule — measured against `docs/PREREG_node_cut_2026-09-18.md` (md5 2af3393070d6c2ded8db3 | 2026-09-18 |
| `NO_READTHROUGH_COUNTERFACTUAL.md` | A world with no readthroughs: it does not rescue the family definition, and on real data it cannot b | 2026-09-18 |
| `NPIP_DISHUCK_TRUTH_CHECK.md` | NPIP subfamily truth checked against Dishuck et al. 2025 | 2026-09-16 |
| `NPIP_IDEAL_EXPRESSION.md` | Agent 2 of 3 — certificates on the idealized synthetic substrate, plus the RNA-structure conjunct | 2026-09-17 |
| `NPIP_MEMBERSHIP_RESCORE.md` | NPIP membership re-run of the shared-junction conjunct t_J (k = 2, δ = 0) | 2026-09-17 |
| `NPIP_NODE_LADDER.md` | NPIP node ladder: annotated DNA nodes → de novo shared-definition nodes (human CHM13, 2026-09-17) | 2026-09-17 |
| `NPIP_NODE_RESCORE.md` | NPIP node arms re-scored: SCORER v2 (one-to-one copy mapping) + GATE v2 (agreement with guided A1) | 2026-09-17 |
| `NPIP_PER_MEMBER_METRICS.md` | NPIP per-member precision / recall / bipartite matching — human CHM13 (2026-09-16) | 2026-09-16 |
| `NPIP_PER_MEMBER_WHOLECHR16.md` | NPIP per-member metrics on whole-chr16 catalogs, with single-exon mapping variants — human CHM13 (20 | 2026-09-17 |
| `NPIP_READ_GROUP_NODES.md` | NPIP read-group nodes vs shipped de novo nodes (human CHM13, 2026-09-17) | 2026-09-17 |
| `NPIP_SIM_CEILING.md` | NPIP algorithmic-ceiling simulation — vs `docs/PREREG_npip_sim_2026-09-18.md` (md5 17be6b031092b79b3 | 2026-09-18 |
| `PERFORMANCE_AND_IO.md` | Performance And Io (consolidated) | 2026-07-09 |
| `PRECISION_LEVERS_CHR20.md` | Getting more transcripts, or better precision — measured on chr20; and what `-R`/`-Q` actually do | 2026-09-19 |
| `READ_ISOFORM_LOCUS.md` | Read-isoform locus and node-admission floor, pre-registered — agent 1 of 2 | 2026-09-17 |
| `RNA_LEVEL_DEFINITION.md` | An RNA-level family definition scored against RNA-DERIVED truth | 2026-09-18 |
| `SOTO_AS_A_REFINEMENT.md` | Are Soto's families a REFINEMENT of ours? (each Soto family inside one of ours) | 2026-09-18 |
| `SOTO_BENCHMARK.md` | Soto benchmark — RNA copy recovery on the human segmental-duplication catalog | 2026-07-16 |
| `SOTO_FAMILY_LOCUS_FIDELITY_2026-09-22.md` | How faithfully does each Soto family's genes get reproduced by our de novo loci? | 2026-09-22 |
| `SQANTI3_POLISH.md` | SQANTI3 on the polished `--assemble-only` output (§6q5, 2026-09-19) | 2026-09-19 |
| `STRAND_ASYMMETRY_AND_REALISTIC_SIM.md` | Strand fix for the shared-definition node set: four variants on two substrates (agent 1 of 2) | 2026-09-18 |
| `STRAND_V4_SCORER_V3.md` | SCORER v3 + genome-wide family panel: V1 / V3 / V4 on S-IDEAL, U00, U20, U40 (agent 1 of 2) | 2026-09-18 |
| `TBC1D3_GUITART_TRUTH_CORRECTION.md` | TBC1D3 subfamily truth corrected from Guitart 2024 Fig 6B/6C | 2026-09-16 |
| `THEORY.md` | Theory (consolidated) | 2026-09-01 |
| `TPM_AND_LOCUS_BED.md` | TPM in the GTF, and a locus BED that shows how close our loci are to the annotation (§6r5, 2026-09-1 | 2026-09-19 |
| `V4_GORILLA_TIEBREAK.md` | V4 on real gorilla data — the tiebreak (agent 1 of 2) | 2026-09-18 |
| `V4_HELDOUT_NOISE_CALIBRATED.md` | V4 on held-back substrates, with the leader rule's noise floor measured and the family clause calibr | 2026-09-18 |
| `V5_RETIRE_PLACEHOLDER.md` | Verification (independent recompute) — agent 2 of 2 | 2026-09-18 |
| `VALIDATION_AND_STATUS.md` | Validation, Reviews & Objective Status (consolidated) | 2026-09-01 |
| `mechanism/consolidation_divergences.md` | Consolidation divergences (deliverable B) | 2026-07-21 |
| `o1_expanded_family_audit/README.md` | Expanded O1 known-family purity audit | 2026-08-19 |
| `o1_fresh_emission_validation/README.md` | O1 fresh-emission validation | 2026-08-19 |
| `o1_gene_family_audit/README.md` | O1 typed gene-family audit | 2026-08-19 |
| `o1_golga2_subfamily_audit.md` | GOLGA2 versus GOLGA6/8: family or false merge? | 2026-08-26 |
| `o1_outgroup_rooting_poc/README.md` | O1 single-outgroup rooting proof of concept | 2026-08-19 |
| `o1_provenance_witness_prototype/README.md` | O1 duplication-provenance witness prototype | 2026-08-19 |
| `soto/MEETING_2026_08_04.md` | O1 status: what is established, what was refuted | 2026-08-04 |
| `soto/bam_tie_signals.md` | BAM signals for multimapper ties: what AS can't do, what `de` fixes, and what's still unused | 2026-08-14 |
| `soto/candidate18_classification.md` | SD classification of the 18 "no-Soto-homology" RNA candidates | 2026-08-22 |
| `soto/coordinate_version_check.md` | Coordinate-version check: Soto BED vs the A119b BAM (CHM13 v1 vs v2) | 2026-07-23 |
| `soto/dna_ceiling_rna_subset.md` | The T2T DNA graph is the 100% ceiling; RNA is a faithful subset | 2026-07-23 |
| `soto/dna_vs_rna_mode.md` | DNA vs RNA, one engine: reproducing Soto's families from the genome vs from transcripts | 2026-07-27 |
| `soto/exon_homogenized_nonexon_signal.md` | Exon-homogenized copies carry recoverable non-exon signal — PPIAL4 (ID_431) | 2026-07-21 |
| `soto/gate_robustness.md` | Soto gate robustness — do the readthrough / mis-chain gates distort detection? (2026-07-23) | 2026-09-01 |
| `soto/merge_quality_analysis.md` | Over- and under-merging in the DNA family partition — diagnosis and levers | 2026-09-01 |
| `soto/nonexon_rescue.md` | Non-exon-signal rescue POC — 24 K=0-bearing Soto families | 2026-07-21 |
| `soto/precision_recall_audit.md` | Soto precision/recall audit — current binary (2026-07-23) | 2026-07-23 |
| `soto/soto_sensitivity_precision.md` | Soto per-family sensitivity / precision — updated 2026-07-21 | 2026-07-21 |
| `soto/tied_seed_eval.md` | Tied-secondary seeding — Phase-1 benchmark result (2026-07-22) | 2026-07-23 |
