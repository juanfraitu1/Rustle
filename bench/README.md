# bench/ — analysis scripts and per-topic reports

Regenerated 2026-09-23 (wave 4). Every `.md` here is a per-topic report the ledger (`docs/o1_ledger.md`) or the
register cites; scripts are the ones a report, `REPRODUCE.md` or a test names, plus their import closure.
Retired reports and scripts are at the git tags `notebook-2026-09-19/20/23` and in `~/Desktop/Rustle_attic/`.

## Reports

| report | title | last change |
|---|---|---|
| `ASSEMBLER_WIDENING.md` | Assembler read-isoform widening + single-exon strand — vs `docs/PREREG_assembler_widening_2026-09-18.md` (md5  | 2026-09-18 |
| `ASSEMBLE_ONLY_MODE.md` | `--assemble-only`: the assembler product, with none of the all-vs-all work | 2026-09-19 |
| `ASSEMBLY_POLISH.md` | Assembly polish: matching StringTie in `--assemble-only` mode (§6p8, 2026-09-19) | 2026-09-23 |
| `CHIMERA_POLICY.md` | `--chimera-policy` — measured against `docs/PREREG_chimera_policy_2026-09-18.md` (md5 57f578ff0c4e2e4829d7c72b | 2026-09-18 |
| `CHR16_JUNCTION_MAJORITY_ARM.md` | The chr16 arm `build_spliced_seq_with` demands before `RUSTLE_JUNCTION_MAJORITY` can be a default | 2026-09-18 |
| `CHR20_ASSEMBLER_COMPARISON.md` | Chr20 assembler comparison: ours vs StringTie vs FLAIR (gffcompare + SQANTI3) | 2026-09-16 |
| `CLUSTERING_OPERATOR_BAKEOFF.md` | Clustering-operator bakeoff on the L2 copy graph — DESCRIPTIVE, development families only | 2026-09-18 |
| `COPY_ASSIGNMENT_AND_GATE.md` | Copy Assignment And Gate (consolidated) | 2026-09-01 |
| `CROSS_SPECIES_NPIP_CONJUNCT.md` | Agent 3 of 4 — the certificates | 2026-09-17 |
| `DEFINITIONS_FORMAL.md` | Five Concepts, Four Baselines: Paralog, Segmental Duplication, Multi-Copy Gene Family, Expansion, and Referenc | 2026-08-14 |
| `DENOVO_PIPELINE.md` | Denovo Pipeline (consolidated) | 2026-09-01 |
| `FALSE_NEGATIVES.md` | False negatives: what the pipeline misses, and why | 2026-08-26 |
| `FAMILY_CERTIFICATES_NPIP_TBC1D3.md` | Family certificates on NPIP and TBC1D3: where each family is an exact connected component, and with how much m | 2026-09-17 |
| `FAMILY_DEF.md` | Family Definition (consolidated) | 2026-09-01 |
| `FAMILY_LEVELS_AND_RELATED.md` | Family Levels (RNA/DNA/Protein) & Related Methods (consolidated) | 2026-09-01 |
| `FEX_SWEEP_LORO.md` | Raising `--min-shared-exon-frac` from 0.30 to 0.60 — leave-one-region-out validated | 2026-09-18 |
| `GAP_CLOSED_FRACTION.md` | Gap-closed fraction: read-isoform widening and evidence-based admission floors, measured relative to the annot | 2026-09-17 |
| `GATE_CENSUS_NPIP.md` | Where pass-1 skeletons die, and why the "ceiling" was not a depth limit | 2026-09-18 |
| `GENOME_WIDE_BAKEOFF_2026-09-22.md` | Genome-wide `--assemble-only` bakeoff: ours vs StringTie vs FLAIR | 2026-09-23 |
| `GFFCOMPARE_CHR20_2026_09_19.md` | gffcompare on chr20: ours vs StringTie vs FLAIR, with today's assembly settings | 2026-09-19 |
| `GTF_SECONDARY_POOL.md` | `RUSTLE_GTF_SECONDARY` — admitting secondary alignments into the `--gtf` read pool | 2026-09-18 |
| `IDEAL_WIDENING_K5.md` | Read-isoform widening at k = 5 on the idealized simulated substrate (agent 1 of 2) | 2026-09-17 |
| `ISOSEQ_FLAIR_MECHANISMS.md` | isoseq collapse and FLAIR: as comparison arms, and as mechanisms (§6q6, 2026-09-19) | 2026-09-19 |
| `JUNCTION_AND_READTHROUGH_RULES.md` | Adopting the two proposed definition changes: what each form costs, and why none of them is adoptable | 2026-09-17 |
| `JUNCTION_MAJORITY_CHR16.md` | Splitting the `build_spliced_seq` bucket — and the chr16 arm the code asked for | 2026-09-18 |
| `LAB_DATASET_BAKEOFF.md` | Against the lab's own isoseq / StringTie / FLAIR runs (§6q7, 2026-09-19) | 2026-09-19 |
| `LATTICE_RULE_STRENGTHENERS.md` | Strengthening the DNA levels of the nested edge-test lattice: five strengtheners, swept, on NPIP and TBC1D3 | 2026-09-17 |
| `LAYER_ORDER_NPIP_TBC1D3.md` | Layer order for NPIP and TBC1D3 (human T2T-CHM13): protein P, DNA catalog E1 ("D"), subfamily clades C, RNA op | 2026-09-16 |
| `LOCUS_ASSEMBLY_NPIP.md` | Locus assembly at NPIP — measured against `docs/PREREG_locus_assembly_2026-09-18.md` (md5 02237f9d6da31ec694a8 | 2026-09-18 |
| `LOCUS_WIDTH_GAP.md` | Locus formation (node width and composition) vs annotated gene records — agent 1 of 2 | 2026-09-17 |
| `MERGED_LOCI_LAYER.md` | The MERGED-LOCUS layer: real fusions recorded as dual membership, outside the partition | 2026-09-19 |
| `NESTED_LATTICE_NPIP_TBC1D3.md` | Nested edge-test lattice on NPIP and TBC1D3 (human T2T-CHM13) | 2026-09-16 |
| `NODE_CUT_RULE.md` | The NODE CUT rule — measured against `docs/PREREG_node_cut_2026-09-18.md` (md5 2af3393070d6c2ded8db3cd89c0a6dc | 2026-09-18 |
| `NO_READTHROUGH_COUNTERFACTUAL.md` | A world with no readthroughs: it does not rescue the family definition, and on real data it cannot be built wi | 2026-09-18 |
| `NPIP_DISHUCK_TRUTH_CHECK.md` | NPIP subfamily truth checked against Dishuck et al. 2025 | 2026-09-16 |
| `NPIP_IDEAL_EXPRESSION.md` | Agent 2 of 3 — certificates on the idealized synthetic substrate, plus the RNA-structure conjunct | 2026-09-17 |
| `NPIP_MEMBERSHIP_RESCORE.md` | NPIP membership re-run of the shared-junction conjunct t_J (k = 2, δ = 0) | 2026-09-17 |
| `NPIP_NODE_LADDER.md` | NPIP node ladder: annotated DNA nodes → de novo shared-definition nodes (human CHM13, 2026-09-17) | 2026-09-17 |
| `NPIP_NODE_RESCORE.md` | NPIP node arms re-scored: SCORER v2 (one-to-one copy mapping) + GATE v2 (agreement with guided A1) | 2026-09-17 |
| `NPIP_PER_MEMBER_METRICS.md` | NPIP per-member precision / recall / bipartite matching — human CHM13 (2026-09-16) | 2026-09-16 |
| `NPIP_PER_MEMBER_WHOLECHR16.md` | NPIP per-member metrics on whole-chr16 catalogs, with single-exon mapping variants — human CHM13 (2026-09-16) | 2026-09-17 |
| `NPIP_READ_GROUP_NODES.md` | NPIP read-group nodes vs shipped de novo nodes (human CHM13, 2026-09-17) | 2026-09-17 |
| `NPIP_SIM_CEILING.md` | NPIP algorithmic-ceiling simulation — vs `docs/PREREG_npip_sim_2026-09-18.md` (md5 17be6b031092b79b343421fa0dd | 2026-09-18 |
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
| `TPM_AND_LOCUS_BED.md` | TPM in the GTF, and a locus BED that shows how close our loci are to the annotation (§6r5, 2026-09-19) | 2026-09-19 |
| `V4_GORILLA_TIEBREAK.md` | V4 on real gorilla data — the tiebreak (agent 1 of 2) | 2026-09-18 |
| `V4_HELDOUT_NOISE_CALIBRATED.md` | V4 on held-back substrates, with the leader rule's noise floor measured and the family clause calibrated to it | 2026-09-18 |
| `V5_RETIRE_PLACEHOLDER.md` | Verification (independent recompute) — agent 2 of 2 | 2026-09-18 |
| `VALIDATION_AND_STATUS.md` | Validation, Reviews & Objective Status (consolidated) | 2026-09-01 |
| `o1_golga2_subfamily_audit.md` | GOLGA2 versus GOLGA6/8: family or false merge? | 2026-08-26 |

## Scripts

173 top-level scripts (`bench/*.py`), plus subfolders: `crossspecies/`, `layer_order/`, `mechanism/`, `multi_copy_eval/`, `negative_control/`, `o1_expanded_family_audit/`, `o1_fresh_emission_validation/`, `o1_gene_family_audit/`, `o1_outgroup_rooting_poc/`, `o1_provenance_witness_prototype/`, `soto/`, `tandem_attribution/`.
Run any script with `--help` or read its module docstring; `bench/_shared.py` holds the common loaders.
