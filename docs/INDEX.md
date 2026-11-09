# docs/ index — every document, one line each (generated 2026-11-08)

Regenerate with `python3 tools/docs_index.py > docs/INDEX.md` after adding or archiving a document. Archived files keep their names; only the directory changed, so grepping a filename still finds it. Verdicts of archived studies are in the register rows that cite them; pre-registrations are kept verbatim (the commit that introduced each one is the proof it preceded the result).

## Living documents (15, undated, maintained in place)

| file | title | last commit | cited by |
|---|---|---|---|
| [`ADVISOR_QUESTIONS.md`](ADVISOR_QUESTIONS.md) | The advisor's standing questions — and what the codebase can answer | 2026-11-08 | — |
| [`DATA.md`](DATA.md) | Datasets — what they are, where they came from, how to rebuild them | 2026-11-08 | — |
| [`LEDGER_GENERAL.md`](LEDGER_GENERAL.md) | Ledger: general (handoffs, pending lists, register drafts, overviews) | 2026-11-08 | — |
| [`LEDGER_O1.md`](LEDGER_O1.md) | Ledger: O1 (family definition, assembler, locus units) | 2026-11-08 | — |
| [`LEDGER_O2.md`](LEDGER_O2.md) | Ledger: O2 (copy assignment of tied reads) | 2026-11-08 | — |
| [`LEDGER_O3.md`](LEDGER_O3.md) | Ledger: O3 (copies the reference lacks) | 2026-11-08 | — |
| [`MODULE_STATUS.md`](MODULE_STATUS.md) | Module status — what is SHIPPED vs what was TRIED | 2026-11-08 | — |
| [`NEGATIVE_RESULTS_REGISTER.md`](NEGATIVE_RESULTS_REGISTER.md) | REGISTER OF NEGATIVE RESULTS (full: archive/2026-11/NEGATIVE_RESULTS_REGISTER.md:1) | 2026-11-08 | — |
| [`PERFORMANCE.md`](PERFORMANCE.md) | Pipeline performance: where the time goes and the switches that change it | 2026-11-06 | — |
| [`README.md`](README.md) | docs/ — what is where | 2026-11-08 | — |
| [`REFERENCE.md`](REFERENCE.md) | REFERENCE — glossary, data locations, and a worked example | 2026-11-07 | — |
| [`SCRIPT_INVENTORY.md`](SCRIPT_INVENTORY.md) | Script inventory (2026-11-02) | 2026-11-08 | — |
| [`THESIS_OBJECTIVES.md`](THESIS_OBJECTIVES.md) | Thesis objectives — the current three | 2026-11-05 | — |
| [`copy_assignment_definition.md`](copy_assignment_definition.md) | Copy assignment (O2): definition, decision rule, abstention, certificate | 2026-11-05 | — |
| [`seeded_family_definition.md`](seeded_family_definition.md) | The seeded multi-copy family: definition, properties, certificate | 2026-11-08 | — |

## Open studies and current handoff (59)

| file | title | last commit | cited by |
|---|---|---|---|
| [`PREREG_collapse_combined_2026-11-08.md`](PREREG_collapse_combined_2026-11-08.md) | PREREG: combined collapse rule (certificate on the copy's own catalog copy OR S/H pile), held out on HG002 (register 158 | 2026-11-08 | register ×1 |
| [`PREREG_collapse_copyphasing_2026-11-08.md`](PREREG_collapse_copyphasing_2026-11-08.md) | PREREG: does copy_phasing's allele-slot certificate detect collapsed copies on GRCh38? (register 1583, 2026-11-08, befor | 2026-11-08 | register ×1 |
| [`PREREG_crossind_copyphasing_2026-11-08.md`](PREREG_crossind_copyphasing_2026-11-08.md) | PREREG: does copy_phasing certify the cross-individual strong loci? (register 1582, 2026-11-08) | 2026-11-08 | register ×1 |
| [`PREREG_sd_add_edges_2026-11-08.md`](PREREG_sd_add_edges_2026-11-08.md) | PREREG row 1581: add weak SD-only edges on top of SD damping (`--sd-add-edges W`), held out on HG002 | 2026-11-08 | register ×1 |
| [`PREREG_sd_damp_span_2026-11-08.md`](PREREG_sd_damp_span_2026-11-08.md) | PREREG row 1580: damp non-SD edges before MCL (`--sd-damp`, exons + span link rule), held out on HG002 and gorilla | 2026-11-08 | register ×1 |
| [`PREREG_prot_sd_edges_2026-11-08.md`](PREREG_prot_sd_edges_2026-11-08.md) | PREREG row 1590: protein edges between SD-linked loci on top of SD damping (STOP at dev, +0.035) | 2026-11-08 | register ×1 |
| [`PREREG_exon_sharing_2026-11-08.md`](PREREG_exon_sharing_2026-11-08.md) | PREREG row 1591: exon sharing vs exon coverage as the edge score (edge level, HG002 PASS +0.368) | 2026-11-08 | register ×1 |
| [`PREREG_exon_sharing_family_2026-11-08.md`](PREREG_exon_sharing_family_2026-11-08.md) | PREREG row 1592: exon sharing as the family edge rule (STOP at dev, +0.0399) | 2026-11-08 | register ×1 |
| [`PREREG_exon_sharing_lowI_2026-11-09.md`](PREREG_exon_sharing_lowI_2026-11-09.md) | PREREG row 1593: exon sharing >= 0.5 at low inflation (HG002 PASS: plain +0.317, damped +0.145) | 2026-11-09 | register ×1 |
| [`PREREG_rna_copynumber_soto_2026-11-09.md`](PREREG_rna_copynumber_soto_2026-11-09.md) | PREREG row 1594: copy_phasing extra-copy calls vs Soto S1E population parCN classes | 2026-11-09 | pending |
| [`PREREG_as_table_parallel_2026-11-07.md`](PREREG_as_table_parallel_2026-11-07.md) | PREREG row 1567: parallel `as_table` scan (ms columns kept) | 2026-11-07 | register ×1 |
| [`PREREG_as_table_shared_2026-11-07.md`](PREREG_as_table_shared_2026-11-07.md) | PREREG row 1568: one best-AS table per BAM, shared by every run on it | 2026-11-07 | register ×1 |
| [`PREREG_assemble_region_threads_2026-11-07.md`](PREREG_assemble_region_threads_2026-11-07.md) | PREREG row 1564: run the assemble stage's region sweep in parallel (`--region-threads`) | 2026-11-07 | register ×1 |
| [`PREREG_crossind_missing_copies_2026-11-07.md`](PREREG_crossind_missing_copies_2026-11-07.md) | PREREG rows 1578-1579: does a second gorilla's Iso-Seq on Jim's T2T show copies the reference lacks? | 2026-11-08 | register ×2 |
| [`PREREG_locus_pieces_2026-11-07.md`](PREREG_locus_pieces_2026-11-07.md) | PREREG row 1571: locus pieces for the families input (genes with no locus of their own) | 2026-11-07 | register ×1 |
| [`PREREG_locus_pieces_attach_2026-11-07.md`](PREREG_locus_pieces_attach_2026-11-07.md) | PREREG row 1573: attach locus pieces to the unchanged families (no re-clustering) | 2026-11-07 | register ×1 |
| [`PREREG_locus_pieces_confined_2026-11-07.md`](PREREG_locus_pieces_confined_2026-11-07.md) | PREREG row 1572: locus pieces confined to loci with family evidence (follow-up to row 1571) | 2026-11-07 | register ×1 |
| [`PREREG_locus_pieces_heldout_hg002_2026-11-07.md`](PREREG_locus_pieces_heldout_hg002_2026-11-07.md) | PREREG row 1576: held-out of the locus-piece attach rules on HG002 (user option 1) | 2026-11-07 | register ×1 |
| [`PREREG_locus_reps_2026-11-07.md`](PREREG_locus_reps_2026-11-07.md) | PREREG row 1574: more than one representative per locus (option 2 for the 71 absorbed pairs) | 2026-11-07 | register ×1 |
| [`PREREG_locus_reps_disjoint_2026-11-07.md`](PREREG_locus_reps_disjoint_2026-11-07.md) | PREREG row 1575: multiple representatives per locus with exon-disjoint groups (fix of row 1574's grouping) | 2026-11-07 | register ×1 |
| [`PREREG_locus_reps_keepall_2026-11-07.md`](PREREG_locus_reps_keepall_2026-11-07.md) | PREREG row 1577: multiple representatives, every attached group kept (no merge-back) | 2026-11-07 | register ×1 |
| [`PREREG_mono_primed_from_sweep_2026-11-07.md`](PREREG_mono_primed_from_sweep_2026-11-07.md) | PREREG row 1566: `--polish-mono-primed` evidence from the sweep, not a second BAM pass | 2026-11-07 | register ×1 |
| [`PREREG_mono_primed_parallel_2026-11-07.md`](PREREG_mono_primed_parallel_2026-11-07.md) | PREREG row 1565: parallel BAM pass for `--polish-mono-primed` (SnapBam::mono_primed) | 2026-11-07 | register ×1 |
| [`PREREG_structure_view_2026-11-07.md`](PREREG_structure_view_2026-11-07.md) | PREREG row 1569: do O1's sequence families also share exon architecture? (a sequence-free second view) | 2026-11-07 | register ×1 |
| [`PREREG_structure_vs_sequence_2026-11-07.md`](PREREG_structure_vs_sequence_2026-11-07.md) | PREREG row 1570: which family definition agrees better with Compara, sequence (O1) or exon architecture? | 2026-11-07 | register ×1 |
| [`PREREG_annotated_only_2026-11-06.md`](PREREG_annotated_only_2026-11-06.md) | PREREG: annotated-only assign (row 1551, 2026-11-06 21:30, user-proposed) | 2026-11-06 | register ×1 |
| [`PREREG_assemble_profile_2026-11-06.md`](PREREG_assemble_profile_2026-11-06.md) | PREREG row 1560: where does the assemble stage's wall go (A119b, singlenode 72 cores) | 2026-11-06 | register ×1 |
| [`PREREG_bam_candidates_2026-11-06.md`](PREREG_bam_candidates_2026-11-06.md) | PREREG: read-star only on BAM-reached copies (row 1557, 2026-11-07 00:10) | 2026-11-06 | register ×1 |
| [`PREREG_bam_stage_2026-11-06.md`](PREREG_bam_stage_2026-11-06.md) | PREREG row 1547: BAM on NVMe workspace instead of NFS home (2026-11-06 10:15, before the copy and run) | 2026-11-06 | register ×1 |
| [`PREREG_budget_shrink_2026-11-06.md`](PREREG_budget_shrink_2026-11-06.md) | PREREG row 1558: release the read-budget estimate once the fetch count is known (`RUSTLE_BUDGET_SHRINK=1`) | 2026-11-06 | register ×1 |
| [`PREREG_censat_exclude_2026-11-06.md`](PREREG_censat_exclude_2026-11-06.md) | PREREG row 1559: exclude families whose copies lie in centromeric satellite (`RUSTLE_EXCLUDE_BED=censat.bed`) | 2026-11-06 | register ×1 |
| [`PREREG_coverage_flow_2026-11-06.md`](PREREG_coverage_flow_2026-11-06.md) | PREREG: coverage-conserving flow for tie resolution (row 1554, 2026-11-06 22:15, user-requested) | 2026-11-06 | register ×1 |
| [`PREREG_dispatch_defer_2026-11-06.md`](PREREG_dispatch_defer_2026-11-06.md) | PREREG: span-first dispatch + budget deferral (row 1548, 2026-11-06 13:40) | 2026-11-06 | register ×1 |
| [`PREREG_graph_assign_2026-11-06.md`](PREREG_graph_assign_2026-11-06.md) | PREREG: read-to-graph assignment pre-check on MCL1 (row 1544) | 2026-11-06 | register ×1 |
| [`PREREG_graph_chain_2026-11-06.md`](PREREG_graph_chain_2026-11-06.md) | PREREG: chain once on a family graph (row 1552, 2026-11-06 20:40, user-requested parallel lever) | 2026-11-06 | register ×1 |
| [`PREREG_mm2fast_2026-11-06.md`](PREREG_mm2fast_2026-11-06.md) | PREREG: mm2-fast as the read-star aligner (row 1549, 2026-11-06 18:45, user-requested) | 2026-11-06 | register ×1 |
| [`PREREG_psv_read_graph_2026-11-06.md`](PREREG_psv_read_graph_2026-11-06.md) | PREREG: reads as nodes, PSV agreement as edges (row 1555, 2026-11-06 22:50, user's advisor's suggestion) | 2026-11-06 | register ×1 |
| [`PREREG_read_overlap_graph_2026-11-06.md`](PREREG_read_overlap_graph_2026-11-06.md) | PREREG: read-overlap graph on the ambiguous pool (row 1553, 2026-11-06 21:50, user-requested) | 2026-11-06 | register ×1 |
| [`PREREG_satgeom_2026-11-06.md`](PREREG_satgeom_2026-11-06.md) | PREREG row 1561: flag satellite-like families by copy geometry alone (bench/satgeom.py) | 2026-11-06 | register ×1 |
| [`PREREG_satseq2_2026-11-06.md`](PREREG_satseq2_2026-11-06.md) | PREREG row 1563: row 1562's rule against repeat-class truth, plus the assign wall it buys | 2026-11-06 | register ×1 |
| [`PREREG_satseq_2026-11-06.md`](PREREG_satseq_2026-11-06.md) | PREREG row 1562: flag satellite families by copy SEQUENCE alone (bench/satseq.py) | 2026-11-06 | register ×1 |
| [`PREREG_sched_goal_2026-11-06.md`](PREREG_sched_goal_2026-11-06.md) | PREREG: targeted scheduling goal (row 1550, 2026-11-06 20:50, user-set goal) | 2026-11-07 | register ×1 |
| [`PREREG_star_pilot_2026-11-06.md`](PREREG_star_pilot_2026-11-06.md) | PREREG row 1545: per-family read-star pilot with early exit (2026-11-06 08:10, before any build or run) | 2026-11-06 | register ×1 |
| [`PREREG_star_unit_2026-11-06.md`](PREREG_star_unit_2026-11-06.md) | PREREG row 1546: unit-form read-star on the hotspot panel (2026-11-06 09:35, before the run) | 2026-11-06 | register ×1 |
| [`PREREG_unannotated_noise_2026-11-06.md`](PREREG_unannotated_noise_2026-11-06.md) | PREREG: are the unannotated families signal or noise? (row 1556, 2026-11-06 23:15, user-requested) | 2026-11-06 | register ×1 |
| [`PREREG_absent_on_draft_2026-11-05.md`](PREREG_absent_on_draft_2026-11-05.md) | PREREG ABSENT-COPY DETECTOR ON A DRAFT GENOME (register 1539, 2026-11-05, before any number) | 2026-11-05 | register ×1 |
| [`PREREG_o2_guard_default_2026-11-05.md`](PREREG_o2_guard_default_2026-11-05.md) | PREREG O2 PRIMARY-ON-FAMILY GUARD DEFAULT UNDER --tie-ms (register 1542, 2026-11-05, before any run) | 2026-11-05 | register ×1 |
| [`PREREG_o2_shape_allele_guards_2026-11-05.md`](PREREG_o2_shape_allele_guards_2026-11-05.md) | PREREG O2 SHAPE GUARD (register 1540) AND ALLELE-SITE GUARD (register 1541), 2026-11-05, before any number | 2026-11-05 | register ×1 |
| [`PREREG_pipeline_speed_2026-11-05.md`](PREREG_pipeline_speed_2026-11-05.md) | PREREG PIPELINE SPEED, GENOME-WIDE END TO END (register 1543, 2026-11-05, before any new number) | 2026-11-06 | register ×1 |
| [`PREREG_psv_pile_2026-11-05.md`](PREREG_psv_pile_2026-11-05.md) | PREREG PSV PILE AT ABSORBING LOCI (register 1538, 2026-11-05, before any number) | 2026-11-05 | register ×1 |
| [`GOAL_family_assembly_2026-11-04.md`](GOAL_family_assembly_2026-11-04.md) | Goal: family-wise assembly (drafted 2026-11-04 from registers 1469-1535) | 2026-11-04 | — |
| [`PREREG_compat_pool_2026-11-04.md`](PREREG_compat_pool_2026-11-04.md) | PREREG r1526: compatible-chain admission for family members with reads but no two identical chains | 2026-11-04 | register ×1 |
| [`PREREG_graph_threading_2026-11-04.md`](PREREG_graph_threading_2026-11-04.md) | PREREG r1532: family graph threading of multimapper reads: a ceiling measurement (user, 2026-11-04) | 2026-11-04 | register ×1 |
| [`PREREG_junction_snap_panel_2026-11-04.md`](PREREG_junction_snap_panel_2026-11-04.md) | PREREG r1527: read-level junction snap (10 bp) for one-chain family misses | 2026-11-04 | register ×1 |
| [`PREREG_missing_copy_hallmarks_2026-11-04.md`](PREREG_missing_copy_hallmarks_2026-11-04.md) | PREREG MISSING-COPY HALLMARKS ON A NON-T2T GENOME (register 1537, 2026-11-04, before any number) | 2026-11-05 | register ×1 |
| [`PREREG_panel_span_guard_2026-11-04.md`](PREREG_panel_span_guard_2026-11-04.md) | PREREG r1535: panel span guard (scorer correction) + the mono layer measured genome-wide (2026-11-04) | 2026-11-04 | register ×1 |
| [`PREREG_pooled_depth_2026-11-04.md`](PREREG_pooled_depth_2026-11-04.md) | PREREG r1531: depth, not rules: does pooling HG002 Kinnex reads reach the family panel's misses? (2026-11-04) | 2026-11-04 | register ×1 |
| [`PREREG_same_symbol_scorer_2026-11-04.md`](PREREG_same_symbol_scorer_2026-11-04.md) | PREREG r1530: pairscore's "false" pairs include same-symbol CAT duplicate ids; correct the scorer, re-score r1528 (2026- | 2026-11-04 | register ×1 |
| [`PREREG_shared_exon_frac_2026-11-04.md`](PREREG_shared_exon_frac_2026-11-04.md) | PREREG r1528: the shared-exon fraction gate and catalog coverage of family members | 2026-11-04 | register ×1 |
| [`PREREG_spliced_no_model_2026-11-04.md`](PREREG_spliced_no_model_2026-11-04.md) | PREREG r1536: why pass 1 emits no spliced model at the 52 B3 members (step 0 = diagnosis; 2026-11-04) | 2026-11-04 | register ×1 |
| [`PREREG_tied_threading_2026-11-04.md`](PREREG_tied_threading_2026-11-04.md) | PREREG r1533: tied-read loci + graph threading, a prototype on the panel (2026-11-04) | 2026-11-04 | register ×1 |
| [`PREREG_transcript_edge_2026-11-04.md`](PREREG_transcript_edge_2026-11-04.md) | PREREG r1529: a transcript-level edge for sibling loci the shared-exon gate rejects (user goal, 2026-11-04) | 2026-11-04 | register ×1 |
| [`PREREG_unique_primary_misses_2026-11-04.md`](PREREG_unique_primary_misses_2026-11-04.md) | PREREG r1534: what the 177 unique-primary no-model members are (step 0 only; 2026-11-04) | 2026-11-04 | register ×1 |
| [`VAULT_METHODS_AUDIT_2026-11-04.md`](VAULT_METHODS_AUDIT_2026-11-04.md) | ThesisVault methods audit, family-wise assembly — 2026-11-04 | 2026-11-04 | — |

## `archive/2026-11/` — 50 closed documents

| file | title | last commit | cited by |
|---|---|---|---|
| [`PREREG_strand_vote_dup_zmw_2026-11-04.md`](archive/2026-11/PREREG_strand_vote_dup_zmw_2026-11-04.md) | PREREG r1525: ship the strand vote; narrow the duplicate rule to different ZMWs | 2026-11-04 | register ×1 |
| [`PREREG_family_mode_tied_seed_2026-11-04.md`](archive/2026-11/PREREG_family_mode_tied_seed_2026-11-04.md) | PREREG r1523: family-wide assembly mode (tied-seed members + fused-model cuts), judged on the r1521 panel | 2026-11-04 | register ×1 |
| [`PREREG_base_family_misses_2026-11-04.md`](archive/2026-11/PREREG_base_family_misses_2026-11-04.md) | PREREG r1524: recover the base-mode family misses (strand vote, pass-1 drop, catalog-family cut) | 2026-11-04 | register ×1 |
| [`PREREG_tandem_cut_neighbour_2026-11-03.md`](archive/2026-11/PREREG_tandem_cut_neighbour_2026-11-03.md) | PREREG TANDEM CUT, NEIGHBOUR BLOCKS (register 1517, 2026-11-03) | 2026-11-04 | register ×1 |
| [`PREREG_tandem_cut_2026-11-03.md`](archive/2026-11/PREREG_tandem_cut_2026-11-03.md) | PREREG TANDEM CUT (register 1516, 2026-11-03) | 2026-11-04 | register ×1 |
| [`PREREG_seed_unspliced_support_2026-11-03.md`](archive/2026-11/PREREG_seed_unspliced_support_2026-11-03.md) | PREREG UNSPLICED SECONDARIES AS SUPPORT ONLY (register 1520, 2026-11-03) | 2026-11-04 | register ×1 |
| [`PREREG_seed_unspliced_primary_2026-11-03.md`](archive/2026-11/PREREG_seed_unspliced_primary_2026-11-03.md) | PREREG SEED UNSPLICED FROM PRIMARIES (register 1518, 2026-11-03) | 2026-11-04 | register ×1 |
| [`PREREG_secondary_usefulness_2026-11-03.md`](archive/2026-11/PREREG_secondary_usefulness_2026-11-03.md) | PREREG SECONDARY USEFULNESS + MULTIMAPPER FAMILY GRAPH (registers 1511-1512, 2026-11-03) | 2026-11-03 | register ×1 |
| [`PREREG_secondary_pool_2026-11-03.md`](archive/2026-11/PREREG_secondary_pool_2026-11-03.md) | PREREG SECONDARY-POOL CANDIDATE COPIES (register 1513, 2026-11-03) | 2026-11-03 | register ×1 |
| [`PREREG_project_models_v2_2026-11-03.md`](archive/2026-11/PREREG_project_models_v2_2026-11-03.md) | PREREG r1522: projected models v2, judged set gated (hotspot excluded / transcribed target) | 2026-11-04 | register ×1 |
| [`PREREG_project_models_2026-11-03.md`](archive/2026-11/PREREG_project_models_2026-11-03.md) | PREREG PROJECTED MODELS, MULTI-COPY LAYER (register 1519, 2026-11-03) | 2026-11-04 | register ×1 |
| [`PREREG_o2_primary_on_family_2026-11-03.md`](archive/2026-11/PREREG_o2_primary_on_family_2026-11-03.md) | PREREG O2 PRIMARY-ON-FAMILY GUARD (register 1514, 2026-11-03) | 2026-11-03 | register ×1 |
| [`PREREG_o2_ms_tie_proof_2026-11-03.md`](archive/2026-11/PREREG_o2_ms_tie_proof_2026-11-03.md) | PREREG O2 MS-TIE PROOF (register 1509, 2026-11-03) | 2026-11-03 | register ×1 |
| [`PREREG_o2_ms_tie_heldout_2026-11-03.md`](archive/2026-11/PREREG_o2_ms_tie_heldout_2026-11-03.md) | PREREG O2 MS-TIE HELD-OUT (register 1510, 2026-11-03) | 2026-11-03 | register ×1 |
| [`PREREG_o2_contested_census_2026-11-03.md`](archive/2026-11/PREREG_o2_contested_census_2026-11-03.md) | PREREG O2 CONTESTED CENSUS (register 1508, 2026-11-03) | 2026-11-03 | register ×1 |
| [`PREREG_isocon_only_2026-11-03.md`](archive/2026-11/PREREG_isocon_only_2026-11-03.md) | PREREG ISOCON-ONLY TRANSCRIPTS (register 1515, 2026-11-03) | 2026-11-04 | register ×1 |
| [`PREREG_family_panel_tools_2026-11-03.md`](archive/2026-11/PREREG_family_panel_tools_2026-11-03.md) | PREREG FAMILY PANEL ACROSS TOOLS (register 1521, 2026-11-03; descriptive, predictions stated) | 2026-11-04 | register ×1 |
| [`PREREG_truth_split_collapsed_2026-11-02.md`](archive/2026-11/PREREG_truth_split_collapsed_2026-11-02.md) | PREREG TRUTH SPLIT COLLAPSED (register 1494, 2026-11-02) | 2026-11-03 | register ×1 |
| [`PREREG_sqanti_filter_2026-11-02.md`](archive/2026-11/PREREG_sqanti_filter_2026-11-02.md) | PREREG SQANTI3 RULES FILTER (register 1499, 2026-11-02, before the run) | 2026-11-03 | register ×1 |
| [`PREREG_shared_exon_split_2026-11-02.md`](archive/2026-11/PREREG_shared_exon_split_2026-11-02.md) | PREREG SHARED-EXON SPLIT (register 1496, 2026-11-02) | 2026-11-03 | register ×2 |
| [`PREREG_readthrough_detect_2026-11-02.md`](archive/2026-11/PREREG_readthrough_detect_2026-11-02.md) | PREREG READTHROUGH DETECT (register 1492, 2026-11-02) | 2026-11-03 | register ×1 |
| [`PREREG_multi_sample_2026-11-02.md`](archive/2026-11/PREREG_multi_sample_2026-11-02.md) | PREREG MULTI-SAMPLE POOL (register 1503, 2026-11-02, before the run) | 2026-11-03 | register ×1 |
| [`PREREG_missing_copy_nont2t_2026-11-02.md`](archive/2026-11/PREREG_missing_copy_nont2t_2026-11-02.md) | PREREG MISSING COPIES ON A NON-T2T GENOME (register 1504, 2026-11-02, before any alignment or number) | 2026-11-03 | register ×1 |
| [`PREREG_fused_loci_2026-11-02.md`](archive/2026-11/PREREG_fused_loci_2026-11-02.md) | PREREG FUSED LOCI (register 1495, 2026-11-02) | 2026-11-03 | register ×1 |
| [`PREREG_final_subchain_2026-11-02.md`](archive/2026-11/PREREG_final_subchain_2026-11-02.md) | PREREG FINAL SUB-CHAIN LAYER (register 1498, 2026-11-02) | 2026-11-03 | register ×1 |
| [`PREREG_end_artifacts_shared_2026-11-02.md`](archive/2026-11/PREREG_end_artifacts_shared_2026-11-02.md) | PREREG END ARTIFACTS, SHARED-EVIDENCE EXEMPTION (register 1501, 2026-11-02, before code or numbers) | 2026-11-03 | register ×1 |
| [`PREREG_end_artifacts_pas_2026-11-02.md`](archive/2026-11/PREREG_end_artifacts_pas_2026-11-02.md) | PREREG END ARTIFACTS, PAS EVIDENCE (register 1502, 2026-11-02, before code or numbers) | 2026-11-03 | register ×1 |
| [`PREREG_end_artifacts_2026-11-02.md`](archive/2026-11/PREREG_end_artifacts_2026-11-02.md) | PREREG END ARTIFACTS (register 1500, 2026-11-02, before computing any arm) | 2026-11-03 | register ×1 |
| [`PREREG_candidates_parallel_families_2026-11-02.md`](archive/2026-11/PREREG_candidates_parallel_families_2026-11-02.md) | PREREG CANDIDATES PARALLEL FAMILIES (register 1506, 2026-11-02) | 2026-11-03 | register ×1 |
| [`PREREG_assign_window_sweep_2026-11-02.md`](archive/2026-11/PREREG_assign_window_sweep_2026-11-02.md) | PREREG ASSIGN WINDOW SWEEP (register 1507, 2026-11-02) | 2026-11-03 | register ×1 |
| [`PREREG_assign_window_fetch_2026-11-02.md`](archive/2026-11/PREREG_assign_window_fetch_2026-11-02.md) | PREREG ASSIGN WINDOW FETCH (register 1505, 2026-11-02) | 2026-11-03 | register ×1 |
| [`ANNOT_CAT_VS_REFSEQ_2026-11-02.md`](archive/2026-11/ANNOT_CAT_VS_REFSEQ_2026-11-02.md) | CHM13 CAT/Liftoff vs NCBI RefSeq (native CHM13): why they differ (2026-11-02, descriptive, no decision rule) | 2026-11-03 | — |
| [`seeded_family_definition.md`](archive/2026-11/seeded_family_definition.md) | The seeded multi-copy family: definition, properties, certificate | 2026-11-01 | — |
| [`o1_ledger.md`](archive/2026-11/o1_ledger.md) | O1 LEDGER — what has and has not worked | 2026-11-01 | — |
| [`copy_assignment_definition.md`](archive/2026-11/copy_assignment_definition.md) | Copy assignment (O2): definition, decision rule, abstention, certificate | 2026-11-01 | — |
| [`THESIS_OBJECTIVES.md`](archive/2026-11/THESIS_OBJECTIVES.md) | Thesis objectives — the current three | 2026-11-01 | — |
| [`REFERENCE.md`](archive/2026-11/REFERENCE.md) | REFERENCE — glossary, data locations, and a worked example | 2026-11-01 | — |
| [`PREREG_sqanti_single_copy_2026-11-01.md`](archive/2026-11/PREREG_sqanti_single_copy_2026-11-01.md) | Pre-registration: SQANTI3 on single-copy genes, our plain assembly vs StringTie, FLAIR and IsoCon (A119b chr20, 2026-11- | 2026-11-03 | — |
| [`PREREG_assembly_ends_2026-11-01.md`](archive/2026-11/PREREG_assembly_ends_2026-11-01.md) | PREREG ASSEMBLY ENDS (2026-11-01, register 1485 / 1486) | 2026-11-03 | register ×3 |
| [`O3_STATUS.md`](archive/2026-11/O3_STATUS.md) | O3 — reference-absent copies: the evidence, and what it does and does not support | 2026-11-01 | — |
| [`NEGATIVE_RESULTS_REGISTER.md`](archive/2026-11/NEGATIVE_RESULTS_REGISTER.md) | REGISTER OF NEGATIVE RESULTS | 2026-11-01 | register ×1720 |
| [`MEMORY_DIGEST.md`](archive/2026-11/MEMORY_DIGEST.md) | Memory digest — the long form of the memory index | 2026-11-01 | — |
| [`LEDGER_O3.md`](archive/2026-11/LEDGER_O3.md) | Ledger: O3 (copies the reference lacks) | 2026-11-01 | — |
| [`LEDGER_O2.md`](archive/2026-11/LEDGER_O2.md) | Ledger: O2 (copy assignment of tied reads) | 2026-11-01 | — |
| [`LEDGER_O1.md`](archive/2026-11/LEDGER_O1.md) | Ledger: O1 (family definition, assembler, locus units) | 2026-11-01 | — |
| [`LEDGER_GENERAL.md`](archive/2026-11/LEDGER_GENERAL.md) | Ledger: general (handoffs, pending lists, register drafts, overviews) | 2026-11-01 | — |
| [`CHAPTER_O1_O2_COMPOSITION_DRAFT.md`](archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md) | O1 → O2 by composition: families of transcribed loci, and assignment among exactly those loci (draft, 2026-09-05) | 2026-11-01 | — |
| [`ADVISOR_QUESTIONS.md`](archive/2026-11/ADVISOR_QUESTIONS.md) | The advisor's standing questions — and what the codebase can answer | 2026-11-01 | — |
| [`ACTIVE_WORKING_SET.md`](archive/2026-11/ACTIVE_WORKING_SET.md) | Active working set — what is actually in use (2026-09-19; wave 4 applied 2026-09-23) | 2026-11-01 | — |
| [`STATUS_O2_O3_2026-10-24.md`](archive/2026-11/STATUS_O2_O3_2026-10-24.md) | O2 / O3 status, 2026-10-24 (main; o3-fast merged) | 2026-11-01 | — |

## `archive/2026-10/` — 37 closed documents

| file | title | last commit | cited by |
|---|---|---|---|
| [`MERGED_REMOVED.md`](archive/2026-10/MERGED_REMOVED.md) | Removed 2026-11-08: the 164 merged originals | 2026-11-08 | — |
| [`SPLICED_COPY_SUPPORT_2026-10-04.md`](archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md) | What counts as a FOUND copy — spliced support at the 25 human NPIP copies (and gorilla's 25): the page's "own node" over | 2026-10-07 | — |
| [`PREREG_spliced_copy_support_2026-10-04.md`](archive/2026-10/PREREG_spliced_copy_support_2026-10-04.md) | PREREG — what counts as a FOUND copy: spliced support, not subregion overlap (written before any run, 2026-10-04) | 2026-10-07 | — |
| [`PREREG_locus_representative_rule_2026-10-04.md`](archive/2026-10/PREREG_locus_representative_rule_2026-10-04.md) | PREREG — the locus representative must carry the locus's structure: most-junctions vs most-reads (written before any run | 2026-10-07 | — |
| [`NPIP_CHAIN_COMPARISON_2026-10-04.md`](archive/2026-10/NPIP_CHAIN_COMPARISON_2026-10-04.md) | Is the annotation what is expressed at NPIP? Chain-level comparison of the reads with the CAT and RefSeq models at the 2 | 2026-10-07 | — |
| [`LRPAP1_FAMILY_2026-10-04.md`](archive/2026-10/LRPAP1_FAMILY_2026-10-04.md) | The LRPAP1 family in gorilla: eight copies, one annotated name — what O1, O2 and O3 say about it (2026-10-04, for the ad | 2026-10-07 | — |
| [`LOCUS_REPRESENTATIVE_RULE_2026-10-04.md`](archive/2026-10/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md) | The locus representative, most-junctions (R_J) vs most-reads (R_M), on pre-f1v2 de novo GTFs: R_J carries the copy's str | 2026-10-07 | — |
| [`PSV_CEILING_2026-10-03.md`](archive/2026-10/PSV_CEILING_2026-10-03.md) | PSV ceiling — how many copies are too similar for any PSV to place a read (2026-10-03) | 2026-10-07 | — |
| [`PREREG_o2_tie_outside_2026-10-03.md`](archive/2026-10/PREREG_o2_tie_outside_2026-10-03.md) | PREREG — O2 tie-outside registration at the molecule's best score (register row 1220), written before any run | 2026-10-07 | — |
| [`O3_CANDIDATES_CONTROL_A14_2026-10-03.md`](archive/2026-10/O3_CANDIDATES_CONTROL_A14_2026-10-03.md) | The no-deletion control of `o3_candidates` (Amendment 14): C1' FAILS (35/53), C2' HOLDS (0.54%) — the default flip 1f49d | 2026-10-07 | — |
| [`O3_CANDIDATES_CONSENSUS_DEFECT_2026-10-03.md`](archive/2026-10/O3_CANDIDATES_CONSENSUS_DEFECT_2026-10-03.md) | The consensus defect behind Amendment 14's false flags — read-only diagnosis (Task 5 of the A13 plan), 2026-10-03 | 2026-10-07 | — |
| [`O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md`](archive/2026-10/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md) | `o3_candidates` re-run on the 53-family held-out (Amendment 13 + 13b-13e): A13-1 PASSES, A13-2 PASSES, A13-3 PASSES — de | 2026-10-07 | — |
| [`NPIP_READ_SUPPORT_2026-10-03.md`](archive/2026-10/NPIP_READ_SUPPORT_2026-10-03.md) | Do the reads at an NPIP copy support its transcript? Junction-level check of the A119b primaries at the 25 CAT/Liftoff c | 2026-10-07 | — |
| [`YAG_CHAIN_GORILLA_2026-10-02.md`](archive/2026-10/YAG_CHAIN_GORILLA_2026-10-02.md) | The reference-absent-copy chain on the Y ampliconic genes — HELD-OUT: gorilla OR6737 testis on mGorGor1 `_pri`, 2026-10- | 2026-10-07 | — |
| [`VAULT_METHODS_AUDIT_2026-10-02.md`](archive/2026-10/VAULT_METHODS_AUDIT_2026-10-02.md) | ThesisVault methods audit, addendum — 2026-10-02 | 2026-10-07 | — |
| [`RNA_ALLELE_INHOUSE_CONSENSUS_2026-10-02.md`](archive/2026-10/RNA_ALLELE_INHOUSE_CONSENSUS_2026-10-02.md) | The in-house consensus in place of IsoCon (Amendment 11): FAILS, and why — 2026-10-02 | 2026-10-07 | — |
| [`O3_CANDIDATES_ACCEPTANCE_2026-10-02.md`](archive/2026-10/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md) | `o3_candidates` on the 53-family held-out (Amendment 12): A12-1 FAILS, A12-2 FAILS, A12-3 PASSES — 2026-10-02 | 2026-10-07 | — |
| [`FAMSIM_DESIGN_2026-10-02.md`](archive/2026-10/FAMSIM_DESIGN_2026-10-02.md) | famsim — controlled gene-family simulations that prove the condition they test | 2026-10-07 | — |
| [`YAG_CHAIN_HUMAN_2026-10-01.md`](archive/2026-10/YAG_CHAIN_HUMAN_2026-10-01.md) | The reference-absent-copy chain on the human Y ampliconic genes (DEV: A119b on CHM13v2.0), 2026-10-01 | 2026-10-07 | — |
| [`SQANTI3_CAT_2026-10-01.md`](archive/2026-10/SQANTI3_CAT_2026-10-01.md) | SQANTI3 and the gffcompare scorecard on CAT/Liftoff v2.0 (CAT re-run step 4, 2026-10-01) | 2026-10-07 | — |
| [`RNA_ALLELE_TRUTH_FROZEN_2026-10-01.md`](archive/2026-10/RNA_ALLELE_TRUTH_FROZEN_2026-10-01.md) | Truth frozen for PREREG_rna_allele_haplotype_count_2026-10-01 (before any RNA call) | 2026-10-07 | — |
| [`RNA_ALLELE_REFABSENT_2026-10-01.md`](archive/2026-10/RNA_ALLELE_REFABSENT_2026-10-01.md) | Real reference-absent copies of KB3781: does the RNA-only chain flag them? (Amendment 10), 2026-10-01 | 2026-10-07 | — |
| [`RNA_ALLELE_MERGE_2026-10-01.md`](archive/2026-10/RNA_ALLELE_MERGE_2026-10-01.md) | Merging a missing copy's new-copy transcripts into one candidate copy, without truth (Amendment 8), 2026-10-01 | 2026-10-07 | — |
| [`RNA_ALLELE_LINKING_2026-10-01.md`](archive/2026-10/RNA_ALLELE_LINKING_2026-10-01.md) | Linking IsoCon transcripts to their source locus, held-out on 53 multi-copy families (Amendment 7), 2026-10-01 | 2026-10-07 | — |
| [`RNA_ALLELE_ISOCON_2026-10-01.md`](archive/2026-10/RNA_ALLELE_ISOCON_2026-10-01.md) | IsoCon (reference-free) on KB3781: simulation and real fibroblast reads (Amendments 1-3), 2026-10-01 | 2026-10-07 | — |
| [`RNA_ALLELE_HELDOUT_2026-10-01.md`](archive/2026-10/RNA_ALLELE_HELDOUT_2026-10-01.md) | Held-out test: IsoCon's reference-absent transcripts as extra copies (Amendment 6), 2026-10-01 | 2026-10-07 | — |
| [`RNA_ALLELE_EXCISION_2026-10-01.md`](archive/2026-10/RNA_ALLELE_EXCISION_2026-10-01.md) | Reference-absent copies from RNA: unmapped reads, the maternal-only NPIP copy, and the NPIPA2 deletion test (Amendment 4 | 2026-10-07 | — |
| [`RNA_ALLELE_CONTROL_2026-10-01.md`](archive/2026-10/RNA_ALLELE_CONTROL_2026-10-01.md) | The no-deletion control of the IsoCon chain (Amendment 9), 2026-10-01 | 2026-10-07 | — |
| [`PREREG_yag_isocon_chain_2026-10-01.md`](archive/2026-10/PREREG_yag_isocon_chain_2026-10-01.md) | PREREG — the RNA-only reference-absent-copy chain on the Y ampliconic genes (written before any run on these families) | 2026-10-07 | — |
| [`PREREG_rna_allele_haplotype_count_2026-10-01.md`](archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md) | Pre-registration: counting haplotype copies from RNA alone ("transcripts that are alleles"), gorilla KB3781, 2026-10-01 | 2026-10-07 | — |
| [`PREREG_npip_read_pool_2026-10-01.md`](archive/2026-10/PREREG_npip_read_pool_2026-10-01.md) | Pre-registration: which alignments build the loci, shown on NPIP (human chr16, A119b), 2026-10-01 | 2026-10-07 | — |
| [`NPIP_READ_POOL_2026-10-01.md`](archive/2026-10/NPIP_READ_POOL_2026-10-01.md) | Which alignments build the loci, on NPIP (human chr16, A119b): result, 2026-10-01 | 2026-10-07 | — |
| [`HELDOUT_FAMILIES_CAT_2026-10-01.md`](archive/2026-10/HELDOUT_FAMILIES_CAT_2026-10-01.md) | Held-out family test re-run on CAT/Liftoff v2.0 (2026-10-01) | 2026-10-07 | — |
| [`COPY_RECOVERY_TOOLS_CAT_2026-10-01.md`](archive/2026-10/COPY_RECOVERY_TOOLS_CAT_2026-10-01.md) | Per-copy recovery of NPIP and TBC1D3, re-scored on the CAT/Liftoff v2.0 truth (human, 2026-10-01) | 2026-10-07 | — |
| [`CAT_RERUN_STEP5_PLAN_2026-10-01.md`](archive/2026-10/CAT_RERUN_STEP5_PLAN_2026-10-01.md) | CAT re-run, steps 5 and 6: plan (written before any CAT result was computed) | 2026-10-07 | — |
| [`CAT_RERUN_PROTOCOL_2026-10-01.md`](archive/2026-10/CAT_RERUN_PROTOCOL_2026-10-01.md) | Re-running the RefSeq benchmarks on CAT/Liftoff v2.0: protocol (2026-10-01) | 2026-10-07 | — |
| [`ANNOTATION_CAT_DEFAULT_2026-10-01.md`](archive/2026-10/ANNOTATION_CAT_DEFAULT_2026-10-01.md) | Default human annotation: CAT/Liftoff on CHM13 v2.0 (2026-10-01) | 2026-10-07 | — |

## `archive/2026-09/` — 201 closed documents

| file | title | last commit | cited by |
|---|---|---|---|
| [`ROADMAP_O1_O2.md`](archive/2026-09/ROADMAP_O1_O2.md) | Roadmap: making O1 slide into O2 (2026-09-05) | 2026-10-07 | — |
| [`O1_O2_LOOSE_ENDS.md`](archive/2026-09/O1_O2_LOOSE_ENDS.md) | O1 → O2 loose ends — the ordered queue (2026-09-05) | 2026-10-07 | — |
| [`O1_O2_COMPOSITION.md`](archive/2026-09/O1_O2_COMPOSITION.md) | O1 and O2 as a composition — the statement of the method | 2026-10-07 | — |
| [`O1_DEFINITION_SWITCH.md`](archive/2026-09/O1_DEFINITION_SWITCH.md) | O1 — the two definitions, side by side, and the record of why a switch is warranted (2026-09-05) | 2026-10-07 | — |
| [`ISOSEQ_COMPARISON_RECIPE.md`](archive/2026-09/ISOSEQ_COMPARISON_RECIPE.md) | The fair comparison: `isoseq collapse` vs ours — cluster recipe | 2026-10-07 | — |
| [`SOTO_VS_OURS_MEETING_2026-09-30.md`](archive/2026-09/SOTO_VS_OURS_MEETING_2026-09-30.md) | Soto 2025 and our families: what is measured (meeting sheet, 2026-09-30) | 2026-10-07 | — |
| [`SOTO_UNIT_COVER_CN_2026-09-30.md`](archive/2026-09/SOTO_UNIT_COVER_CN_2026-09-30.md) | KEY=unitcovercn | 2026-10-07 | — |
| [`SOTO_UNIT_COVER_2026-09-30.md`](archive/2026-09/SOTO_UNIT_COVER_2026-09-30.md) | KEY=unitcover | 2026-10-07 | — |
| [`SOTO_UNION_EDGES_2026-09-30.md`](archive/2026-09/SOTO_UNION_EDGES_2026-09-30.md) | KEY=unionedges | 2026-10-07 | — |
| [`SOTO_SEDEF_EXON_IDENTITY_2026-09-30.md`](archive/2026-09/SOTO_SEDEF_EXON_IDENTITY_2026-09-30.md) |  | 2026-10-07 | — |
| [`SOTO_OUR_HOMOLOGY_2026-09-30.md`](archive/2026-09/SOTO_OUR_HOMOLOGY_2026-09-30.md) | KEY=ourhomology | 2026-10-07 | — |
| [`SOTO_M2_MEETING_EVIDENCE_2026-09-30.md`](archive/2026-09/SOTO_M2_MEETING_EVIDENCE_2026-09-30.md) | Soto vs our families: reproduced and extended (2026-09-30) | 2026-10-07 | — |
| [`PREREG_unit_cover_cn_2026-09-30.md`](archive/2026-09/PREREG_unit_cover_cn_2026-09-30.md) | Pre-registration: finer units for the family cover — duplicon x copy-number class, sliver-free usage (2026-09-30, KEY=un | 2026-10-07 | — |
| [`PREREG_unit_cover_2026-09-30.md`](archive/2026-09/PREREG_unit_cover_2026-09-30.md) | Pre-registration: families as a partition of duplicon units, genes as paths — does it predict Soto's multi-family genes? | 2026-10-07 | — |
| [`PREREG_soto_union_edges_2026-09-30.md`](archive/2026-09/PREREG_soto_union_edges_2026-09-30.md) | Pre-registration: does adding SEDEF links to the exon map-back graph reproduce more of Soto's families? (2026-09-30, KEY | 2026-10-07 | — |
| [`PREREG_soto_sedef_exon_identity_2026-09-30.md`](archive/2026-09/PREREG_soto_sedef_exon_identity_2026-09-30.md) | Pre-registration: SEDEF links restricted to exon pairs that are themselves >= 98% identical (2026-09-30, KEY=sedefexon) | 2026-10-07 | — |
| [`PREREG_soto_cn_vs_duplicon_2026-09-30.md`](archive/2026-09/PREREG_soto_cn_vs_duplicon_2026-09-30.md) | Pre-registration: does Soto's copy-number cut follow duplicon boundaries? (2026-09-30, KEY=cnduplicon) | 2026-10-07 | — |
| [`PREREG_quasiclique_soto_2026-09-30.md`](archive/2026-09/PREREG_quasiclique_soto_2026-09-30.md) | Pre-registration: does our genome-only self-alignment + quasi-clique grouping reconstruct Soto's families? (2026-09-30,  | 2026-10-07 | — |
| [`PREREG_our_homology_s1c_cn_2026-09-30.md`](archive/2026-09/PREREG_our_homology_s1c_cn_2026-09-30.md) | Pre-registration: does OUR homology method find Soto's families when given Soto's copy numbers? (2026-09-30, KEY=ourhomo | 2026-10-07 | — |
| [`PREREG_npip_fusion_duplicon_2026-09-30.md`](archive/2026-09/PREREG_npip_fusion_duplicon_2026-09-30.md) | Pre-registration: are NPIP fusion transcripts duplicon-boundary crossings inside co-duplicated blocks? (2026-09-30, KEY= | 2026-10-07 | — |
| [`PREREG_container_v2_heldout_2026-09-30.md`](archive/2026-09/PREREG_container_v2_heldout_2026-09-30.md) | Pre-registration: container v2 (majority core) on RESERVED substrates: (a) partner-base accessory calls, (b) the RELATIO | 2026-10-07 | — |
| [`PREREG_container_units_v2_dev_2026-09-30.md`](archive/2026-09/PREREG_container_units_v2_dev_2026-09-30.md) | Pre-registration: the UNIT split with the assembler's OWN locus rule, and real detectors in the loop (dev only) | 2026-10-07 | — |
| [`PREREG_container_units_mechanism_2026-09-30.md`](archive/2026-09/PREREG_container_units_mechanism_2026-09-30.md) | Pre-registration: a UNIT split of fusion transcripts for the families input only — the ceiling of the execution under an | 2026-10-07 | — |
| [`PREREG_container_units_heldout_human_2026-09-30.md`](archive/2026-09/PREREG_container_units_heldout_human_2026-09-30.md) | Pre-registration: the UNITS execution of the fused-locus container on RESERVED human substrates (A119b chr1-12, 14, 15,  | 2026-10-07 | — |
| [`PREREG_container_units_heldout_apes_2026-09-30.md`](archive/2026-09/PREREG_container_units_heldout_apes_2026-09-30.md) | Pre-registration: the UNITS execution of the fused-locus container on RESERVED APE substrates (gorilla OR6737 and KB3781 | 2026-10-07 | — |
| [`PREREG_container_units_definition_2026-09-30.md`](archive/2026-09/PREREG_container_units_definition_2026-09-30.md) | Pre-registration: candidate DEFINITIONS of a unit boundary inside fused transcripts (definition study; no pipeline chang | 2026-10-07 | — |
| [`PREREG_container_csm_2026-09-30.md`](archive/2026-09/PREREG_container_csm_2026-09-30.md) | Pre-registration: family-relative CORE-SUPPORT MEMBERSHIP (container v2), a route that needs no fusion detector (dev onl | 2026-10-07 | — |
| [`HANDOFF_SECOND_MACHINE_2026-09-30.md`](archive/2026-09/HANDOFF_SECOND_MACHINE_2026-09-30.md) | Second machine (small disk): reproduce and extend the Soto evidence, then merge (2026-09-30) | 2026-10-07 | — |
| [`CONTAINER_HEADROOM_2026-09-30.md`](archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md) | Container headroom under the current defaults: what a fused-locus container can still find (2026-09-30) | 2026-10-07 | — |
| [`PREREG_soto_reconciliation_2026-09-29.md`](archive/2026-09/PREREG_soto_reconciliation_2026-09-29.md) | Pre-registration: reconciling our pipeline with Soto 2025 — which named choices close the gap (attribution) | 2026-10-07 | register ×1 |
| [`PREREG_soto_parcn_assembly_2026-09-29.md`](archive/2026-09/PREREG_soto_parcn_assembly_2026-09-29.md) | Pre-registration: exact assembly k-mer parCN for Soto 2025's SD98 genes, and their human-vs-ape calls | 2026-10-07 | register ×1 |
| [`PREREG_soto_losses_2026-09-29.md`](archive/2026-09/PREREG_soto_losses_2026-09-29.md) | Pre-registration: attacking the named losses of the Soto 2025 replication, held out by Soto family | 2026-10-07 | register ×1 |
| [`PREREG_soto_famcn_allwssd_2026-09-29.md`](archive/2026-09/PREREG_soto_famcn_allwssd_2026-09-29.md) | Pre-registration: our own WSSD famCN from all local SGDP tracks in Soto's reconciled family recipe | 2026-10-07 | register ×1 |
| [`PREREG_o2_presence_gate_2026-09-29.md`](archive/2026-09/PREREG_o2_presence_gate_2026-09-29.md) | Pre-registration: a DNA PRESENCE GATE for O2 on a matched individual (HG002) | 2026-10-07 | — |
| [`PREREG_o2_power_2026-09-29.md`](archive/2026-09/PREREG_o2_power_2026-09-29.md) | Pre-registration: a per-read power analysis of O2's abstention (KEY=o2_power) | 2026-10-07 | — |
| [`PREREG_o1_cover_growth_2026-09-29.md`](archive/2026-09/PREREG_o1_cover_growth_2026-09-29.md) | Pre-registration: COVER-BY-GROWTH families for O1 (fusion loci may belong to two families), held out genome-wide | 2026-10-07 | register ×1 |
| [`PREREG_lrc_metric_2026-09-29.md`](archive/2026-09/PREREG_lrc_metric_2026-09-29.md) | Pre-registration: %LRC (LRGASP long-read coverage) as a reporting metric — 2026-09-29 | 2026-10-07 | — |
| [`PREREG_ggo_tbc1d3_holdout_2026-09-29.md`](archive/2026-09/PREREG_ggo_tbc1d3_holdout_2026-09-29.md) | Pre-registration — gorilla TBC1D3 as a HELD-OUT family for the NPIP loss hypothesis | 2026-10-07 | — |
| [`PREREG_ggo_npip_variant_sim_2026-09-29.md`](archive/2026-09/PREREG_ggo_npip_variant_sim_2026-09-29.md) | Pre-registration — do the gorilla individuals' own NPIP copy sequences and/or 5′-truncated reads reproduce the real NPIP | 2026-10-07 | — |
| [`PREREG_f1v2_readshare_2026-09-29.md`](archive/2026-09/PREREG_f1v2_readshare_2026-09-29.md) | Pre-registration: F1v2 (F1 + a read-share MINORITY condition on each bridge) held out on two HUMAN libraries | 2026-10-07 | register ×1 |
| [`PREREG_copy_recovery_tools_2026-09-29.md`](archive/2026-09/PREREG_copy_recovery_tools_2026-09-29.md) | Pre-registration — per-copy NPIP and TBC1D3 recovery: our assembler vs the StringTie, FLAIR and isoseq collapse baseline | 2026-10-07 | — |
| [`COPY_RECOVERY_TOOLS_2026-09-29.md`](archive/2026-09/COPY_RECOVERY_TOOLS_2026-09-29.md) | Per-copy recovery of NPIP and TBC1D3: our assembler and the StringTie, FLAIR and isoseq baselines (2026-09-29) | 2026-10-07 | — |
| [`VAULT_METHODS_AUDIT_2026-09-28.md`](archive/2026-09/VAULT_METHODS_AUDIT_2026-09-28.md) | ThesisVault methods audit — 2026-09-28 | 2026-10-07 | — |
| [`SOTO_REPLICATION_STATUS_2026-09-28.md`](archive/2026-09/SOTO_REPLICATION_STATUS_2026-09-28.md) | Soto 2025 replication: status, gaps, and the "our families are different from theirs" evidence | 2026-10-26 | — |
| [`PREREG_rg3_npip_block_2026-09-28.md`](archive/2026-09/PREREG_rg3_npip_block_2026-09-28.md) | Pre-registration: RG3 (regroup after polish) on the NPIP block of a held-out human library, with ape no-change controls | 2026-10-07 | — |
| [`PREREG_fusion_container_sim_2026-09-28.md`](archive/2026-09/PREREG_fusion_container_sim_2026-09-28.md) | Pre-registration: a container for the extra pieces of fused family members, tested on simulated fusions | 2026-10-07 | — |
| [`PREREG_f1_bridge_locus_2026-09-28.md`](archive/2026-09/PREREG_f1_bridge_locus_2026-09-28.md) | Pre-registration: F1 (bridge-aware regrouping) held out on two gorilla samples | 2026-10-07 | register ×1 |
| [`PREREG_completeness_heldout_2026-09-27.md`](archive/2026-09/PREREG_completeness_heldout_2026-09-27.md) | Pre-registration: held-out completeness of `--assemble-only` — `--polish-tes pas-end` judged; tools, `--polish-subchain  | 2026-10-07 | — |
| [`PREREG_complete_transcripts_2026-09-27.md`](archive/2026-09/PREREG_complete_transcripts_2026-09-27.md) | Pre-registration: compatible-containment collapse (arm A), fewer partial transcripts from `--assemble-only` | 2026-10-07 | — |
| [`COMPLETE_TRANSCRIPTS_AND_ENDS_2026-09-27.md`](archive/2026-09/COMPLETE_TRANSCRIPTS_AND_ENDS_2026-09-27.md) | Complete transcripts and read-proven transcript ends: findings (2026-09-27) | 2026-10-07 | — |
| [`PREREG_readthrough_v3_2026-09-26.md`](archive/2026-09/PREREG_readthrough_v3_2026-09-26.md) | Pre-registration: readthrough junction filter v3 (arm R3): v2 with a reads-only majority guard on tier B | 2026-10-07 | — |
| [`PREREG_readthrough_v2_2026-09-26.md`](archive/2026-09/PREREG_readthrough_v2_2026-09-26.md) | Pre-registration: readthrough junction filter v2 (arm R2): an alternative-last-exon exemption and an own-promoter tier | 2026-10-07 | — |
| [`PREREG_locus_representatives_2026-09-26.md`](archive/2026-09/PREREG_locus_representatives_2026-09-26.md) | Pre-registration: locus representatives, one arm: regroup after polish (RG, rule RG3) | 2026-10-07 | — |
| [`READTHROUGH_G50K_AND_LAST_EXON_2026-09-25.md`](archive/2026-09/READTHROUGH_G50K_AND_LAST_EXON_2026-09-25.md) | Readthroughs: minimap2 `-G 50k` vs a "reaches past the last exon" rule | 2026-10-07 | — |
| [`PREREG_tied_read_loop_2026-09-25.md`](archive/2026-09/PREREG_tied_read_loop_2026-09-25.md) | Pre-registration: the closed loop (families -> assign tied reads -> give each read to its copy -> re-assemble) | 2026-10-07 | — |
| [`PREREG_readthrough_ends_representatives_2026-09-25.md`](archive/2026-09/PREREG_readthrough_ends_representatives_2026-09-25.md) | Pre-registration: a read-end readthrough junction filter before locus formation (better locus representatives) | 2026-10-07 | — |
| [`PREREG_protein_attach_2026-09-25.md`](archive/2026-09/PREREG_protein_attach_2026-09-25.md) | Pre-registration: protein attachment, the manual extra-sensitive step of the de novo families | 2026-10-07 | — |
| [`PREREG_liftoff_loci_2026-09-25.md`](archive/2026-09/PREREG_liftoff_loci_2026-09-25.md) | Pre-registration: loci in the Liftoff framework (Figure 8) | 2026-10-15 | — |
| [`PREREG_guided_transcript_comparison_2026-09-25.md`](archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md) | Pre-registration: annotation-guided transcript comparison (figures 1-3, guided path) — 2026-09-25 | 2026-10-07 | — |
| [`PREREG_genome_wide_families_2026-09-25.md`](archive/2026-09/PREREG_genome_wide_families_2026-09-25.md) | Pre-registration: genome-wide family definition on every sample (Figures 6 and 7) | 2026-10-07 | — |
| [`PREREG_genome_wide_copy_assignment_2026-09-25.md`](archive/2026-09/PREREG_genome_wide_copy_assignment_2026-09-25.md) | Pre-registration: genome-wide copy assignment on every sample (Figures 4 and 5) | 2026-10-07 | — |
| [`PREREG_families_copy_table_2026-09-25.md`](archive/2026-09/PREREG_families_copy_table_2026-09-25.md) | Pre-registration: the families stage writes the copy table that copy assignment consumes (parity checks) | 2026-10-07 | — |
| [`O3_WGS_TRIO_CN_2026-09-25.md`](archive/2026-09/O3_WGS_TRIO_CN_2026-09-25.md) | O3 — between-individual copy number: RNA cannot see it, DNA can (2026-09-25) | 2026-10-07 | — |
| [`ADVISOR_JACCARD_ANCHORS_2026-09-25.md`](archive/2026-09/ADVISOR_JACCARD_ANCHORS_2026-09-25.md) | Advisor's proposal: read-Jaccard + multimapper chaining anchors instead of exon sums | 2026-10-07 | — |
| [`PREREG_tandem_copy_sim_2026-09-24.md`](archive/2026-09/PREREG_tandem_copy_sim_2026-09-24.md) | Pre-registration — tandem-copy simulations: what minimap2 and the pipeline do with near-identical adjacent copies | 2026-10-07 | — |
| [`PREREG_single_edge_rule_2026-09-24.md`](archive/2026-09/PREREG_single_edge_rule_2026-09-24.md) | Pre-registration — one edge rule instead of three tiers, and the 60–90% band | 2026-10-07 | — |
| [`PREREG_identity_spectrum_2026-09-24.md`](archive/2026-09/PREREG_identity_spectrum_2026-09-24.md) | Pre-registration — the identity spectrum: which alignment tier finds which paralogue pairs, against Ensembl Compara | 2026-10-07 | — |
| [`PREREG_family_scoped_pool_2026-09-24.md`](archive/2026-09/PREREG_family_scoped_pool_2026-09-24.md) | Pre-registration — family-wise assembly as a SCOPED secondary pool (§6zl), written before any arm is assembled | 2026-10-07 | — |
| [`PREREG_cov_longer_catalog_2026-09-24.md`](archive/2026-09/PREREG_cov_longer_catalog_2026-09-24.md) | Pre-registration — the longer-side coverage floor on the shipped catalog, against external truths (§6zj) | 2026-10-07 | — |
| [`PREREG_streaming_assembly_2026-09-23.md`](archive/2026-09/PREREG_streaming_assembly_2026-09-23.md) | Pre-registration — a streaming `--assemble-only` that a laptop can run genome-wide | 2026-10-07 | — |
| [`PREREG_single_read_and_noncanonical_2026-09-23.md`](archive/2026-09/PREREG_single_read_and_noncanonical_2026-09-23.md) | Pre-registration — can single-read chains and non-canonical junctions be admitted by a rule that holds precision? | 2026-10-07 | — |
| [`PREREG_quant_parity_2026-09-23.md`](archive/2026-09/PREREG_quant_parity_2026-09-23.md) | Pre-registration — is our per-transcript quantification comparable to StringTie / FLAIR / isoseq collapse? | 2026-10-07 | — |
| [`PREREG_polish_refit_true_counts_2026-09-23.md`](archive/2026-09/PREREG_polish_refit_true_counts_2026-09-23.md) | Pre-registration — re-fit the assembly polish on TRUE read counts (`--keep-coordinate-duplicates`) | 2026-10-07 | — |
| [`PREREG_o3_rna_only_2026-09-23.md`](archive/2026-09/PREREG_o3_rna_only_2026-09-23.md) | Pre-registration — O3 as far as RNA can carry it: flag, characterise, screen, hand to DNA | 2026-10-07 | — |
| [`PREREG_o3_reference_bias_2026-09-23.md`](archive/2026-09/PREREG_o3_reference_bias_2026-09-23.md) | Pre-registration — O3: do copies absent from the reference hide as unmapped or multi-mapping reads? | 2026-10-07 | — |
| [`PREREG_o2_read_truth_2026-09-23.md`](archive/2026-09/PREREG_o2_read_truth_2026-09-23.md) | Pre-registration — O2 accuracy with read-level truth on the real family catalogs | 2026-10-07 | — |
| [`PREREG_assembly_precision_levers_2026-09-23.md`](archive/2026-09/PREREG_assembly_precision_levers_2026-09-23.md) | Pre-registration — two precision levers for `--assemble-only`, and whether they become the defaults | 2026-10-07 | — |
| [`PREREG_union_truth_npip_2026-09-22.md`](archive/2026-09/PREREG_union_truth_npip_2026-09-22.md) | Pre-registration — does the guided arm's lead survive a UNION truth, and how much rides on SEDEF? | 2026-10-07 | — |
| [`PREREG_trie_edge_construction_2026-09-22.md`](archive/2026-09/PREREG_trie_edge_construction_2026-09-22.md) | Pre-registration — can the intron-chain certificate CREATE edges minimap2 never proposes? | 2026-10-07 | — |
| [`PREREG_split_graph_admission_2026-09-22.md`](archive/2026-09/PREREG_split_graph_admission_2026-09-22.md) | Pre-registration — do the §6w0 split triggers recover the 406 COLLATERAL evictions? | 2026-10-07 | — |
| [`PREREG_spliced_guided_mode_2026-09-22.md`](archive/2026-09/PREREG_spliced_guided_mode_2026-09-22.md) | Pre-registration — should guided mode have SPLICED and UNSPLICED variants? | 2026-10-07 | — |
| [`PREREG_single_exon_admission_2026-09-22.md`](archive/2026-09/PREREG_single_exon_admission_2026-09-22.md) | Pre-registration — admitting SINGLE-EXON loci whose reads fall entirely inside one exon of a family member | 2026-10-07 | — |
| [`PREREG_semiguided_npip_2026-09-22.md`](archive/2026-09/PREREG_semiguided_npip_2026-09-22.md) | Pre-registration — does a SEMI-GUIDED (DNA-region) node set beat de novo and guided on NPIP? | 2026-10-07 | — |
| [`PREREG_sedef_core_refine_2026-09-22.md`](archive/2026-09/PREREG_sedef_core_refine_2026-09-22.md) | Pre-registration — does consuming SEDEF improve the family definition in the de novo and guided modes? | 2026-10-07 | — |
| [`PREREG_protein_sensitive_mode_2026-09-22.md`](archive/2026-09/PREREG_protein_sensitive_mode_2026-09-22.md) | Pre-registration — a SENSITIVE MODE: protein edges from translated de novo locus representatives | 2026-10-07 | — |
| [`PREREG_protein_false_merge_2026-09-22.md`](archive/2026-09/PREREG_protein_false_merge_2026-09-22.md) | Pre-registration — the FALSE-MERGE measurement that blocks r906's adoption | 2026-10-07 | — |
| [`PREREG_primary_dedupe_2026-09-22.md`](archive/2026-09/PREREG_primary_dedupe_2026-09-22.md) | Pre-registration — `copy_assign` collapses distinct molecules with identical coordinates; fix it and price it | 2026-10-07 | — |
| [`PREREG_overmerge_ceiling_2026-09-22.md`](archive/2026-09/PREREG_overmerge_ceiling_2026-09-22.md) | Pre-registration — what is the CEILING of repairing node-level over-merge? | 2026-10-07 | — |
| [`PREREG_locus_read_pool_2026-09-22.md`](archive/2026-09/PREREG_locus_read_pool_2026-09-22.md) | Pre-registration — which alignments define the INITIAL loci: primaries only, all alignments, or primaries + near-tied se | 2026-10-07 | — |
| [`PREREG_guarded_containment_admission_2026-09-22.md`](archive/2026-09/PREREG_guarded_containment_admission_2026-09-22.md) | Pre-registration — guarded CONTAINMENT at the admission gate (ceiling measurement) | 2026-10-07 | — |
| [`PREREG_full_protein_homology_2026-09-22.md`](archive/2026-09/PREREG_full_protein_homology_2026-09-22.md) | Pre-registration — does requiring FULL protein homology rescue the sensitive mode? | 2026-10-07 | — |
| [`PREREG_external_tool_bakeoff_2026-09-22.md`](archive/2026-09/PREREG_external_tool_bakeoff_2026-09-22.md) | Pre-registration — is our `--assemble-only` step as good as isoseq collapse, StringTie and FLAIR, on the lab's own genom | 2026-10-07 | — |
| [`PREREG_depth_adaptive_k_2026-09-22.md`](archive/2026-09/PREREG_depth_adaptive_k_2026-09-22.md) | Pre-registration — a depth-adaptive locus boundary (fixed QUANTILE instead of fixed RANK) | 2026-10-07 | — |
| [`PREREG_cov_shorter_conjunct_2026-09-22.md`](archive/2026-09/PREREG_cov_shorter_conjunct_2026-09-22.md) | Pre-registration — adapting `RUSTLE_ER_COVERAGE_LONGER_FLOOR` to the CURRENT definition | 2026-10-07 | — |
| [`PREREG_cov_shorter_adoption_2026-09-22.md`](archive/2026-09/PREREG_cov_shorter_adoption_2026-09-22.md) | Pre-registration — does the containment escape improve FAMILIES, not just admission? | 2026-10-07 | — |
| [`PREREG_corroborated_union_rep_2026-09-22.md`](archive/2026-09/PREREG_corroborated_union_rep_2026-09-22.md) | Pre-registration — a CORROBORATED union rep: widening only where the evidence is real | 2026-10-07 | — |
| [`PREREG_containment_poset_2026-09-22.md`](archive/2026-09/PREREG_containment_poset_2026-09-22.md) | Pre-registration — is a POSET a better family object than a partition of an undirected graph? | 2026-10-07 | — |
| [`NODE_FORMATION_2026-09-22.md`](archive/2026-09/NODE_FORMATION_2026-09-22.md) | Node formation — empirical rules for FP/FN (§6z4, 2026-09-22) | 2026-10-07 | — |
| [`NODE_DEFINITION_2026-09-22.md`](archive/2026-09/NODE_DEFINITION_2026-09-22.md) | Node definition — what limits it, what was fixed, and what is closed (§6x3–§6x7) | 2026-10-07 | — |
| [`READ_JUNCTIONS_RESULT_2026-09-21.md`](archive/2026-09/READ_JUNCTIONS_RESULT_2026-09-21.md) | Read-derived junctions for the junction family rule — coverage 5×, F worse | 2026-10-07 | — |
| [`PREREG_read_junctions_2026-09-21.md`](archive/2026-09/PREREG_read_junctions_2026-09-21.md) | Pre-registration — READ-DERIVED junctions for the junction family rule | 2026-10-07 | — |
| [`PREREG_read_bridged_node_merge_2026-09-21.md`](archive/2026-09/PREREG_read_bridged_node_merge_2026-09-21.md) | PREREG — read-bridged node merge (RBM) | 2026-10-07 | — |
| [`PREREG_node_split_2026-09-21.md`](archive/2026-09/PREREG_node_split_2026-09-21.md) | PREREG — read-level node splitting, two triggers tested in parallel | 2026-10-07 | — |
| [`PREREG_neighbourhood_jaccard_2026-09-21.md`](archive/2026-09/PREREG_neighbourhood_jaccard_2026-09-21.md) | Pre-registration — a Jaccard-like metric on NODE CONNECTIVITY, with the size trap controlled | 2026-10-07 | — |
| [`PREREG_multilabel_locus_resolver_2026-09-21.md`](archive/2026-09/PREREG_multilabel_locus_resolver_2026-09-21.md) | PREREG — multi-label locus-to-gene resolution for scoring | 2026-10-07 | — |
| [`PREREG_ideal_chromosome_sim_2026-09-21.md`](archive/2026-09/PREREG_ideal_chromosome_sim_2026-09-21.md) | PREREG — ideal-scenario chromosome simulation | 2026-10-07 | — |
| [`PREREG_five_prime_deficit_2026-09-21.md`](archive/2026-09/PREREG_five_prime_deficit_2026-09-21.md) | Pre-registration — an empirical rule for the missing 5′ of a constructed locus | 2026-10-07 | — |
| [`PREREG_false_negative_rules_2026-09-21.md`](archive/2026-09/PREREG_false_negative_rules_2026-09-21.md) | Pre-registration — empirical rules against FALSE NEGATIVES | 2026-10-07 | — |
| [`PREREG_cover_and_jn_weight_2026-09-21.md`](archive/2026-09/PREREG_cover_and_jn_weight_2026-09-21.md) | PREREG — two candidate improvements to the family definition | 2026-10-07 | — |
| [`PREREG_colinearity_conjunct_2026-09-21.md`](archive/2026-09/PREREG_colinearity_conjunct_2026-09-21.md) | Pre-registration — does block colinearity separate the 76 rejected containment≥0.90 pairs? | 2026-10-07 | — |
| [`NODE_GRAPH_ADMISSION_2026-09-21.md`](archive/2026-09/NODE_GRAPH_ADMISSION_2026-09-21.md) | What is missing for the node definition: it is not the boundary, it is graph admission | 2026-10-07 | — |
| [`NODE_CONSTRUCTION_OVERMERGE_2026-09-21.md`](archive/2026-09/NODE_CONSTRUCTION_OVERMERGE_2026-09-21.md) | Node construction: the de novo defect is OVER-merge, not fragmentation | 2026-10-07 | — |
| [`NEIGHBOURHOOD_JACCARD_RESULT_2026-09-21.md`](archive/2026-09/NEIGHBOURHOOD_JACCARD_RESULT_2026-09-21.md) | Neighbourhood Jaccard — a real connectivity metric, and exactly where it is blind | 2026-10-07 | — |
| [`FALSE_NEGATIVE_RULES_RESULT_2026-09-21.md`](archive/2026-09/FALSE_NEGATIVE_RULES_RESULT_2026-09-21.md) | False negatives: where they are, and why the obvious rules do not transfer | 2026-10-07 | — |
| [`EICHLER_COMPARISON_2026-09-21.md`](archive/2026-09/EICHLER_COMPARISON_2026-09-21.md) | Eichler's AS-margin rule, computed alongside ours | 2026-10-07 | — |
| [`COVER_AND_JN_DEFINITION_2026-09-21.md`](archive/2026-09/COVER_AND_JN_DEFINITION_2026-09-21.md) | Two candidate definitions tested: cover prediction, and J_N as an MCL weight | 2026-10-07 | — |
| [`VG_CONTAINMENT_RESULT_2026-09-20.md`](archive/2026-09/VG_CONTAINMENT_RESULT_2026-09-20.md) | Would a variation graph solve the containment problem? — No, and the number says why | 2026-10-07 | — |
| [`SLIM_REGRESSION_2026-09-20.md`](archive/2026-09/SLIM_REGRESSION_2026-09-20.md) | Does the slimmed codebase still produce the recorded numbers? | 2026-10-07 | — |
| [`PROTEIN_EDGES_RESULT_2026-09-20.md`](archive/2026-09/PROTEIN_EDGES_RESULT_2026-09-20.md) | Protein edges vs the §6o8 no-edge gap — result | 2026-10-07 | — |
| [`PREREG_vg_containment_2026-09-20.md`](archive/2026-09/PREREG_vg_containment_2026-09-20.md) | Pre-registration — does a variation-graph signal separate real containments from repeat/domain ones? | 2026-10-07 | — |
| [`PREREG_protein_edges_2026-09-20.md`](archive/2026-09/PREREG_protein_edges_2026-09-20.md) | Pre-registration — PROTEIN EDGES: can they close §6o8's no-edge gap? | 2026-10-07 | — |
| [`PREREG_operator_heldout_2026-09-20.md`](archive/2026-09/PREREG_operator_heldout_2026-09-20.md) | Pre-registration — the clustering-operator bakeoff, on HELD-OUT substrate | 2026-10-07 | — |
| [`PREREG_junction_family_2026-09-20.md`](archive/2026-09/PREREG_junction_family_2026-09-20.md) | Pre-registration — a family definition based on SHARED SPLICE JUNCTIONS | 2026-10-07 | — |
| [`PREREG_jaccard_sweep_2026-09-20.md`](archive/2026-09/PREREG_jaccard_sweep_2026-09-20.md) | Pre-registration — SIMPLE JACCARD SWEEP vs the shipped rule | 2026-10-07 | — |
| [`PREREG_intron_placeholder_2026-09-20.md`](archive/2026-09/PREREG_intron_placeholder_2026-09-20.md) | Pre-registration — intron placeholders: does preserving exon SPACING separate real containments? | 2026-10-07 | — |
| [`PREREG_heldout_families_2026-09-20.md`](archive/2026-09/PREREG_heldout_families_2026-09-20.md) | Pre-registration — HELD-OUT FAMILIES: does the shipped O1 rule generalise, or is it NPIP-shaped? | 2026-10-07 | — |
| [`PREREG_chained_jaccard_2026-09-20.md`](archive/2026-09/PREREG_chained_jaccard_2026-09-20.md) | Pre-registration — DP CHAINING as the pair score, vs global Jaccard | 2026-10-07 | — |
| [`PREREG_asymmetry_metric_2026-09-20.md`](archive/2026-09/PREREG_asymmetry_metric_2026-09-20.md) | Pre-registration — a pair metric that tolerates length asymmetry without saturating | 2026-10-07 | — |
| [`OPERATOR_BAKEOFF_HELDOUT_2026-09-20.md`](archive/2026-09/OPERATOR_BAKEOFF_HELDOUT_2026-09-20.md) | Cliques, trusses, communities — the operator bakeoff on held-out substrate | 2026-10-07 | — |
| [`O1_FUSION_GAP_2026-09-20.md`](archive/2026-09/O1_FUSION_GAP_2026-09-20.md) | Where O1 is incomplete, and what fusions specifically cost it | 2026-10-07 | — |
| [`JUNCTION_FAMILY_RESULT_2026-09-20.md`](archive/2026-09/JUNCTION_FAMILY_RESULT_2026-09-20.md) | A family definition from splice junctions — it exists, it works, and it cannot be the definition | 2026-10-07 | — |
| [`JACCARD_SWEEP_RESULT_2026-09-20.md`](archive/2026-09/JACCARD_SWEEP_RESULT_2026-09-20.md) | "A simple Jaccard sweep should work" — tested head-to-head | 2026-10-07 | — |
| [`INTRON_PLACEHOLDER_RESULT_2026-09-20.md`](archive/2026-09/INTRON_PLACEHOLDER_RESULT_2026-09-20.md) | Intron placeholders — tested, and the reason they fail closes the containment question | 2026-10-07 | — |
| [`HELDOUT_FAMILIES_RESULT_2026-09-20.md`](archive/2026-09/HELDOUT_FAMILIES_RESULT_2026-09-20.md) | Held-out family test — result | 2026-10-07 | — |
| [`DOMINANT_GAP_RESOLVED_2026-09-20.md`](archive/2026-09/DOMINANT_GAP_RESOLVED_2026-09-20.md) | The dominant gap, resolved: §6o8's 0.052 ceiling is a property of the TRUTH's divergence, not of O1 | 2026-10-07 | — |
| [`CHAINED_JACCARD_RESULT_2026-09-20.md`](archive/2026-09/CHAINED_JACCARD_RESULT_2026-09-20.md) | DP chaining as the pair score — tested, and it cannot help at this level | 2026-10-07 | — |
| [`ASYMMETRY_METRIC_RESULT_2026-09-20.md`](archive/2026-09/ASYMMETRY_METRIC_RESULT_2026-09-20.md) | A pair metric that tolerates length asymmetry — bakeoff result | 2026-10-07 | — |
| [`ASYMMETRY_BLIND_SPOT_RESOLVED_2026-09-20.md`](archive/2026-09/ASYMMETRY_BLIND_SPOT_RESOLVED_2026-09-20.md) | The length-asymmetry blind spot: real, characterised, and correctly rejected | 2026-10-07 | — |
| [`PREREG_assembly_polish_2026-09-19.md`](archive/2026-09/PREREG_assembly_polish_2026-09-19.md) | Pre-registration — assembly-only precision polish (2026-09-19) | 2026-10-07 | — |
| [`PIPELINE_STATE_TEST_2026-09-19.md`](archive/2026-09/PIPELINE_STATE_TEST_2026-09-19.md) | Does the pipeline, as it stands today, fulfil O1, O2 and the advisor's standing questions? | 2026-10-07 | — |
| [`PREREG_npip_sim_2026-09-18.md`](archive/2026-09/PREREG_npip_sim_2026-09-18.md) | Pre-registration — NPIP algorithmic-ceiling simulation (canonical transcripts, distinguishable reads) | 2026-10-07 | — |
| [`PREREG_node_cut_2026-09-18.md`](archive/2026-09/PREREG_node_cut_2026-09-18.md) | Pre-registration — the NODE CUT rule (NC): split a known chimeric record at its parent boundary | 2026-10-07 | — |
| [`PREREG_locus_assembly_2026-09-18.md`](archive/2026-09/PREREG_locus_assembly_2026-09-18.md) | Pre-registration — LOCUS ASSEMBLY: reconstruct full-length NPIP copies from reads, O2 for the ambiguous ones | 2026-10-07 | — |
| [`PREREG_chimera_policy_2026-09-18.md`](archive/2026-09/PREREG_chimera_policy_2026-09-18.md) | Pre-registration — `--chimera-policy`: a truth-side lever for known readthrough/chimeric records | 2026-10-07 | — |
| [`PREREG_assembler_widening_2026-09-18.md`](archive/2026-09/PREREG_assembler_widening_2026-09-18.md) | Pre-registration — read-isoform widening in the assembler, single-exon strand, all-alignment pool | 2026-10-07 | — |
| [`PENDING_2026-09-17.md`](archive/2026-09/PENDING_2026-09-17.md) | Pending work — resume list (written 2026-09-17, end of session) | 2026-10-07 | — |
| [`PREREG_junction_fuzz_2026-09-15.md`](archive/2026-09/PREREG_junction_fuzz_2026-09-15.md) | PREREG — what junction-jitter tolerance should RUSTLE_JUNCTION_FUZZ_BP use? (2026-09-15, written before any gffcompare/S | 2026-10-07 | — |
| [`PREREG_excision_extension_2026-09-15.md`](archive/2026-09/PREREG_excision_extension_2026-09-15.md) | PREREG — EXCISION EXTENSION: does the "missing copy" signature generalize beyond NPIP/ZNF875? (2026-09-15) | 2026-10-07 | — |
| [`PREREG_core_definition_2026-09-12.md`](archive/2026-09/PREREG_core_definition_2026-09-12.md) | PREREG — which contiguous-core definition should `confirm_edge` use? (2026-09-12, written before any score below is comp | 2026-10-07 | — |
| [`PREREG_locus_growth_extent_2026-09-10.md`](archive/2026-09/PREREG_locus_growth_extent_2026-09-10.md) | PREREG: junction-aware greedy locus growth for de novo (reads-only) O1 boundary inference | 2026-10-07 | — |
| [`PREREG_identity_weighted_density_2026-09-10.md`](archive/2026-09/PREREG_identity_weighted_density_2026-09-10.md) | PREREG — validating identity-weighted density (advisor criticism 5) genome-wide, priced on a held-out family (2026-09-10 | 2026-10-07 | — |
| [`PREREG_boundary_outlier_2026-09-10.md`](archive/2026-09/PREREG_boundary_outlier_2026-09-10.md) | PREREG — `--min-boundary-fraction`: flag a locus's most-distal boundary when it sits alone, far past the rest (2026-09-1 | 2026-10-07 | — |
| [`PREREG_tie_partner_candidates_2026-09-09.md`](archive/2026-09/PREREG_tie_partner_candidates_2026-09-09.md) | PREREG — candidates are ALL tied placements (user's rule, 2026-09-09; written before the code change) | 2026-10-07 | — |
| [`PREREG_origin_drop_indels_2026-09-09.md`](archive/2026-09/PREREG_origin_drop_indels_2026-09-09.md) | PREREG — `--origin-drop-indels`: exclude indels from the origin certificate's edit count (2026-09-09) | 2026-10-07 | — |
| [`PREREG_modality_2026-09-09.md`](archive/2026-09/PREREG_modality_2026-09-09.md) | PREREG — a whole-distribution test of subfamily structure (2026-09-09, before any run) | 2026-10-07 | — |
| [`PREREG_known_subclusters_2026-09-09.md`](archive/2026-09/PREREG_known_subclusters_2026-09-09.md) | PREREG — do our families recover the KNOWN subcluster structure of NPIP and TBC1D3? (advisor, 2026-09-09) | 2026-10-07 | — |
| [`PREREG_isoform_pool_2026-09-09.md`](archive/2026-09/PREREG_isoform_pool_2026-09-09.md) | PREREG — the ISOFORM-level certificate: pooling a chain's reads before the pairwise test (2026-09-09, before the prototy | 2026-10-07 | — |
| [`PREREG_isoform_phantoms_2026-09-09.md`](archive/2026-09/PREREG_isoform_phantoms_2026-09-09.md) | PREREG — isoform-level phantoms in the shipped GTF: the same isoform emitted at several copies on abstaining reads (2026 | 2026-10-07 | — |
| [`PREREG_isoform_fraction_2026-09-09.md`](archive/2026-09/PREREG_isoform_fraction_2026-09-09.md) | PREREG — `--min-isoform-fraction`: demote low-relative-depth transcripts, StringTie-style (2026-09-09, before the code c | 2026-10-07 | — |
| [`PREREG_isoform_bakeoff_2026-09-09.md`](archive/2026-09/PREREG_isoform_bakeoff_2026-09-09.md) | PREREG — the ISOFORM-level (lift-aware) hard-locus bakeoff (2026-09-09, before the metric runs) | 2026-10-07 | — |
| [`PREREG_indel_psv_2026-09-09.md`](archive/2026-09/PREREG_indel_psv_2026-09-09.md) | PREREG — indel PSV columns in read-star (`--indel-psv`, 2026-09-09, before the code change and before any run) | 2026-10-07 | — |
| [`PREREG_holdout_mcl58_2026-09-09.md`](archive/2026-09/PREREG_holdout_mcl58_2026-09-09.md) | PREREG — the held-back family for the O2 refactor: gorilla MCL58 (2026-09-09, before any run) | 2026-10-07 | — |
| [`PREREG_hard_locus_bakeoff_2026-09-09.md`](archive/2026-09/PREREG_hard_locus_bakeoff_2026-09-09.md) | PREREG — the HARD-LOCUS bakeoff (advisor's metric, Part 0h; 2026-09-09, before any restricted number is computed) | 2026-10-07 | — |
| [`PREREG_gtf_copy_set_2026-09-09.md`](archive/2026-09/PREREG_gtf_copy_set_2026-09-09.md) | PREREG — emit the GTF O2 believes: evidence-placed transcripts and copy SETS (2026-09-09, before the post-processor runs | 2026-10-07 | — |
| [`PREREG_exon_sum_width_bias_2026-09-09.md`](archive/2026-09/PREREG_exon_sum_width_bias_2026-09-09.md) | PREREG — the de novo exon-sum width deficit: how much, which end, and does it correct? (2026-09-09, before the measureme | 2026-10-07 | — |
| [`PREREG_excision_sweep_2026-09-09.md`](archive/2026-09/PREREG_excision_sweep_2026-09-09.md) | PREREG — excision sweep over every human MCL0 copy (2026-09-09, before any run) | 2026-10-07 | — |
| [`PREREG_excision_odi_2026-09-09.md`](archive/2026-09/PREREG_excision_odi_2026-09-09.md) | PREREG — excision control on the origin-drop-indels assignments (2026-09-09, before the runs) | 2026-10-07 | — |
| [`PREREG_best_by_psv_score_2026-09-09.md`](archive/2026-09/PREREG_best_by_psv_score_2026-09-09.md) | PREREG — the read-star "best" must be scored the way the verdict is (2026-09-09, before the code change) | 2026-10-07 | — |
| [`PREREG_best_by_duel_2026-09-09.md`](archive/2026-09/PREREG_best_by_duel_2026-09-09.md) | PREREG — `--best-by-duel`: the read-star best chosen by PAIRWISE duels, not the column-count score (2026-09-09, before t | 2026-10-07 | — |
| [`PREREG_as_tied_only_2026-09-09.md`](archive/2026-09/PREREG_as_tied_only_2026-09-09.md) | PREREG — O2 admits AS-TIED MULTIMAPPERS ONLY (user correction, 2026-09-09) | 2026-10-07 | — |
| [`PREREG_aligner_disagreement_2026-09-09.md`](archive/2026-09/PREREG_aligner_disagreement_2026-09-09.md) | PREREG — `--admit-aligner-disagreement`: widen the gate to aligner self-disagreement (2026-09-09, before any run) | 2026-10-07 | — |
| [`O1_O2_OPEN_2026-09-09.md`](archive/2026-09/O1_O2_OPEN_2026-09-09.md) | O1 / O2 — the open queue as of 2026-09-09 | 2026-10-07 | — |
| [`PSEUDOCODE_2026-09-08.md`](archive/2026-09/PSEUDOCODE_2026-09-08.md) | How the pipeline goes from reads to gene families | 2026-10-07 | — |
| [`PREREG_tool_bakeoff_2026-09-08.md`](archive/2026-09/PREREG_tool_bakeoff_2026-09-08.md) | PREREG — three-tool bakeoff on one multi-copy family, scored by ONE derived rule (2026-09-08) | 2026-10-07 | — |
| [`PREREG_cross_family_overlap_default_2026-09-08.md`](archive/2026-09/PREREG_cross_family_overlap_default_2026-09-08.md) | PREREG — flipping the cross-family exon-overlap rule ON by default (2026-09-08) | 2026-10-07 | — |
| [`PREREG_yags_daz_2026-09-07.md`](archive/2026-09/PREREG_yags_daz_2026-09-07.md) | PREREG — Y ampliconic gene families, DAZ first (2026-09-07) | 2026-10-07 | — |
| [`PREREG_size_reference_2026-09-07.md`](archive/2026-09/PREREG_size_reference_2026-09-07.md) | PREREG — the SIZE reference for bipartite matching, and the exonic-core bias test (2026-09-07) | 2026-10-07 | — |
| [`PREREG_readthrough_object_2026-09-07.md`](archive/2026-09/PREREG_readthrough_object_2026-09-07.md) | PREREG — the conjoined read-through as an O1 object (2026-09-07) | 2026-10-07 | — |
| [`PREREG_npip_human_2026-09-07.md`](archive/2026-09/PREREG_npip_human_2026-09-07.md) | PREREG — the SD-core definition on HUMAN NPIP (2026-09-07) | 2026-10-07 | — |
| [`PREREG_minimap_readthrough_arms_2026-09-07.md`](archive/2026-09/PREREG_minimap_readthrough_arms_2026-09-07.md) | PREREG — can a minimap2 setting remove the read-throughs without destroying the real ones? (2026-09-07) | 2026-10-07 | — |
| [`PREREG_crossspecies_junctions_2026-09-07.md`](archive/2026-09/PREREG_crossspecies_junctions_2026-09-07.md) | PREREG — cross-species replication of the read-through junctions (2026-09-07) | 2026-10-07 | — |
| [`PREREG_core_projection_2026-09-07.md`](archive/2026-09/PREREG_core_projection_2026-09-07.md) | PREREG — recovering unannotated family members by projecting the family's core (2026-09-07) | 2026-10-07 | — |
| [`PREREG_closure_definition_2026-09-07.md`](archive/2026-09/PREREG_closure_definition_2026-09-07.md) | PREREG — the family as a CLOSURE (fixed point of "shares a core with half the family"), 2026-09-07 | 2026-10-07 | — |
| [`PREREG_annotation_ablation_2026-09-07.md`](archive/2026-09/PREREG_annotation_ablation_2026-09-07.md) | PREREG — the annotation-degradation ablation for O1 (2026-09-07) | 2026-10-07 | — |
| [`PREREG_o3_reconstruction_2026-09-06.md`](archive/2026-09/PREREG_o3_reconstruction_2026-09-06.md) | PREREG — reconstruction over every flagged pair (O3 item 2), 2026-09-06 | 2026-10-07 | — |
| [`PREREG_o3_flag_pass_2026-09-06.md`](archive/2026-09/PREREG_o3_flag_pass_2026-09-06.md) | PREREG — the O3 flag pass: every family's rejected reads through the missing-copy detector and the admission prototype ( | 2026-10-07 | — |
| [`PREREG_npip_known_origin_2026-09-06.md`](archive/2026-09/PREREG_npip_known_origin_2026-09-06.md) | PREREG — one sensitivity / one specificity for NPIP on reads of known origin (2026-09-06) | 2026-10-07 | — |
| [`PREREG_npip_all_members_2026-09-06.md`](archive/2026-09/PREREG_npip_all_members_2026-09-06.md) | PREREG — every NPIP member accepted: two definition polishes and the read-through certificate (2026-09-06) | 2026-10-07 | — |
| [`PREREG_worst2_2026-09-05.md`](archive/2026-09/PREREG_worst2_2026-09-05.md) | PREREG — adjudication of the two worst sweep-v2 families (2026-09-05) | 2026-10-07 | — |
| [`PREREG_soto_two_definitions_2026-09-05.md`](archive/2026-09/PREREG_soto_two_definitions_2026-09-05.md) | PREREG — the two O1 definitions on ONE truth: Soto 2025 families, HUMAN A119b/CHM13 v2.0 (2026-09-05) | 2026-10-07 | — |
| [`PREREG_roster_admission_2026-09-05.md`](archive/2026-09/PREREG_roster_admission_2026-09-05.md) | PREREG — roster admission by read origin (O2 → O1 feedback, §6fc), 2026-09-05 | 2026-10-07 | — |
| [`PREREG_q9_npipAB_2026-09-05.md`](archive/2026-09/PREREG_q9_npipAB_2026-09-05.md) | PREREG — Q9: are NPIPA and NPIPB distinct subfamilies in the gorilla catalog? (2026-09-05) | 2026-10-07 | — |
| [`PREREG_loose_ends_L1_L2_2026-09-05.md`](archive/2026-09/PREREG_loose_ends_L1_L2_2026-09-05.md) | PREREG — L1 dropped members as O2 candidates + L2 the locus extent owned by O1 (2026-09-05) | 2026-10-07 | — |
| [`PREREG_loose_end_L3_2026-09-05.md`](archive/2026-09/PREREG_loose_end_L3_2026-09-05.md) | PREREG — L3: single-candidate ties under the pairwise rule (`-p 0`), 2026-09-05 | 2026-10-07 | — |
| [`PREREG_genomic_read_star_2026-09-05.md`](archive/2026-09/PREREG_genomic_read_star_2026-09-05.md) | PREREG — genomic read-star: the read against each candidate's LOCUS, splice-aware (O2-9d, §6fd), 2026-09-05 | 2026-10-07 | — |
| [`PREREG_g3_thresholds_2026-09-05.md`](archive/2026-09/PREREG_g3_thresholds_2026-09-05.md) | PREREG — G3: sensitivity of the catalog to the edge thresholds identity 0.70 / cov_longer 0.30 / min_bp 300 (2026-09-05) | 2026-10-07 | — |
| [`PREREG_fold_within_clusters_2026-09-05.md`](archive/2026-09/PREREG_fold_within_clusters_2026-09-05.md) | PREREG — root-cause fix for MCL "fragmentation" (row 694): fold within clusters + exon-less span (2026-09-05) | 2026-10-07 | — |
| [`PREREG_excision_2026-09-05.md`](archive/2026-09/PREREG_excision_2026-09-05.md) | PREREG — EXCISION: remove a copy, follow its reads, detect "a copy should be here" (O3 in O2's vocabulary), 2026-09-05 | 2026-10-07 | — |
| [`PREREG_d3_molecule_observation_2026-09-05.md`](archive/2026-09/PREREG_d3_molecule_observation_2026-09-05.md) | PREREG — D3 (register 691): a molecule is one observation, scored on its SEQUENCE against every copy (2026-09-05) | 2026-10-07 | — |
| [`FRAMING_AUDIT_2026-09-05.md`](archive/2026-09/FRAMING_AUDIT_2026-09-05.md) | Framing audit — does "gene family, not duplication block" serve the objectives and answer the advisor? (2026-09-05) | 2026-10-07 | — |
| [`PREREG_rna_corroboration_2026-09-04.md`](archive/2026-09/PREREG_rna_corroboration_2026-09-04.md) | PRE-REGISTRATION — RNA corroboration of the MCL catalog (written 2026-09-04, BEFORE the run completed) | 2026-10-07 | — |
