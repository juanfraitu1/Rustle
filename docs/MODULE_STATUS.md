# Module status — what is SHIPPED vs what was TRIED

Generated 2026-09-03 by a 4-agent reachability survey of every `src/rustle/vg_family/` module
(ledger §6dj). Each tag was assigned by **tracing callers up to a binary in `src/bin/`**, never by
reading the module's own header.

⚠⚠ **THE HEADLINE: 29 of 52 module headers DISAGREE with reachability.** A `//!` header is a
CLAIM. Several modules describe themselves in the present tense as live analysis stages while having
**zero production callers**. Verify a header against its callers before quoting it — this project has
also left *retracted findings* in shipped docstrings (the `RUSTLE_LOCUS_EXON_UNION` doc still asserts it
"broke up the component fusing 40 of 83 Soto families", which memory records as *"a 40-family blob that
does not exist in pipeline output; real worst fusion: 2"*).

⭐ **Only 14 of 52 modules are reachable at defaults.** `TEST-ONLY` is not a
criticism — an idea that was built and left switched off is a legitimate outcome. The point is to be able
to tell which is which without re-deriving it.

Enforced by `vg_family::module_status_tests`: every module must carry a `//! **STATUS:** <TAG>` line, and
this file must list exactly the module set.

| tag | count |
|---|---|
| **SHIPPED-DEFAULT** | 14 |
| **OPT-IN** | 12 |
| **OTHER-BINARY** | 11 |
| **REFUTED** | 1 |
| **TEST-ONLY** | 12 |
| **INFRASTRUCTURE** | 2 |

## SHIPPED-DEFAULT (14)

Reachable from a shipped binary with **no env var and no non-default flag**. This is the method.

| module | gate | deciding evidence |
|---|---|---|
| `catalog_input.rs` | - | MEASURED: catalog_input::exon_blocks_str is called at src/bin/gw_family_catalog.rs:317 (inside fn exon_blocks, :313), used at :484 to write the copies.tsv exon-blocks column inside emit_catalog (:320), which main calls uncondition |
| `copy_assign.rs` | - | MEASURED: copy_assign::assign_read is called at copy_assign_pipeline.rs:719 (fn assign_one_read), reached via assign_family_detailed, which detect_and_assign calls unconditionally at denovo_pipeline.rs:2229; detect_and_assign is i |
| `copy_assign_pipeline.rs` | - | MEASURED: assign_family_detailed is called unconditionally at denovo_pipeline.rs:2229 (Stage-1 of detect_and_assign, the function src/bin/copy_assign.rs:1567 drives); the module is imported wholesale at denovo_pipeline.rs:23-25 (b |
| `copy_split.rs` | - | MEASURED: split_locus_copies is called unconditionally at denovo_pipeline.rs:2346 (the collapsed_copies count inside detect_and_assign, no enclosing flag check — contrast the flagged uses at :1393 under recover_collapsed_candidate |
| `denovo_assemble.rs` | - | MEASURED: imported unconditionally by both flagship binaries — src/bin/copy_assign.rs:27-30 (assemble_gate, pass1_skeletons, reads_in_region, BamIndexCache, BamRead, GATE_MIN_READS) and used at src/bin/gw_family_catalog.rs:1006 (r |
| `denovo_pipeline.rs` | - | MEASURED: it IS the driver both flagships call — src/bin/gw_family_catalog.rs:19 imports it and src/bin/copy_assign.rs:34 imports detect_and_assign/catalog_overlaps/DenovoConfig/FamilyAssignment, with detect_and_assign invoked at  |
| `family_detect.rs` | - | MEASURED: denovo_pipeline.rs:3690 calls `collapse_loci_span_aware(&transcripts, &cfg.detect)` in the final unconditional `else` branch of the rep-selection chain (the three earlier branches are the env-gated RUSTLE_LOCUS_EXON_UNIO |
| `family_graph.rs` | - | MEASURED: contiguous_core_coverage_bounded is imported unconditionally at denovo_pipeline.rs:27 and family_detect.rs:32 (the default POA edge criterion); poa_msa_with_costs is called at copy_assign_pipeline.rs:428 and :539; family |
| `family_rescue.rs` | - | MEASURED: denovo_pipeline.rs:2176 calls `rescue_thin_loci_iterative(&loci, &members, &member_spans, genome, &RescueParams::default(), 3)` inside the `for cf in colocated` loop of detect_and_assign. The only guard is `let rescued = |
| `family_split.rs` | - | MEASURED: denovo_pipeline.rs:3348 calls `family_split::gamma_quasi_clique_partition(reps.len(), &edges3, gamma)` unconditionally inside homology_blocks — the RUSTLE_ER_WEIGHTED_PARTITION read at :3335 only chooses whether edge wei |
| `mosaic.rs` | - | MEASURED: `detect_mosaic` is called unconditionally at src/rustle/vg_family/copy_assign_pipeline.rs:1541 inside `assign_family_detailed_once` (copy_assign_pipeline.rs:1421), which is reached from `assign_family_detailed` at src/ru |
| `read_conflict.rs` | - (no gate; opt-OUT only via cfg.homology_primary. Tuning env RUSTLE_CONFLICT_SIG / RUSTLE_CONFLICT_MIN_READS  | MEASURED: locus_unique_mapper_counts runs unconditionally inside detect_and_assign at src/rustle/vg_family/denovo_pipeline.rs:2000, and conflict_edges/conflict_families is the DEFAULT membership oracle at denovo_pipeline.rs:2004-2 |
| `readonly_copy_number.rs` | - for the chi_h leg. The depth_cn leg is gated by `--lambda-global` (src/bin/copy_assign.rs:322 doc, consumed  | MEASURED: src/bin/copy_assign.rs:1983 calls chi_h_with_junctions in the unconditional famcn_rows.push loop; src/bin/copy_assign.rs:1354 comments the table as "always emitted" and copy_assign.rs:2071 writes <out>.famcn_readonly.tsv |
| `rescue_pipeline.rs` | - (no flag, no env var). Suppressed only when copy_assign is run with --families (src/bin/copy_assign.rs:452,  | MEASURED: thin_loci at src/rustle/vg_family/denovo_pipeline.rs:2175 and rescue_thin_loci_iterative at :2176, inside detect_and_assign's per-family `for cf in colocated` loop (prod; test mod starts at denovo_pipeline.rs:7362); dete |

## OPT-IN (12)

Built and wired, but behind a flag that **defaults off**. An arm, not the method — always name the flag when reporting a result from one.

| module | gate | deciding evidence |
|---|---|---|
| `absent_copy.rs` | --absent-copies (src/bin/copy_assign.rs:255-256, default_value_t = false); second route --vg-realign (src/bin/ | MEASURED: the only production call to absent_copy::admit_candidate is src/rustle/vg_family/denovo_pipeline.rs:2252, guarded by `if absent_copies {` at denovo_pipeline.rs:2233; the second entry admit_novel_pools (denovo_pipeline.rs |
| `collapse_enumerate.rs` | --collapse-enumerate (src/bin/gw_family_catalog.rs:177-178, default_value_t = false) or env RUSTLE_COLLAPSE_EN | MEASURED: `if cfg.collapse_enumerate {` at denovo_pipeline.rs:3852 guards the readmit_locus call at :3853; the enclosing branch at :3841 requires collapse_enumerate // collapse_expressed // dna_family_fallback, and DenovoConfig::d |
| `copy_graph.rs` | --phase (src/bin/copy_assign.rs:235-236, default_value_t = false) | MEASURED: the only production construction sites are build_copy_graph at src/bin/copy_assign.rs:1958 and build_exon_graph at :1973, both inside the `if args.phase {` block opened at src/bin/copy_assign.rs:1902; the <out>.exon.gfa  |
| `em_copy_assign.rs` | `--em` (src/bin/copy_assign.rs:310-311, `#[arg(long, default_value_t = false)]`) OR `--vg-realign` (src/bin/co | MEASURED: em_assign_family has exactly two production call sites — src/bin/copy_assign.rs:1853, wrapped in `if args.em {` at copy_assign.rs:1852; and denovo_pipeline.rs:1586 via recompute_realign_abundance, whose sole production c |
| `from_genome.rs` | `--from-genome <BED>` (src/bin/gw_family_catalog.rs:38-39, `#[arg(long)] from_genome: Option<String>`, default | MEASURED: the sole call site of genome_reps/GenomeRepParams in all of src/ is src/bin/gw_family_catalog.rs:628, inside `if let Some(win_bed) = args.from_genome.as_deref() {` at gw_family_catalog.rs:627. The three other grep hits ( |
| `genome_projection.rs` | `--enumerate-copies` (src/bin/gw_family_catalog.rs:172-173, default false) or `--min-identity 0.98` (gw_family | MEASURED: in gw_family_catalog the projection is fenced by `let enumerate = (args.enumerate_copies // args.min_identity == Some(0.98)) && o1_homology;` at gw_family_catalog.rs:930, with the two project_families_batch calls at :949 |
| `hidden_copy.rs` | `--collapse-enumerate` (src/bin/gw_family_catalog.rs:177-178, `#[arg(long, default_value_t = false)]`) | MEASURED: detect_hidden_copy's only production caller in src/ is collapse_enumerate.rs:97 inside `pub fn readmit_locus` (collapse_enumerate.rs:91); readmit_locus's only production caller is denovo_pipeline.rs:3853, inside `if cfg. |
| `linearize.rs` | --linearize / --linearize-gate, both `#[arg(long, default_value_t = false)]` at src/bin/copy_assign.rs:280-281 | MEASURED: the only production call, `super::linearize::linearize_certificate(...)` at src/rustle/vg_family/denovo_pipeline.rs:1838, is reached only through `linearize_cert_if_enabled` (denovo_pipeline.rs:1846) inside the `if do_li |
| `project_all.rs` | --project-all-families, `#[arg(long, default_value_t = false)]` at src/bin/gw_family_catalog.rs:195-196; equiv | MEASURED: the module's only importer is `use rustle::vg_family::project_all::{CopyIn, all_copy_consensuses, known_from_fams, dedup_overlapping, overlaps_any, format_allproj_row};` at src/bin/gw_family_catalog.rs:991, which sits in |
| `seed_projection.rs` | --seed (src/bin/gw_family_catalog.rs:223-224, `#[arg(long)] seed: Vec<String>`, default = empty vec; project_s | MEASURED: the module's only importer is the `use rustle::vg_family::seed_projection::{...}` inside fn project_seeds at src/bin/gw_family_catalog.rs:238 (prod; that file's test mod starts at :1024); project_seeds is called at gw_fa |
| `single_copy.rs` | --single-copy-baseline (src/bin/gw_family_catalog.rs:189-190, `#[arg(long, default_value_t = false)] single_co | MEASURED: single_copy_loci's only production caller is src/rustle/vg_family/denovo_pipeline.rs:2704 inside detect_single_copy_baseline_genome_wide, and that function's only caller in the whole tree is src/bin/gw_family_catalog.rs: |
| `vg_realign.rs` | --vg-realign or --vg-realign-correct (src/bin/copy_assign.rs:340-341 and :346-347, both `default_value_t = fal | MEASURED: the correction leg is guarded by `if cfg.vg_realign` at src/rustle/vg_family/denovo_pipeline.rs:2356 (call at :2376), and DenovoConfig's default is `vg_realign: false` / `vg_realign_admit: false` at denovo_pipeline.rs:14 |

## OTHER-BINARY (11)

Live, but only from a binary other than `gw_family_catalog` / `copy_assign`.

| module | gate | deciding evidence |
|---|---|---|
| `allele_specific_junctions.rs` | - | MEASURED: src/bin/asj.rs:10-12 imports bh_fdr/scan_gene_asj/scan_gene_asj_multisnp/scan_gene_copy_specific_junctions/AsjParams and src/bin/asj_verify.rs:18 imports anchor_junction_dist/is_transversion; no reference from gw_family_ |
| `annotation_families.rs` | - | MEASURED: sole caller is src/bin/mcl_families.rs:18-20 (build_clusters, graph_from_paf, mcl, Cluster, GeneKey, GraphParams). The one other hit, family_detect.rs:207, is a doc comment (`/// Iterative path-halving union-find (matche |
| `bridge_detector.rs` | - | MEASURED: production callers are driver.rs:39/:397/:402 (load_skeletons, load_strand, ExonFetcher), recombinant_split.rs:47 and multi_repeat_bridge.rs:55 — and all three roll up to driver::build_catalog, whose only caller is src/b |
| `driver.rs` | - | MEASURED: driver::build_catalog is called at src/bin/family_define.rs:144 and driver::write_outputs at :186 (module imported at family_define.rs:18); no other src/ file references `driver::`. Binary: family_define. |
| `edge_oracle.rs` | - | MEASURED: the only production consumers of load_pair_core_str/load_universal_aln/load_allele/demote_gene are driver.rs:302-304 and driver.rs:468, and driver.rs's only caller in all of src/bin/ is src/bin/family_define.rs:149 (`dri |
| `family_definition.rs` | - | MEASURED: distinct_loci's only production call is driver.rs:231; refine_families' only production call is driver.rs:386 — both reachable only from src/bin/family_define.rs:149. refine_component/induced_density have one further pro |
| `family_loaders.rs` | - | MEASURED: every production call is in driver.rs — load_meta :290, load_annot :291, load_raw_families :292, load_edges_str :294, build_genes_dict :296, components_from_edges :371 — and driver.rs is reached only from src/bin/family_ |
| `multi_repeat_bridge.rs` | family_define; default-ON there (driver.rs:143 `repeat_bridge_gate: true`), opt-OUT via --no-repeat-bridge-gat | MEASURED: the only production callers are src/rustle/vg_family/driver.rs:434 (`load_node_mult`) and driver.rs:444 (`split_families_repeat_bridge`), inside `if opt.repeat_bridge_gate` (driver.rs:432); `driver::build_catalog` is inv |
| `parcn.rs` | binary `parcn` (Cargo.toml:116-118, path src/bin/parcn.rs); no flag inside it gates the module — the whole bin | MEASURED: the only importer is src/bin/parcn.rs:15-18 (`use rustle::vg_family::parcn::{assign_locus, dedup_loci, format_family_row, format_parcn_row, parse_copies_fa, sun_positions, tabulate, Assignment, CopySun, Locus};`); no oth |
| `recombinant_split.rs` | family_define only; DEFAULT-ON there, opt out with --no-split-recombinants or RUSTLE_NO_SPLIT_RECOMBINANTS=1 ( | MEASURED: sole production caller is `recombinant_split::split_block` at src/rustle/vg_family/driver.rs:274 (driver.rs test mod starts at :531), inside split_families_recombinant which runs under `if opt.split_recombinants` at driv |
| `repeat_catalog.rs` | family_define only. `load_skeletons` runs under `if opt.repeat_bridge_gate` (driver.rs:432; DEFAULT-ON, opt ou | MEASURED: the only production entry points are `repeat_catalog::load_skeletons` at src/rustle/vg_family/driver.rs:433 and `repeat_catalog::dn_exons` at src/rustle/vg_family/multi_repeat_bridge.rs:61 — both on the bin/family_define |

## REFUTED (1)

Implemented, **measured**, and the measurement went against it. Kept deliberately as a re-runnable instrument — a negative result you can reproduce is worth more than one taken on trust.

| module | gate | deciding evidence |
|---|---|---|
| `collapse_gate.rs` | --collapse-gate (src/bin/copy_assign.rs:393-394, default_value_t = false); DenovoConfig::default sets collapse | MEASURED (the refutation is a recorded measurement in-tree): collapse_gate.rs:17-21 — 'DEFAULT OFF. The instrument is not what this module's name claims, and a control proved it… Run genome-wide, the gate fires on EEF1A1 … and rep |

## TEST-ONLY (12)

**No non-test callers anywhere in `src/`.** Dead in every shipped binary. Not deleted, but nothing it claims is in effect.

| module | gate | deciding evidence |
|---|---|---|
| `asj_genetic_core.rs` | - | MEASURED: grep for `genetic_core`/`GeneticCore`/`parse_verified`/`parse_strandbias` across src/ and tests/ returns nothing outside asj_genetic_core.rs except the mod.rs:36 declaration; its own tests start at asj_genetic_core.rs:47 |
| `asj_strand_bias.rs` | - | MEASURED: the only src/ references outside its own file are doc comments — asj_verify.rs:5 and :256 (`///`) and asj_genetic_core.rs:11/:58 (`//!`, `///`). Zero call sites for StrandBiasEngine/sor/strand_table/fisher_strand_p; own  |
| `asj_verify.rs` | - | MEASURED: zero callers of AsjVerifyEngine/parse_calls/verified_line anywhere in src/ or tests/ outside the module; own tests at asj_verify.rs:369. Critically, the like-named binary src/bin/asj_verify.rs does NOT use it — it import |
| `consensus.rs` | - | MEASURED: grep for `consensus::` / consensus_vetoes / map_junctions_to_edges / family_consensus_vetoes over all of src/ and tests/ returns ZERO hits outside consensus.rs itself. Its three pub fns are consensus.rs:36, :67, :105; th |
| `diagnostic.rs` | - | MEASURED: classify_internal, classify_external and cigar_has_long_indel have ZERO production callers — the only src/ references are the re-export at mod.rs:66 and the field TYPE at types.rs:109 (`pub rescue_class: Option<crate::vg |
| `multi_repeat_bridge_tests.rs` | - | MEASURED: the file is pulled in only by `#[cfg(test)] #[path = "multi_repeat_bridge_tests.rs"] mod tests;` at src/rustle/vg_family/multi_repeat_bridge.rs:640-642 — it is not declared in mod.rs and has no other reference in the tre |
| `o2_columns.rs` | - | MEASURED: the sole non-test consumer of `column_alleles` (o2_columns.rs:83) is `use super::o2_columns::column_alleles;` at src/rustle/vg_family/o2_materialize.rs:44, and o2_materialize itself is imported by no binary (o2_materiali |
| `o2_margin_gate.rs` | - | MEASURED: the only non-test import of `assign_read_margin` (o2_margin_gate.rs:75) is `use super::o2_margin_gate::{assign_read_margin, BTOL, ERR, JW, MARGIN};` at src/rustle/vg_family/o2_materialize.rs:45 (used at o2_materialize.rs |
| `o2_materialize.rs` | - | MEASURED: `grep -rn o2_materialize src/` yields, outside the file itself, only mod.rs:42 and four DOC-COMMENT mentions in src/bin/copy_assign.rs:129, :130, :131, :1414 (the `--read-cap` help text) plus the runtime warning string a |
| `positional.rs` | - | MEASURED: `scan_genome_for_all_families` (positional.rs:214) and `scan_genome_for_family_loci` (positional.rs:370) have ZERO callers in src/ or tests/ — `grep -rn 'positional::/scan_genome_for' src/ tests/` returns exactly one lin |
| `recombinant_abstain.rs` | RUSTLE_NO_RECOMBINANT_ABSTAIN (recombinant_abstain.rs:72, read at :81, documented DEFAULT-ON opt-out) — VACUOU | MEASURED: the module's only non-test caller is apply_abstain_to_vg at src/rustle/vg_family/o2_materialize.rs:866 (prod; that file's test mod starts at :1133) — but o2_materialize is imported by ZERO binaries. grep for "o2_material |
| `segdup.rs` | RUSTLE_VG_SEGDUP_WINDOW / _MIN_ID / _MIN_FLANK / _MIN_EACH / _BAND (segdup.rs:109-117, inside SegdupParams::fr | MEASURED: grep for flank_homology_extent, flank_homology_extent_banded, SegdupParams, SegdupExtent and call_segdup_extent across src/ and tests/ returns ZERO hits outside src/rustle/vg_family/segdup.rs itself, whose only remaining |

## INFRASTRUCTURE (2)

Shared utility with no independent objective claim.

| module | gate | deciding evidence |
|---|---|---|
| `minimizers.rs` | - | MEASURED: no src/bin/*.rs file references `minimizers`; its only three consumers are libraries — src/rustle/vg_family/repeat_catalog.rs:39, src/rustle/vg_family/multi_repeat_bridge.rs:60, src/rustle/vg_family/vg_realign.rs:27 — an |
| `seq_utils.rs` | - | MEASURED: `pub(crate) fn reverse_complement` (seq_utils.rs:9) is the file's only item and is used on the default path — denovo_assemble.rs:22 (import) with production uses at denovo_assemble.rs:1414 and :1547 (test mod starts at : |

## ⚠ Header / reachability mismatches (29)

Each describes itself as doing something its callers do not support.

| module | tag | the mismatch |
|---|---|---|
| `minimizers.rs` | INFRASTRUCTURE | INFERRED: minimizers.rs:3-9 calls itself "the FOUNDATION + byte-parity crux of the O1 over-merge-gate migration" and says "Every gate that rests on the repeat catalog ... consumes the output of this ONE function" — true, but every one of those gates is reachab |
| `absent_copy.rs` | OPT-IN | Header (absent_copy.rs:1-21) documents the five gates as if the module were always live and only mentions that gate 1's floor is env-overridable; it never says the module itself is unreachable at defaults. |
| `copy_graph.rs` | OPT-IN | Header (copy_graph.rs:1-3) presents it as 'Copy-graph objects (v1)' — the vehicle for making 'a reference-absent copy visibly an arm the reference does not take' — with no indication that nothing builds one unless --phase is passed. |
| `em_copy_assign.rs` | OPT-IN | TWO mismatches. (1) The //! header (em_copy_assign.rs:20-22) says 'the coupling to Task 1's ReadEvidence.logl arrives in Task 3' — stale; em_assign_family (em_copy_assign.rs:269) already does that coupling. (2) The docstring at em_copy_assign.rs:258 asserts 'I |
| `hidden_copy.rs` | OPT-IN | YES, at the module-index level: src/rustle/vg_family/mod.rs:4 describes the crate as one that 'drives structural detectors (mosaic, segdup, hidden_copy, positional)', which reads as core machinery. hidden_copy is unreachable in any default run — it is a leaf o |
| `seed_projection.rs` | OPT-IN | - (the header's own thesis — that --seed is a QUERY over the emitted catalog, not a term in the definition — is exactly what the code does: it reads `fams` AFTER emit_catalog and never feeds the node set). |
| `single_copy.rs` | OPT-IN | MEASURED (mild): the header calls this "the λ_global baseline that calibrates depth_cn = E_fam / λ_global", which implies a live pipeline coupling. There is none — the coupling is by FILE: gw_family_catalog writes <out>.lambda_global.tsv (gw_family_catalog.rs: |
| `vg_realign.rs` | OPT-IN | - (the header states "Default OFF => every output byte-identical" at vg_realign.rs:14, which matches the code). Worth noting only that the header's first line, "significance-gated (correct + discover)", describes two legs behind two DIFFERENT defaults-off swit |
| `edge_oracle.rs` | OTHER-BINARY | Soft mismatch worth flagging: the //! header (edge_oracle.rs:1) calls this 'the E_r RNA-homology edge oracle'. The E_r edge oracle that the SHIPPED catalog actually runs is a different implementation (family_detect::confirm_edge / detect_edges, denovo_pipeline |
| `family_definition.rs` | OTHER-BINARY | THE BIG ONE. family_definition.rs:1 announces itself as 'O1 multi-copy family-definition predicate `distinct_loci`' and denovo_pipeline.rs:6 states outright that 'the ≥2-distinct-loci certificate is `family_definition::distinct_loci`'. It is not. The shipped g |
| `repeat_catalog.rs` | OTHER-BINARY | MEASURED, and load-bearing: repeat_catalog.rs:5-10 claims "This is the piece the repeat-bridge GATE consumes at runtime ... this module builds exactly those multiplicities" and repeat_catalog.rs:44 defines M_OP=5 as the gate threshold. The gate does NOT build  |
| `collapse_gate.rs` | REFUTED | The module NAME and mod.rs:9's one-liner ('admit a COLLAPSED single-rep locus as a multi-copy family') both still assert collapse detection; the header body itself retracts that ('detects unresolvable PARALOGY, not collapse'). Name and claim disagree, header a |
| `catalog_input.rs` | SHIPPED-DEFAULT | The header (catalog_input.rs:1-15) is entirely about the O1→O2 FILE CONTRACT (parse_copies_tsv/parse_copies_fa/group_families/to_colocated). That half is OPT-IN: it runs only under `--families` (src/bin/copy_assign.rs:452, Option<String> default None), consume |
| `family_graph.rs` | SHIPPED-DEFAULT | YES — the header describes the DEAD half. family_graph.rs:1-8 defines the module as 'union of per-copy splice graphs ... Nodes are exon-equivalence classes'. That object has no production caller: FamilyGraph (family_graph.rs:48), JunctionEdge (:40), extract_co |
| `mosaic.rs` | SHIPPED-DEFAULT | MEASURED — THE BIGGEST ONE IN THIS SLICE: mosaic.rs:14-15 states "Default-OFF in the pipeline (RUSTLE_VG_MOSAIC_ON)". `grep -rn 'RUSTLE_VG_MOSAIC_ON\/MOSAIC_ON' src/ tests/` returns exactly ONE hit — that docstring line itself. The env var is never read; the d |
| `read_conflict.rs` | SHIPPED-DEFAULT | MEASURED: the header at src/rustle/vg_family/read_conflict.rs:22-23 says "The remaining integration is plumbing per-locus secondary placements (`secondary_index` / `tied_secondary_reads_in_region`) into the detection stage" — i.e. it presents the module as NOT |
| `readonly_copy_number.rs` | SHIPPED-DEFAULT | MEASURED (minor, but it is a severed claim): readonly_copy_number.rs:10 is a dangling fragment — "//!  families e.g. `chi_H=1` on a locus whose true copy number is ~11." — the sentence it belonged to is gone, so the stated lower-bound caveat reads as a floatin |
| `rescue_pipeline.rs` | SHIPPED-DEFAULT | - (header calls it "integration stage 4b", which matches). Scope note only: gw_family_catalog does NOT call detect_and_assign (it imports detect_conflict_catalog_genome_wide* / detect_homology_catalog_genome_wide at gw_family_catalog.rs:19-24), so rescue is de |
| `asj_genetic_core.rs` | TEST-ONLY | Header (asj_genetic_core.rs:1-4) calls this 'the O3 DELIVERABLE' producing the 54-call genetic core, and mod.rs:36 repeats 'O3 ASJ DELIVERABLE' — but no shipped binary can produce that file; only #[cfg(test)] runs the funnel. |
| `asj_strand_bias.rs` | TEST-ONLY | Header (asj_strand_bias.rs:1-10) describes the SOR filter as running 'ON TOP of the shipped asj_calls.tsv', implying a live analysis stage; nothing in src/bin/ invokes it. |
| `asj_verify.rs` | TEST-ONLY | Header (asj_verify.rs:1-16) says it 'writes bench/asj_calls_verified.tsv' as the SECOND analysis layer over the shipped calls; the shipped binary of the same name writes <prefix>.asj_verified.tsv from its own inline code and never touches this module. |
| `consensus.rs` | TEST-ONLY | Header (consensus.rs:1-5) and mod.rs:26 both bill it as 'the one SUBTRACTIVE precision lever' that 'VETOES a low-coverage copy's off-consensus junction' — present tense, as if in the pipeline. No shipped binary can reach it. |
| `diagnostic.rs` | TEST-ONLY | Header (diagnostic.rs:5-10) says classify_internal is 'fast, always runs on every rescued read' and that classify_external runs 'when config.vg_rescue_diagnostic == true'. Neither ever runs: the first has no caller at all, and no code path assigns rescue_class |
| `o2_columns.rs` | TEST-ONLY | INFERRED: o2_columns.rs:5-9 says this "is what feeds `super::o2_margin_gate`" as part of a live pipeline; in Rust it feeds a chain that terminates in an `#[ignore]`d test. |
| `o2_margin_gate.rs` | TEST-ONLY | MEASURED: o2_margin_gate.rs:5-7 describes itself as "the gate that `bench/o2_vg_visualization.py::materialize_family` (the O2 VG-materialization pipeline) actually calls" — true of the PYTHON, but the Rust port has no Rust binary above it. Reading the header a |
| `o2_materialize.rs` | TEST-ONLY | MEASURED: copy_assign advertises a `--read-cap` CLI flag whose help (src/bin/copy_assign.rs:129-131) names `o2_materialize::READ_CAP` / `MaterializeConfig::read_cap`, which makes the module look wired from the flag list; copy_assign.rs:1414-1418 then admits it |
| `positional.rs` | TEST-ONLY | MEASURED: positional.rs:13-15 claims "The resulting candidates are then injected into the family graph as 'ghost' copies in Phase 2 so HMM scoring can target them" — there is no Phase-2 injection; the field meant to carry them (types.rs:937) is written by nobo |
| `recombinant_abstain.rs` | TEST-ONLY | MEASURED, and this is the sharpest one in the slice: recombinant_abstain.rs:18 calls apply_abstain_to_vg "the DEFAULT-ON gate leg", and :30 documents RUSTLE_NO_RECOMBINANT_ABSTAIN as its opt-out. Nothing on any binary's path calls it, so "default-ON" is true o |
| `segdup.rs` | TEST-ONLY | MEASURED: segdup.rs:1-8 says the signal "comes from the GENOME sequence (which the VG already loads for family discovery), anchored at the gene the family discovered" — phrasing that reads as wired into discovery. Nothing in discovery, or anywhere else, calls  |
