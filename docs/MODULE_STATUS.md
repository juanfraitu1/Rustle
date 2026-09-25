# Module status — what is SHIPPED vs what was TRIED

> **2026-09-20 (§6s0) — 11 modules in this table were REMOVED, not just marked.** The survey below
> had already MEASURED them as having zero production callers; they were still compiling. Removed:
> asj_genetic_core.rs, asj_strand_bias.rs, asj_verify.rs (the dropped ASJ objective — a closed island
> whose like-named *binaries* never imported it), consensus.rs, segdup.rs, positional.rs,
> diagnostic.rs, and util's bitset/bitvec/coord/hard_counters. positional and diagnostic were pinned
> only by two *type-only* struct fields (`Bundle.rescue_class`, `RunConfig.vg_candidate_loci`) and one
> re-export — 3 lines, now gone. Also removed: `apply_compat_preset`, `apply_stringtie_exact_overrides`,
> `stringtie_exact()` and 26 dead `vg_*` config fields (all uncalled StringTie/VG-HMM-era residue).
> **Recover any of it from tag `retired-modules-2026-09-20`.** Test suite 931 passed / 0 failed; the
> 46 tests lost were tests OF the removed dead modules. multi_repeat_bridge_tests.rs (retired 2026-09-23 with its module) was deliberately
> KEPT — it covers the live `multi_repeat_bridge` O1 gate.

Generated 2026-09-03 by a 4-agent reachability survey of every `src/rustle/vg_family/` module
(ledger §6dj). Each tag was assigned by **tracing callers up to a binary in `src/bin/`**, never by
reading the module's own header.

⚠⚠ **THE HEADLINE: 29 of 52 module headers DISAGREE with reachability.** A `//!` header is a
CLAIM. Several modules describe themselves in the present tense as live analysis stages while having
**zero production callers**. Verify a header against its callers before quoting it — this project has
also left *retracted findings* in shipped docstrings (the `RUSTLE_LOCUS_EXON_UNION` doc still asserts it
"broke up the component fusing 40 of 83 Soto families", which memory records as *"a 40-family blob that
does not exist in pipeline output; real worst fusion: 2"*).

⭐ **Only 14 of 53 modules are reachable at defaults.** `TEST-ONLY` is not a
criticism — an idea that was built and left switched off is a legitimate outcome. The point is to be able
to tell which is which without re-deriving it.

Enforced by `vg_family::module_status_tests`: every module must carry a `//! **STATUS:** <TAG>` line, and
this file must list exactly the module set.

| tag | count |
|---|---|
| **SHIPPED-DEFAULT** | 15 |
| **OPT-IN** | 15 |
| **OTHER-BINARY** | 3 |
| **REFUTED** | 1 |
| **TEST-ONLY** | 0 |
| **INFRASTRUCTURE** | 1 |

> **2026-09-23 wave 6 (consolidation):** the binaries `asj`, `asj_verify`, `debug_poa`, `bam_null_probe`, `index_bam`,
> `bam_header`, `filter_bam_by_as`, `gamma_refine`, `mcl_refine` and the legacy `family_define` were retired (tag
> `notebook-2026-09-23c`, `~/Desktop/Rustle_attic/2026-09-23c/`), and with them the 12 modules only they reached
> (`driver`, `edge_oracle`, `family_definition`, `family_loaders`, `multi_repeat_bridge(+_tests)`, `o2_columns`,
> `o2_margin_gate`, `o2_materialize`, `recombinant_abstain`, `recombinant_split`, `allele_specific_junctions`; the
> `lgamma` helper `missing_copy_flag_pass` used is inlined there). Counts below were recomputed from the remaining rows.

> **2026-09-24 wave 7 (dead-code cleanup):** the modules minimizers, bridge_detector and repeat_catalog were REMOVED — no
> binary reached them once wave 6 retired driver / multi_repeat_bridge / recombinant_split. Their live survivors moved:
> IndexedFasta to the crate-root genome module (not a vg_family module); revcomp (renamed `revcomp_keep_case`),
> hw_distance and aln_id to seq_utils. The dead FamilyGraph half of family_graph went too (its header now describes the
> live contiguous-core kernel, so its mismatch row is gone). Outside vg_family: the crate-root types module shrank to the
> three hash aliases, the bam module to `open_bam` + `exons_from_cigar`, and util was deleted. Recover any of it from tag
> `notebook-2026-09-24`. Counts below were recomputed from the remaining rows.

## SHIPPED-DEFAULT (15)

Reachable from a shipped binary with **no env var and no non-default flag**. This is the method.

| module | gate | deciding evidence |
|---|---|---|
| `catalog_input.rs` | - | MEASURED: catalog_input::exon_blocks_str is called at src/bin/gw_family_catalog.rs:317 (inside fn exon_blocks, :313), used at :484 to write the copies.tsv exon-blocks column inside emit_catalog (:320), which main calls uncondition |
| `copy_assign.rs` | - | MEASURED 2026-09-24: copy_assign::assign_read_editing is called by copy_assign_pipeline.rs `assign_family_detailed_once`, reached via `assign_family_detailed`, which `detect_and_assign` calls unconditionally. (The earlier citation, `assign_read` via `assign_one_read`, is a test-only driver, now `#[cfg(test)]`.) |
| `copy_assign_pipeline.rs` | - | MEASURED: assign_family_detailed is called unconditionally at denovo_pipeline.rs:2229 (Stage-1 of detect_and_assign, the function src/bin/copy_assign.rs:1567 drives); the module is imported wholesale at denovo_pipeline.rs:23-25 (b |
| `copy_split.rs` | - | MEASURED: split_locus_copies is called unconditionally at denovo_pipeline.rs:2346 (the collapsed_copies count inside detect_and_assign, no enclosing flag check — contrast the flagged uses at :1393 under recover_collapsed_candidate |
| `denovo_assemble.rs` | - | MEASURED: imported unconditionally by both flagship binaries — src/bin/copy_assign.rs:27-30 (assemble_gate, pass1_skeletons, reads_in_region, BamIndexCache, BamRead, GATE_MIN_READS) and used at src/bin/gw_family_catalog.rs:1006 (r |
| `denovo_pipeline.rs` | - | MEASURED: it IS the driver both flagships call — src/bin/gw_family_catalog.rs:19 imports it and src/bin/copy_assign.rs:34 imports detect_and_assign/catalog_overlaps/DenovoConfig/FamilyAssignment, with detect_and_assign invoked at  |
| `family_detect.rs` | - | MEASURED: denovo_pipeline.rs:3690 calls `collapse_loci_span_aware(&transcripts, &cfg.detect)` in the final unconditional `else` branch of the rep-selection chain (the three earlier branches are the env-gated RUSTLE_LOCUS_EXON_UNIO |
| `shared_definition.rs` | OPT-IN | MEASURED: reached only from denovo_pipeline.rs when `RUSTLE_SHARED_DEFINITION` is set (shared_definition::enabled()); unset leaves the catalog byte-identical. Read-isoform widening (ISOFORM_MIN_READS = 5) is ON inside that path by default, off with `RUSTLE_SD_READ_ISOFORM=0` (ledger 6m0). |
| `family_graph.rs` | - | MEASURED: contiguous_core_coverage_bounded is imported unconditionally at denovo_pipeline.rs:27 and family_detect.rs:32 (the default POA edge criterion); poa_msa_with_costs is called at copy_assign_pipeline.rs:428 and :539; family |
| `family_rescue.rs` | - | MEASURED: denovo_pipeline.rs:2176 calls `rescue_thin_loci_iterative(&loci, &members, &member_spans, genome, &RescueParams::default(), 3)` inside the `for cf in colocated` loop of detect_and_assign. The only guard is `let rescued = |
| `family_split.rs` | - | MEASURED: denovo_pipeline.rs:3348 calls `family_split::gamma_quasi_clique_partition(reps.len(), &edges3, gamma)` unconditionally inside homology_blocks — the RUSTLE_ER_WEIGHTED_PARTITION read at :3335 only chooses whether edge wei |
| `mosaic.rs` | - | MEASURED: `detect_mosaic` is called unconditionally at src/rustle/vg_family/copy_assign_pipeline.rs:1541 inside `assign_family_detailed_once` (copy_assign_pipeline.rs:1421), which is reached from `assign_family_detailed` at src/ru |
| `read_conflict.rs` | - (no gate; opt-OUT only via cfg.homology_primary. Tuning env RUSTLE_CONFLICT_SIG / RUSTLE_CONFLICT_MIN_READS  | MEASURED: locus_unique_mapper_counts runs unconditionally inside detect_and_assign at src/rustle/vg_family/denovo_pipeline.rs:2000, and conflict_edges/conflict_families is the DEFAULT membership oracle at denovo_pipeline.rs:2004-2 |
| `readonly_copy_number.rs` | - for the chi_h leg. The depth_cn leg is gated by `--lambda-global` (src/bin/copy_assign.rs:322 doc, consumed  | MEASURED: src/bin/copy_assign.rs:1983 calls chi_h_with_junctions in the unconditional famcn_rows.push loop; src/bin/copy_assign.rs:1354 comments the table as "always emitted" and copy_assign.rs:2071 writes <out>.famcn_readonly.tsv |
| `rescue_pipeline.rs` | - (no flag, no env var). Suppressed only when copy_assign is run with --families (src/bin/copy_assign.rs:452,  | MEASURED: thin_loci at src/rustle/vg_family/denovo_pipeline.rs:2175 and rescue_thin_loci_iterative at :2176, inside detect_and_assign's per-family `for cf in colocated` loop (prod; test mod starts at denovo_pipeline.rs:7362); dete |

## OPT-IN (15)

Built and wired, but behind a flag that **defaults off**. An arm, not the method — always name the flag when reporting a result from one.

| module | gate | deciding evidence |
|---|---|---|
| `run_cache.rs` | RUSTLE_CACHE_DIR (unset = nothing read or written) | MEASURED 2026-09-24: `gw_family_catalog` caches the collapsed representatives (`reps/<key>/`) and the E_r all-vs-all PAF (`paf/<key>/`); a warm run is byte-identical to a cold run and to a run without the cache (gorilla NC_073244.2, human chr16). |
| `absent_copy.rs` | --absent-copies (src/bin/copy_assign.rs:255-256, default_value_t = false); second route --vg-realign (src/bin/ | MEASURED: the only production call to absent_copy::admit_candidate is src/rustle/vg_family/denovo_pipeline.rs:2252, guarded by `if absent_copies {` at denovo_pipeline.rs:2233; the second entry admit_novel_pools (denovo_pipeline.rs |
| `collapse_enumerate.rs` | --collapse-enumerate (src/bin/gw_family_catalog.rs:177-178, default_value_t = false) or env RUSTLE_COLLAPSE_EN | MEASURED: `if cfg.collapse_enumerate {` at denovo_pipeline.rs:3852 guards the readmit_locus call at :3853; the enclosing branch at :3841 requires collapse_enumerate // collapse_expressed // dna_family_fallback, and DenovoConfig::d |
| `copy_discovery.rs` | --discover-copies (src/bin/copy_assign.rs:352-353, `#[arg(long, default_value_t = false)] discover_copies: bool`) | MEASURED: every production reference to `copy_discovery::` sits inside an `if args.discover_copies {` block in src/bin/copy_assign.rs -- `tie_partner_placements` at copy_assign.rs:2917 and `discover_copies_for_family` (copy_assign.rs:1442, which calls `cluster_tie_partners`) at copy_assign.rs:2918, both under the gate opened at copy_assign.rs:2912; the `<out>.discovered_copies.tsv` writer at copy_assign.rs:4503-4521 is gated by the same flag. Unset, the block short-circuits to `Vec::new()` (copy_assign.rs:2919-2921), no file is written, and output is byte-identical (regression: `discover_copies_off_by_default_is_byte_identical`, tests/copy_assign_families.rs). |
| `copy_graph.rs` | --phase (src/bin/copy_assign.rs:235-236, default_value_t = false) | MEASURED: the only production construction sites are build_copy_graph at src/bin/copy_assign.rs:1958 and build_exon_graph at :1973, both inside the `if args.phase {` block opened at src/bin/copy_assign.rs:1902; the <out>.exon.gfa  |
| `em_copy_assign.rs` | `--em` (src/bin/copy_assign.rs:310-311, `#[arg(long, default_value_t = false)]`) OR `--vg-realign` (src/bin/co | MEASURED: em_assign_family has exactly two production call sites — src/bin/copy_assign.rs:1853, wrapped in `if args.em {` at copy_assign.rs:1852; and denovo_pipeline.rs:1586 via recompute_realign_abundance, whose sole production c |
| `from_genome.rs` | `--from-genome <BED>` (src/bin/gw_family_catalog.rs:38-39, `#[arg(long)] from_genome: Option<String>`, default | MEASURED: the sole call site of genome_reps/GenomeRepParams in all of src/ is src/bin/gw_family_catalog.rs:628, inside `if let Some(win_bed) = args.from_genome.as_deref() {` at gw_family_catalog.rs:627. The three other grep hits ( |
| `genome_projection.rs` | `--enumerate-copies` (src/bin/gw_family_catalog.rs:172-173, default false) or `--min-identity 0.98` (gw_family | MEASURED: in gw_family_catalog the projection is fenced by `let enumerate = (args.enumerate_copies // args.min_identity == Some(0.98)) && o1_homology;` at gw_family_catalog.rs:930, with the two project_families_batch calls at :949 |
| `hidden_copy.rs` | `--collapse-enumerate` (src/bin/gw_family_catalog.rs:177-178, `#[arg(long, default_value_t = false)]`) | MEASURED: detect_hidden_copy's only production caller in src/ is collapse_enumerate.rs:97 inside `pub fn readmit_locus` (collapse_enumerate.rs:91); readmit_locus's only production caller is denovo_pipeline.rs:3853, inside `if cfg. |
| `linearize.rs` | --linearize / --linearize-gate, both `#[arg(long, default_value_t = false)]` at src/bin/copy_assign.rs:280-281 | MEASURED: the only production call, `super::linearize::linearize_certificate(...)` at src/rustle/vg_family/denovo_pipeline.rs:1838, is reached only through `linearize_cert_if_enabled` (denovo_pipeline.rs:1846) inside the `if do_li
| `missing_copy_flag_pass.rs` | --flag-missing-copies (src/bin/copy_assign.rs:222, `#[arg(long, default_value_t = false)]`); also requires --families (copy_assign.rs:1400-1402 bails if not) | MEASURED: every production call of `poisson_tail`/`finalize_flags`/`detect_missing_copy_pairs`/`classify_orphan_locus` in `src/bin/copy_assign.rs` sits inside an `if args.flag_missing_copies { ... }` block (copy_assign.rs:1927/1945/1960/2462/4173/4211) — none unconditional. `flag_missing_copies` defaults `false`; unset, the whole block short-circuits to `(Vec::new(), Vec::new())` (copy_assign.rs:2462) and no `o3_*` column or `<out>.missing_copy_loci.tsv` file is ever produced, so the default path is untouched. |
| `project_all.rs` | --project-all-families, `#[arg(long, default_value_t = false)]` at src/bin/gw_family_catalog.rs:195-196; equiv | MEASURED: the module's only importer is `use rustle::vg_family::project_all::{CopyIn, all_copy_consensuses, known_from_fams, dedup_overlapping, overlaps_any, format_allproj_row};` at src/bin/gw_family_catalog.rs:991, which sits in |
| `seed_projection.rs` | --seed (src/bin/gw_family_catalog.rs:223-224, `#[arg(long)] seed: Vec<String>`, default = empty vec; project_s | MEASURED: the module's only importer is the `use rustle::vg_family::seed_projection::{...}` inside fn project_seeds at src/bin/gw_family_catalog.rs:238 (prod; that file's test mod starts at :1024); project_seeds is called at gw_fa |
| `single_copy.rs` | --single-copy-baseline (src/bin/gw_family_catalog.rs:189-190, `#[arg(long, default_value_t = false)] single_co | MEASURED: single_copy_loci's only production caller is src/rustle/vg_family/denovo_pipeline.rs:2704 inside detect_single_copy_baseline_genome_wide, and that function's only caller in the whole tree is src/bin/gw_family_catalog.rs: |
| `vg_realign.rs` | --vg-realign or --vg-realign-correct (src/bin/copy_assign.rs:340-341 and :346-347, both `default_value_t = fal | MEASURED: the correction leg is guarded by `if cfg.vg_realign` at src/rustle/vg_family/denovo_pipeline.rs:2356 (call at :2376), and DenovoConfig's default is `vg_realign: false` / `vg_realign_admit: false` at denovo_pipeline.rs:14 |

## OTHER-BINARY (3)

Live, but only from a binary other than `gw_family_catalog` / `copy_assign`.

| module | gate | deciding evidence |
|---|---|---|
| `annotation_families.rs` | - | MEASURED: sole caller is src/bin/mcl_families.rs:18-20 (build_clusters, graph_from_paf, mcl, Cluster, GeneKey, GraphParams). The one other hit, family_detect.rs:207, is a doc comment (`/// Iterative path-halving union-find (matche |
| `missing_copy.rs` | - | MEASURED: sole caller is src/bin/missing_copy_flag.rs (`use rustle::vg_family::missing_copy::*`, §6ze RNA-only O3 chain: two_means/split_pile/consistency/patched_consensus/home/verdict/confirmed); no reference from gw_family_catalog or copy_assign. |
| `parcn.rs` | binary `parcn` (Cargo.toml:116-118, path src/bin/parcn.rs); no flag inside it gates the module — the whole bin | MEASURED: the only importer is src/bin/parcn.rs:15-18 (`use rustle::vg_family::parcn::{assign_locus, dedup_loci, format_family_row, format_parcn_row, parse_copies_fa, sun_positions, tabulate, Assignment, CopySun, Locus};`); no oth |

## REFUTED (1)

Implemented, **measured**, and the measurement went against it. Kept deliberately as a re-runnable instrument — a negative result you can reproduce is worth more than one taken on trust.

| module | gate | deciding evidence |
|---|---|---|
| `collapse_gate.rs` | --collapse-gate (src/bin/copy_assign.rs:393-394, default_value_t = false); DenovoConfig::default sets collapse | MEASURED (the refutation is a recorded measurement in-tree): collapse_gate.rs:17-21 — 'DEFAULT OFF. The instrument is not what this module's name claims, and a control proved it… Run genome-wide, the gate fires on EEF1A1 … and rep |

## TEST-ONLY (0)

**No non-test callers anywhere in `src/`.** Dead in every shipped binary. Not deleted, but nothing it claims is in effect.

| module | gate | deciding evidence |
|---|---|---|

## INFRASTRUCTURE (1)

Shared utility with no independent objective claim.

| module | gate | deciding evidence |
|---|---|---|
| `seq_utils.rs` | - | MEASURED 2026-09-24: `reverse_complement` is used on the default path (denovo_assemble.rs, family_detect.rs, family_rescue.rs, denovo_pipeline.rs); `revcomp_keep_case` by denovo_pipeline.rs `outside_pseudo_copy` and by vg_realign.rs; `aln_id`/`hw_distance` by vg_realign.rs (opt-in). The last three moved here from the removed `bridge_detector.rs`. |

## ⚠ Header / reachability mismatches (13)

Each describes itself as doing something its callers do not support.

| module | tag | the mismatch |
|---|---|---|
| `absent_copy.rs` | OPT-IN | Header (absent_copy.rs:1-21) documents the five gates as if the module were always live and only mentions that gate 1's floor is env-overridable; it never says the module itself is unreachable at defaults. |
| `copy_graph.rs` | OPT-IN | Header (copy_graph.rs:1-3) presents it as 'Copy-graph objects (v1)' — the vehicle for making 'a reference-absent copy visibly an arm the reference does not take' — with no indication that nothing builds one unless --phase is passed. |
| `em_copy_assign.rs` | OPT-IN | TWO mismatches. (1) The //! header (em_copy_assign.rs:20-22) says 'the coupling to Task 1's ReadEvidence.logl arrives in Task 3' — stale; em_assign_family (em_copy_assign.rs:269) already does that coupling. (2) The docstring at em_copy_assign.rs:258 asserts 'I |
| `hidden_copy.rs` | OPT-IN | YES, at the module-index level: src/rustle/vg_family/mod.rs:4 describes the crate as one that 'drives structural detectors (mosaic, segdup, hidden_copy, positional)', which reads as core machinery. hidden_copy is unreachable in any default run — it is a leaf o |
| `seed_projection.rs` | OPT-IN | - (the header's own thesis — that --seed is a QUERY over the emitted catalog, not a term in the definition — is exactly what the code does: it reads `fams` AFTER emit_catalog and never feeds the node set). |
| `single_copy.rs` | OPT-IN | MEASURED (mild): the header calls this "the λ_global baseline that calibrates depth_cn = E_fam / λ_global", which implies a live pipeline coupling. There is none — the coupling is by FILE: gw_family_catalog writes <out>.lambda_global.tsv (gw_family_catalog.rs: |
| `vg_realign.rs` | OPT-IN | - (the header states "Default OFF => every output byte-identical" at vg_realign.rs:14, which matches the code). Worth noting only that the header's first line, "significance-gated (correct + discover)", describes two legs behind two DIFFERENT defaults-off swit |
| `collapse_gate.rs` | REFUTED | The module NAME and mod.rs:9's one-liner ('admit a COLLAPSED single-rep locus as a multi-copy family') both still assert collapse detection; the header body itself retracts that ('detects unresolvable PARALOGY, not collapse'). Name and claim disagree, header a |
| `catalog_input.rs` | SHIPPED-DEFAULT | The header (catalog_input.rs:1-15) is entirely about the O1→O2 FILE CONTRACT (parse_copies_tsv/parse_copies_fa/group_families/to_colocated). That half is OPT-IN: it runs only under `--families` (src/bin/copy_assign.rs:452, Option<String> default None), consume |
| `mosaic.rs` | SHIPPED-DEFAULT | MEASURED — THE BIGGEST ONE IN THIS SLICE: mosaic.rs:14-15 states "Default-OFF in the pipeline (RUSTLE_VG_MOSAIC_ON)". `grep -rn 'RUSTLE_VG_MOSAIC_ON\/MOSAIC_ON' src/ tests/` returns exactly ONE hit — that docstring line itself. The env var is never read; the d |
| `read_conflict.rs` | SHIPPED-DEFAULT | MEASURED: the header at src/rustle/vg_family/read_conflict.rs:22-23 says "The remaining integration is plumbing per-locus secondary placements (`secondary_index` / `tied_secondary_reads_in_region`) into the detection stage" — i.e. it presents the module as NOT |
| `readonly_copy_number.rs` | SHIPPED-DEFAULT | MEASURED (minor, but it is a severed claim): readonly_copy_number.rs:10 is a dangling fragment — "//!  families e.g. `chi_H=1` on a locus whose true copy number is ~11." — the sentence it belonged to is gone, so the stated lower-bound caveat reads as a floatin |
| `rescue_pipeline.rs` | SHIPPED-DEFAULT | - (header calls it "integration stage 4b", which matches). Scope note only: gw_family_catalog does NOT call detect_and_assign (it imports detect_conflict_catalog_genome_wide* / detect_homology_catalog_genome_wide at gw_family_catalog.rs:19-24), so rescue is de |
