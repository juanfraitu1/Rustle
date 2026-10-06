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

⭐ **Only 14 of 53 modules are reachable at defaults** (2026-09-03; 16 rows after the 2026-09-29 flip below). `TEST-ONLY` is not a
criticism — an idea that was built and left switched off is a legitimate outcome. The point is to be able
to tell which is which without re-deriving it.

Enforced by `vg_family::module_status_tests`: every module must carry a `//! **STATUS:** <TAG>` line, and
this file must list exactly the module set.

| tag | count |
|---|---|
| **SHIPPED-DEFAULT** | 12 |
| **OPT-IN** | 10 |
| **OTHER-BINARY** | 3 |
| **REFUTED** | 0 |
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

> **2026-09-29 default flips (the user's decision; `bench/ASSEMBLY_POLISH.md` addendum 3, `REPRODUCE.md`):** `bridge_regroup.rs` moved from
> OPT-IN to SHIPPED-DEFAULT (`copy_assign --assemble-only` now runs `--bridge-regroup f1v2`; `--bridge-regroup off` = the old
> products), and `mcl_families --min-cov-shorter` (the §6x4 containment escape inside `annotation_families.rs`, OTHER-BINARY)
> defaults to 0.70 instead of 0 (`--min-cov-shorter 0` = the old edge weights; known regressions: NPIP-guided Soto F .833 -> .800,
> semi-guided precision .973 -> .833, register 1007/1009). Both old behaviours were proven byte-identical on human_testis.

> **2026-09-30 units port (opt-in only, nothing flipped):** `bridge_regroup.rs` gained the arm `--bridge-regroup f1units` and
> `--bridge-units-list` (still SHIPPED-DEFAULT for its default `f1v2`; both new switches are OPT-IN and leave every default
> `bench/ASSEMBLY_POLISH.md` 2026-09-30 addendum, `docs/PREREG_container_units_v2_dev_2026-09-30.md` Part C.

## SHIPPED-DEFAULT (17)

Reachable from a shipped binary with **no env var and no non-default flag**. This is the method.

| module | gate | deciding evidence |
|---|---|---|
| `catalog_input.rs` | - | MEASURED: catalog_input::exon_blocks_str is called at src/bin/gw_family_catalog.rs:317 (inside fn exon_blocks, :313), used at :484 to write the copies.tsv exon-blocks column inside emit_catalog (:320), which main calls uncondition |
| `copy_assign.rs` (absorbed `em_copy_assign` + `copy_assign_pipeline` 2026-10-05) | - | MEASURED 2026-09-24: copy_assign::assign_read_editing is called by copy_assign_pipeline.rs `assign_family_detailed_once`, reached via `assign_family_detailed`, which `detect_and_assign` calls unconditionally. (The earlier citation, `assign_read` via `assign_one_read`, is a test-only driver, now `#[cfg(test)]`.) |
| `copy_split.rs` | - | MEASURED: split_locus_copies is called unconditionally at denovo_pipeline.rs:2346 (the collapsed_copies count inside detect_and_assign, no enclosing flag check — contrast the flagged uses at :1393 under recover_collapsed_candidate |
| `denovo_assemble.rs` | - | MEASURED: imported unconditionally by both flagship binaries — src/bin/copy_assign.rs:27-30 (assemble_gate, pass1_skeletons, reads_in_region, BamIndexCache, BamRead, GATE_MIN_READS) and used at src/bin/gw_family_catalog.rs:1006 (r |
| `bridge_regroup.rs` | - (default `f1v2` under `--assemble-only` since 2026-09-29: src/bin/copy_assign.rs `bridge_regroup: Option<String>`, unset resolves to f1v2 in `resolve_bridge_mode`; `--bridge-regroup off` = the 2026-09-25 products; driver `RUSTLE_BRIDGE_REGROUP`, unset = f1v2) | MEASURED 2026-09-29: the only production calls (`bridge_regroup::run`, the evidence hooks in `stream_pass1_region` and `bridge_evidence_region`) sit behind `bridge_mode.is_some()` / `acc.bridge`, now `Some(F1v2)` on every `--assemble-only` run without `--families` and `None` (no file written, products byte-identical to the previous binary) with `off` or outside `--assemble-only`. A port of the frozen `bench/f1_bridge.py` (37ee8e77) and `f1v2.py` (b4e788ad): GTF, families GTF and side tables byte-identical to their held-out outputs (`docs/PREREG_f1_bridge_locus_2026-09-28.md`, `docs/PREREG_f1v2_readshare_2026-09-29.md`). **Default flipped 2026-09-29 by the user's decision** (was OPT-IN, `default_value = "off"`): F1v2 EFFECTIVE held out on both human libraries, and its families beat BASE on every Compara metric on both human substrates (`docs/PREREG_o1_cover_growth_2026-09-29.md` Outcome, the side result). Proven on human_testis: driver `RUSTLE_BRIDGE_REGROUP=off` = the 3007c3d4 build's default products byte for byte; the new default = the port's `f1v2` products. **2026-09-30: two OPT-IN switches of the same module, both off unless named:** `--bridge-regroup f1units` (driver `RUSTLE_BRIDGE_REGROUP=f1units`: F1's bridge junctions without the read-share rule, each bridge transcript cut into UNITS in `<out>.families.gtf` only, scoped native regroup of the gene_ids that hold a cut; `<out>.gtf` and `bridge_junctions.tsv` are `f1`'s) and `--bridge-units-list FILE` (driver `RUSTLE_BRIDGE_UNITS_LIST`: the cuts named by any detector, F1's evidence not read). The only production calls (`run` with `Mode::F1Units`, `run_list`, `build_units`) sit behind `bridge_mode == Some(F1Units)`, never the default; the default products are byte-identical to the e163d955 build (human_testis, gorilla OR6737). Units, gene_ids and tables byte-identical to the dev prototype `units2.py --regroup scoped` (9fea69a1) on gorilla S f = .5 / 1 and human A119b chr16 (R arm and the annotation-oracle list); the production helpers `native_components` + `best_rep` are property-tested equal to `family_detect::collapse_loci_groups`. Under `f1units` `<out>.gtf` is f1's file: it tags exactly the transcripts that were cut. An empty list (header only) is an error that says to run without `--bridge-units-list`. |
| `denovo_pipeline.rs` | - | MEASURED: it IS the driver both flagships call — src/bin/gw_family_catalog.rs:19 imports it and src/bin/copy_assign.rs:34 imports detect_and_assign/catalog_overlaps/DenovoConfig/FamilyAssignment, with detect_and_assign invoked at  |
| `family_detect.rs` (absorbed `family_graph` + `family_split` 2026-10-05) | - | MEASURED: denovo_pipeline.rs:3690 calls `collapse_loci_span_aware(&transcripts, &cfg.detect)` in the final unconditional `else` branch of the rep-selection chain (the three earlier branches are the env-gated RUSTLE_LOCUS_EXON_UNIO |
| `shared_definition.rs` | OPT-IN | MEASURED: reached only from denovo_pipeline.rs when `RUSTLE_SHARED_DEFINITION` is set (shared_definition::enabled()); unset leaves the catalog byte-identical. Read-isoform widening (ISOFORM_MIN_READS = 5) is ON inside that path by default, off with `RUSTLE_SD_READ_ISOFORM=0` (ledger 6m0). |
| `family_rescue.rs` (absorbed `rescue_pipeline` 2026-10-05) | - | MEASURED: denovo_pipeline.rs:2176 calls `rescue_thin_loci_iterative(&loci, &members, &member_spans, genome, &RescueParams::default(), 3)` inside the `for cf in colocated` loop of detect_and_assign. The only guard is `let rescued = |
| `fam_from_gtf.rs` (absorbed `family_container` + `family_relations` 2026-10-05) | - | 2026-10-04: extracted verbatim from `src/bin/mcl_families.rs` (the `--from-gtf` loci, all-vs-all + PAF cache, and `write_locus_rep_copies` copy table) into the library; `mcl_families` re-imports every moved function unchanged. Byte-identity of its products cmp-checked against the pre-extraction binary (bench/MERGED_PIPELINE.md). |
| `mosaic.rs` | - | MEASURED: `detect_mosaic` is called unconditionally at src/rustle/vg_family/copy_assign_pipeline.rs:1541 inside `assign_family_detailed_once` (copy_assign_pipeline.rs:1421), which is reached from `assign_family_detailed` at src/ru |
| `read_conflict.rs` | - (no gate; opt-OUT only via cfg.homology_primary. Tuning env RUSTLE_CONFLICT_SIG / RUSTLE_CONFLICT_MIN_READS  | MEASURED: locus_unique_mapper_counts runs unconditionally inside detect_and_assign at src/rustle/vg_family/denovo_pipeline.rs:2000, and conflict_edges/conflict_families is the DEFAULT membership oracle at denovo_pipeline.rs:2004-2 |
| `readonly_copy_number.rs` | - for the chi_h leg. The depth_cn leg is gated by `--lambda-global` (src/bin/copy_assign.rs:322 doc, consumed  | MEASURED: src/bin/copy_assign.rs:1983 calls chi_h_with_junctions in the unconditional famcn_rows.push loop; src/bin/copy_assign.rs:1354 comments the table as "always emitted" and copy_assign.rs:2071 writes <out>.famcn_readonly.tsv |

## OPT-IN (17)

Built and wired, but behind a flag that **defaults off**. An arm, not the method — always name the flag when reporting a result from one.

| module | gate | deciding evidence |
|---|---|---|
| `run_cache.rs` | RUSTLE_CACHE_DIR (unset = nothing read or written) | MEASURED 2026-09-24: `gw_family_catalog` caches the collapsed representatives (`reps/<key>/`) and the E_r all-vs-all PAF (`paf/<key>/`); a warm run is byte-identical to a cold run and to a run without the cache (gorilla NC_073244.2, human chr16). 2026-09-28: the families PAF (`mcl_families --from-gtf`, key `rustle families paf v2`) is keyed on a 128-bit hash of the loci FASTA taken while it is written and replayed by hard link from a PINNED entry (mtime/inode/sampled-content pins; `RUSTLE_CACHE_VERIFY=1` full re-hash); no cache / cold / warm byte-identical (gorilla NC_073244.2, human_testis, chimp_PTR). |
| `absent_copy.rs` | --absent-copies (src/bin/copy_assign.rs:255-256, default_value_t = false); second route --vg-realign (src/bin/ | MEASURED: the only production call to absent_copy::admit_candidate is src/rustle/vg_family/denovo_pipeline.rs:2252, guarded by `if absent_copies {` at denovo_pipeline.rs:2233; the second entry admit_novel_pools (denovo_pipeline.rs |
| `collapse_enumerate.rs` (absorbed `hidden_copy` + the REFUTED `collapse_gate` instrument 2026-10-05) | --collapse-enumerate (src/bin/gw_family_catalog.rs:177-178, default_value_t = false) or env RUSTLE_COLLAPSE_EN | MEASURED: `if cfg.collapse_enumerate {` at denovo_pipeline.rs:3852 guards the readmit_locus call at :3853; the enclosing branch at :3841 requires collapse_enumerate // collapse_expressed // dna_family_fallback, and DenovoConfig::d |
| `copy_graph.rs` (absorbed `copy_discovery` 2026-10-05) | --phase (src/bin/copy_assign.rs:235-236, default_value_t = false) | MEASURED: the only production construction sites are build_copy_graph at src/bin/copy_assign.rs:1958 and build_exon_graph at :1973, both inside the `if args.phase {` block opened at src/bin/copy_assign.rs:1902; the <out>.exon.gfa  |
| `from_genome.rs` (absorbed `project_all` + `seed_projection` 2026-10-05) | `--from-genome <BED>` (src/bin/gw_family_catalog.rs:38-39, `#[arg(long)] from_genome: Option<String>`, default | MEASURED: the sole call site of genome_reps/GenomeRepParams in all of src/ is src/bin/gw_family_catalog.rs:628, inside `if let Some(win_bed) = args.from_genome.as_deref() {` at gw_family_catalog.rs:627. The three other grep hits ( |
| `genome_projection.rs` | `--enumerate-copies` (src/bin/gw_family_catalog.rs:172-173, default false) or `--min-identity 0.98` (gw_family | MEASURED: in gw_family_catalog the projection is fenced by `let enumerate = (args.enumerate_copies // args.min_identity == Some(0.98)) && o1_homology;` at gw_family_catalog.rs:930, with the two project_families_batch calls at :949 |
| `linearize.rs` | --linearize / --linearize-gate, both `#[arg(long, default_value_t = false)]` at src/bin/copy_assign.rs:280-281 | MEASURED: the only production call, `super::linearize::linearize_certificate(...)` at src/rustle/vg_family/denovo_pipeline.rs:1838, is reached only through `linearize_cert_if_enabled` (denovo_pipeline.rs:1846) inside the `if do_li
| `single_copy.rs` | --single-copy-baseline (src/bin/gw_family_catalog.rs:189-190, `#[arg(long, default_value_t = false)] single_co | MEASURED: single_copy_loci's only production caller is src/rustle/vg_family/denovo_pipeline.rs:2704 inside detect_single_copy_baseline_genome_wide, and that function's only caller in the whole tree is src/bin/gw_family_catalog.rs: |
| `vg_realign.rs` | --vg-realign or --vg-realign-correct (src/bin/copy_assign.rs:340-341 and :346-347, both `default_value_t = fal | MEASURED: the correction leg is guarded by `if cfg.vg_realign` at src/rustle/vg_family/denovo_pipeline.rs:2356 (call at :2376), and DenovoConfig's default is `vg_realign: false` / `vg_realign_admit: false` at denovo_pipeline.rs:14 |

## OTHER-BINARY (4)

Live, but only from a binary other than `gw_family_catalog` / `copy_assign`.

| module | gate | deciding evidence |
|---|---|---|
| `annotation_families.rs` | - | MEASURED: sole caller is src/bin/mcl_families.rs:18-20 (build_clusters, graph_from_paf, mcl, Cluster, GeneKey, GraphParams). The one other hit, family_detect.rs:207, is a doc comment (`/// Iterative path-halving union-find (matche |
| `missing_copy.rs` (absorbed `missing_copy_flag_pass` 2026-10-05) | - | MEASURED: sole caller is src/bin/missing_copy_flag.rs (`use rustle::vg_family::missing_copy::*`, §6ze RNA-only O3 chain: two_means/split_pile/consistency/patched_consensus/home/verdict/confirmed); no reference from gw_family_catalog or copy_assign. |
| `parcn.rs` | binary `parcn` (Cargo.toml:116-118, path src/bin/parcn.rs); no flag inside it gates the module — the whole bin | MEASURED: the only importer is src/bin/parcn.rs:15-18 (`use rustle::vg_family::parcn::{assign_locus, dedup_loci, format_family_row, format_parcn_row, parse_copies_fa, sun_positions, tabulate, Assignment, CopySun, Locus};`); no oth |
| `o3_candidates.rs` | binary `o3_candidates` (Cargo.toml `[[bin]] o3_candidates`, src/bin/o3_candidates.rs); no flag inside it gates the module — the whole binary is the O3 candidates stage (spec docs/superpowers/specs/2026-10-02-o3-candidates-design.md) | MEASURED 2026-10-03: the only importer is src/bin/o3_candidates.rs (`use rustle::vg_family::o3_candidates::{...}`: the pass-B attribution of prereg Amendment 13b (`attribute_by_hits`, `is_poorly_placed`, `MM2_ATTRIB`; the k-mer index of 2026-10-02 is retired), PAF/cs parsing, cluster_reads, the structural template (`structural_template`, `refined_template`, Amendments 13d/13e), consensus_from_template, refine_cluster, variant_is_real, classify, components, is_flagged, union_sequence_with_note, the cached minimap2 runner and the writers); no reference from gw_family_catalog or copy_assign. The pipeline driver's `candidates` stage (plan task 9) runs that binary; the stage is OPT-IN since 2026-10-02 (`--candidates`; ruling R14: its pre-registered acceptance, Amendment 12, FAILED — docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md, its re-run A13 passed — docs/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md; and Amendment 14's no-deletion control failed — docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md), so no default run reaches the module. Was TEST-ONLY until the binary landed (plan task 8). |

## REFUTED (1)

Implemented, **measured**, and the measurement went against it. Kept deliberately as a re-runnable instrument — a negative result you can reproduce is worth more than one taken on trust.

| module | gate | deciding evidence |
|---|---|---|

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
| `single_copy.rs` | OPT-IN | MEASURED (mild): the header calls this "the λ_global baseline that calibrates depth_cn = E_fam / λ_global", which implies a live pipeline coupling. There is none — the coupling is by FILE: gw_family_catalog writes <out>.lambda_global.tsv (gw_family_catalog.rs: |
| `vg_realign.rs` | OPT-IN | - (the header states "Default OFF => every output byte-identical" at vg_realign.rs:14, which matches the code). Worth noting only that the header's first line, "significance-gated (correct + discover)", describes two legs behind two DIFFERENT defaults-off swit |
| `catalog_input.rs` | SHIPPED-DEFAULT | The header (catalog_input.rs:1-15) is entirely about the O1→O2 FILE CONTRACT (parse_copies_tsv/parse_copies_fa/group_families/to_colocated). That half is OPT-IN: it runs only under `--families` (src/bin/copy_assign.rs:452, Option<String> default None), consume |
| `mosaic.rs` | SHIPPED-DEFAULT | MEASURED — THE BIGGEST ONE IN THIS SLICE: mosaic.rs:14-15 states "Default-OFF in the pipeline (RUSTLE_VG_MOSAIC_ON)". `grep -rn 'RUSTLE_VG_MOSAIC_ON\/MOSAIC_ON' src/ tests/` returns exactly ONE hit — that docstring line itself. The env var is never read; the d |
| `read_conflict.rs` | SHIPPED-DEFAULT | MEASURED: the header at src/rustle/vg_family/read_conflict.rs:22-23 says "The remaining integration is plumbing per-locus secondary placements (`secondary_index` / `tied_secondary_reads_in_region`) into the detection stage" — i.e. it presents the module as NOT |
| `readonly_copy_number.rs` | SHIPPED-DEFAULT | MEASURED (minor, but it is a severed claim): readonly_copy_number.rs:10 is a dangling fragment — "//!  families e.g. `chi_H=1` on a locus whose true copy number is ~11." — the sentence it belonged to is gone, so the stated lower-bound caveat reads as a floatin |
