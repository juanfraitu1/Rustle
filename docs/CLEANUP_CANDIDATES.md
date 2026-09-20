# Cleanup candidates — likely dead / likely superseded files

Generated 2026-09-19 at `dna-from-genome@174ae65a` by `tools/audit_cleanup_candidates.py` (re-run it; this file is overwritten). **Read-only: nothing was moved, edited or deleted.** Full per-file table: `docs/cleanup_candidates.tsv` (filter on `class` and `confidence`).

⚠ A mark is a *candidate*, not a verdict. Before deleting anything: (1) grep the path once more, (2) check the `ledger_sections` / `anchor_citers` columns, (3) prefer `git mv` into an archive directory over `rm` for anything tracked, (4) remember `bench/` data can be slow to regenerate (AGENTS.md §2).

## How a file is kept

A file is **anchored** if a path-like token naming it (full path, unique path suffix, unique basename, output prefix, a glob matching exactly 1 file, an enclosing directory of ≤40 files, or a Python `import` resolved against the importing script's directory and its parents) appears in an anchor source: `docs/o1_ledger.md`, `docs/NEGATIVE_RESULTS_REGISTER.md`, `docs/PREREG_*.md`, top-level `docs/*.md`, `docs/experiments/*.md`, `README.md`, the auto-memory directory, Rust `src/`/`tests/`, or `Cargo.toml`. Anchoring then propagates: a file named by an anchored script or anchored markdown write-up is anchored, index sidecars (`.fai`, `.bai`) follow their parent, and a script that names anchored data (its generator — the naming line looks like a write) is anchored; a script that merely reads anchored data is not. An ambiguous bare basename (e.g. `reads.fa`, found in many directories) never anchors on its own. Citations from `docs/archive/`, `docs/superpowers/`, `AGENTS.md` and un-anchored scripts are *weak*: recorded, not counted.

## Rules (applied in this order; stale cutoff = 2026-08-20)

| class | confidence | rule | files | bytes |
|---|---|---|---:|---:|
| **PROTECTED** | - | build/config files, Rust sources and tests (Rust reachability is docs/MODULE_STATUS.md's job, enforced by module_status_tests), and the anchor docs themselves. Never a candidate. | 158 | 6M |
| **TEMP** | high | Python/pytest caches; untracked or git-ignored files at the repo root; untracked *.log / *err* / *out* / *.patch.txt / checkpoint files that nothing cites. | 73 | 17M |
| **REFUTED-MODULE** | medium | Rust module whose `//! **STATUS:**` header is REFUTED and that no other file names (a REFUTED module that is still imported, e.g. collapse_gate.rs, stays PROTECTED). | 0 | 0B |
| **SUPERSEDED-PORTED** | medium | Python script that a Rust source line declares it ports ('Port of', 'Faithful Rust port of', 'Mirrors', 'migration'); 'low' when only a function or part is ported (`x.py::f`, `x.py loaders`). The Python may still serve as a parity oracle or golden-fixture generator -- check tests before deleting. | 23 | 580K |
| **SUPERSEDED-VERSION** | medium | Older member of a version series in the same directory (_v1.._vN, foo/foo2/foo3, dated _YYYY-MM-DD copies, foo vs foo_fix/_final/_new) that nothing anchors. 'high' when the newest member IS anchored. | 26 | 17M |
| **SUPERSEDED-CITED** | low | Older member of a version series that IS anchored: provenance for a recorded result -- archive, don't delete. | 10 | 633K |
| **LEGACY-ASSEMBLER** | medium | Not anchored, and its path or first 200 lines name StringTie-era assembler machinery (bundle/transfrag/parity/gffcompare/...) -- the assembler layer was retired (docs/RETIREMENT_AND_MIGRATION.md). | 52 | 482K |
| **AMBIGUOUS-CITE** | low | Not anchored; an anchor source names it only by a bare basename shared by several files, a directory too large (>40 files) or a glob matching several files -- may or may not mean this copy. | 617 | 31M |
| **ORPHAN** | medium | Not anchored; named only by files that are themselves not anchored (e.g. a figure named only by its un-cited plotting script). | 258 | 141M |
| **UNCITED-STALE** | medium | Not named by anything (wide globs / big directories / scripts' ambiguous basenames ignored), last touched more than 30 days ago. | 393 | 118M |
| **PROBABLE-PROVENANCE** | - | Not cited by name, but it sits in an experiment directory (below bench/, docs/, ...) holding anchored files, or under a folder whose anchored README/write-up covers it -- usually an output of that experiment written under a computed name and cited as a folder. Verification judged 5/6 such files provenance: not a candidate. | 105 | 283K |
| **UNCITED-RECENT** | low | Not cited, touched within 30 days: may be work in progress. | 28 | 1M |
| **KEEP-CITED** | - | Anchored directly or transitively. Not a candidate. | 1184 | 96M |

**1480 candidates** of 2927 files (78 high, 724 medium, 678 low confidence).

## Candidates by directory

| directory | TEMP | REFUTED-MODULE | SUPERSEDED-PORTED | SUPERSEDED-VERSION | SUPERSEDED-CITED | LEGACY-ASSEMBLER | AMBIGUOUS-CITE | ORPHAN | UNCITED-STALE | UNCITED-RECENT | kept |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `bench/` | 39 |  | 23 | 14 | 2 | 16 | 152 | 107 | 157 | 5 | 690 |
| `bench/soto/` | 2 |  |  | 2 |  |  | 322 |  |  |  | 131 |
| `analysis/family_graphs/` |  |  |  |  |  |  | 24 | 72 | 46 |  | 62 |
| `bench/psv_split_stringtie/` |  |  |  |  |  |  |  |  | 80 |  | 0 |
| `scripts/` |  |  |  |  |  | 22 | 1 | 11 | 27 |  | 3 |
| `bench/graph_align_experiment_structural/` |  |  |  |  |  |  | 50 |  |  |  | 0 |
| `tools/demo/` |  |  |  |  |  | 6 | 15 | 22 | 7 |  | 1 |
| `docs/` |  |  |  | 6 | 8 |  | 11 |  | 1 | 15 | 112 |
| `(root)` | 23 |  |  | 4 |  |  |  | 2 | 4 |  | 7 |
| `bench/o1_gene_family_audit/` |  |  |  |  |  |  | 5 |  | 24 |  | 3 |
| `bench/fixtures/` |  |  |  |  |  |  | 24 |  |  |  | 2 |
| `bench/slides/` |  |  |  |  |  |  |  | 16 |  |  | 3 |
| `bench/ggo19_needy_top15_refs/` |  |  |  |  |  |  |  | 11 |  |  | 0 |
| `bench/o1_provenance_witness_prototype/` |  |  |  |  |  |  | 1 |  | 10 |  | 13 |
| `figures/` |  |  |  |  |  |  | 1 |  | 10 |  | 0 |
| `tools/trace_analysis/` |  |  |  |  |  | 7 |  |  | 3 |  | 0 |
| `bench/o1_fresh_emission_validation/` | 2 |  |  |  |  |  | 1 |  | 4 |  | 15 |
| `bench/sim/` |  |  |  |  |  |  | 4 |  | 3 |  | 0 |
| `tools/` |  |  |  |  |  |  |  | 3 | 4 |  | 3 |
| `wf3/ledger/` |  |  |  |  |  |  |  |  |  | 7 | 0 |
| `bench/identifiability_boundary/` |  |  |  |  |  |  |  | 5 | 1 |  | 0 |
| `bench/ggo19_needy_top5_refs/` |  |  |  |  |  |  |  | 5 |  |  | 0 |
| `tests/mechanism/` | 1 |  |  |  |  |  |  |  | 4 |  | 0 |
| `bench/o1_expanded_family_audit/` |  |  |  |  |  |  | 4 |  |  |  | 75 |
| `bench/tandem_attribution/` | 1 |  |  |  |  |  | 1 | 1 | 1 |  | 3 |
| `bench/compara_e/` |  |  |  |  |  |  |  | 1 | 1 |  | 0 |
| `bench/te_bridge_check/` |  |  |  |  |  |  |  |  | 2 |  | 0 |
| `test_data/vg_hmm/` |  |  |  |  |  |  |  | 2 |  |  | 0 |
| `bench/crossspecies/` | 1 |  |  |  |  |  |  |  |  |  | 27 |
| `bench/family_filter_chr19/` |  |  |  |  |  |  |  |  | 1 |  | 0 |
| `bench/family_filter_full/` |  |  |  |  |  |  |  |  | 1 |  | 0 |
| `bench/hidden_copy/` |  |  |  |  |  |  |  |  | 1 |  | 0 |
| `bench/jaccard_sweep/` |  |  |  |  |  |  |  |  | 1 |  | 0 |
| `bench/mechanism/` | 1 |  |  |  |  |  |  |  |  |  | 26 |
| `bench/mosaic_discriminator/` | 1 |  |  |  |  |  |  |  |  |  | 3 |
| `bench/o1_outgroup_rooting_poc/` |  |  |  |  |  |  | 1 |  |  |  | 20 |
| `bench/psv_sizing/` | 1 |  |  |  |  |  |  |  |  |  | 7 |
| `examples/` |  |  |  |  |  |  |  |  |  | 1 | 0 |
| `tests/regression/` |  |  |  |  |  | 1 |  |  |  |  | 2 |
| `bench/__pycache__/` | 1 |  |  |  |  |  |  |  |  |  | 0 |

## TEMP (73)

- `1` [ignored, high]
- `TMPPATH` [ignored, high]
- `as_margin.py` [ignored, high]
- `bench/__pycache__/` [ignored, high] (cache-dir:10-files)
- `bench/crossspecies/__pycache__/` [ignored, high] (cache-dir:3-files)
- `bench/err.log` [ignored, high] (match:glob-wide)
- `bench/family_level_pr_current.log` [ignored, high] (match:glob-wide)
- `bench/filter_families_by_psv.log` [ignored, high] (match:glob-wide)
- `bench/ggo19_needy_top15_run_STRG_125.stderr.log` [ignored, high] (match:glob-wide)
- `bench/ggo19_needy_top15_run_STRG_210.stderr.log` [ignored, high] (match:glob-wide)
- `bench/ggo19_needy_top15_run_STRG_300.stderr.log` [ignored, high] (match:glob-wide)
- `bench/ggo19_needy_top15_run_STRG_52.stderr.log` [ignored, high] (match:glob-wide)
- `bench/ggo19_needy_top15_run_STRG_95.stderr.log` [ignored, high] (match:glob-wide)
- `bench/jfp_study_st_out.gtf` [ignored, high] (match:glob-wide)
- `bench/mechanism/__pycache__/` [ignored, high] (cache-dir:10-files)
- `bench/mosaic_discriminator/__pycache__/` [ignored, high] (cache-dir:2-files)
- `bench/o1_fresh_emission_validation/GGO.fresh.run.log` [ignored, high]
- `bench/o1_fresh_emission_validation/HSA.fresh.run.log` [ignored, high]
- `bench/psv_sizing/__pycache__/` [ignored, high] (cache-dir:1-files)
- `bench/sn_sweep2.log` [ignored, high] (match:glob-wide)
- `bench/soto/.pytest_cache/` [ignored, high] (cache-dir:4-files)
- `bench/soto/__pycache__/` [ignored, high] (cache-dir:23-files)
- `bench/tandem_attribution/__pycache__/` [ignored, high] (cache-dir:2-files)
- `bench/vg_family_prototype.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_antisense.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_cdhit.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_cdhit90.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_cdhit_eval.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_eval.log` [ignored, high] (byte-identical-to:bench/vg_family_prototype_o1vg_cdhit_eval.log;match:glob-wide)
- `bench/vg_family_prototype_exact.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_fp_characterize.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_id95.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_members80.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_members80_eval.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_min2.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_mmseqs.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_mmseqs_eval.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_mmseqs_pairs.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_mmseqs_pairs85.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_mmseqs_pairs85_eval.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_mmseqs_pairs_eval.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_mmseqs_single.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_o1vg.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_o1vg_cdhit.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_o1vg_cdhit_eval.log` [ignored, high] (byte-identical-to:bench/vg_family_prototype_eval.log;match:glob-wide)
- `bench/vg_family_prototype_o1vg_exact_eval.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_protcov_eval.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_protpure_eval.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_repeatgate.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_repeatgate_eval.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_seeded.log` [ignored, high] (match:glob-wide)
- `bench/vg_family_prototype_vsearch.log` [ignored, high] (match:glob-wide)
- `blockcmp.py` [ignored, high]
- `build_err.txt` [ignored, high] (dead-path-refs:28/43)
- `chry.py` [ignored, high]
- `dbg.py` [ignored, high]
- `exon_id2.py` [ignored, high]
- `final_tally.py` [ignored, high]
- `nm_margin.py` [ignored, high]
- `nm_values.py` [ignored, high]
- `probe_ends4.py` [ignored, high]
- `probe_err.txt` [ignored, high] (dead-path-refs:28/30)
- `probe_final.py` [ignored, high]
- `probe_flank_div.py` [ignored, high]
- `probe_k0_confirm.py` [ignored, high]
- `probe_out.txt` [ignored, high]
- `probe_pair0_tss.py` [ignored, high]
- `probe_pair1.py` [ignored, high]
- `psv_panel.py` [ignored, high]
- `ref_psv.py` [ignored, high]
- `resolve_reads.py` [ignored, high]
- `shared_exon_dump.patch.txt` [ignored, high]
- `tests/mechanism/__pycache__/` [ignored, high] (cache-dir:4-files)

## SUPERSEDED-PORTED (23)

- `bench/allele_specific_junctions.py` [tracked, low] → superseded by `src/rustle/vg_family/allele_specific_junctions.rs` (still-used-by:bench/asj_aggregate.py;asj-objective-dropped;match:exact+glob-wide+suffix)
- `bench/asj_genetic_core.py` [tracked, low] → superseded by `src/rustle/vg_family/asj_genetic_core.rs,src/rustle/vg_family/mod.rs` (still-used-by:bench/gen_asj_genetic_core_fixture.py;asj-objective-dropped;match:exact+glob-wide+import+suffix)
- `bench/asj_strand_bias.py` [tracked, low] → superseded by `src/rustle/vg_family/asj_strand_bias.rs,src/rustle/vg_family/mod.rs` (still-used-by:bench/asj_aggregate.py;asj-objective-dropped;match:exact+glob-wide+import+suffix)
- `bench/asj_verify.py` [tracked, low] → superseded by `src/rustle/vg_family/asj_verify.rs,src/rustle/vg_family/mod.rs` (still-used-by:bench/asj_fig.py;asj-objective-dropped;match:exact+glob-wide+import+suffix)
- `bench/copy_assign.py` [tracked, low] → superseded by `src/rustle/vg_family/copy_assign.rs,src/rustle/vg_family/copy_assign_pipeline.rs,src/rustle/vg_family/denovo_pipeline.rs,src/rustle/vg_family/mod.rs,src/rustle/vg_family/o2_margin_gate.rs` (partial-port;still-used-by:bench/align_error_dna_test.py;match:exact+glob-wide+import+suffix)
- `bench/denovo_assemble_gate.py` [tracked, low] → superseded by `src/rustle/vg_family/denovo_assemble.rs` (still-used-by:bench/genome_rna_overlay_readcontent.py;match:exact+glob-wide+suffix)
- `bench/denovo_shared_def.py` [tracked, low] → superseded by `src/rustle/vg_family/shared_definition.rs` (still-used-by:bench/layer_order/lattice_check_c2.py;named-beside-stale-language-by:docs/ACTIVE_WORKING_SET.md;match:exact+glob-wide+import+suffix)
- `bench/family_copy_number.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs` (partial-port;still-used-by:bench/rna_copy_number_depth.py;match:exact+glob-wide+suffix)
- `bench/family_er_pr.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs` (partial-port;still-used-by:bench/divergence_floor.py;match:exact+glob-wide+import+suffix)
- `bench/family_rescue.py` [tracked, medium] → superseded by `src/rustle/vg_family/family_rescue.rs,src/rustle/vg_family/rescue_pipeline.rs` (match:exact+glob-wide+suffix)
- `bench/family_rna_refine.py` [tracked, low] → superseded by `src/bin/family_define.rs,src/rustle/vg_family/driver.rs,src/rustle/vg_family/mod.rs` (still-used-by:bench/divergence_floor.py;match:exact+glob-wide+import+suffix)
- `bench/genome_family_def.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs` (still-used-by:bench/colinear_multiexon_gate.py;match:exact+glob-wide+import+suffix)
- `bench/multi_repeat_bridge_gate.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/multi_repeat_bridge.rs` (still-used-by:bench/divergence_floor.py;match:exact+glob-wide+import+suffix)
- `bench/o2_vg_visualization.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/o2_materialize.rs` (partial-port;still-used-by:bench/gen_o2_margin_gate_fixture.py;match:exact+glob-wide+import+suffix)
- `bench/o3_flag_pass.py` [tracked, medium] → superseded by `src/bin/copy_assign.rs,src/rustle/vg_family/o3_flag_pass.rs` (named-beside-stale-language-by:docs/ACTIVE_WORKING_SET.md;match:exact+glob-wide+suffix)
- `bench/poa_family_definition.py` [tracked, low] → superseded by `src/rustle/vg_family/family_graph.rs` (still-used-by:bench/candidate_generation_recall.py;match:exact+glob-wide+import+suffix)
- `bench/psv_graph_genomewide.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/o2_columns.rs` (partial-port;still-used-by:bench/a1_read_sda_smoketest.py;named-beside-stale-language-by:memory:project_psv_aware_vg.md;match:exact+glob-wide+import+suffix)
- `bench/recombinant_abstain.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/recombinant_abstain.rs` (still-used-by:bench/gen_recombinant_abstain_fixture.py;match:exact+glob-wide+import+suffix)
- `bench/recombinant_split.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/recombinant_split.rs` (still-used-by:bench/divergence_floor.py;match:exact+glob-wide+import+suffix)
- `bench/recombination_bridge_detector.py` [tracked, low] → superseded by `src/rustle/vg_family/bridge_detector.rs,src/rustle/vg_family/mod.rs` (still-used-by:bench/colinear_multiexon_gate.py;match:exact+glob-wide+import+suffix)
- `bench/rna_only_edge_oracle.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs` (still-used-by:bench/colinear_multiexon_gate.py;match:exact+glob-wide+import+suffix)
- `bench/twopass_denovo_gw_pass1.py` [tracked, medium] → superseded by `src/rustle/vg_family/denovo_assemble.rs` (match:exact+glob-wide+suffix)
- `bench/vg_repeat_catalog.py` [tracked, low] → superseded by `src/rustle/vg_family/minimizers.rs,src/rustle/vg_family/mod.rs,src/rustle/vg_family/repeat_catalog.rs` (still-used-by:bench/family_rna_refine.py;named-beside-stale-language-by:docs/ACTIVE_WORKING_SET.md;match:exact+glob-wide+import+suffix)

## SUPERSEDED-VERSION (26)

- `bench/analyze_jfps.py` [tracked, medium] → superseded by `bench/analyze_jfps_v5.py` (series:v;match:glob-wide)
- `bench/analyze_jfps_v2.py` [tracked, medium] → superseded by `bench/analyze_jfps_v5.py` (series:v;match:glob-wide)
- `bench/analyze_jfps_v3.py` [tracked, medium] → superseded by `bench/analyze_jfps_v5.py` (series:v;match:glob-wide)
- `bench/analyze_jfps_v4.py` [tracked, medium] → superseded by `bench/analyze_jfps_v5.py` (series:v;match:glob-wide+suffix)
- `bench/family_def_refute_edgebetw.py` [tracked, medium] → superseded by `bench/family_def_refute_edgebetw3.py` (series:n;match:glob-wide)
- `bench/family_def_refute_edgebetw2.py` [tracked, medium] → superseded by `bench/family_def_refute_edgebetw3.py` (series:n;match:glob-wide)
- `bench/ggo19_tlen_norm.gtf` [ignored, medium] → superseded by `bench/ggo19_tlen_norm_final.gtf` (series:fix;match:glob-wide)
- `bench/gw_rebuild_v3.sh` [tracked, medium] → superseded by `bench/gw_rebuild_v4.sh` (series:v;match:glob-wide)
- `bench/jfp_study_denovo.gtf` [ignored, medium] → superseded by `bench/jfp_study_denovo_fix.gtf` (series:fix;match:exact+glob-wide)
- `bench/no_hs_se_deplete.gtf` [ignored, medium] → superseded by `bench/no_hs_se_deplete2.gtf` (series:n;match:glob-wide)
- `bench/no_hs_spare_denovo.gtf` [ignored, medium] → superseded by `bench/no_hs_spare_denovo2.gtf` (series:n;match:glob-wide)
- `bench/refute_cov_min.py` [tracked, medium] → superseded by `bench/refute_cov_min2.py` (series:n;match:glob-wide)
- `bench/se_last_v2.gtf` [ignored, medium] → superseded by `bench/se_last_v3.gtf` (series:v;match:glob-wide)
- `bench/sn_sweep.log` [ignored, medium] → superseded by `bench/sn_sweep2.log` (series:n;match:glob-wide+prefix)
- `bench/soto/member_attribution.py` [tracked, medium] → superseded by `bench/soto/member_attribution_final.py` (series:fix;match:dir-large)
- `bench/soto/member_attribution.tsv` [tracked, medium] → superseded by `bench/soto/member_attribution_final.tsv` (series:fix;match:dir-large+suffix)
- `docs/sweep_v1_families_2026-09-05.tsv` [tracked, high] → superseded by `docs/sweep_v10_families_2026-09-05.tsv` (series:v)
- `docs/sweep_v1_units_2026-09-05.tsv` [tracked, medium] → superseded by `docs/sweep_v2_units_2026-09-05.tsv` (series:v)
- `docs/sweep_v2_families_2026-09-05.tsv` [tracked, high] → superseded by `docs/sweep_v10_families_2026-09-05.tsv` (series:v)
- `docs/sweep_v3_families_2026-09-05.tsv` [tracked, high] → superseded by `docs/sweep_v10_families_2026-09-05.tsv` (series:v)
- `docs/sweep_v8_families_2026-09-05.tsv` [tracked, high] → superseded by `docs/sweep_v10_families_2026-09-05.tsv` (series:v)
- `docs/sweep_v9_families_2026-09-05.tsv` [tracked, high] → superseded by `docs/sweep_v10_families_2026-09-05.tsv` (series:v)
- `exon_id.py` [ignored, medium] → superseded by `exon_id2.py` (series:n)
- `probe_ends.py` [ignored, medium] → superseded by `probe_ends4.py` (series:n)
- `probe_ends2.py` [ignored, medium] → superseded by `probe_ends4.py` (series:n)
- `probe_ends3.py` [ignored, medium] → superseded by `probe_ends4.py` (series:n)

## Overlay: files tied to the dropped ASJ objective (31)

Not classified as dead by this flag alone — ASJ was dropped as an objective (memory, 2026-08-07), but the binaries still build. Decide as a scope question, not a cleanup one.

- `bench/ASJ.md` — ORPHAN
- `bench/allele_specific_junctions.py` — SUPERSEDED-PORTED
- `bench/allele_specific_junctions_multisnp.py` — KEEP-CITED
- `bench/asj_aggregate.py` — KEEP-CITED
- `bench/asj_calls.tsv` — KEEP-CITED
- `bench/asj_calls_strandbias.tsv` — KEEP-CITED
- `bench/asj_calls_verified.tsv` — KEEP-CITED
- `bench/asj_evidence.py` — KEEP-CITED
- `bench/asj_fig.py` — KEEP-CITED
- `bench/asj_findings.png` — KEEP-CITED
- `bench/asj_genetic_core.py` — SUPERSEDED-PORTED
- `bench/asj_genetic_core.tsv` — KEEP-CITED
- `bench/asj_motif_check.py` — KEEP-CITED
- `bench/asj_strand_bias.py` — SUPERSEDED-PORTED
- `bench/asj_verify.py` — SUPERSEDED-PORTED
- `bench/asjm_aggregate.py` — KEEP-CITED
- `bench/asjm_calls.tsv` — KEEP-CITED
- `bench/asjm_fig.py` — KEEP-CITED
- `bench/asjm_findings.png` — KEEP-CITED
- `bench/gen_asj_genetic_core_fixture.py` — KEEP-CITED
- `bench/gen_asj_strand_bias_fixture.py` — KEEP-CITED
- `bench/gen_asj_verify_fixture.py` — KEEP-CITED
- `src/bin/asj.rs` — PROTECTED
- `src/bin/asj_verify.rs` — PROTECTED
- `src/rustle/vg_family/allele_specific_junctions.rs` — PROTECTED
- `src/rustle/vg_family/asj_genetic_core.rs` — PROTECTED
- `src/rustle/vg_family/asj_strand_bias.rs` — PROTECTED
- `src/rustle/vg_family/asj_verify.rs` — PROTECTED
- `src/rustle/vg_family/testdata/asj_genetic_core_fixture.json` — KEEP-CITED
- `src/rustle/vg_family/testdata/asj_strand_bias_fixture.json` — KEEP-CITED
- `src/rustle/vg_family/testdata/asj_verify_fixture.json` — KEEP-CITED

## Known blind spots

- Paths built at run time (`f"bench/{name}.tsv"`, shell loops) are invisible to a token scan: such outputs land in UNCITED-IN-ANCHORED-DIR or UNCITED-*, never KEEP.
- A citation proves a file was *named*, not that the naming text is still true; an anchor in a retracted ledger section still anchors (see `ledger_sections`).
- Data outside the repo (`/mnt/linuxdisk`, BAMs) is not audited; `.worktrees/`, `tools/stringtie` (submodule), `.remember/`, `.superpowers/` are excluded.
- Rust module reachability is not re-derived here — see `docs/MODULE_STATUS.md`.
