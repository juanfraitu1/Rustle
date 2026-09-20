# Cleanup candidates — likely dead / likely superseded files

Generated 2026-09-19 at `dna-from-genome@14d7c874` by `tools/audit_cleanup_candidates.py` (re-run it; this file is overwritten). **Read-only: nothing was moved, edited or deleted.** Full per-file table: `docs/cleanup_candidates.tsv` (filter on `class` and `confidence`).

⚠ A mark is a *candidate*, not a verdict. Before deleting anything: (1) grep the path once more, (2) check the `ledger_sections` / `anchor_citers` columns, (3) prefer `git mv` into an archive directory over `rm` for anything tracked, (4) remember `bench/` data can be slow to regenerate (AGENTS.md §2).

## How a file is kept

A file is **anchored** if a path-like token naming it (full path, unique path suffix, unique basename, output prefix, a glob matching exactly 1 file, an enclosing directory of ≤40 files, or a Python `import` resolved against the importing script's directory and its parents) appears in an anchor source: `docs/o1_ledger.md`, `docs/NEGATIVE_RESULTS_REGISTER.md`, `docs/PREREG_*.md`, top-level `docs/*.md`, `docs/experiments/*.md`, `README.md`, the auto-memory directory, Rust `src/`/`tests/`, or `Cargo.toml`. Anchoring then propagates: a file named by an anchored script or anchored markdown write-up is anchored, index sidecars (`.fai`, `.bai`) follow their parent, and a script that names anchored data (its generator — the naming line looks like a write) is anchored; a script that merely reads anchored data is not. An ambiguous bare basename (e.g. `reads.fa`, found in many directories) never anchors on its own. Citations from `docs/archive/`, `docs/superpowers/`, `AGENTS.md` and un-anchored scripts are *weak*: recorded, not counted.

## Rules (applied in this order; stale cutoff = 2026-08-20)

| class | confidence | rule | files | bytes |
|---|---|---|---:|---:|
| **PROTECTED** | - | build/config files, Rust sources and tests (Rust reachability is docs/MODULE_STATUS.md's job, enforced by module_status_tests), and the anchor docs themselves. Never a candidate. | 158 | 6M |
| **TEMP** | high | Python/pytest caches; untracked or git-ignored files at the repo root; untracked *.log / *err* / *out* / *.patch.txt / checkpoint files that nothing cites. | 0 | 0B |
| **REFUTED-MODULE** | medium | Rust module whose `//! **STATUS:**` header is REFUTED and that no other file names (a REFUTED module that is still imported, e.g. collapse_gate.rs, stays PROTECTED). | 0 | 0B |
| **SUPERSEDED-PORTED** | medium | Python script that a Rust source line declares it ports ('Port of', 'Faithful Rust port of', 'Mirrors', 'migration'); 'low' when only a function or part is ported (`x.py::f`, `x.py loaders`). The Python may still serve as a parity oracle or golden-fixture generator -- check tests before deleting. | 23 | 580K |
| **SUPERSEDED-VERSION** | medium | Older member of a version series in the same directory (_v1.._vN, foo/foo2/foo3, dated _YYYY-MM-DD copies, foo vs foo_fix/_final/_new) that nothing anchors. 'high' when the newest member IS anchored. | 9 | 75K |
| **SUPERSEDED-CITED** | low | Older member of a version series that IS anchored: provenance for a recorded result -- archive, don't delete. | 12 | 540K |
| **LEGACY-ASSEMBLER** | medium | Not anchored, and its path or first 200 lines name StringTie-era assembler machinery (bundle/transfrag/parity/gffcompare/...) -- the assembler layer was retired (docs/RETIREMENT_AND_MIGRATION.md). | 52 | 482K |
| **AMBIGUOUS-CITE** | low | Not anchored; an anchor source names it only by a bare basename shared by several files, a directory too large (>40 files) or a glob matching several files -- may or may not mean this copy. | 620 | 31M |
| **ORPHAN** | medium | Not anchored; named only by files that are themselves not anchored (e.g. a figure named only by its un-cited plotting script). | 256 | 38M |
| **UNCITED-STALE** | medium | Not named by anything (wide globs / big directories / scripts' ambiguous basenames ignored), last touched more than 30 days ago. | 393 | 118M |
| **PROBABLE-PROVENANCE** | - | Not cited by name, but it sits in an experiment directory (below bench/, docs/, ...) holding anchored files, or under a folder whose anchored README/write-up covers it -- usually an output of that experiment written under a computed name and cited as a folder. Verification judged 5/6 such files provenance: not a candidate. | 107 | 294K |
| **UNCITED-RECENT** | low | Not cited, touched within 30 days: may be work in progress. | 30 | 1M |
| **KEEP-CITED** | - | Anchored directly or transitively. Not a candidate. | 1218 | 222M |

**1395 candidates** of 2878 files (0 high, 710 medium, 685 low confidence).

## Candidates by directory

| directory | TEMP | REFUTED-MODULE | SUPERSEDED-PORTED | SUPERSEDED-VERSION | SUPERSEDED-CITED | LEGACY-ASSEMBLER | AMBIGUOUS-CITE | ORPHAN | UNCITED-STALE | UNCITED-RECENT | kept |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `bench/` |  |  |  |  |  |  | 152 | 106 | 157 | 5 | 691 |
| `bench/soto/` |  |  |  |  |  |  | 322 |  |  |  | 131 |
| `analysis/family_graphs/` |  |  |  |  |  |  | 24 | 72 | 46 |  | 62 |
| `bench/psv_split_stringtie/` |  |  |  |  |  |  |  |  | 80 |  | 0 |
| `archive/bench/` |  |  | 23 | 5 | 2 | 16 | 4 |  |  | 1 | 0 |
| `bench/graph_align_experiment_structural/` |  |  |  |  |  |  | 50 |  |  |  | 0 |
| `tools/demo/` |  |  |  |  |  |  | 15 | 22 | 7 |  | 1 |
| `scripts/` |  |  |  |  |  |  | 1 | 11 | 27 |  | 3 |
| `bench/o1_gene_family_audit/` |  |  |  |  |  |  | 5 |  | 24 |  | 3 |
| `docs/` |  |  |  |  |  |  | 10 |  | 1 | 16 | 112 |
| `bench/fixtures/` |  |  |  |  |  |  | 24 |  |  |  | 2 |
| `archive/scripts/` |  |  |  |  |  | 22 |  |  |  |  | 0 |
| `bench/slides/` |  |  |  |  |  |  |  | 16 |  |  | 3 |
| `archive/tools/` |  |  |  |  |  | 13 |  |  |  |  | 0 |
| `bench/ggo19_needy_top15_refs/` |  |  |  |  |  |  |  | 11 |  |  | 0 |
| `bench/o1_provenance_witness_prototype/` |  |  |  |  |  |  | 1 |  | 10 |  | 13 |
| `figures/` |  |  |  |  |  |  | 1 |  | 10 |  | 0 |
| `archive/docs/` |  |  |  | 4 | 6 |  |  |  |  |  | 4 |
| `bench/sim/` |  |  |  |  |  |  | 4 |  | 3 |  | 0 |
| `tools/` |  |  |  |  |  |  |  | 3 | 4 |  | 4 |
| `wf3/ledger/` |  |  |  |  |  |  |  |  |  | 7 | 0 |
| `bench/identifiability_boundary/` |  |  |  |  |  |  |  | 5 | 1 |  | 0 |
| `(root)` |  |  |  |  |  |  |  | 1 | 4 |  | 8 |
| `bench/ggo19_needy_top5_refs/` |  |  |  |  |  |  |  | 5 |  |  | 0 |
| `bench/o1_fresh_emission_validation/` |  |  |  |  |  |  | 1 |  | 4 |  | 15 |
| `archive/untracked/` |  |  |  |  | 4 |  |  |  |  |  | 29 |
| `bench/o1_expanded_family_audit/` |  |  |  |  |  |  | 4 |  |  |  | 75 |
| `tests/mechanism/` |  |  |  |  |  |  |  |  | 4 |  | 0 |
| `bench/tandem_attribution/` |  |  |  |  |  |  | 1 | 1 | 1 |  | 3 |
| `tools/trace_analysis/` |  |  |  |  |  |  |  |  | 3 |  | 0 |
| `bench/compara_e/` |  |  |  |  |  |  |  | 1 | 1 |  | 0 |
| `bench/te_bridge_check/` |  |  |  |  |  |  |  |  | 2 |  | 0 |
| `test_data/vg_hmm/` |  |  |  |  |  |  |  | 2 |  |  | 0 |
| `archive/tests/` |  |  |  |  |  | 1 |  |  |  |  | 0 |
| `bench/family_filter_chr19/` |  |  |  |  |  |  |  |  | 1 |  | 0 |
| `bench/family_filter_full/` |  |  |  |  |  |  |  |  | 1 |  | 0 |
| `bench/hidden_copy/` |  |  |  |  |  |  |  |  | 1 |  | 0 |
| `bench/jaccard_sweep/` |  |  |  |  |  |  |  |  | 1 |  | 0 |
| `bench/o1_outgroup_rooting_poc/` |  |  |  |  |  |  | 1 |  |  |  | 20 |
| `examples/` |  |  |  |  |  |  |  |  |  | 1 | 0 |

## SUPERSEDED-PORTED (23)

- `archive/bench/allele_specific_junctions.py` [tracked, low] → superseded by `src/rustle/vg_family/allele_specific_junctions.rs` (still-used-by:archive/bench/asj_strand_bias.py;asj-objective-dropped;match:glob-wide+suffix)
- `archive/bench/asj_genetic_core.py` [tracked, low] → superseded by `src/rustle/vg_family/asj_genetic_core.rs,src/rustle/vg_family/mod.rs` (still-used-by:bench/gen_asj_genetic_core_fixture.py;asj-objective-dropped;match:glob-wide+import+suffix)
- `archive/bench/asj_strand_bias.py` [tracked, low] → superseded by `src/rustle/vg_family/asj_strand_bias.rs,src/rustle/vg_family/mod.rs` (still-used-by:bench/asj_aggregate.py;asj-objective-dropped;match:glob-wide+import+suffix)
- `archive/bench/asj_verify.py` [tracked, low] → superseded by `src/rustle/vg_family/asj_verify.rs,src/rustle/vg_family/mod.rs` (still-used-by:bench/asj_fig.py;asj-objective-dropped;match:glob-wide+import+suffix)
- `archive/bench/copy_assign.py` [tracked, low] → superseded by `src/rustle/vg_family/copy_assign.rs,src/rustle/vg_family/copy_assign_pipeline.rs,src/rustle/vg_family/denovo_pipeline.rs,src/rustle/vg_family/mod.rs,src/rustle/vg_family/o2_margin_gate.rs` (partial-port;still-used-by:archive/bench/o2_vg_visualization.py;match:glob-wide+import+suffix)
- `archive/bench/denovo_assemble_gate.py` [tracked, low] → superseded by `src/rustle/vg_family/denovo_assemble.rs` (still-used-by:bench/genome_rna_overlay_readcontent.py;match:glob-wide+suffix)
- `archive/bench/denovo_shared_def.py` [tracked, low] → superseded by `src/rustle/vg_family/shared_definition.rs` (still-used-by:bench/layer_order/lattice_check_c2.py;match:exact+glob-wide+import+suffix)
- `archive/bench/family_copy_number.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs` (partial-port;still-used-by:bench/rna_copy_number_depth.py;match:glob-wide+suffix)
- `archive/bench/family_er_pr.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs` (partial-port;still-used-by:archive/bench/family_rna_refine.py;match:glob-wide+import+suffix)
- `archive/bench/family_rescue.py` [tracked, medium] → superseded by `src/rustle/vg_family/family_rescue.rs,src/rustle/vg_family/rescue_pipeline.rs` (match:glob-wide+suffix)
- `archive/bench/family_rna_refine.py` [tracked, low] → superseded by `src/bin/family_define.rs,src/rustle/vg_family/driver.rs,src/rustle/vg_family/mod.rs` (still-used-by:archive/bench/multi_repeat_bridge_gate.py;match:glob-wide+import+suffix)
- `archive/bench/genome_family_def.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs` (still-used-by:archive/bench/family_er_pr.py;match:glob-wide+import+suffix)
- `archive/bench/multi_repeat_bridge_gate.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/multi_repeat_bridge.rs` (still-used-by:archive/bench/family_rna_refine.py;match:glob-wide+import+suffix)
- `archive/bench/o2_vg_visualization.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/o2_materialize.rs` (partial-port;still-used-by:bench/gen_o2_margin_gate_fixture.py;match:glob-wide+import+suffix)
- `archive/bench/o3_flag_pass.py` [tracked, medium] → superseded by `src/bin/copy_assign.rs,src/rustle/vg_family/o3_flag_pass.rs` (match:exact+glob-wide+suffix)
- `archive/bench/poa_family_definition.py` [tracked, low] → superseded by `src/rustle/vg_family/family_graph.rs` (still-used-by:archive/bench/family_rescue.py;match:glob-wide+import+suffix)
- `archive/bench/psv_graph_genomewide.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/o2_columns.rs` (partial-port;still-used-by:archive/bench/family_copy_number.py;named-beside-stale-language-by:memory:project_psv_aware_vg.md;match:glob-wide+import+suffix)
- `archive/bench/recombinant_abstain.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/recombinant_abstain.rs` (still-used-by:archive/bench/o2_vg_visualization.py;match:glob-wide+import+suffix)
- `archive/bench/recombinant_split.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs,src/rustle/vg_family/recombinant_split.rs` (still-used-by:archive/bench/family_rna_refine.py;match:glob-wide+import+suffix)
- `archive/bench/recombination_bridge_detector.py` [tracked, low] → superseded by `src/rustle/vg_family/bridge_detector.rs,src/rustle/vg_family/mod.rs` (still-used-by:archive/bench/multi_repeat_bridge_gate.py;match:glob-wide+import+suffix)
- `archive/bench/rna_only_edge_oracle.py` [tracked, low] → superseded by `src/rustle/vg_family/mod.rs` (still-used-by:archive/bench/family_rna_refine.py;match:glob-wide+import+suffix)
- `archive/bench/twopass_denovo_gw_pass1.py` [tracked, medium] → superseded by `src/rustle/vg_family/denovo_assemble.rs` (match:glob-wide+suffix)
- `archive/bench/vg_repeat_catalog.py` [tracked, low] → superseded by `src/rustle/vg_family/minimizers.rs,src/rustle/vg_family/mod.rs,src/rustle/vg_family/repeat_catalog.rs` (still-used-by:archive/bench/family_rna_refine.py;match:exact+glob-wide+import+suffix)

## SUPERSEDED-VERSION (9)

- `archive/bench/analyze_jfps.py` [tracked, medium] → superseded by `archive/bench/analyze_jfps_v5.py` (series:v;match:glob-wide)
- `archive/bench/analyze_jfps_v2.py` [tracked, medium] → superseded by `archive/bench/analyze_jfps_v5.py` (series:v;match:glob-wide)
- `archive/bench/analyze_jfps_v3.py` [tracked, medium] → superseded by `archive/bench/analyze_jfps_v5.py` (series:v;match:glob-wide)
- `archive/bench/analyze_jfps_v4.py` [tracked, medium] → superseded by `archive/bench/analyze_jfps_v5.py` (series:v;match:glob-wide+suffix)
- `archive/bench/family_def_refute_edgebetw.py` [tracked, medium] → superseded by `archive/bench/family_def_refute_edgebetw2.py` (series:n;match:glob-wide)
- `archive/docs/sweep_v1_families_2026-09-05.tsv` [tracked, medium] → superseded by `archive/docs/sweep_v9_families_2026-09-05.tsv` (series:v)
- `archive/docs/sweep_v2_families_2026-09-05.tsv` [tracked, medium] → superseded by `archive/docs/sweep_v9_families_2026-09-05.tsv` (series:v)
- `archive/docs/sweep_v3_families_2026-09-05.tsv` [tracked, medium] → superseded by `archive/docs/sweep_v9_families_2026-09-05.tsv` (series:v)
- `archive/docs/sweep_v8_families_2026-09-05.tsv` [tracked, medium] → superseded by `archive/docs/sweep_v9_families_2026-09-05.tsv` (series:v)

## Overlay: files tied to the dropped ASJ objective (31)

Not classified as dead by this flag alone — ASJ was dropped as an objective (memory, 2026-08-07), but the binaries still build. Decide as a scope question, not a cleanup one.

- `archive/bench/allele_specific_junctions.py` — SUPERSEDED-PORTED
- `archive/bench/asj_genetic_core.py` — SUPERSEDED-PORTED
- `archive/bench/asj_strand_bias.py` — SUPERSEDED-PORTED
- `archive/bench/asj_verify.py` — SUPERSEDED-PORTED
- `bench/ASJ.md` — ORPHAN
- `bench/allele_specific_junctions_multisnp.py` — KEEP-CITED
- `bench/asj_aggregate.py` — KEEP-CITED
- `bench/asj_calls.tsv` — KEEP-CITED
- `bench/asj_calls_strandbias.tsv` — KEEP-CITED
- `bench/asj_calls_verified.tsv` — KEEP-CITED
- `bench/asj_evidence.py` — KEEP-CITED
- `bench/asj_fig.py` — KEEP-CITED
- `bench/asj_findings.png` — KEEP-CITED
- `bench/asj_genetic_core.tsv` — KEEP-CITED
- `bench/asj_motif_check.py` — KEEP-CITED
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
