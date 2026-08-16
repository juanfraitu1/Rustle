# Rustle retirement + migration plan

Status: **MARKING PASS COMPLETE** (retirement side: comments only — no logic changed, no files removed).
**Migration side: TIER 1 root item #1 (`genome_family_def.py` → `vg_family/family_definition.rs`) MIGRATED**
— `distinct_loci` + `refine_families`/R, option B deterministic splitter; see §5.1.
Branch: `vg/flow-capacity-apportionment`. Machine-readable map: `bench/RETIREMENT_MAP.tsv`.
**2026-07-13 audit pass:** safe dead code removed + `--vg-realign` docs corrected + intra-thesis
duplication inventoried — see **§6**. §4's 12-symbol tendril list re-confirmed byte-exact (grep 2026-07-13).

Rustle began as a Rust port of the StringTie long-read transcript **assembler**. The thesis has
since moved to **multi-copy gene FAMILIES**:

- **O1** — family DEFINITION (γ-quasi-clique homology components; copy number χ(H))
- **O2** — copy ASSIGNMENT: whether the evidence warrants assigning a read to a copy at all (PSV +
  junction + SUN, assign-or-abstain). ⚠ contested set is **alignment-score near-ties**, not MAPQ-0 —
  restated 2026-08-15, see `docs/copy_assignment_definition.md` §0
- **O3** — allele-specific junctions (ASJ)
- **O4** — reference-absent (collapsed / divergent) copies

The assembler + network-flow layer is **LEGACY / DEAD CODE**. This document is the removal plan and
the Python→Rust migration plan. Nothing here has been deleted yet; dead modules carry a top-of-file
`//! DEAD CODE --` banner and `lib.rs` carries a section banner. This is a *mark, don't remove* pass.

---

## 1. Thesis roots (what actually ships)

The four thesis binaries import **only** `rustle::genome` + `rustle::vg_family::*`:

- `src/bin/copy_assign.rs` (O2)
- `src/bin/gw_family_catalog.rs` (O1 catalog)
- `src/bin/asj.rs` (O3)
- `src/bin/asj_verify.rs` (O3 verification)

Everything else in `src/rustle/` is reachable only via (a) the thesis layer, (b) the legacy `rustle`
assembler bin (`main.rs` → `pipeline::run`), or (c) the four SHARED-tendril files.

---

## 2. Classification (dep-graph verified)

| Class | Meaning | Action |
|---|---|---|
| **KEEP-thesis** | `vg_family/*` (30 submodules), `genome.rs` | keep; build new work here |
| **KEEP-foundational** | `types.rs`, `bam.rs`, `util/*` — used by BOTH layers | keep; carry a *reverse tendril* to relocate |
| **SHARED-tendril** | `vg.rs`, `graph.rs`, `path_extract.rs`, `bundle.rs` — assembler files the thesis still imports from | EXTRACT the named export first, THEN delete |
| **DEAD-flow** | pure network flow — `max_flow.rs`, `global_flow.rs`, `parity/flow_iter_dump.rs` | **remove FIRST** |
| **DEAD-stringtie** | the rest of the assembler stack (63 files) | remove after tendrils extracted |
| **DEAD-orphan** | present but never declared as a module (not compiled): `strand_resolve.rs`, `assembly_constants.rs`, `junction_out.rs` | delete anytime (safe) |

**Marked this pass: 74 files** (3 DEAD-flow, 63 DEAD-stringtie, 3 DEAD-orphan, 4 SHARED-tendril,
1 legacy bin `src/bin/rustle.rs`) + the `lib.rs` section banner.
KEEP modules (`types`, `bam`, `genome`, `util/*`, `vg_family/*`) and the four thesis bins were NOT touched.

Build after marking: `cargo check` = **exit 0**; `cargo build --bin copy_assign --bin gw_family_catalog
--bin asj` = **exit 0** (27 pre-existing naming warnings, unchanged — comments add none).

### DEAD-flow (retire FIRST — user priority #1)
- `max_flow.rs` (110 KB, Edmonds-Karp path-capacity max flow). Users: junction_correction,
  path_extract, killed_junctions, parse_trflong_st, transfrag_process, pipeline, parity/flow_iter_dump.
  **Zero thesis use.**
- `global_flow.rs` (92 KB, greedy flow decomposition, `RUSTLE_GREEDY_DECOMPOSE`). Users: lib.rs +
  pipeline.rs only. **Zero thesis use.**
- `parity/flow_iter_dump.rs` (per-iteration flow dump; flow-coupled).

### DEAD-stringtie (63 files)
pipeline, main, transcript_filter, transfrag_process, graph_build, killed_junctions, map_reads,
bundle_builder, longtrim, junction_graph_st, parse_trflong_st, junction_correction, junctions,
junction_graph, coverage_trim, read_boundaries, cross_strand_predcluster, parallel_predprune,
single_exon_pileup, predcluster_st, stringtie_parity, terminal_oracle, ml_filter, gene_abundance,
merge_transcripts, merge_mode, nascent, polya, tss_tts, hard_boundaries, exon_merge, nodecov, bpcov,
report_losses, guide_path_fixer, per_bnode_graph, assembly_mode, treepat, gtf, reference_gtf,
vg_runtime, vg_ingestion, family_manifest, vg_eval, annotation_families, psv_fasta, graph_comparison,
debug_stage, snapshot, futuretr, junction_trace, ballgown, tracing/{mod,events,pipeline,reference},
parity/{mod,decisions,graph_edges_dump,junction_dump,partition_dump,shadow,trace_dump},
bin `src/bin/rustle.rs`.

> Tendril-support subset (`map_reads`, `bundle_builder`, `assembly_mode`, `junction_correction`,
> `parse_trflong_st`, `annotation_families`, `psv_fasta`): DEAD, but must keep compiling until the
> SHARED-tendril file that pulls them in (`vg` / `bundle` / `path_extract`) is extracted.

---

## 3. Removal ORDER (network-flow first)

1. **Delete the network-flow layer FIRST.**
   Remove `max_flow.rs`, `global_flow.rs`, `parity/flow_iter_dump.rs`. Their `pub mod` lines in
   `lib.rs` and every `use crate::max_flow`/`use crate::global_flow` in DEAD modules go with them.
   (No thesis code references either — the flow layer is the cleanest cut.)

2. **Extract the SHARED-tendril exports into a kept module** (see §4). Nothing else can be deleted
   until the thesis stops importing from `vg` / `graph` / `path_extract` / `bundle`.

3. **Relocate the two REVERSE tendrils** (KEEP → DEAD calls; see §4) so the KEEP modules stop
   referencing DEAD ones.

4. **Drop the legacy assembler bin + stack.** Remove `src/bin/rustle.rs`, `main.rs`,
   `pipeline.rs`, and the rest of DEAD-stringtie (bottom-up: leaves first, `pipeline` last), plus the
   crate-root `pub use` lines that re-export from them (`graph::*`, `gtf::write_gtf`,
   `path_extract::*`, `pipeline::run`). Delete the three DEAD-orphans anytime.

Re-run `cargo check` + the three thesis-bin builds after each step; each step is independently green.

---

## 4. SHARED-tendril EXTRACT-FIRST list (every assembler symbol the thesis imports)

> **⚠ 2026-07-14 CORRECTION (dep-trace verified, supersedes the 12-symbol list below).** Only **ONE**
> symbol is a genuine thesis *runtime* tendril: **`reverse_complement`** (used by `denovo_assemble`,
> `family_detect`, `family_rescue`, `denovo_pipeline`). **DONE** — relocated verbatim to
> `vg_family/seq_utils.rs` (commit `a14b988`), GSTM regression byte-identical. The other 11 symbols are
> **false tendrils**: `FamilyGroup`/`Graph`/`Transcript`/`discover_family_groups_layer2`/`graph_exon_seqs`/
> `exon_kmer_similarity`/`enumerate_diagnostic_sites`/`build_multimap_index_from_secondary_index`/
> `recompute_junction_stats`/`collect_secondary_index_from_bam`/`fnv1a64` are used **only** by
> legacy-reachable code — the `vg_family` modules `layer2`, `rescue`, `tandem`, `secondary_index`,
> `psv_linkage` (each imported by ZERO thesis-reachable module) and the legacy functions
> `family_graph::build_family_graph` / `build_family_graph_from_layer1_graphs` (callers: `pipeline`/`vg`/
> `layer2`/tests only). They are **deleted WITH the island**, not extracted.
>
> **Remaining carve = one atomic deletion** (gate each: 3 thesis bins build + full `--lib` tests +
> GSTM/PCDHB/MAGEA/DAZ byte-identical): (a) sever 2 reverse tendrils inline while deleting their source —
> `bam`→`polya::detect_polya_aligned_unaligned` (+`anchored_run_meets`,`poly_window_meets`) and
> `types`→`stringtie_parity::stringtie_exact`; (b) prune `family_graph`'s `build_family_graph*` (+
> `make_layer1_graph` test helper) so nothing references `graph::Graph`/`vg::FamilyGroup`; (c) drop the 5
> legacy `vg_family` modules from `mod.rs`; (d) remove every assembler `pub mod` in `lib.rs` (all but
> `util`/`types`/`bam`/`genome`/`vg_family`) + the crate-root re-exports (`graph::*`, `gtf::write_gtf`,
> `path_extract::*`, `pipeline::run`); (e) delete `src/bin/rustle.rs`; (f) delete the files. DEAD-orphans
> already removed (commit `84f0a82`).

The thesis imports exactly these symbols from to-be-deleted modules. Extract each into a kept module
(new `vg_family/` submodule, or a small `vg_family/legacy_support.rs`) before deleting the source.

**From `vg.rs`** (heaviest — also calls back into `vg_family` ~47×, so split the *family-discovery
core* out and mark the EM-reweighting / `-G2` assembler body dead):
1. `reverse_complement` (fn) — vg_family/{denovo_assemble,family_detect,family_rescue}.rs
2. `FamilyGroup` (struct) — vg_family/{family_graph,layer2,rescue}.rs
3. `fnv1a64` (fn)
4. `graph_exon_seqs` (fn)
5. `exon_kmer_similarity` (fn)
6. `enumerate_diagnostic_sites` (fn)
7. `discover_family_groups_layer2` (fn)
8. `build_multimap_index_from_secondary_index` (fn) — vg_family/secondary_index.rs

**From `graph.rs`:**
9. `Graph` (struct) + `Graph::new` — vg_family/{family_graph,layer2}.rs. Deps: util, types only (clean).

**From `path_extract.rs`** (551 KB; thesis needs ONE type):
10. `Transcript` (struct) — vg_family/layer2.rs

**From `bundle.rs`** (101 KB; thesis needs TWO fns):
11. `recompute_junction_stats` (fn)
12. `collect_secondary_index_from_bam` (fn) — re-exported at vg_family/secondary_index.rs:279

**REVERSE tendrils (KEEP → DEAD; relocate the small fn before deleting the DEAD module):**
- A. `bam.rs` → `polya::detect_polya_aligned_unaligned` (1 call). Relocate before deleting `polya.rs`.
- B. `types.rs` → `stringtie_parity::stringtie_exact()` (2 calls, an env-flag reader). Inline/relocate
  before deleting `stringtie_parity.rs`.
- C. (DEAD → DEAD, non-blocking) `junctions.rs` → `junction_correction` + `stringtie_parity`; resolves
  itself when `junctions` is removed.

Symbols imported from KEEP modules (`types::{Bundle, BundleRead, DetHashMap, DetHashSet, RunConfig,
Junction, JunctionStats, JunctionStat}`, `bam::{open_bam, exons_from_cigar, record_to_bundle_read,
extract_mismatches_vs_fasta}`, `genome::GenomeIndex`) are NOT extract-first — those modules stay.

---

## 5. Python → Rust migration plan

The shipped **O1 family DEFINITION** is currently **Python-only** (the Rust `vg_family::denovo_pipeline`
/ `family_detect` / `family_split` is an *earlier* read-conflict `E_c` detector; the shipped
γ-quasi-clique **homology** definition lives in `bench/*.py`). O2/O3/O4 already have shipped Rust.
Migrate bottom-up the import DAG so each Rust module compiles against already-ported deps.

### 5.1 MIGRATED so far (TIER 1 root, item #1) — `genome_family_def.py` → `vg_family/family_definition.rs`

The **root of the DAG** is done. Two commits:

- **`distinct_loci` + `LOCUS_OVERLAP`** (union-find ≥2-loci predicate): byte-parity 867/867 against a
  golden fixture from the shipped Python build (`#[test] distinct_loci_parity`; earlier commit e6771e1).
- **`refine_families` + refinement operator R** (`_refine_component` / `_induced_density` / `_split_once`
  / `biotype_purity` / `GAMMA=0.20` / `PURITY=0.50`): ported **option B — DETERMINISTIC, networkx-free,
  RNG-free splitter** (below).

**Option B — deterministic splitter (replaces the seeded networkx Louvain `_split_once`).**
`split_once` is a 3-stage RNG-free, networkx-free splitter, each stage guaranteeing ≥2 non-empty parts
for n>2 so `refine_component`'s density-gated recursion provably terminates with every leaf a
γ-quasi-clique (or |C|≤2):
1. **connected components** of the induced subgraph (deterministic BFS in ascending node order) —
   matches Python's `nx.connected_components` branch exactly;
2. else **deterministic greedy-modularity agglomeration** (Clauset-Newman-Moore): merge the community
   pair with greatest positive gain `dQ = k_xy/m − D_x·D_y/(2m²)`, community id = min node, `dQ` ties
   break to the lexicographically-smallest community pair (sorted keys) — natural communities, preferred
   over crude halving where it applies (this is what actually fires on the real blob);
3. else **sanctioned deterministic halving** by sorted node order (the doc's guaranteed-progress last
   resort; never fires on the real sparse graph).

**KEY MEASURED FACT reproduced.** The reconstructed shipped RNA graph (default gates: core≥0.19 ∧
aln≥0.24 ∧ ¬repeat-hub) has **1933 components, exactly 1 invoking the splitter** (a 213-node blob,
749 edges, density 0.033 < γ). The other 1932 return WHOLE via the density gate → the deterministic
core ports **byte-identically for 1932/1933**.

**Verification (`#[test]`s in `family_definition.rs`, against golden fixture
`vg_family/testdata/refine_families_fixture.json`, generated by `bench/gen_refine_families_fixture.py`,
md5 `05c7fab72931b65b93a7bd822925a2bf`):**
- (1) **SET byte-parity EXACT** on all 1932 clean (non-splitter) components — Rust refined family-set ==
  Python's for every one (584 clean families total);
- (2) **blob certificate VALID** — every Rust family from the 213-node blob is a γ-quasi-clique
  (induced density ≥ 0.20 or size ≤ 2) with ≥2 distinct loci. **Rust 25 vs Python 26** families = a
  different but equally-valid witness of the NP-hard max-γ-quasi-clique partition (the doc's stated
  "non-unique witness"; certificate PROPERTY, not count, is what the port guarantees);
- (3) **determinism** — full `refine_families` over the reconstructed graph reruns byte-identical;
- (4) **networkx-free / RNG-free** — no `rand`/`rng`/`networkx`/`louvain`/`shuffle` in the port code.

**Biotype note.** The RNA `genes` dict does carry `biotype` (real biotype when the projection is kept,
`"NA"` when floored); it is non-constant, so `biotype_purity` and the 3-term sort
`(purity<0.50, −len, name)` are ported faithfully. RNA's downstream `family_rna_refine.write_outputs`
re-sorts the catalog by member tuple, so this internal sort is only observable in the
`genome_family_def` standalone usage — hence the SET-based parity test for the RNA fixture.

**Remaining TIER 1:** #2 loaders, #3 refine-sweep driver, #4 `E_r` edge oracle (crux), #5 driver.

### TIER 1 — thesis-critical, attains O1, shipped/default, NO Rust equivalent
| # | Python (bench/) | KB | Role | Target Rust module | Tractability |
|---|---|---|---|---|---|
| 1 | `genome_family_def.py` | 32 | γ-quasi-clique `refine_families` / `distinct_loci` / GAMMA / SEED — **root of the DAG** | ✅ **MIGRATED** `vg_family/family_definition.rs` | done (option B, deterministic splitter — see §5.1) |
| 2 | `family_er_pr.py` | 36 | loaders + gene projection (`gene_of_dn`); the shipped metric | `vg_family/family_definition/loaders.rs` | straightforward |
| 3 | `graph_def_refine_sweep.py` | 40 | `_refine_component` + modularity-split sweep | `vg_family/family_definition/refine.rs` | medium (networkx `louvain_communities` / greedy-modularity → port) |
| 4 | `rna_only_edge_oracle.py` | 77 | the `E_r` RNA-homology edge oracle | `vg_family/family_definition/edge_oracle.rs` | **hardest** (scipy.stats + POA-substring); highest value |
| 5 | `family_rna_refine.py` | 51 | DRIVER wiring 1–4 + the two gates | extend `gw_family_catalog` bin (or new `family_define` bin) | do LAST (depends on 1–8) |

### TIER 2 — O1 gates + O1 copy-number, shipped/default context
| # | Python (bench/) | Role | Target Rust module | Notes |
|---|---|---|---|---|
| 6 | `recombinant_split.py` + `recombination_bridge_detector.py` | Alu-bridge over-merge split gate | `vg_family/family_gates.rs` (or extend `read_conflict.rs`) | medium |
| 7 | `multi_repeat_bridge_gate.py` + `vg_repeat_catalog.py` | repeat-bridge gate | `vg_family/family_gates.rs` | heavy |
| 8 | `family_copy_number.py` + `sun_identifiability.py` | O1 χ(H) copy-number cardinality + SUN ladder | `vg_family/copy_number.rs` | SUN logic partly overlaps existing `psv_linkage` |
| 9 | `recombinant_abstain.py` | **O2** abstain gate | fold into `vg_family/copy_assign.rs` | small, standalone (no local imports) — **easy quick win** |

### LEAVE IN PYTHON — already have shipped Rust (Python is an analysis mirror)
- `copy_assign.py` (O2) → `src/bin/copy_assign.rs` + `vg_family/copy_assign*`
- `allele_specific_junctions.py` (O3) → `src/bin/asj.rs` + `vg_family/allele_specific_junctions.rs`
- `absent_copy_*.py`, `reference_absent_catalog.py` (O4) → `vg_family/absent_copy.rs`

### LEAVE IN PYTHON — analysis-only / one-off (bulk of the ~335 scripts)
All `*_fig.py`, `make_*_slides.py`, `*_validate.py`, `*_check.py`, `sim_*`, `a1_/a2_/a4_*` overlays,
`soto_*`, `diploid_cn_*`, `dna_*`, `protein_family_def.py` (E_p formal-object companion, not the default
RNA catalog), theory checkers (`copy_assignment_theory_checks.py`, `bridge_theorem_check.py`,
`mwca_approximation_check.py`). These produce findings/figures, not the shipped pipeline.

### Migration ORDER
`9` (trivial, O2 quick win) → ~~`1`~~ ✅ → `2` → `3` → `4` → `5` → `6`/`7` → `8`.
Rationale: bottom-up the import DAG so each Rust module compiles against ported deps;
`rna_only_edge_oracle` (#4) is the tractability crux (POA + stats) and gates the driver (#5).
**Item #1 (`genome_family_def.py` root: `distinct_loci` + `refine_families`/R) is MIGRATED** — see §5.1.

---

## 6. Intra-thesis duplication (2026-07-13 simplification audit)

The assembler carve (§2–§4) is the big lever, but a 46-agent audit of `vg_family/*` itself found
duplication *inside* the thesis code. This is a separate axis from §5 (Python→Rust) and §3 (assembler
removal): here both sides are already LIVE Rust, so the work is unification, not porting.

### 6.1 Completed this pass (committed, `cargo build` + 94 module tests green)
- **`--vg-realign` doc truthfulness.** The `DenovoConfig::vg_realign` field doc and the `copy_assign`
  CLI help both claimed "REPORT-ONLY … does not alter the copy set / assignment / any emitted field".
  False: the flag runs `apply_realign_patch` (corrects assignments), `admit_novel_pools` (widens the
  copy set), and `recompute_realign_abundance` (recomputes EM abundance) — the intended Task-3
  end-to-end behavior. Docs corrected; **no behavior change** (OFF is still byte-identical).
- **Dead code removed:** `copies_overlap` (superseded by `catalog_overlaps`) + its test;
  `detect_families_from_files` (0 callers); the `RUSTLE_PSV_BAND_DEBUG` / `RUSTLE_PSV_COST_AUDIT`
  eprintln diagnostic blocks (87 lines, no callers — functional POASTA/BAND/MINIMAP2 escapes kept);
  `util::constants::KMER` (unreferenced).
- **Root scratch pruned:** 8 StringTie-parity-era `PRECISION_*`/`SESSION_SUMMARY`/`GUIDED_MODE` files
  (kept the LIVE `bench/PRECISION_RECALL_FRONTIER.*` O1 artifacts referenced by `family_define.rs`).

### 6.2 Duplicate logic to unify (all verified LIVE on both sides — medium effort)
| # | Duplication | Where | Note |
|---|---|---|---|
| a | **Two soft-EM abundance engines** | `copy_assign_pipeline::soft_quantify_em` (PSV-only, fixed ε=0.01, 100 iters, no convergence break) vs factored `em_copy_assign` (e_step/m_step/loglik) | both write the same `copy_abundance` column; unify on the factored engine |
| b | **Two γ-quasi-clique refiners** | `family_definition::refine_component` (family_define path) vs `family_split::gamma_quasi_clique_partition` (gw path) | identical density-gated skeleton, each with a *private* `induced_density` + `split_once` |
| c | **Two "≥2 distinct loci" operators, DIVERGENT policy** | `family_definition::distinct_loci` (reciprocal ≥50%) vs `denovo_pipeline::distinct_locus_reps` (any same-chrom overlap collapses) | two binaries → two definitions of "distinct locus"; a *semantic* inconsistency — pick one |
| d | **Two `E_r` edge-oracle sourcing paths** | `driver.rs` LOADS python-precomputed TSVs (`edge_oracle.rs`) vs `family_detect`/`denovo_pipeline` COMPUTE in-engine | load path reachable only via `family_define`; retires when `family_define` folds into `gw` |
| e | **~11 hand-rolled union-finds** | 3 of them (`uf_find`/`uf_union` in `family_split`, `family_detect`, `denovo_pipeline`) are byte-identical | hoist the identical trio to `util`; leave the struct-embedded ones |

Items **b/c/d** are the same story as §5: there are two O1-definition lineages — the `E_c` read-conflict
detector (`denovo_pipeline`/`family_detect`/`family_split`) and the γ-quasi-clique homology definition
(`family_definition`, migrated from `genome_family_def.py`). The unification and the Python→Rust
migration should converge on **one** O1 path (the `gw_family_catalog` genome-wide catalog).

### 6.3 Keep-or-decide (not dead, not obviously wrong)
- **`collapse_gate.rs`** — built-then-REFUTED (fires χ(H)=7 on single-copy EEF1A1), wired but
  permanently default-OFF, 10 unit tests. Keep as a documented negative result, or remove.
- **`o2_materialize` + `o2_columns` + `o2_margin_gate`** (2549 lines) — parity-tested O2 deliverable
  wired to **no** binary. Decide: expose it in a binary, or shelve it.
