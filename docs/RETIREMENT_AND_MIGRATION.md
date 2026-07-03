# Rustle retirement + migration plan

Status: **MARKING PASS COMPLETE** (comments only — no logic changed, no files removed).
Branch: `vg/flow-capacity-apportionment`. Machine-readable map: `bench/RETIREMENT_MAP.tsv`.

Rustle began as a Rust port of the StringTie long-read transcript **assembler**. The thesis has
since moved to **multi-copy gene FAMILIES**:

- **O1** — family DEFINITION (γ-quasi-clique homology components; copy number χ(H))
- **O2** — copy ASSIGNMENT under MAPQ-0 ambiguity (PSV + junction + SUN, assign-or-abstain)
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

### TIER 1 — thesis-critical, attains O1, shipped/default, NO Rust equivalent
| # | Python (bench/) | KB | Role | Target Rust module | Tractability |
|---|---|---|---|---|---|
| 1 | `genome_family_def.py` | 32 | γ-quasi-clique `refine_families` / `distinct_loci` / GAMMA / SEED — **root of the DAG** | new `vg_family/family_definition.rs` | graph + Louvain; petgraph or hand-rolled, deterministic seed |
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
`9` (trivial, O2 quick win) → `1` → `2` → `3` → `4` → `5` → `6`/`7` → `8`.
Rationale: bottom-up the import DAG so each Rust module compiles against ported deps;
`rna_only_edge_oracle` (#4) is the tractability crux (POA + stats) and gates the driver (#5).
