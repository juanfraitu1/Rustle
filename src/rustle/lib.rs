//! Rustle — Long-read transcript assembler (Rust port of the original algorithm pipeline)
//!
//! ## Pipeline stages (in execution order)
//!
//! 1. **Input** — BAM parsing, reference GTF, genome, coordinate utilities
//! 2. **Bundle detection** — read grouping, junction extraction, coverage, boundaries
//! 3. **Graph construction** — node/edge graph from bundle data, node coverage
//! 4. **Read mapping** — mapping reads/transfrags onto the graph
//! 5. **Flow computation** — global flow, max-flow path capacity
//! 6. **Path extraction** — seeded path search, EM, pattern matching, trimming
//! 7. **Transcript filtering** — isofrac, pairwise, runoff, dedup, merge, nascent
//! 8. **Output** — GTF writing, junction output, gene/ballgown abundance tables
//! 9. **Top-level orchestration** — pipeline, assembly modes

// ═══════════════════════════════════════════════════════════════════════════════
//  RETIREMENT NOTICE — see docs/RETIREMENT_AND_MIGRATION.md
//  ---------------------------------------------------------------------------
//  The multi-copy-family THESIS (O1 family-def, O2 copy-assign, O3 ASJ, O4
//  reference-absent copies) lives ONLY in `vg_family` + `genome` + the
//  foundational IO/type modules (`util`, `types`, `bam`).  EVERYTHING ELSE
//  declared below is LEGACY StringTie assembler + NETWORK FLOW — DEAD CODE,
//  slated for removal.  Do NOT extend the dead modules; build new thesis work
//  in `vg_family`.
//
//  Removal ORDER:
//    (1) NETWORK FLOW FIRST: `max_flow`, `global_flow` (+ parity/flow_iter_dump).
//    (2) EXTRACT the SHARED-tendril exports (`vg`, `graph`, `path_extract`,
//        `bundle`) — the only assembler symbols the thesis still imports — into
//        a kept module; then relocate the two REVERSE tendrils
//        (`polya::detect_polya_aligned_unaligned` used by bam.rs,
//        `stringtie_parity::stringtie_exact` used by types.rs).
//    (3) Drop the `rustle` bin (main.rs) + the rest of the assembler stack.
//
//  Inline tags below:  [KEEP] = thesis / foundational (stays);
//  [SHARED-tendril] = assembler file the thesis still imports from (extract 1st);
//  [DEAD] = StringTie assembler;  [DEAD-FLOW] = network flow (retire first).
// ═══════════════════════════════════════════════════════════════════════════════

// ── Infrastructure ────────────────────────────────────────────────────────────
pub mod util;  // [KEEP] bitset, bitvec, constants, coord, hard_counters (crate-wide low-level utilities)
pub mod types; // [KEEP] RunConfig, Bundle, all shared data types (REVERSE tendril -> stringtie_parity::stringtie_exact)

// ── Stage 1: Input ────────────────────────────────────────────────────────────
pub mod bam; // [KEEP] BAM parsing / read ingestion (REVERSE tendril -> polya::detect_polya_aligned_unaligned)
pub mod genome; // [KEEP] genome / chromosome metadata (thesis: GenomeIndex)
pub mod gtf; // [DEAD] GTF reading and writing
pub mod reference_gtf; // [DEAD] reference annotation loading for guided assembly

// ═══ LEGACY StringTie assembler + NETWORK FLOW — DEAD CODE, slated for removal ═══
// (Stages 2–9 below, plus the vg-mode legacy modules, are all assembler-era.
//  The only exceptions are [SHARED-tendril] files and `vg_family` [KEEP].)

// ── Stage 2: Bundle detection ─────────────────────────────────────────────────
pub mod bpcov; // [DEAD] base-pair coverage prefix sums
mod bundle; // [SHARED-tendril] extract recompute_junction_stats + collect_secondary_index_from_bam before removal
pub mod bundle_builder; // [DEAD] Sub-bundle building (public API)
pub mod hard_boundaries; // [DEAD]
pub mod junction_correction; // [DEAD] sserror-based junction position correction
pub mod junction_graph; // [DEAD] exon-endpoint clustering for junction normalization
pub mod junctions; // [DEAD] splice junction data structures (thesis uses types::Junction, not this)
pub mod junction_graph_st; // [DEAD] StringTie-style bundlenode graph inference
pub mod killed_junctions; // [DEAD] good_junc gate, junction splice aggregation
pub mod read_boundaries; // [DEAD] read boundary inference for bundles

// ── Stage 3: Graph construction ───────────────────────────────────────────────
pub mod graph; // [SHARED-tendril] extract Graph (+ Graph::new) before removal
pub mod graph_build; // [DEAD] build graph from bundle nodes
pub mod guide_path_fixer; // [DEAD] verify and reconstruct guide transcript paths
pub mod nodecov; // [DEAD] node coverage accumulation (nodecov, longcov)

// ── Stage 4: Read mapping ─────────────────────────────────────────────────────
pub mod map_reads; // [DEAD] map reads/transfrags onto the graph
pub mod transfrag_process; // [DEAD] transfrag creation and processing

// ── Stage 5: Flow computation ─────────────────────────────────────────────────
pub mod max_flow; // [DEAD-FLOW] max-flow seeded path capacity — RETIRE FIRST
pub mod global_flow; // [DEAD-FLOW] greedy flow decomposition (RUSTLE_GREEDY_DECOMPOSE=1) — RETIRE FIRST

// ── Stage 6: Path extraction ──────────────────────────────────────────────────
pub mod coverage_trim; // [DEAD]
pub mod longtrim; // [DEAD] long-read boundary trimming
pub mod path_extract; // [SHARED-tendril] extract the `Transcript` type before removal
pub mod treepat; // [DEAD] tree pattern matching for path patterns

// ── Stage 7: Transcript filtering ────────────────────────────────────────────
pub mod exon_merge; // [DEAD] exon boundary merging
pub mod polya; // [DEAD] poly-A site detection (REVERSE tendril target: relocate detect_polya_aligned_unaligned)
pub mod transcript_filter; // [DEAD] isofrac, pairwise, runoff, dedup, merge, lowintron
pub mod ml_filter;         // [DEAD] ML-based isofrac alternative (depth-3 decision tree)
pub mod cross_strand_predcluster; // [DEAD] post-pipeline cross-strand KRI (gated, default off)
pub mod single_exon_pileup; // [DEAD] polyA-driven SE detector (gated, opt-in)
pub mod parallel_predprune; // [DEAD] post-pipeline m/k-class pruning (gated, default off)
pub mod parse_trflong_st; // [DEAD] byte-faithful port of ST's parse_trflong + update_abundance (scaffold)
pub mod tss_tts; // [DEAD] TSS / TTS refinement

// ── Stage 8: Output ───────────────────────────────────────────────────────────
pub mod ballgown; // [DEAD] Ballgown table output
pub mod gene_abundance; // [DEAD] gene-level abundance summary
pub mod report_losses; // [DEAD] diagnostic loss reporting

// ── StringTie-parity scaffold ────────────────────────────────────────────────
pub mod stringtie_parity; // [DEAD] RUSTLE_STRINGTIE_EXACT meta-flag (REVERSE tendril target: relocate stringtie_exact())
pub mod terminal_oracle;  // [DEAD] RUSTLE_TERMINAL_ORACLE: inject ST terminal edges (analysis-only, default OFF)
pub mod predcluster_st;   // [DEAD] RUSTLE_PREDCLUSTER_ST: ST-faithful prediction selection (default OFF)
pub mod parity;           // [DEAD] parity dumps + decision log (parity/flow_iter_dump is [DEAD-FLOW])

// ── Stage 9: Orchestration ────────────────────────────────────────────────────
pub mod assembly_mode; // [DEAD] assembly mode dispatch
pub mod merge_mode; // [DEAD] merge mode (multi-sample)
pub mod debug_stage; // [DEAD]
pub mod vg_runtime; // [DEAD] process-global "--vg active" flag for the anchored-coverage channel
pub mod pipeline; // [DEAD] main long-read assembly pipeline (rustle::run entry point)

// ── Per-bundlenode graph processing ───────────────────────────────────────────
pub mod per_bnode_graph; // [DEAD] junction-expanded per-bundlenode graphs

// ── Variation graph (gene family) mode ───────────────────────────────────────
pub mod annotation_families; // [DEAD] legacy -G2 annotation -> paralog families (NOT the thesis O1 catalog)
pub mod vg_eval; // [DEAD] legacy -G2 reports
pub mod vg; // [SHARED-tendril] 8 symbols imported by thesis; split family-discovery core out before removal
pub mod vg_ingestion; // [DEAD] GTF/GFF ingestion mode for template-based family assembly
pub mod family_manifest; // [DEAD] TSV-manifest ingestion: R-annotated multi-locus families
pub mod vg_family; // [KEEP] THESIS LAYER: O1 family detect/split, O2 copy_assign, O3 ASJ, O4 absent-copy
pub mod psv_fasta; // [DEAD] opt-in PSV-FASTA emission (RUSTLE_VG_PSV_FASTA)
pub mod graph_comparison; // [DEAD] Compare de novo graphs against reference GTF for diagnosis

// ── JSON snapshots (debug) ────────────────────────────────────────────────────
pub mod snapshot; // [DEAD] per-bundle JSONL snapshots for 1:1 comparisons

// ── Tracing / diagnostics ─────────────────────────────────────────────────────
pub mod futuretr; // [DEAD]
pub mod tracing; // [DEAD] events, stage-level TSV dumps, per-reference transcript fate tracing
pub mod junction_trace; // [DEAD] Per-junction pipeline tracing for bottleneck diagnosis

// ── Public re-exports ─────────────────────────────────────────────────────────
// NOTE: these re-export from [SHARED-tendril]/[DEAD] modules; they disappear with
// those modules (thesis bins do NOT depend on them — they use vg_family + genome).
pub use graph::{CGraphnode, CMTransfrag, CPath, CTransfrag, Graph, GraphNode, GraphTransfrag};
pub use gtf::write_gtf;
pub use path_extract::{cpredictions_to_transcripts, transcripts_to_cpredictions};
pub use pipeline::run;
pub use types::{
    cjunctions_to_junction_stats, junction_stats_to_cjunctions, AssemblyMode, Bundle, BundleData,
    CBundlenode, CExon, CJunction, CMaxIntv, CPred, CPrediction, GArray, GEdge, GPVec, GVec,
    Junction, JunctionStats, RunConfig,
};
