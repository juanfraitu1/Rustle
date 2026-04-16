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

// ── Infrastructure ────────────────────────────────────────────────────────────
pub mod bitset; // SmallBitset: zero-alloc u64 bitset for ≤64 elements
pub mod bitvec; // GBitVec type alias → SmallBitset
pub mod constants; // numeric constants (BSIZE, KMER, FLOW_EPSILON, …)
pub mod coord;
pub mod types; // RunConfig, Bundle, all shared data types // coordinate / interval utilities

// ── Stage 1: Input ────────────────────────────────────────────────────────────
pub mod bam; // BAM file parsing and read ingestion
pub mod genome; // genome / chromosome metadata
pub mod gtf; // GTF reading and writing
pub mod reference_gtf; // reference annotation loading for guided assembly

// ── Stage 2: Bundle detection ─────────────────────────────────────────────────
pub mod bpcov; // base-pair coverage prefix sums
mod bundle; // Internal bundle implementation (legacy)
pub mod bundle_builder; // Sub-bundle building (public API)
pub mod hard_boundaries;
pub mod junction_correction; // sserror-based junction position correction
pub mod junctions; // splice junction data structures
pub mod killed_junctions; // good_junc gate, junction splice aggregation
pub mod read_boundaries; // read boundary inference // hard bundle boundary detection

// ── Stage 3: Graph construction ───────────────────────────────────────────────
pub mod graph; // GraphNode, Graph, adjacency, reachability
pub mod graph_build; // build graph from bundle nodes
pub mod nodecov; // node coverage accumulation (nodecov, longcov)

// ── Stage 4: Read mapping ─────────────────────────────────────────────────────
pub mod map_reads; // map reads/transfrags onto the graph
pub mod transfrag_process; // transfrag creation and processing

// ── Stage 5: Flow computation ─────────────────────────────────────────────────
pub mod max_flow; // max-flow seeded path capacity

// ── Stage 6: Path extraction ──────────────────────────────────────────────────
pub mod coverage_trim;
pub mod longtrim; // long-read boundary trimming
pub mod path_extract; // seeded forward/back path search, checktrf, EM
pub mod treepat; // tree pattern matching for path patterns // coverage-based transcript trimming

// ── Stage 7: Transcript filtering ────────────────────────────────────────────
pub mod exon_merge; // exon boundary merging
pub mod polya; // poly-A site detection
pub mod transcript_filter; // isofrac, pairwise, runoff, dedup, merge, lowintron
pub mod tss_tts; // TSS / TTS refinement

// ── Stage 8: Output ───────────────────────────────────────────────────────────
pub mod ballgown; // Ballgown table output
pub mod gene_abundance; // gene-level abundance summary
pub mod report_losses; // diagnostic loss reporting

// ── Stage 9: Orchestration ────────────────────────────────────────────────────
pub mod assembly_mode; // assembly mode dispatch
pub mod merge_mode; // merge mode (multi-sample)
pub mod debug_stage;
pub mod pipeline; // main long-read assembly pipeline

// ── Per-bundlenode graph processing ───────────────────────────────────────────
pub mod per_bnode_graph; // junction-expanded per-bundlenode graphs

// ── Variation graph (gene family) mode ───────────────────────────────────────
pub mod vg; // family group discovery, EM reweighting
pub mod vg_mflp; // MFLP (minimum flow linear program) solver

// ── JSON snapshots (debug) ────────────────────────────────────────────────────
pub mod snapshot; // per-bundle JSONL snapshots for 1:1 comparisons

// ── Tracing / diagnostics ─────────────────────────────────────────────────────
pub mod futuretr;
pub mod trace_events; // debug trace helpers
pub mod trace_reference; // per-reference transcript fate tracing // future transcript placeholders (diagnostic)

// ── Public re-exports ─────────────────────────────────────────────────────────
pub use graph::{CGraphnode, CMTransfrag, CPath, CTransfrag, Graph, GraphNode, GraphTransfrag};
pub use gtf::write_gtf;
pub use path_extract::{cpredictions_to_transcripts, transcripts_to_cpredictions};
pub use pipeline::run;
pub use types::{
    cjunctions_to_junction_stats, junction_stats_to_cjunctions, AssemblyMode, Bundle, BundleData,
    CBundlenode, CExon, CJunction, CMaxIntv, CPred, CPrediction, GArray, GEdge, GPVec, GVec,
    Junction, JunctionStats, RunConfig,
};
