//! Parity-with-StringTie scaffolding.
//!
//! All modules here exist to support cross-tool A/B diffs against StringTie
//! during the long-read parity work — they produce structured TSV/JSONL dumps
//! at various pipeline stages so behaviour can be diffed mechanically rather
//! than visually.  None of these are load-bearing for transcript output; they
//! are gated by env vars or `--debug-stage` flags and are off by default.
//!
//! See `project_parity_decisions_*.md` memos for the original motivation and
//! how each dump corresponds to a StringTie stage.

pub mod decisions;        // structured per-decision JSONL log for cross-tool diffing
pub mod shadow;           // layer-by-layer shadow parity logging
pub mod partition_dump;   // canonical partition geometry TSV
pub mod junction_dump;    // multi-stage junction-set TSV
pub mod graph_edges_dump; // splice-graph edges TSV (env-gated)
pub mod flow_iter_dump;   // per-iteration max-flow augmenting-path TSV
pub mod trace_dump;       // cgroup + graph topology + optional read-exon trace TSV
