//! DEAD CODE -- StringTie-era assembler; slated for removal.
//! NOT part of the multi-copy-family thesis (O1 family-def / O2 copy-assign / O3 ASJ / O4 absent-copy).
//! See docs/RETIREMENT_AND_MIGRATION.md. Do not extend.
//! Pipeline tracing / diagnostic emission.
//!
//! All modules here are debug-only helpers for following pipeline behaviour:
//! - `events`    — generic trace event helpers used throughout the pipeline.
//! - `pipeline`  — stage-level TSV dumps for diff-based bisection.
//! - `reference` — per-reference transcript fate tracing.
//!
//! All are gated by env vars (`RUSTLE_TRACE_*`) and produce no output unless
//! enabled.

pub mod events;
pub mod pipeline;
pub mod reference;
