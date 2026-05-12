//! Crate-wide low-level utilities — small data structures and numeric helpers
//! that don't belong to any particular pipeline stage.
//!
//! - `bitset`         — SmallBitset: zero-alloc u64 bitset for ≤ 64 elements.
//! - `bitvec`         — GBitVec type alias → SmallBitset.
//! - `constants`      — numeric constants (BSIZE, KMER, FLOW_EPSILON, …).
//! - `coord`          — coordinate / interval utilities.
//! - `hard_counters`  — global counters for diagnostic accounting.

pub mod bitset;
pub mod bitvec;
pub mod constants;
pub mod coord;
pub mod hard_counters;
