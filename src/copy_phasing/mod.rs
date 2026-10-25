//! Copy phasing (2026-10-24): `bench/o3_graph/o3port.py run` in Rust, byte for byte (plan docs/superpowers/plans/2026-10-24-o3port-rust.md).
//! The Python is the specification; every function here carries its Python name and keeps its iteration and tie-break order.
//! Python dicts iterate in insertion order, so the maps that reach a product or a tie-break are `IndexMap`s; read names are
//! interned as `ReadId`s in sorted-name order, so comparing ids compares names (`sorted(reads)` = sorted ids).

pub mod cluster;
pub mod family;
pub mod pile;
pub mod putative;
pub mod pyjson;
pub mod slots;
pub mod stats;
pub mod tied_reads;

use indexmap::IndexMap;
use std::collections::{BTreeMap, BTreeSet};

pub const MIN_READS: usize = 2; // the project's ">= 2 independent observations"
pub const MIN_Q: u8 = 10; // base quality below which a base counts as missing
pub const AS_TIE: i32 = 10; // O2's tie rule (Eichler): alignments within 10 AS of the read's best
pub const ALPHA: f64 = 0.01; // IsoCon-style test, Bonferroni over the tested columns
pub const MAX_ROUNDS: usize = 100; // label propagation stops earlier when no label changes
pub const CLIP_BP: u32 = 50; // a soft clip this long that is not a poly-A / poly-T tail
pub const CLIP_BIN: i64 = 20; // clip boundaries in the same CLIP_BIN-bp bin on the same side count as one boundary
pub const MIN_SUPPORT: usize = 6; // reads a haplotype needs to take a slot or be certified
pub const SITE_FLANK: i64 = 5; // a distinguishing site needs this many gap-free aligned bases on each side in both alignments
/// Canonical splice-site dinucleotide pairs (donor, acceptor), both strands.
pub const CANONICAL: [(&str, &str); 6] = [("GT", "AG"), ("GC", "AG"), ("AT", "AC"), ("CT", "AC"), ("CT", "GC"), ("GT", "AT")];

/// A read, as an index into the run's sorted table of read names.
pub type ReadId = u32;

/// One observation at a column: Python's `(read, base, qual[, error rate])`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObsX {
    pub read: ReadId,
    pub base: u8,
    pub qual: u8,
    pub eps: Option<f64>,
}

/// `{column: [(read, base, qual[, eps])]}` in insertion order.
pub type Obs = IndexMap<i64, Vec<ObsX>>;
/// Informative columns: `{column: set of supported alleles}` in insertion order.
pub type Cols = IndexMap<i64, BTreeSet<u8>>;
/// A read's alleles at informative columns (`patterns`); only ever read sorted or by key.
pub type Pattern = BTreeMap<i64, u8>;
/// A haplotype: `{column: allele}` in the order of the columns it was built on.
pub type Hap = IndexMap<i64, u8>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Clusterer {
    Cc,
    Lp,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SlotOrder {
    Divergence,
    Support,
}

/// The run's switches (the Python's module globals and `run_family` keywords); `Default` = the Python defaults.
#[derive(Clone, Debug)]
pub struct Params {
    pub clipping: bool,
    pub anchor: bool,
    pub one_site: bool,
    pub guard: bool,
    pub prefer: bool,
    pub reconstruct: bool,
    pub semi: bool,
    pub coverage: bool,
    pub correct: bool,
    pub clusterer: Clusterer,
    pub slot_order: SlotOrder,
    pub structure: bool,
    pub floor_first: bool,
    pub absorb: bool,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            clipping: true,
            anchor: false,
            one_site: false,
            guard: false,
            prefer: false,
            reconstruct: true,
            semi: false,
            coverage: false,
            correct: false,
            clusterer: Clusterer::Cc,
            slot_order: SlotOrder::Divergence,
            structure: true,
            floor_first: true,
            absorb: false,
        }
    }
}

/// A copy of `--copies`: `(tid, chrom, first exon start, last exon end)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Copy {
    pub tid: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
}
