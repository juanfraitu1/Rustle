//! Centralized constants for ref/

/// Bundle block size used by prefix-sum coverage (`get_cov`).
pub const BSIZE: usize = 10_000;
/// K-mer length used in overlap adjustment logic.
pub const KMER: usize = 31;
/// graph node cap sentinel.
pub const MAX_NODE: i64 = 1_000_000;
/// floating sentinel used in merge DP.
pub const MIN_VAL: f64 = -100_000.0;
/// Long intron threshold header.
pub const LONGINTRON: u64 = 100_000;
/// Mismatch fraction threshold header.
pub const MISMATCHFRAC: f64 = 0.02;
/// Low coverage threshold header.
pub const LOWCOV: f64 = 1.5;
/// Splice-site error tolerance window from globals (common default).
pub const SSERROR: u64 = 25;
/// Global floating epsilon used by `0.000001`.
pub const FLOW_EPSILON: f64 = 1e-6;
