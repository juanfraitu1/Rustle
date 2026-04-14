//! Centralized constants for C++ ref/C++ parity.

/// Bundle block size used by C++ prefix-sum coverage (`get_cov`).
pub const BSIZE: usize = 10_000;
/// K-mer length used in overlap adjustment logic.
pub const KMER: usize = 31;
/// C++ graph node cap sentinel.
pub const MAX_NODE: i64 = 1_000_000;
/// C++ floating sentinel used in merge DP.
pub const MIN_VAL: f64 = -100_000.0;
/// Long intron threshold from C++ header.
pub const LONGINTRON: u64 = 100_000;
/// Mismatch fraction threshold from C++ header.
pub const MISMATCHFRAC: f64 = 0.02;
/// Low coverage threshold from C++ header.
pub const LOWCOV: f64 = 1.5;
/// Splice-site error tolerance window from C++ reference globals (common default).
pub const SSERROR: u64 = 25;
/// Global floating epsilon used by C++ (C++ reference): `0.000001`.
pub const FLOW_EPSILON: f64 = 1e-6;
