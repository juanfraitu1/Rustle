//! Shared deterministic hash-container aliases.
//!
//! The StringTie-era data types that used to live here (`Bundle`, `BundleRead`, `Junction`,
//! `RunConfig`, the `C*` assembler structs) were removed 2026-09-24 with the last of their callers;
//! recover them from tag `notebook-2026-09-24`.

/// Fast deterministic hasher: FxHash is ~3-5x faster than SipHash for integer keys.
/// Deterministic because FxHash has no randomization seed.
pub type FixedBuild = fxhash::FxBuildHasher;
pub type DetHashMap<K, V> = std::collections::HashMap<K, V, FixedBuild>;
pub type DetHashSet<T> = std::collections::HashSet<T, FixedBuild>;
