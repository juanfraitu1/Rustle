/// Fast deterministic hasher: FxHash is ~3-5x faster than SipHash for integer keys.
/// Deterministic because FxHash has no randomization seed.
pub type FixedBuild = fxhash::FxBuildHasher;
pub type DetHashMap<K, V> = std::collections::HashMap<K, V, FixedBuild>;
pub type DetHashSet<T> = std::collections::HashSet<T, FixedBuild>;
