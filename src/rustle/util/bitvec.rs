//! GBitVec — now a type alias over SmallBitset.
//! All call sites continue to compile unchanged.
pub use crate::util::bitset::SmallBitset as GBitVec;
