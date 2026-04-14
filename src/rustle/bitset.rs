//! Zero-allocation bitset: u64 inline for ≤ 64 bits, FixedBitSet fallback for larger.
//! Replaces GBitVec (bit-vec wrapper), Vec<bool> masks, HashSet<usize> for small integers.

use fixedbitset::FixedBitSet;

#[derive(Clone, Debug)]
pub enum SmallBitset {
    Inline(u64),
    Large(Box<FixedBitSet>),
}

impl Default for SmallBitset {
    fn default() -> Self {
        Self::Inline(0)
    }
}

impl PartialEq for SmallBitset {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Inline(a), Self::Inline(b)) => a == b,
            (Self::Large(a), Self::Large(b)) => a == b,
            (Self::Inline(w), Self::Large(bs)) | (Self::Large(bs), Self::Inline(w)) => {
                if *w == 0 {
                    return bs.is_clear();
                }
                let w_ones = w.count_ones() as usize;
                if bs.count_ones(..) != w_ones {
                    return false;
                }
                (0..64usize).all(|i| ((*w >> i) & 1 != 0) == bs.contains(i))
            }
        }
    }
}

impl Eq for SmallBitset {}

impl SmallBitset {
    /// Zero-cost empty bitset.
    #[inline]
    pub fn empty() -> Self {
        Self::Inline(0)
    }

    /// Allocate capacity for `n` bits (Inline if n ≤ 64, Large otherwise).
    pub fn with_capacity(n: usize) -> Self {
        if n <= 64 {
            Self::Inline(0)
        } else {
            Self::Large(Box::new(FixedBitSet::with_capacity(n)))
        }
    }

    /// Capacity of this bitset (64 for Inline, allocated size for Large).
    pub fn capacity(&self) -> usize {
        match self {
            Self::Inline(_) => 64,
            Self::Large(bs) => bs.len(),
        }
    }

    /// Insert bit `i`. For Inline, panics in debug if i >= 64.
    #[inline]
    pub fn insert(&mut self, i: usize) {
        match self {
            Self::Inline(w) => {
                debug_assert!(i < 64, "SmallBitset::Inline overflow: i={}", i);
                *w |= 1u64 << i;
            }
            Self::Large(bs) => {
                if i >= bs.len() {
                    bs.grow(i + 1);
                }
                bs.insert(i);
            }
        }
    }

    /// Insert bit `i`, growing to Large if needed (safe for any i).
    pub fn insert_grow(&mut self, i: usize) {
        match self {
            Self::Inline(w) => {
                if i < 64 {
                    *w |= 1u64 << i;
                } else {
                    let mut bs = FixedBitSet::with_capacity(i + 1);
                    for j in 0..64usize {
                        if (*w >> j) & 1 != 0 {
                            bs.insert(j);
                        }
                    }
                    bs.insert(i);
                    *self = Self::Large(Box::new(bs));
                }
            }
            Self::Large(bs) => {
                if i >= bs.len() {
                    bs.grow(i + 1);
                }
                bs.insert(i);
            }
        }
    }

    #[inline]
    pub fn contains(&self, i: usize) -> bool {
        match self {
            Self::Inline(w) => i < 64 && (*w >> i) & 1 != 0,
            Self::Large(bs) => i < bs.len() && bs.contains(i),
        }
    }

    #[inline]
    pub fn remove(&mut self, i: usize) {
        match self {
            Self::Inline(w) => {
                if i < 64 {
                    *w &= !(1u64 << i);
                }
            }
            Self::Large(bs) => {
                if i < bs.len() {
                    bs.set(i, false);
                }
            }
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        match self {
            Self::Inline(w) => *w == 0,
            Self::Large(bs) => bs.is_clear(),
        }
    }

    #[inline]
    pub fn count_ones(&self) -> usize {
        match self {
            Self::Inline(w) => w.count_ones() as usize,
            Self::Large(bs) => bs.count_ones(..),
        }
    }

    #[inline]
    pub fn clear(&mut self) {
        match self {
            Self::Inline(w) => *w = 0,
            Self::Large(bs) => bs.clear(),
        }
    }

    /// True iff `self` and `other` share at least one set bit.
    #[inline]
    pub fn intersects(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Inline(a), Self::Inline(b)) => a & b != 0,
            (Self::Large(a), Self::Large(b)) => a.ones().any(|i| b.contains(i)),
            (Self::Inline(w), Self::Large(bs)) | (Self::Large(bs), Self::Inline(w)) => {
                bs.ones().any(|i| i < 64 && (*w >> i) & 1 != 0)
            }
        }
    }

    /// self |= other
    pub fn union_with(&mut self, other: &Self) {
        match (self, other) {
            (Self::Inline(a), Self::Inline(b)) => *a |= b,
            (Self::Large(a), Self::Large(b)) => a.union_with(b),
            (Self::Large(a), Self::Inline(b)) => {
                for i in 0..64usize {
                    if (*b >> i) & 1 != 0 {
                        if i >= a.len() {
                            a.grow(i + 1);
                        }
                        a.insert(i);
                    }
                }
            }
            (this @ Self::Inline(_), Self::Large(b)) => {
                let mut bs = b.as_ref().clone();
                if let Self::Inline(w) = this {
                    for i in 0..64usize {
                        if (*w >> i) & 1 != 0 {
                            if i >= bs.len() {
                                bs.grow(i + 1);
                            }
                            bs.insert(i);
                        }
                    }
                }
                *this = Self::Large(Box::new(bs));
            }
        }
    }

    /// self &= other
    pub fn intersect_with(&mut self, other: &Self) {
        match (self, other) {
            (Self::Inline(a), Self::Inline(b)) => *a &= b,
            (Self::Large(a), Self::Large(b)) => a.intersect_with(b),
            (Self::Large(a), Self::Inline(b)) => {
                for i in 0..a.len() {
                    if i >= 64 || (*b >> i) & 1 == 0 {
                        a.set(i, false);
                    }
                }
            }
            (Self::Inline(a), Self::Large(b)) => {
                let mut mask = 0u64;
                for i in 0..64usize {
                    if (*a >> i) & 1 != 0 && b.contains(i) {
                        mask |= 1u64 << i;
                    }
                }
                *a = mask;
            }
        }
    }

    /// True iff every bit set in `other` is also set in `self` (other ⊆ self).
    pub fn contains_all(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Inline(a), Self::Inline(b)) => a & b == *b,
            (Self::Large(a), Self::Large(b)) => b.ones().all(|i| a.contains(i)),
            (Self::Inline(a), Self::Large(b)) => b.ones().all(|i| i < 64 && (*a >> i) & 1 != 0),
            (Self::Large(a), Self::Inline(b)) => {
                (0..64usize).all(|i| (*b >> i) & 1 == 0 || a.contains(i))
            }
        }
    }

    /// Iterator over indices of set bits, in ascending order.
    pub fn ones(&self) -> SmallBitsetOnes<'_> {
        match self {
            Self::Inline(w) => SmallBitsetOnes::Inline { word: *w },
            Self::Large(bs) => SmallBitsetOnes::Large { inner: bs.ones() },
        }
    }

    /// Build from iterator, growing to Large if needed.
    pub fn from_iter(iter: impl IntoIterator<Item = usize>) -> Self {
        let mut bs = Self::empty();
        for i in iter {
            bs.insert_grow(i);
        }
        bs
    }

    /// Extend with items from iterator.
    pub fn extend(&mut self, iter: impl IntoIterator<Item = usize>) {
        for i in iter {
            self.insert_grow(i);
        }
    }

    /// Create with a single bit set.
    pub fn singleton(i: usize) -> Self {
        let mut bs = Self::empty();
        bs.insert_grow(i);
        bs
    }

    // ── GBitVec compatibility API ─────────────────────────────────────────────

    /// Create with capacity `n` (GBitVec::new parity).
    pub fn new(n: usize) -> Self {
        Self::with_capacity(n)
    }

    /// Grow to at least `new_size` bits (GBitVec::grow parity).
    pub fn grow(&mut self, new_size: usize) {
        if new_size > self.capacity() {
            match self {
                Self::Inline(w) => {
                    if new_size > 64 {
                        let mut bs = FixedBitSet::with_capacity(new_size);
                        for i in 0..64usize {
                            if (*w >> i) & 1 != 0 {
                                bs.insert(i);
                            }
                        }
                        *self = Self::Large(Box::new(bs));
                    }
                    // else: Inline always has 64-bit capacity, nothing to do
                }
                Self::Large(bs) => bs.grow(new_size),
            }
        }
    }

    pub fn set_bit(&mut self, pos: usize) {
        self.insert_grow(pos);
    }
    pub fn clear_bit(&mut self, pos: usize) {
        self.remove(pos);
    }
    pub fn get_bit(&self, pos: usize) -> bool {
        self.contains(pos)
    }
    pub fn reset(&mut self) {
        self.clear();
    }
    pub fn len_bits(&self) -> usize {
        self.count_ones()
    }

    /// self |= other (GBitVec::or_assign parity).
    pub fn or_assign(&mut self, other: &Self) {
        self.union_with(other);
    }

    /// True iff every bit in `other` is set in `self` (GBitVec::contains_pattern parity).
    pub fn contains_pattern(&self, other: &Self) -> bool {
        self.contains_all(other)
    }

    /// True iff every bit in 0..n_nodes set in `other` is also set in `self`
    /// (GBitVec::contains_nodes_only parity).
    pub fn contains_nodes_only(&self, other: &Self, n_nodes: usize) -> bool {
        other
            .ones()
            .take_while(|&i| i < n_nodes)
            .all(|i| self.contains(i))
    }
}

// ── Iterator ─────────────────────────────────────────────────────────────────

pub enum SmallBitsetOnes<'a> {
    Inline { word: u64 },
    Large { inner: fixedbitset::Ones<'a> },
}

impl Iterator for SmallBitsetOnes<'_> {
    type Item = usize;
    fn next(&mut self) -> Option<usize> {
        match self {
            SmallBitsetOnes::Inline { word } => {
                if *word == 0 {
                    return None;
                }
                let tz = word.trailing_zeros() as usize;
                *word &= *word - 1;
                Some(tz)
            }
            SmallBitsetOnes::Large { inner } => inner.next(),
        }
    }
}

// ── BitAnd / BitOr (GBitVec operator parity) ─────────────────────────────────

impl std::ops::BitAnd for SmallBitset {
    type Output = SmallBitset;
    fn bitand(mut self, rhs: Self) -> Self {
        self.intersect_with(&rhs);
        self
    }
}

impl std::ops::BitOr for SmallBitset {
    type Output = SmallBitset;
    fn bitor(mut self, rhs: Self) -> Self {
        self.union_with(&rhs);
        self
    }
}

// ── Unit tests ────────────────────────────────────────────────────────────────

/// Type alias for node-tracking bitsets (replaces HashSet<usize>).
/// Uses inline storage for ≤64 nodes, heap only when needed.
pub type NodeSet = SmallBitset;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inline_insert_contains() {
        let mut b = SmallBitset::empty();
        b.insert(0);
        b.insert(7);
        b.insert(63);
        assert!(b.contains(0));
        assert!(b.contains(7));
        assert!(b.contains(63));
        assert!(!b.contains(1));
        assert!(!b.contains(62));
    }

    #[test]
    fn inline_remove() {
        let mut b = SmallBitset::empty();
        b.insert(5);
        b.remove(5);
        assert!(!b.contains(5));
        assert!(b.is_empty());
    }

    #[test]
    fn inline_intersects() {
        let mut a = SmallBitset::empty();
        a.insert(3);
        a.insert(7);
        let mut b = SmallBitset::empty();
        b.insert(7);
        b.insert(9);
        assert!(a.intersects(&b));
        let mut c = SmallBitset::empty();
        c.insert(1);
        assert!(!a.intersects(&c));
    }

    #[test]
    fn inline_union_with() {
        let mut a = SmallBitset::empty();
        a.insert(1);
        let mut b = SmallBitset::empty();
        b.insert(2);
        a.union_with(&b);
        assert!(a.contains(1));
        assert!(a.contains(2));
    }

    #[test]
    fn inline_count_ones() {
        let mut b = SmallBitset::empty();
        b.insert(0);
        b.insert(10);
        b.insert(20);
        assert_eq!(b.count_ones(), 3);
    }

    #[test]
    fn ones_iterator() {
        let mut b = SmallBitset::empty();
        b.insert(2);
        b.insert(5);
        b.insert(11);
        let v: Vec<_> = b.ones().collect();
        assert_eq!(v, vec![2, 5, 11]);
    }

    #[test]
    fn large_basic() {
        let mut b = SmallBitset::with_capacity(128);
        b.insert(65);
        b.insert(100);
        assert!(b.contains(65));
        assert!(b.contains(100));
        assert!(!b.contains(64));
    }

    #[test]
    fn grow_inline_to_large() {
        let mut b = SmallBitset::empty();
        b.insert(3);
        b.insert_grow(70);
        assert!(b.contains(3));
        assert!(b.contains(70));
    }

    #[test]
    fn contains_all() {
        let mut a = SmallBitset::empty();
        a.insert(1);
        a.insert(2);
        a.insert(3);
        let mut b = SmallBitset::empty();
        b.insert(1);
        b.insert(3);
        assert!(a.contains_all(&b));
        b.insert(5);
        assert!(!a.contains_all(&b));
    }

    #[test]
    fn gbvec_compat() {
        let mut b = SmallBitset::new(10);
        b.set_bit(4);
        b.set_bit(9);
        assert!(b.get_bit(4));
        assert!(b.get_bit(9));
        assert!(!b.get_bit(5));
        assert_eq!(b.len_bits(), 2);
        b.reset();
        assert!(b.is_empty());
    }

    #[test]
    fn from_iter_basic() {
        let b = SmallBitset::from_iter([1, 5, 10]);
        assert!(b.contains(1));
        assert!(b.contains(5));
        assert!(b.contains(10));
        assert!(!b.contains(0));
        assert!(!b.contains(2));
    }

    #[test]
    fn from_iter_large() {
        let b = SmallBitset::from_iter([0, 64, 128]);
        assert!(b.contains(0));
        assert!(b.contains(64));
        assert!(b.contains(128));
        assert!(!b.contains(63));
    }

    #[test]
    fn extend_basic() {
        let mut b = SmallBitset::empty();
        b.extend([1, 2, 3]);
        assert_eq!(b.count_ones(), 3);
        b.extend([4, 5]);
        assert_eq!(b.count_ones(), 5);
    }

    #[test]
    fn nodeset_alias() {
        let mut ns: NodeSet = NodeSet::empty();
        ns.insert(0);
        ns.insert(5);
        assert!(ns.contains(0));
        assert!(ns.contains(5));
    }
}
