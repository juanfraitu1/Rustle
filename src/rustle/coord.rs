//! Coordinate helpers for
//!
//! Rustle stores genomic intervals as 0-based half-open: `[start, end)`.
//! the original algorithm commonly reasons in inclusive end coordinates in local logic.
//! These helpers make conversions explicit and keep boundary math consistent.

/// Half-open interval length for `[start, end)`.
#[inline]
pub fn len_half_open(start: u64, end: u64) -> u64 {
    end.saturating_sub(start)
}

/// Half-open overlap predicate: true iff `[a0,a1)` intersects `[b0,b1)`.
#[inline]
pub fn overlaps_half_open(a0: u64, a1: u64, b0: u64, b1: u64) -> bool {
    a0 < b1 && b0 < a1
}

/// Half-open overlap length between `[a0,a1)` and `[b0,b1)`.
#[inline]
pub fn overlap_len_half_open(a0: u64, a1: u64, b0: u64, b1: u64) -> u64 {
    if !overlaps_half_open(a0, a1, b0, b1) {
        return 0;
    }
    a1.min(b1).saturating_sub(a0.max(b0))
}

/// True when a boundary is contiguous in half-open convention: `[a0,a1)` followed by `[a1,b1)`.
#[inline]
pub fn contiguous_half_open(left_end: u64, right_start: u64) -> bool {
    left_end == right_start
}

/// Convert inclusive interval `[start, end_inclusive]` to Rust half-open `[start, end_exclusive)`.
#[inline]
pub fn inclusive_to_half_open(start: u64, end_inclusive: u64) -> (u64, u64) {
    if end_inclusive < start {
        return (start, start);
    }
    (start, end_inclusive.saturating_add(1))
}

/// Convert Rust half-open `[start, end_exclusive)` to inclusive `[start, end_inclusive]`.
#[inline]
pub fn half_open_to_inclusive(start: u64, end_exclusive: u64) -> Option<(u64, u64)> {
    if end_exclusive <= start {
        return None;
    }
    Some((start, end_exclusive - 1))
}

/// inclusive interval length for `[start, end_inclusive]`.
#[inline]
pub fn len_inclusive(start: u64, end_inclusive: u64) -> u64 {
    if end_inclusive < start {
        return 0;
    }
    end_inclusive.saturating_sub(start).saturating_add(1)
}

/// Inclusive overlap length between `[a0,a1]` and `[b0,b1]` (convention).
#[inline]
pub fn overlap_len_inclusive(a0: u64, a1: u64, b0: u64, b1: u64) -> u64 {
    let s = a0.max(b0);
    let e = a1.min(b1);
    len_inclusive(s, e)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inclusive_half_open_length_compat() {
        let cases = [(0, 0), (10, 10), (10, 11), (10, 100), (1000, 1005)];
        for (s, einc) in cases {
            let (hs, he) = inclusive_to_half_open(s, einc);
            assert_eq!(len_inclusive(s, einc), len_half_open(hs, he));
        }
    }

    #[test]
    fn inclusive_half_open_overlap_compat() {
        // `[1,3]` and `[3,5]` overlap at base 3 in inclusive coordinates.
        let inc_ov = overlap_len_inclusive(1, 3, 3, 5);
        let a = inclusive_to_half_open(1, 3);
        let b = inclusive_to_half_open(3, 5);
        let ho_ov = overlap_len_half_open(a.0, a.1, b.0, b.1);
        assert_eq!(inc_ov, ho_ov);

        // Disjoint intervals.
        let inc_ov = overlap_len_inclusive(1, 2, 4, 7);
        let a = inclusive_to_half_open(1, 2);
        let b = inclusive_to_half_open(4, 7);
        let ho_ov = overlap_len_half_open(a.0, a.1, b.0, b.1);
        assert_eq!(inc_ov, ho_ov);
    }

    #[test]
    fn roundtrip_half_open_inclusive() {
        let iv = (25, 30);
        let inc = half_open_to_inclusive(iv.0, iv.1).unwrap();
        let back = inclusive_to_half_open(inc.0, inc.1);
        assert_eq!(iv, back);
    }
}
