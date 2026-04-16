//! PolyA/T tail detection (reference-style).
//! Used for strand and boundary evidence.
//! shorten_terminal_exon: shortenFirstExon/shortenLastExon (trim first/last exon by max_trim bp).

/// Check if a sequence window is mostly A or mostly T (min_frac of length).
/// This helper is intentionally generic; BAM processing uses the base-specific
/// boundary-anchored checks below.
pub fn check_polya_window(seq: &[u8], min_frac: f64) -> bool {
    if seq.is_empty() {
        return false;
    }
    let len = seq.len() as f64;
    let a_count = seq.iter().filter(|&&c| c == b'A' || c == b'a').count() as f64;
    let t_count = seq.iter().filter(|&&c| c == b'T' || c == b't').count() as f64;
    (a_count >= len * min_frac) || (t_count >= len * min_frac)
}

/// CIGAR op: (kind, length). Kind: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P (SAM).
const CIGAR_SOFT_CLIP: u8 = 4;

#[inline]
fn is_poly_base(c: u8, good: u8) -> bool {
    c == good || c == good.to_ascii_lowercase()
}

// poly_window_meets(): progressive threshold over prefixes/suffixes,
// with early success/failure relative to the boundary being tested.
fn poly_window_meets(
    seq: &[u8],
    a: usize,
    b: usize,
    good: u8,
    min_count: usize,
    min_frac: f64,
    from_right: bool,
) -> bool {
    if b <= a || b > seq.len() {
        return false;
    }
    let n = b - a;
    let mut matches = 0usize;
    let mut checked = 0usize;

    if !from_right {
        for &base in &seq[a..b] {
            if is_poly_base(base, good) {
                matches += 1;
            }
            checked += 1;
            let req = min_count.max((min_frac * checked as f64).ceil() as usize);
            if matches >= req {
                return true;
            }
            let remaining = n - checked;
            if matches + remaining < min_count {
                return false;
            }
        }
    } else {
        for &base in seq[a..b].iter().rev() {
            if is_poly_base(base, good) {
                matches += 1;
            }
            checked += 1;
            let req = min_count.max((min_frac * checked as f64).ceil() as usize);
            if matches >= req {
                return true;
            }
            let remaining = n - checked;
            if matches + remaining < min_count {
                return false;
            }
        }
    }

    false
}

fn anchored_run_meets(
    seq: &[u8],
    a: usize,
    b: usize,
    good: u8,
    min_consec: usize,
    from_right: bool,
) -> bool {
    if b <= a || b > seq.len() {
        return false;
    }
    let mut run = 0usize;
    if !from_right {
        for &base in &seq[a..b] {
            if is_poly_base(base, good) {
                run += 1;
            } else {
                break;
            }
        }
    } else {
        for &base in seq[a..b].iter().rev() {
            if is_poly_base(base, good) {
                run += 1;
            } else {
                break;
            }
        }
    }
    run >= min_consec
}

/// Detect polyA/T at 5' and 3' of read from sequence + CIGAR.
/// - 5' soft-clip and/or first window bp
/// - 3' soft-clip and/or last window bp
/// Returns (has_poly_at_5prime, has_poly_at_3prime).
pub fn detect_polya_tail(
    seq: &[u8],
    cigar_first: Option<(u8, u32)>,
    cigar_last: Option<(u8, u32)>,
    min_consec: u32,
    window: usize,
    min_frac: f64,
) -> (bool, bool) {
    let mut has_poly_start = false;
    let mut has_poly_end = false;

    if seq.is_empty() {
        return (false, false);
    }
    let min_count = min_consec as usize;

    // 5' soft-clip
    if let Some((op, len)) = cigar_first {
        if op == CIGAR_SOFT_CLIP && len > 0 {
            let clip_len = (len as usize).min(seq.len());
            let a = clip_len.saturating_sub(window);
            let b = clip_len;
            if anchored_run_meets(seq, a, b, b'T', min_count, true)
                || poly_window_meets(seq, a, b, b'T', min_count, min_frac, true)
            {
                has_poly_start = true;
            }
        }
    }
    if !has_poly_start {
        let a = 0usize;
        let b = window.min(seq.len());
        if anchored_run_meets(seq, a, b, b'T', min_count, false)
            || poly_window_meets(seq, a, b, b'T', min_count, min_frac, false)
        {
            has_poly_start = true;
        }
    }

    // 3' soft-clip
    if let Some((op, len)) = cigar_last {
        if op == CIGAR_SOFT_CLIP && len > 0 {
            let clip_len = (len as usize).min(seq.len());
            let a = seq.len().saturating_sub(clip_len);
            let b = (a + window).min(seq.len());
            if anchored_run_meets(seq, a, b, b'A', min_count, false)
                || poly_window_meets(seq, a, b, b'A', min_count, min_frac, false)
            {
                has_poly_end = true;
            }
        }
    }
    if !has_poly_end {
        let b = seq.len();
        let a = b.saturating_sub(window);
        if anchored_run_meets(seq, a, b, b'A', min_count, true)
            || poly_window_meets(seq, a, b, b'A', min_count, min_frac, true)
        {
            has_poly_end = true;
        }
    }

    (has_poly_start, has_poly_end)
}

/// Separate aligned (in-alignment body) vs unaligned (soft-clip tail) polyA/T detection.
/// ref: aligned_polyA / unaligned_polyA / aligned_polyT / unaligned_polyT.
///
/// Returns `(has_aligned_start, has_unaligned_start, has_aligned_end, has_unaligned_end)`.
///
/// - `has_unaligned_start`: polyT found in the 5' soft-clip region.
/// - `has_aligned_start`: polyT found in the first `window` bases of the ALIGNED portion
///   (after 5' soft-clip).
/// - `has_unaligned_end`: polyA found in the 3' soft-clip region.
/// - `has_aligned_end`: polyA found in the last `window` bases of the ALIGNED portion
///   (before 3' soft-clip).
///
/// The condition `aligned && !unaligned` identifies an RT drop-off artifact: the polyA
/// run is inside the mapped alignment, not in the soft-clipped overhang (which would
/// indicate a real polyadenylation site). Used by shortenFirstExon/shortenLastExon.
pub fn detect_polya_aligned_unaligned(
    seq: &[u8],
    clip_left: u32,
    clip_right: u32,
    min_consec: u32,
    window: usize,
    min_frac: f64,
) -> (bool, bool, bool, bool) {
    let slen = seq.len();
    let clip_l = (clip_left as usize).min(slen);
    let clip_r = (clip_right as usize).min(slen);
    let min_count = min_consec as usize;

    // Unaligned 5': boundary-anchored polyT in the left soft clip.
    let has_unaligned_start = if clip_l > 0 {
        let a = clip_l.saturating_sub(window);
        let b = clip_l;
        anchored_run_meets(seq, a, b, b'T', min_count, true)
            || poly_window_meets(seq, a, b, b'T', min_count, min_frac, true)
    } else {
        false
    };

    // Aligned 5': first window of aligned sequence after left soft-clip.
    let astart_5 = clip_l;
    let aend_5 = (clip_l + window).min(slen.saturating_sub(clip_r));
    let has_aligned_start = anchored_run_meets(seq, astart_5, aend_5, b'T', min_count, false)
        || poly_window_meets(seq, astart_5, aend_5, b'T', min_count, min_frac, false);

    // Unaligned 3': boundary-anchored polyA in the right soft clip.
    let has_unaligned_end = if clip_r > 0 {
        let a = slen.saturating_sub(clip_r);
        let b = (a + window).min(slen);
        anchored_run_meets(seq, a, b, b'A', min_count, false)
            || poly_window_meets(seq, a, b, b'A', min_count, min_frac, false)
    } else {
        false
    };

    // Aligned 3': last window of aligned sequence before right soft-clip.
    let aend_3 = slen.saturating_sub(clip_r);
    let astart_3 = aend_3.saturating_sub(window).max(clip_l);
    let has_aligned_end = anchored_run_meets(seq, astart_3, aend_3, b'A', min_count, true)
        || poly_window_meets(seq, astart_3, aend_3, b'A', min_count, min_frac, true);

    (
        has_aligned_start,
        has_unaligned_start,
        has_aligned_end,
        has_unaligned_end,
    )
}

/// Shorten first or last exon by up to max_trim bp (shortenFirstExon/shortenLastExon).
/// Keeps at least 1 bp per exon. Exons are (start, end); first = trim from start, last = trim from end.
pub fn shorten_terminal_exon(read: &mut crate::types::BundleRead, first: bool, max_trim: u64) {
    if read.exons.is_empty() || max_trim == 0 {
        return;
    }
    if first {
        let (s, e) = read.exons[0];
        let len = e.saturating_sub(s);
        let trim = max_trim.min(len.saturating_sub(1));
        if trim > 0 {
            read.exons[0].0 = s + trim;
        }
    } else {
        let i = read.exons.len() - 1;
        let (s, e) = read.exons[i];
        let len = e.saturating_sub(s);
        let trim = max_trim.min(len.saturating_sub(1));
        if trim > 0 {
            read.exons[i].1 = e - trim;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_polya_window() {
        assert!(check_polya_window(b"AAAAA", 0.8));
        assert!(check_polya_window(b"TTTTT", 0.8));
        assert!(!check_polya_window(b"ACGTACGT", 0.8));
    }

    #[test]
    fn test_detect_polya_aligned_unaligned_boundary_run() {
        let seq = b"TTTTTACCCCCCCCCCCCC";
        let got = detect_polya_aligned_unaligned(seq, 0, 0, 5, 20, 0.8);
        assert_eq!(got, (true, false, false, false));
    }

    #[test]
    fn test_detect_polya_aligned_unaligned_base_specific() {
        let seq = b"AAAAACCCCCCCCCCCCCC";
        let got = detect_polya_aligned_unaligned(seq, 0, 0, 5, 20, 0.8);
        assert_eq!(got, (false, false, false, false));
    }

    #[test]
    fn test_detect_polya_aligned_unaligned_right_clip_boundary() {
        let seq = b"CCCCCAAAAAT";
        let got = detect_polya_aligned_unaligned(seq, 0, 6, 5, 20, 0.8);
        assert_eq!(got, (false, false, false, true));
    }
}
