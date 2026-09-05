//! Segmental-duplication extent / breakpoint calling from GENOME homology.
//!
//! A bare gene paralog shares homology only across the gene (exons/introns). A true
//! SEGMENTAL DUPLICATION copied the gene PLUS flanking sequence, so cross-copy homology
//! extends into the intergenic flanks; the duplication BREAKPOINTS are where that homology
//! ends. IsoSeq reads are spliced mRNA and do not cover intergenic flanks, so this signal
//! comes from the GENOME sequence (which the VG already loads for family discovery), anchored
//! at the gene the family discovered. The core is a pure anchored sliding-window homology scan.
//!
//! **STATUS:** TEST-ONLY  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

/// How far cross-copy homology extends OUTWARD from a gene boundary into the flank. Both
/// slices run in the SAME direction away from the gene: index 0 is the gene boundary, larger
/// indices go further into the flank (callers reverse the upstream flank so this holds).
/// Slides a `window`-bp window; returns the offset at which windowed identity first drops
/// below `min_identity` — the breakpoint distance from the gene boundary. Returns the full
/// compared length if homology never drops (the duplication extends past the fetched window).
pub fn flank_homology_extent(a: &[u8], b: &[u8], window: usize, min_identity: f64) -> usize {
    flank_homology_extent_banded(a, b, window, min_identity, 0)
}

/// Banded variant: tolerate flank INDELS up to `band` bp. The offset-locked scan above breaks
/// at the first indel (the offset drifts and identity collapses) even though the homology
/// continues — so a real segdup with a flank indel is truncated. Here, when the current-offset
/// window drops, we search forward (advancing past a window-straddling indel) and across offsets
/// in ±band for where homology RESUMES, absorbing the indel and continuing. A genuine breakpoint
/// (homology → random) finds no resuming window (random identity ≪ floor), so it is NOT
/// over-extended. `band == 0` reproduces the offset-locked scan exactly.
pub fn flank_homology_extent_banded(
    a: &[u8], b: &[u8], window: usize, min_identity: f64, band: usize,
) -> usize {
    let (na, nb) = (a.len(), b.len());
    if na == 0 || nb == 0 || window == 0 {
        return 0;
    }
    let w = window.min(na).min(nb);
    let thr = (min_identity * w as f64).ceil() as usize;
    // Match count of a[i..i+w] vs b[i+δ..i+δ+w], or None if either window is out of bounds.
    let count = |i: usize, delta: i64| -> Option<usize> {
        if i + w > na {
            return None;
        }
        let bj = i as i64 + delta;
        if bj < 0 || bj as usize + w > nb {
            return None;
        }
        let bj = bj as usize;
        Some((0..w).filter(|&k| a[i + k] == b[bj + k]).count())
    };
    let band = band as i64;
    let mut delta: i64 = 0;
    let mut i = 0usize;
    loop {
        match count(i, delta) {
            // Ran off a comparable window while homology was holding → extends to the end.
            None => return na.min(nb),
            Some(m) if m >= thr => i += 1,
            Some(_) => {
                // Drop. Search forward (s, to clear a window-straddling indel) × offsets
                // (δ', the indel shift) for where homology resumes; take the smallest advance,
                // best-aligning offset.
                let mut found: Option<(usize, i64)> = None;
                if band > 0 {
                    's_loop: for s in 0..=(w + band as usize) {
                        if i + s + w > na {
                            break;
                        }
                        let mut best = (0usize, delta);
                        for d in (delta - band)..=(delta + band) {
                            if let Some(m) = count(i + s, d) {
                                if m > best.0 {
                                    best = (m, d);
                                }
                            }
                        }
                        if best.0 >= thr {
                            found = Some((s, best.1));
                            break 's_loop;
                        }
                    }
                }
                match found {
                    Some((s, d)) => {
                        i += s; // the indel region counts toward the duplicated extent
                        delta = d;
                    }
                    None => return if i == 0 { 0 } else { i + w / 2 },
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct SegdupParams {
    pub window: usize,        // homology window (bp)
    pub min_identity: f64,    // windowed identity floor
    pub min_segdup_flank: u64, // total flank homology to call a segdup vs a bare paralog
    pub min_each_flank: u64,   // min homology on EACH side (a segdup has two breakpoints)
    pub band: usize,           // max flank indel (bp) the homology scan tolerates (0 = none)
}

impl Default for SegdupParams {
    fn default() -> Self {
        SegdupParams { window: 50, min_identity: 0.70, min_segdup_flank: 500, min_each_flank: 50, band: 30 }
    }
}

impl SegdupParams {
    pub fn from_env() -> Self {
        let mut p = SegdupParams::default();
        if let Some(v) = std::env::var("RUSTLE_VG_SEGDUP_WINDOW").ok().and_then(|s| s.parse().ok()) { p.window = v; }
        if let Some(v) = std::env::var("RUSTLE_VG_SEGDUP_MIN_ID").ok().and_then(|s| s.parse().ok()) { p.min_identity = v; }
        if let Some(v) = std::env::var("RUSTLE_VG_SEGDUP_MIN_FLANK").ok().and_then(|s| s.parse().ok()) { p.min_segdup_flank = v; }
        if let Some(v) = std::env::var("RUSTLE_VG_SEGDUP_MIN_EACH").ok().and_then(|s| s.parse().ok()) { p.min_each_flank = v; }
        if let Some(v) = std::env::var("RUSTLE_VG_SEGDUP_BAND").ok().and_then(|s| s.parse().ok()) { p.band = v; }
        p
    }
}

/// A called segmental-duplication extent for one copy pair.
#[derive(Debug, Clone, PartialEq)]
pub struct SegdupExtent {
    pub gene_span: u64,
    pub upstream_extent: u64,   // bp of homology upstream of the gene boundary
    pub downstream_extent: u64, // bp downstream
    pub total_extent: u64,      // upstream + gene + downstream (the duplicated segment span)
    pub is_segdup: bool,        // flank homology exceeds min_segdup_flank → true segdup
}

/// Call the duplicated-segment extent for one copy pair from the anchored flanks. `up_a`/`up_b`
/// are the upstream flanks REVERSED (index 0 = gene start, going outward); `down_a`/`down_b`
/// are the downstream flanks (index 0 = gene end, going outward).
pub fn call_segdup_extent(
    gene_span: u64,
    up_a: &[u8], up_b: &[u8],
    down_a: &[u8], down_b: &[u8],
    p: &SegdupParams,
) -> SegdupExtent {
    let up = flank_homology_extent_banded(up_a, up_b, p.window, p.min_identity, p.band) as u64;
    let down = flank_homology_extent_banded(down_a, down_b, p.window, p.min_identity, p.band) as u64;
    // A true segmental duplication has TWO breakpoints: require genuine homology on BOTH
    // sides (not one big one-sided match, which is a partial dup or a coincidence).
    let is_segdup = (up + down) >= p.min_segdup_flank && up.min(down) >= p.min_each_flank;
    SegdupExtent {
        gene_span,
        upstream_extent: up,
        downstream_extent: down,
        total_extent: up + gene_span + down,
        is_segdup,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Deterministic pseudo-random bytes (no Math.random in scope; vary by seed).
    fn seqr(n: usize, seed: u64) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut x = seed.wrapping_mul(2862933555777941757).wrapping_add(3037000493);
        (0..n).map(|_| { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); bases[((x >> 33) & 3) as usize] }).collect()
    }

    #[test]
    fn identical_flanks_extend_fully() {
        let s = seqr(1000, 1);
        assert_eq!(flank_homology_extent(&s, &s, 50, 0.70), 1000);
    }

    #[test]
    fn homologous_prefix_then_random_breaks_near_boundary() {
        // 600 bp shared homology, then divergent random tails.
        let shared = seqr(600, 7);
        let mut a = shared.clone(); a.extend(seqr(400, 11));
        let mut b = shared.clone(); b.extend(seqr(400, 99));
        let ext = flank_homology_extent(&a, &b, 50, 0.70);
        assert!((550..=650).contains(&ext), "extent {ext} not near the 600 bp boundary");
    }

    #[test]
    fn immediately_divergent_has_tiny_extent() {
        let a = seqr(500, 3);
        let b = seqr(500, 4);
        assert!(flank_homology_extent(&a, &b, 50, 0.70) < 50);
    }

    #[test]
    fn empty_is_zero() {
        assert_eq!(flank_homology_extent(&[], &[], 50, 0.70), 0);
        assert_eq!(flank_homology_extent(b"ACGT", &[], 50, 0.70), 0);
    }

    #[test]
    fn band_absorbs_a_flank_insertion() {
        // a = X(400) + INS(25) + Y(400) + tail;  b = X(400) + Y(400) + tail.
        // Offset-locked (band=0) breaks at the insertion (~400); banded extends through it.
        let x = seqr(400, 21);
        let y = seqr(400, 22);
        let ins = seqr(25, 23);
        let mut a = x.clone(); a.extend(&ins); a.extend(&y); a.extend(seqr(300, 24));
        let mut b = x.clone(); b.extend(&y); b.extend(seqr(300, 25));
        let locked = flank_homology_extent_banded(&a, &b, 50, 0.70, 0);
        let banded = flank_homology_extent_banded(&a, &b, 50, 0.70, 30);
        assert!(locked < 500, "offset-locked should break near the insertion, got {locked}");
        assert!(banded > 700, "banded should extend through the indel into Y, got {banded}");
    }

    #[test]
    fn band_absorbs_a_flank_deletion() {
        // a = X(400) + Y(400) + tail;  b = X(400) + EXTRA(25) + Y(400) + tail (deletion in a).
        let x = seqr(400, 31);
        let y = seqr(400, 32);
        let extra = seqr(25, 33);
        let mut a = x.clone(); a.extend(&y); a.extend(seqr(300, 34));
        let mut b = x.clone(); b.extend(&extra); b.extend(&y); b.extend(seqr(300, 35));
        let banded = flank_homology_extent_banded(&a, &b, 50, 0.70, 30);
        assert!(banded > 700, "banded should absorb the deletion and reach Y, got {banded}");
    }

    #[test]
    fn band_does_not_over_extend_a_genuine_breakpoint() {
        // Adversarial: homology then DIVERGENT random — banding must NOT find a spurious resume.
        let shared = seqr(500, 41);
        let mut a = shared.clone(); a.extend(seqr(500, 42));
        let mut b = shared.clone(); b.extend(seqr(500, 99));
        let banded = flank_homology_extent_banded(&a, &b, 50, 0.70, 30);
        assert!((450..=560).contains(&banded), "must stop at the real breakpoint, got {banded}");
    }

    #[test]
    fn band_zero_equals_offset_locked() {
        let s = seqr(800, 51);
        let mut a = s.clone(); a.extend(seqr(200, 52));
        let mut b = s.clone(); b.extend(seqr(200, 53));
        assert_eq!(
            flank_homology_extent(&a, &b, 50, 0.70),
            flank_homology_extent_banded(&a, &b, 50, 0.70, 0),
        );
    }

    #[test]
    fn segdup_vs_bare_paralog_classification() {
        let p = SegdupParams::default(); // min_segdup_flank = 500
        let shared = seqr(800, 5);
        // SEGDUP: 800 bp shared flank both sides.
        let mut up_a = shared.clone(); up_a.extend(seqr(200, 1));
        let mut up_b = shared.clone(); up_b.extend(seqr(200, 2));
        let seg = call_segdup_extent(10_000, &up_a, &up_b, &up_a, &up_b, &p);
        assert!(seg.is_segdup);
        assert!(seg.total_extent > 10_000 + 1000);
        // BARE PARALOG: flanks immediately divergent → not a segdup.
        let bare = call_segdup_extent(10_000, &seqr(1000, 8), &seqr(1000, 9), &seqr(1000, 10), &seqr(1000, 11), &p);
        assert!(!bare.is_segdup);
    }
}
