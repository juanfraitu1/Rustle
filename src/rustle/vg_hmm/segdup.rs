//! Segmental-duplication extent / breakpoint calling from GENOME homology.
//!
//! A bare gene paralog shares homology only across the gene (exons/introns). A true
//! SEGMENTAL DUPLICATION copied the gene PLUS flanking sequence, so cross-copy homology
//! extends into the intergenic flanks; the duplication BREAKPOINTS are where that homology
//! ends. IsoSeq reads are spliced mRNA and do not cover intergenic flanks, so this signal
//! comes from the GENOME sequence (which the VG already loads for family discovery), anchored
//! at the gene the family discovered. The core is a pure anchored sliding-window homology scan.

/// How far cross-copy homology extends OUTWARD from a gene boundary into the flank. Both
/// slices run in the SAME direction away from the gene: index 0 is the gene boundary, larger
/// indices go further into the flank (callers reverse the upstream flank so this holds).
/// Slides a `window`-bp window; returns the offset at which windowed identity first drops
/// below `min_identity` — the breakpoint distance from the gene boundary. Returns the full
/// compared length if homology never drops (the duplication extends past the fetched window).
pub fn flank_homology_extent(a: &[u8], b: &[u8], window: usize, min_identity: f64) -> usize {
    let n = a.len().min(b.len());
    if n == 0 || window == 0 {
        return 0;
    }
    let w = window.min(n);
    let mut matches = (0..w).filter(|&k| a[k] == b[k]).count();
    let mut i = 0usize;
    loop {
        if (matches as f64) / (w as f64) < min_identity {
            // Window [i, i+w) fell below identity → breakpoint is within it. Report its center
            // for a less window-biased estimate — but if the ANCHOR window (i==0) already fails,
            // there's no flank homology at all; report 0, not a phantom w/2 (bare paralogs would
            // otherwise show a spurious ~w/2 flank extent).
            return if i == 0 { 0 } else { i + w / 2 };
        }
        if i + w >= n {
            return n; // homology held to the end of the fetched flank
        }
        if a[i] == b[i] {
            matches -= 1;
        }
        if a[i + w] == b[i + w] {
            matches += 1;
        }
        i += 1;
    }
}

#[derive(Debug, Clone)]
pub struct SegdupParams {
    pub window: usize,        // homology window (bp)
    pub min_identity: f64,    // windowed identity floor
    pub min_segdup_flank: u64, // total flank homology to call a segdup vs a bare paralog
    pub min_each_flank: u64,   // min homology on EACH side (a segdup has two breakpoints)
}

impl Default for SegdupParams {
    fn default() -> Self {
        SegdupParams { window: 50, min_identity: 0.70, min_segdup_flank: 500, min_each_flank: 50 }
    }
}

impl SegdupParams {
    pub fn from_env() -> Self {
        let mut p = SegdupParams::default();
        if let Some(v) = std::env::var("RUSTLE_VG_SEGDUP_WINDOW").ok().and_then(|s| s.parse().ok()) { p.window = v; }
        if let Some(v) = std::env::var("RUSTLE_VG_SEGDUP_MIN_ID").ok().and_then(|s| s.parse().ok()) { p.min_identity = v; }
        if let Some(v) = std::env::var("RUSTLE_VG_SEGDUP_MIN_FLANK").ok().and_then(|s| s.parse().ok()) { p.min_segdup_flank = v; }
        if let Some(v) = std::env::var("RUSTLE_VG_SEGDUP_MIN_EACH").ok().and_then(|s| s.parse().ok()) { p.min_each_flank = v; }
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
    let up = flank_homology_extent(up_a, up_b, p.window, p.min_identity) as u64;
    let down = flank_homology_extent(down_a, down_b, p.window, p.min_identity) as u64;
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
