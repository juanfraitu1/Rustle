//! VG re-align supplement -- re-align poor-fit/unmapped reads to O1's family copy-paths,
//! significance-gated (correct + discover). Task 1: candidate selection.

/// Thresholds for flagging a read as a re-align CANDIDATE (poor-fit or unmapped).
///
/// A read is a candidate if it is low-MAPQ (ambiguous/multi-mapped), heavily soft/hard-clipped
/// (partial alignment, suggesting the reference copy it landed on isn't its true source), or
/// highly divergent (many mismatches relative to the copy it aligned to). Divergence and clip
/// fraction are computed by the CALLER from the CIGAR/NM at wiring time (later task); this
/// struct only carries the thresholds `is_candidate` applies.
pub struct RealignParams {
    pub max_mapq: u8,
    pub min_clip_frac: f64,
    pub min_div: f64,
    pub min_reads: usize,
}

impl Default for RealignParams {
    fn default() -> Self {
        RealignParams { max_mapq: 20, min_clip_frac: 0.20, min_div: 0.05, min_reads: 3 }
    }
}

/// True iff `(mapq, div, clip_frac)` indicates a poor primary-alignment fit under `p`'s
/// thresholds: low MAPQ (`<= max_mapq`), OR heavy clipping (`>= min_clip_frac`), OR high
/// divergence (`>= min_div`). Pure function of the three scalars + params -- `div` and
/// `clip_frac` are computed by the caller from the CIGAR/NM, not derived here.
pub fn is_candidate(mapq: u8, div: f64, clip_frac: f64, p: &RealignParams) -> bool {
    mapq <= p.max_mapq || clip_frac >= p.min_clip_frac || div >= p.min_div
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_candidate_flags_poor_fit() {
        let p = RealignParams::default();

        // low MAPQ alone -> true
        assert!(is_candidate(5, 0.0, 0.0, &p));
        // high clip alone -> true
        assert!(is_candidate(60, 0.0, 0.30, &p));
        // high divergence alone -> true
        assert!(is_candidate(60, 0.08, 0.0, &p));
        // clean read: high MAPQ, no clip, low div -> false
        assert!(!is_candidate(60, 0.0, 0.0, &p));

        // boundary: mapq == max_mapq is still <= -> true
        assert!(is_candidate(20, 0.0, 0.0, &p));
        // boundary: just under both clip and div thresholds, and mapq above max -> false
        assert!(!is_candidate(21, 0.049, 0.19, &p));
    }

    #[test]
    fn realign_params_defaults() {
        let p = RealignParams::default();
        assert_eq!(p.max_mapq, 20);
        assert_eq!(p.min_clip_frac, 0.20);
        assert_eq!(p.min_div, 0.05);
        assert_eq!(p.min_reads, 3);
    }
}
