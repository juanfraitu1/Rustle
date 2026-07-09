//! VG re-align supplement -- re-align poor-fit/unmapped reads to O1's family copy-paths,
//! significance-gated (correct + discover). Task 1: candidate selection. Task 3: re-align a
//! candidate read to the family's copy-paths (identity-based, DRY on `bridge_detector::aln_id`)
//! and route unmapped reads to candidate families by shared minimizers.

use std::collections::HashSet;

use crate::vg_family::bridge_detector::aln_id;
use crate::vg_family::minimizers::minimizers;

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

/// Minimal identity floor for `realign_to_paths`: a read scoring below this against every
/// copy-path in the family fits none of them well enough to be a re-align candidate at all.
/// (Genuine novel copies pooled from low-fit reads are handled separately, by Task 4 --
/// this floor only gates obvious non-members out of the per-family re-align step.)
pub const MIN_ALN_ID: f64 = 0.5;

/// A candidate read's best fit among a family's copy-paths (Task 3), used by Task 4 to decide
/// whether re-aligning to the family beats the read's existing linear-locus alignment.
pub struct RealignHit {
    /// Index into `copy_seqs` of the best-fitting copy-path (earliest index on ties).
    pub best_copy: usize,
    /// `aln_id(read_seq, copy_seqs[best_copy])` -- the best copy-path identity.
    pub id_best: f64,
    /// `aln_id(read_seq, copy_seqs[linear_copy])` when `linear_copy` is `Some` and in range,
    /// else `0.0` -- the read's fit to the copy it would otherwise be attributed to linearly.
    pub id_linear: f64,
}

/// Re-align `read_seq` to every copy-path in `copy_seqs`, reusing `bridge_detector::aln_id` as
/// the fit score (best infix identity, tries both orientations). Returns the best-fitting copy
/// plus (for Task 4's accept comparison) the fit to `linear_copy`'s sequence, if given.
///
/// Returns `None` when `copy_seqs` is empty, or when the best fit is below `MIN_ALN_ID`: a read
/// that fits no copy-path at all is not a re-align candidate for this family.
pub fn realign_to_paths(
    read_seq: &[u8],
    copy_seqs: &[Vec<u8>],
    linear_copy: Option<usize>,
) -> Option<RealignHit> {
    if copy_seqs.is_empty() {
        return None;
    }

    let mut best_copy = 0usize;
    let mut id_best = f64::NEG_INFINITY;
    for (k, seq) in copy_seqs.iter().enumerate() {
        let id = aln_id(read_seq, seq);
        // Strict `>` keeps the earliest index on ties.
        if id > id_best {
            id_best = id;
            best_copy = k;
        }
    }

    if id_best < MIN_ALN_ID {
        return None;
    }

    let id_linear = match linear_copy {
        Some(lc) if lc < copy_seqs.len() => aln_id(read_seq, &copy_seqs[lc]),
        _ => 0.0,
    };

    Some(RealignHit { best_copy, id_best, id_linear })
}

/// Route an unmapped read to candidate families by shared canonical minimizers: compute the
/// read's minimizer value-set (`minimizers(seq, MINIMIZER_K, MINIMIZER_W)`, values only -- the
/// masked flag doesn't matter for routing) and, for each `(family_id, consensus)` pair, count
/// how many minimizer values the read shares with that family's consensus. Returns the
/// `family_id`s meeting `count >= min_shared`, sorted ascending.
pub fn route_unmapped(
    seq: &[u8],
    family_consensuses: &[(usize, Vec<u8>)],
    min_shared: usize,
) -> Vec<usize> {
    use crate::vg_family::minimizers::{MINIMIZER_K, MINIMIZER_W};

    let read_mins: HashSet<u64> =
        minimizers(seq, MINIMIZER_K, MINIMIZER_W).into_iter().map(|(v, _)| v).collect();

    let mut hits: Vec<usize> = family_consensuses
        .iter()
        .filter_map(|(family_id, consensus)| {
            let cons_mins: HashSet<u64> =
                minimizers(consensus, MINIMIZER_K, MINIMIZER_W).into_iter().map(|(v, _)| v).collect();
            let shared = read_mins.intersection(&cons_mins).count();
            (shared >= min_shared).then_some(*family_id)
        })
        .collect();

    hits.sort_unstable();
    hits
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

    /// Deterministic pseudo-random ACGT sequence generator (xorshift64), so test sequences are
    /// reproducible without hand-typing long strings or depending on real fixture data.
    fn pseudo_seq(seed: u64, len: usize) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        // splitmix64-style mix so nearby/small seeds (1, 2, 3, ...) don't collapse to the same
        // or correlated xorshift states (plain `seed | 1` made seeds 2 and 3 identical).
        let mut state =
            seed.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(0x2545F4914F6CDD1D) | 1;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                bases[(state % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn realign_picks_best_copy_path() {
        let copy0 = pseudo_seq(1, 60);
        let copy1 = pseudo_seq(2, 60);
        let copy2 = pseudo_seq(3, 60);
        let copy_seqs = vec![copy0.clone(), copy1.clone(), copy2.clone()];

        // Read == copy 1 exactly -> best_copy == 1, id_best ~1.0.
        let hit = realign_to_paths(&copy1, &copy_seqs, Some(0))
            .expect("exact copy-path match must be a candidate");
        assert_eq!(hit.best_copy, 1);
        assert!(hit.id_best > 0.99, "id_best = {}", hit.id_best);
        // id_linear (fit to copy 0, the "linear locus" copy) must be strictly lower than the
        // true best-copy fit -- the read really belongs to copy 1, not copy 0.
        assert!(
            hit.id_linear < hit.id_best,
            "id_linear = {} should be < id_best = {}",
            hit.id_linear,
            hit.id_best
        );

        // Read == copy 2 with a couple of substitutions -> still best_copy == 2, id_best < 1.0.
        let mut mutated = copy2.clone();
        // Flip two bases to a base guaranteed different from the original.
        for &pos in &[10usize, 40usize] {
            let orig = mutated[pos];
            mutated[pos] = [b'A', b'C', b'G', b'T'].into_iter().find(|&b| b != orig).unwrap();
        }
        let hit2 = realign_to_paths(&mutated, &copy_seqs, None)
            .expect("near-exact copy-path match must be a candidate");
        assert_eq!(hit2.best_copy, 2);
        assert!(hit2.id_best < 1.0, "id_best = {}", hit2.id_best);
        assert!(hit2.id_best > 0.9, "id_best = {} should still be a strong fit", hit2.id_best);
        // No linear_copy given -> id_linear is the documented 0.0 sentinel.
        assert_eq!(hit2.id_linear, 0.0);
    }

    #[test]
    fn realign_returns_none_for_nonmember() {
        let copy_seqs = vec![pseudo_seq(1, 60), pseudo_seq(2, 60), pseudo_seq(3, 60)];
        // A read unrelated to any copy-path (seed chosen so its best infix identity to all
        // three copies lands below MIN_ALN_ID; random same-length sequences under free-end-gap
        // edit distance land around ~0.5 identity by chance, so this isn't a free lunch).
        let read = pseudo_seq(3130, 60);
        let hit = realign_to_paths(&read, &copy_seqs, None);
        match hit {
            None => {}
            Some(h) => panic!("expected None for a non-member read, got id_best = {}", h.id_best),
        }

        // Empty copy_seqs -> always None regardless of the read.
        assert!(realign_to_paths(&read, &[], None).is_none());
    }

    #[test]
    fn route_unmapped_matches_by_minimizers() {
        // Two family consensuses, long enough (>= MINIMIZER_K + several windows) to produce
        // multiple minimizers.
        let consensus_a = pseudo_seq(10, 150);
        let consensus_b = pseudo_seq(20, 150);

        // A read that's an exact interior substring of family A's consensus -- interior so its
        // minimizer windows are fully contained in windows also present in the full consensus,
        // guaranteeing shared minimizer VALUES (not just positions).
        let read: Vec<u8> = consensus_a[40..100].to_vec();

        let families = vec![(0usize, consensus_a.clone()), (1usize, consensus_b.clone())];

        // Sanity: the read and consensus B should share few/no minimizers (unrelated random
        // sequences), while read and consensus A share several -- confirm before picking
        // min_shared so the test isn't tuned to a lucky threshold.
        use crate::vg_family::minimizers::{MINIMIZER_K, MINIMIZER_W};
        let mins_of = |s: &[u8]| -> HashSet<u64> {
            minimizers(s, MINIMIZER_K, MINIMIZER_W).into_iter().map(|(v, _)| v).collect()
        };
        let read_mins = mins_of(&read);
        let a_mins = mins_of(&consensus_a);
        let b_mins = mins_of(&consensus_b);
        let shared_a = read_mins.intersection(&a_mins).count();
        let shared_b = read_mins.intersection(&b_mins).count();
        assert!(shared_a >= 2, "expected several shared minimizers with A, got {shared_a}");
        assert!(shared_b < shared_a, "B should share fewer minimizers than A ({shared_b} vs {shared_a})");

        let min_shared = shared_b + 1; // strictly above B's count, at or below A's count
        assert!(min_shared <= shared_a, "min_shared {min_shared} must still be reachable by A");

        let routed = route_unmapped(&read, &families, min_shared);
        assert_eq!(routed, vec![0], "expected only family A routed, got {routed:?}");
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
