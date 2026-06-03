//! Intra-bundle tandem decomposition (spec 2026-06-02). Detects a tandem
//! paralog array that collapsed into one bundle and splits it into per-copy
//! clusters so the existing family-graph/EM machinery can resolve the copies.

use std::collections::HashSet;

/// Exact k-mer Jaccard between two DNA sequences (ACGT only; other bases
/// skipped). Returns 0.0 if either has no valid k-mers. Clusters are short
/// (one copy's exonic sequence), so an exact set is fine — no min-hash needed.
pub fn cluster_kmer_jaccard(a: &[u8], b: &[u8], k: usize) -> f64 {
    fn kset(s: &[u8], k: usize) -> HashSet<&[u8]> {
        let mut set = HashSet::new();
        if s.len() >= k {
            for w in s.windows(k) {
                if w.iter().all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
                    set.insert(w);
                }
            }
        }
        set
    }
    let sa = kset(a, k);
    let sb = kset(b, k);
    if sa.is_empty() || sb.is_empty() {
        return 0.0;
    }
    let inter = sa.intersection(&sb).count() as f64;
    let union = sa.union(&sb).count() as f64;
    inter / union
}

/// Cluster `[ (ref_start, ref_end) ]` read footprints into genomically-ordered
/// groups separated by gaps `> min_gap`. Returns, per cluster, the index list
/// into `intervals` plus the cluster's (start,end) span. Genomic-ascending.
pub fn cluster_reads_by_position(
    intervals: &[(u64, u64)],
    min_gap: u64,
) -> Vec<(Vec<usize>, (u64, u64))> {
    if intervals.is_empty() {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..intervals.len()).collect();
    order.sort_by_key(|&i| intervals[i].0);

    let mut clusters: Vec<(Vec<usize>, (u64, u64))> = Vec::new();
    for &i in &order {
        let (s, e) = intervals[i];
        match clusters.last_mut() {
            Some((idxs, span)) if s <= span.1.saturating_add(min_gap) => {
                idxs.push(i);
                span.1 = span.1.max(e);
                span.0 = span.0.min(s);
            }
            _ => clusters.push((vec![i], (s, e))),
        }
    }
    clusters
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── cluster_kmer_jaccard tests ──────────────────────────────────────────

    #[test]
    fn identical_sequences_jaccard_one() {
        let s = b"ACGTACGTACGTACGTACGT";
        assert!((cluster_kmer_jaccard(s, s, 15) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn near_identical_high_jaccard() {
        // 100-base non-repeating sequence: one substitution at position 50 affects
        // only 15 k-mers out of 86 total, leaving 71 shared → Jaccard ≈ 0.70 > 0.2.
        let a: Vec<u8> = (0u8..100)
            .map(|i| [b'A', b'C', b'G', b'T', b'G', b'C', b'A', b'T', b'T', b'G'][i as usize % 10])
            .collect();
        let mut b = a.clone();
        b[50] = if b[50] == b'A' { b'C' } else { b'A' }; // one substitution
        assert!(cluster_kmer_jaccard(&a, &b, 15) > 0.2);
    }

    #[test]
    fn disjoint_sequences_jaccard_zero() {
        let a = b"AAAAAAAAAAAAAAAAAAAA";
        let b = b"CCCCCCCCCCCCCCCCCCCC";
        assert_eq!(cluster_kmer_jaccard(a, b, 15), 0.0);
    }

    // ── cluster_reads_by_position tests ────────────────────────────────────

    #[test]
    fn splits_two_separated_copies() {
        // two copies ~23kb apart, gap threshold 5kb -> 2 clusters.
        let ivs = [(1000, 14000), (1500, 14500), (37000, 50000), (37200, 50200)];
        let c = cluster_reads_by_position(&ivs, 5000);
        assert_eq!(c.len(), 2);
        assert_eq!(c[0].0.len(), 2);
        assert_eq!(c[1].0.len(), 2);
        assert_eq!(c[0].1, (1000, 14500));
        assert_eq!(c[1].1, (37000, 50200));
        // Primary contract: returned indices must map back to the original input
        // slots so callers can recover per-read metadata. A sort permutation that
        // re-ordered indices incorrectly would break all callers silently.
        assert_eq!(c[0].0, vec![0, 1], "first cluster must contain original indices 0 and 1");
        assert_eq!(c[1].0, vec![2, 3], "second cluster must contain original indices 2 and 3");
    }

    #[test]
    fn single_cluster_when_within_gap() {
        let ivs = [(1000, 14000), (16000, 30000)]; // 2kb gap < 5kb
        assert_eq!(cluster_reads_by_position(&ivs, 5000).len(), 1);
    }

    #[test]
    fn empty_input() {
        assert!(cluster_reads_by_position(&[], 5000).is_empty());
    }
}
