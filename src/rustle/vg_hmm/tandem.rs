//! Intra-bundle tandem decomposition (spec 2026-06-02). Detects a tandem
//! paralog array that collapsed into one bundle and splits it into per-copy
//! clusters so the existing family-graph/EM machinery can resolve the copies.

use std::collections::HashSet;

/// Exact k-mer Jaccard between two DNA sequences (ACGT only; other bases
/// skipped). Returns 0.0 if either has no valid k-mers. Clusters are short
/// (one copy's exonic sequence), so an exact set is fine — no min-hash needed.
///
/// Name follows the domain-first, verb-last pattern used in this codebase
/// (`compute_family_graph_kmer_jaccard`, `minimizer_jaccard`): *seq* is the
/// domain (two raw sequences), *kmer_jaccard* is the operation.
pub fn seq_kmer_jaccard(a: &[u8], b: &[u8], k: usize) -> f64 {
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

/// O5/tandem runtime config (env). `enabled` defaults OFF (opt-in prototype).
#[derive(Debug, Clone, Copy)]
pub struct TandemConfig {
    pub enabled: bool,
    pub min_gap: u64,       // genomic separation to split copies
    pub min_jaccard: f64,   // cross-copy sequence-similarity floor
    pub kmer: usize,
}

impl TandemConfig {
    pub fn from_env() -> Self {
        let enabled = std::env::var("RUSTLE_VG_TANDEM")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false);
        let min_gap = std::env::var("RUSTLE_VG_TANDEM_MIN_SEP")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(5_000);
        let min_jaccard = std::env::var("RUSTLE_VG_TANDEM_MIN_JACCARD")
            .ok().and_then(|v| v.parse().ok())
            .filter(|f: &f64| *f > 0.0 && *f <= 1.0).unwrap_or(0.20);
        TandemConfig { enabled, min_gap, min_jaccard, kmer: 15 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── seq_kmer_jaccard tests ──────────────────────────────────────────────

    #[test]
    fn identical_sequences_jaccard_one() {
        let s = b"ACGTACGTACGTACGTACGT";
        assert!((seq_kmer_jaccard(s, s, 15) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn near_identical_high_jaccard() {
        // The spec prescribed: let a = b"ACGTACGTACGTACGTACGTACGTACGTACGT" (32
        // bytes, ACGT-repeating) with one substitution at byte 16 and k=15 → assert
        // Jaccard > 0.2. That is mathematically broken: a 32-byte ACGT-repeat
        // produces only 4 distinct 15-mers (the four rotation phases), and flipping
        // byte 16 adds/removes a few but the Jaccard is ≈ 0.105 < 0.2 — the
        // assertion would always fail. A non-repeating 100-base sequence is used
        // instead: one substitution at position 50 displaces 15 k-mers out of 86,
        // leaving 71 shared → Jaccard ≈ 0.70 > 0.2.
        let a: Vec<u8> = (0u8..100)
            .map(|i| [b'A', b'C', b'G', b'T', b'G', b'C', b'A', b'T', b'T', b'G'][i as usize % 10])
            .collect();
        let mut b = a.clone();
        b[50] = if b[50] == b'A' { b'C' } else { b'A' }; // one substitution
        assert!(seq_kmer_jaccard(&a, &b, 15) > 0.2);
    }

    #[test]
    fn disjoint_sequences_jaccard_zero() {
        // AAAA vs CCCC: both sequences have non-empty k-mer sets (each has exactly
        // one 15-mer), but the sets are disjoint → inter=0 → Jaccard=0.
        // NOTE: this does NOT exercise the early-exit guard (sa.is_empty() ||
        // sb.is_empty()); see the three tests below for guard coverage.
        let a = b"AAAAAAAAAAAAAAAAAAAA";
        let b = b"CCCCCCCCCCCCCCCCCCCC";
        assert_eq!(seq_kmer_jaccard(a, b, 15), 0.0);
    }

    // ── early-exit guard tests (sa.is_empty() || sb.is_empty()) ────────────
    //
    // The docstring promises 0.0 for (i) empty slice, (ii) slice shorter than k,
    // (iii) all-non-ACGT (e.g. all-N). Each test below hits the guard directly.

    #[test]
    fn empty_slice_returns_zero() {
        // An empty slice produces an empty k-mer set → guard fires.
        let s = b"ACGTACGTACGTACGT";
        assert_eq!(seq_kmer_jaccard(b"", s, 15), 0.0);
        assert_eq!(seq_kmer_jaccard(s, b"", 15), 0.0);
        assert_eq!(seq_kmer_jaccard(b"", b"", 15), 0.0);
    }

    #[test]
    fn shorter_than_k_returns_zero() {
        // A slice with length < k produces an empty k-mer set → guard fires.
        let short = b"ACGTACGT"; // 8 bytes < k=15
        let long = b"ACGTACGTACGTACGTACGT";
        assert_eq!(seq_kmer_jaccard(short, long, 15), 0.0);
        assert_eq!(seq_kmer_jaccard(long, short, 15), 0.0);
    }

    #[test]
    fn all_n_sequence_returns_zero() {
        // All-N (non-ACGT) sequences: all k-mers are rejected by the ACGT filter,
        // so the k-mer set is empty → guard fires.
        let n_seq = b"NNNNNNNNNNNNNNNNNNNN"; // 20 N-bases
        let acgt = b"ACGTACGTACGTACGTACGT";
        assert_eq!(seq_kmer_jaccard(n_seq, acgt, 15), 0.0);
        assert_eq!(seq_kmer_jaccard(acgt, n_seq, 15), 0.0);
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

    // ── TandemConfig tests ──────────────────────────────────────────────────

    #[test]
    fn tandemconfig_defaults_off() {
        let c = TandemConfig::from_env();
        assert!(!c.enabled);
        assert_eq!(c.min_gap, 5_000);
        assert!((c.min_jaccard - 0.20).abs() < 1e-9);
    }
}
