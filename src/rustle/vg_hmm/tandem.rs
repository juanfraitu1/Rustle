//! Intra-bundle tandem decomposition (spec 2026-06-02). Detects a tandem
//! paralog array that collapsed into one bundle and splits it into per-copy
//! clusters so the existing family-graph/EM machinery can resolve the copies.

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
