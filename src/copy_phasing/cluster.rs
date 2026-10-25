//! Clustering (o3port.py `label_propagation`, `correlation_clustering`, `correct_patterns`, `cluster_pile`, `haplotype`,
//! `conflicting_set`, `absorb`).

use super::stats::{compatible, counter, informative, most_common, patterns};
use super::{Clusterer, Cols, Hap, Obs, Pattern, ReadId, MAX_ROUNDS, MIN_Q, MIN_READS};
use std::collections::{BTreeMap, BTreeSet};

/// A cluster: a set of reads; lists of clusters are sorted by (-size, sorted reads) as the Python sorts them.
pub type Cluster = BTreeSet<ReadId>;

/// Distinct non-empty patterns, `sorted(groups, key=lambda k: (-len(groups[k]), k))`, with their reads.
fn pattern_groups(pat: &BTreeMap<ReadId, Pattern>) -> Vec<(Vec<(i64, u8)>, Vec<ReadId>)> {
    let mut groups: BTreeMap<Vec<(i64, u8)>, Vec<ReadId>> = BTreeMap::new();
    for (&r, p) in pat {
        if !p.is_empty() {
            groups.entry(p.iter().map(|(&c, &b)| (c, b)).collect()).or_default().push(r);
        }
    }
    let mut keys: Vec<(Vec<(i64, u8)>, Vec<ReadId>)> = groups.into_iter().collect();
    keys.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(&b.0)));
    keys
}

fn sort_clusters(mut cs: Vec<Cluster>) -> Vec<Cluster> {
    cs.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a.iter().cmp(b.iter())));
    cs
}

fn by_label(keys: &[(Vec<(i64, u8)>, Vec<ReadId>)], label: &[usize]) -> Vec<Cluster> {
    let mut clusters: BTreeMap<usize, Cluster> = BTreeMap::new();
    for (i, (_, rs)) in keys.iter().enumerate() {
        clusters.entry(label[i]).or_default().extend(rs.iter().copied());
    }
    sort_clusters(clusters.into_values().collect())
}

pub fn label_propagation(pat: &BTreeMap<ReadId, Pattern>) -> Vec<Cluster> {
    todo!()
}

pub fn correlation_clustering(pat: &BTreeMap<ReadId, Pattern>) -> Vec<Cluster> {
    todo!()
}

pub fn correct_patterns(pat: &BTreeMap<ReadId, Pattern>) -> BTreeMap<ReadId, Pattern> {
    todo!()
}

fn clusterer(method: Clusterer, pat: &BTreeMap<ReadId, Pattern>) -> Vec<Cluster> {
    match method {
        Clusterer::Cc => correlation_clustering(pat),
        Clusterer::Lp => label_propagation(pat),
    }
}

/// `cluster_pile`: (clusters of >= MIN_READS reads, sorted; the pile-level informative columns). The Python's middle product (every
/// column used) is never read and is not returned.
pub fn cluster_pile(reads: &BTreeSet<ReadId>, obs: &Obs, method: Clusterer, correct: bool) -> (Vec<Cluster>, Cols) {
    todo!()
}

pub fn haplotype(cluster: &Cluster, obs: &Obs, cols: &Cols) -> Hap {
    todo!()
}

fn differ(a: &Hap, b: &Hap) -> (bool, usize) {
    let mut shared = false;
    let mut n = 0;
    for (c, x) in a {
        if let Some(y) = b.get(c) {
            shared = true;
            n += (x != y) as usize;
        }
    }
    (shared, n)
}

pub fn conflicting_set(haps: &[Hap]) -> Vec<usize> {
    todo!()
}

pub fn absorb(hs: &[Hap], keep: &[usize]) -> BTreeMap<usize, usize> {
    todo!()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::copy_phasing::ObsX;

    /// The Python tests' `obs_from`, with read ids assigned in sorted-name order; returns (obs, reads, id -> name).
    pub(crate) fn obs_from(haps: &[(&str, &[(i64, u8)])], counts: &[usize]) -> (Obs, BTreeSet<ReadId>, Vec<String>) {
        let mut names: Vec<String> = Vec::new();
        for (k, (name, _)) in haps.iter().enumerate() {
            names.extend((0..counts[k]).map(|i| format!("{name}{i}")));
        }
        let mut sorted = names.clone();
        sorted.sort();
        let id = |n: &str| sorted.binary_search_by(|x| x.as_str().cmp(n)).unwrap() as ReadId;
        let mut obs = Obs::new();
        for (k, (name, h)) in haps.iter().enumerate() {
            for i in 0..counts[k] {
                let r = id(&format!("{name}{i}"));
                for &(c, b) in h.iter() {
                    obs.entry(c).or_default().push(ObsX { read: r, base: b, qual: 30, eps: None });
                }
            }
        }
        let reads = (0..sorted.len() as ReadId).collect();
        (obs, reads, sorted)
    }

    fn hap(xs: &[(i64, u8)]) -> Hap {
        xs.iter().copied().collect()
    }

    fn stems(c: &Cluster, names: &[String]) -> Vec<String> {
        c.iter().map(|&r| names[r as usize].trim_end_matches(|ch: char| ch.is_ascii_digit()).to_string()).collect::<BTreeSet<_>>().into_iter().collect()
    }

    #[test]
    fn three_haplotypes_on_one_pile_including_a_minor_one() {
        let (a, b, c): (&[(i64, u8)], &[(i64, u8)], &[(i64, u8)]) = (
            &[(1, b'A'), (2, b'C'), (3, b'G'), (4, b'A')],
            &[(1, b'T'), (2, b'T'), (3, b'G'), (4, b'A')],
            &[(1, b'A'), (2, b'T'), (3, b'C'), (4, b'G')],
        );
        let (obs, reads, names) = obs_from(&[("A", a), ("B", b), ("C", c)], &[45, 45, 9]);
        let (clusters, top) = cluster_pile(&reads, &obs, Clusterer::Cc, false);
        let mut got: Vec<Vec<String>> = clusters.iter().map(|c| stems(c, &names)).collect();
        got.sort();
        assert_eq!(got, vec![vec!["A"], vec!["B"], vec!["C"]]);
        let hs: Vec<Hap> = clusters.iter().map(|c| haplotype(c, &obs, &top)).collect();
        assert_eq!(conflicting_set(&hs).len(), 3);
    }

    #[test]
    fn bridge_reads_do_not_merge_haplotypes() {
        let (obs0, _, names0) = obs_from(&[("A", &[(1, b'A'), (2, b'C')]), ("C", &[(1, b'A'), (2, b'T')])], &[20, 20]);
        // the six bridge reads sort before "C*": rebuild with them in the table
        let mut names = names0.clone();
        names.extend((0..6).map(|i| format!("bridge{i}")));
        names.sort();
        let id = |n: &str| names.binary_search_by(|x| x.as_str().cmp(n)).unwrap() as ReadId;
        let mut obs = Obs::new();
        for (c, xs) in &obs0 {
            obs.insert(*c, xs.iter().map(|x| ObsX { read: id(&names0[x.read as usize]), ..*x }).collect());
        }
        for i in 0..6 {
            obs[&1].push(ObsX { read: id(&format!("bridge{i}")), base: b'A', qual: 30, eps: None });
        }
        let reads: BTreeSet<ReadId> = (0..names.len() as ReadId).collect();
        let (clusters, _) = cluster_pile(&reads, &obs, Clusterer::Lp, false);
        let big: Vec<&Cluster> = clusters.iter().filter(|c| c.len() >= 10).collect();
        assert_eq!(big.len(), 2);
        for c in big {
            let s: BTreeSet<String> = c.iter().map(|&r| names[r as usize].clone()).filter(|n| !n.starts_with("bridge")).map(|n| n[..1].to_string()).collect();
            assert_eq!(s.len(), 1);
        }
    }

    #[test]
    fn fragments_of_one_haplotype_do_not_conflict() {
        assert_eq!(conflicting_set(&[hap(&[(1, b'A'), (2, b'C')]), hap(&[(7, b'G'), (8, b'T')])]).len(), 1);
    }

    #[test]
    fn a_single_differing_column_is_not_a_conflict() {
        assert_eq!(conflicting_set(&[hap(&[(1, b'A'), (2, b'C'), (3, b'G')]), hap(&[(1, b'A'), (2, b'C'), (3, b'T')])]).len(), 1);
        assert_eq!(conflicting_set(&[hap(&[(1, b'A'), (2, b'C'), (3, b'G')]), hap(&[(1, b'A'), (2, b'T'), (3, b'T')])]).len(), 2);
    }

    #[test]
    fn absorb_joins_a_fragment_to_its_one_compatible_haplotype() {
        let hs = [
            hap(&[(1, b'A'), (2, b'C'), (3, b'G')]),
            hap(&[(1, b'T'), (2, b'T'), (3, b'C')]),
            hap(&[(2, b'C'), (3, b'G')]),
            hap(&[(1, b'A'), (2, b'T')]),
            hap(&[(9, b'G')]),
        ];
        assert_eq!(absorb(&hs, &[0, 1]), BTreeMap::from([(2, 0)]));
    }

    fn pat_of(xs: &[(i64, u8)]) -> Pattern {
        xs.iter().copied().collect()
    }

    #[test]
    fn correlation_clustering_three_haplotypes_and_errors() {
        // names: A0..A29, Aerr, B0..B29, C0..C7, bridge (sorted ids)
        let mut names: Vec<String> = (0..30).map(|i| format!("A{i}")).chain((0..30).map(|i| format!("B{i}"))).chain((0..8).map(|i| format!("C{i}"))).collect();
        names.push("Aerr".into());
        names.push("bridge".into());
        names.sort();
        let id = |n: &str| names.binary_search_by(|x| x.as_str().cmp(n)).unwrap() as ReadId;
        let mut pat = BTreeMap::new();
        for i in 0..30 {
            pat.insert(id(&format!("A{i}")), pat_of(&[(1, b'A'), (2, b'C'), (3, b'G'), (4, b'A')]));
            pat.insert(id(&format!("B{i}")), pat_of(&[(1, b'T'), (2, b'T'), (3, b'G'), (4, b'C')]));
        }
        for i in 0..8 {
            pat.insert(id(&format!("C{i}")), pat_of(&[(1, b'A'), (2, b'T'), (3, b'C'), (4, b'G')]));
        }
        pat.insert(id("Aerr"), pat_of(&[(1, b'A'), (2, b'C'), (3, b'G'), (4, b'G')]));
        pat.insert(id("bridge"), pat_of(&[(1, b'A')]));
        let got = correlation_clustering(&pat);
        let st: Vec<Vec<String>> = got.iter().map(|c| stems(c, &names)).collect();
        assert!(st[0] == vec!["A", "Aerr"] || st[0] == vec!["A", "Aerr", "bridge"]);
        assert!(st.contains(&vec!["B".to_string()]));
        assert!(st.iter().any(|n| n == &vec!["C"] || n == &vec!["C", "bridge"]));
    }

    #[test]
    fn correlation_clustering_keeps_a_one_column_difference_together() {
        let pat: BTreeMap<ReadId, Pattern> =
            (0..5).map(|i| (i, pat_of(&[(1, b'A'), (2, b'C')]))).chain((5..10).map(|i| (i, pat_of(&[(1, b'A'), (2, b'T')])))).collect();
        assert_eq!(correlation_clustering(&pat).len(), 1);
    }

    #[test]
    fn correct_patterns_replaces_a_lone_error() {
        // 4 reads AC at (1, 2), one read AT: its neighbours (agree at column 1) carry C at column 2, none T -> corrected to C
        let mut pat: BTreeMap<ReadId, Pattern> = (0..4).map(|i| (i, pat_of(&[(1, b'A'), (2, b'C')]))).collect();
        pat.insert(4, pat_of(&[(1, b'A'), (2, b'T')]));
        let f = correct_patterns(&pat);
        assert_eq!(f[&4], pat_of(&[(1, b'A'), (2, b'C')]));
        assert_eq!(f[&0], pat_of(&[(1, b'A'), (2, b'C')]));
    }
}
