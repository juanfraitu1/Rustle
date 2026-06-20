//! Read-conflict (mutual-mappability) family criterion — the OPERATIONAL definition of a multi-copy family at
//! the RNA level (`bench/family_def_readconflict.md`).
//!
//! The de-novo pipeline currently defines a family by SEQUENCE similarity (`family_detect`: contiguous-core
//! `>= T_CORE`). That conflates true paralogs with DOMAIN-SHARERS (genes sharing one exon) and rests on an
//! arbitrary threshold. The operational definition instead asks the question the copy-assignment problem
//! actually cares about: **do reads cross-map between these loci?** Two loci are linked iff some read has a
//! placement in BOTH with TIED alignment scores (a genuine alternative placement — the multimapping conflict,
//! Canzar's conflict graph). A family is a connected component of that graph.
//!
//! Why this is the right unit: (1) no tuned similarity threshold — the boundary is the alignment-score tie, a
//! property of the data; (2) reads never cross-map outside their component, so the assignment problem
//! decomposes EXACTLY across families with no information lost; (3) it never groups domain-sharers (a read
//! over a shared exon maps to one locus — no alternative placement), and it picks out exactly the families
//! where assignment is needed (validated on Compara labels: 0 conflict on 7/7 domain-sharers, fires on
//! RABL2/APOBEC3, silent on the resolvable RFPL).
//!
//! This is the portable KERNEL: `conflict_edges` (read placements → weighted edges) + `conflict_families`
//! (edges → connected-component families). The remaining integration is plumbing per-locus secondary
//! placements (`secondary_index` / `tied_secondary_reads_in_region`) into the detection stage.

/// One read's placement on a candidate locus: the locus index plus the signals the conflict criterion uses —
/// `de` (gap-compressed divergence, the tie discriminant), `mapq` (both-0 = genuine-multimapper corroboration),
/// and `as_score` (kept only to log the `de-tie ⊆ AS-tie` invariant).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Placement {
    pub locus: usize,
    pub de: f32,
    pub mapq: u8,
    pub as_score: i32,
}

/// One read's placements over the family's candidate loci. Built from the BAM by the adapter.
pub type ReadPlacements = Vec<Placement>;

/// Tunables for the read-conflict criterion (de-tie). Defaults are the bake-off operating point
/// (`bench/family_criterion_bakeoff.md`); env-overridable via `RUSTLE_CONFLICT_DE_DELTA/DE_MAX/MIN_READS`.
#[derive(Clone, Copy, Debug)]
pub struct ConflictParams {
    /// Two placements conflict iff their divergences are within `delta` AND both `<= de_max` (both fit).
    pub delta: f64,
    pub de_max: f64,
    /// Minimum conflicting-read count for an edge (guards the noise floor).
    pub min_reads: usize,
}

impl Default for ConflictParams {
    fn default() -> Self {
        ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 3 }
    }
}

impl ConflictParams {
    /// Read overrides from `RUSTLE_CONFLICT_DE_DELTA`, `RUSTLE_CONFLICT_DE_MAX`, `RUSTLE_CONFLICT_MIN_READS`.
    pub fn from_env() -> Self {
        let d = Self::default();
        let f = |k: &str, v: f64| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(v);
        let u = |k: &str, v: usize| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(v);
        ConflictParams {
            delta: f("RUSTLE_CONFLICT_DE_DELTA", d.delta),
            de_max: f("RUSTLE_CONFLICT_DE_MAX", d.de_max),
            min_reads: u("RUSTLE_CONFLICT_MIN_READS", d.min_reads),
        }
    }
}

/// de-tie: `|de_a − de_b| <= delta` AND `max(de_a, de_b) <= de_max` (the read fits both copies, comparably).
fn de_tied(a: &Placement, b: &Placement, p: &ConflictParams) -> bool {
    let (da, db) = (a.de as f64, b.de as f64);
    (da - db).abs() <= p.delta && da.max(db) <= p.de_max
}

/// `min(a,b) >= as_tie * max(a,b)` — the legacy AS-tie predicate, kept only for the audit edge-set.
fn as_tied(a: i32, b: i32, as_tie: f64) -> bool {
    let (hi, lo) = (a.max(b), a.min(b));
    hi > 0 && (lo as f64) >= as_tie * (hi as f64)
}

/// Build the read-conflict edges over `n_loci` under the **de-tie** criterion. For each read, every pair of its
/// placements that de-ties contributes one conflict observation to that locus pair. Returns `(i, j, weight)`
/// with `i < j` for pairs whose count `>= p.min_reads`, sorted. Self-pairs (same locus) ignored.
pub fn conflict_edges(n_loci: usize, reads: &[ReadPlacements], p: &ConflictParams) -> Vec<(usize, usize, usize)> {
    use std::collections::BTreeMap;
    let mut weight: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for placements in reads {
        for a in 0..placements.len() {
            for b in (a + 1)..placements.len() {
                let (pa, pb) = (&placements[a], &placements[b]);
                if pa.locus == pb.locus || pa.locus >= n_loci || pb.locus >= n_loci {
                    continue;
                }
                if de_tied(pa, pb, p) {
                    let key = (pa.locus.min(pb.locus), pa.locus.max(pb.locus));
                    *weight.entry(key).or_insert(0) += 1;
                }
            }
        }
    }
    weight.into_iter().filter(|&(_, w)| w >= p.min_reads).map(|((i, j), w)| (i, j, w)).collect()
}

/// AS-tie edge node-pairs over `n_loci` — the legacy criterion, used only to log that the de-tie edge set is a
/// subset of the AS-tie edge set (`de-tie ⊆ AS-tie`, the portable regression invariant from the bake-off).
pub fn as_tie_edges(n_loci: usize, reads: &[ReadPlacements], as_tie: f64, min_reads: usize) -> std::collections::BTreeSet<(usize, usize)> {
    use std::collections::BTreeMap;
    let mut weight: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for placements in reads {
        for a in 0..placements.len() {
            for b in (a + 1)..placements.len() {
                let (pa, pb) = (&placements[a], &placements[b]);
                if pa.locus == pb.locus || pa.locus >= n_loci || pb.locus >= n_loci {
                    continue;
                }
                if as_tied(pa.as_score, pb.as_score, as_tie) {
                    *weight.entry((pa.locus.min(pb.locus), pa.locus.max(pb.locus))).or_insert(0) += 1;
                }
            }
        }
    }
    weight.into_iter().filter(|&(_, w)| w >= min_reads).map(|((i, j), _)| (i, j)).collect()
}

/// Connected-component families over the conflict edges (union-find). Returns components of size `>= 2`
/// (a locus with no conflict needs no resolution — it is not a family), each sorted ascending, the list
/// sorted by first member (deterministic).
pub fn conflict_families(n_loci: usize, edges: &[(usize, usize, usize)]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..n_loci).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            r = parent[r];
        }
        let mut c = x;
        while parent[c] != r {
            let next = parent[c];
            parent[c] = r;
            c = next;
        }
        r
    }
    for &(a, b, _) in edges {
        if a < n_loci && b < n_loci {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            if ra != rb {
                parent[ra.max(rb)] = ra.min(rb);
            }
        }
    }
    let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
    for x in 0..n_loci {
        let r = find(&mut parent, x);
        groups.entry(r).or_default().push(x);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().filter(|g| g.len() >= 2).collect();
    for g in &mut out {
        g.sort_unstable();
    }
    out.sort_by_key(|g| g[0]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    fn p(locus: usize, de: f32) -> Placement { Placement { locus, de, mapq: 0, as_score: 100 } }

    #[test]
    fn de_tied_placements_make_an_edge_and_a_family() {
        // both copies fit comparably (de 0.010 vs 0.012, both < 0.05) -> tie -> edge -> family.
        let reads = vec![vec![p(0, 0.010), p(1, 0.012)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert_eq!(edges, vec![(0, 1, 1)]);
        assert_eq!(conflict_families(2, &edges), vec![vec![0, 1]]);
    }

    #[test]
    fn divergence_gap_beyond_delta_is_not_a_conflict() {
        // read fits copy 0 (de 0.001) far better than copy 1 (de 0.020): |Δ|=0.019 > 0.005 -> resolvable.
        let reads = vec![vec![p(0, 0.001), p(1, 0.020)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert!(edges.is_empty());
    }

    #[test]
    fn both_high_divergence_blocked_by_ceiling() {
        // de_a 0.06 ~ de_b 0.061 (tied within delta) but both exceed de_max 0.05 -> read fits neither.
        let reads = vec![vec![p(0, 0.060), p(1, 0.061)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert!(edges.is_empty());
    }

    #[test]
    fn single_placement_read_is_a_singleton_not_a_family() {
        let reads = vec![vec![p(0, 0.01)], vec![p(1, 0.01)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert!(edges.is_empty());
        assert!(conflict_families(2, &edges).is_empty());
    }

    #[test]
    fn min_reads_threshold_drops_thin_conflicts() {
        let one = vec![vec![p(0, 0.01), p(1, 0.012)]];
        let pr = ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 3 };
        assert!(conflict_edges(2, &one, &pr).is_empty());
        let three = vec![one[0].clone(), one[0].clone(), one[0].clone()];
        assert_eq!(conflict_edges(2, &three, &pr), vec![(0, 1, 3)]);
    }

    #[test]
    fn transitive_conflict_closes_into_one_family() {
        let reads = vec![vec![p(0, 0.010), p(1, 0.012)], vec![p(1, 0.010), p(2, 0.013)]];
        let edges = conflict_edges(3, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert_eq!(conflict_families(3, &edges), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn disjoint_conflicts_form_separate_families() {
        let reads = vec![
            vec![p(0, 0.01), p(1, 0.012)],
            vec![p(2, 0.01), p(3, 0.012)],
            vec![p(4, 0.01)],
        ];
        let edges = conflict_edges(5, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert_eq!(conflict_families(5, &edges), vec![vec![0, 1], vec![2, 3]]);
    }

    #[test]
    fn default_params_are_the_operating_point() {
        let d = ConflictParams::default();
        assert!((d.delta - 0.005).abs() < 1e-9);
        assert!((d.de_max - 0.05).abs() < 1e-9);
        assert_eq!(d.min_reads, 3);
    }

    #[test]
    fn as_tie_edges_superset_of_de_edges() {
        // AS ties two placements that de SPLITS (de 0.001 vs 0.020): AS-edge exists, de-edge does not.
        let reads = vec![vec![
            Placement { locus: 0, de: 0.001, mapq: 0, as_score: 500 },
            Placement { locus: 1, de: 0.020, mapq: 0, as_score: 498 },
        ]];
        let de_edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        let as_edges = as_tie_edges(2, &reads, 0.9, 1);
        assert!(de_edges.is_empty());
        assert_eq!(as_edges, std::collections::BTreeSet::from([(0, 1)]));
    }
}
