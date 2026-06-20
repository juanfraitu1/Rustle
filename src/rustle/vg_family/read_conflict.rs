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

/// One read's placements over the family's candidate loci: `(locus index, alignment score)` for the read's
/// primary + secondary alignments that landed on a candidate locus. Built from the BAM by the adapter.
pub type ReadPlacements = Vec<(usize, i32)>;

/// Tunables for the read-conflict criterion.
#[derive(Clone, Copy, Debug)]
pub struct ConflictParams {
    /// A pair of placements is a conflict iff `min(AS) >= as_tie * max(AS)` (the alignment-score tie that
    /// defines genuine ambiguity). `1.0` = exact tie; `0.9` admits a 10% margin.
    pub as_tie: f64,
    /// Minimum number of conflicting reads for an edge to be kept (guards against a lone spurious secondary).
    pub min_reads: usize,
}

impl Default for ConflictParams {
    fn default() -> Self {
        ConflictParams { as_tie: 0.9, min_reads: 1 }
    }
}

/// `min(a,b) >= as_tie * max(a,b)`, with non-positive scores treated as no-tie.
fn as_tied(a: i32, b: i32, as_tie: f64) -> bool {
    let (hi, lo) = (a.max(b), a.min(b));
    hi > 0 && (lo as f64) >= as_tie * (hi as f64)
}

/// Build the read-conflict edges over `n_loci`. For each read, every pair of its placements whose alignment
/// scores are TIED contributes one conflict observation to that locus pair. Returns `(i, j, weight)` with
/// `i < j` for pairs whose conflict count `>= p.min_reads`, sorted (deterministic). A placement appearing
/// twice on the same locus is ignored (no self-edge). `O(reads * placements^2)`.
pub fn conflict_edges(n_loci: usize, reads: &[ReadPlacements], p: &ConflictParams) -> Vec<(usize, usize, usize)> {
    use std::collections::BTreeMap;
    let mut weight: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for placements in reads {
        for a in 0..placements.len() {
            for b in (a + 1)..placements.len() {
                let (la, sa) = placements[a];
                let (lb, sb) = placements[b];
                if la == lb || la >= n_loci || lb >= n_loci {
                    continue;
                }
                if as_tied(sa, sb, p.as_tie) {
                    let key = (la.min(lb), la.max(lb));
                    *weight.entry(key).or_insert(0) += 1;
                }
            }
        }
    }
    weight
        .into_iter()
        .filter(|&(_, w)| w >= p.min_reads)
        .map(|((i, j), w)| (i, j, w))
        .collect()
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

    #[test]
    fn tied_placements_make_an_edge_and_a_family() {
        // one read cross-maps loci 0 and 1 with tied scores -> one edge -> one 2-member family.
        let reads = vec![vec![(0usize, 100i32), (1usize, 98i32)]];
        let edges = conflict_edges(2, &reads, &ConflictParams::default());
        assert_eq!(edges, vec![(0, 1, 1)]);
        assert_eq!(conflict_families(2, &edges), vec![vec![0, 1]]);
    }

    #[test]
    fn untied_placement_is_not_a_conflict() {
        // the secondary is far below the primary (read resolves uniquely) -> no edge -> no family.
        let reads = vec![vec![(0usize, 100i32), (1usize, 50i32)]];
        let edges = conflict_edges(2, &reads, &ConflictParams::default());
        assert!(edges.is_empty());
        assert!(conflict_families(2, &edges).is_empty());
    }

    #[test]
    fn single_placement_read_is_a_singleton_not_a_family() {
        // a read over a shared exon maps to ONE locus (the domain-sharer case) -> no alternative placement.
        let reads = vec![vec![(0usize, 100i32)], vec![(1usize, 100i32)]];
        let edges = conflict_edges(2, &reads, &ConflictParams::default());
        assert!(edges.is_empty());
        assert!(conflict_families(2, &edges).is_empty());
    }

    #[test]
    fn min_reads_threshold_drops_a_lone_conflict() {
        let reads = vec![vec![(0usize, 100i32), (1usize, 99i32)]];
        let p = ConflictParams { as_tie: 0.9, min_reads: 2 };
        assert!(conflict_edges(2, &reads, &p).is_empty());
        // two conflicting reads clear it.
        let reads2 = vec![reads[0].clone(), vec![(0usize, 100i32), (1usize, 99i32)]];
        assert_eq!(conflict_edges(2, &reads2, &p), vec![(0, 1, 2)]);
    }

    #[test]
    fn transitive_conflict_closes_into_one_family() {
        // A~B tied and B~C tied (no read directly links A~C) -> all three are ONE family (the closure).
        let reads = vec![vec![(0usize, 100i32), (1usize, 98i32)], vec![(1usize, 100i32), (2usize, 97i32)]];
        let edges = conflict_edges(3, &reads, &ConflictParams::default());
        assert_eq!(conflict_families(3, &edges), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn disjoint_conflicts_form_separate_families() {
        // {0,1} and {2,3} cross-map internally but never to each other -> two families, 4 is a singleton.
        let reads = vec![
            vec![(0usize, 100i32), (1usize, 99i32)],
            vec![(2usize, 100i32), (3usize, 99i32)],
            vec![(4usize, 100i32)],
        ];
        let edges = conflict_edges(5, &reads, &ConflictParams::default());
        assert_eq!(conflict_families(5, &edges), vec![vec![0, 1], vec![2, 3]]);
    }

    #[test]
    fn deterministic_under_placement_order() {
        let r1 = vec![vec![(2usize, 100i32), (0usize, 98i32)], vec![(0usize, 100i32), (1usize, 99i32)]];
        let r2 = vec![vec![(0usize, 98i32), (2usize, 100i32)], vec![(1usize, 99i32), (0usize, 100i32)]];
        let p = ConflictParams::default();
        assert_eq!(conflict_families(3, &conflict_edges(3, &r1, &p)),
                   conflict_families(3, &conflict_edges(3, &r2, &p)));
    }
}
