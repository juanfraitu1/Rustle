//! Read-conflict (mutual-mappability) graph — the operational SCOPE for copy-resolution and the SEED for the
//! family catalog. NOT a rival definition of a family (`bench/family_def_readconflict.md`).
//!
//! DEFINITION vs SCOPE. The family DEFINITION is homology: a γ-quasi-clique component of the transcribed-
//! homology graph over ≥ 2 loci (`family_definition.rs`, the E_r oracle; `docs/RETIREMENT_AND_MIGRATION.md`,
//! `DEFINITIONS_FORMAL.md`). This module builds a DIFFERENT graph — over the SAME loci but with a different
//! edge — that answers the question copy-assignment actually cares about: **do reads cross-map between these
//! loci?** Two loci are linked iff some read has a placement in BOTH with TIED alignment scores (a genuine
//! alternative placement — the multimapping conflict). Homology says which loci ARE one family; read-conflict
//! says how many COPIES they hide and how the assignment decomposes. So it is the SCOPE (which loci must be
//! co-resolved) and the SEED for the catalog, not the definition of membership.
//!
//! Why this is the right SCOPE: (1) no tuned similarity threshold — the boundary is the alignment-score tie
//! (with `RUSTLE_CONFLICT_SIG`, the SAME significance level α the assignment gate uses), a property of the
//! data; (2) reads never cross-map outside their component, so the assignment problem decomposes EXACTLY
//! across families with no information lost; (3) it never groups domain-sharers (a read over a shared exon
//! maps to one locus — no alternative placement), and it picks out exactly the families where assignment is
//! needed (validated on Compara labels: 0 conflict on 7/7 domain-sharers, fires on RABL2/APOBEC3, silent on
//! the resolvable RFPL).
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
    /// Aligned-block length (number of M/=/X columns) of this placement. Used only by the SIGNIFICANCE
    /// de-tie criterion to convert the divergence rate `de` into a mismatch COUNT (`m = de * aln_len`);
    /// the legacy delta-based `de_tied` ignores it, so existing behaviour is unaffected when sig is off.
    pub aln_len: u32,
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
    /// SIGNIFICANCE de-tie (the UNIFICATION with the copy-assignment gate): when `Some((eps, alpha))`, a read
    /// de-ties between two loci iff it CANNOT significantly distinguish them under the SAME IsoCon real-vs-error
    /// test the assignment gate uses — `eps^delta >= alpha`, where `delta = |m_a - m_b|` is the excess mismatch
    /// count (the per-read distinguishing-column proxy, mirroring Theorem 4's `min_p = eps^delta`). The
    /// `de_max` quality floor still applies. `None` (default) = the legacy fixed-`delta` `de_tied`, so OFF is
    /// byte-identical. Replaces the arbitrary `delta=0.005` with the error-model-derived tie threshold.
    pub sig: Option<(f64, f64)>,
}

impl Default for ConflictParams {
    fn default() -> Self {
        ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 3, sig: None }
    }
}

impl ConflictParams {
    /// Read overrides from `RUSTLE_CONFLICT_DE_DELTA`, `RUSTLE_CONFLICT_DE_MAX`, `RUSTLE_CONFLICT_MIN_READS`.
    pub fn from_env() -> Self {
        let d = Self::default();
        let f = |k: &str, v: f64| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(v);
        let u = |k: &str, v: usize| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(v);
        // Significance de-tie is DEFAULT ON: the conflict edge uses the SAME IsoCon real-vs-error test (and the
        // SAME level `alpha`) as the assignment gate, so no hand-set score-gap `delta` decides the conflict
        // scope. eps = per-distinguishing-column error proxy (e/3 ~ 0.001 HiFi), alpha = significance.
        // RUSTLE_CONFLICT_SIG=0 reverts to the legacy `delta` tie-width (A/B comparison / legacy reproduction).
        let sig = if std::env::var("RUSTLE_CONFLICT_SIG").ok().as_deref() == Some("0") {
            None
        } else {
            Some((f("RUSTLE_CONFLICT_EPS", 0.001), f("RUSTLE_CONFLICT_ALPHA", 1e-3)))
        };
        ConflictParams {
            delta: f("RUSTLE_CONFLICT_DE_DELTA", d.delta),
            de_max: f("RUSTLE_CONFLICT_DE_MAX", d.de_max),
            min_reads: u("RUSTLE_CONFLICT_MIN_READS", d.min_reads),
            sig,
        }
    }
}

/// de-tie: `|de_a − de_b| <= delta` AND `max(de_a, de_b) <= de_max` (the read fits both copies, comparably).
fn de_tied(a: &Placement, b: &Placement, p: &ConflictParams) -> bool {
    let (da, db) = (a.de as f64, b.de as f64);
    (da - db).abs() <= p.delta && da.max(db) <= p.de_max
}

/// SIGNIFICANCE de-tie — the unification with the copy-assignment gate. A read counts as conflict evidence
/// between two loci iff it CANNOT significantly distinguish them: with `m_x = de_x * aln_len_x` the mismatch
/// count to locus `x`, the excess `delta = |m_a - m_b|` is the per-read distinguishing-column proxy, and the
/// read is tied iff `eps^delta >= alpha` — exactly the assignment gate's `min_p >= alpha` (Theorem 4). The
/// `de_max` quality floor still applies (both alignments must genuinely fit). No arbitrary `delta` constant.
fn sig_tied(a: &Placement, b: &Placement, de_max: f64, eps: f64, alpha: f64) -> bool {
    let (da, db) = (a.de as f64, b.de as f64);
    if da.max(db) > de_max {
        return false;
    }
    let ma = da * a.aln_len as f64;
    let mb = db * b.aln_len as f64;
    let delta_cols = (ma - mb).abs();
    eps.powf(delta_cols) >= alpha
}

/// Whether a read's two placements conflict (de-tie), under either criterion (significance if `p.sig` is set).
fn tied(a: &Placement, b: &Placement, p: &ConflictParams) -> bool {
    match p.sig {
        Some((eps, alpha)) => sig_tied(a, b, p.de_max, eps, alpha),
        None => de_tied(a, b, p),
    }
}

/// A read's alignment-score evidence: the best and runner-up `AS:i` over its placements, both raw and
/// normalized by aligned length.
///
/// AS is REPORTED, never decisive — `de` decides. Raw AS is an absolute score that grows with aligned
/// length, so a genuine multimapper whose second placement is a PARTIAL alignment scores far lower there
/// even at identical per-base quality. Measured on GGO Iso-Seq: median runner-up/best AS ratio 0.713 while
/// the aligned-length ratio is 0.897. That length confound is why `de` (a rate) replaced AS (a total), and
/// why `de-tie ⊆ AS-tie` does NOT hold on real data. `per_base` divides out the confound and is the value
/// to compare across placements.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AsEvidence {
    pub best: i32,
    pub second: Option<i32>,
    pub best_per_base: f32,
    pub second_per_base: Option<f32>,
}

impl AsEvidence {
    /// `best - second`; `None` when the read has a single placement (nothing to be ambiguous between).
    pub fn margin(&self) -> Option<i32> {
        self.second.map(|s| self.best - s)
    }
}

/// Best and runner-up alignment score over one read's placements, as `(as_score, aligned_len)` pairs.
/// Ranks on RAW `AS` (that is the quantity people quote and Eichler's `AS >= 10` rule uses); the per-base
/// values are carried alongside for the length-fair comparison. `None` if the read has no placements.
/// Zero-length placements get a per-base score of 0.0 rather than a division by zero.
pub fn as_evidence(placements: &[(i32, u32)]) -> Option<AsEvidence> {
    let per_base = |(s, l): (i32, u32)| if l == 0 { 0.0 } else { s as f32 / l as f32 };
    let mut sorted: Vec<(i32, u32)> = placements.to_vec();
    sorted.sort_by(|a, b| b.0.cmp(&a.0));
    let &first = sorted.first()?;
    let second = sorted.get(1).copied();
    Some(AsEvidence {
        best: first.0,
        second: second.map(|s| s.0),
        best_per_base: per_base(first),
        second_per_base: second.map(per_base),
    })
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
                if tied(pa, pb, p) {
                    let key = (pa.locus.min(pb.locus), pa.locus.max(pb.locus));
                    *weight.entry(key).or_insert(0) += 1;
                }
            }
        }
    }
    weight.into_iter().filter(|&(_, w)| w >= p.min_reads).map(|((i, j), w)| (i, j, w)).collect()
}

/// AS-tie edge node-pairs over `n_loci` — the legacy criterion, kept only as a logged comparison edge-set.
///
/// ⚠ `de-tie ⊆ AS-tie` holds only when both placements are FULL-LENGTH (the unit tests and the planted sims).
/// On real Iso-Seq it is FALSE: secondary placements are partial, raw AS scales with aligned length, so a
/// de-tied pair can fail `as_tied` outright. Every real GGO region logs `de⊆AS=false`. Do not treat the
/// audit line as a regression invariant on real data — it is a diagnostic.
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

/// Count reads supporting a conflict family and how many of those reads have BOTH placements at mapq==0
/// (the genuine-multimapper corroboration). A read is counted when it contributes at least one de-tied pair
/// whose two loci are BOTH in `family`; of those, a `both_mapq0` read has mapq==0 on BOTH placements in
/// that pair. Returns `(supporting_reads, both_mapq0_reads)`. Log-only: does NOT gate any edge.
///
/// For a read with placements on >=3 family loci the mapq0 check is applied to the FIRST qualifying pair by
/// iteration order, so the `both_mapq0` count is CONSERVATIVE for multi-copy families (it can undercount but
/// never overstate multimapper evidence) — exact for the common 2-locus family.
pub fn family_mapq0_support(reads: &[ReadPlacements], family: &[usize], p: &ConflictParams) -> (usize, usize) {
    let fset: std::collections::BTreeSet<usize> = family.iter().copied().collect();
    let mut support = 0usize;
    let mut both_mapq0 = 0usize;
    'read: for placements in reads {
        // Scan every pair within this read; stop at the first pair that fires (count once per read).
        for ai in 0..placements.len() {
            for bi in (ai + 1)..placements.len() {
                let (pa, pb) = (&placements[ai], &placements[bi]);
                if !fset.contains(&pa.locus) || !fset.contains(&pb.locus) {
                    continue;
                }
                if pa.locus == pb.locus {
                    continue;
                }
                if tied(pa, pb, p) {
                    support += 1;
                    if pa.mapq == 0 && pb.mapq == 0 {
                        both_mapq0 += 1;
                    }
                    continue 'read;
                }
            }
        }
    }
    (support, both_mapq0)
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
    fn p(locus: usize, de: f32) -> Placement { Placement { locus, de, mapq: 0, as_score: 100, aln_len: 2000 } }

    #[test]
    fn de_tied_placements_make_an_edge_and_a_family() {
        // both copies fit comparably (de 0.010 vs 0.012, both < 0.05) -> tie -> edge -> family.
        let reads = vec![vec![p(0, 0.010), p(1, 0.012)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert_eq!(edges, vec![(0, 1, 1)]);
        assert_eq!(conflict_families(2, &edges), vec![vec![0, 1]]);
    }

    #[test]
    fn sig_criterion_ties_ambiguous_resolves_distinguishing() {
        // The UNIFICATION: significance edge mirrors the assignment gate's min_p>=alpha (eps^delta>=alpha).
        let sig = ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: Some((1e-3, 1e-3)) };
        // tau = floor(ln(alpha)/ln(eps)) = floor(ln(1e-3)/ln(1e-3)) = 1: tied iff excess-mismatches <= 1.
        // AMBIGUOUS read: equal divergence on both copies (delta_cols = 0) -> tied -> edge.
        let amb = vec![vec![
            Placement { locus: 0, de: 0.0050, mapq: 0, as_score: 100, aln_len: 2000 },
            Placement { locus: 1, de: 0.0050, mapq: 0, as_score: 100, aln_len: 2000 },
        ]];
        assert_eq!(conflict_edges(2, &amb, &sig), vec![(0, 1, 1)]);
        // DISTINGUISHING read: m_a=0.0005*2000=1, m_b=0.0050*2000=10 -> delta_cols=9 >> tau -> NO edge
        // (the boundary is now eps^delta>=alpha, not a hand-set 0.005).
        let dist = vec![vec![
            Placement { locus: 0, de: 0.0005, mapq: 0, as_score: 100, aln_len: 2000 },
            Placement { locus: 1, de: 0.0050, mapq: 0, as_score: 100, aln_len: 2000 },
        ]];
        assert!(conflict_edges(2, &dist, &sig).is_empty());
    }

    #[test]
    fn sig_edge_is_a_refinement_of_de_tied_equal_length() {
        // Exhaustive over a divergence grid (equal aligned length, default eps/alpha): every SIG-tie is also a
        // de-tie, so the significance edge set is a SUBSET of the de-tie edge set -> SIG can only shrink/split
        // families, never invent them (the rigorous refinement claim behind the 81->71 catalog narrowing).
        let de_p = ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None };
        let (eps, alpha, l) = (1e-3f64, 1e-3f64, 2000u32);
        for ia in 0..=60u32 {
            for ib in 0..=60u32 {
                let (da, db) = (ia as f32 * 0.001, ib as f32 * 0.001);
                let a = Placement { locus: 0, de: da, mapq: 0, as_score: 0, aln_len: l };
                let b = Placement { locus: 1, de: db, mapq: 0, as_score: 0, aln_len: l };
                if sig_tied(&a, &b, de_p.de_max, eps, alpha) {
                    assert!(de_tied(&a, &b, &de_p), "sig-tie not a de-tie at de=({da},{db})");
                }
            }
        }
    }

    #[test]
    fn sig_off_default_is_byte_identical_to_de_tied() {
        // With sig: None (default), `tied` == `de_tied` exactly: same edges as the legacy criterion.
        let reads = vec![
            vec![p(0, 0.010), p(1, 0.012)],            // tied under delta=0.005
            vec![p(0, 0.001), p(1, 0.020)],            // resolved under delta=0.005
        ];
        let legacy = ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None };
        assert_eq!(conflict_edges(2, &reads, &legacy), vec![(0, 1, 1)]);   // only the first read ties
        assert_eq!(ConflictParams::default().sig, None);                    // default ships OFF
    }

    #[test]
    fn divergence_gap_beyond_delta_is_not_a_conflict() {
        // read fits copy 0 (de 0.001) far better than copy 1 (de 0.020): |Δ|=0.019 > 0.005 -> resolvable.
        let reads = vec![vec![p(0, 0.001), p(1, 0.020)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert!(edges.is_empty());
    }

    #[test]
    fn both_high_divergence_blocked_by_ceiling() {
        // de_a 0.06 ~ de_b 0.061 (tied within delta) but both exceed de_max 0.05 -> read fits neither.
        let reads = vec![vec![p(0, 0.060), p(1, 0.061)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert!(edges.is_empty());
    }

    #[test]
    fn single_placement_read_is_a_singleton_not_a_family() {
        let reads = vec![vec![p(0, 0.01)], vec![p(1, 0.01)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert!(edges.is_empty());
        assert!(conflict_families(2, &edges).is_empty());
    }

    #[test]
    fn min_reads_threshold_drops_thin_conflicts() {
        let one = vec![vec![p(0, 0.01), p(1, 0.012)]];
        let pr = ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 3, sig: None };
        assert!(conflict_edges(2, &one, &pr).is_empty());
        let three = vec![one[0].clone(), one[0].clone(), one[0].clone()];
        assert_eq!(conflict_edges(2, &three, &pr), vec![(0, 1, 3)]);
    }

    #[test]
    fn transitive_conflict_closes_into_one_family() {
        let reads = vec![vec![p(0, 0.010), p(1, 0.012)], vec![p(1, 0.010), p(2, 0.013)]];
        let edges = conflict_edges(3, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert_eq!(conflict_families(3, &edges), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn disjoint_conflicts_form_separate_families() {
        let reads = vec![
            vec![p(0, 0.01), p(1, 0.012)],
            vec![p(2, 0.01), p(3, 0.012)],
            vec![p(4, 0.01)],
        ];
        let edges = conflict_edges(5, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
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
    fn deterministic_under_placement_order() {
        // Shuffling placement order within a read and across reads must not change the edge/family output.
        let pr = ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None };
        // two reads, each with placements on loci {0,1,2} in different orders.
        let forward = vec![
            vec![p(0, 0.010), p(1, 0.012), p(2, 0.030)],
            vec![p(1, 0.011), p(0, 0.013)],
        ];
        let reversed = vec![
            vec![p(2, 0.030), p(1, 0.012), p(0, 0.010)],
            vec![p(0, 0.013), p(1, 0.011)],
        ];
        let reads_swapped = vec![forward[1].clone(), forward[0].clone()];
        let e_fwd = conflict_edges(3, &forward, &pr);
        let e_rev = conflict_edges(3, &reversed, &pr);
        let e_swp = conflict_edges(3, &reads_swapped, &pr);
        assert_eq!(e_fwd, e_rev, "placement order within reads must not change edges");
        assert_eq!(e_fwd, e_swp, "read order must not change edges");
        assert_eq!(conflict_families(3, &e_fwd), conflict_families(3, &e_rev));
        assert_eq!(conflict_families(3, &e_fwd), conflict_families(3, &e_swp));
    }

    #[test]
    fn de_max_boundary_exactly_at_threshold_fires_just_over_is_blocked() {
        // Safely inside: de 0.049 vs 0.049 — max(de)=0.049 < 0.05, |Δ|=0 ≤ 0.005 → should fire.
        let pr = ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None };
        let inside = vec![vec![p(0, 0.049), p(1, 0.049)]];
        let e_inside = conflict_edges(2, &inside, &pr);
        assert_eq!(e_inside, vec![(0, 1, 1)], "de=0.049 <= de_max=0.05 must fire");
        // NOTE: exactly 0.05f32 widens above 0.05f64 after f32→f64 cast → on the precision boundary;
        // we test the clearly-over case (0.051) which is blocked regardless.
        let over = vec![vec![p(0, 0.051), p(1, 0.051)]];
        let e_over = conflict_edges(2, &over, &pr);
        assert!(e_over.is_empty(), "de=0.051 > de_max=0.05 must be blocked");
    }

    #[test]
    fn family_mapq0_support_counts_correctly() {
        let pr = ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None };
        let family = vec![0usize, 1];
        // read A: both loci de-tied, both mapq==0 → counts in support AND both_mapq0.
        let read_a: ReadPlacements = vec![
            Placement { locus: 0, de: 0.010, mapq: 0, as_score: 100, aln_len: 2000 },
            Placement { locus: 1, de: 0.012, mapq: 0, as_score: 100, aln_len: 2000 },
        ];
        // read B: both loci de-tied, but mapq>0 on one → counts in support but NOT both_mapq0.
        let read_b: ReadPlacements = vec![
            Placement { locus: 0, de: 0.010, mapq: 60, as_score: 100, aln_len: 2000 },
            Placement { locus: 1, de: 0.012, mapq: 0, as_score: 100, aln_len: 2000 },
        ];
        // read C: de NOT tied (gap too large) → not counted at all.
        let read_c: ReadPlacements = vec![
            Placement { locus: 0, de: 0.001, mapq: 0, as_score: 100, aln_len: 2000 },
            Placement { locus: 1, de: 0.020, mapq: 0, as_score: 100, aln_len: 2000 },
        ];
        let reads = vec![read_a, read_b, read_c];
        let (support, mapq0) = family_mapq0_support(&reads, &family, &pr);
        assert_eq!(support, 2, "read_a and read_b both contribute a de-tied pair in the family");
        assert_eq!(mapq0, 1, "only read_a has both placements mapq==0");
    }

    #[test]
    fn as_tie_edges_superset_of_de_edges() {
        // AS ties two placements that de SPLITS (de 0.001 vs 0.020): AS-edge exists, de-edge does not.
        let reads = vec![vec![
            Placement { locus: 0, de: 0.001, mapq: 0, as_score: 500, aln_len: 2000 },
            Placement { locus: 1, de: 0.020, mapq: 0, as_score: 498, aln_len: 2000 },
        ]];
        let de_edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        let as_edges = as_tie_edges(2, &reads, 0.9, 1);
        assert!(de_edges.is_empty());
        assert_eq!(as_edges, std::collections::BTreeSet::from([(0, 1)]));
    }

    #[test]
    fn as_evidence_none_without_placements() {
        assert_eq!(as_evidence(&[]), None);
    }

    #[test]
    fn as_evidence_single_placement_has_no_second_or_margin() {
        let e = as_evidence(&[(1842, 916)]).unwrap();
        assert_eq!(e.best, 1842);
        assert_eq!(e.second, None);
        assert_eq!(e.second_per_base, None);
        assert_eq!(e.margin(), None, "one placement => nothing to be ambiguous between");
    }

    #[test]
    fn as_evidence_ranks_best_and_runner_up_over_three_placements() {
        let e = as_evidence(&[(1310, 800), (1842, 916), (900, 500)]).unwrap();
        assert_eq!(e.best, 1842);
        assert_eq!(e.second, Some(1310));
        assert_eq!(e.margin(), Some(532));
        assert!((e.best_per_base - 1842.0 / 916.0).abs() < 1e-6);
        assert!((e.second_per_base.unwrap() - 1310.0 / 800.0).abs() < 1e-6);
    }

    #[test]
    fn as_evidence_zero_length_placement_does_not_divide_by_zero() {
        let e = as_evidence(&[(50, 0)]).unwrap();
        assert_eq!(e.best_per_base, 0.0);
    }

    /// The measured real-data confound, pinned: two placements of EQUAL per-base quality where the runner-up
    /// is a partial alignment. Raw AS says "not a tie" (ratio 0.60 < 0.90) while per-base AS says they are
    /// indistinguishable. This is why AS is reported and `de` decides.
    #[test]
    fn raw_as_misjudges_a_partial_placement_that_per_base_as_calls_equal() {
        let (full, partial) = ((2000, 1000), (1200, 600)); // both exactly 2.0 AS per aligned base
        let e = as_evidence(&[full, partial]).unwrap();
        assert_eq!(e.margin(), Some(800), "raw AS shows a large margin...");
        assert!(!as_tied(full.0, partial.0, 0.9), "...so the legacy AS-tie predicate rejects the pair");
        assert_eq!(
            e.best_per_base,
            e.second_per_base.unwrap(),
            "...yet per-aligned-base they are identical: the margin is pure length confound"
        );
    }
}
