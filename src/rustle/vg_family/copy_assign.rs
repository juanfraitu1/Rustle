//! Copy ASSIGNMENT (resolution): assign ONE read to a specific paralog COPY, given the family's known
//! copies, using PSV bases + copy-specific JUNCTIONS, behind an IDENTIFIABILITY gate.
//!
//! This is distinct from `copy_split` (which DISCOVERS the copies). Here the copy set is known (from the
//! family layer) and we RESOLVE which copy a read came from -- the "pick one mapping" step for hard
//! multimappers minimap2 leaves at MAPQ 0. A read carries its true copy's alleles regardless of which
//! paralog it aligned to, so matching its observed features against every copy's feature vector resolves it.
//!
//! Identifiability theorem (operational): a read is RESOLVABLE iff it spans >= 1 feature (a PSV column or
//! a junction boundary) where the candidate copies DIFFER. Reads spanning 0 such features are genuinely
//! TIED ("N equally good places"). Among resolvable reads, a log-likelihood-ratio margin over the
//! runner-up calls ASSIGNED vs AMBIGUOUS. Mirrors `bench/copy_assign.py::assign_read`.

use super::copy_split::{allele_at, intron_chain_of, AlignedRead};

/// Per-copy feature profile over the family's PSV columns + intron-boundary set.
#[derive(Clone, Debug)]
pub struct CopyProfile {
    pub copy_id: usize,
    /// Allele base at each family PSV column (index = column). `None` = this copy has a gap there.
    /// A copy that MATCHES the reference allele must carry that base here, not `None` (else the
    /// likelihood would favour copies it shares no columns with).
    pub alleles: Vec<Option<u8>>,
    /// Intron-boundary offsets in transcription-spliced space (a copy's copy-specific junctions).
    pub junctions: Vec<i64>,
}

/// One read's observed features in the family's column/boundary space.
#[derive(Clone, Debug)]
pub struct ReadFeatures {
    /// Observed base per PSV column (`None` = the read does not span that column).
    pub psv_obs: Vec<Option<u8>>,
    /// The read's intron-boundary offsets (transcription-spliced space).
    pub junctions: Vec<i64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AssignStatus {
    /// Resolvable and the best copy clears the likelihood-ratio margin.
    Assigned,
    /// Resolvable (spans a decisive feature) but the margin over the runner-up is too small.
    Ambiguous,
    /// Spans no feature where the copies differ -> genuinely unassignable.
    Tied,
}

#[derive(Clone, Debug)]
pub struct Assignment {
    pub best_copy: usize,
    pub log_lr_margin: f64,
    /// Number of features (PSV columns + junction boundaries) the read observes where copies differ.
    pub n_decisive: usize,
    pub resolvable: bool,
    pub status: AssignStatus,
}

/// Tunable assignment parameters (defaults mirror the Python prototype: HiFi error, jw=5, tol=4, margin=2).
#[derive(Clone, Copy, Debug)]
pub struct AssignParams {
    /// Per-base PSV error rate used in the likelihood (HiFi ~ 0.003).
    pub error_rate: f64,
    /// Per-junction log-odds: a matching copy-specific junction is +this, a non-match is -this.
    pub junction_weight: f64,
    /// Junction-boundary match tolerance (bp) for splice-site jitter.
    pub boundary_tol: i64,
    /// Minimum log-likelihood-ratio over the runner-up to call ASSIGNED.
    pub margin: f64,
}

impl Default for AssignParams {
    fn default() -> Self {
        AssignParams { error_rate: 0.003, junction_weight: 5.0, boundary_tol: 4, margin: 2.0 }
    }
}

fn boundary_present(jb: i64, junctions: &[i64], tol: i64) -> bool {
    junctions.iter().any(|&x| (jb - x).abs() <= tol)
}

/// Assign a read to its most likely copy. Returns `None` only if `copies` is empty.
/// Deterministic: copies are scored in slice order and ties resolve to the earliest (lowest index).
pub fn assign_read(read: &ReadFeatures, copies: &[CopyProfile], p: &AssignParams) -> Option<Assignment> {
    let n = copies.len();
    if n == 0 {
        return None;
    }
    let lp_match = (1.0 - p.error_rate).ln();
    let lp_mis = (p.error_rate / 3.0).ln();
    let mut logl = vec![0.0f64; n];
    let mut n_decisive = 0usize;

    // --- PSV term ---
    let n_cols = read.psv_obs.len();
    for j in 0..n_cols {
        let obs = match read.psv_obs[j] {
            Some(b) => b,
            None => continue, // read does not span this column
        };
        // decisive iff the copies disagree at this column (among those with a defined allele)
        let mut seen: Option<u8> = None;
        let mut differ = false;
        for c in copies {
            if let Some(a) = c.alleles.get(j).copied().flatten() {
                match seen {
                    None => seen = Some(a),
                    Some(s) => {
                        if s != a {
                            differ = true;
                        }
                    }
                }
            }
        }
        if differ {
            n_decisive += 1;
        }
        for (ci, c) in copies.iter().enumerate() {
            if let Some(a) = c.alleles.get(j).copied().flatten() {
                logl[ci] += if obs == a { lp_match } else { lp_mis };
            }
        }
    }

    // --- junction (copy-specific splicing) term ---
    for &jb in &read.junctions {
        let present: Vec<bool> = copies
            .iter()
            .map(|c| boundary_present(jb, &c.junctions, p.boundary_tol))
            .collect();
        let np = present.iter().filter(|&&x| x).count();
        if np > 0 && np < n {
            n_decisive += 1; // some copies have this junction, others lack it
        }
        for ci in 0..n {
            logl[ci] += if present[ci] { p.junction_weight } else { -p.junction_weight };
        }
    }

    // argmax (earliest on ties) + runner-up
    let mut best = 0usize;
    for i in 1..n {
        if logl[i] > logl[best] {
            best = i;
        }
    }
    let mut second = f64::NEG_INFINITY;
    for i in 0..n {
        if i != best && logl[i] > second {
            second = logl[i];
        }
    }
    let margin = if n > 1 { logl[best] - second } else { f64::INFINITY };
    let resolvable = n_decisive >= 1;
    let status = if !resolvable {
        AssignStatus::Tied
    } else if margin >= p.margin {
        AssignStatus::Assigned
    } else {
        AssignStatus::Ambiguous
    };
    Some(Assignment { best_copy: copies[best].copy_id, log_lr_margin: margin, n_decisive, resolvable, status })
}

/// Build the PSV-observation vector for a read aligned in the frame of `psv_positions` (genomic, 0-based,
/// in the mapped copy's coordinate frame). Reuses the `copy_split::allele_at` CIGAR bridge. Bases are
/// forward-genome; the caller reverse-complements for a minus-strand copy before comparing (so allele
/// space matches the copies' transcription-strand alleles). Junction boundaries are computed by the
/// integration layer (they need per-copy exon context) and supplied separately.
pub fn read_psv_obs(read: &AlignedRead, psv_positions: &[u64]) -> Vec<Option<u8>> {
    psv_positions.iter().map(|&p| allele_at(read, p)).collect()
}

/// The read's intron donors/acceptors (genomic), straight from its CIGAR N ops -- the raw input the
/// integration layer maps into transcription-spliced boundary offsets.
pub fn read_introns(read: &AlignedRead) -> Vec<(u64, u64)> {
    intron_chain_of(read)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cp(copy_id: usize, alleles: &[Option<u8>], junctions: &[i64]) -> CopyProfile {
        CopyProfile { copy_id, alleles: alleles.to_vec(), junctions: junctions.to_vec() }
    }
    fn rf(psv: &[Option<u8>], junctions: &[i64]) -> ReadFeatures {
        ReadFeatures { psv_obs: psv.to_vec(), junctions: junctions.to_vec() }
    }
    const A: u8 = b'A';
    const C: u8 = b'C';
    const G: u8 = b'G';
    const T: u8 = b'T';

    #[test]
    fn two_copies_one_psv_resolves() {
        let copies = [cp(0, &[Some(A)], &[]), cp(1, &[Some(C)], &[])];
        let r = rf(&[Some(A)], &[]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 0);
        assert!(a.resolvable);
        assert_eq!(a.n_decisive, 1);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn no_features_is_tied() {
        let copies = [cp(0, &[Some(A)], &[100]), cp(1, &[Some(C)], &[200])];
        let r = rf(&[None], &[]); // spans no PSV column, no junctions
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert!(!a.resolvable);
        assert_eq!(a.n_decisive, 0);
        assert_eq!(a.status, AssignStatus::Tied);
    }

    #[test]
    fn shared_allele_not_decisive() {
        // both copies carry A at the only spanned column -> not decisive -> tied
        let copies = [cp(0, &[Some(A)], &[]), cp(1, &[Some(A)], &[])];
        let r = rf(&[Some(A)], &[]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.n_decisive, 0);
        assert_eq!(a.status, AssignStatus::Tied);
    }

    #[test]
    fn junction_only_resolves() {
        // no PSV spanned; copy0 has a copy-specific junction at 100, copy1 does not
        let copies = [cp(0, &[None], &[100]), cp(1, &[None], &[500])];
        let r = rf(&[None], &[100]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 0);
        assert!(a.resolvable);
        assert_eq!(a.n_decisive, 1);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn psv_error_tolerated() {
        // read matches copy0 at 2 of 3 columns (1 sporadic error) -> still copy0, assigned
        let copies = [
            cp(0, &[Some(A), Some(A), Some(A)], &[]),
            cp(1, &[Some(C), Some(C), Some(C)], &[]),
        ];
        let r = rf(&[Some(A), Some(C), Some(A)], &[]); // middle base is an error
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 0);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn five_copies_two_psv_unique() {
        // sim5x K=2 analogue: base-4 allele vectors over 2 columns give 5 distinct copies.
        let bases = [A, C, G, T];
        let allele = |c: usize, j: usize| bases[(c / 4usize.pow(j as u32)) % 4];
        let copies: Vec<CopyProfile> = (0..5)
            .map(|c| cp(c, &[Some(allele(c, 0)), Some(allele(c, 1))], &[]))
            .collect();
        // a read carrying copy 3's alleles, spanning both columns -> uniquely copy 3
        let r = rf(&[Some(allele(3, 0)), Some(allele(3, 1))], &[]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 3);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn five_copies_one_column_ambiguous_for_collision() {
        // K=1 collision: copies 0 and 4 share allele A at column 0 (base-4: 0%4==0, 4%4==0).
        let bases = [A, C, G, T];
        let allele = |c: usize, j: usize| bases[(c / 4usize.pow(j as u32)) % 4];
        let copies: Vec<CopyProfile> = (0..5)
            .map(|c| cp(c, &[Some(allele(c, 0))], &[])) // only column 0
            .collect();
        // a read from copy 4 spanning only column 0 -> allele A matches copy0 AND copy4 -> ambiguous
        let r = rf(&[Some(allele(4, 0))], &[]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert!(a.resolvable); // the column IS decisive (copies 1,2,3 differ)
        assert_eq!(a.status, AssignStatus::Ambiguous); // but copy0/copy4 tie -> margin 0 < 2
    }

    #[test]
    fn junction_augments_decisive_count() {
        // 1 decisive PSV column + 1 decisive junction -> n_decisive == 2
        let copies = [cp(0, &[Some(A)], &[100]), cp(1, &[Some(C)], &[500])];
        let r = rf(&[Some(A)], &[100]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.n_decisive, 2);
        assert_eq!(a.best_copy, 0);
    }

    #[test]
    fn boundary_tolerance_matches_jitter() {
        // read junction at 102 vs copy junction at 100, tol=4 -> present
        let copies = [cp(0, &[None], &[100]), cp(1, &[None], &[900])];
        let r = rf(&[None], &[102]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 0);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn copy_matching_ref_allele_inherits_it() {
        // the Python bug guard: a copy that matches the ref allele must carry that base, not None.
        // copy0=(A,A), copy1=(C,A): a read (A,C) must NOT be mis-assigned to a copy with a None at col1.
        let copies = [cp(0, &[Some(A), Some(A)], &[]), cp(1, &[Some(A), Some(C)], &[])];
        let r = rf(&[Some(A), Some(C)], &[]); // matches copy1 at col1 (the decisive one)
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 1);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn deterministic_tiebreak_to_lowest_index() {
        // identical profiles -> not decisive -> tied, and best resolves to the earliest copy.
        let copies = [cp(7, &[Some(A)], &[]), cp(3, &[Some(A)], &[])];
        let r = rf(&[Some(A)], &[]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 7); // first in slice order
        assert_eq!(a.status, AssignStatus::Tied);
    }

    #[test]
    fn empty_copies_returns_none() {
        let r = rf(&[Some(A)], &[]);
        assert!(assign_read(&r, &[], &AssignParams::default()).is_none());
    }

    #[test]
    fn read_psv_obs_reuses_cigar_bridge() {
        // 5M2N5M over ref_start 10: positions 10..15 and 17..22 matched; 15,16 inside the intron.
        let read = AlignedRead { ref_start: 10, cigar: vec![('M', 5), ('N', 2), ('M', 5)], seq: b"AAAAACCCCC".to_vec() };
        let obs = read_psv_obs(&read, &[12, 16, 18]);
        assert_eq!(obs, vec![Some(b'A'), None, Some(b'C')]);
    }
}
