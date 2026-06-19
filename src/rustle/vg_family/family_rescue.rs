//! Family-aware copy RESCUE (borrow strength / partial pooling): recover an
//! under-ASSEMBLED copy that falls below the general-assembly read gate but is
//! sequence-homologous to a CONFIRMED de-novo family.
//!
//! Rationale (mirrors `bench/family_rescue.py`): a single confident
//! canonical-junction read forming a multi-exon chain that POA-confirms against
//! an existing family is strong evidence of a real copy -- the family PRIOR
//! overcomes thin read support. This recovers expression-limited copies (e.g.
//! RFPL2, a single read) WITHOUT lowering the genome-wide `>= 3`-read gate
//! (fuzzy / low-gate experiments only added noise; see the pipeline record).
//!
//! This module is the testable DECISION CORE. The BAM-neighbourhood scan that
//! builds the thin candidate loci (intron-chain collapse, member-span exclusion,
//! canonical-junction strand check, spliced-sequence construction) is the
//! integration layer (`gen2off`-style boundary mapping + BAM/FASTA), kept out of
//! here so the rescue logic is unit-testable against `contiguous_core_coverage`.
//!
//! The decision, per candidate thin locus:
//!   1. A CANONICAL exact base-4 k-mer pre-filter (`canonical_kmer_set`,
//!      `KMER = 18`, strand-symmetric) selects the SINGLE best family member by
//!      shared-k-mer overlap (`>= K_RESCUE`). This bounds POA to true candidates
//!      and never decides membership -- POA does.
//!   2. POA-confirm the thin locus against that one member via the validated
//!      `contiguous_core_coverage` primitive, trying BOTH orientations (forward,
//!      then a reverse-complement fallback for copies assembled on the opposite
//!      strand). Rescue iff the contiguous-core coverage `>= T_CORE = 0.13`.

use super::family_graph::contiguous_core_coverage;
use crate::types::DetHashSet;
use crate::vg::reverse_complement;

/// Exact k-mer length for the canonical pre-filter. `KMER = 18` gives a
/// `4^18 ~ 6.9e10` space (negligible coincidental sharing) and fits in `u64`
/// (36 bits). Matches `bench/denovo_families.py::KMER`.
pub const KMER: usize = 18;

/// POA contiguous-core coverage threshold to confirm a rescued copy. A true
/// recent-duplicate copy shares one long homologous core (`>= 0.13` of the
/// shorter sequence); a domain-sharer co-aligns over only a short block.
/// Matches `bench/family_rescue.py::T_CORE`.
pub const T_CORE: f64 = 0.13;

/// Minimum number of shared canonical k-mers a thin locus must have with a
/// family member to even reach POA (real copies share many; common-domain
/// sharers share few). Matches `bench/family_rescue.py::K_RESCUE`.
pub const K_RESCUE: usize = 20;

/// Length cap (bp) on the POA inputs (cost guard; POA is O(L^2)). Matches
/// `bench/family_rescue.py::LEN_CAP`.
pub const LEN_CAP: usize = 9000;

/// 2-bit base code A=0, C=1, G=2, T=3 (lowercase accepted), matching the python
/// `_CODE` table. Returns `None` for any non-ACGT base (so a window touching it
/// is dropped, as in the python `badwin` mask).
fn base_code(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// CANONICAL base-4 code of a single window = `min(forward, reverse-complement)`
/// so a window and its reverse-complement collapse to the SAME code (strand
/// symmetry). Returns `None` if the window touches any non-ACGT base.
///
/// Faithful to `bench/denovo_families.py::kmer_hashes`: forward Horner code
/// `Σ base_t · 4^(k-1-t)`; reverse-complement code substitutes `3 - base` and
/// reverses the position order.
///
/// `pub(crate)` so the family-detection layer can reuse the exact same canonical
/// encoder for its position-aware k-mer signatures.
pub(crate) fn window_canon_code(window: &[u8]) -> Option<u64> {
    let mut fwd: u64 = 0;
    let mut rc: u64 = 0;
    for (i, &b) in window.iter().enumerate() {
        let c = base_code(b)?; // None drops a window touching a non-ACGT base
        fwd = (fwd << 2) | c as u64;
        // reverse-complement code: complement (3 - c), placed in reverse order.
        rc |= ((3 - c) as u64) << (2 * i);
    }
    Some(fwd.min(rc))
}

/// CANONICAL exact base-4 k-mer set (`KMER`-mers, `min(fwd, rc)` per window;
/// windows touching a non-ACGT base dropped). A sequence and its reverse
/// complement yield the SAME set, so rescue is strand-symmetric. Mirrors
/// `set(kmer_hashes(seq))` in the python pipeline.
pub fn canonical_kmer_set(seq: &[u8]) -> DetHashSet<u64> {
    let mut set = DetHashSet::default();
    if seq.len() < KMER {
        return set;
    }
    for w in seq.windows(KMER) {
        if let Some(code) = window_canon_code(w) {
            set.insert(code);
        }
    }
    set
}

/// Number of shared canonical k-mers between two sets (the pre-filter overlap).
pub fn kmer_overlap(a: &DetHashSet<u64>, b: &DetHashSet<u64>) -> usize {
    // Iterate the smaller set for speed; the count is order-independent.
    let (small, large) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    small.iter().filter(|k| large.contains(k)).count()
}

/// A confirmed de-novo family member the rescue compares thin loci against. Its
/// canonical k-mer set is precomputed once (the python pipeline caches
/// `memkmers[t]`).
#[derive(Clone, Debug)]
pub struct FamilyMember {
    /// Member transcript id (e.g. a de-novo transcript id).
    pub tid: String,
    /// The family this member belongs to.
    pub family_id: String,
    /// The member's spliced sequence in transcription orientation.
    pub seq: Vec<u8>,
    /// Precomputed canonical k-mer set of `seq`.
    pub kmers: DetHashSet<u64>,
}

impl FamilyMember {
    /// Build a member, computing its canonical k-mer set from `seq`.
    pub fn new(tid: String, family_id: String, seq: Vec<u8>) -> Self {
        let kmers = canonical_kmer_set(&seq);
        FamilyMember { tid, family_id, seq, kmers }
    }
}

/// Which orientation of the member the contiguous core was found in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orientation {
    /// The member's stored orientation aligned the core directly.
    Forward,
    /// The reverse-complement of the member aligned the core (opposite-strand
    /// assembly).
    ReverseComplement,
}

/// Tunable rescue parameters (defaults mirror `bench/family_rescue.py`).
#[derive(Clone, Copy, Debug)]
pub struct RescueParams {
    /// POA contiguous-core coverage required to confirm a copy.
    pub t_core: f64,
    /// Minimum shared canonical k-mers with the best member to reach POA.
    pub k_rescue: usize,
    /// Length cap on POA inputs.
    pub len_cap: usize,
}

impl Default for RescueParams {
    fn default() -> Self {
        RescueParams { t_core: T_CORE, k_rescue: K_RESCUE, len_cap: LEN_CAP }
    }
}

/// A rescued copy: which family it joins, the member that confirmed it, the POA
/// contiguous-core coverage, and the orientation the core aligned in.
#[derive(Clone, Debug)]
pub struct RescueOutcome {
    pub family_id: String,
    pub best_member: String,
    pub core_recip: f64,
    pub orientation: Orientation,
}

/// Decide whether a thin-locus spliced sequence is a rescued copy of one of the
/// confirmed families represented by `members`.
///
/// `members` is the candidate neighbourhood the integration layer already
/// windowed (the python `WIN = 1 Mb` member set near the locus). The function
/// pre-filters to the SINGLE best member by canonical-k-mer overlap
/// (`>= p.k_rescue`), then POA-confirms via `contiguous_core_coverage` in the
/// forward orientation, falling back to the reverse complement if forward is
/// below `p.t_core`. Returns the rescue iff the best coverage reaches
/// `p.t_core`. Deterministic: members are scanned in slice order and the best
/// overlap ties to the earliest member (matching the python `if ov > best_ov`).
pub fn rescue_thin_locus(
    thin_seq: &[u8],
    members: &[FamilyMember],
    p: &RescueParams,
) -> Option<RescueOutcome> {
    let kset = canonical_kmer_set(thin_seq);
    if kset.is_empty() {
        return None;
    }

    // Pre-filter: the SINGLE best member by canonical-k-mer overlap. `best_ov`
    // starts at `k_rescue - 1` so a member must clear `>= k_rescue` to be picked
    // (python `best_ov = K_RESCUE - 1; if ov > best_ov`). Strict `>` ties to the
    // earliest member in slice order (determinism).
    let mut best: Option<&FamilyMember> = None;
    let mut best_ov = p.k_rescue.saturating_sub(1);
    for m in members {
        // POA cost guard (python `_poa_rescue`: skip if min(len) > LEN_CAP).
        if thin_seq.len().min(m.seq.len()) > p.len_cap {
            continue;
        }
        let ov = kmer_overlap(&kset, &m.kmers);
        if ov > best_ov {
            best_ov = ov;
            best = Some(m);
        }
    }
    let m = best?;

    // POA-confirm against that one member. The k-mer pre-filter is
    // case-insensitive (like python's `_CODE` table), so uppercase the POA
    // operands here to match it: python's `poa_pair_stats` uppercases its inputs,
    // and `reverse_complement` maps lowercase -> N, which would otherwise silently
    // void the RC fallback on soft-masked input. Forward first; reverse-complement
    // fallback only if forward is below threshold, keeping the max (mirrors the
    // python `_poa_rescue` RC retry).
    let thin_up = thin_seq.to_ascii_uppercase();
    let mem_up = m.seq.to_ascii_uppercase();
    let mut core_recip = contiguous_core_coverage(&thin_up, &mem_up);
    let mut orientation = Orientation::Forward;
    if core_recip < p.t_core {
        let rc = contiguous_core_coverage(&thin_up, &reverse_complement(&mem_up));
        if rc > core_recip {
            core_recip = rc;
            orientation = Orientation::ReverseComplement;
        }
    }

    if core_recip >= p.t_core {
        Some(RescueOutcome {
            family_id: m.family_id.clone(),
            best_member: m.tid.clone(),
            core_recip,
            orientation,
        })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Deterministic DNA (SplitMix64), mirroring family_graph's core-coverage tests.
    struct SplitMix64(u64);
    impl SplitMix64 {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
    }
    fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut rng = SplitMix64(seed);
        const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
        (0..n).map(|_| B[(rng.next_u64() % 4) as usize]).collect()
    }
    fn cat(parts: &[&[u8]]) -> Vec<u8> {
        parts.iter().flat_map(|p| p.iter().copied()).collect()
    }
    fn member(tid: &str, fid: &str, seq: Vec<u8>) -> FamilyMember {
        FamilyMember::new(tid.to_string(), fid.to_string(), seq)
    }

    // ---- canonical k-mer encoding (faithful to python base-4 codes) ----

    #[test]
    fn window_canon_code_is_base4_and_strand_canonical() {
        // A=0,C=1,G=2,T=3 ; fwd("AC")=0*4+1=1 ; RC("AC")="GT"=2*4+3=11 ; min=1.
        assert_eq!(window_canon_code(b"AC"), Some(1));
        // "GT" is the reverse-complement of "AC" -> SAME canonical code.
        assert_eq!(window_canon_code(b"GT"), Some(1));
        // a window containing a non-ACGT base is dropped.
        assert_eq!(window_canon_code(b"AN"), None);
    }

    #[test]
    fn canonical_kmer_set_is_strand_symmetric() {
        let s = rand_seq(200, 0x5EED_0001);
        let rc = reverse_complement(&s);
        assert_eq!(canonical_kmer_set(&s), canonical_kmer_set(&rc));
        assert!(!canonical_kmer_set(&s).is_empty());
    }

    #[test]
    fn canonical_kmer_set_drops_only_windows_touching_n() {
        // 60 bp random -> 60-18+1 = 43 windows; random 18-mers are distinct.
        let s = rand_seq(60, 0x5EED_0002);
        let full = canonical_kmer_set(&s);
        assert_eq!(full.len(), 60 - KMER + 1, "distinct random 18-mers");
        let mut withn = s.clone();
        withn[30] = b'N';
        let dropped = canonical_kmer_set(&withn);
        // windows with start in [13, 30] cover index 30 -> 18 windows removed.
        assert_eq!(full.len() - dropped.len(), 18);
    }

    #[test]
    fn kmer_overlap_counts_shared_canonical_kmers() {
        let s = rand_seq(60, 0x7001);
        let a = canonical_kmer_set(&s);
        let b = canonical_kmer_set(&reverse_complement(&s));
        assert_eq!(kmer_overlap(&a, &b), a.len()); // identical canonical sets
        let c = canonical_kmer_set(&rand_seq(60, 0x7002));
        assert!(kmer_overlap(&a, &c) < a.len()); // independent seqs share ~none
    }

    // ---- rescue decision ----

    #[test]
    fn rescue_confirms_homologous_copy() {
        // thin locus and a member share a 400 bp identical core (>> K_RESCUE
        // shared 18-mers, contiguous-core coverage >= 0.13), divergent flanks.
        let core = rand_seq(400, 0xC0FE_2001);
        let thin = cat(&[&rand_seq(80, 0xAAAA_2001), &core, &rand_seq(80, 0xAAAA_2002)]);
        let mseq = cat(&[&rand_seq(80, 0xBBBB_2001), &core, &rand_seq(80, 0xBBBB_2002)]);
        let members = [member("M1", "FAM7", mseq)];
        let out = rescue_thin_locus(&thin, &members, &RescueParams::default())
            .expect("homologous copy should be rescued");
        assert_eq!(out.family_id, "FAM7");
        assert_eq!(out.best_member, "M1");
        assert_eq!(out.orientation, Orientation::Forward);
        assert!(out.core_recip >= T_CORE, "core_recip {} >= {}", out.core_recip, T_CORE);
    }

    #[test]
    fn rescue_rejects_domain_sharer_below_core_threshold() {
        // shares a 40 bp block (40-18+1 = 23 shared 18-mers >= K_RESCUE so the
        // PRE-FILTER passes) but in long otherwise-independent sequences, so the
        // POA contiguous-core coverage is < T_CORE -> NOT rescued (POA decides).
        let block = rand_seq(40, 0xD0D0_3001);
        let thin = cat(&[&rand_seq(350, 0xAAAA_3001), &block, &rand_seq(350, 0xAAAA_3002)]);
        let mseq = cat(&[&rand_seq(350, 0xBBBB_3001), &block, &rand_seq(350, 0xBBBB_3002)]);
        let members = [member("M1", "FAM3", mseq)];
        // precondition: the pre-filter DOES pass (>= K_RESCUE shared k-mers).
        let thin_k = canonical_kmer_set(&thin);
        assert!(kmer_overlap(&thin_k, &members[0].kmers) >= K_RESCUE,
            "pre-filter precondition: shared k-mers >= K_RESCUE");
        assert!(rescue_thin_locus(&thin, &members, &RescueParams::default()).is_none());
    }

    #[test]
    fn rescue_pre_filter_rejects_too_few_shared_kmers() {
        // shares only a 25 bp block -> 25-18+1 = 8 shared 18-mers < K_RESCUE
        // -> NO member clears the pre-filter -> no POA, no rescue.
        let block = rand_seq(25, 0xD0D0_4001);
        let thin = cat(&[&rand_seq(300, 0xAAAA_4001), &block, &rand_seq(300, 0xAAAA_4002)]);
        let mseq = cat(&[&rand_seq(300, 0xBBBB_4001), &block, &rand_seq(300, 0xBBBB_4002)]);
        let members = [member("M1", "FAM4", mseq)];
        assert!(rescue_thin_locus(&thin, &members, &RescueParams::default()).is_none());
    }

    #[test]
    fn rescue_detects_reverse_complement_copy() {
        // A member stored on the OPPOSITE strand: canonical k-mers still match
        // (strand-symmetric), forward POA is LOW, RC fallback POA is HIGH.
        let core = rand_seq(400, 0xC0FE_5001);
        let thin = cat(&[&rand_seq(80, 0xAAAA_5001), &core, &rand_seq(80, 0xAAAA_5002)]);
        let mfwd = cat(&[&rand_seq(80, 0xBBBB_5001), &core, &rand_seq(80, 0xBBBB_5002)]);
        let mseq = reverse_complement(&mfwd); // member assembled on the other strand
        let members = [member("M1", "FAM5", mseq)];
        let out = rescue_thin_locus(&thin, &members, &RescueParams::default())
            .expect("RC homologous copy should be rescued via the RC fallback");
        assert_eq!(out.orientation, Orientation::ReverseComplement);
        assert!(out.core_recip >= T_CORE, "core_recip {} >= {}", out.core_recip, T_CORE);
    }

    #[test]
    fn rescue_picks_best_kmer_matching_member() {
        // thin shares its full 400 bp core with M2 (many k-mers) and nothing with
        // M1 (unrelated); POA must confirm against M2 and report it.
        let core = rand_seq(400, 0xC0FE_6001);
        let thin = cat(&[&rand_seq(80, 0xAAAA_6001), &core, &rand_seq(80, 0xAAAA_6002)]);
        let m1 = rand_seq(660, 0x1111_6001); // fully unrelated
        let m2 = cat(&[&rand_seq(80, 0xBBBB_6001), &core, &rand_seq(80, 0xBBBB_6002)]);
        let members = [member("M1", "FAMx", m1), member("M2", "FAM6", m2)];
        let out = rescue_thin_locus(&thin, &members, &RescueParams::default())
            .expect("should rescue via M2");
        assert_eq!(out.best_member, "M2");
        assert_eq!(out.family_id, "FAM6");
    }

    #[test]
    fn rescue_empty_members_returns_none() {
        let thin = rand_seq(300, 0x9001);
        assert!(rescue_thin_locus(&thin, &[], &RescueParams::default()).is_none());
    }

    #[test]
    fn rescue_thin_seq_shorter_than_kmer_returns_none() {
        let thin = rand_seq(10, 0x9002); // < KMER -> empty k-mer set -> None
        let m = member("M1", "F", rand_seq(400, 0x9003));
        assert!(rescue_thin_locus(&thin, &[m], &RescueParams::default()).is_none());
    }

    #[test]
    fn rescue_detects_reverse_complement_copy_lowercase() {
        // Soft-masked (lowercase) inputs must behave IDENTICALLY to uppercase. The
        // k-mer pre-filter is case-insensitive, so a lowercase opposite-strand copy
        // reaches POA; the RC fallback must not be silently killed by a
        // case-sensitive reverse_complement (which maps lowercase -> N). Same
        // geometry as rescue_detects_reverse_complement_copy, just lowercased.
        let core = rand_seq(400, 0xC0FE_5001);
        let thin = cat(&[&rand_seq(80, 0xAAAA_5001), &core, &rand_seq(80, 0xAAAA_5002)]);
        let mfwd = cat(&[&rand_seq(80, 0xBBBB_5001), &core, &rand_seq(80, 0xBBBB_5002)]);
        let thin_lc = thin.to_ascii_lowercase();
        let mseq_lc = reverse_complement(&mfwd).to_ascii_lowercase();
        let members = [member("M1", "FAM5", mseq_lc)];
        let out = rescue_thin_locus(&thin_lc, &members, &RescueParams::default())
            .expect("lowercase RC homologous copy should still be rescued via the RC fallback");
        assert_eq!(out.orientation, Orientation::ReverseComplement);
        assert!(out.core_recip >= T_CORE, "core_recip {} >= {}", out.core_recip, T_CORE);
    }

    #[test]
    fn rescue_tie_break_picks_earliest_member() {
        // Two members with the IDENTICAL sequence (so provably EQUAL k-mer overlap
        // with the thin locus); the strict `>` update must keep the FIRST in slice
        // order. (Distinct random flanks do NOT guarantee equal overlap — a flank
        // k-mer can coincidentally collide — so identical seqs pin the tie exactly.)
        let core = rand_seq(400, 0xC0FE_7001);
        let thin = cat(&[&rand_seq(80, 0xAAAA_7001), &core, &rand_seq(80, 0xAAAA_7002)]);
        let mseq = cat(&[&rand_seq(80, 0x1111_7001), &core, &rand_seq(80, 0x1111_7002)]);
        let members = [member("FIRST", "FAM7", mseq.clone()), member("SECOND", "FAM7", mseq)];
        // precondition: the two overlaps are genuinely equal (the tie under test).
        let tk = canonical_kmer_set(&thin);
        assert_eq!(
            kmer_overlap(&tk, &members[0].kmers),
            kmer_overlap(&tk, &members[1].kmers),
            "members must have equal overlap for this to test the tie-break"
        );
        let out = rescue_thin_locus(&thin, &members, &RescueParams::default())
            .expect("should rescue");
        assert_eq!(out.best_member, "FIRST", "equal overlap must tie to the earliest member");
    }

    #[test]
    fn rescue_threshold_boundary_is_inclusive() {
        // Pin the `core_recip >= t_core` semantics: at t_core == cr the copy is
        // ACCEPTED (inclusive), just above cr it is REJECTED. `cr` is read off the
        // SAME deterministic primitive the function uses (forward orientation).
        use crate::vg_family::family_graph::contiguous_core_coverage;
        let core = rand_seq(400, 0xC0FE_8001);
        let thin = cat(&[&rand_seq(80, 0xAAAA_8001), &core, &rand_seq(80, 0xAAAA_8002)]);
        let mseq = cat(&[&rand_seq(80, 0xBBBB_8001), &core, &rand_seq(80, 0xBBBB_8002)]);
        let cr = contiguous_core_coverage(&thin, &mseq);
        let members = [member("M1", "FAM8", mseq)];
        let at = RescueParams { t_core: cr, ..RescueParams::default() };
        assert!(rescue_thin_locus(&thin, &members, &at).is_some(),
            "t_core == cr ({cr}) must be ACCEPTED (inclusive >=)");
        let above = RescueParams { t_core: cr + 1e-6, ..RescueParams::default() };
        assert!(rescue_thin_locus(&thin, &members, &above).is_none(),
            "t_core just above cr ({cr}) must be REJECTED");
    }
}

