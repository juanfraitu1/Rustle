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
//!
//! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use crate::types::DetHashSet;
use crate::vg_family::seq_utils::reverse_complement;

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
    use crate::vg_family::family_graph::{contiguous_core_coverage_bounded_with, EDGE_CONFIRM_ASTAR};
    let thin_up = crate::vg_family::family_graph::upper_cow(thin_seq);
    let mem_up = crate::vg_family::family_graph::upper_cow(&m.seq);
    let mut core_recip =
        contiguous_core_coverage_bounded_with(&thin_up, &mem_up, p.len_cap, EDGE_CONFIRM_ASTAR);
    let mut orientation = Orientation::Forward;
    if core_recip < p.t_core {
        let rc = contiguous_core_coverage_bounded_with(
            &thin_up,
            &reverse_complement(&mem_up),
            p.len_cap,
            EDGE_CONFIRM_ASTAR,
        );
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

// ---- merged 2026-10-05: was `vg_family/rescue_pipeline.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod rescue_pipeline {
//! Family-aware RESCUE thin-locus scan (integration stage 4b) — the BAM-neighbourhood orchestration of
//! `bench/family_rescue.py`.
//!
//! Recover an under-ASSEMBLED copy that falls below the general-assembly `>= 3`-read gate but is
//! sequence-homologous to a CONFIRMED de-novo family (borrow strength / partial pooling). The decision core
//! (`family_rescue::rescue_thin_locus`: canonical-k-mer pre-filter + both-orientation POA) is already
//! ported; this builds the thin candidate loci from reads in the family neighbourhood:
//!   1. group reads by exact intron chain, collapse OVERLAPPING chains into loci (best-supported chain wins);
//!   2. drop loci overlapping an already-assembled family member span;
//!   3. build each thin locus's spliced sequence (canonical-junction gate, reverse-complement for `-`);
//!   4. POA-confirm against the family members; dedup by locus, keep the best core coverage.
//!
//! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use std::collections::{BTreeMap, BTreeSet};

use crate::vg_family::denovo_assemble::{build_spliced_seq, PrimaryRead};
use super::{rescue_thin_locus, FamilyMember, RescueOutcome, RescueParams};
use crate::genome::GenomeIndex;

/// Minimum read support for a thin candidate locus (`family_rescue.py::MIN_SUPPORT`).
pub const RESCUE_MIN_SUPPORT: u32 = 1;
/// Minimum spliced length of a thin locus to attempt rescue.
pub const RESCUE_MIN_LEN: usize = 200;

/// A thin candidate locus: a (collapsed) intron chain with read support, below the general-assembly gate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ThinLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub support: u32,
    pub introns: Vec<(u64, u64)>,
}

/// An already-assembled family member's genomic span (used to EXCLUDE overlapping thin loci).
#[derive(Clone, Debug)]
pub struct MemberSpan {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

/// A rescued copy: the thin locus, the family/member that confirmed it (+ core coverage + orientation), the
/// transcription strand, and the spliced sequence (so the rescued copy is assignable).
#[derive(Clone, Debug)]
pub struct RescuedCopy {
    pub locus: ThinLocus,
    pub outcome: RescueOutcome,
    pub strand: char,
    pub seq: Vec<u8>,
}

/// Group reads by exact intron chain, then collapse OVERLAPPING chains (by span) into loci, keeping the
/// best-supported chain (by support, then span) as each locus's representative. Multi-exon only. Mirrors
/// `family_rescue.py`'s per-interval locus collapse. Deterministic (sorted by start, ties by chain).
pub fn thin_loci(reads: &[PrimaryRead], min_support: u32) -> Vec<ThinLocus> {
    // group reads by (chrom, exact intron chain) -> (support, min_start, max_end); multi-exon only.
    let mut chains: BTreeMap<(&str, Vec<(u64, u64)>), (u32, u64, u64)> = BTreeMap::new();
    for r in reads {
        if r.introns.is_empty() {
            continue;
        }
        let e = chains.entry((r.chrom.as_str(), r.introns.clone())).or_insert((0, u64::MAX, 0));
        e.0 += 1;
        e.1 = e.1.min(r.ref_start);
        e.2 = e.2.max(r.ref_end);
    }
    // per chrom: collapse overlapping chains (single-linkage by span); the rep is the best (support, span).
    let mut by_chrom: BTreeMap<&str, Vec<(Vec<(u64, u64)>, u32, u64, u64)>> = BTreeMap::new();
    for ((chrom, chain), (sup, s, e)) in chains {
        by_chrom.entry(chrom).or_default().push((chain, sup, s, e));
    }
    struct Loc {
        s: u64,
        e: u64,
        sup: u32,
        s2: u64,
        e2: u64,
        intr: Vec<(u64, u64)>,
    }
    let mut out = Vec::new();
    for (chrom, mut group) in by_chrom {
        group.sort_by_key(|&(_, _, s, _)| s); // stable: ties keep (chrom,chain) order
        let mut loci: Vec<Loc> = Vec::new();
        for (intr, sup, s, e) in group {
            if sup < min_support {
                continue;
            }
            let mut merged = false;
            for l in loci.iter_mut() {
                if s <= l.e && e >= l.s {
                    l.s = l.s.min(s);
                    l.e = l.e.max(e);
                    if (sup, e - s) > (l.sup, l.e2 - l.s2) {
                        l.sup = sup;
                        l.s2 = s;
                        l.e2 = e;
                        l.intr = intr.clone();
                    }
                    merged = true;
                    break;
                }
            }
            if !merged {
                loci.push(Loc { s, e, sup, s2: s, e2: e, intr });
            }
        }
        for l in loci {
            out.push(ThinLocus {
                chrom: chrom.to_string(),
                start: l.s2,
                end: l.e2,
                support: l.sup,
                introns: l.intr,
            });
        }
    }
    out
}

/// Rescue under-assembled copies: for each thin locus NOT overlapping a member span, build its spliced
/// sequence and POA-confirm it against the family `members` (`rescue_thin_locus`). Dedup by locus, keeping
/// the best `core_recip`. Mirrors `family_rescue.py`'s scan + dedup.
pub fn rescue_thin_loci(
    loci: &[ThinLocus],
    members: &[FamilyMember],
    member_spans: &[MemberSpan],
    genome: &GenomeIndex,
    p: &RescueParams,
) -> Vec<RescuedCopy> {
    let mut by_key: BTreeMap<(String, u64, u64), RescuedCopy> = BTreeMap::new();
    for locus in loci {
        // exclude a thin locus that overlaps an already-assembled family member span.
        if member_spans
            .iter()
            .any(|m| m.chrom == locus.chrom && locus.start < m.end && locus.end > m.start)
        {
            continue;
        }
        let (seq, strand) = match build_spliced_seq(genome, &locus.chrom, locus.start, locus.end, &locus.introns, None) {
            Some(v) => v,
            None => continue,
        };
        if seq.len() < RESCUE_MIN_LEN || seq.len() > p.len_cap {
            continue;
        }
        if let Some(outcome) = rescue_thin_locus(&seq, members, p) {
            let key = (locus.chrom.clone(), locus.start, locus.end);
            let better = by_key
                .get(&key)
                .map_or(true, |prev| outcome.core_recip > prev.outcome.core_recip);
            if better {
                by_key.insert(key, RescuedCopy { locus: locus.clone(), outcome, strand, seq });
            }
        }
    }
    by_key.into_values().collect()
}

/// ITERATIVE family-aware rescue (borrow strength across passes). A rescued copy is itself a new family
/// member that can bridge to OTHER under-assembled copies the first pass couldn't reach (homologous to the
/// rescued copy but not the original family). So: rescue, fold the rescued copies in as members + spans,
/// rescue the remaining loci again, until a pass recovers nothing new (or `max_iters`). Deterministic.
pub fn rescue_thin_loci_iterative(
    loci: &[ThinLocus],
    members: &[FamilyMember],
    member_spans: &[MemberSpan],
    genome: &GenomeIndex,
    p: &RescueParams,
    max_iters: usize,
) -> Vec<RescuedCopy> {
    let mut all_members: Vec<FamilyMember> = members.to_vec();
    let mut all_spans: Vec<MemberSpan> = member_spans.to_vec();
    let mut remaining: Vec<ThinLocus> = loci.to_vec();
    let mut rescued_all: Vec<RescuedCopy> = Vec::new();
    for _ in 0..max_iters {
        let rescued = rescue_thin_loci(&remaining, &all_members, &all_spans, genome, p);
        if rescued.is_empty() {
            break;
        }
        // each rescued copy becomes a member (so it can bridge) and a span (so it is not re-rescued).
        let done: BTreeSet<(String, u64, u64)> =
            rescued.iter().map(|rc| (rc.locus.chrom.clone(), rc.locus.start, rc.locus.end)).collect();
        for rc in &rescued {
            all_members.push(FamilyMember::new(
                format!("RC_{}_{}", rc.locus.chrom, rc.locus.start),
                rc.outcome.family_id.clone(),
                rc.seq.clone(),
            ));
            all_spans.push(MemberSpan {
                chrom: rc.locus.chrom.clone(),
                start: rc.locus.start,
                end: rc.locus.end,
            });
        }
        remaining.retain(|l| !done.contains(&(l.chrom.clone(), l.start, l.end)));
        rescued_all.extend(rescued);
        if remaining.is_empty() {
            break;
        }
    }
    rescued_all
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn read(chrom: &str, s: u64, e: u64, introns: &[(u64, u64)]) -> PrimaryRead {
        PrimaryRead { chrom: chrom.into(), ref_start: s, ref_end: e, introns: introns.to_vec(), reverse: false }
    }

    /// Genome with one thin gene: exon1 [0,200), canonical intron [200,220) GT..AG, exon2 [220,420);
    /// spliced = flank(50) + `core`(300) + flank(50). Returns the genome.
    fn thin_gene_genome(core: &[u8]) -> GenomeIndex {
        let mut g = vec![b'A'; 500];
        g[0..50].copy_from_slice(&rand_seq(50, 0x71));
        g[50..200].copy_from_slice(&core[0..150]);
        g[200] = b'G';
        g[201] = b'T';
        g[218] = b'A';
        g[219] = b'G';
        g[220..370].copy_from_slice(&core[150..300]);
        g[370..420].copy_from_slice(&rand_seq(50, 0x72));
        GenomeIndex::from_seqs(&[("c1", &g)])
    }

    // ---- thin_loci ----

    #[test]
    fn thin_loci_groups_and_keeps_multi_exon() {
        let reads = [
            read("c1", 0, 400, &[(100, 200)]),
            read("c1", 5, 410, &[(100, 200)]),
            read("c1", 0, 400, &[]), // single-exon -> ignored
        ];
        let loci = thin_loci(&reads, 1);
        assert_eq!(loci.len(), 1);
        assert_eq!(loci[0].support, 2);
        assert_eq!(loci[0].introns, vec![(100, 200)]);
    }

    #[test]
    fn thin_loci_collapses_overlapping_chains_best_supported() {
        // chain A [0,400) support 2, chain B [100,500) support 1 overlap -> one locus, rep = chain A.
        let reads = [
            read("c1", 0, 400, &[(100, 200)]),
            read("c1", 0, 400, &[(100, 200)]),
            read("c1", 100, 500, &[(150, 250)]),
        ];
        let loci = thin_loci(&reads, 1);
        assert_eq!(loci.len(), 1, "overlapping chains collapse to one locus");
        assert_eq!(loci[0].support, 2, "best-supported chain wins");
        assert_eq!(loci[0].introns, vec![(100, 200)]);
        assert_eq!((loci[0].start, loci[0].end), (0, 400));
    }

    #[test]
    fn thin_loci_separates_nonoverlapping() {
        let reads = [
            read("c1", 0, 300, &[(100, 200)]),
            read("c1", 1000, 1300, &[(1100, 1200)]),
        ];
        assert_eq!(thin_loci(&reads, 1).len(), 2);
    }

    // ---- rescue_thin_loci ----

    #[test]
    fn rescue_recovers_thin_locus_homologous_to_family() {
        let core = rand_seq(300, 0xC0FE_F00D);
        let genome = thin_gene_genome(&core);
        // a family member sharing the same core (distinct flanks)
        let mseq = cat(&[&rand_seq(50, 0x81), &core, &rand_seq(50, 0x82)]);
        let members = [FamilyMember::new("M1".into(), "FAM1".into(), mseq)];
        // a SINGLE read (support 1, below the >=3 gate) forming the thin locus
        let reads = [read("c1", 0, 420, &[(200, 220)])];
        let loci = thin_loci(&reads, RESCUE_MIN_SUPPORT);
        assert_eq!(loci.len(), 1);
        let rescued = rescue_thin_loci(&loci, &members, &[], &genome, &RescueParams::default());
        assert_eq!(rescued.len(), 1, "the thin copy is rescued into the family");
        assert_eq!(rescued[0].outcome.family_id, "FAM1");
        assert!(rescued[0].outcome.core_recip >= 0.13);
        assert_eq!(rescued[0].seq.len(), 400);
    }

    #[test]
    fn rescue_excludes_loci_overlapping_a_member_span() {
        let core = rand_seq(300, 0xC0FE_F00D);
        let genome = thin_gene_genome(&core);
        let mseq = cat(&[&rand_seq(50, 0x81), &core, &rand_seq(50, 0x82)]);
        let members = [FamilyMember::new("M1".into(), "FAM1".into(), mseq)];
        let reads = [read("c1", 0, 420, &[(200, 220)])];
        let loci = thin_loci(&reads, RESCUE_MIN_SUPPORT);
        // a member already assembled across the locus span -> the thin locus is NOT a new copy.
        let spans = [MemberSpan { chrom: "c1".into(), start: 100, end: 300 }];
        let rescued = rescue_thin_loci(&loci, &members, &spans, &genome, &RescueParams::default());
        assert!(rescued.is_empty(), "locus overlapping an assembled member is excluded");
    }

    #[test]
    fn rescue_rejects_non_homologous_thin_locus() {
        let genome = thin_gene_genome(&rand_seq(300, 0xC0FE_F00D));
        // a family member with a DIFFERENT, unrelated core -> no rescue.
        let mseq = cat(&[&rand_seq(50, 0x81), &rand_seq(300, 0xDEAD_BEEF), &rand_seq(50, 0x82)]);
        let members = [FamilyMember::new("M1".into(), "FAM1".into(), mseq)];
        let reads = [read("c1", 0, 420, &[(200, 220)])];
        let loci = thin_loci(&reads, RESCUE_MIN_SUPPORT);
        let rescued = rescue_thin_loci(&loci, &members, &[], &genome, &RescueParams::default());
        assert!(rescued.is_empty(), "a non-homologous thin locus is not rescued");
    }

    #[test]
    fn iterative_rescue_recovers_a_bridged_copy() {
        // M1 = flank + core1. L1 spliced = core1 + core2 (rescued pass 1 via core1, shared with M1).
        // L2 spliced = core2 + flank (homologous to L1's core2 but NOT to M1) -> rescued ONLY in pass 2,
        // once L1 is a member. Single-pass recovers 1; iterative recovers 2.
        let core1 = rand_seq(200, 0xC0DE_0001);
        let core2 = rand_seq(200, 0xC0DE_0002);
        let mut g = vec![b'A'; 1600];
        // L1 = core1 | core2 ; L2 = flankL | core2 ; M1 = core1 | flankM. Shared cores sit at the SAME
        // relative position (no offset) so POA anchors them cleanly.
        g[0..200].copy_from_slice(&core1);
        g[200] = b'G';
        g[201] = b'T';
        g[218] = b'A';
        g[219] = b'G';
        g[220..420].copy_from_slice(&core2);
        g[1000..1200].copy_from_slice(&rand_seq(200, 0xF00D)); // flankL
        g[1200] = b'G';
        g[1201] = b'T';
        g[1218] = b'A';
        g[1219] = b'G';
        g[1220..1420].copy_from_slice(&core2);
        let genome = GenomeIndex::from_seqs(&[("c1", &g)]);
        let m1 = FamilyMember::new("M1".into(), "FAM1".into(), cat(&[&core1, &rand_seq(200, 0xAA)]));
        let loci = [
            ThinLocus { chrom: "c1".into(), start: 0, end: 420, support: 1, introns: vec![(200, 220)] },
            ThinLocus { chrom: "c1".into(), start: 1000, end: 1420, support: 1, introns: vec![(1200, 1220)] },
        ];
        let single = rescue_thin_loci(&loci, std::slice::from_ref(&m1), &[], &genome, &RescueParams::default());
        assert_eq!(single.len(), 1, "single-pass rescues only the directly-homologous locus");
        let iter = rescue_thin_loci_iterative(&loci, &[m1], &[], &genome, &RescueParams::default(), 5);
        assert_eq!(iter.len(), 2, "iterative recovers the bridged locus via the first rescued copy");
    }
}
}
