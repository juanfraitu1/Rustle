//! VG re-align supplement -- re-align poor-fit/unmapped reads to O1's family copy-paths,
//! significance-gated (correct + discover). Task 1: candidate selection. Task 3: re-align a
//! candidate read to the family's copy-paths (identity-based, DRY on `bridge_detector::aln_id`)
//! and route unmapped reads to candidate families by shared minimizers. Task 4: gate the
//! re-align correction behind a min_p significance certificate (same `epsilon^delta` form as
//! `copy_assign::read_copy_evidence`), and greedily pool reads that fit no existing copy into
//! candidate novel-copy clusters. Task 5: wire Tasks 1/3/4 into a per-family driver
//! (`run_family_realign`) and, behind `DenovoConfig::vg_realign`, into the pipeline -- REPORT-ONLY
//! (see that fn's doc for the exact scope).
//!
//! VG re-align END-TO-END plan, Task 1: `align_traceback` + `path_obs_at`. There is no `edlib`
//! crate; `bridge_detector::hw_distance` is a hand-rolled 2-row DP that gives the HW/infix edit
//! DISTANCE only, no alignment path. To re-extract a read's base at a copy-path's PSV columns
//! (follow-up (c) in `bench/VG_REALIGN.md`) we need the actual traceback, so this keeps a full DP
//! + backtrack matrix (not the rolling 2-row form) and reconstructs the aligned columns.

use std::collections::HashSet;

use crate::vg_family::bridge_detector::{aln_id, hw_distance, revcomp};
use crate::vg_family::copy_assign_pipeline::best_overlap_copy;
use crate::vg_family::denovo_assemble::BamRead;
use crate::vg_family::family_detect::DenovoTranscript;
use crate::vg_family::minimizers::minimizers;

/// Backtrack pointer for one DP cell of `align_traceback`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Trace {
    /// Row 0 (free leading gap on `target`) -- backtrack terminates here.
    Start,
    /// From `dp[i-1][j-1]`: consumes `query[i-1]` AND `target[j-1]` (match or mismatch).
    Diag,
    /// From `dp[i-1][j]`: consumes `query[i-1]` only -- a gap in `target`.
    Up,
    /// From `dp[i][j-1]`: consumes `target[j-1]` only -- a gap in `query`.
    Left,
}

/// HW/infix alignment of `query` against `target`, WITH the traceback path (unlike
/// `bridge_detector::hw_distance`, which only returns the distance via a rolling 2-row DP).
///
/// Same semantics as `hw_distance`: `query` is rows, `target` is columns; row 0 is `0` across
/// every target column (free leading gap on `target`); the alignment ends at the min-cost cell in
/// the LAST query row (free trailing gap on `target`) and is backtracked to query row 0. Match
/// cost 0, substitution/indel cost 1.
///
/// Returns the aligned columns in order from the start of the alignment to its end: `(Some(qi),
/// Some(ti))` for a match/mismatch, `(Some(qi), None)` for a gap in `target`, `(None, Some(ti))`
/// for a gap in `query`. Every `query` index `0..query.len()` appears exactly once (query has no
/// free end-gaps under HW); `target` indices outside the aligned span (the free leading/trailing
/// gap) never appear at all.
pub fn align_traceback(query: &[u8], target: &[u8]) -> Vec<(Option<usize>, Option<usize>)> {
    let lq = query.len();
    let lt = target.len();
    let cols = lt + 1;

    let mut dp: Vec<usize> = vec![0; (lq + 1) * cols];
    let mut back: Vec<Trace> = vec![Trace::Start; (lq + 1) * cols];

    // Row 0: free leading gap on target -> cost 0 everywhere; Start marks the backtrack terminus.
    // (dp[0][j] is already 0 from the `vec!` initializer; `back` is already `Trace::Start`.)

    for i in 1..=lq {
        dp[i * cols] = i;
        back[i * cols] = Trace::Up; // dp[i][0] = dp[i-1][0] + 1 (gap in target, column stays 0).
        let qi = query[i - 1];
        for j in 1..=lt {
            let sub = dp[(i - 1) * cols + (j - 1)] + if qi == target[j - 1] { 0 } else { 1 };
            let up = dp[(i - 1) * cols + j] + 1;
            let left = dp[i * cols + (j - 1)] + 1;

            let (best, tr) = if sub <= up && sub <= left {
                (sub, Trace::Diag)
            } else if up <= left {
                (up, Trace::Up)
            } else {
                (left, Trace::Left)
            };
            dp[i * cols + j] = best;
            back[i * cols + j] = tr;
        }
    }

    // Free trailing gap on target: the alignment ends at the min-cost cell of the last query row.
    let last_row = lq * cols;
    let mut best_j = 0usize;
    let mut best_cost = dp[last_row];
    for j in 1..=lt {
        let c = dp[last_row + j];
        if c < best_cost {
            best_cost = c;
            best_j = j;
        }
    }

    let mut pairs: Vec<(Option<usize>, Option<usize>)> = Vec::new();
    let mut i = lq;
    let mut j = best_j;
    while i > 0 {
        match back[i * cols + j] {
            Trace::Diag => {
                pairs.push((Some(i - 1), Some(j - 1)));
                i -= 1;
                j -= 1;
            }
            Trace::Up => {
                pairs.push((Some(i - 1), None));
                i -= 1;
            }
            Trace::Left => {
                pairs.push((None, Some(j - 1)));
                j -= 1;
            }
            Trace::Start => break, // unreachable for i > 0, but avoid looping forever if it were.
        }
    }
    pairs.reverse();
    pairs
}

/// For each PSV column `t` (a position in `target`, the copy-path/consensus that `align_map` was
/// computed against), the read's base observed there: `Some(query[qi])` when `align_map` has an
/// aligned pair `(Some(qi), Some(t))`, or `None` when the read gaps at that column or the
/// alignment doesn't span it at all (`t` falls in the free leading/trailing target gap, or lands
/// on a `(None, Some(t))` query-gap column). Output length matches
/// `psv_positions_in_consensus.len()`, in the same order.
pub fn path_obs_at(
    align_map: &[(Option<usize>, Option<usize>)],
    psv_positions_in_consensus: &[usize],
    query: &[u8],
) -> Vec<Option<u8>> {
    psv_positions_in_consensus
        .iter()
        .map(|&t| {
            align_map
                .iter()
                .find(|&&(_, ti)| ti == Some(t))
                .and_then(|&(qi, _)| qi.map(|q| query[q]))
        })
        .collect()
}

/// C1 fix: orient `read_seq` to match `copy_seq`'s coordinate frame before `align_traceback`.
///
/// Invariant (see `copy_assign_pipeline::fill_psv_obs`/`FamilyProfiles::strand`): a copy's
/// `DenovoTranscript::seq` (and hence `copy_seqs[k]` here) is in TRANSCRIPTION-strand orientation,
/// while a `BamRead`'s `read.seq` is FORWARD-GENOME orientation (as SAM/BAM already store it). For
/// a `+`-strand copy the two frames coincide; for a `-`-strand copy the copy's consensus is the
/// REVERSE COMPLEMENT of the forward-genome sequence at that locus, so aligning `read.seq` against
/// it literally (as `align_traceback` does, byte-for-byte, no orientation search of its own --
/// unlike `aln_id`, which tries both orientations for its identity SCORE) produces a ~0.5-identity
/// garbage alignment: the wrong or `None` `path_obs` bases that then feed the EM as bogus PSV
/// evidence.
///
/// Fix: try both `read_seq` and `revcomp(read_seq)` against `copy_seq` via `hw_distance` (the same
/// distance `aln_id` uses internally) and return whichever orientation fits better (ties keep the
/// forward/as-given orientation). Callers must align (and extract `path_obs_at` from) THIS
/// returned, oriented sequence -- not the raw `read_seq` -- so the observed bases end up in the
/// copy's own transcription-strand frame, directly comparable to `copy_psv_alleles`.
pub fn orient_for_copy(read_seq: &[u8], copy_seq: &[u8]) -> Vec<u8> {
    let d_fwd = hw_distance(read_seq, copy_seq);
    let d_rev = hw_distance(&revcomp(read_seq), copy_seq);
    if d_rev < d_fwd {
        revcomp(read_seq)
    } else {
        read_seq.to_vec()
    }
}

/// Thresholds for flagging a read as a re-align CANDIDATE (poor-fit or unmapped).
///
/// A read is a candidate if it is low-MAPQ (ambiguous/multi-mapped), heavily soft/hard-clipped
/// (partial alignment, suggesting the reference copy it landed on isn't its true source), or
/// highly divergent (many mismatches relative to the copy it aligned to). Divergence and clip
/// fraction are computed by the CALLER from the CIGAR/NM at wiring time (later task); this
/// struct only carries the thresholds `is_candidate` applies.
pub struct RealignParams {
    pub max_mapq: u8,
    pub min_clip_frac: f64,
    pub min_div: f64,
    pub min_reads: usize,
}

impl Default for RealignParams {
    fn default() -> Self {
        RealignParams { max_mapq: 20, min_clip_frac: 0.20, min_div: 0.05, min_reads: 3 }
    }
}

/// True iff `(mapq, div, clip_frac)` indicates a poor primary-alignment fit under `p`'s
/// thresholds: low MAPQ (`<= max_mapq`), OR heavy clipping (`>= min_clip_frac`), OR high
/// divergence (`>= min_div`). Pure function of the three scalars + params -- `div` and
/// `clip_frac` are computed by the caller from the CIGAR/NM, not derived here.
pub fn is_candidate(mapq: u8, div: f64, clip_frac: f64, p: &RealignParams) -> bool {
    mapq <= p.max_mapq || clip_frac >= p.min_clip_frac || div >= p.min_div
}

/// Minimal identity floor for `realign_to_paths`: a read scoring below this against every
/// copy-path in the family fits none of them well enough to be a re-align candidate at all.
/// (Genuine novel copies pooled from low-fit reads are handled separately, by Task 4 --
/// this floor only gates obvious non-members out of the per-family re-align step.)
pub const MIN_ALN_ID: f64 = 0.5;

/// A candidate read's best fit among a family's copy-paths (Task 3), used by Task 4 to decide
/// whether re-aligning to the family beats the read's existing linear-locus alignment.
pub struct RealignHit {
    /// Index into `copy_seqs` of the best-fitting copy-path (earliest index on ties).
    pub best_copy: usize,
    /// `aln_id(read_seq, copy_seqs[best_copy])` -- the best copy-path identity.
    pub id_best: f64,
    /// `aln_id(read_seq, copy_seqs[linear_copy])` when `linear_copy` is `Some` and in range,
    /// else `0.0` -- the read's fit to the copy it would otherwise be attributed to linearly.
    pub id_linear: f64,
}

/// Re-align `read_seq` to every copy-path in `copy_seqs`, reusing `bridge_detector::aln_id` as
/// the fit score (best infix identity, tries both orientations). Returns the best-fitting copy
/// plus (for Task 4's accept comparison) the fit to `linear_copy`'s sequence, if given.
///
/// Returns `None` when `copy_seqs` is empty, or when the best fit is below `MIN_ALN_ID`: a read
/// that fits no copy-path at all is not a re-align candidate for this family.
pub fn realign_to_paths(
    read_seq: &[u8],
    copy_seqs: &[Vec<u8>],
    linear_copy: Option<usize>,
) -> Option<RealignHit> {
    if copy_seqs.is_empty() {
        return None;
    }

    let mut best_copy = 0usize;
    let mut id_best = f64::NEG_INFINITY;
    for (k, seq) in copy_seqs.iter().enumerate() {
        let id = aln_id(read_seq, seq);
        // Strict `>` keeps the earliest index on ties.
        if id > id_best {
            id_best = id;
            best_copy = k;
        }
    }

    if id_best < MIN_ALN_ID {
        return None;
    }

    let id_linear = match linear_copy {
        Some(lc) if lc < copy_seqs.len() => aln_id(read_seq, &copy_seqs[lc]),
        _ => 0.0,
    };

    Some(RealignHit { best_copy, id_best, id_linear })
}

/// Route an unmapped read to candidate families by shared canonical minimizers: compute the
/// read's minimizer value-set (`minimizers(seq, MINIMIZER_K, MINIMIZER_W)`, values only -- the
/// masked flag doesn't matter for routing) and, for each `(family_id, consensus)` pair, count
/// how many minimizer values the read shares with that family's consensus. Returns the
/// `family_id`s meeting `count >= min_shared`, sorted ascending.
pub fn route_unmapped(
    seq: &[u8],
    family_consensuses: &[(usize, Vec<u8>)],
    min_shared: usize,
) -> Vec<usize> {
    use crate::vg_family::minimizers::{MINIMIZER_K, MINIMIZER_W};

    let read_mins: HashSet<u64> =
        minimizers(seq, MINIMIZER_K, MINIMIZER_W).into_iter().map(|(v, _)| v).collect();

    let mut hits: Vec<usize> = family_consensuses
        .iter()
        .filter_map(|(family_id, consensus)| {
            let cons_mins: HashSet<u64> =
                minimizers(consensus, MINIMIZER_K, MINIMIZER_W).into_iter().map(|(v, _)| v).collect();
            let shared = read_mins.intersection(&cons_mins).count();
            (shared >= min_shared).then_some(*family_id)
        })
        .collect();

    hits.sort_unstable();
    hits
}

/// Task 4's verdict on a candidate re-alignment: either correct the read's copy attribution to
/// `Reassign(best_copy)`, or `Reject` the correction (keep whatever attribution the caller
/// already had -- linear, or none). Admission of genuinely novel copies (no existing attribution
/// at all) is a separate concern, handled by `pool_novel` here and `absent_copy::admit_candidate`
/// at wiring time -- this enum only covers correcting an EXISTING attribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RealignAction {
    Reassign(usize),
    Reject,
}

/// Decide whether `hit` (a candidate read's re-alignment result from `realign_to_paths`) beats
/// its existing linear-locus attribution `linear_copy` significantly enough to correct it.
///
/// Mirrors `copy_assign::read_copy_evidence`'s `min_p` certificate: `n_decisive` is the number of
/// read positions (out of `read_len`) that support the best copy-path over the linear locus
/// (`(id_best - id_linear) * read_len`, rounded), and `min_p = (error_rate / 3)^n_decisive` is the
/// probability that all of those decisive differences arose by sequencing error alone (an
/// `epsilon^delta` bound: each independent error has probability `error_rate / 3` of landing on
/// the specific alternate base that agrees with the best copy-path). `min_p < alpha` certifies the
/// correction; otherwise the evidence isn't strong enough to overturn the existing attribution.
///
/// No correction is needed (and none is offered) when the read's best copy-path already IS its
/// linear attribution, or when there's no decisive evidence (`n_decisive < 1`) at all.
///
/// When `linear_copy` is `None` this always `Reject`s. The `min_p` certificate here certifies a
/// *correction* -- that the best copy-path beats an existing linear attribution significantly
/// enough to overturn it. `realign_to_paths` fills `id_linear` with a `0.0` sentinel when there is
/// no linear copy to compare against, which is not a real identity and cannot serve as a
/// baseline: certifying against it would accept any moderately-well-fitting `id_best` (even a
/// near-random ~0.5 identity) as a "correction" of nothing. A read with no linear attribution at
/// all isn't a correction case -- it's handled by the separate novel-copy path (`pool_novel`).
pub fn accept_realignment(
    hit: &RealignHit,
    linear_copy: Option<usize>,
    read_len: usize,
    error_rate: f64,
    alpha: f64,
) -> RealignAction {
    let Some(linear_copy) = linear_copy else {
        return RealignAction::Reject;
    };

    if linear_copy == hit.best_copy {
        return RealignAction::Reject;
    }

    let n_decisive = ((hit.id_best - hit.id_linear) * read_len as f64).round() as i64;
    if n_decisive < 1 {
        return RealignAction::Reject;
    }

    let min_p = (error_rate / 3.0).powi(n_decisive as i32);
    if min_p < alpha {
        RealignAction::Reassign(hit.best_copy)
    } else {
        RealignAction::Reject
    }
}

/// Greedily single-linkage cluster `unfit` reads (those `realign_to_paths` matched to NO existing
/// copy -- candidate reference-absent/novel-copy material) by pairwise `aln_id >= min_id`.
///
/// Each read joins the first existing cluster whose FIRST member (the cluster's representative)
/// it matches at `>= min_id`; if it matches no cluster's representative, it starts a new
/// singleton cluster. This is a cheap O(n * clusters) pass, not full correlation clustering --
/// good enough to pool obviously-related novel-copy candidates for the Task-5 wiring, which is
/// where the actual `absent_copy::admit_candidate` admission gate (needing the genome + remap)
/// runs. Returns only clusters with `>= min_reads` members, as index vectors into `unfit`.
pub fn pool_novel(unfit: &[(String, Vec<u8>)], min_id: f64, min_reads: usize) -> Vec<Vec<usize>> {
    let mut clusters: Vec<Vec<usize>> = Vec::new();

    for (i, (_, seq)) in unfit.iter().enumerate() {
        let mut joined = false;
        for cluster in clusters.iter_mut() {
            let rep = cluster[0];
            if aln_id(seq, &unfit[rep].1) >= min_id {
                cluster.push(i);
                joined = true;
                break;
            }
        }
        if !joined {
            clusters.push(vec![i]);
        }
    }

    clusters.retain(|c| c.len() >= min_reads);
    clusters
}

/// One per-read decision from the per-family re-align supplement (Task 5), emitted verbatim as a row of
/// `<out>.vg_realign.tsv` by the `copy_assign` binary.
#[derive(Debug, Clone, PartialEq)]
pub struct RealignRecord {
    pub read_name: String,
    /// `"reassigned"` (a significant correction to a different copy), `"rejected"` (a candidate that
    /// re-aligned but didn't clear Task 4's significance certificate), or `"novel-candidate"` (fits no
    /// existing copy-path at all -- `Task 4`'s `pool_novel`/`absent_copy` admission gate is the separate,
    /// out-of-scope next step for these).
    pub action: String,
    /// The copy index the record concerns: `best_copy` when `action == "reassigned"`, else `-1` (no
    /// correction target -- `"rejected"` keeps the read's existing linear attribution, `"novel-candidate"`
    /// has no copy-path fit at all).
    pub target_copy: i64,
    /// `RealignHit::id_best` (0.0 for `"novel-candidate"`, where `realign_to_paths` found no fit at all so
    /// no `id_best` was computed).
    pub id_best: f64,
    /// The copy index the read's own linear (BAM-coordinate) alignment placed it on within this family, or
    /// `-1` if none of the family's copies overlap the read's aligned span (`best_overlap_copy` returned
    /// `None`).
    pub linear_copy: i64,
}

/// VG re-align END-TO-END plan, Task 2: `apply_realign`'s output.
///
/// `admitted` reference-absent copies are NOT produced here -- see the field doc on
/// `novel_pools` and `apply_realign`'s doc for why (this stays genome-free/testable; admission
/// is the follow-up wiring task's job).
pub struct RealignApply {
    /// `read_index -> (new_copy_idx, path_obs)`: a significant correction (Task 4's
    /// `RealignAction::Reassign`) to a DIFFERENT copy than the read's existing linear
    /// attribution, plus the read's base at each of the new copy's PSV columns (from
    /// `align_traceback` + `path_obs_at` against the new copy's consensus).
    pub corrected: std::collections::HashMap<usize, (usize, Vec<Option<u8>>)>,
    /// Clusters (each `>= rp.min_reads` members) of read INDICES (into `bam_reads`) that fit no
    /// existing copy-path at all (`realign_to_paths` returned `None`) but are mutually similar
    /// enough (`pool_novel`, `min_id ~= 0.9`) to be candidate novel/reference-absent copies. Not
    /// yet admitted -- the wiring task turns each pool into a `CollapsedCandidate` and runs it
    /// through `absent_copy::admit_candidate` with the real genome + remap.
    pub novel_pools: Vec<Vec<usize>>,
    /// One record per candidate read processed (`"reassigned"`, `"rejected"`, or
    /// `"novel-candidate"`), same shape as `run_family_realign`'s output.
    pub records: Vec<RealignRecord>,
}

/// Task 5: run the VG re-align supplement (Tasks 1/3/4) over ONE co-located family's reads.
///
/// For every non-supplementary read in `bam_reads`: compute `mapq` (`read.mapq`), `clip_frac` (total
/// soft-clipped CIGAR bases / read length -- `read.seq` keeps soft-clips per
/// `aligned_read_from_record`'s doc, so this is exact), and `div` (`read.de`, minimap2's `de:f` gap-
/// compressed per-base divergence tag -- already parsed onto `BamRead`, so no NM/CIGAR recomputation is
/// needed). `linear_copy` is the copy index the read's own linear alignment overlaps most in this family
/// (`best_overlap_copy`, the SAME greatest-ref-overlap rule `assign_family_detailed` uses internally to
/// seed each read's `mapped_copy`).
///
/// Reads failing `is_candidate` produce NO record at all (clean primary fit -- nothing to reconsider).
/// Candidates are re-aligned to every copy's spliced consensus (`DenovoTranscript::seq`, already spliced
/// by the assembly stage -- no re-splicing needed) via `realign_to_paths`; `None` (fits no copy-path at
/// all) is tagged `"novel-candidate"`. A hit is run through `accept_realignment`'s significance
/// certificate: `Reassign` -> `"reassigned"`, `Reject` -> `"rejected"`.
///
/// REPORT-ONLY / ADDITIVE: this function only classifies reads into a decision log. It does not mutate
/// `bam_reads`/`copies`, does not feed a `"reassigned"` verdict back into the EM/PSV assignment, and does
/// not admit `"novel-candidate"` reads into the copy set -- that deeper wiring (`pool_novel` +
/// `absent_copy::admit_candidate`) is an explicit follow-up, out of scope here.
pub fn run_family_realign(
    bam_reads: &[BamRead],
    copies: &[DenovoTranscript],
    params: &RealignParams,
    error_rate: f64,
    alpha: f64,
) -> Vec<RealignRecord> {
    let copy_seqs: Vec<Vec<u8>> = copies.iter().map(|c| c.seq.clone()).collect();
    let copy_refs: Vec<&DenovoTranscript> = copies.iter().collect();

    let mut out = Vec::new();
    for br in bam_reads {
        if br.is_supplementary {
            continue;
        }
        let read_len = br.read.seq.len();
        let clip: u64 = br.read.cigar.iter().filter(|&&(op, _)| op == 'S').map(|&(_, n)| n).sum();
        let clip_frac = if read_len > 0 { clip as f64 / read_len as f64 } else { 0.0 };
        let div = br.de as f64;

        if !is_candidate(br.mapq, div, clip_frac, params) {
            continue;
        }

        let linear_copy = best_overlap_copy(&br.read, &copy_refs);
        let linear_copy_i64 = linear_copy.map(|c| c as i64).unwrap_or(-1);

        let record = match realign_to_paths(&br.read.seq, &copy_seqs, linear_copy) {
            None => RealignRecord {
                read_name: br.name.clone(),
                action: "novel-candidate".to_string(),
                target_copy: -1,
                id_best: 0.0,
                linear_copy: linear_copy_i64,
            },
            Some(hit) => {
                let id_best = hit.id_best;
                match accept_realignment(&hit, linear_copy, read_len, error_rate, alpha) {
                    RealignAction::Reassign(best_copy) => RealignRecord {
                        read_name: br.name.clone(),
                        action: "reassigned".to_string(),
                        target_copy: best_copy as i64,
                        id_best,
                        linear_copy: linear_copy_i64,
                    },
                    RealignAction::Reject => RealignRecord {
                        read_name: br.name.clone(),
                        action: "rejected".to_string(),
                        target_copy: -1,
                        id_best,
                        linear_copy: linear_copy_i64,
                    },
                }
            }
        };
        out.push(record);
    }
    out
}

/// VG re-align END-TO-END plan, Task 2: apply the per-read decisions (Tasks 1/3/4) into a
/// ready-to-consume correction map + novel-copy candidate pools, over reads spanning potentially
/// SEVERAL families at once (unlike `run_family_realign`'s one-family driver).
///
/// `copies`/`copy_seqs` are parallel (one spliced consensus per copy, `copy_seqs[k] ==
/// copies[k].seq` is the expected caller invariant but only `copy_seqs` is actually read here --
/// `copies` is carried for callers/future use, e.g. locus metadata alongside the correction map).
/// `psv_pos_per_copy[k]` are the family's PSV positions in `copy_seqs[k]`'s consensus coordinates.
/// `linear_copy_of[i]` is read `i`'s existing linear-locus copy attribution (or `None`), parallel
/// to `bam_reads`.
///
/// CORRECTIONS (must-have): a candidate read (`is_candidate`) that re-aligns to a copy-path
/// (`realign_to_paths`) and clears `accept_realignment`'s significance certificate
/// (`RealignAction::Reassign`) is entered into `corrected[read_index] = (new_copy, path_obs)`,
/// where `path_obs` is the read's base at each of the new copy's PSV columns
/// (`align_traceback` + `path_obs_at` against `copy_seqs[new_copy]`).
///
/// ADMISSIONS (best-effort/mechanism-only): a candidate read that fits NO copy-path at all
/// (`realign_to_paths` returns `None`) is pooled with other such "unfit" reads by
/// `pool_novel`; `novel_pools` holds the resulting clusters (mapped back to indices into
/// `bam_reads`) with `>= rp.min_reads` members. This is genome-free and does NOT run
/// `absent_copy::admit_candidate` -- turning a pool into an admitted reference-absent copy needs
/// the real genome + remap, which is the follow-up wiring task's job, not this one's. Real yield
/// here is data-limited (the O4 divergent frontier): most families will produce zero pools.
///
/// Reads failing `is_candidate` are skipped entirely (no record, no correction, no pooling) --
/// a clean primary fit has nothing to reconsider. Supplementary alignments (`is_supplementary`)
/// are also skipped, mirroring `run_family_realign`.
pub fn apply_realign(
    bam_reads: &[BamRead],
    copies: &[DenovoTranscript],
    copy_seqs: &[Vec<u8>],
    psv_pos_per_copy: &[Vec<usize>],
    linear_copy_of: &[Option<usize>],
    rp: &RealignParams,
    error_rate: f64,
    alpha: f64,
) -> RealignApply {
    let _ = copies; // parallel to copy_seqs; not read directly here (see doc).

    let mut corrected = std::collections::HashMap::new();
    let mut records = Vec::new();
    let mut unfit: Vec<(String, Vec<u8>)> = Vec::new();
    let mut unfit_idx: Vec<usize> = Vec::new();

    for (i, br) in bam_reads.iter().enumerate() {
        if br.is_supplementary {
            continue;
        }
        let read_len = br.read.seq.len();
        let clip: u64 = br.read.cigar.iter().filter(|&&(op, _)| op == 'S').map(|&(_, n)| n).sum();
        let clip_frac = if read_len > 0 { clip as f64 / read_len as f64 } else { 0.0 };
        let div = br.de as f64;

        if !is_candidate(br.mapq, div, clip_frac, rp) {
            continue;
        }

        let linear_copy = linear_copy_of.get(i).copied().flatten();
        let linear_copy_i64 = linear_copy.map(|c| c as i64).unwrap_or(-1);

        match realign_to_paths(&br.read.seq, copy_seqs, linear_copy) {
            None => {
                unfit.push((br.name.clone(), br.read.seq.clone()));
                unfit_idx.push(i);
                records.push(RealignRecord {
                    read_name: br.name.clone(),
                    action: "novel-candidate".to_string(),
                    target_copy: -1,
                    id_best: 0.0,
                    linear_copy: linear_copy_i64,
                });
            }
            Some(hit) => {
                let id_best = hit.id_best;
                match accept_realignment(&hit, linear_copy, read_len, error_rate, alpha) {
                    RealignAction::Reassign(best_copy) => {
                        // C1: orient the read to the copy's own (transcription-strand) frame
                        // before the traceback -- a `-`-strand copy's consensus is the reverse
                        // complement of forward-genome, so a literal forward alignment here would
                        // produce garbage/`None` `path_obs`.
                        let oriented = orient_for_copy(&br.read.seq, &copy_seqs[best_copy]);
                        let map = align_traceback(&oriented, &copy_seqs[best_copy]);
                        let obs = path_obs_at(&map, &psv_pos_per_copy[best_copy], &oriented);
                        corrected.insert(i, (best_copy, obs));
                        records.push(RealignRecord {
                            read_name: br.name.clone(),
                            action: "reassigned".to_string(),
                            target_copy: best_copy as i64,
                            id_best,
                            linear_copy: linear_copy_i64,
                        });
                    }
                    RealignAction::Reject => {
                        records.push(RealignRecord {
                            read_name: br.name.clone(),
                            action: "rejected".to_string(),
                            target_copy: -1,
                            id_best,
                            linear_copy: linear_copy_i64,
                        });
                    }
                }
            }
        }
    }

    let clusters = pool_novel(&unfit, 0.9, rp.min_reads);
    let novel_pools: Vec<Vec<usize>> =
        clusters.into_iter().map(|c| c.into_iter().map(|j| unfit_idx[j]).collect()).collect();

    RealignApply { corrected, novel_pools, records }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_family::bridge_detector::hw_distance;
    use crate::vg_family::copy_split::AlignedRead;

    #[test]
    fn is_candidate_flags_poor_fit() {
        let p = RealignParams::default();

        // low MAPQ alone -> true
        assert!(is_candidate(5, 0.0, 0.0, &p));
        // high clip alone -> true
        assert!(is_candidate(60, 0.0, 0.30, &p));
        // high divergence alone -> true
        assert!(is_candidate(60, 0.08, 0.0, &p));
        // clean read: high MAPQ, no clip, low div -> false
        assert!(!is_candidate(60, 0.0, 0.0, &p));

        // boundary: mapq == max_mapq is still <= -> true
        assert!(is_candidate(20, 0.0, 0.0, &p));
        // boundary: just under both clip and div thresholds, and mapq above max -> false
        assert!(!is_candidate(21, 0.049, 0.19, &p));
    }

    /// Deterministic pseudo-random ACGT sequence generator (xorshift64), so test sequences are
    /// reproducible without hand-typing long strings or depending on real fixture data.
    fn pseudo_seq(seed: u64, len: usize) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        // splitmix64-style mix so nearby/small seeds (1, 2, 3, ...) don't collapse to the same
        // or correlated xorshift states (plain `seed | 1` made seeds 2 and 3 identical).
        let mut state =
            seed.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(0x2545F4914F6CDD1D) | 1;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                bases[(state % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn realign_picks_best_copy_path() {
        let copy0 = pseudo_seq(1, 60);
        let copy1 = pseudo_seq(2, 60);
        let copy2 = pseudo_seq(3, 60);
        let copy_seqs = vec![copy0.clone(), copy1.clone(), copy2.clone()];

        // Read == copy 1 exactly -> best_copy == 1, id_best ~1.0.
        let hit = realign_to_paths(&copy1, &copy_seqs, Some(0))
            .expect("exact copy-path match must be a candidate");
        assert_eq!(hit.best_copy, 1);
        assert!(hit.id_best > 0.99, "id_best = {}", hit.id_best);
        // id_linear (fit to copy 0, the "linear locus" copy) must be strictly lower than the
        // true best-copy fit -- the read really belongs to copy 1, not copy 0.
        assert!(
            hit.id_linear < hit.id_best,
            "id_linear = {} should be < id_best = {}",
            hit.id_linear,
            hit.id_best
        );

        // Read == copy 2 with a couple of substitutions -> still best_copy == 2, id_best < 1.0.
        let mut mutated = copy2.clone();
        // Flip two bases to a base guaranteed different from the original.
        for &pos in &[10usize, 40usize] {
            let orig = mutated[pos];
            mutated[pos] = [b'A', b'C', b'G', b'T'].into_iter().find(|&b| b != orig).unwrap();
        }
        let hit2 = realign_to_paths(&mutated, &copy_seqs, None)
            .expect("near-exact copy-path match must be a candidate");
        assert_eq!(hit2.best_copy, 2);
        assert!(hit2.id_best < 1.0, "id_best = {}", hit2.id_best);
        assert!(hit2.id_best > 0.9, "id_best = {} should still be a strong fit", hit2.id_best);
        // No linear_copy given -> id_linear is the documented 0.0 sentinel.
        assert_eq!(hit2.id_linear, 0.0);
    }

    #[test]
    fn realign_returns_none_for_nonmember() {
        let copy_seqs = vec![pseudo_seq(1, 60), pseudo_seq(2, 60), pseudo_seq(3, 60)];
        // A read unrelated to any copy-path (seed chosen so its best infix identity to all
        // three copies lands below MIN_ALN_ID; random same-length sequences under free-end-gap
        // edit distance land around ~0.5 identity by chance, so this isn't a free lunch).
        let read = pseudo_seq(3130, 60);
        let hit = realign_to_paths(&read, &copy_seqs, None);
        match hit {
            None => {}
            Some(h) => panic!("expected None for a non-member read, got id_best = {}", h.id_best),
        }

        // Empty copy_seqs -> always None regardless of the read.
        assert!(realign_to_paths(&read, &[], None).is_none());
    }

    #[test]
    fn route_unmapped_matches_by_minimizers() {
        // Two family consensuses, long enough (>= MINIMIZER_K + several windows) to produce
        // multiple minimizers.
        let consensus_a = pseudo_seq(10, 150);
        let consensus_b = pseudo_seq(20, 150);

        // A read that's an exact interior substring of family A's consensus -- interior so its
        // minimizer windows are fully contained in windows also present in the full consensus,
        // guaranteeing shared minimizer VALUES (not just positions).
        let read: Vec<u8> = consensus_a[40..100].to_vec();

        let families = vec![(0usize, consensus_a.clone()), (1usize, consensus_b.clone())];

        // Sanity: the read and consensus B should share few/no minimizers (unrelated random
        // sequences), while read and consensus A share several -- confirm before picking
        // min_shared so the test isn't tuned to a lucky threshold.
        use crate::vg_family::minimizers::{MINIMIZER_K, MINIMIZER_W};
        let mins_of = |s: &[u8]| -> HashSet<u64> {
            minimizers(s, MINIMIZER_K, MINIMIZER_W).into_iter().map(|(v, _)| v).collect()
        };
        let read_mins = mins_of(&read);
        let a_mins = mins_of(&consensus_a);
        let b_mins = mins_of(&consensus_b);
        let shared_a = read_mins.intersection(&a_mins).count();
        let shared_b = read_mins.intersection(&b_mins).count();
        assert!(shared_a >= 2, "expected several shared minimizers with A, got {shared_a}");
        assert!(shared_b < shared_a, "B should share fewer minimizers than A ({shared_b} vs {shared_a})");

        let min_shared = shared_b + 1; // strictly above B's count, at or below A's count
        assert!(min_shared <= shared_a, "min_shared {min_shared} must still be reachable by A");

        let routed = route_unmapped(&read, &families, min_shared);
        assert_eq!(routed, vec![0], "expected only family A routed, got {routed:?}");
    }

    #[test]
    fn realign_params_defaults() {
        let p = RealignParams::default();
        assert_eq!(p.max_mapq, 20);
        assert_eq!(p.min_clip_frac, 0.20);
        assert_eq!(p.min_div, 0.05);
        assert_eq!(p.min_reads, 3);
    }

    #[test]
    fn accept_significant_reassigns() {
        // id diff 0.14 * read_len 1000 -> n_decisive = 140 decisive positions favoring copy 2
        // over the read's current linear attribution (copy 0). min_p = (0.003/3)^140 is
        // astronomically small (<< alpha = 1e-3) -- certifies the correction.
        let hit = RealignHit { best_copy: 2, id_best: 0.99, id_linear: 0.85 };
        let n_decisive = ((hit.id_best - hit.id_linear) * 1000.0).round() as i64;
        assert_eq!(n_decisive, 140);

        let action = accept_realignment(&hit, Some(0), 1000, 0.003, 1e-3);
        assert_eq!(action, RealignAction::Reassign(2));
    }

    #[test]
    fn accept_marginal_rejects() {
        // id diff 0.001 * read_len 1000 -> n_decisive = 1 (rounds exactly to 1). With a
        // deliberately high error_rate = 0.05, min_p = (0.05/3)^1 = 0.01666... which is >= alpha
        // = 1e-3 -- a single decisive position isn't enough evidence to overturn the read's
        // existing linear attribution (copy 0) under this noisy an error model.
        let hit = RealignHit { best_copy: 2, id_best: 0.90, id_linear: 0.899 };
        let n_decisive = ((hit.id_best - hit.id_linear) * 1000.0).round() as i64;
        assert_eq!(n_decisive, 1);
        let min_p = (0.05_f64 / 3.0).powi(1);
        assert!(min_p >= 1e-3, "min_p = {min_p} should be >= alpha");

        let action = accept_realignment(&hit, Some(0), 1000, 0.05, 1e-3);
        assert_eq!(action, RealignAction::Reject);
    }

    #[test]
    fn accept_zero_decisive_rejects() {
        // id_best == id_linear -> n_decisive = 0 -> no decisive evidence at all -> Reject,
        // regardless of how permissive alpha/error_rate are.
        let hit = RealignHit { best_copy: 2, id_best: 0.95, id_linear: 0.95 };
        let action = accept_realignment(&hit, Some(0), 1000, 0.003, 1.0);
        assert_eq!(action, RealignAction::Reject);
    }

    #[test]
    fn accept_best_equals_linear_rejects() {
        // best_copy already IS the read's linear attribution -- no correction needed even
        // though id_best/id_linear here would otherwise look decisive.
        let hit = RealignHit { best_copy: 2, id_best: 0.99, id_linear: 0.10 };
        let action = accept_realignment(&hit, Some(2), 1000, 0.003, 1e-3);
        assert_eq!(action, RealignAction::Reject);
    }

    #[test]
    fn accept_none_linear_rejects() {
        // No existing linear attribution at all (unmapped read routed by Task 3) -- id_linear's
        // 0.0 sentinel is not a real baseline, so there is nothing to "correct" against. Even a
        // high id_best (0.99) must Reject here, not Reassign off the meaningless zero baseline;
        // genuinely unattributed reads are handled by the separate novel-copy path.
        let hit = RealignHit { best_copy: 1, id_best: 0.99, id_linear: 0.0 };
        let action = accept_realignment(&hit, None, 1000, 0.003, 1e-3);
        assert_eq!(action, RealignAction::Reject);
    }

    #[test]
    fn pool_novel_clusters_unfit() {
        // 3 mutually near-identical reads (one exact copy + two 1-base mutants of it) plus 1
        // unrelated random read. min_id = 0.9, min_reads = 3.
        let base = pseudo_seq(100, 80);
        let mut mut1 = base.clone();
        mut1[5] = [b'A', b'C', b'G', b'T'].into_iter().find(|&b| b != mut1[5]).unwrap();
        let mut mut2 = base.clone();
        mut2[60] = [b'A', b'C', b'G', b'T'].into_iter().find(|&b| b != mut2[60]).unwrap();
        // Unrelated: different seed, same length, chosen so its identity to `base` lands well
        // under 0.9 (random same-length sequences under free-end-gap identity hover ~0.5).
        let unrelated = pseudo_seq(9999, 80);
        assert!(
            aln_id(&base, &unrelated) < 0.9,
            "fixture assumption broken: unrelated read too similar to base"
        );

        let unfit: Vec<(String, Vec<u8>)> = vec![
            ("r0".to_string(), base),
            ("r1".to_string(), mut1),
            ("r2".to_string(), mut2),
            ("r3".to_string(), unrelated),
        ];

        let clusters = pool_novel(&unfit, 0.9, 3);
        assert_eq!(clusters.len(), 1, "expected exactly one cluster to survive min_reads, got {clusters:?}");
        let mut got = clusters[0].clone();
        got.sort_unstable();
        assert_eq!(got, vec![0, 1, 2], "expected the 3 mutually-similar reads clustered together");
    }

    #[test]
    fn pool_novel_below_min_reads_dropped() {
        let a = pseudo_seq(1, 60);
        let b = a.clone();
        let unfit: Vec<(String, Vec<u8>)> = vec![("a".to_string(), a), ("b".to_string(), b)];
        // Only 2 mutually-identical reads, but min_reads = 3 -> nothing survives.
        let clusters = pool_novel(&unfit, 0.9, 3);
        assert!(clusters.is_empty(), "expected no clusters to meet min_reads = 3, got {clusters:?}");
    }

    /// Minimal `DenovoTranscript` builder for the `run_family_realign` tests -- unspliced (no introns),
    /// mirroring how `copy_assign_pipeline`'s own tests construct copies.
    fn transcript(tid: &str, chrom: &str, start: u64, seq: Vec<u8>) -> DenovoTranscript {
        let end = start + seq.len() as u64;
        DenovoTranscript { tid: tid.to_string(), chrom: chrom.to_string(), start, end, n_reads: 10, strand: '+', introns: vec![], seq }
    }

    fn bam_read(name: &str, chrom: &str, ref_start: u64, seq: Vec<u8>, mapq: u8, de: f32) -> BamRead {
        let len = seq.len() as u64;
        BamRead {
            chrom: chrom.to_string(),
            read: AlignedRead { ref_start, cigar: vec![('M', len)], seq, qual: vec![] },
            mapq,
            name: name.to_string(),
            as_score: 0,
            de,
            is_supplementary: false,
            is_secondary: false,
        }
    }

    #[test]
    fn run_family_realign_reassigns_misplaced_low_mapq_read() {
        let copy0_seq = pseudo_seq(1, 200);
        let copy1_seq = pseudo_seq(2, 200);
        // A substitution-variant of copy 1 (2 flipped bases) -- still >99% identical to copy 1, and only
        // ~50% identical to the unrelated copy 0 (random same-length sequences).
        let mut read_seq = copy1_seq.clone();
        for &pos in &[20usize, 100usize] {
            let orig = read_seq[pos];
            read_seq[pos] = [b'A', b'C', b'G', b'T'].into_iter().find(|&b| b != orig).unwrap();
        }

        // Copy 0 at chr1:0-200, copy 1 at chr1:5000-5200 -- distinct, non-overlapping loci.
        let copy0 = transcript("copy0", "chr1", 0, copy0_seq);
        let copy1 = transcript("copy1", "chr1", 5000, copy1_seq);
        let copies = vec![copy0, copy1];

        // The read is "linearly placed" (BAM ref_start) inside copy 0's span -- so its coordinate-overlap
        // locus (linear_copy) is copy 0 -- but its SEQUENCE is really copy 1's, and its MAPQ is low
        // (ambiguous multimapper), so it's a Task-1 candidate.
        let br = bam_read("readA", "chr1", 0, read_seq, 3, 0.0);

        let records = run_family_realign(&[br], &copies, &RealignParams::default(), 0.003, 1e-3);
        assert_eq!(records.len(), 1, "expected exactly one record, got {records:?}");
        let rec = &records[0];
        assert_eq!(rec.read_name, "readA");
        assert_eq!(rec.linear_copy, 0, "read's linear (BAM-coordinate) placement must be copy 0");
        assert_eq!(rec.action, "reassigned");
        assert_eq!(rec.target_copy, 1, "the read's true best-fit copy is copy 1");
        assert!(rec.id_best > 0.9, "id_best = {} should be a strong fit to copy 1", rec.id_best);
    }

    #[test]
    fn run_family_realign_clean_read_produces_no_record() {
        let copy0_seq = pseudo_seq(1, 200);
        let copy1_seq = pseudo_seq(2, 200);
        let copy0 = transcript("copy0", "chr1", 0, copy0_seq.clone());
        let copy1 = transcript("copy1", "chr1", 5000, copy1_seq);
        let copies = vec![copy0, copy1];

        // A read that is exactly copy 0's sequence, high MAPQ, no clipping, on copy 0's own locus --
        // a clean primary fit, not a Task-1 candidate at all.
        let br = bam_read("readB", "chr1", 0, copy0_seq, 60, 0.0);

        let records = run_family_realign(&[br], &copies, &RealignParams::default(), 0.003, 1e-3);
        assert!(records.is_empty(), "a clean high-MAPQ read must produce no record, got {records:?}");
    }

    // -----------------------------------------------------------------------------------------
    // Task 1 (end-to-end plan): align_traceback + path_obs_at
    // -----------------------------------------------------------------------------------------

    /// Sum of mismatched-diagonal columns + gap columns (either side) in an alignment map --
    /// this must equal `hw_distance`'s edit distance for any valid traceback of that DP.
    fn edits_in_map(align_map: &[(Option<usize>, Option<usize>)], query: &[u8], target: &[u8]) -> usize {
        align_map
            .iter()
            .filter(|&&(qi, ti)| match (qi, ti) {
                (Some(q), Some(t)) => query[q] != target[t],
                (Some(_), None) | (None, Some(_)) => true,
                (None, None) => panic!("align_traceback must never emit a (None, None) column"),
            })
            .count()
    }

    #[test]
    fn traceback_edit_distance_matches_hw() {
        let cases: Vec<(&[u8], &[u8])> = vec![
            // exact infix: ACGT occurs verbatim inside TTACGTGG -> dist 0.
            (b"ACGT", b"TTACGTGG"),
            // 1 substitution: ACGT vs an infix that differs at one base (ACCT inside TT-ACCT-GG).
            (b"ACGT", b"TTACCTGG"),
            // 1 insertion (relative to target): query has an extra base not in any target infix.
            (b"ACGGT", b"TTACGTGG"),
            // 1 deletion (relative to target): query is missing a base present in the target infix.
            (b"ACT", b"TTACGTGG"),
            // longer query with 2 edits scattered through it.
            (b"ACGTACGTAC", b"TTACGTACCTAGGG"),
        ];

        for (query, target) in cases {
            let expected = hw_distance(query, target);
            let map = align_traceback(query, target);
            let got = edits_in_map(&map, query, target);
            assert_eq!(
                got, expected,
                "query={:?} target={:?}: traceback edits {got} != hw_distance {expected}",
                std::str::from_utf8(query).unwrap(),
                std::str::from_utf8(target).unwrap()
            );

            // Sanity: the aligned columns must walk query positions 0..query.len() in order (every
            // query base consumed exactly once, monotonically), since HW gives free end-gaps only on
            // the TARGET, not the query.
            let q_positions: Vec<usize> =
                map.iter().filter_map(|&(qi, _)| qi).collect();
            let expected_q_positions: Vec<usize> = (0..query.len()).collect();
            assert_eq!(
                q_positions, expected_q_positions,
                "every query position must appear exactly once, in order"
            );
        }
    }

    // -----------------------------------------------------------------------------------------
    // VG re-align END-TO-END plan, Task 2: apply_realign orchestrator
    // -----------------------------------------------------------------------------------------

    #[test]
    fn apply_corrects_mismapped_read() {
        // Two distinct-consensus copies; copy 1 has a distinguishing base at consensus position
        // `p` relative to copy 0 (guaranteed distinct at that column by construction).
        let copy0_seq = pseudo_seq(1, 200);
        let mut copy1_seq = pseudo_seq(2, 200);
        let p = 50usize;
        // Force column p to differ between the two copies (pseudo_seq draws from independent
        // seeds so this already usually holds, but force it so the test isn't seed-lucky).
        if copy1_seq[p] == copy0_seq[p] {
            copy1_seq[p] = [b'A', b'C', b'G', b'T'].into_iter().find(|&b| b != copy0_seq[p]).unwrap();
        }
        let copy0 = transcript("copy0", "chr1", 0, copy0_seq.clone());
        let copy1 = transcript("copy1", "chr1", 5000, copy1_seq.clone());
        let copies = vec![copy0, copy1];
        let copy_seqs = vec![copy0_seq, copy1_seq.clone()];
        let psv_pos_per_copy = vec![vec![p], vec![p]];

        // Read == copy 1's consensus exactly, but low MAPQ and linearly attributed to copy 0
        // (its BAM placement) -- a Task-1 candidate whose true source is copy 1.
        let br = bam_read("readA", "chr1", 0, copy1_seq.clone(), 3, 0.0);
        let linear_copy_of = vec![Some(0usize)];

        let out = apply_realign(
            &[br],
            &copies,
            &copy_seqs,
            &psv_pos_per_copy,
            &linear_copy_of,
            &RealignParams::default(),
            0.003,
            1e-3,
        );

        assert_eq!(out.corrected.len(), 1, "expected exactly one correction, got {:?}", out.corrected);
        let (new_copy, obs) = out.corrected.get(&0).expect("read 0 must be corrected");
        assert_eq!(*new_copy, 1, "read's true best-fit copy is copy 1");
        assert_eq!(obs, &vec![Some(copy1_seq[p])], "obs at the PSV column must equal copy 1's base");
        assert!(out.novel_pools.is_empty());
        assert_eq!(
            out.records.iter().filter(|r| r.action == "reassigned").count(),
            1,
            "expected one reassigned record, got {:?}",
            out.records
        );
    }

    /// C1: a `-`-strand copy's `copy_seqs[k]` is TRANSCRIPTION strand, i.e. the reverse complement
    /// of the forward-genome sequence at that locus. A read is always FORWARD-GENOME (`read.seq`,
    /// per BAM convention), so a read that is truly this copy's source material arrives as
    /// `revcomp(copy_seq)`, not `copy_seq` itself. `realign_to_paths`/`aln_id` already handle this
    /// (they try both orientations for the identity SCORE), so the correction is still detected --
    /// but before the C1 fix, `apply_realign`'s `align_traceback`/`path_obs_at` aligned the raw
    /// forward `read.seq` literally against `copy_seq`, producing a ~0.5-identity garbage
    /// alignment and a wrong/`None` `path_obs` at the PSV column. This test pins that the
    /// corrected `path_obs` instead equals copy 1's own TRANSCRIPTION-strand allele -- exactly what
    /// `copy_assign_pipeline::fill_psv_obs`'s per-base `rc_base` for `-` copies would produce.
    #[test]
    fn apply_realign_strand_orients_path_obs_for_minus_copy() {
        let copy0_seq = pseudo_seq(1, 200);
        let mut copy1_seq = pseudo_seq(2, 200); // TRANSCRIPTION-strand consensus of a '-'-strand copy
        let p = 50usize;
        if copy1_seq[p] == copy0_seq[p] {
            copy1_seq[p] = [b'A', b'C', b'G', b'T'].into_iter().find(|&b| b != copy0_seq[p]).unwrap();
        }
        let copy0 = transcript("copy0", "chr1", 0, copy0_seq.clone());
        let mut copy1 = transcript("copy1", "chr1", 5000, copy1_seq.clone());
        copy1.strand = '-';
        let copies = vec![copy0, copy1];
        let copy_seqs = vec![copy0_seq, copy1_seq.clone()];
        let psv_pos_per_copy = vec![vec![p], vec![p]];

        // The read is FORWARD-GENOME: the reverse complement of copy 1's transcription-strand
        // consensus. Low MAPQ (a Task-1 candidate) and linearly misattributed to copy 0.
        let read_seq = crate::vg_family::bridge_detector::revcomp(&copy1_seq);
        let br = bam_read("readA", "chr1", 0, read_seq, 3, 0.0);
        let linear_copy_of = vec![Some(0usize)];

        let out = apply_realign(
            &[br],
            &copies,
            &copy_seqs,
            &psv_pos_per_copy,
            &linear_copy_of,
            &RealignParams::default(),
            0.003,
            1e-3,
        );

        assert_eq!(out.corrected.len(), 1, "expected exactly one correction, got {:?}", out.corrected);
        let (new_copy, obs) = out.corrected.get(&0).expect("read 0 must be corrected");
        assert_eq!(*new_copy, 1, "read's true best-fit copy is copy 1 (the '-'-strand copy)");
        assert_eq!(
            obs,
            &vec![Some(copy1_seq[p])],
            "path_obs at the PSV column must equal copy 1's TRANSCRIPTION-strand allele, not a \
             garbage/wrong-orientation base"
        );
    }

    #[test]
    fn apply_pools_novel_unfit() {
        // 3 mutually near-identical reads that match NEITHER copy (id_best < MIN_ALN_ID so
        // realign_to_paths returns None for all of them) -- candidate novel-copy material.
        let copy0_seq = pseudo_seq(1, 80);
        let copy1_seq = pseudo_seq(2, 80);
        let copy0 = transcript("copy0", "chr1", 0, copy0_seq.clone());
        let copy1 = transcript("copy1", "chr1", 5000, copy1_seq.clone());
        let copies = vec![copy0, copy1];
        let copy_seqs = vec![copy0_seq.clone(), copy1_seq.clone()];
        let psv_pos_per_copy = vec![vec![], vec![]];

        let novel_base = pseudo_seq(69, 80);
        assert!(
            aln_id(&novel_base, &copy0_seq) < MIN_ALN_ID && aln_id(&novel_base, &copy1_seq) < MIN_ALN_ID,
            "fixture assumption broken: novel_base must fit neither existing copy"
        );
        // Exact clones of `novel_base` (not single-base mutants): with a random ~80bp sequence
        // sitting right at the ~0.5 "no better than chance" identity floor against the existing
        // copies, even a 1-base mutation can nudge `aln_id` across the `MIN_ALN_ID` boundary in
        // either direction (edit-distance realignment isn't strictly monotone in Hamming
        // distance). Cloning keeps every read's fit to the existing copies IDENTICAL to the
        // already-asserted `novel_base` fit, while still being "mutually near-identical"
        // (identity 1.0) for `pool_novel`'s `>= min_id` clustering.
        let novel1 = novel_base.clone();
        let novel2 = novel_base.clone();

        // A clean read on its correct copy (high MAPQ, no clip, low div) -- not a candidate at
        // all, must produce no record and not enter the pool.
        let clean = bam_read("clean", "chr1", 0, copy0_seq, 60, 0.0);
        let n0 = bam_read("n0", "chr1", 0, novel_base, 3, 0.0);
        let n1 = bam_read("n1", "chr1", 0, novel1, 3, 0.0);
        let n2 = bam_read("n2", "chr1", 0, novel2, 3, 0.0);

        let bam_reads = vec![clean, n0, n1, n2];
        let linear_copy_of = vec![Some(0usize), None, None, None];

        let out = apply_realign(
            &bam_reads,
            &copies,
            &copy_seqs,
            &psv_pos_per_copy,
            &linear_copy_of,
            &RealignParams::default(),
            0.003,
            1e-3,
        );

        assert!(out.corrected.is_empty(), "no corrections expected, got {:?}", out.corrected);
        assert_eq!(out.novel_pools.len(), 1, "expected exactly one novel pool, got {:?}", out.novel_pools);
        let mut got = out.novel_pools[0].clone();
        got.sort_unstable();
        assert_eq!(got, vec![1, 2, 3], "expected read indices 1,2,3 (the 3 novel reads) pooled");

        let novel_records: Vec<_> = out.records.iter().filter(|r| r.action == "novel-candidate").collect();
        assert_eq!(novel_records.len(), 3, "expected one novel-candidate record per unfit read");
        assert!(
            out.records.iter().all(|r| r.read_name != "clean"),
            "the clean read must produce no record at all"
        );
    }

    #[test]
    fn apply_clean_family_noop() {
        let copy0_seq = pseudo_seq(1, 200);
        let copy1_seq = pseudo_seq(2, 200);
        let copy0 = transcript("copy0", "chr1", 0, copy0_seq.clone());
        let copy1 = transcript("copy1", "chr1", 5000, copy1_seq.clone());
        let copies = vec![copy0, copy1];
        let copy_seqs = vec![copy0_seq.clone(), copy1_seq.clone()];
        let psv_pos_per_copy = vec![vec![], vec![]];

        let r0 = bam_read("r0", "chr1", 0, copy0_seq, 60, 0.0);
        let r1 = bam_read("r1", "chr1", 5000, copy1_seq, 60, 0.0);
        let bam_reads = vec![r0, r1];
        let linear_copy_of = vec![Some(0usize), Some(1usize)];

        let out = apply_realign(
            &bam_reads,
            &copies,
            &copy_seqs,
            &psv_pos_per_copy,
            &linear_copy_of,
            &RealignParams::default(),
            0.003,
            1e-3,
        );

        assert!(out.corrected.is_empty(), "expected no corrections, got {:?}", out.corrected);
        assert!(out.novel_pools.is_empty(), "expected no novel pools, got {:?}", out.novel_pools);
        assert!(out.records.is_empty(), "expected no records at all, got {:?}", out.records);
    }

    #[test]
    fn path_obs_reads_bases_at_psv_columns() {
        // T-idx 5 = 'C', T-idx 10 = 'G' (0-based).
        let target: Vec<u8> = b"AAAAACAAAAGAAAAA".to_vec();
        assert_eq!(target[5], b'C');
        assert_eq!(target[10], b'G');

        // Exact full-length match -> both PSV columns are spanned and read verbatim.
        let query = target.clone();
        let map = align_traceback(&query, &target);
        let obs = path_obs_at(&map, &[5, 10], &query);
        assert_eq!(obs, vec![Some(b'C'), Some(b'G')]);

        // A short read that only covers target[..8] -- doesn't reach column 12 at all.
        let short_query: Vec<u8> = target[..8].to_vec();
        let map2 = align_traceback(&short_query, &target);
        let obs2 = path_obs_at(&map2, &[12], &short_query);
        assert_eq!(obs2, vec![None]);
    }
}
