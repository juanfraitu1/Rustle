//! Stage 3: per-block refinement -- strand, exact boundaries, identity.
//!
//! Stage 2 blocks are quantised to the window grid and their identity
//! estimate measures k-mer mass, not collinearity, so mosaic hubs inflate
//! it. This stage re-examines each block's two sequences directly:
//!
//! 1. **Anchors** -- shared canonical k-mers between A and B (canonical =
//!    strand-agnostic, so they serve both orientations), with a
//!    too-common filter like stage 1.
//! 2. **Anchor chaining in BOTH strand spaces** -- the same anchor list is
//!    walked once with B in forward coordinates and once with B
//!    reverse-complemented; each run must advance both coordinates by at
//!    most `max_gap` bp with per-step diagonal drift <= `diag_bin`. The
//!    winning space declares the strand -- collinearity, not voting.
//! 3. **Boundaries + identity** -- the best run's extent gives bp-exact
//!    coordinates on both sides; identity is anchor-interval coverage
//!    averaged over both spans, which only collinear matching sequence
//!    satisfies (mosaic blocks deflate to their real core).
//!
//! Cost per block is O(lenA + lenB) hashing -- no alignment is run, so the
//! whole stage stays laptop-cheap.

use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::HashMap;

use crate::sedef::anchors::BlockAnchors;
use crate::sedef::chain::Block;
use crate::sedef::kmer::for_each_kmer;

pub struct RefineParams {
    pub k: usize,
    /// Max bp gap between consecutive anchors inside a run.
    pub max_gap: u64,
    /// Per-step diagonal drift tolerance in bp (also the bin width).
    pub diag_bin: i64,
    /// Blocks refining to a shorter span than this are flagged.
    pub min_span: u64,
    /// Max anchor pairs emitted per k-mer hash; higher-multiplicity hashes
    /// are deterministically subsampled instead of truncated (truncation
    /// starves late positions in satellite arrays).
    pub max_pairs: u64,
    /// Scan v2 (SEDEF's final filter): a block's ok pairs together need at least this many uppercase (unmasked) bases on
    /// each side and this many aligned matches between uppercase bases; otherwise they are `repeat` records. 0 disables
    /// the filter.
    pub min_uppercase: u64,
    /// Scan v2: a chain whose alignment identity is below this (percent) is not a pair and shadows no other chain.
    pub min_identity: f64,
}

pub struct Refined {
    pub chr_a: String,
    pub ra0: u64,
    pub ra1: u64,
    pub chr_b: String,
    pub rb0: u64,
    pub rb1: u64,
    pub strand: char,
    pub anchors: usize,
    pub est_id: f64,
    /// Refined coordinates are trustworthy (enough anchors/span).
    pub ok: bool,
    /// A and B intervals are (nearly) the same locus -- tandem-overlap
    /// artifact rather than an interlocus duplication.
    pub self_pair: bool,
    /// Scan v2: an otherwise ok pair below `min_uppercase` (made of soft-masked repeats).
    pub repeat: bool,
    /// Stage-2 identity for comparison.
    pub stage2_id: f64,
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' | b'U' | b'u' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Union length of half-open intervals.
pub fn covered_len(intervals: &[(u64, u64)]) -> u64 {
    if intervals.is_empty() {
        return 0;
    }
    let mut v: Vec<(u64, u64)> = intervals.to_vec();
    v.sort_unstable();
    let mut total = 0u64;
    let (mut cs, mut ce) = v[0];
    for &(s, e) in &v[1..] {
        if s <= ce {
            ce = ce.max(e);
        } else {
            total += ce - cs;
            cs = s;
            ce = e;
        }
    }
    total + (ce - cs)
}

fn get_seq<'a>(seqs: &'a HashMap<String, Vec<u8>>, chr: &str) -> Option<&'a [u8]> {
    seqs.get(chr).map(|v| v.as_slice()).or_else(|| {
        // Block names are full FASTA headers; fall back to the first token.
        let first = chr.split_whitespace().next()?;
        seqs.get(first).map(|v| v.as_slice())
    })
}

/// One anchor run: consecutive collinear shared k-mers.
pub struct Run {
    pub pts: Vec<(u32, u32)>, // (posA, posB-in-current-strand-space)
    pub cov_a: u64,
    pub cov_b: u64,
    pub span_a: u64,
    pub span_b: u64,
}

/// Extract shared canonical k-mer anchor pairs between two sequences.
/// Emission per hash is capped by interleaved subsampling so late
/// positions in satellite arrays are not starved.
pub fn extract_anchors(sa: &[u8], sb: &[u8], k: usize, max_pairs: u64) -> Vec<(u32, u32)> {
    let mut pos: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
    for_each_kmer(sa, k, |pa, h| {
        pos.entry(h).or_default().push(pa as u32);
    });
    let mut cnt_b: FxHashMap<u64, u64> = FxHashMap::default();
    for_each_kmer(sb, k, |_, h| {
        if pos.contains_key(&h) {
            *cnt_b.entry(h).or_insert(0) += 1;
        }
    });
    let mut seen_b: FxHashMap<u64, u64> = FxHashMap::default();
    let mut anchors: Vec<(u32, u32)> = Vec::new();
    for_each_kmer(sb, k, |pb, h| {
        let Some(v) = pos.get(&h) else { return };
        let i = *seen_b.entry(h).or_insert(0);
        seen_b.insert(h, i + 1);
        let prod = v.len() as u64 * cnt_b[&h];
        let stride = prod.div_ceil(max_pairs).max(1);
        for (j, &pa) in v.iter().enumerate() {
            if (i * v.len() as u64 + j as u64) % stride == 0 {
                anchors.push((pa, pb as u32));
            }
        }
    });
    anchors
}

/// Union length of k-mer intervals given run points in ascending `pa`
/// order (runs are pushed one point per A-position, sorted).
fn cov_pa(pts: &[(u32, i64)], k: usize) -> u64 {
    let mut total = 0u64;
    let mut cs = 0u64;
    let mut ce = 0u64;
    let mut first = true;
    for &(pa, _) in pts {
        let (s, e) = (pa as u64, pa as u64 + k as u64);
        if first {
            cs = s;
            ce = e;
            first = false;
        } else if s <= ce {
            ce = ce.max(e);
        } else {
            total += ce - cs;
            cs = s;
            ce = e;
        }
    }
    if !first {
        total += ce - cs;
    }
    total
}

/// Union length of intervals [s, s+k) from arbitrary start positions.
fn cov_starts(starts: &mut Vec<u64>, k: u64) -> u64 {
    starts.sort_unstable();
    let mut total = 0u64;
    let mut cs = 0u64;
    let mut ce = 0u64;
    let mut first = true;
    for &s in starts.iter() {
        let e = s + k;
        if first {
            cs = s;
            ce = e;
            first = false;
        } else if s <= ce {
            ce = ce.max(e);
        } else {
            total += ce - cs;
            cs = s;
            ce = e;
        }
    }
    if !first {
        total += ce - cs;
    }
    total
}

fn flush_run(
    run: &mut Vec<(u32, i64)>,
    best_pts: &mut Vec<(u32, u32)>,
    best_cov: &mut u64,
    k: usize,
) {
    if run.len() < 3 {
        run.clear();
        return;
    }
    let cov = cov_pa(run, k);
    if cov > *best_cov {
        *best_cov = cov;
        *best_pts = run.iter().map(|&(pa, pb)| (pa, pb as u32)).collect();
    }
    run.clear();
}

/// The run of `runs` with the most A coverage (the first on ties): scan v1's refined core, and decompose's peeling step.
pub fn best_run(
    anchors: &[(u32, u32)],
    k: usize,
    max_gap: u64,
    diag_bin: i64,
    pb_map: impl Fn(u32) -> i64,
) -> Option<Run> {
    let mut best_pts: Vec<(u32, u32)> = Vec::new();
    let mut best_cov = 0u64;
    for mut run in runs(anchors, k, max_gap, diag_bin, pb_map) {
        flush_run(&mut run, &mut best_pts, &mut best_cov, k);
    }
    if best_pts.is_empty() {
        return None;
    }
    let pa_min = best_pts.first().unwrap().0 as u64;
    let pa_max = best_pts.last().unwrap().0 as u64 + k as u64;
    let mut scratch: Vec<u64> = best_pts.iter().map(|&(_, pb)| pb as u64).collect();
    let pb_min = *scratch.iter().min().unwrap();
    let pb_max = scratch.iter().max().unwrap() + k as u64;
    let cov_a = cov_pa(
        &best_pts
            .iter()
            .map(|&(pa, pb)| (pa, pb as i64))
            .collect::<Vec<_>>(),
        k,
    );
    let cov_b = cov_starts(&mut scratch, k as u64);
    Some(Run {
        pts: best_pts,
        cov_a,
        cov_b,
        span_a: pa_max - pa_min,
        span_b: pb_max - pb_min,
    })
}

/// Every run of the walker below with at least 3 points, in walk order; each run is (posA, posB in the chaining strand
/// space) in ascending posA.
///
/// The walk: greedy chaining of anchors within diagonal bins, with diagonal locking: candidates are grouped per
/// A-position and the walk always continues with the candidate closest to the run's current diagonal. This ignores
/// combinatorial junk pairs (multiple B positions for the same k-mer in repeat hubs) instead of letting them break or
/// pollute the run. `pb_map` converts anchor B positions into the chaining strand space.
///
/// Implementation: one sort by (bin, pa, pb) followed by a linear scan -- no hash maps, no per-bin vectors, no repeated
/// sorting.
pub fn runs(
    anchors: &[(u32, u32)],
    _k: usize,
    max_gap: u64,
    diag_bin: i64,
    pb_map: impl Fn(u32) -> i64,
) -> Vec<Vec<(u32, i64)>> {
    let mut out: Vec<Vec<(u32, i64)>> = Vec::new();
    if anchors.len() < 3 {
        return out;
    }
    let mut cells: Vec<(i64, u32, i64)> = Vec::with_capacity(anchors.len());
    for &(pa, pb) in anchors {
        let pb_m = pb_map(pb);
        let d = pb_m - pa as i64;
        cells.push((d.div_euclid(diag_bin), pa, pb_m));
    }
    cells.sort_unstable();

    let mut run: Vec<(u32, i64)> = Vec::new();
    let mut diag: i64 = 0; // current diagonal pb_m - pa
    let mut end_run = |run: &mut Vec<(u32, i64)>| {
        if run.len() >= 3 {
            out.push(std::mem::take(run));
        } else {
            run.clear();
        }
    };

    let mut i = 0;
    let n = cells.len();
    while i < n {
        let bin = cells[i].0;
        let mut j = i + 1;
        while j < n && cells[j].0 == bin {
            j += 1;
        }
        // Process [i, j): groups of equal pa (candidates sorted by pb_m).
        let mut g = i;
        while g < j {
            let pa = cells[g].1;
            let mut h = g + 1;
            while h < j && cells[h].1 == pa {
                h += 1;
            }
            let cands = &cells[g..h];
            let pick = if run.is_empty() {
                // Start a run on this group's median-ish candidate.
                cands[cands.len() / 2].2
            } else {
                // Continue with the candidate closest to the current
                // diagonal; break the run if none fits.
                let pa0 = run[run.len() - 1].0;
                let gap_ok = (pa - pa0) as u64 <= max_gap;
                let best = cands
                    .iter()
                    .min_by_key(|c| (c.2 - pa as i64 - diag).abs())
                    .map(|c| c.2)
                    .unwrap();
                let drift = (best - pa as i64 - diag).abs();
                if gap_ok && drift <= diag_bin {
                    best
                } else {
                    end_run(&mut run);
                    cands[cands.len() / 2].2
                }
            };
            diag = pick - pa as i64;
            run.push((pa, pick));
            g = h;
        }
        end_run(&mut run);
        i = j;
    }
    out
}

/// Scan v2: join runs of one strand space into chains. Runs are taken in ascending start; a run extends the chain whose
/// last point it follows with a step of at most `max_gap` on both axes (the walker's own gap rule, now allowing the
/// diagonal to move: an indel or a diverged stretch), closest diagonal first. Points of a run that do not advance past
/// the chain's last point on both axes are dropped (a run split at a diagonal-bin edge overlaps its neighbour).
pub fn chain_runs(mut runs: Vec<Vec<(u32, i64)>>, max_gap: u64) -> Vec<Vec<(u32, i64)>> {
    runs.sort_by_key(|r| (r[0].0, r[0].1));
    let mut chains: Vec<Vec<(u32, i64)>> = Vec::new();
    let mut active: Vec<usize> = Vec::new(); // chains whose last point may still be followed
    for run in runs {
        let start = run[0].0;
        active.retain(|&c| chains[c].last().unwrap().0 as u64 + max_gap >= start as u64);
        let mut best: Option<(i64, usize, usize)> = None; // (diagonal change, chain, first kept point of the run)
        for &c in &active {
            let (lpa, lpb) = *chains[c].last().unwrap();
            let Some(f) = run.iter().position(|&(pa, pb)| pa > lpa && pb > lpb) else { continue };
            let (pa, pb) = run[f];
            if (pa - lpa) as u64 > max_gap || (pb - lpb) as u64 > max_gap {
                continue;
            }
            let dd = ((pb - pa as i64) - (lpb - lpa as i64)).abs();
            if best.is_none_or(|b| (dd, c) < (b.0, b.1)) {
                best = Some((dd, c, f));
            }
        }
        match best {
            Some((_, c, f)) => chains[c].extend_from_slice(&run[f..]),
            None => {
                chains.push(run);
                active.push(chains.len() - 1);
            }
        }
    }
    chains
}

/// BISER's refine (`align.codon`, `refine`): a chain scores MATCH per matched base and -MISMATCH per other spanned base;
/// joining chain p to a later chain c costs GAPOPEN + MISMATCH x the shorter gap + GAP x (longer - shorter), for gaps
/// below MAX_GAP; at most MAX_ITER predecessors are tried per chain. Values are BISER's.
const BISER_MATCH: f64 = 10.0;
const BISER_MISMATCH: f64 = 1.0;
const BISER_GAP: f64 = 0.5;
const BISER_GAPOPEN: f64 = 100.0;
const BISER_MAX_GAP: i64 = 10 * 1024;
const BISER_MAX_ITER: usize = 500;

/// Scan v2: BISER's second-level chaining over the chains of one strand space. Chains sorted by end; each takes the
/// best-scoring predecessor (ending before it on both axes, the longer gap < MAX_GAP, the shorter <= `max_gap` so the
/// stretch both copies have is aligned, not assumed: BISER aligns it afterwards); paths are read from the best scores
/// down, each chain used once, and concatenated into one chain.
pub fn join_chains(mut chains: Vec<Vec<(u32, i64)>>, k: usize, max_gap: u64) -> Vec<Vec<(u32, i64)>> {
    if chains.len() < 2 {
        return chains;
    }
    let end = |c: &Vec<(u32, i64)>| {
        let l = c[c.len() - 1];
        (l.0 as i64 + k as i64, l.1 + k as i64)
    };
    chains.sort_by_key(|c| (end(c), c[0].0, c[0].1));
    let score: Vec<f64> = chains
        .iter()
        .map(|c| {
            let cov = cov_pa(c, k) as f64;
            let span = (end(c).0 - c[0].0 as i64) as f64;
            BISER_MATCH * cov - BISER_MISMATCH * (span - cov)
        })
        .collect();
    let n = chains.len();
    let mut dp = score.clone();
    let mut prev = vec![usize::MAX; n];
    for c in 0..n {
        let (cs_a, cs_b) = (chains[c][0].0 as i64, chains[c][0].1);
        let mut tried = 0;
        for p in (0..c).rev() {
            let (pe_a, pe_b) = end(&chains[p]);
            if pe_a + BISER_MAX_GAP <= cs_a {
                break; // sorted by end: every earlier chain ends even further away
            }
            // every examined predecessor counts (BISER counts only the feasible ones; a tandem block has thousands of
            // overlapping chains)
            tried += 1;
            if tried > BISER_MAX_ITER {
                break;
            }
            if pe_a > cs_a || pe_b > cs_b {
                continue;
            }
            let (ga, gb) = (cs_a - pe_a, cs_b - pe_b);
            let (mi, ma) = (ga.min(gb), ga.max(gb));
            if ma >= BISER_MAX_GAP || mi > max_gap as i64 {
                continue;
            }
            let sc = dp[p] + score[c] - BISER_MISMATCH * mi as f64 - BISER_GAPOPEN - BISER_GAP * (ma - mi) as f64;
            if sc >= dp[c] {
                dp[c] = sc;
                prev[c] = p;
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| dp[y].total_cmp(&dp[x]).then(x.cmp(&y)));
    let mut used = vec![false; n];
    let mut out = Vec::new();
    for top in order {
        if used[top] {
            continue;
        }
        let mut path = Vec::new();
        let mut i = top;
        while i != usize::MAX && !used[i] {
            used[i] = true;
            path.push(i);
            i = prev[i];
        }
        path.reverse();
        let mut pts = Vec::new();
        for i in path {
            pts.extend_from_slice(&chains[i]);
        }
        out.push(pts);
    }
    out
}

/// Unit-cost global alignment of two short sequences (case-insensitive): (matches, mismatches, matches of two uppercase
/// bases) on one optimal path, preferring the diagonal on ties.
pub fn align_counts(a: &[u8], b: &[u8]) -> (u64, u64, u64) {
    let (n, m) = (a.len(), b.len());
    if n == 0 || m == 0 {
        return (0, 0, 0);
    }
    let w = m + 1;
    let mut d = vec![0u32; (n + 1) * w];
    for j in 0..=m {
        d[j] = j as u32;
    }
    for i in 1..=n {
        d[i * w] = i as u32;
        let ai = a[i - 1].to_ascii_uppercase();
        for j in 1..=m {
            let sub = d[(i - 1) * w + j - 1] + (ai != b[j - 1].to_ascii_uppercase()) as u32;
            let del = d[(i - 1) * w + j] + 1;
            let ins = d[i * w + j - 1] + 1;
            d[i * w + j] = sub.min(del).min(ins);
        }
    }
    let (mut i, mut j) = (n, m);
    let (mut mat, mut mis, mut up) = (0u64, 0u64, 0u64);
    while i > 0 && j > 0 {
        let same = a[i - 1].to_ascii_uppercase() == b[j - 1].to_ascii_uppercase();
        if d[i * w + j] == d[(i - 1) * w + j - 1] + (!same) as u32 {
            if same {
                mat += 1;
                up += (a[i - 1].is_ascii_uppercase() && b[j - 1].is_ascii_uppercase()) as u64;
            } else {
                mis += 1;
            }
            i -= 1;
            j -= 1;
        } else if d[i * w + j] == d[(i - 1) * w + j] + 1 {
            i -= 1;
        } else {
            j -= 1;
        }
    }
    (mat, mis, up)
}

/// Scan v2 identity of a chain: anchors are k-mer matches, the stretches between consecutive anchors are aligned;
/// (100 x matches / (matches + mismatches), matches of two uppercase bases), indels left out of the identity (SEDEF's
/// fracMatch). `b` is in the chain's strand space.
/// A gap is aligned only up to `max_gap`^2 DP cells, the largest a run's own steps make; a longer gap of a BISER join
/// counts as an indel plus its shorter side as mismatches.
pub fn chain_identity(pts: &[(u32, i64)], a: &[u8], b: &[u8], k: usize, max_gap: u64) -> (f64, u64) {
    let both_up = |i: usize, j: usize| a[i].is_ascii_uppercase() && b[j].is_ascii_uppercase();
    // uppercase matches of `n` matched bases starting at a[i], b[j]
    let up_run = |i: usize, j: usize, n: usize| (0..n).filter(|&t| both_up(i + t, j + t)).count() as u64;
    let (p0a, p0b) = (pts[0].0 as usize, pts[0].1 as usize);
    let mut mat = k as u64;
    let mut up = up_run(p0a, p0b, k);
    let mut mis = 0u64;
    for w in pts.windows(2) {
        let ((pa, pb), (qa, qb)) = (w[0], w[1]);
        let (da, db) = ((qa - pa) as i64, qb - pb);
        if da == db && da as usize <= k {
            mat += da as u64;
            up += up_run(pa as usize + k, pb as usize + k, da as usize);
            continue;
        }
        let (ga, gb) = (pa as usize + k, (pb + k as i64) as usize);
        let cells = (qa as usize).saturating_sub(ga) * (qb.max(0) as usize).saturating_sub(gb);
        if qa as usize >= ga && qb >= 0 && qb as usize >= gb && (qb as usize) <= b.len() && cells as u64 > max_gap * max_gap {
            // a long gap between joined chains, not aligned (BISER aligns it): its longer side's excess is an indel, and
            // the shorter side's bases are unverified, so they count as mismatches (identity stays a lower bound)
            mis += ((qa as usize - ga).min(qb as usize - gb)) as u64;
            mat += k as u64;
            up += up_run(qa as usize, qb as usize, k);
        } else if qa as usize >= ga && qb >= 0 && qb as usize >= gb && (qb as usize) <= b.len() {
            let (m2, x2, u2) = align_counts(&a[ga..qa as usize], &b[gb..qb as usize]);
            mat += m2 + k as u64;
            mis += x2;
            up += u2 + up_run(qa as usize, qb as usize, k);
        } else {
            // the two k-mers overlap off one diagonal (a short indel inside them): count the new A bases both advance
            mat += da.min(db).max(0) as u64;
        }
    }
    if mat + mis == 0 {
        return (0.0, 0);
    }
    (100.0 * mat as f64 / (mat + mis) as f64, up)
}

/// Scan v2: every chain of the block (`chain_runs` in both strand spaces) whose two spans reach `min_span` is a pair,
/// in ascending A start; a block with none gives its best chain (or, without any run, the block itself) as one weak
/// record, as scan v1 did.
fn refine_block(b: &Block, seqs: &HashMap<String, Vec<u8>>, p: &RefineParams) -> Vec<(Refined, BlockAnchors)> {
    let (Some(sa_full), Some(sb_full)) = (get_seq(seqs, &b.chr_a), get_seq(seqs, &b.chr_b)) else {
        return Vec::new();
    };
    let (a0, a1) = (
        (b.a0 as usize).min(sa_full.len()),
        (b.a1 as usize).min(sa_full.len()),
    );
    let (b0, b1) = (
        (b.b0 as usize).min(sb_full.len()),
        (b.b1 as usize).min(sb_full.len()),
    );
    if a1 <= a0 || b1 <= b0 {
        return Vec::new();
    }
    let sa = &sa_full[a0..a1];
    let sb = &sb_full[b0..b1];
    let k = p.k;
    let sb_len = sb.len() as i64;

    // ---- anchors: shared canonical k-mers (strand-agnostic) -------------
    let anchors = extract_anchors(sa, sb, k, p.max_pairs);

    // ---- chains in both strand spaces -----------------------------------
    let mut chains: Vec<(char, Vec<(u32, i64)>)> = Vec::new();
    for strand in ['+', '-'] {
        let rs = if strand == '+' {
            runs(&anchors, k, p.max_gap, p.diag_bin, |pb| pb as i64)
        } else {
            runs(&anchors, k, p.max_gap, p.diag_bin, |pb| sb_len - pb as i64 - k as i64)
        };
        chains.extend(join_chains(chain_runs(rs, p.max_gap), k, p.max_gap).into_iter().map(|c| (strand, c)));
    }
    if chains.is_empty() {
        let rec = Refined {
            chr_a: b.chr_a.clone(),
            ra0: b.a0,
            ra1: b.a1,
            chr_b: b.chr_b.clone(),
            rb0: b.b0,
            rb1: b.b1,
            strand: '+',
            anchors: 0,
            est_id: b.est_id,
            ok: false,
            self_pair: false,
            repeat: false,
            stage2_id: b.est_id,
        };
        // No refined core exists; decompose only consumes ok records,
        // so store an empty anchor list.
        let ba = BlockAnchors {
            chr_a: b.chr_a.clone(),
            a0: b.a0,
            a1: b.a1,
            chr_b: b.chr_b.clone(),
            b0: b.b0,
            b1: b.b1,
            status: 0,
            anchors: Vec::new(),
        };
        return vec![(rec, ba)];
    }
    let sb_rc = if chains.iter().any(|(s, _)| *s == '-') { revcomp(sb) } else { Vec::new() };

    // ---- refined coordinates per chain ----------------------------------
    // (strand, points, ra0, ra1, rb0, rb1, A coverage) for every chain
    let mut cands: Vec<(char, &Vec<(u32, i64)>, u64, u64, u64, u64, u64)> = Vec::new();
    for (strand, pts) in &chains {
        let pa_min = pts.first().unwrap().0 as u64;
        let pa_max = pts.last().unwrap().0 as u64 + k as u64;
        let (ra0, ra1) = (a0 as u64 + pa_min, a0 as u64 + pa_max);
        // chain pb coords are in the strand space: original sb (+) or revcomp(sb) (-), where position q corresponds to
        // the sb interval [L-q-k, L-q).
        let pb_min = pts.iter().map(|&(_, pb)| pb as u64).min().unwrap();
        let pb_max = pts.iter().map(|&(_, pb)| pb as u64).max().unwrap() + k as u64;
        let (rb0, rb1) = if *strand == '-' {
            (b0 as u64 + sb_len as u64 - pb_max, b0 as u64 + sb_len as u64 - pb_min)
        } else {
            (b0 as u64 + pb_min, b0 as u64 + pb_max)
        };
        cands.push((*strand, pts, ra0, ra1, rb0, rb1, cov_pa(pts, k)));
    }
    // Most A coverage first; a chain is kept when its A x B rectangle meets no kept chain's (shifted diagonals of a
    // tandem array, or the other strand of a palindrome, chain the same two stretches again). Chains below min_span are
    // only the fallback weak record.
    cands.sort_by_key(|c| (std::cmp::Reverse(c.6), c.2, c.4, c.0));
    let is_ok = |c: &(char, &Vec<(u32, i64)>, u64, u64, u64, u64, u64)| {
        c.3 - c.2 >= p.min_span && c.5 - c.4 >= p.min_span && c.1.len() >= 3
    };
    let ident = |i: usize| {
        let (strand, pts, ..) = cands[i];
        chain_identity(pts, sa, if strand == '-' { &sb_rc } else { sb }, k, p.max_gap)
    };
    // (candidate, identity, uppercase matches) of the kept chains
    let mut chosen: Vec<(usize, f64, u64)> = Vec::new();
    for (i, c) in cands.iter().enumerate() {
        if !is_ok(c) {
            continue;
        }
        let meets = chosen.iter().any(|&(j, ..)| {
            let d = &cands[j];
            c.2 < d.3 && d.2 < c.3 && c.4 < d.5 && d.4 < c.5
        });
        if meets {
            continue;
        }
        let (id, up) = ident(i);
        if id >= p.min_identity {
            chosen.push((i, id, up));
        }
    }
    let mut fallback = false;
    if chosen.is_empty() {
        let (id, up) = ident(0);
        chosen.push((0, id, up)); // the best chain, as a weak record
        fallback = true;
    }
    let mut up_stats: Vec<(u64, u64, u64)> = Vec::new(); // per kept record: unmasked matches, unmasked bases A, B
    let mut kept: Vec<Refined> = chosen
        .into_iter()
        .map(|(i, est_id, up_matches)| {
            let (strand, pts, ra0, ra1, rb0, rb1, _) = cands[i];
            let ok = is_ok(&cands[i]) && !fallback;
            let upper = |seq: &[u8], x0: u64, x1: u64| seq[x0 as usize..x1 as usize].iter().filter(|c| c.is_ascii_uppercase()).count() as u64;
            up_stats.push((up_matches, upper(sa_full, ra0, ra1), upper(sb_full, rb0, rb1)));
            // "Self" only when the two sides land on (nearly) the SAME locus --
            // a refinement collapse artifact. Shifted copies (offset >= 1 kb)
            // are genuine tandem duplications.
            let self_pair = b.chr_a == b.chr_b && ra0.abs_diff(rb0) < 1000 && ra1.abs_diff(rb1) < 1000;
            Refined {
                chr_a: b.chr_a.clone(),
                ra0,
                ra1,
                chr_b: b.chr_b.clone(),
                rb0,
                rb1,
                strand,
                anchors: pts.len(),
                est_id,
                ok,
                self_pair,
                repeat: false,
                stage2_id: b.est_id,
            }
        })
        .collect();
    // SEDEF's final filter, over the block's ok pairs together (they are pieces of one duplication): unmasked matches and
    // unmasked bases per side, summed
    if p.min_uppercase > 0 {
        let (mut um, mut ua, mut ub) = (0u64, 0u64, 0u64);
        for (r, &(m, a, b)) in kept.iter().zip(&up_stats) {
            if r.ok && !r.self_pair {
                (um, ua, ub) = (um + m, ua + a, ub + b);
            }
        }
        if um < p.min_uppercase || ua < p.min_uppercase || ub < p.min_uppercase {
            for r in kept.iter_mut().filter(|r| r.ok && !r.self_pair) {
                r.repeat = true;
            }
        }
    }
    kept.sort_by_key(|r| (r.ra0, r.ra1, r.rb0, r.rb1, r.strand));

    // Re-extract anchors on each refined core (not the wider stage-2
    // block) so the stored set is exactly what `decompose` would extract
    // from the refined intervals. Filtering the block-frame anchors is
    // NOT equivalent: subsampling under the max_pairs cap depends on
    // per-hash multiplicities, which change with the interval.
    kept.into_iter()
        .map(|r| {
            let status: u8 = if !r.ok {
                0
            } else if r.self_pair {
                2
            } else if r.repeat {
                3
            } else {
                1
            };
            let ra1_c = (r.ra1 as usize).min(sa_full.len());
            let rb1_c = (r.rb1 as usize).min(sb_full.len());
            let anchors = extract_anchors(&sa_full[r.ra0 as usize..ra1_c], &sb_full[r.rb0 as usize..rb1_c], k, p.max_pairs);
            let ba = BlockAnchors {
                chr_a: r.chr_a.clone(),
                a0: r.ra0,
                a1: r.ra1,
                chr_b: r.chr_b.clone(),
                b0: r.rb0,
                b1: r.rb1,
                status,
                anchors,
            };
            (r, ba)
        })
        .collect()
}

pub fn refine_blocks(
    blocks: &[Block],
    seqs: &HashMap<String, Vec<u8>>,
    p: &RefineParams,
) -> Vec<(Refined, BlockAnchors)> {
    blocks
        .par_iter()
        .flat_map_iter(|b| refine_block(b, seqs, p))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // Deterministic pseudo-DNA so tests are reproducible.
    fn synth(seed: u64, n: usize) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 62) as usize]
            })
            .collect()
    }

    /// Mutate `rate` percent of bases (deterministic, no indels).
    fn mutate(seq: &[u8], rate: u32, seed: u64) -> Vec<u8> {
        let mut x = seed;
        seq.iter()
            .map(|&b| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                if (x >> 32) % 100 < rate as u64 {
                    b"ACGT"[((x >> 40) % 4) as usize]
                } else {
                    b
                }
            })
            .collect()
    }

    fn params() -> RefineParams {
        RefineParams {
            k: 16,
            max_gap: 200,
            diag_bin: 50,
            min_span: 500,
            max_pairs: 1024,
            min_uppercase: 0,
            min_identity: 90.0,
        }
    }

    fn block(a0: u64, a1: u64, b0: u64, b1: u64) -> Block {
        Block {
            chr_a: "c".into(),
            a0,
            a1,
            chr_b: "c".into(),
            b0,
            b1,
            cells: 1,
            score: 100,
            est_id: 99.0,
            upper_seed: true,
        }
    }

    #[test]
    fn refines_direct_duplication() {
        let core = mutate(&synth(7, 20_000), 5, 99); // ~95% identical copy
        let mut seq = synth(1, 60_000);
        seq[10_000..30_000].copy_from_slice(&core);
        seq[40_000..60_000].copy_from_slice(&core);
        let mut seqs = HashMap::new();
        seqs.insert("c".to_string(), seq);

        let r = refine_block(&block(9_000, 31_000, 39_000, 61_000), &seqs, &params()).swap_remove(0).0;
        assert_eq!(r.strand, '+');
        assert!(r.ok, "expected ok refinement");
        assert!(r.ra0.abs_diff(10_000) <= 50, "ra0={}", r.ra0);
        assert!(r.ra1.abs_diff(30_000) <= 50, "ra1={}", r.ra1);
        assert!(r.rb0.abs_diff(40_000) <= 50, "rb0={}", r.rb0);
        assert!(r.rb1.abs_diff(60_000) <= 50, "rb1={}", r.rb1);
        assert!(r.est_id > 88.0 && r.est_id <= 100.0, "est_id={}", r.est_id);
    }

    #[test]
    fn detects_inverted_duplication() {
        let core = mutate(&synth(7, 20_000), 5, 99);
        let mut seq = synth(1, 60_000);
        seq[10_000..30_000].copy_from_slice(&core);
        let rc = revcomp(&core);
        seq[40_000..60_000].copy_from_slice(&rc);
        let mut seqs = HashMap::new();
        seqs.insert("c".to_string(), seq);

        let r = refine_block(&block(9_000, 31_000, 39_000, 61_000), &seqs, &params()).swap_remove(0).0;
        assert_eq!(r.strand, '-');
        assert!(r.ok);
        assert!(r.rb0.abs_diff(40_000) <= 50, "rb0={}", r.rb0);
        assert!(r.rb1.abs_diff(60_000) <= 50, "rb1={}", r.rb1);
        assert!(r.est_id > 88.0, "est_id={}", r.est_id);
    }

    #[test]
    fn mosaic_hub_deflates_identity() {
        // A = X ++ Y, B = X ++ Z  (shares only half with A)
        let x = mutate(&synth(11, 10_000), 5, 5);
        let y = synth(13, 10_000);
        let z = synth(17, 10_000);
        let mut seq = synth(1, 60_000);
        seq[10_000..20_000].copy_from_slice(&x);
        seq[20_000..30_000].copy_from_slice(&y);
        seq[40_000..50_000].copy_from_slice(&x);
        seq[50_000..60_000].copy_from_slice(&z);
        let mut seqs = HashMap::new();
        seqs.insert("c".to_string(), seq);

        let r = refine_block(&block(9_000, 31_000, 39_000, 61_000), &seqs, &params()).swap_remove(0).0;
        // stage-2-style estimate would be ~99 (k-mer mass of shared X);
        // refinement must collapse to the ~10 kb real shared core.
        assert!(r.ra1 - r.ra0 < 12_000, "span={}", r.ra1 - r.ra0);
        assert!(r.est_id > 90.0, "est_id={}", r.est_id);
    }

    #[test]
    fn stored_anchors_match_fresh_extraction_on_refined_coords() {
        // The invariant behind `decompose --anchors`: the rebased anchor
        // list stored by refine must equal freshly extracted anchors on
        // the refined intervals -- for BOTH strands.
        let core = mutate(&synth(7, 20_000), 5, 99);
        let rc = revcomp(&core);

        for (want_strand, b_copy) in [('+', &core[..]), ('-', &rc[..])] {
            let mut a_seq = synth(1, 40_000);
            a_seq[10_000..30_000].copy_from_slice(&core);
            let mut b_seq = synth(2, 80_000);
            b_seq[40_000..60_000].copy_from_slice(b_copy);
            let mut seqs = HashMap::new();
            seqs.insert("c1".to_string(), a_seq);
            seqs.insert("c2".to_string(), b_seq);

            let blk = Block {
                chr_a: "c1".into(),
                a0: 9_000,
                a1: 31_000,
                chr_b: "c2".into(),
                b0: 39_000,
                b1: 61_000,
                cells: 1,
                score: 100,
                est_id: 99.0,
                upper_seed: true,
            };
            let (rec, ba) = refine_block(&blk, &seqs, &params()).swap_remove(0);
            assert_eq!(rec.strand, want_strand);
            assert_eq!((ba.a0, ba.a1), (rec.ra0, rec.ra1));
            assert_eq!((ba.b0, ba.b1), (rec.rb0, rec.rb1));

            let sa = &seqs["c1"][ba.a0 as usize..ba.a1 as usize];
            let sb = &seqs["c2"][ba.b0 as usize..ba.b1 as usize];
            let mut got = ba.anchors.clone();
            got.sort_unstable();
            let mut want = extract_anchors(sa, sb, 16, 1024);
            want.sort_unstable();
            assert_eq!(got, want, "strand {}", want_strand);
        }
    }

    fn ok_records(b: &Block, seqs: &HashMap<String, Vec<u8>>) -> Vec<Refined> {
        refine_block(b, seqs, &params()).into_iter().map(|(r, _)| r).filter(|r| r.ok).collect()
    }

    #[test]
    fn an_indel_does_not_cut_the_pair() {
        // ~95% copy with a 150 bp insertion in the B copy halfway: one pair over the whole copy.
        let core = mutate(&synth(7, 20_000), 5, 99);
        let mut b_copy = core[..10_000].to_vec();
        b_copy.extend_from_slice(&synth(21, 150));
        b_copy.extend_from_slice(&core[10_000..]);
        let mut seq = synth(1, 70_000);
        seq[10_000..30_000].copy_from_slice(&core);
        seq[40_000..60_150].copy_from_slice(&b_copy);
        let mut seqs = HashMap::new();
        seqs.insert("c".to_string(), seq);
        let recs = ok_records(&block(9_000, 31_000, 39_000, 61_000), &seqs);
        assert!(recs.iter().any(|r| r.ra1 - r.ra0 >= 19_500 && r.rb1 - r.rb0 >= 19_600), "spans {:?}",
            recs.iter().map(|r| (r.ra0, r.ra1, r.rb0, r.rb1)).collect::<Vec<_>>());
    }

    #[test]
    fn a_long_insertion_inside_a_copy_is_bridged() {
        // a 2 kb insertion in the B copy: BISER's gap-scored join keeps one pair over both halves
        let core = mutate(&synth(7, 10_000), 3, 9);
        let mut b_copy = core[..5_000].to_vec();
        b_copy.extend_from_slice(&synth(23, 2_000));
        b_copy.extend_from_slice(&core[5_000..]);
        let mut seq = synth(1, 60_000);
        seq[10_000..20_000].copy_from_slice(&core);
        seq[30_000..42_000].copy_from_slice(&b_copy);
        let mut seqs = HashMap::new();
        seqs.insert("c".to_string(), seq);
        let recs = ok_records(&block(9_000, 21_000, 29_000, 43_000), &seqs);
        assert!(recs.iter().any(|r| r.ra1 - r.ra0 >= 9_900 && r.rb1 - r.rb0 >= 11_900),
            "{:?}", recs.iter().map(|r| (r.ra0, r.ra1, r.rb0, r.rb1)).collect::<Vec<_>>());
    }

    #[test]
    fn two_duplications_in_one_block_are_two_pairs() {
        // X and Y are both duplicated, with 5 kb between them on A and 7 kb on B: two pairs, not the better one only.
        let x = mutate(&synth(31, 5_000), 5, 1);
        let y = mutate(&synth(37, 5_000), 5, 2);
        let mut seq = synth(1, 70_000);
        seq[10_000..15_000].copy_from_slice(&x);
        seq[20_000..25_000].copy_from_slice(&y);
        seq[40_000..45_000].copy_from_slice(&x);
        seq[52_000..57_000].copy_from_slice(&y);
        let mut seqs = HashMap::new();
        seqs.insert("c".to_string(), seq);
        let recs = ok_records(&block(9_000, 26_000, 39_000, 58_000), &seqs);
        let near = |r: &Refined, a0: u64, b0: u64| r.ra0.abs_diff(a0) <= 50 && r.rb0.abs_diff(b0) <= 50;
        assert!(recs.iter().any(|r| near(r, 10_000, 40_000)), "X missing");
        assert!(recs.iter().any(|r| near(r, 20_000, 52_000)), "Y missing");
    }

    #[test]
    fn a_tandem_array_gives_disjoint_pairs_not_one_per_diagonal() {
        // 20 tandem copies of a 1 kb unit on A and on B: every shifted diagonal chains, but the kept pairs must not
        // overlap each other on both axes.
        let unit = synth(51, 1_000);
        let mut seq = synth(1, 80_000);
        for i in 0..20 {
            seq[5_000 + i * 1_000..6_000 + i * 1_000].copy_from_slice(&unit);
            seq[45_000 + i * 1_000..46_000 + i * 1_000].copy_from_slice(&mutate(&unit, 3, i as u64));
        }
        let mut seqs = HashMap::new();
        seqs.insert("c".to_string(), seq);
        let recs = ok_records(&block(4_000, 26_000, 44_000, 66_000), &seqs);
        assert!(!recs.is_empty());
        for (i, r) in recs.iter().enumerate() {
            for q in &recs[i + 1..] {
                let over = |a0: u64, a1: u64, b0: u64, b1: u64| a0 < b1 && b0 < a1;
                assert!(!(over(r.ra0, r.ra1, q.ra0, q.ra1) && over(r.rb0, r.rb1, q.rb0, q.rb1)), "overlapping pairs");
            }
        }
    }

    #[test]
    fn identity_is_the_alignment_identity() {
        // mutate(rate) substitutes rate% of bases by a random base, a quarter of them by the same base:
        // the true identity is 100 - 0.75 * rate.
        for rate in [4u32, 8, 10] {
            let core = mutate(&synth(41, 10_000), rate, 7);
            let mut seq = synth(1, 40_000);
            seq[5_000..15_000].copy_from_slice(&synth(41, 10_000));
            seq[25_000..35_000].copy_from_slice(&core);
            let mut seqs = HashMap::new();
            seqs.insert("c".to_string(), seq);
            let recs = ok_records(&block(4_000, 16_000, 24_000, 36_000), &seqs);
            let best = recs.iter().max_by_key(|r| r.ra1 - r.ra0).expect("a pair");
            let want = 100.0 - 0.75 * rate as f64;
            assert!((best.est_id - want).abs() <= 1.0, "rate {rate}: est_id {} want {want}", best.est_id);
        }
    }

    #[test]
    fn combinatorial_junk_doesnt_form_long_runs() {
        // A and B both contain 30 copies of a 300 bp repeat, plus a single
        // 2 kb shared core. Only the core should chain; the combinatorial
        // repeat pairs must not inflate coverage.
        let core = mutate(&synth(3, 2_000), 4, 42);
        let rep = synth(9, 300);
        let mut seq = synth(1, 40_000);
        for i in 0..30 {
            let s = 1_000 + i * 500;
            seq[s..s + 300].copy_from_slice(&rep);
        }
        seq[5_000..7_000].copy_from_slice(&core);
        for i in 0..30 {
            let s = 20_000 + i * 500;
            seq[s..s + 300].copy_from_slice(&rep);
        }
        seq[25_000..27_000].copy_from_slice(&core);
        let mut seqs = HashMap::new();
        seqs.insert("c".to_string(), seq);

        let r = refine_block(&block(0, 18_000, 19_000, 37_000), &seqs, &params()).swap_remove(0).0;
        assert_eq!(r.strand, '+');
        assert!(r.ra1 - r.ra0 < 4_000, "span={}", r.ra1 - r.ra0);
        assert!(r.est_id > 90.0, "est_id={}", r.est_id);
    }
}
