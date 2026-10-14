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

/// Greedy chaining of anchors within diagonal bins, with diagonal
/// locking: candidates are grouped per A-position and the walk always
/// continues with the candidate closest to the run's current diagonal.
/// This ignores combinatorial junk pairs (multiple B positions for the
/// same k-mer in repeat hubs) instead of letting them break or pollute
/// the run. `pb_map` converts anchor B positions into the chaining
/// strand space.
///
/// Implementation: one sort by (bin, pa, pb) followed by a linear scan --
/// no hash maps, no per-bin vectors, no repeated sorting.
pub fn best_run(
    anchors: &[(u32, u32)],
    k: usize,
    max_gap: u64,
    diag_bin: i64,
    pb_map: impl Fn(u32) -> i64,
) -> Option<Run> {
    if anchors.len() < 3 {
        return None;
    }
    let mut cells: Vec<(i64, u32, i64)> = Vec::with_capacity(anchors.len());
    for &(pa, pb) in anchors {
        let pb_m = pb_map(pb);
        let d = pb_m - pa as i64;
        cells.push((d.div_euclid(diag_bin), pa, pb_m));
    }
    cells.sort_unstable();

    let mut best_pts: Vec<(u32, u32)> = Vec::new();
    let mut best_cov = 0u64;
    let mut run: Vec<(u32, i64)> = Vec::new();
    let mut diag: i64 = 0; // current diagonal pb_m - pa

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
                    flush_run(&mut run, &mut best_pts, &mut best_cov, k);
                    cands[cands.len() / 2].2
                }
            };
            diag = pick - pa as i64;
            run.push((pa, pick));
            g = h;
        }
        flush_run(&mut run, &mut best_pts, &mut best_cov, k);
        i = j;
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

fn refine_block(
    b: &Block,
    seqs: &HashMap<String, Vec<u8>>,
    p: &RefineParams,
) -> Option<(Refined, BlockAnchors)> {
    let sa_full = get_seq(seqs, &b.chr_a)?;
    let sb_full = get_seq(seqs, &b.chr_b)?;
    let (a0, a1) = (
        (b.a0 as usize).min(sa_full.len()),
        (b.a1 as usize).min(sa_full.len()),
    );
    let (b0, b1) = (
        (b.b0 as usize).min(sb_full.len()),
        (b.b1 as usize).min(sb_full.len()),
    );
    if a1 <= a0 || b1 <= b0 {
        return None;
    }
    let sa = &sa_full[a0..a1];
    let sb = &sb_full[b0..b1];
    let k = p.k;
    let sb_len = sb.len() as i64;

    // ---- anchors: shared canonical k-mers (strand-agnostic) -------------
    let anchors = extract_anchors(sa, sb, k, p.max_pairs);

    // ---- chain in both strand spaces, keep the better run ---------------
    let fwd = best_run(&anchors, k, p.max_gap, p.diag_bin, |pb| pb as i64);
    let rev = best_run(&anchors, k, p.max_gap, p.diag_bin, |pb| sb_len - pb as i64 - k as i64);
    let (strand, run) = match (fwd, rev) {
        (Some(f), Some(r)) => {
            let score = |x: &Run| x.cov_a.min(x.cov_b);
            if score(&f) >= score(&r) {
                ('+', f)
            } else {
                ('-', r)
            }
        }
        (Some(f), None) => ('+', f),
        (None, Some(r)) => ('-', r),
        (None, None) => {
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
            return Some((rec, ba));
        }
    };

    // ---- refined coordinates --------------------------------------------
    let pa_min = run.pts.first().unwrap().0 as u64;
    let pa_max = run.pts.last().unwrap().0 as u64 + k as u64;
    let (ra0, ra1) = (a0 as u64 + pa_min, a0 as u64 + pa_max);

    // run pb coords are already in the winning strand space; the space is
    // either original sb (fwd) or revcomp(sb) (rev). For rev, position q
    // in revcomp(sb) corresponds to sb interval [L-q-k, L-q).
    let pb_min = run.pts.iter().map(|&(_, pb)| pb as u64).min().unwrap();
    let pb_max = run.pts.iter().map(|&(_, pb)| pb as u64).max().unwrap() + k as u64;
    let (rb0, rb1) = if strand == '-' {
        (
            b0 as u64 + sb_len as u64 - pb_max,
            b0 as u64 + sb_len as u64 - pb_min,
        )
    } else {
        (b0 as u64 + pb_min, b0 as u64 + pb_max)
    };

    let span_a = ra1 - ra0;
    let span_b = rb1 - rb0;
    let mean_span = (span_a + span_b) as f64 / 2.0;
    let est_id = if mean_span > 0.0 {
        (100.0 * (run.cov_a + run.cov_b) as f64 / (2.0 * mean_span)).min(100.0)
    } else {
        0.0
    };
    let n_anchors = run.pts.len();
    let ok = span_a >= p.min_span && span_b >= p.min_span && n_anchors >= 3;
    // "Self" only when the two sides land on (nearly) the SAME locus --
    // a refinement collapse artifact. Shifted copies (offset >= 1 kb)
    // are genuine tandem duplications.
    let self_pair = b.chr_a == b.chr_b
        && ra0.abs_diff(rb0) < 1000
        && ra1.abs_diff(rb1) < 1000;
    let status: u8 = if !ok {
        0
    } else if self_pair {
        2
    } else {
        1
    };

    // Re-extract anchors on the refined core (not the wider stage-2
    // block) so the stored set is exactly what `decompose` would extract
    // from the refined intervals. Filtering the block-frame anchors is
    // NOT equivalent: subsampling under the max_pairs cap depends on
    // per-hash multiplicities, which change with the interval.
    let ra1_c = (ra1 as usize).min(sa_full.len());
    let rb1_c = (rb1 as usize).min(sb_full.len());
    let kept = extract_anchors(&sa_full[ra0 as usize..ra1_c], &sb_full[rb0 as usize..rb1_c], k, p.max_pairs);
    let rec = Refined {
        chr_a: b.chr_a.clone(),
        ra0,
        ra1,
        chr_b: b.chr_b.clone(),
        rb0,
        rb1,
        strand,
        anchors: n_anchors,
        est_id,
        ok,
        self_pair,
        stage2_id: b.est_id,
    };
    let ba = BlockAnchors {
        chr_a: b.chr_a.clone(),
        a0: ra0,
        a1: ra1,
        chr_b: b.chr_b.clone(),
        b0: rb0,
        b1: rb1,
        status,
        anchors: kept,
    };
    Some((rec, ba))
}

pub fn refine_blocks(
    blocks: &[Block],
    seqs: &HashMap<String, Vec<u8>>,
    p: &RefineParams,
) -> Vec<(Refined, BlockAnchors)> {
    blocks
        .par_iter()
        .filter_map(|b| refine_block(b, seqs, p))
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

        let r = refine_block(&block(9_000, 31_000, 39_000, 61_000), &seqs, &params()).unwrap().0;
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

        let r = refine_block(&block(9_000, 31_000, 39_000, 61_000), &seqs, &params()).unwrap().0;
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

        let r = refine_block(&block(9_000, 31_000, 39_000, 61_000), &seqs, &params()).unwrap().0;
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
            };
            let (rec, ba) = refine_block(&blk, &seqs, &params()).unwrap();
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

        let r = refine_block(&block(0, 18_000, 19_000, 37_000), &seqs, &params()).unwrap().0;
        assert_eq!(r.strand, '+');
        assert!(r.ra1 - r.ra0 < 4_000, "span={}", r.ra1 - r.ra0);
        assert!(r.est_id > 90.0, "est_id={}", r.est_id);
    }
}
