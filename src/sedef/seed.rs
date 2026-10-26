//! Stage 1 of the SD-detection cascade: seeding + cheap filtering.
//!
//! The genome is cut into non-overlapping windows. Each window is
//! summarised by a "bottom-S" sketch: the S smallest canonical k-mer
//! hashes it contains. An inverted index maps each sketch value to the
//! windows that contain it; window pairs that co-occur on at least
//! `min_shared` sketch values become SD candidates. This mirrors the
//! seeding/filtering philosophy of SEDEF/BISER while staying simple
//! enough to audit.

use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::HashMap;

use crate::sedef::kmer::for_each_kmer;

pub struct SeedParams {
    pub k: usize,
    pub window: usize,
    /// Scan v2: windows start every `stride` bp (scan v1: `stride` = `window`, a non-overlapping grid; v2: half a
    /// window, so a copy is never more than a quarter window off a window of the other copy).
    pub stride: usize,
    pub sketch: usize,
    /// Windows with fewer distinct k-mers than this are treated as pure
    /// satellite repeats and skipped (SDFinder/BISER do analogous repeat
    /// masking upstream).
    pub min_kmers: usize,
    /// Sketch values occurring in more windows than this are discarded
    /// (they carry no pairing signal -- classic too-common filter).
    pub max_postings: usize,
    pub min_shared: usize,
    /// Scan v2 (SEDEF's seed rule, `--uppercase` = k): every k-mer counts toward the shared mass, and a candidate is
    /// marked `upper_seed` when a shared sketch k-mer is all uppercase (unmasked) in both windows; the chain stage keeps
    /// only blocks with such a seed, so common repeats seed nothing while a mostly-masked duplication is extended from
    /// its unmasked part.
    pub uppercase_seed: bool,
    /// Vendored change (2026-10-13): pairs are counted in this many shards (a pair belongs to the shard of its lower
    /// window index modulo `count_shards`), so at most one shard map per thread is alive. Every value gives the same
    /// candidates; 1 is the sandbox's single map.
    pub count_shards: usize,
}

pub struct SeqInfo {
    pub name: String,
    pub len: usize,
}

pub struct WindowSketch {
    pub seq: u32,
    pub start: u64,
    pub end: u64,
    /// Bottom-S distinct canonical k-mer hashes, sorted ascending.
    pub sketch: Vec<u64>,
    /// With `uppercase_seed`: per sketch hash, whether one of its k-mers in the window is all uppercase (else empty).
    pub upper: Vec<bool>,
}

pub struct Candidate {
    pub a: usize,
    pub b: usize,
    pub shared: u32,
    /// A shared sketch k-mer is all uppercase in both windows (always false without `uppercase_seed`).
    pub upper_seed: bool,
}

/// Cut all sequences into window sketches. Returned in genomic order.
pub fn build_windows(
    seqs: &[(String, Vec<u8>)],
    p: &SeedParams,
) -> (Vec<WindowSketch>, Vec<SeqInfo>) {
    let infos: Vec<SeqInfo> = seqs
        .iter()
        .map(|(name, s)| SeqInfo {
            name: name.clone(),
            len: s.len(),
        })
        .collect();

    let mut all: Vec<WindowSketch> = Vec::new();
    for (si, (_, seq)) in seqs.iter().enumerate() {
        let stride = p.stride.max(1);
        let n_win = if seq.len() < p.window { 0 } else { (seq.len() - p.window) / stride + 1 };
        let mut part: Vec<Option<WindowSketch>> = Vec::with_capacity(n_win);
        (0..n_win)
            .into_par_iter()
            .map(|wi| {
                let s = wi * stride;
                let e = (s + p.window).min(seq.len());
                let w = &seq[s..e];
                // (hash, the k-mer is all uppercase)
                let mut hashes: Vec<(u64, bool)> = Vec::with_capacity(p.window);
                if p.uppercase_seed {
                    // up[i] = uppercase bases in w[..i]
                    let mut up = Vec::with_capacity(w.len() + 1);
                    up.push(0u32);
                    for &c in w {
                        up.push(up[up.len() - 1] + c.is_ascii_uppercase() as u32);
                    }
                    for_each_kmer(w, p.k, |i, h| hashes.push((h, (up[i + p.k] - up[i]) as usize == p.k)));
                } else {
                    for_each_kmer(w, p.k, |_, h| hashes.push((h, false)));
                }
                // ascending hash, an all-uppercase occurrence first, so dedup keeps it
                hashes.sort_unstable_by_key(|&(h, u)| (h, !u));
                hashes.dedup_by_key(|x| x.0);
                if hashes.len() < p.min_kmers {
                    return None;
                }
                hashes.truncate(p.sketch);
                let upper: Vec<bool> = if p.uppercase_seed {
                    let mut u = Vec::with_capacity(hashes.len());
                    u.extend(hashes.iter().map(|x| x.1));
                    u
                } else {
                    Vec::new()
                };
                // a fresh exact-size buffer (collecting in place would keep the k-mer buffer's capacity)
                let mut sketch = Vec::with_capacity(hashes.len());
                sketch.extend(hashes.iter().map(|x| x.0));
                let hashes = sketch;
                Some(WindowSketch {
                    seq: si as u32,
                    start: s as u64,
                    end: e as u64,
                    sketch: hashes,
                    upper,
                })
            })
            .collect_into_vec(&mut part);
        all.extend(part.into_iter().flatten());
    }
    (all, infos)
}

/// Inverted index over sketch values + pair counting.
///
/// Implementation: flatten (hash, window) entries and parallel-sort -- group iteration then streams linearly
/// through memory instead of chasing one allocated Vec per hash, and `par_sort_unstable` uses all cores.
///
/// Vendored changes (2026-10-13): each window's sketch is released once its entries are written (nothing reads it
/// again), and the pairs are counted in `p.count_shards` independent shards, so a whole genome's pair map is never
/// held at once. Candidates come back in shard order; the caller sorts them.
pub fn find_candidates(wins: &mut [WindowSketch], p: &SeedParams) -> Vec<Candidate> {
    // (sequence, start) per window: overlapping windows of one sequence (stride < window) share sequence, not a copy
    let at: Vec<(u32, u64)> = wins.iter().map(|w| (w.seq, w.start)).collect();
    let mut entries: Vec<(u64, u32, bool)> = Vec::with_capacity(wins.iter().map(|w| w.sketch.len()).sum());
    for (id, w) in wins.iter_mut().enumerate() {
        for (j, &h) in w.sketch.iter().enumerate() {
            entries.push((h, id as u32, w.upper.get(j).copied().unwrap_or(false)));
        }
        w.sketch = Vec::new();
        w.upper = Vec::new();
    }
    entries.par_sort_unstable();
    let shards = p.count_shards.max(1) as u32;
    let per_shard: Vec<Vec<Candidate>> =
        (0..shards).into_par_iter().map(|s| count_shard(&entries, &at, s, shards, p)).collect();
    per_shard.into_iter().flatten().collect()
}

/// The candidates whose lower window index is `shard` modulo `shards`. Inside one hash group the entries are sorted
/// by window index and a window holds a hash once (its sketch is deduplicated), so `a < b` whenever x < y: the
/// shard test can be made once per x.
fn count_shard(entries: &[(u64, u32, bool)], at: &[(u32, u64)], shard: u32, shards: u32, p: &SeedParams) -> Vec<Candidate> {
    let mut counts: FxHashMap<(u32, u32), u32> =
        FxHashMap::with_capacity_and_hasher(entries.len() / 64 / shards as usize, Default::default());
    let n = entries.len();
    let mut i = 0;
    while i < n {
        let h = entries[i].0;
        let mut j = i + 1;
        while j < n && entries[j].0 == h {
            j += 1;
        }
        let m = j - i;
        if m >= 2 && m <= p.max_postings {
            for x in i..j {
                let a = entries[x].1;
                if a % shards != shard {
                    continue;
                }
                for y in (x + 1)..j {
                    let b = entries[y].1;
                    debug_assert!(a < b);
                    let ((sa, xa), (sb, xb)) = (at[a as usize], at[b as usize]);
                    if sa == sb && xa.abs_diff(xb) < p.window as u64 {
                        continue;
                    }
                    // low 31 bits: shared hashes; bit 31: one of them is all uppercase in both windows
                    let both_up = ((entries[x].2 && entries[y].2) as u32) << 31;
                    let c = counts.entry((a, b)).or_insert(0);
                    *c = (*c + 1) | both_up;
                }
            }
        }
        i = j;
    }
    counts
        .into_iter()
        .filter(|(_, c)| (c & 0x7fff_ffff) >= p.min_shared as u32)
        .map(|((a, b), c)| Candidate { a: a as usize, b: b as usize, shared: c & 0x7fff_ffff, upper_seed: c >> 31 == 1 })
        .collect()
}

/// Rough Jaccard estimate from the number of shared bottom-S elements.
pub fn est_jaccard(shared: u32, sketch: usize) -> f64 {
    shared as f64 / (2.0 * sketch as f64 - shared as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(shards: usize) -> SeedParams {
        SeedParams { k: 16, window: 2000, sketch: 8, min_kmers: 1, max_postings: 3, min_shared: 2, count_shards: shards, stride: 2000, uppercase_seed: false }
    }

    fn win(id: u64, sketch: &[u64]) -> WindowSketch {
        WindowSketch { seq: 0, start: id * 2000, end: id * 2000 + 2000, sketch: sketch.to_vec(), upper: Vec::new() }
    }

    /// Hash 9 sits in 4 windows (> max_postings 3) and pairs nothing. Shared counts: (0,1) = {1,2,3} = 3;
    /// (0,2) = {3} = 1; (1,2) = {3,4} = 2; (3,4) = {5} = 1. With min_shared 2: (0,1) and (1,2).
    fn windows() -> Vec<WindowSketch> {
        vec![win(0, &[1, 2, 3, 9]), win(1, &[1, 2, 3, 4, 9]), win(2, &[3, 4, 9]), win(3, &[5, 9]), win(4, &[5, 6])]
    }

    fn sorted(c: Vec<Candidate>) -> Vec<(usize, usize, u32)> {
        let mut v: Vec<(usize, usize, u32)> = c.into_iter().map(|c| (c.a, c.b, c.shared)).collect();
        v.sort();
        v
    }

    #[test]
    fn sharded_counting_gives_the_one_map_candidates() {
        let one = sorted(find_candidates(&mut windows(), &params(1)));
        assert_eq!(one, vec![(0, 1, 3), (1, 2, 2)]);
        for s in [2, 3, 7, 64, 256] {
            assert_eq!(sorted(find_candidates(&mut windows(), &params(s))), one, "{s} shards");
        }
    }

    #[test]
    fn counting_releases_every_sketch() {
        let mut w = windows();
        find_candidates(&mut w, &params(4));
        assert!(w.iter().all(|x| x.sketch.is_empty()));
    }

    #[test]
    fn window_sketches_keep_only_their_hashes() {
        // vendored change (2026-10-13): the sandbox kept each window's 2,000-hash buffer (16 KB) behind a 256-hash
        // sketch, ~25 GB on a human genome
        let mut x = 7u64;
        let seq: Vec<u8> = (0..20_000)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 62) as usize]
            })
            .collect();
        let p = SeedParams { k: 16, window: 2000, sketch: 256, min_kmers: 32, max_postings: 100, min_shared: 40, count_shards: 1, stride: 2000, uppercase_seed: false };
        let (wins, _) = build_windows(&[("s".to_string(), seq)], &p);
        assert_eq!(wins.len(), 10);
        assert!(wins.iter().all(|w| w.sketch.len() == 256 && w.sketch.capacity() == 256));
    }

    fn synth(seed: u64, n: usize) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 62) as usize]
            })
            .collect()
    }

    #[test]
    fn a_copy_half_a_window_off_the_grid_is_seeded_with_half_window_stride() {
        // a 4 kb segment at 10 kb (on the 2 kb grid) and a ~91% copy at 41 kb (1 kb off it)
        let seg = synth(5, 4_000);
        let mut y = 99u64;
        let copy: Vec<u8> = seg
            .iter()
            .map(|&b| {
                y = y.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                if (y >> 32) % 100 < 12 { b"ACGT"[((y >> 40) % 4) as usize] } else { b }
            })
            .collect();
        let mut seq = synth(1, 60_000);
        seq[10_000..14_000].copy_from_slice(&seg);
        seq[41_000..45_000].copy_from_slice(&copy);
        let seqs = [("s".to_string(), seq)];
        let pairs_copy = |stride: usize, max_postings: usize| {
            let p = SeedParams { k: 16, window: 2000, stride, sketch: 256, min_kmers: 32, max_postings, min_shared: 40, count_shards: 1, uppercase_seed: false };
            let (mut wins, _) = build_windows(&seqs, &p);
            let c = find_candidates(&mut wins, &p);
            c.iter()
                .filter(|c| wins[c.a].start < 14_000 && wins[c.a].end > 10_000 && wins[c.b].start < 45_000 && wins[c.b].end > 41_000)
                .count()
        };
        assert_eq!(pairs_copy(2000, 100), 0, "scan v1's grid misses it (the reason for v2)");
        assert!(pairs_copy(1000, 200) > 0, "the half-window grid seeds it");
    }

    #[test]
    fn overlapping_windows_of_one_locus_are_not_candidates() {
        let p = SeedParams { k: 16, window: 2000, stride: 1000, sketch: 256, min_kmers: 32, max_postings: 200, min_shared: 40, count_shards: 1, uppercase_seed: false };
        let (mut wins, _) = build_windows(&[("s".to_string(), synth(3, 20_000))], &p);
        assert_eq!(wins.len(), 19);
        assert!(find_candidates(&mut wins, &p).is_empty());
    }

    #[test]
    fn only_windows_sharing_an_unmasked_kmer_are_upper_seeds() {
        // a masked 2 kb unit twice (windows 0 and 5) and an unmasked one twice (windows 10 and 15)
        let masked = synth(61, 2_000).to_ascii_lowercase();
        let unmasked = synth(62, 2_000);
        let mut seq = synth(8, 34_000);
        seq[0..2_000].copy_from_slice(&masked);
        seq[10_000..12_000].copy_from_slice(&masked);
        seq[20_000..22_000].copy_from_slice(&unmasked);
        seq[30_000..32_000].copy_from_slice(&unmasked);
        for (upper, want) in [(false, vec![]), (true, vec![(20_000u64, 30_000u64)])] {
            let p = SeedParams { k: 16, window: 2000, stride: 2000, sketch: 256, min_kmers: 32, max_postings: 100, min_shared: 40, count_shards: 1, uppercase_seed: upper };
            let (mut wins, _) = build_windows(&[("s".to_string(), seq.clone())], &p);
            let cands = find_candidates(&mut wins, &p);
            let mut all: Vec<(u64, u64)> = cands.iter().map(|c| (wins[c.a].start, wins[c.b].start)).collect();
            all.sort();
            assert_eq!(all, vec![(0, 10_000), (20_000, 30_000)], "both units are candidates");
            let mut got: Vec<(u64, u64)> =
                cands.iter().filter(|c| c.upper_seed).map(|c| (wins[c.a].start, wins[c.b].start)).collect();
            got.sort();
            assert_eq!(got, want, "uppercase_seed {upper}");
        }
    }

}
