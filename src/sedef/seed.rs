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
    pub sketch: usize,
    /// Windows with fewer distinct k-mers than this are treated as pure
    /// satellite repeats and skipped (SDFinder/BISER do analogous repeat
    /// masking upstream).
    pub min_kmers: usize,
    /// Sketch values occurring in more windows than this are discarded
    /// (they carry no pairing signal -- classic too-common filter).
    pub max_postings: usize,
    pub min_shared: usize,
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
}

pub struct Candidate {
    pub a: usize,
    pub b: usize,
    pub shared: u32,
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
        let n_win = seq.len() / p.window;
        let mut part: Vec<Option<WindowSketch>> = Vec::with_capacity(n_win);
        (0..n_win)
            .into_par_iter()
            .map(|wi| {
                let s = wi * p.window;
                let e = (s + p.window).min(seq.len());
                let mut hashes: Vec<u64> = Vec::with_capacity(p.window);
                for_each_kmer(&seq[s..e], p.k, |_, h| hashes.push(h));
                hashes.sort_unstable();
                hashes.dedup();
                if hashes.len() < p.min_kmers {
                    return None;
                }
                hashes.truncate(p.sketch);
                // vendored change: release the k-mer buffer's spare capacity (the sandbox kept 16 KB per window)
                hashes.shrink_to_fit();
                Some(WindowSketch {
                    seq: si as u32,
                    start: s as u64,
                    end: e as u64,
                    sketch: hashes,
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
    let mut entries: Vec<(u64, u32)> = Vec::with_capacity(wins.iter().map(|w| w.sketch.len()).sum());
    for (id, w) in wins.iter_mut().enumerate() {
        for &h in &w.sketch {
            entries.push((h, id as u32));
        }
        w.sketch = Vec::new();
    }
    entries.par_sort_unstable();
    let shards = p.count_shards.max(1) as u32;
    let per_shard: Vec<Vec<Candidate>> =
        (0..shards).into_par_iter().map(|s| count_shard(&entries, s, shards, p)).collect();
    per_shard.into_iter().flatten().collect()
}

/// The candidates whose lower window index is `shard` modulo `shards`. Inside one hash group the entries are sorted
/// by window index and a window holds a hash once (its sketch is deduplicated), so `a < b` whenever x < y: the
/// shard test can be made once per x.
fn count_shard(entries: &[(u64, u32)], shard: u32, shards: u32, p: &SeedParams) -> Vec<Candidate> {
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
                    *counts.entry((a, b)).or_insert(0) += 1;
                }
            }
        }
        i = j;
    }
    counts
        .into_iter()
        .filter(|(_, c)| *c >= p.min_shared as u32)
        .map(|((a, b), c)| Candidate { a: a as usize, b: b as usize, shared: c })
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
        SeedParams { k: 16, window: 2000, sketch: 8, min_kmers: 1, max_postings: 3, min_shared: 2, count_shards: shards }
    }

    fn win(id: u64, sketch: &[u64]) -> WindowSketch {
        WindowSketch { seq: 0, start: id * 2000, end: id * 2000 + 2000, sketch: sketch.to_vec() }
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
        let p = SeedParams { k: 16, window: 2000, sketch: 256, min_kmers: 32, max_postings: 100, min_shared: 40, count_shards: 1 };
        let (wins, _) = build_windows(&[("s".to_string(), seq)], &p);
        assert_eq!(wins.len(), 10);
        assert!(wins.iter().all(|w| w.sketch.len() == 256 && w.sketch.capacity() == 256));
    }
}
