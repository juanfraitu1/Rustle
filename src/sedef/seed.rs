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
/// Implementation: flatten (hash, window) entries and parallel-sort --
/// group iteration then streams linearly through memory instead of
/// chasing one allocated Vec per hash, and `par_sort_unstable` uses all
/// cores.
pub fn find_candidates(wins: &[WindowSketch], p: &SeedParams) -> Vec<Candidate> {
    let mut entries: Vec<(u64, u32)> = Vec::with_capacity(wins.len() * p.sketch);
    for (id, w) in wins.iter().enumerate() {
        for &h in &w.sketch {
            entries.push((h, id as u32));
        }
    }
    entries.par_sort_unstable();

    let mut counts: FxHashMap<(u32, u32), u32> =
        FxHashMap::with_capacity_and_hasher(entries.len() / 64, Default::default());
    let mut i = 0;
    let n = entries.len();
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
                for y in (x + 1)..j {
                    let b = entries[y].1;
                    let key = if a < b { (a, b) } else { (b, a) };
                    match counts.get_mut(&key) {
                        Some(c) => *c += 1,
                        None => {
                            counts.insert(key, 1);
                        }
                    }
                }
            }
        }
        i = j;
    }

    counts
        .into_iter()
        .filter(|(_, c)| *c >= p.min_shared as u32)
        .map(|((a, b), c)| Candidate {
            a: a as usize,
            b: b as usize,
            shared: c,
        })
        .collect()
}

/// Rough Jaccard estimate from the number of shared bottom-S elements.
pub fn est_jaccard(shared: u32, sketch: usize) -> f64 {
    shared as f64 / (2.0 * sketch as f64 - shared as f64)
}
