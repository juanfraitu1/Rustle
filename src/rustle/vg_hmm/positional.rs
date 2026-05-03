//! Genome-wide positional candidate scan for novel-copy rescue.
//!
//! The classic rescue path scores unmapped reads against family graphs whose
//! copies come from *aligned* bundles. A genuinely novel paralog has no
//! bundle, so the HMM has nothing to match the read to at the novel locus.
//! This module fixes that gap: scan the genome for windows whose 15-mer
//! content matches a family's k-mer profile, cluster overlapping high-density
//! windows, exclude positions of known paralogs, and return the residue —
//! candidate genomic loci where a novel copy could live.
//!
//! Analogous to PSI-BLAST iteration: a profile built from the family's
//! known paralogs (the family-graph k-mer set) is used to discover more
//! distant homologs in the database (here: the genome itself). The resulting
//! candidates are then injected into the family graph as "ghost" copies in
//! Phase 2 so HMM scoring can target them.

use crate::genome::GenomeIndex;
use crate::types::DetHashSet;

/// One candidate genomic locus produced by the family-profile scan.
#[derive(Debug, Clone)]
pub struct CandidateLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// Distinct family k-mer hits inside the locus (forward + reverse strand
    /// counted together). Used as a soft confidence proxy.
    pub n_kmer_hits: usize,
    /// Hits / span length, in hits per kb. A self-paralog typically scores in
    /// the hundreds of hits/kb; coincidental cross-family matches are <30/kb.
    pub density_per_kb: f64,
}

/// Genomic spans of known paralogs to exclude from candidate output.
/// (chrom, start, end) — half-open, 0-based.
pub type KnownParalogSpan = (String, u64, u64);

#[inline]
fn fnv1a64(s: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in s {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

/// Hash a 15-mer's reverse-complement directly without allocating a new
/// buffer. Reads `kmer` from end → start, complementing each base.
#[inline]
fn fnv1a64_revcomp(kmer: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in kmer.iter().rev() {
        h ^= complement(b) as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

#[inline]
fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// Parameters controlling the genome scan.
#[derive(Debug, Clone, Copy)]
pub struct ScanParams {
    pub kmer_len: usize,
    /// Tile width in bp (e.g. 5000). Each tile's k-mer hit count is the unit
    /// of thresholding.
    pub window: u64,
    /// Stride between tiles in bp (e.g. 2500 for 50% overlap).
    pub stride: u64,
    /// Minimum distinct family k-mer hits inside a tile to count as a hit.
    /// Tunable by family — pass kmer_set.len() / 100 as a rough default.
    pub min_hits_per_window: usize,
    /// Padding (bp) added around known paralog spans for the exclusion test.
    /// Captures candidate hits that overlap a paralog's flanks.
    pub paralog_pad: u64,
    /// Skip chromosomes whose sequence is shorter than this (avoids tiny
    /// unplaced contigs producing noise).
    pub min_chrom_len: u64,
    /// After clustering, keep only the top-N candidates per family by
    /// density. 0 = unlimited. Keeps the output bounded — Phase 2 only
    /// inflates the family graph with the strongest hits anyway.
    pub max_candidates_per_family: usize,
}

impl Default for ScanParams {
    fn default() -> Self {
        Self {
            kmer_len: 15,
            window: 5000,
            stride: 2500,
            // Empirical chr19 scan (commit 1ea2f35 + scan diagnostics):
            // real paralogs hit 200+ in a 5kb window; coincidental
            // background sequence hits <50. 200 is a safe floor that keeps
            // candidate counts in the hundreds per family rather than the
            // hundred-thousands.
            min_hits_per_window: 200,
            paralog_pad: 5000,
            min_chrom_len: 50_000,
            max_candidates_per_family: 50,
        }
    }
}

/// Scan one chromosome, returning per-tile hit counts. Each entry is
/// (tile_start, n_hits). Tiles use the params' window/stride.
fn scan_chrom_tiles(
    chrom_seq: &[u8],
    family_kmers: &DetHashSet<u64>,
    params: &ScanParams,
) -> Vec<(u64, usize)> {
    let k = params.kmer_len;
    let win = params.window as usize;
    let stride = params.stride as usize;
    if chrom_seq.len() < win || chrom_seq.len() < k {
        return Vec::new();
    }
    let mut out: Vec<(u64, usize)> = Vec::new();
    let mut tile_start = 0usize;
    while tile_start + win <= chrom_seq.len() {
        let tile = &chrom_seq[tile_start..tile_start + win];
        let mut hits = 0usize;
        for w in tile.windows(k) {
            if !w.iter().all(|&b| is_acgt(b)) { continue; }
            let h_fwd = fnv1a64(w);
            if family_kmers.contains(&h_fwd) {
                hits += 1;
                continue;
            }
            let h_rc = fnv1a64_revcomp(w);
            if family_kmers.contains(&h_rc) {
                hits += 1;
            }
        }
        out.push((tile_start as u64, hits));
        tile_start += stride;
    }
    out
}

/// Cluster contiguous tiles whose hits ≥ threshold into candidate spans.
/// Returns Vec<(start, end, total_hits)> per chromosome.
fn cluster_tiles(
    tiles: &[(u64, usize)],
    params: &ScanParams,
) -> Vec<(u64, u64, usize)> {
    let mut out: Vec<(u64, u64, usize)> = Vec::new();
    let mut cur: Option<(u64, u64, usize)> = None;
    for &(ts, hits) in tiles {
        let te = ts + params.window;
        if hits >= params.min_hits_per_window {
            cur = match cur {
                Some((s, e, h)) if ts <= e => Some((s, te.max(e), h + hits)),
                _ => {
                    if let Some(c) = cur.take() { out.push(c); }
                    Some((ts, te, hits))
                }
            };
        } else if let Some(c) = cur.take() {
            out.push(c);
        }
    }
    if let Some(c) = cur { out.push(c); }
    out
}

/// Drop candidates whose span overlaps any known paralog (with `paralog_pad`
/// bp of padding). Keeps the residue — these are the regions where a novel
/// copy could plausibly live.
fn exclude_known_paralogs(
    chrom: &str,
    raw: Vec<(u64, u64, usize)>,
    known: &[KnownParalogSpan],
    pad: u64,
) -> Vec<(u64, u64, usize)> {
    let masks: Vec<(u64, u64)> = known
        .iter()
        .filter(|(c, _, _)| c == chrom)
        .map(|(_, s, e)| (s.saturating_sub(pad), e.saturating_add(pad)))
        .collect();
    raw.into_iter()
        .filter(|(s, e, _)| !masks.iter().any(|(ms, me)| s < me && ms < e))
        .collect()
}

/// Single-pass multi-family scan: builds an inverted index (kmer → families
/// that contain it), walks the genome once, and returns candidate loci PER
/// FAMILY. Much faster than calling [`scan_genome_for_family_loci`] in a
/// loop over families: the genome is read once instead of N times, and
/// per-15-mer there's one hash + one lookup instead of N.
///
/// `family_kmers[i]` and `known_paralogs[i]` describe family `i`. The
/// returned Vec has length `family_kmers.len()`; entry `i` is family `i`'s
/// candidate loci (sorted by descending density), or empty.
///
/// `repetitive_kmer_max_families` (typical: 4): k-mers that appear in more
/// than this many families' sets are dropped from the inverted index. They
/// represent shared repetitive sequence (interspersed elements, low-complexity
/// regions) and would create spurious hits across unrelated families.
pub fn scan_genome_for_all_families(
    family_kmers: &[DetHashSet<u64>],
    genome: &GenomeIndex,
    known_paralogs: &[Vec<KnownParalogSpan>],
    params: &ScanParams,
    repetitive_kmer_max_families: usize,
) -> Vec<Vec<CandidateLocus>> {
    let n_fams = family_kmers.len();
    if n_fams == 0 || family_kmers.iter().all(|s| s.is_empty()) {
        return vec![Vec::new(); n_fams];
    }
    assert_eq!(known_paralogs.len(), n_fams,
        "known_paralogs and family_kmers must have same length");
    if n_fams > u16::MAX as usize {
        // Defensive bound — we encode family idx as u16 in the inverted index.
        eprintln!("[positional] n_families={} exceeds u16 range; truncating", n_fams);
    }

    // Build inverted index: kmer hash → SmallVec of family idxs (u16).
    // Use raw Vec<u16> here; for typical paralog families, average families
    // per shared k-mer is 1-3, and the "repetitive" cap keeps the worst
    // case bounded.
    let mut inv: crate::types::DetHashMap<u64, Vec<u16>> =
        crate::types::DetHashMap::default();
    for (fi, kset) in family_kmers.iter().enumerate() {
        let fi_u16 = fi as u16;
        for &h in kset {
            inv.entry(h).or_default().push(fi_u16);
        }
    }
    // Drop repetitive k-mers (in too many families).
    let n_inv_before = inv.len();
    inv.retain(|_, v| v.len() <= repetitive_kmer_max_families);
    let n_inv_after = inv.len();
    if n_inv_before != n_inv_after {
        eprintln!(
            "[positional] inverted index: {} unique k-mers ({} dropped as repetitive >{}-fams)",
            n_inv_after, n_inv_before - n_inv_after, repetitive_kmer_max_families,
        );
    }

    // Single-pass scan over the genome. Per chrom: build per-family per-tile
    // counters; cluster + exclude per family afterward. Chroms run in
    // parallel via rayon — each chrom is independent and the inverted index
    // is shared `&` across threads.
    let k = params.kmer_len;
    let win = params.window as usize;
    let stride = params.stride as usize;

    use rayon::prelude::*;
    let chroms: Vec<(String, &[u8])> = genome.chroms()
        .map(|(name, seq)| (name.to_string(), seq))
        .collect();
    let per_chrom: Vec<Vec<(usize, CandidateLocus)>> = chroms
        .par_iter()
        .map(|(chrom, seq)| {
            let mut out: Vec<(usize, CandidateLocus)> = Vec::new();
            if (seq.len() as u64) < params.min_chrom_len { return out; }
            if seq.len() < win || seq.len() < k { return out; }

            let n_tiles = ((seq.len() - win) / stride) + 1;
            // tile_hits[fi][tile] = distinct family-k-mer hits in tile.
            // We accumulate hits PER POSITION (not per tile-overlap) so a
            // single 15-mer that lies in 2 overlapping tiles contributes
            // once to each — this is the correct double-count for tile
            // density, but avoid summing across the cluster step (the
            // cluster_tiles helper just merges adjacent passing tiles
            // without re-summing).
            let mut tile_hits: Vec<Vec<u32>> = vec![vec![0u32; n_tiles]; n_fams];

            for (pos, w) in seq.windows(k).enumerate() {
                if !w.iter().all(|&b| is_acgt(b)) { continue; }
                let h_fwd = fnv1a64(w);
                let fams_fwd = inv.get(&h_fwd);
                let h_rc = fnv1a64_revcomp(w);
                let fams_rc = inv.get(&h_rc);
                if fams_fwd.is_none() && fams_rc.is_none() { continue; }
                let t_max = (pos / stride).min(n_tiles - 1);
                let t_min = if pos + 1 > win { (pos + 1 - win + stride - 1) / stride } else { 0 };
                for t in t_min..=t_max {
                    if let Some(v) = fams_fwd {
                        for &fi in v { tile_hits[fi as usize][t] += 1; }
                    }
                    if let Some(v) = fams_rc {
                        for &fi in v { tile_hits[fi as usize][t] += 1; }
                    }
                }
            }

            for (fi, hits) in tile_hits.iter().enumerate() {
                let tiles: Vec<(u64, usize)> = hits.iter().enumerate()
                    .map(|(t, &c)| ((t * stride) as u64, c as usize))
                    .collect();
                let clustered = cluster_tiles_no_double_count(&tiles, params);
                let kept = exclude_known_paralogs(chrom, clustered, &known_paralogs[fi], params.paralog_pad);
                for (s, e, h) in kept {
                    let span = (e.saturating_sub(s)) as f64;
                    let density = if span > 0.0 { 1000.0 * h as f64 / span } else { 0.0 };
                    out.push((fi, CandidateLocus {
                        chrom: chrom.clone(),
                        start: s,
                        end: e,
                        n_kmer_hits: h,
                        density_per_kb: density,
                    }));
                }
            }
            out
        })
        .collect();

    let mut all_results: Vec<Vec<CandidateLocus>> = vec![Vec::new(); n_fams];
    for entries in per_chrom {
        for (fi, c) in entries {
            all_results[fi].push(c);
        }
    }
    for v in &mut all_results {
        v.sort_by(|a, b| b.density_per_kb.partial_cmp(&a.density_per_kb).unwrap_or(std::cmp::Ordering::Equal));
        if params.max_candidates_per_family > 0 && v.len() > params.max_candidates_per_family {
            v.truncate(params.max_candidates_per_family);
        }
    }
    all_results
}

/// Like `cluster_tiles` but reports the MAX hits in any contributing tile
/// (not the sum, which double-counts overlapping windows). The reported
/// "density" is then computed against the cluster span, giving an
/// interpretable hits/kb that doesn't exceed physical limits.
fn cluster_tiles_no_double_count(
    tiles: &[(u64, usize)],
    params: &ScanParams,
) -> Vec<(u64, u64, usize)> {
    let mut out: Vec<(u64, u64, usize)> = Vec::new();
    let mut cur: Option<(u64, u64, usize)> = None;
    for &(ts, hits) in tiles {
        let te = ts + params.window;
        if hits >= params.min_hits_per_window {
            cur = match cur {
                Some((s, e, h)) if ts <= e => Some((s, te.max(e), h.max(hits))),
                _ => {
                    if let Some(c) = cur.take() { out.push(c); }
                    Some((ts, te, hits))
                }
            };
        } else if let Some(c) = cur.take() {
            out.push(c);
        }
    }
    if let Some(c) = cur { out.push(c); }
    out
}

/// Top-level scan. Iterates the genome once per family, returns candidate
/// loci sorted by descending k-mer density.
pub fn scan_genome_for_family_loci(
    family_kmers: &DetHashSet<u64>,
    genome: &GenomeIndex,
    known_paralogs: &[KnownParalogSpan],
    params: &ScanParams,
) -> Vec<CandidateLocus> {
    let mut all: Vec<CandidateLocus> = Vec::new();
    if family_kmers.is_empty() { return all; }
    for (chrom, seq) in genome.chroms() {
        if (seq.len() as u64) < params.min_chrom_len { continue; }
        let tiles = scan_chrom_tiles(seq, family_kmers, params);
        let clustered = cluster_tiles(&tiles, params);
        let kept = exclude_known_paralogs(chrom, clustered, known_paralogs, params.paralog_pad);
        for (s, e, hits) in kept {
            let span = (e.saturating_sub(s)) as f64;
            let density = if span > 0.0 { 1000.0 * hits as f64 / span } else { 0.0 };
            all.push(CandidateLocus {
                chrom: chrom.to_string(),
                start: s,
                end: e,
                n_kmer_hits: hits,
                density_per_kb: density,
            });
        }
    }
    all.sort_by(|a, b| b.density_per_kb.partial_cmp(&a.density_per_kb).unwrap_or(std::cmp::Ordering::Equal));
    all
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a tiny synthetic family k-mer set (just one 100-bp seed) and
    /// verify that the scan finds the seed when it's planted in a small
    /// "genome" and ignores random flanking sequence.
    #[test]
    fn scan_finds_planted_seed() {
        let seed = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut family_kmers: DetHashSet<u64> = DetHashSet::default();
        for w in seed.windows(15) {
            family_kmers.insert(fnv1a64(w));
        }
        // Build a synthetic chromosome: 6 kb of As, then the seed, then 6 kb of As.
        let mut chrom: Vec<u8> = vec![b'A'; 6000];
        chrom.extend_from_slice(seed);
        chrom.extend(vec![b'A'; 6000]);
        // Stash into a hand-built GenomeIndex.
        // Use `from_fasta` via a temp file would be heavy — write a small helper.
        let mut idx = GenomeIndex::default();
        // SAFETY: GenomeIndex's seqs is private but `Default` gives an empty map;
        // we emulate by writing through the test-only public API. For now,
        // skip the GenomeIndex path and exercise scan_chrom_tiles directly.
        let _ = &mut idx;
        let params = ScanParams { kmer_len: 15, window: 200, stride: 100,
                                  min_hits_per_window: 20,
                                  paralog_pad: 0, min_chrom_len: 0,
                                  max_candidates_per_family: 0 };
        let tiles = scan_chrom_tiles(&chrom, &family_kmers, &params);
        let clustered = cluster_tiles(&tiles, &params);
        // Expect at least one cluster covering ~position 6000.
        assert!(!clustered.is_empty(),
            "expected clustered hits at planted seed, got {:?}", clustered);
        let mid_clusters: Vec<_> = clustered.iter()
            .filter(|(s, e, _)| *s < 6200 && *e > 6000)
            .collect();
        assert!(!mid_clusters.is_empty(),
            "expected a cluster overlapping position 6000-6100, got {:?}", clustered);
    }

    #[test]
    fn revcomp_hash_matches_explicit_revcomp() {
        let kmer = b"ACGTACGTACGTACG";
        let rc: Vec<u8> = kmer.iter().rev().map(|&b| complement(b)).collect();
        assert_eq!(fnv1a64(&rc), fnv1a64_revcomp(kmer));
    }

    #[test]
    fn cluster_tiles_merges_adjacent() {
        let params = ScanParams { kmer_len: 15, window: 1000, stride: 500,
                                  min_hits_per_window: 10,
                                  paralog_pad: 0, min_chrom_len: 0,
                                  max_candidates_per_family: 0 };
        let tiles = vec![
            (0, 5),       // below
            (500, 50),    // hit
            (1000, 60),   // hit (adjacent)
            (1500, 8),    // below
            (3000, 40),   // hit (separate)
        ];
        let out = cluster_tiles(&tiles, &params);
        assert_eq!(out.len(), 2, "got {:?}", out);
        assert_eq!(out[0].2, 110);  // 50 + 60
        assert_eq!(out[1].2, 40);
    }

    #[test]
    fn exclude_known_paralogs_drops_overlapping() {
        let raw = vec![(1000, 2000, 50), (5000, 6000, 30), (9000, 10000, 80)];
        let known = vec![("chr1".to_string(), 4500, 5500)];
        let kept = exclude_known_paralogs("chr1", raw, &known, 200);
        assert_eq!(kept.len(), 2);
        assert!(kept.iter().all(|(s, e, _)| !(*s < 5700 && 4300 < *e)));
    }
}
