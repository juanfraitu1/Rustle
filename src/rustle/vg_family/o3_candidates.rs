//! O3 candidate copies (spec docs/superpowers/specs/2026-10-02-o3-candidates-design.md): the reference-absent-copy chain of
//! PREREG_rna_allele_haplotype_count Amendments 7-11 without IsoCon. Pure functions here; BAM passes and minimap2 calls in the binary.
//!
//! **STATUS:** TEST-ONLY — no caller yet (the `o3_candidates` binary and the driver stage land in later tasks of the plan; flips to OPT-IN then)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use std::collections::HashMap;

pub const KMER_K: usize = 31;
pub const SKETCH_W: usize = 5;
pub const MIN_UNMAPPED_LEN: usize = 300;
const ATTRIB_MIN_FRAC: f64 = 0.30;
const ATTRIB_MIN_RATIO: f64 = 2.0;

fn code(b: u8) -> Option<u64> { match b { b'A' | b'a' => Some(0), b'C' | b'c' => Some(1), b'G' | b'g' => Some(2), b'T' | b't' => Some(3), _ => None } }

/// Canonical (min of forward / reverse-complement) 2-bit k-mers, k <= 32; windows with N are skipped.
pub fn canonical_kmers(seq: &[u8], k: usize) -> Vec<u64> {
    let mask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    let mut out = Vec::with_capacity(seq.len().saturating_sub(k) + 1);
    for &b in seq {
        match code(b) {
            Some(c) => { fw = ((fw << 2) | c) & mask; rv = (rv >> 2) | ((3 - c) << (2 * (k - 1))); valid += 1; }
            None => { valid = 0; fw = 0; rv = 0; }
        }
        if valid >= k { out.push(fw.min(rv)); }
    }
    out
}

fn mix(x: u64) -> u64 { let mut z = x.wrapping_add(0x9E3779B97F4A7C15); z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }

/// (k, w) minimizers over the canonical k-mers: the smallest hashed k-mer of every window of w consecutive k-mers, deduplicated
/// consecutively; sorted and deduplicated so sketches compare by merge.
pub fn minimizer_sketch(seq: &[u8], k: usize, w: usize) -> Vec<u64> {
    let km: Vec<u64> = canonical_kmers(seq, k).into_iter().map(mix).collect();
    let mut out = Vec::new();
    if km.len() < w { out.extend(km.iter().copied().min()); }
    for win in km.windows(w) { let m = *win.iter().min().unwrap(); if out.last() != Some(&m) { out.push(m); } }
    out.sort_unstable(); out.dedup(); out
}

/// Shared minimizers / size of the smaller sketch (both sorted, deduplicated).
pub fn sketch_share(a: &[u64], b: &[u64]) -> f64 {
    let (mut i, mut j, mut shared) = (0, 0, 0usize);
    while i < a.len() && j < b.len() { if a[i] == b[j] { shared += 1; i += 1; j += 1; } else if a[i] < b[j] { i += 1; } else { j += 1; } }
    shared as f64 / a.len().min(b.len()).max(1) as f64
}

/// Canonical 31-mers of every copy sequence -> the families carrying them (k-mers in > max_families families dropped as repeats).
pub struct FamilyKmerIndex { map: HashMap<u64, Vec<String>>, k: usize }
impl FamilyKmerIndex {
    pub fn build(copies: &[(String, Vec<u8>)], k: usize, max_families: usize) -> Self {
        let mut map: HashMap<u64, Vec<String>> = HashMap::new();
        for (fam, seq) in copies {
            for km in canonical_kmers(seq, k) { let v = map.entry(km).or_default(); if !v.contains(fam) { v.push(fam.clone()); } }
        }
        map.retain(|_, v| v.len() <= max_families);
        FamilyKmerIndex { map, k }
    }
    /// The family with the most k-mer hits when >= 30% of the read's k-mers hit it and it leads the runner-up >= 2x; reads < 300 bp: None.
    pub fn attribute(&self, read: &[u8]) -> Option<String> {
        if read.len() < MIN_UNMAPPED_LEN { return None; }
        let kms = canonical_kmers(read, self.k);
        let mut hits: HashMap<&str, usize> = HashMap::new();
        for km in &kms { if let Some(f) = self.map.get(km) { for fam in f { *hits.entry(fam).or_insert(0) += 1; } } }
        let mut v: Vec<(&str, usize)> = hits.into_iter().collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
        let (best, n) = *v.first()?;
        let second = v.get(1).map(|x| x.1).unwrap_or(0);
        if (n as f64) < ATTRIB_MIN_FRAC * kms.len() as f64 { return None; }
        if second > 0 && (n as f64) < ATTRIB_MIN_RATIO * second as f64 { return None; }
        Some(best.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed; (0..n).map(|_| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b"ACGT"[(x % 4) as usize] }).collect()
    }
    fn mutate(s: &[u8], rate_per_kb: usize, seed: u64) -> Vec<u8> {
        // substitutes every (1000 / rate_per_kb)-th base; the new base is always different from the original
        let mut v = s.to_vec(); let mut x = seed.max(1);
        for i in (0..v.len()).step_by(1000 / rate_per_kb.max(1)) {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            let orig = b"ACGT".iter().position(|&c| c == v[i]).expect("ACGT input");
            v[i] = b"ACGT"[(orig + 1 + (x % 3) as usize) % 4];   // index + 1..=3 (mod 4) never returns the original index
        }
        v
    }
    #[test]
    fn canonical_kmers_are_strand_free() {
        let s = rand_seq(500, 7);
        let rc = crate::vg_family::seq_utils::reverse_complement(&s);
        let mut a = canonical_kmers(&s, 31); let mut b = canonical_kmers(&rc, 31);
        a.sort(); b.sort();
        assert_eq!(a, b);
    }
    #[test]
    fn sketch_share_separates_copies_from_errors() {
        let copy_a = rand_seq(3000, 1);
        let same_copy_read = mutate(&copy_a, 2, 3);        // 0.2%: HiFi-like (every 500th base: 6 substitutions in 3 kb)
        let copy_b = mutate(&copy_a, 25, 5);               // 2.5% diverged paralog (every 40th base: 75 substitutions in 3 kb)
        let real_subs = |x: &[u8]| x.iter().zip(&copy_a).filter(|(p, q)| p != q).count();
        assert_eq!((real_subs(&same_copy_read), real_subs(&copy_b)), (6, 75), "every substitution must be real: exactly 0.2% and 2.5% of 3000 bp");
        let sa = minimizer_sketch(&copy_a, 31, 5);
        assert!(sketch_share(&sa, &minimizer_sketch(&same_copy_read, 31, 5)) > 0.8);
        assert!(sketch_share(&sa, &minimizer_sketch(&copy_b, 31, 5)) < 0.6);
    }
    #[test]
    fn attribution_needs_a_clear_winner_and_ignores_short_reads() {
        let fam1 = rand_seq(2000, 11); let fam2 = rand_seq(2000, 12);
        let idx = FamilyKmerIndex::build(&[("F1".into(), fam1.clone()), ("F2".into(), fam2.clone())], 31, 8);
        assert_eq!(idx.attribute(&mutate(&fam1, 2, 2)).as_deref(), Some("F1"));
        assert_eq!(idx.attribute(&rand_seq(1500, 99)), None);               // nothing hits
        let half: Vec<u8> = fam1[..1000].iter().chain(fam2[..1000].iter()).copied().collect();
        assert_eq!(idx.attribute(&half), None);                             // no 2x winner
        assert_eq!(idx.attribute(&fam1[..69]), None);                       // short_unmapped_reads_are_ignored
    }
}
