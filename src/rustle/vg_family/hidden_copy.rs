//! Detect gene-family copies PRESENT in the reads but ABSENT from the reference genome
//! (collapsed segdup / assembly gap / CNV / private duplication).
//!
//! A hidden copy's reads have no correct home, so they mismap to the closest sibling reference
//! copy carrying their PRIVATE SNPs — a COHERENT second haplotype the reference (one copy at
//! that locus) cannot explain. This detector finds that second haplotype among the locus's
//! PRIMARY alignments and FLAGS the discrepancy. Per the DAZ3 discipline it DETECTS and reports
//! evidence ("the reads imply ≥2 copies; the reference models 1") and ABSTAINS from placing or
//! fabricating the missing copy — it never manufactures a copy.
//!
//! Design = synthesis of an independent design panel (statistical / algorithmic / honesty
//! lenses), which converged on: PRIMARY-alignments-only matrix (the paralog-bleed firewall, since an
//! in-reference paralog's reads are primary at THEIR locus), candidate columns where the
//! non-reference allele frequency sits in a balanced band (excludes 0.5% sequencing error and
//! fixed differences), and a co-segregation/block test (a hidden copy's alt columns co-occur on
//! ONE read subset — distinguishing it from scattered heterozygous SNPs by requiring many).

/// One primary alignment's observation at the locus: its covered span [start, end) and the positions
/// where it carries a non-reference allele (a mismatch).
#[derive(Debug, Clone)]
pub struct ReadObs {
    pub start: u64,
    pub end: u64,
    pub alts: Vec<u64>,
}

#[derive(Debug, Clone)]
pub struct HiddenCopyParams {
    pub balanced_lo: f64,        // min alt-allele fraction for a candidate column (≫ error rate)
    pub balanced_hi: f64,        // max alt-allele fraction (above = fixed diff / ref error)
    pub min_depth: usize,        // min coverage at a candidate column
    pub min_alt_positions: usize,// min candidate columns to call a hidden copy (≫ a few hets)
    pub min_alt_reads: usize,    // min reads in the alt haplotype (the hidden copy's depth)
    pub share_hi: f64,           // a read joins H if alt at ≥ this fraction of candidate cols it covers
}

impl Default for HiddenCopyParams {
    fn default() -> Self {
        HiddenCopyParams {
            balanced_lo: 0.20,
            balanced_hi: 0.60,
            min_depth: 8,
            min_alt_positions: 12, // het firewall: a diploid het is 1-2 columns; a copy is dozens
            min_alt_reads: 5,
            share_hi: 0.60,
        }
    }
}

impl HiddenCopyParams {
    pub fn from_env() -> Self {
        let mut p = HiddenCopyParams::default();
        let getf = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<f64>().ok());
        let getu = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<usize>().ok());
        if let Some(v) = getf("RUSTLE_VG_HIDDEN_ALT_LO") { p.balanced_lo = v; }
        if let Some(v) = getf("RUSTLE_VG_HIDDEN_ALT_HI") { p.balanced_hi = v; }
        if let Some(v) = getu("RUSTLE_VG_HIDDEN_MIN_DEPTH") { p.min_depth = v; }
        if let Some(v) = getu("RUSTLE_VG_HIDDEN_MIN_POSITIONS") { p.min_alt_positions = v; }
        if let Some(v) = getu("RUSTLE_VG_HIDDEN_MIN_READS") { p.min_alt_reads = v; }
        p
    }
}

/// Evidence for a copy not in the reference. DETECT + FLAG only — no placement, no sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct HiddenCopyEvidence {
    pub n_primary_reads: usize,
    pub n_alt_positions: usize, // coherent second-haplotype columns
    pub n_alt_reads: usize,     // reads in the alt haplotype (the hidden copy's apparent depth)
    pub alt_read_fraction: f64, // n_alt_reads / n_primary_reads
    pub flagged: bool,          // evidence of an unmodeled copy at this locus
}

/// Pure detector over PRIMARY alignments at one reference-copy locus. The caller MUST pass primary
/// reads only (the paralog-bleed firewall). Deterministic; no I/O.
pub fn detect_hidden_copy(reads: &[ReadObs], p: &HiddenCopyParams) -> HiddenCopyEvidence {
    let n = reads.len();
    let none = HiddenCopyEvidence {
        n_primary_reads: n, n_alt_positions: 0, n_alt_reads: 0,
        alt_read_fraction: 0.0, flagged: false,
    };
    if n < p.min_depth {
        return none;
    }

    // Per-position alt-read count (only positions some read calls alt are candidates).
    let mut alt_count: crate::types::DetHashMap<u64, usize> = Default::default();
    for r in reads {
        for &pos in &r.alts {
            *alt_count.entry(pos).or_insert(0) += 1;
        }
    }

    // Candidate columns: balanced alt fraction over sufficient depth. Error (~0.5%) never
    // reaches balanced_lo; a fixed difference / reference error sits above balanced_hi.
    let mut candidates: Vec<u64> = Vec::new();
    for (&pos, &ac) in &alt_count {
        if ac < 2 {
            continue; // a singleton is error, not a haplotype
        }
        let cov = reads.iter().filter(|r| r.start <= pos && pos < r.end).count();
        if cov < p.min_depth {
            continue;
        }
        let frac = ac as f64 / cov as f64;
        if frac >= p.balanced_lo && frac <= p.balanced_hi {
            candidates.push(pos);
        }
    }
    candidates.sort_unstable();
    let n_alt_positions = candidates.len();

    // Co-segregation: a hidden COPY's candidate columns co-occur on ONE read subset (H). Partition
    // reads by their alt-share over the candidate columns they cover; H = the alt haplotype.
    let cand_set: crate::types::DetHashSet<u64> = candidates.iter().copied().collect();
    let n_alt_reads = reads.iter().filter(|r| {
        let covered = candidates.iter().filter(|&&pos| r.start <= pos && pos < r.end).count();
        if covered == 0 {
            return false;
        }
        let alt_at = r.alts.iter().filter(|pos| cand_set.contains(pos)).count();
        (alt_at as f64 / covered as f64) >= p.share_hi
    }).count();

    // Flag only with MANY co-segregating positions (≫ a het) AND a real alt-haplotype read group.
    let flagged = n_alt_positions >= p.min_alt_positions && n_alt_reads >= p.min_alt_reads;

    HiddenCopyEvidence {
        n_primary_reads: n,
        n_alt_positions,
        n_alt_reads,
        alt_read_fraction: if n > 0 { n_alt_reads as f64 / n as f64 } else { 0.0 },
        flagged,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p() -> HiddenCopyParams { HiddenCopyParams::default() }

    // n reads all spanning [0, span); `hap` reads carry alt at every position in `shared`,
    // the rest carry alt at `noise` random-but-distinct positions each (sequencing error).
    fn reads(n: usize, span: u64, hap: usize, shared: &[u64], noise_per_read: u64) -> Vec<ReadObs> {
        (0..n).map(|r| {
            let mut alts: Vec<u64> = if r < hap { shared.to_vec() } else { Vec::new() };
            // distinct error positions per read (no cross-read coherence)
            for k in 0..noise_per_read {
                alts.push(span - 1 - (r as u64 * 17 + k)); // deterministic, scattered, unique-ish
            }
            ReadObs { start: 0, end: span, alts }
        }).collect()
    }

    #[test]
    fn hidden_copy_is_flagged() {
        // 60 reads, 30 carry a 20-position shared alt haplotype → a hidden copy.
        let shared: Vec<u64> = (100..120).collect();
        let rs = reads(60, 2000, 30, &shared, 3);
        let e = detect_hidden_copy(&rs, &p());
        assert!(e.flagged);
        assert_eq!(e.n_alt_positions, 20);
        assert_eq!(e.n_alt_reads, 30);
        assert!((e.alt_read_fraction - 0.5).abs() < 1e-9);
    }

    #[test]
    fn sequencing_error_is_not_flagged() {
        // No shared haplotype, just ~9 random errors per read → no candidate columns.
        let rs = reads(60, 2000, 0, &[], 9);
        let e = detect_hidden_copy(&rs, &p());
        assert!(!e.flagged);
        assert_eq!(e.n_alt_positions, 0);
    }

    #[test]
    fn heterozygous_snp_is_not_flagged() {
        // 30 of 60 reads share just 2 alt positions (a diploid het) → below min_alt_positions.
        let rs = reads(60, 2000, 30, &[500, 900], 3);
        let e = detect_hidden_copy(&rs, &p());
        assert!(!e.flagged, "a het (2 positions) must not be called a hidden copy");
        assert!(e.n_alt_positions < p().min_alt_positions);
    }

    #[test]
    fn fixed_difference_above_band_is_not_a_candidate() {
        // A position alt in ~all reads (fixed diff / reference error) is above balanced_hi.
        let rs = reads(40, 2000, 40, &[700], 0); // all 40 alt at 700 → frac 1.0 > 0.60
        let e = detect_hidden_copy(&rs, &p());
        assert_eq!(e.n_alt_positions, 0);
        assert!(!e.flagged);
    }

    #[test]
    fn low_depth_abstains() {
        let shared: Vec<u64> = (100..120).collect();
        let rs = reads(4, 2000, 2, &shared, 0); // below min_depth
        assert!(!detect_hidden_copy(&rs, &p()).flagged);
    }
}
