//! End-artifact flags for a final assembly (register 1500): reverse-transcriptase template switching (RTS) at a junction and
//! intrapriming (genomic A-run after the 3' end). Ports of SQANTI3 6.1 `utilities/rt_switching.py` as `sqanti3_qc.py` calls it
//! (min_match 8, `-a` = one mismatch allowed, wiggle 1, all junctions) and of the `perc_A_downstream_TTS` window of
//! `classification_classifiers.py` (20 bp, strand-aware), so SQANTI3 is a reference for the rule, not a dependency.
//! Rule E (PREREG docs/PREREG_end_artifacts_2026-11-02.md): flag a model when [a junction is RTS OR the A share > 60%] AND its
//! `reads` <= 2.

use crate::genome::IndexedFasta;

/// rt_switching.py PATSEQLEN: sequence searched per side, without the wiggle.
const PATSEQLEN: i64 = 10;
pub const RTS_MIN_MATCH: usize = 8;
pub const RTS_WIGGLE: i64 = 1;
pub const POLYA_WINDOW: i64 = 20;
pub const POLYA_MAX_PCT: f64 = 60.0;
pub const MAX_READS: u64 = 2;

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            x => x,
        })
        .collect()
}

fn upper(v: Vec<u8>) -> Vec<u8> {
    v.into_iter().map(|b| b.to_ascii_uppercase()).collect()
}

/// rt_switching.seq_match: equal length and identical, or (allow_mismatch) at most one mismatch.
fn seq_match(a: &[u8], b: &[u8], allow_mismatch: bool) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a == b || (allow_mismatch && a.iter().zip(b).filter(|(x, y)| x != y).count() <= 1)
}

/// rt_switching.checkForRepeatPat: a direct repeat of `min_match` bases (seed min_match/2 exact, extended exactly, then
/// completed right or left with <= 1 mismatch) between the exon-side and intron-side sequences.
pub fn has_repeat(exon: &[u8], intron: &[u8], min_match: usize, allow_mismatch: bool) -> bool {
    let seed = min_match / 2;
    let n = exon.len();
    if n < seed {
        return false;
    }
    for i in 0..=(n - seed) {
        let s = &exon[i..i + seed];
        let mut offset = 0;
        while offset + seed <= intron.len() {
            let Some(p) = intron[offset..].windows(seed).position(|w| w == s) else {
                break;
            };
            let j = offset + p;
            offset = j + 1;
            let mut k = seed;
            while i + k < n && j + k < intron.len() && exon[i + k] == intron[j + k] {
                k += 1;
            }
            // m may be negative in Python (k > min_match): slices then run backwards and are empty vs non-empty
            let m = min_match as i64 - k as i64;
            let (ii, jj, kk) = (i as i64, j as i64, k as i64);
            if ii + kk + m <= n as i64 && jj + kk + m <= intron.len() as i64 {
                let (a0, a1) = ((ii + kk) as usize, (ii + kk + m).max(ii + kk) as usize);
                let (b0, b1) = ((jj + kk) as usize, (jj + kk + m).max(jj + kk) as usize);
                if seq_match(&exon[a0..a1.min(n)], &intron[b0..b1.min(intron.len())], allow_mismatch) {
                    return true;
                }
            }
            if ii - m >= 0 && jj - m >= 0 {
                let (a0, b0) = ((ii - m) as usize, (jj - m) as usize);
                let a = if a0 <= i { &exon[a0..i] } else { &exon[0..0] };
                let b = if b0 <= j { &intron[b0..j] } else { &intron[0..0] };
                if seq_match(a, b, allow_mismatch) {
                    return true;
                }
            }
        }
    }
    false
}

/// Python slice [s, e) of a contig: negative s counts from the end (only reachable at contig starts; treated as empty).
fn slice(fa: &IndexedFasta, chrom: &str, s: i64, e: i64) -> Vec<u8> {
    if s < 0 {
        return Vec::new();
    }
    upper(fa.fetch(chrom, s, e).unwrap_or_default())
}

/// RTS at one junction; `strpos` = first intron base, `endpos` = last intron base, both 1-based (SQANTI junctions.txt).
pub fn rts_junction(fa: &IndexedFasta, chrom: &str, strand: char, strpos: i64, endpos: i64) -> bool {
    let w = RTS_WIGGLE;
    let cnt = PATSEQLEN + 2 * w;
    let (exon, intron) = if strand == '+' {
        let s = strpos - cnt + w - 1;
        let e = endpos - cnt + w;
        (slice(fa, chrom, s, s + cnt), slice(fa, chrom, e, e + cnt))
    } else {
        let e = strpos - w - 1;
        let s = endpos - w;
        (revcomp(&slice(fa, chrom, s, s + cnt)), revcomp(&slice(fa, chrom, e, e + cnt)))
    };
    !exon.is_empty() && !intron.is_empty() && has_repeat(&exon, &intron, RTS_MIN_MATCH, true)
}

/// Percent A of the POLYA_WINDOW genomic bases after the 3' end (transcript orientation); `start0`/`end` = 0-based start and
/// exclusive end of the model.
pub fn a_share(fa: &IndexedFasta, chrom: &str, strand: char, start0: i64, end: i64) -> f64 {
    let seq = if strand == '+' {
        slice(fa, chrom, end, end + POLYA_WINDOW)
    } else {
        revcomp(&slice(fa, chrom, start0 - POLYA_WINDOW, start0))
    };
    seq.iter().filter(|&&b| b == b'A').count() as f64 / POLYA_WINDOW as f64 * 100.0
}

/// One model's flags: (any RTS junction, A share).
pub fn model_flags(fa: &IndexedFasta, chrom: &str, strand: char, exons: &[(i64, i64)]) -> (bool, f64) {
    // exons: 0-based half-open, sorted
    let rts = exons
        .windows(2)
        .any(|w| rts_junction(fa, chrom, strand, w[0].1 + 1, w[1].0));
    let a = a_share(fa, chrom, strand, exons[0].0, exons[exons.len() - 1].1);
    (rts, a)
}

/// Rule E.
pub fn rule_e(rts: bool, a_pct: f64, reads: u64) -> bool {
    (rts || a_pct > POLYA_MAX_PCT) && reads <= MAX_READS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repeat_of_eight_exact_is_found() {
        assert!(has_repeat(b"TTTTACGTACGTCC", b"GGGACGTACGTGGG", 8, false));
    }

    #[test]
    fn one_mismatch_only_when_allowed() {
        // seed ACGT + exact extension 6, completed by 2 bases with one mismatch
        let (e, i) = (b"CCACGTACGAACC".as_slice(), b"GGACGTACGTAGG".as_slice());
        assert!(!has_repeat(e, i, 8, false));
        assert!(has_repeat(e, i, 8, true));
    }

    #[test]
    fn unrelated_sequences_have_no_repeat() {
        assert!(!has_repeat(b"AAAAAAAAAAAA", b"CCCCCCCCCCCC", 8, true));
    }

    #[test]
    fn rule_e_needs_low_reads() {
        assert!(rule_e(true, 0.0, 2));
        assert!(rule_e(false, 65.0, 1));
        assert!(!rule_e(false, 60.0, 1));
        assert!(!rule_e(true, 90.0, 3));
    }
}
