//! Failure-mode classifier for VG-HMM novel-copy rescue.
//!
//! ## Two-stage pipeline
//!
//! **Stage 1 – Internal** (`classify_internal`): fast, always runs on every
//! rescued read using information already in memory.
//!
//! **Stage 2 – External** (`classify_external`): slower minimap2-based
//! verification; only runs when the internal stage returns
//! `NeedsExternalVerification` **and** `config.vg_rescue_diagnostic == true`.
//!
//! ## Buckets
//!
//! | Variant                   | Meaning                                                  |
//! |---------------------------|----------------------------------------------------------|
//! | `BelowThresholdChain`     | Too few k-mer hits → chain score below 40.0              |
//! | `SeedMasked`              | >30 % of seed positions land in repeat-masked regions    |
//! | `Divergent`               | Maps with minimap2 but has no long indels (SNV/small)    |
//! | `Structural`              | Maps with minimap2 and has ≥1 indel ≥ 50 bp in CIGAR    |
//! | `ReferenceAbsent`         | Cannot be aligned even with maximally-relaxed minimap2   |
//! | `NeedsExternalVerification` | Internal stage inconclusive; requires minimap2 check   |

use anyhow::Result;
use std::path::Path;

// ── Enum ─────────────────────────────────────────────────────────────────────

/// Classification of why a read was or was not rescued.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RescueClass {
    /// Chain score below threshold (too few k-mer hits).
    BelowThresholdChain,
    /// More than 30 % of seed windows fall in repeat-masked (N) positions.
    SeedMasked,
    /// Aligned by minimap2 but CIGAR has no indel ≥ 50 bp.
    Divergent,
    /// Aligned by minimap2 and CIGAR contains at least one indel ≥ 50 bp.
    Structural,
    /// No primary alignment found even with maximally-relaxed minimap2 parameters.
    ReferenceAbsent,
    /// Internal stage inconclusive; requires external minimap2 verification.
    NeedsExternalVerification,
}

impl RescueClass {
    /// Short snake_case label suitable for GTF attributes and TSV columns.
    pub fn to_str(self) -> &'static str {
        match self {
            RescueClass::BelowThresholdChain => "below_threshold_chain",
            RescueClass::SeedMasked => "seed_masked",
            RescueClass::Divergent => "divergent",
            RescueClass::Structural => "structural",
            RescueClass::ReferenceAbsent => "reference_absent",
            RescueClass::NeedsExternalVerification => "needs_external_verification",
        }
    }
}

impl std::fmt::Display for RescueClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.to_str())
    }
}

// ── Stage 1: Internal classifier ─────────────────────────────────────────────

/// Fast internal failure-mode classification using in-memory statistics.
///
/// # Parameters
/// - `n_kmer_hits`: number of distinct k-mer hits in the read against the
///   family k-mer set (from `prefilter_read`).
/// - `k`: k-mer length used (typically 15).
/// - `read_len`: length of the read in bases.
/// - `masked_fraction`: fraction of the read's k-mer windows that overlap
///   at least one `N` base (proxy for repeat-masked regions).
///
/// # Decision tree
/// 1. `approx_chain_score = n_kmer_hits * k`; if < 40.0 → `BelowThresholdChain`.
/// 2. If `masked_fraction > 0.30` → `SeedMasked`.
/// 3. Otherwise → `NeedsExternalVerification`.
pub fn classify_internal(
    n_kmer_hits: usize,
    k: usize,
    _read_len: usize,
    masked_fraction: f64,
) -> RescueClass {
    let approx_chain_score = n_kmer_hits as f64 * k as f64;
    if approx_chain_score < 40.0 {
        return RescueClass::BelowThresholdChain;
    }
    if masked_fraction > 0.30 {
        return RescueClass::SeedMasked;
    }
    RescueClass::NeedsExternalVerification
}

// ── Stage 2: External minimap2 classifier ────────────────────────────────────

/// Parse a CIGAR string for any `I`, `D`, or `N` operation of length ≥
/// `threshold`.
///
/// Returns `true` if at least one such operation is found.
pub fn cigar_has_long_indel(cigar: &str, threshold: usize) -> bool {
    let mut num_buf = 0usize;
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_buf = num_buf * 10 + (ch as usize - '0' as usize);
        } else {
            let op = ch;
            let len = num_buf;
            num_buf = 0;
            if matches!(op, 'I' | 'D' | 'N') && len >= threshold {
                return true;
            }
        }
    }
    false
}

/// Determine whether a SAM output line (non-header) represents a primary
/// alignment.  A primary alignment has flag bit 4 (unmapped) clear AND flag
/// bit 256 (secondary) clear AND flag bit 2048 (supplementary) clear.
fn is_primary_alignment(sam_line: &str) -> bool {
    let mut fields = sam_line.splitn(12, '\t');
    let _qname = fields.next();
    let flag_str = match fields.next() {
        Some(f) => f,
        None => return false,
    };
    let flag: u16 = match flag_str.parse() {
        Ok(f) => f,
        Err(_) => return false,
    };
    // bit 4 = unmapped, bit 256 = secondary, bit 2048 = supplementary
    (flag & 4) == 0 && (flag & 256) == 0 && (flag & 2048) == 0
}

/// Extract the CIGAR string (field index 5) from a SAM line.
fn extract_cigar(sam_line: &str) -> Option<&str> {
    sam_line.splitn(12, '\t').nth(5)
}

/// Try to align `read_seq` to `ref_fasta` using minimap2 with `extra_args`.
///
/// Preferred minimap2 binary path (checked first; falls back to PATH lookup).
const MINIMAP2_PREFERRED_PATH: &str = "/home/juanfra/miniforge3/bin/minimap2";

/// Resolve the minimap2 command: use `MINIMAP2_PREFERRED_PATH` if it exists,
/// otherwise fall back to `"minimap2"` (PATH lookup).
fn minimap2_cmd() -> std::process::Command {
    if std::path::Path::new(MINIMAP2_PREFERRED_PATH).exists() {
        std::process::Command::new(MINIMAP2_PREFERRED_PATH)
    } else {
        std::process::Command::new("minimap2")
    }
}

/// Returns `(has_primary, primary_cigar)`.
fn try_minimap2(
    ref_fasta: &Path,
    read_seq: &[u8],
    extra_args: &[&str],
) -> Result<(bool, Option<String>)> {
    use std::io::Write;
    use std::process::Stdio;

    let mut child = minimap2_cmd()
        .args(extra_args)
        .arg("-a") // SAM output
        .arg(ref_fasta)
        .arg("-") // read from stdin
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()?;

    // Write a single FASTA record to stdin.
    if let Some(stdin) = child.stdin.as_mut() {
        stdin.write_all(b">query\n")?;
        stdin.write_all(read_seq)?;
        stdin.write_all(b"\n")?;
    }

    let output = child.wait_with_output()?;
    let stdout = String::from_utf8_lossy(&output.stdout);

    for line in stdout.lines() {
        if line.starts_with('@') {
            continue;
        }
        if is_primary_alignment(line) {
            let cigar = extract_cigar(line).map(|s| s.to_owned());
            return Ok((true, cigar));
        }
    }
    Ok((false, None))
}

/// External minimap2-based classification.
///
/// Tries progressively-relaxed minimap2 invocations until either a primary
/// alignment is found or all attempts are exhausted.
///
/// When a primary alignment is found, the CIGAR is inspected:
/// - Any `I`/`D`/`N` operation ≥ 50 bp → `Structural`.
/// - Otherwise → `Divergent`.
///
/// If no primary alignment is found after all attempts → `ReferenceAbsent`.
///
/// # Minimap2 parameter ladder
/// 1. Default (long-read ONT preset `-x map-ont`)
/// 2. `-x map-ont -s 20`
/// 3. `-x map-ont -f 0.0001`
/// 4. `-x map-ont -N 50 --secondary=yes`
/// 5. `-s 10 -f 0 -N 100` (most permissive, no preset)
pub fn classify_external(ref_fasta: &Path, read_seq: &[u8]) -> Result<RescueClass> {
    let ladders: &[&[&str]] = &[
        &["-x", "map-ont"],
        &["-x", "map-ont", "-s", "20"],
        &["-x", "map-ont", "-f", "0.0001"],
        &["-x", "map-ont", "-N", "50", "--secondary=yes"],
        &["-s", "10", "-f", "0", "-N", "100"],
    ];

    for args in ladders {
        match try_minimap2(ref_fasta, read_seq, args) {
            Ok((true, cigar_opt)) => {
                let cigar = cigar_opt.as_deref().unwrap_or("*");
                if cigar_has_long_indel(cigar, 50) {
                    return Ok(RescueClass::Structural);
                } else {
                    return Ok(RescueClass::Divergent);
                }
            }
            Ok((false, _)) => continue,
            Err(_) => continue,
        }
    }

    Ok(RescueClass::ReferenceAbsent)
}
