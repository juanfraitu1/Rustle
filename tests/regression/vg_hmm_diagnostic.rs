//! Regression tests for `vg_hmm::diagnostic` — internal and external classifiers.

use rustle::vg_hmm::diagnostic::{
    cigar_has_long_indel, classify_internal, RescueClass,
};

// ── Task 6.1: Internal classifier ────────────────────────────────────────────

/// Many short islands (small n_kmer_hits, k=15) → approx_chain_score < 40 →
/// BelowThresholdChain.
#[test]
fn test_internal_below_threshold_chain() {
    // 2 hits × 15 = 30.0 < 40.0  → BelowThresholdChain
    let class = classify_internal(2, 15, 1000, 0.0);
    assert_eq!(
        class,
        RescueClass::BelowThresholdChain,
        "expected BelowThresholdChain, got {:?}",
        class
    );
    // Boundary: exactly 2 hits → still below
    let class2 = classify_internal(2, 15, 500, 0.05);
    assert_eq!(class2, RescueClass::BelowThresholdChain);
}

/// High hit count but high masked fraction → SeedMasked.
#[test]
fn test_internal_seed_masked() {
    // 10 hits × 15 = 150.0 ≥ 40.0; masked_fraction = 0.50 > 0.30 → SeedMasked
    let class = classify_internal(10, 15, 2000, 0.50);
    assert_eq!(
        class,
        RescueClass::SeedMasked,
        "expected SeedMasked, got {:?}",
        class
    );
}

/// High hits and low mask → NeedsExternalVerification.
#[test]
fn test_internal_needs_external() {
    // 5 hits × 15 = 75.0 ≥ 40.0; masked_fraction = 0.10 ≤ 0.30
    let class = classify_internal(5, 15, 1500, 0.10);
    assert_eq!(
        class,
        RescueClass::NeedsExternalVerification,
        "expected NeedsExternalVerification, got {:?}",
        class
    );
}

// ── CIGAR helper ─────────────────────────────────────────────────────────────

#[test]
fn test_cigar_long_indel_detection() {
    assert!(cigar_has_long_indel("10M50I5M", 50));
    assert!(cigar_has_long_indel("5M100D3M", 50));
    assert!(cigar_has_long_indel("2M200N8M", 50));
    assert!(!cigar_has_long_indel("10M5I2M3D", 50));
    assert!(!cigar_has_long_indel("20M", 50));
    assert!(!cigar_has_long_indel("*", 50));
    // Exactly at threshold is a match.
    assert!(cigar_has_long_indel("5M50I", 50));
    // One below threshold is not.
    assert!(!cigar_has_long_indel("5M49I", 50));
}

// ── Task 6.2: External classifier (skip-if-missing) ──────────────────────────

#[test]
fn test_external_classify_skip_if_no_minimap2() {
    use std::process::Command;
    use std::io::Write;

    if Command::new("minimap2").arg("--version").output().is_err() {
        eprintln!("skipping: minimap2 not on PATH");
        return;
    }

    // Build a tiny reference FASTA and a maximally divergent read in a temp dir.
    let dir = tempfile::tempdir().expect("tempdir");
    let ref_path = dir.path().join("ref.fa");
    let mut ref_file = std::fs::File::create(&ref_path).unwrap();
    // 200 bp of A/T reference.
    let ref_seq = b"AAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTT";
    writeln!(ref_file, ">ref").unwrap();
    ref_file.write_all(ref_seq).unwrap();
    writeln!(ref_file).unwrap();
    drop(ref_file);

    // A random GC-rich read that is maximally divergent from the AT reference.
    let read_seq = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";

    let result = rustle::vg_hmm::diagnostic::classify_external(&ref_path, read_seq);
    let class = result.expect("classify_external should not error");
    assert!(
        matches!(class, RescueClass::Divergent | RescueClass::ReferenceAbsent),
        "expected Divergent or ReferenceAbsent, got {:?}",
        class
    );
}
