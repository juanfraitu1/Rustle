use rustle::types::RunConfig;

#[test]
fn runconfig_g2_defaults_off() {
    let c = RunConfig::default();
    assert!(c.guide2_path.is_none(), "guide2_path must default to None");
    assert!((c.family_exon_similarity - rustle::vg_family::family_graph::family_merge_jaccard()).abs() < 1e-9,
        "family_exon_similarity default must equal family_merge_jaccard()");
}

use rustle::annotation_families::{cluster_families, load_copies, CopyStructure};

#[test]
fn copy_sequence_concatenates_exons() {
    // Build genome chr1 = "ACGTACGTAC" via a temp FASTA file + GenomeIndex::from_fasta.
    let dir = tempfile::tempdir().unwrap();
    let fa = dir.path().join("test.fa");
    std::fs::write(&fa, ">chr1\nACGTACGTAC\n").unwrap();
    let genome = rustle::genome::GenomeIndex::from_fasta(fa.to_str().unwrap()).unwrap();

    let copy = CopyStructure {
        copy_id: "A".into(),
        chrom: "chr1".into(),
        strand: '+',
        exons: vec![(0, 4), (6, 10)],
    };
    // exons[0]: bytes[0..4] = "ACGT"
    // exons[1]: bytes[6..10] = "GTAC"
    // concatenated: "ACGTGTAC"
    let seq = rustle::annotation_families::copy_sequence(&copy, &genome).unwrap();
    assert_eq!(seq, b"ACGTGTAC");
}

#[test]
fn load_copies_parses_two_transcripts() {
    let gtf = "chr1\ttest\texon\t101\t200\t.\t+\t.\ttranscript_id \"A\";\n\
               chr5\ttest\texon\t301\t400\t.\t+\t.\ttranscript_id \"B\";\n";
    let path = std::env::temp_dir().join("g2_load.gtf");
    std::fs::write(&path, gtf).unwrap();
    let copies = load_copies(&path).unwrap();
    assert_eq!(copies.len(), 2);
    assert_eq!(copies[0].copy_id, "A");
    assert_eq!(copies[0].chrom, "chr1");
    assert_eq!(copies[1].chrom, "chr5");
    // 1-based inclusive [101,200] -> 0-based half-open [100,200)
    assert_eq!(copies[0].exons, vec![(100, 200)]);
}

// ── Task 4: cluster_families ──────────────────────────────────────────────────

fn copy(id: &str, chrom: &str, exons: Vec<(u64, u64)>) -> CopyStructure {
    CopyStructure { copy_id: id.into(), chrom: chrom.into(), strand: '+', exons }
}

// 60 bp sequences (2× repeat of the 30 bp motifs) so that k=15, w=10 produces
// enough minimizers for the Jaccard thresholds to behave as asserted.
// With 60 bp: n = 60-15+1 = 46 positions, ~36 windows → rich minimizer sets.

/// Two sequences that are nearly identical (1-base difference at position 59):
/// high Jaccard >= 0.5, used for the trans-chromosomal and far-apart tests.
const SEQ_A: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTAC";
const SEQ_A2: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTAG";

/// Completely different sequence: Jaccard with SEQ_A should be ~0.0.
const SEQ_DIFF: &[u8] = b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAATTTTGGGGCCCCAAAATTTTGGGGCCCCAA";

/// Partially similar to SEQ_A: differs in 15 central bases (positions 22-36).
/// Jaccard expected to be moderate: >= 0.5 (clusters at T=0.5) but < 0.99 (splits at T=0.99).
const SEQ_PART: &[u8] = b"ACGTACGTACGTACGTACGTACTTTTTTTTTTTTTTACGTACGTACGTACGTACGTACGT";

#[test]
fn trans_chromosomal_copies_form_one_family() {
    let seqs = vec![
        (copy("A", "chr1", vec![(0, 60)]), SEQ_A.to_vec()),
        (copy("B", "chr5", vec![(0, 60)]), SEQ_A2.to_vec()),
    ];
    let fams = cluster_families(seqs, 0.5);
    assert_eq!(fams.len(), 1, "trans-chromosomal homologs must form ONE family");
    assert_eq!(fams[0].copies.len(), 2);
    let chroms: Vec<_> = fams[0].copies.iter().map(|c| c.chrom.as_str()).collect();
    assert!(chroms.contains(&"chr1") && chroms.contains(&"chr5"));
}

#[test]
fn far_apart_same_chrom_form_one_family() {
    let seqs = vec![
        (copy("A", "chr1", vec![(0, 60)]), SEQ_A.to_vec()),
        (copy("B", "chr1", vec![(20_000_000, 20_000_060)]), SEQ_A.to_vec()),
    ];
    assert_eq!(cluster_families(seqs, 0.5).len(), 1, ">10Mb-apart homologs must form ONE family");
}

#[test]
fn dissimilar_copies_stay_separate_and_threshold_is_real() {
    let seqs = vec![
        (copy("A", "chr1", vec![(0, 60)]), SEQ_A.to_vec()),
        (copy("B", "chr5", vec![(0, 60)]), SEQ_DIFF.to_vec()),
    ];
    assert_eq!(cluster_families(seqs.clone(), 0.5).len(), 2, "dissimilar -> 2 families");
    let sim = vec![
        (copy("A", "chr1", vec![(0, 60)]), SEQ_A.to_vec()),
        (copy("B", "chr5", vec![(0, 60)]), SEQ_PART.to_vec()),
    ];
    assert_eq!(cluster_families(sim.clone(), 0.99).len(), 2, "T above achieved sim -> split");
}
