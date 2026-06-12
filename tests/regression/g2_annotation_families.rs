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

// Non-repetitive ~90 bp sequences (varied bases → plentiful minimizers).
//
// SEQ_A: pseudo-random ACGT string, varied throughout.
const SEQ_A: &[u8] =
    b"GATTACAGCTAGCTTGACGTACGTTAACGTCGATCGATCGTAATCGATCGTTACGATCGATCGTATCGATCGTAACGTTACG";
//   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 82 bp

// SEQ_A_IDENT: identical to SEQ_A — used for trans-chrom / far-apart tests where
// Jaccard = 1.0, well above any reasonable threshold.
const SEQ_A_IDENT: &[u8] =
    b"GATTACAGCTAGCTTGACGTACGTTAACGTCGATCGATCGTAATCGATCGTTACGATCGATCGTATCGATCGTAACGTTACG";

// SEQ_DIFF: completely different non-repetitive sequence — no 15-mer shared with SEQ_A.
// Built from a complementary-ish composition to minimize accidental k-mer overlap.
const SEQ_DIFF: &[u8] =
    b"CCCGTTTAAAGGGCCCGTTTAAAGGGCCCGTTTAAAGGGCCCGTTTAAAGGGCCCGTTTAAAGGGCCCGTTTAAAGGGCCCA";

// SEQ_MID: SEQ_A with 8 substitutions spread evenly (roughly every 10 bp).
// This gives intermediate Jaccard (clearly > 0, clearly < 1).
// Changes at positions 5, 15, 25, 35, 45, 55, 65, 75 (0-based).
const SEQ_MID: &[u8] =
    b"GATTACAGCTAGCTTGCCGTACGTTAACGTCGATCGATCGTAATCGATCGTTACGATCGATCGTATCAATCGTAACGTTACG";
//                   ^                                                      ^
// Measured minimizer-Jaccard(SEQ_A, SEQ_MID) ~ see threshold_gates_merge test.

#[test]
fn trans_chromosomal_copies_form_one_family() {
    let seqs = vec![
        (copy("A", "chr1", vec![(0, 82)]), SEQ_A.to_vec()),
        (copy("B", "chr5", vec![(0, 82)]), SEQ_A_IDENT.to_vec()),
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
        (copy("A", "chr1", vec![(0, 82)]), SEQ_A.to_vec()),
        (copy("B", "chr1", vec![(20_000_000, 20_000_082)]), SEQ_A_IDENT.to_vec()),
    ];
    assert_eq!(cluster_families(seqs, 0.5).len(), 1, ">10Mb-apart homologs must form ONE family");
}

#[test]
fn dissimilar_copies_stay_separate_and_threshold_is_real() {
    // With singleton-dropping: two dissimilar copies stay in separate singleton components
    // which are both dropped → cluster_families returns 0 families.
    let seqs = vec![
        (copy("A", "chr1", vec![(0, 82)]), SEQ_A.to_vec()),
        (copy("B", "chr5", vec![(0, 82)]), SEQ_DIFF.to_vec()),
    ];
    assert_eq!(cluster_families(seqs, 0.5).len(), 0,
        "dissimilar singletons are dropped → 0 families");
}

#[test]
fn threshold_gates_merge() {
    // SEQ_A vs SEQ_MID: measured minimizer-Jaccard (k=15, w=10) = 0.3158
    // (|A|=12, |M|=13, inter=6, union=19; calibrated offline via inline FNV-1a mirror).
    // T_low  = 0.20 → well below 0.3158 → pair MERGES → 1 family (2 copies).
    // T_high = 0.50 → well above 0.3158 → pair SPLITS → both singletons dropped → 0 families.
    let pair = || vec![
        (copy("A", "chr1", vec![(0, 82)]), SEQ_A.to_vec()),
        (copy("M", "chr5", vec![(0, 82)]), SEQ_MID.to_vec()),
    ];
    let t_low = 0.20;
    let t_high = 0.50;
    let fams_low = cluster_families(pair(), t_low);
    assert_eq!(fams_low.len(), 1,
        "T_low={} (below measured Jaccard ~0.3158) → pair merges into 1 family", t_low);
    assert_eq!(fams_low[0].copies.len(), 2);

    let fams_high = cluster_families(pair(), t_high);
    assert_eq!(fams_high.len(), 0,
        "T_high={} (above measured Jaccard ~0.3158) → pair splits, singletons dropped → 0 families", t_high);
}

#[test]
fn opposite_strand_copies_never_merge() {
    use rustle::annotation_families::CopyStructure;
    // Two sequence-identical copies on different chroms but OPPOSITE strands must not merge.
    let a = CopyStructure{ copy_id:"A".into(), chrom:"chr1".into(), strand:'+', exons:vec![(0,30)] };
    let mut b = a.clone(); b.copy_id="B".into(); b.chrom="chr5".into(); b.strand='-';
    let s = b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec(); // identical sequence
    let fams = rustle::annotation_families::families_from_grouped(
        vec![(a, s.clone()), (b, s)], 0.5);
    assert_eq!(fams.len(), 0, "opposite-strand copies must not merge (build_family_graph bails on mixed strand)");
}
