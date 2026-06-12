use rustle::types::RunConfig;

#[test]
fn runconfig_g2_defaults_off() {
    let c = RunConfig::default();
    assert!(c.guide2_path.is_none(), "guide2_path must default to None");
    assert!((c.family_exon_similarity - rustle::vg_family::family_graph::family_merge_jaccard()).abs() < 1e-9,
        "family_exon_similarity default must equal family_merge_jaccard()");
}

use rustle::annotation_families::{load_copies, CopyStructure};

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
