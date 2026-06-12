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

#[test]
fn family_ids_unique_across_strands() {
    use rustle::annotation_families::{families_from_grouped, CopyStructure};
    // Reuse two similar copies per strand so each strand forms a family.
    let plus_a = CopyStructure{copy_id:"PA".into(),chrom:"chr1".into(),strand:'+',exons:vec![(0,80)]};
    let plus_b = CopyStructure{copy_id:"PB".into(),chrom:"chr2".into(),strand:'+',exons:vec![(0,80)]};
    let minus_a = CopyStructure{copy_id:"MA".into(),chrom:"chr3".into(),strand:'-',exons:vec![(0,80)]};
    let minus_b = CopyStructure{copy_id:"MB".into(),chrom:"chr4".into(),strand:'-',exons:vec![(0,80)]};
    // Identical sequence within each strand -> each strand merges into one family.
    let s = SEQ_A.to_vec(); // reuse the non-repetitive 82bp SEQ_A constant from earlier cluster tests
    let fams = families_from_grouped(vec![
        (plus_a, s.clone()), (plus_b, s.clone()), (minus_a, s.clone()), (minus_b, s.clone())], 0.5);
    assert_eq!(fams.len(), 2, "one family per strand");
    let ids: std::collections::HashSet<_> = fams.iter().map(|f| f.family_id.clone()).collect();
    assert_eq!(ids.len(), 2, "family_ids must be globally unique across strands");
}

#[test]
fn vg_families_report_shows_per_copy_chroms() {
    use rustle::annotation_families::{AnnotationFamily, CopyStructure};
    let fam = AnnotationFamily {
        family_id: "FAM0".into(),
        copies: vec![
            CopyStructure{copy_id:"A".into(),chrom:"chr1".into(),strand:'+',exons:vec![(0,30)]},
            CopyStructure{copy_id:"B".into(),chrom:"chr5".into(),strand:'+',exons:vec![(0,30)]},
        ],
        achieved_min_sim: 0.93, achieved_mean_sim: 0.93,
    };
    let tsv = rustle::vg_eval::render_vg_families(&[fam], 0.5);
    assert!(tsv.contains("chr1;chr5"), "copy_chroms must list both contigs; got:\n{}", tsv);
    assert!(tsv.contains("0.500"), "must record merge_threshold_T");
    assert!(tsv.contains("0.930"), "must record achieved similarity");
}

#[test]
fn scorecard_verdicts() {
    use rustle::vg_eval::verdict;
    assert_eq!(verdict(true, false), "WIN");        // in_vg, not in stringtie
    assert_eq!(verdict(true, true), "tie");
    assert_eq!(verdict(false, false), "miss");
    assert_eq!(verdict(false, true), "regression"); // stringtie had it, vg lost it
}

#[test]
fn chain_match_within_tolerance() {
    use rustle::vg_eval::chain_in_set;
    let mut set = std::collections::HashSet::new();
    set.insert(vec![(100u64, 200u64), (300, 400)]);
    // exact match
    assert!(chain_in_set(&[(100,200),(300,400)], &set));
    // within +/-2 bp per coordinate
    assert!(chain_in_set(&[(101,199),(302,400)], &set));
    // off by >2 -> no match
    assert!(!chain_in_set(&[(105,200),(300,400)], &set));
    // different length -> no match
    assert!(!chain_in_set(&[(100,200)], &set));
}

#[test]
fn render_vg_eval_summary_counts() {
    use rustle::vg_eval::render_vg_eval;
    let entries = vec![
        ("c1".to_string(), "FAM0".to_string(), true, true),   // in_st, in_vg -> tie
        ("c2".to_string(), "FAM0".to_string(), false, true),  // !in_st, in_vg -> WIN
        ("c3".to_string(), "FAM1".to_string(), true, false),  // in_st, !in_vg -> regression
        ("c4".to_string(), "FAM1".to_string(), false, false), // miss
    ];
    let tsv = render_vg_eval(&entries);
    assert!(tsv.contains("c2\tFAM0\tyes\tfalse\ttrue\tWIN"), "WIN row; got:\n{}", tsv);
    assert!(tsv.contains("WIN=1"));
    assert!(tsv.contains("tie=1"));
    assert!(tsv.contains("miss=1"));
    assert!(tsv.contains("regression=1"));
    assert!(tsv.contains("mode=guided"));
}
