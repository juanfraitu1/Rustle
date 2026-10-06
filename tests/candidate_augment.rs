//! Integration tests for `utilities candidate-augment`, the Rust port of the retired Python script
//! (kept for provenance at `tools/candidate_augment.py.legacy`).

use std::process::{Command, Output};
use tempfile::TempDir;

const BIN: &str = env!("CARGO_BIN_EXE_utilities");

fn write_fixture(dir: &TempDir) {
    std::fs::write(
        dir.path().join("genome.fa"),
        ">chr1\nACGTACGTACGT\n>chr2\nTGCATGCATGCA\n",
    )
    .unwrap();
    std::fs::write(
        dir.path().join("copies.tsv"),
        "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\tsource\tgene_id\tcore_hull\tsd_depth\tcore_bp\trep_frac\tmember_status\tlocus_start\tlocus_end\n\
         F1\t0\tchr1_A\tchr1\t100\t200\t1\t+\t10\t100-200\t0.99\tfixture\tgA\tNA\t1\t100\t0.000\tkept\t100\t200\n\
         F1\t1\tchr1_B\tchr1\t500\t600\t1\t+\t8\t500-600\t0.98\tfixture\tgA\tNA\t1\t100\t0.000\tkept\t500\t600\n\
         F2\t0\tchr2_C\tchr2\t1000\t1100\t1\t+\t5\t1000-1100\t0.97\tfixture\tgB\tNA\t1\t100\t0.000\tkept\t1000\t1100\n",
    )
    .unwrap();
    std::fs::write(
        dir.path().join("copies.fa"),
        ">F1|0|chr1_A:100-200|+|nexon=1\nACGTACGTACGTACGTACGT\n\
         >F1|1|chr1_B:500-600|+|nexon=1\nTGCATGCATGCATGCATGCA\n\
         >F2|0|chr2_C:1000-1100|+|nexon=1\nAAAATTTTCCCCGGGGAAAA\n",
    )
    .unwrap();
    std::fs::write(
        dir.path().join("regions.txt"),
        "F1\tchr1:90-210\nF1\tchr1:490-610\nF2\tchr2:990-1110\n",
    )
    .unwrap();
    std::fs::write(
        dir.path().join("cand.candidates.tsv"),
        "family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used\n\
         F1\tcand_F1_0\t1\t6\t1\t30\tchr1:300-330\t0.01000\t6\t6\n\
         F2\tcand_F2_0\t1\t7\t1\t25\tchr2:2000-2025\t0.02000\t7\t7\n\
         F1\tcand_F1_1\t1\t4\t0\t28\tchr1:400-428\t0.03000\t4\t4\n",
    )
    .unwrap();
    std::fs::write(
        dir.path().join("cand.contigs.fa"),
        ">cand_F1_0\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\
         >cand_F2_0\nCCCCCCCCCCCCCCCCCCCCCCCCC\n\
         >cand_F1_1\nGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n",
    )
    .unwrap();
}

fn run(dir: &TempDir, extra_args: &[&str]) -> Output {
    let mut c = Command::new(BIN);
    c.arg("candidate-augment")
        .arg("--fasta")
        .arg(dir.path().join("genome.fa"))
        .arg("--copies")
        .arg(dir.path().join("copies.tsv"))
        .arg("--copies-fa")
        .arg(dir.path().join("copies.fa"))
        .arg("--regions")
        .arg(dir.path().join("regions.txt"))
        .arg("--cand")
        .arg(dir.path().join("cand"))
        .arg("--out")
        .arg(dir.path().join("out.aug"));
    c.args(extra_args).output().unwrap()
}

#[test]
fn augments_all_products_for_flagged_candidates() {
    let dir = TempDir::new().unwrap();
    write_fixture(&dir);
    let o = run(&dir, &[]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    let fa = std::fs::read_to_string(dir.path().join("out.aug.fa")).unwrap();
    assert!(fa.starts_with(">chr1\nACGTACGTACGT\n>chr2\nTGCATGCATGCA\n>cand_F1_0\n"));
    assert!(fa.contains(">cand_F2_0\nCCCCCCCCCCCCCCCCCCCCCCCCC\n"));

    let copies = std::fs::read_to_string(dir.path().join("out.aug.copies.tsv")).unwrap();
    assert!(copies.contains("F1\t2\tcand_F1_0\tcand_F1_0\t0\t30"));
    assert!(copies.contains("F2\t1\tcand_F2_0\tcand_F2_0\t0\t25"));
    assert!(!copies.contains("cand_F1_1"));

    let copies_fa = std::fs::read_to_string(dir.path().join("out.aug.copies.fa")).unwrap();
    assert!(copies_fa.contains(">F1|2|cand_F1_0:0-30|+|nexon=1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"));
    assert!(copies_fa.contains(">F2|1|cand_F2_0:0-25|+|nexon=1\nCCCCCCCCCCCCCCCCCCCCCCCCC\n"));

    let regions = std::fs::read_to_string(dir.path().join("out.aug.regions.txt")).unwrap();
    assert_eq!(
        regions.lines().collect::<Vec<_>>(),
        vec![
            "chr1:90-210",
            "chr1:490-610",
            "chr2:990-1110",
            "cand_F1_0:0-30",
            "cand_F2_0:0-25"
        ]
    );

    let families = std::fs::read_to_string(dir.path().join("out.aug.families.txt")).unwrap();
    assert_eq!(families, "F1\nF2\n");

    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(stderr.contains("candidate_augment: 2 candidate contig(s) of 2 families"));
}

#[test]
fn respects_existing_fai_for_name_collision_check() {
    let dir = TempDir::new().unwrap();
    write_fixture(&dir);
    // samtools faidx writes a .fai. The binary should use it for genome names.
    let status = Command::new("samtools")
        .args(["faidx", dir.path().join("genome.fa").to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let o = run(&dir, &[]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
}

#[test]
fn refuses_genome_name_collision_and_exits_2() {
    let dir = TempDir::new().unwrap();
    write_fixture(&dir);
    std::fs::write(
        dir.path().join("cand.candidates.tsv"),
        "family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used\n\
         F1\tchr1\t1\t6\t1\t30\tchr1:300-330\t0.01000\t6\t6\n",
    )
    .unwrap();
    std::fs::write(
        dir.path().join("cand.contigs.fa"),
        ">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
    )
    .unwrap();
    let o = run(&dir, &[]);
    assert!(!o.status.success());
    assert_eq!(o.status.code(), Some(2));
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("candidate chr1 already names a sequence"));
    // Nothing should be written.
    for suffix in [
        "fa",
        "copies.tsv",
        "copies.fa",
        "regions.txt",
        "families.txt",
    ] {
        assert!(
            !dir.path().join(format!("out.aug.{suffix}")).exists(),
            "{suffix} exists"
        );
    }
}

#[test]
fn refuses_missing_family_in_copies_table_and_exits_2() {
    let dir = TempDir::new().unwrap();
    write_fixture(&dir);
    std::fs::write(
        dir.path().join("cand.candidates.tsv"),
        "family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used\n\
         F3\tcand_F3_0\t1\t6\t1\t30\tchr1:300-330\t0.01000\t6\t6\n",
    )
    .unwrap();
    std::fs::write(
        dir.path().join("cand.contigs.fa"),
        ">cand_F3_0\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
    )
    .unwrap();
    let o = run(&dir, &[]);
    assert!(!o.status.success());
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("family F3 has no row"));
}

#[test]
fn refuses_missing_region_and_exits_2() {
    let dir = TempDir::new().unwrap();
    write_fixture(&dir);
    std::fs::write(dir.path().join("regions.txt"), "F1\tchr1:90-210\n").unwrap();
    let o = run(&dir, &[]);
    assert!(!o.status.success());
    assert_eq!(o.status.code(), Some(2));
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("has no region for the candidate family F2")
    );
}

#[test]
fn refuses_no_flagged_candidates_and_exits_2() {
    let dir = TempDir::new().unwrap();
    write_fixture(&dir);
    std::fs::write(
        dir.path().join("cand.candidates.tsv"),
        "family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used\n\
         F1\tcand_F1_0\t1\t6\t0\t30\tchr1:300-330\t0.01000\t6\t6\n",
    )
    .unwrap();
    let o = run(&dir, &[]);
    assert!(!o.status.success());
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("no flagged candidate"));
}

#[test]
fn refuses_gzipped_genome_and_exits_2() {
    let dir = TempDir::new().unwrap();
    write_fixture(&dir);
    let gz = dir.path().join("genome.fa.gz");
    let status = Command::new("gzip")
        .args(["-c", dir.path().join("genome.fa").to_str().unwrap()])
        .stdout(std::fs::File::create(&gz).unwrap())
        .status()
        .unwrap();
    assert!(status.success());
    let mut c = Command::new(BIN);
    let o = c
        .arg("candidate-augment")
        .arg("--fasta")
        .arg(&gz)
        .arg("--copies")
        .arg(dir.path().join("copies.tsv"))
        .arg("--copies-fa")
        .arg(dir.path().join("copies.fa"))
        .arg("--regions")
        .arg(dir.path().join("regions.txt"))
        .arg("--cand")
        .arg(dir.path().join("cand"))
        .arg("--out")
        .arg(dir.path().join("out.aug"))
        .output()
        .unwrap();
    assert!(!o.status.success());
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("is compressed"));
}
