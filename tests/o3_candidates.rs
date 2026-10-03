use std::process::Command;
#[test]
fn flags_the_deleted_copy_and_assign_places_its_reads() {
    let dir = tempfile::tempdir().unwrap();
    let fx = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/o3_candidates");
    let mmi = dir.path().join("genome.splice.mmi");
    assert!(Command::new(std::env::var("RUSTLE_MINIMAP2").unwrap_or("minimap2".into())).args(["-x", "splice", "-d"]).arg(&mmi).arg(format!("{fx}/genome.fa")).status().unwrap().success());
    let out = dir.path().join("t.cand");
    let o = Command::new(env!("CARGO_BIN_EXE_o3_candidates")).args(["--bam", &format!("{fx}/reads.bam"), "--fasta", &format!("{fx}/genome.fa"), "--copies", &format!("{fx}/copies.tsv"), "--copies-fa", &format!("{fx}/copies.fa"), "--index", mmi.to_str().unwrap(), "--out", out.to_str().unwrap(), "--threads", "2"]).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let cands = std::fs::read_to_string(format!("{}.candidates.tsv", out.display())).unwrap();
    let flagged: Vec<&str> = cands.lines().skip(1).filter(|l| l.split('\t').nth(4) == Some("1")).collect();
    assert_eq!(flagged.len(), 1, "{cands}");
    let fa = std::fs::read_to_string(format!("{}.contigs.fa", out.display())).unwrap();
    assert!(fa.starts_with(">cand_MCL0_0\n"));
    let union_len = fa.lines().nth(1).unwrap().len();
    assert!((850..=950).contains(&union_len), "union {union_len} bp, expected the 900-bp spliced copy");
}
