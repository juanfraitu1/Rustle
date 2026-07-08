use std::process::Command;

#[test]
fn same_chrom_supplement_recovers_same_chrom_family() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out";
    let status = Command::new(bin)
        .args([
            "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
            "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
            "--out", out,
            "--cross-chrom",
            "--same-chrom-supplement-win", "1000",
        ])
        .status()
        .expect("gw_family_catalog failed to run");
    assert!(status.success());

    let families = std::fs::read_to_string(format!("{}.families.tsv", out)).unwrap();
    let copies = std::fs::read_to_string(format!("{}.copies.tsv", out)).unwrap();
    let sidecar = std::fs::read_to_string(format!("{}.same_chrom_supplement.tsv", out)).unwrap();

    // With the fixture we expect at least one family from the supplement.
    let supp_count = sidecar.lines().count().saturating_sub(1);
    assert!(supp_count >= 1, "expected at least one same-chrom supplement family, got {}\n{}", supp_count, sidecar);

    // The main catalog should contain the cross-chrom family as well.
    assert!(families.lines().count() > 1, "expected at least one family\n{}", families);
    assert!(copies.lines().count() > 1, "expected at least one copy\n{}", copies);
}
