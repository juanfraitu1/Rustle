use std::process::Command;

#[test]
fn default_cross_chrom_output_is_unchanged() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out_default";
    let status = Command::new(bin)
        .args([
            "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
            "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
            "--out", out,
            "--cross-chrom",
        ])
        .status()
        .expect("gw_family_catalog failed to run");
    assert!(status.success());

    let families = std::fs::read_to_string(format!("{}.families.tsv", out)).unwrap();
    let copies = std::fs::read_to_string(format!("{}.copies.tsv", out)).unwrap();

    // Without the supplement flag, no sidecar should be written and the catalog
    // should only contain the cross-chrom family.
    assert!(
        std::fs::metadata(format!("{}.same_chrom_supplement.tsv", out)).is_err(),
        "sidecar must not be created when flag is absent"
    );
    assert!(families.lines().count() > 1);
    assert!(copies.lines().count() > 1);
}
