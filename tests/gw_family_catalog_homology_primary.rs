use std::process::Command;
#[test]
fn homology_primary_emits_families_on_fixture() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out_hom";
    let status = Command::new(bin).args([
        "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
        "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
        "--out", out, "--homology-primary",
    ]).status().unwrap();
    assert!(status.success());
    let fams = std::fs::read_to_string(format!("{}.families.tsv", out)).unwrap();
    assert!(fams.lines().count() > 1, "expected >=1 homology family\n{}", fams);
}

#[test]
fn homology_primary_writes_famcn_when_enumerating() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out_famcn";
    let status = Command::new(bin).args([
        "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
        "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
        "--out", out, "--homology-primary", "--enumerate-copies",
    ]).status().unwrap();
    assert!(status.success());
    assert!(std::fs::metadata(format!("{}.famcn.tsv", out)).is_ok(), "famcn.tsv must be written");
}
