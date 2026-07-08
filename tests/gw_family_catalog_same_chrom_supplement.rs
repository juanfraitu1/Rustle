use std::process::Command;

/// Same-chromosome multi-copy families must appear in the `--cross-chrom` O1 catalog directly (via the
/// chrom/strand-agnostic Louvain emission), NOT through a separate lossy `colocated_families` supplement.
/// The fixture carries a same-chrom pair (c1:B + c1:C) and a cross-chrom pair (c1:A + c2:X); plain
/// `--cross-chrom` must emit the same-chrom family. (The removed `--same-chrom-supplement-win` path routed
/// this through `colocated_families`, whose `(chrom,strand)` + `win` filters dropped inverted-duplication
/// and distant same-chrom paralogs.)
#[test]
fn cross_chrom_catalog_emits_same_chrom_family() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out";
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

    // At least one SAME-chromosome multi-copy family (n_chroms == 1, n_copies >= 2) must be present.
    // families.tsv columns: family_id  n_copies  n_chroms  chroms  cross_chrom  avg_reads
    let same_chrom_families: Vec<&str> = families
        .lines()
        .skip(1)
        .filter(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f.len() >= 3 && f[2] == "1" && f[1].parse::<usize>().map(|n| n >= 2).unwrap_or(false)
        })
        .collect();
    assert!(
        !same_chrom_families.is_empty(),
        "the --cross-chrom catalog must contain a same-chromosome multi-copy family\n{}",
        families
    );

    // The cross-chromosome family (c1:A + c2:X) must also still be present.
    let has_cross_chrom = families
        .lines()
        .skip(1)
        .any(|l| l.split('\t').nth(2) == Some("2"));
    assert!(has_cross_chrom, "the cross-chromosome family must remain in the catalog\n{}", families);

    assert!(copies.lines().count() > 1, "expected at least one copy row\n{}", copies);

    // The same-chrom supplement feature (and its sidecar) was removed: no sidecar must be written.
    assert!(
        std::fs::metadata(format!("{}.same_chrom_supplement.tsv", out)).is_err(),
        "the .same_chrom_supplement.tsv sidecar must no longer be created"
    );
}
