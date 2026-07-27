// Integration test for `gw_family_catalog --from-genome`: runs the binary on the real subset fixture
// (3 near-identical NCF1 copies + 2 unrelated decoys) and asserts the family is recovered from the
// genome alone (>=2-copy family emitted, all 3 NCF1 copies present) and the decoys are NOT grouped
// (precision control). Fast — the fixture is ~68 kb, so one minimap2 self-align in seconds.
use std::process::Command;

#[test]
fn from_genome_recovers_family_and_excludes_decoys() {
    if Command::new("minimap2").arg("--version").output().is_err() {
        eprintln!("minimap2 absent; skip");
        return;
    }
    let fa = "tests/fixtures/from_genome/subset.fa";
    let win = "tests/fixtures/from_genome/windows.bed";
    if std::fs::metadata(fa).is_err() || std::fs::metadata(win).is_err() {
        eprintln!("fixture absent; skip");
        return;
    }
    let dir = std::env::temp_dir().join("rustle_from_genome_it");
    let _ = std::fs::create_dir_all(&dir);
    let out = dir.join("dna");
    let out_s = out.to_str().unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_gw_family_catalog"))
        .args([
            "--from-genome", win,
            "--fasta", fa,
            "--min-identity", "0.90",
            "--out", out_s,
        ])
        .status()
        .expect("run gw_family_catalog --from-genome");
    assert!(status.success(), "binary exited non-zero");

    let fams = std::fs::read_to_string(format!("{out_s}.families.tsv")).expect("families.tsv");
    let copies = std::fs::read_to_string(format!("{out_s}.copies.tsv")).expect("copies.tsv");

    // >=1 family with n_copies >= 2 (column 2 = n_copies).
    let has_multicopy = fams
        .lines()
        .skip(1)
        .any(|l| l.split('\t').nth(1).and_then(|n| n.parse::<usize>().ok()).map_or(false, |n| n >= 2));
    assert!(has_multicopy, "expected a >=2-copy family from the genome; families.tsv:\n{fams}");

    // the three NCF1 copies must all appear as copies (column 4 = chrom/contig).
    let copy_contigs: Vec<&str> = copies
        .lines()
        .skip(1)
        .filter_map(|l| l.split('\t').nth(3))
        .collect();
    for want in ["NCF1", "NCF1B", "NCF1C"] {
        assert!(copy_contigs.contains(&want), "family copy {want} missing; copies.tsv:\n{copies}");
    }
    // precision: the unrelated decoys must NOT be grouped into any family.
    for decoy in ["DECOY1", "DECOY2"] {
        assert!(!copy_contigs.contains(&decoy), "decoy {decoy} wrongly grouped; copies.tsv:\n{copies}");
    }
}
