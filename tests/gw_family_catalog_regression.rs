use std::process::Command;

/// Byte-for-byte golden test: the default `--cross-chrom` catalog must reproduce the committed baseline in
/// `tests/fixtures/same_chrom_supplement/out_default.*`.
///
/// The run writes into a scratch dir under `target/` (`CARGO_TARGET_TMPDIR`, gitignored). It must NEVER
/// write to the golden path. The previous version of this test did exactly that -- it pointed `--out` at the
/// tracked golden files, so every run OVERWROTE the baseline it was meant to compare against and then
/// asserted only `lines().count() > 1` on what it had just written. It could not fail by construction, and
/// it silently masked a real regression: the same-chrom family `GWFAM1` (c1:250-460 + c1:380-550) stopped
/// being emitted once `--refine` became the default in 121b7ea, because `distinct_locus_reps` collapsed the
/// two overlapping copies -- `distinguishing_uniq` was never populated on the `--cross-chrom` path, so the
/// chi(H) read-evidence guard was inert. The committed golden still contained GWFAM1 and would have caught
/// it on the first real comparison.
#[test]
fn default_cross_chrom_output_is_unchanged() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let golden = "tests/fixtures/same_chrom_supplement/out_default";
    let scratch = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("gwcat_default");
    let out = scratch.to_str().expect("scratch path is valid UTF-8").to_string();

    let status = Command::new(bin)
        .arg("--bam")
        .arg("tests/fixtures/same_chrom_supplement/reads.bam")
        .arg("--fasta")
        .arg("tests/fixtures/same_chrom_supplement/genome.fa")
        .arg("--out")
        .arg(&out)
        .arg("--cross-chrom")
        .status()
        .expect("gw_family_catalog failed to run");
    assert!(status.success());

    // No sidecar without the (removed) supplement flag.
    assert!(
        std::fs::metadata(format!("{}.same_chrom_supplement.tsv", out)).is_err(),
        "sidecar must not be created when flag is absent"
    );

    for ext in ["families.tsv", "copies.tsv", "copies.fa"] {
        let want = std::fs::read_to_string(format!("{}.{}", golden, ext))
            .unwrap_or_else(|e| panic!("read committed golden {golden}.{ext}: {e}"));
        let got = std::fs::read_to_string(format!("{}.{}", out, ext))
            .unwrap_or_else(|e| panic!("read produced {out}.{ext}: {e}"));
        assert_eq!(
            got, want,
            "{ext} differs from the committed golden.\n--- committed ---\n{want}\n--- produced ---\n{got}"
        );
    }

    // The golden must contain BOTH families -- the cross-chrom one and the same-chrom one. Asserted
    // explicitly so that a future regeneration of the baseline cannot quietly bless a catalog that lost
    // the same-chrom family (which is precisely how the regression above went unnoticed).
    let families = std::fs::read_to_string(format!("{}.families.tsv", golden)).unwrap();
    let same_chrom = families
        .lines()
        .skip(1)
        .filter(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f.len() >= 3 && f[2] == "1" && f[1].parse::<usize>().map(|n| n >= 2).unwrap_or(false)
        })
        .count();
    assert_eq!(same_chrom, 1, "golden must retain the same-chrom family\n{families}");
}
