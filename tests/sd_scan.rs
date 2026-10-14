//! `utilities sd-scan` (src/sedef/scan.rs) from the command line: the toy genome's tables, `--version-key`, and the
//! arguments a scan needs (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 4).

use std::path::PathBuf;
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_utilities");

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("sd_scan").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn sd_scan_writes_every_table_of_the_toy_genome() {
    let d = scratch("toy");
    let fa = d.join("toy.fa");
    std::fs::write(&fa, rustle::sedef::scan::toy_genome_fasta(true)).unwrap();
    let o = Command::new(BIN)
        .args(["sd-scan", "--genome"])
        .arg(&fa)
        .arg("--out")
        .arg(d.join("toy"))
        .args(["--threads", "2"])
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for ext in ["candidates.tsv", "blocks.tsv", "pairs.tsv", "anchors.bin", "segments.tsv", "families.tsv", "duplicons.tsv", "params.tsv"] {
        assert!(d.join(format!("toy.sd.{ext}")).exists(), "{ext}");
    }
    let params = std::fs::read_to_string(d.join("toy.sd.params.tsv")).unwrap();
    assert!(params.contains("\nfamilies\t5\n") && params.contains("\nduplicon_rows\t10\n"), "{params}");
}

#[test]
fn version_key_needs_no_genome_and_a_scan_needs_the_genome() {
    let o = Command::new(BIN).args(["sd-scan", "--version-key"]).output().unwrap();
    assert!(o.status.success());
    assert_eq!(
        String::from_utf8_lossy(&o.stdout).trim(),
        rustle::sedef::scan::ScanParams::default().version_key()
    );
    let o = Command::new(BIN).args(["sd-scan", "--out", "x"]).output().unwrap();
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("--genome"));
}
