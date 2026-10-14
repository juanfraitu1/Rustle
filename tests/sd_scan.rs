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

/// The toy genome's tables through the hierarchy's readers (final-review fix m2): `parse_sd` reads sd-scan's pairs.tsv
/// by its shape (11 columns, last `ok` / `self` / `weak`) and `DupIndex::from_text` reads duplicons.tsv, so a change of
/// either file's format on the writer's side shows up here. The expected numbers are the toy genome's golden output
/// (src/sedef/scan.rs, PLANTED_*), which only a SCAN_VERSION bump may change.
#[test]
fn the_toy_genomes_tables_are_read_back_by_the_hierarchy_readers() {
    let d = scratch("cross_format");
    let fa = d.join("toy.fa");
    std::fs::write(&fa, rustle::sedef::scan::toy_genome_fasta(true)).unwrap();
    let o = Command::new(BIN).args(["sd-scan", "--genome"]).arg(&fa).arg("--out").arg(d.join("toy")).args(["--threads", "2"]).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    let pairs_text = std::fs::read_to_string(d.join("toy.sd.pairs.tsv")).unwrap();
    // columns by their header names, so the check does not repeat parse_sd's positions
    let header: Vec<&str> = pairs_text.lines().next().unwrap().trim_start_matches('#').split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap_or_else(|| panic!("no column {name}"));
    let (est, status) = (col("est_id"), col("status"));
    let ok_rows: Vec<Vec<&str>> = pairs_text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|f| f[status] == "ok")
        .collect();
    assert_eq!(ok_rows.len(), 5, "the toy genome has 5 ok pairs");

    let parsed = rustle::family_hierarchy::parse_sd(&pairs_text);
    assert_eq!(parsed.len(), 5, "parse_sd keeps the 5 ok pairs");
    for (p, f) in parsed.iter().zip(&ok_rows) {
        let n = |i: usize| f[i].parse::<u64>().unwrap();
        assert_eq!(p.a, (f[col("chr_a")].to_string(), n(col("start_a")), n(col("end_a"))));
        assert_eq!(p.b, (f[col("chr_b")].to_string(), n(col("start_b")), n(col("end_b"))));
        let want = f[est].parse::<f64>().unwrap() / 100.0;
        assert!((p.identity.unwrap() - want).abs() < 1e-12, "identity {:?} is est_id / 100 = {want}", p.identity);
    }

    // chrA:50000-51999 is the first elementary segment of the inverted copy, family F1 (src/sedef/scan.rs PLANTED_DUPLICONS):
    // [50000, 51999) is 1,999 bases, and F2 starts at 52000, outside the query [50000, 52000)
    let dup = rustle::family_hierarchy::DupIndex::from_text(&std::fs::read_to_string(d.join("toy.sd.duplicons.tsv")).unwrap());
    let got = dup.id_bases("chrA", 50_000, 52_000);
    assert_eq!(got.into_iter().collect::<Vec<_>>(), vec![("F1".to_string(), 1999u64)]);
}
