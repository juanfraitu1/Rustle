//! `tools/rustle_pipeline.sh sd` (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 5; the cache
//! marker is the final-review fix I2) against a STUB `utilities` that records each sd-scan call and writes the eight
//! files the real `utilities sd-scan --out P` writes (`P.sd.*`): what is tested is what the driver decides (the cache
//! key, reuse, the DONE marker, the prefix file, the arguments the stage needs).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const DRIVER: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tools/rustle_pipeline.sh");

/// The eight files `utilities sd-scan --out P` writes as `P.sd.<name>` (src/sedef/scan.rs), plus the driver's key.tsv:
/// the nine files of a cache entry, as `genome.sd.<name>` and `key.tsv`.
const ENTRY_FILES: [&str; 9] = [
    "genome.sd.anchors.bin",
    "genome.sd.blocks.tsv",
    "genome.sd.candidates.tsv",
    "genome.sd.duplicons.tsv",
    "genome.sd.families.tsv",
    "genome.sd.pairs.tsv",
    "genome.sd.params.tsv",
    "genome.sd.segments.tsv",
    "key.tsv",
];

/// The stub: `--version-key` prints STUB_KEY (default `stub key 1`); a scan logs its arguments to calls.txt, exits 3
/// if STUB_FAIL is set (after writing nothing), else writes the eight files with fixed contents.
const STUB: &str = r##"#!/bin/bash
if [ "$1 $2" = 'sd-scan --version-key' ]; then echo "${STUB_KEY:-stub key 1}"; exit 0; fi
echo "$*" >> "$(dirname "$0")/calls.txt"
[ -n "${STUB_FAIL:-}" ] && exit 3
while [ $# -gt 0 ]; do case "$1" in --out) OUT=$2; shift 2;; *) shift;; esac; done
for t in candidates blocks segments duplicons params; do printf '#%s\n' "$t" > "$OUT.sd.$t.tsv"; done
printf 'ANCHORS\000\001\002\n' > "$OUT.sd.anchors.bin"
printf '#h\nchr1\t0\t2000\tchr2\t0\t2000\t+\t9\t99.0\t99.0\tok\n' > "$OUT.sd.pairs.tsv"
printf '#h\nF0\t2\t4000\tchr1;chr2\n' > "$OUT.sd.families.tsv"
"##;

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_sd").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("bin")).unwrap();
    let stub = d.join("bin/utilities");
    std::fs::write(&stub, STUB).unwrap();
    Command::new("chmod").args(["+x", stub.to_str().unwrap()]).status().unwrap();
    std::fs::write(d.join("g.fa"), ">chr1\nACGT\n").unwrap();
    d
}

fn sd(dir: &Path, env: &[(&str, &str)]) -> Output {
    sd_with(dir, &[], env)
}

/// The `sd` stage with extra driver arguments (e.g. `--no-cache`).
fn sd_with(dir: &Path, extra: &[&str], env: &[(&str, &str)]) -> Output {
    let mut c = Command::new("bash");
    c.arg(DRIVER)
        .args(["sd", "--fasta"])
        .arg(dir.join("g.fa"))
        .arg("--out")
        .arg(dir.join("run"))
        .arg("--bin")
        .arg(dir.join("bin"))
        .args(extra)
        .env_remove("RUSTLE_SD_CACHE_DIR")
        .env_remove("STUB_KEY")
        .env_remove("STUB_FAIL")
        .env("TMPDIR", dir);
    for (k, v) in env {
        c.env(k, v);
    }
    c.output().expect("bash failed to spawn")
}

fn calls(dir: &Path) -> Vec<String> {
    std::fs::read_to_string(dir.join("bin/calls.txt")).unwrap_or_default().lines().map(String::from).collect()
}

fn prefix(dir: &Path) -> String {
    std::fs::read_to_string(dir.join("run.sd.prefix")).unwrap().trim().to_string()
}

fn err(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

/// The cache entry directory a run's `run.sd.prefix` names (`<cache>/sd/<key>`).
fn entry(dir: &Path) -> PathBuf {
    Path::new(&prefix(dir)).parent().unwrap().to_path_buf()
}

/// `DONE` as `(name, bytes)` lines; panics on a line that is not `name<TAB>bytes`.
fn done_lines(entry: &Path) -> Vec<(String, u64)> {
    std::fs::read_to_string(entry.join("DONE"))
        .unwrap()
        .lines()
        .map(|l| {
            let (n, b) = l.split_once('\t').unwrap_or_else(|| panic!("DONE line {l:?} is not name<TAB>bytes"));
            (n.to_string(), b.parse().unwrap_or_else(|_| panic!("DONE line {l:?}: bytes is not a number")))
        })
        .collect()
}

/// Cut one byte off `f`; returns its former length.
fn cut_one_byte(f: &Path) -> u64 {
    let len = std::fs::metadata(f).unwrap().len();
    std::fs::OpenOptions::new().write(true).open(f).unwrap().set_len(len - 1).unwrap();
    len
}

/// `<cache>/sd/<32 hex>-<8 hex>/genome`
fn assert_key_dir(p: &str, cache: &Path) {
    let rest = p
        .strip_prefix(&format!("{}/sd/", cache.display()))
        .unwrap_or_else(|| panic!("{p} is not under {}/sd/", cache.display()));
    let (key, tail) = rest.split_once('/').unwrap();
    assert_eq!(tail, "genome");
    let (md5, ver) = key.split_once('-').unwrap();
    assert!(md5.len() == 32 && ver.len() == 8, "{key}");
    assert!(md5.chars().chain(ver.chars()).all(|c| c.is_ascii_hexdigit()), "{key}");
}

#[test]
fn sd_runs_sd_scan_once_per_genome_without_a_bam() {
    let d = scratch("once");
    let o = sd(&d, &[]);
    assert!(o.status.success(), "{}", err(&o));
    let c = calls(&d);
    assert_eq!(c.len(), 1);
    let p = prefix(&d);
    assert_key_dir(&p, &d.join("run.cache"));
    assert_eq!(c[0], format!("sd-scan --genome {} --out {p} --threads 4", d.join("g.fa").display()));
    let dir = Path::new(&p).parent().unwrap();
    assert!(dir.join("DONE").exists());
    assert!(std::fs::read_to_string(dir.join("key.tsv")).unwrap().starts_with("sd\t"));
    assert!(err(&o).contains("1 SD pairs (ok), 1 duplication families"), "{}", err(&o));
    let o = sd(&d, &[]);
    assert!(o.status.success(), "{}", err(&o));
    assert_eq!(calls(&d).len(), 1, "a cached genome is not scanned again");
    assert!(err(&o).contains("sd: cached"), "{}", err(&o));
    assert_eq!(prefix(&d), p);
}

#[test]
fn a_changed_genome_or_version_key_is_scanned_again() {
    let d = scratch("changed");
    assert!(sd(&d, &[]).status.success());
    let first = prefix(&d);
    std::fs::write(d.join("g.fa"), ">chr1\nACGTA\n").unwrap();
    assert!(sd(&d, &[]).status.success());
    let second = prefix(&d);
    assert_ne!(first, second);
    assert!(sd(&d, &[("STUB_KEY", "stub key 2")]).status.success());
    assert_ne!(prefix(&d), second);
    assert_eq!(calls(&d).len(), 3);
}

#[test]
fn a_failed_scan_is_not_cached() {
    let d = scratch("failed");
    let o = sd(&d, &[("STUB_FAIL", "1")]);
    assert!(!o.status.success());
    assert!(err(&o).contains("sd-scan failed"), "{}", err(&o));
    assert!(sd(&d, &[]).status.success());
    assert_eq!(calls(&d).len(), 2, "the failed scan left no DONE, so the next call scans");
}

#[test]
fn sd_needs_an_existing_fasta_and_takes_the_cache_dir_from_the_environment() {
    let d = scratch("args");
    let o = Command::new("bash").arg(DRIVER).args(["sd", "--out"]).arg(d.join("run")).arg("--bin").arg(d.join("bin")).output().unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("sd needs --fasta and --out"), "{}", err(&o));
    let o = Command::new("bash")
        .arg(DRIVER)
        .args(["sd", "--fasta"])
        .arg(d.join("missing.fa"))
        .arg("--out")
        .arg(d.join("run"))
        .arg("--bin")
        .arg(d.join("bin"))
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("is missing or empty"), "{}", err(&o));
    let cache = d.join("shared_cache");
    assert!(sd(&d, &[("RUSTLE_SD_CACHE_DIR", cache.to_str().unwrap())]).status.success());
    assert_key_dir(&prefix(&d), &cache);
}

#[test]
fn done_lists_every_file_of_the_entry_with_its_size() {
    let d = scratch("done_format");
    let o = sd(&d, &[]);
    assert!(o.status.success(), "{}", err(&o));
    let e = entry(&d);
    let mut lines = done_lines(&e);
    lines.sort();
    assert_eq!(lines.iter().map(|(n, _)| n.as_str()).collect::<Vec<_>>(), ENTRY_FILES, "one line per file of the entry");
    for (n, b) in &lines {
        assert_eq!(std::fs::metadata(e.join(n)).unwrap().len(), *b, "{n}");
        assert!(*b > 0, "{n} is empty");
    }
    assert!(!e.join("DONE.tmp").exists(), "DONE.tmp is renamed to DONE");
}

#[test]
fn a_truncated_cached_file_is_scanned_again() {
    let d = scratch("truncated");
    assert!(sd(&d, &[]).status.success());
    let (p, e) = (prefix(&d), entry(&d));
    let pairs = e.join("genome.sd.pairs.tsv");
    let len = cut_one_byte(&pairs);
    let o = sd(&d, &[]);
    assert!(o.status.success(), "{}", err(&o));
    assert_eq!(calls(&d).len(), 2, "a cached file shorter than its DONE size is a miss: the next call scans");
    assert!(!err(&o).contains("sd: cached"), "{}", err(&o));
    assert_eq!(std::fs::metadata(&pairs).unwrap().len(), len, "the rescan restored the file");
    assert_eq!(prefix(&d), p);
    assert!(sd(&d, &[]).status.success());
    assert_eq!(calls(&d).len(), 2, "the entry is marked done again, so the next call hits");
}

#[test]
fn an_empty_done_a_missing_file_or_a_longer_file_is_a_miss() {
    // an empty DONE is what the stage wrote before the marker listed the files
    let d = scratch("miss_empty");
    assert!(sd(&d, &[]).status.success());
    std::fs::write(entry(&d).join("DONE"), "").unwrap();
    assert!(sd(&d, &[]).status.success());
    assert_eq!(calls(&d).len(), 2, "an empty DONE is a miss");
    assert!(!done_lines(&entry(&d)).is_empty(), "and the rescan re-marks the entry");

    let d = scratch("miss_missing");
    assert!(sd(&d, &[]).status.success());
    std::fs::remove_file(entry(&d).join("genome.sd.segments.tsv")).unwrap();
    assert!(sd(&d, &[]).status.success());
    assert_eq!(calls(&d).len(), 2, "a file DONE lists but the entry lacks is a miss");
    assert!(entry(&d).join("genome.sd.segments.tsv").exists());

    let d = scratch("miss_longer");
    assert!(sd(&d, &[]).status.success());
    let anchors = entry(&d).join("genome.sd.anchors.bin");
    let mut bytes = std::fs::read(&anchors).unwrap();
    bytes.push(0);
    std::fs::write(&anchors, bytes).unwrap();
    assert!(sd(&d, &[]).status.success());
    assert_eq!(calls(&d).len(), 2, "a file longer than its DONE size is a miss");

    let d = scratch("miss_garbled");
    assert!(sd(&d, &[]).status.success());
    std::fs::write(entry(&d).join("DONE"), "genome.sd.pairs.tsv\n").unwrap();
    assert!(sd(&d, &[]).status.success());
    assert_eq!(calls(&d).len(), 2, "a DONE line without a size is a miss");
}

#[test]
fn a_failed_rescan_leaves_neither_prefix_nor_done() {
    let d = scratch("failed_rescan");
    assert!(sd(&d, &[]).status.success());
    let e = entry(&d);
    cut_one_byte(&e.join("genome.sd.pairs.tsv"));
    let o = sd(&d, &[("STUB_FAIL", "1")]);
    assert_eq!(o.status.code(), Some(1), "{}", err(&o));
    assert!(err(&o).contains("sd-scan failed"), "{}", err(&o));
    assert_eq!(calls(&d).len(), 2, "the damaged entry was scanned again");
    assert!(!d.join("run.sd.prefix").exists(), "a failed run leaves no prefix");
    assert!(!e.join("DONE").exists(), "the entry is not marked done");
    assert!(sd(&d, &[]).status.success());
    assert_eq!(calls(&d).len(), 3, "without DONE the next call scans");
    assert_eq!(done_lines(&e).len(), 9);
    assert_eq!(prefix(&d), e.join("genome").to_str().unwrap());
}

#[test]
fn a_failed_scan_of_another_genome_removes_the_previous_prefix() {
    let d = scratch("stale_prefix");
    assert!(sd(&d, &[]).status.success());
    assert!(d.join("run.sd.prefix").exists());
    std::fs::write(d.join("g.fa"), ">chr1\nACGTA\n").unwrap();
    let o = sd(&d, &[("STUB_FAIL", "1")]);
    assert_eq!(o.status.code(), Some(1), "{}", err(&o));
    assert!(!d.join("run.sd.prefix").exists(), "the first genome's prefix must not outlive a failed run");
}

#[test]
fn no_cache_scans_a_warm_entry_again_and_marks_it_again() {
    let d = scratch("no_cache");
    assert!(sd(&d, &[]).status.success());
    let (p, e) = (prefix(&d), entry(&d));
    let o = sd_with(&d, &["--no-cache"], &[]);
    assert!(o.status.success(), "{}", err(&o));
    assert_eq!(calls(&d).len(), 2, "--no-cache scans a warm entry");
    assert!(!err(&o).contains("sd: cached"), "{}", err(&o));
    assert_eq!(prefix(&d), p, "into the same entry");
    assert_eq!(done_lines(&e).len(), 9, "which is marked done again");
    assert!(sd(&d, &[]).status.success());
    assert_eq!(calls(&d).len(), 2, "and a plain call hits it");
    assert!(sd_with(&d, &["--cache"], &[]).status.success());
    assert_eq!(calls(&d).len(), 2, "--cache is the default");
}
