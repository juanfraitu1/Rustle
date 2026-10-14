//! `tools/rustle_pipeline.sh sd` (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 5) against a
//! STUB `utilities` that records each sd-scan call and writes the two tables the stage summarises: what is tested is
//! what the driver decides (the cache key, reuse, the prefix file, the arguments the stage needs).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const DRIVER: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tools/rustle_pipeline.sh");

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_sd").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("bin")).unwrap();
    let stub = d.join("bin/utilities");
    std::fs::write(
        &stub,
        "#!/bin/bash\n\
         if [ \"$1 $2\" = 'sd-scan --version-key' ]; then echo \"${STUB_KEY:-stub key 1}\"; exit 0; fi\n\
         echo \"$*\" >> \"$(dirname \"$0\")/calls.txt\"\n\
         [ -n \"${STUB_FAIL:-}\" ] && exit 3\n\
         while [ $# -gt 0 ]; do case \"$1\" in --out) OUT=$2; shift 2;; *) shift;; esac; done\n\
         printf '#h\\nchr1\\t0\\t2000\\tchr2\\t0\\t2000\\t+\\t9\\t99.0\\t99.0\\tok\\n' > \"$OUT.sd.pairs.tsv\"\n\
         printf '#h\\nF0\\t2\\t4000\\tchr1;chr2\\n' > \"$OUT.sd.families.tsv\"\n",
    )
    .unwrap();
    Command::new("chmod").args(["+x", stub.to_str().unwrap()]).status().unwrap();
    std::fs::write(d.join("g.fa"), ">chr1\nACGT\n").unwrap();
    d
}

fn sd(dir: &Path, env: &[(&str, &str)]) -> Output {
    let mut c = Command::new("bash");
    c.arg(DRIVER)
        .args(["sd", "--fasta"])
        .arg(dir.join("g.fa"))
        .arg("--out")
        .arg(dir.join("run"))
        .arg("--bin")
        .arg(dir.join("bin"))
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
