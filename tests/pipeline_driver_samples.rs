//! Multi-sample runs (docs/superpowers/specs/2026-11-02-multi-sample-design.md): `tools/rustle_pipeline.sh` hands a comma
//! `--bam` to `tools/rustle_samples.sh`. Tested here: the input checks, the pool and its reuse, and that one BAM never enters
//! the multi-sample path. Needs samtools on PATH (the pixi rustle env).
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const DRIVER: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tools/rustle_pipeline.sh");

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_samples").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A sorted, indexed BAM at `path` with reads `names` on contig c1 (or `sq` as the @SQ line).
fn bam(path: &Path, names: &[&str], sq: &str) {
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut s = format!("@HD\tVN:1.6\tSO:coordinate\n{sq}\n");
    for (i, n) in names.iter().enumerate() {
        s.push_str(&format!("{n}\t0\tc1\t{}\t60\t50M\t*\t0\t0\t*\t*\n", 100 + i * 10));
    }
    let sam = path.with_extension("sam");
    std::fs::write(&sam, s).unwrap();
    assert!(Command::new("samtools").args(["view", "-b", "-o"]).arg(path).arg(&sam).status().unwrap().success());
    assert!(Command::new("samtools").arg("index").arg(path).status().unwrap().success());
}

const SQ: &str = "@SQ\tSN:c1\tLN:10000";

fn run(stage: &str, bams: &str, extra: &[&str], out: &Path) -> Output {
    Command::new("bash")
        .arg(DRIVER)
        .arg(stage)
        .args(["--bam", bams, "--fasta", "g.fa", "--out"])
        .arg(out)
        .args(extra)
        .env("TMPDIR", out.parent().unwrap())
        .output()
        .unwrap()
}

fn err(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

fn two(d: &Path) -> String {
    bam(&d.join("a.bam"), &["r1", "r2"], SQ);
    bam(&d.join("b.bam"), &["s1"], SQ);
    format!("{},{}", d.join("a.bam").display(), d.join("b.bam").display())
}

#[test]
fn sample_names_must_match_the_bams_and_be_safe() {
    let d = scratch("names");
    let b = two(&d);
    let o = run("pool", &b, &["--samples", "A"], &d.join("run"));
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("--samples names 1 samples for 2 BAMs"), "{}", err(&o));
    let o = run("pool", &b, &["--samples", "A,b c"], &d.join("run"));
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("must match"), "{}", err(&o));
}

#[test]
fn same_basenames_ask_for_samples() {
    let d = scratch("basenames");
    bam(&d.join("x/s.bam"), &["r1"], SQ);
    bam(&d.join("y/s.bam"), &["s1"], SQ);
    let b = format!("{},{}", d.join("x/s.bam").display(), d.join("y/s.bam").display());
    let o = run("pool", &b, &[], &d.join("run"));
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("pass --samples"), "{}", err(&o));
}

#[test]
fn bams_must_share_the_genome() {
    let d = scratch("genome");
    bam(&d.join("a.bam"), &["r1"], SQ);
    bam(&d.join("b.bam"), &["s1"], "@SQ\tSN:c1\tLN:9999");
    let b = format!("{},{}", d.join("a.bam").display(), d.join("b.bam").display());
    let o = run("pool", &b, &[], &d.join("run"));
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("b.bam' (b) has other @SQ lines"), "{}", err(&o));
}

#[test]
fn a_shared_read_name_names_both_samples() {
    let d = scratch("dup");
    bam(&d.join("a.bam"), &["r1", "r2"], SQ);
    bam(&d.join("b.bam"), &["r2"], SQ);
    let b = format!("{},{}", d.join("a.bam").display(), d.join("b.bam").display());
    let o = run("pool", &b, &[], &d.join("run"));
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("read name r2 in a and b"), "{}", err(&o));
}

#[test]
fn v1_refuses_candidates_and_as_table() {
    let d = scratch("v1");
    let b = two(&d);
    for f in [&["--candidates"][..], &["--as-table", "t.tsv"][..]] {
        let o = run("all", &b, f, &d.join("run"));
        assert_eq!(o.status.code(), Some(2));
        assert!(err(&o).contains("not supported in multi-sample runs"), "{}", err(&o));
    }
}

#[test]
fn pool_merges_once_and_is_reused() {
    let d = scratch("pool with space");
    let b = two(&d);
    let out = d.join("run");
    let o = run("pool", &b, &["--samples", "A,B"], &out);
    assert!(o.status.success(), "{}", err(&o));
    let pooled = d.join("run.pooled.bam");
    let n = Command::new("samtools").args(["view", "-c"]).arg(&pooled).output().unwrap();
    assert_eq!(String::from_utf8_lossy(&n.stdout).trim(), "3");
    let t0 = std::fs::metadata(&pooled).unwrap().modified().unwrap();
    std::thread::sleep(std::time::Duration::from_millis(1100));
    let o = run("pool", &b, &["--samples", "A,B"], &out);
    assert!(o.status.success() && err(&o).contains("reused"), "{}", err(&o));
    assert_eq!(std::fs::metadata(&pooled).unwrap().modified().unwrap(), t0);
}

#[test]
fn one_bam_never_enters_the_multi_sample_path() {
    let d = scratch("single");
    bam(&d.join("a.bam"), &["r1"], SQ);
    // `pool` is not a stage of the single-sample driver
    let o = run("pool", d.join("a.bam").to_str().unwrap(), &[], &d.join("run"));
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("unknown stage pool"), "{}", err(&o));
    assert!(!d.join("run.samples.log").exists());
}

#[test]
fn a_shared_read_name_in_big_bams_still_gets_the_message() {
    // > a pipe buffer of names: the first duplicate must not end the check by SIGPIPE (review finding 1)
    let d = scratch("dupbig");
    let a: Vec<String> = (0..20000).map(|i| format!("read_{i:06}_with_a_long_name_to_fill_the_pipe")).collect();
    let ar: Vec<&str> = a.iter().map(|s| s.as_str()).collect();
    bam(&d.join("a.bam"), &ar, SQ);
    bam(&d.join("b.bam"), &[ar[0]], SQ);
    let b = format!("{},{}", d.join("a.bam").display(), d.join("b.bam").display());
    let o = run("pool", &b, &[], &d.join("run"));
    assert_eq!(o.status.code(), Some(2), "{}", err(&o));
    assert!(err(&o).contains("in a and b"), "{}", err(&o));
}

#[test]
fn the_copy_table_keeps_columns_when_a_sample_has_no_quant() {
    // review finding 2: a missing quant file must give that sample zeros, not shift the next sample's numbers
    let d = scratch("copytable");
    std::fs::create_dir_all(d.join("run.samples")).unwrap();
    let h = "family_id\tcopy_index\tcopy_tid\tcopy_chrom\tabundance\tn_reads_hard\n";
    std::fs::write(d.join("run.samples/C.assign.quant.tsv"), format!("{h}F1\t0\tt0\tc1\t0.75\t9\n")).unwrap();
    let o = Command::new("bash")
        .arg(concat!(env!("CARGO_MANIFEST_DIR"), "/tools/rustle_samples.sh"))
        .args(["copy-table", "--samples", "A,B,C", "--out"])
        .arg(d.join("run"))
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", err(&o));
    let t = std::fs::read_to_string(d.join("run.copy_counts.tsv")).unwrap();
    assert_eq!(
        t,
        "family_id\tcopy_index\tcopy_tid\tA_reads\tA_abundance\tB_reads\tB_abundance\tC_reads\tC_abundance\n\
         F1\t0\tt0\t0\t0\t0\t0\t9\t0.75\n"
    );
}
