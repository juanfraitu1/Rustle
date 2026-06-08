//! End-to-end smoke on the chr19 BAM. Skipped when BAM/genome missing.
//!
//! The chr19 BAM has its unmapped reads stripped, so HMM rescue produces 0
//! candidates on it. The smoke goal is "no crash + pipeline runs to completion
//! without error". Real rescue testing needs full GGO.bam (Task 7.4).
//!
//! Note: env!("CARGO_BIN_EXE_rustle") requires that Cargo pre-build the rustle
//! binary before running tests. With `cargo test --profile quick` that happens
//! automatically because the [[test]] depends on [[bin]]. If running the test
//! binary directly, build with `cargo build --profile dev-opt` first.

use std::process::Command;

const BAM: &str = "/mnt/c/Users/jfris/Desktop/GGO_19.bam";
const GENOME: &str = "/mnt/c/Users/jfris/Desktop/GGO.fasta";

#[test]
fn hmm_mode_smoke_on_chr19() {
    if !std::path::Path::new(BAM).exists() {
        eprintln!("skipping: {} missing", BAM);
        return;
    }
    if !std::path::Path::new(GENOME).exists() {
        eprintln!("skipping: {} missing", GENOME);
        return;
    }

    let out_gtf = std::env::temp_dir().join("rustle_hmm_smoke.gtf");
    let status = Command::new(env!("CARGO_BIN_EXE_rustle"))
        .args([
            "-L",
            "--vg",
            "--vg-discover-novel",
            "--vg-discover-novel-mode",
            "hmm",
            "--genome-fasta",
            GENOME,
            "-o",
            out_gtf.to_str().unwrap(),
            BAM,
        ])
        .status()
        .expect("rustle failed to spawn");
    assert!(status.success(), "rustle exited with {:?}", status);

    let gtf = std::fs::read_to_string(&out_gtf).expect("read gtf");
    assert!(gtf.lines().count() > 0, "empty gtf");
}
