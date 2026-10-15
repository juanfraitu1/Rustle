//! `tools/rustle_pipeline.sh`: which stages a call runs, and in which order.
//!
//! The test runs the driver's own final stage dispatch (the last `case "$STAGE" in` .. `esac` at column 0, the block
//! figures/samples.py's driver_stage_code also treats as the dispatch) with every `stage_*` function replaced by a stub
//! that appends its name to a log. What is checked is the decision the dispatch makes, not what a stage does, so nothing
//! here needs a genome, a BAM or a heavy tool.
//!
//! 2026-10-06: the merge f8ae8a63 left two sets of `assemble) .. flag)` and `all)` arms in the one `case`. Bash takes the
//! first arm that matches, so `all` always ran the legacy catalog stage and never the candidates stage, whatever
//! `--candidates` and `--legacy-catalog` said.

use std::path::PathBuf;
use std::process::Command;

const DRIVER: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tools/rustle_pipeline.sh");
const STAGES: [&str; 7] = ["assemble", "families", "candidates", "catalog", "assign", "flag", "merged"];

/// The driver's final dispatch: from the last line that starts `case "$STAGE" in` to the next line that is `esac`.
fn dispatch_block() -> String {
    let text = std::fs::read_to_string(DRIVER).expect("read the driver");
    let lines: Vec<&str> = text.lines().collect();
    let start = lines
        .iter()
        .rposition(|l| l.starts_with("case \"$STAGE\" in"))
        .expect("the driver has no stage dispatch");
    let len = lines[start..]
        .iter()
        .position(|l| *l == "esac")
        .expect("the stage dispatch has no esac at column 0");
    lines[start..=start + len].join("\n")
}

/// A top-level function of the driver: from the line `NAME() {` to the next line that is `}`.
fn function_block(name: &str) -> String {
    let text = std::fs::read_to_string(DRIVER).expect("read the driver");
    let lines: Vec<&str> = text.lines().collect();
    let head = format!("{name}() {{");
    let start = lines.iter().position(|l| *l == head).unwrap_or_else(|| panic!("the driver has no {name}()"));
    let len = lines[start..].iter().position(|l| *l == "}").expect("function has no closing brace at column 0");
    lines[start..=start + len].join("\n")
}

/// The stages the dispatch ran for `stage` (in order) and its exit status. `tag` names the log, so tests do not share one.
/// `all` runs the driver's real `run_all`, which re-invokes the driver once per stage (`bash "$0" STAGE ..`): `$0` is a stub
/// that logs the stage, so the plan and its order are what is checked. An `--index` is given, so no index stage is planned.
fn dispatch(tag: &str, stage: &str, legacy_catalog: bool, candidates: bool) -> (Vec<String>, i32) {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch");
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let log = dir.join(format!("{tag}.log"));
    let _ = std::fs::remove_file(&log);
    let stub_driver = dir.join(format!("{tag}.driver.sh"));
    std::fs::write(&stub_driver, "echo \"$1\" >> \"$LOG\"\n").expect("write the stub driver");
    let stubs: Vec<String> = STAGES
        .iter()
        .map(|s| format!("stage_{s}() {{ echo {s} >> \"$LOG\"; }}"))
        .collect();
    let script = format!(
        "set -euo pipefail\n\
         export LOG=$1; STAGE=$2; LEGACY_CATALOG=$3; CANDIDATES=$4\n\
         OUT={out}; RESUME=0; NO_FLAG=0; INDEX=given.mmi; LOW_MEM=0; FASTA=genome.fa; ORIG_ARGS=()\n\
         say() {{ :; }}\n\
         all_env() {{ echo env; }}\n\
         {}\n\
         {}\n\
         {}\n",
        stubs.join("\n"),
        function_block("run_all"),
        dispatch_block(),
        out = dir.join(tag).display()
    );
    let out = Command::new("bash")
        .arg("-c")
        .arg(&script)
        .arg(&stub_driver)
        .arg(&log)
        .arg(stage)
        .arg(if legacy_catalog { "1" } else { "0" })
        .arg(if candidates { "1" } else { "0" })
        .output()
        .expect("bash failed to spawn");
    let ran = std::fs::read_to_string(&log)
        .unwrap_or_default()
        .lines()
        .map(String::from)
        .collect();
    (ran, out.status.code().unwrap_or(-1))
}

#[test]
fn all_runs_assemble_families_assign_flag_and_neither_catalog_nor_candidates() {
    let (ran, code) = dispatch("all_default", "all", false, false);
    assert_eq!(code, 0);
    assert_eq!(ran, ["assemble", "families", "assign", "flag"]);
}

#[test]
fn all_with_candidates_runs_the_candidates_stage_between_families_and_assign() {
    let (ran, code) = dispatch("all_candidates", "all", false, true);
    assert_eq!(code, 0);
    assert_eq!(ran, ["assemble", "families", "candidates", "assign", "flag"]);
}

#[test]
fn all_with_legacy_catalog_runs_the_catalog_stage_between_families_and_assign() {
    let (ran, code) = dispatch("all_legacy", "all", true, false);
    assert_eq!(code, 0);
    assert_eq!(ran, ["assemble", "families", "catalog", "assign", "flag"]);
}

#[test]
fn a_named_stage_runs_only_itself() {
    for s in STAGES {
        let (ran, code) = dispatch(&format!("only_{s}"), s, false, false);
        assert_eq!(code, 0, "stage {s}");
        assert_eq!(ran, [s], "stage {s}");
    }
}

#[test]
fn an_unknown_stage_is_refused_with_status_2_and_runs_nothing() {
    let (ran, code) = dispatch("unknown", "no-such-stage", false, false);
    assert_eq!(code, 2);
    assert!(ran.is_empty(), "ran {ran:?}");
}
