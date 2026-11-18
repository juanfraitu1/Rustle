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
/// `all` runs the driver's real `stage_all`, which re-invokes the driver once per stage (`bash "$0" STAGE ..`): `$0` is a stub
/// that logs the stage, so the plan and its order are what is checked. An `--index` is given, so no index stage is planned.
fn dispatch(tag: &str, stage: &str, legacy_catalog: bool, candidates: bool) -> (Vec<String>, i32) {
    let (ran, code, _) = dispatch_env(tag, stage, legacy_catalog, candidates, &[]);
    (ran, code)
}

/// `dispatch` with extra environment variables; also returns the driver's stderr.
fn dispatch_env(tag: &str, stage: &str, legacy_catalog: bool, candidates: bool, env: &[(&str, &str)]) -> (Vec<String>, i32, String) {
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
         export RUSTLE_MINIMAP2=true  # stage_all checks for minimap2 up front
         OUT={out}; RESUME=0; NO_FLAG=0; INDEX=given.mmi; LOW_MEM=0; FASTA=genome.fa; ORIG_ARGS=(); FROM_STAGE=; TO_STAGE=\n\
         say() {{ :; }}\n\
         stage_all_env() {{ echo env; }}\n\
         {}\n\
         {}\n\
         {}\n",
        stubs.join("\n"),
        function_block("stage_all"),
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
        .envs(env.iter().copied())
        .output()
        .expect("bash failed to spawn");
    // `flag_scan` runs in the background from the end of assemble (round-39 audit), so its line may land anywhere after
    // assemble's: the plan's order is checked without it, and its single run by `all_runs_..._flag_scan_once`
    let ran = std::fs::read_to_string(&log)
        .unwrap_or_default()
        .lines()
        .filter(|l| *l != "flag_scan" && *l != "flag_bg")
        .map(String::from)
        .collect();
    (ran, out.status.code().unwrap_or(-1), String::from_utf8_lossy(&out.stderr).into_owned())
}

/// Audit round 52 (2026-11-11): assign's skip lists (RUSTLE_EXCLUDE_BED, RUSTLE_ANNOTATED_ONLY) are built by a Python
/// helper. Without python3 `all` used to run assemble and families (13 min on A119b) and then stop in assign with a bare
/// exit 2. It now refuses at the start and names what is missing; without a skip list python3 is not needed at all.
#[test]
fn all_refuses_at_the_start_when_a_skip_list_needs_a_python_that_is_not_there() {
    let none = [("RUSTLE_PYTHON", "no-such-python-xyz")];
    let (ran, code, _) = dispatch_env("py_unneeded", "all", false, false, &[none[0], ("RUSTLE_FAMREPS", "0")]);
    assert_eq!((code, ran.first().map(String::as_str)), (0, Some("assemble")), "no skip list, no famreps: python is not needed");
    // row 1662: families' extra representatives (default) are built by Python helpers
    let (ran, code, err) = dispatch_env("py_famreps", "all", false, false, &none);
    assert!(code == 2 && ran.is_empty() && err.contains("RUSTLE_FAMREPS"), "{code} {ran:?} {err}");
    for var in ["RUSTLE_EXCLUDE_BED", "RUSTLE_EXCLUDE_SOFTMASK", "RUSTLE_EXCLUDE_FASTA", "RUSTLE_ANNOTATED_ONLY"] {
        let (ran, code, err) = dispatch_env(&format!("py_{var}"), "all", false, false, &[none[0], ("RUSTLE_FAMREPS", "0"), (var, "x.bed")]);
        assert_eq!(code, 2, "{var}: {err}");
        assert!(ran.is_empty(), "{var}: ran {ran:?}");
        assert!(err.contains("no-such-python-xyz") && err.contains(var), "{var}: {err}");
    }
}

#[test]
fn all_runs_assemble_families_assign_flag_and_neither_catalog_nor_candidates() {
    let (ran, code) = dispatch("all_default", "all", false, false);
    assert_eq!(code, 0);
    // round 55: flag (scan, align, verdict) runs in the background from the end of assemble, as `flag_bg`
    assert_eq!(ran, ["assemble", "families", "assign"]);
}

#[test]
fn all_runs_the_flag_scan_once_in_the_background_after_assemble() {
    let (_, code) = dispatch("all_scan", "all", false, false);
    assert_eq!(code, 0);
    let log = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("all_scan.log");
    let ran: Vec<String> = std::fs::read_to_string(log).unwrap().lines().map(String::from).collect();
    // round 55: the whole flag stage is that background job (`flag_bg` = scan, then align + verdict); no second scan or flag
    assert_eq!(ran.iter().filter(|l| *l == "flag_bg").count(), 1, "{ran:?}");
    assert!(!ran.iter().any(|l| l == "flag" || l == "flag_scan"), "{ran:?}");
    let pos = |s: &str| ran.iter().position(|l| l == s).unwrap();
    assert!(pos("assemble") < pos("flag_bg"), "{ran:?}");
    // with --candidates the verdict reads the candidates table: only the scan runs early, flag stays last
    let (ran, code) = dispatch("all_scan_cand", "all", false, true);
    assert_eq!((code, ran.last().map(String::as_str)), (0, Some("flag")), "{ran:?}");
    let log = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("all_scan_cand.log");
    let all = std::fs::read_to_string(log).unwrap();
    assert!(all.lines().any(|l| l == "flag_scan") && !all.lines().any(|l| l == "flag_bg"), "{all}");
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
    assert_eq!(ran, ["assemble", "families", "catalog", "assign"]);
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

/// Round-10 audit (2026-11-07): `all`'s resume. Each stage keeps its own settings record (PREFIX.<stage>.env), written only
/// when it succeeds. A failed stage is redone without redoing the ones before it (one shared record, removed before every
/// stage, used to send the next resume back to assemble). A product newer than its stage's time record (a direct stage call
/// rewrote it) is redone. A settings change redoes everything.
#[test]
fn all_resumes_at_the_failed_stage_and_redoes_a_product_rewritten_outside_it() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("resume");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let out = dir.join("run");
    let log = dir.join("stages.log");
    let stub_driver = dir.join("driver.sh");
    std::fs::write(
        &stub_driver,
        "echo \"$1\" >> \"$LOG\"\n\
         [ \"${FAIL_STAGE:-}\" = \"$1\" ] && exit 3\n\
         case \"$1\" in assemble) p=gtf;; earlyends) p=earlyends.gtf;; capstarts) p=capstarts.gtf;; families) p=fam.copies.tsv;; assign) p=assign.assignments.tsv;; flag|flag_bg) p=flag.missing_copy.tsv;; esac\n\
         echo x > \"$OUT.$p\"\n",
    )
    .unwrap();
    let script = format!(
        "set -euo pipefail\n\
         export LOG=$1 OUT={out} RUSTLE_MINIMAP2=true; STAGE=all; LEGACY_CATALOG=0; CANDIDATES=0\n\
         RESUME=1; NO_FLAG=0; INDEX=given.mmi; LOW_MEM=0; FASTA=genome.fa; ORIG_ARGS=(); FROM_STAGE=; TO_STAGE=\n\
         say() {{ :; }}\n\
         stage_all_env() {{ echo \"env $SETTING\"; }}\n\
         {}\n\
         {}\n",
        function_block("stage_all"),
        dispatch_block(),
        out = out.display()
    );
    let run = |setting: &str, fail: &str| -> (Vec<String>, i32) {
        let _ = std::fs::remove_file(&log);
        let o = Command::new("bash")
            .arg("-c")
            .arg(&script)
            .arg(&stub_driver)
            .arg(&log)
            .env("SETTING", setting)
            .env("FAIL_STAGE", fail)
            .output()
            .expect("bash failed to spawn");
        let ran = std::fs::read_to_string(&log).unwrap_or_default().lines().filter(|l| *l != "flag_scan" && *l != "flag_bg").map(String::from).collect();
        (ran, o.status.code().unwrap_or(-1))
    };
    assert_eq!(run("a", ""), (vec!["assemble".into(), "families".into(), "assign".into(), "earlyends".into(), "capstarts".into()], 0), "flag ran in the background; the early-end and cap-start assemblies (default side products, rows 1717-1722) are written last");
    assert_eq!(run("a", ""), (vec![], 0), "nothing changed: every stage is skipped");
    // a direct `families` call rewrote its product after the recorded run: its time record is older than the product
    let earlier = std::time::SystemTime::now() - std::time::Duration::from_secs(60);
    std::fs::File::options().write(true).open(format!("{}.families.time", out.display())).unwrap().set_modified(earlier).unwrap();
    assert_eq!(run("a", ""), (vec!["families".into(), "assign".into(), "flag".into()], 0));
    // new settings, assign fails: then the same settings again resume AT assign
    let (ran, code) = run("b", "assign");
    assert_eq!((ran, code), (vec!["assemble".into(), "families".into(), "assign".into()], 3));
    // the early-end assembly of the first run is older than the assembly just rewritten (set by hand: the runs here share a second)
    for side in ["earlyends", "capstarts"] {
        std::fs::File::options().write(true).open(format!("{}.{side}.gtf", out.display())).unwrap().set_modified(earlier).unwrap();
    }
    assert_eq!(run("b", ""), (vec!["assign".into(), "flag".into(), "earlyends".into(), "capstarts".into()], 0), "assemble and families are not redone; the side assemblies follow the new assembly");
}

/// Round-41 audit: when a stage fails, `all` stops its background flag scan (the scan writes PREFIX.flag_scan.* in place, and a
/// rerun's own scan would race it). The stub scan sleeps, then logs `scan_end`: after the failed run it must never get there.
#[test]
fn a_failed_stage_stops_the_background_flag_scan() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch_bg_kill");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let log = dir.join("stages.log");
    let stub_driver = dir.join("driver.sh");
    std::fs::write(
        &stub_driver,
        "echo \"$1\" >> \"$LOG\"\n\
         [ \"$1\" = flag_bg ] && { ( ( sleep 3; echo grand_end >> \"$LOG\" ) & wait ); sleep 3; echo scan_end >> \"$LOG\"; }\n\
         [ \"$1\" = assign ] && exit 3\n\
         exit 0\n",
    )
    .unwrap();
    let script = format!(
        "set -euo pipefail\n\
         export LOG=$1 OUT={out} RUSTLE_MINIMAP2=true; STAGE=all; LEGACY_CATALOG=0; CANDIDATES=0\n\
         RESUME=0; NO_FLAG=0; INDEX=given.mmi; LOW_MEM=0; FASTA=genome.fa; ORIG_ARGS=(); FROM_STAGE=; TO_STAGE=\n\
         say() {{ :; }}\n\
         stage_all_env() {{ echo env; }}\n\
         {}\n\
         {}\n",
        function_block("stage_all"),
        dispatch_block(),
        out = dir.join("run").display()
    );
    let o = Command::new("bash").arg("-c").arg(&script).arg(&stub_driver).arg(&log).output().expect("bash failed to spawn");
    assert_eq!(o.status.code(), Some(3), "assign's exit status");
    std::thread::sleep(std::time::Duration::from_secs(5));
    let ran = std::fs::read_to_string(&log).unwrap();
    assert!(ran.lines().any(|l| l == "flag_bg"), "the scan started: {ran}");
    assert!(!ran.lines().any(|l| l == "scan_end"), "the scan outlived the failed run: {ran}");
    // round 60: flag's align step has a minimap2 child of its own; the whole tree is stopped, not only the first level
    assert!(!ran.lines().any(|l| l == "grand_end"), "a grandchild of the background job outlived the failed run: {ran}");
}

/// Audit round 46: `all --guide` writes PREFIX.guided.gtf as a side product of assemble. A resumed run that skips assemble
/// must still write a missing one, and writing it must not make the later stages rerun.
#[test]
fn all_with_a_guide_writes_the_guided_assembly_even_when_assemble_is_skipped() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("guided");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, log, stub_driver, guide) = (dir.join("run"), dir.join("stages.log"), dir.join("driver.sh"), dir.join("guide.gtf"));
    std::fs::write(&guide, "x\n").unwrap();
    std::fs::write(
        &stub_driver,
        "echo \"$1\" >> \"$LOG\"\n\
         case \"$1\" in assemble) p=gtf;; earlyends) p=earlyends.gtf;; capstarts) p=capstarts.gtf;; guided) p=guided.gtf;; families) p=fam.copies.tsv;; assign) p=assign.assignments.tsv;; flag|flag_bg) p=flag.missing_copy.tsv;; esac\n\
         echo x > \"$OUT.$p\"\n",
    )
    .unwrap();
    let script = format!(
        "set -euo pipefail\n\
         export LOG=$1 OUT={out} RUSTLE_MINIMAP2=true; STAGE=all; LEGACY_CATALOG=0; CANDIDATES=0; GUIDE={guide}\n\
         RESUME=1; NO_FLAG=0; INDEX=given.mmi; LOW_MEM=0; FASTA=genome.fa; ORIG_ARGS=(); FROM_STAGE=; TO_STAGE=\n\
         say() {{ :; }}\n\
         stage_all_env() {{ echo env; }}\n\
         {}\n\
         {}\n",
        function_block("stage_all"),
        dispatch_block(),
        out = out.display(),
        guide = guide.display()
    );
    let run = || -> (Vec<String>, i32) {
        let _ = std::fs::remove_file(&log);
        let o = Command::new("bash").arg("-c").arg(&script).arg(&stub_driver).arg(&log).output().expect("bash failed to spawn");
        let ran = std::fs::read_to_string(&log).unwrap_or_default().lines().filter(|l| *l != "flag_scan" && *l != "flag_bg").map(String::from).collect();
        (ran, o.status.code().unwrap_or(-1))
    };
    let all: Vec<String> = ["assemble", "guided", "families", "assign", "earlyends", "capstarts"].iter().map(|s| s.to_string()).collect();
    assert_eq!(run(), (all, 0));
    assert_eq!(run(), (vec![], 0), "nothing changed: the guided assembly is not rewritten");
    std::fs::remove_file(format!("{}.guided.gtf", out.display())).unwrap();
    assert_eq!(run(), (vec!["guided".into()], 0), "only the missing guided assembly is written");
}

/// Audit round 54 (2026-11-11): assign reads 24 regions at once under a 16M-record budget when the node is large (>= 48
/// cores, >= 128 GB) AND the satellite exclusion is set, the configuration measured on A119b (286 s against 484 s at 8 / 4M,
/// tables identical, peak RSS 40 GB against 34 GB). Anywhere else the r1543 defaults stay: without the exclusion they
/// peaked at 199 GB on the satellite families.
#[test]
fn assign_reads_more_regions_at_once_only_on_the_measured_configuration() {
    let run = |args: &str| {
        let script = format!("set -euo pipefail\n{}\nassign_parallel_defaults {args}\n", function_block("assign_parallel_defaults"));
        let out = Command::new("bash").arg("-c").arg(&script).output().expect("bash failed to spawn");
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        String::from_utf8_lossy(&out.stdout).trim().to_string()
    };
    // cores, memory in kB, exclusion BED, low-mem
    assert_eq!(run("72 527502548 censat.bed 0"), "24 16000000");
    assert_eq!(run("48 134217728 censat.bed 0"), "24 16000000");
    assert_eq!(run("72 527502548 '' 0"), "8 4000000", "no exclusion: the satellite families set the memory");
    assert_eq!(run("32 527502548 censat.bed 0"), "8 4000000", "few cores");
    assert_eq!(run("72 65000000 censat.bed 0"), "8 4000000", "little memory");
    assert_eq!(run("72 527502548 censat.bed 1"), "8 4000000", "--low-mem sets its own");
}

/// Audit round 58: assemble sweeps 24 contigs at once under an 80M-record budget on a large node (A119b: 291 s at 8 / 20M,
/// 237 s at 16 / 40M, 212 s at 25 / 80M, the three GTFs identical, 9 / 11 / 13 GB per process); elsewhere 8 / 20M as before.
#[test]
fn assemble_sweeps_more_contigs_at_once_only_on_a_large_node() {
    let run = |args: &str| {
        let script = format!("set -euo pipefail\n{}\nassemble_parallel_defaults {args}\n", function_block("assemble_parallel_defaults"));
        let out = Command::new("bash").arg("-c").arg(&script).output().expect("bash failed to spawn");
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        String::from_utf8_lossy(&out.stdout).trim().to_string()
    };
    // cores, memory in kB, low-mem
    assert_eq!(run("72 527502548 0"), "24 80000000");
    assert_eq!(run("48 134217728 0"), "24 80000000");
    assert_eq!(run("32 527502548 0"), "8 20000000");
    assert_eq!(run("72 65000000 0"), "8 20000000");
    assert_eq!(run("72 527502548 1"), "8 20000000");
}

/// Audit round 66: families needs PREFIX.families.gtf, not the layered PREFIX.final.gtf, whose run ends about 25 s after the
/// base run on A119b. Under `all` (not `--low-mem`) assemble leaves the layered run to finish in the background
/// (`RP_LAYERS_BG=1`, pid in PREFIX.final.pid); `all` waits for it before it ends and stops it when a stage fails.
#[test]
fn all_lets_the_layered_assembly_finish_under_families() {
    let run = |tag: &str, low_mem: u8, fail: &str| -> (String, i32, PathBuf) {
        let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch_layers").join(tag);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create scratch dir");
        let (log, stub_driver, out) = (dir.join("stages.log"), dir.join("driver.sh"), dir.join("run"));
        std::fs::write(
            &stub_driver,
            "echo \"$1\" >> \"$LOG\"\n\
             if [ \"$1\" = assemble ]; then\n\
               if [ \"${RP_LAYERS_BG:-0}\" = 1 ]; then ( sleep 2; echo x > \"$OUT.final.gtf\" ) > /dev/null 2>&1 & echo $! > \"$OUT.final.pid\"\n\
               else echo layers_in_stage >> \"$LOG\"; fi\n\
             fi\n\
             if [ \"$1\" = families ] && [ -e \"$OUT.final.gtf\" ]; then echo families_saw_final >> \"$LOG\"; fi\n\
             [ \"$1\" = \"$FAIL\" ] && exit 3\n\
             exit 0\n",
        )
        .unwrap();
        let script = format!(
            "set -euo pipefail\n\
             export LOG=$1 OUT={out} RUSTLE_MINIMAP2=true FAIL={fail}; STAGE=all; LEGACY_CATALOG=0; CANDIDATES=0\n\
             RESUME=0; NO_FLAG=0; INDEX=given.mmi; LOW_MEM={low_mem}; FASTA=genome.fa; ORIG_ARGS=(); FROM_STAGE=; TO_STAGE=\n\
             say() {{ :; }}\n\
             layered_report() {{ :; }}\n\
             stage_all_env() {{ echo env; }}\n\
             {}\n\
             {}\n",
            function_block("stage_all"),
            dispatch_block(),
            out = out.display()
        );
        let o = Command::new("bash").arg("-c").arg(&script).arg(&stub_driver).arg(&log).output().expect("bash failed to spawn");
        let ran = std::fs::read_to_string(&log).unwrap_or_default();
        (format!("{ran}{}", String::from_utf8_lossy(&o.stderr)), o.status.code().unwrap_or(-1), out)
    };
    let fin = |out: &PathBuf| PathBuf::from(format!("{}.final.gtf", out.display())).exists();
    let (ran, code, out) = run("bg", 0, "none");
    assert_eq!(code, 0, "{ran}");
    assert!(!ran.contains("layers_in_stage") && !ran.contains("families_saw_final"), "families did not wait: {ran}");
    assert!(fin(&out), "`all` ended before the layered assembly was written");
    assert!(!PathBuf::from(format!("{}.final.pid", out.display())).exists());
    let (ran, code, _) = run("lowmem", 1, "none");
    assert!(code == 0 && ran.contains("layers_in_stage"), "--low-mem keeps the layered run inside assemble: {ran}");
    let (ran, code, out) = run("fail", 0, "families");
    assert_eq!(code, 3, "{ran}");
    std::thread::sleep(std::time::Duration::from_secs(3));
    assert!(!fin(&out), "the layered run outlived the failed run");
}

/// Register 1628 (opt-in, user decision 2026-11-11): `siblings` writes PREFIX.sibling.gtf from the assembly, the family copy
/// table and the best-AS table; it refuses before running anything when one of them is missing.
#[test]
fn siblings_needs_the_assembly_and_the_families_and_then_writes_the_sibling_gtf() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("siblings");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, log, py) = (dir.join("run"), dir.join("py.log"), dir.join("py.sh"));
    // the stub interpreter logs the script it was given; the `emit` call writes its last argument
    std::fs::write(&py, "#!/bin/bash\nbasename \"$1\" >> \"$LOG\"\nif [ \"$2\" = emit ]; then echo x > \"${@: -1}\"; fi\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport LOG={log} RUSTLE_PYTHON='bash {py}'\nOUT={out}; BAM=r.bam; FASTA=g.fa; THREADS=2; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\nstage_siblings\n",
        function_block("stage_siblings"),
        log = log.display(),
        py = py.display(),
        out = out.display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("siblings needs") && !log.exists());
    for p in ["gtf", "fam.copies.tsv", "molecules.tsv"] {
        std::fs::write(format!("{}.{p}", out.display()), "x\n").unwrap();
    }
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(&log).unwrap(), "project_models.py\nsibling_chains.py\n");
    assert!(std::path::Path::new(&format!("{}.sibling.gtf", out.display())).exists());
}

/// Register 1640 E5 (opt-in, user decision 2026-11-12): `flairlike` refuses without a BAM and genome, then runs the first pass
/// and the realignment and writes PREFIX.flairlike.gtf.
#[test]
fn flairlike_needs_the_bam_and_genome_and_then_writes_the_flairlike_gtf() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("flairlike");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, log, py) = (dir.join("run"), dir.join("py.log"), dir.join("py.sh"));
    std::fs::write(&py, "#!/bin/bash\necho \"$2\" >> \"$LOG\"\nif [ \"$2\" = gtf ]; then cat > /dev/null; echo x > \"${@: -1}\"; fi\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport LOG={log} RUSTLE_PYTHON='bash {py}'\nOUT={out}; BAM={bam}; FASTA={bam}; THREADS=2; RP_DRIVER_DIR=tools\nsay() {{ :; }}\nsamtools() {{ :; }}\nminimap2() {{ cat; }}\n{}\nstage_flairlike\n",
        function_block("stage_flairlike"),
        log = log.display(),
        py = py.display(),
        out = out.display(),
        bam = dir.join("in").display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("flairlike needs") && !log.exists());
    std::fs::write(dir.join("in"), "x\n").unwrap();
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(&log).unwrap(), "firstpass\ngtf\n");
    assert!(std::path::Path::new(&format!("{}.flairlike.gtf", out.display())).exists());
}

/// Opt-in `altends` (2026-11-12): refuses without the assembly, then writes PREFIX.altends.gtf and leaves PREFIX.gtf alone.
#[test]
fn altends_needs_the_assembly_and_then_writes_the_altends_gtf() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("altends");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py, inp) = (dir.join("run"), dir.join("py.sh"), dir.join("in"));
    std::fs::write(&py, "#!/bin/bash\necho annotated > \"$5\"\necho done\n").unwrap();
    std::fs::write(&inp, "x\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}'\nOUT={out}; BAM={inp}; FASTA={inp}; THREADS=2; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\nstage_altends\n",
        function_block("stage_altends"),
        py = py.display(),
        out = out.display(),
        inp = inp.display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("altends needs"));
    let gtf = format!("{}.gtf", out.display());
    std::fs::write(&gtf, "model\n").unwrap();
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(format!("{}.altends.gtf", out.display())).unwrap(), "annotated\n");
    assert_eq!(std::fs::read_to_string(&gtf).unwrap(), "model\n");
}

/// Opt-in `relatives` (register 1676 / 1737): refuses without the families output or without a SEDEF BED, then writes the table.
#[test]
fn relatives_needs_the_families_output_and_a_sedef_bed() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("relatives");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py, fa) = (dir.join("run"), dir.join("py.sh"), dir.join("g.fa"));
    // stand-in for relative_links.py: `--sedef BED PREFIX OUT`
    std::fs::write(&py, "#!/bin/bash\necho links > \"$5\"; echo done\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}'\nOUT={out}; FASTA={fa}; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\nstage_relatives\n",
        function_block("stage_relatives"),
        py = py.display(),
        out = out.display(),
        fa = fa.display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).env_remove("RUSTLE_SEDEF_BED").output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("relatives needs"));
    for f in ["copies.tsv", "loci.gff3", "loci.gated.paf"] {
        std::fs::write(format!("{}.fam.{f}", out.display()), "x\n").unwrap();
    }
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("SEDEF BED"));
    std::fs::write(format!("{}.sedef.bed", fa.display()), "x\n").unwrap();
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(format!("{}.fam.relatives.tsv", out.display())).unwrap(), "links\n");
}

/// Opt-in `sdcore` (register 1756): refuses without the copies table or without a SEDEF BED, then writes the label table.
#[test]
fn sdcore_needs_the_copies_table_and_a_sedef_bed() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("sdcore");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py, fa) = (dir.join("run"), dir.join("py.sh"), dir.join("g.fa"));
    // stand-in for sdcore.py COPIES SEDEF OUT
    std::fs::write(&py, "#!/bin/bash\necho labels > \"$4\"; echo done\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}'\nOUT={out}; FASTA={fa}; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\nstage_sdcore\n",
        function_block("stage_sdcore"),
        py = py.display(),
        out = out.display(),
        fa = fa.display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("sdcore needs"));
    std::fs::write(format!("{}.fam.copies.tsv", out.display()), "x\n").unwrap();
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("SEDEF BED"));
    std::fs::write(format!("{}.sedef.bed", fa.display()), "x\n").unwrap();
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(format!("{}.fam.sdcore.tsv", out.display())).unwrap(), "labels\n");
}

/// Opt-in `exonmap` (register 1595): refuses without the families output, then runs the table and the graph view.
#[test]
fn exonmap_needs_the_families_output_and_then_writes_table_and_graphs() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("exonmap");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py) = (dir.join("run"), dir.join("py.sh"));
    // stand-in for exonmap.py: `gfa PREFIX COLUMNS DIR` makes DIR, the table call writes OUT_PREFIX.columns.tsv
    std::fs::write(&py, "#!/bin/bash\nif [ \"$2\" = gfa ]; then mkdir \"$5\"; echo graphs; else echo table > \"$5.columns.tsv\"; fi\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}'\nOUT={out}; THREADS=2; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\nstage_exonmap\n",
        function_block("stage_exonmap"),
        py = py.display(),
        out = out.display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("exonmap needs"));
    for f in ["clusters.tsv", "loci.gff3", "loci.gated.paf"] {
        std::fs::write(format!("{}.fam.{f}", out.display()), "x\n").unwrap();
    }
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(format!("{}.fam.exonmap.columns.tsv", out.display())).unwrap(), "table\n");
    assert!(std::path::Path::new(&format!("{}.fam.exonmap.gfa", out.display())).is_dir());
}

/// Opt-in `repeatloci` (registers 1698 / 1702): refuses without the assembly, then writes PREFIX.repeat_loci.tsv only.
#[test]
fn repeatloci_needs_the_assembly_and_then_writes_the_table() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("repeatloci");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py, inp) = (dir.join("run"), dir.join("py.sh"), dir.join("in"));
    std::fs::write(&py, "#!/bin/bash\necho table > \"$6.tsv\"\necho seq > \"$6.fa\"\necho done\n").unwrap();
    std::fs::write(&inp, "x\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}' RUSTLE_REPEAT_FASTA={inp}\nOUT={out}; BAM={inp}; FASTA={inp}; THREADS=2; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\nstage_repeatloci\n",
        function_block("stage_repeatloci"),
        py = py.display(),
        out = out.display(),
        inp = inp.display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("repeatloci needs"));
    std::fs::write(format!("{}.gtf", out.display()), "model\n").unwrap();
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(format!("{}.repeat_loci.tsv", out.display())).unwrap(), "table\n");
    assert!(!std::path::Path::new(&format!("{}.repeat_loci.tmp.fa", out.display())).exists());
    assert_eq!(std::fs::read_to_string(format!("{}.gtf", out.display())).unwrap(), "model\n");
}

/// Opt-in `trunc5` (register 1677): refuses without the assembly, then writes PREFIX.trunc5.gtf and leaves PREFIX.gtf alone.
#[test]
fn trunc5_needs_the_assembly_and_then_writes_the_trunc5_gtf() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("trunc5");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py, inp) = (dir.join("run"), dir.join("py.sh"), dir.join("in"));
    std::fs::write(&py, "#!/bin/bash\necho annotated > \"$4\"\necho done\n").unwrap();
    std::fs::write(&inp, "x\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}'\nOUT={out}; BAM={inp}; FASTA={inp}; THREADS=2; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\nstage_trunc5\n",
        function_block("stage_trunc5"),
        py = py.display(),
        out = out.display(),
        inp = inp.display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("trunc5 needs"));
    let gtf = format!("{}.gtf", out.display());
    std::fs::write(&gtf, "model\n").unwrap();
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(format!("{}.trunc5.gtf", out.display())).unwrap(), "annotated\n");
    assert_eq!(std::fs::read_to_string(&gtf).unwrap(), "model\n");
}

/// Opt-in `chainscore` (register 1761): refuses without the assembly, then writes PREFIX.chainscore.gtf and leaves PREFIX.gtf alone.
#[test]
fn chainscore_needs_the_assembly_and_then_writes_the_chainscore_gtf() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("chainscore");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py, inp) = (dir.join("run"), dir.join("py.sh"), dir.join("in"));
    std::fs::write(&py, "#!/bin/bash\necho annotated > \"$7\"\necho done\n").unwrap();
    std::fs::write(&inp, "x\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}'\nOUT={out}; BAM={inp}; FASTA={inp}; THREADS=2; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\n{}\nstage_chainscore\n",
        function_block("chainscore_model"),
        function_block("stage_chainscore"),
        py = py.display(),
        out = out.display(),
        inp = inp.display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("chainscore needs"));
    let gtf = format!("{}.gtf", out.display());
    std::fs::write(&gtf, "model\n").unwrap();
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(format!("{}.chainscore.gtf", out.display())).unwrap(), "annotated\n");
    assert_eq!(std::fs::read_to_string(&gtf).unwrap(), "model\n");
}

/// Opt-in `shortexons` (register 1753 / 1754 / 1757): refuses without the reads, then rewrites them, assembles the rewritten BAM through the
/// driver itself and leaves PREFIX.gtf alone; the rewritten BAM is removed.
#[test]
fn shortexons_needs_the_reads_and_then_writes_the_shortexons_gtf() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("shortexons");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py, me, inp) = (dir.join("run"), dir.join("py.sh"), dir.join("driver.sh"), dir.join("in"));
    std::fs::write(&py, "#!/bin/bash\necho reads > \"$4\"\necho 3 units\n").unwrap();
    // the stub driver: `assemble --bam B --fasta G --out P ..` writes P.gtf from the BAM it was given
    std::fs::write(&me, "#!/bin/bash\n[ \"$1\" = assemble ] && cat \"$3\" > \"$7.gtf\"\n").unwrap();
    let gtf = format!("{}.gtf", out.display());
    std::fs::write(&gtf, "model\n").unwrap();
    let script = |bam: &str| {
        format!(
            "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}'\nOUT={out}; BAM={bam}; FASTA={inp}; THREADS=2; BIN=bin; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\nstage_shortexons\n",
            function_block("stage_shortexons"),
            py = py.display(),
            out = out.display(),
            inp = inp.display()
        )
    };
    let run = |bam: &str| std::process::Command::new("bash").arg("-c").arg(script(bam)).arg(&me).output().expect("run bash");
    std::fs::write(&inp, "x\n").unwrap();
    let o = run(&dir.join("missing").display().to_string());
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("shortexons needs"));
    let o = run(&inp.display().to_string());
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(format!("{}.shortexons.gtf", out.display())).unwrap(), "reads\n");
    assert_eq!(std::fs::read_to_string(&gtf).unwrap(), "model\n");
    assert!(!PathBuf::from(format!("{}.shortexons.d", out.display())).exists());
}

/// Opt-in `shiftdrop` (register 1727 / 1730 / 1732): refuses without the assembly, then writes PREFIX.shiftdrop.gtf and leaves PREFIX.gtf alone.
#[test]
fn shiftdrop_needs_the_assembly_and_then_writes_the_shiftdrop_gtf() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("shiftdrop");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py, inp) = (dir.join("run"), dir.join("py.sh"), dir.join("in"));
    std::fs::write(&py, "#!/bin/bash\necho annotated > \"$4\"\necho done\n").unwrap();
    std::fs::write(&inp, "x\n").unwrap();
    let script = format!(
        "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}'\nOUT={out}; BAM={inp}; FASTA={inp}; THREADS=2; RP_DRIVER_DIR=tools\nsay() {{ :; }}\n{}\nstage_shiftdrop\n",
        function_block("stage_shiftdrop"),
        py = py.display(),
        out = out.display(),
        inp = inp.display()
    );
    let run = || std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
    let o = run();
    assert_eq!(o.status.code(), Some(2), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("shiftdrop needs"));
    let gtf = format!("{}.gtf", out.display());
    std::fs::write(&gtf, "model\n").unwrap();
    let o = run();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read_to_string(format!("{}.shiftdrop.gtf", out.display())).unwrap(), "annotated\n");
    assert_eq!(std::fs::read_to_string(&gtf).unwrap(), "model\n");
}

/// Register 1761, default since 2026-11-16: the layered assembly is filtered by the chain score; RUSTLE_CHAINSCORE=0 turns it
/// off, and a failing filter leaves PREFIX.final.gtf as the layered run wrote it.
#[test]
fn the_final_gtf_is_filtered_by_the_chain_score_unless_switched_off_or_failing() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_dispatch").join("final_chainscore");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let (out, py) = (dir.join("run"), dir.join("py.sh"));
    // stand-in for chainscore.py apply MODEL GTF BAM FASTA OUT THREADS GTF2 OUT2
    std::fs::write(&py, "#!/bin/bash\n[ -z \"${FAIL:-}\" ] || exit 1\necho \"$3\" > \"$7.model\"; echo scored > \"$7\"; echo filtered > \"${10}\"; echo 'x: dropped 1'\n").unwrap();
    let fin = format!("{}.final.gtf", out.display());
    let run = |env: &str| {
        std::fs::write(&fin, "layered\n").unwrap();
        let script = format!(
            "set -euo pipefail\nexport RUSTLE_PYTHON='bash {py}' {env}\nOUT={out}; BAM=b; FASTA=f; THREADS=2; RP_DRIVER_DIR=tools\nsay() {{ echo \"$*\"; }}\n{}\n{}\nfinal_chainscore\n",
            function_block("chainscore_model"),
            function_block("final_chainscore"),
            py = py.display(),
            out = out.display()
        );
        let o = std::process::Command::new("bash").arg("-c").arg(&script).output().expect("run bash");
        assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
        (std::fs::read_to_string(&fin).unwrap(), String::from_utf8_lossy(&o.stdout).into_owned())
    };
    assert_eq!(run("").0, "filtered\n");
    assert_eq!(std::fs::read_to_string(format!("{}.chainscore.gtf", out.display())).unwrap(), "scored\n");
    // the weight file: r1769 (default since 2026-11-18) unless RUSTLE_CHAINSCORE_MODEL names another one (r1677, the earlier default) or gives a path
    let model = || std::fs::read_to_string(format!("{}.chainscore.gtf.tmp.model", out.display())).unwrap();
    assert_eq!(model(), "tools/data/chainscore_r1769.json\n");
    run("RUSTLE_CHAINSCORE_MODEL=r1677");
    assert_eq!(model(), "tools/data/chainscore_r1677.json\n");
    run("RUSTLE_CHAINSCORE_MODEL=/some/where/m.json");
    assert_eq!(model(), "/some/where/m.json\n");
    assert_eq!(run("RUSTLE_CHAINSCORE=0").0, "layered\n");
    let (kept, said) = run("FAIL=1");
    assert_eq!(kept, "layered\n");
    assert!(said.contains("WARNING the chain score filter failed"), "{said}");
}
