//! `tools/rustle_pipeline.sh`: which alignments seed the assembly (the seeding pool), and what the switch passes down.
//!
//! The whole driver runs its `assemble` stage against stub `copy_assign` / `as_table` binaries (`--bin`) that record the
//! three environment variables the Rust side reads (`RUSTLE_GTF_SECONDARY`, `RUSTLE_GTF_SECONDARY_AS_RATIO`,
//! `RUSTLE_GTF_SECONDARY_AS_TABLE`) and the argument list. What is checked is the decision the driver makes, so nothing here
//! needs a genome or a BAM.
//!
//! 2026-10-07: the pool became a switch, `--seed-pool primary|good|all` (+ `--seed-as-ratio R` for `good`), so the effect of
//! the choice can be shown on one command line. Before, `good` was the only non-primary pool the driver could ask for and
//! `all` needed a hand-written `RUSTLE_GTF_SECONDARY=1` call of `copy_assign`. The default (`good`, ratio 0.98, the AS
//! table) and both legacy flags must keep their exact meaning.

use std::path::{Path, PathBuf};
use std::process::Command;

const DRIVER: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tools/rustle_pipeline.sh");

const COPY_ASSIGN_STUB: &str = r#"#!/bin/bash
out=""; args=("$@")
while [ $# -gt 0 ]; do case "$1" in --out) out=$2; shift 2;; *) shift;; esac; done
printf 'copy_assign sec=%s ratio=%s table=%s\n' "${RUSTLE_GTF_SECONDARY-unset}" "${RUSTLE_GTF_SECONDARY_AS_RATIO-unset}" "${RUSTLE_GTF_SECONDARY_AS_TABLE-unset}" >> "$STUB_LOG"
printf 'argv %s\n' "${args[*]}" >> "$STUB_LOG"
printf 'chr1\tstub\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' > "$out.gtf"
"#;

const AS_TABLE_STUB: &str = r#"#!/bin/bash
bam=""; out=""
while [ $# -gt 0 ]; do case "$1" in --bam) bam=$2; shift 2;; --out) out=$2; shift 2;; *) shift;; esac; done
echo as_table >> "$STUB_LOG"
printf 'bam=%s\tstub\n' "$(readlink -f "$bam")" > "$out"
"#;

struct Run {
    code: i32,
    stdout: String,
    stderr: String,
    /// the stub's records, one per line
    log: Vec<String>,
    /// the three seeding variables and nothing else, as the stub `copy_assign` saw them
    seen: String,
    argv: String,
    /// whether `as_table` ran
    as_table_ran: bool,
}

fn write_exec(path: &Path, text: &str) {
    std::fs::write(path, text).expect("write stub");
    use std::os::unix::fs::PermissionsExt;
    std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).expect("chmod stub");
}

/// The FASTA index every case gets unless it asks for none: two contigs, the second with the length of CHM13 chr16.
const FAI: &str = "chr1\t248387328\t7\t80\t81\nchr16\t96330374\t252000000\t80\t81\n";

/// The scratch directory of a case (the helper empties it first).
fn case_dir(tag: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_seed_pool").join(tag)
}

/// Run `assemble` with `flags`. `env` is added to an otherwise empty environment (PATH only), so an exported
/// `RUSTLE_*` of the machine cannot reach the driver unless the test puts it there.
fn assemble(tag: &str, flags: &[&str], env: &[(&str, &str)]) -> Run {
    assemble_with(tag, flags, env, Some(FAI), &|_| {})
}

/// `assemble` with the FASTA index (None: no index at all) and a hook that runs once the case directory exists.
fn assemble_with(tag: &str, flags: &[&str], env: &[(&str, &str)], fai: Option<&str>, setup: &dyn Fn(&Path)) -> Run {
    let dir = case_dir(tag);
    let _ = std::fs::remove_dir_all(&dir);
    let bin = dir.join("bin");
    std::fs::create_dir_all(&bin).expect("create scratch dir");
    write_exec(&bin.join("copy_assign"), COPY_ASSIGN_STUB);
    write_exec(&bin.join("as_table"), AS_TABLE_STUB);
    let (bam, fasta, out) = (dir.join("reads.bam"), dir.join("g.fa"), dir.join("run"));
    std::fs::write(&bam, "").unwrap();
    std::fs::write(&fasta, "").unwrap();
    if let Some(fai) = fai {
        std::fs::write(dir.join("g.fa.fai"), fai).unwrap();
    }
    setup(&dir);
    let stub_log = dir.join("stub.log");
    let mut cmd = Command::new("bash");
    cmd.arg(DRIVER)
        .arg("assemble")
        .args(["--bam", bam.to_str().unwrap(), "--fasta", fasta.to_str().unwrap(), "--out", out.to_str().unwrap()])
        .args(["--bin", bin.to_str().unwrap(), "--no-cache"])
        .args(flags)
        .env_clear()
        .env("PATH", std::env::var("PATH").unwrap_or_default())
        .env("TMPDIR", &dir)
        .env("STUB_LOG", &stub_log);
    for (k, v) in env {
        cmd.env(k, v);
    }
    let o = cmd.output().expect("bash failed to spawn");
    let log: Vec<String> = std::fs::read_to_string(&stub_log).unwrap_or_default().lines().map(String::from).collect();
    let seen = log.iter().find(|l| l.starts_with("copy_assign ")).cloned().unwrap_or_default();
    let argv = log.iter().find(|l| l.starts_with("argv ")).cloned().unwrap_or_default();
    Run {
        code: o.status.code().unwrap_or(-1),
        stdout: String::from_utf8_lossy(&o.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&o.stderr).into_owned(),
        as_table_ran: log.iter().any(|l| l == "as_table"),
        log,
        seen: seen.replace(&dir.to_string_lossy().to_string(), "<dir>"),
        argv: argv.replace(&dir.to_string_lossy().to_string(), "<dir>"),
    }
}

const GOOD: &str = "copy_assign sec=1 ratio=0.98 table=<dir>/run.molecules.tsv";
const PRIMARY: &str = "copy_assign sec=unset ratio=unset table=unset";
const ALL: &str = "copy_assign sec=1 ratio=unset table=unset";

// ---- the default and the legacy flags keep their meaning (these passed before the switch existed) ----

#[test]
fn the_default_seeds_with_good_secondaries_ratio_098_and_the_as_table() {
    let r = assemble("default", &[], &[]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, GOOD);
    assert!(r.as_table_ran, "the good pool needs the genome-wide best-AS table");
}

#[test]
fn no_seed_secondaries_is_the_primary_pool() {
    let r = assemble("legacy_primary", &["--no-seed-secondaries"], &[]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, PRIMARY);
    assert!(!r.as_table_ran, "the primary pool needs no AS table");
}

#[test]
fn seed_secondaries_is_the_good_pool() {
    let r = assemble("legacy_good", &["--no-seed-secondaries", "--seed-secondaries"], &[]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, GOOD);
}

// ---- the switch ----

#[test]
fn seed_pool_primary_sets_nothing_and_skips_the_as_table() {
    let r = assemble("pool_primary", &["--seed-pool", "primary"], &[]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, PRIMARY);
    assert!(!r.as_table_ran);
}

#[test]
fn seed_pool_all_admits_every_secondary_with_no_ratio_and_no_table() {
    let r = assemble("pool_all", &["--seed-pool", "all"], &[]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, ALL);
    assert!(!r.as_table_ran, "the all pool filters nothing, so it needs no AS table");
}

#[test]
fn seed_pool_good_is_the_default_spelled_out() {
    let r = assemble("pool_good", &["--seed-pool", "good"], &[]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, GOOD);
    assert!(r.as_table_ran);
}

#[test]
fn seed_as_ratio_sets_the_width_of_the_good_pool() {
    let r = assemble("ratio", &["--seed-pool", "good", "--seed-as-ratio", "0.90"], &[]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, "copy_assign sec=1 ratio=0.90 table=<dir>/run.molecules.tsv");
    let r1 = assemble("ratio_one", &["--seed-as-ratio", "1"], &[]);
    assert_eq!(r1.code, 0, "{}", r1.stderr);
    assert_eq!(r1.seen, "copy_assign sec=1 ratio=1 table=<dir>/run.molecules.tsv", "the ratio alone keeps the default pool (good)");
}

#[test]
fn only_the_environment_the_pool_needs_differs_the_argument_list_is_the_same() {
    let base = assemble("argv_default", &[], &[]).argv;
    assert!(base.contains("--assemble-only") && base.contains("--bridge-regroup"), "{base}");
    for (tag, flags) in [
        ("argv_primary", vec!["--seed-pool", "primary"]),
        ("argv_all", vec!["--seed-pool", "all"]),
        ("argv_ratio", vec!["--seed-as-ratio", "0.9"]),
        ("argv_legacy", vec!["--no-seed-secondaries"]),
    ] {
        assert_eq!(assemble(tag, &flags, &[]).argv, base, "{flags:?}");
    }
}

#[test]
fn an_exported_secondary_variable_cannot_leak_into_the_primary_pool() {
    let polluted = [
        ("RUSTLE_GTF_SECONDARY", "1"),
        ("RUSTLE_GTF_SECONDARY_AS_RATIO", "0.5"),
        ("RUSTLE_GTF_SECONDARY_AS_TABLE", "/nonexistent.tsv"),
    ];
    let r = assemble("leak_primary", &["--seed-pool", "primary"], &polluted);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, PRIMARY, "the pool the flag names is the pool the binary gets");
    let r = assemble("leak_all", &["--seed-pool", "all"], &polluted);
    assert_eq!(r.seen, ALL, "an exported ratio or table must not narrow the all pool");
}

#[test]
fn the_pool_can_come_from_the_environment_and_a_flag_wins() {
    let r = assemble("env_pool", &[], &[("RUSTLE_SEED_POOL", "all")]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, ALL);
    let r = assemble("flag_over_env", &["--seed-pool", "primary"], &[("RUSTLE_SEED_POOL", "all")]);
    assert_eq!(r.seen, PRIMARY);
    let r = assemble("env_ratio", &[], &[("RUSTLE_SEED_AS_RATIO", "0.95")]);
    assert_eq!(r.seen, "copy_assign sec=1 ratio=0.95 table=<dir>/run.molecules.tsv");
    // an exported ratio is a default for the good pool, not an error for the others
    let r = assemble("env_ratio_primary", &["--seed-pool", "primary"], &[("RUSTLE_SEED_AS_RATIO", "0.95")]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert_eq!(r.seen, PRIMARY);
}

#[test]
fn the_pool_is_named_in_the_log() {
    for (tag, flags, want) in [
        ("log_default", vec![], "seeding pool: good (AS >= 0.98 x the molecule's genome-wide best)"),
        ("log_primary", vec!["--seed-pool", "primary"], "seeding pool: primary"),
        ("log_all", vec!["--seed-pool", "all"], "seeding pool: all"),
        ("log_ratio", vec!["--seed-as-ratio", "0.9"], "seeding pool: good (AS >= 0.9 x the molecule's genome-wide best)"),
    ] {
        let r = assemble(tag, &flags, &[]);
        assert!(r.stderr.contains(want), "{tag}: wanted '{want}' in\n{}", r.stderr);
    }
}

// ---- refusals: a wrong spelling must not quietly run another pool ----

#[test]
fn a_bad_pool_or_ratio_is_refused_with_status_2_and_nothing_runs() {
    let cases: [(&str, Vec<&str>, Vec<(&str, &str)>, &str); 8] = [
        ("bad_pool", vec!["--seed-pool", "secondary"], vec![], "must be primary, good or all"),
        ("bad_env_pool", vec![], vec![("RUSTLE_SEED_POOL", "Good")], "must be primary, good or all"),
        ("ratio_zero", vec!["--seed-as-ratio", "0"], vec![], "number in (0, 1]"),
        ("ratio_big", vec!["--seed-as-ratio", "1.5"], vec![], "number in (0, 1]"),
        ("ratio_text", vec!["--seed-as-ratio", "abc"], vec![], "number in (0, 1]"),
        ("ratio_env_bad", vec![], vec![("RUSTLE_SEED_AS_RATIO", "2")], "number in (0, 1]"),
        ("ratio_with_primary", vec!["--seed-pool", "primary", "--seed-as-ratio", "0.9"], vec![], "only the good pool"),
        ("ratio_with_all", vec!["--seed-as-ratio", "0.9", "--seed-pool", "all"], vec![], "only the good pool"),
    ];
    for (tag, flags, env, want) in cases {
        let r = assemble(tag, &flags, &env);
        assert_eq!(r.code, 2, "{tag}: {}", r.stderr);
        assert!(r.stderr.contains(want), "{tag}: wanted '{want}' in\n{}", r.stderr);
        assert!(r.log.is_empty(), "{tag}: ran {:?}", r.log);
        assert!(r.stdout.is_empty(), "{tag}: {}", r.stdout);
    }
}

// ---- the resume key of `merged` must tell the pools apart (and keep the default's key as it was) ----

/// The driver's own `merged_env` function: from the line that starts `merged_env()` to the first line that ends it.
fn merged_env_block() -> String {
    let text = std::fs::read_to_string(DRIVER).expect("read the driver");
    let lines: Vec<&str> = text.lines().collect();
    let start = lines.iter().position(|l| l.starts_with("merged_env()")).expect("the driver has no merged_env");
    // a one-line definition ends on its own line; a longer one at the first line that is a lone `}`
    let len = if lines[start].trim_end().ends_with('}') {
        0
    } else {
        lines[start..].iter().position(|l| *l == "}").expect("merged_env has no closing brace at column 0")
    };
    lines[start..=start + len].join("\n")
}

fn merged_env(pool: &str, ratio: &str) -> String {
    merged_env_for(pool, ratio, "")
}

fn merged_env_for(pool: &str, ratio: &str, contig: &str) -> String {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_seed_pool").join("merged_env");
    std::fs::create_dir_all(&dir).unwrap();
    let (bam, fasta) = (dir.join("reads.bam"), dir.join("g.fa"));
    std::fs::write(&bam, "").unwrap();
    std::fs::write(&fasta, "").unwrap();
    let script = format!(
        "set -euo pipefail\nBAM=$1; FASTA=$2; THREADS=4; SEED_POOL=$3; SEED_RATIO=$4; CONTIG=$5; BRIDGE_MODE=f1v2\n{}\nmerged_env\n",
        merged_env_block()
    );
    let o = Command::new("bash")
        .arg("-c")
        .arg(script)
        .arg("merged_env")
        .arg(&bam)
        .arg(&fasta)
        .arg(pool)
        .arg(ratio)
        .arg(contig)
        .output()
        .expect("bash failed to spawn");
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    String::from_utf8_lossy(&o.stdout).replace(&dir.to_string_lossy().to_string(), "<dir>")
}

#[test]
fn the_merged_resume_key_of_the_default_and_of_primary_is_unchanged() {
    // byte for byte what the driver wrote before the switch, so a finished `merged` run is still recognised
    assert_eq!(merged_env("good", "0.98"), "bam=<dir>/reads.bam\nfasta=<dir>/g.fa\nthreads=4\nseed_sec=1\nbridge=f1v2\n");
    assert_eq!(merged_env("primary", "0.98"), "bam=<dir>/reads.bam\nfasta=<dir>/g.fa\nthreads=4\nseed_sec=0\nbridge=f1v2\n");
}

#[test]
fn the_merged_resume_key_tells_all_and_other_ratios_apart() {
    let default = merged_env("good", "0.98");
    let all = merged_env("all", "0.98");
    let wide = merged_env("good", "0.90");
    assert_ne!(all, default, "a finished good-pool run must not stand in for an all-pool run");
    assert!(all.contains("seed_pool=all"), "{all}");
    assert_ne!(wide, default, "a finished 0.98 run must not stand in for a 0.90 run");
    assert!(wide.contains("seed_ratio=0.90"), "{wide}");
    assert_ne!(all, wide);
}

// ---- the scope of the assemble stage and the table of the good pool (2026-10-07, for one-contig pool comparisons) ----

#[test]
fn without_a_contig_the_assembly_is_genome_wide() {
    let r = assemble("scope_default", &[], &[]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert!(r.argv.contains(" --genome-wide ") && !r.argv.contains("--region"), "{}", r.argv);
}

#[test]
fn contig_scopes_the_assembly_to_one_region_from_the_fasta_index() {
    let r = assemble("scope_contig", &["--contig", "chr16"], &[]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert!(r.argv.contains(" --region chr16:0-96330374 ") && !r.argv.contains("--genome-wide"), "{}", r.argv);
    assert_eq!(r.seen, GOOD, "the scope does not change the pool");
    assert!(r.stderr.contains("contig chr16 only"), "{}", r.stderr);
}

#[test]
fn an_unknown_contig_or_a_missing_fasta_index_is_refused_before_anything_runs() {
    let r = assemble("scope_unknown", &["--contig", "chrZ"], &[]);
    assert_eq!(r.code, 2, "{}", r.stderr);
    assert!(r.stderr.contains("chrZ") && r.stderr.contains("g.fa.fai"), "{}", r.stderr);
    assert!(r.log.is_empty(), "ran {:?}", r.log);
    let r = assemble_with("scope_no_fai", &["--contig", "chr16"], &[], None, &|_| {});
    assert_eq!(r.code, 2, "{}", r.stderr);
    assert!(r.stderr.contains("g.fa.fai"), "{}", r.stderr);
    assert!(r.log.is_empty(), "ran {:?}", r.log);
}

#[test]
fn as_table_names_the_table_the_good_pool_reads_and_builds_it_there_when_absent() {
    let table = case_dir("table_absent").join("shared").join("mol.tsv");
    let t = table.to_str().unwrap().to_string();
    let r = assemble_with("table_absent", &["--as-table", &t], &[], Some(FAI), &|d| std::fs::create_dir_all(d.join("shared")).unwrap());
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert!(r.as_table_ran, "absent table: build it");
    assert_eq!(r.seen, format!("copy_assign sec=1 ratio=0.98 table={}", t.replace(&case_dir("table_absent").to_string_lossy().to_string(), "<dir>")));
    assert!(table.exists(), "the table is written where --as-table says");
    assert!(!case_dir("table_absent").join("run.molecules.tsv").exists(), "and not at the default place");
}

#[test]
fn an_as_table_built_for_this_bam_is_reused_and_one_for_another_bam_is_rebuilt() {
    let write = |name: &'static str, bam_of: fn(&Path) -> PathBuf| {
        move |d: &Path| {
            std::fs::create_dir_all(d.join("shared")).unwrap();
            std::fs::write(d.join("shared").join(name), format!("#as_table\tbam={}\trecords=1\tmolecules=1\n", bam_of(d).display())).unwrap();
        }
    };
    let same = |d: &Path| d.join("reads.bam");
    let other = |_: &Path| PathBuf::from("/some/other.bam");
    let t = case_dir("table_same").join("shared").join("mol.tsv");
    let r = assemble_with("table_same", &["--as-table", t.to_str().unwrap()], &[], Some(FAI), &write("mol.tsv", same));
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert!(!r.as_table_ran, "a table made from this BAM is reused");
    let t = case_dir("table_other").join("shared").join("mol.tsv");
    let r = assemble_with("table_other", &["--as-table", t.to_str().unwrap()], &[], Some(FAI), &write("mol.tsv", other));
    assert!(r.as_table_ran, "a table made from another BAM is rebuilt");
}

#[test]
fn the_other_pools_never_read_or_build_a_table_even_when_one_is_named() {
    for pool in ["primary", "all"] {
        let r = assemble(&format!("table_unused_{pool}"), &["--seed-pool", pool, "--as-table", "/nonexistent/mol.tsv"], &[]);
        assert_eq!(r.code, 0, "{pool}: {}", r.stderr);
        assert!(!r.as_table_ran, "{pool}");
    }
}

#[test]
fn the_merged_resume_key_records_the_contig() {
    let whole = merged_env_for("good", "0.98", "");
    assert_eq!(whole, merged_env("good", "0.98"));
    let one = merged_env_for("good", "0.98", "chr16");
    assert_ne!(one, whole, "a finished genome-wide run must not stand in for a one-contig run");
    assert!(one.contains("contig=chr16"), "{one}");
}

#[test]
fn several_contigs_become_a_regions_file_one_line_per_contig() {
    let r = assemble_with("scope_several", &["--contig", "chr1,chr16"], &[], Some(FAI), &|_| {});
    assert_eq!(r.code, 0, "{}", r.stderr);
    assert!(r.argv.contains(" --regions <dir>/run.assemble_regions.txt ") && !r.argv.contains("--region ") && !r.argv.contains("--genome-wide"), "{}", r.argv);
    let text = std::fs::read_to_string(case_dir("scope_several").join("run.assemble_regions.txt")).expect("the regions file");
    assert_eq!(text, "chr1:0-248387328\nchr16:0-96330374\n");
    assert!(r.stderr.contains("contigs chr1,chr16 only"), "{}", r.stderr);
}

#[test]
fn one_unknown_contig_in_a_list_refuses_the_whole_call() {
    let r = assemble("scope_several_bad", &["--contig", "chr1,chrZ"], &[]);
    assert_eq!(r.code, 2, "{}", r.stderr);
    assert!(r.stderr.contains("chrZ"), "{}", r.stderr);
    assert!(r.log.is_empty(), "ran {:?}", r.log);
}
