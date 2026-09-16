//! `copy_assign --families` — the O1 → O2 FILE CONTRACT.
//!
//! O1 (`gw_family_catalog`) and O2 (`copy_assign`) already shared one node type, one edge engine and one
//! admission primitive BY FUNCTION CALL, and nothing BY FILE: each binary re-derived its own families from
//! the BAM, and their family ids (`GWFAM{i}` vs `CAFAM{i}`) were assigned independently, so the two tables
//! had no join key. These tests pin the three things that fixes:
//!
//! 1. the FLAG — `--families` makes O2 consume the catalog instead of detecting;
//! 2. the JOIN KEY — the emitted rows carry the catalog's own `family_id`/`tid`;
//! 3. the LOUD-FAILURE contract — every way a supplied copy could go missing is an error, never a filter.
//!
//! All of these use the committed `same_chrom_supplement` fixture, whose `out_default.copies.tsv` /
//! `.copies.fa` are a real `gw_family_catalog` output (see `tests/gw_family_catalog_regression.rs`).

use std::path::PathBuf;
use std::process::{Command, Output};

const FIX: &str = "tests/fixtures/same_chrom_supplement";
/// The fixture's one SAME-CHROMOSOME catalog family: c1:250-460 + c1:380-550.
const GWFAM1_TSV: &str = "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n\
GWFAM1\t0\tDN_c1_250_2\tc1\t250\t460\t2\t+\t6\t250-299,371-460\n\
GWFAM1\t1\tDN_c1_380_2\tc1\t380\t550\t2\t+\t3\t380-419,451-550\n";

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("copy_assign_families").join(name);
    std::fs::create_dir_all(&d).expect("create scratch dir");
    d
}

fn write(dir: &PathBuf, name: &str, body: &str) -> String {
    let p = dir.join(name);
    std::fs::write(&p, body).expect("write fixture");
    p.to_str().expect("utf-8 path").to_string()
}

/// Run `copy_assign` over the fixture region with `extra` args appended. Never asserts success — the
/// contract tests need the failing runs.
fn run(dir: &PathBuf, extra: &[&str]) -> (Output, String) {
    let out = dir.join("o");
    let out_s = out.to_str().expect("utf-8 path").to_string();
    let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
        .args(["--region", "c1:200-600", "--out", &out_s])
        .args(extra)
        .output()
        .expect("copy_assign failed to spawn");
    (o, out_s)
}

fn read(out: &str, ext: &str) -> String {
    std::fs::read_to_string(format!("{out}.{ext}")).unwrap_or_else(|e| panic!("read {out}.{ext}: {e}"))
}

/// Column `col` of every non-header row.
fn col(text: &str, col: usize) -> Vec<String> {
    text.lines().skip(1).filter(|l| !l.trim().is_empty()).map(|l| l.split('\t').nth(col).unwrap_or("").to_string()).collect()
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

// ---- 1. the FLAG -----------------------------------------------------------------------------------

/// The supplied catalog IS the copy set: both catalog copies come back, and only those.
#[test]
fn families_flag_assigns_exactly_the_supplied_copy_set() {
    let d = scratch("flag");
    let fam = write(&d, "cat.copies.tsv", GWFAM1_TSV);
    let (o, out) = run(&d, &["--families", &fam, "--copies-fa", &format!("{FIX}/out_default.copies.fa")]);
    assert!(o.status.success(), "run failed:\n{}", stderr(&o));

    let quant = read(&out, "quant.tsv");
    let mut tids = col(&quant, 2);
    tids.sort();
    assert_eq!(tids, vec!["DN_c1_250_2", "DN_c1_380_2"], "the copy set must be the supplied catalog's, exactly");

    // and family CONSTRUCTION must not have run: no detection, no refine, no rescue.
    let err = stderr(&o);
    assert!(err.contains("assigned AS GIVEN"), "expected the as-given banner in:\n{err}");
    assert!(!err.contains("[detect_and_assign] refine:"), "refine must not run on a supplied roster:\n{err}");
    let rescued: Vec<String> = col(&read(&out, "families.tsv"), 3); // rescued_copies
    assert!(rescued.iter().all(|v| v == "0"), "rescue must not widen a supplied roster: {rescued:?}");
}

/// Without `--copies-fa`, the sequences are rebuilt from `--fasta` at the catalog's own exon coordinates —
/// and that must reach the SAME copy set (this is the fallback documented on the flag).
#[test]
fn families_without_copies_fa_rebuilds_the_sequences_from_the_genome() {
    let d = scratch("rebuild");
    let fam = write(&d, "cat.copies.tsv", GWFAM1_TSV);
    let (o, out) = run(&d, &["--families", &fam]);
    assert!(o.status.success(), "run failed:\n{}", stderr(&o));
    assert!(stderr(&o).contains("rebuilt at the catalog's exon coordinates"), "{}", stderr(&o));
    let mut tids = col(&read(&out, "quant.tsv"), 2);
    tids.sort();
    assert_eq!(tids, vec!["DN_c1_250_2", "DN_c1_380_2"]);
}

/// `--discover-copies` is opt-in and REPORT ONLY (Task 5): with the flag unset, the binary must never
/// even know the feature exists -- every other output file this invocation unconditionally produces
/// (`assignments.tsv`, `families.tsv`, `quant.tsv`, plus the two always-written files `famcn_readonly.tsv`
/// and `params.tsv`, see `src/bin/copy_assign.rs:4286` and `:4839`) must come out byte-for-byte identical
/// to a run with the flag added. This is the single most important untested claim from the copy-discovery
/// feature itself (Tasks 1-4), so it is checked directly against the real binary, not the library code.
#[test]
fn discover_copies_off_by_default_is_byte_identical() {
    // Two separate scratch dirs give the two runs distinct --out prefixes (the `run()` helper always
    // writes to `<dir>/o`), mirroring `without_families_ids_are_minted_and_no_join_file_is_written`'s
    // plain (non-`--families`) invocation above.
    let d_off = scratch("discover_off");
    let (o_off, out_off) = run(&d_off, &["--no-refine"]);
    assert!(o_off.status.success(), "flag-off run failed:\n{}", stderr(&o_off));

    let d_on = scratch("discover_on");
    let (o_on, out_on) = run(&d_on, &["--no-refine", "--discover-copies"]);
    assert!(o_on.status.success(), "flag-on run failed:\n{}", stderr(&o_on));

    for ext in ["assignments.tsv", "families.tsv", "quant.tsv", "famcn_readonly.tsv", "params.tsv"] {
        let off_path = format!("{out_off}.{ext}");
        let on_path = format!("{out_on}.{ext}");
        let a = std::fs::read(&off_path).unwrap_or_else(|e| panic!("read {off_path}: {e}"));
        let b = std::fs::read(&on_path).unwrap_or_else(|e| panic!("read {on_path}: {e}"));
        assert_eq!(a, b, "--discover-copies must not perturb .{ext} (unset vs set)");
    }

    // The flag's own report is additive, not a rename of an existing file: absent when unset, present
    // (even if empty of candidate rows) when set.
    assert!(
        std::fs::metadata(format!("{out_off}.discovered_copies.tsv")).is_err(),
        "discovered_copies.tsv must not exist without --discover-copies"
    );
    assert!(
        std::fs::metadata(format!("{out_on}.discovered_copies.tsv")).is_ok(),
        "discovered_copies.tsv must exist with --discover-copies"
    );
}

// ---- 2. the JOIN KEY -------------------------------------------------------------------------------

/// The whole point: rows must carry the CATALOG's `family_id`, and `<out>.family_join.tsv` must name the
/// catalog row behind every assigned copy.
#[test]
fn families_flag_emits_the_catalog_id_as_the_join_key() {
    let d = scratch("join");
    let fam = write(&d, "cat.copies.tsv", GWFAM1_TSV);
    let (o, out) = run(&d, &["--families", &fam, "--copies-fa", &format!("{FIX}/out_default.copies.fa")]);
    assert!(o.status.success(), "run failed:\n{}", stderr(&o));

    let fams = read(&out, "families.tsv");
    assert_eq!(col(&fams, 0), vec!["GWFAM1"], "family_id must BE the catalog id, not a minted CAFAM id");

    let join = read(&out, "family_join.tsv");
    assert!(join.starts_with("family_id\tcopy_index\tcopy_tid\tcatalog_family_id\tcatalog_copy_idx\t"), "{join}");
    let rows: Vec<Vec<&str>> = join.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 2, "one join row per assigned copy: {join}");
    for r in &rows {
        assert_eq!(r[3], "GWFAM1", "catalog_family_id");
        // the join must be able to look the copy back up in the catalog table it came from
        let want = format!("GWFAM1\t{}\t{}\t", r[4], r[2]);
        assert!(GWFAM1_TSV.contains(&want), "row {r:?} does not join back to copies.tsv (looked for {want:?})");
    }
}

/// Control: WITHOUT `--families` the binary still mints its own ids and writes no join file — i.e. the
/// join key is what the flag ADDS, not something the test would have seen anyway.
#[test]
fn without_families_ids_are_minted_and_no_join_file_is_written() {
    let d = scratch("noflag");
    let (o, out) = run(&d, &["--no-refine"]);
    assert!(o.status.success(), "run failed:\n{}", stderr(&o));
    assert!(
        std::fs::metadata(format!("{out}.family_join.tsv")).is_err(),
        "the join file must only exist under --families"
    );
    for id in col(&read(&out, "families.tsv"), 0) {
        assert!(id.starts_with("CAFAM") || id.starts_with("DSFAM") || id.starts_with("TSFAM"), "{id}");
    }
}

// ---- 3. the LOUD-FAILURE contract ------------------------------------------------------------------

/// 2026-09-15: a cross-chromosome catalog family is no longer refused. `copy_assign` gathers its reads
/// directly from every one of its copies' own chromosomes and pools them, so the AS-tied certificate
/// compares a read against the family's FULL copy set — never one truncated to a single region's contig.
/// The fixture's `read_cross_0/1/2` are built exactly for this: identical sequence, identical AS score
/// (100), one placed at MAPQ 60 on `c1:1` and a second, equally-scoring placement at MAPQ 0 on `c2:1` — a
/// real tie ACROSS chromosomes that only a cross-chromosome-aware certificate can see as one molecule.
#[test]
fn a_cross_chrom_family_is_assigned_not_refused() {
    let d = scratch("xchrom");
    // GWFAM0 of the committed catalog: c1:0-260 + c2:0-260.
    let tsv = std::fs::read_to_string(format!("{FIX}/out_default.copies.tsv")).unwrap();
    let only0: String = tsv
        .lines()
        .filter(|l| l.starts_with("family_id") || l.starts_with("GWFAM0\t"))
        .collect::<Vec<_>>()
        .join("\n");
    let fam = write(&d, "x.copies.tsv", &format!("{only0}\n"));
    // The committed `run()` helper hardcodes a c1-only `--region`; this family also needs c2 covered, so
    // sweep both contigs directly via `--regions`.
    let regions = write(&d, "regions.txt", "c1:0-600\nc2:0-320\n");
    let out = d.join("o");
    let out_s = out.to_str().expect("utf-8 path").to_string();
    let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
        .args(["--regions", &regions, "--out", &out_s])
        .args(["--families", &fam, "--copies-fa", &format!("{FIX}/out_default.copies.fa")])
        .output()
        .expect("copy_assign failed to spawn");
    assert!(o.status.success(), "a cross-chrom family must now be assignable:\n{}", stderr(&o));
    let e = stderr(&o);
    assert!(e.contains("spans 2 chromosomes"), "expected the cross-chrom banner in:\n{e}");
    assert!(e.contains("cross-chromosome pass"), "{e}");

    let quant = read(&out_s, "quant.tsv");
    let mut tids = col(&quant, 2);
    tids.sort();
    assert_eq!(tids, vec!["DN_c1_0_2", "DN_c2_0_2"], "both cross-chrom copies must be assignable: {quant}");

    // The property this fix actually guarantees: read_cross_0/1/2's records on BOTH c1 and c2 survive the
    // per-region read-gathering and dedup (before the fix a same-name, same-offset record on a SECOND
    // chromosome was silently collapsed onto the first one — the exact truncation this feature exists to
    // avoid) and the run completes rather than refusing the family outright.
    let assignments = read(&out_s, "assignments.tsv");
    let cross_rows: Vec<&str> = assignments.lines().filter(|l| l.contains("read_cross_")).collect();
    assert_eq!(cross_rows.len(), 3, "all 3 read_cross_* molecules must appear in the assignment output:\n{assignments}");
    // ⚠ NOT ASSERTED HERE (a known, separate limitation, not something this fix touches): whether a
    // molecule genuinely AS-tied ACROSS two chromosomes is scored as `n_candidates == 2` by the deeper
    // PSV/mosaic certificate (`assign_family_detailed_once` / `best_overlap_copy` in
    // `copy_assign_pipeline.rs`) depends on `AlignedRead`, which carries NO chromosome field at all — its
    // overlap math compares bare numeric ranges. For a family whose copies sit on DIFFERENT chromosomes
    // but at OVERLAPPING numeric coordinates (as `c1:0-260` and `c2:0-260` deliberately do here), that
    // layer can pick the wrong "best overlap" copy for a record having nothing to do with its real
    // chromosome. This is safe for a family whose copies' coordinates do not numerically coincide across
    // chromosomes (checked by hand for the real target this feature was built for), but is not a general
    // guarantee — fixing it means threading chromosome through `AlignedRead` and every PSV/mosaic call
    // site that compares positions, a much larger change than the read-gathering fix here.
}

/// A supplied copy outside every swept region would never have its reads read. Loud, not skipped.
#[test]
fn a_copy_outside_the_swept_regions_aborts() {
    let d = scratch("outside");
    let tsv = "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n\
GWFAM1\t0\tDN_c1_0_1\tc1\t0\t60\t1\t+\t6\t0-60\n\
GWFAM1\t1\tDN_c1_380_2\tc1\t380\t550\t2\t+\t3\t380-419,451-550\n";
    let fam = write(&d, "o.copies.tsv", tsv);
    let (o, _out) = run(&d, &["--families", &fam]); // region is c1:200-600; the family starts at 0
    assert!(!o.status.success(), "a family outside the swept regions must abort");
    assert!(stderr(&o).contains("lies outside every --region"), "{}", stderr(&o));
}

/// A supplied copy with no reads in the region cannot be assigned; dropping it silently would understate
/// the family's copy count and loosen the Bonferroni certificate over the survivors.
#[test]
fn a_copy_with_no_reads_aborts_rather_than_being_dropped() {
    let d = scratch("noreads");
    // c1 is 600 bp and the fixture's reads stop at 550: nothing overlaps 555-600.
    let tsv = "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n\
GWFAM7\t0\tDN_c1_555_1\tc1\t555\t575\t1\t+\t3\t555-575\n\
GWFAM7\t1\tDN_c1_578_1\tc1\t578\t598\t1\t+\t3\t578-598\n";
    let fam = write(&d, "n.copies.tsv", tsv);
    let (o, _out) = run(&d, &["--families", &fam]);
    assert!(!o.status.success(), "a read-less supplied copy must abort");
    let e = stderr(&o);
    assert!(e.contains("has NO reads"), "{e}");
    assert!(e.contains("subset BAM"), "the message should name the recurring cause: {e}");
}

/// `--copies-fa` that does not cover every supplied copy is a mismatched pair of files, not a licence to
/// fall back to the genome for the missing ones.
#[test]
fn a_copies_fa_missing_a_record_aborts() {
    let d = scratch("nofa");
    let fam = write(&d, "cat.copies.tsv", GWFAM1_TSV);
    // only the FIRST of GWFAM1's two copies
    let fa = write(&d, "partial.copies.fa", ">GWFAM1|0|c1:250-460|+|nexon=2\nACGTACGTAC\n");
    let (o, _out) = run(&d, &["--families", &fam, "--copies-fa", &fa]);
    assert!(!o.status.success(), "an incomplete --copies-fa must abort");
    assert!(stderr(&o).contains("has no record for GWFAM1 copy 1"), "{}", stderr(&o));
}

/// A `--copies-fa` record whose header disagrees with its `copies.tsv` row means the two files are from
/// different runs. Checked, not trusted.
#[test]
fn a_copies_fa_record_disagreeing_with_the_tsv_aborts() {
    let d = scratch("fadisagree");
    let fam = write(&d, "cat.copies.tsv", GWFAM1_TSV);
    let fa = write(
        &d,
        "wrong.copies.fa",
        ">GWFAM1|0|c1:250-461|+|nexon=2\nACGT\n>GWFAM1|1|c1:380-550|+|nexon=2\nACGT\n",
    );
    let (o, _out) = run(&d, &["--families", &fam, "--copies-fa", &fa]);
    assert!(!o.status.success());
    assert!(stderr(&o).contains("but copies.tsv says"), "{}", stderr(&o));
}

/// Flags that would CHANGE the roster are refused up front, so `--families` can never quietly assign
/// against a copy set that is not the supplied one.
#[test]
fn roster_changing_flags_are_refused_under_families() {
    let d = scratch("incompat");
    let fam = write(&d, "cat.copies.tsv", GWFAM1_TSV);
    for flag in ["--absent-copies", "--vg-realign", "--iterative-prune", "--collapse-gate", "--tied-seed", "--recover-copies"] {
        let (o, _out) = run(&d, &["--families", &fam, flag]);
        assert!(!o.status.success(), "{flag} must be refused under --families");
        assert!(stderr(&o).contains(&format!("--families is incompatible with {flag}")), "{}", stderr(&o));
    }
}

/// `--copies-fa` without `--families` is a user error, not a silent no-op.
#[test]
fn copies_fa_without_families_is_an_error() {
    let d = scratch("orphanfa");
    let (o, _out) = run(&d, &["--copies-fa", &format!("{FIX}/out_default.copies.fa")]);
    assert!(!o.status.success());
    assert!(stderr(&o).contains("only meaningful with --families"), "{}", stderr(&o));
}
