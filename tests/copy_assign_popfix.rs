//! `copy_assign --psv-population FILE [--psv-fixed-only]` end to end (spec 2026-10-12-population-fixed-psvs; final
//! review, controller ruling R30). Its own test binary: the population table is a process-wide `OnceLock` in the
//! library (each run here is a separate `copy_assign` process anyway).
//!
//! FIXTURE. `tests/fixtures/o2_origin_resolve` (the one `tests/o2_origin_resolve.rs` uses): one family `FAM1`, two
//! two-exon copies on `c1` that differ at 24 PSV positions, and the molecule `MOL_PRIMARY_WRONG` (AS-tied between
//! them) that the read-star certificate ASSIGNS to copy 1 on 15 columns. The brief's first choice,
//! `same_chrom_supplement` with the GWFAM1 catalog of `tests/copy_assign_families.rs`, makes no certificate-decided
//! Assigned call (its three molecules tie on 0 columns), so it cannot exercise a strict Assigned call; the precondition
//! below pins that this fixture can.
//!
//! What is pinned:
//! 1. baseline (`--dump-star`) has >= 1 `assigned` row with a read-star proof (PRECONDITION);
//! 2. a table that makes every read-star genome position FIXED: the non-strict run's assignments.tsv is byte-identical
//!    to the baseline and every baseline assigned certificate row is PORTABLE; the strict run (`--psv-fixed-only`)
//!    assigns exactly the baseline's certificate-decided Assigned set, and that set is not empty;
//! 3. a table with no FIXED read-star column (every row POLYMORPHIC plus one FIXED decoy inside copy 0's span at no
//!    column — a table with no FIXED row at all is refused, see 4): strict assigns nothing among tier rows;
//!    non-strict assignments are byte-identical and never PORTABLE;
//! 4. table sanity: no FIXED row at all is an error; FIXED sites on no catalog copy's span warn (non-strict) / abort
//!    (strict);
//! 5. guards: `--psv-population` without `--families`, or with `--absent-copies`, exits non-zero naming the flag.

use std::collections::{BTreeSet, HashMap};
use std::path::PathBuf;
use std::process::{Command, Output};

const FIX: &str = "tests/fixtures/o2_origin_resolve";
const POP_HEADER: &str = "chrom\tpos\tfamily_id\tcol_index\tstatus\tmethod\tn_men\tn_cover\tn_carry\tn_mismatch\n";

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join("copy_assign_popfix")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("create scratch dir");
    d
}

fn write(dir: &PathBuf, name: &str, body: &str) -> String {
    let p = dir.join(name);
    std::fs::write(&p, body).expect("write fixture");
    p.to_str().expect("utf-8 path").to_string()
}

/// Run `copy_assign` over the fixture with `extra` args appended, output prefix `<dir>/<name>`. Never asserts
/// success — the guard tests need the failing runs.
fn run_raw(dir: &PathBuf, name: &str, extra: &[&str]) -> (Output, String) {
    let out_s = dir.join(name).to_str().expect("utf-8 path").to_string();
    let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args([
            "--bam",
            &format!("{FIX}/reads.bam"),
            "--fasta",
            &format!("{FIX}/genome.fa"),
        ])
        .args(["--region", "c1:0-3500", "--out", &out_s])
        .args(extra)
        .output()
        .expect("copy_assign failed to spawn");
    (o, out_s)
}

/// `run_raw` with the fixture's catalog (`--families` + `--copies-fa`, as `tests/o2_origin_resolve.rs`) and `--dump-star`.
fn run(dir: &PathBuf, name: &str, extra: &[&str]) -> (Output, String) {
    let mut a: Vec<&str> = Vec::new();
    let fam = format!("{FIX}/copies.tsv");
    let fa = format!("{FIX}/copies.fa");
    a.extend(["--families", fam.as_str(), "--copies-fa", fa.as_str(), "--dump-star"]);
    a.extend_from_slice(extra);
    run_raw(dir, name, &a)
}

fn read(out: &str, ext: &str) -> String {
    std::fs::read_to_string(format!("{out}.{ext}")).unwrap_or_else(|e| panic!("read {out}.{ext}: {e}"))
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

/// Rows of a TSV with a header, as column-name maps.
fn rows(text: &str) -> Vec<HashMap<String, String>> {
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().expect("tsv header").split('\t').collect();
    lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| head.iter().zip(l.split('\t')).map(|(h, v)| (h.to_string(), v.to_string())).collect())
        .collect()
}

type Key = (String, String); // (read_name, family_id)

fn key(r: &HashMap<String, String>) -> Key {
    (r["read_name"].clone(), r["family_id"].clone())
}

/// (read, family, assigned copy) of every `assigned` row of an assignments.tsv whose (read, family) is in `within`.
fn assigned_set(assignments: &str, within: &BTreeSet<Key>) -> BTreeSet<(String, String, String)> {
    rows(assignments)
        .iter()
        .filter(|r| r["status"] == "assigned" && within.contains(&key(r)))
        .map(|r| (r["read_name"].clone(), r["family_id"].clone(), r["assigned_copy"].clone()))
        .collect()
}

/// Every (chrom, 0-based pos) the baseline star dump's `positions` column names: each candidate's own genome position
/// at each column, its chrom = its copy's in `family_join.tsv` (by `(family_id, copy_index)`).
fn star_sites(out: &str) -> BTreeSet<(String, u64)> {
    let chrom: HashMap<(String, String), String> = rows(&read(out, "family_join.tsv"))
        .iter()
        .map(|r| ((r["family_id"].clone(), r["copy_index"].clone()), r["chrom"].clone()))
        .collect();
    let mut sites = BTreeSet::new();
    for r in rows(&read(out, "star_reads.tsv")) {
        let cand: Vec<&str> = if r["candidates"].is_empty() { Vec::new() } else { r["candidates"].split(',').collect() };
        if r["positions"].is_empty() {
            continue;
        }
        for col in r["positions"].split(',') {
            for (k, p) in col.split('/').enumerate() {
                if p == "." {
                    continue;
                }
                let c = chrom
                    .get(&(r["family_id"].clone(), cand[k].to_string()))
                    .unwrap_or_else(|| panic!("no family_join row for {} copy {}", r["family_id"], cand[k]));
                sites.insert((c.clone(), p.parse().expect("position")));
            }
        }
    }
    sites
}

/// A `psv_population.tsv` (bench/popfix/table.py header) with every site at `status`, plus `extra` FIXED sites.
fn pop_table(sites: &BTreeSet<(String, u64)>, status: &str, extra_fixed: &[(String, u64)]) -> String {
    let mut s = POP_HEADER.to_string();
    let all = sites.iter().map(|x| (x, status)).chain(extra_fixed.iter().map(|x| (x, "FIXED")));
    for (i, ((c, p), st)) in all.enumerate() {
        s.push_str(&format!("{c}\t{p}\tFAM1\t{i}\t{st}\ttest\t0\t0\t0\t0\n"));
    }
    s
}

/// (read, family) -> its population_tier.tsv row.
fn tiers(out: &str) -> HashMap<Key, HashMap<String, String>> {
    rows(&read(out, "population_tier.tsv")).into_iter().map(|r| (key(&r), r)).collect()
}

struct Baseline {
    out: String,
    assignments: String,
    sites: BTreeSet<(String, u64)>,
    /// (read, family) of the baseline's `assigned` rows that carry a read-star proof
    assigned_cert: BTreeSet<Key>,
}

fn baseline(d: &PathBuf) -> Baseline {
    let (o, out) = run(d, "base", &[]);
    assert!(o.status.success(), "baseline run failed:\n{}", stderr(&o));
    let assignments = read(&out, "assignments.tsv");
    let proofs: BTreeSet<Key> = rows(&read(&out, "star_reads.tsv")).iter().map(key).collect();
    let assigned_cert: BTreeSet<Key> =
        rows(&assignments).iter().filter(|r| r["status"] == "assigned").map(key).filter(|k| proofs.contains(k)).collect();
    // PRECONDITION (R30): the strict arm below must be able to make an Assigned call at all
    assert!(
        !assigned_cert.is_empty(),
        "fixture precondition: no `assigned` row with a read-star proof in\n{assignments}"
    );
    let sites = star_sites(&out);
    assert!(!sites.is_empty(), "the baseline star dump carries no genome positions");
    Baseline { out, assignments, sites, assigned_cert }
}

/// 2. Every read-star position FIXED: non-strict changes no call and marks every assigned certificate row PORTABLE;
/// strict reproduces the baseline's certificate-decided Assigned set, and that set is not empty.
#[test]
fn popfix_all_fixed_table_reproduces_the_baseline_and_strict_assigns() {
    let d = scratch("all_fixed");
    let b = baseline(&d);
    let table = write(&d, "all_fixed.tsv", &pop_table(&b.sites, "FIXED", &[]));

    let (o, ns) = run(&d, "ns", &["--psv-population", &table]);
    assert!(o.status.success(), "non-strict run failed:\n{}", stderr(&o));
    assert_eq!(read(&ns, "assignments.tsv"), b.assignments, "non-strict --psv-population must not change assignments.tsv");
    let t_ns = tiers(&ns);
    for k in &b.assigned_cert {
        let t = t_ns.get(k).unwrap_or_else(|| panic!("no tier row for baseline assigned certificate row {k:?}"));
        assert_eq!(t["tier"], "PORTABLE", "{k:?}: every column FIXED, so the fixed-only call must agree: {t:?}");
    }

    let (o, st) = run(&d, "strict", &["--psv-population", &table, "--psv-fixed-only"]);
    assert!(o.status.success(), "strict run failed:\n{}", stderr(&o));
    // certificate-decided rows = rows with a tier that the sole-candidate promotion did not decide (R24: strict never
    // promotes a sole candidate, so such a row is not a certificate call in the baseline)
    let base_sole: BTreeSet<Key> =
        rows(&b.assignments).iter().filter(|r| r["sole_candidate"] == "1").map(key).collect();
    let decided: BTreeSet<Key> =
        tiers(&st).into_iter().filter(|(_, t)| t["status_fixed"] != "NA").map(|(k, _)| k).filter(|k| !base_sole.contains(k)).collect();
    let want = assigned_set(&b.assignments, &decided);
    let got = assigned_set(&read(&st, "assignments.tsv"), &decided);
    assert!(!got.is_empty(), "the strict run must make >= 1 Assigned call (the point of this test)");
    assert_eq!(got, want, "all-FIXED strict must assign exactly the baseline's certificate-decided Assigned set");
    // and every baseline assigned certificate row is among them
    for k in &b.assigned_cert {
        assert!(decided.contains(k), "baseline assigned certificate row {k:?} has no strict tier");
    }
}

/// 3. No FIXED read-star column: strict assigns nothing among tier rows; non-strict keeps every call and none is
/// PORTABLE. The table keeps every site POLYMORPHIC and adds one FIXED decoy inside copy 0's span at a position no
/// column uses, so the table passes the sanity check (4) while no column is FIXED.
#[test]
fn popfix_no_fixed_column_strict_assigns_nothing_and_nonstrict_is_never_portable() {
    let d = scratch("no_fixed");
    let b = baseline(&d);
    let join = rows(&read(&b.out, "family_join.tsv"));
    let c0 = join.iter().find(|r| r["copy_index"] == "0").expect("copy 0 in family_join");
    let (c0_chrom, c0_start, c0_end): (String, u64, u64) =
        (c0["chrom"].clone(), c0["start"].parse().unwrap(), c0["end"].parse().unwrap());
    let decoy = (c0_start..c0_end)
        .map(|p| (c0_chrom.clone(), p))
        .find(|s| !b.sites.contains(s))
        .expect("a position of copy 0 that no column uses");
    let table = write(&d, "no_fixed.tsv", &pop_table(&b.sites, "POLYMORPHIC", &[decoy]));

    let (o, ns) = run(&d, "ns", &["--psv-population", &table]);
    assert!(o.status.success(), "non-strict run failed:\n{}", stderr(&o));
    assert_eq!(read(&ns, "assignments.tsv"), b.assignments, "non-strict --psv-population must not change assignments.tsv");
    let t_ns = tiers(&ns);
    assert!(!t_ns.is_empty(), "no tier rows at all: the assertions below would be vacuous");
    for (k, t) in &t_ns {
        assert_ne!(t["tier"], "PORTABLE", "{k:?}: no FIXED column, nothing can be PORTABLE: {t:?}");
    }

    let (o, st) = run(&d, "strict", &["--psv-population", &table, "--psv-fixed-only"]);
    assert!(o.status.success(), "strict run failed:\n{}", stderr(&o));
    let with_tier: BTreeSet<Key> = tiers(&st).into_keys().collect();
    assert!(with_tier.is_superset(&b.assigned_cert), "every baseline assigned certificate row must have a strict tier");
    let got = assigned_set(&read(&st, "assignments.tsv"), &with_tier);
    assert!(got.is_empty(), "no FIXED column: strict must not assign any certificate row, got {got:?}");
}

/// 4. Table sanity: no FIXED row at all is an error naming the file; FIXED sites on no catalog copy's span (here: the
/// right positions under another chrom name) warn in a non-strict run and abort a strict one.
#[test]
fn popfix_table_sanity_empty_table_errors_and_off_span_sites_warn_or_abort() {
    let d = scratch("sanity");
    let b = baseline(&d);

    let empty = write(&d, "empty.tsv", &pop_table(&b.sites, "POLYMORPHIC", &[]));
    let (o, _) = run(&d, "empty", &["--psv-population", &empty]);
    assert!(!o.status.success(), "a table with no FIXED site must be refused");
    assert!(stderr(&o).contains(&format!("no FIXED sites in {empty}")), "{}", stderr(&o));

    let renamed: BTreeSet<(String, u64)> = b.sites.iter().map(|(c, p)| (format!("chr_{c}"), *p)).collect();
    let off = write(&d, "off_span.tsv", &pop_table(&renamed, "FIXED", &[]));
    let (o, ns) = run(&d, "off_ns", &["--psv-population", &off]);
    assert!(o.status.success(), "a non-strict run with off-span sites still runs:\n{}", stderr(&o));
    let err = stderr(&o);
    assert!(err.contains("WARNING") && err.contains("0 of") && err.contains("chrom naming"), "{err}");
    assert_eq!(read(&ns, "assignments.tsv"), b.assignments);
    let (o, _) = run(&d, "off_strict", &["--psv-population", &off, "--psv-fixed-only"]);
    assert!(!o.status.success(), "a strict run with no FIXED site on any catalog copy must abort");
    assert!(stderr(&o).contains("catalog copy's span"), "{}", stderr(&o));

    // the count is printed on a normal run
    let ok = write(&d, "ok.tsv", &pop_table(&b.sites, "FIXED", &[]));
    let (o, _) = run(&d, "ok", &["--psv-population", &ok]);
    assert!(o.status.success(), "{}", stderr(&o));
    let n = b.sites.len();
    assert!(stderr(&o).contains(&format!("{n} FIXED sites, {n} inside a catalog copy's span")), "{}", stderr(&o));
}

/// 5. Guards through the binary: they run before any read is touched and name the flag.
#[test]
fn popfix_guards_need_families_and_refuse_absent_copies() {
    let d = scratch("guards");
    let table = write(&d, "t.tsv", &format!("{POP_HEADER}c1\t520\tFAM1\t0\tFIXED\ttest\t0\t0\t0\t0\n"));

    let (o, _) = run_raw(&d, "nofam", &["--psv-population", &table]);
    assert!(!o.status.success(), "--psv-population without --families must be refused");
    let err = stderr(&o);
    assert!(err.contains("--psv-population") && err.contains("--families"), "{err}");

    let (o, _) = run(&d, "absent", &["--psv-population", &table, "--absent-copies"]);
    assert!(!o.status.success(), "--psv-population with --absent-copies must be refused");
    let err = stderr(&o);
    // the population guard itself (the --families roster guard would name --absent-copies but not --psv-population)
    assert!(err.contains("--psv-population") && err.contains("--absent-copies"), "{err}");
}

/// R31: under `--no-as-tied-only` the placement fallback assigns molecules (sole candidates Tied by R24 among them) that
/// no FIXED column backs, so `--psv-fixed-only` there is refused unless `--no-placement-assign` turns placement off.
/// Non-strict `--psv-population` under `--no-as-tied-only` is unchanged.
#[test]
fn popfix_fixed_only_refuses_no_as_tied_only_placement() {
    let d = scratch("placement");
    let table = write(&d, "t.tsv", &format!("{POP_HEADER}c1\t520\tFAM1\t0\tFIXED\ttest\t0\t0\t0\t0\n"));

    let (o, _) = run(&d, "strict", &["--psv-population", &table, "--psv-fixed-only", "--no-as-tied-only"]);
    assert!(!o.status.success(), "--psv-fixed-only with --no-as-tied-only (placement on) must be refused");
    let err = stderr(&o);
    assert!(err.contains("--no-as-tied-only") && err.contains("--no-placement-assign"), "{err}");
    let (o, _) = run(
        &d,
        "strict_first",
        &["--psv-population", &table, "--psv-fixed-only", "--no-as-tied-only", "--placement-first"],
    );
    assert!(!o.status.success(), "--placement-first is still a placement path");
    assert!(stderr(&o).contains("--no-as-tied-only"), "{}", stderr(&o));

    // controls: placement off, and the non-strict table, both run
    let (o, _) = run(
        &d,
        "strict_noplace",
        &["--psv-population", &table, "--psv-fixed-only", "--no-as-tied-only", "--no-placement-assign"],
    );
    assert!(o.status.success(), "strict with --no-placement-assign must run:\n{}", stderr(&o));
    let (o, _) = run(&d, "nonstrict", &["--psv-population", &table, "--no-as-tied-only"]);
    assert!(o.status.success(), "non-strict --psv-population under --no-as-tied-only is unchanged:\n{}", stderr(&o));
}
