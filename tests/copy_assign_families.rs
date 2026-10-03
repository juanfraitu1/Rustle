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
/// and `params.tsv`, see `src/bin/copy_assign.rs:4286` and `:4886`) must come out byte-for-byte identical
/// to a run with the flag added. This is the single most important untested claim from the copy-discovery
/// feature itself (Tasks 1-4), so it is checked directly against the real binary, not the library code.
///
/// ⚠ Runs under `--families` (final whole-branch review, Finding 5). The original version of this test
/// used the PLAIN invocation, which on this fixture yields header-only `assignments.tsv`/`families.tsv`/
/// `quant.tsv` in BOTH arms: the byte-parity assertion was real but VACUOUS -- it compared empty files and
/// never entered the per-family discovery path at all. With `--families GWFAM1` + `--copies-fa`, both arms
/// produce non-empty assignment/quant output AND `--discover-copies` actually runs its per-family
/// clustering, so the parity being asserted is parity of a real result.
#[test]
fn discover_copies_off_by_default_is_byte_identical() {
    // Two separate scratch dirs give the two runs distinct --out prefixes (the `run()` helper always
    // writes to `<dir>/o`), but the CATALOG must be one shared file: `params.tsv` records the `--families`
    // path verbatim, so two per-dir copies of the same table would differ there for a reason that has
    // nothing to do with this flag.
    let d_cat = scratch("discover_cat");
    let fam = write(&d_cat, "cat.copies.tsv", GWFAM1_TSV);

    let d_off = scratch("discover_off");
    let (o_off, out_off) =
        run(&d_off, &["--families", &fam, "--copies-fa", &format!("{FIX}/out_default.copies.fa")]);
    assert!(o_off.status.success(), "flag-off run failed:\n{}", stderr(&o_off));

    let d_on = scratch("discover_on");
    let (o_on, out_on) = run(
        &d_on,
        &["--families", &fam, "--copies-fa", &format!("{FIX}/out_default.copies.fa"), "--discover-copies"],
    );
    assert!(o_on.status.success(), "flag-on run failed:\n{}", stderr(&o_on));

    // The parity comparison is only meaningful if these files actually carry rows.
    assert!(!col(&read(&out_off, "quant.tsv"), 2).is_empty(), "the flag-off arm must produce real quant rows");
    assert!(
        !col(&read(&out_off, "assignments.tsv"), 0).is_empty(),
        "the flag-off arm must produce real assignment rows"
    );

    for ext in ["assignments.tsv", "families.tsv", "quant.tsv", "famcn_readonly.tsv", "params.tsv", "family_join.tsv"] {
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
    let report = read(&out_on, "discovered_copies.tsv");
    assert_eq!(
        report.lines().next(),
        Some("family_id\tchrom\tstart\tend\tstrand\tn_supporting_reads\tread_names\tnearest_copy_tid\tnearest_copy_distance"),
        "the report header must carry the strand column (final whole-branch review, Finding 4): {report}"
    );
}

/// ⚠ THE CROSS-FAMILY POOLING BUG, at the real-binary level (final whole-branch review, Finding 1).
///
/// Pre-fix, `cluster_tie_partners` was handed the WHOLE region's AS-tied read list once per family, so one
/// out-of-catalog site was emitted once per family in the region, with the identical `read_names` list
/// under a different `family_id` -- confirmed on real data, where one site came out under 3-8 ids.
///
/// The catalog here is built to make that visible with the committed fixture BAM. It holds THREE families
/// over one region set:
///   * `GWFAM1` (c1:250-460 + c1:380-550) -- owns `read_same_0/1/2`, whose two AS-tied placements both sit
///     inside its own copies, so it legitimately discovers nothing;
///   * `GWFAM0` (c1:0-260 + c2:0-260) -- also considers `read_same_*` (their primary overlaps c1:0-260),
///     and their OTHER max-AS placement at c1:380-550 is outside every GWFAM0 copy: the one legitimate
///     candidate row, present both before and after the fix;
///   * `FAMX` (c1:0-59 + c1:100-150) -- a DECOY at the far end of the same region. No `read_same_*`
///     placement comes near it, so `FAMX` never considers those reads (it has no `assignments.tsv` row for
///     them at all) and must report nothing. Pre-fix it reported `c1:250-550` naming
///     `read_same_0,read_same_1,read_same_2` -- reads that were never its to reason about. VERIFIED by
///     re-introducing the bug against this exact catalog: the buggy binary emits that extra `FAMX` row and
///     this test fails; the fixed binary emits only the `GWFAM0` row.
#[test]
fn discovered_copies_are_never_pooled_across_families() {
    let d = scratch("discover_xfam");
    // FAMX (decoy) + GWFAM1 (same-chrom) + GWFAM0 (cross-chrom), swept over both contigs. No
    // `--copies-fa`: FAMX is synthetic, so the sequences are rebuilt from the genome at the catalog's own
    // exon coordinates (the documented fallback, pinned by
    // `families_without_copies_fa_rebuilds_the_sequences_from_the_genome` above).
    const FAMX_TSV: &str = "FAMX\t0\tX0\tc1\t0\t59\t1\t+\t3\t0-59\n\
FAMX\t1\tX1\tc1\t100\t150\t1\t+\t3\t100-150\n";
    let tsv = std::fs::read_to_string(format!("{FIX}/out_default.copies.tsv")).unwrap();
    let fam0: String = tsv.lines().filter(|l| l.starts_with("GWFAM0\t")).collect::<Vec<_>>().join("\n");
    // GWFAM1_TSV already carries the header row; FAMX/GWFAM0 rows append to it.
    let (hdr, gwfam1) = GWFAM1_TSV.split_at(GWFAM1_TSV.find('\n').unwrap() + 1);
    let fam = write(&d, "three.copies.tsv", &format!("{hdr}{FAMX_TSV}{gwfam1}{fam0}\n"));
    let regions = write(&d, "regions.txt", "c1:0-600\nc2:0-320\n");
    let out = d.join("o");
    let out_s = out.to_str().expect("utf-8 path").to_string();
    let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
        .args(["--regions", &regions, "--out", &out_s])
        .args(["--families", &fam])
        .arg("--discover-copies")
        .output()
        .expect("copy_assign failed to spawn");
    assert!(o.status.success(), "run failed:\n{}", stderr(&o));

    // All three families must really be in play -- otherwise "no cross-family attribution" is vacuous.
    let mut fams: Vec<String> = col(&read(&out_s, "families.tsv"), 0);
    fams.sort();
    fams.dedup();
    assert_eq!(
        fams,
        vec!["FAMX".to_string(), "GWFAM0".to_string(), "GWFAM1".to_string()],
        "all three catalog families must be assigned"
    );

    let report = read(&out_s, "discovered_copies.tsv");
    // The decoy considered none of the region's AS-tied reads, so it must claim nothing.
    assert!(
        !col(&report, 0).iter().any(|f| f == "FAMX"),
        "FAMX considered no AS-tied read yet claims a candidate -- cross-family pooling:\n{report}"
    );
    // Non-vacuous: this fixture really does yield a candidate (GWFAM0's own tied reads have a max-AS
    // placement at c1:380-550, outside every GWFAM0 copy). Without a row, everything below is empty-set
    // true and the test would prove nothing.
    assert!(report.lines().skip(1).any(|l| !l.trim().is_empty()), "expected at least one candidate:\n{report}");

    // (1) THE DIRECT INVARIANT: every read named by a discovered row must be a read the REPORTING family
    // actually considered. `assignments.tsv` is that ground truth, emitted from the same `fa.assignments`
    // `discover_copies_for_family` now restricts on. Pre-fix, a family was handed the whole region's tied
    // reads, so a row could name reads that appear nowhere under its own family_id.
    let assignments = read(&out_s, "assignments.tsv");
    let mut considered: std::collections::HashMap<String, std::collections::HashSet<String>> =
        std::collections::HashMap::new();
    for l in assignments.lines().skip(1).filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = l.split('\t').collect();
        considered.entry(f[1].to_string()).or_default().insert(f[0].to_string());
    }

    // (2) and no two families may claim the SAME site with the SAME supporting reads.
    let mut by_site: std::collections::HashMap<(String, String, String, String), Vec<String>> =
        std::collections::HashMap::new();
    for l in report.lines().skip(1).filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = l.split('\t').collect();
        assert_eq!(f.len(), 9, "row must have 9 columns (strand included): {l}");
        let (fid, names) = (f[0].to_string(), f[6]);
        let mine = considered.get(&fid).cloned().unwrap_or_default();
        for n in names.split(',') {
            assert!(
                mine.contains(n),
                "{fid} reports read {n}, which it never considered (not in its assignments.tsv rows) \
                 -- cross-family pooling:\n{report}\n{assignments}"
            );
        }
        by_site
            .entry((f[1].to_string(), f[2].to_string(), f[3].to_string(), names.to_string()))
            .or_default()
            .push(fid);
    }
    for (site, ids) in &by_site {
        assert_eq!(ids.len(), 1, "site {site:?} reported under {ids:?} -- cross-family pooling:\n{report}");
    }
    // Strand must be a real call, never blank or a placeholder character.
    for s in col(&report, 4) {
        assert!(s == "+" || s == "-", "strand column must be + or -, got {s:?}:\n{report}");
    }
    // And the distance column must never be a raw `u64::MAX` sentinel (final whole-branch review, Minor 6).
    for d in col(&report, 8) {
        assert_ne!(d, "18446744073709551615", "absent nearest copy must print NA:\n{report}");
    }
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
    // Not asserted here: whether a molecule AS-tied ACROSS the two chromosomes is scored as `n_candidates == 2`.
    // The certificate now takes each read's chromosome (`read_chroms`) and `detect_and_assign` hands the family
    // its reads on BOTH c1 and c2, so a read overlaps only the copies on its own chromosome even where `c1:0-260`
    // and `c2:0-260` coincide numerically (pinned in `denovo_pipeline.rs` by
    // `a_cross_chromosome_family_is_assigned_the_reads_on_every_chromosome_it_has_a_copy_on`).
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

/// `--only-families` / `--skip-families` (spec 2026-10-02 §7, the driver's two-run O2 split over the families with an O3
/// candidate copy): each run assigns exactly its selection of the supplied table and the two runs cover the full run's
/// families; an unknown id or a selection that keeps nothing aborts, and so does a list without `--families`.
#[test]
fn only_and_skip_families_split_the_supplied_roster() {
    let d = scratch("split");
    let regions = write(&d, "regions.txt", "c1:0-600\nc2:0-320\n");
    let go = |name: &str, extra: &[&str]| -> (Output, String) {
        let out = d.join(name).to_str().expect("utf-8 path").to_string();
        let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
            .args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
            .args(["--regions", &regions, "--out", &out])
            .args(["--families", &format!("{FIX}/out_default.copies.tsv"), "--copies-fa", &format!("{FIX}/out_default.copies.fa")])
            .args(extra)
            .output()
            .expect("copy_assign failed to spawn");
        (o, out)
    };
    let list = write(&d, "gwfam0.txt", "GWFAM0\n");
    let (o, all) = go("all", &[]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(!stderr(&o).contains("families selected"), "no list: nothing selected, nothing said");
    let (o, only) = go("only", &["--only-families", &list]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(stderr(&o).contains("1 of 2 families selected (--only-families"), "{}", stderr(&o));
    let (o, skip) = go("skip", &["--skip-families", &list]);
    assert!(o.status.success(), "{}", stderr(&o));
    let mut fams_all = col(&read(&all, "families.tsv"), 0);
    fams_all.sort();
    assert_eq!(fams_all, vec!["GWFAM0", "GWFAM1"]);
    assert_eq!(col(&read(&only, "families.tsv"), 0), vec!["GWFAM0"]);
    assert_eq!(col(&read(&skip, "families.tsv"), 0), vec!["GWFAM1"]);
    for (out, fam) in [(&only, "GWFAM0"), (&skip, "GWFAM1")] {
        for t in ["assignments.tsv", "quant.tsv", "family_join.tsv", "famcn_readonly.tsv"] {
            let ids = col(&read(out, t), 0 + (t == "assignments.tsv") as usize);
            assert!(ids.iter().all(|f| f == fam), "{out}.{t} holds another family: {ids:?}");
        }
    }
    // the run certificate names the list it was given, and only then
    assert!(read(&only, "params.tsv").contains(&format!("only_families\t{list}\n")));
    assert!(read(&skip, "params.tsv").contains(&format!("skip_families\t{list}\n")));
    assert!(!read(&all, "params.tsv").contains("_families\t"), "no list, no row");
    let typo = write(&d, "typo.txt", "GWFAM0\nGWFAM9\n");
    let (o, _) = go("typo", &["--only-families", &typo]);
    assert!(!o.status.success() && stderr(&o).contains("GWFAM9"), "{}", stderr(&o));
    let both = write(&d, "both.txt", "GWFAM0\nGWFAM1\n");
    let (o, _) = go("none", &["--skip-families", &both]);
    assert!(!o.status.success() && stderr(&o).contains("leave no family"), "{}", stderr(&o));
    let (o, _) = run(&d, &["--only-families", &list]);
    assert!(!o.status.success() && stderr(&o).contains("only meaningful with --families"), "{}", stderr(&o));
}

/// Task 9 review (spec 2026-10-02 §7): a family made cross-chromosome by a zero-read copy on another contig, which is
/// what the driver's augmentation does to every family with an O3 candidate, is bound to its `~xchrom~` key, so the real
/// regions holding it are swept with NO family bound. The §6gz tie-outside test there had no unit to be inside of: it
/// registered every AS-tied molecule as tied outside (a process-wide registry), and the family's own verdicts on those
/// molecules then read `tie_outside_catalog 1` and could never be `assigned`. With the extra copy, the rows of the
/// single-chromosome run keep their status and `tie_outside_catalog`. The molecules the c2 copy pools (`read_cross_*`,
/// AS-tied between c1:0-260 and c2:0-260) still read 1: their c2 placement lies outside every unit of the family, so the
/// family's own sweep still flags a real outside competitor.
#[test]
fn a_zero_read_copy_on_another_contig_does_not_mark_the_family_tied_outside() {
    let d = scratch("xchrom_zero_read_copy");
    let regions = write(&d, "regions.txt", "c1:0-600\nc2:0-320\n");
    let only = write(&d, "only.txt", "GWFAM1\n");
    let tsv = std::fs::read_to_string(format!("{FIX}/out_default.copies.tsv")).unwrap();
    let fa = std::fs::read_to_string(format!("{FIX}/out_default.copies.fa")).unwrap();
    let genome = std::fs::read_to_string(format!("{FIX}/genome.fa")).unwrap();
    let c2: String = genome.split('>').find(|r| r.starts_with("c2")).expect("c2 in genome.fa").lines().skip(1).collect();
    let tsv_b = write(&d, "b.copies.tsv", &format!("{tsv}GWFAM1\t2\tCAND_c2\tc2\t270\t320\t1\t+\t0\t270-320\t0\n"));
    let fa_b = write(&d, "b.copies.fa", &format!("{fa}>GWFAM1|2|c2:270-320|+|nexon=1\n{}\n", &c2[270..320]));
    let go = |name: &str, tsv: &str, fa: &str| -> (Output, String) {
        let out = d.join(name).to_str().expect("utf-8 path").to_string();
        let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
            .args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
            .args(["--regions", &regions, "--families", tsv, "--copies-fa", fa, "--only-families", &only, "--out", &out])
            .output()
            .expect("copy_assign failed to spawn");
        (o, out)
    };
    let (o, a) = go("a", &format!("{FIX}/out_default.copies.tsv"), &format!("{FIX}/out_default.copies.fa"));
    assert!(o.status.success(), "{}", stderr(&o));
    let (o, b) = go("b", &tsv_b, &fa_b);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(stderr(&o).contains("GWFAM1 spans 2 chromosomes"), "the extra copy must make GWFAM1 cross-chromosome:\n{}", stderr(&o));
    // read -> (status, tie_outside_catalog), columns by name
    let rows = |out: &str| -> std::collections::BTreeMap<String, (String, String)> {
        let t = read(out, "assignments.tsv");
        let head: Vec<&str> = t.lines().next().unwrap().split('\t').collect();
        let at = |name: &str| head.iter().position(|c| *c == name).unwrap_or_else(|| panic!("no {name} column"));
        let (i_read, i_status, i_out) = (at("read_name"), at("status"), at("tie_outside_catalog"));
        t.lines()
            .skip(1)
            .map(|l| l.split('\t').collect::<Vec<_>>())
            .map(|f| (f[i_read].to_string(), (f[i_status].to_string(), f[i_out].to_string())))
            .collect()
    };
    let (ra, rb) = (rows(&a), rows(&b));
    assert!(!ra.is_empty(), "the single-chromosome run must have rows to compare");
    for (read, v) in &ra {
        assert_eq!(v.1, "0", "{read}: tied only between GWFAM1's own copies in the single-chromosome run");
        assert_eq!(rb.get(read), Some(v), "{read}: status / tie_outside_catalog changed by the zero-read c2 copy");
    }
    let cross: Vec<(&String, &(String, String))> = rb.iter().filter(|(r, _)| r.starts_with("read_cross_")).collect();
    assert_eq!(cross.len(), 3, "the c2 copy pools the 3 read_cross molecules: {rb:?}");
    assert!(cross.iter().all(|(_, v)| v.1 == "1"), "their c2 placement is outside every GWFAM1 unit: {cross:?}");
}
