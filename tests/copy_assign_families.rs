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

// ---- 4. RUSTLE_JUNCTION_FUZZ_BP regression --------------------------------------------------------

/// WHAT THIS PROVES: `RUSTLE_JUNCTION_FUZZ_BP` unset and explicitly `0` are the same "off" state -- both
/// must produce byte-identical `--gtf` output, since both resolve to `tolerance_bp == 0`
/// (`merge_fuzzy_skeletons`'s own proven no-op) before the merge function is ever called. That is a real,
/// valid off-state parity contract and this test genuinely exercises it.
///
/// WHAT THIS DOES NOT PROVE, and a correction of an earlier claim here (confirmed 2026-09-16, Task 5
/// fix round 1): this comment used to say `--gtf` on the `same_chrom_supplement` fixture emits a
/// 0-BYTE `.gtf` file. That was true of the pre-2026-09-16 dedup behaviour (still reproducible today
/// via `RUSTLE_LEGACY_PLACEMENT_DEDUP=1`, and on the committed pre-fix `copy_assign` binary) -- both
/// measured 0 `.gtf` lines on this exact invocation. This branch's cross-window dedup fix changed
/// that: the CURRENT binary emits 9 `.gtf` rows on this fixture. Re-measuring the second claim against
/// that current binary: `RUSTLE_JUNCTION_FUZZ_BP=200` vs unset on this fixture is NO LONGER
/// byte-identical (9 rows unset vs 6 rows at `FUZZ_BP=200` -- the merge now has real chains to act on).
/// The test below only asserts unset == explicit-`0`, which stays true and is a real no-op check
/// either way (`tolerance_bp == 0` short-circuits before `merge_fuzzy_skeletons` is ever called,
/// independent of how many rows exist), so this test's own assertions are unaffected by the dedup fix
/// -- but it still cannot and does not exercise the merge's ON-state (tolerance > 0) at all. The
/// positive-control evidence that the merge fires correctly (and, at its pre-registered tolerance, harmfully)
/// on real `--gtf` output lives elsewhere: `bench/CHR20_ASSEMBLER_COMPARISON.md`'s real chr20 acceptance run
/// (976 -> 836 transcripts, matching intron chains 345 -> 284, a real measured effect), and
/// `merge_fuzzy_skeletons`'s own unit tests (`denovo_assemble.rs`, `fuzzy_merge_*`), which exercise the merge
/// logic directly on synthetic skeletons.
#[test]
fn junction_fuzz_off_state_parity_unset_vs_explicit_zero() {
    let d_unset = scratch("junction_fuzz_unset");
    let (o_unset, out_unset) = run(&d_unset, &["--no-refine", "--gtf"]);
    assert!(o_unset.status.success(), "unset run failed:\n{}", stderr(&o_unset));

    let d_zero = scratch("junction_fuzz_zero");
    let out_s_zero = d_zero.join("o").to_str().expect("utf-8 path").to_string();
    let o_zero = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
        .args(["--region", "c1:200-600", "--out", &out_s_zero])
        .args(["--no-refine", "--gtf"])
        .env("RUSTLE_JUNCTION_FUZZ_BP", "0")
        .output()
        .expect("copy_assign failed to spawn");
    assert!(o_zero.status.success(), "explicit-zero run failed:\n{}", String::from_utf8_lossy(&o_zero.stderr));

    for ext in ["assignments.tsv", "families.tsv", "quant.tsv", "gtf"] {
        let a = std::fs::read(format!("{out_unset}.{ext}")).unwrap_or_else(|e| panic!("read {out_unset}.{ext}: {e}"));
        let b = std::fs::read(format!("{out_s_zero}.{ext}")).unwrap_or_else(|e| panic!("read {out_s_zero}.{ext}: {e}"));
        assert_eq!(a, b, "RUSTLE_JUNCTION_FUZZ_BP unset vs explicit 0 must be byte-identical for .{ext}");
    }
}

/// `--gtf-refine` is validated up front: it needs `--gtf`, and an unknown component is an error. Its ON-state
/// behaviour is not exercised here: with the current (2026-09-16 dedup-fixed) binary this fixture's `--gtf`
/// emits 9 rows (measured via `cmp`), but none of them trigger any `--gtf-refine` rule -- `--gtf-refine all`
/// produces byte-identical output (every `.tsv`/`.gtf`) to unset on this exact invocation. So the third
/// assertion below (`--gtf-refine all` succeeds) proves only accept-and-run, not that any rule actually
/// fired. ON-state evidence lives elsewhere: the `gtf_refine` unit tests (`src/rustle/vg_family/gtf_refine.rs`)
/// and the chr20 fidelity anchors (docs/superpowers/plans/2026-09-16-gtf-refine-and-dedup-fix.md, Task 6).
#[test]
fn gtf_refine_requires_gtf_and_known_components() {
    let d = scratch("gtf_refine_cli");
    let (o, _) = run(&d, &["--no-refine", "--gtf-refine", "all"]);
    assert!(!o.status.success(), "--gtf-refine without --gtf must fail");
    assert!(stderr(&o).contains("--gtf-refine is only meaningful with --gtf"), "{}", stderr(&o));
    let (o, _) = run(&d, &["--no-refine", "--gtf", "--gtf-refine", "bogus"]);
    assert!(!o.status.success(), "an unknown component must fail");
    let (o, _) = run(&d, &["--no-refine", "--gtf", "--gtf-refine", "all"]);
    assert!(o.status.success(), "{}", stderr(&o));
}
