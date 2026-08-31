//! `RUSTLE_XFAM_RECONCILE` — CROSS-FAMILY reconciliation of a molecule assigned in TWO families.
//!
//! ## Why there was no such test before
//!
//! Every pre-existing `copy_assign` e2e fixture supplies exactly ONE family (`copy_assign_families.rs`'s
//! `GWFAM1`, and `tests/fixtures/same_chrom_supplement`'s catalog has one same-chrom family), and every
//! in-tree assignment-count assertion in `denovo_pipeline.rs` is single-family. A defect whose necessary
//! condition is TWO families' significance gates accepting the same molecule was therefore invisible to
//! the whole suite — which is why it shipped. `tests/fixtures/xfam_conflict` is the first two-family
//! e2e fixture; see its `make_bam.py` for the layout and the three planted molecules.
//!
//! ## The three strata, one planted molecule each
//!
//! * `MOL_CONTRA` — TWO records, one primary at `XFA` copy0 (`x1:1000-1600`) and one SECONDARY at
//!   `XFB` copy1 (`x1:10000-10600`). Disjoint copies, DIFFERENT records ⇒ `cross_family_contradiction`.
//! * `MOL_READTHROUGH` — ONE record, `x1:3000` `600M4400N600M`, whose N-gap spans `XFA` copy1 to
//!   `XFB` copy0. Disjoint copies, SAME record ⇒ `readthrough_span`.
//! * `MOL_SHARED` — ONE record inside `x1:15000-15600`, the interval `XFA` copy2 and `XFB` copy2 BOTH
//!   claim ⇒ `shared_locus` (an O1 partition artifact, shared with the 8 support reads at that locus).
//!
//! Only the first abstains. The other two are reported and left assigned, on purpose.

use std::path::PathBuf;
use std::process::{Command, Output};

const FIX: &str = "tests/fixtures/xfam_conflict";

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("copy_assign_xfam").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("create scratch dir");
    d
}

/// Run the fixture with `RUSTLE_XFAM_RECONCILE` set to `mode` (`None` = leave it UNSET).
/// `--phase`/`--posterior`/`--dump-psv` are on so every status-carrying output is produced.
fn run_mode(dir: &PathBuf, mode: Option<&str>) -> (Output, String) {
    let out = dir.join("o");
    let out_s = out.to_str().expect("utf-8 path").to_string();
    let mut c = Command::new(env!("CARGO_BIN_EXE_copy_assign"));
    c.args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
        .args(["--families", &format!("{FIX}/copies.tsv"), "--copies-fa", &format!("{FIX}/copies.fa")])
        .args(["--region", "x1:0-20000", "--dump-psv", "--phase", "--posterior"])
        .args(["--out", &out_s]);
    match mode {
        Some(m) => {
            c.env("RUSTLE_XFAM_RECONCILE", m);
        }
        None => {
            c.env_remove("RUSTLE_XFAM_RECONCILE");
        }
    }
    let o = c.output().expect("copy_assign failed to spawn");
    assert!(o.status.success(), "run ({mode:?}) failed:\n{}", String::from_utf8_lossy(&o.stderr));
    (o, out_s)
}

fn read(out: &str, ext: &str) -> String {
    std::fs::read_to_string(format!("{out}.{ext}")).unwrap_or_else(|e| panic!("read {out}.{ext}: {e}"))
}

fn md5_of(out: &str, ext: &str) -> String {
    // No md5 crate in-tree; a content hash is all the byte-identity assertions need.
    let b = std::fs::read(format!("{out}.{ext}")).unwrap_or_else(|e| panic!("read {out}.{ext}: {e}"));
    format!("{:016x}:{}", fxhash(&b), b.len())
}

fn fxhash(b: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &x in b {
        h ^= x as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// `(read_name, family_id) -> status` from a TSV whose status lives at 0-based column `col`.
fn statuses(text: &str, col: usize) -> std::collections::BTreeMap<(String, String), String> {
    text.lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            ((f[0].to_string(), f[1].to_string()), f[col].to_string())
        })
        .collect()
}

/// Rows of `<out>.xfam_conflicts.tsv` as `(read_name, stratum, demoted)`.
fn conflicts(text: &str) -> Vec<(String, String, String)> {
    text.lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[1].to_string(), f[14].to_string())
        })
        .collect()
}

// ---- the rule ---------------------------------------------------------------------------------------

/// THE CLAUSE: two independent placements naming disjoint loci ⇒ the molecule abstains in EVERY family.
/// A molecule that SPANS both loci (one record) and one inside an interval both families claim do NOT.
#[test]
fn abstain_demotes_only_the_two_record_contradiction() {
    let d = scratch("abstain");
    let (_, out) = run_mode(&d, Some("abstain"));
    let st = statuses(&read(&out, "assignments.tsv"), 3);
    let g = |r: &str, f: &str| st.get(&(r.to_string(), f.to_string())).cloned().unwrap_or_default();

    // CLAUSE.
    for fam in ["XFA", "XFB"] {
        assert_eq!(
            g("MOL_CONTRA", fam),
            "ambiguous",
            "two independent placements naming disjoint loci => the molecule abstains in every family \
             (MOL_CONTRA in {fam})"
        );
        assert_eq!(
            g("MOL_READTHROUGH", fam),
            "assigned",
            "a molecule whose single record SPANS both loci really does come from both, so \
             'a molecule cannot come from two disjoint loci' does not apply to it (MOL_READTHROUGH in {fam})"
        );
        assert_eq!(
            g("MOL_SHARED", fam),
            "assigned",
            "two families claiming ONE overlapping interval is an O1 partition artifact, not an O2 \
             assignment error; demoting it would charge that defect to O2's abstention rate \
             (MOL_SHARED in {fam})"
        );
    }
}

/// THE VALUE, pinned separately from the clause: the exact per-status row counts on this fixture.
#[test]
fn abstain_row_counts_are_exactly_two_demotions() {
    let d = scratch("abstain_counts");
    let (_, out) = run_mode(&d, Some("abstain"));
    let a = read(&out, "assignments.tsv");
    let mut n = std::collections::BTreeMap::new();
    for v in statuses(&a, 3).values() {
        *n.entry(v.clone()).or_insert(0usize) += 1;
    }
    assert_eq!(n.get("assigned").copied().unwrap_or(0), 52, "assigned rows: {n:?}");
    assert_eq!(n.get("ambiguous").copied().unwrap_or(0), 2, "ambiguous rows: {n:?}");
    assert_eq!(a.lines().count() - 1, 54, "total rows must be INVARIANT across modes (no row is deleted)");
}

// ---- the OFF arm ------------------------------------------------------------------------------------

/// Unset, `off` and the DEFAULT must be one and the same byte stream, and `off` must write no side file.
#[test]
fn off_arm_is_byte_identical() {
    let du = scratch("off_unset");
    let (_, u) = run_mode(&du, None);
    let df = scratch("off_explicit");
    let (_, f) = run_mode(&df, Some("off"));
    for ext in [
        "assignments.tsv", "families.tsv", "quant.tsv", "psv_reads.tsv", "psv_copies.tsv", "psv_cols.tsv",
        "posterior.tsv", "phased_reads.tsv", "phase_blocks.tsv", "phased_haplotypes.tsv", "phase.gfa",
        "exon.gfa", "famcn_readonly.tsv", "family_join.tsv",
    ] {
        assert_eq!(md5_of(&u, ext), md5_of(&f, ext), "unset and off must be byte-identical in {ext}");
    }
    for out in [&u, &f] {
        assert!(
            !std::path::Path::new(&format!("{out}.xfam_conflicts.tsv")).exists(),
            "off must write NO side file: {out}.xfam_conflicts.tsv exists"
        );
    }
}

/// `report` changes no pre-existing byte, and its side file lists exactly the three strata.
#[test]
fn report_arm_changes_no_existing_output() {
    let doff = scratch("rep_off");
    let (_, off) = run_mode(&doff, Some("off"));
    let drep = scratch("rep_on");
    let (o, rep) = run_mode(&drep, Some("report"));
    for ext in [
        "assignments.tsv", "families.tsv", "quant.tsv", "psv_reads.tsv", "psv_copies.tsv", "psv_cols.tsv",
        "posterior.tsv", "phased_reads.tsv", "phase_blocks.tsv", "phased_haplotypes.tsv", "phase.gfa",
        "exon.gfa", "famcn_readonly.tsv", "family_join.tsv",
    ] {
        assert_eq!(md5_of(&off, ext), md5_of(&rep, ext), "report must not change {ext}");
    }
    let c = conflicts(&read(&rep, "xfam_conflicts.tsv"));

    // CLAUSE: the three strata are present, and only the contradiction is marked demoted.
    let by: std::collections::BTreeMap<&str, &str> =
        c.iter().map(|(r, s, _)| (r.as_str(), s.as_str())).collect();
    assert_eq!(by.get("MOL_CONTRA").copied(), Some("cross_family_contradiction"));
    assert_eq!(by.get("MOL_READTHROUGH").copied(), Some("readthrough_span"));
    assert_eq!(by.get("MOL_SHARED").copied(), Some("shared_locus"));
    assert!(
        c.iter().all(|(_, s, d)| (s == "cross_family_contradiction") == (d == "true")),
        "under `report` the `demoted` column still marks the stratum the abstain arm WOULD demote, \
         and only that stratum: {c:?}"
    );
    // and the stderr summary must carry its denominators (a rate without one is not reportable here).
    let err = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(err.contains("[xfam] RUSTLE_XFAM_RECONCILE=report"), "{err}");
    assert!(err.contains("distinct molecules"), "the summary must quote its denominators:\n{err}");

    // VALUE, pinned separately: 1 contradiction + 1 readthrough + 9 shared (MOL_SHARED plus the 8
    // support reads at the interval both families claim).
    let mut n = std::collections::BTreeMap::new();
    for (_, s, _) in &c {
        *n.entry(s.clone()).or_insert(0usize) += 1;
    }
    assert_eq!(n.get("cross_family_contradiction").copied().unwrap_or(0), 1, "{n:?}");
    assert_eq!(n.get("readthrough_span").copied().unwrap_or(0), 1, "{n:?}");
    assert_eq!(n.get("shared_locus").copied().unwrap_or(0), 9, "{n:?}");
    assert_eq!(c.len(), 11, "total contested pairs: {n:?}");
}

// ---- the regression that catches a post-pass ---------------------------------------------------------

/// A demotion applied as a post-pass over `assign_rows` would fix `.assignments.tsv` and leave
/// `.psv_reads.tsv`, `.posterior.tsv` and `.phased_reads.tsv` disagreeing with it. All four are emitted
/// from `fa.assignments` inside the drain, so they must agree row-for-row on `(read_name, family_id)`.
#[test]
fn status_consistency_across_outputs() {
    let d = scratch("consistency");
    let (_, out) = run_mode(&d, Some("abstain"));
    let a = statuses(&read(&out, "assignments.tsv"), 3);
    for (ext, col) in [("psv_reads.tsv", 3), ("posterior.tsv", 2), ("phased_reads.tsv", 5)] {
        let b = statuses(&read(&out, ext), col);
        assert_eq!(a.len(), b.len(), "{ext} has a different row set than assignments.tsv");
        for (k, v) in &a {
            assert_eq!(
                b.get(k),
                Some(v),
                "every status emit site must see the EFFECTIVE status: {ext} disagrees at {k:?}"
            );
        }
    }
    // and the --phase haplotype must follow the effective status, not the raw one.
    let ph = read(&out, "phased_reads.tsv");
    for l in ph.lines().skip(1).filter(|l| l.starts_with("MOL_CONTRA\t")) {
        let f: Vec<&str> = l.split('\t').collect();
        assert_eq!(f[2], "-1", "an abstaining read is unphaseable: {l}");
    }
}

// ---- honest scope ------------------------------------------------------------------------------------

/// `.quant.tsv` is DELIBERATELY unmoved: `n_reads_hard` has no status filter and `abundance`/`ci95` come
/// from the per-family EM, whose observations ignore status. This fails loudly if someone silently
/// "improves" it — the abundance double-count is a separate, two-pass decision.
#[test]
fn quant_is_unmoved_by_demotion() {
    let dr = scratch("quant_rep");
    let (_, rep) = run_mode(&dr, Some("report"));
    let da = scratch("quant_abs");
    let (_, abs) = run_mode(&da, Some("abstain"));
    assert_eq!(
        md5_of(&rep, "quant.tsv"),
        md5_of(&abs, "quant.tsv"),
        "quant.tsv must be byte-identical between report and abstain: a demotion changes STATUS only"
    );
    assert_eq!(md5_of(&rep, "families.tsv"), md5_of(&abs, "families.tsv"), "families.tsv counters too");
}

// ---- the M2 guard ------------------------------------------------------------------------------------

/// The params certificate is what makes an ON arm distinguishable from an OFF arm from its outputs alone.
/// Without it the two arms of an env flag are identical on disk — the recurring M2 defect.
#[test]
fn params_certificate_distinguishes_the_arms() {
    let mut seen = Vec::new();
    for (name, mode) in [("p_unset", None), ("p_off", Some("off")), ("p_rep", Some("report")), ("p_abs", Some("abstain"))] {
        let d = scratch(name);
        let (_, out) = run_mode(&d, mode);
        let p = read(&out, "params.tsv");
        let v = p
            .lines()
            .find_map(|l| l.strip_prefix("xfam_reconcile\t"))
            .unwrap_or_else(|| panic!("params.tsv has no xfam_reconcile row:\n{p}"))
            .to_string();
        assert!(p.starts_with("key\tvalue\n"), "params.tsv header: {p}");
        assert!(p.contains("posterior_prior\t"), "RUSTLE_POSTERIOR_PRIOR must be recorded too:\n{p}");
        seen.push(v);
    }
    assert_eq!(seen, vec!["off", "off", "report", "abstain"], "the EFFECTIVE mode must be recorded");
}

/// An unrecognized value is an ERROR, not a silent `off`: a typo'd flag that quietly disables the pass is
/// how an ON arm gets reported as OFF.
#[test]
fn an_unknown_mode_is_refused() {
    let d = scratch("bad_mode");
    let out = d.join("o");
    let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
        .args(["--families", &format!("{FIX}/copies.tsv"), "--copies-fa", &format!("{FIX}/copies.fa")])
        .args(["--region", "x1:0-20000", "--out", out.to_str().unwrap()])
        .env("RUSTLE_XFAM_RECONCILE", "reprot")
        .output()
        .expect("spawn");
    assert!(!o.status.success(), "a typo'd mode must abort, not run as off");
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("off|report|abstain"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
}

// ---- determinism -------------------------------------------------------------------------------------

/// Same input ⇒ byte-identical output, in every mode and including the new side file; and the
/// region-parallel sweep must equal the serial one (the new pass iterates `works` in the flat region
/// order and uses BTree containers only, so nothing depends on thread scheduling).
#[test]
fn determinism_across_repeats_and_region_threads() {
    for mode in [None, Some("off"), Some("report"), Some("abstain")] {
        let tag = mode.unwrap_or("unset");
        let d1 = scratch(&format!("det_a_{tag}"));
        let (_, a) = run_mode(&d1, mode);
        let d2 = scratch(&format!("det_b_{tag}"));
        let (_, b) = run_mode(&d2, mode);
        for ext in ["assignments.tsv", "families.tsv", "quant.tsv", "psv_reads.tsv", "params.tsv", "phase.gfa"] {
            assert_eq!(md5_of(&a, ext), md5_of(&b, ext), "{tag}: two identical runs differ in {ext}");
        }
        if mode.is_some() && mode != Some("off") {
            assert_eq!(
                md5_of(&a, "xfam_conflicts.tsv"),
                md5_of(&b, "xfam_conflicts.tsv"),
                "{tag}: the side file must be deterministic too"
            );
        }
    }
    // --region-threads > 1 over two regions covering the same fixture.
    let ds = scratch("det_serial");
    let dp = scratch("det_par");
    let mut outs = Vec::new();
    for (d, threads) in [(&ds, "1"), (&dp, "4")] {
        let out = d.join("o").to_str().unwrap().to_string();
        let mut rf = d.clone();
        rf.push("regions.txt");
        std::fs::write(&rf, "x1:0-9000\nx1:9000-20000\n").unwrap();
        let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
            .args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
            .args(["--regions", rf.to_str().unwrap(), "--region-threads", threads])
            .args(["--dump-psv", "--out", &out])
            .env("RUSTLE_XFAM_RECONCILE", "report")
            .output()
            .expect("spawn");
        assert!(o.status.success(), "region-threads={threads} failed:\n{}", String::from_utf8_lossy(&o.stderr));
        outs.push(out);
    }
    for ext in ["assignments.tsv", "families.tsv", "quant.tsv", "xfam_conflicts.tsv"] {
        assert_eq!(
            md5_of(&outs[0], ext),
            md5_of(&outs[1], ext),
            "--region-threads>1 must equal the serial sweep in {ext}"
        );
    }
}
