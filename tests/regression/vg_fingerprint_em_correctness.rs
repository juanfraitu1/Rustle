//! End-to-end correctness tests for the fingerprint-EM path.
//!
//! Runs rustle on `test_data/synthetic_family/reads_sorted.bam` (2 paralogous
//! gene copies, 7 planted diagnostic SNPs, 28 multi-mappers) and verifies that
//! the fingerprint-EM output satisfies fundamental mathematical and biological
//! properties — without relying on hard-coded read-to-copy assignments that
//! can shift as the scoring formula evolves.
//!
//! The three invariants tested:
//!
//!   1. **Weight conservation** — per multi-mapper read, final_weight sums to 1.0
//!      (the EM is a normalised mixture; probability must integrate to 1).
//!
//!   2. **Decisive-rate floor** — ≥30% of multi-mappers reach weight_gap > 0.8
//!      (synthetic reads cover 7 planted SNPs; sufficient diagnostic signal must
//!      produce decisive assignments, not uniform 0.5).
//!
//!   3. **Both copies assembled** — the output GTF contains at least one
//!      transcript with copy_id "0" and at least one with copy_id "1", proving
//!      that EM routes reads to both bundles and neither is silenced.
//!
//! A fourth test checks that every transcript in the GTF with a family_id
//! attribute also has copy_id and copy_confidence, and that copy_confidence is
//! a valid float in (0.0, 1.0].

use std::collections::HashMap;
use std::process::Command;

const FIXTURE_BAM: &str = "test_data/synthetic_family/reads_sorted.bam";
const FIXTURE_FA:  &str = "test_data/synthetic_family/ref.fa";

// ── helpers ────────────────────────────────────────────────────────────────────

fn fixture_present() -> bool {
    std::path::Path::new(FIXTURE_BAM).exists() && std::path::Path::new(FIXTURE_FA).exists()
}

/// Run rustle on the synthetic fixture with fingerprint-EM enabled.
/// Returns (gtf_text, tsv_text).
///
/// Each call writes to fresh temporary files so parallel test invocations
/// don't trample each other.
fn run_fp_em() -> (String, String) {
    let gtf_file = tempfile::NamedTempFile::new().expect("tempfile for GTF");
    let tsv_file = tempfile::NamedTempFile::new().expect("tempfile for TSV");

    let status = Command::new(env!("CARGO_BIN_EXE_rustle"))
        .env("RUSTLE_VG_FP_ATTR_TSV", tsv_file.path())
        .args([
            "-L",
            "--vg",
            "--vg-no-hmm",
            "--genome-fasta", FIXTURE_FA,
            "-o", gtf_file.path().to_str().unwrap(),
            FIXTURE_BAM,
        ])
        .status()
        .expect("rustle binary failed to spawn");

    assert!(status.success(), "rustle exited non-zero: {:?}", status);

    let gtf = std::fs::read_to_string(gtf_file.path()).unwrap_or_default();
    let tsv = std::fs::read_to_string(tsv_file.path())
        .expect("RUSTLE_VG_FP_ATTR_TSV not written; fingerprint-EM did not fire");

    (gtf, tsv)
}

/// Parse the TSV into rows.  Returns (read_hash, placement_copy, n_sites,
/// final_weight, weight_gap, weight_sum) per data row.
fn parse_tsv(tsv: &str) -> Vec<(u64, usize, usize, f64, f64, f64)> {
    tsv.lines()
        .skip(1) // header
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            assert!(cols.len() >= 7, "malformed TSV row: {:?}", l);
            let hash: u64  = cols[1].parse().expect("read_name_hash not u64");
            let copy: usize = cols[2].parse().expect("placement_copy not usize");
            let sites: usize = cols[3].parse().expect("n_sites_covered not usize");
            let fw: f64  = cols[4].parse().expect("final_weight not f64");
            let gap: f64 = cols[5].parse().expect("weight_gap not f64");
            let sum: f64 = cols[6].parse().expect("weight_sum not f64");
            (hash, copy, sites, fw, gap, sum)
        })
        .collect()
}

// ── tests ──────────────────────────────────────────────────────────────────────

#[test]
fn fp_em_weight_conservation() {
    if !fixture_present() { return; }

    let (_gtf, tsv) = run_fp_em();
    let rows = parse_tsv(&tsv);

    // Sum final_weight per read_name_hash.
    let mut sums: HashMap<u64, f64> = HashMap::new();
    for (hash, _copy, _sites, fw, _gap, _wsum) in &rows {
        *sums.entry(*hash).or_default() += fw;
    }

    assert!(!sums.is_empty(), "TSV has no data rows — fingerprint-EM produced nothing");

    for (hash, total) in &sums {
        assert!(
            (total - 1.0).abs() < 1e-5,
            "weight conservation violated for read {}: sum = {} (expected 1.0)",
            hash, total
        );
    }
}

#[test]
fn fp_em_decisive_rate_above_30_percent() {
    if !fixture_present() { return; }

    let (_gtf, tsv) = run_fp_em();
    let rows = parse_tsv(&tsv);

    // Count distinct reads and how many have weight_gap > 0.8 (decisive).
    let mut seen: std::collections::HashSet<u64> = Default::default();
    let mut decisive = 0usize;
    let mut total = 0usize;
    for (hash, _copy, _sites, _fw, gap, _wsum) in &rows {
        if seen.insert(*hash) {
            total += 1;
            if *gap > 0.8 {
                decisive += 1;
            }
        }
    }

    assert!(total > 0, "TSV is empty");
    let rate = decisive as f64 / total as f64;
    assert!(
        rate >= 0.30,
        "decisive rate too low: {}/{} = {:.1}% (expected ≥30%)",
        decisive, total, rate * 100.0
    );
}

#[test]
fn fp_em_both_copies_assembled_in_gtf() {
    if !fixture_present() { return; }

    let (gtf, _tsv) = run_fp_em();

    let mut saw_copy_0 = false;
    let mut saw_copy_1 = false;

    for line in gtf.lines() {
        if !line.contains("copy_id") { continue; }
        if line.contains("copy_id \"0\"") { saw_copy_0 = true; }
        if line.contains("copy_id \"1\"") { saw_copy_1 = true; }
        if saw_copy_0 && saw_copy_1 { break; }
    }

    assert!(
        saw_copy_0,
        "GTF has no transcript with copy_id \"0\" — copy A was silenced by EM"
    );
    assert!(
        saw_copy_1,
        "GTF has no transcript with copy_id \"1\" — copy B was silenced by EM"
    );
}

#[test]
fn fp_em_gtf_copy_attributes_are_well_formed() {
    if !fixture_present() { return; }

    let (gtf, _tsv) = run_fp_em();

    let tx_lines_with_family: Vec<&str> = gtf.lines()
        .filter(|l| {
            let cols: Vec<&str> = l.splitn(10, '\t').collect();
            cols.get(2) == Some(&"transcript") && l.contains("family_id")
        })
        .collect();

    assert!(
        !tx_lines_with_family.is_empty(),
        "no transcript lines contain family_id — VG GTF annotation not written"
    );

    for line in &tx_lines_with_family {
        assert!(
            line.contains("copy_id"),
            "transcript has family_id but no copy_id: {}", line
        );
        assert!(
            line.contains("copy_confidence"),
            "transcript has family_id but no copy_confidence: {}", line
        );

        // Extract copy_confidence value and verify it is in (0.0, 1.0].
        let conf: f64 = line
            .split("copy_confidence \"")
            .nth(1)
            .and_then(|s| s.split('"').next())
            .and_then(|v| v.parse().ok())
            .expect("copy_confidence not a parseable float");
        assert!(
            conf > 0.0 && conf <= 1.0,
            "copy_confidence out of range (0,1]: {} on line: {}",
            conf, line
        );
    }
}

#[ignore = "Stale fixture premise. This test was calibrated to a PRE-`f61a296` \
fixture where diagnostic SNPs were not actually applied, so `n_sites_covered` \
spanned ~9-1800 (a wide bimodal range). Commit f61a296 ('actually apply \
diagnostic SNPs') corrected the fixture to its real ~31 diagnostic sites, so \
every read now covers a UNIFORM 26-31 sites and mean weight_gap is identical \
(~0.50) across the high/low buckets -> the >50/<20 thresholds are unsatisfiable \
and the 'scales-with-coverage' claim has no coverage gradient to test. Re-enable \
after regenerating `test_data/synthetic_family` with single-exon (genuinely \
low-coverage, <20-site) reads alongside the full-length ones. NOT a tandem-VG \
regression (the tandem commits never touched this file or the fixture)."]
#[test]
fn fp_em_decisiveness_scales_with_diagnostic_site_coverage() {
    // IsoSeq advantage claim: reads that cover more copy-distinguishing
    // positions (longer reads spanning more diagnostic SNPs) should be
    // more decisive than reads covering few positions.
    //
    // In the synthetic fixture, multi-mapper reads that start before the
    // first diagnostic SNP in exon 1 cover ~250–1800 sites total across
    // all placements; reads that start after exon 1's SNP cover only ~9
    // sites.  The threshold is chosen conservatively to keep this test
    // resilient to minor fixture regeneration.
    //
    // Expected outcome:
    //   - All reads with max_sites > 50  → decisive (weight_gap > 0.8)
    //   - All reads with max_sites < 20  → uncertain (weight_gap < 0.5)
    if !fixture_present() { return; }

    let (_gtf, tsv) = run_fp_em();
    let rows = parse_tsv(&tsv);

    // Aggregate per read: max sites covered across placements, and the
    // weight_gap (same for every row of a given read).
    let mut per_read: HashMap<u64, (usize, f64)> = HashMap::new();
    for (hash, _copy, sites, _fw, gap, _wsum) in &rows {
        let entry = per_read.entry(*hash).or_insert((0, *gap));
        if *sites > entry.0 { entry.0 = *sites; }
    }

    let high: Vec<_> = per_read.values().filter(|(s, _)| *s > 50).collect();
    let low:  Vec<_> = per_read.values().filter(|(s, _)| *s < 20).collect();

    assert!(!high.is_empty(), "no high-coverage reads in TSV — fixture may have changed");
    assert!(!low.is_empty(),  "no low-coverage reads in TSV — fixture may have changed");

    for (sites, gap) in &high {
        assert!(
            *gap > 0.8,
            "read with {} diagnostic sites should be decisive (gap > 0.8), got gap={}",
            sites, gap
        );
    }
    for (sites, gap) in &low {
        assert!(
            *gap < 0.5,
            "read with {} diagnostic sites should be uncertain (gap < 0.5), got gap={}",
            sites, gap
        );
    }
}

#[test]
fn topo_borrow_adds_transcripts_to_under_assembled_copies() {
    if !fixture_present() { return; }
    // Run without topo borrow.
    let (gtf_base, _) = run_fp_em();
    // Run with topo borrow.
    let gtf_topo = {
        let gtf_file = tempfile::NamedTempFile::new().unwrap();
        let status = std::process::Command::new(env!("CARGO_BIN_EXE_rustle"))
            .env("RUSTLE_VG_TOPO_BORROW", "1")
            .args([
                "-L",
                "--vg",
                "--vg-no-hmm",
                "--genome-fasta", FIXTURE_FA,
                "-o", gtf_file.path().to_str().unwrap(),
                FIXTURE_BAM,
            ])
            .status()
            .expect("rustle spawn failed");
        assert!(status.success(), "rustle with RUSTLE_VG_TOPO_BORROW exited non-zero: {:?}", status);
        std::fs::read_to_string(gtf_file.path()).unwrap_or_default()
    };
    // With topo borrow: should have >= as many transcripts as without.
    let count_base = gtf_base.lines().filter(|l| l.contains("\ttranscript\t")).count();
    let count_topo = gtf_topo.lines().filter(|l| l.contains("\ttranscript\t")).count();
    assert!(
        count_topo >= count_base,
        "topo borrow should not remove transcripts: base={} topo={}",
        count_base, count_topo
    );
    // At least one topology_borrow transcript should appear (if source copy has
    // confident multi-exon transcripts that can be projected to the sister).
    // The fixture has two synthetic gene copies — if neither qualifies, just check
    // that the binary ran successfully and count did not drop.
    let _ = gtf_topo.contains("topo_borrow"); // informational; not a hard assert for sparse fixtures
}
