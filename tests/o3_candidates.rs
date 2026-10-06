use std::collections::BTreeMap;
use std::path::Path;
use std::process::{Command, Output};

const FX: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/o3_candidates");

/// The stage on the fixture's genome, copies and splice index (built in `dir`), with `bam` as the reads; the products' prefix is `dir/t.cand`.
fn run_stage(dir: &Path, bam: &str) -> (Output, String) {
    let mmi = dir.join("genome.splice.mmi");
    assert!(
        Command::new(std::env::var("RUSTLE_MINIMAP2").unwrap_or("minimap2".into()))
            .args(["-x", "splice", "-d"])
            .arg(&mmi)
            .arg(format!("{FX}/genome.fa"))
            .status()
            .unwrap()
            .success()
    );
    let out = dir.join("t.cand");
    let o = Command::new(env!("CARGO_BIN_EXE_o3_candidates"))
        .args([
            "--bam",
            bam,
            "--fasta",
            &format!("{FX}/genome.fa"),
            "--copies",
            &format!("{FX}/copies.tsv"),
            "--copies-fa",
            &format!("{FX}/copies.fa"),
            "--index",
            mmi.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--threads",
            "2",
        ])
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    (o, out.display().to_string())
}

/// The one flagged candidate's union: exactly one flagged row in candidates.tsv, and contigs.fa holds `cand_MCL0_0` alone.
fn flagged_union_len(out: &str) -> usize {
    let cands = std::fs::read_to_string(format!("{out}.candidates.tsv")).unwrap();
    let flagged: Vec<&str> = cands
        .lines()
        .skip(1)
        .filter(|l| l.split('\t').nth(4) == Some("1"))
        .collect();
    assert_eq!(flagged.len(), 1, "{cands}");
    let fa = std::fs::read_to_string(format!("{out}.contigs.fa")).unwrap();
    assert!(fa.starts_with(">cand_MCL0_0\n"), "{fa}");
    assert_eq!(fa.lines().count(), 2, "{fa}");
    fa.lines().nth(1).unwrap().len()
}

#[test]
fn flags_the_deleted_copy_and_assign_places_its_reads() {
    let dir = tempfile::tempdir().unwrap();
    let (_, out) = run_stage(dir.path(), &format!("{FX}/reads.bam"));
    let union_len = flagged_union_len(&out);
    assert!(
        (850..=950).contains(&union_len),
        "union {union_len} bp, expected the 900-bp spliced copy"
    );
}

fn samtools(args: &[&str]) -> Output {
    let o = Command::new("samtools").args(args).output().unwrap();
    assert!(
        o.status.success(),
        "samtools {args:?}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    o
}

/// The fixture BAM with pass B's attribution set in it (built here, the SAM text rewritten in Rust; only samtools is needed): the 60 copy-B
/// reads (`B_*`, make_fixture.py's names) as UNMAPPED records (flag 4; RNAME, POS, MAPQ, CIGAR and the tags dropped; SEQ/QUAL as stored,
/// which is as sequenced: every fixture record is a forward primary, asserted), plus four reads under the 300-bp floor (Amendment 13c): `U_00`,
/// the first 200 bases of B_05, unmapped; `P_00`..`P_02`, the first 250 bases of B_00..B_02, mapped far from the copy (chrT:40001, MAPQ 0,
/// de 0.05), so in no net and poorly placed. Sorted and indexed in `dir`.
fn bam_with_unmapped_copy_b(dir: &Path) -> String {
    let sam =
        String::from_utf8(samtools(&["view", "-h", &format!("{FX}/reads.bam")]).stdout).unwrap();
    let (mut text, mut b_seq) = (String::new(), BTreeMap::new());
    for line in sam.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if line.starts_with('@') || !f[0].starts_with("B_") {
            text += line;
            text += "\n";
            continue;
        }
        assert_eq!(
            f[1], "0",
            "{}: SEQ is kept as stored, so every copy-B record must be a forward primary",
            f[0]
        );
        text += &format!("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\n", f[0], f[9], f[10]);
        b_seq.insert(f[0].to_string(), f[9].to_string());
    }
    assert_eq!(b_seq.len(), 60, "the fixture holds 60 copy-B reads");
    assert!(
        b_seq.values().all(|s| s.len() >= 300),
        "every copy-B read is above the floor"
    );
    text += &format!(
        "U_00\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t*\n",
        &b_seq["B_05"][..200]
    );
    for k in 0..3 {
        text += &format!(
            "P_{k:02}\t0\tchrT\t40001\t0\t250M\t*\t0\t0\t{}\t*\tde:f:0.0500\n",
            &b_seq[&format!("B_{k:02}")][..250]
        );
    }
    let (sam_path, bam) = (dir.join("variant.sam"), dir.join("variant.bam"));
    std::fs::write(&sam_path, text).unwrap();
    samtools(&[
        "sort",
        "-o",
        bam.to_str().unwrap(),
        sam_path.to_str().unwrap(),
    ]);
    samtools(&["index", bam.to_str().unwrap()]);
    bam.display().to_string()
}

#[test]
fn pass_b_attributes_the_unmapped_copy_b_reads_and_keeps_reads_under_the_floor_out() {
    // prereg Amendments 13 / 13b / 13c at the BAM level: the deleted copy's 60 reads are unmapped, so only pass B's attribution (map-ont
    // against the net reads and the copy) can bring them into MCL0's net; the four reads under 300 bp (one unmapped, three poorly placed)
    // are never written to the attribution set, the poorly placed ones counted as below the floor
    let dir = tempfile::tempdir().unwrap();
    let bam = bam_with_unmapped_copy_b(dir.path());
    let (o, out) = run_stage(dir.path(), &bam);
    let log = String::from_utf8_lossy(&o.stderr);
    let pass_b = log
        .lines()
        .find(|l| l.starts_with("[o3_candidates] pass B: unmapped"))
        .unwrap_or_else(|| panic!("no pass B line: {log}"));
    for want in [
        "unmapped >= 300 bp 60;",
        "poorly placed >= 300 bp 0 (",
        "3 below the floor",
        "aligned 60: unmapped 60, poorly placed 0;",
        "attributed 60 (read coverage >= 0.5, de <= 0.20): unmapped 60, poorly placed 0;",
        "joined this run's families 60",
    ] {
        assert!(pass_b.contains(want), "`{want}` not in: {pass_b}");
    }
    // the net is the 60 copy-A reads of pass A and the 60 attributed copy-B reads, nothing under the floor
    let families = std::fs::read_to_string(format!("{out}.families.tsv")).unwrap();
    assert_eq!(
        families
            .lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .take(3)
            .collect::<Vec<_>>(),
        ["MCL0", "120", "120"],
        "{families}"
    );
    let union_len = flagged_union_len(&out);
    assert!(
        (850..=950).contains(&union_len),
        "union {union_len} bp, expected the 900-bp spliced copy"
    );
}
