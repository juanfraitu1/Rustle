//! `utilities sample-counts` on a small BAM built with samtools: only primary mapped reads count, strand is not a filter, and
//! each sample's column sums to its compatible reads.
use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_utilities"))
}

fn bam(dir: &PathBuf, name: &str, recs: &[&str]) -> PathBuf {
    let sam = dir.join(format!("{name}.sam"));
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:10000\n");
    for r in recs {
        s.push_str(r);
        s.push('\n');
    }
    std::fs::write(&sam, s).unwrap();
    let b = dir.join(format!("{name}.bam"));
    assert!(Command::new("samtools").args(["view", "-b", "-o"]).arg(&b).arg(&sam).status().unwrap().success());
    b
}

#[test]
fn only_primary_reads_count_and_columns_sum_to_compatible() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("sample_counts_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // t1: exons 101-200, 301-500 (1-based); t2 single exon 1001-1500
    std::fs::write(
        dir.join("f.gtf"),
        "c1\tr\ttranscript\t101\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
         c1\tr\texon\t101\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
         c1\tr\texon\t301\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
         c1\tr\ttranscript\t1001\t1500\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\";\n\
         c1\tr\texon\t1001\t1500\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\";\n",
    )
    .unwrap();
    let q = "*\t0\t0\t*\t*";
    let a = bam(&dir, "a", &[
        &format!("r1\t0\tc1\t121\t60\t80M100N100M\t{q}"),        // t1, primary
        &format!("r2\t16\tc1\t121\t60\t80M100N100M\t{q}"),       // t1, reverse strand: still counts
        &format!("r3\t256\tc1\t121\t60\t80M100N100M\t{q}"),      // secondary: never
        &format!("r4\t2048\tc1\t121\t60\t80M100N100M\t{q}"),     // supplementary: never
        &format!("r6\t0\tc1\t1101\t60\t200M\t{q}"),              // t2
        &format!("r5\t4\t*\t0\t0\t*\t{q}"),                      // unmapped: never
    ]);
    let b = bam(&dir, "b", &[&format!("s1\t0\tc1\t5001\t60\t50M\t{q}")]); // no transcript
    let out = dir.join("counts.tsv");
    let o = Command::new(bin())
        .args(["sample-counts", "--gtf"]).arg(dir.join("f.gtf"))
        .arg("--bam").arg(format!("{},{}", a.display(), b.display()))
        .args(["--samples", "A,B", "--out"]).arg(&out)
        .output().unwrap();
    let log = String::from_utf8_lossy(&o.stderr);
    assert!(o.status.success(), "{log}");
    assert!(log.contains("A: seen 3, compatible 3, unassigned 0"), "{log}");
    assert!(log.contains("B: seen 1, compatible 0, unassigned 1"), "{log}");
    let t = std::fs::read_to_string(&out).unwrap();
    assert_eq!(t, "transcript_id\tgene_id\tA\tB\nt1\tg1\t2.00\t0.00\nt2\tg2\t1.00\t0.00\n");
}
