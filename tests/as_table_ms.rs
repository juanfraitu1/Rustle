//! `as_table` reports, beside the best/second AS, the ms (max-scoring segment) scores, the primary's MAPQ and
//! splice shape, and the best-AS rival's AS, ms, splice shape and chromosome: the columns that size the O2
//! contested populations (near-tie on AS or on ms, spliced-vs-unspliced) without a second BAM pass.
use std::path::PathBuf;
use std::process::Command;

fn rec(name: &str, flag: u32, chrom: &str, pos: u64, mapq: u8, cigar: &str, as_: i32, ms: i32) -> String {
    let seq: String = cigar
        .split_inclusive(|c: char| !c.is_ascii_digit())
        .filter(|op| op.ends_with('M') || op.ends_with('S'))
        .map(|op| "A".repeat(op[..op.len() - 1].parse::<usize>().unwrap()))
        .collect();
    let qual = "I".repeat(seq.len());
    format!("{name}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tAS:i:{as_}\tms:i:{ms}")
}

#[test]
fn as_table_reports_ms_mapq_splice_and_rival_columns() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("as_table_ms_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let sam = dir.join("t.sam");
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:10000\n@SQ\tSN:c2\tLN:10000\n");
    for r in [
        // spliced primary on c1 (AS 300, ms 320, MAPQ 3); unspliced rival on c2 with higher AS but lower ms:
        // the processed-copy class. A weaker third secondary must not become the rival.
        rec("pseudo", 0, "c1", 100, 3, "50M500N50M", 300, 320),
        rec("pseudo", 256, "c2", 100, 0, "100M", 310, 300),
        rec("pseudo", 256, "c2", 5000, 0, "100M", 250, 250),
        // single placement: no rival at all
        rec("unique", 0, "c1", 2000, 60, "100M", 400, 400),
        // exact tie on both scores, both unspliced; supplementary record ignored
        rec("tie", 0, "c1", 3000, 0, "100M", 200, 200),
        rec("tie", 256, "c1", 6000, 0, "100M", 200, 200),
        rec("tie", 2048, "c1", 8000, 0, "20M", 50, 50),
    ] {
        s.push_str(&r);
        s.push('\n');
    }
    std::fs::write(&sam, s).unwrap();
    let bam = dir.join("t.bam");
    assert!(Command::new("samtools").args(["view", "-b", "-o"]).arg(&bam).arg(&sam).status().unwrap().success());
    let out = dir.join("t.molecules.tsv");
    let st = Command::new(env!("CARGO_BIN_EXE_as_table"))
        .env("RUSTLE_AS_TABLE_SIDECAR", "0")
        .arg("--bam")
        .arg(&bam)
        .arg("--out")
        .arg(&out)
        .status()
        .unwrap();
    assert!(st.success());
    let text = std::fs::read_to_string(&out).unwrap();
    let rows: Vec<Vec<&str>> = text.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').collect()).collect();
    // name best_as second_as n primary_chrom primary_as | best_ms second_ms primary_ms primary_mapq primary_spliced
    // rival_as rival_ms rival_spliced rival_chrom
    let want: Vec<Vec<&str>> = vec![
        vec!["pseudo", "310", "300", "3", "c1", "300", "320", "300", "320", "3", "1", "310", "300", "0", "c2"],
        vec!["tie", "200", "200", "2", "c1", "200", "200", "200", "200", "0", "0", "200", "200", "0", "c1"],
        vec!["unique", "400", "-1", "1", "c1", "400", "400", "-1", "400", "60", "0", "NA", "NA", "NA", "NA"],
    ];
    assert_eq!(rows, want, "table:\n{text}");
}
