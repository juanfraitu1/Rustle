//! r1507: `BamIndexCache::sweep_windows` (one query per chromosome, records bucketed into windows) returns, per
//! window, exactly what `reads_in_region_mode` returns for that window — same records, same order — in every mode.
use rustle::denovo_assemble::{BamIndexCache, ReadMode};
use rustle::types::DetHashSet;
use std::path::PathBuf;
use std::process::Command;

fn bam(dir: &PathBuf, recs: &[&str]) -> String {
    let sam = dir.join("w.sam");
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:10000\n");
    for r in recs {
        s.push_str(r);
        s.push('\n');
    }
    std::fs::write(&sam, s).unwrap();
    let b = dir.join("w.bam");
    assert!(Command::new("samtools").args(["view", "-b", "-o"]).arg(&b).arg(&sam).status().unwrap().success());
    assert!(Command::new("samtools").arg("index").arg(&b).status().unwrap().success());
    b.to_string_lossy().into_owned()
}

fn rec(name: &str, flag: u32, pos: u64, cigar: &str) -> String {
    let n = cigar.split(|c: char| !c.is_ascii_digit()).filter(|x| !x.is_empty()).count();
    let seq: String = cigar
        .split_inclusive(|c: char| !c.is_ascii_digit())
        .filter(|op| op.ends_with('M') || op.ends_with('S'))
        .map(|op| "A".repeat(op[..op.len() - 1].parse::<usize>().unwrap()))
        .collect();
    let _ = n;
    let (seq, qual) = if cigar == "*" { ("*".to_string(), "*".to_string()) } else { (seq.clone(), "I".repeat(seq.len())) };
    format!("{name}\t{flag}\tc1\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tAS:i:90\tNM:i:1\tde:f:0.01")
}

#[test]
fn sweep_windows_matches_the_per_window_query_in_every_mode() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("window_sweep_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // windows (0-based half-open): W0 [1000,2000), W1 [1500,3000), W2 [2000,2500) (W1 holds it; starts where W0 ends),
    // W3 = W0 again; given out of order
    let path = bam(
        &dir,
        &[
            &rec("r1_before", 0, 901, "50M"),                 // ends before every window
            &rec("r2_ends_at_lo", 0, 951, "50M"),             // 1-based 951-1000: ends exactly at W0's lo -> none
            &rec("r4_spliced", 0, 1200, "100M1500N100M"),     // 1200-2899: W0, W1, W2 (intron over W2)
            &rec("r11_no_cigar", 0, 1600, "*"),               // placed, no alignment
            &rec("r6_secondary", 256, 1600, "50M"),
            &rec("r7_supplementary", 2048, 1700, "50M"),
            &rec("r3_ends_at_hi", 0, 1951, "50M"),            // 1951-2000: W0 (hi 2000), W1; not W2 (lo 2000 -> 2001)
            &rec("r12_starts_on_w0_last_base", 0, 2000, "50M"), // 2000-2049: W0 (its last base), W1, W2
            &rec("r5_after_w0", 0, 2001, "100M"),             // W1, W2
            &rec("r9_crosses_w2_end", 0, 2490, "50M"),        // 2490-2539: W1, W2
            &rec("r10_reverse_last", 16, 2999, "50M"),        // 2999-3048: W1 only (hi 3000)
        ],
    );
    let cache = BamIndexCache::open(&path).unwrap();
    let only: DetHashSet<String> = ["r4_spliced", "r9_crosses_w2_end"].iter().map(|s| s.to_string()).collect();
    let wins = [(2000u64, 2500u64), (1000, 2000), (1500, 3000), (1000, 2000)];
    for mode in [ReadMode::Light, ReadMode::Full, ReadMode::Only(&only)] {
        let spec: Vec<(u64, u64, ReadMode<'_>)> = wins.iter().map(|w| (w.0, w.1, mode)).collect();
        let swept = cache.sweep_windows(&path, "c1", &spec).unwrap();
        assert_eq!(swept.len(), wins.len());
        let mut n_total = 0;
        for (i, w) in wins.iter().enumerate() {
            let (p, b) = cache.reads_in_region_mode(&path, "c1", w.0, w.1, mode).unwrap();
            assert_eq!(swept[i].0, p, "primaries of window {i} {w:?}");
            assert_eq!(format!("{:?}", swept[i].1), format!("{:?}", b), "bam reads of window {i} {w:?}");
            n_total += b.len();
        }
        assert!(n_total > 0, "the fixture yields records");
    }
    // the per-window query itself holds the expected sets (the fixture is as described)
    let names = |lo, hi| -> Vec<String> {
        cache.reads_in_region_mode(&path, "c1", lo, hi, ReadMode::Full).unwrap().1.into_iter().map(|b| b.name).collect()
    };
    assert_eq!(names(1000, 2000), ["r4_spliced", "r6_secondary", "r7_supplementary", "r3_ends_at_hi", "r12_starts_on_w0_last_base"]);
    assert_eq!(names(2000, 2500), ["r4_spliced", "r12_starts_on_w0_last_base", "r5_after_w0", "r9_crosses_w2_end"]);
    let w1 = names(1500, 3000);
    assert!(w1.contains(&"r10_reverse_last".to_string()) && !w1.contains(&"r2_ends_at_lo".to_string()));
}
