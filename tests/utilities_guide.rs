//! `utilities guide` on a two-contig BAM built with samtools: a guide chain with a full-chain primary alignment and no
//! de novo model is rescued, an assembled one is not duplicated, secondary alignments do not count, and the output does
//! not depend on `--threads`.
use std::path::PathBuf;
use std::process::Command;

#[test]
fn rescues_guide_chains_with_a_full_chain_primary_alignment_at_any_thread_count() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("guide_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let ex = |c: &str, t: &str, a: u32, b: u32| format!("{c}\tr\texon\t{a}\t{b}\t.\t+\t.\tgene_id \"{t}\"; transcript_id \"{t}\";\n");
    // tG1 assembled (d1), tG2 only in reads, tG3 no read, tG4 only a secondary alignment
    let guide = [ex("c1", "tG1", 101, 200), ex("c1", "tG1", 301, 500), ex("c2", "tG2", 101, 200), ex("c2", "tG2", 401, 600),
        ex("c1", "tG3", 1001, 1100), ex("c1", "tG3", 1301, 1500), ex("c2", "tG4", 2001, 2100), ex("c2", "tG4", 2301, 2400)].concat();
    std::fs::write(dir.join("guide.gtf"), guide).unwrap();
    let d1 = "c1\trustle\ttranscript\t101\t500\t.\t+\t.\tgene_id \"d\"; transcript_id \"d1\"; reads \"5\";\n".to_string()
        + &ex("c1", "d1", 101, 200) + &ex("c1", "d1", 301, 500);
    std::fs::write(dir.join("denovo.gtf"), &d1).unwrap();
    let q = "*\t0\t0\t*\t*";
    let sam = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:10000\n@SQ\tSN:c2\tLN:10000\n\
         r1\t0\tc1\t121\t60\t80M100N100M\t{q}\nr2\t0\tc2\t151\t60\t50M200N50M\t{q}\nr3\t16\tc2\t151\t60\t50M200N50M\t{q}\n\
         r4\t256\tc2\t2051\t0\t50M200N50M\t{q}\n"
    );
    std::fs::write(dir.join("a.sam"), sam).unwrap();
    let bam = dir.join("a.bam");
    assert!(Command::new("samtools").args(["view", "-b", "-o"]).arg(&bam).arg(dir.join("a.sam")).status().unwrap().success());
    assert!(Command::new("samtools").arg("index").arg(&bam).status().unwrap().success());
    let run = |threads: &str| {
        let out = dir.join(format!("g{threads}.gtf"));
        let st = Command::new(env!("CARGO_BIN_EXE_utilities"))
            .args(["guide", "--threads", threads, "--gtf"]).arg(dir.join("denovo.gtf")).arg("--guide").arg(dir.join("guide.gtf"))
            .arg("--bam").arg(&bam).arg("--out").arg(&out).status().unwrap();
        assert!(st.success());
        std::fs::read_to_string(out).unwrap()
    };
    let g1 = run("1");
    assert_eq!(g1, run("4"));
    assert!(g1.starts_with(&d1), "the de novo model is kept verbatim");
    assert!(g1.contains("c2\tguide\ttranscript\t101\t600\t.\t+\t.\tgene_id \"tG2\"; transcript_id \"tG2\"; reads \"2\"; guide \"rescued\";"));
    for absent in ["\"tG1\"", "\"tG3\"", "\"tG4\""] {
        assert!(!g1.contains(absent), "{absent} must not be written");
    }
    let _ = std::fs::remove_dir_all(&dir);
}
