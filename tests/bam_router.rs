//! r1600 task 3: `bam_router::route` gives every window exactly its own noodles query's records, in file order, at any
//! worker count, converting each routed record once.
use noodles_sam::alignment::RecordBuf;
use rustle::bam_router::route;
use rustle::denovo_assemble::BamIndexCache;
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

fn dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("bam_router_{name}_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A coordinate-sorted, indexed BAM from SAM lines (samtools, as tests/window_sweep.rs).
fn bam(d: &PathBuf, recs: &[String]) -> String {
    let sam = d.join("r.sam");
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:100000\n@SQ\tSN:c2\tLN:100000\n");
    for r in recs {
        s.push_str(r);
        s.push('\n');
    }
    std::fs::write(&sam, s).unwrap();
    let b = d.join("r.bam");
    assert!(Command::new("samtools").args(["sort", "-o"]).arg(&b).arg(&sam).status().unwrap().success());
    assert!(Command::new("samtools").arg("index").arg(&b).status().unwrap().success());
    b.to_string_lossy().into_owned()
}

fn rec(name: &str, flag: u32, chrom: &str, pos: u64, cigar: &str) -> String {
    let len: usize = cigar
        .split_inclusive(|c: char| !c.is_ascii_digit())
        .filter(|op| op.ends_with('M') || op.ends_with('S') || op.ends_with('I'))
        .map(|op| op[..op.len() - 1].parse::<usize>().unwrap())
        .sum();
    let (seq, qual) = if cigar == "*" { ("ACGT".to_string(), "IIII".to_string()) } else { ("ACGT".repeat(len / 4 + 1)[..len].to_string(), "I".repeat(len)) };
    format!("{name}\t{flag}\t{chrom}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tAS:i:90")
}

fn render(cache: &BamIndexCache, r: &noodles_bam::Record) -> String {
    format!("{:?}", RecordBuf::try_from_alignment_record(cache.header(), r).unwrap())
}

/// Each window's own noodles query, rendered.
fn queried(path: &str, cache: &BamIndexCache, windows: &[(String, u64, u64)]) -> Vec<Vec<String>> {
    windows
        .iter()
        .map(|(c, lo, hi)| {
            let mut r = noodles_bam::io::reader::Builder::default().build_from_path(path).unwrap();
            let region: noodles_core::Region = format!("{c}:{}-{hi}", lo + 1).parse().unwrap();
            r.query(cache.header(), cache.index(), &region).unwrap().map(|x| render(cache, &x.unwrap())).collect()
        })
        .collect()
}

fn routed(path: &str, cache: &BamIndexCache, windows: &[(String, u64, u64)], workers: usize) -> Vec<Vec<String>> {
    route(path, cache, windows, workers, &|r: &noodles_bam::Record, _w: &[usize]| Ok(Some(render(cache, r)))).unwrap()
}

fn w(c: &str, lo: u64, hi: u64) -> (String, u64, u64) {
    (c.to_string(), lo, hi)
}

fn edge_bam(d: &PathBuf) -> String {
    bam(
        d,
        &[
            rec("a_1000_1999", 0, "c1", 1001, "1000M"),     // 0-based [1000, 2000)
            rec("b_2000_2999", 0, "c1", 2001, "1000M"),     // [2000, 3000)
            rec("c_spliced", 0, "c1", 1501, "100M7300N100M"), // [1500, 9000)
            rec("d_secondary", 256, "c1", 1801, "300M"),
            rec("e_unmapped_placed", 4, "c1", 1501, "*"),
            rec("f_other_contig", 0, "c2", 1501, "200M"),
        ],
    )
}

fn edge_windows() -> Vec<(String, u64, u64)> {
    vec![
        w("c1", 999, 1000),   // the base before a_1000_1999
        w("c1", 1000, 1001),  // a's first base
        w("c1", 1999, 2000),  // a's last base
        w("c1", 2000, 2001),  // b's first base
        w("c1", 3000, 3001),  // the base after b: only c_spliced (intron)
        w("c1", 0, 20000),
        w("c2", 1500, 1501),
        w("c1", 1500, 1500),  // empty
    ]
}

#[test]
fn routed_equals_query_on_edge_windows() {
    let d = dir("edge");
    let path = edge_bam(&d);
    let cache = BamIndexCache::open(&path).unwrap();
    let wins = edge_windows();
    let want = queried(&path, &cache, &wins);
    assert!(want.iter().map(|v| v.len()).sum::<usize>() >= 8, "fixture too thin: {want:?}");
    assert_eq!(routed(&path, &cache, &wins, 1), want);
}

#[test]
fn routed_is_identical_at_1_2_7_and_8_workers() {
    let d = dir("workers");
    // 5,000 records, fixed seed (LCG), over both contigs: ~80 BGZF blocks of 64 KiB uncompressed
    let mut x: u64 = 0x2545F4914F6CDD1D;
    let mut next = || {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        x >> 33
    };
    // 30 long spliced records at c1's start sit in a top-level BAI bin, so a window far downstream reads two separate
    // file ranges (that bin's chunk at the file start, then its local bins): its list spans worker groups
    let mut recs: Vec<String> = (0..30).map(|i| rec(&format!("long{i}"), 0, "c1", 1 + i, "100M89800N100M")).collect();
    recs.extend((0..5000)
        .map(|i| {
            let c = if next() % 3 == 0 { "c2" } else { "c1" };
            let pos = 1 + next() % 95000;
            let cig = if next() % 4 == 0 { format!("{}M{}N{}M", 50 + next() % 400, 100 + next() % 3000, 50 + next() % 400) } else { format!("{}M", 100 + next() % 900) };
            rec(&format!("r{i}"), if next() % 5 == 0 { 256 } else { 0 }, c, pos, &cig)
        }));
    let path = bam(&d, &recs);
    let cache = BamIndexCache::open(&path).unwrap();
    let mut wins: Vec<(String, u64, u64)> = (0..60).map(|_| {
        let c = if next() % 3 == 0 { "c2" } else { "c1" };
        let lo = next() % 98000;
        w(c, lo, lo + 1 + next() % 2000)
    }).collect();
    wins.push(w("c1", 89900, 89950)); // the long records' second exons and local records
    let want = queried(&path, &cache, &wins);
    assert!(want.iter().map(|v| v.len()).sum::<usize>() > 1000);
    for k in [1, 2, 7, 8] {
        assert_eq!(routed(&path, &cache, &wins, k), want, "{k} workers");
    }
}

#[test]
fn duplicate_windows_each_get_the_full_list() {
    let d = dir("dup");
    let path = edge_bam(&d);
    let cache = BamIndexCache::open(&path).unwrap();
    let wins = vec![w("c1", 1500, 2500), w("c1", 1500, 2500)];
    let got = routed(&path, &cache, &wins, 2);
    assert!(!got[0].is_empty());
    assert_eq!(got[0], got[1]);
}

#[test]
fn unknown_contig_fails_naming_the_window() {
    let d = dir("unknown");
    let path = edge_bam(&d);
    let cache = BamIndexCache::open(&path).unwrap();
    let e = route(&path, &cache, &[w("c9", 0, 100)], 1, &|_r: &noodles_bam::Record, _w: &[usize]| Ok(Some(()))).unwrap_err();
    assert!(format!("{e:#}").contains("c9"), "{e:#}");
}

#[test]
fn convert_is_called_once_per_record() {
    let d = dir("once");
    let path = edge_bam(&d);
    let cache = BamIndexCache::open(&path).unwrap();
    let wins = edge_windows();
    let calls = AtomicUsize::new(0);
    let got = route(&path, &cache, &wins, 3, &|r: &noodles_bam::Record, ws: &[usize]| {
        calls.fetch_add(1, Ordering::Relaxed);
        assert!(!ws.is_empty() && ws.windows(2).all(|p| p[0] < p[1]), "windows sorted, non-empty: {ws:?}");
        Ok(Some(render(&cache, r)))
    })
    .unwrap();
    let mut distinct: Vec<&String> = got.iter().flatten().collect();
    distinct.sort();
    distinct.dedup();
    assert_eq!(calls.load(Ordering::Relaxed), distinct.len());
}

/// One window whose records lie in two separate file ranges: 400 long spliced records at c1's start (a 128-kb BAI bin,
/// at the file start) and 30 local records at 89,800, with ~3,000 short unrelated records between them in the file.
/// At 2 workers the two ranges go to different groups, so the window's list is joined across groups in file order.
#[test]
fn a_window_split_across_groups_keeps_file_order() {
    let d = dir("split");
    // random bases and qualities: incompressible, so bins keep real compressed sizes (htslib folds a bin under 64 kb into
    // its parent, which would make the window one range) and the long records' range outweighs the local one, so 2
    // workers cut between the two ranges
    let mut x: u64 = 7;
    let mut rnd = |n: usize, alpha: &[u8]| -> String {
        (0..n).map(|_| { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); alpha[(x >> 33) as usize % alpha.len()] as char }).collect()
    };
    const Q: &[u8] = b"!#%&()*+-/0123456789:;<=>?@ABCDEFGHI";
    let mut recs: Vec<String> = (0..400)
        .map(|i| {
            let (seq, qual) = (rnd(400, b"ACGT"), rnd(400, Q));
            format!("long{i}\t0\tc1\t{}\t60\t200M89500N200M\t*\t0\t0\t{seq}\t{qual}\tAS:i:90", 1 + i)
        })
        .collect();
    // short records kept inside 16-kb tiles (no 128-kb bin entries), between the long records and the window
    recs.extend((0..3000).map(|i| {
        let i = i as u64;
        let pos = 16384 * (2 + i % 3) + 101 + (i % 50) * 200; // offset + 300 < 16384: never crosses a tile
        let (seq, qual) = (rnd(300, b"ACGT"), rnd(300, Q));
        format!("mid{i}\t0\tc1\t{pos}\t60\t300M\t*\t0\t0\t{seq}\t{qual}\tAS:i:90")
    }));
    recs.extend((0..30).map(|i| rec(&format!("local{i}"), 0, "c1", 89801 + i, "400M")));
    let path = bam(&d, &recs);
    let cache = BamIndexCache::open(&path).unwrap();
    let wins = vec![w("c1", 89900, 89950)];
    let region: noodles_core::Region = "c1:89901-89950".parse().unwrap();
    use noodles_csi::BinningIndex;
    let chunks = cache.index().query(0, region.interval()).unwrap();
    assert!(chunks.len() >= 2, "the window must read separate ranges: {chunks:?}");
    let want = queried(&path, &cache, &wins);
    assert_eq!(want[0].len(), 429); // 399 long (long0 ends at 89,900, one base before the window) + 30 local
    for k in [1, 2, 3] {
        assert_eq!(routed(&path, &cache, &wins, k), want, "{k} workers");
    }
}
