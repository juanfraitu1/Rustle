//! r1600 task 2: `BamIndexCache::convert_record` applied to a window's query records gives exactly what
//! `reads_in_region_mode` gives for that window, in every mode, on the fixture BAMs.
use rustle::denovo_assemble::{BamIndexCache, ReadMode};
use rustle::types::DetHashSet;

const FIXTURES: [&str; 4] = ["xfam_conflict", "o2_origin_resolve", "candidate_copies", "same_chrom_supplement"];

#[test]
fn convert_record_equals_reads_in_region() {
    let mut compared = 0usize;
    for fx in FIXTURES {
        let path = format!("{}/tests/fixtures/{fx}/reads.bam", env!("CARGO_MANIFEST_DIR"));
        let cache = BamIndexCache::open(&path).unwrap();
        let contigs: Vec<(String, u64)> = cache
            .header()
            .reference_sequences()
            .iter()
            .map(|(n, m)| (n.to_string(), m.length().get() as u64))
            .collect();
        for (chrom, len) in &contigs {
            for (lo, hi) in [(0u64, *len), (100, 2000), (0, 1)] {
                let first_two: DetHashSet<String> = {
                    let mut r = noodles_bam::io::reader::Builder::default().build_from_path(&path).unwrap();
                    let region: noodles_core::Region = format!("{chrom}:{}-{hi}", lo + 1).parse().unwrap();
                    r.query(cache.header(), cache.index(), &region)
                        .unwrap()
                        .take(2)
                        .map(|x| String::from_utf8_lossy(x.unwrap().name().unwrap().as_ref()).into_owned())
                        .collect()
                };
                for mode in [ReadMode::Full, ReadMode::Light, ReadMode::Only(&first_two)] {
                    let want = cache.reads_in_region_mode(&path, chrom, lo, hi, mode).unwrap();
                    let mut r = noodles_bam::io::reader::Builder::default().build_from_path(&path).unwrap();
                    let region: noodles_core::Region = format!("{chrom}:{}-{hi}", lo + 1).parse().unwrap();
                    let (mut p, mut b) = (Vec::new(), Vec::new());
                    for rec in r.query(cache.header(), cache.index(), &region).unwrap() {
                        let (pr, br) = cache.convert_record(&rec.unwrap(), chrom, mode).unwrap();
                        p.extend(pr);
                        b.extend(br);
                    }
                    compared += want.1.len();
                    assert_eq!(format!("{:?}", (p, b)), format!("{want:?}"), "{fx} {chrom}:{lo}-{hi}");
                }
            }
        }
    }
    assert!(compared > 100, "only {compared} records compared");
    eprintln!("compared {compared} records");
}
#[test]
/// The query path's output on the fixtures, pinned before the split (digest computed on 5276a03e): the split changed nothing.
fn reads_in_region_output_is_unchanged_by_the_split() {
    let mut h: u64 = 0xcbf29ce484222325;
    for fx in FIXTURES {
        let path = format!("{}/tests/fixtures/{fx}/reads.bam", env!("CARGO_MANIFEST_DIR"));
        let cache = BamIndexCache::open(&path).unwrap();
        let mut r = noodles_bam::io::reader::Builder::default().build_from_path(&path).unwrap();
        let header = r.read_header().unwrap();
        for (n, m) in header.reference_sequences() {
            let (chrom, len) = (n.to_string(), m.length().get() as u64);
            for (lo, hi) in [(0u64, len), (100, 2000), (0, 1)] {
                let only: DetHashSet<String> = cache.reads_in_region_mode(&path, &chrom, lo, hi, ReadMode::Full).unwrap().1.iter().take(2).map(|b| b.name.clone()).collect();
                for mode in [ReadMode::Full, ReadMode::Light, ReadMode::Only(&only)] {
                    for b in format!("{:?}", cache.reads_in_region_mode(&path, &chrom, lo, hi, mode).unwrap()).bytes() {
                        h = (h ^ b as u64).wrapping_mul(0x100000001b3);
                    }
                }
            }
        }
    }
    assert_eq!(h, 0xf84303da1d238089, "reads_in_region_mode output changed: {h:#x}");
}
