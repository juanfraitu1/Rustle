//! `family_gfa` end to end on a synthetic families-stage run: FAM1 = A (+), B (+, one SNP in exon 1), C (-, the
//! reverse complement of A); FAM2 = D, E identical; F (chr6) is in no family.

use std::collections::HashSet;
use std::path::Path;
use std::process::{Command, Output};

/// Deterministic pseudo-random DNA.
fn dna(n: usize, seed: u64) -> String {
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 62) as usize] as char
        })
        .collect()
}

fn rc(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            o => o,
        })
        .collect()
}

pub struct Loc {
    chrom: &'static str,
    start: u64,
    strand: char,
    /// exons as 1-based inclusive offsets into `seq`
    exons: Vec<(u64, u64)>,
    seq: String,
    gid: &'static str,
    /// family id, or "-" for a locus in no family
    fam: &'static str,
}

impl Loc {
    fn key(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start, self.start + self.seq.len() as u64 - 1)
    }
}

fn fixture() -> Vec<Loc> {
    let s = dna(60, 7);
    let mut b = s.clone().into_bytes();
    b[4] = if b[4] == b'A' { b'C' } else { b'A' };
    let b = String::from_utf8(b).unwrap();
    let t = dna(30, 11);
    vec![
        Loc { chrom: "chr1", start: 1001, strand: '+', exons: vec![(1, 20), (41, 60)], seq: s.clone(), gid: "GA", fam: "FAM1" },
        Loc { chrom: "chr2", start: 2001, strand: '+', exons: vec![(1, 20), (41, 60)], seq: b, gid: "GB", fam: "FAM1" },
        Loc { chrom: "chr3", start: 3001, strand: '-', exons: vec![(1, 20), (41, 60)], seq: rc(&s), gid: "GC", fam: "FAM1" },
        Loc { chrom: "chr4", start: 4001, strand: '+', exons: vec![(1, 30)], seq: t.clone(), gid: "GD", fam: "FAM2" },
        Loc { chrom: "chr5", start: 5001, strand: '+', exons: vec![(1, 30)], seq: t, gid: "GE", fam: "FAM2" },
        Loc { chrom: "chr6", start: 6001, strand: '+', exons: vec![(1, 25)], seq: dna(25, 13), gid: "GF", fam: "-" },
    ]
}

/// Writes PREFIX.{clusters.tsv, loci.gff3, loci.fa, loci.paf} into `dir`; returns PREFIX.
fn write_inputs(dir: &Path, locs: &[Loc], paf: &[String]) -> String {
    let p = dir.join("t.fam").to_string_lossy().to_string();
    let mut cl = String::from("cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n");
    let mut g = String::from("##gff-version 3\n");
    let mut fa = String::new();
    for l in locs {
        let end = l.start + l.seq.len() as u64 - 1;
        if l.fam != "-" {
            cl += &format!("{}\t0\t1\t0\tNA\t{}\t{}\t{end}\n", l.fam, l.chrom, l.start);
        }
        g += &format!("{}\t.\tgene\t{}\t{end}\t.\t{}\t.\tID=gene-{};Name={}\n", l.chrom, l.start, l.strand, l.gid, l.gid);
        for &(a, b) in &l.exons {
            g += &format!(
                "{}\t.\texon\t{}\t{}\t.\t{}\t.\tParent=gene-{};gene={}\n",
                l.chrom,
                l.start + a - 1,
                l.start + b - 1,
                l.strand,
                l.gid,
                l.gid
            );
        }
        fa += &format!(">{}\n{}\n", l.key(), l.seq.to_lowercase());
    }
    std::fs::write(format!("{p}.clusters.tsv"), cl).unwrap();
    std::fs::write(format!("{p}.loci.gff3"), g).unwrap();
    std::fs::write(format!("{p}.loci.fa"), fa).unwrap();
    std::fs::write(format!("{p}.loci.paf"), paf.join("")).unwrap();
    p
}

/// One full-length PAF line between two loci of equal length.
fn paf_full(q: &Loc, t: &Loc, strand: char) -> String {
    let n = q.seq.len();
    format!("{}\t{n}\t0\t{n}\t{strand}\t{}\t{n}\t0\t{n}\t{n}\t{n}\t60\tNM:i:0\tcg:Z:{n}M\n", q.key(), t.key())
}

/// All-vs-all records within FAM1 (C against A and B is reverse-strand) and FAM2, in both directions.
fn all_pafs(l: &[Loc]) -> Vec<String> {
    vec![
        paf_full(&l[0], &l[1], '+'),
        paf_full(&l[1], &l[0], '+'),
        paf_full(&l[0], &l[2], '-'),
        paf_full(&l[2], &l[0], '-'),
        paf_full(&l[1], &l[2], '-'),
        paf_full(&l[2], &l[1], '-'),
        paf_full(&l[3], &l[4], '+'),
        paf_full(&l[4], &l[3], '+'),
    ]
}

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_family_gfa")).args(args).output().unwrap()
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

#[test]
fn loads_a_valid_run() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let out = d.path().join("out");
    let o = run(&["--prefix", &p, "--outdir", out.to_str().unwrap()]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(stderr(&o).contains("2 families selected, 5 loci, 8 alignment records"), "{}", stderr(&o));
}

#[test]
fn a_member_missing_from_the_gff3_is_an_error_naming_it() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let mut cl = std::fs::read_to_string(format!("{p}.clusters.tsv")).unwrap();
    cl += "FAM3\t0\t1\t0\tNA\tchr9\t1\t10\n";
    std::fs::write(format!("{p}.clusters.tsv"), cl).unwrap();
    let o = run(&["--prefix", &p, "--outdir", d.path().join("o").to_str().unwrap()]);
    assert!(!o.status.success());
    assert!(stderr(&o).contains("chr9:1-10") && stderr(&o).contains("FAM3"), "{}", stderr(&o));
}

#[test]
fn an_unknown_family_is_an_error_naming_it() {
    let d = tempfile::tempdir().unwrap();
    let l = fixture();
    let p = write_inputs(d.path(), &l, &all_pafs(&l));
    let o = run(&["--prefix", &p, "--outdir", d.path().join("o").to_str().unwrap(), "--families", "FAM1,NOPE"]);
    assert!(!o.status.success());
    assert!(stderr(&o).contains("family NOPE is not in"), "{}", stderr(&o));
}

#[allow(dead_code)]
fn unused_helpers_kept_for_later_tasks(_: HashSet<String>) {}
