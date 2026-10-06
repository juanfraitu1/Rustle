//! `o2_origin_resolve` — end-to-end proof that O2 PSV-resolution overrides the aligner primary flag.
//!
//! Fixture layout (`tests/fixtures/o2_origin_resolve/`):
//!   * one contig `c1` with two two-exon gene copies A and B of the same family (`FAM1`);
//!   * copies A and B share a 600 bp spliced cDNA backbone but differ at 24 PSV positions;
//!   * 4 support reads per copy, with different 5' offsets so they survive assembler de-duplication
//!     and collapse into one isoform per copy;
//!   * one molecule `MOL_PRIMARY_WRONG` whose primary record aligns to copy A but whose sequence
//!     matches copy B, plus a secondary record at copy B with the same AS score.
//!
//! Assertions:
//!   1. `copy_assign` assigns `MOL_PRIMARY_WRONG` to copy_idx 1 (copy B) under the default AS-tied gate.
//!   2. The emitted `*.gtf` contains a transcript whose exons lie inside copy B's genomic span.
//!   3. `missing_copy_flag --assignments` attributes `MOL_PRIMARY_WRONG` to copy B: the copy B locus
//!     has more reads than the unambiguous copy B support set, while the copy A locus does NOT
//!     count the molecule.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Command;

const FIX: &str = "tests/fixtures/o2_origin_resolve";
const CHROM: &str = "c1";
const A_START: u64 = 500;
const A_END: u64 = 1200;
const B_START: u64 = 2600;
const B_END: u64 = 3300;
const N_SUPPORT: usize = 4;

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join("o2_origin_resolve")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("create scratch dir");
    d
}

fn read(out: &PathBuf, ext: &str) -> String {
    let p = out.with_extension(ext);
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {:?}: {e}", p))
}

/// Parse a TSV with a header row into rows as field vectors.
fn parse_tsv(text: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let mut lines = text.lines();
    let header = lines
        .next()
        .expect("tsv header")
        .split('\t')
        .map(|s| s.to_string())
        .collect();
    let rows: Vec<Vec<String>> = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').map(|s| s.to_string()).collect())
        .collect();
    (header, rows)
}

fn col_idx(header: &[String], name: &str) -> usize {
    header
        .iter()
        .position(|h| h == name)
        .unwrap_or_else(|| panic!("no `{name}` column in {header:?}"))
}

/// Parse `<out>.gtf` into a map transcript_id -> (chrom, start, end, exon intervals, attrs).
fn parse_gtf(text: &str) -> BTreeMap<String, (String, u64, u64, Vec<(u64, u64)>, String)> {
    let mut exons: BTreeMap<String, Vec<(u64, u64)>> = BTreeMap::new();
    let mut transcripts: BTreeMap<String, (String, u64, u64, String)> = BTreeMap::new();
    for l in text.lines().filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let chrom = f[0].to_string();
        let start: u64 = f[3].parse().expect("gtf start");
        let end: u64 = f[4].parse().expect("gtf end");
        let attr = f[8].to_string();
        let tid = attr
            .split(';')
            .find_map(|kv| {
                let kv = kv.trim();
                kv.strip_prefix("transcript_id \"")
                    .and_then(|v| v.strip_suffix('"'))
                    .map(|s| s.to_string())
            })
            .expect("transcript_id in gtf row");
        match f[2] {
            "transcript" => {
                transcripts.insert(tid, (chrom, start, end, attr));
            }
            "exon" => {
                exons.entry(tid).or_default().push((start - 1, end));
            }
            _ => {}
        }
    }
    transcripts
        .into_iter()
        .map(|(tid, (c, s, e, a))| {
            let mut ee = exons.remove(&tid).unwrap_or_default();
            ee.sort();
            (tid, (c, s - 1, e, ee, a))
        })
        .collect()
}

#[test]
fn o2_resolves_origin_in_assignment_gtf_and_missing_copy_flag() {
    let d = scratch("e2e");
    let out = d.join("o");

    // ---- 1. copy_assign ----
    let ca = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args([
            "--bam",
            &format!("{FIX}/reads.bam"),
            "--fasta",
            &format!("{FIX}/genome.fa"),
            "--families",
            &format!("{FIX}/copies.tsv"),
            "--copies-fa",
            &format!("{FIX}/copies.fa"),
            "--region",
            "c1:0-3500",
            "--gtf",
            "--out",
            out.to_str().unwrap(),
        ])
        .output()
        .expect("copy_assign spawn");
    assert!(
        ca.status.success(),
        "copy_assign failed:\n{}",
        String::from_utf8_lossy(&ca.stderr)
    );

    // 1a. assignment: MOL_PRIMARY_WRONG -> copy_idx 1 (copy B).
    let (a_hdr, a_rows) = parse_tsv(&read(&out, "assignments.tsv"));
    let i_name = col_idx(&a_hdr, "read_name");
    let i_copy = col_idx(&a_hdr, "assigned_copy");
    let i_status = col_idx(&a_hdr, "status");
    let wrong_row = a_rows
        .iter()
        .find(|r| r[i_name] == "MOL_PRIMARY_WRONG")
        .expect("MOL_PRIMARY_WRONG in assignments.tsv");
    assert_eq!(
        wrong_row[i_copy], "1",
        "O2 must resolve the molecule to copy B"
    );
    assert_eq!(
        wrong_row[i_status], "assigned",
        "MOL_PRIMARY_WRONG must be assigned, not ambiguous/tied"
    );

    // 1b. GTF: a transcript exists at copy B, not copy A.
    let gtf = parse_gtf(&read(&out, "gtf"));
    assert!(
        !gtf.is_empty(),
        "copy_assign --gtf must emit at least one transcript"
    );
    let b_txs: Vec<_> = gtf
        .values()
        .filter(|(chrom, _s, _e, exons, _attrs)| {
            chrom == CHROM
                && exons.iter().all(|(es, ee)| B_START < *ee && B_END > *es)
                && !exons.iter().any(|(es, ee)| A_START < *ee && A_END > *es)
        })
        .collect();
    assert!(
        !b_txs.is_empty(),
        "the GTF must contain a transcript whose exons lie inside copy B and not copy A"
    );

    // ---- 2. missing_copy_flag with O2 assignments ----
    let mm_i = d.join("genome.mmi");
    let mm = Command::new("minimap2")
        .args([
            "-x",
            "splice:hq",
            "-d",
            mm_i.to_str().unwrap(),
            &format!("{FIX}/genome.fa"),
        ])
        .output()
        .expect("minimap2 index spawn");
    assert!(
        mm.status.success(),
        "minimap2 index failed:\n{}",
        String::from_utf8_lossy(&mm.stderr)
    );

    // missing_copy_flag's copy-key parser needs locus labels that encode
    // family_id/copy_index.  A small BED derived from the catalog is the
    // cleanest way to give each copy a keyed locus.
    let loci = d.join("loci.bed");
    std::fs::write(
        &loci,
        format!(
            "{CHROM}\t{A_START}\t{A_END}\tFAM1_copy_0\n\
             {CHROM}\t{B_START}\t{B_END}\tFAM1_copy_1\n"
        ),
    )
    .unwrap();

    let mcf_out = d.join("mcf");
    let mcf = Command::new(env!("CARGO_BIN_EXE_missing_copy_flag"))
        .args([
            "--bam",
            &format!("{FIX}/reads.bam"),
            "--fasta",
            &format!("{FIX}/genome.fa"),
            "--loci",
            loci.to_str().unwrap(),
            "--index",
            mm_i.to_str().unwrap(),
            "--out",
            mcf_out.to_str().unwrap(),
            "--assignments",
            &format!("{}.assignments.tsv", out.to_str().unwrap()),
            "--min-reads",
            "1",
        ])
        .output()
        .expect("missing_copy_flag spawn");
    assert!(
        mcf.status.success(),
        "missing_copy_flag failed:\n{}",
        String::from_utf8_lossy(&mcf.stderr)
    );

    // 2a. copy B locus must include the O2-resolved MOL_PRIMARY_WRONG read.
    let (m_hdr, m_rows) = parse_tsv(&read(&mcf_out, "missing_copy.tsv"));
    let i_locus = col_idx(&m_hdr, "locus");
    let i_reads = col_idx(&m_hdr, "n_reads");
    let by_locus: BTreeMap<String, usize> = m_rows
        .iter()
        .map(|r| (r[i_locus].clone(), r[i_reads].parse().unwrap()))
        .collect();

    let a_count = by_locus
        .get("FAM1_copy_0")
        .copied()
        .expect("FAM1_copy_0 row in missing_copy.tsv");
    let b_count = by_locus
        .get("FAM1_copy_1")
        .copied()
        .expect("FAM1_copy_1 row in missing_copy.tsv");

    assert_eq!(
        a_count, N_SUPPORT,
        "copy A locus must contain only the unambiguous copy-A support reads; \
         MOL_PRIMARY_WRONG (assigned to copy B) must NOT be attributed here"
    );
    assert_eq!(
        b_count,
        N_SUPPORT + 1,
        "copy B locus must contain the 4 copy-B support reads PLUS MOL_PRIMARY_WRONG, \
         attributed via --assignments despite its primary alignment being at copy A"
    );

    // 2b. Sanity: without --assignments the molecule would be counted at copy A by primary overlap.
    let mcf_no = Command::new(env!("CARGO_BIN_EXE_missing_copy_flag"))
        .args([
            "--bam",
            &format!("{FIX}/reads.bam"),
            "--fasta",
            &format!("{FIX}/genome.fa"),
            "--loci",
            loci.to_str().unwrap(),
            "--index",
            mm_i.to_str().unwrap(),
            "--out",
            d.join("mcf_no").to_str().unwrap(),
            "--min-reads",
            "1",
        ])
        .output()
        .expect("missing_copy_flag (no assignments) spawn");
    assert!(
        mcf_no.status.success(),
        "missing_copy_flag without --assignments failed:\n{}",
        String::from_utf8_lossy(&mcf_no.stderr)
    );
    let (m_hdr2, m_rows2) = parse_tsv(&read(&d.join("mcf_no"), "missing_copy.tsv"));
    let i_locus2 = col_idx(&m_hdr2, "locus");
    let i_reads2 = col_idx(&m_hdr2, "n_reads");
    let by_no: BTreeMap<String, usize> = m_rows2
        .iter()
        .map(|r| (r[i_locus2].clone(), r[i_reads2].parse().unwrap()))
        .collect();
    assert_eq!(
        by_no.get("FAM1_copy_0").copied().unwrap_or(0),
        N_SUPPORT + 1,
        "without --assignments the primary-at-A molecule must be counted at copy A"
    );
    assert_eq!(
        by_no.get("FAM1_copy_1").copied().unwrap_or(0),
        N_SUPPORT,
        "without --assignments the molecule's secondary at copy B must be ignored"
    );
}
