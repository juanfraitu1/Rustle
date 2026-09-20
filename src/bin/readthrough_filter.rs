//! Secondary-dominated readthrough filter (§6n9) — Rust port of
//! `bench/readthrough_secondary_filter.py`, byte-identical on stdout and stderr.
//!
//! WHAT IT IS. A readthrough record is FLAGGED when the reads carrying its FUSION junction are mostly
//! SECONDARY alignments, i.e. the record's defining junction exists at that locus mainly as a
//! multimapping echo of a molecule whose primary placement is elsewhere.
//!
//! WHY IT IS EVIDENCE-BASED AND NOT A BLANKET RULE. Register 844 refuted the blunt "drop every
//! readthrough" node rule because it deletes `PKD1P6-NPIPP1`, a genuine NPIP member with 110 primary
//! MAPQ-60 reads over a canonical GT-AG junction (§6m1 addendum). This filter keeps that record and
//! keeps `PDXDC2P-NPIPB14P` (726 primary / 1 secondary) while flagging `PKD1P4-NPIPA8` (10 / 1,150)
//! and `PKD1P3-NPIPA1` (214 / 628).
//!
//! ⚠ SCOPE AND A TRAP (§6n8). The MEDIAN readthrough in every class is purely primary-supported, and
//! only 7 of 124 paralog-joining records are secondary-majority, so this flags a small minority BY
//! DESIGN. Do not read it as "readthroughs are artifacts": supplementary support is 0.1%, so these are
//! contiguous alignments, not split ones. It flags where the EVIDENCE IS BORROWED, nothing more.
//!
//! DEFAULT IS OFF. `--max-secondary-frac` defaults to 1.0, which flags nothing and is the explicit
//! no-op. 0.50 is the "secondary-majority" setting §6n8 reports.
//!
//! §6s4: ports away the Python's two fresh-clone blockers — the `samtools view` subprocess (now a
//! noodles indexed fetch) and `sys.path.insert("/mnt/linuxdisk/.../family_cert/dna"); import dna_cert`
//! (now `--nodes` / `--exonless`, which default to that path but can point anywhere).

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};

const DEFAULT_DNA_DIR: &str = "/mnt/linuxdisk/home/juanfraitu/family_cert/dna";

#[derive(Debug, Clone)]
struct Node {
    chrom: String,
    /// 0-based half-open exon blocks, as `nodes.tsv` stores them.
    exons: Vec<(i64, i64)>,
}

/// Pure rule. Returns `(flagged, secondary_fraction)`.
///
/// A record with NO fusion-junction support is never flagged — absence of evidence is not evidence of
/// borrowing, and its fraction is NaN. At `max_secondary_frac >= 1.0` nothing is ever flagged. The
/// comparison is STRICT, so a record exactly at the threshold survives.
fn classify(primary: u64, secondary: u64, supplementary: u64, max_secondary_frac: f64) -> (bool, f64) {
    let tot = primary + secondary + supplementary;
    if tot == 0 {
        return (false, f64::NAN);
    }
    let frac = secondary as f64 / tot as f64;
    if max_secondary_frac >= 1.0 {
        return (false, frac);
    }
    (frac > max_secondary_frac, frac)
}

/// The annotated intron spanning the cut point, as a 1-based `(first intron base, first base of next
/// exon)`. `exons` are 0-based half-open and need not be sorted.
fn fusion_junction(exons: &[(i64, i64)], cut_coord: i64) -> Option<(i64, i64)> {
    let mut ex = exons.to_vec();
    ex.sort_unstable();
    for w in ex.windows(2) {
        if w[0].1 <= cut_coord && cut_coord <= w[1].0 {
            return Some((w[0].1 + 1, w[1].0 + 1));
        }
    }
    None
}

/// Count reads whose CIGAR carries exactly `junction` as an `N` gap, split by alignment class.
/// Mirrors the Python's `samtools view <bam> <chrom>:<lo>-<hi>` walk one op at a time.
fn junction_counts(
    bam_path: &str,
    chrom: &str,
    lo: i64,
    hi: i64,
    junction: (i64, i64),
) -> Result<(u64, u64, u64)> {
    let bai_path = format!("{bam_path}.bai");
    anyhow::ensure!(
        std::path::Path::new(&bai_path).exists(),
        "no .bai index beside {bam_path}"
    );
    let mut reader = noodles_bam::io::reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let index = noodles_bam::bai::read(&bai_path)?;
    let region: noodles_core::Region = format!("{chrom}:{lo}-{hi}").parse()?;

    let (mut p, mut s, mut u) = (0u64, 0u64, 0u64);
    let query = match reader.query(&header, &index, &region) {
        Ok(q) => q,
        // samtools view on a reference the BAM does not know prints nothing and exits 0.
        Err(_) => return Ok((0, 0, 0)),
    };
    for result in query {
        let rb = RecordBuf::try_from_alignment_record(&header, &result?)?;
        let Some(start) = rb.alignment_start() else {
            continue;
        };
        let mut pos = usize::from(start) as i64; // 1-based, as SAM POS
        for op in rb.cigar().as_ref() {
            use noodles_sam::alignment::record::cigar::op::Kind;
            let n = op.len() as i64;
            match op.kind() {
                Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    pos += n;
                }
                Kind::Skip => {
                    if (pos, pos + n) == junction {
                        let f = rb.flags();
                        if f.is_supplementary() {
                            u += 1;
                        } else if f.is_secondary() {
                            s += 1;
                        } else {
                            p += 1;
                        }
                    }
                    pos += n;
                }
                _ => {}
            }
        }
    }
    Ok((p, s, u))
}

fn read_tsv(path: &str) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    let fh = std::fs::File::open(path).with_context(|| format!("open {path}"))?;
    let mut lines = BufReader::new(fh).lines();
    let hdr: Vec<String> = lines
        .next()
        .transpose()?
        .unwrap_or_default()
        .split('\t')
        .map(str::to_string)
        .collect();
    let mut rows = Vec::new();
    for line in lines {
        let line = line?;
        rows.push(line.split('\t').map(str::to_string).collect());
    }
    Ok((hdr, rows))
}

fn col<'a>(hdr: &[String], row: &'a [String], name: &str) -> Option<&'a str> {
    hdr.iter()
        .position(|h| h == name)
        .and_then(|i| row.get(i))
        .map(String::as_str)
}

fn parse_blocks(s: &str) -> Vec<(i64, i64)> {
    s.split(',')
        .filter_map(|b| {
            let (a, z) = b.split_once('-')?;
            Some((a.parse().ok()?, z.parse().ok()?))
        })
        .collect()
}

/// `nodes.tsv` (exon blocks) plus `exonless.tsv` (span used as the single block).
fn load_nodes(nodes_tsv: &str, exonless_tsv: &str) -> Result<HashMap<String, Node>> {
    let mut nodes = HashMap::new();
    let (hdr, rows) = read_tsv(nodes_tsv)?;
    for r in &rows {
        let (Some(id), Some(chrom), Some(exons)) = (
            col(&hdr, r, "node_id"),
            col(&hdr, r, "chrom"),
            col(&hdr, r, "exons"),
        ) else {
            continue;
        };
        nodes.insert(
            id.to_string(),
            Node { chrom: chrom.to_string(), exons: parse_blocks(exons) },
        );
    }
    if std::path::Path::new(exonless_tsv).exists() {
        let (hdr, rows) = read_tsv(exonless_tsv)?;
        for r in &rows {
            let (Some(id), Some(chrom), Some(s0), Some(e0)) = (
                col(&hdr, r, "node_id"),
                col(&hdr, r, "chrom"),
                col(&hdr, r, "span_start0"),
                col(&hdr, r, "span_end"),
            ) else {
                continue;
            };
            let (Ok(s), Ok(e)) = (s0.parse::<i64>(), e0.parse::<i64>()) else {
                continue;
            };
            nodes.insert(
                id.to_string(),
                Node { chrom: chrom.to_string(), exons: vec![(s, e)] },
            );
        }
    }
    Ok(nodes)
}

/// Python prints a float with `repr`, so 1.0 renders as "1.0", not Rust's "1".
fn py_float(v: f64) -> String {
    if v == v.trunc() && v.is_finite() {
        format!("{v:.1}")
    } else {
        format!("{v}")
    }
}

fn usage() -> ! {
    eprintln!(
        "usage: readthrough_filter --bam BAM --cuts cuts.tsv [--max-secondary-frac F] [--out T]\n\
                                   [--nodes nodes.tsv] [--exonless exonless.tsv] | --self-test"
    );
    std::process::exit(2)
}

fn main() -> Result<()> {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    if argv.iter().any(|a| a == "--self-test") {
        self_test();
        println!("self-test OK");
        return Ok(());
    }
    let mut bam = String::new();
    let mut cuts = String::new();
    let mut out = "-".to_string();
    let mut max_secondary_frac = 1.0f64;
    let mut nodes_tsv = format!("{DEFAULT_DNA_DIR}/nodes.tsv");
    let mut exonless_tsv = format!("{DEFAULT_DNA_DIR}/exonless.tsv");
    let mut i = 0;
    while i < argv.len() {
        let need = |i: usize| -> String {
            argv.get(i + 1).cloned().unwrap_or_else(|| usage())
        };
        match argv[i].as_str() {
            "--bam" => bam = need(i),
            "--cuts" => cuts = need(i),
            "--out" => out = need(i),
            "--nodes" => nodes_tsv = need(i),
            "--exonless" => exonless_tsv = need(i),
            "--max-secondary-frac" => max_secondary_frac = need(i).parse()?,
            "--gff" => {} // accepted and ignored, as in the Python
            _ => usage(),
        }
        i += 2;
    }
    if bam.is_empty() || cuts.is_empty() {
        usage();
    }

    let nodes = load_nodes(&nodes_tsv, &exonless_tsv)?;
    let (hdr, rows) = read_tsv(&cuts)?;
    let mut recs: Vec<(String, String, u64, u64, u64, f64, bool)> = Vec::new();
    for r in &rows {
        if col(&hdr, r, "status") != Some("cut") {
            continue;
        }
        let (Some(node_id), Some(name), Some(cut)) = (
            col(&hdr, r, "node"),
            col(&hdr, r, "name"),
            col(&hdr, r, "cut_coord"),
        ) else {
            continue;
        };
        let Some(n) = nodes.get(node_id) else { continue };
        let mut ex = n.exons.clone();
        ex.sort_unstable();
        let Some(j) = fusion_junction(&ex, cut.parse()?) else {
            continue;
        };
        let (lo, hi) = (ex[0].0, ex[ex.len() - 1].1);
        let (p, s, u) = junction_counts(&bam, &n.chrom, lo, hi, j)?;
        let (flagged, frac) = classify(p, s, u, max_secondary_frac);
        recs.push((name.to_string(), n.chrom.clone(), p, s, u, frac, flagged));
    }

    // Python: sorted(rows, key=lambda x: -x[3]) — descending secondary, stable on ties.
    recs.sort_by(|a, b| b.3.cmp(&a.3));

    let mut w: Box<dyn Write> = if out == "-" {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::fs::File::create(&out)?)
    };
    writeln!(w, "name\tchrom\tprimary\tsecondary\tsupplementary\tsecondary_frac\tflagged")?;
    for (nm, ch, p, s, u, fr, fl) in &recs {
        let frac = if fr.is_nan() { String::new() } else { format!("{fr:.4}") };
        writeln!(w, "{nm}\t{ch}\t{p}\t{s}\t{u}\t{frac}\t{}", u8::from(*fl))?;
    }
    w.flush()?;
    let n_fl = recs.iter().filter(|r| r.6).count();
    eprintln!(
        "[readthrough-filter] {} records with fusion-junction support; flagged {} at max_secondary_frac={}",
        recs.len(),
        n_fl,
        py_float(max_secondary_frac)
    );
    Ok(())
}

/// The Python's `_self_test()`, case for case.
fn self_test() {
    assert!(classify(10, 1150, 0, 0.50).0);
    assert!(!classify(726, 1, 0, 0.50).0);
    assert!(
        !classify(110, 0, 0, 0.50).0,
        "PKD1P6-NPIPP1 must survive (register 844)"
    );
    assert!(!classify(0, 0, 0, 0.50).0, "no support is not evidence");
    assert!(!classify(1, 999, 0, 1.0).0, "1.0 must be the no-op");
    // exactly at the threshold is NOT flagged (strict >)
    assert!(!classify(50, 50, 0, 0.50).0);
    assert_eq!(fusion_junction(&[(10, 20), (40, 50)], 30), Some((21, 41)));
    assert_eq!(fusion_junction(&[(10, 20), (40, 50)], 5), None);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_self_test_cases_hold() {
        self_test();
    }

    #[test]
    fn no_support_is_nan_and_unflagged() {
        let (flagged, frac) = classify(0, 0, 0, 0.5);
        assert!(!flagged);
        assert!(frac.is_nan(), "absence of evidence must not render as 0.0000");
    }

    #[test]
    fn supplementary_counts_toward_the_denominator() {
        // 1 primary, 1 secondary, 2 supplementary -> 1/4, not 1/2.
        let (_, frac) = classify(1, 1, 2, 0.5);
        assert!((frac - 0.25).abs() < 1e-12);
    }

    #[test]
    fn threshold_is_strict() {
        assert!(!classify(50, 50, 0, 0.50).0, "== threshold survives");
        assert!(classify(49, 51, 0, 0.50).0, "just over is flagged");
    }

    #[test]
    fn fusion_junction_sorts_its_input() {
        assert_eq!(fusion_junction(&[(40, 50), (10, 20)], 30), Some((21, 41)));
    }

    #[test]
    fn fusion_junction_accepts_the_closed_interval_ends() {
        assert_eq!(fusion_junction(&[(10, 20), (40, 50)], 20), Some((21, 41)));
        assert_eq!(fusion_junction(&[(10, 20), (40, 50)], 40), Some((21, 41)));
    }

    #[test]
    fn py_float_matches_python_repr_for_the_cli_values() {
        assert_eq!(py_float(1.0), "1.0");
        assert_eq!(py_float(0.5), "0.5");
        assert_eq!(py_float(0.25), "0.25");
    }

    #[test]
    fn parse_blocks_reads_the_nodes_tsv_exon_format() {
        assert_eq!(parse_blocks("10-20,40-50"), vec![(10, 20), (40, 50)]);
    }
}
