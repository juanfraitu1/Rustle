//! GTF merge mode: combine transcript annotations from one or more GTF files.
//! This is a structural merge (no abundance estimation), intended to mirror the original algorithm-style
//! merge workflows where transcript models are consolidated into a non-redundant set.

use crate::types::DetHashMap as HashMap;
use anyhow::{Context, Result};
use std::path::Path;

use crate::gtf::write_gtf;
use crate::path_extract::Transcript;
use crate::transcript_filter::{dedup_equal_intron_chains, filter_min_transcript_length};
use crate::types::RunConfig;

#[derive(Debug, Clone)]
struct GtfTx {
    chrom: String,
    strand: char,
    source: String,
    exons: Vec<(u64, u64)>,
}

fn parse_gtf_attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        if !part.starts_with(key) {
            continue;
        }
        let v = part
            .strip_prefix(key)?
            .trim()
            .trim_start_matches(' ')
            .strip_prefix('"')?;
        let end = v.find('"')?;
        return Some(v[..end].to_string());
    }
    None
}

fn parse_gtf_transcripts<P: AsRef<Path>>(path: P) -> Result<Vec<GtfTx>> {
    let p = path.as_ref();
    let s = std::fs::read_to_string(p)
        .with_context(|| format!("failed to read GTF: {}", p.display()))?;

    let mut by_tid: HashMap<String, GtfTx> = Default::default();
    for line in s.lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.splitn(9, '\t');
        let chrom = fields.next().unwrap_or("").to_string();
        let source = fields.next().unwrap_or("").to_string();
        let feature = fields.next().unwrap_or("");
        let start_1b: u64 = fields.next().and_then(|x| x.parse().ok()).unwrap_or(0);
        let end_1b: u64 = fields.next().and_then(|x| x.parse().ok()).unwrap_or(0);
        let _score = fields.next();
        let strand_c = fields.next().unwrap_or(".").chars().next().unwrap_or('.');
        let _frame = fields.next();
        let attrs = fields.next().unwrap_or("");

        if feature != "exon" || chrom.is_empty() || start_1b == 0 || end_1b == 0 {
            continue;
        }
        let tid = parse_gtf_attr(attrs, "transcript_id")
            .unwrap_or_else(|| format!("{}:{}-{}", chrom, start_1b, end_1b));
        let strand = if strand_c == '-' { '-' } else { '+' };
        let exon = (start_1b.saturating_sub(1), end_1b);
        let entry = by_tid.entry(tid).or_insert_with(|| GtfTx {
            chrom: chrom.clone(),
            strand,
            source: source.clone(),
            exons: Vec::new(),
        });
        entry.exons.push(exon);
    }

    let mut out = Vec::with_capacity(by_tid.len());
    for (_tid, mut tx) in by_tid {
        tx.exons.sort_by_key(|e| e.0);
        tx.exons.dedup();
        if tx.exons.is_empty() {
            continue;
        }
        out.push(tx);
    }
    out.sort_by(|a, b| {
        (
            a.chrom.as_str(),
            a.strand,
            a.exons.first().map(|e| e.0).unwrap_or(0),
            a.exons.last().map(|e| e.1).unwrap_or(0),
        )
            .cmp(&(
                b.chrom.as_str(),
                b.strand,
                b.exons.first().map(|e| e.0).unwrap_or(0),
                b.exons.last().map(|e| e.1).unwrap_or(0),
            ))
    });
    Ok(out)
}

fn dedup_exact_structure(transcripts: Vec<Transcript>) -> Vec<Transcript> {
    let mut seen: HashMap<(String, char, Vec<(u64, u64)>), usize> = Default::default();
    let mut out: Vec<Transcript> = Vec::new();
    for tx in transcripts {
        let key = (tx.chrom.clone(), tx.strand, tx.exons.clone());
        if seen.contains_key(&key) {
            continue;
        }
        let idx = out.len();
        out.push(tx);
        seen.insert(key, idx);
    }
    out
}

fn intron_chain(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if exons.len() < 2 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(exons.len() - 1);
    for i in 0..exons.len() - 1 {
        out.push((exons[i].1, exons[i + 1].0));
    }
    out
}

fn is_chain_subset(a: &[(u64, u64)], b: &[(u64, u64)]) -> bool {
    // Strictly ordered subsequence check.
    if a.is_empty() || a.len() > b.len() {
        return false;
    }
    let mut i = 0usize;
    let mut j = 0usize;
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            i += 1;
        }
        j += 1;
    }
    i == a.len()
}

fn chains_near_equal(a: &[(u64, u64)], b: &[(u64, u64)], tol: u64) -> bool {
    if a.len() != b.len() || a.is_empty() {
        return false;
    }
    a.iter()
        .zip(b.iter())
        .all(|(&(ad, aa), &(bd, ba))| ad.abs_diff(bd) <= tol && aa.abs_diff(ba) <= tol)
}

fn exon_span(tx: &Transcript) -> Option<(u64, u64)> {
    Some((tx.exons.first()?.0, tx.exons.last()?.1))
}

fn tx_contains(a: &Transcript, b: &Transcript) -> bool {
    let (a0, a1) = match exon_span(a) {
        Some(v) => v,
        None => return false,
    };
    let (b0, b1) = match exon_span(b) {
        Some(v) => v,
        None => return false,
    };
    a0 <= b0 && b1 <= a1
}

fn merge_contained_transcripts(mut txs: Vec<Transcript>) -> Vec<Transcript> {
    if txs.len() < 2 {
        return txs;
    }
    txs.sort_by(|a, b| {
        (
            a.chrom.as_str(),
            a.strand,
            a.exons.first().map(|e| e.0).unwrap_or(0),
            a.exons.last().map(|e| e.1).unwrap_or(0),
            a.exons.len(),
        )
            .cmp(&(
                b.chrom.as_str(),
                b.strand,
                b.exons.first().map(|e| e.0).unwrap_or(0),
                b.exons.last().map(|e| e.1).unwrap_or(0),
                b.exons.len(),
            ))
    });

    let n = txs.len();
    let chains: Vec<Vec<(u64, u64)>> = txs.iter().map(|t| intron_chain(&t.exons)).collect();
    let mut remove = vec![false; n];
    for i in 0..n {
        if remove[i] {
            continue;
        }
        for j in 0..n {
            if i == j || remove[i] {
                continue;
            }
            let a = &txs[i];
            let b = &txs[j];
            if a.chrom != b.chrom || a.strand != b.strand {
                continue;
            }
            // Single-exon containment.
            if a.exons.len() == 1 && b.exons.len() == 1 {
                if tx_contains(b, a)
                    && (b.exons[0].1 - b.exons[0].0) >= (a.exons[0].1 - a.exons[0].0)
                {
                    remove[i] = true;
                }
                continue;
            }
            // Prefer structurally richer transcript when intron chain is subset and bounds are contained.
            if chains[i].len() < chains[j].len()
                && is_chain_subset(&chains[i], &chains[j])
                && tx_contains(b, a)
            {
                remove[i] = true;
            }
        }
    }

    txs.into_iter()
        .enumerate()
        .filter(|(i, _)| !remove[*i])
        .map(|(_, t)| t)
        .collect()
}

fn collapse_near_equal_intron_chains(txs: Vec<Transcript>, tol: u64) -> Vec<Transcript> {
    if txs.len() < 2 {
        return txs;
    }
    let n = txs.len();
    let chains: Vec<Vec<(u64, u64)>> = txs.iter().map(|t| intron_chain(&t.exons)).collect();
    let mut remove = vec![false; n];
    for i in 0..n {
        if remove[i] || chains[i].is_empty() {
            continue;
        }
        for j in (i + 1)..n {
            if remove[j] || txs[i].chrom != txs[j].chrom || txs[i].strand != txs[j].strand {
                continue;
            }
            if !chains_near_equal(&chains[i], &chains[j], tol) {
                continue;
            }
            let len_i: u64 = txs[i].exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let len_j: u64 = txs[j].exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            // Keep transcript with tighter boundaries around the intron chain.
            let keep_i = if len_i != len_j {
                len_i > len_j
            } else if let (Some((i0, i1)), Some((j0, j1))) =
                (exon_span(&txs[i]), exon_span(&txs[j]))
            {
                (i0 <= j0 && i1 >= j1) || (i0 == j0 && i1 == j1)
            } else {
                true
            };
            if keep_i {
                remove[j] = true;
            } else {
                remove[i] = true;
                break;
            }
        }
    }
    txs.into_iter()
        .enumerate()
        .filter(|(i, _)| !remove[*i])
        .map(|(_, t)| t)
        .collect()
}

pub fn run_merge<P: AsRef<Path>>(
    input_gtfs: &[String],
    output_gtf: P,
    config: &RunConfig,
) -> Result<()> {
    if input_gtfs.is_empty() {
        anyhow::bail!("merge mode requires at least one input GTF");
    }

    let mut all: Vec<Transcript> = Vec::new();
    for path in input_gtfs {
        let txs = parse_gtf_transcripts(path)?;
        all.extend(txs.into_iter().map(|t| Transcript {
            chrom: t.chrom,
            strand: t.strand,
            exons: t.exons.clone(),
            coverage: 1.0,
            exon_cov: vec![1.0; t.exons.len()],
            tpm: 0.0,
            fpkm: 0.0,
            source: Some(t.source),
            is_longread: false,
            longcov: 1.0,
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: false,
            hardend: false,
                    alt_tts_end: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, intron_low: Vec::new(), synthetic: false,
        }));
    }

    all = dedup_exact_structure(all);
    all = dedup_equal_intron_chains(all, config.verbose);
    // bundledist = 250 in merge mode (default for merge mode)
    all = collapse_near_equal_intron_chains(all, 250);
    all = merge_contained_transcripts(all);
    all = filter_min_transcript_length(all, config.min_transcript_length, config.verbose);
    all.sort_by(|a, b| {
        (
            a.chrom.as_str(),
            a.strand,
            a.exons.first().map(|e| e.0).unwrap_or(0),
            a.exons.last().map(|e| e.1).unwrap_or(0),
            a.exons.len(),
        )
            .cmp(&(
                b.chrom.as_str(),
                b.strand,
                b.exons.first().map(|e| e.0).unwrap_or(0),
                b.exons.last().map(|e| e.1).unwrap_or(0),
                b.exons.len(),
            ))
    });

    let mut f = std::fs::File::create(output_gtf.as_ref())
        .with_context(|| format!("failed to create {}", output_gtf.as_ref().display()))?;
    write_gtf(&all, &mut f, &config.label)?;

    if config.verbose {
        eprintln!(
            "original algorithm: merge mode wrote {} merged transcripts from {} input GTF(s)",
            all.len(),
            input_gtfs.len()
        );
    }
    Ok(())
}
