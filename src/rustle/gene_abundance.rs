//! Gene-level abundance summary output (parity scaffold for C++ reference geneabundance mode).

use crate::types::DetHashMap as HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::bpcov::BpcovStranded;
use crate::path_extract::Transcript;

/// Strand-weighted bpcov contribution from a genomic window.
///
/// Port of C++ reference switch(s) block.
/// `localcov[0]` = rev strand coverage, `localcov[1]` = total coverage, `localcov[2]` = fwd strand coverage.
/// - '+' strand (s=2): fwd + unstranded × fwd/(fwd+rev) — if rev > 0, else total
/// - '-' strand (s=0): rev + unstranded × rev/(fwd+rev) — if fwd > 0, else total
/// - '.' unstranded (s=1): total − fwd − rev
pub fn strand_weighted_cov(fwd: f64, total: f64, rev: f64, strand: char) -> f64 {
    let unstranded = total - fwd - rev;
    let sum = fwd + rev;
    match strand {
        '+' => {
            // C++: if(localcov[0]) — condition is rev > 0
            if rev > 0.0 {
                fwd + unstranded * fwd / sum
            } else {
                total
            }
        }
        '-' => {
            // C++: if(localcov[2]) — condition is fwd > 0
            if fwd > 0.0 {
                rev + unstranded * rev / sum
            } else {
                total
            }
        }
        _ => unstranded, // C++: cov+=localcov[1]-localcov[2]-localcov[0]
    }
}

/// Compute strand-weighted bpcov sum over all exon bases for a transcript.
///
/// Port of C++ reference per-transcript bpcov accumulation used in gene_abundance output.
/// Sums `strand_weighted_cov(fwd, total, rev, tx.strand)` over each base position in the transcript's exons.
/// Returns 0.0 if bpcov is empty or exons are empty.
pub fn compute_transcript_bpcov(tx: &Transcript, bpcov: &BpcovStranded) -> f64 {
    let mut sum = 0.0f64;
    for &(s, e) in &tx.exons {
        let si = bpcov.plus.idx(s);
        let ei = bpcov.plus.idx(e);
        let end = ei.min(bpcov.plus.cov.len());
        if end <= si {
            continue;
        }
        for i in si..end {
            let fwd = bpcov.plus.cov[i];
            let rev = bpcov.minus.cov.get(i).copied().unwrap_or(0.0);
            // In Rust BpcovStranded, plus+minus = all stranded reads; no separate unstranded bucket.
            // total = fwd + rev; unstranded contribution = 0 for fully-stranded data.
            let total = fwd + rev;
            sum += strand_weighted_cov(fwd, total, rev, tx.strand);
        }
    }
    sum
}

#[derive(Debug, Clone)]
struct GeneAgg {
    chrom: String,
    strand: char,
    start: u64,
    end: u64,
    n_tx: usize,
    cov_sum: f64,
    tpm_sum: f64,
    fpkm_sum: f64,
}

fn tx_bounds(tx: &Transcript) -> Option<(u64, u64)> {
    Some((tx.exons.first()?.0, tx.exons.last()?.1))
}

/// Exon-overlap Union-Find gene clustering matching write_gtf assign_gene_tx_numbers.
fn assign_gene_numbers(transcripts: &[Transcript]) -> Vec<usize> {
    let n = transcripts.len();
    let mut out = vec![0usize; n];
    if n == 0 {
        return out;
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let ta = &transcripts[a];
        let tb = &transcripts[b];
        let sa = ta.exons.first().map(|e| e.0).unwrap_or(0);
        let sb = tb.exons.first().map(|e| e.0).unwrap_or(0);
        ta.chrom
            .cmp(&tb.chrom)
            .then_with(|| ta.strand.cmp(&tb.strand))
            .then_with(|| sa.cmp(&sb))
    });

    let mut parent: Vec<usize> = (0..n).collect();
    let bounds: Vec<(u64, u64)> = order
        .iter()
        .map(|&idx| tx_bounds(&transcripts[idx]).unwrap_or((0, 0)))
        .collect();

    let find = |parent: &mut Vec<usize>, mut x: usize| -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    };

    for i in 0..n {
        let (si, ei) = bounds[i];
        let ti = &transcripts[order[i]];
        for j in (0..i).rev() {
            let (sj, ej) = bounds[j];
            let tj = &transcripts[order[j]];
            if ej <= si || sj >= ei {
                continue;
            }
            if tj.chrom != ti.chrom || tj.strand != ti.strand {
                break;
            }
            // Check actual exon overlap.
            let overlap = ti
                .exons
                .iter()
                .any(|&(as_, ae)| tj.exons.iter().any(|&(bs, be)| as_ < be && bs < ae));
            if overlap {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
        }
    }

    let mut root_to_gene: HashMap<usize, usize> = Default::default();
    let mut gene_no = 0usize;
    for i in 0..n {
        let root = find(&mut parent, i);
        if !root_to_gene.contains_key(&root) {
            gene_no += 1;
            root_to_gene.insert(root, gene_no);
        }
        out[order[i]] = root_to_gene[&find(&mut parent, i)];
    }
    out
}

/// Write gene-level abundance table:
/// gene_id, chr, strand, start_1b, end_1b, num_transcripts, cov_sum, cov_mean, tpm_sum, fpkm_sum
pub fn write_gene_abundance<P: AsRef<Path>>(
    transcripts: &[Transcript],
    out_path: P,
    label: &str,
) -> std::io::Result<()> {
    let mut agg: HashMap<usize, GeneAgg> = Default::default();
    let gno = assign_gene_numbers(transcripts);
    for (i, tx) in transcripts.iter().enumerate() {
        let gid = gno.get(i).copied().unwrap_or(0);
        if gid == 0 {
            continue;
        }
        let Some((s, e)) = tx_bounds(tx) else {
            continue;
        };
        let eagg = agg.entry(gid).or_insert_with(|| GeneAgg {
            chrom: tx.chrom.clone(),
            strand: tx.strand,
            start: s,
            end: e,
            n_tx: 0,
            cov_sum: 0.0,
            tpm_sum: 0.0,
            fpkm_sum: 0.0,
        });
        eagg.start = eagg.start.min(s);
        eagg.end = eagg.end.max(e);
        eagg.n_tx += 1;
        // Use strand-weighted bpcov sum if computed; fall back to path abundance.
        // C++ reference uses per-base bpcov for gene-level coverage.
        eagg.cov_sum += if tx.bpcov_cov > 0.0 {
            tx.bpcov_cov
        } else {
            tx.coverage.max(0.0)
        };
        eagg.tpm_sum += tx.tpm.max(0.0);
        eagg.fpkm_sum += tx.fpkm.max(0.0);
    }

    let mut ids: Vec<usize> = agg.keys().copied().collect();
    ids.sort_unstable();

    let mut w = BufWriter::new(File::create(out_path)?);
    writeln!(
        w,
        "gene_id\tchr\tstrand\tstart\tend\tnum_transcripts\tcov_sum\tcov_mean\ttpm_sum\tfpkm_sum"
    )?;
    for gid in ids {
        let a = &agg[&gid];
        let cov_mean = if a.n_tx > 0 {
            a.cov_sum / a.n_tx as f64
        } else {
            0.0
        };
        writeln!(
            w,
            "{}.{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            label,
            gid,
            a.chrom,
            a.strand,
            a.start + 1,
            a.end,
            a.n_tx,
            a.cov_sum,
            cov_mean,
            a.tpm_sum,
            a.fpkm_sum
        )?;
    }
    w.flush()
}
