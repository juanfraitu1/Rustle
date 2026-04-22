//! Gene-level abundance summary output (compatibility scaffold for geneabundance mode).

use crate::types::DetHashMap as HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::bpcov::BpcovStranded;
use crate::path_extract::Transcript;

/// Strand-weighted bpcov contribution from a genomic window.
///
/// Port of switch(s) block.
/// `localcov[0]` = rev strand coverage, `localcov[1]` = total coverage, `localcov[2]` = fwd strand coverage.
/// - '+' strand (s=2): fwd + unstranded × fwd/(fwd+rev) — if rev > 0, else total
/// - '-' strand (s=0): rev + unstranded × rev/(fwd+rev) — if fwd > 0, else total
/// - '.' unstranded (s=1): total − fwd − rev
pub fn strand_weighted_cov(fwd: f64, total: f64, rev: f64, strand: char) -> f64 {
    let unstranded = total - fwd - rev;
    let sum = fwd + rev;
    match strand {
        '+' => {
            // if(localcov[0]) — condition is rev > 0
            if rev > 0.0 {
                fwd + unstranded * fwd / sum
            } else {
                total
            }
        }
        '-' => {
            // if(localcov[2]) — condition is fwd > 0
            if fwd > 0.0 {
                rev + unstranded * rev / sum
            } else {
                total
            }
        }
        _ => unstranded, // cov+=localcov[1]-localcov[2]-localcov[0]
    }
}

/// Compute strand-weighted bpcov sum over all exon bases for a transcript.
///
/// Port of per-transcript bpcov accumulation used in gene_abundance output.
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

/// Per-intron low-coverage classifier. Returns a Vec<bool> of length
/// `exons.len() - 1` where `true` indicates the intron at that index has
/// average bpcov (all-strand) below `ERROR_PERC * flanking_exon_avg` — the
/// signal StringTie's retainedintron() uses via `lowintron[n1][i-1]`.
/// Returns empty when fewer than 2 exons.
pub fn compute_transcript_intron_low(tx: &Transcript, bpcov: &BpcovStranded) -> Vec<bool> {
    if tx.exons.len() < 2 {
        return Vec::new();
    }
    // StringTie rlink.cpp:18780-18793: compares LOCAL 25bp (longintronanchor) windows
    // at the donor-side — last 25bp of upstream exon vs first 25bp of intron.
    // If that fails (introncov >= exoncov * intronfrac), also tries acceptor-side
    // right-drop: last 25bp of intron vs first 25bp of next exon.
    // intronfrac = ERROR_PERC = 0.1 for long reads.
    const LONGINTRONANCHOR: u64 = 25;
    let intronfrac: f64 = std::env::var("RUSTLE_INTRON_LOW_RATIO")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.1);
    let mut out = Vec::with_capacity(tx.exons.len() - 1);
    let avg_cov = |lo: usize, hi: usize| -> f64 {
        if hi <= lo {
            return 0.0;
        }
        let mut sum = 0.0f64;
        let hi = hi.min(bpcov.plus.cov.len());
        if hi <= lo {
            return 0.0;
        }
        for i in lo..hi {
            let fwd = bpcov.plus.cov[i];
            let rev = bpcov.minus.cov.get(i).copied().unwrap_or(0.0);
            sum += fwd + rev;
        }
        sum / (hi - lo) as f64
    };
    for i in 1..tx.exons.len() {
        let prev = tx.exons[i - 1];
        let curr = tx.exons[i];
        let intron_start = prev.1; // last exon base + 1 (exclusive)... Rustle uses half-open
        let intron_end = curr.0;   // first exon base (exclusive end of intron)
        if intron_end <= intron_start {
            out.push(false);
            continue;
        }

        // StringTie parity (rlink.cpp:19334-19361): three-tier check.
        //
        // Tier 1: intron cov < 1 → LOW. Short-circuit for essentially-empty introns.
        // Tier 2: full-intron vs full-flanking-exon avg cov. intron < 0.1 * exon → LOW.
        // Tier 3: anchor-window (25bp) donor-side or acceptor-side drop. Either → LOW.
        //
        // Rustle previously only implemented Tier 3, which missed cases where the
        // whole intron is depleted but anchor windows at the boundaries have enough
        // reads to fail the drop contrast (e.g., alt-donor residue bleed).
        // Opt-out via RUSTLE_INTRON_LOW_STRINGTIE_OFF=1.
        let stringtie_mode = std::env::var_os("RUSTLE_INTRON_LOW_STRINGTIE_OFF").is_none();
        let trace_il = std::env::var_os("RUSTLE_INTRON_LOW_TRACE").is_some();
        let intron_full_cov = avg_cov(bpcov.plus.idx(intron_start), bpcov.plus.idx(intron_end));
        if stringtie_mode && intron_full_cov < 1.0 {
            if trace_il {
                eprintln!(
                    "[IL] intron={}-{} tier=1 full_cov={:.4} RESULT=LOW",
                    intron_start, intron_end, intron_full_cov
                );
            }
            out.push(true);
            continue;
        }
        if stringtie_mode {
            let prev_len = prev.1.saturating_sub(prev.0);
            let curr_len = curr.1.saturating_sub(curr.0);
            if prev_len > 0 && curr_len > 0 {
                let prev_sum = avg_cov(bpcov.plus.idx(prev.0), bpcov.plus.idx(prev.1))
                    * prev_len as f64;
                let curr_sum = avg_cov(bpcov.plus.idx(curr.0), bpcov.plus.idx(curr.1))
                    * curr_len as f64;
                let exon_full_avg = (prev_sum + curr_sum) / (prev_len + curr_len) as f64;
                if trace_il {
                    eprintln!(
                        "[IL] intron={}-{} tier=2 intron_full={:.4} exon_full={:.4} threshold={:.4}",
                        intron_start, intron_end, intron_full_cov,
                        exon_full_avg, exon_full_avg * intronfrac
                    );
                }
                if intron_full_cov < exon_full_avg * intronfrac {
                    if trace_il {
                        eprintln!("[IL]   tier=2 RESULT=LOW");
                    }
                    out.push(true);
                    continue;
                }
            }
        }

        // Tier 3 (existing): anchor-window donor/acceptor drop.
        let anchor = LONGINTRONANCHOR.min(intron_end - intron_start);
        let exon_left_start = prev.1.saturating_sub(LONGINTRONANCHOR).max(prev.0);
        let exon_left_end = prev.1;
        let intron_left_start = intron_start;
        let intron_left_end = intron_start + anchor;
        let exon_l = avg_cov(bpcov.plus.idx(exon_left_start), bpcov.plus.idx(exon_left_end));
        let intron_l = avg_cov(bpcov.plus.idx(intron_left_start), bpcov.plus.idx(intron_left_end));
        if intron_l < exon_l * intronfrac {
            out.push(true);
            continue;
        }
        let intron_right_end = intron_end;
        let intron_right_start = intron_end.saturating_sub(anchor);
        let exon_right_start = curr.0;
        let exon_right_end = curr.0 + LONGINTRONANCHOR.min(curr.1.saturating_sub(curr.0));
        let intron_r = avg_cov(bpcov.plus.idx(intron_right_start), bpcov.plus.idx(intron_right_end));
        let exon_r = avg_cov(bpcov.plus.idx(exon_right_start), bpcov.plus.idx(exon_right_end));
        out.push(intron_r < exon_r * intronfrac);
    }
    out
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
        // uses per-base bpcov for gene-level coverage.
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
