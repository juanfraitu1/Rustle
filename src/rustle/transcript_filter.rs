//! Post-extraction transcript filters (containment, min exon length, isofrac, TPM, merge micro exons).
//!
//! Advanced filters below mirror C++ reference `filter_predictions` and related logic:
//! - Pairwise / included: C++ reference ~18855-18952 (included_pred, single-exon containment, intronic).
//! - Retained intron: C++ reference 7558-7568 `has_retained_intron`, 7742-7743.
//! - Polymerase runoff: C++ reference 19247-19276 (single-exon between multi-exon, cov < exoncov+singlethr).
//! - Short terminal exon: C++ header SMALL_EXON 35 (C++ reference 10090, 14846).

use crate::assembly_mode::LONGINTRONANCHOR;
use crate::bitset::SmallBitset;
use crate::bpcov::Bpcov;
use crate::constants::FLOW_EPSILON;
use crate::coord::{len_half_open, overlaps_half_open};
use crate::path_extract::Transcript;
use crate::reference_gtf::RefTranscript;
use crate::trace_reference::debug_ref_stage;
use crate::types::{
    CExon, CMaxIntv, CPrediction, DetHashMap as HashMap, DetHashSet as HashSet, RunConfig,
};
use std::sync::OnceLock;

const EPSILON: f64 = crate::constants::FLOW_EPSILON;

fn predcluster_trace_enabled() -> bool {
    std::env::var_os("RUSTLE_TRACE_PREDCLUSTER").is_some()
}

fn relax_guide_lowcov_lr() -> bool {
    std::env::var_os("RUSTLE_RELAX_GUIDE_LOWCOV_LR").is_some()
}

fn relax_pairwise_guide_equiv() -> bool {
    std::env::var_os("RUSTLE_RELAX_PAIRWISE_GUIDE_EQUIV").is_some()
}

/// C++ reference `print_predcluster` readthr uses `pred[n]->cov` only (e.g. rawreads ~20090; long-read
/// longunder uses the same field at ~18991–18999). Set `RUSTLE_READTHR_LONGCOV_FALLBACK` to also
/// accept `longcov >= threshold` when `coverage` is low (legacy Rustle heuristic).
fn readthr_allow_longcov_fallback() -> bool {
    std::env::var_os("RUSTLE_READTHR_LONGCOV_FALLBACK").is_some()
}

/// Machine-readable per-transcript coverage dump for diffing against C++ reference / the reference assembler reasoning.
/// Enable with `RUSTLE_BUNDLE_COV_DUMP=1` and scope with `RUSTLE_TRACE_LOCUS=start-end` (same locus
/// filter as `[BUNDLE_TRACE]`).
fn bundle_cov_dump_enabled() -> bool {
    std::env::var_os("RUSTLE_BUNDLE_COV_DUMP").is_some() && trace_locus_range().is_some()
}

fn bundle_cov_dump_stage(stage: &str, txs: &[Transcript]) {
    if !bundle_cov_dump_enabled() {
        return;
    }
    static HEADER: std::sync::Once = std::sync::Once::new();
    HEADER.call_once(|| {
        eprintln!(
            "# RUSTLE_BUNDLE_COV_TSV v1  (set RUSTLE_TRACE_LOCUS + RUSTLE_BUNDLE_COV_DUMP=1)"
        );
        eprintln!(
            "# stage\tchrom\tstart\tend\tstrand\texons\tcov\tlongcov\tbpcov_cov\tabs_tlen_x_cov\thard_s\thard_e\tis_lr\tsource"
        );
    });
    for t in txs {
        if !tx_in_trace_locus(t) {
            continue;
        }
        let start = t.exons.first().map(|e| e.0).unwrap_or(0);
        let end = t.exons.last().map(|e| e.1).unwrap_or(0);
        let score = tx_score(t);
        let src = t.source.as_deref().unwrap_or("-");
        eprintln!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.9}\t{:.9}\t{:.9}\t{:.9}\t{}\t{}\t{}\t{}",
            stage,
            t.chrom,
            start,
            end,
            t.strand,
            t.exons.len(),
            t.coverage,
            t.longcov,
            t.bpcov_cov,
            score,
            t.hardstart as u8,
            t.hardend as u8,
            t.is_longread as u8,
            src
        );
    }
}

/// Ref/trace assist transcripts that historically bypassed predcluster because they carried
/// placeholder abundance. Enumerated `junction_path` transcripts now use **node bottleneck**
/// coverage and run through the full predcluster chain for precision; set
/// `RUSTLE_PROTECT_JUNCTION_PATH` to restore the old bypass (diagnostics only).
#[inline]
fn predcluster_protected_source(src: Option<&str>) -> bool {
    if std::env::var_os("RUSTLE_PROTECT_JUNCTION_PATH").is_some() {
        return matches!(
            src,
            Some(
                "junction_path"
                    | "junction_chain"
                    | "ref_intron_chain"
                    | "ref_chain"
                    | "ref_chain_nearmiss"
                    | "ref_chain_junction_witness",
            ),
        );
    }
    matches!(
        src,
        Some(
            "junction_chain"
                | "ref_intron_chain"
                | "ref_chain"
                | "ref_chain_nearmiss"
                | "ref_chain_junction_witness",
        ),
    )
}

fn included_trace_enabled() -> bool {
    std::env::var_os("RUSTLE_TRACE_INCLUDED").is_some()
}

fn pairwise_target_trace_enabled() -> bool {
    std::env::var_os("RUSTLE_TRACE_PAIRWISE_TARGET").is_some()
}

fn parse_pairwise_target_spec(spec: &str) -> Option<(u64, u64, usize)> {
    let spec = spec.trim();
    if spec.is_empty() {
        return None;
    }
    let (span, exons) = match spec.split_once(':') {
        Some((span, exons)) => (span.trim(), exons.trim().parse().ok()?),
        None => (spec, 0),
    };
    let (start, end) = span.split_once('-')?;
    let start: u64 = start.trim().parse().ok()?;
    let end: u64 = end.trim().parse().ok()?;
    Some((start, end, exons))
}

fn pairwise_trace_targets() -> &'static Vec<(u64, u64, usize)> {
    static TARGETS: OnceLock<Vec<(u64, u64, usize)>> = OnceLock::new();
    TARGETS.get_or_init(|| {
        let mut parsed = std::env::var("RUSTLE_TRACE_PAIRWISE_TARGET")
            .ok()
            .map(|v| {
                v.split(',')
                    .filter_map(parse_pairwise_target_spec)
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        if !parsed.is_empty() {
            return parsed;
        }
        parsed = vec![
            (23_628_384, 23_783_928, 41),
            (23_628_384, 23_783_928, 42),
            (23_628_384, 23_783_371, 41),
            (23_628_384, 23_783_371, 42),
            (23_640_804, 23_783_928, 34),
            (23_660_260, 23_783_928, 24),
            (23_662_604, 23_783_371, 24),
            (23_666_572, 23_783_371, 22),
            (23_666_572, 23_783_371, 23),
        ];
        parsed
    })
}

fn trace_locus_range() -> Option<(u64, u64)> {
    static LOCUS: OnceLock<Option<(u64, u64)>> = OnceLock::new();
    *LOCUS.get_or_init(|| {
        let value = std::env::var("RUSTLE_TRACE_LOCUS").ok()?;
        let (start, end) = value.split_once('-')?;
        let start: u64 = start.trim().parse().ok()?;
        let end: u64 = end.trim().parse().ok()?;
        Some((start.min(end), start.max(end)))
    })
}

/// Check if transcript overlaps with the trace locus range.
fn tx_in_trace_locus(t: &Transcript) -> bool {
    let Some((locus_start, locus_end)) = trace_locus_range() else {
        return false;
    };
    let tx_start = t.exons.first().map(|e| e.0).unwrap_or(0);
    let tx_end = t.exons.last().map(|e| e.1).unwrap_or(0);
    tx_start <= locus_end && tx_end >= locus_start
}

/// Log detailed transcript info for bundle tracing.
fn trace_tx_detail(label: &str, t: &Transcript, extra: Option<&str>) {
    if !tx_in_trace_locus(t) {
        return;
    }
    let start = t.exons.first().map(|e| e.0).unwrap_or(0);
    let end = t.exons.last().map(|e| e.1).unwrap_or(0);
    let introns: Vec<String> = t
        .exons
        .windows(2)
        .map(|w| format!("{}-{}", w[0].1, w[1].0))
        .collect();
    let exon_lens: Vec<u64> = t.exons.iter().map(|(s, e)| e - s).collect();
    let tlen: u64 = t.exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
    let score = (tlen.max(1) as f64).abs() * t.coverage;
    let src = t.source.as_deref().unwrap_or("-");
    eprintln!(
        "[BUNDLE_TRACE] {} {}-{} exons={} cov={:.6} longcov={:.6} bpcov_cov={:.3} tlen={} score={:.3} is_lr={} src={} exon_lens=[{}] introns=[{}] {}",
        label,
        start,
        end,
        t.exons.len(),
        t.coverage,
        t.longcov,
        t.bpcov_cov,
        tlen,
        score,
        t.is_longread,
        src,
        exon_lens.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(","),
        introns.join(","),
        extra.unwrap_or("")
    );
}

fn is_rescue_protected(t: &Transcript) -> bool {
    // Guide-matched transcripts are protected from dedup/filtering.
    if t.source.as_deref().map_or(false, |s| s.starts_with("guide:")) {
        return true;
    }
    // High-confidence rescued variants are protected.
    if matches!(
        t.source.as_deref(),
        Some("terminal_alt_acceptor") | Some("micro_exon_rescue") | Some("terminal_stop_variant")
    ) {
        return true;
    }
    // Verified complete transcripts (both boundaries found in reads/graph) are also protected.
    t.hardstart && t.hardend
}

fn pred_fate_trace_enabled() -> bool {
    std::env::var_os("RUSTLE_TRACE_PRED_FATE").is_some()
}

fn parse_trace_locus() -> Option<(u64, u64)> {
    let val = std::env::var("RUSTLE_TRACE_LOCUS").ok()?;
    let (a, b) = val.split_once('-')?;
    let start = a.trim().parse::<u64>().ok()?;
    let end = b.trim().parse::<u64>().ok()?;
    Some((start.min(end), start.max(end)))
}

fn tx_overlaps_trace_locus(t: &Transcript) -> bool {
    let Some((lo, hi)) = parse_trace_locus() else {
        return true;
    };
    let start = t.exons.first().map(|e| e.0).unwrap_or(0);
    let end = t.exons.last().map(|e| e.1).unwrap_or(0);
    start <= hi && end >= lo
}

fn tx_summary(t: &Transcript) -> String {
    let start = t.exons.first().map(|e| e.0).unwrap_or(0);
    let end = t.exons.last().map(|e| e.1).unwrap_or(0);
    let guide = t
        .source
        .as_deref()
        .and_then(|s| s.strip_prefix("guide:"))
        .unwrap_or("novel");
    format!(
        "{}-{} cov={:.4} longcov={:.1} strand={} exons={} guide=({})",
        start,
        end,
        t.coverage,
        t.longcov,
        t.strand,
        t.exons.len(),
        guide
    )
}

fn emit_pred_entries(stage: &str, txs: &[Transcript]) {
    eprintln!("print_predcluster: {} count={}", stage, txs.len());
    for (i, t) in txs.iter().enumerate() {
        if !tx_overlaps_trace_locus(t) {
            continue;
        }
        eprintln!("  pred[{}]: {}", i, tx_summary(t));
    }
}

fn pairwise_trace_target(a: &Transcript, b: &Transcript) -> bool {
    let Some((as_, ae)) = tx_span(a) else {
        return false;
    };
    let Some((bs, be)) = tx_span(b) else {
        return false;
    };
    let a_key = (as_, ae, a.exons.len());
    let b_key = (bs, be, b.exons.len());
    let interesting = pairwise_trace_targets();
    let matches = |key: (u64, u64, usize)| {
        interesting
            .iter()
            .any(|&(s, e, exons)| s == key.0 && e == key.1 && (exons == 0 || exons == key.2))
    };
    matches(a_key) && matches(b_key)
}

fn predcluster_exonic_overlap_len(a: &[(u64, u64)], b: &[(u64, u64)]) -> u64 {
    let mut total = 0u64;
    for &(as_, ae) in a {
        for &(bs, be) in b {
            if overlaps_half_open(as_, ae, bs, be) {
                let start = as_.max(bs);
                let end = ae.min(be);
                total += len_half_open(start, end);
            }
        }
    }
    total
}

fn predcluster_trace_dump(stage: &str, ref_tx: &RefTranscript, txs: &[Transcript]) {
    if !predcluster_trace_enabled() {
        return;
    }
    eprintln!(
        "[TRACE_PREDCLUSTER] ref={} stage={} txs={}",
        ref_tx.id,
        stage,
        txs.len()
    );
    for (idx, tx) in txs.iter().enumerate() {
        if tx.chrom != ref_tx.chrom || (ref_tx.strand != '.' && tx.strand != ref_tx.strand) {
            continue;
        }
        let Some((start, end)) = tx_span(tx) else {
            continue;
        };
        let overlap = predcluster_exonic_overlap_len(&ref_tx.exons, &tx.exons);
        if overlap == 0 {
            continue;
        }
        let first = tx.exons.first().copied().unwrap_or((0, 0));
        let last = tx.exons.last().copied().unwrap_or((0, 0));
        let source = tx.source.as_deref().unwrap_or("-");
        eprintln!(
            "[TRACE_PREDCLUSTER]   idx={} span={}-{} cov={:.4} exons={} first={}-{} last={}-{} source={} overlap_bp={}",
            idx,
            start,
            end,
            tx.coverage,
            tx.exons.len(),
            first.0,
            first.1,
            last.0,
            last.1,
            source,
            overlap
        );
    }
}

/// Add single-exon prediction y to prediction x (C++ reference add_pred): extend first or last exon of x and merge coverage.
pub fn add_pred(x: &mut Transcript, y: &Transcript, cov: f64) {
    if y.exons.len() != 1 || x.exons.is_empty() {
        return;
    }
    let (ys, ye) = y.exons[0];
    let y_len = len_half_open(ys, ye) as f64;
    let x_len: u64 = x.exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
    let x_len_f = x_len as f64;

    if ys < x.exons[0].1 {
        let addlen: i64 = if ye < x.exons[0].0 {
            x.exons[0].0 as i64 - ye as i64 - 1
        } else {
            0
        };
        let denom = x_len_f + addlen.max(0) as f64 + y_len;
        if denom > 0.0 {
            x.coverage = (x.coverage * x_len_f + cov * y_len) / denom;
        }
        let first = &mut x.exons[0];
        if ys < first.0 {
            first.0 = ys;
        }
    } else {
        let last_idx = x.exons.len() - 1;
        let addlen: i64 = if x.exons[last_idx].1 < ys {
            ys as i64 - x.exons[last_idx].1 as i64 - 1
        } else {
            0
        };
        let denom = x_len_f + addlen.max(0) as f64 + y_len;
        if denom > 0.0 {
            x.coverage = (x.coverage * x_len_f + cov * y_len) / denom;
        }
        let last = &mut x.exons[last_idx];
        if ye > last.1 {
            last.1 = ye;
        }
    }
}

/// Compute TPM and FPKM using standard output formulas:
/// TPM_i = cov_i * 1e6 / sum(cov)
/// FPKM_i = cov_i * 1e9 / (num_fragments * mean_frag_len * transcript_len)
/// where mean_frag_len = frag_len_sum / num_fragments, so denominator is
/// equivalent to `frag_len_sum * transcript_len`.
pub fn compute_tpm_fpkm(transcripts: &mut [Transcript], num_fragments: f64, frag_len_sum: f64) {
    if transcripts.is_empty() {
        return;
    }
    let cov_sum: f64 = transcripts.iter().map(|t| t.coverage.max(0.0)).sum();
    let cov_sum = if cov_sum < EPSILON { 1.0 } else { cov_sum };
    for tx in transcripts.iter_mut() {
        let tcov = tx.coverage.max(0.0);
        let tlen: u64 = tx.exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
        tx.tpm = tcov * 1e6 / cov_sum;
        tx.fpkm = if num_fragments > EPSILON && frag_len_sum > EPSILON && tlen > 0 {
            tcov * 1e9 / (frag_len_sum * tlen as f64)
        } else {
            0.0
        };
    }
}

/// Backward-compatible wrapper.
pub fn compute_tpm(transcripts: &mut [Transcript]) {
    compute_tpm_fpkm(transcripts, 0.0, 0.0);
}

#[derive(Debug, Clone, Default)]
pub struct PairwiseKillSummary {
    pub included_kill: usize,
    pub secontained_kill: usize,
    pub intronic_kill: usize,
    pub bettercov_kill: usize,
    pub other_kill: usize,
    pub exact_reason_note: String,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct IsofracKillSummary {
    pub longunder_kill: usize,
}

fn classify_pairwise_reason(reason: &str, summary: &mut PairwiseKillSummary) {
    if reason.contains("included") || reason == "guide_lowcov" {
        summary.included_kill += 1;
    } else if reason.contains("se_contained") {
        summary.secontained_kill += 1;
    } else if reason.contains("intronic") {
        summary.intronic_kill += 1;
    } else if reason.contains("bettercov_conflict") {
        summary.bettercov_kill += 1;
    } else {
        summary.other_kill += 1;
    }
}

fn format_reason_histogram(reason_counts: &std::collections::BTreeMap<&'static str, usize>) -> String {
    reason_counts
        .iter()
        .map(|(reason, count)| format!("{reason}={count}"))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Remove transcripts that have any exon shorter than min_exon_length (bp).
pub fn filter_min_exon_length(
    transcripts: Vec<Transcript>,
    min_exon_length: u64,
    verbose: bool,
) -> Vec<Transcript> {
    if min_exon_length == 0 {
        return transcripts;
    }
    let mut kept = Vec::with_capacity(transcripts.len());
    let mut removed = 0;
    for t in transcripts {
        let any_short = t
            .exons
            .iter()
            .any(|(s, e)| len_half_open(*s, *e) < min_exon_length);
        if any_short {
            removed += 1;
        } else {
            kept.push(t);
        }
    }
    if verbose && removed > 0 {
        eprintln!(
            "    Min exon length filter: removed {} transcripts (exon < {}bp)",
            removed, min_exon_length
        );
    }
    kept
}

/// Isofrac filter: remove transcripts whose coverage is < isofrac fraction of max coverage
/// in the same locus (overlapping, same chrom/strand). Keep if coverage >= keep_min_reads.
///
/// Dynamic scaling for long-read mode (C++ reference printResults 19041-19049):
/// - cov > 100 (CHI_WIN): isofrac *= ERROR_PERC * DROP = 0.05  (200x more permissive)
/// - cov > 50  (CHI_THR): isofrac *= DROP             = 0.5    (2x more permissive)
/// - cov <= 50:           isofrac unchanged
///
/// Also for multi-exon: only filter if cov < 5.0 (DROP/ERROR_PERC absolute minimum; C++ reference 19043).
const CHI_WIN: f64 = 100.0;
const CHI_THR: f64 = 50.0;
const ISOFRAC_DROP: f64 = 0.5;
const ISOFRAC_ERROR_PERC: f64 = 0.1;
const ISOFRAC_ABS_MIN: f64 = 5.0; // DROP/ERROR_PERC = 0.5/0.1

pub fn isofrac_filter(
    transcripts: Vec<Transcript>,
    isofrac: f64,
    _keep_min_reads: f64,
    single_exon_stricter: bool,
    verbose: bool,
) -> Vec<Transcript> {
    if transcripts.is_empty() || isofrac <= 0.0 {
        return transcripts;
    }
    let single_exon_frac = if single_exon_stricter {
        (isofrac * 10.0).min(1.0)
    } else {
        isofrac
    };
    let mut by_strand: HashMap<(String, char), Vec<Transcript>> = Default::default();
    for tx in transcripts {
        by_strand
            .entry((tx.chrom.clone(), tx.strand))
            .or_default()
            .push(tx);
    }
    let total_in: usize = by_strand.values().map(|v| v.len()).sum();
    let mut out = Vec::new();
    for (_, mut strand_txs) in by_strand {
        strand_txs.sort_by_key(|tx| tx.exons.first().map(|e| e.0).unwrap_or(0));
        let mut loci: Vec<Vec<Transcript>> = Vec::new();
        if let Some(first) = strand_txs.first().cloned() {
            let mut current_end = first.exons.last().map(|e| e.1).unwrap_or(0);
            let mut current = vec![first];
            for tx in strand_txs.into_iter().skip(1) {
                let start = tx.exons.first().map(|e| e.0).unwrap_or(0);
                if start <= current_end {
                    current_end = current_end.max(tx.exons.last().map(|e| e.1).unwrap_or(0));
                    current.push(tx);
                } else {
                    loci.push(std::mem::take(&mut current));
                    current_end = tx.exons.last().map(|e| e.1).unwrap_or(0);
                    current = vec![tx];
                }
            }
            loci.push(current);
        }
        for locus in loci {
            let max_cov: f64 = locus.iter().map(|t| t.coverage).fold(0.0, f64::max);
            let multi_max_cov: f64 = locus
                .iter()
                .filter(|t| t.exons.len() > 1)
                .map(|t| t.coverage)
                .fold(0.0, f64::max);
            for tx in locus {
                let is_single = tx.exons.len() == 1;

                // Dynamic isofrac scaling (C++ reference 19041-19049): relax for high-coverage transcripts
                let base_frac = if is_single { single_exon_frac } else { isofrac };
                let effective_frac = if tx.coverage > CHI_WIN {
                    base_frac * ISOFRAC_ERROR_PERC * ISOFRAC_DROP
                } else if tx.coverage > CHI_THR {
                    base_frac * ISOFRAC_DROP
                } else {
                    base_frac
                };

                if is_single {
                    // Single-exon: filter if cov < effective_frac * locus_max (C++ reference 19050)
                    if tx.coverage >= effective_frac * max_cov {
                        out.push(tx);
                    }
                } else {
                    // Multi-exon long-read isofrac (C++ reference 19043-19047):
                    // Filter if: (multicov==0 AND cov < frac*usedcov AND cov < ABS_MIN)
                    //         OR (multicov>0  AND cov < frac*multicov)
                    let longunder = if multi_max_cov <= 0.0 {
                        tx.coverage < effective_frac * max_cov && tx.coverage < ISOFRAC_ABS_MIN
                    } else {
                        tx.coverage < effective_frac * multi_max_cov
                    };
                    if !longunder {
                        out.push(tx);
                    }
                }
            }
        }
    }
    if verbose && out.len() < total_in {
        eprintln!(
            "    Isofrac filter: removed {} transcripts (isofrac={})",
            total_in - out.len(),
            isofrac
        );
    }
    out
}

/// Merge consecutive exons when gap <= max_gap (micro-intron merge). In-place per transcript.
pub fn merge_micro_intron_exons(
    transcripts: Vec<Transcript>,
    max_gap: u64,
    verbose: bool,
) -> Vec<Transcript> {
    if max_gap == 0 {
        return transcripts;
    }
    let mut merged_count = 0;
    let out: Vec<Transcript> = transcripts
        .into_iter()
        .map(|mut tx| {
            if tx.exons.len() < 2 {
                return tx;
            }
            let mut sorted: Vec<(u64, u64)> = tx.exons;
            sorted.sort_by_key(|e| e.0);
            let mut new_exons = Vec::new();
            let mut i = 0;
            while i < sorted.len() {
                let (start, mut end) = sorted[i];
                let mut j = i + 1;
                while j < sorted.len() {
                    let (n_start, n_end) = sorted[j];
                    let gap = n_start.saturating_sub(end);
                    if gap > max_gap {
                        break;
                    }
                    end = end.max(n_end);
                    j += 1;
                }
                new_exons.push((start, end));
                if j - i > 1 {
                    merged_count += 1;
                }
                i = j;
            }
            tx.exons = new_exons;
            tx
        })
        .collect();
    if verbose && merged_count > 0 {
        eprintln!(
            "    Merge micro-intron exons: merged {} transcript(s) (gap <={}bp)",
            merged_count, max_gap
        );
    }
    out
}

/// True if exon (s, e) is contained in some exon of container (container exons must cover [s,e]).
fn exon_contained_in((s, e): (u64, u64), container_exons: &[(u64, u64)]) -> bool {
    container_exons.iter().any(|&(s2, e2)| s2 <= s && e <= e2)
}

/// True if transcript a is fully contained in transcript b (same chrom/strand, each exon of a in b).
fn transcript_contained_in(a: &Transcript, b: &Transcript) -> bool {
    if a.chrom != b.chrom || a.strand != b.strand {
        return false;
    }
    if a.exons.len() > b.exons.len() {
        return false;
    }
    a.exons.iter().all(|&ex| exon_contained_in(ex, &b.exons))
}

/// Remove transcripts whose total length (sum of exon lengths) is less than min_length.
pub fn filter_min_transcript_length(
    transcripts: Vec<Transcript>,
    min_length: u64,
    verbose: bool,
) -> Vec<Transcript> {
    if min_length == 0 {
        return transcripts;
    }
    let mut kept = Vec::with_capacity(transcripts.len());
    let mut removed = 0;
    for t in transcripts {
        let length: u64 = t.exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
        if length < min_length {
            removed += 1;
        } else {
            kept.push(t);
        }
    }
    if verbose && removed > 0 {
        eprintln!(
            "    Min transcript length filter: removed {} transcripts (length < {}bp)",
            removed, min_length
        );
    }
    kept
}

/// Remove transcripts that are fully contained in another with at least as much coverage.
/// Keeps the container; drops the contained. When coverage is equal, keeps first in list.
pub fn filter_contained_transcripts(
    transcripts: Vec<Transcript>,
    verbose: bool,
) -> Vec<Transcript> {
    let n = transcripts.len();
    let mut drop = SmallBitset::with_capacity(n.min(64));
    for i in 0..n {
        if drop.contains(i) {
            continue;
        }
        let a = &transcripts[i];
        for j in 0..n {
            if i == j || drop.contains(j) {
                continue;
            }
            let b = &transcripts[j];
            if transcript_contained_in(a, b) && b.coverage >= a.coverage {
                drop.insert_grow(i);
                break;
            }
        }
    }
    let removed = drop.count_ones();
    if verbose && removed > 0 {
        eprintln!(
            "    Transcript containment filter: removed {} contained transcripts",
            removed
        );
    }
    transcripts
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !drop.contains(*i))
        .map(|(_, t)| t)
        .collect()
}

// =============================================================================
// Advanced filters (post-assembly; C++ reference filter_predictions + reference script)
// =============================================================================
// C++ header: DROP=0.5, ERROR_PERC=0.1, SMALL_EXON=35

type Junction = (u64, u64);

fn exons_to_junction_chain(exons: &[(u64, u64)]) -> Vec<Junction> {
    (0..exons.len().saturating_sub(1))
        .map(|k| (exons[k].1, exons[k + 1].0))
        .collect()
}

fn tx_exonic_len(tx: &Transcript) -> u64 {
    tx.exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum()
}

fn tx_span(tx: &Transcript) -> Option<(u64, u64)> {
    let s = tx.exons.first().map(|e| e.0)?;
    let e = tx.exons.last().map(|e| e.1)?;
    Some((s, e))
}

fn tx_score(tx: &Transcript) -> f64 {
    // C++ reference: sort key for isofrac filter is abs(tlen)*cov (predordCmp), not bare cov.
    // Line 18464 uses bare cov only for CMaxIntv seeding (separate from filter sort).
    (tx_exonic_len(tx).max(1) as f64) * tx.coverage
}

fn guide_equiv_id(tx: &Transcript) -> Option<&str> {
    // C++ parity: guide-ness follows pred->t_eq, not only a textual source tag.
    tx.ref_transcript_id
        .as_deref()
        .or_else(|| tx.source.as_deref().and_then(|s| s.strip_prefix("guide:")))
}

fn is_guide_pair(tx: &Transcript) -> bool {
    guide_equiv_id(tx).is_some()
}

fn guide_id_pair(tx: &Transcript) -> Option<&str> {
    guide_equiv_id(tx)
}

fn same_intron_chain(a: &Transcript, b: &Transcript) -> bool {
    if a.strand != b.strand || a.exons.len() < 2 || b.exons.len() < 2 {
        return false;
    }
    exons_to_junction_chain(&a.exons) == exons_to_junction_chain(&b.exons)
}

fn guide_chain_equivalent_pair(a: &Transcript, b: &Transcript) -> bool {
    if is_guide_pair(a) == is_guide_pair(b) {
        return false;
    }
    same_intron_chain(a, b)
}

/// C++ reference included_pred parity (coverage-free variant):
/// checks whether the smaller transcript is structurally included in the larger one.
fn included_pred(
    txs: &[Transcript],
    n1: usize,
    n2: usize,
    longreads: bool,
    bpcov: Option<&Bpcov>,
    singlethr: f64,
    error_perc: f64,
    reference: bool,
) -> bool {
    let trace = included_trace_enabled();
    let trace_pair = |reason: &str, a: &Transcript, b: &Transcript| {
        if !trace {
            return;
        }
        let Some((as_, ae)) = tx_span(a) else {
            return;
        };
        let Some((bs, be)) = tx_span(b) else {
            return;
        };
        if let Some((lo, hi)) = trace_locus_range() {
            let a_overlaps = as_ <= hi && ae >= lo;
            let b_overlaps = bs <= hi && be >= lo;
            if !a_overlaps && !b_overlaps {
                return;
            }
        } else if ae < 23_620_000 || as_ > 23_790_000 || be < 23_620_000 || bs > 23_790_000 {
            return;
        }
        eprintln!(
            "[TRACE_INCLUDED] reason={} a={}-{} exons={} cov={:.4} b={}-{} exons={} cov={:.4}",
            reason,
            as_,
            ae,
            a.exons.len(),
            a.coverage,
            bs,
            be,
            b.exons.len(),
            b.coverage
        );
    };
    if n1 >= txs.len() || n2 >= txs.len() {
        return false;
    }
    let p1 = &txs[n1];
    let p2 = &txs[n2];
    let Some((p1s, p1e)) = tx_span(p1) else {
        trace_pair("p1_no_span", p1, p2);
        return false;
    };
    let Some((p2s, p2e)) = tx_span(p2) else {
        trace_pair("p2_no_span", p1, p2);
        return false;
    };
    if !overlaps_half_open(p1s, p1e, p2s, p2e) {
        trace_pair("no_overlap", p1, p2);
        return false;
    }

    let mut big = n1;
    let mut small = n2;
    if txs[n1].exons.len() < txs[n2].exons.len() {
        big = n2;
        small = n1;
    }
    let b = &txs[big];
    let s = &txs[small];
    if b.exons.is_empty() || s.exons.is_empty() {
        trace_pair("empty_exons", b, s);
        return false;
    }

    if reference && is_guide_pair(s) {
        if is_guide_pair(b) && guide_id_pair(s) != guide_id_pair(b) {
            trace_pair("guide_id_mismatch", b, s);
            return false;
        }
        if s.exons.len() != b.exons.len() {
            trace_pair("guide_exon_count", b, s);
            return false;
        }
    }

    let intronextension = if longreads {
        100u64 // CHI_WIN
    } else {
        2 * LONGINTRONANCHOR
    };

    let mut bex = 0usize;
    while bex < b.exons.len() {
        // C++ reference uses inclusive exon ends here:
        //   small.start > big.end  <=>  small.start >= big.end_exclusive
        //   small.end   < big.start <=> small.end_exclusive <= big.start
        // Rust stores half-open exons, so abutting intervals must not be treated as overlap.
        if s.exons[0].0 >= b.exons[bex].1 {
            bex += 1;
            continue;
        }
        if s.exons[0].1 <= b.exons[bex].0 {
            trace_pair("first_exon_no_overlap", b, s);
            return false;
        }
        let mut sex = 0usize;
        while sex < s.exons.len() && bex < b.exons.len() {
            if sex + 1 == s.exons.len() {
                if bex + 1 == b.exons.len() {
                    // C++ parity: If small starts before big and they share the same 3' end,
                    // small is a 5' extension, not an inclusion
                    if s.exons[0].0 < b.exons[0].0 && s.exons[sex].1 == b.exons[bex].1 {
                        trace_pair("5prime_extension_same_3prime", b, s);
                        return false;
                    }
                    trace_pair("success_same_tail", b, s);
                    return true;
                }
                // C++ parity: Single-exon transcript that starts before big's first exon
                // is a 5' extension, not an inclusion, even if it overlaps an internal exon
                if s.exons.len() == 1 && b.exons.len() > 1 && s.exons[0].0 < b.exons[0].0 {
                    trace_pair("single_exon_5prime_extension", b, s);
                    return false;
                }
                if sex > 0 && s.exons[sex].1 > b.exons[bex + 1].0 {
                    trace_pair("last_exon_cross_next_start", b, s);
                    return false;
                }
                if sex > 0 && s.exons[sex].1 > b.exons[bex].1.saturating_add(intronextension) {
                    if let Some(bp) = bpcov {
                        let endin = s.exons[sex].1.min(b.exons[bex + 1].0);
                        let covintron = cov_avg_all(bp, b.exons[bex].1, endin);
                        if covintron > singlethr {
                            let startex = s.exons[sex].0.max(b.exons[bex].0);
                            let covexon = cov_avg_all(bp, startex, b.exons[bex].1);
                            if covintron > error_perc * covexon {
                                let exondrop = cov_avg_all(
                                    bp,
                                    b.exons[bex].1.saturating_sub(LONGINTRONANCHOR),
                                    b.exons[bex].1,
                                );
                                let introndrop = cov_avg_all(
                                    bp,
                                    b.exons[bex].1,
                                    b.exons[bex].1.saturating_add(LONGINTRONANCHOR),
                                );
                                if introndrop > error_perc * exondrop {
                                    trace_pair("last_exon_intron_drop", b, s);
                                    return false;
                                }
                            }
                        }
                    } else {
                        trace_pair("last_exon_no_bpcov", b, s);
                        return false;
                    }
                }
                trace_pair("success_last_exon", b, s);
                return true;
            }
            if bex + 1 == b.exons.len() {
                trace_pair("small_past_big", b, s);
                return false;
            }
            // C++ parity: Special case for single-exon transcripts that could be 5' terminal exons
            // If small is single-exon and starts before big, and big has multiple exons,
            // the small transcript represents a 5' extension, not an inclusion
            if s.exons.len() == 1 && b.exons.len() > 1 && s.exons[0].0 < b.exons[0].0 {
                trace_pair("single_exon_5prime_extension", b, s);
                return false;
            }
            if s.exons[sex].1 != b.exons[bex].1 {
                trace_pair("exon_end_mismatch", b, s);
                return false;
            }
            if sex == 0
                && bex > 0
                && s.exons[sex].0.saturating_add(intronextension) < b.exons[bex].0
            {
                if let Some(bp) = bpcov {
                    let startin = s.exons[sex].0.max(b.exons[bex - 1].1);
                    let covintron = cov_avg_all(bp, startin, b.exons[bex].0);
                    if covintron > singlethr {
                        let endex = s.exons[sex].1.min(b.exons[bex].1);
                        let covexon = cov_avg_all(bp, b.exons[bex].0, endex);
                        if covintron > error_perc * covexon {
                            let introndrop = cov_avg_all(
                                bp,
                                b.exons[bex].0.saturating_sub(LONGINTRONANCHOR),
                                b.exons[bex].0,
                            );
                            let exondrop = cov_avg_all(
                                bp,
                                b.exons[bex].0,
                                b.exons[bex].0.saturating_add(LONGINTRONANCHOR),
                            );
                            if introndrop > error_perc * exondrop {
                                trace_pair("first_exon_intron_drop", b, s);
                                return false;
                            }
                        }
                    }
                } else {
                    trace_pair("first_exon_no_bpcov", b, s);
                    return false;
                }
            }
            bex += 1;
            sex += 1;
            if sex >= s.exons.len() || bex >= b.exons.len() {
                trace_pair("ran_off_after_advance", b, s);
                return false;
            }
            if s.exons[sex].0 != b.exons[bex].0 {
                trace_pair("next_exon_start_mismatch", b, s);
                return false;
            }
        }
        trace_pair("fell_out_inner_loop", b, s);
        return false;
    }
    trace_pair("fell_out_outer_loop", b, s);
    false
}

#[cfg(test)]
mod tests {
    use super::{included_pred, retainedintron_like};
    use crate::path_extract::Transcript;

    fn test_tx_cov(exons: &[(u64, u64)], coverage: f64) -> Transcript {
        Transcript {
            chrom: "chrTest".to_string(),
            strand: '+',
            exons: exons.to_vec(),
            coverage,
            exon_cov: vec![coverage; exons.len()],
            tpm: 0.0,
            fpkm: 0.0,
            source: None,
            is_longread: true,
            longcov: 1.0,
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: false,
            hardend: false,
        }
    }

    fn test_tx(exons: &[(u64, u64)]) -> Transcript {
        test_tx_cov(exons, 1.0)
    }

    #[test]
    fn included_pred_treats_abutting_single_exons_as_non_overlapping() {
        let txs = vec![test_tx(&[(120, 130)]), test_tx(&[(110, 120)])];
        assert!(!included_pred(&txs, 0, 1, true, None, 0.0, 0.1, true));
    }

    #[test]
    fn included_pred_skips_abutting_big_exon_before_alignment() {
        let txs = vec![
            test_tx(&[(100, 110), (120, 130), (140, 150)]),
            test_tx(&[(110, 130), (140, 150)]),
        ];
        assert!(included_pred(&txs, 0, 1, true, None, 0.0, 0.1, true));
    }

    #[test]
    fn retainedintron_skips_exon_ending_at_next_exon_start() {
        let txs = vec![
            test_tx_cov(&[(100, 110), (120, 130)], 10.0),
            test_tx_cov(&[(80, 90), (95, 120), (140, 150)], 1.0),
        ];
        let lowintron = vec![vec![true], vec![false, false]];
        assert!(!retainedintron_like(&txs, 0, 1, &lowintron, 0.1));
    }

    #[test]
    fn retainedintron_requires_overlap_before_prev_exon_end() {
        let txs = vec![
            test_tx_cov(&[(100, 110), (120, 130)], 10.0),
            test_tx_cov(&[(80, 90), (110, 125), (140, 150)], 1.0),
        ];
        let lowintron = vec![vec![true], vec![false, false]];
        assert!(!retainedintron_like(&txs, 0, 1, &lowintron, 0.1));
    }

}

fn cov_avg_all(bpcov: &Bpcov, start: u64, end: u64) -> f64 {
    let span = len_half_open(start, end);
    if span == 0 {
        return 0.0;
    }
    let i0 = bpcov.idx(start).min(bpcov.cov.len());
    let i1 = bpcov.idx(end).min(bpcov.cov.len());
    if i1 <= i0 {
        return 0.0;
    }
    bpcov.get_cov_range(i0, i1) / (span as f64)
}

fn cov_sum_all(bpcov: &Bpcov, start: u64, end: u64) -> f64 {
    if end <= start {
        return 0.0;
    }
    let i0 = bpcov.idx(start).min(bpcov.cov.len());
    let i1 = bpcov.idx(end).min(bpcov.cov.len());
    if i1 <= i0 {
        return 0.0;
    }
    bpcov.get_cov_range(i0, i1)
}

#[inline]
fn cpred_exonic_len(pred: &CPrediction) -> u64 {
    pred.exons
        .iter()
        .map(|(s, e)| len_half_open(*s, *e))
        .sum::<u64>()
}

#[inline]
fn cpred_span(pred: &CPrediction) -> Option<(u64, u64)> {
    let s = pred.exons.first().map(|e| e.0)?;
    let e = pred.exons.last().map(|e| e.1)?;
    Some((s, e))
}

#[inline]
fn cpred_exon_cov(pred: &CPrediction, exon_idx: usize) -> f64 {
    pred.exoncov.get(exon_idx).copied().unwrap_or(pred.cov)
}

#[inline]
fn cpred_score(pred: &CPrediction) -> f64 {
    let tlen = pred.tlen.unsigned_abs() as u64;
    let tlen = if tlen == 0 {
        cpred_exonic_len(pred).max(1)
    } else {
        tlen
    };
    (tlen as f64) * pred.cov
}

/// C++-style significance check for overlapping exon pair (`update_overlap` boundary logic).
fn significant_overlap_pair_cpred(
    p: &CPrediction,
    e: usize,
    pi: &CPrediction,
    ei: usize,
    error_perc: f64,
) -> bool {
    if p.exons.is_empty() || pi.exons.is_empty() {
        return false;
    }
    let Some((p_start, p_end)) = cpred_span(p) else {
        return false;
    };
    let Some((pi_start, pi_end)) = cpred_span(pi) else {
        return false;
    };
    let plen = cpred_exonic_len(p).max(1) as f64;
    let pilen = cpred_exonic_len(pi).max(1) as f64;
    let mut sig = true;

    if e == 0 && p.exons[0].1 >= pi_end {
        let len = len_half_open(p_start, pi_end) as f64;
        if len < error_perc * plen {
            sig = false;
        }
    }
    if sig && ei == 0 && pi.exons[0].1 >= p_end {
        let len = len_half_open(pi_start, p_end) as f64;
        if len < error_perc * pilen {
            sig = false;
        }
    }
    if sig && e + 1 == p.exons.len() && p.exons[e].0 <= pi_start {
        let len = len_half_open(pi_start, p_end) as f64;
        if len < error_perc * plen {
            sig = false;
        }
    }
    if sig && ei + 1 == pi.exons.len() && pi.exons[ei].0 <= p_start {
        let len = len_half_open(p_start, pi_end) as f64;
        if len < error_perc * pilen {
            sig = false;
        }
    }
    sig
}

/// Build C++-style interval partition over prediction exons (`CMaxIntv` nodes).
fn build_cmaxintv_segments(preds: &[CPrediction]) -> Vec<CMaxIntv> {
    let mut breakpoints: Vec<u64> = Vec::new();
    for pred in preds {
        for &(s, e) in &pred.exons {
            if e > s {
                breakpoints.push(s);
                breakpoints.push(e);
            }
        }
    }
    if breakpoints.len() < 2 {
        return Vec::new();
    }
    breakpoints.sort_unstable();
    breakpoints.dedup();

    // Build all elementary segments [bp[i], bp[i+1]) once, then add each exon
    // to the segment range it spans. This avoids O(segments * exons) scans.
    let mut out: Vec<CMaxIntv> = breakpoints
        .windows(2)
        .filter_map(|w| {
            if w[1] > w[0] {
                Some(CMaxIntv::with_node(Vec::new(), w[0], w[1], 0.0))
            } else {
                None
            }
        })
        .collect();
    if out.is_empty() {
        return out;
    }

    for (predno, pred) in preds.iter().enumerate() {
        for (exonno, &(es, ee)) in pred.exons.iter().enumerate() {
            if ee <= es {
                continue;
            }
            let si = breakpoints.partition_point(|&x| x < es);
            let ei = breakpoints.partition_point(|&x| x < ee);
            if si >= ei || si >= out.len() {
                continue;
            }
            let stop = ei.min(out.len());
            let excov = cpred_exon_cov(pred, exonno);
            for seg in out.iter_mut().take(stop).skip(si) {
                seg.node.push(CExon::new(predno, exonno, excov));
                seg.cov += excov;
            }
        }
    }

    out.retain(|s| !s.node.is_empty());
    out
}

/// Recompute per-prediction coverage from `CMaxIntv` partitions.
/// Returns number of predictions with changed overall coverage.
fn recompute_cpred_cov_from_cmaxintv(preds: &mut [CPrediction], bpcov: Option<&Bpcov>) -> usize {
    if preds.len() < 2 {
        return 0;
    }
    let maxint = build_cmaxintv_segments(preds);
    if maxint.is_empty() {
        return 0;
    }

    let mut exon_mass: Vec<Vec<f64>> = preds.iter().map(|p| vec![0.0; p.exons.len()]).collect();
    let mut exon_bases: Vec<Vec<f64>> = preds.iter().map(|p| vec![0.0; p.exons.len()]).collect();

    for cur in &maxint {
        let ilen = len_half_open(cur.start, cur.end) as f64;
        if ilen > 0.0 {
            let intv_cov = if let Some(bp) = bpcov {
                cov_avg_all(bp, cur.start, cur.end)
            } else if !cur.node.is_empty() {
                cur.cov / (cur.node.len() as f64)
            } else {
                0.0
            };
            if intv_cov.is_finite() && intv_cov >= 0.0 {
                for cx in &cur.node {
                    if cx.predno >= preds.len() || cx.exonno >= preds[cx.predno].exons.len() {
                        continue;
                    }
                    // Segments are built on exon boundaries, so membership implies full coverage.
                    let olen = ilen;
                    exon_mass[cx.predno][cx.exonno] += intv_cov * olen;
                    exon_bases[cx.predno][cx.exonno] += olen;
                }
            }
        }
    }

    let mut changed = 0usize;
    for (i, pred) in preds.iter_mut().enumerate() {
        if pred.exons.is_empty() {
            continue;
        }
        let mut new_exoncov = vec![pred.cov.max(0.0); pred.exons.len()];
        let mut total_len = 0.0f64;
        let mut total_mass = 0.0f64;
        for ei in 0..pred.exons.len() {
            let elen = len_half_open(pred.exons[ei].0, pred.exons[ei].1).max(1) as f64;
            let cov = if exon_bases[i][ei] > 0.0 {
                exon_mass[i][ei] / exon_bases[i][ei]
            } else {
                cpred_exon_cov(pred, ei)
            };
            let cov = if cov.is_finite() { cov.max(0.0) } else { 0.0 };
            new_exoncov[ei] = cov;
            total_len += elen;
            total_mass += cov * elen;
        }
        let new_cov = if total_len > 0.0 {
            total_mass / total_len
        } else {
            pred.cov.max(0.0)
        };
        if (pred.cov - new_cov).abs() > EPSILON {
            changed += 1;
        }
        pred.cov = if new_cov.is_finite() {
            new_cov.max(0.0)
        } else {
            0.0
        };
        pred.exoncov = new_exoncov;
    }
    changed
}

/// Build significant-overlap matrix using C++-style `CMaxIntv` interval nodes.
fn build_significant_overlap_adj_cpred(preds: &[CPrediction], error_perc: f64) -> Vec<Vec<usize>> {
    let n = preds.len();
    let mut adj = vec![Vec::<usize>::new(); n];
    let mut seen_pairs: HashSet<(usize, usize)> = Default::default();
    let maxint = build_cmaxintv_segments(preds);
    for cur in &maxint {
        for i in 0..cur.node.len() {
            for j in (i + 1)..cur.node.len() {
                let a = cur.node[i];
                let b = cur.node[j];
                if a.predno == b.predno || a.predno >= n || b.predno >= n {
                    continue;
                }
                let (x, y) = if a.predno < b.predno {
                    (a.predno, b.predno)
                } else {
                    (b.predno, a.predno)
                };
                if !seen_pairs.insert((x, y)) {
                    continue;
                }
                if significant_overlap_pair_cpred(
                    &preds[a.predno],
                    a.exonno,
                    &preds[b.predno],
                    b.exonno,
                    error_perc,
                ) {
                    adj[a.predno].push(b.predno);
                    adj[b.predno].push(a.predno);
                }
            }
        }
    }
    adj
}

fn overlap_matrix_to_adj(overlaps: &[crate::bitset::SmallBitset]) -> Vec<Vec<usize>> {
    let n = overlaps.len();
    let mut adj = vec![Vec::<usize>::new(); n];
    for i in 0..n {
        for j in overlaps[i].ones().filter(|&j| j > i && j < n) {
            adj[i].push(j);
            adj[j].push(i);
        }
    }
    adj
}

fn build_lowintron_flags(
    txs: &[Transcript],
    bpcov: &Bpcov,
    longreads: bool,
    singlethr: f64,
    error_perc: f64,
    drop: f64,
) -> Vec<Vec<bool>> {
    let intronfrac = if longreads {
        error_perc
    } else {
        drop - error_perc
    };
    let mut out: Vec<Vec<bool>> = Vec::with_capacity(txs.len());
    for tx in txs {
        if tx.exons.len() < 2 {
            out.push(Vec::new());
            continue;
        }
        let mut flags = vec![false; tx.exons.len() - 1];
        for j in 1..tx.exons.len() {
            let left = tx.exons[j - 1];
            let right = tx.exons[j];
            let intr_s = left.1;
            let intr_e = right.0;
            if intr_e <= intr_s {
                flags[j - 1] = true;
                continue;
            }
            let introncov = cov_avg_all(bpcov, intr_s, intr_e);
            // C++ parity: if(introncov) { ... } — zero intron coverage means
            // no reads in intron, so don't flag as "low intron" (C++ reference).
            if introncov == 0.0 {
                continue;
            }
            if introncov < 1.0 || (!longreads && introncov < singlethr) {
                flags[j - 1] = true;
                continue;
            }
            // C++ parity: length-weighted average (C++ reference)
            let left_len = len_half_open(left.0, left.1).max(1) as f64;
            let right_len = len_half_open(right.0, right.1).max(1) as f64;
            let left_sum = cov_sum_all(bpcov, left.0, left.1);
            let right_sum = cov_sum_all(bpcov, right.0, right.1);
            let exoncov = (left_sum + right_sum) / (left_len + right_len);
            if introncov < exoncov * intronfrac {
                flags[j - 1] = true;
                continue;
            }
            let w = LONGINTRONANCHOR;
            let l0 = left.1.saturating_sub(w);
            let l1 = left.1;
            let r0 = right.0;
            let r1 = right.0.saturating_add(w);
            let left_exon_drop = cov_avg_all(bpcov, l0, l1);
            let left_intron_drop = cov_avg_all(bpcov, left.1, left.1.saturating_add(w));
            if left_intron_drop < left_exon_drop * intronfrac {
                flags[j - 1] = true;
                continue;
            }
            let right_intron_drop = cov_avg_all(bpcov, right.0.saturating_sub(w), r0);
            let right_exon_drop = cov_avg_all(bpcov, r0, r1);
            if right_intron_drop < right_exon_drop * intronfrac {
                flags[j - 1] = true;
            }
        }
        out.push(flags);
    }
    out
}

fn retainedintron_like(
    txs: &[Transcript],
    n1: usize,
    n2: usize,
    lowintron: &[Vec<bool>],
    frac: f64,
) -> bool {
    if n1 >= txs.len() || n2 >= txs.len() {
        return false;
    }
    let a = &txs[n1];
    let b = &txs[n2];
    if a.exons.len() < 2 || b.exons.is_empty() {
        return false;
    }
    let mut j = 0usize;
    for i in 1..a.exons.len() {
        if j > b.exons.len().saturating_sub(1) {
            return false;
        }
        if !lowintron
            .get(n1)
            .and_then(|v| v.get(i - 1))
            .copied()
            .unwrap_or(false)
        {
            continue;
        }
        if j == b.exons.len() - 1
            && b.coverage < frac * a.coverage
            && b.exons[j].0 < a.exons[i - 1].1
        {
            return true;
        }
        // C++ reference uses inclusive exon ends:
        //   b.end < a.next.start  <=>  b.end_exclusive <= a.next.start
        //   b.start <= a.prev.end <=>  b.start < a.prev.end_exclusive
        while j < b.exons.len() && b.exons[j].1 <= a.exons[i].0 {
            j += 1;
        }
        if j == 0 && b.coverage < frac * a.coverage {
            return true;
        }
        if j < b.exons.len() && b.exons[j].0 < a.exons[i - 1].1 {
            if j > 0 && j < b.exons.len() - 1 {
                return true;
            }
            if b.coverage < frac * a.coverage {
                return true;
            }
        }
    }
    false
}

/// Approximate C++ reference update_overlap() significance for one overlapping exon pair.
fn significant_overlap_pair(
    p: &Transcript,
    e: usize,
    pi: &Transcript,
    ei: usize,
    error_perc: f64,
) -> bool {
    if p.exons.is_empty() || pi.exons.is_empty() {
        return false;
    }
    let p_start = p.exons[0].0;
    let p_end = p.exons[p.exons.len() - 1].1;
    let pi_start = pi.exons[0].0;
    let pi_end = pi.exons[pi.exons.len() - 1].1;
    let plen = tx_exonic_len(p).max(1) as f64;
    let pilen = tx_exonic_len(pi).max(1) as f64;
    let mut sig = true;

    if e == 0 && p.exons[0].1 >= pi_end {
        let len = len_half_open(p_start, pi_end) as f64;
        if len < error_perc * plen {
            sig = false;
        }
    }
    if sig && ei == 0 && pi.exons[0].1 >= p_end {
        let len = len_half_open(pi_start, p_end) as f64;
        if len < error_perc * pilen {
            sig = false;
        }
    }
    if sig && e + 1 == p.exons.len() && p.exons[e].0 <= pi_start {
        let len = len_half_open(pi_start, p_end) as f64;
        if len < error_perc * plen {
            sig = false;
        }
    }
    if sig && ei + 1 == pi.exons.len() && pi.exons[ei].0 <= p_start {
        let len = len_half_open(p_start, pi_end) as f64;
        if len < error_perc * pilen {
            sig = false;
        }
    }
    sig
}

/// Build significant-overlap matrix via active-interval sweep (C++ OvlTracker-like).
fn build_significant_overlap_matrix(
    txs: &[Transcript],
    error_perc: f64,
) -> Vec<crate::bitset::SmallBitset> {
    let n = txs.len();
    let mut ov: Vec<crate::bitset::SmallBitset> = (0..n)
        .map(|_| crate::bitset::SmallBitset::with_capacity(n.min(64)))
        .collect();
    for i in 0..n {
        ov[i].insert_grow(i);
    }
    // Sweep events over exon intervals and evaluate pairwise significance
    // only when exon intervals overlap in genomic space.
    let mut all_exons: Vec<(u64, u64, usize, usize)> = Vec::new(); // (start,end,tx,exon_idx)
    for (ti, t) in txs.iter().enumerate() {
        for (ei, &(s, e)) in t.exons.iter().enumerate() {
            if e > s {
                all_exons.push((s, e, ti, ei));
            }
        }
    }
    if all_exons.len() < 2 {
        return ov;
    }

    #[derive(Clone, Copy)]
    struct Ev {
        pos: u64,
        is_start: bool,
        idx: usize,
    }
    let mut evs: Vec<Ev> = Vec::with_capacity(all_exons.len() * 2);
    for (idx, &(s, e, _, _)) in all_exons.iter().enumerate() {
        evs.push(Ev {
            pos: s,
            is_start: true,
            idx,
        });
        evs.push(Ev {
            pos: e,
            is_start: false,
            idx,
        });
    }
    evs.sort_by(|a, b| {
        a.pos
            .cmp(&b.pos)
            // End events before start events at same coordinate to preserve half-open interval semantics.
            .then_with(|| a.is_start.cmp(&b.is_start))
    });

    let mut active: Vec<usize> = Vec::new(); // indices into all_exons
    for ev in evs {
        if !ev.is_start {
            active.retain(|&k| k != ev.idx);
            continue;
        }
        let (_s, _e, ti, ei) = all_exons[ev.idx];
        for &other in &active {
            let (_os, _oe, tj, ej) = all_exons[other];
            if ti == tj || ov[ti].contains(tj) {
                continue;
            }
            if significant_overlap_pair(&txs[ti], ei, &txs[tj], ej, error_perc) {
                ov[ti].insert_grow(tj);
                ov[tj].insert_grow(ti);
            }
        }
        active.push(ev.idx);
    }
    ov
}

/// C++ reference intronic(m, M): m is mostly intronic to M with only small boundary overlap.
fn intronic_soft(m: &Transcript, m_big: &Transcript, error_perc: f64) -> bool {
    if m.exons.is_empty() || m_big.exons.len() < 2 {
        return false;
    }
    let (m_start, m_end) = (
        m.exons.first().map(|e| e.0).unwrap_or(0),
        m.exons.last().map(|e| e.1).unwrap_or(0),
    );
    let mut i = 0usize;
    while i < m_big.exons.len() && m_start >= m_big.exons[i].0 {
        i += 1;
    }
    if i == 0 || i == m_big.exons.len() {
        return false;
    }
    let left = m_big.exons[i - 1];
    let right = m_big.exons[i];
    let m_first_end = m.exons[0].1;
    let m_last_start = m.exons[m.exons.len() - 1].0;
    if m_first_end < left.1 || m_last_start > right.0 {
        return false;
    }
    let mlen = tx_exonic_len(m).max(1) as f64;
    if m_start <= left.1 {
        let len = len_half_open(m_start, left.1) as f64;
        if len > error_perc * mlen {
            return false;
        }
    }
    if m_end >= right.0 {
        let len = len_half_open(right.0, m_end) as f64;
        if len > error_perc * mlen {
            return false;
        }
    }
    true
}

/// C++ reference transcript_overlap parity:
/// must overlap genomically, and overlap cannot be only first-vs-last terminal exon crossing.
fn transcript_overlap_sig(a: &Transcript, b: &Transcript) -> bool {
    if a.exons.is_empty() || b.exons.is_empty() {
        return false;
    }
    let (as_, ae) = (a.exons[0].0, a.exons[a.exons.len() - 1].1);
    let (bs, be) = (b.exons[0].0, b.exons[b.exons.len() - 1].1);
    if !overlaps_half_open(as_, ae, bs, be) {
        return false;
    }
    let a_last_start = a.exons[a.exons.len() - 1].0;
    let b_last_start = b.exons[b.exons.len() - 1].0;
    let a_first_end = a.exons[0].1;
    let b_first_end = b.exons[0].1;
    !(a_last_start <= b_first_end || b_last_start <= a_first_end)
}

fn junction_chain_is_subset(a: &[Junction], b: &[Junction]) -> bool {
    if a.is_empty() || a.len() >= b.len() {
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

/// Long-read isofrac filter matching C++ reference print_predcluster() lines 18961-19017.
/// For each maximal overlap interval, sorts transcripts by abs(tlen)*cov (descending),
/// then filters secondary transcripts whose coverage is below a dynamic isofrac threshold
/// relative to the accumulated coverage of higher-ranking transcripts.
/// With flow enabled, tx.coverage is flow-distributed (like C++ pred->cov).
fn isofrac(
    transcripts: Vec<Transcript>,
    isofrac: f64,
    error_perc: f64,
    drop: f64,
    _bpcov: Option<&Bpcov>,
    verbose: bool,
    keep_min_abundance: f64,
) -> Vec<Transcript> {
    isofrac_with_summary(
        transcripts,
        isofrac,
        error_perc,
        drop,
        _bpcov,
        verbose,
        keep_min_abundance,
    )
    .0
}

fn isofrac_with_summary(
    transcripts: Vec<Transcript>,
    isofrac: f64,
    error_perc: f64,
    drop: f64,
    _bpcov: Option<&Bpcov>,
    verbose: bool,
    keep_min_abundance: f64,
) -> (Vec<Transcript>, IsofracKillSummary) {
    if std::env::var_os("RUSTLE_SKIP_UNDERTHRESHOLD").is_some() {
        return (transcripts, IsofracKillSummary::default());
    }
    if transcripts.len() <= 1 {
        return (transcripts, IsofracKillSummary::default());
    }
    let txs = transcripts;
    let n = txs.len();
    let mut dead = SmallBitset::with_capacity(n.min(64));
    let mut summary = IsofracKillSummary::default();

    let mut starts: std::collections::BTreeMap<u64, Vec<usize>> = std::collections::BTreeMap::new();
    let mut ends: std::collections::BTreeMap<u64, Vec<usize>> = std::collections::BTreeMap::new();
    for (ti, t) in txs.iter().enumerate() {
        for &(s, e) in &t.exons {
            if e <= s {
                continue;
            }
            starts.entry(s).or_default().push(ti);
            ends.entry(e).or_default().push(ti);
        }
    }
    let mut breakpoints: Vec<u64> = starts.keys().chain(ends.keys()).copied().collect();
    breakpoints.sort_unstable();
    breakpoints.dedup();
    let mut active = crate::bitset::SmallBitset::with_capacity(n.min(64));
    let mut removed = 0usize;
    for w in breakpoints.windows(2) {
        let pos = w[0];
        let next = w[1];
        if let Some(v) = ends.get(&pos) {
            for &ti in v {
                active.remove(ti);
            }
        }
        if let Some(v) = starts.get(&pos) {
            for &ti in v {
                active.insert_grow(ti);
            }
        }
        if next <= pos {
            continue;
        }
        // C++ parity: skip intervals beyond bpcov bounds (C++ reference)
        if let Some(bpc) = _bpcov {
            if pos < bpc.bundle_start || (pos - bpc.bundle_start) as usize >= bpc.cov.len() {
                continue;
            }
        }
        let mut uniq: Vec<usize> = active.ones().collect();
        uniq.retain(|&k| !dead.contains(k));
        if uniq.is_empty() {
            continue;
        }

        // C++: sort by abs(tlen)*cov descending (predordCmp)
        uniq.sort_by(|&a, &b| {
            tx_score(&txs[b])
                .partial_cmp(&tx_score(&txs[a]))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let first = uniq[0];
        let first_cov = txs[first].coverage;
        let mut usedcov = [0.0f64; 2];
        let mut multicov = [0.0f64; 2];
        let fs = if txs[first].strand == '+' {
            1usize
        } else {
            0usize
        };
        usedcov[fs] = txs[first].coverage;
        if txs[first].exons.len() > 1 {
            multicov[fs] = usedcov[fs];
        }

        for &k in uniq.iter().skip(1) {
            if dead.contains(k) {
                continue;
            }
            let sidx = if txs[k].strand == '+' { 1usize } else { 0usize };
            let cov = txs[k].coverage;

            // C++: skip isofrac check for guide-matched transcripts (pred->t_eq)
            if is_guide_pair(&txs[k]) || is_rescue_protected(&txs[k]) {
                usedcov[sidx] += cov;
                if txs[k].exons.len() > 1 {
                    multicov[sidx] += cov;
                }
                continue;
            }

            // C++ parity: CHI_WIN / CHI_THR scaling uses `pred[n]->cov` only (C++ reference).
            // That is flow-derived `coverage` here, not `longcov` (pre-flow abundance).
            // `raw_cov` below is only for `--transcript-isofrac-keep-min` (abundance floor).
            let raw_cov = txs[k].longcov.max(cov);
            let mut isofraclong = isofrac;
            if cov > 100.0 {
                isofraclong = isofrac * error_perc * drop;
            } else if cov > 50.0 {
                isofraclong = isofrac * drop;
            }

            // C++: multi-exon vs single-exon threshold check (use effective_cov for comparison)
            let mut longunder = if txs[k].exons.len() > 1 {
                (multicov[sidx] <= 0.0
                    && cov < isofraclong * usedcov[sidx]
                    && cov < drop / error_perc)
                    || cov < isofraclong * multicov[sidx]
            } else {
                cov < isofraclong * usedcov[sidx]
            };
            // Optional floor: keep isoforms with at least this read-abundance (longcov/cov max)
            // even when longunder would drop them (CLI: --transcript-isofrac-keep-min).
            if longunder && keep_min_abundance > 0.0 && raw_cov >= keep_min_abundance {
                longunder = false;
            }

            if std::env::var("RUSTLE_FILTER_DEBUG").is_ok()
                || std::env::var("RUSTLE_PARITY_DEBUG").is_ok()
            {
                if trace_locus_range().is_some() && !tx_in_trace_locus(&txs[k]) {
                    continue;
                }
                let first_exon_start = txs[k].exons.first().map(|(s, _)| *s).unwrap_or(0);
                let tlen = tx_exonic_len(&txs[k]);
                let cov_per_base = cov / tlen.max(1) as f64;
                eprintln!("PARITY_LONGUNDER win={}-{} firstcov={:.4} start={} exons={} cov={:.4} cov_per_base={:.4} tlen={} usedcov={:.4} multicov={:.4} isofraclong={:.6} killed={}",
                    pos, next, first_cov,
                    first_exon_start, txs[k].exons.len(), cov, cov_per_base, tlen,
                    usedcov[sidx], multicov[sidx], isofraclong, longunder);
            }

            if longunder {
                if std::env::var_os("RUSTLE_UNDERTHRESHOLD_DIAG").is_some() {
                    if trace_locus_range().is_some() && !tx_in_trace_locus(&txs[k]) {
                        continue;
                    }
                    let first_exon_start = txs[k].exons.first().map(|e| e.0).unwrap_or(0);
                    eprintln!(
                        "UNDERTHRESHOLD_KILL win={}-{} firstcov={:.4} start={} exons={} cov={:.4} usedcov={:.4} multicov={:.4} isofraclong={:.6} strand={}",
                        pos, next, first_cov, first_exon_start, txs[k].exons.len(), cov, usedcov[sidx], multicov[sidx], isofraclong, txs[k].strand
                    );
                }
                dead.insert_grow(k);
                removed += 1;
                summary.longunder_kill += 1;
            } else {
                usedcov[sidx] += cov;
                if txs[k].exons.len() > 1 {
                    multicov[sidx] += cov;
                }
            }
        }
    }

    if verbose && removed > 0 {
        eprintln!(
            "    Longread interval threshold filter: removed {} transcript(s)",
            removed
        );
    }
    (
        txs.into_iter()
            .enumerate()
            .filter(|(i, _)| !dead.contains(*i))
            .map(|(_, t)| t)
            .collect(),
        summary,
    )
}

fn eliminate_incomplete_vs_guides(transcripts: Vec<Transcript>, verbose: bool) -> Vec<Transcript> {
    if transcripts.len() < 2 {
        return transcripts;
    }
    let txs = transcripts;
    let n = txs.len();
    let mut dead = SmallBitset::with_capacity(n.min(64));
    let mut guide_introns: HashSet<(u64, char, u64)> = Default::default();
    let mut guide_chains: Vec<(char, Vec<Junction>)> = Vec::new();
    for t in &txs {
        if !is_guide_pair(t) || t.exons.len() < 2 {
            continue;
        }
        let chain = exons_to_junction_chain(&t.exons);
        for &(jstart, jend) in &chain {
            guide_introns.insert((jstart, t.strand, jend));
        }
        if !chain.is_empty() {
            guide_chains.push((t.strand, chain));
        }
    }
    if guide_introns.is_empty() {
        return txs;
    }
    let mut removed = 0usize;
    for (i, t) in txs.iter().enumerate() {
        if is_guide_pair(t) || t.exons.len() < 2 {
            continue;
        }
        let mut all_in_guides = true;
        for w in t.exons.windows(2) {
            if !guide_introns.contains(&(w[0].1, t.strand, w[1].0)) {
                all_in_guides = false;
                break;
            }
        }
        if all_in_guides {
            let tx_chain = exons_to_junction_chain(&t.exons);
            let exact_guide_chain_match = guide_chains
                .iter()
                .any(|(strand, gchain)| *strand == t.strand && *gchain == tx_chain);
            if exact_guide_chain_match {
                continue;
            }
            let proper_subset_of_guide = guide_chains.iter().any(|(strand, gchain)| {
                *strand == t.strand && junction_chain_is_subset(&tx_chain, gchain)
            });
            if proper_subset_of_guide {
                dead.insert_grow(i);
                removed += 1;
            }
        }
    }
    if verbose && removed > 0 {
        eprintln!(
            "    Incomplete-vs-guide filter: removed {} transcript(s)",
            removed
        );
    }
    txs.into_iter()
        .enumerate()
        .filter(|(i, _)| !dead.contains(*i))
        .map(|(_, t)| t)
        .collect()
}

/// Remove multi-exon transcripts whose first or last exon is shorter than min_terminal_exon_len (bp).
/// Single-exon transcripts are kept.
/// Corresponds to C++ header SMALL_EXON 35; reference script filter_short_terminal_exons.
pub fn filter_short_terminal_exons(
    transcripts: Vec<Transcript>,
    min_terminal_exon_len: u64,
    verbose: bool,
) -> Vec<Transcript> {
    if min_terminal_exon_len == 0 {
        return transcripts;
    }
    let mut kept = Vec::with_capacity(transcripts.len());
    let mut removed = 0;
    for tx in transcripts {
        if tx.exons.len() < 2 {
            kept.push(tx);
            continue;
        }
        let first_len = tx.exons[0].1.saturating_sub(tx.exons[0].0);
        let last_len = tx.exons[tx.exons.len() - 1]
            .1
            .saturating_sub(tx.exons[tx.exons.len() - 1].0);
        if first_len < min_terminal_exon_len || last_len < min_terminal_exon_len {
            removed += 1;
        } else {
            kept.push(tx);
        }
    }
    if verbose && removed > 0 {
        eprintln!(
            "    Short terminal exon filter (<{}bp): removed {} transcripts",
            min_terminal_exon_len, removed
        );
    }
    kept
}

/// Pairwise overlap filter: junction inclusion (n2 junctions subset of n1) + coverage ratio,
/// single-exon containment, and intronic containment. Runs on transcripts already in same bundle.
/// Mirrors C++ reference filter_predictions: included_pred (17552), single-exon elimination (18908-18919),
/// intronic elimination (18936-18938). Uses isofrac, lowisofrac (C++ header), error_perc (ERROR_PERC 0.1), drop (DROP 0.5).
pub fn pairwise_overlap_filter(
    transcripts: Vec<Transcript>,
    isofrac: f64,
    lowisofrac: f64,
    error_perc: f64,
    drop: f64,
    singlethr: f64,
    longreads: bool,
    mixed_mode: bool,
    bpcov: Option<&Bpcov>,
    verbose: bool,
) -> Vec<Transcript> {
    pairwise_overlap_filter_with_summary(
        transcripts,
        isofrac,
        lowisofrac,
        error_perc,
        drop,
        singlethr,
        longreads,
        mixed_mode,
        bpcov,
        verbose,
    )
    .0
}

pub fn pairwise_overlap_filter_with_summary(
    transcripts: Vec<Transcript>,
    isofrac: f64,
    lowisofrac: f64,
    error_perc: f64,
    drop: f64,
    singlethr: f64,
    longreads: bool,
    mixed_mode: bool,
    bpcov: Option<&Bpcov>,
    verbose: bool,
) -> (Vec<Transcript>, PairwiseKillSummary) {
    if transcripts.len() <= 1 {
        return (transcripts, PairwiseKillSummary::default());
    }
    let mut txs = transcripts;
    let mut pairwise_summary = PairwiseKillSummary::default();
    let mut reason_counts: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();
    let mut cpreds: Vec<CPrediction> = txs
        .iter()
        .enumerate()
        .map(|(i, t)| t.to_cprediction(i as i32))
        .collect();
    let mut cov_updates = 0usize;
    // C++ parity: long-read path does not run CMaxIntv coverage redistribution.
    // Keep CMaxIntv overlap/coverage recompute strictly for non-long-read paths.
    let use_cmaxintv = cpreds.len() > 1 && !longreads;
    if use_cmaxintv {
        // C++ parity: run interval recomputation in CMaxIntv space before overlap elimination.
        // Default to convergence (with a safety cap); env var can tighten/expand the cap.
        let cmax_iters = std::env::var("RUSTLE_CMAXINTV_ITERS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(32)
            .max(1);
        for _ in 0..cmax_iters {
            let changed = recompute_cpred_cov_from_cmaxintv(&mut cpreds, bpcov);
            cov_updates += changed;
            if changed == 0 {
                break;
            }
        }
        for (i, p) in cpreds.iter().enumerate() {
            txs[i].coverage = p.cov;
            if p.exoncov.len() == txs[i].exons.len() {
                txs[i].exon_cov = p.exoncov.clone();
            }
        }
    }
    let cpred_scores: Vec<f64> = cpreds.iter().map(cpred_score).collect();
    let mut ord: Vec<usize> = (0..txs.len()).collect();
    if use_cmaxintv {
        ord.sort_by(|&a, &b| {
            cpred_scores[b]
                .partial_cmp(&cpred_scores[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    } else {
        ord.sort_by(|&a, &b| {
            tx_score(&txs[b])
                .partial_cmp(&tx_score(&txs[a]))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    let n = txs.len();
    let overlaps_adj = if use_cmaxintv {
        build_significant_overlap_adj_cpred(&cpreds, error_perc)
    } else {
        overlap_matrix_to_adj(&build_significant_overlap_matrix(&txs, error_perc))
    };
    let overlap_sets: Vec<HashSet<usize>> = overlaps_adj
        .iter()
        .map(|v| v.iter().copied().collect())
        .collect();
    let mut ord_pos = vec![usize::MAX; n];
    for (pos, &idx) in ord.iter().enumerate() {
        ord_pos[idx] = pos;
    }
    let mut ord_neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        if overlaps_adj[i].is_empty() {
            continue;
        }
        let mut nb = overlaps_adj[i].clone();
        nb.sort_by_key(|&j| ord_pos[j]);
        ord_neighbors[i] = nb;
    }
    let mut dead = SmallBitset::with_capacity(n.min(64));
    let mut kill_reason: Vec<&str> = vec![""; n];
    let mut kill_by: Vec<usize> = vec![usize::MAX; n];
    let mut bettercov: Vec<Vec<usize>> = vec![Vec::new(); n];
    let lowintron =
        bpcov.map(|bp| build_lowintron_flags(&txs, bp, longreads, singlethr, error_perc, drop));
    let parity_debug = std::env::var("RUSTLE_PARITY_DEBUG").is_ok();
    let predcluster_trace = std::env::var_os("RUSTLE_PARITY_PRED_FINAL").is_some() || parity_debug;

    for oi in 0..ord.len().saturating_sub(1) {
        let n1 = ord[oi];
        if dead.contains(n1) {
            continue;
        }
        for &n2 in &ord_neighbors[n1] {
            if ord_pos[n2] <= oi {
                continue;
            }
            if dead.contains(n2) {
                continue;
            }

            let t1 = &txs[n1];
            let t2 = &txs[n2];
            let Some((s1, e1)) = tx_span(t1) else {
                continue;
            };
            let Some((s2, e2)) = tx_span(t2) else {
                continue;
            };
            if pairwise_target_trace_enabled() && pairwise_trace_target(t1, t2) {
                eprintln!(
                    "[TRACE_PAIRWISE_TARGET] compare n1={} {}-{} exons={} cov={:.4} n2={} {}-{} exons={} cov={:.4}",
                    n1,
                    s1,
                    e1,
                    t1.exons.len(),
                    t1.coverage,
                    n2,
                    s2,
                    e2,
                    t2.exons.len(),
                    t2.coverage
                );
            }
            let n1g = is_guide_pair(t1);
            let n2g = is_guide_pair(t2);
            if relax_pairwise_guide_equiv() && guide_chain_equivalent_pair(t1, t2) {
                continue;
            }

            // Macro-like closure to record kill reason for parity debug
            macro_rules! kill {
                ($idx:expr, $killer:expr, $reason:expr) => {
                    dead.insert_grow($idx);
                    let vtx = &txs[$idx];
                    // Bundle tracing: log killed transcripts
                    if tx_in_trace_locus(vtx) {
                        let extra = format!("KILLED_BY={} REASON={}", $killer, $reason);
                        trace_tx_detail("PAIRWISE_KILL", vtx, Some(&extra));
                    }
                    if predcluster_trace {
                        let vstart = vtx.exons.first().map(|(s, _)| *s).unwrap_or(0);
                        let vend = vtx.exons.last().map(|(_, e)| *e).unwrap_or(0);
                        if $killer < txs.len() {
                            let ktx = &txs[$killer];
                            let kstart = ktx.exons.first().map(|(s, _)| *s).unwrap_or(0);
                            let kend = ktx.exons.last().map(|(_, e)| *e).unwrap_or(0);
                            eprintln!(
                                "--- print_predcluster: FILTER pred[{}] {}-{} cov={:.4} strand={} exons={} KILLED_BY pred[{}] {}-{} cov={:.4} reason={}",
                                $idx,
                                vstart,
                                vend,
                                vtx.coverage,
                                vtx.strand,
                                vtx.exons.len(),
                                $killer,
                                kstart,
                                kend,
                                ktx.coverage,
                                $reason.to_ascii_uppercase()
                            );
                        } else {
                            eprintln!(
                                "--- print_predcluster: FILTER pred[{}] {}-{} cov={:.4} strand={} exons={} KILLED_BY pred[{}] reason={}",
                                $idx,
                                vstart,
                                vend,
                                vtx.coverage,
                                vtx.strand,
                                vtx.exons.len(),
                                $killer,
                                $reason.to_ascii_uppercase()
                            );
                        }
                    }
                    if parity_debug {
                        kill_reason[$idx] = $reason;
                        kill_by[$idx] = $killer;
                    }
                    classify_pairwise_reason($reason, &mut pairwise_summary);
                    *reason_counts.entry($reason).or_insert(0) += 1;
                };
            }

            if !n2g {
                if n1g && t2.coverage < error_perc * t1.coverage {
                    let kill_lr_multiexon = if relax_guide_lowcov_lr()
                        && t2.is_longread
                        && t2.exons.len() > 1
                    {
                        // Optional high-sensitivity mode: in long-read runs,
                        // only apply guide_lowcov when the LR transcript is
                        // structurally included in the guide.
                        included_pred(
                            &txs, n1, n2, longreads, bpcov, singlethr, error_perc, true,
                        )
                    } else {
                        true
                    };
                    if kill_lr_multiexon {
                        kill!(n2, n1, "guide_lowcov");
                        continue;
                    }
                }

                // C++ reference — relax anchor from 25 to 10 for long-read predictions with cov > readthr
                let anchor = if t2.is_longread && t2.coverage > 1.0 {
                    10u64
                } else {
                    LONGINTRONANCHOR
                };
                let first_ok = t2
                    .exons
                    .first()
                    .map(|x| x.1.saturating_sub(x.0) >= anchor)
                    .unwrap_or(false);
                let last_ok = t2
                    .exons
                    .last()
                    .map(|x| x.1.saturating_sub(x.0) >= anchor)
                    .unwrap_or(false);
                if !(first_ok && last_ok) {
                    kill!(n2, n1, "short_terminal_exon");
                    continue;
                }
                // Retained-intron filter is primarily a short-read artifact guard.
                // In long-read mode it can incorrectly eliminate real long-read isoforms (e.g. STRG.320.*),
                // so skip it when both transcripts are long-read.
                if lowintron.is_some()
                    && !(t1.is_longread && t2.is_longread)
                    && retainedintron_like(
                        &txs,
                        n1,
                        n2,
                        lowintron.as_ref().unwrap(),
                        error_perc,
                    )
                {
                    kill!(n2, n1, "retained_intron");
                    continue;
                }
                if t1.strand != t2.strand {
                    // C++ reference — single-exon OR long-read n1 kills low-cov short-read n2
                    if t2.exons.len() == 1
                        || (t1.is_longread && !t2.is_longread && t2.coverage < 1.0 / error_perc)
                    {
                        kill!(n2, n1, "strand_conflict_se_or_lr");
                    } else if !n1g && t1.exons.len() == 1 && t1.coverage < t2.coverage + singlethr {
                        kill!(n1, n2, "strand_conflict_n1_se_lowcov");
                        break;
                    } else if t2.exons.len() <= t1.exons.len()
                        && included_pred(
                            &txs, n1, n2, longreads, bpcov, singlethr, error_perc, true,
                        )
                    {
                        kill!(n2, n1, "strand_conflict_included");
                    } else if !t2.is_longread && !mixed_mode && bpcov.is_some() {
                        let mut lowcov = true;
                        let bp = bpcov.unwrap();
                        for (eidx, &(es, ee)) in t2.exons.iter().enumerate() {
                            let len = len_half_open(es, ee).max(1) as f64;
                            let totalcov = cov_sum_all(bp, es, ee);
                            let pred_exoncov =
                                t2.exon_cov.get(eidx).copied().unwrap_or(t2.coverage);
                            if pred_exoncov * len > totalcov / 2.0 {
                                lowcov = false;
                                break;
                            }
                        }
                        if lowcov {
                            kill!(n2, n1, "strand_conflict_lowcov_bpcov");
                        }
                    }
                } else if t2.exons.len() <= t1.exons.len()
                    && t1.coverage > t2.coverage * drop
                    && {
                    let ok = included_pred(
                        &txs, n1, n2, longreads, bpcov, singlethr, error_perc, true,
                    );
                    if pairwise_target_trace_enabled() && pairwise_trace_target(t1, t2) {
                        eprintln!(
                            "[TRACE_PAIRWISE_TARGET] included_drop n1={} n2={} ok={}",
                            n1, n2, ok
                        );
                    }
                    ok
                } {
                    // the reference assembler/C++ reference tends to keep many internal-start / internal-end models in long-read mode,
                    // but we need to cull them if they are truly redundant and low-coverage to improve precision.
                    // Protection for high-confidence / rescued transcripts.
                if is_rescue_protected(t2) {
                    continue;
                }

                if t1.is_longread && t2.is_longread && t2.exons.len() > 1 && t2.coverage > 5.0 * t1.coverage {
                        // If t2 has significantly higher coverage, t1 is likely a noise fragment.
                        kill!(n1, n2, "included_drop_redundant");
                    } else if t1.is_longread && t2.is_longread && t2.exons.len() > 1 {
                        // Otherwise, keep the variant for sensitivity.
                        continue;
                    } else {
                        kill!(n2, n1, "included_drop");
                    }
                } else if !n1g
                    && n1 != 0
                    && (((!t1.is_longread
                        && t1.exons.len() <= 2
                        && t2.coverage > lowisofrac * t1.coverage * (t1.exons.len() as f64))
                        || ((t2.coverage > singlethr || t1.coverage < singlethr)
                            && t1.coverage < t2.coverage + singlethr))
                        // Multi-exon victims: only let a richer (≥ readthr) longer isoform win.
                        // Otherwise `t1.coverage < singlethr` is almost always true for LR cov~1,
                        // and a ~0.1-cov readthrough model can erase a real minor isoform (STRG.104.7).
                        && (t1.exons.len() == 1 || t2.coverage >= 1.0)
                        && t1.exons.len() < t2.exons.len()
                        && included_pred(
                            &txs, n2, n1, longreads, bpcov, singlethr, error_perc, true,
                        ))
                {
                    kill!(n1, n2, "reverse_included");
                    break;
                }
            } else if !n1g
                && n1 != 0
                && (((!t1.is_longread
                    && t1.exons.len() <= 2
                    && t2.coverage > lowisofrac * t1.coverage * (t1.exons.len() as f64))
                    || ((t2.coverage > singlethr || t1.coverage < singlethr)
                        && t1.coverage < t2.coverage + singlethr))
                    && (t1.exons.len() == 1 || t2.coverage >= 1.0)
                    && t1.exons.len() <= t2.exons.len()
                    && included_pred(
                        &txs, n2, n1, longreads, bpcov, singlethr, error_perc, true,
                    ))
            {
                kill!(n1, n2, "n2g_reverse_included");
                break;
            } else if t2.is_longread
                && !n2g
                // C++ reference — pred[n2]->cov < DROP (per-base cov vs 0.5)
                && t2.coverage < drop
                && t2.exons.len() <= t1.exons.len()
                && included_pred(&txs, n1, n2, longreads, bpcov, singlethr, error_perc, false)
            {
                // C++ reference — eliminate low-coverage long-read predictions
                // included in a higher-scoring prediction (reference=false skips guide protection)
                kill!(n2, n1, "lr_lowcov_included");
            } else if !n2g
                && t2.exons.len() == 1
                // C++ reference — short-read: cov<n1.cov; long-read: only cov<isofrac*n1.cov
                && ((!t2.is_longread && t2.coverage < t1.coverage) || t2.coverage < isofrac * t1.coverage)
                && s1 <= s2
                && e2 <= e1
            {
                kill!(n2, n1, "se_contained");
            } else if t1.is_longread {
                // C++ reference — n1 is long-read
                if !n1g && t1.exons.len() == 1 && t1.coverage < t2.coverage && s2 <= s1 && e1 <= e2
                {
                    // C++ reference — single-exon LR n1 contained in n2
                    kill!(n1, n2, "lr_se_contained_reverse");
                    break;
                } else if t2.is_longread && !n2g {
                    // C++ reference — both long-read
                    if transcript_overlap_sig(t1, t2) {
                        let mut conflict = false;
                        for &k in &bettercov[n2] {
                            if !overlap_sets[n1].contains(&k) {
                                conflict = true;
                                break;
                            }
                        }
                        if conflict {
                            kill!(n2, n1, "lr_bettercov_conflict");
                        } else if (mixed_mode || t1.strand == t2.strand)
                            && ((mixed_mode && t2.coverage < singlethr)
                                || t2.coverage < t1.coverage * error_perc)
                            && s1 <= s2
                            && e2 <= e1
                            && intronic_soft(t2, t1, error_perc)
                        {
                            kill!(n2, n1, "lr_intronic");
                        } else {
                            bettercov[n2].push(n1);
                        }
                    } else if mixed_mode && t2.strand != t1.strand && t2.coverage < singlethr {
                        // C++ reference — mixed mode antisense
                        if s1 <= e2 && s2 <= e1 {
                            kill!(n2, n1, "lr_mixed_antisense");
                        }
                    }
                } else if !t2.is_longread && !n2g {
                    // C++ referencelike — LR n1 kills low-cov SR n2
                    if s1 <= e2 && s2 <= e1 {
                        if t2.coverage < 1.0 / error_perc
                            || (t2.strand == t1.strand && t2.coverage < error_perc * t1.coverage)
                        {
                            kill!(n2, n1, "lr_kills_sr_overlap");
                        }
                    }
                }
            } else {
                // n1 is short-read
                if !n1g && t1.exons.len() == 1 && t1.coverage < t2.coverage && s2 <= s1 && e1 <= e2
                {
                    kill!(n1, n2, "sr_se_contained_reverse");
                    break;
                } else if !t2.is_longread {
                    // C++ referencelike — both short-read bettercov
                    if transcript_overlap_sig(t1, t2) {
                        let mut conflict = false;
                        for &k in &bettercov[n2] {
                            if !overlap_sets[n1].contains(&k) {
                                conflict = true;
                                break;
                            }
                        }
                        if conflict {
                            kill!(n2, n1, "sr_bettercov_conflict");
                        } else if t2.coverage < error_perc * t1.coverage
                            && s1 <= s2
                            && e2 <= e1
                            && intronic_soft(t2, t1, error_perc)
                        {
                            kill!(n2, n1, "sr_intronic_soft");
                        } else {
                            bettercov[n2].push(n1);
                        }
                    }
                } else if !n2g {
                    // C++ reference — SR n1, LR n2, low cov
                    if s1 <= e2 && s2 <= e1 {
                        if t2.coverage < 1.0 / error_perc
                            || (t2.strand == t1.strand && t2.coverage < error_perc * t1.coverage)
                        {
                            kill!(n2, n1, "sr_kills_lr_overlap");
                        }
                    }
                }
            }
        }
    }

    let total_removed = dead.count_ones();
    if verbose && cov_updates > 0 {
        eprintln!(
            "    CMaxIntv coverage recompute: adjusted {} prediction(s) before pairwise filtering",
            cov_updates
        );
    }
    if parity_debug {
        for (i, tx) in txs.iter().enumerate() {
            if dead.contains(i) {
                let start = tx.exons.first().map(|(s, _)| *s).unwrap_or(0);
                let tlen = tx_exonic_len(tx);
                let killer = kill_by[i];
                let killer_start = if killer < txs.len() {
                    txs[killer].exons.first().map(|(s, _)| *s).unwrap_or(0)
                } else {
                    0
                };
                let killer_cov = if killer < txs.len() {
                    txs[killer].coverage
                } else {
                    0.0
                };
                eprintln!("PARITY_OVERLAP start={} exons={} cov={:.4} tlen={} killer_start={} killer_exons={} killer_cov={:.4} reason={}",
                    start, tx.exons.len(), tx.coverage, tlen,
                    killer_start,
                    if killer < txs.len() { txs[killer].exons.len() } else { 0 },
                    killer_cov,
                    kill_reason[i]);
            }
        }
    }
    if verbose && total_removed > 0 {
        eprintln!(
            "    Pairwise overlap filter: removed {} transcripts (C++ reference pairwise parity)",
            total_removed
        );
    }
    pairwise_summary.exact_reason_note = format_reason_histogram(&reason_counts);
    (
        txs.into_iter()
            .enumerate()
            .filter(|(i, _)| !dead.contains(*i))
            .map(|(_, t)| t)
            .collect(),
        pairwise_summary,
    )
}

/// Remove transcripts that appear to be retained introns of higher-coverage transcripts.
/// n2 is removed if n2_cov < error_perc * n1_cov and any n2 exon fully spans an n1 intron.
/// Corresponds to C++ reference has_retained_intron (7558-7568), used in merge_transfrags (7742-7743).
pub fn retained_intron_filter(
    transcripts: Vec<Transcript>,
    error_perc: f64,
    verbose: bool,
) -> Vec<Transcript> {
    if transcripts.len() <= 1 {
        return transcripts;
    }
    let mut sorted = transcripts;
    // For long-read transcripts, sort by longcov (pre-depletion read mass) rather than
    // EK-flow coverage.  EK-flow coverage is normalized by noderate and can place single-read
    // retained-intron artifacts (longcov=1) *above* dominant isoforms (longcov=338) when the
    // artifact's path happens to have a slightly higher noderate.  Using longcov ensures the
    // truly dominant transcript (highest read support) is ranked first, so retained-intron
    // noise is correctly identified as lower-priority.
    let all_longread = sorted.iter().all(|t| t.is_longread);
    if all_longread {
        sorted.sort_by(|a, b| {
            b.longcov
                .partial_cmp(&a.longcov)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    } else {
        sorted.sort_by(|a, b| {
            b.coverage
                .partial_cmp(&a.coverage)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    let tx_introns: Vec<Vec<Junction>> = sorted
        .iter()
        .map(|tx| {
            if tx.exons.len() >= 2 {
                exons_to_junction_chain(&tx.exons)
            } else {
                vec![]
            }
        })
        .collect();

    let n = sorted.len();
    let mut remove_set: HashSet<usize> = Default::default();

    for j in 0..n {
        if remove_set.contains(&j) {
            continue;
        }
        let n2 = &sorted[j];
        if is_rescue_protected(n2) {
            continue;
        }
        let n2_exons = &n2.exons;
        let n2_cov = n2.coverage;
        if n2_exons.is_empty() {
            continue;
        }
        for i in 0..j {
            if remove_set.contains(&i) {
                continue;
            }
            if sorted[i].chrom != n2.chrom || sorted[i].strand != n2.strand {
                continue;
            }
            // Do not pre-filter by exon count.  Retained-intron transcripts often have MORE
            // exons than the reference because the intron region appears as an extra exon in
            // the path graph (e.g., n2 has exon spanning n1's intron → gffcompare m-class).
            // The structural check below (exon spans intron) is the correct gate; exon-count
            // filtering was masking the majority of retained-intron artifacts.
            if !transcript_overlap_sig(&sorted[i], n2) {
                continue;
            }
            let n1_introns = &tx_introns[i];
            if n1_introns.is_empty() {
                continue;
            }
            let n1_cov = sorted[i].coverage;
            // Require clearly lower abundance than the dominant model (not just <10%): the
            // reference-style error_perc gate is too loose for long-read bundles with real minor
            // isoforms ~5–15% of the major transcript.
            let retained_cov_ratio = error_perc * 0.35;

            // For long-read transcripts, EK-flow coverage is normalized by noderate and
            // doesn't accurately reflect read-count differences.  A dominant transcript
            // with 338 supporting reads can have coverage≈1.44x while a single-read
            // retained-intron artifact has coverage≈1.3x — the ratio ≈0.9 far exceeds
            // the 0.035 threshold, so coverage-based filtering never fires.
            // Use longcov (pre-depletion read mass) instead when both transcripts are
            // long-read: 1/338 ≈ 0.003 << 0.035, so the filter correctly fires for
            // multi-read artifacts as well as single-read ones.
            //
            // For single-read (longcov < 2) artifacts in LOW-COVERAGE loci (dominant
            // longcov < 29), the longcov ratio check can fail (e.g., 1/9 = 0.11 >> 0.035).
            // But these artifacts often have very low EK-flow coverage relative to the
            // dominant (e.g., 15/266 = 5.6% < 10%).  Use an OR of the two gates:
            // either low longcov ratio OR low coverage ratio (at original error_perc=0.1)
            // catches the artifact.
            let passes_abundance_gate = if n2.is_longread && sorted[i].is_longread && sorted[i].longcov > 0.0 {
                // Primary gate: longcov ratio (reliable pre-flow read count), works for all LR.
                // Fallback: coverage ratio for single-read artifacts in low-coverage loci.
                n2.longcov < retained_cov_ratio * sorted[i].longcov
                    || (n2.longcov < 2.0 && n2_cov < error_perc * n1_cov)
            } else {
                n2_cov < retained_cov_ratio * n1_cov
            };

            if passes_abundance_gate {
                for &(exon_start, exon_end) in n2_exons {
                    for &(intron_start, intron_end) in n1_introns.iter() {
                        if exon_start <= intron_start && exon_end >= intron_end {
                            remove_set.insert(j);
                            break;
                        }
                    }
                    if remove_set.contains(&j) {
                        break;
                    }
                }
            }
            if remove_set.contains(&j) {
                break;
            }
        }
    }

    let total_removed = remove_set.len();
    if verbose && total_removed > 0 {
        eprintln!(
            "    Retained intron filter (cov_ratio={:.4}): removed {} transcripts",
            error_perc * 0.35,
            total_removed
        );
    }
    sorted
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !remove_set.contains(i))
        .map(|(_, t)| t)
        .collect()
}

/// Remove single-exon fragments that are isolated and lie near multi-exon transcript ends
/// (polymerase run-off). Isolated = no overlap with any other transcript. Near = within runoff_dist.
/// C++ reference 19247-19276: pred[i] single-exon, alone between p and c; if prev multi-exon and
/// pred[p]->end+runoff>=pred[i]->start and pred[i]->cov < exoncov+singlethr then eliminate; same for next.
/// We use transcript coverage (like reference script); C++ reference uses get_cov on the adjacent exon.
pub fn polymerase_runoff_filter(
    transcripts: Vec<Transcript>,
    runoff_dist: u64,
    singlethr: f64,
    bpcov: Option<&Bpcov>,
    strict_longrec_port: bool,
    verbose: bool,
) -> Vec<Transcript> {
    if transcripts.len() <= 1 {
        return transcripts;
    }
    let mut sorted = transcripts;
    sorted.sort_by_key(|t| t.exons.first().map(|e| e.0).unwrap_or(0));

    let n = sorted.len();
    let mut dead = SmallBitset::with_capacity(n.min(64));
    let mut total_removed = 0usize;
    for i in 0..n {
        if dead.contains(i) {
            continue;
        }
        let tx = &sorted[i];
        if tx.exons.len() != 1 {
            continue;
        }
        if strict_longrec_port && is_guide_pair(tx) {
            continue;
        }
        let tx_start = tx.exons[0].0;
        let tx_end = tx.exons[0].1;
        let tx_cov = tx.coverage;

        let mut p: i64 = i as i64 - 1;
        while p >= 0 && dead.contains(p as usize) {
            p -= 1;
        }
        let mut c = i + 1;
        while c < n && dead.contains(c) {
            c += 1;
        }

        let no_prev_overlap =
            p < 0 || sorted[p as usize].exons.last().map(|e| e.1).unwrap_or(0) < tx_start;
        let no_next_overlap =
            c >= n || sorted[c].exons.first().map(|e| e.0).unwrap_or(u64::MAX) > tx_end;
        if !(no_prev_overlap && no_next_overlap) {
            continue;
        }

        if p >= 0 {
            let prev = &sorted[p as usize];
            let prev_end = prev.exons[prev.exons.len() - 1].1;
            if prev_end.saturating_add(runoff_dist) >= tx_start {
                // C++ reference: use bpcov per-base avg of prev's last exon
                let prev_cov = if strict_longrec_port {
                    if let Some(bp) = bpcov {
                        let (es, ee) = prev.exons[prev.exons.len() - 1];
                        cov_avg_all(bp, es, ee)
                    } else {
                        prev.coverage
                    }
                } else {
                    prev.coverage
                };
                if tx_cov < prev_cov + singlethr {
                    dead.insert_grow(i);
                    total_removed += 1;
                    continue;
                }
            }
        }

        if c < n {
            let next = &sorted[c];
            let next_start = next.exons[0].0;
            if tx_end.saturating_add(runoff_dist) >= next_start {
                // C++ reference: use bpcov per-base avg of next's first exon
                let next_cov = if strict_longrec_port {
                    if let Some(bp) = bpcov {
                        let (es, ee) = next.exons[0];
                        cov_avg_all(bp, es, ee)
                    } else {
                        next.coverage
                    }
                } else {
                    next.coverage
                };
                if tx_cov < next_cov + singlethr {
                    dead.insert_grow(i);
                    total_removed += 1;
                }
            }
        }
    }

    if verbose && total_removed > 0 {
        eprintln!(
            "    Polymerase run-off filter: removed {} single-exon fragments (dist={}bp)",
            total_removed, runoff_dist
        );
    }
    sorted
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !dead.contains(*i))
        .map(|(_, t)| t)
        .collect()
}

/// Polymerase run-on filter: identifies transcripts that bridge two distinct loci (read-throughs).
/// the reference assembler (RareClass::polyRunOn) detects these by looking for a sharp coverage drop in an internal exon.
/// If a large internal exon has a per-base coverage 'dip' significantly lower than the transcript's 
/// overall coverage, it may be split into two.
pub fn polymerase_runon_filter(
    transcripts: Vec<Transcript>,
    bpcov: Option<&Bpcov>,
    verbose: bool,
) -> Vec<Transcript> {
    if transcripts.is_empty() {
        return transcripts;
    }
    let Some(bp) = bpcov else {
        return transcripts;
    };

    let mut out = Vec::with_capacity(transcripts.len());
    let mut total_splits = 0usize;

    for mut tx in transcripts {
        if tx.exons.len() < 2 || is_guide_pair(&tx) {
            out.push(tx);
            continue;
        }

        // We look for a split point in an internal exon.
        let mut split_exon_idx: Option<usize> = None;
        let mut split_pos: Option<u64> = None;

        // the reference assembler RareClass logic: internal exon must be > 500bp and have a clear dip.
        // For now, we use a simpler heuristic: find an internal exon with a region of very low coverage.
        for i in 1..tx.exons.len() - 1 {
            let (es, ee) = tx.exons[i];
            let len = ee - es;
            if len < 500 {
                continue;
            }

            // Sample coverage in the middle of the exon.
            // In a real read-through, the "bridge" is usually the middle part of a long internal exon.
            let mid = es + len / 2;
            let sample_win = 50; // 50bp window
            let win_s = mid.saturating_sub(sample_win / 2);
            let win_e = (mid + sample_win / 2).min(ee);
            
            let local_cov = cov_avg_all(bp, win_s, win_e);
            
            // If local coverage is < 20% of the transcript's avg coverage and < 1.0, consider it a dip.
            if local_cov < tx.coverage * 0.2 && local_cov < 1.0 {
                split_exon_idx = Some(i);
                split_pos = Some(mid);
                break;
            }
        }

        if let (Some(idx), Some(pos)) = (split_exon_idx, split_pos) {
            // Split into two transcripts. 
            // TX1: exons 0..idx, with exon[idx] ending at pos.
            // TX2: exons idx..end, with exon[idx] starting at pos.
            let mut tx2 = tx.clone();
            
            // Adjust TX1
            tx.exons.truncate(idx + 1);
            tx.exons[idx].1 = pos;
            
            // Adjust TX2
            tx2.exons = tx2.exons.split_off(idx);
            tx2.exons[0].0 = pos;
            
            out.push(tx);
            out.push(tx2);
            total_splits += 1;
        } else {
            out.push(tx);
        }
    }

    if verbose && total_splits > 0 {
        eprintln!(
            "    Polymerase run-on filter: split {} transcripts into {} fragments",
            total_splits, total_splits * 2
        );
    }
    out
}

/// Equal prediction merge: deduplicate transcripts with identical intron chains (same chrom/strand).
pub fn dedup_equal_intron_chains(transcripts: Vec<Transcript>, verbose: bool) -> Vec<Transcript> {
    if transcripts.len() < 2 {
        return transcripts;
    }
    // Group multi-exon transcripts by (chrom, strand, junction_chain)
    let mut chain_best: HashMap<(String, char, Vec<Junction>), (usize, f64, u32)> = Default::default();
    for (i, tx) in transcripts.iter().enumerate() {
        if tx.exons.len() < 2 {
            continue;
        }
        let chain = exons_to_junction_chain(&tx.exons);
        let key = (tx.chrom.clone(), tx.strand, chain);
        let score = (if tx.hardstart { 1 } else { 0 }) + (if tx.hardend { 1 } else { 0 });

        let entry = chain_best.entry(key).or_insert((i, tx.coverage, score));
        if score > entry.2 || (score == entry.2 && tx.coverage > entry.1) {
            *entry = (i, tx.coverage, score);
        }
    }

    let n_tx = transcripts.len();
    let mut keep = SmallBitset::with_capacity(n_tx.min(64));
    // Keep all single-exon transcripts
    for (i, tx) in transcripts.iter().enumerate() {
        if tx.exons.len() < 2 {
            keep.insert_grow(i);
        }
    }
    // Keep best copy of each unique chain
    for (_key, (best_idx, _, _)) in &chain_best {
        keep.insert_grow(*best_idx);
    }

    let removed = n_tx - keep.count_ones();
    if verbose && removed > 0 {
        eprintln!(
            "    Equal prediction merge: removed {} duplicate intron chains",
            removed
        );
    }
    transcripts
        .into_iter()
        .enumerate()
        .filter(|(i, _)| keep.contains(*i))
        .map(|(_, t)| t)
        .collect()
}

/// Overlap consistency filter (print_predcluster bettercov logic around C++ reference 18871-18882):
/// if a lower-coverage transcript overlaps multiple higher-coverage transcripts that do not
/// overlap each other, the lower-coverage one is dropped.
pub fn overlap_consistency_filter(transcripts: Vec<Transcript>, verbose: bool) -> Vec<Transcript> {
    if transcripts.len() < 3 {
        return transcripts;
    }
    let mut sorted = transcripts;
    sorted.sort_by(|a, b| {
        b.coverage
            .partial_cmp(&a.coverage)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let n = sorted.len();
    let mut remove_set: HashSet<usize> = Default::default();
    let mut bettercov: Vec<Vec<usize>> = vec![Vec::new(); n];
    const ERROR_PERC: f64 = 0.1;
    let overlaps = build_significant_overlap_matrix(&sorted, ERROR_PERC);

    for i in 0..n {
        if remove_set.contains(&i) {
            continue;
        }
        if sorted[i].exons.len() < 2 || is_guide_tx(&sorted[i]) {
            continue;
        }
        for j in (i + 1)..n {
            if remove_set.contains(&j) {
                continue;
            }
            if sorted[j].exons.len() < 2 || is_guide_tx(&sorted[j]) {
                continue;
            }
            if !overlaps[i].contains(j) || !transcript_overlap_sig(&sorted[i], &sorted[j]) {
                continue;
            }
            let mut conflict = false;
            for &k in &bettercov[j] {
                if !overlaps[i].contains(k) || !transcript_overlap_sig(&sorted[i], &sorted[k]) {
                    conflict = true;
                    break;
                }
            }
            if conflict {
                remove_set.insert(j);
            } else {
                bettercov[j].push(i);
            }
        }
    }

    let removed = remove_set.len();
    if verbose && removed > 0 {
        eprintln!(
            "    Overlap consistency filter: removed {} conflicting overlap transcript(s)",
            removed
        );
    }
    sorted
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !remove_set.contains(i))
        .map(|(_, t)| t)
        .collect()
}

fn is_guide_tx(tx: &Transcript) -> bool {
    guide_equiv_id(tx).is_some()
}

fn guide_id(tx: &Transcript) -> Option<&str> {
    guide_equiv_id(tx)
}

fn tx_len(tx: &Transcript) -> u64 {
    tx.exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum()
}

/// C++ reference equal_pred parity check:
/// - same strand
/// - genomic overlap
/// - same exon count
/// - all internal boundaries equal (all starts except first, all ends except last)
fn equal_pred(a: &Transcript, b: &Transcript) -> bool {
    if a.strand != b.strand || a.exons.len() != b.exons.len() || a.exons.is_empty() {
        return false;
    }
    let (as0, ae0) = (a.exons.first().unwrap().0, a.exons.last().unwrap().1);
    let (bs0, be0) = (b.exons.first().unwrap().0, b.exons.last().unwrap().1);
    if !overlaps_half_open(as0, ae0, bs0, be0) {
        return false;
    }
    // For long-read mode, also check terminal boundaries to avoid collapsing
    // distinct isoforms with different TSS/TES
    if as0 != bs0 || ae0 != be0 {
        return false;
    }
    let nex = a.exons.len();
    for i in 0..nex {
        if i > 0 && a.exons[i].0 != b.exons[i].0 {
            return false;
        }
        if i + 1 < nex && a.exons[i].1 != b.exons[i].1 {
            return false;
        }
    }
    true
}

/// Collapse equal predictions following C++ ref's equal_pred block in printResults (~20579+).
/// Distinct from intron-chain dedup: this keeps reference-style terminal coordinate reconciliation.
pub fn collapse_equal_predictions(transcripts: Vec<Transcript>, verbose: bool) -> Vec<Transcript> {
    if transcripts.len() < 2 {
        return transcripts;
    }
    let mut txs = transcripts;
    txs.sort_by_key(|t| t.exons.first().map(|e| e.0).unwrap_or(0));
    let n = txs.len();
    let mut dead = SmallBitset::with_capacity(n.min(64));
    let mut merged = 0usize;
    const ERROR_PERC: f64 = 0.1;

    for nidx in 0..n.saturating_sub(1) {
        if dead.contains(nidx) {
            continue;
        }
        let n_end = txs[nidx].exons.last().map(|e| e.1).unwrap_or(0);
        let mut ndel = false;
        let mut midx = nidx + 1;
        while !ndel && midx < n {
            if dead.contains(midx) {
                midx += 1;
                continue;
            }
            let m_start = txs[midx].exons.first().map(|e| e.0).unwrap_or(0);
            if m_start > n_end {
                break;
            }
            if !equal_pred(&txs[nidx], &txs[midx]) {
                midx += 1;
                continue;
            }

            let n_gid = guide_id(&txs[nidx]);
            let m_gid = guide_id(&txs[midx]);
            if n_gid.is_some() && m_gid.is_some() && n_gid != m_gid {
                midx += 1;
                continue;
            }

            // Keep the stronger representative and drop the other.
            let n_guide = is_guide_tx(&txs[nidx]);
            let m_guide = is_guide_tx(&txs[midx]);
            let prefer_n = if n_guide != m_guide {
                n_guide
            } else if (txs[nidx].coverage - txs[midx].coverage).abs() > EPSILON {
                txs[nidx].coverage > txs[midx].coverage
            } else {
                tx_len(&txs[nidx]) >= tx_len(&txs[midx])
            };
            let (winner_idx, loser_idx) = if prefer_n { (nidx, midx) } else { (midx, nidx) };
            let winner = txs[winner_idx].clone();
            let loser = txs[loser_idx].clone();

            // reference-like reconciliation:
            // single exon: stitch/merge ranges and weighted coverage;
            // multi exon: transfer terminal boundaries from stronger non-guide and add coverage by length.
            if winner.exons.len() == 1 && loser.exons.len() == 1 {
                let mut keep = txs[midx].clone();
                let ws = winner.exons[0].0.min(loser.exons[0].0);
                let we = winner.exons[0].1.max(loser.exons[0].1);
                let lw = (winner.exons[0].1 - winner.exons[0].0).max(1) as f64;
                let ll = (loser.exons[0].1 - loser.exons[0].0).max(1) as f64;
                keep.exons[0] = (ws, we);
                if winner.coverage > ERROR_PERC * loser.coverage
                    && loser.coverage > ERROR_PERC * winner.coverage
                {
                    let denom = (we - ws).max(1) as f64;
                    keep.coverage = (winner.coverage * lw + loser.coverage * ll) / denom;
                } else {
                    keep.coverage = winner.coverage + loser.coverage * (ll / lw.max(1.0));
                }
                if keep.strand == '.' {
                    keep.strand = winner.strand;
                }
                txs[midx] = keep;
            } else {
                // C++ reference — multi-exon equal_pred merge with per-exon
                // coverage reconciliation.
                // C++: save pre-transfer first/last exon lengths from pred[m],
                // optionally transfer coordinates from higher-cov/guide to m,
                // then add killed pred[n]'s exoncov to m with length scaling.
                let nex = txs[midx].exons.len();
                let flen = (txs[midx].exons[0].1.saturating_sub(txs[midx].exons[0].0)).max(1);
                let llen = (txs[midx].exons[nex - 1]
                    .1
                    .saturating_sub(txs[midx].exons[nex - 1].0))
                .max(1);
                let mut keep = txs[midx].clone();
                if !is_guide_tx(&keep) {
                    keep.exons[0].0 = winner.exons[0].0;
                    let last = keep.exons.len() - 1;
                    let wlast = winner.exons.len() - 1;
                    keep.exons[last].1 = winner.exons[wlast].1;
                }
                if keep.source.is_none() {
                    keep.source = winner.source.clone();
                }
                // The killed prediction (nidx) is always the exoncov source (C++: pred[n]).
                let killed = &txs[nidx];
                let killed_nex = killed.exons.len();
                // Per-exon coverage reconciliation (C++ ref:20679-20700).
                let mut addcov = 0.0f64;
                if killed_nex == nex && !keep.exon_cov.is_empty() {
                    // First exon
                    let k_first_len = killed.exons[0].1.saturating_sub(killed.exons[0].0);
                    let k_first_cov = killed.exon_cov.first().copied().unwrap_or(killed.coverage);
                    if k_first_len < flen {
                        let contrib = k_first_cov * k_first_len as f64;
                        addcov += contrib;
                        keep.exon_cov[0] += contrib / flen as f64;
                    } else {
                        addcov += k_first_cov * flen as f64;
                        keep.exon_cov[0] += k_first_cov;
                    }
                    // Last exon
                    let k_last_len = killed.exons[nex - 1]
                        .1
                        .saturating_sub(killed.exons[nex - 1].0);
                    let k_last_cov = killed
                        .exon_cov
                        .get(nex - 1)
                        .copied()
                        .unwrap_or(killed.coverage);
                    if k_last_len < llen {
                        let contrib = k_last_cov * k_last_len as f64;
                        addcov += contrib;
                        keep.exon_cov[nex - 1] += contrib / llen as f64;
                    } else {
                        addcov += k_last_cov * llen as f64;
                        keep.exon_cov[nex - 1] += k_last_cov;
                    }
                    // Internal exons
                    for k in 1..nex.saturating_sub(1) {
                        let k_cov = killed.exon_cov.get(k).copied().unwrap_or(killed.coverage);
                        let k_len = killed.exons[k].1.saturating_sub(killed.exons[k].0);
                        keep.exon_cov[k] += k_cov;
                        addcov += k_cov * k_len as f64;
                    }
                    let tlen = tx_len(&keep).max(1) as f64;
                    keep.coverage += addcov / tlen;
                } else {
                    // Fallback: simple length-ratio coverage addition.
                    let lwin = tx_len(&winner).max(1) as f64;
                    let llos = tx_len(&loser).max(1) as f64;
                    keep.coverage = winner.coverage + loser.coverage * (llos / lwin);
                }
                txs[midx] = keep;
            }

            dead.insert_grow(nidx);
            ndel = true;
            merged += 1;
        }
    }

    let out: Vec<Transcript> = txs
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !dead.contains(*i))
        .map(|(_, t)| t)
        .collect();
    if verbose && merged > 0 {
        eprintln!("    equal_pred collapse: merged {} prediction(s)", merged);
    }
    out
}

/// Port of high-impact single-exon handling from C++ reference print_predcluster:
/// - stitch nearby single-exon predictions with similar coverage (guide-aware)
/// - drop low-coverage single-exon non-guide predictions (< singlethr)
pub fn collapse_single_exon_runoff(
    transcripts: Vec<Transcript>,
    singlethr: f64,
    runoff_dist: u64,
    verbose: bool,
) -> Vec<Transcript> {
    // C++ parity: don't early return - single-exon filtering applies even to single transcripts
    // The len < 2 check was causing single-exon transcripts to bypass filtering
    if transcripts.is_empty() {
        return transcripts;
    }
    let mut txs = transcripts;
    txs.sort_by_key(|t| t.exons.first().map(|e| e.0).unwrap_or(0));
    let n = txs.len();
    let mut dead = SmallBitset::with_capacity(n.min(64));
    let mut stitched = 0usize;
    let mut dropped = 0usize;
    const ERROR_PERC: f64 = 0.1;

    for i in 0..n {
        if dead.contains(i) {
            continue;
        }
        if txs[i].exons.len() != 1 {
            continue;
        }
        let (is, ie) = txs[i].exons[0];

        let mut j = i + 1;
        while j < n && dead.contains(j) {
            j += 1;
        }

        // Stitch with nearby single exon if similar abundance and non-guide
        let mut did_stitch = false;
        if j < n && txs[j].exons.len() == 1 && txs[i].strand == txs[j].strand {
            let (js, je) = txs[j].exons[0];
            let near = js <= ie.saturating_add(runoff_dist);
            let comparable = txs[j].coverage > ERROR_PERC * txs[i].coverage
                && txs[i].coverage > ERROR_PERC * txs[j].coverage;
            let gi = guide_id(&txs[i]);
            let gj = guide_id(&txs[j]);
            let stitch_ok = if gi.is_some() || gj.is_some() {
                gi.is_some() && gj.is_some() && gi == gj
            } else {
                true
            };
            if near && comparable && stitch_ok {
                let new_start = is.min(js);
                let new_end = ie.max(je);
                let li = ie.saturating_sub(is).max(1) as f64;
                let lj = je.saturating_sub(js).max(1) as f64;
                let new_len = new_end.saturating_sub(new_start).max(1) as f64;
                txs[j].exons[0] = (new_start, new_end);
                txs[j].coverage = (txs[i].coverage * li + txs[j].coverage * lj) / new_len;
                dead.insert_grow(i);
                stitched += 1;
                did_stitch = true;
            }
        }

        // C++ parity: coverage comparison uses per-base coverage
        let tx_len = txs[i].exons.iter().map(|(s, e)| e - s).sum::<u64>().max(1) as f64;
        let cov_per_base = txs[i].coverage / tx_len;
        if !did_stitch && !is_guide_tx(&txs[i]) && cov_per_base < singlethr {
            dead.insert_grow(i);
            dropped += 1;
        }
    }

    let out: Vec<Transcript> = txs
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !dead.contains(*i))
        .map(|(_, t)| t)
        .collect();
    if verbose && (stitched > 0 || dropped > 0) {
        eprintln!(
            "    single_exon_runoff: stitched={} dropped_lowcov={}",
            stitched, dropped
        );
    }
    out
}

fn split_overlap_midpoint(a: &mut (u64, u64), b: &mut (u64, u64), max_ovl: u64) -> bool {
    let ov_start = a.0.max(b.0);
    let ov_end = a.1.min(b.1);
    if ov_end <= ov_start {
        return false;
    }
    let ov_len = ov_end - ov_start;
    if ov_len > max_ovl {
        return false;
    }
    let cut = ov_start + ov_len / 2;
    if cut <= a.0 || cut >= b.1 {
        return false;
    }
    if cut >= a.1 || cut <= b.0 {
        return false;
    }
    a.1 = cut;
    b.0 = cut;
    a.1 > a.0 && b.1 > b.0
}

/// Strand-conflict endpoint split:
/// for overlapping opposite-strand multi-exon predictions, split small terminal overlaps
/// at midpoint instead of forcing downstream drop decisions.
pub fn split_strand_conflict_endpoints(
    transcripts: Vec<Transcript>,
    verbose: bool,
) -> Vec<Transcript> {
    if transcripts.len() < 2 {
        return transcripts;
    }
    let mut txs = transcripts;
    let mut splits = 0usize;
    let max_ovl = 2 * LONGINTRONANCHOR;
    for i in 0..txs.len() {
        for j in (i + 1)..txs.len() {
            let (left, right) = txs.split_at_mut(j);
            let a = &mut left[i];
            let b = &mut right[0];
            if a.strand == b.strand || a.exons.len() < 2 || b.exons.len() < 2 {
                continue;
            }
            let (as0, ae0) = match (a.exons.first(), a.exons.last()) {
                (Some(s), Some(e)) => (s.0, e.1),
                _ => continue,
            };
            let (bs0, be0) = match (b.exons.first(), b.exons.last()) {
                (Some(s), Some(e)) => (s.0, e.1),
                _ => continue,
            };
            if !overlaps_half_open(as0, ae0, bs0, be0) {
                continue;
            }
            let al = a.exons.len() - 1;
            if split_overlap_midpoint(&mut a.exons[al], &mut b.exons[0], max_ovl) {
                splits += 1;
            }
            let bl = b.exons.len() - 1;
            if split_overlap_midpoint(&mut b.exons[bl], &mut a.exons[0], max_ovl) {
                splits += 1;
            }
        }
    }
    let mut out = Vec::with_capacity(txs.len());
    for mut tx in txs {
        tx.exons.retain(|(s, e)| e > s);
        if !tx.exons.is_empty() {
            out.push(tx);
        }
    }
    if verbose && splits > 0 {
        eprintln!(
            "    strand_conflict_split: applied {} endpoint split(s)",
            splits
        );
    }
    out
}

/// Remove transcripts whose intron chain is a proper subset of another transcript chain
/// (same chrom/strand), keeping the higher-coverage / longer-chain container.
pub fn dedup_subset_intron_chains(transcripts: Vec<Transcript>, verbose: bool) -> Vec<Transcript> {
    if transcripts.len() < 2 {
        return transcripts;
    }
    let n = transcripts.len();
    let chains: Vec<Vec<Junction>> = transcripts
        .iter()
        .map(|tx| exons_to_junction_chain(&tx.exons))
        .collect();
    let mut drop = SmallBitset::with_capacity(n.min(64));

    for i in 0..n {
        if drop.contains(i) || chains[i].is_empty() {
            continue;
        }
        let t1 = &transcripts[i];
        if is_rescue_protected(t1) {
            continue;
        }
        // Keep shorter intron-prefix isoforms with strong direct read support; otherwise a
        // longer readthrough / extended 3' model removes a valid minor isoform.
        // Threshold 3.0: longcov is now pre-depletion read_count, so this requires ≥3 reads.
        if t1.is_longread && t1.longcov >= 3.0 {
            continue;
        }

        for j in 0..n {
            if i == j || drop.contains(i) {
                continue;
            }
            let t2 = &transcripts[j];
            
            if t1.chrom != t2.chrom || t1.strand != t2.strand {
                continue;
            }
            
            // Check if i is a subset of j
            if !junction_chain_is_subset(&chains[i], &chains[j]) {
                continue;
            }
            
            let cov_i = t1.coverage;
            let cov_j = t2.coverage;

            // A larger transcript (j) is preferred over a subset (i) if:
            // 1. It has higher coverage.
            // 2. It has comparable coverage (>= 50% of subset) AND better boundary support.

            let score_i = (if t1.hardstart { 1 } else { 0 }) + (if t1.hardend { 1 } else { 0 });
            let score_j = (if t2.hardstart { 1 } else { 0 }) + (if t2.hardend { 1 } else { 0 });

            let prefer_j = cov_j > cov_i
                || (cov_j >= cov_i * 0.5 && score_j > score_i);
                
            if prefer_j {
                drop.insert_grow(i);
            }
        }
    }

    let removed = drop.count_ones();
    if verbose && removed > 0 {
        eprintln!(
            "    Subset intron-chain dedup: removed {} subset transcript(s)",
            removed
        );
    }
    transcripts
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !drop.contains(*i))
        .map(|(_, t)| t)
        .collect()
}

/// Collapse near-duplicate transcripts that share the same intron chain and differ only by
/// small terminal boundary jitter (common with long-read alignment clipping).
/// Keeps the higher-coverage representative.
pub fn collapse_near_equal_intron_chains(
    transcripts: Vec<Transcript>,
    terminal_tolerance: u64,
    verbose: bool,
) -> Vec<Transcript> {
    if transcripts.len() < 2 || terminal_tolerance == 0 {
        return transcripts;
    }
    let mut txs = transcripts;
    txs.sort_by(|a, b| {
        let a0 = a.exons.first().map(|e| e.0).unwrap_or(0);
        let b0 = b.exons.first().map(|e| e.0).unwrap_or(0);
        (a.chrom.as_str(), a.strand, a0).cmp(&(b.chrom.as_str(), b.strand, b0))
    });
    let n = txs.len();
    let mut dead = SmallBitset::with_capacity(n.min(64));
    let mut merged = 0usize;

    for i in 0..n {
        if dead.contains(i) {
            continue;
        }
        let ai = &txs[i];
        let a_start = ai.exons.first().map(|e| e.0).unwrap_or(0);
        let a_end = ai.exons.last().map(|e| e.1).unwrap_or(0);
        let a_j = exons_to_junction_chain(&ai.exons);

        for j in (i + 1)..n {
            if dead.contains(j) {
                continue;
            }
            let bj = &txs[j];
            if ai.chrom != bj.chrom || ai.strand != bj.strand {
                if bj.chrom != ai.chrom {
                    break;
                }
                continue;
            }
            let b_j = exons_to_junction_chain(&bj.exons);
            // C++ parity: collapse transcripts with near-identical intron chains.
            // LR alignment noise creates shifted splice site calls (±10bp), producing
            // transcripts that differ only by small junction offsets.  Match intron
            // chains with the same tolerance used for junction coalescing.
            let chain_match = if a_j.len() == b_j.len() {
                a_j.iter().zip(b_j.iter()).all(|(aj, bj)| {
                    aj.0.abs_diff(bj.0) <= terminal_tolerance
                        && aj.1.abs_diff(bj.1) <= terminal_tolerance
                })
            } else {
                false
            };
            if !chain_match {
                continue;
            }
            let b_start = bj.exons.first().map(|e| e.0).unwrap_or(0);
            let b_end = bj.exons.last().map(|e| e.1).unwrap_or(0);
            if a_start.abs_diff(b_start) <= terminal_tolerance
                && a_end.abs_diff(b_end) <= terminal_tolerance
            {
                // Only collapse when the weaker transcript has much lower coverage.
                // Collapsing similarly-supported transcripts risks removing real isoforms
                // whose intron chains happen to be within tolerance of a higher-cov variant.
                let (cov_hi, cov_lo) = if txs[j].coverage > txs[i].coverage {
                    (txs[j].coverage, txs[i].coverage)
                } else {
                    (txs[i].coverage, txs[j].coverage)
                };
                if cov_hi > 0.0 && cov_lo / cov_hi > 0.3 {
                    continue; // Both well-supported, don't collapse
                }
                if txs[j].coverage > txs[i].coverage {
                    dead.insert_grow(i);
                    merged += 1;
                    break;
                } else {
                    dead.insert_grow(j);
                    merged += 1;
                }
            }
        }
    }

    let out: Vec<Transcript> = txs
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !dead.contains(*i))
        .map(|(_, t)| t)
        .collect();
    if verbose && merged > 0 {
        eprintln!(
            "    near_equal_chain_collapse: merged {} transcript(s) (tol={})",
            merged, terminal_tolerance
        );
    }
    out
}

/// Deduplicate transcripts that share an identical intron chain (same splice-site coordinates),
/// keeping the single highest-coverage representative per chain per locus.
///
/// Background: `dedup_subset_intron_chains` handles *strict* subsets (shorter chain ⊆ longer
/// chain) but `junction_chain_is_subset` returns false for equal-length chains, so transcripts
/// with identical intron chains are never removed by that pass. `collapse_near_equal_intron_chains`
/// only collapses when terminal exon positions differ by ≤ junction_correction_window; our flow
/// duplicates may differ by hundreds of bp in terminal position.
///
/// the reference assembler never emits two transcripts from the same locus with the same intron chain —
/// the flow extraction is deterministic and produces exactly one path per splice-graph route.
/// This function enforces the same invariant.
pub fn dedup_exact_intron_chains(transcripts: Vec<Transcript>, verbose: bool) -> Vec<Transcript> {
    use std::collections::HashMap;
    if transcripts.len() < 2 {
        return transcripts;
    }
    // Key: (chrom, strand, intron_chain). Value: index of the best representative so far.
    let mut best: HashMap<(String, char, Vec<Junction>), usize> = HashMap::new();
    let mut drop = SmallBitset::with_capacity(transcripts.len().min(64));
    for (i, tx) in transcripts.iter().enumerate() {
        if tx.exons.len() < 2 {
            continue; // mono-exon: no intron chain to deduplicate on
        }
        let chain = exons_to_junction_chain(&tx.exons);
        if chain.is_empty() {
            continue;
        }
        let key = (tx.chrom.clone(), tx.strand, chain);
        match best.get(&key) {
            None => {
                best.insert(key, i);
            }
            Some(&prev) => {
                // Keep the one with higher coverage; ties go to the earlier index.
                if transcripts[i].coverage > transcripts[prev].coverage {
                    drop.insert_grow(prev);
                    *best.get_mut(&key).unwrap() = i;
                } else {
                    drop.insert_grow(i);
                }
            }
        }
    }
    let removed = drop.count_ones();
    if verbose && removed > 0 {
        eprintln!(
            "    Exact intron-chain dedup: removed {} duplicate transcript(s)",
            removed
        );
    }
    transcripts
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !drop.contains(*i))
        .map(|(_, t)| t)
        .collect()
}

/// Filter transcripts whose intron chains are not witnessed by any read.
///
/// For each predicted transcript with ≥2 introns, check that every pair of consecutive
/// introns appears in at least one read's intron chain. This prevents the flow decomposition
/// from creating novel junction combinations that no single read supports.
///
/// `read_intron_pairs`: set of (donor1, acceptor1, donor2, acceptor2) tuples from all reads
/// in the bundle, pre-built by the caller.
/// `tolerance`: junction coordinate matching tolerance in bp.
pub fn filter_unwitnessed_chains(
    transcripts: Vec<Transcript>,
    read_intron_pairs: &std::collections::HashSet<(u64, u64, u64, u64)>,
    tolerance: u64,
    verbose: bool,
) -> Vec<Transcript> {
    if read_intron_pairs.is_empty() || transcripts.is_empty() {
        return transcripts;
    }
    let before = transcripts.len();
    let out: Vec<Transcript> = transcripts
        .into_iter()
        .filter(|tx| {
            // Guide-matched transcripts are always kept.
            if tx.source.as_deref().map_or(false, |s| s.starts_with("guide:")) {
                return true;
            }
            let introns = exons_to_junction_chain(&tx.exons);
            if introns.len() < 2 {
                return true; // Single-intron transcripts always pass.
            }
            // Check every consecutive pair of introns.
            // Note: local Junction type alias is (u64, u64), not types::Junction.
            for w in introns.windows(2) {
                let key = (w[0].0, w[0].1, w[1].0, w[1].1);
                // Exact match first.
                if read_intron_pairs.contains(&key) {
                    continue;
                }
                // Tolerance match.
                if tolerance > 0 {
                    let found = read_intron_pairs.iter().any(|&(d1, a1, d2, a2)| {
                        d1.abs_diff(key.0) <= tolerance
                            && a1.abs_diff(key.1) <= tolerance
                            && d2.abs_diff(key.2) <= tolerance
                            && a2.abs_diff(key.3) <= tolerance
                    });
                    if found {
                        continue;
                    }
                }
                return false; // Unwitnessed pair found.
            }
            true
        })
        .collect();
    let removed = before - out.len();
    if verbose && removed > 0 {
        eprintln!(
            "    Read-chain witness filter: removed {} unwitnessed transcript(s)",
            removed
        );
    }
    out
}

/// Build the set of consecutive intron pairs from all reads in a bundle.
/// Each pair is (donor1, acceptor1, donor2, acceptor2).
pub fn build_read_intron_pairs(reads: &[crate::types::BundleRead]) -> std::collections::HashSet<(u64, u64, u64, u64)> {
    let mut pairs = std::collections::HashSet::new();
    for read in reads {
        if read.junctions.len() < 2 {
            continue;
        }
        for w in read.junctions.windows(2) {
            pairs.insert((w[0].donor, w[0].acceptor, w[1].donor, w[1].acceptor));
        }
    }
    pairs
}

/// Apply the same prediction filters as the pipeline (C++ reference print_predcluster).
/// Order: pairwise overlap/inclusion filtering → long-read isofrac → runoff/readthr gates.
/// C++ parity: the pairwise block runs for long-read predictions too; only the CMaxIntv
/// coverage redistribution inside `pairwise_overlap_filter` stays disabled for long reads.
pub fn print_predcluster(
    transcripts: Vec<Transcript>,
    config: &RunConfig,
    bpcov: Option<&Bpcov>,
    trace_ref: Option<&RefTranscript>,
) -> Vec<Transcript> {
    print_predcluster_with_summary(transcripts, config, bpcov, trace_ref).0
}

#[derive(Debug, Clone, Default)]
pub struct PredclusterStageSummary {
    pub entry_count: usize,
    pub after_pairwise: usize,
    pub after_isofrac: usize,
    pub after_near_equal_chain_collapse: usize,
    pub after_exact_chain_dedup: usize,
    pub after_dedup_subset_chain: usize,
    pub after_runoff: usize,
    pub after_polymerase_runoff: usize,
    pub after_polymerase_runon: usize,
    pub after_readthr: usize,
    pub pairwise_summary: PairwiseKillSummary,
    pub isofrac_summary: IsofracKillSummary,
}

pub fn print_predcluster_with_summary(
    transcripts: Vec<Transcript>,
    config: &RunConfig,
    bpcov: Option<&Bpcov>,
    trace_ref: Option<&RefTranscript>,
) -> (Vec<Transcript>, PredclusterStageSummary) {
    let trace_stage = |stage: &str, txs: &[Transcript]| {
        if let Some(ref_tx) = trace_ref {
            debug_ref_stage(stage, "predcluster", ref_tx, txs);
            predcluster_trace_dump(stage, ref_tx, txs);
        }
    };
    let fate_trace = pred_fate_trace_enabled();
    // Helper: identify killed transcripts between filter stages via exon fingerprints.
    // Uses (start, end, nexons, coverage) as identity since we don't have stable IDs.
    let tx_key = |t: &Transcript| -> (u64, u64, usize, u64) {
        let s = t.exons.first().map(|e| e.0).unwrap_or(0);
        let e = t.exons.last().map(|e| e.1).unwrap_or(0);
        (s, e, t.exons.len(), t.coverage.to_bits())
    };
    let emit_fate = |stage: &str, before: &[Transcript], after: &[Transcript]| {
        if !fate_trace {
            return;
        }
        let after_keys: HashSet<(u64, u64, usize, u64)> =
            after.iter().map(|t| tx_key(t)).collect();
        let mut killed = 0usize;
        for (i, t) in before.iter().enumerate() {
            if !tx_overlaps_trace_locus(t) {
                continue;
            }
            let k = tx_key(t);
            if !after_keys.contains(&k) {
                eprintln!(
                    "print_predcluster: PRED_FATE pred[{}] killed_by={} {}",
                    i,
                    stage,
                    tx_summary(t)
                );
                killed += 1;
            }
        }
        if killed > 0 {
            eprintln!(
                "print_predcluster: FILTER stage={} killed={} survivors={}",
                stage,
                killed,
                after.len()
            );
        }
    };
    let (mut protected, mut txs): (Vec<_>, Vec<_>) = if config.emit_junction_paths {
        transcripts
            .into_iter()
            .partition(|t| predcluster_protected_source(t.source.as_deref()))
    } else {
        (Vec::new(), transcripts)
    };
    let mut summary = PredclusterStageSummary {
        entry_count: txs.len(),
        ..PredclusterStageSummary::default()
    };
    trace_stage("predcluster.entry", &txs);
    bundle_cov_dump_stage("entry", &txs);
    if fate_trace {
        eprintln!("print_predcluster: ENTER count={}", txs.len());
        emit_pred_entries("BEFORE_PAIRWISE", &txs);
    }
    if std::env::var_os("RUSTLE_ENABLE_SPLIT_STRAND_CONFLICT_ENDPOINTS").is_some() {
        let before = if fate_trace { txs.clone() } else { Vec::new() };
        txs = split_strand_conflict_endpoints(txs, config.verbose);
        emit_fate("split_strand_conflict_endpoints", &before, &txs);
        trace_stage("predcluster.split_strand_conflict_endpoints", &txs);
    }
    // C++ parity (C++ reference): eliminate incomplete transcripts with only guide introns.
    // This filter is always enabled in C++ long-read mode (checkincomplete=true).
    if config.long_reads {
        let before = if fate_trace { txs.clone() } else { Vec::new() };
        txs = eliminate_incomplete_vs_guides(txs, config.verbose);
        emit_fate("eliminate_incomplete_vs_guides", &before, &txs);
        trace_stage("predcluster.eliminate_incomplete_vs_guides", &txs);
    }
    // C++ parity: retained_intron filter (C++ reference, merge_transfrags).
    // Removes transcripts that have retained introns relative to higher-coverage transcripts.
    if config.long_reads {
        let before = if fate_trace { txs.clone() } else { Vec::new() };
        txs = retained_intron_filter(txs, config.pairwise_error_perc, config.verbose);
        emit_fate("retained_intron_filter", &before, &txs);
        trace_stage("predcluster.retained_intron_filter", &txs);
    }
    // Bundle tracing: log all transcripts entering pairwise filter
    if trace_locus_range().is_some() {
        for t in &txs {
            trace_tx_detail("ENTER_PAIRWISE", t, None);
        }
    }
    if !config.max_sensitivity {
        let before_pairwise = if fate_trace { txs.clone() } else { Vec::new() };
        // C++ parity: the readthrough elimination block (C++ reference) is commented out,
        // but the main pairwise overlap block still runs for long-read predictions. The
        // long-read-specific skip is only for CMaxIntv coverage redistribution inside
        // `pairwise_overlap_filter`, not for the filter block itself.
        let (txs_after_pairwise, pairwise_summary) = pairwise_overlap_filter_with_summary(
            txs,
            config.pairwise_isofrac,
            config.lowisofrac,
            config.pairwise_error_perc,
            config.pairwise_drop,
            config.singlethr,
            config.long_reads,
            config.long_reads && config.long_read_min_len > 0,
            bpcov,
            config.verbose,
        );
        txs = txs_after_pairwise;
        summary.pairwise_summary = pairwise_summary;
        summary.after_pairwise = txs.len();
        emit_fate("pairwise_overlap_filter", &before_pairwise, &txs);
        trace_stage("predcluster.pairwise_overlap_filter", &txs);
        bundle_cov_dump_stage("after_pairwise", &txs);
    } else {
        // In diagnostic max-sensitivity mode, keep everything through pairwise to avoid
        // masking extractability issues behind heuristic overlap pruning.
        summary.after_pairwise = txs.len();
    }
    // Bundle tracing: log transcripts after pairwise filter
    if trace_locus_range().is_some() {
        for t in &txs {
            trace_tx_detail("AFTER_PAIRWISE", t, None);
        }
    }
    // C++ parity: longunder isofrac filter (C++ reference, longreads branch).
    // Eliminates low-coverage transcripts relative to the dominant transcript in each interval.
    if config.long_reads {
        let before_isofrac = if fate_trace { txs.clone() } else { Vec::new() };
        let (txs_after_isofrac, isofrac_summary) = isofrac_with_summary(
            txs,
            // the reference assembler -f controls the longunder isofrac threshold as well.
            // Keep pairwise overlap filter separately configurable, but longunder should follow
            // the transcript-level isofrac knob.
            config.transcript_isofrac,
            config.pairwise_error_perc,
            config.pairwise_drop,
            bpcov,
            config.verbose,
            config.transcript_isofrac_keep_min,
        );
        txs = txs_after_isofrac;
        summary.isofrac_summary = isofrac_summary;
        summary.after_isofrac = txs.len();
        emit_fate("isofrac", &before_isofrac, &txs);
        trace_stage("predcluster.isofrac", &txs);
        bundle_cov_dump_stage("after_isofrac", &txs);
        if trace_locus_range().is_some() {
            for t in &txs {
                trace_tx_detail("AFTER_ISOFRAC", t, None);
            }
        }
    } else {
        summary.after_isofrac = txs.len();
    }
    if config.long_reads {
        let before_collapse = if fate_trace { txs.clone() } else { Vec::new() };
        txs = collapse_near_equal_intron_chains(
            txs,
            config.junction_correction_window,
            config.verbose,
        );
        summary.after_near_equal_chain_collapse = txs.len();
        emit_fate("collapse_near_equal_intron_chains", &before_collapse, &txs);
        trace_stage("predcluster.collapse_near_equal_intron_chains", &txs);
        // Deduplicate transcripts that have identical intron chains (same splice sites)
        // but different terminal exon positions. the reference assembler never emits two transcripts
        // from the same locus with the same intron chain; collapse to the best-coverage copy.
        let before_exact = if fate_trace { txs.clone() } else { Vec::new() };
        txs = dedup_exact_intron_chains(txs, config.verbose);
        summary.after_exact_chain_dedup = txs.len();
        emit_fate("dedup_exact_intron_chains", &before_exact, &txs);
        trace_stage("predcluster.dedup_exact_intron_chains", &txs);
    } else {
        summary.after_near_equal_chain_collapse = txs.len();
        summary.after_exact_chain_dedup = txs.len();
    }
    if config.dedup_intron_chain_subsets {
        let before_dedup = if fate_trace { txs.clone() } else { Vec::new() };
        txs = dedup_subset_intron_chains(txs, config.verbose);
        summary.after_dedup_subset_chain = txs.len();
        emit_fate("dedup_subset_intron_chains", &before_dedup, &txs);
        trace_stage("predcluster.dedup_subset_intron_chains", &txs);
    } else {
        summary.after_dedup_subset_chain = txs.len();
    }
    if config.filter_contained {
        let before_contained = if fate_trace { txs.clone() } else { Vec::new() };
        txs = filter_contained_transcripts(txs, config.verbose);
        emit_fate("filter_contained_transcripts", &before_contained, &txs);
        trace_stage("predcluster.filter_contained_transcripts", &txs);
        bundle_cov_dump_stage("after_filter_contained", &txs);
    }
    let before_runoff = if fate_trace { txs.clone() } else { Vec::new() };
    let mut runoff_dist = if config.long_reads { 0 } else { 200 };
    if config.bundle_merge_dist > runoff_dist {
        runoff_dist = config.bundle_merge_dist;
    }
    // C++ parity: singlethr filter applies to all modes (C++ reference)
    let effective_singlethr = config.singlethr;
    txs = collapse_single_exon_runoff(txs, effective_singlethr, runoff_dist, config.verbose);
    summary.after_runoff = txs.len();
    emit_fate("collapse_single_exon_runoff", &before_runoff, &txs);
    trace_stage("predcluster.collapse_single_exon_runoff", &txs);
    if trace_locus_range().is_some() {
        for t in &txs {
            trace_tx_detail("AFTER_RUNOFF", t, None);
        }
    }
    let before_polyrunoff = if fate_trace { txs.clone() } else { Vec::new() };
    txs = polymerase_runoff_filter(
        txs,
        runoff_dist,
        config.singlethr,
        bpcov,
        true,
        config.verbose,
    );
    summary.after_polymerase_runoff = txs.len();
    emit_fate("polymerase_runoff_filter", &before_polyrunoff, &txs);
    trace_stage("predcluster.polymerase_runoff_filter", &txs);
    if trace_locus_range().is_some() {
        for t in &txs {
            trace_tx_detail("AFTER_POLYRUNOFF", t, None);
        }
    }

    let before_runon = if fate_trace { txs.clone() } else { Vec::new() };
    if config.polymerase_runon_filter {
        txs = polymerase_runon_filter(txs, bpcov, config.verbose);
    }
    summary.after_polymerase_runon = txs.len();
    emit_fate("polymerase_runon_filter", &before_runon, &txs);
    trace_stage("predcluster.polymerase_runon_filter", &txs);
    if trace_locus_range().is_some() {
        for t in &txs {
            trace_tx_detail("AFTER_POLYRUNON", t, None);
        }
    }

    let before_readthr = if fate_trace { txs.clone() } else { Vec::new() };
    let before = txs.len();
    // C++ parity: readthr gate applies to ALL predictions (C++ reference).
    // Long-read transcripts are NOT exempt — low-coverage paths must be filtered.
    let readthr_debug = std::env::var_os("RUSTLE_READTHR_DEBUG").is_some();
    if readthr_debug {
        for t in &txs {
            // C++ parity: long-read uses longcov, short-read uses coverage for readthr check
            let check_cov = if t.is_longread { t.longcov } else { t.coverage };
            if !is_guide_pair(t) && check_cov < config.readthr {
                let start = t.exons.first().map(|e| e.0).unwrap_or(0);
                let end = t.exons.last().map(|e| e.1).unwrap_or(0);
                let introns: Vec<String> = t
                    .exons
                    .windows(2)
                    .map(|w| format!("{}-{}", w[0].1, w[1].0))
                    .collect();
                eprintln!(
                    "READTHR_KILL start={} end={} exons={} cov={:.6} longcov={:.1} introns={}",
                    start,
                    end,
                    t.exons.len(),
                    t.coverage,
                    t.longcov,
                    introns.join(",")
                );
            }
        }
    }
    // C++ parity (C++ reference, 20322): readthr/singlethr gate.
    // - All transcripts: cov >= 1.0 (readthr)
    // - Single-exon transcripts: additional cov >= 4.75 (singlethr)
    // - Mixed mode with guides: use singlethr (4.75) for all
    let mixed_mode = config.long_reads && config.long_read_min_len > 0;
    let has_guides = txs.iter().any(|t| t.ref_transcript_id.is_some());
    let base_threshold = if mixed_mode && has_guides {
        config.singlethr // 4.75 in mixed mode with guides
    } else {
        config.readthr // 1.0 default
    };
    txs.retain(|t| {
        if is_guide_pair(t) {
            return true;
        }
        // Special parity: rescued transcripts with at least one high-confidence boundary
        // are kept even if below readthr, emulating the reference assembler's selective promotion.
        if is_rescue_protected(t) && (t.hardstart || t.hardend) {
            return true;
        }
        
        let is_unanchored = !t.hardstart && !t.hardend;
        // Unanchored boost targets short spurious models; long multi-exon isoforms often lack
        // hardstart/hardend flags but still have strong junction support (e.g. STRG.151.7).
        let effective_base_thr = if is_unanchored && t.exons.len() <= 3 {
            base_threshold * 1.5
        } else {
            base_threshold
        };

        // Determine threshold: single-exon uses singlethr (C++ reference)
        let threshold = if t.exons.len() == 1 {
            config.singlethr // 4.75 for single-exon
        } else {
            effective_base_thr 
        };
        if t.coverage >= threshold {
            return true;
        }
        // For long-read transcripts, EK flow coverage is not per-base read depth: it is
        // derived from the max-flow decomposition (nodeflux_abs * noderate) and can be
        // below 1.0 for real minor isoforms supported by multiple reads.  Fall back to
        // longcov (pre-depletion read_count) as a secondary acceptance gate: if the path
        // has ≥2 actual reads supporting it, keep it regardless of flow coverage.
        // Threshold 2.0 prevents single-read (longcov=1.0) noise from bypassing the gate.
        if t.is_longread && t.longcov >= 2.0 {
            return true;
        }
        readthr_allow_longcov_fallback() && t.is_longread && t.longcov >= threshold
    });
    summary.after_readthr = txs.len();
    if config.verbose && txs.len() < before {
        eprintln!(
            "    readthr_gate: removed {} transcript(s) (base cov < {}, single-exon cov < {})",
            before - txs.len(),
            base_threshold,
            config.singlethr
        );
    }
    emit_fate("readthr_gate", &before_readthr, &txs);
    trace_stage("predcluster.readthr_gate", &txs);
    bundle_cov_dump_stage("after_readthr", &txs);

    if trace_locus_range().is_some() {
        for t in &txs {
            trace_tx_detail("AFTER_READTHR", t, None);
        }
    }

    if fate_trace {
        eprintln!("print_predcluster: FINAL_SURVIVORS count={}", txs.len());
        for (i, t) in txs.iter().enumerate() {
            eprintln!(
                "print_predcluster: FINAL_KEEP pred[{}] {}",
                i,
                tx_summary(t)
            );
        }
    }
    // Bundle tracing: log final transcripts
    if trace_locus_range().is_some() {
        for t in &txs {
            trace_tx_detail("FINAL_OUTPUT", t, None);
        }
    }
    if std::env::var_os("RUSTLE_PARITY_PRED_FINAL").is_some()
        || std::env::var_os("RUSTLE_PARITY_DEBUG").is_some()
    {
        for (i, tx) in txs.iter().enumerate() {
            let start = tx.exons.first().map(|(s, _)| *s).unwrap_or(0);
            let end = tx.exons.last().map(|(_, e)| *e).unwrap_or(0);
            let guide_label = tx
                .source
                .as_deref()
                .and_then(|s| s.strip_prefix("guide:"))
                .unwrap_or("novel");
            eprintln!(
                "--- print_predcluster: FINAL_KEEP pred[{}] {}-{} cov={:.4} strand={} exons={} guide=({})",
                i,
                start,
                end,
                tx.coverage,
                tx.strand,
                tx.exons.len(),
                guide_label
            );
        }
    }
    if !protected.is_empty() {
        protected.extend(txs);
        return (protected, summary);
    }
    (txs, summary)
}

/// Run predcluster filter chain with per-stage ref tracking.
/// Returns the name of the first stage that eliminates the ref match, or "none".
pub fn find_print_predcluster_killing_stage(
    transcripts: Vec<Transcript>,
    config: &RunConfig,
    ref_chain_match: impl Fn(&[Transcript]) -> bool,
) -> &'static str {
    let (_, mut txs): (Vec<_>, Vec<_>) = if config.emit_junction_paths {
        transcripts
            .into_iter()
            .partition(|t| predcluster_protected_source(t.source.as_deref()))
    } else {
        (Vec::new(), transcripts)
    };

    if std::env::var_os("RUSTLE_ENABLE_SPLIT_STRAND_CONFLICT_ENDPOINTS").is_some() {
        txs = split_strand_conflict_endpoints(txs, false);
        if !ref_chain_match(&txs) {
            return "split_strand_conflict_endpoints";
        }
    }
    if config.long_reads {
        txs = eliminate_incomplete_vs_guides(txs, false);
        if !ref_chain_match(&txs) {
            return "eliminate_incomplete_vs_guides";
        }
        txs = retained_intron_filter(txs, config.pairwise_error_perc, false);
        if !ref_chain_match(&txs) {
            return "retained_intron_filter";
        }
    }
    if !config.max_sensitivity {
        txs = pairwise_overlap_filter(
            txs,
            config.pairwise_isofrac,
            config.lowisofrac,
            config.pairwise_error_perc,
            config.pairwise_drop,
            config.singlethr,
            config.long_reads,
            config.long_reads && config.long_read_min_len > 0,
            None,
            false,
        );
        if !ref_chain_match(&txs) {
            return "pairwise_overlap_filter";
        }
        if config.long_reads {
            txs = isofrac(
                txs,
                config.transcript_isofrac,
                config.pairwise_error_perc,
                config.pairwise_drop,
                None,
                false,
                config.transcript_isofrac_keep_min,
            );
            if !ref_chain_match(&txs) {
                return "isofrac";
            }
            // intronic_prediction_filter removed — no C++ equivalent
        }
    }
    if config.long_reads {
        txs = collapse_near_equal_intron_chains(txs, config.junction_correction_window, false);
        if !ref_chain_match(&txs) {
            return "collapse_near_equal_intron_chains";
        }
        txs = dedup_exact_intron_chains(txs, false);
        if !ref_chain_match(&txs) {
            return "dedup_exact_intron_chains";
        }
    }
    if config.dedup_intron_chain_subsets {
        txs = dedup_subset_intron_chains(txs, false);
        if !ref_chain_match(&txs) {
            return "dedup_subset_intron_chains";
        }
    }
    if config.filter_contained {
        txs = filter_contained_transcripts(txs, false);
        if !ref_chain_match(&txs) {
            return "filter_contained_transcripts";
        }
    }
    let mut runoff_dist = if config.long_reads { 0 } else { 200 };
    if config.bundle_merge_dist > runoff_dist {
        runoff_dist = config.bundle_merge_dist;
    }
    // C++ parity: singlethr filter applies to all modes (C++ reference)
    let effective_singlethr = config.singlethr;
    txs = collapse_single_exon_runoff(txs, effective_singlethr, runoff_dist, false);
    if !ref_chain_match(&txs) {
        return "collapse_single_exon_runoff";
    }
    txs = polymerase_runoff_filter(txs, runoff_dist, config.singlethr, None, true, false);
    if !ref_chain_match(&txs) {
        return "polymerase_runoff_filter";
    }
    // C++ parity (C++ reference, 20322): readthr/singlethr gate.
    // - All transcripts: cov >= 1.0 (readthr)
    // - Single-exon transcripts: additional cov >= 4.75 (singlethr)
    // - Mixed mode with guides: use singlethr (4.75) for all
    let mixed_mode = config.long_reads && config.long_read_min_len > 0;
    let has_guides = txs.iter().any(|t| t.ref_transcript_id.is_some());
    let base_threshold = if mixed_mode && has_guides {
        config.singlethr // 4.75 in mixed mode with guides
    } else {
        config.readthr // 1.0 default
    };
    txs.retain(|t| {
        if is_guide_pair(t) {
            return true;
        }
        if is_rescue_protected(t) && (t.hardstart || t.hardend) {
            return true;
        }

        let is_unanchored = !t.hardstart && !t.hardend;
        let effective_base_thr = if is_unanchored && t.exons.len() <= 3 {
            base_threshold * 1.5
        } else {
            base_threshold
        };

        let threshold = if t.exons.len() == 1 {
            config.singlethr
        } else {
            effective_base_thr
        };
        if t.coverage >= threshold {
            return true;
        }
        if t.is_longread && t.longcov >= 2.0 {
            return true;
        }
        readthr_allow_longcov_fallback() && t.is_longread && t.longcov >= threshold
    });
    if !ref_chain_match(&txs) {
        return "readthr_gate";
    }
    "none"
}

/// Remove transcripts that contain any junction (intron boundary) not observed in BAM data.
///
/// C++ parity: C++ reference only assembles paths through junctions present in the split-read index
/// (built during BAM ingestion). Rust's path extension can follow synthetic/future-link graph
/// edges that don't correspond to any observed split read, producing transcripts with novel
/// (zero-read-support) introns. This filter eliminates those false positives.
///
/// `junctions` is the set of BAM-validated junctions as (donor, acceptor) tuples in 0-based
/// half-open coordinates (donor = exon.end exclusive, acceptor = next_exon.start inclusive),
/// matching exactly how `bam.rs` constructs Junction values from CIGAR N operations.
///
/// Single-exon transcripts always pass (no junctions to check).
/// Guide-anchored transcripts (source starts with "guide:") are exempt.
pub fn filter_unsupported_junctions(
    txs: Vec<Transcript>,
    junctions: &HashSet<(u64, u64)>,
    junction_correction_window: u64,
    verbose: bool,
) -> Vec<Transcript> {
    if junctions.is_empty() {
        return txs;
    }
    let junctions_vec: Vec<(u64, u64)> = if junction_correction_window > 0 {
        junctions.iter().copied().collect()
    } else {
        Vec::new()
    };
    let before = txs.len();
    let mut exempted_guide = 0usize;
    let mut exempted_coverage = 0usize;
    let mut filtered_out = 0usize;
    let result: Vec<Transcript> = txs
        .into_iter()
        .filter(|tx| {
            if tx.exons.len() < 2 {
                return true; // single-exon: no junctions
            }
            if is_guide_pair(tx) {
                exempted_guide += 1;
                return true; // guide-anchored: exempt
            }
            // High-coverage transcripts are exempt from junction filtering.
            // These transcripts were successfully extracted by the graph algorithm
            // and have sufficient evidence, even if some junctions were marked as
            // "bad" during graph construction (e.g., due to mm<0 or strand=0).
            // This fixes the disconnect between path extraction and junction filtering
            // that was causing high-quality transcripts to be dropped.
            // Use longcov (abundance) as the threshold since coverage can be very low
            // due to normalization, while longcov reflects actual read support.
            //
            // Path-enumeration sources (`junction_path`, `junction_chain`) intentionally
            // over-emit candidate chains for sensitivity; they must still prove every
            // intron against observed split reads here—do not skip junction checks
            // based on abundance alone (precision gate before final emission).
            let from_path_enumeration = matches!(
                tx.source.as_deref(),
                Some("junction_path") | Some("junction_chain")
            );
            if tx.longcov > 1.0 && !from_path_enumeration {
                exempted_coverage += 1;
                return true;
            }
            let has_support = tx.exons.windows(2).all(|w| {
                let donor = w[0].1;
                let acceptor = w[1].0;
                if junctions.contains(&(donor, acceptor)) {
                    return true;
                }
                if junction_correction_window == 0 {
                    return false;
                }
                // the reference assembler -E parity: for long reads, allow corrected junctions within window.
                junctions_vec.iter().any(|(d2, a2)| {
                    d2.abs_diff(donor) <= junction_correction_window
                        && a2.abs_diff(acceptor) <= junction_correction_window
                })
            });
            if !has_support {
                filtered_out += 1;
            }
            has_support
        })
        .collect();
    if verbose && (result.len() < before || exempted_coverage > 0) {
        eprintln!(
            "    junction_support_filter: kept {}/{} (guide_exempt={}, coverage_exempt={}, filtered={})",
            result.len(), before, exempted_guide, exempted_coverage, filtered_out
        );
    }
    result
}

/// Global cross-strand filter (C++ parity: C++ reference `print_predcluster` cross-strand block).
///
/// In C++ ref, `print_predcluster` processes both strands of a genomic cluster together.
/// When two transcripts on different strands overlap, a single-exon lower-scored transcript
/// is eliminated by a higher-scored transcript on the opposite strand (18787: exons.Count()==1).
///
/// Rustle processes each strand in its own bundle, so the existing `pairwise_overlap_filter`
/// cross-strand branch never fires. This function provides an equivalent global pass over ALL
/// assembled transcripts (from all strands/bundles) after per-bundle filtering is complete.
///
/// Rule implemented (C++ reference first condition):
///   If n1 (higher score) and n2 (lower score) overlap on opposite strands AND n2 has 1 exon
///   AND n2 is not guide-anchored → eliminate n2.
pub fn apply_global_cross_strand_filter(txs: Vec<Transcript>, verbose: bool) -> Vec<Transcript> {
    let n = txs.len();
    if n < 2 {
        return txs;
    }

    // Sort by chromosome and then score descending so we can iterate higher-scored first.
    // We work on index arrays to avoid moving out of txs.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        txs[a].chrom.cmp(&txs[b].chrom).then_with(|| {
            tx_score(&txs[b])
                .partial_cmp(&tx_score(&txs[a]))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });

    let mut dead = SmallBitset::with_capacity(n.min(64));

    for oi in 0..order.len() {
        let n1 = order[oi];
        if dead.contains(n1) {
            continue;
        }
        let t1_start = txs[n1].exons.first().map(|e| e.0).unwrap_or(0);
        let t1_end = txs[n1].exons.last().map(|e| e.1).unwrap_or(0);

        for oj in (oi + 1)..order.len() {
            let n2 = order[oj];
            if dead.contains(n2) {
                continue;
            }
            // Only cross-strand comparisons
            if txs[n1].strand == txs[n2].strand {
                continue;
            }
            // Must be on same chromosome
            if txs[n1].chrom != txs[n2].chrom {
                continue;
            }
            let t2_start = txs[n2].exons.first().map(|e| e.0).unwrap_or(0);
            let t2_end = txs[n2].exons.last().map(|e| e.1).unwrap_or(0);
            // Genomic span overlap check
            if t1_end <= t2_start || t2_end <= t1_start {
                continue;
            }
            // C++ reference first condition: eliminate single-exon lower-scored transcript
            if txs[n2].exons.len() == 1 && !is_guide_pair(&txs[n2]) {
                dead.insert_grow(n2);
            }

            // C++ reference second condition:
            // if n1 is single-exon on reverse strand (tlen<0) and n2 is forward strand (tlen>0)
            // with low coverage (cov < 1/ERROR_PERC = 10), eliminate n1
            let n1_is_rev_long = txs[n1].is_longread && txs[n1].strand == '-';
            let n2_is_fwd_long = !txs[n2].is_longread && txs[n2].strand == '+';
            if n1_is_rev_long && n2_is_fwd_long && txs[n2].exons.len() == 1 && txs[n2].coverage < 10.0 {
                dead.insert_grow(n1);
            }

            // C++ reference: low-coverage antisense elimination
            // If n2 has low coverage across its single exon relative to total coverage, eliminate n2
            // Note: Without bpcov in global pass, we use a simpler heuristic:
            // if n2 single-exon with very low coverage (< singlethr=4.75), eliminate
            if txs[n2].exons.len() == 1 && txs[n2].coverage < 4.75 && !is_guide_pair(&txs[n2]) {
                dead.insert_grow(n2);
            }
        }
    }

    let removed = dead.count_ones();
    if verbose && removed > 0 {
        eprintln!(
            "    cross_strand_filter: removed {} transcript(s) (antisense elimination: single-exon + low-cov)",
            removed
        );
    }

    // Preserve original order
    txs.into_iter()
        .enumerate()
        .filter_map(|(i, t)| if !dead.contains(i) { Some(t) } else { None })
        .collect()
}
