//! Trace why reference transcripts are missing: run filter stages and report blocker per ref.

use crate::path_extract::{SeedOutcome, Transcript};
use crate::reference_gtf::RefTranscript;
use crate::transcript_filter::find_print_predcluster_killing_stage;
use crate::types::{Bundle, DetHashMap as HashMap, RunConfig};

/// Default tolerance (bp) when comparing intron chains (donor/acceptor).
const INTRON_TOL_DEFAULT: u64 = 5;

fn intron_tol() -> u64 {
    std::env::var("RUSTLE_TRACE_INTRON_TOL")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(INTRON_TOL_DEFAULT)
}

/// Intron chain from exons (start, end) with end exclusive: introns are (exon_i.end, exon_{i+1}.start).
fn intron_chain(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    for i in 0..exons.len().saturating_sub(1) {
        out.push((exons[i].1, exons[i + 1].0));
    }
    out
}

fn intron_chains_match(a: &[(u64, u64)], b: &[(u64, u64)]) -> bool {
    let tol = intron_tol();
    if a.len() != b.len() {
        return false;
    }
    a.iter()
        .zip(b.iter())
        .all(|((d1, ac1), (d2, ac2))| d1.abs_diff(*d2) <= tol && ac1.abs_diff(*ac2) <= tol)
}

fn intron_eq(a: (u64, u64), b: (u64, u64)) -> bool {
    let tol = intron_tol();
    a.0.abs_diff(b.0) <= tol && a.1.abs_diff(b.1) <= tol
}

fn intron_diff_summary(ref_chain: &[(u64, u64)], tx_chain: &[(u64, u64)]) -> String {
    let mut missing = 0usize;
    for &ri in ref_chain {
        if !tx_chain.iter().copied().any(|ti| intron_eq(ri, ti)) {
            missing += 1;
        }
    }

    let mut extra = 0usize;
    for &ti in tx_chain {
        if !ref_chain.iter().copied().any(|ri| intron_eq(ri, ti)) {
            extra += 1;
        }
    }

    let first_mismatch = ref_chain
        .iter()
        .copied()
        .zip(tx_chain.iter().copied())
        .enumerate()
        .find(|(_, (ri, ti))| !intron_eq(*ri, *ti));

    if let Some((idx, (ri, ti))) = first_mismatch {
        return format!(
            "introns ref={} tx={} missing={} extra={} first_mismatch#{} ref={}-{} tx={}-{}",
            ref_chain.len(),
            tx_chain.len(),
            missing,
            extra,
            idx + 1,
            ri.0,
            ri.1,
            ti.0,
            ti.1
        );
    }

    if ref_chain.len() != tx_chain.len() {
        return format!(
            "introns ref={} tx={} missing={} extra={} prefix_match=true",
            ref_chain.len(),
            tx_chain.len(),
            missing,
            extra
        );
    }

    format!(
        "introns ref={} tx={} missing={} extra={} match=true",
        ref_chain.len(),
        tx_chain.len(),
        missing,
        extra
    )
}

fn exons_exact_match(a: &[(u64, u64)], b: &[(u64, u64)]) -> bool {
    a == b
}

fn exonic_overlap_len(a: &[(u64, u64)], b: &[(u64, u64)]) -> u64 {
    let mut i = 0usize;
    let mut j = 0usize;
    let mut total = 0u64;
    while i < a.len() && j < b.len() {
        let (a0, a1) = a[i];
        let (b0, b1) = b[j];
        let start = a0.max(b0);
        let end = a1.min(b1);
        if start < end {
            total += end - start;
        }
        if a1 <= b1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    total
}

fn find_exact_transcript(ref_tx: &RefTranscript, txs: &[Transcript]) -> Option<usize> {
    for (i, t) in txs.iter().enumerate() {
        if t.chrom != ref_tx.chrom {
            continue;
        }
        if ref_tx.strand != '.' && t.strand != ref_tx.strand {
            continue;
        }
        if exons_exact_match(&ref_tx.exons, &t.exons) {
            return Some(i);
        }
    }
    None
}

fn trace_target_ref_id() -> Option<String> {
    std::env::var("RUSTLE_TRACE_REF_ID")
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
}

fn boundary_distance(ref_tx: &RefTranscript, tx: &Transcript) -> u64 {
    let (Some(ref_first), Some(ref_last), Some(tx_first), Some(tx_last)) = (
        ref_tx.exons.first(),
        ref_tx.exons.last(),
        tx.exons.first(),
        tx.exons.last(),
    ) else {
        return u64::MAX;
    };
    ref_first.0.abs_diff(tx_first.0)
        + ref_first.1.abs_diff(tx_first.1)
        + ref_last.0.abs_diff(tx_last.0)
        + ref_last.1.abs_diff(tx_last.1)
}

fn summarize_exons(exons: &[(u64, u64)]) -> String {
    match (exons.first(), exons.last()) {
        (Some(first), Some(last)) => format!(
            "{} exons first={}-{} last={}-{}",
            exons.len(),
            first.0,
            first.1,
            last.0,
            last.1
        ),
        _ => "0 exons".to_string(),
    }
}

/// Emit focused stage-by-stage debug for a single reference transcript selected by
/// `RUSTLE_TRACE_REF_ID`. This reports whether the exact exon structure is still present
/// and which assembled transcript is currently the closest match.
pub fn debug_ref_stage(stage: &str, context: &str, ref_tx: &RefTranscript, txs: &[Transcript]) {
    let exact_idx = find_exact_transcript(ref_tx, txs);
    let chain_idx = find_matching_transcript(ref_tx, txs);
    eprintln!(
        "[TRACE_REF] ref={} {} stage={} txs={} exact={} chain={} ref={}",
        ref_tx.id,
        context,
        stage,
        txs.len(),
        exact_idx
            .map(|idx| idx.to_string())
            .unwrap_or_else(|| "-".to_string()),
        chain_idx
            .map(|idx| idx.to_string())
            .unwrap_or_else(|| "-".to_string()),
        summarize_exons(&ref_tx.exons)
    );

    let ref_chain = intron_chain(&ref_tx.exons);
    let mut ranked: Vec<(usize, u64, bool, bool, u64)> = txs
        .iter()
        .enumerate()
        .filter(|(_, tx)| {
            tx.chrom == ref_tx.chrom && (ref_tx.strand == '.' || tx.strand == ref_tx.strand)
        })
        .map(|(idx, tx)| {
            (
                idx,
                exonic_overlap_len(&ref_tx.exons, &tx.exons),
                exons_exact_match(&ref_tx.exons, &tx.exons),
                intron_chains_match(&ref_chain, &intron_chain(&tx.exons)),
                boundary_distance(ref_tx, tx),
            )
        })
        .filter(|(_, overlap, _, _, _)| *overlap > 0)
        .collect();
    ranked.sort_by(|a, b| {
        b.2.cmp(&a.2)
            .then_with(|| b.3.cmp(&a.3))
            .then_with(|| b.1.cmp(&a.1))
            .then_with(|| a.4.cmp(&b.4))
            .then_with(|| b.0.cmp(&a.0))
    });

    if ranked.is_empty() {
        eprintln!("[TRACE_REF]   no overlapping transcripts");
        return;
    }

    for (rank, (idx, overlap, exact, chain, boundary)) in ranked.into_iter().take(3).enumerate() {
        let tx = &txs[idx];
        let source = tx.source.as_deref().unwrap_or("-");
        let tx_chain = intron_chain(&tx.exons);
        eprintln!(
            "[TRACE_REF]   cand#{} idx={} exact={} chain={} overlap_bp={} boundary_delta={} cov={:.4} src={} {}",
            rank + 1,
            idx,
            exact,
            chain,
            overlap,
            boundary,
            tx.coverage,
            source,
            summarize_exons(&tx.exons)
        );
        eprintln!(
            "[TRACE_REF]     {}",
            intron_diff_summary(&ref_chain, &tx_chain)
        );
    }
}

pub fn find_traced_ref_in_bundle<'a>(
    bundle: &Bundle,
    ref_transcripts: &'a [RefTranscript],
) -> Option<&'a RefTranscript> {
    let target_id = trace_target_ref_id()?;
    ref_transcripts
        .iter()
        .find(move |r| r.id == target_id && ref_overlaps_bundle(r, bundle))
}

pub fn debug_target_ref_stage(
    stage: &str,
    bundle: &Bundle,
    ref_transcripts: &[RefTranscript],
    txs: &[Transcript],
) {
    let Some(ref_tx) = find_traced_ref_in_bundle(bundle, ref_transcripts) else {
        return;
    };
    let bundle_region = format!("bundle={}:{}-{}", bundle.chrom, bundle.start, bundle.end);
    debug_ref_stage(stage, &bundle_region, ref_tx, txs);
}

/// True if ref transcript overlaps bundle (same chrom, strand, and span overlaps).
pub fn ref_overlaps_bundle(ref_tx: &RefTranscript, bundle: &Bundle) -> bool {
    if ref_tx.chrom != bundle.chrom {
        return false;
    }
    if ref_tx.strand != '.' && ref_tx.strand != bundle.strand {
        return false;
    }
    let (ref_start, ref_end) = match (ref_tx.exons.first(), ref_tx.exons.last()) {
        (Some(&(s, _)), Some(&(_, e))) => (s, e),
        _ => return false,
    };
    ref_start < bundle.end.saturating_add(1) && ref_end > bundle.start
}

/// True if ref transcript is fully contained within bundle span (same chrom/strand).
///
/// This is important for trace mode: bundles are frequently split into sub-bundles/components,
/// and "overlap" alone causes us to compare a full transcript intron chain to a tiny sub-bundle
/// junction set, producing misleading `missing_junctions(0/N)` reports.
fn ref_contained_in_bundle(ref_tx: &RefTranscript, bundle: &Bundle) -> bool {
    if ref_tx.chrom != bundle.chrom {
        return false;
    }
    if ref_tx.strand != '.' && ref_tx.strand != bundle.strand {
        return false;
    }
    let (ref_start, ref_end) = match (ref_tx.exons.first(), ref_tx.exons.last()) {
        (Some(&(s, _)), Some(&(_, e))) => (s, e),
        _ => return false,
    };
    // Bundle end is treated as inclusive throughout the pipeline (trace printing uses end+1).
    // RefTranscript exon end is exclusive. Require the ref end to be <= bundle.end+1.
    ref_start >= bundle.start && ref_end <= bundle.end.saturating_add(1)
}

/// True when `ref_tx` is a single-exon annotation and `t` covers the same outer span (Rustle may
/// split the interval into sub-exons; gffcompare still treats them as one body).
fn single_exon_ref_outer_span_match(ref_tx: &RefTranscript, t: &Transcript, tol: u64) -> bool {
    if ref_tx.exons.len() != 1 || t.exons.is_empty() {
        return false;
    }
    let (r0, r1) = ref_tx.exons[0];
    let q0 = t.exons[0].0;
    let q1 = t.exons.last().map(|e| e.1).unwrap_or(q0);
    q0.abs_diff(r0) <= tol && q1.abs_diff(r1) <= tol
}

/// Find index of first transcript in `txs` that has the same intron chain as `ref_tx`.
pub fn find_matching_transcript(ref_tx: &RefTranscript, txs: &[Transcript]) -> Option<usize> {
    let ref_chain = intron_chain(&ref_tx.exons);
    let boundary_tol = intron_tol().max(10);
    if ref_tx.exons.len() == 1 {
        for (i, t) in txs.iter().enumerate() {
            if t.chrom != ref_tx.chrom {
                continue;
            }
            if ref_tx.strand != '.' && t.strand != ref_tx.strand {
                continue;
            }
            if single_exon_ref_outer_span_match(ref_tx, t, boundary_tol) {
                return Some(i);
            }
        }
        return None;
    }
    for (i, t) in txs.iter().enumerate() {
        if t.chrom != ref_tx.chrom {
            continue;
        }
        if ref_tx.strand != '.' && t.strand != ref_tx.strand {
            continue;
        }
        if intron_chains_match(&ref_chain, &intron_chain(&t.exons)) {
            return Some(i);
        }
    }
    None
}

/// Subcategory for NotExtracted misses.
#[derive(Debug, Clone)]
pub enum NotExtractedReason {
    /// Single-exon ref (no intron chain to check).
    SingleExon,
    /// Some required junctions are missing from the bundle's junction set.
    MissingJunctions { present: usize, required: usize },
    /// All junctions present but path extraction failed (flow/seed issue).
    JunctionsPresent,
    /// Unknown (no junction data available).
    Unknown,
}

/// Blocker reason for a reference transcript.
#[derive(Debug, Clone)]
pub enum Blocker {
    /// No transcript with matching intron chain was extracted (no transfrag or below readthr/singlethr/min_length).
    NotExtracted(NotExtractedReason),
    /// A filter removed it.
    Filter(&'static str),
}

/// Classify NotExtracted subcategory by checking junction availability.
fn classify_not_extracted(
    ref_tx: &RefTranscript,
    bundle: &Bundle,
    bundle_junctions: Option<&[(u64, u64)]>,
) -> NotExtractedReason {
    let tol = intron_tol();
    let ref_chain_all = intron_chain(&ref_tx.exons);
    // Only require introns that fall within the traced bundle span. Trace mode often runs on
    // component-split sub-bundles; without this restriction we can report spurious "missing"
    // junctions that lie completely outside the sub-bundle.
    let ref_chain: Vec<(u64, u64)> = ref_chain_all
        .into_iter()
        .filter(|&(d, a)| {
            // Bundle end is treated as inclusive throughout the pipeline.
            d >= bundle.start && a <= bundle.end.saturating_add(1)
        })
        .collect();
    if ref_chain.is_empty() {
        if ref_tx.exons.len() == 1 {
            return NotExtractedReason::SingleExon;
        }
        // Multi-exon transcript but no introns fall inside this bundle region.
        return NotExtractedReason::Unknown;
    }
    let Some(junctions) = bundle_junctions else {
        return NotExtractedReason::Unknown;
    };
    let required = ref_chain.len();
    let present = ref_chain
        .iter()
        .filter(|&&(d, a)| {
            junctions
                .iter()
                .any(|&(jd, ja)| d.abs_diff(jd) <= tol && a.abs_diff(ja) <= tol)
        })
        .count();
    if present == required {
        NotExtractedReason::JunctionsPresent
    } else {
        // Log which specific junctions are missing to aid diagnosis.
        let missing: Vec<_> = ref_chain
            .iter()
            .filter(|&&(d, a)| {
                !junctions
                    .iter()
                    .any(|&(jd, ja)| d.abs_diff(jd) <= tol && a.abs_diff(ja) <= tol)
            })
            .collect();
        eprintln!(
            "[TRACE_MISSING_JUNC] ref={} present={}/{} missing_coords={:?}",
            ref_tx.id, present, required, missing
        );
        NotExtractedReason::MissingJunctions { present, required }
    }
}

/// Run the same filter chain as the pipeline, recording which stage removes the ref match.
/// `txs` = transcripts after extraction for this bundle. Returns blocker (or Filter("none") if kept).
pub fn trace_ref_in_bundle(
    ref_tx: &RefTranscript,
    bundle: &Bundle,
    txs: Vec<Transcript>,
    bundle_junctions: Option<&[(u64, u64)]>,
    config: &RunConfig,
) -> Blocker {
    if !ref_overlaps_bundle(ref_tx, bundle) {
        return Blocker::NotExtracted(classify_not_extracted(
            ref_tx,
            bundle,
            bundle_junctions,
        ));
    }
    if find_matching_transcript(ref_tx, &txs).is_none() {
        return Blocker::NotExtracted(classify_not_extracted(
            ref_tx,
            bundle,
            bundle_junctions,
        ));
    }
    // Run the actual predcluster filter chain (same as pipeline) with per-stage tracking.
    let ref_tx_clone = ref_tx.clone();
    let killing_stage = find_print_predcluster_killing_stage(txs, config, |txs| {
        find_matching_transcript(&ref_tx_clone, txs).is_some()
    });
    if killing_stage != "none" {
        return Blocker::Filter(killing_stage);
    }
    Blocker::Filter("none")
}

/// (Bundle, transcripts right after extraction, junction set as (donor,acceptor) pairs, seed outcomes) for trace.
pub type BundleRawTx = (
    Bundle,
    Vec<Transcript>,
    Vec<(u64, u64)>,
    Vec<(usize, SeedOutcome)>,
);

/// Aggregate blockers for all refs: ref_id -> list of (bundle_region, blocker).
/// Result of tracing: per-ref blockers + per-ref seed diagnosis strings for junctions_present.
pub struct TraceResult {
    pub by_ref: HashMap<String, Vec<(String, Blocker)>>,
    /// Per ref_id: diagnosis string for junctions_present failures.
    pub seed_diag: HashMap<String, String>,
}

pub fn trace_refs_in_bundles(
    ref_transcripts: &[RefTranscript],
    bundle_raw_txs: &[BundleRawTx],
    config: &RunConfig,
) -> TraceResult {
    let mut by_ref: HashMap<String, Vec<(String, Blocker)>> = Default::default();
    let mut seed_diag: HashMap<String, String> = Default::default();
    for ref_tx in ref_transcripts {
        // Candidate bundles:
        // 1) all bundles that fully contain the ref transcript (preferred),
        // 2) if none, all overlapping bundles.
        let containing_indices: Vec<usize> = bundle_raw_txs
            .iter()
            .enumerate()
            .filter(|(_, (bundle, _, _, _))| ref_contained_in_bundle(ref_tx, bundle))
            .map(|(i, _)| i)
            .collect();
        let candidate_indices: Vec<usize> = if containing_indices.is_empty() {
            bundle_raw_txs
                .iter()
                .enumerate()
                .filter(|(_, (bundle, _, _, _))| ref_overlaps_bundle(ref_tx, bundle))
                .map(|(i, _)| i)
                .collect()
        } else {
            containing_indices
        };

        if candidate_indices.is_empty() {
            continue;
        }

        // Evaluate all candidates. If any candidate keeps the ref ("none"), treat it as resolved.
        let mut has_any_match = false;
        let mut failures: Vec<(usize, String, Blocker, Option<String>)> = Vec::new();
        for i in candidate_indices {
            let (bundle, txs, junctions, seed_outcomes) = &bundle_raw_txs[i];
            let region = format!("{}:{}-{}", bundle.chrom, bundle.start, bundle.end);
            let blocker = trace_ref_in_bundle(ref_tx, bundle, txs.clone(), Some(junctions), config);
            match &blocker {
                Blocker::Filter("none") => {
                    has_any_match = true;
                    break;
                }
                Blocker::NotExtracted(NotExtractedReason::JunctionsPresent) => {
                    let diag = diagnose_junctions_present(ref_tx, junctions, seed_outcomes);
                    failures.push((i, region, blocker, Some(diag)));
                }
                _ => {
                    failures.push((i, region, blocker, None));
                }
            }
        }

        if has_any_match || failures.is_empty() {
            continue;
        }

        // Keep a single primary blocker per ref to avoid duplicate cross-sub-bundle noise.
        let mut best_idx = 0usize;
        for i in 1..failures.len() {
            let (_, _, ref cand_blocker, _) = failures[i];
            let (_, _, ref best_blocker, _) = failures[best_idx];
            let cand_pri = blocker_priority(cand_blocker);
            let best_pri = blocker_priority(best_blocker);
            if cand_pri > best_pri {
                best_idx = i;
                continue;
            }
            if cand_pri == best_pri {
                let cand_span = {
                    let bundle = &bundle_raw_txs[failures[i].0].0;
                    bundle.end.saturating_sub(bundle.start)
                };
                let best_span = {
                    let bundle = &bundle_raw_txs[failures[best_idx].0].0;
                    bundle.end.saturating_sub(bundle.start)
                };
                if cand_span < best_span {
                    best_idx = i;
                }
            }
        }

        let (_bundle_idx, region, blocker, diag_opt) = failures.swap_remove(best_idx);
        if let Some(diag) = diag_opt {
            seed_diag.insert(ref_tx.id.clone(), diag);
        }
        by_ref
            .entry(ref_tx.id.clone())
            .or_default()
            .push((region, blocker));
    }
    TraceResult { by_ref, seed_diag }
}

/// Remove blockers for refs that are present in the final assembled transcript set.
/// This prevents extraction-stage-only misses from being reported when a later bundle
/// or fallback path still emits the matching intron chain.
pub fn drop_resolved_blockers_by_final_matches(
    trace_result: &mut TraceResult,
    ref_transcripts: &[RefTranscript],
    final_txs: &[Transcript],
) {
    if trace_result.by_ref.is_empty() || final_txs.is_empty() {
        return;
    }
    let ref_by_id: HashMap<&str, &RefTranscript> = ref_transcripts
        .iter()
        .map(|r| (r.id.as_str(), r))
        .collect();
    let mut resolved_ids: Vec<String> = Vec::new();
    for ref_id in trace_result.by_ref.keys() {
        if let Some(ref_tx) = ref_by_id.get(ref_id.as_str()) {
            if find_matching_transcript(ref_tx, final_txs).is_some() {
                resolved_ids.push(ref_id.clone());
            }
        }
    }
    for ref_id in resolved_ids {
        trace_result.by_ref.remove(&ref_id);
        trace_result.seed_diag.remove(&ref_id);
    }
}

fn blocker_priority(blocker: &Blocker) -> u8 {
    match blocker {
        Blocker::Filter("none") => 0,
        Blocker::Filter(_) => 5,
        Blocker::NotExtracted(NotExtractedReason::JunctionsPresent) => 4,
        Blocker::NotExtracted(NotExtractedReason::MissingJunctions { .. }) => 3,
        Blocker::NotExtracted(NotExtractedReason::SingleExon) => 2,
        Blocker::NotExtracted(NotExtractedReason::Unknown) => 1,
    }
}

fn blocker_label(blocker: &Blocker) -> String {
    match blocker {
        Blocker::NotExtracted(reason) => match reason {
            NotExtractedReason::SingleExon => "not_extracted:single_exon".to_string(),
            NotExtractedReason::MissingJunctions { present, required } => {
                format!("not_extracted:missing_junctions({}/{})", present, required)
            }
            NotExtractedReason::JunctionsPresent => "not_extracted:junctions_present".to_string(),
            NotExtractedReason::Unknown => "not_extracted:unknown".to_string(),
        },
        Blocker::Filter(n) => n.to_string(),
    }
}

fn blocker_category(blocker: &Blocker) -> &'static str {
    match blocker {
        Blocker::NotExtracted(NotExtractedReason::SingleExon) => "not_extracted:single_exon",
        Blocker::NotExtracted(NotExtractedReason::MissingJunctions { .. }) => {
            "not_extracted:missing_junctions"
        }
        Blocker::NotExtracted(NotExtractedReason::JunctionsPresent) => {
            "not_extracted:junctions_present"
        }
        Blocker::NotExtracted(NotExtractedReason::Unknown) => "not_extracted:unknown",
        Blocker::Filter(n) => n,
    }
}

/// Write blocker report: one line per ref with blocker(s), summary by blocker type, and by strand.
pub fn write_blocker_report(
    by_ref: &HashMap<String, Vec<(String, Blocker)>>,
    seed_diag: &HashMap<String, String>,
    ref_transcripts: &[RefTranscript],
    out: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    writeln!(out, "# reference_id\tblocker\tbundle_region\tseed_diag")?;
    let mut counts: HashMap<String, u64> = Default::default();
    let ref_strand: HashMap<&str, char> = ref_transcripts
        .iter()
        .map(|r| (r.id.as_str(), r.strand))
        .collect();
    let mut blocked_by_strand: HashMap<char, u64> = Default::default();
    let mut ref_ids: Vec<_> = by_ref.keys().collect();
    ref_ids.sort();
    for id in ref_ids {
        let strand = ref_strand.get(id.as_str()).copied().unwrap_or('.');
        *blocked_by_strand.entry(strand).or_insert(0) += 1;
        let diag = seed_diag.get(id.as_str()).map(|s| s.as_str()).unwrap_or("");
        for (region, blocker) in &by_ref[id] {
            let label = blocker_label(blocker);
            let cat = blocker_category(blocker);
            *counts.entry(cat.to_string()).or_insert(0) += 1;
            writeln!(out, "{}\t{}\t{}\t{}", id, label, region, diag)?;
        }
    }
    writeln!(out, "# summary:")?;
    let mut sum: Vec<_> = counts.into_iter().collect();
    sum.sort_by(|a, b| b.1.cmp(&a.1));
    for (name, count) in &sum {
        writeln!(out, "#   {}: {}", name, count)?;
    }
    // Aggregate seed_diag sub-causes for junctions_present refs
    let jp_count = sum
        .iter()
        .find(|(n, _)| n == "not_extracted:junctions_present")
        .map(|(_, c)| *c)
        .unwrap_or(0);
    if jp_count > 0 {
        let mut sub_counts: HashMap<String, u64> = Default::default();
        for (id, blockers) in by_ref {
            let is_jp = blockers.iter().any(|(_, b)| {
                matches!(
                    b,
                    Blocker::NotExtracted(NotExtractedReason::JunctionsPresent)
                )
            });
            if !is_jp {
                continue;
            }
            if let Some(diag) = seed_diag.get(id.as_str()) {
                for token in diag.split_whitespace() {
                    if let Some(eq_pos) = token.find('=') {
                        let key = &token[..eq_pos];
                        if key == "ref_introns" || key == "total_seeds" {
                            continue;
                        }
                        if let Ok(val) = token[eq_pos + 1..].parse::<u64>() {
                            *sub_counts.entry(key.to_string()).or_insert(0) += val;
                        }
                    }
                }
            }
        }
        if !sub_counts.is_empty() {
            let mut sub_sorted: Vec<_> = sub_counts.into_iter().collect();
            sub_sorted.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
            let parts: Vec<String> = sub_sorted
                .iter()
                .map(|(k, v)| format!("{}={}", k, v))
                .collect();
            writeln!(out, "#     seed_diag sub-causes: {}", parts.join(", "))?;
        }
    }
    let total_by_strand: HashMap<char, u64> =
        ref_transcripts
            .iter()
            .fold(Default::default(), |mut m: HashMap<char, u64>, r| {
                *m.entry(r.strand).or_insert(0) += 1;
                m
            });
    writeln!(
        out,
        "# blocked by strand (refs with at least one bundle where not extracted or filtered):"
    )?;
    for strand in ['+', '-'] {
        let blocked = blocked_by_strand.get(&strand).copied().unwrap_or(0);
        let total = total_by_strand.get(&strand).copied().unwrap_or(0);
        let pct = if total > 0 {
            100.0 * blocked as f64 / total as f64
        } else {
            0.0
        };
        writeln!(
            out,
            "#   strand {}: {} blocked / {} total ({:.1}%)",
            strand, blocked, total, pct
        )?;
    }
    Ok(())
}

/// For a junctions_present ref, summarize what happened to seeds whose junctions
/// overlap the ref's intron chain. Returns a diagnostic string.
fn diagnose_junctions_present(
    ref_tx: &RefTranscript,
    _bundle_junctions: &[(u64, u64)],
    seed_outcomes: &[(usize, SeedOutcome)],
) -> String {
    if seed_outcomes.is_empty() {
        return "no_seed_outcomes_collected".to_string();
    }

    // Summarize outcome distribution across ALL seeds (we don't have per-seed junction info
    // at this level — transfrags are consumed during extraction).  The most useful signal is
    // the aggregate outcome distribution: how many seeds hit each failure mode.
    let ref_chain = intron_chain(&ref_tx.exons);
    let n_introns = ref_chain.len();

    let mut counts: HashMap<&'static str, usize> = Default::default();
    for (_, outcome) in seed_outcomes {
        let label = match outcome {
            SeedOutcome::Skipped(reason) => *reason,
            SeedOutcome::BackToSourceFail => "back_to_source_fail",
            SeedOutcome::FwdToSinkFail => "fwd_to_sink_fail",
            SeedOutcome::UnwitnessedSplice => "unwitnessed_splice",
            SeedOutcome::HardBoundaryMismatch => "hard_boundary_mismatch",
            SeedOutcome::ZeroFlux => "zero_flux",
            SeedOutcome::LowCoverage(_) => "low_coverage",
            SeedOutcome::EonlyNonGuide => "eonly_non_guide",
            SeedOutcome::TooShort => "too_short",
            SeedOutcome::ChecktrfReadthr => "checktrf_readthr",
            SeedOutcome::ChecktrfEonlySkip => "checktrf_eonly_skip",
            SeedOutcome::ChecktrfRedistributed => "checktrf_redistributed",
            SeedOutcome::ChecktrfRescued => "checktrf_rescued",
            SeedOutcome::ChecktrfIncomplete => "checktrf_incomplete",
            SeedOutcome::ChecktrfRescueFail => "checktrf_rescue_fail",
            SeedOutcome::Stored(_) => "stored",
        };
        *counts.entry(label).or_insert(0) += 1;
    }

    // Build sorted summary string.
    let mut pairs: Vec<_> = counts.into_iter().collect();
    pairs.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
    let summary: Vec<String> = pairs.iter().map(|(k, v)| format!("{}={}", k, v)).collect();
    format!(
        "ref_introns={} total_seeds={} {}",
        n_introns,
        seed_outcomes.len(),
        summary.join(" ")
    )
}
