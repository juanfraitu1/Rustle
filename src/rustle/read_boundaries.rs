//! Read boundaries for longtrim (lstart/lend; original algorithm collect_read_boundaries).
//! CPred in has predno (position) and cov; used by longtrim to split nodes at read starts/ends.

use crate::bpcov::Bpcov;
use crate::tss_tts::{
    add_cpas_trimpoint, cluster_positions_with_counts, CPAS_MIN_SUPPORT, CPAS_POS_BIN,
};
use crate::types::{AssemblyMode, BundleRead, CBundlenode, CPred, ReadBoundary};

const CHI_THR: u64 = 50;
const CHI_WIN: u64 = 100;
const ERROR_PERC: f64 = 0.1;
const LONGINTRONANCHOR: u64 = 25;

fn collect_raw_boundaries(
    reads: &[BundleRead],
) -> (
    crate::types::DetHashMap<u64, f64>,
    crate::types::DetHashMap<u64, f64>,
) {
    use crate::types::DetHashMap as HashMap;
    let mut start_counts: HashMap<u64, f64> = Default::default();
    let mut end_counts: HashMap<u64, f64> = Default::default();

    for r in reads {
        let w = r.weight;
        if let Some(&s) = r.exons.first().map(|e| &e.0) {
            *start_counts.entry(s).or_insert(0.0) += w;
        }
        if let Some(e) = r.exons.last().map(|x| x.1) {
            let last_base = e.saturating_sub(1);
            *end_counts.entry(last_base).or_insert(0.0) += w;
        }
    }

    (start_counts, end_counts)
}

fn finalize_boundaries(
    start_counts: crate::types::DetHashMap<u64, f64>,
    end_counts: crate::types::DetHashMap<u64, f64>,
) -> (Vec<ReadBoundary>, Vec<ReadBoundary>) {
    let mut lstart: Vec<ReadBoundary> = start_counts
        .into_iter()
        .map(|(pos, cov)| ReadBoundary { pos, cov })
        .collect();
    lstart.sort_by_key(|b| b.pos);

    let mut lend: Vec<ReadBoundary> = end_counts
        .into_iter()
        .map(|(pos, cov)| ReadBoundary { pos, cov })
        .collect();
    lend.sort_by_key(|b| b.pos);

    (lstart, lend)
}

#[inline]
fn is_mixed_long_read(read: &BundleRead, long_read_min_len: u64) -> bool {
    if long_read_min_len == 0 {
        return true;
    }
    read.query_length.map_or(true, |q| q >= long_read_min_len)
}

#[inline]
fn cov_at(bpcov: &Bpcov, pos: u64) -> f64 {
    let idx = bpcov.idx(pos);
    if idx >= bpcov.cov.len() {
        0.0
    } else {
        bpcov.cov[idx]
    }
}

fn push_boundary(boundaries: &mut Vec<ReadBoundary>, pos: u64, cov: f64) {
    if let Some(last) = boundaries.last_mut() {
        if last.pos.abs_diff(pos) < CHI_WIN {
            if last.cov > 0.0 && cov > last.cov {
                last.pos = pos;
                last.cov = cov;
            }
            return;
        }
    }
    boundaries.push(ReadBoundary { pos, cov });
}

/// Direct port of the original algorithm's longtrim boundary scan from `create_graph()`.
///
/// This follows the `diffval`/`sumstartleft`/`sumstartright` logic in
/// `original algorithm_debug/` and returns candidate start/end boundaries for
/// the half-open node span `[span_start, span_end)`.
///
/// `start_features`/`end_features` are optional external hard boundaries. When
/// their `cov` is negative, they behave like the original algorithm's `CPred(..., -1)` point
/// features and suppress nearby diff-based candidates within `CHI_WIN`.
pub fn collect_longtrim_boundaries_in_span(
    bpcov: &Bpcov,
    span_start: u64,
    span_end: u64,
    start_features: &[ReadBoundary],
    end_features: &[ReadBoundary],
) -> (Vec<ReadBoundary>, Vec<ReadBoundary>) {
    let mut lstart: Vec<ReadBoundary> = start_features
        .iter()
        .copied()
        .filter(|b| b.pos >= span_start && b.pos < span_end)
        .collect();
    let mut lend: Vec<ReadBoundary> = end_features
        .iter()
        .copied()
        .filter(|b| b.pos >= span_start && b.pos < span_end)
        .collect();
    lstart.sort_by_key(|b| b.pos);
    lend.sort_by_key(|b| b.pos);

    if span_end <= span_start {
        return (lstart, lend);
    }

    let endbundle = span_end.saturating_sub(1);
    if endbundle.saturating_sub(span_start) <= CHI_WIN + CHI_THR {
        return (lstart, lend);
    }

    let mut diffval: Vec<f64> = Vec::with_capacity(CHI_WIN as usize);
    let mut sumstartleft = 0.0;
    let mut sumendleft = 0.0;
    let mut sumstartright = 0.0;
    let mut sumendright = 0.0;

    let mut lastcov = cov_at(
        bpcov,
        span_start
            .saturating_add(LONGINTRONANCHOR)
            .saturating_sub(1),
    );

    for pos in span_start.saturating_add(LONGINTRONANCHOR)
        ..span_start.saturating_add(CHI_THR + LONGINTRONANCHOR)
    {
        let icov = cov_at(bpcov, pos);
        let diff = icov - lastcov;
        diffval.push(diff);
        if diff > 0.0 {
            sumstartleft += diff;
        } else {
            sumendleft -= diff;
        }
        lastcov = icov;
    }

    for pos in span_start.saturating_add(CHI_THR + LONGINTRONANCHOR)
        ..span_start.saturating_add(CHI_WIN + LONGINTRONANCHOR)
    {
        let icov = cov_at(bpcov, pos);
        let diff = icov - lastcov;
        diffval.push(diff);
        if diff > 0.0 {
            sumstartright += diff;
        } else {
            sumendright -= diff;
        }
        lastcov = icov;
    }

    let max_cov_pos = bpcov.bundle_end.saturating_sub(1);
    let covend = endbundle.min(max_cov_pos).saturating_sub(LONGINTRONANCHOR);
    let start_scan = span_start.saturating_add(CHI_WIN + LONGINTRONANCHOR);
    if covend <= start_scan {
        return (lstart, lend);
    }

    for pos in start_scan..covend {
        let m = ((pos
            .saturating_sub(CHI_THR)
            .saturating_sub(LONGINTRONANCHOR)
            .saturating_sub(span_start))
            % CHI_WIN) as usize;
        if diffval[m] > 1.0 / ERROR_PERC && sumstartleft < sumstartright * ERROR_PERC {
            let istart = pos.saturating_sub(CHI_THR);
            push_boundary(&mut lstart, istart, diffval[m]);
        }

        let icov = cov_at(bpcov, pos);
        let p = (m + CHI_THR as usize) % CHI_WIN as usize;
        let diff = icov - lastcov;
        if diffval[p] > 0.0 {
            sumstartleft -= diffval[p];
        } else {
            sumendleft += diffval[p];
        }
        if diffval[m] > 0.0 {
            sumstartleft += diffval[m];
            sumstartright -= diffval[m];
        } else {
            sumendleft -= diffval[m];
            sumendright += diffval[m];
        }
        diffval[p] = diff;
        if diff > 0.0 {
            sumstartright += diff;
        } else {
            sumendright -= diff;
        }
        lastcov = icov;

        if diffval[m] < -1.0 / ERROR_PERC && sumendleft * ERROR_PERC > sumendright {
            let iend = pos.saturating_sub(CHI_THR).saturating_sub(1);
            push_boundary(&mut lend, iend, -diffval[m]);
        }
    }

    (lstart, lend)
}

pub type LongtrimBoundaryMap =
    crate::types::DetHashMap<usize, (Vec<ReadBoundary>, Vec<ReadBoundary>)>;

/// Build one the original algorithm-style longtrim boundary stream per original bundlenode.
///
/// the original algorithm computes `lstart/lend` once for the full bundlenode span, then reuses
/// that stream across inline `longtrim()` calls while create_graph keeps splitting
/// the same bundlenode. This map preserves that bundlenode-scoped behavior.
pub fn collect_longtrim_boundary_map(
    bpcov: &Bpcov,
    bundlenodes: Option<&CBundlenode>,
    lstart: &[ReadBoundary],
    lend: &[ReadBoundary],
    min_boundary_cov: f64,
) -> LongtrimBoundaryMap {
    let mut out: LongtrimBoundaryMap = Default::default();
    let mut cur = bundlenodes;
    while let Some(bn) = cur {
        // the original algorithm longtrim uses read-boundary peaks (lstart/lend) as candidates, plus optional
        // negative-cov "feature points" that suppress nearby diff-based candidates.
        // Keep this list compact by pre-filtering to points that can actually trigger a split.
        let keep = |b: &ReadBoundary| b.cov < 0.0 || b.cov >= min_boundary_cov;
        let start_features: Vec<ReadBoundary> = lstart
            .iter()
            .copied()
            .filter(keep)
            .filter(|b| b.pos >= bn.start && b.pos < bn.end)
            .collect();
        let end_features: Vec<ReadBoundary> = lend
            .iter()
            .copied()
            .filter(keep)
            .filter(|b| b.pos >= bn.start && b.pos < bn.end)
            .collect();
        out.insert(
            bn.bid,
            collect_longtrim_boundaries_in_span(
                bpcov,
                bn.start,
                bn.end,
                &start_features,
                &end_features,
            ),
        );
        cur = bn.next.as_deref();
    }
    out
}

/// Collect (lstart, lend): sorted lists of (position, coverage) for read starts and ends.
/// Simple read-count approach: sums weighted read start/end positions.
/// Used by short-read coverage_trim (trimnode_all).
pub fn collect_read_boundaries(reads: &[BundleRead]) -> (Vec<ReadBoundary>, Vec<ReadBoundary>) {
    let (start_counts, end_counts) = collect_raw_boundaries(reads);
    finalize_boundaries(start_counts, end_counts)
}

/// Endpoint transport (`CPred`) for pipeline plumbing.
pub fn collect_read_cpreds(reads: &[BundleRead]) -> (Vec<CPred>, Vec<CPred>) {
    let (lstart, lend) = collect_read_boundaries(reads);
    (
        lstart.into_iter().map(CPred::from).collect(),
        lend.into_iter().map(CPred::from).collect(),
    )
}

/// Collect read boundaries and inject clustered long-read CPAS points.
///
/// Parity intent with 
/// - CPAS from long reads is gathered from unaligned poly tail evidence.
/// - plus-strand CPAS contributes to `lend` (sink-side boundaries).
/// - minus-strand CPAS contributes to `lstart` (source-side boundaries).
/// - unstranded (`.`) reads contribute both sides when evidence exists.
pub fn collect_read_boundaries_with_cpas(
    reads: &[BundleRead],
    mode: AssemblyMode,
    long_read_min_len: u64,
    bundle_start: u64,
    bundle_end: u64,
) -> (Vec<ReadBoundary>, Vec<ReadBoundary>) {
    let (mut start_counts, mut end_counts) = collect_raw_boundaries(reads);

    if !matches!(mode, AssemblyMode::LongRead | AssemblyMode::Mixed) {
        return finalize_boundaries(start_counts, end_counts);
    }

    let mut raw_plus: Vec<(u64, f64)> = Vec::new();
    let mut raw_minus: Vec<(u64, f64)> = Vec::new();

    for r in reads {
        if mode == AssemblyMode::Mixed && !is_mixed_long_read(r, long_read_min_len) {
            continue;
        }
        let start = match r.exons.first() {
            Some(e) => e.0,
            None => continue,
        };
        let end = match r.exons.last() {
            Some(e) => e.1.saturating_sub(1),
            None => continue,
        };
        // Match CPAS weighting: use integer unaligned tail evidence counts.
        let w_plus = r.unaligned_poly_a as f64;
        let w_minus = r.unaligned_poly_t as f64;

        match r.strand {
            '+' => {
                if w_plus > 0.0 {
                    raw_plus.push((end, w_plus));
                }
            }
            '-' => {
                if w_minus > 0.0 {
                    raw_minus.push((start, w_minus));
                }
            }
            _ => {
                if w_plus > 0.0 {
                    raw_plus.push((end, w_plus));
                }
                if w_minus > 0.0 {
                    raw_minus.push((start, w_minus));
                }
            }
        }
    }

    let plus_clusters = cluster_positions_with_counts(&raw_plus, CPAS_POS_BIN, CPAS_MIN_SUPPORT);
    let minus_clusters = cluster_positions_with_counts(&raw_minus, CPAS_POS_BIN, CPAS_MIN_SUPPORT);
    if !plus_clusters.is_empty() || !minus_clusters.is_empty() {
        let refend = bundle_end.saturating_sub(1).max(bundle_start);
        let mut cpas_starts: Vec<(u64, f64)> = Vec::new();
        let mut cpas_ends: Vec<(u64, f64)> = Vec::new();

        for (pos, support) in plus_clusters {
            let p = pos.clamp(bundle_start, refend);
            add_cpas_trimpoint(&mut cpas_ends, p, support, CPAS_POS_BIN);
        }
        for (pos, support) in minus_clusters {
            let p = pos.clamp(bundle_start, refend);
            add_cpas_trimpoint(&mut cpas_starts, p, support, CPAS_POS_BIN);
        }

        for (pos, cov) in cpas_starts {
            *start_counts.entry(pos).or_insert(0.0) += cov;
        }
        for (pos, cov) in cpas_ends {
            *end_counts.entry(pos).or_insert(0.0) += cov;
        }
    }

    finalize_boundaries(start_counts, end_counts)
}

/// CPAS-augmented endpoints exported as `CPred` entries.
pub fn collect_read_cpreds_with_cpas(
    reads: &[BundleRead],
    mode: AssemblyMode,
    long_read_min_len: u64,
    bundle_start: u64,
    bundle_end: u64,
) -> (Vec<CPred>, Vec<CPred>) {
    let (lstart, lend) =
        collect_read_boundaries_with_cpas(reads, mode, long_read_min_len, bundle_start, bundle_end);
    (
        lstart.into_iter().map(CPred::from).collect(),
        lend.into_iter().map(CPred::from).collect(),
    )
}
