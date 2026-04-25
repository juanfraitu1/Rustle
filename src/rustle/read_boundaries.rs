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

    // RUSTLE_LONGTRIM_CURATED_ONLY=1: skip the bpcov-derivative scan and
    // use only the externally-curated start/end features. Set when the
    // caller has already produced ST-faithful trim points (via
    // find_all_trims_for_bundle) and we don't want this function adding
    // additional bpcov-derivative candidates that over-mark hardstart.
    if std::env::var_os("RUSTLE_LONGTRIM_CURATED_ONLY").is_some() {
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
    _mode: AssemblyMode,
    _long_read_min_len: u64,
    bundle_start: u64,
    bundle_end: u64,
) -> (Vec<ReadBoundary>, Vec<ReadBoundary>) {
    let (mut start_counts, mut end_counts) = collect_raw_boundaries(reads);

    // Long-read mode only: always collect CPAS boundaries

    let mut raw_plus: Vec<(u64, f64)> = Vec::new();
    let mut raw_minus: Vec<(u64, f64)> = Vec::new();

    for r in reads {
        // long-read only: all reads contribute
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

/// Port of StringTie's `find_all_trims` (rlink.cpp:2686-2900) as a
/// pre-filter on `lstart`/`lend`. Walks bpcov over each bundlenode
/// region, accepting only positions where the local coverage drop ratio
/// passes ST's localdrop threshold. CPAS cluster positions force-add
/// regardless of drop (mirroring `lastdrop=0` at rlink.cpp:2712).
///
/// Constants mirror rlink.cpp: CHI_THR=50, CHI_WIN=100, ERROR_PERC=0.1,
/// DROP=0.5, longintronanchor=25, localdrop=ERROR_PERC/DROP=0.2.
pub fn find_all_trims_for_bundle(
    bpcov: &Bpcov,
    bundlenodes_head: Option<&CBundlenode>,
    reads: &[BundleRead],
    bundle_start: u64,
    bundle_end: u64,
    _bundle_strand: char,
) -> (Vec<ReadBoundary>, Vec<ReadBoundary>) {
    // Flatten the linked list of bundlenodes for iteration.
    let mut bundlenodes: Vec<(u64, u64)> = Vec::new();
    let mut cur = bundlenodes_head;
    while let Some(bn) = cur {
        bundlenodes.push((bn.start, bn.end));
        cur = bn.next.as_deref();
    }
    const CHI_WIN: u64 = 100;
    const CHI_THR: u64 = 50;
    const DROP: f64 = 0.5;
    const ERROR_PERC: f64 = 0.1;
    const LONGINTRONANCHOR: u64 = 25;
    // Configurable localdrop threshold (default ST-faithful 0.2 = 5×
    // drop required for source-trim). Higher = more permissive =
    // more trim points emitted. Override via RUSTLE_FAT_LOCALDROP=N.
    let env_localdrop: Option<f64> = std::env::var("RUSTLE_FAT_LOCALDROP")
        .ok().and_then(|v| v.parse().ok());

    // Build CPAS clusters from polyA-bearing read endpoints (same as
    // collect_read_boundaries_with_cpas).
    let mut raw_plus: Vec<(u64, f64)> = Vec::new();
    let mut raw_minus: Vec<(u64, f64)> = Vec::new();
    for r in reads {
        let s_p = match r.exons.first() { Some(e) => e.0, None => continue };
        let e_p = match r.exons.last() { Some(e) => e.1.saturating_sub(1), None => continue };
        let w_plus = r.unaligned_poly_a as f64;
        let w_minus = r.unaligned_poly_t as f64;
        match r.strand {
            '+' => if w_plus > 0.0 { raw_plus.push((e_p, w_plus)); },
            '-' => if w_minus > 0.0 { raw_minus.push((s_p, w_minus)); },
            _ => {
                if w_plus > 0.0 { raw_plus.push((e_p, w_plus)); }
                if w_minus > 0.0 { raw_minus.push((s_p, w_minus)); }
            }
        }
    }
    let cpas_starts_in =
        cluster_positions_with_counts(&raw_minus, CPAS_POS_BIN, CPAS_MIN_SUPPORT);
    let cpas_ends_in =
        cluster_positions_with_counts(&raw_plus, CPAS_POS_BIN, CPAS_MIN_SUPPORT);
    let refend = bundle_end.saturating_sub(1).max(bundle_start);
    let mut cpas_starts: Vec<(u64, f64)> = Vec::new();
    let mut cpas_ends: Vec<(u64, f64)> = Vec::new();
    for (pos, sup) in cpas_starts_in {
        let p = pos.clamp(bundle_start, refend);
        add_cpas_trimpoint(&mut cpas_starts, p, sup, CPAS_POS_BIN);
    }
    for (pos, sup) in cpas_ends_in {
        let p = pos.clamp(bundle_start, refend);
        add_cpas_trimpoint(&mut cpas_ends, p, sup, CPAS_POS_BIN);
    }
    cpas_starts.sort_by_key(|&(p, _)| p);
    cpas_ends.sort_by_key(|&(p, _)| p);

    let cov_avg = |a: u64, b: u64| -> f64 {
        if b < a { return 0.0; }
        let ai = bpcov.idx(a);
        if ai >= bpcov.cov.len() { return 0.0; }
        let bi = bpcov.idx(b).min(bpcov.cov.len() - 1);
        let len = (bi - ai + 1) as f64;
        bpcov.get_cov_range(ai, bi + 1) / len.max(1.0)
    };

    let mut lstart_out: Vec<ReadBoundary> = Vec::new();
    let mut lend_out: Vec<ReadBoundary> = Vec::new();

    fn push_or_replace(
        out: &mut Vec<ReadBoundary>,
        pos: u64,
        cov: f64,
        lastdrop: &mut f64,
        thisdrop: f64,
    ) {
        const CHI_THR_LOCAL: u64 = 50;
        if let Some(last) = out.last() {
            if pos.saturating_sub(last.pos) <= CHI_THR_LOCAL {
                if thisdrop < *lastdrop {
                    let l = out.last_mut().unwrap();
                    l.pos = pos;
                    l.cov = cov;
                    *lastdrop = thisdrop;
                }
                return;
            }
        }
        out.push(ReadBoundary { pos, cov });
        *lastdrop = thisdrop;
    }

    fn match_cpas(i: u64, cpas: &[(u64, f64)]) -> Option<f64> {
        for &(p, sup) in cpas {
            if p.abs_diff(i) <= CPAS_POS_BIN {
                return Some(sup);
            }
        }
        None
    }

    for &(bn_start, bn_end) in &bundlenodes {
        let start = bn_start.max(bundle_start);
        let end = bn_end.min(refend);
        if end <= start { continue; }
        let len = end.saturating_sub(start).saturating_add(1);
        if len < CHI_THR { continue; }

        let local_s: Vec<(u64, f64)> = cpas_starts
            .iter().copied().filter(|&(p, _)| p >= start && p <= end).collect();
        let local_e: Vec<(u64, f64)> = cpas_ends
            .iter().copied().filter(|&(p, _)| p >= start && p <= end).collect();

        let mut localdrop: f64 = env_localdrop.unwrap_or(ERROR_PERC / DROP);
        if env_localdrop.is_none()
            && len < 2 * (CHI_WIN + CHI_THR) + 1
            && len < CHI_WIN
        {
            localdrop = ERROR_PERC / (10.0 * DROP);
        }
        let mut lastdrop_s = localdrop;
        let mut lastdrop_e = localdrop;

        let walk = |i: u64,
                    cov_l: f64,
                    cov_r: f64,
                    lstart_out: &mut Vec<ReadBoundary>,
                    lend_out: &mut Vec<ReadBoundary>,
                    lastdrop_s: &mut f64,
                    lastdrop_e: &mut f64,
                    drop_scale: f64| {
            // CPAS force.
            if let Some(sup) = match_cpas(i, &local_s) {
                push_or_replace(lstart_out, i + 1, sup, lastdrop_s, 0.0);
                return;
            }
            if let Some(sup) = match_cpas(i, &local_e) {
                push_or_replace(lend_out, i, sup, lastdrop_e, 0.0);
                return;
            }
            if cov_l < cov_r {
                let thisdrop = cov_l / cov_r.max(f64::EPSILON);
                if thisdrop < localdrop {
                    push_or_replace(
                        lstart_out, i + 1, (cov_r - cov_l) / drop_scale,
                        lastdrop_s, thisdrop,
                    );
                }
            } else if cov_l > cov_r {
                let thisdrop = cov_r / cov_l.max(f64::EPSILON);
                if thisdrop < localdrop {
                    push_or_replace(
                        lend_out, i, (cov_l - cov_r) / drop_scale,
                        lastdrop_e, thisdrop,
                    );
                }
            }
        };

        // Short region (rlink.cpp:2698-2754).
        if len < 2 * (CHI_WIN + CHI_THR) + 1 {
            let lo = start.saturating_add(LONGINTRONANCHOR);
            let hi = end.saturating_sub(LONGINTRONANCHOR);
            for i in lo..hi {
                let cov_l = cov_avg(start, i.saturating_sub(1));
                let cov_r = cov_avg(i, end);
                walk(i, cov_l, cov_r,
                     &mut lstart_out, &mut lend_out,
                     &mut lastdrop_s, &mut lastdrop_e, DROP);
            }
            continue;
        }

        // Phase 1: ramp-up window (rlink.cpp:2762-2811).
        let winlen = CHI_WIN + CHI_THR;
        for i in (start + CHI_THR - 1)..(start + winlen - 1) {
            let cov_l = cov_avg(start, i);
            let cov_r = cov_avg(i + 1, i + winlen);
            walk(i, cov_l, cov_r,
                 &mut lstart_out, &mut lend_out,
                 &mut lastdrop_s, &mut lastdrop_e, DROP);
        }
        // Phase 2: full sliding window (rlink.cpp:2818-2868).
        let mid_end = end.saturating_sub(winlen);
        for i in (start + winlen - 1)..mid_end {
            let cov_l = cov_avg(i + 1 - winlen, i);
            let cov_r = cov_avg(i + 1, i + winlen);
            walk(i, cov_l, cov_r,
                 &mut lstart_out, &mut lend_out,
                 &mut lastdrop_s, &mut lastdrop_e, DROP * winlen as f64);
        }
        // Phase 3: ramp-down window (rlink.cpp:2872+).
        let phase3_end = end.saturating_sub(CHI_THR + 1);
        for i in mid_end..phase3_end {
            let cov_l = cov_avg(i + 1 - winlen, i);
            let cov_r = cov_avg(i + 1, end);
            walk(i, cov_l, cov_r,
                 &mut lstart_out, &mut lend_out,
                 &mut lastdrop_s, &mut lastdrop_e, DROP);
        }

        // Per-position derivative-magnitude trim — port of
        // collect_longtrim_boundaries_in_span:120-215. Catches sharp
        // single-position coverage transitions that the window-average
        // ratio gates miss. Adds candidates at positions where:
        //   diffval[m] > 1/ERROR_PERC (≥10-read jump) AND
        //   sumstartleft < sumstartright * ERROR_PERC  → lstart trim
        //   diffval[m] < -1/ERROR_PERC AND
        //   sumendleft * ERROR_PERC > sumendright → lend trim
        // Default ON; opt-out via RUSTLE_FAT_DIFFVAL_OFF=1.
        if std::env::var_os("RUSTLE_FAT_DIFFVAL_OFF").is_none()
            && end.saturating_sub(start) > CHI_WIN + CHI_THR
        {
            let mut diffval: Vec<f64> = Vec::with_capacity(CHI_WIN as usize);
            let mut sumstartleft = 0.0f64;
            let mut sumendleft = 0.0f64;
            let mut sumstartright = 0.0f64;
            let mut sumendright = 0.0f64;

            // Use bpcov.cov directly via cov_at (per-position, not avg).
            let cov_at_pos = |p: u64| -> f64 {
                let idx = bpcov.idx(p);
                if idx >= bpcov.cov.len() { 0.0 } else { bpcov.cov[idx] }
            };

            let mut lastcov = cov_at_pos(
                start.saturating_add(LONGINTRONANCHOR).saturating_sub(1),
            );
            // Phase A: prime sumstartleft/sumendleft over first CHI_THR
            for pos in start.saturating_add(LONGINTRONANCHOR)
                ..start.saturating_add(CHI_THR + LONGINTRONANCHOR)
            {
                let icov = cov_at_pos(pos);
                let diff = icov - lastcov;
                diffval.push(diff);
                if diff > 0.0 { sumstartleft += diff; } else { sumendleft -= diff; }
                lastcov = icov;
            }
            // Phase B: prime sumstartright/sumendright over next CHI_WIN
            for pos in start.saturating_add(CHI_THR + LONGINTRONANCHOR)
                ..start.saturating_add(CHI_WIN + LONGINTRONANCHOR)
            {
                let icov = cov_at_pos(pos);
                let diff = icov - lastcov;
                diffval.push(diff);
                if diff > 0.0 { sumstartright += diff; } else { sumendright -= diff; }
                lastcov = icov;
            }
            // Main loop: walk pos through scan range, emit at
            // diffval-magnitude conditions.
            let max_cov_pos = bpcov.bundle_end.saturating_sub(1);
            let covend = end.min(max_cov_pos).saturating_sub(LONGINTRONANCHOR);
            let start_scan = start.saturating_add(CHI_WIN + LONGINTRONANCHOR);
            if covend > start_scan {
                let inv_err = 1.0 / ERROR_PERC; // 10.0
                for pos in start_scan..covend {
                    let m = ((pos
                        .saturating_sub(CHI_THR)
                        .saturating_sub(LONGINTRONANCHOR)
                        .saturating_sub(start))
                        % CHI_WIN) as usize;
                    if m >= diffval.len() { break; }
                    if diffval[m] > inv_err
                        && sumstartleft < sumstartright * ERROR_PERC
                    {
                        let istart = pos.saturating_sub(CHI_THR);
                        push_boundary(&mut lstart_out, istart, diffval[m]);
                    }

                    let icov = cov_at_pos(pos);
                    let p = (m + CHI_THR as usize) % CHI_WIN as usize;
                    if p >= diffval.len() { break; }
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

                    if diffval[m] < -inv_err
                        && sumendleft * ERROR_PERC > sumendright
                    {
                        let iend = pos.saturating_sub(CHI_THR).saturating_sub(1);
                        push_boundary(&mut lend_out, iend, -diffval[m]);
                    }
                }
            }
        }
    }

    lstart_out.sort_by_key(|b| b.pos);
    lend_out.sort_by_key(|b| b.pos);
    if std::env::var_os("RUSTLE_FAT_EMIT_TRACE").is_some() {
        eprintln!(
            "[FAT_EMIT] bundle={}-{} bnodes={} lstart_out={} lend_out={}",
            bundle_start, bundle_end,
            bundlenodes.len(), lstart_out.len(), lend_out.len()
        );
    }
    (lstart_out, lend_out)
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
