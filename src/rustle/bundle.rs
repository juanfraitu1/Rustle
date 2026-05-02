//! Bundle detection: group reads by chromosome/strand and junction connectivity.
//! Optional CGroup-based bundlenode construction (union-find grouping).

use anyhow::Result;
use std::path::Path;
use std::sync::Arc;

use crate::assembly_mode::LONGINTRONANCHOR;
use crate::bam::{junctions_from_exons, open_bam, record_ref_span};
use crate::junction_correction::mismatch_anchor;
use crate::types::{
    Bundle, BundleRead, CBundlenode, DetHashMap as HashMap, DetHashSet as HashSet, Junction,
    JunctionStat, JunctionStats, RunConfig,
};

// header: CHI_THR=50, DROP=0.5 (CGroup small-exon logic)
const CGROUP_CHI_THR: u64 = 50;
const CGROUP_DROP: f64 = 0.5;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct PairKey {
    name: Arc<str>,
    lo_start: u64,
    hi_start: u64,
    hi_tag: u32,
}

#[inline]
fn trace_log_style_active() -> bool {
    std::env::var_os("RUSTLE_TRACE_LOG_STYLE").is_some()
}

#[inline]
fn trace_strand_i8(strand: char) -> i8 {
    match strand {
        '-' => -1,
        '+' => 1,
        _ => 0,
    }
}

fn format_trace_exons(exons: &[(u64, u64)]) -> String {
    exons
        .iter()
        .map(|(s, e)| format!("{}-{}", s.saturating_add(1), e))
        .collect::<Vec<_>>()
        .join(" ")
}

fn format_trace_junctions(junctions: &[Junction]) -> String {
    junctions
        .iter()
        .map(|j| format!("{}-{}", j.donor, j.acceptor))
        .collect::<Vec<_>>()
        .join(" ")
}

fn format_deljunc_offsets(exons: &[(u64, u64)], junctions: &[Junction]) -> String {
    junctions
        .iter()
        .enumerate()
        .filter_map(|(i, j)| {
            let left = exons.get(i)?;
            let right = exons.get(i + 1)?;
            let start_off = left.1.saturating_sub(j.donor);
            let end_off = j.acceptor.saturating_sub(right.0);
            Some(format!("{}:({},{})", i, start_off, end_off))
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[derive(Debug, Clone)]
struct DebugBundleTarget {
    chrom: String,
    start: u64,
    end: u64,
}

#[derive(Debug, Default)]
struct EarlyBundleDebug {
    parsed_reads: usize,
    accepted_reads: usize,
    collapsed_reads: usize,
    stored_reads: usize,
    collapse_same_start_candidates: usize,
    collapse_exon_miss: usize,
    collapse_deljunc_miss: usize,
    collapse_exon_miss_logged: usize,
    collapse_deljunc_miss_logged: usize,
    drop_secondary_supp: usize,
    drop_poly_artifact: usize,
    drop_trim_empty: usize,
    pretrim_raw: HashMap<(u64, u64, i8), f64>,
    pretrim_del: HashMap<(u64, u64, i8), f64>,
    posttrim_raw: HashMap<(u64, u64, i8), f64>,
    posttrim_del: HashMap<(u64, u64, i8), f64>,
    stored_raw: HashMap<(u64, u64, i8), f64>,
    stored_del: HashMap<(u64, u64, i8), f64>,
}

fn parse_debug_bundle_target(config: &RunConfig) -> Option<DebugBundleTarget> {
    if !config.debug_junctions {
        return None;
    }
    let spec = config.debug_bundle.as_ref()?;
    let (chrom, rest) = spec.split_once(':')?;
    let (start, end) = rest.split_once('-')?;
    let start = start.replace(',', "").parse::<u64>().ok()?;
    let end = end.replace(',', "").parse::<u64>().ok()?;
    Some(DebugBundleTarget {
        chrom: chrom.to_string(),
        start,
        end,
    })
}

#[inline]
fn intervals_overlap(start_a: u64, end_a: u64, start_b: u64, end_b: u64) -> bool {
    start_a <= end_b && start_b <= end_a
}

#[inline]
fn debug_target_overlaps(
    target: Option<&DebugBundleTarget>,
    chrom: &str,
    start: u64,
    end: u64,
) -> bool {
    target
        .map(|t| t.chrom == chrom && intervals_overlap(start, end, t.start, t.end))
        .unwrap_or(false)
}

fn accumulate_stranded_junction_counts(
    counts: &mut HashMap<(u64, u64, i8), f64>,
    junctions: &[Junction],
    strand: char,
    weight: f64,
) {
    let strand_i8 = trace_strand_i8(strand);
    for j in junctions {
        *counts
            .entry((j.donor, j.acceptor, strand_i8))
            .or_insert(0.0) += weight;
    }
}

fn dump_stranded_junction_counts(label: &str, counts: &HashMap<(u64, u64, i8), f64>) {
    let mut list: Vec<((u64, u64, i8), f64)> = counts.iter().map(|(k, c)| (*k, *c)).collect();
    list.sort_unstable_by_key(|((donor, acceptor, strand), _)| (*donor, *acceptor, *strand));
    eprintln!("[EARLY_DEBUG] {}: {} junctions", label, list.len());
    for ((donor, acceptor, strand), count) in list {
        eprintln!(
            "  {}-{} strand={} count={:.2}",
            donor, acceptor, strand, count
        );
    }
}

fn add_pair_link(reads: &mut [BundleRead], a: usize, b: usize, weight: f64) {
    if a == b || a >= reads.len() || b >= reads.len() || weight <= 0.0 {
        return;
    }

    let mut found = false;
    for i in 0..reads[a].pair_idx.len() {
        if reads[a].pair_idx[i] == b {
            reads[a].pair_count[i] += weight;
            found = true;
            break;
        }
    }
    if !found {
        reads[a].pair_idx.push(b);
        reads[a].pair_count.push(weight);
    }

    found = false;
    for i in 0..reads[b].pair_idx.len() {
        if reads[b].pair_idx[i] == a {
            reads[b].pair_count[i] += weight;
            found = true;
            break;
        }
    }
    if !found {
        reads[b].pair_idx.push(a);
        reads[b].pair_count.push(weight);
    }
}

fn pop_pending_pair(
    pending_by_chrom: &mut HashMap<String, HashMap<PairKey, Vec<(usize, f64)>>>,
    chrom: &str,
    key: &PairKey,
) -> Option<(usize, f64)> {
    let chrom_map = pending_by_chrom.get_mut(chrom)?;
    let vec = chrom_map.get_mut(key)?;
    let v = vec.pop();
    if vec.is_empty() {
        chrom_map.remove(key);
    }
    v
}

fn push_pending_pair(
    pending_by_chrom: &mut HashMap<String, HashMap<PairKey, Vec<(usize, f64)>>>,
    chrom: &str,
    key: PairKey,
    read_idx: usize,
    weight: f64,
) {
    pending_by_chrom
        .entry(chrom.to_string())
        .or_default()
        .entry(key)
        .or_default()
        .push((read_idx, weight));
}

/// One coverage group from union-find grouping (CGroup / STGroup).
#[derive(Debug, Clone)]
struct CGroup {
    start: u64,
    end: u64,
    color: usize,
    grid: usize,
    cov_sum: f64,
    multi: f64,
    next_gr_idx: Option<usize>,
}

fn find_color_root(mut color: usize, eqcol: &[usize]) -> usize {
    while color < eqcol.len() && eqcol[color] != color {
        color = eqcol[color];
    }
    color
}

fn merge_colors(color1: usize, color2: usize, eqcol: &mut [usize]) {
    let mut root1 = color1;
    while eqcol[root1] != root1 {
        root1 = eqcol[root1];
    }
    let mut root2 = color2;
    while eqcol[root2] != root2 {
        root2 = eqcol[root2];
    }
    if root1 < root2 {
        eqcol[root2] = root1;
    } else if root2 < root1 {
        eqcol[root1] = root2;
    }
}

/// Merge group2 into group1 (forward merge along list). merge_fwd_groups.
fn merge_fwd_groups(
    group: &mut [CGroup],
    g1_idx: usize,
    g2_idx: usize,
    merge: &mut [usize],
    eqcol: &mut [usize],
) {
    let c1 = find_color_root(group[g1_idx].color, eqcol);
    let c2 = find_color_root(group[g2_idx].color, eqcol);
    if c1 < c2 {
        eqcol[c2] = c1;
    } else if c2 < c1 {
        eqcol[c1] = c2;
        group[g1_idx].color = c2;
    }
    group[g1_idx].end = group[g2_idx].end;
    group[g1_idx].cov_sum += group[g2_idx].cov_sum;
    group[g1_idx].multi += group[g2_idx].multi;
    group[g1_idx].next_gr_idx = group[g2_idx].next_gr_idx;
    merge[group[g2_idx].grid] = g1_idx;
}

/// Keep exon in CGroup building when groups already exist on this strand
/// (merge_read_to_group, currgroup!=NULL branch).
fn cgroup_keep_exon_existing(
    exon_len: u64,
    exon_idx: usize,
    n_exons: usize,
    junc_supported: &[bool],
    read_len: u64,
    junction_support: u64,
    longreads: bool,
) -> bool {
    let small = exon_len < junction_support
        || (longreads
            && exon_len < CGROUP_CHI_THR
            && (read_len == 0 || (exon_len as f64) < CGROUP_DROP * (read_len as f64)));
    if !small {
        return true;
    }
    if n_exons <= 1 {
        return true;
    }
    if exon_idx == 0 {
        if !junc_supported[0] {
            return false;
        }
        return true;
    }
    if exon_idx < n_exons - 1 {
        if !junc_supported[exon_idx] && !junc_supported[exon_idx - 1] {
            return false;
        }
        return true;
    }
    if n_exons > 1 && !junc_supported[exon_idx - 1] {
        return false;
    }
    true
}

/// Keep exon in CGroup building when no groups yet exist on this strand
/// (merge_read_to_group, currgroup==NULL branch).
fn cgroup_keep_exon_new(
    exon_len: u64,
    exon_idx: usize,
    n_exons: usize,
    junc_supported: &[bool],
    junction_support: u64,
) -> bool {
    if exon_len >= junction_support {
        return true;
    }
    if n_exons <= 1 {
        return true;
    }
    if exon_idx == 0 {
        return n_exons > 1 && junc_supported[0];
    }
    if exon_idx == n_exons - 1 {
        return n_exons > 1 && junc_supported[exon_idx - 1];
    }
    junc_supported[exon_idx - 1] || junc_supported[exon_idx]
}

/// Merge groups within bundledist ( st_merge_close_groups).
///
/// checks boundaryleft/boundaryright before merging
/// If guides are present, merging is prevented across junction boundaries.
fn merge_close_groups(
    group: &mut [CGroup],
    startgroup_idx: Option<usize>,
    bundledist: u64,
    eqcol: &mut [usize],
    merge: &mut [usize],
    boundary_left: &HashSet<u64>,
    boundary_right: &HashSet<u64>,
) {
    let mut lastgroup_idx: Option<usize> = None;
    let mut procgroup_idx = startgroup_idx;
    while let Some(pi) = procgroup_idx {
        if pi >= group.len() {
            break;
        }
        if let Some(li) = lastgroup_idx {
            let gap = group[pi].start.saturating_sub(group[li].end);
            if gap <= bundledist {
                // check boundaryleft/boundaryright
                // Don't merge if there's a boundary at lastgroup.end or procgroup.start.
                let has_boundary = (!boundary_left.is_empty()
                    && boundary_left.contains(&group[li].end))
                    || (!boundary_right.is_empty() && boundary_right.contains(&group[pi].start));
                if !has_boundary {
                    merge_colors(group[li].color, group[pi].color, eqcol);
                    if group[pi].end > group[li].end {
                        group[li].end = group[pi].end;
                    }
                    group[li].cov_sum += group[pi].cov_sum;
                    group[li].next_gr_idx = group[pi].next_gr_idx;
                    merge[group[pi].grid] = li;
                    procgroup_idx = group[pi].next_gr_idx;
                    continue;
                }
            }
        }
        lastgroup_idx = Some(pi);
        procgroup_idx = group[pi].next_gr_idx;
    }
}

#[inline]
fn read_is_long_class(read: &BundleRead, config: &RunConfig) -> bool {
    if !config.long_reads {
        return false;
    }
    if config.long_read_min_len == 0 {
        return true;
    }
    let rlen = read.query_length.unwrap_or_else(|| {
        read.exons
            .iter()
            .map(|(s, e)| e.saturating_sub(*s))
            .sum::<u64>()
    });
    rlen >= config.long_read_min_len
}

/// processRead(701-745): trim terminal polyA/T artifact exons before
/// duplicate matching/pair linking. Returns false when the read should be discarded.
fn apply_processread_terminal_poly_trim(read: &mut BundleRead) -> bool {
    if read.exons.len() < 2 {
        read.has_last_exon_polya = false;
        read.has_first_exon_polyt = false;
        return true;
    }

    let mut rm_last = read.has_last_exon_polya;
    let mut rm_first = read.has_first_exon_polyt;
    if rm_last && rm_first {
        if read.exons.len() == 2 {
            return false;
        }
        if read.strand == '+' {
            rm_first = false;
        } else if read.strand == '-' {
            rm_last = false;
        }
    }

    let mut trimmed = false;
    if rm_last {
        if !read.junctions.is_empty() {
            read.junctions.pop();
        }
        if !read.junctions_raw.is_empty() {
            read.junctions_raw.pop();
        }
        if !read.junctions_del.is_empty() {
            read.junctions_del.pop();
        }
        read.exons.pop();
        read.has_poly_end_aligned = false;
        if !read.has_poly_end_unaligned {
            read.has_poly_end_unaligned = true;
            if read.unaligned_poly_a == 0 {
                read.unaligned_poly_a = 1;
            }
        }
        trimmed = true;
    }

    if rm_first {
        if !read.junctions.is_empty() {
            read.junctions.remove(0);
        }
        if !read.junctions_raw.is_empty() {
            read.junctions_raw.remove(0);
        }
        if !read.junctions_del.is_empty() {
            read.junctions_del.remove(0);
        }
        read.exons.remove(0);
        read.has_poly_start_aligned = false;
        if !read.has_poly_start_unaligned {
            read.has_poly_start_unaligned = true;
            if read.unaligned_poly_t == 0 {
                read.unaligned_poly_t = 1;
            }
        }
        trimmed = true;
    }

    if read.exons.is_empty() {
        return false;
    }

    if trimmed {
        let expected = read.exons.len().saturating_sub(1);
        if read.junctions.len() != expected {
            read.junctions = junctions_from_exons(&read.exons);
        }
        if read.junctions_raw.len() != expected {
            read.junctions_raw = junctions_from_exons(&read.exons);
        }
        if read.junctions_del.len() != expected {
            read.junctions_del = read.junctions.clone();
        }
        if let Some((s, _)) = read.exons.first().copied() {
            read.ref_start = s;
        }
        if let Some((_, e)) = read.exons.last().copied() {
            read.ref_end = e;
        }
        // processRead: infer unknown strand from poly evidence when unambiguous.
        if read.strand == '.' {
            let plus = read.has_poly_end_unaligned || read.has_poly_end_aligned;
            let minus = read.has_poly_start_unaligned || read.has_poly_start_aligned;
            if plus && !minus {
                read.strand = '+';
            } else if minus && !plus {
                read.strand = '-';
            }
        }
    }

    read.has_poly_start = read.has_poly_start_aligned || read.has_poly_start_unaligned;
    read.has_poly_end = read.has_poly_end_aligned || read.has_poly_end_unaligned;
    read.has_last_exon_polya = false;
    read.has_first_exon_polyt = false;
    true
}

/// processRead mismatch gate for one alignment event.
/// Returns rdcount contribution that should be added to junction.nm.
fn processread_mismatch_weight(
    read: &BundleRead,
    this_long: bool,
    current_bundle_start: u64,
    junction_support: u64,
    mismatchfrac: f64,
) -> f64 {
    if read.junctions.is_empty() {
        return 0.0;
    }
    let mut mismatch = false;
    if this_long {
        mismatch = true;
    } else {
        let read_len = read
            .query_length
            .unwrap_or_else(|| read.exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum())
            .max(1);
        if (read.nm as f64) / (read_len as f64) > mismatchfrac {
            mismatch = true;
        } else if read.nm > 0 {
            let first_len = read
                .exons
                .first()
                .map(|(s, e)| e.saturating_sub(*s))
                .unwrap_or(0);
            let last_len = read
                .exons
                .last()
                .map(|(s, e)| e.saturating_sub(*s))
                .unwrap_or(0);
            if (read.clip_left > 0
                && first_len < junction_support.saturating_add(read.clip_left as u64))
                || (read.clip_right > 0
                    && last_len < junction_support.saturating_add(read.clip_right as u64))
                || mismatch_anchor(read, junction_support, current_bundle_start)
            {
                mismatch = true;
            }
        }
    }

    if mismatch || read.nh > 2 {
        read.weight
    } else {
        0.0
    }
}

#[inline]
fn processread_pair_weight(current_rdcount: f64, current_nh: u32, mate_nh: u32) -> f64 {
    if mate_nh > current_nh {
        1.0 / (mate_nh as f64).max(1.0)
    } else {
        current_rdcount
    }
}

#[inline]
fn exonmatch(prev: &[(u64, u64)], exons: &[(u64, u64)]) -> bool {
    prev == exons
}

#[inline]
fn deljuncmatch(prev: &BundleRead, jd: &[Junction]) -> bool {
    if prev.junctions_del.len() != jd.len() {
        return false;
    }
    if jd.is_empty() {
        return true;
    }
    if prev.exons.len() < 2 {
        return false;
    }
    for i in 0..jd.len() {
        let p = prev.junctions_del[i];
        let left = prev.exons[i];
        let right = prev.exons[i + 1];
        let prev_startj = left.1.saturating_sub(p.donor);
        let prev_endj = p.acceptor.saturating_sub(right.0);
        let cur_startj = left.1.saturating_sub(jd[i].donor);
        let cur_endj = jd[i].acceptor.saturating_sub(right.0);
        if prev_startj != cur_startj || prev_endj != cur_endj {
            return false;
        }
    }
    true
}

/// Same CGroup bundling pass, but also returns the exact per-read bundlenode ids implied by the
/// CGroup -> group2bundle flow, so multigraph read mapping can reuse standard readgroup order.
/// Third element: CGroup color root per bundlenode — used for bundle splitting.
pub fn build_bundlenodes_and_readgroups_from_cgroups(
    reads: &[BundleRead],
    good_junctions: &HashSet<Junction>,
    killed_junctions: &HashSet<Junction>,
    junction_support: u64,
    bundledist: u64,
    longreads: bool,
    boundary_left: &HashSet<u64>,
    boundary_right: &HashSet<u64>,
) -> (Option<CBundlenode>, Vec<Vec<usize>>, Vec<usize>) {
    if reads.is_empty() {
        return (None, Vec::new(), Vec::new());
    }

    type ReadSeg = (usize, u64, u64, Vec<(u64, u64)>, Vec<bool>, f64);
    let mut reads_sorted: Vec<ReadSeg> = Vec::new();
    for (read_idx, r) in reads.iter().enumerate() {
        if r.exons.is_empty() {
            continue;
        }
        let start = r.exons[0].0;
        let end = r.exons[r.exons.len() - 1].1;
        let mut junc_supported = Vec::with_capacity(r.exons.len());
        for (i, ex) in r.exons.iter().enumerate() {
            let mut supported = false;
            if i + 1 < r.exons.len() {
                let j = Junction::new(ex.1, r.exons[i + 1].0);
                if good_junctions.contains(&j) {
                    supported = true;
                }
            }
            if i > 0 {
                let j = Junction::new(r.exons[i - 1].1, ex.0);
                if good_junctions.contains(&j) {
                    supported = true;
                }
            }
            junc_supported.push(supported);
        }
        reads_sorted.push((
            read_idx,
            start,
            end,
            r.exons.clone(),
            junc_supported,
            r.weight,
        ));
    }
    // sort by (start, end) to match GSeg sort / coordStableCmp.
    // For longreads mode:readlist.Sort() which uses GSeg comparison (start, then end).
    // For mixed mode: uses coordStableCmp stable sort (start, then end, then pre-sort index).
    // Rust sort_by_key is stable, so equal (start, end) pairs preserve insertion order — matching both.
    reads_sorted.sort_unstable_by_key(|r| (r.1, r.2));

    let mut group: Vec<CGroup> = Vec::new();
    let mut merge: Vec<usize> = Vec::new();
    let mut eqcol: Vec<usize> = Vec::new();
    let mut read_groups: Vec<Vec<usize>> = vec![Vec::new(); reads.len()];
    let mut currgroup_idx: Option<usize> = None;
    let mut startgroup_idx: Option<usize> = None;
    let localdist = if longreads {
        0u64
    } else {
        bundledist.saturating_add(25)
    };

    // RUSTLE_DEBUG_BUNDLE env var enables detailed bundle/bundlenode debug output
    let debug_bundle = std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some();

    for (read_idx, start, end, exons, junc_supported, readcov) in reads_sorted {
        let had_existing_groups = currgroup_idx.is_some();
        let mut read_color = eqcol.len();
        eqcol.push(read_color);
        let read_len = end.saturating_sub(start).saturating_add(1);
        let mut lastgroup_idx: Option<usize> = None;
        let mut thisgroup_idx = currgroup_idx;

        for (ei, &(seg_start, seg_end)) in exons.iter().enumerate() {
            // break color chain on bad junctions (~1477).
            // Color breaks when juncs[i-1]->strand == 0 (junction was killed).
            // We check if the junction is in the killed_junctions set (strand == Some(0)).
            let bad_junction_before = if ei > 0 {
                let j = Junction::new(exons[ei - 1].1, exons[ei].0);
                killed_junctions.contains(&j)
            } else {
                false
            };
            let seg_len = seg_end.saturating_sub(seg_start).saturating_add(1);
            let keep = if had_existing_groups {
                cgroup_keep_exon_existing(
                    seg_len,
                    ei,
                    exons.len(),
                    &junc_supported,
                    read_len,
                    junction_support,
                    longreads,
                )
            } else {
                cgroup_keep_exon_new(seg_len, ei, exons.len(), &junc_supported, junction_support)
            };
            if !keep {
                if had_existing_groups && ei == 0 {
                    if let Some(li) = lastgroup_idx {
                        thisgroup_idx = Some(li);
                    }
                }
                if !had_existing_groups && ei == 0 {
                    // currgroup==NULL branch still creates the first group even when the
                    // first exon is not kept, but does not attach it to readgroup.
                } else {
                    continue;
                }
            }

            while let Some(ti) = thisgroup_idx {
                if ti >= group.len() {
                    break;
                }
                if seg_start > group[ti].end {
                    lastgroup_idx = Some(ti);
                    thisgroup_idx = group[ti].next_gr_idx;
                } else {
                    break;
                }
            }

            let overlaps = thisgroup_idx.map_or(false, |ti| {
                ti < group.len() && seg_end.saturating_add(localdist) >= group[ti].start
            });

            if overlaps {
                let ti = thisgroup_idx.unwrap();
                if read_groups[read_idx].last() != Some(&ti) {
                    read_groups[read_idx].push(ti);
                }
                if seg_start < group[ti].start {
                    group[ti].start = seg_start;
                }
                let mut next_idx = group[ti].next_gr_idx;
                while let Some(ni) = next_idx {
                    if ni >= group.len() {
                        break;
                    }
                    if seg_end >= group[ni].start {
                        merge_fwd_groups(&mut group, ti, ni, &mut merge, &mut eqcol);
                        next_idx = group[ti].next_gr_idx;
                    } else {
                        break;
                    }
                }
                if seg_end > group[ti].end {
                    group[ti].end = seg_end;
                }
                group[ti].cov_sum += (seg_len as f64) * readcov;
                group[ti].multi += (seg_len as f64) * readcov;

                // on bad junction, adopt group's color WITHOUT union-find merge
                // (~1428-1431: readcol = thisgroup->color). This keeps colors
                // on opposite sides of a bad junction separate → different bundles.
                if bad_junction_before {
                    let group_root = find_color_root(group[ti].color, &eqcol);
                    if debug_bundle {
                        eprintln!(
                            "  COLOR_ADOPT read={} exon={} group_root={} old_read_color={}",
                            read_idx, ei, group_root, read_color
                        );
                    }
                    read_color = group_root;
                } else {
                    let this_root = find_color_root(group[ti].color, &eqcol);
                    let read_root = find_color_root(read_color, &eqcol);
                    if this_root != read_root {
                        if debug_bundle {
                            eprintln!(
                                "  COLOR_MERGE read={} exon={} this_root={} read_root={} -> merge",
                                read_idx, ei, this_root, read_root
                            );
                        }
                        if read_root < this_root {
                            eqcol[this_root] = read_root;
                            group[ti].color = read_root;
                        } else {
                            eqcol[read_root] = this_root;
                        }
                    }
                }
            } else {
                // (merge_read_to_group):
                // create a new color on bad-junction boundaries only when this exon
                // starts a new group (non-overlap path), not on every bad junction.
                if bad_junction_before {
                    let new_color = eqcol.len();
                    if debug_bundle {
                        eprintln!(
                            "  COLOR_BREAK_NEWGROUP read={} exon={} junction={}-{} new_color={}",
                            read_idx,
                            ei,
                            exons[ei - 1].1,
                            exons[ei].0,
                            new_color
                        );
                    }
                    eqcol.push(new_color);
                    read_color = new_color;
                }
                let ngroup = group.len();
                let next_gr_idx = thisgroup_idx;
                group.push(CGroup {
                    start: seg_start,
                    end: seg_end,
                    color: read_color, // uses readcol
                    grid: ngroup,
                    cov_sum: (seg_len as f64) * readcov,
                    multi: (seg_len as f64) * readcov,
                    next_gr_idx,
                });
                merge.push(ngroup);
                if let Some(li) = lastgroup_idx {
                    group[li].next_gr_idx = Some(ngroup);
                } else {
                    currgroup_idx = Some(ngroup);
                    if startgroup_idx.is_none() {
                        startgroup_idx = Some(ngroup);
                    }
                }
                if keep {
                    read_groups[read_idx].push(ngroup);
                }
                lastgroup_idx = Some(ngroup);
                thisgroup_idx = Some(ngroup);
            }
        }
    }

    if bundledist > 0 {
        merge_close_groups(
            &mut group,
            startgroup_idx,
            bundledist,
            &mut eqcol,
            &mut merge,
            boundary_left,
            boundary_right,
        );
    }

    fn canonical_grid(g: usize, merge: &[usize]) -> usize {
        let mut g = g;
        while merge[g] != g {
            g = merge[g];
        }
        g
    }

    let mut all_indices: Vec<usize> = Vec::new();
    let mut idx = startgroup_idx;
    while let Some(i) = idx {
        if i < group.len() {
            all_indices.push(i);
            idx = group[i].next_gr_idx;
        } else {
            break;
        }
    }

    fn find_color_root_eqcol(c: usize, eqcol: &[usize]) -> usize {
        let mut r = c;
        while r < eqcol.len() && eqcol[r] != r {
            r = eqcol[r];
        }
        r
    }

    let mut canonical: Vec<(usize, u64, u64, f64, usize)> = Vec::new();
    for &i in &all_indices {
        if canonical_grid(i, &merge) == i {
            let col_root = find_color_root_eqcol(group[i].color, &eqcol);
            canonical.push((i, group[i].start, group[i].end, group[i].cov_sum, col_root));
        }
    }
    canonical.sort_unstable_by_key(|c| c.1);

    let localdist_bnode = bundledist;
    // RUSTLE_DEBUG_BUNDLE: detailed CGroup/bundlenode state for comparison
    let debug_detail = std::env::var_os("RUSTLE_DEBUG_DETAIL").is_some();
    if debug_detail || debug_bundle {
        let mut color_set = crate::bitset::SmallBitset::with_capacity(canonical.len().min(64));
        for &(_, _, _, _, col) in &canonical {
            color_set.insert_grow(col);
        }
        eprintln!(
            "DEBUG_CGROUP canonical={} colors={} localdist_bnode={}",
            canonical.len(),
            color_set.count_ones(),
            localdist_bnode
        );
        for &(root, start, end, cov, col) in &canonical {
            eprintln!(
                "  CGROUP root={} {}-{} cov={:.1} color={}",
                root, start, end, cov, col
            );
        }
        // Emit DEBUG_BG_GRP lines matching format for cross-comparison
        eprintln!("DEBUG_BG_GROUPS ngroups={}", canonical.len());
        eprintln!("DEBUG_BG_DISTINCT_COLORS total={}", color_set.count_ones());
        for &(root, start, end, cov, col) in &canonical {
            eprintln!(
                "DEBUG_BG_GRP grid={} start={} end={} color={} cov={:.1}",
                root,
                start + 1,
                end,
                col,
                cov
            );
        }
        // Print eqcol state for debugging color union-find
        eprintln!("DEBUG_EQCOL len={}", eqcol.len());
        for (i, &root) in eqcol.iter().enumerate() {
            if i != root {
                eprintln!("  EQCOL color={} -> root={}", i, root);
            }
        }
    }
    let mut bnode_list: Vec<CBundlenode> = Vec::new();
    let mut bnode_colors: Vec<usize> = Vec::new(); // CGroup color root per bundlenode
    let mut root_to_bid: Vec<usize> = vec![usize::MAX; group.len()];
    let mut last_by_color: HashMap<usize, usize> = Default::default();
    for (root, start, end, cov_sum, col_root) in canonical {
        if let Some(&i) = last_by_color.get(&col_root) {
            if start <= bnode_list[i].end.saturating_add(localdist_bnode) {
                bnode_list[i].end = bnode_list[i].end.max(end);
                bnode_list[i].cov += cov_sum;
                root_to_bid[root] = i;
                continue;
            }
        }
        let bid = bnode_list.len();
        bnode_list.push(CBundlenode::new(start, end, cov_sum, bid));
        bnode_colors.push(col_root);
        root_to_bid[root] = bid;
        last_by_color.insert(col_root, bid);
    }
    if bnode_list.is_empty() {
        return (None, read_groups, Vec::new());
    }

    let mut read_bnodes: Vec<Vec<usize>> = vec![Vec::new(); reads.len()];
    for (read_idx, groups_for_read) in read_groups.iter().enumerate() {
        for &g in groups_for_read {
            if g >= merge.len() {
                continue;
            }
            let root = canonical_grid(g, &merge);
            if root >= root_to_bid.len() {
                continue;
            }
            let bid = root_to_bid[root];
            if bid == usize::MAX {
                continue;
            }
            if read_bnodes[read_idx].last() != Some(&bid) {
                read_bnodes[read_idx].push(bid);
            }
        }
    }

    // RUSTLE_DEBUG_BUNDLE: detailed bundlenode output
    if debug_bundle {
        eprintln!("DEBUG_BNODE_LIST count={}", bnode_list.len());
        for (i, bn) in bnode_list.iter().enumerate() {
            let color = bnode_colors.get(i).copied().unwrap_or(0);
            eprintln!(
                "  BNODE bid={} start={} end={} cov={:.1} color={}",
                i,
                bn.start + 1,
                bn.end,
                bn.cov,
                color
            );
        }
        // Show color -> bundlenode mapping
        let mut color_to_bnodes: HashMap<usize, Vec<usize>> = Default::default();
        for (i, &color) in bnode_colors.iter().enumerate() {
            color_to_bnodes.entry(color).or_default().push(i);
        }
        eprintln!("DEBUG_COLOR_MAP unique_colors={}", color_to_bnodes.len());
        for (color, bnodes) in color_to_bnodes.iter() {
            eprintln!("  COLOR {} -> bnodes {:?}", color, bnodes);
        }
    }

    let mut head = bnode_list[0].clone();
    let mut curr = &mut head;
    for (i, bn) in bnode_list.iter().enumerate().skip(1) {
        let mut next = bn.clone();
        next.bid = i;
        curr.next = Some(Box::new(next));
        curr = curr.next.as_deref_mut().unwrap();
    }
    (Some(head), read_bnodes, bnode_colors)
}

/// 3-strand CGroup/bundle construction using group-color flow
/// (`add_read_to_group`, `set_strandcol`, `eqnegcol/eqposcol`).
///
/// Returns bundlenodes for `target_strand` only, but computes color propagation using
/// all `region_reads` strands together.
pub fn build_bundlenodes_and_readgroups_from_cgroups_3strand(
    region_reads: &[BundleRead],
    target_strand: char,
    good_junctions: &HashSet<Junction>,
    killed_junctions: &HashSet<Junction>,
    junction_support: u64,
    bundledist: u64,
    longreads: bool,
    _boundary_left: &HashSet<u64>,
    _boundary_right: &HashSet<u64>,
    _junction_stats: Option<&crate::types::JunctionStats>,
    _junction_redirect_map: &crate::types::DetHashMap<Junction, Junction>,
) -> (Option<CBundlenode>, Vec<Vec<usize>>, Vec<usize>, Vec<f64>) {
    if region_reads.is_empty() {
        return (None, Vec::new(), Vec::new(), Vec::new());
    }

    let mut cfg = RunConfig::default();
    cfg.long_reads = longreads;
    cfg.long_read_min_len = 0;
    cfg.bundle_merge_dist = bundledist;
    cfg.junction_support = junction_support;

    let mut sub_bundles = match crate::bundle_builder::build_sub_bundles(
        region_reads,
        &cfg,
        good_junctions,
        killed_junctions,
        None, None, 0,
    ) {
        Ok(v) => v,
        Err(e) => {
            if std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() {
                eprintln!(
                    "DEBUG_3STRAND_ERROR strand={} reads={} err={}",
                    target_strand,
                    region_reads.len(),
                    e
                );
            }
            return (
                None,
                vec![Vec::new(); region_reads.len()],
                Vec::new(),
                vec![0.0; region_reads.len()],
            );
        }
    };

    // for stranded targets, include the unstranded (sno=1) lane so
    // rprop-weighted neutral evidence can contribute to +/- graph construction.
    // Keep an explicit opt-out for A/B diagnostics.
    let include_dot_for_stranded =
        target_strand != '.' && std::env::var_os("RUSTLE_DISABLE_3STRAND_INCLUDE_DOT").is_none();
    // Diagnostic debug mode: when processing an unstranded ('.') bundle in Rust's
    // per-bundle pipeline, include projected +/- strand bundle outputs from the
    // 3-strand sweep as well. processes all three bundle arrays together,
    // while Rust currently invokes this per target bundle.
    let include_projected_for_dot =
        target_strand == '.' && std::env::var_os("RUSTLE_3STRAND_DOT_INCLUDE_STRANDED").is_some();
    let diag_strand_count = std::env::var_os("RUSTLE_3STRAND_STRANDCOUNT").is_some();
    let mut n_neg = 0usize;
    let mut n_dot = 0usize;
    let mut n_pos = 0usize;
    if diag_strand_count {
        for b in &sub_bundles {
            match b.strand {
                '-' => n_neg += 1,
                '+' => n_pos += 1,
                _ => n_dot += 1,
            }
        }
    }
    // Pick SubBundleResults for this pipeline bundle's strand. For unstranded ('.')
    // bundles, `build_sub_bundles` may emit only stranded ('-' / '+') TB rows (common on
    // sparse long-read loci). The strict `strand == '.'` filter would then drop everything,
    // leaving no bundlenodes / empty partition geometry vs StringTie. When the primary
    // filter matches nothing, fall back to all stranded sub-bundles (same effect as
    // `RUSTLE_3STRAND_DOT_INCLUDE_STRANDED=1`, but automatic).
    let primary_indices: Vec<usize> = sub_bundles
        .iter()
        .enumerate()
        .filter(|(_, b)| {
            if include_projected_for_dot {
                return true;
            }
            if include_dot_for_stranded {
                return b.strand == target_strand || b.strand == '.';
            }
            b.strand == target_strand
        })
        .map(|(i, _)| i)
        .collect();
    let chosen_indices: Vec<usize> = if primary_indices.is_empty() && target_strand == '.' {
        sub_bundles
            .iter()
            .enumerate()
            .filter(|(_, b)| b.strand == '-' || b.strand == '+')
            .map(|(i, _)| i)
            .collect()
    } else {
        primary_indices
    };
    let mut chosen_sorted = chosen_indices;
    chosen_sorted.sort_unstable();
    let mut strand_bundles: Vec<crate::bundle_builder::SubBundleResult> = Vec::new();
    for &i in chosen_sorted.iter().rev() {
        strand_bundles.push(sub_bundles.swap_remove(i));
    }
    strand_bundles.reverse();
    if diag_strand_count {
        let mut n_sel_neg = 0usize;
        let mut n_sel_dot = 0usize;
        let mut n_sel_pos = 0usize;
        for b in &strand_bundles {
            match b.strand {
                '-' => n_sel_neg += 1,
                '+' => n_sel_pos += 1,
                _ => n_sel_dot += 1,
            }
        }
        eprintln!(
            "DEBUG_3STRAND_COUNTS target={} include_dot_for_stranded={} include_dot_proj={} all(-/.+/)={}/{}/{} selected(-/.+/)={}/{}/{} region_reads={}",
            target_strand,
            include_dot_for_stranded,
            include_projected_for_dot,
            n_neg,
            n_dot,
            n_pos,
            n_sel_neg,
            n_sel_dot,
            n_sel_pos,
            region_reads.len()
        );
    }
    if strand_bundles.is_empty() {
        return (
            None,
            vec![Vec::new(); region_reads.len()],
            Vec::new(),
            vec![0.0; region_reads.len()],
        );
    }
    strand_bundles.sort_unstable_by_key(|b| b.start);

    let mut merged_nodes_raw: Vec<(u64, u64, f64, usize)> = Vec::new(); // (start,end,cov,old_bid)
    let mut merged_colors_raw: Vec<usize> = Vec::new(); // one color root per merged node
    let mut read_bnodes_old: Vec<Vec<usize>> = vec![Vec::new(); region_reads.len()];
    let mut read_scales: Vec<f64> = vec![0.0; region_reads.len()];

    // Preserve per-bundlenode cgroup colors from each `SubBundleResult` (see `bnode_colors`).
    // Pack (sub_bundle_index, local_color) so disjoint sub-bundles never share a root id.
    const COLOR_SPACE: usize = 1_000_000_000;
    for (sub_bundle_idx, bres) in strand_bundles.into_iter().enumerate() {
        let mut local_nodes: Vec<(u64, u64, f64)> = Vec::new();
        let mut cur = bres.bnode_head.as_ref();
        while let Some(n) = cur {
            local_nodes.push((n.start, n.end, n.cov));
            cur = n.next.as_deref();
        }
        let bid_offset = merged_nodes_raw.len();
        for (li, (s, e, c)) in local_nodes.into_iter().enumerate() {
            merged_nodes_raw.push((s, e, c, bid_offset + li));
            let local_col = bres
                .bnode_colors
                .get(li)
                .copied()
                .unwrap_or(sub_bundle_idx);
            let global_col = sub_bundle_idx
                .saturating_mul(COLOR_SPACE)
                .saturating_add(local_col);
            merged_colors_raw.push(global_col);
        }
        for (ri, bids) in bres.read_to_bnodes.iter().enumerate() {
            if ri >= read_bnodes_old.len() {
                break;
            }
            for &lbid in bids {
                let gbid = bid_offset.saturating_add(lbid);
                read_bnodes_old[ri].push(gbid);
            }
        }
        for (ri, &scale) in bres.read_scale.iter().enumerate() {
            if ri < read_scales.len() && scale > read_scales[ri] {
                read_scales[ri] = scale;
            }
        }
    }

    // Deduplicate per-read old bnode assignments in one pass (faster than
    // repeated Vec::contains checks during accumulation).
    for bids in &mut read_bnodes_old {
        bids.sort_unstable();
        bids.dedup();
    }

    if merged_nodes_raw.is_empty() {
        return (
            None,
            vec![Vec::new(); region_reads.len()],
            Vec::new(),
            read_scales,
        );
    }

    let mut order: Vec<usize> = (0..merged_nodes_raw.len()).collect();
    order.sort_unstable_by_key(|&i| (merged_nodes_raw[i].0, merged_nodes_raw[i].1));
    let mut old_to_new: Vec<usize> = vec![0; merged_nodes_raw.len()];
    let mut sorted_nodes: Vec<(u64, u64, f64)> = Vec::with_capacity(merged_nodes_raw.len());
    let mut sorted_colors: Vec<usize> = Vec::with_capacity(merged_nodes_raw.len());
    for (new_bid, old_bid) in order.into_iter().enumerate() {
        old_to_new[old_bid] = new_bid;
        sorted_nodes.push((
            merged_nodes_raw[old_bid].0,
            merged_nodes_raw[old_bid].1,
            merged_nodes_raw[old_bid].2,
        ));
        sorted_colors.push(merged_colors_raw[old_bid]);
    }

    let mut read_bnodes_new: Vec<Vec<usize>> = vec![Vec::new(); read_bnodes_old.len()];
    for (ri, bids) in read_bnodes_old.iter().enumerate() {
        for &old_bid in bids {
            if old_bid >= old_to_new.len() {
                continue;
            }
            let nb = old_to_new[old_bid];
            read_bnodes_new[ri].push(nb);
        }
        read_bnodes_new[ri].sort_unstable();
        read_bnodes_new[ri].dedup();
    }

    let mut head = CBundlenode::new(sorted_nodes[0].0, sorted_nodes[0].1, sorted_nodes[0].2, 0);
    let mut cur = &mut head;
    for (i, &(s, e, c)) in sorted_nodes.iter().enumerate().skip(1) {
        cur.next = Some(Box::new(CBundlenode::new(s, e, c, i)));
        cur = cur.next.as_deref_mut().unwrap();
    }

    (Some(head), read_bnodes_new, sorted_colors, read_scales)
}

/// Detect bundles from BAM and optionally use a reference FASTA for SNP extraction.
pub fn detect_bundles_from_bam_with_snp<P: AsRef<Path>>(
    bam_path: P,
    config: &RunConfig,
    chrom_filter: Option<&str>,
    genome: Option<&crate::genome::GenomeIndex>,
) -> Result<Vec<Bundle>> {
    let min_intron = config.min_intron_length;
    let runoffdist = config.bundle_runoff_dist;
    let trace_log_style = trace_log_style_active();

    let mut reader = open_bam(bam_path, config.threads)?;
    let header = reader.read_header()?;
    let mut record = noodles_sam::alignment::RecordBuf::default();

    // Boundary detection: all alignments (incl. supp/sec). Regions 0-based (start_incl, end_incl).
    let mut regions: Vec<(String, u64, u64)> = Vec::new();
    let mut current_chrom: Option<String> = None;
    let mut current_start: u64 = 0;
    let mut current_end_incl: u64 = 0;
    let mut has_current = false;

    // readlist ingestion should include all mapped alignments by default.
    // Keep a diagnostic env switch to force primary-only behavior for A/B tests.
    let mut reads_by_chrom: HashMap<String, Vec<BundleRead>> = Default::default();
    let mut pair_pending_by_chrom: HashMap<String, HashMap<PairKey, Vec<(usize, f64)>>> =
        Default::default();
    let boundary_primary_only = std::env::var_os("RUSTLE_BOUNDARY_PRIMARY_ONLY").is_some();
    let reads_primary_only = std::env::var_os("RUSTLE_READS_PRIMARY_ONLY").is_some();
    let stream_coupled_boundary = std::env::var_os("RUSTLE_STREAM_COUPLED_BOUNDARY").is_some();
    let debug_target = parse_debug_bundle_target(config);
    let mut early_debug = EarlyBundleDebug::default();

    while reader.read_record_buf(&header, &mut record)? > 0 {
        let name = record
            .reference_sequence(&header)
            .transpose()
            .ok()
            .flatten()
            .map(|(name_ref, _)| String::from_utf8_lossy(name_ref).into_owned())
            .unwrap_or_else(|| "?".to_string());

        if let Some(cf) = chrom_filter {
            if name != cf {
                continue;
            }
        }

        // LOO experiment: skip reads whose alignment overlaps a mask region.
        // The same read's sequence is collected by the HMM rescue path
        // (vg_hmm/rescue.rs) so it can be recovered as a "novel" copy.
        if !config.vg_mask_regions.is_empty() {
            if let Some((rs, re_excl)) = crate::bam::record_ref_span(&record) {
                let masked = config.vg_mask_regions.iter().any(|(c, ms, me)| {
                    c == &name && rs < *me && *ms < re_excl
                });
                if masked {
                    continue;
                }
            }
        }

        let is_primary = !record.flags().is_secondary() && !record.flags().is_supplementary();
        let mut parsed_read = crate::bam::record_to_bundle_read_with_snp(
            &record,
            Some(&name),
            genome,
        );

        let mut boundary_span: Option<(u64, u64)> = None;
        if stream_coupled_boundary {
            // Experimental stream-coupled boundary mode:
            // update bundle regions from the same parsed alignment stream used by ingest,
            // with processRead's long-read sec/supp gate applied before boundary growth.
            if let Some(read) = parsed_read.as_ref() {
                let this_long = read_is_long_class(read, config);
                let is_secondary_or_supp =
                    record.flags().is_secondary() || record.flags().is_supplementary();
                let drop_sec_supp_for_mode = is_secondary_or_supp && config.long_reads && this_long;
                if (!boundary_primary_only || is_primary) && !drop_sec_supp_for_mode {
                    boundary_span = Some((read.ref_start, read.ref_end.saturating_sub(1)));
                }
            } else if !boundary_primary_only || is_primary {
                if let Some((ref_start, ref_end_excl)) = record_ref_span(&record) {
                    boundary_span = Some((ref_start, ref_end_excl.saturating_sub(1)));
                }
            }
        } else if !boundary_primary_only || is_primary {
            // Default mode: preserve historical region construction from full alignment stream.
            let huge_intron_cap: Option<u64> = std::env::var("RUSTLE_BOUNDARY_SPLIT_HUGE_INTRON")
                .ok()
                .and_then(|v| v.parse().ok())
                .filter(|&v| v > 0);
            let span_opt = if let Some(cap) = huge_intron_cap {
                crate::bam::record_ref_span_capped(&record, cap)
            } else {
                record_ref_span(&record)
            };
            if let Some((ref_start, ref_end_excl)) = span_opt {
                boundary_span = Some((ref_start, ref_end_excl.saturating_sub(1)));
            }
        }
        if let Some((ref_start, end_incl)) = boundary_span {
            let new_chrom = current_chrom.as_ref().map(|c| c.as_str()) != Some(name.as_str());
            let runoff_gap = has_current && ref_start > current_end_incl.saturating_add(runoffdist);

            if new_chrom {
                if has_current {
                    regions.push((
                        current_chrom.take().unwrap(),
                        current_start,
                        current_end_incl,
                    ));
                }
                current_chrom = Some(name.clone());
                current_start = ref_start;
                current_end_incl = end_incl;
                has_current = true;
            } else if runoff_gap {
                if let Some(ref c) = current_chrom {
                    regions.push((c.clone(), current_start, current_end_incl));
                }
                current_start = ref_start;
                current_end_incl = end_incl;
            } else if has_current {
                current_end_incl = current_end_incl.max(end_incl);
            }
        }

        // processRead:
        // - drop secondary/supplementary for long-read-class alignments
        // - singleton long-read artifact screening
        // - terminal poly exon trimming before duplicate matching / pair linking
        // - per-event mismatch accumulation for junction.nm
        if let Some(mut read) = parsed_read.take() {
            let debug_hit = debug_target_overlaps(
                debug_target.as_ref(),
                &name,
                read.ref_start,
                read.ref_end.saturating_sub(1),
            );
            if debug_hit {
                early_debug.parsed_reads += 1;
                accumulate_stranded_junction_counts(
                    &mut early_debug.pretrim_raw,
                    &read.junctions_raw,
                    read.strand,
                    read.weight,
                );
                accumulate_stranded_junction_counts(
                    &mut early_debug.pretrim_del,
                    &read.junctions_del,
                    read.strand,
                    read.weight,
                );
            }
            let is_secondary_or_supp =
                record.flags().is_secondary() || record.flags().is_supplementary();
            let this_long = read_is_long_class(&read, config);
            // VG flow-based multi-map redirection: in VG mode, retain secondary
            // alignments as per-copy evidence. Defaults to on whenever --vg is
            // set so the downstream EM solver can replace the uniform 1/NH
            // weighting with evidence-based (junction-compatibility + context)
            // weights per copy. Override with RUSTLE_VG_DROP_SECONDARY=1 to
            // restore StringTie-equivalent behaviour.
            let vg_include_secondary = (config.vg_mode
                || std::env::var_os("RUSTLE_VG_INCLUDE_SECONDARY").is_some())
                && std::env::var_os("RUSTLE_VG_DROP_SECONDARY").is_none()
                && is_secondary_or_supp
                && config.long_reads
                && this_long;
            let drop_sec_supp_for_mode =
                is_secondary_or_supp && config.long_reads && this_long && !vg_include_secondary;
            if trace_log_style && drop_sec_supp_for_mode {
                eprintln!(
                    "--- processRead: DROP name={} reason=SECONDARY_SUPPLEMENTARY flags=0x{:x}",
                    read.read_name,
                    record.flags().bits()
                );
            }
            if (!reads_primary_only || is_primary) && !drop_sec_supp_for_mode {
                // processRead: for de novo long-read flow (ovlpguide=false),
                // discard <=2 exon long reads with aligned poly artifact.
                if this_long
                    && read.exons.len() <= 2
                    && (read.has_poly_start_aligned || read.has_poly_end_aligned)
                {
                    if debug_hit {
                        early_debug.drop_poly_artifact += 1;
                    }
                    if trace_log_style {
                        eprintln!(
                            "--- processRead: DROP name={} reason=POLY_ARTIFACT flags=0x{:x}",
                            read.read_name,
                            record.flags().bits()
                        );
                    }
                    continue;
                }
                let pretrim_exons = if debug_hit {
                    Some(read.exons.clone())
                } else {
                    None
                };
                let pretrim_junctions_raw = if debug_hit {
                    Some(read.junctions_raw.clone())
                } else {
                    None
                };
                let pretrim_junctions_del = if debug_hit {
                    Some(read.junctions_del.clone())
                } else {
                    None
                };
                if this_long && !apply_processread_terminal_poly_trim(&mut read) {
                    if debug_hit {
                        early_debug.drop_trim_empty += 1;
                    }
                    if trace_log_style {
                        eprintln!(
                            "--- processRead: DROP name={} reason=TERMINAL_POLY_TRIM_EMPTY flags=0x{:x}",
                            read.read_name,
                            record.flags().bits()
                        );
                    }
                    continue;
                }
                if debug_hit {
                    if let (Some(before_exons), Some(before_raw), Some(before_del)) = (
                        pretrim_exons.as_ref(),
                        pretrim_junctions_raw.as_ref(),
                        pretrim_junctions_del.as_ref(),
                    ) {
                        if *before_exons != read.exons
                            || *before_raw != read.junctions_raw
                            || *before_del != read.junctions_del
                        {
                            eprintln!(
                                "[EARLY_DEBUG] trim_delta read={} before_exons=\"{}\" after_exons=\"{}\" before_raw=\"{}\" after_raw=\"{}\" before_del=\"{}\" after_del=\"{}\"",
                                read.read_name,
                                format_trace_exons(before_exons),
                                format_trace_exons(&read.exons),
                                format_trace_junctions(before_raw),
                                format_trace_junctions(&read.junctions_raw),
                                format_trace_junctions(before_del),
                                format_trace_junctions(&read.junctions_del),
                            );
                        }
                    }
                }
                if debug_hit {
                    early_debug.accepted_reads += 1;
                    accumulate_stranded_junction_counts(
                        &mut early_debug.posttrim_raw,
                        &read.junctions_raw,
                        read.strand,
                        read.weight,
                    );
                    accumulate_stranded_junction_counts(
                        &mut early_debug.posttrim_del,
                        &read.junctions_del,
                        read.strand,
                        read.weight,
                    );
                }
                if trace_log_style {
                    eprintln!(
                        "Process read {} with strand={} and exons: {}",
                        read.read_name,
                        trace_strand_i8(read.strand),
                        format_trace_exons(&read.exons)
                    );
                }

                // PARITY READ DUMP: emit read_name + final exon coords for cross-tool diff
                if let Ok(dump_path) = std::env::var("RUSTLE_PARITY_READ_DUMP_TSV") {
                    if !dump_path.is_empty() {
                        use std::io::Write;
                        use std::sync::Mutex;
                        use std::sync::OnceLock;
                        static DUMP_MUTEX: OnceLock<Mutex<bool>> = OnceLock::new();
                        let mtx = DUMP_MUTEX.get_or_init(|| Mutex::new(false));
                        if let Ok(mut hdr_written) = mtx.lock() {
                            if let Ok(mut f) = std::fs::OpenOptions::new()
                                .create(true)
                                .append(true)
                                .open(&dump_path)
                            {
                                if !*hdr_written {
                                    let _ = writeln!(f, "source\tchrom\tread_name\tstart\tend\tstrand\texons");
                                    *hdr_written = true;
                                }
                                let mut exons_str = String::new();
                                for (i, (s, e)) in read.exons.iter().enumerate() {
                                    if i > 0 { exons_str.push(','); }
                                    // StringTie uses 1-based inclusive coords; convert
                                    // Rustle's 0-based half-open [s, e) to 1-based [s+1, e].
                                    exons_str.push_str(&format!("{}-{}", s + 1, *e));
                                }
                                let strand_i = match read.strand {
                                    '+' => 1i32,
                                    '-' => -1i32,
                                    _ => 0i32,
                                };
                                let _ = writeln!(f, "rustle\t{}\t{}\t{}\t{}\t{}\t{}",
                                    &name,
                                    &*read.read_name,
                                    read.ref_start + 1,
                                    read.ref_end,
                                    strand_i,
                                    exons_str);
                            }
                        }
                    }
                }

                // mismatch gate is per alignment event (before collapsing into readlist entry).
                let anchor_start = if has_current {
                    current_start
                } else {
                    read.ref_start
                };
                read.junc_mismatch_weight = processread_mismatch_weight(
                    &read,
                    this_long,
                    anchor_start,
                    config.junction_support,
                    config.mismatchfrac,
                );

                let event_weight = read.weight;
                let chrom_reads = reads_by_chrom.entry(name.clone()).or_default();
                // exonmatch/deljuncmatch: if an immediately preceding
                // same-start read has identical exon structure (and deletion-aware
                // junction offsets for long/mixed), aggregate into that read record.
                let mut matched_idx: Option<usize> = None;
                for idx in (0..chrom_reads.len()).rev() {
                    let prev = &chrom_reads[idx];
                    if prev.ref_start != read.ref_start {
                        break;
                    }
                    if debug_hit {
                        early_debug.collapse_same_start_candidates += 1;
                    }
                    if prev.strand != read.strand {
                        continue;
                    }
                    if read_is_long_class(prev, config) != this_long {
                        continue;
                    }
                    if !exonmatch(&prev.exons, &read.exons) {
                        if debug_hit {
                            early_debug.collapse_exon_miss += 1;
                            if early_debug.collapse_exon_miss_logged < 20 {
                                eprintln!(
                                    "[EARLY_DEBUG] collapse_exon_miss read={} prev_read={} ref_start={} prev_exons=\"{}\" read_exons=\"{}\"",
                                    read.read_name,
                                    prev.read_name,
                                    read.ref_start,
                                    format_trace_exons(&prev.exons),
                                    format_trace_exons(&read.exons),
                                );
                                early_debug.collapse_exon_miss_logged += 1;
                            }
                        }
                        continue;
                    }
                    if config.long_reads && !deljuncmatch(prev, &read.junctions_del) {
                        if debug_hit {
                            early_debug.collapse_deljunc_miss += 1;
                            if early_debug.collapse_deljunc_miss_logged < 20 {
                                eprintln!(
                                    "[EARLY_DEBUG] collapse_deljunc_miss read={} prev_read={} ref_start={} exons=\"{}\" prev_del=\"{}\" read_del=\"{}\" prev_jdel=\"{}\" read_jdel=\"{}\"",
                                    read.read_name,
                                    prev.read_name,
                                    read.ref_start,
                                    format_trace_exons(&read.exons),
                                    format_deljunc_offsets(&prev.exons, &prev.junctions_del),
                                    format_deljunc_offsets(&read.exons, &read.junctions_del),
                                    format_trace_junctions(&prev.junctions_del),
                                    format_trace_junctions(&read.junctions_del),
                                );
                                early_debug.collapse_deljunc_miss_logged += 1;
                            }
                        }
                        continue;
                    }
                    matched_idx = Some(idx);
                    break;
                }
                if let Some(mi) = matched_idx {
                    if debug_hit {
                        early_debug.collapsed_reads += 1;
                    }
                    if let Some(prev) = chrom_reads.get_mut(mi) {
                        prev.weight += read.weight;
                        prev.countfrag_len += read.countfrag_len;
                        prev.countfrag_num += read.countfrag_num;
                        prev.junc_mismatch_weight += read.junc_mismatch_weight;
                        // retain smallest NH across redundant alignments.
                        if read.nh < prev.nh {
                            prev.nh = read.nh;
                        }
                        // poly evidence is accumulated only for long-read entries.
                        if this_long {
                            prev.has_poly_start |= read.has_poly_start;
                            prev.has_poly_end |= read.has_poly_end;
                            prev.has_poly_start_aligned |= read.has_poly_start_aligned;
                            prev.has_poly_start_unaligned |= read.has_poly_start_unaligned;
                            prev.has_poly_end_aligned |= read.has_poly_end_aligned;
                            prev.has_poly_end_unaligned |= read.has_poly_end_unaligned;
                            // Redundant long-read alignments contribute additional tail-evidence
                            // counts on the same read entry (uint16 saturated).
                            prev.unaligned_poly_t =
                                prev.unaligned_poly_t.saturating_add(read.unaligned_poly_t);
                            prev.unaligned_poly_a =
                                prev.unaligned_poly_a.saturating_add(read.unaligned_poly_a);
                        }
                    }
                    if let Some(prev) = chrom_reads.get_mut(mi) {
                        if prev.read_name.is_empty() {
                            prev.read_name = read.read_name.clone();
                        }
                        if prev.ref_id.is_none() {
                            prev.ref_id = read.ref_id;
                        }
                        if prev.mate_ref_id.is_none() {
                            prev.mate_ref_id = read.mate_ref_id;
                        }
                        if prev.mate_start.is_none() {
                            prev.mate_start = read.mate_start;
                        }
                        if prev.hi == 0 {
                            prev.hi = read.hi;
                        }
                    }

                    if record.flags().is_segmented() && !record.flags().is_mate_unmapped() {
                        if let (Some(rid), Some(mrid), Some(mstart)) =
                            (read.ref_id, read.mate_ref_id, read.mate_start)
                        {
                            if rid == mrid {
                                let key = PairKey {
                                    name: read.read_name.clone(),
                                    lo_start: read.ref_start.min(mstart),
                                    hi_start: read.ref_start.max(mstart),
                                    hi_tag: read.hi,
                                };
                                let maybe_pair =
                                    pop_pending_pair(&mut pair_pending_by_chrom, &name, &key);
                                if let Some((mate_idx, _mate_weight)) = maybe_pair {
                                    let mate_nh =
                                        chrom_reads.get(mate_idx).map(|r| r.nh).unwrap_or(read.nh);
                                    let pair_w =
                                        processread_pair_weight(event_weight, read.nh, mate_nh);
                                    add_pair_link(chrom_reads, mi, mate_idx, pair_w);
                                } else {
                                    push_pending_pair(
                                        &mut pair_pending_by_chrom,
                                        &name,
                                        key,
                                        mi,
                                        event_weight,
                                    );
                                }
                            }
                        }
                    }
                    if trace_log_style {
                        eprintln!(
                            "Add read {} with strand={} and exons: {}",
                            read.read_name,
                            trace_strand_i8(read.strand),
                            format_trace_exons(&read.exons)
                        );
                    }
                } else {
                    chrom_reads.push(read);
                    let new_idx = chrom_reads.len() - 1;
                    let inserted = &chrom_reads[new_idx];
                    if trace_log_style {
                        eprintln!(
                            "Add read {} with strand={} and exons: {}",
                            inserted.read_name,
                            trace_strand_i8(inserted.strand),
                            format_trace_exons(&inserted.exons)
                        );
                    }
                    if record.flags().is_segmented() && !record.flags().is_mate_unmapped() {
                        if let (Some(rid), Some(mrid), Some(mstart)) =
                            (inserted.ref_id, inserted.mate_ref_id, inserted.mate_start)
                        {
                            if rid == mrid {
                                let key = PairKey {
                                    name: inserted.read_name.clone(),
                                    lo_start: inserted.ref_start.min(mstart),
                                    hi_start: inserted.ref_start.max(mstart),
                                    hi_tag: inserted.hi,
                                };
                                let maybe_pair =
                                    pop_pending_pair(&mut pair_pending_by_chrom, &name, &key);
                                if let Some((mate_idx, _mate_weight)) = maybe_pair {
                                    let mate_nh = chrom_reads
                                        .get(mate_idx)
                                        .map(|r| r.nh)
                                        .unwrap_or(inserted.nh);
                                    let pair_w =
                                        processread_pair_weight(event_weight, inserted.nh, mate_nh);
                                    add_pair_link(chrom_reads, new_idx, mate_idx, pair_w);
                                } else {
                                    push_pending_pair(
                                        &mut pair_pending_by_chrom,
                                        &name,
                                        key,
                                        new_idx,
                                        event_weight,
                                    );
                                }
                            }
                        }
                    }
                }
            } else if debug_hit && drop_sec_supp_for_mode {
                early_debug.drop_secondary_supp += 1;
            }
        }
    }

    if has_current {
        if let Some(c) = current_chrom {
            regions.push((c, current_start, current_end_incl));
        }
    }

    // Optional: post-pass that splits regions at zero-coverage stretches >= min_gap.
    // Uses exon-level coverage from primary reads (introns don't contribute coverage).
    // Activated by RUSTLE_BPCOV_REGION_SPLIT=<min_gap_bp> (e.g. 200).
    if let Some(min_gap) = std::env::var("RUSTLE_BPCOV_REGION_SPLIT")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .filter(|&v| v > 0)
    {
        let mut new_regions: Vec<(String, u64, u64)> = Vec::with_capacity(regions.len());
        for (chrom, start, end_incl) in regions.drain(..) {
            let span = end_incl.saturating_sub(start).saturating_add(1) as usize;
            if span <= min_gap as usize {
                new_regions.push((chrom, start, end_incl));
                continue;
            }
            let Some(reads) = reads_by_chrom.get(&chrom) else {
                new_regions.push((chrom, start, end_incl));
                continue;
            };
            // Exon-level coverage bitmap for this region.
            let mut covered = vec![false; span];
            for r in reads {
                if r.ref_end <= start || r.ref_start > end_incl {
                    continue;
                }
                for (s, e) in &r.exons {
                    let es = (*s).max(start);
                    let ee = (*e).min(end_incl);
                    if es > ee {
                        continue;
                    }
                    let a = (es - start) as usize;
                    let b = (ee - start) as usize;
                    for idx in a..=b.min(span - 1) {
                        covered[idx] = true;
                    }
                }
            }
            // Scan for zero-cov stretches >= min_gap.
            let mut cur_start = start;
            let mut last_cov: Option<u64> = None;
            let mut in_gap_start: Option<usize> = None;
            for i in 0..span {
                if covered[i] {
                    if let Some(gs) = in_gap_start.take() {
                        let gap_len = (i - gs) as u64;
                        if gap_len >= min_gap {
                            if let Some(lc) = last_cov {
                                new_regions.push((chrom.clone(), cur_start, lc));
                                cur_start = start + i as u64;
                            }
                        }
                    }
                    last_cov = Some(start + i as u64);
                } else if in_gap_start.is_none() {
                    in_gap_start = Some(i);
                }
            }
            new_regions.push((chrom, cur_start, end_incl));
        }
        regions = new_regions;
    }

    if let Some(target) = debug_target.as_ref() {
        if let Some(reads) = reads_by_chrom.get(&target.chrom) {
            for r in reads.iter().filter(|r| {
                intervals_overlap(
                    r.ref_start,
                    r.ref_end.saturating_sub(1),
                    target.start,
                    target.end,
                )
            }) {
                early_debug.stored_reads += 1;
                accumulate_stranded_junction_counts(
                    &mut early_debug.stored_raw,
                    &r.junctions_raw,
                    r.strand,
                    r.weight,
                );
                accumulate_stranded_junction_counts(
                    &mut early_debug.stored_del,
                    &r.junctions_del,
                    r.strand,
                    r.weight,
                );
            }
        }
        eprintln!(
            "[EARLY_DEBUG] ingest {}:{}-{} parsed_reads={} accepted_reads={} collapsed_reads={} stored_reads={} collapse_same_start_candidates={} collapse_exon_miss={} collapse_deljunc_miss={} drop_secondary_supp={} drop_poly_artifact={} drop_trim_empty={}",
            target.chrom,
            target.start,
            target.end,
            early_debug.parsed_reads,
            early_debug.accepted_reads,
            early_debug.collapsed_reads,
            early_debug.stored_reads,
            early_debug.collapse_same_start_candidates,
            early_debug.collapse_exon_miss,
            early_debug.collapse_deljunc_miss,
            early_debug.drop_secondary_supp,
            early_debug.drop_poly_artifact,
            early_debug.drop_trim_empty
        );
        dump_stranded_junction_counts("pretrim_raw", &early_debug.pretrim_raw);
        dump_stranded_junction_counts("pretrim_del", &early_debug.pretrim_del);
        dump_stranded_junction_counts("posttrim_raw", &early_debug.posttrim_raw);
        dump_stranded_junction_counts("posttrim_del", &early_debug.posttrim_del);
        dump_stranded_junction_counts("stored_raw", &early_debug.stored_raw);
        dump_stranded_junction_counts("stored_del", &early_debug.stored_del);
    }

    // Create per-strand bundles does
    // Each strand gets its own bundle with its own junction_stats.
    // Cross-strand color mapping (eqposcol/eqnegcol) happens at the CGroup level.
    let mut bundles: Vec<Bundle> = Vec::new();
    for (chrom, start, end_incl) in regions {
        let reads = match reads_by_chrom.get(&chrom) {
            Some(r) => r,
            None => continue,
        };
        let bundle_read_entries: Vec<(usize, BundleRead)> = reads
            .iter()
            .enumerate()
            .filter(|(_, r)| r.ref_start < end_incl.saturating_add(1) && r.ref_end > start)
            .map(|(i, r)| (i, r.clone()))
            .collect();
        let mut orig_to_local: HashMap<usize, usize> = Default::default();
        for (local_idx, (orig_idx, _)) in bundle_read_entries.iter().enumerate() {
            orig_to_local.insert(*orig_idx, local_idx);
        }
        let mut bundle_reads: Vec<BundleRead> =
            bundle_read_entries.into_iter().map(|(_, r)| r).collect();
        for r in bundle_reads.iter_mut() {
            let old_idx = std::mem::take(&mut r.pair_idx);
            let old_count = std::mem::take(&mut r.pair_count);
            for (i, pidx) in old_idx.into_iter().enumerate() {
                if let Some(&mapped) = orig_to_local.get(&pidx) {
                    r.pair_idx.push(mapped);
                    r.pair_count.push(old_count.get(i).copied().unwrap_or(0.0));
                }
            }
        }
        if bundle_reads.is_empty() {
            continue;
        }

        // Create separate bundles for each strand (
        // sno: 0 = negative, 1 = unstranded, 2 = positive
        for strand in ['-', '.', '+'] {
            let strand_reads: Vec<BundleRead> = bundle_reads
                .iter()
                .filter(|r| {
                    if strand == '.' {
                        r.strand == strand || r.strand == '?' || r.strand == ' '
                    } else {
                        r.strand == strand
                    }
                })
                .cloned()
                .collect();

            if strand_reads.is_empty() {
                continue;
            }

            // Compute junction_stats for this strand
            let mut junction_stats: JunctionStats = Default::default();

            for r in &strand_reads {
                // compute cumulative max exon lengths from left and right
                // (count_good_junctions: leftsup[i] = max seg length from start to i,
                // rightsup[nex-i-1] = max seg length from end to nex-i).
                let nex = r.exons.len();
                let mut leftsup = vec![0u64; nex];
                let mut rightsup = vec![0u64; nex];
                if nex > 0 {
                    let mut max_left = 0u64;
                    let mut max_right = 0u64;
                    for ei in 0..nex {
                        let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                        if seg_len > max_left {
                            max_left = seg_len;
                        }
                        leftsup[ei] = max_left;
                    }
                    for ei in (0..nex).rev() {
                        let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                        if seg_len > max_right {
                            max_right = seg_len;
                        }
                        rightsup[ei] = max_right;
                    }
                }

                for i in 0..r.junctions.len() {
                    let j = r.junctions[i];
                    let intron_len = j.acceptor.saturating_sub(j.donor);
                    if intron_len < min_intron {
                        continue;
                    }
                    if !(j.donor >= start
                        && j.donor <= end_incl
                        && j.acceptor >= start
                        && j.acceptor <= end_incl)
                    {
                        continue;
                    }
                    // leftsup[i-1] → cumulative max exon length from exon 0..=i
                    // rightsup[nex-i-1] → cumulative max exon length from exon nex-1..=i+1
                    let left_anchor = leftsup.get(i).copied().unwrap_or(0);
                    let right_anchor = rightsup.get(i + 1).copied().unwrap_or(0);
                    let st = junction_stats
                        .entry(j)
                        .or_insert_with(JunctionStat::default);
                    st.mrcount += r.weight;
                    let mut anchor = config.junction_support;
                    if intron_len > config.longintron && anchor < LONGINTRONANCHOR {
                        anchor = LONGINTRONANCHOR;
                    }
                    // reference-style anchor witnesses: for long introns, raise the minimum anchor
                    // to LONGINTRONANCHOR before counting left/right support.
                    if left_anchor >= anchor {
                        st.leftsupport += r.weight;
                        if right_anchor >= anchor {
                            st.rightsupport += r.weight;
                            // count_good_junctions: nreads_good increments only when both sides pass anchor.
                            st.nreads_good += r.weight;
                        }
                    } else if right_anchor >= anchor {
                        st.rightsupport += r.weight;
                    }
                    // processRead: nm is accumulated per event (`mismatch || nh>2`),
                    // then carried by merged read entries as junc_mismatch_weight.
                    st.nm += r.junc_mismatch_weight;
                    // mm counts reads where both flanking exons > LONGINTRONANCHOR
                    if leftsup.get(i).copied().unwrap_or(0) > LONGINTRONANCHOR
                        && rightsup.get(i + 1).copied().unwrap_or(0) > LONGINTRONANCHOR
                    {
                        st.mm += r.weight;
                    }
                    if st.strand.is_none() {
                        // Keep historical Rust behavior for now: unresolved read strands are
                        // carried as +1 in initial bundle stats so junctions are not dropped
                        // before downstream correction/resolution passes.
                        st.strand = Some(if r.strand == '-' { -1 } else { 1 });
                    }
                    if std::env::var_os("RUSTLE_JUNC_TRACE").is_some()
                        && j.donor >= 20648000
                        && j.acceptor <= 20733000
                    {
                        eprintln!(
                            "JUNC_ADD read_strand={} jstart={} jend={} obj_strand={:?}",
                            r.strand, j.donor, j.acceptor, st.strand
                        );
                    }
                }
            }

            if debug_target_overlaps(debug_target.as_ref(), &chrom, start, end_incl) {
                let mut raw_counts: HashMap<(u64, u64, i8), f64> = Default::default();
                let mut del_counts: HashMap<(u64, u64, i8), f64> = Default::default();
                for r in &strand_reads {
                    accumulate_stranded_junction_counts(
                        &mut raw_counts,
                        &r.junctions_raw,
                        r.strand,
                        r.weight,
                    );
                    accumulate_stranded_junction_counts(
                        &mut del_counts,
                        &r.junctions_del,
                        r.strand,
                        r.weight,
                    );
                }
                eprintln!(
                    "[EARLY_DEBUG] strand_bundle {}:{}-{} strand={} reads={} initial_jstats={}",
                    chrom,
                    start,
                    end_incl,
                    strand,
                    strand_reads.len(),
                    junction_stats.len()
                );
                dump_stranded_junction_counts("strand_bundle_raw", &raw_counts);
                dump_stranded_junction_counts("strand_bundle_del", &del_counts);
            }

            bundles.push(Bundle {
                chrom: chrom.clone(),
                start,
                end: end_incl,
                strand,
                reads: strand_reads,
                junction_stats,
                bundlenodes: None,
                read_bnodes: None,
                bnode_colors: None,
                synthetic: false,
                rescue_class: None,
            });
        } // end for strand
    } // end for regions

    bundles.sort_unstable_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
    Ok(bundles)
}

/// Parsed region from the original algorithm bundle log: >bundle CHROM:START-END. Log is 1-based inclusive.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BundleLogRegion {
    pub chrom: String,
    /// 0-based inclusive start (start_1b - 1)
    pub start: u64,
    /// 0-based inclusive end (end_1b - 1)
    pub end: u64,
}

/// Parse the original algorithm -v log for bundle lines. Returns 0-based inclusive (start, end).
pub fn parse_bundle_log<P: AsRef<Path>>(log_path: P) -> Result<Vec<BundleLogRegion>> {
    let s = std::fs::read_to_string(log_path)?;
    let mut out = Vec::new();
    for line in s.lines() {
        let Some(after) = line.find(">bundle ") else {
            continue;
        };
        let rest = line[after + ">bundle ".len()..].trim();
        let region = rest.split_whitespace().next().unwrap_or("");
        let mut parts = region.split(':');
        let chrom = parts.next().unwrap_or("").to_string();
        let range = parts.next().unwrap_or("");
        let mut se = range.split('-');
        let start_1b: u64 = se.next().and_then(|x| x.parse().ok()).unwrap_or(0);
        let end_1b: u64 = se.next().and_then(|x| x.parse().ok()).unwrap_or(0);
        if chrom.is_empty() || end_1b < start_1b {
            continue;
        }
        out.push(BundleLogRegion {
            chrom,
            start: start_1b.saturating_sub(1),
            end: end_1b.saturating_sub(1),
        });
    }
    Ok(out)
}

/// Compare detected bundles to log regions. Our bundles are per (chrom, start, end, strand); log has (chrom, start, end).
/// Reports match count and only-in-Rust / only-in-log.
pub fn compare_bundles_to_log(bundles: &[Bundle], log_regions: &[BundleLogRegion]) {
    let ours: HashSet<(String, u64, u64)> = bundles
        .iter()
        .map(|b| (b.chrom.clone(), b.start, b.end))
        .collect();
    let theirs: HashSet<(String, u64, u64)> = log_regions
        .iter()
        .map(|r| (r.chrom.clone(), r.start, r.end))
        .collect();
    let matched = ours.intersection(&theirs).count();
    let only_rust = ours.difference(&theirs).count();
    let only_log = theirs.difference(&ours).count();
    eprintln!(
        "[bundle compare] ours={} regions, log={} regions | match={}, only_rust={}, only_log={}",
        ours.len(),
        theirs.len(),
        matched,
        only_rust,
        only_log
    );

    let mut only_rust_regions: Vec<(String, u64, u64)> =
        ours.difference(&theirs).cloned().collect();
    let mut only_log_regions: Vec<(String, u64, u64)> = theirs.difference(&ours).cloned().collect();
    only_rust_regions.sort_unstable();
    only_log_regions.sort_unstable();
    let max_show = std::env::var("RUSTLE_BUNDLE_DIFF_MAX")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(20);
    if !only_rust_regions.is_empty() {
        eprintln!(
            "[bundle compare] only_rust_regions (showing up to {}):",
            max_show
        );
        for (chrom, start, end) in only_rust_regions.into_iter().take(max_show) {
            eprintln!(
                "  rust {}:{}-{}",
                chrom,
                start.saturating_add(1),
                end.saturating_add(1)
            );
        }
    }
    if !only_log_regions.is_empty() {
        eprintln!(
            "[bundle compare] only_log_regions (showing up to {}):",
            max_show
        );
        for (chrom, start, end) in only_log_regions.into_iter().take(max_show) {
            eprintln!(
                "  log  {}:{}-{}",
                chrom,
                start.saturating_add(1),
                end.saturating_add(1)
            );
        }
    }
}

/// Rebuild a bundle's `junction_stats` from a filtered view of its `reads`,
/// applying the same anchor/intron-length logic used at ingest time.
///
/// `primary_only`: when true, junction stats are built only from primary
/// alignments. Used to clean cross-mapping noise in --vg mode where secondary
/// alignments from other paralogs leave wrong-intron junctions in the splice
/// graph. The reads vec itself is left untouched — only `junction_stats` is
/// rewritten.
pub fn recompute_junction_stats(bundle: &mut Bundle, config: &RunConfig) {
    recompute_junction_stats_inner(bundle, config, /*primary_only*/ false);
}

/// Rebuild `junction_stats` using only primary alignments — secondaries are
/// dropped entirely. Used in --vg mode AFTER family discovery to clean
/// cross-mapping junction noise from non-family bundles where secondaries
/// have no useful role.
pub fn recompute_junction_stats_primary_only(bundle: &mut Bundle, config: &RunConfig) {
    recompute_junction_stats_inner(bundle, config, /*primary_only*/ true);
}

/// Rebuild `junction_stats` using ALL reads, but only count junctions that
/// at least one PRIMARY read supports. This is the right knob for family
/// bundles: it drops "wrong-intron" cross-mapping junctions (introns that
/// no primary read at this locus uses — they came from a sister paralog's
/// structure when its read was secondary-aligned here) while preserving
/// real junctions where primaries AND secondaries agree. The
/// secondary-confirmed counts on real junctions stay intact so EM and
/// flow-based assembly still see the cross-mapping evidence.
///
/// This is what the GOLGA6L10 fix needs WITHOUT the LOC129523543 collateral
/// damage: GOLGA6L10's wrong-intron secondaries (8228N/858N/4557N/210N)
/// have no primary support → dropped. LOC129523543's primary+secondary
/// supported junctions keep their full counts → paralog still recovered.
pub fn recompute_junction_stats_primary_supported(bundle: &mut Bundle, config: &RunConfig) {
    let primary_supported: std::collections::HashSet<crate::types::Junction> = bundle
        .reads
        .iter()
        .filter(|r| r.is_primary_alignment)
        .flat_map(|r| r.junctions.iter().copied())
        .collect();
    recompute_junction_stats_filtered(bundle, config, &primary_supported);
}

fn recompute_junction_stats_filtered(
    bundle: &mut Bundle,
    config: &RunConfig,
    keep_set: &std::collections::HashSet<crate::types::Junction>,
) {
    let min_intron = config.min_intron_length;
    let start = bundle.start;
    let end_incl = bundle.end;
    let mut new_stats: JunctionStats = Default::default();
    for r in &bundle.reads {
        let nex = r.exons.len();
        let mut leftsup = vec![0u64; nex];
        let mut rightsup = vec![0u64; nex];
        if nex > 0 {
            let mut max_left = 0u64;
            let mut max_right = 0u64;
            for ei in 0..nex {
                let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                if seg_len > max_left { max_left = seg_len; }
                leftsup[ei] = max_left;
            }
            for ei in (0..nex).rev() {
                let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                if seg_len > max_right { max_right = seg_len; }
                rightsup[ei] = max_right;
            }
        }
        for i in 0..r.junctions.len() {
            let j = r.junctions[i];
            if !keep_set.contains(&j) { continue; }
            let intron_len = j.acceptor.saturating_sub(j.donor);
            if intron_len < min_intron { continue; }
            if !(j.donor >= start && j.donor <= end_incl
                && j.acceptor >= start && j.acceptor <= end_incl) { continue; }
            let left_anchor = leftsup.get(i).copied().unwrap_or(0);
            let right_anchor = rightsup.get(i + 1).copied().unwrap_or(0);
            let st = new_stats.entry(j).or_insert_with(JunctionStat::default);
            st.mrcount += r.weight;
            let mut anchor = config.junction_support;
            if intron_len > config.longintron && anchor < LONGINTRONANCHOR {
                anchor = LONGINTRONANCHOR;
            }
            if left_anchor >= anchor {
                st.leftsupport += r.weight;
                if right_anchor >= anchor {
                    st.rightsupport += r.weight;
                    st.nreads_good += r.weight;
                }
            } else if right_anchor >= anchor {
                st.rightsupport += r.weight;
            }
            st.nm += r.junc_mismatch_weight;
            if leftsup.get(i).copied().unwrap_or(0) > LONGINTRONANCHOR
                && rightsup.get(i + 1).copied().unwrap_or(0) > LONGINTRONANCHOR
            {
                st.mm += r.weight;
            }
            if st.strand.is_none() {
                st.strand = Some(if r.strand == '-' { -1 } else { 1 });
            }
        }
    }
    bundle.junction_stats = new_stats;
}

fn recompute_junction_stats_inner(bundle: &mut Bundle, config: &RunConfig, primary_only: bool) {
    let min_intron = config.min_intron_length;
    let start = bundle.start;
    let end_incl = bundle.end;

    let mut new_stats: JunctionStats = Default::default();
    for r in &bundle.reads {
        if primary_only && !r.is_primary_alignment {
            continue;
        }
        let nex = r.exons.len();
        let mut leftsup = vec![0u64; nex];
        let mut rightsup = vec![0u64; nex];
        if nex > 0 {
            let mut max_left = 0u64;
            let mut max_right = 0u64;
            for ei in 0..nex {
                let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                if seg_len > max_left {
                    max_left = seg_len;
                }
                leftsup[ei] = max_left;
            }
            for ei in (0..nex).rev() {
                let seg_len = r.exons[ei].1.saturating_sub(r.exons[ei].0);
                if seg_len > max_right {
                    max_right = seg_len;
                }
                rightsup[ei] = max_right;
            }
        }
        for i in 0..r.junctions.len() {
            let j = r.junctions[i];
            let intron_len = j.acceptor.saturating_sub(j.donor);
            if intron_len < min_intron {
                continue;
            }
            if !(j.donor >= start
                && j.donor <= end_incl
                && j.acceptor >= start
                && j.acceptor <= end_incl)
            {
                continue;
            }
            let left_anchor = leftsup.get(i).copied().unwrap_or(0);
            let right_anchor = rightsup.get(i + 1).copied().unwrap_or(0);
            let st = new_stats
                .entry(j)
                .or_insert_with(JunctionStat::default);
            st.mrcount += r.weight;
            let mut anchor = config.junction_support;
            if intron_len > config.longintron && anchor < LONGINTRONANCHOR {
                anchor = LONGINTRONANCHOR;
            }
            if left_anchor >= anchor {
                st.leftsupport += r.weight;
                if right_anchor >= anchor {
                    st.rightsupport += r.weight;
                    st.nreads_good += r.weight;
                }
            } else if right_anchor >= anchor {
                st.rightsupport += r.weight;
            }
            st.nm += r.junc_mismatch_weight;
            if leftsup.get(i).copied().unwrap_or(0) > LONGINTRONANCHOR
                && rightsup.get(i + 1).copied().unwrap_or(0) > LONGINTRONANCHOR
            {
                st.mm += r.weight;
            }
            if st.strand.is_none() {
                st.strand = Some(if r.strand == '-' { -1 } else { 1 });
            }
        }
    }
    bundle.junction_stats = new_stats;
}
