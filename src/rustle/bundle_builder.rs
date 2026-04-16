//! Faithful port of the original algorithm bundling logic from 
//! This module implements build_graphs bundle creation with exact

use anyhow::Result;
use hashbrown::{HashMap as HbHashMap, HashSet as HbHashSet};
use roaring::RoaringBitmap;
use std::time::Instant;
use std::sync::Arc;

use crate::types::{BundleRead, CBundlenode, DetHashSet as HashSet, Junction, RunConfig};

/// CGroup structure matching header:145-155 exactly
/// next_gr is the index of the next group in the linked list (like pointer)
#[derive(Debug, Clone)]
struct CGroup {
    start: u64,
    end: u64,
    color: usize,
    grid: usize,
    cov_sum: f64,
    multi: f64,
    neg_prop: f64, // Proportion of negative strand overlap (for unstranded groups)
    next_gr: Option<usize>, // Index of next group in linked list, not Box
    hardstart: bool,
    hardend: bool,
}

impl CGroup {
    fn new(start: u64, end: u64, color: usize, grid: usize, cov_sum: f64, multi: f64) -> Self {
        Self {
            start,
            end,
            color,
            grid,
            cov_sum,
            multi,
            neg_prop: 0.0,
            next_gr: None,
            hardstart: false,
            hardend: false,
        }
    }
}

/// CBundle structure matching header:260-268
#[derive(Debug, Clone)]
struct CBundle {
    _len: u64,
    _cov: f64,
    _multi: f64,
    _startnode: i32,  // id of first bundlenode
    _lastnodeid: i32, // id of last bundlenode added
}

impl CBundle {
    #[allow(dead_code)]
    fn new(len: u64, cov: f64, multi: f64, startnode: i32, lastnodeid: i32) -> Self {
        Self {
            _len: len,
            _cov: cov,
            _multi: multi,
            _startnode: startnode,
            _lastnodeid: lastnodeid,
        }
    }
}

/// Per-strand bundle data (sno: 0=negative, 1=unstranded, 2=positive)
struct StrandBundles {
    groups: Vec<CGroup>,       // All groups for this strand
    _bundles: Vec<CBundle>,    // All bundles for this strand
    _bnodes: Vec<CBundlenode>, // All bundlenodes for this strand
    startgroup: Option<usize>, // Index of first group in linked list
    currgroup: Option<usize>,  // Index of current group
}

impl StrandBundles {
    fn new() -> Self {
        Self {
            groups: Vec::new(),
            _bundles: Vec::new(),
            _bnodes: Vec::new(),
            startgroup: None,
            currgroup: None,
        }
    }
}

/// Result of Original-style bundle creation
pub struct SubBundleResult {
    pub strand: char,
    pub start: u64,
    pub end: u64,
    pub bnode_head: Option<CBundlenode>,
    pub read_to_bnodes: Vec<Vec<usize>>, // Maps read index to bundlenode IDs
    pub bnode_colors: Vec<usize>,        // Color for each bundlenode (by bid)
    pub read_scale: Vec<f64>,            // Per-read strand proportion (rprop)
}

/// Union-find find operation
fn find_color_root(color: usize, eqcol: &[usize]) -> usize {
    let mut c = color;
    let mut steps = 0usize;
    while c < eqcol.len() && eqcol[c] != c {
        let next = eqcol[c];
        if next >= eqcol.len() {
            break;
        }
        c = next;
        steps += 1;
        if steps > eqcol.len() {
            // Cycle guard for corrupted union-find links.
            return color.min(eqcol.len().saturating_sub(1));
        }
    }
    c
}

#[inline]
fn find_color_root_mut(color: usize, eqcol: &mut [usize]) -> usize {
    if color >= eqcol.len() {
        return color;
    }
    let mut root = color;
    let mut steps = 0usize;
    while root < eqcol.len() && eqcol[root] != root {
        let next = eqcol[root];
        if next >= eqcol.len() {
            break;
        }
        root = next;
        steps += 1;
        if steps > eqcol.len() {
            // Cycle guard for corrupted links.
            root = color;
            eqcol[color] = color;
            break;
        }
    }
    let mut cur = color;
    let mut compress_steps = 0usize;
    while cur < eqcol.len() && eqcol[cur] != cur {
        let next = eqcol[cur];
        eqcol[cur] = root;
        if next >= eqcol.len() {
            break;
        }
        cur = next;
        compress_steps += 1;
        if compress_steps > eqcol.len() {
            break;
        }
    }
    root
}

#[inline]
fn find_merge_root(idx: usize, merge: &mut [usize]) -> usize {
    if idx >= merge.len() {
        return idx;
    }
    let mut root = idx;
    let mut steps = 0usize;
    while root < merge.len() && merge[root] != root {
        let next = merge[root];
        if next >= merge.len() {
            break;
        }
        root = next;
        steps += 1;
        if steps > merge.len() {
            // Cycle guard for corrupted union-find links.
            root = idx;
            merge[idx] = idx;
            break;
        }
    }
    let mut cur = idx;
    let mut compress_steps = 0usize;
    while cur < merge.len() && merge[cur] != cur {
        let next = merge[cur];
        merge[cur] = root;
        if next >= merge.len() {
            break;
        }
        cur = next;
        compress_steps += 1;
        if compress_steps > merge.len() {
            break;
        }
    }
    root
}

#[inline]
fn find_merge_root_const(idx: usize, merge: &[usize]) -> usize {
    let mut root = idx;
    let mut steps = 0usize;
    while root < merge.len() && merge[root] != root {
        let next = merge[root];
        if next >= merge.len() {
            break;
        }
        root = next;
        steps += 1;
        if steps > merge.len() {
            return idx.min(merge.len().saturating_sub(1));
        }
    }
    root
}

fn group_color_from_grid(
    grid: usize,
    grid_to_group: &[(usize, usize)],
    strand_data: &mut [StrandBundles; 3],
    eqcol: &mut [usize],
) -> Option<usize> {
    let &(sno, gidx) = grid_to_group.get(grid)?;
    let group = strand_data.get_mut(sno)?.groups.get_mut(gidx)?;
    let c = find_color_root_mut(group.color, eqcol);
    group.color = c;
    Some(c)
}

/// Get strand index: 0=negative('-'), 1=unstranded('.'), 2=positive('+')
fn strand_idx(strand: char) -> usize {
    match strand {
        '-' => 0,
        '+' => 2,
        _ => 1, // '.' or any other
    }
}

/// Get next group index (strand with minimum start coordinate)
/// Returns (strand_idx, group_idx) for the strand with smallest start
/// Returns None if no groups remain
fn get_min_start(
    currgroup: &[Option<usize>], // [sno] -> group index in strand_data[sno].groups
    strand_data: &[StrandBundles; 3],
) -> Option<(usize, usize)> {
    let mut min_start = u64::MAX;
    let mut min_sno: Option<usize> = None;

    for sno in 0..3 {
        if let Some(gidx) = currgroup[sno] {
            if gidx < strand_data[sno].groups.len() {
                let start = strand_data[sno].groups[gidx].start;
                if start < min_start {
                    min_start = start;
                    min_sno = Some(sno);
                }
            }
        }
    }

    min_sno.map(|sno| (sno, currgroup[sno].unwrap()))
}

fn count_reachable_groups(sb: &StrandBundles) -> (usize, bool) {
    let mut seen: HbHashSet<usize> = HbHashSet::new();
    let mut cur = sb.startgroup;
    while let Some(idx) = cur {
        if idx >= sb.groups.len() {
            return (seen.len(), true);
        }
        if !seen.insert(idx) {
            // cycle in linked list
            return (seen.len(), true);
        }
        cur = sb.groups[idx].next_gr;
    }
    (seen.len(), false)
}

/// Set cross-strand color equivalence (set_strandcol 
/// `eqmap` is either eqnegcol or eqposcol and is keyed by the overlapped stranded color.
fn set_strandcol(
    prev_color: usize,
    group: &mut CGroup,
    grcol: usize,
    eqmap: &mut Vec<i64>,
    equalcolor: &mut [usize],
) {
    if prev_color >= eqmap.len() || grcol >= equalcolor.len() {
        return;
    }
    let zerocol = eqmap[prev_color];
    if zerocol > -1 {
        let mut z = zerocol as usize;
        while eqmap[z] != -1 && eqmap[z] != z as i64 {
            z = eqmap[z] as usize;
        }
        let tmpcol = z;
        while equalcolor[z] != z {
            z = equalcolor[z];
        }
        eqmap[prev_color] = z as i64;
        eqmap[tmpcol] = z as i64;

        if z < grcol {
            equalcolor[grcol] = z;
            group.color = z;
        } else if grcol < z {
            equalcolor[z] = grcol;
            eqmap[prev_color] = grcol as i64;
        }
    } else {
        eqmap[prev_color] = grcol as i64;
    }
}

/// CHI_THR constant for small exon filtering (header:58)
const CHI_THR: u64 = 50;
/// DROP constant for long read small exon filtering (constants)
const DROP: f64 = 0.5;
/// epsilon-like threshold when deciding whether a single fragment count remains.
const EPSILON: f64 = 1e-9;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct PairKey {
    name: Arc<str>,
    lo_start: u64,
    hi_start: u64,
    hi_tag: u32,
}

fn add_local_pair_pairing(out: &mut [Vec<(usize, f64)>], i: usize, j: usize, w: f64) {
    if i == j || i >= out.len() || j >= out.len() || w <= 0.0 {
        return;
    }
    let mut found = false;
    for (idx, cw) in out[i].iter_mut() {
        if *idx == j {
            *cw += w;
            found = true;
            break;
        }
    }
    if !found {
        out[i].push((j, w));
    }

    found = false;
    for (idx, cw) in out[j].iter_mut() {
        if *idx == i {
            *cw += w;
            found = true;
            break;
        }
    }
    if !found {
        out[j].push((i, w));
    }
}

fn build_local_pair_links(reads: &[BundleRead]) -> Vec<Vec<(usize, f64)>> {
    let mut out: Vec<Vec<(usize, f64)>> = vec![Vec::new(); reads.len()];

    // First, use any pre-computed pair_idx/pair_count that already map to this read vector.
    for (i, r) in reads.iter().enumerate() {
        for (k, &j) in r.pair_idx.iter().enumerate() {
            if j >= reads.len() || i >= j {
                continue;
            }
            let w = r.pair_count.get(k).copied().unwrap_or(0.0);
            if w > 0.0 {
                add_local_pair_pairing(&mut out, i, j, w);
            }
        }
    }

    // Then add standard mate links from read metadata to ensure local-index validity.
    let mut pending: HbHashMap<PairKey, Vec<(usize, f64)>> = HbHashMap::new();
    for (i, r) in reads.iter().enumerate() {
        let (Some(rid), Some(mrid), Some(mstart)) = (r.ref_id, r.mate_ref_id, r.mate_start) else {
            continue;
        };
        if rid != mrid || r.read_name.is_empty() {
            continue;
        }
        let key = PairKey {
            name: r.read_name.clone(),
            lo_start: r.ref_start.min(mstart),
            hi_start: r.ref_start.max(mstart),
            hi_tag: r.hi,
        };

        let maybe_prev = {
            let slot = pending.entry(key.clone()).or_default();
            let prev = slot.pop();
            if slot.is_empty() {
                pending.remove(&key);
            }
            prev
        };

        if let Some((j, _wprev)) = maybe_prev {
            let w = r.weight.max(0.0);
            add_local_pair_pairing(&mut out, i, j, w);
        } else {
            pending.entry(key).or_default().push((i, r.weight.max(0.0)));
        }
    }

    out
}

/// Determine if exon should be kept (keep exon logic from 
/// For long reads, also checks CHI_THR and DROP thresholds
fn should_keep_exon(
    seg_len: u64,
    exon_idx: usize,
    total_exons: usize,
    has_good_left_junc: bool,
    has_good_right_junc: bool,
    junction_support: u64,
    is_long_read: bool,
    read_len: u64,
) -> bool {
    // 
    // seg_len < junctionsupport || (longread && seg_len < CHI_THR && seg_len < DROP * read_len)
    let is_small = seg_len < junction_support
        || (is_long_read && seg_len < CHI_THR && (seg_len as f64) < DROP * (read_len as f64));

    if is_small {
        // First exon
        if exon_idx == 0 {
            if !has_good_right_junc {
                return false;
            }
        }
        // Last exon
        else if exon_idx == total_exons - 1 {
            if !has_good_left_junc {
                return false;
            }
        }
        // Middle exon
        else {
            if !has_good_left_junc && !has_good_right_junc {
                return false;
            }
        }
    }
    true
}

#[derive(Debug, Clone, Copy)]
struct ProcessReadResult {
    next_color: usize,
}

fn first_exon_valid(
    read: &BundleRead,
    junction_support: u64,
    killed_junctions: &HashSet<Junction>,
) -> bool {
    let first_len = read
        .exons
        .first()
        .map(|(s, e)| e.saturating_sub(*s).saturating_add(1))
        .unwrap_or(0);
    if first_len >= junction_support {
        return true;
    }
    if let Some(j) = read.junctions.first() {
        return !killed_junctions.contains(j);
    }
    false
}

fn last_exon_valid(
    read: &BundleRead,
    junction_support: u64,
    killed_junctions: &HashSet<Junction>,
) -> bool {
    let last_len = read
        .exons
        .last()
        .map(|(s, e)| e.saturating_sub(*s).saturating_add(1))
        .unwrap_or(0);
    if last_len >= junction_support {
        return true;
    }
    if let Some(j) = read.junctions.last() {
        return !killed_junctions.contains(j);
    }
    false
}

fn fragment_pair_joinable(
    left: &BundleRead,
    right: &BundleRead,
    localdist: u64,
    junction_support: u64,
    killed_junctions: &HashSet<Junction>,
) -> bool {
    left.ref_end < right.ref_start
        && left.ref_end.saturating_add(localdist) > right.ref_start
        && last_exon_valid(left, junction_support, killed_junctions)
        && first_exon_valid(right, junction_support, killed_junctions)
}

/// Main Original-style bundle building function
/// Returns bundles per strand
pub fn build_sub_bundles(
    reads: &[BundleRead],
    config: &RunConfig,
    good_junctions: &HashSet<Junction>,
    killed_junctions: &HashSet<Junction>,
) -> Result<Vec<SubBundleResult>> {
    if reads.is_empty() {
        return Ok(Vec::new());
    }
    let profile_3strand = std::env::var_os("RUSTLE_PROFILE_3STRAND").is_some();
    let t_start = Instant::now();
    if profile_3strand {
        eprintln!("PROFILE_3STRAND begin reads={}", reads.len());
    }
    // Initialize per-strand data structures (3 strands: -, ., +)
    let mut strand_data: [StrandBundles; 3] = [
        StrandBundles::new(),
        StrandBundles::new(),
        StrandBundles::new(),
    ];

    let mut next_color: usize = 0;

    // Union-find structures
    let mut eqcol: Vec<usize> = Vec::new();
    let mut merge: Vec<usize> = Vec::new();

    // Cross-strand color mapping for unknown strand reads 
    // pre-sizes both arrays to equalcolor.Count() before phase 1.
    // Initialized after eqcol is fully built (phase 1 sweep below).
    let mut eqnegcol: Vec<i64>; // color -> equivalent negative strand color
    let mut eqposcol: Vec<i64>; // color -> equivalent positive strand color

    // Track which global-group IDs each read belongs to (readgroup[n] stores group->grid IDs).
    let mut readgroups: Vec<Vec<usize>> = vec![Vec::new(); reads.len()];
    let pair_links = build_local_pair_links(reads);
    if profile_3strand {
        let mut max_links = 0usize;
        let mut max_i = 0usize;
        let mut total_links = 0usize;
        for (i, v) in pair_links.iter().enumerate() {
            total_links += v.len();
            if v.len() > max_links {
                max_links = v.len();
                max_i = i;
            }
        }
        eprintln!(
            "PROFILE_3STRAND pair_links total={} max={} at_read={}",
            total_links, max_links, max_i
        );
    }
    let mut grid_to_group: Vec<(usize, usize)> = Vec::new();
    let mut next_group_grid: usize = 0;

    // Sort reads by (start, end)
    let mut read_order: Vec<usize> = (0..reads.len()).collect();
    read_order.sort_unstable_by_key(|&i| {
        let r = &reads[i];
        (r.ref_start, r.ref_end)
    });

    // Process each read (add_read_to_group style: pair_count first, then single_count).
    for (ord_pos, &read_idx) in read_order.iter().enumerate() {
        let read = &reads[read_idx];
        if read.exons.is_empty() {
            continue;
        }
        if profile_3strand && ord_pos % 100 == 0 {
            eprintln!(
                "PROFILE_3STRAND ingest pos={} elapsed_ms={} groups=[{},{},{}] colors={}",
                ord_pos,
                t_start.elapsed().as_millis(),
                strand_data[0].groups.len(),
                strand_data[1].groups.len(),
                strand_data[2].groups.len(),
                eqcol.len()
            );
        }
        if profile_3strand && reads.len() > 800 && (430..=520).contains(&ord_pos) {
            eprintln!(
                "PROFILE_3STRAND read_start pos={} read_idx={} exons={} pair_links={} existing_readgroups={}",
                ord_pos,
                read_idx,
                read.exons.len(),
                pair_links[read_idx].len(),
                readgroups[read_idx].len()
            );
        }

        let mut usedcol_read: [i64; 3] = [-1, -1, -1];
        let mut single_count = read.weight.max(0.0);
        let localdist_pair = if config.long_reads {
            0u64
        } else {
            config.bundle_merge_dist.saturating_add(25)
        };

        for &(pair_idx, pair_cov_raw) in &pair_links[read_idx] {
            if pair_idx >= reads.len() {
                continue;
            }
            let pair_cov = pair_cov_raw.max(0.0);
            if pair_cov <= 0.0 {
                continue;
            }

            // merge_read_to_group early-return path when pair was already consumed
            // as a contiguous fragment in the earlier read.
            if pair_idx < read_idx
                && fragment_pair_joinable(
                    &reads[pair_idx],
                    read,
                    localdist_pair,
                    config.junction_support,
                    killed_junctions,
                )
            {
                single_count -= pair_cov;
                if single_count < 0.0 {
                    single_count = 0.0;
                }
                continue;
            }

            let mut sno = strand_idx(read.strand);
            let snop = strand_idx(reads[pair_idx].strand);
            if sno != snop {
                if sno == 1 {
                    sno = snop;
                } else if snop != 1 {
                    continue;
                }
            }

            single_count -= pair_cov;
            if single_count < 0.0 {
                single_count = 0.0;
            }

            let mut readcol = usedcol_read[sno];
            if pair_idx < read_idx {
                if let Some(&first_pair_grid) = readgroups[pair_idx].first() {
                    if merge.len() <= first_pair_grid {
                        let old_len = merge.len();
                        merge.resize(first_pair_grid + 1, 0);
                        for i in old_len..merge.len() {
                            merge[i] = i;
                        }
                    }
                    let root_grid = find_merge_root(first_pair_grid, &mut merge);
                    if let Some(pc) = group_color_from_grid(
                        root_grid,
                        &grid_to_group,
                        &mut strand_data,
                        &mut eqcol,
                    ) {
                        readcol = pc as i64;
                    }
                }
            } else if readcol < 0 {
                readcol = next_color as i64;
                usedcol_read[sno] = readcol;
                eqcol.push(next_color);
                next_color += 1;
            }

            if readcol < 0 {
                readcol = next_color as i64;
                usedcol_read[sno] = readcol;
                eqcol.push(next_color);
                next_color += 1;
            }

            let r1 = process_read_for_group(
                read_idx,
                read,
                reads,
                sno,
                readcol as usize,
                pair_cov,
                &mut strand_data[sno],
                &mut readgroups,
                &mut eqcol,
                &mut merge,
                &mut usedcol_read,
                &mut next_group_grid,
                &mut grid_to_group,
                next_color,
                pair_idx < read_idx,
                if pair_idx > read_idx {
                    Some(pair_idx)
                } else {
                    None
                },
                config,
                good_junctions,
                killed_junctions,
            )?;
            next_color = r1.next_color;
        }

        if single_count > EPSILON {
            let sno = strand_idx(read.strand);
            let mut readcol = usedcol_read[sno];
            if readcol < 0 {
                readcol = next_color as i64;
                usedcol_read[sno] = readcol;
                eqcol.push(next_color);
                next_color += 1;
            }

            let r = process_read_for_group(
                read_idx,
                read,
                reads,
                sno,
                readcol as usize,
                single_count,
                &mut strand_data[sno],
                &mut readgroups,
                &mut eqcol,
                &mut merge,
                &mut usedcol_read,
                &mut next_group_grid,
                &mut grid_to_group,
                next_color,
                false,
                None,
                config,
                good_junctions,
                killed_junctions,
            )?;
            next_color = r.next_color;
        }
        if profile_3strand && reads.len() > 800 && (430..=520).contains(&ord_pos) {
            eprintln!(
                "PROFILE_3STRAND read_done pos={} read_idx={} elapsed_ms={} groups=[{},{},{}] colors={}",
                ord_pos,
                read_idx,
                t_start.elapsed().as_millis(),
                strand_data[0].groups.len(),
                strand_data[1].groups.len(),
                strand_data[2].groups.len(),
                eqcol.len()
            );
        }
    }
    if profile_3strand {
        eprintln!(
            "PROFILE_3STRAND ingest_done elapsed_ms={} groups=[{},{},{}] colors={}",
            t_start.elapsed().as_millis(),
            strand_data[0].groups.len(),
            strand_data[1].groups.len(),
            strand_data[2].groups.len(),
            eqcol.len()
        );
        for sno in 0..3 {
            let (reachable, cycle_or_invalid) = count_reachable_groups(&strand_data[sno]);
            eprintln!(
                "PROFILE_3STRAND chain_pre sno={} startgroup={:?} reachable={} total={} broken={}",
                sno,
                strand_data[sno].startgroup,
                reachable,
                strand_data[sno].groups.len(),
                cycle_or_invalid
            );
        }
    }

    // only runs this pass when bundledist>0 (or guide-driven in short-read mode).
    // For long-read de novo with bundledist=0, this pass is skipped.
    let t_merge = Instant::now();
    if config.bundle_merge_dist > 0 {
        for sno in 0..3 {
            merge_close_groups(
                &mut strand_data[sno],
                &mut merge,
                &mut eqcol,
                config.bundle_merge_dist,
            )?;
        }
    }
    if profile_3strand {
        eprintln!(
            "PROFILE_3STRAND merge_close elapsed_ms={} total_ms={}",
            t_merge.elapsed().as_millis(),
            t_start.elapsed().as_millis()
        );
        for sno in 0..3 {
            let (reachable, cycle_or_invalid) = count_reachable_groups(&strand_data[sno]);
            eprintln!(
                "PROFILE_3STRAND chain_post sno={} startgroup={:?} reachable={} total={} broken={}",
                sno,
                strand_data[sno].startgroup,
                reachable,
                strand_data[sno].groups.len(),
                cycle_or_invalid
            );
        }
    }

    // Phase 1: Cross-strand color assignment (
    // Process all 3 strands together in coordinate order
    // long-read path (~16165+) uses strict overlap (no bundledist)
    // for cross-strand color projection of unstranded groups.
    let overlap_dist = if config.long_reads {
        0u64
    } else {
        config.bundle_merge_dist
    };
    eqnegcol = vec![-1; eqcol.len()];
    eqposcol = vec![-1; eqcol.len()];
    let mut currgroup: [Option<usize>; 3] = [None, None, None];
    let mut prevgroup: [Option<usize>; 3] = [None, None, None];

    // Initialize current groups
    for sno in 0..3 {
        currgroup[sno] = strand_data[sno].startgroup;
    }

    let t_phase1 = Instant::now();
    // Process groups in coordinate order.
    while currgroup[0].is_some() || currgroup[1].is_some() || currgroup[2].is_some() {
        let Some((nextgr, gidx)) = get_min_start(&currgroup, &strand_data) else {
            break;
        };
        if gidx >= strand_data[nextgr].groups.len() {
            currgroup[nextgr] = None;
            continue;
        }

        // Canonicalize current group's color root.
        let mut grcol = strand_data[nextgr].groups[gidx].color;
        while grcol < eqcol.len() && eqcol[grcol] != grcol {
            grcol = eqcol[grcol];
        }
        strand_data[nextgr].groups[gidx].color = grcol;

        if nextgr == 1 {
            // Unknown-strand group: project to overlapping negative/positive contexts.
            let ug_start = strand_data[1].groups[gidx].start;
            let ug_end = strand_data[1].groups[gidx].end;
            strand_data[1].groups[gidx].neg_prop = 0.0;

            // Previous negative overlap.
            if let Some(pg0) = prevgroup[0] {
                if pg0 < strand_data[0].groups.len() {
                    let (ng_start, ng_end, ng_cov, ng_color, ng_len) = {
                        let ng = &strand_data[0].groups[pg0];
                        (
                            ng.start,
                            ng.end,
                            ng.cov_sum,
                            ng.color,
                            ng.end.saturating_sub(ng.start).saturating_add(1).max(1),
                        )
                    };
                    if ug_start <= ng_end.saturating_add(overlap_dist) {
                        let ug_col = strand_data[1].groups[gidx].color;
                        // previous-group overlap passes raw prevgroup->color
                        // into set_strandcol (no pre-canonicalization at call site).
                        {
                            let ng = &mut strand_data[0].groups[pg0];
                            set_strandcol(ug_col, ng, ng_color, &mut eqnegcol, &mut eqcol);
                        }

                        let ov_start = ug_start.max(ng_start);
                        let ov_end = ug_end.min(ng_end);
                        let ov_len = ov_end.saturating_sub(ov_start).saturating_add(1).max(1);
                        strand_data[1].groups[gidx].neg_prop +=
                            ng_cov * (ov_len as f64) / (ng_len as f64);
                    }
                }
            }

            // Current negative overlaps (advance pointer).
            while let Some(cg0) = currgroup[0] {
                if cg0 >= strand_data[0].groups.len() {
                    currgroup[0] = None;
                    break;
                }
                let (ng_start, ng_end, ng_cov, ng_color, ng_next, ng_len) = {
                    let ng = &strand_data[0].groups[cg0];
                    (
                        ng.start,
                        ng.end,
                        ng.cov_sum,
                        ng.color,
                        ng.next_gr,
                        ng.end.saturating_sub(ng.start).saturating_add(1).max(1),
                    )
                };
                let overlaps = ug_start <= ng_end.saturating_add(overlap_dist)
                    && ng_start <= ug_end.saturating_add(overlap_dist);
                if !overlaps {
                    break;
                }

                let ng_col = find_color_root(ng_color, &eqcol);
                strand_data[0].groups[cg0].neg_prop = 1.0;
                let ug_col = strand_data[1].groups[gidx].color;
                {
                    let ng = &mut strand_data[0].groups[cg0];
                    ng.color = ng_col;
                    set_strandcol(ug_col, ng, ng_col, &mut eqnegcol, &mut eqcol);
                }
                let ov_start = ug_start.max(ng_start);
                let ov_end = ug_end.min(ng_end);
                let ov_len = ov_end.saturating_sub(ov_start).saturating_add(1).max(1);
                strand_data[1].groups[gidx].neg_prop += ng_cov * (ov_len as f64) / (ng_len as f64);

                prevgroup[0] = Some(cg0);
                currgroup[0] = ng_next;
            }

            // Positive overlaps.
            let mut pos_prop = 0.0f64;
            if let Some(pg2) = prevgroup[2] {
                if pg2 < strand_data[2].groups.len() {
                    let (pg_start, pg_end, pg_cov, pg_color, pg_len) = {
                        let pg = &strand_data[2].groups[pg2];
                        (
                            pg.start,
                            pg.end,
                            pg.cov_sum,
                            pg.color,
                            pg.end.saturating_sub(pg.start).saturating_add(1).max(1),
                        )
                    };
                    if ug_start <= pg_end.saturating_add(overlap_dist) {
                        let ug_col = strand_data[1].groups[gidx].color;
                        // previous-group overlap passes raw prevgroup->color
                        // into set_strandcol (no pre-canonicalization at call site).
                        {
                            let pg = &mut strand_data[2].groups[pg2];
                            set_strandcol(ug_col, pg, pg_color, &mut eqposcol, &mut eqcol);
                        }
                        if strand_data[1].groups[gidx].neg_prop > 0.0 {
                            let ov_start = ug_start.max(pg_start);
                            let ov_end = ug_end.min(pg_end);
                            let ov_len = ov_end.saturating_sub(ov_start).saturating_add(1).max(1);
                            pos_prop += pg_cov * (ov_len as f64) / (pg_len as f64);
                        }
                    }
                }
            }

            while let Some(cg2) = currgroup[2] {
                if cg2 >= strand_data[2].groups.len() {
                    currgroup[2] = None;
                    break;
                }
                let (pg_start, pg_end, pg_cov, pg_color, pg_next, pg_len) = {
                    let pg = &strand_data[2].groups[cg2];
                    (
                        pg.start,
                        pg.end,
                        pg.cov_sum,
                        pg.color,
                        pg.next_gr,
                        pg.end.saturating_sub(pg.start).saturating_add(1).max(1),
                    )
                };
                let overlaps = ug_start <= pg_end.saturating_add(overlap_dist)
                    && pg_start <= ug_end.saturating_add(overlap_dist);
                if !overlaps {
                    break;
                }

                let ug_col = strand_data[1].groups[gidx].color;
                let pg_col = find_color_root(pg_color, &eqcol);
                {
                    let pg = &mut strand_data[2].groups[cg2];
                    pg.color = pg_col;
                    set_strandcol(ug_col, pg, pg_col, &mut eqposcol, &mut eqcol);
                }
                if strand_data[1].groups[gidx].neg_prop > 0.0 {
                    let ov_start = ug_start.max(pg_start);
                    let ov_end = ug_end.min(pg_end);
                    let ov_len = ov_end.saturating_sub(ov_start).saturating_add(1).max(1);
                    pos_prop += pg_cov * (ov_len as f64) / (pg_len as f64);
                }

                prevgroup[2] = Some(cg2);
                currgroup[2] = pg_next;
            }

            if pos_prop > 0.0 {
                let np = strand_data[1].groups[gidx].neg_prop;
                strand_data[1].groups[gidx].neg_prop = np / (np + pos_prop);
            } else if strand_data[1].groups[gidx].neg_prop > 0.0 {
                strand_data[1].groups[gidx].neg_prop = 1.0;
            }
        } else if nextgr == 0 {
            strand_data[0].groups[gidx].neg_prop = 1.0;
        }

        prevgroup[nextgr] = Some(gidx);
        currgroup[nextgr] = strand_data[nextgr].groups[gidx].next_gr;
    }
    if profile_3strand {
        eprintln!(
            "PROFILE_3STRAND phase1 elapsed_ms={} total_ms={}",
            t_phase1.elapsed().as_millis(),
            t_start.elapsed().as_millis()
        );
    }
    if std::env::var_os("RUSTLE_DEBUG_BUNDLE").is_some() {
        for sno in 0..3 {
            eprintln!(
                "DEBUG_BG_GROUPS sno={} ngroups={}",
                sno,
                strand_data[sno].groups.len()
            );
        }
        let mut color_set: HbHashSet<usize> = HbHashSet::new();
        let mut rows: Vec<(usize, usize, u64, u64, usize, f64)> = Vec::new();
        for sno in 0..3 {
            for g in &strand_data[sno].groups {
                let root = find_color_root(g.color, &eqcol);
                color_set.insert(root);
                rows.push((sno, g.grid, g.start, g.end, root, g.cov_sum));
            }
        }
        eprintln!("DEBUG_BG_DISTINCT_COLORS total={}", color_set.len());
        rows.sort_unstable_by_key(|r| (r.0, r.2, r.3, r.1));
        for (sno, grid, start, end, color, cov) in rows {
            eprintln!(
                "DEBUG_BG_GRP sno={} grid={} start={} end={} color={} cov={:.1}",
                sno,
                grid,
                start.saturating_add(1),
                end,
                color,
                cov
            );
        }
        if std::env::var_os("RUSTLE_DEBUG_BUNDLE_DEEP").is_some() {
            eprintln!(
                "DEBUG_EQCOL_TABLE eqneg_size={} eqpos_size={} equalcolor_size={}",
                eqnegcol.len(),
                eqposcol.len(),
                eqcol.len()
            );
            let lim = eqcol.len().max(eqnegcol.len()).max(eqposcol.len());
            for c in 0..lim {
                let neg = if c < eqnegcol.len() { eqnegcol[c] } else { -1 };
                let pos = if c < eqposcol.len() { eqposcol[c] } else { -1 };
                if neg != -1 || pos != -1 {
                    eprintln!("  EQCOL c={} negcol={} poscol={}", c, neg, pos);
                }
            }
        }
    }

    // Phase 2: Create bundles in one coordinate-ordered 3-strand sweep.
    let t_phase2 = Instant::now();
    let out = create_bundles_global(
        &strand_data,
        &eqcol,
        &eqnegcol,
        &eqposcol,
        reads,
        &readgroups,
        &merge,
        config,
    );
    if profile_3strand {
        let rg_nonempty = readgroups.iter().filter(|g| !g.is_empty()).count();
        eprintln!(
            "PROFILE_3STRAND readgroups_nonempty={}/{}",
            rg_nonempty,
            reads.len()
        );
        if let Ok(ref bundles) = out {
            let mut mapped_reads: usize = 0;
            let mut seen = vec![false; reads.len()];
            for b in bundles {
                for (ri, ids) in b.read_to_bnodes.iter().enumerate() {
                    if !ids.is_empty() && ri < seen.len() && !seen[ri] {
                        seen[ri] = true;
                        mapped_reads += 1;
                    }
                }
            }
            eprintln!(
                "PROFILE_3STRAND output_bundles={} mapped_reads={}/{}",
                bundles.len(),
                mapped_reads,
                reads.len()
            );
        }
        eprintln!(
            "PROFILE_3STRAND phase2 elapsed_ms={} total_ms={}",
            t_phase2.elapsed().as_millis(),
            t_start.elapsed().as_millis()
        );
    }
    out
}

/// Process a single read for group assignment 
fn process_read_for_group(
    read_idx: usize,
    read: &BundleRead,
    reads: &[BundleRead],
    sno: usize,
    initial_color: usize,
    readcov: f64,
    strand_data: &mut StrandBundles,
    readgroups: &mut [Vec<usize>],
    eqcol: &mut Vec<usize>,
    merge: &mut Vec<usize>,
    usedcol: &mut [i64; 3],
    next_group_grid: &mut usize,
    grid_to_group: &mut Vec<(usize, usize)>,
    mut next_color: usize,
    pair_precedes: bool,
    pair_forward: Option<usize>,
    config: &RunConfig,
    good_junctions: &HashSet<Junction>,
    killed_junctions: &HashSet<Junction>,
) -> Result<ProcessReadResult> {
    let localdist = if config.long_reads {
        0u64
    } else {
        config.bundle_merge_dist.saturating_add(25)
    };

    let mut readcol = initial_color;
    let mut active_idx = read_idx;
    let mut active_read = read;
    let mut active_pair_precedes = pair_precedes;
    let mut pending_forward = pair_forward;

    // merge_read_to_group pointer flow
    let mut currgroup = strand_data.currgroup;

    if currgroup.is_some() {
        // set currgroup first
        let mut lastgroup: Option<usize> = None;
        while let Some(cg) = currgroup {
            if cg >= strand_data.groups.len() {
                currgroup = None;
                break;
            }
            if read.ref_start > strand_data.groups[cg].end {
                lastgroup = Some(cg);
                currgroup = strand_data.groups[cg].next_gr;
            } else {
                break;
            }
        }

        // Keep currgroup near this read's first segment.
        if currgroup.is_none()
            || currgroup
                .map(|cg| read.exons[0].1 < strand_data.groups[cg].start)
                .unwrap_or(false)
        {
            currgroup = lastgroup;
        }

        let mut first = true;
        let mut thisgroup = currgroup;
        let mut lastpushedgroup: Option<usize> = None;
        let mut ei = 0usize;
        let mut ncoord = active_read.exons.len();
        while ei < ncoord {
            let (seg_start, seg_end) = active_read.exons[ei];
            let seg_len = seg_end.saturating_sub(seg_start).saturating_add(1);
            let read_len = active_read.query_length.unwrap_or_else(|| {
                active_read
                    .exons
                    .iter()
                    .map(|(s, e)| e.saturating_sub(*s).saturating_add(1))
                    .sum()
            });

            let has_good_left = ei > 0 && {
                let j = Junction::new(active_read.exons[ei - 1].1, seg_start);
                good_junctions.contains(&j)
            };
            let has_good_right = ei + 1 < active_read.exons.len() && {
                let j = Junction::new(seg_end, active_read.exons[ei + 1].0);
                good_junctions.contains(&j)
            };
            let _has_bad_left = ei > 0 && {
                let j = Junction::new(active_read.exons[ei - 1].1, seg_start);
                killed_junctions.contains(&j)
            };

            // in the currgroup!=NULL branch,
            // keep-exon includes the long-read CHI_THR/DROP clause.
            let keep = should_keep_exon(
                seg_len,
                ei,
                active_read.exons.len(),
                has_good_left,
                has_good_right,
                config.junction_support,
                config.long_reads,
                read_len,
            );

            if !keep {
                if ei == 0 && lastgroup.is_some() {
                    currgroup = lastgroup;
                }
                ei += 1;
                continue;
            }

            while let Some(ti) = thisgroup {
                if ti >= strand_data.groups.len() {
                    thisgroup = None;
                    break;
                }
                if seg_start > strand_data.groups[ti].end {
                    lastgroup = Some(ti);
                    thisgroup = strand_data.groups[ti].next_gr;
                } else {
                    break;
                }
            }

            let overlaps = thisgroup
                .map(|ti| {
                    ti < strand_data.groups.len()
                        && seg_end.saturating_add(localdist) >= strand_data.groups[ti].start
                })
                .unwrap_or(false);

            if overlaps {
                let ti = thisgroup.unwrap();
                if ei == 0 {
                    let this_grid = strand_data.groups[ti].grid;
                    if merge.len() <= this_grid {
                        let old_len = merge.len();
                        merge.resize(this_grid + 1, 0);
                        for i in old_len..merge.len() {
                            merge[i] = i;
                        }
                    }
                    let root_grid = find_merge_root(this_grid, merge);
                    for &rg in &readgroups[active_idx] {
                        if rg < merge.len() && find_merge_root(rg, merge) == root_grid {
                            first = false;
                            break;
                        }
                    }
                }

                if active_read.nh > 1 {
                    strand_data.groups[ti].multi += readcov * (seg_len as f64);
                }
                if seg_start < strand_data.groups[ti].start {
                    strand_data.groups[ti].start = seg_start;
                }

                let mut nextgroup = strand_data.groups[ti].next_gr;
                while let Some(ni) = nextgroup {
                    if ni >= strand_data.groups.len() {
                        break;
                    }
                    if seg_end < strand_data.groups[ni].start {
                        break;
                    }
                    merge_fwd_groups(strand_data, ti, ni, merge, eqcol)?;
                    nextgroup = strand_data.groups[ti].next_gr;
                }
                if seg_end > strand_data.groups[ti].end {
                    strand_data.groups[ti].end = seg_end;
                }

                if active_read.ref_start == seg_start && ei == 0 {
                    strand_data.groups[ti].hardstart = true;
                }
                if active_read.ref_end == seg_end && ei + 1 == ncoord {
                    strand_data.groups[ti].hardend = true;
                }

                let group_color = find_color_root_mut(strand_data.groups[ti].color, eqcol);
                strand_data.groups[ti].color = group_color;

                // Pair color didn't reach this group at first exon.
                if ei == 0 && active_pair_precedes && group_color != readcol {
                    let mut rc = usedcol[sno];
                    if rc < 0 {
                        usedcol[sno] = next_color as i64;
                        rc = next_color as i64;
                        eqcol.push(next_color);
                        next_color += 1;
                    }
                    readcol = rc as usize;
                }

                if readcol != group_color {
                    // non-canonical junction (strand==0)
                    // prevents color propagation — readcol takes group's color instead of merging.
                    // Equivalent in rustle: junction not in good_junctions (!has_good_left).
                    if ei > 0 && !has_good_left {
                        readcol = group_color;
                    } else if readcol < group_color {
                        eqcol[group_color] = readcol;
                        strand_data.groups[ti].color = readcol;
                    } else {
                        eqcol[readcol] = group_color;
                        readcol = group_color;
                    }
                }

                let this_grid = strand_data.groups[ti].grid;
                if lastpushedgroup != Some(this_grid) {
                    if first {
                        readgroups[active_idx].push(this_grid);
                    }
                    lastpushedgroup = Some(this_grid);
                }
                strand_data.groups[ti].cov_sum += (seg_len as f64) * readcov;
            } else {
                // First exon split-pair: no overlap, color did not reach this position.
                if ei == 0 && active_pair_precedes {
                    let mut rc = usedcol[sno];
                    if rc < 0 {
                        usedcol[sno] = next_color as i64;
                        rc = next_color as i64;
                        eqcol.push(next_color);
                        next_color += 1;
                    }
                    readcol = rc as usize;
                } else if ei > 0 && !has_good_left {
                    // non-canonical junction → fresh color.
                    // Breaks color chain so groups across non-canonical introns are separate.
                    usedcol[sno] = next_color as i64;
                    readcol = next_color;
                    eqcol.push(next_color);
                    next_color += 1;
                }

                let ngroup = strand_data.groups.len();
                let ngroup_grid = *next_group_grid;
                let multi = if active_read.nh > 1 {
                    readcov * (seg_len as f64)
                } else {
                    0.0
                };
                let mut newgroup = CGroup::new(
                    seg_start,
                    seg_end,
                    readcol,
                    ngroup_grid,
                    (seg_len as f64) * readcov,
                    multi,
                );
                if active_read.ref_start == seg_start && ei == 0 {
                    newgroup.hardstart = true;
                }
                if active_read.ref_end == seg_end && ei + 1 == ncoord {
                    newgroup.hardend = true;
                }
                newgroup.next_gr = thisgroup;

                if let Some(li) = lastgroup {
                    strand_data.groups[li].next_gr = Some(ngroup);
                }

                strand_data.groups.push(newgroup);
                if merge.len() <= ngroup_grid {
                    let old_len = merge.len();
                    merge.resize(ngroup_grid + 1, 0);
                    for i in old_len..merge.len() {
                        merge[i] = i;
                    }
                } else {
                    merge[ngroup_grid] = ngroup_grid;
                }
                if grid_to_group.len() <= ngroup_grid {
                    grid_to_group.resize(ngroup_grid + 1, (0, 0));
                }
                grid_to_group[ngroup_grid] = (sno, ngroup);
                *next_group_grid += 1;

                readgroups[active_idx].push(ngroup_grid);
                lastpushedgroup = Some(ngroup_grid);
                thisgroup = Some(ngroup);
            }

            // pair-fragment continuation inside merge_read_to_group:
            // process forward pair in-place at the end of current read.
            if keep && ei + 1 == ncoord {
                if let Some(np) = pending_forward {
                    if np < reads.len()
                        && fragment_pair_joinable(
                            active_read,
                            &reads[np],
                            localdist,
                            config.junction_support,
                            killed_junctions,
                        )
                    {
                        active_idx = np;
                        active_read = &reads[np];
                        pending_forward = None;
                        active_pair_precedes = false;
                        ncoord = active_read.exons.len();
                        first = true;
                        if let Some(ti) = thisgroup {
                            if ti < strand_data.groups.len()
                                && active_read.ref_start > strand_data.groups[ti].end
                            {
                                let mut nextgroup = strand_data.groups[ti].next_gr;
                                while let Some(ni) = nextgroup {
                                    if ni >= strand_data.groups.len() {
                                        break;
                                    }
                                    if active_read.ref_start < strand_data.groups[ni].start {
                                        break;
                                    }
                                    merge_fwd_groups(strand_data, ti, ni, merge, eqcol)?;
                                    nextgroup = strand_data.groups[ti].next_gr;
                                }
                                if active_read.ref_start > strand_data.groups[ti].end {
                                    strand_data.groups[ti].end = active_read.ref_start;
                                }
                            }
                        }
                        lastpushedgroup = None;
                        ei = 0;
                        continue;
                    }
                }
            }

            ei += 1;
        }

        strand_data.currgroup = currgroup;
        if strand_data.startgroup.is_none() {
            strand_data.startgroup = currgroup;
        }
    } else {
        // currgroup == NULL: create new chain of groups for this read.
        if active_pair_precedes {
            let mut rc = usedcol[sno];
            if rc < 0 {
                usedcol[sno] = next_color as i64;
                rc = next_color as i64;
                eqcol.push(next_color);
                next_color += 1;
            }
            readcol = rc as usize;
        }

        let mut first_created_idx: Option<usize> = None;
        let mut lastgroup: Option<usize> = None;

        let mut ei = 0usize;
        let mut ncoord = active_read.exons.len();
        while ei < ncoord {
            let (seg_start, seg_end) = active_read.exons[ei];
            let seg_len = seg_end.saturating_sub(seg_start).saturating_add(1);
            let read_len = active_read.query_length.unwrap_or_else(|| {
                active_read
                    .exons
                    .iter()
                    .map(|(s, e)| e.saturating_sub(*s).saturating_add(1))
                    .sum()
            });
            let has_good_left = ei > 0 && {
                let j = Junction::new(active_read.exons[ei - 1].1, seg_start);
                good_junctions.contains(&j)
            };
            let has_good_right = ei + 1 < active_read.exons.len() && {
                let j = Junction::new(seg_end, active_read.exons[ei + 1].0);
                good_junctions.contains(&j)
            };
            let _has_bad_left = ei > 0 && {
                let j = Junction::new(active_read.exons[ei - 1].1, seg_start);
                killed_junctions.contains(&j)
            };
            let keep = should_keep_exon(
                seg_len,
                ei,
                active_read.exons.len(),
                has_good_left,
                has_good_right,
                config.junction_support,
                false,
                read_len,
            );

            if ei > 0 && !keep {
                ei += 1;
                continue;
            }

            if ei > 0 && !has_good_left {
                // non-canonical junction → fresh color.
                // When creating entirely new groups, a non-canonical junction breaks color chain.
                usedcol[sno] = next_color as i64;
                readcol = next_color;
                eqcol.push(next_color);
                next_color += 1;
            }

            let ngroup = strand_data.groups.len();
            let ngroup_grid = *next_group_grid;
            let multi = if active_read.nh > 1 {
                readcov * (seg_len as f64)
            } else {
                0.0
            };
            let mut newgroup = CGroup::new(
                seg_start,
                seg_end,
                readcol,
                ngroup_grid,
                (seg_len as f64) * readcov,
                multi,
            );
            if active_read.ref_start == seg_start && ei == 0 {
                newgroup.hardstart = true;
            }
            if active_read.ref_end == seg_end && ei + 1 == ncoord {
                newgroup.hardend = true;
            }

            strand_data.groups.push(newgroup);
            if let Some(li) = lastgroup {
                strand_data.groups[li].next_gr = Some(ngroup);
            } else {
                first_created_idx = Some(ngroup);
            }
            lastgroup = Some(ngroup);

            if keep {
                readgroups[active_idx].push(ngroup_grid);
            }
            if merge.len() <= ngroup_grid {
                let old_len = merge.len();
                merge.resize(ngroup_grid + 1, 0);
                for i in old_len..merge.len() {
                    merge[i] = i;
                }
            } else {
                merge[ngroup_grid] = ngroup_grid;
            }
            if grid_to_group.len() <= ngroup_grid {
                grid_to_group.resize(ngroup_grid + 1, (0, 0));
            }
            grid_to_group[ngroup_grid] = (sno, ngroup);
            *next_group_grid += 1;

            // pair-fragment continuation for currgroup == NULL branch.
            if keep && ei + 1 == ncoord {
                if let Some(np) = pending_forward {
                    if np < reads.len()
                        && fragment_pair_joinable(
                            active_read,
                            &reads[np],
                            localdist,
                            config.junction_support,
                            killed_junctions,
                        )
                    {
                        active_idx = np;
                        active_read = &reads[np];
                        pending_forward = None;
                        ncoord = active_read.exons.len();
                        if let Some(li) = lastgroup {
                            if li < strand_data.groups.len() {
                                let first_pair_end = active_read
                                    .exons
                                    .first()
                                    .map(|(_, e)| *e)
                                    .unwrap_or(active_read.ref_end);
                                if first_pair_end > strand_data.groups[li].end {
                                    strand_data.groups[li].end = first_pair_end;
                                    readgroups[active_idx].push(ngroup_grid);
                                    let first_pair_len = active_read
                                        .exons
                                        .first()
                                        .map(|(s, e)| e.saturating_sub(*s).saturating_add(1))
                                        .unwrap_or(0);
                                    strand_data.groups[li].cov_sum +=
                                        (first_pair_len as f64) * readcov;
                                }
                            }
                        }
                        ei = 1;
                        continue;
                    }
                }
            }

            ei += 1;
        }

        strand_data.currgroup = first_created_idx;
        if strand_data.startgroup.is_none() {
            strand_data.startgroup = first_created_idx;
        }
    }

    Ok(ProcessReadResult { next_color })
}

/// Merge two groups forward 
fn merge_fwd_groups(
    strand_data: &mut StrandBundles,
    g1_idx: usize,
    g2_idx: usize,
    merge: &mut Vec<usize>,
    eqcol: &mut Vec<usize>,
) -> Result<()> {
    if g1_idx >= strand_data.groups.len() || g2_idx >= strand_data.groups.len() {
        return Ok(());
    }

    let g1_grid = strand_data.groups[g1_idx].grid;
    let g2_grid = strand_data.groups[g2_idx].grid;

    // Ensure merge array is sized.
    let max_grid = g1_grid.max(g2_grid);
    if merge.len() <= max_grid {
        let old_len = merge.len();
        merge.resize(max_grid + 1, 0);
        for i in old_len..merge.len() {
            merge[i] = i;
        }
    }

    // :
    // - group1 is the merge target
    // - group2 grid maps directly to group1 grid
    // - color union does not use merge roots
    // - group1 end is set to group2 end
    strand_data.groups[g1_idx].end = strand_data.groups[g2_idx].end;

    let c1 = find_color_root_mut(strand_data.groups[g1_idx].color, eqcol);
    strand_data.groups[g1_idx].color = c1;
    let c2 = find_color_root_mut(strand_data.groups[g2_idx].color, eqcol);
    strand_data.groups[g2_idx].color = c2;
    let max_color = c1.max(c2);
    if eqcol.len() <= max_color {
        let old_len = eqcol.len();
        eqcol.resize(max_color + 1, 0);
        for i in old_len..eqcol.len() {
            eqcol[i] = i;
        }
    }
    if c1 < c2 {
        eqcol[c2] = c1;
    } else if c1 > c2 {
        eqcol[c1] = c2;
        strand_data.groups[g1_idx].color = c2;
    }

    strand_data.groups[g1_idx].cov_sum += strand_data.groups[g2_idx].cov_sum;
    strand_data.groups[g1_idx].next_gr = strand_data.groups[g2_idx].next_gr;
    merge[g2_grid] = g1_grid;
    strand_data.groups[g1_idx].multi += strand_data.groups[g2_idx].multi;

    Ok(())
}

/// Merge close groups pass 
fn merge_close_groups(
    strand_data: &mut StrandBundles,
    merge: &mut Vec<usize>,
    eqcol: &mut Vec<usize>,
    _bundledist: u64,
) -> Result<()> {
    // Ensure merge array is sized
    let max_grid = strand_data.groups.iter().map(|g| g.grid).max().unwrap_or(0);
    if merge.len() <= max_grid {
        for i in merge.len()..=max_grid {
            merge.push(i);
        }
    }

    // for each strand, merge groups within bundledist
    let mut lastgroup_idx = strand_data.startgroup;
    let mut outer_loops = 0;

    while let Some(li) = lastgroup_idx {
        outer_loops += 1;
        if outer_loops > strand_data.groups.len() * 2 {
            break;
        }
        let mut procgroup_idx = strand_data.groups[li].next_gr;
        let mut inner_loops = 0;

        while let Some(pi) = procgroup_idx {
            inner_loops += 1;
            if inner_loops > strand_data.groups.len() * 2 {
                break;
            }
            if pi >= strand_data.groups.len() {
                break;
            }

            let gap = strand_data.groups[pi]
                .start
                .saturating_sub(strand_data.groups[li].end);

            if gap <= _bundledist {
                merge_fwd_groups(strand_data, li, pi, merge, eqcol)?;
                procgroup_idx = strand_data.groups[li].next_gr;
                continue;
            }

            // Move to next procgroup
            procgroup_idx = strand_data.groups[pi].next_gr;
        }

        lastgroup_idx = strand_data.groups[li].next_gr;
    }

    Ok(())
}

fn add_group_to_tmp_bundle(nodes: &mut Vec<CBundlenode>, group: &CGroup, localdist: u64) -> usize {
    if nodes.is_empty() {
        let mut n = CBundlenode::new(group.start, group.end, group.cov_sum, 0);
        n.hardstart = group.hardstart;
        n.hardend = group.hardend;
        nodes.push(n);
        return 0;
    }
    let last_idx = nodes.len().saturating_sub(1);
    let last = nodes.last_mut().unwrap();
    if group.start > last.end.saturating_add(localdist) {
        let bid = nodes.len();
        let mut n = CBundlenode::new(group.start, group.end, group.cov_sum, bid);
        n.hardstart = group.hardstart;
        n.hardend = group.hardend;
        nodes.push(n);
        bid
    } else {
        if group.end > last.end {
            last.end = group.end;
        }
        last.cov += group.cov_sum;
        if group.hardstart && group.start == last.start {
            last.hardstart = true;
        }
        if group.hardend && group.end == last.end {
            last.hardend = true;
        }
        last_idx
    }
}

/// Original-style bundle creation pass
fn create_bundles_global(
    strand_data: &[StrandBundles; 3],
    eqcol: &[usize],
    eqnegcol: &[i64],
    eqposcol: &[i64],
    reads: &[BundleRead],
    readgroups: &[Vec<usize>],
    merge: &[usize],
    config: &RunConfig,
) -> Result<Vec<SubBundleResult>> {
    let localdist = config.bundle_merge_dist;
    let mut currgroup: [Option<usize>; 3] = [
        strand_data[0].startgroup,
        strand_data[1].startgroup,
        strand_data[2].startgroup,
    ];

    #[derive(Default, Clone)]
    struct TmpBundle {
        nodes: Vec<CBundlenode>,
        node_colors: Vec<usize>,
        color_root: Option<usize>,
        // group2bundle semantics: merge-root grid -> lastnodeid when group is added.
        rootgrid_to_node: HbHashMap<usize, usize>,
        // rprop support: merge-root grid -> strand proportion for this bundle.
        rootgrid_to_scale: HbHashMap<usize, f64>,
    }

    let mut bundles: [Vec<TmpBundle>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    let mut bundlecol: Vec<i64> = vec![-1; eqcol.len()];
    let trace_bundle = std::env::var_os("RUSTLE_DEBUG_BUNDLE_DEEP").is_some();
    let trace_bnode = std::env::var_os("RUSTLE_BNODE_TRACE").is_some();

    while currgroup[0].is_some() || currgroup[1].is_some() || currgroup[2].is_some() {
        let Some((nextgr, gidx)) = get_min_start(&currgroup, strand_data) else {
            break;
        };
        if gidx >= strand_data[nextgr].groups.len() {
            currgroup[nextgr] = None;
            continue;
        }
        let group = &strand_data[nextgr].groups[gidx];
        let mut grcol = find_color_root(group.color, eqcol);
        if grcol >= bundlecol.len() {
            bundlecol.resize(grcol + 1, -1);
        }

        let mut place_group =
            |target_sno: usize, color_root: usize, scale: f64| -> Option<(usize, usize)> {
                if scale <= 0.0 {
                    return None;
                }
                if color_root >= bundlecol.len() {
                    bundlecol.resize(color_root + 1, -1);
                }
                let was_new_bundle = bundlecol[color_root] < 0;
                let mut bno = bundlecol[color_root];
                if bno < 0 {
                    bno = bundles[target_sno].len() as i64;
                    bundlecol[color_root] = bno;
                    bundles[target_sno].push(TmpBundle::default());
                }
                let bidx = bno as usize;
                if bidx >= bundles[target_sno].len() {
                    bundles[target_sno].resize_with(bidx + 1, TmpBundle::default);
                }

                let scaled = CGroup {
                    start: group.start,
                    end: group.end,
                    color: color_root,
                    grid: group.grid,
                    cov_sum: group.cov_sum * scale,
                    multi: group.multi * scale,
                    neg_prop: group.neg_prop,
                    next_gr: group.next_gr,
                    hardstart: group.hardstart,
                    hardend: group.hardend,
                };
                let tb = &mut bundles[target_sno][bidx];
                if tb.color_root.is_none() {
                    tb.color_root = Some(color_root);
                }
                let nodes_before = tb.nodes.len();
                let node_idx = add_group_to_tmp_bundle(&mut tb.nodes, &scaled, localdist);
                let action = if was_new_bundle {
                    "create"
                } else if tb.nodes.len() > nodes_before {
                    "newnode"
                } else {
                    "extend"
                };
                if node_idx >= tb.node_colors.len() {
                    tb.node_colors.resize(node_idx + 1, color_root);
                }
                if tb.node_colors[node_idx] == 0 && color_root != 0 {
                    tb.node_colors[node_idx] = color_root;
                }
                let root_grid = find_merge_root_const(group.grid, merge);
                tb.rootgrid_to_node.insert(root_grid, node_idx);
                match tb.rootgrid_to_scale.get_mut(&root_grid) {
                    Some(s) => {
                        if scale > *s {
                            *s = scale;
                        }
                    }
                    None => {
                        tb.rootgrid_to_scale.insert(root_grid, scale);
                    }
                }
                if trace_bnode {
                    eprintln!(
                        "BNODE_ASSIGN sno={} bno={} grid={} start={} end={} cov={:.4} color={} action={} lastnodeid={}",
                        target_sno,
                        bidx,
                        group.grid,
                        group.start.saturating_add(1), // 1-based
                        group.end,
                        group.cov_sum * scale,
                        color_root,
                        action,
                        node_idx
                    );
                }
                Some((bidx, node_idx))
            };

        let has_neg_mapping = grcol < eqnegcol.len() && eqnegcol[grcol] >= 0;
        let has_pos_mapping = grcol < eqposcol.len() && eqposcol[grcol] >= 0;

        if nextgr == 0 || nextgr == 2 || (nextgr == 1 && !has_neg_mapping && !has_pos_mapping) {
            if let Some((bidx, node_idx)) = place_group(nextgr, grcol, 1.0) {
                if trace_bundle {
                    eprintln!(
                        "GROUP_TO_BNODE sno={} grid={} {}-{} color={} neg_prop={:.4} bno={} lastnodeid={}",
                        nextgr,
                        group.grid,
                        group.start,
                        group.end,
                        grcol,
                        group.neg_prop,
                        bidx,
                        node_idx
                    );
                }
            }
        } else {
            if has_neg_mapping {
                let negcol = find_color_root(eqnegcol[grcol] as usize, eqcol);
                let scale = if nextgr == 1 { group.neg_prop } else { 1.0 };
                if let Some((bidx, node_idx)) = place_group(0, negcol, scale) {
                    if trace_bundle {
                        eprintln!(
                            "GROUP_TO_BNODE_NEG grid={} {}-{} neg_prop={:.4} bno={} lastnodeid={}",
                            group.grid, group.start, group.end, group.neg_prop, bidx, node_idx
                        );
                    }
                }
            }
            if has_pos_mapping {
                let poscol = find_color_root(eqposcol[grcol] as usize, eqcol);
                let scale = if nextgr == 1 {
                    1.0 - group.neg_prop
                } else {
                    1.0
                };
                if let Some((bidx, node_idx)) = place_group(2, poscol, scale) {
                    if trace_bundle {
                        eprintln!(
                            "GROUP_TO_BNODE_POS grid={} {}-{} neg_prop={:.4} bno={} lastnodeid={}",
                            group.grid, group.start, group.end, group.neg_prop, bidx, node_idx
                        );
                    }
                }
            }
        }

        currgroup[nextgr] = group.next_gr;
        grcol = find_color_root(grcol, eqcol);
        let _ = grcol;
    }

    let mut results = Vec::new();
    let profile_detail = std::env::var_os("RUSTLE_PROFILE_3STRAND_DETAIL").is_some();

    // Pre-index read membership by merge-root group id so bundle assembly
    // does not need to rescan every readgroup for every bundle.
    let mut root_to_reads: HbHashMap<usize, Vec<usize>> = HbHashMap::new();
    for (ri, rgs) in readgroups.iter().enumerate() {
        if rgs.is_empty() {
            continue;
        }
        let mut roots_u32 = RoaringBitmap::new();
        let mut roots_overflow: Option<HbHashSet<usize>> = None;
        for &rg in rgs {
            let root = find_merge_root_const(rg, merge);
            if root <= u32::MAX as usize {
                roots_u32.insert(root as u32);
            } else {
                roots_overflow.get_or_insert_with(Default::default).insert(root);
            }
        }
        for root in roots_u32.iter().map(|v| v as usize) {
            root_to_reads.entry(root).or_default().push(ri);
        }
        if let Some(extra_roots) = roots_overflow {
            for root in extra_roots {
                root_to_reads.entry(root).or_default().push(ri);
            }
        }
    }

    for reads in root_to_reads.values_mut() {
        // Keep stable-ish local ordering for deterministic downstream diagnostics.
        reads.sort_unstable();
        reads.dedup();
    }

    for sno in 0..3 {
        let strand = match sno {
            0 => '-',
            2 => '+',
            _ => '.',
        };
        for tb in &bundles[sno] {
            if tb.nodes.is_empty() {
                continue;
            }
            if profile_detail && reads.len() >= 500 {
                let mut mapped_ids: HbHashSet<usize> = HbHashSet::new();
                let mut intersect_keys = 0usize;
                for &root_rg in tb.rootgrid_to_node.keys() {
                    if let Some(read_ids) = root_to_reads.get(&root_rg) {
                        intersect_keys += 1;
                        for &ri in read_ids {
                            mapped_ids.insert(ri);
                        }
                    }
                }
                eprintln!(
                    "PROFILE_3STRAND_DETAIL strand={} reads={} root_to_reads_keys={} tb_rootnode_keys={} tb_rootscale_keys={} intersect_keys={} mapped_ids={}",
                    strand,
                    reads.len(),
                    root_to_reads.len(),
                    tb.rootgrid_to_node.len(),
                    tb.rootgrid_to_scale.len(),
                    intersect_keys,
                    mapped_ids.len()
                );
            }
            let mut nodes = tb.nodes.clone();
            for (i, n) in nodes.iter_mut().enumerate() {
                n.bid = i;
                n.next = None;
            }
            let mut head = nodes[0].clone();
            let mut cur = &mut head;
            for n in nodes.iter().skip(1) {
                cur.next = Some(Box::new(n.clone()));
                cur = cur.next.as_deref_mut().unwrap();
            }
            let start = nodes.first().map(|n| n.start).unwrap_or(0);
            let end = nodes.last().map(|n| n.end).unwrap_or(0);
            results.push(SubBundleResult {
                strand,
                start,
                end,
                bnode_head: Some(head),
                read_to_bnodes: {
                    let mut out: Vec<Vec<usize>> = vec![Vec::new(); reads.len()];
                    for (&root_rg, &nid) in tb.rootgrid_to_node.iter() {
                        if let Some(read_ids) = root_to_reads.get(&root_rg) {
                            for &ri in read_ids {
                                out[ri].push(nid);
                            }
                        }
                    }
                    for dst in out.iter_mut() {
                        dst.sort_unstable();
                        dst.dedup();
                    }
                    out
                },
                bnode_colors: if tb.node_colors.len() == nodes.len() {
                    tb.node_colors.clone()
                } else {
                    vec![tb.color_root.unwrap_or(0); nodes.len()]
                },
                read_scale: {
                    let mut out = vec![0.0f64; reads.len()];
                    for (ri, rgs) in readgroups.iter().enumerate() {
                        if rgs.is_empty() {
                            continue;
                        }
                        let mut sum = 0.0f64;
                        let mut n = 0usize;
                        for &rg in rgs {
                            let root_rg = find_merge_root_const(rg, merge);
                            if let Some(&s) = tb.rootgrid_to_scale.get(&root_rg) {
                                sum += s;
                                n += 1;
                            }
                        }
                        if n > 0 {
                            out[ri] = (sum / n as f64).clamp(0.0, 1.0);
                        }
                    }
                    out
                },
            });
        }
    }
    results.sort_unstable_by_key(|r| r.start);
    Ok(results)
}
