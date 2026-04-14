//! Junction filtering: keep junctions above threshold; weak and nearby-weak filters.

use crate::junction_correction::correct_bundle_junctions_higherr;
use crate::types::{
    DetHashMap as HashMap, DetHashSet as HashSet, Junction, JunctionStat, JunctionStats,
};

/// Filter to junctions with mrcount >= min_junction_reads.
/// C++ parity: only keep junctions with non-negative `nreads_good` at threshold
/// and positive mrcount (actual read support).  Junctions with mrcount=0 are
/// counting artifacts — they have nreads_good > 0 from anchor-based accounting but
/// zero actual read coverage.  These create spurious graph node boundaries and
/// generate excess transfrags (7.7x over the reference assembler in the worst loci).
pub fn filter_junctions(stats: &JunctionStats, min_reads: f64) -> Vec<Junction> {
    stats
        .iter()
        // IMPORTANT: `mm` can be used as a temporary BAD_MM_NEG marker in higherr/long-read
        // parity paths (C++ uses `mm=-1` as a marker, not an immediate deletion).
        // Dropping junctions purely because `mm < 0` can remove real junctions prematurely.
        .filter(|(_, s)| s.strand != Some(0) && s.nreads_good >= min_reads && s.mrcount > 0.0)
        .map(|(j, _)| *j)
        .collect()
}


/// Remove junctions with mrcount < min_reads from stats (in-place). Returns set of removed junctions.
/// Demoted left/right junctions must drop here before graph construction.
pub fn filter_weak_junctions(
    stats: &mut JunctionStats,
    min_reads: f64,
    verbose: bool,
) -> HashSet<Junction> {
    let to_remove: Vec<Junction> = stats
        .iter()
        // See `filter_junctions`: `mm < 0` is not an unconditional removal criterion.
        .filter(|(_, s)| s.nreads_good < min_reads || s.mrcount < 0.0)
        .map(|(j, _)| *j)
        .collect();
    let n = to_remove.len();
    for j in &to_remove {
        stats.remove(j);
    }
    if verbose && n > 0 {
        eprintln!(
            "    Weak junction filter removed {} junctions (nreads_good < {})",
            n, min_reads
        );
    }
    to_remove.into_iter().collect()
}

const NEARBY_WEAK_WINDOW: u64 = 25;
/// Minimum read-count fraction relative to the strongest nearby junction.
/// C++ default is ~5% but long-read aligners produce shifted splice site calls
/// (same junction offset by 1-30bp) that create redundant graph node boundaries
/// and 3.3x transfrag proliferation.  10% threshold removes the weakest shifted
/// copies while preserving genuine alternative splice sites.
const NEARBY_WEAK_RATIO: f64 = 0.10;

/// Remove junctions that are within window bp (donor and acceptor) of a much stronger junction
/// and have reads ratio < ratio_threshold vs that stronger one (reference-style nearby weak filter).
pub fn filter_nearby_weak_junctions(
    stats: &JunctionStats,
    window: u64,
    ratio_threshold: f64,
    verbose: bool,
) -> HashSet<Junction> {
    if stats.is_empty() {
        return Default::default();
    }
    let junction_list: Vec<Junction> = {
        let mut v: Vec<_> = stats.keys().copied().collect();
        v.sort_by_key(|j| (j.donor, j.acceptor));
        v
    };
    let n = junction_list.len();
    let reads: Vec<f64> = junction_list
        .iter()
        .map(|j| stats.get(j).map(|s| s.mrcount).unwrap_or(0.0))
        .collect();

    let mut to_remove: HashSet<Junction> = Default::default();
    for i in 0..n {
        let ji = junction_list[i];
        let reads_i = reads[i];
        let mut max_nearby_reads = reads_i;
        for j in 0..n {
            if i == j {
                continue;
            }
            let jj = junction_list[j];
            let donor_ok = ji.donor.abs_diff(jj.donor) <= window;
            let acceptor_ok = ji.acceptor.abs_diff(jj.acceptor) <= window;
            if donor_ok && acceptor_ok {
                if reads[j] > max_nearby_reads {
                    max_nearby_reads = reads[j];
                }
            }
        }
        if max_nearby_reads > reads_i && reads_i / max_nearby_reads < ratio_threshold {
            to_remove.insert(ji);
        }
    }

    if verbose && !to_remove.is_empty() {
        eprintln!(
            "    Nearby weak junction filter: removed {} junctions (within {}bp of stronger, ratio < {})",
            to_remove.len(),
            window,
            ratio_threshold
        );
    }
    to_remove
}

/// Canonicalize junctions: merge those within tolerance bp (donor and acceptor);
/// keep the junction with highest mrcount per cluster (reference-style boundary noise reduction).
pub fn canonicalize_junctions(
    stats: &JunctionStats,
    tolerance: u64,
    verbose: bool,
) -> JunctionStats {
    if stats.is_empty() {
        return Default::default();
    }
    let mut junction_list: Vec<(Junction, JunctionStat)> =
        stats.iter().map(|(j, s)| (*j, s.clone())).collect();
    junction_list.sort_by_key(|(j, _)| (j.donor, j.acceptor));

    let mut visited: HashSet<Junction> = Default::default();
    let mut canonical_map: HashMap<Junction, Junction> = Default::default();
    let n = junction_list.len();

    for i in 0..n {
        let (ji, _) = junction_list[i];
        if visited.contains(&ji) {
            continue;
        }
        let mut cluster: Vec<(Junction, JunctionStat)> = vec![junction_list[i].clone()];
        visited.insert(ji);

        for j in (i + 1)..n {
            let (jj, _) = junction_list[j];
            if jj.donor.saturating_sub(ji.donor) > tolerance {
                break;
            }
            if ji.acceptor.abs_diff(jj.acceptor) <= tolerance {
                cluster.push(junction_list[j].clone());
                visited.insert(jj);
            }
        }

        let canonical = cluster
            .iter()
            .max_by(|a, b| {
                // Priority 1: Guide match
                if a.1.guide_match != b.1.guide_match {
                    return a.1.guide_match.cmp(&b.1.guide_match);
                }
                // Priority 2: Canonical splice sites (GT-AG etc)
                let a_canonical = a.1.consleft == 1 && a.1.consright == 1;
                let b_canonical = b.1.consleft == 1 && b.1.consright == 1;
                if a_canonical != b_canonical {
                    return a_canonical.cmp(&b_canonical);
                }
                // Priority 3: Coverage (mrcount)
                a.1.mrcount
                    .partial_cmp(&b.1.mrcount)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(j, _)| *j)
            .unwrap_or(ji);
        for (j, _) in &cluster {
            canonical_map.insert(*j, canonical);
        }
    }

    let mut new_stats: JunctionStats = Default::default();
    for (j, stat) in stats.iter() {
        let canonical = canonical_map.get(j).copied().unwrap_or(*j);
        let entry = new_stats.entry(canonical).or_insert_with(Default::default);
        entry.mrcount += stat.mrcount;
        entry.nreads_good += stat.nreads_good;
        entry.rcount += stat.rcount;
        entry.nm += stat.nm;
        entry.mm += stat.mm;
        entry.leftsupport += stat.leftsupport;
        entry.rightsupport += stat.rightsupport;
        if entry.strand.is_none() {
            entry.strand = stat.strand;
        }
    }

    let merged = canonical_map.iter().filter(|(a, b)| a != b).count();
    if verbose && merged > 0 {
        eprintln!(
            "    Junction canonicalization: merged {} boundary variants (tolerance={}bp)",
            merged, tolerance
        );
    }
    new_stats
}

/// Coalesce junctions: merge junctions within tolerance (donor and acceptor) by combining their support.
/// This is more aggressive than canonicalization as it actively collapses similar junctions 
/// to reduce noise from alignment artifacts (small deletions near splice sites).
pub fn coalesce_junctions(
    stats: JunctionStats,
    tolerance: u64,
    verbose: bool,
) -> JunctionStats {
    if stats.is_empty() || tolerance == 0 {
        return stats;
    }
    let mut junction_list: Vec<(Junction, JunctionStat)> =
        stats.into_iter().map(|(j, s)| (j, s)).collect();
    junction_list.sort_by_key(|(j, _)| (j.donor, j.acceptor));

    let mut visited: HashSet<usize> = Default::default();
    let mut new_stats: JunctionStats = Default::default();
    let n = junction_list.len();

    for i in 0..n {
        if visited.contains(&i) {
            continue;
        }
        let (mut ji, mut si) = junction_list[i].clone();
        visited.insert(i);

        let mut cluster_indices = vec![i];
        for j in (i + 1)..n {
            if visited.contains(&j) {
                continue;
            }
            let (jj, _) = &junction_list[j];
            if jj.donor.saturating_sub(ji.donor) > tolerance {
                break;
            }
            if ji.acceptor.abs_diff(jj.acceptor) <= tolerance {
                cluster_indices.push(j);
                visited.insert(j);
            }
        }

        if cluster_indices.len() > 1 {
            // Find best representative in cluster
            let mut best_idx = i;
            for &idx in &cluster_indices {
                let (_j_cand, s_cand) = &junction_list[idx];
                let (_j_best, s_best) = &junction_list[best_idx];
                
                let cand_score = (if s_cand.guide_match { 2 } else { 0 }) + 
                                (if s_cand.consleft == 1 && s_cand.consright == 1 { 1 } else { 0 });
                let best_score = (if s_best.guide_match { 2 } else { 0 }) + 
                                (if s_best.consleft == 1 && s_best.consright == 1 { 1 } else { 0 });
                
                if cand_score > best_score || (cand_score == best_score && s_cand.mrcount > s_best.mrcount) {
                    best_idx = idx;
                }
            }
            
            ji = junction_list[best_idx].0;
            // Sum support into the best representative
            let mut total_stat = junction_list[best_idx].1.clone();
            for &idx in &cluster_indices {
                if idx == best_idx { continue; }
                let (_, s) = &junction_list[idx];
                total_stat.mrcount += s.mrcount;
                total_stat.nreads_good += s.nreads_good;
                total_stat.rcount += s.rcount;
                total_stat.nm += s.nm;
                total_stat.mm += s.mm;
                total_stat.leftsupport += s.leftsupport;
                total_stat.rightsupport += s.rightsupport;
            }
            si = total_stat;
        }
        
        new_stats.insert(ji, si);
    }

    let merged = n - new_stats.len();
    if verbose && merged > 0 {
        eprintln!(
            "    Junction coalescing: merged {} near-identical junctions (tolerance={}bp)",
            merged, tolerance
        );
    }
    new_stats
}

/// Per-splice-site isofrac filter (C++ ref): remove junctions whose mrcount is < isofrac% of
/// the total mrcount at that donor site. isofrac is in percent (e.g. 1.0 = 1%).
pub fn filter_per_splice_site_isofrac(
    stats: &JunctionStats,
    isofrac_percent: f64,
    verbose: bool,
) -> HashSet<Junction> {
    if stats.is_empty() || isofrac_percent <= 0.0 {
        return Default::default();
    }
    let threshold = isofrac_percent / 100.0;
    let mut donor_sum: HashMap<u64, f64> = Default::default();
    for (j, s) in stats.iter() {
        *donor_sum.entry(j.donor).or_insert(0.0) += s.mrcount;
    }
    let to_remove: Vec<Junction> = stats
        .iter()
        .filter(|(j, s)| {
            let total = donor_sum.get(&j.donor).copied().unwrap_or(0.0);
            total > 0.0 && (s.mrcount / total) < threshold
        })
        .map(|(j, _)| *j)
        .collect();
    if verbose && !to_remove.is_empty() {
        eprintln!(
            "    Per-splice-site filter removed {} junctions (isofrac={}%)",
            to_remove.len(),
            isofrac_percent
        );
    }
    to_remove.into_iter().collect()
}

/// Apply both weak and nearby-weak filters; returns updated stats (owned).
pub fn apply_junction_filters(
    mut stats: JunctionStats,
    min_junction_reads: f64,
    verbose: bool,
) -> JunctionStats {
    // strand=0 marks killed junctions (good_junc/good_merge_junc parity); remove from active set.
    let killed: Vec<Junction> = stats
        .iter()
        .filter(|(_, s)| s.strand == Some(0))
        .map(|(j, _)| *j)
        .collect();
    for j in &killed {
        stats.remove(j);
    }
    if verbose && !killed.is_empty() {
        eprintln!(
            "    Removed {} killed junction(s) before weak filters",
            killed.len()
        );
    }
    filter_weak_junctions(&mut stats, min_junction_reads, verbose);
    let to_remove =
        filter_nearby_weak_junctions(&stats, NEARBY_WEAK_WINDOW, NEARBY_WEAK_RATIO, verbose);
    for j in to_remove {
        stats.remove(&j);
    }
    stats
}

/// Apply filters, optional per-splice-site isofrac, then canonicalize (full junction preprocessing).
pub fn apply_junction_filters_and_canonicalize(
    stats: JunctionStats,
    min_junction_reads: f64,
    canonical_tolerance: u64,
    per_splice_site_isofrac: f64,
    verbose: bool,
) -> JunctionStats {
    let mut stats = coalesce_junctions(stats, canonical_tolerance, verbose);
    stats = apply_junction_filters(stats, min_junction_reads, verbose);
    if per_splice_site_isofrac > 0.0 {
        let to_remove = filter_per_splice_site_isofrac(&stats, per_splice_site_isofrac, verbose);
        for j in to_remove {
            stats.remove(&j);
        }
    }
    canonicalize_junctions(&stats, canonical_tolerance, verbose)
}

/// Run correction-only stage and return correction map.
/// C++ only corrects bad junctions via higherr (C++ reference).
pub fn correct_junctions_with_map(
    stats: &JunctionStats,
    verbose: bool,
) -> (JunctionStats, HashMap<Junction, Junction>) {
    let (stats1, higherr_map) = correct_bundle_junctions_higherr(stats, 25, 0.9, verbose);
    if verbose && !higherr_map.is_empty() {
        eprintln!("    total correction remap entries: {}", higherr_map.len());
    }
    (stats1, higherr_map)
}
