//! Killed junctions and witness logic (good_junc strand=0; original algorithm V99 killed_junction_pairs,
//! killed_junction_orphan; has_lr_witness_two_splices; witness left/right = leftsupport/rightsupport).

use crate::bpcov::{BpcovStranded, BPCOV_STRAND_ALL, BPCOV_STRAND_MINUS, BPCOV_STRAND_PLUS};
use crate::graph::{Graph, GraphTransfrag};
use crate::types::{
    CJunction, DetHashMap as HashMap, DetHashSet as HashSet, Junction, JunctionStats,
};

#[inline]
fn goodjunc_trace_active() -> bool {
    std::env::var_os("RUSTLE_GOODJUNC_DEBUG").is_some()
        || std::env::var_os("RUSTLE_TRACE_LOG_STYLE").is_some()
}

#[inline]
fn trace_junction_acceptor(jn: &Junction) -> u64 {
    jn.acceptor.saturating_add(1)
}

#[inline]
fn trace_cjunction_acceptor(cj: &CJunction) -> u64 {
    cj.end.saturating_add(1)
}

/// Junctions with strand == 0 (good_junc), map-based variant.
/// Key = (donor, acceptor) same as Junction.
pub fn compute_killed_junction_pairs_stats(junction_stats: &JunctionStats) -> HashSet<Junction> {
    junction_stats
        .iter()
        .filter(|(_, stat)| {
            stat.strand == Some(0) || stat.mm < 0.0 || stat.mrcount < 0.0 || stat.nreads_good < 0.0
        })
        .map(|(j, _)| *j)
        .collect()
}

/// Second-stage good_junc validation for read-loop use (build_graphs).
/// 
/// good_junc is called TWICE:
/// 1. During initial junction filtering (before build_graphs)
/// 2. AGAIN inside the read loop during build_graphs
/// 
/// This second-stage validation catches junctions that passed initial filtering
/// but fail when re-evaluated with current coverage context.
/// 
/// Unlike the full good_junc, this only checks coverage witness for edge cases
/// where the junction might have passed initial filtering due to aggregated
/// statistics but fails when checked against the actual read coverage context.
/// 
/// Returns true if the junction passes validation, false if it should be killed.
pub fn good_junc_second_stage(
    junction: Junction,
    strand: i8,
    stats: &JunctionStats,
    bpcov: &BpcovStranded,
    refstart: u64,
    junction_thr: f64,
    _long_intron: u64,
    longreads: bool,
) -> bool {
    let Some(st) = stats.get(&junction) else {
        return true; // No stats = pass through (will be caught elsewhere)
    };
    
    // Already marked as killed - pass through
    if st.strand == Some(0) || st.mm < 0.0 || st.mrcount < 0.0 || st.nreads_good < 0.0 {
        return false;
    }
    
    // Guide-matched junctions always pass
    if st.guide_match {
        return true;
    }
    
    // Check strand consistency
    if st.strand != Some(strand) {
        return false;
    }
    
    const ERROR_PERC: f64 = 0.1;
    
    // Only check coverage witness for "all-bad" junctions (nm ~ nreads)
    // These are the junctions most likely to be artifacts
    let all_bad = st.nm > 0.0 && (st.nm - st.mrcount).abs() < 1.0;

    // Also check junctions with marginal support
    let marginal = st.nreads_good < 1.25 * junction_thr;

    if !all_bad && !marginal {
        return true;
    }

    let mismatch = all_bad;

    // Coverage witness check for edge-case junctions
    let mut mult = 1.0 / ERROR_PERC;
    if longreads {
        mult /= ERROR_PERC;
    }
    let bw = 5u64;
    let sno = if strand > 0 {
        BPCOV_STRAND_PLUS
    } else {
        BPCOV_STRAND_MINUS
    };
    
    // Donor side check
    let j_left = junction.donor.saturating_sub(refstart) as i64;
    if j_left >= bw as i64 && (j_left + bw as i64 + 1) < bpcov.plus.cov.len() as i64 {
        let a = (j_left - bw as i64 + 1) as usize;
        let b = (j_left + 1) as usize;
        let c = (j_left + 1) as usize;
        let d = (j_left + bw as i64 + 1) as usize;
        let lleftcov = cov_sign_window(bpcov, sno, a, b, longreads);
        let lrightcov = cov_sign_window(bpcov, sno, c, d, longreads);
        
        if lleftcov > 1.0 / ERROR_PERC
            && st.leftsupport * mult < ERROR_PERC * lleftcov
            && (mismatch || lrightcov > lleftcov * (1.0 - ERROR_PERC))
        {
            return false;
        }
    }
    
    // Acceptor side check - only for edge cases
    let j_right = junction.acceptor.saturating_sub(refstart).saturating_sub(1) as i64;
    if j_right >= bw as i64 && (j_right + bw as i64 + 1) < bpcov.plus.cov.len() as i64 {
        let a = (j_right - bw as i64 + 1) as usize;
        let b = (j_right + 1) as usize;
        let c = (j_right + 1) as usize;
        let d = (j_right + bw as i64 + 1) as usize;
        let rleftcov = cov_sign_window(bpcov, sno, a, b, longreads);
        let rrightcov = cov_sign_window(bpcov, sno, c, d, longreads);
        
        if rrightcov > 1.0 / ERROR_PERC
            && st.rightsupport * mult < ERROR_PERC * rrightcov
            && (mismatch || rleftcov > rrightcov * (1.0 - ERROR_PERC))
        {
            return false;
        }
    }
    
    // Step 6: don't keep if splice fraction is extremely low
    if st.mrcount * 10.0 < ERROR_PERC * st.leftsupport
        || st.mrcount * 10.0 < ERROR_PERC * st.rightsupport
    {
        return false;
    }

    true
}



/// True if the path uses no killed junction (has non–killed-junction only). Used to require
/// transcripts to use only good junctions when validating.
pub fn has_non_kjo(
    path_node_ids: &[usize],
    graph: &Graph,
    killed_junction_pairs: &HashSet<Junction>,
) -> bool {
    if path_node_ids.len() < 2 || killed_junction_pairs.is_empty() {
        return true;
    }
    let nodes = &graph.nodes;
    for i in 0..path_node_ids.len() - 1 {
        let (a, b) = (path_node_ids[i], path_node_ids[i + 1]);
        if a >= nodes.len() || b >= nodes.len() {
            continue;
        }
        let donor = nodes[a].end; // 0-based: end of exon segment
        let acceptor = nodes[b].start;
        let j = Junction::new(donor, acceptor);
        if killed_junction_pairs.contains(&j) {
            return false;
        }
    }
    true
}

/// Witness left: support on donor (left) side of splice (leftsupport).
#[inline]
pub fn witness_left(junction: Junction, junction_stats: &JunctionStats) -> f64 {
    junction_stats
        .get(&junction)
        .map(|s| s.leftsupport)
        .unwrap_or(0.0)
}

/// Witness right: support on acceptor (right) side of splice (rightsupport).
#[inline]
pub fn witness_right(junction: Junction, junction_stats: &JunctionStats) -> f64 {
    junction_stats
        .get(&junction)
        .map(|s| s.rightsupport)
        .unwrap_or(0.0)
}

/// Extract junctions (donor, acceptor) from a path of node IDs using graph coordinates.
fn path_junctions(path: &[usize], graph: &Graph) -> Vec<Junction> {
    let mut out = Vec::new();
    let nodes = &graph.nodes;
    for i in 0..path.len().saturating_sub(1) {
        let (a, b) = (path[i], path[i + 1]);
        if a >= nodes.len() || b >= nodes.len() {
            continue;
        }
        let donor = nodes[a].end;
        let acceptor = nodes[b].start;
        if donor != acceptor {
            out.push(Junction::new(donor, acceptor));
        }
    }
    out
}

/// At least one long-read transfrag contains both junctions (in order) within tolerance.
/// has_lr_witness_two_splices; original algorithm: consecutive junction pairs must be connected.
pub fn has_lr_witness_two_splices(
    junctions: &[Junction],
    transfrags: &[GraphTransfrag],
    graph: &Graph,
    tolerance: u64,
) -> bool {
    if junctions.len() < 2 {
        return true;
    }
    let tol = tolerance as i64;
    let j_match = |j1: Junction, j2: Junction| {
        (j1.donor as i64 - j2.donor as i64).abs() <= tol
            && (j1.acceptor as i64 - j2.acceptor as i64).abs() <= tol
    };
    for i in 0..junctions.len() - 1 {
        let j_curr = junctions[i];
        let j_next = junctions[i + 1];
        let mut pair_connected = false;
        for tf in transfrags {
            if !tf.longread {
                continue;
            }
            let tf_js = path_junctions(&tf.node_ids, graph);
            let has_curr = tf_js.iter().any(|&j| j_match(j, j_curr));
            let has_next = tf_js.iter().any(|&j| j_match(j, j_next));
            if has_curr && has_next {
                pair_connected = true;
                break;
            }
        }
        if !pair_connected {
            return false;
        }
    }
    true
}

const ERROR_PERC: f64 = 0.1;
const CHI_WIN: f64 = 100.0;

#[inline]
pub fn cov_sign_window(
    bpcov: &BpcovStranded,
    strand_idx: usize,
    start_idx: usize,
    end_idx: usize,
    longreads: bool,
) -> f64 {
    let base = bpcov.get_cov_range(BPCOV_STRAND_ALL, start_idx, end_idx);
    if !longreads {
        return base;
    }
    let opp = if strand_idx == BPCOV_STRAND_PLUS {
        BPCOV_STRAND_MINUS
    } else {
        BPCOV_STRAND_PLUS
    };
    base - bpcov.get_cov_range(opp, start_idx, end_idx)
}

#[inline]
fn cj_ok_to_demote(cur: &CJunction, better: &CJunction) -> bool {
    let nb = better.nreads.max(0.0);
    if nb <= 0.0 {
        return false;
    }
    let na = cur.nreads.max(0.0);
    let thr = if na <= 25.0 {
        1.0
    } else if na <= 100.0 {
        0.50
    } else if na <= 200.0 {
        0.25
    } else {
        0.10
    };
    (na / nb) <= thr
}

#[inline]
fn cj_higherr_candidate(cj: &CJunction) -> bool {
    cj.strand != 0 && !cj.guide_match && cj.nm > 0.0 && cj.nm + 1e-9 >= cj.nreads
}

#[inline]
fn cj_reliable(cj: &CJunction) -> bool {
    cj.guide_match || cj.nm + 1e-9 < cj.nreads
}

fn resolved_cj_index(_order: &[usize], slot_neg: f64) -> Option<usize> {
    // Redirect pointers now store the cjunctions ARRAY index directly
    // (as -(idx+1)), so we just decode without going through the order array.
    let idx = (-slot_neg).round() as isize - 1;
    if idx < 0 {
        return None;
    }
    Some(idx as usize)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bpcov::Bpcov;

    #[test]
    fn cov_sign_window_longreads_subtracts_opposite_strand() {
        let plus = Bpcov::from_cov(vec![5.0, 5.0, 5.0, 0.0], 0, 4);
        let minus = Bpcov::from_cov(vec![2.0, 2.0, 2.0, 0.0], 0, 4);
        let bpcov = BpcovStranded {
            minus,
            plus,
            all: Some(vec![7.0, 7.0, 7.0, 0.0]),
        };

        // Long-read: signed window = all - opposite strand.
        let plus_signed = cov_sign_window(&bpcov, BPCOV_STRAND_PLUS, 0, 3, true);
        let minus_signed = cov_sign_window(&bpcov, BPCOV_STRAND_MINUS, 0, 3, true);
        assert!((plus_signed - 15.0).abs() < 1e-9);
        assert!((minus_signed - 6.0).abs() < 1e-9);

        // Non-longread path: no subtraction, use all coverage directly.
        let non_long = cov_sign_window(&bpcov, BPCOV_STRAND_PLUS, 0, 3, false);
        assert!((non_long - 21.0).abs() < 1e-9);
    }

    #[test]
    fn higherr_demotes_acceptor_variant() {
        let mut cjunctions = vec![
            CJunction {
                start: 23044188,
                end: 23044337,
                strand: 1,
                nreads: 15.0,
                nreads_good: 15.0,
                leftsupport: 18.0,
                rightsupport: 15.0,
                nm: 15.0,
                mm: 15.0,
                guide_match: false,
                consleft: -1,
                consright: -1,
                rcount: 0,
            },
            CJunction {
                start: 23044188,
                end: 23044355,
                strand: 1,
                nreads: 2.0,
                nreads_good: 2.0,
                leftsupport: 18.0,
                rightsupport: 2.0,
                nm: 2.0,
                mm: 2.0,
                guide_match: false,
                consleft: -1,
                consright: -1,
                rcount: 0,
            },
            CJunction {
                start: 23044188,
                end: 23044839,
                strand: 1,
                nreads: 1.0,
                nreads_good: 1.0,
                leftsupport: 18.0,
                rightsupport: 21.0,
                nm: 1.0,
                mm: 1.0,
                guide_match: false,
                consleft: -1,
                consright: -1,
                rcount: 0,
            },
        ];

        // Use a higher junction threshold so this low-support junction is considered demotable.
        apply_higherr_demotions(&mut cjunctions, 25, 5.0, None, 0);

        assert!(cjunctions[1].nreads_good < 0.0);
        assert_eq!(cjunctions[1].nreads, 2.0);
    }

    #[test]
    fn higherr_demotes_donor_variant() {
        let mut cjunctions = vec![
            CJunction {
                start: 23050320,
                end: 23050466,
                strand: 1,
                nreads: 87.0,
                nreads_good: 87.0,
                leftsupport: 87.0,
                rightsupport: 88.0,
                nm: 87.0,
                mm: 87.0,
                guide_match: false,
                consleft: -1,
                consright: -1,
                rcount: 0,
            },
            CJunction {
                start: 23050324,
                end: 23050466,
                strand: 1,
                nreads: 1.0,
                nreads_good: 1.0,
                leftsupport: 1.0,
                rightsupport: 88.0,
                nm: 1.0,
                mm: 1.0,
                guide_match: false,
                consleft: -1,
                consright: -1,
                rcount: 0,
            },
        ];

        apply_higherr_demotions(&mut cjunctions, 25, 1.0, None, 0);

        assert!(cjunctions[1].nreads < 0.0);
        assert_eq!(cjunctions[1].nreads_good, 1.0);
    }
}

/// Aggregate leftsupport by donor site and rightsupport by acceptor site (per strand), map-based variant.
/// Matches  before good_junc gating, accumulates support
/// from all junctions sharing the same splice site so the witness checks (step 4)
/// and low_splice_frac (step 6) use the aggregate.
///
/// quirk: the last donor group (left) and last acceptor group (right) are NOT
/// propagated — the loop ends without a trailing propagation pass.
pub fn aggregate_splice_site_support_stats(stats: &mut JunctionStats) {
    if stats.len() < 2 {
        return;
    }

    let gjd = goodjunc_trace_active();

    // Part A: leftsupport aggregation by donor
    // junction array is sorted by start (donor). Group by donor, accumulate
    // leftsupport per strand, propagate to all group members when entering next group.
    {
        let mut keys: Vec<Junction> = stats.keys().copied().collect();
        keys.sort_by_key(|j| (j.donor, j.acceptor));

        let mut current_donor: u64 = u64::MAX;
        let mut group_start_idx: usize = 0;
        let mut leftsupport: [f64; 2] = [0.0, 0.0];

        for idx in 0..keys.len() {
            let jn = keys[idx];
            let strand = stats.get(&jn).map(|s| s.strand.unwrap_or(0)).unwrap_or(0);

            if jn.donor != current_donor {
                // Propagate to previous group
                for j in group_start_idx..idx {
                    let prev_jn = keys[j];
                    if let Some(prev_st) = stats.get_mut(&prev_jn) {
                        let ps = prev_st.strand.unwrap_or(0);
                        if ps != 0 {
                            let si = ((1 + ps as i32) / 2) as usize;
                            if gjd && (prev_st.leftsupport - leftsupport[si]).abs() > 0.01 {
                                eprintln!(
                                    "AGG_LEFT {}-{} strand={} old={:.1} new={:.1}",
                                    prev_jn.donor,
                                    prev_jn.acceptor,
                                    ps,
                                    prev_st.leftsupport,
                                    leftsupport[si]
                                );
                            }
                            prev_st.leftsupport = leftsupport[si];
                        }
                    }
                }

                // Reset and init new group
                leftsupport = [0.0, 0.0];
                if strand != 0 {
                    let si = ((1 + strand as i32) / 2) as usize;
                    leftsupport[si] = stats.get(&jn).map(|s| s.leftsupport).unwrap_or(0.0);
                }
                current_donor = jn.donor;
                group_start_idx = idx;
            } else if strand != 0 {
                // Same donor group — accumulate
                let si = ((1 + strand as i32) / 2) as usize;
                leftsupport[si] += stats.get(&jn).map(|s| s.leftsupport).unwrap_or(0.0);
            }
        }
        // quirk: last donor group is NOT propagated (loop ends at 14533)
    }

    // Part B: rightsupport aggregation by acceptor
    // ejunction is sorted by end (acceptor). Same algorithm for rightsupport.
    {
        let mut keys: Vec<Junction> = stats.keys().copied().collect();
        keys.sort_by_key(|j| (j.acceptor, j.donor));

        let mut current_acceptor: u64 = u64::MAX;
        let mut group_start_idx: usize = 0;
        let mut rightsupport: [f64; 2] = [0.0, 0.0];

        for idx in 0..keys.len() {
            let jn = keys[idx];
            let strand = stats.get(&jn).map(|s| s.strand.unwrap_or(0)).unwrap_or(0);

            if jn.acceptor != current_acceptor {
                // Propagate to previous group
                for j in group_start_idx..idx {
                    let prev_jn = keys[j];
                    if let Some(prev_st) = stats.get_mut(&prev_jn) {
                        let ps = prev_st.strand.unwrap_or(0);
                        if ps != 0 {
                            let si = ((1 + ps as i32) / 2) as usize;
                            if gjd && (prev_st.rightsupport - rightsupport[si]).abs() > 0.01 {
                                eprintln!(
                                    "AGG_RIGHT {}-{} strand={} old={:.1} new={:.1}",
                                    prev_jn.donor,
                                    prev_jn.acceptor,
                                    ps,
                                    prev_st.rightsupport,
                                    rightsupport[si]
                                );
                            }
                            prev_st.rightsupport = rightsupport[si];
                        }
                    }
                }

                // Reset and init new group
                rightsupport = [0.0, 0.0];
                if strand != 0 {
                    let si = ((1 + strand as i32) / 2) as usize;
                    rightsupport[si] = stats.get(&jn).map(|s| s.rightsupport).unwrap_or(0.0);
                }
                current_acceptor = jn.acceptor;
                group_start_idx = idx;
            } else if strand != 0 {
                // Same acceptor group — accumulate
                let si = ((1 + strand as i32) / 2) as usize;
                rightsupport[si] += stats.get(&jn).map(|s| s.rightsupport).unwrap_or(0.0);
            }
        }
        // quirk: last acceptor group is NOT propagated
    }

    // Part C: Same-position opposite-strand conflict resolution
    // Not applicable in Rust: Junction key is (donor, acceptor) without strand,
    // so only one entry per coordinate pair exists.

    // Part D: Consensus splice-site kills ( 14493)
    // Skipped: genome sequence not available in -L mode without -G genome.fa.
    // The if(bdata->gseq) guard means also skips this.
}

/// Strict good_junc/good_merge_junc-style gating for junctions, map-based variant.
/// Marks rejected junctions with strand=0; downstream killed-junction handling uses this bit.
pub fn good_junc_stats(
    stats: &mut JunctionStats,
    bpcov: &BpcovStranded,
    refstart: u64,
    junction_thr: f64,
    long_intron: u64,
    isofrac_percent: f64,
    longreads: bool,
    eonly: bool,
) {
    if stats.is_empty() {
        return;
    }

    let mut keys: Vec<Junction> = stats.keys().copied().collect();
    keys.sort_by_key(|j| (j.donor, j.acceptor));

    let gjd = goodjunc_trace_active();
    for jn in &keys {
        let Some(st) = stats.get_mut(jn) else {
            continue;
        };
        let strand = st.strand.unwrap_or(0);
        if strand == 0 {
            continue;
        }
        if gjd {
            eprintln!(
                "consider junction:{}-{}:{} support={:.6},{:.6} nreads={:.6}",
                jn.donor,
                trace_junction_acceptor(jn),
                strand,
                st.leftsupport,
                st.rightsupport,
                st.mrcount
            );
        }

        // eonly kills all non-guide junctions.
        if eonly && !st.guide_match {
            if gjd {
                eprintln!(
                    "GJ_KILL {}-{} reason=eonly mrcount={:.1} nreads_good={:.1}",
                    jn.donor,
                    trace_junction_acceptor(jn),
                    st.mrcount,
                    st.nreads_good
                );
            }
            st.strand = Some(0);
            continue;
        }

        // guide-matched junctions always pass.
        if st.guide_match {
            continue;
        }

        // good_junc step 3: minimum support.
        if st.nreads_good < junction_thr {
            if gjd {
                eprintln!(
                    "GJ_KILL {}-{} reason=min_support nreads_good={:.1} thr={:.1}",
                    jn.donor,
                    trace_junction_acceptor(jn),
                    st.nreads_good,
                    junction_thr
                );
            }
            st.strand = Some(0);
            continue;
        }

        // higherr (+ 15365):
        // all-bad junctions (nm >= nreads) under 1.25*junctionthr are marked with mm=-1
        // (mm=-1 there, then strand=0 in read loop only if good_junc fails).
        // This gate is critical for long-read weak all-bad junction splitting behavior.
        let all_bad = st.nm > 0.0 && st.nm + 1e-9 >= st.mrcount;
        if all_bad && st.nreads_good < 1.25 * junction_thr {
            if gjd {
                eprintln!(
                    "GJ_MARK_BAD_MM {}-{} reason=higherr_low_support_bad nm={:.1} mrcount={:.1} good={:.1} thr={:.2}",
                    jn.donor,
                    trace_junction_acceptor(jn),
                    st.nm,
                    st.mrcount,
                    st.nreads_good,
                    junction_thr
                );
            }
            st.mm = -1.0;  // Mark for BAD_MM_NEG stage, don't kill yet
            // Don't continue - let junction survive to apply_bad_mm_neg_stage
        }

        // good_junc step 5: all-bad long intron with low count ( jd.nm >= jd.nreads).
        let mismatch = all_bad;
        if mismatch
            && jn.acceptor.saturating_sub(jn.donor) > long_intron
            && st.mrcount < CHI_WIN * ERROR_PERC
        {
            if gjd {
                eprintln!(
                    "GJ_KILL {}-{} reason=bad_long_intron nm={:.1} mrcount={:.1}",
                    jn.donor,
                    trace_junction_acceptor(jn),
                    st.nm,
                    st.mrcount
                );
            }
            st.strand = Some(0);
            continue;
        }

        // good_junc step 4: local coverage witness at donor/acceptor.
        let mut mult = 1.0 / ERROR_PERC;
        if longreads {
            mult /= ERROR_PERC;
        }
        let bw = 5u64;
        let sno = if strand > 0 {
            BPCOV_STRAND_PLUS
        } else {
            BPCOV_STRAND_MINUS
        };

        let j_left = jn.donor.saturating_sub(refstart) as i64;
        let mut lleftcov = 0.0;
        let mut lrightcov = 0.0;
        if j_left >= bw as i64 && (j_left + bw as i64 + 1) < bpcov.plus.cov.len() as i64 {
            let a = (j_left - bw as i64 + 1) as usize;
            let b = (j_left + 1) as usize;
            let c = (j_left + 1) as usize;
            let d = (j_left + bw as i64 + 1) as usize;
            lleftcov = cov_sign_window(bpcov, sno, a, b, longreads);
            lrightcov = cov_sign_window(bpcov, sno, c, d, longreads);
        }
        if lleftcov > 1.0 / ERROR_PERC
            && st.leftsupport * mult < ERROR_PERC * lleftcov
            && (mismatch || lrightcov > lleftcov * (1.0 - ERROR_PERC))
        {
            if gjd {
                eprintln!("GJ_KILL {}-{} reason=left_witness lsup={:.1} lleftcov={:.1} lrightcov={:.1} mult={:.1} mm={}", jn.donor, trace_junction_acceptor(jn), st.leftsupport, lleftcov, lrightcov, mult, mismatch);
            }
            st.strand = Some(0);
            continue;
        }

        let j_right = jn.acceptor.saturating_sub(refstart).saturating_sub(1) as i64;
        let mut rleftcov = 0.0;
        let mut rrightcov = 0.0;
        if j_right >= bw as i64 && (j_right + bw as i64 + 1) < bpcov.plus.cov.len() as i64 {
            let a = (j_right - bw as i64 + 1) as usize;
            let b = (j_right + 1) as usize;
            let c = (j_right + 1) as usize;
            let d = (j_right + bw as i64 + 1) as usize;
            rleftcov = cov_sign_window(bpcov, sno, a, b, longreads);
            rrightcov = cov_sign_window(bpcov, sno, c, d, longreads);
        }
        if rrightcov > 1.0 / ERROR_PERC
            && st.rightsupport * mult < ERROR_PERC * rrightcov
            && (mismatch || rleftcov > rrightcov * (1.0 - ERROR_PERC))
        {
            if gjd {
                eprintln!("GJ_KILL {}-{} reason=right_witness rsup={:.1} rleftcov={:.1} rrightcov={:.1} mult={:.1} mm={}", jn.donor, trace_junction_acceptor(jn), st.rightsupport, rleftcov, rrightcov, mult, mismatch);
            }
            st.strand = Some(0);
            continue;
        }

        // good_junc step 6: extremely low splice fraction relative to side support.
        // uses aggregated values here (runs before good_junc).
        if st.mrcount * 10.0 < ERROR_PERC * st.leftsupport
            || st.mrcount * 10.0 < ERROR_PERC * st.rightsupport
        {
            if gjd {
                eprintln!(
                    "GJ_KILL {}-{} reason=low_splice_frac mrcount={:.1} lsup={:.1} rsup={:.1}",
                    jn.donor,
                    trace_junction_acceptor(jn),
                    st.mrcount,
                    st.leftsupport,
                    st.rightsupport
                );
            }
            st.strand = Some(0);
            continue;
        }

        if gjd {
            eprintln!(
                "--- good_junc: ACCEPTED {}-{}:{} nreads={:.1} left={:.1} right={:.1}",
                jn.donor,
                trace_junction_acceptor(jn),
                strand,
                st.mrcount,
                st.leftsupport,
                st.rightsupport
            );
        }
    }

    // good_merge_junc donor-level isofrac guard.
    if isofrac_percent > 0.0 {
        let mut donor_sum: HashMap<u64, f64> = Default::default();
        for (j, s) in stats.iter() {
            if s.strand != Some(0) {
                *donor_sum.entry(j.donor).or_insert(0.0) += s.nreads_good;
            }
        }
        for jn in keys {
            let Some(st) = stats.get_mut(&jn) else {
                continue;
            };
            if st.strand == Some(0) {
                continue;
            }
            let sum = donor_sum.get(&jn.donor).copied().unwrap_or(0.0);
            if sum > 0.0 && st.nreads_good * 100.0 / sum < isofrac_percent {
                st.strand = Some(0);
            }
        }
    }
}

/// Junctions that will be treated as bad during the per-read pass.
/// Besides explicit strand=0 kills, the original algorithm also splits reads on junctions
/// demoted on either side (`nreads<0` / `nreads_good<0`) and on `mm<0`.
pub fn compute_killed_junction_pairs(cjunctions: &[CJunction]) -> HashSet<Junction> {
    cjunctions
        .iter()
        .filter(|cj| cj.strand == 0 || cj.mm < 0.0 || cj.nreads < 0.0 || cj.nreads_good < 0.0)
        .map(|cj| Junction::new(cj.start, cj.end))
        .collect()
}

#[inline]
fn cjunc_strand_index(strand: i8) -> Option<usize> {
    match strand {
        -1 => Some(0),
        1 => Some(1),
        _ => None,
    }
}

/// Aggregate splice-site support using CJunction transport (pre-good_junc pass).
pub fn aggregate_splice_site_support(cjunctions: &mut Vec<CJunction>) {
    if cjunctions.len() < 2 {
        return;
    }
    let gjd = goodjunc_trace_active();

    // uses the start-sorted junction vector, with a second end-sorted copy.
    let mut donor_order: Vec<usize> = (0..cjunctions.len()).collect();
    donor_order.sort_by_key(|&i| (cjunctions[i].start, cjunctions[i].end, cjunctions[i].strand));

    let mut current_start = u64::MAX;
    let mut group_start = 0usize;
    let mut leftsupport = [0.0f64; 2];
    for ord in 0..donor_order.len() {
        let idx = donor_order[ord];
        let start = cjunctions[idx].start;
        let strand = cjunctions[idx].strand;
        if start != current_start {
            for prev_ord in group_start..ord {
                let prev_idx = donor_order[prev_ord];
                let prev = &mut cjunctions[prev_idx];
                if let Some(si) = cjunc_strand_index(prev.strand) {
                    if gjd && (prev.leftsupport - leftsupport[si]).abs() > 0.01 {
                        eprintln!(
                            "AGG_LEFT {}-{} strand={} old={:.1} new={:.1}",
                            prev.start, prev.end, prev.strand, prev.leftsupport, leftsupport[si]
                        );
                    }
                    prev.leftsupport = leftsupport[si];
                }
            }

            // Same-position opposite-strand conflict resolution.
            if strand != 0 {
                let mut j = ord + 1;
                while j < donor_order.len() {
                    let other_idx = donor_order[j];
                    let other = cjunctions[other_idx];
                    if other.start != start || other.end != cjunctions[idx].end {
                        break;
                    }
                    if other.strand != 0 && other.strand != strand {
                        if cjunctions[idx].nreads > other.nreads
                            && cjunctions[idx].nreads_good > other.nreads_good
                        {
                            cjunctions[other_idx].strand = 0;
                        } else if cjunctions[idx].nreads < other.nreads
                            && cjunctions[idx].nreads_good < other.nreads_good
                        {
                            cjunctions[idx].strand = 0;
                        }
                    }
                    j += 1;
                }
            }

            leftsupport = [0.0, 0.0];
            if let Some(si) = cjunc_strand_index(cjunctions[idx].strand) {
                leftsupport[si] = cjunctions[idx].leftsupport;
            }
            current_start = start;
            group_start = ord;
        } else if let Some(si) = cjunc_strand_index(strand) {
            leftsupport[si] += cjunctions[idx].leftsupport;
        }
    }
    // quirk: final donor group is not propagated.

    let mut acceptor_order: Vec<usize> = (0..cjunctions.len()).collect();
    acceptor_order.sort_by_key(|&i| (cjunctions[i].end, cjunctions[i].start, cjunctions[i].strand));

    let mut current_end = u64::MAX;
    let mut group_start = 0usize;
    let mut rightsupport = [0.0f64; 2];
    for ord in 0..acceptor_order.len() {
        let idx = acceptor_order[ord];
        let end = cjunctions[idx].end;
        let strand = cjunctions[idx].strand;
        if end != current_end {
            for prev_ord in group_start..ord {
                let prev_idx = acceptor_order[prev_ord];
                let prev = &mut cjunctions[prev_idx];
                if let Some(si) = cjunc_strand_index(prev.strand) {
                    if gjd && (prev.rightsupport - rightsupport[si]).abs() > 0.01 {
                        eprintln!(
                            "AGG_RIGHT {}-{} strand={} old={:.1} new={:.1}",
                            prev.start, prev.end, prev.strand, prev.rightsupport, rightsupport[si]
                        );
                    }
                    prev.rightsupport = rightsupport[si];
                }
            }

            rightsupport = [0.0, 0.0];
            if let Some(si) = cjunc_strand_index(strand) {
                rightsupport[si] = cjunctions[idx].rightsupport;
            }
            current_end = end;
            group_start = ord;
        } else if let Some(si) = cjunc_strand_index(strand) {
            rightsupport[si] += cjunctions[idx].rightsupport;
        }
    }
    // quirk: final acceptor group is not propagated.
}

/// good_junc/good_merge_junc gate using CJunction transport.
/// Demote "run-through" junctions whose donor coincides with a strong read-END
/// (LEND) peak AND acceptor coincides with a strong read-START (LSTART) peak.
/// These are typically chimeric reads bridging the TTS of one gene and TSS of
/// the next gene — NOT real splice sites. StringTie-faithful bundles (per
/// GGO_19.log analysis) split STRG.29/STRG.31 across such boundaries; we
/// replicate by killing the bridging junction before graph construction.
///
/// Activated by `RUSTLE_RUNTHROUGH_DEMOTE=<thr>` where `thr` is the minimum
/// LEND and LSTART count each (default unset = disabled).
pub fn demote_runthrough_junctions(
    cjunctions: &mut Vec<crate::types::CJunction>,
    bpcov: &BpcovStranded,
    refstart: u64,
) {
    // NOTE: Infrastructure kept but DISABLED by default. Coverage-based discrimination
    // cannot distinguish run-through junctions from real introns (both have low intron cov
    // and high flanking cov). On GGO_19 the criterion catastrophically regresses matches.
    // The real STRG.29/STRG.31 divergence lives in path extraction / max_flow, not here.
    // Gate: RUSTLE_RUNTHROUGH_DEMOTE=<ratio> — kept for future experiments.
    let ratio_thr: f64 = match std::env::var("RUSTLE_RUNTHROUGH_DEMOTE")
        .ok()
        .and_then(|v| v.parse().ok())
    {
        Some(v) if v > 0.0 => v,
        _ => return,
    };
    let min_flank_cov: f64 = std::env::var("RUSTLE_RUNTHROUGH_MIN_FLANK")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(10.0);
    let min_intron_bp: u64 = std::env::var("RUSTLE_RUNTHROUGH_MIN_INTRON_BP")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(100);
    if cjunctions.is_empty() {
        return;
    }
    let gjd = goodjunc_trace_active();
    let cov_len = bpcov.plus.cov.len();
    for cj in cjunctions.iter_mut() {
        if cj.strand == 0 {
            continue;
        }
        // Intron range in Rustle convention: cj.start = last exon base, cj.end = last intron base
        // Actually from our trace: cj.start=17451950 cj.end=17452250, corresponding to samtools junction
        // 17451951-17452250 (intron: first base = 17451951, last base = 17452250).
        // So cj.start = last exon base; intron = (cj.start+1, cj.end); acceptor = cj.end + 1.
        let intron_start = cj.start.saturating_add(1);
        let intron_end = cj.end;
        if intron_end <= intron_start || intron_end - intron_start < min_intron_bp {
            continue;
        }
        if cj.start < refstart || intron_end < refstart {
            continue;
        }
        let donor_idx = (cj.start - refstart) as usize;
        let acceptor_idx = (cj.end.saturating_add(1) - refstart) as usize;
        let intron_len = (intron_end - intron_start + 1) as usize;
        let mid_idx = ((intron_start as i64 + intron_end as i64) / 2 - refstart as i64) as usize;
        if mid_idx >= cov_len || donor_idx >= cov_len || acceptor_idx >= cov_len {
            continue;
        }
        // Flanking cov: 50bp windows on exon sides.
        let win: usize = 50;
        let donor_lo = donor_idx.saturating_sub(win.saturating_sub(1));
        let donor_hi = donor_idx;
        let acceptor_lo = acceptor_idx;
        let acceptor_hi = (acceptor_idx + win - 1).min(cov_len - 1);
        let donor_cov = bpcov.get_cov_range(BPCOV_STRAND_ALL, donor_lo, donor_hi + 1)
            / ((donor_hi - donor_lo + 1) as f64);
        let acceptor_cov = bpcov.get_cov_range(BPCOV_STRAND_ALL, acceptor_lo, acceptor_hi + 1)
            / ((acceptor_hi - acceptor_lo + 1) as f64);
        // Intron middle cov: average over middle 50% of intron.
        let inner_start = mid_idx.saturating_sub(intron_len / 4);
        let inner_end = (mid_idx + intron_len / 4).min(cov_len - 1);
        if inner_end <= inner_start {
            continue;
        }
        let intron_cov = bpcov.get_cov_range(BPCOV_STRAND_ALL, inner_start, inner_end + 1)
            / ((inner_end - inner_start + 1) as f64);
        if donor_cov < min_flank_cov || acceptor_cov < min_flank_cov {
            continue;
        }
        let ratio_d = donor_cov / intron_cov.max(0.5);
        let ratio_a = acceptor_cov / intron_cov.max(0.5);
        if ratio_d >= ratio_thr && ratio_a >= ratio_thr {
            if gjd {
                eprintln!(
                    "RUNTHROUGH_DEMOTE {}-{} donor_cov={:.1} intron_cov={:.2} acceptor_cov={:.1} ratio_d={:.1} ratio_a={:.1}",
                    cj.start, cj.end, donor_cov, intron_cov, acceptor_cov, ratio_d, ratio_a
                );
            }
            cj.strand = 0;
        }
    }
}

pub fn good_junc(
    cjunctions: &mut Vec<CJunction>,
    bpcov: &BpcovStranded,
    refstart: u64,
    junction_thr: f64,
    long_intron: u64,
    isofrac_percent: f64,
    longreads: bool,
    eonly: bool,
) {
    if cjunctions.is_empty() {
        return;
    }
    let mut order: Vec<usize> = (0..cjunctions.len()).collect();
    order.sort_by_key(|&i| (cjunctions[i].start, cjunctions[i].end, cjunctions[i].strand));

    let gjd = goodjunc_trace_active();
    for &idx in &order {
        let cj = &mut cjunctions[idx];
        let strand = cj.strand;
        if strand == 0 {
            continue;
        }
        if gjd {
            eprintln!(
                "consider junction:{}-{}:{} support={:.6},{:.6} nreads={:.6}",
                cj.start,
                trace_cjunction_acceptor(cj),
                strand,
                cj.leftsupport,
                cj.rightsupport,
                cj.nreads
            );
        }

        if eonly && !cj.guide_match {
            if gjd {
                eprintln!(
                    "GJ_KILL {}-{} reason=eonly mrcount={:.1} nreads_good={:.1}",
                    cj.start,
                    trace_cjunction_acceptor(cj),
                    cj.nreads,
                    cj.nreads_good
                );
            }
            cj.strand = 0;
            continue;
        }

        if cj.guide_match {
            continue;
        }

        // `higherr` stores negative "demotion pointers" in `nreads` / `nreads_good`.
        // Those are not counts and must not participate in support-based killing here.
        // Such junctions are later treated as "bad" during read splitting and graph building.
        if cj.nreads < 0.0 || cj.nreads_good < 0.0 {
            continue;
        }

        if cj.nreads_good < junction_thr {
            if gjd {
                eprintln!(
                    "GJ_KILL {}-{} reason=min_support nreads_good={:.1} thr={:.1}",
                    cj.start,
                    trace_cjunction_acceptor(cj),
                    cj.nreads_good,
                    junction_thr
                );
            }
            cj.strand = 0;
            continue;
        }

        let all_bad = cj.nm > 0.0 && cj.nm + 1e-9 >= cj.nreads;
        if all_bad && cj.nreads_good < 1.25 * junction_thr {
            // Optional rescue: if bpcov is DIS-continuous at the donor (right << left),
            // this looks like a real splice site despite low mismatch-free support — skip
            // the mm=-1 demotion. Mirrors StringTie rlink.cpp:14958 good_junc() rescue.
            let mut rescued = false;
            if std::env::var_os("RUSTLE_RESCUE_LOWSUPPORT_BAD").is_some() && cj.start > refstart {
                let donor_idx = (cj.start - refstart) as usize;
                let len = bpcov.plus.cov.len();
                if donor_idx > 0 && donor_idx + 1 < len {
                    let leftcov = bpcov.get_cov_range(
                        crate::bpcov::BPCOV_STRAND_ALL,
                        donor_idx.saturating_sub(1),
                        donor_idx,
                    );
                    let rightcov = bpcov.get_cov_range(
                        crate::bpcov::BPCOV_STRAND_ALL,
                        donor_idx,
                        donor_idx + 1,
                    );
                    if leftcov > 0.0 && rightcov < 0.9 * leftcov {
                        rescued = true;
                    }
                }
            }
            if !rescued {
                cj.mm = -1.0;
                if gjd {
                    eprintln!(
                        "GJ_DEFER {}-{} reason=higherr_low_support_bad mm=-1.0",
                        cj.start,
                        trace_cjunction_acceptor(cj)
                    );
                }
            } else if gjd {
                eprintln!(
                    "GJ_RESCUE {}-{} reason=discontinuous_bpcov",
                    cj.start,
                    trace_cjunction_acceptor(cj)
                );
            }
            // Don't continue; let the junction survive to be processed by BAD_MM_NEG logic
        }

        let mismatch = all_bad;
        if mismatch
            && cj.end.saturating_sub(cj.start) > long_intron
            && cj.nreads < CHI_WIN * ERROR_PERC
        {
            if gjd {
                eprintln!(
                    "GJ_KILL {}-{} reason=bad_long_intron nm={:.1} mrcount={:.1}",
                    cj.start,
                    trace_cjunction_acceptor(cj),
                    cj.nm,
                    cj.nreads
                );
            }
            cj.strand = 0;
            continue;
        }

        let mut mult = 1.0 / ERROR_PERC;
        if longreads {
            mult /= ERROR_PERC;
        }
        let bw = 5u64;
        let sno = if strand > 0 {
            BPCOV_STRAND_PLUS
        } else {
            BPCOV_STRAND_MINUS
        };

        let j_left = cj.start.saturating_sub(refstart) as i64;
        let mut lleftcov = 0.0;
        let mut lrightcov = 0.0;
        if j_left >= bw as i64 && (j_left + bw as i64 + 1) < bpcov.plus.cov.len() as i64 {
            let a = (j_left - bw as i64 + 1) as usize;
            let b = (j_left + 1) as usize;
            let c = (j_left + 1) as usize;
            let d = (j_left + bw as i64 + 1) as usize;
            lleftcov = cov_sign_window(bpcov, sno, a, b, longreads);
            lrightcov = cov_sign_window(bpcov, sno, c, d, longreads);
        }
        if lleftcov > 1.0 / ERROR_PERC
            && cj.leftsupport * mult < ERROR_PERC * lleftcov
            && (mismatch || lrightcov > lleftcov * (1.0 - ERROR_PERC))
        {
            if gjd {
                eprintln!(
                    "GJ_KILL {}-{} reason=left_witness lsup={:.1} lleftcov={:.1} lrightcov={:.1} mult={:.1} mm={}",
                    cj.start,
                    trace_cjunction_acceptor(cj),
                    cj.leftsupport,
                    lleftcov,
                    lrightcov,
                    mult,
                    mismatch
                );
            }
            cj.strand = 0;
            continue;
        }

        let j_right = cj.end.saturating_sub(refstart).saturating_sub(1) as i64;
        let mut rleftcov = 0.0;
        let mut rrightcov = 0.0;
        if j_right >= bw as i64 && (j_right + bw as i64 + 1) < bpcov.plus.cov.len() as i64 {
            let a = (j_right - bw as i64 + 1) as usize;
            let b = (j_right + 1) as usize;
            let c = (j_right + 1) as usize;
            let d = (j_right + bw as i64 + 1) as usize;
            rleftcov = cov_sign_window(bpcov, sno, a, b, longreads);
            rrightcov = cov_sign_window(bpcov, sno, c, d, longreads);
        }
        if rrightcov > 1.0 / ERROR_PERC
            && cj.rightsupport * mult < ERROR_PERC * rrightcov
            && (mismatch || rleftcov > rrightcov * (1.0 - ERROR_PERC))
        {
            if gjd {
                eprintln!(
                    "GJ_KILL {}-{} reason=right_witness rsup={:.1} rleftcov={:.1} rrightcov={:.1} mult={:.1} mm={}",
                    cj.start,
                    trace_cjunction_acceptor(cj),
                    cj.rightsupport,
                    rleftcov,
                    rrightcov,
                    mult,
                    mismatch
                );
            }
            cj.strand = 0;
            continue;
        }

        if cj.nreads * 10.0 < ERROR_PERC * cj.leftsupport
            || cj.nreads * 10.0 < ERROR_PERC * cj.rightsupport
        {
            // if all_bad, don't kill here; mark for BAD_MM_NEG instead.
            if all_bad {
                cj.mm = -1.0;
                if gjd {
                    eprintln!(
                        "GJ_DEFER {}-{} reason=low_splice_frac_all_bad mm=-1.0",
                        cj.start,
                        trace_cjunction_acceptor(cj)
                    );
                }
            } else {
                if gjd {
                    eprintln!(
                        "GJ_KILL {}-{} reason=low_splice_frac mrcount={:.1} lsup={:.1} rsup={:.1}",
                        cj.start,
                        trace_cjunction_acceptor(cj),
                        cj.nreads,
                        cj.leftsupport,
                        cj.rightsupport
                    );
                }
                cj.strand = 0;
                continue;
            }
        }

        if gjd {
            eprintln!(
                "--- good_junc: ACCEPTED {}-{}:{} nreads={:.1} left={:.1} right={:.1}",
                cj.start,
                trace_cjunction_acceptor(cj),
                strand,
                cj.nreads,
                cj.leftsupport,
                cj.rightsupport
            );
        }
    }

    if isofrac_percent > 0.0 {
        let mut donor_sum: HashMap<u64, f64> = Default::default();
        for cj in cjunctions.iter() {
            // Exclude demotion pointers and mm<0 markers from the donor sum, matching the
            // behavior of using valid support counts only.
            if cj.strand != 0 && cj.nreads_good >= 0.0 && cj.mm >= 0.0 {
                *donor_sum.entry(cj.start).or_insert(0.0) += cj.nreads_good;
            }
        }
        for &idx in &order {
            let cj = &mut cjunctions[idx];
            if cj.strand == 0 {
                continue;
            }
            let sum = donor_sum.get(&cj.start).copied().unwrap_or(0.0);
            if sum > 0.0 && cj.nreads_good * 100.0 / sum < isofrac_percent {
                cj.strand = 0;
            }
        }
    }
}

/// the original algorithm `higherr`: donor/acceptor demotion pointers on all-bad junctions.
pub fn apply_higherr_demotions(
    cjunctions: &mut [CJunction],
    sserror: u64,
    junction_thr: f64,
    bpcov: Option<&BpcovStranded>,
    refstart: u64,
) {
    if cjunctions.len() < 2 {
        return;
    }

    let gjd = goodjunc_trace_active();
    let tolerance = 1.0 - ERROR_PERC;

    // Run-through junction detection (StringTie rlink.cpp:15015-15026 analog) was
    // tested but disabled. The STRG.120 case (2-read junction outcompeted by 70+
    // run-through reads) has a real coverage drop at the donor position
    // (leftcov=32, rightcov=20, ratio 0.62) — it IS a weak splice site, not a
    // run-through. Discriminating the two via bpcov alone produces false kills.
    //
    // StringTie likely uses per-read cumulative `good_junc` in build_graphs
    // (line 15276) to kill these via incremental state mutation — a mechanism
    // we couldn't cleanly port (see project_color_break_root_cause memory).
    let _ = (bpcov, refstart);

    let mut donor_order: Vec<usize> = (0..cjunctions.len()).collect();
    donor_order.sort_by_key(|&i| (cjunctions[i].start, cjunctions[i].end, cjunctions[i].strand));

    for ord_i in 0..donor_order.len() {
        let idx_i = donor_order[ord_i];
        let cur = cjunctions[idx_i];
        if !cj_higherr_candidate(&cur) {
            continue;
        }
        // StringTie rlink.cpp:15008-15026: for all-bad junctions (nm>=nreads) with
        // nreads_good >= 1.25*junctionthr, StringTie marks mm=-1 if bpcov is
        // continuous across the donor (rightcov > 0.9 * leftcov).
        //
        // Tested but gated behind RUSTLE_HIGHERR_CONT — bugs in bpcov indexing
        // across bundle/subbundle contexts cause false kills (e.g., a 195-read
        // junction was killed because leftcov/rightcov came back as 2,2 from
        // wrong offsets). Also STRG.120's 2-read junction is NOT killed by this
        // filter in Rustle even when the check runs — suggesting our cur.start
        // offset doesn't correspond to the same position StringTie reads.
        // StringTie rlink.cpp:15015-15026: for all-bad junctions with enough support,
        // check bpcov continuity at the donor. If coverage is continuous across the
        // "junction" (rightcov > 0.9 * leftcov), it's a run-through → mark mm=-1.
        // Enabled by default (kills 86 spurious transcripts for -1 match, +3.6% precision).
        // Disable with RUSTLE_HIGHERR_CONT_OFF=1 if it causes regression on your data.
        if cur.nreads_good >= 0.0
            && cur.nreads_good >= 1.25 * junction_thr
            && std::env::var_os("RUSTLE_HIGHERR_CONT_OFF").is_none()
        {
            if let Some(bp) = bpcov {
                if cur.start > refstart {
                    let donor_idx = (cur.start - refstart) as usize;
                    let len = bp.plus.cov.len();
                    if donor_idx > 0 && donor_idx + 1 < len {
                        // StringTie: point=jd.start-refstart, point+1=jd.start+1-refstart
                        // In Rustle, cur.start is 0-based first intron base, so:
                        //   leftcov at last exon base = donor_idx - 1
                        //   rightcov at first intron base = donor_idx
                        let leftcov = bp.get_cov_range(
                            crate::bpcov::BPCOV_STRAND_ALL,
                            donor_idx.saturating_sub(1),
                            donor_idx,
                        );
                        let rightcov = bp.get_cov_range(
                            crate::bpcov::BPCOV_STRAND_ALL,
                            donor_idx,
                            donor_idx + 1,
                        );
                        // Sanity check: if bpcov at the donor is much smaller than the junction's
                        // read support, the bpcov doesn't reflect this junction's context
                        // (e.g., the junction was accumulated from reads in a different subbundle).
                        // Skip to avoid spurious kills of real junctions.
                        if leftcov < cur.nreads / 2.0 {
                            continue;
                        }
                        if leftcov > 0.0 && rightcov > tolerance * leftcov {
                            cjunctions[idx_i].mm = -1.0;
                            if gjd {
                                eprintln!(
                                    "HE_CONT_DEMOTE {}-{} nreads={:.1} leftcov={:.0} rightcov={:.0} ratio={:.2}",
                                    cur.start,
                                    trace_cjunction_acceptor(&cur),
                                    cur.nreads,
                                    leftcov,
                                    rightcov,
                                    rightcov / leftcov
                                );
                            }
                        }
                    }
                }
            }
            // FALL THROUGH to demotion search (StringTie rlink.cpp:15030 onwards
            // runs demotion for ALL higherr candidates, independently of whether
            // the continuity check marked mm=-1).
        }
        let _ = tolerance;

        let mut support = 0.0;
        let mut search = true;
        let mut reliable = false;
        let mut j = ord_i as isize - 1;
        while j >= 0 {
            let idx_j = donor_order[j as usize];
            let cand = cjunctions[idx_j];
            if cur.start.saturating_sub(cand.start) >= sserror {
                break;
            }
            if cand.strand == cur.strand {
                if cand.start == cur.start {
                    if cand.nreads < 0.0 {
                        if let Some(point_idx) = resolved_cj_index(&donor_order, cand.nreads) {
                            if cj_ok_to_demote(&cjunctions[idx_i], &cjunctions[point_idx]) {
                                cjunctions[idx_i].nreads = cand.nreads;
                                search = false;
                            }
                        }
                    }
                    break;
                } else if cj_reliable(&cand) {
                    if cand.leftsupport > cur.leftsupport * tolerance
                        && cj_ok_to_demote(&cjunctions[idx_i], &cand)
                    {
                        reliable = true;
                        // Store the CJUNCTIONS ARRAY index (not sorted order position)
                        // so the per-read redirect can decode it directly.
                        cjunctions[idx_i].nreads = -(idx_j as f64 + 1.0);
                        if gjd {
                            eprintln!(
                                "HE_DEMOTE_LEFT {}-{} -> {}-{} left={:.1}->{:.1}",
                                cur.start,
                                trace_cjunction_acceptor(&cur),
                                cand.start,
                                trace_cjunction_acceptor(&cand),
                                cur.leftsupport,
                                cand.leftsupport
                            );
                        }
                        break;
                    }
                } else if cand.leftsupport > support
                    && cur.start.saturating_sub(cand.start) < sserror
                    && cand.leftsupport * tolerance > cur.leftsupport
                    && cj_ok_to_demote(&cjunctions[idx_i], &cand)
                {
                    cjunctions[idx_i].nreads = -(idx_j as f64 + 1.0);
                    support = cand.leftsupport;
                    if gjd {
                        eprintln!(
                            "HE_DEMOTE_LEFT {}-{} -> {}-{} left={:.1}->{:.1}",
                            cur.start,
                            trace_cjunction_acceptor(&cur),
                            cand.start,
                            trace_cjunction_acceptor(&cand),
                            cur.leftsupport,
                            cand.leftsupport
                        );
                    }
                }
            }
            j -= 1;
        }

        if !search {
            continue;
        }

        let mut dist = sserror as i64;
        if cjunctions[idx_i].nreads < 0.0 {
            if let Some(point_idx) = resolved_cj_index(&donor_order, cjunctions[idx_i].nreads) {
                dist = cjunctions[idx_i].start as i64 - cjunctions[point_idx].start as i64;
            }
        }

        let mut j = ord_i + 1;
        while j < donor_order.len() {
            let idx_j = donor_order[j];
            let cand = cjunctions[idx_j];
            let d = cand.start as i64 - cur.start as i64;
            if d >= sserror as i64 {
                break;
            }
            if cand.strand == cur.strand && cand.start != cur.start {
                if cj_reliable(&cand) {
                    if (d < dist || (d == dist && cand.leftsupport > support))
                        && cand.leftsupport > cur.leftsupport * tolerance
                        && cj_ok_to_demote(&cjunctions[idx_i], &cand)
                    {
                        cjunctions[idx_i].nreads = -(idx_j as f64 + 1.0);
                        if gjd {
                            eprintln!(
                                "HE_DEMOTE_LEFT {}-{} -> {}-{} left={:.1}->{:.1}",
                                cur.start,
                                trace_cjunction_acceptor(&cur),
                                cand.start,
                                trace_cjunction_acceptor(&cand),
                                cur.leftsupport,
                                cand.leftsupport
                            );
                        }
                        break;
                    }
                } else if !reliable
                    && cand.leftsupport > support
                    && d < sserror as i64
                    && cand.leftsupport * tolerance > cur.leftsupport
                    && cj_ok_to_demote(&cjunctions[idx_i], &cand)
                {
                    cjunctions[idx_i].nreads = -(idx_j as f64 + 1.0);
                    support = cand.leftsupport;
                    if gjd {
                        eprintln!(
                            "HE_DEMOTE_LEFT {}-{} -> {}-{} left={:.1}->{:.1}",
                            cur.start,
                            trace_cjunction_acceptor(&cur),
                            cand.start,
                            trace_cjunction_acceptor(&cand),
                            cur.leftsupport,
                            cand.leftsupport
                        );
                    }
                }
            }
            j += 1;
        }
    }

    let mut acceptor_order: Vec<usize> = (0..cjunctions.len()).collect();
    acceptor_order.sort_by_key(|&i| (cjunctions[i].end, cjunctions[i].start, cjunctions[i].strand));

    for ord_i in 0..acceptor_order.len() {
        let idx_i = acceptor_order[ord_i];
        let cur = cjunctions[idx_i];
        if !cj_higherr_candidate(&cur) {
            continue;
        }
        // Port of StringTie rlink.cpp:15112-15120 acceptor-side continuity check.
        // Check bpcov at last intron base vs first downstream exon base: if coverage
        // is continuous across the acceptor (leftcov > 0.9 * rightcov), it's a
        // run-through → mm=-1. Symmetric to donor-side check in start loop above.
        if cur.nreads_good >= 0.0
            && cur.mm >= 0.0
            && cur.nreads_good >= 1.25 * junction_thr
            && std::env::var_os("RUSTLE_HIGHERR_CONT_OFF").is_none()
        {
            if let Some(bp) = bpcov {
                if cur.end > refstart + 1 {
                    // cur.end = 0-based first exon base after intron.
                    // last intron base idx = cur.end - 1 - refstart.
                    let last_intron_idx = (cur.end - 1 - refstart) as usize;
                    let len = bp.plus.cov.len();
                    if last_intron_idx + 2 < len {
                        let leftcov = bp.get_cov_range(
                            crate::bpcov::BPCOV_STRAND_ALL,
                            last_intron_idx,
                            last_intron_idx + 1,
                        );
                        let rightcov = bp.get_cov_range(
                            crate::bpcov::BPCOV_STRAND_ALL,
                            last_intron_idx + 1,
                            last_intron_idx + 2,
                        );
                        // Sanity: bpcov must reflect this junction's context.
                        if rightcov >= cur.nreads / 2.0
                            && rightcov > 0.0
                            && leftcov > tolerance * rightcov
                        {
                            cjunctions[idx_i].mm = -1.0;
                            if gjd {
                                eprintln!(
                                    "HE_CONT_R_DEMOTE {}-{} nreads={:.1} leftcov={:.0} rightcov={:.0} ratio={:.2}",
                                    cur.start,
                                    trace_cjunction_acceptor(&cur),
                                    cur.nreads,
                                    leftcov,
                                    rightcov,
                                    leftcov / rightcov
                                );
                            }
                        }
                    }
                }
            }
            // FALL THROUGH to demotion search — StringTie rlink.cpp:15123 onwards
            // runs neighbor demotion for all higherr candidates regardless of
            // whether the continuity check marked mm=-1.
        }

        let mut support = 0.0;
        let mut search = true;
        let mut reliable = false;
        let mut j = ord_i as isize - 1;
        while j >= 0 {
            let idx_j = acceptor_order[j as usize];
            let cand = cjunctions[idx_j];
            if cur.end.saturating_sub(cand.end) >= sserror {
                break;
            }
            if cand.strand == cur.strand {
                if cand.end == cur.end {
                    if cand.nreads_good < 0.0 {
                        if let Some(point_idx) =
                            resolved_cj_index(&acceptor_order, cand.nreads_good)
                        {
                            if cj_ok_to_demote(&cjunctions[idx_i], &cjunctions[point_idx]) {
                                cjunctions[idx_i].nreads_good = cand.nreads_good;
                                search = false;
                            }
                        }
                    }
                    break;
                } else if cj_reliable(&cand) {
                    if cand.rightsupport > cur.rightsupport * tolerance
                        && cj_ok_to_demote(&cjunctions[idx_i], &cand)
                    {
                        reliable = true;
                        cjunctions[idx_i].nreads_good = -(idx_j as f64 + 1.0);
                        if gjd {
                            eprintln!(
                                "HE_DEMOTE_RIGHT {}-{} -> {}-{} right={:.1}->{:.1}",
                                cur.start,
                                trace_cjunction_acceptor(&cur),
                                cand.start,
                                trace_cjunction_acceptor(&cand),
                                cur.rightsupport,
                                cand.rightsupport
                            );
                        }
                        break;
                    }
                } else if cand.rightsupport > support
                    && cur.end.saturating_sub(cand.end) < sserror
                    && cand.rightsupport * tolerance > cur.rightsupport
                    && cj_ok_to_demote(&cjunctions[idx_i], &cand)
                    && !(std::env::var_os("RUSTLE_ALT_ACCEPTOR_PRESERVE").is_some()
                        && {
                            let d_up = cur.end.saturating_sub(cand.end);
                            d_up >= 4 && d_up < 9 && cur.nreads >= 5.0
                        })
                {
                    cjunctions[idx_i].nreads_good = -(idx_j as f64 + 1.0);
                    support = cand.rightsupport;
                    if gjd {
                        eprintln!(
                            "HE_DEMOTE_RIGHT {}-{} -> {}-{} right={:.1}->{:.1}",
                            cur.start,
                            trace_cjunction_acceptor(&cur),
                            cand.start,
                            trace_cjunction_acceptor(&cand),
                            cur.rightsupport,
                            cand.rightsupport
                        );
                    }
                }
            }
            j -= 1;
        }

        if !search {
            continue;
        }

        let mut dist = sserror as i64;
        if cjunctions[idx_i].nreads_good < 0.0 {
            if let Some(point_idx) =
                resolved_cj_index(&acceptor_order, cjunctions[idx_i].nreads_good)
            {
                dist = cjunctions[point_idx].end as i64 - cjunctions[idx_i].end as i64;
            }
        }

        let mut j = ord_i + 1;
        while j < acceptor_order.len() {
            let idx_j = acceptor_order[j];
            let cand = cjunctions[idx_j];
            let d = cand.end as i64 - cur.end as i64;
            if d >= sserror as i64 {
                break;
            }
            if cand.strand == cur.strand && cand.end != cur.end {
                if cj_reliable(&cand) {
                    if (d < dist || (d == dist && cand.rightsupport > support))
                        && cand.rightsupport > cur.rightsupport * tolerance
                        && cj_ok_to_demote(&cjunctions[idx_i], &cand)
                    {
                        cjunctions[idx_i].nreads_good = -(idx_j as f64 + 1.0);
                        if gjd {
                            eprintln!(
                                "HE_DEMOTE_RIGHT {}-{} -> {}-{} right={:.1}->{:.1}",
                                cur.start,
                                trace_cjunction_acceptor(&cur),
                                cand.start,
                                trace_cjunction_acceptor(&cand),
                                cur.rightsupport,
                                cand.rightsupport
                            );
                        }
                        break;
                    }
                } else if ((!reliable
                    && cand.rightsupport > support
                    && d < sserror as i64
                    && cand.rightsupport * tolerance > cur.rightsupport)
                    || (d < dist && cj_reliable(&cand)))
                    && cj_ok_to_demote(&cjunctions[idx_i], &cand)
                    && !(std::env::var_os("RUSTLE_ALT_ACCEPTOR_PRESERVE").is_some()
                        && d >= 4
                        && d < 9
                        && cur.nreads >= 5.0)
                {
                    cjunctions[idx_i].nreads_good = -((j + 1) as f64);
                    support = cand.rightsupport;
                    if gjd {
                        eprintln!(
                            "HE_DEMOTE_RIGHT {}-{} -> {}-{} right={:.1}->{:.1}",
                            cur.start,
                            trace_cjunction_acceptor(&cur),
                            cand.start,
                            trace_cjunction_acceptor(&cand),
                            cur.rightsupport,
                            cand.rightsupport
                        );
                    }
                }
            }
            j += 1;
        }
    }
}

/// Back-compat alias: CJunction transport wrapper for splice-site support aggregation.
pub fn aggregate_splice_site_support_cjuncs(cjunctions: &mut Vec<CJunction>) {
    aggregate_splice_site_support(cjunctions);
}

/// Back-compat alias: CJunction transport wrapper for good-junction gating.
pub fn apply_good_junction_gate_cjuncs(
    cjunctions: &mut Vec<CJunction>,
    bpcov: &BpcovStranded,
    refstart: u64,
    junction_thr: f64,
    long_intron: u64,
    isofrac_percent: f64,
    longreads: bool,
    eonly: bool,
) {
    good_junc(
        cjunctions,
        bpcov,
        refstart,
        junction_thr,
        long_intron,
        isofrac_percent,
        longreads,
        eonly,
    );
}

/// Back-compat alias for map-based function name.
pub fn apply_good_junction_gate(
    stats: &mut JunctionStats,
    bpcov: &BpcovStranded,
    refstart: u64,
    junction_thr: f64,
    long_intron: u64,
    isofrac_percent: f64,
    longreads: bool,
    eonly: bool,
) {
    good_junc_stats(
        stats,
        bpcov,
        refstart,
        junction_thr,
        long_intron,
        isofrac_percent,
        longreads,
        eonly,
    );
}
