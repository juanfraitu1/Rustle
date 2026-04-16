//! Coverage-based node splitting (find_trims_wsign + trimnode_all).
//! Splits graph nodes at positions where coverage drops significantly and adds
//! interior source/sink edges plus synthetic transfrags.

use crate::bpcov::Bpcov;
use crate::graph::{Graph, GraphTransfrag};
use crate::types::ReadBoundary;
use std::collections::VecDeque;

// header
const CHI_WIN: u64 = 100;
const CHI_THR: u64 = 50;
const DROP: f64 = 0.5;
const ERROR_PERC: f64 = 0.1;
const LONGINTRONANCHOR: u64 = 25;
const TRTHR: f64 = 1.0;

// helpers
#[allow(dead_code)]
fn compute_chi(winleft: &[f64], winright: &[f64], sumleft: f64, sumright: f64) -> f64 {
    let n = winleft.len().min(winright.len());
    if n == 0 {
        return 0.0;
    }
    let mut chi = 0.0;
    for j in 0..n {
        let denom = sumleft + sumright;
        if denom <= 0.0 {
            continue;
        }
        let mut mul = (winleft[j] + winright[j]) / denom;
        let mut mur = mul;
        mul *= sumleft;
        mur *= sumright;
        if mul > 0.0 {
            chi += (winleft[j] - mul) * (winleft[j] - mul) / mul;
        }
        if mur > 0.0 {
            chi += (winright[j] - mur) * (winright[j] - mur) / mur;
        }
    }
    chi
}

#[allow(dead_code)]
fn compute_chi2(winleft: &[f64], winright: &[f64], sumleft: f64, sumright: f64) -> f64 {
    let n = winleft.len().min(winright.len());
    if n == 0 {
        return 0.0;
    }
    let mut basecost = 0.0;
    let mut leftcost = 0.0;
    let mut rightcost = 0.0;
    let mubase = (sumleft + sumright) / (2.0 * n as f64);
    let muleft = sumleft / n as f64;
    let muright = sumright / n as f64;
    for j in 0..n {
        basecost += (winleft[j] - mubase) * (winleft[j] - mubase)
            + (winright[j] - mubase) * (winright[j] - mubase);
        leftcost += (winleft[j] - muleft) * (winleft[j] - muleft);
        rightcost += (winright[j] - muright) * (winright[j] - muright);
    }
    basecost / (2.0 * n as f64) - (leftcost + rightcost) / (n as f64)
}

#[allow(dead_code)]
fn compute_cost_range(
    wincov: &[f64],
    sumleft: f64,
    sumright: f64,
    leftstart: usize,
    rightstart: usize,
    rightend: usize,
) -> f64 {
    if rightend < leftstart || rightstart <= leftstart || rightstart > rightend {
        return 0.0;
    }
    let len = (rightend - leftstart + 1) as f64;
    let leftlen = (rightstart - leftstart) as f64;
    let rightlen = (rightend - rightstart + 1) as f64;
    if leftlen <= 0.0 || rightlen <= 0.0 {
        return 0.0;
    }
    let mubase = (sumleft + sumright) / len;
    let muleft = sumleft / leftlen;
    let muright = sumright / rightlen;
    let mut basecost = 0.0;
    let mut leftcost = 0.0;
    let mut rightcost = 0.0;
    for j in leftstart..rightstart {
        if j >= wincov.len() {
            break;
        }
        basecost += (wincov[j] - mubase).abs();
        leftcost += (wincov[j] - muleft).abs();
    }
    for j in rightstart..=rightend {
        if j >= wincov.len() {
            break;
        }
        basecost += (wincov[j] - mubase).abs();
        rightcost += (wincov[j] - muright).abs();
    }
    basecost / len - leftcost / leftlen - rightcost / rightlen
}

#[allow(dead_code)]
fn compute_cost_split(wincov: &[f64], split: usize, sumleft: f64, sumright: f64) -> f64 {
    let len = wincov.len();
    if len == 0 || split == 0 || split >= len {
        return 0.0;
    }
    let mut basecost = 0.0;
    let mut leftcost = 0.0;
    let mut rightcost = 0.0;
    let mubase = (sumleft + sumright) / len as f64;
    let muleft = sumleft / split as f64;
    let muright = sumright / (len - split) as f64;
    for (j, &v) in wincov.iter().enumerate() {
        basecost += (v - mubase) * (v - mubase);
        if j < split {
            leftcost += (v - muleft) * (v - muleft);
        } else {
            rightcost += (v - muright) * (v - muright);
        }
    }
    basecost / len as f64 - leftcost / split as f64 - rightcost / (len - split) as f64
}

#[derive(Debug, Clone)]
struct TrimPoint {
    pos: u64,
    start: bool,
    abundance: f64,
    feature: bool,
}

fn push_or_replace(
    trimpoints: &mut Vec<TrimPoint>,
    tp: TrimPoint,
    lastdrop: &mut f64,
    thisdrop: f64,
) {
    if trimpoints.is_empty() {
        trimpoints.push(tp);
        *lastdrop = thisdrop;
        return;
    }
    let last = trimpoints.last().unwrap();
    if (!last.start && tp.start && tp.pos.saturating_sub(last.pos) > CHI_THR)
        || (last.start && !tp.start && tp.pos.saturating_sub(last.pos) > CHI_THR)
    {
        trimpoints.push(tp);
        *lastdrop = thisdrop;
    } else if thisdrop < *lastdrop {
        if let Some(last_mut) = trimpoints.last_mut() {
            *last_mut = tp;
            *lastdrop = thisdrop;
        }
    } else if (thisdrop - *lastdrop).abs() <= f64::EPSILON {
        // Tie-break by preserving larger abundance trimpoint.
        if let Some(last_mut) = trimpoints.last_mut() {
            if tp.abundance > last_mut.abundance {
                *last_mut = tp;
            }
        }
    }
}

fn multi_trim_discovery(
    bpcov: &Bpcov,
    start: u64,
    end: u64,
    mixed_mode: bool,
    feature_points: &[(u64, bool, f64)],
) -> Vec<TrimPoint> {
    let len = end.saturating_sub(start).saturating_add(1);
    if len < CHI_THR {
        return Vec::new();
    }
    let mut trimpoints: Vec<TrimPoint> = Vec::new();
    let mut localdrop = ERROR_PERC / DROP;
    let mut lastdrop = localdrop;

    if len < 2 * (CHI_WIN + CHI_THR) + 1 {
        if len < CHI_WIN {
            localdrop = ERROR_PERC / (10.0 * DROP);
        } else {
            localdrop = ERROR_PERC;
        }
        lastdrop = localdrop;
        let left = start.saturating_add(LONGINTRONANCHOR);
        let right = end.saturating_sub(LONGINTRONANCHOR);
        for i in left..right {
            let covleft = bpcov.get_cov_range(bpcov.idx(start), bpcov.idx(i))
                / (i.saturating_sub(start).max(1) as f64);
            let covright = bpcov.get_cov_range(bpcov.idx(i), bpcov.idx(end).saturating_add(1))
                / (end.saturating_sub(i).saturating_add(1).max(1) as f64);
            if covleft < covright && covright > 0.0 {
                let thisdrop = covleft / covright;
                if thisdrop < localdrop {
                    push_or_replace(
                        &mut trimpoints,
                        TrimPoint {
                            pos: i + 1,
                            start: true,
                            abundance: (covright - covleft) / DROP,
                            feature: false,
                        },
                        &mut lastdrop,
                        thisdrop,
                    );
                }
            } else if covright < covleft && covleft > 0.0 {
                let thisdrop = covright / covleft;
                if thisdrop < localdrop {
                    push_or_replace(
                        &mut trimpoints,
                        TrimPoint {
                            pos: i,
                            start: false,
                            abundance: (covleft - covright) / DROP,
                            feature: false,
                        },
                        &mut lastdrop,
                        thisdrop,
                    );
                }
            }
        }
    } else {
        let winlen = CHI_WIN + CHI_THR;
        let phase1_end = start + winlen - 1;
        for i in (start + CHI_THR - 1)..phase1_end {
            let covleft = bpcov.get_cov_range(bpcov.idx(start), bpcov.idx(i).saturating_add(1))
                / (i.saturating_sub(start).max(1) as f64);
            let covright = bpcov
                .get_cov_range(bpcov.idx(i + 1), bpcov.idx(i + winlen).saturating_add(1))
                / winlen as f64;
            if covleft < covright && covright > 0.0 {
                let thisdrop = covleft / covright;
                if thisdrop < localdrop {
                    push_or_replace(
                        &mut trimpoints,
                        TrimPoint {
                            pos: i + 1,
                            start: true,
                            abundance: (covright - covleft) / DROP,
                            feature: false,
                        },
                        &mut lastdrop,
                        thisdrop,
                    );
                }
            } else if covright < covleft && covleft > 0.0 {
                let thisdrop = covright / covleft;
                if thisdrop < localdrop {
                    push_or_replace(
                        &mut trimpoints,
                        TrimPoint {
                            pos: i,
                            start: false,
                            abundance: (covleft - covright) / DROP,
                            feature: false,
                        },
                        &mut lastdrop,
                        thisdrop,
                    );
                }
            }
        }

        localdrop = if mixed_mode {
            DROP * DROP
        } else {
            ERROR_PERC / DROP
        };
        if trimpoints.is_empty() {
            lastdrop = localdrop;
        }
        let mid_end = end.saturating_sub(winlen);
        for i in phase1_end..mid_end {
            let covleft = bpcov.get_cov_range(
                bpcov.idx(i.saturating_sub(winlen).saturating_add(1)),
                bpcov.idx(i).saturating_add(1),
            );
            let covright =
                bpcov.get_cov_range(bpcov.idx(i + 1), bpcov.idx(i + winlen).saturating_add(1));
            if covleft < covright && covright > 0.0 {
                let thisdrop = covleft / covright;
                if thisdrop < localdrop {
                    push_or_replace(
                        &mut trimpoints,
                        TrimPoint {
                            pos: i + 1,
                            start: true,
                            abundance: (covright - covleft) / (DROP * winlen as f64),
                            feature: false,
                        },
                        &mut lastdrop,
                        thisdrop,
                    );
                }
            } else if covright < covleft && covleft > 0.0 {
                let thisdrop = covright / covleft;
                if thisdrop < localdrop {
                    push_or_replace(
                        &mut trimpoints,
                        TrimPoint {
                            pos: i,
                            start: false,
                            abundance: (covleft - covright) / (DROP * winlen as f64),
                            feature: false,
                        },
                        &mut lastdrop,
                        thisdrop,
                    );
                }
            }
        }

        localdrop = ERROR_PERC;
        for i in mid_end..end.saturating_sub(CHI_THR + 1) {
            let covleft = bpcov.get_cov_range(
                bpcov.idx(i.saturating_sub(winlen).saturating_add(1)),
                bpcov.idx(i).saturating_add(1),
            ) / winlen as f64;
            let covright = bpcov.get_cov_range(bpcov.idx(i + 1), bpcov.idx(end).saturating_add(1))
                / (end.saturating_sub(i).max(1) as f64);
            if covleft < covright && covright > 0.0 {
                let thisdrop = covleft / covright;
                if thisdrop < localdrop {
                    push_or_replace(
                        &mut trimpoints,
                        TrimPoint {
                            pos: i + 1,
                            start: true,
                            abundance: (covright - covleft) / DROP,
                            feature: false,
                        },
                        &mut lastdrop,
                        thisdrop,
                    );
                }
            } else if covright < covleft && covleft > 0.0 {
                let thisdrop = covright / covleft;
                if thisdrop < localdrop {
                    push_or_replace(
                        &mut trimpoints,
                        TrimPoint {
                            pos: i,
                            start: false,
                            abundance: (covleft - covright) / DROP,
                            feature: false,
                        },
                        &mut lastdrop,
                        thisdrop,
                    );
                }
            }
        }
    }

    // Inject feature trim points (tstartend behavior) and protect them from global elimination.
    for &(pos, is_start, cov) in feature_points {
        if pos <= start || pos >= end {
            continue;
        }
        trimpoints.push(TrimPoint {
            pos,
            start: is_start,
            abundance: (cov.abs() + 1.0) / DROP,
            feature: true,
        });
    }
    trimpoints.sort_unstable_by_key(|t| t.pos);

    // global threshold gate with backward invalidation of earlier points
    // when a downstream point fails.
    let localdrop_global = ERROR_PERC / DROP;
    let mut keep: Vec<bool> = vec![true; trimpoints.len()];
    let mut laststart = start;

    for i in 0..trimpoints.len() {
        let tp = &trimpoints[i];
        if tp.feature {
            laststart = if tp.start {
                tp.pos
            } else {
                tp.pos.saturating_add(1)
            };
            continue;
        }
        if tp.pos <= start || tp.pos >= end {
            keep[i] = false;
            continue;
        }

        let mid = if tp.start {
            tp.pos.saturating_sub(1)
        } else {
            tp.pos
        };
        let mut leftstart = laststart;
        if mid.saturating_sub(leftstart) > 3 * CHI_WIN {
            leftstart = mid.saturating_sub(3 * CHI_WIN);
        }
        let mut endpos = if i + 1 < trimpoints.len() {
            trimpoints[i + 1].pos
        } else {
            end
        };
        if endpos.saturating_sub(mid) > 3 * CHI_WIN {
            endpos = mid.saturating_add(3 * CHI_WIN);
        }
        if endpos <= mid {
            keep[i] = false;
            continue;
        }
        let covleft = bpcov.get_cov_range(bpcov.idx(leftstart), bpcov.idx(mid).saturating_add(1))
            / (mid.saturating_sub(leftstart).saturating_add(1) as f64);
        let covright = bpcov.get_cov_range(
            bpcov.idx(mid.saturating_add(1)),
            bpcov.idx(endpos).saturating_add(1),
        ) / (endpos.saturating_sub(mid) as f64);

        let pass = (tp.start && covleft < localdrop_global * covright)
            || (covright < localdrop_global * covleft);
        if pass {
            laststart = tp.pos.saturating_add(1);
            continue;
        }

        keep[i] = false;
        // Backward re-evaluation (skip feature points) as in reference impl.
        let mut k = i as isize - 1;
        while k >= 0 {
            let ku = k as usize;
            if !keep[ku] {
                k -= 1;
                continue;
            }
            if trimpoints[ku].feature {
                break;
            }

            let k_tp = &trimpoints[ku];
            let k_mid = if k_tp.start {
                k_tp.pos.saturating_sub(1)
            } else {
                k_tp.pos
            };
            let mut left = start;
            let mut j = k - 1;
            while j >= 0 {
                let ju = j as usize;
                if keep[ju] {
                    left = trimpoints[ju].pos;
                    break;
                }
                j -= 1;
            }
            if k_mid.saturating_sub(left) > 3 * CHI_WIN {
                left = k_mid.saturating_sub(3 * CHI_WIN);
            }
            let mut right = endpos;
            if right.saturating_sub(k_mid) > 3 * CHI_WIN {
                right = k_mid.saturating_add(3 * CHI_WIN);
            }
            if right <= k_mid {
                keep[ku] = false;
                k -= 1;
                continue;
            }
            let lcov = bpcov.get_cov_range(bpcov.idx(left), bpcov.idx(k_mid).saturating_add(1))
                / (k_mid.saturating_sub(left).saturating_add(1) as f64);
            let rcov = bpcov.get_cov_range(
                bpcov.idx(k_mid.saturating_add(1)),
                bpcov.idx(right).saturating_add(1),
            ) / (right.saturating_sub(k_mid) as f64);
            let k_pass =
                (k_tp.start && lcov < localdrop_global * rcov) || (rcov < localdrop_global * lcov);
            if k_pass {
                break;
            }
            keep[ku] = false;
            k -= 1;
        }
    }

    trimpoints
        .into_iter()
        .enumerate()
        .filter_map(|(i, tp)| if keep[i] { Some(tp) } else { None })
        .collect()
}

fn find_all_trims(
    bpcov: &Bpcov,
    start: u64,
    end: u64,
    mixed_mode: bool,
    feature_points: &[(u64, bool, f64)],
) -> Vec<TrimPoint> {
    multi_trim_discovery(bpcov, start, end, mixed_mode, feature_points)
}

/// Public wrapper for ref_map: (sourcestart, sinkend, source_abundance, sink_abundance).
pub fn find_trims_wsign_export(
    bpcov: &Bpcov,
    start: u64,
    end: u64,
) -> (Option<u64>, Option<u64>, f64, f64) {
    let trims = find_all_trims(bpcov, start, end, false, &[]);
    let mut ss = None;
    let mut se = None;
    let mut sa = 0.0;
    let mut ska = 0.0;
    for t in trims {
        if t.start && ss.is_none() {
            ss = Some(t.pos);
            sa = t.abundance;
        } else if !t.start && se.is_none() {
            se = Some(t.pos);
            ska = t.abundance;
        }
        if ss.is_some() && se.is_some() {
            break;
        }
    }
    (ss, se, sa, ska)
}

/// trimnode_all (trimnode_all + original algorithm apply_trimnode_all): for each node (length >= CHI_THR),
/// run find_trims_wsign; if source and/or sink trim, split node, add edges and synthetic transfrags.
/// Returns synthetic transfrags to be appended to the main transfrag list.
/// Does NOT set hardstart/hardend on coverage-trimmed nodes (interior edges only; original algorithm V68).
pub fn apply_coverage_trim(
    graph: &mut Graph,
    bpcov: &Bpcov,
    lstart: &[ReadBoundary],
    lend: &[ReadBoundary],
    _bundle_start: u64,
    _strand: char,
    min_abundance: f64,
    mixed_mode: bool,
    verbose: bool,
) -> Vec<GraphTransfrag> {
    let mut synthetic = Vec::new();
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let original_n = graph.n_nodes;
    // trimnode_all: newly split nodes can expose additional split points.
    // Process a dynamic queue so one original node may be split multiple times in one pass.
    let mut nodes_to_process: VecDeque<usize> = (1..original_n)
        .filter(|&i| i != source_id && i != sink_id)
        .collect();
    let mut split_ops = 0usize;
    let split_budget = original_n.saturating_mul(8).max(64);

    while let Some(nid) = nodes_to_process.pop_front() {
        let node_len = match graph.nodes.get(nid) {
            Some(n) => n.end.saturating_sub(n.start),
            None => continue,
        };
        if node_len < CHI_THR {
            continue;
        }

        let node = graph.nodes[nid].clone();
        let (start, end) = (node.start, node.end);
        let mut features: Vec<(u64, bool, f64)> = Vec::new();
        for rb in lstart {
            if rb.pos > start && rb.pos < end {
                features.push((rb.pos, true, rb.cov));
            }
        }
        for rb in lend {
            if rb.pos > start && rb.pos < end {
                features.push((rb.pos, false, rb.cov));
            }
        }
        let trims = find_all_trims(bpcov, start, end, mixed_mode, &features);
        if trims.is_empty() {
            continue;
        }

        let mut current = nid;
        for tp in trims {
            let mut ab = tp.abundance + TRTHR;
            if ab < min_abundance {
                ab = min_abundance;
            }
            let (cur_start, cur_end) = {
                let n = &graph.nodes[current];
                (n.start, n.end)
            };
            if tp.pos <= cur_start || tp.pos >= cur_end {
                continue;
            }

            let children: Vec<usize> = graph.nodes[current].children.ones().collect();
            for c in children.iter().copied() {
                graph.remove_edge(current, c);
            }

            if tp.start {
                graph.nodes[current].end = tp.pos - 1;
                let new_id = {
                    let n = graph.add_node(tp.pos, cur_end);
                    n.node_id
                };
                graph.add_edge(source_id, new_id);
                graph.add_edge(current, new_id);
                for c in children.iter().copied() {
                    graph.add_edge(new_id, c);
                }
                let mut stf = GraphTransfrag::new(vec![source_id, new_id], graph.pattern_size());
                stf.abundance = ab;
                stf.longstart = tp.pos;
                stf.longend = tp.pos + 1;
                synthetic.push(stf);
                current = new_id;
                split_ops += 1;
                if split_ops < split_budget {
                    nodes_to_process.push_back(new_id);
                }
            } else {
                graph.nodes[current].end = tp.pos;
                let new_id = {
                    let n = graph.add_node(tp.pos + 1, cur_end);
                    n.node_id
                };
                graph.add_edge(current, new_id);
                for c in children.iter().copied() {
                    graph.add_edge(new_id, c);
                }
                // futuretr: avoid redundant near-sink link creation.
                if !graph.has_sink_within_anchor_downstream(new_id, LONGINTRONANCHOR) {
                    graph.add_edge(current, sink_id);
                    let mut stf = GraphTransfrag::new(vec![current, sink_id], graph.pattern_size());
                    stf.abundance = ab;
                    stf.longstart = cur_start;
                    stf.longend = tp.pos + 1;
                    synthetic.push(stf);
                }
                current = new_id;
                split_ops += 1;
                if split_ops < split_budget {
                    nodes_to_process.push_back(new_id);
                }
            }
        }
    }

    if !synthetic.is_empty() {
        graph.compute_reachability();
    }
    if verbose && !synthetic.is_empty() {
        eprintln!(
            "      [Rustle] coverage_trim: {} synthetic transfrags from node splits",
            synthetic.len()
        );
    }
    synthetic
}
