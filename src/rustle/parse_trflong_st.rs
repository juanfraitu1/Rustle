//! Byte-faithful port of StringTie's `parse_trflong` and `update_abundance`.
//!
//! Goal: produce a Rust translation of the corresponding rlink.cpp logic that,
//! given identical inputs (graph, transfrags, gpos, no2gnode, nodecov, etc.),
//! emits the same predictions and updates state identically. Once verified,
//! flip the default behind `RUSTLE_PARSE_TRFLONG_ST=1`.
//!
//! ## Status
//!
//! **Scaffold only.** None of the slices below are implemented yet; they are
//! stubs that fall back to the existing rustle path when not enabled. This
//! file documents ST's structure so the port can proceed slice-by-slice.
//!
//! ## ST source map (rlink.cpp)
//!
//! | Slice | rlink.cpp range | Lines | Description |
//! |---|---|---|---|
//! | `parse_trflong` (main loop) | 10080–10518 | 439 | iterates trflong[] backward, extracts paths, emits CPredictions |
//! | `update_abundance` | 4681–4870 | 189 | adjusts long-read transfrag nodes; merges or creates transfrag in tr2no |
//! | `back_to_source_fast_long` | ~9540 | ~150 | walk pathpat back to source via best parents |
//! | `fwd_to_sink_fast_long` | ~9700 | ~150 | walk forward to sink |
//! | `long_max_flow` | ~7700 | ~200 | per-path flow assignment with capacity subtraction |
//! | `has_lr_witness_two_splices` | ~9920 | ~50 | LR witness gate for splice consecutivity |
//!
//! ## ST's parse_trflong main-loop pseudocode (paraphrased)
//!
//! ```text
//! for f in (trflong.len()-1 .. 0):           // REVERSE order
//!     t = trflong[f]
//!     if nasc-mode-mismatch: continue
//!     longcov = transfrag[t].abundance
//!     path = []
//!     pathpat = transfrag[t].pattern
//!     mark source/sink edges in pathpat
//!     maxi = transfrag[t].nodes[0]
//!     path.push(maxi)
//!     pathpat[maxi] = 1
//!     istranscript = empty
//!     istranscript[t] = 1
//!     flux = 0
//!     nodeflux = []
//!     tocheck = true
//!     if back_to_source_fast_long(maxi, path, ...):
//!         path.push(0); path.reverse()
//!         if fwd_to_sink_fast_long(maxi, path, ...):
//!             thardstart = no2gnode[transfrag[t].nodes[0]].hardstart
//!             thardend = no2gnode[transfrag[t].nodes.last()].hardend
//!             startnode = 1; lastnode = path.len() - 2
//!             checkpath = true
//!             // poly_end_unaligned trim (per-strand path adjustment) ...
//!             // hardstart/hardend boundary trim ...
//!             // LR-witness-of-consecutive-splices gate ...
//!             if checkpath and lastnode >= startnode:
//!                 if !lr_witness_holds: continue
//!             flux = long_max_flow(...)
//!             if flux > 0 or transfrag[t].guide:
//!                 tocheck = false
//!                 build exons from path
//!                 store CPrediction
//!                 update_abundance(...) for each consumed transfrag
//!     if tocheck:
//!         conditional checktrf.Add(t)
//! ```
//!
//! ## Comparison harness
//!
//! When `RUSTLE_PARSE_TRFLONG_ST_COMPARE=1` is set alongside the existing
//! pipeline, both old and new paths run. The old path is canonical (we use
//! its output). The new path's intermediate decisions (seed order, path
//! built, flux computed, prediction span) are emitted as
//! `parse_trflong_st_*` parity_decision events for diff against the old
//! path's `path_emit` events. Once equality is established for a slice, we
//! flip the slice's gate to use the new code as canonical.

#![allow(dead_code)]

use crate::graph::{Graph, GraphTransfrag};

/// Top-level switch: when true, the ST-faithful path is used as canonical.
/// Until we've verified slice-by-slice equivalence, this remains false.
pub fn enabled() -> bool {
    std::env::var_os("RUSTLE_PARSE_TRFLONG_ST").is_some()
}

/// Compare-mode: run both old and new paths, log the new path's decisions
/// as parity_decisions events for offline diff. Old path is canonical.
pub fn compare_mode() -> bool {
    std::env::var_os("RUSTLE_PARSE_TRFLONG_ST_COMPARE").is_some()
}

// ───────────────────────────────────────────────────────────────────────────
// Slice 1: seed iteration order
// ───────────────────────────────────────────────────────────────────────────
//
// ST's `parse_trflong` iterates `trflong[]` in REVERSE: `for f in count-1..0`.
// The `trflong` vector itself is built earlier in `process_transfrags`
// (rlink.cpp ~5970+) sorted by abundance and node count (descending). So
// the iteration order is: lowest-abundance / shortest first → highest-
// abundance / longest last. Why backward? Because ST ADDS to trflong over
// time (via `update_abundance` creating new transfrags); reversing visits
// the original entries in the same order they were ADDED.
//
// Rustle's current `parse_trflong` (in path_extract.rs ~5425) just sorts
// the seed indices and returns them. The order depends on a knob:
//   - default: sort by `usepath` (insertion order)
//   - RUSTLE_PARSE_TRFLONG_STRINGTIE_SORT: sort by abundance DESC, then
//     node count DESC
//
// **Slice 1 spec.** Produce a seed iteration order that matches ST's
// `for f=count-1..0` over the post-process_transfrags trflong[]. This
// requires:
//   1. Reproduce ST's trflong[] construction order (sort criteria + tie-
//      breaks).
//   2. Iterate the resulting Vec backward.
//
// Until then, this fn is a stub that falls back to the existing impl.

/// ST-faithful seed iteration order (slice 1, NOT YET IMPLEMENTED).
///
/// Returns Some(order) when ST mode is enabled and the slice is complete;
/// None to indicate the caller should use the existing rustle order.
pub fn seed_order_st(transfrags: &[GraphTransfrag]) -> Option<Vec<usize>> {
    if !enabled() {
        return None;
    }
    let _ = transfrags;
    // TODO: implement ST's trflong[] construction order, then return its
    // reverse iteration as Vec<usize>.
    None
}

// ───────────────────────────────────────────────────────────────────────────
// Slice 2: back_to_source / fwd_to_sink path extension
// ───────────────────────────────────────────────────────────────────────────
//
// ST source: `fwd_to_sink_fast_long` (rlink.cpp:8051-8259, ~208 lines)
//            `back_to_source_fast_long` (rlink.cpp:8261-8482, ~221 lines)
//            helpers: `replace_transfrag` (rlink.cpp:8025-8049),
//                     `onpath_long` (rlink.cpp:7979 — already ported)
//
// Rustle has bloated equivalents in path_extract.rs (3735-4561 + 4563-5415
// = 1678 lines, 4× ST). The bloat is from accumulated env-knob features.
//
// **Approach.** Port ST's lean version into this module. Run side-by-side
// with rustle's existing fn via `RUSTLE_PARSE_TRFLONG_ST_COMPARE=1`. Each
// call logs (input_seed, output_path) on both sides. Once equivalent on
// chr19, flip the canonical to ST-faithful via `RUSTLE_PARSE_TRFLONG_ST=1`.
//
// **Status (2026-05-05): partial port.**
//   - `replace_transfrag_st`: ported (~30 lines).
//   - `fwd_to_sink_fast_long_st`: signature in place; body TODO.
//   - `back_to_source_fast_long_st`: not yet started.
//   - Comparison harness: not yet wired.

use crate::bitvec::GBitVec;
use std::cell::RefCell;

thread_local! {
    /// Reusable pathpat scratch buffer for the comparison harness.
    /// Allocated once per thread, reset in place per call to avoid
    /// the per-call GBitVec.clone() that drove RSS to 100GB+ on chr19
    /// when running with `RUSTLE_PARSE_TRFLONG_ST_GATE=1`.
    static SCRATCH_PATHPAT: RefCell<GBitVec> = RefCell::new(GBitVec::default());
    /// Reusable path scratch (Vec<usize>) — same purpose.
    static SCRATCH_PATH: RefCell<Vec<usize>> = RefCell::new(Vec::with_capacity(256));
}

const DROP: f64 = 0.5;
const ERROR_PERC: f64 = 0.1;

/// Edge bit accessor mirroring ST's `gpos[edge(from,to,gno)]` lookup +
/// pathpat[edge_pos] read. Returns true when the edge bit is set in pathpat.
fn edge_pat_bit(graph: &Graph, pathpat: &GBitVec, from: usize, to: usize) -> bool {
    match graph.edge_bit_index(from, to) {
        Some(idx) => pathpat.get_bit(idx),
        None => false,
    }
}

fn edge_pat_set(graph: &Graph, pathpat: &mut GBitVec, from: usize, to: usize, value: bool) {
    if let Some(idx) = graph.edge_bit_index(from, to) {
        if value {
            pathpat.set_bit(idx);
        } else {
            pathpat.clear_bit(idx);
        }
    }
}

/// Port of StringTie's `replace_transfrag` (rlink.cpp:8025-8049).
///
/// Picks the better candidate between `t` and `tmax` (a -1 sentinel for
/// "no current pick"). "Better" = non-weak preferred; among equal weak
/// status, higher abundance preferred. Mutates `tmax` in place.
///
/// **Rustle weak semantics differ from ST.** ST's `weak` is `int` with
/// tri-state {-1 unset, 0 strong, 1 weak}; the -1 triggers a lazy
/// `compute_weak` resolution. Rustle's `weak` is `u8` with binary
/// {0 strong, 1 weak}; there is no unset state, since rustle marks
/// `weak=1` for ABSORBED transfrags during process_transfrags and
/// otherwise leaves it at 0. We treat `weak==0` as ST's "strong" and
/// `weak!=0` as ST's "weak", skipping the lazy compute branch.
pub fn replace_transfrag_st(
    t: i32,
    tmax: &mut i32,
    transfrags: &[GraphTransfrag],
) {
    if *tmax == -1 {
        *tmax = t;
        return;
    }
    let tmax_idx = *tmax as usize;
    let t_idx = t as usize;
    let tmax_weak = transfrags[tmax_idx].weak;
    let t_weak = transfrags[t_idx].weak;
    let tmax_ab = transfrags[tmax_idx].abundance;
    let t_ab = transfrags[t_idx].abundance;
    if tmax_weak == 0 {
        if t_weak == 0 && t_ab > tmax_ab {
            *tmax = t;
        }
    } else {
        if t_weak == 0 || t_ab > tmax_ab {
            *tmax = t;
        }
    }
}

/// ST-faithful `fwd_to_sink_fast_long` (rlink.cpp:8051-8259).
///
/// Walks forward from node `start_i` to sink, selecting the highest-cov
/// child at each step subject to pathpat compatibility, reachability of
/// `maxpath`, and an exclusion rule for adjacent low-cov children.
///
/// Returns true when sink is reached, false otherwise. Mutates `path`,
/// `pathpat`, `minpath`, `maxpath` in place.
///
/// **Differences from C++ ST**:
///   - Loop instead of tail recursion (semantically equivalent).
///   - Skips the `cnode->trf.Delete(j)` in-place optimization (depleted
///     transfrags still get filtered by the `abundance < epsilon` test
///     each iteration; the deletion is just a perf optimization for
///     subsequent calls). May be added later if it matters for parity.
///   - Uses rustle's bitset-iter order for children (= ascending node
///     ID) which mirrors ST's `inode->child[c]` order in normal flow
///     graphs (children appended in coordinate order during build).
#[allow(clippy::too_many_arguments)]
pub fn fwd_to_sink_fast_long_st(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    nodecov: &[f64],
    pathpat: &mut GBitVec,
    path: &mut Vec<usize>,
    minpath: &mut usize,
    maxpath: &mut usize,
    start_i: usize,
) -> bool {
    const EPS: f64 = 1e-9;
    let mut i = start_i;
    let sink_id = graph.sink_id;
    let _ = nodecov;
    let max_iters = graph.nodes.len() + 4; // bounded by graph size; +4 slack
    let mut iters = 0usize;

    loop {
        iters += 1;
        if iters > max_iters {
            return false; // safety: prevent infinite loop on pathological inputs
        }
        let inode = match graph.nodes.get(i) {
            Some(n) => n,
            None => return false,
        };

        // ST line 8056: i<maxpath && !inode->childpat[maxpath] => unreachable
        if i < *maxpath {
            let reachable = inode
                .childpat
                .as_ref()
                .map(|cp| cp.contains(*maxpath))
                .unwrap_or(false);
            if !reachable {
                return false;
            }
        }

        // ST line 8068: no children => sink reached
        if inode.children.ones().next().is_none() {
            return true;
        }

        let mut maxc: i64 = -1;
        let mut maxcov: f64 = 0.0;
        let mut tmax: i32 = -1;
        let mut exclude = false;
        let mut nextnode: usize = 0;
        let mut reach = *maxpath <= i;

        // ST line 8081: pathpat[i+1] short-circuit
        if i + 1 < graph.nodes.len() && pathpat.contains(i + 1) {
            maxc = (i + 1) as i64;
            tmax = -1;
            reach = true;
        } else {
            for c in inode.children.ones() {
                let childonpath = pathpat.contains(c);

                // ST line 8092: pathpat[edge_pos]
                let pos = graph.edge_bit_index(i, c);
                if let Some(p) = pos {
                    if pathpat.get_bit(p) {
                        maxc = c as i64;
                        tmax = -1;
                        reach = true;
                        break;
                    }
                }

                // ST line 8101: reachability check via nextnode
                if *maxpath > i {
                    if nextnode == 0 {
                        let mut j = i + 2;
                        while j <= *maxpath {
                            if pathpat.contains(j) {
                                nextnode = j;
                                break;
                            }
                            j += 1;
                        }
                    }
                    if c != nextnode {
                        let cn_reach = graph
                            .nodes
                            .get(c)
                            .and_then(|n| n.childpat.as_ref())
                            .map(|cp| cp.contains(nextnode))
                            .unwrap_or(false);
                        if !cn_reach {
                            continue;
                        }
                    }
                    reach = true;
                }

                let mut childcov: f64 = 0.0;
                let mut tchild: i32 = -1;
                let endpath = (*maxpath).max(c);

                let cnode = match graph.nodes.get(c) {
                    Some(n) => n,
                    None => continue,
                };

                // ST line 8134: adjacent + low-cov-drop "exclude" rule
                let cnode_len = cnode.length() as f64;
                let inode_end = inode.end;
                let nodecov_i = nodecov[i];
                let nodecov_c = if c < nodecov.len() { nodecov[c] } else { 0.0 };
                if c == i + 1
                    && i < graph.nodes.len().saturating_sub(2)
                    && inode_end == cnode.start
                    && cnode_len > 0.0
                    && nodecov_c / cnode_len < 1000.0
                    && nodecov_i * (DROP + ERROR_PERC) > nodecov_c
                {
                    exclude = true;
                } else {
                    // Mark candidate node + edge in pathpat
                    pathpat.insert_grow(c);
                    if let Some(p) = pos {
                        pathpat.set_bit(p);
                    }

                    // Iterate transfrags through cnode
                    let trf_count = cnode.trf_ids.len();
                    for ti in 0..trf_count {
                        let t = cnode.trf_ids[ti];
                        if t >= transfrags.len() {
                            continue;
                        }
                        let abund = transfrags[t].abundance;
                        if abund < EPS {
                            continue; // depleted
                        }
                        let tf_longread = transfrags[t].longread;
                        if !tf_longread {
                            continue;
                        }
                        let tf_first = transfrags[t].node_ids.first().copied().unwrap_or(0);
                        let tf_last = transfrags[t].node_ids.last().copied().unwrap_or(0);
                        if tf_first == graph.source_id {
                            continue; // ST line 8150: nodes[0] != 0 (= not source)
                        }
                        if c == sink_id {
                            // child is sink: only accept transfrags ending at i
                            if tf_first == i && *maxpath <= i {
                                childcov += abund;
                                if tchild == -1
                                    || abund > transfrags[tchild as usize].abundance
                                {
                                    tchild = t as i32;
                                }
                            }
                        } else if tf_first <= i
                            && tf_last >= c
                            && crate::path_extract::onpath_long_pub(
                                &transfrags[t].pattern,
                                &transfrags[t].node_ids,
                                pathpat,
                                *minpath,
                                endpath,
                                graph,
                            )
                        {
                            childcov += abund;
                            replace_transfrag_st(t as i32, &mut tchild, transfrags);
                        }
                    }

                    // ST line 8168: pick best child by cov + weak preference
                    if childcov > maxcov {
                        let take = if tmax == -1 {
                            true
                        } else {
                            transfrags[tmax as usize].weak > 0
                                || (tchild >= 0 && transfrags[tchild as usize].weak == 0)
                        };
                        if take {
                            maxcov = childcov;
                            maxc = c as i64;
                            tmax = tchild;
                        }
                    } else if maxc != -1 && childcov > maxcov - EPS {
                        let take = if tmax == -1 {
                            tchild >= 0 && transfrags[tchild as usize].weak == 0
                        } else {
                            transfrags[tmax as usize].weak > 0
                                || (tchild >= 0 && transfrags[tchild as usize].weak == 0)
                        };
                        if take {
                            let nodecov_maxc =
                                if (maxc as usize) < nodecov.len() { nodecov[maxc as usize] } else { 0.0 };
                            if nodecov_maxc < nodecov_c {
                                maxc = c as i64;
                                tmax = tchild;
                            }
                        }
                    }

                    // Roll back pathpat marks (we'll re-set the winner later)
                    if let Some(p) = pos {
                        pathpat.clear_bit(p);
                    }
                    if childonpath {
                        break;
                    }
                    pathpat.remove(c);
                }
            }
        }

        if !reach {
            return false;
        }

        if maxc == -1 {
            // ST line 8195: try the excluded i+1 path if available
            if exclude && i + 1 < nodecov.len() && nodecov[i + 1] > 0.0 {
                let cidx = i + 1;
                let cnode = match graph.nodes.get(cidx) {
                    Some(n) => n,
                    None => return false,
                };
                let mut childcov: f64 = 0.0;
                let mut tmax_local: i32 = -1;
                pathpat.insert_grow(cidx);
                let pos = graph.edge_bit_index(i, cidx);
                if let Some(p) = pos {
                    pathpat.set_bit(p);
                }
                let endpath = (*maxpath).max(cidx);
                let trf_count = cnode.trf_ids.len();
                for ti in 0..trf_count {
                    let t = cnode.trf_ids[ti];
                    if t >= transfrags.len() {
                        continue;
                    }
                    let abund = transfrags[t].abundance;
                    if abund < EPS {
                        continue;
                    }
                    let tf_longread = transfrags[t].longread;
                    if !tf_longread {
                        continue;
                    }
                    let tf_first = transfrags[t].node_ids.first().copied().unwrap_or(0);
                    let tf_last = transfrags[t].node_ids.last().copied().unwrap_or(0);
                    if tf_first <= i
                        && tf_last >= cidx
                        && crate::path_extract::onpath_long_pub(
                            &transfrags[t].pattern,
                            &transfrags[t].node_ids,
                            pathpat,
                            *minpath,
                            endpath,
                            graph,
                        )
                    {
                        childcov += abund;
                        replace_transfrag_st(t as i32, &mut tmax_local, transfrags);
                    }
                }
                pathpat.remove(cidx);
                if let Some(p) = pos {
                    pathpat.clear_bit(p);
                }
                if childcov > 0.0 {
                    maxc = cidx as i64;
                    tmax = tmax_local;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }

        // Commit maxc to path
        let maxc_idx = maxc as usize;
        path.push(maxc_idx);
        pathpat.insert_grow(maxc_idx);
        if let Some(p) = graph.edge_bit_index(i, maxc_idx) {
            pathpat.set_bit(p);
        }
        if tmax >= 0 {
            // pathpat |= transfrag[tmax].pattern
            let tf_pat = transfrags[tmax as usize].pattern.clone();
            pathpat.union_with(&tf_pat);
            let tf_first = transfrags[tmax as usize]
                .node_ids
                .first()
                .copied()
                .unwrap_or(0);
            let tf_last = transfrags[tmax as usize]
                .node_ids
                .last()
                .copied()
                .unwrap_or(0);
            if tf_first < *minpath {
                *minpath = tf_first;
            }
            if tf_last > *maxpath {
                *maxpath = tf_last;
            }
        }

        // Tail recursion → loop
        if maxc_idx == sink_id {
            return true;
        }
        i = maxc_idx;
    }
}

/// ST-faithful `back_to_source_fast_long` (rlink.cpp:8261-8482).
///
/// Symmetric to `fwd_to_sink_fast_long_st` — walks from `start_i` back
/// to source by selecting the highest-cov parent at each step.
///
/// Returns true when source is reached, false otherwise. ST's behavior
/// for `if(maxp) path.Add(maxp)` (line 8467) translates to: only add
/// the chosen parent to `path` when it isn't `source_id` (matching
/// ST's `0 = source` convention via the truthy check).
#[allow(clippy::too_many_arguments)]
pub fn back_to_source_fast_long_st(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    nodecov: &[f64],
    pathpat: &mut GBitVec,
    path: &mut Vec<usize>,
    minpath: &mut usize,
    maxpath: &mut usize,
    start_i: usize,
) -> bool {
    const EPS: f64 = 1e-9;
    let mut i = start_i;
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let _ = nodecov;
    let max_iters = graph.nodes.len() + 4;
    let mut iters = 0usize;

    loop {
        iters += 1;
        if iters > max_iters {
            return false;
        }
        let inode = match graph.nodes.get(i) {
            Some(n) => n,
            None => return false,
        };

        // ST line 8267: minpath<i && !inode->parentpat[minpath] => unreachable
        if *minpath < i {
            let reachable = inode
                .parentpat
                .as_ref()
                .map(|pp| pp.contains(*minpath))
                .unwrap_or(false);
            if !reachable {
                return false;
            }
        }

        // ST line 8279: no parents => source reached
        if inode.parents.ones().next().is_none() {
            return true;
        }

        let mut maxp: i64 = -1;
        let mut maxcov: f64 = 0.0;
        let mut tmax: i32 = -1;
        let mut exclude = false;
        let n_nodes = graph.nodes.len();
        // ST sets `nextnode=gno` as sentinel "not found"; rustle uses n_nodes.
        let mut nextnode: usize = n_nodes;
        let mut reach = *minpath >= i;

        // ST line 8298: pathpat[i-1] short-circuit
        if i > 0 && pathpat.contains(i - 1) {
            maxp = (i - 1) as i64;
            tmax = -1;
            reach = true;
        } else {
            for p in inode.parents.ones() {
                let parentonpath = pathpat.contains(p);

                // ST line 8312: edge already in pathpat
                let pos = graph.edge_bit_index(p, i);
                if let Some(ep) = pos {
                    if pathpat.get_bit(ep) {
                        maxp = p as i64;
                        tmax = -1;
                        reach = true;
                        break;
                    }
                }

                // ST line 8321: reachability via nextnode
                if *minpath < i {
                    if nextnode == n_nodes {
                        // ST: j=i-2; while(j>=minpath) { if(pathpat[j]) {nextnode=j; break} j++ }
                        // Note ST's `j++` looks like a bug-but-faithful: should be j-- to walk
                        // backward. Following ST exactly preserves byte-faithfulness even if
                        // the loop never advances meaningfully. We'll diff against ST and
                        // adjust if needed.
                        let mut j = i.saturating_sub(2);
                        while j >= *minpath && j != usize::MAX {
                            if pathpat.contains(j) {
                                nextnode = j;
                                break;
                            }
                            j = j.wrapping_add(1); // ST's j++ — preserved verbatim
                            // Guard against runaway: bail if j exceeds bounds
                            if j > i || j == usize::MAX {
                                break;
                            }
                        }
                    }
                    if p != nextnode {
                        let pn_reach = graph
                            .nodes
                            .get(p)
                            .and_then(|n| n.parentpat.as_ref())
                            .map(|pp| pp.contains(nextnode))
                            .unwrap_or(false);
                        if !pn_reach {
                            continue;
                        }
                    }
                    reach = true;
                }

                let mut parentcov: f64 = 0.0;
                let mut tpar: i32 = -1;
                let startpath = (*minpath).min(p);

                let pnode = match graph.nodes.get(p) {
                    Some(n) => n,
                    None => continue,
                };

                // ST line 8355: adjacent + low-cov-drop "exclude" rule
                let pnode_len = pnode.length() as f64;
                let inode_start = inode.start;
                let nodecov_i = nodecov[i];
                let nodecov_p = if p < nodecov.len() { nodecov[p] } else { 0.0 };
                if p == i.saturating_sub(1)
                    && i > 1
                    && inode_start == pnode.end
                    && pnode_len > 0.0
                    && nodecov_p / pnode_len < 1000.0
                    && nodecov_i * (DROP + ERROR_PERC) > nodecov_p
                {
                    exclude = true;
                } else {
                    pathpat.insert_grow(p);
                    if let Some(ep) = pos {
                        pathpat.set_bit(ep);
                    }

                    let trf_count = pnode.trf_ids.len();
                    for ti in 0..trf_count {
                        let t = pnode.trf_ids[ti];
                        if t >= transfrags.len() {
                            continue;
                        }
                        let abund = transfrags[t].abundance;
                        if abund < EPS {
                            continue;
                        }
                        let tf_longread = transfrags[t].longread;
                        if !tf_longread {
                            continue;
                        }
                        let tf_first = transfrags[t].node_ids.first().copied().unwrap_or(0);
                        let tf_last = transfrags[t].node_ids.last().copied().unwrap_or(0);
                        // ST line 8371: skip transfrags that include sink
                        // (transfrag[t]->nodes.Last() != gno-1)
                        if tf_last == sink_id {
                            continue;
                        }
                        if p == source_id {
                            // ST line 8372: parent is source; need tf ending at i
                            if tf_last == i && *minpath >= i {
                                parentcov += abund;
                                if tpar == -1
                                    || abund > transfrags[tpar as usize].abundance
                                {
                                    tpar = t as i32;
                                }
                            }
                        } else if tf_first <= p
                            && tf_last >= i
                            && crate::path_extract::onpath_long_pub(
                                &transfrags[t].pattern,
                                &transfrags[t].node_ids,
                                pathpat,
                                startpath,
                                *maxpath,
                                graph,
                            )
                        {
                            parentcov += abund;
                            replace_transfrag_st(t as i32, &mut tpar, transfrags);
                        }
                    }

                    // ST line 8390: pick best parent by cov + weak preference
                    if parentcov > maxcov {
                        let take = if tmax == -1 {
                            true
                        } else {
                            transfrags[tmax as usize].weak > 0
                                || (tpar >= 0 && transfrags[tpar as usize].weak == 0)
                        };
                        if take {
                            maxcov = parentcov;
                            maxp = p as i64;
                            tmax = tpar;
                        }
                    } else if maxp != -1 && parentcov > maxcov - EPS {
                        let take = if tmax == -1 {
                            tpar >= 0 && transfrags[tpar as usize].weak == 0
                        } else {
                            transfrags[tmax as usize].weak > 0
                                || (tpar >= 0 && transfrags[tpar as usize].weak == 0)
                        };
                        if take {
                            let nodecov_maxp =
                                if (maxp as usize) < nodecov.len() { nodecov[maxp as usize] } else { 0.0 };
                            if nodecov_maxp < nodecov_p {
                                maxp = p as i64;
                                tmax = tpar;
                            }
                        }
                    }

                    if let Some(ep) = pos {
                        pathpat.clear_bit(ep);
                    }
                    if parentonpath {
                        break;
                    }
                    pathpat.remove(p);
                }
            }
        }

        if !reach {
            return false;
        }

        if maxp == -1 {
            // ST line 8417: try the excluded i-1 path
            if exclude && i > 0 && nodecov[i - 1] > 0.0 {
                let pidx = i - 1;
                let pnode = match graph.nodes.get(pidx) {
                    Some(n) => n,
                    None => return false,
                };
                let mut parentcov: f64 = 0.0;
                let mut tmax_local: i32 = -1;
                pathpat.insert_grow(pidx);
                let pos = graph.edge_bit_index(pidx, i);
                if let Some(ep) = pos {
                    pathpat.set_bit(ep);
                }
                let startpath = (*minpath).min(pidx);
                let trf_count = pnode.trf_ids.len();
                for ti in 0..trf_count {
                    let t = pnode.trf_ids[ti];
                    if t >= transfrags.len() {
                        continue;
                    }
                    let abund = transfrags[t].abundance;
                    if abund < EPS {
                        continue;
                    }
                    let tf_longread = transfrags[t].longread;
                    if !tf_longread {
                        continue;
                    }
                    let tf_first = transfrags[t].node_ids.first().copied().unwrap_or(0);
                    let tf_last = transfrags[t].node_ids.last().copied().unwrap_or(0);
                    if tf_first <= pidx
                        && tf_last >= i
                        && crate::path_extract::onpath_long_pub(
                            &transfrags[t].pattern,
                            &transfrags[t].node_ids,
                            pathpat,
                            startpath,
                            *maxpath,
                            graph,
                        )
                    {
                        parentcov += abund;
                        replace_transfrag_st(t as i32, &mut tmax_local, transfrags);
                    }
                }
                pathpat.remove(pidx);
                if let Some(ep) = pos {
                    pathpat.clear_bit(ep);
                }
                if parentcov > 0.0 {
                    maxp = pidx as i64;
                    tmax = tmax_local;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }

        let maxp_idx = maxp as usize;
        // ST line 8467: only add to path if not source.
        if maxp_idx != source_id {
            path.push(maxp_idx);
        }
        pathpat.insert_grow(maxp_idx);
        if let Some(ep) = graph.edge_bit_index(maxp_idx, i) {
            pathpat.set_bit(ep);
        }
        if tmax >= 0 {
            let tf_pat = transfrags[tmax as usize].pattern.clone();
            pathpat.union_with(&tf_pat);
            let tf_first = transfrags[tmax as usize]
                .node_ids
                .first()
                .copied()
                .unwrap_or(0);
            let tf_last = transfrags[tmax as usize]
                .node_ids
                .last()
                .copied()
                .unwrap_or(0);
            if tf_first < *minpath {
                *minpath = tf_first;
            }
            if tf_last > *maxpath {
                *maxpath = tf_last;
            }
        }

        if maxp_idx == source_id {
            return true;
        }
        i = maxp_idx;
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Slice 3: long_max_flow capacity assignment
// ───────────────────────────────────────────────────────────────────────────
//
// ST's `long_max_flow` distributes per-path flux subject to nodecov caps,
// then SUBTRACTS the consumed nodecov from the global nodecov vector. This
// is what makes subsequent seeds compete for diminishing capacity. The
// rustle equivalent is in `path_extract.rs` flow_flux + nodecov updates.

// stubs intentionally elided until slice 3 begins

// ───────────────────────────────────────────────────────────────────────────
// Slice 4: update_abundance node-trim + accumulate
// ───────────────────────────────────────────────────────────────────────────
//
// ST's `update_abundance` (rlink.cpp:4681) is the second target. For each
// long-read seed, it:
//   1. (is_lr branch) Trims the first/last nodes of the transfrag if the
//      read's start/end falls *inside* a flanking node and there's a
//      coverage drop signaling a different transcript boundary.
//   2. Looks up the (possibly trimmed) node sequence in `tr2no` (a per-
//      bucket tree-pattern index). If found → bump existing transfrag's
//      abundance. If not → create a new transfrag with abundance.
//
// Rustle's equivalent is scattered across `transfrag_process.rs` and
// `path_extract.rs`. A faithful port would centralize the trim + lookup
// logic here.

// stubs intentionally elided until slice 4 begins

// ───────────────────────────────────────────────────────────────────────────
// Comparison harness: run BOTH the existing rustle path-extension fns and
// the new ST-faithful ones on cloned inputs, log divergences for offline diff.
// ───────────────────────────────────────────────────────────────────────────

/// Outcome of one path-extension call.
#[derive(Debug, Clone, PartialEq)]
pub struct PathExtendOutcome {
    pub returned_true: bool,
    pub path_len: usize,
    pub path_first: Option<usize>,
    pub path_last: Option<usize>,
    pub minpath: usize,
    pub maxpath: usize,
    /// Hash of pathpat bits (cheap divergence signal).
    pub pathpat_hash: u64,
}

fn pathpat_hash(p: &GBitVec) -> u64 {
    use std::hash::Hasher;
    let mut h = std::collections::hash_map::DefaultHasher::new();
    let bits: Vec<usize> = p.ones().collect();
    h.write_usize(bits.len());
    for b in bits {
        h.write_usize(b);
    }
    h.finish()
}

/// Reset `dst` so it has the same set bits as `src`. Reuses `dst`'s
/// allocation; equivalent to `*dst = src.clone()` but avoids realloc
/// when `dst` already has enough capacity.
fn pathpat_copy_from(src: &GBitVec, dst: &mut GBitVec) {
    dst.clear();
    dst.union_with(src);
}

/// Run the ST-faithful `fwd_to_sink_fast_long_st` on cloned inputs and
/// produce an `PathExtendOutcome`. Caller is responsible for passing the
/// SAME initial state that rustle's existing fwd_to_sink_fast_long sees.
///
/// Uses thread-local scratch for `pathpat`/`path` to amortize the
/// allocation cost across calls (otherwise per-call GBitVec.clone() ran
/// rustle to 100GB+ RSS on chr19 with the gate enabled).
#[allow(clippy::too_many_arguments)]
pub fn run_fwd_to_sink_st_on_clone(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    nodecov: &[f64],
    pathpat: &GBitVec,
    path: &[usize],
    minpath: usize,
    maxpath: usize,
    i: usize,
) -> PathExtendOutcome {
    SCRATCH_PATHPAT.with(|spp| {
        SCRATCH_PATH.with(|sp| {
            let mut pp = spp.borrow_mut();
            let mut path_buf = sp.borrow_mut();
            pathpat_copy_from(pathpat, &mut pp);
            path_buf.clear();
            path_buf.extend_from_slice(path);
            let mut minp = minpath;
            let mut maxp = maxpath;
            let returned_true = fwd_to_sink_fast_long_st(
                graph,
                transfrags,
                nodecov,
                &mut pp,
                &mut path_buf,
                &mut minp,
                &mut maxp,
                i,
            );
            PathExtendOutcome {
                returned_true,
                path_len: path_buf.len(),
                path_first: path_buf.first().copied(),
                path_last: path_buf.last().copied(),
                minpath: minp,
                maxpath: maxp,
                pathpat_hash: pathpat_hash(&pp),
            }
        })
    })
}

/// Symmetric for `back_to_source_fast_long_st`.
#[allow(clippy::too_many_arguments)]
pub fn run_back_to_source_st_on_clone(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    nodecov: &[f64],
    pathpat: &GBitVec,
    path: &[usize],
    minpath: usize,
    maxpath: usize,
    i: usize,
) -> PathExtendOutcome {
    SCRATCH_PATHPAT.with(|spp| {
        SCRATCH_PATH.with(|sp| {
            let mut pp = spp.borrow_mut();
            let mut path_buf = sp.borrow_mut();
            pathpat_copy_from(pathpat, &mut pp);
            path_buf.clear();
            path_buf.extend_from_slice(path);
            let mut minp = minpath;
            let mut maxp = maxpath;
            let returned_true = back_to_source_fast_long_st(
                graph,
                transfrags,
                nodecov,
                &mut pp,
                &mut path_buf,
                &mut minp,
                &mut maxp,
                i,
            );
            PathExtendOutcome {
                returned_true,
                path_len: path_buf.len(),
                path_first: path_buf.first().copied(),
                path_last: path_buf.last().copied(),
                minpath: minp,
                maxpath: maxp,
                pathpat_hash: pathpat_hash(&pp),
            }
        })
    })
}

/// Capture rustle's outcome from the SAME state for comparison.
pub fn outcome_from_rustle(
    returned_true: bool,
    pathpat: &GBitVec,
    path: &[usize],
    minpath: usize,
    maxpath: usize,
) -> PathExtendOutcome {
    PathExtendOutcome {
        returned_true,
        path_len: path.len(),
        path_first: path.first().copied(),
        path_last: path.last().copied(),
        minpath,
        maxpath,
        pathpat_hash: pathpat_hash(pathpat),
    }
}

/// Emit a divergence event when the two outcomes differ.
pub fn emit_diff_if_diverges(
    rustle: &PathExtendOutcome,
    st: &PathExtendOutcome,
    direction: &str,
    seed_idx: usize,
    start_node: usize,
) {
    if rustle == st {
        return;
    }
    if !crate::parity_decisions::is_enabled() {
        return;
    }
    let r_last = rustle.path_last.map(|x| x as i64).unwrap_or(-1);
    let s_last = st.path_last.map(|x| x as i64).unwrap_or(-1);
    let payload = format!(
        r#""direction":"{}","seed":{},"start":{},"r_ret":{},"st_ret":{},"r_path_len":{},"st_path_len":{},"r_path_last":{},"st_path_last":{},"r_minp":{},"st_minp":{},"r_maxp":{},"st_maxp":{},"r_pphash":{},"st_pphash":{}"#,
        direction,
        seed_idx,
        start_node,
        rustle.returned_true,
        st.returned_true,
        rustle.path_len,
        st.path_len,
        r_last,
        s_last,
        rustle.minpath,
        st.minpath,
        rustle.maxpath,
        st.maxpath,
        rustle.pathpat_hash,
        st.pathpat_hash,
    );
    crate::parity_decisions::emit("path_extend_diff", None, 0, 0, '.', &payload);
}

/// Whether the comparison harness is enabled. Off by default; set
/// `RUSTLE_PARSE_TRFLONG_ST_COMPARE=1` to activate at every fwd_to_sink
/// and back_to_source callsite.
pub fn comparison_active() -> bool {
    std::env::var_os("RUSTLE_PARSE_TRFLONG_ST_COMPARE").is_some()
}

/// Canonical mode: use ST-faithful back_to_source + fwd_to_sink as the
/// actual path builders (not just as a gate/comparator).
///
/// When active, path_extract.rs skips rustle's own back/fwd calls entirely
/// and calls back_to_source_fast_long_st / fwd_to_sink_fast_long_st on the
/// real pathpat+path. This fixes same-seed-wrong-path cases where both tools
/// extend successfully but pick different branches.
///
/// Enable: `RUSTLE_PARSE_TRFLONG_ST_CANONICAL=1`
pub fn canonical_active() -> bool {
    std::env::var_os("RUSTLE_PARSE_TRFLONG_ST_CANONICAL").is_some()
}

// ───────────────────────────────────────────────────────────────────────────
// Tests for the ported helpers
// ───────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{GraphNode, GraphTransfrag};

    fn make_tf(node_ids: Vec<usize>, abund: f64, longread: bool, weak: u8) -> GraphTransfrag {
        let pat_size = 64;
        let mut tf = GraphTransfrag::new(node_ids, pat_size);
        tf.abundance = abund;
        tf.longread = longread;
        tf.weak = weak;
        tf
    }

    #[test]
    fn replace_transfrag_st_picks_first_when_tmax_is_neg1() {
        let transfrags = vec![make_tf(vec![1, 2], 5.0, true, 0)];
        let mut tmax: i32 = -1;
        replace_transfrag_st(0, &mut tmax, &transfrags);
        assert_eq!(tmax, 0);
    }

    // (rustle's weak is u8 binary; no unset/-1 case to resolve, unlike ST.)

    #[test]
    fn replace_transfrag_st_prefers_strong_over_weak_regardless_of_abundance() {
        let mut transfrags = vec![
            make_tf(vec![1, 2], 100.0, true, 1), // strong-zero=false (weak=1)
            make_tf(vec![1, 2], 1.0, true, 0),   // strong (weak=0)
        ];
        // Start with the weak one as tmax
        let mut tmax: i32 = 0;
        replace_transfrag_st(1, &mut tmax, &mut transfrags);
        // tmax was weak; t is strong → take t even though abundance is lower
        assert_eq!(tmax, 1);
    }

    #[test]
    fn replace_transfrag_st_among_strong_picks_higher_abundance() {
        let mut transfrags = vec![
            make_tf(vec![1, 2], 5.0, true, 0),
            make_tf(vec![1, 2], 10.0, true, 0),
        ];
        let mut tmax: i32 = 0;
        replace_transfrag_st(1, &mut tmax, &mut transfrags);
        assert_eq!(tmax, 1);
    }

    #[test]
    fn replace_transfrag_st_among_strong_keeps_tmax_when_abundance_lower() {
        let mut transfrags = vec![
            make_tf(vec![1, 2], 10.0, true, 0),
            make_tf(vec![1, 2], 5.0, true, 0),
        ];
        let mut tmax: i32 = 0;
        replace_transfrag_st(1, &mut tmax, &mut transfrags);
        assert_eq!(tmax, 0);
    }

    /// Trivial fwd_to_sink test: a 4-node graph (source → A → B → sink) with
    /// a long-read transfrag spanning A → B. The walk from A should pick B
    /// as the next child and continue to sink.
    #[test]
    fn fwd_to_sink_st_minimal_graph() {
        let mut graph = crate::graph::Graph::new();
        let _ = graph.add_node(0, 0);     // node 0 = source
        let _ = graph.add_node(100, 200); // node 1 = exon A
        let _ = graph.add_node(300, 400); // node 2 = exon B
        let _ = graph.add_node(0, 0);     // node 3 = sink
        graph.source_id = 0;
        graph.sink_id = 3;
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(2, 3);
        graph.compute_reachability();
        let pat_size = graph.pattern_size();

        // ST convention: long-read transfrag does NOT include source/sink.
        // Use [1, 2] as the internal path; pattern includes the edges to
        // source/sink as in ST's parse_trflong setup.
        let mut tf = GraphTransfrag::new(vec![1, 2], pat_size);
        tf.abundance = 10.0;
        tf.longread = true;
        tf.weak = 0;
        // Set node bits + edge bits in pattern (ST convention)
        tf.pattern.insert_grow(1);
        tf.pattern.insert_grow(2);
        graph.set_pattern_edges_for_path(&mut tf.pattern, &[1, 2]);
        let mut transfrags = vec![tf];
        graph.nodes[1].trf_ids.push(0);
        graph.nodes[2].trf_ids.push(0);

        let mut nodecov = vec![0.0; graph.nodes.len()];
        nodecov[1] = 10.0;
        nodecov[2] = 10.0;
        // Mirror ST's parse_trflong setup: pathpat starts as the seed
        // transfrag's pattern (includes all its node bits + edges), then
        // source-to-minp and maxp-to-sink edge bits get added.
        let mut pathpat = transfrags[0].pattern.clone();
        // ST line 10110-10114: add gpos[edge(0,minp)] and gpos[edge(maxp,gno-1)]
        if let Some(idx) = graph.edge_bit_index(graph.source_id, 1) {
            pathpat.set_bit(idx);
        }
        if let Some(idx) = graph.edge_bit_index(2, graph.sink_id) {
            pathpat.set_bit(idx);
        }
        let mut path = vec![1usize];
        let mut minpath = 1usize;
        let mut maxpath = 2usize;

        let reached = fwd_to_sink_fast_long_st(
            &graph,
            &mut transfrags,
            &mut nodecov,
            &mut pathpat,
            &mut path,
            &mut minpath,
            &mut maxpath,
            1,
        );
        assert!(reached, "expected to reach sink");
        assert_eq!(
            path.last().copied(),
            Some(graph.sink_id),
            "expected last node to be sink"
        );
        // Path should be [1, 2, sink].
        assert!(path.contains(&2), "expected to pick exon B (node 2)");
    }

    /// Symmetric back_to_source test: 4-node graph (source → A → B → sink),
    /// long-read transfrag [1, 2]. Walk from B back to source.
    #[test]
    fn back_to_source_st_minimal_graph() {
        let mut graph = crate::graph::Graph::new();
        let _ = graph.add_node(0, 0);
        let _ = graph.add_node(100, 200);
        let _ = graph.add_node(300, 400);
        let _ = graph.add_node(0, 0);
        graph.source_id = 0;
        graph.sink_id = 3;
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(2, 3);
        graph.compute_reachability();
        let pat_size = graph.pattern_size();

        let mut tf = GraphTransfrag::new(vec![1, 2], pat_size);
        tf.abundance = 10.0;
        tf.longread = true;
        tf.weak = 0;
        tf.pattern.insert_grow(1);
        tf.pattern.insert_grow(2);
        graph.set_pattern_edges_for_path(&mut tf.pattern, &[1, 2]);
        let mut transfrags = vec![tf];
        graph.nodes[1].trf_ids.push(0);
        graph.nodes[2].trf_ids.push(0);

        let mut nodecov = vec![0.0; graph.nodes.len()];
        nodecov[1] = 10.0;
        nodecov[2] = 10.0;

        let mut pathpat = transfrags[0].pattern.clone();
        if let Some(idx) = graph.edge_bit_index(graph.source_id, 1) {
            pathpat.set_bit(idx);
        }
        if let Some(idx) = graph.edge_bit_index(2, graph.sink_id) {
            pathpat.set_bit(idx);
        }
        let mut path = vec![2usize];
        let mut minpath = 1usize;
        let mut maxpath = 2usize;

        // Walk from node 2 back to source (= node 0).
        let reached = back_to_source_fast_long_st(
            &graph,
            &mut transfrags,
            &mut nodecov,
            &mut pathpat,
            &mut path,
            &mut minpath,
            &mut maxpath,
            2,
        );
        assert!(reached, "expected to reach source from node 2");
        // Path before walk: [2]. Walk should add [1] (not source — ST's `if(maxp)` rule).
        assert!(path.contains(&1), "expected to pick node 1 on the way back");
        assert!(
            !path.contains(&graph.source_id),
            "ST convention: source not added to path"
        );
    }
}
